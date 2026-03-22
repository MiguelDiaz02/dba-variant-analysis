"""
domain_analysis.py
==================
Análisis standalone de dominios proteicos en variantes DBA.

Pregunta biológica principal: ¿caen las variantes patogénicas dentro de dominios
anotados de las proteínas ribosomales afectadas en DBA?

Tres sub-análisis:
    1. Cohorte local sola (14 variantes, 12 con posición resuelta)
    2. Cohorte local + ClinVar patogénicas (si variant_summary.txt.gz disponible)
    3. Tabla comparativa: % intra-dominio local vs ClinVar por gen

Figuras generadas:
    graphical_results2/domain_architecture_all.{png,svg}   — Fig 1 global (7 proteínas)
    graphical_results2/domain_track_{GENE}.{png,svg}        — Fig 2 por gen × 7
    graphical_results2/domain_lollipop_{GENE}.{png,svg}     — Fig 3 por gen × 7

Reporte:
    reports/domain_analysis.md

Uso:
    python domain_analysis.py
    python domain_analysis.py --skip-clinvar --verbose
    python domain_analysis.py --genes RPL5 RPL11
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import plotly.graph_objects as go
import yaml

sys.path.insert(0, str(Path(__file__).parent))
import annotator
import data_loader
import visualizer

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------

# Tipos de feature a incluir en el análisis (excluye estructura secundaria)
_DOMAIN_FEATURE_TYPES = frozenset({
    "Domain", "Region", "Modified residue", "Motif", "Binding site", "Active site",
})

# Tipos que cuentan para intra/inter-dominio (excluye PTMs)
_OVERLAP_FEATURE_TYPES = frozenset({
    "Domain", "Region", "Motif", "Binding site", "Active site",
})

# PTMs — solo se marcan cuando una variante local coincide exactamente
_PTM_FEATURE_TYPES = frozenset({"Modified residue"})

# Prioridad para resolver solapamientos múltiples (menor = mayor prioridad)
_FEATURE_PRIORITY: dict[str, int] = {
    "Domain": 0, "Region": 1, "Motif": 2, "Binding site": 3, "Active site": 4,
}

# Bandas de fondo por subunidad
_SUBUNIT_BAND_COLORS: dict[str, str] = {
    "40S": "rgba(173, 216, 230, 0.18)",
    "60S": "rgba(144, 238, 144, 0.18)",
}

# Símbolo según ACMG para Fig 3 (lollipop)
_ACMG_SYMBOL: dict[str, str] = {
    "P":   "circle",
    "LP":  "circle-open",
    "VUS": "diamond-open",
}

_NAVY = "#0A2342"


# ---------------------------------------------------------------------------
# Logging (copiado de pipeline.py para mantener standalone)
# ---------------------------------------------------------------------------

def _setup_logging(verbose: bool, log_file: Path) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    fmt = "%(asctime)s [%(levelname)-8s] %(name)s — %(message)s"
    stream_handler = logging.StreamHandler(sys.stdout)
    if hasattr(stream_handler.stream, "reconfigure"):
        stream_handler.stream.reconfigure(encoding="utf-8")
    handlers: list[logging.Handler] = [
        stream_handler,
        logging.FileHandler(log_file, encoding="utf-8"),
    ]
    for h in handlers:
        h.setFormatter(logging.Formatter(fmt, "%H:%M:%S"))
        h.setLevel(level)
    root = logging.getLogger()
    root.setLevel(level)
    for h in root.handlers[:]:
        root.removeHandler(h)
    for h in handlers:
        root.addHandler(h)


# ---------------------------------------------------------------------------
# Configuración (copiado de pipeline.py)
# ---------------------------------------------------------------------------

def _load_config(config_path: Path) -> dict:
    if not config_path.exists():
        raise FileNotFoundError(f"Config no encontrado: {config_path}")
    with config_path.open(encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _resolve_paths(cfg: dict, base: Path) -> dict[str, Path]:
    paths = cfg.get("paths", {})
    log_base = (base / paths.get("log_file", "../pipeline.log")).resolve()
    return {
        "csv":       (base / paths["variants_csv"]).resolve(),
        "clinvar":   (base / paths["clinvar_gz"]).resolve(),
        "out_viz":   (base / paths.get("graphical_results",  "../graphical_results")).resolve(),
        "out_viz2":  (base / paths.get("graphical_results2", "../graphical_results2")).resolve(),
        "out_rep":   (base / paths.get("reports", "../reports")).resolve(),
        "log_file":  log_base,
        "log_domain": log_base.parent / "domain_analysis.log",
    }


def _build_gene_dicts(cfg: dict) -> tuple[list[str], dict, dict, dict]:
    genes_list, prot_len, subunit_map, uniprot_ids = [], {}, {}, {}
    for entry in cfg.get("genes", []):
        name = entry["name"]
        genes_list.append(name)
        prot_len[name]    = entry["protein_length"]
        subunit_map[name] = entry["subunit"]
        uniprot_ids[name] = entry.get("uniprot_id", "")
    return genes_list, prot_len, subunit_map, uniprot_ids


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="domain_analysis.py",
        description="Análisis de dominios proteicos en variantes DBA — script standalone",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", type=Path,
                   default=Path(__file__).parent / "config.yaml")
    p.add_argument("--skip-clinvar", action="store_true",
                   help="Omitir carga de ClinVar.")
    p.add_argument("--genes", nargs="+", metavar="GEN",
                   help="Analizar solo estos genes.")
    p.add_argument("--verbose", "-v", action="store_true")
    return p


# ---------------------------------------------------------------------------
# Carga de dominios UniProt
# ---------------------------------------------------------------------------

def _load_all_domains(
    genes_list: list[str],
    uniprot_ids: dict[str, str],
) -> dict[str, list[dict]]:
    """
    Obtiene features UniProt para cada gen y filtra a tipos relevantes.
    Devuelve {gene: [feat_dict, ...]} donde feat_dict tiene type, description, start, end.
    """
    features_by_gene: dict[str, list[dict]] = {}
    for gene in genes_list:
        uid = uniprot_ids.get(gene, "")
        if not uid:
            log.warning("%s: sin UniProt ID en config.yaml", gene)
            features_by_gene[gene] = []
            continue
        raw = annotator.get_uniprot_features(uid)
        filtered = [f for f in raw if f["type"] in _DOMAIN_FEATURE_TYPES]
        log.info("%s (%s): %d features relevantes (de %d raw)", gene, uid, len(filtered), len(raw))
        features_by_gene[gene] = filtered
    return features_by_gene


# ---------------------------------------------------------------------------
# Mapeo variante → dominio
# ---------------------------------------------------------------------------

def _map_variants_to_domains(
    variants_df: pd.DataFrame,
    features_by_gene: dict[str, list[dict]],
) -> pd.DataFrame:
    """
    Añade columnas domain_context, domain_type y ptm_overlap al DataFrame.

    - domain_context: nombre del dominio solapante (prioridad Domain > Region > ...)
                      o "Inter-domain" si no hay solapamiento.
                      "Sin posición" para variantes sin posición resuelta.
    - domain_type:    tipo de feature solapante, o "" si inter-domain / sin posición.
    - ptm_overlap:    True si la posición coincide exactamente con un Modified residue.
    """
    out = variants_df.copy()
    contexts, types, ptm_flags = [], [], []

    for _, row in out.iterrows():
        pos = row.get("position")
        gene = row.get("Gene", "")
        features = features_by_gene.get(gene, [])

        if pd.isna(pos) or pos is None:
            contexts.append("Sin posición")
            types.append("")
            ptm_flags.append(False)
            continue

        pos = int(pos)

        # Verificar solapamiento con PTM (posición exacta)
        ptm_hit = any(
            f["type"] in _PTM_FEATURE_TYPES and f["start"] <= pos <= f["end"]
            for f in features
        )

        # Buscar features de solapamiento (no PTM)
        overlapping = [
            f for f in features
            if f["type"] in _OVERLAP_FEATURE_TYPES and f["start"] <= pos <= f["end"]
        ]

        if overlapping:
            best = min(overlapping, key=lambda f: _FEATURE_PRIORITY.get(f["type"], 99))
            desc = (best.get("description") or "").strip() or best["type"]
            contexts.append(desc)
            types.append(best["type"])
        else:
            contexts.append("Inter-domain")
            types.append("")

        ptm_flags.append(ptm_hit)

    out["domain_context"] = contexts
    out["domain_type"] = types
    out["ptm_overlap"] = ptm_flags
    return out


def _map_positions_to_domains(
    positions: list[int],
    features: list[dict],
) -> tuple[int, int]:
    """
    Cuenta cuántas posiciones ClinVar caen dentro de dominios (_OVERLAP_FEATURE_TYPES).
    Devuelve (n_intra, n_inter).
    """
    n_intra = sum(
        1 for pos in positions
        if any(
            f["type"] in _OVERLAP_FEATURE_TYPES and f["start"] <= pos <= f["end"]
            for f in features
        )
    )
    return n_intra, len(positions) - n_intra


# ---------------------------------------------------------------------------
# ClinVar (carga con fallback gracioso)
# ---------------------------------------------------------------------------

def _load_clinvar_data(
    resolved: dict[str, Path],
    cfg: dict,
    skip_clinvar: bool,
) -> Optional[pd.DataFrame]:
    if skip_clinvar:
        log.info("ClinVar omitido (--skip-clinvar)")
        return None
    try:
        df = annotator.load_clinvar_file(resolved["clinvar"])
        return annotator.filter_clinvar_pathogenic(df)
    except FileNotFoundError as exc:
        log.warning("Archivo ClinVar no encontrado — análisis solo con cohorte local: %s", exc)
        return None
    except Exception as exc:
        log.error("Error inesperado cargando ClinVar: %s", exc)
        return None


# ---------------------------------------------------------------------------
# Estadísticas de dominios
# ---------------------------------------------------------------------------

def _compute_domain_stats(
    annotated_df: pd.DataFrame,
    features_by_gene: dict[str, list[dict]],
    clinvar_positions_by_gene: dict[str, list[int]],
    prot_len: dict[str, int],
    subunit_map: dict[str, str],
) -> dict:
    """
    Calcula estadísticas de solapamiento dominio por gen y por subunidad.

    Devuelve dict con claves:
        per_gene   → {gene: stats_dict}
        per_subunit → {"40S": ..., "60S": ...}
    """
    per_gene: dict[str, dict] = {}

    for gene, features in features_by_gene.items():
        gene_df = annotated_df[annotated_df["Gene"] == gene].copy()
        with_pos = gene_df[gene_df["domain_context"] != "Sin posición"]
        n_total   = len(gene_df)
        n_with_pos = len(with_pos)
        n_intra = int((with_pos["domain_context"] != "Inter-domain").sum())
        pct_intra = round(n_intra / n_with_pos * 100, 1) if n_with_pos > 0 else 0.0

        # PTM overlaps
        ptm_variants = gene_df[gene_df["ptm_overlap"] == True][["AA_change", "domain_context"]].to_dict("records")

        # Breakdown por dominio
        domain_breakdown: dict[str, int] = {}
        for _, row in with_pos.iterrows():
            ctx = row["domain_context"]
            if ctx not in ("Inter-domain", "Sin posición"):
                domain_breakdown[ctx] = domain_breakdown.get(ctx, 0) + 1

        # ClinVar
        clinvar_pos = clinvar_positions_by_gene.get(gene, [])
        n_cv = len(clinvar_pos)
        if n_cv > 0 and features:
            cv_intra, cv_inter = _map_positions_to_domains(clinvar_pos, features)
        else:
            cv_intra = 0
        pct_cv_intra = round(cv_intra / n_cv * 100, 1) if n_cv > 0 else None

        per_gene[gene] = {
            "subunit":          subunit_map.get(gene, "?"),
            "n_local_total":    n_total,
            "n_local_with_pos": n_with_pos,
            "n_local_intra":    n_intra,
            "pct_local_intra":  pct_intra,
            "n_clinvar_total":  n_cv,
            "n_clinvar_intra":  cv_intra,
            "pct_clinvar_intra": pct_cv_intra,
            "ptm_variants":     ptm_variants,
            "domain_breakdown": domain_breakdown,
            "n_features":       len(features),
        }

    # Por subunidad
    per_subunit: dict[str, dict] = {}
    for su in ("40S", "60S"):
        su_genes = [g for g, d in per_gene.items() if d["subunit"] == su]
        n_wp  = sum(per_gene[g]["n_local_with_pos"] for g in su_genes)
        n_in  = sum(per_gene[g]["n_local_intra"]    for g in su_genes)
        n_cv  = sum(per_gene[g]["n_clinvar_total"]  for g in su_genes)
        cv_in = sum(per_gene[g]["n_clinvar_intra"]  for g in su_genes)
        per_subunit[su] = {
            "n_local_with_pos": n_wp,
            "n_local_intra":    n_in,
            "pct_local_intra":  round(n_in / n_wp * 100, 1) if n_wp > 0 else 0.0,
            "n_clinvar_total":  n_cv,
            "n_clinvar_intra":  cv_in,
            "pct_clinvar_intra": round(cv_in / n_cv * 100, 1) if n_cv > 0 else None,
        }

    return {"per_gene": per_gene, "per_subunit": per_subunit}


# ---------------------------------------------------------------------------
# Figura 3: Lollipop con fondo de dominios (por gen)
# ---------------------------------------------------------------------------

def _fig_domain_lollipop(
    gene: str,
    gene_df: pd.DataFrame,
    features: list[dict],
    protein_length: int,
) -> go.Figure:
    """
    Lollipop plot con rectángulos de dominio como fondo semitransparente.
    Color del marcador = variant_type; símbolo = ACMG.
    """
    fig = go.Figure()

    # Rectángulos de dominio (fondo, layer=below)
    # Dibujar por prioridad inversa: Domain encima
    sorted_feats = sorted(
        [f for f in features if f["type"] in _OVERLAP_FEATURE_TYPES],
        key=lambda f: _FEATURE_PRIORITY.get(f["type"], 99),
        reverse=True,
    )
    for feat in sorted_feats:
        fcolor = visualizer._FEATURE_TYPE_COLORS.get(feat["type"], visualizer._FEATURE_TYPE_DEFAULT)
        fig.add_shape(
            type="rect",
            x0=feat["start"], x1=feat["end"],
            y0=-0.3, y1=1.9,
            fillcolor=fcolor, opacity=0.18,
            line_width=0, layer="below",
        )
        mid = (feat["start"] + feat["end"]) / 2.0
        desc = (feat.get("description") or feat["type"])[:20]
        fig.add_annotation(
            x=mid, y=1.82,
            text=f"<span style='font-size:8px'>{desc}</span>",
            showarrow=False,
            font=dict(size=8, color=fcolor),
            bgcolor="rgba(255,255,255,0.7)",
            borderpad=1,
        )

    # Línea base de la proteína
    fig.add_shape(
        type="line", x0=0, x1=protein_length, y0=0, y1=0,
        line=dict(color="black", width=2.5),
    )

    # Variantes: stems + heads agrupados por (variant_type, ACMG)
    df_local = gene_df.dropna(subset=["position"]).copy()
    for vtype, vcolor in visualizer._VARIANT_TYPE_COLORS.items():
        sub_vtype = df_local[df_local["variant_type"] == vtype]
        if sub_vtype.empty:
            continue
        for acmg, symbol in _ACMG_SYMBOL.items():
            sub = sub_vtype[sub_vtype["ACMG"] == acmg]
            if sub.empty:
                continue
            positions = sub["position"].tolist()
            aa_changes = sub["AA_change"].fillna("").tolist()
            domains_ctx = sub["domain_context"].tolist() if "domain_context" in sub.columns else ["?"] * len(sub)

            for pos in positions:
                fig.add_shape(
                    type="line", x0=pos, x1=pos, y0=0, y1=1,
                    line=dict(color=vcolor, width=1.5),
                )

            fill = vcolor if "open" not in symbol else "rgba(255,255,255,0)"
            fig.add_trace(go.Scatter(
                x=positions, y=[1] * len(positions),
                mode="markers",
                marker=dict(
                    symbol=symbol, size=12, color=fill,
                    line=dict(color=vcolor, width=2),
                ),
                name=f"{vtype} / {acmg}",
                legendgroup=vtype,
                customdata=list(zip(aa_changes, domains_ctx, [acmg]*len(positions))),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "Tipo: " + vtype + "<br>"
                    "ACMG: %{customdata[2]}<br>"
                    "Dominio: %{customdata[1]}"
                    "<extra></extra>"
                ),
            ))

    n_local = int(df_local.shape[0])
    fig.update_layout(
        title=dict(
            text=f"<b>{gene}</b> — Variantes en contexto de dominio (n={n_local})",
            font=dict(size=15, color=_NAVY),
        ),
        xaxis=dict(
            title="Posición aminoacídica",
            range=[-2, protein_length + 3],
            showgrid=True, gridcolor="rgba(200,200,200,0.4)", zeroline=False,
        ),
        yaxis=dict(visible=False, range=[-0.5, 2.2]),
        plot_bgcolor="white", paper_bgcolor="white",
        height=480,
        legend=dict(
            title="Tipo / ACMG (símbolo = ACMG)",
            orientation="v", x=1.02, y=1,
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor="lightgray", borderwidth=1,
        ),
        margin=dict(r=220, t=70, b=60),
    )
    fig.update_xaxes(showgrid=True, gridcolor="lightgray", zeroline=False)
    return fig


# ---------------------------------------------------------------------------
# Figura 2: Track por gen (reutiliza visualizer._track_fig)
# ---------------------------------------------------------------------------

def _fig_domain_track(
    gene: str,
    gene_df: pd.DataFrame,
    clinvar_positions: list[int],
    features: list[dict],
    protein_length: int,
    acmg_colors: dict,
) -> go.Figure:
    """
    Figura tipo 'genomic track' standalone.
    Pasa TODOS los features filtrados como 'domains' para mostrarlos completos.
    """
    return visualizer._track_fig(
        gene_df=gene_df,
        clinvar_positions=clinvar_positions,
        hotspots=[],          # no hay hotspots en este análisis standalone
        domains=features,     # todos los features, no solo los solapantes con hotspots
        protein_length=protein_length,
        gene=gene,
        acmg_colors=acmg_colors,
    )


# ---------------------------------------------------------------------------
# Figura 1: Diagrama de arquitectura global (todos los genes)
# ---------------------------------------------------------------------------

def _fig_domain_architecture(
    annotated_df: pd.DataFrame,
    clinvar_positions_by_gene: dict[str, list[int]],
    features_by_gene: dict[str, list[dict]],
    prot_len: dict[str, int],
    subunit_map: dict[str, str],
    ordered_genes: list[str],
) -> go.Figure:
    """
    Figura global con todas las proteínas como barras horizontales.
    Eje X normalizado [0, 1] (N-terminus → C-terminus).
    Cada fila = un gen, ordenadas 40S arriba → 60S abajo.
    """
    fig = go.Figure()
    n_genes = len(ordered_genes)
    row_gap = 1.5

    # --- Bandas de fondo por subunidad ---
    for su in ("40S", "60S"):
        su_indices = [i for i, g in enumerate(ordered_genes) if subunit_map.get(g) == su]
        if not su_indices:
            continue
        # y_center para gen i: (n_genes - 1 - i) * row_gap
        y_top = (n_genes - 1 - min(su_indices)) * row_gap + 0.8
        y_bot = (n_genes - 1 - max(su_indices)) * row_gap - 0.8
        fig.add_shape(
            type="rect", xref="paper", yref="y",
            x0=0, x1=1, y0=y_bot, y1=y_top,
            fillcolor=_SUBUNIT_BAND_COLORS[su],
            line_width=0, layer="below",
        )
        # Etiqueta subunidad
        y_mid = (y_top + y_bot) / 2
        fig.add_annotation(
            xref="paper", x=-0.11, y=y_mid,
            text=f"<b>{su}</b>", showarrow=False,
            font=dict(size=13, color=_NAVY),
            xanchor="center",
        )

    # Rastrear qué tipos de feature y variant_type aparecen (para leyenda)
    seen_feat_types: set[str] = set()
    seen_vtypes: set[str] = set()

    for i, gene in enumerate(ordered_genes):
        y_center = (n_genes - 1 - i) * row_gap
        plen = prot_len[gene]
        features = features_by_gene.get(gene, [])
        gene_df = annotated_df[annotated_df["Gene"] == gene].copy()

        # --- Barra de proteína (fondo gris) ---
        fig.add_shape(
            type="rect", x0=0, x1=1,
            y0=y_center - 0.14, y1=y_center + 0.14,
            fillcolor="#DDDDDD", line_color="#BDBDBD", line_width=0.8,
        )

        # --- Rectángulos de dominio ---
        sorted_feats = sorted(
            [f for f in features if f["type"] in _OVERLAP_FEATURE_TYPES],
            key=lambda f: _FEATURE_PRIORITY.get(f["type"], 99),
            reverse=True,
        )
        for feat in sorted_feats:
            fcolor = visualizer._FEATURE_TYPE_COLORS.get(feat["type"], visualizer._FEATURE_TYPE_DEFAULT)
            x0n = feat["start"] / plen
            x1n = feat["end"] / plen
            fig.add_shape(
                type="rect", x0=x0n, x1=x1n,
                y0=y_center - 0.24, y1=y_center + 0.24,
                fillcolor=fcolor, opacity=0.78,
                line_color="white", line_width=0.5,
            )
            seen_feat_types.add(feat["type"])
            # Tooltip: trace invisible en el centro del dominio
            mid_n = (x0n + x1n) / 2
            desc = (feat.get("description") or feat["type"])[:30]
            fig.add_trace(go.Scatter(
                x=[mid_n], y=[y_center],
                mode="markers",
                marker=dict(color=fcolor, size=1, opacity=0),
                text=[f"{feat['type']}: {desc}<br>aa {feat['start']}–{feat['end']}"],
                hovertemplate="%{text}<extra></extra>",
                showlegend=False, name="",
            ))

        # --- Etiqueta del gen (izquierda) ---
        fig.add_annotation(
            x=-0.02, y=y_center,
            text=f"<b>{gene}</b>",
            showarrow=False, xanchor="right",
            font=dict(size=11, color=_NAVY),
        )
        # Longitud de proteína (derecha)
        fig.add_annotation(
            x=1.02, y=y_center,
            text=f"<span style='font-size:9px; color:#888'>{plen} aa</span>",
            showarrow=False, xanchor="left",
            font=dict(size=9),
        )

        # --- ClinVar: ticks debajo de la barra ---
        clinvar_pos = clinvar_positions_by_gene.get(gene, [])
        if clinvar_pos:
            cv_norm = [p / plen for p in clinvar_pos]
            fig.add_trace(go.Scatter(
                x=cv_norm,
                y=[y_center - 0.40] * len(cv_norm),
                mode="markers",
                marker=dict(symbol="line-ns", size=8, color="#9E9E9E",
                            line=dict(width=1, color="#9E9E9E"), opacity=0.35),
                name="ClinVar (P/LP)",
                legendgroup="clinvar",
                showlegend=(i == 0),
                hovertemplate=f"ClinVar {gene}<extra></extra>",
            ))

        # --- Variantes locales: líneas + marcadores encima de la barra ---
        df_local = gene_df.dropna(subset=["position"])
        for vtype, vcolor in visualizer._VARIANT_TYPE_COLORS.items():
            sub = df_local[df_local["variant_type"] == vtype]
            if sub.empty:
                continue
            seen_vtypes.add(vtype)
            pos_norm = (sub["position"] / plen).tolist()
            aa_list = sub["AA_change"].fillna("").tolist()
            acmg_list = sub["ACMG"].tolist()
            ctx_list = sub["domain_context"].tolist() if "domain_context" in sub.columns else ["?"] * len(sub)

            # Stems (líneas verticales hacia arriba)
            for pn in pos_norm:
                fig.add_shape(
                    type="line", x0=pn, x1=pn,
                    y0=y_center + 0.24, y1=y_center + 0.50,
                    line=dict(color=vcolor, width=1.8),
                )

            fig.add_trace(go.Scatter(
                x=pos_norm,
                y=[y_center + 0.50] * len(pos_norm),
                mode="markers",
                marker=dict(symbol="triangle-down", size=9, color=vcolor,
                            line=dict(color=vcolor, width=1)),
                name=vtype,
                legendgroup=f"vtype_{vtype}",
                showlegend=(i == 0),
                customdata=list(zip(aa_list, acmg_list, ctx_list)),
                hovertemplate=(
                    f"<b>%{{customdata[0]}}</b> ({gene})<br>"
                    "ACMG: %{customdata[1]}<br>"
                    "Tipo: " + vtype + "<br>"
                    "Dominio: %{customdata[2]}"
                    "<extra></extra>"
                ),
            ))

        # --- Sitios PTM con variante local (marcador diamante en la barra) ---
        ptm_sub = gene_df[gene_df.get("ptm_overlap", pd.Series(False, index=gene_df.index)) == True] if "ptm_overlap" in gene_df.columns else pd.DataFrame()
        if not ptm_sub.empty:
            ptm_pos = ptm_sub.dropna(subset=["position"])
            if not ptm_pos.empty:
                ptm_norm = (ptm_pos["position"] / plen).tolist()
                ptm_aa = ptm_pos["AA_change"].fillna("").tolist()
                fig.add_trace(go.Scatter(
                    x=ptm_norm, y=[y_center] * len(ptm_norm),
                    mode="markers",
                    marker=dict(symbol="diamond", size=8, color="#E53935",
                                line=dict(color="darkred", width=1.5)),
                    name="PTM site (variante local)",
                    legendgroup="ptm",
                    showlegend=(i == 0),
                    text=ptm_aa,
                    hovertemplate="PTM: <b>%{text}</b><extra></extra>",
                ))

    # --- Leyenda fantasma para tipos de feature (rectángulos de dominio) ---
    for ftype in sorted(seen_feat_types, key=lambda t: _FEATURE_PRIORITY.get(t, 99)):
        fcolor = visualizer._FEATURE_TYPE_COLORS.get(ftype, visualizer._FEATURE_TYPE_DEFAULT)
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode="markers",
            marker=dict(symbol="square", size=12, color=fcolor, opacity=0.78),
            name=ftype,
            legendgroup=f"feat_{ftype}",
            showlegend=True,
        ))

    fig.update_layout(
        title=dict(
            text="<b>Arquitectura de dominios — Proteínas ribosomales DBA</b>",
            font=dict(size=16, color=_NAVY),
        ),
        xaxis=dict(
            title="Posición normalizada (N-terminus → C-terminus)",
            range=[-0.04, 1.06],
            tickvals=[0, 0.25, 0.5, 0.75, 1.0],
            ticktext=["0%", "25%", "50%", "75%", "100%"],
            showgrid=True, gridcolor="rgba(200,200,200,0.4)", zeroline=False,
        ),
        yaxis=dict(
            visible=False,
            range=[(n_genes - 1 - (n_genes - 1)) * row_gap - 1.0,
                   (n_genes - 1) * row_gap + 1.0],
        ),
        plot_bgcolor="white", paper_bgcolor="white",
        height=max(600, 140 * n_genes + 100),
        width=1400,
        legend=dict(
            title=dict(text="Leyenda", font=dict(size=11)),
            orientation="v", x=1.06, y=1,
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="lightgray", borderwidth=1,
            tracegroupgap=6,
        ),
        margin=dict(l=130, r=240, t=60, b=60),
    )
    return fig


# ---------------------------------------------------------------------------
# Reporte Markdown
# ---------------------------------------------------------------------------

def _pct(n: int, d: int) -> str:
    return f"{round(n / d * 100, 1):.1f}%" if d > 0 else "N/A"


def _generate_report(
    stats: dict,
    features_by_gene: dict[str, list[dict]],
    annotated_df: pd.DataFrame,
    output_path: Path,
    has_clinvar: bool,
    figure_paths: list[Path],
) -> None:
    per_gene = stats["per_gene"]
    per_sub  = stats["per_subunit"]
    lines: list[str] = []

    # ── Encabezado ──────────────────────────────────────────────────────────
    lines += [
        "# Domain Analysis: Variant–Domain Mapping in DBA Ribosomal Proteins",
        "",
        "## 1. Biological Rationale",
        "",
        "Ribosomal proteins in Diamond-Blackfan anemia are organized into defined "
        "structural and functional domains that mediate critical interactions with "
        "ribosomal RNA and assembly factors. Pathogenic variants in DBA preferentially "
        "disrupt protein function, suggesting that mutations may concentrate within "
        "structurally or functionally indispensable regions.",
        "",
        "For the 60S compartment, RPL5 and RPL11 participate in the 5S ribonucleoprotein "
        "(5S RNP) complex that stabilizes p53 in response to ribosomal stress via MDM2 "
        "inhibition. The uL18 fold of RPL5 — responsible for 5S rRNA binding — represents "
        "a domain where missense variants may abrogate function without requiring truncation. "
        "For 40S proteins, RNA-contact surfaces at the platform and beak regions of the "
        "small subunit are candidate hotspot domains.",
        "",
        "This analysis maps each variant to its domain context using UniProt annotations, "
        "addressing whether pathogenic mutations preferentially co-localize with known "
        "structural or functional domains.",
        "",
    ]

    # ── Métodos ─────────────────────────────────────────────────────────────
    overlap_feature_types_str = ", ".join(f'"{t}"' for t in sorted(_OVERLAP_FEATURE_TYPES))
    all_feature_types_str = ", ".join(f'"{t}"' for t in sorted(_DOMAIN_FEATURE_TYPES))
    lines += [
        "## 2. Methods",
        "",
        "### 2.1 Data Source",
        "",
        "Protein domain annotations were obtained from the UniProt Knowledge Base REST API "
        f"(https://rest.uniprot.org/uniprotkb/{{UniProtID}}.json) for each of the seven "
        "genes analyzed. UniProt consolidates domain information from Pfam, SMART, "
        "SUPERFAMILY, and PANTHER, providing a curated, non-redundant annotation set.",
        "",
        "### 2.2 Feature Types Included",
        "",
        f"All features of the following types were retrieved: {all_feature_types_str}. "
        "Secondary structure elements (Helix, Beta strand, Turn) were explicitly excluded "
        "to focus on structurally and functionally defined regions.",
        "",
        f"For intra/inter-domain classification, the following types were used: "
        f"{overlap_feature_types_str}. "
        '"Modified residue" (PTM) sites were loaded but are only reported when a local '
        "variant falls at the exact annotated residue.",
        "",
        "### 2.3 Variant–Domain Mapping",
        "",
        "For each variant with a resolved amino acid position, domain overlap was assessed "
        "as: `feature_start ≤ variant_position ≤ feature_end`. When multiple features "
        "overlapped, the highest-priority type was retained "
        "(Domain > Region > Motif > Binding site > Active site). Variants without a resolved "
        'position (e.g., splice variants) are reported as "Sin posición".',
        "",
        "### 2.4 Normalization",
        "",
        "The domain architecture figure (Figure 1) uses normalized positions "
        "(`position / protein_length ∈ [0, 1]`) to enable visual comparison across "
        "proteins of different lengths.",
        "",
    ]

    # ── Resultados por gen ──────────────────────────────────────────────────
    lines += ["## 3. Per-Gene Results", ""]
    for gene, d in per_gene.items():
        subunit = d["subunit"]
        n_feats = d["n_features"]
        lines += [
            f"### {gene} ({subunit}, {per_gene[gene].get('subunit','?')})",
            "",
        ]
        if n_feats == 0:
            lines += [
                f"> **Note:** No annotated features of relevant types were returned by "
                f"UniProt for {gene}. All variants will be classified as Inter-domain.",
                "",
            ]

        gene_df = annotated_df[annotated_df["Gene"] == gene].copy()
        if gene_df.empty:
            lines += ["No variants found for this gene in the local cohort.", ""]
            continue

        lines += [
            "| Variant (AA change) | Position | ACMG | Type | Domain context | Domain type | PTM overlap |",
            "|---------------------|----------|------|------|----------------|-------------|-------------|",
        ]
        for _, row in gene_df.iterrows():
            pos_str = str(int(row["position"])) if not pd.isna(row.get("position")) else "—"
            ctx = row.get("domain_context", "?")
            dtype = row.get("domain_type", "") or "—"
            ptm = "Yes ⚠" if row.get("ptm_overlap") else "No"
            lines.append(
                f"| {row.get('AA_change','?')} | {pos_str} | {row.get('ACMG','?')} "
                f"| {row.get('variant_type','?')} | {ctx} | {dtype} | {ptm} |"
            )
        lines.append("")

    # ── Tabla resumen ────────────────────────────────────────────────────────
    lines += [
        "## 4. Summary Statistics",
        "",
        "### 4.1 Per-Gene",
        "",
    ]
    if has_clinvar:
        lines += [
            "| Gene | Subunit | Local (n) | % Intra-domain (local) | ClinVar (n) | % Intra-domain (ClinVar) |",
            "|------|---------|-----------|------------------------|-------------|--------------------------|",
        ]
        for gene, d in per_gene.items():
            cv_pct = f"{d['pct_clinvar_intra']:.1f}%" if d["pct_clinvar_intra"] is not None else "N/A"
            lines.append(
                f"| {gene} | {d['subunit']} | {d['n_local_with_pos']} | "
                f"{d['pct_local_intra']:.1f}% | {d['n_clinvar_total']} | {cv_pct} |"
            )
    else:
        lines += [
            "| Gene | Subunit | Local (n with position) | Intra-domain (n) | % Intra-domain |",
            "|------|---------|------------------------|------------------|----------------|",
        ]
        for gene, d in per_gene.items():
            lines.append(
                f"| {gene} | {d['subunit']} | {d['n_local_with_pos']} | "
                f"{d['n_local_intra']} | {d['pct_local_intra']:.1f}% |"
            )
    lines.append("")

    lines += ["### 4.2 Per-Subunit", ""]
    if has_clinvar:
        lines += [
            "| Subunit | Local (n) | % Intra-domain | ClinVar (n) | % Intra-domain (ClinVar) |",
            "|---------|-----------|----------------|-------------|--------------------------|",
        ]
        for su, d in per_sub.items():
            cv_pct = f"{d['pct_clinvar_intra']:.1f}%" if d["pct_clinvar_intra"] is not None else "N/A"
            lines.append(
                f"| {su} | {d['n_local_with_pos']} | {d['pct_local_intra']:.1f}% | "
                f"{d['n_clinvar_total']} | {cv_pct} |"
            )
    else:
        lines += [
            "| Subunit | Local (n with position) | Intra-domain (n) | % Intra-domain |",
            "|---------|------------------------|------------------|----------------|",
        ]
        for su, d in per_sub.items():
            lines.append(
                f"| {su} | {d['n_local_with_pos']} | {d['n_local_intra']} | "
                f"{d['pct_local_intra']:.1f}% |"
            )
    lines.append("")

    # ── PTMs ─────────────────────────────────────────────────────────────────
    lines += ["## 5. PTM Co-localization Findings", ""]
    ptm_found = False
    for gene, d in per_gene.items():
        if d["ptm_variants"]:
            ptm_found = True
            for pv in d["ptm_variants"]:
                lines.append(
                    f"- **{gene}** — `{pv['AA_change']}` co-localizes with a known "
                    f"post-translational modification site ({pv['domain_context']})."
                )
    if not ptm_found:
        lines.append(
            "No local variants co-localize with known post-translational modification sites "
            "annotated in UniProt for these ribosomal proteins."
        )
    lines.append("")

    # ── Comparación Local vs ClinVar ─────────────────────────────────────────
    lines += ["## 6. Local vs ClinVar Comparison", ""]
    if not has_clinvar:
        lines += [
            "ClinVar analysis was not performed (--skip-clinvar flag or file unavailable). "
            "Run with `variant_summary.txt.gz` present to enable this comparison.",
            "",
        ]
    else:
        lines += [
            "The table in Section 4 compares the proportion of intra-domain variants between "
            "the local cohort and ClinVar pathogenic variants. Concordance between the two "
            "sources would support the hypothesis that domain co-localization is a general "
            "property of pathogenic DBA mutations rather than a sampling artefact.",
            "",
            "> **Caveat:** The local cohort is small (n = 12–14 with resolved positions). "
            "Percentage comparisons are descriptive only.",
            "",
        ]

    # ── Recomendaciones de figuras ────────────────────────────────────────────
    lines += [
        "## 7. Figure Recommendations for Article Inclusion",
        "",
        "Three figure styles were generated to support different narrative contexts:",
        "",
        "### Figure 1 — Domain Architecture (`domain_architecture_all`)",
        "**Recommended for the article main body.** This figure provides an integrated "
        "cross-protein overview showing the relationship between annotated domain structure "
        "and variant positions for all seven proteins simultaneously. The normalized X-axis "
        "enables direct visual comparison despite heterogeneous protein lengths. Suitable "
        "as a full-width figure in a clinical genetics manuscript.",
        "",
        "### Figure 3 — Lollipop with Domain Background (`domain_lollipop_{GENE}`)",
        "**Recommended as supplementary figures.** The per-gene lollipop plots provide "
        "variant-level resolution with domain context as background shading. "
        "Particularly informative for RPL5 (N-terminal clustering within the uL18 domain) "
        "and for genes with multiple domain types.",
        "",
        "### Figure 2 — Track-style (`domain_track_{GENE}`)",
        "**Recommended as supplementary figures or for conference presentations.** "
        "The multi-track layout integrates domain structure, ClinVar density, and local "
        "variant positions in a single vertical panel per gene. Most information-dense "
        "format; requires more reader familiarity with genome-browser-style figures.",
        "",
        "### Suggested framing for manuscript:",
        "",
        '> *"Mapping of local cohort variants onto annotated UniProt domains revealed that '
        "[X]% of variants with resolved amino acid positions co-localize with known "
        "structural or functional domains (Figure N). This proportion is consistent with "
        "[ClinVar concordance statement if applicable], supporting the view that pathogenic "
        'DBA mutations preferentially target functionally constrained regions of the ribosomal protein."*',
        "",
    ]

    # ── Archivos generados ────────────────────────────────────────────────────
    lines += ["## 8. Generated Files", ""]
    for p in figure_paths:
        lines.append(f"- `{p.name}`")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("*DBA Variant Analysis Pipeline — `scripts/domain_analysis.py`*")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    log.info("Reporte escrito: %s", output_path)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    try:
        cfg = _load_config(args.config)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    base = args.config.parent
    resolved = _resolve_paths(cfg, base)
    _setup_logging(args.verbose, resolved["log_domain"])

    log.info("=" * 60)
    log.info("DBA Domain Analysis — script standalone")
    log.info("=" * 60)

    all_genes, prot_len, subunit_map, uniprot_ids = _build_gene_dicts(cfg)
    genes = args.genes if args.genes else all_genes

    unknown = set(genes) - set(all_genes)
    if unknown:
        log.error("Gen(es) desconocidos: %s | Disponibles: %s", sorted(unknown), all_genes)
        sys.exit(1)

    vis_cfg   = cfg.get("visualization", {})
    acmg_colors = vis_cfg.get("acmg_colors", {"P": "#d73027", "LP": "#fc8d59", "VUS": "#fee08b"})
    output_cfg = cfg.get("output", {})
    formats  = output_cfg.get("formats", ["png", "svg"])
    width    = output_cfg.get("width", 1400)
    scale    = output_cfg.get("dpi_scale", 3)
    out_dir  = resolved["out_viz2"]
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Cargar variantes locales ─────────────────────────────────────────
    log.info("─" * 40)
    log.info("Paso 1/6: Cargando variantes locales")
    valid_acmg = list(acmg_colors.keys())
    raw_df  = data_loader.load_variants(resolved["csv"])
    full_df = data_loader.prepare_dataframe(raw_df, valid_acmg)
    full_df = full_df[full_df["Gene"].isin(genes)].copy()
    log.info("Variantes cargadas: %d (genes: %s)", len(full_df), genes)

    # ── 2. Obtener dominios UniProt ─────────────────────────────────────────
    log.info("─" * 40)
    log.info("Paso 2/6: Obteniendo features UniProt")
    features_by_gene = _load_all_domains(genes, uniprot_ids)

    # ── 3. Mapear variantes → dominios ──────────────────────────────────────
    log.info("─" * 40)
    log.info("Paso 3/6: Mapeando variantes a dominios")
    annotated_df = _map_variants_to_domains(full_df, features_by_gene)

    # ── 4. ClinVar (opcional) ───────────────────────────────────────────────
    log.info("─" * 40)
    log.info("Paso 4/6: Cargando ClinVar")
    clinvar_path_df = _load_clinvar_data(resolved, cfg, args.skip_clinvar)
    has_clinvar = clinvar_path_df is not None

    clinvar_positions_by_gene: dict[str, list[int]] = {}
    for gene in genes:
        if has_clinvar:
            clinvar_positions_by_gene[gene] = annotator.get_clinvar_positions(
                gene=gene,
                protein_length=prot_len[gene],
                cfg=cfg,
                clinvar_path_df=clinvar_path_df,
            )
        else:
            clinvar_positions_by_gene[gene] = []

    # ── 5. Estadísticas ─────────────────────────────────────────────────────
    log.info("─" * 40)
    log.info("Paso 5/6: Calculando estadísticas")
    stats = _compute_domain_stats(
        annotated_df, features_by_gene,
        clinvar_positions_by_gene, prot_len, subunit_map,
    )

    for gene, d in stats["per_gene"].items():
        log.info(
            "  %s: %d/%d variantes intra-dominio (%.1f%%)",
            gene, d["n_local_intra"], d["n_local_with_pos"], d["pct_local_intra"],
        )

    # ── 6. Figuras ──────────────────────────────────────────────────────────
    log.info("─" * 40)
    log.info("Paso 6/6: Generando figuras")
    saved_paths: list[Path] = []

    # Orden: 40S primero, 60S después (igual que subunit_comparison)
    ordered_genes = [g for g in genes if subunit_map.get(g) == "40S"] + \
                    [g for g in genes if subunit_map.get(g) == "60S"]

    # --- Figura 1: Arquitectura global ---
    log.info("  Figura 1: domain_architecture_all")
    fig1 = _fig_domain_architecture(
        annotated_df, clinvar_positions_by_gene,
        features_by_gene, prot_len, subunit_map, ordered_genes,
    )
    arch_height = max(600, 140 * len(ordered_genes) + 100)
    saved_paths += visualizer._save_figure(
        fig1, out_dir / "domain_architecture_all",
        formats, width=1500, height=arch_height, scale=scale,
    )

    # --- Figuras 2 y 3: Por gen ---
    for gene in ordered_genes:
        gene_df  = annotated_df[annotated_df["Gene"] == gene].copy()
        features = features_by_gene.get(gene, [])
        cv_pos   = clinvar_positions_by_gene.get(gene, [])

        # Figura 2: Track
        log.info("  Figura 2 track: %s", gene)
        fig2 = _fig_domain_track(gene, gene_df, cv_pos, features, prot_len[gene], acmg_colors)
        n_tracks = 1 + (1 if cv_pos else 0) + (1 if not gene_df.dropna(subset=["position"]).empty else 0)
        saved_paths += visualizer._save_figure(
            fig2, out_dir / f"domain_track_{gene}",
            formats, width=width, height=max(380, 130 * n_tracks + 120), scale=scale,
        )

        # Figura 3: Lollipop con dominios
        log.info("  Figura 3 lollipop: %s", gene)
        fig3 = _fig_domain_lollipop(gene, gene_df, features, prot_len[gene])
        saved_paths += visualizer._save_figure(
            fig3, out_dir / f"domain_lollipop_{gene}",
            formats, width=width, height=480, scale=scale,
        )

    # ── Reporte ─────────────────────────────────────────────────────────────
    report_path = resolved["out_rep"] / "domain_analysis.md"
    _generate_report(stats, features_by_gene, annotated_df, report_path, has_clinvar, saved_paths)

    log.info("=" * 60)
    log.info("Análisis de dominios completado.")
    log.info("Figuras: %d archivos en %s", len(saved_paths), out_dir)
    log.info("Reporte: %s", report_path)
    log.info("=" * 60)


if __name__ == "__main__":
    main()
