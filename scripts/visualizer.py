"""
visualizer.py
=============
Responsabilidad única: selección inteligente de visualizaciones y generación
de gráficos con Plotly (graph_background/plotly.py).

El selector decide qué tipos de gráficos son apropiados para cada gen y para
el conjunto multi-gen, basándose en las características de los datos disponibles.
Todas las decisiones se registran en el log con su justificación.

Funciones exportadas:
    run_gene_visualizations(gene, gene_df, clinvar_positions, hotspots, cfg, out_dir)
    run_global_visualizations(all_data, protein_lengths, subunit_map, cfg, out_dir)
    select_gene_plots(gene, n_local, n_clinvar, hotspots)   -> list[str]
    select_global_plots(all_data)                           -> list[str]
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import gaussian_kde

log = logging.getLogger(__name__)

# Nombre de símbolo plotly equivalente al marcador matplotlib
_MARKER_MAP = {
    "o": "circle",
    "s": "square",
    "D": "diamond",
    "^": "triangle-up",
    "X": "x",
}

_EXPORT_FORMATS = ("png", "svg")  # sobreescrito por config si se especifica


# ---------------------------------------------------------------------------
# Selector inteligente de visualizaciones
# ---------------------------------------------------------------------------

def select_gene_plots(
    gene: str,
    n_local: int,
    n_clinvar: int,
    hotspots: list[dict],
) -> list[str]:
    """
    Decide qué tipos de gráficos generar para un gen, basándose en los datos.

    Criterios:
        - lollipop:        siempre si hay variantes locales con posición (n_local > 0)
        - kde:             si ClinVar aporta ≥ 30 posiciones (suficiente para densidad)
        - histogram:       si ClinVar aporta entre 5 y 29 posiciones
        - hotspot_overlay: si se detectaron hotspots estadísticamente significativos

    Args:
        gene:      Nombre del gen.
        n_local:   Número de variantes locales con posición proteica válida.
        n_clinvar: Número de posiciones ClinVar disponibles.
        hotspots:  Lista de hotspots detectados.

    Returns:
        Lista ordenada de identificadores de gráfico a generar.
    """
    selected: list[str] = []

    if n_local > 0:
        selected.append("lollipop")
        log.info("[VIZ][%s] ✓ lollipop — %d variantes locales con posición", gene, n_local)
    else:
        log.info("[VIZ][%s] ✗ lollipop — sin variantes locales con posición", gene)

    if n_clinvar >= 30:
        selected.append("kde")
        log.info("[VIZ][%s] ✓ KDE — %d posiciones ClinVar (≥30 para estimación continua)", gene, n_clinvar)
    elif n_clinvar >= 5:
        selected.append("histogram")
        log.info("[VIZ][%s] ✓ histogram — %d posiciones ClinVar (5–29, suficiente para frecuencias)", gene, n_clinvar)
    else:
        log.info("[VIZ][%s] ✗ densidad — %d posiciones ClinVar insuficientes (<5)", gene, n_clinvar)

    if hotspots:
        selected.append("hotspot_overlay")
        log.info("[VIZ][%s] ✓ hotspot_overlay — %d hotspot(s) detectado(s)", gene, len(hotspots))

    return selected


def select_global_plots(all_data: dict) -> list[str]:
    """
    Decide qué visualizaciones multi-gen generar.

    Criterios:
        - multi_overview: si ≥ 2 genes tienen variantes locales con posición
        - kde_multi:      si la suma de posiciones ClinVar de todos los genes es ≥ 30
        - heatmap:        si ≥ 3 genes tienen posiciones ClinVar
        - violin:         si ≥ 3 genes tienen posiciones ClinVar

    Args:
        all_data: {gene: {"local_df": df, "clinvar_positions": [], "hotspots": []}}

    Returns:
        Lista de identificadores de gráfico global a generar.
    """
    selected: list[str] = []

    genes_with_local = [
        g for g, d in all_data.items()
        if d.get("local_df") is not None and d["local_df"]["position"].notna().any()
    ]
    genes_with_clinvar = [
        g for g, d in all_data.items()
        if d.get("clinvar_positions")
    ]
    total_clinvar = sum(len(d.get("clinvar_positions", [])) for d in all_data.values())

    if len(genes_with_local) >= 2:
        selected.append("multi_overview")
        log.info("[VIZ][global] ✓ multi_overview — %d genes con variantes locales", len(genes_with_local))

    if total_clinvar >= 30:
        selected.append("kde_multi")
        log.info("[VIZ][global] ✓ KDE multi-gen — %d posiciones ClinVar totales", total_clinvar)

    if len(genes_with_clinvar) >= 3:
        selected.append("heatmap")
        selected.append("violin")
        log.info(
            "[VIZ][global] ✓ heatmap + violin — %d genes con datos ClinVar",
            len(genes_with_clinvar),
        )

    return selected


# ---------------------------------------------------------------------------
# Gráficos por gen
# ---------------------------------------------------------------------------

def _lollipop_fig(
    gene_df: pd.DataFrame,
    gene: str,
    protein_length: int,
    acmg_colors: dict,
    variant_markers: dict,
    hotspots: Optional[list[dict]] = None,
) -> go.Figure:
    """Genera la figura Plotly del lollipop plot para un gen."""
    fig = go.Figure()

    # Línea base de la proteína
    fig.add_shape(
        type="line", x0=0, x1=protein_length, y0=0, y1=0,
        line=dict(color="black", width=2),
    )

    # Hotspots como regiones coloreadas de fondo
    for hs in (hotspots or []):
        fig.add_vrect(
            x0=hs["start"], x1=hs["end"],
            fillcolor="lightyellow", opacity=0.5,
            layer="below", line_width=0,
            annotation_text=f"Hotspot\np={hs['p_value']:.2e}",
            annotation_position="top left",
            annotation_font_size=8,
        )

    # Variantes: stems + heads agrupados por ACMG × variant_type
    for acmg, color in acmg_colors.items():
        for vtype, mpl_marker in variant_markers.items():
            symbol = _MARKER_MAP.get(mpl_marker, "circle")
            group = gene_df[
                (gene_df["ACMG"] == acmg) & (gene_df["variant_type"] == vtype)
            ]
            if group.empty:
                continue

            positions = group["position"].tolist()
            reported_flags = (
                group["Reported"].str.upper().str.strip() == "YES"
                if "Reported" in group.columns
                else [True] * len(group)
            )

            # Stems (líneas verticales con shapes)
            for pos in positions:
                fig.add_shape(
                    type="line", x0=pos, x1=pos, y0=0, y1=1,
                    line=dict(color=color, width=1.5),
                )

            # Heads: separados en reportados (sólido) y no reportados (contorno)
            reported_pos = [p for p, r in zip(positions, reported_flags) if r]
            unreported_pos = [p for p, r in zip(positions, reported_flags) if not r]

            hover_texts = group["AA_change"].fillna("").tolist()

            if reported_pos:
                fig.add_trace(go.Scatter(
                    x=reported_pos, y=[1] * len(reported_pos),
                    mode="markers",
                    marker=dict(symbol=symbol, size=10, color=color,
                                line=dict(color=color, width=1.5)),
                    name=f"{acmg} {vtype} (reportada)",
                    legendgroup=f"{acmg}_{vtype}",
                    hovertemplate="%{text}<extra></extra>",
                    text=[t for t, r in zip(hover_texts, reported_flags) if r],
                ))

            if unreported_pos:
                fig.add_trace(go.Scatter(
                    x=unreported_pos, y=[1] * len(unreported_pos),
                    mode="markers",
                    marker=dict(symbol=f"{symbol}-open", size=10, color=color,
                                line=dict(color=color, width=1.5)),
                    name=f"{acmg} {vtype} (no reportada)",
                    legendgroup=f"{acmg}_{vtype}",
                    showlegend=not reported_pos,  # evita entrada duplicada
                    hovertemplate="%{text}<extra></extra>",
                    text=[t for t, r in zip(hover_texts, reported_flags) if not r],
                ))

    n = len(gene_df)
    fig.update_layout(
        title=dict(text=f"<b>{gene}</b>  (n={n})", font=dict(size=16)),
        xaxis=dict(title="Posición en la proteína (aa)", range=[0, protein_length + 1]),
        yaxis=dict(visible=False, range=[-0.5, 1.8]),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(
            orientation="v", x=1.02, y=1,
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="lightgray", borderwidth=1,
        ),
        margin=dict(r=200),
    )
    fig.update_xaxes(showgrid=True, gridcolor="lightgray", zeroline=False)
    return fig


def _histogram_fig(
    positions: list[int],
    gene: str,
    protein_length: int,
    hotspots: Optional[list[dict]] = None,
) -> go.Figure:
    """Histograma de frecuencias de posiciones ClinVar para un gen."""
    fig = go.Figure()

    for hs in (hotspots or []):
        fig.add_vrect(
            x0=hs["start"], x1=hs["end"],
            fillcolor="lightyellow", opacity=0.5,
            layer="below", line_width=0,
        )

    fig.add_trace(go.Histogram(
        x=positions,
        nbinsx=min(30, protein_length // 4),
        marker_color="#4a90d9",
        marker_line=dict(color="white", width=0.5),
        name="Variantes ClinVar",
    ))

    fig.update_layout(
        title=dict(text=f"<b>{gene}</b> — Distribución de variantes ClinVar patogénicas",
                   font=dict(size=14)),
        xaxis=dict(title="Posición (aa)", range=[0, protein_length + 1]),
        yaxis=dict(title="Frecuencia"),
        bargap=0.1,
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    return fig


def _kde_fig(
    positions: list[int],
    gene: str,
    protein_length: int,
    hotspots: Optional[list[dict]] = None,
) -> go.Figure:
    """KDE de densidad de variantes ClinVar para un gen."""
    if len(positions) < 3:
        return _histogram_fig(positions, gene, protein_length, hotspots)

    x_range = np.linspace(1, protein_length, 300)
    try:
        kde = gaussian_kde(positions, bw_method="scott")
        y_kde = kde(x_range)
    except Exception as exc:
        log.warning("[VIZ][%s] KDE falló (%s) — usando histograma", gene, exc)
        return _histogram_fig(positions, gene, protein_length, hotspots)

    fig = go.Figure()

    for hs in (hotspots or []):
        fig.add_vrect(
            x0=hs["start"], x1=hs["end"],
            fillcolor="lightyellow", opacity=0.5,
            layer="below", line_width=0,
        )

    fig.add_trace(go.Scatter(
        x=x_range, y=y_kde,
        mode="lines",
        line=dict(color="#d73027", width=2),
        fill="tozeroy",
        fillcolor="rgba(215,48,39,0.15)",
        name="Densidad KDE",
    ))

    fig.add_trace(go.Scatter(
        x=positions, y=[0] * len(positions),
        mode="markers",
        marker=dict(symbol="line-ns", color="#333", size=8, line=dict(width=1)),
        name="Posiciones observadas",
        hovertemplate="aa %{x}<extra></extra>",
    ))

    fig.update_layout(
        title=dict(text=f"<b>{gene}</b> — Densidad de variantes patogénicas (KDE)",
                   font=dict(size=14)),
        xaxis=dict(title="Posición (aa)", range=[0, protein_length + 1]),
        yaxis=dict(title="Densidad"),
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    return fig


# ---------------------------------------------------------------------------
# Gráficos globales (multi-gen)
# ---------------------------------------------------------------------------

def _multi_overview_fig(
    all_data: dict,
    protein_lengths: dict,
    subunit_map: dict,
    acmg_colors: dict,
    variant_markers: dict,
) -> go.Figure:
    """Panel multi-gen: lollipop normalizado de todos los genes por subunidad."""
    genes_40s = [g for g in all_data if subunit_map.get(g) == "40S" and all_data[g].get("local_df") is not None]
    genes_60s = [g for g in all_data if subunit_map.get(g) == "60S" and all_data[g].get("local_df") is not None]
    ordered_genes = genes_40s + genes_60s

    if not ordered_genes:
        return go.Figure()

    fig = make_subplots(
        rows=len(ordered_genes), cols=1,
        shared_xaxes=True,
        vertical_spacing=0.02,
        row_titles=ordered_genes,
    )

    for row_idx, gene in enumerate(ordered_genes, start=1):
        gene_df = all_data[gene]["local_df"]
        prot_len = protein_lengths.get(gene, 300)

        gene_df = gene_df.dropna(subset=["position"]).copy()
        if gene_df.empty:
            continue

        norm_pos = gene_df["position"] / prot_len

        # Stems
        for pos in norm_pos:
            fig.add_shape(
                type="line", x0=pos, x1=pos, y0=0, y1=1,
                line=dict(color="lightgray", width=1),
                row=row_idx, col=1,
            )

        colors = gene_df["ACMG"].map(acmg_colors).fillna("gray")
        symbols = gene_df["variant_type"].map(
            {k: _MARKER_MAP.get(v, "circle") for k, v in variant_markers.items()}
        ).fillna("circle")

        fig.add_trace(
            go.Scatter(
                x=norm_pos, y=[1] * len(norm_pos),
                mode="markers",
                marker=dict(
                    color=colors.tolist(),
                    symbol=symbols.tolist(),
                    size=9,
                    line=dict(width=1, color="gray"),
                ),
                text=gene_df["AA_change"].fillna("").tolist(),
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
                name=gene,
            ),
            row=row_idx, col=1,
        )

        fig.update_yaxes(visible=False, row=row_idx, col=1)

    # Línea separadora 40S / 60S
    if genes_40s and genes_60s:
        sep_y = 1 - len(genes_40s) / len(ordered_genes)
        fig.add_hline(y=sep_y, line_dash="dash", line_color="steelblue",
                      annotation_text="40S ↑  60S ↓", annotation_position="top right")

    fig.update_layout(
        title=dict(text="<b>Variantes DBA — visión global (posición normalizada)</b>",
                   font=dict(size=15)),
        xaxis=dict(title="Posición relativa en la proteína (0–1)", range=[-0.02, 1.05]),
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=120 * len(ordered_genes) + 100,
    )
    return fig


def _kde_multi_fig(
    all_data: dict,
    protein_lengths: dict,
) -> go.Figure:
    """KDE normalizado superpuesto para todos los genes con datos ClinVar."""
    fig = go.Figure()
    palette = ["#d73027", "#fc8d59", "#fee08b", "#91bfdb", "#4575b4", "#1a9850", "#762a83"]

    for i, (gene, data) in enumerate(all_data.items()):
        positions = data.get("clinvar_positions", [])
        if len(positions) < 3:
            continue
        prot_len = protein_lengths.get(gene, 300)
        norm = [p / prot_len for p in positions]
        x_range = np.linspace(0, 1, 300)
        try:
            kde = gaussian_kde(norm, bw_method="scott")
            y_kde = kde(x_range)
        except Exception:
            continue

        color = palette[i % len(palette)]
        fig.add_trace(go.Scatter(
            x=x_range, y=y_kde,
            mode="lines",
            name=gene,
            line=dict(color=color, width=2),
        ))

    fig.update_layout(
        title=dict(text="<b>Densidad de variantes patogénicas ClinVar — todos los genes</b>",
                   font=dict(size=14)),
        xaxis=dict(title="Posición normalizada en la proteína (0 = N-term, 1 = C-term)"),
        yaxis=dict(title="Densidad"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(title="Gen"),
    )
    return fig


def _heatmap_fig(
    all_data: dict,
    protein_lengths: dict,
    n_bins: int = 50,
) -> go.Figure:
    """Heatmap de densidad de mutaciones: genes × posición normalizada."""
    genes = [g for g, d in all_data.items() if d.get("clinvar_positions")]
    if not genes:
        return go.Figure()

    z = np.zeros((len(genes), n_bins))
    x_labels = np.round(np.linspace(0, 1, n_bins), 2).tolist()

    for i, gene in enumerate(genes):
        positions = all_data[gene]["clinvar_positions"]
        prot_len = protein_lengths.get(gene, 300)
        norm = [p / prot_len for p in positions]
        counts, _ = np.histogram(norm, bins=n_bins, range=(0.0, 1.0))
        z[i] = counts

    fig = go.Figure(go.Heatmap(
        z=z,
        x=x_labels,
        y=genes,
        colorscale="Magma",
        colorbar=dict(title="Frecuencia"),
        hovertemplate="Gen: %{y}<br>Posición rel: %{x}<br>Conteo: %{z}<extra></extra>",
    ))

    fig.update_layout(
        title=dict(text="<b>Mapa de densidad de mutaciones — genes DBA</b>",
                   font=dict(size=14)),
        xaxis=dict(title="Posición normalizada (0–1)"),
        yaxis=dict(title="Gen"),
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    return fig


def _violin_fig(
    all_data: dict,
    protein_lengths: dict,
) -> go.Figure:
    """Violin plot de distribución de variantes ClinVar (posición normalizada)."""
    fig = go.Figure()
    palette = ["#d73027", "#fc8d59", "#fee08b", "#91bfdb", "#4575b4", "#1a9850", "#762a83"]

    for i, (gene, data) in enumerate(all_data.items()):
        positions = data.get("clinvar_positions", [])
        if len(positions) < 5:
            continue
        prot_len = protein_lengths.get(gene, 300)
        norm = [p / prot_len for p in positions]

        fig.add_trace(go.Violin(
            y=norm,
            name=gene,
            box_visible=True,
            meanline_visible=True,
            fillcolor=palette[i % len(palette)],
            line_color="black",
            opacity=0.7,
            hoverinfo="y+name",
        ))

    fig.update_layout(
        title=dict(text="<b>Distribución de variantes patogénicas — comparación entre genes</b>",
                   font=dict(size=14)),
        yaxis=dict(title="Posición normalizada (0–1)"),
        xaxis=dict(title="Gen"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        violingap=0.2,
        violinmode="overlay",
    )
    return fig


# ---------------------------------------------------------------------------
# Exportación
# ---------------------------------------------------------------------------

def _save_figure(
    fig: go.Figure,
    output_path: Path,
    formats: list[str],
    width: int = 1400,
    height: int = 700,
    scale: int = 3,
) -> list[Path]:
    """
    Exporta una figura Plotly en los formatos especificados.

    Args:
        fig:         Figura Plotly.
        output_path: Ruta base sin extensión (p. ej. .../RPS19_lollipop).
        formats:     Lista de formatos: "png", "svg", "pdf", "html".
        width:       Ancho en píxeles (para PNG/SVG/PDF).
        height:      Alto en píxeles.
        scale:       Factor de escala para PNG (3 ≈ 300 DPI).

    Returns:
        Lista de rutas de archivos guardados.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []

    for fmt in formats:
        dest = output_path.with_suffix(f".{fmt}")
        try:
            if fmt == "html":
                fig.write_html(str(dest), include_plotlyjs="cdn")
            else:
                fig.write_image(str(dest), width=width, height=height,
                                scale=scale, format=fmt)
            log.info("Guardado: %s", dest)
            saved.append(dest)
        except Exception as exc:
            log.warning("No se pudo exportar %s: %s", dest, exc)
            # Fallback a HTML si falta kaleido
            if fmt in ("png", "svg", "pdf"):
                html_dest = output_path.with_suffix(".html")
                try:
                    fig.write_html(str(html_dest), include_plotlyjs="cdn")
                    log.info("Fallback HTML guardado: %s (instala kaleido para %s)", html_dest, fmt)
                    saved.append(html_dest)
                except Exception:
                    pass

    return saved


# ---------------------------------------------------------------------------
# Puntos de entrada públicos
# ---------------------------------------------------------------------------

def run_gene_visualizations(
    gene: str,
    gene_df: pd.DataFrame,
    clinvar_positions: list[int],
    hotspots: list[dict],
    protein_length: int,
    cfg: dict,
    out_dir: Path,
) -> list[Path]:
    """
    Selecciona y genera todas las visualizaciones apropiadas para un gen.

    Args:
        gene:              Nombre del gen.
        gene_df:           DataFrame de variantes locales del gen (ya con 'position').
        clinvar_positions: Lista de posiciones ClinVar patogénicas.
        hotspots:          Lista de hotspots detectados.
        protein_length:    Longitud de la proteína.
        cfg:               Configuración completa.
        out_dir:           Directorio de salida para las figuras.

    Returns:
        Lista de rutas de archivos generados.
    """
    vis = cfg.get("visualization", {})
    acmg_colors   = vis.get("acmg_colors", {})
    variant_markers = vis.get("variant_markers", {})
    output_cfg = cfg.get("output", {})
    formats = output_cfg.get("formats", list(_EXPORT_FORMATS))
    width   = output_cfg.get("width", 1400)
    height  = output_cfg.get("height", 700)
    scale   = output_cfg.get("dpi_scale", 3)

    n_local = int(gene_df["position"].notna().sum()) if not gene_df.empty else 0
    plot_types = select_gene_plots(gene, n_local, len(clinvar_positions), hotspots)

    saved_paths: list[Path] = []

    for ptype in plot_types:
        if ptype == "lollipop":
            fig = _lollipop_fig(gene_df.dropna(subset=["position"]), gene,
                                protein_length, acmg_colors, variant_markers, hotspots)
            saved_paths += _save_figure(fig, out_dir / f"{gene}_lollipop",
                                        formats, width, height, scale)

        elif ptype == "kde":
            fig = _kde_fig(clinvar_positions, gene, protein_length, hotspots)
            saved_paths += _save_figure(fig, out_dir / f"{gene}_kde",
                                        formats, width, height, scale)

        elif ptype == "histogram":
            fig = _histogram_fig(clinvar_positions, gene, protein_length, hotspots)
            saved_paths += _save_figure(fig, out_dir / f"{gene}_histogram",
                                        formats, width, height, scale)

        elif ptype == "hotspot_overlay":
            # Si ya se generó lollipop o KDE con overlay, no duplicar
            log.debug("[VIZ][%s] hotspot_overlay ya incluido en lollipop/KDE", gene)

    return saved_paths


def run_global_visualizations(
    all_data: dict,
    protein_lengths: dict,
    subunit_map: dict,
    cfg: dict,
    out_dir: Path,
) -> list[Path]:
    """
    Selecciona y genera todas las visualizaciones multi-gen.

    Args:
        all_data:        {gene: {"local_df": df, "clinvar_positions": [], "hotspots": []}}
        protein_lengths: {gene: int}
        subunit_map:     {gene: "40S"|"60S"}
        cfg:             Configuración completa.
        out_dir:         Directorio de salida.

    Returns:
        Lista de rutas de archivos generados.
    """
    vis = cfg.get("visualization", {})
    acmg_colors   = vis.get("acmg_colors", {})
    variant_markers = vis.get("variant_markers", {})
    output_cfg = cfg.get("output", {})
    formats = output_cfg.get("formats", list(_EXPORT_FORMATS))
    width   = output_cfg.get("width", 1400)
    height  = output_cfg.get("height", 700)
    scale   = output_cfg.get("dpi_scale", 3)

    plot_types = select_global_plots(all_data)
    saved_paths: list[Path] = []

    for ptype in plot_types:
        if ptype == "multi_overview":
            fig = _multi_overview_fig(all_data, protein_lengths, subunit_map,
                                      acmg_colors, variant_markers)
            saved_paths += _save_figure(
                fig, out_dir / "global_multi_overview",
                formats, width, max(700, 120 * len(all_data) + 100), scale,
            )

        elif ptype == "kde_multi":
            fig = _kde_multi_fig(all_data, protein_lengths)
            saved_paths += _save_figure(fig, out_dir / "global_kde_multi",
                                        formats, width, height, scale)

        elif ptype == "heatmap":
            fig = _heatmap_fig(all_data, protein_lengths)
            saved_paths += _save_figure(fig, out_dir / "global_heatmap",
                                        formats, width, 500, scale)

        elif ptype == "violin":
            fig = _violin_fig(all_data, protein_lengths)
            saved_paths += _save_figure(fig, out_dir / "global_violin",
                                        formats, width, height, scale)

    return saved_paths


# =============================================================================
# SECCIÓN V2 — Visualizaciones alternativas orientadas a publicación científica
# =============================================================================
#
# Justificación general: los lollipop plots muestran solo la posición y ACMG
# de cada variante. Las visualizaciones v2 explotan dimensiones adicionales
# de los datos (tipo de variante, severidad relativa, densidad poblacional
# ClinVar, contexto funcional UniProt, comparación entre genes/subunidades)
# que el lollipop no puede representar sin saturación visual.
#
# Figuras por gen:
#   {GENE}_acmg_sunburst    — composición ACMG × tipo de variante (sunburst)
#   {GENE}_position_scatter — posición vs severidad ACMG (scatter 3-dim)
#   {GENE}_clinvar_overlay  — densidad ClinVar + cohorte local superpuesta
#   {GENE}_track            — mapa de anotación integrado (genomic-track style)
#
# Figuras globales:
#   global_acmg_comparison  — barras agrupadas ACMG × gen
#   global_variant_spectrum — heatmap espectro mutacional
#   global_subunit_comparison — violin 40S vs 60S
#   global_data_overview    — bubble chart panorama de datos
# =============================================================================

# Constantes v2 ---------------------------------------------------------------

_VARIANT_TYPE_COLORS: dict[str, str] = {
    "Missense":   "#2196F3",
    "Nonsense":   "#F44336",
    "Frameshift": "#FF9800",
    "Splice":     "#9C27B0",
    "Other":      "#607D8B",
}

_ACMG_SEVERITY: dict[str, int] = {"P": 3, "LP": 2, "VUS": 1}

_FEATURE_TYPE_COLORS: dict[str, str] = {
    "Domain":       "#1565C0",
    "Region":       "#2E7D32",
    "Motif":        "#7B1FA2",
    "Zinc finger":  "#AD1457",
    "Binding site": "#E65100",
    "Active site":  "#BF360C",
    "Helix":        "#A5D6A7",
    "Beta strand":  "#90CAF9",
    "Turn":         "#CE93D8",
}
_FEATURE_TYPE_DEFAULT = "#78909C"


# ---------------------------------------------------------------------------
# Figuras v2 por gen
# ---------------------------------------------------------------------------

def _acmg_sunburst_fig(
    gene_df: pd.DataFrame,
    gene: str,
    acmg_colors: dict,
) -> go.Figure:
    """
    Sunburst de dos anillos: ACMG (interior) × tipo de variante (exterior).

    Revela: la composición de patogenicidad Y el espectro mutacional en una
    sola figura. El área de cada segmento es proporcional al número de
    variantes — imposible de transmitir en un lollipop sin saturación visual.
    """
    acmg_order = list(acmg_colors.keys())
    variant_types = list(_VARIANT_TYPE_COLORS.keys())
    root_id = f"{gene}_root"
    ids, labels, parents, values, colors = [], [], [], [], []

    ids.append(root_id)
    labels.append(gene)
    parents.append("")
    values.append(0)
    colors.append("white")

    for acmg in acmg_order:
        total = int((gene_df["ACMG"] == acmg).sum())
        if total == 0:
            continue
        acmg_id = f"{gene}_{acmg}"
        ids.append(acmg_id)
        labels.append(f"<b>{acmg}</b><br>n={total}")
        parents.append(root_id)
        values.append(total)
        colors.append(acmg_colors.get(acmg, "#aaa"))

        for vtype in variant_types:
            n = int(
                ((gene_df["ACMG"] == acmg) & (gene_df["variant_type"] == vtype)).sum()
            )
            if n == 0:
                continue
            ids.append(f"{gene}_{acmg}_{vtype}")
            labels.append(f"{vtype}<br>n={n}")
            parents.append(acmg_id)
            values.append(n)
            colors.append(_VARIANT_TYPE_COLORS.get(vtype, "#aaa"))

    if len(ids) <= 1:
        return go.Figure()

    fig = go.Figure(go.Sunburst(
        ids=ids,
        labels=labels,
        parents=parents,
        values=values,
        branchvalues="total",
        marker=dict(colors=colors, line=dict(color="white", width=1.5)),
        hovertemplate=(
            "<b>%{label}</b><br>"
            "%{value} variantes (%{percentParent:.0%} del grupo)"
            "<extra></extra>"
        ),
        textfont=dict(size=11),
    ))
    fig.update_layout(
        title=dict(
            text=(
                f"<b>{gene}</b> — Clasificacion ACMG x Espectro mutacional"
                f" (n={len(gene_df)})"
            ),
            font=dict(size=14),
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(t=60, l=10, r=10, b=10),
    )
    return fig


def _position_scatter_fig(
    gene_df: pd.DataFrame,
    gene: str,
    protein_length: int,
    hotspots: Optional[list[dict]],
    acmg_colors: dict,
) -> go.Figure:
    """
    Scatter: posición en la proteína (X) × severidad ACMG ordinal (Y).
    Color = tipo de variante; símbolo abierto = no publicada en literatura.

    Revela: concentración posicional de variantes clasificadas por severidad.
    Codifica 3 dimensiones (posición, severidad, tipo) sin saturación visual.
    Ventaja vs lollipop: el lollipop colapsa toda la severidad en el color;
    aquí la posición en Y separa P/LP/VUS permitiendo leer los tres ejes.
    """
    df = gene_df.dropna(subset=["position"]).copy()
    if df.empty:
        return go.Figure()

    df["severity"] = df["ACMG"].map(_ACMG_SEVERITY).fillna(0).astype(float)
    rng = np.random.default_rng(42)
    df["y_j"] = df["severity"] + rng.uniform(-0.18, 0.18, size=len(df))

    if "Reported" in df.columns:
        df["is_reported"] = df["Reported"].str.upper().str.strip() == "YES"
    else:
        df["is_reported"] = True

    fig = go.Figure()
    fig.add_shape(
        type="line", x0=1, x1=protein_length, y0=0.4, y1=0.4,
        line=dict(color="#BDBDBD", width=3),
    )
    for hs in (hotspots or []):
        fig.add_vrect(
            x0=hs["start"], x1=hs["end"],
            fillcolor="rgba(249,168,37,0.18)", layer="below", line_width=0,
        )

    hover_tmpl = (
        "<b>%{text}</b><br>"
        "Posicion: %{x} aa<br>"
        "ACMG: %{customdata[0]}<br>"
        "Tipo: %{customdata[1]}<extra></extra>"
    )
    for vtype, color in _VARIANT_TYPE_COLORS.items():
        sub = df[df["variant_type"] == vtype]
        if sub.empty:
            continue
        for mask, symbol, show_legend, suffix in [
            (sub["is_reported"],  "circle",      True,  ""),
            (~sub["is_reported"], "circle-open", False, " (no pub.)"),
        ]:
            grp = sub[mask]
            if grp.empty:
                continue
            fig.add_trace(go.Scatter(
                x=grp["position"].tolist(),
                y=grp["y_j"].tolist(),
                mode="markers",
                name=f"{vtype}{suffix}",
                legendgroup=vtype,
                showlegend=show_legend,
                marker=dict(
                    color=color, symbol=symbol, size=12,
                    line=dict(color=color, width=1.5),
                ),
                text=grp["AA_change"].fillna("").tolist(),
                customdata=grp[["ACMG", "variant_type"]].values.tolist(),
                hovertemplate=hover_tmpl,
            ))

    fig.update_layout(
        title=dict(
            text=f"<b>{gene}</b> — Distribucion posicional por severidad ACMG",
            font=dict(size=14),
        ),
        xaxis=dict(
            title="Posicion en la proteina (aa)",
            range=[0, protein_length + 2],
            showgrid=True, gridcolor="#F5F5F5",
        ),
        yaxis=dict(
            title="Severidad ACMG",
            tickvals=[1, 2, 3],
            ticktext=["VUS", "LP", "P"],
            range=[0.2, 3.8],
            gridcolor="#F0F0F0",
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(title="Tipo de variante"),
    )
    return fig


def _clinvar_local_overlay_fig(
    gene_df: pd.DataFrame,
    clinvar_positions: list[int],
    gene: str,
    protein_length: int,
    acmg_colors: dict,
    hotspots: Optional[list[dict]],
) -> go.Figure:
    """
    KDE de variantes ClinVar patogénicas (área sombreada) con variantes locales
    superpuestas como marcas rug coloreadas por ACMG.

    Revela: correspondencia o divergencia entre el patrón poblacional (ClinVar)
    y la cohorte estudiada. Regiones con alta densidad ClinVar pero sin
    variantes locales pueden reflejar sesgo de diagnóstico o diferencias
    genéticas de la población.
    Ventaja vs lollipop: integra epidemiología global y datos propios en el
    mismo eje posicional sin necesidad de comparar dos figuras separadas.
    """
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    for hs in (hotspots or []):
        fig.add_vrect(
            x0=hs["start"], x1=hs["end"],
            fillcolor="rgba(249,168,37,0.18)", layer="below", line_width=0,
        )

    if len(clinvar_positions) >= 3:
        x_range = np.linspace(1, protein_length, 400)
        try:
            kde = gaussian_kde(clinvar_positions, bw_method="scott")
            y_kde = kde(x_range)
            fig.add_trace(go.Scatter(
                x=x_range, y=y_kde,
                mode="lines",
                fill="tozeroy",
                line=dict(color="rgba(183,28,28,0.85)", width=2),
                fillcolor="rgba(183,28,28,0.10)",
                name=f"ClinVar patogenicas (n={len(clinvar_positions)})",
            ), secondary_y=False)
        except Exception as exc:
            log.warning("[VIZ2][%s] KDE overlay fallo: %s", gene, exc)

    df_local = gene_df.dropna(subset=["position"]).copy()
    for acmg, color in acmg_colors.items():
        sub = df_local[df_local["ACMG"] == acmg]
        if sub.empty:
            continue
        fig.add_trace(go.Scatter(
            x=sub["position"].tolist(),
            y=[0.0] * len(sub),
            mode="markers",
            marker=dict(
                symbol="line-ns", size=14, color=color,
                line=dict(width=2.5, color=color),
            ),
            name=f"Local {acmg} (n={len(sub)})",
            text=sub["AA_change"].fillna("").tolist(),
            hovertemplate="<b>%{text}</b><br>pos. %{x} aa<extra></extra>",
        ), secondary_y=True)

    fig.update_xaxes(title_text="Posicion en la proteina (aa)", range=[0, protein_length + 1])
    fig.update_yaxes(title_text="Densidad KDE (ClinVar)", secondary_y=False, rangemode="tozero")
    fig.update_yaxes(visible=False, secondary_y=True)
    fig.update_layout(
        title=dict(
            text=f"<b>{gene}</b> — Densidad ClinVar + variantes de la cohorte local",
            font=dict(size=14),
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(orientation="h", y=-0.18, xanchor="center", x=0.5),
    )
    return fig


def _track_fig(
    gene_df: pd.DataFrame,
    clinvar_positions: list[int],
    hotspots: Optional[list[dict]],
    domains: list[dict],
    protein_length: int,
    gene: str,
    acmg_colors: dict,
) -> go.Figure:
    """
    Figura tipo 'genomic track': pistas apiladas con eje X compartido.

    Pistas incluidas según datos disponibles:
        1. Dominios UniProt — estructura y función de la proteína
        2. Densidad ClinVar — frecuencia de variantes por región (bar chart)
        3. Variantes locales — rug marks coloreados por ACMG
        4. Hotspots — regiones rectangulares con p-valor anotado

    Revela: contexto funcional de cada variante. Identifica si las variantes
    caen en dominios conocidos, si los hotspots solapan con regiones
    funcionalmente relevantes, y cómo se compara la cohorte local con ClinVar.
    Ventaja vs lollipop: única figura que integra estructura proteica,
    epidemiología (ClinVar), datos propios e inferencia estadística en un
    solo panel con eje posicional compartido.
    """
    has_domains  = bool(domains)
    has_clinvar  = bool(clinvar_positions)
    has_hotspots = bool(hotspots)
    df_local = gene_df.dropna(subset=["position"]).copy()
    has_local = not df_local.empty

    track_names = []
    if has_domains:
        track_names.append("domains")
    if has_clinvar:
        track_names.append("clinvar")
    if has_local:
        track_names.append("local")
    if has_hotspots:
        track_names.append("hotspots")

    n_rows = len(track_names)
    if n_rows == 0:
        return go.Figure()

    _rh = {"domains": 0.15, "clinvar": 0.45, "local": 0.25, "hotspots": 0.15}
    raw_h = [_rh[t] for t in track_names]
    total_h = sum(raw_h)
    row_heights = [h / total_h for h in raw_h]

    fig = make_subplots(
        rows=n_rows, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.05,
        row_heights=row_heights,
        subplot_titles=[t.capitalize() for t in track_names],
    )

    for row_idx, track in enumerate(track_names, start=1):
        if track == "domains":
            fig.add_shape(
                type="rect", x0=1, x1=protein_length, y0=0.3, y1=0.7,
                fillcolor="#EEEEEE", line_color="#BDBDBD", line_width=1,
                row=row_idx, col=1,
            )
            for feat in domains:
                fcolor = _FEATURE_TYPE_COLORS.get(feat["type"], _FEATURE_TYPE_DEFAULT)
                fig.add_shape(
                    type="rect",
                    x0=feat["start"], x1=feat["end"], y0=0.1, y1=0.9,
                    fillcolor=fcolor, opacity=0.75,
                    line_color="white", line_width=1,
                    row=row_idx, col=1,
                )
                mid = (feat["start"] + feat["end"]) / 2.0
                desc = feat.get("description", feat.get("type", ""))[:28]
                fig.add_trace(go.Scatter(
                    x=[mid], y=[0.5],
                    mode="markers",
                    marker=dict(color=fcolor, size=1, opacity=0),
                    text=[f"{feat['type']}: {desc}<br>aa {feat['start']}-{feat['end']}"],
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False, name="",
                ), row=row_idx, col=1)
            fig.update_yaxes(visible=False, range=[0, 1], row=row_idx, col=1)

        elif track == "clinvar":
            n_bins = max(10, min(40, protein_length // 5))
            counts, bin_edges = np.histogram(
                clinvar_positions, bins=n_bins, range=(1, protein_length)
            )
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
            bar_w = (protein_length / n_bins) * 0.85
            fig.add_trace(go.Bar(
                x=bin_centers, y=counts,
                marker_color="#C62828", opacity=0.80,
                name=f"ClinVar (n={len(clinvar_positions)})",
                width=bar_w,
                hovertemplate="aa %{x:.0f}<br>n=%{y}<extra></extra>",
            ), row=row_idx, col=1)
            fig.update_yaxes(title_text="n", row=row_idx, col=1)

        elif track == "local":
            for acmg, color in acmg_colors.items():
                sub = df_local[df_local["ACMG"] == acmg]
                if sub.empty:
                    continue
                fig.add_trace(go.Scatter(
                    x=sub["position"].tolist(),
                    y=[0.5] * len(sub),
                    mode="markers",
                    marker=dict(
                        symbol="line-ns", size=16, color=color,
                        line=dict(width=2.5, color=color),
                    ),
                    name=f"Local {acmg}",
                    text=sub["AA_change"].fillna("").tolist(),
                    hovertemplate="<b>%{text}</b><extra></extra>",
                ), row=row_idx, col=1)
            fig.update_yaxes(visible=False, range=[0, 1], row=row_idx, col=1)

        elif track == "hotspots":
            fig.add_shape(
                type="line", x0=1, x1=protein_length, y0=0.5, y1=0.5,
                line=dict(color="#BDBDBD", width=1),
                row=row_idx, col=1,
            )
            for hs in hotspots:
                mid = (hs["start"] + hs["end"]) / 2.0
                fig.add_shape(
                    type="rect",
                    x0=hs["start"], x1=hs["end"], y0=0.2, y1=0.8,
                    fillcolor="#F9A825", opacity=0.85,
                    line_color="#F57F17", line_width=1,
                    row=row_idx, col=1,
                )
                fig.add_trace(go.Scatter(
                    x=[mid], y=[0.5],
                    mode="markers",
                    marker=dict(color="#F9A825", size=1, opacity=0),
                    text=[
                        f"Hotspot aa {hs['start']}-{hs['end']}<br>"
                        f"n={hs['count']}<br>p={hs['p_value']:.2e}"
                    ],
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False, name="",
                ), row=row_idx, col=1)
            fig.update_yaxes(visible=False, range=[0, 1], row=row_idx, col=1)

        fig.update_xaxes(range=[0, protein_length + 1], row=row_idx, col=1)
        if row_idx < n_rows:
            fig.update_xaxes(showticklabels=False, row=row_idx, col=1)
        else:
            fig.update_xaxes(
                title_text="Posicion en la proteina (aa)", row=row_idx, col=1
            )

    fig.update_layout(
        title=dict(
            text=f"<b>{gene}</b> — Mapa de anotacion proteica integrado",
            font=dict(size=14),
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=max(450, 130 * n_rows + 120),
        showlegend=True,
        legend=dict(orientation="h", y=-0.12, xanchor="center", x=0.5),
        bargap=0,
    )
    return fig


# ---------------------------------------------------------------------------
# Figuras v2 globales (multi-gen)
# ---------------------------------------------------------------------------

def _acmg_comparison_bar(
    all_data: dict,
    acmg_colors: dict,
) -> go.Figure:
    """
    Barras agrupadas: genes (eje X) × clasificación ACMG (grupos de barras).

    Revela: carga total de variantes y composición de patogenicidad por gen.
    Permite comparar directamente qué genes tienen mayor proporción de P vs
    LP vs VUS — respuesta imposible mirando lollipops individuales.
    """
    genes_with_local = [
        g for g, d in all_data.items()
        if d.get("local_df") is not None and not d["local_df"].empty
    ]
    if not genes_with_local:
        return go.Figure()

    totals = {g: len(all_data[g]["local_df"]) for g in genes_with_local}
    genes_sorted = sorted(genes_with_local, key=lambda g: totals[g], reverse=True)

    fig = go.Figure()
    for acmg, color in acmg_colors.items():
        counts = [
            int((all_data[g]["local_df"]["ACMG"] == acmg).sum())
            for g in genes_sorted
        ]
        fig.add_trace(go.Bar(
            name=acmg,
            x=genes_sorted,
            y=counts,
            marker_color=color,
            text=counts,
            textposition="outside",
            hovertemplate=f"<b>{acmg}</b><br>%{{x}}: %{{y}} variantes<extra></extra>",
        ))

    fig.update_layout(
        barmode="group",
        title=dict(
            text="<b>Carga de variantes por gen — Clasificacion ACMG</b>",
            font=dict(size=14),
        ),
        xaxis=dict(title="Gen", showgrid=False),
        yaxis=dict(title="Numero de variantes", showgrid=True, gridcolor="#F0F0F0"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(title="ACMG"),
    )
    return fig


def _variant_spectrum_heatmap_v2(all_data: dict) -> go.Figure:
    """
    Heatmap: genes (filas) × tipo de variante (columnas), valor = proporción.
    Las celdas muestran n y % simultáneamente.

    Revela: si existe una 'firma mutacional' por gen (ej. RPS19 con alta
    proporción de Nonsense vs RPL5 dominado por Missense). La comparación
    del espectro entre genes requiere ver todos a la vez.
    """
    variant_types = list(_VARIANT_TYPE_COLORS.keys())
    genes_with_local = [
        g for g, d in all_data.items()
        if d.get("local_df") is not None and not d["local_df"].empty
    ]
    if len(genes_with_local) < 2:
        return go.Figure()

    z_raw, text_cells = [], []
    for gene in genes_with_local:
        df = all_data[gene]["local_df"]
        counts = [int((df["variant_type"] == vt).sum()) for vt in variant_types]
        total = sum(counts) or 1
        pcts = [c / total for c in counts]
        z_raw.append(pcts)
        text_cells.append([
            f"n={counts[j]}<br>{pcts[j]*100:.0f}%"
            for j in range(len(variant_types))
        ])

    fig = go.Figure(go.Heatmap(
        z=z_raw,
        x=variant_types,
        y=genes_with_local,
        text=text_cells,
        texttemplate="%{text}",
        textfont=dict(size=11),
        colorscale="Blues",
        colorbar=dict(title="Proporcion", tickformat=".0%"),
        zmin=0, zmax=1,
        hovertemplate="<b>%{y} — %{x}</b><br>%{text}<extra></extra>",
    ))
    fig.update_layout(
        title=dict(
            text="<b>Espectro mutacional por gen — Proporcion de cada tipo de variante</b>",
            font=dict(size=14),
        ),
        xaxis=dict(title="Tipo de variante"),
        yaxis=dict(title="Gen"),
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    return fig


def _subunit_comparison_fig(
    all_data: dict,
    protein_lengths: dict,
    subunit_map: dict,
) -> go.Figure:
    """
    Violin plot (con puntos individuales jitter) de posición ClinVar normalizada,
    agrupado por subunidad ribosomal (40S vs 60S).

    Revela: si las subunidades 40S y 60S muestran patrones diferentes de
    concentración mutacional (N-terminal vs C-terminal, unimodal vs
    multimodal). Hipótesis biológica: diferentes restricciones evolutivas
    según el rol funcional de cada subunidad en el ribosoma.
    """
    _colors = {"40S": "#1565C0", "60S": "#B71C1C"}
    fig = go.Figure()

    for subunit in ["40S", "60S"]:
        genes_in_subunit = [g for g, s in subunit_map.items() if s == subunit]
        all_norm: list[float] = []
        for gene in genes_in_subunit:
            positions = all_data.get(gene, {}).get("clinvar_positions", [])
            if not positions:
                continue
            prot_len = protein_lengths.get(gene, 300)
            all_norm.extend([p / prot_len for p in positions])

        if not all_norm:
            continue

        fig.add_trace(go.Violin(
            y=all_norm,
            x=[subunit] * len(all_norm),
            name=subunit,
            box_visible=True,
            meanline_visible=True,
            points="all",
            jitter=0.35,
            pointpos=0,
            fillcolor=_colors.get(subunit, "#607D8B"),
            line_color="black",
            marker=dict(color=_colors.get(subunit, "#607D8B"), size=4, opacity=0.5),
            opacity=0.55,
            hovertemplate=(
                f"<b>{subunit}</b><br>pos. normalizada: %{{y:.3f}}<extra></extra>"
            ),
        ))

    if not fig.data:
        return go.Figure()

    fig.update_layout(
        title=dict(
            text="<b>Distribucion de variantes ClinVar — Subunidades 40S vs 60S</b>",
            font=dict(size=14),
        ),
        yaxis=dict(
            title="Posicion normalizada (0 = N-term, 1 = C-term)",
            range=[-0.05, 1.05],
        ),
        xaxis=dict(title="Subunidad ribosomal"),
        violinmode="group",
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    return fig


def _data_overview_bubble(
    all_data: dict,
    protein_lengths: dict,
    subunit_map: dict,
) -> go.Figure:
    """
    Bubble chart: longitud proteína (X) × variantes locales (Y),
    tamaño burbuja = n ClinVar patogénicas, color = subunidad.

    Revela: panorama completo de disponibilidad de datos por gen. Identifica
    genes con alta cobertura ClinVar pero pocos datos locales (priorización
    diagnóstica) o viceversa (genes recientemente implicados en DBA).
    """
    _colors = {"40S": "#1565C0", "60S": "#B71C1C", "Unknown": "#607D8B"}
    fig = go.Figure()

    for subunit in ["40S", "60S"]:
        gene_list = [g for g in all_data if subunit_map.get(g) == subunit]
        if not gene_list:
            continue

        x_prot  = [protein_lengths.get(g, 0) for g in gene_list]
        df_list = [all_data[g].get("local_df") for g in gene_list]
        y_local = [
            int(d["position"].notna().sum()) if d is not None and not d.empty else 0
            for d in df_list
        ]
        n_clinvar  = [len(all_data[g].get("clinvar_positions", [])) for g in gene_list]
        n_hotspots = [len(all_data[g].get("hotspots", [])) for g in gene_list]

        fig.add_trace(go.Scatter(
            x=x_prot,
            y=y_local,
            mode="markers+text",
            name=f"Subunidad {subunit}",
            marker=dict(
                size=[max(16, nc * 0.25) for nc in n_clinvar],
                color=_colors[subunit],
                opacity=0.72,
                line=dict(color="white", width=1.5),
                sizemode="diameter",
            ),
            text=gene_list,
            textposition="top center",
            customdata=list(zip(n_clinvar, n_hotspots)),
            hovertemplate=(
                "<b>%{text}</b><br>"
                "Longitud proteina: %{x} aa<br>"
                "Variantes locales: %{y}<br>"
                "ClinVar patogenicas: %{customdata[0]}<br>"
                "Hotspots detectados: %{customdata[1]}<extra></extra>"
            ),
        ))

    fig.update_layout(
        title=dict(
            text=(
                "<b>Panorama de datos por gen DBA</b><br>"
                "<sub>X: longitud proteina | Y: variantes locales | "
                "Tamano burbuja: n ClinVar patogenicas</sub>"
            ),
            font=dict(size=14),
        ),
        xaxis=dict(title="Longitud de la proteina (aa)"),
        yaxis=dict(title="Variantes locales con posicion"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(title="Subunidad"),
    )
    return fig


# ---------------------------------------------------------------------------
# Selectores inteligentes v2
# ---------------------------------------------------------------------------

def select_gene_plots_v2(
    gene: str,
    n_local: int,
    n_clinvar: int,
    hotspots: list[dict],
    domains: list[dict],
) -> list[str]:
    """
    Decide qué visualizaciones v2 generar para un gen dado sus datos.

    Criterios:
        - acmg_sunburst:         n_local >= 2
        - position_scatter:      n_local >= 2
        - clinvar_local_overlay: n_local > 0 y n_clinvar >= 10
        - track:                 hay dominios, o (hotspots + clinvar), o (clinvar + local)
    """
    selected: list[str] = []

    if n_local >= 2:
        selected.append("acmg_sunburst")
        log.info(
            "[VIZ2][%s] acmg_sunburst — composicion ACMG x tipo (%d variantes)",
            gene, n_local,
        )
        selected.append("position_scatter")
        log.info(
            "[VIZ2][%s] position_scatter — posicion x severidad x tipo (3 dim.)",
            gene,
        )

    if n_local > 0 and n_clinvar >= 10:
        selected.append("clinvar_local_overlay")
        log.info(
            "[VIZ2][%s] clinvar_local_overlay — ClinVar(n=%d) + local(n=%d)",
            gene, n_clinvar, n_local,
        )

    if domains or (hotspots and n_clinvar > 0) or (n_clinvar > 0 and n_local > 0):
        selected.append("track")
        log.info(
            "[VIZ2][%s] track — mapa integrado (dominios=%d hotspots=%d clinvar=%d)",
            gene, len(domains), len(hotspots), n_clinvar,
        )

    if not selected:
        log.info("[VIZ2][%s] sin datos suficientes para visualizaciones v2", gene)
    return selected


def select_global_plots_v2(
    all_data: dict,
    subunit_map: dict,
) -> list[str]:
    """
    Decide qué visualizaciones globales v2 generar.

    Criterios:
        - acmg_comparison:    >= 2 genes con variantes locales
        - variant_spectrum:   >= 3 genes con variantes locales
        - subunit_comparison: datos ClinVar en ambas subunidades (40S y 60S)
        - data_overview:      >= 3 genes en total
    """
    selected: list[str] = []

    genes_with_local = [
        g for g, d in all_data.items()
        if d.get("local_df") is not None
        and not d["local_df"].empty
        and d["local_df"]["position"].notna().any()
    ]
    genes_with_clinvar = [g for g, d in all_data.items() if d.get("clinvar_positions")]
    subunits_with_clinvar = {subunit_map.get(g) for g in genes_with_clinvar}

    if len(genes_with_local) >= 2:
        selected.append("acmg_comparison")
        log.info(
            "[VIZ2][global] acmg_comparison — %d genes con datos locales",
            len(genes_with_local),
        )

    if len(genes_with_local) >= 3:
        selected.append("variant_spectrum")
        log.info(
            "[VIZ2][global] variant_spectrum — espectro mutacional %d genes",
            len(genes_with_local),
        )

    if len(subunits_with_clinvar) >= 2 and len(genes_with_clinvar) >= 3:
        selected.append("subunit_comparison")
        log.info(
            "[VIZ2][global] subunit_comparison — 40S vs 60S (%d genes ClinVar)",
            len(genes_with_clinvar),
        )

    if len(all_data) >= 3:
        selected.append("data_overview")
        log.info("[VIZ2][global] data_overview — bubble chart %d genes", len(all_data))

    return selected


# ---------------------------------------------------------------------------
# Puntos de entrada públicos v2
# ---------------------------------------------------------------------------

def run_gene_visualizations_v2(
    gene: str,
    gene_df: pd.DataFrame,
    clinvar_positions: list[int],
    hotspots: list[dict],
    domains: list[dict],
    protein_length: int,
    cfg: dict,
    out_dir: Path,
) -> list[Path]:
    """
    Genera las visualizaciones alternativas (v2) para un gen individual.

    Selección automática basada en datos disponibles. Guarda en out_dir
    (graphical_results2/) sin tocar graphical_results/.

    Args:
        gene:              Nombre del gen.
        gene_df:           DataFrame con variantes locales (columna 'position' incluida).
        clinvar_positions: Posiciones ClinVar patogénicas del gen.
        hotspots:          Hotspots estadísticamente significativos.
        domains:           Dominios UniProt solapantes con hotspots.
        protein_length:    Longitud de la proteína en aa.
        cfg:               Configuración completa (visualization, output).
        out_dir:           Directorio de salida (graphical_results2/).

    Returns:
        Lista de rutas de archivos generados.
    """
    vis = cfg.get("visualization", {})
    acmg_colors = vis.get("acmg_colors", {})
    output_cfg  = cfg.get("output", {})
    formats = output_cfg.get("formats", list(_EXPORT_FORMATS))
    width   = output_cfg.get("width", 1400)
    scale   = output_cfg.get("dpi_scale", 3)

    n_local = int(gene_df["position"].notna().sum()) if not gene_df.empty else 0
    plot_types = select_gene_plots_v2(
        gene, n_local, len(clinvar_positions), hotspots, domains
    )

    saved_paths: list[Path] = []
    out_dir.mkdir(parents=True, exist_ok=True)

    for ptype in plot_types:
        try:
            if ptype == "acmg_sunburst":
                df = gene_df.dropna(subset=["position"])
                fig = _acmg_sunburst_fig(df, gene, acmg_colors)
                saved_paths += _save_figure(
                    fig, out_dir / f"{gene}_acmg_sunburst", formats, 800, 700, scale
                )

            elif ptype == "position_scatter":
                fig = _position_scatter_fig(
                    gene_df, gene, protein_length, hotspots, acmg_colors
                )
                saved_paths += _save_figure(
                    fig, out_dir / f"{gene}_position_scatter", formats, width, 550, scale
                )

            elif ptype == "clinvar_local_overlay":
                fig = _clinvar_local_overlay_fig(
                    gene_df, clinvar_positions, gene,
                    protein_length, acmg_colors, hotspots,
                )
                saved_paths += _save_figure(
                    fig, out_dir / f"{gene}_clinvar_overlay", formats, width, 500, scale
                )

            elif ptype == "track":
                fig = _track_fig(
                    gene_df, clinvar_positions, hotspots,
                    domains, protein_length, gene, acmg_colors,
                )
                n_tracks = sum([
                    bool(domains), bool(clinvar_positions),
                    not gene_df.dropna(subset=["position"]).empty,
                    bool(hotspots),
                ])
                saved_paths += _save_figure(
                    fig, out_dir / f"{gene}_track",
                    formats, width, max(450, 130 * n_tracks + 120), scale,
                )

        except Exception as exc:
            log.warning("[VIZ2][%s] Error en '%s': %s", gene, ptype, exc)

    return saved_paths


def run_global_visualizations_v2(
    all_data: dict,
    protein_lengths: dict,
    subunit_map: dict,
    cfg: dict,
    out_dir: Path,
) -> list[Path]:
    """
    Genera las visualizaciones globales alternativas (v2) multi-gen.

    Args:
        all_data:        {gene: {"local_df", "clinvar_positions", "hotspots", "domains"}}
        protein_lengths: {gene: int}
        subunit_map:     {gene: "40S"|"60S"}
        cfg:             Configuración completa.
        out_dir:         Directorio de salida (graphical_results2/).

    Returns:
        Lista de rutas de archivos generados.
    """
    vis = cfg.get("visualization", {})
    acmg_colors = vis.get("acmg_colors", {})
    output_cfg  = cfg.get("output", {})
    formats = output_cfg.get("formats", list(_EXPORT_FORMATS))
    width   = output_cfg.get("width", 1400)
    scale   = output_cfg.get("dpi_scale", 3)

    plot_types = select_global_plots_v2(all_data, subunit_map)
    saved_paths: list[Path] = []
    out_dir.mkdir(parents=True, exist_ok=True)

    for ptype in plot_types:
        try:
            if ptype == "acmg_comparison":
                fig = _acmg_comparison_bar(all_data, acmg_colors)
                saved_paths += _save_figure(
                    fig, out_dir / "global_acmg_comparison", formats, width, 560, scale
                )

            elif ptype == "variant_spectrum":
                fig = _variant_spectrum_heatmap_v2(all_data)
                saved_paths += _save_figure(
                    fig, out_dir / "global_variant_spectrum", formats, 960, 580, scale
                )

            elif ptype == "subunit_comparison":
                fig = _subunit_comparison_fig(all_data, protein_lengths, subunit_map)
                saved_paths += _save_figure(
                    fig, out_dir / "global_subunit_comparison", formats, 820, 650, scale
                )

            elif ptype == "data_overview":
                fig = _data_overview_bubble(all_data, protein_lengths, subunit_map)
                saved_paths += _save_figure(
                    fig, out_dir / "global_data_overview", formats, width, 600, scale
                )

        except Exception as exc:
            log.warning("[VIZ2][global] Error en '%s': %s", ptype, exc)

    return saved_paths
