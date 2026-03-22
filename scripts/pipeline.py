"""
pipeline.py
===========
Punto de entrada del sistema de análisis de variantes DBA.

Orquesta la ejecución secuencial de los módulos:
    1. data_loader  — carga y limpieza del CSV local
    2. annotator    — integración con ClinVar y UniProt
    3. statistics   — detección de hotspots
    4. visualizer   — generación de gráficos con selección automática
    5. reporter     — reporte científico en Markdown

Uso:
    python pipeline.py
    python pipeline.py --config /ruta/config.yaml
    python pipeline.py --genes RPS19 RPL5
    python pipeline.py --skip-hotspots --skip-report
    python pipeline.py --skip-clinvar-file  # usa API Entrez en lugar del archivo local
    python pipeline.py --verbose

Salidas:
    graphical_results/   — figuras PNG/SVG (o HTML si falta kaleido)
    reports/             — reporte científico en Markdown
    pipeline.log         — log completo de la ejecución
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import yaml

# Importaciones de módulos del pipeline (en el mismo directorio)
sys.path.insert(0, str(Path(__file__).parent))
import annotator
import data_loader
import reporter
import statistics as stats_mod
import visualizer


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def _setup_logging(verbose: bool, log_file: Path) -> None:
    """Configura el sistema de logging con salida a consola y a archivo."""
    level = logging.DEBUG if verbose else logging.INFO
    fmt = "%(asctime)s [%(levelname)-8s] %(name)s — %(message)s"
    datefmt = "%H:%M:%S"

    # StreamHandler con UTF-8 explícito para evitar UnicodeEncodeError en Windows (CP1252)
    stream_handler = logging.StreamHandler(sys.stdout)
    if hasattr(stream_handler.stream, "reconfigure"):
        stream_handler.stream.reconfigure(encoding="utf-8")  # Python 3.7+

    handlers: list[logging.Handler] = [
        stream_handler,
        logging.FileHandler(log_file, encoding="utf-8"),
    ]
    for h in handlers:
        h.setFormatter(logging.Formatter(fmt, datefmt))
        h.setLevel(level)

    root = logging.getLogger()
    root.setLevel(level)
    for h in root.handlers[:]:
        root.removeHandler(h)
    for h in handlers:
        root.addHandler(h)


log = logging.getLogger("pipeline")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pipeline.py",
        description="Pipeline de análisis de variantes DBA — Diamond-Blackfan Anemia",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--config", type=Path,
        default=Path(__file__).parent / "config.yaml",
        help="Ruta al archivo de configuración YAML.",
    )
    p.add_argument(
        "--genes", nargs="+", metavar="GEN",
        help="Genes a analizar. Si se omite, se usan todos los de config.yaml.",
    )
    p.add_argument(
        "--skip-clinvar-file", action="store_true",
        help="Usar la API de NCBI en lugar del archivo ClinVar local.",
    )
    p.add_argument(
        "--skip-hotspots", action="store_true",
        help="Omitir el análisis estadístico de hotspots.",
    )
    p.add_argument(
        "--skip-viz", action="store_true",
        help="Omitir la generación de visualizaciones.",
    )
    p.add_argument(
        "--skip-report", action="store_true",
        help="Omitir la generación del reporte científico.",
    )
    p.add_argument(
        "--verbose", "-v", action="store_true",
        help="Activar logging en modo DEBUG.",
    )
    return p


# ---------------------------------------------------------------------------
# Carga de configuración
# ---------------------------------------------------------------------------

def _load_config(config_path: Path) -> dict:
    if not config_path.exists():
        raise FileNotFoundError(
            f"Archivo de configuración no encontrado: {config_path}\n"
            "Asegúrate de que config.yaml esté en el mismo directorio que pipeline.py."
        )
    with config_path.open(encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    return cfg


def _resolve_paths(cfg: dict, base: Path) -> dict[str, Path]:
    """Resuelve todas las rutas del config relativas al directorio del config."""
    paths = cfg.get("paths", {})
    return {
        "csv":       (base / paths["variants_csv"]).resolve(),
        "clinvar":   (base / paths["clinvar_gz"]).resolve(),
        "out_viz":   (base / paths.get("graphical_results",  "../graphical_results")).resolve(),
        "out_viz2":  (base / paths.get("graphical_results2", "../graphical_results2")).resolve(),
        "out_rep":   (base / paths.get("reports",   "../reports")).resolve(),
        "log_file":  (base / paths.get("log_file",  "../pipeline.log")).resolve(),
    }


def _build_gene_dicts(cfg: dict) -> tuple[list[str], dict, dict, dict]:
    """Construye listas y diccionarios de genes a partir del bloque 'genes' del config."""
    genes_list, prot_len, subunit_map, uniprot_ids = [], {}, {}, {}
    for entry in cfg.get("genes", []):
        name = entry["name"]
        genes_list.append(name)
        prot_len[name]    = entry["protein_length"]
        subunit_map[name] = entry["subunit"]
        uniprot_ids[name] = entry.get("uniprot_id", "")
    return genes_list, prot_len, subunit_map, uniprot_ids


# ---------------------------------------------------------------------------
# Pipeline principal
# ---------------------------------------------------------------------------

def run_pipeline(args: argparse.Namespace) -> int:
    """
    Ejecuta el pipeline completo en el orden establecido.

    Returns:
        0 si completó sin errores críticos, 1 si hubo al menos un error.
    """
    # ------------------------------------------------------------------
    # Configuración
    # ------------------------------------------------------------------
    try:
        cfg = _load_config(args.config)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    base = args.config.parent
    resolved = _resolve_paths(cfg, base)
    _setup_logging(args.verbose, resolved["log_file"])

    log.info("=" * 60)
    log.info("DBA Variant Analysis Pipeline")
    log.info("Config: %s", args.config)
    log.info("=" * 60)

    all_genes, prot_len, subunit_map, uniprot_ids = _build_gene_dicts(cfg)
    genes = args.genes if args.genes else all_genes

    # Validar genes solicitados
    unknown = set(genes) - set(all_genes)
    if unknown:
        log.error("Gen(es) no encontrado(s) en config.yaml: %s", sorted(unknown))
        log.error("Genes disponibles: %s", all_genes)
        return 1

    vis_cfg = cfg.get("visualization", {})
    valid_acmg = list(vis_cfg.get("acmg_colors", {}).keys())

    had_error = False

    # ------------------------------------------------------------------
    # MÓDULO 1: Carga y limpieza de datos locales
    # ------------------------------------------------------------------
    log.info("─" * 40)
    log.info("MÓDULO 1/5: Carga de variantes locales")
    try:
        raw_df = data_loader.load_variants(resolved["csv"])
        full_df = data_loader.prepare_dataframe(raw_df, valid_acmg)
    except (FileNotFoundError, ValueError) as exc:
        log.error("Error en carga de datos: %s", exc)
        return 1

    # ------------------------------------------------------------------
    # MÓDULO 2: ClinVar — carga del archivo (si no se usa API)
    # ------------------------------------------------------------------
    log.info("─" * 40)
    log.info("MÓDULO 2/5: Integración con ClinVar")
    clinvar_path_df = None

    if not args.skip_hotspots and not args.skip_clinvar_file:
        try:
            clinvar_df = annotator.load_clinvar_file(resolved["clinvar"])
            clinvar_path_df = annotator.filter_clinvar_pathogenic(clinvar_df)
        except FileNotFoundError as exc:
            log.warning("%s", exc)
            log.warning("Se intentará usar la API de NCBI Entrez como alternativa.")
        except Exception as exc:
            log.error("Error inesperado al cargar ClinVar: %s", exc)
            had_error = True

    # ------------------------------------------------------------------
    # MÓDULOS 3 (stats) y 4 (viz) por gen
    # ------------------------------------------------------------------
    log.info("─" * 40)
    log.info("MÓDULOS 3-4/5: Análisis estadístico y visualización por gen")

    all_results: dict[str, dict] = {}
    all_figure_paths: list[Path] = []

    for gene in genes:
        log.info("┌── Procesando: %s", gene)
        protein_length = prot_len[gene]
        gene_df = full_df[full_df["Gene"] == gene].copy()
        n_local = int(gene_df["position"].notna().sum())
        log.info("│   Variantes locales con posición: %d", n_local)

        # 3a. Posiciones ClinVar
        clinvar_positions: list[int] = []
        if not args.skip_hotspots:
            try:
                clinvar_positions = annotator.get_clinvar_positions(
                    gene=gene,
                    protein_length=protein_length,
                    cfg=cfg,
                    clinvar_path_df=clinvar_path_df,
                )
            except Exception as exc:
                log.warning("│   ClinVar falló para %s: %s — se omite", gene, exc)
                had_error = True

        # 3b. Hotspots
        hotspots: list[dict] = []
        if not args.skip_hotspots and len(clinvar_positions) >= cfg.get("analysis", {}).get("min_variants_for_hotspot", 20):
            try:
                analysis_cfg = cfg.get("analysis", {})
                hotspots = stats_mod.find_hotspots(
                    positions=clinvar_positions,
                    protein_length=protein_length,
                    window=analysis_cfg.get("hotspot_window", 10),
                    alpha=analysis_cfg.get("alpha", 0.05),
                )
                summary = stats_mod.summarize_hotspots(gene, hotspots, protein_length)
                if hotspots:
                    log.info("│   Hotspots: %s", summary["hotspot_regions"])
                else:
                    log.info("│   Sin hotspots significativos")
            except Exception as exc:
                log.warning("│   Error en detección de hotspots para %s: %s", gene, exc)
                had_error = True
        elif not args.skip_hotspots:
            log.info("│   Hotspots omitidos (%d posiciones ClinVar < mínimo)", len(clinvar_positions))

        # 3c. Anotación UniProt
        domains: list[dict] = []
        if hotspots and uniprot_ids.get(gene):
            try:
                features = annotator.get_uniprot_features(uniprot_ids[gene])
                domains = annotator.find_overlapping_domains(hotspots, features)
                if domains:
                    log.info("│   Dominios solapantes: %d", len(domains))
            except Exception as exc:
                log.warning("│   UniProt falló para %s: %s", gene, exc)
                had_error = True

        # 4. Visualizaciones por gen (lollipop + KDE/histograma — graphical_results/)
        gene_figs: list[Path] = []
        if not args.skip_viz:
            try:
                resolved["out_viz"].mkdir(parents=True, exist_ok=True)
                gene_figs = visualizer.run_gene_visualizations(
                    gene=gene,
                    gene_df=gene_df,
                    clinvar_positions=clinvar_positions,
                    hotspots=hotspots,
                    protein_length=protein_length,
                    cfg=cfg,
                    out_dir=resolved["out_viz"],
                )
                all_figure_paths.extend(gene_figs)
            except Exception as exc:
                log.warning("│   Error en visualizacion para %s: %s", gene, exc)
                had_error = True

        # 4b. Visualizaciones v2 (sunburst, scatter, overlay, track — graphical_results2/)
        if not args.skip_viz:
            try:
                gene_figs_v2 = visualizer.run_gene_visualizations_v2(
                    gene=gene,
                    gene_df=gene_df,
                    clinvar_positions=clinvar_positions,
                    hotspots=hotspots,
                    domains=domains,
                    protein_length=protein_length,
                    cfg=cfg,
                    out_dir=resolved["out_viz2"],
                )
                all_figure_paths.extend(gene_figs_v2)
            except Exception as exc:
                log.warning("│   Error en visualizaciones v2 para %s: %s", gene, exc)
                had_error = True

        all_results[gene] = {
            "n_local":           n_local,
            "local_df":          gene_df,
            "clinvar_positions": clinvar_positions,
            "hotspots":          hotspots,
            "domains":           domains,
            "figure_paths":      gene_figs,
        }
        log.info("└── %s completado", gene)

    # ------------------------------------------------------------------
    # Visualizaciones globales v1 (multi-overview, KDE multi, heatmap, violin)
    # ------------------------------------------------------------------
    if not args.skip_viz and len(genes) > 1:
        log.info("─" * 40)
        log.info("Generando visualizaciones globales (v1)...")
        try:
            global_figs = visualizer.run_global_visualizations(
                all_data=all_results,
                protein_lengths=prot_len,
                subunit_map=subunit_map,
                cfg=cfg,
                out_dir=resolved["out_viz"],
            )
            all_figure_paths.extend(global_figs)
        except Exception as exc:
            log.warning("Error en visualizaciones globales v1: %s", exc)
            had_error = True

    # ------------------------------------------------------------------
    # Visualizaciones globales v2 (ACMG comparison, espectro, 40S/60S, bubble)
    # ------------------------------------------------------------------
    if not args.skip_viz and len(genes) > 1:
        log.info("Generando visualizaciones globales (v2)...")
        try:
            global_figs_v2 = visualizer.run_global_visualizations_v2(
                all_data=all_results,
                protein_lengths=prot_len,
                subunit_map=subunit_map,
                cfg=cfg,
                out_dir=resolved["out_viz2"],
            )
            all_figure_paths.extend(global_figs_v2)
        except Exception as exc:
            log.warning("Error en visualizaciones globales v2: %s", exc)
            had_error = True

    # ------------------------------------------------------------------
    # MÓDULO 5: Reporte científico
    # ------------------------------------------------------------------
    if not args.skip_report:
        log.info("─" * 40)
        log.info("MÓDULO 5/5: Generando reporte científico…")
        try:
            report_path = resolved["out_rep"] / "dba_variant_analysis_report.md"
            reporter.generate_report(
                results=all_results,
                cfg=cfg,
                output_path=report_path,
                figure_paths=all_figure_paths,
            )
        except Exception as exc:
            log.error("Error al generar el reporte: %s", exc)
            had_error = True

    # ------------------------------------------------------------------
    # Resumen final
    # ------------------------------------------------------------------
    log.info("=" * 60)
    log.info("RESUMEN FINAL")
    log.info("Genes procesados: %s", genes)
    log.info(
        "Figuras generadas: %d (v1: %s | v2: %s)",
        len(all_figure_paths), resolved["out_viz"], resolved["out_viz2"],
    )
    if not args.skip_report:
        log.info("Reporte: %s", resolved["out_rep"] / "dba_variant_analysis_report.md")
    log.info("Log completo: %s", resolved["log_file"])

    hotspot_summary = {
        g: [(h["start"], h["end"]) for h in r["hotspots"]]
        for g, r in all_results.items() if r["hotspots"]
    }
    if hotspot_summary:
        log.info("Hotspots detectados:")
        for gene, regions in hotspot_summary.items():
            log.info("  %-10s → %s", gene, regions)
    else:
        log.info("Sin hotspots estadísticamente significativos en los genes analizados.")

    log.info("Estado: %s", "COMPLETADO CON ADVERTENCIAS" if had_error else "COMPLETADO SIN ERRORES")
    log.info("=" * 60)

    return 1 if had_error else 0


# ---------------------------------------------------------------------------
# Entrypoint
# ---------------------------------------------------------------------------

def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()
    sys.exit(run_pipeline(args))


if __name__ == "__main__":
    main()
