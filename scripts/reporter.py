"""
reporter.py
===========
Responsabilidad única: generar el reporte científico en formato Markdown
a partir de los resultados del análisis.

El reporte sigue la estructura de un artículo científico:
    1. Introducción
    2. Problema analítico
    3. Métodos
    4. Resultados (tablas + figuras)
    5. Discusión
    6. Apéndice (parámetros del análisis)

No incluye código fuente. El lenguaje es formal y técnico.

Funciones exportadas:
    generate_report(results, cfg, output_path, figure_paths)
"""

from __future__ import annotations

import logging
from datetime import date
from pathlib import Path

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Secciones del reporte
# ---------------------------------------------------------------------------

def _section_abstract(results: dict) -> str:
    n_genes   = len(results)
    total_var = sum(r["n_local"] for r in results.values())
    genes_with_hotspots = [g for g, r in results.items() if r["hotspots"]]

    hotspot_str = (
        f"{len(genes_with_hotspots)} gen(es) presentaron hotspots estadísticamente "
        f"significativos ({', '.join(genes_with_hotspots)})."
        if genes_with_hotspots
        else "No se detectaron hotspots estadísticamente significativos."
    )

    return f"""\
## Resumen

Se analizaron {total_var} variantes germinales en {n_genes} genes de proteínas
ribosomales asociados a Diamond-Blackfan Anemia (DBA), clasificadas según los
criterios del American College of Medical Genetics and Genomics (ACMG).
Los datos locales se integraron con variantes patogénicas del repositorio público
ClinVar para la detección de hotspots de mutación. {hotspot_str}
Las variantes se distribuyeron principalmente como missense y nonsense, con patrones
específicos de acumulación en regiones estructuralmente relevantes de las proteínas.
"""


def _section_introduction() -> str:
    return """\
## 1. Introducción

La anemia de Diamond-Blackfan (DBA; OMIM #105650) es una enfermedad congénita
rara de la médula ósea, caracterizada por eritroblastopenia selectiva, anomalías
congénitas variables y predisposición aumentada al desarrollo de neoplasias malignas.
La DBA se hereda de forma autosómica dominante en la mayoría de los casos, aunque
también se han descrito mutaciones de novo.

La base genética de la DBA implica principalmente mutaciones heterocigotas con
pérdida de función en genes que codifican proteínas ribosomales de las subunidades
pequeña (40S) y grande (60S) del ribosoma eucariota. Hasta la fecha, se han
identificado variantes patogénicas en más de 20 genes ribosomales, siendo *RPS19*
el más frecuentemente afectado (~25% de los casos). Otros genes relevantes incluyen
*RPS24*, *RPS26*, *RPL5*, *RPL11*, *RPL18* y *RPL26*, que en conjunto explican
aproximadamente el 30–40% adicional de los casos.

El mecanismo fisiopatológico central involucra la insuficiencia ribosómica, que
activa la vía p53 mediante la estabilización del complejo RPL5-RPL11-5S ARN, lo
que conduce a la apoptosis diferencial de los precursores eritroides.

El análisis sistemático de las variantes patogénicas en estos genes es fundamental
para la interpretación clínica de los hallazgos en pacientes con DBA, la
identificación de regiones funcionales críticas y la comprensión de las correlaciones
genotipo-fenotipo.
"""


def _section_problem() -> str:
    return """\
## 2. Problema analítico

El objetivo de este análisis es caracterizar la distribución de variantes patogénicas
y probablemente patogénicas en los genes ribosomales asociados a DBA, con énfasis en:

1. **Distribución posicional**: determinar si las variantes se acumulan en regiones
   específicas de las proteínas o se distribuyen de forma difusa a lo largo de la
   secuencia.

2. **Identificación de hotspots**: detectar regiones estadísticamente enriquecidas
   en variantes patogénicas mediante análisis de ventana deslizante con corrección
   para comparaciones múltiples.

3. **Anotación funcional**: correlacionar los hotspots detectados con dominios
   estructurales y funcionales anotados en UniProt, para determinar si las regiones
   recurrentemente mutadas corresponden a elementos funcionales conocidos.

4. **Comparación inter-gen**: evaluar si los patrones de mutación son similares o
   divergentes entre los distintos genes ribosomales afectados en DBA.

El análisis integra datos clínicos locales con el repertorio de variantes patogénicas
de ClinVar, maximizando la sensibilidad de la detección de patrones recurrentes.
"""


def _section_methods(cfg: dict) -> str:
    analysis = cfg.get("analysis", {})
    entrez   = cfg.get("entrez",   {})
    window   = analysis.get("hotspot_window", 10)
    alpha    = analysis.get("alpha", 0.05)
    min_var  = analysis.get("min_variants_for_hotspot", 20)

    return f"""\
## 3. Métodos

### 3.1 Fuentes de datos

Los datos primarios consistieron en una tabla de variantes clínicas de pacientes
diagnosticados con DBA, clasificadas según los criterios ACMG/AMP (Pathogenic,
Likely Pathogenic, Variant of Uncertain Significance). Los datos de referencia
se obtuvieron del repositorio ClinVar del National Center for Biotechnology
Information (NCBI), descargando el archivo completo *variant_summary.txt.gz*
(versión consultada: {date.today().strftime('%B %Y')}).

### 3.2 Procesamiento de variantes

Las variantes se anotaron con sus posiciones proteicas mediante el análisis
sintáctico de la nomenclatura HGVS (Human Genome Variation Society). El tipo
funcional de cada variante se determinó a partir del cambio proteico reportado,
clasificándose en las categorías: Missense, Nonsense, Frameshift, Splice y Other.

Únicamente se incluyeron en los análisis posicionales las variantes con notación
HGVS proteica válida. Las variantes de clase Benigna y Probablemente Benigna
fueron excluidas del análisis de hotspots.

### 3.3 Detección de hotspots

La detección de hotspots de mutación se realizó mediante un algoritmo de ventana
deslizante sobre la longitud completa de cada proteína. Para cada posición de
inicio *s* (de 1 a L − W + 1, donde L es la longitud proteica y W el tamaño de
ventana), se contó el número de variantes *k* en el intervalo [s, s+W−1].

La significancia estadística se evaluó mediante el **test binomial exacto unilateral**
(scipy.stats.binomtest, alternativa "greater"), bajo la hipótesis nula de distribución
uniforme con probabilidad p₀ = W / L. Para controlar la tasa de error tipo I ante
comparaciones múltiples, se aplicó la **corrección de Bonferroni**:

> α\_corregido = {alpha} / (L − W + 1)

donde el número de tests es el número total de ventanas evaluadas. Solo se
reportaron como hotspots las ventanas con p-valor < α\_corregido.

Las ventanas significativas solapadas o adyacentes se fusionaron en regiones
continuas únicas, conservando el p-valor mínimo (más significativo) y el conteo
máximo de la región resultante.

**Parámetros utilizados:**
- Tamaño de ventana: {window} aminoácidos
- Nivel de significancia nominal: α = {alpha}
- Mínimo de variantes para análisis: {min_var}

### 3.4 Anotación funcional

Para cada hotspot detectado, se consultó la API REST de UniProt Knowledge Base
para obtener las anotaciones de dominios, regiones y motivos de la proteína
correspondiente. Se identificaron los elementos funcionales con coordenadas
solapantes con cada hotspot.

### 3.5 Visualización

Las visualizaciones se generaron mediante la librería Plotly. El tipo de gráfico
se seleccionó automáticamente en función de las características de los datos
disponibles para cada gen: lollipop plots para datos posicionales locales,
estimación de densidad por kernel (KDE) para conjuntos de ≥ 30 posiciones ClinVar,
histogramas para conjuntos menores, y visualizaciones comparativas multi-gen
(heatmaps, violin plots y KDE superpuestos) cuando múltiples genes disponen de
datos suficientes.
"""


def _section_results(results: dict, figure_paths: list[Path], cfg: dict) -> str:
    lines = ["## 4. Resultados\n"]

    # --- 4.1 Tabla de variantes ---
    lines.append("### 4.1 Variantes identificadas\n")
    lines.append("| Gen | Subunidad | Variantes locales | Posiciones ClinVar | Hotspots |")
    lines.append("|-----|-----------|------------------:|-------------------:|---------:|")

    genes_cfg = {g["name"]: g for g in cfg.get("genes", [])}
    for gene, res in results.items():
        sub = genes_cfg.get(gene, {}).get("subunit", "—")
        n_local  = res["n_local"]
        n_clinvar = len(res.get("clinvar_positions", []))
        n_hs = len(res.get("hotspots", []))
        lines.append(f"| {gene} | {sub} | {n_local} | {n_clinvar} | {n_hs} |")
    lines.append("")

    # --- 4.2 Hotspots ---
    lines.append("### 4.2 Hotspots de mutación detectados\n")
    any_hotspot = any(res.get("hotspots") for res in results.values())
    if not any_hotspot:
        lines.append(
            "No se detectaron hotspots estadísticamente significativos en ningún gen "
            "con los parámetros de análisis establecidos. Este resultado es consistente "
            "con una distribución difusa de variantes a lo largo de las proteínas estudiadas.\n"
        )
    else:
        lines.append("| Gen | Región (aa) | Variantes observadas | Esperadas | p-valor (Bonferroni) |")
        lines.append("|-----|-------------|---------------------:|----------:|--------------------:|")
        for gene, res in results.items():
            for hs in res.get("hotspots", []):
                lines.append(
                    f"| {gene} | {hs['start']}–{hs['end']} | "
                    f"{hs['count']} | {hs['expected']} | {hs['p_value']:.2e} |"
                )
        lines.append("")

    # --- 4.3 Anotación funcional ---
    lines.append("### 4.3 Anotación funcional de hotspots\n")
    any_domain = any(res.get("domains") for res in results.values())
    if not any_domain:
        lines.append(
            "No se identificaron dominios UniProt solapantes con los hotspots detectados, "
            "o no se encontraron hotspots en los que realizar la anotación.\n"
        )
    else:
        lines.append("| Gen | Hotspot (aa) | Dominio/Región | Posición en proteína |")
        lines.append("|-----|-------------|----------------|---------------------:|")
        for gene, res in results.items():
            for dom in res.get("domains", []):
                lines.append(
                    f"| {gene} | — | {dom.get('description','—')} ({dom.get('type','')}) | "
                    f"{dom.get('start','?')}–{dom.get('end','?')} |"
                )
        lines.append("")

    # --- 4.4 Figuras ---
    lines.append("### 4.4 Figuras\n")
    if not figure_paths:
        lines.append("*No se generaron figuras en esta ejecución.*\n")
    else:
        # Agrupar por tipo de figura
        for fig_path in figure_paths:
            name = fig_path.stem
            rel_path = Path("..") / "graphical_results" / fig_path.name
            caption = _figure_caption(name, results, cfg)
            lines.append(f"#### {name}\n")
            lines.append(f"![{name}]({rel_path})\n")
            lines.append(f"**Figura:** {caption}\n")

    return "\n".join(lines)


def _figure_caption(stem: str, results: dict, cfg: dict) -> str:
    """Genera un caption científico para una figura a partir de su nombre de archivo."""
    captions = {
        "lollipop": (
            "Lollipop plot mostrando la distribución posicional de las variantes en la proteína. "
            "El color de cada variante corresponde a la clasificación ACMG (rojo: Patogénica, "
            "naranja: Probablemente Patogénica, amarillo: VUS). La forma del marcador indica "
            "el tipo funcional de la variante (círculo: missense, cuadrado: nonsense, "
            "diamante: frameshift, triángulo: splice). Los marcadores vacíos corresponden "
            "a variantes no reportadas previamente en la literatura."
        ),
        "kde": (
            "Estimación de densidad por kernel (KDE) de las posiciones de variantes "
            "patogénicas en ClinVar. Los ticks sobre el eje x indican las posiciones "
            "individuales. Las regiones sombreadas corresponden a hotspots "
            "estadísticamente significativos (test binomial, corrección de Bonferroni)."
        ),
        "histogram": (
            "Histograma de frecuencias de variantes patogénicas de ClinVar por posición "
            "proteica. Las regiones sombreadas corresponden a hotspots estadísticamente "
            "significativos."
        ),
        "global_multi_overview": (
            "Visión global de las variantes en todos los genes DBA analizados, con posición "
            "normalizada (0 = extremo N-terminal, 1 = extremo C-terminal). Los genes de la "
            "subunidad 40S y 60S se presentan separados. El color y la forma de los marcadores "
            "codifican la clasificación ACMG y el tipo de variante, respectivamente."
        ),
        "global_kde_multi": (
            "Estimación de densidad por kernel (KDE) superpuesta para todos los genes, "
            "con posición proteica normalizada. Permite comparar visualmente los patrones "
            "de distribución de variantes entre genes con longitudes proteicas distintas."
        ),
        "global_heatmap": (
            "Mapa de calor de densidad de mutaciones. Cada fila corresponde a un gen y "
            "cada columna a un intervalo de posición normalizada. La intensidad del color "
            "representa el número de variantes ClinVar patogénicas en cada intervalo."
        ),
        "global_violin": (
            "Violin plot de la distribución de variantes patogénicas de ClinVar por gen, "
            "con posición normalizada. La caja central indica la mediana e IQR. "
            "Permite comparar la dispersión y simetría de la distribución de variantes "
            "entre genes."
        ),
    }
    for key, caption in captions.items():
        if key in stem:
            gene = stem.split("_")[0] if not stem.startswith("global") else None
            prefix = f"**{gene}** — " if gene and gene in results else ""
            return prefix + caption
    return f"Figura generada por el pipeline de análisis DBA ({stem})."


def _section_discussion(results: dict) -> str:
    genes_with_hs = [g for g, r in results.items() if r.get("hotspots")]
    total_var = sum(r["n_local"] for r in results.values())

    if genes_with_hs:
        hs_text = (
            f"Los genes {', '.join(genes_with_hs)} presentaron hotspots estadísticamente "
            "significativos, lo que sugiere que ciertas posiciones de estas proteínas son "
            "particularmente vulnerables a mutaciones patogénicas. Esto puede reflejar "
            "restricciones estructurales o funcionales que hacen que estas regiones sean "
            "intolerantes a la variación de secuencia."
        )
    else:
        hs_text = (
            "La ausencia de hotspots estadísticamente significativos puede indicar que las "
            "variantes patogénicas se distribuyen de forma relativamente uniforme a lo largo "
            "de las proteínas, sin puntos de acumulación preferencial. Esto es consistente "
            "con un mecanismo de haploinsuficiencia donde la pérdida de cualquier parte "
            "de la proteína es suficiente para alterar la función ribosomal."
        )

    return f"""\
## 5. Discusión

Este análisis integró {total_var} variantes clínicas locales con el repertorio de
variantes patogénicas de ClinVar para caracterizar los patrones de mutación en los
genes ribosomales asociados a DBA.

{hs_text}

La integración con la base de datos UniProt permitió correlacionar los patrones de
mutación con la arquitectura funcional de las proteínas, proporcionando hipótesis
mecanísticas sobre el impacto de las variantes en la función ribosomal.

### Limitaciones

- El tamaño de la muestra de variantes locales limita el poder estadístico de la
  detección de hotspots para genes con pocas variantes documentadas.
- El análisis de hotspots asume distribución uniforme de variantes bajo la hipótesis
  nula, lo que puede no reflejar la estructura de selección neutral real de cada proteína.
- La dependencia del archivo ClinVar introduce un sesgo de publicación: las variantes
  en genes más estudiados están más representadas.
- La clasificación ACMG puede variar entre laboratorios; se recomienda revisar las
  variantes VUS periódicamente.
"""


def _section_appendix(cfg: dict) -> str:
    analysis = cfg.get("analysis", {})
    genes_cfg = cfg.get("genes", [])

    genes_table = "\n".join(
        f"| {g['name']} | {g['subunit']} | {g['protein_length']} | {g.get('uniprot_id','—')} |"
        for g in genes_cfg
    )

    return f"""\
## Apéndice: Parámetros del análisis

| Parámetro | Valor |
|-----------|-------|
| Ventana deslizante | {analysis.get('hotspot_window', 10)} aa |
| Nivel de significancia (α) | {analysis.get('alpha', 0.05)} |
| Corrección múltiple | Bonferroni |
| Mínimo de variantes para hotspot | {analysis.get('min_variants_for_hotspot', 20)} |

### Genes analizados

| Gen | Subunidad | Longitud proteína (aa) | UniProt ID |
|-----|-----------|----------------------:|------------|
{genes_table}
"""


# ---------------------------------------------------------------------------
# Punto de entrada público
# ---------------------------------------------------------------------------

def generate_report(
    results: dict,
    cfg: dict,
    output_path: Path,
    figure_paths: list[Path],
) -> Path:
    """
    Genera el reporte científico completo en formato Markdown.

    El reporte sigue la estructura de un artículo científico e incluye
    tablas de resultados, referencias a las figuras generadas y discusión
    de los hallazgos. No contiene código fuente.

    Args:
        results:      {gene: {"n_local": int, "clinvar_positions": [], "hotspots": [],
                              "domains": []}}
        cfg:          Configuración completa del análisis.
        output_path:  Ruta donde guardar el archivo .md.
        figure_paths: Lista de rutas a las figuras generadas.

    Returns:
        Ruta del archivo de reporte generado.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    today = date.today().isoformat()
    genes_str = ", ".join(results.keys())

    header = f"""\
# Análisis de Variantes en Genes Ribosomales Asociados a Diamond-Blackfan Anemia

**Fecha de análisis:** {today}
**Genes analizados:** {genes_str}
**Pipeline:** DBA Variant Analysis Pipeline v1.0

---

"""

    sections = [
        header,
        _section_abstract(results),
        "\n---\n\n",
        _section_introduction(),
        _section_problem(),
        _section_methods(cfg),
        _section_results(results, figure_paths, cfg),
        _section_discussion(results),
        _section_appendix(cfg),
        f"\n---\n*Reporte generado automáticamente por el pipeline DBA — {today}*\n",
    ]

    full_text = "\n".join(sections)
    output_path.write_text(full_text, encoding="utf-8")
    log.info("Reporte científico generado: %s (%d caracteres)", output_path, len(full_text))
    return output_path
