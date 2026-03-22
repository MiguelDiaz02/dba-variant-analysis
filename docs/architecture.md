# Arquitectura del sistema — DBA Variant Analysis Pipeline

## Visión general

El pipeline está organizado como un conjunto de módulos Python independientes
orquestados por un punto de entrada central (`pipeline.py`). Cada módulo tiene
una responsabilidad única y no mantiene estado global.

## Flujo de datos

```
config.yaml
    │
    ▼
pipeline.py (orquestador + CLI)
    │
    ├─── data_loader.py ──────► DataFrame limpio (variantes locales)
    │                                │
    ├─── annotator.py ───────► ClinVar positions + UniProt features
    │    ├── variant_summary.txt.gz (local)
    │    └── NCBI Entrez API (fallback)
    │                                │
    ├─── statistics.py ──────► Hotspots detectados (list[dict])
    │                                │
    ├─── visualizer.py ──────► Figuras exportadas (graphical_results/)
    │    └── Plotly (graph_background/plotly.py)
    │                                │
    └─── reporter.py ────────► Reporte Markdown (reports/)
```

## Módulos

### `data_loader.py`
- Lee el CSV de variantes con doble cabecera y encoding BOM
- Normaliza columnas usando búsqueda case-insensitive + parcial
- Exporta `extract_position()` y `classify_variant()` usados por otros módulos
- Define `_HGVS_POSITION_RE` — el único regex HGVS del sistema

### `annotator.py`
- Importa `extract_position` de `data_loader` para consistencia de regex
- Soporta dos modos de ClinVar: archivo local (preferido) o API Entrez
- Devuelve listas de posiciones — no DataFrames — para facilitar el consumo

### `statistics.py`
- No importa ningún módulo interno (sin dependencias circulares)
- `find_hotspots()` aplica corrección de Bonferroni internamente
- `merge_overlapping_hotspots()` es llamada automáticamente desde `find_hotspots()`

### `visualizer.py`
- Selección automática de plots documentada en log con justificación explícita
- Usa `plotly.graph_objects` para control granular de cada traza
- Fallback a HTML si `kaleido` no está instalado (notificado en log)
- Stems de lollipop implementados como `fig.add_shape()` (no traces) para colores independientes

### `reporter.py`
- No importa módulos internos — recibe todos los datos como argumentos
- Genera Markdown con estructura de artículo científico
- Los captions de figuras se generan automáticamente del nombre de archivo

### `pipeline.py`
- Único punto de entrada: `python pipeline.py [opciones]`
- Cada módulo puede saltarse con flags `--skip-*`
- Errores no fatales (p. ej. API down) se registran como WARNING y continúa

## Decisiones de diseño

### ¿Por qué Plotly y no matplotlib?
El directorio `graph_background/plotly.py` es el repositorio fuente de Plotly.
Usar Plotly permite:
- Interactividad (HTML) además de exportación estática (PNG/SVG/PDF)
- Mejor control de colores y símbolos por traza
- Formato publication-ready con `kaleido`

### ¿Por qué ClinVar en archivo y no solo API?
El archivo local (~414 MB comprimido) es más rápido (segundos vs. minutos),
más reproducible (versión fija) y no depende de conectividad. La API es el
fallback para cuando no se dispone del archivo.

### ¿Por qué Bonferroni y no FDR?
La corrección de Bonferroni es conservadora pero transparente e interpretable.
Para hotspot detection en proteínas cortas (100–300 aa), el número de tests
es manejable (< 300), por lo que el costo en potencia estadística es razonable.
FDR (Benjamini-Hochberg) puede considerarse en análisis con muchos genes simultáneos.

### ¿Por qué el regex HGVS está en `data_loader.py`?
Es la única fuente de verdad para extraer posiciones proteicas de strings HGVS.
`annotator.py` lo importa directamente para garantizar consistencia. El notebook
original tenía tres variantes distintas del mismo regex, causando resultados
inconsistentes.

## Limitaciones conocidas

1. **Proteínas no en config.yaml**: cualquier gen no listado en el config es ignorado.
2. **ClinVar columnas**: si NCBI cambia el esquema del TSV, `filter_clinvar_pathogenic`
   y `get_pathogenic_positions_from_clinvar` necesitarán actualización.
3. **kaleido**: sin él, la exportación PNG/SVG no funciona. El pipeline hace fallback
   a HTML, que no es adecuado para publicación impresa.
4. **Ventana fija**: el tamaño de ventana del hotspot es el mismo para todas las
   proteínas. Proteínas muy cortas (< 50 aa) o muy largas podrían requerir ajuste.

## Extensión del pipeline

### Añadir un nuevo gen
1. Añadir entrada en `config.yaml → genes` con name, subunit, protein_length, uniprot_id
2. No se requiere ningún cambio en el código

### Añadir una nueva fuente de variantes
1. Crear función en `annotator.py` que devuelva `list[int]`
2. Añadir opción en `pipeline.py` (flag CLI + llamada condicional)

### Añadir un nuevo tipo de gráfico
1. Implementar función `_nuevo_fig(...)` en `visualizer.py`
2. Añadir criterio de selección en `select_gene_plots()` o `select_global_plots()`
3. Añadir `elif ptype == "nuevo"` en `run_gene_visualizations()` o `run_global_visualizations()`
