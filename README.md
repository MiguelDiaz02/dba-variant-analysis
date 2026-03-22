# DBA Variant Analysis Pipeline

Análisis sistemático de variantes patogénicas en genes de proteínas ribosomales
asociados a **Diamond-Blackfan Anemia (DBA)**.

## Descripción

Este pipeline integra datos clínicos locales con el repositorio público ClinVar
para caracterizar la distribución posicional de variantes patogénicas en 7 genes
ribosomales (RPS19, RPS24, RPS26, RPL5, RPL11, RPL18, RPL26), detectar hotspots
de mutación estadísticamente significativos y generar un reporte científico listo
para publicación.

## Contexto científico

La DBA es una anemia congénita rara (~7 casos por millón de nacidos vivos) causada
principalmente por mutaciones heterocigotas con pérdida de función en genes de
proteínas ribosomales. El análisis de la distribución posicional de las variantes
permite identificar regiones funcionales críticas y establecer correlaciones
genotipo-fenotipo relevantes para el diagnóstico y consejo genético.

## Requisitos

- Python 3.9+
- Archivo ClinVar: [`variant_summary.txt.gz`](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/) (opcional, ~414 MB comprimido)

## Instalación

```bash
# Clona el repositorio
git clone <url-del-repositorio>
cd paula_datos

# Instala las dependencias
pip install -r requirements.txt

# (Opcional) Para exportar PNG/SVG en alta resolución
pip install kaleido
```

## Uso

Todos los comandos se ejecutan desde la carpeta `scripts/`:

```bash
cd scripts/

# Pipeline completo (requiere variant_summary.txt.gz en additional_files/)
python pipeline.py

# Solo lollipop plots, sin análisis estadístico ni reporte
python pipeline.py --skip-hotspots --skip-report

# Analizar genes específicos
python pipeline.py --genes RPS19 RPL5 RPL11

# Usar la API de NCBI en lugar del archivo ClinVar local
python pipeline.py --skip-clinvar-file

# Modo verbose (log detallado)
python pipeline.py --verbose

# Config personalizado
python pipeline.py --config /ruta/mi_config.yaml
```

### Opciones disponibles

| Opción | Descripción |
|--------|-------------|
| `--config PATH` | Archivo YAML de configuración (por defecto: `config.yaml`) |
| `--genes GEN [GEN ...]` | Genes a analizar (por defecto: todos en config) |
| `--skip-clinvar-file` | Usar API Entrez en lugar del archivo local |
| `--skip-hotspots` | Omitir análisis estadístico de hotspots |
| `--skip-viz` | Omitir generación de visualizaciones |
| `--skip-report` | Omitir generación del reporte científico |
| `--verbose`, `-v` | Logging en modo DEBUG |

## Pipeline

```
CSV local → limpieza → ClinVar → hotspot detection → visualización → reporte
```

1. **Carga de datos** (`data_loader.py`): limpia el CSV de doble cabecera, normaliza columnas, parsea notación HGVS
2. **Integración ClinVar** (`annotator.py`): extrae posiciones patogénicas del archivo local o vía API NCBI
3. **Estadística** (`statistics.py`): ventana deslizante + test binomial + corrección Bonferroni
4. **Visualización** (`visualizer.py`): selección automática del tipo de gráfico según los datos disponibles
5. **Reporte** (`reporter.py`): documento Markdown con estructura de artículo científico

## Estructura del repositorio

```
paula_datos/
├── scripts/
│   ├── pipeline.py       ← Punto de entrada principal
│   ├── data_loader.py    ← Carga y limpieza de datos
│   ├── annotator.py      ← ClinVar + UniProt
│   ├── statistics.py     ← Detección de hotspots
│   ├── visualizer.py     ← Visualizaciones (Plotly)
│   ├── reporter.py       ← Reporte científico
│   └── config.yaml       ← Configuración (editar rutas aquí)
├── graphical_results/    ← Figuras generadas (PNG, SVG)
├── reports/              ← Reporte Markdown generado
├── docs/
│   └── architecture.md   ← Documentación técnica del sistema
├── csv_files/            ← Datos locales (no incluidos en el repo)
├── additional_files/     ← ClinVar .gz (no incluido — demasiado grande)
├── DBA-context/          ← Papers de referencia
├── notebooks/            ← Notebook original (referencia histórica)
├── graph_background/     ← Código fuente de Plotly (dependencia)
├── requirements.txt
├── .gitignore
└── LICENSE
```

## Resultados

- **Figuras** → `graphical_results/` (PNG de alta resolución + SVG vectorial)
  - `{GEN}_lollipop.png` — distribución posicional de variantes locales
  - `{GEN}_kde.png` o `{GEN}_histogram.png` — densidad de variantes ClinVar
  - `global_multi_overview.png` — todos los genes en un único panel
  - `global_heatmap.png`, `global_violin.png`, `global_kde_multi.png`

- **Reporte** → `reports/dba_variant_analysis_report.md`
  Documento con introducción, métodos, resultados tabulados, figuras y discusión.

- **Log** → `pipeline.log`
  Registro completo con tiempos, decisiones de visualización justificadas y advertencias.

## Configuración

Edita `scripts/config.yaml` para adaptar el pipeline a tu entorno:

```yaml
paths:
  variants_csv:      "../csv_files/TABLA_VARIANTES_IBMFS_DBA1.csv"
  clinvar_gz:        "../additional_files/variant_summary.txt.gz"
  graphical_results: "../graphical_results"
  reports:           "../reports"

entrez:
  email: "tu.email@ejemplo.com"   # requerido por NCBI

analysis:
  hotspot_window: 10    # aminoácidos
  alpha: 0.05
```

Para añadir un gen, agrega una entrada en la sección `genes` del YAML — no se
requieren cambios en el código.

## Notas sobre exportación de figuras

- Se requiere [`kaleido`](https://github.com/plotly/Kaleido) para exportar PNG/SVG/PDF.
- Sin kaleido, el pipeline hace fallback automático a HTML interactivo y notifica en el log.
- Para publicación científica se recomienda SVG (vectorial, sin pérdida de calidad).

## Licencia

MIT — ver [LICENSE](LICENSE).
