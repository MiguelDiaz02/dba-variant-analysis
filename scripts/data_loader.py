"""
data_loader.py
==============
Responsabilidad única: ingesta, limpieza y normalización de los datos de variantes DBA.

Funciones exportadas:
    load_variants(filepath)          -> pd.DataFrame
    prepare_dataframe(df, valid_acmg) -> pd.DataFrame
    extract_position(aa_change)       -> int | None
    classify_variant(protein_change)  -> str
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Optional

import pandas as pd

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constante compartida: regex único para extraer posición proteica HGVS
# Cubre: p.Ser4Ter, p.Cys74Tyr, p.R119*, (p.Ala25Val), NM_xxx:p.Leu12fs
# ---------------------------------------------------------------------------
_HGVS_POSITION_RE = re.compile(r"p\.[A-Za-z]{1,3}(\d+)", re.IGNORECASE)

# Candidatos de columnas por concepto (búsqueda case-insensitive + parcial)
_COLUMN_ALIASES: dict[str, list[str]] = {
    "Gene":     ["gen", "gene"],
    "ACMG":     ["acmg"],
    "AA_change": ["cambio en aminoácido", "aminoacido", "protein_change", "aa_change", "aa"],
    "Reported": ["reportada en la literatura", "reportada", "literatura", "reported"],
    "cDNA":     ["cambio a nivel cdna", "cdna", "hgvs_c"],
}


# ---------------------------------------------------------------------------
# Carga y limpieza
# ---------------------------------------------------------------------------

def load_variants(filepath: Path) -> pd.DataFrame:
    """
    Lee el CSV de variantes DBA y devuelve un DataFrame con columnas normalizadas.

    El archivo tiene dos filas de cabecera que se fusionan antes de procesar
    los registros. Las columnas vacías y las de código de familia (FM_) se
    eliminan automáticamente.

    Args:
        filepath: Ruta al archivo TABLA_VARIANTES_IBMFS_DBA1.csv.

    Returns:
        DataFrame con columnas estandarizadas: Gene, ACMG, AA_change, Reported.

    Raises:
        FileNotFoundError: si el archivo no existe.
        ValueError: si el archivo tiene menos de 3 líneas (cabecera + datos).
    """
    if not filepath.exists():
        raise FileNotFoundError(
            f"Archivo de variantes no encontrado: {filepath}\n"
            "Verifica la ruta en config.yaml → paths.variants_csv"
        )

    try:
        raw_lines = filepath.read_text(encoding="utf-8-sig").splitlines()
    except UnicodeDecodeError:
        raw_lines = filepath.read_text(encoding="latin-1").splitlines()
        log.warning("Fallback a encoding latin-1 para %s", filepath.name)

    if len(raw_lines) < 3:
        raise ValueError(
            f"{filepath.name} tiene {len(raw_lines)} líneas — "
            "se esperan al menos 3 (2 cabeceras + 1 fila de datos)."
        )

    def _split(line: str) -> list[str]:
        return [c.strip() for c in line.split(",")]

    row0, row1 = _split(raw_lines[0]), _split(raw_lines[1])
    headers = _merge_headers(row0, row1)

    records = []
    for raw in raw_lines[2:]:
        values = _split(raw)
        if not any(v for v in values):
            continue
        records.append(dict(zip(headers, values)))

    df = pd.DataFrame(records)

    # Eliminar columnas vacías y columnas FM_
    df = df.loc[:, df.columns.str.strip().astype(bool)]
    df = df.loc[:, ~df.columns.str.startswith("FM_")]

    # Renombrar columnas a nombres canónicos
    df = _rename_columns(df)

    # Limpieza de Gene
    df["Gene"] = df["Gene"].str.strip()
    df = df[df["Gene"].notna() & (df["Gene"] != "")]

    # Limpieza de ACMG
    if "ACMG" in df.columns:
        df["ACMG"] = df["ACMG"].str.strip()

    log.info(
        "Variantes cargadas: %d filas — genes: %s",
        len(df),
        sorted(df["Gene"].unique().tolist()),
    )
    return df


def prepare_dataframe(df: pd.DataFrame, valid_acmg: list[str]) -> pd.DataFrame:
    """
    Filtra el DataFrame a las clasificaciones ACMG de interés y añade
    las columnas derivadas 'position' y 'variant_type'.

    Args:
        df:          DataFrame resultado de load_variants().
        valid_acmg:  Lista de códigos ACMG a conservar (p. ej. ["P","LP","VUS"]).

    Returns:
        DataFrame filtrado con columnas 'position' (int|None) y 'variant_type' (str).
    """
    out = df[df["ACMG"].isin(valid_acmg)].copy()
    out["position"] = out["AA_change"].apply(extract_position)
    out["variant_type"] = out["AA_change"].apply(classify_variant)
    log.debug("prepare_dataframe: %d filas después de filtro ACMG=%s", len(out), valid_acmg)
    return out


# ---------------------------------------------------------------------------
# Extracción y clasificación de variantes
# ---------------------------------------------------------------------------

def extract_position(aa_change) -> Optional[int]:
    """
    Extrae la posición numérica del residuo de un cambio proteico en notación HGVS.

    Ejemplos:
        "p.Ser4Ter"    → 4
        "p.Cys74Tyr"   → 74
        "(p.Ala12Val)" → 12
        NaN / ""       → None

    Args:
        aa_change: Cadena con cambio proteico HGVS, o NaN/None.

    Returns:
        Número entero de la posición, o None si no se puede extraer.
    """
    if pd.isna(aa_change) or str(aa_change).strip() == "":
        return None
    match = _HGVS_POSITION_RE.search(str(aa_change))
    return int(match.group(1)) if match else None


def classify_variant(protein_change) -> str:
    """
    Clasifica el tipo funcional de una variante a partir del cambio proteico HGVS.

    Criterios (en orden de prioridad):
        - Frameshift: contiene "fs"
        - Nonsense:   contiene "*" (codón stop prematuro)
        - Missense:   empieza con "p." y no contiene "("
        - Splice:     cadena vacía o NaN (sin cambio proteico → variante de splicing)
        - Other:      cualquier otro caso

    Args:
        protein_change: Cadena con cambio proteico, o NaN/None.

    Returns:
        Uno de: "Frameshift", "Nonsense", "Missense", "Splice", "Other".
    """
    val = "" if pd.isna(protein_change) else str(protein_change).strip()
    if "fs" in val:
        return "Frameshift"
    if "*" in val:
        return "Nonsense"
    if val.startswith("p.") and "(" not in val:
        return "Missense"
    if val in ("", "nan"):
        return "Splice"
    return "Other"


# ---------------------------------------------------------------------------
# Utilidades internas
# ---------------------------------------------------------------------------

def _merge_headers(row0: list[str], row1: list[str]) -> list[str]:
    """Fusiona dos filas de cabecera en una sola lista de nombres únicos."""
    headers = []
    seen: dict[str, int] = {}
    for a, b in zip(row0, row1):
        name = f"{a}_{b}".strip("_") if b else a
        name = name.strip() or "col"
        # Deduplica añadiendo sufijo numérico si el nombre ya existe
        if name in seen:
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 0
        headers.append(name)
    return headers


def _find_column(df: pd.DataFrame, aliases: list[str]) -> Optional[str]:
    """
    Busca en las columnas de df alguna que coincida (case-insensitive, parcial)
    con cualquiera de los alias dados.

    Returns:
        El nombre real de la columna, o None si no se encuentra.
    """
    lower_map = {c.lower(): c for c in df.columns}
    for alias in aliases:
        key = alias.lower()
        # Coincidencia exacta primero
        if key in lower_map:
            return lower_map[key]
        # Coincidencia parcial
        for col_lower, col_orig in lower_map.items():
            if key in col_lower:
                return col_orig
    return None


def _rename_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Renombra las columnas del DataFrame a nombres canónicos según _COLUMN_ALIASES."""
    rename_map = {}
    for canonical, aliases in _COLUMN_ALIASES.items():
        found = _find_column(df, aliases)
        if found and found != canonical:
            rename_map[found] = canonical
        elif found is None:
            log.warning(
                "Columna '%s' no encontrada. Alias buscados: %s. "
                "Columnas disponibles: %s",
                canonical, aliases, list(df.columns),
            )
    return df.rename(columns=rename_map)
