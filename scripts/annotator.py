"""
annotator.py
============
Responsabilidad única: integración con fuentes externas de anotación funcional.

Fuentes soportadas:
    - ClinVar: archivo local variant_summary.txt.gz o API NCBI Entrez
    - UniProt: REST API para dominios, regiones y motivos proteicos

Funciones exportadas:
    load_clinvar_file(path)                         -> pd.DataFrame
    filter_clinvar_pathogenic(clinvar_df)           -> pd.DataFrame
    get_clinvar_positions(source, gene, prot_len, cfg) -> list[int]
    get_uniprot_features(uniprot_id)                -> list[dict]
    find_overlapping_domains(hotspots, features)    -> list[dict]
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

from data_loader import _HGVS_POSITION_RE, extract_position

log = logging.getLogger(__name__)

_UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"


# ---------------------------------------------------------------------------
# ClinVar — archivo local
# ---------------------------------------------------------------------------

def load_clinvar_file(path: Path) -> pd.DataFrame:
    """
    Carga el archivo variant_summary.txt.gz descargado desde ClinVar FTP.

    El archivo se lee completo en memoria. Para proteínas típicas en DBA,
    el filtrado posterior reduce el tamaño a miles de filas.

    Args:
        path: Ruta al archivo variant_summary.txt.gz.

    Returns:
        DataFrame completo del archivo ClinVar.

    Raises:
        FileNotFoundError: si el archivo no existe, con instrucciones de descarga.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Archivo ClinVar no encontrado: {path}\n"
            "Descarga desde:\n"
            "  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        )
    log.info("Cargando ClinVar desde %s (archivo grande — puede tardar ~30s)…", path.name)
    try:
        df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    except Exception as exc:
        raise RuntimeError(f"Error al leer {path.name}: {exc}") from exc
    log.info("ClinVar cargado: %d variantes, %d columnas", len(df), df.shape[1])
    return df


def filter_clinvar_pathogenic(clinvar_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filtra el DataFrame de ClinVar a variantes con clasificación patogénica.

    Conserva registros cuyo campo 'ClinicalSignificance' contiene la cadena
    'Pathogenic' (case-insensitive), incluyendo 'Likely pathogenic' y
    'Pathogenic/Likely pathogenic'.

    Args:
        clinvar_df: DataFrame completo de ClinVar.

    Returns:
        Subconjunto filtrado.
    """
    mask = clinvar_df["ClinicalSignificance"].str.contains(
        "Pathogenic", case=False, na=False
    )
    filtered = clinvar_df[mask].copy()
    log.info("ClinVar patogénicas: %d variantes (de %d totales)", len(filtered), len(clinvar_df))
    return filtered


# ---------------------------------------------------------------------------
# ClinVar — API NCBI Entrez
# ---------------------------------------------------------------------------

def _fetch_clinvar_via_api(
    gene: str,
    protein_length: int,
    email: str,
    max_results: int = 500,
    sleep_s: float = 0.15,
) -> list[int]:
    """
    Consulta la API de NCBI Entrez ClinVar para obtener posiciones proteicas.

    Hace una búsqueda esearch para obtener IDs, luego efetch por cada ID
    y extrae posiciones mediante regex HGVS del XML de respuesta.

    Args:
        gene:           Nombre del gen (p. ej. "RPS19").
        protein_length: Longitud proteica para filtrar posiciones válidas.
        email:          Email requerido por NCBI para uso de la API.
        max_results:    Máximo de registros ClinVar a recuperar.
        sleep_s:        Segundos de espera entre peticiones (respeta rate limit).

    Returns:
        Lista de posiciones proteicas enteras en rango [1, protein_length].
    """
    try:
        from Bio import Entrez  # import local para no romper si BioPython no está instalado
    except ImportError:
        log.error("BioPython no instalado. Instala con: pip install biopython")
        return []

    Entrez.email = email
    log.info("ClinVar API → %s (max %d resultados)…", gene, max_results)

    try:
        handle = Entrez.esearch(db="clinvar", term=f"{gene}[gene]", retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
    except Exception as exc:
        log.error("Entrez.esearch falló para %s: %s", gene, exc)
        return []

    ids = record.get("IdList", [])
    log.debug("%s: %d IDs obtenidos de ClinVar", gene, len(ids))

    positions: list[int] = []
    for clinvar_id in ids:
        time.sleep(sleep_s)
        try:
            fh = Entrez.efetch(db="clinvar", id=clinvar_id, rettype="vcv", retmode="xml")
            content = fh.read().decode("utf-8", errors="ignore")
            fh.close()
            for match in _HGVS_POSITION_RE.finditer(content):
                pos = int(match.group(1))
                if 1 <= pos <= protein_length:
                    positions.append(pos)
        except Exception as exc:
            log.warning("Error en efetch para ClinVar ID %s: %s", clinvar_id, exc)

    log.info("%s: %d posiciones recuperadas vía API Entrez", gene, len(positions))
    return positions


def _extract_positions_from_clinvar_df(
    clinvar_path_df: pd.DataFrame,
    gene: str,
    protein_length: int,
) -> list[int]:
    """
    Extrae posiciones proteicas de un DataFrame de ClinVar ya filtrado.

    Args:
        clinvar_path_df: ClinVar filtrado a variantes patogénicas.
        gene:            Nombre del gen.
        protein_length:  Longitud de la proteína para filtrar posiciones.

    Returns:
        Lista de posiciones enteras válidas.
    """
    gene_df = clinvar_path_df[
        clinvar_path_df["GeneSymbol"].str.upper() == gene.upper()
    ]
    positions: list[int] = []
    for name in gene_df["Name"].dropna():
        pos = extract_position(str(name))
        if pos is not None and 1 <= pos <= protein_length:
            positions.append(pos)
    log.info("%s: %d posiciones extraídas del archivo ClinVar local", gene, len(positions))
    return positions


def get_clinvar_positions(
    gene: str,
    protein_length: int,
    cfg: dict,
    clinvar_path_df: Optional[pd.DataFrame] = None,
) -> list[int]:
    """
    Punto de entrada unificado para obtener posiciones ClinVar patogénicas.

    Si `clinvar_path_df` está disponible, usa el archivo local.
    Si no, recurre a la API de NCBI Entrez.

    Args:
        gene:            Nombre del gen.
        protein_length:  Longitud de la proteína.
        cfg:             Configuración completa (bloque 'entrez').
        clinvar_path_df: DataFrame ClinVar pre-filtrado (opcional).

    Returns:
        Lista de posiciones proteicas enteras.
    """
    if clinvar_path_df is not None:
        return _extract_positions_from_clinvar_df(clinvar_path_df, gene, protein_length)

    entrez = cfg.get("entrez", {})
    return _fetch_clinvar_via_api(
        gene=gene,
        protein_length=protein_length,
        email=entrez.get("email", ""),
        max_results=entrez.get("max_results", 500),
        sleep_s=entrez.get("sleep_between_requests", 0.15),
    )


# ---------------------------------------------------------------------------
# UniProt — anotación funcional
# ---------------------------------------------------------------------------

def get_uniprot_features(uniprot_id: str) -> list[dict]:
    """
    Obtiene dominios, regiones y motivos de una proteína desde UniProt REST API.

    Consulta el endpoint JSON de UniProt Knowledge Base y extrae las anotaciones
    estructurales y funcionales con sus coordenadas en la secuencia.

    Args:
        uniprot_id: Identificador UniProt (p. ej. "P39019" para RPS19 humana).

    Returns:
        Lista de dicts con claves: type, description, start, end.
        Lista vacía si la petición falla o el ID no existe.
    """
    if not uniprot_id:
        return []

    url = f"{_UNIPROT_BASE}/{uniprot_id}.json"
    try:
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
    except requests.exceptions.HTTPError as exc:
        log.warning("UniProt HTTP error para %s: %s", uniprot_id, exc)
        return []
    except requests.exceptions.RequestException as exc:
        log.warning("UniProt request failed para %s: %s", uniprot_id, exc)
        return []

    features: list[dict] = []
    for feat in resp.json().get("features", []):
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end   = loc.get("end",   {}).get("value")
        if start is not None and end is not None:
            features.append({
                "type":        feat.get("type", ""),
                "description": feat.get("description", ""),
                "start":       int(start),
                "end":         int(end),
            })

    log.debug("UniProt %s: %d features", uniprot_id, len(features))
    return features


# ---------------------------------------------------------------------------
# Solapamiento hotspot ↔ dominio
# ---------------------------------------------------------------------------

def find_overlapping_domains(
    hotspots: list[dict],
    features: list[dict],
) -> list[dict]:
    """
    Identifica dominios UniProt que se solapan con al menos un hotspot detectado.

    Dos regiones se solapan si comparten al menos un aminoácido:
        max(start1, start2) <= min(end1, end2)

    Args:
        hotspots: Lista de hotspots (con claves 'start' y 'end').
        features: Lista de features UniProt (con claves 'start' y 'end').

    Returns:
        Subconjunto de features que solapan con algún hotspot.
    """
    if not hotspots or not features:
        return []

    overlapping = []
    for feat in features:
        for hs in hotspots:
            if max(feat["start"], hs["start"]) <= min(feat["end"], hs["end"]):
                overlapping.append(feat)
                break  # ya solapó con al menos un hotspot

    log.debug(
        "find_overlapping_domains: %d features, %d solapan con hotspots",
        len(features), len(overlapping),
    )
    return overlapping
