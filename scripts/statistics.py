"""
statistics.py
=============
Responsabilidad única: análisis estadístico de hotspots de mutación.

Implementa detección de regiones de mutación recurrente mediante ventana
deslizante con test binomial (scipy.stats.binomtest) y corrección de
Bonferroni para comparaciones múltiples.

Funciones exportadas:
    find_hotspots(positions, protein_length, window, alpha) -> list[dict]
    merge_overlapping_hotspots(hotspots)                    -> list[dict]
    summarize_hotspots(gene, hotspots, protein_length)      -> dict
"""

from __future__ import annotations

import logging
from typing import TypedDict

from scipy.stats import binomtest

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Tipos
# ---------------------------------------------------------------------------

class Hotspot(TypedDict):
    """Región de la proteína con enriquecimiento estadísticamente significativo."""
    start:    int    # posición inicial (inclusive)
    end:      int    # posición final (inclusive)
    count:    int    # número de variantes observadas en la ventana
    expected: float  # número esperado bajo distribución uniforme
    p_value:  float  # p-valor corregido por Bonferroni


# ---------------------------------------------------------------------------
# Detección de hotspots
# ---------------------------------------------------------------------------

def find_hotspots(
    positions: list[int],
    protein_length: int,
    window: int = 10,
    alpha: float = 0.05,
) -> list[Hotspot]:
    """
    Detecta hotspots de mutación mediante ventana deslizante y test binomial
    con corrección de Bonferroni.

    Hipótesis nula: las variantes se distribuyen uniformemente a lo largo de
    la proteína, con probabilidad p = window / protein_length de caer en
    cualquier ventana de tamaño `window`.

    Se realizan (protein_length - window + 1) tests en total. El umbral de
    significancia se corrige dividiendo alpha por ese número de tests
    (corrección de Bonferroni conservadora).

    Args:
        positions:      Lista de posiciones proteicas observadas (enteros ≥1).
        protein_length: Longitud total de la proteína en aminoácidos.
        window:         Tamaño de la ventana deslizante en aminoácidos.
        alpha:          Nivel de significancia nominal (antes de corrección).

    Returns:
        Lista de Hotspot ordenada por posición de inicio, con regiones
        solapadas ya fusionadas. Lista vacía si no se detectan hotspots.

    Raises:
        ValueError: si protein_length < window, o si window < 1.
    """
    if window < 1:
        raise ValueError(f"window debe ser ≥ 1, recibido: {window}")
    if protein_length < window:
        raise ValueError(
            f"protein_length ({protein_length}) debe ser ≥ window ({window})"
        )
    if not positions:
        log.debug("find_hotspots: lista de posiciones vacía → sin hotspots")
        return []

    n_total = len(positions)
    p_uniform = window / protein_length
    n_tests = protein_length - window + 1
    alpha_corrected = alpha / n_tests

    log.debug(
        "find_hotspots: n=%d, prot_len=%d, window=%d, α=%.4f, α_Bonf=%.2e",
        n_total, protein_length, window, alpha, alpha_corrected,
    )

    raw_hotspots: list[Hotspot] = []
    positions_set = positions  # list para contar con sum()

    for start in range(1, n_tests + 1):
        end = start + window - 1
        count = sum(start <= pos <= end for pos in positions_set)
        if count == 0:
            continue

        result = binomtest(count, n_total, p_uniform, alternative="greater")
        if result.pvalue < alpha_corrected:
            raw_hotspots.append(
                Hotspot(
                    start=start,
                    end=end,
                    count=count,
                    expected=round(n_total * p_uniform, 3),
                    p_value=result.pvalue,
                )
            )

    merged = merge_overlapping_hotspots(raw_hotspots)
    log.debug(
        "find_hotspots: %d ventanas significativas → %d hotspot(s) después de fusión",
        len(raw_hotspots), len(merged),
    )
    return merged


def merge_overlapping_hotspots(hotspots: list[Hotspot]) -> list[Hotspot]:
    """
    Fusiona hotspots solapados o adyacentes en regiones continuas únicas.

    Cuando dos hotspots se fusionan, el resultado conserva:
        - start: el menor de los dos
        - end:   el mayor de los dos
        - count: el máximo (ventana más enriquecida)
        - p_value: el mínimo (ventana más significativa)

    Args:
        hotspots: Lista de Hotspot, en cualquier orden.

    Returns:
        Lista de Hotspot fusionados, ordenada por 'start'.
    """
    if not hotspots:
        return []

    sorted_hs = sorted(hotspots, key=lambda h: h["start"])
    merged: list[Hotspot] = [dict(sorted_hs[0])]  # type: ignore[assignment]

    for current in sorted_hs[1:]:
        last = merged[-1]
        if current["start"] <= last["end"]:
            last["end"]     = max(last["end"], current["end"])
            last["count"]   = max(last["count"], current["count"])
            last["p_value"] = min(last["p_value"], current["p_value"])
        else:
            merged.append(dict(current))  # type: ignore[arg-type]

    return merged  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# Resumen estadístico
# ---------------------------------------------------------------------------

def summarize_hotspots(
    gene: str,
    hotspots: list[Hotspot],
    protein_length: int,
) -> dict:
    """
    Genera un resumen estadístico de los hotspots detectados para un gen.

    Args:
        gene:           Nombre del gen.
        hotspots:       Lista de hotspots detectados.
        protein_length: Longitud de la proteína.

    Returns:
        Dict con claves: gene, n_hotspots, total_aa_affected, hotspot_coverage_pct,
        most_significant (p_value más bajo), hotspot_regions (lista de strings).
    """
    if not hotspots:
        return {
            "gene": gene,
            "n_hotspots": 0,
            "total_aa_affected": 0,
            "hotspot_coverage_pct": 0.0,
            "most_significant_p": None,
            "hotspot_regions": [],
        }

    total_aa = sum(h["end"] - h["start"] + 1 for h in hotspots)
    best_p = min(h["p_value"] for h in hotspots)
    regions = [f"aa {h['start']}–{h['end']} (n={h['count']})" for h in hotspots]

    return {
        "gene": gene,
        "n_hotspots": len(hotspots),
        "total_aa_affected": total_aa,
        "hotspot_coverage_pct": round(100 * total_aa / protein_length, 1),
        "most_significant_p": best_p,
        "hotspot_regions": regions,
    }
