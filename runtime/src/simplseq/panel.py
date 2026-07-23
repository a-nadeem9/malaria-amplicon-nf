"""Describe recovered loci relative to the published SIMPLseq panel."""

from __future__ import annotations

import re
from typing import Iterable


SIMPLSEQ_PANEL = {
    "CSP": {"csp", "pf3d70304600", "0304600"},
    "TRAP": {"trap", "pf3d71335900", "1335900"},
    "WDCP": {"wdcp", "pf3d71410300", "1410300"},
    "KELT": {"kelt", "pf3d71475900", "1475900"},
    "SERA8": {"sera8", "pf3d70207300", "0207300"},
    "SURFIN4.2": {"surfin42", "surfin4.2", "pf3d70424400", "0424400"},
}


def normalize_locus(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").strip().lower())


def canonical_simplseq_locus(value: object) -> str:
    token = normalize_locus(value)
    for canonical, aliases in SIMPLSEQ_PANEL.items():
        if token in {normalize_locus(alias) for alias in aliases}:
            return canonical
    return ""


def panel_profile(loci: Iterable[object]) -> dict[str, object]:
    raw_loci = sorted({str(value).strip() for value in loci if str(value).strip()})
    detected_core: set[str] = set()
    additional: list[str] = []
    for locus in raw_loci:
        canonical = canonical_simplseq_locus(locus)
        if canonical:
            detected_core.add(canonical)
        else:
            additional.append(locus)
    ordered_core = [locus for locus in SIMPLSEQ_PANEL if locus in detected_core]
    missing_core = [locus for locus in SIMPLSEQ_PANEL if locus not in detected_core]
    return {
        "schema_version": 1,
        "reference_panel": "SIMPLseq six-locus panel",
        "reference_loci": list(SIMPLSEQ_PANEL),
        "detected_reference_loci": ordered_core,
        "missing_reference_loci": missing_core,
        "additional_loci": additional,
        "detected_reference_count": len(ordered_core),
        "reference_locus_count": len(SIMPLSEQ_PANEL),
        "all_detected_loci": raw_loci,
        "is_complete_simplseq_panel": len(ordered_core) == len(SIMPLSEQ_PANEL),
        "kelt_sentinel_detected": "KELT" in detected_core,
        "interpretation": (
            "SIMPLseq panel or subset with additional loci"
            if ordered_core and additional
            else "SIMPLseq panel or subset"
            if ordered_core
            else "Custom or legacy panel"
        ),
    }
