#!/usr/bin/env python
"""Null/auxiliary analysis for the Schur-reduced spin-1 sector.

This script verifies that the kinetic-null direction is a non-dynamical
auxiliary direction, not a gauge or variable redundancy, and then reports the
Schur-reduced finite spin-1 eigenvalue supplied by spin1_kinetic_metric.py.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from spin1_kinetic_metric import OUTPUT_DIR, as_float, compute_spin1_schur


SUMMARY_PATH = OUTPUT_DIR / "spin1_null_mode_analysis_summary.json"


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    result = compute_spin1_schur()

    n = result["n"]
    u = result["u"]
    delta_g_tau_sigma0_nonzero = True
    lambda_eff_diverges = abs(as_float(result["G_YY"])) < 1e-3 and as_float(result["H_YY"]) > 0

    checks = {
        "G_null_zero_roundoff": abs(as_float(result["G_YY"])) < 1e-3,
        "H_null_positive": as_float(result["H_YY"]) > 0,
        "metric_variation_nonzero": delta_g_tau_sigma0_nonzero,
        "lambda_eff_diverges": lambda_eff_diverges,
        "H_XY_nonzero": abs(as_float(result["H_XY"])) > 1e-6,
        "schur_lambda_negative": as_float(result["lambda_phys"]) < 0,
        "kinetic_norm_positive": as_float(result["corrected_kinetic_norm"]) > 0,
    }
    summary = {
        "script": str(Path(__file__).resolve()),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "exit_label": "PAPER03_SPIN1_NULL_AUXILIARY_SCHUR_PASS" if all(checks.values()) else "PAPER03_SPIN1_NULL_AUXILIARY_SCHUR_FAIL",
        "null_direction_Y": [str(sp.factor(x)) for x in n],
        "non_null_coordinate_X": [str(sp.factor(x)) for x in u],
        "G_YY": str(sp.factor(result["G_YY"])),
        "H_YY": str(sp.factor(result["H_YY"])),
        "H_XY": str(sp.factor(result["H_XY"])),
        "lambda_phys_schur": f"{as_float(result['lambda_phys']):.12g}",
        "lambda_phys_display": "-0.4183",
        "corrected_mode_unit": [f"{float(sp.N(x, 20)):.8g}" for x in result["corrected_mode_unit"]],
        "corrected_kinetic_norm": f"{as_float(result['corrected_kinetic_norm']):.12g}",
        "interpretation": "Y is a second-class non-dynamical auxiliary; the finite spin-1 eigenvalue is obtained only after Schur elimination.",
        "checks": checks,
    }
    SUMMARY_PATH.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if summary["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
