#!/usr/bin/env python
"""Reproduce the projected spin-1 quotient used in Appendix I.

This script is intentionally limited to the pre-Schur projected quotient.  It
reproduces lambda_proj and verifies that the projected direction is not a
finite generalized eigenvector of the rank-deficient (H, G) pencil.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path
from typing import Any

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
PAPER_ROOT = SCRIPT_DIR.parents[2]
WITNESS_CANDIDATES = [
    SCRIPT_DIR / "spin1_auxiliary_exclusion.py",
    PAPER_ROOT / "script" / "scripts" / "paper03" / "spin1_auxiliary_exclusion.py",
]
WITNESS_SOURCE = next((path for path in WITNESS_CANDIDATES if path.exists()), WITNESS_CANDIDATES[-1])
OUTPUT_DIR = PAPER_ROOT / "data"
SUMMARY_PATH = OUTPUT_DIR / "spin1_appI_projection_artifact_summary.json"


class LiteralEvaluator:
    def __init__(self) -> None:
        self.env: dict[str, Any] = {}

    def eval(self, node: ast.AST) -> Any:
        if isinstance(node, ast.Constant):
            if isinstance(node.value, (int, float)):
                return sp.Rational(str(node.value))
            return node.value
        if isinstance(node, ast.Name):
            return self.env[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -self.eval(node.operand)
        if isinstance(node, ast.BinOp):
            left = self.eval(node.left)
            right = self.eval(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
        if isinstance(node, (ast.List, ast.Tuple)):
            return [self.eval(item) for item in node.elts]
        if isinstance(node, ast.Call):
            name = ast.unparse(node.func)
            if name.endswith("array") and node.args:
                return sp.Matrix(self.eval(node.args[0]))
        raise ValueError(f"unsupported witness expression: {ast.dump(node)}")

    def assign(self, target: ast.expr, value: Any) -> None:
        if isinstance(target, ast.Name):
            self.env[target.id] = value


def load_projection_witness() -> dict[str, Any]:
    wanted = {"G_oo", "G_dd", "G_od", "mu2_oo", "mu2_dd", "m_od", "u_vec"}
    tree = ast.parse(WITNESS_SOURCE.read_text(encoding="utf-8"), filename=str(WITNESS_SOURCE))
    evaluator = LiteralEvaluator()
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id in wanted:
                evaluator.assign(target, evaluator.eval(node.value))
    missing = wanted - set(evaluator.env)
    if missing:
        raise RuntimeError(f"missing witness assignments: {sorted(missing)}")
    env = evaluator.env
    return {
        "G": sp.Matrix([[env["G_oo"], env["G_od"]], [env["G_od"], env["G_dd"]]]),
        "H": sp.Matrix([[env["mu2_oo"], env["m_od"]], [env["m_od"], env["mu2_dd"]]]),
        "u": sp.Matrix(env["u_vec"]),
    }


def as_float(value: object) -> float:
    return float(sp.N(value, 30))


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_projection_witness()
    G: sp.Matrix = data["G"]
    H: sp.Matrix = data["H"]
    u: sp.Matrix = data["u"]

    lambda_proj = sp.simplify((u.T * H * u)[0] / (u.T * G * u)[0])
    det_at_projection = sp.simplify((H - lambda_proj * G).det())
    residual = sp.simplify(H * u - lambda_proj * G * u)
    residual_norm2 = sp.simplify((residual.T * residual)[0])

    checks = {
        "lambda_proj_positive": as_float(lambda_proj) > 0,
        "det_at_projection_nonzero": abs(as_float(det_at_projection)) > 1e-6,
        "projected_direction_not_eigenvector": abs(as_float(residual_norm2)) > 1e-6,
    }
    summary = {
        "script": str(Path(__file__).resolve()),
        "witness_source": str(WITNESS_SOURCE),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "exit_label": "PAPER03_APP_I_PROJECTION_QUOTIENT_REPRODUCED_PASS" if all(checks.values()) else "PAPER03_APP_I_PROJECTION_QUOTIENT_REPRODUCED_FAIL",
        "lambda_proj": f"{as_float(lambda_proj):.12g}",
        "det_H_minus_lambda_proj_G": f"{as_float(det_at_projection):.12g}",
        "residual_norm2": f"{as_float(residual_norm2):.12g}",
        "interpretation": "The projected quotient is reproduced, but the projected direction does not solve the finite generalized eigenvalue equation.",
        "checks": checks,
    }
    SUMMARY_PATH.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if summary["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
