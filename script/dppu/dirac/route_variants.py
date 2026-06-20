"""RBD route variant audit: Routh-before-Dirac elimination variants.

Pure compute module — no I/O, no logging, no side effects on import.

RBD route variant audit helpers (math unchanged).
The ``C`` module-object argument previously used to call back into the
DBR constraint-closure module is replaced by direct imports from
dppu.dirac.constraint_closure.
"""

from __future__ import annotations

from typing import Any, Dict, List

import sympy as sp

from dppu.dirac.constraint_closure import (
    sample_geometry,
    build_dirac_for_second_order_lagrangian,
)


def routh_eliminate(local: Dict[str, Any], also_K0: bool) -> Dict[str, Any]:
    """Routh-eliminate the spatial omega (and optionally K0) from the LOCAL Lagrangian.
    Returns L_red plus the dynamical coframe velocity map."""
    N, q = local["N"], local["q"]
    lp = local["lorentz_pairs"]
    omega_list = [local["omega_syms"][(a, b, i)] for a, b in lp for i in range(1, 4)]
    dom_list = [local["dot_omega_pos_syms"][(a, b, i)] for a, b in lp for i in range(1, 4)]
    K0_list = [local["K0_syms"][(a, b)] for a, b in lp]
    L_kin, L_static = local["L_kin_A"], local["L_static_A"]
    s, qA, X = local["s"], local["qA_vars"], local["X_vars"]
    dyn = [s] + list(qA) + list(X)
    ddyn = {v: sp.Symbol(f"d_{v}") for v in dyn}
    sub1 = {N: sp.S.One, q: sp.S.One}
    y = omega_list + (K0_list if also_K0 else [])
    coeff = [sp.diff(L_kin, dom) for dom in dom_list]
    dotC = [sp.expand(sum(sp.diff(c, dv) * ddyn[dv] for dv in dyn)) for c in coeff]
    yzero = {v: sp.S.Zero for v in y}
    M = sp.Matrix([[sp.diff(L_static, y1, y2).subs(yzero) for y2 in y] for y1 in y]).subs(sub1)
    rank = int(M.rank())
    if rank < len(y):
        return {"status": "BLOCK_SINGULAR", "rank": rank, "dim": len(y)}
    b = sp.Matrix([sp.diff(L_static, yv).subs(yzero) for yv in y]).subs(sub1)
    dotC_ext = [sp.expand(c.subs(sub1)) for c in dotC] + ([sp.S.Zero] * len(K0_list) if also_K0 else [])
    u = sp.Matrix([dotC_ext[i] - b[i] for i in range(len(y))])
    L0 = sp.expand(L_static.subs(yzero).subs(sub1))
    L_red = sp.expand(L0 - sp.Rational(1, 2) * (u.T * M.inv() * u)[0])
    return {"status": "REDUCED", "rank": rank, "dim": len(y), "L_red": L_red,
            "dyn": dyn, "ddyn": ddyn, "K0_eliminated": also_K0}


def project_reduced_dirac(local: Dict[str, Any], routh: Dict[str, Any], sample_l: int) -> Dict[str, Any]:
    geom = sample_geometry(sample_l, "x")
    sl, V, H = geom["sqrt_lambda"], geom["V"], geom["H"]
    inc_long, inc_E, inc_X = sample_l >= 1, sample_l >= 2, sample_l >= 1
    keep_K0 = not routh["K0_eliminated"]
    s, qA, X = local["s"], local["qA_vars"], local["X_vars"]
    ddyn = routh["ddyn"]
    lp = local["lorentz_pairs"]
    n_l = sp.Symbol(f"n_l{sample_l}"); s_l = sp.Symbol(f"s_l{sample_l}")
    B_l = sp.Symbol(f"B_l{sample_l}"); E_l = sp.Symbol(f"E_l{sample_l}"); X_l = sp.Symbol(f"X_l{sample_l}")
    ds_l = sp.Symbol(f"ds_l{sample_l}"); dE_l = sp.Symbol(f"dE_l{sample_l}"); dX_l = sp.Symbol(f"dX_l{sample_l}")
    K0_l = {(a, b): sp.Symbol(f"K0_{a}{b}_l{sample_l}") for a, b in lp}
    subs: Dict[sp.Symbol, sp.Expr] = {local["n_pert"]: n_l, s: s_l, ddyn[s]: ds_l}
    qA00, qA11, qA01, qA02, qA12 = qA
    emap = {qA00: H[(0, 0)], qA11: H[(1, 1)], qA01: H[(0, 1)], qA02: H[(0, 2)], qA12: H[(1, 2)]}
    for qv, h in emap.items():
        subs[qv] = E_l * h if inc_E else sp.S.Zero
        subs[ddyn[qv]] = dE_l * h if inc_E else sp.S.Zero
    for idx, xv in enumerate(X):
        subs[xv] = sl * X_l * V[idx] if inc_X else sp.S.Zero
        subs[ddyn[xv]] = sl * dX_l * V[idx] if inc_X else sp.S.Zero
    for j in range(3):
        subs[local["NADM"][j]] = sl * B_l * V[j] if inc_long else sp.S.Zero
    if keep_K0:
        for a, b in lp:
            subs[local["K0_syms"][(a, b)]] = K0_l[(a, b)]
    L_red_p = sp.expand(routh["L_red"].subs(subs))
    q_red: List[sp.Symbol] = [n_l, s_l]
    v_red: List[sp.Symbol] = [sp.Symbol(f"dn_aux_l{sample_l}"), ds_l]
    if inc_long:
        q_red.append(B_l); v_red.append(sp.Symbol(f"dB_aux_l{sample_l}"))
    if inc_X:
        q_red.append(X_l); v_red.append(dX_l)
    if inc_E:
        q_red.append(E_l); v_red.append(dE_l)
    if keep_K0:
        for a, b in lp:
            q_red.append(K0_l[(a, b)]); v_red.append(sp.Symbol(f"dK0_{a}{b}_aux_l{sample_l}"))
    res = build_dirac_for_second_order_lagrangian(L_red_p, q_red, v_red)
    res["N_config"] = res.get("N_config", len(q_red))
    return res
