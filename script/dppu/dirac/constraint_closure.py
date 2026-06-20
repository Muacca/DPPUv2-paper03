"""Dirac constraint-closure engine for Branch-A first-order Palatini (DBR route).

Pure compute module — no I/O, no logging, no side effects on import.

DBR constraint-closure engine for the Branch-A first-order EC+NY system.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any, Dict, List, Sequence, Tuple

import sympy as sp


def eps3(i: int, j: int, k: int) -> int:
    arr = [i, j, k]
    if sorted(arr) != [1, 2, 3]:
        return 0
    sign = 1
    for p in range(3):
        while arr[p] != p + 1:
            q_ = arr[p] - 1
            arr[p], arr[q_] = arr[q_], arr[p]
            sign = -sign
    return sign


def lc4(a: int, b: int, c: int, d: int) -> int:
    return int(sp.LeviCivita(a, b, c, d))


def build_L_conn_alg(
    N0_sym: sp.Expr,
    q0_sym: sp.Expr,
    K_syms: Dict[Tuple[int, int, int], sp.Symbol],
) -> sp.Expr:
    eta_diag = {0: -1, 1: 1, 2: 1, 3: 1}

    def e_bg_local(a: int, mu: int) -> sp.Expr:
        if a == 0 and mu == 0:
            return N0_sym
        if a == mu and a > 0:
            return q0_sym
        return sp.S.Zero

    def K_up(c: int, e: int, rho: int) -> sp.Expr:
        if c == e:
            return sp.S.Zero
        key = (min(c, e), max(c, e), rho)
        sign = 1 if c < e else -1
        return sign * K_syms.get(key, sp.S.Zero)

    L = sp.S.Zero
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    eps_abcd = lc4(a, b, c, d)
                    if eps_abcd == 0:
                        continue
                    for e in range(4):
                        eta_ee = eta_diag[e]
                        for mu in range(4):
                            ea = e_bg_local(a, mu)
                            if ea == sp.S.Zero:
                                continue
                            for nu in range(4):
                                eb = e_bg_local(b, nu)
                                if eb == sp.S.Zero:
                                    continue
                                for rho in range(4):
                                    Kce = K_up(c, e, rho)
                                    if Kce == sp.S.Zero:
                                        continue
                                    for sigma in range(4):
                                        Ked = K_up(e, d, sigma)
                                        if Ked == sp.S.Zero:
                                            continue
                                        eps_pos = lc4(mu, nu, rho, sigma)
                                        if eps_pos == 0:
                                            continue
                                        L += eps_abcd * eta_ee * eps_pos * ea * eb * Kce * Ked
    return sp.expand(L)


def build_branchA_local_L2() -> Dict[str, Any]:
    N, q = sp.symbols("N q", positive=True)

    def e_bg(a: int, mu: int) -> sp.Expr:
        if a == 0 and mu == 0:
            return N
        if a == mu and a > 0:
            return q
        return sp.S.Zero

    s = sp.Symbol("s")
    n_pert = sp.Symbol("n_pert")
    NADM = [sp.Symbol(f"NADM_{i}") for i in range(1, 4)]
    qA00, qA11, qA01, qA02, qA12 = sp.symbols("qA_00 qA_11 qA_01 qA_02 qA_12")
    qA22 = -(qA00 + qA11)
    qA_comp = {
        (0, 0): qA00,
        (1, 1): qA11,
        (2, 2): qA22,
        (0, 1): qA01,
        (1, 0): qA01,
        (0, 2): qA02,
        (2, 0): qA02,
        (1, 2): qA12,
        (2, 1): qA12,
    }
    X01, X02, X12 = sp.symbols("X_01 X_02 X_12")
    X_comp = {
        (0, 1): X01,
        (1, 0): -X01,
        (0, 2): X02,
        (2, 0): -X02,
        (1, 2): X12,
        (2, 1): -X12,
        (0, 0): sp.S.Zero,
        (1, 1): sp.S.Zero,
        (2, 2): sp.S.Zero,
    }

    lorentz_pairs = list(combinations(range(4), 2))
    K0_syms = {(a, b): sp.Symbol(f"K0_{a}{b}") for a, b in lorentz_pairs}
    omega_syms = {
        (a, b, i): sp.Symbol(f"om_{a}{b}_{i}")
        for a, b in lorentz_pairs
        for i in range(1, 4)
    }
    dot_omega_pos_syms = {
        (a, b, i): sp.Symbol(f"domega_{a}{b}_{i}")
        for a, b in lorentz_pairs
        for i in range(1, 4)
    }
    dot_omega: Dict[Tuple[int, int, int], sp.Expr] = {}
    for a in range(4):
        for b in range(4):
            for i in range(1, 4):
                if a < b:
                    dot_omega[(a, b, i)] = dot_omega_pos_syms[(a, b, i)]
                elif a > b:
                    dot_omega[(a, b, i)] = -dot_omega_pos_syms[(b, a, i)]
                else:
                    dot_omega[(a, b, i)] = sp.S.Zero

    de_A: Dict[Tuple[int, int], sp.Expr] = {}
    for j in range(1, 4):
        de_A[(0, j)] = sp.S.Zero
    for a in range(1, 4):
        for j in range(1, 4):
            ia, ij = a - 1, j - 1
            val = qA_comp.get((ia, ij), sp.S.Zero) + X_comp.get((ia, ij), sp.S.Zero)
            if ia == ij:
                val += s / 3
            de_A[(a, j)] = val

    de_ext_A = dict(de_A)
    de_ext_A[(0, 0)] = n_pert
    for j in range(1, 4):
        de_ext_A[(j, 0)] = q * NADM[j - 1]

    def omega_lc(c: int, e: int, rho: int) -> sp.Expr:
        if c == 0 or e == 0 or rho == 0:
            return sp.S.Zero
        return sp.Rational(1, 1) / q * eps3(c, e, rho)

    def K_sp_up(e: int, d: int, sigma: int) -> sp.Expr:
        if e == d or sigma == 0:
            return sp.S.Zero
        if e < d:
            return omega_syms.get((e, d, sigma), sp.S.Zero)
        return -omega_syms.get((d, e, sigma), sp.S.Zero)

    L_kin_A = sp.S.Zero
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    lc = lc4(a, b, c, d)
                    if lc == 0:
                        continue
                    for i in range(1, 4):
                        for j in range(1, 4):
                            for k in range(1, 4):
                                eijk = eps3(i, j, k)
                                if eijk == 0:
                                    continue
                                de_val = de_A.get((d, k), sp.S.Zero)
                                if de_val == sp.S.Zero:
                                    continue
                                L_kin_A += -2 * lc * eijk * e_bg(c, j) * de_val * dot_omega[(a, b, i)]
    L_kin_A = sp.expand(L_kin_A)

    K_build = {
        (a, b, mu): sp.Symbol(f"K_build_{a}{b}_{mu}")
        for a, b in lorentz_pairs
        for mu in range(4)
    }
    N0_b = sp.Symbol("N0_build", positive=True)
    q0_b = sp.Symbol("q0_build", positive=True)
    L_conn_alg_raw = build_L_conn_alg(N0_b, q0_b, K_build)
    K_subs: Dict[sp.Symbol, sp.Expr] = {N0_b: N + n_pert, q0_b: q}
    for a, b in lorentz_pairs:
        K_subs[K_build[(a, b, 0)]] = K0_syms[(a, b)]
        for i in range(1, 4):
            K_subs[K_build[(a, b, i)]] = omega_syms[(a, b, i)]
    L_conn_alg = sp.expand(L_conn_alg_raw.subs(K_subs))

    L_wK_A = sp.S.Zero
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    eps_abcd = lc4(a, b, c, d)
                    if eps_abcd == 0:
                        continue
                    for e in range(4):
                        for mu in range(4):
                            dea = de_ext_A.get((a, mu), sp.S.Zero)
                            if dea == sp.S.Zero:
                                continue
                            for nu in range(4):
                                ebg_b = e_bg(b, nu)
                                if ebg_b == sp.S.Zero:
                                    continue
                                for rho in range(4):
                                    olc = omega_lc(c, e, rho)
                                    if olc == sp.S.Zero:
                                        continue
                                    for sigma in range(4):
                                        Ksp = K_sp_up(e, d, sigma)
                                        if Ksp == sp.S.Zero:
                                            continue
                                        eps_pos = lc4(mu, nu, rho, sigma)
                                        if eps_pos == 0:
                                            continue
                                        L_wK_A += 2 * eps_abcd * eps_pos * dea * ebg_b * olc * Ksp
    L_wK_A = sp.expand(L_wK_A)

    L_de2_A = sp.S.Zero
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    eps_abcd = lc4(a, b, c, d)
                    if eps_abcd == 0:
                        continue
                    for e in range(4):
                        for mu in range(4):
                            dea = de_ext_A.get((a, mu), sp.S.Zero)
                            if dea == sp.S.Zero:
                                continue
                            for nu in range(4):
                                deb = de_ext_A.get((b, nu), sp.S.Zero)
                                if deb == sp.S.Zero:
                                    continue
                                for rho in range(4):
                                    olc1 = omega_lc(c, e, rho)
                                    if olc1 == sp.S.Zero:
                                        continue
                                    for sigma in range(4):
                                        olc2 = omega_lc(e, d, sigma)
                                        if olc2 == sp.S.Zero:
                                            continue
                                        eps_pos = lc4(mu, nu, rho, sigma)
                                        if eps_pos == 0:
                                            continue
                                        L_de2_A += eps_abcd * eps_pos * dea * deb * olc1 * olc2
    L_de2_A = sp.expand(L_de2_A)

    return {
        "N": N,
        "q": q,
        "n_pert": n_pert,
        "NADM": NADM,
        "s": s,
        "qA_vars": [qA00, qA11, qA01, qA02, qA12],
        "X_vars": [X01, X02, X12],
        "lorentz_pairs": lorentz_pairs,
        "K0_syms": K0_syms,
        "omega_syms": omega_syms,
        "dot_omega_pos_syms": dot_omega_pos_syms,
        "L_kin_A": L_kin_A,
        "L_static_A": sp.expand(L_conn_alg + L_wK_A + L_de2_A),
        "L_total_A": sp.expand(L_kin_A + L_conn_alg + L_wK_A + L_de2_A),
    }


def poisson(f: sp.Expr, g: sp.Expr, q_vars: Sequence[sp.Symbol], p_vars: Sequence[sp.Symbol]) -> sp.Expr:
    result = sp.S.Zero
    for Q, P in zip(q_vars, p_vars):
        result += sp.diff(f, Q) * sp.diff(g, P) - sp.diff(f, P) * sp.diff(g, Q)
    return sp.simplify(result)


def coeff_vector(expr: sp.Expr, vars_all: Sequence[sp.Symbol]) -> List[sp.Expr]:
    expr = sp.expand(expr)
    return [sp.simplify(sp.diff(expr, var)) for var in vars_all]


def is_independent(expr: sp.Expr, existing: List[sp.Expr], vars_all: Sequence[sp.Symbol]) -> bool:
    expr = sp.expand(expr)
    if expr == sp.S.Zero:
        return False
    rows_old = [coeff_vector(item, vars_all) for item in existing]
    rows_new = rows_old + [coeff_vector(expr, vars_all)]
    rank_old = sp.Matrix(rows_old).rank() if rows_old else 0
    rank_new = sp.Matrix(rows_new).rank()
    return rank_new > rank_old


def harmonic_lambda(sample_l: int) -> sp.Expr:
    return sp.Integer(sample_l * (sample_l + 2))


def sample_geometry(sample_l: int, variant: str = "x") -> Dict[str, Any]:
    lam = harmonic_lambda(sample_l)
    sqrt_lam = sp.sqrt(lam)
    if variant == "x":
        V = [sp.S.One, sp.S.Zero, sp.S.Zero]
        H_diag = [sp.sqrt(sp.Rational(1, 6)), sp.sqrt(sp.Rational(1, 6)), -2 * sp.sqrt(sp.Rational(1, 6))]
    elif variant == "y":
        V = [sp.S.Zero, sp.S.One, sp.S.Zero]
        H_diag = [sp.sqrt(sp.Rational(1, 6)), -2 * sp.sqrt(sp.Rational(1, 6)), sp.sqrt(sp.Rational(1, 6))]
    elif variant == "diag":
        V = [sp.sqrt(sp.Rational(1, 3)), sp.sqrt(sp.Rational(1, 3)), sp.sqrt(sp.Rational(1, 3))]
        H_diag = [sp.sqrt(sp.Rational(1, 6)), sp.sqrt(sp.Rational(1, 6)), -2 * sp.sqrt(sp.Rational(1, 6))]
    else:
        raise ValueError(f"unknown projection variant: {variant}")
    H = {
        (0, 0): H_diag[0],
        (1, 1): H_diag[1],
        (2, 2): H_diag[2],
        (0, 1): sp.S.Zero,
        (1, 0): sp.S.Zero,
        (0, 2): sp.S.Zero,
        (2, 0): sp.S.Zero,
        (1, 2): sp.S.Zero,
        (2, 1): sp.S.Zero,
    }
    return {
        "lambda": lam,
        "sqrt_lambda": sqrt_lam,
        "V": V,
        "H": H,
        "variant": variant,
        "V_norm": sp.simplify(sum(item * item for item in V)),
        "H_trace": sp.simplify(H[(0, 0)] + H[(1, 1)] + H[(2, 2)]),
        "H_norm": sp.simplify(sum(H[(i, j)] * H[(i, j)] for i in range(3) for j in range(3))),
        "vector_projection": f"variant={variant}; unit representative for normalized D_iY_l/sqrt(lambda_l)",
        "tensor_projection": f"variant={variant}; unit traceless representative for normalized scalar Hessian",
    }


def project_lagrangian(local: Dict[str, Any], sample_l: int, l_sector: str, projection_variant: str = "x") -> Dict[str, Any]:
    geom = sample_geometry(sample_l, projection_variant)
    lam = geom["lambda"]
    sqrt_lam = geom["sqrt_lambda"]
    V = geom["V"]
    H = geom["H"]
    include_longitudinal = sample_l >= 1
    include_homogeneous_connection = sample_l == 0
    include_E = sample_l >= 2
    include_X = sample_l >= 1

    n_l = sp.Symbol(f"n_l{sample_l}")
    s_l = sp.Symbol(f"s_l{sample_l}")
    B_l = sp.Symbol(f"B_l{sample_l}")
    E_l = sp.Symbol(f"E_l{sample_l}")
    X_l = sp.Symbol(f"X_l{sample_l}")

    lorentz_pairs = local["lorentz_pairs"]
    K0_l = {(a, b): sp.Symbol(f"K0_{a}{b}_l{sample_l}") for a, b in lorentz_pairs}
    k_l = {(a, b): sp.Symbol(f"k_{a}{b}_l{sample_l}") for a, b in lorentz_pairs}
    dk_l = {(a, b): sp.Symbol(f"dk_{a}{b}_l{sample_l}") for a, b in lorentz_pairs}
    omega0_l = {
        (a, b, i): sp.Symbol(f"om_{a}{b}_{i}_l0")
        for a, b in lorentz_pairs
        for i in range(1, 4)
    }
    domega0_l = {
        (a, b, i): sp.Symbol(f"dom_{a}{b}_{i}_l0")
        for a, b in lorentz_pairs
        for i in range(1, 4)
    }

    q_cfg: List[sp.Symbol] = [n_l, s_l] + [K0_l[(a, b)] for a, b in lorentz_pairs]
    sectors: Dict[sp.Symbol, str] = {n_l: "lapse", s_l: "trace_scalar"}
    for a, b in lorentz_pairs:
        sectors[K0_l[(a, b)]] = "temporal_connection_K0"

    if include_longitudinal:
        q_cfg.insert(1, B_l)
        sectors[B_l] = "scalar_shift_longitudinal"
        q_cfg.append(X_l)
        sectors[X_l] = "antisymmetric_frame_scalar_longitudinal"
        for a, b in lorentz_pairs:
            q_cfg.append(k_l[(a, b)])
            sectors[k_l[(a, b)]] = "spatial_connection_longitudinal"
    if include_homogeneous_connection:
        for a, b in lorentz_pairs:
            for i in range(1, 4):
                q_cfg.append(omega0_l[(a, b, i)])
                sectors[omega0_l[(a, b, i)]] = "homogeneous_spatial_connection"
    if include_E:
        q_cfg.insert(3 if include_longitudinal else 2, E_l)
        sectors[E_l] = "traceless_scalar_E"

    velocity_by_q = {k_l[(a, b)]: dk_l[(a, b)] for a, b in lorentz_pairs} if include_longitudinal else {}
    if include_homogeneous_connection:
        velocity_by_q.update({omega0_l[(a, b, i)]: domega0_l[(a, b, i)] for a, b in lorentz_pairs for i in range(1, 4)})

    subs: Dict[sp.Symbol, sp.Expr] = {
        local["n_pert"]: n_l,
        local["s"]: s_l,
    }
    for j in range(3):
        subs[local["NADM"][j]] = sqrt_lam * B_l * V[j] if include_longitudinal else sp.S.Zero

    qA00, qA11, qA01, qA02, qA12 = local["qA_vars"]
    qA_subs = {
        qA00: E_l * H[(0, 0)] if include_E else sp.S.Zero,
        qA11: E_l * H[(1, 1)] if include_E else sp.S.Zero,
        qA01: E_l * H[(0, 1)] if include_E else sp.S.Zero,
        qA02: E_l * H[(0, 2)] if include_E else sp.S.Zero,
        qA12: E_l * H[(1, 2)] if include_E else sp.S.Zero,
    }
    subs.update(qA_subs)

    X01, X02, X12 = local["X_vars"]
    X_components = [X01, X02, X12]
    for idx, x_var in enumerate(X_components):
        subs[x_var] = sqrt_lam * X_l * V[idx] if include_X else sp.S.Zero

    for a, b in lorentz_pairs:
        subs[local["K0_syms"][(a, b)]] = K0_l[(a, b)]
        for i in range(1, 4):
            if include_homogeneous_connection:
                subs[local["omega_syms"][(a, b, i)]] = omega0_l[(a, b, i)]
                subs[local["dot_omega_pos_syms"][(a, b, i)]] = domega0_l[(a, b, i)]
            elif include_longitudinal:
                subs[local["omega_syms"][(a, b, i)]] = sqrt_lam * k_l[(a, b)] * V[i - 1]
                subs[local["dot_omega_pos_syms"][(a, b, i)]] = sqrt_lam * dk_l[(a, b)] * V[i - 1]
            else:
                subs[local["omega_syms"][(a, b, i)]] = sp.S.Zero
                subs[local["dot_omega_pos_syms"][(a, b, i)]] = sp.S.Zero

    L_proj = sp.expand(local["L_total_A"].subs(subs))
    L_static = sp.expand(local["L_static_A"].subs(subs))
    L_kin = sp.expand(local["L_kin_A"].subs(subs))

    absent = []
    if not include_longitudinal:
        absent += ["B_l", "X_l", "k_ab_l"]
    if include_homogeneous_connection:
        absent.append("longitudinal_connection_ansatz_replaced_by_homogeneous_omega_ab_i")
    if not include_E:
        absent.append("E_l")

    return {
        "l_sector": l_sector,
        "sample_l": sample_l,
        "lambda": lam,
        "q_cfg": q_cfg,
        "velocity_by_q": velocity_by_q,
        "sectors": sectors,
        "L_proj": L_proj,
        "L_static": L_static,
        "L_kin": L_kin,
        "subs": subs,
        "absent": absent,
        "projection": geom,
        "retained_homogeneous_connection": include_homogeneous_connection,
    }


def build_dirac_for_projection(proj: Dict[str, Any], local: Dict[str, Any]) -> Dict[str, Any]:
    q_cfg: List[sp.Symbol] = proj["q_cfg"]
    p_cfg = [sp.Symbol(f"p_{q}") for q in q_cfg]
    p_by_q = dict(zip(q_cfg, p_cfg))
    velocity_by_q: Dict[sp.Symbol, sp.Symbol] = proj["velocity_by_q"]
    vars_all = q_cfg + p_cfg
    subs_rank = {local["N"]: sp.S.One, local["q"]: sp.S.One}

    def linear_bracket_with_H(expr: sp.Expr, H_expr: sp.Expr) -> sp.Expr:
        # H_c has no momentum dependence. For a linear constraint f,
        # {f,H_c} = - sum_Q (df/dp_Q) dH_c/dQ.
        return sp.expand(-sum(sp.diff(expr, p_by_q[Q]) * sp.diff(H_expr, Q) for Q in q_cfg))

    primary_rows: List[Dict[str, Any]] = []
    primaries: List[sp.Expr] = []
    for Q in q_cfg:
        vel = velocity_by_q.get(Q)
        momentum_expr = sp.diff(proj["L_proj"], vel) if vel is not None else sp.S.Zero
        phi = sp.expand(p_by_q[Q] - momentum_expr)
        primaries.append(phi)
        primary_rows.append(
            {
                "constraint_id": f"PHI_PRIMARY_{Q}",
                "constraint_kind": "primary",
                "expression": sp.sstr(phi),
                "sector": proj["sectors"].get(Q, "unknown"),
                "independent": "YES",
                "notes": "p - dL/dvelocity" if vel is not None else "p=0 because velocity is absent",
            }
        )

    H_c = -proj["L_static"]
    secondaries_raw = [linear_bracket_with_H(phi, H_c) for phi in primaries]

    independent = list(primaries)
    secondary_rows: List[Dict[str, Any]] = []
    secondaries_ind: List[sp.Expr] = []
    for idx, sec in enumerate(secondaries_raw):
        sec_num = sp.expand(sec.subs(subs_rank))
        independent_flag = is_independent(sec_num, [sp.expand(c.subs(subs_rank)) for c in independent], vars_all)
        if independent_flag:
            independent.append(sec)
            secondaries_ind.append(sec)
        Q = q_cfg[idx]
        secondary_rows.append(
            {
                "constraint_id": f"ZETA_{Q}",
                "constraint_kind": "secondary",
                "expression": sp.sstr(sp.expand(sec)),
                "sector": proj["sectors"].get(Q, "unknown"),
                "independent": "YES" if independent_flag else "NO",
                "notes": "primary preservation {phi,Hc}; independence tested after q=N=1 rank substitution",
            }
        )

    constraints = primaries + secondaries_ind
    n_constraints = len(constraints)
    constraints_num = [sp.expand(c.subs(subs_rank)) for c in constraints]
    coeff_mat = sp.Matrix([coeff_vector(c, vars_all) for c in constraints_num])
    J = sp.zeros(2 * len(q_cfg), 2 * len(q_cfg))
    n_q = len(q_cfg)
    for idx in range(n_q):
        J[idx, n_q + idx] = 1
        J[n_q + idx, idx] = -1
    C_num = coeff_mat * J * coeff_mat.T
    rank_C = int(C_num.rank())
    kernel_dim = n_constraints - rank_C
    N_config = len(q_cfg)
    N_phase = 2 * N_config
    N_first = n_constraints - rank_C
    N_second = rank_C
    N_phys = sp.Rational(N_phase - 2 * N_first - N_second, 2)

    tertiary_count = 0
    tertiary_details: List[str] = []
    H_c_num = sp.expand(H_c.subs(subs_rank))
    if kernel_dim > 0:
        nullspace = C_num.nullspace()
        for null_idx, vec in enumerate(nullspace):
            combo = sp.expand(sum(coeff * constraint for coeff, constraint in zip(vec, constraints_num)))
            consistency_num = sp.expand(linear_bracket_with_H(combo, H_c_num))
            adds = is_independent(consistency_num, constraints_num, vars_all)
            if adds:
                tertiary_count += 1
                tertiary_details.append(f"null_{null_idx}: {sp.sstr(consistency_num)}")

    return {
        "q_cfg": q_cfg,
        "p_cfg": p_cfg,
        "primary_rows": primary_rows,
        "secondary_rows": secondary_rows,
        "constraints": constraints,
        "C_num": C_num,
        "rank_C": rank_C,
        "kernel_dim": kernel_dim,
        "N_config": N_config,
        "N_phase": N_phase,
        "N_primary": len(primaries),
        "N_secondary": len(secondaries_ind),
        "N_constraints": n_constraints,
        "N_first": N_first,
        "N_second": N_second,
        "N_phys": N_phys,
        "tertiary_count": tertiary_count,
        "tertiary_details": tertiary_details,
    }


def build_dirac_for_second_order_lagrangian(L_red: sp.Expr, q_red: List[sp.Symbol], v_red: List[sp.Symbol]) -> Dict[str, Any]:
    p_red = [sp.Symbol(f"P_{Q}") for Q in q_red]
    vars_all = q_red + p_red
    momenta = [sp.diff(L_red, vel) for vel in v_red]
    W = sp.Matrix([[sp.diff(mom, vel) for vel in v_red] for mom in momenta])
    W_rank = int(W.rank())
    _, pivots = W.rref()
    pivot_indices = list(pivots)
    pivot_velocities = [v_red[idx] for idx in pivot_indices]
    equations = [sp.Eq(p_red[idx] - momenta[idx], 0) for idx in pivot_indices]
    solutions = sp.solve(equations, pivot_velocities, dict=True, simplify=False)
    if not solutions and pivot_velocities:
        return {"status": "LEGENDRE_SOLVE_FAILED", "W_rank": W_rank}
    sol = solutions[0] if solutions else {}
    nonpivot_indices = [idx for idx in range(len(q_red)) if idx not in pivot_indices]
    nonpivot_zero = {v_red[idx]: sp.S.Zero for idx in nonpivot_indices}

    primaries = [
        sp.expand((p_red[idx] - momenta[idx]).subs(sol).subs(nonpivot_zero))
        for idx in nonpivot_indices
    ]
    H_c = sp.expand(sum(P * vel for P, vel in zip(p_red, v_red)) - L_red).subs(sol).subs(nonpivot_zero)

    independent = list(primaries)
    secondaries_ind: List[sp.Expr] = []
    for phi in primaries:
        sec = sp.expand(poisson(phi, H_c, q_red, p_red))
        if is_independent(sec, independent, vars_all):
            independent.append(sec)
            secondaries_ind.append(sec)

    C = sp.Matrix([[poisson(c1, c2, q_red, p_red) for c2 in independent] for c1 in independent])
    rank_C = int(C.rank()) if independent else 0
    N_constraints = len(independent)
    N_first = N_constraints - rank_C
    N_second = rank_C
    N_phys = sp.Rational(2 * len(q_red) - 2 * N_first - N_second, 2)

    return {
        "status": "COMPUTED",
        "W_rank": W_rank,
        "N_config": len(q_red),
        "N_primary": len(primaries),
        "N_secondary": len(secondaries_ind),
        "N_constraints": N_constraints,
        "rank_C": rank_C,
        "N_first": N_first,
        "N_second": N_second,
        "N_phys": N_phys,
        "H_c": H_c,
    }
