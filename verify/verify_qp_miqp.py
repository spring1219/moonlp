"""
Independent mathematical cross-verification for moonlp v1.1.0 QP & MIQP models.
Uses pure NumPy / SciPy (SLSQP + exact combinatorial subset enumeration for MIQP)
to cross-check the recorded objective values produced by moonlp.
"""
import itertools
import math
import sys

try:
    import numpy as np
    from scipy.optimize import minimize
except ImportError:
    print("Install numpy & scipy to run QP/MIQP verification: pip install numpy scipy")
    sys.exit(1)

MU = np.array([0.11, 0.10, 0.06, 0.035, 0.085, 0.055])
SIGMA = np.array([
    [0.0400, 0.0240, 0.0040, -0.0020, 0.0180, 0.0020],
    [0.0240, 0.0484, 0.0035, -0.0015, 0.0160, 0.0050],
    [0.0040, 0.0035, 0.0064, 0.0025, 0.0045, 0.0015],
    [-0.0020, -0.0015, 0.0025, 0.0025, -0.0010, 0.0010],
    [0.0180, 0.0160, 0.0045, -0.0010, 0.0324, 0.0030],
    [0.0020, 0.0050, 0.0015, 0.0010, 0.0030, 0.0225],
])
NAMES = ["US_Eq", "Intl_Eq", "Corp_Bond", "Gov_Bond", "REIT", "Gold"]


def check(label, got, expected, tol=1e-6):
    err = abs(got - expected)
    status = "PASS" if err <= tol else "FAIL"
    print(f"[{status}] {label}: scipy={got:.12f}, moonlp={expected:.12f}, |diff|={err:.2e}")
    assert err <= tol, f"Mismatch on {label}: {got} vs {expected}"


def verify_markowitz_qp():
    res = minimize(
        lambda w: 0.5 * float(w @ SIGMA @ w),
        np.ones(6) / 6.0,
        jac=lambda w: SIGMA @ w,
        bounds=[(0.0, 0.5)] * 6,
        constraints=[
            {"type": "eq", "fun": lambda w: np.sum(w) - 1.0, "jac": lambda w: np.ones(6)},
            {"type": "ineq", "fun": lambda w: float(MU @ w) - 0.08, "jac": lambda w: MU},
        ],
        method="SLSQP",
        options={"ftol": 1e-14, "maxiter": 500},
    )
    assert res.success
    check("Markowitz 6-Asset QP (r>=8%)", res.fun, 0.004958677971393078, tol=1e-7)


def verify_cardinality_miqp():
    # Exact brute-force enumeration over all subsets of size 1..3
    best_obj = float("inf")
    best_subset = None
    best_w = None
    for k in range(1, 4):
        for comb in itertools.combinations(range(6), k):
            idx = list(comb)
            sub_mu = MU[idx]
            sub_sig = SIGMA[np.ix_(idx, idx)]
            res = minimize(
                lambda w: 0.5 * float(w @ sub_sig @ w),
                np.ones(k) / k,
                jac=lambda w: sub_sig @ w,
                bounds=[(0.1, 0.5)] * k,
                constraints=[
                    {"type": "eq", "fun": lambda w: np.sum(w) - 1.0, "jac": lambda w: np.ones(k)},
                    {"type": "ineq", "fun": lambda w: float(sub_mu @ w) - 0.08, "jac": lambda w: sub_mu},
                ],
                method="SLSQP",
                options={"ftol": 1e-14, "maxiter": 500},
            )
            if res.success and abs(np.sum(res.x) - 1.0) < 1e-7 and float(sub_mu @ res.x) >= 0.08 - 1e-7:
                if res.fun < best_obj:
                    best_obj = res.fun
                    best_subset = [NAMES[i] for i in idx]
                    best_w = res.x
    print(f"    Exact MIQP subset = {list(zip(best_subset, [round(float(v), 6) for v in best_w]))}")
    check("Cardinality-Constrained Markowitz MIQP (K<=3)", best_obj, 0.0052019200000437155, tol=1e-7)


def verify_linear_mpc_qp():
    N = 8
    dt = 0.5
    half_dt2 = 0.5 * dt * dt
    # Vars: p_1..p_N (0..N-1), v_1..v_N (N..2N-1), u_0..u_{N-1} (2N..3N-1)
    q_diag = np.zeros(3 * N)
    for t in range(N):
        q_diag[t] = 20.0 if t == N - 1 else 10.0
        q_diag[N + t] = 4.0 if t == N - 1 else 2.0
        q_diag[2 * N + t] = 0.5

    bounds = [(-1e9, 1e9)] * N + [(-1.5, 1.5)] * N + [(-1.0, 1.0)] * N
    A_eq = []
    b_eq = []
    for t in range(N):
        rp = np.zeros(3 * N)
        rv = np.zeros(3 * N)
        if t == 0:
            rp[0] = 1.0
            rp[2 * N] = -half_dt2
            rv[N] = 1.0
            rv[2 * N] = -dt
            A_eq.extend([rp, rv])
            b_eq.extend([2.0, 0.0])
        else:
            rp[t] = 1.0
            rp[t - 1] = -1.0
            rp[N + t - 1] = -dt
            rp[2 * N + t] = -half_dt2
            rv[N + t] = 1.0
            rv[N + t - 1] = -1.0
            rv[2 * N + t] = -dt
            A_eq.extend([rp, rv])
            b_eq.extend([0.0, 0.0])
    A_eq = np.array(A_eq)
    b_eq = np.array(b_eq)

    res = minimize(
        lambda x: 0.5 * float(np.sum(q_diag * x * x)),
        np.zeros(3 * N),
        jac=lambda x: q_diag * x,
        bounds=bounds,
        constraints=[{"type": "eq", "fun": lambda x: A_eq @ x - b_eq, "jac": lambda x: A_eq}],
        method="SLSQP",
        options={"ftol": 1e-14, "maxiter": 1000},
    )
    assert res.success
    check("Constrained Linear MPC QP (N=8)", res.fun, 38.79579681548938, tol=1e-6)


def verify_svm_dual_qp():
    xs = np.array([
        [2.0, 2.5], [2.5, 1.8], [3.0, 3.0], [1.6, 2.2],
        [-2.0, -2.2], [-1.5, -2.6], [-2.8, -1.8], [-1.2, -1.5],
    ])
    ys = np.array([1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0])
    K = np.outer(ys, ys) * (xs @ xs.T)
    res = minimize(
        lambda a: 0.5 * float(a @ K @ a) - float(np.sum(a)),
        np.zeros(8),
        jac=lambda a: K @ a - np.ones(8),
        bounds=[(0.0, 10.0)] * 8,
        constraints=[{"type": "eq", "fun": lambda a: float(ys @ a), "jac": lambda a: ys}],
        method="SLSQP",
        options={"ftol": 1e-14, "maxiter": 500},
    )
    assert res.success
    check("Soft-Margin SVM Dual QP (rank-2 PSD)", res.fun, -0.09289363677146076, tol=1e-7)


def verify_unit_commitment_miqp():
    p_min = np.array([30.0, 20.0, 15.0, 10.0, 5.0])
    p_max = np.array([80.0, 60.0, 50.0, 40.0, 30.0])
    a = np.array([0.008, 0.012, 0.018, 0.025, 0.032])
    b = np.array([18.0, 22.0, 26.0, 31.0, 35.0])
    c = np.array([220.0, 140.0, 90.0, 55.0, 35.0])
    best = float("inf")
    best_u = None
    for mask in range(1, 32):
        idx = [g for g in range(5) if (mask >> g) & 1]
        if np.sum(p_max[idx]) < 160.0 or np.sum(p_min[idx]) > 140.0:
            continue
        k = len(idx)
        ai, bi, ci = a[idx], b[idx], float(np.sum(c[idx]))
        res = minimize(
            lambda p: float(np.sum(ai * p * p + bi * p)) + ci,
            p_min[idx],
            jac=lambda p: 2.0 * ai * p + bi,
            bounds=[(float(p_min[g]), float(p_max[g])) for g in idx],
            constraints=[{"type": "eq", "fun": lambda p: np.sum(p) - 140.0, "jac": lambda p: np.ones(k)}],
            method="SLSQP",
            options={"ftol": 1e-14, "maxiter": 500},
        )
        if res.success and abs(np.sum(res.x) - 140.0) < 1e-6 and res.fun < best:
            best = res.fun
            best_u = (idx, res.x)
    print(f"    Exact UC-MIQP units = {[(i, round(float(v), 4)) for i, v in zip(best_u[0], best_u[1])]}")
    check("Thermal Unit Commitment MIQP (140MW+20MW)", best, 3308.3000000000166, tol=1e-5)


def verify_best_subset_ridge_miqp():
    X = np.array([
        [1.0, 0.8, -0.3, 0.5, 0.4],
        [2.0, 1.5, 0.4, -1.0, -0.2],
        [-1.5, -1.1, 0.9, 1.2, 0.7],
        [0.5, 0.6, -1.2, 2.0, -0.5],
        [-2.0, -1.7, 0.2, -1.5, 0.9],
        [1.8, 1.2, 1.1, 0.8, -0.8],
        [-0.7, -0.4, -0.8, 1.7, 0.3],
        [1.2, 1.1, 0.5, -2.2, -0.6],
    ])
    y = np.array([1.62, 6.78, -5.88, -2.37, -2.31, 3.09, -4.83, 6.95])
    lam = 0.05
    best = float("inf")
    best_comb = None
    for k in range(1, 3):
        for comb in itertools.combinations(range(5), k):
            idx = list(comb)
            Xs = X[:, idx]
            beta = np.linalg.solve(Xs.T @ Xs + lam * np.eye(k), Xs.T @ y)
            loss = 0.5 * float(np.sum((Xs @ beta - y) ** 2)) + 0.5 * lam * float(np.sum(beta ** 2))
            if loss < best:
                best = loss
                best_comb = (idx, beta)
    print(f"    Exact Best-Subset Ridge support = {best_comb[0]}, beta = {[round(float(v), 6) for v in best_comb[1]]}")
    check("L0 Best-Subset Ridge Regression MIQP (K<=2)", best, 0.23809506861601154, tol=1e-7)


def verify_hybrid_mpc_miqp():
    N = 5
    dt = 0.5
    h2 = 0.125
    best = float("inf")
    best_seq = None
    for mask in range(32):
        z = [(mask >> t) & 1 for t in range(N)]
        if sum(z) > 3:
            continue
        def eval_cost(u_act, active_idx, z_vec):
            u = np.zeros(N)
            for k, t in enumerate(active_idx):
                u[t] = u_act[k]
            p, v = 2.0, 0.0
            cost = 0.2 * sum(z_vec) + 0.25 * float(np.sum(u * u))
            for t in range(N):
                p = p + dt * v - h2 * u[t]
                v = v - dt * u[t]
                if abs(v) > 2.0 + 1e-7:
                    return 1e12
                qp = 20.0 if t == N - 1 else 10.0
                qv = 4.0 if t == N - 1 else 2.0
                cost += 0.5 * (qp * p * p + qv * v * v)
            return cost
        active = [t for t in range(N) if z[t] == 1]
        if not active:
            val = eval_cost(np.array([]), [], z)
            if val < best:
                best = val
                best_seq = (z, [])
        else:
            res = minimize(
                lambda ua: eval_cost(ua, active, z),
                np.full(len(active), 0.6),
                bounds=[(0.3, 1.0)] * len(active),
                method="SLSQP",
                options={"ftol": 1e-14, "maxiter": 500},
            )
            if res.success and res.fun < best:
                best = res.fun
                best_seq = (z, [round(float(v), 4) for v in res.x])
    print(f"    Exact Hybrid MPC firing = z={best_seq[0]}, u={best_seq[1]}")
    check("Hybrid Switched-Thruster Deadzone MPC MIQP", best, 41.2281250000146, tol=1e-5)


if __name__ == "__main__":
    verify_markowitz_qp()
    verify_cardinality_miqp()
    verify_linear_mpc_qp()
    verify_svm_dual_qp()
    verify_unit_commitment_miqp()
    verify_best_subset_ridge_miqp()
    verify_hybrid_mpc_miqp()
    print("\nALL 7 QP/MIQP MODELS VERIFIED AGAINST SCIPY + EXACT COMBINATORIAL ENUMERATION!")
