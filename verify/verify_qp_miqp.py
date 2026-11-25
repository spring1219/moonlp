"""
Independent mathematical cross-verification for moonlp v1.0.0 QP & MIQP models.
Uses pure NumPy / SciPy (SLSQP + exact combinatorial subset enumeration for MIQP)
to independently verify every objective value reported by moonlp.
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


if __name__ == "__main__":
    verify_markowitz_qp()
    verify_cardinality_miqp()
    verify_linear_mpc_qp()
    verify_svm_dual_qp()
    print("\nALL 4 QP/MIQP MODELS VERIFIED TO <= 1e-6 AGAINST SCIPY + EXACT COMBINATORIAL ENUMERATION!")
