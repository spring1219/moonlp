"""
Independent mathematical cross-verification for moonlp v1.1.0 SOCP & MISOCP models.
Uses SciPy SLSQP (with exact analytical Jacobians for portable convergence across SciPy versions)
+ exact combinatorial enumeration over all binary subsets to cross-check the recorded objective
values produced by moonlp's Nesterov-Todd Conic IPM and MISOCP Branch-and-Bound solver.
"""
import itertools
import math
import numpy as np
from scipy.optimize import minimize

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


def check(label, got, expected, tol=1e-5):
    err = abs(got - expected)
    status = "PASS" if err <= tol else "FAIL"
    print(f"[{status}] {label}: ref={got:.10f}, moonlp={expected:.10f}, |diff|={err:.2e}")
    assert err <= tol, f"Mismatch on {label}: {got} vs {expected}"


def verify_smallest_enclosing_disk():
    # Analytical circumcircle of (0,0), (4,0), (2,3.5)
    exact_r = 16.25 / 7.0
    check("Smallest Enclosing Disk SOCP (6 Lorentz cones)", exact_r, 2.3214285786404494, tol=1e-6)


def verify_robust_portfolio_socp():
    kappa = 1.5

    def obj(w):
        vol = math.sqrt(max(float(w @ SIGMA @ w), 1e-18))
        return -(float(MU @ w) - kappa * vol)

    def obj_jac(w):
        vol = math.sqrt(max(float(w @ SIGMA @ w), 1e-18))
        return -MU + kappa * (SIGMA @ w) / vol

    res = minimize(
        obj,
        np.ones(6) / 6.0,
        jac=obj_jac,
        bounds=[(0.0, 0.5)] * 6,
        constraints=[{"type": "eq", "fun": lambda w: np.sum(w) - 1.0, "jac": lambda w: np.ones(6)}],
        method="SLSQP",
        options={"ftol": 1e-10, "maxiter": 500},
    )
    assert res.success, res.message
    check("Robust Ellipsoidal Portfolio SOCP (kappa=1.5)", -res.fun, -0.02622927526278558, tol=1e-6)


def verify_vol_capped_cardinality_misocp():
    vol_cap = 0.105
    best_ret = -1e9
    best_sub = None
    for k in range(1, 4):
        for comb in itertools.combinations(range(6), k):
            idx = list(comb)
            sub_mu = MU[idx]
            sub_sig = SIGMA[np.ix_(idx, idx)]
            res = minimize(
                lambda w: -float(sub_mu @ w),
                np.ones(k) / k,
                jac=lambda w: -sub_mu,
                bounds=[(0.1, 0.5)] * k,
                constraints=[
                    {"type": "eq", "fun": lambda w: np.sum(w) - 1.0, "jac": lambda w: np.ones(k)},
                    {
                        "type": "ineq",
                        "fun": lambda w: vol_cap * vol_cap - float(w @ sub_sig @ w),
                        "jac": lambda w: -2.0 * (sub_sig @ w),
                    },
                ],
                method="SLSQP",
                options={"ftol": 1e-10, "maxiter": 500},
            )
            if (
                res.success
                and abs(np.sum(res.x) - 1.0) < 1e-6
                and math.sqrt(max(float(res.x @ sub_sig @ res.x), 0.0)) <= vol_cap + 1e-6
            ):
                ret = float(sub_mu @ res.x)
                if ret > best_ret:
                    best_ret = ret
                    best_sub = ([NAMES[i] for i in idx], [round(float(v), 6) for v in res.x])
    print(f"    Exact Vol-Capped MISOCP subset = {list(zip(best_sub[0], best_sub[1]))}")
    check("Volatility-Capped Cardinality MISOCP (K<=3, vol<=10.5%)", best_ret, 0.0811046240541594, tol=1e-6)


def verify_chance_constrained_socp():
    z_a = 1.644854
    mean = np.array([2.0, 3.0, 1.5])
    std = np.array([0.5, 0.8, 0.4])
    std2 = std * std
    profit = np.array([50.0, 65.0, 45.0])

    def cc_fun(x):
        s = math.sqrt(max(float(np.sum(std2 * x * x)), 1e-18))
        return 100.0 - float(mean @ x) - z_a * s

    def cc_jac(x):
        s = math.sqrt(max(float(np.sum(std2 * x * x)), 1e-18))
        return -mean - z_a * (std2 * x) / s

    res = minimize(
        lambda x: -float(profit @ x),
        np.array([10.0, 10.0, 10.0]),
        jac=lambda x: -profit,
        bounds=[(0.0, None)] * 3,
        constraints=[
            {"type": "ineq", "fun": lambda x: 40.0 - float(np.sum(x)), "jac": lambda x: -np.ones(3)},
            {"type": "ineq", "fun": cc_fun, "jac": cc_jac},
        ],
        method="SLSQP",
        options={"ftol": 1e-10, "maxiter": 500},
    )
    assert res.success, res.message
    sigma_x = math.sqrt(float(np.sum(std2 * res.x * res.x)))
    z_score = (100.0 - float(mean @ res.x)) / sigma_x
    exceed_prob = 0.5 * math.erfc(z_score / math.sqrt(2.0))
    print(f"    95% Chance-Constrained x = {np.round(res.x, 4)}, exceedance prob = {exceed_prob * 100:.4f}% (vs 50.00% for LP)")
    check("95% Gaussian Chance-Constrained SOCP", -res.fun, 2046.9776935197942, tol=1e-5)


def eval_rocket_cost(u_seq, active_z=None):
    N = 4
    dt, h2 = 0.5, 0.125
    px, py, vx, vy = 2.0, 1.5, 0.0, 0.0
    cost = 0.5 * sum(active_z) if active_z is not None else 0.0
    for t in range(N):
        ux, uy = u_seq[2 * t], u_seq[2 * t + 1]
        px = px + dt * vx + h2 * ux
        py = py + dt * vy + h2 * uy
        vx = vx + dt * ux
        vy = vy + dt * uy
        qp = 20.0 if t == N - 1 else 10.0
        qv = 4.0 if t == N - 1 else 2.0
        cost += 0.5 * (qp * (px * px + py * py) + qv * (vx * vx + vy * vy))
        cost += 0.25 * (ux * ux + uy * uy)
    return cost


def verify_rocket_socp_and_misocp():
    # Continuous 4-step circular-thrust SOCP
    cons = [{"type": "ineq", "fun": lambda u, t=t: 1.0 - (u[2 * t] ** 2 + u[2 * t + 1] ** 2)} for t in range(4)]
    res = minimize(
        eval_rocket_cost,
        np.full(8, -0.4),
        bounds=[(-1.0, 1.0)] * 8,
        constraints=cons,
        method="SLSQP",
        options={"ftol": 1e-10, "maxiter": 500},
    )
    assert res.success, res.message
    check("2D Isotropic Circular-Thrust Rocket SOCP (N=4)", res.fun, 69.74245692941984, tol=1e-5)

    # Exact brute-force over all <=2-burn binary schedules (11 schedules)
    best = float("inf")
    best_z = None
    for mask in range(16):
        z = [(mask >> t) & 1 for t in range(4)]
        if sum(z) > 2:
            continue
        bnds = [(-1.0, 1.0) if z[t] else (0.0, 0.0) for t in range(4) for _ in range(2)]
        x0 = np.array([-0.4 if z[t] else 0.0 for t in range(4) for _ in range(2)])
        res_s = minimize(
            lambda u, z=z: eval_rocket_cost(u, z),
            x0,
            bounds=bnds,
            constraints=cons,
            method="SLSQP",
            options={"ftol": 1e-10, "maxiter": 500},
        )
        if res_s.success and res_s.fun < best:
            best = res_s.fun
            best_z = z
    print(f"    Exact Switched 2D Rocket burn schedule = z={best_z}")
    check("2D Switched Circular-Thrust Rocket MISOCP (<=2 burns)", best, 75.20312590949864, tol=1e-5)


if __name__ == "__main__":
    verify_smallest_enclosing_disk()
    verify_robust_portfolio_socp()
    verify_vol_capped_cardinality_misocp()
    verify_chance_constrained_socp()
    verify_rocket_socp_and_misocp()
    print("\nALL 6 SOCP / MISOCP MODELS VERIFIED AGAINST SCIPY + EXACT COMBINATORIAL ENUMERATION!")
