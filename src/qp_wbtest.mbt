///|
/// Comprehensive Convex QP & Mixed-Integer QP (MIQP) Verification Suite (v1.0.0):
/// Validates analytical solutions, KKT + Wolfe Strong Duality certification,
/// non-convex / infeasible / unbounded guards, CPLEX .lp & .qps QUADOBJ round-trips,
/// and 5 concrete real-world optimization models strictly beyond linear LP/MILP solvers.

///|
test "64. [v1.0 QP] Analytical 2D Convex QP with active inequality & Wolfe Strong Duality KKT" {
  // min (x - 2)^2 + (y - 3)^2 = x^2 + y^2 - 4x - 6y + 13
  // s.t. x + y <= 3, x >= 0, y >= 0
  // Analytical optimum: projection of (2, 3) onto x + y = 3 => x = 1, y = 2, obj = 2.0
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("x", obj=-4.0)
  let y = m.add_continuous_var("y", obj=-6.0)
  m.add_quad_obj_term(x, x, 1.0)
  m.add_quad_obj_term(y, y, 1.0)
  m.set_obj_constant(13.0)
  m.add_constraint(
    "c1",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    3.0,
  )

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_approx(sol.primal[0], 1.0, 1.0e-7)
  assert_approx(sol.primal[1], 2.0, 1.0e-7)
  assert_approx(sol.objective, 2.0, 1.0e-7)
  assert_approx(sol.dual[0], -2.0, 1.0e-7)
  assert_true(m.verify_solution(sol, tol=1.0e-7))

  let cert = m.certify_qp_solution(sol, tol=1.0e-7)
  assert_true(cert.certified_optimal)
  assert_true(cert.duality_gap < 1.0e-7)
}

///|
test "65. [v1.0 QP] Concave Maximize QP with off-diagonal coupling & box bounds" {
  // max -(x^2 + xy + y^2) + 4x + 5y  s.t. 0 <= x <= 1.5, 0 <= y <= 2.5, x + y <= 3
  // Unconstrained stationary point: 2x + y = 4, x + 2y = 5 => x = 1, y = 2 (strictly interior except x+y=3!)
  // Objective at (1, 2): -(1 + 2 + 4) + 4 + 10 = 7.0
  let m = Model::new(Maximize)
  let x = m.add_continuous_var("x", lb=0.0, ub=1.5, obj=4.0)
  let y = m.add_continuous_var("y", lb=0.0, ub=2.5, obj=5.0)
  m.add_quad_obj_term(x, x, -1.0)
  m.add_quad_obj_term(x, y, -1.0)
  m.add_quad_obj_term(y, y, -1.0)
  m.add_constraint(
    "cap",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    3.0,
  )

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_approx(sol.primal[0], 1.0, 1.0e-7)
  assert_approx(sol.primal[1], 2.0, 1.0e-7)
  assert_approx(sol.objective, 7.0, 1.0e-7)
  let cert = m.certify_qp_solution(sol, tol=1.0e-7)
  assert_true(cert.certified_optimal)
}

///|
test "66. [v1.0 QP] Unconstrained & pure equality-constrained QP (m == 0 fast path)" {
  // min 0.5*(4 x^2 + 2 x y + 2 y^2) - 2 x - 3 y  s.t. x + y = 2 (free variables!)
  // KKT: 4x + y - 2 + lam = 0, x + 2y - 3 + lam = 0, x + y = 2
  // Subtracting: 3x - y + 1 = 0 => y = 3x + 1 => 4x + 1 = 2 => x = 0.25, y = 1.75
  // obj = 2*(0.0625) + (0.25*1.75) + (1.75^2) - 0.5 - 5.25 = 0.125 + 0.4375 + 3.0625 - 5.75 = -2.125
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("x", lb=-1.0e300, ub=1.0e300, obj=-2.0)
  let y = m.add_continuous_var("y", lb=-1.0e300, ub=1.0e300, obj=-3.0)
  m.set_quad_obj_matrix([x, y], [[4.0, 1.0], [1.0, 2.0]], half=true)
  m.add_constraint(
    "eq1",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Eq,
    2.0,
  )

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_approx(sol.primal[0], 0.25, 1.0e-8)
  assert_approx(sol.primal[1], 1.75, 1.0e-8)
  assert_approx(sol.objective, -2.125, 1.0e-8)
  let cert = m.certify_qp_solution(sol, tol=1.0e-7)
  assert_true(cert.certified_optimal)
}

///|
test "67. [v1.0 QP] NonConvex, Infeasible, and Unbounded PSD QP status guards" {
  // 1) Non-convex minimization: min -x^2 + y^2
  let m_nc = Model::new(Minimize)
  let x1 = m_nc.add_continuous_var("x", lb=0.0, ub=5.0)
  let y1 = m_nc.add_continuous_var("y", lb=0.0, ub=5.0)
  m_nc.add_quad_obj_term(x1, x1, -1.0)
  m_nc.add_quad_obj_term(y1, y1, 1.0)
  assert_eq(m_nc.solve().status, NonConvex)

  // 2) Infeasible convex QP
  let m_inf = Model::new(Minimize)
  let x2 = m_inf.add_continuous_var("x", lb=0.0)
  m_inf.add_quad_obj_term(x2, x2, 1.0)
  m_inf.add_constraint("c1", [{ var_id: x2, coeff: 1.0, }], Le, 1.0)
  m_inf.add_constraint("c2", [{ var_id: x2, coeff: 1.0, }], Ge, 3.0)
  assert_eq(m_inf.solve().status, Infeasible)

  // 3) Unbounded PSD QP: min x^2 - y  s.t. x >= 0, y >= 0 (H = diag(2, 0), slope -1 along y)
  let m_unb = Model::new(Minimize)
  let x3 = m_unb.add_continuous_var("x", lb=0.0)
  let _ = m_unb.add_continuous_var("y", lb=0.0, obj=-1.0)
  m_unb.add_quad_obj_term(x3, x3, 1.0)
  assert_eq(m_unb.solve().status, Unbounded)
}

///|
test "68. [v1.0 QP Model] Markowitz 6-Asset Mean-Variance Portfolio & Efficient Frontier" {
  // Trace three points on the Efficient Frontier: target_return = 0.06, 0.08, 0.10
  let m_low = build_markowitz_portfolio_qp(target_return=0.06)
  let m_mid = build_markowitz_portfolio_qp(target_return=0.08)
  let m_high = build_markowitz_portfolio_qp(target_return=0.10)

  let s_low = m_low.solve()
  let s_mid = m_mid.solve()
  let s_high = m_high.solve()

  assert_eq(s_low.status, Optimal)
  assert_eq(s_mid.status, Optimal)
  assert_eq(s_high.status, Optimal)

  // Strictly certified by independent Quadratic KKT + Wolfe Strong Duality
  assert_true(m_low.certify_qp_solution(s_low, tol=1.0e-6).certified_optimal)
  assert_true(m_mid.certify_qp_solution(s_mid, tol=1.0e-6).certified_optimal)
  assert_true(m_high.certify_qp_solution(s_high, tol=1.0e-6).certified_optimal)

  // Higher target return must require strictly higher minimum variance on the efficient frontier
  assert_true(s_low.objective > 0.0)
  assert_true(s_low.objective < s_mid.objective)
  assert_true(s_mid.objective < s_high.objective)
}

///|
test "69. [v1.0 MIQP Model] Cardinality-Constrained Markowitz Portfolio (max 3 of 6 assets)" {
  let m_qp = build_markowitz_portfolio_qp(target_return=0.08)
  let sol_qp = m_qp.solve()
  assert_eq(sol_qp.status, Optimal)

  let m_miqp = build_cardinality_portfolio_miqp(
    max_assets=3,
    min_weight=0.1,
    max_weight=0.5,
    target_return=0.08,
  )
  let sol_miqp = m_miqp.solve()
  assert_eq(sol_miqp.status, Optimal)
  assert_true(m_miqp.verify_solution(sol_miqp, tol=1.0e-6))
  assert_true(sol_miqp.mip_gap <= 1.0e-5)

  // Verify cardinality constraint: at most 3 binary indicators active, each held asset >= 0.1
  let mut active_count = 0
  let mut weight_sum = 0.0
  for i = 0; i < 6; i = i + 1 {
    let wi = sol_miqp.primal[i]
    let zi = sol_miqp.primal[6 + i]
    weight_sum = weight_sum + wi
    if zi > 0.5 {
      active_count = active_count + 1
      assert_true(wi >= 0.1 - 1.0e-6)
      assert_true(wi <= 0.5 + 1.0e-6)
    } else {
      assert_true(wi.abs() <= 1.0e-6)
    }
  }
  assert_true(active_count <= 3)
  assert_approx(weight_sum, 1.0, 1.0e-6)
  // Cardinality-restricted optimum variance must be >= unrestricted continuous QP variance
  assert_true(sol_miqp.objective >= sol_qp.objective - 1.0e-7)
}

///|
test "70. [v1.0 QP Model] Constrained Linear MPC Trajectory Tracking (8-step horizon)" {
  let m = build_linear_mpc_qp(horizon=8, x0_pos=2.0, x0_vel=0.0, u_max=1.0)
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol, tol=1.0e-6))

  let cert = m.certify_qp_solution(sol, tol=1.0e-6)
  assert_true(cert.certified_optimal)
  assert_true(cert.duality_gap < 1.0e-6)

  // Initial control u_0 (index 16) must decelerate negatively and respect |u_0| <= 1.0
  let u0 = sol.primal[16]
  assert_true(u0 <= -0.99 && u0 >= -1.0 - 1.0e-6)
  // Terminal position p_8 (index 14) must be driven close to origin from p_0 = 2.0
  let p_end = sol.primal[14]
  assert_true(p_end.abs() < 0.1)
}

///|
test "71. [v1.0 QP Model] Soft-Margin SVM Dual Training (rank-2 PSD Gram matrix)" {
  let m = build_svm_dual_qp(c_reg=10.0)
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol, tol=1.0e-6))

  let cert = m.certify_qp_solution(sol, tol=1.0e-6)
  assert_true(cert.certified_optimal)

  // Reconstruct separating hyperplane w = sum alpha_i y_i x_i and verify all 8 points classified cleanly
  let xs = [
    (2.0, 2.5),
    (2.5, 1.8),
    (3.0, 3.0),
    (1.6, 2.2),
    (-2.0, -2.2),
    (-1.5, -2.6),
    (-2.8, -1.8),
    (-1.2, -1.5),
  ]
  let ys = [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0]
  let mut w1 = 0.0
  let mut w2 = 0.0
  for i = 0; i < 8; i = i + 1 {
    let a = sol.primal[i]
    let (x1, x2) = xs[i]
    w1 = w1 + a * ys[i] * x1
    w2 = w2 + a * ys[i] * x2
  }
  // Compute bias b from any free support vector (0 < alpha_i < C)
  let mut bias = 0.0
  let mut sv_cnt = 0
  for i = 0; i < 8; i = i + 1 {
    if sol.primal[i] > 1.0e-4 && sol.primal[i] < 9.999 {
      let (x1, x2) = xs[i]
      bias = bias + (ys[i] - (w1 * x1 + w2 * x2))
      sv_cnt = sv_cnt + 1
    }
  }
  assert_true(sv_cnt >= 2)
  bias = bias / sv_cnt.to_double()

  // Every training point must have functional margin y_i * (w^T x_i + b) >= 1 - 1e-5
  for i = 0; i < 8; i = i + 1 {
    let (x1, x2) = xs[i]
    let margin = ys[i] * (w1 * x1 + w2 * x2 + bias)
    assert_true(margin >= 1.0 - 1.0e-5)
  }
}

///|
test "72. [v1.0 QP Model] Simplex-Constrained Least Squares & CPLEX .lp / .qps QUADOBJ round-trip" {
  let m = build_constrained_lsq_qp()
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.certify_qp_solution(sol, tol=1.0e-6).certified_optimal)

  // CPLEX .lp quadratic format round-trip
  let lp_str = m.to_lp_string()
  let m_lp = parse_lp(lp_str)
  let sol_lp = m_lp.solve()
  assert_eq(sol_lp.status, Optimal)
  assert_approx(sol_lp.objective, sol.objective, 1.0e-7)
  assert_true(m_lp.certify_qp_solution(sol_lp, tol=1.0e-6).certified_optimal)

  // Industry-standard .qps (MPS QUADOBJ) round-trip
  let qps_str = m.to_mps_string(name="CONSTR_LSQ")
  let m_qps = parse_mps(qps_str)
  let sol_qps = m_qps.solve()
  assert_eq(sol_qps.status, Optimal)
  assert_approx(sol_qps.objective, sol.objective, 1.0e-7)
  assert_true(m_qps.certify_qp_solution(sol_qps, tol=1.0e-6).certified_optimal)
}

///|
test "73. [v1.0 MIQP] Pure-integer quadratic assignment / discrete target tracking" {
  // min (x - 2.4)^2 + (y - 3.7)^2 = x^2 + y^2 - 4.8x - 7.4y + 19.45
  // s.t. x + y <= 5, x, y in Z_+
  // Candidates near (2.4, 3.7) with x + y <= 5:
  //   (2, 3): 0.4^2 + 0.7^2 = 0.65
  //   (1, 4): 1.4^2 + 0.3^2 = 2.05
  //   (3, 2): 0.6^2 + 1.7^2 = 3.25
  // Unique global integer optimum: (2, 3) with obj = 0.65
  let m = Model::new(Minimize)
  let x = m.add_int_var("x", lb=0.0, ub=5.0, obj=-4.8)
  let y = m.add_int_var("y", lb=0.0, ub=5.0, obj=-7.4)
  m.add_quad_obj_term(x, x, 1.0)
  m.add_quad_obj_term(y, y, 1.0)
  m.set_obj_constant(19.45)
  m.add_constraint(
    "sum",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    5.0,
  )

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_approx(sol.primal[0], 2.0, 1.0e-7)
  assert_approx(sol.primal[1], 3.0, 1.0e-7)
  assert_approx(sol.objective, 0.65, 1.0e-7)
  assert_true(m.verify_solution(sol, tol=1.0e-7))
}
