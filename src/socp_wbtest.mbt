///|
/// Comprehensive Second-Order Cone (SOCP), Quadratically Constrained (QCQP),
/// and Mixed-Integer Second-Order Cone (MISOCP) Verification Suite (v2.0.0).

///|
test "77. [v2.0 SOCP] Smallest Enclosing Disk (analytical circumcircle + Conic KKT Certificate)" {
  let m = build_smallest_enclosing_disk_socp()
  let (sol, cert) = m.solve_socp_certified()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol, tol=1.0e-6))
  assert_true(cert.certified_optimal)
  assert_true(cert.duality_gap < 1.0e-5)

  // Analytical circumcircle of (0,0), (4,0), (2,3.5):
  // center = (2.0, (3.5^2 - 2^2)/(2*3.5)) = (2.0, 8.25/7.0) = (2.0, 1.1785714285714286)
  // radius = 3.5 - 8.25/7.0 = 16.25/7.0 = 2.3214285714285716
  assert_approx(sol.primal[0], 2.0, 1.0e-6)
  assert_approx(sol.primal[1], 8.25 / 7.0, 1.0e-6)
  assert_approx(sol.primal[2], 16.25 / 7.0, 1.0e-6)
  assert_approx(sol.objective, 16.25 / 7.0, 1.0e-6)
}

///|
test "78. [v2.0 QCQP] Convex quadratic constraint auto-lifted via add_quad_constraint" {
  // max 3 x + 4 y  s.t.  x^2 + y^2 <= 25  (circle of radius 5!)
  // Analytical optimum: (3, 4), obj = 25.0
  let m = Model::new(Maximize)
  let x = m.add_continuous_var("x", lb=-1.0e300, ub=1.0e300, obj=3.0)
  let y = m.add_continuous_var("y", lb=-1.0e300, ub=1.0e300, obj=4.0)
  let ok = m.add_quad_constraint(
    "disk",
    [x, y],
    [[1.0, 0.0], [0.0, 1.0]],
    LinearExpr::zero(),
    25.0,
    half=false,
  )
  assert_true(ok)
  let (sol, cert) = m.solve_socp_certified()
  assert_eq(sol.status, Optimal)
  assert_approx(sol.primal[0], 3.0, 1.0e-6)
  assert_approx(sol.primal[1], 4.0, 1.0e-6)
  assert_approx(sol.objective, 25.0, 1.0e-6)
  assert_true(cert.certified_optimal)

  // Non-convex quadratic constraint must be rejected cleanly by add_quad_constraint
  let bad = m.add_quad_constraint(
    "nonconvex",
    [x, y],
    [[1.0, 0.0], [0.0, -1.0]],
    LinearExpr::zero(),
    1.0,
    half=false,
  )
  assert_false(bad)
}

///|
test "79. [v2.0 SOCP] Disjoint Second-Order Cones detected as Infeasible via Conic Farkas" {
  // ||(x, y) - (0, 0)||_2 <= 1  AND  ||(x, y) - (5, 0)||_2 <= 1  (linearly feasible, conically infeasible!)
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("x", lb=-10.0, ub=10.0, obj=1.0)
  let y = m.add_continuous_var("y", lb=-10.0, ub=10.0)
  m.add_soc_constraint(
    "ball1",
    [LinearExpr::from_var(x), LinearExpr::from_var(y)],
    LinearExpr::zero().add_const(1.0),
  )
  m.add_soc_constraint(
    "ball2",
    [LinearExpr::from_var(x).add_const(-5.0), LinearExpr::from_var(y)],
    LinearExpr::zero().add_const(1.0),
  )
  assert_eq(m.solve().status, Infeasible)
}

///|
test "80. [v2.0 SOCP Model] Robust Ellipsoidal-Uncertainty Portfolio Optimization" {
  let m = build_robust_portfolio_socp(kappa=1.5, max_weight=0.5)
  let (sol, cert) = m.solve_socp_certified()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol, tol=1.0e-6))
  assert_true(cert.certified_optimal)
  assert_true(cert.duality_gap < 1.0e-5)

  // Tighter uncertainty (higher kappa) must yield lower robust worst-case return
  let m_strict = build_robust_portfolio_socp(kappa=2.5, max_weight=0.5)
  let s_strict = m_strict.solve()
  assert_eq(s_strict.status, Optimal)
  assert_true(s_strict.objective < sol.objective)
}

///|
test "81. [v2.0 MISOCP Model] Volatility-Capped Cardinality-Constrained Portfolio (max 3 assets, vol <= 10.5%)" {
  let m = build_vol_capped_cardinality_misocp(
    max_assets=3,
    vol_cap=0.105,
    min_weight=0.1,
    max_weight=0.5,
  )
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol, tol=1.0e-5))
  assert_true(sol.mip_gap <= 1.0e-4)

  let mut active = 0
  let mut wsum = 0.0
  for i = 0; i < 6; i = i + 1 {
    wsum = wsum + sol.primal[i]
    if sol.primal[6 + i] > 0.5 {
      active = active + 1
      assert_true(sol.primal[i] >= 0.1 - 1.0e-5)
    } else {
      assert_true(sol.primal[i].abs() <= 1.0e-5)
    }
  }
  assert_true(active <= 3)
  assert_approx(wsum, 1.0, 1.0e-5)
  assert_true(sol.objective > 0.08)
}

///|
test "82. [v2.0 SOCP Model] 95% Gaussian Chance-Constrained Production Planning" {
  let m = build_chance_constrained_socp(z_alpha=1.644854)
  let (sol, cert) = m.solve_socp_certified()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol, tol=1.0e-5))
  assert_true(cert.certified_optimal)

  // Deterministic (z_alpha = 0) ignores uncertainty and achieves strictly higher optimistic profit
  let m_det = build_chance_constrained_socp(z_alpha=0.0)
  let s_det = m_det.solve()
  assert_eq(s_det.status, Optimal)
  assert_true(sol.objective < s_det.objective)
}

///|
test "83. [v2.0 SOCP & MISOCP Model] 2D Isotropic Circular-Thrust Rocket Guidance (continuous & switched)" {
  // 1) Continuous 2D isotropic circular-thrust MPC (quad objective + 4 Lorentz cones)
  let m_cont = build_isotropic_thrust_mpc_socp(horizon=4, u_max=1.0)
  let s_cont = m_cont.solve()
  assert_eq(s_cont.status, Optimal)
  assert_true(m_cont.verify_solution(s_cont, tol=1.0e-5))
  // Verify exact circular norm ||(ux_t, uy_t)||_2 <= 1.0 at every step
  for t = 0; t < 4; t = t + 1 {
    let ux = s_cont.primal[16 + 2 * t]
    let uy = s_cont.primal[16 + 2 * t + 1]
    assert_true(d_sqrt(ux * ux + uy * uy) <= 1.0 + 1.0e-5)
  }

  // 2) Switched 2D circular-thrust MISOCP (at most 2 engine burns out of 4 steps!)
  let m_sw = build_switched_isotropic_mpc_misocp(
    horizon=4,
    max_burns=2,
    u_max=1.0,
  )
  let s_sw = m_sw.solve()
  assert_eq(s_sw.status, Optimal)
  assert_true(m_sw.verify_solution(s_sw, tol=1.0e-5))
  assert_true(s_sw.mip_gap <= 1.0e-4)
  let mut burns = 0
  for t = 0; t < 4; t = t + 1 {
    let ux = s_sw.primal[16 + 2 * t]
    let uy = s_sw.primal[16 + 2 * t + 1]
    let zt = s_sw.primal[24 + t]
    let nrm = d_sqrt(ux * ux + uy * uy)
    if zt > 0.5 {
      burns = burns + 1
      assert_true(nrm <= 1.0 + 1.0e-5)
    } else {
      assert_true(nrm <= 1.0e-5)
    }
  }
  assert_true(burns <= 2 && burns >= 1)
  // Restricted 2-burn switched cost must exceed unrestricted 4-burn continuous cost
  assert_true(s_sw.objective > s_cont.objective)
}
