///|
/// Automated executable tests proving real interoperability, extension, and exact
/// capability boundaries between `Freon793/moonopt@0.6.0` and `spring1219/moonlp`.

///|
test "1. Exact linear equivalence & MPS handoff between moonopt and moonlp" {
  let mopt = build_moonopt_production_lp()
  let s_mopt = @moonopt.solve(mopt)
  assert_eq(s_mopt.status, @moonopt.SolveStatus::Optimal)

  // In-memory conversion via from_moonopt()
  let mlp_mem = from_moonopt(mopt)
  let s_mem = mlp_mem.solve()
  assert_eq(s_mem.status, @moonlp.Status::Optimal)
  assert_true((s_mopt.objective - s_mem.objective).abs() < 1.0e-7)

  // MPS text handoff via @mopt_format.write_mps -> @moonlp.parse_mps
  let mps = @mopt_format.write_mps(mopt)
  let mlp_mps = @moonlp.parse_mps(mps, sense=Maximize)
  let s_mps = mlp_mps.solve()
  assert_eq(s_mps.status, @moonlp.Status::Optimal)
  assert_true((s_mopt.objective - s_mps.objective).abs() < 1.0e-7)

  // Extending into 95% Gaussian Chance-Constrained SOCP lowers optimistic profit safely
  let (s_cc, cert_cc) = @moonlp.build_chance_constrained_socp(z_alpha=1.644854)
    .solve_socp_certified()
  assert_eq(s_cc.status, @moonlp.Status::Optimal)
  assert_true(cert_cc.certified_optimal)
  assert_true(s_cc.objective < s_mopt.objective)
}

///|
test "2. moonopt rejects .qps QUADOBJ and CPLEX QCP .lp while moonlp solves & certifies both" {
  // 2a. .qps with QUADOBJ section
  let mk = @moonlp.build_markowitz_portfolio_qp(target_return=0.08)
  let qps = mk.to_mps_string(name="MARKOWITZ")
  assert_true(@mopt_format.parse_mps(qps) is Err(_))
  let mk_rt = @moonlp.parse_mps(qps)
  let mk_sol = mk_rt.solve()
  assert_eq(mk_sol.status, @moonlp.Status::Optimal)
  assert_true(mk_rt.certify_qp_solution(mk_sol, tol=1.0e-6).certified_optimal)

  // 2b. CPLEX .lp with quadratic constraint `[ x ^ 2 + y ^ 2 ] <= 25` (QCP)
  let qcp_lp =
    #|Maximize
    #|  obj: 3 x + 4 y
    #|Subject To
    #|  disk: [ x ^ 2 + y ^ 2 ] <= 25
    #|Bounds
    #|  x free
    #|  y free
    #|End
  // moonopt is linear-only and cannot represent or solve the circular QCP to optimum 25.0
  let mopt_cannot_solve_qcp = match @mopt_format.read_lp(qcp_lp) {
    Err(_) => true
    Ok(m_bad) => {
      let s_bad = @moonopt.solve(m_bad)
      s_bad.status != @moonopt.SolveStatus::Optimal ||
        (s_bad.objective - 25.0).abs() > 1.0
    }
  }
  assert_true(mopt_cannot_solve_qcp)
  let qcp_m = @moonlp.parse_lp(qcp_lp)
  let (qcp_sol, qcp_cert) = qcp_m.solve_socp_certified()
  assert_eq(qcp_sol.status, @moonlp.Status::Optimal)
  assert_true(qcp_cert.certified_optimal)
  assert_true((qcp_sol.objective - 25.0).abs() < 1.0e-5)
}

///|
test "3. Direct extension of moonopt Unit Commitment MILP into Quadratic Fuel-Cost MIQP" {
  let mopt_uc = build_moonopt_uc_milp()
  let s_lin = @moonopt.solve(mopt_uc)
  assert_eq(s_lin.status, @moonopt.SolveStatus::Optimal)

  // Extend the exact moonopt MILP in-place with quadratic heat-rate curves a_g * p_g^2
  let uc_miqp = from_moonopt(mopt_uc)
  let quad_a = [0.008, 0.012, 0.018, 0.025, 0.032]
  for g = 0; g < 5; g = g + 1 {
    uc_miqp.add_quad_obj_term({ index: g, }, { index: g, }, quad_a[g])
  }
  let s_miqp = uc_miqp.solve()
  assert_eq(s_miqp.status, @moonlp.Status::Optimal)
  assert_true(uc_miqp.verify_solution(s_miqp, tol=1.0e-5))
  assert_true((s_miqp.objective - 3308.3).abs() < 1.0e-3)
  // True quadratic cost strictly exceeds linearized moonopt MILP cost
  assert_true(s_miqp.objective > s_lin.objective)
}
