///|
/// Official MoonBit `moon bench --release` performance benchmark suite.
/// Measures authentic Netlib/HiGHS LP instances and Cold-Start vs Dual Simplex Warm-Start MILP.

///|
test "bench netlib afiro LP (28x32)" (b : @bench.T) {
  let m = try {
    parse_mps(netlib_afiro_mps())
  } catch {
    _ => abort("afiro parse failed")
  }
  b.bench(fn() { b.keep(m.solve().objective) })
}

///|
test "bench highs avgas LP (11x8)" (b : @bench.T) {
  let m = try {
    parse_mps(highs_avgas_mps())
  } catch {
    _ => abort("avgas parse failed")
  }
  b.bench(fn() { b.keep(m.solve().objective) })
}

///|
test "bench 10-item knapsack MILP: cold vs warm dual simplex B&B" (b : @bench.T) {
  let vals = [15.0, 22.0, 14.0, 26.0, 32.0, 9.0, 16.0, 8.0, 25.0, 19.0]
  let wts = [4.0, 7.0, 5.0, 8.0, 10.0, 3.0, 5.0, 2.0, 8.0, 6.0]
  let m = Model::new(Maximize)
  let terms : Array[Term] = []
  for i = 0; i < 10; i = i + 1 {
    let v = m.add_binary_var("k_" + i.to_string(), obj=vals[i])
    terms.push({ var_id: v, coeff: wts[i], })
  }
  m.add_constraint("cap", terms, Le, 23.0)
  b.bench(name="knapsack10_cold_bb", fn() { b.keep(solve_milp(m).objective) })
  b.bench(name="knapsack10_warm_dual_bb", fn() {
    b.keep(solve_milp_warm(m).objective)
  })
}

///|
test "bench 5-city TSP-MTZ MILP: cold vs warm dual simplex B&B" (b : @bench.T) {
  let m = build_tsp5_mtz_model()
  b.bench(name="tsp5_cold_bb", fn() { b.keep(solve_milp(m).objective) })
  b.bench(name="tsp5_warm_dual_bb", fn() { b.keep(solve_milp_warm(m).objective) })
}

///|
test "bench netlib adlittle LP (56x97, 480 nonzeros)" (b : @bench.T) {
  let m = try {
    parse_mps(netlib_adlittle_mps())
  } catch {
    _ => abort("adlittle parse failed")
  }
  b.bench(fn() { b.keep(m.solve().objective) })
}

///|
test "bench 5-city TSP-MTZ: unified Branch-and-Cut" (b : @bench.T) {
  let m = build_tsp5_mtz_model()
  b.bench(fn() { b.keep(solve_milp_branch_and_cut(m).0.objective) })
}

///|
test "bench Markowitz 6-asset portfolio QP (Primal-Dual IPM)" (b : @bench.T) {
  let m = build_markowitz_portfolio_qp(target_return=0.08)
  b.bench(fn() { b.keep(solve_qp(m).objective) })
}

///|
test "bench Cardinality-constrained Markowitz MIQP (max 3 of 6 assets)" (
  b : @bench.T,
) {
  let m = build_cardinality_portfolio_miqp(
    max_assets=3,
    min_weight=0.1,
    max_weight=0.5,
    target_return=0.08,
  )
  b.bench(fn() { b.keep(solve_miqp(m).objective) })
}

///|
test "bench Constrained Linear MPC QP (8-step horizon, 24 vars, 16 eqs)" (
  b : @bench.T,
) {
  let m = build_linear_mpc_qp(horizon=8, x0_pos=2.0, x0_vel=0.0, u_max=1.0)
  b.bench(fn() { b.keep(solve_qp(m).objective) })
}

///|
test "bench Soft-Margin SVM Dual QP (8x8 rank-2 PSD)" (b : @bench.T) {
  let m = build_svm_dual_qp(c_reg=10.0)
  b.bench(fn() { b.keep(solve_qp(m).objective) })
}
