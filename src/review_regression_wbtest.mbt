///|
test "tiny LP objective preserves bounded gains and unbounded directions" {
  let width = 17592186044416.0
  let c = 1.0 / width
  let box = Model::new(Maximize)
  let _ = box.add_continuous_var("x", obj=c, ub=width)
  let _ = box.add_continuous_var("y", obj=1.0, ub=1.0)
  let best = box.solve()
  assert_eq(best.status, Optimal)
  assert_eq(best.objective, 2.0)
  assert_true(box.certify_lp_solution(best).certified_optimal)
  for sign in [1.0, -1.0] {
    let ray = Model::new(if sign > 0.0 { Maximize } else { Minimize })
    let _ = ray.add_continuous_var("x", obj=sign * c)
    assert_eq(ray.solve().status, Unbounded)
    let fake = {
      ..best,
      primal: [0.0],
      dual: [],
      reduced_costs: [sign * c],
      objective: 0.0,
    }
    assert_false(ray.certify_lp_solution(fake).certified_optimal)
  }
}

///|
test "dual entry handles scaled rows and tiny prices" {
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("x", obj=1.0)
  m.add_constraint("demand", [{ var_id: x, coeff: 1.0e-14, }], Ge, 1.0)
  let r = solve_lp_dual(m)
  assert_eq(r.status, Optimal)
  assert_approx(r.objective / 1.0e14, 1.0, 1.0e-10)
  assert_true(m.certify_lp_solution(r).certified_optimal)
  let tiny = Model::new(Minimize)
  let a = tiny.add_continuous_var("a", obj=1.0e-16)
  let b = tiny.add_continuous_var("b", obj=1.0e-18)
  tiny.add_constraint(
    "choose",
    [{ var_id: a, coeff: 2.0, }, { var_id: b, coeff: 1.0, }],
    Ge,
    1.0,
  )
  for sol in [solve_lp(tiny), solve_lp_dual(tiny)] {
    assert_eq(sol.status, Optimal)
    assert_approx(sol.objective / 1.0e-18, 1.0, 1.0e-10)
    assert_true(tiny.certify_lp_solution(sol).certified_optimal)
  }
}

///|
test "near-zero cut bound retains incumbent and respects remaining nodes" {
  let m = Model::new(Maximize)
  let x = m.add_binary_var("x", obj=1.0)
  m.add_constraint("half", [{ var_id: x, coeff: 2.0, }], Le, 1.0)
  for budget in [1, 2, 24] {
    let cfg = { ..SolverConfig::default(), max_nodes: budget, }
    let (s, stats) = solve_milp_branch_and_cut(m, config=cfg)
    assert_true(stats.cuts_added > 0)
    assert_eq(s.nodes_explored, stats.nodes_explored)
    assert_true(s.nodes_explored <= budget)
    assert_eq(s.iterations, stats.total_pivots)
    assert_eq(s.primal, [0.0])
    assert_eq(s.objective, 0.0)
    if budget == 24 {
      assert_eq(s.status, Optimal)
      assert_eq(s.mip_gap, 0.0)
    } else {
      assert_true(s.status == Optimal || s.status == NodeLimit)
    }
  }
}

///|
test "MPS section keywords survive as column and vector names" {
  for
    name in [
      "NAME", "OBJSENSE", "ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS", "ENDATA",
    ] {
    let m = Model::new(Maximize)
    let _ = m.add_continuous_var(name, ub=5.0, obj=1.0)
    let parsed = parse_mps(m.to_mps_string())
    assert_eq(parsed.vars.length(), 1)
    assert_eq(parsed.vars[0].name, name)
    assert_eq(parsed.solve().objective, 5.0)
  }
  let model = parse_mps(
    "NAME P\nOBJSENSE MAX\nROWS\n N OBJ\n L CAP\nCOLUMNS\nX OBJ 1 CAP 1\nRHS\nRHS CAP 5\nENDATA\n",
  )
  assert_eq(model.solve().objective, 5.0)
}
