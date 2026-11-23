///|
/// Comprehensive verification test suite for moonlp.

///|
fn assert_approx(actual : Double, expected : Double, tol : Double) -> Unit raise {
  if (actual - expected).abs() > tol {
    fail("Expected \{expected}, got \{actual} (tol=\{tol})")
  }
}

///|
test "classic 2D production mix LP maximize with exact dual prices" {
  // Maximize 3*x + 2*y
  // s.t.  x + y <= 4
  //       2*x + y <= 5
  //       x, y >= 0
  // Optimum: x = 1, y = 3, obj = 9
  // Dual shadow prices: y1 = 1, y2 = 1 (strong duality: 4*1 + 5*1 = 9)
  let m = Model::new(Maximize)
  let x = m.add_continuous_var("x", obj=3.0)
  let y = m.add_continuous_var("y", obj=2.0)
  m.add_constraint("labor", [{ var_id: x, coeff: 1.0 }, { var_id: y, coeff: 1.0 }], Le, 4.0)
  m.add_constraint("machine", [{ var_id: x, coeff: 2.0 }, { var_id: y, coeff: 1.0 }], Le, 5.0)

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 9.0, 1.0e-7)
  assert_approx(sol.primal[x.index], 1.0, 1.0e-7)
  assert_approx(sol.primal[y.index], 3.0, 1.0e-7)
  assert_approx(sol.dual[0], 1.0, 1.0e-7)
  assert_approx(sol.dual[1], 1.0, 1.0e-7)
}

///|
test "diet problem LP minimize with >= constraints and strong duality" {
  // Minimize 0.6*x1 + 0.35*x2
  // s.t.  5*x1 + 7*x2 >= 8
  //       4*x1 + 2*x2 >= 15
  //       2*x1 + 1*x2 >= 3
  //       x1, x2 >= 0
  // Binding: 4*x1 + 2*x2 = 15 (x1=3.75, x2=0, obj=2.25)
  let m = Model::new(Minimize)
  let x1 = m.add_continuous_var("x1", obj=0.6)
  let x2 = m.add_continuous_var("x2", obj=0.35)
  m.add_constraint("vitA", [{ var_id: x1, coeff: 5.0 }, { var_id: x2, coeff: 7.0 }], Ge, 8.0)
  m.add_constraint("vitB", [{ var_id: x1, coeff: 4.0 }, { var_id: x2, coeff: 2.0 }], Ge, 15.0)
  m.add_constraint("vitC", [{ var_id: x1, coeff: 2.0 }, { var_id: x2, coeff: 1.0 }], Ge, 3.0)

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 2.25, 1.0e-7)
  assert_approx(sol.primal[x1.index], 3.75, 1.0e-7)
  assert_approx(sol.primal[x2.index], 0.0, 1.0e-7)
  // Strong duality: b^T y == c^T x == 2.25
  let dual_obj = 8.0 * sol.dual[0] + 15.0 * sol.dual[1] + 3.0 * sol.dual[2]
  assert_approx(dual_obj, sol.objective, 1.0e-6)
}

///|
test "free variables and negative lower bounds" {
  // Minimize 2*x - 3*y
  // s.t. x + y == 1
  //      -4 <= x <= 2
  //      y free (-inf .. +inf)
  // Since y = 1 - x, obj = 2*x - 3*(1-x) = 5*x - 3 -> minimized at x = -4, y = 5, obj = -23
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("x", lb=-4.0, ub=2.0, obj=2.0)
  let y = m.add_continuous_var("y", lb=-1.0e300, ub=1.0e300, obj=-3.0)
  m.add_constraint("eq1", [{ var_id: x, coeff: 1.0 }, { var_id: y, coeff: 1.0 }], Eq, 1.0)

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.primal[x.index], -4.0, 1.0e-7)
  assert_approx(sol.primal[y.index], 5.0, 1.0e-7)
  assert_approx(sol.objective, -23.0, 1.0e-7)
}

///|
test "infeasible LP detection" {
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("x", obj=1.0)
  let y = m.add_continuous_var("y", obj=1.0)
  m.add_constraint("c1", [{ var_id: x, coeff: 1.0 }, { var_id: y, coeff: 1.0 }], Le, 2.0)
  m.add_constraint("c2", [{ var_id: x, coeff: 1.0 }, { var_id: y, coeff: 1.0 }], Ge, 5.0)
  let sol = m.solve()
  assert_eq(sol.status, Infeasible)
}

///|
test "unbounded LP detection" {
  let m = Model::new(Maximize)
  let x = m.add_continuous_var("x", obj=2.0)
  let y = m.add_continuous_var("y", obj=1.0)
  m.add_constraint("c1", [{ var_id: x, coeff: -1.0 }, { var_id: y, coeff: 1.0 }], Le, 1.0)
  let sol = m.solve()
  assert_eq(sol.status, Unbounded)
}

///|
test "Beale degenerate cycling LP solved via Bland anti-cycling rule" {
  // Beale's classic cycling example (1955):
  // Minimize -0.75*x1 + 20*x2 - 0.5*x3 + 6*x4
  // s.t.  0.25*x1 - 8*x2 - x3 + 9*x4 <= 0
  //       0.5*x1 - 12*x2 - 0.5*x3 + 3*x4 <= 0
  //       x3 <= 1
  //       x1, x2, x3, x4 >= 0
  // Exact global optimum: obj = -1.25 (x1 = 1, x2 = 0, x3 = 1, x4 = 0)
  let m = Model::new(Minimize)
  let x1 = m.add_continuous_var("x1", obj=-0.75)
  let x2 = m.add_continuous_var("x2", obj=20.0)
  let x3 = m.add_continuous_var("x3", obj=-0.5)
  let x4 = m.add_continuous_var("x4", obj=6.0)
  m.add_constraint(
    "r1",
    [
      { var_id: x1, coeff: 0.25 },
      { var_id: x2, coeff: -8.0 },
      { var_id: x3, coeff: -1.0 },
      { var_id: x4, coeff: 9.0 },
    ],
    Le,
    0.0,
  )
  m.add_constraint(
    "r2",
    [
      { var_id: x1, coeff: 0.5 },
      { var_id: x2, coeff: -12.0 },
      { var_id: x3, coeff: -0.5 },
      { var_id: x4, coeff: 3.0 },
    ],
    Le,
    0.0,
  )
  m.add_constraint("r3", [{ var_id: x3, coeff: 1.0 }], Le, 1.0)

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, -1.25, 1.0e-7)
}

///|
test "0-1 Knapsack MILP via Branch-and-Bound" {
  // 5 items: values = [10, 14, 18, 22, 9], weights = [3, 4, 5, 7, 2], capacity = 11
  // Optimal subset: items 0 (w=3,v=10) + 1 (w=4,v=14) + 3? w=14>11
  // Check: 1 (w=4,v=14) + 2 (w=5,v=18) + 4 (w=2,v=9) -> w=11, v=41!
  let m = Model::new(Maximize)
  let vals = [10.0, 14.0, 18.0, 22.0, 9.0]
  let wts = [3.0, 4.0, 5.0, 7.0, 2.0]
  let terms : Array[Term] = []
  for i = 0; i < 5; i = i + 1 {
    let v = m.add_binary_var("item_" + i.to_string(), obj=vals[i])
    terms.push({ var_id: v, coeff: wts[i] })
  }
  m.add_constraint("capacity", terms, Le, 11.0)

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 41.0, 1.0e-6)
}

///|
test "general integer programming with non-trivial LP relaxation gap" {
  // Maximize 5*x + 8*y
  // s.t. x + y <= 6
  //      5*x + 9*y <= 45
  //      x, y >= 0, integer
  // LP relaxation: x=2.25, y=3.75 -> obj=41.25
  // Integer optimum: x=0, y=5 -> obj=40
  let m = Model::new(Maximize)
  let x = m.add_int_var("x", obj=5.0)
  let y = m.add_int_var("y", obj=8.0)
  m.add_constraint("c1", [{ var_id: x, coeff: 1.0 }, { var_id: y, coeff: 1.0 }], Le, 6.0)
  m.add_constraint("c2", [{ var_id: x, coeff: 5.0 }, { var_id: y, coeff: 9.0 }], Le, 45.0)

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 40.0, 1.0e-6)
  assert_approx(sol.primal[x.index], 0.0, 1.0e-6)
  assert_approx(sol.primal[y.index], 5.0, 1.0e-6)
}

///|
test "CPLEX .lp exporter smoke test" {
  let m = Model::new(Maximize)
  let x = m.add_binary_var("x", obj=3.0)
  let y = m.add_int_var("y", lb=0.0, ub=10.0, obj=2.0)
  m.add_constraint("cap", [{ var_id: x, coeff: 2.0 }, { var_id: y, coeff: 1.0 }], Le, 7.0)
  let lp = m.to_lp_string()
  assert_true(lp.contains("Maximize"))
  assert_true(lp.contains("Subject To"))
  assert_true(lp.contains("Binaries"))
  assert_true(lp.contains("Generals"))
  assert_true(lp.contains("End"))
}
