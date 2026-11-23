///|
/// Comprehensive verification test suite for moonlp (15 tests covering LP, Duals,
/// Strong Duality, Beale Cycling, Klee-Minty 4D, Network Max-Flow, Knapsack,
/// General IP, Assignment, Facility Location, 4x4 Sudoku MILP, CPLEX .lp Parser, and Presolve).

///|
fn assert_approx(
  actual : Double,
  expected : Double,
  tol : Double,
) -> Unit raise {
  if (actual - expected).abs() > tol {
    fail("Expected \{expected}, got \{actual} (tol=\{tol})")
  }
}

///|
test "1. classic 2D production mix LP maximize with exact dual shadow prices" {
  let m = Model::new(Maximize)
  let x = m.add_continuous_var("x", obj=3.0)
  let y = m.add_continuous_var("y", obj=2.0)
  m.add_constraint(
    "labor",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    4.0,
  )
  m.add_constraint(
    "machine",
    [{ var_id: x, coeff: 2.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    5.0,
  )

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
test "2. diet problem LP minimize with >= constraints and strong duality" {
  let m = Model::new(Minimize)
  let x1 = m.add_continuous_var("x1", obj=0.6)
  let x2 = m.add_continuous_var("x2", obj=0.35)
  m.add_constraint(
    "vitA",
    [{ var_id: x1, coeff: 5.0, }, { var_id: x2, coeff: 7.0, }],
    Ge,
    8.0,
  )
  m.add_constraint(
    "vitB",
    [{ var_id: x1, coeff: 4.0, }, { var_id: x2, coeff: 2.0, }],
    Ge,
    15.0,
  )
  m.add_constraint(
    "vitC",
    [{ var_id: x1, coeff: 2.0, }, { var_id: x2, coeff: 1.0, }],
    Ge,
    3.0,
  )

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 2.25, 1.0e-7)
  assert_approx(sol.primal[x1.index], 3.75, 1.0e-7)
  assert_approx(sol.primal[x2.index], 0.0, 1.0e-7)
  let dual_obj = 8.0 * sol.dual[0] + 15.0 * sol.dual[1] + 3.0 * sol.dual[2]
  assert_approx(dual_obj, sol.objective, 1.0e-6)
}

///|
test "3. free variables and negative lower bounds" {
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("x", lb=-4.0, ub=2.0, obj=2.0)
  let y = m.add_continuous_var("y", lb=-1.0e300, ub=1.0e300, obj=-3.0)
  m.add_constraint(
    "eq1",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Eq,
    1.0,
  )

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.primal[x.index], -4.0, 1.0e-7)
  assert_approx(sol.primal[y.index], 5.0, 1.0e-7)
  assert_approx(sol.objective, -23.0, 1.0e-7)
}

///|
test "4. infeasible LP detection" {
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("x", obj=1.0)
  let y = m.add_continuous_var("y", obj=1.0)
  m.add_constraint(
    "c1",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    2.0,
  )
  m.add_constraint(
    "c2",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Ge,
    5.0,
  )
  let sol = m.solve()
  assert_eq(sol.status, Infeasible)
}

///|
test "5. unbounded LP detection" {
  let m = Model::new(Maximize)
  let x = m.add_continuous_var("x", obj=2.0)
  let y = m.add_continuous_var("y", obj=1.0)
  m.add_constraint(
    "c1",
    [{ var_id: x, coeff: -1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    1.0,
  )
  let sol = m.solve()
  assert_eq(sol.status, Unbounded)
}

///|
test "6. Beale 1955 degenerate cycling LP solved via Bland anti-cycling rule" {
  let m = Model::new(Minimize)
  let x1 = m.add_continuous_var("x1", obj=-0.75)
  let x2 = m.add_continuous_var("x2", obj=20.0)
  let x3 = m.add_continuous_var("x3", obj=-0.5)
  let x4 = m.add_continuous_var("x4", obj=6.0)
  m.add_constraint(
    "r1",
    [
      { var_id: x1, coeff: 0.25, },
      { var_id: x2, coeff: -8.0, },
      { var_id: x3, coeff: -1.0, },
      { var_id: x4, coeff: 9.0, },
    ],
    Le,
    0.0,
  )
  m.add_constraint(
    "r2",
    [
      { var_id: x1, coeff: 0.5, },
      { var_id: x2, coeff: -12.0, },
      { var_id: x3, coeff: -0.5, },
      { var_id: x4, coeff: 3.0, },
    ],
    Le,
    0.0,
  )
  m.add_constraint("r3", [{ var_id: x3, coeff: 1.0, }], Le, 1.0)

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, -1.25, 1.0e-7)
}

///|
test "7. Klee-Minty 4D cube worst-case LP torture test" {
  let m = Model::new(Maximize)
  let x0 = m.add_continuous_var("x0", obj=8.0)
  let x1 = m.add_continuous_var("x1", obj=4.0)
  let x2 = m.add_continuous_var("x2", obj=2.0)
  let x3 = m.add_continuous_var("x3", obj=1.0)
  m.add_constraint("km0", [{ var_id: x0, coeff: 1.0, }], Le, 5.0)
  m.add_constraint(
    "km1",
    [{ var_id: x0, coeff: 4.0, }, { var_id: x1, coeff: 1.0, }],
    Le,
    25.0,
  )
  m.add_constraint(
    "km2",
    [
      { var_id: x0, coeff: 8.0, },
      { var_id: x1, coeff: 4.0, },
      { var_id: x2, coeff: 1.0, },
    ],
    Le,
    125.0,
  )
  m.add_constraint(
    "km3",
    [
      { var_id: x0, coeff: 16.0, },
      { var_id: x1, coeff: 8.0, },
      { var_id: x2, coeff: 4.0, },
      { var_id: x3, coeff: 1.0, },
    ],
    Le,
    625.0,
  )

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 625.0, 1.0e-6)
}

///|
test "8. network max-flow formulated as LP (CLRS classic 6-node graph = 23)" {
  // Edges: (s,v1,16), (s,v2,13), (v1,v2,10), (v2,v1,4), (v1,v3,12),
  //        (v3,v2,9), (v2,v4,14), (v4,v3,7), (v3,t,20), (v4,t,4)
  // Exact max flow = 23
  let m = Model::new(Maximize)
  let e_s1 = m.add_continuous_var("s1", ub=16.0, obj=1.0)
  let e_s2 = m.add_continuous_var("s2", ub=13.0, obj=1.0)
  let e_12 = m.add_continuous_var("12", ub=10.0)
  let e_21 = m.add_continuous_var("21", ub=4.0)
  let e_13 = m.add_continuous_var("13", ub=12.0)
  let e_32 = m.add_continuous_var("32", ub=9.0)
  let e_24 = m.add_continuous_var("24", ub=14.0)
  let e_43 = m.add_continuous_var("43", ub=7.0)
  let e_3t = m.add_continuous_var("3t", ub=20.0)
  let e_4t = m.add_continuous_var("4t", ub=4.0)

  // Flow conservation at v1, v2, v3, v4
  m.add_constraint(
    "v1",
    [
      { var_id: e_s1, coeff: 1.0, },
      { var_id: e_21, coeff: 1.0, },
      { var_id: e_12, coeff: -1.0, },
      { var_id: e_13, coeff: -1.0, },
    ],
    Eq,
    0.0,
  )
  m.add_constraint(
    "v2",
    [
      { var_id: e_s2, coeff: 1.0, },
      { var_id: e_12, coeff: 1.0, },
      { var_id: e_32, coeff: 1.0, },
      { var_id: e_21, coeff: -1.0, },
      { var_id: e_24, coeff: -1.0, },
    ],
    Eq,
    0.0,
  )
  m.add_constraint(
    "v3",
    [
      { var_id: e_13, coeff: 1.0, },
      { var_id: e_43, coeff: 1.0, },
      { var_id: e_32, coeff: -1.0, },
      { var_id: e_3t, coeff: -1.0, },
    ],
    Eq,
    0.0,
  )
  m.add_constraint(
    "v4",
    [
      { var_id: e_24, coeff: 1.0, },
      { var_id: e_43, coeff: -1.0, },
      { var_id: e_4t, coeff: -1.0, },
    ],
    Eq,
    0.0,
  )

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 23.0, 1.0e-6)
}

///|
test "9. 0-1 Knapsack MILP via Branch-and-Bound" {
  let m = Model::new(Maximize)
  let vals = [10.0, 14.0, 18.0, 22.0, 9.0]
  let wts = [3.0, 4.0, 5.0, 7.0, 2.0]
  let terms : Array[Term] = []
  for i = 0; i < 5; i = i + 1 {
    let v = m.add_binary_var("item_" + i.to_string(), obj=vals[i])
    terms.push({ var_id: v, coeff: wts[i], })
  }
  m.add_constraint("capacity", terms, Le, 11.0)

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 41.0, 1.0e-6)
}

///|
test "10. general integer programming with LP relaxation gap" {
  let m = Model::new(Maximize)
  let x = m.add_int_var("x", obj=5.0)
  let y = m.add_int_var("y", obj=8.0)
  m.add_constraint(
    "c1",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    6.0,
  )
  m.add_constraint(
    "c2",
    [{ var_id: x, coeff: 5.0, }, { var_id: y, coeff: 9.0, }],
    Le,
    45.0,
  )

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 40.0, 1.0e-6)
  assert_approx(sol.primal[x.index], 0.0, 1.0e-6)
  assert_approx(sol.primal[y.index], 5.0, 1.0e-6)
}

///|
test "11. 3x3 combinatorial assignment MILP" {
  let costs = [[9.0, 2.0, 7.0], [6.0, 4.0, 3.0], [5.0, 8.0, 1.0]]
  let m = Model::new(Minimize)
  let grid : Array[Array[VarId]] = []
  for i = 0; i < 3; i = i + 1 {
    let row : Array[VarId] = []
    for j = 0; j < 3; j = j + 1 {
      row.push(
        m.add_binary_var(
          "x_" + i.to_string() + "_" + j.to_string(),
          obj=costs[i][j],
        ),
      )
    }
    grid.push(row)
  }
  for i = 0; i < 3; i = i + 1 {
    m.add_constraint(
      "row_" + i.to_string(),
      [
        { var_id: grid[i][0], coeff: 1.0, },
        { var_id: grid[i][1], coeff: 1.0, },
        { var_id: grid[i][2], coeff: 1.0, },
      ],
      Eq,
      1.0,
    )
  }
  for j = 0; j < 3; j = j + 1 {
    m.add_constraint(
      "col_" + j.to_string(),
      [
        { var_id: grid[0][j], coeff: 1.0, },
        { var_id: grid[1][j], coeff: 1.0, },
        { var_id: grid[2][j], coeff: 1.0, },
      ],
      Eq,
      1.0,
    )
  }
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 9.0, 1.0e-6)
}

///|
test "12. mixed binary + continuous facility location MILP" {
  let m = Model::new(Minimize)
  let y0 = m.add_binary_var("open0", obj=10.0)
  let y1 = m.add_binary_var("open1", obj=12.0)
  let ship = [[4.0, 6.0, 9.0], [8.0, 3.0, 2.0]]
  let f : Array[Array[VarId]] = []
  for w = 0; w < 2; w = w + 1 {
    let row : Array[VarId] = []
    for c = 0; c < 3; c = c + 1 {
      row.push(
        m.add_continuous_var(
          "f_" + w.to_string() + "_" + c.to_string(),
          ub=1.0,
          obj=ship[w][c],
        ),
      )
    }
    f.push(row)
  }
  for c = 0; c < 3; c = c + 1 {
    m.add_constraint(
      "dem_" + c.to_string(),
      [{ var_id: f[0][c], coeff: 1.0, }, { var_id: f[1][c], coeff: 1.0, }],
      Eq,
      1.0,
    )
  }
  let opens = [y0, y1]
  for w = 0; w < 2; w = w + 1 {
    for c = 0; c < 3; c = c + 1 {
      m.add_constraint(
        "link_" + w.to_string() + "_" + c.to_string(),
        [{ var_id: f[w][c], coeff: 1.0, }, { var_id: opens[w], coeff: -1.0, }],
        Le,
        0.0,
      )
    }
  }
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 25.0, 1.0e-6)
}

///|
test "13. 4x4 Sudoku solved as 64-binary-variable MILP" {
  // Puzzle (0 = empty):
  // 1 0 | 0 4
  // 0 0 | 1 0
  // ----+----
  // 0 1 | 0 0
  // 4 0 | 0 1
  let m = Model::new(Minimize)
  let vars : Array[VarId] = []
  for r = 0; r < 4; r = r + 1 {
    for c = 0; c < 4; c = c + 1 {
      for d = 0; d < 4; d = d + 1 {
        vars.push(
          m.add_binary_var(
            "s_" + r.to_string() + "_" + c.to_string() + "_" + d.to_string(),
          ),
        )
      }
    }
  }
  let idx = fn(r : Int32, c : Int32, d : Int32) -> VarId {
    vars[(r * 4 + c) * 4 + d]
  }
  // Each cell gets exactly one digit
  for r = 0; r < 4; r = r + 1 {
    for c = 0; c < 4; c = c + 1 {
      let terms : Array[Term] = []
      for d = 0; d < 4; d = d + 1 {
        terms.push({ var_id: idx(r, c, d), coeff: 1.0, })
      }
      m.add_constraint("cell", terms, Eq, 1.0)
    }
  }
  // Each row & col has each digit once
  for d = 0; d < 4; d = d + 1 {
    for r = 0; r < 4; r = r + 1 {
      let rt : Array[Term] = []
      for c = 0; c < 4; c = c + 1 {
        rt.push({ var_id: idx(r, c, d), coeff: 1.0, })
      }
      m.add_constraint("row", rt, Eq, 1.0)
    }
    for c = 0; c < 4; c = c + 1 {
      let ct : Array[Term] = []
      for r = 0; r < 4; r = r + 1 {
        ct.push({ var_id: idx(r, c, d), coeff: 1.0, })
      }
      m.add_constraint("col", ct, Eq, 1.0)
    }
    // Each 2x2 box has each digit once
    for br = 0; br < 2; br = br + 1 {
      for bc = 0; bc < 2; bc = bc + 1 {
        let bt : Array[Term] = []
        for dr = 0; dr < 2; dr = dr + 1 {
          for dc = 0; dc < 2; dc = dc + 1 {
            bt.push({ var_id: idx(br * 2 + dr, bc * 2 + dc, d), coeff: 1.0, })
          }
        }
        m.add_constraint("box", bt, Eq, 1.0)
      }
    }
  }
  // Clues: (0,0)=1, (0,3)=4, (1,2)=1, (2,1)=1, (3,0)=4, (3,3)=1
  m.set_var_bounds(idx(0, 0, 0), 1.0, 1.0)
  m.set_var_bounds(idx(0, 3, 3), 1.0, 1.0)
  m.set_var_bounds(idx(1, 2, 0), 1.0, 1.0)
  m.set_var_bounds(idx(2, 1, 0), 1.0, 1.0)
  m.set_var_bounds(idx(3, 0, 3), 1.0, 1.0)
  m.set_var_bounds(idx(3, 3, 0), 1.0, 1.0)

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
}

///|
test "14. CPLEX .lp parser & round-trip solve" {
  let lp_text =
    #|\\ Classic MILP in standard CPLEX .lp format
    #|Maximize
    #|  obj: 5 x + 8 y
    #|Subject To
    #|  c1: x + y <= 6
    #|  c2: 5 x + 9 y <= 45
    #|Bounds
    #|  0 <= x <= 10
    #|  0 <= y <= 10
    #|Generals
    #|  x y
    #|End
  let parsed = parse_lp(lp_text)
  let sol1 = parsed.solve()
  assert_eq(sol1.status, Optimal)
  assert_approx(sol1.objective, 40.0, 1.0e-6)

  // Round-trip: export -> re-parse -> solve again
  let exported = parsed.to_lp_string()
  let reparsed = parse_lp(exported)
  let sol2 = reparsed.solve()
  assert_eq(sol2.status, Optimal)
  assert_approx(sol2.objective, 40.0, 1.0e-6)
}

///|
test "15. presolve singleton bound tightening & infeasibility shortcut" {
  let m = Model::new(Maximize)
  let x = m.add_int_var("x", lb=0.0, ub=100.0, obj=3.0)
  let y = m.add_int_var("y", lb=0.0, ub=100.0, obj=2.0)
  // Singleton constraints that tighten x <= 3.7 -> 3 and y >= 2.2 -> 3
  m.add_constraint("s1", [{ var_id: x, coeff: 1.0, }], Le, 3.7)
  m.add_constraint("s2", [{ var_id: y, coeff: 1.0, }], Ge, 2.2)
  m.add_constraint(
    "sum",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    7.0,
  )

  let (reduced, report) = m.presolve()
  assert_true(!report.infeasible_detected)
  assert_eq(report.bounds_tightened, 2)
  assert_eq(report.redundant_rows_removed, 2)
  assert_eq(reduced.num_constraints(), 1)

  let sol = m.solve_presolved()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  // Max 3*x + 2*y with x <= 3, y >= 3, x + y <= 7 -> x=3, y=4 -> obj = 17
  assert_approx(sol.objective, 17.0, 1.0e-6)
}

///|
test "16. Netlib MPS (.mps) standard format parser & solve (LP)" {
  let mps_text =
    #|* Classic Netlib-style MPS LP example
    #|NAME          TESTPROB
    #|ROWS
    #| N  COST
    #| L  LIM1
    #| G  LIM2
    #| E  MYEQN
    #|COLUMNS
    #|    XONE      COST                 1.0   LIM1                 1.0
    #|    XONE      LIM2                 1.0
    #|    YTWO      COST                 4.0   LIM1                 1.0
    #|    YTWO      MYEQN               -1.0
    #|    ZTHREE    COST                 9.0   LIM2                 1.0
    #|    ZTHREE    MYEQN                1.0
    #|RHS
    #|    RHS1      LIM1                 5.0   LIM2                10.0
    #|    RHS1      MYEQN                7.0
    #|BOUNDS
    #| UP BND1      XONE                 4.0
    #| LO BND1      YTWO                -1.0
    #| UP BND1      YTWO                 1.0
    #|ENDATA
  // Minimize x + 4*y + 9*z
  // s.t. x + y <= 5, x + z >= 10, -y + z == 7
  //      0 <= x <= 4, -1 <= y <= 1, z >= 0
  // Since z = 7 + y and x + z = x + y + 7 >= 10 -> x + y >= 3.
  // Obj = x + 4*y + 9*(7 + y) = 63 + x + 13*y.
  // To minimize, make y as small as possible: y = -1 -> z = 6, x + z >= 10 -> x >= 4 -> x = 4!
  // Check: x=4, y=-1, z=6 -> x+y=3<=5, x+z=10>=10, -y+z=7==7 -> obj = 4 - 4 + 54 = 54.0!
  let m = parse_mps(mps_text)
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 54.0, 1.0e-6)
}

///|
test "17. Netlib MPS Mixed-Integer (INTORG/INTEND) parser & solve" {
  let mps_text =
    #|NAME          MIPDEMO
    #|ROWS
    #| N  OBJ
    #| L  C1
    #| L  C2
    #|COLUMNS
    #|    MARK0     'MARKER'                 'INTORG'
    #|    X         OBJ                 -5.0   C1                   1.0
    #|    X         C2                   5.0
    #|    Y         OBJ                 -8.0   C1                   1.0
    #|    Y         C2                   9.0
    #|    MARK1     'MARKER'                 'INTEND'
    #|RHS
    #|    RHS1      C1                   6.0   C2                  45.0
    #|ENDATA
  // Minimize -5x - 8y s.t. x + y <= 6, 5x + 9y <= 45, x,y >= 0 integer -> opt = -40.0
  let m = parse_mps(mps_text)
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, -40.0, 1.0e-6)
}

///|
test "18. KKT strong duality & complementary slackness optimality certificate" {
  let m = Model::new(Minimize)
  let x1 = m.add_continuous_var("x1", obj=0.6)
  let x2 = m.add_continuous_var("x2", obj=0.35)
  m.add_constraint(
    "vitA",
    [{ var_id: x1, coeff: 5.0, }, { var_id: x2, coeff: 7.0, }],
    Ge,
    8.0,
  )
  m.add_constraint(
    "vitB",
    [{ var_id: x1, coeff: 4.0, }, { var_id: x2, coeff: 2.0, }],
    Ge,
    15.0,
  )
  m.add_constraint(
    "vitC",
    [{ var_id: x1, coeff: 2.0, }, { var_id: x2, coeff: 1.0, }],
    Ge,
    3.0,
  )

  let sol = m.solve()
  let cert = m.certify_lp_solution(sol)
  assert_true(cert.certified_optimal)
  assert_true(cert.primal_feasible)
  assert_true(cert.duality_gap < 1.0e-7)
  assert_true(cert.max_complementarity_error < 1.0e-7)
}

///|
test "19. parser error handling for malformed .lp and .mps inputs" {
  let lp_err = try {
    let _ = parse_lp("Subject To\n c1: x <= 5\nEnd")
    false
  } catch {
    ParseError(_) => true
  }
  assert_true(lp_err)

  let mps_err = try {
    let _ = parse_mps("NAME BAD\nENDATA")
    false
  } catch {
    ParseError(_) => true
  }
  assert_true(mps_err)
}

///|
test "20. multi-period supply chain production & inventory LP" {
  // 3 periods: demand = [10, 20, 15], production cost = [4, 7, 5], holding cost = 1/unit/period
  // Production capacity per period = 20, initial inventory = 0
  // Optimal strategy: produce 20 in p0 (cost 80, hold 10 -> +10), produce 10 in p1 (cost 70, inv 0), produce 15 in p2 (cost 75) -> total = 235!
  let m = Model::new(Minimize)
  let p0 = m.add_continuous_var("p0", ub=20.0, obj=4.0)
  let p1 = m.add_continuous_var("p1", ub=20.0, obj=7.0)
  let p2 = m.add_continuous_var("p2", ub=20.0, obj=5.0)
  let inv0 = m.add_continuous_var("inv0", obj=1.0)
  let inv1 = m.add_continuous_var("inv1", obj=1.0)
  let inv2 = m.add_continuous_var("inv2", obj=1.0)

  m.add_constraint(
    "bal0",
    [{ var_id: p0, coeff: 1.0, }, { var_id: inv0, coeff: -1.0, }],
    Eq,
    10.0,
  )
  m.add_constraint(
    "bal1",
    [
      { var_id: inv0, coeff: 1.0, },
      { var_id: p1, coeff: 1.0, },
      { var_id: inv1, coeff: -1.0, },
    ],
    Eq,
    20.0,
  )
  m.add_constraint(
    "bal2",
    [
      { var_id: inv1, coeff: 1.0, },
      { var_id: p2, coeff: 1.0, },
      { var_id: inv2, coeff: -1.0, },
    ],
    Eq,
    15.0,
  )

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 235.0, 1.0e-6)
}

///|
test "21. Netlib MPS export -> parse_mps -> solve round-trip (LP & MILP)" {
  // Build a mixed integer + continuous model, export to MPS, re-parse and verify identical optimum
  let m = Model::new(Minimize)
  let x = m.add_int_var("X1", lb=1.0, ub=6.0, obj=3.0)
  let y = m.add_binary_var("B1", obj=5.0)
  let z = m.add_continuous_var("C1", lb=0.0, ub=10.0, obj=2.0)
  m.add_constraint(
    "DEM",
    [
      { var_id: x, coeff: 1.0, },
      { var_id: y, coeff: 3.0, },
      { var_id: z, coeff: 1.0, },
    ],
    Ge,
    6.5,
  )
  let sol1 = m.solve()
  assert_eq(sol1.status, Optimal)

  let mps = m.to_mps_string(name="ROUNDTRIP")
  let parsed = parse_mps(mps)
  let sol2 = parsed.solve()
  assert_eq(sol2.status, Optimal)
  assert_true(parsed.verify_solution(sol2))
  assert_approx(sol2.objective, sol1.objective, 1.0e-6)
}

///|
test "22. post-optimality constraint sensitivity & exact shadow price perturbation" {
  // Maximize 3*x + 2*y s.t. x + y <= 4, 2*x + y <= 5, x <= 10 (non-binding)
  // Optimum: x=1, y=3, obj=9; shadow prices: labor=1, machine=1, cap=0
  let m = Model::new(Maximize)
  let x = m.add_continuous_var("x", obj=3.0)
  let y = m.add_continuous_var("y", obj=2.0)
  m.add_constraint(
    "labor",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    4.0,
  )
  m.add_constraint(
    "machine",
    [{ var_id: x, coeff: 2.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    5.0,
  )
  m.add_constraint("extra_cap", [{ var_id: x, coeff: 1.0, }], Le, 10.0)

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  let cs = m.analyze_constraints(sol)
  assert_eq(cs.length(), 3)
  assert_true(cs[0].binding)
  assert_true(cs[1].binding)
  assert_true(!cs[2].binding)
  assert_approx(cs[2].slack, 9.0, 1.0e-6)
  assert_approx(cs[2].shadow_price, 0.0, 1.0e-6)

  // Perturbing labor RHS by +0.5 must increase objective by shadow_price * 0.5 = +0.5
  match m.rhs_perturbation_delta(sol, 0, 0.5) {
    Some(d) => assert_approx(d, 0.5, 1.0e-6)
    None => fail("Expected optimal perturbed LP")
  }
}

///|
test "23. variable reduced-cost and bound-activity sensitivity analysis" {
  let m = Model::new(Minimize)
  let x1 = m.add_continuous_var("x1", obj=0.6)
  let x2 = m.add_continuous_var("x2", obj=0.35)
  m.add_constraint(
    "vitA",
    [{ var_id: x1, coeff: 5.0, }, { var_id: x2, coeff: 7.0, }],
    Ge,
    8.0,
  )
  m.add_constraint(
    "vitB",
    [{ var_id: x1, coeff: 4.0, }, { var_id: x2, coeff: 2.0, }],
    Ge,
    15.0,
  )
  let sol = m.solve()
  let vs = m.analyze_variables(sol)
  assert_eq(vs.length(), 2)
  assert_true(!vs[0].at_lower_bound)
  assert_true(vs[1].at_lower_bound)
  assert_approx(vs[0].reduced_cost, 0.0, 1.0e-6)
  assert_true(vs[1].reduced_cost > 0.0)
}

///|
test "24. 10-item 0-1 Knapsack MILP multi-level Branch-and-Bound stress test" {
  // Exact DP optimum = 77
  let vals = [15.0, 22.0, 14.0, 26.0, 32.0, 9.0, 16.0, 8.0, 25.0, 19.0]
  let wts = [4.0, 7.0, 5.0, 8.0, 10.0, 3.0, 5.0, 2.0, 8.0, 6.0]
  let m = Model::new(Maximize)
  let terms : Array[Term] = []
  for i = 0; i < 10; i = i + 1 {
    let v = m.add_binary_var("k_" + i.to_string(), obj=vals[i])
    terms.push({ var_id: v, coeff: wts[i], })
  }
  m.add_constraint("cap", terms, Le, 23.0)

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 77.0, 1.0e-6)
}

///|
test "25. LP blending & linear expression builder (LinearExpr) end-to-end" {
  // Blend two alloys x, y to meet total weight = 100, min grade >= 45%, max impurity <= 12%
  // Alloy A: cost 4/kg, grade 60%, impurity 10%
  // Alloy B: cost 2.5/kg, grade 30%, impurity 15%
  // From grade: 0.6x + 0.3(100-x) >= 45 => 0.3x >= 15 => x >= 50
  // From impurity: 0.1x + 0.15(100-x) <= 12 => -0.05x <= -3 => x >= 60
  // Optimal: x = 60, y = 40 -> cost = 60*4 + 40*2.5 = 340
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("alloy_A")
  let y = m.add_continuous_var("alloy_B")
  let obj_expr = LinearExpr::from_var(x, coeff=4.0).add_term(y, 2.5)
  m.set_objective(obj_expr)
  m.add_expr_constraint(
    "total_mass",
    LinearExpr::from_var(x).add_term(y, 1.0),
    Eq,
    100.0,
  )
  m.add_expr_constraint(
    "min_grade",
    LinearExpr::from_var(x, coeff=0.6).add_term(y, 0.3),
    Ge,
    45.0,
  )
  m.add_expr_constraint(
    "max_impurity",
    LinearExpr::from_var(x, coeff=0.1).add_term(y, 0.15),
    Le,
    12.0,
  )

  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.primal[x.index], 60.0, 1.0e-6)
  assert_approx(sol.primal[y.index], 40.0, 1.0e-6)
  assert_approx(sol.objective, 340.0, 1.0e-6)
}
