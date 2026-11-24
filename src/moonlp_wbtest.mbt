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

///|
test "26. Dual Simplex direct solve on >= Diet LP without Phase I" {
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

  let sol_dual = solve_lp_dual(m)
  assert_eq(sol_dual.status, Optimal)
  assert_true(m.verify_solution(sol_dual))
  assert_approx(sol_dual.objective, 2.25, 1.0e-7)
  let cert = m.certify_lp_solution(sol_dual)
  assert_true(cert.certified_optimal)
  assert_true(cert.duality_gap < 1.0e-7)
}

///|
test "27. Dual Simplex primal infeasibility detection" {
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("x", obj=1.0)
  let y = m.add_continuous_var("y", obj=2.0)
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
    6.0,
  )
  let sol = solve_lp_dual(m)
  assert_eq(sol.status, Infeasible)
}

///|
test "28. Gomory Mixed-Integer (GMI) cuts + Dual Simplex solve classic IP at root" {
  // Maximize 5*x + 8*y s.t. x + y <= 6, 5*x + 9*y <= 45, x,y >= 0 integer
  // Pure B&B takes 9 nodes; GMI cuts + Dual Simplex solve it right at the root (1 node)!
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

  let (sol, rep) = solve_milp_with_cuts(m)
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 40.0, 1.0e-6)
  assert_true(rep.cuts_added > 0)
  assert_true(rep.solved_at_root)
  assert_eq(sol.nodes_explored, 1)
}

///|
test "29. GMI cuts close large LP-to-IP gap (63.0 -> 55.0) at root in 2 cuts" {
  // Maximize 7*x1 + 9*x2 s.t. -x1 + 3*x2 <= 6, 7*x1 + x2 <= 35, x1,x2 >= 0 integer
  let m = Model::new(Maximize)
  let x1 = m.add_int_var("x1", obj=7.0)
  let x2 = m.add_int_var("x2", obj=9.0)
  m.add_constraint(
    "c1",
    [{ var_id: x1, coeff: -1.0, }, { var_id: x2, coeff: 3.0, }],
    Le,
    6.0,
  )
  m.add_constraint(
    "c2",
    [{ var_id: x1, coeff: 7.0, }, { var_id: x2, coeff: 1.0, }],
    Le,
    35.0,
  )

  let (sol, rep) = solve_milp_with_cuts(m)
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(rep.root_lp_obj, 63.0, 1.0e-6)
  assert_approx(sol.objective, 55.0, 1.0e-6)
  assert_approx(sol.primal[x1.index], 4.0, 1.0e-6)
  assert_approx(sol.primal[x2.index], 3.0, 1.0e-6)
  assert_true(rep.solved_at_root)
  assert_eq(sol.nodes_explored, 1)
}

///|
test "30. GMI cut solves true Mixed Integer-Continuous MILP (x int, y cont) in 1 cut" {
  // Maximize 5*x + 4*y s.t. x + y <= 5, 10*x + 6*y <= 45, x integer, y continuous
  // LP opt = 23.75 (x=3.75, y=1.25); Mixed IP opt = 23.333333 (x=4, y=5/6)
  let m = Model::new(Maximize)
  let x = m.add_int_var("x", obj=5.0)
  let y = m.add_continuous_var("y", obj=4.0)
  m.add_constraint(
    "c1",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    5.0,
  )
  m.add_constraint(
    "c2",
    [{ var_id: x, coeff: 10.0, }, { var_id: y, coeff: 6.0, }],
    Le,
    45.0,
  )

  let (sol, rep) = solve_milp_with_cuts(m)
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(rep.root_lp_obj, 23.75, 1.0e-6)
  assert_approx(sol.objective, 70.0 / 3.0, 1.0e-5)
  assert_approx(sol.primal[x.index], 4.0, 1.0e-6)
  assert_approx(sol.primal[y.index], 5.0 / 6.0, 1.0e-5)
  assert_true(rep.solved_at_root)
  assert_eq(rep.cuts_added, 1)
}

///|
test "31. 3x3 Transportation LP (9 variables) + KKT strong duality" {
  let supply = [20.0, 30.0, 25.0]
  let demand = [15.0, 35.0, 25.0]
  let costs = [[8.0, 6.0, 10.0], [9.0, 12.0, 13.0], [14.0, 9.0, 16.0]]
  let m = Model::new(Minimize)
  let x : Array[Array[VarId]] = []
  for i = 0; i < 3; i = i + 1 {
    let row : Array[VarId] = []
    for j = 0; j < 3; j = j + 1 {
      row.push(
        m.add_continuous_var(
          "x_" + i.to_string() + "_" + j.to_string(),
          obj=costs[i][j],
        ),
      )
    }
    x.push(row)
  }
  for i = 0; i < 3; i = i + 1 {
    m.add_constraint(
      "sup_" + i.to_string(),
      [
        { var_id: x[i][0], coeff: 1.0, },
        { var_id: x[i][1], coeff: 1.0, },
        { var_id: x[i][2], coeff: 1.0, },
      ],
      Le,
      supply[i],
    )
  }
  for j = 0; j < 3; j = j + 1 {
    m.add_constraint(
      "dem_" + j.to_string(),
      [
        { var_id: x[0][j], coeff: 1.0, },
        { var_id: x[1][j], coeff: 1.0, },
        { var_id: x[2][j], coeff: 1.0, },
      ],
      Ge,
      demand[j],
    )
  }
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol))
  assert_approx(sol.objective, 715.0, 1.0e-6)
  let cert = m.certify_lp_solution(sol)
  assert_true(cert.certified_optimal)
}

///|
test "32. 4-period capacitated lot-sizing MILP (setup binary + continuous flow)" {
  let dem = [10.0, 15.0, 20.0, 12.0]
  let ucost = [3.0, 4.0, 3.0, 5.0]
  let cap = 30.0
  let m = Model::new(Minimize)
  let y : Array[VarId] = []
  let p : Array[VarId] = []
  let inv : Array[VarId] = []
  for t = 0; t < 4; t = t + 1 {
    y.push(m.add_binary_var("y" + t.to_string(), obj=50.0))
    p.push(m.add_continuous_var("p" + t.to_string(), ub=cap, obj=ucost[t]))
    inv.push(m.add_continuous_var("I" + t.to_string(), obj=1.0))
  }
  for t = 0; t < 4; t = t + 1 {
    let terms : Array[Term] = [
      { var_id: p[t], coeff: 1.0, },
      { var_id: inv[t], coeff: -1.0, },
    ]
    if t > 0 {
      terms.push({ var_id: inv[t - 1], coeff: 1.0, })
    }
    m.add_constraint("bal" + t.to_string(), terms, Eq, dem[t])
    m.add_constraint(
      "cap" + t.to_string(),
      [{ var_id: p[t], coeff: 1.0, }, { var_id: y[t], coeff: -cap, }],
      Le,
      0.0,
    )
  }
  let (sol_cut, _) = solve_milp_with_cuts(m)
  assert_eq(sol_cut.status, Optimal)
  assert_true(m.verify_solution(sol_cut))
  assert_approx(sol_cut.objective, 302.0, 1.0e-5)
}

///|
test "33. authentic Netlib AFIRO benchmark (28 rows, 32 cols, opt = -464.75314286)" {
  let m = parse_mps(netlib_afiro_mps())
  assert_eq(m.num_vars(), 32)
  assert_eq(m.num_constraints(), 27)
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol, tol=1.0e-5))
  assert_approx(sol.objective, -464.75314285714285, 1.0e-5)
  let cert = m.certify_lp_solution(sol, tol=1.0e-5)
  assert_true(cert.certified_optimal)
}

///|
test "34. authentic HiGHS AVGAS benchmark (11 rows, 8 cols, opt = -7.75)" {
  let m = parse_mps(highs_avgas_mps())
  assert_eq(m.num_vars(), 8)
  assert_eq(m.num_constraints(), 10)
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol, tol=1.0e-5))
  assert_approx(sol.objective, -7.75, 1.0e-6)
}

///|
test "35. Dual Simplex Warm-Start B&B beats Cold-Start B&B on pivot count" {
  let vals = [15.0, 22.0, 14.0, 26.0, 32.0, 9.0, 16.0, 8.0, 25.0, 19.0]
  let wts = [4.0, 7.0, 5.0, 8.0, 10.0, 3.0, 5.0, 2.0, 8.0, 6.0]
  let m = Model::new(Maximize)
  let terms : Array[Term] = []
  for i = 0; i < 10; i = i + 1 {
    let v = m.add_binary_var("k_" + i.to_string(), obj=vals[i])
    terms.push({ var_id: v, coeff: wts[i], })
  }
  m.add_constraint("cap", terms, Le, 23.0)

  let cold = solve_milp(m)
  let warm = solve_milp_warm(m)
  assert_eq(cold.status, Optimal)
  assert_eq(warm.status, Optimal)
  assert_true(m.verify_solution(warm))
  assert_approx(warm.objective, 77.0, 1.0e-6)
  assert_approx(warm.objective, cold.objective, 1.0e-6)
  assert_true(warm.iterations < cold.iterations)
}

///|
test "36. 5-City Traveling Salesman Problem (TSP-MTZ MILP, opt = 80) via Warm-Start B&B" {
  let m = build_tsp5_mtz_model()
  let cold = solve_milp(m)
  let warm = solve_milp_warm(m)
  assert_eq(cold.status, Optimal)
  assert_eq(warm.status, Optimal)
  assert_true(m.verify_solution(warm, tol=1.0e-5))
  assert_approx(warm.objective, 80.0, 1.0e-5)
  assert_true(warm.iterations < cold.iterations)
}

///|
test "37. authentic Netlib ADLITTLE benchmark (56 rows, 97 cols, 480 nonzeros)" {
  let m = parse_mps(netlib_adlittle_mps())
  assert_eq(m.num_vars(), 97)
  assert_eq(m.num_constraints(), 56)
  let sol = m.solve()
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol, tol=1.0e-4))
  assert_approx(sol.objective, 225494.96316238, 1.0e-2)
  let cert = m.certify_lp_solution(sol, tol=1.0e-3)
  assert_true(cert.certified_optimal)
}

///|
test "38. Ruiz L-infinity equilibration scaling slashes dynamic range & preserves exact duals" {
  // Ill-conditioned LP with 1e6 coefficient ratio
  let m = Model::new(Maximize)
  let x = m.add_continuous_var("x", obj=3000.0)
  let y = m.add_continuous_var("y", obj=0.002)
  m.add_constraint(
    "r1",
    [{ var_id: x, coeff: 1000.0, }, { var_id: y, coeff: 0.001, }],
    Le,
    4.0,
  )
  m.add_constraint(
    "r2",
    [{ var_id: x, coeff: 2000.0, }, { var_id: y, coeff: 0.001, }],
    Le,
    5.0,
  )
  let (sol_s, rep) = m.solve_scaled(passes=6)
  assert_eq(sol_s.status, Optimal)
  assert_true(rep.original_ratio >= 1.0e6)
  assert_true(rep.scaled_ratio < 5.0)
  assert_true(m.verify_solution(sol_s, tol=1.0e-5))
  assert_approx(sol_s.objective, 9.0, 1.0e-6)
  let cert = m.certify_lp_solution(sol_s, tol=1.0e-5)
  assert_true(cert.certified_optimal)

  // Also test Ruiz scaling on real Netlib ADLITTLE (56x97)
  let ad = parse_mps(netlib_adlittle_mps())
  let (ad_sol, ad_rep) = ad.solve_scaled(passes=5)
  assert_eq(ad_sol.status, Optimal)
  assert_true(ad_rep.scaled_ratio < ad_rep.original_ratio)
  assert_approx(ad_sol.objective, 225494.96316238, 1.0e-2)
}

///|
test "39. authentic HiGHS SMALL_MIP with MPS RANGES section (opt = 61.5/19)" {
  let m = parse_mps(highs_small_mip_mps())
  // 5 base rows + 2 ranged companion rows = 7 constraints
  assert_eq(m.num_vars(), 8)
  assert_eq(m.num_constraints(), 7)
  let (sol, _) = solve_milp_branch_and_cut(m)
  assert_eq(sol.status, Optimal)
  assert_true(m.verify_solution(sol, tol=1.0e-5))
  assert_approx(sol.objective, 61.5 / 19.0, 1.0e-5)
}

///|
test "40. authentic MIPLIB 3.0 FLUGPL airline scheduling (LP = 1167185.73 + Branch-and-Cut)" {
  let m = parse_mps(miplib_flugpl_mps())
  assert_eq(m.num_vars(), 18)
  assert_eq(m.num_constraints(), 18)
  let lp_sol = solve_lp(m)
  assert_eq(lp_sol.status, Optimal)
  assert_approx(lp_sol.objective, 1167185.73, 0.05)

  // Unified Branch-and-Cut (GMI cuts + Hybrid DFS-Plunge + Dual Warm-Start + Rounding Heuristic)
  let cfg : SolverConfig = { ..SolverConfig::default(), max_nodes: 400, }
  let (mip_sol, stats) = solve_milp_branch_and_cut(m, max_root_cuts=6, config=cfg)
  assert_true(stats.cuts_added > 0)
  // Root GMI cuts must strictly tighten LP bound above 1167185.73
  assert_true(stats.post_cut_lp_obj > stats.root_lp_obj + 1.0)
  if mip_sol.status == Optimal || mip_sol.status == NodeLimit {
    if mip_sol.primal.length() == 18 {
      assert_true(mip_sol.objective >= 1201500.0 - 1.0)
    }
  }
}

///|
test "41. unified Branch-and-Cut exactness across Knapsack, TSP-MTZ, and Mixed IP" {
  // 1) 5-city TSP-MTZ (exact = 80.0)
  let tsp = build_tsp5_mtz_model()
  let (tsp_sol, _) = solve_milp_branch_and_cut(tsp)
  assert_eq(tsp_sol.status, Optimal)
  assert_true(tsp.verify_solution(tsp_sol, tol=1.0e-5))
  assert_approx(tsp_sol.objective, 80.0, 1.0e-5)

  // 2) Mixed Integer-Continuous MILP (exact = 70/3)
  let mx = Model::new(Maximize)
  let x = mx.add_int_var("x", obj=5.0)
  let y = mx.add_continuous_var("y", obj=4.0)
  mx.add_constraint(
    "c1",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    5.0,
  )
  mx.add_constraint(
    "c2",
    [{ var_id: x, coeff: 10.0, }, { var_id: y, coeff: 6.0, }],
    Le,
    45.0,
  )
  let (mx_sol, mx_stats) = solve_milp_branch_and_cut(mx)
  assert_eq(mx_sol.status, Optimal)
  assert_true(mx.verify_solution(mx_sol))
  assert_approx(mx_sol.objective, 70.0 / 3.0, 1.0e-5)
  assert_eq(mx_stats.nodes_explored, 1)
}

///|
test "42. [P0-1] MILP max_iters interruption never misclassified as Infeasible" {
  let m = Model::new(Minimize)
  let x = m.add_binary_var("x", obj=1.0)
  m.add_constraint("fix_x", [{ var_id: x, coeff: 1.0, }], Eq, 1.0)

  let tight_cfg : SolverConfig = { ..SolverConfig::default(), max_iters: 1, }
  let lp_sol = solve_lp(m, config=tight_cfg)
  assert_eq(lp_sol.status, MaxIterations)

  let s1 = solve_milp(m, config=tight_cfg)
  assert_eq(s1.status, MaxIterations)
  let s2 = solve_milp_warm(m, config=tight_cfg)
  assert_eq(s2.status, MaxIterations)
  let (s3, _) = solve_milp_with_cuts(m, config=tight_cfg)
  assert_eq(s3.status, MaxIterations)
  let (s4, _) = solve_milp_branch_and_cut(m, config=tight_cfg)
  assert_eq(s4.status, MaxIterations)

  // Default budget solves cleanly to Optimal (x = 1)
  let ok = m.solve()
  assert_eq(ok.status, Optimal)
  assert_approx(ok.objective, 1.0, 1.0e-6)
  assert_approx(ok.primal[0], 1.0, 1.0e-6)
}

///|
test "43. [P0-2] Unbounded LP relaxation disambiguated between Infeasible and Unbounded MILP" {
  // Case A: LP relaxation unbounded (x = 0.5, y -> +inf), MILP is Infeasible!
  let m_inf = Model::new(Maximize)
  let x = m_inf.add_binary_var("x", obj=0.0)
  let _y = m_inf.add_continuous_var("y", lb=0.0, ub=1.0e300, obj=1.0)
  m_inf.add_constraint("half_x", [{ var_id: x, coeff: 1.0, }], Eq, 0.5)

  assert_eq(solve_lp(m_inf).status, Unbounded)
  assert_eq(solve_milp(m_inf).status, Infeasible)
  assert_eq(solve_milp_warm(m_inf).status, Infeasible)
  let (s3, _) = solve_milp_with_cuts(m_inf)
  assert_eq(s3.status, Infeasible)
  let (s4, _) = solve_milp_branch_and_cut(m_inf)
  assert_eq(s4.status, Infeasible)

  // Case B: LP relaxation unbounded AND integer feasible (x = 0, y -> +inf) => Unbounded
  let m_unb = Model::new(Maximize)
  let xb = m_unb.add_binary_var("x", obj=0.0)
  let _yb = m_unb.add_continuous_var("y", lb=0.0, ub=1.0e300, obj=1.0)
  m_unb.add_constraint("le_x", [{ var_id: xb, coeff: 1.0, }], Le, 1.0)

  assert_eq(solve_milp(m_unb).status, Unbounded)
  assert_eq(solve_milp_warm(m_unb).status, Unbounded)
  let (u3, _) = solve_milp_with_cuts(m_unb)
  assert_eq(u3.status, Unbounded)
  let (u4, _) = solve_milp_branch_and_cut(m_unb)
  assert_eq(u4.status, Unbounded)
}

///|
test "44. [P0-3] Sub-tolerance integer rounding never accepts infeasible or drops optimal node" {
  // max x, x in {0, 1}, 1e8 * x <= 99_999_950
  // LP relaxation gives x = 0.9999995 (within 1e-6 of 1), but x = 1 violates constraint by 50!
  // True MILP optimum is x = 0, obj = 0.
  let m = Model::new(Maximize)
  let x = m.add_binary_var("x", obj=1.0)
  m.add_constraint("big_c", [{ var_id: x, coeff: 1.0e8, }], Le, 99999950.0)

  let s1 = solve_milp(m)
  assert_eq(s1.status, Optimal)
  assert_true(m.verify_solution(s1))
  assert_approx(s1.primal[0], 0.0, 1.0e-9)
  assert_approx(s1.objective, 0.0, 1.0e-9)

  let s2 = solve_milp_warm(m)
  assert_eq(s2.status, Optimal)
  assert_true(m.verify_solution(s2))
  assert_approx(s2.primal[0], 0.0, 1.0e-9)

  let (s3, _) = solve_milp_with_cuts(m)
  assert_eq(s3.status, Optimal)
  assert_true(m.verify_solution(s3))
  assert_approx(s3.primal[0], 0.0, 1.0e-9)

  let (s4, _) = solve_milp_branch_and_cut(m)
  assert_eq(s4.status, Optimal)
  assert_true(m.verify_solution(s4))
  assert_approx(s4.primal[0], 0.0, 1.0e-9)
}

///|
test "45. [P0-4] Strict KKT optimality certificate rejects suboptimal points & NaN/Inf" {
  // Min x, 0 <= x <= 1
  let m = Model::new(Minimize)
  let _x = m.add_continuous_var("x", lb=0.0, ub=1.0, obj=1.0)

  let subopt : Solution = {
    status: Optimal,
    objective: 1.0,
    primal: [1.0],
    dual: [],
    reduced_costs: [1.0],
    iterations: 0,
    nodes_explored: 1,
    mip_gap: 0.0,
  }
  let cert_bad = m.certify_lp_solution(subopt)
  assert_false(cert_bad.certified_optimal)
  assert_false(cert_bad.dual_feasible)
  assert_true(cert_bad.duality_gap >= 0.99)

  let opt_sol = m.solve()
  let cert_ok = m.certify_lp_solution(opt_sol)
  assert_true(cert_ok.certified_optimal)
  assert_true(cert_ok.primal_feasible)
  assert_true(cert_ok.dual_feasible)
  assert_true(cert_ok.duality_gap < 1.0e-9)

  // Max x, 0 <= x <= 1: reject x = 0, certify x = 1
  let m2 = Model::new(Maximize)
  let _x2 = m2.add_continuous_var("x", lb=0.0, ub=1.0, obj=1.0)
  let subopt2 : Solution = { ..subopt, objective: 0.0, primal: [0.0], }
  assert_false(m2.certify_lp_solution(subopt2).certified_optimal)
  assert_true(m2.certify_lp_solution(m2.solve()).certified_optimal)

  // Reject NaN and +Inf in verify_solution and certify_lp_solution
  let nan = 0.0 / 0.0
  let pos_inf = 1.0 / 0.0
  let nan_sol : Solution = { ..subopt, objective: nan, primal: [nan], }
  let inf_sol : Solution = { ..subopt, objective: pos_inf, primal: [pos_inf], }
  assert_false(m.verify_solution(nan_sol))
  assert_false(m.verify_solution(inf_sol))
  assert_false(m.certify_lp_solution(nan_sol).certified_optimal)
  assert_false(m.certify_lp_solution(inf_sol).certified_optimal)
}

///|
test "46. [P1-1] Default solve_lp handles extreme coefficient scales & tight ratio tests" {
  // Case A: max x, 1e-10 * x <= 1 => opt x = 1e10, obj = 1e10 (must NOT report Unbounded)
  let m1 = Model::new(Maximize)
  let x1 = m1.add_continuous_var("x", lb=0.0, obj=1.0)
  m1.add_constraint("tiny", [{ var_id: x1, coeff: 1.0e-10, }], Le, 1.0)
  let s1 = m1.solve()
  assert_eq(s1.status, Optimal)
  assert_true(m1.verify_solution(s1, tol=1.0e-4))
  assert_approx(s1.objective, 1.0e10, 1.0e4)
  assert_approx(s1.primal[0], 1.0e10, 1.0e4)

  // Case B: two-constraint tight ratio test where smaller-pivot row is tighter by 0.5
  let m2 = Model::new(Maximize)
  let x2 = m2.add_continuous_var("x", lb=0.0, obj=1.0)
  m2.add_constraint("tight_small", [{ var_id: x2, coeff: 1.0e-5, }], Le, 0.999995)
  m2.add_constraint("loose_big", [{ var_id: x2, coeff: 1.0, }], Le, 100000.0)
  let s2 = m2.solve()
  assert_eq(s2.status, Optimal)
  assert_true(m2.verify_solution(s2, tol=1.0e-5))
  assert_approx(s2.primal[0], 99999.5, 1.0e-4)
  assert_approx(s2.objective, 99999.5, 1.0e-4)
}

///|
test "47. [P1-2a] CPLEX .lp round-trip preserves positive and negative obj_constant" {
  let m_pos = Model::new(Maximize)
  let x = m_pos.add_continuous_var("x", lb=0.0, ub=4.0, obj=3.0)
  let y = m_pos.add_continuous_var("y", lb=0.0, ub=5.0, obj=2.0)
  m_pos.set_obj_constant(17.5)
  m_pos.add_constraint("c1", [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }], Le, 6.0)

  let lp_pos = m_pos.to_lp_string()
  let p_pos = parse_lp(lp_pos)
  assert_approx(p_pos.obj_constant, 17.5, 1.0e-9)
  let sol_orig = m_pos.solve()
  let sol_rt = p_pos.solve()
  assert_eq(sol_rt.status, Optimal)
  assert_approx(sol_rt.objective, sol_orig.objective, 1.0e-6)

  let m_neg = Model::new(Minimize)
  let a = m_neg.add_continuous_var("a", lb=1.0, ub=5.0, obj=2.0)
  m_neg.set_obj_constant(-12.25)
  m_neg.add_constraint("c1", [{ var_id: a, coeff: 1.0, }], Ge, 2.0)
  let p_neg = parse_lp(m_neg.to_lp_string())
  assert_approx(p_neg.obj_constant, -12.25, 1.0e-9)
  assert_approx(p_neg.solve().objective, m_neg.solve().objective, 1.0e-6)
}

///|
test "48. [P1-2b] Netlib MPS round-trip preserves OBJSENSE MAX, obj_constant & fixed binaries" {
  let m = Model::new(Maximize)
  let b1 = m.add_binary_var("B1", obj=5.0)
  let b2 = m.add_binary_var("B2", obj=3.0)
  let x = m.add_continuous_var("X1", lb=0.0, ub=4.0, obj=2.0)
  m.set_var_bounds(b1, 1.0, 1.0) // fixed binary = 1
  m.set_var_bounds(b2, 0.0, 0.0) // fixed binary = 0
  m.set_obj_constant(10.5)
  m.add_constraint(
    "CAP",
    [{ var_id: b1, coeff: 1.0, }, { var_id: b2, coeff: 1.0, }, { var_id: x, coeff: 1.0, }],
    Le,
    3.0,
  )

  let sol1 = m.solve()
  assert_eq(sol1.status, Optimal)
  assert_approx(sol1.objective, 19.5, 1.0e-6)

  let mps = m.to_mps_string(name="MAXFIX")
  let parsed = parse_mps(mps)
  assert_eq(parsed.sense, Maximize)
  assert_approx(parsed.obj_constant, 10.5, 1.0e-9)
  assert_approx(parsed.vars[0].lb, 1.0, 1.0e-9)
  assert_approx(parsed.vars[0].ub, 1.0, 1.0e-9)
  assert_approx(parsed.vars[1].lb, 0.0, 1.0e-9)
  assert_approx(parsed.vars[1].ub, 0.0, 1.0e-9)

  let sol2 = parsed.solve()
  assert_eq(sol2.status, Optimal)
  assert_true(parsed.verify_solution(sol2))
  assert_approx(sol2.objective, 19.5, 1.0e-6)
  assert_approx(sol2.primal[0], 1.0, 1.0e-6)
  assert_approx(sol2.primal[1], 0.0, 1.0e-6)
  assert_approx(sol2.primal[2], 2.0, 1.0e-6)
}

///|
test "49. [P1-2c] solve_presolved recovers full original-space duals & passes KKT certificate" {
  // Maximize 5*x + 4*y s.t.:
  //   c0 (singleton binding ub): 2*x <= 6      -> x <= 3  (dual = 2.5? Wait: coupled with c2!)
  //   c1 (singleton non-binding): y <= 10      -> dual = 0
  //   c2 (multi-var binding): x + y <= 5       -> x=3, y=2, obj=23; dual(c2)=4, dual(c0)=0.5
  //   c3 (singleton negative coeff binding lb): -2*y <= -2 -> y >= 1 (non-binding at y=2, dual=0)
  let m = Model::new(Maximize)
  let x = m.add_continuous_var("x", lb=0.0, obj=5.0)
  let y = m.add_continuous_var("y", lb=0.0, obj=4.0)
  m.add_constraint("c0_xub", [{ var_id: x, coeff: 2.0, }], Le, 6.0)
  m.add_constraint("c1_yub", [{ var_id: y, coeff: 1.0, }], Le, 10.0)
  m.add_constraint("c2_sum", [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }], Le, 5.0)
  m.add_constraint("c3_ylb", [{ var_id: y, coeff: -2.0, }], Le, -2.0)

  let psol = m.solve_presolved()
  assert_eq(psol.status, Optimal)
  assert_eq(psol.dual.length(), 4)
  assert_approx(psol.objective, 23.0, 1.0e-6)
  assert_approx(psol.dual[0], 0.5, 1.0e-6)
  assert_approx(psol.dual[1], 0.0, 1.0e-6)
  assert_approx(psol.dual[2], 4.0, 1.0e-6)
  assert_approx(psol.dual[3], 0.0, 1.0e-6)
  let cert = m.certify_lp_solution(psol)
  assert_true(cert.certified_optimal)

  // Minimize 3*x + 2*y s.t. 2*x >= 4 (binding singleton lb), y = 3 (binding singleton eq)
  // All constraints eliminated by presolve -> m_kept = 0!
  let m2 = Model::new(Minimize)
  let x2 = m2.add_continuous_var("x", lb=0.0, obj=3.0)
  let y2 = m2.add_continuous_var("y", lb=0.0, obj=2.0)
  m2.add_constraint("s_lb", [{ var_id: x2, coeff: 2.0, }], Ge, 4.0)
  m2.add_constraint("s_eq", [{ var_id: y2, coeff: 1.0, }], Eq, 3.0)
  let psol2 = m2.solve_presolved()
  assert_eq(psol2.status, Optimal)
  assert_eq(psol2.dual.length(), 2)
  assert_approx(psol2.objective, 12.0, 1.0e-6)
  assert_approx(psol2.dual[0], 1.5, 1.0e-6)
  assert_approx(psol2.dual[1], 2.0, 1.0e-6)
  assert_true(m2.certify_lp_solution(psol2).certified_optimal)
}

///|
test "50. [v0.7 P0-1] Degraded incumbent repair never prematurely fathoms fractional node" {
  // max y s.t. y <= 2e7 * x, x + y <= 10, x in {0, 1}, y >= 0
  // Root LP has x ~ 5e-7 (within int_tol=1e-6), repair gives feasible (0,0) obj=0,
  // but node LP bound is ~10, so solver MUST branch and find true optimum (1, 9) obj=9.
  let m = Model::new(Maximize)
  let x = m.add_binary_var("x", obj=0.0)
  let y = m.add_continuous_var("y", lb=0.0, obj=1.0)
  m.add_constraint(
    "link",
    [{ var_id: y, coeff: 1.0, }, { var_id: x, coeff: -2.0e7, }],
    Le,
    0.0,
  )
  m.add_constraint(
    "cap",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: 1.0, }],
    Le,
    10.0,
  )

  let s1 = solve_milp(m)
  assert_eq(s1.status, Optimal)
  assert_true(m.verify_solution(s1))
  assert_approx(s1.objective, 9.0, 1.0e-6)
  assert_approx(s1.primal[0], 1.0, 1.0e-6)
  assert_approx(s1.primal[1], 9.0, 1.0e-6)

  let s2 = solve_milp_warm(m)
  assert_eq(s2.status, Optimal)
  assert_true(m.verify_solution(s2))
  assert_approx(s2.objective, 9.0, 1.0e-6)

  let (s3, _) = solve_milp_with_cuts(m)
  assert_eq(s3.status, Optimal)
  assert_true(m.verify_solution(s3))
  assert_approx(s3.objective, 9.0, 1.0e-6)

  let (s4, _) = solve_milp_branch_and_cut(m)
  assert_eq(s4.status, Optimal)
  assert_true(m.verify_solution(s4))
  assert_approx(s4.objective, 9.0, 1.0e-6)
}

///|
test "51. [v0.7 P0-2] GMI cuts retain small nonzero coefficients & never cut off optimum" {
  // max x - 2.1e-9 * y s.t. x - 2e-9 * y <= 0.1, 0 <= y <= 4.6e8, x in {0, 1}
  // Exact optimum: x = 1, y = 4.5e8, obj = 0.055.
  let m = Model::new(Maximize)
  let x = m.add_binary_var("x", obj=1.0)
  let y = m.add_continuous_var("y", lb=0.0, ub=4.6e8, obj=-2.1e-9)
  m.add_constraint(
    "c1",
    [{ var_id: x, coeff: 1.0, }, { var_id: y, coeff: -2.0e-9, }],
    Le,
    0.1,
  )

  let s1 = solve_milp(m)
  assert_eq(s1.status, Optimal)
  assert_true(m.verify_solution(s1, tol=1.0e-5))
  assert_approx(s1.objective, 0.055, 1.0e-6)

  let s2 = solve_milp_warm(m)
  assert_eq(s2.status, Optimal)
  assert_true(m.verify_solution(s2, tol=1.0e-5))
  assert_approx(s2.objective, 0.055, 1.0e-6)

  let (s3, _) = solve_milp_with_cuts(m)
  assert_eq(s3.status, Optimal)
  assert_true(m.verify_solution(s3, tol=1.0e-5))
  assert_approx(s3.objective, 0.055, 1.0e-6)

  let (s4, _) = solve_milp_branch_and_cut(m)
  assert_eq(s4.status, Optimal)
  assert_true(m.verify_solution(s4, tol=1.0e-5))
  assert_approx(s4.objective, 0.055, 1.0e-6)
}

///|
test "52. [v0.7 P0-3] Relative MIP gap pruning does not prune tiny-scale 1e-6 optimum" {
  // max 1e-6 * x + 9e-7 * y s.t. 1.5*x + y <= 1.5, x, y in {0, 1}
  // Root LP = 1.2333e-6; x=0 branch gives incumbent 9e-7; x=1 branch gives true optimum 1e-6.
  let m = Model::new(Maximize)
  let x = m.add_binary_var("x", obj=1.0e-6)
  let y = m.add_binary_var("y", obj=9.0e-7)
  m.add_constraint(
    "c1",
    [{ var_id: x, coeff: 1.5, }, { var_id: y, coeff: 1.0, }],
    Le,
    1.5,
  )

  let s1 = solve_milp(m)
  assert_eq(s1.status, Optimal)
  assert_approx(s1.objective, 1.0e-6, 1.0e-12)
  assert_approx(s1.primal[0], 1.0, 1.0e-9)

  let s2 = solve_milp_warm(m)
  assert_eq(s2.status, Optimal)
  assert_approx(s2.objective, 1.0e-6, 1.0e-12)

  let (s3, _) = solve_milp_branch_and_cut(m, max_root_cuts=0)
  assert_eq(s3.status, Optimal)
  assert_approx(s3.objective, 1.0e-6, 1.0e-12)

  let (s4, _) = solve_milp_branch_and_cut(m)
  assert_eq(s4.status, Optimal)
  assert_approx(s4.objective, 1.0e-6, 1.0e-12)
}

///|
test "53. [v0.7 P0-4] KKT certificate primal/dual feasibility is invariant to huge obj_constant" {
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("x", lb=0.0, obj=0.0)
  m.set_obj_constant(1.0e9)
  m.add_constraint("zero_cap", [{ var_id: x, coeff: 1.0, }], Le, 0.0)

  let bad : Solution = {
    status: Optimal,
    objective: 1.0e9,
    primal: [1.0],
    dual: [0.0],
    reduced_costs: [0.0],
    iterations: 0,
    nodes_explored: 1,
    mip_gap: 0.0,
  }
  let cert_bad = m.certify_lp_solution(bad)
  assert_false(cert_bad.primal_feasible)
  assert_false(cert_bad.certified_optimal)
  assert_approx(cert_bad.max_primal_violation, 1.0, 1.0e-9)

  let opt = m.solve()
  let cert_ok = m.certify_lp_solution(opt)
  assert_true(cert_ok.primal_feasible)
  assert_true(cert_ok.certified_optimal)
}

///|
test "54. [v0.7 P1-1] Small-coeff boundedness across presolve & MILP relaxations + tight Phase I" {
  // A) max x, 1e-10 * x <= 1 via solve_presolved() must stay Optimal (not Unbounded!)
  let m1 = Model::new(Maximize)
  let x1 = m1.add_continuous_var("x", lb=0.0, obj=1.0)
  m1.add_constraint("tiny", [{ var_id: x1, coeff: 1.0e-10, }], Le, 1.0)
  let p1 = m1.solve_presolved()
  assert_eq(p1.status, Optimal)
  assert_approx(p1.objective, 1.0e10, 1.0e4)
  assert_true(m1.certify_lp_solution(p1, tol=1.0e-4).certified_optimal)

  // B) max x, 1e-14 * x <= 1 with an unrelated binary var across ALL 4 MILP entry points
  let m2 = Model::new(Maximize)
  let x2 = m2.add_continuous_var("x", lb=0.0, obj=1.0)
  let _b = m2.add_binary_var("b", obj=0.0)
  m2.add_constraint("ultra_tiny", [{ var_id: x2, coeff: 1.0e-14, }], Le, 1.0)

  let r1 = solve_milp(m2)
  assert_eq(r1.status, Optimal)
  assert_approx(r1.objective, 1.0e14, 1.0e8)
  let r2 = solve_milp_warm(m2)
  assert_eq(r2.status, Optimal)
  assert_approx(r2.objective, 1.0e14, 1.0e8)
  let (r3, _) = solve_milp_with_cuts(m2)
  assert_eq(r3.status, Optimal)
  assert_approx(r3.objective, 1.0e14, 1.0e8)
  let (r4, _) = solve_milp_branch_and_cut(m2)
  assert_eq(r4.status, Optimal)
  assert_approx(r4.objective, 1.0e14, 1.0e8)

  // C) Contradictory constraints differing by 5e-8 (> config.tol = 1e-9) must be Infeasible
  let m3 = Model::new(Minimize)
  let x3 = m3.add_continuous_var("x", lb=0.0, obj=1.0)
  m3.add_constraint("ub", [{ var_id: x3, coeff: 1.0, }], Le, 1.0)
  m3.add_constraint("lb", [{ var_id: x3, coeff: 1.0, }], Ge, 1.0 + 5.0e-8)
  assert_eq(solve_lp(m3).status, Infeasible)
}

///|
test "55. [v0.7 P1-2] Presolve assigns singleton shadow price to exact tightening row" {
  // max x s.t. c0: x <= 100_000_001, c1: x <= 100_000_000
  let m = Model::new(Maximize)
  let x = m.add_continuous_var("x", lb=0.0, obj=1.0)
  m.add_constraint("c0_loose", [{ var_id: x, coeff: 1.0, }], Le, 100000001.0)
  m.add_constraint("c1_tight", [{ var_id: x, coeff: 1.0, }], Le, 100000000.0)

  let direct = m.solve()
  assert_eq(direct.status, Optimal)
  assert_approx(direct.dual[0], 0.0, 1.0e-9)
  assert_approx(direct.dual[1], 1.0, 1.0e-9)

  let presolved = m.solve_presolved()
  assert_eq(presolved.status, Optimal)
  assert_eq(presolved.dual.length(), 2)
  assert_approx(presolved.dual[0], 0.0, 1.0e-9)
  assert_approx(presolved.dual[1], 1.0, 1.0e-9)
  assert_true(m.certify_lp_solution(presolved).certified_optimal)

  // Swapped wrong dual [1, 0] must be rejected by KKT certificate (violates complementarity by 1.0)
  let wrong_dual : Solution = { ..presolved, dual: [1.0, 0.0], }
  assert_false(m.certify_lp_solution(wrong_dual).certified_optimal)
}

///|
test "56. [v0.7 P1-3] CPLEX .lp parser rejects trailing tokens on constraint RHS & bounds" {
  let bad_rhs = try {
    let _ = parse_lp("Maximize\n x\nSubject To\n c1: x <= 1 + 1\nEnd")
    false
  } catch {
    ParseError(_) => true
  }
  assert_true(bad_rhs)

  let bad_bound = try {
    let _ = parse_lp(
      "Maximize\n x\nSubject To\n c1: x <= 5\nBounds\n x <= 2 + 3\nEnd",
    )
    false
  } catch {
    ParseError(_) => true
  }
  assert_true(bad_bound)
}

///|
test "57. [v0.8 P0-1] MILP candidate & heuristic verification strictly obeys config.tol" {
  // max x, x in {0, 1}, x <= 0.9999995 => violation of x=1 is 5e-7 > config.tol (1e-9).
  // Only feasible integer solution is x = 0.
  let m = Model::new(Maximize)
  let x = m.add_binary_var("x", obj=1.0)
  m.add_constraint("tight_ub", [{ var_id: x, coeff: 1.0, }], Le, 0.9999995)

  let s1 = solve_milp(m)
  assert_eq(s1.status, Optimal)
  assert_true(m.verify_solution(s1, tol=1.0e-9))
  assert_approx(s1.primal[0], 0.0, 1.0e-9)
  assert_approx(s1.objective, 0.0, 1.0e-9)

  let s2 = solve_milp_warm(m)
  assert_eq(s2.status, Optimal)
  assert_true(m.verify_solution(s2, tol=1.0e-9))
  assert_approx(s2.primal[0], 0.0, 1.0e-9)

  let (s3, _) = solve_milp_with_cuts(m)
  assert_eq(s3.status, Optimal)
  assert_true(m.verify_solution(s3, tol=1.0e-9))
  assert_approx(s3.primal[0], 0.0, 1.0e-9)

  let (s4, _) = solve_milp_branch_and_cut(m)
  assert_eq(s4.status, Optimal)
  assert_true(m.verify_solution(s4, tol=1.0e-9))
  assert_approx(s4.primal[0], 0.0, 1.0e-9)
}

///|
test "58. [v0.8 P0-2] Fixed-zero variable coefficient M=0 vs M=1e6 invariance across LP & MILP" {
  // a = 2^-29; max x - 1.0625*a*y
  // s.t. x - y + z + M*w <= 0.125, z - (1-a)*y = 0, 0 <= y <= 0.90625/a, w = 0, x in {0, 1}
  // Exact optimum for both M=0 and M=1e6: x = 1, y = 0.875/a, obj = 0.0703125
  let a = 1.0 / 536870912.0
  let y_ub = 0.90625 / a
  let cfg : SolverConfig = {
    ..SolverConfig::default(),
    max_iters: 100,
    max_nodes: 24,
  }
  let ms = [0.0, 1.0e6]

  for mi = 0; mi < ms.length(); mi = mi + 1 {
    let big_m = ms[mi]
    // Test both w fixed by constraint (w == 0) and w fixed by bounds (lb=ub=0)
    for mode = 0; mode < 2; mode = mode + 1 {
      let m = Model::new(Maximize)
      let x = m.add_binary_var("x", obj=1.0)
      let y = m.add_continuous_var("y", lb=0.0, ub=y_ub, obj=-1.0625 * a)
      let z = m.add_continuous_var("z", lb=0.0, obj=0.0)
      let w = if mode == 0 {
        m.add_continuous_var("w", lb=0.0, ub=0.0, obj=0.0)
      } else {
        m.add_continuous_var("w", lb=0.0, obj=0.0)
      }
      m.add_constraint(
        "c1",
        [
          { var_id: x, coeff: 1.0, },
          { var_id: y, coeff: -1.0, },
          { var_id: z, coeff: 1.0, },
          { var_id: w, coeff: big_m, },
        ],
        Le,
        0.125,
      )
      m.add_constraint(
        "c2",
        [{ var_id: z, coeff: 1.0, }, { var_id: y, coeff: -(1.0 - a), }],
        Eq,
        0.0,
      )
      if mode == 1 {
        m.add_constraint("fix_w", [{ var_id: w, coeff: 1.0, }], Eq, 0.0)
      }

      // Verify fixed x = 1 LP directly via public solve_lp()
      let m_x1 = m.clone_with_bounds(
        [1.0, 0.0, 0.0, 0.0],
        [1.0, y_ub, 1.0e300, if mode == 0 { 0.0 } else { 1.0e300 }],
      )
      let lp_x1 = solve_lp(m_x1, config=cfg)
      assert_eq(lp_x1.status, Optimal)
      assert_approx(lp_x1.objective, 0.0703125, 1.0e-6)

      let s1 = solve_milp(m, config=cfg)
      assert_eq(s1.status, Optimal)
      assert_true(m.verify_solution(s1, tol=1.0e-7))
      assert_approx(s1.objective, 0.0703125, 1.0e-6)

      let s2 = solve_milp_warm(m, config=cfg)
      assert_eq(s2.status, Optimal)
      assert_true(m.verify_solution(s2, tol=1.0e-7))
      assert_approx(s2.objective, 0.0703125, 1.0e-6)

      let (s3, _) = solve_milp_with_cuts(m, config=cfg)
      assert_eq(s3.status, Optimal)
      assert_true(m.verify_solution(s3, tol=1.0e-7))
      assert_approx(s3.objective, 0.0703125, 1.0e-6)

      let (s4, _) = solve_milp_branch_and_cut(m, config=cfg)
      assert_eq(s4.status, Optimal)
      assert_true(m.verify_solution(s4, tol=1.0e-7))
      assert_approx(s4.objective, 0.0703125, 1.0e-6)
    }
  }
}

///|
test "59. [v0.8 P0-3] KKT certificate rejects wrong-sign dual despite unrelated 1e9 obj coeff" {
  // min x + 1e9 * z s.t. 0 <= x <= 1, c0: x <= 1, c1: z = 0
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("x", lb=0.0, ub=1.0, obj=1.0)
  let z = m.add_continuous_var("z", lb=0.0, obj=1.0e9)
  m.add_constraint("x_ub", [{ var_id: x, coeff: 1.0, }], Le, 1.0)
  m.add_constraint("z_zero", [{ var_id: z, coeff: 1.0, }], Eq, 0.0)

  // Suboptimal (1, 0) with wrong-sign dual +1 on Le constraint x <= 1
  let bad_sign : Solution = {
    status: Optimal,
    objective: 1.0,
    primal: [1.0, 0.0],
    dual: [1.0, 1.0e9],
    reduced_costs: [0.0, 0.0],
    iterations: 0,
    nodes_explored: 1,
    mip_gap: 0.0,
  }
  let c1 = m.certify_lp_solution(bad_sign)
  assert_false(c1.dual_feasible)
  assert_false(c1.certified_optimal)
  assert_true(c1.max_dual_violation >= 0.99)

  // Suboptimal (1, 0) with valid-sign dual [0, 1e9] must also fail local stationarity
  let bad_stat : Solution = { ..bad_sign, dual: [0.0, 1.0e9], }
  let c2 = m.certify_lp_solution(bad_stat)
  assert_false(c2.dual_feasible)
  assert_false(c2.certified_optimal)

  // True optimal (0, 0) must be certified
  let opt = m.solve()
  assert_eq(opt.status, Optimal)
  assert_approx(opt.objective, 0.0, 1.0e-9)
  let c_ok = m.certify_lp_solution(opt)
  assert_true(c_ok.primal_feasible)
  assert_true(c_ok.dual_feasible)
  assert_true(c_ok.certified_optimal)
}
