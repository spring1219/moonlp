///|
/// Runnable demonstration of moonlp:
/// 1. Production Planning LP with Shadow Prices & KKT Optimality Certificate
/// 2. CPLEX .lp text parsing & MILP solving
/// 3. Live 4x4 Sudoku solved via 64-binary-variable Branch-and-Bound MILP

///|
fn main {
  println("=== moonlp: Pure MoonBit LP/MILP Optimization Solver Demo ===\n")

  // 1. Production Planning LP + KKT Strong Duality Certificate
  let m = @moonlp.Model::new(Maximize)
  let x = m.add_continuous_var("widget_A", obj=3.0)
  let y = m.add_continuous_var("widget_B", obj=2.0)
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
  let cert = m.certify_lp_solution(sol)
  println("[1] Production Mix LP (Maximize 3*A + 2*B)")
  println(
    "    Optimum = \{sol.objective} | A = \{sol.primal[0]}, B = \{sol.primal[1]}",
  )
  println(
    "    Shadow Prices = [labor: \{sol.dual[0]}, machine: \{sol.dual[1]}]",
  )
  println(
    "    KKT Strong Duality Certified = \{cert.certified_optimal} (gap = \{cert.duality_gap})\n",
  )

  // 2. Parse & Solve a MILP directly from CPLEX .lp text
  let lp_text =
    #|Maximize
    #|  profit: 5 x + 8 y
    #|Subject To
    #|  c1: x + y <= 6
    #|  c2: 5 x + 9 y <= 45
    #|Generals
    #|  x y
    #|End
  let parsed_res = try {
    let pm = @moonlp.parse_lp(lp_text)
    let psol = pm.solve()
    "Optimum = \{psol.objective} (x=\{psol.primal[0]}, y=\{psol.primal[1]}, B&B nodes=\{psol.nodes_explored})"
  } catch {
    @moonlp.ParseError(msg) => "Parse error: \{msg}"
  }
  println("[2] CPLEX .lp Format Parser + Integer Programming")
  println("    \{parsed_res}\n")

  // 3. Solve a 4x4 Sudoku live as a 64-binary-variable MILP
  let sm = @moonlp.Model::new(Minimize)
  let svars : Array[@moonlp.VarId] = []
  for r = 0; r < 4; r = r + 1 {
    for c = 0; c < 4; c = c + 1 {
      for d = 0; d < 4; d = d + 1 {
        svars.push(sm.add_binary_var("s_\{r}_\{c}_\{d}"))
      }
    }
  }
  let idx = fn(r : Int32, c : Int32, d : Int32) -> @moonlp.VarId {
    svars[(r * 4 + c) * 4 + d]
  }
  for r = 0; r < 4; r = r + 1 {
    for c = 0; c < 4; c = c + 1 {
      let terms : Array[@moonlp.Term] = []
      for d = 0; d < 4; d = d + 1 {
        terms.push({ var_id: idx(r, c, d), coeff: 1.0, })
      }
      sm.add_constraint("cell", terms, Eq, 1.0)
    }
  }
  for d = 0; d < 4; d = d + 1 {
    for r = 0; r < 4; r = r + 1 {
      let rt : Array[@moonlp.Term] = []
      for c = 0; c < 4; c = c + 1 {
        rt.push({ var_id: idx(r, c, d), coeff: 1.0, })
      }
      sm.add_constraint("row", rt, Eq, 1.0)
    }
    for c = 0; c < 4; c = c + 1 {
      let ct : Array[@moonlp.Term] = []
      for r = 0; r < 4; r = r + 1 {
        ct.push({ var_id: idx(r, c, d), coeff: 1.0, })
      }
      sm.add_constraint("col", ct, Eq, 1.0)
    }
    for br = 0; br < 2; br = br + 1 {
      for bc = 0; bc < 2; bc = bc + 1 {
        let bt : Array[@moonlp.Term] = []
        for dr = 0; dr < 2; dr = dr + 1 {
          for dc = 0; dc < 2; dc = dc + 1 {
            bt.push({ var_id: idx(br * 2 + dr, bc * 2 + dc, d), coeff: 1.0, })
          }
        }
        sm.add_constraint("box", bt, Eq, 1.0)
      }
    }
  }
  sm.set_var_bounds(idx(0, 0, 0), 1.0, 1.0)
  sm.set_var_bounds(idx(0, 3, 3), 1.0, 1.0)
  sm.set_var_bounds(idx(1, 2, 0), 1.0, 1.0)
  sm.set_var_bounds(idx(2, 1, 0), 1.0, 1.0)
  sm.set_var_bounds(idx(3, 0, 3), 1.0, 1.0)
  sm.set_var_bounds(idx(3, 3, 0), 1.0, 1.0)

  let ssol = sm.solve()
  println(
    "[3] 4x4 Sudoku Solved as 64-Binary-Variable MILP (pivots=\{ssol.iterations}):",
  )
  for r = 0; r < 4; r = r + 1 {
    let line = StringBuilder()
    line.write_string("    ")
    for c = 0; c < 4; c = c + 1 {
      let mut digit = 0
      for d = 0; d < 4; d = d + 1 {
        if ssol.primal[idx(r, c, d).index] > 0.5 {
          digit = d + 1
        }
      }
      line.write_string(digit.to_string())
      line.write_string(" ")
      if c == 1 {
        line.write_string("| ")
      }
    }
    println(line.to_string())
    if r == 1 {
      println("    ----+----")
    }
  }
}
