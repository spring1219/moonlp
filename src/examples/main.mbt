///|
/// Runnable demonstration of moonlp:
/// 1. Production Planning LP with Shadow Prices & KKT Optimality Certificate
/// 2. CPLEX .lp text parsing & MILP solving
/// 3. Live 4x4 Sudoku solved via 64-binary-variable Branch-and-Bound MILP

///|
fn main {
  println(
    "=== moonlp v1.1.0: Pure MoonBit SOCP / MISOCP / QCQP / MIQP / QP / LP / MILP Solver Demo ===\n",
  )

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

  // 4. Gomory Mixed-Integer (GMI) Cuts + Dual Simplex root-node solve
  let gm = @moonlp.Model::new(Maximize)
  let gx1 = gm.add_int_var("x1", obj=7.0)
  let gx2 = gm.add_int_var("x2", obj=9.0)
  gm.add_constraint(
    "c1",
    [{ var_id: gx1, coeff: -1.0, }, { var_id: gx2, coeff: 3.0, }],
    Le,
    6.0,
  )
  gm.add_constraint(
    "c2",
    [{ var_id: gx1, coeff: 7.0, }, { var_id: gx2, coeff: 1.0, }],
    Le,
    35.0,
  )
  let (gsol, grep) = @moonlp.solve_milp_with_cuts(gm)
  println("\n[4] Gomory Mixed-Integer (GMI) Cuts + Dual Simplex (Max 7*x1 + 9*x2)")
  println(
    "    Root LP = \{grep.root_lp_obj} -> Integer Optimum = \{gsol.objective} (x1=\{gsol.primal[0]}, x2=\{gsol.primal[1]})",
  )
  println(
    "    GMI Cuts Added = \{grep.cuts_added} | Dual Pivots = \{grep.dual_pivots} | Solved at Root = \{grep.solved_at_root}",
  )

  // 5. Authentic Netlib AFIRO benchmark + 5-City TSP-MTZ Warm-Start vs Cold-Start B&B
  let afiro_res = try {
    let am = @moonlp.parse_mps(@moonlp.netlib_afiro_mps())
    let asol = am.solve()
    let acert = am.certify_lp_solution(asol, tol=1.0e-5)
    "Optimum = \{asol.objective} (canonical -464.753143) | pivots=\{asol.iterations} | KKT=\{acert.certified_optimal}"
  } catch {
    @moonlp.ParseError(msg) => "Parse error: \{msg}"
  }
  println("\n[5] Authentic Netlib AFIRO Benchmark (28 rows, 32 cols)")
  println("    \{afiro_res}")

  let tsp = @moonlp.build_tsp5_mtz_model()
  let tsp_cold = @moonlp.solve_milp(tsp)
  let tsp_warm = @moonlp.solve_milp_warm(tsp)
  println("\n[6] 5-City TSP (MTZ Subtour Elimination MILP, 20 binary + 4 continuous vars)")
  println(
    "    Optimum Tour = \{tsp_warm.objective} | Cold B&B pivots=\{tsp_cold.iterations} vs Dual Warm-Start pivots=\{tsp_warm.iterations}",
  )

  // 7. Markowitz 6-Asset Mean-Variance Portfolio (Continuous Convex QP) vs
  //    Cardinality-Constrained Portfolio Selection (Mixed-Integer QP, MIQP)
  let mk_qp = @moonlp.build_markowitz_portfolio_qp(target_return=0.08)
  let mk_qsol = mk_qp.solve()
  let mk_qcert = mk_qp.certify_qp_solution(mk_qsol)
  println(
    "\n[7] Markowitz 6-Asset Mean-Variance Portfolio (Convex QP, target return >= 8%)",
  )
  println(
    "    Min 0.5*w^T*Sigma*w = \{mk_qsol.objective} | IPM iters=\{mk_qsol.iterations} | Wolfe KKT Certified=\{mk_qcert.certified_optimal} (gap=\{mk_qcert.duality_gap})",
  )
  let names = ["US_Eq", "Intl_Eq", "Corp_Bond", "Gov_Bond", "REIT", "Gold"]
  let q_weights = StringBuilder()
  for i = 0; i < 6; i = i + 1 {
    if i > 0 {
      q_weights.write_string(", ")
    }
    q_weights.write_string(names[i])
    q_weights.write_string("=")
    q_weights.write_string(mk_qsol.primal[i].to_string())
  }
  println("    Weights = [\{q_weights.to_string()}]")

  let mk_miqp = @moonlp.build_cardinality_portfolio_miqp(
    max_assets=3,
    min_weight=0.1,
    max_weight=0.5,
    target_return=0.08,
  )
  let mk_msol = mk_miqp.solve()
  let picked = StringBuilder()
  for i = 0; i < 6; i = i + 1 {
    if mk_msol.primal[6 + i] > 0.5 {
      if picked.to_string().length() > 0 {
        picked.write_string(", ")
      }
      picked.write_string(names[i])
      picked.write_string("=")
      picked.write_string(mk_msol.primal[i].to_string())
    }
  }
  println(
    "\n[8] Cardinality-Constrained Markowitz Portfolio (MIQP: max 3 of 6 assets, min hold 10%)",
  )
  println(
    "    Global Optimal 0.5*w^T*Sigma*w = \{mk_msol.objective} (vs QP relaxation \{mk_qsol.objective}) | B&B nodes=\{mk_msol.nodes_explored} | mip_gap=\{mk_msol.mip_gap}",
  )
  println("    Selected 3 Assets = [\{picked.to_string()}]")

  // 9. Constrained Linear Model Predictive Control (MPC Trajectory Tracking QP)
  let mpc = @moonlp.build_linear_mpc_qp(
    horizon=8,
    x0_pos=2.0,
    x0_vel=0.0,
    u_max=1.0,
  )
  let mpc_sol = mpc.solve()
  let mpc_cert = mpc.certify_qp_solution(mpc_sol)
  println(
    "\n[9] Constrained Linear MPC Trajectory Tracking (8-step horizon, 24 vars, 16 dynamics eqs)",
  )
  println(
    "    LQR Cost = \{mpc_sol.objective} | u_0 = \{mpc_sol.primal[16]} (saturated at -1.0) | p_8 = \{mpc_sol.primal[14]} | KKT=\{mpc_cert.certified_optimal}",
  )

  // 10. Soft-Margin SVM Dual Training (Rank-2 PSD 8x8 Gram Matrix QP)
  let svm = @moonlp.build_svm_dual_qp(c_reg=10.0)
  let svm_sol = svm.solve()
  let svm_cert = svm.certify_qp_solution(svm_sol)
  let mut sv_count = 0
  for i = 0; i < 8; i = i + 1 {
    if svm_sol.primal[i] > 1.0e-4 {
      sv_count = sv_count + 1
    }
  }
  println(
    "\n[10] Soft-Margin SVM Dual Training (8x8 Rank-2 PSD Gram Matrix QP)",
  )
  println(
    "     Dual Objective = \{svm_sol.objective} | Support Vectors = \{sv_count} | IPM iters=\{svm_sol.iterations} | KKT=\{svm_cert.certified_optimal}",
  )

  // 11. Thermal Unit Commitment & Quadratic Economic Dispatch (MIQP)
  let uc = @moonlp.build_unit_commitment_miqp(demand=140.0, reserve=20.0)
  let uc_sol = uc.solve()
  let uc_names = ["Coal_Base", "CCGT_Mid", "Gas_Flex", "Peaker_1", "Peaker_2"]
  let uc_on = StringBuilder()
  for g = 0; g < 5; g = g + 1 {
    if uc_sol.primal[5 + g] > 0.5 {
      if uc_on.to_string().length() > 0 {
        uc_on.write_string(", ")
      }
      uc_on.write_string(uc_names[g])
      uc_on.write_string("=")
      uc_on.write_string(uc_sol.primal[g].to_string())
      uc_on.write_string("MW")
    }
  }
  println(
    "\n[11] Thermal Unit Commitment & Quadratic Economic Dispatch (MIQP, 140MW + 20MW reserve)",
  )
  println(
    "     Min Total Cost = \{uc_sol.objective} | B&B nodes=\{uc_sol.nodes_explored} | mip_gap=\{uc_sol.mip_gap}",
  )
  println("     Committed Units = [\{uc_on.to_string()}]")

  // 12. L0 Best-Subset Sparse Ridge Regression (MIQP)
  let bsr = @moonlp.build_best_subset_ridge_miqp(
    max_features=2,
    ridge_lambda=0.05,
  )
  let bsr_sol = bsr.solve()
  println(
    "\n[12] L0 Best-Subset Sparse Ridge Regression (MIQP, max 2 of 5 features)",
  )
  println(
    "     Min Regularized Loss = \{bsr_sol.objective} | B&B nodes=\{bsr_sol.nodes_explored} | beta_0=\{bsr_sol.primal[0]}, beta_3=\{bsr_sol.primal[3]}",
  )

  // 13. Hybrid Switched-Thruster Deadzone MPC (MIQP-MPC)
  let hmpc = @moonlp.build_hybrid_mpc_miqp(
    horizon=5,
    max_firings=3,
    x0_pos=2.0,
    x0_vel=0.0,
    u_min=0.3,
    u_max=1.0,
  )
  let hmpc_sol = hmpc.solve()
  let fire_str = StringBuilder()
  for t = 0; t < 5; t = t + 1 {
    if t > 0 {
      fire_str.write_string(", ")
    }
    fire_str.write_string("t\{t}:")
    fire_str.write_string(
      if hmpc_sol.primal[15 + t] > 0.5 {
        hmpc_sol.primal[10 + t].to_string()
      } else {
        "OFF"
      },
    )
  }
  println(
    "\n[13] Hybrid Switched-Thruster Deadzone MPC (MIQP-MPC, u in {0} U [0.3, 1.0], <= 3 firings)",
  )
  println(
    "     Min Hybrid Cost = \{hmpc_sol.objective} | B&B nodes=\{hmpc_sol.nodes_explored} | Thruster = [\{fire_str.to_string()}]",
  )

  // 14. Smallest Enclosing Disk / 5G Base Station Coverage (SOCP)
  let disk = @moonlp.build_smallest_enclosing_disk_socp()
  let (dsol, dcert) = disk.solve_socp_certified()
  println("\n[14] Smallest Enclosing Disk / 5G Coverage Radius (SOCP, 6 Lorentz cones)")
  println(
    "     Min Radius = \{dsol.objective} | Center = (\{dsol.primal[0]}, \{dsol.primal[1]}) | Conic KKT=\{dcert.certified_optimal} (gap=\{dcert.duality_gap})",
  )

  // 15. Robust Ellipsoidal-Uncertainty Portfolio (SOCP)
  let rport = @moonlp.build_robust_portfolio_socp(kappa=1.5, max_weight=0.5)
  let (rpsol, rpcert) = rport.solve_socp_certified()
  println("\n[15] Robust Portfolio under Ellipsoidal Return Uncertainty (SOCP, kappa=1.5)")
  println(
    "     Max Worst-Case Return = \{rpsol.objective} | Risk ||Sigma^0.5 w||_2 = \{rpsol.primal[6]} | Conic KKT=\{rpcert.certified_optimal}",
  )

  // 16. Volatility-Capped Cardinality-Constrained Portfolio (MISOCP / MIQCQP)
  let vmisocp = @moonlp.build_vol_capped_cardinality_misocp(
    max_assets=3,
    vol_cap=0.105,
    min_weight=0.1,
    max_weight=0.5,
  )
  let vsol = vmisocp.solve()
  let vpick = StringBuilder()
  for i = 0; i < 6; i = i + 1 {
    if vsol.primal[6 + i] > 0.5 {
      if vpick.to_string().length() > 0 {
        vpick.write_string(", ")
      }
      vpick.write_string(names[i])
      vpick.write_string("=")
      vpick.write_string(vsol.primal[i].to_string())
    }
  }
  println(
    "\n[16] Volatility-Capped Cardinality Portfolio (MISOCP: max 3 assets, vol <= 10.5%)",
  )
  println(
    "     Max Return = \{vsol.objective} | B&B nodes=\{vsol.nodes_explored} | mip_gap=\{vsol.mip_gap}",
  )
  println("     Selected 3 Assets = [\{vpick.to_string()}]")

  // 17. 95% Gaussian Chance-Constrained Production Planning (SOCP)
  let cc = @moonlp.build_chance_constrained_socp(z_alpha=1.644854)
  let (ccsol, cccert) = cc.solve_socp_certified()
  println("\n[17] 95% Gaussian Chance-Constrained Production Planning (SOCP)")
  println(
    "     Max Safe Profit = \{ccsol.objective} | x = (\{ccsol.primal[0]}, \{ccsol.primal[1]}, \{ccsol.primal[2]}) | Conic KKT=\{cccert.certified_optimal}",
  )

  // 18. 2D Isotropic Circular-Thrust Rocket Guidance (SOCP & Switched 2-Burn MISOCP)
  let r_cont = @moonlp.build_isotropic_thrust_mpc_socp(horizon=4, u_max=1.0)
  let r_csol = r_cont.solve()
  let r_sw = @moonlp.build_switched_isotropic_mpc_misocp(
    horizon=4,
    max_burns=2,
    u_max=1.0,
  )
  let r_ssol = r_sw.solve()
  println(
    "\n[18] 2D Isotropic Circular-Thrust Rocket Guidance (SOCP vs 2-Burn Switched MISOCP)",
  )
  println(
    "     Continuous 4-step SOCP Cost = \{r_csol.objective} | Switched <=2-Burn MISOCP Cost = \{r_ssol.objective} (B&B nodes=\{r_ssol.nodes_explored})",
  )
}
