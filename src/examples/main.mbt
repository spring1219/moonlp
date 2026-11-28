///|
/// Runnable demonstration of moonlp v1.1.0:
/// Leads with MoonBit-ecosystem-unique Conic (SOCP / QCQP / MISOCP) and Mixed-Integer
/// Quadratic (MIQP / QP) models, followed by LP/MILP compatibility baselines.

///|
fn main {
  println(
    "=== moonlp v1.1.0: Pure MoonBit SOCP / MISOCP / QCQP / MIQP Optimization Kernel ===\n",
  )

  let names = ["US_Eq", "Intl_Eq", "Corp_Bond", "Gov_Bond", "REIT", "Gold"]

  // [1] Volatility-Capped Cardinality-Constrained Portfolio (MISOCP / MIQCQP)
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
    "[1] Volatility-Capped Cardinality Portfolio (MISOCP: max 3 assets, ||Sigma^0.5 w||_2 <= 10.5%)",
  )
  println(
    "    Max Return = \{vsol.objective} | B&B nodes=\{vsol.nodes_explored} | mip_gap=\{vsol.mip_gap}",
  )
  println("    Selected 3 Assets = [\{vpick.to_string()}]\n")

  // [2] 2D Isotropic Circular-Thrust Rocket Guidance (SOCP & Switched 2-Burn MISOCP)
  let r_cont = @moonlp.build_isotropic_thrust_mpc_socp(horizon=4, u_max=1.0)
  let (r_csol, r_ccert) = r_cont.solve_socp_certified(cert_tol=1.0e-4)
  let r_sw = @moonlp.build_switched_isotropic_mpc_misocp(
    horizon=4,
    max_burns=2,
    u_max=1.0,
  )
  let r_ssol = r_sw.solve()
  println(
    "[2] 2D Isotropic Circular-Thrust Rocket Guidance (Quad Obj + Lorentz Thrust Cones)",
  )
  println(
    "    Continuous 4-step SOCP = \{r_csol.objective} (Wolfe-Conic KKT=\{r_ccert.certified_optimal}, gap=\{r_ccert.duality_gap})",
  )
  println(
    "    Switched <=2-Burn MISOCP = \{r_ssol.objective} (B&B nodes=\{r_ssol.nodes_explored})\n",
  )

  // [3] Thermal Unit Commitment & Quadratic Economic Dispatch (MIQP)
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
    "[3] Thermal Unit Commitment & Quadratic Economic Dispatch (MIQP, 140MW + 20MW reserve)",
  )
  println(
    "    Min Total Cost = \{uc_sol.objective} | B&B nodes=\{uc_sol.nodes_explored} | mip_gap=\{uc_sol.mip_gap}",
  )
  println("    Committed Units = [\{uc_on.to_string()}]\n")

  // [4] Robust Portfolio under Ellipsoidal Return Uncertainty (SOCP)
  let rport = @moonlp.build_robust_portfolio_socp(kappa=1.5, max_weight=0.5)
  let (rpsol, rpcert) = rport.solve_socp_certified()
  println("[4] Robust Portfolio under Ellipsoidal Return Uncertainty (SOCP, kappa=1.5)")
  println(
    "    Max Worst-Case Return = \{rpsol.objective} | Risk = \{rpsol.primal[6]} | Conic KKT=\{rpcert.certified_optimal}\n",
  )

  // [5] 95% Gaussian Chance-Constrained Production Planning (SOCP)
  let cc = @moonlp.build_chance_constrained_socp(z_alpha=1.644854)
  let (ccsol, cccert) = cc.solve_socp_certified()
  println("[5] 95% Gaussian Chance-Constrained Production Planning (SOCP)")
  println(
    "    Max Safe Profit = \{ccsol.objective} | x = (\{ccsol.primal[0]}, \{ccsol.primal[1]}, \{ccsol.primal[2]}) | Conic KKT=\{cccert.certified_optimal}\n",
  )

  // [6] Smallest Enclosing Disk / 5G Coverage Radius (SOCP)
  let disk = @moonlp.build_smallest_enclosing_disk_socp()
  let (dsol, dcert) = disk.solve_socp_certified()
  println("[6] Smallest Enclosing Disk / 5G Coverage Radius (SOCP, 6 Lorentz cones)")
  println(
    "    Min Radius = \{dsol.objective} | Center = (\{dsol.primal[0]}, \{dsol.primal[1]}) | Conic KKT=\{dcert.certified_optimal} (gap=\{dcert.duality_gap})\n",
  )

  // [7] Cardinality-Constrained Markowitz (MIQP) vs Continuous QP + Wolfe KKT
  let mk_qp = @moonlp.build_markowitz_portfolio_qp(target_return=0.08)
  let mk_qsol = mk_qp.solve()
  let mk_qcert = mk_qp.certify_qp_solution(mk_qsol)
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
  println("[7] Markowitz Portfolio: Continuous QP vs Cardinality-Constrained MIQP (<= 3 assets)")
  println(
    "    Continuous QP 0.5*w^T*Sigma*w = \{mk_qsol.objective} | Wolfe KKT=\{mk_qcert.certified_optimal} (gap=\{mk_qcert.duality_gap})",
  )
  println(
    "    MIQP (<=3 assets) = \{mk_msol.objective} | Selected = [\{picked.to_string()}]\n",
  )

  // [8] L0 Best-Subset Sparse Ridge Regression (MIQP)
  let bsr = @moonlp.build_best_subset_ridge_miqp(
    max_features=2,
    ridge_lambda=0.05,
  )
  let bsr_sol = bsr.solve()
  println("[8] L0 Best-Subset Sparse Ridge Regression (MIQP, max 2 of 5 features)")
  println(
    "    Min Loss = \{bsr_sol.objective} | B&B nodes=\{bsr_sol.nodes_explored} | beta_0=\{bsr_sol.primal[0]}, beta_3=\{bsr_sol.primal[3]}\n",
  )

  // [9] Hybrid Switched-Thruster Deadzone MPC (MIQP-MPC) & Soft-Margin SVM Dual (PSD QP)
  let hmpc = @moonlp.build_hybrid_mpc_miqp(
    horizon=5,
    max_firings=3,
    x0_pos=2.0,
    x0_vel=0.0,
    u_min=0.3,
    u_max=1.0,
  )
  let hmpc_sol = hmpc.solve()
  let svm = @moonlp.build_svm_dual_qp(c_reg=10.0)
  let svm_sol = svm.solve()
  let svm_cert = svm.certify_qp_solution(svm_sol)
  println("[9] Hybrid Deadzone MPC (MIQP-MPC) & Soft-Margin SVM Dual (Rank-2 PSD QP)")
  println(
    "    Hybrid MPC Cost = \{hmpc_sol.objective} (u_0=\{hmpc_sol.primal[10]}, u_1=\{hmpc_sol.primal[11]}, t2..4=OFF)",
  )
  println(
    "    SVM Dual Obj = \{svm_sol.objective} | IPM iters=\{svm_sol.iterations} | Wolfe KKT=\{svm_cert.certified_optimal}\n",
  )

  // [10] Compatibility Baseline: Authentic Netlib AFIRO LP & 5-City TSP-MTZ MILP
  let afiro_res = try {
    let am = @moonlp.parse_mps(@moonlp.netlib_afiro_mps())
    let asol = am.solve()
    let acert = am.certify_lp_solution(asol, tol=1.0e-5)
    "Optimum = \{asol.objective} | KKT=\{acert.certified_optimal}"
  } catch {
    @moonlp.ParseError(msg) => "Parse error: \{msg}"
  }
  let tsp_warm = @moonlp.solve_milp_warm(@moonlp.build_tsp5_mtz_model())
  println("[10] Compatibility Baseline (LP/MILP; for large sparse LP/MILP see interop_moonopt/)")
  println("     Netlib AFIRO LP: \{afiro_res}")
  println(
    "     5-City TSP-MTZ MILP: Optimum Tour = \{tsp_warm.objective} (Dual Warm-Start pivots=\{tsp_warm.iterations})",
  )
}
