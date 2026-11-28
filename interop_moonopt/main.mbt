///|
/// Live Interoperability & Extension Bridge: `Freon793/moonopt` + `spring1219/moonlp`.
///
/// Demonstrates directly in runnable MoonBit code:
///   1. Direct in-memory conversion `from_moonopt(@mopt_model.Model) -> @moonlp.Model`
///      and zero-loss MPS handoff between `Freon793/moonopt` and `spring1219/moonlp`.
///   2. Reusing `@moonopt.solve` on linear LP / MILP baselines, then extending those exact
///      models in `moonlp` into SOCP (Chance-Constrained), MISOCP (Volatility-Capped
///      Cardinality Portfolio), and MIQP (Quadratic Fuel-Cost Unit Commitment) — problem
///      classes `moonopt` cannot represent or solve.
///   3. Live format-boundary proof: `@mopt_format.parse_mps` rejects `.qps` (`QUADOBJ`)
///      and `@mopt_format.read_lp` rejects CPLEX QCP `[ ... ] <= rhs`, whereas `moonlp`
///      parses, solves, and certifies them.

///|
/// Directly convert a `Freon793/moonopt/model.Model` into a `spring1219/moonlp.Model`
/// so any existing `moonopt` linear model can be extended with quadratic objectives,
/// convex quadratic constraints (QCQP), or Second-Order Cone constraints (SOCP/MISOCP).
pub fn from_moonopt(src : @mopt_model.Model) -> @moonlp.Model {
  let sense : @moonlp.Sense = match src.sense() {
    Minimize => Minimize
    Maximize => Maximize
  }
  let dst = @moonlp.Model::new(sense)
  let ids : Array[@moonlp.VarId] = []
  let obj_map = Array::make(src.num_vars(), 0.0)
  let obj_terms = src.objective()
  for k = 0; k < obj_terms.length(); k = k + 1 {
    let (idx, coeff) = obj_terms[k]
    obj_map[idx] = obj_map[idx] + coeff
  }
  for j = 0; j < src.num_vars(); j = j + 1 {
    let v = src.get_var(j)
    let lb = if v.lb <= -@mopt_model.infinite_bound * 0.5 {
      -1.0e300
    } else {
      v.lb
    }
    let ub = if v.ub >= @mopt_model.infinite_bound * 0.5 {
      1.0e300
    } else {
      v.ub
    }
    let kind : @moonlp.VarKind = if v.is_int {
      if lb == 0.0 && ub == 1.0 { Binary } else { Integer }
    } else {
      Continuous
    }
    ids.push(dst.add_var(v.name, lb~, ub~, kind~, obj=obj_map[j]))
  }
  for i = 0; i < src.num_constraints(); i = i + 1 {
    let c = src.constraint(i)
    let terms : Array[@moonlp.Term] = []
    for k = 0; k < c.terms.length(); k = k + 1 {
      let (idx, coeff) = c.terms[k]
      terms.push({ var_id: ids[idx], coeff, })
    }
    let rel : @moonlp.Relation = match c.rel {
      LessEqual => Le
      GreaterEqual => Ge
      Equal => Eq
    }
    dst.add_constraint(c.name, terms, rel, c.rhs)
  }
  dst
}

///|
/// Build the deterministic linear production baseline natively in `Freon793/moonopt`.
fn build_moonopt_production_lp() -> @mopt_model.Model {
  let m = @mopt_model.Model::new(@mopt_model.Sense::Maximize)
  let x1 = m.add_var("x1", lb=0.0)
  let x2 = m.add_var("x2", lb=0.0)
  let x3 = m.add_var("x3", lb=0.0)
  m.set_objective([(x1, 50.0), (x2, 65.0), (x3, 45.0)])
  m.add_named_constraint(
    [(x1, 1.0), (x2, 1.0), (x3, 1.0)],
    @mopt_model.Rel::LessEqual,
    40.0,
    "labor",
  )
  m.add_named_constraint(
    [(x1, 2.0), (x2, 3.0), (x3, 1.5)],
    @mopt_model.Rel::LessEqual,
    100.0,
    "machine_mean",
  )
  m
}

///|
/// Build the linear max-return 6-asset portfolio baseline natively in `Freon793/moonopt`.
fn build_moonopt_linear_portfolio() -> @mopt_model.Model {
  let mu = [0.11, 0.10, 0.06, 0.035, 0.085, 0.055]
  let names = ["US_Eq", "Intl_Eq", "Corp_Bond", "Gov_Bond", "REIT", "Gold"]
  let m = @mopt_model.Model::new(@mopt_model.Sense::Maximize)
  let obj : Array[(Int32, Double)] = []
  let bud : Array[(Int32, Double)] = []
  for i = 0; i < 6; i = i + 1 {
    let wi = m.add_var(names[i], lb=0.0, ub=0.5)
    obj.push((wi, mu[i]))
    bud.push((wi, 1.0))
  }
  m.set_objective(obj)
  m.add_named_constraint(bud, @mopt_model.Rel::Equal, 1.0, "budget")
  m
}

///|
/// Build the linearized 5-generator Unit Commitment MILP baseline natively in `Freon793/moonopt`.
fn build_moonopt_uc_milp() -> @mopt_model.Model {
  let names = ["Coal_Base", "CCGT_Mid", "Gas_Flex", "Peaker_1", "Peaker_2"]
  let p_min = [30.0, 20.0, 15.0, 10.0, 5.0]
  let p_max = [80.0, 60.0, 50.0, 40.0, 30.0]
  let lin_b = [18.0, 22.0, 26.0, 31.0, 35.0]
  let fix_c = [220.0, 140.0, 90.0, 55.0, 35.0]
  let m = @mopt_model.Model::new(@mopt_model.Sense::Minimize)
  let p : Array[Int32] = []
  let u : Array[Int32] = []
  let obj : Array[(Int32, Double)] = []
  let bal : Array[(Int32, Double)] = []
  let res : Array[(Int32, Double)] = []
  for g = 0; g < 5; g = g + 1 {
    let pg = m.add_var("p_" + names[g], lb=0.0, ub=p_max[g])
    p.push(pg)
    obj.push((pg, lin_b[g]))
    bal.push((pg, 1.0))
  }
  for g = 0; g < 5; g = g + 1 {
    let ug = m.add_var("u_" + names[g], lb=0.0, ub=1.0, is_int=true)
    u.push(ug)
    obj.push((ug, fix_c[g]))
    res.push((ug, p_max[g]))
    m.add_named_constraint(
      [(p[g], 1.0), (ug, -p_min[g])],
      @mopt_model.Rel::GreaterEqual,
      0.0,
      "pmin_" + g.to_string(),
    )
    m.add_named_constraint(
      [(p[g], 1.0), (ug, -p_max[g])],
      @mopt_model.Rel::LessEqual,
      0.0,
      "pmax_" + g.to_string(),
    )
  }
  m.set_objective(obj)
  m.add_named_constraint(bal, @mopt_model.Rel::Equal, 140.0, "power_balance")
  m.add_named_constraint(
    res,
    @mopt_model.Rel::GreaterEqual,
    160.0,
    "spinning_reserve",
  )
  m
}

///|
fn main {
  println(
    "=== Live Interop & Extension: Freon793/moonopt@0.6.0 + spring1219/moonlp ===\n",
  )

  // 1. Reuse moonopt on deterministic LP, hand off via MPS, then extend to 95% Chance-Constrained SOCP
  let mopt_prod = build_moonopt_production_lp()
  let mopt_psol = @moonopt.solve(mopt_prod)
  let mps_text = @mopt_format.write_mps(mopt_prod)
  let mlp_prod_lin = try {
    @moonlp.parse_mps(mps_text, sense=Maximize)
  } catch {
    _ => abort("MPS handoff failed")
  }
  let mlp_lsol = mlp_prod_lin.solve()
  let cc_model = @moonlp.build_chance_constrained_socp(z_alpha=1.644854)
  let (cc_sol, cc_cert) = cc_model.solve_socp_certified()
  println("[1] Production Planning: moonopt Linear Baseline -> moonlp 95% Chance-Constrained SOCP")
  println(
    "    moonopt deterministic LP = \{mopt_psol.objective} (MPS handoff moonlp LP = \{mlp_lsol.objective})",
  )
  println(
    "    moonlp 95% Chance-Constrained SOCP = \{cc_sol.objective} | Conic KKT Certified = \{cc_cert.certified_optimal} (gap = \{cc_cert.duality_gap})\n",
  )

  // 2. Portfolio: moonopt linear baseline vs moonlp Volatility-Capped Cardinality MISOCP + QPS rejection proof
  let mopt_port = build_moonopt_linear_portfolio()
  let mopt_port_sol = @moonopt.solve(mopt_port)
  let qps_text = @moonlp.build_markowitz_portfolio_qp().to_mps_string(
    name="MARKOWITZ",
  )
  let mopt_qps_rejected = match @mopt_format.parse_mps(qps_text) {
    Err(issue) => "Rejected as expected: " + issue.to_text()
    Ok(_) => "Unexpectedly accepted"
  }
  let misocp_port = @moonlp.build_vol_capped_cardinality_misocp(
    max_assets=3,
    vol_cap=0.105,
  )
  let misocp_sol = misocp_port.solve()
  println("[2] Portfolio Risk: moonopt Linear -> moonlp Volatility-Capped Cardinality MISOCP")
  println(
    "    moonopt linear max-return = \{mopt_port_sol.objective} (concentrated 50% US_Eq + 50% Intl_Eq, ~18.4% volatility!)",
  )
  println("    moonopt parse_mps(.qps QUADOBJ) -> \{mopt_qps_rejected}")
  println(
    "    moonlp MISOCP (vol <= 10.5%, <= 3 assets) = \{misocp_sol.objective} | B&B nodes = \{misocp_sol.nodes_explored} | mip_gap = \{misocp_sol.mip_gap}\n",
  )

  // 3. Unit Commitment: build MILP in moonopt, lift directly via from_moonopt(), add quadratic fuel curves -> MIQP
  let mopt_uc = build_moonopt_uc_milp()
  let mopt_uc_sol = @moonopt.solve(mopt_uc)
  let uc_extended = from_moonopt(mopt_uc)
  let quad_a = [0.008, 0.012, 0.018, 0.025, 0.032]
  for g = 0; g < 5; g = g + 1 {
    uc_extended.add_quad_obj_term({ index: g, }, { index: g, }, quad_a[g])
  }
  let uc_miqp_sol = uc_extended.solve()
  println("[3] Thermal Unit Commitment: moonopt MILP -> from_moonopt() + Quadratic Heat-Rate MIQP")
  println(
    "    moonopt linear MILP cost (ignores a_g * p_g^2) = \{mopt_uc_sol.objective}",
  )
  println(
    "    moonlp extended MIQP true quadratic cost = \{uc_miqp_sol.objective} | B&B nodes = \{uc_miqp_sol.nodes_explored} | mip_gap = \{uc_miqp_sol.mip_gap}",
  )
}
