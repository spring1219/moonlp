///|
/// Concrete Real-World Second-Order Cone (SOCP), Quadratically Constrained (QCQP),
/// and Mixed-Integer Second-Order Cone (MISOCP) Models (v2.0.0):
/// Every model below requires nonlinear Euclidean/Lorentz cone constraints (and/or
/// integer-coupled conic constraints), which NEITHER pure LP/MILP solvers (moonopt)
/// NOR linearly-constrained continuous QP solvers (moonqp) can directly express or solve.

///|
/// 1. Smallest Enclosing Disk / 5G Base-Station Minimum Coverage Radius (SOCP):
///   min  r
///   s.t. || [x - p_ix, y - p_iy] ||_2 <= r,   i = 0..5
pub fn build_smallest_enclosing_disk_socp() -> Model {
  let pts = [
    (0.0, 0.0),
    (4.0, 0.0),
    (2.0, 3.5),
    (1.0, 1.0),
    (3.0, 1.5),
    (2.0, 0.5),
  ]
  let m = Model::new(Minimize)
  let x = m.add_continuous_var("cx", lb=-1.0e300, ub=1.0e300)
  let y = m.add_continuous_var("cy", lb=-1.0e300, ub=1.0e300)
  let r = m.add_continuous_var("radius", lb=0.0, obj=1.0)
  for i = 0; i < pts.length(); i = i + 1 {
    let (px, py) = pts[i]
    let dx = LinearExpr::from_var(x).add_const(-px)
    let dy = LinearExpr::from_var(y).add_const(-py)
    m.add_soc_constraint(
      "cover_" + i.to_string(),
      [dx, dy],
      LinearExpr::from_var(r),
    )
  }
  m
}

///|
/// 2. Robust Portfolio Optimization under Ellipsoidal Return Uncertainty (SOCP):
///   max  mu^T w - kappa * risk
///   s.t. || Sigma^{1/2} w ||_2 <= risk,  sum(w_i) = 1,  0 <= w_i <= max_weight
/// Maximizes exact worst-case portfolio return over the confidence ellipsoid
/// { tilde_mu : (tilde_mu - mu)^T Sigma^{-1} (tilde_mu - mu) <= kappa^2 }.
pub fn build_robust_portfolio_socp(
  kappa~ : Double = 1.5,
  max_weight~ : Double = 0.5,
) -> Model {
  let (mu, sigma) = markowitz6_data()
  let names = ["US_Eq", "Intl_Eq", "Corp_Bond", "Gov_Bond", "REIT", "Gold"]
  let m = Model::new(Maximize)
  let w : Array[VarId] = []
  let budget : Array[Term] = []
  for i = 0; i < 6; i = i + 1 {
    let wi = m.add_continuous_var(names[i], lb=0.0, ub=max_weight, obj=mu[i])
    w.push(wi)
    budget.push({ var_id: wi, coeff: 1.0, })
  }
  let risk = m.add_continuous_var("risk", lb=0.0, obj=-kappa)
  m.add_constraint("budget", budget, Eq, 1.0)

  match psd_factor_rows(sigma, 6) {
    Some(frows) => {
      let norm_exprs : Array[LinearExpr] = []
      for r = 0; r < frows.length(); r = r + 1 {
        let e = LinearExpr::zero()
        for j = 0; j < 6; j = j + 1 {
          if frows[r][j] != 0.0 {
            let _ = e.add_term(w[j], frows[r][j])
          }
        }
        norm_exprs.push(e)
      }
      m.add_soc_constraint(
        "ellipsoidal_risk",
        norm_exprs,
        LinearExpr::from_var(risk),
      )
    }
    None => ()
  }
  m
}

///|
/// 3. Volatility-Capped Cardinality-Constrained Portfolio Selection (MISOCP / MIQCQP):
///   max  mu^T w
///   s.t. || Sigma^{1/2} w ||_2 <= vol_cap   (nonlinear conic volatility cap!)
///        sum(w_i) = 1
///        min_weight * z_i <= w_i <= max_weight * z_i,  z_i in {0, 1}
///        sum(z_i) <= max_assets
pub fn build_vol_capped_cardinality_misocp(
  max_assets~ : Int32 = 3,
  vol_cap~ : Double = 0.105,
  min_weight~ : Double = 0.1,
  max_weight~ : Double = 0.5,
) -> Model {
  let (mu, sigma) = markowitz6_data()
  let names = ["US_Eq", "Intl_Eq", "Corp_Bond", "Gov_Bond", "REIT", "Gold"]
  let m = Model::new(Maximize)
  let w : Array[VarId] = []
  let z : Array[VarId] = []
  let budget : Array[Term] = []
  let card : Array[Term] = []

  for i = 0; i < 6; i = i + 1 {
    let wi = m.add_continuous_var(
      "w_" + names[i],
      lb=0.0,
      ub=max_weight,
      obj=mu[i],
    )
    w.push(wi)
    budget.push({ var_id: wi, coeff: 1.0, })
  }
  for i = 0; i < 6; i = i + 1 {
    let zi = m.add_binary_var("z_" + names[i])
    z.push(zi)
    card.push({ var_id: zi, coeff: 1.0, })
    m.add_constraint(
      "link_lb_" + i.to_string(),
      [{ var_id: w[i], coeff: 1.0, }, { var_id: zi, coeff: -min_weight, }],
      Ge,
      0.0,
    )
    m.add_constraint(
      "link_ub_" + i.to_string(),
      [{ var_id: w[i], coeff: 1.0, }, { var_id: zi, coeff: -max_weight, }],
      Le,
      0.0,
    )
  }
  m.add_constraint("budget", budget, Eq, 1.0)
  m.add_constraint("cardinality", card, Le, max_assets.to_double())

  match psd_factor_rows(sigma, 6) {
    Some(frows) => {
      let norm_exprs : Array[LinearExpr] = []
      for r = 0; r < frows.length(); r = r + 1 {
        let e = LinearExpr::zero()
        for j = 0; j < 6; j = j + 1 {
          if frows[r][j] != 0.0 {
            let _ = e.add_term(w[j], frows[r][j])
          }
        }
        norm_exprs.push(e)
      }
      m.add_soc_constraint(
        "volatility_cap",
        norm_exprs,
        LinearExpr::zero().add_const(vol_cap),
      )
    }
    None => ()
  }
  m
}

///|
/// 4. Gaussian Chance-Constrained Production Planning (SOCP):
///   max  50 x1 + 65 x2 + 45 x3
///   s.t. Pr(tilde_a^T x <= 100) >= 95%  where tilde_a ~ N([2, 3, 1.5], diag(0.5^2, 0.8^2, 0.4^2))
///        <=> 2 x1 + 3 x2 + 1.5 x3 + 1.644854 * || [0.5 x1, 0.8 x2, 0.4 x3] ||_2 <= 100
///        x1 + x2 + x3 <= 40,  x_i >= 0
pub fn build_chance_constrained_socp(z_alpha~ : Double = 1.644854) -> Model {
  let m = Model::new(Maximize)
  let x1 = m.add_continuous_var("x1", lb=0.0, obj=50.0)
  let x2 = m.add_continuous_var("x2", lb=0.0, obj=65.0)
  let x3 = m.add_continuous_var("x3", lb=0.0, obj=45.0)
  m.add_constraint(
    "labor",
    [
      { var_id: x1, coeff: 1.0, },
      { var_id: x2, coeff: 1.0, },
      { var_id: x3, coeff: 1.0, },
    ],
    Le,
    40.0,
  )
  // || [z_alpha*0.5 x1, z_alpha*0.8 x2, z_alpha*0.4 x3] ||_2 <= 100 - (2 x1 + 3 x2 + 1.5 x3)
  let rhs = LinearExpr::zero()
    .add_const(100.0)
    .add_term(x1, -2.0)
    .add_term(x2, -3.0)
    .add_term(x3, -1.5)
  m.add_soc_constraint(
    "chance_95pct",
    [
      LinearExpr::from_var(x1, coeff=z_alpha * 0.5),
      LinearExpr::from_var(x2, coeff=z_alpha * 0.8),
      LinearExpr::from_var(x3, coeff=z_alpha * 0.4),
    ],
    rhs,
  )
  m
}

///|
/// 5a. 2D Isotropic Circular-Thrust Rocket/Drone Guidance MPC (SOCP + Quadratic Objective):
/// Planar double integrator from (2.0, 1.5) -> (0, 0) with exact Euclidean circular thrust
/// limit || (u_x, u_y) ||_2 <= u_max (unlike box constraints which distort diagonal thrust by 41%).
pub fn build_isotropic_thrust_mpc_socp(
  horizon~ : Int32 = 4,
  u_max~ : Double = 1.0,
) -> Model {
  let dt = 0.5
  let h2 = 0.5 * dt * dt
  let m = Model::new(Minimize)
  let px : Array[VarId] = []
  let py : Array[VarId] = []
  let vx : Array[VarId] = []
  let vy : Array[VarId] = []
  let ux : Array[VarId] = []
  let uy : Array[VarId] = []

  for t = 1; t <= horizon; t = t + 1 {
    let pxt = m.add_continuous_var("px_" + t.to_string(), lb=-1.0e300, ub=1.0e300)
    let pyt = m.add_continuous_var("py_" + t.to_string(), lb=-1.0e300, ub=1.0e300)
    let vxt = m.add_continuous_var("vx_" + t.to_string(), lb=-2.0, ub=2.0)
    let vyt = m.add_continuous_var("vy_" + t.to_string(), lb=-2.0, ub=2.0)
    px.push(pxt)
    py.push(pyt)
    vx.push(vxt)
    vy.push(vyt)
    let qp = if t == horizon { 20.0 } else { 10.0 }
    let qv = if t == horizon { 4.0 } else { 2.0 }
    m.add_quad_obj_term(pxt, pxt, 0.5 * qp)
    m.add_quad_obj_term(pyt, pyt, 0.5 * qp)
    m.add_quad_obj_term(vxt, vxt, 0.5 * qv)
    m.add_quad_obj_term(vyt, vyt, 0.5 * qv)
  }
  for t = 0; t < horizon; t = t + 1 {
    let uxt = m.add_continuous_var("ux_" + t.to_string(), lb=-u_max, ub=u_max)
    let uyt = m.add_continuous_var("uy_" + t.to_string(), lb=-u_max, ub=u_max)
    ux.push(uxt)
    uy.push(uyt)
    m.add_quad_obj_term(uxt, uxt, 0.25)
    m.add_quad_obj_term(uyt, uyt, 0.25)
    // Exact 2D circular thrust cone: || [ux_t, uy_t] ||_2 <= u_max
    m.add_soc_constraint(
      "thrust_circle_" + t.to_string(),
      [LinearExpr::from_var(uxt), LinearExpr::from_var(uyt)],
      LinearExpr::zero().add_const(u_max),
    )
  }
  for t = 0; t < horizon; t = t + 1 {
    if t == 0 {
      m.add_constraint(
        "dx0",
        [{ var_id: px[0], coeff: 1.0, }, { var_id: ux[0], coeff: -h2, }],
        Eq,
        2.0,
      )
      m.add_constraint(
        "dy0",
        [{ var_id: py[0], coeff: 1.0, }, { var_id: uy[0], coeff: -h2, }],
        Eq,
        1.5,
      )
      m.add_constraint(
        "dvx0",
        [{ var_id: vx[0], coeff: 1.0, }, { var_id: ux[0], coeff: -dt, }],
        Eq,
        0.0,
      )
      m.add_constraint(
        "dvy0",
        [{ var_id: vy[0], coeff: 1.0, }, { var_id: uy[0], coeff: -dt, }],
        Eq,
        0.0,
      )
    } else {
      m.add_constraint(
        "dx" + t.to_string(),
        [
          { var_id: px[t], coeff: 1.0, },
          { var_id: px[t - 1], coeff: -1.0, },
          { var_id: vx[t - 1], coeff: -dt, },
          { var_id: ux[t], coeff: -h2, },
        ],
        Eq,
        0.0,
      )
      m.add_constraint(
        "dy" + t.to_string(),
        [
          { var_id: py[t], coeff: 1.0, },
          { var_id: py[t - 1], coeff: -1.0, },
          { var_id: vy[t - 1], coeff: -dt, },
          { var_id: uy[t], coeff: -h2, },
        ],
        Eq,
        0.0,
      )
      m.add_constraint(
        "dvx" + t.to_string(),
        [
          { var_id: vx[t], coeff: 1.0, },
          { var_id: vx[t - 1], coeff: -1.0, },
          { var_id: ux[t], coeff: -dt, },
        ],
        Eq,
        0.0,
      )
      m.add_constraint(
        "dvy" + t.to_string(),
        [
          { var_id: vy[t], coeff: 1.0, },
          { var_id: vy[t - 1], coeff: -1.0, },
          { var_id: uy[t], coeff: -dt, },
        ],
        Eq,
        0.0,
      )
    }
  }
  m
}

///|
/// 5b. Switched 2D Circular-Thrust Rocket Guidance MPC (MISOCP):
/// Extends 5a with binary engine ignition switches z_t in {0, 1}, conic thrust
/// coupling || [ux_t, uy_t] ||_2 <= u_max * z_t, and maximum burn budget sum(z_t) <= max_burns.
pub fn build_switched_isotropic_mpc_misocp(
  horizon~ : Int32 = 4,
  max_burns~ : Int32 = 2,
  u_max~ : Double = 1.0,
) -> Model {
  let dt = 0.5
  let h2 = 0.5 * dt * dt
  let m = Model::new(Minimize)
  let px : Array[VarId] = []
  let py : Array[VarId] = []
  let vx : Array[VarId] = []
  let vy : Array[VarId] = []
  let ux : Array[VarId] = []
  let uy : Array[VarId] = []
  let burn_terms : Array[Term] = []

  for t = 1; t <= horizon; t = t + 1 {
    let pxt = m.add_continuous_var("px_" + t.to_string(), lb=-1.0e300, ub=1.0e300)
    let pyt = m.add_continuous_var("py_" + t.to_string(), lb=-1.0e300, ub=1.0e300)
    let vxt = m.add_continuous_var("vx_" + t.to_string(), lb=-2.0, ub=2.0)
    let vyt = m.add_continuous_var("vy_" + t.to_string(), lb=-2.0, ub=2.0)
    px.push(pxt)
    py.push(pyt)
    vx.push(vxt)
    vy.push(vyt)
    let qp = if t == horizon { 20.0 } else { 10.0 }
    let qv = if t == horizon { 4.0 } else { 2.0 }
    m.add_quad_obj_term(pxt, pxt, 0.5 * qp)
    m.add_quad_obj_term(pyt, pyt, 0.5 * qp)
    m.add_quad_obj_term(vxt, vxt, 0.5 * qv)
    m.add_quad_obj_term(vyt, vyt, 0.5 * qv)
  }
  for t = 0; t < horizon; t = t + 1 {
    let uxt = m.add_continuous_var("ux_" + t.to_string(), lb=-u_max, ub=u_max)
    let uyt = m.add_continuous_var("uy_" + t.to_string(), lb=-u_max, ub=u_max)
    ux.push(uxt)
    uy.push(uyt)
    m.add_quad_obj_term(uxt, uxt, 0.25)
    m.add_quad_obj_term(uyt, uyt, 0.25)
  }
  for t = 0; t < horizon; t = t + 1 {
    let zt = m.add_binary_var("z_" + t.to_string(), obj=0.5)
    burn_terms.push({ var_id: zt, coeff: 1.0, })
    // Conic engine switch: || [ux_t, uy_t] ||_2 <= u_max * z_t
    m.add_soc_constraint(
      "thrust_switch_" + t.to_string(),
      [LinearExpr::from_var(ux[t]), LinearExpr::from_var(uy[t])],
      LinearExpr::from_var(zt, coeff=u_max),
    )
  }
  m.add_constraint("max_burns", burn_terms, Le, max_burns.to_double())

  for t = 0; t < horizon; t = t + 1 {
    if t == 0 {
      m.add_constraint(
        "dx0",
        [{ var_id: px[0], coeff: 1.0, }, { var_id: ux[0], coeff: -h2, }],
        Eq,
        2.0,
      )
      m.add_constraint(
        "dy0",
        [{ var_id: py[0], coeff: 1.0, }, { var_id: uy[0], coeff: -h2, }],
        Eq,
        1.5,
      )
      m.add_constraint(
        "dvx0",
        [{ var_id: vx[0], coeff: 1.0, }, { var_id: ux[0], coeff: -dt, }],
        Eq,
        0.0,
      )
      m.add_constraint(
        "dvy0",
        [{ var_id: vy[0], coeff: 1.0, }, { var_id: uy[0], coeff: -dt, }],
        Eq,
        0.0,
      )
    } else {
      m.add_constraint(
        "dx" + t.to_string(),
        [
          { var_id: px[t], coeff: 1.0, },
          { var_id: px[t - 1], coeff: -1.0, },
          { var_id: vx[t - 1], coeff: -dt, },
          { var_id: ux[t], coeff: -h2, },
        ],
        Eq,
        0.0,
      )
      m.add_constraint(
        "dy" + t.to_string(),
        [
          { var_id: py[t], coeff: 1.0, },
          { var_id: py[t - 1], coeff: -1.0, },
          { var_id: vy[t - 1], coeff: -dt, },
          { var_id: uy[t], coeff: -h2, },
        ],
        Eq,
        0.0,
      )
      m.add_constraint(
        "dvx" + t.to_string(),
        [
          { var_id: vx[t], coeff: 1.0, },
          { var_id: vx[t - 1], coeff: -1.0, },
          { var_id: ux[t], coeff: -dt, },
        ],
        Eq,
        0.0,
      )
      m.add_constraint(
        "dvy" + t.to_string(),
        [
          { var_id: vy[t], coeff: 1.0, },
          { var_id: vy[t - 1], coeff: -1.0, },
          { var_id: uy[t], coeff: -dt, },
        ],
        Eq,
        0.0,
      )
    }
  }
  m
}
