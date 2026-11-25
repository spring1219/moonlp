///|
/// Concrete Real-World Convex QP & Mixed-Integer QP (MIQP) Optimization Models (v1.0.0):
/// Demonstrates new capabilities strictly beyond linear LP/MILP solvers:
///   1. Markowitz Mean-Variance Portfolio Optimization (Continuous QP & Efficient Frontier)
///   2. Cardinality-Constrained Portfolio Selection (Mixed-Integer QP, MIQP)
///   3. Constrained Linear Model Predictive Control (MPC Trajectory Tracking QP)
///   4. Soft-Margin Support Vector Machine Dual Training (Rank-Deficient PSD QP)
///   5. Simplex-Constrained Least Squares / Spectral Unmixing (QP)

///|
/// Return the 6-asset expected return vector and strictly positive definite
/// annual covariance matrix shared by the Markowitz QP and Cardinality MIQP benchmarks.
fn markowitz6_data() -> (Array[Double], Array[Array[Double]]) {
  // Assets: [US_Eq, Intl_Eq, Corp_Bond, Gov_Bond, REIT, Gold]
  let mu = [0.11, 0.10, 0.06, 0.035, 0.085, 0.055]
  let sigma = [
    [0.0400, 0.0240, 0.0040, -0.0020, 0.0180, 0.0020],
    [0.0240, 0.0484, 0.0035, -0.0015, 0.0160, 0.0050],
    [0.0040, 0.0035, 0.0064, 0.0025, 0.0045, 0.0015],
    [-0.0020, -0.0015, 0.0025, 0.0025, -0.0010, 0.0010],
    [0.0180, 0.0160, 0.0045, -0.0010, 0.0324, 0.0030],
    [0.0020, 0.0050, 0.0015, 0.0010, 0.0030, 0.0225],
  ]
  (mu, sigma)
}

///|
/// 1. Markowitz Mean-Variance Portfolio Optimization (Convex QP):
///   min  0.5 * w^T Sigma w
///   s.t. sum(w_i) = 1, mu^T w >= target_return, 0 <= w_i <= max_weight
pub fn build_markowitz_portfolio_qp(
  target_return~ : Double = 0.08,
  max_weight~ : Double = 0.5,
) -> Model {
  let (mu, sigma) = markowitz6_data()
  let names = ["US_Eq", "Intl_Eq", "Corp_Bond", "Gov_Bond", "REIT", "Gold"]
  let m = Model::new(Minimize)
  let w : Array[VarId] = []
  let budget_terms : Array[Term] = []
  let ret_terms : Array[Term] = []
  for i = 0; i < 6; i = i + 1 {
    let wi = m.add_continuous_var(names[i], lb=0.0, ub=max_weight)
    w.push(wi)
    budget_terms.push({ var_id: wi, coeff: 1.0, })
    ret_terms.push({ var_id: wi, coeff: mu[i], })
  }
  m.set_quad_obj_matrix(w, sigma, half=true)
  m.add_constraint("budget", budget_terms, Eq, 1.0)
  m.add_constraint("min_return", ret_terms, Ge, target_return)
  m
}

///|
/// 2. Cardinality-Constrained Markowitz Portfolio Selection (MIQP):
///   min  0.5 * w^T Sigma w
///   s.t. sum(w_i) = 1, mu^T w >= target_return
///        min_weight * z_i <= w_i <= max_weight * z_i,  z_i in {0, 1}
///        sum(z_i) <= max_assets
/// Notice: 12x12 Hessian has rank 6 (zero block on binary indicators z_i).
pub fn build_cardinality_portfolio_miqp(
  max_assets~ : Int32 = 3,
  min_weight~ : Double = 0.1,
  max_weight~ : Double = 0.5,
  target_return~ : Double = 0.08,
) -> Model {
  let (mu, sigma) = markowitz6_data()
  let names = ["US_Eq", "Intl_Eq", "Corp_Bond", "Gov_Bond", "REIT", "Gold"]
  let m = Model::new(Minimize)
  let w : Array[VarId] = []
  let z : Array[VarId] = []
  let budget_terms : Array[Term] = []
  let ret_terms : Array[Term] = []
  let card_terms : Array[Term] = []
  for i = 0; i < 6; i = i + 1 {
    let wi = m.add_continuous_var("w_" + names[i], lb=0.0, ub=max_weight)
    w.push(wi)
    budget_terms.push({ var_id: wi, coeff: 1.0, })
    ret_terms.push({ var_id: wi, coeff: mu[i], })
  }
  for i = 0; i < 6; i = i + 1 {
    let zi = m.add_binary_var("z_" + names[i])
    z.push(zi)
    card_terms.push({ var_id: zi, coeff: 1.0, })
    m.add_constraint(
      "link_ub_" + i.to_string(),
      [{ var_id: w[i], coeff: 1.0, }, { var_id: zi, coeff: -max_weight, }],
      Le,
      0.0,
    )
    m.add_constraint(
      "link_lb_" + i.to_string(),
      [{ var_id: w[i], coeff: 1.0, }, { var_id: zi, coeff: -min_weight, }],
      Ge,
      0.0,
    )
  }
  m.set_quad_obj_matrix(w, sigma, half=true)
  m.add_constraint("budget", budget_terms, Eq, 1.0)
  m.add_constraint("min_return", ret_terms, Ge, target_return)
  m.add_constraint("cardinality", card_terms, Le, max_assets.to_double())
  m
}

///|
/// 3. Constrained Linear Model Predictive Control (MPC, Convex QP):
/// Discrete double-integrator plant (dt = 0.5):
///   p_{t+1} = p_t + 0.5 * v_t + 0.125 * u_t
///   v_{t+1} = v_t + 0.5 * u_t
/// Minimize 0.5 * sum_{t=1..N} (10 p_t^2 + 2 v_t^2) + 0.5 * sum_{t=0..N-1} 0.5 u_t^2
/// s.t. |u_t| <= u_max, |v_t| <= v_max, starting from (x0_pos, x0_vel).
pub fn build_linear_mpc_qp(
  horizon~ : Int32 = 8,
  x0_pos~ : Double = 2.0,
  x0_vel~ : Double = 0.0,
  u_max~ : Double = 1.0,
  v_max~ : Double = 1.5,
) -> Model {
  let dt = 0.5
  let half_dt2 = 0.5 * dt * dt
  let m = Model::new(Minimize)
  let p_vars : Array[VarId] = []
  let v_vars : Array[VarId] = []
  let u_vars : Array[VarId] = []

  // States p_t, v_t for t = 1..horizon
  for t = 1; t <= horizon; t = t + 1 {
    let pt = m.add_continuous_var(
      "p_" + t.to_string(),
      lb=-1.0e300,
      ub=1.0e300,
    )
    let vt = m.add_continuous_var("v_" + t.to_string(), lb=-v_max, ub=v_max)
    p_vars.push(pt)
    v_vars.push(vt)
    let q_p = if t == horizon { 20.0 } else { 10.0 }
    let q_v = if t == horizon { 4.0 } else { 2.0 }
    m.add_quad_obj_term(pt, pt, 0.5 * q_p)
    m.add_quad_obj_term(vt, vt, 0.5 * q_v)
  }

  // Control inputs u_t for t = 0..horizon-1
  for t = 0; t < horizon; t = t + 1 {
    let ut = m.add_continuous_var("u_" + t.to_string(), lb=-u_max, ub=u_max)
    u_vars.push(ut)
    m.add_quad_obj_term(ut, ut, 0.5 * 0.5)
  }

  // Dynamics constraints
  for t = 0; t < horizon; t = t + 1 {
    if t == 0 {
      // p_1 - 0.125 u_0 = x0_pos + 0.5 * x0_vel
      m.add_constraint(
        "dyn_p_0",
        [
          { var_id: p_vars[0], coeff: 1.0, },
          { var_id: u_vars[0], coeff: -half_dt2, },
        ],
        Eq,
        x0_pos + dt * x0_vel,
      )
      // v_1 - 0.5 u_0 = x0_vel
      m.add_constraint(
        "dyn_v_0",
        [{ var_id: v_vars[0], coeff: 1.0, }, { var_id: u_vars[0], coeff: -dt, }],
        Eq,
        x0_vel,
      )
    } else {
      m.add_constraint(
        "dyn_p_" + t.to_string(),
        [
          { var_id: p_vars[t], coeff: 1.0, },
          { var_id: p_vars[t - 1], coeff: -1.0, },
          { var_id: v_vars[t - 1], coeff: -dt, },
          { var_id: u_vars[t], coeff: -half_dt2, },
        ],
        Eq,
        0.0,
      )
      m.add_constraint(
        "dyn_v_" + t.to_string(),
        [
          { var_id: v_vars[t], coeff: 1.0, },
          { var_id: v_vars[t - 1], coeff: -1.0, },
          { var_id: u_vars[t], coeff: -dt, },
        ],
        Eq,
        0.0,
      )
    }
  }
  m
}

///|
/// 4. Soft-Margin Support Vector Machine (SVM) Dual QP:
/// 8 labeled 2D points (y_i in {+1, -1}) with rank-2 positive semi-definite Gram matrix:
///   min  0.5 * sum_{i,j} (y_i y_j x_i^T x_j) alpha_i alpha_j - sum_i alpha_i
///   s.t. sum_i y_i alpha_i = 0,  0 <= alpha_i <= c_reg
pub fn build_svm_dual_qp(c_reg~ : Double = 10.0) -> Model {
  let xs = [
    (2.0, 2.5),
    (2.5, 1.8),
    (3.0, 3.0),
    (1.6, 2.2),
    (-2.0, -2.2),
    (-1.5, -2.6),
    (-2.8, -1.8),
    (-1.2, -1.5),
  ]
  let ys = [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0]
  let n = xs.length()
  let m = Model::new(Minimize)
  let alphas : Array[VarId] = []
  let eq_terms : Array[Term] = []
  for i = 0; i < n; i = i + 1 {
    let ai = m.add_continuous_var(
      "a_" + i.to_string(),
      lb=0.0,
      ub=c_reg,
      obj=-1.0,
    )
    alphas.push(ai)
    eq_terms.push({ var_id: ai, coeff: ys[i], })
  }
  let gram : Array[Array[Double]] = []
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    let (xi1, xi2) = xs[i]
    for j = 0; j < n; j = j + 1 {
      let (xj1, xj2) = xs[j]
      row[j] = ys[i] * ys[j] * (xi1 * xj1 + xi2 * xj2)
    }
    gram.push(row)
  }
  m.set_quad_obj_matrix(alphas, gram, half=true)
  m.add_constraint("sum_y_alpha", eq_terms, Eq, 0.0)
  m
}

///|
/// 5. Simplex-Constrained Least Squares / Spectral Unmixing QP:
///   min 0.5 * ||A x - b||_2^2   s.t.  sum(x_j) = 1, x_j >= 0
pub fn build_constrained_lsq_qp() -> Model {
  // 4 measurements, 3 endmembers (true mixture roughly [0.5, 0.3, 0.2] + noise)
  let a_mat = [
    [1.0, 0.2, 0.1],
    [0.1, 1.0, 0.3],
    [0.2, 0.1, 1.0],
    [0.8, 0.6, 0.4],
  ]
  let b_obs = [0.60, 0.42, 0.34, 0.65]
  let m_rows = a_mat.length()
  let n_cols = 3
  let ata : Array[Array[Double]] = []
  for _i = 0; _i < n_cols; _i = _i + 1 {
    ata.push(Array::make(n_cols, 0.0))
  }
  let atb = Array::make(n_cols, 0.0)
  let mut btb = 0.0
  for r = 0; r < m_rows; r = r + 1 {
    btb = btb + b_obs[r] * b_obs[r]
    for j = 0; j < n_cols; j = j + 1 {
      atb[j] = atb[j] + a_mat[r][j] * b_obs[r]
      for k = 0; k < n_cols; k = k + 1 {
        ata[j][k] = ata[j][k] + a_mat[r][j] * a_mat[r][k]
      }
    }
  }
  let model = Model::new(Minimize)
  let x : Array[VarId] = []
  let simplex_terms : Array[Term] = []
  for j = 0; j < n_cols; j = j + 1 {
    let xj = model.add_continuous_var(
      "x_" + j.to_string(),
      lb=0.0,
      obj=-atb[j],
    )
    x.push(xj)
    simplex_terms.push({ var_id: xj, coeff: 1.0, })
  }
  model.set_quad_obj_matrix(x, ata, half=true)
  model.set_obj_constant(0.5 * btb)
  model.add_constraint("simplex", simplex_terms, Eq, 1.0)
  model
}
