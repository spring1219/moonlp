///|
/// Rigorous Machine-Checkable KKT Optimality Certificate for Linear Programs (v0.6.0):
/// Independently verifies ALL KKT conditions for arbitrary variable bounds & senses:
///   0. Finiteness: no NaN / +/-Inf in primal, dual, or objective
///   1. Primal feasibility (constraints + variable bounds)
///   2. Dual sign feasibility (y_i <= 0 for Min/Le, y_i >= 0 for Min/Ge, reversed for Max)
///   3. Constraint complementary slackness: |y_i * (a_i^T x - b_i)| <= tol
///   4. KKT stationarity & bound complementary slackness on r = c - A^T y:
///        s_j = sense_sign * r_j must satisfy:
///          x_j > lb_j + tol => s_j <= tol
///          x_j < ub_j - tol => s_j >= -tol
///   5. Full strong duality including finite variable bound multipliers:
///        |z_primal - z_dual| <= tol * (1 + |z_primal|)

///|
#derive(Debug)
pub(all) struct OptimalityCertificate {
  primal_feasible : Bool
  dual_feasible : Bool
  max_primal_violation : Double
  max_dual_violation : Double
  max_complementarity_error : Double
  primal_objective : Double
  dual_objective : Double
  duality_gap : Double
  certified_optimal : Bool
}

///|
/// Generate a full, rigorous KKT optimality certificate for an LP solution.
pub fn Model::certify_lp_solution(
  self : Model,
  sol : Solution,
  tol~ : Double = 1.0e-6,
) -> OptimalityCertificate {
  let invalid : OptimalityCertificate = {
    primal_feasible: false,
    dual_feasible: false,
    max_primal_violation: 1.0e300,
    max_dual_violation: 1.0e300,
    max_complementarity_error: 1.0e300,
    primal_objective: 0.0,
    dual_objective: 0.0,
    duality_gap: 1.0e300,
    certified_optimal: false,
  }

  let n = self.vars.length()
  let m = self.constraints.length()
  if sol.status != Optimal ||
    sol.primal.length() != n ||
    sol.dual.length() != m ||
    !is_finite_double(sol.objective) {
    return invalid
  }

  let sense_sign = match self.sense {
    Minimize => 1.0
    Maximize => -1.0
  }

  // 1. Check finiteness & variable bounds, compute linear & total primal objective
  let mut max_pv = 0.0
  let mut p_lin = 0.0
  for j = 0; j < n; j = j + 1 {
    let x = sol.primal[j]
    if !is_finite_double(x) {
      return invalid
    }
    let v = self.vars[j]
    let term_obj = v.obj * x
    p_lin = p_lin + term_obj
    if !is_finite_double(term_obj) || !is_finite_double(p_lin) {
      return invalid
    }
    if x < v.lb {
      let viol = v.lb - x
      if viol > max_pv {
        max_pv = viol
      }
    }
    if x > v.ub {
      let viol = x - v.ub
      if viol > max_pv {
        max_pv = viol
      }
    }
  }
  let p_obj = p_lin + self.obj_constant
  if !is_finite_double(p_obj) ||
    (p_obj - sol.objective).abs() >
    tol * (1.0 + p_lin.abs() + self.obj_constant.abs()) {
    return invalid
  }

  // 2. Check constraints, per-row dual sign feasibility, local complementary slackness,
  //    and accumulate A^T y & local column scales |A|^T |y| (strictly isolated per row/col!)
  let at_y = Array::make(n, 0.0)
  let at_y_abs = Array::make(n, 0.0)
  let mut max_dv = 0.0
  let mut max_comp = 0.0
  let mut dual_ok = true
  let mut comp_ok = true
  let mut d_lin = 0.0

  for i = 0; i < m; i = i + 1 {
    let y = sol.dual[i]
    if !is_finite_double(y) {
      return invalid
    }
    let c = self.constraints[i]
    let mut lhs = 0.0
    for k = 0; k < c.terms.length(); k = k + 1 {
      let t = c.terms[k]
      let idx = t.var_id.index
      let p_term = t.coeff * sol.primal[idx]
      lhs = lhs + p_term
      let d_term = t.coeff * y
      at_y[idx] = at_y[idx] + d_term
      at_y_abs[idx] = at_y_abs[idx] + d_term.abs()
      if !is_finite_double(p_term) ||
        !is_finite_double(lhs) ||
        !is_finite_double(d_term) ||
        !is_finite_double(at_y[idx]) ||
        !is_finite_double(at_y_abs[idx]) {
        return invalid
      }
    }

    let slack = match c.relation {
      Le => c.rhs - lhs
      Ge => lhs - c.rhs
      Eq => -(lhs - c.rhs).abs()
    }
    if !is_finite_double(slack) {
      return invalid
    }
    if slack < 0.0 {
      let viol = -slack
      if viol > max_pv {
        max_pv = viol
      }
    }

    // Strict per-row dual sign feasibility (pure absolute `tol`, never inflated by unrelated vars!):
    // Minimize: Le => y <= 0, Ge => y >= 0, Eq => free
    // Maximize: Le => y >= 0, Ge => y <= 0, Eq => free
    let y_min = sense_sign * y
    match c.relation {
      Le =>
        if y_min > 0.0 {
          if y_min > max_dv {
            max_dv = y_min
          }
          if y_min > tol {
            dual_ok = false
          }
        }
      Ge =>
        if y_min < 0.0 {
          let dv = -y_min
          if dv > max_dv {
            max_dv = dv
          }
          if dv > tol {
            dual_ok = false
          }
        }
      Eq => ()
    }

    let abs_res = (lhs - c.rhs).abs()
    let comp = (y * (lhs - c.rhs)).abs()
    let d_term_rhs = y * c.rhs
    d_lin = d_lin + d_term_rhs
    if !is_finite_double(abs_res) ||
      !is_finite_double(comp) ||
      !is_finite_double(d_term_rhs) ||
      !is_finite_double(d_lin) {
      return invalid
    }
    if comp > max_comp {
      max_comp = comp
    }
    let min_pair = if y.abs() < abs_res { y.abs() } else { abs_res }
    if min_pair > tol || comp > tol * (1.0 + y.abs()) {
      comp_ok = false
    }
  }

  // 3. Per-variable KKT stationarity & bound multipliers (strictly local scale per variable!)
  for j = 0; j < n; j = j + 1 {
    let v = self.vars[j]
    let x = sol.primal[j]
    let s_j = sense_sign * (v.obj - at_y[j]) // reduced cost in minimization space
    let slope_scale = v.obj.abs() + at_y_abs[j]
    let local_scale = 1.0 + slope_scale
    if !is_finite_double(s_j) || !is_finite_double(local_scale) {
      return invalid
    }
    let has_lb = v.lb > -inf_threshold
    let has_ub = v.ub < inf_threshold

    // If s_j > 0, lower bound must be finite and active (x == lb)
    if s_j > 0.0 {
      if !has_lb {
        if s_j > max_dv {
          max_dv = s_j
        }
        if s_j > tol * slope_scale {
          dual_ok = false
        }
      } else {
        let gap_lb = x - v.lb
        if gap_lb > tol {
          if s_j > max_dv {
            max_dv = s_j
          }
          let bc = (s_j * gap_lb).abs()
          if !is_finite_double(bc) {
            return invalid
          }
          if bc > max_comp {
            max_comp = bc
          }
          if s_j > tol * local_scale {
            dual_ok = false
          }
          if bc > tol * local_scale {
            comp_ok = false
          }
        }
        let lb_contrib = sense_sign * s_j * v.lb
        d_lin = d_lin + lb_contrib
        if !is_finite_double(lb_contrib) || !is_finite_double(d_lin) {
          return invalid
        }
      }
    } else if s_j < 0.0 {
      // If s_j < 0, upper bound must be finite and active (x == ub)
      let mag = -s_j
      if !has_ub {
        if mag > max_dv {
          max_dv = mag
        }
        if mag > tol * slope_scale {
          dual_ok = false
        }
      } else {
        let gap_ub = v.ub - x
        if gap_ub > tol {
          if mag > max_dv {
            max_dv = mag
          }
          let bc = (mag * gap_ub).abs()
          if !is_finite_double(bc) {
            return invalid
          }
          if bc > max_comp {
            max_comp = bc
          }
          if mag > tol * local_scale {
            dual_ok = false
          }
          if bc > tol * local_scale {
            comp_ok = false
          }
        }
        let ub_contrib = sense_sign * mag * v.ub
        d_lin = d_lin - ub_contrib
        if !is_finite_double(ub_contrib) || !is_finite_double(d_lin) {
          return invalid
        }
      }
    }
  }

  let d_obj = d_lin + self.obj_constant
  let lin_scale = 1.0 +
    (if p_lin.abs() > d_lin.abs() { p_lin.abs() } else { d_lin.abs() })
  let gap = (p_lin - d_lin).abs()
  if !is_finite_double(d_obj) ||
    !is_finite_double(lin_scale) ||
    !is_finite_double(gap) ||
    !is_finite_double(max_pv) ||
    !is_finite_double(max_dv) ||
    !is_finite_double(max_comp) {
    return invalid
  }
  let primal_ok = max_pv <= tol
  let gap_ok = gap <= tol * lin_scale

  {
    primal_feasible: primal_ok,
    dual_feasible: dual_ok,
    max_primal_violation: max_pv,
    max_dual_violation: max_dv,
    max_complementarity_error: max_comp,
    primal_objective: p_obj,
    dual_objective: d_obj,
    duality_gap: gap,
    certified_optimal: primal_ok && dual_ok && comp_ok && gap_ok,
  }
}

///|
/// Generate a full, rigorous Quadratic KKT + Wolfe Strong Duality optimality
/// certificate for a convex QP solution (also subsuming pure LP when Q == 0).
/// Independently verifies:
///   0. Convexity (sense_sign * Q is positive semi-definite) & finiteness
///   1. Primal feasibility (constraints + variable bounds)
///   2. Per-row dual sign feasibility (y_i <= 0 for Min/Le, y_i >= 0 for Min/Ge)
///   3. Constraint complementary slackness: |y_i * (a_i^T x - b_i)| <= tol
///   4. Quadratic KKT stationarity & bound complementarity on
///        s_j = sense_sign * ((Q x)_j + c_j - (A^T y)_j)
///   5. Wolfe Strong Duality:
///        p_obj = 0.5 x^T Q x + c^T x + c_0
///        d_obj = -0.5 x^T Q x + y^T b + sense_sign * (z_l^T l - z_u^T u) + c_0
pub fn Model::certify_qp_solution(
  self : Model,
  sol : Solution,
  tol~ : Double = 1.0e-6,
) -> OptimalityCertificate {
  if !self.has_quad_obj() {
    return self.certify_lp_solution(sol, tol~)
  }

  let invalid : OptimalityCertificate = {
    primal_feasible: false,
    dual_feasible: false,
    max_primal_violation: 1.0e300,
    max_dual_violation: 1.0e300,
    max_complementarity_error: 1.0e300,
    primal_objective: 0.0,
    dual_objective: 0.0,
    duality_gap: 1.0e300,
    certified_optimal: false,
  }

  let n = self.vars.length()
  let m = self.constraints.length()
  if sol.status != Optimal ||
    sol.primal.length() != n ||
    sol.dual.length() != m ||
    !is_finite_double(sol.objective) {
    return invalid
  }

  let sense_sign = match self.sense {
    Minimize => 1.0
    Maximize => -1.0
  }

  // Verify convexity of Hessian in minimization space on non-fixed variables
  let q_orig = self.build_symmetric_hessian()
  let act_indices : Array[Int32] = []
  for j = 0; j < n; j = j + 1 {
    let v = self.vars[j]
    let fixed = v.lb > -inf_threshold && v.ub < inf_threshold && v.lb >= v.ub
    if !fixed {
      act_indices.push(j)
    }
  }
  let n_act = act_indices.length()
  let h_act : Array[Array[Double]] = []
  for a = 0; a < n_act; a = a + 1 {
    let row = Array::make(n_act, 0.0)
    for b = 0; b < n_act; b = b + 1 {
      row[b] = q_orig[act_indices[a]][act_indices[b]] * sense_sign
    }
    h_act.push(row)
  }
  if !is_symmetric_psd(h_act, n_act) {
    return invalid
  }

  // 1. Check finiteness, variable bounds, linear + quadratic primal objective, and Q x
  let mut max_pv = 0.0
  let mut p_lin = 0.0
  for j = 0; j < n; j = j + 1 {
    let x = sol.primal[j]
    if !is_finite_double(x) {
      return invalid
    }
    let v = self.vars[j]
    let term_obj = v.obj * x
    p_lin = p_lin + term_obj
    if !is_finite_double(term_obj) || !is_finite_double(p_lin) {
      return invalid
    }
    if x < v.lb {
      let viol = v.lb - x
      if viol > max_pv {
        max_pv = viol
      }
    }
    if x > v.ub {
      let viol = x - v.ub
      if viol > max_pv {
        max_pv = viol
      }
    }
  }

  let qx = Array::make(n, 0.0)
  let qx_abs = Array::make(n, 0.0)
  let mut x_q_x = 0.0
  for i = 0; i < n; i = i + 1 {
    let mut sum = 0.0
    let mut sum_abs = 0.0
    for j = 0; j < n; j = j + 1 {
      let qij = q_orig[i][j]
      if qij != 0.0 {
        let term = qij * sol.primal[j]
        sum = sum + term
        sum_abs = sum_abs + term.abs()
      }
    }
    if !is_finite_double(sum) || !is_finite_double(sum_abs) {
      return invalid
    }
    qx[i] = sum
    qx_abs[i] = sum_abs
    x_q_x = x_q_x + sol.primal[i] * sum
    if !is_finite_double(x_q_x) {
      return invalid
    }
  }

  let p_quad = 0.5 * x_q_x
  let p_var = p_lin + p_quad
  let p_obj = p_var + self.obj_constant
  if !is_finite_double(p_obj) ||
    (p_obj - sol.objective).abs() >
    tol * (1.0 + p_var.abs() + self.obj_constant.abs()) {
    return invalid
  }

  // 2. Constraints, dual sign feasibility, complementarity, and A^T y
  let at_y = Array::make(n, 0.0)
  let at_y_abs = Array::make(n, 0.0)
  let mut max_dv = 0.0
  let mut max_comp = 0.0
  let mut dual_ok = true
  let mut comp_ok = true
  let mut d_var = -p_quad

  for i = 0; i < m; i = i + 1 {
    let y = sol.dual[i]
    if !is_finite_double(y) {
      return invalid
    }
    let c = self.constraints[i]
    let mut lhs = 0.0
    for k = 0; k < c.terms.length(); k = k + 1 {
      let t = c.terms[k]
      let idx = t.var_id.index
      let p_term = t.coeff * sol.primal[idx]
      lhs = lhs + p_term
      let d_term = t.coeff * y
      at_y[idx] = at_y[idx] + d_term
      at_y_abs[idx] = at_y_abs[idx] + d_term.abs()
      if !is_finite_double(p_term) ||
        !is_finite_double(lhs) ||
        !is_finite_double(d_term) ||
        !is_finite_double(at_y[idx]) ||
        !is_finite_double(at_y_abs[idx]) {
        return invalid
      }
    }

    let slack = match c.relation {
      Le => c.rhs - lhs
      Ge => lhs - c.rhs
      Eq => -(lhs - c.rhs).abs()
    }
    if !is_finite_double(slack) {
      return invalid
    }
    if slack < 0.0 {
      let viol = -slack
      if viol > max_pv {
        max_pv = viol
      }
    }

    let y_min = sense_sign * y
    match c.relation {
      Le =>
        if y_min > 0.0 {
          if y_min > max_dv {
            max_dv = y_min
          }
          if y_min > tol {
            dual_ok = false
          }
        }
      Ge =>
        if y_min < 0.0 {
          let dv = -y_min
          if dv > max_dv {
            max_dv = dv
          }
          if dv > tol {
            dual_ok = false
          }
        }
      Eq => ()
    }

    let abs_res = (lhs - c.rhs).abs()
    let comp = (y * (lhs - c.rhs)).abs()
    let d_term_rhs = y * c.rhs
    d_var = d_var + d_term_rhs
    if !is_finite_double(abs_res) ||
      !is_finite_double(comp) ||
      !is_finite_double(d_term_rhs) ||
      !is_finite_double(d_var) {
      return invalid
    }
    if comp > max_comp {
      max_comp = comp
    }
    let min_pair = if y.abs() < abs_res { y.abs() } else { abs_res }
    if min_pair > tol || comp > tol * (1.0 + y.abs()) {
      comp_ok = false
    }
  }

  // 3. Quadratic KKT stationarity & bound multipliers
  for j = 0; j < n; j = j + 1 {
    let v = self.vars[j]
    let x = sol.primal[j]
    let s_j = sense_sign * (qx[j] + v.obj - at_y[j])
    let slope_scale = qx_abs[j] + v.obj.abs() + at_y_abs[j]
    let local_scale = 1.0 + slope_scale
    if !is_finite_double(s_j) || !is_finite_double(local_scale) {
      return invalid
    }
    let has_lb = v.lb > -inf_threshold
    let has_ub = v.ub < inf_threshold

    if s_j > 0.0 {
      if !has_lb {
        if s_j > max_dv {
          max_dv = s_j
        }
        if s_j > tol * local_scale {
          dual_ok = false
        }
      } else {
        let gap_lb = x - v.lb
        if gap_lb > tol {
          if s_j > max_dv {
            max_dv = s_j
          }
          let bc = (s_j * gap_lb).abs()
          if !is_finite_double(bc) {
            return invalid
          }
          if bc > max_comp {
            max_comp = bc
          }
          if s_j > tol * local_scale {
            dual_ok = false
          }
          if bc > tol * local_scale {
            comp_ok = false
          }
        }
        let lb_contrib = sense_sign * s_j * v.lb
        d_var = d_var + lb_contrib
        if !is_finite_double(lb_contrib) || !is_finite_double(d_var) {
          return invalid
        }
      }
    } else if s_j < 0.0 {
      let mag = -s_j
      if !has_ub {
        if mag > max_dv {
          max_dv = mag
        }
        if mag > tol * local_scale {
          dual_ok = false
        }
      } else {
        let gap_ub = v.ub - x
        if gap_ub > tol {
          if mag > max_dv {
            max_dv = mag
          }
          let bc = (mag * gap_ub).abs()
          if !is_finite_double(bc) {
            return invalid
          }
          if bc > max_comp {
            max_comp = bc
          }
          if mag > tol * local_scale {
            dual_ok = false
          }
          if bc > tol * local_scale {
            comp_ok = false
          }
        }
        let ub_contrib = sense_sign * mag * v.ub
        d_var = d_var - ub_contrib
        if !is_finite_double(ub_contrib) || !is_finite_double(d_var) {
          return invalid
        }
      }
    }
  }

  let d_obj = d_var + self.obj_constant
  let var_scale = 1.0 +
    (if p_var.abs() > d_var.abs() { p_var.abs() } else { d_var.abs() })
  let gap = (p_var - d_var).abs()
  if !is_finite_double(d_obj) ||
    !is_finite_double(var_scale) ||
    !is_finite_double(gap) ||
    !is_finite_double(max_pv) ||
    !is_finite_double(max_dv) ||
    !is_finite_double(max_comp) {
    return invalid
  }
  let primal_ok = max_pv <= tol
  let gap_ok = gap <= tol * var_scale

  {
    primal_feasible: primal_ok,
    dual_feasible: dual_ok,
    max_primal_violation: max_pv,
    max_dual_violation: max_dv,
    max_complementarity_error: max_comp,
    primal_objective: p_obj,
    dual_objective: d_obj,
    duality_gap: gap,
    certified_optimal: primal_ok && dual_ok && comp_ok && gap_ok,
  }
}

///|
/// Preserve genuine small shadow prices. When roundoff prevents certification,
/// try the zeroed multipliers and use them only if the original LP certifies.
fn clean_dual_roundoff(
  model : Model,
  sol : Solution,
  config : SolverConfig,
) -> Solution {
  if !sol.dual.any(y => y != 0.0 && y.abs() < 1.0e-15) ||
    model.certify_lp_solution(sol, tol=config.tol).certified_optimal {
    return sol
  }
  let dual = sol.dual.map(y => if y.abs() < 1.0e-15 { 0.0 } else { y })
  let cleaned = { ..sol, dual, }
  if model.certify_lp_solution(cleaned, tol=config.tol).certified_optimal {
    cleaned
  } else {
    sol
  }
}
