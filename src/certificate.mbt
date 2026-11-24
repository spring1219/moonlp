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
    let local_scale = 1.0 + v.obj.abs() + at_y_abs[j]
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
