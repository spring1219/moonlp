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
  let mut max_c = 0.0
  for j = 0; j < n; j = j + 1 {
    let x = sol.primal[j]
    if !is_finite_double(x) {
      return invalid
    }
    let v = self.vars[j]
    if v.obj.abs() > max_c {
      max_c = v.obj.abs()
    }
    p_lin = p_lin + v.obj * x
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

  // 2. Check constraints, dual sign feasibility, complementary slackness, accumulate A^T y
  let at_y = Array::make(n, 0.0)
  let mut max_dv = 0.0
  let mut max_comp = 0.0
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
      lhs = lhs + t.coeff * sol.primal[t.var_id.index]
      at_y[t.var_id.index] = at_y[t.var_id.index] + t.coeff * y
    }
    if !is_finite_double(lhs) {
      return invalid
    }

    let slack = match c.relation {
      Le => c.rhs - lhs
      Ge => lhs - c.rhs
      Eq => -(lhs - c.rhs).abs()
    }
    if slack < 0.0 {
      let viol = -slack
      if viol > max_pv {
        max_pv = viol
      }
    }

    // Dual sign feasibility:
    // Minimize: Le => y <= 0, Ge => y >= 0, Eq => free
    // Maximize: Le => y >= 0, Ge => y <= 0, Eq => free
    let y_min = sense_sign * y
    match c.relation {
      Le =>
        if y_min > 0.0 && y_min > max_dv {
          max_dv = y_min
        }
      Ge =>
        if y_min < 0.0 && -y_min > max_dv {
          max_dv = -y_min
        }
      Eq => ()
    }

    let comp = (y * (lhs - c.rhs)).abs()
    if comp > max_comp {
      max_comp = comp
    }
    d_lin = d_lin + y * c.rhs
  }

  // 3. KKT stationarity & bound multipliers (exact general strong duality, invariant to obj_constant!)
  for j = 0; j < n; j = j + 1 {
    let v = self.vars[j]
    let x = sol.primal[j]
    let s_j = sense_sign * (v.obj - at_y[j]) // reduced cost in minimization space
    let has_lb = v.lb > -inf_threshold
    let has_ub = v.ub < inf_threshold

    // If s_j > 0, lower bound must be finite and active (x == lb)
    if s_j > 0.0 {
      if !has_lb {
        if s_j > max_dv {
          max_dv = s_j
        }
      } else {
        let gap_lb = x - v.lb
        if gap_lb > tol {
          if s_j > max_dv {
            max_dv = s_j
          }
          let bc = (s_j * gap_lb).abs()
          if bc > max_comp {
            max_comp = bc
          }
        }
        d_lin = d_lin + sense_sign * s_j * v.lb
      }
    } else if s_j < 0.0 {
      // If s_j < 0, upper bound must be finite and active (x == ub)
      let mag = -s_j
      if !has_ub {
        if mag > max_dv {
          max_dv = mag
        }
      } else {
        let gap_ub = v.ub - x
        if gap_ub > tol {
          if mag > max_dv {
            max_dv = mag
          }
          let bc = (mag * gap_ub).abs()
          if bc > max_comp {
            max_comp = bc
          }
        }
        d_lin = d_lin - sense_sign * mag * v.ub
      }
    }
  }

  // Scale-separated, obj_constant-invariant KKT checks:
  // - Primal feasibility uses pure absolute tolerance `tol` (never inflated by obj_constant)
  // - Dual & complementarity use objective-coefficient scale `1 + max|c_j|`
  // - Strong duality compares linear objective parts `p_lin` vs `d_lin` directly
  let d_obj = d_lin + self.obj_constant
  let c_scale = 1.0 + max_c
  let lin_scale = 1.0 +
    (if p_lin.abs() > d_lin.abs() { p_lin.abs() } else { d_lin.abs() })
  let gap = (p_lin - d_lin).abs()
  let primal_ok = max_pv <= tol
  let dual_ok = max_dv <= tol * c_scale
  let comp_ok = max_comp <= tol * c_scale
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
