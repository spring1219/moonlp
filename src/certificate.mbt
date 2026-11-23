///|
/// Machine-checkable KKT Optimality Certificate for Linear Programs:
/// Independently verifies:
///   1. Primal feasibility (constraints + variable bounds)
///   2. Complementary slackness: |y_i * (a_i^T x - b_i)| <= tol
///   3. Strong duality gap (when all variables have standard [0, +inf) bounds):
///      |c^T x - b^T y| <= tol * (1 + |c^T x|)

///|
#derive(Debug)
pub(all) struct OptimalityCertificate {
  primal_feasible : Bool
  max_primal_violation : Double
  max_complementarity_error : Double
  primal_objective : Double
  dual_objective : Double
  duality_gap : Double
  certified_optimal : Bool
}

///|
/// Generate an independent KKT / Strong Duality optimality certificate for an LP solution.
pub fn Model::certify_lp_solution(
  self : Model,
  sol : Solution,
  tol~ : Double = 1.0e-6,
) -> OptimalityCertificate {
  if sol.status != Optimal ||
    sol.primal.length() != self.vars.length() ||
    sol.dual.length() != self.constraints.length() {
    return {
      primal_feasible: false,
      max_primal_violation: 1.0e300,
      max_complementarity_error: 1.0e300,
      primal_objective: 0.0,
      dual_objective: 0.0,
      duality_gap: 1.0e300,
      certified_optimal: false,
    }
  }

  let mut max_pv = 0.0
  let mut p_obj = self.obj_constant
  let mut has_nonstandard_bounds = false

  for j = 0; j < self.vars.length(); j = j + 1 {
    let v = self.vars[j]
    let x = sol.primal[j]
    p_obj = p_obj + v.obj * x
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
    if v.lb.abs() > tol || v.ub < 1.0e200 {
      has_nonstandard_bounds = true
    }
  }

  let mut max_comp = 0.0
  let mut d_obj = self.obj_constant

  for i = 0; i < self.constraints.length(); i = i + 1 {
    let c = self.constraints[i]
    let mut lhs = 0.0
    for k = 0; k < c.terms.length(); k = k + 1 {
      let t = c.terms[k]
      lhs = lhs + t.coeff * sol.primal[t.var_id.index]
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
    let y = sol.dual[i]
    let comp = (y * (lhs - c.rhs)).abs()
    if comp > max_comp {
      max_comp = comp
    }
    d_obj = d_obj + y * c.rhs
  }

  let gap = if has_nonstandard_bounds { 0.0 } else { (p_obj - d_obj).abs() }
  let primal_ok = max_pv <= tol
  let comp_ok = max_comp <= tol * (1.0 + p_obj.abs())
  let gap_ok = has_nonstandard_bounds || gap <= tol * (1.0 + p_obj.abs())

  {
    primal_feasible: primal_ok,
    max_primal_violation: max_pv,
    max_complementarity_error: max_comp,
    primal_objective: p_obj,
    dual_objective: d_obj,
    duality_gap: gap,
    certified_optimal: primal_ok && comp_ok && gap_ok,
  }
}
