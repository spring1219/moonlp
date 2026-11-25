///|
/// Post-optimality sensitivity analysis for Linear Programs:
/// - Constraint binding status, slack, shadow prices, and local RHS perturbation checks
/// - Non-basic variable reduced-cost analysis

///|
#derive(Debug)
pub(all) struct ConstraintSensitivity {
  name : String
  relation : Relation
  rhs : Double
  lhs : Double
  slack : Double
  binding : Bool
  shadow_price : Double
}

///|
#derive(Debug)
pub(all) struct VarSensitivity {
  name : String
  value : Double
  reduced_cost : Double
  at_lower_bound : Bool
  at_upper_bound : Bool
}

///|
/// Analyze constraint slack, binding status, and shadow prices at an LP optimum.
pub fn Model::analyze_constraints(
  self : Model,
  sol : Solution,
  tol~ : Double = 1.0e-6,
) -> Array[ConstraintSensitivity] {
  let out : Array[ConstraintSensitivity] = []
  if sol.status != Optimal || sol.primal.length() != self.vars.length() {
    return out
  }
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
      Eq => (c.rhs - lhs).abs()
    }
    let binding = slack.abs() <= tol
    let sp = if i < sol.dual.length() { sol.dual[i] } else { 0.0 }
    out.push({
      name: c.name,
      relation: c.relation,
      rhs: c.rhs,
      lhs,
      slack,
      binding,
      shadow_price: sp,
    })
  }
  out
}

///|
/// Analyze variable bound activity and reduced costs at an LP optimum.
pub fn Model::analyze_variables(
  self : Model,
  sol : Solution,
  tol~ : Double = 1.0e-6,
) -> Array[VarSensitivity] {
  let out : Array[VarSensitivity] = []
  if sol.status != Optimal || sol.primal.length() != self.vars.length() {
    return out
  }
  for j = 0; j < self.vars.length(); j = j + 1 {
    let v = self.vars[j]
    let x = sol.primal[j]
    let rc = if j < sol.reduced_costs.length() {
      sol.reduced_costs[j]
    } else {
      0.0
    }
    out.push({
      name: v.name,
      value: x,
      reduced_cost: rc,
      at_lower_bound: (x - v.lb).abs() <= tol,
      at_upper_bound: (x - v.ub).abs() <= tol,
    })
  }
  out
}

///|
/// Perturb constraint `constraint_idx`'s RHS by `delta` and re-solve to measure exact objective change.
pub fn Model::rhs_perturbation_delta(
  self : Model,
  base_sol : Solution,
  constraint_idx : Int32,
  delta : Double,
  config~ : SolverConfig = SolverConfig::default(),
) -> Double? {
  if base_sol.status != Optimal ||
    constraint_idx < 0 ||
    constraint_idx >= self.constraints.length() {
    return None
  }
  let new_constraints : Array[Constraint] = []
  for i = 0; i < self.constraints.length(); i = i + 1 {
    let c = self.constraints[i]
    if i == constraint_idx {
      new_constraints.push({ ..c, rhs: c.rhs + delta, })
    } else {
      new_constraints.push(c)
    }
  }
  let perturbed : Model = {
    sense: self.sense,
    vars: self.vars,
    constraints: new_constraints,
    obj_constant: self.obj_constant,
    quad_terms: self.quad_terms,
  }
  let new_sol = perturbed.solve(config~)
  if new_sol.status != Optimal {
    return None
  }
  Some(new_sol.objective - base_sol.objective)
}
