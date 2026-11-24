///|
/// Algebraic LP/MILP Model builder and solution verifier.

///|
#derive(Debug)
pub struct Model {
  mut sense : Sense
  vars : Array[VarSpec]
  constraints : Array[Constraint]
  mut obj_constant : Double
}

///|
pub fn Model::new(sense : Sense) -> Model {
  { sense, vars: [], constraints: [], obj_constant: 0.0, }
}

///|
pub fn Model::num_vars(self : Model) -> Int32 {
  self.vars.length()
}

///|
pub fn Model::num_constraints(self : Model) -> Int32 {
  self.constraints.length()
}

///|
pub fn Model::set_sense(self : Model, sense : Sense) -> Unit {
  self.sense = sense
}

///|
pub fn Model::set_obj_constant(self : Model, c : Double) -> Unit {
  self.obj_constant = c
}

///|
/// Add a variable with explicit bounds, kind, and objective coefficient.
pub fn Model::add_var(
  self : Model,
  name : String,
  lb~ : Double = 0.0,
  ub~ : Double = 1.0e300,
  kind~ : VarKind = Continuous,
  obj~ : Double = 0.0,
) -> VarId {
  let idx = self.vars.length()
  let eff_lb = match kind {
    Binary => if lb < 0.0 { 0.0 } else { lb }
    _ => lb
  }
  let eff_ub = match kind {
    Binary => if ub > 1.0 { 1.0 } else { ub }
    _ => ub
  }
  self.vars.push({ name, lb: eff_lb, ub: eff_ub, kind, obj, })
  { index: idx, }
}

///|
pub fn Model::add_continuous_var(
  self : Model,
  name : String,
  lb~ : Double = 0.0,
  ub~ : Double = 1.0e300,
  obj~ : Double = 0.0,
) -> VarId {
  self.add_var(name, lb~, ub~, kind=Continuous, obj~)
}

///|
pub fn Model::add_int_var(
  self : Model,
  name : String,
  lb~ : Double = 0.0,
  ub~ : Double = 1.0e300,
  obj~ : Double = 0.0,
) -> VarId {
  self.add_var(name, lb~, ub~, kind=Integer, obj~)
}

///|
pub fn Model::add_binary_var(
  self : Model,
  name : String,
  obj~ : Double = 0.0,
) -> VarId {
  self.add_var(name, lb=0.0, ub=1.0, kind=Binary, obj~)
}

///|
pub fn Model::set_var_obj(self : Model, v : VarId, coeff : Double) -> Unit {
  let old = self.vars[v.index]
  self.vars[v.index] = { ..old, obj: coeff, }
}

///|
pub fn Model::set_var_bounds(
  self : Model,
  v : VarId,
  lb : Double,
  ub : Double,
) -> Unit {
  let old = self.vars[v.index]
  self.vars[v.index] = { ..old, lb, ub, }
}

///|
pub fn Model::set_objective(self : Model, expr : LinearExpr) -> Unit {
  for i = 0; i < self.vars.length(); i = i + 1 {
    let old = self.vars[i]
    self.vars[i] = { ..old, obj: 0.0, }
  }
  for i = 0; i < expr.terms.length(); i = i + 1 {
    let t = expr.terms[i]
    let old = self.vars[t.var_id.index]
    self.vars[t.var_id.index] = { ..old, obj: old.obj + t.coeff, }
  }
  self.obj_constant = expr.constant
}

///|
pub fn Model::add_constraint(
  self : Model,
  name : String,
  terms : Array[Term],
  relation : Relation,
  rhs : Double,
) -> Unit {
  self.constraints.push({ name, terms, relation, rhs, })
}

///|
pub fn Model::add_expr_constraint(
  self : Model,
  name : String,
  expr : LinearExpr,
  relation : Relation,
  rhs : Double,
) -> Unit {
  self.constraints.push({
    name,
    terms: expr.terms,
    relation,
    rhs: rhs - expr.constant,
  })
}

///|
pub fn Model::has_integer_vars(self : Model) -> Bool {
  for i = 0; i < self.vars.length(); i = i + 1 {
    match self.vars[i].kind {
      Continuous => ()
      Integer | Binary => return true
    }
  }
  false
}

///|
/// Clone model with updated variable bounds (used internally by Branch-and-Bound).
pub fn Model::clone_with_bounds(
  self : Model,
  lbs : Array[Double],
  ubs : Array[Double],
) -> Model {
  let new_vars : Array[VarSpec] = []
  for i = 0; i < self.vars.length(); i = i + 1 {
    let v = self.vars[i]
    new_vars.push({ ..v, lb: lbs[i], ub: ubs[i], })
  }
  {
    sense: self.sense,
    vars: new_vars,
    constraints: self.constraints,
    obj_constant: self.obj_constant,
  }
}

///|
/// Independently verify primal feasibility, variable bounds, integrality, and objective value.
pub fn Model::verify_solution(
  self : Model,
  sol : Solution,
  tol~ : Double = 1.0e-5,
) -> Bool {
  if sol.status != Optimal {
    return false
  }
  if sol.primal.length() != self.vars.length() {
    return false
  }
  if !is_finite_double(sol.objective) {
    return false
  }
  let mut obj_val = self.obj_constant
  for i = 0; i < self.vars.length(); i = i + 1 {
    let x = sol.primal[i]
    if !is_finite_double(x) {
      return false
    }
    let v = self.vars[i]
    if x < v.lb - tol || x > v.ub + tol {
      return false
    }
    match v.kind {
      Continuous => ()
      Integer | Binary => {
        let rounded = d_round(x)
        if (x - rounded).abs() > tol {
          return false
        }
      }
    }
    obj_val = obj_val + v.obj * x
  }
  if !is_finite_double(obj_val) ||
    (obj_val - sol.objective).abs() > tol * (1.0 + obj_val.abs()) {
    return false
  }
  for c_idx = 0; c_idx < self.constraints.length(); c_idx = c_idx + 1 {
    let c = self.constraints[c_idx]
    let mut lhs = 0.0
    for j = 0; j < c.terms.length(); j = j + 1 {
      let t = c.terms[j]
      lhs = lhs + t.coeff * sol.primal[t.var_id.index]
    }
    if !is_finite_double(lhs) {
      return false
    }
    match c.relation {
      Le => if lhs > c.rhs + tol { return false }
      Ge => if lhs < c.rhs - tol { return false }
      Eq => if (lhs - c.rhs).abs() > tol { return false }
    }
  }
  true
}

///|
/// High-level entry point: automatically dispatches to LP Simplex or MILP Branch-and-Bound.
pub fn Model::solve(
  self : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  if self.has_integer_vars() {
    solve_milp(self, config~)
  } else {
    solve_lp(self, config~)
  }
}
