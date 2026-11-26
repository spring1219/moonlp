///|
/// Algebraic LP/MILP/QP/MIQP/SOCP/MISOCP Model builder and solution verifier.

///|
#derive(Debug)
pub struct Model {
  mut sense : Sense
  vars : Array[VarSpec]
  constraints : Array[Constraint]
  mut obj_constant : Double
  quad_terms : Array[QuadTerm]
  soc_constraints : Array[SocConstraint]
}

///|
pub fn Model::new(sense : Sense) -> Model {
  {
    sense,
    vars: [],
    constraints: [],
    obj_constant: 0.0,
    quad_terms: [],
    soc_constraints: [],
  }
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
/// Add a quadratic objective term `coeff * x_{var1} * x_{var2}` directly to the objective.
/// Duplicate pairs are accumulated symmetrically when assembling the Hessian.
pub fn Model::add_quad_obj_term(
  self : Model,
  var1 : VarId,
  var2 : VarId,
  coeff : Double,
) -> Unit {
  if coeff != 0.0 {
    self.quad_terms.push({ var1, var2, coeff, })
  }
}

///|
/// Add a symmetric quadratic form to the objective over `vars`:
///   - `half = true` (default): adds `0.5 * sum_{i,j} q_mat[i][j] * x_i * x_j`
///   - `half = false`:          adds `sum_{i,j} q_mat[i][j] * x_i * x_j`
pub fn Model::set_quad_obj_matrix(
  self : Model,
  vars : Array[VarId],
  q_mat : Array[Array[Double]],
  half~ : Bool = true,
) -> Unit {
  self.quad_terms.clear()
  let k = vars.length()
  let factor = if half { 0.5 } else { 1.0 }
  for i = 0; i < k; i = i + 1 {
    if i >= q_mat.length() {
      break
    }
    let row = q_mat[i]
    for j = i; j < k; j = j + 1 {
      if j >= row.length() {
        break
      }
      if i == j {
        let c = factor * row[i]
        if c != 0.0 {
          self.quad_terms.push({ var1: vars[i], var2: vars[i], coeff: c, })
        }
      } else {
        let qji = if j < q_mat.length() && i < q_mat[j].length() {
          q_mat[j][i]
        } else {
          0.0
        }
        let c = factor * (row[j] + qji)
        if c != 0.0 {
          self.quad_terms.push({ var1: vars[i], var2: vars[j], coeff: c, })
        }
      }
    }
  }
}

///|
pub fn Model::clear_quad_obj(self : Model) -> Unit {
  self.quad_terms.clear()
}

///|
/// Check whether the model has any nonzero aggregated quadratic objective term.
pub fn Model::has_quad_obj(self : Model) -> Bool {
  if self.quad_terms.length() == 0 {
    return false
  }
  let n = self.vars.length()
  let q = self.build_symmetric_hessian()
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if q[i][j] != 0.0 {
        return true
      }
    }
  }
  false
}

///|
/// Assemble the symmetric n x n Hessian matrix Q such that the quadratic part of
/// the objective is exactly `0.5 * x^T Q x`.
fn Model::build_symmetric_hessian(self : Model) -> Array[Array[Double]] {
  let n = self.vars.length()
  let q : Array[Array[Double]] = []
  for _i = 0; _i < n; _i = _i + 1 {
    q.push(Array::make(n, 0.0))
  }
  for k = 0; k < self.quad_terms.length(); k = k + 1 {
    let qt = self.quad_terms[k]
    let i = qt.var1.index
    let j = qt.var2.index
    if i >= 0 && i < n && j >= 0 && j < n {
      if i == j {
        q[i][i] = q[i][i] + 2.0 * qt.coeff
      } else {
        q[i][j] = q[i][j] + qt.coeff
        q[j][i] = q[j][i] + qt.coeff
      }
    }
  }
  q
}

///|
/// Evaluate the exact objective value `obj_constant + c^T x + sum q_k x_i x_j`.
pub fn Model::eval_objective(self : Model, x : Array[Double]) -> Double {
  let mut val = self.obj_constant
  let n = self.vars.length()
  for i = 0; i < n; i = i + 1 {
    val = val + self.vars[i].obj * x[i]
  }
  for k = 0; k < self.quad_terms.length(); k = k + 1 {
    let qt = self.quad_terms[k]
    val = val + qt.coeff * x[qt.var1.index] * x[qt.var2.index]
  }
  val
}

///|
/// Pivoted Cholesky factorization of a symmetric positive semi-definite matrix P:
/// returns `Some(rows)` such that `P = sum_k rows[k] * rows[k]^T` (with rank <= n),
/// or `None` if P has a materially negative eigenvalue (non-convex).
fn psd_factor_rows(
  p_mat : Array[Array[Double]],
  n : Int32,
) -> Array[Array[Double]]? {
  if n <= 0 {
    return Some([])
  }
  let a : Array[Array[Double]] = []
  let mut max_abs = 0.0
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      let pij = if i < p_mat.length() && j < p_mat[i].length() {
        p_mat[i][j]
      } else {
        0.0
      }
      let pji = if j < p_mat.length() && i < p_mat[j].length() {
        p_mat[j][i]
      } else {
        0.0
      }
      let sym = 0.5 * (pij + pji)
      if !is_finite_double(sym) {
        return None
      }
      row[j] = sym
      if sym.abs() > max_abs {
        max_abs = sym.abs()
      }
    }
    a.push(row)
  }
  if max_abs == 0.0 {
    return Some([])
  }
  let neg_tol = 1.0e-9 * (1.0 + max_abs)
  let zero_tol = 1.0e-13 * (1.0 + max_abs)
  let rows : Array[Array[Double]] = []
  for _step = 0; _step < n; _step = _step + 1 {
    let mut piv = -1
    let mut best_d = -1.0e308
    for i = 0; i < n; i = i + 1 {
      let d = a[i][i]
      if d < -neg_tol {
        return None
      }
      if d > best_d {
        best_d = d
        piv = i
      }
    }
    if piv == -1 || best_d <= zero_tol {
      // Verify no large residual off-diagonals remain (which would imply indefiniteness)
      for i = 0; i < n; i = i + 1 {
        for j = i + 1; j < n; j = j + 1 {
          if a[i][j].abs() > neg_tol {
            return None
          }
        }
      }
      break
    }
    let root = d_sqrt(best_d)
    let r = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      r[j] = a[piv][j] / root
    }
    rows.push(r)
    for i = 0; i < n; i = i + 1 {
      let ri = r[i]
      if ri != 0.0 {
        let row_i = a[i]
        for j = 0; j < n; j = j + 1 {
          row_i[j] = row_i[j] - ri * r[j]
        }
      }
    }
  }
  Some(rows)
}

///|
/// Evaluate an affine LinearExpr at primal point `x`.
fn eval_linear_expr(expr : LinearExpr, x : Array[Double]) -> Double {
  let mut v = expr.constant
  for k = 0; k < expr.terms.length(); k = k + 1 {
    let t = expr.terms[k]
    v = v + t.coeff * x[t.var_id.index]
  }
  v
}

///|
/// Add a Second-Order Cone (Lorentz) constraint:
///   || [norm_exprs[0](x), ..., norm_exprs[k-1](x)] ||_2 <= rhs_expr(x)
pub fn Model::add_soc_constraint(
  self : Model,
  name : String,
  norm_exprs : Array[LinearExpr],
  rhs_expr : LinearExpr,
) -> Unit {
  self.soc_constraints.push({ name, norm_exprs, rhs_expr, })
}

///|
/// Add a convex quadratic constraint (QCQP) over `vars`:
///   - `half = true` (default): `0.5 * x_sub^T P x_sub + lin_expr(x) <= rhs`
///   - `half = false`:          `x_sub^T P x_sub + lin_expr(x) <= rhs`
/// Automatically factors PSD `P` via pivoted Cholesky and lifts into an exact
/// Second-Order Cone constraint. Returns `false` if `P` is not positive semi-definite.
pub fn Model::add_quad_constraint(
  self : Model,
  name : String,
  vars : Array[VarId],
  p_mat : Array[Array[Double]],
  lin_expr : LinearExpr,
  rhs : Double,
  half~ : Bool = true,
) -> Bool {
  let k = vars.length()
  let scaled_p : Array[Array[Double]] = []
  let factor = if half { 1.0 } else { 2.0 }
  for i = 0; i < k; i = i + 1 {
    let row = Array::make(k, 0.0)
    for j = 0; j < k; j = j + 1 {
      let pij = if i < p_mat.length() && j < p_mat[i].length() {
        p_mat[i][j]
      } else {
        0.0
      }
      row[j] = factor * pij
    }
    scaled_p.push(row)
  }
  match psd_factor_rows(scaled_p, k) {
    None => false
    Some(frows) => {
      if frows.length() == 0 {
        self.add_expr_constraint(name, lin_expr, Le, rhs)
        return true
      }
      // 0.5 * ||L^T x||^2 + u <= 0  <=>  || [L^T x, u + 0.5] ||_2 <= 0.5 - u
      // where u(x) = lin_expr(x) - rhs
      let norm_exprs : Array[LinearExpr] = []
      for r = 0; r < frows.length(); r = r + 1 {
        let e = LinearExpr::zero()
        let fr = frows[r]
        for j = 0; j < k; j = j + 1 {
          if fr[j] != 0.0 {
            let _ = e.add_term(vars[j], fr[j])
          }
        }
        norm_exprs.push(e)
      }
      let u_plus_half = LinearExpr::zero()
        .add_expr(lin_expr)
        .add_const(0.5 - rhs)
      norm_exprs.push(u_plus_half)
      let half_minus_u = lin_expr.scale(-1.0).add_const(0.5 + rhs)
      self.add_soc_constraint(name, norm_exprs, half_minus_u)
      true
    }
  }
}

///|
pub fn Model::has_soc_constraints(self : Model) -> Bool {
  self.soc_constraints.length() > 0
}

///|
pub fn Model::num_soc_constraints(self : Model) -> Int32 {
  self.soc_constraints.length()
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
    quad_terms: self.quad_terms,
    soc_constraints: self.soc_constraints,
  }
}

///|
/// Independently verify primal feasibility, variable bounds, integrality, and objective value
/// (supporting both linear and quadratic objectives).
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
    let obj_term = v.obj * x
    obj_val = obj_val + obj_term
    if !is_finite_double(obj_term) || !is_finite_double(obj_val) {
      return false
    }
  }
  for k = 0; k < self.quad_terms.length(); k = k + 1 {
    let qt = self.quad_terms[k]
    let q_term = qt.coeff * sol.primal[qt.var1.index] * sol.primal[qt.var2.index]
    obj_val = obj_val + q_term
    if !is_finite_double(q_term) || !is_finite_double(obj_val) {
      return false
    }
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
      let p_term = t.coeff * sol.primal[t.var_id.index]
      lhs = lhs + p_term
      if !is_finite_double(p_term) || !is_finite_double(lhs) {
        return false
      }
    }
    match c.relation {
      Le => if lhs > c.rhs + tol { return false }
      Ge => if lhs < c.rhs - tol { return false }
      Eq => if (lhs - c.rhs).abs() > tol { return false }
    }
  }
  // Verify Second-Order Cone constraints: ||norm_exprs(x)||_2 <= rhs_expr(x)
  for s_idx = 0; s_idx < self.soc_constraints.length(); s_idx = s_idx + 1 {
    let sc = self.soc_constraints[s_idx]
    let rhs_v = eval_linear_expr(sc.rhs_expr, sol.primal)
    if !is_finite_double(rhs_v) {
      return false
    }
    let mut sum_sq = 0.0
    for k = 0; k < sc.norm_exprs.length(); k = k + 1 {
      let ek = eval_linear_expr(sc.norm_exprs[k], sol.primal)
      if !is_finite_double(ek) {
        return false
      }
      sum_sq = sum_sq + ek * ek
    }
    let nrm = d_sqrt(sum_sq)
    if !is_finite_double(nrm) || rhs_v < -tol || nrm > rhs_v + tol {
      return false
    }
  }
  true
}

///|
/// Compute the exact maximum primal bound/constraint violation of a continuous solution
/// (ignoring integrality), or 1.0e300 if non-Optimal / non-finite / inconsistent objective.
fn Model::max_lp_violation(
  self : Model,
  sol : Solution,
  obj_tol~ : Double = 1.0e-5,
) -> Double {
  if sol.status != Optimal ||
    sol.primal.length() != self.vars.length() ||
    !is_finite_double(sol.objective) {
    return 1.0e300
  }
  let mut max_v = 0.0
  let mut obj_val = self.obj_constant
  for i = 0; i < self.vars.length(); i = i + 1 {
    let x = sol.primal[i]
    if !is_finite_double(x) {
      return 1.0e300
    }
    let v = self.vars[i]
    if x < v.lb {
      let d = v.lb - x
      if d > max_v {
        max_v = d
      }
    }
    if x > v.ub {
      let d = x - v.ub
      if d > max_v {
        max_v = d
      }
    }
    let term = v.obj * x
    obj_val = obj_val + term
    if !is_finite_double(term) || !is_finite_double(obj_val) {
      return 1.0e300
    }
  }
  for k = 0; k < self.quad_terms.length(); k = k + 1 {
    let qt = self.quad_terms[k]
    let q_term = qt.coeff * sol.primal[qt.var1.index] * sol.primal[qt.var2.index]
    obj_val = obj_val + q_term
    if !is_finite_double(q_term) || !is_finite_double(obj_val) {
      return 1.0e300
    }
  }
  if (obj_val - sol.objective).abs() >
    obj_tol * (1.0 + obj_val.abs() + self.obj_constant.abs()) {
    return 1.0e300
  }
  for c_idx = 0; c_idx < self.constraints.length(); c_idx = c_idx + 1 {
    let c = self.constraints[c_idx]
    let mut lhs = 0.0
    for j = 0; j < c.terms.length(); j = j + 1 {
      let t = c.terms[j]
      let term = t.coeff * sol.primal[t.var_id.index]
      lhs = lhs + term
      if !is_finite_double(term) || !is_finite_double(lhs) {
        return 1.0e300
      }
    }
    let viol = match c.relation {
      Le => if lhs > c.rhs { lhs - c.rhs } else { 0.0 }
      Ge => if lhs < c.rhs { c.rhs - lhs } else { 0.0 }
      Eq => (lhs - c.rhs).abs()
    }
    if !is_finite_double(viol) {
      return 1.0e300
    }
    if viol > max_v {
      max_v = viol
    }
  }
  for s_idx = 0; s_idx < self.soc_constraints.length(); s_idx = s_idx + 1 {
    let sc = self.soc_constraints[s_idx]
    let rhs_v = eval_linear_expr(sc.rhs_expr, sol.primal)
    if !is_finite_double(rhs_v) {
      return 1.0e300
    }
    let mut sum_sq = 0.0
    for k = 0; k < sc.norm_exprs.length(); k = k + 1 {
      let ek = eval_linear_expr(sc.norm_exprs[k], sol.primal)
      if !is_finite_double(ek) {
        return 1.0e300
      }
      sum_sq = sum_sq + ek * ek
    }
    let nrm = d_sqrt(sum_sq)
    if !is_finite_double(nrm) {
      return 1.0e300
    }
    let cv = if nrm > rhs_v { nrm - rhs_v } else { 0.0 }
    if cv > max_v {
      max_v = cv
    }
    if -rhs_v > max_v {
      max_v = -rhs_v
    }
  }
  max_v
}

///|
/// High-level entry point: automatically dispatches across all 6 optimization classes:
///   - LP (Two-Phase Primal Simplex + Ruiz Scaling)
///   - MILP (Best-Bound Branch-and-Bound)
///   - Convex QP (Scale-Normalized Mehrotra Primal-Dual IPM)
///   - MIQP (Branch-and-Bound over Convex QP relaxations)
///   - SOCP / QCQP (Nesterov-Todd Primal-Dual Conic IPM)
///   - MISOCP / MIQCQP (Branch-and-Bound over SOCP relaxations)
pub fn Model::solve(
  self : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  let has_s = self.has_soc_constraints()
  let has_q = self.has_quad_obj()
  let has_i = self.has_integer_vars()
  if has_s && has_i {
    solve_misocp(self, config~)
  } else if has_s {
    solve_socp(self, config~)
  } else if has_q && has_i {
    solve_miqp(self, config~)
  } else if has_q {
    solve_qp(self, config~)
  } else if has_i {
    solve_milp(self, config~)
  } else {
    solve_lp(self, config~)
  }
}
