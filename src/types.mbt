///|
/// Core types for moonlp: Linear Programming (LP) & Mixed-Integer Linear Programming (MILP).

///|
pub(all) enum Sense {
  Minimize
  Maximize
} derive(Show, Eq)

///|
pub(all) enum Relation {
  Le // <=
  Ge // >=
  Eq // ==
} derive(Show, Eq)

///|
pub(all) enum VarKind {
  Continuous
  Integer
  Binary
} derive(Show, Eq)

///|
pub(all) enum Status {
  Optimal
  Infeasible
  Unbounded
  MaxIterations
  NodeLimit
} derive(Show, Eq)

///|
pub(all) struct VarId {
  index : Int
} derive(Show, Eq)

///|
pub(all) struct Term {
  var_id : VarId
  coeff : Double
} derive(Show, Eq)

///|
pub(all) struct VarSpec {
  name : String
  lb : Double // Double::neg_infinity() for unbounded below
  ub : Double // Double::infinity() for unbounded above
  kind : VarKind
  obj : Double
} derive(Show)

///|
pub(all) struct Constraint {
  name : String
  terms : Array[Term]
  relation : Relation
  rhs : Double
} derive(Show)

///|
pub(all) struct SolverConfig {
  tol : Double
  int_tol : Double
  max_iters : Int
  max_nodes : Int
  mip_gap_tol : Double
  bland_stall_limit : Int
} derive(Show)

///|
pub fn SolverConfig::default() -> SolverConfig {
  {
    tol: 1.0e-9,
    int_tol: 1.0e-6,
    max_iters: 20000,
    max_nodes: 5000,
    mip_gap_tol: 1.0e-6,
    bland_stall_limit: 25,
  }
}

///|
pub(all) struct Solution {
  status : Status
  objective : Double
  primal : Array[Double]
  dual : Array[Double]
  reduced_costs : Array[Double]
  iterations : Int
  nodes_explored : Int
  mip_gap : Double
} derive(Show)

///|
pub fn Solution::empty(status : Status, iterations : Int, nodes : Int) -> Solution {
  {
    status,
    objective: 0.0,
    primal: [],
    dual: [],
    reduced_costs: [],
    iterations,
    nodes_explored: nodes,
    mip_gap: 0.0,
  }
}

///|
/// Algebraic linear expression builder for ergonomic constraint & objective construction.
pub struct LinearExpr {
  terms : Array[Term]
  mut constant : Double
} derive(Show)

///|
pub fn LinearExpr::zero() -> LinearExpr {
  { terms: [], constant: 0.0 }
}

///|
pub fn LinearExpr::from_var(v : VarId, coeff? : Double = 1.0) -> LinearExpr {
  { terms: [{ var_id: v, coeff }], constant: 0.0 }
}

///|
pub fn LinearExpr::add_term(
  self : LinearExpr,
  v : VarId,
  coeff : Double,
) -> LinearExpr {
  if coeff != 0.0 {
    self.terms.push({ var_id: v, coeff })
  }
  self
}

///|
pub fn LinearExpr::add_const(self : LinearExpr, c : Double) -> LinearExpr {
  self.constant = self.constant + c
  self
}

///|
pub fn LinearExpr::scale(self : LinearExpr, factor : Double) -> LinearExpr {
  let out : Array[Term] = []
  for i = 0; i < self.terms.length(); i = i + 1 {
    let t = self.terms[i]
    out.push({ var_id: t.var_id, coeff: t.coeff * factor })
  }
  { terms: out, constant: self.constant * factor }
}

///|
pub fn LinearExpr::add_expr(
  self : LinearExpr,
  other : LinearExpr,
) -> LinearExpr {
  let out : Array[Term] = []
  for i = 0; i < self.terms.length(); i = i + 1 {
    out.push(self.terms[i])
  }
  for i = 0; i < other.terms.length(); i = i + 1 {
    out.push(other.terms[i])
  }
  { terms: out, constant: self.constant + other.constant }
}

///|
pub fn LinearExpr::sub_expr(
  self : LinearExpr,
  other : LinearExpr,
) -> LinearExpr {
  self.add_expr(other.scale(-1.0))
}

///|
pub suberror ParseError {
  ParseError(String)
} derive(Show, Eq)
