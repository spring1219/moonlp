///|
/// Presolver & Post-Optimality Sensitivity Analysis for moonlp.

///|
#derive(Debug)
pub(all) struct PresolveReport {
  bounds_tightened : Int32
  redundant_rows_removed : Int32
  infeasible_detected : Bool
}

///|
/// Tighten variable bounds from singleton constraints (a * x <= / >= / == b),
/// round integer/binary bounds, and drop trivially redundant rows.
pub fn Model::presolve(
  self : Model,
  tol~ : Double = 1.0e-9,
) -> (Model, PresolveReport) {
  let new_vars : Array[VarSpec] = []
  for i = 0; i < self.vars.length(); i = i + 1 {
    new_vars.push(self.vars[i])
  }
  let mut tightened = 0
  let mut removed = 0
  let kept_constraints : Array[Constraint] = []

  // First check existing bounds
  for i = 0; i < new_vars.length(); i = i + 1 {
    let v = new_vars[i]
    if v.lb > v.ub + tol {
      return (
        { ..self, vars: new_vars, constraints: kept_constraints, },
        {
          bounds_tightened: tightened,
          redundant_rows_removed: removed,
          infeasible_detected: true,
        },
      )
    }
  }

  for c_idx = 0; c_idx < self.constraints.length(); c_idx = c_idx + 1 {
    let c = self.constraints[c_idx]
    // Aggregate duplicate terms in the row
    let mut nonzeros = 0
    let mut single_var = -1
    let mut single_coeff = 0.0
    for j = 0; j < c.terms.length(); j = j + 1 {
      let t = c.terms[j]
      if t.coeff.abs() > tol {
        if single_var == -1 || single_var == t.var_id.index {
          single_var = t.var_id.index
          single_coeff = single_coeff + t.coeff
          nonzeros = 1
        } else {
          nonzeros = 2
          break
        }
      }
    }

    if nonzeros == 0 {
      // Constant row: 0 <= / >= / == rhs
      let ok = match c.relation {
        Le => 0.0 <= c.rhs + tol
        Ge => 0.0 >= c.rhs - tol
        Eq => c.rhs.abs() <= tol
      }
      if !ok {
        return (
          { ..self, vars: new_vars, constraints: kept_constraints, },
          {
            bounds_tightened: tightened,
            redundant_rows_removed: removed,
            infeasible_detected: true,
          },
        )
      }
      removed = removed + 1
      continue
    }

    if nonzeros == 1 && single_coeff.abs() > tol {
      let v = new_vars[single_var]
      let bound_val = c.rhs / single_coeff
      let mut new_lb = v.lb
      let mut new_ub = v.ub
      let eff_rel = if single_coeff > 0.0 {
        c.relation
      } else {
        match c.relation {
          Le => Ge
          Ge => Le
          Eq => Eq
        }
      }
      match eff_rel {
        Le => if bound_val < new_ub { new_ub = bound_val }
        Ge => if bound_val > new_lb { new_lb = bound_val }
        Eq => {
          if bound_val > new_lb {
            new_lb = bound_val
          }
          if bound_val < new_ub {
            new_ub = bound_val
          }
        }
      }
      // Tighten integer/binary bounds cleanly
      match v.kind {
        Integer | Binary => {
          if new_lb > -1.0e200 {
            new_lb = d_ceil(new_lb - tol)
          }
          if new_ub < 1.0e200 {
            new_ub = d_floor(new_ub + tol)
          }
        }
        Continuous => ()
      }
      if new_lb != v.lb || new_ub != v.ub {
        tightened = tightened + 1
        new_vars[single_var] = { ..v, lb: new_lb, ub: new_ub, }
      }
      if new_lb > new_ub + tol {
        return (
          { ..self, vars: new_vars, constraints: kept_constraints, },
          {
            bounds_tightened: tightened,
            redundant_rows_removed: removed,
            infeasible_detected: true,
          },
        )
      }
      removed = removed + 1
    } else {
      kept_constraints.push(c)
    }
  }

  (
    {
      sense: self.sense,
      vars: new_vars,
      constraints: kept_constraints,
      obj_constant: self.obj_constant,
    },
    {
      bounds_tightened: tightened,
      redundant_rows_removed: removed,
      infeasible_detected: false,
    },
  )
}

///|
/// Solve a Model after applying Presolve bound tightening & redundant-row elimination.
pub fn Model::solve_presolved(
  self : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  let (reduced, report) = self.presolve(tol=config.tol)
  if report.infeasible_detected {
    return Solution::empty(Infeasible, 0, 1)
  }
  reduced.solve(config~)
}
