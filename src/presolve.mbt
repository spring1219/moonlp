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
priv struct RemovedSingleton {
  orig_idx : Int32
  var_idx : Int32
  coeff : Double
  eff_rel : Relation
  bound_val : Double
}

///|
fn Model::presolve_with_map(
  self : Model,
  tol : Double,
) -> (Model, PresolveReport, Array[Int32], Array[RemovedSingleton]) {
  let new_vars : Array[VarSpec] = []
  for i = 0; i < self.vars.length(); i = i + 1 {
    new_vars.push(self.vars[i])
  }
  let mut tightened = 0
  let mut removed = 0
  let kept_constraints : Array[Constraint] = []
  let kept_orig_indices : Array[Int32] = []
  let removed_singletons : Array[RemovedSingleton] = []

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
        kept_orig_indices,
        removed_singletons,
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
          kept_orig_indices,
          removed_singletons,
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
          kept_orig_indices,
          removed_singletons,
        )
      }
      removed_singletons.push({
        orig_idx: c_idx,
        var_idx: single_var,
        coeff: single_coeff,
        eff_rel,
        bound_val,
      })
      removed = removed + 1
    } else {
      kept_constraints.push(c)
      kept_orig_indices.push(c_idx)
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
    kept_orig_indices,
    removed_singletons,
  )
}

///|
/// Tighten variable bounds from singleton constraints (a * x <= / >= / == b),
/// round integer/binary bounds, and drop trivially redundant rows.
pub fn Model::presolve(
  self : Model,
  tol~ : Double = 1.0e-9,
) -> (Model, PresolveReport) {
  let (reduced, report, _, _) = self.presolve_with_map(tol)
  (reduced, report)
}

///|
/// Solve a Model after applying Presolve bound tightening & redundant-row elimination,
/// recovering original-space dual shadow prices and reduced costs for LP models.
pub fn Model::solve_presolved(
  self : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  let (reduced, report, kept_orig_indices, removed_singletons) = self.presolve_with_map(
    config.tol,
  )
  if report.infeasible_detected {
    return Solution::empty(Infeasible, 0, 1)
  }
  let rsol = reduced.solve(config~)
  if rsol.status != Optimal ||
    self.has_integer_vars() ||
    rsol.dual.length() != kept_orig_indices.length() {
    return rsol
  }

  // Postsolve: map kept constraint duals back to original indices & recover singleton duals
  let m_orig = self.constraints.length()
  let n_orig = self.vars.length()
  let orig_dual = Array::make(m_orig, 0.0)
  for k = 0; k < kept_orig_indices.length(); k = k + 1 {
    orig_dual[kept_orig_indices[k]] = rsol.dual[k]
  }

  let at_y = Array::make(n_orig, 0.0)
  for k = 0; k < kept_orig_indices.length(); k = k + 1 {
    let orig_i = kept_orig_indices[k]
    let y_i = orig_dual[orig_i]
    if y_i != 0.0 {
      let c = self.constraints[orig_i]
      for t_idx = 0; t_idx < c.terms.length(); t_idx = t_idx + 1 {
        let t = c.terms[t_idx]
        at_y[t.var_id.index] = at_y[t.var_id.index] + t.coeff * y_i
      }
    }
  }

  let sense_sign = match self.sense {
    Minimize => 1.0
    Maximize => -1.0
  }
  let orig_rc = Array::make(n_orig, 0.0)
  for j = 0; j < n_orig; j = j + 1 {
    let s_j = sense_sign * (self.vars[j].obj - at_y[j])
    let orig_v = self.vars[j]
    let red_v = reduced.vars[j]
    let mut absorbed = false
    if s_j > config.tol && red_v.lb > orig_v.lb + config.tol {
      let btol = 1.0e-7 * (1.0 + red_v.lb.abs())
      for r = 0; r < removed_singletons.length(); r = r + 1 {
        let rs = removed_singletons[r]
        if rs.var_idx == j &&
          (rs.eff_rel == Ge || rs.eff_rel == Eq) &&
          (rs.bound_val - red_v.lb).abs() <= btol {
          orig_dual[rs.orig_idx] = s_j / (sense_sign * rs.coeff)
          absorbed = true
          break
        }
      }
    } else if s_j < -config.tol && red_v.ub < orig_v.ub - config.tol {
      let btol = 1.0e-7 * (1.0 + red_v.ub.abs())
      for r = 0; r < removed_singletons.length(); r = r + 1 {
        let rs = removed_singletons[r]
        if rs.var_idx == j &&
          (rs.eff_rel == Le || rs.eff_rel == Eq) &&
          (rs.bound_val - red_v.ub).abs() <= btol {
          orig_dual[rs.orig_idx] = s_j / (sense_sign * rs.coeff)
          absorbed = true
          break
        }
      }
    }
    orig_rc[j] = if absorbed { 0.0 } else { self.vars[j].obj - at_y[j] }
  }

  { ..rsol, dual: orig_dual, reduced_costs: orig_rc, }
}
