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
fn Model::presolve_with_map(
  self : Model,
  tol : Double,
) -> (
  Model,
  PresolveReport,
  Array[Int32],
  Array[(Int32, Double)],
  Array[(Int32, Double)],
) {
  let n_vars = self.vars.length()
  let new_vars : Array[VarSpec] = []
  let lb_source : Array[(Int32, Double)] = []
  let ub_source : Array[(Int32, Double)] = []
  for i = 0; i < n_vars; i = i + 1 {
    new_vars.push(self.vars[i])
    lb_source.push((-1, 0.0))
    ub_source.push((-1, 0.0))
  }
  let mut tightened = 0
  let mut removed = 0
  let kept_constraints : Array[Constraint] = []
  let kept_orig_indices : Array[Int32] = []

  // First check existing bounds
  for i = 0; i < n_vars; i = i + 1 {
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
        lb_source,
        ub_source,
      )
    }
  }

  for c_idx = 0; c_idx < self.constraints.length(); c_idx = c_idx + 1 {
    let c = self.constraints[c_idx]
    // Aggregate duplicate terms in the row without discarding small nonzero coefficients!
    let mut nonzeros = 0
    let mut single_var = -1
    let mut single_coeff = 0.0
    for j = 0; j < c.terms.length(); j = j + 1 {
      let t = c.terms[j]
      if t.coeff != 0.0 {
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

    if nonzeros == 0 || (nonzeros == 1 && single_coeff == 0.0) {
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
          lb_source,
          ub_source,
        )
      }
      removed = removed + 1
      continue
    }

    if nonzeros == 1 && single_coeff.abs() > 1.0e-30 {
      let v = new_vars[single_var]
      let bound_val = c.rhs / single_coeff
      if !is_finite_double(bound_val) {
        kept_constraints.push(c)
        kept_orig_indices.push(c_idx)
        continue
      }
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
        Le =>
          if bound_val < new_ub {
            new_ub = bound_val
            ub_source[single_var] = (c_idx, single_coeff)
          }
        Ge =>
          if bound_val > new_lb {
            new_lb = bound_val
            lb_source[single_var] = (c_idx, single_coeff)
          }
        Eq => {
          if bound_val > new_lb {
            new_lb = bound_val
            lb_source[single_var] = (c_idx, single_coeff)
          }
          if bound_val < new_ub {
            new_ub = bound_val
            ub_source[single_var] = (c_idx, single_coeff)
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
          lb_source,
          ub_source,
        )
      }
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
    lb_source,
    ub_source,
  )
}

///|
/// Tighten variable bounds from singleton constraints (a * x <= / >= / == b),
/// round integer/binary bounds, and drop trivially redundant rows.
pub fn Model::presolve(
  self : Model,
  tol~ : Double = 1.0e-9,
) -> (Model, PresolveReport) {
  let (reduced, report, _, _, _) = self.presolve_with_map(tol)
  (reduced, report)
}

///|
/// Solve a Model after applying Presolve bound tightening & redundant-row elimination,
/// recovering exact original-space dual shadow prices and reduced costs for LP models.
pub fn Model::solve_presolved(
  self : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  let (reduced, report, kept_orig_indices, lb_source, ub_source) = self.presolve_with_map(
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

  // Postsolve: map kept constraint duals back to original indices & assign exact singleton sources
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
    let mut absorbed = false
    if s_j > config.tol {
      let (src_idx, src_coeff) = lb_source[j]
      if src_idx >= 0 && src_coeff != 0.0 {
        orig_dual[src_idx] = s_j / (sense_sign * src_coeff)
        absorbed = true
      }
    } else if s_j < -config.tol {
      let (src_idx, src_coeff) = ub_source[j]
      if src_idx >= 0 && src_coeff != 0.0 {
        orig_dual[src_idx] = s_j / (sense_sign * src_coeff)
        absorbed = true
      }
    }
    orig_rc[j] = if absorbed { 0.0 } else { self.vars[j].obj - at_y[j] }
  }

  { ..rsol, dual: orig_dual, reduced_costs: orig_rc, }
}
