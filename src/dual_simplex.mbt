///|
/// Dual Simplex Engine for moonlp v0.3.0:
/// - Maintains dual feasibility (all active reduced costs >= -tol) while driving
///   primal infeasibilities (negative RHS entries) to >= 0.
/// - Ideal for warm-start re-optimization after adding Gomory cutting planes or
///   Branch-and-Bound bound cuts (zero Phase I restarts needed).

///|
/// Run Dual Simplex pivots in-place on a dual-feasible minimization tableau.
/// Returns (status, iterations_used): Optimal | Infeasible | MaxIterations.
fn run_dual_simplex(
  tab : Array[Array[Double]],
  basis : Array[Int32],
  col_active : Array[Bool],
  m : Int32,
  num_cols : Int32,
  config : SolverConfig,
) -> (Status, Int32) {
  let rhs_col = num_cols
  let mut iters = 0
  let mut stall_count = 0
  let mut prev_obj = tab[m][rhs_col]
  let nz_cols : Array[Int32] = []
  let nz_vals : Array[Double] = []

  while iters < config.max_iters {
    let use_bland = stall_count >= config.bland_stall_limit

    // 1. Choose leaving row: negative RHS entry (most negative, or Bland smallest basis var on stall)
    let mut leave_row = -1
    let mut most_neg = -config.tol
    for i = 0; i < m; i = i + 1 {
      let rhs_val = tab[i][rhs_col]
      if rhs_val < -config.tol {
        if use_bland {
          if leave_row == -1 || basis[i] < basis[leave_row] {
            leave_row = i
          }
        } else if rhs_val < most_neg {
          most_neg = rhs_val
          leave_row = i
        }
      }
    }

    if leave_row == -1 {
      return (Optimal, iters)
    }

    // 2. Dual minimum ratio test with tight relative tie tolerance (1e-12) and Harris pivot guard
    let pivot_tol = 1.0e-12
    let mut enter_col = -1
    let mut best_ratio = 1.0e308
    let mut best_pivot = 0.0
    let obj_row = tab[m]
    for j = 0; j < num_cols; j = j + 1 {
      if !col_active[j] {
        continue
      }
      let a_lj = tab[leave_row][j]
      if a_lj < -pivot_tol {
        let rc = if obj_row[j] < 0.0 { 0.0 } else { obj_row[j] }
        let ratio = rc / -a_lj
        let tie_tol = 1.0e-12 * best_ratio.abs()
        if enter_col == -1 || ratio < best_ratio - tie_tol {
          best_ratio = ratio
          best_pivot = -a_lj
          enter_col = j
        } else if (ratio - best_ratio).abs() <= tie_tol {
          if use_bland {
            if j < enter_col {
              enter_col = j
              best_pivot = -a_lj
            }
          } else if -a_lj > best_pivot {
            enter_col = j
            best_pivot = -a_lj
          }
        }
      }
    }

    if enter_col == -1 {
      // Dual unbounded => Primal Infeasible
      return (Infeasible, iters)
    }

    // 3. Pivot (zero-allocation sparse update)
    pivot_tableau_buf(
      tab, basis, m, rhs_col, leave_row, enter_col, nz_cols, nz_vals,
    )
    iters = iters + 1

    let cur_obj = tab[m][rhs_col]
    if (cur_obj - prev_obj).abs() <= config.tol {
      stall_count = stall_count + 1
    } else {
      stall_count = 0
      prev_obj = cur_obj
    }
  }

  (MaxIterations, iters)
}

///|
/// Solve an LP directly via Dual Simplex when the initial slack basis is dual-feasible
/// (e.g., Minimize c^T x with c >= 0 and >= / <= constraints, zero Phase I overhead!).
/// Automatically falls back to Two-Phase Primal Simplex if the initial basis is not dual-feasible.
pub fn solve_lp_dual(
  model : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  let n_orig = model.vars.length()
  let m_user = model.constraints.length()
  let sense_sign = match model.sense {
    Minimize => 1.0
    Maximize => -1.0
  }

  // Fast eligibility check: non-negative lower-bounded variables with dual-feasible reduced costs
  // and inequality constraints only (Le / Ge)
  for j = 0; j < n_orig; j = j + 1 {
    let v = model.vars[j]
    if v.lb > v.ub + config.tol {
      return Solution::empty(Infeasible, 0, 1)
    }
    if v.lb < -inf_threshold || v.ub < inf_threshold {
      return solve_lp(model, config~)
    }
    if v.obj * sense_sign < 0.0 {
      return solve_lp(model, config~)
    }
  }
  for i = 0; i < m_user; i = i + 1 {
    match model.constraints[i].relation {
      Eq => return solve_lp(model, config~)
      Le | Ge => ()
    }
  }
  if m_user == 0 || needs_auto_scaling(model) {
    return solve_lp(model, config~)
  }

  // Build dual-feasible initial tableau: convert Ge (a^T x >= b) to (-a^T x + s = -b)
  let num_cols = n_orig + m_user
  let rhs_col = num_cols
  let tab : Array[Array[Double]] = []
  let basis : Array[Int32] = Array::make(m_user, -1)
  let row_signs : Array[Double] = Array::make(m_user, 1.0)
  let mut obj_offset = model.obj_constant * sense_sign

  for j = 0; j < n_orig; j = j + 1 {
    obj_offset = obj_offset + model.vars[j].obj * sense_sign * model.vars[j].lb
  }

  for i = 0; i < m_user; i = i + 1 {
    let c = model.constraints[i]
    let row = Array::make(num_cols + 1, 0.0)
    let rsign = match c.relation {
      Le => 1.0
      Ge => -1.0
      Eq => 1.0
    }
    row_signs[i] = rsign
    let mut rhs = c.rhs
    for k = 0; k < c.terms.length(); k = k + 1 {
      let t = c.terms[k]
      let lb = model.vars[t.var_id.index].lb
      row[t.var_id.index] = row[t.var_id.index] + rsign * t.coeff
      rhs = rhs - t.coeff * lb
    }
    let scol = n_orig + i
    row[scol] = 1.0
    row[rhs_col] = rsign * rhs
    basis[i] = scol
    tab.push(row)
  }

  let obj_row = Array::make(num_cols + 1, 0.0)
  for j = 0; j < n_orig; j = j + 1 {
    obj_row[j] = model.vars[j].obj * sense_sign
  }
  tab.push(obj_row)

  let col_active = Array::make(num_cols, true)
  let (status, iters) = run_dual_simplex(
    tab, basis, col_active, m_user, num_cols, config,
  )
  if status != Optimal {
    return Solution::empty(status, iters, 1)
  }

  let primal = Array::make(n_orig, 0.0)
  for j = 0; j < n_orig; j = j + 1 {
    primal[j] = model.vars[j].lb
  }
  for i = 0; i < m_user; i = i + 1 {
    let bcol = basis[i]
    if bcol >= 0 && bcol < n_orig {
      let val = tab[i][rhs_col]
      primal[bcol] = model.vars[bcol].lb +
        (if val < 0.0 && val > -config.tol { 0.0 } else { val })
    }
  }

  let reduced_costs : Array[Double] = []
  for j = 0; j < n_orig; j = j + 1 {
    reduced_costs.push(tab[m_user][j] * sense_sign)
  }

  let dual : Array[Double] = []
  for i = 0; i < m_user; i = i + 1 {
    let y_std = -tab[m_user][n_orig + i]
    let y_user = y_std * row_signs[i] * sense_sign
    dual.push(y_user)
  }

  let final_obj = (-tab[m_user][rhs_col] + obj_offset) * sense_sign
  let result : Solution = {
    status: Optimal,
    objective: final_obj,
    primal,
    dual,
    reduced_costs,
    iterations: iters,
    nodes_explored: 1,
    mip_gap: 0.0,
  }
  clean_dual_roundoff(model, result, config)
}
