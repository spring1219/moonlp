///|
/// Gomory Mixed-Integer (GMI) Cutting Planes + Dual Simplex Re-Optimization
/// for Mixed-Integer Linear Programming (Branch-and-Cut).

///|
#derive(Debug)
pub(all) struct GmiCutReport {
  cuts_added : Int32
  root_lp_obj : Double
  cut_lp_obj : Double
  dual_pivots : Int32
  solved_at_root : Bool
}

///|
fn is_near_int(x : Double, tol : Double) -> Bool {
  (x - d_round(x)).abs() <= tol
}

///|
/// Solve a MILP using Root Gomory Mixed-Integer (GMI) Cutting Planes with
/// Dual Simplex re-optimization, seamlessly falling back to Best-Bound
/// Branch-and-Bound if fractional variables remain after `max_cuts` rounds.
pub fn solve_milp_with_cuts(
  model : Model,
  max_cuts~ : Int32 = 20,
  config~ : SolverConfig = SolverConfig::default(),
) -> (Solution, GmiCutReport) {
  let empty_report : GmiCutReport = {
    cuts_added: 0,
    root_lp_obj: 0.0,
    cut_lp_obj: 0.0,
    dual_pivots: 0,
    solved_at_root: false,
  }

  if !model.has_integer_vars() {
    let sol = solve_lp(model, config~)
    return (
      sol,
      {
        ..empty_report,
        root_lp_obj: sol.objective,
        cut_lp_obj: sol.objective,
        solved_at_root: sol.status == Optimal,
      },
    )
  }

  let n_orig = model.vars.length()
  let m_user = model.constraints.length()
  let sense_sign = match model.sense {
    Minimize => 1.0
    Maximize => -1.0
  }

  // Step 1: Standardize non-negative shifted variables & track which columns are integer
  let lbs : Array[Double] = []
  let std_obj : Array[Double] = []
  let is_int_col : Array[Bool] = []
  let ub_rows : Array[(Int32, Double, Bool)] = [] // (col, delta, is_int_slack)
  let mut obj_offset = model.obj_constant * sense_sign

  for j = 0; j < n_orig; j = j + 1 {
    let v = model.vars[j]
    let eff_lb = match v.kind {
      Binary => if v.lb > 0.0 { d_ceil(v.lb) } else { 0.0 }
      Integer => if v.lb > -inf_threshold { d_ceil(v.lb) } else { v.lb }
      Continuous => v.lb
    }
    let eff_ub = match v.kind {
      Binary => if v.ub < 1.0 { d_floor(v.ub) } else { 1.0 }
      Integer => if v.ub < inf_threshold { d_floor(v.ub) } else { v.ub }
      Continuous => v.ub
    }
    if eff_lb > eff_ub + config.tol {
      return (Solution::empty(Infeasible, 0, 1), empty_report)
    }
    // For general free/neg-infinite lower bounds, delegate to Branch-and-Bound
    if eff_lb <= -inf_threshold {
      let sol = solve_milp(model, config~)
      return (sol, empty_report)
    }
    let c_j = v.obj * sense_sign
    lbs.push(eff_lb)
    std_obj.push(c_j)
    let is_int = is_integer_var(v.kind)
    is_int_col.push(is_int)
    obj_offset = obj_offset + c_j * eff_lb
    if eff_ub < inf_threshold {
      ub_rows.push((j, eff_ub - eff_lb, is_int))
    }
  }

  let n_std = n_orig
  let m_init = m_user + ub_rows.length()
  if m_init == 0 {
    let sol = solve_milp(model, config~)
    return (sol, empty_report)
  }

  // Step 2: Build transformed constraint rows (RHS >= 0) & detect integral slacks
  let row_coeffs : Array[Array[Double]] = []
  let row_rhs : Array[Double] = []
  let row_rels : Array[Relation] = []
  let row_int_slack : Array[Bool] = []

  for i = 0; i < m_user; i = i + 1 {
    let c = model.constraints[i]
    let coeffs = Array::make(n_std, 0.0)
    let mut rhs = c.rhs
    for k = 0; k < c.terms.length(); k = k + 1 {
      let t = c.terms[k]
      coeffs[t.var_id.index] = coeffs[t.var_id.index] + t.coeff
      rhs = rhs - t.coeff * lbs[t.var_id.index]
    }
    let mut rel = c.relation
    if rhs < 0.0 {
      for j = 0; j < n_std; j = j + 1 {
        coeffs[j] = -coeffs[j]
      }
      rhs = -rhs
      rel = match c.relation {
        Le => Ge
        Ge => Le
        Eq => Eq
      }
    }
    // A slack/surplus variable is provably integral iff RHS is integer and every
    // nonzero term is an integer variable with an integer coefficient.
    let mut slack_int = is_near_int(rhs, config.tol)
    if slack_int {
      for j = 0; j < n_std; j = j + 1 {
        if coeffs[j].abs() > config.tol {
          if !is_int_col[j] || !is_near_int(coeffs[j], config.tol) {
            slack_int = false
            break
          }
        }
      }
    }
    row_coeffs.push(coeffs)
    row_rhs.push(rhs)
    row_rels.push(rel)
    row_int_slack.push(slack_int)
  }

  for i = 0; i < ub_rows.length(); i = i + 1 {
    let (col, delta, is_int_s) = ub_rows[i]
    let coeffs = Array::make(n_std, 0.0)
    coeffs[col] = 1.0
    row_coeffs.push(coeffs)
    row_rhs.push(delta)
    row_rels.push(Le)
    row_int_slack.push(is_int_s && is_near_int(delta, config.tol))
  }

  // Step 3: Allocate slack, surplus, and artificial columns
  let mut num_cols = n_std
  let row_slack_col = Array::make(m_init, -1)
  let row_art_col = Array::make(m_init, -1)
  let is_artificial = Array::make(n_std, false)

  for i = 0; i < m_init; i = i + 1 {
    match row_rels[i] {
      Le => {
        row_slack_col[i] = num_cols
        is_artificial.push(false)
        is_int_col.push(row_int_slack[i])
        num_cols = num_cols + 1
      }
      Ge => {
        row_slack_col[i] = num_cols
        is_artificial.push(false)
        is_int_col.push(row_int_slack[i])
        num_cols = num_cols + 1
        row_art_col[i] = num_cols
        is_artificial.push(true)
        is_int_col.push(false)
        num_cols = num_cols + 1
      }
      Eq => {
        row_art_col[i] = num_cols
        is_artificial.push(true)
        is_int_col.push(false)
        num_cols = num_cols + 1
      }
    }
  }

  let mut m = m_init
  let mut rhs_col = num_cols
  let tab : Array[Array[Double]] = []
  let basis : Array[Int32] = Array::make(m, -1)
  let mut has_art = false

  for i = 0; i < m; i = i + 1 {
    let row = Array::make(num_cols + 1, 0.0)
    for j = 0; j < n_std; j = j + 1 {
      row[j] = row_coeffs[i][j]
    }
    match row_rels[i] {
      Le => {
        let s = row_slack_col[i]
        row[s] = 1.0
        basis[i] = s
      }
      Ge => {
        let s = row_slack_col[i]
        let a = row_art_col[i]
        row[s] = -1.0
        row[a] = 1.0
        basis[i] = a
        has_art = true
      }
      Eq => {
        let a = row_art_col[i]
        row[a] = 1.0
        basis[i] = a
        has_art = true
      }
    }
    row[rhs_col] = row_rhs[i]
    tab.push(row)
  }
  tab.push(Array::make(num_cols + 1, 0.0))

  let col_active = Array::make(num_cols, true)
  let mut total_iters = 0

  // Phase I if needed
  if has_art {
    let obj_row = tab[m]
    for j = 0; j < num_cols; j = j + 1 {
      obj_row[j] = if is_artificial[j] { 1.0 } else { 0.0 }
    }
    obj_row[rhs_col] = 0.0
    for i = 0; i < m; i = i + 1 {
      if is_artificial[basis[i]] {
        for j = 0; j <= rhs_col; j = j + 1 {
          obj_row[j] = obj_row[j] - tab[i][j]
        }
      }
    }
    let (p1_st, p1_it) = run_simplex_phase(tab, basis, col_active, m, num_cols, config)
    total_iters = total_iters + p1_it
    if p1_st != Optimal || -tab[m][rhs_col] > config.tol * 100.0 {
      return (Solution::empty(Infeasible, total_iters, 1), empty_report)
    }
    for i = 0; i < m; i = i + 1 {
      if is_artificial[basis[i]] {
        for j = 0; j < num_cols; j = j + 1 {
          if !is_artificial[j] && tab[i][j].abs() > config.tol {
            pivot_tableau(tab, basis, m, rhs_col, i, j)
            break
          }
        }
      }
    }
    for j = 0; j < num_cols; j = j + 1 {
      if is_artificial[j] {
        col_active[j] = false
      }
    }
  }

  // Phase II root LP solve
  let obj_row = tab[m]
  for j = 0; j <= rhs_col; j = j + 1 {
    obj_row[j] = 0.0
  }
  for j = 0; j < n_std; j = j + 1 {
    obj_row[j] = std_obj[j]
  }
  for i = 0; i < m; i = i + 1 {
    let bc = basis[i]
    let cb = if bc < n_std { std_obj[bc] } else { 0.0 }
    if cb != 0.0 {
      for j = 0; j <= rhs_col; j = j + 1 {
        obj_row[j] = obj_row[j] - cb * tab[i][j]
      }
    }
  }

  let (p2_st, p2_it) = run_simplex_phase(tab, basis, col_active, m, num_cols, config)
  total_iters = total_iters + p2_it
  if p2_st != Optimal {
    return (Solution::empty(p2_st, total_iters, 1), empty_report)
  }

  let root_lp_obj = (-tab[m][rhs_col] + obj_offset) * sense_sign

  // Step 4: Iterative Gomory Mixed-Integer (GMI) Cuts + Dual Simplex Re-Optimization
  let mut cuts_added = 0
  let mut dual_pivots = 0
  let mut is_integral = false

  while cuts_added <= max_cuts {
    // Find basic structural integer variable with most fractional RHS
    let mut best_r = -1
    let mut best_dist = 0.0
    for i = 0; i < m; i = i + 1 {
      let bc = basis[i]
      if bc >= 0 && bc < n_std && is_int_col[bc] {
        let val = tab[i][rhs_col]
        let f = val - d_floor(val)
        if f > config.int_tol && 1.0 - f > config.int_tol {
          let dist = if f < 1.0 - f { f } else { 1.0 - f }
          if dist > best_dist {
            best_dist = dist
            best_r = i
          }
        }
      }
    }

    if best_r == -1 {
      is_integral = true
      break
    }
    if cuts_added == max_cuts {
      break
    }

    // Generate GMI cut from row `best_r`
    let rhs_val = tab[best_r][rhs_col]
    let f0 = rhs_val - d_floor(rhs_val)
    let scale = f0 / (1.0 - f0)
    let new_col = num_cols
    num_cols = num_cols + 1
    rhs_col = num_cols

    // Expand all existing rows + obj row by 1 column before RHS
    for i = 0; i <= m; i = i + 1 {
      let row_i = tab[i]
      let old_rhs = row_i[new_col]
      row_i[new_col] = 0.0
      row_i.push(old_rhs)
    }

    let cut_row = Array::make(num_cols + 1, 0.0)
    let mut cut_has_coeff = false
    for j = 0; j < new_col; j = j + 1 {
      if !col_active[j] || j == basis[best_r] {
        continue
      }
      let a = tab[best_r][j]
      if a.abs() <= config.tol {
        continue
      }
      let gj = if is_int_col[j] {
        let fj = a - d_floor(a)
        if fj <= config.tol || 1.0 - fj <= config.tol {
          0.0
        } else if fj <= f0 {
          fj
        } else {
          scale * (1.0 - fj)
        }
      } else if a >= 0.0 {
        a
      } else {
        -scale * a
      }
      if gj > config.tol {
        cut_row[j] = -gj
        cut_has_coeff = true
      }
    }

    if !cut_has_coeff {
      // Degenerate cut; abort cut loop and fall back to B&B
      let _ = tab[m].pop()
      break
    }

    cut_row[new_col] = 1.0
    cut_row[rhs_col] = -f0

    // Insert cut_row before obj_row
    let saved_obj = tab[m]
    tab[m] = cut_row
    tab.push(saved_obj)
    basis.push(new_col)
    col_active.push(true)
    is_int_col.push(false)
    m = m + 1
    cuts_added = cuts_added + 1

    // Re-optimize immediately via Dual Simplex (NO Phase I!)
    let (dst, dit) = run_dual_simplex(tab, basis, col_active, m, num_cols, config)
    dual_pivots = dual_pivots + dit
    total_iters = total_iters + dit
    if dst == Infeasible {
      return (
        Solution::empty(Infeasible, total_iters, 1),
        {
          cuts_added,
          root_lp_obj,
          cut_lp_obj: root_lp_obj,
          dual_pivots,
          solved_at_root: true,
        },
      )
    }
    if dst != Optimal {
      break
    }
  }

  let cut_lp_obj = (-tab[m][rhs_col] + obj_offset) * sense_sign

  if is_integral {
    let primal = Array::make(n_orig, 0.0)
    for j = 0; j < n_orig; j = j + 1 {
      primal[j] = lbs[j]
    }
    for i = 0; i < m; i = i + 1 {
      let bc = basis[i]
      if bc >= 0 && bc < n_std {
        let raw = lbs[bc] + tab[i][rhs_col]
        primal[bc] = if is_int_col[bc] { d_round(raw) } else { raw }
      }
    }
    let mut exact_obj = model.obj_constant
    for j = 0; j < n_orig; j = j + 1 {
      exact_obj = exact_obj + model.vars[j].obj * primal[j]
    }
    let root_sol : Solution = {
      status: Optimal,
      objective: exact_obj,
      primal,
      dual: [],
      reduced_costs: [],
      iterations: total_iters,
      nodes_explored: 1,
      mip_gap: 0.0,
    }
    // Safety verification: guarantee 100% primal feasibility before accepting root cut solution
    if model.verify_solution(root_sol, tol=1.0e-5) {
      return (
        root_sol,
        {
          cuts_added,
          root_lp_obj,
          cut_lp_obj: exact_obj,
          dual_pivots,
          solved_at_root: true,
        },
      )
    }
  }

  // Fallback to exact Best-Bound Branch-and-Bound if fractionality persists
  let bb_sol = solve_milp(model, config~)
  (
    { ..bb_sol, iterations: bb_sol.iterations + total_iters, },
    {
      cuts_added,
      root_lp_obj,
      cut_lp_obj,
      dual_pivots,
      solved_at_root: false,
    },
  )
}
