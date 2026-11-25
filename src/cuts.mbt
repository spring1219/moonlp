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
/// Sound interval propagation over shifted variables u_j = x_j - lbs[j] >= 0
/// to derive valid finite upper bounds for variables constrained by linear inequalities.
fn compute_implied_shifted_ubs(
  model : Model,
  lbs : Array[Double],
  init_ubs : Array[Double],
) -> Array[Double] {
  let n = lbs.length()
  let ubs = init_ubs.copy()
  for _pass = 0; _pass < 3; _pass = _pass + 1 {
    for i = 0; i < model.constraints.length(); i = i + 1 {
      let c = model.constraints[i]
      let coeffs = Array::make(n, 0.0)
      let mut rhs = c.rhs
      for k = 0; k < c.terms.length(); k = k + 1 {
        let t = c.terms[k]
        let idx = t.var_id.index
        coeffs[idx] = coeffs[idx] + t.coeff
        rhs = rhs - t.coeff * lbs[idx]
      }
      // Check <= direction (for Le and Eq): sum coeffs[j] * u_j <= rhs
      if c.relation == Le || c.relation == Eq {
        for j = 0; j < n; j = j + 1 {
          let aj = coeffs[j]
          if aj > 1.0e-12 {
            let mut min_other = 0.0
            let mut ok = true
            for p = 0; p < n; p = p + 1 {
              if p != j && coeffs[p] < 0.0 {
                if ubs[p] >= inf_threshold {
                  ok = false
                  break
                }
                min_other = min_other + coeffs[p] * ubs[p]
              }
            }
            if ok {
              let cand = (rhs - min_other) / aj
              if is_finite_double(cand) && cand >= 0.0 && cand < ubs[j] {
                ubs[j] = cand
              }
            }
          }
        }
      }
      // Check >= direction (for Ge and Eq): sum (-coeffs[j]) * u_j <= -rhs
      if c.relation == Ge || c.relation == Eq {
        for j = 0; j < n; j = j + 1 {
          let aj = -coeffs[j]
          if aj > 1.0e-12 {
            let mut min_other = 0.0
            let mut ok = true
            for p = 0; p < n; p = p + 1 {
              if p != j && -coeffs[p] < 0.0 {
                if ubs[p] >= inf_threshold {
                  ok = false
                  break
                }
                min_other = min_other - coeffs[p] * ubs[p]
              }
            }
            if ok {
              let cand = (-rhs - min_other) / aj
              if is_finite_double(cand) && cand >= 0.0 && cand < ubs[j] {
                ubs[j] = cand
              }
            }
          }
        }
      }
    }
  }
  ubs
}

///|
/// Compute a valid upper bound on the non-negative slack/surplus variable of a shifted row.
fn compute_row_slack_ub(
  coeffs : Array[Double],
  rhs : Double,
  rel : Relation,
  var_ubs : Array[Double],
) -> Double {
  match rel {
    Le => {
      // s = rhs - sum a_j u_j <= rhs - sum_{a_j < 0} a_j * ub_j
      let mut ub = rhs
      for j = 0; j < coeffs.length(); j = j + 1 {
        if coeffs[j] < 0.0 {
          if var_ubs[j] >= inf_threshold {
            return 1.0e300
          }
          ub = ub - coeffs[j] * var_ubs[j]
        }
      }
      if is_finite_double(ub) && ub >= 0.0 {
        ub
      } else {
        1.0e300
      }
    }
    Ge => {
      // s = sum a_j u_j - rhs <= sum_{a_j > 0} a_j * ub_j - rhs
      let mut ub = -rhs
      for j = 0; j < coeffs.length(); j = j + 1 {
        if coeffs[j] > 0.0 {
          if var_ubs[j] >= inf_threshold {
            return 1.0e300
          }
          ub = ub + coeffs[j] * var_ubs[j]
        }
      }
      if is_finite_double(ub) && ub >= 0.0 {
        ub
      } else {
        1.0e300
      }
    }
    Eq => 0.0
  }
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
  if has_extreme_scale(model) {
    let sol = solve_milp(model, config~)
    return (sol, empty_report)
  }

  let n_orig = model.vars.length()
  let m_user = model.constraints.length()
  let sense_sign = match model.sense {
    Minimize => 1.0
    Maximize => -1.0
  }

  // Step 1: Standardize non-negative shifted variables & track which columns are integer
  let lbs : Array[Double] = []
  let init_ubs : Array[Double] = []
  let is_fixed : Array[Bool] = []
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
    let fixed_j = eff_ub < inf_threshold && eff_lb == eff_ub
    let c_j = v.obj * sense_sign
    lbs.push(eff_lb)
    is_fixed.push(fixed_j)
    init_ubs.push(
      if fixed_j {
        0.0
      } else if eff_ub < inf_threshold {
        eff_ub - eff_lb
      } else {
        1.0e300
      },
    )
    std_obj.push(if fixed_j { 0.0 } else { c_j })
    let is_int = is_integer_var(v.kind)
    is_int_col.push(is_int)
    obj_offset = obj_offset + c_j * eff_lb
    if !fixed_j && eff_ub < inf_threshold {
      ub_rows.push((j, eff_ub - eff_lb, is_int))
    }
  }

  let col_ubs = compute_implied_shifted_ubs(model, lbs, init_ubs)

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
      let idx = t.var_id.index
      if !is_fixed[idx] {
        coeffs[idx] = coeffs[idx] + t.coeff
      }
      rhs = rhs - t.coeff * lbs[idx]
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
    // A slack/surplus variable is provably integral ONLY if RHS is an exact integer
    // and EVERY nonzero term (coeffs[j] != 0.0) is an integer variable with an exact integer coeff.
    let mut slack_int = rhs == d_round(rhs)
    if slack_int {
      for j = 0; j < n_std; j = j + 1 {
        if coeffs[j] != 0.0 {
          if !is_int_col[j] || coeffs[j] != d_round(coeffs[j]) {
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
    row_int_slack.push(is_int_s && delta == d_round(delta))
  }

  // Step 3: Allocate slack, surplus, and artificial columns
  let mut num_cols = n_std
  let row_slack_col = Array::make(m_init, -1)
  let row_art_col = Array::make(m_init, -1)
  let is_artificial = Array::make(n_std, false)

  for i = 0; i < m_init; i = i + 1 {
    let slack_ub = compute_row_slack_ub(
      row_coeffs[i],
      row_rhs[i],
      row_rels[i],
      col_ubs,
    )
    match row_rels[i] {
      Le => {
        row_slack_col[i] = num_cols
        is_artificial.push(false)
        is_int_col.push(row_int_slack[i])
        col_ubs.push(slack_ub)
        num_cols = num_cols + 1
      }
      Ge => {
        row_slack_col[i] = num_cols
        is_artificial.push(false)
        is_int_col.push(row_int_slack[i])
        col_ubs.push(slack_ub)
        num_cols = num_cols + 1
        row_art_col[i] = num_cols
        is_artificial.push(true)
        is_int_col.push(false)
        col_ubs.push(0.0)
        num_cols = num_cols + 1
      }
      Eq => {
        row_art_col[i] = num_cols
        is_artificial.push(true)
        is_int_col.push(false)
        col_ubs.push(0.0)
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
  for j = 0; j < n_std; j = j + 1 {
    if is_fixed[j] {
      col_active[j] = false
    }
  }
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
    if p1_st != Optimal {
      return (Solution::empty(p1_st, total_iters, 1), empty_report)
    }
    if -tab[m][rhs_col] > config.tol {
      return (Solution::empty(Infeasible, total_iters, 1), empty_report)
    }
    drive_out_artificials(
      tab, basis, col_active, is_artificial, m, num_cols, config,
    )
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

  let (p2_st, p2_it) = run_simplex_phase(tab, basis, col_active, m, num_cols, config, costs=std_obj)
  total_iters = total_iters + p2_it
  if p2_st == Unbounded {
    return (
      Solution::empty(
        resolve_unbounded_milp_status(model, config),
        total_iters,
        1,
      ),
      empty_report,
    )
  }
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

    // Generate numerically safe GMI cut from row `best_r` (build & validate BEFORE expanding tableau).
    // Never discard any nonzero coefficient without checking its maximum contribution over `col_ubs`.
    let rhs_val = tab[best_r][rhs_col]
    let f0 = rhs_val - d_floor(rhs_val)
    if f0 <= 1.0e-5 || 1.0 - f0 <= 1.0e-5 {
      break
    }
    let scale = f0 / (1.0 - f0)
    let new_col = num_cols
    let cut_row = Array::make(num_cols + 2, 0.0)
    let mut f0_eff = f0 - 1.0e-11
    let mut cut_has_coeff = false
    let mut cut_safe = true
    let mut max_g = 0.0
    let mut min_g = 1.0e308
    let mut cut_slack_ub = -f0_eff
    let mut cut_slack_finite = true

    for j = 0; j < new_col; j = j + 1 {
      if !col_active[j] || j == basis[best_r] {
        continue
      }
      if j < col_ubs.length() && col_ubs[j] == 0.0 {
        continue
      }
      let a = tab[best_r][j]
      if a == 0.0 {
        continue
      }
      let gj = if is_int_col[j] {
        let fj = a - d_floor(a)
        if fj <= 0.0 || fj >= 1.0 {
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
      if gj > 1.0e-12 {
        cut_row[j] = -gj
        cut_has_coeff = true
        if gj > max_g {
          max_g = gj
        }
        if gj < min_g {
          min_g = gj
        }
        if cut_slack_finite {
          if j < col_ubs.length() && col_ubs[j] < inf_threshold {
            cut_slack_ub = cut_slack_ub + gj * col_ubs[j]
          } else {
            cut_slack_finite = false
          }
        }
      } else if gj > 0.0 {
        // Never silently drop a positive GMI coefficient: only relax RHS if col has finite upper bound
        if j < col_ubs.length() && col_ubs[j] < inf_threshold {
          f0_eff = f0_eff - gj * col_ubs[j]
        } else {
          cut_safe = false
          break
        }
      }
    }

    if !cut_safe ||
      !cut_has_coeff ||
      f0_eff <= 1.0e-6 ||
      max_g <= 1.0e-12 ||
      max_g / min_g > 1.0e12 {
      // Unsafe or ill-conditioned cut; abort cut loop and fall back to exact B&B
      break
    }

    num_cols = num_cols + 1
    rhs_col = num_cols
    for i = 0; i <= m; i = i + 1 {
      let row_i = tab[i]
      let old_rhs = row_i[new_col]
      row_i[new_col] = 0.0
      row_i.push(old_rhs)
    }
    cut_row[new_col] = 1.0
    cut_row[rhs_col] = -f0_eff

    // Insert cut_row before obj_row
    let saved_obj = tab[m]
    tab[m] = cut_row
    tab.push(saved_obj)
    basis.push(new_col)
    col_active.push(true)
    is_int_col.push(false)
    col_ubs.push(
      if cut_slack_finite && is_finite_double(cut_slack_ub) && cut_slack_ub >= 0.0 {
        cut_slack_ub
      } else {
        1.0e300
      },
    )
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
    if dst == MaxIterations {
      return (
        Solution::empty(MaxIterations, total_iters, 1),
        {
          cuts_added,
          root_lp_obj,
          cut_lp_obj: (-tab[m][rhs_col] + obj_offset) * sense_sign,
          dual_pivots,
          solved_at_root: false,
        },
      )
    }
    if dst != Optimal {
      break
    }
  }

  let cut_lp_min = -tab[m][rhs_col] + obj_offset
  let cut_lp_obj = cut_lp_min * sense_sign

  if is_integral {
    let raw_primal = Array::make(n_orig, 0.0)
    for j = 0; j < n_orig; j = j + 1 {
      raw_primal[j] = lbs[j]
    }
    for i = 0; i < m; i = i + 1 {
      let bc = basis[i]
      if bc >= 0 && bc < n_std {
        raw_primal[bc] = lbs[bc] + tab[i][rhs_col]
      }
    }
    match verify_or_repair_int_candidate(model, raw_primal, config) {
      Some(cand) => {
        let cand_min = cand.objective * sense_sign
        // Only declare root optimality if the verified candidate actually attains the LP bound!
        if candidate_fathoms_node(cand_min, cut_lp_min, config.mip_gap_tol) {
          let root_sol : Solution = {
            ..cand,
            status: Optimal,
            iterations: total_iters,
            nodes_explored: 1,
            mip_gap: compute_mip_gap(cand_min, cut_lp_min),
          }
          return (
            root_sol,
            {
              cuts_added,
              root_lp_obj,
              cut_lp_obj: cand.objective,
              dual_pivots,
              solved_at_root: true,
            },
          )
        }
      }
      None => ()
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
