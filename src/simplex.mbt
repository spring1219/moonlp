///|
/// Two-Phase Simplex LP Solver with:
/// - Arbitrary variable bounds (including negative lower bounds & free variables)
/// - Dantzig steepest-reduced-cost pricing + automatic Bland's anti-cycling rule
/// - Exact primal solution, dual shadow prices, and reduced cost extraction
/// - Infeasibility & unboundedness certificates

///|
let inf_threshold : Double = 1.0e200

///|
priv enum VarMapping {
  Shifted(Int32, Double) // x = u + lb
  NegShifted(Int32, Double) // x = ub - v
  Free(Int32, Int32) // x = u_pos - u_neg
  Fixed(Double) // x = lb == ub
}

///|
/// Zero-allocation sparse-aware Gauss-Jordan pivot using reusable scratch buffers.
fn pivot_tableau_buf(
  tab : Array[Array[Double]],
  basis : Array[Int32],
  m : Int32,
  rhs_col : Int32,
  pivot_row : Int32,
  pivot_col : Int32,
  nz_cols : Array[Int32],
  nz_vals : Array[Double],
) -> Unit {
  let pv = tab[pivot_row][pivot_col]
  let inv_pv = 1.0 / pv
  let prow = tab[pivot_row]
  nz_cols.clear()
  nz_vals.clear()
  for j = 0; j <= rhs_col; j = j + 1 {
    let v = prow[j]
    if v != 0.0 {
      let scaled = v * inv_pv
      prow[j] = scaled
      if j != pivot_col {
        nz_cols.push(j)
        nz_vals.push(scaled)
      }
    }
  }
  prow[pivot_col] = 1.0

  let nz_len = nz_cols.length()
  for i = 0; i <= m; i = i + 1 {
    if i != pivot_row {
      let factor = tab[i][pivot_col]
      if factor != 0.0 {
        let row_i = tab[i]
        for k = 0; k < nz_len; k = k + 1 {
          let j = nz_cols[k]
          row_i[j] = row_i[j] - factor * nz_vals[k]
        }
        row_i[pivot_col] = 0.0
      }
    }
  }
  basis[pivot_row] = pivot_col
}

///|
fn pivot_tableau(
  tab : Array[Array[Double]],
  basis : Array[Int32],
  m : Int32,
  rhs_col : Int32,
  pivot_row : Int32,
  pivot_col : Int32,
) -> Unit {
  let nz_cols : Array[Int32] = []
  let nz_vals : Array[Double] = []
  pivot_tableau_buf(
    tab, basis, m, rhs_col, pivot_row, pivot_col, nz_cols, nz_vals,
  )
}

///|
/// Run simplex iterations on the current tableau objective row `m`.
/// Returns (status, iterations_used).
fn run_simplex_phase(
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
    let obj_row = tab[m]

    // 1. Choose entering column
    let mut enter_col = -1
    let mut best_rc = -config.tol
    for j = 0; j < num_cols; j = j + 1 {
      if !col_active[j] {
        continue
      }
      let rc = obj_row[j]
      if rc < -config.tol {
        if use_bland {
          enter_col = j
          break
        } else if rc < best_rc {
          best_rc = rc
          enter_col = j
        }
      }
    }

    if enter_col == -1 {
      return (Optimal, iters)
    }

    // 2. Minimum ratio test with tight relative tie tolerance (1e-12) and Harris pivot guard
    let pivot_tol = 1.0e-12
    let mut leave_row = -1
    let mut best_ratio = 1.0e308
    let mut best_pivot = 0.0
    for i = 0; i < m; i = i + 1 {
      let a_ie = tab[i][enter_col]
      if a_ie > pivot_tol {
        let rhs_val = if tab[i][rhs_col] < 0.0 && tab[i][rhs_col] > -config.tol {
          0.0
        } else {
          tab[i][rhs_col]
        }
        let ratio = rhs_val / a_ie
        let tie_tol = 1.0e-12 * (1.0 + best_ratio.abs())
        if leave_row == -1 || ratio < best_ratio - tie_tol {
          best_ratio = ratio
          best_pivot = a_ie
          leave_row = i
        } else if (ratio - best_ratio).abs() <= tie_tol {
          if use_bland {
            if basis[i] < basis[leave_row] {
              leave_row = i
              best_pivot = a_ie
            }
          } else if a_ie > best_pivot {
            leave_row = i
            best_pivot = a_ie
          }
        }
      }
    }

    if leave_row == -1 {
      return (Unbounded, iters)
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
/// Check whether constraint coefficients exhibit extreme scales (< 1e-4, > 1e4, or ratio > 1e3).
fn needs_auto_scaling(model : Model) -> Bool {
  if model.constraints.length() == 0 || model.vars.length() == 0 {
    return false
  }
  let mut min_v = 1.0e308
  let mut max_v = 0.0
  for i = 0; i < model.constraints.length(); i = i + 1 {
    let terms = model.constraints[i].terms
    for k = 0; k < terms.length(); k = k + 1 {
      let a = terms[k].coeff.abs()
      if a > 0.0 {
        if a < min_v {
          min_v = a
        }
        if a > max_v {
          max_v = a
        }
      }
    }
  }
  if max_v <= 0.0 {
    return false
  }
  min_v < 1.0e-4 || max_v > 1.0e4 || max_v / min_v > 1.0e3
}

///|
/// Solve the LP relaxation of a Model using Two-Phase Simplex with automatic
/// Ruiz equilibration scaling on poorly-scaled matrices.
pub fn solve_lp(
  model : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  if !model.has_integer_vars() && needs_auto_scaling(model) {
    let (sol, _) = model.solve_scaled(passes=6, config~)
    if sol.status == Optimal && model.verify_solution(sol, tol=1.0e-4) {
      return sol
    }
    if sol.status != Optimal {
      return sol
    }
  }
  let raw = solve_lp_unscaled(model, config)
  if raw.status == Optimal &&
    !model.has_integer_vars() &&
    !model.verify_solution(raw, tol=1.0e-4) {
    let (scaled_sol, _) = model.solve_scaled(passes=6, config~)
    if scaled_sol.status == Optimal &&
      model.verify_solution(scaled_sol, tol=1.0e-4) {
      return scaled_sol
    }
  }
  raw
}

///|
fn solve_lp_unscaled(
  model : Model,
  config : SolverConfig,
) -> Solution {
  let n_orig = model.vars.length()
  let m_user = model.constraints.length()
  let sense_sign = match model.sense {
    Minimize => 1.0
    Maximize => -1.0
  }

  // Step 1: Variable bound checks and standard non-negative variable mapping
  let mappings : Array[VarMapping] = []
  let std_obj : Array[Double] = []
  let mut obj_offset = model.obj_constant * sense_sign
  let ub_rows : Array[(Int32, Double)] = [] // (std_col, upper_bound_delta)

  for j = 0; j < n_orig; j = j + 1 {
    let v = model.vars[j]
    if v.lb > v.ub + config.tol {
      return Solution::empty(Infeasible, 0, 1)
    }
    let c_j = v.obj * sense_sign
    let has_lb = v.lb > -inf_threshold
    let has_ub = v.ub < inf_threshold

    if has_lb && has_ub && (v.ub - v.lb).abs() <= config.tol {
      obj_offset = obj_offset + c_j * v.lb
      mappings.push(Fixed(v.lb))
    } else if has_lb {
      let col = std_obj.length()
      std_obj.push(c_j)
      obj_offset = obj_offset + c_j * v.lb
      mappings.push(Shifted(col, v.lb))
      if has_ub {
        ub_rows.push((col, v.ub - v.lb))
      }
    } else if has_ub {
      let col = std_obj.length()
      std_obj.push(-c_j)
      obj_offset = obj_offset + c_j * v.ub
      mappings.push(NegShifted(col, v.ub))
    } else {
      let col_p = std_obj.length()
      std_obj.push(c_j)
      let col_n = std_obj.length()
      std_obj.push(-c_j)
      mappings.push(Free(col_p, col_n))
    }
  }

  let n_std = std_obj.length()
  let m_total = m_user + ub_rows.length()

  // Handle unconstrained trivial case (0 rows)
  if m_total == 0 {
    for j = 0; j < n_std; j = j + 1 {
      if std_obj[j] < -config.tol {
        return Solution::empty(Unbounded, 0, 1)
      }
    }
    let primal : Array[Double] = []
    let reduced_costs : Array[Double] = []
    for j = 0; j < n_orig; j = j + 1 {
      match mappings[j] {
        Shifted(_, lb) => {
          primal.push(lb)
          reduced_costs.push(model.vars[j].obj)
        }
        NegShifted(_, ub) => {
          primal.push(ub)
          reduced_costs.push(model.vars[j].obj)
        }
        Free(_, _) => {
          primal.push(0.0)
          reduced_costs.push(model.vars[j].obj)
        }
        Fixed(val) => {
          primal.push(val)
          reduced_costs.push(0.0)
        }
      }
    }
    return {
      status: Optimal,
      objective: obj_offset * sense_sign,
      primal,
      dual: [],
      reduced_costs,
      iterations: 0,
      nodes_explored: 1,
      mip_gap: 0.0,
    }
  }

  // Step 2: Build transformed constraint rows (all RHS >= 0)
  let row_coeffs : Array[Array[Double]] = []
  let row_rhs : Array[Double] = []
  let row_rels : Array[Relation] = []
  let row_signs : Array[Double] = []

  for i = 0; i < m_user; i = i + 1 {
    let c = model.constraints[i]
    let coeffs = Array::make(n_std, 0.0)
    let mut rhs = c.rhs
    for k = 0; k < c.terms.length(); k = k + 1 {
      let t = c.terms[k]
      match mappings[t.var_id.index] {
        Shifted(col, lb) => {
          coeffs[col] = coeffs[col] + t.coeff
          rhs = rhs - t.coeff * lb
        }
        NegShifted(col, ub) => {
          coeffs[col] = coeffs[col] - t.coeff
          rhs = rhs - t.coeff * ub
        }
        Free(col_p, col_n) => {
          coeffs[col_p] = coeffs[col_p] + t.coeff
          coeffs[col_n] = coeffs[col_n] - t.coeff
        }
        Fixed(val) => rhs = rhs - t.coeff * val
      }
    }
    if rhs < 0.0 {
      for j = 0; j < n_std; j = j + 1 {
        coeffs[j] = -coeffs[j]
      }
      rhs = -rhs
      let flipped = match c.relation {
        Le => Ge
        Ge => Le
        Eq => Eq
      }
      row_coeffs.push(coeffs)
      row_rhs.push(rhs)
      row_rels.push(flipped)
      row_signs.push(-1.0)
    } else {
      row_coeffs.push(coeffs)
      row_rhs.push(rhs)
      row_rels.push(c.relation)
      row_signs.push(1.0)
    }
  }

  // Upper bound rows: u_col <= ub - lb
  for i = 0; i < ub_rows.length(); i = i + 1 {
    let (col, delta) = ub_rows[i]
    let coeffs = Array::make(n_std, 0.0)
    coeffs[col] = 1.0
    row_coeffs.push(coeffs)
    row_rhs.push(delta)
    row_rels.push(Le)
    row_signs.push(1.0)
  }

  // Step 3: Allocate slack, surplus, and artificial columns
  let mut next_col = n_std
  let row_slack_col : Array[Int32] = Array::make(m_total, -1)
  let row_art_col : Array[Int32] = Array::make(m_total, -1)
  let init_unit_col : Array[Int32] = Array::make(m_total, -1)
  let is_artificial : Array[Bool] = Array::make(n_std, false)

  for i = 0; i < m_total; i = i + 1 {
    match row_rels[i] {
      Le => {
        row_slack_col[i] = next_col
        init_unit_col[i] = next_col
        is_artificial.push(false)
        next_col = next_col + 1
      }
      Ge => {
        row_slack_col[i] = next_col // surplus (-1)
        is_artificial.push(false)
        next_col = next_col + 1
        row_art_col[i] = next_col // artificial (+1)
        init_unit_col[i] = next_col
        is_artificial.push(true)
        next_col = next_col + 1
      }
      Eq => {
        row_art_col[i] = next_col // artificial (+1)
        init_unit_col[i] = next_col
        is_artificial.push(true)
        next_col = next_col + 1
      }
    }
  }

  let num_cols = next_col
  let rhs_col = num_cols
  let tab : Array[Array[Double]] = []
  let basis : Array[Int32] = Array::make(m_total, -1)
  let mut has_art = false

  for i = 0; i < m_total; i = i + 1 {
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
  tab.push(Array::make(num_cols + 1, 0.0)) // objective row m_total

  let col_active = Array::make(num_cols, true)
  let mut total_iters = 0

  // Step 4: Phase I (if artificial variables exist)
  if has_art {
    let obj_row = tab[m_total]
    for j = 0; j < num_cols; j = j + 1 {
      obj_row[j] = if is_artificial[j] { 1.0 } else { 0.0 }
    }
    obj_row[rhs_col] = 0.0
    for i = 0; i < m_total; i = i + 1 {
      if is_artificial[basis[i]] {
        for j = 0; j <= rhs_col; j = j + 1 {
          obj_row[j] = obj_row[j] - tab[i][j]
        }
      }
    }

    let (p1_status, p1_iters) = run_simplex_phase(
      tab, basis, col_active, m_total, num_cols, config,
    )
    total_iters = total_iters + p1_iters
    if p1_status != Optimal {
      return Solution::empty(p1_status, total_iters, 1)
    }

    let p1_obj = -tab[m_total][rhs_col]
    if p1_obj > config.tol * 100.0 {
      return Solution::empty(Infeasible, total_iters, 1)
    }

    // Drive degenerate artificial variables out of basis where possible
    for i = 0; i < m_total; i = i + 1 {
      if is_artificial[basis[i]] {
        let mut pivot_col = -1
        for j = 0; j < num_cols; j = j + 1 {
          if !is_artificial[j] && tab[i][j].abs() > config.tol {
            pivot_col = j
            break
          }
        }
        if pivot_col != -1 {
          pivot_tableau(tab, basis, m_total, rhs_col, i, pivot_col)
        }
      }
    }

    // Deactivate all artificial columns for Phase II (columns remain updated for dual extraction)
    for j = 0; j < num_cols; j = j + 1 {
      if is_artificial[j] {
        col_active[j] = false
      }
    }
  }

  // Step 5: Phase II (optimize true objective)
  let obj_row = tab[m_total]
  for j = 0; j <= rhs_col; j = j + 1 {
    obj_row[j] = 0.0
  }
  for j = 0; j < n_std; j = j + 1 {
    obj_row[j] = std_obj[j]
  }
  // Eliminate basic variables from objective row
  for i = 0; i < m_total; i = i + 1 {
    let bcol = basis[i]
    let cb = if bcol < n_std { std_obj[bcol] } else { 0.0 }
    if cb != 0.0 {
      for j = 0; j <= rhs_col; j = j + 1 {
        obj_row[j] = obj_row[j] - cb * tab[i][j]
      }
    }
  }

  let (p2_status, p2_iters) = run_simplex_phase(
    tab, basis, col_active, m_total, num_cols, config,
  )
  total_iters = total_iters + p2_iters
  if p2_status != Optimal {
    return Solution::empty(p2_status, total_iters, 1)
  }

  // Step 6: Reconstruct primal solution, dual shadow prices, and reduced costs
  let std_x = Array::make(n_std, 0.0)
  for i = 0; i < m_total; i = i + 1 {
    let bcol = basis[i]
    if bcol >= 0 && bcol < n_std {
      let val = tab[i][rhs_col]
      std_x[bcol] = if val.abs() < config.tol { 0.0 } else { val }
    }
  }

  let primal : Array[Double] = []
  let reduced_costs : Array[Double] = []
  for j = 0; j < n_orig; j = j + 1 {
    match mappings[j] {
      Shifted(col, lb) => {
        primal.push(std_x[col] + lb)
        reduced_costs.push(tab[m_total][col] * sense_sign)
      }
      NegShifted(col, ub) => {
        primal.push(ub - std_x[col])
        reduced_costs.push(-tab[m_total][col] * sense_sign)
      }
      Free(col_p, col_n) => {
        primal.push(std_x[col_p] - std_x[col_n])
        reduced_costs.push(tab[m_total][col_p] * sense_sign)
      }
      Fixed(val) => {
        primal.push(val)
        reduced_costs.push(0.0)
      }
    }
  }

  // Extract dual shadow prices for original user constraints via initial unit columns
  let dual : Array[Double] = []
  for i = 0; i < m_user; i = i + 1 {
    let ucol = init_unit_col[i]
    let y_std = -tab[m_total][ucol]
    let y_user = y_std * row_signs[i] * sense_sign
    dual.push(if y_user.abs() < config.tol { 0.0 } else { y_user })
  }

  let min_obj = -tab[m_total][rhs_col] + obj_offset
  let final_obj = min_obj * sense_sign

  {
    status: Optimal,
    objective: final_obj,
    primal,
    dual,
    reduced_costs,
    iterations: total_iters,
    nodes_explored: 1,
    mip_gap: 0.0,
  }
}
