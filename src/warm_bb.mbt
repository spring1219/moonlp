///|
/// Dual Simplex Warm-Start Branch-and-Bound MILP Solver (v0.4.0):
/// Instead of rebuilding the LP and restarting Phase I from scratch at every B&B node,
/// each child node inherits its parent's optimal simplex tableau, appends a single
/// bound cut row (eliminating the branching basic variable in O(n)), and re-optimizes
/// directly via Dual Simplex in 1-3 pivots (zero Phase I restarts!).

///|
priv struct WarmBBNode {
  tab : Array[Array[Double]]
  basis : Array[Int32]
  col_active : Array[Bool]
  m : Int32
  num_cols : Int32
  bound : Double // LP relaxation lower bound in minimization space
}

///|
fn clone_tableau(tab : Array[Array[Double]]) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = []
  for i = 0; i < tab.length(); i = i + 1 {
    out.push(tab[i].copy())
  }
  out
}

///|
fn pop_best_warm_node(queue : Array[WarmBBNode]) -> WarmBBNode {
  let mut best_i = 0
  let mut best_b = queue[0].bound
  for i = 1; i < queue.length(); i = i + 1 {
    if queue[i].bound < best_b {
      best_b = queue[i].bound
      best_i = i
    }
  }
  let node = queue[best_i]
  let last = queue.length() - 1
  if best_i != last {
    queue[best_i] = queue[last]
  }
  let _ = queue.pop()
  node
}

///|
/// Append a single bound cut row to a cloned optimal tableau and keep it dual-feasible:
/// Given basic row `r`: x_k + sum_{j != k} a_{rj} x_j = b_r (fractional):
/// - Upper branch (x_k <= fl): -sum_{j != k} a_{rj} x_j + s_new = fl - b_r < 0
/// - Lower branch (x_k >= cl): +sum_{j != k} a_{rj} x_j + s_new = b_r - cl < 0
fn make_warm_child(
  parent : WarmBBNode,
  basic_row : Int32,
  is_upper_branch : Bool,
  bound_rhs : Double,
) -> WarmBBNode {
  let tab = clone_tableau(parent.tab)
  let basis = parent.basis.copy()
  let col_active = parent.col_active.copy()
  let m = parent.m
  let old_num_cols = parent.num_cols
  let new_col = old_num_cols
  let new_num_cols = old_num_cols + 1
  let new_rhs_col = new_num_cols

  // Expand all existing rows + objective row by 1 column before RHS
  for i = 0; i <= m; i = i + 1 {
    let row_i = tab[i]
    let old_rhs = row_i[old_num_cols]
    row_i[old_num_cols] = 0.0
    row_i.push(old_rhs)
  }

  let bcol = basis[basic_row]
  let b_r = tab[basic_row][new_rhs_col]
  let cut_row = Array::make(new_num_cols + 1, 0.0)
  let sign = if is_upper_branch { -1.0 } else { 1.0 }
  for j = 0; j < old_num_cols; j = j + 1 {
    if j != bcol {
      let a = tab[basic_row][j]
      if a != 0.0 {
        cut_row[j] = sign * a
      }
    }
  }
  cut_row[new_col] = 1.0
  cut_row[new_rhs_col] = if is_upper_branch {
    bound_rhs - b_r
  } else {
    b_r - bound_rhs
  }

  let saved_obj = tab[m]
  tab[m] = cut_row
  tab.push(saved_obj)
  basis.push(new_col)
  col_active.push(true)

  {
    tab,
    basis,
    col_active,
    m: m + 1,
    num_cols: new_num_cols,
    bound: parent.bound,
  }
}

///|
/// Solve a MILP using Best-Bound Branch-and-Bound with Dual Simplex Warm-Start.
/// Falls back cleanly to `solve_milp` if any variable has `-inf` lower bound.
pub fn solve_milp_warm(
  model : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  if !model.has_integer_vars() {
    return solve_lp(model, config~)
  }
  if has_extreme_scale(model) {
    return solve_milp(model, config~)
  }

  let n_orig = model.vars.length()
  let m_user = model.constraints.length()
  let sense_sign = match model.sense {
    Minimize => 1.0
    Maximize => -1.0
  }

  // Step 1: Standardize non-negative shifted variables
  let lbs : Array[Double] = []
  let is_fixed : Array[Bool] = []
  let std_obj : Array[Double] = []
  let is_int_col : Array[Bool] = []
  let ub_rows : Array[(Int32, Double)] = []
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
      return Solution::empty(Infeasible, 0, 1)
    }
    if eff_lb <= -inf_threshold {
      return solve_milp(model, config~)
    }
    let fixed_j = eff_ub < inf_threshold && eff_lb == eff_ub
    let c_j = v.obj * sense_sign
    lbs.push(eff_lb)
    is_fixed.push(fixed_j)
    std_obj.push(if fixed_j { 0.0 } else { c_j })
    is_int_col.push(is_integer_var(v.kind))
    obj_offset = obj_offset + c_j * eff_lb
    if !fixed_j && eff_ub < inf_threshold {
      ub_rows.push((j, eff_ub - eff_lb))
    }
  }

  let n_std = n_orig
  let m_total = m_user + ub_rows.length()
  if m_total == 0 {
    return solve_milp(model, config~)
  }

  // Step 2: Build transformed constraint rows (RHS >= 0)
  let row_coeffs : Array[Array[Double]] = []
  let row_rhs : Array[Double] = []
  let row_rels : Array[Relation] = []

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
    row_coeffs.push(coeffs)
    row_rhs.push(rhs)
    row_rels.push(rel)
  }

  for i = 0; i < ub_rows.length(); i = i + 1 {
    let (col, delta) = ub_rows[i]
    let coeffs = Array::make(n_std, 0.0)
    coeffs[col] = 1.0
    row_coeffs.push(coeffs)
    row_rhs.push(delta)
    row_rels.push(Le)
  }

  // Step 3: Allocate slack, surplus, and artificial columns
  let mut num_cols = n_std
  let row_slack_col = Array::make(m_total, -1)
  let row_art_col = Array::make(m_total, -1)
  let is_artificial = Array::make(n_std, false)

  for i = 0; i < m_total; i = i + 1 {
    match row_rels[i] {
      Le => {
        row_slack_col[i] = num_cols
        is_artificial.push(false)
        num_cols = num_cols + 1
      }
      Ge => {
        row_slack_col[i] = num_cols
        is_artificial.push(false)
        num_cols = num_cols + 1
        row_art_col[i] = num_cols
        is_artificial.push(true)
        num_cols = num_cols + 1
      }
      Eq => {
        row_art_col[i] = num_cols
        is_artificial.push(true)
        num_cols = num_cols + 1
      }
    }
  }

  let rhs_col = num_cols
  let tab : Array[Array[Double]] = []
  let basis = Array::make(m_total, -1)
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
  tab.push(Array::make(num_cols + 1, 0.0))

  let col_active = Array::make(num_cols, true)
  for j = 0; j < n_std; j = j + 1 {
    if is_fixed[j] {
      col_active[j] = false
    }
  }
  let mut total_iters = 0

  // Solve Root LP ONCE (Phase I + Phase II)
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
    let (p1_st, p1_it) = run_simplex_phase(
      tab, basis, col_active, m_total, num_cols, config,
    )
    total_iters = total_iters + p1_it
    if p1_st != Optimal {
      return Solution::empty(p1_st, total_iters, 1)
    }
    if -tab[m_total][rhs_col] > config.tol {
      return Solution::empty(Infeasible, total_iters, 1)
    }
    drive_out_artificials(
      tab, basis, col_active, is_artificial, m_total, num_cols, config,
    )
    for j = 0; j < num_cols; j = j + 1 {
      if is_artificial[j] {
        col_active[j] = false
      }
    }
  }

  let obj_row = tab[m_total]
  for j = 0; j <= rhs_col; j = j + 1 {
    obj_row[j] = 0.0
  }
  for j = 0; j < n_std; j = j + 1 {
    obj_row[j] = std_obj[j]
  }
  for i = 0; i < m_total; i = i + 1 {
    let bc = basis[i]
    let cb = if bc < n_std { std_obj[bc] } else { 0.0 }
    if cb != 0.0 {
      for j = 0; j <= rhs_col; j = j + 1 {
        obj_row[j] = obj_row[j] - cb * tab[i][j]
      }
    }
  }

  let (p2_st, p2_it) = run_simplex_phase(
    tab, basis, col_active, m_total, num_cols, config,
  )
  total_iters = total_iters + p2_it
  if p2_st == Unbounded {
    return Solution::empty(
      resolve_unbounded_milp_status(model, config),
      total_iters,
      1,
    )
  }
  if p2_st != Optimal {
    return Solution::empty(p2_st, total_iters, 1)
  }

  // Step 4: Best-Bound Branch-and-Bound with Dual Simplex Warm-Start
  let root_bound = -tab[m_total][rhs_col] + obj_offset
  let queue : Array[WarmBBNode] = [
    { tab, basis, col_active, m: m_total, num_cols, bound: root_bound, },
  ]
  let mut best_min_obj = 1.0e308
  let mut best_pruned_bound = 1.0e308
  let mut best_sol : Solution? = None
  let mut nodes_explored = 0
  let mut hit_max_iters = false
  let mut has_unproved = false

  while queue.length() > 0 && nodes_explored < config.max_nodes {
    let node = pop_best_warm_node(queue)
    if best_sol is Some(_) &&
      can_prune_by_bound(node.bound, best_min_obj, config.mip_gap_tol) {
      if node.bound < best_pruned_bound {
        best_pruned_bound = node.bound
      }
      continue
    }

    // Re-optimize child node with Dual Simplex (root is already optimal, 0 pivots!)
    let (dst, dit) = run_dual_simplex(
      node.tab,
      node.basis,
      node.col_active,
      node.m,
      node.num_cols,
      config,
    )
    nodes_explored = nodes_explored + 1
    total_iters = total_iters + dit

    if dst == MaxIterations {
      hit_max_iters = true
      has_unproved = true
      continue
    }
    if dst != Optimal {
      continue
    }

    let min_obj = -node.tab[node.m][node.num_cols] + obj_offset
    if best_sol is Some(_) &&
      can_prune_by_bound(min_obj, best_min_obj, config.mip_gap_tol) {
      if min_obj < best_pruned_bound {
        best_pruned_bound = min_obj
      }
      continue
    }

    // Pick most-fractional basic integer variable
    let mut branch_row = -1
    let mut best_dist = 0.0
    for i = 0; i < node.m; i = i + 1 {
      let bc = node.basis[i]
      if bc >= 0 && bc < n_std && is_int_col[bc] {
        let val = node.tab[i][node.num_cols]
        let dist = (val - d_round(val)).abs()
        if dist > config.int_tol && dist > best_dist {
          best_dist = dist
          branch_row = i
        }
      }
    }

    if branch_row == -1 {
      let raw_primal = Array::make(n_orig, 0.0)
      for j = 0; j < n_orig; j = j + 1 {
        raw_primal[j] = lbs[j]
      }
      for i = 0; i < node.m; i = i + 1 {
        let bc = node.basis[i]
        if bc >= 0 && bc < n_std {
          raw_primal[bc] = lbs[bc] + node.tab[i][node.num_cols]
        }
      }
      let mut fathomed = false
      match verify_or_repair_int_candidate(model, raw_primal, config) {
        Some(cand) => {
          let cand_min = cand.objective * sense_sign
          if cand_min < best_min_obj {
            best_min_obj = cand_min
            best_sol = Some({
              ..cand,
              iterations: total_iters,
              nodes_explored,
              mip_gap: 0.0,
            })
          }
          if candidate_fathoms_node(cand_min, min_obj, config.mip_gap_tol) {
            fathomed = true
            if min_obj < best_pruned_bound {
              best_pruned_bound = min_obj
            }
          }
        }
        None => ()
      }
      if !fathomed {
        // Sub-tolerance fractionality or degraded repair: branch strictly above dual-simplex tol!
        let mut sub_dist = 0.0
        for i = 0; i < node.m; i = i + 1 {
          let bc = node.basis[i]
          if bc >= 0 && bc < n_std && is_int_col[bc] {
            let val = node.tab[i][node.num_cols]
            let dist = (val - d_round(val)).abs()
            if dist > config.tol && dist > sub_dist {
              sub_dist = dist
              branch_row = i
            }
          }
        }
        if branch_row == -1 {
          has_unproved = true
        }
      }
    }

    if branch_row != -1 {
      let val = node.tab[branch_row][node.num_cols]
      let fl = d_floor(val)
      let cl = d_ceil(val)
      let updated_parent : WarmBBNode = { ..node, bound: min_obj, }
      queue.push(make_warm_child(updated_parent, branch_row, true, fl))
      queue.push(make_warm_child(updated_parent, branch_row, false, cl))
    }
  }

  match best_sol {
    Some(sol) =>
      if queue.length() == 0 && !has_unproved {
        let eff_bound = if best_pruned_bound < best_min_obj {
          best_pruned_bound
        } else {
          best_min_obj
        }
        let final_gap = compute_mip_gap(best_min_obj, eff_bound)
        {
          ..sol,
          status: Optimal,
          iterations: total_iters,
          nodes_explored,
          mip_gap: final_gap,
        }
      } else {
        let mut best_rem = best_pruned_bound
        for i = 0; i < queue.length(); i = i + 1 {
          if queue[i].bound < best_rem {
            best_rem = queue[i].bound
          }
        }
        if best_rem > best_min_obj {
          best_rem = best_min_obj
        }
        let gap = if has_unproved && queue.length() == 0 {
          1.0
        } else {
          compute_mip_gap(best_min_obj, best_rem)
        }
        let st = if !has_unproved && gap <= config.mip_gap_tol {
          Optimal
        } else if hit_max_iters && nodes_explored < config.max_nodes {
          MaxIterations
        } else if nodes_explored >= config.max_nodes {
          NodeLimit
        } else {
          NumericalError
        }
        { ..sol, status: st, iterations: total_iters, nodes_explored, mip_gap: gap, }
      }
    None =>
      if hit_max_iters {
        Solution::empty(MaxIterations, total_iters, nodes_explored)
      } else if nodes_explored >= config.max_nodes {
        Solution::empty(NodeLimit, total_iters, nodes_explored)
      } else if has_unproved {
        Solution::empty(NumericalError, total_iters, nodes_explored)
      } else {
        Solution::empty(Infeasible, total_iters, nodes_explored)
      }
  }
}
