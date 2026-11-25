///|
/// Unified Flagship Branch-and-Cut MILP Solver (v0.5.0):
/// Combines:
///   1. Root Gomory Mixed-Integer (GMI) Cutting Planes + Dual Simplex re-optimization
///   2. Hybrid DFS-Plunge + Best-Bound node selection (finds incumbents fast, prunes aggressively)
///   3. Dual Simplex Warm-Start branching inheriting root cuts into every child node (zero Phase I)
///   4. Objective-weighted fractional branching + nearest-integer dive ordering
///   5. Primal Rounding Heuristic to capture feasible incumbents early

///|
#derive(Debug)
pub(all) struct BranchAndCutStats {
  root_lp_obj : Double
  post_cut_lp_obj : Double
  cuts_added : Int32
  heuristic_hits : Int32
  nodes_explored : Int32
  total_pivots : Int32
}

///|
priv struct BCNode {
  tab : Array[Array[Double]]
  basis : Array[Int32]
  col_active : Array[Bool]
  m : Int32
  num_cols : Int32
  bound : Double
  depth : Int32
}

///|
fn pop_hybrid_bc_node(
  stack : Array[BCNode],
  has_incumbent : Bool,
  nodes_explored : Int32,
) -> BCNode {
  // Plunge DFS (LIFO) until an incumbent is found, then alternate 3 DFS plunges + 1 Best-Bound
  if !has_incumbent || nodes_explored % 4 != 0 {
    match stack.pop() {
      Some(n) => n
      None => abort("unreachable: empty stack")
    }
  } else {
    let mut best_i = 0
    let mut best_b = stack[0].bound
    for i = 1; i < stack.length(); i = i + 1 {
      if stack[i].bound < best_b {
        best_b = stack[i].bound
        best_i = i
      }
    }
    let node = stack[best_i]
    let last = stack.length() - 1
    if best_i != last {
      stack[best_i] = stack[last]
    }
    let _ = stack.pop()
    node
  }
}

///|
fn make_bc_child(
  parent : BCNode,
  basic_row : Int32,
  is_upper_branch : Bool,
  bound_rhs : Double,
  parent_bound : Double,
) -> BCNode {
  let tab = clone_tableau(parent.tab)
  let basis = parent.basis.copy()
  let col_active = parent.col_active.copy()
  let m = parent.m
  let old_num_cols = parent.num_cols
  let new_col = old_num_cols
  let new_num_cols = old_num_cols + 1
  let new_rhs_col = new_num_cols

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
    bound: parent_bound,
    depth: parent.depth + 1,
  }
}

///|
/// Try simple nearest-integer rounding heuristic on a fractional LP node solution.
fn try_rounding_heuristic(
  model : Model,
  lbs : Array[Double],
  tab : Array[Array[Double]],
  basis : Array[Int32],
  m : Int32,
  n_std : Int32,
  rhs_col : Int32,
  is_int_col : Array[Bool],
  sense_sign : Double,
  best_min_obj : Double,
  config : SolverConfig,
) -> Solution? {
  let n_orig = model.vars.length()
  let cand_x = Array::make(n_orig, 0.0)
  for j = 0; j < n_orig; j = j + 1 {
    cand_x[j] = lbs[j]
  }
  for i = 0; i < m; i = i + 1 {
    let bc = basis[i]
    if bc >= 0 && bc < n_std {
      let raw = lbs[bc] + tab[i][rhs_col]
      cand_x[bc] = if is_int_col[bc] { d_round(raw) } else { raw }
    }
  }
  let mut exact_obj = model.obj_constant
  for j = 0; j < n_orig; j = j + 1 {
    exact_obj = exact_obj + model.vars[j].obj * cand_x[j]
  }
  let h_tol = 1.0e-12 * best_min_obj.abs() + 1.0e-15
  if exact_obj * sense_sign >= best_min_obj - h_tol {
    return None
  }
  let sol : Solution = {
    status: Optimal,
    objective: exact_obj,
    primal: cand_x,
    dual: [],
    reduced_costs: [],
    iterations: 0,
    nodes_explored: 0,
    mip_gap: 0.0,
  }
  if model.verify_solution(sol, tol=config.tol) {
    Some(sol)
  } else {
    None
  }
}

///|
/// Unified Flagship Branch-and-Cut MILP Solver.
pub fn solve_milp_branch_and_cut(
  model : Model,
  max_root_cuts~ : Int32 = 8,
  config~ : SolverConfig = SolverConfig::default(),
) -> (Solution, BranchAndCutStats) {
  let empty_stats : BranchAndCutStats = {
    root_lp_obj: 0.0,
    post_cut_lp_obj: 0.0,
    cuts_added: 0,
    heuristic_hits: 0,
    nodes_explored: 0,
    total_pivots: 0,
  }
  if !model.has_integer_vars() {
    let sol = solve_lp(model, config~)
    return (
      sol,
      {
        ..empty_stats,
        root_lp_obj: sol.objective,
        post_cut_lp_obj: sol.objective,
        nodes_explored: 1,
        total_pivots: sol.iterations,
      },
    )
  }
  if has_extreme_scale(model) {
    let sol = solve_milp(model, config~)
    return (
      sol,
      {
        ..empty_stats,
        nodes_explored: sol.nodes_explored,
        total_pivots: sol.iterations,
      },
    )
  }

  let n_orig = model.vars.length()
  let m_user = model.constraints.length()
  let sense_sign = match model.sense {
    Minimize => 1.0
    Maximize => -1.0
  }

  let lbs : Array[Double] = []
  let init_ubs : Array[Double] = []
  let is_fixed : Array[Bool] = []
  let std_obj : Array[Double] = []
  let is_int_col : Array[Bool] = []
  let ub_rows : Array[(Int32, Double, Bool)] = []
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
      return (Solution::empty(Infeasible, 0, 1), empty_stats)
    }
    if eff_lb <= -inf_threshold {
      let sol = solve_milp(model, config~)
      return (
        sol,
        {
          ..empty_stats,
          nodes_explored: sol.nodes_explored,
          total_pivots: sol.iterations,
        },
      )
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
    return (sol, empty_stats)
  }

  // Build standard rows & detect integral slacks
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
  let basis = Array::make(m, -1)
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
  let mut total_pivots = 0

  // Phase I + Phase II root LP
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
    let (p1_st, p1_it) = run_simplex_phase(
      tab, basis, col_active, m, num_cols, config,
    )
    total_pivots = total_pivots + p1_it
    if p1_st != Optimal {
      return (Solution::empty(p1_st, total_pivots, 1), empty_stats)
    }
    if -tab[m][rhs_col] > config.tol {
      return (Solution::empty(Infeasible, total_pivots, 1), empty_stats)
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
  let (p2_st, p2_it) = run_simplex_phase(
    tab, basis, col_active, m, num_cols, config, costs=std_obj,
  )
  total_pivots = total_pivots + p2_it
  if p2_st == Unbounded {
    return (
      Solution::empty(
        resolve_unbounded_milp_status(model, config),
        total_pivots,
        1,
      ),
      empty_stats,
    )
  }
  if p2_st != Optimal {
    return (Solution::empty(p2_st, total_pivots, 1), empty_stats)
  }

  let root_lp_obj = (-tab[m][rhs_col] + obj_offset) * sense_sign

  // Root GMI Cutting Planes (tightening root tableau BEFORE branching!)
  let mut cuts_added = 0
  while cuts_added < max_root_cuts {
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
      break
    }

    let rhs_val = tab[best_r][rhs_col]
    let f0 = rhs_val - d_floor(rhs_val)
    if f0 <= 1.0e-5 || 1.0 - f0 <= 1.0e-5 {
      break
    }
    let scale = f0 / (1.0 - f0)
    let new_col = num_cols
    let mut f0_eff = f0 - 1.0e-11
    let mut cut_has_coeff = false
    let mut cut_safe = true
    let mut max_g = 0.0
    let mut min_g = 1.0e308
    let mut cut_slack_ub = -f0_eff
    let mut cut_slack_finite = true
    let cut_row = Array::make(num_cols + 2, 0.0)
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

    let (dst, dit) = run_dual_simplex(tab, basis, col_active, m, num_cols, config)
    total_pivots = total_pivots + dit
    if dst == Infeasible {
      return (
        Solution::empty(Infeasible, total_pivots, 1),
        {
          root_lp_obj,
          post_cut_lp_obj: root_lp_obj,
          cuts_added,
          heuristic_hits: 0,
          nodes_explored: 1,
          total_pivots,
        },
      )
    }
    if dst == MaxIterations {
      return (
        Solution::empty(MaxIterations, total_pivots, 1),
        {
          root_lp_obj,
          post_cut_lp_obj: (-tab[m][rhs_col] + obj_offset) * sense_sign,
          cuts_added,
          heuristic_hits: 0,
          nodes_explored: 1,
          total_pivots,
        },
      )
    }
    if dst != Optimal {
      break
    }
  }

  let post_cut_lp_obj = (-tab[m][rhs_col] + obj_offset) * sense_sign
  let root_min_bound = -tab[m][rhs_col] + obj_offset

  // Hybrid DFS-Plunge + Best-Bound Dual Warm-Start Branch-and-Cut Tree
  let stack : Array[BCNode] = [
    { tab, basis, col_active, m, num_cols, bound: root_min_bound, depth: 0, },
  ]
  let mut best_min_obj = 1.0e308
  let mut best_pruned_bound = 1.0e308
  let mut best_sol : Solution? = None
  let mut nodes_explored = 0
  let mut heuristic_hits = 0
  let mut hit_max_iters = false
  let mut has_unproved = false
  let mut retry_without_cuts = false

  while stack.length() > 0 && nodes_explored < config.max_nodes {
    let node = pop_hybrid_bc_node(stack, best_sol is Some(_), nodes_explored)
    if best_sol is Some(_) &&
      can_prune_by_bound(node.bound, best_min_obj, config.mip_gap_tol) {
      if node.bound < best_pruned_bound {
        best_pruned_bound = node.bound
      }
      continue
    }

    let (dst, dit) = run_dual_simplex(
      node.tab,
      node.basis,
      node.col_active,
      node.m,
      node.num_cols,
      config,
    )
    nodes_explored = nodes_explored + 1
    total_pivots = total_pivots + dit
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

    // Objective-weighted fractional branching variable selection
    let mut branch_row = -1
    let mut best_score = -1.0
    for i = 0; i < node.m; i = i + 1 {
      let bc = node.basis[i]
      if bc >= 0 && bc < n_std && is_int_col[bc] {
        let val = node.tab[i][node.num_cols]
        let dist = (val - d_round(val)).abs()
        if dist > config.int_tol {
          let score = dist * (1.0 + model.vars[bc].obj.abs())
          if score > best_score {
            best_score = score
            branch_row = i
          }
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
              iterations: total_pivots,
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
        // Rounding failed or degraded objective vs LP bound: branch on any sub-tol fractional var
        // strictly above dual-simplex feasibility tol so child bound cuts force strict progress!
        let mut sub_score = -1.0
        for i = 0; i < node.m; i = i + 1 {
          let bc = node.basis[i]
          if bc >= 0 && bc < n_std && is_int_col[bc] {
            let val = node.tab[i][node.num_cols]
            let dist = (val - d_round(val)).abs()
            if dist > config.tol {
              let score = dist * (1.0 + model.vars[bc].obj.abs())
              if score > sub_score {
                sub_score = score
                branch_row = i
              }
            }
          }
        }
        if branch_row == -1 {
          // A weakened cut may leave a gap too small for another warm pivot.
          // Re-prove the original model with a fresh tree and the remaining budget.
          if cuts_added > 0 && nodes_explored < config.max_nodes {
            retry_without_cuts = true
            break
          }
          has_unproved = true
        }
      }
    }

    if branch_row != -1 {
      // Try fast Primal Rounding Heuristic before branching
      match
        try_rounding_heuristic(
          model,
          lbs,
          node.tab,
          node.basis,
          node.m,
          n_std,
          node.num_cols,
          is_int_col,
          sense_sign,
          best_min_obj,
          config,
        ) {
        Some(hsol) => {
          best_min_obj = hsol.objective * sense_sign
          best_sol = Some(hsol)
          heuristic_hits = heuristic_hits + 1
          if can_prune_by_bound(min_obj, best_min_obj, config.mip_gap_tol) {
            if min_obj < best_pruned_bound {
              best_pruned_bound = min_obj
            }
            continue
          }
        }
        None => ()
      }

      let val = node.tab[branch_row][node.num_cols]
      let fl = d_floor(val)
      let cl = d_ceil(val)
      let c_down = make_bc_child(node, branch_row, true, fl, min_obj)
      let c_up = make_bc_child(node, branch_row, false, cl, min_obj)
      if val - fl <= cl - val {
        stack.push(c_up)
        stack.push(c_down)
      } else {
        stack.push(c_down)
        stack.push(c_up)
      }
    }
  }

  let stats : BranchAndCutStats = {
    root_lp_obj,
    post_cut_lp_obj,
    cuts_added,
    heuristic_hits,
    nodes_explored,
    total_pivots,
  }

  if retry_without_cuts {
    let remaining = { ..config, max_nodes: config.max_nodes - nodes_explored, }
    let cold = solve_milp_core(model, remaining, true, incumbent?=best_sol)
    let nodes = nodes_explored + cold.nodes_explored
    let pivots = total_pivots + cold.iterations
    return (
      { ..cold, nodes_explored: nodes, iterations: pivots, },
      { ..stats, nodes_explored: nodes, total_pivots: pivots, },
    )
  }
  match best_sol {
    Some(sol) =>
      if stack.length() == 0 && !has_unproved {
        let eff_bound = if best_pruned_bound < best_min_obj {
          best_pruned_bound
        } else {
          best_min_obj
        }
        let final_gap = compute_mip_gap(best_min_obj, eff_bound)
        (
          {
            ..sol,
            status: Optimal,
            iterations: total_pivots,
            nodes_explored,
            mip_gap: final_gap,
          },
          stats,
        )
      } else {
        let mut best_rem = best_pruned_bound
        for i = 0; i < stack.length(); i = i + 1 {
          if stack[i].bound < best_rem {
            best_rem = stack[i].bound
          }
        }
        if best_rem > best_min_obj {
          best_rem = best_min_obj
        }
        let gap = if has_unproved && stack.length() == 0 {
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
        (
          { ..sol, status: st, iterations: total_pivots, nodes_explored, mip_gap: gap, },
          stats,
        )
      }
    None =>
      (
        Solution::empty(
          if hit_max_iters {
            MaxIterations
          } else if nodes_explored >= config.max_nodes {
            NodeLimit
          } else if has_unproved {
            NumericalError
          } else {
            Infeasible
          },
          total_pivots,
          nodes_explored,
        ),
        stats,
      )
  }
}
