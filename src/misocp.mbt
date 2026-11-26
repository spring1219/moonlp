///|
/// Mixed-Integer Second-Order Cone Programming (MISOCP / MIQCQP) Solver
/// via Best-Bound Branch-and-Bound over Nesterov-Todd SOCP Relaxations (v2.0.0).
///
/// Combines:
///   1. Continuous SOCP / QCQP relaxations at every Branch-and-Bound node
///   2. Verified integer candidate acceptance with fixed-integer continuous SOCP repair
///   3. Periodic primal rounding + fixed-integer SOCP heuristic to capture incumbents early
///   4. Strict status & MIP gap propagation across all termination modes

///|
fn verify_or_repair_misocp_candidate(
  model : Model,
  raw_primal : Array[Double],
  config : SolverConfig,
) -> Solution? {
  let n = model.vars.length()
  let snapped : Array[Double] = []
  let mut has_continuous = false
  for j = 0; j < n; j = j + 1 {
    if is_integer_var(model.vars[j].kind) {
      snapped.push(d_round(raw_primal[j]))
    } else {
      has_continuous = true
      snapped.push(raw_primal[j])
    }
  }
  let direct : Solution = {
    status: Optimal,
    objective: model.eval_objective(snapped),
    primal: snapped,
    dual: [],
    reduced_costs: [],
    iterations: 0,
    nodes_explored: 0,
    mip_gap: 0.0,
  }
  let vtol = if config.tol * 100.0 > 1.0e-5 {
    config.tol * 100.0
  } else {
    1.0e-5
  }
  if model.verify_solution(direct, tol=vtol) {
    return Some(direct)
  }
  if !has_continuous {
    return None
  }
  let fixed_lbs : Array[Double] = []
  let fixed_ubs : Array[Double] = []
  for j = 0; j < n; j = j + 1 {
    if is_integer_var(model.vars[j].kind) {
      let sv = snapped[j]
      if sv < model.vars[j].lb - config.tol || sv > model.vars[j].ub + config.tol {
        return None
      }
      fixed_lbs.push(sv)
      fixed_ubs.push(sv)
    } else {
      fixed_lbs.push(model.vars[j].lb)
      fixed_ubs.push(model.vars[j].ub)
    }
  }
  let fixed_socp = model.clone_with_bounds(fixed_lbs, fixed_ubs)
  let cont_sol = solve_socp(fixed_socp, config~)
  if cont_sol.status == Optimal && model.verify_solution(cont_sol, tol=vtol) {
    Some(cont_sol)
  } else {
    None
  }
}

///|
/// Solve a Mixed-Integer Second-Order Cone Program (MISOCP) or MIQCQP via
/// Best-First Branch-and-Bound with Nesterov-Todd SOCP node relaxations.
pub fn solve_misocp(
  model : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  if !model.has_soc_constraints() {
    if model.has_quad_obj() {
      return solve_miqp(model, config~)
    } else {
      return solve_milp(model, config~)
    }
  }
  if !model.has_integer_vars() {
    return solve_socp(model, config~)
  }

  let sense_sign = match model.sense {
    Minimize => 1.0
    Maximize => -1.0
  }
  let n = model.vars.length()
  let root_lbs : Array[Double] = []
  let root_ubs : Array[Double] = []
  for j = 0; j < n; j = j + 1 {
    let v = model.vars[j]
    let lb = match v.kind {
      Binary => if v.lb > 0.0 { d_ceil(v.lb - config.tol) } else { 0.0 }
      Integer =>
        if v.lb > -inf_threshold {
          d_ceil(v.lb - config.tol)
        } else {
          v.lb
        }
      Continuous => v.lb
    }
    let ub = match v.kind {
      Binary => if v.ub < 1.0 { d_floor(v.ub + config.tol) } else { 1.0 }
      Integer =>
        if v.ub < inf_threshold {
          d_floor(v.ub + config.tol)
        } else {
          v.ub
        }
      Continuous => v.ub
    }
    if lb > ub + config.tol {
      return Solution::empty(Infeasible, 0, 1)
    }
    root_lbs.push(lb)
    root_ubs.push(ub)
  }

  let queue : Array[BBNode] = [
    { lbs: root_lbs, ubs: root_ubs, bound: -1.0e308, },
  ]
  let mut best_min_obj = 1.0e308
  let mut best_pruned_bound = 1.0e308
  let mut best_sol : Solution? = None
  let mut nodes_explored = 0
  let mut total_iters = 0
  let mut hit_unbounded = false
  let mut hit_max_iters = false
  let mut hit_nonconvex = false
  let mut has_unproved = false

  while queue.length() > 0 && nodes_explored < config.max_nodes {
    let node = pop_best_node(queue)

    if best_sol is Some(_) &&
      can_prune_by_bound(node.bound, best_min_obj, config.mip_gap_tol) {
      if node.bound < best_pruned_bound {
        best_pruned_bound = node.bound
      }
      continue
    }

    let sub_model = model.clone_with_bounds(node.lbs, node.ubs)
    let socp_sol = solve_socp(sub_model, config~)
    nodes_explored = nodes_explored + 1
    total_iters = total_iters + socp_sol.iterations

    match socp_sol.status {
      Infeasible => continue
      Unbounded => {
        hit_unbounded = true
        has_unproved = true
        continue
      }
      NonConvex => {
        hit_nonconvex = true
        has_unproved = true
        continue
      }
      MaxIterations => {
        hit_max_iters = true
        has_unproved = true
        continue
      }
      NodeLimit | InfeasibleOrUnbounded | NumericalError => {
        has_unproved = true
        continue
      }
      Optimal => {
        let min_obj = socp_sol.objective * sense_sign
        if best_sol is Some(_) &&
          can_prune_by_bound(min_obj, best_min_obj, config.mip_gap_tol) {
          if min_obj < best_pruned_bound {
            best_pruned_bound = min_obj
          }
          continue
        }

        let mut branch_idx = pick_branch_var(
          model,
          socp_sol.primal,
          config.int_tol,
        )
        if branch_idx == -1 {
          let mut fathomed = false
          match
            verify_or_repair_misocp_candidate(model, socp_sol.primal, config) {
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
            branch_idx = pick_branch_var(model, socp_sol.primal, 1.0e-12)
            if branch_idx == -1 {
              has_unproved = true
            }
          }
        } else if nodes_explored <= 12 || nodes_explored % 4 == 0 {
          match
            verify_or_repair_misocp_candidate(model, socp_sol.primal, config) {
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
              if can_prune_by_bound(min_obj, best_min_obj, config.mip_gap_tol) {
                if min_obj < best_pruned_bound {
                  best_pruned_bound = min_obj
                }
                continue
              }
            }
            None => ()
          }
        }

        if branch_idx != -1 {
          let val = socp_sol.primal[branch_idx]
          let fl = d_floor(val)
          let cl = d_ceil(val)
          if fl >= node.lbs[branch_idx] - config.tol &&
            fl < node.ubs[branch_idx] + config.tol {
            let left_lbs = node.lbs.copy()
            let left_ubs = node.ubs.copy()
            left_ubs[branch_idx] = fl
            queue.push({ lbs: left_lbs, ubs: left_ubs, bound: min_obj, })
          }
          if cl <= node.ubs[branch_idx] + config.tol &&
            cl > node.lbs[branch_idx] - config.tol {
            let right_lbs = node.lbs.copy()
            let right_ubs = node.ubs.copy()
            right_lbs[branch_idx] = cl
            queue.push({ lbs: right_lbs, ubs: right_ubs, bound: min_obj, })
          }
        }
      }
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
        } else if hit_nonconvex {
          NonConvex
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
      if hit_nonconvex {
        Solution::empty(NonConvex, total_iters, nodes_explored)
      } else if hit_max_iters {
        Solution::empty(MaxIterations, total_iters, nodes_explored)
      } else if hit_unbounded {
        Solution::empty(InfeasibleOrUnbounded, total_iters, nodes_explored)
      } else if nodes_explored >= config.max_nodes {
        Solution::empty(NodeLimit, total_iters, nodes_explored)
      } else if has_unproved {
        Solution::empty(NumericalError, total_iters, nodes_explored)
      } else {
        Solution::empty(Infeasible, total_iters, nodes_explored)
      }
  }
}
