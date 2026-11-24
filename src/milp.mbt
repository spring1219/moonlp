///|
/// Mixed-Integer Linear Programming (MILP) Solver via Best-Bound Branch-and-Bound (v0.6.0).
/// Hardened status propagation, zero-objective unboundedness disambiguation (Meyer's theorem),
/// and verified integer-candidate acceptance with sub-tolerance fallback branching.

///|
priv struct BBNode {
  lbs : Array[Double]
  ubs : Array[Double]
  bound : Double // LP relaxation bound (in minimization space)
}

///|
fn is_integer_var(kind : VarKind) -> Bool {
  match kind {
    Continuous => false
    Integer | Binary => true
  }
}

///|
/// Find the integer/binary variable furthest from an integer beyond `frac_tol`.
/// Returns -1 if all integer/binary variables are within `frac_tol` of an integer.
fn pick_branch_var(
  model : Model,
  primal : Array[Double],
  frac_tol : Double,
) -> Int32 {
  let mut best_idx = -1
  let mut best_frac = 0.0
  for j = 0; j < model.vars.length(); j = j + 1 {
    if !is_integer_var(model.vars[j].kind) {
      continue
    }
    let val = primal[j]
    let rounded = d_round(val)
    let dist = (val - rounded).abs()
    if dist > frac_tol && dist > best_frac {
      best_frac = dist
      best_idx = j
    }
  }
  best_idx
}

///|
/// Pop the node with the best (lowest) minimization LP bound.
fn pop_best_node(queue : Array[BBNode]) -> BBNode {
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
/// Build a zero-objective copy of `model` to test pure integer feasibility.
fn make_zero_obj_model(model : Model) -> Model {
  let new_vars : Array[VarSpec] = []
  for j = 0; j < model.vars.length(); j = j + 1 {
    let v = model.vars[j]
    new_vars.push({ ..v, obj: 0.0, })
  }
  {
    sense: Minimize,
    vars: new_vars,
    constraints: model.constraints,
    obj_constant: 0.0,
  }
}

///|
/// Disambiguate an unbounded MILP LP relaxation via Meyer's theorem:
/// Solve zero-objective integer feasibility model:
///   - Feasible integer point exists => MILP is Unbounded
///   - Proven no integer feasible point => MILP is Infeasible
///   - Budget exhausted => InfeasibleOrUnbounded
fn resolve_unbounded_milp_status(
  model : Model,
  config : SolverConfig,
) -> Status {
  let feas_model = make_zero_obj_model(model)
  let feas_sol = solve_milp_core(feas_model, config, false)
  match feas_sol.status {
    Optimal => Unbounded
    Infeasible => Infeasible
    _ => InfeasibleOrUnbounded
  }
}

///|
/// Try fixing integer variables to rounded values and re-solving continuous variables
/// to verify/repair an integral candidate solution against the original model.
fn verify_or_repair_int_candidate(
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
  let mut exact_obj = model.obj_constant
  for j = 0; j < n; j = j + 1 {
    exact_obj = exact_obj + model.vars[j].obj * snapped[j]
  }
  let direct : Solution = {
    status: Optimal,
    objective: exact_obj,
    primal: snapped,
    dual: [],
    reduced_costs: [],
    iterations: 0,
    nodes_explored: 0,
    mip_gap: 0.0,
  }
  if model.verify_solution(direct, tol=1.0e-5) {
    return Some(direct)
  }
  if !has_continuous {
    return None
  }
  // Re-solve continuous variables with integer variables fixed to snapped values
  let fixed_lbs : Array[Double] = []
  let fixed_ubs : Array[Double] = []
  for j = 0; j < n; j = j + 1 {
    if is_integer_var(model.vars[j].kind) {
      fixed_lbs.push(snapped[j])
      fixed_ubs.push(snapped[j])
    } else {
      fixed_lbs.push(model.vars[j].lb)
      fixed_ubs.push(model.vars[j].ub)
    }
  }
  let fixed_lp = model.clone_with_bounds(fixed_lbs, fixed_ubs)
  let cont_sol = solve_lp(fixed_lp, config~)
  if cont_sol.status == Optimal && model.verify_solution(cont_sol, tol=1.0e-5) {
    Some(cont_sol)
  } else {
    None
  }
}

///|
fn solve_milp_core(
  model : Model,
  config : SolverConfig,
  check_unbounded : Bool,
) -> Solution {
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
  let mut best_sol : Solution? = None
  let mut nodes_explored = 0
  let mut total_iters = 0
  let mut hit_unbounded = false
  let mut hit_max_iters = false
  let mut has_unproved = false

  while queue.length() > 0 && nodes_explored < config.max_nodes {
    let node = pop_best_node(queue)

    if best_sol is Some(_) && node.bound >= best_min_obj - config.mip_gap_tol {
      continue
    }

    let sub_model = model.clone_with_bounds(node.lbs, node.ubs)
    let lp_sol = solve_lp(sub_model, config~)
    nodes_explored = nodes_explored + 1
    total_iters = total_iters + lp_sol.iterations

    match lp_sol.status {
      Infeasible => continue
      Unbounded => {
        hit_unbounded = true
        has_unproved = true
        continue
      }
      MaxIterations => {
        hit_max_iters = true
        has_unproved = true
        continue
      }
      NodeLimit | InfeasibleOrUnbounded => {
        has_unproved = true
        continue
      }
      Optimal => {
        let min_obj = lp_sol.objective * sense_sign
        if best_sol is Some(_) && min_obj >= best_min_obj - config.mip_gap_tol {
          continue
        }

        let mut branch_idx = pick_branch_var(
          model,
          lp_sol.primal,
          config.int_tol,
        )
        if branch_idx == -1 {
          match verify_or_repair_int_candidate(model, lp_sol.primal, config) {
            Some(cand) =>
              if cand.objective * sense_sign < best_min_obj {
                best_min_obj = cand.objective * sense_sign
                best_sol = Some({
                  ..cand,
                  iterations: total_iters,
                  nodes_explored,
                  mip_gap: 0.0,
                })
              }
            None => {
              // Rounding failed original constraint verification (e.g. 1e8 * 0.9999995):
              // branch on any sub-tolerance fractional integer variable instead of dropping node!
              branch_idx = pick_branch_var(model, lp_sol.primal, 1.0e-12)
              if branch_idx == -1 {
                has_unproved = true
              }
            }
          }
        }

        if branch_idx != -1 {
          let val = lp_sol.primal[branch_idx]
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
        { ..sol, status: Optimal, iterations: total_iters, nodes_explored, mip_gap: 0.0, }
      } else {
        let mut best_rem = best_min_obj
        if queue.length() > 0 {
          best_rem = queue[0].bound
          for i = 1; i < queue.length(); i = i + 1 {
            if queue[i].bound < best_rem {
              best_rem = queue[i].bound
            }
          }
        }
        let denom = if best_min_obj.abs() > 1.0e-9 {
          best_min_obj.abs()
        } else {
          1.0
        }
        let gap = if has_unproved && queue.length() == 0 {
          1.0
        } else {
          (best_min_obj - best_rem).abs() / denom
        }
        let st = if !has_unproved && gap <= config.mip_gap_tol {
          Optimal
        } else if hit_max_iters && nodes_explored < config.max_nodes {
          MaxIterations
        } else {
          NodeLimit
        }
        { ..sol, status: st, iterations: total_iters, nodes_explored, mip_gap: gap, }
      }
    None =>
      if hit_max_iters {
        Solution::empty(MaxIterations, total_iters, nodes_explored)
      } else if hit_unbounded {
        let st = if check_unbounded {
          resolve_unbounded_milp_status(model, config)
        } else {
          InfeasibleOrUnbounded
        }
        Solution::empty(st, total_iters, nodes_explored)
      } else if nodes_explored >= config.max_nodes || has_unproved {
        Solution::empty(NodeLimit, total_iters, nodes_explored)
      } else {
        Solution::empty(Infeasible, total_iters, nodes_explored)
      }
  }
}

///|
/// Solve a Mixed-Integer Linear Program using Best-First Branch-and-Bound.
pub fn solve_milp(
  model : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  solve_milp_core(model, config, true)
}
