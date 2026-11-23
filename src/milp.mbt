///|
/// Mixed-Integer Linear Programming (MILP) Solver via Best-Bound Branch-and-Bound.
/// Supports Binary, General Integer, and Continuous mixed variables.

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
/// Find the integer/binary variable with value furthest from an integer (closest to 0.5).
/// Returns -1 if all integer/binary variables are integral within `int_tol`.
fn pick_branch_var(
  model : Model,
  primal : Array[Double],
  int_tol : Double,
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
    if dist > int_tol && dist > best_frac {
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
/// Solve a Mixed-Integer Linear Program using Best-First Branch-and-Bound.
pub fn solve_milp(
  model : Model,
  config~ : SolverConfig = SolverConfig::default(),
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
    match v.kind {
      Binary => {
        root_lbs.push(if v.lb > 0.0 { d_ceil(v.lb) } else { 0.0 })
        root_ubs.push(if v.ub < 1.0 { d_floor(v.ub) } else { 1.0 })
      }
      Integer => {
        root_lbs.push(if v.lb > -1.0e200 { d_ceil(v.lb) } else { v.lb })
        root_ubs.push(if v.ub < 1.0e200 { d_floor(v.ub) } else { v.ub })
      }
      Continuous => {
        root_lbs.push(v.lb)
        root_ubs.push(v.ub)
      }
    }
  }

  let queue : Array[BBNode] = [
    { lbs: root_lbs, ubs: root_ubs, bound: -1.0e308, },
  ]
  let mut best_min_obj = 1.0e308
  let mut best_sol : Solution? = None
  let mut nodes_explored = 0
  let mut total_iters = 0
  let mut hit_unbounded = false

  while queue.length() > 0 && nodes_explored < config.max_nodes {
    let node = pop_best_node(queue)

    // Prune by bound before solving if parent bound already worse than incumbent
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
        continue
      }
      MaxIterations | NodeLimit => continue
      Optimal => {
        let min_obj = lp_sol.objective * sense_sign
        // Prune by bound
        if best_sol is Some(_) && min_obj >= best_min_obj - config.mip_gap_tol {
          continue
        }

        let branch_idx = pick_branch_var(model, lp_sol.primal, config.int_tol)
        if branch_idx == -1 {
          // Integral feasible solution found — snap integer vars to exact integers
          let snapped : Array[Double] = []
          for j = 0; j < n; j = j + 1 {
            if is_integer_var(model.vars[j].kind) {
              snapped.push(d_round(lp_sol.primal[j]))
            } else {
              snapped.push(lp_sol.primal[j])
            }
          }
          // Recompute exact objective on snapped vector
          let mut exact_obj = model.obj_constant
          for j = 0; j < n; j = j + 1 {
            exact_obj = exact_obj + model.vars[j].obj * snapped[j]
          }
          best_min_obj = exact_obj * sense_sign
          best_sol = Some({
            status: Optimal,
            objective: exact_obj,
            primal: snapped,
            dual: lp_sol.dual,
            reduced_costs: lp_sol.reduced_costs,
            iterations: total_iters,
            nodes_explored,
            mip_gap: 0.0,
          })
        } else {
          // Branch on x[branch_idx] <= floor(val) and x[branch_idx] >= ceil(val)
          let val = lp_sol.primal[branch_idx]
          let fl = d_floor(val)
          let cl = d_ceil(val)

          if fl >= node.lbs[branch_idx] - config.tol {
            let left_lbs = node.lbs.copy()
            let left_ubs = node.ubs.copy()
            left_ubs[branch_idx] = fl
            queue.push({ lbs: left_lbs, ubs: left_ubs, bound: min_obj, })
          }
          if cl <= node.ubs[branch_idx] + config.tol {
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
      if queue.length() == 0 {
        { ..sol, iterations: total_iters, nodes_explored, mip_gap: 0.0, }
      } else {
        // Compute remaining best lower bound gap if node limit reached
        let mut best_rem = queue[0].bound
        for i = 1; i < queue.length(); i = i + 1 {
          if queue[i].bound < best_rem {
            best_rem = queue[i].bound
          }
        }
        let denom = if best_min_obj.abs() > 1.0e-9 {
          best_min_obj.abs()
        } else {
          1.0
        }
        let gap = (best_min_obj - best_rem).abs() / denom
        {
          ..sol,
          status: if gap <= config.mip_gap_tol {
            Optimal
          } else {
            NodeLimit
          },
          iterations: total_iters,
          nodes_explored,
          mip_gap: gap,
        }
      }
    None =>
      if hit_unbounded {
        Solution::empty(Unbounded, total_iters, nodes_explored)
      } else if nodes_explored >= config.max_nodes {
        Solution::empty(NodeLimit, total_iters, nodes_explored)
      } else {
        Solution::empty(Infeasible, total_iters, nodes_explored)
      }
  }
}
