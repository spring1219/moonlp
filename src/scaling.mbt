///|
/// Ruiz Iterative L-infinity Matrix Equilibration Scaling (v0.5.0):
/// Dramatically reduces numerical dynamic range (condition ratio max|A_ij| / min|A_ij|)
/// on ill-conditioned industrial LPs (e.g. Netlib benchmarks spanning 1e-4 .. 1e6).

///|
#derive(Debug)
pub(all) struct ScalingReport {
  original_ratio : Double
  scaled_ratio : Double
  passes : Int32
  row_scales : Array[Double]
  col_scales : Array[Double]
}

///|
fn d_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  // Exact power-of-4 range reduction to [1.0, 4.0] so Newton converges to 1e-16
  // across the entire IEEE-754 range (1e-300 .. 1e300).
  let mut y = x
  let mut factor = 1.0
  while y < 1.0 && y > 0.0 {
    y = y * 4.0
    factor = factor * 0.5
  }
  while y > 4.0 {
    y = y * 0.25
    factor = factor * 2.0
  }
  let mut g = 0.5 * (y + 1.0)
  for _i = 0; _i < 16; _i = _i + 1 {
    let next = 0.5 * (g + y / g)
    if (next - g).abs() <= 1.0e-15 * g {
      return next * factor
    }
    g = next
  }
  g * factor
}

///|
/// Compute dynamic range ratio (max|a_ij| / min|a_ij|) of nonzero constraint coefficients.
pub fn Model::matrix_scale_ratio(self : Model) -> Double {
  let mut min_v = 1.0e308
  let mut max_v = 0.0
  for i = 0; i < self.constraints.length(); i = i + 1 {
    let terms = self.constraints[i].terms
    for k = 0; k < terms.length(); k = k + 1 {
      let a = terms[k].coeff.abs()
      if a > 1.0e-30 {
        if a < min_v {
          min_v = a
        }
        if a > max_v {
          max_v = a
        }
      }
    }
  }
  if max_v <= 0.0 || min_v >= 1.0e307 {
    1.0
  } else {
    max_v / min_v
  }
}

///|
/// Apply `passes` iterations of Ruiz L-infinity equilibration scaling to an LP Model.
/// Returns (scaled_model, scaling_report).
/// Exact variable/dual mapping:
///   x_orig[j]  = col_scales[j] * x_scaled[j]
///   y_orig[i]  = row_scales[i] * y_scaled[i]
///   rc_orig[j] = rc_scaled[j] / col_scales[j]
///   obj_orig   = obj_scaled (invariant)
fn Model::ruiz_scale_core(
  self : Model,
  passes : Int32,
  relax_int : Bool,
) -> (Model, ScalingReport) {
  let m = self.constraints.length()
  let n = self.vars.length()
  let orig_ratio = self.matrix_scale_ratio()
  let r_acc = Array::make(m, 1.0)
  let c_acc = Array::make(n, 1.0)

  // Only scale continuous LPs unless explicitly relaxing integrality for LP relaxation
  if (!relax_int && self.has_integer_vars()) || m == 0 || n == 0 {
    return (
      self,
      {
        original_ratio: orig_ratio,
        scaled_ratio: orig_ratio,
        passes: 0,
        row_scales: r_acc,
        col_scales: c_acc,
      },
    )
  }

  // Work on mutable copy of constraint term coefficients
  let work_terms : Array[Array[Term]] = []
  for i = 0; i < m; i = i + 1 {
    work_terms.push(self.constraints[i].terms.copy())
  }

  for _p = 0; _p < passes; _p = _p + 1 {
    let row_max = Array::make(m, 0.0)
    let col_max = Array::make(n, 0.0)
    for i = 0; i < m; i = i + 1 {
      let terms = work_terms[i]
      for k = 0; k < terms.length(); k = k + 1 {
        let t = terms[k]
        let a = t.coeff.abs()
        if a > row_max[i] {
          row_max[i] = a
        }
        if a > col_max[t.var_id.index] {
          col_max[t.var_id.index] = a
        }
      }
    }

    let dr = Array::make(m, 1.0)
    let dc = Array::make(n, 1.0)
    for i = 0; i < m; i = i + 1 {
      if row_max[i] > 1.0e-30 {
        dr[i] = 1.0 / d_sqrt(row_max[i])
        r_acc[i] = r_acc[i] * dr[i]
      }
    }
    for j = 0; j < n; j = j + 1 {
      if col_max[j] > 1.0e-30 {
        dc[j] = 1.0 / d_sqrt(col_max[j])
        c_acc[j] = c_acc[j] * dc[j]
      }
    }

    for i = 0; i < m; i = i + 1 {
      let terms = work_terms[i]
      let dri = dr[i]
      for k = 0; k < terms.length(); k = k + 1 {
        let t = terms[k]
        terms[k] = {
          var_id: t.var_id,
          coeff: t.coeff * dri * dc[t.var_id.index],
        }
      }
    }
  }

  // Construct scaled model
  let scaled_vars : Array[VarSpec] = []
  for j = 0; j < n; j = j + 1 {
    let v = self.vars[j]
    let cj = c_acc[j]
    let new_lb = if v.lb > -inf_threshold { v.lb / cj } else { v.lb }
    let new_ub = if v.ub < inf_threshold { v.ub / cj } else { v.ub }
    scaled_vars.push({
      ..v,
      lb: new_lb,
      ub: new_ub,
      obj: v.obj * cj,
    })
  }
  let scaled_constraints : Array[Constraint] = []
  for i = 0; i < m; i = i + 1 {
    let c = self.constraints[i]
    scaled_constraints.push({
      ..c,
      terms: work_terms[i],
      rhs: c.rhs * r_acc[i],
    })
  }

  let scaled_model : Model = {
    sense: self.sense,
    vars: scaled_vars,
    constraints: scaled_constraints,
    obj_constant: self.obj_constant,
  }
  let scaled_ratio = scaled_model.matrix_scale_ratio()
  (
    scaled_model,
    {
      original_ratio: orig_ratio,
      scaled_ratio,
      passes,
      row_scales: r_acc,
      col_scales: c_acc,
    },
  )
}

///|
/// Apply `passes` iterations of Ruiz L-infinity equilibration scaling to an LP Model.
/// Returns (scaled_model, scaling_report).
pub fn Model::ruiz_scale(
  self : Model,
  passes~ : Int32 = 6,
) -> (Model, ScalingReport) {
  self.ruiz_scale_core(passes, false)
}

///|
/// Solve the LP relaxation (ignoring integrality) with Ruiz L-infinity scaling
/// and exact unscaling back to original variable/dual space.
fn solve_lp_relaxation_scaled(
  model : Model,
  passes~ : Int32 = 6,
  config~ : SolverConfig = SolverConfig::default(),
) -> (Solution, ScalingReport) {
  let (scaled, rep) = model.ruiz_scale_core(passes, true)
  let ssol = solve_lp_unscaled(scaled, config)
  if ssol.status != Optimal {
    return (ssol, rep)
  }
  let primal : Array[Double] = []
  for j = 0; j < ssol.primal.length(); j = j + 1 {
    primal.push(ssol.primal[j] * rep.col_scales[j])
  }
  let dual : Array[Double] = []
  for i = 0; i < ssol.dual.length(); i = i + 1 {
    dual.push(ssol.dual[i] * rep.row_scales[i])
  }
  let rc : Array[Double] = []
  for j = 0; j < ssol.reduced_costs.length(); j = j + 1 {
    rc.push(ssol.reduced_costs[j] / rep.col_scales[j])
  }
  ({ ..ssol, primal, dual, reduced_costs: rc, }, rep)
}

///|
/// Solve an LP with Ruiz L-infinity equilibration scaling and exact unscaling of
/// primal values, dual shadow prices, and reduced costs.
pub fn Model::solve_scaled(
  self : Model,
  passes~ : Int32 = 6,
  config~ : SolverConfig = SolverConfig::default(),
) -> (Solution, ScalingReport) {
  let (scaled, rep) = self.ruiz_scale(passes~)
  let ssol = if scaled.has_integer_vars() {
    solve_milp(scaled, config~)
  } else {
    solve_lp_unscaled(scaled, config)
  }
  if ssol.status != Optimal {
    return (ssol, rep)
  }

  let primal : Array[Double] = []
  for j = 0; j < ssol.primal.length(); j = j + 1 {
    primal.push(ssol.primal[j] * rep.col_scales[j])
  }
  let dual : Array[Double] = []
  for i = 0; i < ssol.dual.length(); i = i + 1 {
    dual.push(ssol.dual[i] * rep.row_scales[i])
  }
  let rc : Array[Double] = []
  for j = 0; j < ssol.reduced_costs.length(); j = j + 1 {
    rc.push(ssol.reduced_costs[j] / rep.col_scales[j])
  }

  (
    { ..ssol, primal, dual, reduced_costs: rc, },
    rep,
  )
}
