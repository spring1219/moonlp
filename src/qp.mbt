///|
/// Convex Quadratic Programming (QP) Solver via Scale-Normalized Mehrotra
/// Predictor-Corrector Primal-Dual Interior-Point Method (v1.0.0).
///
/// Solves general convex quadratic programs:
///   min/max  0.5 * x^T Q x + c^T x + c_0
///   s.t.     A_eq x = b_eq, A_le x <= b_le, A_ge x >= b_ge, l <= x <= u
///
/// Key engineering properties:
/// - Supports arbitrary positive semi-definite Hessians Q >= 0 (including rank-deficient
///   Q from MIQP relaxations with binary indicator variables having zero Hessian blocks).
/// - Validates convexity (PSD check via shifted Cholesky) and rejects non-convex QPs with
///   explicit `NonConvex` status.
/// - Scale-normalized augmented KKT formulation avoids dividing by tiny s_i or z_i near
///   convergence, preserving full IEEE-754 precision down to mu ~ 1e-13.
/// - Exact primal feasibility pre-check & unbounded recession-ray certification.

///|
/// Check whether symmetric n x n matrix H is positive semi-definite (H + eps*I > 0)
/// via in-place Cholesky factorization on a copy.
fn is_symmetric_psd(h : Array[Array[Double]], n : Int32) -> Bool {
  if n <= 0 {
    return true
  }
  let mut h_max = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let a = h[i][j].abs()
      if !is_finite_double(a) {
        return false
      }
      if a > h_max {
        h_max = a
      }
    }
  }
  if h_max == 0.0 {
    return true
  }
  let eps = 1.0e-10 * (1.0 + h_max)
  let l : Array[Array[Double]] = []
  for _i = 0; _i < n; _i = _i + 1 {
    l.push(Array::make(n, 0.0))
  }
  for i = 0; i < n; i = i + 1 {
    for j = 0; j <= i; j = j + 1 {
      let mut sum = h[i][j]
      if i == j {
        sum = sum + eps
      }
      for k = 0; k < j; k = k + 1 {
        sum = sum - l[i][k] * l[j][k]
      }
      if i == j {
        if sum <= 0.0 || !is_finite_double(sum) {
          return false
        }
        l[i][j] = d_sqrt(sum)
      } else {
        l[i][j] = sum / l[j][j]
      }
    }
  }
  true
}

///|
/// Dense LU factorization with partial pivoting: P * A = L * U.
/// Returns Some((lu, perm)) on success, None if numerically singular.
fn lu_factorize(
  mat : Array[Array[Double]],
  n : Int32,
) -> (Array[Array[Double]], Array[Int32])? {
  let lu : Array[Array[Double]] = []
  let perm : Array[Int32] = []
  for i = 0; i < n; i = i + 1 {
    lu.push(mat[i].copy())
    perm.push(i)
  }
  for k = 0; k < n; k = k + 1 {
    let mut piv_row = k
    let mut max_val = lu[k][k].abs()
    for i = k + 1; i < n; i = i + 1 {
      let v = lu[i][k].abs()
      if v > max_val {
        max_val = v
        piv_row = i
      }
    }
    if max_val < 1.0e-15 || !is_finite_double(max_val) {
      return None
    }
    if piv_row != k {
      let tmp_r = lu[k]
      lu[k] = lu[piv_row]
      lu[piv_row] = tmp_r
      let tmp_p = perm[k]
      perm[k] = perm[piv_row]
      perm[piv_row] = tmp_p
    }
    let piv = lu[k][k]
    let inv_piv = 1.0 / piv
    let row_k = lu[k]
    for i = k + 1; i < n; i = i + 1 {
      let row_i = lu[i]
      let f = row_i[k] * inv_piv
      row_i[k] = f
      if f != 0.0 {
        for j = k + 1; j < n; j = j + 1 {
          row_i[j] = row_i[j] - f * row_k[j]
        }
      }
    }
  }
  Some((lu, perm))
}

///|
/// Solve A * x = b using precomputed LU factorization with partial pivoting.
fn lu_solve(
  lu : Array[Array[Double]],
  perm : Array[Int32],
  b : Array[Double],
  n : Int32,
) -> Array[Double] {
  let y = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut sum = b[perm[i]]
    let row_i = lu[i]
    for j = 0; j < i; j = j + 1 {
      sum = sum - row_i[j] * y[j]
    }
    y[i] = sum
  }
  let x = Array::make(n, 0.0)
  let mut i = n - 1
  while i >= 0 {
    let mut sum = y[i]
    let row_i = lu[i]
    for j = i + 1; j < n; j = j + 1 {
      sum = sum - row_i[j] * x[j]
    }
    x[i] = sum / row_i[i]
    i = i - 1
  }
  x
}

///|
/// Solve A * x = b with 1 step of iterative refinement using original matrix A.
fn lu_solve_refined(
  orig : Array[Array[Double]],
  lu : Array[Array[Double]],
  perm : Array[Int32],
  b : Array[Double],
  n : Int32,
) -> Array[Double] {
  let x = lu_solve(lu, perm, b, n)
  let res = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let row_i = orig[i]
    let mut ax = 0.0
    for j = 0; j < n; j = j + 1 {
      ax = ax + row_i[j] * x[j]
    }
    res[i] = b[i] - ax
  }
  let dx = lu_solve(lu, perm, res, n)
  for i = 0; i < n; i = i + 1 {
    if is_finite_double(dx[i]) {
      x[i] = x[i] + dx[i]
    }
  }
  x
}

///|
fn max_positive_step(v : Array[Double], dv : Array[Double]) -> Double {
  let mut alpha = 1.0
  for i = 0; i < v.length(); i = i + 1 {
    if dv[i] < 0.0 {
      let cand = -v[i] / dv[i]
      if cand < alpha {
        alpha = cand
      }
    }
  }
  if alpha < 0.0 {
    0.0
  } else {
    alpha
  }
}

///|
fn vec_inf_norm(v : Array[Double]) -> Double {
  let mut m = 0.0
  for i = 0; i < v.length(); i = i + 1 {
    let a = v[i].abs()
    if a > m {
      m = a
    }
  }
  m
}

///|
/// Solve a Convex Quadratic Program (continuous relaxation if integer variables exist)
/// using a Scale-Normalized Mehrotra Predictor-Corrector Primal-Dual Interior-Point Method.
pub fn solve_qp(
  model : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  if !model.has_quad_obj() {
    return solve_lp(model, config~)
  }

  let n_orig = model.vars.length()
  let m_user = model.constraints.length()
  let sense_sign = match model.sense {
    Minimize => 1.0
    Maximize => -1.0
  }

  // 1. Validate finiteness & bounds, partition active vs fixed variables
  let orig_to_act = Array::make(n_orig, -1)
  let act_to_orig : Array[Int32] = []
  let fixed_val = Array::make(n_orig, 0.0)
  let is_fixed = Array::make(n_orig, false)

  for j = 0; j < n_orig; j = j + 1 {
    let v = model.vars[j]
    if !is_finite_double(v.obj) ||
      (v.lb > -inf_threshold && !is_finite_double(v.lb)) ||
      (v.ub < inf_threshold && !is_finite_double(v.ub)) {
      return Solution::empty(NumericalError, 0, 1)
    }
    if v.lb > v.ub + config.tol {
      return Solution::empty(Infeasible, 0, 1)
    }
    let has_lb = v.lb > -inf_threshold
    let has_ub = v.ub < inf_threshold
    if has_lb && has_ub && v.lb >= v.ub {
      is_fixed[j] = true
      fixed_val[j] = v.lb
    } else {
      orig_to_act[j] = act_to_orig.length()
      act_to_orig.push(j)
    }
  }

  let n = act_to_orig.length()

  // Assemble full symmetric Hessian Q_min = sense_sign * Q
  let q_full = model.build_symmetric_hessian()
  for i = 0; i < n_orig; i = i + 1 {
    for j = 0; j < n_orig; j = j + 1 {
      let v = q_full[i][j] * sense_sign
      if !is_finite_double(v) {
        return Solution::empty(NumericalError, 0, 1)
      }
      q_full[i][j] = v
    }
  }

  // All variables fixed: evaluate feasibility & exact objective directly
  if n == 0 {
    for i = 0; i < m_user; i = i + 1 {
      let c = model.constraints[i]
      let mut lhs = 0.0
      for k = 0; k < c.terms.length(); k = k + 1 {
        let t = c.terms[k]
        lhs = lhs + t.coeff * fixed_val[t.var_id.index]
      }
      let ok = match c.relation {
        Le => lhs <= c.rhs + config.tol
        Ge => lhs >= c.rhs - config.tol
        Eq => (lhs - c.rhs).abs() <= config.tol
      }
      if !ok {
        return Solution::empty(Infeasible, 0, 1)
      }
    }
    let primal = fixed_val.copy()
    let exact_obj = model.eval_objective(primal)
    let dual = Array::make(m_user, 0.0)
    let rc : Array[Double] = []
    for j = 0; j < n_orig; j = j + 1 {
      let mut grad = model.vars[j].obj
      for k = 0; k < n_orig; k = k + 1 {
        grad = grad + q_full[j][k] * sense_sign * primal[k]
      }
      rc.push(grad)
    }
    return {
      status: Optimal,
      objective: exact_obj,
      primal,
      dual,
      reduced_costs: rc,
      iterations: 0,
      nodes_explored: 1,
      mip_gap: 0.0,
    }
  }

  // Extract active Hessian H (n x n) and reduced linear gradient g (n)
  let h_mat : Array[Array[Double]] = []
  let g_vec = Array::make(n, 0.0)
  for a = 0; a < n; a = a + 1 {
    let orig_a = act_to_orig[a]
    let row = Array::make(n, 0.0)
    for b = 0; b < n; b = b + 1 {
      row[b] = q_full[orig_a][act_to_orig[b]]
    }
    h_mat.push(row)
    let mut ga = model.vars[orig_a].obj * sense_sign
    for j = 0; j < n_orig; j = j + 1 {
      if is_fixed[j] {
        ga = ga + q_full[orig_a][j] * fixed_val[j]
      }
    }
    g_vec[a] = ga
  }

  // 2. Convexity check: active Hessian must be positive semi-definite
  if !is_symmetric_psd(h_mat, n) {
    return Solution::empty(NonConvex, 0, 1)
  }

  // 3. Build active equality system A x = b and inequality system G x <= h
  let a_rows : Array[Array[Double]] = []
  let b_vec : Array[Double] = []
  let eq_user_idx : Array[Int32] = []

  let g_rows : Array[Array[Double]] = []
  let h_rhs : Array[Double] = []
  // ineq_user_idx[k] >= 0 maps to user constraint index; -1 means variable bound row
  let ineq_user_idx : Array[Int32] = []
  // ineq_user_sign[k]: -1.0 for Le (user dual = -z), +1.0 for Ge (user dual = +z)
  let ineq_user_sign : Array[Double] = []

  for i = 0; i < m_user; i = i + 1 {
    let c = model.constraints[i]
    let coeffs = Array::make(n, 0.0)
    let mut rhs = c.rhs
    for k = 0; k < c.terms.length(); k = k + 1 {
      let t = c.terms[k]
      if !is_finite_double(t.coeff) {
        return Solution::empty(NumericalError, 0, 1)
      }
      let idx = t.var_id.index
      if is_fixed[idx] {
        rhs = rhs - t.coeff * fixed_val[idx]
      } else {
        let a = orig_to_act[idx]
        coeffs[a] = coeffs[a] + t.coeff
      }
    }
    let mut has_nz = false
    for a = 0; a < n; a = a + 1 {
      if coeffs[a] != 0.0 {
        has_nz = true
        break
      }
    }
    if !has_nz {
      let ok = match c.relation {
        Le => 0.0 <= rhs + config.tol
        Ge => 0.0 >= rhs - config.tol
        Eq => rhs.abs() <= config.tol
      }
      if !ok {
        return Solution::empty(Infeasible, 0, 1)
      }
      continue
    }
    match c.relation {
      Eq => {
        a_rows.push(coeffs)
        b_vec.push(rhs)
        eq_user_idx.push(i)
      }
      Le => {
        g_rows.push(coeffs)
        h_rhs.push(rhs)
        ineq_user_idx.push(i)
        ineq_user_sign.push(-1.0)
      }
      Ge => {
        for a = 0; a < n; a = a + 1 {
          coeffs[a] = -coeffs[a]
        }
        g_rows.push(coeffs)
        h_rhs.push(-rhs)
        ineq_user_idx.push(i)
        ineq_user_sign.push(1.0)
      }
    }
  }

  // Append finite variable bounds on active variables as inequality rows
  for a = 0; a < n; a = a + 1 {
    let v = model.vars[act_to_orig[a]]
    if v.lb > -inf_threshold {
      let row = Array::make(n, 0.0)
      row[a] = -1.0
      g_rows.push(row)
      h_rhs.push(-v.lb)
      ineq_user_idx.push(-1)
      ineq_user_sign.push(0.0)
    }
    if v.ub < inf_threshold {
      let row = Array::make(n, 0.0)
      row[a] = 1.0
      g_rows.push(row)
      h_rhs.push(v.ub)
      ineq_user_idx.push(-1)
      ineq_user_sign.push(0.0)
    }
  }

  let p = a_rows.length()
  let m = g_rows.length()

  // 4. Exact primal feasibility check & feasible warm start via zero-objective LP
  let feas_model = make_zero_obj_model(model)
  let feas_sol = solve_lp(feas_model, config~)
  if feas_sol.status == Infeasible {
    return Solution::empty(Infeasible, feas_sol.iterations, 1)
  }

  let x0 = Array::make(n, 0.0)
  if feas_sol.status == Optimal && feas_sol.primal.length() == n_orig {
    for a = 0; a < n; a = a + 1 {
      x0[a] = feas_sol.primal[act_to_orig[a]]
    }
  } else {
    for a = 0; a < n; a = a + 1 {
      let v = model.vars[act_to_orig[a]]
      if v.lb > -inf_threshold && v.ub < inf_threshold {
        x0[a] = 0.5 * (v.lb + v.ub)
      } else if v.lb > -inf_threshold {
        x0[a] = v.lb + 1.0
      } else if v.ub < inf_threshold {
        x0[a] = v.ub - 1.0
      } else {
        x0[a] = 0.0
      }
    }
  }

  let mut h_norm = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let a = h_mat[i][j].abs()
      if a > h_norm {
        h_norm = a
      }
    }
  }
  let g_norm = vec_inf_norm(g_vec)
  let b_norm = vec_inf_norm(b_vec)
  let h_rhs_norm = vec_inf_norm(h_rhs)

  // 5. Special fast path: m == 0 (unconstrained or pure equality-constrained QP)
  if m == 0 {
    let k_dim = n + p
    let kkt : Array[Array[Double]] = []
    for _i = 0; _i < k_dim; _i = _i + 1 {
      kkt.push(Array::make(k_dim, 0.0))
    }
    let rhs = Array::make(k_dim, 0.0)
    let reg = 1.0e-13 * (1.0 + h_norm)
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        kkt[i][j] = h_mat[i][j]
      }
      kkt[i][i] = kkt[i][i] + reg
      rhs[i] = -g_vec[i]
    }
    for i = 0; i < p; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        let aij = a_rows[i][j]
        kkt[n + i][j] = aij
        kkt[j][n + i] = aij
      }
      kkt[n + i][n + i] = -1.0e-13
      rhs[n + i] = b_vec[i]
    }
    match lu_factorize(kkt, k_dim) {
      None => return Solution::empty(NumericalError, 1, 1)
      Some((lu, perm)) => {
        let sol_v = lu_solve_refined(kkt, lu, perm, rhs, k_dim)
        let x = Array::make(n, 0.0)
        let y = Array::make(p, 0.0)
        for i = 0; i < n; i = i + 1 {
          x[i] = sol_v[i]
        }
        for i = 0; i < p; i = i + 1 {
          y[i] = sol_v[n + i]
        }
        // Check stationarity & possible unboundedness along nullspace of H & A
        let mut max_rd = 0.0
        for i = 0; i < n; i = i + 1 {
          let mut r = g_vec[i]
          for j = 0; j < n; j = j + 1 {
            r = r + h_mat[i][j] * x[j]
          }
          for k = 0; k < p; k = k + 1 {
            r = r + a_rows[k][i] * y[k]
          }
          if r.abs() > max_rd {
            max_rd = r.abs()
          }
        }
        let x_inf = vec_inf_norm(x)
        if !is_finite_double(x_inf) ||
          x_inf > 1.0e11 ||
          max_rd > config.tol * (1.0 + g_norm + h_norm * x_inf) * 100.0 {
          // Verify recession ray d = (x - x0) / ||x - x0||_inf
          if is_finite_double(x_inf) && x_inf > 1.0e6 {
            let d = Array::make(n, 0.0)
            let mut d_norm = 0.0
            for i = 0; i < n; i = i + 1 {
              d[i] = x[i] - x0[i]
              if d[i].abs() > d_norm {
                d_norm = d[i].abs()
              }
            }
            if d_norm > 0.0 {
              for i = 0; i < n; i = i + 1 {
                d[i] = d[i] / d_norm
              }
              let mut hd_inf = 0.0
              let mut gd = 0.0
              for i = 0; i < n; i = i + 1 {
                gd = gd + g_vec[i] * d[i]
                let mut hdi = 0.0
                for j = 0; j < n; j = j + 1 {
                  hdi = hdi + h_mat[i][j] * d[j]
                }
                if hdi.abs() > hd_inf {
                  hd_inf = hdi.abs()
                }
              }
              let mut ad_inf = 0.0
              for i = 0; i < p; i = i + 1 {
                let mut adi = 0.0
                for j = 0; j < n; j = j + 1 {
                  adi = adi + a_rows[i][j] * d[j]
                }
                if adi.abs() > ad_inf {
                  ad_inf = adi.abs()
                }
              }
              if hd_inf <= 1.0e-5 * (1.0 + h_norm) &&
                ad_inf <= 1.0e-5 &&
                gd < -1.0e-7 {
                return Solution::empty(Unbounded, 1, 1)
              }
            }
          }
          return Solution::empty(NumericalError, 1, 1)
        }
        let primal = fixed_val.copy()
        for a = 0; a < n; a = a + 1 {
          primal[act_to_orig[a]] = x[a]
        }
        let dual = Array::make(m_user, 0.0)
        for i = 0; i < p; i = i + 1 {
          dual[eq_user_idx[i]] = -y[i] * sense_sign
        }
        let rc = compute_qp_reduced_costs(model, q_full, sense_sign, primal, dual)
        return {
          status: Optimal,
          objective: model.eval_objective(primal),
          primal,
          dual,
          reduced_costs: rc,
          iterations: 1,
          nodes_explored: 1,
          mip_gap: 0.0,
        }
      }
    }
  }

  // 6. General Primal-Dual Interior-Point Method (m > 0)
  let x = x0.copy()
  let y = Array::make(p, 0.0)
  let s = Array::make(m, 1.0)
  let z = Array::make(m, 1.0)
  let init_z_scale = if g_norm > 1.0 { 1.0 + 0.1 * g_norm } else { 1.0 }
  for i = 0; i < m; i = i + 1 {
    let mut gx = 0.0
    let grow = g_rows[i]
    for j = 0; j < n; j = j + 1 {
      gx = gx + grow[j] * x[j]
    }
    let gap = h_rhs[i] - gx
    s[i] = if gap > 1.0 { gap } else { 1.0 }
    z[i] = init_z_scale
  }

  let k_dim = n + p + m
  let max_ipm_iters = if config.max_iters < 200 { config.max_iters } else { 200 }
  let stop_tol = if config.tol < 1.0e-10 { config.tol } else { 1.0e-10 }
  let mut iters = 0
  let mut converged = false

  let r_d = Array::make(n, 0.0)
  let r_p = Array::make(p, 0.0)
  let r_g = Array::make(m, 0.0)
  let hx = Array::make(n, 0.0)
  let w_scale = Array::make(m, 1.0)
  let z_tilde = Array::make(m, 1.0)
  let s_tilde = Array::make(m, 1.0)
  let kkt : Array[Array[Double]] = []
  for _i = 0; _i < k_dim; _i = _i + 1 {
    kkt.push(Array::make(k_dim, 0.0))
  }
  let rhs_aff = Array::make(k_dim, 0.0)
  let rhs_corr = Array::make(k_dim, 0.0)
  let ds_aff = Array::make(m, 0.0)
  let ds_corr = Array::make(m, 0.0)

  while iters < max_ipm_iters {
    // Compute residuals & complementarity mu
    for i = 0; i < n; i = i + 1 {
      let mut hxi = 0.0
      let hrow = h_mat[i]
      for j = 0; j < n; j = j + 1 {
        hxi = hxi + hrow[j] * x[j]
      }
      hx[i] = hxi
      r_d[i] = hxi + g_vec[i]
    }
    for i = 0; i < p; i = i + 1 {
      let arow = a_rows[i]
      let mut ax = 0.0
      let yi = y[i]
      for j = 0; j < n; j = j + 1 {
        let aij = arow[j]
        if aij != 0.0 {
          ax = ax + aij * x[j]
          r_d[j] = r_d[j] + aij * yi
        }
      }
      r_p[i] = ax - b_vec[i]
    }
    let mut sz_sum = 0.0
    let mut max_min_sz = 0.0
    for i = 0; i < m; i = i + 1 {
      let grow = g_rows[i]
      let mut gx = 0.0
      let zi = z[i]
      for j = 0; j < n; j = j + 1 {
        let gij = grow[j]
        if gij != 0.0 {
          gx = gx + gij * x[j]
          r_d[j] = r_d[j] + gij * zi
        }
      }
      r_g[i] = gx + s[i] - h_rhs[i]
      sz_sum = sz_sum + s[i] * zi
      let min_sz = if s[i] < zi { s[i] } else { zi }
      if min_sz > max_min_sz {
        max_min_sz = min_sz
      }
    }
    let mu = sz_sum / m.to_double()
    if !is_finite_double(mu) {
      return Solution::empty(NumericalError, iters, 1)
    }

    let x_inf = vec_inf_norm(x)
    let rp_inf = vec_inf_norm(r_p)
    let rg_inf = vec_inf_norm(r_g)
    let rd_inf = vec_inf_norm(r_d)
    let hx_inf = vec_inf_norm(hx)
    let p_res_raw = if rp_inf > rg_inf { rp_inf } else { rg_inf }
    let p_denom = 1.0 +
      (if b_norm > h_rhs_norm {
        if b_norm > x_inf { b_norm } else { x_inf }
      } else if h_rhs_norm > x_inf {
        h_rhs_norm
      } else {
        x_inf
      })
    let d_denom = 1.0 + g_norm + hx_inf
    let p_err = p_res_raw / p_denom
    let d_err = rd_inf / d_denom

    if p_err <= stop_tol &&
      d_err <= stop_tol &&
      mu <= stop_tol &&
      (max_min_sz <= 1.0e-8 || mu <= 1.0e-15) {
      converged = true
      break
    }

    // Detect unbounded convex QP (e.g. rank-deficient H with negative linear slope along ray)
    if x_inf > 1.0e11 && p_err < 1.0e-3 {
      let d = Array::make(n, 0.0)
      let mut d_norm = 0.0
      for i = 0; i < n; i = i + 1 {
        d[i] = x[i] - x0[i]
        if d[i].abs() > d_norm {
          d_norm = d[i].abs()
        }
      }
      if d_norm > 0.0 {
        for i = 0; i < n; i = i + 1 {
          d[i] = d[i] / d_norm
        }
        let mut hd_inf = 0.0
        let mut gd = 0.0
        for i = 0; i < n; i = i + 1 {
          gd = gd + g_vec[i] * d[i]
          let mut hdi = 0.0
          for j = 0; j < n; j = j + 1 {
            hdi = hdi + h_mat[i][j] * d[j]
          }
          if hdi.abs() > hd_inf {
            hd_inf = hdi.abs()
          }
        }
        let mut ad_inf = 0.0
        for i = 0; i < p; i = i + 1 {
          let mut adi = 0.0
          for j = 0; j < n; j = j + 1 {
            adi = adi + a_rows[i][j] * d[j]
          }
          if adi.abs() > ad_inf {
            ad_inf = adi.abs()
          }
        }
        let mut max_g_d = -1.0e300
        for i = 0; i < m; i = i + 1 {
          let mut gdi = 0.0
          for j = 0; j < n; j = j + 1 {
            gdi = gdi + g_rows[i][j] * d[j]
          }
          if gdi > max_g_d {
            max_g_d = gdi
          }
        }
        if hd_inf <= 1.0e-5 * (1.0 + h_norm) &&
          ad_inf <= 1.0e-5 &&
          max_g_d <= 1.0e-5 &&
          gd < -1.0e-7 {
          return Solution::empty(Unbounded, iters, 1)
        }
      }
    }

    // Scale-normalized KKT assembly (never divides by tiny s_i or z_i!)
    for i = 0; i < m; i = i + 1 {
      let w = if s[i] > z[i] { s[i] } else { z[i] }
      let ws = if w > 1.0e-30 { w } else { 1.0e-30 }
      w_scale[i] = ws
      z_tilde[i] = z[i] / ws
      s_tilde[i] = s[i] / ws
    }

    let mut reg_p = 1.0e-13 * (1.0 + h_norm)
    let mut reg_d = 1.0e-13
    let mut lu_opt : (Array[Array[Double]], Array[Int32])? = None

    for _attempt = 0; _attempt < 3; _attempt = _attempt + 1 {
      // Clear KKT matrix
      for i = 0; i < k_dim; i = i + 1 {
        let row = kkt[i]
        for j = 0; j < k_dim; j = j + 1 {
          row[j] = 0.0
        }
      }
      // Block (0,0): H + reg_p I
      for i = 0; i < n; i = i + 1 {
        let row = kkt[i]
        let hrow = h_mat[i]
        for j = 0; j < n; j = j + 1 {
          row[j] = hrow[j]
        }
        row[i] = row[i] + reg_p
      }
      // Blocks (0,1) A^T, (1,0) A, (1,1) -reg_d I
      for i = 0; i < p; i = i + 1 {
        let arow = a_rows[i]
        let r_idx = n + i
        for j = 0; j < n; j = j + 1 {
          let aij = arow[j]
          if aij != 0.0 {
            kkt[j][r_idx] = aij
            kkt[r_idx][j] = aij
          }
        }
        kkt[r_idx][r_idx] = -reg_d
      }
      // Blocks (0,2) G^T, (2,0) Z_tilde G, (2,2) -(S_tilde + reg_d I)
      for i = 0; i < m; i = i + 1 {
        let grow = g_rows[i]
        let r_idx = n + p + i
        let zt = z_tilde[i]
        for j = 0; j < n; j = j + 1 {
          let gij = grow[j]
          if gij != 0.0 {
            kkt[j][r_idx] = gij
            kkt[r_idx][j] = zt * gij
          }
        }
        kkt[r_idx][r_idx] = -(s_tilde[i] + reg_d)
      }

      match lu_factorize(kkt, k_dim) {
        Some(f) => {
          lu_opt = Some(f)
          break
        }
        None => {
          reg_p = reg_p * 1.0e4
          reg_d = reg_d * 1.0e4
        }
      }
    }

    let (lu, perm) = match lu_opt {
      Some(f) => f
      None => break
    }

    // 6a. Predictor (Affine) Step
    for i = 0; i < n; i = i + 1 {
      rhs_aff[i] = -r_d[i]
    }
    for i = 0; i < p; i = i + 1 {
      rhs_aff[n + i] = -r_p[i]
    }
    for i = 0; i < m; i = i + 1 {
      rhs_aff[n + p + i] = -z_tilde[i] * r_g[i] + s[i] * z[i] / w_scale[i]
    }
    let sol_aff = lu_solve_refined(kkt, lu, perm, rhs_aff, k_dim)
    for i = 0; i < m; i = i + 1 {
      let grow = g_rows[i]
      let mut gdx = 0.0
      for j = 0; j < n; j = j + 1 {
        gdx = gdx + grow[j] * sol_aff[j]
      }
      ds_aff[i] = -r_g[i] - gdx
    }
    let dz_aff = Array::make(m, 0.0)
    for i = 0; i < m; i = i + 1 {
      dz_aff[i] = sol_aff[n + p + i]
    }

    let a_p_aff = max_positive_step(s, ds_aff)
    let a_d_aff = max_positive_step(z, dz_aff)
    let mut mu_aff = 0.0
    for i = 0; i < m; i = i + 1 {
      let sp = s[i] + a_p_aff * ds_aff[i]
      let zp = z[i] + a_d_aff * dz_aff[i]
      mu_aff = mu_aff + sp * zp
    }
    mu_aff = mu_aff / m.to_double()

    let mut sigma = if mu > 1.0e-30 && mu_aff > 0.0 {
      let ratio = mu_aff / mu
      let r3 = ratio * ratio * ratio
      if r3 < 1.0e-4 {
        1.0e-4
      } else if r3 > 0.999 {
        0.999
      } else {
        r3
      }
    } else {
      0.1
    }
    if mu_aff > mu {
      sigma = 0.5
    }

    // 6b. Mehrotra Corrector Step (reusing same LU factorization!)
    for i = 0; i < n; i = i + 1 {
      rhs_corr[i] = -r_d[i]
    }
    for i = 0; i < p; i = i + 1 {
      rhs_corr[n + i] = -r_p[i]
    }
    for i = 0; i < m; i = i + 1 {
      let r_sz = s[i] * z[i] + ds_aff[i] * dz_aff[i] - sigma * mu
      rhs_corr[n + p + i] = -z_tilde[i] * r_g[i] + r_sz / w_scale[i]
    }
    let sol_corr = lu_solve_refined(kkt, lu, perm, rhs_corr, k_dim)
    for i = 0; i < m; i = i + 1 {
      let grow = g_rows[i]
      let mut gdx = 0.0
      for j = 0; j < n; j = j + 1 {
        gdx = gdx + grow[j] * sol_corr[j]
      }
      ds_corr[i] = -r_g[i] - gdx
    }
    let dz_corr = Array::make(m, 0.0)
    for i = 0; i < m; i = i + 1 {
      dz_corr[i] = sol_corr[n + p + i]
    }

    let raw_ap = max_positive_step(s, ds_corr)
    let raw_ad = max_positive_step(z, dz_corr)
    let a_p = if 0.995 * raw_ap < 1.0 { 0.995 * raw_ap } else { 1.0 }
    let a_d = if 0.995 * raw_ad < 1.0 { 0.995 * raw_ad } else { 1.0 }

    for i = 0; i < n; i = i + 1 {
      x[i] = x[i] + a_p * sol_corr[i]
    }
    for i = 0; i < p; i = i + 1 {
      y[i] = y[i] + a_d * sol_corr[n + i]
    }
    for i = 0; i < m; i = i + 1 {
      let ns = s[i] + a_p * ds_corr[i]
      let nz = z[i] + a_d * dz_corr[i]
      s[i] = if ns > 1.0e-16 { ns } else { 1.0e-16 }
      z[i] = if nz > 1.0e-16 { nz } else { 1.0e-16 }
    }

    iters = iters + 1
  }

  // Reconstruct primal, dual, reduced costs, and verify final accuracy
  let primal = fixed_val.copy()
  for a = 0; a < n; a = a + 1 {
    let orig_j = act_to_orig[a]
    let v = model.vars[orig_j]
    let mut xj = x[a]
    if v.lb > -inf_threshold && xj < v.lb && v.lb - xj <= config.tol {
      xj = v.lb
    }
    if v.ub < inf_threshold && xj > v.ub && xj - v.ub <= config.tol {
      xj = v.ub
    }
    primal[orig_j] = xj
  }

  let dual = Array::make(m_user, 0.0)
  for i = 0; i < p; i = i + 1 {
    let val = -y[i] * sense_sign
    dual[eq_user_idx[i]] = if val.abs() < 1.0e-15 { 0.0 } else { val }
  }
  for i = 0; i < m; i = i + 1 {
    let ui = ineq_user_idx[i]
    if ui >= 0 {
      let val = ineq_user_sign[i] * z[i] * sense_sign
      dual[ui] = if val.abs() < 1.0e-15 { 0.0 } else { val }
    }
  }

  let rc = compute_qp_reduced_costs(model, q_full, sense_sign, primal, dual)
  let exact_obj = model.eval_objective(primal)
  let sol : Solution = {
    status: Optimal,
    objective: exact_obj,
    primal,
    dual,
    reduced_costs: rc,
    iterations: iters,
    nodes_explored: 1,
    mip_gap: 0.0,
  }

  if !converged {
    // Accept if KKT certificate passes within relaxed user tolerance, else report status honestly
    let cert_tol = if config.tol * 100.0 > 1.0e-6 {
      config.tol * 100.0
    } else {
      1.0e-6
    }
    let cert = model.certify_qp_solution(sol, tol=cert_tol)
    if cert.certified_optimal {
      return sol
    }
    if iters >= config.max_iters {
      return Solution::empty(MaxIterations, iters, 1)
    }
    return Solution::empty(NumericalError, iters, 1)
  }

  let viol = model.max_lp_violation(sol)
  if viol > config.tol * 100.0 && viol > 1.0e-6 {
    return Solution::empty(NumericalError, iters, 1)
  }
  sol
}

///|
fn compute_qp_reduced_costs(
  model : Model,
  q_min : Array[Array[Double]],
  sense_sign : Double,
  primal : Array[Double],
  dual : Array[Double],
) -> Array[Double] {
  let n = model.vars.length()
  let at_y = Array::make(n, 0.0)
  for i = 0; i < model.constraints.length(); i = i + 1 {
    let yi = dual[i]
    if yi != 0.0 {
      let c = model.constraints[i]
      for k = 0; k < c.terms.length(); k = k + 1 {
        let t = c.terms[k]
        at_y[t.var_id.index] = at_y[t.var_id.index] + t.coeff * yi
      }
    }
  }
  let rc : Array[Double] = []
  for j = 0; j < n; j = j + 1 {
    let mut qx = 0.0
    for k = 0; k < n; k = k + 1 {
      qx = qx + q_min[j][k] * sense_sign * primal[k]
    }
    let r = qx + model.vars[j].obj - at_y[j]
    rc.push(if r.abs() < 1.0e-15 { 0.0 } else { r })
  }
  rc
}
