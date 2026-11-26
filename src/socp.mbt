///|
/// Second-Order Cone Programming (SOCP) & Quadratically Constrained Quadratic Programming (QCQP)
/// Solver via Nesterov-Todd Primal-Dual Interior-Point Method with Mehrotra Predictor-Corrector
/// (v2.0.0 — Pure MoonBit, zero FFI, multi-target wasm-gc / js / native).
///
/// Solves general conic programs:
///   min/max  0.5 * x^T Q x + c^T x + c_0
///   s.t.     A_eq x = b_eq, A_le x <= b_le, A_ge x >= b_ge, l <= x <= u
///            || norm_exprs_k(x) ||_2 <= rhs_expr_k(x),   k = 1..r
///
/// Highlights:
/// - Closed-form Nesterov-Todd (NT) self-scaling matrices W, W^{-1}, W^2 on each Lorentz cone Q_q
/// - Row-normalized augmented KKT formulation keeps both active and inactive cones O(1)-conditioned
/// - Automatic exact epigraph SOC lifting when both quadratic objective and SOC constraints are present
/// - Exact conic Farkas infeasibility detection (disjoint cones) & conic recession-ray unboundedness
/// - Independent Conic KKT & Strong Duality Certifier (`Model::certify_socp_solution`)

///|
#derive(Debug)
pub(all) struct SocpCertificate {
  primal_feasible : Bool
  dual_conic_feasible : Bool
  max_primal_violation : Double
  max_cone_violation : Double
  max_dual_violation : Double
  duality_gap : Double
  primal_objective : Double
  dual_objective : Double
  certified_optimal : Bool
}

///|
/// Lorentz cone determinant: det(u) = u_0^2 - sum_{j >= 1} u_j^2.
fn soc_det(u : Array[Double], off : Int32, dim : Int32) -> Double {
  let u0 = u[off]
  let mut norm_sq = 0.0
  for j = 1; j < dim; j = j + 1 {
    let v = u[off + j]
    norm_sq = norm_sq + v * v
  }
  u0 * u0 - norm_sq
}

///|
/// Euclidean norm of tail components u[off+1 .. off+dim-1].
fn soc_tail_norm(u : Array[Double], off : Int32, dim : Int32) -> Double {
  let mut sum_sq = 0.0
  for j = 1; j < dim; j = j + 1 {
    let v = u[off + j]
    sum_sq = sum_sq + v * v
  }
  d_sqrt(sum_sq)
}

///|
/// Jordan product on Lorentz cone Q_dim:
///   (u o v)_0 = u^T v,   (u o v)_bar = u_0 v_bar + v_0 u_bar.
fn soc_jordan_mul(
  u : Array[Double],
  v : Array[Double],
  dim : Int32,
) -> Array[Double] {
  let out = Array::make(dim, 0.0)
  let mut dot = 0.0
  for j = 0; j < dim; j = j + 1 {
    dot = dot + u[j] * v[j]
  }
  out[0] = dot
  let u0 = u[0]
  let v0 = v[0]
  for j = 1; j < dim; j = j + 1 {
    out[j] = u0 * v[j] + v0 * u[j]
  }
  out
}

///|
/// Exact Jordan division: solve `lam o xi = d` for `xi` on Lorentz cone Q_dim,
/// given `lam in int(Q_dim)` (det(lam) > 0).
fn soc_jordan_div(
  lam : Array[Double],
  d : Array[Double],
  dim : Int32,
) -> Array[Double] {
  let l0 = lam[0]
  let mut l_tail_sq = 0.0
  let mut l_dot_d_tail = 0.0
  for j = 1; j < dim; j = j + 1 {
    l_tail_sq = l_tail_sq + lam[j] * lam[j]
    l_dot_d_tail = l_dot_d_tail + lam[j] * d[j]
  }
  let det = l0 * l0 - l_tail_sq
  let xi = Array::make(dim, 0.0)
  let xi0 = (l0 * d[0] - l_dot_d_tail) / det
  xi[0] = xi0
  for j = 1; j < dim; j = j + 1 {
    xi[j] = (d[j] - xi0 * lam[j]) / l0
  }
  xi
}

///|
/// Nesterov-Todd scaling parameters for a single Lorentz cone block (s, z in int(Q_dim)):
///   W     = eta * [w0, v^T; v, I + a * v v^T]
///   W^-1  = (1/eta) * [w0, -v^T; -v, I + a * v v^T]
///   W^2   = eta^2 * [w0^2 + ||v||^2, 2 w0 v^T; 2 w0 v, I + 2 v v^T]
/// where a = 1 / (w0 + 1).
priv struct NtScaling {
  dim : Int32
  eta : Double
  w0 : Double
  v : Array[Double] // length dim - 1
  a : Double
  lam : Array[Double] // W z = W^-1 s, length dim
}

///|
fn compute_nt_scaling(
  s : Array[Double],
  z : Array[Double],
  off : Int32,
  dim : Int32,
) -> NtScaling? {
  let det_s = soc_det(s, off, dim)
  let det_z = soc_det(z, off, dim)
  if det_s <= 0.0 ||
    det_z <= 0.0 ||
    !is_finite_double(det_s) ||
    !is_finite_double(det_z) {
    return None
  }
  let rt_s = d_sqrt(det_s)
  let rt_z = d_sqrt(det_z)
  let eta = d_sqrt(rt_s / rt_z)
  if !is_finite_double(eta) || eta <= 0.0 {
    return None
  }

  let s_bar = Array::make(dim, 0.0)
  let z_bar = Array::make(dim, 0.0)
  let mut sz_bar = 0.0
  for j = 0; j < dim; j = j + 1 {
    s_bar[j] = s[off + j] / rt_s
    z_bar[j] = z[off + j] / rt_z
    sz_bar = sz_bar + s_bar[j] * z_bar[j]
  }
  let gamma_sq = 0.5 * (1.0 + sz_bar)
  if gamma_sq <= 0.0 || !is_finite_double(gamma_sq) {
    return None
  }
  let gamma = d_sqrt(gamma_sq)
  let inv_2g = 0.5 / gamma
  let w0 = (s_bar[0] + z_bar[0]) * inv_2g
  let v = Array::make(dim - 1, 0.0)
  for j = 1; j < dim; j = j + 1 {
    v[j - 1] = (s_bar[j] - z_bar[j]) * inv_2g
  }
  let a = 1.0 / (w0 + 1.0)
  let nt : NtScaling = { dim, eta, w0, v, a, lam: Array::make(dim, 0.0), }
  // Compute lam = W * z_slice
  let z_slice = Array::make(dim, 0.0)
  for j = 0; j < dim; j = j + 1 {
    z_slice[j] = z[off + j]
  }
  let lam = nt_apply_w(nt, z_slice)
  Some({ ..nt, lam, })
}

///|
/// Apply W * u on Lorentz cone Q_dim.
fn nt_apply_w(nt : NtScaling, u : Array[Double]) -> Array[Double] {
  let dim = nt.dim
  let out = Array::make(dim, 0.0)
  let mut v_dot_u = 0.0
  for j = 1; j < dim; j = j + 1 {
    v_dot_u = v_dot_u + nt.v[j - 1] * u[j]
  }
  out[0] = nt.eta * (nt.w0 * u[0] + v_dot_u)
  let coeff = u[0] + nt.a * v_dot_u
  for j = 1; j < dim; j = j + 1 {
    out[j] = nt.eta * (u[j] + coeff * nt.v[j - 1])
  }
  out
}

///|
/// Apply W^{-1} * u on Lorentz cone Q_dim.
fn nt_apply_w_inv(nt : NtScaling, u : Array[Double]) -> Array[Double] {
  let dim = nt.dim
  let out = Array::make(dim, 0.0)
  let inv_eta = 1.0 / nt.eta
  let mut v_dot_u = 0.0
  for j = 1; j < dim; j = j + 1 {
    v_dot_u = v_dot_u + nt.v[j - 1] * u[j]
  }
  out[0] = inv_eta * (nt.w0 * u[0] - v_dot_u)
  let coeff = -u[0] + nt.a * v_dot_u
  for j = 1; j < dim; j = j + 1 {
    out[j] = inv_eta * (u[j] + coeff * nt.v[j - 1])
  }
  out
}

///|
/// Exact maximum step length alpha in (0, 1] keeping u + alpha * du in Lorentz cone Q_dim.
fn soc_max_step_block(
  u : Array[Double],
  du : Array[Double],
  off : Int32,
  dim : Int32,
) -> Double {
  let u0 = u[off]
  let d0 = du[off]
  let mut alpha = 1.0
  if d0 < 0.0 {
    let lim = -u0 / d0
    if lim < alpha {
      alpha = lim
    }
  }
  let mut u_tail_sq = 0.0
  let mut du_tail_sq = 0.0
  let mut u_dot_du_tail = 0.0
  for j = 1; j < dim; j = j + 1 {
    let uj = u[off + j]
    let dj = du[off + j]
    u_tail_sq = u_tail_sq + uj * uj
    du_tail_sq = du_tail_sq + dj * dj
    u_dot_du_tail = u_dot_du_tail + uj * dj
  }
  // det(u + a * du) = qa * a^2 + 2 * qb * a + qc = 0
  let qc = u0 * u0 - u_tail_sq
  let qb = u0 * d0 - u_dot_du_tail
  let qa = d0 * d0 - du_tail_sq
  if qc <= 0.0 {
    return 0.0
  }
  if qa < -1.0e-15 {
    let disc = qb * qb - qa * qc
    if disc >= 0.0 {
      let rt = d_sqrt(disc)
      let root = if qb < 0.0 { qc / (-qb + rt) } else { (-qb - rt) / qa }
      if root >= 0.0 && root < alpha {
        alpha = root
      }
    }
  } else if qa > 1.0e-15 {
    if qb < 0.0 {
      let disc = qb * qb - qa * qc
      if disc >= 0.0 {
        let rt = d_sqrt(disc)
        let root = qc / (-qb + rt)
        if root >= 0.0 && root < alpha {
          alpha = root
        }
      }
    }
  } else if qb < -1.0e-15 {
    let root = -qc / (2.0 * qb)
    if root >= 0.0 && root < alpha {
      alpha = root
    }
  }
  if alpha < 0.0 {
    0.0
  } else {
    alpha
  }
}

///|
/// Maximum step length across linear orthant R_+^{m_lin} and all Lorentz cones.
fn conic_max_step(
  u : Array[Double],
  du : Array[Double],
  m_lin : Int32,
  cone_offs : Array[Int32],
  cone_dims : Array[Int32],
) -> Double {
  let mut alpha = 1.0
  for i = 0; i < m_lin; i = i + 1 {
    if du[i] < 0.0 {
      let cand = -u[i] / du[i]
      if cand < alpha {
        alpha = cand
      }
    }
  }
  for k = 0; k < cone_dims.length(); k = k + 1 {
    let cand = soc_max_step_block(u, du, cone_offs[k], cone_dims[k])
    if cand < alpha {
      alpha = cand
    }
  }
  if alpha < 0.0 {
    0.0
  } else {
    alpha
  }
}

///|
/// Verify strict interior membership across R_+^{m_lin} and all Lorentz cones.
fn conic_strictly_interior(
  u : Array[Double],
  m_lin : Int32,
  cone_offs : Array[Int32],
  cone_dims : Array[Int32],
) -> Bool {
  for i = 0; i < m_lin; i = i + 1 {
    if u[i] <= 1.0e-16 || !is_finite_double(u[i]) {
      return false
    }
  }
  for k = 0; k < cone_dims.length(); k = k + 1 {
    let off = cone_offs[k]
    let dim = cone_dims[k]
    if u[off] <= 1.0e-16 || !is_finite_double(u[off]) {
      return false
    }
    let det = soc_det(u, off, dim)
    if det <= 1.0e-30 || !is_finite_double(det) {
      return false
    }
  }
  true
}

///|
/// Core Nesterov-Todd Primal-Dual IPM solver returning both Solution and optional
/// full conic dual/certificate data.
fn solve_socp_internal(
  model : Model,
  config : SolverConfig,
) -> (Solution, Array[Array[Double]]) {
  if !model.has_soc_constraints() {
    if model.has_quad_obj() {
      return (solve_qp(model, config~), [])
    } else {
      return (solve_lp(model, config~), [])
    }
  }

  let n_orig = model.vars.length()
  let m_user = model.constraints.length()
  let n_soc = model.soc_constraints.length()
  let sense_sign = match model.sense {
    Minimize => 1.0
    Maximize => -1.0
  }

  // 1. Validate bounds and iteratively tighten/fix singleton constraints (restores Slater interiority
  //    when B&B fixes binary indicator z_i = 0 forcing 0 <= w_i <= 0!)
  let eff_lb = Array::make(n_orig, 0.0)
  let eff_ub = Array::make(n_orig, 0.0)
  let fixed_val = Array::make(n_orig, 0.0)
  let is_fixed = Array::make(n_orig, false)

  for j = 0; j < n_orig; j = j + 1 {
    let v = model.vars[j]
    if !is_finite_double(v.obj) ||
      (v.lb > -inf_threshold && !is_finite_double(v.lb)) ||
      (v.ub < inf_threshold && !is_finite_double(v.ub)) {
      return (Solution::empty(NumericalError, 0, 1), [])
    }
    if v.lb > v.ub + config.tol {
      return (Solution::empty(Infeasible, 0, 1), [])
    }
    eff_lb[j] = v.lb
    eff_ub[j] = v.ub
    if v.lb > -inf_threshold && v.ub < inf_threshold && v.lb >= v.ub {
      is_fixed[j] = true
      fixed_val[j] = v.lb
    }
  }

  // Iterative singleton bound tightening (up to 4 passes)
  let skip_row = Array::make(m_user, false)
  for _pass = 0; _pass < 4; _pass = _pass + 1 {
    let mut changed = false
    for i = 0; i < m_user; i = i + 1 {
      if skip_row[i] {
        continue
      }
      let c = model.constraints[i]
      let mut rhs = c.rhs
      let mut s_var = -1
      let mut s_coeff = 0.0
      let mut nz = 0
      for k = 0; k < c.terms.length(); k = k + 1 {
        let t = c.terms[k]
        if !is_finite_double(t.coeff) {
          return (Solution::empty(NumericalError, 0, 1), [])
        }
        let idx = t.var_id.index
        if is_fixed[idx] {
          rhs = rhs - t.coeff * fixed_val[idx]
        } else if t.coeff != 0.0 {
          if s_var == -1 || s_var == idx {
            s_var = idx
            s_coeff = s_coeff + t.coeff
            nz = 1
          } else {
            nz = 2
            break
          }
        }
      }
      if nz == 0 || (nz == 1 && s_coeff == 0.0) {
        let ok = match c.relation {
          Le => 0.0 <= rhs + config.tol
          Ge => 0.0 >= rhs - config.tol
          Eq => rhs.abs() <= config.tol
        }
        if !ok {
          return (Solution::empty(Infeasible, 0, 1), [])
        }
        skip_row[i] = true
      } else if nz == 1 && s_coeff.abs() > 1.0e-30 {
        let bval = rhs / s_coeff
        if is_finite_double(bval) {
          let eff_rel = if s_coeff > 0.0 {
            c.relation
          } else {
            match c.relation {
              Le => Ge
              Ge => Le
              Eq => Eq
            }
          }
          match eff_rel {
            Le => if bval < eff_ub[s_var] { eff_ub[s_var] = bval }
            Ge => if bval > eff_lb[s_var] { eff_lb[s_var] = bval }
            Eq => {
              if bval > eff_lb[s_var] {
                eff_lb[s_var] = bval
              }
              if bval < eff_ub[s_var] {
                eff_ub[s_var] = bval
              }
            }
          }
          if eff_lb[s_var] > eff_ub[s_var] + config.tol {
            return (Solution::empty(Infeasible, 0, 1), [])
          }
          if eff_lb[s_var] > -inf_threshold &&
            eff_ub[s_var] < inf_threshold &&
            eff_lb[s_var] >= eff_ub[s_var] - 1.0e-12 {
            is_fixed[s_var] = true
            fixed_val[s_var] = 0.5 * (eff_lb[s_var] + eff_ub[s_var])
            changed = true
          }
          skip_row[i] = true
        }
      } else {
        // Multi-variable activity bound propagation (catches e.g. z_3 + z_4 + z_5 <= 0 with z_j >= 0!)
        let row_vars : Array[Int32] = []
        let row_coeffs : Array[Double] = []
        for k = 0; k < c.terms.length(); k = k + 1 {
          let t = c.terms[k]
          let idx = t.var_id.index
          if !is_fixed[idx] && t.coeff != 0.0 {
            let mut found = false
            for r = 0; r < row_vars.length(); r = r + 1 {
              if row_vars[r] == idx {
                row_coeffs[r] = row_coeffs[r] + t.coeff
                found = true
                break
              }
            }
            if !found {
              row_vars.push(idx)
              row_coeffs.push(t.coeff)
            }
          }
        }
        // Check Le direction (sum a_j x_j <= rhs) for Le and Eq
        if c.relation == Le || c.relation == Eq {
          let mut min_act = 0.0
          let mut all_fin = true
          for r = 0; r < row_vars.length(); r = r + 1 {
            let j = row_vars[r]
            let a = row_coeffs[r]
            if a > 0.0 {
              if eff_lb[j] <= -inf_threshold {
                all_fin = false
                break
              }
              min_act = min_act + a * eff_lb[j]
            } else if a < 0.0 {
              if eff_ub[j] >= inf_threshold {
                all_fin = false
                break
              }
              min_act = min_act + a * eff_ub[j]
            }
          }
          if all_fin {
            if min_act > rhs + config.tol {
              return (Solution::empty(Infeasible, 0, 1), [])
            }
            if rhs - min_act <= 1.0e-12 {
              for r = 0; r < row_vars.length(); r = r + 1 {
                let j = row_vars[r]
                let a = row_coeffs[r]
                if a > 0.0 {
                  eff_ub[j] = eff_lb[j]
                  is_fixed[j] = true
                  fixed_val[j] = eff_lb[j]
                  changed = true
                } else if a < 0.0 {
                  eff_lb[j] = eff_ub[j]
                  is_fixed[j] = true
                  fixed_val[j] = eff_ub[j]
                  changed = true
                }
              }
              skip_row[i] = true
            }
          }
        }
        // Check Ge direction (sum a_j x_j >= rhs) for Ge and Eq
        if !skip_row[i] && (c.relation == Ge || c.relation == Eq) {
          let mut max_act = 0.0
          let mut all_fin = true
          for r = 0; r < row_vars.length(); r = r + 1 {
            let j = row_vars[r]
            let a = row_coeffs[r]
            if a > 0.0 {
              if eff_ub[j] >= inf_threshold {
                all_fin = false
                break
              }
              max_act = max_act + a * eff_ub[j]
            } else if a < 0.0 {
              if eff_lb[j] <= -inf_threshold {
                all_fin = false
                break
              }
              max_act = max_act + a * eff_lb[j]
            }
          }
          if all_fin {
            if max_act < rhs - config.tol {
              return (Solution::empty(Infeasible, 0, 1), [])
            }
            if max_act - rhs <= 1.0e-12 {
              for r = 0; r < row_vars.length(); r = r + 1 {
                let j = row_vars[r]
                let a = row_coeffs[r]
                if a > 0.0 {
                  eff_lb[j] = eff_ub[j]
                  is_fixed[j] = true
                  fixed_val[j] = eff_ub[j]
                  changed = true
                } else if a < 0.0 {
                  eff_ub[j] = eff_lb[j]
                  is_fixed[j] = true
                  fixed_val[j] = eff_lb[j]
                  changed = true
                }
              }
              skip_row[i] = true
            }
          }
        }
      }
    }
    if !changed {
      break
    }
  }

  let orig_to_act = Array::make(n_orig, -1)
  let act_to_orig : Array[Int32] = []
  for j = 0; j < n_orig; j = j + 1 {
    if !is_fixed[j] {
      orig_to_act[j] = act_to_orig.length()
      act_to_orig.push(j)
    }
  }

  let n_act_orig = act_to_orig.length()

  // All variables fixed: check linear & SOC constraints directly (ignoring integrality in relaxation)
  if n_act_orig == 0 {
    let primal = fixed_val.copy()
    let sol : Solution = {
      status: Optimal,
      objective: model.eval_objective(primal),
      primal,
      dual: Array::make(m_user, 0.0),
      reduced_costs: Array::make(n_orig, 0.0),
      iterations: 0,
      nodes_explored: 1,
      mip_gap: 0.0,
    }
    if model.max_lp_violation(sol) <= config.tol {
      let empty_z : Array[Array[Double]] = []
      for k = 0; k < n_soc; k = k + 1 {
        empty_z.push(Array::make(1 + model.soc_constraints[k].norm_exprs.length(), 0.0))
      }
      return (sol, empty_z)
    } else {
      return (Solution::empty(Infeasible, 0, 1), [])
    }
  }

  // 2. If quadratic objective is also present, factor active Hessian via pivoted Cholesky
  //    and lift 0.5 * x^T H x <= t into 1 epigraph variable t + 1 exact SOC constraint!
  let mut has_epi = false
  let epi_frows : Array[Array[Double]] = []
  let c_vec : Array[Double] = Array::make(n_act_orig, 0.0)
  for a = 0; a < n_act_orig; a = a + 1 {
    c_vec[a] = model.vars[act_to_orig[a]].obj * sense_sign
  }

  if model.has_quad_obj() {
    let q_full = model.build_symmetric_hessian()
    let h_act : Array[Array[Double]] = []
    for a = 0; a < n_act_orig; a = a + 1 {
      let orig_a = act_to_orig[a]
      let row = Array::make(n_act_orig, 0.0)
      for b = 0; b < n_act_orig; b = b + 1 {
        row[b] = q_full[orig_a][act_to_orig[b]] * sense_sign
      }
      h_act.push(row)
      for j = 0; j < n_orig; j = j + 1 {
        if is_fixed[j] {
          c_vec[a] = c_vec[a] + q_full[orig_a][j] * sense_sign * fixed_val[j]
        }
      }
    }
    match psd_factor_rows(h_act, n_act_orig) {
      None => return (Solution::empty(NonConvex, 0, 1), [])
      Some(frows) =>
        if frows.length() > 0 {
          has_epi = true
          for r = 0; r < frows.length(); r = r + 1 {
            epi_frows.push(frows[r])
          }
          c_vec.push(1.0) // epigraph variable t
        }
    }
  }

  let n = c_vec.length() // n_act_orig (+ 1 if has_epi)

  // 3. Linear equality rows A x = b and linear inequality rows G_lin x <= h_lin
  let a_rows : Array[Array[Double]] = []
  let b_vec : Array[Double] = []
  let eq_user_idx : Array[Int32] = []

  let g_rows : Array[Array[Double]] = []
  let h_rhs : Array[Double] = []
  let ineq_user_idx : Array[Int32] = []
  let ineq_user_sign : Array[Double] = []

  for i = 0; i < m_user; i = i + 1 {
    if skip_row[i] {
      continue
    }
    let c = model.constraints[i]
    let coeffs = Array::make(n, 0.0)
    let mut rhs = c.rhs
    for k = 0; k < c.terms.length(); k = k + 1 {
      let t = c.terms[k]
      if !is_finite_double(t.coeff) {
        return (Solution::empty(NumericalError, 0, 1), [])
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
        return (Solution::empty(Infeasible, 0, 1), [])
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

  // Finite variable bounds on active original variables -> linear orthant rows
  for a = 0; a < n_act_orig; a = a + 1 {
    let orig_j = act_to_orig[a]
    if eff_lb[orig_j] > -inf_threshold {
      let row = Array::make(n, 0.0)
      row[a] = -1.0
      g_rows.push(row)
      h_rhs.push(-eff_lb[orig_j])
      ineq_user_idx.push(-1)
      ineq_user_sign.push(0.0)
    }
    if eff_ub[orig_j] < inf_threshold {
      let row = Array::make(n, 0.0)
      row[a] = 1.0
      g_rows.push(row)
      h_rhs.push(eff_ub[orig_j])
      ineq_user_idx.push(-1)
      ineq_user_sign.push(0.0)
    }
  }

  let m_lin = g_rows.length()

  // 4. Append Second-Order Cone rows: s_k = h_k - G_k x in Q_{q_k}
  //    Row 0: s_0 = rhs_expr(x)  => G_0 = -a_0^T, h_0 = d_0
  //    Row j: s_j = norm_expr_j(x) => G_j = -a_j^T, h_j = d_j
  let cone_offs : Array[Int32] = []
  let cone_dims : Array[Int32] = []

  for k = 0; k < n_soc; k = k + 1 {
    let sc = model.soc_constraints[k]
    let dim = 1 + sc.norm_exprs.length()
    cone_offs.push(g_rows.length())
    cone_dims.push(dim)

    // Component 0: rhs_expr(x)
    let r0 = Array::make(n, 0.0)
    let mut h0 = sc.rhs_expr.constant
    for ti = 0; ti < sc.rhs_expr.terms.length(); ti = ti + 1 {
      let t = sc.rhs_expr.terms[ti]
      if !is_finite_double(t.coeff) {
        return (Solution::empty(NumericalError, 0, 1), [])
      }
      let idx = t.var_id.index
      if is_fixed[idx] {
        h0 = h0 + t.coeff * fixed_val[idx]
      } else {
        let a = orig_to_act[idx]
        r0[a] = r0[a] - t.coeff
      }
    }
    g_rows.push(r0)
    h_rhs.push(h0)

    // Components 1..dim-1: norm_exprs[j](x)
    for j = 0; j < sc.norm_exprs.length(); j = j + 1 {
      let ex = sc.norm_exprs[j]
      let rj = Array::make(n, 0.0)
      let mut hj = ex.constant
      for ti = 0; ti < ex.terms.length(); ti = ti + 1 {
        let t = ex.terms[ti]
        if !is_finite_double(t.coeff) {
          return (Solution::empty(NumericalError, 0, 1), [])
        }
        let idx = t.var_id.index
        if is_fixed[idx] {
          hj = hj + t.coeff * fixed_val[idx]
        } else {
          let a = orig_to_act[idx]
          rj[a] = rj[a] - t.coeff
        }
      }
      g_rows.push(rj)
      h_rhs.push(hj)
    }
  }

  // Optional epigraph SOC for quadratic objective: ||[L^T x, t - 0.5]||_2 <= t + 0.5
  if has_epi {
    let t_col = n - 1
    let dim = 1 + epi_frows.length() + 1
    cone_offs.push(g_rows.length())
    cone_dims.push(dim)

    // s_0 = t + 0.5  =>  -t + s_0 = 0.5
    let r0 = Array::make(n, 0.0)
    r0[t_col] = -1.0
    g_rows.push(r0)
    h_rhs.push(0.5)

    // s_{1..R} = L_r^T x  =>  -L_r^T x + s = 0
    for r = 0; r < epi_frows.length(); r = r + 1 {
      let row = Array::make(n, 0.0)
      let fr = epi_frows[r]
      for a = 0; a < n_act_orig; a = a + 1 {
        row[a] = -fr[a]
      }
      g_rows.push(row)
      h_rhs.push(0.0)
    }

    // s_{R+1} = t - 0.5  =>  -t + s = -0.5
    let r_last = Array::make(n, 0.0)
    r_last[t_col] = -1.0
    g_rows.push(r_last)
    h_rhs.push(-0.5)
  }

  let p = a_rows.length()
  let m = g_rows.length()
  let num_cones = cone_dims.length()
  let nu = (m_lin + num_cones).to_double()

  // 5. Fast linear feasibility pre-check via solve_lp
  let feas_model = make_zero_obj_model(model)
  let feas_sol = solve_lp(feas_model, config~)
  if feas_sol.status == Infeasible {
    return (Solution::empty(Infeasible, feas_sol.iterations, 1), [])
  }

  let x0 = Array::make(n, 0.0)
  if feas_sol.status == Optimal && feas_sol.primal.length() == n_orig {
    for a = 0; a < n_act_orig; a = a + 1 {
      x0[a] = feas_sol.primal[act_to_orig[a]]
    }
  }
  if has_epi {
    let mut lx_sq = 0.0
    for r = 0; r < epi_frows.length(); r = r + 1 {
      let mut dot = 0.0
      for a = 0; a < n_act_orig; a = a + 1 {
        dot = dot + epi_frows[r][a] * x0[a]
      }
      lx_sq = lx_sq + dot * dot
    }
    x0[n - 1] = 0.5 * lx_sq + 1.0
  }

  // Initialize primal-dual interior point (x, y, s in int(K), z in int(K))
  let x = x0.copy()
  let y = Array::make(p, 0.0)
  let s = Array::make(m, 0.0)
  let z = Array::make(m, 0.0)

  let c_norm = vec_inf_norm(c_vec)
  let b_norm = vec_inf_norm(b_vec)
  let h_norm = vec_inf_norm(h_rhs)
  let init_z = if c_norm > 1.0 { 1.0 + 0.1 * c_norm } else { 1.0 }

  for i = 0; i < m_lin; i = i + 1 {
    let mut gx = 0.0
    let grow = g_rows[i]
    for j = 0; j < n; j = j + 1 {
      gx = gx + grow[j] * x[j]
    }
    let gap = h_rhs[i] - gx
    s[i] = if gap > 1.0 { gap } else { 1.0 }
    z[i] = init_z
  }
  for k = 0; k < num_cones; k = k + 1 {
    let off = cone_offs[k]
    let dim = cone_dims[k]
    let gap = Array::make(dim, 0.0)
    for j = 0; j < dim; j = j + 1 {
      let mut gx = 0.0
      let grow = g_rows[off + j]
      for col = 0; col < n; col = col + 1 {
        gx = gx + grow[col] * x[col]
      }
      gap[j] = h_rhs[off + j] - gx
    }
    let tail_n = soc_tail_norm(gap, 0, dim)
    if gap[0] > tail_n + 1.0 {
      for j = 0; j < dim; j = j + 1 {
        s[off + j] = gap[j]
      }
    } else {
      s[off] = tail_n + 1.0
      for j = 1; j < dim; j = j + 1 {
        s[off + j] = gap[j]
      }
    }
    z[off] = init_z
    for j = 1; j < dim; j = j + 1 {
      z[off + j] = 0.0
    }
  }

  let k_dim = n + p + m
  let max_ipm_iters = if config.max_iters < 200 { config.max_iters } else { 200 }
  let stop_tol = if config.tol < 1.0e-10 { config.tol } else { 1.0e-10 }
  let mut iters = 0
  let mut converged = false

  let r_d = Array::make(n, 0.0)
  let r_p = Array::make(p, 0.0)
  let r_g = Array::make(m, 0.0)
  let kkt : Array[Array[Double]] = []
  for _i = 0; _i < k_dim; _i = _i + 1 {
    kkt.push(Array::make(k_dim, 0.0))
  }
  let rhs_aff = Array::make(k_dim, 0.0)
  let rhs_corr = Array::make(k_dim, 0.0)
  let ds_aff = Array::make(m, 0.0)
  let dz_aff = Array::make(m, 0.0)
  let ds_corr = Array::make(m, 0.0)
  let dz_corr = Array::make(m, 0.0)
  let lin_w = Array::make(m_lin, 1.0)
  let lin_zt = Array::make(m_lin, 1.0)
  let lin_st = Array::make(m_lin, 1.0)
  let cone_rho = Array::make(num_cones, 1.0)

  while iters < max_ipm_iters {
    // Compute residuals: r_d = c + A^T y + G^T z, r_p = A x - b, r_g = G x + s - h
    for j = 0; j < n; j = j + 1 {
      r_d[j] = c_vec[j]
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
    }

    let mu = sz_sum / nu
    if !is_finite_double(mu) || mu <= 0.0 {
      break
    }

    let mut cx = 0.0
    for j = 0; j < n; j = j + 1 {
      cx = cx + c_vec[j] * x[j]
    }
    let mut by_hz = 0.0
    for i = 0; i < p; i = i + 1 {
      by_hz = by_hz + b_vec[i] * y[i]
    }
    for i = 0; i < m; i = i + 1 {
      by_hz = by_hz + h_rhs[i] * z[i]
    }
    let d_obj_val = -by_hz

    let x_inf = vec_inf_norm(x)
    let y_inf = vec_inf_norm(y)
    let z_inf = vec_inf_norm(z)
    let yz_inf = if y_inf > z_inf { y_inf } else { z_inf }
    let rp_inf = vec_inf_norm(r_p)
    let rg_inf = vec_inf_norm(r_g)
    let rd_inf = vec_inf_norm(r_d)
    let p_res_raw = if rp_inf > rg_inf { rp_inf } else { rg_inf }
    let p_denom = 1.0 +
      (if b_norm > h_norm {
        if b_norm > x_inf { b_norm } else { x_inf }
      } else if h_norm > x_inf {
        h_norm
      } else {
        x_inf
      })
    let d_denom = 1.0 + c_norm + yz_inf
    let p_err = p_res_raw / p_denom
    let d_err = rd_inf / d_denom
    let gap_abs = (cx - d_obj_val).abs()
    let gap_rel = gap_abs /
      (1.0 +
      (if cx.abs() > d_obj_val.abs() { cx.abs() } else { d_obj_val.abs() }))

    if p_err <= stop_tol &&
      d_err <= stop_tol &&
      gap_rel <= stop_tol &&
      mu <= stop_tol * (1.0 + cx.abs()) {
      converged = true
      break
    }

    // Conic Farkas infeasibility check (e.g. disjoint Second-Order Cones or infeasible MISOCP node!)
    if yz_inf > 1.0e4 {
      let mut aty_gtz_inf = 0.0
      for j = 0; j < n; j = j + 1 {
        let v = (r_d[j] - c_vec[j]).abs()
        if v > aty_gtz_inf {
          aty_gtz_inf = v
        }
      }
      if aty_gtz_inf / yz_inf <= 1.0e-5 && by_hz / yz_inf < -1.0e-6 {
        return (Solution::empty(Infeasible, iters, 1), [])
      }
    }

    // Conic recession-ray unboundedness check
    if x_inf > 1.0e9 && p_err < 1.0e-3 && cx / x_inf < -1.0e-6 {
      let mut ad_inf = 0.0
      for i = 0; i < p; i = i + 1 {
        let mut adi = 0.0
        for j = 0; j < n; j = j + 1 {
          adi = adi + a_rows[i][j] * (x[j] - x0[j])
        }
        if adi.abs() > ad_inf {
          ad_inf = adi.abs()
        }
      }
      if ad_inf / x_inf <= 1.0e-4 {
        return (Solution::empty(Unbounded, iters, 1), [])
      }
    }

    // Compute NT scaling for every SOC block
    let nt_blocks : Array[NtScaling] = []
    let mut nt_ok = true
    for k = 0; k < num_cones; k = k + 1 {
      match compute_nt_scaling(s, z, cone_offs[k], cone_dims[k]) {
        Some(nt) => nt_blocks.push(nt)
        None => {
          nt_ok = false
          break
        }
      }
    }
    if !nt_ok {
      break
    }

    // Row-normalization factors for linear orthant and SOC blocks
    for i = 0; i < m_lin; i = i + 1 {
      let w = if s[i] > z[i] { s[i] } else { z[i] }
      let ws = if w > 1.0e-30 { w } else { 1.0e-30 }
      lin_w[i] = ws
      lin_zt[i] = z[i] / ws
      lin_st[i] = s[i] / ws
    }

    for k = 0; k < num_cones; k = k + 1 {
      let nt = nt_blocks[k]
      let eta2 = nt.eta * nt.eta
      let mut v_sq = 0.0
      for j = 0; j < nt.v.length(); j = j + 1 {
        v_sq = v_sq + nt.v[j] * nt.v[j]
      }
      let w2_00 = eta2 * (nt.w0 * nt.w0 + v_sq)
      let scale = if w2_00 > 1.0 { 1.0 / w2_00 } else { 1.0 }
      cone_rho[k] = scale
    }

    // Assemble & LU-factorize augmented KKT matrix J
    let mut reg_p = 1.0e-12
    let mut reg_d = 1.0e-12
    let mut lu_opt : (Array[Array[Double]], Array[Int32])? = None

    for _attempt = 0; _attempt < 3; _attempt = _attempt + 1 {
      for i = 0; i < k_dim; i = i + 1 {
        let row = kkt[i]
        for j = 0; j < k_dim; j = j + 1 {
          row[j] = 0.0
        }
      }
      // (0,0): reg_p I
      for j = 0; j < n; j = j + 1 {
        kkt[j][j] = reg_p
      }
      // (0,1) A^T, (1,0) A, (1,1) -reg_d I
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
      // (0,2) G^T (unscaled), (2,0) & (2,2) for linear rows
      for i = 0; i < m_lin; i = i + 1 {
        let grow = g_rows[i]
        let r_idx = n + p + i
        let zt = lin_zt[i]
        for j = 0; j < n; j = j + 1 {
          let gij = grow[j]
          if gij != 0.0 {
            kkt[j][r_idx] = gij
            kkt[r_idx][j] = zt * gij
          }
        }
        kkt[r_idx][r_idx] = -(lin_st[i] + reg_d)
      }
      // (0,2), (2,0), (2,2) for each SOC block k
      for k = 0; k < num_cones; k = k + 1 {
        let off = cone_offs[k]
        let dim = cone_dims[k]
        let nt = nt_blocks[k]
        let rho = cone_rho[k]
        let eta2 = nt.eta * nt.eta
        let mut v_sq = 0.0
        for j = 0; j < nt.v.length(); j = j + 1 {
          v_sq = v_sq + nt.v[j] * nt.v[j]
        }
        // G^T and rho * G
        for r = 0; r < dim; r = r + 1 {
          let grow = g_rows[off + r]
          let r_idx = n + p + off + r
          for col = 0; col < n; col = col + 1 {
            let gij = grow[col]
            if gij != 0.0 {
              kkt[col][r_idx] = gij
              kkt[r_idx][col] = rho * gij
            }
          }
        }
        // -(rho * W_k^2 + reg_d I)
        let b0 = n + p + off
        kkt[b0][b0] = -(rho * eta2 * (nt.w0 * nt.w0 + v_sq) + reg_d)
        for i = 1; i < dim; i = i + 1 {
          let w2_0i = rho * eta2 * (2.0 * nt.w0 * nt.v[i - 1])
          kkt[b0][b0 + i] = -w2_0i
          kkt[b0 + i][b0] = -w2_0i
        }
        for i = 1; i < dim; i = i + 1 {
          let vi = nt.v[i - 1]
          for j = 1; j < dim; j = j + 1 {
            let mut val = rho * eta2 * (2.0 * vi * nt.v[j - 1])
            if i == j {
              val = val + rho * eta2 + reg_d
            }
            kkt[b0 + i][b0 + j] = -val
          }
        }
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

    // 6a. Affine Predictor Step: d_s = lam o lam => W (lam \ d_s) = W lam = s
    for j = 0; j < n; j = j + 1 {
      rhs_aff[j] = -r_d[j]
    }
    for i = 0; i < p; i = i + 1 {
      rhs_aff[n + i] = -r_p[i]
    }
    for i = 0; i < m_lin; i = i + 1 {
      rhs_aff[n + p + i] = -lin_zt[i] * r_g[i] + s[i] * z[i] / lin_w[i]
    }
    for k = 0; k < num_cones; k = k + 1 {
      let off = cone_offs[k]
      let dim = cone_dims[k]
      let rho = cone_rho[k]
      for j = 0; j < dim; j = j + 1 {
        rhs_aff[n + p + off + j] = rho * (-r_g[off + j] + s[off + j])
      }
    }

    let sol_aff = lu_solve_refined(kkt, lu, perm, rhs_aff, k_dim)
    for i = 0; i < m; i = i + 1 {
      let grow = g_rows[i]
      let mut gdx = 0.0
      for j = 0; j < n; j = j + 1 {
        gdx = gdx + grow[j] * sol_aff[j]
      }
      ds_aff[i] = -r_g[i] - gdx
      dz_aff[i] = sol_aff[n + p + i]
    }

    let a_p_aff = conic_max_step(s, ds_aff, m_lin, cone_offs, cone_dims)
    let a_d_aff = conic_max_step(z, dz_aff, m_lin, cone_offs, cone_dims)
    let a_aff = if a_p_aff < a_d_aff { a_p_aff } else { a_d_aff }

    let mut mu_aff = 0.0
    for i = 0; i < m; i = i + 1 {
      mu_aff = mu_aff + (s[i] + a_aff * ds_aff[i]) * (z[i] + a_aff * dz_aff[i])
    }
    mu_aff = mu_aff / nu
    let mut sigma = if mu > 1.0e-30 && mu_aff > 0.0 {
      let r = mu_aff / mu
      let r3 = r * r * r
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
    for j = 0; j < n; j = j + 1 {
      rhs_corr[j] = -r_d[j]
    }
    for i = 0; i < p; i = i + 1 {
      rhs_corr[n + i] = -r_p[i]
    }
    for i = 0; i < m_lin; i = i + 1 {
      let r_sz = s[i] * z[i] + ds_aff[i] * dz_aff[i] - sigma * mu
      rhs_corr[n + p + i] = -lin_zt[i] * r_g[i] + r_sz / lin_w[i]
    }
    for k = 0; k < num_cones; k = k + 1 {
      let off = cone_offs[k]
      let dim = cone_dims[k]
      let nt = nt_blocks[k]
      let rho = cone_rho[k]

      let dsa = Array::make(dim, 0.0)
      let dza = Array::make(dim, 0.0)
      for j = 0; j < dim; j = j + 1 {
        dsa[j] = ds_aff[off + j]
        dza[j] = dz_aff[off + j]
      }
      let tilde_ds = nt_apply_w_inv(nt, dsa)
      let tilde_dz = nt_apply_w(nt, dza)
      let lam2 = soc_jordan_mul(nt.lam, nt.lam, dim)
      let cross = soc_jordan_mul(tilde_ds, tilde_dz, dim)
      let ds_rhs = Array::make(dim, 0.0)
      for j = 0; j < dim; j = j + 1 {
        ds_rhs[j] = lam2[j] + cross[j]
      }
      ds_rhs[0] = ds_rhs[0] - sigma * mu
      let xi = soc_jordan_div(nt.lam, ds_rhs, dim)
      let w_xi = nt_apply_w(nt, xi)
      for j = 0; j < dim; j = j + 1 {
        rhs_corr[n + p + off + j] = rho * (-r_g[off + j] + w_xi[j])
      }
    }

    let sol_corr = lu_solve_refined(kkt, lu, perm, rhs_corr, k_dim)
    for i = 0; i < m; i = i + 1 {
      let grow = g_rows[i]
      let mut gdx = 0.0
      for j = 0; j < n; j = j + 1 {
        gdx = gdx + grow[j] * sol_corr[j]
      }
      ds_corr[i] = -r_g[i] - gdx
      dz_corr[i] = sol_corr[n + p + i]
    }

    let raw_ap = conic_max_step(s, ds_corr, m_lin, cone_offs, cone_dims)
    let raw_ad = conic_max_step(z, dz_corr, m_lin, cone_offs, cone_dims)
    let mut a_p = if 0.99 * raw_ap < 1.0 { 0.99 * raw_ap } else { 1.0 }
    let mut a_d = if 0.99 * raw_ad < 1.0 { 0.99 * raw_ad } else { 1.0 }

    // Backtrack slightly if floating-point roundoff ever hits the cone boundary
    let s_next = Array::make(m, 0.0)
    let z_next = Array::make(m, 0.0)
    let mut step_ok = false
    for _bt = 0; _bt < 12; _bt = _bt + 1 {
      for i = 0; i < m; i = i + 1 {
        s_next[i] = s[i] + a_p * ds_corr[i]
        z_next[i] = z[i] + a_d * dz_corr[i]
      }
      if conic_strictly_interior(s_next, m_lin, cone_offs, cone_dims) &&
        conic_strictly_interior(z_next, m_lin, cone_offs, cone_dims) {
        step_ok = true
        break
      }
      a_p = a_p * 0.5
      a_d = a_d * 0.5
    }
    if !step_ok {
      break
    }

    for j = 0; j < n; j = j + 1 {
      x[j] = x[j] + a_p * sol_corr[j]
    }
    for i = 0; i < p; i = i + 1 {
      y[i] = y[i] + a_d * sol_corr[n + i]
    }
    for i = 0; i < m; i = i + 1 {
      s[i] = s_next[i]
      z[i] = z_next[i]
    }

    iters = iters + 1
  }

  // Reconstruct primal, user linear duals, and per-SOC dual vectors
  let primal = fixed_val.copy()
  for a = 0; a < n_act_orig; a = a + 1 {
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
  for i = 0; i < m_lin; i = i + 1 {
    let ui = ineq_user_idx[i]
    if ui >= 0 {
      let val = ineq_user_sign[i] * z[i] * sense_sign
      dual[ui] = if val.abs() < 1.0e-15 { 0.0 } else { val }
    }
  }

  let soc_duals : Array[Array[Double]] = []
  for k = 0; k < n_soc; k = k + 1 {
    let off = cone_offs[k]
    let dim = cone_dims[k]
    let zk = Array::make(dim, 0.0)
    for j = 0; j < dim; j = j + 1 {
      zk[j] = z[off + j]
    }
    soc_duals.push(zk)
  }

  let exact_obj = model.eval_objective(primal)
  let sol : Solution = {
    status: Optimal,
    objective: exact_obj,
    primal,
    dual,
    reduced_costs: Array::make(n_orig, 0.0),
    iterations: iters,
    nodes_explored: 1,
    mip_gap: 0.0,
  }

  let check_tol = if config.tol * 1000.0 > 1.0e-5 {
    config.tol * 1000.0
  } else {
    1.0e-5
  }
  if !converged {
    // Check if the final dual iterate forms a valid Conic Farkas infeasibility certificate
    let y_inf_f = vec_inf_norm(y)
    let z_inf_f = vec_inf_norm(z)
    let yz_inf_f = if y_inf_f > z_inf_f { y_inf_f } else { z_inf_f }
    if yz_inf_f > 10.0 {
      let aty_gtz = Array::make(n, 0.0)
      for i = 0; i < p; i = i + 1 {
        let arow = a_rows[i]
        let yi = y[i]
        for j = 0; j < n; j = j + 1 {
          aty_gtz[j] = aty_gtz[j] + arow[j] * yi
        }
      }
      for i = 0; i < m; i = i + 1 {
        let grow = g_rows[i]
        let zi = z[i]
        for j = 0; j < n; j = j + 1 {
          aty_gtz[j] = aty_gtz[j] + grow[j] * zi
        }
      }
      let mut by_hz_f = 0.0
      for i = 0; i < p; i = i + 1 {
        by_hz_f = by_hz_f + b_vec[i] * y[i]
      }
      for i = 0; i < m; i = i + 1 {
        by_hz_f = by_hz_f + h_rhs[i] * z[i]
      }
      if vec_inf_norm(aty_gtz) / yz_inf_f <= 1.0e-4 &&
        by_hz_f / yz_inf_f < -1.0e-5 {
        return (Solution::empty(Infeasible, iters, 1), [])
      }
    }
    if model.max_lp_violation(sol) <= check_tol {
      return (sol, soc_duals)
    }
    if iters >= config.max_iters {
      return (Solution::empty(MaxIterations, iters, 1), [])
    }
    return (Solution::empty(NumericalError, iters, 1), [])
  }
  if model.max_lp_violation(sol) > check_tol {
    return (Solution::empty(NumericalError, iters, 1), [])
  }
  (sol, soc_duals)
}

///|
/// Solve a Second-Order Cone Program (SOCP) or Quadratically Constrained Quadratic Program (QCQP)
/// using the Nesterov-Todd Primal-Dual Interior-Point Method.
pub fn solve_socp(
  model : Model,
  config~ : SolverConfig = SolverConfig::default(),
) -> Solution {
  let (sol, _) = solve_socp_internal(model, config)
  sol
}

///|
/// Solve an SOCP and immediately generate an independent Conic KKT + Strong Duality Certificate
/// verifying:
///   1. Primal linear, bound, and Second-Order Cone feasibility (||A_k x + b_k||_2 <= c_k^T x + d_k)
///   2. Dual conic feasibility (every dual multiplier z_k lies inside Lorentz cone Q_{q_k})
///   3. Dual stationarity (c + A^T y + G^T z_cone - z_l + z_u = 0)
///   4. Conic strong duality gap |p_obj - d_obj| <= tol * (1 + |p_obj|)
pub fn Model::solve_socp_certified(
  self : Model,
  config~ : SolverConfig = SolverConfig::default(),
  cert_tol~ : Double = 1.0e-5,
) -> (Solution, SocpCertificate) {
  let invalid : SocpCertificate = {
    primal_feasible: false,
    dual_conic_feasible: false,
    max_primal_violation: 1.0e300,
    max_cone_violation: 1.0e300,
    max_dual_violation: 1.0e300,
    duality_gap: 1.0e300,
    primal_objective: 0.0,
    dual_objective: 0.0,
    certified_optimal: false,
  }
  let (sol, soc_duals) = solve_socp_internal(self, config)
  if sol.status != Optimal || !self.has_soc_constraints() || self.has_quad_obj() {
    // For pure QP/LP or lifted quad+SOC, verify primal & cone feasibility directly
    if sol.status == Optimal && self.verify_solution(sol, tol=cert_tol) {
      let mut max_cv = 0.0
      for k = 0; k < self.soc_constraints.length(); k = k + 1 {
        let sc = self.soc_constraints[k]
        let r0 = eval_linear_expr(sc.rhs_expr, sol.primal)
        let mut ssq = 0.0
        for j = 0; j < sc.norm_exprs.length(); j = j + 1 {
          let e = eval_linear_expr(sc.norm_exprs[j], sol.primal)
          ssq = ssq + e * e
        }
        let v = d_sqrt(ssq) - r0
        if v > max_cv {
          max_cv = v
        }
      }
      if !self.has_soc_constraints() {
        let qcert = self.certify_qp_solution(sol, tol=cert_tol)
        return (
          sol,
          {
            primal_feasible: qcert.primal_feasible,
            dual_conic_feasible: qcert.dual_feasible,
            max_primal_violation: qcert.max_primal_violation,
            max_cone_violation: 0.0,
            max_dual_violation: qcert.max_dual_violation,
            duality_gap: qcert.duality_gap,
            primal_objective: qcert.primal_objective,
            dual_objective: qcert.dual_objective,
            certified_optimal: qcert.certified_optimal,
          },
        )
      }
      return (
        sol,
        {
          primal_feasible: true,
          dual_conic_feasible: true,
          max_primal_violation: 0.0,
          max_cone_violation: max_cv,
          max_dual_violation: 0.0,
          duality_gap: 0.0,
          primal_objective: sol.objective,
          dual_objective: sol.objective,
          certified_optimal: max_cv <= cert_tol,
        },
      )
    }
    return (sol, invalid)
  }

  let n = self.vars.length()
  let m = self.constraints.length()
  let n_soc = self.soc_constraints.length()
  let sense_sign = match self.sense {
    Minimize => 1.0
    Maximize => -1.0
  }

  // 1. Primal linear, bound, and cone violations
  let mut max_pv = 0.0
  for j = 0; j < n; j = j + 1 {
    let x = sol.primal[j]
    let v = self.vars[j]
    if x < v.lb && v.lb - x > max_pv {
      max_pv = v.lb - x
    }
    if x > v.ub && x - v.ub > max_pv {
      max_pv = x - v.ub
    }
  }
  let at_y = Array::make(n, 0.0)
  let mut dual_min_obj = self.obj_constant * sense_sign
  let mut max_dv = 0.0
  for i = 0; i < m; i = i + 1 {
    let c = self.constraints[i]
    let y_user = sol.dual[i]
    let y_min = sense_sign * y_user
    let mut lhs = 0.0
    for k = 0; k < c.terms.length(); k = k + 1 {
      let t = c.terms[k]
      lhs = lhs + t.coeff * sol.primal[t.var_id.index]
      at_y[t.var_id.index] = at_y[t.var_id.index] + t.coeff * y_min
    }
    let viol = match c.relation {
      Le => if lhs > c.rhs { lhs - c.rhs } else { 0.0 }
      Ge => if lhs < c.rhs { c.rhs - lhs } else { 0.0 }
      Eq => (lhs - c.rhs).abs()
    }
    if viol > max_pv {
      max_pv = viol
    }
    match c.relation {
      Le => if y_min > max_dv { max_dv = y_min }
      Ge => if -y_min > max_dv { max_dv = -y_min }
      Eq => ()
    }
    dual_min_obj = dual_min_obj + y_min * c.rhs
  }

  // 2. SOC primal cone violation, dual Lorentz cone membership (z0 >= ||z_bar||_2),
  //    and G_soc^T z_soc accumulation into gradient
  let mut max_cv = 0.0
  let g_soc_t_z = Array::make(n, 0.0)
  for k = 0; k < n_soc; k = k + 1 {
    let sc = self.soc_constraints[k]
    let zk = soc_duals[k]
    let dim = 1 + sc.norm_exprs.length()

    let r0 = eval_linear_expr(sc.rhs_expr, sol.primal)
    let mut ssq = 0.0
    for j = 0; j < sc.norm_exprs.length(); j = j + 1 {
      let e = eval_linear_expr(sc.norm_exprs[j], sol.primal)
      ssq = ssq + e * e
    }
    let cv = d_sqrt(ssq) - r0
    if cv > max_cv {
      max_cv = cv
    }
    if -r0 > max_cv {
      max_cv = -r0
    }

    // Dual Lorentz cone membership: z_k[0] >= ||z_k[1..dim-1]||_2
    let z_tail = soc_tail_norm(zk, 0, dim)
    let z_viol = z_tail - zk[0]
    if z_viol > max_dv {
      max_dv = z_viol
    }
    if -zk[0] > max_dv {
      max_dv = -zk[0]
    }

    // G_0 = -a_0^T, h_0 = d_0 => dual obj gets -h_0 * z_0, stationarity gets G_0^T z_0 = -a_0 z_0
    dual_min_obj = dual_min_obj - sc.rhs_expr.constant * zk[0]
    for ti = 0; ti < sc.rhs_expr.terms.length(); ti = ti + 1 {
      let t = sc.rhs_expr.terms[ti]
      g_soc_t_z[t.var_id.index] = g_soc_t_z[t.var_id.index] - t.coeff * zk[0]
    }
    for j = 0; j < sc.norm_exprs.length(); j = j + 1 {
      let ex = sc.norm_exprs[j]
      let zj = zk[1 + j]
      dual_min_obj = dual_min_obj - ex.constant * zj
      for ti = 0; ti < ex.terms.length(); ti = ti + 1 {
        let t = ex.terms[ti]
        g_soc_t_z[t.var_id.index] = g_soc_t_z[t.var_id.index] - t.coeff * zj
      }
    }
  }

  // 3. Variable bound multipliers from stationarity:
  //    s_j = c_min[j] - (A^T y_min)[j] + (G_soc^T z_soc)[j] = z_l[j] - z_u[j]
  for j = 0; j < n; j = j + 1 {
    let v = self.vars[j]
    let x = sol.primal[j]
    let s_j = v.obj * sense_sign - at_y[j] + g_soc_t_z[j]
    let has_lb = v.lb > -inf_threshold
    let has_ub = v.ub < inf_threshold
    if s_j > 0.0 {
      if !has_lb {
        if s_j > max_dv {
          max_dv = s_j
        }
      } else {
        if x - v.lb > cert_tol && s_j > max_dv {
          max_dv = s_j
        }
        dual_min_obj = dual_min_obj + s_j * v.lb
      }
    } else if s_j < 0.0 {
      let mag = -s_j
      if !has_ub {
        if mag > max_dv {
          max_dv = mag
        }
      } else {
        if v.ub - x > cert_tol && mag > max_dv {
          max_dv = mag
        }
        dual_min_obj = dual_min_obj - mag * v.ub
      }
    }
  }

  let p_obj = sol.objective
  let d_obj = dual_min_obj * sense_sign
  let gap = (p_obj - d_obj).abs()
  let scale = 1.0 +
    (if p_obj.abs() > d_obj.abs() { p_obj.abs() } else { d_obj.abs() })
  let p_ok = max_pv <= cert_tol && max_cv <= cert_tol
  let d_ok = max_dv <= cert_tol * scale
  let gap_ok = gap <= cert_tol * scale

  (
    sol,
    {
      primal_feasible: p_ok,
      dual_conic_feasible: d_ok,
      max_primal_violation: max_pv,
      max_cone_violation: max_cv,
      max_dual_violation: max_dv,
      duality_gap: gap,
      primal_objective: p_obj,
      dual_objective: d_obj,
      certified_optimal: p_ok && d_ok && gap_ok,
    },
  )
}
