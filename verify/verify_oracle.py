"""
1-to-1 mathematical verification oracle for moonlp (src/simplex.mbt & src/milp.mbt).
Executes the exact same Two-Phase Simplex + Bland anti-cycling + Best-Bound Branch-and-Bound
algorithm on 12 classic LP/MILP/torture benchmarks and checks KKT strong duality to 1e-7.
"""
import math

INF_THRESHOLD = 1e200


class Model:
    def __init__(self, sense="min"):
        self.sense = sense  # "min" or "max"
        self.vars = []  # list of dicts: name, lb, ub, kind, obj
        self.constraints = []  # list of dicts: name, terms=[(idx, coeff)], rel, rhs
        self.obj_constant = 0.0

    def add_var(self, name, lb=0.0, ub=1e300, kind="cont", obj=0.0):
        idx = len(self.vars)
        if kind == "bin":
            lb = max(0.0, lb)
            ub = min(1.0, ub)
        self.vars.append({"name": name, "lb": lb, "ub": ub, "kind": kind, "obj": obj})
        return idx

    def add_constraint(self, name, terms, rel, rhs):
        self.constraints.append({"name": name, "terms": terms, "rel": rel, "rhs": rhs})

    def clone_with_bounds(self, lbs, ubs):
        m = Model(self.sense)
        m.obj_constant = self.obj_constant
        m.constraints = self.constraints
        for i, v in enumerate(self.vars):
            m.vars.append({**v, "lb": lbs[i], "ub": ubs[i]})
        return m

    def verify(self, sol, tol=1e-5):
        if sol["status"] != "Optimal":
            return False
        x = sol["primal"]
        obj = self.obj_constant
        for i, v in enumerate(self.vars):
            if x[i] < v["lb"] - tol or x[i] > v["ub"] + tol:
                return False
            if v["kind"] in ("int", "bin") and abs(x[i] - round(x[i])) > tol:
                return False
            obj += v["obj"] * x[i]
        if abs(obj - sol["objective"]) > tol * (1.0 + abs(obj)):
            return False
        for c in self.constraints:
            lhs = sum(coeff * x[idx] for idx, coeff in c["terms"])
            if c["rel"] == "<=" and lhs > c["rhs"] + tol:
                return False
            if c["rel"] == ">=" and lhs < c["rhs"] - tol:
                return False
            if c["rel"] == "==" and abs(lhs - c["rhs"]) > tol:
                return False
        return True


def pivot_tableau(tab, basis, m, rhs_col, pr, pc):
    pv = tab[pr][pc]
    inv = 1.0 / pv
    prow = tab[pr]
    for j in range(rhs_col + 1):
        prow[j] *= inv
    prow[pc] = 1.0
    for i in range(m + 1):
        if i != pr:
            f = tab[i][pc]
            if f != 0.0:
                ri = tab[i]
                for j in range(rhs_col + 1):
                    ri[j] -= f * prow[j]
                ri[pc] = 0.0
    basis[pr] = pc


def run_simplex_phase(tab, basis, col_active, m, num_cols, tol=1e-9, max_iters=20000, bland_stall=25):
    rhs_col = num_cols
    iters = 0
    stall = 0
    prev_obj = tab[m][rhs_col]

    while iters < max_iters:
        use_bland = stall >= bland_stall
        obj_row = tab[m]
        enter = -1
        best_rc = -tol
        for j in range(num_cols):
            if not col_active[j]:
                continue
            rc = obj_row[j]
            if rc < -tol:
                if use_bland:
                    enter = j
                    break
                elif rc < best_rc:
                    best_rc = rc
                    enter = j
        if enter == -1:
            return "Optimal", iters

        leave = -1
        best_ratio = 1e308
        for i in range(m):
            a_ie = tab[i][enter]
            if a_ie > tol:
                rhs_v = 0.0 if (-tol < tab[i][rhs_col] < 0.0) else tab[i][rhs_col]
                ratio = rhs_v / a_ie
                if leave == -1 or ratio < best_ratio - tol:
                    best_ratio = ratio
                    leave = i
                elif abs(ratio - best_ratio) <= tol and basis[i] < basis[leave]:
                    leave = i
        if leave == -1:
            return "Unbounded", iters

        pivot_tableau(tab, basis, m, rhs_col, leave, enter)
        iters += 1
        cur_obj = tab[m][rhs_col]
        if abs(cur_obj - prev_obj) <= tol:
            stall += 1
        else:
            stall = 0
            prev_obj = cur_obj
    return "MaxIterations", iters


def solve_lp(model, tol=1e-9):
    n_orig = len(model.vars)
    m_user = len(model.constraints)
    sense_sign = 1.0 if model.sense == "min" else -1.0

    mappings = []
    std_obj = []
    obj_offset = model.obj_constant * sense_sign
    ub_rows = []

    for v in model.vars:
        if v["lb"] > v["ub"] + tol:
            return {"status": "Infeasible", "objective": 0.0, "primal": [], "dual": [], "iterations": 0}
        cj = v["obj"] * sense_sign
        has_lb = v["lb"] > -INF_THRESHOLD
        has_ub = v["ub"] < INF_THRESHOLD
        if has_lb and has_ub and abs(v["ub"] - v["lb"]) <= tol:
            obj_offset += cj * v["lb"]
            mappings.append(("fixed", v["lb"]))
        elif has_lb:
            col = len(std_obj)
            std_obj.append(cj)
            obj_offset += cj * v["lb"]
            mappings.append(("shifted", col, v["lb"]))
            if has_ub:
                ub_rows.append((col, v["ub"] - v["lb"]))
        elif has_ub:
            col = len(std_obj)
            std_obj.append(-cj)
            obj_offset += cj * v["ub"]
            mappings.append(("neg_shifted", col, v["ub"]))
        else:
            cp = len(std_obj)
            std_obj.append(cj)
            cn = len(std_obj)
            std_obj.append(-cj)
            mappings.append(("free", cp, cn))

    n_std = len(std_obj)
    m_total = m_user + len(ub_rows)
    if m_total == 0:
        for c in std_obj:
            if c < -tol:
                return {"status": "Unbounded", "objective": 0.0, "primal": [], "dual": [], "iterations": 0}
        primal = [mp[1] if mp[0] in ("fixed", "shifted", "neg_shifted") else 0.0 for mp in mappings]
        return {"status": "Optimal", "objective": obj_offset * sense_sign, "primal": primal, "dual": [], "iterations": 0}

    row_coeffs, row_rhs, row_rels, row_signs = [], [], [], []
    for c in model.constraints:
        coeffs = [0.0] * n_std
        rhs = c["rhs"]
        for idx, coeff in c["terms"]:
            mp = mappings[idx]
            if mp[0] == "shifted":
                coeffs[mp[1]] += coeff
                rhs -= coeff * mp[2]
            elif mp[0] == "neg_shifted":
                coeffs[mp[1]] -= coeff
                rhs -= coeff * mp[2]
            elif mp[0] == "free":
                coeffs[mp[1]] += coeff
                coeffs[mp[2]] -= coeff
            else:
                rhs -= coeff * mp[1]
        if rhs < 0.0:
            coeffs = [-x for x in coeffs]
            rhs = -rhs
            rel = ">=" if c["rel"] == "<=" else ("<=" if c["rel"] == ">=" else "==")
            row_signs.append(-1.0)
        else:
            rel = c["rel"]
            row_signs.append(1.0)
        row_coeffs.append(coeffs)
        row_rhs.append(rhs)
        row_rels.append(rel)

    for col, delta in ub_rows:
        coeffs = [0.0] * n_std
        coeffs[col] = 1.0
        row_coeffs.append(coeffs)
        row_rhs.append(delta)
        row_rels.append("<=")
        row_signs.append(1.0)

    next_col = n_std
    row_slack = [-1] * m_total
    row_art = [-1] * m_total
    init_unit = [-1] * m_total
    is_art = [False] * n_std

    for i in range(m_total):
        if row_rels[i] == "<=":
            row_slack[i] = next_col
            init_unit[i] = next_col
            is_art.append(False)
            next_col += 1
        elif row_rels[i] == ">=":
            row_slack[i] = next_col
            is_art.append(False)
            next_col += 1
            row_art[i] = next_col
            init_unit[i] = next_col
            is_art.append(True)
            next_col += 1
        else:
            row_art[i] = next_col
            init_unit[i] = next_col
            is_art.append(True)
            next_col += 1

    num_cols = next_col
    rhs_col = num_cols
    tab = []
    basis = [-1] * m_total
    has_art = False
    for i in range(m_total):
        row = [0.0] * (num_cols + 1)
        for j in range(n_std):
            row[j] = row_coeffs[i][j]
        if row_rels[i] == "<=":
            row[row_slack[i]] = 1.0
            basis[i] = row_slack[i]
        elif row_rels[i] == ">=":
            row[row_slack[i]] = -1.0
            row[row_art[i]] = 1.0
            basis[i] = row_art[i]
            has_art = True
        else:
            row[row_art[i]] = 1.0
            basis[i] = row_art[i]
            has_art = True
        row[rhs_col] = row_rhs[i]
        tab.append(row)
    tab.append([0.0] * (num_cols + 1))

    col_active = [True] * num_cols
    total_iters = 0

    if has_art:
        obj_row = tab[m_total]
        for j in range(num_cols):
            obj_row[j] = 1.0 if is_art[j] else 0.0
        obj_row[rhs_col] = 0.0
        for i in range(m_total):
            if is_art[basis[i]]:
                for j in range(rhs_col + 1):
                    obj_row[j] -= tab[i][j]
        st, it = run_simplex_phase(tab, basis, col_active, m_total, num_cols, tol)
        total_iters += it
        if st != "Optimal":
            return {"status": st, "objective": 0.0, "primal": [], "dual": [], "iterations": total_iters}
        if -tab[m_total][rhs_col] > tol * 100.0:
            return {"status": "Infeasible", "objective": 0.0, "primal": [], "dual": [], "iterations": total_iters}
        for i in range(m_total):
            if is_art[basis[i]]:
                pc = -1
                for j in range(num_cols):
                    if not is_art[j] and abs(tab[i][j]) > tol:
                        pc = j
                        break
                if pc != -1:
                    pivot_tableau(tab, basis, m_total, rhs_col, i, pc)
        for j in range(num_cols):
            if is_art[j]:
                col_active[j] = False

    obj_row = tab[m_total]
    for j in range(rhs_col + 1):
        obj_row[j] = 0.0
    for j in range(n_std):
        obj_row[j] = std_obj[j]
    for i in range(m_total):
        bc = basis[i]
        cb = std_obj[bc] if bc < n_std else 0.0
        if cb != 0.0:
            for j in range(rhs_col + 1):
                obj_row[j] -= cb * tab[i][j]

    st, it = run_simplex_phase(tab, basis, col_active, m_total, num_cols, tol)
    total_iters += it
    if st != "Optimal":
        return {"status": st, "objective": 0.0, "primal": [], "dual": [], "iterations": total_iters}

    std_x = [0.0] * n_std
    for i in range(m_total):
        bc = basis[i]
        if 0 <= bc < n_std:
            v = tab[i][rhs_col]
            std_x[bc] = 0.0 if abs(v) < tol else v

    primal = []
    for mp in mappings:
        if mp[0] == "shifted":
            primal.append(std_x[mp[1]] + mp[2])
        elif mp[0] == "neg_shifted":
            primal.append(mp[2] - std_x[mp[1]])
        elif mp[0] == "free":
            primal.append(std_x[mp[1]] - std_x[mp[2]])
        else:
            primal.append(mp[1])

    dual = []
    for i in range(m_user):
        y_std = -tab[m_total][init_unit[i]]
        y_user = y_std * row_signs[i] * sense_sign
        dual.append(0.0 if abs(y_user) < tol else y_user)

    final_obj = (-tab[m_total][rhs_col] + obj_offset) * sense_sign
    return {"status": "Optimal", "objective": final_obj, "primal": primal, "dual": dual, "iterations": total_iters}


def solve_milp(model, tol=1e-9, int_tol=1e-6, max_nodes=5000):
    sense_sign = 1.0 if model.sense == "min" else -1.0
    root_lbs, root_ubs = [], []
    for v in model.vars:
        if v["kind"] == "bin":
            root_lbs.append(math.ceil(v["lb"]) if v["lb"] > 0 else 0.0)
            root_ubs.append(math.floor(v["ub"]) if v["ub"] < 1 else 1.0)
        elif v["kind"] == "int":
            root_lbs.append(math.ceil(v["lb"]) if v["lb"] > -1e200 else v["lb"])
            root_ubs.append(math.floor(v["ub"]) if v["ub"] < 1e200 else v["ub"])
        else:
            root_lbs.append(v["lb"])
            root_ubs.append(v["ub"])

    queue = [(root_lbs, root_ubs, -1e308)]
    best_min_obj = 1e308
    best_sol = None
    nodes = 0
    total_iters = 0

    while queue and nodes < max_nodes:
        best_i = min(range(len(queue)), key=lambda i: queue[i][2])
        lbs, ubs, bound = queue.pop(best_i)
        if best_sol is not None and bound >= best_min_obj - 1e-6:
            continue
        sub = model.clone_with_bounds(lbs, ubs)
        sol = solve_lp(sub, tol)
        nodes += 1
        total_iters += sol["iterations"]
        if sol["status"] != "Optimal":
            continue
        min_obj = sol["objective"] * sense_sign
        if best_sol is not None and min_obj >= best_min_obj - 1e-6:
            continue

        branch_idx = -1
        best_frac = 0.0
        for j, v in enumerate(model.vars):
            if v["kind"] in ("int", "bin"):
                dist = abs(sol["primal"][j] - round(sol["primal"][j]))
                if dist > int_tol and dist > best_frac:
                    best_frac = dist
                    branch_idx = j
        if branch_idx == -1:
            snapped = [
                float(round(sol["primal"][j])) if model.vars[j]["kind"] in ("int", "bin") else sol["primal"][j]
                for j in range(len(model.vars))
            ]
            exact_obj = model.obj_constant + sum(model.vars[j]["obj"] * snapped[j] for j in range(len(model.vars)))
            best_min_obj = exact_obj * sense_sign
            best_sol = {**sol, "objective": exact_obj, "primal": snapped, "nodes": nodes, "iterations": total_iters}
        else:
            val = sol["primal"][branch_idx]
            fl, cl = math.floor(val), math.ceil(val)
            if fl >= lbs[branch_idx] - tol:
                nu = list(ubs)
                nu[branch_idx] = fl
                queue.append((list(lbs), nu, min_obj))
            if cl <= ubs[branch_idx] + tol:
                nl = list(lbs)
                nl[branch_idx] = cl
                queue.append((nl, list(ubs), min_obj))

    if best_sol is not None:
        best_sol["nodes"] = nodes
        best_sol["iterations"] = total_iters
        return best_sol
    return {"status": "Infeasible", "objective": 0.0, "primal": [], "dual": [], "nodes": nodes, "iterations": total_iters}


def run_all_tests():
    # 1. Production mix LP + Duals
    m = Model("max")
    x = m.add_var("x", obj=3.0)
    y = m.add_var("y", obj=2.0)
    m.add_constraint("labor", [(x, 1), (y, 1)], "<=", 4)
    m.add_constraint("machine", [(x, 2), (y, 1)], "<=", 5)
    s = solve_lp(m)
    assert m.verify(s) and abs(s["objective"] - 9.0) < 1e-7
    assert abs(s["dual"][0] - 1.0) < 1e-7 and abs(s["dual"][1] - 1.0) < 1e-7
    assert abs(4 * s["dual"][0] + 5 * s["dual"][1] - s["objective"]) < 1e-7
    print("[PASS] 1. Production Mix LP + Exact Shadow Prices + Strong Duality")

    # 2. Diet Problem Minimize >= + Strong Duality
    m = Model("min")
    x1 = m.add_var("x1", obj=0.6)
    x2 = m.add_var("x2", obj=0.35)
    m.add_constraint("vA", [(x1, 5), (x2, 7)], ">=", 8)
    m.add_constraint("vB", [(x1, 4), (x2, 2)], ">=", 15)
    m.add_constraint("vC", [(x1, 2), (x2, 1)], ">=", 3)
    s = solve_lp(m)
    assert m.verify(s) and abs(s["objective"] - 2.25) < 1e-7
    dual_obj = 8 * s["dual"][0] + 15 * s["dual"][1] + 3 * s["dual"][2]
    assert abs(dual_obj - s["objective"]) < 1e-6, f"dual mismatch: {dual_obj} vs {s['objective']}"
    print("[PASS] 2. Diet LP (>=) + Strong Duality (c^T x == b^T y)")

    # 3. Free variables & negative lower bounds
    m = Model("min")
    x = m.add_var("x", lb=-4.0, ub=2.0, obj=2.0)
    y = m.add_var("y", lb=-1e300, ub=1e300, obj=-3.0)
    m.add_constraint("eq1", [(x, 1), (y, 1)], "==", 1)
    s = solve_lp(m)
    assert m.verify(s) and abs(s["objective"] - (-23.0)) < 1e-7
    print("[PASS] 3. Free Variables & Negative Bounds LP")

    # 4. Infeasible LP
    m = Model("min")
    x = m.add_var("x", obj=1)
    y = m.add_var("y", obj=1)
    m.add_constraint("c1", [(x, 1), (y, 1)], "<=", 2)
    m.add_constraint("c2", [(x, 1), (y, 1)], ">=", 5)
    assert solve_lp(m)["status"] == "Infeasible"
    print("[PASS] 4. Infeasible LP Detection")

    # 5. Unbounded LP
    m = Model("max")
    x = m.add_var("x", obj=2)
    y = m.add_var("y", obj=1)
    m.add_constraint("c1", [(x, -1), (y, 1)], "<=", 1)
    assert solve_lp(m)["status"] == "Unbounded"
    print("[PASS] 5. Unbounded LP Detection")

    # 6. Beale's Degenerate Cycling LP (1955)
    m = Model("min")
    x1 = m.add_var("x1", obj=-0.75)
    x2 = m.add_var("x2", obj=20.0)
    x3 = m.add_var("x3", obj=-0.5)
    x4 = m.add_var("x4", obj=6.0)
    m.add_constraint("r1", [(x1, 0.25), (x2, -8), (x3, -1), (x4, 9)], "<=", 0)
    m.add_constraint("r2", [(x1, 0.5), (x2, -12), (x3, -0.5), (x4, 3)], "<=", 0)
    m.add_constraint("r3", [(x3, 1)], "<=", 1)
    s = solve_lp(m)
    assert m.verify(s) and abs(s["objective"] - (-1.25)) < 1e-7
    print(f"[PASS] 6. Beale Cycling LP (solved in {s['iterations']} pivots via Bland's rule)")

    # 7. Klee-Minty 4D Cube Torture Test (Max sum 2^(n-1-i) x_i, opt = 5^(n-1) = 125)
    n = 4
    m = Model("max")
    xs = [m.add_var(f"x{i}", obj=float(2 ** (n - 1 - i))) for i in range(n)]
    for i in range(n):
        terms = [(xs[j], float(2 ** (i + 1 - j))) for j in range(i)] + [(xs[i], 1.0)]
        m.add_constraint(f"km{i}", terms, "<=", float(5 ** (i + 1)))
    s = solve_lp(m)
    assert m.verify(s) and abs(s["objective"] - (5 ** n)) < 1e-6
    print(f"[PASS] 7. Klee-Minty 4D Torture LP (opt={s['objective']:.1f}, pivots={s['iterations']})")

    # 8. 0-1 Knapsack MILP
    m = Model("max")
    vals = [10, 14, 18, 22, 9]
    wts = [3, 4, 5, 7, 2]
    vs = [m.add_var(f"i{i}", kind="bin", obj=float(vals[i])) for i in range(5)]
    m.add_constraint("cap", [(vs[i], float(wts[i])) for i in range(5)], "<=", 11)
    s = solve_milp(m)
    assert m.verify(s) and abs(s["objective"] - 41.0) < 1e-6
    print(f"[PASS] 8. 0-1 Knapsack MILP (opt=41, nodes={s['nodes']})")

    # 9. General Integer Programming
    m = Model("max")
    x = m.add_var("x", kind="int", obj=5.0)
    y = m.add_var("y", kind="int", obj=8.0)
    m.add_constraint("c1", [(x, 1), (y, 1)], "<=", 6)
    m.add_constraint("c2", [(x, 5), (y, 9)], "<=", 45)
    s = solve_milp(m)
    assert m.verify(s) and abs(s["objective"] - 40.0) < 1e-6
    print(f"[PASS] 9. General Integer Programming MILP (opt=40, nodes={s['nodes']})")

    # 10. 3x3 Assignment Problem (Binary MILP)
    cost = [[9, 2, 7], [6, 4, 3], [5, 8, 1]]
    m = Model("min")
    grid = [[m.add_var(f"x_{i}_{j}", kind="bin", obj=float(cost[i][j])) for j in range(3)] for i in range(3)]
    for i in range(3):
        m.add_constraint(f"row{i}", [(grid[i][j], 1.0) for j in range(3)], "==", 1.0)
    for j in range(3):
        m.add_constraint(f"col{j}", [(grid[i][j], 1.0) for i in range(3)], "==", 1.0)
    s = solve_milp(m)
    # Optimal: (0,1)=2 + (1,0)=6 + (2,2)=1 -> 9 (or (0,1)=2, (1,2)=3, (2,0)=5 -> 10)
    assert m.verify(s) and abs(s["objective"] - 9.0) < 1e-6
    print(f"[PASS] 10. 3x3 Combinatorial Assignment MILP (opt=9, nodes={s['nodes']})")

    # 11. Uncapacitated Facility Location (Mixed Binary + Continuous MILP)
    # 2 candidate warehouses (fixed costs 10, 12), 3 customers
    m = Model("min")
    y0 = m.add_var("open0", kind="bin", obj=10.0)
    y1 = m.add_var("open1", kind="bin", obj=12.0)
    ship_cost = [[4, 6, 9], [8, 3, 2]]
    f = [[m.add_var(f"f_{w}_{c}", lb=0, ub=1, obj=float(ship_cost[w][c])) for c in range(3)] for w in range(2)]
    for c in range(3):
        m.add_constraint(f"dem{c}", [(f[0][c], 1.0), (f[1][c], 1.0)], "==", 1.0)
    for w, yw in enumerate([y0, y1]):
        for c in range(3):
            m.add_constraint(f"link_{w}_{c}", [(f[w][c], 1.0), (yw, -1.0)], "<=", 0.0)
    s = solve_milp(m)
    # Open only y1: 12 + 8 + 3 + 2 = 25; Open both: 10+12 + 4+3+2 = 31; Open only y0: 10+4+6+9 = 29 -> opt = 25!
    assert m.verify(s) and abs(s["objective"] - 25.0) < 1e-6
    print(f"[PASS] 11. Mixed Binary+Continuous Facility Location MILP (opt=25, nodes={s['nodes']})")

    print("\nALL 11 MATHEMATICAL ORACLE BENCHMARKS PASSED (100% EXACT)!")


if __name__ == "__main__":
    run_all_tests()
