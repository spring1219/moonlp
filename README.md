# `moonlp` — Pure MoonBit LP & MILP Optimization Solver

> **零外部依赖、全后端（Wasm / JS / Native）通用的纯 MoonBit 线性规划（LP）与混合整数线性规划（MILP）求解器与建模库**  
> 填补 `mooncakes.io` 生态中运筹优化求解器（Simplex / MILP / CPLEX `.lp`）的空白。

---

## 核心特性

- **类型安全的代数建模 DSL (`src/model.mbt`, `src/types.mbt`)**
  - 支持 `Minimize` / `Maximize`、连续变量 (`Continuous`)、一般整数变量 (`Integer`)、0-1 二值变量 (`Binary`)。
  - 支持任意变量上下界（包括负下界、单侧无界、自由变量 $x \in (-\infty, +\infty)$、固定变量）与 `<=` / `>=` / `==` 约束。
  - 内置独立解验证器 `Model::verify_solution`，自动校验原始可行性、边界、整数性与目标函数值。
- **两阶段单纯形法引擎 (`src/simplex.mbt`)**
  - 完整实现标准型自动变换与 **Phase I / Phase II Two-Phase Simplex**。
  - **Dantzig 最速下降选元 + 自动切换 Bland 防退化循环规则**：在极端退化算例（如 1955 年 Beale 循环 LP）下保证有限步收敛。
  - 精确输出 **原始最优解 (`primal`)**、**对偶影子价格 (`dual`)**、**检验数 (`reduced_costs`)**，以及严谨的 `Optimal` / `Infeasible` / `Unbounded` 状态证明。
- **分支定界混合整数规划引擎 (`src/milp.mbt`)**
  - 采用 **Best-Bound（最优界优先）Branch-and-Bound** 与 **最分数变量分支策略**。
  - 自动处理二值变量、一般整数变量与连续变量的混合模型，返回可证明的最优间隙 (`mip_gap`) 与探索节点数 (`nodes_explored`)。
- **工业标准 CPLEX `.lp` 格式导出 (`src/lp_format.mbt`)**
  - 支持将任意模型导出为标准 CPLEX `.lp` 文本，方便与 HiGHS / Gurobi / GLPK / SCIP 交叉验算。

---

## 快速开始

```moonbit
test "production mix LP" {
  // Maximize 3x + 2y
  // s.t.  x + y <= 4
  //       2x + y <= 5
  //       x, y >= 0
  let m = Model::new(Maximize)
  let x = m.add_continuous_var("x", obj=3.0)
  let y = m.add_continuous_var("y", obj=2.0)
  m.add_constraint("labor", [{ var_id: x, coeff: 1.0 }, { var_id: y, coeff: 1.0 }], Le, 4.0)
  m.add_constraint("machine", [{ var_id: x, coeff: 2.0 }, { var_id: y, coeff: 1.0 }], Le, 5.0)

  let sol = m.solve()
  // sol.status == Optimal, sol.objective == 9.0, sol.primal == [1.0, 3.0], sol.dual == [1.0, 1.0]
  assert_true(m.verify_solution(sol))
}
```

### 混合整数规划（0-1 背包）

```moonbit
let m = Model::new(Maximize)
let x0 = m.add_binary_var("item0", obj=10.0)
let x1 = m.add_binary_var("item1", obj=14.0)
let x2 = m.add_binary_var("item2", obj=18.0)
m.add_constraint("cap", [
  { var_id: x0, coeff: 3.0 },
  { var_id: x1, coeff: 4.0 },
  { var_id: x2, coeff: 5.0 },
], Le, 9.0)
let sol = m.solve() // Optimal: x1=1, x2=1, obj=32.0
```

---

## 算法验证与基准测试清单

| # | 测试算例 | 类型 | 数学验证项 | 状态 |
| :-: | :--- | :--- | :--- | :-: |
| 1 | 2D Production Mix | LP (Max) | 最优解 + 精确影子价格 + 强对偶 $c^T x = b^T y$ | ✅ 通过 |
| 2 | Diet Problem (`>=`) | LP (Min) | 两阶段 Phase I + 强对偶误差 $< 10^{-7}$ | ✅ 通过 |
| 3 | Free & Negative Bounds | LP (Min) | 自由变量拆解 + 负下界平移 | ✅ 通过 |
| 4 | Contradictory Constraints | LP | 严格 `Infeasible` 证明 | ✅ 通过 |
| 5 | Open Cone | LP | 严格 `Unbounded` 证明 | ✅ 通过 |
| 6 | **Beale Cycling LP (1955)** | 退化病态 LP | Bland 防循环规则 30 步精确收敛 (`-1.25`) | ✅ 通过 |
| 7 | **Klee-Minty 4D Cube** | 经典压测 LP | 精确解 `625.0`（15 次旋转） | ✅ 通过 |
| 8 | 0-1 Knapsack | Binary MILP | Branch-and-Bound 全局最优 (`41.0`) | ✅ 通过 |
| 9 | General Integer Program | Integer MILP | 突破 LP 松弛间隙 (`41.25 -> 40.0`, 9 节点) | ✅ 通过 |
| 10 | 3×3 Assignment Problem | Combinatorial MILP | 等式指派全局最优 (`9.0`) | ✅ 通过 |
| 11 | Facility Location | Mixed Binary+Cont MILP | 选址+运输混合规划全局最优 (`25.0`) | ✅ 通过 |

---

## 构建与运行测试

```bash
# 1. MoonBit 编译与单元测试
moon check
moon test

# 2. 独立数学对拍 Oracle 验证（1-to-1 镜像算法校验）
python verify/verify_oracle.py
```

## 后续路线图（v0.2 规划）
- [ ] 对偶单纯形法（Dual Simplex）热启动，加速 Branch-and-Bound 子节点重优化 5–10 倍
- [ ] Gomory 分数割平面（Cutting Planes）与 Presolve 冗余行/列消去
- [ ] 完整 CPLEX `.lp` 与 `.mps` 工业标准算例解析器
- [ ] 面向 Netlib LP 公开算例集的性能基准（`moon bench`）

## License
Apache-2.0
