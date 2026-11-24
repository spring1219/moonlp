# `moonlp` — Pure MoonBit LP & MILP Optimization Solver

[![MoonBit CI](https://github.com/spring1219/moonlp/actions/workflows/ci.yml/badge.svg)](https://github.com/spring1219/moonlp/actions/workflows/ci.yml)
![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)
![Backends: wasm-gc | js | native](https://img.shields.io/badge/targets-wasm--gc%20%7C%20js%20%7C%20native-brightgreen.svg)

> **零外部依赖、全后端（`wasm-gc` / `js` / `native`）通用的纯 MoonBit 线性规划（LP）与混合整数线性规划（MILP）求解器与建模库**  
> 填补 `mooncakes.io` 生态中独立运筹优化求解器的空白。

---

## 核心特性（v0.3.0 · 全后端三目标 CI 验证）

- **类型安全的代数建模 DSL** (`src/types.mbt`, `src/model.mbt`)
  - 支持 `Minimize` / `Maximize`、连续 (`Continuous`)、一般整数 (`Integer`)、0-1 二值 (`Binary`) 变量、任意上下界与自由变量，内置独立可行性校验器 `Model::verify_solution`。
- **两阶段原始单纯形 + 对偶单纯形引擎** (`src/simplex.mbt`, `src/dual_simplex.mbt`)
  - **Two-Phase Primal Simplex**：标准型自动变换 + Dantzig 选元 + **Bland 防退化循环规则**，精确输出原始解、对偶影子价格 (`dual`) 与检验数 (`reduced_costs`)。
  - **Dual Simplex (`solve_lp_dual` / `run_dual_simplex`)**：面向对偶可行初基实现零 Phase I 求解，并作为割平面与分支定界的热启动重优化引擎。
- **Gomory 混合整数割平面（GMI Cuts）+ 分支定界 MILP 引擎** (`src/cuts.mbt`, `src/milp.mbt`)
  - **`solve_milp_with_cuts`**：自动识别整数松弛行，从最优 LP 单纯形表生成 **Gomory Mixed-Integer (GMI) 割平面**，通过 **Dual Simplex** 1–2 步热启动重优化；对纯整数和连续+整数混合模型均可显著收紧根节点松弛界（典型算例 1–2 个割平面直接在根节点收敛，无需展开分支树），未完全收敛时无缝衔接 Best-Bound 分支定界保证全局最优。
- **Presolve 预处理器、KKT 强对偶证明器与后验灵敏度分析** (`src/presolve.mbt`, `src/certificate.mbt`, `src/sensitivity.mbt`)
  - 单变量约束界收紧与冗余行消除；独立机检原始可行性误差、互补松弛误差与强对偶间隙 $|c^T x - b^T y|$；约束紧致度/松弛量/影子价格与 RHS 扰动分析。
- **双工业标准格式完整双向读写** (`src/lp_format.mbt`, `src/mps_format.mbt`)
  - **CPLEX `.lp`**：`parse_lp` + `to_lp_string`
  - **Netlib / MIPLIB `.mps`**：`parse_mps` + `to_mps_string`（含 `N/L/G/E` 行、`LO/UP/FX/FR/BV/LI/UI` 边界与 `'MARKER' 'INTORG'/'INTEND'` 整数段）。

---

## 快速开始

```bash
moon check
moon test
moon run src/examples

# 跨后端全量验证
moon test --target wasm-gc
moon test --target js
moon test --target native
```

### 代码示例：Gomory 混合整数割平面（GMI Cuts）+ Dual Simplex

```moonbit
// Maximize 7*x1 + 9*x2 s.t. -x1 + 3*x2 <= 6, 7*x1 + x2 <= 35 (x1, x2 integer)
// LP relaxation = 63.0 -> GMI cuts + Dual Simplex solve directly at root to 55.0!
let m = @moonlp.Model::new(Maximize)
let x1 = m.add_int_var("x1", obj=7.0)
let x2 = m.add_int_var("x2", obj=9.0)
m.add_constraint("c1", [{ var_id: x1, coeff: -1.0 }, { var_id: x2, coeff: 3.0 }], Le, 6.0)
m.add_constraint("c2", [{ var_id: x1, coeff: 7.0 }, { var_id: x2, coeff: 1.0 }], Le, 35.0)

let (sol, rep) = @moonlp.solve_milp_with_cuts(m)
// sol.objective == 55.0, rep.solved_at_root == true, rep.cuts_added == 2
```

---

## 32 项客观数学验证测试清单（`wasm-gc` / `js` / `native` 全绿）

| 类别 | 覆盖测试算例 |
| :--- | :--- |
| **原始 & 对偶单纯形 LP** | 2D 生产规划 + 影子价格、Diet `>=` 最小化（Primal 与 Dual Simplex 双校验）、自由变量与负下界、不可行与无界证明（Primal & Dual）、KKT 强对偶与互补松弛证书、RHS/检验数灵敏度分析 |
| **病态压测与工业 LP** | **Beale 1955 退化循环 LP**（Bland 防循环）、**Klee-Minty 4D 立方体**、CLRS 网络最大流（=23）、3×3 运输问题 LP（=715）、多周期供应链排产 LP（=235）、合金混料 `LinearExpr` 端到端测试 |
| **GMI 割平面 & MILP** | **GMI 割平面根节点收敛（41.25→40、63→55）**、**真·混合整数连续 MILP（1 刀 GMI 割平面 23.75→23.333）**、5 物品与 10 物品 0-1 背包、3×3 指派、混合设施选址、4 周期批量排产 Lot-Sizing MILP（=302）、**64 变量 4×4 数独实时 MILP** |
| **标准格式与预处理** | CPLEX `.lp` 往返解析导出、Netlib `.mps` LP/MILP 双向解析导出往返闭环、非法输入 `ParseError` 容错、Presolve 单变量界收紧与冗余行消除 |

---

## 目录结构

```text
moonlp/
├── moon.mod                    # MoonBit 1.0 模块配置 (v0.3.0)
├── .github/workflows/ci.yml    # wasm-gc / js / native 三后端 CI 矩阵
├── src/
│   ├── moon.pkg                # 主库包配置
│   ├── types.mbt               # 核心类型、LinearExpr 与 IEEE-754 数值工具
│   ├── model.mbt               # 代数建模 DSL 与解验证器
│   ├── simplex.mbt             # Two-Phase Primal Simplex + Bland 防循环
│   ├── dual_simplex.mbt        # Dual Simplex 引擎（零 Phase I & 热启动）
│   ├── milp.mbt                # Best-Bound Branch-and-Bound MILP 引擎
│   ├── cuts.mbt                # Gomory Mixed-Integer (GMI) 割平面 + Branch-and-Cut
│   ├── presolve.mbt            # Presolve 单变量界收紧与冗余行消除
│   ├── certificate.mbt         # KKT 强对偶与互补松弛最优性证明器
│   ├── sensitivity.mbt         # 后验约束/变量灵敏度与 RHS 扰动分析
│   ├── lp_format.mbt           # CPLEX .lp 双向解析器与导出器
│   ├── mps_format.mbt          # Netlib / MIPLIB .mps 双向解析器与导出器
│   ├── moonlp_test.mbt         # 32 项全量单元测试集
│   └── examples/main.mbt       # CLI 综合演示（LP + .lp + 64 变量数独 + GMI 割平面）
└── verify/                     # 独立 Python 数学对拍校验脚本
```

## License
Apache-2.0
