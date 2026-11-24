# `moonlp` — Pure MoonBit LP & MILP Optimization Solver

[![MoonBit CI](https://github.com/spring1219/moonlp/actions/workflows/ci.yml/badge.svg)](https://github.com/spring1219/moonlp/actions/workflows/ci.yml)
![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)
![Backends: wasm-gc | js | native](https://img.shields.io/badge/targets-wasm--gc%20%7C%20js%20%7C%20native-brightgreen.svg)

> **零外部依赖、全后端（`wasm-gc` / `js` / `native`）通用的纯 MoonBit 线性规划（LP）与混合整数线性规划（MILP）求解器与建模库**  
> 已通过真实国际标准 **Netlib `AFIRO`**、**HiGHS `AVGAS`** 与 36 项全后端测试及 `moon bench --release` 基准验证。

---

## 核心特性（v0.4.0 · 三后端 CI + `moon bench` 全绿验证）

- **类型安全的代数建模 DSL** (`src/types.mbt`, `src/model.mbt`)
  - 支持 `Minimize` / `Maximize`、连续 (`Continuous`)、一般整数 (`Integer`)、0-1 二值 (`Binary`) 变量、任意上下界与自由变量，内置独立可行性校验器 `Model::verify_solution`。
- **稀疏感知原始单纯形 + 对偶单纯形引擎** (`src/simplex.mbt`, `src/dual_simplex.mbt`)
  - **稀疏感知 Gauss-Jordan 旋转**：先收集主元行非零列索引，跳过稀疏矩阵中的大量零元素更新，大幅降低工业算例旋转开销。
  - **Harris 最大主元稳定性选元 + Bland 防退化循环规则**：常规迭代中在最小比值平局时优先选取最大主元 $|a_{re}|$ 提升数值条件数，检测到停滞后自动切换 Bland 规则保证有限步收敛。
  - **Dual Simplex (`solve_lp_dual` / `run_dual_simplex`)**：零 Phase I 直接求解对偶可行 LP，并驱动割平面与分支定界子节点热启动。
- **三大 MILP 求解引擎（支持同算例横向对拍）** (`src/warm_bb.mbt`, `src/cuts.mbt`, `src/milp.mbt`)
  - **`solve_milp_warm`（Dual Simplex 热启动分支定界）**：子节点直接继承父节点最优单纯形表，追加单行界割并用 Dual Simplex 1–3 步重优化，彻底避免每节点重跑 Phase I（在 10 物品背包与 5 城 TSP-MTZ 算例上旋转次数显著低于冷启动 B&B）。
  - **`solve_milp_with_cuts`（Gomory 混合整数割平面 GMI Cuts + Branch-and-Cut）**：自动识别整数松弛行，对纯整数及连续+整数混合模型生成有效 GMI 割平面并用 Dual Simplex 重优化，多个经典算例在根节点 1 个节点直接收敛。
  - **`solve_milp`（Best-Bound 分支定界基线）**：支持任意含无下界自由变量的通用 MILP。
- **Presolve 预处理器、KKT 强对偶证明器与后验灵敏度分析** (`src/presolve.mbt`, `src/certificate.mbt`, `src/sensitivity.mbt`)
  - 单变量约束界收紧与冗余行消除；独立机检原始可行性误差、互补松弛误差与强对偶间隙 $|c^T x - b^T y|$；约束紧致度/松弛量/影子价格与 RHS 扰动分析。
- **双工业标准格式完整双向读写** (`src/lp_format.mbt`, `src/mps_format.mbt`)
  - **CPLEX `.lp`**：`parse_lp` + `to_lp_string`
  - **Netlib / MIPLIB `.mps`**：`parse_mps` + `to_mps_string`（兼容 `.301`、`-.4`、`310.` 等真实 Netlib 数值格式与 `'MARKER' 'INTORG'/'INTEND'` 整数段）。

---

## 快速开始

```bash
moon check
moon test
moon run src/examples
moon bench --release

# 跨后端全量验证（与 GitHub Actions CI 一致）
moon test --target wasm-gc && moon bench --release --target wasm-gc
moon test --target js      && moon bench --release --target js
moon test --target native  && moon bench --release --target native
```

### 示例：求解真实国际标准 Netlib `AFIRO` 与 5 城旅行商 TSP-MTZ

```moonbit
// 1. 解析并求解真实 Stanford SOL / Netlib AFIRO 基准算例 (28 rows, 32 cols)
let afiro = @moonlp.parse_mps(@moonlp.netlib_afiro_mps())
let sol = afiro.solve()              // -464.75314286 (canonical Netlib optimum)
let cert = afiro.certify_lp_solution(sol, tol=1.0e-5) // certified_optimal == true

// 2. 使用 Dual Simplex 热启动分支定界求解 5 城 TSP (MTZ 消去子环路 MILP)
let tsp = @moonlp.build_tsp5_mtz_model()
let tsp_sol = @moonlp.solve_milp_warm(tsp) // Optimal tour = 80.0
```

---

## 36 项客观数学验证测试清单（`wasm-gc` / `js` / `native` 全绿）

| 类别 | 覆盖测试算例 |
| :--- | :--- |
| **国际标准工业基准** | **真实 Netlib `AFIRO`**（28×32，最优值 `-464.75314286` + KKT 强对偶认证）、**真实 HiGHS `AVGAS`**（11×8，最优值 `-7.75`） |
| **原始 & 对偶单纯形 LP** | 2D 生产规划 + 影子价格、Diet `>=` 最小化（Primal 与 Dual Simplex 双校验）、自由变量与负下界、不可行与无界证明（Primal & Dual）、KKT 强对偶证书、RHS/检验数灵敏度分析 |
| **病态压测与图论 LP** | **Beale 1955 退化循环 LP**（Bland 防循环）、**Klee-Minty 4D 立方体**、CLRS 网络最大流（=23）、3×3 运输问题 LP（=715）、多周期供应链排产 LP（=235）、合金混料 `LinearExpr` |
| **热启动 B&B / GMI 割平面 MILP** | **5 城旅行商 TSP-MTZ MILP（=80，验证热启动旋转数严格优于冷启动）**、**GMI 割平面根节点收敛（41.25→40、63→55）**、**真·混合整数连续 MILP（1 刀 GMI 割平面 23.75→23.333）**、5/10 物品 0-1 背包、3×3 指派、混合设施选址、4 周期批量排产 Lot-Sizing MILP（=302）、**64 变量 4×4 数独实时 MILP** |
| **标准格式与预处理** | CPLEX `.lp` 往返解析导出、Netlib `.mps` LP/MILP 双向解析导出往返闭环、非法输入 `ParseError` 容错、Presolve 界收紧与冗余行消除 |

---

## 目录结构

```text
moonlp/
├── moon.mod                      # MoonBit 1.0 模块配置 (v0.4.0)
├── .github/workflows/ci.yml      # wasm-gc / js / native 三后端测试 + moon bench 矩阵
├── src/
│   ├── moon.pkg                  # 包配置（含 wbtest @bench 依赖）
│   ├── types.mbt                 # 核心类型、LinearExpr 与 IEEE-754 数值工具
│   ├── model.mbt                 # 代数建模 DSL 与解验证器
│   ├── simplex.mbt               # 稀疏感知 Two-Phase Primal Simplex + Harris/Bland
│   ├── dual_simplex.mbt          # Dual Simplex 引擎（零 Phase I & 热启动）
│   ├── milp.mbt                  # Best-Bound Branch-and-Bound 基线引擎
│   ├── warm_bb.mbt               # Dual Simplex 热启动 Branch-and-Bound 引擎
│   ├── cuts.mbt                  # Gomory Mixed-Integer (GMI) 割平面 + Branch-and-Cut
│   ├── presolve.mbt              # Presolve 单变量界收紧与冗余行消除
│   ├── certificate.mbt           # KKT 强对偶与互补松弛最优性证明器
│   ├── sensitivity.mbt           # 后验约束/变量灵敏度与 RHS 扰动分析
│   ├── lp_format.mbt             # CPLEX .lp 双向解析器与导出器
│   ├── mps_format.mbt            # Netlib / MIPLIB .mps 双向解析器与导出器
│   ├── benchmarks.mbt            # 真实 Netlib AFIRO、HiGHS AVGAS 与 TSP-MTZ 算例集
│   ├── moonlp_test.mbt           # 36 项全量单元测试集
│   ├── moonlp_bench_wbtest.mbt   # 官方 moon bench --release 性能基准套件
│   └── examples/main.mbt         # 6 合 1 CLI 演示（含 Netlib AFIRO 与 5 城 TSP-MTZ）
└── verify/                       # 独立 Python 数学对拍校验脚本
```

## License
Apache-2.0
