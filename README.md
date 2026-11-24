# `moonlp` — Pure MoonBit LP & MILP Optimization Solver

[![MoonBit CI](https://github.com/spring1219/moonlp/actions/workflows/ci.yml/badge.svg)](https://github.com/spring1219/moonlp/actions/workflows/ci.yml)
![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)
![Backends: wasm-gc | js | native](https://img.shields.io/badge/targets-wasm--gc%20%7C%20js%20%7C%20native-brightgreen.svg)

> **零外部依赖、全后端（`wasm-gc` / `js` / `native`）通用的纯 MoonBit 线性规划（LP）与混合整数线性规划（MILP）求解器与建模库**  
> 已通过真实国际标准 **Netlib `ADLITTLE`（56×97）**、**Netlib `AFIRO`（28×32）**、**HiGHS `AVGAS`/`SMALL_MIP`**、**MIPLIB 3.0 `FLUGPL`** 及 **56 项** `--deny-warn` 全后端正确性/回归测试与 `moon bench --release` 验证。

---

## 核心特性（v0.7.0 · MILP 最优性/割平面/数值尺度全面加固 · 三后端 `--deny-warn` 全绿）

- **类型安全的代数建模 DSL 与严格校验** (`src/types.mbt`, `src/model.mbt`)
  - 支持 `Minimize` / `Maximize`、连续 (`Continuous`)、一般整数 (`Integer`)、0-1 二值 (`Binary`) 变量、任意上下界与自由变量；内置拒绝 `NaN`/`±Inf` 的独立可行性校验器 `Model::verify_solution`。
- **稀疏感知原始单纯形 + 对偶单纯形 + 自动 Ruiz 平衡缩放** (`src/simplex.mbt`, `src/dual_simplex.mbt`, `src/scaling.mbt`)
  - **自动 Ruiz 迭代 $L_\infty$ 平衡缩放 (`solve_lp` / `solve_scaled` / `ruiz_scale`)**：当系数跨度超界时默认求解路径自动触发 Ruiz 平衡缩放与原尺度复核，稳健求解 $10^{-10} x \le 1$ 等极端尺度 LP，并精确反缩放原始解、对偶影子价格与检验数。
  - **稀疏感知 Gauss-Jordan 旋转 + $10^{-12}$ 紧凑最小比率测试 + Harris/Bland 防循环**：将主元阈值与比率平局容差独立解耦，避免大尺度差异下的错误出基行。
  - **Dual Simplex (`solve_lp_dual` / `run_dual_simplex`)**：零 Phase I 直接求解对偶可行 LP，并驱动割平面与分支定界子节点热启动。
- **统一状态语义的 Branch-and-Cut 与四大 MILP 引擎** (`src/branch_and_cut.mbt`, `src/warm_bb.mbt`, `src/cuts.mbt`, `src/milp.mbt`)
  - **`solve_milp_branch_and_cut`（旗舰统一引擎）**：融合根节点 **GMI 混合整数割平面** + **混合 DFS 潜行/Best-Bound 搜索** + **继承割平面的 Dual Simplex 热启动分支** + **目标加权分支** + **原始舍入启发式（Primal Rounding Heuristic）**。
  - **严格状态传播与整数候选验收（v0.6.0 加固）**：四大 MILP 入口统一区分 `Optimal` / `Infeasible` / `Unbounded` / `InfeasibleOrUnbounded` / `MaxIterations` / `NodeLimit`——不再将迭代耗尽误判为不可行；利用零目标整数可行性模型（Meyer 定理）区分"LP 松弛无界但整数不可行"与"真正无界"；取整候选统一在原模型上验证/修复连续变量，大系数次容差分数解自动回退细分分支。
- **Presolve 双向对偶恢复、完整 KKT 最优性证书与灵敏度分析** (`src/presolve.mbt`, `src/certificate.mbt`, `src/sensitivity.mbt`)
  - **Presolve + Postsolve (`solve_presolved`)**：单变量约束界收紧与冗余行消除后，自动将保留约束与被消除单变量紧约束的对偶影子价格完整映射回原模型空间。
  - **严格 KKT 最优性证明器 (`certify_lp_solution`)**：完整检查有限性、原始可行性、对偶符号可行性、约束互补松弛、含有限上下界乘子的 KKT 驻点与完整强对偶间隙（彻底杜绝次优点误签）。
- **双工业标准格式无损往返** (`src/lp_format.mbt`, `src/mps_format.mbt`)
  - **CPLEX `.lp`**：`parse_lp` + `to_lp_string`（完整保留目标常数项 `obj_constant`）。
  - **Netlib / MIPLIB `.mps`**：`parse_mps` + `to_mps_string`（完整支持 `OBJSENSE MIN/MAX`、`RHS` 目标常数偏置、固定二值变量 `BV`+`FX` 导出、`N/L/G/E` 行、`RANGES` 区间约束、`LO/UP/FX/FR/BV/LI/UI` 边界与 `'MARKER' 'INTORG'/'INTEND'` 整数段）。

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

## 实测性能基准（`moon bench --release --target wasm-gc`）

以下为 Windows 本机实测耗时（`moon bench --release --target wasm-gc`，10 轮统计均值 ± 标准差）：

| 基准算例 | 规模 / 类型 | 实测均值 (`mean ± σ`) | 区间 (`min … max`) |
| :--- | :--- | ---: | ---: |
| **HiGHS `AVGAS`** | 11×8 LP（38 非零元） | **`4.56 µs`** ± 34.11 ns | `4.51 µs … 4.59 µs` |
| **Netlib `AFIRO`** | 28×32 LP（88 非零元） | **`10.12 µs`** ± 422.24 ns | `9.83 µs … 11.14 µs` |
| **Netlib `ADLITTLE`** | 56×97 LP（480 非零元） | **`654.07 µs`** ± 19.83 µs | `638.61 µs … 704.44 µs` |
| **5 城 TSP-MTZ（冷启动 B&B）** | 20 二值 + 4 连续 MILP | `55.79 µs` ± 382.78 ns | `55.38 µs … 56.39 µs` |
| **5 城 TSP-MTZ（Dual 热启动）** | 20 二值 + 4 连续 MILP | `53.76 µs` ± 1.13 µs | `52.31 µs … 55.19 µs` |
| **5 城 TSP-MTZ（Branch-and-Cut）** | 融合 GMI + 热启动 + 启发式 | **`52.23 µs`** ± 422.26 ns | `51.73 µs … 53.03 µs` |
| **10 物品 0-1 背包（冷启动 B&B）** | 10 二值 MILP | **`62.96 µs`** ± 878.77 ns | `61.75 µs … 64.48 µs` |

---

## 56 项客观数学与正确性回归测试清单（`wasm-gc` / `js` / `native` `--deny-warn` 全绿）

| 类别 | 覆盖测试算例 |
| :--- | :--- |
| **国际标准工业基准** | **真实 Netlib `ADLITTLE`**（56×97，Ruiz 平衡缩放 + KKT 认证）、**真实 Netlib `AFIRO`**（28×32，最优值 `-464.75314286` + KKT 强对偶认证）、**真实 HiGHS `AVGAS` & `SMALL_MIP`**、**真实 MIPLIB 3.0 `FLUGPL`**（根节点 GMI 割平面收紧） |
| **v0.7.0 最优性/割平面回归集（7 项）** | 候选修复不提前关闭松弛未达界节点（$y \le 2\times10^7 x$）、GMI 保留有意义小系数不切掉真最优解（$0.055$）、纯相对 MIP gap 剪枝不剪掉 $10^{-6}$ 微尺度真最优、KKT 证书对 $10^9$ 目标常数尺度不变、$10^{-10}$/$10^{-14}$ 跨入口有界性与紧 Phase I、Presolve 精确来源行对偶分配、`.lp` 严格消费完整 RHS/Bounds |
| **v0.6.0 正确性回归集（8 项）** | `max_iters=1` 状态传播（4 入口统一返回 `MaxIterations`）、LP 无界 vs MILP 不可行/无界消歧、$10^8 x \le 99{,}999{,}950$ 次容差取整细分分支、严格 KKT 拒绝次优点与 `NaN`/`Inf`、$10^{-10} x \le 1$ 自动缩放与紧比率测试、`.lp` 目标常数往返、`.mps` `OBJSENSE MAX`+目标常数+固定二值变量往返、`solve_presolved` 原空间对偶恢复与 KKT 认证 |
| **原始 & 对偶单纯形 LP** | 2D 生产规划 + 影子价格、Diet `>=` 最小化（Primal 与 Dual Simplex 双校验）、自由变量与负下界、不可行与无界证明（Primal & Dual）、KKT 强对偶证书、RHS/检验数灵敏度分析 |
| **病态压测与图论 LP** | **Beale 1955 退化循环 LP**（Bland 防循环）、**Klee-Minty 4D 立方体**、CLRS 网络最大流（=23）、3×3 运输问题 LP（=715）、多周期供应链排产 LP（=235）、合金混料 `LinearExpr` |
| **热启动 B&B / GMI 割平面 MILP** | **5 城旅行商 TSP-MTZ MILP（=80）**、**GMI 割平面根节点收敛（41.25→40、63→55）**、**真·混合整数连续 MILP（1 刀 GMI 割平面 23.75→23.333）**、5/10 物品 0-1 背包、3×3 指派、混合设施选址、4 周期批量排产 Lot-Sizing MILP（=302）、**64 变量 4×4 数独实时 MILP** |
| **标准格式与预处理** | CPLEX `.lp` 往返解析导出、Netlib `.mps` LP/MILP 双向解析导出往返闭环、非法输入 `ParseError` 容错、Presolve 界收紧与冗余行消除 |

---

## 目录结构

```text
moonlp/
├── moon.mod                      # MoonBit 1.0 模块配置 (v0.7.0)
├── .github/workflows/ci.yml      # wasm-gc / js / native 三后端测试 + moon bench 矩阵
├── src/
│   ├── moon.pkg                  # 包配置（含 wbtest @bench 依赖）
│   ├── types.mbt                 # 核心类型、LinearExpr 与 IEEE-754 数值/有限性工具
│   ├── model.mbt                 # 代数建模 DSL 与解验证器（拒绝 NaN/Inf）
│   ├── simplex.mbt               # 稀疏感知 Two-Phase Primal Simplex + 自动 Ruiz 缩放 + Harris/Bland
│   ├── dual_simplex.mbt          # Dual Simplex 引擎（零 Phase I & 热启动）
│   ├── scaling.mbt               # Ruiz 迭代 L-infinity 矩阵平衡缩放与精确反缩放
│   ├── milp.mbt                  # Best-Bound B&B 基线引擎 + 候选修复 + 无界消歧
│   ├── warm_bb.mbt               # Dual Simplex 热启动 Branch-and-Bound 引擎
│   ├── cuts.mbt                  # Gomory Mixed-Integer (GMI) 割平面引擎
│   ├── branch_and_cut.mbt        # 旗舰统一 Branch-and-Cut 引擎（GMI + 热启动 + 舍入启发式）
│   ├── presolve.mbt              # Presolve 单变量界收紧 + Postsolve 原空间对偶恢复
│   ├── certificate.mbt           # 严格 KKT 最优性证明器（含对偶可行性与有限界乘子）
│   ├── sensitivity.mbt           # 后验约束/变量灵敏度与 RHS 扰动分析
│   ├── lp_format.mbt             # CPLEX .lp 双向解析器与导出器（含目标常数）
│   ├── mps_format.mbt            # Netlib / MIPLIB .mps 双向解析器与导出器（含 OBJSENSE/FX/常数）
│   ├── benchmarks.mbt            # 真实 Netlib AFIRO、HiGHS AVGAS 与 TSP-MTZ 算例集
│   ├── benchmarks_netlib_large.mbt # 真实 Netlib ADLITTLE、HiGHS SMALL_MIP 与 MIPLIB FLUGPL
│   ├── moonlp_wbtest.mbt         # 56 项全量单元测试与正确性回归集（--deny-warn 全绿）
│   ├── moonlp_bench_wbtest.mbt   # 官方 moon bench --release 性能基准套件
│   └── examples/main.mbt         # 6 合 1 CLI 演示（含 Netlib AFIRO 与 5 城 TSP-MTZ）
└── verify/                       # 独立 Python 数学对拍校验脚本
```

## License
Apache-2.0
