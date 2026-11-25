# `moonlp` — Pure MoonBit MIQP / Convex QP & LP / MILP Optimization Solver

[![MoonBit CI](https://github.com/spring1219/moonlp/actions/workflows/ci.yml/badge.svg)](https://github.com/spring1219/moonlp/actions/workflows/ci.yml)
![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)
![Backends: wasm-gc | js | native](https://img.shields.io/badge/targets-wasm--gc%20%7C%20js%20%7C%20native-brightgreen.svg)

> **零外部依赖、全后端（`wasm-gc` / `js` / `native`）通用的纯 MoonBit 混合整数二次规划（MIQP）、凸二次规划（QP）与 LP/MILP 求解及证明库**  
> 精准填补 MoonBit 生态在 **混合整数二次规划（MIQP）**、**二阶原始-对偶内点法精确凸 QP（含不可行/无界判定）**、**二次 KKT / Wolfe 强对偶独立证明** 与 **`.qps` (`QUADOBJ`/`QMATRIX`) / CPLEX 二次 `.lp`** 上的空白。

---

## 与 MoonBit 生态已有优化库（`moonopt` / `moonqp`）的分工与实质差异

详细全生态排重报告见 [docs/comparison_with_moonopt.md](docs/comparison_with_moonopt.md)：

| 能力维度 | `Freon793/moonopt`（线性稀疏 LP/MILP） | `yujieMamba/moonqp`（连续 ADMM QP） | **`spring1219/moonlp` v1.0.0（本项目）** |
| :--- | :---: | :---: | :---: |
| **混合整数二次规划（MIQP）** | 不支持 | **明确不支持** | **原生支持**（`solve_miqp`：B&B + 凸 QP 松弛 + 固定整数修复） |
| **凸 QP 求解内核** | 不支持 | 一阶 ADMM（稠密数组） | **二阶尺度归一化 Mehrotra 原始-对偶内点法**（7~15 步至 $10^{-11}$ 精度） |
| **QP 不可行 / 无界 / 非凸精确判定** | 仅 LP | **不支持**（仅 `IterationLimit`） | **原生支持**（`Infeasible` / `Unbounded` / `NonConvex` / `Optimal`） |
| **独立 Wolfe 强对偶 + 二次 KKT 证明** | 仅线性 | 无 Wolfe 强对偶证明 | **原生支持** `certify_qp_solution` |
| **`.qps` (`QUADOBJ`/`QMATRIX`) & CPLEX 二次 `.lp`** | 仅线性 MPS/LP | **不支持** | **原生支持双向无损往返**（并与 `moonopt` 标准格式互通） |
| **大规模稀疏纯线性 LP/MILP** | **强项（推荐使用）** | 不支持 | 兼容支持（详见对比文档） |

---

## 核心新增具体模型与客观验证结果（均经 SciPy + 全组合穷举对拍）

以下模型全部内置于 [`src/benchmarks_qp.mbt`](src/benchmarks_qp.mbt)，经 [`verify/verify_qp_miqp.py`](verify/verify_qp_miqp.py) 独立对拍，全部误差 $\le 1.55\times 10^{-11}$：

### 1. 混合整数二次规划（MIQP）模型 —— `moonopt` 与 `moonqp` 均无法直接求解

| 具体 MIQP 模型 | 实际问题与结构特征 | `moonlp` 全局最优结果 | 穷举+SciPy 误差 | `wasm-gc` 耗时 |
| :--- | :--- | :--- | :---: | ---: |
| **电力系统机组组合经济调度** (`build_unit_commitment_miqp`) | 5 机组二值启停 + 旋转备用 + **二次煤耗曲线** $a_g p_g^2 + b_g p_g$ | 总成本 **`3308.30`**（9 节点，`mip_gap=0`；**Coal 80MW + CCGT 55MW + Peaker2 5MW**） | **`1.55e-11`** | **`1.58 ms`** |
| **基数约束 Markowitz 资产配置** (`build_cardinality_portfolio_miqp`) | 6 选 $\le 3$ 资产 + 最小持仓 10% + $12\times 12$ 秩 6 半正定协方差 | 全局最优方差 **`0.00520192`**（9 节点，精选 **US_Eq 41% / Corp_Bond 49% / Gold 10%**） | **`4.37e-14`** | **`4.34 ms`** |
| **$L_0$ 基数约束稀疏岭回归** (`build_best_subset_ridge_miqp`) | 5 含噪相关特征限选 $\le 2$ 个非零系数（严格 $L_0$ 约束） | 损失 **`0.23809507`**（7 节点，精准恢复真实支撑集 **`{f0, f3}`**） | **`1.59e-13`** | **`1.82 ms`** |
| **推力死区混杂模型预测控制** (`build_hybrid_mpc_miqp`) | 推进器死区 $u_t \in \{0\}\cup[0.3,1.0]$ + 最多点火 3 次（对比连续 MPC） | 最优代价 **`41.228125`**，开关序列 **`[t0:1.0, t1:1.0, t2..4:OFF]`** | **`1.46e-11`** | **`580.05 µs`** |

### 2. 连续凸 QP 模型与 Wolfe 强对偶证明 / 格式往返

| 具体 QP 模型 | `moonlp` 实测最优结果 | Wolfe KKT / SciPy 误差 | `wasm-gc` 耗时 |
| :--- | :--- | :---: | ---: |
| **Markowitz 6 资产有效前沿** ($r \ge 8\%$) | $\frac{1}{2}w^T\Sigma w = \mathbf{0.0049586780}$（8 次内点迭代） | Wolfe=`true` (gap `5.89e-11`)；**`1.47e-12`** | **`56.13 µs`** |
| **8 步饱和/速度约束线性 MPC**（24 变量 16 等式） | LQR 代价 **`38.79579682`**（$u_0=-1.0$ 限幅，$p_8=-0.011$） | Wolfe=`true`；**`5.76e-13`** | **`859.64 µs`** |
| **软间隔 SVM 对偶训练**（$8\times 8$ 秩 2 半正定） | 对偶最优 **`-0.09289364`**（7 次迭代，2 支持向量） | Wolfe=`true`；**`1.44e-11`** | **`73.01 µs`** |
| **单纯形最小二乘 + `.lp` / `.qps` (`QUADOBJ`) 往返** | 原模型、CPLEX 二次 `.lp`、`.qps` 三路径一致 | Wolfe=`true`（三路径全绿） | — |

---

## 快速开始

```bash
moon check --deny-warn
moon test --deny-warn
moon run src/examples
moon bench --release

# 独立 Python / SciPy + 全组合穷举对拍验证全部 7 个 QP/MIQP 模型
python verify/verify_qp_miqp.py

# 跨后端全量验证（与 GitHub Actions CI 一致）
moon test --deny-warn --target wasm-gc && moon bench --release --target wasm-gc
moon test --deny-warn --target js      && moon bench --release --target js
moon test --deny-warn --target native  && moon bench --release --target native
```

### 示例：电力机组组合（MIQP）与推力死区混杂控制（MIQP-MPC）

```moonbit
// 1. 电力系统火电机组组合 + 二次煤耗经济调度 (MIQP)
let uc = @moonlp.build_unit_commitment_miqp(demand=140.0, reserve=20.0)
let uc_sol = uc.solve()        // 3308.30 (Coal 80MW, CCGT 55MW, Peaker2 5MW), ~1.58 ms

// 2. 带最小冲量死区与点火预算的混杂 MPC (MIQP-MPC)
let hmpc = @moonlp.build_hybrid_mpc_miqp(horizon=5, max_firings=3, u_min=0.3, u_max=1.0)
let h_sol = hmpc.solve()       // 41.228125 ([1.0, 1.0, OFF, OFF, OFF]), ~580 us on wasm-gc

// 3. Markowitz 连续投资组合 (QP) + 独立 Wolfe 强对偶 KKT 证明
let qp = @moonlp.build_markowitz_portfolio_qp(target_return=0.08)
let q_sol = qp.solve()         // 0.00495868 (~56 us)
let cert = qp.certify_qp_solution(q_sol) // certified_optimal == true
```

---

## 全量实测性能基准（`moon bench --release --target wasm-gc`）

以下为 Windows 本机实测耗时（`moon bench --release --target wasm-gc`，10 轮统计均值 ± 标准差）：

| 基准算例 | 规模 / 类型 | 实测均值 (`mean ± σ`) | 区间 (`min … max`) |
| :--- | :--- | ---: | ---: |
| **Markowitz 6 资产投资组合（QP）** | 6 变量正定协方差 QP（IPM） | **`56.13 µs`** ± 739.62 ns | `55.28 µs … 57.51 µs` |
| **软间隔 SVM 对偶训练（QP）** | 8×8 秩 2 半正定 Gram QP（IPM） | **`73.01 µs`** ± 2.47 µs | `71.87 µs … 79.85 µs` |
| **推力死区混杂 MPC（MIQP-MPC）** | 5 步死区 + $\le 3$ 次点火 MIQP | **`580.05 µs`** ± 18.52 µs | `566.03 µs … 631.23 µs` |
| **8 步饱和线性 MPC（QP）** | 24 变量 + 16 动力学等式 QP（IPM） | **`859.64 µs`** ± 23.80 µs | `839.58 µs … 913.50 µs` |
| **电力机组组合经济调度（UC-MIQP）** | 5 机组二值启停 + 二次煤耗 MIQP | **`1.58 ms`** ± 47.41 µs | `1.54 ms … 1.70 ms` |
| **$L_0$ 基数约束稀疏岭回归（MIQP）** | 5 选 $\le 2$ 特征子集选择 MIQP | **`1.82 ms`** ± 94.15 µs | `1.75 ms … 2.06 ms` |
| **基数约束 Markowitz 选股（MIQP）** | 6 选 $\le 3$ 资产、秩 6 半正定 MIQP | **`4.34 ms`** ± 46.39 µs | `4.27 ms … 4.41 ms` |
| **HiGHS `AVGAS`（LP）** | 11×8 LP（38 非零元） | **`5.34 µs`** ± 144.16 ns | `5.18 µs … 5.56 µs` |
| **Netlib `AFIRO`（LP）** | 28×32 LP（88 非零元） | **`17.82 µs`** ± 1.49 µs | `16.21 µs … 20.73 µs` |
| **Netlib `ADLITTLE`（LP）** | 56×97 LP（480 非零元） | **`658.86 µs`** ± 23.01 µs | `633.45 µs … 698.54 µs` |
| **5 城 TSP-MTZ（Branch-and-Cut MILP）** | 20 二值 + 4 连续 MILP | **`67.68 µs`** ± 550.37 ns | `66.85 µs … 68.67 µs` |

---

## 核心特性与测试覆盖（80 项单元测试 + 13 项性能基准，`--deny-warn` 三后端全绿）

- **混合整数二次规划（MIQP）与凸二次规划（QP）引擎** (`src/qp.mbt`, `src/miqp.mbt`, `src/benchmarks_qp.mbt`, `src/qp_wbtest.mbt`)
  - 尺度归一化 Mehrotra 预测-校正原始-对偶内点法、平移 Cholesky 半正定校验与 `NonConvex` 显式状态、无界退化射线证明、Best-Bound B&B MIQP 与固定整数连续 QP 修复启发式。
  - **13 项 v1.0.0 QP/MIQP 专项测试**：解析 2D 凸 QP、凹目标最大化与弱有效约束退化精度、无约束/纯等式快速路径、非凸/不可行/无界半正定拦截、Markowitz 有效前沿、基数约束选股 MIQP、8 步饱和线性 MPC、秩 2 半正定 SVM 对偶、单纯形最小二乘与 `.lp`/`.qps` `QUADOBJ` 无损往返、纯整数二次追踪 MIQP、电力机组组合 UC-MIQP、$L_0$ 稀疏岭回归 MIQP、推力死区混杂 MPC MIQP。
- **LP / MILP 完整求解与数值加固（67 项回归集）** (`src/simplex.mbt`, `src/dual_simplex.mbt`, `src/scaling.mbt`, `src/milp.mbt`, `src/warm_bb.mbt`, `src/cuts.mbt`, `src/branch_and_cut.mbt`, `src/presolve.mbt`)
  - 自动 Ruiz 平衡缩放、定价容差与可行性容差解耦、Dual Simplex 热启动、GMI 割平面、Presolve/Postsolve 对偶恢复，通过 Netlib `ADLITTLE`/`AFIRO`、HiGHS `AVGAS`/`SMALL_MIP`、MIPLIB `FLUGPL` 与多轮对抗数值回归。
- **独立最优性证明与双工业标准格式** (`src/certificate.mbt`, `src/lp_format.mbt`, `src/mps_format.mbt`)
  - `certify_lp_solution` & `certify_qp_solution`（线性/二次 KKT + Wolfe 强对偶）。
  - CPLEX `.lp`（含常数项与 `[ ... ] / 2` 二次目标）及 Netlib / Maros-Mészáros `.mps`/`.qps`（含 `OBJSENSE`、`RANGES`、`QUADOBJ`、`QMATRIX`）。

---

## 目录结构

```text
moonlp/
├── moon.mod                      # MoonBit 1.0 模块配置 (v1.0.0)
├── docs/
│   └── comparison_with_moonopt.md # 全生态排重对比（moonopt/moonqp）、差异矩阵与复审说明
├── .github/workflows/ci.yml      # wasm-gc / js / native 三后端 --deny-warn + bench 矩阵
├── src/
│   ├── types.mbt                 # 核心类型（含 QuadTerm、NonConvex 状态、LinearExpr）
│   ├── model.mbt                 # LP/MILP/QP/MIQP 统一代数建模 DSL 与独立验证器
│   ├── qp.mbt                    # 尺度归一化 Mehrotra 原始-对偶内点法凸 QP 求解器
│   ├── miqp.mbt                  # Best-Bound Branch-and-Bound MIQP 引擎
│   ├── benchmarks_qp.mbt         # 8 个具体 QP/MIQP 应用模型（含 4 个真实工业 MIQP）
│   ├── qp_wbtest.mbt             # 13 项 QP/MIQP 专项与格式往返测试
│   ├── certificate.mbt           # 线性 KKT + 二次 KKT / Wolfe 强对偶独立证明器
│   ├── lp_format.mbt             # CPLEX .lp 双向解析/导出（支持 [ ... ] / 2 二次目标）
│   ├── mps_format.mbt            # MPS / QPS 双向解析/导出（支持 QUADOBJ / QMATRIX）
│   ├── simplex.mbt               # Two-Phase Primal Simplex + 自动 Ruiz 缩放 + 定价解耦
│   ├── dual_simplex.mbt          # Dual Simplex 引擎
│   ├── scaling.mbt               # Ruiz 迭代 L-infinity 矩阵平衡缩放
│   ├── milp.mbt / warm_bb.mbt / cuts.mbt / branch_and_cut.mbt # 四大 MILP 引擎
│   ├── presolve.mbt / sensitivity.mbt # Presolve 对偶恢复 & 灵敏度分析
│   ├── benchmarks.mbt / benchmarks_netlib_large.mbt # Netlib/HiGHS/MIPLIB 算例集
│   ├── moonlp_wbtest.mbt / review_regression_wbtest.mbt # 67 项 LP/MILP 单元测试与回归集
│   ├── moonlp_bench_wbtest.mbt   # 13 项官方 moon bench --release 性能基准（含 4 MIQP + 3 QP）
│   └── examples/main.mbt         # 13 合 1 CLI 演示（含 4 个 MIQP 与 3 个 QP 实测输出）
└── verify/
    └── verify_qp_miqp.py         # 独立 Python/SciPy + 全组合穷举对拍脚本（7 个 QP/MIQP 模型）
```

## License
Apache-2.0
