# `moonlp` — Pure MoonBit Convex QP / MIQP & LP / MILP Optimization Solver

[![MoonBit CI](https://github.com/spring1219/moonlp/actions/workflows/ci.yml/badge.svg)](https://github.com/spring1219/moonlp/actions/workflows/ci.yml)
![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)
![Backends: wasm-gc | js | native](https://img.shields.io/badge/targets-wasm--gc%20%7C%20js%20%7C%20native-brightgreen.svg)

> **零外部依赖、全后端（`wasm-gc` / `js` / `native`）通用的纯 MoonBit 凸二次规划（QP）、混合整数二次规划（MIQP）及线性/混合整数线性规划（LP/MILP）求解与证明库**  
> 在兼容标准 LP/MILP 与 `.mps`/`.lp` 的基础上，补齐 MoonBit 生态在 **凸二次规划（Primal-Dual 内点法）**、**混合整数二次规划（MIQP B&B）**、**二次 KKT / Wolfe 强对偶认证** 与 **`.qps` (`QUADOBJ`) / CPLEX 二次 `.lp`** 上的能力空白。

---

## 与生态已有 `Freon793/moonopt` 的定位分工与实质差异

MoonBit 生态中已有优秀的纯 MoonBit 稀疏修正单纯形 LP/MILP 求解器 [`Freon793/moonopt`](https://github.com/Freon793/moonopt)。`moonlp` v1.0.0 与其形成清晰互补（详细对比见 [docs/comparison_with_moonopt.md](docs/comparison_with_moonopt.md)）：

- **纯大规模稀疏 LP/MILP**：推荐使用 `Freon793/moonopt`（稀疏修正单纯形 + 稀疏 LU 分解）；两库通过标准 `.mps` / `.lp` 无损互通。
- **凸二次规划（QP）与混合整数二次规划（MIQP）**：`moonopt` 仅支持线性目标 $c^T x$，无法表达或求解 $\frac{1}{2}x^T Q x + c^T x$。`moonlp` v1.0.0 原生提供：
  1. **尺度归一化 Mehrotra 预测-校正原始-对偶内点法（`solve_qp`）**：支持任意正半定 Hessian $Q \succeq 0$（包括秩亏 SVM Gram 矩阵与含二值指示变量的 MIQP 松弛）、任意等式/不等式/盒约束/自由变量，显式拦截非凸输入（`NonConvex`）与无界退化射线（`Unbounded`）；
  2. **混合整数二次规划引擎（`solve_miqp`）**：Best-Bound Branch-and-Bound + 节点凸 QP 松弛 + 固定整数连续 QP 修复启发式 + 严格状态与 `mip_gap` 传播；
  3. **独立二次 KKT 与 Wolfe 强对偶证明器（`certify_qp_solution`）**：独立检验 $Q \succeq 0$、原始可行性、对偶符号、互补松弛、梯度驻点 $Qx + c - A^T y - z_l + z_u = 0$ 与 Wolfe 强对偶间隙；
  4. **工业二次格式无损往返**：支持 `.qps`（MPS `QUADOBJ` / `QMATRIX` 段）与 CPLEX 二次 `.lp`（`+ [ ... ] / 2`）。

---

## 核心新增具体模型与客观验证结果（`moonopt` 无法直接求解）

以下具体模型全部内置于 [`src/benchmarks_qp.mbt`](src/benchmarks_qp.mbt)，经 [`verify/verify_qp_miqp.py`](verify/verify_qp_miqp.py)（SciPy SLSQP + MIQP 全子集穷举）独立对拍验证，误差均 $\le 1.44\times 10^{-11}$：

| 具体优化模型 | 问题类别与特征 | `moonlp` 实测最优结果 | 独立验证 / SciPy 对拍误差 | `wasm-gc` 实测耗时 |
| :--- | :--- | :--- | :--- | ---: |
| **Markowitz 6 资产投资组合** (`build_markowitz_portfolio_qp`, $r \ge 8\%$) | 连续凸 QP（$6\times 6$ 协方差矩阵 $\Sigma$ + 有效前沿） | $\frac{1}{2}w^T\Sigma w = \mathbf{0.004958677971}$（8 次迭代，5 资产分散配置） | Wolfe KKT=`true` (gap `5.89e-11`)；vs SciPy **`1.47e-12`** | **`57.35 µs`** |
| **基数约束 Markowitz 选股** (`build_cardinality_portfolio_miqp`, 6 选 $\le 3$，最小持仓 10%) | **混合整数二次规划 MIQP**（12 变量、秩 6 半正定 Hessian + 6 二值变量） | 全局最优 **`0.005201920000`**（9 个 B&B 节点，`mip_gap=2.16e-12`；精选 **US_Eq 41% + Corp_Bond 49% + Gold 10%**） | vs 全 41 子集穷举精确解误差 **`4.37e-14`** | **`4.40 ms`** |
| **带饱和约束线性 MPC 轨迹跟踪** (`build_linear_mpc_qp`, $N=8$) | 连续凸 QP（24 变量 + 16 动力学等式 + $\|u_t\| \le 1.0$, $\|v_t\| \le 1.5$） | LQR 代价 **`38.7957968155`**；$u_0=-1.0$ 精准限幅，$p_0=2.0 \to p_8=-0.011$ | Wolfe KKT=`true`；vs SciPy **`5.76e-13`** | **`857.16 µs`**（>1,100 Hz） |
| **软间隔 SVM 对偶训练** (`build_svm_dual_qp`, $C=10$) | **秩亏半正定 QP**（$8\times 8$ **秩 2** Gram 矩阵 + 盒约束） | 对偶最优 **`-0.0928936368`**（7 次迭代，2 支持向量，全样本间隔 $\ge 1.0$） | Wolfe KKT=`true`；vs SciPy **`1.44e-11`** | **`74.12 µs`** |
| **单纯形约束最小二乘 / 光谱解混** (`build_constrained_lsq_qp`) | 连续凸 QP + `.lp` 二次 / `.qps` `QUADOBJ` 格式往返 | 原模型、`.lp` 往返、`.qps` 往返最优值与 KKT 证书完全一致 | Wolfe KKT=`true`（三路径全绿） | — |

---

## 快速开始

```bash
moon check --deny-warn
moon test --deny-warn
moon run src/examples
moon bench --release

# 独立 Python / SciPy + 穷举对拍验证 QP & MIQP 模型
python verify/verify_qp_miqp.py

# 跨后端全量验证（与 GitHub Actions CI 一致）
moon test --deny-warn --target wasm-gc && moon bench --release --target wasm-gc
moon test --deny-warn --target js      && moon bench --release --target js
moon test --deny-warn --target native  && moon bench --release --target native
```

### 示例 1：基数约束 Markowitz 投资组合优化（MIQP）与 Wolfe KKT 认证（QP）

```moonbit
// 1. 连续 Markowitz 均值-方差投资组合（Convex QP）+ Wolfe 强对偶 KKT 证明
let qp = @moonlp.build_markowitz_portfolio_qp(target_return=0.08)
let qsol = qp.solve()                               // 0.00495868 (8 IPM iters, ~57 us)
let qcert = qp.certify_qp_solution(qsol)            // certified_optimal == true

// 2. 基数约束投资组合（MIQP：6 资产限选 <= 3 只、单只持仓 10%~50%）
let miqp = @moonlp.build_cardinality_portfolio_miqp(
  max_assets=3, min_weight=0.1, max_weight=0.5, target_return=0.08
)
let msol = miqp.solve()                             // 0.00520192 (US_Eq 41%, Corp_Bond 49%, Gold 10%)
```

### 示例 2：带执行器饱和约束的实时线性 MPC（Convex QP）

```moonbit
let mpc = @moonlp.build_linear_mpc_qp(horizon=8, x0_pos=2.0, x0_vel=0.0, u_max=1.0)
let sol = mpc.solve()                               // LQR = 38.7958, u_0 = -1.0, ~857 us on wasm-gc
let cert = mpc.certify_qp_solution(sol)             // certified_optimal == true
```

---

## 全量实测性能基准（`moon bench --release --target wasm-gc`）

以下为 Windows 本机实测耗时（`moon bench --release --target wasm-gc`，10 轮统计均值 ± 标准差）：

| 基准算例 | 规模 / 类型 | 实测均值 (`mean ± σ`) | 区间 (`min … max`) |
| :--- | :--- | ---: | ---: |
| **Markowitz 6 资产投资组合（QP）** | 6 变量正定协方差 QP（IPM） | **`57.35 µs`** ± 4.76 µs | `54.42 µs … 67.02 µs` |
| **软间隔 SVM 对偶训练（QP）** | 8×8 秩 2 半正定 Gram QP（IPM） | **`74.12 µs`** ± 2.27 µs | `72.60 µs … 79.90 µs` |
| **8 步约束线性 MPC（QP）** | 24 变量 + 16 动力学等式 QP（IPM） | **`857.16 µs`** ± 6.55 µs | `845.24 µs … 868.28 µs` |
| **基数约束 Markowitz 选股（MIQP）** | 6 连续 + 6 二值、秩 6 半正定 MIQP | **`4.40 ms`** ± 65.32 µs | `4.34 ms … 4.52 ms` |
| **HiGHS `AVGAS`（LP）** | 11×8 LP（38 非零元） | **`5.00 µs`** ± 119.00 ns | `4.93 µs … 5.33 µs` |
| **Netlib `AFIRO`（LP）** | 28×32 LP（88 非零元） | **`15.89 µs`** ± 653.70 ns | `15.34 µs … 17.42 µs` |
| **Netlib `ADLITTLE`（LP）** | 56×97 LP（480 非零元） | **`623.70 µs`** ± 3.89 µs | `619.11 µs … 632.12 µs` |
| **5 城 TSP-MTZ（Branch-and-Cut MILP）** | 20 二值 + 4 连续 MILP | **`67.42 µs`** ± 1.92 µs | `64.79 µs … 70.17 µs` |
| **10 物品 0-1 背包（B&B MILP）** | 10 二值 MILP | **`70.53 µs`** ± 277.37 ns | `70.28 µs … 71.21 µs` |

---

## 核心特性与测试覆盖（77 项单元测试 + 10 项性能基准，`--deny-warn` 三后端全绿）

- **凸二次规划（QP）与混合整数二次规划（MIQP）引擎** (`src/qp.mbt`, `src/miqp.mbt`, `src/benchmarks_qp.mbt`, `src/qp_wbtest.mbt`)
  - 尺度归一化 Mehrotra 预测-校正原始-对偶内点法（避免 $s_i, z_i \to 0$ 除零精度衰减）、平移 Cholesky 半正定校验与 `NonConvex` 显式状态、无界退化射线验证、固定整数连续 QP 修复启发式。
  - **10 项 v1.0.0 QP/MIQP 专项测试**：解析 2D 凸 QP、凹目标最大化与弱有效约束退化精度、无约束/纯等式快速路径、非凸/不可行/无界半正定拦截、Markowitz 有效前沿单调性、基数约束 MIQP 全局最优、8 步约束线性 MPC、秩 2 半正定 SVM 对偶训练、单纯形最小二乘与 `.lp`/`.qps` `QUADOBJ` 无损往返、纯整数二次追踪 MIQP。
- **LP / MILP 完整求解与数值加固（63 项回归集）** (`src/simplex.mbt`, `src/dual_simplex.mbt`, `src/scaling.mbt`, `src/milp.mbt`, `src/warm_bb.mbt`, `src/cuts.mbt`, `src/branch_and_cut.mbt`, `src/presolve.mbt`)
  - 自动 Ruiz 平衡缩放、定价容差与可行性容差解耦、Dual Simplex 热启动、GMI 割平面、Presolve/Postsolve 对偶恢复，通过 Netlib `ADLITTLE`/`AFIRO`、HiGHS `AVGAS`/`SMALL_MIP`、MIPLIB `FLUGPL` 与 4 轮对抗数值回归。
- **独立最优性证明与双工业标准格式** (`src/certificate.mbt`, `src/lp_format.mbt`, `src/mps_format.mbt`)
  - `certify_lp_solution` & `certify_qp_solution`（线性/二次 KKT + Wolfe 强对偶）。
  - CPLEX `.lp`（含常数项与 `[ ... ] / 2` 二次目标）及 Netlib / Maros-Mészáros `.mps`/`.qps`（含 `OBJSENSE`、`RANGES`、`QUADOBJ`、`QMATRIX`）。

---

## 目录结构

```text
moonlp/
├── moon.mod                      # MoonBit 1.0 模块配置 (v1.0.0)
├── docs/
│   └── comparison_with_moonopt.md # 与 Freon793/moonopt 的定位分工、差异对照与实测验证
├── .github/workflows/ci.yml      # wasm-gc / js / native 三后端 --deny-warn + bench 矩阵
├── src/
│   ├── types.mbt                 # 核心类型（含 QuadTerm、NonConvex 状态、LinearExpr）
│   ├── model.mbt                 # LP/MILP/QP/MIQP 统一代数建模 DSL 与独立验证器
│   ├── qp.mbt                    # 尺度归一化 Mehrotra 原始-对偶内点法凸 QP 求解器
│   ├── miqp.mbt                  # Best-Bound Branch-and-Bound MIQP 引擎
│   ├── benchmarks_qp.mbt         # 5 个具体 QP/MIQP 应用模型（Markowitz/MIQP/MPC/SVM/LSQ）
│   ├── qp_wbtest.mbt             # 10 项 QP/MIQP 专项与格式往返测试
│   ├── certificate.mbt           # 线性 KKT + 二次 KKT / Wolfe 强对偶独立证明器
│   ├── lp_format.mbt             # CPLEX .lp 双向解析/导出（支持 [ ... ] / 2 二次目标）
│   ├── mps_format.mbt            # MPS / QPS 双向解析/导出（支持 QUADOBJ / QMATRIX）
│   ├── simplex.mbt               # Two-Phase Primal Simplex + 自动 Ruiz 缩放 + 定价解耦
│   ├── dual_simplex.mbt          # Dual Simplex 引擎
│   ├── scaling.mbt               # Ruiz 迭代 L-infinity 矩阵平衡缩放
│   ├── milp.mbt / warm_bb.mbt / cuts.mbt / branch_and_cut.mbt # 四大 MILP 引擎
│   ├── presolve.mbt / sensitivity.mbt # Presolve 对偶恢复 & 灵敏度分析
│   ├── benchmarks.mbt / benchmarks_netlib_large.mbt # Netlib/HiGHS/MIPLIB 算例集
│   ├── moonlp_wbtest.mbt         # 63 项 LP/MILP 单元测试与 4 轮对抗回归集
│   ├── moonlp_bench_wbtest.mbt   # 10 项官方 moon bench --release 性能基准（含 QP/MIQP）
│   └── examples/main.mbt         # 10 合 1 CLI 演示（含 Markowitz QP/MIQP、MPC、SVM）
└── verify/
    └── verify_qp_miqp.py         # 独立 Python/SciPy + 全子集穷举对拍验证脚本
```

## License
Apache-2.0
