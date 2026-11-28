# `spring1219/moonlp` v1.1.0 全生态排重对比、实质差异与新增模型验证报告

## 1. MoonBit 生态已有优化项目全面排重（截至 2026-11-26）

针对初审意见，我们对 GitHub 与 `mooncakes.io` 上的全部 MoonBit 优化相关仓库进行了全面排重，客观梳理各已有项目的能力边界与真实生态空白：

| 生态已有项目 | 最近维护 | 核心定位与算法 | 明确边界 / 不支持项 |
| :--- | :---: | :--- | :--- |
| **`Freon793/moonopt`** | 2026-11-13 | 纯线性 **LP / MILP**（稀疏修正单纯形 + 稀疏 LU + Presolve + Farkas/最优证书 + B&B + MPS/LP） | **仅支持线性目标与线性约束**；不支持 QP、MIQP、QCQP、SOCP 或 MISOCP |
| **`yujieMamba/moonqp`** | 2026-10-13 | **纯连续稠密凸 QP**（一阶 ADMM，仅线性区间约束 $l \le Ax \le u$，数组接口） | README 明文声明在支持范围之外：**不支持任何非线性/二次/二阶锥约束（QCQP/SOCP）**、**不支持任何整数变量（MIQP/MISOCP）**、**不支持不可行/无界判定与证明**、不支持 `.qps`/`.lp` 格式与代数建模 DSL |
| **`Luna-Flow/linear-program`** | 2026-11-12 | 轻量级 LP 建模与基础求解 | 无 MILP / QP / MIQP / SOCP |
| **`paipai-Studio/BlackBoxOptimization`** | 2026-11-12 | 黑盒元启发式 / 群智能优化 | 无精确数学规划、无约束全局最优保证、无 KKT/锥对偶证明 |

### 本项目 `spring1219/moonlp` v1.1.0 的定位与生态净贡献

> **面向 `moonopt` 线性模型的纯 MoonBit 混合整数锥与二次优化扩展库：在保持零外部依赖与 `wasm-gc`/`js`/`native` 三后端一致性的前提下，将 MoonBit 数学规划求解能力从纯线性（`Freon793/moonopt@0.6.0`）与纯连续线性约束 QP（`yujieMamba/moonqp`）扩展至凸二次约束规划（QCQP）、二阶锥规划（SOCP / 旋转二阶锥）及其混合整数分支定界（MIQP / MISOCP）与 Wolfe-Conic 强对偶证明。**

---

## 2. 核心能力矩阵对比（`moonopt@0.6.0` vs `moonqp` vs `moonlp` v1.1.0）

| 能力维度 | `Freon793/moonopt@0.6.0` | `yujieMamba/moonqp` | **`spring1219/moonlp` v1.1.0** |
| :--- | :---: | :---: | :---: |
| **二阶锥规划（SOCP）** $\|A_k x + b_k\|_2 \le c_k^T x + d_k$ | 不支持 | 不支持 | **原生支持**（Nesterov-Todd 自缩放 + Mehrotra 预测-校正原始-对偶锥内点法 `solve_socp`） |
| **凸二次约束规划（QCQP）** $\frac{1}{2}x^T P_k x + q_k^T x \le r_k$ | 不支持 | 不支持 | **原生支持**（`add_quad_constraint` 主元 Cholesky 半正定分解自动无损提升为 Lorentz 锥） |
| **混合整数二阶锥规划（MISOCP / MIQCQP）** | 不支持 | **明确不支持** | **原生支持**（`solve_misocp`：Best-Bound B&B + 节点 SOCP 松弛 + 活动界传播 + 固定整数修复） |
| **混合整数二次规划（MIQP）** | 不支持 | **明确不支持** | **原生支持**（`solve_miqp`：Best-Bound B&B + 节点凸 QP 松弛 + 固定整数 QP 修复） |
| **锥不可行 / 验证射线无界 / 非凸诊断** | 仅 LP | **不支持** | **原生支持**（Conic Farkas `Infeasible`、仅在确认锥衰退射线时返回 `Unbounded`、非凸返回 `NonConvex`） |
| **独立锥 KKT 与 Wolfe 强对偶证明** | 仅线性 | 无强对偶证明 | **原生支持** `solve_socp_certified`（验证对偶 Lorentz 锥归属 $z_{k,0}\ge\|\bar{z}_k\|_2$ 与强对偶间隙）+ `certify_qp_solution` |
| **`.qps` (`QUADOBJ`/`QMATRIX`) & CPLEX 二次/QCP `.lp`** | 仅线性 MPS/LP | **不支持** | **支持 LP/MILP/QP/MIQP 及凸 QCP/常右端球约束双向无损往返**；一般仿射右端二阶锥经 `to_lp_string_checked`/`to_mps_string_checked` 显式报错拒绝 |
| **大规模稀疏纯线性 LP/MILP** | **强项（推荐）** | 不支持 | 兼容支持（文档明确推荐大规模纯线性稀疏问题优先使用 `moonopt`） |

---

## 3. 与 `Freon793/moonopt@0.6.0` 的真实代码级复用与原地扩展实测（`interop_moonopt/`）

除算法差异外，仓库通过 [`moon.work`](../moon.work) 在独立互通模块 [`interop_moonopt/`](../interop_moonopt/) 中**直接引入了 `Freon793/moonopt@0.6.0`**，实现 `from_moonopt(@mopt_model.Model) -> @moonlp.Model` 内存桥接与 MPS 无损交接。三组端到端实验均在 `from_moonopt()` 转换所得的同一模型上原地追加新约束/二次项，并在 CI 三后端自动运行验证：

1. **生产规划（主决策案例：复用 `moonopt` LP → 原地扩展 95% 机会约束 SOCP）**：`moonopt` 求得确定性 LP = `2333.33`（MPS 交接同解 `2333.33`，但在模型高斯不确定性假设下产能超限概率达 `50%`）；经 `from_moonopt()` 原地追加 95% 高斯机会约束 Lorentz 锥后得到安全最优 **`2046.98`**（产量 `(11.495, 9.475, 19.030)`，独立回算超限概率严格降至 **`5.00%`**，Conic KKT 认证通过，`gap = 1.36e-7`）。
2. **投资组合（实证 `.qps`/QCP 边界 → 原地扩展波动率上限 MISOCP）**：`moonopt` 线性基线得 `10.50%`（波动率 `18.4%`），且实测 `@mopt_format.parse_mps` 对 `.qps` `QUADOBJ` 直接报错拒绝（`token 'QUADOBJ'`）、`@mopt_format.read_lp` 无法求解 CPLEX QCP `[ x^2 + y^2 ] <= 25`；经 `from_moonopt()` 原地追加波动率上限 $\le 10.5\%$ 与限选 3 资产约束后求得全局最优 **`8.1105%`**（11 节点，`mip_gap = 0`）。
3. **电力机组组合（MILP → MIQP 成本精确核算）**：`moonopt` 线性化 MILP 得 `3220.00`（选中 3 机组出力 `80 / 55 / 5 MW`）；经 `from_moonopt()` 原地导入并叠加 5 台机组真实二次热耗系数 $a_g p_g^2$，在 `moonlp` 精确核算真实二次总成本 **`3308.30`**（补全 `+88.30` 二次热耗，`mip_gap = 0`）。

---

## 4. 具体优化模型与新增可用结果（13 个非线性/锥/混合整数模型，全部经 SciPy + 穷举独立对拍）

以下 **10 个核心代表模型**均无法被 `moonopt` 或 `moonqp` 直接建模求解（详见 [`src/benchmarks_socp.mbt`](../src/benchmarks_socp.mbt) 与 [`src/benchmarks_qp.mbt`](../src/benchmarks_qp.mbt)），可通过 `moon run src/examples`、`python verify/verify_socp_misocp.py` 与 `python verify/verify_qp_miqp.py` 一键复现。

### A. 二阶锥规划（SOCP / QCQP）与混合整数二阶锥规划（MISOCP）—— 全生态 0 覆盖

| # | 具体模型 | 锥/离散结构特征 | `moonlp` 实测最优结果 | 独立证明 / SciPy+穷举误差 | `wasm-gc` 耗时 |
| :---: | :--- | :--- | :--- | :---: | ---: |
| **1** | **波动率硬上限基数约束投资组合（MISOCP / MIQCQP）** (`build_vol_capped_cardinality_misocp`, 6 选 $\le 3$ 资产, $\sigma \le 10.5\%$) | **非线性波动率锥约束** $\|\Sigma^{1/2}w\|_2 \le 0.105$ + 6 二值选股变量 + 半连续持仓 | 最大收益 **`8.11046%`**（11 节点，`mip_gap=7.83e-10`；精选 **US_Eq 43.21% / Corp_Bond 46.79% / Gold 10.00%**） | vs 41 子集穷举精确最优 **`7.14e-11`** | **`15.03 ms`** |
| **2** | **2D 圆推力矢量火箭制导与开关控制（SOCP & MISOCP）** (`build_isotropic_thrust_mpc_socp` / `build_switched_isotropic_mpc_misocp`) | 二次 LQR 目标（自动 epigraph 提升）+ **精确圆推力锥** $\sqrt{u_{x,t}^2+u_{y,t}^2}\le u_{\max} z_t$ + $\le 2$ 次点火 | 4 步连续 SOCP 代价 **`69.742457`**；限 2 次点火 MISOCP 代价 **`75.203126`**（点火序列 `[1,1,0,0]`） | vs SciPy+全序列穷举 **`3.29e-8` / `9.10e-7`** | **`40.19 ms`** |
| **3** | **椭球不确定集鲁棒投资组合（Robust SOCP）** (`build_robust_portfolio_socp`, $\kappa=1.5$) | 精确最大化置信椭球最坏情形收益 $\max \mu^T w - \kappa \|\Sigma^{1/2} w\|_2$ | 最坏情形收益 **`-0.02622928`**（风险 `5.350%`） | Conic KKT=`true`；vs SciPy **`1.16e-12`** | **`165.27 µs`** |
| **4** | **95% 高斯机会约束生产规划（Chance-Constrained SOCP）** (`build_chance_constrained_socp`) | 概率约束 $\Pr(\tilde{a}^T x \le 100)\ge 95\%$ 精确等价二阶锥 | 安全最大利润 **`2046.9777`**，$x=(11.495, 9.475, 19.030)$ | Conic KKT=`true`；vs SciPy **`7.93e-8`** | **`62.31 µs`** |
| **5** | **最小包围圆 / 5G 基站最小覆盖半径（SOCP）** (`build_smallest_enclosing_disk_socp`) | 6 个 Lorentz 锥 $\|x - p_i\|_2 \le r$ | 最小半径 **`2.32142858`**，圆心 `(2.000, 1.1786)` | Conic KKT=`true` (gap `1.52e-8`)；vs 解析解 **`7.21e-9`** | **`2.00 ms`** |

### B. 混合整数二次规划（MIQP）—— `moonopt` 与 `moonqp` 均无法求解

| # | 具体 MIQP 模型 | `moonlp` 全局最优结果 | 穷举+SciPy 误差 | `wasm-gc` 耗时 |
| :---: | :--- | :--- | :---: | ---: |
| **6** | **电力系统机组组合与二次经济调度（UC-MIQP）**（5 机组启停 + 140 MW + 20 MW 备用） | 总成本 **`3308.30`**（`Coal 80MW + CCGT 55MW + Peaker2 5MW`） | **`1.55e-11`** | **`1.78 ms`** |
| **7** | **基数约束 Markowitz 资产配置（MIQP）**（6 选 $\le 3$ 资产 + 最小持仓 10%） | 最小方差 **`0.00520192`**（`US_Eq 41% / Corp_Bond 49% / Gold 10%`） | **`4.37e-14`** | **`5.17 ms`** |
| **8** | **$L_0$ 基数约束稀疏岭回归（MIQP）**（5 选 $\le 2$ 特征） | 损失 **`0.23809507`**，精准恢复真实支撑集 **`{f0, f3}`** | **`1.59e-13`** | **`2.11 ms`** |
| **9** | **推力最小冲量死区混杂 MPC（MIQP-MPC）**（$u_t \in \{0\}\cup[0.3,1.0]$，$\le 3$ 次点火） | 代价 **`41.228125`**，开关 **`[1.0, 1.0, OFF, OFF, OFF]`** | **`1.46e-11`** | **`619.16 µs`** |
| **10** | **Markowitz 有效前沿 QP / 8 步饱和 MPC QP / 秩 2 SVM 对偶 QP + `.qps`/`.lp` 往返** | Wolfe KKT 强对偶认证全通过，误差 $\le 1.44\times 10^{-11}$ | **`5.76e-13..1.44e-11`** | **`59.9 µs..996 µs`** |

---

## 5. 复审申报表可直接引用的补充说明（精简版）

> 针对初审意见，我们对截至 2026-11-26 的 MoonBit 生态全部已有优化库（纯线性 LP/MILP 求解器 `Freon793/moonopt@0.6.0` 与纯连续线性约束 ADMM QP 库 `yujieMamba/moonqp`）完成了边界核查，并将项目定位明确升级为 **`v1.1.0` — 面向 `moonopt` 线性模型的纯 MoonBit 混合整数锥与二次优化扩展库**：
>
> 1. **与 `moonopt` 的真实代码级复用与扩展，以及相对已有生态的新增内核**：
>    - **直接复用与原地扩展 `moonopt`**：通过 `moon.work` 在伴随模块 `interop_moonopt/` 中直接依赖 `Freon793/moonopt@0.6.0`，实现 `from_moonopt(@mopt_model.Model) -> @moonlp.Model` 与 MPS 无损交接——三组端到端案例（生产规划、投资组合、火电机组组合）均先由 `moonopt` 求解线性基线，再直接在 `from_moonopt()` 转换所得的同一模型上原地追加二阶锥/二次目标/离散约束，并实测证明 `moonopt` 会拒绝 `.qps (QUADOBJ)` 且无法求解 CPLEX QCP `.lp`。
>    - **相对 `moonopt` 与 `moonqp` 的核心新增能力**：`moonopt` 专注纯线性 LP/MILP；`yujieMamba/moonqp` 仅支持连续线性约束 ADMM QP（其 README 明确排除整数变量、非线性/二阶锥约束、不可行/无界证书与 `.qps`/`.lp` 格式）。本项目零 FFI 三后端新增实现：**Nesterov-Todd 自缩放 Mehrotra 原始-对偶二阶锥内点法（`solve_socp`）**、**旋转二阶锥（`add_rotated_soc_constraint`）与主元 Cholesky 二次约束自动提升（`add_quad_constraint`）**、**在整数分支定界节点求解非线性凸松弛的 MISOCP / MIQP 引擎（`solve_misocp` / `solve_miqp`）**、**Wolfe-Conic 强对偶 KKT 独立证明器（`solve_socp_certified` / `certify_qp_solution`）**，以及 **`.qps (QUADOBJ)` 与 CPLEX 二次目标/凸 QCP `.lp` 双向无损往返**（对不支持导出的仿射右端一般二阶锥显式返回 `Err`）。
> 2. **具体优化模型与新增可用结果（13 个模型均经 SciPy + 全组合穷举脚本交叉核验并纳入 CI）**：
>    - **主决策案例 — 95% 高斯机会约束生产规划（SOCP）**：`moonopt` 确定性 LP 利润 `2333.33` 在模型高斯假设下产能超限概率高达 `50%`；经 `from_moonopt()` 原地追加 95% 机会约束二阶锥后求得安全最大利润 **`2046.98`**、产量 `(11.495, 9.475, 19.030)`，独立回算超限概率严格降至 **`5.00%`**（`4.999996%`，Conic KKT 认证通过，`gap=1.36e-7`）；
>    - **风险与离散联合案例 — 波动率硬上限基数约束投资组合（MISOCP）**：`moonopt` 线性基线收益 `10.50%` 但波动率高达 `18.4%`；经 `from_moonopt()` 原地施加 $\|\Sigma^{1/2}w\|_2 \le 10.5\%$ 锥约束与限选 3 只资产约束求得全局最大收益 **`8.11046%`**（`US_Eq 43.21% / Corp_Bond 46.79% / Gold 10.00%`，11 节点，`mip_gap=0`，与 41 子集穷举误差 `7.14e-11`）；
>    - **代码复用与成本核算案例 — 电力机组组合（MILP → MIQP）**：经 `from_moonopt()` 保留全部机组变量与约束并原地追加二次热耗曲线 $a_g p_g^2$，在同样选中 `80 / 55 / 5 MW` 机组组合下，精确核算真实二次总成本 **`3308.30`**（补全线性化 MILP `3220.00` 忽略的 `+88.30` 二次热耗，`mip_gap=0`，误差 `1.55e-11`）；
>    - **其余工业模型**：2D 圆推力矢量火箭制导 SOCP（`69.7425`，Wolfe-Conic KKT `gap=2.02e-6`）与 2 次点火 MISOCP（`75.2031`）、椭球鲁棒投资组合 SOCP（`-0.026229`）、最小包围圆 5G 基站覆盖 SOCP（`2.321429`）、$L_0$ 稀疏岭回归 MIQP（精准恢复 `{f0, f3}`）与推力死区混杂 MPC MIQP（`41.2281`）。
> 3. **工程验收**：全部 **95 项工作区测试**（92 项核心 + 3 项 `Freon793/moonopt@0.6.0` 真实互通扩展测试）、**18 项 `moon bench --release` 基准**以及 **Python SciPy 对拍脚本**（`verify/verify_qp_miqp.py`、`verify/verify_socp_misocp.py`，将独立 SciPy/穷举解与记录的 MoonBit 结果交叉验证）均在 GitHub Actions `wasm-gc` / `js` / `native` + `verify-scipy` CI 全绿通过。
