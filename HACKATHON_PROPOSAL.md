# MoonBit 2026 开发者黑客松（11月赛季）复审项目申报书

> **仓库地址**：https://github.com/spring1219/moonlp （`v1.1.0`，GitHub ID：`spring1219`）
> **补充对比与验证文档**：[`docs/comparison_with_moonopt.md`](docs/comparison_with_moonopt.md) · 实跑互通模块：[`interop_moonopt/`](interop_moonopt/)

---

### 1. 项目名称
**`moonlp` — 面向 `moonopt` 线性模型的纯 MoonBit 混合整数锥与二次优化扩展库（QP / MIQP / QCQP / SOCP / MISOCP）**

### 2. 参赛赛道
**新项目方向 — 开源基础工具库 / 科学计算与运筹优化**

### 3. 针对初审意见的定位重构：生态边界、复用方式与净贡献

截至 2026-11-26，MoonBit 公开包生态与 GitHub 上已有数学优化项目的能力边界如下：

| 生态已有项目 | 核查版本 / 最近更新 | 核心能力 | 明确边界（无法覆盖的需求） |
| :--- | :---: | :--- | :--- |
| **`Freon793/moonopt`** | `0.6.0` (2026-11-13) | 稀疏修正单纯形 LP/MILP、分支定界、线性证书、MPS/LP | 仅支持线性目标与线性约束；**不支持二次目标（QP/MIQP）及二阶锥/二次约束（QCQP/SOCP/MISOCP）** |
| **`yujieMamba/moonqp`** | `main` (2026-10-13) | 纯连续稠密凸 QP（一阶 ADMM，线性区间约束，数组接口） | 其 README 明确排除：**整数变量（MIQP/MISOCP）**、**非线性/二阶锥约束（QCQP/SOCP）**、不可行/无界证书、`.qps`/`.lp` 文件格式与代数 DSL |
| **`Luna-Flow/linear-program`** | `main` (2026-11-12) | 轻量级 LP 求解 | 不支持 MILP / QP / MIQP / SOCP |
| **`paipai-Studio/BlackBoxOptimization`** | `main` (2026-11-12) | 黑盒元启发式 / 群智能搜索 | 无精确数学规划松弛与 KKT 强对偶证明 |

**独立建库的必要性与净贡献**：
- 我们完全认可 `Freon793/moonopt` 在纯线性稀疏 LP/MILP 上的成熟实现，不以重复实现纯线性单纯形作为本项目价值主张（大规模纯线性稀疏问题在文档中明确推荐优先使用 `moonopt`）。
- **`moonlp` v1.1.0 的实质增量工作量与不可替代贡献**在于将 MoonBit 数学规划求解内核从纯线性推进至**凸二次规划（QP）、凸二次约束规划（QCQP）、二阶锥规划（SOCP/旋转二阶锥）及其混合整数分支定界（MIQP / MISOCP）**：
  1. **Nesterov-Todd 自缩放 Mehrotra 原始-对偶二阶锥内点法**（`src/socp.mbt`）与**尺度归一化 Mehrotra 凸 QP 内点法**（`src/qp.mbt`）；
  2. **在整数分支定界节点上直接求解非线性凸 QP 与二阶锥松弛**的 `solve_miqp`（`src/miqp.mbt`）与 `solve_misocp`（`src/misocp.mbt`），配套主元 Cholesky 二次约束转锥（`add_quad_constraint`）与旋转二阶锥 DSL（`add_rotated_soc_constraint`）；
  3. **独立 Wolfe-Conic 强对偶 KKT 证明器**（`solve_socp_certified` / `certify_qp_solution`）、Conic Farkas 不可行判定、基于验证衰退射线的无界判定与非凸检测（`NonConvex`）；
  4. **工业非线性格式支持**：`.qps`（`QUADOBJ`/`QMATRIX`）与 CPLEX 二次目标/凸二次约束 `.lp`（QCP）解析与无损导出（对一般仿射右端二阶锥通过 `to_lp_string_checked` / `to_mps_string_checked` 显式报错而非静默丢约束）。

### 4. 与 `Freon793/moonopt@0.6.0` 的真实代码级复用与扩展（`interop_moonopt/`）

项目通过 MoonBit Workspace（`moon.work`）设立独立伴随模块 `interop_moonopt/`（直接依赖 `Freon793/moonopt@0.6.0` 与 `spring1219/moonlp@1.1.0`，同时保持核心库零第三方依赖），提供内存模型转换函数 `from_moonopt(@mopt_model.Model) -> @moonlp.Model` 及 MPS 无损交接。三组端到端扩展演示均已在 `wasm-gc` / `js` / `native` 三后端 CI 自动运行：

1. **主案例（决策变化最直观）— 生产规划：`moonopt` 确定性 LP → 原地扩展 95% 高斯机会约束 SOCP**
   - 在 `moonopt` 中建立确定性生产规划 LP，最优利润为 **`2333.33`**；其方案把机器期望工时用到上限，在模型的高斯不确定性假设下产能超限概率高达 **`50%`**。
   - 调用 `from_moonopt()` 无损保留原变量、目标与约束，并原地追加 95% 机会约束 Lorentz 二阶锥后，`moonlp` 求得安全最优利润 **`2046.98`**，最优产量 `(11.495, 9.475, 19.030)`，独立回算超限概率严格降至 **`5.00%`**（`4.999996%`），Conic KKT 强对偶认证通过（`gap = 1.36e-7`）。用户可直接量化为换取 95% 履约可靠性需让出的利润边际。
2. **风险与离散联合决策 — 投资组合：`moonopt` 线性基线 → 原地扩展波动率上限基数约束 MISOCP**
   - `moonopt` 线性最大化收益得 **`10.50%`**（集中持仓两类股票，年化波动率高达 `18.4%`），且实测 `@mopt_format.parse_mps` 解析 `.qps` `QUADOBJ` 直接报错拒绝（`token 'QUADOBJ'`）。
   - 调用 `from_moonopt()` 原地叠加 6 个 0-1 选股变量、最多选 3 只资产基数约束、10% 最小持仓与波动率二阶锥 $\|\Sigma^{1/2}w\|_2 \le 10.5\%$ 后，`moonlp` 求得全局最优收益 **`8.1105%`**（精选 `US_Eq 43.21% / Corp_Bond 46.79% / Gold 10.00%`，11 个 B&B 节点，`mip_gap = 0`，与 41 子集全组合穷举误差 `7.14e-11`）。
3. **MILP → MIQP 成本精确核算 — 火电机组组合：`moonopt` 线性化 MILP → 原地追加二次煤耗 MIQP**
   - 在 `moonopt` 中建立 5 机组启停与旋转备用 MILP，求解得线性化成本 **`3220.00`**（选中 3 台机组出力 `80 / 55 / 5 MW`）。
   - 调用 `from_moonopt()` 原地保留全部整数/连续变量与约束，直接追加 5 台机组真实二次热耗系数 $a_g p_g^2$ 后，`moonlp` 求解 MIQP 得精确总成本 **`3308.30`**（`mip_gap = 0`，与穷举误差 `1.55e-11`），在不改变机组启停与出力的前提下准确核算了被线性化忽略的 **`+88.30`** 二次热耗成本。

### 5. 其余新增工业基准模型（共 13 组非线性/锥/混合整数模型）

- **2D 圆推力矢量火箭制导（SOCP & 2 次点火 MISOCP）**：二次 LQR 目标自动 epigraph 提升 + 精确圆推力锥 $\sqrt{u_x^2+u_y^2} \le u_{\max} z_t$，连续 SOCP 代价 **`69.7425`**（Wolfe-Conic KKT `gap = 2.02e-6`）、限 2 次点火 MISOCP 代价 **`75.2031`**。
- **椭球不确定集鲁棒投资组合（SOCP）**（`-0.026229`，Conic KKT 认证通过）、**最小包围圆 5G 基站覆盖（6 个 Lorentz 锥 SOCP）**（`r = 2.321429`）。
- **$L_0$ 基数约束稀疏岭回归（MIQP）**（精准恢复真实支撑集 `{f0, f3}`）、**推力死区混杂 MPC（MIQP）**（`41.2281`）、**基数约束 Markowitz 资产配置（MIQP）**（`0.005202`）。

### 6. 客观工程验收与一键复现

- **自动化测试与基准**：全工作区 **95 项测试**（92 项核心测试 + 3 项 `Freon793/moonopt@0.6.0` 实跑互通测试）与 **18 项 `moon bench --release` 基准**在 `wasm-gc` / `js` / `native` 三后端 `--deny-warn` CI 全绿通过。
- **独立 Python / SciPy 数值对拍 CI**：仓库提供 `verify/verify_qp_miqp.py` 与 `verify/verify_socp_misocp.py`（使用带解析雅可比的 SciPy SLSQP + 精确组合穷举，与记录的 MoonBit 求解结果交叉核验），已集成至 `.github/workflows/ci.yml` 的 `verify-scipy` 任务持续验收。
- **一键复现命令**：
  ```bash
  moon update && moon check --deny-warn && moon test --deny-warn
  moon run interop_moonopt
  moon run src/examples
  python verify/verify_qp_miqp.py && python verify/verify_socp_misocp.py
  ```
