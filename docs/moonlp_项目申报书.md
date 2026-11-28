# `moonlp` 项目申报书

## 1. 项目名称
**moonlp**（面向 `moonopt` 线性模型的纯 MoonBit 混合整数锥与二次优化扩展库）

## 2. 项目简介
`moonlp`（`v1.1.0`）是一个零外部依赖、全后端（`wasm-gc` / `js` / `native`）通用的纯 MoonBit 数学优化求解与证明扩展库，聚焦于**凸二次规划（QP）、混合整数二次规划（MIQP）、凸二次约束规划（QCQP）、二阶锥规划（SOCP / 旋转二阶锥）与混合整数二阶锥规划（MISOCP）**。

项目在工作区伴随模块 `interop_moonopt/` 中直接集成 `Freon793/moonopt@0.6.0`，提供内存模型转换器 `from_moonopt(@mopt_model.Model) -> @moonlp.Model` 与 MPS 无损交接，使已有 `moonopt` 线性/MILP 模型可直接原地追加二次目标、二阶锥或凸二次约束并升级求解，同时配套独立 Wolfe-Conic 强对偶 KKT 证明器与 `.qps` / CPLEX QCP `.lp` 双向格式支持。

## 3. 项目方向与通用性说明
- **项目方向**：基础工具库 / 科学计算与运筹优化。
- **生态定位与差异化**：截至 2026-11-26，MoonBit 生态已有纯线性稀疏 LP/MILP 求解器 `Freon793/moonopt@0.6.0` 与纯连续线性约束 ADMM QP 库 `yujieMamba/moonqp`（其文档明确排除整数变量、非线性/二阶锥约束、不可行/无界证书与 `.qps`/`.lp` 格式）。`moonlp` 不与 `moonopt` 竞争大规模纯线性稀疏场景（文档明确推荐纯线性稀疏问题优先使用 `moonopt`），而是向上补齐**凸二次与二阶锥内点松弛、整数分支定界非线性子节点求解（MIQP/MISOCP）及 Wolfe-Conic 强对偶验证**。
- **跨平台通用性**：核心库零平台 FFI、零外部依赖，可编译至 `wasm-gc`、`js`、`native` 三后端，在浏览器、边缘设备与服务端均可离线运行。

## 4. 预期使用场景
1. **不确定需求下的机会约束生产规划（SOCP）**：在 `moonopt` 确定性生产规划 LP（利润 `2333.33`，在模型高斯假设下产能超限概率高达 `50%`）基础上，经 `from_moonopt()` 原地追加 95% 高斯机会约束二阶锥，求得安全最优利润 **`2046.98`**、产量 `(11.495, 9.475, 19.030)`，将产能超限概率严格降至 **`5.00%`**（独立回算 `4.999996%`，Conic KKT 认证通过，`gap = 1.36e-7`），帮助企业直接量化可靠性与利润的权衡代价。
2. **风险可控的基数约束资产配置（MISOCP / MIQCQP）**：在 `moonopt` 线性组合基线（最大收益 `10.50%`，年化波动率高达 `18.4%`）基础上，经 `from_moonopt()` 原地叠加限选 3 只资产、10% 最小持仓与波动率硬上限二阶锥 $\|\Sigma^{1/2}w\|_2 \le 10.5\%$，求得全局最优收益 **`8.1105%`**（精选 `US_Eq 43.21% / Corp_Bond 46.79% / Gold 10.00%`，11 节点，`mip_gap = 0`，与 41 子集穷举误差 `7.14e-11`）。
3. **电力系统机组组合二次热耗核算（MILP → MIQP）**：在 `moonopt` 5 机组线性化 MILP（成本 `3220.00`，选中 3 台机组出力 `80 / 55 / 5 MW`）基础上，经 `from_moonopt()` 保留全部变量与约束并原地追加二次煤耗曲线 $a_g p_g^2$，精确核算真实二次总成本 **`3308.30`**（补全线性化忽略的 `+88.30` 二次热耗，`mip_gap = 0`）。
4. **航天/机器人圆推力矢量控制与嵌入式边缘优化**：支持 2D 圆推力约束火箭制导 SOCP（`69.7425`，Wolfe-Conic KKT `gap = 2.02e-6`）与 2 次点火开关 MISOCP（`75.2031`）、$L_0$ 稀疏岭回归 MIQP 与推力死区混杂 MPC，编译为 Wasm 后可在前端或边缘节点本地求解。

## 5. 拟实现的核心功能
- **统一代数建模 DSL**：支持 `Minimize`/`Maximize`、连续/一般整数/0-1 变量、二次目标、标准二阶锥（`add_soc_constraint`）、旋转二阶锥（`add_rotated_soc_constraint`）与凸二次约束主元 Cholesky 自动转锥（`add_quad_constraint`）。
- **凸二次与二阶锥内点法引擎**：尺度归一化 Mehrotra 原始-对偶凸 QP 内点法（`solve_qp`）与 Nesterov-Todd 自缩放 Mehrotra 二阶锥内点法（`solve_socp`，自动对二次目标+锥约束执行 epigraph 提升）。
- **混合整数非线性分支定界引擎**：Best-Bound `solve_miqp` 与 `solve_misocp`（在整数搜索节点上求解凸 QP / SOCP 松弛，配套多变量活动界传播与整数固定节点修复）。
- **独立数学证明与数值诊断**：输出线性/二次/Wolfe-Conic 强对偶 KKT 证明（`solve_socp_certified` / `certify_qp_solution`）、Conic Farkas 不可行判定、基于验证锥衰退射线的无界判定与非凸检测（`NonConvex`）。
- **工业标准格式与 `moonopt` 互通**：支持 `.qps`（`QUADOBJ`/`QMATRIX`）与 CPLEX 二次目标/凸 QCP `.lp` 双向无损往返（对一般仿射右端二阶锥通过 `to_lp_string_checked` / `to_mps_string_checked` 显式报错而非静默丢失）；通过 `interop_moonopt/` 直接集成 `Freon793/moonopt@0.6.0` 实现内存与 MPS 互通扩展。

## 6. 原创 / 移植说明
**原创项目**（Apache-2.0），核心求解算法（Mehrotra QP IPM、Nesterov-Todd Conic IPM、MIQP/MISOCP Branch-and-Bound、Wolfe-Conic KKT 证明器）及格式解析导出器均从零以纯 MoonBit 编写，非已有代码库直接翻译移植。

## 7. 参考开源项目说明
本项目非移植项目。在生态互操作、接口风格与标准算例对拍上参考或集成了以下公开项目：
- `Freon793/moonopt@0.6.0`（https://github.com/Freon793/moonopt ）：在独立工作区伴随模块 `interop_moonopt/` 中直接作为依赖集成，用于演示线性基线复用、`from_moonopt()` 原地扩展及格式边界实测；
- `microlp`（https://github.com/Specy/microlp ，MIT 协议）与 `good_lp`（https://github.com/rust-or/good_lp ，MIT 协议）：参考其轻量化纯代码建模表达式组织风格；
- SciPy（`scipy.optimize.minimize`）：在 `verify/` 目录中用于对 13 组新增 QP / MIQP / SOCP / MISOCP 算例执行独立数值交叉验证。

## 8. GitHub 仓库链接
https://github.com/spring1219/moonlp （`v1.1.0`，全工作区 **95 项测试**（含 92 项核心测试 + 3 项 `Freon793/moonopt@0.6.0` 实时互通扩展测试）、**18 项性能基准**与 **13 组 Python/SciPy 交叉验证算例**均在 GitHub Actions `wasm-gc` / `js` / `native` 及 `verify-scipy` 四项 CI 全绿通过，对比文档见 `docs/comparison_with_moonopt.md`）
