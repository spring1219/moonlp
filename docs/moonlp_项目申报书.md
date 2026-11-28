# `moonlp` 项目申报书（v1.1.0 复审版）

## 1. 项目名称
**`moonlp` — 面向 `moonopt` 线性模型的纯 MoonBit 混合整数锥与二次优化扩展库**

## 2. 项目简介
`moonlp`（`v1.1.0`）是一个零外部依赖、支持 `wasm-gc` / `js` / `native` 全后端的纯 MoonBit **凸二次规划（QP）、混合整数二次规划（MIQP）、凸二次约束规划（QCQP）、二阶锥规划（SOCP / 旋转二阶锥）与混合整数二阶锥规划（MISOCP）**求解与独立强对偶证明内核。

针对 MoonBit 生态已有成熟的纯线性稀疏 LP/MILP 求解器 `Freon793/moonopt@0.6.0` 与纯连续线性约束 ADMM QP 库 `yujieMamba/moonqp`，本项目明确以**向上扩展非线性锥与混合整数求解能力**为核心定位，并在工作区伴随模块 `interop_moonopt/` 中真实引入 `Freon793/moonopt@0.6.0`，提供 `from_moonopt(@mopt_model.Model) -> @moonlp.Model` 内存桥接与 MPS 无损交接，使已有 `moonopt` 线性/MILP 模型可直接原地升级为 SOCP、MIQP 与 MISOCP。

## 3. 生态边界与独立建库必要性说明
截至 2026-11-26，MoonBit 公开生态已有相关仓库的明确边界：
- **`Freon793/moonopt@0.6.0`**：强于稀疏修正单纯形 LP/MILP、分支定界与线性证明，但不支持二次目标（QP/MIQP）或二阶锥/二次约束（QCQP/SOCP/MISOCP）。
- **`yujieMamba/moonqp`**：纯连续稠密凸 QP（一阶 ADMM，仅线性区间约束），其 README 明确声明不支持整数变量（MIQP/MISOCP）、非线性/二阶锥约束（QCQP/SOCP）、不可行/无界证书与 `.qps`/`.lp` 格式。
- **`moonlp` v1.1.0 的核心净贡献**：提供零 FFI 纯 MoonBit 的 **Nesterov-Todd 原始-对偶二阶锥内点法**、**Mehrotra 凸 QP 内点法**、**基于非线性凸松弛的 MIQP / MISOCP 分支定界**以及 **Wolfe-Conic 强对偶 KKT 独立验证器**。对大规模纯线性稀疏 LP/MILP，文档明确推荐继续使用 `moonopt`。

## 4. 代表性落地场景与新增可用结果
1. **不确定需求下的机会约束生产规划（主决策案例，SOCP）**：在 `moonopt` 确定性生产规划 LP（利润 `2333.33`，在模型高斯假设下产能超限概率 `50%`）基础上，经 `from_moonopt()` 原地追加 95% 机会约束二阶锥，求得安全最优利润 **`2046.98`**、产量 `(11.495, 9.475, 19.030)`，将产能超限概率降至 **`5.00%`**（独立回算 `4.999996%`，Conic KKT 认证通过，`gap = 1.36e-7`）。
2. **风险可控的基数约束资产配置（MISOCP / MIQCQP）**：在 `moonopt` 线性组合基线（收益 `10.50%`，年化波动率 `18.4%`）基础上，经 `from_moonopt()` 原地叠加限选 3 只资产、10% 最小持仓与波动率硬上限二阶锥 $\|\Sigma^{1/2}w\|_2 \le 10.5\%$，求得全局最优收益 **`8.1105%`**（`mip_gap = 0`，与 41 子集精确穷举误差 `7.14e-11`）。
3. **电力系统机组组合二次热耗核算（MILP → MIQP 复用证据）**：在 `moonopt` 5 机组线性化 MILP（成本 `3220.00`，选中 3 机组出力 `80 / 55 / 5 MW`）基础上，经 `from_moonopt()` 原地追加二次煤耗曲线 $a_g p_g^2$，精确核算真实二次总成本 **`3308.30`**（补全 `+88.30` 二次热耗，`mip_gap = 0`）。
4. **航天/机器人圆推力矢量控制与嵌入式边缘优化**：支持 2D 圆推力约束火箭制导 SOCP（`69.7425`，Wolfe-Conic KKT `gap = 2.02e-6`）与 2 次点火开关 MISOCP（`75.2031`），纯 MoonBit 编译为 Wasm 后可在浏览器与边缘节点离线运行。

## 5. 已实现的核心技术能力（v1.1.0）
- **统一代数建模 DSL**：支持 `Minimize`/`Maximize`、连续/一般整数/0-1 变量、二次目标、标准二阶锥（`add_soc_constraint`）、旋转二阶锥（`add_rotated_soc_constraint`）与凸二次约束自动 Cholesky 转锥（`add_quad_constraint`）。
- **非线性凸与锥内点法核心**：尺度归一化 Mehrotra 凸 QP 内点法（`src/qp.mbt`）与 Nesterov-Todd 自缩放 Mehrotra 二阶锥内点法（`src/socp.mbt`），自动对二次目标+锥约束执行 epigraph 提升。
- **混合整数非线性分支定界**：Best-Bound `solve_miqp` 与 `solve_misocp`（带多变量活动界传播与整数固定节点修复）。
- **数学证明与诊断**：独立线性/二次/Wolfe-Conic 强对偶 KKT 证明器、Conic Farkas 不可行判定、基于验证衰退射线的无界判定与非凸检测（`NonConvex`）。
- **标准格式与 `moonopt` 互通**：支持 `.qps`（`QUADOBJ`/`QMATRIX`）与 CPLEX 二次目标/凸 QCP `.lp` 双向无损往返（不支持一般仿射右端锥约束时通过 `to_lp_string_checked` / `to_mps_string_checked` 显式返回 `Err` 而非静默丢失）；伴随模块 `interop_moonopt/` 实测验证 `moonopt` 复用与格式边界。

## 6. 原创与工程验证说明
- **原创项目**（Apache-2.0），核心求解算法、证明器与格式解析器均由纯 MoonBit 原创实现。
- **自动化验证**：GitHub 仓库 `https://github.com/spring1219/moonlp` 包含 **95 项工作区测试**（92 项核心 + 3 项 `Freon793/moonopt@0.6.0` 实时互通测试）、**18 项官方性能基准**及 **13 组非线性/锥/混合整数模型的独立 SciPy + 穷举验证脚本**（`verify/verify_qp_miqp.py`、`verify/verify_socp_misocp.py`），全部接入 GitHub Actions 三后端（`wasm-gc` / `js` / `native`）与 `verify-scipy` CI 持续验证。
