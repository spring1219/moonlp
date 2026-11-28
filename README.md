# `moonlp` v1.1.0 — 面向 `moonopt` 线性模型的纯 MoonBit 混合整数锥与二次优化扩展库

[![MoonBit CI](https://github.com/spring1219/moonlp/actions/workflows/ci.yml/badge.svg)](https://github.com/spring1219/moonlp/actions/workflows/ci.yml)
![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)
![Backends: wasm-gc | js | native](https://img.shields.io/badge/targets-wasm--gc%20%7C%20js%20%7C%20native-brightgreen.svg)

> **零外部依赖、全后端（`wasm-gc` / `js` / `native`）通用的纯 MoonBit 二阶锥规划（SOCP / 旋转二阶锥）、混合整数二阶锥规划（MISOCP）、凸二次约束规划（QCQP）与混合整数二次规划（MIQP）求解及强对偶证明扩展库**（可复用 `Freon793/moonopt@0.6.0` 线性模型并原地扩展，同时兼容 `.mps` / `.qps` / `.lp` 工业格式）。

---

## 与 `Freon793/moonopt` 的真实代码级互通与原地扩展（`interop_moonopt/`）

本项目通过 MoonBit Workspace（[`moon.work`](moon.work)）在独立伴随模块 [`interop_moonopt/`](interop_moonopt/) 中**直接引入 `Freon793/moonopt@0.6.0`**（核心库 `spring1219/moonlp` 本身仍保持零外部依赖），提供内存转换器 `from_moonopt(@mopt_model.Model) -> @moonlp.Model` 与 MPS 无损交接。以下三组端到端实验均直接对 `from_moonopt()` 转换所得的同一模型原地追加非线性/锥/离散约束，并随 CI 三后端自动执行 `moon run interop_moonopt`：

| 实测场景 | `Freon793/moonopt@0.6.0` 线性基线实测 | `spring1219/moonlp` v1.1.0 原地扩展后新增可用结果 |
| :--- | :--- | :--- |
| **1. 生产规划 → 95% 高斯机会约束 SOCP（主决策案例）** | 确定性 LP = `2333.33`（MPS 交接同解 `2333.33`；在模型高斯假设下产能超限概率 **`50%`**） | 经 `from_moonopt()` 原地追加 95% 机会约束 Lorentz 锥 → 安全最优 **`2046.98`**，超限概率严格降至 **`5.00%`**（Conic KKT=`true`，`gap=1.36e-7`） |
| **2. 投资组合 → 波动率上限基数约束 MISOCP** | 线性最大收益 = `10.50%`（波动率 `18.4%`；`parse_mps(.qps QUADOBJ)` 直接报错拒绝） | 经 `from_moonopt()` 原地追加波动率锥 $\|\Sigma^{1/2}w\|_2\le 10.5\%$ 与限选 3 资产约束 → 全局最优 **`8.1105%`**（11 节点，`mip_gap=0`） |
| **3. 电力机组组合 → 二次煤耗成本精确核算 MIQP** | 线性化 MILP 成本 = `3220.00`（选中 3 机组出力 `80 / 55 / 5 MW`，忽略 $a_g p_g^2$） | 经 `from_moonopt()` 原地追加二次热耗曲线 → 精确核算总成本 **`3308.30`**（补全 `+88.30` 二次热耗，`mip_gap=0`） |

---

## 与 MoonBit 生态已有优化库（截至 2026-11-26）的边界对比

详细全生态排重对比见 [docs/comparison_with_moonopt.md](docs/comparison_with_moonopt.md)：

| 能力维度 | `Freon793/moonopt@0.6.0`（线性稀疏 LP/MILP） | `yujieMamba/moonqp`（连续 ADMM QP） | **`spring1219/moonlp` v1.1.0（本项目）** |
| :--- | :---: | :---: | :---: |
| **二阶锥规划（SOCP）** $\|A_k x+b_k\|_2 \le c_k^T x+d_k$ | 不支持 | 不支持 | **原生支持**（Nesterov-Todd 自缩放 Mehrotra 原始-对偶锥内点法 `solve_socp`） |
| **凸二次约束规划（QCQP）** $\frac{1}{2}x^T P x + q^T x \le r$ | 不支持 | 不支持 | **原生支持**（`add_quad_constraint` 主元 Cholesky 自动无损二阶锥提升） |
| **混合整数二阶锥（MISOCP / MIQCQP）** | 不支持 | **明确不支持** | **原生支持**（`solve_misocp`：B&B + 节点 SOCP 松弛 + 活动界传播） |
| **混合整数二次规划（MIQP）** | 不支持 | **明确不支持** | **原生支持**（`solve_miqp`：B&B + 节点凸 QP 松弛 + 固定整数修复） |
| **锥不可行 / 验证射线无界 / 非凸诊断** | 仅 LP | **不支持** | **原生支持**（Conic Farkas `Infeasible` / 验证锥衰退射线 `Unbounded` / 非凸 `NonConvex`；证据不足时不误判无界） |
| **独立锥 KKT & Wolfe 强对偶证明** | 仅线性 | 无强对偶证明 | **原生支持** `solve_socp_certified` & `certify_qp_solution` |
| **`.qps` (`QUADOBJ`) & CPLEX 二次/QCP `.lp`** | 仅线性 MPS/LP | **不支持** | **支持 LP/MILP/QP/MIQP 及凸 QCP/常右端球约束双向无损往返**；一般仿射右端二阶锥经 `to_lp_string_checked` / `to_mps_string_checked` 显式报错而非静默丢失 |
| **大规模稀疏纯线性 LP/MILP** | **强项（推荐）** | 不支持 | 兼容支持（大规模纯线性稀疏问题推荐优先使用 `moonopt`） |

---

## 核心具体模型与客观对拍结果（13 个模型均经 SciPy + 全组合穷举独立验证）

### 1. 二阶锥规划（SOCP / QCQP）与混合整数二阶锥规划（MISOCP）— [benchmarks_socp.mbt](src/benchmarks_socp.mbt)

全部经 [verify/verify_socp_misocp.py](verify/verify_socp_misocp.py) 对拍通过：

| 具体模型 | 锥/离散结构特征 | `moonlp` 实测最优结果 | 独立证明 / SciPy+穷举误差 | `wasm-gc` 耗时 |
| :--- | :--- | :--- | :---: | ---: |
| **波动率硬上限基数约束投资组合（MISOCP）** | $\|\Sigma^{1/2}w\|_2 \le 10.5\%$ 锥约束 + 6 选 $\le 3$ 资产 + 最小持仓 10% | 最大收益 **`8.11046%`**（11 节点；**US_Eq 43.21% / Corp_Bond 46.79% / Gold 10.00%**） | vs 41 子集穷举 **`7.14e-11`** | **`15.03 ms`** |
| **2D 圆推力矢量火箭制导（SOCP & MISOCP）** | 二次 LQR 目标（自动 epigraph 提升）+ $\sqrt{u_x^2+u_y^2} \le u_{\max} z_t$ + $\le 2$ 次点火 | 连续 SOCP **`69.742457`**；2 次点火 MISOCP **`75.203126`** (`[1,1,0,0]`) | vs SciPy+穷举 **`3.29e-8` / `9.10e-7`** | **`40.19 ms`** |
| **椭球不确定集鲁棒投资组合（Robust SOCP）** | $\max \mu^T w - 1.5 \|\Sigma^{1/2} w\|_2$ | 最坏情形收益 **`-0.02622928`**（风险 `5.350%`） | Conic KKT=`true`；**`1.16e-12`** | **`165.27 µs`** |
| **95% 高斯机会约束生产规划（SOCP）** | $\Pr(\tilde{a}^T x \le 100) \ge 95\%$ 精确等价二阶锥 | 安全最大利润 **`2046.9777`** | Conic KKT=`true`；**`7.93e-8`** | **`62.31 µs`** |
| **最小包围圆 / 5G 基站覆盖半径（SOCP）** | 6 个 Lorentz 锥 $\|x - p_i\|_2 \le r$ | 最小半径 **`2.32142858`**，圆心 `(2.000, 1.1786)` | Conic KKT=`true`；**`7.21e-9`** | **`2.00 ms`** |

### 2. 混合整数二次规划（MIQP）与凸 QP — [benchmarks_qp.mbt](src/benchmarks_qp.mbt)

全部经 [verify/verify_qp_miqp.py](verify/verify_qp_miqp.py) 对拍通过：

| 具体模型 | 类型与结构特征 | `moonlp` 全局最优结果 | 穷举+SciPy 误差 | `wasm-gc` 耗时 |
| :--- | :--- | :--- | :---: | ---: |
| **电力机组组合与二次经济调度（UC-MIQP）** | 5 机组启停 + 140 MW 负荷 + 20 MW 备用 + 二次煤耗 | 总成本 **`3308.30`**（`Coal 80MW + CCGT 55MW + Peaker2 5MW`） | **`1.55e-11`** | **`1.78 ms`** |
| **基数约束 Markowitz 资产配置（MIQP）** | 6 选 $\le 3$ 资产 + 最小持仓 10% + 秩 6 半正定协方差 | 最小方差 **`0.00520192`**（`US_Eq 41% / Corp_Bond 49% / Gold 10%`） | **`4.37e-14`** | **`5.17 ms`** |
| **$L_0$ 基数约束稀疏岭回归（MIQP）** | 5 含噪相关特征限选 $\le 2$ 个非零系数 | 损失 **`0.23809507`**，精准恢复真实支撑集 **`{f0, f3}`** | **`1.59e-13`** | **`2.11 ms`** |
| **推力最小冲量死区混杂 MPC（MIQP-MPC）** | $u_t \in \{0\}\cup[0.3,1.0]$ + $\le 3$ 次点火 | 代价 **`41.228125`**，开关 **`[1.0, 1.0, OFF, OFF, OFF]`** | **`1.46e-11`** | **`619.16 µs`** |
| **Markowitz 有效前沿 / 8 步饱和 MPC / 秩 2 SVM 对偶（QP）** | 凸 QP + Wolfe 强对偶 KKT + `.qps`/`QUADOBJ`/`.lp` 往返 | Wolfe KKT 全认证通过 | **`5.76e-13..1.44e-11`** | **`59.9 µs..996 µs`** |

---

## 快速开始

```bash
moon update
moon check --deny-warn
moon test --deny-warn
moon run src/examples
moon run interop_moonopt
moon bench --release

# 独立 Python / SciPy (numpy>=1.26, scipy>=1.13) + 全组合穷举对拍验证（将 SciPy 参考解与记录的 MoonBit 结果交叉核验，已接入 CI `verify-scipy`）
python verify/verify_socp_misocp.py
python verify/verify_qp_miqp.py

# 跨后端全量验证（与 GitHub Actions CI 一致）
moon test --deny-warn --target wasm-gc && moon run interop_moonopt --target wasm-gc && moon bench --release --target wasm-gc
moon test --deny-warn --target js      && moon run interop_moonopt --target js      && moon bench --release --target js
moon test --deny-warn --target native  && moon run interop_moonopt --target native  && moon bench --release --target native
```

### 示例：波动率硬上限基数约束投资组合（MISOCP）与最小包围圆（SOCP）

```moonbit
// 1. 波动率硬上限 (vol <= 10.5%) + 最多选 3 只资产的混合整数二阶锥投资组合 (MISOCP)
let misocp = @moonlp.build_vol_capped_cardinality_misocp(max_assets=3, vol_cap=0.105)
let msol = misocp.solve()       // Max Return = 8.11046% (US_Eq 43.2%, Corp_Bond 46.8%, Gold 10.0%)

// 2. 最小包围圆 (SOCP) + 独立二阶锥 KKT 强对偶认证
let disk = @moonlp.build_smallest_enclosing_disk_socp()
let (dsol, dcert) = disk.solve_socp_certified() // r = 2.321429, Conic KKT certified_optimal == true

// 3. 2D 圆推力矢量火箭制导 (二次 LQR 目标 + 精确圆推力锥 ||u_t||_2 <= u_max * z_t, MISOCP)
let rocket = @moonlp.build_switched_isotropic_mpc_misocp(horizon=4, max_burns=2, u_max=1.0)
let rsol = rocket.solve()       // Optimal 2-burn cost = 75.203126
```

---

## 测试与工程验收（95 项工作区测试 + 18 项性能基准，`--deny-warn` 三后端与 SciPy CI 全绿）

- **`Freon793/moonopt@0.6.0` 实时互通与原地扩展测试（3 项，`interop_moonopt/`）**：`from_moonopt()` 内存转换与 MPS 无损交接等价校验、生产规划 LP 原地扩展 95% 机会约束 SOCP、投资组合线性基线原地扩展波动率上限 MISOCP、`moonopt` 拒绝 `.qps QUADOBJ` 与无法求解 CPLEX QCP `.lp` 边界校验、`moonopt` MILP 原地扩展二次热耗 MIQP 验证。
- **SOCP / 旋转二阶锥 / QCQP / MISOCP 引擎（12 项专项测试）**：最小包围圆 + 独立锥 KKT 认证、`add_quad_constraint` 自动二阶锥提升与非凸拦截、分离二阶锥 Conic Farkas 不可行证明、椭球鲁棒投资组合 SOCP、波动率硬上限基数约束 MISOCP、95% 高斯机会约束规划 SOCP、2D 圆推力火箭制导 SOCP & 2 次点火 MISOCP、旋转二阶锥双曲几何均值（`add_rotated_soc_constraint`）、二次目标+二阶锥真实 Wolfe-Conic 强对偶认证（`gap=2.02e-6`）、CPLEX QCP `.lp` 解析与双向无损往返 + 一般二阶锥导出显式拒收（`to_lp_string_checked` / `to_mps_string_checked`）、大半径有界锥（`||(x,y)||_2 <= 1e10`）不误报无界 + 验证锥衰退射线无界测试。
- **QP / MIQP 引擎（13 项专项测试）**：解析 2D 凸 QP、凹目标最大化与弱有效约束退化、纯等式快速路径、非凸/不可行/无界拦截、Markowitz 有效前沿、基数约束选股 MIQP、8 步饱和 MPC、秩 2 半正定 SVM 对偶、单纯形最小二乘与 `.lp`/`.qps` `QUADOBJ` 往返、纯整数二次追踪 MIQP、电力机组组合 UC-MIQP、$L_0$ 稀疏岭回归 MIQP、推力死区混杂 MPC MIQP。
- **LP / MILP 兼容基线（67 项回归测试）**：Netlib `ADLITTLE`/`AFIRO`、HiGHS `AVGAS`/`SMALL_MIP`、MIPLIB `FLUGPL`、GMI 割平面、Ruiz 平衡缩放与多轮对抗数值回归。

---

## 目录结构

```text
moonlp/
├── moon.mod                      # MoonBit 1.0 模块配置 (v1.1.0)
├── docs/
│   └── comparison_with_moonopt.md # 全生态排重对比（moonopt/moonqp）、差异矩阵与复审说明
├── .github/workflows/ci.yml      # wasm-gc / js / native 三后端 --deny-warn + bench 矩阵
├── src/
│   ├── types.mbt                 # 核心类型（SocConstraint、QuadTerm、NonConvex、LinearExpr）
│   ├── model.mbt                 # 6 类问题统一代数建模 DSL、主元 Cholesky 与独立验证器
│   ├── socp.mbt                  # Nesterov-Todd 自缩放 Mehrotra 原始-对偶二阶锥内点法 + Conic KKT 证明器
│   ├── misocp.mbt                # Best-Bound Branch-and-Bound MISOCP / MIQCQP 引擎
│   ├── benchmarks_socp.mbt       # 6 个 SOCP / QCQP / MISOCP 真实模型（鲁棒组合/波动率MISOCP/机会约束/圆推力火箭/最小包围圆）
│   ├── socp_wbtest.mbt           # 12 项 SOCP / QCQP / MISOCP / QCP 往返与边界专项测试
│   ├── qp.mbt / miqp.mbt         # 尺度归一化 Mehrotra 凸 QP 内核 & Best-Bound MIQP 引擎
│   ├── benchmarks_qp.mbt         # 8 个 QP / MIQP 模型（含 UC-MIQP、L0 岭回归、死区 MPC、Markowitz MIQP）
│   ├── qp_wbtest.mbt             # 13 项 QP / MIQP 与 .qps/.lp 格式往返测试
│   ├── certificate.mbt           # 线性 KKT + 二次 KKT / Wolfe 强对偶独立证明器
│   ├── lp_format.mbt / mps_format.mbt # CPLEX 二次 .lp & MPS / QPS (QUADOBJ/QMATRIX) 双向往返
│   ├── simplex.mbt / dual_simplex.mbt / scaling.mbt / milp.mbt / warm_bb.mbt / cuts.mbt / branch_and_cut.mbt / presolve.mbt / sensitivity.mbt
│   ├── moonlp_wbtest.mbt / review_regression_wbtest.mbt # 67 项 LP/MILP 单元测试与回归集
│   ├── moonlp_bench_wbtest.mbt   # 18 项官方 moon bench --release 性能基准（SOCP/MISOCP/QP/MIQP/LP/MILP）
│   └── examples/main.mbt         # 18 合 1 CLI 演示
└── verify/
    ├── verify_socp_misocp.py     # 独立 SciPy + 全组合穷举对拍脚本（6 个 SOCP/MISOCP 模型）
    └── verify_qp_miqp.py         # 独立 SciPy + 全组合穷举对拍脚本（7 个 QP/MIQP 模型）
```

## License
Apache-2.0
