# `moonlp` v1.1.0 — Pure MoonBit SOCP / MISOCP / QCQP & MIQP Optimization Kernel

[![MoonBit CI](https://github.com/spring1219/moonlp/actions/workflows/ci.yml/badge.svg)](https://github.com/spring1219/moonlp/actions/workflows/ci.yml)
![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)
![Backends: wasm-gc | js | native](https://img.shields.io/badge/targets-wasm--gc%20%7C%20js%20%7C%20native-brightgreen.svg)

> **零外部依赖、全后端（`wasm-gc` / `js` / `native`）通用的纯 MoonBit 二阶锥规划（SOCP）、混合整数二阶锥规划（MISOCP）、二次约束规划（QCQP）与混合整数二次规划（MIQP）求解及强对偶证明内核**（同时兼容标准 LP / MILP 与 `.mps` / `.qps` / `.lp` 工业格式）。

---

## 与 MoonBit 生态已有优化库（`moonopt` / `moonqp`）的本质差异

详细全生态排重对比见 [docs/comparison_with_moonopt.md](docs/comparison_with_moonopt.md)：

| 能力维度 | `Freon793/moonopt`（线性稀疏 LP/MILP） | `yujieMamba/moonqp`（连续 ADMM QP） | **`spring1219/moonlp` v1.1.0（本项目）** |
| :--- | :---: | :---: | :---: |
| **二阶锥规划（SOCP）** $\|A_k x+b_k\|_2 \le c_k^T x+d_k$ | 不支持 | 不支持 | **原生支持**（Nesterov-Todd 自缩放 Mehrotra 原始-对偶锥内点法 `solve_socp`） |
| **凸二次约束规划（QCQP）** $\frac{1}{2}x^T P x + q^T x \le r$ | 不支持 | 不支持 | **原生支持**（`add_quad_constraint` 主元 Cholesky 自动无损二阶锥提升） |
| **混合整数二阶锥（MISOCP / MIQCQP）** | 不支持 | **明确不支持** | **原生支持**（`solve_misocp`：B&B + 节点 SOCP 松弛 + 活动界传播） |
| **混合整数二次规划（MIQP）** | 不支持 | **明确不支持** | **原生支持**（`solve_miqp`：B&B + 节点凸 QP 松弛 + 固定整数修复） |
| **锥不可行（Conic Farkas）/ 无界 / 非凸精确判定** | 仅 LP | **不支持** | **原生支持**（`Infeasible` / `Unbounded` / `NonConvex` / `Optimal`） |
| **独立锥 KKT & Wolfe 强对偶证明** | 仅线性 | 无强对偶证明 | **原生支持** `solve_socp_certified` & `certify_qp_solution` |
| **`.qps` (`QUADOBJ`/`QMATRIX`) & CPLEX 二次 `.lp`** | 仅线性 MPS/LP | **不支持** | **原生支持双向无损往返**（与 `moonopt` 标准 `.mps`/`.lp` 互通） |
| **大规模稀疏纯线性 LP/MILP** | **强项（推荐）** | 不支持 | 兼容支持（详见对比文档） |

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
moon check --deny-warn
moon test --deny-warn
moon run src/examples
moon bench --release

# 独立 Python / SciPy + 全组合穷举对拍验证（全部 13 个 SOCP/MISOCP/QP/MIQP 模型）
python verify/verify_socp_misocp.py
python verify/verify_qp_miqp.py

# 跨后端全量验证（与 GitHub Actions CI 一致）
moon test --deny-warn --target wasm-gc && moon bench --release --target wasm-gc
moon test --deny-warn --target js      && moon bench --release --target js
moon test --deny-warn --target native  && moon bench --release --target native
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

## 测试与工程验收（87 项单元测试 + 18 项性能基准，`--deny-warn` 三后端全绿）

- **SOCP / QCQP / MISOCP 引擎（7 项专项测试）**：最小包围圆 + 独立锥 KKT 认证、`add_quad_constraint` 自动二阶锥提升与非凸拦截、分离二阶锥 Conic Farkas 不可行证明、椭球鲁棒投资组合 SOCP、波动率硬上限基数约束 MISOCP、95% 高斯机会约束规划 SOCP、2D 圆推力火箭制导 SOCP & 2 次点火 MISOCP。
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
│   ├── socp_wbtest.mbt           # 7 项 SOCP / QCQP / MISOCP 专项测试
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
