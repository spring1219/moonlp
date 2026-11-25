# `spring1219/moonlp` v1.0.0 全生态排重对比、实质差异与新增模型验证报告

## 1. MoonBit 生态已有优化项目全面排重（截至 2026-11-26）

针对初审意见，我们对 GitHub 与 `mooncakes.io` 上的全部 MoonBit 优化相关仓库进行了系统排重，客观梳理各项目能力边界与真实生态空白：

| 生态已有项目 | 最近维护 | 核心定位与算法 | 明确边界 / 不支持项 |
| :--- | :---: | :--- | :--- |
| **`Freon793/moonopt`** | 2026-11-13 | 纯线性 **LP / MILP**（稀疏修正单纯形 + 稀疏 LU + Presolve + Farkas/最优证书 + B&B + MPS/LP） | **仅支持线性目标**，不支持二次目标 $\frac{1}{2}x^T Q x$、凸 QP 或混合整数二次规划（MIQP） |
| **`yujieMamba/moonqp`** | 2026-10-13 | **纯连续稠密凸 QP**（一阶 ADMM 分裂，稠密数组矩阵接口 `QpProblem::new`，含连续 MPC/估计/分配示例） | 其 README 明文声明以下均在支持范围之外：**不支持混合整数规划（MIQP）**、**不支持不可行/无界判定与证明**（`IterationLimit` 不代表不可行）、不支持 `.qps`/`.lp` 工业格式、不支持代数命名建模 DSL、无 Wolfe 强对偶自动认证 |
| **`Luna-Flow/linear-program`** | 2026-11-12 | 轻量级 LP 建模与基础线性规划工具 | 无 MILP、无 QP / MIQP |
| **`paipai-Studio/BlackBoxOptimization`** | 2026-11-12 | 黑盒元启发式 / 群智能优化框架 | 无精确数学规划、无约束全局最优保证、无 KKT 证明 |

### 本项目 `spring1219/moonlp` v1.0.0 的不可替代净贡献

> **填补整个 MoonBit 生态在「混合整数二次规划（MIQP）+ 二阶原始-对偶内点法精确凸 QP（含不可行/无界判定与 Wolfe 强对偶 KKT 认证）+ 工业 `.qps` (`QUADOBJ`/`QMATRIX`) & CPLEX 二次 `.lp` 无损格式互操作」上的完整空白。**

---

## 2. 核心能力矩阵对比（`moonopt` vs `moonqp` vs `moonlp` v1.0.0）

| 能力维度 | `Freon793/moonopt` | `yujieMamba/moonqp` | `spring1219/moonlp` v1.0.0 |
| :--- | :---: | :---: | :---: |
| **混合整数二次规划（MIQP）** | 不支持 | **明确不支持**（README 声明 outside contract） | **原生支持**（`solve_miqp`：Best-Bound B&B + 节点凸 QP 松弛 + 固定整数连续 QP 修复 + `mip_gap`） |
| **连续凸 QP 核心算法** | 不支持 | 一阶 ADMM（稠密数组） | **二阶尺度归一化 Mehrotra 预测-校正原始-对偶内点法**（7~15 步收敛至 $10^{-11}$ 精度） |
| **半正定秩亏 Hessian（$Q \succeq 0$）与非凸拦截** | 不支持 | 支持（Jacobi PSD 校验） | **支持**（平移 Cholesky 校验 + 显式 `NonConvex` 状态 + 秩亏 Gram/MIQP 松弛稳健收敛） |
| **QP 不可行与无界退化射线精确判定** | 仅 LP | **不支持**（仅返回 `IterationLimit`） | **原生支持**（精确区分 `Infeasible` / `Unbounded` / `NonConvex` / `Optimal`） |
| **独立二次 KKT 与 Wolfe 强对偶证明** | 仅线性 KKT/Farkas | 仅提供残差评估器（无 Wolfe 强对偶证明） | **原生支持** `certify_qp_solution`（验证 $Q \succeq 0$、驻点、互补松弛与 Wolfe 强对偶间隙） |
| **工业二次格式 `.qps` (`QUADOBJ`/`QMATRIX`) & CPLEX QP `.lp`** | 仅线性 MPS/LP | **不支持** | **原生支持双向无损解析与导出** |
| **命名代数建模 DSL** | 支持（纯线性） | 不支持（原始稠密数组） | **统一支持 LP / MILP / QP / MIQP**，并通过标准 `.mps`/`.lp` 与 `moonopt` 互通 |
| **大规模稀疏纯线性 LP/MILP** | **强项（稀疏 LU 修正单纯形）** | 不支持 | 兼容支持；文档明确推荐大规模纯线性稀疏问题优先使用 `moonopt` |

---

## 3. 具体优化模型与新增可用结果（均经 SciPy + 全组合穷举独立对拍验证）

以下 **8 个具体优化模型** 全部内置于 [`src/benchmarks_qp.mbt`](../src/benchmarks_qp.mbt)，其中 **4 个为全生态任何已有库（包括 `moonopt` 与 `moonqp`）均无法直接求解的真实工业 MIQP 模型**，另外 4 个展示二阶内点法的高精度 Wolfe 强对偶认证与 `.qps`/`.lp` 格式互操作。全部结果可通过 `moon run src/examples`、`moon test --deny-warn` 与独立脚本 [`verify/verify_qp_miqp.py`](../verify/verify_qp_miqp.py) 一键复现。

### A. 混合整数二次规划（MIQP）模型 —— `moonopt` 与 `moonqp` 均无法求解

| # | 具体 MIQP 模型 | 实际业务背景与离散+二次结构 | `moonlp` 实测全局最优结果 | 独立全组合穷举+SciPy 对拍误差 | `wasm-gc` 实测耗时 |
| :---: | :--- | :--- | :--- | :---: | ---: |
| **1** | **电力系统机组组合与二次经济调度** (`build_unit_commitment_miqp`, 140 MW 负荷 + 20 MW 旋转备用) | 5 台火电机组二值启停 $u_g \in \{0,1\}$ + 固定开机费 + 最小/最大技术出力 + **真实二次煤耗曲线** $a_g p_g^2 + b_g p_g$ | 最优总成本 **`3308.30`**（9 节点，`mip_gap=0`），启机方案 **`Coal_Base 80MW + CCGT_Mid 55MW + Peaker_2 5MW`** | vs 全 31 种启停组合穷举误差 **`1.55e-11`** | **`1.58 ms`** |
| **2** | **基数约束 Markowitz 资产配置** (`build_cardinality_portfolio_miqp`, 6 选 $\le 3$ 只 + 最小持仓 10%) | 6 资产协方差矩阵 $\Sigma$ + 二值选股变量 $z_i \in \{0,1\}$ + 半连续持仓约束 $0.1 z_i \le w_i \le 0.5 z_i$（$12\times 12$ 秩 6 半正定） | 全局最优方差 **`0.00520192`**（9 节点，`mip_gap=2.16e-12`；对比连续松弛 5 资产 `0.00495868`），精选 **US_Eq 41% / Corp_Bond 49% / Gold 10%** | vs 全 41 个子集穷举误差 **`4.37e-14`** | **`4.34 ms`** |
| **3** | **$L_0$ 基数约束稀疏岭回归** (`build_best_subset_ridge_miqp`, 5 选 $\le 2$ 特征) | 含相关干扰特征的回归设计矩阵 $X$，施加严格 $L_0$ 约束 $\sum z_j \le 2$ 与 Big-M 开关 | 损失 **`0.23809507`**（7 节点），精准恢复真实支撑集 **`{f0, f3}`**（$\beta_0=2.493, \beta_3=-1.792$，干扰项全为 0） | vs 全特征子集精确线性求解误差 **`1.59e-13`** | **`1.82 ms`** |
| **4** | **带最小冲量死区与点火预算的混杂 MPC** (`build_hybrid_mpc_miqp`, 5 步 $\le 3$ 次点火) | 对比 `moonqp` 纯连续 MPC：推进器满足离散死区 $u_t \in \{0\} \cup [0.3, 1.0]$ + 阀门点火成本与次数上限 | 全局最优代价 **`41.228125`**，精准输出开关序列 **`[t0:1.0, t1:1.0, t2:OFF, t3:OFF, t4:OFF]`** | vs 全 32 种点火序列穷举误差 **`1.46e-11`** | **`580.05 µs`** |

### B. 连续凸 QP 模型与工业格式 / Wolfe KKT 证明

| # | 具体 QP 模型 | `moonlp` 实测最优结果 | Wolfe KKT 认证 / SciPy 对拍误差 | `wasm-gc` 耗时 |
| :---: | :--- | :--- | :---: | ---: |
| **5** | **Markowitz 6 资产均值-方差有效前沿** ($r \ge 8\%$) | $\frac{1}{2}w^T\Sigma w = \mathbf{0.0049586780}$（8 次内点迭代） | Wolfe KKT=`true` (gap `5.89e-11`)；vs SciPy **`1.47e-12`** | **`56.13 µs`** |
| **6** | **8 步带饱和/速度约束线性 MPC**（24 变量 16 等式） | LQR 代价 **`38.79579682`**（$u_0=-1.0$ 限幅，$p_8=-0.011$） | Wolfe KKT=`true`；vs SciPy **`5.76e-13`** | **`859.64 µs`** |
| **7** | **软间隔 SVM 对偶训练**（$8\times 8$ 秩 2 半正定 Gram） | 对偶最优 **`-0.09289364`**（7 次迭代，2 支持向量） | Wolfe KKT=`true`；vs SciPy **`1.44e-11`** | **`73.01 µs`** |
| **8** | **单纯形约束最小二乘 + `.lp` / `.qps` (`QUADOBJ`) 往返** | 原模型、CPLEX 二次 `.lp`、`.qps` `QUADOBJ` 三路径一致 | Wolfe KKT=`true`（三路径全绿） | — |

---

## 4. 复审申报表可直接引用的补充说明（精简版）

> 针对初审意见中指出的与生态已有 `Freon793/moonopt` 职责重叠的问题，我们对 MoonBit 生态全部优化库（包括纯线性 LP/MILP 求解器 `Freon793/moonopt` 与纯连续 ADMM QP 库 `yujieMamba/moonqp`）进行了全面排重，并将本项目升级至 **`v1.0.0`**，精准聚焦于 **现有生态明确缺失的混合整数二次规划（MIQP）与高精度二阶内点法凸 QP（含精确不可行/无界判定、Wolfe 强对偶认证及 `.qps`/`.lp` 工业二次格式）**：
>
> 1. **与 `moonopt`、`moonqp` 的清晰边界与实质差异**：
>    - `Freon793/moonopt` 专注大规模稀疏线性 LP/MILP，不支持二次目标；`yujieMamba/moonqp` 仅支持一阶 ADMM 连续稠密 QP，其 README 明文声明不支持混合整数模型（MIQP）、不支持不可行/无界判定、不支持 `.qps`/`.lp` 格式与代数建模 DSL。
>    - 本项目在零 FFI、三后端（`wasm-gc`/`js`/`native`）下新增：**二阶尺度归一化 Mehrotra 预测-校正原始-对偶内点法凸 QP 内核（`solve_qp`）**、**Best-Bound Branch-and-Bound 混合整数二次规划引擎（`solve_miqp`）**、**独立二次 KKT 与 Wolfe 强对偶证明器（`certify_qp_solution`）**、**`.qps` (`QUADOBJ`/`QMATRIX`) 与 CPLEX 二次 `.lp` (`[ ... ] / 2`) 无损往返**；并与 `moonopt` 保持 `.mps`/`.lp` 标准格式互通，在文档中明确推荐大规模纯线性稀疏问题优先使用 `moonopt`。
> 2. **4 个全生态无人能解的具体 MIQP 模型与新增可用结果（全部经独立 SciPy + 全组合穷举对拍，误差 $\le 1.55\times 10^{-11}$）**：
>    - **电力系统机组组合与二次经济调度（UC-MIQP，140 MW + 20 MW 备用）**：5 台火电机组二值启停 + 二次燃料热耗曲线，9 个 B&B 节点求得全局最优总成本 `3308.30`（`Coal_Base 80MW + CCGT_Mid 55MW + Peaker_2 5MW`，`mip_gap=0`），与全 31 种启停组合穷举误差 `1.55e-11`，`wasm-gc` 实测 **`1.58 ms`**；
>    - **基数约束 Markowitz 资产配置（MIQP，6 选 $\le 3$ 资产 + 最小持仓 10%）**：处理 $12\times 12$ 秩 6 半正定 Hessian + 6 二值变量，9 节点求得全局最优方差 `0.00520192`（`US_Eq 41% / Corp_Bond 49% / Gold 10%`），与全 41 子集穷举误差 `4.37e-14`，`wasm-gc` 实测 **`4.34 ms`**；
>    - **$L_0$ 基数约束稀疏岭回归（MIQP，5 选 $\le 2$ 特征）**：7 节点精准恢复真实支撑集 `{f0, f3}`（损失 `0.23809507`，误差 `1.59e-13`），`wasm-gc` 实测 **`1.82 ms`**；
>    - **带推进器最小冲量死区（$u_t \in \{0\}\cup[0.3,1.0]$）与点火次数上限的混杂 MPC（MIQP-MPC）**：精准求得全局最优代价 `41.228125` 与最优开关序列 `[1.0, 1.0, OFF, OFF, OFF]`（与全 32 种序列穷举误差 `1.46e-11`），`wasm-gc` 实测 **`580.05 µs`**；
>    - 此外还包含 Markowitz 连续有效前沿（`56.13 µs`，Wolfe 对偶间隙 `5.89e-11`）、8 步饱和线性 MPC（`859.64 µs`）、秩 2 半正定 SVM 对偶训练（`73.01 µs`）及 `.qps`/`.lp` 往返验证。
> 3. **工程验收**：全部 **80 项单元测试**（67 项 LP/MILP + 13 项 QP/MIQP）与 **13 项 `moon bench --release` 性能基准**在 `wasm-gc` / `js` / `native` 三后端 `--deny-warn` CI 全绿通过，可运行 `moon run src/examples` 与 `python verify/verify_qp_miqp.py` 一键复现。
