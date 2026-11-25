# `spring1219/moonlp` v1.0.0 与 `Freon793/moonopt` 生态定位、实质差异与新增模型验证报告

## 1. 背景与定位澄清

MoonBit 生态中已有近期持续维护的优秀纯 MoonBit 线性与混合整数线性规划求解器 **`Freon793/moonopt`**（基于稀疏修正单纯形 + 稀疏 LU 基分解 + Presolve/Postsolve + Farkas/最优独立证书 + B&B）。

在充分尊重并认可 `moonopt` 在**大规模稀疏纯线性 LP/MILP** 领域优势的基础上，`spring1219/moonlp` v1.0.0 明确将核心演进重心与生态净贡献定位于 **`moonopt` 明确不覆盖的非线性二次规划领域**：

> **补齐 MoonBit 生态在纯 MoonBit（零 FFI、全后端 `wasm-gc`/`js`/`native`）凸二次规划（Convex QP）、混合整数二次规划（MIQP）、二次 KKT / Wolfe 强对偶认证及 `.qps` / CPLEX 二次 `.lp` 格式互操作上的空白。**

---

## 2. 能力边界与实质差异对照

| 能力维度 | `Freon793/moonopt`（已有生态基线） | `spring1219/moonlp` v1.0.0（本项目净增量） | 生态关系 |
| :--- | :--- | :--- | :--- |
| **线性规划（LP）与 MILP** | 支持（稀疏修正单纯形 + 稀疏 LU，适合大规模稀疏线性模型） | 兼容支持（稠密稀疏感知单纯形 + Ruiz 缩放 + GMI 割平面 + B&B） | 互补：纯大规模稀疏 LP/MILP 推荐使用 `moonopt`；两者通过标准 `.mps` / `.lp` 无损互通 |
| **凸二次规划（Convex QP）** $\min \frac{1}{2}x^T Q x + c^T x$ | **不支持**（仅线性目标） | **原生支持**（尺度归一化 Mehrotra 预测-校正原始-对偶内点法 `solve_qp`） | **实质新增能力** |
| **半正定秩亏 Hessian（$Q \succeq 0$）与非凸拦截** | **不支持** | **原生支持**（平移 Cholesky 半正定校验、显式 `NonConvex` 状态拒绝非凸问题、秩亏 Gram/MIQP 松弛稳健收敛、无界退化射线证明） | **实质新增能力** |
| **混合整数二次规划（MIQP）** | **不支持** | **原生支持**（Best-Bound B&B + 节点凸 QP 松弛 + 固定整数连续 QP 修复启发式 `solve_miqp`） | **实质新增能力** |
| **二次 KKT 与 Wolfe 强对偶独立证明** | 仅线性 KKT / Farkas 证书 | **原生支持** `certify_qp_solution`（验证 $Q \succeq 0$、$\nabla f(x) = Qx + c - A^T y - z_l + z_u = 0$ 与 Wolfe 强对偶间隙 $-\frac{1}{2}x^TQx + y^Tb + z_l^Tl - z_u^Tu$） | **实质新增能力** |
| **工业二次格式互操作（`.qps` / CPLEX QP `.lp`）** | 仅线性 MPS / LP | **原生支持** MPS `QUADOBJ` / `QMATRIX` 与 CPLEX `[ ... ] / 2` 双向无损往返 | **实质新增能力** |

### 为何保留内部 LP 内核而非完全依赖外部线性后端？
1. **零依赖全后端自包含**：`solve_qp` 的可行性预检与热启动初值生成、以及 `.qps`/`.lp` 混合模型的统一入口 `Model::solve()` 可在零 FFI、零第三方依赖下同时编译到 `wasm-gc`、`js` 和 `native`。
2. **与 `moonopt` 的格式级互通与分工**：所有线性约束与标准模型均可通过 `to_mps_string()` / `to_lp_string()` 与 `Freon793/moonopt` 无缝交换；项目文档明确建议大规模纯线性稀疏问题优先使用 `moonopt`，而涉及二次目标、方差风险、LQR/MPC 控制、SVM 核训练或基数约束 MIQP 的场景使用 `moonlp`。

---

## 3. 具体优化模型与新增可用结果（均经独立 SciPy / 穷举对拍验证）

以下 5 个具体优化模型均包含二次目标项 $\frac{1}{2}x^T Q x$（其中模型 2 同时包含二值整数变量，属于 MIQP），属于 `moonopt` 无法直接建模与求解的问题族。所有代码位于 [`src/benchmarks_qp.mbt`](../src/benchmarks_qp.mbt)，可通过 `moon run src/examples` 与 `moon test --deny-warn` 一键复现，并经独立脚本 [`verify/verify_qp_miqp.py`](../verify/verify_qp_miqp.py)（SciPy SLSQP + 全组合穷举）验证。

| # | 具体优化模型 | 类型与结构特征 | `moonlp` 实测最优结果 | 独立证明 / SciPy 对拍误差 | `wasm-gc` 耗时 (`moon bench --release`) |
| :---: | :--- | :--- | :--- | :--- | ---: |
| **1** | **Markowitz 6 资产均值-方差投资组合** (`build_markowitz_portfolio_qp`, $r \ge 8\%$) | 连续凸 QP（$6\times 6$ 正定协方差 $\Sigma$，预算等式 + 收益下界 + 持仓上限） | $\frac{1}{2}w^T\Sigma w = \mathbf{0.004958677971}$（年化波动率 $9.959\%$，8 次内点迭代，分散持有 5 只资产） | Wolfe KKT = `true`（gap `5.89e-11`）；vs SciPy 误差 **`1.47e-12`** | **`57.35 µs`** |
| **2** | **基数约束真实投资组合选股** (`build_cardinality_portfolio_miqp`, 6 选 $\le 3$ 资产 + 最小持仓 10%) | **混合整数二次规划 MIQP**（6 连续权重 + 6 二值选股变量，$12\times 12$ **秩 6 半正定** Hessian） | 全局最优 $\frac{1}{2}w^T\Sigma w = \mathbf{0.005201920000}$（9 个 B&B 节点，`mip_gap = 2.16e-12`），精准选出 **US_Eq 41% + Corp_Bond 49% + Gold 10%** | vs 全 41 个子集穷举精确最优误差 **`4.37e-14`** | **`4.40 ms`** |
| **3** | **带输入/速度硬约束的线性模型预测控制 MPC** (`build_linear_mpc_qp`, $N=8$) | 连续凸 QP（24 决策变量 + 16 动力学等式约束 + 执行器饱和 $\|u_t\| \le 1.0$ + 速度安全界 $\|v_t\| \le 1.5$） | LQR 代价 $= \mathbf{38.7957968155}$；$u_0 = -1.0$（饱和限幅生效），终端位置由 $p_0=2.0$ 平滑收敛至 $p_8 = -0.0111$ | Wolfe KKT = `true`；vs SciPy 误差 **`5.76e-13`** | **`857.16 µs`**（>1,100 Hz 实时闭环） |
| **4** | **软间隔支持向量机（SVM）对偶训练** (`build_svm_dual_qp`, $C=10$) | **秩亏半正定 QP**（$8\times 8$ **秩 2** Gram 矩阵 $y_i y_j x_i^T x_j$ + 等式 $\sum y_i \alpha_i = 0$ + 盒约束 $0 \le \alpha_i \le C$） | 对偶最优 $= \mathbf{-0.0928936368}$（7 次内点迭代，检出 2 个支持向量，全部样本间隔 $y_i(w^Tx_i+b) \ge 1.0$） | Wolfe KKT = `true`；vs SciPy 误差 **`1.44e-11`** | **`74.12 µs`** |
| **5** | **单纯形约束最小二乘 / 光谱解混** (`build_constrained_lsq_qp`) | 连续凸 QP + `.qps` (`QUADOBJ`) / CPLEX `.lp` 二次格式往返验证 | 原模型、`.lp` 导出再解析、`.qps` 导出再解析三者最优值完全一致（误差 $< 10^{-7}$） | Wolfe KKT = `true`（三路径全通过） | — |

---

## 4. 复审申报表可直接引用的补充说明（精简版）

> 针对初审意见中指出的与生态已有 `Freon793/moonopt` 的职责重叠问题，我们对项目定位与核心能力进行了实质性拓展（升级至 `v1.0.0`），明确与 `moonopt` 的互补边界，补齐了 MoonBit 生态目前空白的 **纯 MoonBit 凸二次规划（QP）与混合整数二次规划（MIQP）** 求解及证明能力：
>
> 1. **与 `moonopt` 的清晰分工与实质差异**：`Freon793/moonopt` 专注于稀疏修正单纯形 LP/MILP 求解，不支持二次目标函数 $\frac{1}{2}x^T Q x$。本项目新增了零 FFI、全后端（`wasm-gc`/`js`/`native`）的 **尺度归一化 Mehrotra 预测-校正原始-对偶内点法凸 QP 求解器（`solve_qp`）**、**Best-Bound Branch-and-Bound 混合整数二次规划引擎（`solve_miqp`）**、**独立二次 KKT 与 Wolfe 强对偶证明器（`certify_qp_solution`）**，以及 **`.qps`（MPS `QUADOBJ`/`QMATRIX`）与 CPLEX 二次 `.lp`（`[ ... ] / 2`）** 双向无损解析与导出；同时通过标准 `.mps`/`.lp` 与 `moonopt` 保持格式级互通，并明确在文档中推荐纯大规模稀疏 LP/MILP 使用 `moonopt`。
> 2. **具体优化模型与新增可用结果（全部经独立 SciPy / 穷举对拍验证，误差 $\le 1.44\times 10^{-11}$）**：
>    - **Markowitz 6 资产均值-方差投资组合（QP）**：求得目标收益 $\ge 8\%$ 的最小方差组合 $\frac{1}{2}w^T\Sigma w = 0.00495868$，8 次内点迭代收敛，Wolfe 强对偶间隙 `5.89e-11`，`wasm-gc` 实测 **`57.35 µs`**；
>    - **基数约束 Markowitz 选股优化（MIQP，6 选 $\le 3$ 资产且最小持仓 10%）**：处理 $12\times 12$ 秩 6 半正定 Hessian + 6 二值选股变量，9 个 B&B 节点求得全局最优 `0.00520192`（`US_Eq 41% / Corp_Bond 49% / Gold 10%`，`mip_gap = 2.16e-12`），与全 41 个子集穷举最优误差 `4.37e-14`，`wasm-gc` 实测 **`4.40 ms`**；
>    - **带执行器饱和与速度约束的线性 MPC 控制（QP，8 步 24 变量 16 动力学等式）**：求得最优 LQR 代价 `38.795797`（$u_0 = -1.0$ 饱和限幅、终端位置由 $2.0$ 收敛至 $-0.011$），`wasm-gc` 实测 **`857.16 µs`**（>1,100 Hz 实时频率）；
>    - **软间隔 SVM 对偶训练（秩 2 半正定 $8\times 8$ Gram 矩阵 QP）**：7 次迭代求得最优对偶值 `-0.09289364`，`wasm-gc` 实测 **`74.12 µs`**，独立 KKT 认证通过。
> 3. **工程验收**：全部 77 项单元测试（含 67 项 LP/MILP 回归 + 10 项 QP/MIQP 专项）与 10 项性能基准在 `wasm-gc` / `js` / `native` 三后端 `--deny-warn` CI 全绿通过，可运行 `moon run src/examples` 与 `python verify/verify_qp_miqp.py` 复现全部结果。
