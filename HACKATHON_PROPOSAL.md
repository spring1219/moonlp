# MoonBit 2026 开发者黑客松（11月赛季）项目申报表文案

> 可直接复制填入官方飞书申报表单（官网入口：`moonbitlang.github.io` / `moonbitlang.cn` 赛事页面）

---

### 1. 项目名称
**`moonlp` — 纯 MoonBit 线性规划（LP）与混合整数规划（MILP）求解器与建模库**

### 1.1 开源仓库地址
**https://github.com/spring1219/moonlp** （GitHub ID：`spring1219`）

### 2. 参赛赛道
**新项目方向 — 开源基础工具库 / 科学计算与运筹优化**

### 3. 项目定位与生态空白分析
经检索 `mooncakes.io` 现有包生态：
- MoonBit 社区已具备线性代数基础库（`moonmcl`、`Luna-Flow/linear-algebra`）、SAT 求解器（`moonbit-sat`）、向量与全文检索引擎（`moonhnsw`、`moonseek`）与策略引擎（`moonrego`、`moon_cel`）。
- 但**独立的运筹优化求解器（Linear Programming / Mixed-Integer Linear Programming）目前在 `mooncakes.io` 上仍属空白**（类似 Rust 生态的 `microlp` / `good_lp` / `highs`）。
- 依托 MoonBit 编译至 Wasm / JS / Native 的极小产物与高性能优势，零外部依赖的 `moonlp` 可直接部署于浏览器、边缘节点与服务端，服务于排班调度、供应链排产、网络流优化、组合指派、游戏数值平衡与编译器寄存器分配等真实工程场景。

### 4. 核心技术方案与已验收成果（v0.2.0）
本项目全部采用纯 MoonBit 1.0 实现（零外部依赖，`moon check` 0 error 0 warning，`moon test` 20/20 全绿）：
1. **类型安全代数建模 DSL**：支持 `Minimize`/`Maximize`、连续/一般整数/0-1 二值变量、任意上下界与自由变量、内置独立可行性校验器。
2. **Two-Phase Simplex 核心引擎**：标准型自动变换 + Phase I/II 两阶段单纯形法 + Dantzig 选元与 **Bland 防退化循环规则**，精确输出原始解、对偶影子价格与检验数。
3. **Best-Bound Branch-and-Bound MILP 引擎**：支持二值/整数/连续混合规划、最分数分支、界剪枝与 MIP Gap 计算。
4. **Presolve 预处理器 + 可机检 KKT 强对偶证明器**：自动收紧单变量约束界并消除冗余行；对 LP 最优解独立校验互补松弛与强对偶间隙 $|c^T x - b^T y|$。
5. **双工业标准格式解析与导出**：完整支持 **CPLEX `.lp` 解析/导出** 与 **Netlib/MIPLIB `.mps` 解析**（含 `INTORG`/`INTEND` 整数段）。
6. **20 项客观测试集与可运行 CLI 示例**：涵盖 Beale 1955 退化循环 LP、Klee-Minty 4D 压测、CLRS 网络最大流、供应链排产、0-1 背包、设施选址、64 变量 4×4 数独 MILP 实时求解演示（`moon run src/examples` 76 次旋转秒解）。

### 5. 后续迭代计划（冲刺季度赛）
- 实现 Dual Simplex 热启动，显著加速大规模 MILP 分支定界子节点重优化。
- 引入稀疏矩阵 LU 分解与 Gomory 分数割平面（Cutting Planes）。
- 对接 Netlib LP 公开算例集补充 `moon bench` 跨后端（Wasm / JS / Native）性能基准并发布至 `mooncakes.io`。
