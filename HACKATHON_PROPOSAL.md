# MoonBit 2026 开发者黑客松（11月赛季）项目申报表文案

> 可直接复制填入官方飞书申报表单（官网入口：`moonbitlang.github.io` / `moonbitlang.cn` 赛事页面）

---

### 1. 项目名称
**`moonlp` — 纯 MoonBit 线性规划（LP）与混合整数规划（MILP）求解器与建模库**

### 2. 参赛赛道
**新项目方向 — 开源基础工具库 / 科学计算与运筹优化**

### 3. 项目定位与生态空白分析
经检索 `mooncakes.io` 现有生态：
- MoonBit 社区目前已拥有较完善的线性代数基础库（如 `moonmcl`、`Luna-Flow/linear-algebra`）、SAT 求解器（`moonbit-sat`）、向量与全文检索引擎（`moonhnsw`、`moonseek`）与策略引擎（`moonrego`、`moon_cel`）。
- 但**独立的运筹优化求解器（Linear Programming / Mixed-Integer Linear Programming）在 `mooncakes.io` 上尚属空白**（类似 Rust 生态中的 `microlp` / `good_lp` / `highs`）。
- 借助 MoonBit 编译到 Wasm / JS / Native 的极小体积与高性能优势，一个零外部依赖的纯 MoonBit LP/MILP 求解器可直接运行在浏览器端、边缘节点与服务端，服务于排班调度、资源分配、路径规划、组合优化、游戏数值平衡与编译器寄存器分配等真实工业场景。

### 4. 核心技术方案与交付内容
1. **类型安全代数建模 DSL**：支持 `Minimize`/`Maximize`、连续/整数/二值变量、任意上下界与自由变量、独立可行性与目标函数验证器。
2. **Two-Phase Simplex 核心引擎**：标准型自动变换 + Phase I/II 两阶段单纯形法 + Dantzig 选元与 Bland 防退化循环规则，输出精确原始解、对偶影子价格（Shadow Prices）与检验数。
3. **Best-Bound Branch-and-Bound MILP 引擎**：支持二值/一般整数/连续混合规划、最分数变量分支、边界剪枝与可证明的 MIP Gap 报告。
4. **工业标准互操作**：支持 CPLEX `.lp` 格式导出（后续迭代完整 `.lp`/`.mps` 解析器、Dual Simplex 热启动与 Gomory 割平面）。
5. **全量客观数学验证集**：覆盖强对偶校验（$c^T x = b^T y$）、Beale 1955 退化循环 LP、Klee-Minty 4D 立方体、0-1 背包、指派问题、设施选址等 11+ 经典算例。

### 5. 里程碑计划
- **M1（已完成核心原型）**：完成 DSL、Two-Phase Simplex + Bland 规则、Branch-and-Bound MILP、CPLEX `.lp` 导出器及 11 个经典数学算例全量通过。
- **M2（11月底验收版）**：在本地 MoonBit 工具链完成 `moon check` / `moon test` / `moon bench` 全绿验证，补充 Dual Simplex 热启动、CPLEX `.lp` 解析器、完整文档与多场景 Example，发布至 GitHub 与 `mooncakes.io`。
- **M3（季度赛冲刺版）**：引入稀疏矩阵 LU 分解（Forrest-Tomlin 更新）、Gomory 割平面与 Netlib 标准测试集基准跑分。
