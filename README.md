# `moonlp` — Pure MoonBit LP & MILP Optimization Solver

> **零外部依赖、全后端（Wasm / JS / Native）通用的纯 MoonBit 线性规划（LP）与混合整数线性规划（MILP）求解器与建模库**  
> 填补 `mooncakes.io` 生态中运筹优化求解器（Simplex / MILP / KKT 强对偶证明 / CPLEX `.lp` / Netlib `.mps`）的空白。

---

## 核心特性（v0.2.0 · MoonBit 1.0 全绿验证）

- **类型安全的代数建模 DSL** (`src/types.mbt`, `src/model.mbt`)
  - 支持 `Minimize` / `Maximize`、连续变量 (`Continuous`)、一般整数变量 (`Integer`)、0-1 二值变量 (`Binary`)。
  - 支持任意变量上下界（负下界、单侧无界、自由变量 $x \in (-\infty, +\infty)$、固定变量）与 `<=` / `>=` / `==` 约束。
  - 内置独立解验证器 `Model::verify_solution`，自动检验原始可行性、边界、整数性与目标函数值。
- **两阶段单纯形法核心引擎** (`src/simplex.mbt`)
  - 标准型自动变换 + **Phase I / Phase II Two-Phase Simplex**。
  - **Dantzig 最速下降选元 + 自动切换 Bland 防退化循环规则**：在极端病态退化算例（Beale 1955 循环 LP）下保证有限步收敛。
  - 精确输出原始最优解 (`primal`)、**对偶影子价格 (`dual`)**、检验数 (`reduced_costs`)，以及严谨的 `Optimal` / `Infeasible` / `Unbounded` 判定。
- **分支定界 MILP 引擎** (`src/milp.mbt`)
  - **Best-Bound（最优界优先）Branch-and-Bound** + 最分数变量分支 + 边界剪枝，支持二值/一般整数/连续混合规划与可证明 `mip_gap`。
- **Presolve 预处理器** (`src/presolve.mbt`)
  - 单变量约束自动收紧上下界、整数变量上下界取整、冗余行剔除与不可行提前短路。
- **可机检 KKT 强对偶最优性证书** (`src/certificate.mbt`)
  - 独立校验原始可行性误差、互补松弛误差 $|y_i(a_i^T x - b_i)|$ 与强对偶间隙 $|c^T x - b^T y|$。
- **双工业标准文件格式支持** (`src/lp_format.mbt`, `src/mps_format.mbt`)
  - **CPLEX `.lp`**：支持完整双向解析 (`parse_lp`) 与导出 (`to_lp_string`)。
  - **Netlib / MIPLIB `.mps`**：支持完整解析 (`parse_mps`)，含 `N/L/G/E` 行、`LO/UP/FX/FR/BV/LI/UI` 边界与 `'MARKER' 'INTORG'/'INTEND'` 整数段。

---

## 快速开始

### 运行内置示例与 20 项全量测试

```bash
moon check
moon test
moon run src/examples
```

示例实测输出：
```text
=== moonlp: Pure MoonBit LP/MILP Optimization Solver Demo ===

[1] Production Mix LP (Maximize 3*A + 2*B)
    Optimum = 9 | A = 1, B = 3
    Shadow Prices = [labor: 1, machine: 1]
    KKT Strong Duality Certified = true (gap = 0)

[2] CPLEX .lp Format Parser + Integer Programming
    Optimum = 40 (x=0, y=5, B&B nodes=9)

[3] 4x4 Sudoku Solved as 64-Binary-Variable MILP (pivots=76):
    1 2 | 3 4
    3 4 | 1 2
    ----+----
    2 1 | 4 3
    4 3 | 2 1
```

### 代码示例

```moonbit
// 1. 代数建模 + KKT 强对偶验证
let m = @moonlp.Model::new(Maximize)
let x = m.add_continuous_var("x", obj=3.0)
let y = m.add_continuous_var("y", obj=2.0)
m.add_constraint("labor", [{ var_id: x, coeff: 1.0 }, { var_id: y, coeff: 1.0 }], Le, 4.0)
m.add_constraint("machine", [{ var_id: x, coeff: 2.0 }, { var_id: y, coeff: 1.0 }], Le, 5.0)

let sol = m.solve()
let cert = m.certify_lp_solution(sol) // cert.certified_optimal == true, gap == 0.0

// 2. 直接解析工业标准 CPLEX .lp 或 Netlib .mps 文件求解
let model_from_lp = @moonlp.parse_lp(lp_text)
let model_from_mps = @moonlp.parse_mps(mps_text)
```

---

## 20 项客观数学验证测试清单（`moon test`: 20 passed, 0 failed）

| 类别 | 覆盖测试算例 |
| :--- | :--- |
| **基础与对偶 LP** | 2D 生产规划 + 影子价格、Diet `>=` 最小化 + 强对偶、自由变量与负下界、不可行与无界证明、KKT 强对偶与互补松弛证书 |
| **经典压测与图论 LP** | **Beale 1955 退化循环 LP**（Bland 规则收敛）、**Klee-Minty 4D 立方体**、CLRS 网络最大流 LP（精确解 `23`）、多周期供应链排产与库存 LP |
| **混合整数规划 MILP** | 0-1 背包、一般整数规划（突破 LP 松弛间隙）、3×3 组合指派、混合二值+连续设施选址、**64 变量 4×4 数独 MILP** |
| **工程与标准格式** | CPLEX `.lp` 解析与导出往返测试、Netlib `.mps` LP 与 `INTORG/INTEND` MILP 解析、非法格式 `ParseError` 容错、Presolve 界收紧与冗余行消除 |

---

## 目录结构

```text
moonlp/
├── moon.mod                  # MoonBit 1.0 模块配置
├── src/
│   ├── moon.pkg              # 主库包配置
│   ├── types.mbt             # 核心数据结构、LinearExpr 与 IEEE-754 数值工具
│   ├── model.mbt             # 代数建模 DSL 与解验证器
│   ├── simplex.mbt           # Two-Phase Simplex + Bland 防循环 + 对偶提取
│   ├── milp.mbt              # Best-Bound Branch-and-Bound MILP 引擎
│   ├── presolve.mbt          # Presolve 单变量界收紧与冗余行消除
│   ├── certificate.mbt       # KKT 强对偶与互补松弛最优性证明器
│   ├── lp_format.mbt         # CPLEX .lp 双向解析器与导出器
│   ├── mps_format.mbt        # Netlib / MIPLIB .mps 解析器
│   ├── moonlp_test.mbt       # 20 项全量单元测试集
│   └── examples/
│       ├── moon.pkg          # 可运行示例包配置
│       └── main.mbt          # CLI 综合演示（LP + .lp + 64 变量数独 MILP）
├── verify/verify_oracle.py   # 独立 Python 数学对拍校验脚本
└── HACKATHON_PROPOSAL.md     # MoonBit 2026 黑客松申报材料
```

## License
Apache-2.0
