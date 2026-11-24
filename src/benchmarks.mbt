///|
/// Authentic Netlib / HiGHS industrial MPS benchmarks and combinatorial MILP builders.

///|
/// Authentic Stanford SOL / Netlib AFIRO benchmark (28 rows, 32 cols, 88 nonzeros).
/// Canonical published optimal objective value: -464.75314285714285.
pub fn netlib_afiro_mps() -> String {
  #|NAME          AFIRO   SIZE: N=32, M=28, NZ=115
  #|ROWS
  #| E  R09
  #| E  R10
  #| L  X05
  #| L  X21
  #| E  R12
  #| E  R13
  #| L  X17
  #| L  X18
  #| L  X19
  #| L  X20
  #| E  R19
  #| E  R20
  #| L  X27
  #| L  X44
  #| E  R22
  #| E  R23
  #| L  X40
  #| L  X41
  #| L  X42
  #| L  X43
  #| L  X45
  #| L  X46
  #| L  X47
  #| L  X48
  #| L  X49
  #| L  X50
  #| L  X51
  #| N  COST
  #|COLUMNS
  #|    X01       X48               .301   R09                -1.
  #|    X01       R10              -1.06   X05                 1.
  #|    X02       X21                -1.   R09                 1.
  #|    X02       COST               -.4
  #|    X03       X46                -1.   R09                 1.
  #|    X04       X50                 1.   R10                 1.
  #|    X06       X49               .301   R12                -1.
  #|    X06       R13              -1.06   X17                 1.
  #|    X07       X49               .313   R12                -1.
  #|    X07       R13              -1.06   X18                 1.
  #|    X08       X49               .313   R12                -1.
  #|    X08       R13               -.96   X19                 1.
  #|    X09       X49               .326   R12                -1.
  #|    X09       R13               -.86   X20                 1.
  #|    X10       X45              2.364   X17                -1.
  #|    X11       X45              2.386   X18                -1.
  #|    X12       X45              2.408   X19                -1.
  #|    X13       X45              2.429   X20                -1.
  #|    X14       X21                1.4   R12                 1.
  #|    X14       COST              -.32
  #|    X15       X47                -1.   R12                 1.
  #|    X16       X51                 1.   R13                 1.
  #|    X22       X46               .109   R19                -1.
  #|    X22       R20               -.43   X27                 1.
  #|    X23       X44                -1.   R19                 1.
  #|    X23       COST               -.6
  #|    X24       X48                -1.   R19                 1.
  #|    X25       X45                -1.   R19                 1.
  #|    X26       X50                 1.   R20                 1.
  #|    X28       X47               .109   R22               -.43
  #|    X28       R23                 1.   X40                 1.
  #|    X29       X47               .108   R22               -.43
  #|    X29       R23                 1.   X41                 1.
  #|    X30       X47               .108   R22               -.39
  #|    X30       R23                 1.   X42                 1.
  #|    X31       X47               .107   R22               -.37
  #|    X31       R23                 1.   X43                 1.
  #|    X32       X45              2.191   X40                -1.
  #|    X33       X45              2.219   X41                -1.
  #|    X34       X45              2.249   X42                -1.
  #|    X35       X45              2.279   X43                -1.
  #|    X36       X44                1.4   R23                -1.
  #|    X36       COST              -.48
  #|    X37       X49                -1.   R23                 1.
  #|    X38       X51                 1.   R22                 1.
  #|    X39       R23                 1.   COST               10.
  #|RHS
  #|    B         X50               310.   X51               300.
  #|    B         X05                80.   X17                80.
  #|    B         X27               500.   R23                44.
  #|    B         X40               500.
  #|ENDATA
}

///|
/// Authentic HiGHS AVGAS aviation gasoline LP benchmark (11 rows, 8 cols, 38 nonzeros).
/// Canonical published optimal objective value: -7.75.
pub fn highs_avgas_mps() -> String {
  #|NAME          AVGAS SIZE: N=8, M=11, NZ=38
  #|ROWS
  #| G  R1
  #| G  R2
  #| G  R3
  #| G  R4
  #| G  R5
  #| G  R6
  #| G  R7
  #| G  R8
  #| G  R9
  #| G  R10
  #| N  COST
  #|COLUMNS
  #|    C1        R1           -1.0        R5         -1.0
  #|    C1        R7            2.0        R8          5.0
  #|    C2        R1           -1.0        R6         -1.0
  #|    C2        R9            1.0        R10         1.0
  #|    C2        COST         -2.0
  #|    C3        R2           -1.0        R5         -1.0
  #|    C3        R7            1.0        R8          3.0
  #|    C3        COST         -1.0
  #|    C4        R2           -1.0        R6         -1.0
  #|    C4        R9           -1.0
  #|    C4        COST         -3.0
  #|    C5        R3           -1.0        R5         -1.0
  #|    C5        R8           -3.0
  #|    C5        COST         -2.0
  #|    C6        R3           -1.0        R6         -1.0
  #|    C6        R9           -3.0        R10        -3.0
  #|    C6        COST         -4.0
  #|    C7        R4           -1.0        R5         -1.0
  #|    C7        R7           -1.0        R8         -1.0
  #|    C7        COST         -3.0
  #|    C8        R4           -1.0        R6         -1.0
  #|    C8        R9           -5.0        R10        -2.0
  #|    C8        COST         -5.0
  #|RHS
  #|    DEMANDS   R1           -1.0        R2         -1.0
  #|    DEMANDS   R3           -1.0        R4         -1.0
  #|    DEMANDS   R5           -2.0        R6         -2.0
  #|BOUNDS
  #| UP SERVINGS  C1            1.0
  #| UP SERVINGS  C2            1.0
  #| UP SERVINGS  C3            1.0
  #| UP SERVINGS  C4            1.0
  #| UP SERVINGS  C5            1.0
  #| UP SERVINGS  C6            1.0
  #| UP SERVINGS  C7            1.0
  #| UP SERVINGS  C8            1.0
  #|ENDATA
}

///|
/// Build a 5-city Traveling Salesman Problem (TSP) using Miller-Tucker-Zemlin (MTZ)
/// subtour elimination constraints (20 binary edge vars + 4 continuous order vars).
/// Exact optimal tour length = 80.0 (0 -> 1 -> 3 -> 4 -> 2 -> 0).
pub fn build_tsp5_mtz_model() -> Model {
  let dist = [
    [0.0, 10.0, 15.0, 20.0, 12.0],
    [10.0, 0.0, 35.0, 25.0, 18.0],
    [15.0, 35.0, 0.0, 30.0, 14.0],
    [20.0, 25.0, 30.0, 0.0, 16.0],
    [12.0, 18.0, 14.0, 16.0, 0.0],
  ]
  let n = 5
  let m = Model::new(Minimize)
  let x : Array[Array[VarId]] = []
  for i = 0; i < n; i = i + 1 {
    let row : Array[VarId] = []
    for j = 0; j < n; j = j + 1 {
      if i == j {
        row.push({ index: -1, })
      } else {
        row.push(
          m.add_binary_var(
            "x_" + i.to_string() + "_" + j.to_string(),
            obj=dist[i][j],
          ),
        )
      }
    }
    x.push(row)
  }
  let u : Array[VarId] = [{ index: -1, }]
  for i = 1; i < n; i = i + 1 {
    u.push(m.add_continuous_var("u_" + i.to_string(), lb=1.0, ub=4.0))
  }
  for i = 0; i < n; i = i + 1 {
    let out_terms : Array[Term] = []
    let in_terms : Array[Term] = []
    for j = 0; j < n; j = j + 1 {
      if i != j {
        out_terms.push({ var_id: x[i][j], coeff: 1.0, })
        in_terms.push({ var_id: x[j][i], coeff: 1.0, })
      }
    }
    m.add_constraint("out_" + i.to_string(), out_terms, Eq, 1.0)
    m.add_constraint("in_" + i.to_string(), in_terms, Eq, 1.0)
  }
  for i = 1; i < n; i = i + 1 {
    for j = 1; j < n; j = j + 1 {
      if i != j {
        m.add_constraint(
          "mtz_" + i.to_string() + "_" + j.to_string(),
          [
            { var_id: u[i], coeff: 1.0, },
            { var_id: u[j], coeff: -1.0, },
            { var_id: x[i][j], coeff: 4.0, },
          ],
          Le,
          3.0,
        )
      }
    }
  }
  m
}
