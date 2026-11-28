///|
/// CPLEX .lp format exporter & parser for industry-standard interoperability
/// (compatible with HiGHS, Gurobi, GLPK, SCIP, and COIN-OR CBC).

///|
fn format_linear_terms(model : Model, terms : Array[Term]) -> String {
  if terms.length() == 0 {
    return "0"
  }
  let buf = StringBuilder()
  let mut emitted = false
  for i = 0; i < terms.length(); i = i + 1 {
    let t = terms[i]
    let c = t.coeff
    if c == 0.0 {
      continue
    }
    let name = model.vars[t.var_id.index].name
    if !emitted {
      emitted = true
      if c < 0.0 {
        buf.write_string("- ")
        if -c != 1.0 {
          buf.write_string((-c).to_string())
          buf.write_string(" ")
        }
      } else if c != 1.0 {
        buf.write_string(c.to_string())
        buf.write_string(" ")
      }
      buf.write_string(name)
    } else if c < 0.0 {
      buf.write_string(" - ")
      if -c != 1.0 {
        buf.write_string((-c).to_string())
        buf.write_string(" ")
      }
      buf.write_string(name)
    } else {
      buf.write_string(" + ")
      if c != 1.0 {
        buf.write_string(c.to_string())
        buf.write_string(" ")
      }
      buf.write_string(name)
    }
  }
  if !emitted {
    "0"
  } else {
    buf.to_string()
  }
}

///|
fn format_quad_block(model : Model) -> String {
  let q = model.build_symmetric_hessian()
  let n = model.vars.length()
  let buf = StringBuilder()
  let mut emitted = false
  for i = 0; i < n; i = i + 1 {
    for j = i; j < n; j = j + 1 {
      let c = if i == j { q[i][i] } else { 2.0 * q[i][j] }
      if c == 0.0 {
        continue
      }
      if !emitted {
        emitted = true
        if c < 0.0 {
          buf.write_string("- ")
          if -c != 1.0 {
            buf.write_string((-c).to_string())
            buf.write_string(" ")
          }
        } else if c != 1.0 {
          buf.write_string(c.to_string())
          buf.write_string(" ")
        }
      } else if c < 0.0 {
        buf.write_string(" - ")
        if -c != 1.0 {
          buf.write_string((-c).to_string())
          buf.write_string(" ")
        }
      } else {
        buf.write_string(" + ")
        if c != 1.0 {
          buf.write_string(c.to_string())
          buf.write_string(" ")
        }
      }
      buf.write_string(model.vars[i].name)
      if i == j {
        buf.write_string(" ^ 2")
      } else {
        buf.write_string(" * ")
        buf.write_string(model.vars[j].name)
      }
    }
  }
  buf.to_string()
}

///|
/// Export a Model (LP/MILP/QP/MIQP) to standard CPLEX `.lp` text format.
pub fn Model::to_lp_string(self : Model) -> String {
  let buf = StringBuilder()
  match self.sense {
    Minimize => buf.write_string("Minimize\n")
    Maximize => buf.write_string("Maximize\n")
  }
  let obj_terms : Array[Term] = []
  for j = 0; j < self.vars.length(); j = j + 1 {
    if self.vars[j].obj != 0.0 {
      obj_terms.push({ var_id: { index: j, }, coeff: self.vars[j].obj, })
    }
  }
  buf.write_string("  obj: ")
  buf.write_string(format_linear_terms(self, obj_terms))
  if self.obj_constant != 0.0 {
    if self.obj_constant > 0.0 {
      buf.write_string(" + ")
      buf.write_string(self.obj_constant.to_string())
    } else {
      buf.write_string(" - ")
      buf.write_string((-self.obj_constant).to_string())
    }
  }
  if self.has_quad_obj() {
    let qstr = format_quad_block(self)
    if qstr.length() > 0 {
      buf.write_string(" + [ ")
      buf.write_string(qstr)
      buf.write_string(" ] / 2")
    }
  }
  buf.write_string("\nSubject To\n")

  for i = 0; i < self.constraints.length(); i = i + 1 {
    let c = self.constraints[i]
    let cname = if c.name.length() > 0 { c.name } else { "c" + i.to_string() }
    let rel_str = match c.relation {
      Le => "<="
      Ge => ">="
      Eq => "="
    }
    buf.write_string("  ")
    buf.write_string(cname)
    buf.write_string(": ")
    buf.write_string(format_linear_terms(self, c.terms))
    buf.write_string(" ")
    buf.write_string(rel_str)
    buf.write_string(" ")
    buf.write_string(c.rhs.to_string())
    buf.write_string("\n")
  }

  buf.write_string("Bounds\n")
  let generals : Array[String] = []
  let binaries : Array[String] = []
  for j = 0; j < self.vars.length(); j = j + 1 {
    let v = self.vars[j]
    let has_lb = v.lb > -1.0e200
    let has_ub = v.ub < 1.0e200
    if !has_lb && !has_ub {
      buf.write_string("  ")
      buf.write_string(v.name)
      buf.write_string(" free\n")
    } else if has_lb && has_ub {
      buf.write_string("  ")
      buf.write_string(v.lb.to_string())
      buf.write_string(" <= ")
      buf.write_string(v.name)
      buf.write_string(" <= ")
      buf.write_string(v.ub.to_string())
      buf.write_string("\n")
    } else if has_lb {
      buf.write_string("  ")
      buf.write_string(v.lb.to_string())
      buf.write_string(" <= ")
      buf.write_string(v.name)
      buf.write_string("\n")
    } else {
      buf.write_string("  -inf <= ")
      buf.write_string(v.name)
      buf.write_string(" <= ")
      buf.write_string(v.ub.to_string())
      buf.write_string("\n")
    }
    match v.kind {
      Integer => generals.push(v.name)
      Binary => binaries.push(v.name)
      Continuous => ()
    }
  }

  if generals.length() > 0 {
    buf.write_string("Generals\n ")
    for i = 0; i < generals.length(); i = i + 1 {
      buf.write_string(" ")
      buf.write_string(generals[i])
    }
    buf.write_string("\n")
  }

  if binaries.length() > 0 {
    buf.write_string("Binaries\n ")
    for i = 0; i < binaries.length(); i = i + 1 {
      buf.write_string(" ")
      buf.write_string(binaries[i])
    }
    buf.write_string("\n")
  }

  buf.write_string("End\n")
  buf.to_string()
}

///|
/// ---- CPLEX .lp Format Parser ----

///|
priv enum LpSection {
  None
  Objective
  SubjectTo
  Bounds
  Generals
  Binaries
}

///|
fn is_space_char(c : Char) -> Bool {
  c == ' ' || c == '\t' || c == '\r'
}

///|
fn is_digit_char(c : Char) -> Bool {
  c >= '0' && c <= '9'
}

///|
fn is_ident_char(c : Char) -> Bool {
  (c >= 'a' && c <= 'z') ||
  (c >= 'A' && c <= 'Z') ||
  (c >= '0' && c <= '9') ||
  c == '_' ||
  c == '.' ||
  c == '$' ||
  c == '#'
}

///|
fn to_lower_ascii(s : String) -> String {
  let buf = StringBuilder()
  for ch in s {
    if ch >= 'A' && ch <= 'Z' {
      let lower = (ch.to_int() + 32).unsafe_to_char()
      buf.write_char(lower)
    } else {
      buf.write_char(ch)
    }
  }
  buf.to_string()
}

///|
fn strip_comment_and_trim(line : String) -> String {
  let chars : Array[Char] = []
  for ch in line {
    if ch == '\\' {
      break
    }
    chars.push(ch)
  }
  let mut start = 0
  let mut end = chars.length()
  while start < end && is_space_char(chars[start]) {
    start = start + 1
  }
  while end > start && is_space_char(chars[end - 1]) {
    end = end - 1
  }
  let buf = StringBuilder()
  for i = start; i < end; i = i + 1 {
    buf.write_char(chars[i])
  }
  buf.to_string()
}

///|
fn split_lines(text : String) -> Array[String] {
  let out : Array[String] = []
  let mut buf = StringBuilder()
  for ch in text {
    if ch == '\n' {
      out.push(buf.to_string())
      buf = StringBuilder()
    } else if ch != '\r' {
      buf.write_char(ch)
    }
  }
  let rem = buf.to_string()
  if rem.length() > 0 {
    out.push(rem)
  }
  out
}

///|
/// Tokenize a single LP expression / bound line into clean tokens.
fn tokenize_lp_line(line : String) -> Array[String] {
  let chars : Array[Char] = []
  for ch in line {
    chars.push(ch)
  }
  let n = chars.length()
  let tokens : Array[String] = []
  let mut i = 0
  while i < n {
    let c = chars[i]
    if is_space_char(c) {
      i = i + 1
      continue
    }
    if c == ':' ||
      c == '+' ||
      c == '-' ||
      c == '[' ||
      c == ']' ||
      c == '*' ||
      c == '^' ||
      c == '/' {
      let b = StringBuilder()
      b.write_char(c)
      tokens.push(b.to_string())
      i = i + 1
      continue
    }
    if c == '<' || c == '>' || c == '=' {
      let b = StringBuilder()
      b.write_char(c)
      if i + 1 < n && chars[i + 1] == '=' {
        b.write_char('=')
        i = i + 2
      } else {
        i = i + 1
      }
      tokens.push(b.to_string())
      continue
    }
    // Number token (including scientific notation like 1.5e-3)
    if is_digit_char(c) ||
      (c == '.' && i + 1 < n && is_digit_char(chars[i + 1])) {
      let b = StringBuilder()
      while i < n && (is_digit_char(chars[i]) || chars[i] == '.') {
        b.write_char(chars[i])
        i = i + 1
      }
      if i < n && (chars[i] == 'e' || chars[i] == 'E') {
        b.write_char(chars[i])
        i = i + 1
        if i < n && (chars[i] == '+' || chars[i] == '-') {
          b.write_char(chars[i])
          i = i + 1
        }
        while i < n && is_digit_char(chars[i]) {
          b.write_char(chars[i])
          i = i + 1
        }
      }
      tokens.push(b.to_string())
      continue
    }
    // Identifier token
    if is_ident_char(c) {
      let b = StringBuilder()
      while i < n && is_ident_char(chars[i]) {
        b.write_char(chars[i])
        i = i + 1
      }
      tokens.push(b.to_string())
      continue
    }
    i = i + 1
  }
  tokens
}

///|
/// Parse an ASCII floating-point number (including optional sign and scientific exponent).
fn try_parse_double(tok : String) -> Double? {
  let low = to_lower_ascii(tok)
  if low == "inf" || low == "infinity" || low == "+inf" || low == "+infinity" {
    return Some(1.0e300)
  }
  if low == "-inf" || low == "-infinity" {
    return Some(-1.0e300)
  }
  let chars : Array[Char] = []
  for ch in tok {
    chars.push(ch)
  }
  let n = chars.length()
  if n == 0 {
    return None
  }
  let mut i = 0
  let mut sign = 1.0
  if chars[i] == '+' {
    i = i + 1
  } else if chars[i] == '-' {
    sign = -1.0
    i = i + 1
  }
  let mut int_part = 0.0
  let mut frac_part = 0.0
  let mut frac_scale = 0.1
  let mut has_digits = false
  while i < n && is_digit_char(chars[i]) {
    has_digits = true
    int_part = int_part * 10.0 + (chars[i].to_int() - '0'.to_int()).to_double()
    i = i + 1
  }
  if i < n && chars[i] == '.' {
    i = i + 1
    while i < n && is_digit_char(chars[i]) {
      has_digits = true
      frac_part = frac_part +
        (chars[i].to_int() - '0'.to_int()).to_double() * frac_scale
      frac_scale = frac_scale * 0.1
      i = i + 1
    }
  }
  if !has_digits {
    return None
  }
  let mut exp_val = 0
  let mut exp_sign = 1
  if i < n && (chars[i] == 'e' || chars[i] == 'E') {
    i = i + 1
    if i < n && chars[i] == '+' {
      i = i + 1
    } else if i < n && chars[i] == '-' {
      exp_sign = -1
      i = i + 1
    }
    let mut has_exp_digits = false
    while i < n && is_digit_char(chars[i]) {
      has_exp_digits = true
      exp_val = exp_val * 10 + (chars[i].to_int() - '0'.to_int())
      i = i + 1
    }
    if !has_exp_digits {
      return None
    }
  }
  if i != n {
    return None
  }
  let mut val = int_part + frac_part
  let mut p = 0
  while p < exp_val {
    if exp_sign > 0 {
      val = val * 10.0
    } else {
      val = val * 0.1
    }
    p = p + 1
  }
  Some(sign * val)
}

///|
fn get_or_create_var(
  model : Model,
  name_to_id : Map[String, VarId],
  name : String,
) -> VarId {
  match name_to_id.get(name) {
    Some(id) => id
    None => {
      let id = model.add_continuous_var(name)
      name_to_id.set(name, id)
      id
    }
  }
}

///|
fn parse_linear_tokens(
  tokens : Array[String],
  start_idx : Int32,
  end_idx : Int32,
  model : Model,
  name_to_id : Map[String, VarId],
  allow_const~ : Bool = false,
) -> (Array[Term], Double) raise ParseError {
  let terms : Array[Term] = []
  let mut const_sum = 0.0
  let mut i = start_idx
  let mut sign = 1.0
  while i < end_idx {
    let tok = tokens[i]
    if tok == "+" {
      sign = 1.0
      i = i + 1
      continue
    }
    if tok == "-" {
      sign = -1.0
      i = i + 1
      continue
    }
    match try_parse_double(tok) {
      Some(num) => {
        if i + 1 < end_idx {
          let next_tok = tokens[i + 1]
          if next_tok != "+" && next_tok != "-" {
            let vid = get_or_create_var(model, name_to_id, next_tok)
            terms.push({ var_id: vid, coeff: sign * num, })
            sign = 1.0
            i = i + 2
            continue
          }
        }
        if allow_const {
          const_sum = const_sum + sign * num
        } else if num != 0.0 {
          raise ParseError(
            "Unexpected standalone constant in linear expression: " + tok,
          )
        }
        sign = 1.0
        i = i + 1
      }
      None => {
        let vid = get_or_create_var(model, name_to_id, tok)
        terms.push({ var_id: vid, coeff: sign, })
        sign = 1.0
        i = i + 1
      }
    }
  }
  (terms, const_sum)
}

///|
fn parse_signed_number(
  tokens : Array[String],
  idx : Int32,
) -> (Double, Int32) raise ParseError {
  if idx >= tokens.length() {
    raise ParseError("Expected number in LP bound/RHS")
  }
  if tokens[idx] == "+" || tokens[idx] == "-" {
    if idx + 1 >= tokens.length() {
      raise ParseError("Dangling sign in LP bound/RHS")
    }
    let sign = if tokens[idx] == "-" { -1.0 } else { 1.0 }
    match try_parse_double(tokens[idx + 1]) {
      Some(v) => (sign * v, idx + 2)
      None => raise ParseError("Invalid number after sign: " + tokens[idx + 1])
    }
  } else {
    match try_parse_double(tokens[idx]) {
      Some(v) => (v, idx + 1)
      None => raise ParseError("Invalid number: " + tokens[idx])
    }
  }
}

///|
fn parse_quad_tokens(
  tokens : Array[String],
  model : Model,
  name_to_id : Map[String, VarId],
  scale : Double,
) -> Unit raise ParseError {
  let n = tokens.length()
  let mut i = 0
  let mut sign = 1.0
  while i < n {
    let tok = tokens[i]
    if tok == "+" {
      sign = 1.0
      i = i + 1
      continue
    }
    if tok == "-" {
      sign = -1.0
      i = i + 1
      continue
    }
    let mut coeff = 1.0
    match try_parse_double(tok) {
      Some(num) => {
        coeff = num
        i = i + 1
        if i >= n {
          raise ParseError("Dangling coefficient in LP quadratic objective")
        }
      }
      None => ()
    }
    let v1_name = tokens[i]
    let v1 = get_or_create_var(model, name_to_id, v1_name)
    i = i + 1
    if i >= n {
      raise ParseError("Incomplete quadratic term after variable: " + v1_name)
    }
    if tokens[i] == "^" {
      if i + 1 >= n || tokens[i + 1] != "2" {
        raise ParseError("Only exponent ^ 2 is supported in QP objective")
      }
      model.add_quad_obj_term(v1, v1, sign * coeff * scale)
      i = i + 2
    } else if tokens[i] == "*" {
      if i + 1 >= n {
        raise ParseError("Missing second variable after * in QP objective")
      }
      let v2 = get_or_create_var(model, name_to_id, tokens[i + 1])
      model.add_quad_obj_term(v1, v2, sign * coeff * scale)
      i = i + 2
    } else {
      raise ParseError(
        "Expected '^ 2' or '* var' in QP objective near: " + tokens[i],
      )
    }
    sign = 1.0
  }
}

///|
fn parse_objective_tokens(
  tokens : Array[String],
  model : Model,
  name_to_id : Map[String, VarId],
) -> Unit raise ParseError {
  let n = tokens.length()
  let lin_tokens : Array[String] = []
  let mut i = 0
  while i < n {
    let tok = tokens[i]
    if tok == "[" {
      // Check if previous token in lin_tokens is a standalone outer sign (+ or -)
      let mut outer_sign = 1.0
      if lin_tokens.length() > 0 {
        let last = lin_tokens[lin_tokens.length() - 1]
        if last == "+" {
          let _ = lin_tokens.pop()
          outer_sign = 1.0
        } else if last == "-" {
          let _ = lin_tokens.pop()
          outer_sign = -1.0
        }
      }
      let q_tokens : Array[String] = []
      i = i + 1
      let mut closed = false
      while i < n {
        if tokens[i] == "]" {
          closed = true
          i = i + 1
          break
        }
        q_tokens.push(tokens[i])
        i = i + 1
      }
      if !closed {
        raise ParseError("Unclosed '[' in LP quadratic objective")
      }
      let mut scale = outer_sign
      if i + 1 < n && tokens[i] == "/" {
        match try_parse_double(tokens[i + 1]) {
          Some(denom) =>
            if denom == 0.0 {
              raise ParseError("Division by zero in LP quadratic objective")
            } else {
              scale = outer_sign / denom
              i = i + 2
            }
          None =>
            raise ParseError(
              "Invalid divisor after '/' in LP quadratic objective: " +
              tokens[i + 1],
            )
        }
      }
      parse_quad_tokens(q_tokens, model, name_to_id, scale)
    } else {
      lin_tokens.push(tok)
      i = i + 1
    }
  }
  if lin_tokens.length() > 0 {
    let (terms, obj_c) = parse_linear_tokens(
      lin_tokens,
      0,
      lin_tokens.length(),
      model,
      name_to_id,
      allow_const=true,
    )
    if obj_c != 0.0 {
      model.set_obj_constant(model.obj_constant + obj_c)
    }
    for k = 0; k < terms.length(); k = k + 1 {
      let t = terms[k]
      let old = model.vars[t.var_id.index]
      model.vars[t.var_id.index] = { ..old, obj: old.obj + t.coeff, }
    }
  }
}

///|
/// Parse standard CPLEX `.lp` text (LP/MILP/QP/MIQP) into a `Model`.
pub fn parse_lp(text : String) -> Model raise ParseError {
  let lines = split_lines(text)
  let model = Model::new(Minimize)
  let name_to_id : Map[String, VarId] = Map([])
  let mut section : LpSection = None
  let mut has_sense = false
  let mut auto_cname = 0
  let obj_all_tokens : Array[String] = []
  let mut obj_header_stripped = false
  let mut obj_parsed = false

  for li = 0; li < lines.length(); li = li + 1 {
    let clean = strip_comment_and_trim(lines[li])
    if clean.length() == 0 {
      continue
    }
    let low = to_lower_ascii(clean)

    // Detect section headers
    if low == "minimize" || low == "min" || low == "minimum" {
      model.set_sense(Minimize)
      section = Objective
      has_sense = true
      continue
    }
    if low == "maximize" || low == "max" || low == "maximum" {
      model.set_sense(Maximize)
      section = Objective
      has_sense = true
      continue
    }
    if low == "subject to" || low == "s.t." || low == "st" || low == "such that" {
      if !obj_parsed && obj_all_tokens.length() > 0 {
        parse_objective_tokens(obj_all_tokens, model, name_to_id)
        obj_parsed = true
      }
      section = SubjectTo
      continue
    }
    if low == "bounds" || low == "bound" {
      if !obj_parsed && obj_all_tokens.length() > 0 {
        parse_objective_tokens(obj_all_tokens, model, name_to_id)
        obj_parsed = true
      }
      section = Bounds
      continue
    }
    if low == "generals" ||
      low == "general" ||
      low == "gen" ||
      low == "integers" {
      if !obj_parsed && obj_all_tokens.length() > 0 {
        parse_objective_tokens(obj_all_tokens, model, name_to_id)
        obj_parsed = true
      }
      section = Generals
      continue
    }
    if low == "binaries" || low == "binary" || low == "bin" {
      if !obj_parsed && obj_all_tokens.length() > 0 {
        parse_objective_tokens(obj_all_tokens, model, name_to_id)
        obj_parsed = true
      }
      section = Binaries
      continue
    }
    if low == "end" {
      break
    }

    let tokens = tokenize_lp_line(clean)
    if tokens.length() == 0 {
      continue
    }

    match section {
      None =>
        raise ParseError(
          "Unexpected content before Minimize/Maximize: " + clean,
        )
      Objective => {
        let start = if !obj_header_stripped &&
          tokens.length() >= 2 &&
          tokens[1] == ":" {
          obj_header_stripped = true
          2
        } else {
          obj_header_stripped = true
          0
        }
        for k = start; k < tokens.length(); k = k + 1 {
          obj_all_tokens.push(tokens[k])
        }
      }
      SubjectTo => {
        let mut start = 0
        let mut cname = "c" + auto_cname.to_string()
        if tokens.length() >= 2 && tokens[1] == ":" {
          cname = tokens[0]
          start = 2
        } else {
          auto_cname = auto_cname + 1
        }
        // Locate relational operator
        let mut rel_idx = -1
        let mut rel : Relation = Le
        for k = start; k < tokens.length(); k = k + 1 {
          if tokens[k] == "<=" || tokens[k] == "<" {
            rel_idx = k
            rel = Le
            break
          }
          if tokens[k] == ">=" || tokens[k] == ">" {
            rel_idx = k
            rel = Ge
            break
          }
          if tokens[k] == "=" || tokens[k] == "==" {
            rel_idx = k
            rel = Eq
            break
          }
        }
        if rel_idx == -1 {
          raise ParseError(
            "Missing relational operator in constraint: " + clean,
          )
        }
        // Check for CPLEX QCP quadratic constraint `[ ... ] <= / >= rhs`
        let mut has_q_bracket = false
        for k = start; k < rel_idx; k = k + 1 {
          if tokens[k] == "[" {
            has_q_bracket = true
            break
          }
        }
        if has_q_bracket {
          let (rhs, rhs_end) = parse_signed_number(tokens, rel_idx + 1)
          if rhs_end != tokens.length() {
            raise ParseError(
              "Unexpected trailing tokens after QCP constraint RHS: " + clean,
            )
          }
          let lin_toks : Array[String] = []
          let q_pairs : Array[(VarId, VarId, Double)] = []
          let mut k = start
          while k < rel_idx {
            if tokens[k] == "[" {
              let mut outer_s = 1.0
              if lin_toks.length() > 0 {
                let last = lin_toks[lin_toks.length() - 1]
                if last == "+" {
                  let _ = lin_toks.pop()
                } else if last == "-" {
                  let _ = lin_toks.pop()
                  outer_s = -1.0
                }
              }
              k = k + 1
              let q_toks : Array[String] = []
              let mut closed = false
              while k < rel_idx {
                if tokens[k] == "]" {
                  closed = true
                  k = k + 1
                  break
                }
                q_toks.push(tokens[k])
                k = k + 1
              }
              if !closed {
                raise ParseError("Unclosed '[' in QCP constraint: " + clean)
              }
              // Parse quadratic terms inside `[ ... ]` directly into q_pairs
              let mut qi = 0
              let mut qs = 1.0
              while qi < q_toks.length() {
                if q_toks[qi] == "+" {
                  qs = 1.0
                  qi = qi + 1
                  continue
                }
                if q_toks[qi] == "-" {
                  qs = -1.0
                  qi = qi + 1
                  continue
                }
                let mut qc = 1.0
                match try_parse_double(q_toks[qi]) {
                  Some(num) => {
                    qc = num
                    qi = qi + 1
                    if qi >= q_toks.length() {
                      raise ParseError("Dangling coeff in QCP constraint")
                    }
                  }
                  None => ()
                }
                let v1 = get_or_create_var(model, name_to_id, q_toks[qi])
                qi = qi + 1
                if qi >= q_toks.length() {
                  raise ParseError("Incomplete quad term in QCP constraint")
                }
                if q_toks[qi] == "^" {
                  if qi + 1 >= q_toks.length() || q_toks[qi + 1] != "2" {
                    raise ParseError("Only ^ 2 supported in QCP constraint")
                  }
                  q_pairs.push((v1, v1, outer_s * qs * qc))
                  qi = qi + 2
                } else if q_toks[qi] == "*" {
                  if qi + 1 >= q_toks.length() {
                    raise ParseError("Missing var after * in QCP constraint")
                  }
                  let v2 = get_or_create_var(model, name_to_id, q_toks[qi + 1])
                  q_pairs.push((v1, v2, outer_s * qs * qc))
                  qi = qi + 2
                } else {
                  raise ParseError("Invalid quad term syntax in QCP constraint")
                }
                qs = 1.0
              }
            } else {
              lin_toks.push(tokens[k])
              k = k + 1
            }
          }
          let (lterms, _) = if lin_toks.length() > 0 {
            parse_linear_tokens(
              lin_toks,
              0,
              lin_toks.length(),
              model,
              name_to_id,
            )
          } else {
            ([], 0.0)
          }
          // Convert Ge to Le by negating everything; reject Eq on nonlinear quadratic constraints
          let dir_sign = match rel {
            Le => 1.0
            Ge => -1.0
            Eq =>
              raise ParseError(
                "Convex QCP requires <= or >= quadratic constraint: " + clean,
              )
          }
          let n_cur = model.vars.length()
          let q_vars : Array[VarId] = []
          for vi = 0; vi < n_cur; vi = vi + 1 {
            q_vars.push({ index: vi, })
          }
          let p_mat : Array[Array[Double]] = []
          for _vi = 0; _vi < n_cur; _vi = _vi + 1 {
            p_mat.push(Array::make(n_cur, 0.0))
          }
          for qp = 0; qp < q_pairs.length(); qp = qp + 1 {
            let (v1, v2, coeff) = q_pairs[qp]
            let c_eff = dir_sign * coeff
            if v1.index == v2.index {
              p_mat[v1.index][v1.index] = p_mat[v1.index][v1.index] + c_eff
            } else {
              p_mat[v1.index][v2.index] = p_mat[v1.index][v2.index] +
                0.5 * c_eff
              p_mat[v2.index][v1.index] = p_mat[v2.index][v1.index] +
                0.5 * c_eff
            }
          }
          let lexpr = LinearExpr::zero()
          for lt = 0; lt < lterms.length(); lt = lt + 1 {
            let _ = lexpr.add_term(
              lterms[lt].var_id,
              dir_sign * lterms[lt].coeff,
            )
          }
          let ok_qcp = model.add_quad_constraint(
            cname,
            q_vars,
            p_mat,
            lexpr,
            dir_sign * rhs,
            half=false,
          )
          if !ok_qcp {
            raise ParseError("Non-convex quadratic constraint in LP: " + clean)
          }
        } else {
          let (terms, _) = parse_linear_tokens(
            tokens, start, rel_idx, model, name_to_id,
          )
          let (rhs, rhs_end) = parse_signed_number(tokens, rel_idx + 1)
          if rhs_end != tokens.length() {
            raise ParseError(
              "Unexpected trailing tokens after constraint RHS: " + clean,
            )
          }
          model.add_constraint(cname, terms, rel, rhs)
        }
      }
      Bounds => {
        // Case 1: `x free`
        if tokens.length() == 2 && to_lower_ascii(tokens[1]) == "free" {
          let vid = get_or_create_var(model, name_to_id, tokens[0])
          model.set_var_bounds(vid, -1.0e300, 1.0e300)
          continue
        }
        // Find all relational operators in bound line
        let rel_positions : Array[Int32] = []
        for k = 0; k < tokens.length(); k = k + 1 {
          if tokens[k] == "<=" ||
            tokens[k] == ">=" ||
            tokens[k] == "=" ||
            tokens[k] == "<" ||
            tokens[k] == ">" {
            rel_positions.push(k)
          }
        }
        if rel_positions.length() == 2 {
          // Format: lb <= x <= ub
          let (lb, lb_end) = parse_signed_number(tokens, 0)
          if lb_end != rel_positions[0] ||
            rel_positions[0] + 2 != rel_positions[1] {
            raise ParseError("Invalid range bound declaration: " + clean)
          }
          let var_tok = tokens[rel_positions[0] + 1]
          let (ub, ub_end) = parse_signed_number(tokens, rel_positions[1] + 1)
          if ub_end != tokens.length() {
            raise ParseError("Unexpected trailing tokens in bound: " + clean)
          }
          let vid = get_or_create_var(model, name_to_id, var_tok)
          model.set_var_bounds(vid, lb, ub)
        } else if rel_positions.length() == 1 {
          let rp = rel_positions[0]
          let op = tokens[rp]
          // Check whether LHS is a variable name or a number
          let lhs_is_var = rp == 1 && try_parse_double(tokens[0]) is None
          if lhs_is_var {
            let vid = get_or_create_var(model, name_to_id, tokens[0])
            let (val, val_end) = parse_signed_number(tokens, rp + 1)
            if val_end != tokens.length() {
              raise ParseError("Unexpected trailing tokens in bound: " + clean)
            }
            let old = model.vars[vid.index]
            if op == "<=" || op == "<" {
              model.set_var_bounds(vid, old.lb, val)
            } else if op == ">=" || op == ">" {
              model.set_var_bounds(vid, val, old.ub)
            } else {
              model.set_var_bounds(vid, val, val)
            }
          } else {
            // Format: lb <= x
            let (val, val_end) = parse_signed_number(tokens, 0)
            if val_end != rp || rp + 2 != tokens.length() {
              raise ParseError("Invalid bound declaration: " + clean)
            }
            let var_tok = tokens[rp + 1]
            let vid = get_or_create_var(model, name_to_id, var_tok)
            let old = model.vars[vid.index]
            if op == "<=" || op == "<" {
              model.set_var_bounds(vid, val, old.ub)
            } else if op == ">=" || op == ">" {
              model.set_var_bounds(vid, old.lb, val)
            } else {
              model.set_var_bounds(vid, val, val)
            }
          }
        } else {
          raise ParseError("Invalid bound declaration: " + clean)
        }
      }
      Generals =>
        for k = 0; k < tokens.length(); k = k + 1 {
          let vid = get_or_create_var(model, name_to_id, tokens[k])
          let old = model.vars[vid.index]
          model.vars[vid.index] = { ..old, kind: Integer, }
        }
      Binaries =>
        for k = 0; k < tokens.length(); k = k + 1 {
          let vid = get_or_create_var(model, name_to_id, tokens[k])
          let old = model.vars[vid.index]
          let lb = if old.lb < 0.0 { 0.0 } else { old.lb }
          let ub = if old.ub > 1.0 { 1.0 } else { old.ub }
          model.vars[vid.index] = { ..old, lb, ub, kind: Binary, }
        }
    }
  }

  if !obj_parsed && obj_all_tokens.length() > 0 {
    parse_objective_tokens(obj_all_tokens, model, name_to_id)
  }
  if !has_sense {
    raise ParseError("LP input missing Minimize/Maximize section")
  }
  model
}
