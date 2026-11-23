///|
/// CPLEX .lp format exporter for industry-standard interoperability
/// (compatible with HiGHS, Gurobi, GLPK, SCIP, and COIN-OR CBC).

///|
fn format_linear_terms(
  model : Model,
  terms : Array[Term],
) -> String {
  if terms.length() == 0 {
    return "0"
  }
  let buf = StringBuilder::new()
  for i = 0; i < terms.length(); i = i + 1 {
    let t = terms[i]
    let c = t.coeff
    if c == 0.0 {
      continue
    }
    let name = model.vars[t.var_id.index].name
    if i == 0 {
      if c < 0.0 {
        buf.write_string("- ")
        if (-c) != 1.0 {
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
      if (-c) != 1.0 {
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
  buf.to_string()
}

///|
/// Export a Model to standard CPLEX `.lp` text format.
pub fn Model::to_lp_string(self : Model) -> String {
  let buf = StringBuilder::new()
  match self.sense {
    Minimize => buf.write_string("Minimize\n")
    Maximize => buf.write_string("Maximize\n")
  }
  let obj_terms : Array[Term] = []
  for j = 0; j < self.vars.length(); j = j + 1 {
    if self.vars[j].obj != 0.0 {
      obj_terms.push({ var_id: { index: j }, coeff: self.vars[j].obj })
    }
  }
  buf.write_string("  obj: ")
  buf.write_string(format_linear_terms(self, obj_terms))
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
