///|
/// Industry-standard MPS (.mps) format parser (compatible with Netlib LP & MIPLIB).
/// Supports NAME, ROWS (N/L/G/E), COLUMNS (including 'MARKER' 'INTORG'/'INTEND'),
/// RHS, BOUNDS (LO/UP/FX/FR/BV/LI/UI/MI/PL), and ENDATA.

///|
priv enum MpsSection {
  None
  Rows
  Columns
  Rhs
  Ranges
  Bounds
}

///|
fn split_whitespace(s : String) -> Array[String] {
  let out : Array[String] = []
  let mut buf = StringBuilder()
  for ch in s {
    if ch == ' ' || ch == '\t' || ch == '\r' {
      let tok = buf.to_string()
      if tok.length() > 0 {
        out.push(tok)
        buf = StringBuilder()
      }
    } else {
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
fn is_mps_comment(line : String) -> Bool {
  for ch in line {
    return ch == '*'
  }
  false
}

///|
/// Parse an industry-standard free-format MPS string into a `Model` (default Minimize).
pub fn parse_mps(
  text : String,
  sense~ : Sense = Minimize,
) -> Model raise ParseError {
  let lines = split_lines(text)
  let model = Model::new(sense)
  let mut section : MpsSection = None
  let mut obj_row_name = ""
  let row_name_to_idx : Map[String, Int32] = {}
  let row_names : Array[String] = []
  let row_rels : Array[Relation] = []
  let row_terms : Array[Array[Term]] = []
  let row_rhs : Array[Double] = []
  let row_range : Array[Double?] = []

  let var_name_to_id : Map[String, VarId] = {}
  let mut in_int_block = false
  let mut saw_rows = false

  for li = 0; li < lines.length(); li = li + 1 {
    let raw = lines[li]
    if raw.length() == 0 || is_mps_comment(raw) {
      continue
    }
    let fields = split_whitespace(raw)
    if fields.length() == 0 {
      continue
    }
    let head = to_lower_ascii(fields[0])

    if head == "name" {
      continue
    }
    if head == "rows" {
      section = Rows
      saw_rows = true
      continue
    }
    if head == "columns" {
      section = Columns
      continue
    }
    if head == "rhs" {
      section = Rhs
      continue
    }
    if head == "ranges" {
      section = Ranges
      continue
    }
    if head == "bounds" {
      section = Bounds
      continue
    }
    if head == "objsense" {
      continue
    }
    if head == "endata" {
      break
    }

    match section {
      None => raise ParseError("MPS content before ROWS section: " + raw)
      Rows => {
        if fields.length() < 2 {
          raise ParseError("Invalid MPS ROWS entry: " + raw)
        }
        let rtype = to_lower_ascii(fields[0])
        let rname = fields[1]
        if rtype == "n" {
          if obj_row_name.length() == 0 {
            obj_row_name = rname
          }
        } else {
          let rel = if rtype == "l" {
            Le
          } else if rtype == "g" {
            Ge
          } else if rtype == "e" {
            Eq
          } else {
            raise ParseError("Unknown MPS row type: " + fields[0])
          }
          let ridx = row_names.length()
          row_name_to_idx.set(rname, ridx)
          row_names.push(rname)
          row_rels.push(rel)
          row_terms.push([])
          row_rhs.push(0.0)
          row_range.push(None)
        }
      }
      Columns => {
        // Check for integer block markers: <name> 'MARKER' 'INTORG' / 'INTEND'
        if fields.length() >= 3 && to_lower_ascii(fields[1]) == "'marker'" {
          let mkind = to_lower_ascii(fields[2])
          if mkind == "'intorg'" {
            in_int_block = true
          } else if mkind == "'intend'" {
            in_int_block = false
          }
          continue
        }
        if fields.length() < 3 {
          raise ParseError("Invalid MPS COLUMNS entry: " + raw)
        }
        let vname = fields[0]
        let vid = match var_name_to_id.get(vname) {
          Some(id) => id
          None => {
            let kind : VarKind = if in_int_block { Integer } else { Continuous }
            let id = model.add_var(vname, kind~)
            var_name_to_id.set(vname, id)
            id
          }
        }
        let mut p = 1
        while p + 1 < fields.length() {
          let rname = fields[p]
          let coeff = match try_parse_double(fields[p + 1]) {
            Some(v) => v
            None =>
              raise ParseError("Invalid MPS coefficient: " + fields[p + 1])
          }
          if rname == obj_row_name {
            let old = model.vars[vid.index]
            model.vars[vid.index] = { ..old, obj: old.obj + coeff, }
          } else {
            match row_name_to_idx.get(rname) {
              Some(ridx) => row_terms[ridx].push({ var_id: vid, coeff, })
              None => ()
            }
          }
          p = p + 2
        }
      }
      Rhs => {
        if fields.length() < 3 {
          raise ParseError("Invalid MPS RHS entry: " + raw)
        }
        let mut p = 1
        while p + 1 < fields.length() {
          let rname = fields[p]
          let val = match try_parse_double(fields[p + 1]) {
            Some(v) => v
            None => raise ParseError("Invalid MPS RHS number: " + fields[p + 1])
          }
          match row_name_to_idx.get(rname) {
            Some(ridx) => row_rhs[ridx] = val
            None => ()
          }
          p = p + 2
        }
      }
      Ranges => {
        if fields.length() < 3 {
          raise ParseError("Invalid MPS RANGES entry: " + raw)
        }
        let mut p = 1
        while p + 1 < fields.length() {
          let rname = fields[p]
          let val = match try_parse_double(fields[p + 1]) {
            Some(v) => v
            None =>
              raise ParseError("Invalid MPS RANGES number: " + fields[p + 1])
          }
          match row_name_to_idx.get(rname) {
            Some(ridx) => row_range[ridx] = Some(val)
            None => ()
          }
          p = p + 2
        }
      }
      Bounds => {
        if fields.length() < 3 {
          raise ParseError("Invalid MPS BOUNDS entry: " + raw)
        }
        let btype = to_lower_ascii(fields[0])
        let vname = fields[2]
        let vid = match var_name_to_id.get(vname) {
          Some(id) => id
          None => {
            let id = model.add_continuous_var(vname)
            var_name_to_id.set(vname, id)
            id
          }
        }
        let old = model.vars[vid.index]
        if btype == "fr" {
          model.set_var_bounds(vid, -1.0e300, 1.0e300)
        } else if btype == "mi" {
          model.set_var_bounds(vid, -1.0e300, old.ub)
        } else if btype == "pl" {
          model.set_var_bounds(vid, old.lb, 1.0e300)
        } else if btype == "bv" {
          model.vars[vid.index] = { ..old, lb: 0.0, ub: 1.0, kind: Binary, }
        } else {
          if fields.length() < 4 {
            raise ParseError("Missing bound value in MPS BOUNDS: " + raw)
          }
          let bval = match try_parse_double(fields[3]) {
            Some(v) => v
            None => raise ParseError("Invalid MPS bound value: " + fields[3])
          }
          if btype == "lo" {
            model.set_var_bounds(vid, bval, old.ub)
          } else if btype == "up" {
            model.set_var_bounds(vid, old.lb, bval)
          } else if btype == "fx" {
            model.set_var_bounds(vid, bval, bval)
          } else if btype == "li" {
            model.vars[vid.index] = { ..old, lb: bval, kind: Integer, }
          } else if btype == "ui" {
            model.vars[vid.index] = { ..old, ub: bval, kind: Integer, }
          } else {
            raise ParseError("Unknown MPS bound type: " + fields[0])
          }
        }
      }
    }
  }

  if !saw_rows {
    raise ParseError("MPS input missing ROWS section")
  }

  for i = 0; i < row_names.length(); i = i + 1 {
    match row_range[i] {
      None =>
        model.add_constraint(row_names[i], row_terms[i], row_rels[i], row_rhs[i])
      Some(rval) => {
        let abs_r = rval.abs()
        match row_rels[i] {
          Ge => {
            // b <= a^T x <= b + |r|
            model.add_constraint(row_names[i], row_terms[i], Ge, row_rhs[i])
            model.add_constraint(
              row_names[i] + "_ub",
              row_terms[i],
              Le,
              row_rhs[i] + abs_r,
            )
          }
          Le => {
            // b - |r| <= a^T x <= b
            model.add_constraint(row_names[i], row_terms[i], Le, row_rhs[i])
            model.add_constraint(
              row_names[i] + "_lb",
              row_terms[i],
              Ge,
              row_rhs[i] - abs_r,
            )
          }
          Eq =>
            if rval >= 0.0 {
              model.add_constraint(row_names[i], row_terms[i], Ge, row_rhs[i])
              model.add_constraint(
                row_names[i] + "_ub",
                row_terms[i],
                Le,
                row_rhs[i] + rval,
              )
            } else {
              model.add_constraint(
                row_names[i] + "_lb",
                row_terms[i],
                Ge,
                row_rhs[i] + rval,
              )
              model.add_constraint(row_names[i], row_terms[i], Le, row_rhs[i])
            }
        }
      }
    }
  }
  model
}

///|
/// Export a Model to industry-standard free-format Netlib/MIPLIB MPS (.mps) text.
pub fn Model::to_mps_string(self : Model, name~ : String = "MOONLP") -> String {
  let buf = StringBuilder()
  buf.write_string("NAME          ")
  buf.write_string(name)
  buf.write_string("\nROWS\n")
  buf.write_string(" N  OBJ\n")

  let row_names : Array[String] = []
  for i = 0; i < self.constraints.length(); i = i + 1 {
    let c = self.constraints[i]
    let rname = if c.name.length() > 0 { c.name } else { "R" + i.to_string() }
    row_names.push(rname)
    let rtype = match c.relation {
      Le => " L  "
      Ge => " G  "
      Eq => " E  "
    }
    buf.write_string(rtype)
    buf.write_string(rname)
    buf.write_string("\n")
  }

  // Build column-major entries: for each variable j, list of (row_name, coeff)
  let col_entries : Array[Array[(String, Double)]] = []
  let sense_sign = match self.sense {
    Minimize => 1.0
    Maximize => -1.0
  }
  for j = 0; j < self.vars.length(); j = j + 1 {
    let entries : Array[(String, Double)] = []
    if self.vars[j].obj != 0.0 {
      entries.push(("OBJ", self.vars[j].obj * sense_sign))
    }
    col_entries.push(entries)
  }
  for i = 0; i < self.constraints.length(); i = i + 1 {
    let c = self.constraints[i]
    let rname = row_names[i]
    for k = 0; k < c.terms.length(); k = k + 1 {
      let t = c.terms[k]
      if t.coeff != 0.0 {
        col_entries[t.var_id.index].push((rname, t.coeff))
      }
    }
  }

  buf.write_string("COLUMNS\n")
  let mut in_int = false
  let mut mark_id = 0
  for j = 0; j < self.vars.length(); j = j + 1 {
    let v = self.vars[j]
    let need_int = match v.kind {
      Continuous => false
      Integer | Binary => true
    }
    if need_int && !in_int {
      buf.write_string("    MARK")
      buf.write_string(mark_id.to_string())
      buf.write_string("     'MARKER'                 'INTORG'\n")
      mark_id = mark_id + 1
      in_int = true
    } else if !need_int && in_int {
      buf.write_string("    MARK")
      buf.write_string(mark_id.to_string())
      buf.write_string("     'MARKER'                 'INTEND'\n")
      mark_id = mark_id + 1
      in_int = false
    }
    let entries = col_entries[j]
    if entries.length() == 0 {
      buf.write_string("    ")
      buf.write_string(v.name)
      buf.write_string("  OBJ  0.0\n")
    } else {
      for e = 0; e < entries.length(); e = e + 1 {
        let (rname, coeff) = entries[e]
        buf.write_string("    ")
        buf.write_string(v.name)
        buf.write_string("  ")
        buf.write_string(rname)
        buf.write_string("  ")
        buf.write_string(coeff.to_string())
        buf.write_string("\n")
      }
    }
  }
  if in_int {
    buf.write_string("    MARK")
    buf.write_string(mark_id.to_string())
    buf.write_string("     'MARKER'                 'INTEND'\n")
  }

  buf.write_string("RHS\n")
  for i = 0; i < self.constraints.length(); i = i + 1 {
    if self.constraints[i].rhs != 0.0 {
      buf.write_string("    RHS1  ")
      buf.write_string(row_names[i])
      buf.write_string("  ")
      buf.write_string(self.constraints[i].rhs.to_string())
      buf.write_string("\n")
    }
  }

  buf.write_string("BOUNDS\n")
  for j = 0; j < self.vars.length(); j = j + 1 {
    let v = self.vars[j]
    let has_lb = v.lb > -1.0e200
    let has_ub = v.ub < 1.0e200
    match v.kind {
      Binary => {
        buf.write_string(" BV BND1  ")
        buf.write_string(v.name)
        buf.write_string("\n")
      }
      _ =>
        if !has_lb && !has_ub {
          buf.write_string(" FR BND1  ")
          buf.write_string(v.name)
          buf.write_string("\n")
        } else if has_lb && has_ub && (v.ub - v.lb).abs() <= 1.0e-12 {
          buf.write_string(" FX BND1  ")
          buf.write_string(v.name)
          buf.write_string("  ")
          buf.write_string(v.lb.to_string())
          buf.write_string("\n")
        } else {
          if !has_lb {
            buf.write_string(" MI BND1  ")
            buf.write_string(v.name)
            buf.write_string("\n")
          } else if v.lb != 0.0 {
            buf.write_string(" LO BND1  ")
            buf.write_string(v.name)
            buf.write_string("  ")
            buf.write_string(v.lb.to_string())
            buf.write_string("\n")
          }
          if has_ub {
            buf.write_string(" UP BND1  ")
            buf.write_string(v.name)
            buf.write_string("  ")
            buf.write_string(v.ub.to_string())
            buf.write_string("\n")
          }
        }
    }
  }

  buf.write_string("ENDATA\n")
  buf.to_string()
}
