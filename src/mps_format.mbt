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
    model.add_constraint(row_names[i], row_terms[i], row_rels[i], row_rhs[i])
  }
  model
}
