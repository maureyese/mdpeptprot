#!/usr/bin/env python3
"""
protonate_from_pqr.py INPUT.pdb PDB2PQR.pqr OUTPUT.pdb AMINOACIDS.rtp REPORT.tsv PH

Reads the protonation states chosen by PDB2PQR/PROPKA (from residue names AND
from which hydrogens it placed) and writes them into the residue names CHARMM36
understands: HSD/HSE/HSP, ASPP, GLUP, LSN. Coordinates and atom names of INPUT
are untouched.
"""
import sys
from collections import Counter

STD = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
       "MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
AMBER2CHARMM = {"HID": "HSD", "HIE": "HSE", "HIP": "HSP",
                "ASH": "ASPP", "GLH": "GLUP", "LYN": "LSN"}
KNOWN = STD | set(AMBER2CHARMM) | {"CYX", "CYM"}

def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

def norm(n):                       # strip AMBER terminal prefixes (NILE -> ILE)
    n = n.upper()
    if len(n) == 4 and n[0] in "NC" and n[1:] in KNOWN:
        return n[1:]
    return n

def assign(old, st, a):
    """Return (new_resname, warning_or_None) for one residue."""
    if st in AMBER2CHARMM:
        return AMBER2CHARMM[st], None
    if old == "HIS":
        if "HD1" in a and "HE2" in a: return "HSP", None
        if "HE2" in a:                return "HSE", None
        if "HD1" in a:                return "HSD", None
        return "HSD", "no ring hydrogen found; using HSD"
    if old == "ASP" and "HD2" in a:
        return "ASPP", None
    if old == "GLU" and "HE2" in a:
        return "GLUP", None
    if old == "LYS" and "HZ1" in a and "HZ3" not in a:
        return "LSN", None
    if old == "CYS" and st == "CYM":
        return "CYM", None
    if old == "TYR" and "HH" not in a:
        fail("tyrosine predicted deprotonated by PROPKA, but the force field has no "
             "deprotonated TYR residue; cannot build this system at this pH")
    return old, None

def main():
    inp, pqr, out, rtp, report, ph = sys.argv[1:7]

    states, atoms = {}, {}          # (chain, resid) -> name / set of atom names
    for l in open(pqr):
        if l.startswith("ATOM"):
            f = l.split()
            key = (f[4], f[5])
            states[key] = norm(f[3])
            atoms.setdefault(key, set()).add(f[2])

    rtp_names = set()
    for l in open(rtp):
        s = l.strip()
        if s.startswith("[") and s.endswith("]"):
            rtp_names.add(s[1:-1].strip())

    res_order, res_old, res_new, warns = [], {}, {}, []
    out_lines = []
    for l in open(inp).read().splitlines():
        if l.startswith("ATOM  "):
            key = (l[21], l[22:27].strip())
            if key not in res_old:
                old = l[17:20].strip()
                new, w = old, None
                if old in ("HIS", "ASP", "GLU", "LYS", "CYS", "TYR"):
                    new, w = assign(old, states.get(key), atoms.get(key, set()))
                if w:
                    warns.append(f"{key}: {w}")
                res_order.append(key)
                res_old[key], res_new[key] = old, new
            l = l[:17] + f"{res_new[key]:<4}" + l[21:]
        out_lines.append(l)

    n_unmatched = sum(1 for k in res_order if k not in states)
    if n_unmatched > 0.05 * len(res_order):
        fail(f"could not match {n_unmatched}/{len(res_order)} residues between the input and the pdb2pqr output")

    missing = sorted({n for n in res_new.values() if n not in rtp_names})
    if missing:
        fail(f"residue names {missing} are not defined in the force field's aminoacids.rtp")

    open(out, "w").write("\n".join(out_lines) + "\n")

    with open(report, "w") as f:
        f.write("chain\tresid\toriginal\tassigned\n")
        for k in res_order:
            if res_old[k] == "HIS" or res_old[k] != res_new[k]:
                f.write(f"{k[0]}\t{k[1]}\t{res_old[k]}\t{res_new[k]}\n")

    for w in warns[:10]:
        print("WARNING:", w)
    if len(warns) > 10:
        print(f"WARNING: ... and {len(warns) - 10} more")
    c = Counter(res_new[k] for k in res_order
                if res_new[k] in ("HSD", "HSE", "HSP", "ASPP", "GLUP", "LSN"))
    print(f"Protonation at pH {ph}: " + (", ".join(f"{n} x{c[n]}" for n in sorted(c)) or "default states only"))

main()
