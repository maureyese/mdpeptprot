#!/usr/bin/env bash
# =============================================================================
# 00_prepare_input.sh
# Validates the input complex (1 protein chain + 1 short peptide chain) and
# writes input_topology.pdb with the protein as chain A and the peptide as
# chain B. Fails early, with a clear message, on anything the later stages
# cannot handle.
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/common.sh"

STAGE="00_prepare_input"
stage_start
log "=== Stage: $STAGE ==="

require_file "$INPUT_PDB"

RAW_COPY="$RUN_DIR/00_raw_input.pdb"
cp "$INPUT_PDB" "$RAW_COPY"
cd "$RUN_DIR"

python3 - "$RAW_COPY" "$RUN_DIR/input_topology.pdb" "$RUN_DIR/00_input_report.json" \
          "${PEPTIDE_MIN_LEN:-3}" "${PEPTIDE_MAX_LEN:-6}" <<'PYEOF'
import sys, json, math

src, dst, report_path = sys.argv[1], sys.argv[2], sys.argv[3]
min_len, max_len = int(sys.argv[4]), int(sys.argv[5])

AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
      "MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}

def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

lines = open(src).read().splitlines()

n_models = sum(1 for l in lines if l.startswith("MODEL"))
if n_models > 1:
    fail(f"input has {n_models} models; provide a single model")

chains, hetero = {}, set()
for l in lines:
    if l.startswith("ENDMDL"):
        break
    if not l.startswith(("ATOM  ", "HETATM")):
        continue
    res = l[17:20].strip()
    if res not in AA:
        hetero.add((l[21], res))
        continue
    if l[16] not in (" ", "A"):          # first alternate location only
        continue
    chains.setdefault(l[21], []).append(l)

if hetero:
    fail(f"non-amino-acid residues found (chain, name): {sorted(hetero)}. "
         "Water, ligands, sugars and cofactors are not supported: remove them first.")

def nres(ls): return len({l[22:27] for l in ls})
def xyz(l): return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
def heavy(l):
    el = l[76:78].strip() or l[12:16].strip().lstrip("0123456789")[:1]
    return el.upper() != "H"

if len(chains) != 2:
    fail(f"expected exactly 2 chains (protein + peptide), found {len(chains)}: {sorted(chains)}")

pep, prot = sorted(chains, key=lambda c: nres(chains[c]))   # shorter chain = peptide
pep_l, prot_l = chains[pep], chains[prot]
n_pep, n_prot = nres(pep_l), nres(prot_l)
if not (min_len <= n_pep <= max_len):
    fail(f"peptide chain {pep} has {n_pep} residues; expected {min_len}-{max_len}")
if n_prot < 50:
    fail(f"protein chain {prot} has only {n_prot} residues; expected a full protein")

def residues_in_order(ls):
    order, atoms = [], {}
    for l in ls:
        k = l[22:27]
        if k not in atoms:
            atoms[k] = {}
            order.append(k)
        atoms[k][l[12:16].strip()] = xyz(l)
    return order, atoms

breaks = {}
for name, ls in (("protein", prot_l), ("peptide", pep_l)):
    order, atoms = residues_in_order(ls)
    bad = [k.strip() for k in order if not {"N", "CA", "C"} <= set(atoms[k])]
    if bad:
        fail(f"{name} chain: residues missing backbone atoms N/CA/C: {bad[:10]}")
    brk = [(a.strip(), b.strip()) for a, b in zip(order, order[1:])
           if math.dist(atoms[a]["C"], atoms[b]["N"]) > 2.0]
    breaks[name] = brk
    if brk:
        print(f"WARNING: {name} chain has {len(brk)} chain break(s) (residue pairs): {brk[:10]}")

seq, seen = [], None
for l in pep_l:
    if l[22:27] != seen:
        seq.append(l[17:20])
        seen = l[22:27]
sequence = "-".join(seq)

pc = [xyz(l) for l in prot_l if heavy(l)]
n_contacts, dmin2 = 0, 1e9
for l in pep_l:
    if not heavy(l):
        continue
    a = xyz(l)
    for b in pc:
        d2 = (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2
        if d2 < 20.25:                    # 4.5 A
            n_contacts += 1
        if d2 < dmin2:
            dmin2 = d2
dmin = math.sqrt(dmin2)
if dmin < 1.5:
    fail(f"severe clash: peptide and protein heavy atoms are only {dmin:.2f} A apart")
if n_contacts == 0:
    fail(f"peptide is not in contact with the protein (closest heavy atoms {dmin:.2f} A apart)")
if dmin < 2.0:
    print(f"WARNING: close contact between peptide and protein ({dmin:.2f} A)")

with open(dst, "w") as f:
    for new, ls in (("A", prot_l), ("B", pep_l)):
        for l in ls:
            f.write("ATOM  " + l[6:21] + new + l[22:] + "\n")
        f.write("TER\n")
    f.write("END\n")

json.dump({"protein_chain_original": prot, "peptide_chain_original": pep,
           "protein_residues": n_prot, "peptide_residues": n_pep,
           "peptide_sequence": sequence, "contacts_4.5A": n_contacts,
           "min_distance_A": round(dmin, 2), "chain_breaks": breaks},
          open(report_path, "w"), indent=2)

print(f"Protein chain {prot} ({n_prot} res) + peptide chain {pep} ({sequence}); "
      f"{n_contacts} contacts < 4.5 A, min distance {dmin:.2f} A")
print("Relabelled in output: protein = chain A, peptide = chain B")
PYEOF

require_file "$RUN_DIR/input_topology.pdb"
log "Input prepared: $RUN_DIR/input_topology.pdb"
mark_done "$STAGE"
