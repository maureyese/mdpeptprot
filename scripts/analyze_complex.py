#!/usr/bin/env python3
"""
analyze_complex.py - peptide-protein MD analysis for one run folder.

    python analyze_complex.py <run_dir> [--skip-frac 0.2] [--cutoff 4.5] [--stride 1]

Reads <run_dir>/00_input_report.json and every <run_dir>/rep*/ (protein.tpr + md.xtc).
Writes <run_dir>/analysis/ and a movie-ready <rep>/md_clean.xtc (+ complex_start.pdb).
Units: Angstrom, ps.
"""
import argparse
import datetime
import glob
import shutil
import subprocess
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore", message="Found no information for attr")

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.lib.distances import capped_distance, distance_array, calc_angles, minimize_vectors


def donor_hydrogen_pairs(ag):
    """(donor heavy-atom indices, hydrogen indices) for every N-H / O-H in ag."""
    don, hyd = [], []
    for h in ag.select_atoms("name H*"):
        for b in h.bonded_atoms:
            if b.name[:1] in ("N", "O"):
                don.append(b.index)
                hyd.append(h.index)
    return np.array(don, dtype=int), np.array(hyd, dtype=int)


def frame_hbonds(pos, don, hyd, acc, dist, angle):
    """Set of (donor_idx, acceptor_idx) pairs that satisfy the distance and angle rules."""
    if len(don) == 0 or len(acc) == 0:
        return set()
    pairs = capped_distance(pos[don], pos[acc], dist, return_distances=False)
    if len(pairs) == 0:
        return set()
    di, ai = pairs[:, 0], pairs[:, 1]
    ang = np.degrees(calc_angles(pos[don[di]], pos[hyd[di]], pos[acc[ai]]))
    ok = ang >= angle
    return set(zip(don[di][ok].tolist(), acc[ai][ok].tolist()))


_NORM = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
         "ASPP": "ASP", "GLUP": "GLU", "LSN": "LYS", "CYS2": "CYS", "CYM": "CYS"}


def load_numbering(path):
    """[(resname, resnum), ...] in file order, from a PDB with the original residue numbers."""
    out, prev = [], None
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21], line[22:27])
            if key != prev:
                out.append((line[17:20].strip(), int(line[22:26])))
                prev = key
    return out


def analyze_rep(rep_dir, out_dir, n_prot, n_pep, a):
    name = os.path.basename(rep_dir.rstrip("/"))
    xtc_in = os.path.join(rep_dir, "md.xtc")
    whole = os.path.join(rep_dir, ".md_whole.xtc")
    gmx = os.environ.get("GMX_BIN") or shutil.which("gmx")
    made_whole = False
    if gmx and not a.no_gmx_whole:       # fast: GROMACS makes the molecules whole for the whole trajectory
        r = subprocess.run([gmx, "trjconv", "-s", os.path.join(rep_dir, "protein.tpr"), "-f", xtc_in,
                            "-o", whole, "-pbc", "whole"], input="0\n", text=True, capture_output=True)
        made_whole = r.returncode == 0 and os.path.isfile(whole)
        if not made_whole:
            print("  warning: gmx trjconv failed, using the slow MDAnalysis unwrap\n" + r.stderr[-300:])
    u = mda.Universe(os.path.join(rep_dir, "protein.tpr"), whole if made_whole else xtc_in)
    n_res = len(u.residues)
    if n_res != n_prot + n_pep:
        sys.exit(f"ERROR {name}: {n_res} residues in the topology, expected {n_prot}+{n_pep}")
    if a.numbering_list is not None:
        nl = a.numbering_list
        if len(nl) != n_res:
            sys.exit(f"ERROR: numbering file has {len(nl)} residues, topology has {n_res}")
        bad = [(i, t, o[0]) for i, (t, o) in enumerate(zip(u.residues.resnames, nl)) if _NORM.get(t, t) != o[0]]
        if bad:
            sys.exit(f"ERROR: residue names disagree with the numbering file, first cases: {bad[:5]}")
        u.residues.resids = np.array([o[1] for o in nl])
        u.residues.resnames = np.array([o[0] for o in nl])
    prot = u.select_atoms(f"resindex 0:{n_prot - 1}")
    pep = u.select_atoms(f"resindex {n_prot}:{n_res - 1}")
    if prot.indices.max() >= pep.indices.min():
        sys.exit(f"ERROR {name}: atom order is not protein-then-peptide")

    prot_heavy = prot.select_atoms("not name H*")
    pep_heavy = pep.select_atoms("not name H*")
    fit_atoms = prot.select_atoms("name CA")
    pep_bb = pep.select_atoms("name N CA C O OT1 OT2")
    don_p, hyd_p = donor_hydrogen_pairs(pep)
    don_P, hyd_P = donor_hydrogen_pairs(prot)
    acc_p = pep.select_atoms("name O* ND1 NE2").indices
    acc_P = prot.select_atoms("name O* ND1 NE2").indices

    chain = np.where(np.arange(u.atoms.n_atoms) < prot.n_atoms, "A", "B")
    try:
        u.add_TopologyAttr("chainIDs", chain)
    except Exception:
        pass

    writer = None if a.no_clean_traj else mda.Writer(os.path.join(rep_dir, "md_clean.xtc"), u.atoms.n_atoms)

    ref_fit = ref_bb = ref_com = None
    rows, hits, hbs, coords = [], [], [], []
    pair_frames = []
    for ts in u.trajectory[::a.stride]:
        # 1) whole chains
        if not made_whole:
            prot.unwrap(compound="fragments")
            pep.unwrap(compound="fragments")
        # 2) peptide into the image nearest to the protein (small vector, so no ambiguity)
        com_p = pep.center_of_mass()
        pp = prot_heavy.positions
        j = int(np.argmin(distance_array(com_p[None, :], pp, box=ts.dimensions)))
        w = minimize_vectors((com_p - pp[j])[None, :].astype(np.float32), ts.dimensions)[0]
        pep.translate(pp[j] + w - com_p)
        # 3) center on the protein and fit on its CA atoms (frame 0 is the reference)
        pf = fit_atoms.positions
        cf = pf.mean(axis=0)
        if ref_fit is None:
            ref_fit = pf - cf
        R, _ = rotation_matrix(pf - cf, ref_fit)
        u.atoms.translate(-cf)
        u.atoms.rotate(R)

        # --- metrics (molecules are whole and adjacent now: no PBC needed) ---
        P = pep_bb.positions.copy()
        com = pep.center_of_mass()
        if ref_bb is None:
            ref_bb, ref_com = P, com
            try:
                u.atoms.write(os.path.join(rep_dir, "complex_start.pdb"))
            except Exception as e:
                print(f"  warning: could not write complex_start.pdb ({e})")
        rmsd = float(np.sqrt(((P - ref_bb) ** 2).sum(axis=1).mean()))
        shift = float(np.linalg.norm(com - ref_com))
        d = distance_array(pep_heavy.positions, prot_heavy.positions)
        hit = np.unique(prot_heavy.resindices[np.where(d < a.cutoff)[1]])
        rr, cc = np.where(d < a.cutoff)
        pair_frames.append(sorted({(int(pep_heavy.resindices[r_]), int(prot_heavy.resindices[c_]))
                                   for r_, c_ in zip(rr, cc)}))
        pos = u.atoms.positions
        keys = (frame_hbonds(pos, don_p, hyd_p, acc_P, a.hb_dist, a.hb_angle)
                | frame_hbonds(pos, don_P, hyd_P, acc_p, a.hb_dist, a.hb_angle))
        rows.append([ts.time, rmsd, shift, float(d.min()), len(hit), len(keys)])
        hits.append(hit)
        hbs.append(keys)
        coords.append(pep_heavy.positions.copy())
        if writer is not None:
            writer.write(u.atoms)
    if writer is not None:
        writer.close()
    if made_whole:
        os.remove(whole)

    ts_arr = np.array(rows)
    # exact H-bond atom pairs per frame (PDB serial numbers = atom index + 1), used to draw the movie
    json.dump([[[int(d) + 1, int(ac) + 1] for d, ac in sorted(k)] for k in hbs],
              open(os.path.join(out_dir, f"{name}_hbond_pairs.json"), "w"))
    # protein residues in contact with the peptide, per frame (labels with the original numbering)
    json.dump([[f"{u.residues[int(ri)].resname}{u.residues[int(ri)].resid}" for ri in h] for h in hits],
              open(os.path.join(out_dir, f"{name}_contact_residues.json"), "w"))
    # residue-residue contacts per frame: [peptide residue, protein residue] (original numbering)
    def _lab(ri):
        r = u.residues[int(ri)]
        return f"{r.resname}{r.resid}"
    json.dump([[[_lab(p_), _lab(q_)] for p_, q_ in fr] for fr in pair_frames],
              open(os.path.join(out_dir, f"{name}_contact_pairs.json"), "w"))
    n = len(ts_arr)
    skip = min(int(a.skip_frac * n), n - 1)
    sel = ts_arr[skip:]
    np.savetxt(os.path.join(out_dir, f"{name}_timeseries.tsv"), ts_arr, delimiter="\t", fmt="%.3f",
               header="time_ps\tpep_rmsd_A\tpep_com_shift_A\tmin_dist_A\tcontact_residues\thbonds", comments="")

    # per-residue contact occupancy
    cnt = np.zeros(n_prot)
    for h in hits[skip:]:
        cnt[h] += 1
    occ = cnt / max(len(hits[skip:]), 1)
    top = []
    with open(os.path.join(out_dir, f"{name}_contacts.tsv"), "w") as f:
        f.write("residue\toccupancy\n")
        for ri in np.argsort(-occ):
            if occ[ri] <= 0:
                break
            r = u.residues[int(ri)]
            f.write(f"{r.resname}{r.resid}\t{occ[ri]:.3f}\n")
            if len(top) < 5:
                top.append([f"{r.resname}{r.resid}", round(float(occ[ri]), 3)])

    # per-pair hydrogen-bond occupancy
    hbc = {}
    for keys in hbs[skip:]:
        for k in keys:
            hbc[k] = hbc.get(k, 0) + 1
    nfr = max(len(hbs[skip:]), 1)

    def alab(i):
        at = u.atoms[i]
        return f"{'pep' if i >= prot.n_atoms else 'prot'}:{at.resname}{at.resid}:{at.name}"

    with open(os.path.join(out_dir, f"{name}_hbonds.tsv"), "w") as f:
        f.write("donor\tacceptor\toccupancy\n")
        for (dn, ac), c in sorted(hbc.items(), key=lambda kv: -kv[1]):
            f.write(f"{alab(dn)}\t{alab(ac)}\t{c / nfr:.3f}\n")

    # peptide flexibility per residue
    X = np.array(coords)[skip:]
    fl = np.sqrt(((X - X.mean(axis=0)) ** 2).sum(axis=2).mean(axis=0))
    rmsf = {}
    for ri in np.unique(pep_heavy.resindices):
        r = u.residues[int(ri)]
        rmsf[f"{r.resname}{r.resid}"] = round(float(fl[pep_heavy.resindices == ri].mean()), 3)

    bound = float((sel[:, 3] < a.cutoff).mean())
    shift_max = float(sel[:, 2].max())
    return {
        "frames": n, "t_end_ps": float(ts_arr[-1, 0]), "analysed_from_ps": float(ts_arr[skip, 0]),
        "pep_rmsd_mean": float(sel[:, 1].mean()), "pep_rmsd_sd": float(sel[:, 1].std()),
        "pep_rmsd_max": float(sel[:, 1].max()),
        "com_shift_mean": float(sel[:, 2].mean()), "com_shift_max": shift_max,
        "contact_residues_mean": float(sel[:, 4].mean()), "hbonds_mean": float(sel[:, 5].mean()),
        "bound_fraction": bound,
        "status": "bound" if (bound >= a.min_bound and shift_max <= a.diss_shift) else "LEFT SITE",
        "top_contacts": top, "peptide_rmsf_A": rmsf,
    }


def _read_tsv(path):
    if not os.path.isfile(path):
        return []
    return [ln.split("\t") for ln in open(path, encoding="utf-8").read().splitlines()[1:] if ln.strip()]


def _consensus(L, title, names, tables, n_top, label_cols):
    """Markdown table: one row per contact / H-bond, occupancy per replica and the mean."""
    occ = {n: {" -> ".join(r[:label_cols]): float(r[label_cols]) for r in tables[n]} for n in names}
    keys = sorted({k for d in occ.values() for k in d}, key=lambda k: -np.mean([occ[n].get(k, 0.0) for n in names]))
    L += ["", title, "", "| " + ("residue" if label_cols == 1 else "donor -> acceptor") + " | " + " | ".join(names) + " | mean |",
          "|---|" + "---|" * (len(names) + 1)]
    for k in keys[:n_top]:
        v = [occ[n].get(k, 0.0) for n in names]
        L.append(f"| {k} | " + " | ".join(f"{x:.2f}" for x in v) + f" | {np.mean(v):.2f} |")


def write_report(out_dir, run_name, rep, summary, overall, a):
    names = list(summary)
    pep = rep.get("peptide_sequence", "peptide")
    ov = overall
    t_end = max(s["t_end_ps"] for s in summary.values()) / 1000
    t_from = min(s["analysed_from_ps"] for s in summary.values()) / 1000
    pm = lambda k: f"{ov[k]['mean']:.2f} ± {ov[k]['sd_between_replicas']:.2f}"
    n_bound = len(names) - ov["n_left_site"]
    L = [f"# MD analysis report: {run_name}", "",
         f"- Peptide **{pep}** | {len(names)} replicas | {t_end:.2f} ns simulated each | analysed from {t_from:.2f} ns (first {a.skip_frac:.0%} dropped as settling time)",
         f"- Contact = any heavy atoms closer than {a.cutoff} Å. Hydrogen bond = donor-acceptor ≤ {a.hb_dist} Å and D-H···A angle ≥ {a.hb_angle:.0f}°.",
         f"- 'bound' = in contact in ≥ {a.min_bound:.0%} of frames and peptide centre of mass moved ≤ {a.diss_shift} Å from its start pose.",
         "- Residue numbers follow the original PDB numbering (peptide residues as numbered in the input).",
         "", "## Per replica", "",
         "| replica | frames | RMSD mean (Å) | RMSD max (Å) | COM shift max (Å) | contact residues | H-bonds | bound fraction | status |",
         "|---|---|---|---|---|---|---|---|---|"]
    for n, s in summary.items():
        L.append(f"| {n} | {s['frames']} | {s['pep_rmsd_mean']:.2f} | {s['pep_rmsd_max']:.2f} | {s['com_shift_max']:.2f} | "
                 f"{s['contact_residues_mean']:.1f} | {s['hbonds_mean']:.2f} | {s['bound_fraction']:.2f} | {s['status']} |")
    L.append(f"| **mean ± sd** | | {pm('pep_rmsd_mean')} | | {pm('com_shift_max')} | {pm('contact_residues_mean')} | "
             f"{pm('hbonds_mean')} | {pm('bound_fraction')} | {n_bound}/{len(names)} bound |")

    _consensus(L, "## Contacts held (fraction of analysed frames; top 15 by mean)", names,
               {n: _read_tsv(os.path.join(out_dir, f"{n}_contacts.tsv")) for n in names}, 15, 1)
    _consensus(L, "## Hydrogen bonds (fraction of analysed frames; top 10 by mean)", names,
               {n: _read_tsv(os.path.join(out_dir, f"{n}_hbonds.tsv")) for n in names}, 10, 2)

    resn = list(next(iter(summary.values()))["peptide_rmsf_A"])
    L += ["", "## Peptide flexibility (RMSF, Å)", "", "| replica | " + " | ".join(resn) + " |", "|---|" + "---|" * len(resn)]
    for n, s in summary.items():
        L.append(f"| {n} | " + " | ".join(f"{s['peptide_rmsf_A'][r]:.2f}" for r in resn) + " |")

    eq = os.path.join(os.path.dirname(out_dir), "05_equilibration_report.json")
    if os.path.isfile(eq):
        try:
            L += ["", "## Equilibration check", ""] + [f"- {k}: {v}" for k, v in json.load(open(eq)).items()
                                                      if isinstance(v, (int, float, str, bool))]
        except Exception:
            pass
    job = os.path.join(os.path.dirname(os.path.dirname(out_dir)), f"{run_name}.job.yaml")
    if os.path.isfile(job):
        L += ["", "## Job file", "", "```yaml", open(job, encoding="utf-8").read().rstrip(), "```"]
    L += ["", f"_Generated {datetime.datetime.now():%Y-%m-%d %H:%M}. Per-frame tables: rep*_timeseries.tsv, contacts: rep*_contacts.tsv, hydrogen bonds: rep*_hbonds.tsv._"]
    with open(os.path.join(out_dir, "report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")

    cols = [("frames", "frames"), ("pep_rmsd_mean", "rmsd_mean_A"), ("pep_rmsd_max", "rmsd_max_A"),
            ("com_shift_max", "com_shift_max_A"), ("contact_residues_mean", "contact_residues"),
            ("hbonds_mean", "hbonds"), ("bound_fraction", "bound_fraction")]
    with open(os.path.join(out_dir, "summary.tsv"), "w", encoding="utf-8") as f:
        f.write("\t".join(["run", "peptide", "replica"] + [c[1] for c in cols] + ["status"]) + "\n")
        for n, s in summary.items():
            f.write("\t".join([run_name, pep, n] + [f"{s[c[0]]:.3f}" for c in cols] + [s["status"]]) + "\n")
        f.write("\t".join([run_name, pep, "mean"] + [f"{np.mean([s[c[0]] for s in summary.values()]):.3f}" for c in cols]
                         + [f"{n_bound}/{len(names)} bound"]) + "\n")
    print(f"[analyze] report: {os.path.join(out_dir, 'report.md')}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--skip-frac", type=float, default=0.2, help="fraction of frames dropped as settling time")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--cutoff", type=float, default=4.5, help="contact distance, heavy atoms (A)")
    ap.add_argument("--hb-dist", type=float, default=3.5)
    ap.add_argument("--hb-angle", type=float, default=150.0)
    ap.add_argument("--min-bound", type=float, default=0.9, help="min fraction of frames in contact")
    ap.add_argument("--diss-shift", type=float, default=10.0, help="peptide COM shift (A) that counts as leaving")
    ap.add_argument("--no-clean-traj", action="store_true")
    ap.add_argument("--no-gmx-whole", action="store_true", help="use the (slow) MDAnalysis unwrap")
    ap.add_argument("--numbering", default=None, help="PDB with the original residue numbers (default: <run_dir>/input_topology.pdb)")
    a = ap.parse_args()

    rep = json.load(open(os.path.join(a.run_dir, "00_input_report.json")))
    n_prot, n_pep = int(rep["protein_residues"]), int(rep["peptide_residues"])
    out_dir = os.path.join(a.run_dir, "analysis")
    num_path = a.numbering or os.path.join(a.run_dir, "input_topology.pdb")
    a.numbering_list = load_numbering(num_path) if os.path.isfile(num_path) else None
    print(f"[analyze] residue numbering: {num_path if a.numbering_list else 'topology numbering (no numbering file found)'}")
    os.makedirs(out_dir, exist_ok=True)

    reps = sorted(glob.glob(os.path.join(a.run_dir, "rep*", "protein.tpr")))
    if not reps:
        sys.exit("ERROR: no rep*/protein.tpr found")
    summary = {}
    for tpr in reps:
        rd = os.path.dirname(tpr)
        name = os.path.basename(rd)
        print(f"[analyze] {name} ...")
        summary[name] = analyze_rep(rd, out_dir, n_prot, n_pep, a)

    keys = ["pep_rmsd_mean", "com_shift_max", "contact_residues_mean", "hbonds_mean", "bound_fraction"]
    overall = {k: {"mean": float(np.mean([s[k] for s in summary.values()])),
                   "sd_between_replicas": float(np.std([s[k] for s in summary.values()]))} for k in keys}
    overall["n_left_site"] = sum(s["status"] != "bound" for s in summary.values())
    write_report(out_dir, os.path.basename(os.path.abspath(a.run_dir)), rep, summary, overall, a)
    json.dump({"peptide": rep.get("peptide_sequence"), "replicas": summary, "overall": overall},
              open(os.path.join(out_dir, "summary.json"), "w"), indent=2)

    print(f"\nPeptide {rep.get('peptide_sequence')}   (first {a.skip_frac:.0%} of frames dropped)")
    print(f"{'rep':<6}{'frames':>7}{'rmsd':>8}{'rmsd_max':>10}{'shift_max':>11}{'contacts':>10}{'hbonds':>8}{'bound':>7}  status")
    for n, s in summary.items():
        print(f"{n:<6}{s['frames']:>7}{s['pep_rmsd_mean']:>8.2f}{s['pep_rmsd_max']:>10.2f}{s['com_shift_max']:>11.2f}"
              f"{s['contact_residues_mean']:>10.1f}{s['hbonds_mean']:>8.2f}{s['bound_fraction']:>7.2f}  {s['status']}")
    print(f"\nWrote {out_dir}/ (summary.json, rep*_timeseries.tsv, rep*_contacts.tsv, rep*_hbonds.tsv)")


if __name__ == "__main__":
    main()
