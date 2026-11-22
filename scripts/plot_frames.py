#!/usr/bin/env python3
"""plot_frames.py <run_dir> <rep_number> <movie_dir>
Per rendered 3D frame writes lbl_NNNN.png (3D frame + title + frosted residue labels + contacts table)
and p_NNNN.png (plot panel)."""
import json
import math
import os
import re
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from PIL import Image, ImageDraw, ImageFilter, ImageFont
import MDAnalysis as mda

warnings.filterwarnings("ignore")
BLUR, WHITE_ALPHA = 4, 0.60                 # frosted-glass look of labels and table
LABEL_FONT, TABLE_FONT = 14, 12             # text sizes (px)

run_dir, rep, mdir = sys.argv[1], int(sys.argv[2]), sys.argv[3]
rd = os.path.join(run_dir, f"rep{rep}")
an = os.path.join(run_dir, "analysis")
pep = json.load(open(os.path.join(run_dir, "00_input_report.json"))).get("peptide_sequence", "peptide")
_AA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
       "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
       "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
if "-" in pep:                                   # three-letter names -> one-letter code
    pep = "".join(_AA.get(x, "X") for x in pep.split("-"))

ts = np.genfromtxt(os.path.join(an, f"rep{rep}_timeseries.tsv"), skip_header=1)
states = [int(x) for x in open(os.path.join(mdir, "states.txt"))]
lab = json.load(open(os.path.join(mdir, "labels.json")))
hbpairs = json.load(open(os.path.join(an, f"rep{rep}_hbond_pairs.json")))
contacts = json.load(open(os.path.join(an, f"rep{rep}_contact_pairs.json")))   # per frame: [[pep res, prot res], ...]


def nice(label):
    """ILE1 -> Ile1"""
    m = re.match(r"([A-Za-z]+)(-?\d+)", label)
    return m.group(1).capitalize() + m.group(2) if m else label


def entry(pair):
    return f"Pept-{nice(pair[0])}:Prot-{nice(pair[1])}"


# atom serial -> (chain, residue label), to know which residue pairs are H-bonded
serial = {}
for line in open(os.path.join(rd, "complex_start.pdb")):
    if line.startswith(("ATOM", "HETATM")):
        serial[int(line[6:11])] = (line[21], line[17:20].strip() + line[22:26].strip())


def hb_pair(pr):
    b = [serial[x][1] for x in pr if serial.get(x, ("", ""))[0] == "B"]
    a_ = [serial[x][1] for x in pr if serial.get(x, ("", ""))[0] == "A"]
    return (b[0], a_[0]) if b and a_ else None


hb_res = [{p for p in (hb_pair(pr) for pr in fr) if p} for fr in hbpairs]

u = mda.Universe(os.path.join(rd, "complex_start.pdb"), os.path.join(rd, "md_clean.xtc"))
times = np.array([t.time for t in u.trajectory])
tsel = times[np.array(states) - 1]
idx = np.abs(ts[:, 0][None, :] - tsel[:, None]).argmin(axis=1)
tn = ts[:, 0] / 1000.0

# ---------------- plot panels ----------------
plt.style.use("dark_background")
panels = [("Peptide RMSD (Å)", 1, "#ff6ec7"), ("Min. distance to protein (Å)", 3, "#00e5ff"),
          ("Contact residues", 4, "#b388ff"), ("Hydrogen bonds", 5, "#ffe600")]
fig, axes = plt.subplots(4, 1, figsize=(6.4, 7.2), dpi=100, sharex=True)
artists = []
for ax, (label, col, color) in zip(axes, panels):
    y = ts[:, col]
    ax.plot(tn, y, color=color, alpha=0.18, lw=1)
    (line,) = ax.plot([], [], color=color, lw=1.6)
    (dot,) = ax.plot([], [], "o", color=color, ms=5)
    pad = 0.1 * (y.max() - y.min()) + 1e-6
    ax.set_xlim(tn[0], tn[-1])
    ax.set_ylim(y.min() - pad, y.max() + pad)
    ax.set_ylabel(label, fontsize=9)
    ax.tick_params(labelsize=8)
    artists.append((line, dot, y))
axes[-1].set_xlabel("time (ns)", fontsize=9)
fig.tight_layout()

# ---------------- overlay helpers ----------------
FONT_PATH = fm.findfont(fm.FontProperties(family="DejaVu Sans", weight="bold"))


def F(size):
    return ImageFont.truetype(FONT_PATH, size)


def fit_font(draw, text, max_w, size=26, min_size=12):
    while size > min_size:
        if draw.textlength(text, font=F(size)) <= max_w:
            return F(size)
        size -= 2
    return F(min_size)


def frosted(img, box, radius=8):
    x0, y0, x1, y1 = [int(round(v)) for v in box]
    x0, y0, x1, y1 = max(x0, 0), max(y0, 0), min(x1, img.width), min(y1, img.height)
    if x1 <= x0 or y1 <= y0:
        return
    region = img.crop((x0, y0, x1, y1)).filter(ImageFilter.GaussianBlur(BLUR))
    region = Image.blend(region, Image.new("RGB", region.size, (255, 255, 255)), WHITE_ALPHA)
    mask = Image.new("L", region.size, 0)
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, region.width - 1, region.height - 1), radius=radius, fill=255)
    img.paste(region, (x0, y0), mask)


def overlap(a, b):
    return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])


def make_projector(view, fov_deg, W, H):
    R = np.array(view[0:9]).reshape(3, 3)
    pos = np.array(view[9:12])
    origin = np.array(view[12:15])
    f = 1.0 / math.tan(math.radians(fov_deg) / 2.0)

    def project(xyz):
        c = R.T @ (np.asarray(xyz) - origin) + pos
        z = -c[2]
        if z <= 1e-3:
            return None
        return ((f * c[0] / z / (W / H) * 0.5 + 0.5) * W, (0.5 - f * c[1] / z * 0.5) * H)
    return project


project = make_projector(lab["view"], lab["fov"], lab["W"], lab["H"])

# contacts table: fixed size for the whole movie, so it does not jump around
probe = ImageDraw.Draw(Image.new("RGB", (10, 10)))
f_row, f_title, f_lab, f_leg = F(TABLE_FONT), F(TABLE_FONT + 3), F(LABEL_FONT), F(TABLE_FONT - 1)
all_entries = sorted({entry(p) for fr in contacts for p in fr}) or ["Pept-Xxx0:Prot-Xxx000"]
colw = int(max(probe.textlength(e, font=f_row) for e in all_entries)) + 20
max_n = max([len(contacts[s - 1]) for s in states if s - 1 < len(contacts)] + [1])
ncols = 1 if max_n <= 8 else 2
rows = math.ceil(max_n / ncols)
PAD, TITLE_H, ROW_H, LEGEND_H = 10, 26, TABLE_FONT + 7, 22
panel_w = max(2 * PAD + ncols * colw, 2 * PAD + int(probe.textlength("Residue contacts (00)", font=f_title)))
panel_h = PAD + TITLE_H + rows * ROW_H + LEGEND_H + PAD // 2

LH, LP = LABEL_FONT + 8, 6                           # label box height and text padding
offsets = lambda w: [(10, -28), (10, 8), (-10 - w, -28), (-10 - w, 8), (10, -52), (10, 34),
                     (-10 - w, -52), (-10 - w, 34), (28, -11), (-28 - w, -11)]

# ---------------- per-frame loop ----------------
for i, k in enumerate(idx):
    for line, dot, y in artists:
        line.set_data(tn[:k + 1], y[:k + 1])
        dot.set_data([tn[k]], [y[k]])
    fig.savefig(os.path.join(mdir, "p_%04d.png" % i), dpi=100)

    st = states[i]
    img = Image.open(os.path.join(mdir, "f_%04d.png" % i)).convert("RGB")
    d = ImageDraw.Draw(img)
    W_, H_ = img.size

    # title (one-letter sequence, shrinks / abbreviates to fit) + time
    seq_txt = pep if len(pep) <= 40 else pep[:18] + "..." + pep[-18:]
    f1 = fit_font(d, seq_txt, W_ - 32)
    d.text((16, 12), seq_txt, fill=(255, 255, 255), font=f1)
    d.text((16, 12 + int(f1.size) + 8), f"rep{rep} | t = {tsel[i] / 1000:.2f} ns", fill=(255, 255, 255), font=F(22))

    # contacts table, bottom-left: one peptide-residue : protein-residue pair per line
    cur = contacts[st - 1] if st - 1 < len(contacts) else []
    hb_cur = hb_res[st - 1] if st - 1 < len(hb_res) else set()
    px0, py0 = 16, H_ - 16 - panel_h
    panel = (px0, py0, px0 + panel_w, py0 + panel_h)
    frosted(img, panel, radius=10)
    d = ImageDraw.Draw(img)
    d.text((px0 + PAD, py0 + 6), f"Residue contacts ({len(cur)})", fill=(15, 15, 15), font=f_title)
    for j, pair in enumerate(cur[:rows * ncols]):
        col, row = divmod(j, rows)
        x, y = px0 + PAD + col * colw, py0 + PAD + TITLE_H + row * ROW_H
        text = entry(pair)
        if tuple(pair) in hb_cur:
            d.rounded_rectangle((x - 3, y - 1, x + probe.textlength(text, font=f_row) + 3, y + ROW_H - 4),
                                radius=4, fill=(255, 225, 0))
        d.text((x, y), text, fill=(15, 15, 15), font=f_row)
    ly = py0 + panel_h - LEGEND_H + 2
    d.rounded_rectangle((px0 + PAD, ly + 2, px0 + PAD + 18, ly + 14), radius=3, fill=(255, 225, 0))
    d.text((px0 + PAD + 25, ly), "hydrogen-bonded now", fill=(15, 15, 15), font=f_leg)

    # frosted residue labels for the current hydrogen bonds
    placed = [panel, (0, 0, W_, 80)]
    items = []
    for text, chain, x, y, z in sorted(lab["frames"][i], key=lambda t: t[3]):
        p = project((x, y, z))
        if p is None or not (0 <= p[0] < W_ and 0 <= p[1] < H_):
            continue
        w, h = int(d.textlength(text, font=f_lab)) + 2 * LP, LH
        box = None
        for dx, dy in offsets(w):
            b = (int(p[0] + dx), int(p[1] + dy), int(p[0] + dx) + w, int(p[1] + dy) + h)
            if b[0] >= 0 and b[1] >= 0 and b[2] <= W_ and b[3] <= H_ and not any(overlap(b, o) for o in placed):
                box = b
                break
        if box is None:
            bx = min(max(int(p[0]) + 10, 0), W_ - w)
            by = min(max(int(p[1]) - 28, 0), H_ - h)
            box = (bx, by, bx + w, by + h)
        placed.append(box)
        items.append((text, chain, p, box))
    for text, chain, p, box in items:
        frosted(img, box)
    d = ImageDraw.Draw(img)
    for text, chain, p, box in items:
        cx, cy = min(max(p[0], box[0]), box[2]), min(max(p[1], box[1]), box[3])
        d.line((p[0], p[1], cx, cy), fill=(255, 255, 255), width=1)
        d.ellipse((p[0] - 3, p[1] - 3, p[0] + 3, p[1] + 3), fill=(255, 225, 0), outline=(0, 0, 0))
        d.text((box[0] + LP, box[1] + 3), text, font=f_lab, fill=(0, 95, 115) if chain == "A" else (140, 0, 140))
    img.save(os.path.join(mdir, "lbl_%04d.png" % i))
print(f"{len(idx)} plot panels and labelled frames written")
