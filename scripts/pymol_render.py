"""Render frames of the peptide-protein complex. Run with: pymol -cq pymol_render.py
Environment: MV_PDB, MV_XTC, MV_OUT, MV_TARGET, MV_W, MV_H, MV_SMOOTH, MV_THREADS, MV_HB (H-bond pairs json, optional).
Writes f_NNNN.png, states.txt and labels.json (H-bond residue anchors + camera, used by plot_frames.py)."""
import json
import os

import numpy as np
from pymol import cmd
from pymol.cgo import CYLINDER

pdb, xtc, out = os.environ["MV_PDB"], os.environ["MV_XTC"], os.environ["MV_OUT"]
target = int(os.environ.get("MV_TARGET", "300"))
W, H = int(os.environ.get("MV_W", "960")), int(os.environ.get("MV_H", "720"))
smooth = int(os.environ.get("MV_SMOOTH", "5"))
os.makedirs(out, exist_ok=True)
hb_path = os.environ.get("MV_HB", "")
hb_frames = json.load(open(hb_path)) if hb_path and os.path.isfile(hb_path) else None


def dashed(p, q, radius=0.07, dash=0.28, gap=0.2, color=(1.0, 0.9, 0.1)):
    """Yellow dashed line from p to q as a list of CGO cylinders."""
    v = q - p
    length = float(np.linalg.norm(v))
    if length < 1e-6:
        return []
    u = v / length
    obj, pos = [], 0.0
    while pos < length:
        end = min(pos + dash, length)
        a, b = p + u * pos, p + u * end
        obj += [CYLINDER, float(a[0]), float(a[1]), float(a[2]), float(b[0]), float(b[1]), float(b[2]),
                radius, color[0], color[1], color[2], color[0], color[1], color[2]]
        pos += dash + gap
    return obj


cmd.set("max_threads", int(os.environ.get("MV_THREADS", "8")))
cmd.load(pdb, "cx")
cmd.load_traj(xtc, "cx", state=1)
n = cmd.count_states("cx")
if n < 2:
    raise SystemExit("ERROR: trajectory did not load (only %d state)" % n)
if cmd.count_atoms("chain A") == 0 or cmd.count_atoms("chain B") == 0:
    raise SystemExit("ERROR: expected chain A (protein) and chain B (peptide) in the PDB")
if smooth > 1:
    cmd.smooth("cx", 1, smooth)          # light running average, removes thermal jitter

cmd.remove("name H*")
cmd.hide("everything", "cx")
cmd.bg_color("black")
cmd.frame(1)
cmd.select("site", "byres (chain A within 5 of chain B)")
cmd.show("cartoon", "chain A")
cmd.color("gray60", "chain A")
cmd.set("cartoon_transparency", 0.35, "chain A")
cmd.show("sticks", "site or chain B")
cmd.set("stick_radius", 0.15, "site")
cmd.set("stick_radius", 0.28, "chain B")
cmd.color("cyan", "site and name C*")
cmd.color("magenta", "chain B and name C*")
cmd.color("blue", "(site or chain B) and name N*")
cmd.color("red", "(site or chain B) and name O*")
cmd.set("ray_shadow", 0)
cmd.set("antialias", 2)
cmd.set("ray_opaque_background", 1)
cmd.set("depth_cue", 0)
cmd.orient("site or chain B")
cmd.zoom("site or chain B", 6)
cmd.clip("slab", 45)                      # show only the region around the pocket
view0 = cmd.get_view()                    # the camera framing, restored before every frame

stride = max(1, n // target)
states = list(range(1, n + 1, stride))
frames_labels = []
with open(os.path.join(out, "states.txt"), "w") as f:
    for i, s in enumerate(states):
        cmd.frame(s)
        cmd.delete("hb")
        labs = []
        if hb_frames is not None and s - 1 < len(hb_frames):
            obj, seen = [], set()
            for d_id, a_id in hb_frames[s - 1]:
                xyz = {}
                for atom_id in (d_id, a_id):
                    c = cmd.get_coords("id %d" % atom_id, state=s)
                    xyz[atom_id] = None if c is None else c[0]
                if xyz[d_id] is not None and xyz[a_id] is not None:
                    obj += dashed(xyz[d_id], xyz[a_id])
                for atom_id in (d_id, a_id):
                    if xyz[atom_id] is None:
                        continue
                    at = cmd.get_model("id %d" % atom_id).atom[0]
                    key = (at.chain, at.resi)
                    if key in seen:
                        continue
                    seen.add(key)
                    x = xyz[atom_id]
                    labs.append([at.resn + at.resi, at.chain, float(x[0]), float(x[1]), float(x[2])])
            if obj:
                cmd.load_cgo(obj, "hb", zoom=0)
        frames_labels.append(labs)
        cmd.set_view(view0)
        cmd.png(os.path.join(out, "f_%04d.png" % i), W, H, ray=1)
        f.write("%d\n" % s)

json.dump({"view": list(view0), "fov": float(cmd.get("field_of_view")), "W": W, "H": H, "frames": frames_labels},
          open(os.path.join(out, "labels.json"), "w"))
print("rendered %d frames (every %d-th of %d states)" % (len(states), stride, n))
cmd.quit()
