#!/usr/bin/env python3
"""
run_md.py - run the mdpeptprot pipeline from a YAML job file.

    python run_md.py [--dry-run] [--force] job.yaml [job2.yaml ...]

The YAML is the only thing that changes between runs/machines. Values left out
of the YAML keep their defaults from scripts/config.sh.
"""
import argparse
import datetime
import os
import re
import shutil
import subprocess
import sys
import time

try:
    import yaml
except ImportError:
    sys.exit("ERROR: PyYAML is missing (pip install pyyaml).")

ROOT = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(ROOT, "scripts")
ALL_STAGES = ["prepare", "equilibrate", "production", "analyze"]
SCHEMA = {
    "name": None, "input": None, "stages": None,
    "conditions": {"ph", "temperature_k", "salt_molar"},
    "production": {"replicas", "ns"},
    "hardware": {"gpu_id", "cores", "threads_per_replica"},
}
# settings this runner controls: cleared from the caller's environment first,
# so nothing leaks in from the shell.
MANAGED = ["RUN_NAME", "INPUT_PDB", "PH", "TEMPERATURE_K", "ION_CONC", "N_REPLICAS",
           "PRODUCTION_NS", "GPU_ID", "TOTAL_CORES", "MAX_THREADS_PER_REPLICA",
           "WAIT", "FORCE_RELAUNCH"]


class JobError(Exception):
    pass


def _section(job, key):
    sec = job.get(key) or {}
    if not isinstance(sec, dict):
        raise JobError(f"'{key}' must be a mapping")
    extra = set(sec) - SCHEMA[key]
    if extra:
        raise JobError(f"unknown key(s) in '{key}': {sorted(extra)} (allowed: {sorted(SCHEMA[key])})")
    return sec


def _num(sec, section, key, lo, hi, integer=False):
    if sec.get(key) is None:
        return None
    v = sec[key]
    ok = isinstance(v, int) if integer else isinstance(v, (int, float))
    if isinstance(v, bool) or not ok:
        raise JobError(f"{section}.{key} must be {'an integer' if integer else 'a number'}, got {v!r}")
    if not lo <= v <= hi:
        raise JobError(f"{section}.{key} = {v} is outside the supported range {lo}..{hi}")
    return v


def load_job(path):
    with open(path) as f:
        job = yaml.safe_load(f)
    if not isinstance(job, dict):
        raise JobError("the job file must be a YAML mapping")
    extra = set(job) - set(SCHEMA)
    if extra:
        raise JobError(f"unknown top-level key(s): {sorted(extra)} (allowed: {sorted(SCHEMA)})")

    name = job.get("name")
    if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9._-]+", name):
        raise JobError("'name' is required (letters, digits, . _ - only)")
    inp = job.get("input")
    if not isinstance(inp, str):
        raise JobError("'input' is required (path to the protein+peptide PDB)")
    inp_abs = inp if os.path.isabs(inp) else os.path.join(ROOT, inp)
    if not os.path.isfile(inp_abs):
        raise JobError(f"input file not found: {inp}")

    cond = _section(job, "conditions")
    prod = _section(job, "production")
    hw = _section(job, "hardware")

    env = {"RUN_NAME": name, "INPUT_PDB": inp_abs}
    settings = [
        ("PH", _num(cond, "conditions", "ph", 5, 9)),
        ("TEMPERATURE_K", _num(cond, "conditions", "temperature_k", 273, 400)),
        ("ION_CONC", _num(cond, "conditions", "salt_molar", 0, 1)),
        ("N_REPLICAS", _num(prod, "production", "replicas", 1, 8, True)),
        ("PRODUCTION_NS", _num(prod, "production", "ns", 0.01, 1000)),
        ("GPU_ID", _num(hw, "hardware", "gpu_id", 0, 15, True)),
        ("TOTAL_CORES", _num(hw, "hardware", "cores", 1, 1024, True)),
        ("MAX_THREADS_PER_REPLICA", _num(hw, "hardware", "threads_per_replica", 1, 64, True)),
    ]
    for key, val in settings:
        if val is not None:
            env[key] = str(val)

    stages = job.get("stages", ALL_STAGES)
    if not isinstance(stages, list) or not stages or any(s not in ALL_STAGES for s in stages):
        raise JobError(f"'stages' must be a non-empty list drawn from {ALL_STAGES}")
    return name, env, [s for s in ALL_STAGES if s in stages]


def plan(stages):
    """Script codes to run, in order. 01 always runs (idempotent, refreshes the force-field path)."""
    seq = []
    if "prepare" in stages:
        seq.append("00")
    seq.append("01")
    if "prepare" in stages:
        seq += ["02", "03"]
    if "equilibrate" in stages:
        seq += ["04", "05"]
    if "production" in stages:
        seq += ["06", "08"]
    if "analyze" in stages:
        seq.append("11")
    return seq


def script_for(code):
    hits = sorted(f for f in os.listdir(SCRIPTS) if f.startswith(code + "_") and f.endswith(".sh"))
    if len(hits) != 1:
        raise JobError(f"expected exactly one scripts/{code}_*.sh, found {hits}")
    return os.path.join(SCRIPTS, hits[0])


def build_env(extra):
    env = {k: v for k, v in os.environ.items() if k not in MANAGED}
    env.update(extra)
    return env


def is_done(script, env):
    stage = os.path.basename(script)[:-3]
    cmd = 'source scripts/config.sh >/dev/null 2>&1; source scripts/common.sh; stage_is_done "$1"'
    r = subprocess.run(["bash", "-c", cmd, "_", stage], cwd=ROOT, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return r.returncode == 0


def run_stage(script, env, logf):
    p = subprocess.Popen(["bash", script], cwd=ROOT, env=env, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in p.stdout:
        sys.stdout.write(line)
        logf.write(line)
        logf.flush()
    return p.wait()


def run_job(path, dry, force):
    name, extra, stages = load_job(path)
    env = build_env(extra)
    codes = plan(stages)
    scripts = [(c, script_for(c)) for c in codes]

    print(f"\n=== Job: {name}   ({path}) ===")
    for k in MANAGED:
        if k in extra:
            print(f"  {k:<24}= {extra[k]}")
    print(f"  stages: {', '.join(stages)}")

    if dry:
        for code, script in scripts:
            base = os.path.basename(script)
            if code == "01":
                state = "run (always)"
            elif force:
                state = "run (--force)"
            else:
                state = "skip (already done)" if is_done(script, env) else "run"
            print(f"  {base:<32} {state}")
        return True

    os.makedirs(os.path.join(ROOT, "runs"), exist_ok=True)
    shutil.copy(path, os.path.join(ROOT, "runs", f"{name}.job.yaml"))
    t_job = time.time()
    with open(os.path.join(ROOT, "runs", f"{name}.run_md.log"), "a") as logf:
        logf.write(f"\n##### {datetime.datetime.now().isoformat(timespec='seconds')}  {path}\n")
        for code, script in scripts:
            base = os.path.basename(script)
            if code != "01" and not force and is_done(script, env):
                print(f"[run_md] {base}: already done, skipping")
                continue
            stage_env = dict(env)
            if code == "08":
                stage_env["WAIT"] = "1"      # block until every replica really finished
            print(f"[run_md] >>> {base}")
            t0 = time.time()
            rc = run_stage(script, stage_env, logf)
            if rc != 0:
                print(f"[run_md] FAILED at {base} (exit {rc}). Fix it and run again: finished stages are skipped.")
                return False
            print(f"[run_md] <<< {base} ok ({time.time() - t0:.0f} s)")
    print(f"[run_md] Job {name} finished in {(time.time() - t_job) / 60:.1f} min")
    return True


def main():
    sys.stdout.reconfigure(line_buffering=True)
    ap = argparse.ArgumentParser(description="Run the mdpeptprot pipeline from YAML job files.")
    ap.add_argument("jobs", nargs="+", help="one or more job YAML files, run one after another")
    ap.add_argument("--dry-run", action="store_true", help="show the plan, run nothing")
    ap.add_argument("--force", action="store_true", help="re-run stages even if already done")
    args = ap.parse_args()

    failed = []
    for path in args.jobs:
        try:
            ok = run_job(path, args.dry_run, args.force)
        except JobError as e:
            print(f"ERROR in {path}: {e}", file=sys.stderr)
            ok = False
        if not ok:
            failed.append(path)
    if failed:
        print(f"\n{len(failed)} job(s) failed: {', '.join(failed)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
