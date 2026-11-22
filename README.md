# mdpeptprot: molecular dynamics filtering of peptide-protein complexes

Workflow:

```
germinatedPeptides -> probabilisticDocking -> mdpeptprot
```

- `germinatedPeptides`: physicochemical filtering
- `probabilisticDocking`: biophysical filtering (docking)
- `mdpeptprot`: molecular dynamics filtering (this repository)

Input: one docked complex (a protein chain + a 3-6 residue peptide, no PTMs).
Output: whether the peptide **stays bound** in explicit water at physiological conditions, **how** it is held
(which residues, which hydrogen bonds), whether independent replicas **agree**, and two videos of the binding.

---

## Quick start

```bash
conda activate mdpeptprot
python3 run_md.py --dry-run jobs/1WCY.yaml     # show the plan, run nothing
python3 run_md.py jobs/1WCY.yaml               # run every stage listed in the job file
python3 run_md.py jobs/a.yaml jobs/b.yaml      # several jobs, one after another
python3 run_md.py --force jobs/1WCY.yaml       # re-run stages even if already finished
```

`run_md.py` reads the YAML job, checks every value, sets the settings that `scripts/config.sh` understands, and
runs the stages in order. Finished stages are skipped, so a stopped job resumes where it left off; a failed stage
stops the job and names the stage. After changing parameters of a finished run, use `--force` or a new `name`.

### Job file (YAML): the only thing that changes between runs and machines

| key | meaning | default |
|---|---|---|
| `name` | run folder name, `runs/<name>/` | required |
| `input` | PDB with the protein and the peptide chain | required |
| `conditions.ph` | pH; drives all protonation | 7.4 |
| `conditions.temperature_k` | temperature, K | 310 |
| `conditions.salt_molar` | NaCl concentration, M | 0.15 |
| `production.replicas` | independent replicas (different random seeds), run concurrently on one GPU | 3 |
| `production.ns` | length of each replica, ns | `config.sh` |
| `hardware.gpu_id` | GPU to use | 0 |
| `hardware.cores` / `hardware.threads_per_replica` | optional: total CPU cores / threads per replica (default: detected, max 4 per replica) | auto |
| `stages` | any of `prepare`, `equilibrate`, `production`, `analyze` | all four |

Unknown keys are rejected, so a typo cannot be silently ignored.

---

## What each script does

Scripts live in `scripts/` and are run in this order by `run_md.py`. Each can also be run alone with
`RUN_NAME=<name> bash scripts/<script>`. Stage names refer to the `stages:` list of the job file.

| script | stage | what it does | main outputs |
|---|---|---|---|
| `config.sh` | | All settings (paths, conditions, box, equilibration and production lengths, hardware). Every value can be overridden from the environment. | |
| `common.sh` | | Shared helpers: logging, checks, stage markers, replica progress and completion tests. | `stage_timings.tsv` |
| `00_prepare_input.sh` | prepare | Validates the input: a single model, only standard amino acids (no water, sugars or ligands), exactly two chains, peptide of 3-6 residues, protein of at least 50, backbone atoms complete, peptide in contact with the protein (< 4.5 A), no severe clashes (< 1.5 A). The shorter chain is the peptide; the output always has protein = chain A, peptide = chain B. | `input_topology.pdb`, `00_input_report.json` |
| `01_setup_forcefield.sh` | prepare | Extracts the CHARMM36 (jul2022) force field once into `forcefields/_extracted/` and checks that TIP3P water and the chosen ions exist. Always runs first. | `.gmxlib_env` |
| `02_pdb2gmx.sh` + `protonate_from_pqr.py` | prepare | Protonates every residue **at the chosen pH** (PDB2PQR + PROPKA), converts the result to CHARMM names (HSD/HSE/HSP, ASPP, GLUP, LSN), and builds the topology with `pdb2gmx`. Termini are charged (NH3+ / COO-). His, Asp, Glu, Lys and Cys follow PROPKA (a Cys predicted deprotonated becomes CYM). A Tyr predicted deprotonated stops the run with an error, because the force field has no such residue. Nothing is hard-coded per residue. | `topol.top`, `topol_*.itp`, `02_protonation.tsv` |
| `03_solvate_ionize.sh` | prepare | Dodecahedral box (1.0 nm margin), TIP3P water, NaCl at the chosen concentration, neutralized. Fails if the net charge is not zero. | `ionized.gro` |
| `04_minimize.sh` | equilibrate | Energy minimization (steepest descent). Fails if it did not converge. | `minim.gro` |
| `05_equilibrate.sh` | equilibrate | NVT then NPT, with position restraints on the solute, at the target temperature. A **gate** then checks the result: mean temperature within 3 K of the target, density between 950 and 1100 kg/m3, density drift between the last two blocks < 0.3%. | `npt.gro`, `05_equilibration_report.json` |
| `06_setup_production.sh` | production | Builds one production input per replica (new random velocities and seed each), plus `protein.tpr`, the solute-only topology used by the analysis. | `rep*/md.tpr`, `rep*/protein.tpr` |
| `07_benchmark.sh` | (manual) | Optional. Measures speed for different thread counts and prints the longest run that fits a time budget. Run once per machine, not per complex. | `07_benchmark.tsv` |
| `08_launch_production.sh` + `run_replica.sh` | production | Starts all replicas at once on the same GPU (GPU-resident: non-bonded, PME, bonded and update on the GPU, 4 CPU threads pinned per replica). Safe to re-run: finished replicas are skipped, running ones left alone, stopped ones resume from their checkpoint. A replica counts as finished only when its log reaches the target time. `WAIT=1` blocks until all finish (used by `run_md.py`). | `rep*/md.xtc`, `md.edr`, `md.log`, `md.cpt`, `.production.done` |
| `09_monitor.sh` | | Read-only. Progress, speed (ns/day) and ETA of each replica, plus GPU status. Use `watch -n 60 "RUN_NAME=<name> bash scripts/09_monitor.sh"`. | |
| `10_stop_production.sh` | | Stops all replicas gracefully (same as one Ctrl-C: finish the step, write a clean checkpoint). Never use `kill -9` on `mdrun`. Resume with `08`. | |
| `11_analyze.sh` + `analyze_complex.py` | analyze | Analyzes every replica (see "Analysis outputs"). Refuses unfinished replicas unless `ANALYZE_PARTIAL=1`. | `analysis/`, `rep*/md_clean.xtc` |
| `12_movies.sh` + `pymol_render.py` + `plot_frames.py` | (manual) | Two videos per replica (see "Movies"). Needs the separate `mdmovie` conda env (PyMOL + ffmpeg). | `movies/` |
| `progress.sh` | | Older progress helper, superseded by `09_monitor.sh`. | |

---

## Run folder layout

```
runs/<name>/
  00_input_report.json  input_topology.pdb     checked input (protein = chain A, peptide = chain B)
  02_protonation.tsv                           protonation state chosen for every residue at this pH
  topol.top, topol_*.itp, posre_*.itp          topology
  ionized.gro, minim.*, nvt.*, npt.*           system building and equilibration
  05_equilibration_report.json                 equilibration gate result
  rep1 ... repN/
    md.xtc          trajectory: protein + peptide only (no water), one frame every 10 ps
    md.edr, md.log  energies and the GROMACS log (final speed is at the end of md.log)
    md.cpt          checkpoint, used to resume or extend
    protein.tpr     solute-only topology, used by the analysis
    md_clean.xtc    trajectory made whole, centred and fitted on the protein (for viewing and movies)
    complex_start.pdb   start structure for md_clean.xtc (chains A and B, original residue numbers)
  analysis/                                    see below
  movies/               rep*_structure.mp4, rep*_structure_plots.mp4
runs/<name>.job.yaml                           copy of the job file that produced the run
runs/<name>.run_md.log                         everything the stages printed
```

---

## Analysis outputs (`runs/<name>/analysis/`)

`report.md` is the readable summary. The same numbers are in `summary.tsv` (one row per replica plus a mean row,
easy to stack across peptides) and `summary.json`. Per replica: `rep*_timeseries.tsv` (every frame),
`rep*_contacts.tsv`, `rep*_hbonds.tsv`. `rep*_hbond_pairs.json`, `rep*_contact_pairs.json` and
`rep*_contact_residues.json` hold per-frame details used by the movies.

The first 20% of every trajectory is dropped from all averages as settling time (`--skip-frac`).
Residue numbers follow the **original PDB numbering**; peptide residues are numbered as in the input.

| quantity | what it is | how to read it |
|---|---|---|
| **Peptide RMSD (A)** | Movement of the peptide backbone (N, CA, C, O) away from its start pose, after superposing the protein's CA atoms. The reference is the first production frame. | Low and flat: the peptide holds its pose. A jump of ~1.5 A in the first 10 ps is normal (equilibration held it with restraints). Steadily rising or large: the peptide is moving or leaving. |
| **COM shift (A)** | How far the peptide's centre of mass has moved from its start position (protein fitted). | Small: stays put. The report lists the maximum. Above 10 A counts as having left. |
| **Min. distance to protein (A)** | Shortest distance between any peptide heavy atom and any protein heavy atom. | About 2.5-2.8 A while hydrogen-bonded. Above 4.5 A means no contact in that frame. |
| **Contact residues** | Number of protein residues with any heavy atom within 4.5 A of any peptide heavy atom. | A stable number means a stable interface; a falling number means the peptide is letting go. |
| **Hydrogen bonds** | Number of peptide-protein hydrogen bonds per frame: donor (N-H or O-H) to acceptor distance <= 3.5 A and D-H...A angle >= 150 degrees. | Water is not saved in `md.xtc`, so water-mediated bridges are not counted. |
| **Contacts held** | For each protein residue, the fraction of analysed frames in contact with the peptide. | 1.00 = always. Compare across replicas: residues held in all replicas are the reliable interactions. |
| **Hydrogen-bond occupancy** | For each donor-acceptor pair, the fraction of analysed frames it exists. | The strongest pairs are the anchors of the binding mode. A pair present in one replica only is a sign of an alternative pose. |
| **Peptide flexibility (RMSF, A)** | Per peptide residue, how much its heavy atoms fluctuate around their average position. | Higher = more mobile end or side chain. |
| **bound fraction / status** | `bound` if the peptide is in contact in at least 90% of analysed frames and its centre of mass never moved more than 10 A; otherwise `LEFT SITE`. | These thresholds are initial values and are **not calibrated**; run a negative control (a poorly docked peptide) to set them. |
| **mean +- sd** | Mean of the replica means, and the standard deviation **between replicas**. | A large sd means the replicas disagree: extend the run before trusting it. |

`rep*_timeseries.tsv` columns: `time_ps`, `pep_rmsd_A`, `pep_com_shift_A`, `min_dist_A`, `contact_residues`, `hbonds`.

---

## Movies

`RUN_NAME=<name> REP=1 bash scripts/12_movies.sh` writes, for one replica:

- `rep<N>_structure.mp4`: the 3D binding only.
- `rep<N>_structure_plots.mp4`: the same 3D view on the left; on the right, four plots that draw themselves as the movie
  plays (peptide RMSD, minimum distance to protein, contact residues, hydrogen bonds), each with a dot at the current time.

Legend of the 3D view: grey cartoon = protein; cyan sticks = protein residues within 5 A of the peptide at the start;
magenta sticks = the peptide; **yellow dashes = hydrogen bonds present in that frame**, each with a frosted label for the
residue at both ends (teal = protein, purple = peptide). The top-left shows the peptide sequence (one-letter code),
the replica and the time. The table in the bottom-left corner lists the residue pairs in contact in that frame,
as `Pept-Ile1:Prot-Glu205` (any heavy atoms within 4.5 A), updating whenever a contact forms or breaks; pairs that are
also hydrogen-bonded at that moment are highlighted in yellow.

Options (environment variables): `REP` (replica), `MOVIE_FRAMES` (frames to render, default 300), `MOVIE_FPS`
(default 10; lower = slower video), `MOVIE_CPUS` (restrict to cores, e.g. `16-31` while a simulation runs),
`KEEP_FRAMES=1` (keep the PNG frames), `REUSE_FRAMES=1` (redo labels, table, plots and encoding without re-rendering),
`MOVIE_ENV` (conda env with PyMOL and ffmpeg, default `mdmovie`).

---

## Simulation settings

| | |
|---|---|
| Engine | GROMACS 2026.3 (CUDA, thread-MPI) |
| Force field / water | CHARMM36 (jul2022), TIP3P, explicit solvent |
| Conditions | pH 7.4 (PROPKA, drives all protonation including termini), 310 K, 0.15 M NaCl |
| Box | dodecahedron, 1.0 nm solute-box distance; ~87,000 atoms for DPP-IV + tripeptide |
| Electrostatics / coupling | PME; V-rescale thermostat; C-rescale barostat; LINCS on H-bonds; 2 fs time step |
| Equilibration | NVT 100 ps + NPT 1000 ps with position restraints, followed by the gate described above |
| Production | independent replicas with different seeds, 10 ps between saved frames |
| Speed (RTX 5060, 3 replicas at once) | about 45 ns/day per replica, ~135 ns/day combined |

## Requirements

- An NVIDIA GPU and driver.
- **With Docker (recommended):** nothing else; everything except the movies is inside the image.
- **Without Docker (conda):** two environments, created from the files in the repository:
  ```bash
  conda env create -f environment.yml          # pipeline: GROMACS 2026.3 (CUDA), PDB2PQR, PROPKA, MDAnalysis, ...
  conda env create -f environment-movies.yml   # movies only: PyMOL and ffmpeg
  conda activate mdpeptprot
  ```
  The CUDA build of GROMACS is solved only on a machine where the NVIDIA driver is visible to conda.

## Docker

```bash
docker build -t mdpeptprot:dev .
docker run --rm --gpus all --user "$(id -u):$(id -g)" \
  -v "$PWD/inputs:/app/inputs" -v "$PWD/jobs:/app/jobs" -v "$PWD/runs:/app/runs" \
  mdpeptprot:dev jobs/1WCY.yaml
```

The image contains GROMACS, the force field, PDB2PQR/PROPKA, MDAnalysis and the scripts; inputs, job files and results
live in the mounted folders. Movies are not part of the image. Validated: GPU access and the complete four-stage YAML job
(prepare, equilibrate, production, analyze) inside the container, 3 replicas x 0.2 ns in 18 minutes on an RTX 5060.

## Running a new complex on a new machine

### Once per machine

1. A Linux machine with an NVIDIA GPU. Install the NVIDIA driver, Docker Engine and the NVIDIA Container Toolkit.
   Check that Docker sees the GPU (it should print the GPU table):
   ```bash
   docker run --rm --gpus all nvidia/cuda:12.6.0-base-ubuntu22.04 nvidia-smi
   ```
2. Get the repository and build the image (a few minutes, once):
   ```bash
   git clone <repository-url> mdpeptprot && cd mdpeptprot
   docker build -t mdpeptprot:dev .
   ```

Nothing else is installed on the host: GROMACS, PDB2PQR, PROPKA and the force field are inside the image.
The CHARMM36 force field ships in the repository as `forcefields/charmm36-jul2022.ff.tgz`. There is no manual unpacking:
stage 01 extracts it automatically into `forcefields/_extracted/` the first time a job runs (also inside the container).

### For each complex

1. **Check the PDB.** Exactly two chains: the protein and a peptide of 3-6 residues. Only the 20 standard amino acids,
   so remove waters, ligands, ions and sugars first. Stage 00 rejects anything else and says why.
2. **Put it in `inputs/`.**
3. **Write the job file.** Copy `jobs/1WCY.yaml` to `jobs/<name>.yaml` and change `name` and `input`. Change pH,
   temperature, salt or `production.ns` only if needed. This file is the only thing that differs between runs and machines.
4. **Dry run** (shows the plan, runs nothing):
   ```bash
   docker run --rm --user "$(id -u):$(id -g)" \
     -v "$PWD/inputs:/app/inputs" -v "$PWD/jobs:/app/jobs" -v "$PWD/runs:/app/runs" \
     mdpeptprot:dev --dry-run jobs/<name>.yaml
   ```
5. **Run in the background:**
   ```bash
   docker run -d --init --name md_job --rm --gpus all --user "$(id -u):$(id -g)" \
     -v "$PWD/inputs:/app/inputs" -v "$PWD/jobs:/app/jobs" -v "$PWD/runs:/app/runs" \
     mdpeptprot:dev jobs/<name>.yaml
   ```
   The terminal is free immediately. Everything the stages print is also written to `runs/<name>.run_md.log`.
6. **Watch it:**
   ```bash
   tail -f runs/<name>.run_md.log
   watch -n 30 "docker exec -e RUN_NAME=<name> md_job bash scripts/09_monitor.sh"   # once production has started
   ```
   Stop gracefully with `docker exec -e RUN_NAME=<name> md_job bash scripts/10_stop_production.sh`.
   Do not use `kill -9`. Running the same job file again resumes: finished stages are skipped.
7. **Read the results:** `runs/<name>/analysis/report.md` (and `summary.tsv` to compare complexes). Copy `runs/<name>/analysis/`
   off the machine when done.
8. **Movies (optional).** They need the `mdmovie` conda env, so they run outside the container on a machine that has it:
   `RUN_NAME=<name> REP=1 bash scripts/12_movies.sh`.

Not yet implemented: extending the best peptides to longer runs.

### Many complexes (batch)

Give one job file per complex. Put the PDBs in `inputs/batch/` and generate the job files from a template:

```bash
mkdir -p jobs/batch
for f in inputs/batch/*.pdb; do
  n=$(basename "$f" .pdb)
  sed -E "s|^name:.*|name: $n|; s|^input:.*|input: inputs/batch/$n.pdb|" jobs/1WCY.yaml > jobs/batch/$n.yaml
done
```

Then run them all in one go (one after another, because each complex already uses the whole GPU):

```bash
docker run -d --init --name md_job --rm --gpus all --user "$(id -u):$(id -g)" \
  -v "$PWD/inputs:/app/inputs" -v "$PWD/jobs:/app/jobs" -v "$PWD/runs:/app/runs" \
  mdpeptprot:dev jobs/batch/*.yaml
```

If a job fails (a bad YAML or a stage that stops), the failed stage is printed and the next job still runs. At the end the
failed jobs are listed and the exit code is 1. Run the same command again to resume: finished stages are skipped.
Compare the results with `runs/*/analysis/summary.tsv`.

## Limitations

- Explicit solvent only; no post-translational modifications; peptides of 3-6 residues.
- CHARMM36 has no deprotonated Tyr, so a system where PROPKA predicts one (very high pH or an unusual pKa) is rejected
  with an error rather than built with the wrong charge.
- 5 ns per replica is a **screening** length. It shows whether a peptide stays in the site and how consistently across
  replicas, not binding free energies. Extend the promising peptides before drawing conclusions.
- The "bound" thresholds are uncalibrated (see above).
- Water is not saved, so water-mediated hydrogen bonds are not analysed.

---

# END