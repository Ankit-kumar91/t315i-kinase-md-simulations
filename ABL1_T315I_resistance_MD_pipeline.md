# Probing Imatinib Resistance in BCR-ABL1: An MD Study of the T315I Gatekeeper Mutation

**One-line premise:** Use all-atom molecular dynamics + MM-PBSA to reproduce, from first principles, why imatinib fails against the T315I mutant of ABL1 while ponatinib survives it — a textbook drug-resistance mechanism turned into a reproducible computational result.

---

## 1. Scientific question & hypothesis

**Question:** Can physics-based MD distinguish a resistance-conferring mutation from a benign one, and correctly rank a first-generation drug (imatinib) against a later-generation drug (ponatinib) designed to overcome that mutation?

**Hypothesis:**
- The T315I mutation destabilizes the imatinib complex (higher ligand RMSD, loss of the Thr315 H-bond, unfavorable ΔΔG_bind).
- The same mutation has little effect on ponatinib (low ligand RMSD, near-zero ΔΔG_bind).
- Therefore `ΔΔG_bind(imatinib, WT→T315I)  >>  ΔΔG_bind(ponatinib, WT→T315I)`.

**Why this system is strong for a portfolio:** T315I is *the* clinically famous "gatekeeper" mutation in CML; ponatinib (AP24534) was rationally designed to beat it. Recruiters in oncology/comp-chem recognize it instantly, and there is abundant structural + clinical literature to validate against.

---

## 2. Experimental design (the 2×2)

| System | Protein | Ligand | Expected outcome | Role |
|--------|---------|--------|------------------|------|
| S1 | ABL1 WT | imatinib | stable, strong binding | baseline |
| S2 | ABL1 T315I | imatinib | destabilized, weak binding | **resistance signal** |
| S3 | ABL1 WT | ponatinib | stable, strong binding | control |
| S4 | ABL1 T315I | ponatinib | stable, retained binding | **rescue signal** |

The science lives in two differences: **S2 − S1** (imatinib loses binding) vs **S4 − S3** (ponatinib does not).

---

## 3. Required starting structures

Pull these from the RCSB PDB. **Verify the exact codes on rcsb.org before downloading** — codes below are the canonical ones but confirm the contents (correct protein, correct ligand, DFG-out inactive conformation):

| Structure | PDB | Notes |
|-----------|-----|-------|
| ABL1 WT + imatinib | **1IEP** | classic Abl–imatinib (STI-571) complex; imatinib resname `STI`; gatekeeper = THR A 315 (1b numbering) |
| ABL1 WT + ponatinib | **3OXZ** | WT Abl + ponatinib; confirmed THR A 315 in downloaded file |
| ABL1 T315I + ponatinib | **3IK3** | T315I mutant + AP24534 (ponatinib, resname `0LI`); PDB title confirms T315I; 1.90 Å resolution — **note:** residue 315 in 3IK3 is labelled ILE in the file (different chain numbering than 1b), gatekeeper is not at position 315 in this file's numbering |
| ABL1 T315I + imatinib | *build in silico* | no co-crystal exists (that's the point — it doesn't bind well); mutate THR315→ILE from prepared 1IEP REC.pdb using ChimeraX `swapaa` |

So three of the four complexes come straight from crystal structures, and only S2 (T315I + imatinib) is built by mutation. That keeps the binding poses experimentally grounded.

> **Critical numbering caveat:** ABL1 has two isoform numbering schemes (1a vs 1b) that differ by ~19 residues. The clinical "T315I" uses 1b numbering. **1IEP uses 1b numbering** — gatekeeper = THR A 315. **3IK3 uses a different numbering** — the residue at position 315 in the file is a naturally-occurring ILE (not the gatekeeper); the actual gatekeeper in 3IK3 is at a different label. For S2 mutation: always mutate from the **1IEP-derived REC.pdb** where THR 315 is the confirmed gatekeeper.

---

## 4. Software stack

| Purpose | Tool |
|---------|------|
| Structure prep / mutation | **ChimeraX** (`swapaa` command); PyMOL Mutagenesis Wizard also works |
| Protonation states (pH 7.4) | PROPKA / H++ |
| Ligand parameters (CHARMM) | **CHARMm 27** (recommended) or CGenFF server + `cgenff_charmm2gmx.py` |
| System building | **CHARMM-GUI Solution Builder** (recommended) or `pdb2gmx`/`editconf`/`solvate`/`genion` |
| MD engine | GROMACS (CUDA build) |
| Free-energy | gmx_MMPBSA (AmberTools MMPBSA.py port for GROMACS) |
| Analysis | MDAnalysis (Python), gmx tools, matplotlib/seaborn, pandas |

**Strong recommendation:** since you're CHARMM-trained, use **CHARMM-GUI** for steps 5–6. Its *Solution Builder* with a bound ligand handles CGenFF parametrization, protonation, solvation, ionization, and emits GROMACS-ready topology + `.mdp` files in one workflow. It removes ~80% of the setup pain and the most common path errors. Below I give the manual command path too so you understand what it's doing.

---

## 5. Stage 0 — Protein prep & the mutation

```bash
# 1. Clean the PDB: keep one chain, remove crystallographic waters/buffers,
#    keep the ligand. Inspect in ChimeraX first.
# 2. Assign protonation states at pH 7.4 (PROPKA or H++). Pay attention to
#    histidines and the catalytic residues.

# Build the T315I mutant from the prepared 1IEP REC.pdb (ChimeraX):
open /path/to/1IEP_MD/Replica1/REC.pdb
swapaa ILE :315          # picks best rotamer by Dunbrack score + least clash
save /path/to/structures/1IEP_T315I.pdb format pdb models #1

# Verify:
# grep " 315 " 1IEP_T315I.pdb | grep "^ATOM"   → should show ILE A 315

# GUI alternative: Tools > Structure Editing > Rotamers
#   → select THR A 315 → change to ILE → pick lowest-clash rotamer → save

# NOTE: GROMACS EM (steepest descent) will relax the introduced Ile315
# sidechain adequately before production; a separate PyMOL minimization
# is not required when using the standard NVT/NPT equilibration with
# position restraints.
```

**Ligand protonation (do not skip):** both imatinib and ponatinib carry an *N*-methylpiperazine that is **protonated (cationic, +1) at physiological pH** (piperazine pKa ≈ 8). Model them in the protonated state — it materially changes the electrostatics with the pocket. Set this explicitly when you generate ligand parameters.

---

## 6. Stage 1 — Ligand parametrization (swissparam)

Both drugs are drug-like small molecules; swissparam will parametrize them but **check the penalty scores** it reports — high penalties (>~50) mean poorly-covered chemistry that may need validation or QM refinement.

- **Imatinib:** generally well-behaved in swissparam.
- **Ponatinib:** watch the **alkyne (C≡C) linker** — swissparam has triple-bond atom types (`CG1T1`/`CG1T2`) but penalties around the alkyne and the imidazo-pyridazine can be elevated. Note any high-penalty atoms in your README as a documented limitation.

Manual path:
```bash
# Upload ligand mol2/sdf (with correct protonation + formal charge) to the
ParamChem server -> download .str
# Convert to GROMACS:
python cgenff_charmm2gmx.py LIG ligand.mol2 ligand.str charmm36.ff
```
CHARMM-GUI's Ligand Reader does the equivalent and merges it into the system automatically.

---

## 7. Stage 2 — System building

Manual GROMACS path (CHARMM36m force field, TIP3P water):

```bash
# Generate protein topology
gmx pdb2gmx -f protein.pdb -o processed.gro -water tip3p -ff charmm36-jul2022

# Merge protein + ligand topology (insert ligand .itp + coordinates)

# Box: rhombic dodecahedron saves ~30% waters vs cubic -> cheaper sims
gmx editconf -f complex.gro -o box.gro -bt dodecahedron -d 1.2 -c

# Solvate
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

# Add ions: neutralize + 0.15 M NaCl
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o ions.gro -p topol.top -neutral -conc 0.15 -pname NA -nname CL
```

Expected system size: ABL kinase domain (~290 residues) in a 1.2 nm dodecahedral box ≈ **40,000–55,000 atoms**. Throughput on a single T4/P100 ≈ **30–55 ns/day**.

---

## 8. Stage 3 — MD protocol

Standard CHARMM36 GROMACS protocol. The CHARMM-specific cutoff settings below matter — don't use default Amber-style cutoffs.

**Critical `.mdp` settings for CHARMM36 (apply to NVT/NPT/production):**
```
; CHARMM force-switch van der Waals
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
rlist                   = 1.2
; Electrostatics
coulombtype             = PME
rcoulomb                = 1.2
; Constraints / timestep
constraints             = h-bonds
constraint-algorithm    = lincs
dt                      = 0.002
```

**Pipeline:**
1. **Energy minimization** — steepest descent, until Fmax < 1000 kJ/mol/nm.
2. **NVT equilibration** — 200–500 ps, 300 K (or 310 K), V-rescale thermostat, position restraints (1000 kJ/mol/nm²) on protein heavy atoms **and** ligand.
3. **NPT equilibration** — 1 ns, 1 bar, C-rescale or Parrinello-Rahman barostat, restraints retained.
4. **(Recommended) restraint release** — a second 1 ns NPT with restraints reduced/removed, so the mutated S2 system relaxes gently.
5. **Production** — 50–100 ns, no restraints, save coordinates every 10–20 ps.

**Replicas (essential here):** run each system with **independent initial velocities** (`gen-vel = yes`, different `gen-seed`). The resistance effect can be within single-trajectory noise — replicas are what let you put error bars on ΔΔG and claim significance.

**Restart after Kaggle/spot interruptions:**
```bash
gmx mdrun -deffnm prod -cpi prod.cpt -append   # resumes from checkpoint
```

---

## 9. Compute budget

| Plan | Design | Total sim | ~GPU-days (40 ns/day) | ~Cost (GCP T4 spot) |
|------|--------|-----------|----------------------|---------------------|
| Budget | 4 sys × 2 rep × 50 ns | 400 ns | ~10 | **~$80–100** |
| Recommended | imatinib pair 3 rep, ponatinib pair 2 rep, 100 ns | 1000 ns | ~25 | **~$200–250** |
| Gold | 4 sys × 3 rep × 100 ns | 1200 ns | ~30 | **~$290** (risky vs $300) |

Use **spot/preemptible** instances, checkpoint, and **stop instances when idle**. Prioritize replicas on the imatinib pair (S1/S2) — that's the must-have comparison; ponatinib (S3/S4) is the bonus that completes the narrative.

---

## 10. Stage 4 — Trajectory analysis

Run the same battery on every replica of every system; report mean ± SD across replicas.

| Metric | What it tells you | Expected result |
|--------|-------------------|-----------------|
| Protein backbone RMSD | convergence/stability sanity check | all plateau |
| **Ligand RMSD** (fit on protein) | **headline stability metric** | imatinib drifts in T315I; stable elsewhere |
| RMSF per residue | flexibility of P-loop / activation loop | local changes near pocket |
| **H-bond occupancy to gatekeeper (315)** | the mechanistic core for imatinib | present in S1/S3 WT, lost in S2 T315I+imatinib |
| **Gatekeeper hydrophobic contact (315)** | the mechanistic core for ponatinib | S4 (T315I+ponatinib): ILE315 maintains hydrophobic contact; do NOT use H-bond threshold (<3.5 Å) for ponatinib — it bypasses the gatekeeper H-bond by design; use contact fraction with cutoff ~4.5 Å instead |
| Ligand–residue-315 min distance / contacts | shows the steric/polar switch | imatinib: H-bond lost in T315I; ponatinib: hydrophobic contact retained |
| Pocket contact map | overall interaction pattern | imatinib loses contacts in mutant; ponatinib contact map largely preserved |
| (optional) ligand pose clustering | visualize the pose shift | imatinib re-poses in T315I |

**Example: ligand RMSD with MDAnalysis**
```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms

u = mda.Universe("prod.tpr", "prod_nojump.xtc")
R = rms.RMSD(u, u,
             select="protein and name CA",      # superpose on protein
             groupselections=["resname STI"])    # measure ligand drift
R.run()
# R.results.rmsd[:, 3] = ligand RMSD over time (after protein fit)
```

**Example: gatekeeper contact / H-bond tracking**
```python
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA

hbonds = HBA(universe=u,
             between=["resname STI", "resid <GATEKEEPER_RESID>"])
hbonds.run()
# occupancy = fraction of frames the Thr315–ligand H-bond is present
```

---

## 11. Stage 5 — MM-PBSA & the headline ΔΔG

Use **gmx_MMPBSA** on the equilibrated portion of each trajectory (single-trajectory protocol: extract complex/receptor/ligand from the *same* trajectory — this cancels intramolecular errors and is the standard approach). You don't need the full 100 ns for the bound-state estimate; sample evenly across the converged window (e.g. every 10 ps over the last 30–50 ns).

```bash
gmx_MMPBSA -O -i mmpbsa.in -cs prod.tpr -ci index.ndx \
           -cg <receptor_group> <ligand_group> \
           -ct prod_nojump.xtc -cp topol.top
```

Report, per system, mean ± SD of ΔG_bind across replicas, then compute the two key deltas:

- `ΔΔG_imatinib = ΔG(S2) − ΔG(S1)`  → expect **large, unfavorable** (positive, several kcal/mol)
- `ΔΔG_ponatinib = ΔG(S4) − ΔG(S3)`  → expect **near zero**

Add **per-residue decomposition** focused on the gatekeeper — showing residue 315's energetic contribution flip between imatinib and ponatinib is a beautiful, depth-signaling figure.

**Optional, cheap, impressive:** add an **Interaction Entropy** estimate (gmx_MMPBSA supports it) so you're not silently dropping the entropy term — and say so.

### Headline result table (fill in)

| System | ΔG_bind (kcal/mol) | | ΔΔG vs WT |
|--------|--------------------|--|-----------|
| S1 WT + imatinib | −__ ± __ | | — |
| S2 T315I + imatinib | −__ ± __ | | **+__** (resistance) |
| S3 WT + ponatinib | −__ ± __ | | — |
| S4 T315I + ponatinib | −__ ± __ | | **~0** (retained) |

If that table shows a big positive ΔΔG for imatinib and ~0 for ponatinib, your single figure tells the entire clinical story.

---

## 12. Simulation progress (updated 2026-06-24)

| System | PDB | Replicas | Status |
|--------|-----|----------|--------|
| S1 — WT + imatinib | 1IEP | 3 × 100 ns | ✅ Complete — all analysis done |
| S2 — T315I + imatinib | *in silico* (from 1IEP) | 0 | ❌ Structure not yet built — next priority |
| S3 — WT + ponatinib | 3OXZ | in progress | 🔄 Running |
| S4 — T315I + ponatinib | 3IK3 | Rep2 done (100 ns); Rep1 + Rep3 pending | 🔄 In progress |

**S1 key findings (cross-replica, 3 × 100 ns):**
- Backbone RMSD: 2.26 ± 0.43 Å (stable; transient >3 Å excursions in R1/R2)
- Ligand RMSD: 1.12 ± 0.11 Å (tightly bound throughout)
- THR315 gatekeeper contact: **>95% occupancy in all replicas** (mean 3.10 ± 0.06 Å)
- Dominant H-bonds: MET318(N)→imatinib 98.1%, ASP381(N)→imatinib 94.2%
- Core binding shell (≥80% in all 3 replicas): ALA269, ALA380, ASP381, GLU286, GLU316, ILE313, ILE360, LEU248, LEU370, LYS271, MET290, MET318, PHE317, PHE382, THR315, VAL289, VAL299

**S4 preliminary findings (3IK3 Rep2, 1 × 100 ns):**
- ILE315 contact: 36.4% at <3.5 Å — **expected; do not interpret as failure** (ponatinib uses hydrophobic not H-bond contact; metric threshold is wrong for this system)
- Ligand RMSD: 1.32 Å (stable binding despite gatekeeper mutation)
- Larger binding pocket (25 residues ≥80% vs 17 for imatinib) — consistent with ponatinib's broader contacts

---

## 13. System-specific pitfalls (read before you start)

1. **Residue numbering** — the #1 mistake. Confirm the gatekeeper threonine's identity structurally before mutating (see §3).
2. **DFG-out conformation** — imatinib only binds the *inactive* DFG-out state. Start from and preserve it; do not use an active-state structure.
3. **Ligand protonation** — model the *N*-methylpiperazine as cationic (+1) at pH 7.4 for both drugs.
4. **CGenFF penalties** — flag high-penalty atoms (esp. ponatinib's alkyne) as a documented limitation.
5. **Relax the mutant gently** — minimize/equilibrate S2 with restraints before unrestrained production, or the introduced Ile315 clash can crash the run.
6. **MM-PBSA caveats** — no configurational entropy by default, sensitive to internal dielectric and salt; single-trajectory protocol; treat as a *ranking* tool, not absolute ΔG. State all of this.
7. **Convergence & replicas** — the resistance effect can be within noise; replicas + error bars are non-negotiable for the claim.

---

## 14. Suggested repository structure

```
abl1-t315i-resistance-md/
├── README.md                 # mini-paper: motivation, methods, results, limitations
├── structures/               # raw + cleaned PDBs, mutant build notes
├── ligands/                  # CGenFF .str/.itp, protonation notes, penalty report
├── systems/                  # built systems S1–S4 (gro/top, mdp files)
├── mdp/                       # em / nvt / npt / prod .mdp
├── runs/                      # per-system, per-replica run dirs (large files gitignored)
├── analysis/
│   ├── rmsd.py               # MDAnalysis ligand/protein RMSD
│   ├── hbonds.py             # gatekeeper contact/H-bond tracking
│   ├── mmpbsa/               # gmx_MMPBSA inputs + per-residue decomposition
│   └── figures.ipynb        # publication-quality plots
├── results/                  # ΔG table, figures
└── environment.yml           # reproducible env
```

---

## 15. What turns this from "homework" into a portfolio piece

- **Replicas + error bars** on every reported number.
- A **README written like a mini-paper** — explicit methods (force field, water model, box, protocol, .mdp settings), embedded figures, and an honest **limitations** section.
- **Full reproducibility**: mdp files, environment spec, analysis notebook, run instructions.
- The **per-residue decomposition** flipping at residue 315 — your one figure that ties physics to clinical fact.
- **Stretch (optional):** add a benign mutation as a negative control, or feed the per-residue interaction energies / dynamic descriptors into the QSAR/ML layer from your other project to bridge MD → ML in one story.

---

*Targets and structures referenced here are well-established in the literature; verify all PDB codes and residue numbering against RCSB before running. This is a study-design blueprint, not clinical or therapeutic guidance.*
