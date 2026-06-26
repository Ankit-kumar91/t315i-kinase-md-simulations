# MD Analysis Pipeline — User Manual

Automated GROMACS molecular dynamics analysis pipeline.  
Produces plots, tables, and a summary report from a completed MD simulation.

---

## Files

| File | Purpose |
|---|---|
| `md_analysis.py` | Main analysis script — never needs to be edited |
| `md_config.ini` | Configuration file — edit this for each simulation |
| `install.py` | One-time dependency installer |

---

## Step 1 — Install Dependencies

Run this **once** on any new computer before your first analysis.

```
python install.py
```

This installs:

| Package | Used for |
|---|---|
| NumPy | Numerical calculations |
| Pandas | Data tables and CSV export |
| Matplotlib | Plotting |
| Seaborn | Plot styling |
| MDAnalysis | Reading GROMACS trajectories (RMSD, RMSF, H-bonds, contacts) |

GROMACS must also be available on your system (for PBC correction, energy extraction, and radius of gyration). On HPC clusters, load it with:

```
module load gromacs
```

---

## Step 2 — Configure md_config.ini

Open `md_config.ini` in any text editor.

**The only required change is `sim_dir`** — set it to the folder containing your GROMACS output files (`.tpr`, `.xtc`, `.edr`, `.gro`, `.ndx`).

```ini
[paths]
sim_dir = /home/user/projects/1IEP_MD
out_dir = /home/user/projects/1IEP_MD/results
```

Use the **full absolute path**.  
`out_dir` is created automatically — you do not need to make it yourself.

### File names are auto-detected

The script scans `sim_dir` and picks the right files automatically.  
You only need to set file names manually if detection fails or you have multiple files of the same type:

```ini
[paths]
tpr = md_production.tpr
xtc = md_production.xtc
```

### Common optional settings

**Change the time axis to nanoseconds:**
```ini
[analysis]
time_unit = ns
```

**Change the gatekeeper residue** (set to 0 to skip):
```ini
[system]
gatekeeper_res = 315
```

**Add a label to your plot titles:**
```ini
[system]
sim_label = ABL1 WT + Imatinib
```

**Set the ligand residue name** (if auto-detection picks the wrong molecule):
```ini
[system]
lig_resname = LIG
```

**Adjust plot appearance:**
```ini
[figure]
dpi            = 300        # 72=screen  150=standard  300=publication
width          = 12
height         = 5
title_fontsize = 14
label_fontsize = 12
ligand_color   = #ff6600
```

---

## Step 3 — Run the Analysis

```
python md_analysis.py --config /path/to/md_config.ini
```

Or pass the simulation folder directly (the script looks for `md_config.ini` inside it):

```
python md_analysis.py /path/to/1IEP_MD
```

### Running multiple simulations

Because `md_analysis.py` is never modified, you can keep one copy of the script and run it against any simulation by pointing to that simulation's config file:

```
python md_analysis.py --config /data/1IEP_MD/md_config.ini
python md_analysis.py --config /data/3IK3_MD/md_config.ini
python md_analysis.py --config /data/T315I_MD/md_config.ini
```

Each `md_config.ini` has its own `sim_dir` and `out_dir`.

---

## What the Script Does

The pipeline runs ten steps in order. Each step is independent — if one fails, the rest continue.

### Step 1 — PBC Correction
Fixes periodic-boundary artefacts using GROMACS `trjconv`.  
Produces a corrected trajectory (`md_center.xtc`) used by all subsequent steps.  
Skipped automatically if the corrected trajectory already exists.

### Step 2 — Thermodynamic Properties
Extracts seven properties from the `.edr` energy file using `gmx energy`:

- Potential energy (kJ/mol)
- Kinetic energy (kJ/mol)
- Total energy (kJ/mol)
- Temperature (K)
- Pressure (bar)
- Volume (nm³)
- Density (kg/m³)

Each property is saved as an individual plot with mean ± 1σ shading.

### Step 3 — RMSD
Root-mean-square deviation of:
- **Protein backbone Cα** — measures overall structural stability
- **Ligand heavy atoms** (fit on backbone) — measures binding pose stability

A value consistently below 2–3 Å indicates a stable simulation.

### Step 4 — RMSF
Per-residue Cα root-mean-square fluctuation — shows which parts of the protein are flexible.  
The gatekeeper residue is highlighted with a vertical marker.

### Step 5 — Radius of Gyration
Protein compactness over time using `gmx gyrate`.  
A stable Rg indicates the protein maintains its folded conformation.

### Step 6 — Protein–Ligand H-bonds
Counts hydrogen bonds between protein and ligand at every frame.  
Reports per-pair occupancy (% of frames each H-bond is present).

### Step 7 — Gatekeeper Distance
Minimum distance between the ligand and the gatekeeper residue (e.g. THR315 in ABL1) over the trajectory.  
Shows whether the ligand maintains contact with this key residue.

### Step 8 — Binding-Pocket Contact Map
Identifies all protein residues within the contact cutoff (default 4.5 Å) of the ligand and reports their occupancy.  
Residues below `min_plot_contact_pct` are hidden from the bar chart.

### Step 9 — Summary Figure
A four-panel figure combining backbone RMSD, ligand RMSD, radius of gyration, and per-residue RMSF.

### Step 10 — Text Report
A plain-text summary of all key numerical results saved to `tables/analysis_report.txt`.

---

## Output Files

All output is written to `out_dir` (default: `results/` inside the simulation folder).

```
results/
├── energy_potential.png
├── energy_kinetic.png
├── energy_total_energy.png
├── energy_temperature.png
├── energy_pressure.png
├── energy_volume.png
├── energy_density.png
├── rmsd.png
├── rmsf.png
├── gyration.png
├── hbonds.png
├── gatekeeper_distance.png
├── pocket_contacts.png
├── summary_panel.png
└── tables/
    ├── rmsd.csv
    ├── rmsf.csv
    ├── hbonds_count.csv
    ├── hbonds_pairs.csv
    ├── gatekeeper_distance.csv
    ├── pocket_contacts.csv
    └── analysis_report.txt
```

---

## Troubleshooting

**"GROMACS not found"**  
Either GROMACS is not installed, or it is not on your PATH.  
On HPC: `module load gromacs`  
If the command is `gmx_mpi`: set `gromacs_cmd = gmx_mpi` in `[system]`

**"Missing files in sim_dir"**  
The script lists what it found. Either `sim_dir` is wrong, or your files have non-standard names.  
Fix by setting explicit file names in `[paths]` of `md_config.ini`.

**"Multiple non-protein residues found"**  
The script picked the first one. If it is wrong, set the correct name:
```ini
[system]
lig_resname = LIG
```

**"Protein_LIG group NOT FOUND"**  
The script will try to create this group automatically with `gmx make_ndx`.  
If that also fails, add the group number manually:
```ini
[groups]
Protein_LIG = 22
```
(The group numbers are printed in the startup table — use the number from the `Protein_LIG` row.)

**MDAnalysis-dependent steps are skipped**  
Run `python install.py` to install missing packages.

---

## Quick-Reference — md_config.ini Settings

| Section | Key | Default | Description |
|---|---|---|---|
| `[paths]` | `sim_dir` | *(required)* | Folder with GROMACS files |
| `[paths]` | `out_dir` | `results` | Output folder |
| `[system]` | `sim_label` | `simulation` | Label in plot titles |
| `[system]` | `gromacs_cmd` | `gmx` | GROMACS executable name |
| `[system]` | `lig_resname` | auto | Ligand residue name |
| `[system]` | `gatekeeper_res` | `315` | Gatekeeper residue number (0 = skip) |
| `[analysis]` | `time_unit` | `ps` | Time axis unit: `ps` or `ns` |
| `[analysis]` | `contact_cutoff_A` | `4.5` | Contact distance cutoff (Å) |
| `[analysis]` | `hbond_distance_A` | `3.5` | H-bond donor–acceptor cutoff (Å) |
| `[analysis]` | `hbond_angle_deg` | `120` | H-bond angle cutoff (°) |
| `[analysis]` | `min_plot_contact_pct` | `10` | Min occupancy % shown in contact map |
| `[figure]` | `dpi` | `150` | Image resolution |
| `[figure]` | `width` / `height` | `10` / `4` | Figure size (inches) |
| `[figure]` | `title_fontsize` | `13` | Plot title font size |
| `[figure]` | `label_fontsize` | `11` | Axis label font size |
| `[figure]` | `font_scale` | `1.1` | Overall text scale multiplier |
| `[figure]` | `*_color` | various | Colours for each trace |
| `[figure]` | `ligand_marker_size` | `60` | Marker size on RMSF plot |
