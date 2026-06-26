#!/usr/bin/env python3
"""
md_analysis.py — Automated GROMACS MD analysis pipeline

Settings are read from md_config.yaml (searched next to this script by default).
Every setting can be overridden on the command line — see --help.

Usage
─────
    python md_analysis.py
    python md_analysis.py --config /path/to/md_config.yaml
    python md_analysis.py --sim-dir /data/my_sim --out-dir /data/results
    python md_analysis.py --dpi 300 --ligand-color "#ff6600" --ligand-size 80
    python md_analysis.py --width 14 --height 5
    python md_analysis.py --help

Analyses performed
──────────────────
1.  PBC correction (gmx trjconv)
2.  Thermodynamic properties (gmx energy) — potential, kinetic, total energy,
    temperature, pressure, volume, density
3.  Protein backbone RMSD
4.  Ligand (LIG / imatinib) RMSD
5.  Per-residue Cα RMSF
6.  Radius of gyration (gmx gyrate)
7.  H-bond occupancy — protein↔ligand
8.  Ligand–gatekeeper (THR315) minimum distance
9.  Binding-pocket contact map (residues within cutoff Å of ligand)
10. Summary report (text + multi-panel figure)

Dependencies
────────────
    pip install MDAnalysis matplotlib seaborn pandas numpy pyyaml
    GROMACS ≥ 2021 must be on your PATH (or set gromacs_cmd in the config)
"""

import sys
import subprocess

# ──────────────────────────────────────────────────────────────────────────────
# AUTO-INSTALL MISSING DEPENDENCIES
# Runs before any third-party import so missing packages are handled gracefully.
# ──────────────────────────────────────────────────────────────────────────────
def _bootstrap() -> None:
    """Check for required packages and install any that are missing via pip."""
    import importlib.util

    # (import_name, pip_package_name)
    required = [
        ("numpy",       "numpy"),
        ("pandas",      "pandas"),
        ("matplotlib",  "matplotlib"),
        ("seaborn",     "seaborn"),
        ("MDAnalysis",  "MDAnalysis"),
    ]

    missing_pkgs = [pkg for mod, pkg in required
                    if importlib.util.find_spec(mod) is None]

    if not missing_pkgs:
        return

    print("=" * 62)
    print("  [bootstrap] Missing Python packages detected:")
    for p in missing_pkgs:
        print(f"              • {p}")
    print("  [bootstrap] Installing now via pip …")
    print("=" * 62)

    try:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "--quiet"] + missing_pkgs
        )
        print("\n  [bootstrap] Installation complete — restarting script …\n")
    except subprocess.CalledProcessError as exc:
        print(f"\n  [bootstrap] pip install failed (exit code {exc.returncode}).")
        print("  Please install manually and re-run:")
        print(f"      pip install {' '.join(missing_pkgs)}")
        sys.exit(1)

    # Restart so the newly installed modules are importable in a fresh process
    result = subprocess.run([sys.executable] + sys.argv)
    sys.exit(result.returncode)


_bootstrap()
# ──────────────────────────────────────────────────────────────────────────────

import argparse
import re
import textwrap
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

warnings.filterwarnings("ignore")

import configparser

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA
    HAS_MDA = True
except ImportError:
    HAS_MDA = False
    print("WARNING: MDAnalysis not available — MDAnalysis-dependent steps skipped.\n")


# ──────────────────────────────────────────────────────────────────────────────
# GLOBAL CONFIG VARIABLES — populated by _configure() before any analysis runs
# ──────────────────────────────────────────────────────────────────────────────
SIM_DIR = OUT_DIR = TABLES_DIR = TPR = XTC = EDR = GRO = NDX = None
XTC_NOJUMP = XTC_CENTER = None
LIG_RESNAME    = "LIG"
GATEKEEPER_RES = 315
#SIM_LABEL      = "simulation"
GROMACS        = "gmx"
GRP            = {}
COLORS         = {}

PLOT_DPI         = 150
PLOT_W           = 10
PLOT_H           = 4
LIG_MARKER_SIZE  = 60
CONTACT_CUTOFF   = 4.5
HBOND_DIST       = 3.5
HBOND_ANGLE      = 120
MIN_CONTACT_PCT  = 10


# ──────────────────────────────────────────────────────────────────────────────
# CONFIG LOADING
# ──────────────────────────────────────────────────────────────────────────────
def _load_config(config_path: Path) -> dict:
    """
    Load an INI config file (md_config.ini) using Python's built-in
    configparser and return it as a nested dict {section: {key: value}}.
    No extra packages needed — configparser is part of standard Python.
    """
    if not config_path.exists():
        print(f"[ERROR] Config file not found: {config_path}")
        print( "        Expected location:  " + str(config_path))
        print( "        Make sure md_config.ini is in the same folder as md_analysis.py")
        sys.exit(1)

    cp = configparser.RawConfigParser(
        comment_prefixes=("#", ";"),   # whole-line comments with # or ;
        inline_comment_prefixes=None,  # disabled so hex colours like #d7191c work fine
        strict=False,                  # allows the same section to appear more than once
    )
    cp.read(config_path, encoding="utf-8")
    # Collapse duplicate sections: later keys override earlier ones within the same section
    merged: dict = {}
    for section in cp.sections():
        merged.setdefault(section, {}).update(dict(cp.items(section)))
    return merged


def _resolve(base: Path, raw: str) -> Path:
    """Resolve a path that may be absolute, relative, or use ~."""
    p = Path(raw).expanduser()
    return p if p.is_absolute() else (base / p).resolve()


# ──────────────────────────────────────────────────────────────────────────────
# LIGAND AUTO-DETECTION
# ──────────────────────────────────────────────────────────────────────────────
_STANDARD_RESIDUES = frozenset({
    # Amino acids (all standard + common protonation variants)
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","HIE","HID","HIP","HSE","HSD","HSP","HSD",
    "ILE","LEU","LYS","MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL",
    # Water
    "WAT","SOL","HOH","TIP","SPC","TIP3","TIP4","T3P",
    # Common ions
    "NA","CL","K","MG","CA","ZN","FE","CU","MN",
    "BR","IOD","F","LI","RB","CS","PO4","SO4",
    # Capping groups and cofactors sometimes mistaken for ligands
    "ACE","NME","NHE","FOR",
})


def _detect_ligand(gro_path: Path) -> tuple:
    """
    Parse a GRO file and return the unique non-protein/water/ion residue name.
    Returns (result, status) where status is one of:
      "found"         – exactly one candidate, result is the name string
      "multiple"      – >1 candidates, result is a sorted list
      "not_found"     – no candidates
      "gro_missing"   – GRO file does not exist
      "read_error"    – file exists but could not be parsed
    """
    if not gro_path.exists():
        return None, "gro_missing"
    try:
        with open(gro_path) as fh:
            lines = fh.readlines()
        if len(lines) < 3:
            return None, "read_error"
        # GRO format: line[0]=title, line[1]=natoms, line[2:-1]=atoms, line[-1]=box
        candidates: set = set()
        for line in lines[2:-1]:
            if len(line) >= 15:
                resname = line[5:10].strip().upper()
                if resname and not resname.isdigit() and resname not in _STANDARD_RESIDUES:
                    candidates.add(resname)
    except Exception as exc:
        return None, f"read_error:{exc}"

    if len(candidates) == 1:
        return candidates.pop(), "found"
    if len(candidates) > 1:
        return sorted(candidates), "multiple"
    return None, "not_found"


# ──────────────────────────────────────────────────────────────────────────────
# INPUT VALIDATION
# ──────────────────────────────────────────────────────────────────────────────
def _validate_inputs() -> bool:
    """
    Check GROMACS availability and all required simulation files.
    Prints specific, actionable error messages.
    Returns True if everything is in order.
    """
    ok = True

    # ── GROMACS ──────────────────────────────────────────────────────────────
    gmx_check = subprocess.run(
        [GROMACS, "--version"], capture_output=True, text=True
    )
    if gmx_check.returncode != 0:
        print(f"\n[ERROR] GROMACS command '{GROMACS}' not found.")
        print( "  Fixes:")
        print( "    HPC/cluster  →  module load gromacs")
        print( "    conda        →  conda install -c bioconda gromacs")
        print(f"    custom name  →  set  gromacs_cmd = gmx_mpi  in md_config.ini [system]")
        ok = False

    # ── Simulation files ──────────────────────────────────────────────────────
    required = {
        "TPR (topology + coordinates)": TPR,
        "XTC (trajectory)":             XTC,
        "EDR (energy)":                 EDR,
        "GRO (structure)":              GRO,
        "NDX (index)":                  NDX,
    }
    missing_files = {lbl: p for lbl, p in required.items() if not p.exists()}

    if missing_files:
        print(f"\n[ERROR] Missing files in {SIM_DIR}:")
        for lbl, p in missing_files.items():
            print(f"  {lbl}: '{p.name}'  <-- NOT FOUND")

        # Show what IS present to help the user
        found = sorted(
            list(SIM_DIR.glob("*.tpr")) + list(SIM_DIR.glob("*.xtc")) +
            list(SIM_DIR.glob("*.gro")) + list(SIM_DIR.glob("*.edr")) +
            list(SIM_DIR.glob("*.ndx"))
        )
        if found:
            print(f"\n  Files found in {SIM_DIR}:")
            for f in found:
                print(f"    {f.name}")
            print(f"\n  If your file names differ, update [paths] in md_config.ini:")
            print( "    [paths]")
            for lbl, p in missing_files.items():
                ext = p.suffix.lstrip(".")
                print(f"    {ext} = <your-actual-file-name>.{ext}")
        else:
            print(f"\n  No GROMACS files found in {SIM_DIR}")
            print( "  Make sure you copied .tpr  .xtc  .edr  .gro  .ndx")
            print(f"  into the same folder as md_analysis.py  ({SIM_DIR})")
        ok = False

    return ok


# ──────────────────────────────────────────────────────────────────────────────
# INDEX GROUP AUTO-DETECTION
# ──────────────────────────────────────────────────────────────────────────────
def _parse_ndx(ndx_path: Path) -> dict:
    """
    Parse a GROMACS index.ndx file.
    Returns {group_name: group_number_str} in declaration order (0-based).
    """
    groups = {}
    idx = -1
    with open(ndx_path) as fh:
        for line in fh:
            m = re.match(r'\[\s*(.+?)\s*\]', line)
            if m:
                idx += 1
                groups[m.group(1)] = str(idx)
    return groups


def _auto_groups(ndx_path: Path, lig_resname: str,
                 gro_path: Path, gmx_cmd: str) -> dict:
    """
    Read index.ndx and map each required group key to its group number.
    If the Protein_LIG combined group is missing, tries to create it via
    gmx make_ndx.  Prints a summary table so the user can verify.
    Returns {internal_key: group_number_str}.
    """
    if not ndx_path.exists():
        print(f"  [groups] index.ndx not found at {ndx_path} — using fallback numbers")
        return {}

    raw = _parse_ndx(ndx_path)

    # Candidate names for each required role (checked in priority order)
    candidates = {
        "System":         ["System"],
        "Protein":        ["Protein"],
        "C_alpha":        ["C-alpha", "CA", "Calpha", "C_alpha"],
        "Backbone":       ["Backbone"],
        "LIG":            [lig_resname, "LIG", "MOL", "UNK", "DRG"],
        "Water_and_ions": ["Water_and_ions", "W_and_ions", "SOL_Ion",
                           "Water&Ions", "non-Protein"],
        "Protein_LIG":    [f"Protein_{lig_resname}", "Protein_LIG",
                           f"Protein_{lig_resname.lower()}"],
    }

    grp = {}
    for key, cands in candidates.items():
        for cand in cands:
            if cand in raw:
                grp[key] = raw[cand]
                break

    # If Protein_LIG is missing, try creating it with gmx make_ndx
    if "Protein_LIG" not in grp and "Protein" in grp and "LIG" in grp:
        print(f"  [groups] Protein_LIG group not found — attempting to create it ...")
        try:
            cmd = [gmx_cmd, "make_ndx",
                   "-f", str(gro_path),
                   "-n", str(ndx_path),
                   "-o", str(ndx_path)]
            stdin = f"{grp['Protein']} | {grp['LIG']}\nq\n"
            r = subprocess.run(cmd, input=stdin, text=True,
                               capture_output=True, cwd=str(ndx_path.parent))
            if r.returncode == 0:
                raw2 = _parse_ndx(ndx_path)
                # The new group is whatever appeared last in the updated file
                new_keys = [k for k in raw2 if k not in raw]
                if new_keys:
                    grp["Protein_LIG"] = raw2[new_keys[-1]]
                    print(f"  [groups] Created '{new_keys[-1]}' "
                          f"→ group {grp['Protein_LIG']}")
                else:
                    print("  [groups] make_ndx ran but no new group was written.")
            else:
                print(f"  [groups] make_ndx failed (return code {r.returncode}); "
                      "set Protein_LIG manually in md_config.yaml → groups")
        except FileNotFoundError:
            print(f"  [groups] '{gmx_cmd}' not found on PATH; "
                  "set Protein_LIG manually in md_config.yaml → groups")
        except Exception as exc:
            print(f"  [groups] Could not create Protein_LIG: {exc}")

    # Print summary table
    print(f"\n  {'Role':<16}  {'NDX name':<22}  Group#")
    print(f"  {'-'*16}  {'-'*22}  {'------'}")
    for key, cands in candidates.items():
        num = grp.get(key, "NOT FOUND")
        # find the actual name that matched
        matched = next((c for c in cands if c in raw or
                        (key == "Protein_LIG" and num != "NOT FOUND")), "—")
        flag = "  ← WARNING: missing" if num == "NOT FOUND" else ""
        print(f"  {key:<16}  {matched:<22}  {num}{flag}")
    print()

    return grp


def _configure(cfg: dict, cli: argparse.Namespace, config_dir: Path) -> None:
    """
    Populate all global config variables from the YAML dict,
    then apply any CLI overrides on top.
    """
    global SIM_DIR, OUT_DIR, TABLES_DIR, TPR, XTC, EDR, GRO, NDX
    global XTC_NOJUMP, XTC_CENTER
    global LIG_RESNAME, GATEKEEPER_RES, SIM_LABEL, GROMACS
    global GRP, COLORS
    global PLOT_DPI, PLOT_W, PLOT_H, LIG_MARKER_SIZE
    global CONTACT_CUTOFF, HBOND_DIST, HBOND_ANGLE, MIN_CONTACT_PCT

    # ── paths ────────────────────────────────────────────────────────────────
    paths = cfg.get("paths", {})
    # Default sim_dir = "." means files live in the same folder as the config/script
    sim_dir_raw = getattr(cli, "sim_dir", None) or paths.get("sim_dir", ".")
    out_dir_raw = getattr(cli, "out_dir", None) or paths.get("out_dir", "results")

    SIM_DIR = _resolve(config_dir, sim_dir_raw)
    OUT_DIR = _resolve(config_dir, out_dir_raw)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    TABLES_DIR = OUT_DIR / "tables"
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    def _sim(name, default):
        return SIM_DIR / paths.get(name, default)

    TPR = _sim("tpr", "md.tpr")
    XTC = _sim("xtc", "md.xtc")
    EDR = _sim("edr", "md.edr")
    GRO = _sim("gro", "md.gro")
    NDX = _sim("ndx", "index.ndx")

    XTC_NOJUMP = OUT_DIR / "md_nojump.xtc"
    XTC_CENTER = OUT_DIR / "md_center.xtc"

    # ── system ───────────────────────────────────────────────────────────────
    sys_cfg = cfg.get("system", {})
    SIM_LABEL = sys_cfg.get("sim_label", "simulation")
    GROMACS   = sys_cfg.get("gromacs_cmd", "gmx")

    # ── ligand auto-detection ─────────────────────────────────────────────────
    lig_from_cfg = sys_cfg.get("lig_resname", "").strip()
    if lig_from_cfg:
        LIG_RESNAME = lig_from_cfg
        print(f"\n  [ligand] Using config value: '{LIG_RESNAME}'")
    else:
        detected, status = _detect_ligand(GRO)
        if status == "found":
            LIG_RESNAME = detected
            print(f"\n  [ligand] Auto-detected from {GRO.name}: '{LIG_RESNAME}'")
        elif status == "multiple":
            LIG_RESNAME = detected[0]
            others = ", ".join(detected[1:])
            print(f"\n  [ligand] Multiple non-protein residues found: {', '.join(detected)}")
            print(f"           Using '{LIG_RESNAME}'. If wrong, set in md_config.ini:")
            print(f"             [system]")
            print(f"             lig_resname = {others.split(',')[0].strip()}")
        elif status == "gro_missing":
            LIG_RESNAME = "LIG"   # will be caught by _validate_inputs
        else:
            LIG_RESNAME = "LIG"
            print(f"\n  [ligand] Could not auto-detect — defaulting to 'LIG'")
            print(f"           If wrong, set  lig_resname = YOUR_LIG  in md_config.ini [system]")

    # ── gatekeeper (optional) ─────────────────────────────────────────────────
    gk_raw = sys_cfg.get("gatekeeper_res", "0").strip()
    try:
        GATEKEEPER_RES = int(gk_raw) if gk_raw else 0
    except ValueError:
        GATEKEEPER_RES = 0

    if GATEKEEPER_RES == 0:
        print("  [gatekeeper] Not configured — gatekeeper analysis will be skipped")
        print("               To enable: set  gatekeeper_res = 315  in md_config.ini [system]")
    else:
        print(f"  [gatekeeper] Tracking residue {GATEKEEPER_RES}")

    # ── index groups — auto-detected, then overridden by config if present ────
    print("\n── Auto-detecting GROMACS index groups ─────────────────────────")
    auto = _auto_groups(NDX, LIG_RESNAME, GRO, GROMACS)
    GRP.clear()
    GRP.update(auto)
    overrides = {k: str(v) for k, v in cfg.get("groups", {}).items()}
    if overrides:
        GRP.update(overrides)
        print(f"  [groups] Config overrides applied: {list(overrides.keys())}")

    # ── figure / plot ─────────────────────────────────────────────────────────
    fig_cfg = cfg.get("figure", {})
    PLOT_DPI        = int(getattr(cli, "dpi",          None) or fig_cfg.get("dpi",          150))
    PLOT_W          = float(getattr(cli, "width",      None) or fig_cfg.get("width",         10))
    PLOT_H          = float(getattr(cli, "height",     None) or fig_cfg.get("height",         4))
    LIG_MARKER_SIZE = int(getattr(cli, "ligand_size", None) or fig_cfg.get("ligand_marker_size", 60))
    font_scale      = float(fig_cfg.get("font_scale", 1.1))

    COLORS.clear()
    COLORS.update({
        "backbone":   fig_cfg.get("backbone_color",   "#2c7bb6"),
        "ligand":     getattr(cli, "ligand_color", None) or fig_cfg.get("ligand_color", "#d7191c"),
        "hbond":      fig_cfg.get("hbond_color",      "#1a9641"),
        "gatekeeper": fig_cfg.get("gatekeeper_color", "#fd8d3c"),
        "rmsf":       fig_cfg.get("rmsf_color",       "#756bb1"),
        "energy":     fig_cfg.get("energy_color",     "#636363"),
    })

    sns.set_theme(style="whitegrid", palette="muted", font_scale=font_scale)

    # ── analysis params ───────────────────────────────────────────────────────
    ana_cfg = cfg.get("analysis", {})
    CONTACT_CUTOFF  = float(ana_cfg.get("contact_cutoff_A", 4.5))
    HBOND_DIST      = float(ana_cfg.get("hbond_distance_A", 3.5))
    HBOND_ANGLE     = float(ana_cfg.get("hbond_angle_deg",  120))
    MIN_CONTACT_PCT = float(ana_cfg.get("min_plot_contact_pct", 10))


# ──────────────────────────────────────────────────────────────────────────────
# UTILITIES
# ──────────────────────────────────────────────────────────────────────────────
def run_gmx(args: list, stdin: str = "", label: str = "") -> tuple:
    cmd = [GROMACS] + args
    tag = label or args[0]
    print(f"  [gmx] {' '.join(cmd[:4])} ...")
    r = subprocess.run(cmd, input=stdin, text=True,
                       capture_output=True, cwd=str(SIM_DIR))
    if r.returncode != 0:
        print(f"  [WARN] gmx {tag} exited {r.returncode}")
        for line in r.stderr.splitlines()[-6:]:
            print(f"    {line}")
    return r.stdout, r.stderr, r.returncode


def parse_xvg(path: Path) -> pd.DataFrame:
    """Read a GROMACS .xvg file; column headers extracted from @ lines."""
    headers, rows = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("@ s") and "legend" in line:
                label = line.split('"')[1]
                headers.append(label)
            elif not line.startswith(("#", "@")) and line:
                rows.append([float(x) for x in line.split()])
    arr = np.array(rows)
    cols = ["time_ps"] + (headers if len(headers) == arr.shape[1] - 1
                          else [f"col{i}" for i in range(1, arr.shape[1])])
    return pd.DataFrame(arr, columns=cols)


def save_fig(fig, name: str):
    p = OUT_DIR / name
    fig.savefig(p, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {p}")


# ──────────────────────────────────────────────────────────────────────────────
# STEP 1 — PBC CORRECTION
# ──────────────────────────────────────────────────────────────────────────────
def preprocess_trajectory() -> Path:
    """
    Fix periodic-boundary artefacts in two passes:
      pass 1: -pbc nojump removes molecules jumping across box faces
      pass 2: -center -pbc mol re-centres on Protein_LIG, wraps molecules whole
    """
    print("\n── Step 1: PBC correction ──────────────────────────────────────")

    if XTC_CENTER.exists():
        print(f"  Centred trajectory already exists: {XTC_CENTER.name}  (skipping)")
        return XTC_CENTER

    run_gmx(
        ["trjconv", "-s", str(TPR), "-f", str(XTC),
         "-o", str(XTC_NOJUMP), "-pbc", "nojump", "-n", str(NDX)],
        stdin=f"{GRP['System']}\n",
        label="trjconv-nojump",
    )

    run_gmx(
        ["trjconv", "-s", str(TPR), "-f", str(XTC_NOJUMP),
         "-o", str(XTC_CENTER), "-center", "-pbc", "mol", "-n", str(NDX)],
        stdin=f"{GRP['Protein_LIG']}\n{GRP['System']}\n",
        label="trjconv-center",
    )

    if XTC_CENTER.exists():
        print(f"  Centred trajectory written: {XTC_CENTER.name}")
    else:
        print("  [ERROR] Centred trajectory not created — check GROMACS output above")
    return XTC_CENTER


# ──────────────────────────────────────────────────────────────────────────────
# STEP 2 — THERMODYNAMIC PROPERTIES
# ──────────────────────────────────────────────────────────────────────────────
def energy_analysis() -> dict:
    print("\n── Step 2: Energy analysis ─────────────────────────────────────")

    terms = {
        "potential":    ("12", "Potential energy (kJ/mol)"),
        "kinetic":      ("13", "Kinetic energy (kJ/mol)"),
        "total_energy": ("14", "Total energy (kJ/mol)"),
        "temperature":  ("16", "Temperature (K)"),
        "pressure":     ("17", "Pressure (bar)"),
        "volume":       ("22", "Volume (nm³)"),
        "density":      ("23", "Density (kg/m³)"),
    }

    data = {}
    for name, (idx, ylabel) in terms.items():
        xvg = TABLES_DIR / f"energy_{name}.xvg"
        if not xvg.exists():
            run_gmx(
                ["energy", "-f", str(EDR), "-o", str(xvg)],
                stdin=f"{idx}\n0\n",
                label=f"energy-{name}",
            )
        if xvg.exists():
            df = parse_xvg(xvg)
            df.columns = ["time_ps", name]
            data[name] = df
            mu  = df[name].mean()
            std = df[name].std()
            print(f"  {name:15s}: {mu:>12.3f} ± {std:.3f}  ({ylabel.split('(')[1][:-1]})")

    # Save each thermodynamic property as its own figure
    for name, (_, ylabel) in terms.items():
        if name not in data:
            continue
        df  = data[name]
        mu  = df[name].mean()
        std = df[name].std()

        fig, ax = plt.subplots(figsize=(PLOT_W, PLOT_H))
        ax.plot(df["time_ps"], df[name], color=COLORS["energy"], lw=1.0)
        ax.fill_between(df["time_ps"], mu - std, mu + std,
                        color="red", alpha=0.08, label=f"±1σ  ({std:.3f})")
        ax.axhline(mu, color="red", ls="--", lw=1.0,
                   label=f"mean = {mu:.3f}")
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{name.replace('_', ' ').title()}")
        ax.legend(fontsize=9)
        plt.tight_layout()
        save_fig(fig, f"energy_{name}.png")

    return data


# ──────────────────────────────────────────────────────────────────────────────
# STEP 3 — RMSD
# ──────────────────────────────────────────────────────────────────────────────
def rmsd_analysis(traj: Path) -> pd.DataFrame:
    print("\n── Step 3: RMSD analysis ───────────────────────────────────────")
    if not HAS_MDA:
        print("  [SKIP] MDAnalysis not available"); return pd.DataFrame()

    u   = mda.Universe(str(GRO), str(traj))
    ref = mda.Universe(str(GRO), str(traj))

    ca_sel  = "backbone and name CA"
    lig_sel = f"resname {LIG_RESNAME} and not name H*"

    R = rms.RMSD(u, ref, select=ca_sel, groupselections=[lig_sel], ref_frame=0)
    R.run()

    rmsd_arr = R.results.rmsd
    df = pd.DataFrame({
        "time_ps":       rmsd_arr[:, 1],
        "backbone_RMSD": rmsd_arr[:, 2],
        "ligand_RMSD":   rmsd_arr[:, 3],
    })

    for col in ("backbone_RMSD", "ligand_RMSD"):
        mu, std = df[col].mean(), df[col].std()
        print(f"  {col}: {mu:.3f} ± {std:.3f} Å")

    df.to_csv(TABLES_DIR / "rmsd.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(PLOT_W, PLOT_H * 1.75), sharex=True)

    axes[0].plot(df["time_ps"], df["backbone_RMSD"],
                 color=COLORS["backbone"], lw=1.2, label="Backbone Cα")
    axes[0].set_ylabel("RMSD (Å)")
    axes[0].set_title("Protein Backbone RMSD")
    axes[0].legend()
    axes[0].axhline(df["backbone_RMSD"].mean(), color="k", ls="--", lw=0.8)

    axes[1].plot(df["time_ps"], df["ligand_RMSD"],
                 color=COLORS["ligand"], lw=1.2, label=f"Ligand ({LIG_RESNAME})")
    axes[1].set_xlabel("Time (ps)")
    axes[1].set_ylabel("RMSD (Å)")
    axes[1].set_title("Ligand RMSD (fit on backbone)")
    axes[1].legend()
    axes[1].axhline(df["ligand_RMSD"].mean(), color="k", ls="--", lw=0.8)

    fig.suptitle(f"RMSD", fontsize=13)
    plt.tight_layout()
    save_fig(fig, "rmsd.png")
    return df


# ──────────────────────────────────────────────────────────────────────────────
# STEP 4 — RMSF per Cα residue
# ──────────────────────────────────────────────────────────────────────────────
def rmsf_analysis(traj: Path) -> pd.DataFrame:
    print("\n── Step 4: RMSF analysis ───────────────────────────────────────")
    if not HAS_MDA:
        print("  [SKIP] MDAnalysis not available"); return pd.DataFrame()

    u = mda.Universe(str(GRO), str(traj))

    aligner = align.AlignTraj(u, u, select="backbone and name CA", in_memory=True)
    aligner.run()

    ca = u.select_atoms("backbone and name CA")
    rmsf_calc = rms.RMSF(ca)
    rmsf_calc.run()

    df = pd.DataFrame({
        "resid":   ca.resids,
        "resname": ca.resnames,
        "rmsf_A":  rmsf_calc.results.rmsf,
    })
    df.to_csv(TABLES_DIR / "rmsf.csv", index=False)

    top5 = df.nlargest(5, "rmsf_A")[["resid", "resname", "rmsf_A"]]
    print("  Top-5 most flexible residues:")
    print(textwrap.indent(top5.to_string(index=False), "    "))

    fig, ax = plt.subplots(figsize=(PLOT_W * 1.4, PLOT_H))
    ax.fill_between(df["resid"], df["rmsf_A"], alpha=0.35, color=COLORS["rmsf"])
    ax.plot(df["resid"], df["rmsf_A"], color=COLORS["rmsf"], lw=0.9)

    gk_row = df[df["resid"] == GATEKEEPER_RES]
    if not gk_row.empty:
        gk_val = gk_row["rmsf_A"].values[0]
        ax.axvline(GATEKEEPER_RES, color="red", ls="--", lw=1.0,
                   label=f"Gatekeeper ({gk_row['resname'].values[0]}{GATEKEEPER_RES})")
        ax.scatter([GATEKEEPER_RES], [gk_val],
                   color=COLORS["ligand"], zorder=5, s=LIG_MARKER_SIZE)
        ax.legend()

    ax.set_xlabel("Residue ID")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title(f"Per-residue Cα RMSF")
    plt.tight_layout()
    save_fig(fig, "rmsf.png")
    return df


# ──────────────────────────────────────────────────────────────────────────────
# STEP 5 — RADIUS OF GYRATION
# ──────────────────────────────────────────────────────────────────────────────
def gyration_analysis() -> pd.DataFrame:
    print("\n── Step 5: Radius of gyration ──────────────────────────────────")
    xvg = TABLES_DIR / "gyrate.xvg"

    if not xvg.exists():
        run_gmx(
            ["gyrate", "-s", str(TPR), "-f", str(XTC_CENTER),
             "-o", str(xvg), "-n", str(NDX)],
            stdin=f"{GRP['Protein']}\n",
            label="gyrate",
        )

    if not xvg.exists():
        print("  [SKIP] gyrate.xvg not produced"); return pd.DataFrame()

    df = parse_xvg(xvg)
    df.columns = ["time_ps", "Rg", "Rg_X", "Rg_Y", "Rg_Z"][:df.shape[1]]
    print(f"  Mean Rg: {df['Rg'].mean():.3f} ± {df['Rg'].std():.3f} nm")

    fig, ax = plt.subplots(figsize=(PLOT_W, PLOT_H))
    ax.plot(df["time_ps"], df["Rg"] * 10, color=COLORS["energy"], lw=1.2)
    ax.axhline(df["Rg"].mean() * 10, color="k", ls="--", lw=0.8,
               label=f"mean = {df['Rg'].mean()*10:.2f} Å")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Radius of Gyration (Å)")
    ax.set_title(f"Protein Radius of Gyration")
    ax.legend()
    plt.tight_layout()
    save_fig(fig, "gyration.png")
    return df


# ──────────────────────────────────────────────────────────────────────────────
# STEP 6 — H-BOND ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────
def hbond_analysis(traj: Path) -> pd.DataFrame:
    print("\n── Step 6: H-bond analysis ─────────────────────────────────────")
    if not HAS_MDA:
        print("  [SKIP] MDAnalysis not available"); return pd.DataFrame()

    u = mda.Universe(str(GRO), str(traj))

    prot_sel  = "protein"
    lig_sel2  = f"resname {LIG_RESNAME}"
    donors_sel    = (f"({prot_sel} and (name N O)) or "
                     f"({lig_sel2} and (name N* O*))")
    acceptors_sel = (f"({prot_sel} and (name O N)) or "
                     f"({lig_sel2} and (name N* O*))")

    hba = HBA(
        universe=u,
        donors_sel=donors_sel,
        hydrogens_sel="name H*",
        acceptors_sel=acceptors_sel,
        between=[lig_sel2, prot_sel],
        d_a_cutoff=HBOND_DIST,
        d_h_a_angle_cutoff=HBOND_ANGLE,
        update_selections=False,
    )
    hba.run()

    counts = hba.count_by_time()
    u2 = mda.Universe(str(GRO), str(traj))
    times = [ts.time for ts in u2.trajectory]
    time_arr = np.array(times) if len(times) == len(counts) else np.arange(len(counts))

    df_count = pd.DataFrame({"time_ps": time_arr, "n_hbonds": counts})
    df_count.to_csv(TABLES_DIR / "hbonds_count.csv", index=False)
    print(f"  Mean H-bonds/frame: {df_count['n_hbonds'].mean():.2f} "
          f"± {df_count['n_hbonds'].std():.2f}")

    pair_occ = pd.DataFrame()
    if len(hba.results.hbonds) > 0:
        hb_df = pd.DataFrame(
            hba.results.hbonds,
            columns=["frame", "donor_ix", "H_ix", "acceptor_ix",
                     "distance_A", "angle_deg"],
        )
        all_atoms = u.atoms
        hb_df["donor_resid"]   = hb_df["donor_ix"].astype(int).map(lambda ix: all_atoms[ix].resid)
        hb_df["donor_resname"] = hb_df["donor_ix"].astype(int).map(lambda ix: all_atoms[ix].resname)
        hb_df["donor_name"]    = hb_df["donor_ix"].astype(int).map(lambda ix: all_atoms[ix].name)
        hb_df["acc_resid"]     = hb_df["acceptor_ix"].astype(int).map(lambda ix: all_atoms[ix].resid)
        hb_df["acc_resname"]   = hb_df["acceptor_ix"].astype(int).map(lambda ix: all_atoms[ix].resname)
        hb_df["acc_name"]      = hb_df["acceptor_ix"].astype(int).map(lambda ix: all_atoms[ix].name)

        n_frames = len(u.trajectory)
        pair_occ = (
            hb_df.groupby(["donor_resid", "donor_resname", "donor_name",
                            "acc_resid", "acc_resname", "acc_name"])
            .size()
            .reset_index(name="n_frames")
        )
        pair_occ["occupancy_%"] = (pair_occ["n_frames"] / n_frames * 100).round(1)
        pair_occ = pair_occ.sort_values("occupancy_%", ascending=False)
        pair_occ.to_csv(TABLES_DIR / "hbonds_pairs.csv", index=False)

        print("  Top H-bond pairs (occupancy %):")
        print(textwrap.indent(pair_occ.head(10).to_string(index=False), "    "))

        if GATEKEEPER_RES > 0:
            gk_hb = pair_occ[
                (pair_occ["donor_resid"] == GATEKEEPER_RES) |
                (pair_occ["acc_resid"]   == GATEKEEPER_RES)
            ]
            if gk_hb.empty:
                print(f"  [NOTE] No H-bond with gatekeeper residue {GATEKEEPER_RES}")
            else:
                print(f"\n  Gatekeeper {GATEKEEPER_RES} H-bond occupancy:")
                print(textwrap.indent(gk_hb.to_string(index=False), "    "))
    else:
        print("  [NOTE] No protein↔ligand H-bonds detected in this trajectory")

    fig, ax = plt.subplots(figsize=(PLOT_W, PLOT_H))
    ax.fill_between(df_count["time_ps"], df_count["n_hbonds"],
                    alpha=0.4, color=COLORS["hbond"])
    ax.plot(df_count["time_ps"], df_count["n_hbonds"],
            color=COLORS["hbond"], lw=1.0)
    ax.axhline(df_count["n_hbonds"].mean(), color="k", ls="--", lw=0.8,
               label=f"mean = {df_count['n_hbonds'].mean():.2f}")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Number of H-bonds")
    ax.set_title(f"Protein–Ligand H-bonds")
    ax.legend()
    plt.tight_layout()
    save_fig(fig, "hbonds.png")

    return pair_occ


# ──────────────────────────────────────────────────────────────────────────────
# STEP 7 — GATEKEEPER DISTANCE
# ──────────────────────────────────────────────────────────────────────────────
def gatekeeper_distance(traj: Path) -> pd.DataFrame:
    print("\n── Step 7: Gatekeeper distance ─────────────────────────────────")
    if GATEKEEPER_RES == 0:
        print("  [SKIP] gatekeeper_res not set in md_config.ini")
        print("         Set  gatekeeper_res = 315  (or your residue number) to enable")
        return pd.DataFrame()
    if not HAS_MDA:
        print("  [SKIP] MDAnalysis not available"); return pd.DataFrame()

    u         = mda.Universe(str(GRO), str(traj))
    lig_atoms = u.select_atoms(f"resname {LIG_RESNAME} and not name H*")
    gk_atoms  = u.select_atoms(f"resid {GATEKEEPER_RES} and not name H*")

    if gk_atoms.n_atoms == 0:
        print(f"  [WARN] Residue {GATEKEEPER_RES} not found — check gatekeeper_res in config")
        return pd.DataFrame()

    gk_resname = gk_atoms.resnames[0]
    print(f"  Gatekeeper: {gk_resname}{GATEKEEPER_RES} ({gk_atoms.n_atoms} heavy atoms)")

    from MDAnalysis.lib.distances import distance_array

    times, dists = [], []
    for ts in u.trajectory:
        d = distance_array(lig_atoms.positions, gk_atoms.positions, box=ts.dimensions)
        times.append(ts.time)
        dists.append(d.min())

    df = pd.DataFrame({"time_ps": times, "min_dist_A": dists})
    df.to_csv(TABLES_DIR / "gatekeeper_distance.csv", index=False)
    print(f"  Mean min-distance to {gk_resname}{GATEKEEPER_RES}: "
          f"{df['min_dist_A'].mean():.2f} ± {df['min_dist_A'].std():.2f} Å")

    close_frac = (df["min_dist_A"] < HBOND_DIST).mean() * 100
    print(f"  Frames < {HBOND_DIST} Å: {close_frac:.1f}%")

    fig, ax = plt.subplots(figsize=(PLOT_W, PLOT_H))
    ax.fill_between(df["time_ps"], df["min_dist_A"],
                    alpha=0.3, color=COLORS["gatekeeper"])
    ax.plot(df["time_ps"], df["min_dist_A"], color=COLORS["gatekeeper"], lw=1.2)
    ax.axhline(HBOND_DIST, color="k", ls="--", lw=0.9,
               label=f"H-bond cutoff ({HBOND_DIST} Å)")
    ax.axhline(df["min_dist_A"].mean(), color="red", ls=":", lw=0.9,
               label=f"mean = {df['min_dist_A'].mean():.2f} Å")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Min distance (Å)")
    ax.set_title(
        f"Ligand – {gk_resname}{GATEKEEPER_RES} (gatekeeper) Min Distance")
    ax.legend()
    plt.tight_layout()
    save_fig(fig, "gatekeeper_distance.png")
    return df


# ──────────────────────────────────────────────────────────────────────────────
# STEP 8 — BINDING-POCKET CONTACT MAP
# ──────────────────────────────────────────────────────────────────────────────
def pocket_contacts(traj: Path) -> pd.DataFrame:
    cutoff = CONTACT_CUTOFF
    print(f"\n── Step 8: Pocket contact map (cutoff {cutoff} Å) ──────────────")
    if not HAS_MDA:
        print("  [SKIP] MDAnalysis not available"); return pd.DataFrame()

    from MDAnalysis.lib.distances import distance_array

    u    = mda.Universe(str(GRO), str(traj))
    lig  = u.select_atoms(f"resname {LIG_RESNAME} and not name H*")
    prot = u.select_atoms("protein and not name H*")
    resids   = np.unique(prot.resids)
    resnames = {r: prot.select_atoms(f"resid {r}").resnames[0] for r in resids}

    contact_count = {r: 0 for r in resids}
    n_frames = 0

    for ts in u.trajectory:
        d = distance_array(lig.positions, prot.positions, box=ts.dimensions)
        in_contact = d.min(axis=0) < cutoff
        for r in resids:
            mask = prot.resids == r
            if in_contact[mask].any():
                contact_count[r] += 1
        n_frames += 1

    df = pd.DataFrame({
        "resid":        list(contact_count.keys()),
        "resname":      [resnames[r] for r in contact_count],
        "contact_frac": [v / n_frames for v in contact_count.values()],
    })
    df = df[df["contact_frac"] > 0].sort_values("contact_frac", ascending=False)
    df["contact_%"] = (df["contact_frac"] * 100).round(1)
    df.to_csv(TABLES_DIR / "pocket_contacts.csv", index=False)

    print(f"  Residues in contact (> 0 frames): {len(df)}")
    print(f"  Residues with > 80% occupancy:    {(df['contact_%'] > 80).sum()}")
    print("\n  Top-15 contact residues:")
    print(textwrap.indent(
        df.head(15)[["resid", "resname", "contact_%"]].to_string(index=False), "    "))

    df_plot = df[df["contact_%"] >= MIN_CONTACT_PCT].copy()
    df_plot["label"] = df_plot["resname"] + df_plot["resid"].astype(str)

    bar_colors = [
        COLORS["gatekeeper"] if r == GATEKEEPER_RES else COLORS["backbone"]
        for r in df_plot["resid"]
    ]

    fig, ax = plt.subplots(figsize=(max(PLOT_W, len(df_plot) * 0.45), PLOT_H * 1.25))
    ax.bar(df_plot["label"], df_plot["contact_%"],
           color=bar_colors, edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Residue")
    ax.set_ylabel("Contact occupancy (%)")
    ax.set_title(
        f"Binding-pocket Contacts (≥ {MIN_CONTACT_PCT}% occupancy)")
    ax.set_ylim(0, 105)
    plt.xticks(rotation=45, ha="right", fontsize=8)

    from matplotlib.patches import Patch
    legend_els = [Patch(facecolor=COLORS["backbone"], label="Pocket residue")]
    if GATEKEEPER_RES > 0:
        legend_els.append(
            Patch(facecolor=COLORS["gatekeeper"], label=f"Gatekeeper ({GATEKEEPER_RES})")
        )
    ax.legend(handles=legend_els, loc="upper right")
    plt.tight_layout()
    save_fig(fig, "pocket_contacts.png")
    return df


# ──────────────────────────────────────────────────────────────────────────────
# STEP 9 — SUMMARY PANEL FIGURE
# ──────────────────────────────────────────────────────────────────────────────
def summary_figure(rmsd_df: pd.DataFrame, rg_df: pd.DataFrame,
                   rmsf_df: pd.DataFrame) -> None:
    print("\n── Step 9: Summary figure ──────────────────────────────────────")

    fig = plt.figure(figsize=(PLOT_W * 1.4, PLOT_H * 2.5))
    gs  = gridspec.GridSpec(2, 2, hspace=0.4, wspace=0.35)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    if not rmsd_df.empty and "backbone_RMSD" in rmsd_df.columns:
        ax1.plot(rmsd_df["time_ps"], rmsd_df["backbone_RMSD"],
                 color=COLORS["backbone"], lw=1.2)
        ax1.axhline(rmsd_df["backbone_RMSD"].mean(), color="k", ls="--", lw=0.8)
        ax1.set_xlabel("Time (ps)")
        ax1.set_ylabel("RMSD (Å)")
        ax1.set_title("Backbone Cα RMSD")
    else:
        ax1.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax1.transAxes)

    if not rmsd_df.empty and "ligand_RMSD" in rmsd_df.columns:
        ax2.plot(rmsd_df["time_ps"], rmsd_df["ligand_RMSD"],
                 color=COLORS["ligand"], lw=1.2)
        ax2.axhline(rmsd_df["ligand_RMSD"].mean(), color="k", ls="--", lw=0.8)
        ax2.set_xlabel("Time (ps)")
        ax2.set_ylabel("RMSD (Å)")
        ax2.set_title(f"Ligand ({LIG_RESNAME}) RMSD")
    else:
        ax2.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax2.transAxes)

    if not rg_df.empty:
        ax3.plot(rg_df["time_ps"], rg_df["Rg"] * 10,
                 color=COLORS["energy"], lw=1.2)
        ax3.axhline(rg_df["Rg"].mean() * 10, color="k", ls="--", lw=0.8,
                    label=f"mean = {rg_df['Rg'].mean()*10:.2f} Å")
        ax3.set_xlabel("Time (ps)")
        ax3.set_ylabel("Radius of Gyration (Å)")
        ax3.set_title("Protein Radius of Gyration")
        ax3.legend(fontsize=8)
    else:
        ax3.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax3.transAxes)

    if not rmsf_df.empty:
        ax4.fill_between(rmsf_df["resid"], rmsf_df["rmsf_A"],
                         alpha=0.35, color=COLORS["rmsf"])
        ax4.plot(rmsf_df["resid"], rmsf_df["rmsf_A"], color=COLORS["rmsf"], lw=0.9)
        ax4.set_xlabel("Residue ID")
        ax4.set_ylabel("RMSF (Å)")
        ax4.set_title("Per-residue Cα RMSF")
    else:
        ax4.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax4.transAxes)

    fig.suptitle(f"MD Analysis Summary", fontsize=14, fontweight="bold")
    save_fig(fig, "summary_panel.png")


# ──────────────────────────────────────────────────────────────────────────────
# STEP 10 — TEXT REPORT
# ──────────────────────────────────────────────────────────────────────────────
def write_report(energy_data: dict, rmsd_df: pd.DataFrame,
                 rmsf_df: pd.DataFrame, gk_df: pd.DataFrame,
                 contacts_df: pd.DataFrame) -> None:
    print("\n── Step 10: Writing summary report ─────────────────────────────")

    lines = [
        "=" * 70,
        f"  MD ANALYSIS REPORT",
        "=" * 70,
        "",
        "THERMODYNAMIC AVERAGES",
    ]

    labels = {
        "potential":    ("Potential energy",    "kJ/mol"),
        "kinetic":      ("Kinetic energy",      "kJ/mol"),
        "total_energy": ("Total energy",        "kJ/mol"),
        "temperature":  ("Temperature",         "K"),
        "pressure":     ("Pressure",            "bar"),
        "volume":       ("Volume",              "nm³"),
        "density":      ("Density",             "kg/m³"),
    }
    for key, (name, unit) in labels.items():
        if key in energy_data:
            vals = energy_data[key].iloc[:, 1]
            lines.append(
                f"  {name:<22}: {vals.mean():>12.3f} ± {vals.std():.3f} {unit}")

    lines += ["", "─" * 70, "STRUCTURAL STABILITY"]
    if not rmsd_df.empty:
        bk = rmsd_df["backbone_RMSD"]
        lg = rmsd_df["ligand_RMSD"]
        lines += [
            f"  Backbone Cα RMSD  : {bk.mean():.3f} ± {bk.std():.3f} Å  (max {bk.max():.3f} Å)",
            f"  Ligand RMSD       : {lg.mean():.3f} ± {lg.std():.3f} Å  (max {lg.max():.3f} Å)",
        ]
        if lg.max() > 3.0:
            lines.append("  [FLAG] Ligand RMSD > 3 Å — possible pose instability")
        if bk.max() > 3.0:
            lines.append("  [FLAG] Backbone RMSD > 3 Å — consider longer equilibration")

    if not rmsf_df.empty:
        lines += ["", "  Top-5 flexible residues (Cα RMSF):"]
        for _, row in rmsf_df.nlargest(5, "rmsf_A").iterrows():
            lines.append(f"    {row['resname']}{int(row['resid']):<6}: {row['rmsf_A']:.3f} Å")

    if GATEKEEPER_RES > 0:
        lines += ["", "─" * 70, "GATEKEEPER CONTACT"]
    if GATEKEEPER_RES > 0 and not gk_df.empty:
        mu    = gk_df["min_dist_A"].mean()
        mn    = gk_df["min_dist_A"].min()
        close = (gk_df["min_dist_A"] < HBOND_DIST).mean() * 100
        lines += [
            f"  Gatekeeper residue: {GATEKEEPER_RES}",
            f"  Mean min-distance : {mu:.2f} Å  (min = {mn:.2f} Å)",
            f"  Frames < {HBOND_DIST} Å   : {close:.1f}%  (H-bond contact threshold)",
        ]
        if close > 50:
            lines.append("  [OK]  Gatekeeper contact maintained (> 50% occupancy)")
        else:
            lines.append("  [NOTE] Gatekeeper contact < 50% — check H-bond stability")

    lines += ["", "─" * 70, "BINDING-POCKET CONTACTS"]
    if not contacts_df.empty:
        high = contacts_df[contacts_df["contact_%"] >= 80]
        lines += [
            f"  Total residues in contact  : {len(contacts_df)}",
            f"  Residues ≥ 80% occupancy   : {len(high)}",
            "",
            "  High-occupancy contacts (≥ 80%):",
        ]
        for _, r in high.iterrows():
            lines.append(f"    {r['resname']}{int(r['resid']):<6}: {r['contact_%']:.1f}%")

    lines += [
        "",
        "─" * 70,
        "OUTPUT FILES",
        f"  Figures  →  {OUT_DIR}/",
        f"    energy_potential.png   energy_kinetic.png   energy_total_energy.png",
        f"    energy_temperature.png energy_pressure.png  energy_volume.png",
        f"    energy_density.png     rmsd.png             rmsf.png",
        f"    gyration.png           hbonds.png           gatekeeper_distance.png",
        f"    pocket_contacts.png    summary_panel.png",
        f"",
        f"  Tables   →  {TABLES_DIR}/",
        f"    rmsd.csv  rmsf.csv  hbonds_count.csv  hbonds_pairs.csv",
        f"    gatekeeper_distance.csv  pocket_contacts.csv",
        f"    analysis_report.txt",
        "",
        "─" * 70,
        "NEXT STEPS",
        "  1. Extend to 50–100 ns production run for statistical validity.",
        "  2. Repeat for T315I mutant — compare ligand RMSD & gatekeeper contact loss.",
        "  3. Run gmx_MMPBSA on the converged portion for ΔG_bind.",
        "  4. Add replicas; report mean ± SD across replicates.",
        "=" * 70,
    ]

    report_path = TABLES_DIR / "analysis_report.txt"
    report_path.write_text("\n".join(lines))
    print(f"  → {report_path}")
    print()
    print("\n".join(lines))


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="md_analysis.py",
        description=(
            "Automated GROMACS MD analysis pipeline.\n"
            "Settings are read from md_config.ini; all can be overridden here."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            examples:
              python md_analysis.py
              python md_analysis.py --config /data/project/md_config.ini
              python md_analysis.py --sim-dir /data/my_sim --out-dir /data/results
              python md_analysis.py --dpi 300 --width 14 --height 5
              python md_analysis.py --ligand-color "#ff6600" --ligand-size 80
        """),
    )

    p.add_argument(
        "--config", metavar="FILE",
        help=(
            "Path to config file (.ini).  "
            "Default: md_config.ini next to this script."
        ),
    )

    grp_paths = p.add_argument_group("path overrides")
    grp_paths.add_argument(
        "--sim-dir", metavar="DIR",
        help="Directory containing GROMACS output files (.tpr, .xtc, .edr, …).",
    )
    grp_paths.add_argument(
        "--out-dir", metavar="DIR",
        help="Directory for analysis output (created if absent). Default: results/",
    )

    grp_fig = p.add_argument_group("figure overrides")
    grp_fig.add_argument(
        "--dpi", type=int, metavar="N",
        help="Image resolution in DPI (72=screen, 150=default, 300=publication).",
    )
    grp_fig.add_argument(
        "--width", type=float, metavar="INCHES",
        help="Base figure width in inches (default: 10).",
    )
    grp_fig.add_argument(
        "--height", type=float, metavar="INCHES",
        help="Base figure height in inches (default: 4).",
    )
    grp_fig.add_argument(
        "--ligand-color", metavar="COLOR",
        help="Colour for ligand traces, e.g. \"#d7191c\" or \"red\".",
    )
    grp_fig.add_argument(
        "--ligand-size", type=int, metavar="PT2",
        help="Marker size (pt²) for ligand scatter points (default: 60).",
    )

    return p


def main():
    parser = _build_parser()
    args   = parser.parse_args()

    # Locate config file — prefer .ini, fall back to legacy .yaml if present
    script_dir = Path(__file__).resolve().parent
    if args.config:
        config_path = Path(args.config).resolve()
    else:
        config_path = next(
            (script_dir / f for f in ("md_config.ini", "md_config.yaml")
             if (script_dir / f).exists()),
            script_dir / "md_config.ini",   # will trigger a clear error below
        )

    if config_path.exists():
        print(f"[config] Loading: {config_path}")
        cfg = _load_config(config_path)
        config_dir = config_path.parent
    else:
        if args.config:
            print(f"[ERROR] Config file not found: {config_path}")
            sys.exit(1)
        print("[config] md_config.ini not found — using built-in defaults")
        cfg        = {}
        config_dir = script_dir

    _configure(cfg, args, config_dir)

    print("=" * 70)
    print(f"  MD Analysis Pipeline")
    print(f"  Ligand   : {LIG_RESNAME}")
    print(f"  Gatekeeper: {'residue ' + str(GATEKEEPER_RES) if GATEKEEPER_RES else 'not set (skipped)'}")
    print(f"  Figures  : {OUT_DIR}/")
    print(f"  Tables   : {TABLES_DIR}/")
    print("=" * 70)

    # ── Validate inputs before doing anything ────────────────────────────────
    if not _validate_inputs():
        sys.exit(1)

    if not HAS_MDA:
        print("\n[WARNING] MDAnalysis not installed — structural analyses will be skipped.")
        print("  Fix:  python install.py    or    pip install MDAnalysis\n")

    # ── Run pipeline — each step is independent; failures are non-fatal ──────
    traj = None
    try:
        traj = preprocess_trajectory()
        if traj is None or not traj.exists():
            raise RuntimeError("PBC-corrected trajectory was not produced")
    except Exception as exc:
        print(f"\n[ERROR] Step 1 (PBC correction) failed: {exc}")
        print("  This step requires GROMACS. Check that gmx is on your PATH.")
        print("  Structural analyses that need the corrected trajectory will be skipped.\n")

    energy_d = {}
    try:
        energy_d = energy_analysis()
    except Exception as exc:
        print(f"  [WARN] Energy analysis failed: {exc}")
        print("         Check that the EDR file exists and is not corrupted.\n")

    rmsd_df = pd.DataFrame()
    if traj and traj.exists():
        try:
            rmsd_df = rmsd_analysis(traj)
        except Exception as exc:
            print(f"  [WARN] RMSD analysis failed: {exc}")
            print(f"         Check that '{LIG_RESNAME}' exists in your topology.\n")

    rmsf_df = pd.DataFrame()
    if traj and traj.exists():
        try:
            rmsf_df = rmsf_analysis(traj)
        except Exception as exc:
            print(f"  [WARN] RMSF analysis failed: {exc}\n")

    rg_df = pd.DataFrame()
    try:
        rg_df = gyration_analysis()
    except Exception as exc:
        print(f"  [WARN] Gyration analysis failed: {exc}\n")

    contacts_df = pd.DataFrame()
    if traj and traj.exists():
        try:
            contacts_df = pocket_contacts(traj)
        except Exception as exc:
            print(f"  [WARN] Contact map failed: {exc}\n")

    if traj and traj.exists():
        try:
            hbond_analysis(traj)
        except Exception as exc:
            print(f"  [WARN] H-bond analysis failed: {exc}\n")

    gk_df = pd.DataFrame()
    if traj and traj.exists():
        try:
            gk_df = gatekeeper_distance(traj)
        except Exception as exc:
            print(f"  [WARN] Gatekeeper distance failed: {exc}\n")

    try:
        summary_figure(rmsd_df, rg_df, rmsf_df)
    except Exception as exc:
        print(f"  [WARN] Summary figure failed: {exc}\n")

    try:
        write_report(energy_d, rmsd_df, rmsf_df, gk_df, contacts_df)
    except Exception as exc:
        print(f"  [WARN] Report generation failed: {exc}\n")

    print(f"\n[DONE] Figures  →  {OUT_DIR}/")
    print(f"       Tables   →  {TABLES_DIR}/\n")


if __name__ == "__main__":
    main()
