#!/usr/bin/env python3
# md_analysis.py — GROMACS MD analysis pipeline, fully configured via md_config.ini
# Usage:  python md_analysis.py /path/to/sim_folder
#         python md_analysis.py --config /path/to/sim_folder/md_config.ini

import sys, subprocess

def _bootstrap():
    import importlib.util
    needed = [("numpy","numpy"),("pandas","pandas"),("matplotlib","matplotlib"),
              ("seaborn","seaborn"),("MDAnalysis","MDAnalysis")]
    miss = [p for m,p in needed if importlib.util.find_spec(m) is None]
    if not miss:
        return
    print(f"[bootstrap] Installing: {', '.join(miss)} ...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "--quiet"] + miss)
    except subprocess.CalledProcessError:
        print(f"Install failed. Run:  pip install {' '.join(miss)}")
        sys.exit(1)
    sys.exit(subprocess.run([sys.executable] + sys.argv).returncode)

_bootstrap()

import argparse, re, textwrap, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams['axes.xmargin'] = 0
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import configparser
warnings.filterwarnings("ignore")

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA
    HAS_MDA = True
except ImportError:
    HAS_MDA = False
    print("WARNING: MDAnalysis not available — structural analyses skipped.\n")


# Globals — populated by _configure() before analysis
SIM_DIR = OUT_DIR = TABLES_DIR = TPR = XTC = EDR = GRO = NDX = None
XTC_NOJUMP = XTC_CENTER = None
LIG_RESNAME    = "LIG"
GATEKEEPER_RES = 315
SIM_LABEL      = "simulation"
GROMACS        = "gmx"
GRP            = {}
COLORS         = {}
PLOT_DPI        = 150
PLOT_W          = 10
PLOT_H          = 4
LIG_MARKER_SIZE = 60
TITLE_FONTSIZE  = 13
LABEL_FONTSIZE  = 11
CONTACT_CUTOFF  = 4.5
HBOND_DIST      = 3.5
HBOND_ANGLE     = 120
MIN_CONTACT_PCT = 10
TIME_UNIT  = "ps"
TIME_DIV   = 1.0
TIME_LABEL = "Time (ps)"


def t_axis(s):
    return s / TIME_DIV


_STANDARD_RESIDUES = frozenset({
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","HIE","HID","HIP","HSE","HSD","HSP",
    "ILE","LEU","LYS","MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL",
    "WAT","SOL","HOH","TIP","SPC","TIP3","TIP4","T3P",
    "NA","CL","K","MG","CA","ZN","FE","CU","MN",
    "BR","IOD","F","LI","RB","CS","PO4","SO4",
    "ACE","NME","NHE","FOR",
})


def _detect_ligand(gro_path):
    if not gro_path.exists():
        return None, "gro_missing"
    try:
        lines = gro_path.read_text().splitlines()
        if len(lines) < 3:
            return None, "read_error"
        candidates = set()
        for line in lines[2:-1]:
            if len(line) >= 15:
                resname = line[5:10].strip().upper()
                if resname and not resname.isdigit() and resname not in _STANDARD_RESIDUES:
                    candidates.add(resname)
    except Exception as e:
        return None, f"read_error:{e}"
    if len(candidates) == 1:
        return candidates.pop(), "found"
    if len(candidates) > 1:
        return sorted(candidates), "multiple"
    return None, "not_found"


def _load_config(path):
    cp = configparser.RawConfigParser(
        comment_prefixes=("#", ";"),
        inline_comment_prefixes=None,
        strict=False,
    )
    cp.read(path, encoding="utf-8")
    merged = {}
    for section in cp.sections():
        merged.setdefault(section, {}).update(dict(cp.items(section)))
    return merged


def _resolve(base, raw):
    p = Path(raw).expanduser()
    return p if p.is_absolute() else (base / p).resolve()


def _parse_ndx(ndx_path):
    groups, idx = {}, -1
    for line in ndx_path.read_text().splitlines():
        m = re.match(r'\[\s*(.+?)\s*\]', line)
        if m:
            idx += 1
            groups[m.group(1)] = str(idx)
    return groups


def _auto_groups(ndx_path, lig_resname, gro_path, gmx_cmd):
    if not ndx_path.exists():
        print("  [groups] index.ndx not found — using fallback numbers")
        return {}

    raw = _parse_ndx(ndx_path)
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
        for c in cands:
            if c in raw:
                grp[key] = raw[c]
                break

    if "Protein_LIG" not in grp and "Protein" in grp and "LIG" in grp:
        print("  [groups] Protein_LIG missing — attempting to create it ...")
        try:
            r = subprocess.run(
                [gmx_cmd, "make_ndx", "-f", str(gro_path),
                 "-n", str(ndx_path), "-o", str(ndx_path)],
                input=f"{grp['Protein']} | {grp['LIG']}\nq\n",
                text=True, capture_output=True, cwd=str(ndx_path.parent),
            )
            if r.returncode == 0:
                raw2 = _parse_ndx(ndx_path)
                new = [k for k in raw2 if k not in raw]
                if new:
                    grp["Protein_LIG"] = raw2[new[-1]]
                    print(f"  [groups] Created '{new[-1]}' → group {grp['Protein_LIG']}")
        except Exception as e:
            print(f"  [groups] make_ndx failed: {e}")

    print(f"\n  {'Role':<16}  {'NDX name':<22}  Group#")
    print(f"  {'-'*16}  {'-'*22}  {'------'}")
    for key, cands in candidates.items():
        num = grp.get(key, "NOT FOUND")
        matched = next((c for c in cands if c in raw), "—")
        flag = "  ← WARNING" if num == "NOT FOUND" else ""
        print(f"  {key:<16}  {matched:<22}  {num}{flag}")
    print()
    return grp


def _configure(cfg, config_dir):
    global SIM_DIR, OUT_DIR, TABLES_DIR, TPR, XTC, EDR, GRO, NDX
    global XTC_NOJUMP, XTC_CENTER
    global LIG_RESNAME, GATEKEEPER_RES, SIM_LABEL, GROMACS
    global GRP, COLORS
    global PLOT_DPI, PLOT_W, PLOT_H, LIG_MARKER_SIZE, TITLE_FONTSIZE, LABEL_FONTSIZE
    global CONTACT_CUTOFF, HBOND_DIST, HBOND_ANGLE, MIN_CONTACT_PCT
    global TIME_UNIT, TIME_DIV, TIME_LABEL

    paths   = cfg.get("paths",    {})
    sys_cfg = cfg.get("system",   {})
    fig_cfg = cfg.get("figure",   {})
    ana_cfg = cfg.get("analysis", {})

    SIM_DIR = _resolve(config_dir, paths.get("sim_dir", "."))
    OUT_DIR = _resolve(config_dir, paths.get("out_dir", "results"))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR = OUT_DIR / "tables"
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    _SKIP = {"md_nojump.xtc", "md_center.xtc"}

    def _find(key, ext, fallback):
        name = paths.get(key, "")
        if name and (SIM_DIR / name).exists():
            return SIM_DIR / name
        hits = [h for h in sorted(SIM_DIR.glob(f"*.{ext}")) if h.name not in _SKIP]
        if not hits:
            return SIM_DIR / (name or fallback)
        md = [h for h in hits if "md" in h.stem.lower()]
        return (md or hits)[0]

    TPR = _find("tpr", "tpr", "md.tpr")
    XTC = _find("xtc", "xtc", "md.xtc")
    EDR = _find("edr", "edr", "md.edr")
    GRO = _find("gro", "gro", "md.gro")
    NDX = _find("ndx", "ndx", "index.ndx")
    print(f"  [files] tpr={TPR.name}  xtc={XTC.name}  edr={EDR.name}  gro={GRO.name}  ndx={NDX.name}")

    XTC_NOJUMP = OUT_DIR / "md_nojump.xtc"
    XTC_CENTER = OUT_DIR / "md_center.xtc"

    SIM_LABEL = sys_cfg.get("sim_label", "simulation")
    GROMACS   = sys_cfg.get("gromacs_cmd", "gmx")

    lig_raw = sys_cfg.get("lig_resname", "").strip()
    if lig_raw:
        LIG_RESNAME = lig_raw
        print(f"  [ligand] Config: '{LIG_RESNAME}'")
    else:
        detected, status = _detect_ligand(GRO)
        if status == "found":
            LIG_RESNAME = detected
            print(f"  [ligand] Auto-detected: '{LIG_RESNAME}'")
        elif status == "multiple":
            LIG_RESNAME = detected[0]
            print(f"  [ligand] Multiple found: {', '.join(detected)} — using '{LIG_RESNAME}'")
            print("           Set lig_resname in [system] to override.")
        else:
            LIG_RESNAME = "LIG"
            print("  [ligand] Could not detect — defaulting to 'LIG'")

    try:
        GATEKEEPER_RES = int(sys_cfg.get("gatekeeper_res", "0").strip() or "0")
    except ValueError:
        GATEKEEPER_RES = 0
    if GATEKEEPER_RES == 0:
        print("  [gatekeeper] Not set — gatekeeper analysis skipped")
    else:
        print(f"  [gatekeeper] Tracking residue {GATEKEEPER_RES}")

    print("\n── Auto-detecting GROMACS index groups ─────────────────────────")
    auto = _auto_groups(NDX, LIG_RESNAME, GRO, GROMACS)
    GRP.clear(); GRP.update(auto)
    overrides = {k: str(v) for k, v in cfg.get("groups", {}).items()}
    if overrides:
        GRP.update(overrides)
        print(f"  [groups] Overrides: {list(overrides.keys())}")

    PLOT_DPI        = int(fig_cfg.get("dpi",                 150))
    PLOT_W          = float(fig_cfg.get("width",              10))
    PLOT_H          = float(fig_cfg.get("height",              4))
    LIG_MARKER_SIZE = int(fig_cfg.get("ligand_marker_size",   60))
    TITLE_FONTSIZE  = int(fig_cfg.get("title_fontsize",       13))
    LABEL_FONTSIZE  = int(fig_cfg.get("label_fontsize",       11))
    font_scale      = float(fig_cfg.get("font_scale",        1.1))

    COLORS.clear()
    COLORS.update({
        "backbone":   fig_cfg.get("backbone_color",   "#2c7bb6"),
        "ligand":     fig_cfg.get("ligand_color",     "#d7191c"),
        "hbond":      fig_cfg.get("hbond_color",      "#1a9641"),
        "gatekeeper": fig_cfg.get("gatekeeper_color", "#fd8d3c"),
        "rmsf":       fig_cfg.get("rmsf_color",       "#756bb1"),
        "energy":     fig_cfg.get("energy_color",     "#636363"),
    })

    sns.set_theme(style="whitegrid", palette="muted", font_scale=font_scale)

    CONTACT_CUTOFF  = float(ana_cfg.get("contact_cutoff_A",     4.5))
    HBOND_DIST      = float(ana_cfg.get("hbond_distance_A",     3.5))
    HBOND_ANGLE     = float(ana_cfg.get("hbond_angle_deg",      120))
    MIN_CONTACT_PCT = float(ana_cfg.get("min_plot_contact_pct", 10))

    tu = ana_cfg.get("time_unit", "ps").strip().lower()
    if tu == "ns":
        TIME_UNIT, TIME_DIV, TIME_LABEL = "ns", 1000.0, "Time (ns)"
    else:
        TIME_UNIT, TIME_DIV, TIME_LABEL = "ps", 1.0,    "Time (ps)"


def _validate_inputs():
    ok = True
    r = subprocess.run([GROMACS, "--version"], capture_output=True, text=True)
    if r.returncode != 0:
        print(f"\n[ERROR] GROMACS '{GROMACS}' not found.")
        print("  module load gromacs  or  set gromacs_cmd in [system]")
        ok = False
    missing = {k: v for k, v in
               {"TPR": TPR, "XTC": XTC, "EDR": EDR, "GRO": GRO, "NDX": NDX}.items()
               if not v.exists()}
    if missing:
        print(f"\n[ERROR] Missing files in {SIM_DIR}:")
        for k, v in missing.items():
            print(f"  {k}: {v.name}")
        found = [f.name for ext in ("tpr","xtc","gro","edr","ndx")
                 for f in sorted(SIM_DIR.glob(f"*.{ext}"))]
        if found:
            print(f"  Found: {found}")
            print("  Set file names in [paths] of md_config.ini if they differ.")
        ok = False
    return ok


def run_gmx(args, stdin="", label=""):
    cmd = [GROMACS] + args
    print(f"  [gmx] {' '.join(cmd[:4])} ...")
    r = subprocess.run(cmd, input=stdin, text=True, capture_output=True, cwd=str(SIM_DIR))
    if r.returncode != 0:
        print(f"  [WARN] gmx {label or args[0]} exited {r.returncode}")
        for line in r.stderr.splitlines()[-6:]:
            print(f"    {line}")
    return r.stdout, r.stderr, r.returncode


def parse_xvg(path):
    headers, rows = [], []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith("@ s") and "legend" in line:
            headers.append(line.split('"')[1])
        elif not line.startswith(("#", "@")) and line:
            rows.append([float(x) for x in line.split()])
    arr = np.array(rows)
    cols = ["time_ps"] + (headers if len(headers) == arr.shape[1]-1
                          else [f"col{i}" for i in range(1, arr.shape[1])])
    return pd.DataFrame(arr, columns=cols)


def save_fig(fig, name):
    p = OUT_DIR / name
    fig.savefig(p, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {p}")


def preprocess_trajectory():
    print("\n── Step 1: PBC correction")
    if XTC_CENTER.exists():
        print(f"  Already exists: {XTC_CENTER.name}  (skipping)")
        return XTC_CENTER
    run_gmx(
        ["trjconv", "-s", str(TPR), "-f", str(XTC),
         "-o", str(XTC_NOJUMP), "-pbc", "nojump", "-n", str(NDX)],
        stdin=f"{GRP['System']}\n", label="trjconv-nojump",
    )
    run_gmx(
        ["trjconv", "-s", str(TPR), "-f", str(XTC_NOJUMP),
         "-o", str(XTC_CENTER), "-center", "-pbc", "mol", "-n", str(NDX)],
        stdin=f"{GRP['Protein_LIG']}\n{GRP['System']}\n", label="trjconv-center",
    )
    if XTC_CENTER.exists():
        print(f"  Written: {XTC_CENTER.name}")
    else:
        print("  [ERROR] Centred trajectory not created — check GROMACS output above")
    return XTC_CENTER


def energy_analysis():
    print("\n── Step 2: Energy analysis")
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
            run_gmx(["energy", "-f", str(EDR), "-o", str(xvg)],
                    stdin=f"{idx}\n0\n", label=f"energy-{name}")
        if xvg.exists():
            df = parse_xvg(xvg)
            df.columns = ["time_ps", name]
            data[name] = df
            print(f"  {name:15s}: {df[name].mean():>12.3f} ± {df[name].std():.3f}")

    for name, (_, ylabel) in terms.items():
        if name not in data:
            continue
        df  = data[name]
        mu  = df[name].mean()
        std = df[name].std()
        fig, ax = plt.subplots(figsize=(PLOT_W, PLOT_H))
        ax.plot(t_axis(df["time_ps"]), df[name], color=COLORS["energy"], lw=1.0)
        ax.fill_between(t_axis(df["time_ps"]), mu-std, mu+std,
                        color="red", alpha=0.08, label=f"±1σ  ({std:.3f})")
        ax.axhline(mu, color="red", ls="--", lw=1.0, label=f"mean = {mu:.3f}")
        ax.set_xlabel(TIME_LABEL, fontsize=LABEL_FONTSIZE)
        ax.set_ylabel(ylabel, fontsize=LABEL_FONTSIZE)
        ax.set_title(f"{SIM_LABEL} — {name.replace('_',' ').title()}",
                     fontsize=TITLE_FONTSIZE)
        ax.legend(fontsize=9)
        plt.tight_layout()
        save_fig(fig, f"energy_{name}.png")
    return data


def rmsd_analysis(traj):
    print("\n── Step 3: RMSD analysis")
    if not HAS_MDA:
        print("  [SKIP] MDAnalysis not available"); return pd.DataFrame()

    u   = mda.Universe(str(GRO), str(traj))
    ref = mda.Universe(str(GRO), str(traj))
    R = rms.RMSD(u, ref,
                 select="backbone and name CA",
                 groupselections=[f"resname {LIG_RESNAME} and not name H*"],
                 ref_frame=0)
    R.run()

    arr = R.results.rmsd
    df = pd.DataFrame({
        "time_ps":       arr[:, 1],
        "backbone_RMSD": arr[:, 2],
        "ligand_RMSD":   arr[:, 3],
    })
    for col in ("backbone_RMSD", "ligand_RMSD"):
        print(f"  {col}: {df[col].mean():.3f} ± {df[col].std():.3f} Å")
    df.to_csv(TABLES_DIR / "rmsd.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(PLOT_W, PLOT_H*1.75), sharex=True)
    axes[0].plot(t_axis(df["time_ps"]), df["backbone_RMSD"],
                 color=COLORS["backbone"], lw=1.2, label="Backbone Cα")
    axes[0].axhline(df["backbone_RMSD"].mean(), color="k", ls="--", lw=0.8)
    axes[0].set_ylabel("RMSD (Å)", fontsize=LABEL_FONTSIZE)
    axes[0].set_title(f"{SIM_LABEL} — Protein Backbone RMSD", fontsize=TITLE_FONTSIZE)
    axes[0].legend()

    axes[1].plot(t_axis(df["time_ps"]), df["ligand_RMSD"],
                 color=COLORS["ligand"], lw=1.2, label=f"Ligand ({LIG_RESNAME})")
    axes[1].axhline(df["ligand_RMSD"].mean(), color="k", ls="--", lw=0.8)
    axes[1].set_xlabel(TIME_LABEL, fontsize=LABEL_FONTSIZE)
    axes[1].set_ylabel("RMSD (Å)", fontsize=LABEL_FONTSIZE)
    axes[1].set_title(f"{SIM_LABEL} — Ligand RMSD", fontsize=TITLE_FONTSIZE)
    axes[1].legend()
    plt.tight_layout()
    save_fig(fig, "rmsd.png")
    return df


def rmsf_analysis(traj):
    print("\n── Step 4: RMSF analysis")
    if not HAS_MDA:
        print("  [SKIP] MDAnalysis not available"); return pd.DataFrame()

    u = mda.Universe(str(GRO), str(traj))
    align.AlignTraj(u, u, select="backbone and name CA", in_memory=True).run()
    ca = u.select_atoms("backbone and name CA")
    rmsf_calc = rms.RMSF(ca)
    rmsf_calc.run()

    df = pd.DataFrame({
        "resid":   ca.resids,
        "resname": ca.resnames,
        "rmsf_A":  rmsf_calc.results.rmsf,
    })
    df.to_csv(TABLES_DIR / "rmsf.csv", index=False)
    print("  Top-5 most flexible residues:")
    print(textwrap.indent(df.nlargest(5, "rmsf_A")[["resid","resname","rmsf_A"]].to_string(index=False), "    "))

    fig, ax = plt.subplots(figsize=(PLOT_W*1.4, PLOT_H))
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
    ax.set_xlabel("Residue ID", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("RMSF (Å)", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"{SIM_LABEL} — Per-residue Cα RMSF", fontsize=TITLE_FONTSIZE)
    plt.tight_layout()
    save_fig(fig, "rmsf.png")
    return df


def gyration_analysis():
    print("\n── Step 5: Radius of gyration")
    xvg = TABLES_DIR / "gyrate.xvg"
    if not xvg.exists():
        run_gmx(["gyrate", "-s", str(TPR), "-f", str(XTC_CENTER),
                 "-o", str(xvg), "-n", str(NDX)],
                stdin=f"{GRP['Protein']}\n", label="gyrate")
    if not xvg.exists():
        print("  [SKIP] gyrate.xvg not produced"); return pd.DataFrame()

    df = parse_xvg(xvg)
    df.columns = ["time_ps", "Rg", "Rg_X", "Rg_Y", "Rg_Z"][:df.shape[1]]
    print(f"  Mean Rg: {df['Rg'].mean():.3f} ± {df['Rg'].std():.3f} nm")

    fig, ax = plt.subplots(figsize=(PLOT_W, PLOT_H))
    ax.plot(t_axis(df["time_ps"]), df["Rg"]*10, color=COLORS["energy"], lw=1.2)
    ax.axhline(df["Rg"].mean()*10, color="k", ls="--", lw=0.8,
               label=f"mean = {df['Rg'].mean()*10:.2f} Å")
    ax.set_xlabel(TIME_LABEL, fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Radius of Gyration (Å)", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"{SIM_LABEL} — Protein Radius of Gyration", fontsize=TITLE_FONTSIZE)
    ax.legend()
    plt.tight_layout()
    save_fig(fig, "gyration.png")
    return df


def hbond_analysis(traj):
    print("\n── Step 6: H-bond analysis")
    if not HAS_MDA:
        print("  [SKIP] MDAnalysis not available"); return pd.DataFrame()

    u = mda.Universe(str(GRO), str(traj))
    prot = "protein"
    lig  = f"resname {LIG_RESNAME}"
    hba  = HBA(
        universe=u,
        donors_sel=f"({prot} and (name N O)) or ({lig} and (name N* O*))",
        hydrogens_sel="name H*",
        acceptors_sel=f"({prot} and (name O N)) or ({lig} and (name N* O*))",
        between=[lig, prot],
        d_a_cutoff=HBOND_DIST,
        d_h_a_angle_cutoff=HBOND_ANGLE,
        update_selections=False,
    )
    hba.run()

    counts = hba.count_by_time()
    times  = [ts.time for ts in u.trajectory]
    time_arr = np.array(times) if len(times) == len(counts) else np.arange(len(counts))
    df_count = pd.DataFrame({"time_ps": time_arr, "n_hbonds": counts})
    df_count.to_csv(TABLES_DIR / "hbonds_count.csv", index=False)
    print(f"  Mean H-bonds/frame: {df_count['n_hbonds'].mean():.2f} ± {df_count['n_hbonds'].std():.2f}")

    pair_occ = pd.DataFrame()
    if len(hba.results.hbonds) > 0:
        hb_df = pd.DataFrame(hba.results.hbonds,
                             columns=["frame","donor_ix","H_ix","acceptor_ix",
                                      "distance_A","angle_deg"])
        atoms = u.atoms
        for side, ix_col in [("donor","donor_ix"),("acc","acceptor_ix")]:
            for attr, name in [("resid","resid"),("resname","resname"),("name","name")]:
                hb_df[f"{side}_{name}"] = (hb_df[ix_col].astype(int)
                                           .map(lambda i: getattr(atoms[i], attr)))
        n_frames = len(u.trajectory)
        pair_occ = (
            hb_df.groupby(["donor_resid","donor_resname","donor_name",
                            "acc_resid","acc_resname","acc_name"])
            .size().reset_index(name="n_frames")
        )
        pair_occ["occupancy_%"] = (pair_occ["n_frames"] / n_frames * 100).round(1)
        pair_occ = pair_occ.sort_values("occupancy_%", ascending=False)
        pair_occ.to_csv(TABLES_DIR / "hbonds_pairs.csv", index=False)
        print("  Top H-bond pairs (occupancy %):")
        print(textwrap.indent(pair_occ.head(10).to_string(index=False), "    "))
        if GATEKEEPER_RES > 0:
            gk_hb = pair_occ[(pair_occ["donor_resid"] == GATEKEEPER_RES) |
                             (pair_occ["acc_resid"]   == GATEKEEPER_RES)]
            if gk_hb.empty:
                print(f"  [NOTE] No H-bond with gatekeeper residue {GATEKEEPER_RES}")
            else:
                print(f"\n  Gatekeeper {GATEKEEPER_RES} H-bond occupancy:")
                print(textwrap.indent(gk_hb.to_string(index=False), "    "))
    else:
        print("  [NOTE] No protein↔ligand H-bonds detected")

    fig, ax = plt.subplots(figsize=(PLOT_W, PLOT_H))
    ax.fill_between(t_axis(df_count["time_ps"]), df_count["n_hbonds"],
                    alpha=0.4, color=COLORS["hbond"])
    ax.plot(t_axis(df_count["time_ps"]), df_count["n_hbonds"],
            color=COLORS["hbond"], lw=1.0)
    ax.axhline(df_count["n_hbonds"].mean(), color="k", ls="--", lw=0.8,
               label=f"mean = {df_count['n_hbonds'].mean():.2f}")
    ax.set_xlabel(TIME_LABEL, fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Number of H-bonds", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"{SIM_LABEL} — Protein–Ligand H-bonds", fontsize=TITLE_FONTSIZE)
    ax.legend()
    plt.tight_layout()
    save_fig(fig, "hbonds.png")
    return pair_occ


def gatekeeper_distance(traj):
    print("\n── Step 7: Gatekeeper distance")
    if GATEKEEPER_RES == 0:
        print("  [SKIP] gatekeeper_res not set in md_config.ini [system]")
        return pd.DataFrame()
    if not HAS_MDA:
        print("  [SKIP] MDAnalysis not available"); return pd.DataFrame()

    from MDAnalysis.lib.distances import distance_array as mda_dist
    u        = mda.Universe(str(GRO), str(traj))
    lig_atm  = u.select_atoms(f"resname {LIG_RESNAME} and not name H*")
    gk_atm   = u.select_atoms(f"resid {GATEKEEPER_RES} and not name H*")
    if gk_atm.n_atoms == 0:
        print(f"  [WARN] Residue {GATEKEEPER_RES} not found")
        return pd.DataFrame()

    gk_resname = gk_atm.resnames[0]
    print(f"  Gatekeeper: {gk_resname}{GATEKEEPER_RES} ({gk_atm.n_atoms} heavy atoms)")

    times, dists = [], []
    for ts in u.trajectory:
        d = mda_dist(lig_atm.positions, gk_atm.positions, box=ts.dimensions)
        times.append(ts.time); dists.append(d.min())

    df = pd.DataFrame({"time_ps": times, "min_dist_A": dists})
    df.to_csv(TABLES_DIR / "gatekeeper_distance.csv", index=False)
    close_frac = (df["min_dist_A"] < HBOND_DIST).mean() * 100
    print(f"  Mean min-distance: {df['min_dist_A'].mean():.2f} ± {df['min_dist_A'].std():.2f} Å")
    print(f"  Frames < {HBOND_DIST} Å: {close_frac:.1f}%")

    fig, ax = plt.subplots(figsize=(PLOT_W, PLOT_H))
    ax.fill_between(t_axis(df["time_ps"]), df["min_dist_A"],
                    alpha=0.3, color=COLORS["gatekeeper"])
    ax.plot(t_axis(df["time_ps"]), df["min_dist_A"],
            color=COLORS["gatekeeper"], lw=1.2)
    ax.axhline(HBOND_DIST, color="k", ls="--", lw=0.9,
               label=f"H-bond cutoff ({HBOND_DIST} Å)")
    ax.axhline(df["min_dist_A"].mean(), color="red", ls=":", lw=0.9,
               label=f"mean = {df['min_dist_A'].mean():.2f} Å")
    ax.set_xlabel(TIME_LABEL, fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Min distance (Å)", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"{SIM_LABEL} — Ligand–{gk_resname}{GATEKEEPER_RES} Distance",
                 fontsize=TITLE_FONTSIZE)
    ax.legend()
    plt.tight_layout()
    save_fig(fig, "gatekeeper_distance.png")
    return df


def pocket_contacts(traj):
    print(f"\n── Step 8: Pocket contact map (cutoff {CONTACT_CUTOFF} Å)")
    if not HAS_MDA:
        print("  [SKIP] MDAnalysis not available"); return pd.DataFrame()

    from MDAnalysis.lib.distances import distance_array as mda_dist
    u    = mda.Universe(str(GRO), str(traj))
    lig  = u.select_atoms(f"resname {LIG_RESNAME} and not name H*")
    prot = u.select_atoms("protein and not name H*")
    resids      = np.unique(prot.resids)
    resnames    = {r: prot.select_atoms(f"resid {r}").resnames[0] for r in resids}
    contact_count = {r: 0 for r in resids}
    prot_resid_arr = prot.resids          # precompute once; shape (n_atoms,)
    n_frames = 0
    for ts in u.trajectory:
        d = mda_dist(lig.positions, prot.positions, box=ts.dimensions)
        in_contact = d.min(axis=0) < CONTACT_CUTOFF
        for r in np.unique(prot_resid_arr[in_contact]):  # only touched residues
            contact_count[r] += 1
        n_frames += 1

    df = pd.DataFrame({
        "resid":        list(contact_count.keys()),
        "resname":      [resnames[r] for r in contact_count],
        "contact_frac": [v/n_frames for v in contact_count.values()],
    })
    df = df[df["contact_frac"] > 0].sort_values("contact_frac", ascending=False)
    df["contact_%"] = (df["contact_frac"] * 100).round(1)
    df.to_csv(TABLES_DIR / "pocket_contacts.csv", index=False)
    print(f"  Residues in contact: {len(df)}")
    print(f"  Residues > 80%: {(df['contact_%'] > 80).sum()}")
    print("\n  Top-15:")
    print(textwrap.indent(df.head(15)[["resid","resname","contact_%"]].to_string(index=False), "    "))

    df_plot = df[df["contact_%"] >= MIN_CONTACT_PCT].copy()
    df_plot["label"] = df_plot["resname"] + df_plot["resid"].astype(str)
    bar_colors = [COLORS["gatekeeper"] if r == GATEKEEPER_RES else COLORS["backbone"]
                  for r in df_plot["resid"]]

    fig, ax = plt.subplots(figsize=(max(PLOT_W, len(df_plot)*0.45), PLOT_H*1.25))
    ax.bar(df_plot["label"], df_plot["contact_%"],
           color=bar_colors, edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Residue", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Contact occupancy (%)", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"{SIM_LABEL} — Binding-pocket Contacts (≥{MIN_CONTACT_PCT}%)",
                 fontsize=TITLE_FONTSIZE)
    ax.set_ylim(0, 105)
    plt.xticks(rotation=45, ha="right", fontsize=8)
    from matplotlib.patches import Patch
    legend_els = [Patch(facecolor=COLORS["backbone"], label="Pocket residue")]
    if GATEKEEPER_RES > 0:
        legend_els.append(Patch(facecolor=COLORS["gatekeeper"],
                                label=f"Gatekeeper ({GATEKEEPER_RES})"))
    ax.legend(handles=legend_els, loc="upper right")
    plt.tight_layout()
    save_fig(fig, "pocket_contacts.png")
    return df


def summary_figure(rmsd_df, rg_df, rmsf_df):
    print("\n── Step 9: Summary figure")
    fig = plt.figure(figsize=(PLOT_W*1.4, PLOT_H*2.5))
    gs  = gridspec.GridSpec(2, 2, hspace=0.4, wspace=0.35)
    ax1, ax2 = fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1])
    ax3, ax4 = fig.add_subplot(gs[1,0]), fig.add_subplot(gs[1,1])

    def _no_data(ax):
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)

    if not rmsd_df.empty and "backbone_RMSD" in rmsd_df.columns:
        ax1.plot(t_axis(rmsd_df["time_ps"]), rmsd_df["backbone_RMSD"],
                 color=COLORS["backbone"], lw=1.2)
        ax1.axhline(rmsd_df["backbone_RMSD"].mean(), color="k", ls="--", lw=0.8)
        ax1.set_xlabel(TIME_LABEL, fontsize=LABEL_FONTSIZE)
        ax1.set_ylabel("RMSD (Å)", fontsize=LABEL_FONTSIZE)
        ax1.set_title("Backbone Cα RMSD", fontsize=TITLE_FONTSIZE)
    else:
        _no_data(ax1)

    if not rmsd_df.empty and "ligand_RMSD" in rmsd_df.columns:
        ax2.plot(t_axis(rmsd_df["time_ps"]), rmsd_df["ligand_RMSD"],
                 color=COLORS["ligand"], lw=1.2)
        ax2.axhline(rmsd_df["ligand_RMSD"].mean(), color="k", ls="--", lw=0.8)
        ax2.set_xlabel(TIME_LABEL, fontsize=LABEL_FONTSIZE)
        ax2.set_ylabel("RMSD (Å)", fontsize=LABEL_FONTSIZE)
        ax2.set_title(f"Ligand ({LIG_RESNAME}) RMSD", fontsize=TITLE_FONTSIZE)
    else:
        _no_data(ax2)

    if not rg_df.empty:
        ax3.plot(t_axis(rg_df["time_ps"]), rg_df["Rg"]*10,
                 color=COLORS["energy"], lw=1.2)
        ax3.axhline(rg_df["Rg"].mean()*10, color="k", ls="--", lw=0.8,
                    label=f"mean = {rg_df['Rg'].mean()*10:.2f} Å")
        ax3.set_xlabel(TIME_LABEL, fontsize=LABEL_FONTSIZE)
        ax3.set_ylabel("Radius of Gyration (Å)", fontsize=LABEL_FONTSIZE)
        ax3.set_title("Protein Radius of Gyration", fontsize=TITLE_FONTSIZE)
        ax3.legend(fontsize=8)
    else:
        _no_data(ax3)

    if not rmsf_df.empty:
        ax4.fill_between(rmsf_df["resid"], rmsf_df["rmsf_A"],
                         alpha=0.35, color=COLORS["rmsf"])
        ax4.plot(rmsf_df["resid"], rmsf_df["rmsf_A"], color=COLORS["rmsf"], lw=0.9)
        ax4.set_xlabel("Residue ID", fontsize=LABEL_FONTSIZE)
        ax4.set_ylabel("RMSF (Å)", fontsize=LABEL_FONTSIZE)
        ax4.set_title("Per-residue Cα RMSF", fontsize=TITLE_FONTSIZE)
    else:
        _no_data(ax4)

    fig.suptitle(f"{SIM_LABEL} — MD Analysis Summary",
                 fontsize=TITLE_FONTSIZE+2, fontweight="bold")
    save_fig(fig, "summary_panel.png")


def write_report(energy_data, rmsd_df, rmsf_df, gk_df, contacts_df):
    print("\n── Step 10: Summary report")
    lines = [
        "=" * 70,
        "  MD ANALYSIS REPORT",
        f"  Simulation: {SIM_LABEL}",
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
            lines.append(f"  {name:<22}: {vals.mean():>12.3f} ± {vals.std():.3f} {unit}")

    lines += ["", "─"*70, "STRUCTURAL STABILITY"]
    if not rmsd_df.empty:
        bk = rmsd_df["backbone_RMSD"]
        lg = rmsd_df["ligand_RMSD"]
        lines += [
            f"  Backbone Cα RMSD  : {bk.mean():.3f} ± {bk.std():.3f} Å  (max {bk.max():.3f} Å)",
            f"  Ligand RMSD       : {lg.mean():.3f} ± {lg.std():.3f} Å  (max {lg.max():.3f} Å)",
        ]
        if lg.mean() > 2.0 or lg.max() > 3.0:
            lines.append("  [FLAG] Ligand RMSD elevated — possible pose instability")
        if bk.mean() > 3.5:
            lines.append("  [FLAG] Mean backbone RMSD > 3.5 Å — consider longer equilibration or staged restraint release")
        elif bk.max() > 5.0:
            lines.append("  [FLAG] Backbone RMSD spike > 5 Å — transient instability, check trajectory")

    if not rmsf_df.empty:
        lines += ["", "  Top-5 flexible residues (Cα RMSF):"]
        for _, row in rmsf_df.nlargest(5, "rmsf_A").iterrows():
            lines.append(f"    {row['resname']}{int(row['resid']):<6}: {row['rmsf_A']:.3f} Å")

    if GATEKEEPER_RES > 0:
        lines += ["", "─"*70, "GATEKEEPER CONTACT"]
    if GATEKEEPER_RES > 0 and not gk_df.empty:
        mu    = gk_df["min_dist_A"].mean()
        mn    = gk_df["min_dist_A"].min()
        close = (gk_df["min_dist_A"] < HBOND_DIST).mean() * 100
        lines += [
            f"  Gatekeeper residue: {GATEKEEPER_RES}",
            f"  Mean min-distance : {mu:.2f} Å  (min = {mn:.2f} Å)",
            f"  Frames < {HBOND_DIST} Å   : {close:.1f}%",
        ]
        lines.append("  [OK]  Contact > 50%" if close > 50
                     else "  [NOTE] Contact < 50% — check H-bond stability")

    lines += ["", "─"*70, "BINDING-POCKET CONTACTS"]
    if not contacts_df.empty:
        high = contacts_df[contacts_df["contact_%"] >= 80]
        lines += [
            f"  Total residues in contact: {len(contacts_df)}",
            f"  Residues ≥ 80% occupancy : {len(high)}",
            "", "  High-occupancy contacts (≥ 80%):",
        ]
        for _, r in high.iterrows():
            lines.append(f"    {r['resname']}{int(r['resid']):<6}: {r['contact_%']:.1f}%")

    lines += [
        "", "─"*70, "OUTPUT FILES",
        f"  Figures → {OUT_DIR}/",
        f"  Tables  → {TABLES_DIR}/",
        "=" * 70,
    ]

    report_path = TABLES_DIR / "analysis_report.txt"
    report_path.write_text("\n".join(lines))
    print(f"  → {report_path}")
    print("\n" + "\n".join(lines))


def main():
    p = argparse.ArgumentParser(
        prog="md_analysis.py",
        description="GROMACS MD analysis pipeline — all settings via md_config.ini",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            examples:
              python md_analysis.py /path/to/1IEP_MD
              python md_analysis.py /path/to/T315I_MD
              python md_analysis.py --config /path/to/1IEP_MD/md_config.ini
        """),
    )
    p.add_argument("sim_path", nargs="?", metavar="SIM_FOLDER_OR_CONFIG",
                   help="Path to simulation folder (must contain md_config.ini) "
                        "or path to an md_config.ini file directly.")
    p.add_argument("--config", metavar="FILE",
                   help="Explicit path to md_config.ini (alternative to positional).")
    args = p.parse_args()

    # Resolve config path
    if args.config:
        config_path = Path(args.config).resolve()
    elif args.sim_path:
        sp = Path(args.sim_path).resolve()
        config_path = sp if sp.is_file() else sp / "md_config.ini"
    else:
        # Fall back to script directory
        config_path = Path(__file__).resolve().parent / "md_config.ini"

    if not config_path.exists():
        print(f"[ERROR] Config not found: {config_path}")
        print("  Place md_config.ini in the simulation folder and pass that folder as argument.")
        sys.exit(1)

    print(f"[config] {config_path}")
    cfg = _load_config(config_path)
    _configure(cfg, config_path.parent)

    print("=" * 70)
    print(f"  MD Analysis Pipeline")
    print(f"  Label      : {SIM_LABEL}")
    print(f"  Sim dir    : {SIM_DIR}")
    print(f"  Out dir    : {OUT_DIR}")
    print(f"  Ligand     : {LIG_RESNAME}")
    print(f"  Gatekeeper : {'residue ' + str(GATEKEEPER_RES) if GATEKEEPER_RES else 'skipped'}")
    print(f"  Time unit  : {TIME_UNIT}")
    print("=" * 70)

    if not _validate_inputs():
        sys.exit(1)
    if not HAS_MDA:
        print("\n[WARNING] MDAnalysis not installed — structural analyses will be skipped.\n")

    traj = None
    try:
        traj = preprocess_trajectory()
        if not traj.exists():
            raise RuntimeError("PBC-corrected trajectory not produced")
    except Exception as e:
        print(f"\n[ERROR] Step 1 (PBC correction) failed: {e}\n")

    energy_d = {}
    try:   energy_d  = energy_analysis()
    except Exception as e: print(f"  [WARN] Energy analysis failed: {e}\n")

    rmsd_df = pd.DataFrame()
    if traj and traj.exists():
        try:   rmsd_df = rmsd_analysis(traj)
        except Exception as e: print(f"  [WARN] RMSD analysis failed: {e}\n")

    rmsf_df = pd.DataFrame()
    if traj and traj.exists():
        try:   rmsf_df = rmsf_analysis(traj)
        except Exception as e: print(f"  [WARN] RMSF analysis failed: {e}\n")

    rg_df = pd.DataFrame()
    try:   rg_df = gyration_analysis()
    except Exception as e: print(f"  [WARN] Gyration analysis failed: {e}\n")

    if traj and traj.exists():
        try:   hbond_analysis(traj)
        except Exception as e: print(f"  [WARN] H-bond analysis failed: {e}\n")

    gk_df = pd.DataFrame()
    if traj and traj.exists():
        try:   gk_df = gatekeeper_distance(traj)
        except Exception as e: print(f"  [WARN] Gatekeeper distance failed: {e}\n")

    contacts_df = pd.DataFrame()
    if traj and traj.exists():
        try:   contacts_df = pocket_contacts(traj)
        except Exception as e: print(f"  [WARN] Contact map failed: {e}\n")

    try:   summary_figure(rmsd_df, rg_df, rmsf_df)
    except Exception as e: print(f"  [WARN] Summary figure failed: {e}\n")

    try:   write_report(energy_d, rmsd_df, rmsf_df, gk_df, contacts_df)
    except Exception as e: print(f"  [WARN] Report failed: {e}\n")

    print(f"\n[DONE]  Figures → {OUT_DIR}/")
    print(f"        Tables  → {TABLES_DIR}/\n")


if __name__ == "__main__":
    main()
