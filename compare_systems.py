#!/usr/bin/env python3
# compare_systems.py — Compare the 4-system 2x2 resistance panel, 3 replicates each:
#   S1 1IEP  = WT   + imatinib   (baseline)
#   S2 T315I = T315I + imatinib   (resistance)
#   S3 3OXZ  = WT   + ponatinib  (control)
#   S4 3IK3  = T315I + ponatinib  (rescue)
# Usage:  python compare_systems.py

from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.signal import find_peaks
import warnings; warnings.filterwarnings("ignore")

ROOT = Path(__file__).parent

SYSTEMS = {
    "WT+imatinib (1IEP)": [
        ROOT / "1IEP_MD/Replica1/analysis/results/tables",
        ROOT / "1IEP_MD/Replica2/analysis/results/tables",
        ROOT / "1IEP_MD/Replica3/analysis/results/tables",
    ],
    "T315I+imatinib": [
        ROOT / "T315I_MD/Rep1/results/tables",
        ROOT / "T315I_MD/Rep2/results/tables",
        ROOT / "T315I_MD/Rep3/results/tables",
    ],
    "WT+ponatinib (3OXZ)": [
        ROOT / "3OXZ_MD/Repli1/results/tables",
        ROOT / "3OXZ_MD/Repli2/results/tables",
        ROOT / "3OXZ_MD/Repli3/results/tables",
    ],
    "T315I+ponatinib (3IK3)": [
        ROOT / "3IK3_MD/Rep1/results/tables",
        ROOT / "3IK3_MD/Rep2/results/tables",
        ROOT / "3IK3_MD/Rep3/tables",   # note: Rep3 has tables/ directly, not results/tables
    ],
}

# imatinib pair = warm (blue WT / red mutant), ponatinib pair = cool (green WT / purple mutant)
COLORS = {
    "WT+imatinib (1IEP)":     "#2c7bb6",
    "T315I+imatinib":         "#d7191c",
    "WT+ponatinib (3OXZ)":    "#1a9641",
    "T315I+ponatinib (3IK3)": "#7b3294",
}

# the two scientific contrasts (mutant − WT) for each drug
CONTRASTS = [
    ("imatinib", "WT+imatinib (1IEP)",  "T315I+imatinib"),
    ("ponatinib", "WT+ponatinib (3OXZ)", "T315I+ponatinib (3IK3)"),
]

# plot each drug's WT-vs-T315I pair on its own figure instead of all 4 systems
# crammed together — that was unreadable.
PAIRS = {
    "imatinib":  {k: v for k, v in SYSTEMS.items() if "imatinib" in k},
    "ponatinib": {k: v for k, v in SYSTEMS.items() if "ponatinib" in k},
}
ALPHA_TRACE = 0.25
ALPHA_MEAN  = 0.85

OUT = ROOT / "comparison"
OUT.mkdir(exist_ok=True)


# ── helpers ─────────────────────────────────────────────────────────────────

def load_csv(tables: Path, name: str) -> pd.DataFrame | None:
    p = tables / name
    return pd.read_csv(p) if p.exists() else None


def mean_sd_traces(dfs: list[pd.DataFrame], x_col: str, y_col: str, scale=1.0):
    """Interpolate each replicate onto a common x grid, return mean±SD."""
    x_min = max(df[x_col].min() for df in dfs)
    x_max = min(df[x_col].max() for df in dfs)
    xs = np.linspace(x_min, x_max, 500)
    interp = [np.interp(xs, df[x_col].values, df[y_col].values * scale) for df in dfs]
    arr = np.array(interp)
    return xs, arr.mean(0), arr.std(0), interp


def ps_to_ns(df: pd.DataFrame, col: str) -> pd.DataFrame:
    df = df.copy()
    if df[col].max() > 1000:
        df[col] = df[col] / 1000
    return df


# ── per-metric comparison plots ─────────────────────────────────────────────

def plot_rmsd(drug, systems):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=False)
    titles = ["Backbone Cα RMSD", "Ligand RMSD"]
    cols   = ["backbone_RMSD", "ligand_RMSD"]

    for ax, col, title in zip(axes, cols, titles):
        for sys, tdirs in systems.items():
            dfs = []
            for td in tdirs:
                df = load_csv(td, "rmsd.csv")
                if df is None: continue
                df = ps_to_ns(df, df.columns[0])
                dfs.append(df)
            if not dfs: continue
            xs, mean, sd, traces = mean_sd_traces(dfs, dfs[0].columns[0], col)
            c = COLORS[sys]
            for t in traces:
                ax.plot(xs, t, color=c, alpha=ALPHA_TRACE, lw=0.8)
            ax.plot(xs, mean, color=c, lw=2, alpha=ALPHA_MEAN, label=sys)
            ax.fill_between(xs, mean - sd, mean + sd, color=c, alpha=0.15)

        ax.set_xlabel("Time (ns)"); ax.set_ylabel("RMSD (Å)")
        ax.set_title(title, fontweight="bold"); ax.legend()
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"Structural Stability: WT vs T315I ({drug})", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / f"rmsd_comparison_{drug}.png", dpi=150)
    plt.close(fig)
    print(f"  Saved rmsd_comparison_{drug}.png")


def plot_gatekeeper(drug, systems):
    fig, ax = plt.subplots(figsize=(8, 4))
    for sys, tdirs in systems.items():
        dfs = []
        for td in tdirs:
            df = load_csv(td, "gatekeeper_distance.csv")
            if df is None: continue
            df = ps_to_ns(df, df.columns[0])
            dfs.append(df)
        if not dfs: continue
        xs, mean, sd, traces = mean_sd_traces(dfs, dfs[0].columns[0], dfs[0].columns[1])
        c = COLORS[sys]
        for t in traces:
            ax.plot(xs, t, color=c, alpha=ALPHA_TRACE, lw=0.8)
        ax.plot(xs, mean, color=c, lw=2, alpha=ALPHA_MEAN, label=sys)
        ax.fill_between(xs, mean - sd, mean + sd, color=c, alpha=0.15)

    ax.axhline(3.5, color="grey", ls="--", lw=1, label="3.5 Å threshold")
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Min distance to ligand (Å)")
    ax.set_title(f"Gatekeeper Residue 315 – Ligand Distance ({drug})", fontweight="bold")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / f"gatekeeper_comparison_{drug}.png", dpi=150)
    plt.close(fig)
    print(f"  Saved gatekeeper_comparison_{drug}.png")


def plot_hbonds(drug, systems):
    fig, ax = plt.subplots(figsize=(8, 4))
    for sys, tdirs in systems.items():
        dfs = []
        for td in tdirs:
            df = load_csv(td, "hbonds_count.csv")
            if df is None: continue
            df = ps_to_ns(df, df.columns[0])
            dfs.append(df)
        if not dfs: continue
        xs, mean, sd, traces = mean_sd_traces(dfs, dfs[0].columns[0], "n_hbonds")
        c = COLORS[sys]
        for t in traces:
            ax.plot(xs, t, color=c, alpha=ALPHA_TRACE, lw=0.8)
        ax.plot(xs, mean, color=c, lw=2, alpha=ALPHA_MEAN, label=sys)
        ax.fill_between(xs, mean - sd, mean + sd, color=c, alpha=0.15)

    ax.set_xlabel("Time (ns)"); ax.set_ylabel("# H-bonds (ligand–protein)")
    ax.set_title(f"Ligand–Protein Hydrogen Bonds ({drug})", fontweight="bold")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / f"hbonds_comparison_{drug}.png", dpi=150)
    plt.close(fig)
    print(f"  Saved hbonds_comparison_{drug}.png")


def plot_rmsf(drug, systems):
    fig, ax = plt.subplots(figsize=(11, 4))
    for sys, tdirs in systems.items():
        dfs = []
        for td in tdirs:
            df = load_csv(td, "rmsf.csv")
            if df is None: continue
            dfs.append(df)
        if not dfs: continue
        res_col = dfs[0].columns[0]
        rmsf_col = [c for c in dfs[0].columns if "rmsf" in c.lower()][0]
        # use union of residues present in all reps
        residues = sorted(set.intersection(*[set(df[res_col]) for df in dfs]))
        vals = []
        for df in dfs:
            sub = df[df[res_col].isin(residues)].set_index(res_col)[rmsf_col]
            vals.append(sub.reindex(residues).values)
        arr = np.array(vals)
        mean, sd = arr.mean(0), arr.std(0)
        c = COLORS[sys]
        ax.plot(residues, mean, color=c, lw=1.5, alpha=ALPHA_MEAN, label=sys)
        ax.fill_between(residues, mean - sd, mean + sd, color=c, alpha=0.15)

        # label the top 3 interior flexibility peaks with their residue number
        # (chain termini are excluded — their high RMSF is a free-end artifact, not a real peak)
        residues_arr = np.asarray(residues)
        interior = slice(5, -5)
        peaks, _ = find_peaks(mean[interior], prominence=0.3, distance=5)
        peaks = peaks + 5
        top_peaks = peaks[np.argsort(mean[peaks])[-3:]] if len(peaks) else []
        for p in top_peaks:
            ax.annotate(str(residues_arr[p]), xy=(residues_arr[p], mean[p]),
                        xytext=(0, 6), textcoords="offset points", ha="center",
                        fontsize=8, color=c, fontweight="bold")

    ax.set_xlabel("Residue number"); ax.set_ylabel("RMSF (Å)")
    ax.set_title(f"Per-residue Flexibility (Cα RMSF, {drug})", fontweight="bold")
    ax.legend(); ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(OUT / f"rmsf_comparison_{drug}.png", dpi=150)
    plt.close(fig)
    print(f"  Saved rmsf_comparison_{drug}.png")


def plot_pocket_contacts(drug, systems):
    """Bar chart of mean contact occupancy for residues seen in either system of the pair."""
    sys_contacts = {}
    for sys, tdirs in systems.items():
        all_df = []
        for td in tdirs:
            df = load_csv(td, "pocket_contacts.csv")
            if df is None: continue
            df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
            all_df.append(df)
        if not all_df: continue
        merged = {}
        for df in all_df:
            res_col2 = df.columns[0]
            occ_col2 = [c for c in df.columns if "occ" in c or "%" in c][0]
            for _, row in df.iterrows():
                k = f"{row['resname']}{row[res_col2]}" if 'resname' in df.columns else str(row[res_col2])
                merged.setdefault(k, []).append(float(row[occ_col2]))
        sys_contacts[sys] = {k: (np.mean(v), np.std(v)) for k, v in merged.items()}

    # residues in either system with mean ≥ 50% in at least one
    all_res = set()
    for d in sys_contacts.values():
        all_res |= {k for k, (m, _) in d.items() if m >= 50}
    all_res = sorted(all_res)

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(all_res))
    width = 0.38
    for i, (sys, d) in enumerate(sys_contacts.items()):
        means = [d.get(r, (0, 0))[0] for r in all_res]
        sds   = [d.get(r, (0, 0))[1] for r in all_res]
        offset = (i - 0.5) * width
        ax.bar(x + offset, means, width, yerr=sds, capsize=3,
               color=COLORS[sys], alpha=0.8, label=sys, error_kw={"elinewidth": 1})

    ax.set_xticks(x); ax.set_xticklabels(all_res, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mean contact occupancy (%)"); ax.set_ylim(0, 115)
    ax.set_title(f"Binding-pocket Contact Occupancy ({drug})", fontweight="bold")
    ax.axhline(80, color="grey", ls="--", lw=0.8, alpha=0.6, label="80% threshold")
    ax.legend(); ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / f"pocket_contacts_comparison_{drug}.png", dpi=150)
    plt.close(fig)
    print(f"  Saved pocket_contacts_comparison_{drug}.png")


# ── summary panel ────────────────────────────────────────────────────────────

def summary_bar(label, per_sys_values, ax, unit=""):
    """per_sys_values: {system_name: [per-replicate scalar, ...]}."""
    names  = list(per_sys_values.keys())
    means  = [np.mean(v) if v else np.nan for v in per_sys_values.values()]
    sds    = [np.std(v)  if v else 0.0    for v in per_sys_values.values()]
    colors = [COLORS[n] for n in names]
    bars = ax.bar(range(len(names)), means, yerr=sds, capsize=6, color=colors, alpha=0.85,
                  error_kw={"elinewidth": 1.5})
    ax.set_title(label, fontweight="bold")
    ax.set_ylabel(unit)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([n.split(" (")[0] for n in names], rotation=25, ha="right", fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    for bar, m, s in zip(bars, means, sds):
        if np.isnan(m): continue
        ax.text(bar.get_x() + bar.get_width()/2, m + s + (ax.get_ylim()[1] * 0.02),
                f"{m:.2f}", ha="center", va="bottom", fontsize=8)


def plot_summary(drug, systems):
    # per-replicate scalar mean for a metric, for every system in the pair
    def collect_all(csv_name, col):
        out = {}
        for sys, tdirs in systems.items():
            vals = []
            for td in tdirs:
                df = load_csv(td, csv_name)
                if df is None: continue
                vals.append(float(df[col].mean()))
            out[sys] = vals
        return out

    def collect_gk():
        out = {}
        for sys, tdirs in systems.items():
            vals = []
            for td in tdirs:
                df = load_csv(td, "gatekeeper_distance.csv")
                if df is None: continue
                vals.append(float(df.iloc[:, 1].mean()))
            out[sys] = vals
        return out

    metrics = [
        ("Backbone RMSD",    collect_all("rmsd.csv", "backbone_RMSD"), "Å"),
        ("Ligand RMSD",      collect_all("rmsd.csv", "ligand_RMSD"),   "Å"),
        ("Gatekeeper dist.", collect_gk(),                             "Å"),
        ("H-bonds (mean)",   collect_all("hbonds_count.csv", "n_hbonds"), "count"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(14, 4.5))
    for ax, (label, per_sys, unit) in zip(axes, metrics):
        summary_bar(label, per_sys, ax, unit)
    fig.suptitle(f"Key Metrics Summary — WT vs T315I ({drug}, mean ± SD across 3 replicates)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / f"summary_comparison_{drug}.png", dpi=150)
    plt.close(fig)
    print(f"  Saved summary_comparison_{drug}.png")


# ── text report ──────────────────────────────────────────────────────────────

def write_report():
    def stats(tdirs, csv_name, col):
        vals = []
        for td in tdirs:
            df = load_csv(td, csv_name)
            if df is None: continue
            if col not in df.columns: continue
            vals.extend(df[col].dropna().tolist())
        if not vals: return float("nan"), float("nan")
        return float(np.mean(vals)), float(np.std(vals))

    def gk_stats(tdirs):
        means, pcts = [], []
        for td in tdirs:
            df = load_csv(td, "gatekeeper_distance.csv")
            if df is None: continue
            col = df.columns[1]
            means.append(df[col].mean())
            pcts.append((df[col] < 3.5).mean() * 100)
        if not means: return float("nan"), float("nan"), float("nan")
        return float(np.mean(means)), float(np.std(means)), float(np.mean(pcts))

    # per-system scalar table
    rows = {}
    for sys, tdirs in SYSTEMS.items():
        bb_m, bb_s   = stats(tdirs, "rmsd.csv", "backbone_RMSD")
        lig_m, lig_s = stats(tdirs, "rmsd.csv", "ligand_RMSD")
        hb_m, hb_s   = stats(tdirs, "hbonds_count.csv", "n_hbonds")
        gk_m, gk_s, gk_pct = gk_stats(tdirs)
        rows[sys] = dict(bb_m=bb_m, bb_s=bb_s, lig_m=lig_m, lig_s=lig_s,
                         hb_m=hb_m, hb_s=hb_s, gk_m=gk_m, gk_s=gk_s, gk_pct=gk_pct)

    lines = []
    a = lines.append
    sep = "=" * 78

    a(sep)
    a("  4-SYSTEM COMPARISON REPORT — ABL1 T315I resistance panel")
    a("  S1 WT+imatinib (1IEP) | S2 T315I+imatinib | S3 WT+ponatinib (3OXZ) | S4 T315I+ponatinib (3IK3)")
    a(f"  Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}")
    a("  NOTE: MD ensembles generated with CHARMM27 (legacy FF); endpoint MM-GBSA/PBSA")
    a("        not reported here — FF basis inconsistent, values unreliable (see README).")
    a(sep)
    a("")
    a("PER-SYSTEM METRICS (mean ± SD across 3 replicates, all frames)")
    a(f"  {'System':<26}{'bbRMSD(Å)':>12}{'ligRMSD(Å)':>12}{'GKdist(Å)':>12}{'GK<3.5Å%':>10}{'H-bonds':>10}")
    a("  " + "-" * 82)
    for sys, r in rows.items():
        a(f"  {sys:<26}{r['bb_m']:>7.2f}±{r['bb_s']:<4.2f}{r['lig_m']:>7.2f}±{r['lig_s']:<4.2f}"
          f"{r['gk_m']:>7.2f}±{r['gk_s']:<4.2f}{r['gk_pct']:>9.1f}{r['hb_m']:>7.2f}±{r['hb_s']:<3.2f}")
    a("")
    a("─" * 78)
    a("THE TWO CONTRASTS  (mutant − WT;  Δ = destabilisation caused by T315I)")
    a("  The resistance hypothesis: T315I destabilises IMATINIB but not PONATINIB,")
    a("  i.e. Δ(imatinib) should be large, Δ(ponatinib) near zero.")
    a("")
    for drug, wt_key, mut_key in CONTRASTS:
        w, m = rows[wt_key], rows[mut_key]
        a(f"  {drug.upper()}   ({wt_key.split(' (')[0]}  →  {mut_key.split(' (')[0]})")
        a(f"    Δ ligand RMSD     = {m['lig_m']-w['lig_m']:+.2f} Å")
        a(f"    Δ gatekeeper dist = {m['gk_m']-w['gk_m']:+.2f} Å   "
          f"(contact% {w['gk_pct']:.0f} → {m['gk_pct']:.0f})")
        a(f"    Δ H-bonds/frame   = {m['hb_m']-w['hb_m']:+.2f}")
        a(f"    Δ backbone RMSD   = {m['bb_m']-w['bb_m']:+.2f} Å")
        a("")
    a("─" * 78)
    a("RESISTANCE SIGNAL  (does imatinib destabilise MORE than ponatinib?)")
    im = rows["T315I+imatinib"]["lig_m"] - rows["WT+imatinib (1IEP)"]["lig_m"]
    po = rows["T315I+ponatinib (3IK3)"]["lig_m"] - rows["WT+ponatinib (3OXZ)"]["lig_m"]
    a(f"  Δ ligand RMSD:  imatinib {im:+.2f} Å   vs   ponatinib {po:+.2f} Å")
    a(f"  Differential (imatinib − ponatinib) = {im - po:+.2f} Å")
    a("  Positive & larger for imatinib ⇒ structural signature of T315I resistance"
      " reproduced.")
    a("  (Interpret against per-replicate scatter in summary_comparison_<drug>.png.)")
    a("")
    a("FIGURES GENERATED")
    for f in sorted(OUT.glob("*.png")):
        a(f"  {f.name}")
    a(sep)

    report_path = OUT / "comparison_report.txt"
    report_path.write_text("\n".join(lines))
    print("  Saved comparison_report.txt")
    print("\n" + "\n".join(lines))


# ── main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print(f"Output → {OUT}/")
    for drug, systems in PAIRS.items():
        plot_rmsd(drug, systems)
        plot_gatekeeper(drug, systems)
        plot_hbonds(drug, systems)
        plot_rmsf(drug, systems)
        plot_pocket_contacts(drug, systems)
        plot_summary(drug, systems)
    write_report()
    print(f"\nDone. All files in: {OUT}/")
