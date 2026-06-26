#!/usr/bin/env python3
# generate_manual.py — produces USER_MANUAL.pdf

import sys
from pathlib import Path
from fpdf import FPDF, XPos, YPos

OUT      = Path(__file__).resolve().parent / "USER_MANUAL.pdf"
FONT_DIR = Path(sys.executable).resolve().parent.parent / "fonts"

# colours
C_DARK    = (30,  30,  30)
C_BLUE1   = (26,  82, 118)
C_BLUE2   = (40, 116, 166)
C_CODE_FG = (44,  44,  44)
C_CODE_BG = (245, 245, 245)
C_THEAD   = (26,  82, 118)
C_TROW_A  = (240, 245, 250)
C_TROW_B  = (255, 255, 255)
C_RULE    = (180, 180, 180)
C_WHITE   = (255, 255, 255)


class PDF(FPDF):
    def __init__(self):
        super().__init__()
        self.set_auto_page_break(auto=True, margin=22)
        self.set_margins(20, 20, 20)
        self._add_fonts()

    def _add_fonts(self):
        fd = FONT_DIR
        self.add_font("Body",  style="",  fname=str(fd / "Ubuntu-R.ttf"))
        self.add_font("Body",  style="B", fname=str(fd / "Ubuntu-B.ttf"))
        self.add_font("Body",  style="I", fname=str(fd / "Ubuntu-RI.ttf"))
        self.add_font("Body",  style="BI",fname=str(fd / "Ubuntu-BI.ttf"))
        self.add_font("Mono",  style="",  fname=str(fd / "SourceCodePro-Regular.ttf"))
        self.add_font("Mono",  style="B", fname=str(fd / "SourceCodePro-Bold.ttf"))

    def header(self):
        if self.page_no() == 1:
            return
        self.set_font("Body", "I", 8)
        self.set_text_color(*C_RULE)
        self.cell(0, 5, "MD Analysis Pipeline — User Manual",
                  new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_draw_color(*C_RULE)
        self.line(20, self.get_y(), 190, self.get_y())
        self.ln(3)
        self.set_text_color(*C_DARK)

    def _f(self, name, style, size):
        self.set_font(name, style, size)

    def h1(self, text):
        self.ln(5)
        self._f("Body", "B", 17)
        self.set_text_color(*C_BLUE1)
        self.cell(0, 10, text, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_draw_color(*C_BLUE1)
        self.set_line_width(0.6)
        self.line(20, self.get_y(), 190, self.get_y())
        self.set_line_width(0.2)
        self.ln(4)
        self.set_text_color(*C_DARK)

    def h2(self, text):
        self.ln(4)
        self._f("Body", "B", 12)
        self.set_text_color(*C_BLUE2)
        self.cell(0, 7, text, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(1)
        self.set_text_color(*C_DARK)

    def h3(self, text):
        self.ln(3)
        self._f("Body", "B", 10)
        self.set_text_color(*C_DARK)
        self.cell(0, 6, text, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    def body(self, text):
        self._f("Body", "", 10)
        self.set_text_color(*C_DARK)
        self.multi_cell(0, 5.5, text)
        self.ln(1)

    def bullet(self, items):
        self._f("Body", "", 10)
        self.set_text_color(*C_DARK)
        for item in items:
            x = self.get_x()
            self.set_x(26)
            self.cell(5, 5.5, "•")
            self.multi_cell(0, 5.5, item)
        self.ln(1)

    def code(self, text):
        self.ln(1)
        lines = text.strip().split("\n")
        h = len(lines) * 5.5 + 5
        y = self.get_y()
        self.set_fill_color(*C_CODE_BG)
        self.set_draw_color(*C_RULE)
        self.rect(20, y, 170, h, style="FD")
        self.set_xy(25, y + 2.5)
        self._f("Mono", "", 8.5)
        self.set_text_color(*C_CODE_FG)
        for line in lines:
            self.set_x(25)
            self.cell(0, 5.5, line, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(2)
        self.set_text_color(*C_DARK)

    def hrule(self):
        self.ln(3)
        self.set_draw_color(*C_RULE)
        self.line(20, self.get_y(), 190, self.get_y())
        self.ln(4)

    def table(self, headers, rows, col_widths):
        self.ln(2)
        # header row
        self._f("Body", "B", 9)
        self.set_fill_color(*C_THEAD)
        self.set_text_color(*C_WHITE)
        for h, w in zip(headers, col_widths):
            self.cell(w, 7, h, border=1, fill=True)
        self.ln()
        # data rows
        self._f("Body", "", 9)
        for i, row in enumerate(rows):
            fill_col = C_TROW_A if i % 2 == 0 else C_TROW_B
            self.set_fill_color(*fill_col)
            self.set_text_color(*C_DARK)
            # measure max height needed
            x0, y0 = self.get_x(), self.get_y()
            max_h = 6
            for cell, w in zip(row, col_widths):
                lines_n = self.get_string_width(cell) / (w - 2) + 1
                max_h = max(max_h, int(lines_n) * 6)
            x = x0
            for cell, w in zip(row, col_widths):
                self.set_xy(x, y0)
                self.multi_cell(w, 6, cell, border=1, fill=True)
                x += w
            self.set_xy(x0, y0 + max_h)
        self.ln(3)


def build():
    pdf = PDF()

    # ── Cover ────────────────────────────────────────────────────────────────
    pdf.add_page()
    pdf.set_fill_color(*C_BLUE1)
    pdf.rect(0, 0, 210, 75, style="F")
    pdf.set_y(20)
    pdf.set_font("Body", "B", 26)
    pdf.set_text_color(*C_WHITE)
    pdf.cell(0, 13, "MD Analysis Pipeline", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
    pdf.set_font("Body", "", 16)
    pdf.cell(0, 10, "User Manual", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
    pdf.set_font("Body", "I", 10)
    pdf.cell(0, 8, "Automated GROMACS MD Analysis - Installation, Configuration, and Usage",
             new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")

    pdf.set_y(85)
    pdf.set_text_color(*C_DARK)
    pdf.set_font("Body", "B", 11)
    pdf.set_text_color(*C_BLUE1)
    pdf.cell(0, 8, "Files included", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(1)
    pdf.table(
        ["File", "Purpose"],
        [
            ["md_analysis.py",  "Main analysis script - never needs to be edited"],
            ["md_config.ini",   "Configuration file - edit for each simulation"],
            ["install.py",      "One-time dependency installer"],
            ["USER_MANUAL.pdf", "This document"],
        ],
        [55, 115],
    )
    pdf.set_text_color(*C_DARK)
    pdf.set_font("Body", "B", 11)
    pdf.set_text_color(*C_BLUE1)
    pdf.cell(0, 8, "Quick Start", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(1)
    pdf.set_text_color(*C_DARK)
    pdf.table(
        ["Step", "Command"],
        [
            ["1. Install dependencies (once)", "python install.py"],
            ["2. Edit md_config.ini",          "Set sim_dir = /path/to/your/sim/folder"],
            ["3. Run analysis",                "python md_analysis.py /path/to/sim/folder"],
        ],
        [75, 95],
    )

    # ── Page 2 ───────────────────────────────────────────────────────────────
    pdf.add_page()

    pdf.h1("1.  Install Dependencies")
    pdf.body(
        "Run this once on any new computer before your first analysis. "
        "It installs all required Python packages automatically."
    )
    pdf.code("python install.py")
    pdf.table(
        ["Package", "Used for"],
        [
            ["NumPy",       "Numerical calculations"],
            ["Pandas",      "Data tables and CSV export"],
            ["Matplotlib",  "Plotting"],
            ["Seaborn",     "Plot styling"],
            ["MDAnalysis",  "Reading trajectories - RMSD, RMSF, H-bonds, contacts"],
            ["fpdf2",       "PDF generation (this manual)"],
        ],
        [40, 130],
    )
    pdf.body(
        "GROMACS must also be available on your system. "
        "On HPC clusters, load it before running the analysis:"
    )
    pdf.code("module load gromacs")

    pdf.hrule()

    pdf.h1("2.  Configure md_config.ini")
    pdf.body(
        "Open md_config.ini in any text editor (Notepad on Windows, "
        "TextEdit on Mac, gedit or nano on Linux). "
        "The only required change for each simulation is sim_dir."
    )

    pdf.h2("Required - set the simulation folder")
    pdf.body("Change sim_dir and out_dir to point at your GROMACS files:")
    pdf.code(
        "[paths]\n"
        "sim_dir = /home/user/projects/1IEP_MD\n"
        "out_dir = /home/user/projects/1IEP_MD/results"
    )
    pdf.body(
        "Use the full absolute path. "
        "out_dir is created automatically if it does not exist."
    )

    pdf.h2("File names are auto-detected")
    pdf.body(
        "The script scans sim_dir and picks .tpr, .xtc, .edr, .gro, and .ndx files "
        "automatically, preferring files whose name contains 'md' (skipping em, npt, nvt). "
        "Uncomment and set file names manually only if detection picks the wrong file:"
    )
    pdf.code(
        "[paths]\n"
        "tpr = md_production.tpr\n"
        "xtc = md_production.xtc"
    )

    pdf.h2("Common optional settings")

    pdf.h3("Change the time axis to nanoseconds")
    pdf.code("[analysis]\ntime_unit = ns")

    pdf.h3("Add a label to all plot titles")
    pdf.code("[system]\nsim_label = ABL1 WT + Imatinib 100 ns")

    pdf.h3("Set the gatekeeper residue (0 = skip)")
    pdf.code("[system]\ngatekeeper_res = 315")

    pdf.h3("Set ligand residue name (if auto-detection is wrong)")
    pdf.code("[system]\nlig_resname = LIG")

    pdf.h3("Use a different GROMACS command (e.g. on HPC clusters)")
    pdf.code("[system]\ngromacs_cmd = gmx_mpi")

    pdf.h3("Adjust plot appearance")
    pdf.code(
        "[figure]\n"
        "dpi            = 300       # 72=screen  150=standard  300=publication\n"
        "width          = 12\n"
        "height         = 5\n"
        "title_fontsize = 14\n"
        "label_fontsize = 12\n"
        "ligand_color   = #ff6600\n"
        "time_unit      = ns"
    )

    pdf.add_page()

    pdf.h1("3.  Run the Analysis")
    pdf.body(
        "Pass the simulation folder as an argument. "
        "The script looks for md_config.ini inside that folder automatically:"
    )
    pdf.code("python md_analysis.py /path/to/1IEP_MD")
    pdf.body("Or pass the config file directly:")
    pdf.code("python md_analysis.py --config /path/to/1IEP_MD/md_config.ini")

    pdf.h2("Running multiple simulations")
    pdf.body(
        "Keep one copy of md_analysis.py. Place a copy of md_config.ini in each "
        "simulation folder and change only sim_dir in each one:"
    )
    pdf.code(
        "python md_analysis.py /data/1IEP_MD\n"
        "python md_analysis.py /data/3IK3_MD\n"
        "python md_analysis.py /data/T315I_MD"
    )
    pdf.body(
        "Each run writes its results to its own out_dir, so the outputs never mix."
    )

    pdf.hrule()

    pdf.h1("4.  What the Script Does")
    pdf.body(
        "The pipeline runs ten steps in order. "
        "Each step is independent - if one fails, the rest continue."
    )

    steps = [
        ("Step 1 - PBC Correction",
         "Fixes periodic-boundary artefacts using GROMACS trjconv. Produces a corrected "
         "trajectory (md_center.xtc) used by all subsequent steps. "
         "Skipped automatically if the corrected trajectory already exists."),
        ("Step 2 - Thermodynamic Properties",
         "Extracts potential energy, kinetic energy, total energy, temperature, pressure, "
         "volume, and density from the .edr file using gmx energy. "
         "Each property is saved as an individual plot with mean +/- 1 sigma shading."),
        ("Step 3 - RMSD",
         "Root-mean-square deviation of the protein backbone Ca atoms (structural "
         "stability) and the ligand heavy atoms fitted on the backbone (binding pose "
         "stability). Values consistently below 2-3 Angstroms indicate a stable simulation."),
        ("Step 4 - RMSF",
         "Per-residue Ca root-mean-square fluctuation. Shows which parts of the protein "
         "are flexible. The gatekeeper residue is highlighted with a vertical marker."),
        ("Step 5 - Radius of Gyration",
         "Protein compactness over time using gmx gyrate. "
         "A stable Rg indicates the protein maintains its folded conformation."),
        ("Step 6 - Protein-Ligand H-bonds",
         "Counts hydrogen bonds between protein and ligand at every frame. "
         "Reports per-pair occupancy (% of frames each H-bond is present)."),
        ("Step 7 - Gatekeeper Distance",
         "Minimum distance between the ligand and the gatekeeper residue (e.g. THR315 "
         "in ABL1) over the trajectory. Shows whether the ligand maintains contact with "
         "this key residue throughout the simulation."),
        ("Step 8 - Binding-Pocket Contact Map",
         "Identifies all protein residues within the contact cutoff (default 4.5 A) "
         "of the ligand and reports their occupancy as a bar chart. Residues below "
         "min_plot_contact_pct are hidden."),
        ("Step 9 - Summary Figure",
         "A four-panel figure combining backbone RMSD, ligand RMSD, radius of gyration, "
         "and per-residue RMSF in one image."),
        ("Step 10 - Text Report",
         "A plain-text summary of all key numerical results saved to "
         "tables/analysis_report.txt."),
    ]
    for title, desc in steps:
        pdf.h3(title)
        pdf.body(desc)

    pdf.add_page()

    pdf.h1("5.  Output Files")
    pdf.body(
        "All output is written to out_dir (default: results/ inside the simulation folder)."
    )
    pdf.table(
        ["File", "Contents"],
        [
            ["energy_potential.png",            "Potential energy vs time"],
            ["energy_kinetic.png",              "Kinetic energy vs time"],
            ["energy_total_energy.png",         "Total energy vs time"],
            ["energy_temperature.png",          "Temperature vs time"],
            ["energy_pressure.png",             "Pressure vs time"],
            ["energy_volume.png",               "Volume vs time"],
            ["energy_density.png",              "Density vs time"],
            ["rmsd.png",                        "Backbone and ligand RMSD (two panels)"],
            ["rmsf.png",                        "Per-residue Ca RMSF"],
            ["gyration.png",                    "Radius of gyration"],
            ["hbonds.png",                      "H-bond count vs time"],
            ["gatekeeper_distance.png",         "Ligand-gatekeeper minimum distance"],
            ["pocket_contacts.png",             "Binding-pocket contact occupancy bar chart"],
            ["summary_panel.png",               "Four-panel summary figure"],
            ["tables/rmsd.csv",                 "RMSD data (time, backbone, ligand)"],
            ["tables/rmsf.csv",                 "RMSF data (resid, resname, rmsf_A)"],
            ["tables/hbonds_count.csv",         "H-bond count per frame"],
            ["tables/hbonds_pairs.csv",         "Per-pair H-bond occupancy table"],
            ["tables/gatekeeper_distance.csv",  "Gatekeeper min-distance per frame"],
            ["tables/pocket_contacts.csv",      "Contact occupancy per residue"],
            ["tables/analysis_report.txt",      "Full numerical summary report"],
        ],
        [78, 92],
    )

    pdf.hrule()

    pdf.h1("6.  Troubleshooting")
    problems = [
        ("GROMACS not found",
         "Either GROMACS is not installed or not on your PATH.\n"
         "On HPC: run  module load gromacs  before the analysis.\n"
         "Custom name: set  gromacs_cmd = gmx_mpi  in [system] of md_config.ini."),
        ("Missing files in sim_dir",
         "Check that sim_dir is the correct folder. The script prints what it found.\n"
         "If your file names are non-standard, set them explicitly in [paths]."),
        ("Multiple non-protein residues found",
         "The script used the first one. If wrong, set the correct residue name:\n"
         "    [system]\n    lig_resname = LIG"),
        ("Protein_LIG group NOT FOUND",
         "The script tries to create this group automatically with gmx make_ndx.\n"
         "If that also fails, find the correct group number in the startup table and add:\n"
         "    [groups]\n    Protein_LIG = 22"),
        ("MDAnalysis steps are skipped",
         "MDAnalysis is not installed. Run  python install.py  to fix this."),
        ("Plots show 'No data'",
         "That step's analysis failed silently. Check the terminal output for "
         "a [WARN] message on that step for the specific error."),
    ]
    for title, desc in problems:
        pdf.h3(title)
        pdf.body(desc)

    pdf.add_page()

    pdf.h1("7.  Quick-Reference - All Settings")
    pdf.table(
        ["Section", "Key", "Default", "Description"],
        [
            ["[paths]",    "sim_dir",              "(required)", "Folder with GROMACS files"],
            ["[paths]",    "out_dir",              "results",    "Output folder (auto-created)"],
            ["[system]",   "sim_label",            "simulation", "Label in plot titles"],
            ["[system]",   "gromacs_cmd",          "gmx",        "GROMACS executable name"],
            ["[system]",   "lig_resname",          "auto",       "Ligand residue name"],
            ["[system]",   "gatekeeper_res",       "315",        "Gatekeeper residue (0=skip)"],
            ["[analysis]", "time_unit",            "ps",         "Time axis unit: ps or ns"],
            ["[analysis]", "contact_cutoff_A",     "4.5",        "Contact distance cutoff (A)"],
            ["[analysis]", "hbond_distance_A",     "3.5",        "H-bond donor-acceptor cutoff (A)"],
            ["[analysis]", "hbond_angle_deg",      "120",        "H-bond angle cutoff (degrees)"],
            ["[analysis]", "min_plot_contact_pct", "10",         "Min occupancy % shown in contact map"],
            ["[figure]",   "dpi",                  "150",        "Image resolution (DPI)"],
            ["[figure]",   "width / height",       "10 / 4",     "Figure size in inches"],
            ["[figure]",   "title_fontsize",       "13",         "Plot title font size (pt)"],
            ["[figure]",   "label_fontsize",       "11",         "Axis label font size (pt)"],
            ["[figure]",   "font_scale",           "1.1",        "Overall text scale multiplier"],
            ["[figure]",   "backbone_color",       "#2c7bb6",    "Protein backbone trace colour"],
            ["[figure]",   "ligand_color",         "#d7191c",    "Ligand trace colour"],
            ["[figure]",   "hbond_color",          "#1a9641",    "H-bond trace colour"],
            ["[figure]",   "gatekeeper_color",     "#fd8d3c",    "Gatekeeper trace colour"],
            ["[figure]",   "rmsf_color",           "#756bb1",    "RMSF trace colour"],
            ["[figure]",   "energy_color",         "#636363",    "Energy trace colour"],
            ["[figure]",   "ligand_marker_size",   "60",         "Marker size on RMSF plot (pt2)"],
        ],
        [28, 44, 26, 72],
    )

    pdf.output(str(OUT))
    print(f"Written: {OUT}")


if __name__ == "__main__":
    build()
