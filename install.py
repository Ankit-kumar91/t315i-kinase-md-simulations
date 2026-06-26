#!/usr/bin/env python3
"""
install.py  --  Install all packages needed for the MD analysis pipeline.

Run this ONCE before your first analysis:

    python install.py

After it finishes, run the analysis with:

    python analyze_md.py
"""

import sys
import subprocess

PACKAGES = [
    "numpy",
    "pandas",
    "matplotlib",
    "seaborn",
    "MDAnalysis",
    "fpdf2",          # used by generate_manual.py to create the PDF guide
]

print()
print("=" * 55)
print("  MD Analysis Pipeline  --  Installing packages")
print("=" * 55)
print(f"  Python: {sys.executable}")
print()

all_ok = True
for pkg in PACKAGES:
    print(f"  Installing {pkg} ...", end=" ", flush=True)
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", "--quiet", pkg],
        capture_output=True, text=True,
    )
    if result.returncode == 0:
        print("OK")
    else:
        print("FAILED")
        print(f"    {result.stderr.strip()}")
        all_ok = False

print()
print("=" * 55)
if all_ok:
    print("  All packages installed successfully!")
    print()
    print("  Next step:")
    print("    python analyze_md.py")
else:
    print("  Some packages failed. Check the errors above.")
    print("  Try running:  pip install <package_name>")
print("=" * 55)
print()
