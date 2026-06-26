#!/usr/bin/env bash
# Run NVT -> NPT -> 100ns production for Replica1.
# Resumable: re-running this script skips steps already finished and
# continues an interrupted production run from its .cpt checkpoint.
set -euo pipefail
cd "$(dirname "$0")"

log() { echo "[$(date +%T)] $*"; }

# --- NVT ---
if [ -f nvt.gro ]; then
  log "NVT already done, skipping"
else
  log "Running NVT"
  gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
  gmx mdrun -deffnm nvt -nb gpu -pme gpu -bonded gpu
fi

# --- NPT ---
if [ -f npt.gro ]; then
  log "NPT already done, skipping"
else
  log "Running NPT"
  gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
  gmx mdrun -deffnm npt -nb gpu -pme gpu -bonded gpu
fi

# --- Production (100 ns) ---
if [ -f md.gro ]; then
  log "Production already done"
elif [ -f md.cpt ]; then
  log "Resuming production from checkpoint"
  gmx mdrun -deffnm md -cpi md.cpt -nb gpu -pme gpu -bonded gpu
else
  log "Starting production"
  gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
  gmx mdrun -deffnm md -nb gpu -pme gpu -bonded gpu
fi

log "Replica1 done"
