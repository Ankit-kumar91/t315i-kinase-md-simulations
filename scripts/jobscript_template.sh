#!/usr/bin/env bash
#SBATCH --job-name=MY_JOB_NAME
#SBATCH --time=00:20:00
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --output=sl_%j.out
#SBATCH --error=sl_%j.err
#SBATCH --gres=gpu:a100:1

module reset
module load bio/GROMACS/2025.2-foss-2024a-CUDA-12.6.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=true

srun gmx mdrun -s topol.tpr -nb gpu
