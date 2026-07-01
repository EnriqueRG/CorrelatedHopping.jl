#!/bin/bash
#SBATCH -J ch_fig8
#SBATCH -A naissYYYY-X-XX
#SBATCH -p fat
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --mem=256G
#SBATCH -o logs/fig8-%j.out
#SBATCH -e logs/fig8-%j.err

set -euo pipefail

# Edit the module names after checking `ml avail julia` on Arrhenius.
ml purge
ml julia/1.12.1

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs examples/output

export JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="${PWD}/.julia_depot:${JULIA_DEPOT_PATH:-}"
export FIG8_PROFILE=cluster

srun julia --project=examples --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
srun julia --project=examples --startup-file=no examples/fig8_generate_data.jl
