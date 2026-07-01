#!/bin/bash
#SBATCH -J ch_fig8
#SBATCH -A naiss2026-3-465-cpu
#SBATCH -p fat
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --mem=256G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=enriquerozasgarcia@gmail.com
#SBATCH -o logs/fig8-%j.out
#SBATCH -e logs/fig8-%j.err

set -euo pipefail

ml purge
ml Julia/1.12.2-bdist

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs examples/output

export JULIA_NUM_THREADS=1
export JULIA_DEPOT_PATH="${PWD}/.julia_depot:${JULIA_DEPOT_PATH:-}"
export FIG8_PROFILE=cluster

srun julia --project=examples --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
srun julia --project=examples --startup-file=no examples/fig8_generate_data.jl
