#!/bin/bash
#SBATCH -J ch_fig6
#SBATCH -A naiss2026-3-465-cpu
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=enriquerozasgarcia@gmail.com
#SBATCH -o logs/fig6-%j.out
#SBATCH -e logs/fig6-%j.err

set -euo pipefail

ml purge
ml Julia/1.12.2-bdist

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs examples/output

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export JULIA_DEPOT_PATH="${PWD}/.julia_depot:${JULIA_DEPOT_PATH:-}"

srun julia --project=examples --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
srun julia --threads="${JULIA_NUM_THREADS}" --project=examples --startup-file=no examples/fig6_generate_data.jl
