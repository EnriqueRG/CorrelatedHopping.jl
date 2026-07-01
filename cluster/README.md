# Arrhenius Cluster Jobs

These Slurm scripts are meant to be submitted from the package root, the
directory containing `Project.toml` and `examples/`.

Before submitting, edit `#SBATCH -A naissYYYY-X-XX` in each script to your NAISS
project allocation, and adjust the Julia module line if Arrhenius reports a
different version:

```sh
ml avail julia
```

Set up the repository once on the cluster:

```sh
git clone <repo-url>
cd <repo-directory>
```

If you want Fig. 8 to resume from data already produced locally, copy the
current CSV before submitting:

```sh
mkdir -p examples/output
# copy fig8_data.csv into examples/output/fig8_data.csv
```

Submit the data jobs:

```sh
sbatch cluster/fig6_data.sbatch
sbatch cluster/fig8_data.sbatch
```

Monitor them with:

```sh
squeue -u "$USER"
tail -f logs/fig6-<jobid>.out
tail -f logs/fig8-<jobid>.out
```

`fig6_data.sbatch` uses multiple Julia threads through `SLURM_CPUS_PER_TASK`.
`fig8_data.sbatch` is single-threaded but requests a large-memory node and sets
`FIG8_PROFILE=cluster`, so it includes the additional cluster-only Fig. 8 data
points.
