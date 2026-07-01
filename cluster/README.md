# Arrhenius Cluster Jobs

These Slurm scripts are meant to be submitted from the package root, the
directory containing `Project.toml` and `examples/`.

The scripts currently use the Arrhenius CPU allocation
`naiss2026-3-465-cpu`, load `Julia/1.12.2-bdist`, and request mail on job
completion or failure:

```sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=enriquerozasgarcia@gmail.com
```

If Arrhenius changes the available Julia module, check the current names with:

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
