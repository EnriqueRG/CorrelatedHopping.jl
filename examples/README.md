# Examples

These scripts are runnable plot-generating versions of selected workflows from
`legacy/notebooks/PaperPlots.ipynb`.

They use reduced sizes compared with the full paper workflow, but several still
run for seconds to minutes depending on the machine. Before increasing `L`,
`n_ens`, or `n_ensemble`, expect the runtime to grow quickly.

From the repository root, set up the example environment with:

```julia
using Pkg
Pkg.activate("examples")
Pkg.develop(path=".")
Pkg.instantiate()
```

Then run a figure script with:

```sh
julia --project=examples examples/fig2.jl
```

Generated PDFs are written to `examples/output/`.

## Figure Examples

- `fig2.jl`: larger PaperPlots density-evolution figure port. Uses `L = 2^12`
  and `n_traj = 10`.
- `fig4.jl`: PaperPlots final-time histogram with Gumbel fits. Uses
  `L_vals = [2^4, 2^8, 2^12]` and `n_ens = 100`.
- `fig5.jl`: PaperPlots finite-size scaling comparison with full and reduced
  `ell=12` theory curves.
- `fig6.jl`: Small-system version of the joint fast/slow-rate
  `\langle \hat{T}_L\rangle` panel.
- `fig7.jl`: PaperPlots diffusive final-time scaling figure. Uses
  `L_values = [2^4, 2^8, 2^12]` and `n_ens = 100`.
