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
- `fig4_generate_data.jl` and `fig4_plot.jl`: PaperPlots final-time histogram
  with Gumbel fits. Uses `L_values = [2^4, 2^8, 2^12]` and
  `ensemble_size = 500`. Generate the data first, preferably with multiple
  Julia threads:

  ```sh
  julia --threads auto --project=examples examples/fig4_generate_data.jl
  julia --project=examples examples/fig4_plot.jl
  ```
- `fig5.jl`: PaperPlots finite-size scaling comparison with full and reduced
  `ell=12` theory curves.
- `fig6_generate_data.jl` and `fig6_plot.jl`: Small-system version of the joint
  fast/slow-rate `\langle \hat{T}_L\rangle` panel. Generate the data first,
  preferably with multiple Julia threads:

  ```sh
  julia --threads auto --project=examples examples/fig6_generate_data.jl
  julia --project=examples examples/fig6_plot.jl
  ```

- `fig7_generate_data.jl` and `fig7_plot.jl`: PaperPlots diffusive final-time
  scaling figure. Uses `L_values = [2^4, 2^8, 2^12]` and
  `ensemble_size = 200`. Generate the data first, preferably with multiple
  Julia threads:

  ```sh
  julia --threads auto --project=examples examples/fig7_generate_data.jl
  julia --project=examples examples/fig7_plot.jl
  ```
