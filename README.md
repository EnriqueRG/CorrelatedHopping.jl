<h1 align="center">
  <img src="./logo.jpg" alt="CorrelatedHopping.jl" width="720">
</h1>

Julia tools for simulating and analyzing constrained reaction-diffusion models on
periodic lattices.

This package is extracted from the research notebooks in this repository. The
notebooks remain useful as provenance and figure-generation workflows, while
the package in `src/` is the reusable API.

## Quick Start

```julia
using CorrelatedHopping
using Random

rng = MersenneTwister(1)
L = 16
rho0 = 0.5
initial_state = Int.(rand(rng, L) .< rho0)

sys = initialize_system(
    L,
    initial_state,
    1.0,
    1.0;
    reaction = Reaction(2, 0),
    dynamics = CorrelatedHoppingDynamics(),
)

times, particles = simulate!(
    sys,
    (s, _t) -> is_final_binary(s);
    rng,
)
```

For larger production runs, increase `L` and the number of ensemble samples in
your scripts or notebooks. The figure examples use reduced settings compared
with the paper workflow, but some still run for tens of seconds.

## Main Concepts

- `Reaction(n, m)` defines a finite-rate reaction `nA -> mA`.
- `InstantaneousReaction(n, m)` enforces `nA -> mA` immediately after dynamics.
- `CorrelatedHoppingDynamics()` implements the constrained move `0110 <-> 1001`.
- `StandardDiffusion()` implements independent nearest-neighbor diffusion.
- `build_generator` constructs a symmetry-reduced finite-state Markov generator
  for small systems.

## Repository Layout

- `src/`: reusable Julia package code.
- `test/`: automatic correctness checks.
- `examples/`: runnable examples ported from selected `PaperPlots.ipynb` cells.
- `legacy/notebooks/`: research notebooks and figure-generation provenance.
- `legacy/figures/`: generated PDFs from the notebook workflow.
- `legacy/data/`: serialized `.jls` data and raw `.bin` outputs.

## Development

From this directory:

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

The notebooks are not part of the automatic test suite because several cells
are intentionally large research runs.

## Examples

The plotting examples have their own environment so the package itself does not
depend on plotting libraries:

```julia
using Pkg
Pkg.activate("examples")
Pkg.develop(path=".")
Pkg.instantiate()
```

Then run an example from the repository root:

```sh
julia --project=examples examples/fig2.jl
```

- `examples/fig2.jl`: larger PaperPlots density-evolution figure port.
- `examples/fig4.jl`: PaperPlots final-time histogram with Gumbel fits.
- `examples/fig5.jl`: PaperPlots finite-size scaling comparison.
- `examples/fig6.jl`: small-system joint fast/slow-rate final-time panel.
- `examples/fig7.jl`: PaperPlots diffusive final-time scaling figure.
