<h1 align="center">
  <img src="./logo.jpg" alt="CorrelatedHopping.jl" width="720">
</h1>

Julia tools for simulating and analyzing constrained reaction-diffusion models on
periodic lattices.

This package contains reusable Julia code developped for the simulations in [Physical Review Research 6, L032071 (2024)](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.L032071)
and a companion manuscript.

## Installation

From Julia, install the package from GitHub with:

```julia
using Pkg
Pkg.add(url="https://github.com/EnriqueRG/CorrelatedHopping.jl.git")
```

For local development from this repository, activate the project and run the
test suite:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.test()
```

## Quick Start

```julia
using CorrelatedHopping
using Random

# Parameters
rng = MersenneTwister(1)
L = 16
initial_density = 0.5
hopping_rate  = 1.0
reaction_rate = 1.0

# Bernoulli initial state
initial_state = Int.(rand(rng, L) .< initial_density)

# Simulation
sys = initialize_system(
    L,
    initial_state,
    hopping_rate,
    reaction_rate;
    reaction = Reaction(2, 0),
    dynamics = CorrelatedHoppingDynamics(),
)

simulation = simulate!(
    sys,
    (s, _t) -> is_final_binary(s);
    rng,
)

times = simulation.recorded_times
particles = simulation.history
```

This example uses the Gillespie algorithm and records `times` and the
corresponding total particle count `particles`. The absorption time is available
as `simulation.final_time`.

For larger production runs, increase `L` and the number of ensemble samples in
your scripts. The figure examples use reduced settings compared with the paper
workflow, but some still run for tens of seconds.

## Main Concepts

- `Reaction(n, m)` defines a finite-rate reaction `nA -> mA` with `n >= 1`
  and `m >= 0`; finite reactions may remove or create particles.
- `InstantaneousReaction(n, m)` enforces `nA -> mA` immediately after dynamics
  and requires `n >= 1`, `m >= 0`, and `m < n`.
- `CorrelatedHoppingDynamics()` implements the constrained move `0110 <-> 1001`.
- `StandardDiffusion()` implements independent nearest-neighbor diffusion.
- `build_pair_annihilation_generator` constructs a symmetry-reduced finite-state
  Markov generator for the small-system `2A -> 0` correlated-hopping model.

## Repository Layout

- `src/`: reusable Julia package code.
- `test/`: automatic correctness checks.
- `examples/`: runnable plotting examples with reduced settings.

## Examples

The examples reproduce the figures in the companion manuscript. The plotting examples have their own environment so the package itself does not depend on plotting libraries:

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
