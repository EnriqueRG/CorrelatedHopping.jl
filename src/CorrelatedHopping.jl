module CorrelatedHopping

using Arpack
using Printf
using Random
using SparseArrays

export 
    Reaction,
    InstantaneousReaction,
    CorrelatedHoppingDynamics,
    StandardDiffusion,
    DirectLatticeSystem,
    initialize_system,
    simulate!,
    is_final_binary,
    build_pair_annihilation_generator,
    find_spectral_gap,
    full_decay_coefficients,
    run_ensemble,
    active_sites

include("gillespie.jl")
include("spectral.jl")
include("fragments.jl")
include("krylov.jl")

end
