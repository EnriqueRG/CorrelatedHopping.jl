module CorrelatedHopping

using Arpack
using LinearAlgebra
using Random
using SparseArrays
using Statistics

export 
    AbstractReactionModel,
    Reaction,
    InstantaneousReaction,
    AbstractDynamicsModel,
    CorrelatedHoppingDynamics,
    StandardDiffusion,
    DirectLatticeSystem,
    initialize_system,
    update_local_rates!,
    simulate!,
    is_final_binary,
    build_generator,
    find_spectral_gap,
    full_gamma_coefficients,
    calculate_coefficient,
    run_ensemble

include("gillespie.jl")
include("spectral.jl")
include("fragments.jl")

end
