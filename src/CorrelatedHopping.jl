module CorrelatedHopping

using Arpack
using LinearAlgebra
using Printf
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
    run_ensemble,
    initialize_active_sites,
    search_active_sites,
    get_sandboxes,
    sandbox_search,
    sandbox_search_history,
    symmetry_sector_counts,
    largest_symmetry_sector,
    generate_symmetry_sector_states,
    analyze_krylov_sectors!,
    largest_krylov_sector

include("gillespie.jl")
include("spectral.jl")
include("fragments.jl")
include("krylov.jl")

end
