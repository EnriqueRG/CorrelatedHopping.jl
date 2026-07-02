using CorrelatedHopping
using Dates
using Random
using Serialization
using Statistics

###########################################################
# Fig. 6 data generation                                  #
# Mean final times across density and rate ratios.        #
###########################################################

function fig6_parameters(;
    rng_seed = 6,
    L = 1000,
    hop_rate = 1.0,
    ensemble_size = 10000,
    rho0_values = LinRange(0, 1, 19),
    reaction_over_hop_powers = collect(-3:3),
    n_display_ratios = 4,
)
    reaction_over_hop_powers = collect(reaction_over_hop_powers)
    reaction_over_hop_values = collect(10.0 .^ reaction_over_hop_powers)

    return (;
        rng_seed,
        L,
        hop_rate,
        ensemble_size,
        rho0_values = collect(rho0_values),
        reaction_over_hop_powers,
        reaction_over_hop_values,
        n_display_ratios,
    )
end

"""
Obtain an ensemble of final times for a grid of densities and rate ratios.

Each grid point gets its own RNG stream so the output is reproducible
independently of thread scheduling.
"""
function simulate_rate_grid_parallel(params; verbose = true)
    rho0_values = params.rho0_values
    reaction_over_hop_values = params.reaction_over_hop_values

    mean_times = zeros(length(reaction_over_hop_values), length(rho0_values))
    stderr_times = zeros(length(reaction_over_hop_values), length(rho0_values))
    reaction_rate_values = params.hop_rate .* reaction_over_hop_values

    grid_points = collect(CartesianIndices(mean_times))
    progress_lock = ReentrantLock()
    completed = Ref(0)

    Threads.@threads for k in eachindex(grid_points)
        idx = grid_points[k]
        i, j = Tuple(idx)
        reaction_rate = reaction_rate_values[i]
        rho0 = rho0_values[j]
        cell_seed = params.rng_seed + 1_000_003 * i + 9_176 * j
        rng = MersenneTwister(cell_seed)

        samples = run_ensemble(
            params.L,
            params.ensemble_size,
            rho0;
            hop_rate = params.hop_rate,
            reaction_rate,
            reaction = Reaction(2, 0),
            dynamics = CorrelatedHoppingDynamics(),
            rng,
        )

        mean_times[i, j] = mean(samples)
        stderr_times[i, j] = params.ensemble_size == 1 ? 0.0 : std(samples) / sqrt(params.ensemble_size)

        if verbose
            lock(progress_lock)
            try
                completed[] += 1
                println(
                    "Completed ",
                    completed[],
                    "/",
                    length(grid_points),
                    ": reaction_rate/hop_rate=",
                    reaction_over_hop_values[i],
                    ", rho0=",
                    rho0,
                )
                flush(stdout)
            finally
                unlock(progress_lock)
            end
        end
    end

    return (; mean_times, stderr_times, reaction_rate_values)
end

"""
Obtain the amplitude c1 of the slowest decay mode as a function of density
rho0 in the reaction_rate -> infinity limit.
"""
function c1_curve()
    L = 12
    hop_rate = 1.0
    reaction_rate = 10.0 # Close enough to the reaction_rate -> infinity limit for L = 12.
    rho0_values = collect(LinRange(0, 1, 100))
    W, representatives, multiplicities = build_pair_annihilation_generator(
        L;
        hop_rate,
        reaction_rate,
        lump_final = true,
    )
    _, right_vectors, left_vectors = find_spectral_gap(W; num_eigenvalues = 2)
    c1_values = [
        CorrelatedHopping.calculate_coefficient(
            rho0,
            1,
            right_vectors,
            left_vectors,
            representatives,
            multiplicities,
            L,
        ) for rho0 in rho0_values
    ]
    return (; rho0_values, c1_values, L, hop_rate, reaction_rate)
end

function generate_fig6_data(;
    verbose = true,
    kwargs...,
)
    params = fig6_parameters(; kwargs...)
    metadata = (;
        generated_at = string(now()),
        julia_version = string(VERSION),
        nthreads = Threads.nthreads(),
        reaction = "Reaction(2, 0)",
        dynamics = "CorrelatedHoppingDynamics()",
    )

    verbose && println(
        "Generating fig6 data with ",
        Threads.nthreads(),
        " Julia thread(s); L=",
        params.L,
        ", ensemble_size=",
        params.ensemble_size,
    )
    simulation = simulate_rate_grid_parallel(params; verbose)
    verbose && println("Calculating inset c1 curve.")
    inset = c1_curve()

    data = (; metadata, parameters = params, simulation, inset)
    output_dir = joinpath(@__DIR__, "output")
    mkpath(output_dir)
    path = joinpath(output_dir, "fig6_data.jls")
    serialize(path, data)
    println("Saved fig6 data to ", path)
    return data
end

generate_fig6_data();
