using CorrelatedHopping
using Dates
using Random
using Serialization

###########################################################
# Fig. 7 data generation                                  #
# Diffusive final-time scaling histograms.                #
###########################################################

function fig7_parameters(;
    rng_seed = 7,
    rho0 = 0.5,
    hop_rate = 1.0,
    reaction_rate = 1.0,
    ensemble_size = 10_000, # Reduce this parameter for a quick test run
    L_values = [2^4, 2^8, 2^12],
)
    return (;
        rng_seed,
        rho0,
        hop_rate,
        reaction_rate,
        ensemble_size,
        L_values = collect(L_values),
    )
end

function chunk_sizes(n, n_chunks)
    base, remainder = divrem(n, n_chunks)
    return [base + (i <= remainder ? 1 : 0) for i in 1:n_chunks if base + (i <= remainder ? 1 : 0) > 0]
end

function run_ensemble_parallel(
    L,
    ensemble_size,
    rho0;
    hop_rate,
    reaction_rate,
    rng_seed,
    reaction = Reaction(2, 0),
    dynamics = StandardDiffusion(),
    condition_even = true,
    stop_on_reaction_only = true,
    stop = nothing,
)
    sizes = chunk_sizes(ensemble_size, min(Threads.nthreads(), ensemble_size))
    chunks = Vector{Vector{Float64}}(undef, length(sizes))

    Threads.@threads for chunk_index in eachindex(sizes)
        rng = MersenneTwister(rng_seed + 1_000_003 * chunk_index)
        chunks[chunk_index] = run_ensemble(
            L,
            sizes[chunk_index],
            rho0;
            hop_rate,
            reaction_rate,
            reaction,
            dynamics,
            condition_even,
            stop_on_reaction_only,
            stop,
            rng,
        )
    end

    return vcat(chunks...)
end

function generate_fig7_data(; verbose = true, kwargs...)
    params = fig7_parameters(; kwargs...)
    verbose && println(
        "Generating fig7 data with ",
        Threads.nthreads(),
        " Julia thread(s); ensemble_size=",
        params.ensemble_size,
    )

    results = [
        begin
            final_times = run_ensemble_parallel(
                L,
                params.ensemble_size,
                params.rho0;
                hop_rate = params.hop_rate,
                reaction_rate = params.reaction_rate,
                rng_seed = params.rng_seed + 100_003 * i,
                reaction = Reaction(2, 0),
                dynamics = StandardDiffusion(),
                condition_even = true,
                stop = (s, _t) -> sum(s.occupations) < 2,
                stop_on_reaction_only = true,
            )
            verbose && println("Completed L=", L)
            (; L, final_times)
        end
        for (i, L) in enumerate(params.L_values)
    ]

    metadata = (;
        generated_at = string(now()),
        julia_version = string(VERSION),
        nthreads = Threads.nthreads(),
        reaction = "Reaction(2, 0)",
        dynamics = "StandardDiffusion()",
    )
    data = (; metadata, parameters = params, results)

    output_dir = joinpath(@__DIR__, "output")
    mkpath(output_dir)
    path = joinpath(output_dir, "fig7_data.jls")
    serialize(path, data)
    println("Saved fig7 data to ", path)
    return data
end

generate_fig7_data();
