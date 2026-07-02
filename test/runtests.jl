using CorrelatedHopping
using LinearAlgebra
using Random
using SparseArrays
using Test

@testset "reactions and final states" begin
    sys = initialize_system(
        8,
        [2, 0, 0, 0, 0, 0, 0, 0],
        1.0,
        1.0;
        reaction = InstantaneousReaction(2, 0),
    )
    @test sys.occupations == zeros(Int, 8)
    @test is_final_binary([1, 0, 0, 0, 0, 1, 0, 0])
    @test !is_final_binary([1, 1, 1, 0, 0, 0, 0, 0])
    @test !is_final_binary([2, 0, 0, 0, 0, 0, 0, 0])
end

@testset "system initialization" begin
    sys = initialize_system(
        8,
        Int8[1, 0, 0, 0, 0, 0, 0, 0],
        1,
        1;
        dynamics = :pairwise,
    )
    @test sys.occupations isa Vector{Int}
    @test sum(sys.occupations) == 1

    diffusion_sys = initialize_system(
        4,
        [1, 0, 0, 0],
        2,
        1;
        dynamics = :diffusion,
    )
    @test diffusion_sys.site_rates == [4.0, 0.0, 0.0, 0.0]
    @test diffusion_sys.tree[1] == 4.0

    @test_throws ArgumentError initialize_system(3, zeros(Int, 3); dynamics = CorrelatedHoppingDynamics())
    @test_throws ArgumentError initialize_system(8, [-1, 0, 0, 0, 0, 0, 0, 0])
    @test_throws ArgumentError initialize_system(8, zeros(Int, 8), -1, 1)
end

@testset "simulation" begin
    rng = MersenneTwister(1)
    sys = initialize_system(
        12,
        [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        1.0,
        1.0;
        reaction = Reaction(2, 0),
        dynamics = CorrelatedHoppingDynamics(),
    )
    times, particles = simulate!(
        sys,
        (s, _t) -> is_final_binary(s);
        rng,
        stop_on_reaction_only = true,
    )
    @test first(times) == 0.0
    @test issorted(times)
    @test is_final_binary(sys)
    @test last(particles) == sum(sys.occupations)

    diffusion_sys = initialize_system(
        4,
        [1, 0, 0, 0],
        1.0,
        1.0;
        reaction = Reaction(2, 0),
        dynamics = StandardDiffusion(),
    )
    event_kinds = Symbol[]
    event_count = Ref(0)
    simulate!(
        diffusion_sys,
        (_sys, _t) -> event_count[] >= 3;
        stop_on_reaction_only = false,
        rng = MersenneTwister(3),
        event_callback = (_sys, kind, _site, _time, _change) -> begin
            push!(event_kinds, kind)
            event_count[] += 1
        end,
    )
    @test event_kinds == [:dynamics, :dynamics, :dynamics]
end

@testset "generator" begin
    W, reps, mults = build_generator(4, 1.0, 1.0; max_total = 2)
    @test size(W, 1) == length(reps) == length(mults)
    @test norm(vec(sum(W, dims = 1)), Inf) < 1e-10
end

@testset "ensemble" begin
    rng = MersenneTwister(2)
    results = run_ensemble(
        8,
        2,
        0.0,
        1.0,
        1.0;
        reaction = InstantaneousReaction(2, 0),
        rng,
    )
    @test length(results) == 2
    @test all(results .>= 0)
    @test_throws ArgumentError run_ensemble(8, 0, 0.0, 1.0, 1.0)
    @test_throws ArgumentError run_ensemble(8, 1, 1.5, 1.0, 1.0)
end

@testset "active sites API" begin
    state = [1, 1, 0, 0, 0, 0, 0, 0]
    state_activity = active_sites(state; max_states = 1000)
    sandbox_state_activity = active_sites(state; algorithm = :sandbox, max_states = 1000)
    history = reshape(state, 1, length(state))
    history_activity = active_sites(history; max_states = 1000)

    @test state_activity isa Vector{Bool}
    @test history_activity isa Matrix{Bool}
    @test size(history_activity) == size(history)
    @test sandbox_state_activity == state_activity
    @test history_activity == reshape(state_activity, 1, length(state_activity))
    @test active_sites(state; reaction = Reaction(2, 0), max_states = 1000) isa Vector{Bool}
    @test_throws ArgumentError active_sites(state; algorithm = :unknown)

    non_final_state = [1, 1, 1, 0, 0, 0, 0, 0]
    montecarlo_state_activity = active_sites(
        non_final_state;
        algorithm = :montecarlo,
        reaction = InstantaneousReaction(2, 0),
        n_sims = 4,
        max_steps = 8,
        rng = MersenneTwister(3),
    )
    montecarlo_history_activity = active_sites(
        history;
        algorithm = :montecarlo,
        reaction = InstantaneousReaction(2, 0),
        n_sims = 2,
        max_steps = 4,
        rng = MersenneTwister(4),
    )

    @test montecarlo_state_activity isa Vector{Bool}
    @test montecarlo_history_activity isa Matrix{Bool}
    @test size(montecarlo_history_activity) == size(history)
    @test any(montecarlo_state_activity)
    @test !any(active_sites(
        zeros(Int, length(state));
        algorithm = :montecarlo,
        n_sims = 1,
        rng = MersenneTwister(5),
    ))
    @test !any(active_sites(
        [1, 0, 0, 0, 0, 1, 0, 0];
        algorithm = :montecarlo,
        n_sims = 1,
        final_time_horizon = nothing,
        rng = MersenneTwister(6),
    ))
    final_active_state = [1, 1, 0, 0, 0, 0, 0, 0]
    @test is_final_binary(final_active_state)
    @test !any(active_sites(
        final_active_state;
        algorithm = :montecarlo,
        n_sims = 1,
        max_steps = 1,
        final_time_horizon = nothing,
        rng = MersenneTwister(7),
    ))
    @test any(active_sites(
        final_active_state;
        algorithm = :montecarlo,
        n_sims = 1,
        max_steps = 1,
        rng = MersenneTwister(7),
    ))
    @test_throws ArgumentError active_sites(state; algorithm = :montecarlo, n_sims = 0)
    @test_throws ArgumentError active_sites(state; algorithm = :montecarlo, max_steps = 0)
    @test_throws ArgumentError active_sites(state; algorithm = :montecarlo, final_time_horizon = -1.0)
end

@testset "krylov sectors" begin
    L, N = 6, 3
    counts = symmetry_sector_counts(L, N)
    @test sum(values(counts)) == binomial(big(N + L - 1), L - 1)

    symmetry_sector, symmetry_sector_size, number_of_symmetry_sectors =
        largest_symmetry_sector(L, N)
    states = generate_symmetry_sector_states(
        L,
        N,
        symmetry_sector[1],
        symmetry_sector[2],
        symmetry_sector_size;
        verbose = false,
    )
    @test length(states) == symmetry_sector_size

    krylov_sector_result = analyze_krylov_sectors!(copy(states), L)
    @test krylov_sector_result.certified
    @test krylov_sector_result.largest <= symmetry_sector_size
    @test krylov_sector_result.largest >= krylov_sector_result.remaining

    result = largest_krylov_sector(
        L,
        N;
        target_symmetry_sector = symmetry_sector,
        symmetry_sector_size,
        number_of_symmetry_sectors,
        verbose = false,
    )
    @test result.certified
    @test result.largest == krylov_sector_result.largest
    @test result.ratio == krylov_sector_result.largest / symmetry_sector_size

    @test_throws ArgumentError symmetry_sector_counts(5, N)
end
