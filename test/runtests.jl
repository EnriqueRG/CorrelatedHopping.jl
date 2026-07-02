using CorrelatedHopping
using Random
using Test

@testset "public exports" begin
    exported_names = names(CorrelatedHopping)
    @test :simulate! in exported_names
    @test :run_ensemble in exported_names
    @test :build_pair_annihilation_generator in exported_names
    @test :update_local_rates! ∉ exported_names
    @test :calculate_coefficient ∉ exported_names
    @test :largest_krylov_sector ∉ exported_names
end

@testset "reactions and final states" begin
    @test Reaction(2, 0) isa Reaction{2,0}
    @test Reaction(1, 2) isa Reaction{1,2}
    @test InstantaneousReaction(2, 0) isa InstantaneousReaction{2,0}
    @test_throws ArgumentError Reaction(0, 0)
    @test_throws ArgumentError Reaction(1, -1)
    @test_throws ArgumentError Reaction{0,0}()
    @test_throws ArgumentError Reaction{1,-1}()
    @test_throws ArgumentError InstantaneousReaction(0, 0)
    @test_throws ArgumentError InstantaneousReaction(2, -1)
    @test_throws ArgumentError InstantaneousReaction(2, 2)
    @test_throws ArgumentError InstantaneousReaction(2, 3)
    @test_throws ArgumentError InstantaneousReaction{2,2}()
    @test_throws ArgumentError InstantaneousReaction{2,3}()

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
    result = simulate!(
        sys,
        (s, _t) -> is_final_binary(s);
        rng,
        stop_on_reaction_only = true,
    )
    @test first(result.recorded_times) == 0.0
    @test issorted(result.recorded_times)
    @test is_final_binary(sys)
    @test result.final_time >= last(result.recorded_times)
    @test result.final_state == sys.occupations
    @test last(result.history) == sum(sys.occupations)

    clock_sys = initialize_system(
        4,
        [1, 0, 0, 0],
        1.0,
        0.0;
        reaction = Reaction(2, 0),
        dynamics = StandardDiffusion(),
    )
    clock_result = simulate!(
        clock_sys,
        (_sys, t) -> t > 0;
        stop_on_reaction_only = false,
        rng = MersenneTwister(2),
    )
    @test clock_result.final_time > last(clock_result.recorded_times)
    @test clock_result.recorded_times == [0.0]
    @test clock_result.history == [1]

    full_record_sys = initialize_system(
        4,
        [1, 0, 0, 0],
        1.0,
        0.0;
        reaction = Reaction(2, 0),
        dynamics = StandardDiffusion(),
    )
    full_record_result = simulate!(
        full_record_sys,
        (_sys, t) -> t > 0;
        record = :all,
        stop_on_reaction_only = false,
        rng = MersenneTwister(2),
    )
    @test length(full_record_result.recorded_times) == size(full_record_result.history, 1)
    @test size(full_record_result.history, 2) == 4

    no_record_sys = initialize_system(
        4,
        [1, 0, 0, 0],
        1.0,
        0.0;
        reaction = Reaction(2, 0),
        dynamics = StandardDiffusion(),
    )
    no_record_result = simulate!(
        no_record_sys,
        (_sys, t) -> t > 0;
        record = :none,
        stop_on_reaction_only = false,
        rng = MersenneTwister(2),
    )
    @test no_record_result.final_time > 0.0
    @test isempty(no_record_result.recorded_times)
    @test no_record_result.history === nothing
    @test_throws ArgumentError simulate!(no_record_sys, (_sys, _t) -> true; record = :unknown)

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

@testset "pair annihilation generator" begin
    W, reps, mults = build_pair_annihilation_generator(
        4;
        hop_rate = 1.0,
        reaction_rate = 1.0,
        max_total = 2,
    )
    @test size(W, 1) == length(reps) == length(mults)
    @test maximum(abs, vec(sum(W, dims = 1))) < 1e-10

    W_rates, reps_rates, _ = build_pair_annihilation_generator(
        5;
        hop_rate = 2.0,
        reaction_rate = 5.0,
        max_total = 2,
    )
    rep_index(state) = findfirst(==(CorrelatedHopping._representative(state, false)), reps_rates)

    empty_idx = rep_index(zeros(Int, 5))
    reaction_source_idx = rep_index([2, 0, 0, 0, 0])
    @test W_rates[empty_idx, reaction_source_idx] ≈ 5.0

    hop_source_idx = rep_index([0, 1, 1, 0, 0])
    hop_dest_idx = rep_index([1, 0, 0, 1, 0])
    @test W_rates[hop_dest_idx, hop_source_idx] ≈ 2.0

    @test_throws ArgumentError build_pair_annihilation_generator(4; hop_rate = -1.0)
    @test_throws ArgumentError build_pair_annihilation_generator(4; reaction_rate = -1.0)
end

@testset "ensemble" begin
    rng = MersenneTwister(2)
    results = run_ensemble(
        8,
        2,
        0.0;
        hop_rate = 1.0,
        reaction_rate = 1.0,
        reaction = InstantaneousReaction(2, 0),
        rng,
    )
    @test length(results) == 2
    @test all(results .>= 0)
    @test_throws ArgumentError run_ensemble(8, 0, 0.0; hop_rate = 1.0, reaction_rate = 1.0)
    @test_throws ArgumentError run_ensemble(8, 1, 1.5; hop_rate = 1.0, reaction_rate = 1.0)
    @test_throws ArgumentError run_ensemble(8, 1, 0.0; hop_rate = -1.0, reaction_rate = 1.0)
    @test_throws ArgumentError run_ensemble(8, 1, 0.0; hop_rate = 1.0, reaction_rate = -1.0)
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

    search_rng = MersenneTwister(11)
    for L in 4:8, _ in 1:12
        random_state = Int.(rand(search_rng, L) .< 0.45)
        full_activity = CorrelatedHopping._search_active_sites(
            random_state;
            reaction = InstantaneousReaction(2, 0),
            max_states = 10_000,
        )
        sandbox_activity = CorrelatedHopping._sandbox_search(
            random_state;
            reaction = InstantaneousReaction(2, 0),
            max_states = 10_000,
        )
        @test sandbox_activity == full_activity
    end

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
    counts = CorrelatedHopping.symmetry_sector_counts(L, N)
    @test sum(values(counts)) == binomial(big(N + L - 1), L - 1)

    symmetry_sector, symmetry_sector_size, number_of_symmetry_sectors =
        CorrelatedHopping.largest_symmetry_sector(L, N)
    states = CorrelatedHopping.generate_symmetry_sector_states(
        L,
        N,
        symmetry_sector[1],
        symmetry_sector[2],
        symmetry_sector_size;
        verbose = false,
    )
    @test length(states) == symmetry_sector_size

    krylov_sector_result = CorrelatedHopping.analyze_krylov_sectors!(copy(states), L)
    @test krylov_sector_result.certified
    @test krylov_sector_result.largest <= symmetry_sector_size
    @test krylov_sector_result.largest >= krylov_sector_result.remaining

    result = CorrelatedHopping.largest_krylov_sector(
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

    @test_throws ArgumentError CorrelatedHopping.symmetry_sector_counts(5, N)
end
