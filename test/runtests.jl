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
