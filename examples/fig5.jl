using CorrelatedHopping
using LaTeXStrings
using LinearAlgebra
using Measures
using Plots
using Random
using Roots
using Statistics

###########################################################
# Fig. 5                                                  #
# Extreme-value scaling of survival times.                #
###########################################################

"""
Computes the decay rates γk and amplitudes ck of the cdf tail
for the reduced ℓ=12 model, which captures the evolution of a small number of states.
"""
function reducedℓ12_ckγk(λ, γ, ρ0)
    function slow_mode(M, p0)
        F = eigen(M)
        order = sortperm(real.(F.values), rev = true)
        k = order[2]              # Index of the slowest non-zero mode
        v = F.vectors[:, k]       # Right eigenvector
        u = inv(F.vectors)[k, :]  # Left eigenvector
        initialStateOverlap = sum(u .* p0) / sum(u .* v)
        return -real(F.values[k]), -real(initialStateOverlap * v[1])
    end

    # Slowest mode
    M1 = [
        0.0 λ 0.0 0.0 0.0 0.0
        0.0 -λ-2γ γ 0.0 0.0 0.0
        0.0 2γ -2γ γ 0.0 0.0
        0.0 0.0 γ -2γ γ 0.0
        0.0 0.0 0.0 γ -2γ γ
        0.0 0.0 0.0 0.0 γ -γ
    ]
    w1 = 26 * ρ0^4 * (1 - ρ0)^(12 - 4)
    γ1, c1 = slow_mode(M1, [0.0, 0.0, w1, w1, w1, w1])

    # Second slowest mode
    M2 = [
        0.0 λ 0.0 0.0 0.0
        0.0 -λ-2γ γ 0.0 0.0
        0.0 2γ -2γ γ 0.0
        0.0 0.0 γ -2γ γ
        0.0 0.0 0.0 γ -γ
    ]
    w2 = 26 * ρ0^3 * (1 - ρ0)^(12 - 3)
    γ2, c2 = slow_mode(M2, [0.0, 0.0, w2, w2, w2])

    return [γ1, γ2], [c1, c2]
end


"""
Computes the extreme-value scaling constants (an, bn) 
given the decay rates γk and amplitudes ck of the cdf tail.
"""
function anbn_from_ckγk(n, ck, γk)
    cdf_tail(t) = sum(c * exp(-γ * t) for (c, γ) in zip(ck, γk))
    pdf_tail(t) = sum(c * γ * exp(-γ * t) for (c, γ) in zip(ck, γk))

    if cdf_tail(0.0) <= 1 / n
        # No solution is possible
        # This situation can arise if n is not large enough, so we are not in the asymptotic regime where the theory applies
        return (NaN, NaN)
    else
        # When the solution exists, it is unique and guaranteed to be in this interval
        bracket = (0.0, log(n * sum(ck)) / first(γk))

        bn = find_zero(t -> cdf_tail(t) - 1 / n, bracket)
        an = 1 / (n * pdf_tail(bn))
        return an, bn
    end
end


"""
Computes the expected value of T_L as the maximum of n i.i.d. random variables 
drawn from a distribution with tail decay rates γk and amplitudes ck.
"""
function meanT_evs(n_range, ck, γk)
    return map(n_range) do n
        an, bn = anbn_from_ckγk(n, ck, γk)
        0.5772156649 * an + bn
    end
end

# Parameters
rng = MersenneTwister(5)
ρ0_values = [0.25, 0.75]
λ = 1.0
Γ = 1.0
ℓ = 12
ensemble_size = 1000
L_values = ℓ * round.(Int, 2 .^ range(0, 10, 10))
n_replicas = 2 .^ range(3, 10, 100)

# Run simulations and compute theory values
results = map(ρ0_values) do ρ0
    samples = [
        run_ensemble(
            L,
            ensemble_size,
            ρ0,
            Γ,
            λ;
            reaction = Reaction(2, 0),
            dynamics = CorrelatedHoppingDynamics(),
            rng,
        )
        for L in L_values
    ]

    μT_simulation = mean.(samples)
    σT_simulation = std.(samples) ./ sqrt(ensemble_size)

    full_gammas, full_coefficients = full_gamma_coefficients(ℓ, Γ, λ, ρ0; nmax = 2)
    μT_theory_full = meanT_evs(n_replicas, full_coefficients, full_gammas)

    reduced_gammas, reduced_coefficients = reducedℓ12_ckγk(λ, Γ, ρ0)
    μT_theory_rℓ12 = meanT_evs(n_replicas, reduced_coefficients, reduced_gammas)

    (; ρ0, μT_simulation, σT_simulation, μT_theory_full, μT_theory_rℓ12)
end

# Set up plot
max_μT = maximum(maximum(result.μT_simulation) for result in results)
ytick_step = 5 * floor(Int, max_μT / (3 * 5))
yticks = ytick_step .* (0:3)
ymax = 1.1 * max_μT
plt = plot(
    xlabel = L"L",
    left_margin = 7.5mm,
    ylims = (0, ymax),
    yticks = yticks,
    xscale = :log10,
    legend = (0.2, 0.9),
    foreground_color_legend = :transparent,
    size = (400, 275),
)
annotate!(plt, 4.0, 0.5*max_μT, text(L"\frac{\langle \hat{T}_L\rangle}{\Gamma^{-1}}", 11))

# Plot simulation and theory results
for (i, result) in enumerate(results)
    scatter!(
        plt,
        L_values,
        result.μT_simulation;
        yerror = result.σT_simulation,
        c = i,
        label = "ρ₀=$(result.ρ0)",
    )
    plot!(plt, ℓ .* n_replicas, result.μT_theory_full; ls = :dash, lw = 2, c = i, label = nothing)
    plot!(plt, ℓ .* n_replicas, result.μT_theory_rℓ12; ls = :dot, lw = 2, c = i, label = nothing)
end

# Dummy lines for legend
plot!(plt, [NaN], [NaN]; ls = :dash, lw = 2, c = :black, label = "Numerical ℓ=$ℓ")
plot!(plt, [NaN], [NaN]; ls = :dot, lw = 2, c = :black, label = "Reduced ℓ=$ℓ")

# Save plot
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)
path = joinpath(output_dir, "fig5.pdf")
savefig(plt, path)
println("Saved fig5.pdf to ", path)
