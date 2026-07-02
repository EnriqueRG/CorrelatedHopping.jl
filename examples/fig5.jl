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
# Extreme-value scaling of final times.                   #
###########################################################

"""
Computes the decay rates and tail amplitudes for the reduced ell=12 model,
which captures the evolution of a small number of states.
"""
function reduced_ell12_decay_tail(reaction_rate, hop_rate, rho0)
    function slow_mode(M, p0)
        F = eigen(M)
        order = sortperm(real.(F.values), rev = true)
        k = order[2]              # Index of the slowest non-zero mode
        v = F.vectors[:, k]       # Right eigenvector
        u = inv(F.vectors)[k, :]  # Left eigenvector
        initial_state_overlap = sum(u .* p0) / sum(u .* v)
        return -real(F.values[k]), -real(initial_state_overlap * v[1])
    end

    # Slowest mode
    M1 = [
        0.0 reaction_rate 0.0 0.0 0.0 0.0
        0.0 -(reaction_rate + 2 * hop_rate) hop_rate 0.0 0.0 0.0
        0.0 (2 * hop_rate) -(2 * hop_rate) hop_rate 0.0 0.0
        0.0 0.0 hop_rate -(2 * hop_rate) hop_rate 0.0
        0.0 0.0 0.0 hop_rate -(2 * hop_rate) hop_rate
        0.0 0.0 0.0 0.0 hop_rate -hop_rate
    ]
    w1 = 26 * rho0^4 * (1 - rho0)^(12 - 4)
    decay_rate_1, coefficient_1 = slow_mode(M1, [0.0, 0.0, w1, w1, w1, w1])

    # Second slowest mode
    M2 = [
        0.0 reaction_rate 0.0 0.0 0.0
        0.0 -(reaction_rate + 2 * hop_rate) hop_rate 0.0 0.0
        0.0 (2 * hop_rate) -(2 * hop_rate) hop_rate 0.0
        0.0 0.0 hop_rate -(2 * hop_rate) hop_rate
        0.0 0.0 0.0 hop_rate -hop_rate
    ]
    w2 = 26 * rho0^3 * (1 - rho0)^(12 - 3)
    decay_rate_2, coefficient_2 = slow_mode(M2, [0.0, 0.0, w2, w2, w2])

    return [decay_rate_1, decay_rate_2], [coefficient_1, coefficient_2]
end


"""
Computes the extreme-value scaling constants (an, bn) 
given the decay rates and tail amplitudes of the cdf tail.
"""
function extreme_value_constants(n, tail_coefficients, decay_rates)
    cdf_tail(t) = sum(c * exp(-decay_rate * t) for (c, decay_rate) in zip(tail_coefficients, decay_rates))
    pdf_tail(t) = sum(c * decay_rate * exp(-decay_rate * t) for (c, decay_rate) in zip(tail_coefficients, decay_rates))

    if cdf_tail(0.0) <= 1 / n
        # No solution is possible
        # This situation can arise if n is not large enough, so we are not in the asymptotic regime where the theory applies
        return (NaN, NaN)
    else
        # When the solution exists, it is unique and guaranteed to be in this interval
        bracket = (0.0, log(n * sum(tail_coefficients)) / first(decay_rates))

        bn = find_zero(t -> cdf_tail(t) - 1 / n, bracket)
        an = 1 / (n * pdf_tail(bn))
        return an, bn
    end
end


"""
Computes the expected value of T_L as the maximum of n i.i.d. random variables 
drawn from a distribution with given tail decay rates and amplitudes.
"""
function mean_extreme_value_times(n_range, tail_coefficients, decay_rates)
    return map(n_range) do n
        an, bn = extreme_value_constants(n, tail_coefficients, decay_rates)
        0.5772156649 * an + bn
    end
end

# Parameters
rng = MersenneTwister(5)
rho0_values = [0.25, 0.75]
reaction_rate = 1.0
hop_rate = 1.0
ell = 12
ensemble_size = 1000
L_values = ell * round.(Int, 2 .^ range(0, 10, 10))
n_replicas = 2 .^ range(3, 10, 100)

# Run simulations and compute theory values
results = map(rho0_values) do rho0
    samples = [
        run_ensemble(
            L,
            ensemble_size,
            rho0;
            hop_rate,
            reaction_rate,
            reaction = Reaction(2, 0),
            dynamics = CorrelatedHoppingDynamics(),
            rng,
        )
        for L in L_values
    ]

    mean_times_simulation = mean.(samples)
    stderr_times_simulation = std.(samples) ./ sqrt(ensemble_size)

    full_decay_rates, full_coefficients = full_decay_coefficients(
        ell,
        rho0;
        hop_rate,
        reaction_rate,
        nmax = 2,
    )
    mean_times_theory_full = mean_extreme_value_times(n_replicas, full_coefficients, full_decay_rates)

    reduced_decay_rates, reduced_coefficients = reduced_ell12_decay_tail(reaction_rate, hop_rate, rho0)
    mean_times_theory_reduced_ell = mean_extreme_value_times(n_replicas, reduced_coefficients, reduced_decay_rates)

    (; rho0, mean_times_simulation, stderr_times_simulation, mean_times_theory_full, mean_times_theory_reduced_ell)
end

# Set up plot
max_mean_time = maximum(maximum(result.mean_times_simulation) for result in results)
ytick_step = 5 * floor(Int, max_mean_time / (3 * 5))
yticks = ytick_step .* (0:3)
ymax = 1.1 * max_mean_time
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
annotate!(plt, 4.0, 0.5 * max_mean_time, text(L"\frac{\langle \hat{T}_L\rangle}{\Gamma^{-1}}", 11))

# Plot simulation and theory results
for (i, result) in enumerate(results)
    scatter!(
        plt,
        L_values,
        result.mean_times_simulation;
        yerror = result.stderr_times_simulation,
        c = i,
        label = L"\rho_0=%$(result.rho0)",
    )
    plot!(plt, ell .* n_replicas, result.mean_times_theory_full; ls = :dash, lw = 2, c = i, label = nothing)
    plot!(plt, ell .* n_replicas, result.mean_times_theory_reduced_ell; ls = :dot, lw = 2, c = i, label = nothing)
end

# Dummy lines for legend
plot!(plt, [NaN], [NaN]; ls = :dash, lw = 2, c = :black, label = L"Numerical \ell=%$ell")
plot!(plt, [NaN], [NaN]; ls = :dot, lw = 2, c = :black, label = L"Reduced \ell=%$ell")

# Save plot
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)
path = joinpath(output_dir, "fig5.pdf")
savefig(plt, path)
println("Saved fig5.pdf to ", path)
