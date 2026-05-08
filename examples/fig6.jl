using CorrelatedHopping
using LaTeXStrings
using Measures
using Plots
using Random
using Statistics

###########################################################
# Fig. 6                                                  #
# Mean survival times across density and rate ratios.     #
###########################################################

"""
Helper function which obtains an ensemble of escape times 
for a grid of density and rate ratios.
"""
function simulate_rate_grid(L, ensemble_size, ρ0_values, λ_over_Γ_values, Γ, rng)
    μT = zeros(length(λ_over_Γ_values), length(ρ0_values))
    σT = zeros(length(λ_over_Γ_values), length(ρ0_values))
    λ_values = zeros(length(λ_over_Γ_values))

    for (i, λ_over_Γ) in enumerate(λ_over_Γ_values)
        λ = λ_over_Γ * Γ
        λ_values[i] = λ
        for (j, ρ0) in enumerate(ρ0_values)
            samples = run_ensemble(
                L,
                ensemble_size,
                ρ0,
                Γ,
                λ;
                reaction = Reaction(2, 0),
                dynamics = CorrelatedHoppingDynamics(),
                rng,
            )
            μT[i, j] = mean(samples)
            σT[i, j] = std(samples) / sqrt(ensemble_size)
        end
    end

    return (; μT, σT, λ_values)
end

"""
Helper function which obtains the amplitude c1 of the slowest decay mode
as a function of density ρ0 in the λ → ∞ limit.
"""
function c1_curve()
    L = 12
    Γ = 1.0
    λ = 10.0 # A close enough approximation of the λ → ∞ limit for L = 12 # Eventually update so that build_generator works in the direct limit
    ρ0_range = LinRange(0, 1, 100)
    W, representatives, multiplicities = build_generator(L, Γ, λ; lump_final = true)
    _, right_vectors, left_vectors = find_spectral_gap(W; num_eigenvalues = 2)
    c1_values = [
        calculate_coefficient(
            ρ0,
            1,
            right_vectors,
            left_vectors,
            representatives,
            multiplicities,
            L,
        ) for ρ0 in ρ0_range
    ]
    return (; ρ0_range, c1_values)
end

# Parameters
rng = MersenneTwister(6)
L = 2^9
Γ = 1.0
ensemble_size = 10
ρ0_values = LinRange(0, 1, 20)
λ_over_Γ_powers = -3:3
λ_over_Γ_values = 10.0 .^ λ_over_Γ_powers
n_display_ratios = 4
slow_reaction_indices = firstindex(λ_over_Γ_values):(firstindex(λ_over_Γ_values) + n_display_ratios - 1)
fast_reaction_indices = (lastindex(λ_over_Γ_values) - n_display_ratios + 1):lastindex(λ_over_Γ_values)

# Calculate average escape time for the parameter range
simulation = simulate_rate_grid(L, ensemble_size, ρ0_values, λ_over_Γ_values, Γ, rng)
# Calculate c1 curve for the inset plot
inset = c1_curve()

# Set up left plot: fast reactions
plt1 = plot(
    xlabel = L"\rho_0",
    left_margin = 8.5mm,
    bottom_margin = 4mm,
    right_margin = 15mm,
    xlims = (0, 1.02),
    ylims = (0, 53),
    legend = (1.17, 0.66),
    foreground_color_legend = :transparent,
    size = (400, 275),
)
annotate!(plt1, -0.17, 25, text(L"\frac{\langle \hat{T}_L\rangle}{\Gamma^{-1}}", 11))
annotate!(plt1, -0.17, 50, text(L"(a)", 11))

# Plot simulation
colors = cgrad([:blue, :red])[range(0.0, 1.0, length = length(λ_over_Γ_values))]
for i in fast_reaction_indices
    scatter!(
        plt1,
        ρ0_values,
        Γ .* simulation.μT[i, :];
        yerror = Γ .* simulation.σT[i, :],
        c = colors[i],
        label = nothing,
    )
end

# Plot inset
plot!(
    plt1,
    inset.ρ0_range,
    inset.c1_values;
    inset = (1, bbox(0.465, 0.66, 0.3, 0.2)),
    subplot = 2,
    label = nothing,
    lw = 2,
    ylims = (0, 0.10),
    yticks = 0:0.05:0.1,
    xlims = (0, 1),
    xticks = 0:0.5:1.0,
    framestyle = :box,
    guidefontsize = 7,
    tickfontsize = 6,
    c = :red,
    xtickfontvalign = :bottom,
    ytickfonthalign = :left,
)
annotate!(plt1, 0.37, 12.5, text(L"c_1", 9))
annotate!(plt1, 0.64, 3, text(L"\rho_0", 9))

# Global legend
labels = [L"\lambda/\Gamma = 10^{%$p}" for p in λ_over_Γ_powers]
for i in reverse(eachindex(labels))
    scatter!(plt1, [NaN], [NaN]; c = colors[i], label = labels[i])
end
annotate!(plt1, 1.25, 42, text("Fast\nreactions", 8, color = :red))
annotate!(plt1, 1.25, 6, text("Slow\nreactions", 8, color = :blue))

# Set up right plot: slow reactions
plt2 = plot(
    xlabel = L"\rho_0",
    left_margin = 23.5mm,
    bottom_margin = 4mm,
    xlims = (0, 1.02),
    ylims = (0, 53),
    legend = (0.6, 0.27),
    foreground_color_legend = :transparent,
    size = (400, 275),
)
annotate!(plt2, -0.17, 25, text(L"\frac{\langle \hat{T}_L\rangle}{\lambda^{-1}}", 11))
annotate!(plt2, -0.17, 50, text(L"(b)", 11))

# Plot simulation and mean-field curves
ρ0_mf = LinRange(0.15, 1.0, 200)
for i in slow_reaction_indices
    μT_λ_units = simulation.λ_values[i] .* simulation.μT[i, :]
    σT_λ_units = simulation.λ_values[i] .* simulation.σT[i, :]

    plot!(
        plt2,
        ρ0_mf,
        (last(μT_λ_units) + 1) .- 1 ./ ρ0_mf; # Mean-field curve
        c = colors[i],
        lw = 2,
        ls = :dash,
        label = nothing,
    )
    scatter!(
        plt2,
        ρ0_values,
        μT_λ_units;
        yerror = σT_λ_units,
        c = colors[i],
        label = nothing,
    )
end

# Save figure
plt = plot(plt1, plt2; layout = (1, 2), size = (860, 275))
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)
path = joinpath(output_dir, "fig6.pdf")
savefig(plt, path)
println("Saved fig6.pdf to ", path)
