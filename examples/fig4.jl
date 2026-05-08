using CorrelatedHopping
using Distributions
using LaTeXStrings
using Measures
using Plots
using Random
using Statistics

###########################################################
# Fig. 4                                                  #
# Convergence of survival times to a Gumbel distribution. #
###########################################################

# Parameters
rng = MersenneTwister(4)
ρ0 = 0.5
λ = 1.0
Γ = 1.0
ensemble_size = 500
L_values = [2^4, 2^8, 2^12]

# Run the ensemble simulations
results = [
    begin
        survival_times = run_ensemble(
            L,
            ensemble_size,
            ρ0,
            Γ,
            λ;
            reaction = Reaction(2, 0),
            dynamics = CorrelatedHoppingDynamics(),
            rng,
        )

        (; L, survival_times)
    end
    for L in L_values
]

# Set up plot
plt = plot(
    xlabel = L"t",
    foreground_color_legend = :transparent,
    left_margin = 7.5mm,
    xlims = (0, 100),
    ylims = (0, 0.11),
    size = (400, 275),
    yticks = 0:0.03:0.12,
)

# Plot histograms of the survival times and Gumbel distributions from empirical mean and standard deviation.
for (i, result) in enumerate(results)
    survival_times = result.survival_times

    theta_est = std(survival_times) * sqrt(6) / pi
    mu_est = mean(survival_times) - theta_est * 0.5772156649015329
    fitted_gumbel = Gumbel(mu_est, theta_est)
    x_vals = range(0.001, maximum(survival_times); length = 500)
    gumbel_vals = pdf.(fitted_gumbel, x_vals)

    clr = palette(:lisbon10)[i + 1]
    stephist!(
        plt,
        survival_times;
        bins = 100,
        normalize = :pdf,
        lw = 2,
        c = clr,
        label = nothing,
    )
    plot!(plt, x_vals, gumbel_vals; lw = 2, c = clr, label = "L=$(result.L)")
end
annotate!(plt, -15, 0.045, text(L"P(t)", 11))

# Save figure
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)
path = joinpath(output_dir, "fig4.pdf")
savefig(plt, path)
println("Saved fig4.pdf to ", path)
