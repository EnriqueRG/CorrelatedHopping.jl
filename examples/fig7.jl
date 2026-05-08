using CorrelatedHopping
using LaTeXStrings
using Measures
using Plots
using Random

# Port of the PaperPlots diffusive final-time figure.

rng = MersenneTwister(7)
rho0 = 0.5
gamma = 1.0
lambda = 1.0
n_ens = 200
L_values = [2^4, 2^8, 2^12]

simulation_data = Dict{Int,Vector{Float64}}()
for L in L_values
    simulation_data[L] = run_ensemble(
        L,
        n_ens,
        rho0,
        gamma,
        lambda;
        reaction = Reaction(2, 0),
        dynamics = StandardDiffusion(),
        condition_even = true,
        stop = (s, _t) -> sum(s.occupations) < 2,
        stop_on_reaction_only = true,
        rng,
    )
end

plt = plot(
    xlabel = L"(t-\mu)/\sigma",
    foreground_color_legend = :transparent,
    left_margin = 7.5mm,
    xlims = (-1.2, 2.2),
    ylims = (0, 1.05),
    size = (400, 275),
)

for (i, L) in enumerate(L_values)
    mu = L^2 / (2 * pi^2 * gamma)
    sigma = mu
    t_ens = (simulation_data[L] .- mu) ./ sigma

    clr = palette(:lisbon10)[i + 1]
    stephist!(
        plt,
        t_ens;
        bins = 80,
        normalize = :pdf,
        lw = 2,
        c = clr,
        label = "L=$L",
    )
end

x_range = range(-1.0, 2.2; length = 500)
y_vals = exp.(-x_range .- 1)
plot!(
    plt,
    x_range,
    y_vals;
    lw = 2,
    ls = :dash,
    c = :black,
    label = "Exponential distribution",
)

annotate!(plt, -1.75, 0.5, text(L"P(t)", 11))

output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)
path = joinpath(output_dir, "fig7.pdf")
savefig(plt, path)
println("saved ", path)
