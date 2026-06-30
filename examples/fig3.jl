using CorrelatedHopping
using LaTeXStrings
using Measures
using Plots
using Random

###########################################################
# Fig. 3                                                  #
# Sample trajectories for the particle density rho(t)     #
# using different values of λ/Γ and initial density ρ0.   #
###########################################################

# Parameters
rng = MersenneTwister(2)
L = 2^12
Γ = 1.0
cases = [
    (ρ0 = 0.25, λ = 10.0, label = L"\lambda/\Gamma=10" ),
    (ρ0 = 0.50, λ = 1.00, label = L"\lambda/\Gamma=1"  ),
    (ρ0 = 1.00, λ = 0.10, label = L"\lambda/\Gamma=0.1"),
]
n_trajectories = 10

default_blue = palette(:default)[1]
seafoam = RGB(0.36, 0.75, 0.59)
julia_green = RGB(0.22, 0.60, 0.15)
amber = RGB(0.90, 0.62, 0.00)
lambda_over_gamma_values = [case.λ / Γ for case in cases]
gradient_colors = cgrad([default_blue, seafoam, julia_green, amber, :red])[range(0.0, 1.0, length = length(cases))]
case_colors = fill(default_blue, length(cases))
for (rank, case_index) in enumerate(sortperm(lambda_over_gamma_values))
    case_colors[case_index] = gradient_colors[rank]
end

# Simulate trajectories and compute mean-field curves
trajectories = [
    begin
        initial_state = Int.(rand(rng, L) .< case.ρ0)
        sys = initialize_system(
            L,
            initial_state,
            Γ,
            case.λ;
            reaction = Reaction(2, 0),
            dynamics = CorrelatedHoppingDynamics(),
        )
        times, particle_counts = simulate!(
            sys,
            (sys, _t) -> is_final_binary(sys);
            rng,
        )
        times = Γ .* times .+ 1e-10
        density = particle_counts ./ L
        theory_times = times[times .< 27]
        mean_field_density = 1 ./ (case.λ .* theory_times ./ Γ .+ case.ρ0^(-1))

        (; case_index, ρ0 = case.ρ0, λ = case.λ, times, density, theory_times, mean_field_density)
    end
    for _ in 1:n_trajectories
    for (case_index, case) in enumerate(cases)
]

# Set up plot
plt = plot(
    xlabel = L"t/\Gamma^{-1}",
    foreground_color_legend = :transparent,
    left_margin = 7.5mm,
    label = nothing,
    xscale = :log10,
    xlims = (1e-1, 1e3),
    ylims = (0.1, 1.0),
    yticks = 0:0.2:1.0,
    size = (400, 275),
)
annotate!(plt, 0.03, 0.5, text(L"\rho(t)", 11))

# Plot simulation and theory results
for trajectory in trajectories
    plot!(plt, trajectory.times, trajectory.density; c = case_colors[trajectory.case_index], lw = 2, alpha = 0.3, label = nothing)
    scatter!(
        plt,
        [last(trajectory.times)],
        [last(trajectory.density)];
        c = case_colors[trajectory.case_index],
        ms = 3,
        label = nothing,
    )
end

# Text annotations
annotation_y = [0.25, 0.5, 0.8]
for (i, case) in enumerate(cases)
    annotate!(
        plt,
        0.14,
        annotation_y[i],
        text(case.label, 11, color = case_colors[i], halign = :left, valign = :bottom),
    )
end

# Inset
lens!(
    plt,
    [10, 1000],
    [0.14, 0.20];
    inset = (1, bbox(0.6, 0.0, 0.4, 0.4)),
    alpha = 1.0,
    lw = 2,
)
plot!(
    plt[2];
    xscale = :log10,
    xlabel = nothing,
    ylabel = nothing,
    framestyle = :box,
    xticks = :auto,
    yticks = 0.14:0.03:0.2,
)
x, xcap, y1, y2, xlabel = 75, 65, 0.155, 0.172, 90
plot!(plt[2], [x, x], [y1, y2]; c=:black, lw=1.2, label=nothing)
plot!(plt[2], [xcap, x], [y1, y1]; c=:black, lw=1.2, label=nothing)
plot!(plt[2], [xcap, x], [y2, y2]; c=:black, lw=1.2, label=nothing)
annotate!(plt[2], xlabel, (y1+y2)/2.03, text("Gaussian", 8, halign = :left, valign = :center))
y, ycap, x1, x2, ylabel = 0.15, 0.153, 15, 58, 0.149
plot!(plt[2], [x1, x2], [y, y]; c=:black, lw=1.2, label=nothing)
plot!(plt[2], [x1, x1], [y, ycap]; c=:black, lw=1.2, label=nothing)
plot!(plt[2], [x2, x2], [y, ycap]; c=:black, lw=1.2, label=nothing)
annotate!(plt[2], sqrt(x1*x2), ylabel, text("Gumbel", 8, halign = :center, valign = :top))

# Save figure
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)
path = joinpath(output_dir, "fig3.pdf")
savefig(plt, path)
println("Saved fig3.pdf to ", path)
