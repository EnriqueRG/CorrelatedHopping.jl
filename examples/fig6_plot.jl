using LaTeXStrings
using Measures
using Plots
using Serialization

###########################################################
# Fig. 6 plotting                                         #
# Load precomputed data and render the figure.            #
###########################################################

function load_fig6_data()
    output_dir = joinpath(@__DIR__, "output")
    path = joinpath(output_dir, "fig6_data.jls")
    isfile(path) || error("Missing fig6 data file: $(path). Run fig6_generate_data.jl first.")
    return deserialize(path)
end

function interpolate_at_density(rho0_values, values, rho_ref)
    exact_index = findfirst(rho0 -> isapprox(rho0, rho_ref; atol = 1e-12, rtol = 0), rho0_values)
    if !isnothing(exact_index)
        return values[exact_index]
    end

    upper_index = searchsortedfirst(rho0_values, rho_ref)
    if upper_index <= firstindex(rho0_values) || upper_index > lastindex(rho0_values)
        error("Reference density $(rho_ref) is outside the simulated density range.")
    end

    lower_index = upper_index - 1
    lower_rho = rho0_values[lower_index]
    upper_rho = rho0_values[upper_index]
    weight = (rho_ref - lower_rho) / (upper_rho - lower_rho)
    return (1 - weight) * values[lower_index] + weight * values[upper_index]
end

function plot_fig6(data)
    params = data.parameters
    simulation = data.simulation
    inset = data.inset

    rho0_values = params.rho0_values
    hop_rate = params.hop_rate
    reaction_over_hop_powers = params.reaction_over_hop_powers
    reaction_over_hop_values = params.reaction_over_hop_values
    n_display_ratios = params.n_display_ratios
    slow_reaction_indices = firstindex(reaction_over_hop_values):(firstindex(reaction_over_hop_values) + n_display_ratios - 1)
    fast_reaction_indices = (lastindex(reaction_over_hop_values) - n_display_ratios + 1):lastindex(reaction_over_hop_values)

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
    default_blue = palette(:default)[1]
    seafoam = RGB(0.36, 0.75, 0.59)
    julia_green = RGB(0.22, 0.60, 0.15)
    amber = RGB(0.90, 0.62, 0.00)
    colors = cgrad([default_blue, seafoam, julia_green, amber, :red])[range(0.0, 1.0, length = length(reaction_over_hop_values))]
    for i in fast_reaction_indices
        scatter!(
            plt1,
            rho0_values,
            hop_rate .* simulation.mean_times[i, :];
            yerror = hop_rate .* simulation.stderr_times[i, :],
            c = colors[i],
            label = nothing,
        )
    end

    # Plot inset
    plot!(
        plt1,
        inset.rho0_values,
        inset.c1_values;
        inset = (1, bbox(0.465, 0.66, 0.3, 0.2)),
        subplot = 2,
        label = nothing,
        lw = 2,
        ylims = (0, 0.10),
        yticks = 0:0.05:0.1,
        xlims = (0, 1),
        xticks = 0.5:0.5:1.0,
        framestyle = :box,
        guidefontsize = 7,
        tickfontsize = 6,
        c = :red,
        xtickfontvalign = :bottom,
        ytickfonthalign = :left,
    )
    annotate!(plt1, 0.36, 12.5, text(L"c_1", 9))
    annotate!(plt1, 0.64, 2, text(L"\rho_0", 9))

    # Global legend
    labels = [L"\lambda/\Gamma = 10^{%$p}" for p in reaction_over_hop_powers]
    for i in reverse(eachindex(labels))
        scatter!(plt1, [NaN], [NaN]; c = colors[i], label = labels[i])
    end
    annotate!(plt1, 1.25, 42, text("Fast\nreactions", 8, color = :red))
    annotate!(plt1, 1.25, 6, text("Slow\nreactions", 8, color = default_blue))

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
    rho0_mf = LinRange(0.15, 1.0, 200)
    rho_ref = 1 / 2
    for i in reverse(slow_reaction_indices)
        mean_times_reaction_units = simulation.reaction_rate_values[i] .* simulation.mean_times[i, :]
        stderr_times_reaction_units = simulation.reaction_rate_values[i] .* simulation.stderr_times[i, :]
        reference_mean_time = interpolate_at_density(rho0_values, mean_times_reaction_units, rho_ref)

        plot!(
            plt2,
            rho0_mf,
            reference_mean_time .+ 1 / rho_ref .- 1 ./ rho0_mf; # Mean-field curve
            c = colors[i],
            lw = 2,
            ls = :dash,
            label = nothing,
        )
        scatter!(
            plt2,
            rho0_values,
            mean_times_reaction_units;
            yerror = stderr_times_reaction_units,
            c = colors[i],
            label = nothing,
        )
    end

    plt = plot(plt1, plt2; layout = (1, 2), size = (860, 275))
    output_dir = joinpath(@__DIR__, "output")
    mkpath(output_dir)
    path = joinpath(output_dir, "fig6.pdf")
    savefig(plt, path)
    println("Saved fig6.pdf to ", path)
    return plt
end


data = load_fig6_data()
plot_fig6(data)
