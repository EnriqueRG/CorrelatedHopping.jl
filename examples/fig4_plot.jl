using Distributions
using LaTeXStrings
using Measures
using Plots
using Serialization
using Statistics

###########################################################
# Fig. 4 plotting                                         #
# Load precomputed data and render the figure.            #
###########################################################

function load_fig4_data()
    output_dir = joinpath(@__DIR__, "output")
    path = joinpath(output_dir, "fig4_data.jls")
    isfile(path) || error("Missing fig4 data file: $(path). Run fig4_generate_data.jl first.")
    return deserialize(path)
end

function plot_fig4(data)
    plt = plot(
        xlabel = L"t",
        foreground_color_legend = :transparent,
        left_margin = 7.5mm,
        xlims = (0, 100),
        ylims = (0, 0.11),
        size = (400, 275),
        yticks = 0:0.03:0.12,
    )

    green_colors = [
        RGB(0.70, 0.86, 0.64),
        RGB(0.22, 0.60, 0.15),
        RGB(0.04, 0.30, 0.06),
    ]
    for (i, result) in enumerate(data.results)
        final_times = result.final_times

        theta_est = std(final_times) * sqrt(6) / pi
        mu_est = mean(final_times) - theta_est * 0.5772156649015329
        fitted_gumbel = Gumbel(mu_est, theta_est)
        x_vals = range(0.001, maximum(final_times); length = 500)
        gumbel_vals = pdf.(fitted_gumbel, x_vals)

        clr = green_colors[i]
        stephist!(
            plt,
            final_times;
            bins = 100,
            normalize = :pdf,
            lw = 2,
            c = clr,
            label = nothing,
        )
        plot!(plt, x_vals, gumbel_vals; lw = 2, c = clr, label = "L=$(result.L)")
    end
    annotate!(plt, -15, 0.045, text(L"P(t)", 11))

    output_dir = joinpath(@__DIR__, "output")
    mkpath(output_dir)
    path = joinpath(output_dir, "fig4.pdf")
    savefig(plt, path)
    println("Saved fig4.pdf to ", path)
    return plt
end

data = load_fig4_data()
plot_fig4(data)
