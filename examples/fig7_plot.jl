using LaTeXStrings
using Measures
using Plots
using Serialization

###########################################################
# Fig. 7 plotting                                         #
# Load precomputed data and render the figure.            #
###########################################################

function load_fig7_data()
    output_dir = joinpath(@__DIR__, "output")
    path = joinpath(output_dir, "fig7_data.jls")
    isfile(path) || error("Missing fig7 data file: $(path). Run fig7_generate_data.jl first.")
    return deserialize(path)
end

function plot_fig7(data)
    params = data.parameters

    plt = plot(
        xlabel = L"(t-\mu)/\sigma",
        foreground_color_legend = :transparent,
        left_margin = 7.5mm,
        xlims = (-1.2, 2.2),
        ylims = (0, 1.05),
        size = (400, 275),
    )

    green_colors = [
        RGB(0.70, 0.86, 0.64),
        RGB(0.22, 0.60, 0.15),
        RGB(0.04, 0.30, 0.06),
    ]
    for (i, result) in enumerate(data.results)
        L = result.L
        mu = L^2 / (2 * pi^2 * params.gamma)
        sigma = mu
        t_ens = (result.survival_times .- mu) ./ sigma

        clr = green_colors[i]
        stephist!(
            plt,
            t_ens;
            bins = 100,
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
    println("Saved fig7.pdf to ", path)
    return plt
end

data = load_fig7_data()
plot_fig7(data)
