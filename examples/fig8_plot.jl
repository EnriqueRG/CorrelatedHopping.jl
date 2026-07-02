using LaTeXStrings
using Measures
using Plots

###########################################################
# Fig. 8 plotting                                         #
# Load certified Krylov-sector data and render.           #
###########################################################

output_dir = joinpath(@__DIR__, "output")
data_path = joinpath(output_dir, "fig8_data.csv")
pdf_path = joinpath(output_dir, "fig8.pdf")

function read_fig8_data(path)
    isfile(path) || error("Missing fig8 data file: $(path). Run fig8_generate_data.jl first.")
    rows = NamedTuple[]
    for (line_number, line) in enumerate(eachline(path))
        line_number == 1 && continue
        isempty(strip(line)) && continue
        fields = split(line, ','; limit = 6)
        length(fields) == 6 || error("Expected 6 CSV fields on line $line_number, got $(length(fields))")
        push!(
            rows,
            (
                series = fields[1],
                L = parse(Int, fields[2]),
                ratio = parse(Float64, fields[5]),
            ),
        )
    end
    return rows
end

function plot_fig8(rows)
    series_definitions = [
        (name = "below half", label = L"\rho = 1/2 - 1/L", marker = :circle),
        (name = "half filling", label = L"\rho = 1/2", marker = :utriangle),
        (name = "above half", label = L"\rho = 1/2 + 1/L", marker = :dtriangle),
    ]

    plt = plot(
        xlabel = L"L",
        ylabel = nothing,
        left_margin = 15mm,
        bottom_margin = 0mm,
        right_margin = 5mm,
        xlims = (15, 27),
        ylims = (0, 1.05),
        xticks = 16:2:26,
        yticks = 0:0.25:1.0,
        legend = (0.65, 0.5),
        foreground_color_legend = :transparent,
        size = (420, 285),
    )
    annotate!(plt, 12, 0.525, text(L"\frac{|\mathcal{K}_{\max}|}{|\mathcal{C}|}", 11))

    for series in series_definitions
        selected = sort([row for row in rows if row.series == series.name]; by = row -> row.L)
        isempty(selected) && continue
        x = [row.L for row in selected]
        y = [row.ratio for row in selected]

        plot!(
            plt,
            x,
            y;
            lw = 2,
            ls = :solid,
            c = :black,
            label = nothing,
        )
        scatter!(
            plt,
            x,
            y;
            markershape = series.marker,
            markersize = 5.5,
            markercolor = :black,
            markerstrokecolor = :black,
            markerstrokewidth = 1.4,
            c = :black,
            label = series.label,
        )
    end

    mkpath(output_dir)
    savefig(plt, pdf_path)
    println("Saved fig8.pdf to ", pdf_path)
    return plt
end

rows = read_fig8_data(data_path)
plot_fig8(rows)
