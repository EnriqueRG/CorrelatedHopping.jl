using CorrelatedHopping
using LaTeXStrings
using Measures
using Plots
using Random

###########################################################
# Fig. 2                                                  #
# Active sites for a small system.                        #
###########################################################

function particle_offsets(count::Integer)
    d = 0.18
    if count == 1 || count >= 7
        return [(0.0, 0.0)]
    elseif count == 2
        return [(-d, -d), (d, d)]
    elseif count == 3
        y_shift = 0.06
        dr = 1.15 * d
        return [
            (0.0, -2 * dr * sqrt(3) / 3 + y_shift),
            (-dr, dr * sqrt(3) / 3 + y_shift),
            (dr, dr * sqrt(3) / 3 + y_shift),
        ]
    elseif count == 4
        return [(-d, -d), (d, -d), (-d, d), (d, d)]
    elseif count == 5
        return [(-d, -d), (d, -d), (0.0, 0.0), (-d, d), (d, d)]
    elseif count == 6
        return [(-d, -d), (d, -d), (-d, 0.0), (d, 0.0), (-d, d), (d, d)]
    else
        return Tuple{Float64,Float64}[]
    end
end

function plot_active_intervals(
    state_history::AbstractMatrix{<:Integer}, 
    active_sites_history::AbstractMatrix{Bool}
)
    t, L = size(state_history)

    # The heatmap natively maps matrix rows to the y-axis and columns to the x-axis.
    h = heatmap(
        active_sites_history,
        c = [:white, :lightgray],
        yflip = true, # Orient time to flow downwards
        aspect_ratio = :equal,
        ticks = false,
        colorbar = false,
        axis = false,
        # Scale both width and height based on the lattice size and time steps
        size = (max(360, Int(L * 26)), max(90, Int(t * 26))),
        xlims = (-1.0, L + 0.5),
        ylims = (0.5, t + 0.5),
    )

    # Generate horizontal grid lines for each time step
    for y in 0.5:1:(t + 0.5)
        plot!(h, [0.5, L + 0.5], [y, y], color = :black, lw = 0.5, label = "")
    end
    # Generate vertical grid lines for each spatial site
    for x in 0.5:1:(L + 0.5)
        plot!(h, [x, x], [0.5, t + 0.5], color = :black, lw = 0.5, label = "")
    end

    single_xs = Float64[]
    single_ys = Float64[]
    pip_xs = Float64[]
    pip_ys = Float64[]
    big_xs = Float64[]
    big_ys = Float64[]

    for row in 1:t, col in 1:L
        count = state_history[row, col]
        count == 0 && continue

        if count == 1
            push!(single_xs, col)
            push!(single_ys, row)
        elseif count >= 7
            push!(big_xs, col)
            push!(big_ys, row)
        else
            for (dx, dy) in particle_offsets(count)
                push!(pip_xs, col + dx)
                push!(pip_ys, row + dy)
            end
        end
    end

    scatter!(
        h,
        single_xs,
        single_ys,
        color = :black,
        markersize = 5,
        markerstrokewidth = 0,
        legend = false,
    )
    scatter!(
        h,
        pip_xs,
        pip_ys,
        color = :black,
        markersize = 5,
        markerstrokewidth = 0,
        legend = false,
    )
    scatter!(
        h,
        big_xs,
        big_ys,
        color = :black,
        markersize = 8,
        markerstrokewidth = 0,
        legend = false,
    )

    # Arrow
    plot!(h, [-0.2, -0.2], [t, 1], arrow = :closed, color = :black, lw = 2, label = "")
    annotate!(h, -0.7, t / 2, text("Time", 12, :black, :center, rotation = 90))

    return h
end


# Simulate the evolution
rng = MersenneTwister(540)
L = 30
rho0 = 0.6
initial_state = Int.(rand(rng, L) .< rho0)
hop_rate = 1.0
reaction_rate = 1.0
reaction_model = Reaction(2, 0)
sys = CorrelatedHopping.initialize_system(
    L,
    initial_state,
    hop_rate,
    reaction_rate;
    reaction = reaction_model,
)
simulation = CorrelatedHopping.simulate!(
    sys,
    (sys, _t) -> is_final_binary(sys);
    record = :all,
    rng,
)
state_history = simulation.history

# Evaluate active sites and plot
active_sites_history = active_sites(
    state_history;
    algorithm = :montecarlo,
    reaction = reaction_model,
    hop_rate = hop_rate,
    reaction_rate = reaction_rate,
)
plt = plot_active_intervals(state_history, active_sites_history)

# Save figure
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)
path = joinpath(output_dir, "fig2.pdf")
savefig(plt, path)
println("Saved fig2.pdf to ", path)
