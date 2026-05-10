using CorrelatedHopping
using LaTeXStrings
using Measures
using Plots
using Random

###########################################################
# Fig. 3                                                  #
# Active sites for a small system.                        #
###########################################################

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

    # Findall on a matrix returns a vector of CartesianIndex(row, col)
    particle_indices = findall(x -> x != 0, state_history)
    
    # Extract spatial coordinates (columns) and temporal coordinates (rows)
    xs = [idx[2] for idx in particle_indices]
    ys = [idx[1] for idx in particle_indices]

    scatter!(
        h,
        xs,
        ys,
        color = :black,
        markersize = 5,
        markerstrokewidth = 0,
        legend = false,
    )

    # Arrow
    plot!(h, [-0.2, -0.2], [t, 1], arrow = :closed, color = :black, lw = 2, label = "")
    annotate!(h, -0.7, t / 2, text("Time", 12, :black, :center, rotation = 90))

    return h
end


# Simulate the evolution
rng = MersenneTwister(3)
L = 30
ρ0 = 0.6
initial_state = Int.(rand(rng, L) .< ρ0)
sys = CorrelatedHopping.initialize_system(L, initial_state, reaction = InstantaneousReaction(2, 0))
times, particles_history = CorrelatedHopping.simulate!(
    sys,
    (sys, _t) -> is_final_binary(sys);
    full_history = true,
    rng,
)

# Evaluate active sites and plot
active_sites_history = sandbox_search_history(particles_history);
plt = plot_active_intervals(particles_history, active_sites_history)

# Save figure
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)
path = joinpath(output_dir, "fig3.pdf")
savefig(plt, path)
println("Saved fig3.pdf to ", path)
