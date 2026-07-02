using CorrelatedHopping
using Printf

###########################################################
# Fig. 8 data generation                                  #
# Krylov-sector connectivity near half filling.           #
###########################################################

output_dir = joinpath(@__DIR__, "output")
data_path = joinpath(output_dir, "fig8_data.csv")

function transition_density_tasks(profile::AbstractString = "local")
    profile in ("local", "cluster") ||
        throw(ArgumentError("profile must be \"local\" or \"cluster\"."))

    local_tasks = [
        (series = "below half", L = 16, N = 7),
        (series = "below half", L = 18, N = 8),
        (series = "below half", L = 20, N = 9),
        (series = "below half", L = 22, N = 10),
        (series = "below half", L = 24, N = 11),
        (series = "below half", L = 26, N = 12),

        (series = "half filling", L = 16, N = 8),
        (series = "half filling", L = 18, N = 9),
        (series = "half filling", L = 20, N = 10),
        (series = "half filling", L = 22, N = 11),
        (series = "half filling", L = 24, N = 12),

        (series = "above half", L = 16, N = 9),
        (series = "above half", L = 18, N = 10),
        (series = "above half", L = 20, N = 11),
        (series = "above half", L = 22, N = 12),
    ]

    profile == "local" && return local_tasks

    return [
        local_tasks...,
        (series = "half filling", L = 26, N = 13),
        (series = "above half", L = 24, N = 13),
        (series = "above half", L = 26, N = 14),
    ]
end

function generate_fig8_data(; profile::AbstractString = "local", path::String = data_path)
    mkpath(dirname(path))

    completed = Set{Tuple{String,Int,Int}}()
    if isfile(path)
        for (line_number, line) in enumerate(eachline(path))
            line_number == 1 && continue
            isempty(strip(line)) && continue
            fields = split(line, ','; limit = 4)
            length(fields) >= 3 || error("Could not read CSV key on line $line_number.")
            push!(completed, (String(fields[1]), parse(Int, fields[2]), parse(Int, fields[3])))
        end
        @printf("Loaded %d completed row(s) from %s\n", length(completed), path)
    else
        open(path, "w") do io
            println(io, "series,L,N,rho,ratio,note")
        end
    end

    open(path, "a") do io
        for task in transition_density_tasks(profile)
            key = (task.series, task.L, task.N)
            if key in completed
                @printf("Skipping %s L=%d N=%d; already in %s\n", task.series, task.L, task.N, path)
                continue
            end

            result = CorrelatedHopping.largest_krylov_sector(task.L, task.N)
            note = "symmetry_sector=(P=$(result.P); A=$(result.A)); symmetry_sector_states=$(result.symmetry_sector_size); number_of_symmetry_sectors=$(result.number_of_symmetry_sectors); krylov_sectors_seen=$(result.krylov_sectors_seen); remaining=$(result.remaining); exhausted=$(result.exhausted)"

            @printf(
                io,
                "%s,%d,%d,%.17g,%.17g,%s\n",
                task.series,
                task.L,
                task.N,
                result.rho,
                result.ratio,
                note,
            )
            flush(io)
            push!(completed, key)
            @printf("Saved %s L=%d N=%d to %s\n", task.series, task.L, task.N, path)
        end
    end

    return path
end

generate_fig8_data(profile = get(ENV, "FIG8_PROFILE", "local"))
