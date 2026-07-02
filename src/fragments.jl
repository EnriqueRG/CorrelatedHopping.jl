"""
    active_sites(
        state_or_history;
        algorithm = :sandbox,
        reaction = InstantaneousReaction(2, 0),
        max_states = 100000,
        n_sims = 1000,
        max_steps = nothing,
        time_horizon = nothing,
        hop_rate = 1.0,
        final_time_horizon = 5 / hop_rate,
        activity_threshold = 0.0,
        reaction_rate = 1.0,
        rng = Random.default_rng(),
    )

Find active sites for a single state or an evolution history.

A single state is treated as a one-step history and returns a vector. A history
matrix with one row per time step returns a matrix with the same shape. The
sandbox backend propagates active sites backwards in time when possible to avoid
repeated searches. Currently supported algorithms are:

- `:sandbox`: sandboxed state-space search.
- `:montecarlo`: Monte Carlo estimate from repeated stochastic dynamics.

Arguments:

- `state_or_history`: a vector of site occupations for one state, or a matrix
  whose rows are successive states in an evolution history.
- `algorithm`: activity estimator. Use `:sandbox` for deterministic sandboxed
  state-space search, or `:montecarlo` for sampled stochastic activity.
- `reaction`: reaction model used by the estimator. The sandbox backend enforces
  instantaneous reactions during searched hops; the Monte Carlo backend passes
  the model to the stochastic simulator.
- `max_states`: maximum number of states explored by each sandbox search before
  warning and returning the activity found so far. Only used by `:sandbox`.
- `n_sims`: number of stochastic simulations run from each state. Only used by
  `:montecarlo`.
- `max_steps`: optional maximum number of Gillespie steps per stochastic
  simulation. `nothing` means no artificial step cap. Only used by
  `:montecarlo`.
- `time_horizon`: optional physical-time horizon for each stochastic simulation.
  `nothing` means simulations stop at `is_final_binary`. Only used by
  `:montecarlo`.
- `hop_rate`: correlated-hopping rate used in stochastic simulations. It also
  sets the default final-state search window through `final_time_horizon =
  5 / hop_rate`. Only used by `:montecarlo`.
- `final_time_horizon`: optional physical-time horizon used only when a
  Monte Carlo simulation starts from a final binary state and `time_horizon` is
  `nothing`. This allows final states to be checked for residual hopping
  activity. The default is five microscopic hopping times, `5 / hop_rate`; pass
  `nothing` to keep the immediate final-state stop. Only used by `:montecarlo`.
- `activity_threshold`: minimum average hop count needed to mark a site active.
  Only used by `:montecarlo`.
- `reaction_rate`: finite reaction rate used in stochastic simulations. This has
  no effect for instantaneous reactions. Only used by `:montecarlo`.
- `rng`: random number generator used by stochastic simulations. Only used by
  `:montecarlo`.
"""
function active_sites(
    initial_state::AbstractVector{<:Integer};
    algorithm::Symbol = :sandbox,
    reaction::AbstractReactionModel = InstantaneousReaction(2, 0),
    max_states::Int = 100000,
    n_sims::Integer = 1000,
    max_steps::Union{Nothing,Integer} = nothing,
    time_horizon::Union{Nothing,Real} = nothing,
    hop_rate::Real = 1.0,
    final_time_horizon::Union{Nothing,Real} = 5 / hop_rate,
    activity_threshold::Real = 0.0,
    reaction_rate::Real = 1.0,
    rng::AbstractRNG = Random.default_rng(),
)
    state_history = reshape(initial_state, 1, length(initial_state))
    return vec(active_sites(
        state_history;
        algorithm = algorithm,
        reaction = reaction,
        max_states = max_states,
        n_sims = n_sims,
        max_steps = max_steps,
        time_horizon = time_horizon,
        hop_rate = hop_rate,
        final_time_horizon = final_time_horizon,
        activity_threshold = activity_threshold,
        reaction_rate = reaction_rate,
        rng = rng,
    ))
end

function active_sites(
    state_history::AbstractMatrix{<:Integer};
    algorithm::Symbol = :sandbox,
    reaction::AbstractReactionModel = InstantaneousReaction(2, 0),
    max_states::Int = 100000,
    n_sims::Integer = 1000,
    max_steps::Union{Nothing,Integer} = nothing,
    time_horizon::Union{Nothing,Real} = nothing,
    hop_rate::Real = 1.0,
    final_time_horizon::Union{Nothing,Real} = 5 / hop_rate,
    activity_threshold::Real = 0.0,
    reaction_rate::Real = 1.0,
    rng::AbstractRNG = Random.default_rng(),
)
    if algorithm === :sandbox
        return _sandbox_active_sites(
            state_history;
            reaction = reaction,
            max_states = max_states,
        )
    elseif algorithm === :montecarlo
        return _montecarlo_active_sites(
            state_history;
            reaction = reaction,
            n_sims = n_sims,
            max_steps = max_steps,
            time_horizon = time_horizon,
            hop_rate = hop_rate,
            final_time_horizon = final_time_horizon,
            activity_threshold = activity_threshold,
            reaction_rate = reaction_rate,
            rng = rng,
        )
    else
        throw(ArgumentError("Unsupported active-sites algorithm: $algorithm. Use :sandbox or :montecarlo."))
    end
end


# Sandbox state-space search

"""
    _sandbox_active_sites(state_history; reaction, max_states)

Evaluate active sites for every row of a history using the sandboxed
state-space search. The last row is searched directly, then activity is
propagated backwards when adjacent rows have the same particle number or when
all later sites are already active.
"""
function _sandbox_active_sites(
    state_history::AbstractMatrix{<:Integer};
    reaction::AbstractReactionModel,
    max_states::Int,
)
    t, L = size(state_history)
    active_sites_history = Matrix{Bool}(undef, t, L)
    active_sites_history[t, :] .= _sandbox_search(
        state_history[t, :];
        reaction = reaction,
        max_states = max_states,
    )

    # Propagate backwards in time
    for i in (t-1):-1:1
        # If the number of particles did not change, no need to recompute
        if sum(@view state_history[i, :]) == sum(@view state_history[i+1, :])
            active_sites_history[i, :] .= @view active_sites_history[i+1, :]
        # If all sites are active, no need to recompute
        elseif all(@view active_sites_history[i+1, :])
             active_sites_history[i, :] .= 1
        # Otherwise, recompute
        else
            active_sites_history[i, :] .= _search_active_sites(
                state_history[i, :], 
                active_sites_history[i+1, :];
                reaction = reaction,
                max_states = max_states,
            )
        end
    end

    return active_sites_history
end


"""
    _initialize_active_sites(state)

Return the sites that are active in the input state before any state-space
exploration. A site is marked active when it belongs to a four-site pattern that
can immediately participate in a correlated hop.
"""
function _initialize_active_sites(state::AbstractVector{<:Integer})
    L = length(state)
    @assert L >= 4 "Correlated hopping requires L >= 4."

    is_occupied = BitVector(state .> 0)
    active_sites = falses(L)

    # Push all possible hops into the list of active sites.
    for i1 in 1:L
        i2 = right(i1, L)
        i3 = right(i2, L)
        i4 = right(i3, L)
        if (is_occupied[i2] && is_occupied[i3]) || (is_occupied[i1] && is_occupied[i4])
            active_sites[i1] = true
            active_sites[i2] = true
            active_sites[i3] = true
            active_sites[i4] = true
        end
    end
    
    return active_sites
end


"""
    _search_active_sites(initial_state, active_sites=_initialize_active_sites(initial_state); reaction, max_states)

Breadth-first search over states reachable from `initial_state` under correlated
hopping and instantaneous reaction enforcement. The supplied `active_sites`
vector is updated as new hopping supports are discovered. Search stops early if
all sites are active or if more than `max_states` states have been visited.
"""
function _search_active_sites(
    initial_state::AbstractVector{<:Integer},
    active_sites::AbstractVector{Bool} = _initialize_active_sites(initial_state); 
    reaction::AbstractReactionModel = InstantaneousReaction(2, 0), 
    max_states::Int = 100000
    )

    L = length(initial_state)
    
    # Track visited states
    visited = Set{Vector{Int64}}()
    push!(visited, initial_state)
    
    # Queue for exploration
    queue = [initial_state]

    while !isempty(queue)
        # Recursion depth break, exploring everything is generally impossible
        if length(visited) > max_states
            @warn "Reached maximum state limit; active sites may be underestimated." max_states = max_states
            break
        end

        # Early stop condition is all sites are active
        if all(active_sites .== 1)
            break
        end
        
        current_state = popfirst!(queue)
        
        for i in 1:L
            i1 = i
            i2 = right(i1, L)
            i3 = right(i2, L)
            i4 = right(i3, L)
            
            # Hop out
            if current_state[i2] > 0 && current_state[i3] > 0
                # Mark sites as active
                active_sites[[i1, i2, i3, i4]] .= 1
                
                # Move to new state
                new_state = copy(current_state)
                new_state[[i1, i2, i3, i4]] .+= (1, -1, -1, 1)
                
                # Enforce reactions depending on type
                CorrelatedHopping.enforce_instantaneous!(reaction, new_state, i1)
                CorrelatedHopping.enforce_instantaneous!(reaction, new_state, i4)
                
                # If it's a new state, queue it for exploration
                if new_state ∉ visited
                    push!(visited, new_state)
                    push!(queue, new_state)
                end
            end
            
            # Hop in
            if current_state[i1] > 0 && current_state[i4] > 0
                active_sites[[i1, i2, i3, i4]] .= 1
                new_state = copy(current_state)
                new_state[[i1, i2, i3, i4]] .+= (-1, 1, 1, -1)
                CorrelatedHopping.enforce_instantaneous!(reaction, new_state, i2)
                CorrelatedHopping.enforce_instantaneous!(reaction, new_state, i3)
                if new_state ∉ visited
                    push!(visited, new_state)
                    push!(queue, new_state)
                end
            end
        end
    end
    
    return active_sites
end


"""
    _get_sandboxes(active_sites)

Build the sandbox regions used by `_sandbox_search`. A sandbox is a contiguous
periodic block obtained by dilating active or occupied sites to include nearby
boundary effects; disconnected sandboxes can then be searched independently.
"""
function _get_sandboxes(active_sites::AbstractVector{Bool})
    sandbox_mask = copy(active_sites)
    
    # Apply boolean dilation via sequential circular shifts
    for j in 1:3
        sandbox_mask .|= circshift(active_sites, j)
    end
    
    indices = findall(sandbox_mask)
    isempty(indices) && return Vector{Vector{Int}}()

    # Extract contiguous blocks
    breaks = [0; findall(diff(indices) .> 1); length(indices)]
    regions = [indices[breaks[i]+1 : breaks[i+1]] for i in 1:(length(breaks)-1)]

    # Enforce periodic boundaries on the extracted arrays
    if length(regions) > 1 && sandbox_mask[begin] && sandbox_mask[end]
        prepend!(regions[begin], pop!(regions))
    end

    return regions
end


"""
    _sandbox_search(initial_state, active_sites=_initialize_active_sites(initial_state); reaction, max_states)

Search for active sites by repeatedly splitting the state into local sandboxes.
Each sandbox is explored with `_search_active_sites`; when any sandbox expands
the global activity mask, the sandbox decomposition is rebuilt and the process
continues.
"""
function _sandbox_search(
    initial_state::AbstractVector{<:Integer},
    active_sites::AbstractVector{Bool} = _initialize_active_sites(initial_state); 
    reaction::AbstractReactionModel = InstantaneousReaction(2, 0), 
    max_states::Int = 100000
)
    
    L = length(initial_state)
    is_occupied = BitVector(initial_state .> 0)
    while true
        sandboxes = _get_sandboxes(active_sites .| is_occupied)
        state_changed = false

        for sandbox_sites in sandboxes
            sandbox = zeros(eltype(initial_state), L)
            sandbox[sandbox_sites] .= @view initial_state[sandbox_sites]
            old_active_sites = copy(active_sites)

            new_active_sites = _search_active_sites(
                sandbox, 
                copy(active_sites);
                reaction = reaction, 
                max_states = max_states
            )

            # Check for changes
            if new_active_sites != old_active_sites
                active_sites = new_active_sites
                state_changed = true
                break
            end
        end
        if !state_changed
            break
        end
    end
    
    return active_sites
end


# Monte Carlo stochastic activity search

"""
    _montecarlo_active_sites(state_history; reaction, n_sims, max_steps, time_horizon,
                             final_time_horizon, activity_threshold,
                             hop_rate, reaction_rate, rng)

Estimate active sites by Monte Carlo sampling. This converts the raw average
hop counts from `_montecarlo_activity` into a boolean matrix by marking sites
whose sampled activity is greater than `activity_threshold`.
"""
function _montecarlo_active_sites(
    state_history::AbstractMatrix{<:Integer};
    reaction::AbstractReactionModel,
    n_sims::Integer,
    max_steps::Union{Nothing,Integer},
    time_horizon::Union{Nothing,Real},
    final_time_horizon::Union{Nothing,Real},
    activity_threshold::Real,
    hop_rate::Real,
    reaction_rate::Real,
    rng::AbstractRNG,
)
    activity_threshold >= 0 ||
        throw(ArgumentError("activity_threshold must be nonnegative."))
    activity_history = _montecarlo_activity(
        state_history;
        reaction = reaction,
        n_sims = n_sims,
        max_steps = max_steps,
        time_horizon = time_horizon,
        final_time_horizon = final_time_horizon,
        hop_rate = hop_rate,
        reaction_rate = reaction_rate,
        rng = rng,
    )
    return Matrix{Bool}(activity_history .> activity_threshold)
end

"""
    _montecarlo_activity(state_history; reaction, n_sims, max_steps, time_horizon,
                         final_time_horizon, hop_rate, reaction_rate, rng)

Return the average correlated-hop count at each site for each row of a history.
For every input row, the function runs `n_sims` stochastic trajectories from
that row and averages the per-site hop counts. If `final_time_horizon` is set,
trajectories that start from a final binary state can still run for that time
window when no general `time_horizon` is supplied.
"""
function _montecarlo_activity(
    state_history::AbstractMatrix{<:Integer};
    reaction::AbstractReactionModel,
    n_sims::Integer,
    max_steps::Union{Nothing,Integer},
    time_horizon::Union{Nothing,Real},
    final_time_horizon::Union{Nothing,Real},
    hop_rate::Real,
    reaction_rate::Real,
    rng::AbstractRNG,
)
    n_sims > 0 || throw(ArgumentError("n_sims must be positive."))
    if max_steps !== nothing
        max_steps > 0 || throw(ArgumentError("max_steps must be positive."))
    end
    hop_rate >= 0 || throw(ArgumentError("hop_rate must be nonnegative."))
    reaction_rate >= 0 || throw(ArgumentError("reaction_rate must be nonnegative."))
    if time_horizon !== nothing
        time_horizon >= 0 || throw(ArgumentError("time_horizon must be nonnegative."))
    end
    if final_time_horizon !== nothing
        final_time_horizon >= 0 || throw(ArgumentError("final_time_horizon must be nonnegative."))
    end

    t, L = size(state_history)
    activity_history = zeros(Float64, t, L)

    for row in 1:t
        current_state = Int.(state_history[row, :])
        local_hops = zeros(Float64, L)

        for _ in 1:n_sims
            local_hops .+= _montecarlo_hop_counts(
                current_state;
                reaction = reaction,
                max_steps = max_steps,
                time_horizon = time_horizon,
                final_time_horizon = final_time_horizon,
                hop_rate = hop_rate,
                reaction_rate = reaction_rate,
                rng = rng,
            )
        end

        activity_history[row, :] .= local_hops ./ n_sims
    end

    return activity_history
end

"""
    _montecarlo_hop_counts(initial_state; reaction, max_steps, time_horizon,
                           final_time_horizon, hop_rate, reaction_rate, rng)

Run one stochastic correlated-hopping trajectory from `initial_state` with
`simulate!` and count how many sampled hop events involve each site. The
trajectory stops when there are no events left, when the optional `max_steps`
limit is reached, when `time_horizon` has been reached, or by default when the
state is final binary. When the initial state is already final binary and
`time_horizon === nothing`, `final_time_horizon` can provide one last
time-based search window for residual hopping activity.
"""
function _montecarlo_hop_counts(
    initial_state::AbstractVector{<:Integer};
    reaction::AbstractReactionModel,
    max_steps::Union{Nothing,Integer},
    time_horizon::Union{Nothing,Real},
    final_time_horizon::Union{Nothing,Real},
    hop_rate::Real,
    reaction_rate::Real,
    rng::AbstractRNG,
)
    L = length(initial_state)
    sys = initialize_system(
        L,
        Int.(initial_state),
        Float64(hop_rate),
        Float64(reaction_rate);
        reaction = reaction,
        dynamics = CorrelatedHoppingDynamics(),
    )

    hop_counts = zeros(Float64, L)
    step_count = Ref(0)
    simulation_horizon = time_horizon
    if simulation_horizon === nothing && is_final_binary(sys)
        simulation_horizon = final_time_horizon
    end

    function count_hop!(sys, event_kind, site_idx, _time, _particle_change)
        step_count[] += 1
        if event_kind === :dynamics
            _count_hop_sites!(hop_counts, sys.dynamics, site_idx, sys.L)
        end
        return nothing
    end

    function should_stop(_sys, t)
        if max_steps !== nothing && step_count[] >= max_steps
            return true
        end
        if simulation_horizon !== nothing
            return t >= simulation_horizon
        end
        return is_final_binary(_sys)
    end

    simulate!(
        sys,
        should_stop;
        record = :none,
        stop_on_reaction_only = false,
        rng = rng,
        event_callback = count_hop!,
    )

    return hop_counts
end

"""
    _count_hop_sites!(hop_counts, dynamics, site_idx, L)

Increment the four sites touched by a correlated hop anchored at `site_idx`.
This records participation in the hop, independent of whether the hop direction
is outward or inward.
"""
function _count_hop_sites!(
    hop_counts::AbstractVector{<:Real},
    ::CorrelatedHoppingDynamics,
    site_idx::Int,
    L::Int,
)
    i1 = site_idx
    i2 = right(i1, L)
    i3 = right(i2, L)
    i4 = right(i3, L)
    hop_counts[i1] += 1
    hop_counts[i2] += 1
    hop_counts[i3] += 1
    hop_counts[i4] += 1
    return hop_counts
end
