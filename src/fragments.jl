"""
Given an array of occupations (state), initializes the active sites as a boolean array based on hopping participation.
"""
function initialize_active_sites(state::AbstractVector{<:Integer})
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
Applies the rules from the correlated hopping model to find all active sites starting from an initial state.
"""
function search_active_sites(
    initial_state::AbstractVector{<:Integer},
    active_sites::AbstractVector{Bool} = initialize_active_sites(initial_state); 
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
Helper function. A sanbox is defined as a contiguous block of active sites, 
extended by one site on either side to capture boundary effects.
The idea is that sandboxes that do not overlap can be simulated independently, which reduces the state space.
"""
function get_sandboxes(active_sites::AbstractVector{Bool})
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
Sandbox-based search for active sites. Each sandbox is simulated independently.
Results are combined to update the global active sites and sandboxes.
"""
function sandbox_search(
    initial_state::AbstractVector{<:Integer},
    active_sites::AbstractVector{Bool} = initialize_active_sites(initial_state); 
    reaction::AbstractReactionModel = InstantaneousReaction(2, 0), 
    max_states::Int = 100000
)
    
    L = length(initial_state)
    is_occupied = BitVector(initial_state .> 0)
    while true
        sandboxes = get_sandboxes(active_sites .| is_occupied)
        state_changed = false

        for sandbox_sites in sandboxes
            sandbox = zeros(eltype(initial_state), L)
            sandbox[sandbox_sites] .= @view initial_state[sandbox_sites]

            new_active_sites = search_active_sites(
                sandbox, 
                active_sites; 
                reaction = reaction, 
                max_states = max_states
            )

            # Check for changes
            if new_active_sites != active_sites
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


"""
For evolution histories, we can propagate active sites backwards in time.
This saves some repeated computation.
"""
function sandbox_search_history(state_history::AbstractMatrix{<:Integer})
    t, L = size(state_history)
    active_sites_history = Matrix{Bool}(undef, t, L)
    active_sites_history[t, :] .= sandbox_search(state_history[t, :])

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
            active_sites_history[i, :] .= search_active_sites(
                state_history[i, :], 
                active_sites_history[i+1, :]
            )
        end
    end

    return active_sites_history
end