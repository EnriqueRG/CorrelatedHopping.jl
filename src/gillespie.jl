# Gillespie simulation core: reaction models, dynamics models, system initialization, and simulation

abstract type AbstractReactionModel end

function _validate_reaction_parameters(n, m; instantaneous::Bool)
    n isa Integer || throw(ArgumentError("n must be an integer."))
    m isa Integer || throw(ArgumentError("m must be an integer."))
    n >= 1 || throw(ArgumentError("n must be at least 1."))
    m >= 0 || throw(ArgumentError("m must be nonnegative."))
    if instantaneous && m >= n
        throw(ArgumentError("instantaneous reactions require m < n."))
    end
    return nothing
end

"""
    Reaction(n, m)

Finite-rate reaction `nA -> mA`. Finite reactions require `n >= 1` and
`m >= 0`; they may either remove particles (`m < n`) or create particles
(`m > n`).
"""
struct Reaction{N,M} <: AbstractReactionModel
    function Reaction{N,M}() where {N,M}
        _validate_reaction_parameters(N, M; instantaneous = false)
        return new{N,M}()
    end
end

# Convenience constructor that stores `n` and `m` as type parameters.
Reaction(n::Int, m::Int) = Reaction{n,m}()

"""
    InstantaneousReaction(n, m)

Infinite-rate reaction `nA -> mA`, enforced immediately after a dynamics event.
Instantaneous reactions require `n >= 1`, `m >= 0`, and `m < n` so that
reaction enforcement strictly reduces local occupation.
"""
struct InstantaneousReaction{N,M} <: AbstractReactionModel
    function InstantaneousReaction{N,M}() where {N,M}
        _validate_reaction_parameters(N, M; instantaneous = true)
        return new{N,M}()
    end
end

# Convenience constructor that stores `n` and `m` as type parameters.
InstantaneousReaction(n::Int, m::Int) = InstantaneousReaction{n,m}()

"""
    get_reaction_rate(reaction, rate, occ)

Return the local reaction propensity for one site with occupation `occ`.
"""
@inline function get_reaction_rate(::Reaction{N,M}, rate::Float64, occ::Int) where {N,M}
    occ < N && return 0.0
    combinations = 1.0
    for i in 0:(N - 1)
        combinations *= occ - i
    end
    return rate * combinations / factorial(N)
end

# Instantaneous reactions are enforced after motion, not sampled as events.
@inline get_reaction_rate(::InstantaneousReaction, ::Float64, ::Int) = 0.0

"""
    execute_reaction!(reaction, occupations, i)

Apply one finite-rate reaction event at site `i` and return the particle change.
"""
@inline function execute_reaction!(::Reaction{N,M}, occupations::Vector{Int}, i::Int) where {N,M}
    change = M - N
    occupations[i] += change
    return change
end

"""
    enforce_instantaneous!(reaction, occupations, i)

Apply instantaneous reactions at site `i` after a dynamics event.
"""
# Finite-rate reactions have nothing to enforce immediately after hopping.
@inline enforce_instantaneous!(::Reaction, occupations, i) = 0

# Collapse an over-occupied site according to the instantaneous reaction rule.
@inline function enforce_instantaneous!(
    ::InstantaneousReaction{N,M},
    occupations,
    i,
) where {N,M}
    loss = 0
    if occupations[i] >= N
        old = occupations[i]
        occupations[i] -= N - M
        loss -= N - M

        if occupations[i] >= N
            new = M + (old - M) % (N - M)
            loss += new - (old + loss)
            occupations[i] = new
        end
    end
    return loss
end

# Dynamics models

abstract type AbstractDynamicsModel end

"""
    CorrelatedHoppingDynamics()

Constrained four-site dynamics `0110 <-> 1001` on a periodic 1D lattice.
"""
struct CorrelatedHoppingDynamics <: AbstractDynamicsModel end

"""
    StandardDiffusion()

Independent nearest-neighbor diffusion on a periodic 1D lattice.
"""
struct StandardDiffusion <: AbstractDynamicsModel end

# Periodic one-dimensional index helpers.
@inline left(i, L) = i == 1 ? L : i - 1
@inline right(i, L) = i == L ? 1 : i + 1
@inline wrap(i, L) = mod1(i, L)

"""
    get_dynamics_rate(dynamics, sys, i)

Return the total hopping propensity associated with anchor site `i`.
"""
# Nearest-neighbor diffusion can hop left or right from the same site.
@inline function get_dynamics_rate(::StandardDiffusion, sys, i)
    return 2.0 * sys.hop_rate * sys.occupations[i]
end

"""
    execute_dynamics!(dynamics, sys, site_idx, target)

Sample and apply one dynamics event anchored at `site_idx`.
"""
# Use `target` to choose left versus right diffusion from `site_idx`.
@inline function execute_dynamics!(::StandardDiffusion, sys, site_idx, target)
    rate_left = sys.hop_rate * sys.occupations[site_idx]
    sys.occupations[site_idx] -= 1
    target_idx = target < rate_left ? left(site_idx, sys.L) : right(site_idx, sys.L)
    sys.occupations[target_idx] += 1
    return enforce_instantaneous!(sys.reaction, sys.occupations, target_idx)
end

"""
    update_neighborhood!(dynamics, sys, site_idx)

Refresh local rates affected by an event at `site_idx`.
"""
# Diffusion changes rates only near the departure and arrival sites.
@inline function update_neighborhood!(::StandardDiffusion, sys, site_idx)
    for k in (site_idx - 1):(site_idx + 1)
        update_local_rates!(sys, k)
    end
end

# Correlated hopping sums outward and inward pair-hop propensities.
@inline function get_dynamics_rate(::CorrelatedHoppingDynamics, sys, i)
    idx_p1 = right(i, sys.L)
    idx_p2 = right(idx_p1, sys.L)
    idx_p3 = right(idx_p2, sys.L)

    r_exp = sys.hop_rate * sys.occupations[idx_p1] * sys.occupations[idx_p2]
    r_con = sys.hop_rate * sys.occupations[i] * sys.occupations[idx_p3]
    return r_exp + r_con
end

# Use `target` to choose between the two correlated-hop directions.
@inline function execute_dynamics!(::CorrelatedHoppingDynamics, sys, site_idx, target)
    idx_p1 = right(site_idx, sys.L)
    idx_p2 = right(idx_p1, sys.L)
    idx_p3 = right(idx_p2, sys.L)
    r_exp = sys.hop_rate * sys.occupations[idx_p1] * sys.occupations[idx_p2]

    change = 0
    if target < r_exp
        sys.occupations[site_idx] += 1
        sys.occupations[idx_p1] -= 1
        sys.occupations[idx_p2] -= 1
        sys.occupations[idx_p3] += 1
        change += enforce_instantaneous!(sys.reaction, sys.occupations, site_idx)
        change += enforce_instantaneous!(sys.reaction, sys.occupations, idx_p3)
    else
        sys.occupations[site_idx] -= 1
        sys.occupations[idx_p1] += 1
        sys.occupations[idx_p2] += 1
        sys.occupations[idx_p3] -= 1
        change += enforce_instantaneous!(sys.reaction, sys.occupations, idx_p1)
        change += enforce_instantaneous!(sys.reaction, sys.occupations, idx_p2)
    end
    return change
end

# A four-site hop can affect rates up to three anchors away.
@inline function update_neighborhood!(::CorrelatedHoppingDynamics, sys, site_idx)
    for k in (site_idx - 3):(site_idx + 3)
        update_local_rates!(sys, k)
    end
end

# System initialization

"""
    DirectLatticeSystem

Mutable simulation state for a periodic one-dimensional lattice. The struct
stores occupations, rates, the reaction and dynamics models, and the binary
tree used to sample Gillespie events.
"""
struct DirectLatticeSystem{
    R<:AbstractReactionModel,
    D<:AbstractDynamicsModel,
}
    L::Int
    occupations::Vector{Int}
    hop_rate::Float64
    reaction_rate::Float64
    reaction::R
    dynamics::D
    site_rates::Vector{Float64}
    tree_capacity::Int
    tree::Vector{Float64}
end

"""
    _dynamics_from_symbol(dynamics)

Convert a user-facing dynamics symbol into its model object.
"""
function _dynamics_from_symbol(dynamics::Symbol)
    dynamics == :pairwise && return CorrelatedHoppingDynamics()
    dynamics == :diffusion && return StandardDiffusion()
    error("Unknown dynamics type: $dynamics. Use :pairwise or :diffusion.")
end

"""
    _validate_dynamics_size(dynamics, L)

Check the minimum lattice size required by a dynamics model.
"""
function _validate_dynamics_size(::CorrelatedHoppingDynamics, L::Int)
    L >= 4 || throw(ArgumentError("Pairwise dynamics require L >= 4."))
end

# Standard diffusion needs at least two sites on a periodic chain.
function _validate_dynamics_size(::StandardDiffusion, L::Int)
    L >= 2 || throw(ArgumentError("Standard diffusion requires L >= 2."))
end

"""
    initialize_system(L, initial_sites, hop_rate=1.0, reaction_rate=1.0; reaction, dynamics)

Create a periodic one-dimensional lattice system for Gillespie simulation.
`initial_sites` gives the integer occupation number at each site. `hop_rate`
sets the dynamics rate and `reaction_rate` sets the local reaction rate.

`reaction` defaults to `Reaction(2, 0)` and `dynamics` defaults to
`CorrelatedHoppingDynamics()`. The `dynamics` keyword also accepts `:pairwise` and
`:diffusion`.
"""
function initialize_system(
    L::Int,
    initial_sites::AbstractVector{<:Integer},
    hop_rate::Real = 1.0,
    reaction_rate::Real = 1.0;
    reaction::AbstractReactionModel = Reaction(2, 0),
    dynamics::Union{AbstractDynamicsModel,Symbol} = CorrelatedHoppingDynamics(),
)
    L > 0 || throw(ArgumentError("L must be positive."))
    length(initial_sites) == L ||
        throw(ArgumentError("initial_sites must have length L."))
    hop_rate >= 0 || throw(ArgumentError("hop_rate must be nonnegative."))
    reaction_rate >= 0 || throw(ArgumentError("reaction_rate must be nonnegative."))
    any(<(0), initial_sites) &&
        throw(ArgumentError("initial_sites must contain nonnegative occupations."))

    dynamics_model = dynamics isa Symbol ? _dynamics_from_symbol(dynamics) : dynamics
    _validate_dynamics_size(dynamics_model, L)
    initial_state = Int.(initial_sites)
    tree_capacity = nextpow(2, L)
    tree = zeros(Float64, 2 * tree_capacity)
    site_rates = zeros(Float64, L)

    sys = DirectLatticeSystem(
        L,
        initial_state,
        Float64(hop_rate),
        Float64(reaction_rate),
        reaction,
        dynamics_model,
        site_rates,
        tree_capacity,
        tree,
    )

    for i in 1:L
        enforce_instantaneous!(sys.reaction, sys.occupations, i)
    end
    for i in 1:L
        update_local_rates!(sys, i)
    end
    return sys
end

"""
    update_local_rates!(sys, i)

Recompute the total event rate at site `i` and update the binary rate tree.
Indices are periodic, so values outside `1:sys.L` are wrapped.
"""
function update_local_rates!(sys::DirectLatticeSystem, i::Int)
    i = wrap(i, sys.L)
    r_react = get_reaction_rate(sys.reaction, sys.reaction_rate, sys.occupations[i])
    total_site_rate = r_react + get_dynamics_rate(sys.dynamics, sys, i)
    sys.site_rates[i] = total_site_rate

    tree_idx = sys.tree_capacity + i - 1
    diff = total_site_rate - sys.tree[tree_idx]
    while tree_idx >= 1
        sys.tree[tree_idx] += diff
        tree_idx >>= 1
    end
    return sys
end

# Simulation

"""
    simulate!(
        sys,
        should_stop;
        record = :particle_changes,
        stop_on_reaction_only = true,
        rng = Random.default_rng(),
        event_callback = nothing,
    )

Run a Gillespie simulation in place until `should_stop(sys, t)` is true or no
events remain. Returns a named tuple with `final_time`, `recorded_times`,
`history`, and `final_state`.

The `record` keyword controls stored history:

- `:particle_changes`: record total particle number initially and after
  particle-changing events.
- `:all`: record the full lattice state initially and after every event.
- `:none`: do not store trajectory history.

When `stop_on_reaction_only=true`, the stopping rule is checked initially and
after particle-changing events. Set it to `false` for stopping rules that depend
on diffusion-only moves or elapsed time.

If `event_callback` is provided, it is called after every sampled event has been
applied and local rates have been refreshed:
`event_callback(sys, event_kind, site_idx, t, particle_change)`. `event_kind` is
`:reaction` for a finite reaction event and `:dynamics` for a hopping event.
"""
function simulate!(
    sys::DirectLatticeSystem,
    should_stop::Function;
    record::Symbol = :particle_changes,
    stop_on_reaction_only::Bool = true,
    rng::AbstractRNG = Random.default_rng(),
    event_callback::Union{Nothing,Function} = nothing,
)
    record in (:particle_changes, :all, :none) ||
        throw(ArgumentError("record must be :particle_changes, :all, or :none."))

    t = 0.0
    recorded_times = record === :none ? Float64[] : Float64[0.0]
    history = if record === :all
        Vector{Int}[copy(sys.occupations)]
    elseif record === :particle_changes
        Int[sum(sys.occupations)]
    else
        nothing
    end
    check_stop = true

    while true
        check_stop && should_stop(sys, t) && break

        total_rate = sys.tree[1]
        total_rate <= 1e-12 && break

        t += randexp(rng) / total_rate

        target = rand(rng) * total_rate
        idx = 1
        while idx < sys.tree_capacity
            left_child = 2 * idx
            left_child > length(sys.tree) && break
            val_left = sys.tree[left_child]
            if target <= val_left
                idx = left_child
            else
                target -= val_left
                idx = left_child + 1
            end
        end
        site_idx = idx - sys.tree_capacity + 1

        r_react = get_reaction_rate(
            sys.reaction,
            sys.reaction_rate,
            sys.occupations[site_idx],
        )
        local_target = rand(rng) * sys.site_rates[site_idx]
        particle_change = 0
        event_kind = :reaction

        if local_target < r_react
            particle_change += execute_reaction!(sys.reaction, sys.occupations, site_idx)
        else
            event_kind = :dynamics
            shifted_target = local_target - r_react
            particle_change += execute_dynamics!(sys.dynamics, sys, site_idx, shifted_target)
        end

        update_neighborhood!(sys.dynamics, sys, site_idx)
        if event_callback !== nothing
            event_callback(sys, event_kind, site_idx, t, particle_change)
        end

        if record === :all
            push!(recorded_times, t)
            push!(history, copy(sys.occupations))
        elseif record === :particle_changes && particle_change != 0
            push!(recorded_times, t)
            push!(history, history[end] + particle_change)
        end

        check_stop = !stop_on_reaction_only || particle_change != 0
    end

    recorded_history = record === :all ? stack(history, dims = 1) : history
    return (;
        final_time = t,
        recorded_times,
        history = recorded_history,
        final_state = copy(sys.occupations),
    )
end

# Fragment and final-state analysis: binary patterns and ensemble sampling

const FINAL_BINARY_FORBIDDEN_PATTERNS = (
    (1, 1, 1),
    (1, 1, 0, 1),
    (1, 0, 1, 1),
    (1, 1, 0, 0, 1),
    (1, 1, 0, 0, 0, 1, 1),
    (1, 0, 0, 1, 1),
    (1, 0, 0, 1, 0, 1),
    (1, 0, 1, 0, 0, 1),
    (1, 1, 0, 0, 0, 1, 0, 1),
    (1, 0, 1, 0, 0, 0, 1, 1),
)

"""
    is_final_binary(sys_or_state)

Return `true` when a binary periodic state contains none of the forbidden
motifs for the constrained pairwise dynamics. States with occupation greater
than one are not final. The exact forbidden-pattern check assumes `L > 7`.
"""
function is_final_binary(sys::DirectLatticeSystem)
    (isa(sys.reaction, Reaction{2}) || isa(sys.reaction, InstantaneousReaction{2})) || throw(ArgumentError("is_final_binary requires reactions with two reactants."))
    isa(sys.dynamics, CorrelatedHoppingDynamics) || throw(ArgumentError("is_final_binary requires CorrelatedHoppingDynamics dynamics."))
    return is_final_binary(sys.occupations)
end

# Raw-vector overload used by simulations, spectral lumping, and tests.
function is_final_binary(state::AbstractVector{<:Integer})
    L = length(state)
    L > 7 || throw(ArgumentError("Exact binary final-state patterns assume L > 7."))

    @inbounds for i in 1:L
        state[i] > 1 && return false

        for pattern in FINAL_BINARY_FORBIDDEN_PATTERNS
            match = true
            for j in 1:length(pattern)
                idx = i + j - 1
                idx = idx > L ? idx - L : idx
                if state[idx] != pattern[j]
                    match = false
                    break
                end
            end
            match && return false
        end
    end
    return true
end

# Ensemble sampling

"""
    run_ensemble(L, n_ensemble, rho0; hop_rate=1.0, reaction_rate=1.0, reaction, dynamics, condition_even=false, stop=nothing, rng=Random.default_rng())

Run `n_ensemble` independent simulations for a single system size `L`.
Initial states are Bernoulli binary states with density `rho0`. Returns a
vector of final times. `hop_rate` sets the dynamics rate and `reaction_rate`
sets the local reaction rate.

By default, simulations stop at `is_final_binary`. Pass `stop=(sys, t) -> ...`
to use another stopping rule, such as ordinary diffusion stopping when fewer
than two particles remain.
"""
function run_ensemble(
    L::Integer,
    n_ensemble::Int,
    rho0;
    hop_rate::Real = 1.0,
    reaction_rate::Real = 1.0,
    reaction::AbstractReactionModel = Reaction(2, 0),
    dynamics::Union{AbstractDynamicsModel,Symbol} = CorrelatedHoppingDynamics(),
    condition_even::Bool = false,
    stop_on_reaction_only::Bool = true,
    stop = nothing,
    rng::AbstractRNG = Random.default_rng(),
)
    n_ensemble > 0 || throw(ArgumentError("n_ensemble must be positive."))
    0 <= rho0 <= 1 || throw(ArgumentError("rho0 must be between 0 and 1."))
    hop_rate >= 0 || throw(ArgumentError("hop_rate must be nonnegative."))
    reaction_rate >= 0 || throw(ArgumentError("reaction_rate must be nonnegative."))
    stop_rule = isnothing(stop) ? ((s, _t) -> is_final_binary(s)) : stop

    L > 0 || throw(ArgumentError("L must be positive."))
    times = zeros(Float64, n_ensemble)
    for j in 1:n_ensemble
        initial_state = Int.(rand(rng, Int(L)) .< rho0)
        if condition_even && isodd(sum(initial_state))
            idx = rand(rng, 1:Int(L))
            initial_state[idx] = 1 - initial_state[idx]
        end

        sys = initialize_system(
            Int(L),
            initial_state,
            Float64(hop_rate),
            Float64(reaction_rate);
            reaction,
            dynamics,
        )
        result = simulate!(
            sys,
            stop_rule;
            record = :none,
            stop_on_reaction_only,
            rng,
        )
        times[j] = result.final_time
    end

    return times
end
