# Symmetry-sector enumeration and Krylov-sector connectivity tools for the
# four-site hopping rule 0110 <-> 1001 on an even periodic chain.

# To save memory, states are packed as a UInt128 using 4 bits per lattice site.
# Thus: 
# each site occupation must fit in 4 bits, so n_x <= 15
# the whole chain must fit in UInt128, so 4L <= 128
const KRYLOV_STATE_BITS = 4
const KRYLOV_STATE_MASK = UInt128(0x0f)


"""
    _krylov_getocc(state, site)

Return the occupation stored at `site` in a packed Krylov state.
"""
@inline function _krylov_getocc(state::UInt128, site::Int)
    return Int((state >> (KRYLOV_STATE_BITS * site)) & KRYLOV_STATE_MASK)
end


"""
    _krylov_addocc(state, site, delta)

Return a packed state with `delta` added to the occupation at `site`.
"""
@inline function _krylov_addocc(state::UInt128, site::Int, delta::Int)
    shift = KRYLOV_STATE_BITS * site
    old = Int((state >> shift) & KRYLOV_STATE_MASK)
    new = old + delta
    @assert 0 <= new <= Int(KRYLOV_STATE_MASK)
    return (state & ~(KRYLOV_STATE_MASK << shift)) | (UInt128(new) << shift)
end


"""
    _validate_krylov_size(L, N)

Check that the requested system fits the current `UInt128` state encoding.
"""
function _validate_krylov_size(L::Int, N::Int)
    iseven(L) || throw(ArgumentError("Symmetry-sector diagnostics require even L."))
    N <= Int(KRYLOV_STATE_MASK) ||
        throw(ArgumentError("Increase the state representation before using N > 15."))
    KRYLOV_STATE_BITS * L <= 128 ||
        throw(ArgumentError("Increase the packed state representation before using this L."))
    return nothing
end


"""
    symmetry_sector_counts(L, N)

Count all symmetry sectors `(P, A)` at fixed particle number `N` for an even periodic chain of length `L`.
Here `P` is the dipole modulo `L`, and `A = N_even - N_odd`.
"""
function symmetry_sector_counts(L::Int, N::Int)
    _validate_krylov_size(L, N)

    dp = Dict{Tuple{Int,Int,Int},BigInt}((0, 0, 0) => big(1))
    for x in 0:(L - 1)
        next = Dict{Tuple{Int,Int,Int},BigInt}()
        sign = iseven(x) ? 1 : -1
        for ((used, P, A), count) in dp
            for added in 0:(N - used)
                key = (used + added, mod(P + x * added, L), A + sign * added)
                next[key] = get(next, key, big(0)) + count
            end
        end
        dp = next
    end

    return Dict((P, A) => count for ((used, P, A), count) in dp if used == N)
end


"""
    largest_symmetry_sector(L, N)

Return `((P, A), symmetry_sector_size, number_of_symmetry_sectors)` for a largest symmetry sector.
Ties are broken by picking the sector with the smallest `P`, then `A`, values.
"""
function largest_symmetry_sector(L::Int, N::Int)
    counts = symmetry_sector_counts(L, N)
    best_key = minimum(keys(counts))
    best_count = counts[best_key]
    for key in sort(collect(keys(counts)))
        count = counts[key]
        if count > best_count
            best_key = key
            best_count = count
        end
    end
    return best_key, Int(best_count), length(counts)
end


"""
    _symmetry_suffix_counts(L, N)

Build pruning counts for completing a partial configuration to each symmetry sector.
"""
function _symmetry_suffix_counts(L::Int, N::Int)
    suffix = Vector{Dict{Tuple{Int,Int,Int},Int64}}(undef, L + 1)
    suffix[L + 1] = Dict{Tuple{Int,Int,Int},Int64}((0, 0, 0) => 1)

    for pos in (L - 1):-1:0
        next = suffix[pos + 2]
        current = Dict{Tuple{Int,Int,Int},Int64}()
        sign = iseven(pos) ? 1 : -1
        for ((used, P, A), count) in next
            for added in 0:(N - used)
                key = (used + added, mod(P + pos * added, L), A + sign * added)
                current[key] = get(current, key, 0) + count
            end
        end
        suffix[pos + 1] = current
    end
    return suffix
end

"""
    generate_symmetry_sector_states(L, N, target_P, target_A, target_count; verbose=true)

Enumerate exactly the symmetry sector `(target_P, target_A)` at fixed `L,N`.
States are packed in a `UInt128` with four bits per site.
"""
function generate_symmetry_sector_states(
    L::Int,
    N::Int,
    target_P::Int,
    target_A::Int,
    target_count::Int;
    verbose::Bool = true,
)
    _validate_krylov_size(L, N)
    suffix = _symmetry_suffix_counts(L, N)
    states = Set{UInt128}()
    sizehint!(states, target_count)
    generated = 0
    t0 = time()

    function rec(pos::Int, remaining::Int, P::Int, A::Int, packed::UInt128)
        if pos == L
            if remaining == 0 && P == target_P && A == target_A
                push!(states, packed)
                generated += 1
                if verbose && generated % 2_000_000 == 0
                    @printf("  generated=%d elapsed=%.1fs\n", generated, time() - t0)
                    flush(stdout)
                end
            end
            return
        end

        sign = iseven(pos) ? 1 : -1
        next_suffix = suffix[pos + 2]
        for value in 0:remaining
            next_remaining = remaining - value
            next_P = mod(P + pos * value, L)
            next_A = A + sign * value
            residual_key = (
                next_remaining,
                mod(target_P - next_P, L),
                target_A - next_A,
            )
            if get(next_suffix, residual_key, 0) > 0
                rec(
                    pos + 1,
                    next_remaining,
                    next_P,
                    next_A,
                    packed | (UInt128(value) << (KRYLOV_STATE_BITS * pos)),
                )
            end
        end
    end

    rec(0, N, 0, 0, UInt128(0))
    @assert length(states) == target_count
    return states
end

"""
    _push_krylov_neighbors!(neighbors, state, L)

Fill `neighbors` with states reachable by one `0110 <-> 1001` hop.
"""
function _push_krylov_neighbors!(neighbors::Vector{UInt128}, state::UInt128, L::Int)
    empty!(neighbors)
    for x in 0:(L - 1)
        i0 = x
        i1 = mod(x + 1, L)
        i2 = mod(x + 2, L)
        i3 = mod(x + 3, L)

        n0 = _krylov_getocc(state, i0)
        n1 = _krylov_getocc(state, i1)
        n2 = _krylov_getocc(state, i2)
        n3 = _krylov_getocc(state, i3)

        if n1 > 0 && n2 > 0
            new_state = state
            new_state = _krylov_addocc(new_state, i0, 1)
            new_state = _krylov_addocc(new_state, i1, -1)
            new_state = _krylov_addocc(new_state, i2, -1)
            new_state = _krylov_addocc(new_state, i3, 1)
            push!(neighbors, new_state)
        end

        if n0 > 0 && n3 > 0
            new_state = state
            new_state = _krylov_addocc(new_state, i0, -1)
            new_state = _krylov_addocc(new_state, i1, 1)
            new_state = _krylov_addocc(new_state, i2, 1)
            new_state = _krylov_addocc(new_state, i3, -1)
            push!(neighbors, new_state)
        end
    end
    return neighbors
end

"""
    analyze_krylov_sectors!(unvisited, L)

Explore Krylov sectors under the four-site hopping rule inside one symmetry
sector. The input set is mutated. The search stops once the largest Krylov
sector found is at least as large as the entire unvisited remainder, which
certifies that it is the true largest Krylov sector.
"""
function analyze_krylov_sectors!(unvisited::Set{UInt128}, L::Int)
    neighbors = UInt128[]
    component_sizes = Int[]
    largest = 0
    second = 0
    t0 = time()

    while !isempty(unvisited)
        start = first(unvisited)
        delete!(unvisited, start)
        queue = UInt128[start]
        head = 1

        while head <= length(queue)
            current = queue[head]
            head += 1
            _push_krylov_neighbors!(neighbors, current, L)
            for neighbor in neighbors
                if neighbor in unvisited
                    delete!(unvisited, neighbor)
                    push!(queue, neighbor)
                end
            end
        end

        component_size = length(queue)
        push!(component_sizes, component_size)
        if component_size > largest
            second = largest
            largest = component_size
        elseif component_size > second
            second = component_size
        end

        isempty(unvisited) && break

        if largest >= length(unvisited)
            return (
                krylov_sectors_seen = length(component_sizes),
                largest = largest,
                second_seen = second,
                remaining = length(unvisited),
                certified = true,
                exhausted = false,
                elapsed = time() - t0,
            )
        end
    end

    return (
        krylov_sectors_seen = length(component_sizes),
        largest = largest,
        second_seen = second,
        remaining = 0,
        certified = true,
        exhausted = true,
        elapsed = time() - t0,
    )
end

"""
    largest_krylov_sector(L, N; target_symmetry_sector=nothing, verbose=true)

Generate one symmetry sector and return its largest Krylov-sector ratio.
"""
function largest_krylov_sector(
    L::Int,
    N::Int;
    target_symmetry_sector = nothing,
    symmetry_sector_size = nothing,
    number_of_symmetry_sectors = nothing,
    verbose::Bool = true,
)
    if isnothing(target_symmetry_sector)
        (P, A), symmetry_sector_size, number_of_symmetry_sectors = largest_symmetry_sector(L, N)
    else
        P, A = target_symmetry_sector
        if isnothing(symmetry_sector_size) || isnothing(number_of_symmetry_sectors)
            counts = symmetry_sector_counts(L, N)
            symmetry_sector_size = Int(counts[(P, A)])
            number_of_symmetry_sectors = length(counts)
        else
            symmetry_sector_size = Int(symmetry_sector_size)
            number_of_symmetry_sectors = Int(number_of_symmetry_sectors)
        end
    end

    if verbose
        @printf(
            "\nSymmetry sector: L=%d N=%d rho=%.6f sector=(P=%d,A=%d) states=%d\n",
            L,
            N,
            N / L,
            P,
            A,
            symmetry_sector_size,
        )
    end

    states = generate_symmetry_sector_states(L, N, P, A, symmetry_sector_size; verbose)
    krylov_result = analyze_krylov_sectors!(states, L)
    ratio = krylov_result.largest / symmetry_sector_size

    if verbose
        @printf(
            "  largest=%d ratio=%.12f remaining=%d elapsed=%.1fs\n",
            krylov_result.largest,
            ratio,
            krylov_result.remaining,
            krylov_result.elapsed,
        )
    end

    return (;
        L,
        N,
        rho = N / L,
        P,
        A,
        symmetry_sector_size,
        number_of_symmetry_sectors,
        ratio,
        krylov_result...,
    )
end
