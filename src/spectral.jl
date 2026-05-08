function _representative(state::Vector{Int}, lump_final::Bool)
    L = length(state)
    if lump_final && is_final_binary(state)
        return zeros(Int, L)
    end

    best = copy(state)
    for shift in 0:(L - 1)
        shifted = circshift(state, shift)
        shifted < best && (best = shifted)
        reflected = reverse(shifted)
        reflected < best && (best = reflected)
    end
    return best
end

"""
    build_generator(L, gamma, lambda; max_total=L, lump_final=false)

Construct the symmetry-reduced Markov generator for the finite pairwise model
on a periodic lattice. Returns `(W, representatives, multiplicities)`, where
columns of `W` sum to zero and each representative stands for all translations
and reflections of a state.

Set `lump_final=true` to merge all exact binary final states into one absorbing
representative.
"""
function build_generator(
    L::Int,
    gamma::Float64,
    lambda::Float64;
    max_total::Int = L,
    lump_final::Bool = false,
)
    L >= 4 || throw(ArgumentError("Pairwise dynamics require L >= 4."))
    lump_final && L <= 7 &&
        throw(ArgumentError("Exact binary final-state patterns assume L > 7."))

    representatives = Vector{Vector{Int}}()
    multiplicities = Int[]
    rep_to_idx = Dict{Vector{Int},Int}()

    function generate_states!(current_state, site, current_sum)
        if site > L
            rep = _representative(current_state, lump_final)
            idx = get(rep_to_idx, rep, 0)
            if idx == 0
                push!(representatives, copy(rep))
                push!(multiplicities, 1)
                rep_to_idx[rep] = length(representatives)
            else
                multiplicities[idx] += 1
            end
            return
        end

        for occ in 0:(max_total - current_sum)
            current_state[site] = occ
            generate_states!(current_state, site + 1, current_sum + occ)
        end
    end

    generate_states!(zeros(Int, L), 1, 0)
    n_states = length(representatives)
    row_idx = Int[]
    col_idx = Int[]
    values = Float64[]

    function add_transition!(source_idx, new_state, rate)
        new_rep = _representative(new_state, lump_final)
        dest_idx = get(rep_to_idx, new_rep, 0)
        dest_idx == 0 && return

        push!(row_idx, dest_idx)
        push!(col_idx, source_idx)
        push!(values, rate)

        push!(row_idx, source_idx)
        push!(col_idx, source_idx)
        push!(values, -rate)
    end

    for (source_idx, state) in enumerate(representatives)
        lump_final && is_final_binary(state) && continue

        for i in 1:L
            if state[i] >= 2
                new_state = copy(state)
                new_state[i] -= 2
                rate = gamma * state[i] * (state[i] - 1) / 2.0
                add_transition!(source_idx, new_state, rate)
            end

            i1 = i
            i2 = i % L + 1
            i3 = (i + 1) % L + 1
            i4 = (i + 2) % L + 1

            if state[i2] >= 1 && state[i3] >= 1
                new_state = copy(state)
                new_state[i1] += 1
                new_state[i2] -= 1
                new_state[i3] -= 1
                new_state[i4] += 1
                add_transition!(source_idx, new_state, lambda * state[i2] * state[i3])
            end

            if state[i1] >= 1 && state[i4] >= 1
                new_state = copy(state)
                new_state[i1] -= 1
                new_state[i2] += 1
                new_state[i3] += 1
                new_state[i4] -= 1
                add_transition!(source_idx, new_state, lambda * state[i1] * state[i4])
            end
        end
    end

    W = sparse(row_idx, col_idx, values, n_states, n_states)
    return W, representatives, multiplicities
end

"""
    find_spectral_gap(W; num_eigenvalues=3, max_iter=3000, krylov_dim=30)

Compute the leading right and left eigenvectors of a generator using ARPACK.
Returns eigenvalues sorted by decreasing real part, followed by right and left
eigenvector matrices in the same order.
"""
function find_spectral_gap(
    W;
    num_eigenvalues::Int = 3,
    max_iter::Int = 3000,
    krylov_dim::Int = 30,
)
    n = size(W, 1)
    nev = min(num_eigenvalues, n - 1)
    ncv = min(krylov_dim, n - 1)
    ncv = max(ncv, min(n - 1, 2 * nev + 1))

    vals_r, vecs_r = eigs(
        W;
        nev,
        which = :LR,
        maxiter = max_iter,
        ncv,
    )
    vals_l, vecs_l = eigs(
        copy(W');
        nev,
        which = :LR,
        maxiter = max_iter,
        ncv,
    )

    perm_r = sortperm(real.(vals_r), rev = true)
    perm_l = sortperm(real.(vals_l), rev = true)
    return vals_r[perm_r], vecs_r[:, perm_r], vecs_l[:, perm_l]
end

"""
    calculate_coefficient(rho0, k, right_vectors, left_vectors, representatives, multiplicities, L)

Calculate the initial-state overlap coefficient of mode `k` for a Bernoulli
binary initial ensemble with density `rho0`.
"""
function calculate_coefficient(
    rho0,
    k,
    right_vectors,
    left_vectors,
    representatives,
    multiplicities,
    L,
)
    eig_idx = k + 1
    v_k = right_vectors[:, eig_idx]
    u_k = left_vectors[:, eig_idx]

    overlap_numerator = 0.0
    for idx in eachindex(representatives)
        if maximum(representatives[idx]) < 2
            n = sum(representatives[idx])
            p0_state = idx == 1 ? 0.0 : multiplicities[idx] * rho0^n * (1 - rho0)^(L - n)
            overlap_numerator += u_k[idx] * p0_state
        end
    end

    overlap = overlap_numerator / sum(u_k .* v_k)
    return -real(overlap * v_k[1])
end

"""
    full_gamma_coefficients(L, gamma, lambda, rho0; nmax=2)

Build a lumped generator and return the first `nmax` decay rates and overlap
coefficients for the Bernoulli initial density `rho0`.
"""
function full_gamma_coefficients(L, gamma, lambda, rho0; nmax = 2)
    W, representatives, multiplicities = build_generator(
        L,
        Float64(gamma),
        Float64(lambda);
        lump_final = true,
    )
    vals, right_vectors, left_vectors = find_spectral_gap(W; num_eigenvalues = nmax + 1)
    gammas = -real.(vals[2:end])
    coefficients = [
        calculate_coefficient(
            rho0,
            k,
            right_vectors,
            left_vectors,
            representatives,
            multiplicities,
            L,
        ) for k in 1:nmax
    ]
    return gammas, coefficients
end
