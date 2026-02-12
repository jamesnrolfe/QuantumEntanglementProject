function get_squared_sorted_schmidt_spectrum_from_mps(psi::MPS, cutoff=cutoff)
    N = length(psi)

    # MIDDLE bond
    middle_bond = N ÷ 2

    orthogonalize!(psi, middle_bond)
    link_idx = linkind(psi, middle_bond)
    _, S, _ = svd(psi[middle_bond], (siteind(psi, middle_bond), link_idx), cutoff=cutoff)

    schmidt_dim = dim(S, 1)

    schmidt_vals = []
    for n in 1:schmidt_dim
        push!(schmidt_vals, S[n, n]) # steal diagonals
    end

    # sort desc and square
    return sort(schmidt_vals .^ 2, rev=true)
end