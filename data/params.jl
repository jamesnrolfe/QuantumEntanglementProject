function get_params()::Dict
    return Dict{String, Any}(
        "J" => 1.0,
        "Δ" => 1.0,
        "N_vals" => [6],
        "σ_vals" => [0.000, 0.001, 0.002],
        "μ" => 1.0,
        "NUM_SWEEPS" => 30,
        "MAX_BOND_DIM" => 1000,
        "ACC" => 1e-10
    )
end