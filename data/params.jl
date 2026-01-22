function get_system_params()::Dict{String, Any}
    return Dict{String, Any}(
        "J" => 1.0,
        "Δ" => 1.0,
        "μ" => 1.0,
        "NUM_SWEEPS" => 30,
        "MAX_BOND_DIM" => 1000,
        "ACC" => 1e-10
    )
end