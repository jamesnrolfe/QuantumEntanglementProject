using Printf

function get_system_params()::Dict{String,Any}
    return Dict{String,Any}(
        "J" => 1.0,
        "Δ" => 1.0,
        "μ" => 1.0,
        "NUM_SWEEPS" => 10,
        "MAX_BOND_DIM" => 1000,
        "ACC" => 1e-10
    )
end

"""Get a standardised filename for a set of system parameters."""
function extract_filename_from_system_params(params::Dict{String,Any})::String
    string = ""
    for (key, value) in params
        if !isempty(string)
            string *= "_"
        end
        # format floats in scientific notation with 2 decimal places
        formatted_value = value isa AbstractFloat ? @sprintf("%.2e", value) : value
        string *= "$(key)=$(formatted_value)"
    end
    return string
end