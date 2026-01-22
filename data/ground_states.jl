include("helpers/mps.jl")
include("helpers/saving.jl")
include("params.jl")

system_params::Dict{String, Any} = get_system_params()

params = Dict{String, Any}(
    "N" => 4,
    "σ" => 0.001
)
mps = find_ground_state_mps(params["N"], params["σ"], system_params)
save_mps_with_params("data/mps/test_save.hd5", mps, params)

psi, params = load_mps_with_params("data/mps/test_save.hd5")
println(psi)
println(params)