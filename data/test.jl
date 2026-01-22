include("helpers/dmrg.jl")
include("helpers/saving.jl")
include("params.jl")

system_params::Dict{String,Any} = get_system_params()

run_params = Dict{String,Any}(
    "N" => 4,
    "σ" => 0.001
)
mps = find_ground_state_mps(run_params, system_params)
save_mps_with_params(
    "data/mps/test_save.hd5",
    mps,
    system_params,
    run_params
)

psi, sys, run, = load_mps_with_params("data/mps/test_save.hd5")
println(psi)
println(sys)
println(run)
