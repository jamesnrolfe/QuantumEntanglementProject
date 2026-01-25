using Distributed

const NUM_WORKERS = 8

# Add workers if needed
current_workers = nprocs() - 1
if current_workers < NUM_WORKERS
    addprocs(NUM_WORKERS - current_workers)
end

# Load necessary modules on all workers
@everywhere begin
    using Pkg
    Pkg.add(["ITensors", "ITensorMPS", "HDF5"])
    using ITensors, ITensorMPS, HDF5
    include("helpers/dmrg.jl")
end

# Load on main process
include("helpers/dmrg.jl")
include("helpers/saving.jl")
include("params.jl")

# Worker function: compute ground state MPS for given parameters
@everywhere function compute_ground_state(N::Int, σ::Float64, system_params::Dict{String,Any})
    println("Started task: N=$N, σ=$σ")
    run_params = Dict{String,Any}("N" => N, "σ" => σ)
    @time ψ_gs = find_ground_state_mps(run_params, system_params)
    println("Finished task: N=$N, σ=$σ")
    return (run_params, ψ_gs)
end

function main()
    save_path = "data/mps/default_params.hd5"
    N_vals = 10:5:80
    sigma_vals = [0.000, 0.001, 0.002]

    # Get system parameters
    system_params = get_system_params()

    # Create list of tasks
    tasks = [(N, σ) for N in N_vals for σ in sigma_vals]

    @info "Computing $(length(tasks)) ground states using $(nprocs()-1) workers..."

    # Distribute computation across workers
    results = pmap(tasks) do (N, σ)
        compute_ground_state(N, σ, system_params)
    end

    @info "Computation complete. Saving results to $save_path..."

    # Save all results to file
    for (run_params, ψ) in results
        save_mps_with_params(save_path, ψ, system_params, run_params)
        @info "Saved: N=$(run_params["N"]), σ=$(run_params["σ"])"
    end

    @info "All simulations complete and saved."

    # Clean up workers
    rmprocs(workers())
end

main()
