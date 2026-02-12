current_dir = @__DIR__

include("../helpers/dmrg.jl")
include("../helpers/saving.jl")
include("../params.jl")

function main()
    N_vals = 40
    sigma_vals = [0.002]
    accuracies = [10^-16]
    repeats = 1

    num_tasks = repeats * length(N_vals) * length(sigma_vals) * length(accuracies)
    cur_task = 1
    for r in 1:repeats
        @info "Starting wave $r"
        for N in N_vals
            for accuracy in accuracies
                for sigma in sigma_vals
                    @info "Computing results for N=$N, σ=$sigma and acc=$accuracy"

                    system_params = Dict{String,Any}(
                    "J" => -1.0,
                    "Δ" => -1.0,
                    "μ" => 1.0,
                    "NUM_SWEEPS" => 30,
                    "MAX_BOND_DIM" => 21,
                    "ACC" => accuracy
                    )
                    run_params = Dict{String,Any}(
                    "N" => N,
                    "σ" => sigma
                    )

                    filename = extract_filename_from_system_params(system_params)
                    save_path = "data/storage/$(filename).hd5"

                    @time psi = find_ground_state_mps(run_params, system_params)
                    @info "Computation complete, saving..."

                    max_dim = maxlinkdim(psi)
                    @info "DMRG reported max link dimension: $max_dim"

                    save_mps_with_params(save_path, psi, system_params, run_params)

                    @info "Saved results for N=$N, σ=$sigma and acc=$accuracy"
                    percent_complete = (cur_task / num_tasks) * 100
                    @info "=== $percent_complete% complete ==="
                    cur_task += 1

                    # garbage collector has troublue with these large objects
                    # so we expiclity clear it here
                    GC.gc()
                    psi = nothing  # ... and large objects (mps)
                end
            end
        end
    end

end

main()
