
CURRENT_DIR = @__DIR__
include("../helpers/dmrg.jl")
include("../helpers/saving.jl")
include("../params.jl")

function main()
    N_vals = [40]
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
                    "NUM_SWEEPS" => 10,
                    "MAX_BOND_DIM" => 1000,
                    "ACC" => accuracy
                    )
                    clean_run_params = Dict{String,Any}(
                    "N" => N,
                    "σ" => 0.00
                    )
                    disorder_run_params = Dict{String,Any}(
                    "N" => N,
                    "σ" => sigma
                    )

                    filename = "perturb_$(extract_filename_from_system_params(system_params))"
                    save_path = "data/storage/$(filename).hd5"

                    @info "Creating adjacency matrices..."

                    A_cln = generate_fully_connected_wam(N, 0.00, system_params["μ"])
                    A_dis = generate_fully_connected_wam(N, sigma, system_params["μ"])

                    @info "Created adjacency matrices"

                    @info "Creating sites..."

                    sites = siteinds("S=1/2", N; conserve_qns=true)
                    psi_init = MPS(sites, [isodd(i) ? "Up" : "Dn" for i = 1:N])

                    @info "Created sites"

                    @info "Creating MPOs..."

                    H_cln = create_xxz_hamiltonian_mpo(N, A_cln, system_params["J"], system_params["Δ"], sites)
                    H_dis = create_xxz_hamiltonian_mpo(N, A_dis, system_params["J"], system_params["Δ"], sites)

                    @info "Created MPOs"

                    NUM_SWEEPS = system_params["NUM_SWEEPS"]
                    MAX_BOND_DIM = system_params["MAX_BOND_DIM"]
                    ACC = system_params["ACC"]

                    @info "Computing ground state..."

                    _, gs_cln = solve_xxz_hamiltonian_dmrg(H_cln, psi_init, NUM_SWEEPS, MAX_BOND_DIM, ACC)
                    _, gs_dis = solve_xxz_hamiltonian_dmrg(H_dis, psi_init, NUM_SWEEPS, MAX_BOND_DIM, ACC)

                    @info "Computation complete, saving..."

                    save_mps_with_params(save_path, gs_cln, system_params, clean_run_params)
                    save_mps_with_params(save_path, gs_dis, system_params, disorder_run_params)

                    @info "Saved results for N=$N, σ=$sigma and acc=$accuracy"
                    percent_complete = (cur_task / num_tasks) * 100
                    @info "=== $percent_complete% complete ==="
                    cur_task += 1

                    GC.gc()
                    gs_cln = nothing
                    gs_dis = nothing

                end
            end
        end
    end

end

main()