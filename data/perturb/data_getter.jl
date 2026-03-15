using Pkg; Pkg.add(["ProgressMeter"])
using ProgressMeter


CURRENT_DIR = @__DIR__
include("../helpers/dmrg.jl")
include("../helpers/saving.jl")
include("../params.jl")

function main()
    # N_vals = vcat(2:1:19, 20:2:40, [50, 60])
    N_vals = [40]
    sigma_vals = [0.002]
    accuracies = [1e-16]
    repeats = 1

    J = Δ = -1.0
    μ = 1.0

    num_tasks = repeats * length(N_vals) * length(sigma_vals) * length(accuracies)

    prog = Progress(num_tasks; desc="DMRG Main Loop: ", dt=0.5)

    for r in 1:repeats
        for N in N_vals
            A_cln = generate_fully_connected_wam(N, 0.00, μ)
            sites = siteinds("S=1/2", N; conserve_qns=true)
            psi_init = MPS(sites, [isodd(i) ? "Up" : "Dn" for i = 1:N])
            H_cln = create_xxz_hamiltonian_mpo(N, A_cln, J, Δ, sites)

            for accuracy in accuracies
                system_params = Dict{String,Any}(
                    "J" => J,
                    "Δ" => Δ,
                    "μ" => μ,
                    "NUM_SWEEPS" => 25,
                    "MAX_BOND_DIM" => 2000,
                    "ACC" => accuracy
                )
                clean_run_params = Dict{String,Any}("N" => N, "σ" => 0.00)

                NUM_SWEEPS = system_params["NUM_SWEEPS"]
                MAX_BOND_DIM = system_params["MAX_BOND_DIM"]
                ACC = system_params["ACC"]

                _, gs_cln = solve_xxz_hamiltonian_dmrg(H_cln, psi_init, NUM_SWEEPS, MAX_BOND_DIM, ACC)

                for sigma in sigma_vals
                    disorder_run_params = Dict{String,Any}("N" => N, "σ" => sigma)

                    filename = "perturb_sigma_$(sigma)_$(extract_filename_from_system_params(system_params))"
                    save_path = "data/storage/$(filename).hd5"

                    A_dis = generate_fully_connected_wam(N, sigma, system_params["μ"])
                    H_dis = create_xxz_hamiltonian_mpo(N, A_dis, system_params["J"], system_params["Δ"], sites)

                    _, gs_dis = solve_xxz_hamiltonian_dmrg(H_dis, psi_init, NUM_SWEEPS, MAX_BOND_DIM, ACC)

                    save_mps_with_params(save_path, gs_cln, system_params, clean_run_params)
                    save_mps_with_params(save_path, gs_dis, system_params, disorder_run_params)

                    next!(prog; showvalues = [
                        (:Wave, r),
                        (:N, N),
                        (:Sigma, sigma),
                        (:Acc, accuracy),
                        (:Status, "Completed N=$N, σ=$sigma")
                    ])

                    gs_dis = nothing
                end
                gs_cln = nothing
                GC.gc()
            end
        end
    end
end

main()