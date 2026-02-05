current_dir = @__DIR__

include("../helpers/dmrg.jl")
include("../helpers/saving.jl")
include("../params.jl")

function main()
  N_vals = 10:10:80
  sigma_vals = [0.00, 0.00001, 0.0001, 0.001]
  accuracies = 10.0 .^ -(1:16)

  for N in N_vals
    for accuracy in accuracies
      for sigma in sigma_vals
        @info "Computing results for N=$N, σ=$sigma and acc=$accuracy"

        system_params = Dict{String,Any}(
          "J" => 1.0,
          "Δ" => 1.0,
          "μ" => 1.0,
          "NUM_SWEEPS" => 10,
          "MAX_BOND_DIM" => 1000,
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

        save_mps_with_params(save_path, psi, system_params, run_params)

        @info "Saved results for N=$N, σ=$sigma and acc=$accuracy"
      end
    end
  end

end

main()
