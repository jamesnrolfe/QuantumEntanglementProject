include("../helpers/dmrg.jl")
include("../helpers/saving.jl")

function main()
  N = 30
  sigma_vals = [0.00, 0.0001, 0.001, 0.01]
  accuracies = 10.0 .^ -(1:16)

  for accuracy in accuracies
    for sigma in sigma_vals
      @info "Computing results for σ=$sigma and acc=$accuracy"
      system_params = Dict{String,Any}(
        "J" => 1.0,
        "Δ" => 1.0,
        "μ" => 1.0,
        "NUM_SWEEPS" => 10,
        "MAX_BOND_DIM" => 1000,
        "ACC" => accuracy
      )
      current_dir = @__DIR__
      save_path = "$current_dir/mps_data_at_acc_$accuracy.hd5"

      @time run_params, psi = compute_ground_state(N, sigma, system_params)
      @info "Computation complete, saving..."

      save_mps_with_params(save_path, psi, system_params, run_params)

      @info "Saved results for σ=$sigma and acc=$acc"
    end
  end

end

main()
