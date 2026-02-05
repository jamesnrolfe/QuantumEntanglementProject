using ITensors, ITensorMPS, HDF5, Plots, LinearAlgebra

include("../../data/params.jl")
include("../../data/helpers/saving.jl")   # must provide `find_runs_by_params`
include("../../data/helpers/dmrg.jl")

CURRENT_DIR = @__DIR__

base_system_params = Dict{String,Any}(
  "J" => 1.0,
  "Δ" => 1.0,
  "μ" => 1.0,
  "NUM_SWEEPS" => 10,
  "MAX_BOND_DIM" => 1000,
)

N = 30
sigma_vals = [0.00, 0.00001, 0.0001, 0.001]

# accuracies to scan
accuracies = 10.0 .^ -(1:16)

SIGMA_DATA = Dict{String,Dict{String,Any}}()

for sigma in sigma_vals
  SIGMA_DATA["$sigma"] = Dict{String,Any}()

  for acc in accuracies
    base_system_params["ACC"] = acc
    fn = extract_filename_from_system_params(base_system_params)
    path = joinpath(CURRENT_DIR, "../../data/storage/$(fn).hd5")

    query::Dict{String, Any} = Dict("N" => N, "σ" => sigma)

    matches = nothing
    try
      @info "Loading file at: $path with $query"
      matches = find_runs_by_params(path, query; load_psi=true, instance_selection=:latest)
    catch e
      @warn "Querying file failed" file=path exception=(e,catch_backtrace())
      matches = []
    end

    if isempty(matches)
      continue
    end

    # take the first match (should be only one per file with instance_selection=:latest)
    m = matches[1]
    if m.psi === nothing
      @warn "Found matching run but psi was not loaded for $(path) (acc=$(acc) σ=$(sigma))"
      continue
    end

    SIGMA_DATA["$sigma"]["$acc"] = m.psi
  end
end

ENERGY_DATA = Dict{String,Dict{String,Float64}}()

for sigma in sigma_vals
  ENERGY_DATA["$sigma"] = Dict{String,Float64}()

  @info "Collecting energy data for sigma=$sigma"

  sigma_data = get(SIGMA_DATA, "$sigma", nothing)
  if sigma_data === nothing
    continue
  end

  for acc in accuracies
    @info "Generating ground state energy for acc=$acc"
    mps = get(sigma_data, "$acc", nothing)
    if mps === nothing
      continue
    end
    sites = siteinds(mps)
    local A = generate_fully_connected_wam(N, sigma, base_system_params["μ"])
    local H = create_xxz_hamiltonian_mpo(N, A, base_system_params["J"], base_system_params["Δ"], sites)
    energy = inner(mps', H, mps)
    ENERGY_DATA["$sigma"]["$acc"] = energy
  end
end

@info "Plotting..."

p = plot(title="N=$N", xscale=:log10, xlabel="DMRG Cutoff Acc", ylabel="Ground state energy", legend=:top)

for sigma in sigma_vals
  acc_vals = Float64[]
  energy_vals = Float64[]

  energy_data = get(ENERGY_DATA, "$sigma", nothing)
  if energy_data === nothing
    continue
  end

  for acc in sort(accuracies)
    energy = get(energy_data, "$acc", nothing)
    if energy === nothing
      continue
    end

    push!(acc_vals, acc)
    push!(energy_vals, energy)
  end

  if !isempty(acc_vals)
    plot!(acc_vals, energy_vals,
      marker=:circle,
      markersize=6,
      linewidth=2,
      label="σ=$sigma")
  end
end

display(p)

out_path = joinpath(CURRENT_DIR, "energies_N$(N).png")
savefig(p, out_path)
@info "Saved figure to $out_path"
