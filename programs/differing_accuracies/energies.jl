using ITensors, ITensorMPS, HDF5, Plots, LinearAlgebra

include("../../data/params.jl")
include("../../data/helpers/saving.jl")
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

accuracies = 10.0 .^ -(1:16)

DATA = Dict{String,Any}()
for acc in accuracies
    base_system_params["ACC"] = acc
    fn = extract_filename_from_system_params(base_system_params)
    path = joinpath(CURRENT_DIR, "../../data/storage/$(fn).hd5")
    try
        system_params, runs = load_all_mps_from_file(path)
        DATA["$acc"] = runs
    catch
        @error "Failed to load data from file $path"
    end
end

SIGMA_DATA = Dict{String, Dict{String, Any}}()
for sigma in sigma_vals
    SIGMA_DATA["$sigma"] = Dict{String, Any}()
    for acc in accuracies
        runs = get(DATA, "$acc", nothing)
        if runs === nothing continue end
        for run in runs
            if run.params["σ"] ≈ sigma && run.params["N"] == N
                SIGMA_DATA["$sigma"]["$acc"] = run.psi
            end
        end

    end
end

ENERGY_DATA = Dict{String, Dict{String, Float64}}()



for sigma in sigma_vals
    ENERGY_DATA["$sigma"] = Dict{String, Float64}()

    sigma_data = get(SIGMA_DATA, "$sigma", nothing)
    if sigma_data === nothing continue end

    for acc in accuracies
        mps = get(sigma_data, "$acc", nothing)
        if mps === nothing continue end
        # reconstruct the Hamiltonian
        sites = siteinds(mps)
        local A = generate_fully_connected_wam(N, sigma, base_system_params["μ"])
        local H = create_xxz_hamiltonian_mpo(N, A, base_system_params["J"], base_system_params["Δ"], sites)
        # energy is expectation value of the hamiltonian on the state
        energy = inner(mps', H, mps)
        ENERGY_DATA["$sigma"]["$acc"] = energy
    end
end

p = plot(title="N=$N", xscale=:log10, xlabel="DMRG Cutoff Acc", ylabel="Ground state energy", legend=:topright)

for sigma in sigma_vals
    acc_vals = Float64[]
    energy_vals = Float64[]

    energy_data = get(ENERGY_DATA, "$sigma", nothing)
    if energy_data === nothing continue end

    for acc in sort(accuracies)
        energy = get(energy_data, "$acc", nothing)
        if energy === nothing continue end

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
