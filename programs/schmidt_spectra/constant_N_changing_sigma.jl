# the idea here is that if we keep the level of disorder constant
# we can see how the tail scales with system size.
# we expect to see behaviour where the tail gets longer and flatter,
# essentially holding the same amount of total "area"
# or, we might see something completely different

using ITensors, ITensorMPS, Plots
include("../../data/params.jl")
include("../../data/helpers/saving.jl")
include("./get_spectrum.jl")

@info "Setting paramters..."
# now want to load in the mps states we want
σ_vals = [0.00, 0.0001, 0.001, 0.002]
N = 40
system_params = Dict{String,Any}(
    "J" => -1.0,
    "Δ" => -1.0,
    "μ" => 1.0,
    "NUM_SWEEPS" => 30,
    "MAX_BOND_DIM" => 1000,
    "ACC" => 1e-16
)

FILENAME = extract_filename_from_system_params(system_params)
CURRENT_DIR = @__DIR__
MPS_FILEPATH = path = joinpath(CURRENT_DIR, "../../data/storage/$(FILENAME).hd5")
@info "Searching filepath '$MPS_FILEPATH'"


@info "Creating Schmidt spectra..."
# store schmidt spectra in a map relating to the system size
sigma_to_SS_map = Dict{Float64, Vector{Float64}}()
for sigma in σ_vals
    @info "Starting sigma=$sigma"
    run_params = Dict{String, Any}(
        "N" => N,
        "σ" => sigma
    )
    matches = nothing
    try
        @info "Finding N=$N, σ=$sigma"
        matches = find_runs_by_params(MPS_FILEPATH, run_params; load_psi=true, instance_selection=:all)
    catch e
        @warn "Could not find MPS data for N=$N, σ=$σ"
        matches = []
    end

    if isempty(matches) continue end



    num_runs = min(3, length(matches))
    @info "Using latest $num_runs runs for averaging."

    all_spectra = []
    for i in 1:num_runs
        m = matches[i]
        if m.psi === nothing
            @warn "Could not find psi for run $i"
            continue end

        @info "Run $i: Calculating Schmidt Spectra"
        schmidt_spectra = get_squared_sorted_schmidt_spectrum_from_mps(m.psi)
        push!(all_spectra, schmidt_spectra)
    end

    if isempty(all_spectra)
        @warn "No valid spectra for N=$N"
        continue end

    max_len = maximum(length.(all_spectra))

    averaged_spectrum = zeros(max_len)
    for spec in all_spectra
        for (i,val) in enumerate(spec)
            averaged_spectrum[i] += val
        end
    end

    @debug "Len all spectra $length(all_spectra)"
    averaged_spectrum ./ length(all_spectra)

    # remove trailing zeros if any
    last_nonzero = findlast(x -> x > 0, averaged_spectrum)
    if last_nonzero !== nothing
        averaged_spectrum = averaged_spectrum[1:last_nonzero]
    end

    @info "Averaged $(length(all_spectra)) spectra for sigma=$sigma"

    sum_ss = sum(averaged_spectrum)
    @info " Average Schmidt Spectra sum: Σλ_i^2 = $sum_ss"

    sigma_to_SS_map[sigma] = averaged_spectrum
end

acc = system_params["ACC"]

p = plot(title="Schmidt Spectra, N=$N; cutoff=$acc", ylabel="Schmidt Values Squared", legend=:topright,
    yscale=:log,
    # xscale=:log
)

for sigma in σ_vals
    if !haskey(sigma_to_SS_map, sigma)
        continue
    end

    spec = sigma_to_SS_map[sigma]

    plot!(p, 1:length(spec), spec, label="σ=$sigma")
end

savefig(p, joinpath(CURRENT_DIR, "schmidt_spectra_constant_size_N=$(N)_acc=$(acc).png"))
@info "Plot saved!"
