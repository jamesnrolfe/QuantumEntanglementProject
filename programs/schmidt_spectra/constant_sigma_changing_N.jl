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
σ = 0.001
N_vals = 20:20:80
system_params = Dict{String,Any}(
    "J" => -1.0,
    "Δ" => -1.0,
    "μ" => 1.0,
    "NUM_SWEEPS" => 30,
    "MAX_BOND_DIM" => 1000,
    "ACC" => 1e-10 # use a relatively low accuracy
)

FILENAME = extract_filename_from_system_params(system_params)
CURRENT_DIR = @__DIR__
MPS_FILEPATH = path = joinpath(CURRENT_DIR, "../../data/storage/$(FILENAME).hd5")
@info "Searching filepath '$MPS_FILEPATH'"


@info "Creating Schmidt spectra..."
# store schmidt spectra in a map relating to the system size
N_to_SS_map = Dict{Int, Vector{Float64}}()
for N in N_vals
    @info "Starting N=$N"
    run_params = Dict{String, Any}(
        "N" => N,
        "σ" => σ
    )
    matches = nothing
    try
        @info "Finding N=$N, σ=$σ"
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

    averaged_spectrum ./ length(all_spectra)

    # remove trailing zeros if any
    last_nonzero = findlast(x -> x > 0, averaged_spectrum)
    if last_nonzero !== nothing
        averaged_spectrum = averaged_spectrum[1:last_nonzero]
    end

    @info "Averaged $(length(all_spectra)) spectra for N=$N"

    sum_ss = sum(averaged_spectrum)
    @info " Average Schmidt Spectra sum: Σλ_i^2 = $sum_ss"

    N_to_SS_map[N] = averaged_spectrum
end

acc = system_params["ACC"]

p = plot(title="Schmidt Spectra, σ=$σ; cutoff=$acc", ylabel="Schmidt Values Squared", legend=:topright,
    yscale=:log,
    # xscale=:log
)

for N in N_vals
    if !haskey(N_to_SS_map, N)
        continue
    end

    spec = N_to_SS_map[N]

    plot!(p, 1:length(spec), spec, label="N=$N")
end

savefig(p, joinpath(CURRENT_DIR, "schmidt_spectra_constant_disorder_σ=$(σ)_acc=$acc.png"))
@info "Plot saved!"
