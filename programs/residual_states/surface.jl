using Pkg; Pkg.add(["ITensors", "ITensorMPS", "Kneedle", "Plots", "LsqFit", "Statistics", "LinearAlgebra", "ProgressMeter", "LaTeXStrings", "SavitzkyGolay", "PlotlyJS", "DelimitedFiles"])
using ITensors, ITensorMPS, LsqFit, Statistics, LinearAlgebra, ProgressMeter, LaTeXStrings, SavitzkyGolay, Kneedle, DelimitedFiles, Plots

include("../../data/params.jl")
include("../../data/helpers/saving.jl")
include("../schmidt_spectra/get_spectrum.jl")
include("../schmidt_spectra/deviation/deviation.jl")

function fetch_states(
    N::Int, σ::Float64, system_params::Dict{String, Any};
    max::Union{Nothing, Int}=nothing, type::Symbol=:latest
)::Tuple{Vector{MPS}, Vector{MPS}}
    CURRENT_DIR = @__DIR__
    FP = extract_filename_from_system_params(system_params)
    MPS_FP = joinpath(CURRENT_DIR, "../../data/storage/perturb_sigma_$(σ)_$(FP).hd5")
    matches_cln = find_runs_by_params(MPS_FP, Dict{String, Any}("N"=>N, "σ"=>0.00); load_psi=true, instance_selection=type)
    matches_dis = find_runs_by_params(MPS_FP, Dict{String, Any}("N"=>N, "σ"=>σ); load_psi=true, instance_selection=type)

    cln_mpss = MPS[]
    dis_mpss = MPS[]

    for (cln, dis) in zip(matches_cln, matches_dis)
        if (cln.psi === nothing || dis.psi === nothing)
            continue
        end

        push!(cln_mpss, cln.psi)
        push!(dis_mpss, dis.psi)

        if max !== nothing && length(cln_mpss) >= max
            break
        end
    end

    return cln_mpss, dis_mpss
end

system_params = Dict{String,Any}(
                    "J" => -1.0,
                    "Δ" => -1.0,
                    "μ" => 1.0,
                    "NUM_SWEEPS" => 25,
                    "MAX_BOND_DIM" => 2000,
                    "ACC" => 1e-16
                )

Ns = 2:1:80
σs = [5e-7, 6e-7, 7e-7, 8e-7, 9e-7, 1e-6, 5e-6]

get_β(ψ0, ψσ) = sqrt(1 - min(1, abs(inner(ψ0, ψσ))^2))


plotlyjs()

Z = fill(NaN, length(Ns), length(σs))

total_iters = length(σs) * length(Ns)
prog = Progress(total_iters; desc="Total Progress")

for (j, σ) in enumerate(σs)
    for (i, N) in enumerate(Ns)
        ψ0s, ψσs = fetch_states(N, σ, system_params; max=3, type=:all)

        if !isempty(ψ0s) && !isempty(ψσs)
            # Standardizing the inner product to avoid DomainErrors with sqrt
            loc_βs = [get_β(p0, ps) for (p0, ps) in zip(ψ0s, ψσs)]
            Z[i, j] = mean(loc_βs)
        end

        # Manual memory management for MPS objects
        ψ0s = nothing
        ψσs = nothing
        next!(prog)
    end

    # Apply Savitzky-Golay smoothing per sigma column
    column_data = Z[:, j]
    mask = .!isnan.(column_data)
    if sum(mask) >= 7
        try
            # window=7, poly=2
            smoothed = savitzky_golay(column_data[mask], 7, 2).y
            Z[mask, j] .= smoothed
        catch e
            @warn "SG smoothing failed for σ=$σ: $e"
        end
    end
    GC.gc()
end

p_surface = surface(Ns, σs, Z',
    title="Infidelity Surface: sqrt(1-ξ^2)",
    xlabel="System Size N",
    ylabel="Sigma σ",
    zlabel="Infidelity",
    colorbar=true,
    alpha=0.9,
    colormap=:viridis,
    size=(900, 700)
)

savefig(p_surface, joinpath(@__DIR__, "infidelity_surface.html"))
display(p_surface)