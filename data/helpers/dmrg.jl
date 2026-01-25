using Pkg
Pkg.add(["LinearAlgebra", "ITensors", "ITensorMPS", "Statistics"])

using LinearAlgebra
using ITensors
using ITensorMPS
using Statistics

const _sites_cache = Dict{Tuple{Int,Bool},Any}()

"""
    Generate a randomly weighted adjacency matrix for a fully connected graph of N nodes.
    """
function generate_fully_connected_wam(N::Int, σ::Float64, μ::Float64)::Matrix{Float64}
    A = zeros(Float64, N, N)
    for i in 1:N
        for j in (i+1):N
            weight = μ + σ * randn() # weight from ND with mean μ and std σ
            A[i, j] = weight
            A[j, i] = weight
        end # j loop
    end # i loop
    return A
end # function

"""Create the XXZ Hamiltonian as an MPO given an adjacency matrix."""
function create_xxz_hamiltonian_mpo(N::Int, A::Matrix{Float64}, J::Float64, Δ::Float64, sites::Vector{Index{Vector{Pair{QN,Int64}}}})::MPO
    mpo = OpSum()
    for i = 1:N-1
        for j = i+1:N
            weight = A[i, j]
            if weight != 0.0
                # if the weight is zero, then we shouldn't add a connection
                # XX and YY terms: S+S- + S-S+ = 2(SxSx + SySy)
                # So to get J(SxSx + SySy), we need J/2 * (S+S- + S-S+)
                mpo += weight * J / 2, "S+", i, "S-", j
                mpo += weight * J / 2, "S-", i, "S+", j
                # ZZ term
                mpo += weight * J * Δ, "Sz", i, "Sz", j
            end # weight conditional
        end # j loop
    end # i loop
    H = MPO(mpo, sites)
    return H
end # function

"""Apply the DMRG to a Hamiltonian."""
function solve_xxz_hamiltonian_dmrg(H::MPO, ψ0::MPS, num_sweeps::Int, bond_dim::Int, cutoff::Float64)::Tuple{Float64,MPS}
    local sweeps = Sweeps(num_sweeps)
    setmaxdim!(sweeps, bond_dim)
    setcutoff!(sweeps, cutoff)
    E, ψ = dmrg(H, ψ0, sweeps; outputlevel=0) # output level 0 to make it quieter
    return E, ψ
end # function

"""Create a random MPS for a spin-1/2 graph of size N."""
function create_mps(N::Int; conserve_qns::Bool=true)::Tuple{MPS,Vector{Index{Vector{Pair{QN,Int64}}}}}
    # create a site set for a spin-1/2 system

    key = (N, conserve_qns)
    sites = get(_sites_cache, key, nothing)
    if sites === nothing
        sites = siteinds("S=1/2", N; conserve_qns=conserve_qns)
        _sites_cache[key] = sites
    end

    # create a random MPS
    return MPS(sites, [isodd(i) ? "Up" : "Dn" for i = 1:N]), sites
end # function

"""Helper function to find the ground state MPS for a given N and σ."""
function find_ground_state_mps(run_params::Dict{String,Any}, system_params::Dict{String,Any})::MPS
    N = run_params["N"]
    σ = run_params["σ"]

    J = system_params["J"]
    Δ = system_params["Δ"]
    μ = system_params["μ"]
    NUM_SWEEPS = system_params["NUM_SWEEPS"]
    MAX_BOND_DIM = system_params["MAX_BOND_DIM"]
    ACC = system_params["ACC"]

    ψ, sites = create_mps(N)
    A = generate_fully_connected_wam(N, σ, μ)
    H = create_xxz_hamiltonian_mpo(N, A, J, Δ, sites)
    _, ψ_gs = solve_xxz_hamiltonian_dmrg(H, ψ, NUM_SWEEPS, MAX_BOND_DIM, ACC)
    return ψ_gs
end # function
