using Pkg
Pkg.add(["LinearAlgebra", "ITensors", "ITensorMPS", "Statistics"])

using LinearAlgebra
using ITensors
using ITensorMPS
using Statistics

# PARAMS
include("../params.jl")
params = get_params()

J = params["J"]
Δ = params["Δ"]
N_vals = params["N_vals"]
σ_vals = params["σ_vals"]
μ = params["μ"]
NUM_SWEEPS = params["NUM_SWEEPS"]
MAX_BOND_DIM = params["MAX_BOND_DIM"]
ACC = params["ACC"]

function generate_fully_connected_wam(N::Int, σ::Float64, μ::Float64)::Matrix{Float64}
    """
    Generate a randomly weighted adjacency matrix for a fully connected graph of N nodes.
    """
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

function create_xxz_hamiltonian_mpo(N::Int, A::Matrix{Float64}, J::Float64, Δ::Float64, sites::Vector{Index{Vector{Pair{QN, Int64}}}})::MPO
    """Create the XXZ Hamiltonian as an MPO given an adjacency matrix."""
    mpo = OpSum()
    for i = 1:N-1
        for j = i+1:N
            weight = A[i, j]
            if weight != 0.0
                # if the weight is zero, then we shouldn't add a connection
                # XX and YY terms: S+S- + S-S+ = 2(SxSx + SySy)
                # So to get J(SxSx + SySy), we need J/2 * (S+S- + S-S+)
                mpo += weight * J/2, "S+", i, "S-", j
                mpo += weight * J/2, "S-", i, "S+", j
                # ZZ term
                mpo += weight * J * Δ, "Sz", i, "Sz", j
            end # weight conditional
          end # j loop
    end # i loop
    H = MPO(mpo, sites)
    return H
end # function

function solve_xxz_hamiltonian_dmrg(H::MPO, ψ0::MPS, num_sweeps::Int, bond_dim::Int, cutoff::Float64)::tuple{Float64, MPS}
    """Apply the DMRG to a Hamiltonian."""
    local sweeps = Sweeps(num_sweeps)
    setmaxdim!(sweeps, bond_dim)
    setcutoff!(sweeps, cutoff)
    E, ψ = dmrg(H, ψ0, sweeps; outputlevel = 0) # output level 0 to make it quieter
    return E, ψ
end # function

function create_mps(N::Int; conserve_qns::Bool=true)::Tuple{MPS, Vector{Index{Vector{Pair{QN, Int64}}}}}
    """Create a random MPS for a spin-1/2 graph of size N."""
    # create a site set for a spin-1/2 system
    sites::Vector{Index{Vector{Pair{QN, Int64}}}} = siteinds("S=1/2", N; conserve_qns=conserve_qns)
    # create a random MPS
    return MPS(sites, [isodd(i) ? "Up" : "Dn" for i = 1:N]), sites
end # function

function find_ground_state_mps(N::Int, σ::Float64)::MPS
    """Helper function to find the ground state MPS for a given N and σ."""
    ψ, sites = create_mps(N)
    A = generate_fully_connected_wam(N, σ, μ)
    H = create_xxz_hamiltonian_mpo(N, A, J, Δ, sites)
    _, ψ_gs = (H, ψ, NUM_SWEEPS, MAX_BOND_DIM, ACC)
    return ψ_gs
end # function

println(find_ground_state_mps(6, 0.001))