using Pkg; Pkg.add(["HDF5", "ITensors"])
using HDF5, ITensors

function save_mps_with_params(filepath::String, ψ::MPS, params::Dict{String, Any})
    h5open(filepath, "w") do file
        write(file, "psi", ψ)

        g = create_group(file, "params")
        for (key, value) in params
            g[key] = value
        end # expand params
    end # hd5open
    println("Saved MPS to $filepath.")
end # function

function load_mps_with_params(filepath::String)
    params = Dict{String, Any}()

    h5open(filepath, "r") do file
        # Load the MPS
        ψ = read(file, "psi", MPS)

        # Load parameters
        g = file["params"]
        for key in keys(g)
            params[key] = read(g, key)
        end

        return ψ, params
    end
end