using Pkg;
Pkg.add(["HDF5", "ITensors"]);
using HDF5, ITensors
using Dates, Logging

"""
Save an MPS with its system parameters and its run parameters.

Storage Layout
- /system_params        (group)     : holds system-wide parameters
- /runs/<run_id>/psi    (dataset)   : the MPS for this run
- /runs/<run_id>/params (dataset)   : run-specific parameters
"""
function save_mps_with_params(
    filepath::String,
    ψ::MPS,
    system_params::Dict{String,Any},
    run_params::Dict{String,Any};
    run_id::Union{String,Nothing}=nothing,
    param_safety::Bool=true
)
    mode = isfile(filepath) ? "r+" : "w"
    h5open(filepath, mode) do file
        # system params - create only if not present
        if !haskey(file, "system_params")
            system_group = create_group(file, "system_params")
            for (k, v) in system_params
                try
                    system_group[k] = v
                catch e
                    @warn "Could not write system_params[$k] as HDF5; saving as string" exception = (e, catch_backtrace())
                    system_group[k] = string(v)
                end # try
            end # k,v loop
        else
            # if system params exist, do not overwrite - warn if new keys present
            existing = file["system_params"]
            for (k, v) in system_params
                if !haskey(existing, k)
                    @warn "system_parrams in file is missing key $k. This function will not overwrite."
                    if param_safety
                        error("param_safety is enabled, program exiting.
                            If you wish to force exit the program, disable this parameter (NOT RECOMMENDED).")
                    else
                        @warn "Data will be written, but its system data might be incorred. Consider writing to a different file."
                    end
                end # haskey
            end # k,v loop
        end # haskey

        # prepare runs group
        runs_group = haskey(file, "runs") ? file["runs"] : create_group(file, "runs")

        # determine run_id
        if run_id === nothing
            # Find numeric run names and pick next integer id
            existing_names = collect(keys(runs_group))
            nums = Int[]
            for nm in existing_names
                n = tryparse(Int, string(nm))
                if n !== nothing
                    push!(nums, n)
                end
            end
            next_num = isempty(nums) ? 1 : maximum(nums) + 1
            run_id = string(next_num)
        else
            # ensure run_id is a string and unique; if it collides, append suffixes
            run_id = string(run_id)
            if haskey(runs_group, run_id)
                base = run_id
                counter = 1
                while haskey(runs_group, "$(base)_$(counter)")
                    counter += 1
                end
                run_id = "$(base)_$(counter)"
            end
        end
        # create run group
        run_group = create_group(runs_group, run_id)

        # write the mps in this run_group
        try
            write(run_group, "psi", ψ)
        catch e
            @error "Failed to write MPS to HDF5" exception = (e, catch_backtrace())
            rethrow(e)
        end

        # write run specific params under runs/<run_id>/params
        run_params_group = create_group(run_group, "params")
        for (k, v) in run_params
            try
                run_params_group[k] = v
            catch e
                @warn "Could not write run_params[$k] as HDF5; saving as string" exception = (e, catch_backtrace())
            end
        end

        @info "Saved run '$run_id' to $filepath."
    end
end # function


"""
Load an MPS together with system and run params.

If `run_id` is ommitted, uses the most recently created run.

Returns ψ::MPS, system_params::Dict{String, Any}, run_params::Dict{String, Any}.
"""
function load_mps_with_params(
    filepath::String;
    run_id::Union{String,Nothing}=nothing
)::Tuple{MPS,Dict{String,Any},Dict{String,Any}}

    system_params = Dict{String,Any}()
    run_params = Dict{String,Any}()

    result = h5open(filepath, "r") do file
        if haskey(file, "system_params")
            system_group = file["system_params"]
            for key in keys(system_group)
                try
                    system_params[string(key)] = read(system_group, key)
                catch
                    # read as string as fallback
                    system_params[string(key)] = string(system_group[key])
                end
            end
        else
            @warn "No `system_params` group found in $filepath."
        end

        # choose run
        if !haskey(file, "runs")
            error("No 'runs' group found in $filepath")
        end

        runs_group = file["runs"]
        run_names = collect(keys(runs_group))
        if isempty(run_names)
            error("No runs found in 'runs' group of $filepath")
        end

        if run_id === nothing
            # Prefer numeric run ids: choose the run with largest integer name.
            nums_map = Dict{Int,String}()
            for nm in run_names
                n = tryparse(Int, string(nm))
                if n !== nothing
                    nums_map[n] = string(nm)
                end
            end
            if !isempty(nums_map)
                run_id = nums_map[maximum(keys(nums_map))]
            else
                # fallback to lexicographic last if there are no purely numeric run names
                run_id = sort(run_names)[end]
            end
        elseif !in(string(run_id), run_names)
            error("Requested run_id `$(run_id)` not found in file. Available runs: $(run_names).")
        else
            # normalize to the stored string name
            run_id = string(run_id)
        end

        run_group = runs_group[run_id]

        # read MPS
        psi = nothing
        try
            # Prefer reading with a target type so HDF5.jl can convert to MPS
            psi = read(run_group, "psi", MPS)
        catch e
            # If conversion fails, read raw object and print diagnostics
            @warn "Direct read into MPS failed; attempting raw read for diagnostics."
            raw = try
                read(run_group, "psi")
            catch inner
                @error "Raw read also failed" exception = (inner, catch_backtrace())
                rethrow(inner)
            end
            @info "Raw read type: $(typeof(raw))"
            if isa(raw, Dict)
                @info "Raw read keys: $(collect(keys(raw)))"
            end
            rethrow(e)
        end

        # read run params
        if haskey(run_group, "params")
            run_group_params = run_group["params"]
            for key in keys(run_group_params)
                try
                    run_params[string(key)] = read(run_group_params, key)
                catch
                    run_params[string(key)] = string(run_group_params[key])
                end
            end
        end

        (psi, system_params, run_params)
    end

    return result
end

"""
Return the system parameters in a specific file.
"""
function get_system_params_of_file(filepath::String)::Union{Dict{String,Any},Nothing}
    system_params = Dict{String,Any}()

    h5open(filepath, "r") do file
        if haskey(file, "system_params")
            system_group = file["system_params"]
            for key in keys(system_group)
                try
                    system_params[string(key)] = read(system_group, key)
                catch
                    # read as string as fallback
                    system_params[string(key)] = string(system_group[key])
                end
            end
        else
            @warn "No `system_params` group found in $filepath."
            return nothing
        end
    end

    return system_params
end

"""
Load all the run data from a specific file.

Will return all runs associated with the set of system parameters defining the file.
"""
function load_all_mps_from_file(filepath::String)
    # Try to read system params (may return nothing)
    system_params = get_system_params_of_file(filepath)
    if system_params === nothing
        system_params = Dict{String,Any}()
    end

    runs = Vector{NamedTuple{(:run_id, :psi, :params),Tuple{String,MPS,Dict{String,Any}}}}()

    h5open(filepath, "r") do file
        if !haskey(file, "runs")
            @warn "No 'runs' group found in $filepath"
            return (system_params, runs)
        end

        runs_group = file["runs"]
        run_names = collect(keys(runs_group))
        if isempty(run_names)
            @warn "No runs found in 'runs' group of $filepath"
            return (system_params, runs)
        end

        # Determine ordering: prefer numeric ids in ascending order if any are purely numeric,
        # otherwise use lexicographic ordering.
        nums_map = Dict{Int,String}()
        for nm in run_names
            n = tryparse(Int, string(nm))
            if n !== nothing
                nums_map[n] = string(nm)
            end
        end

        ordered_run_names = if !isempty(nums_map)
            [nums_map[k] for k in sort(collect(keys(nums_map)))]
        else
            sort(string.(run_names))
        end

        for rn in ordered_run_names
            run_group = runs_group[rn]

            # Read MPS (prefer direct conversion to MPS, otherwise provide diagnostics and rethrow)
            psi = nothing
            try
                psi = read(run_group, "psi", MPS)
            catch e
                @warn "Direct read into MPS failed for run $rn; attempting raw read for diagnostics."
                raw = try
                    read(run_group, "psi")
                catch inner
                    @error "Raw read also failed for run $rn" exception = (inner, catch_backtrace())
                    rethrow(inner)
                end
                @info "Raw read for run $rn type: $(typeof(raw))"
                if isa(raw, Dict)
                    @info "Raw read keys for run $rn: $(collect(keys(raw)))"
                end
                rethrow(e)
            end

            # Read run params if present
            run_params = Dict{String,Any}()
            if haskey(run_group, "params")
                run_params_group = run_group["params"]
                for key in keys(run_params_group)
                    try
                        run_params[string(key)] = read(run_params_group, key)
                    catch
                        run_params[string(key)] = string(run_params_group[key])
                    end
                end
            end

            push!(runs, (run_id=string(rn), psi=psi, params=run_params))
        end
    end

    return (system_params, runs)
end
