using Kneedle
using SavitzkyGolay

function find_deviation_idx(schmidt_spectrum::Vector{Float64}; window_size=11, poly_order=2)::Int64

    log_spec = log10.(max.(schmidt_spectrum, 1e-18))
    sg = savitzky_golay(log_spec, window_size, poly_order) # smooth to get better results
    smoothed_log_y = sg.y
    x = 1:length(smoothed_log_y)
    res = kneedle(x, smoothed_log_y)

    deviation_index = knees(res)
    return deviation_index[1]
end