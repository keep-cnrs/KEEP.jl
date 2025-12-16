module TestUtils

using LinearAlgebra: norm

export pretty_string, rel_neg_log_norm_diff, should_verbose
using Logging

function should_verbose()
    return Logging.min_enabled_level(current_logger()) <= Info
end

"""Pretty string of an object"""
function pretty_string(obj)
    io = IOBuffer()
    show(io, "text/plain", obj)
    String(take!(io))
end

"""arr has samples in columns and variables in rows. Returns the relative negative log norm difference between each pair of samples"""
function rel_neg_log_norm_diff(arr)
    norms = map(norm, eachslice(arr, dims=2))
    arr1 = reshape(arr, size(arr)..., 1)
    arr2 = permutedims(arr1, (1, 3, 2))
    diff_norms = map(norm, eachslice(arr1 .- arr2, dims=(2, 3)))
    @. -log10(diff_norms / max(norms, norms'))
end

end # module
