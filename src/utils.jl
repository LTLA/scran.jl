function transform_factor(f::Nothing, n::Integer, msg::String)
    if isnothing(f)
        return false, Vector{Int32}(), Vector{String}()
    end
end

function transform_factor(f::Vector{T}, n::Integer, msg::String) where T
    if length(f) != n
        throw(ErrorException(msg))
    end

    keys = unique(f)
    sort!(keys)

    # Short cut if it's already properly formatted.
    if T == Int32 && length(keys) && keys[1] == 0 && keys[length(keys)] == length(keys) - 1
        return f
    end

    levels = Dict{T, Int32}()
    for i in eachindex(keys)
        levels[keys[i]] = i - 1 # Get to 0-based indices
    end

    ids = Vector{Int32}(undef, n)
    for i in eachindex(f)
        ids[i] = levels[f[i]]
    end

    return true, ids, keys
end

function cast_to_logical(v::Vector{T}, NR::Integer) where T <: Bool
    if length(v) != NR
        throw(ErrorException("subsetting vector is longer than the dimension extent"))
    end
    return UInt8.(v)
end

function cast_to_logical(v::Vector{T}, NR::Integer) where T <: Integer
    output = zeros(UInt8, NR)
    for i in v
        if i <= 0 || i > NR
            throw(ErrorException("subsetting indices are out of range of the dimension extent"))
        end
        output[i] = 1
    end
    return output
end
