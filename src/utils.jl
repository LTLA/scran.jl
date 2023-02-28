function transform_factor(f, n::Integer, msg::String)
    if isnothing(f)
        return false, Vector{Int32}(), Vector{String}()
    end

    if length(f) != n
        throw(ErrorException(msg))
    end

    counter = 0
    levels = Dict{String, Int32}()
    ids = Vector{Int32}(undef, n)

    for i in eachindex(f)
        l = string(f[i])
        if haskey(levels, l)
            ids[i] = levels[l]
        else
            levels[l] = counter
            ids[i] = counter
            counter +=1
        end
    end

    olevels = Vector{String}(undef, length(levels))
    for (k, v) in levels
        olevels[v + 1] = k
    end

    return true, ids, olevels
end

function cast_to_logical(v::Vector{T}, NR::Integer) where T <: Bool
    if length(v) != NR
        throw(ErrorException("subsetting vector is longer than the dimension extent"))
    end
    return UInt8.(v)
end

function cast_to_logical(v::Vector{T}, NR::Integer) where T <: Integer
    output = zeros(UInt8, NR)
    for i in eachindex(v)
        if i <= 0 || i > NR
            throw(ErrorException("subsetting indices are out of range of the dimension extent"))
        end
        output[i] = 1
    end
    return output
end
