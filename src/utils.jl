function transform_factor(f, n::Integer, msg::String)
    if isnothing(f)
        return false, Vector{Int32}(), Vector{String}()
    end

    if length(f) != n
        throw(ErrorException(msg))
    end

    counter = 0
    levels = Dict{String, Int32}
    ids = Vector{Int32}(undef, n)

    for i in eachindex(f)
        l = string(f[i])
        if haskey(levels, l)
            ids[i] = levels[l]
        else
            levels[l] = counter
            counter +=1
        end
    end

    olevels = Vector{String}(undef, length(levels))
    for (k, v) in levels
        olevels[v] = k
    end

    return true, ids, olevels
end

