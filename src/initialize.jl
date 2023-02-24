using SparseArrays

"""
    initializesparsematrix(x; nocopy = false, forceint = true)

Initialize a `ScranMatrix` from the matrix `x`.
Currently `x` is expected to be a `SparseMatrixCSC` instance from the SparseArrays package.

If `nocopy = true`, the non-zero values in `x` will not be copied into the output `ScranMatrix`.
This improves memory efficiency but assumes that `x` will not be mutated or garbage collected.

If `forceint = true`, non-integer values in `x` will be coerced to integer.
This is set to `true` for typical applications involving count matrices.

# Examples
```jldoctest

julia> using SparseArrays

julia> x = sprand(Int32, 8, 12, 0.2)

julia> using scran

julia> mat = scran.initializesparsematrix(sparsemat)

julia> mat isa scran.ScranMatrix
true
"""
function initializesparsematrix(x::SparseMatrixCSC{Tv, Ti}; nocopy = false, forceint = true) where {Tv, Ti}
    if Ti == Int32
        if Tv == Int32
            return initialize_from_memory_int32_int32(x.nzval, x.rowval, x.colptr, x.m, x.n, nocopy, false, forceint)
        elseif Tv == Int64
            return initialize_from_memory_int64_int32(x.nzval, x.rowval, x.colptr, x.m, x.n, nocopy, false, forceint)
        elseif Tv == Float32 
            return initialize_from_memory_float32_int32(x.nzval, x.rowval, x.colptr, x.m, x.n, nocopy, false, forceint)
        elseif Tv == Float64
            return initialize_from_memory_float64_int32(x.nzval, x.rowval, x.colptr, x.m, x.n, nocopy, false, forceint)
        end
    elseif Ti == Int64
        if Tv == Int32
            return initialize_from_memory_int32_int64(x.nzval, x.rowval, x.colptr, x.m, x.n, nocopy, false, forceint)
        elseif Tv == Int64
            return initialize_from_memory_int64_int64(x.nzval, x.rowval, x.colptr, x.m, x.n, nocopy, false, forceint)
        elseif Tv == Float32 
            return initialize_from_memory_float32_int64(x.nzval, x.rowval, x.colptr, x.m, x.n, nocopy, false, forceint)
        elseif Tv == Float64
            return initialize_from_memory_float64_int64(x.nzval, x.rowval, x.colptr, x.m, x.n, nocopy, false, forceint)
        end
    end

    throw(ErrorException("unsupported types for ScranMatrix initialization in SparseMatrixCSC 'x'"))
end

"""
    size(x)

Get the dimensions of a `ScranMatrix`, returning a 2-tuple containing the rows and columns.

# Examples
```jldoctest

julia> using SparseArrays

julia> x = sprand(Int64, 8, 12, 0.2)

julia> using scran

julia> mat = scran.initializesparsematrix(x)

julia> size(mat)
(8, 12)
```
"""
function Base.size(x::ScranMatrix)
    return (num_rows(x), num_columns(x))
end

"""
    extractrow(x, r)

Extract the values of the row `r` from the `ScranMatrix` `x`, returning a vector of double-precision floats of length equal to the number of columns.
This is mostly intended for use in testing; looping over this and calling a function on each row will be rather inefficient. 

# Examples
```jldoctest

julia> using SparseArrays

julia> x = sprand(Float64, 8, 12, 0.2) .* 10

julia> using scran

julia> mat = scran.initializesparsematrix(x; forceint = false)

julia> extractrow(mat, 1)
```
"""
function extractrow(x::ScranMatrix, r::Integer) 
    buffer = Vector{Float64}(undef, num_columns(x))
    row(x, r - 1, buffer)
    return buffer
end

"""
    extractcolumn(x, c)

Extract the values of the column `c` from the `ScranMatrix` `x`, returning a vector of double-precision floats of length equal to the number of rows.
This is mostly intended for use in testing; looping over this and calling a function on each column will be rather inefficient. 

# Examples
```jldoctest

julia> using SparseArrays

julia> x = sprand(Float64, 8, 12, 0.2) .* 10

julia> using scran

julia> mat = scran.initializesparsematrix(x)

julia> extractcolumn(mat, 1)
```
"""
function extractcolumn(x::ScranMatrix, c::Integer) 
    buffer = Vector{Float64}(undef, num_rows(x))
    column(x, c - 1, buffer)
    return buffer
end
