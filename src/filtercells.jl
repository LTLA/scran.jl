"""
    filtercells(x::ScranMatrix, discard::Vector{Bool})

Filters the columns of `x` to remove low-quality cells, marked as truthy values in `discard`.
This returns a `ScranMatrix` containing a subset of columns of `x`.

# Examples
```jldoctest

julia> using SparseArrays

julia> x = abs.(sprand(Int32, 50, 10, 0.2))

julia> using scran

julia> mat = scran.initializesparsematrix(x)

julia> metrics = scran.percellrnaqcmetrics(mat)

julia> suggested = scran.suggestrnaqcfilters(
           metrics["sums"], 
           metrics["detected"], 
           metrics["subsetproportions"]
       )

julia> filtered = scran.filtercells(mat, suggested["filter"])

julia> size(filtered, 2)
"""   
function filtercells(x::ScranMatrix, discard::Vector{Bool})
    if length(discard) != size(x, 2)
        throw(ErrorException("length of 'discard' should be equal to the number of columns in 'x'"))
    end
    return filter_cells(x, UInt8.(discard))
end
