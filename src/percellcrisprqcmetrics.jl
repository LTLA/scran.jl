"""
    percellcrisprqcmetrics(x::ScranMatrix; numthreads = 1)

Compute per-cell QC metrics from a `ScranMatrix` `x` of guide counts for CRISPR abundance data.
This returns a dictionary containing:

- `sums`: a vector of doubles containing the sum of counts within each cell.
- `detected`: a vector of integers containing the number of detected (i.e., non-zero) tags in each cell.
- `maxproportion`: a vector of doubles containing the proportion of counts in the most abundant guide for each cell.
- `maxindex`: a vector of integers specifying the row index of the most abundant guide for each cell.

The calculation may be parallelized by setting `numthreads` accordingly.

# Examples
```jldoctest

julia> using SparseArrays

julia> x = abs.(sprand(Int32, 50, 10, 0.2))

julia> using scran

julia> mat = scran.initializesparsematrix(x)

julia> scran.percellcrisprqcmetrics(mat)
```
"""
function percellcrisprqcmetrics(x::ScranMatrix; numthreads = 1)
    NR = size(x, 1)
    NC = size(x, 2)

    sums = Vector{Float64}(undef, NC)
    detected = Vector{Int32}(undef, NC)
    maxprop = Vector{Float64}(undef, NC)
    maxindex = Vector{Int32}(undef, NC)
    per_cell_crispr_qc_metrics(x, sums, detected, maxprop, maxindex, numthreads)

    maxindex .+= 1 # get back to 1-based indexing.

    return Dict{String, Any}(
        "sums" => sums, 
        "detected" => detected,
        "maxproportion" => maxprop,
        "maxindex" => maxindex
    )
end

