"""
    percellrnaqcmetrics(x::ScranMatrix; subsets = nothing, numthreads = 1)

Compute per-cell QC metrics from a `ScranMatrix` `x` of gene counts for RNA expression data.
This returns a dictionary containing:

- `sums`, a vector of doubles containing the sum of counts within each cell.
- `detected`, a vector of integers containing the number of detected (i.e., non-zero) counts in each cell.
- `subsetproportions`, an inner dictionary where each key is the name of a feature subset.
  Each value is a vector of doubles containing the proportion of counts in the subset.

If `subsets` is provided, it should be a dictionary where each key is a string containing the name of a feature subset.
If the value is a vector of booleans, it should be of length equal to the number of rows in `x`, where truthy values represent membership of the corresponding row in the subset.
If the value is a vector of integers, those integers are treated as row indices of `x` specifying the features that are members of the subset.

The calculation may be parallelized by setting `numthreads` accordingly.

# Examples
```jldoctest

julia> using SparseArrays

julia> x = abs.(sprand(Int32, 50, 10, 0.2))

julia> using scran

julia> mat = scran.initializesparsematrix(x)

julia> scran.percellrnaqcmetrics(mat)

julia> sub = Dict{String,Vector{Int}}("mito" => Vector{Int}([1,2,3,4,5]))

julia> percellrnaqcmetrics(mat; subsets = sub)
```
"""
function percellrnaqcmetrics(x::ScranMatrix; subsets = nothing, numthreads = 1)
    NR = size(x, 1)
    NC = size(x, 2)

    subs = Vector{Any}()
    subout = Vector{Any}()
    subnames = Vector{String}()
    if !isnothing(subsets)
        for (k, v) in subsets
            push!(subs, cast_to_logical(v, NR))
            push!(subnames, String(k))
            push!(subout, Vector{Float64}(undef, NC))
        end
    end

    sums = Vector{Float64}(undef, NC)
    detected = Vector{Int32}(undef, NC)
    per_cell_rna_qc_metrics(x, subs, sums, detected, subout, numthreads)

    output = Dict{String, Any}("sums" => sums, "detected" => detected)
    subdict = Dict{String, Vector{Float64}}()
    for i in eachindex(subnames)
        subdict[subnames[i]] = subout[i]
    end
    output["subsetproportions"] = subdict

    return output
end

