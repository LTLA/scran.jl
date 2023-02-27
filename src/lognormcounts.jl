include("utils.jl")

"""
    lognormcounts(x; 
        sizefactors = nothing, 
        block = nothing, 
        blockmethod = "lowest",
        center = true, 
        allowzeros = false, 
        numthreads = 1)

Compute log-normalized expression values from a count matrix in a `ScranMatrix` `x`.
This divides all counts by a cell-specific size factor before applying a log-2 transformation.
It returns a new `ScranMatrix` containing the log-normalized matrix.

`sizefactors` should be a vector of reals of length equal to the number of columns of `x`.
This should contain non-negative size factors to be applied to each cell in `x`.
If not provided, the size factor for each cell is derived from its library size.

`block` should be a vector of length equal to the number of columns of `x`.
This should specify the block assignment for each cell.
If not provided, all cells are assumed to belong to the same block.

`blockmethod` should describe how multi-block normalization should be performed.
If `"lowest"`, we downscale all batches to the coverage of the lowest batch.
If `"perblock"`, we scale each batch to a mean of 1.
This option is only relevant when `block` is provided.

If `center = true`, size factors are centered before use in normalization.
This should only be set to `false` if pre-centered `sizefactors` are supplied.

If `allowzeros = false`, an error is raised upon encountering size factors of zero.
If `true`, size factors of zero will automatically be be set to the smallest non-zero size factor after centering (or 1, if all size factors are zero).

The number of threads used is defined by `numthreads`.
This is only used when `sizefactors` is not supplied and we compute the library size for each cell.

# Examples
```jldoctest

julia> using SparseArrays

julia> x = abs.(sprand(Int32, 20, 10, 0.2))

julia> using scran

julia> mat = scran.initializesparsematrix(x)

julia> normed = scran.lognormcounts(mat)

julia> extractcolumn(normed, 1)
```
"""
function lognormcounts(x; sizefactors = nothing, block = nothing, blockmethod = "lowest", center = true, allowzeros = false, numthreads = 1)
    NC = size(x)[1]

    use_sf = false;
    if isnothing(sizefactors)
        sizefactors = Vector{Float64}()
    else
        if length(sizefactors) != NC
            throw(ErrorException("length of 'sizefactors' vector should be equal to number of columns in 'x'"))
        end

        use_sf = true
        if !(sizefactors isa Vector{Float64})
            sizefactors = Float64.(sizefactors)
        end
    end

    use_block, block_ids, _ = transform_factor(block, NC, "length of 'block' vector should be equal to number of columns in 'x'")

    return log_norm_counts(x, use_block, block_ids, blockmethod, use_sf, sizefactors, center, allowzeros, numthreads)
end
