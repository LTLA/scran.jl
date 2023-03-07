"""
    runpca(
        x::ScranMatrix; 
        components = 50, 
        subset = nothing, 
        block = nothing,
        scale = false, 
        numthreads = 1, 
        rotation = false
    )

Perform a principal components analysis on a log-normalized expression matrix `x`, returning the top `components` components that explain the most variance.
This returns a dictionary containing the following fields:

- `principalcomponents`: a matrix where columns are cells (corresponding to columns of `x`) and rows are principal components.
  The number of reported components is no greater than `compoennts`.
- `variancexplained`: a vector of length equal to the number of reported components,
  containing the proportion of variance explained by each component.
- `rotation`: a matrix where columns are components and rows are genes, containing the rotation vector corresponding to each component.
  If `subset = nothing`, the number of genes is equal to the number of rows of `x`;
  otherwise, if `subset` is supplied, the number of genes is equal to the size of the subset.
  Only reported if `rotation = true`.

Users can perform the PCA on a subset of features (e.g., highly variable genes) by setting `subset`.
This can either be a vector of booleans of length equal to the number of rows in `x`, where a true value indicates that the corresponding feature should be used;
or a vector of integers containing the row indices of the features of interest in `x`.
(In the latter case, it is recommended that the vector be strictly increasing for easier matching to rows of the `rotation` matrix.)
If `nothing`, all features are used.

If `block` is supplied and `blockmethod = "block"`, a PCA is performed after regressing out the blocking factor from each row of `x`.
This effectively performs the pCA on the residuals, which eliminates any linear block-to-block differences. 

If `scale = true`, genes are scaled to unit variance prior to PCA.
If `block` is supplied and `blockmethod = "block"`, this scaling is done on the variance of the residuals.

Users can set `numthreads` to enable parallelization.

# Examples
```jldoctest
julia> using SparseArrays

julia> x = abs.(sprand(Int32, 500, 100, 0.2))

julia> using scran

julia> mat = scran.initializesparsematrix(x)

julia> norm = scran.lognormcounts(mat)

julia> runpca(norm)
```
"""
function runpca(x::ScranMatrix; components = 50, subset = nothing, block = nothing, blockmethod = "block", scale = false, numthreads = 1, rotation = false)
    use_subset = !isnothing(subset)
    feats = Vector{UInt8}()
    if use_subset
        feats = cast_to_logical(subset, size(x, 1))
    end

    if isnothing(block)
        res = run_pca(x, components, use_subset, feats, scale, numthreads)
    else
        use_block, block_ids, _ = transform_factor(block, size(x, 2), "length of 'block' vector should be equal to number of columns in 'x'")
        if blockmethod == "block"
            res = run_blocked_pca(x, components, use_subset, feats, block_ids, scale, numthreads) 
        end
    end

    pcs = Array{Float64}(undef, num_dims(res), num_obs(res))
    pcs[:,:] = principal_components(res)
    output = Dict{String, Any}(
        "principalcomponents" => pcs, 
        "varianceexplained" => variance_explained(res) ./ total_variance(res)
    )

    if rotation
        rot = Array{Float64}(undef, num_genes(res), num_dims(res))
        rot[:,:] = scran.rotation(res)
        output["rotation"] = rot
    end

    return output
end
