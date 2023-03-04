"""
    modelgenevar(x::ScranMatrix; block = nothing, span = 0.3, numthreads = 1)

Model the per-gene variances from a (typically log-normalized) expression matrix `x` by fitting a trend to account for the mean-variance relationship.
This returns a dictionary containing the following keys:

- `statistics`: a dictionary where each value is a double-precision vector of statistics across all genes.
  - `means`: the mean expression of each gene.
  - `variances`: the variance of each gene.
  - `fitted`: the fitted value of the mean-variance trend for each gene.
  - `residuals`: the residual from the trend for each gene.

`block` should be a vector of length equal to the number of columns of `x`.
This should specify the block assignment for each cell.
Variance modelling is then performed within each block to avoid inflated variances due to block-to-block differences.
If not provided, all cells are assumed to belong to the same block.

If `block` is supplied, the returned `statistics` are defined as the averages of the corresponding statistics across blocks.
The returned dictionary additionally contains `perblock`, a dictionary where each key is the name of a block 
and each value contains the per-block statistics in the same format as `statistics`.

# Examples
```jldoctest
julia> using SparseArrays

julia> x = abs.(sprand(Int32, 50, 10, 0.2))

julia> using scran

julia> mat = scran.initializesparsematrix(x)

julia> norm = scran.lognormcounts(mat)

julia> scran.modelgenevar(norm)
```
"""
function modelgenevar(x::ScranMatrix; block = nothing, span = 0.3, numthreads = 1)
    NR = size(x, 1)
    NC = size(x, 2)

    use_block, block_ids, block_levels = transform_factor(block, NC, "length of 'block' vector should be equal to number of columns in 'x'")
    res = model_gene_var(x, use_block, block_ids, span, numthreads)

    output = Dict{String, Any}()

    if !use_block
        output["statistics"] = Dict{String, Vector{Float64}}(
            "means" => copy(means(res, 0)),
            "variances" => copy(variances(res, 0)),
            "fitted" => copy(fitted(res, 0)),
            "residuals" => copy(residuals(res, 0))
        )
    else
        perblock = Dict{eltype(block_levels), Dict{String, Vector{Float64}}}()
        avemeans = zeros(Float64, NR)
        avevariances = zeros(Float64, NR)
        avefitted = zeros(Float64, NR)
        averesiduals = zeros(Float64, NR)

        for i in eachindex(block_levels)
            m = copy(means(res, i - 1))
            v = copy(variances(res, i - 1))
            f = copy(fitted(res, i - 1))
            r = copy(residuals(res, i - 1))

            perblock[block_levels[i]] = Dict{String, Vector{Float64}}("means" => m, "variances" => v, "fitted" => f, "residuals" => r)
            avemeans .+= m
            avevariances .+= v
            avefitted .+= f
            averesiduals .+= r
        end
        output["perblock"] = perblock

        avemeans ./= length(block_levels)
        avevariances ./= length(block_levels)
        avefitted ./= length(block_levels)
        averesiduals ./= length(block_levels)
        output["statistics"] = Dict{String, Vector{Float64}}("means" => avemeans, "variances" => avevariances, "fitted" => avefitted, "residuals" => averesiduals)
    end

    return output
end
