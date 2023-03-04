"""
    suggestadtqcfilters(
        detected::Vector{Int32}, 
        subsettotals::Dict{String,Vector{Float64}}; 
        block = nothing, 
        mindetecteddrop = 0.1,
        nummads = 3)

Suggest suitable filtering thresholds for the ADT-derived QC metrics.
This is based on the assumption that most cells are of high quality, such that outliers for the QC metrics (detected with a MAD-based algorithm) are potentially problematic and should be removed.
The function accepts the following arguments, typically generated by `percelladtqcmetrics`:

- `detected`, a vector containing the number of detected tags per cell.
- `subsettotals`, a dictionary containing the total count in each tag subset.

All vectors should have length equal to the number of cells.
This function then returns a dictionary containing the following keys:

- `thresholds`: a dictionary containing:
  - `detected`: the filtering threshold applied to the detected number of genes.
  - `subsettotals`: a dictionary where each key is a subset name and each value is the filtering threshold applied to the subset totals.
- `filter`: a vector of booleans where a truthy value indicates that the corresponding cell failed to pass at least one of the thresholds.
  Such cells are deemed to be of low-quality and should be discarded.

`block` should be a vector of length equal to the number of columns of `x`.
This should specify the block assignment for each cell.
Filtering thresholds are computed within each block, which avoids inflated MADs when there are large block-to-block differences.
If `block` is provided, all thresholds are instead dictionaries where each key is the block name and each value is the filtering threshold in that block.
If not provided, all cells are assumed to belong to the same block.

The number of MADs used to define the outlier thresholds is defined by `nummads`.
Larger values yield more relaxed filters.

For the number of detected tags, the MAD-based threshold is supplemented by `mindetecteddrop`.
This defines the minimum relative drop in the number of detected tags from the median across all cells;
by default, cells must have less than 90% of the median in order to fail the threshold.
The aim is to avoid overly aggressive filters when the MAD is low.

# Examples
```jldoctest

julia> using SparseArrays

julia> x = abs.(sprand(Int32, 50, 10, 0.2))

julia> using scran

julia> mat = scran.initializesparsematrix(x)

julia> sub = Dict{String,Vector{Int}}("IgG" => Vector{Int}([1,2,3,4,5]))

julia> metrics = scran.percelladtqcmetrics(mat; subsets = sub)

julia> scran.suggestadtqcfilters(metrics["detected"], metrics["subsettotals"])
```
"""
function suggestadtqcfilters(detected::Vector{Int32}, subsettotals::Dict{String,Vector{Float64}}; block = nothing, mindetectedrop=0.1, nummads = 3)
    ncells = length(detected)
    subnames = Vector{String}()
    props = Vector{Any}()
    for (k, v) in subsettotals 
        if length(v) != ncells
            throw(ErrorException("each entry of 'subsettotals' should be of the same length as 'sums'"))
        end
        push!(props, v)
        push!(subnames, k)
    end

    use_block, block_ids, block_levels = transform_factor(block, ncells, "length of 'block' vector should be equal to length of 'detected'")

    survivors = Vector{UInt8}(undef, ncells)
    thresholds = suggest_adt_qc_filters(detected, props, use_block, block_ids, survivors, mindetectedrop, nummads)

    prop = Dict{String, Any}()
    for i in eachindex(subnames)
        prop[subnames[i]] = transform_blocked_thresholds(thresholds_totals(thresholds, i - 1), use_block, block_levels)
    end

    output = Dict{String, Any}(
        "thresholds" => Dict{String, Any}(
            "detected" => transform_blocked_thresholds(thresholds_detected(thresholds), use_block, block_levels),
            "subsettotals" => prop
        ),
        "filter" => Vector{Bool}(survivors)
    )

    return output
end
