"""
    groupedsizefactors(
        x::ScranMatrix, 
        group; 
        center = true, 
        priorcount = 10, 
        reference = nothing, 
        numthreads=1
    )

Compute size factors to remove composition biases between groups of related cells in a `ScranMatrix` `x`.
This constructs pseudo-bulk profiles for each group, and then computes median-based size factors to normalize those profiles between groups,
taking advantage of the higher pseudo-bulk counts to obtain a stable median.

`group` should be a vector of length equal to the number of columns in `x`.
Each entry should specify the group that the corresponding cell is assigned to.
Cells in the same group should be biologically similar (e.g., derived from clustering on data in `x`) to remove composition biases between biological subpopulations.

If `center = true`, the output size factors are centered.

`priorcount` is a positive number specifying the prior count to add to the pseudo-bulk profiles prior to normalization.
Specifically, a scaled version of the reference profile is added to each profile, 
where the scaling is defined so that the increase in the total count from this addition is equal to `priorcount` on average.
Larger values improve stability for near-zero size factors at the cost of accurate removal of composition biases.

`reference` specifies the identity of the group in `group` to be used as a reference profile.
This should ideally be a group with many cells, such that the pseudo-bulk profile has few zeros;
and be not too dissimilar in its expression profile when compared to other groups.
If `reference = nothing`, a high-coverage group is automatically chosen.

The number of threads used is defined by `numthreads`.

# Examples
```jldoctest

julia> using SparseArrays

julia> x = abs.(sprand(Int32, 50, 10, 0.2))

julia> using scran

julia> mat = scran.initializesparsematrix(x)

julia> group = Vector{String}(undef, size(x, 2))

julia> levels = [ "A", "B", "C" ]

julia> for i in 1:size(x, 2)
           group[i] = levels[i % length(levels) + 1]
       end

julia> factors = scran.groupedsizefactors(mat, group)
"""
function groupedsizefactors(x::ScranMatrix, group; center = true, priorcount=10, reference = nothing, numthreads=1)
    NC = size(x, 2)
    use_group, group_ids, group_levels = transform_factor(group, NC, "length of 'group' vector should be equal to number of columns in 'x'")

    ref = -1
    if !isnothing(reference)
        for i in 1:length(group_levels)
            if group_levels[i]== reference
                ref = i - 1
            end
        end
        if ref == -1
            throw(ErrorException("reference level '" + String(reference) + " not found in 'group'"))
        end
    end

    output = Vector{Float64}(undef, NC)
    grouped_size_factors(
        x,
        group_ids,
        output,
        center, 
        priorcount, 
        ref,
        numthreads
    )

    return output
end
