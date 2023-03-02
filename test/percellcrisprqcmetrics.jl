using SparseArrays
using scran

@testset "Per-cell CRISPR QC metrics testing" begin
    x = round.(abs.(sprand(Int32, 1000, 20, 0.2)) / 1e6)
    colsums = sum(x, dims=1)[1,:]
    mat = scran.initializesparsematrix(x)

    metrics = scran.percellcrisprqcmetrics(mat)
    @test isapprox(metrics["sums"], colsums)

    detected = sum(x .> 0, dims=1)[1,:]
    @test detected == metrics["detected"]

    refmaxprop = Vector{Float64}(undef, size(x, 2))
    refmaxprop2 = Vector{Float64}(undef, size(x, 2))
    counter = 1
    for col in eachcol(x)
        best = argmax(col)
        refmaxprop[counter] = col[best] / colsums[counter]
        refmaxprop2[counter] = col[metrics["maxindex"][counter]] / colsums[counter]
        counter += 1
    end
    @test isapprox(metrics["maxproportion"], refmaxprop)
    @test refmaxprop == refmaxprop2 # don't compare indices directly in case of ties.
end
