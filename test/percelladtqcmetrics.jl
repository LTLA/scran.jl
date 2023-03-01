using SparseArrays
using scran

@testset "Per-cell ADT QC metrics testing" begin
    x = round.(abs.(sprand(Int32, 1000, 20, 0.2)) / 1e6)
    colsums = sum(x, dims=1)[1,:]
    mat = scran.initializesparsematrix(x)

    @testset "no subsets" begin
        metrics = scran.percelladtqcmetrics(mat)
        @test isapprox(metrics["sums"], colsums)

        detected = sum(x .> 0, dims=1)[1,:]
        @test detected == metrics["detected"]

        @test length(metrics["subsettotals"]) == 0
    end

    @testset "some subsets" begin
        sub1 = Vector{Int64}([5, 10, 20, 50])
        sub2 = Vector{Bool}(undef, size(x, 1))
        for i in eachindex(sub2)
            sub2[i] = rand() < 0.05
        end

        subsets = Dict{String,Any}("SUB1" => sub1, "SUB2" => sub2)
        metrics = scran.percelladtqcmetrics(mat; subsets = subsets)
        @test length(metrics["subsettotals"]) == 2

        totals1 = sum(x[sub1,:], dims=1)[1,:]
        totals2 = sum(x[sub2,:], dims=1)[1,:]
        @test isapprox(metrics["subsettotals"]["SUB1"], totals1)
        @test isapprox(metrics["subsettotals"]["SUB2"], totals2)
    end
end
