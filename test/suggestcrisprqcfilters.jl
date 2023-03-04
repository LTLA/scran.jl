using SparseArrays
using scran

@testset "Testing the CRISPR QC filter suggestions" begin
    x = round.(abs.(sprand(Int32, 1000, 20, 0.2)) / 1e6)
    colsums = sum(x, dims=1)[1,:]
    mat = scran.initializesparsematrix(x)
    metrics = scran.percellcrisprqcmetrics(mat)

    @testset "no blocks" begin
        suggested = scran.suggestcrisprqcfilters(metrics["sums"], metrics["maxproportion"])

        @test suggested["thresholds"]["maxcount"] > 0

        maxcounts = metrics["sums"] .* metrics["maxproportion"]
        refabove = (maxcounts .< suggested["thresholds"]["maxcount"]) 
        @test refabove == suggested["filter"]
    end

    @testset "now with blocks" begin
        block = Vector{String}(undef, size(mat, 2))
        slices = Dict{String, Vector{Int32}}("a" => Vector{Int32}(), "b" => Vector{Int32}())
        for i in eachindex(block)
            chosen = mod(i, 2) == 1 ? "a" : "b"
            block[i] = chosen
            push!(slices[chosen], i)
        end

        suggested = scran.suggestcrisprqcfilters(metrics["sums"], metrics["maxproportion"]; block = block)
        @test length(suggested["thresholds"]["maxcount"]) == 2

        ref_a = scran.suggestcrisprqcfilters(metrics["sums"][slices["a"]], metrics["maxproportion"][slices["a"]])
        @test suggested["thresholds"]["maxcount"]["a"] == ref_a["thresholds"]["maxcount"]

        ref_b = scran.suggestcrisprqcfilters(metrics["sums"][slices["b"]], metrics["maxproportion"][slices["b"]])
        @test suggested["thresholds"]["maxcount"]["b"] == ref_b["thresholds"]["maxcount"]

        combined = zeros(Bool, size(mat, 2))
        combined[slices["a"]] = ref_a["filter"]
        combined[slices["b"]] = ref_b["filter"]
        @test combined == suggested["filter"]
    end
end
