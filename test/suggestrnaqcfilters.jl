using SparseArrays
using scran

@testset "Testing the RNA QC filter suggestions" begin
    x = round.(abs.(sprand(Int32, 1000, 20, 0.2)) / 1e6)
    colsums = sum(x, dims=1)[1,:]
    mat = scran.initializesparsematrix(x)
    subsets = Dict{String,Any}("mito" => Vector{Int64}([5, 10, 20, 50]))
    metrics = scran.percellrnaqcmetrics(mat; subsets = subsets)

    function reference_filter(
        sums::Vector{Float64}, 
        detected::Vector{Int32}, 
        proportions::Dict{String,Vector{Float64}},
        thresholds_sums::Float64, 
        thresholds_detected::Float64, 
        thresholds_proportions::Dict{String,Float64}
    )
        output = zeros(Bool, length(sums))

        for i in eachindex(sums)
            if sums[i] < thresholds_sums
                output[i] = 1
                continue
            end
            if detected[i] < thresholds_detected
                output[i] = 1
                continue
            end
        end

        for (k, v) in proportions
            threshold = thresholds_proportions[k]
            for i in eachindex(v)
                if v[i] > threshold
                    output[i] = 1
                end
            end
        end

        return output
    end

    @testset "no subsets, no blocks" begin
        proportions = Dict{String,Vector{Float64}}()
        suggested = scran.suggestrnaqcfilters(metrics["sums"], metrics["detected"], proportions)

        @test suggested["thresholds"]["sums"] > 0
        @test suggested["thresholds"]["detected"] > 0
        @test length(suggested["thresholds"]["subsetproportions"]) == 0

        ref = reference_filter(
            metrics["sums"], 
            metrics["detected"], 
            proportions, 
            suggested["thresholds"]["sums"], 
            suggested["thresholds"]["detected"], 
            Dict{String,Float64}()
        )
        @test ref == suggested["filter"]
    end

    @testset "some subsets, no blocks" begin
        suggested = scran.suggestrnaqcfilters(metrics["sums"], metrics["detected"], metrics["subsetproportions"])

        @test suggested["thresholds"]["sums"] > 0
        @test suggested["thresholds"]["detected"] > 0
        @test length(suggested["thresholds"]["subsetproportions"]) == length(subsets)

        ref = reference_filter(
            metrics["sums"], 
            metrics["detected"], 
            metrics["subsetproportions"], 
            suggested["thresholds"]["sums"], 
            suggested["thresholds"]["detected"], 
            Dict{String,Float64}(suggested["thresholds"]["subsetproportions"])
        )
        @test ref == suggested["filter"]
    end

    @testset "now with blocks" begin
        block = Vector{String}(undef, size(mat, 2))
        slices = Dict{String, Vector{Int32}}("a" => Vector{Int32}(), "b" => Vector{Int32}())
        for i in eachindex(block)
            chosen = mod(i, 2) == 1 ? "a" : "b"
            block[i] = chosen
            push!(slices[chosen], i)
        end

        suggested = scran.suggestrnaqcfilters(metrics["sums"], metrics["detected"], metrics["subsetproportions"]; block = block)
        @test length(suggested["thresholds"]["sums"]) == 2
        @test length(suggested["thresholds"]["detected"]) == 2
        @test length(suggested["thresholds"]["subsetproportions"]) == length(subsets)

        prop_a = Dict{String, Vector{Float64}}()
        prop_b = Dict{String, Vector{Float64}}()
        for (k, v) in metrics["subsetproportions"]
            prop_a[k] = v[slices["a"]]
            prop_b[k] = v[slices["b"]]
        end

        ref_a = scran.suggestrnaqcfilters(metrics["sums"][slices["a"]], metrics["detected"][slices["a"]], prop_a)
        @test suggested["thresholds"]["sums"]["a"] == ref_a["thresholds"]["sums"]
        @test suggested["thresholds"]["detected"]["a"] == ref_a["thresholds"]["detected"]
        @test suggested["thresholds"]["subsetproportions"]["mito"]["a"] == ref_a["thresholds"]["subsetproportions"]["mito"]

        ref_b = scran.suggestrnaqcfilters(metrics["sums"][slices["b"]], metrics["detected"][slices["b"]], prop_b)
        @test suggested["thresholds"]["sums"]["b"] == ref_b["thresholds"]["sums"]
        @test suggested["thresholds"]["detected"]["b"] == ref_b["thresholds"]["detected"]
        @test suggested["thresholds"]["subsetproportions"]["mito"]["b"] == ref_b["thresholds"]["subsetproportions"]["mito"]

        combined = zeros(Bool, size(mat, 2))
        combined[slices["a"]] = ref_a["filter"]
        combined[slices["b"]] = ref_b["filter"]
        @test combined == suggested["filter"]
    end
end
