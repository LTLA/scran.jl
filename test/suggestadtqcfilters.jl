using SparseArrays
using scran

@testset "Testing the ADT QC filter suggestions" begin
    x = round.(abs.(sprand(Int32, 1000, 20, 0.2)) / 1e6)
    colsums = sum(x, dims=1)[1,:]
    mat = scran.initializesparsematrix(x)
    subsets = Dict{String,Any}("mito" => Vector{Int64}([5, 10, 20, 50]))
    metrics = scran.percelladtqcmetrics(mat; subsets = subsets)

    function reference_filter(
        detected::Vector{Int32}, 
        subsettotals::Dict{String,Vector{Float64}},
        thresholds_detected::Float64, 
        thresholds_subsettotals::Dict{String,Float64}
    )
        output = zeros(Bool, length(detected))

        for i in eachindex(detected)
            if detected[i] < thresholds_detected
                output[i] = 1
                continue
            end
        end

        for (k, v) in subsettotals
            threshold = thresholds_subsettotals[k]
            for i in eachindex(v)
                if v[i] > threshold
                    output[i] = 1
                end
            end
        end

        return output
    end

    @testset "no subsets, no blocks" begin
        subsettotals = Dict{String,Vector{Float64}}()
        suggested = scran.suggestadtqcfilters(metrics["detected"], subsettotals)

        @test suggested["thresholds"]["detected"] > 0
        @test length(suggested["thresholds"]["subsettotals"]) == 0

        ref = reference_filter(metrics["detected"], subsettotals, suggested["thresholds"]["detected"], Dict{String,Float64}())
        @test ref == suggested["filter"]
    end

    @testset "some subsets, no blocks" begin
        suggested = scran.suggestadtqcfilters(metrics["detected"], metrics["subsettotals"])

        @test suggested["thresholds"]["detected"] > 0
        @test length(suggested["thresholds"]["subsettotals"]) == length(subsets)

        ref = reference_filter(metrics["detected"], metrics["subsettotals"], suggested["thresholds"]["detected"], Dict{String,Float64}(suggested["thresholds"]["subsettotals"]))
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

        suggested = scran.suggestadtqcfilters(metrics["detected"], metrics["subsettotals"]; block = block)
        @test length(suggested["thresholds"]["detected"]) == 2
        @test length(suggested["thresholds"]["subsettotals"]) == length(subsets)

        prop_a = Dict{String, Vector{Float64}}()
        prop_b = Dict{String, Vector{Float64}}()
        for (k, v) in metrics["subsettotals"]
            prop_a[k] = v[slices["a"]]
            prop_b[k] = v[slices["b"]]
        end

        ref_a = scran.suggestadtqcfilters(metrics["detected"][slices["a"]], prop_a)
        @test suggested["thresholds"]["detected"]["a"] == ref_a["thresholds"]["detected"]
        @test suggested["thresholds"]["subsettotals"]["mito"]["a"] == ref_a["thresholds"]["subsettotals"]["mito"]

        ref_b = scran.suggestadtqcfilters(metrics["detected"][slices["b"]], prop_b)
        @test suggested["thresholds"]["detected"]["b"] == ref_b["thresholds"]["detected"]
        @test suggested["thresholds"]["subsettotals"]["mito"]["b"] == ref_b["thresholds"]["subsettotals"]["mito"]

        combined = zeros(Bool, size(mat, 2))
        combined[slices["a"]] = ref_a["filter"]
        combined[slices["b"]] = ref_b["filter"]
        @test combined == suggested["filter"]
    end
end
