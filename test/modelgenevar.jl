using SparseArrays
using scran

@testset "Variance modelling tests" begin
    x = round.(abs.(sprand(Int32, 1000, 20, 0.2)) / 1e6)
    NR = size(x, 1)
    NC = size(x, 2)

    mat = scran.initializesparsematrix(x)
    norm = scran.lognormcounts(mat)

    @testset "Basic usage" begin
        res = scran.modelgenevar(norm)
        stats = res["statistics"]

        refmeans = Vector{Float64}(undef, NR)
        refvars = Vector{Float64}(undef, NR)
        for i = 1:NR
            stuff = extractrow(norm, i)
            refmeans[i] = sum(stuff) / NC
            refvars[i] = sum((stuff .- refmeans[i]).^2) / (NC - 1)
        end

        @test isapprox(res["statistics"]["means"], refmeans)
        @test isapprox(res["statistics"]["variances"], refvars)
        @test isapprox(res["statistics"]["variances"], res["statistics"]["fitted"] .+ res["statistics"]["residuals"])
    end

    @testset "Plus blocks" begin
        block = Vector{String}(undef, NC)
        slices = Dict{String, Vector{Int32}}("A" => Vector{Int32}(), "B" => Vector{Int32}())
        for i in eachindex(block)
            chosen = mod(i, 2) == 1 ? "A" : "B"
            block[i] = chosen
            push!(slices[chosen], i)
        end

        res = scran.modelgenevar(norm, block = block)
        stats = res["statistics"]
        perblock = res["perblock"]

        @test isapprox(stats["means"], (perblock["A"]["means"] .+ perblock["B"]["means"]) ./ 2)
        @test isapprox(stats["variances"], (perblock["A"]["variances"] .+ perblock["B"]["variances"]) ./ 2)
        @test isapprox(stats["fitted"], (perblock["A"]["fitted"] .+ perblock["B"]["fitted"]) ./ 2)
        @test isapprox(stats["residuals"], (perblock["A"]["residuals"] .+ perblock["B"]["residuals"]) ./ 2)

        for (k, v) in slices
            refmeans = Vector{Float64}(undef, NR)
            refvars = Vector{Float64}(undef, NR)
            for i = 1:NR
                stuff = extractrow(norm, i)[v]
                refmeans[i] = sum(stuff) / length(v)
                refvars[i] = sum((stuff .- refmeans[i]).^2) / (length(v) - 1)
            end

            @test isapprox(perblock[k]["means"], refmeans)
            @test isapprox(perblock[k]["variances"], refvars)
            @test isapprox(perblock[k]["variances"], perblock[k]["fitted"] .+ perblock[k]["residuals"])
        end
    end
end
