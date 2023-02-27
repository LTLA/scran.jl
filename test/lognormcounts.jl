using SparseArrays
using scran

@testset "Log-normalization tests" begin
    x = round.(abs.(sprand(Int32, 1000, 20, 0.2)) / 1e6)
    colsums = sum(x, dims=1)[1,:]
    mat = scran.initializesparsematrix(x)

    @testset "Basic usage" begin
        normed = scran.lognormcounts(mat)
        refsf = colsums / (sum(colsums) / length(colsums))
        @test isapprox(scran.extractrow(normed, 2), log2.(Vector{Float64}(x[2,:]) ./ refsf .+ 1))
        @test isapprox(scran.extractcolumn(normed, 9), log2.(Vector{Float64}(x[:,9]) / refsf[9] .+ 1))
    end

    @testset "Presupplied size factors" begin
        rawsf = Vector{Float64}(undef, size(x, 2))
        for i in eachindex(rawsf)
            rawsf[i] = rand() * 2
        end

        @testset "With automatic centering" begin
            normed = scran.lognormcounts(mat; sizefactors = rawsf)
            sf = rawsf / (sum(rawsf) / length(rawsf))
            @test isapprox(scran.extractrow(normed, 900), log2.(Vector{Float64}(x[900,:]) ./ sf .+ 1))
            @test isapprox(scran.extractcolumn(normed, 15), log2.(Vector{Float64}(x[:,15]) / sf[15] .+ 1))
        end

        @testset "Without automatic centering" begin
            normed = scran.lognormcounts(mat; sizefactors = rawsf, center = false)
            @test isapprox(scran.extractrow(normed, 800), log2.(Vector{Float64}(x[800,:]) ./ rawsf .+ 1))
            @test isapprox(scran.extractcolumn(normed, 12), log2.(Vector{Float64}(x[:,12]) / rawsf[12] .+ 1))
        end

        @testset "Fails if not of the right length" begin
            @test_throws "number of columns" scran.lognormcounts(mat; sizefactors = rawsf[1:10])
        end
    end

    @testset "Handle size factors of zero" begin
        sf = zeros(size(x, 2))
        @test_throws "size factors should be positive" scran.lognormcounts(mat; sizefactors = sf)

        # Doesn't bother normalizing in this case.
        normed = scran.lognormcounts(mat; sizefactors = sf, allowzeros = true)
        @test isapprox(scran.extractrow(normed, 100), log2.(scran.extractrow(mat, 100) .+ 1))
        @test isapprox(scran.extractcolumn(normed, 10), log2.(scran.extractcolumn(mat, 10) .+ 1))
    end

    @testset "Handle blocking factors" begin
        block = Vector{String}(undef, size(mat, 2))
        totals = Dict{String, Vector{Float64}}("A" => [0, 0], "B" => [0, 0])
        for i in eachindex(block)
            chosen = mod(i, 2) == 1 ? "A" : "B"
            block[i] = chosen
            current = totals[chosen]
            current[1] += colsums[i]
            current[2] += 1
        end

        for (k, v) in totals
            v[1] /= v[2]
        end

        @testset "Works with per-block" begin
            moresf = deepcopy(colsums)
            for i in eachindex(moresf)
                moresf[i] /= totals[block[i]][1]
            end

            normed = scran.lognormcounts(mat; block = block, blockmethod = "perblock")
            @test isapprox(scran.extractrow(normed, 200), log2.(scran.extractrow(mat, 200) ./ moresf .+ 1))
            @test isapprox(scran.extractcolumn(normed, 12), log2.(scran.extractcolumn(mat, 12) / moresf[12] .+ 1))
        end

        @testset "Works with per-block" begin
            lowest = min(totals["A"][1], totals["B"][1])
            moresf = deepcopy(colsums)
            for i in eachindex(moresf)
                moresf[i] /= lowest
            end

            normed = scran.lognormcounts(mat; block = block)
            @test isapprox(scran.extractrow(normed, 200), log2.(scran.extractrow(mat, 200) ./ moresf .+ 1))
            @test isapprox(scran.extractcolumn(normed, 12), log2.(scran.extractcolumn(mat, 12) / moresf[12] .+ 1))
        end

        @testset "Fail if block is not of the right size" begin
            @test_throws "number of columns" scran.lognormcounts(mat; block = block[1:10])
        end
    end
end
