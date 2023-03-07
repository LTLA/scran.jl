using SparseArrays
using scran

function define_subset_indices(x::scran.ScranMatrix)
    sub = Vector{Int}()
    for i in 1:size(x, 1)
        if i % 2 == 1 
            push!(sub, i)
        end
    end
    return sub
end

function define_block(x::scran.ScranMatrix)
    block = Vector{String}(undef, size(x, 2))
    for i in eachindex(block)
        block[i] = (i % 2 == 1 ? "A" : "B")
    end
    return block
end

@testset "PCA tests" begin
    x = round.(abs.(sprand(Int32, 500, 200, 0.2)) / 1e6)
    mat = scran.initializesparsematrix(x)
    norm = scran.lognormcounts(mat)
    res = runpca(norm; components = 50)

    @testset "Unblocked, no subset" begin
        @test size(res["principalcomponents"], 1) == 50
        @test size(res["principalcomponents"], 2) == size(x, 2)
        @test length(res["varianceexplained"]) == 50
        @test issorted(-res["varianceexplained"]) 

        res = runpca(norm; components = 50, rotation = true)
        @test size(res["rotation"], 1) == size(x, 1)
        @test size(res["rotation"], 2) == 50

        # Parallel execution gives the same result.
        parallel = runpca(norm, components = 50, numthreads = 2)
        @test parallel["principalcomponents"] == res["principalcomponents"]

        # Scaling has some kind of effect.
        scaled = runpca(norm, components = 50, scale = true)
        @test scaled["principalcomponents"] != res["principalcomponents"]
    end

    @testset "Unblocked with subset" begin
        sub = define_subset_indices(norm)
        res = runpca(norm; subset = sub, rotation = true) 
        @test size(res["rotation"], 1) == length(sub)
        @test size(res["rotation"], 2) == 50

        # Same with a boolean.
        sub = Vector{Bool}(undef, size(x, 1))
        for i in eachindex(sub)
            sub[i] = i % 2 == 1 
        end
        res2 = runpca(norm; subset = sub, rotation = true) 
        @test res["rotation"] == res2["rotation"]
    end

    @testset "Blocked with regression" begin
        block = define_block(norm)
        bres = runpca(norm; components = 50, block = block)
        @test size(bres["principalcomponents"]) == size(res["principalcomponents"])
        @test bres["principalcomponents"] != res["principalcomponents"]
        @test bres["varianceexplained"] != res["varianceexplained"]

        # Combined with feature subsetting.
        sub = define_subset_indices(norm)
        bres = runpca(norm; subset = sub, block = block, rotation = true) 
        @test size(bres["rotation"], 1) == length(sub)
    end
end

