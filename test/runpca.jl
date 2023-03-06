using SparseArrays
using scran

@testset "PCA tests" begin
    x = round.(abs.(sprand(Int32, 500, 200, 0.2)) / 1e6)
    mat = scran.initializesparsematrix(x)
    norm = scran.lognormcounts(mat)
    
    @testset "Unblocked, no subset" begin
        res = runpca(norm; components = 50)
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
        sub = Vector{Int}()
        for i in 1:size(x, 1)
            if i % 2 == 1 
                push!(sub, i)
            end
        end

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
end

