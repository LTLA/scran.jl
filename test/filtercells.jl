using SparseArrays
using scran

@testset "Cell filtering tests" begin
    x = round.(abs.(sprand(Int32, 1000, 20, 0.2)) / 1e6)
    colsums = sum(x, dims=1)[1,:]
    mat = scran.initializesparsematrix(x)

    NC = size(x, 2)
    discard = Vector{Bool}(undef, NC)
    survivors = Vector{Int}()
    for i = 1:NC
        lose = i % 3 == 0
        discard[i] = lose
        if !lose
            push!(survivors, i)
        end
    end

    filtered = scran.filtercells(mat, discard)
    @test size(filtered, 2) == length(survivors)
    @test extractcolumn(mat, survivors[1]) == extractcolumn(filtered, 1)
    @test extractcolumn(mat, survivors[length(survivors)]) == extractcolumn(filtered, length(survivors))
end
