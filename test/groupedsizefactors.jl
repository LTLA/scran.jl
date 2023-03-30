using SparseArrays
using scran

@testset "Grouped size factor tests" begin
    x = round.(abs.(sprand(Int32, 1000, 20, 0.2)) / 1e6)
    mat = scran.initializesparsematrix(x)

    group = Vector{String}(undef, size(x, 2))
    levels = [ "A", "B", "C" ]
    for i in 1:size(x, 2)
        group[i] = levels[i % length(levels) + 1]
    end

    factors = scran.groupedsizefactors(mat, group)
    @test isapprox(sum(factors) / length(factors), 1)

    pfactors = scran.groupedsizefactors(mat, group; numthreads=2)
    @test factors == pfactors

    factors2 = scran.groupedsizefactors(mat, group; center = false)
    @test factors != factors2

    factors2 = scran.groupedsizefactors(mat, group; priorcount = 5)
    @test factors != factors2

    factors2 = scran.groupedsizefactors(mat, group; reference = "B")
    @test isapprox(sum(factors2) / length(factors2), 1)
end
