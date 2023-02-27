using Test
@testset "scran.jl" begin
using scran, Documenter
#doctest(scran)
include("initializesparsematrix.jl")
include("lognormcounts.jl")
end
