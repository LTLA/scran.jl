using Test
@testset "scran.jl" begin
using scran, Documenter
#doctest(scran)
include("initializesparsematrix.jl")
include("lognormcounts.jl")
include("percellrnaqcmetrics.jl")
include("suggestrnaqcfilters.jl")
include("percelladtqcmetrics.jl")
include("suggestadtqcfilters.jl")
include("percellcrisprqcmetrics.jl")
include("suggestcrisprqcfilters.jl")
include("filtercells.jl")
include("modelgenevar.jl")
include("runpca.jl")
end
