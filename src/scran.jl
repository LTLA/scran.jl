# Load the module and generate the functions
module scran

using CxxWrap
@wrapmodule(joinpath(@__DIR__, "cpp/build/lib", "libscrancxx"), :define_julia_module)

function __init__()
    @initcxx
end

include("utils.jl")

include("initialize.jl")
export initializesparsematrix, size, extractrow, extractcolumn

include("lognormcounts.jl")
export lognormcounts

include("groupedsizefactors.jl")
export groupedsizefactors

include("percellrnaqcmetrics.jl")
export percellrnaqcmetrics

include("suggestrnaqcfilters.jl")
export suggestrnaqcfilters

include("percelladtqcmetrics.jl")
export percelladtqcmetrics

include("suggestadtqcfilters.jl")
export suggestadtqcfilters

include("percellcrisprqcmetrics.jl")
export percellcrisprqcmetrics

include("suggestcrisprqcfilters.jl")
export suggestcrisprqcfilters

include("filtercells.jl")
export filtercells

include("modelgenevar.jl")
export modelgenevar

include("runpca.jl")
export runpca

end
