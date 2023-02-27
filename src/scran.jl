# Load the module and generate the functions
module scran

using CxxWrap
@wrapmodule(joinpath(@__DIR__, "cpp/build/lib", "libscrancxx"), :define_julia_module)

function __init__()
    @initcxx
end

include("initialize.jl")
export initializesparsematrix, size, extractrow, extractcolumn

include("lognormcounts.jl")
export lognormcounts

end