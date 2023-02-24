# Load the module and generate the functions
module scran

using CxxWrap
@wrapmodule(joinpath(@__DIR__, "cpp/build/lib", "libscrancxx"), :define_module_initialize_from_memory)

function __init__()
    @initcxx
end

include("initialize.jl")
export initializesparsematrix, size, extractrow, extractcolumn

end
