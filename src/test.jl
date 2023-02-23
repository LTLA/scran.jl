# Load the module and generate the functions
module scran
  using CxxWrap
  @wrapmodule(joinpath("cpp/build/lib", "libscrancxx"), :define_module_initialize_from_memory)

  function __init__()
    @initcxx
  end
end

# Call greet and show the result
using SparseArrays

i = Vector{Int32}([0, 3, 4, 5, 7, 7, 2, 3, 4, 5, 2, 5, 0, 2, 4, 6, 1, 4, 6, 1, 5, 7, 0, 3, 7, 7, 1, 6, 0, 5]);
x = Vector{Int32}([7, 1, 8, 9, 10, 2, 12, 4, 2, 3, 8, 36, 5, 2, 24, 12, 1, 6, 4, 7, 1, 14, 6, 18, 3, 19, 10, 3, 12, 17])
p = Vector{Int64}([0, 5, 6, 10, 12, 16, 19, 22, 25, 26, 28, 29, 30])

mat = scran.initialize_from_memory_integer(x, i, p, 8, 12, false, false, true)
@show mat

row = Vector{Float64}(undef, 12)
scran.row(mat, 2, row)
@show row
