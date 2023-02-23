# Load the module and generate the functions
module CppHello
  using CxxWrap
  @wrapmodule(joinpath("cpp/build/lib", "libscrancxx"))

  function __init__()
    @initcxx
  end
end

# Call greet and show the result
@show CppHello.greet()
