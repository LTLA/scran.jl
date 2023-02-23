#include <string>
#include "jlcxx/jlcxx.hpp"

struct MyStruct { MyStruct() {} };

std::string greet() {
    return "Hello";
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.add_type<MyStruct>("MyStruct");
  mod.method("greet", &greet);
}
