cd src/cpp
if [ ! -e build ]
then
    PREFIX_PATH=$(julia -e "using CxxWrap; println(CxxWrap.prefix_path())")
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$PREFIX_PATH
fi
cmake --build build --config Release
