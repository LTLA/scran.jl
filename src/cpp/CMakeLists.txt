project(scranlib)

cmake_minimum_required(VERSION 3.14)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

find_package(JlCxx)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")

set(CMAKE_CXX_STANDARD 17)

add_library(scrancxx SHARED 
    src/julia.cpp
    src/initialize_from_memory.cpp
    src/log_norm_counts.cpp
    src/per_cell_rna_qc_metrics.cpp
    src/suggest_rna_qc_filters.cpp
    src/per_cell_adt_qc_metrics.cpp
    src/suggest_adt_qc_filters.cpp
    src/per_cell_crispr_qc_metrics.cpp
    src/suggest_crispr_qc_filters.cpp
    src/filter_cells.cpp
    src/model_gene_var.cpp
)

add_subdirectory(extern)

target_link_libraries(scrancxx JlCxx::cxxwrap_julia)
target_link_libraries(scrancxx scran mnncorrect qdtsne umappp)

install(TARGETS scrancxx
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION lib)

