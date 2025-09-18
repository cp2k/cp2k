add_library(CP2K_OPENPMD_BINDINGS STATIC openPMD.cpp)

find_package(openPMD 0.16.1 REQUIRED)
target_link_libraries(CP2K_OPENPMD_BINDINGS PRIVATE openPMD::openPMD)

find_package(MPI REQUIRED)
target_link_libraries(CP2K_OPENPMD_BINDINGS PRIVATE MPI::MPI_C)

target_include_directories(CP2K_OPENPMD_BINDINGS PRIVATE .)
