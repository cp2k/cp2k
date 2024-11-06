add_library(CP2K_OPENPMD_BINDINGS STATIC
    src/binding/openPMD/openPMD.cpp
)

target_link_libraries(CP2K_OPENPMD_BINDINGS PRIVATE openPMD::openPMD)
target_include_directories(
    CP2K_OPENPMD_BINDINGS
    PRIVATE src/binding/openPMD
)

add_executable(cp2k_test_openPMD
    src/binding/openPMD/testing/main.c)
target_link_libraries(cp2k_test_openPMD PRIVATE CP2K_OPENPMD_BINDINGS)
target_include_directories(
    cp2k_test_openPMD
    PRIVATE
    $<TARGET_PROPERTY:CP2K_OPENPMD_BINDINGS,INCLUDE_DIRECTORIES>
)
