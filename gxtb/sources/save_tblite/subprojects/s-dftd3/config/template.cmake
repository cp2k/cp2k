@PACKAGE_INIT@

set("@PROJECT_NAME@_WITH_API" @WITH_API@)
set("@PROJECT_NAME@_WITH_OpenMP" @WITH_OpenMP@)
set("@PROJECT_NAME@_WITH_BLAS" @WITH_BLAS@)
set(
  "@PROJECT_NAME@_INCLUDE_DIRS"
  "@CMAKE_INSTALL_PREFIX@/@CMAKE_INSTALL_INCLUDEDIR@"
  "@CMAKE_INSTALL_PREFIX@/@CMAKE_INSTALL_INCLUDEDIR@/@module-dir@"
)

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")

  include(CMakeFindDependencyMacro)

  if(NOT TARGET "OpenMP::OpenMP_Fortran" AND "@PROJECT_NAME@_WITH_OpenMP")
    find_dependency("OpenMP")
  endif()

  if(NOT TARGET "BLAS::BLAS" AND "@PROJECT_NAME@_WITH_BLAS")
    find_dependency("BLAS")
  endif()

  if(NOT TARGET "mctc-lib::mctc-lib")
    find_dependency("mctc-lib")
  endif()
endif()
