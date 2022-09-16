# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindLAPACK
----------

Find Linear Algebra PACKage (LAPACK) library

This module finds an installed Fortran library that implements the
LAPACK linear-algebra interface (see http://www.netlib.org/lapack/).

Input Variables
^^^^^^^^^^^^^^^

The following variables may be set to influence this module's behavior:

``BLAS_STATIC``
  if ``ON`` use static linkage

  ``BLAS_THREADS``
  if set, check for threading support. List of possible options
  * ``sequential``
  * ``thread``
  * ``gnu-thread``
  * ``openmp``
  * ``tbb``

  tbb is specific to Intel MKL. thread and gnu-thread are synonimous in most
  systems, openmp is used mostly in openblas, sequential is defined for all
  libraries

``BLAS_VENDOR``
  If set, checks only the specified vendor, if not set checks all the
  possibilities.  List of vendors valid in this module:

  * ``OpenBLAS``
  * ``FLAME``
  * ``MKL``
  * ``ACML``
  * ``Generic``

``BLAS_F95``
  if ``ON`` tries to find the BLAS95/LAPACK95 interfaces

Result Variables
^^^^^^^^^^^^^^^^

This module defines the following variables:

``LAPACK_FOUND``
  library implementing the LAPACK interface is found
``LAPACK_LINKER_FLAGS``
  uncached list of required linker flags (excluding ``-l`` and ``-L``).
``LAPACK_LIBRARIES``
  uncached list of libraries (using full path name) to link against
  to use LAPACK
``LAPACK95_LIBRARIES``
  uncached list of libraries (using full path name) to link against
  to use LAPACK95
``LAPACK95_FOUND``
  library implementing the LAPACK95 interface is found

.. note::

  C, CXX or Fortran must be enabled to detect a BLAS/LAPACK library.
  C or CXX must be enabled to use Intel Math Kernel Library (MKL).

  For example, to use Intel MKL libraries and/or Intel compiler:

  .. code-block:: cmake

  set(BLAS_VENDOR "Generic")
  find_package(LAPACK)
#]=======================================================================]

find_package(PkgConfig)

# check for blas first. Most of the vendor libraries bundle lapack and blas in
# the same library. (MKL, and OPENBLAS do)

find_package(PkgConfig)

find_package(Blas REQUIRED)

if(CP2K_BLAS_FOUND)
  # LAPACK in the Intel MKL 10+ library?
  if(CP2K_BLAS_VENDOR MATCHES "MKL"
      OR CP2K_BLAS_VENDOR MATCHES "OpenBLAS"
      OR CP2K_BLAS_VENDOR MATCHES "Arm"
      OR CP2K_BLAS_VENDOR MATCHES "SCI")
    # we just need to create the interface that's all
    set(CP2K_LAPACK_FOUND TRUE)
    get_target_property(CP2K_LAPACK_INCLUDE_DIRS CP2K_BLAS::blas
      INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(CP2K_LAPACK_LIBRARIES CP2K_BLAS::blas INTERFACE_LINK_LIBRARIES)
  else()

    # we might get lucky to find a pkgconfig package for lapack (fedora provides
    # one for instance)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(CP2K_LAPACK lapack)
    endif()

    if(NOT CP2K_LAPACK_FOUND)
      find_library(
        CP2K_LAPACK_LIBRARIES
        NAMES "lapack" "lapack64"
        PATH_SUFFIXES "openblas" "openblas64" "openblas-pthread"
                      "openblas-openmp" "lib" "lib64"
        NO_DEFAULT_PATH)
    endif()
  endif()
endif()

# check if found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Lapack REQUIRED_VARS CP2K_LAPACK_LIBRARIES)

if(CP2K_LAPACK_FOUND)
  if(NOT TARGET CP2K_LAPACK::lapack)
    add_library(CP2K_LAPACK::lapack INTERFACE IMPORTED)
    set_property(TARGET CP2K_LAPACK::lapack PROPERTY INTERFACE_LINK_LIBRARIES
      ${CP2K_LAPACK_LIBRARIES})
    if (CP2K_LAPACK_INCLUDE_DIRS)
      set_property(TARGET CP2K_LAPACK::lapack PROPERTY INTERFACE_INCLUDE_DIRECTORIES
        ${CP2K_LAPACK_INCLUDE_DIRS})
    endif()
    add_library(CP2K_LAPACK::LAPACK ALIAS CP2K_LAPACK::lapack)
  endif()
endif()

# prevent clutter in cache
mark_as_advanced(CP2K_LAPACK_FOUND CP2K_LAPACK_LIBRARIES CP2K_LAPACK_INCLUDE_DIRS)
