#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

set(CP2K_C_COMPILER_LIST "GNU;Intel;IntelLLVM;NAG;Cray;PGI;Clang;AppleClang")
set(CP2K_Fortran_COMPILER_LIST "GNU;Intel;IntelLLVM;NAG;Cray;PGI")

if(NOT CMAKE_C_COMPILER_ID IN_LIST CP2K_C_COMPILER_LIST)
  message(
    WARNING
      "Unknown C compiler, trying without any additional (optimization) flags.\n"
      "You may have to specify flags for C11 support manually.\n Please open an\n"
      "issue at https://github.com/cp2k/cp2k/issues with the reported compiler \n"
      "name and the required flags.")

  message("-- CMAKE_C_COMPILER_ID: ${CMAKE_C_COMPILER_ID}\n"
          "-- CMAKE_C_COMPILER full path: ${CMAKE_C_COMPILER}\n")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID IN_LIST CP2K_Fortran_COMPILER_LIST)
  message(
    WARNING
      "Unknown Fortran compiler, trying without any additional (optimization) flags.\n"
      "You will most likely have to specify extra options for free-form Fortran 2008\n"
      "for your compiler! Please open an issue at https://github.com/cp2k/cp2k/issues\n"
      "with the reported compiler name and the required flags.")
  message("-- CMAKE_Fortran_COMPILER_ID: ${CMAKE_Fortran_COMPILER_ID}\n"
          "-- CMAKE_Fortran_COMPILER full path: ${CMAKE_Fortran_COMPILER}\n")
endif()

# ================================ GNU Compilers ===============================

# Baseline
add_compile_options(
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2008;-ffree-form;-fimplicit-none>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-g;-fno-omit-frame-pointer;-fbacktrace>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:$<$<VERSION_GREATER:$<Fortran_COMPILER_VERSION>,10>:-fallow-argument-mismatch>>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-deprecated-declarations;-Wno-maybe-uninitialized;-Wuninitialized;-Wuse-without-only>"
)
add_compile_options(
  "$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-g;-fno-omit-frame-pointer>"
  "$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-deprecated-declarations;-Wno-vla-parameter>"
)
add_compile_options(
  "$<$<COMPILE_LANG_AND_ID:C,GNU>:-g;-fno-omit-frame-pointer>"
  "$<$<COMPILE_LANG_AND_ID:C,GNU>:-Wno-deprecated-declarations;-Wno-vla-parameter>"
)

# -- Apple Silicon + GCC: -march=native expands internally to -march=apple-m1
# (invalid). Use -mcpu=native instead.
set(_CP2K_GNU_NATIVE_TUNE "-march=native;-mtune=native")
if(APPLE
   AND CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64"
   AND (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
        OR CMAKE_C_COMPILER_ID STREQUAL "GNU"
        OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU"))
  set(_CP2K_GNU_NATIVE_TUNE "-mcpu=native")
endif()
if(APPLE)
  add_definitions(-D__MACOSX)
endif()
if(APPLE OR BSD)
  add_definitions(-D__NO_STATM_ACCESS)
endif()

# Release
add_compile_options(
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-O3;${_CP2K_GNU_NATIVE_TUNE};-funroll-loops>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:CXX,GNU>>:-O3;${_CP2K_GNU_NATIVE_TUNE};-funroll-loops>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,GNU>>:-O3;${_CP2K_GNU_NATIVE_TUNE};-funroll-loops>"
)

# Generic
add_compile_options(
  "$<$<AND:$<CONFIG:GENERIC>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-O3;-mtune=generic;-funroll-loops>"
  "$<$<AND:$<CONFIG:GENERIC>,$<COMPILE_LANG_AND_ID:CXX,GNU>>:-O3;-mtune=generic;-funroll-loops>"
  "$<$<AND:$<CONFIG:GENERIC>,$<COMPILE_LANG_AND_ID:C,GNU>>:-O3;-mtune=generic;-funroll-loops>"
)

# Debug
add_compile_options(
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-O1;${_CP2K_GNU_NATIVE_TUNE}>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:CXX,GNU>>:-O1;${_CP2K_GNU_NATIVE_TUNE}>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,GNU>>:-O1;${_CP2K_GNU_NATIVE_TUNE};-Wall;-Wextra;-Werror>"
)
add_compile_options(
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-fsanitize=leak;-Werror=realloc-lhs>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-fcheck=all,no-array-temps;-finline-matmul-limit=0>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-ffpe-trap=invalid,zero,overflow>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-finit-derived;-finit-real=snan;-finit-integer=-42>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-Werror=aliasing;-Werror=ampersand;-Werror=c-binding-type>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-Werror=character-truncation;-Werror=intrinsic-shadow>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-Werror=intrinsics-std;-Werror=line-truncation;-Werror=tabs;-Werror=target-lifetime>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-Werror=underflow;-Werror=unused-but-set-variable;-Werror=unused-variable>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-Werror=unused-dummy-argument;-Werror=unused-parameter>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-Werror=unused-label;-Werror=conversion;-Werror=zerotrip>"
)
add_compile_definitions("$<$<CONFIG:DEBUG>:__HAS_IEEE_EXCEPTIONS;__CHECK_DIAG>")
add_link_options(
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-fsanitize=leak>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:CXX,GNU>>:-fsanitize=leak>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,GNU>>:-fsanitize=leak>")

# Conventions
add_compile_options(
  "$<$<CONFIG:CONVENTIONS>:-O1;${_CP2K_GNU_NATIVE_TUNE}>"
  "$<$<AND:$<CONFIG:CONVENTIONS>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-Wno-pedantic;-Wall;-Wextra;-Wsurprising>"
  "$<$<AND:$<CONFIG:CONVENTIONS>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-Warray-temporaries;-Wconversion-extra;-Wimplicit-interface>"
  "$<$<AND:$<CONFIG:CONVENTIONS>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-Wimplicit-procedure;-Wreal-q-constant;-Walign-commons>"
  "$<$<AND:$<CONFIG:CONVENTIONS>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-Wfunction-elimination;-Wrealloc-lhs;-Wcompare-reals;-Wzerotrip>"
  # Writes a lot of output. Make sure to use:
  # -DCMAKE_Fortran_COMPILER_LAUNCHER="redirect_gfortran_output.py"
  "$<$<AND:$<CONFIG:CONVENTIONS>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-fdump-fortran-original>"
)

# Coverage
add_compile_options(
  "$<$<CONFIG:COVERAGE>:-coverage;-fkeep-static-functions;-O1;${_CP2K_GNU_NATIVE_TUNE}>"
)
add_compile_definitions("$<$<CONFIG:COVERAGE>:__NO_ABORT>")

# Address Sanitizer
add_compile_options(
  "$<$<CONFIG:ASAN>:-fsanitize=address;-no-pie;-O3;${_CP2K_GNU_NATIVE_TUNE};-funroll-loops>"
)
add_link_options("$<$<CONFIG:ASAN>:-fsanitize=address>")

# =============================== Other Compilers ==============================

# Baseline
add_compile_options(
  "$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-free;-stand f18;-fpp;-qopenmp;-heap-arrays;-D__MAX_CONTR=4>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-free;-fpp;-qopenmp;-D__MAX_CONTR=4>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,PGI>:-Mfreeform;-Mextend;-Mallocatable=03>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-f2008;-free;-Warn=reallocation;-Warn=subnormal>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,Cray>:-f;free;-M3105;-ME7212;-hnoacc;-M1234>"
)
add_compile_options(
  "$<$<COMPILE_LANG_AND_ID:C,Cray>:-hnoacc;-h;nomessage=1234>")

# Release
add_compile_options(
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,PGI>>:-fast>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:-O3;-g;-D__HAS_IEEE_EXCEPTIONS>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>>:-O3;-g;-D__HAS_IEEE_EXCEPTIONS>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,Cray>>:-O2;-G2>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,NAG>>:-gline>")
add_compile_options(
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,PGI>>:-fast>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,Intel>>:-O3;-g>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,IntelLLVM>>:-O3;-g>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,Cray>>:-O3>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,NAG>>:-gline>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,Clang>>:-O3;-funroll-loops>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,AppleClang>>:-O3;-funroll-loops>"
  "$<$<NOT:$<BOOL:OpenMP_FOUND>>:$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,NAG>>:-gline>>"
)

# Debug
add_compile_options(
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,PGI>>:-g>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:-O2;-debug>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>>:-O2;-debug>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,Cray>>:-G2>")
add_compile_options(
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,PGI>>:-fast>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,Intel>>:-O2;-g>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,IntelLLVM>>:-O2;-g>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,Cray>>:-G2>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,NAG>>:-g;-C>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,Clang>>:-O1;-g;-fno-omit-frame-pointer>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,AppleClang>>:-O0;-g>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,NAG>>:-C=all>")

# =================================== Tweaks ===================================
# Workaround https://gitlab.kitware.com/cmake/cmake/-/issues/27231
if(TARGET MPI::MPI_Fortran)
  get_target_property(opts MPI::MPI_Fortran INTERFACE_COMPILE_OPTIONS)
  set_target_properties(
    MPI::MPI_Fortran PROPERTIES INTERFACE_COMPILE_OPTIONS
                                "$<$<COMPILE_LANGUAGE:Fortran>:${opts}>")
  unset(opts)
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
  set(CMAKE_Fortran_MODOUT_FLAG "-ef") # override to get lower-case module file
                                       # names
  set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "-lgcov") # Apple's Clang needs an extra
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  set(CMAKE_Fortran_MODOUT_FLAG "-ef") # override to get lower-case module file
                                       # names
endif()
if(CMAKE_C_COMPILER_ID STREQUAL "Cray" AND CMAKE_C_COMPILER_VERSION
                                           VERSION_LESS 9)
  # prevent deallocation failures due to tcmalloc's free with glibc's
  # aligned_alloc, see https://bugzilla.redhat.com/show_bug.cgi?id=1569391
  add_compile_options("$<$<LANGUAGE:C>:-h;system_alloc>"
                      "$<$<LANGUAGE:Fortran>:-h;system_alloc>")
  # since the detection of the implicitly linked libraries occurs before we can
  # intervene, filter them out again
  list(FILTER CMAKE_C_IMPLICIT_LINK_LIBRARIES EXCLUDE REGEX "tcmalloc")
  list(FILTER CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES EXCLUDE REGEX "tcmalloc")
endif()
