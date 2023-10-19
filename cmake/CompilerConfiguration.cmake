#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

set(CP2K_C_COMPILER_LIST "GNU;Intel;NAG;Cray;PGI;Clang;AppleClang")
set(CP2K_Fortran_COMPILER_LIST "GNU;Intel;NAG;Cray;PGI")

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

# OpenACC support with CCE is EOL: causes
# https://github.com/cp2k/dbcsr/issues/261 eventually check compiler version
# (similar to -h system_alloc)
add_compile_options(
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-mtune=native;-ffree-line-length-512;-ffree-form;-std=f2008;-fimplicit-none;-Werror=aliasing;-Werror=ampersand;-Werror=c-binding-type;-Werror=conversion;-Werror=intrinsic-shadow;-Werror=intrinsics-std;-Werror=line-truncation;-Werror=tabs;-Werror=target-lifetime;-Werror=underflow;-Werror=unused-but-set-variable;-Werror=unused-variable>"
  "$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<VERSION_GREATER_EQUAL:${CMAKE_Fortran_COMPILER_VERSION},11>>:-fallow-argument-mismatch>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-free -stand=f18 -fpp -heap-arrays>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,PGI>:-Mfreeform -Mextend -Mallocatable=03>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-f2008 -free -Warn=reallocation -Warn=subnormal>"
  "$<$<COMPILE_LANG_AND_ID:C,Cray>:-hnoacc -h nomessage=1234>"
  "$<$<COMPILE_LANG_AND_ID:Fortran,Cray>:-f free -M3105 -ME7212  -hnoacc -M1234>"
)

add_compile_options(
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-O3;-g;-funroll-loops>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,PGI>>:-fast>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:-O3;-g>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,Cray>>:-O2;-G2>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,NAG>>:-gline>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,GNU>>:-O3;-g;-funroll-loops;-Wall>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,PGI>>:-fast>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,Intel>>:-O3;-g>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,Cray>>:-O3>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,NAG>>:-gline>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,Clang>>:-O3;-funroll-loops>"
  "$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:C,AppleClang>>:-O3;-funroll-loops>"
  "$<$<NOT:$<BOOL:OpenMP_FOUND>>:$<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANG_AND_ID:Fortran,NAG>>:-gline>>"
)

add_compile_options(
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-fbacktrace;-ffree-form;-fimplicit-none;-std=f2008;-fsanitize=leak;-fcheck=all,no-array-temps;-ffpe-trap=invalid,zero,overflow;-finit-derived;-finit-real=snan;-finit-integer=-42;-Werror=realloc-lhs;-finline-matmul-limit=0>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,PGI>>:-g>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,Intel>>:-O2;-debug>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,Cray>>:-G2>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,GNU>>:-O2;-ggdb;-Wall;-fsanitize=undefined;-fsanitize=address;-fsanitize-recover=all;-Wall;-Wextra;-Werror;-Wno-vla-parameter;-Wno-deprecated-declarations>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,PGI>>:-fast>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,Intel>>:-O3;-g>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,Cray>>:-G2>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,NAG>>:-g;-C>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,Clang>>:-O1;-g;-fno-omit-frame-pointer>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:C,AppleClang>>:-O0;-g>"
  "$<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANG_AND_ID:Fortran,NAG>>:-C=all>")

add_compile_options(
  "$<$<AND:$<CONFIG:COVERAGE>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-O0;-g;--coverage;-fno-omit-frame-pointer;-fcheck=all;-ffpe-trap=invalid,zero,overflow;-fbacktrace;-finit-real=snan;-finit-integer=-42;-finit-derived;-Werror=realloc-lhs;-finline-matmul-limit=0>"
  "$<$<AND:$<CONFIG:COVERAGE>,$<COMPILE_LANG_AND_ID:C,GNU>>:-O0;-g;--coverage;-Wall>"
)

if(NOT CP2K_USE_MPI OR NOT "${MPI_Fortran_LIBRARY_VERSION_STRING}" MATCHES "Open
  MPI")
  add_compile_options(
    "$<$<AND:$<CONFIG:COVERAGE>,$<COMPILE_LANG_AND_ID:Fortran,GNU>>:-fsanitize=leak>"
    "$<$<AND:$<CONFIG:COVERAGE>,$<COMPILE_LANG_AND_ID:C,GNU>>:-O0;-g;--coverage;-Wall>"
  )
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
