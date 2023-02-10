#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS
      "-mtune=native -ffree-line-length-512 -ffree-form -std=f2008 -fimplicit-none -Werror=aliasing -Werror=ampersand -Werror=c-binding-type -Werror=intrinsic-shadow -Werror=intrinsics-std -Werror=tabs -Werror=target-lifetime -Werror=underflow -Werror=unused-but-set-variable -Werror=unused-variable  -Werror=conversion -Werror=zerotrip -Wuninitialized -Wno-maybe-uninitialized -Wunused-parameter"
  )

  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch"
    )# gcc 10+ has this automatically
  endif()

  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -g -funroll-loops")
  set(CMAKE_Fortran_FLAGS_COVERAGE
      "-O0 -g --coverage -fno-omit-frame-pointer -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-real=snan -finit-integer=-42 -finit-derived -Werror=realloc-lhs -finline-matmul-limit=0"
  )
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-fbacktrace -ffree-form -fimplicit-none -std=f2008 -fsanitize=leak -fcheck=all,no-array-temps -ffpe-trap=invalid,zero,overflow -finit-derived -finit-real=snan -finit-integer=-42 -Werror=realloc-lhs -finline-matmul-limit=0"
  )
  if((NOT (USE_MPI)) OR (NOT ("${MPI_Fortran_LIBRARY_VERSION_STRING}" MATCHES
                              "Open MPI")))
    set(CMAKE_Fortran_FLAGS_COVERAGE
        "${CMAKE_Fortran_FLAGS_COVERAGE} -fsanitize=leak")
    set(CMAKE_Fortran_FLAGS_DEBUG
        "${CMAKE_Fortran_FLAGS_DEBUG} -fsanitize=leak")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "-free -stand=f18 -fpp -heap-arrays")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -g")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O2 -debug")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  set(CMAKE_Fortran_FLAGS "-Mfreeform -Mextend -Mallocatable=03"
  )# -Mallocatable=03: enable F2003+ assignment semantics
  set(CMAKE_Fortran_FLAGS_RELEASE "-fast")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
  set(CMAKE_Fortran_FLAGS "-f2008 -free -Warn=reallocation -Warn=subnormal")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -C")
  if(NOT OpenMP_FOUND)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -gline"
    )# -gline is only supported without OpenMP
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -C=all"
    )# some checks are not available with OpenMP
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  set(CMAKE_Fortran_FLAGS "-f free -M3105 -ME7212") # -M3105: hide a
                                                    # false-positive warning
                                                    # about modified loop
                                                    # variables due to loop
                                                    # fusing, promote warning
                                                    # 7212 to an error
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -G2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-G2")
  set(CMAKE_Fortran_MODOUT_FLAG "-ef") # override to get lower-case module file
                                       # names
else()
  message(
    WARNING
      "\
Unknown Fortran compiler, trying without any additional (optimization) flags.\n\
You will most likely have to specify extra options for free-form Fortran 2008 for your compiler!\n\
Please open an issue at https://github.com/cp2k/dbcsr/issues with the reported compiler name and the required flags."
  )
  message("-- CMAKE_Fortran_COMPILER_ID: " ${CMAKE_Fortran_COMPILER_ID})
  message("-- CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -g -funroll-loops -Wall")
  set(CMAKE_C_FLAGS_COVERAGE "-O0 -g --coverage -Wall")
  set(CMAKE_C_FLAGS_DEBUG
      "-O2 -ggdb -Wall -fsanitize=undefined -fsanitize=address -fsanitize-recover=all -Wall -Wextra -Werror -Wno-vla-parameter -Wno-deprecated-declarations"
  )
  if((NOT (USE_MPI)) OR (NOT ("${MPI_Fortran_LIBRARY_VERSION_STRING}" MATCHES
                              "Open MPI")))
    set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_C_FLAGS_COVERAGE} -fsanitize=leak")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=leak")
  endif()
elseif(CMAKE_C_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -funroll-loops")
  set(CMAKE_C_FLAGS_COVERAGE "-O0 -g --coverage")
  set(CMAKE_C_FLAGS_DEBUG "-O1 -g -fno-omit-frame-pointer")
elseif(CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -funroll-loops")
  set(CMAKE_C_FLAGS_COVERAGE "-O0 -g --coverage")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
  set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "-lgcov") # Apple's Clang needs an extra
elseif(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -g")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -debug")
elseif(CMAKE_C_COMPILER_ID STREQUAL "PGI")
  set(CMAKE_C_FLAGS_RELEASE "-fast")
  set(CMAKE_C_FLAGS_DEBUG "-g")
elseif(CMAKE_C_COMPILER_ID STREQUAL "Cray")
  set(CMAKE_C_FLAGS_RELEASE "-O3")
  set(CMAKE_C_FLAGS_DEBUG "-G2")
  if(CMAKE_C_COMPILER_VERSION VERSION_LESS 9)
    # prevent deallocation failures due to tcmalloc's free with glibc's
    # aligned_alloc, see https://bugzilla.redhat.com/show_bug.cgi?id=1569391
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -h system_alloc")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -h system_alloc")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -h system_alloc")
    # since the detection of the implicitly linked libraries occurs before we
    # can intervene, filter them out again
    list(FILTER CMAKE_C_IMPLICIT_LINK_LIBRARIES EXCLUDE REGEX "tcmalloc")
    list(FILTER CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES EXCLUDE REGEX "tcmalloc")
  endif()
  # OpenACC support with CCE is EOL: causes
  # https://github.com/cp2k/dbcsr/issues/261 eventually check compiler version
  # (similar to -h system_alloc)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -hnoacc -h nomessage=1234")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -hnoacc -h nomessage=1234")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -hnoacc -M1234")
else()
  message(
    WARNING
      "\
Unknown C++ compiler, trying without any additional (optimization) flags.\n\
You may have to specify flags for C++11/14 support manually.\n\
Please open an issue at https://github.com/cp2k/dbcsr/issues with the reported compiler name and the required flags."
  )
  message("-- CMAKE_C_COMPILER_ID: " ${CMAKE_C_COMPILER_ID})
  message("-- CMAKE_C_COMPILER full path: " ${CMAKE_C_COMPILER})
endif()
