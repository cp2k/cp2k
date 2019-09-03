# -*- Autoconf -*-

#####
#
# SYNOPSIS
#
# AX_CHECK_CUDA
#
# DESCRIPTION
#
# This macro search for the presence of the CUDA toolkit and set up several
# environement variables accordingly. The macro also search for clang and allows
# the user to use it instead of nvcc. In that case, it sets up additional
# variables that can be used in makefiles to compile cuda kernels.

# Note that cu files are not natively handled by autotools and libtool so a
# compilation rule should be defined accordingly. For instance, the following
# rule should do the trick in most cases

# if USE_CLANG_FOR_CUDA
# .cu.o :
# $(LIBTOOL) --tag=CC --mode=compile $(CLANG)  $(CUDA_FLAGS) \
# -I$(top_srcdir)/src/lib -I$(top_srcdir)/src/lib/Model -c $< -o $@
# .cu.lo :
# $(LIBTOOL) --tag=CC --mode=compile $(CLANG)  $(CUDA_FLAGS) -fPIC -I$(top_srcdir)/src/lib -c $< -o $@
# else
# .cu.o :
# $(LIBTOOL) --tag=CC --mode=compile $(NVCC)  $(NVCC_FLAGS) -I$(top_srcdir) $(NVCC_ADDITIONAL_FLAGS) -c $< -o $@
# .cu.lo :
# $(LIBTOOL) --tag=CC --mode=compile $(NVCC)  $(NVCC_FLAGS) -I$(top_srcdir) $(NVCC_ADDITIONAL_FLAGS) -c $< -o $@
# endif

# the macro defines the following variables
#
# general variables
#
# CUDA_CFLAGS : compilation parameters for the cuda libraries, include files, etc...
# CUDA_LDFLAGS : linker flags for the libraries
# CUDA_LIBS : flags for the libraries
#
# cuda kernels
#
# NVCC : location of nvcc
# CLANG : location of clang if found
# CUDA_FLAGS : flags for compiling cuda kernels CUDA_ADDIOTNAL_FLAGS : extra
# flags for indicating for instance the compiler to use to compile the C/C++
# parts of cu files.

# USE_CLANG_FOR_CUDA : indicate which compiler to choose for compiling cuda kernels
#
# config.h also contains a #define HAVE_CUDA if cuda is found
#
#
# the macro also defines an additional variable ac_have_cuda_ that can be used in
# configure.ac or others macros.


# LICENCE GPL
#
# AUTHOR
# Mathieu Taillefumier
#
#####

AU_ALIAS([ACX_CUDA], [AX_CUDA])
AC_DEFUN([AX_CUDA],
         [
         ac_have_cuda_=no

         NVCC=""
         NVCC_FLAGS=""
         CUDA_LDFLAGS=""
         # compiler used in combination with nvcc
         NVCC_ADDITIONAL_FLAGS=""
         CUDA_CFLAGS=""
         CUDA_LIBS=""

         # Provide your CUDA path with this
         AC_ARG_WITH(cuda, [AS_HELP_STRING([--with-cuda=yes/no/PATH],
                                           [compile with cuda support [default yes]])],
                     [cuda_prefix=$withval],
                     [cuda_prefix="/usr/local/cuda"])


         AS_IF([test "$with_cuda" != "no"], [
         # Setting the prefix to the default if only --with-cuda was given

         AS_IF([test "$cuda_prefix" == "yes"], [cuda_prefix="/usr/local/cuda"])

         ac_save_path_=$PATH
         PATH="$PATH:$cuda_prefix/bin"

         #search for nvcc
         AC_PATH_PROG(nvcc_compiler, nvcc, [no])
         AS_IF([test "$nvcc_compiler" != "no"], [cuda_prefix="${ac_cv_path_nvcc_compiler%%"/bin/nvcc"}"])
         CUDA_CFLAGS="-I$cuda_prefix/include"
         AC_ARG_WITH(cuda-compiler,
                     [AS_HELP_STRING([--with-cuda-compiler=nvcc/clang],
                                     [set which compiler should be used to compile cuda code [nvcc]])],
                     [cuda_compiler=$withval],
                     [cuda_compiler=nvcc])

        AC_ARG_WITH(cuda-compiler-options,
                    [AS_HELP_STRING([--with-cuda-compiler-options=FLAGS],
                                    [compilation parameters for cuda code [default -O2 -g]])],
                    [cuda_compiler_options=$withval],
                    [cuda_compiler_options="-O2 -g"])

        AC_ARG_WITH(cuda-arch,
                    [AS_HELP_STRING([--with-cuda-arch=sm_36],
                                    [GPU architecture : sm_20, sm_30, ... or P100, V100, K80 [default sm_35 (NVidia K80)]])],
                    [cuda_arch=$withval],
                    [cuda_arch="sm_37"])

        AC_ARG_WITH(cuda-cc-compiler,
                    [AS_HELP_STRING([--with-cuda-cc-compiler=compiler],
                                    [compiler used for compiling c/c++ codes included in cu files [default $CC]])],
                    [cuda_cc_compiler=$withval],
                    [cuda_cc_compiler=$CC])

        #check that the compiler used in combination with nvcc is present
        AS_IF([test "$cuda_cc_compiler" != "$CC"],
              [AC_CHECK_PROG(NVCC_CC_COMPILER,
                             $cuda_cc_compiler,
                             [yes], [no])])
       AS_IF([test "$NVCC_CC_COMPILER" == "yes"],
             [NVCC_ADDITIONAL_FLAGS="-ccbin $with_cuda_cc_compiler"])

       # check for cuda environment

       # first include files and libraries

       #check for 64 bits versions of the library
       AC_CHECK_FILE([$cuda_prefix/lib64],[cuda_lib64_found=yes],[cuda_lib64_found=no cuda_found=no])
       #check for 32 bits versions of the library
       AC_CHECK_FILE([$cuda_prefix/lib],[cuda_lib32_found=yes],[cuda_lib32_found=no])

       AC_CHECK_FILE($cuda_prefix/include/cuda.h,
                      [ac_have_cuda_=yes],
                      [ac_have_cuda_=no])

       AC_CHECK_LIB([cuda],
                    [cuInit],
                    [ac_have_cuda_=yes],
                    [ac_have_cuda_=no])

       # Saving the current flags
       ax_save_CFLAGS="${CFLAGS}"
       ax_save_LDFLAGS="${LDFLAGS}"

       AS_IF([test "$cuda_lib64_found" == "yes"],
             [CUDA_CFLAGS+=" -m64"
             CUDA_LDFLAGS="-L$cuda_prefix/lib64 -Wl,-rpath,$cuda_prefix/lib64"
             CUDA_LIBS="-lcublas -lcufft -lcurand -lcusparse -lcudart -L/usr/lib64/nvidia -lnvrtc -lcuda"],
             [test "$$cuda_lib32_found" == "yes"],
             [CUDA_LDFLAGS="-L$cuda_prefix/lib -Wl,-rpath,$cuda_prefix/lib"
             CUDA_LIBS="-lcublas -lcufft -lcurand -lcusparse -lcudart -L/usr/lib/nvidia -lnvrtc -lcuda"],
             [ac_have_cuda_=no])

             AS_IF([test "$cuda_compiler" == "nvcc" -o "$cuda_compiler" == "clang"],
             [ac_have_cuda_=yes],
             [ac_have_cuda_=no
            AC_MSG_ERROR([invalid parameter for the option --with-cuda-compiler. Please use --help for valid values])])

        # setup the default options for compiling gpu code
        NVCC_FLAGS="$cuda_compiler_options"

        # search for clang
        AC_PATH_PROG(clang_compiler, clang, [no])

        # check if nvcc or clang are found. one of them should be present to compile cuda code
        AS_IF([test "$nvcc_compiler" != "no" -o "$clang_compiler" = "yes"],
              [ac_have_cuda_=yes],
              [AC_MSG_ERROR("could not find clang or nvcc"); ac_have_cuda_=no])

        AS_CASE($cuda_arch,
                ["K80"], [GPUVER="K80"; cuda_arch="sm_37"],
                ["P100"], [GPUVER="P100"; cuda_arch="sm_60"],
                ["P40"], [GPUVER="P40"; cuda_arch="sm_61"],
                ["P6"], [GPUVER="P6"; cuda_arch="sm_61"],
                ["P4"], [GPUVER="P4"; cuda_arch="sm_61"],
                ["V100"], [GPUVER="V100"; cuda_arch="sm_70"],
                ["T4"], [GPUVER="T4"; cuda_arch="sm_72"],
                ["compute_30"], [],
                ["compute_32"], [],
                ["compute_35"], [],
                ["compute_37"], [GPUVER="K80"],
                ["compute_50"], [],
                ["compute_52"], [],
                ["compute_53"], [],
                ["compute_60"], [GPUVER="P100"],
                ["compute_61"], [GPUVER="P6"],
                ["compute_62"], [],
                ["compute_70"], [GPUVER="V100"],
                ["compute_72"], [],
                ["compute_75"], [],
                ["compute_80"], [],
                ["sm_30"], [],
                ["sm_32"], [],
                ["sm_35"], [],
                ["sm_37"], [GPUVER="K80"],
                ["sm_50"], [],
                ["sm_52"], [],
                ["sm_53"], [],
                ["sm_60"], [GPUVER="P100"],
                ["sm_61"], [GPUVER="P6"],
                ["sm_62"], [],
                ["sm_70"], [GPUVER="V100"],
                ["sm_72"], [],
                ["sm_75"], [],
                ["sm_80"], [],
                [AC_MSG_WARN([Unknown GPU architecture. will default to sm_36])
                GPUVER="K80"])


        AS_CASE($cuda_compiler,
                [nvcc], [NVCC=$cuda_prefix/bin/nvcc
                         NVCC_FLAGS+=" -arch $cuda_arch"],
                [clang], [ NVCC=clang dnl clang can compile cuda code but it needs a few extra options by default
                           CLANG=clang
                           NVCC_FLAGS+=" --cuda-gpu-arch=$cuda_arch --cuda-path=$cuda_prefix"],
                [AC_MSG_ERROR("Invalid value for the cuda compiler. Should be nvcc or clang")
                ac_have_cuda_=no])

        AS_IF([test "$ac_have_cuda_" == "yes"],
              [AC_DEFINE(HAVE_CUDA, 1, [cuda support enabled])])

    ])

    dnl everything for cuda libraries
    AC_SUBST([CUDA_LIBS])
    AC_SUBST([CUDA_LDFLAGS])
    AC_SUBST([CUDA_CFLAGS])

    dnl everything for compiling kernels
    AC_SUBST([NVCC])
    AC_SUBST([CLANG])
    AC_SUBST([NVCC_FLAGS])
    AC_SUBST([NVCC_ADDITIONAL_FLAGS])
    AM_CONDITIONAL([USE_CUDA], [test "$ac_have_cuda_" == "yes"])
    AM_CONDITIONAL([USE_CLANG_FOR_CUDA], [test "$cuda_compiler" == "clang"])

    AS_IF([test "$ac_have_cuda_" == "yes"], [$1], [$2])
])
