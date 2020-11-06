#!/bin/bash -e

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

# ------------------------------------------------------------------------
# generate arch file for compiling cp2k
# ------------------------------------------------------------------------

echo "==================== generating arch files ===================="
echo "arch files can be found in the ${INSTALLDIR}/arch subdirectory"
! [ -f "${INSTALLDIR}/arch" ] && mkdir -p ${INSTALLDIR}/arch
cd ${INSTALLDIR}/arch

# -------------------------
# set compiler flags
# -------------------------

# need to switch between FC and MPICC etc in arch file, but cannot use
# same variable names, so use _arch suffix
CC_arch="$CC"
CXX_arch="$CXX"
FC_arch="IF_MPI(${MPIFC}|${FC})"
LD_arch="IF_MPI(${MPIFC}|${FC})"

# we always want good line information and backtraces
BASEFLAGS="-march=native -mtune=native -fno-omit-frame-pointer -g ${TSANFLAGS}"
OPT_FLAGS="-O3 -funroll-loops"
NOOPT_FLAGS="-O1"

# those flags that do not influence code generation are used always, the others if debug
FCDEB_FLAGS="-fbacktrace -ffree-form -fimplicit-none -std=f2008"
FCDEB_FLAGS_DEBUG="-fsanitize=leak -fcheck=all,no-array-temps -ffpe-trap=invalid,zero,overflow -finit-derived -finit-real=snan -finit-integer=-42 -Werror=realloc-lhs -finline-matmul-limit=0"

# code coverage generation flags
COVERAGE_FLAGS="-O1 -coverage -fkeep-static-functions"
COVERAGE_DFLAGS="-D__NO_ABORT"

# profile based optimization, see https://www.cp2k.org/howto:pgo
PROFOPT_FLAGS="\$(PROFOPT)"

# special flags for gfortran
# https://gcc.gnu.org/onlinedocs/gfortran/Error-and-Warning-Options.html
# we error out for these warnings (-Werror=uninitialized -Wno-maybe-uninitialized -> error on variables that must be used uninitialized)
WFLAGS_ERROR="-Werror=aliasing -Werror=ampersand -Werror=c-binding-type -Werror=intrinsic-shadow -Werror=intrinsics-std -Werror=line-truncation -Werror=tabs -Werror=target-lifetime -Werror=underflow -Werror=unused-but-set-variable -Werror=unused-variable -Werror=unused-dummy-argument -Werror=conversion -Werror=zerotrip -Wno-maybe-uninitialized"
# we just warn for those (that eventually might be promoted to WFLAGSERROR). It is useless to put something here with 100s of warnings.
WFLAGS_WARN="-Wuninitialized -Wuse-without-only"
# while here we collect all other warnings, some we'll ignore
# TODO: -Wpedantic with -std2008 requires an upgrade of the MPI interfaces from mpi to mpi_f08
WFLAGS_WARNALL="-Wno-pedantic -Wall -Wextra -Wsurprising -Wunused-parameter -Warray-temporaries -Wcharacter-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wreal-q-constant -Walign-commons -Wfunction-elimination -Wrealloc-lhs -Wcompare-reals -Wzerotrip"

# IEEE_EXCEPTIONS dependency
IEEE_EXCEPTIONS_DFLAGS="-D__HAS_IEEE_EXCEPTIONS"

# check all of the above flags, filter out incompatible flags for the
# current version of gcc in use
BASEFLAGS=$(allowed_gfortran_flags         $BASEFLAGS)
OPT_FLAGS=$(allowed_gfortran_flags         $OPT_FLAGS)
NOOPT_FLAGS=$(allowed_gfortran_flags       $NOOPT_FLAGS)
FCDEB_FLAGS=$(allowed_gfortran_flags       $FCDEB_FLAGS)
FCDEB_FLAGS_DEBUG=$(allowed_gfortran_flags $FCDEB_FLAGS_DEBUG)
COVERAGE_FLAGS=$(allowed_gfortran_flags    $COVERAGE_FLAGS)
WFLAGS_ERROR=$(allowed_gfortran_flags      $WFLAGS_ERROR)
WFLAGS_WARN=$(allowed_gfortran_flags       $WFLAGS_WARN)
WFLAGS_WARNALL=$(allowed_gfortran_flags    $WFLAGS_WARNALL)

# check if ieee_exeptions module is available for the current version
# of gfortran being used
if ! (check_gfortran_module ieee_exceptions) ; then
    IEEE_EXCEPTIONS_DFLAGS=""
fi

# concatenate the above flags into WFLAGS, FCDEBFLAGS, DFLAGS and
# finally into FCFLAGS and CFLAGS
WFLAGS="$WFLAGS_ERROR $WFLAGS_WARN IF_WARNALL(${WFLAGS_WARNALL}|)"
FCDEBFLAGS="$FCDEB_FLAGS IF_DEBUG($FCDEB_FLAGS_DEBUG|)"
DFLAGS="${CP_DFLAGS} IF_DEBUG($IEEE_EXCEPTIONS_DFLAGS -D__CHECK_DIAG|) IF_COVERAGE($COVERAGE_DFLAGS|)"
# language independent flags
# valgrind with avx can lead to spurious out-of-bound results
G_CFLAGS="$BASEFLAGS IF_VALGRIND(-mno-avx -mno-avx2|)"
G_CFLAGS="$G_CFLAGS IF_COVERAGE($COVERAGE_FLAGS|IF_DEBUG($NOOPT_FLAGS|$OPT_FLAGS))"
G_CFLAGS="$G_CFLAGS IF_DEBUG(|$PROFOPT_FLAGS)"
G_CFLAGS="$G_CFLAGS $CP_CFLAGS"
# FCFLAGS, for gfortran
FCFLAGS="$G_CFLAGS \$(FCDEBFLAGS) \$(WFLAGS) \$(DFLAGS)"
# CFLAGS, special flags for gcc
CFLAGS="$G_CFLAGS -std=c99 -Wall -Wextra -Werror \$(DFLAGS)"

# Linker flags
LDFLAGS="\$(FCFLAGS) ${CP_LDFLAGS}"

# Library flags
# add standard libs
LIBS="${CP_LIBS} -lstdc++"

# CUDA handling
CUDA_LIBS="-lcudart -lnvrtc -lcuda -lcufft -lcublas -lrt IF_DEBUG(-lnvToolsExt|)"
CUDA_DFLAGS="-D__COLLOCATE_GPU -D__ACC -D__DBCSR_ACC -D__PW_CUDA IF_DEBUG(-D__CUDA_PROFILING|)"
if [ "${ENABLE_CUDA}" = __TRUE__ ] && [ "${GPUVER}" != no ] ; then
    LIBS="${LIBS} IF_CUDA(${CUDA_LIBS}|)"
    DFLAGS="IF_CUDA(${CUDA_DFLAGS}|) ${DFLAGS}"
    NVFLAGS="-arch sm_${ARCH_NUM} -O3 -Xcompiler='-fopenmp' --std=c++11 \$(DFLAGS)"
    check_command nvcc "cuda"
    check_lib -lcudart "cuda"
    check_lib -lnvrtc "cuda"
    check_lib -lcuda "cuda"
    check_lib -lcufft "cuda"
    check_lib -lcublas "cuda"

    # Set include flags
    CUDA_CFLAGS=''
    add_include_from_paths CUDA_CFLAGS "cuda.h" $INCLUDE_PATHS
    export CUDA_CFLAGS="${CUDA_CFLAGS}"
    CFLAGS+=" ${CUDA_CFLAGS}"

    # Set LD-flags
    CUDA_LDFLAGS=''
    add_lib_from_paths CUDA_LDFLAGS "libcudart.*" $LIB_PATHS
    add_lib_from_paths CUDA_LDFLAGS "libnvrtc.*" $LIB_PATHS
    add_lib_from_paths CUDA_LDFLAGS "libcuda.*" $LIB_PATHS
    add_lib_from_paths CUDA_LDFLAGS "libcufft.*" $LIB_PATHS
    add_lib_from_paths CUDA_LDFLAGS "libcublas.*" $LIB_PATHS
    export CUDA_LDFLAGS="${CUDA_LDFLAGS}"
    LDFLAGS+=" ${CUDA_LDFLAGS}"
fi

# -------------------------
# generate the arch files
# -------------------------

# generator for CP2K ARCH files
gen_arch_file() {
    # usage: gen_arch_file file_name flags
    #
    # If the flags are present they are assumed to be on, otherwise
    # they switched off
    require_env ARCH_FILE_TEMPLATE
    local __filename=$1
    shift
    local __flags=$@
    local __full_flag_list="MPI DEBUG CUDA WARNALL VALGRIND COVERAGE"
    local __flag=''
    for __flag in $__full_flag_list ; do
        eval "local __${__flag}=off"
    done
    for __flag in $__flags ; do
        eval "__${__flag}=on"
    done
    # generate initial arch file
    cat $ARCH_FILE_TEMPLATE > $__filename
    # add additional parts
    if [ "$__CUDA" = "on" ] ; then
      cat <<EOF >> $__filename
#
CXX         = \${CC}
CXXFLAGS    = \${CXXFLAGS} -I\\\${CUDA_PATH}/include -std=c++11 -fopenmp
GPUVER      = \${GPUVER}
NVCC        = \${NVCC}
NVFLAGS     = \${NVFLAGS}
EOF
    fi
    if [ "$__WARNALL" = "on" ] ; then
        cat <<EOF >> $__filename
#
FCLOGPIPE   = 2>&1 | tee \\\$(notdir \\\$<).warn
EOF
    fi
    if [ "$with_gcc" != "__DONTUSE__" ] ; then
        cat <<EOF >> $__filename
#
FYPPFLAGS   = -n --line-marker-format=gfortran5
EOF
    fi
    # replace variable values in output file using eval
    local __TMPL=$(cat $__filename)
    eval "printf \"${__TMPL}\n\"" > $__filename
    # pass this to parsers to replace all of the IF_XYZ statements
    "${SCRIPTDIR}/parse_if.py" -i -f "${__filename}" $__flags
    echo "Wrote ${INSTALLDIR}/arch/$__filename"
}

rm -f ${INSTALLDIR}/arch/local*
# normal production arch files
gen_arch_file "local.ssmp"
gen_arch_file "local.sdbg" DEBUG
arch_vers="ssmp sdbg"

if [ "$MPI_MODE" != no ] ; then
    gen_arch_file "local.psmp"           MPI
    gen_arch_file "local.pdbg"           MPI DEBUG
    gen_arch_file "local_warn.psmp"      MPI WARNALL
    gen_arch_file "local_coverage.pdbg"  MPI COVERAGE
    arch_vers="${arch_vers} psmp pdbg"
fi

# cuda enabled arch files
if [ "$ENABLE_CUDA" = __TRUE__ ] ; then
    gen_arch_file "local_cuda.ssmp"              CUDA
    gen_arch_file "local_cuda.sdbg"              CUDA DEBUG
    if [ "$MPI_MODE" != no ] ; then
        gen_arch_file "local_cuda.psmp"          CUDA MPI
        gen_arch_file "local_cuda.pdbg"          CUDA MPI DEBUG
        gen_arch_file "local_cuda_warn.psmp"     CUDA MPI WARNALL
        gen_arch_file "local_coverage_cuda.pdbg" CUDA MPI COVERAGE
    fi
fi

# valgrind enabled arch files
if [ "$ENABLE_VALGRIND" = __TRUE__ ] ; then
    gen_arch_file "local_valgrind.ssmp"          VALGRIND
    if [ "$MPI_MODE" != no ] ; then
        gen_arch_file "local_valgrind.psmp"      VALGRIND MPI
    fi
fi

cd "${ROOTDIR}"

# -------------------------
# print out user instructions
# -------------------------

cat <<EOF
========================== usage =========================
Done!
Now copy:
  cp ${INSTALLDIR}/arch/* to the cp2k/arch/ directory
To use the installed tools and libraries and cp2k version
compiled with it you will first need to execute at the prompt:
  source ${SETUPFILE}
To build CP2K you should change directory:
  cd cp2k/
  make -j ${NPROCS} ARCH=local VERSION="${arch_vers}"

arch files for GPU enabled CUDA versions are named "local_cuda.*"
arch files for valgrind versions are named "local_valgrind.*"
arch files for coverage versions are named "local_coverage.*"

Note that these pre-built arch files are for the GNU compiler, users have to adapt them for other compilers.
It is possible to use the provided CP2K arch files as guidance.
EOF

#EOF
