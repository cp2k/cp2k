#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

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
CC_arch="IF_MPI(${MPICC}|${CC})"
CXX_arch="IF_MPI(${MPICXX}|${CXX})"
FC_arch="IF_MPI(${MPIFC}|${FC})"
LD_arch="IF_MPI(${MPIFC}|${FC})"

# we always want good line information and backtraces
if [ "${with_intel}" != "__DONTUSE__" ]; then
  if [ "${TARGET_CPU}" = "native" ]; then
    BASEFLAGS="-fPIC -fp-model=precise -g -qopenmp -traceback -xHost"
  elif [ "${TARGET_CPU}" = "generic" ]; then
    BASEFLAGS="-fPIC -fp-model=precise -g -mtune=${TARGET_CPU} -qopenmp -traceback"
  else
    BASEFLAGS="-fPIC -fp-model=precise -g -march=${TARGET_CPU} -mtune=${TARGET_CPU} -qopenmp -traceback"
  fi
  OPT_FLAGS="-O2 -funroll-loops"
  LDFLAGS_C="-nofor-main"
else
  if [ "${TARGET_CPU}" = "generic" ]; then
    BASEFLAGS="-fno-omit-frame-pointer -fopenmp -g -mtune=${TARGET_CPU} IF_ASAN(-fsanitize=address|)"
  else
    BASEFLAGS="-fno-omit-frame-pointer -fopenmp -g -march=${TARGET_CPU} -mtune=${TARGET_CPU} IF_ASAN(-fsanitize=address|)"
  fi
  OPT_FLAGS="-O3 -funroll-loops"
  LDFLAGS_C=""
fi

NOOPT_FLAGS="-O1"

# those flags that do not influence code generation are used always, the others if debug
if [ "${with_intel}" != "__DONTUSE__" ]; then
  FCDEB_FLAGS=""
  FCDEB_FLAGS_DEBUG=""
else
  FCDEB_FLAGS="-fbacktrace -ffree-form -fimplicit-none -std=f2008"
  FCDEB_FLAGS_DEBUG="-fsanitize=leak -fcheck=all,no-array-temps -ffpe-trap=invalid,zero,overflow -finit-derived -finit-real=snan -finit-integer=-42 -Werror=realloc-lhs -finline-matmul-limit=0"
fi

# code coverage generation flags
COVERAGE_FLAGS="-O1 -coverage -fkeep-static-functions"
COVERAGE_DFLAGS="-D__NO_ABORT"

# profile based optimization, see https://www.cp2k.org/howto:pgo
PROFOPT_FLAGS="\$(PROFOPT)"

# special flags for gfortran
# https://gcc.gnu.org/onlinedocs/gfortran/Error-and-Warning-Options.html
# we error out for these warnings (-Werror=uninitialized -Wno-maybe-uninitialized -> error on variables that must be used uninitialized)
WFLAGS_ERROR="-Werror=aliasing -Werror=ampersand -Werror=c-binding-type -Werror=intrinsic-shadow -Werror=intrinsics-std -Werror=line-truncation -Werror=tabs -Werror=target-lifetime -Werror=underflow -Werror=unused-but-set-variable -Werror=unused-variable -Werror=unused-dummy-argument -Werror=unused-parameter -Werror=unused-label -Werror=conversion -Werror=zerotrip -Wno-maybe-uninitialized"
# we just warn for those (that eventually might be promoted to WFLAGSERROR). It is useless to put something here with 100s of warnings.
WFLAGS_WARN="-Wuninitialized -Wuse-without-only"
# while here we collect all other warnings, some we'll ignore
# TODO: -Wpedantic with -std2008 requires us to drop the old MPI-90 interface entirely and SIRIUS (-Wuninitialized) and DBCSR (default initializers) to initialize their types
WFLAGS_WARNALL="-Wno-pedantic -Wall -Wextra -Wsurprising -Warray-temporaries -Wcharacter-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wreal-q-constant -Walign-commons -Wfunction-elimination -Wrealloc-lhs -Wcompare-reals -Wzerotrip"

# IEEE_EXCEPTIONS dependency
IEEE_EXCEPTIONS_DFLAGS="-D__HAS_IEEE_EXCEPTIONS"

# check all of the above flags, filter out incompatible flags for the
# current version of gcc in use
if [ "${with_intel}" == "__DONTUSE__" ]; then
  OPT_FLAGS=$(allowed_gfortran_flags $OPT_FLAGS)
  NOOPT_FLAGS=$(allowed_gfortran_flags $NOOPT_FLAGS)
  FCDEB_FLAGS=$(allowed_gfortran_flags $FCDEB_FLAGS)
  FCDEB_FLAGS_DEBUG=$(allowed_gfortran_flags $FCDEB_FLAGS_DEBUG)
  COVERAGE_FLAGS=$(allowed_gfortran_flags $COVERAGE_FLAGS)
  WFLAGS_ERROR=$(allowed_gfortran_flags $WFLAGS_ERROR)
  WFLAGS_WARN=$(allowed_gfortran_flags $WFLAGS_WARN)
  WFLAGS_WARNALL=$(allowed_gfortran_flags $WFLAGS_WARNALL)
else
  WFLAGS_ERROR=""
  WFLAGS_WARN=""
  WFLAGS_WARNALL=""
fi

# check if ieee_exeptions module is available for the current version
# of gfortran being used
if ! (check_gfortran_module ieee_exceptions); then
  IEEE_EXCEPTIONS_DFLAGS=""
fi

# concatenate the above flags into WFLAGS, FCDEBFLAGS, DFLAGS and
# finally into FCFLAGS and CFLAGS
WFLAGS="$WFLAGS_ERROR $WFLAGS_WARN IF_WARNALL(${WFLAGS_WARNALL}|)"
FCDEBFLAGS="$FCDEB_FLAGS IF_DEBUG($FCDEB_FLAGS_DEBUG|)"
DFLAGS="${CP_DFLAGS} IF_DEBUG($IEEE_EXCEPTIONS_DFLAGS -D__CHECK_DIAG|) IF_COVERAGE($COVERAGE_DFLAGS|)"
# language independent flags
G_CFLAGS="$BASEFLAGS"
G_CFLAGS="$G_CFLAGS IF_COVERAGE($COVERAGE_FLAGS|IF_DEBUG($NOOPT_FLAGS|$OPT_FLAGS))"
G_CFLAGS="$G_CFLAGS IF_DEBUG(|$PROFOPT_FLAGS)"
G_CFLAGS="$G_CFLAGS $CP_CFLAGS"
if [ "${with_intel}" == "__DONTUSE__" ]; then
  # FCFLAGS, for gfortran
  FCFLAGS="$G_CFLAGS \$(FCDEBFLAGS) \$(WFLAGS) \$(DFLAGS)"
  FCFLAGS+=" IF_MPI($(allowed_gfortran_flags "-fallow-argument-mismatch")|)"
else
  FCFLAGS="$G_CFLAGS \$(FCDEBFLAGS) \$(WFLAGS) \$(DFLAGS)"
fi
# CFLAGS, special flags for gcc

# TODO: Remove -Wno-vla-parameter after upgrade to gcc 11.3.
# https://gcc.gnu.org/bugzilla//show_bug.cgi?id=101289
if [ "${with_intel}" == "__DONTUSE__" ]; then
  CFLAGS="$G_CFLAGS -std=c11 -Wall -Wextra -Werror -Wno-vla-parameter -Wno-deprecated-declarations \$(DFLAGS)"
else
  if [ "${with_ifx}" == "no" ]; then
    CC_arch+=" IF_MPI(-cc=${I_MPI_CC}|)"
    CXX_arch+=" IF_MPI(-cxx=${I_MPI_CXX}|)"
    FC_arch+=" IF_MPI(-fc=${I_MPI_FC}|)"
    LD_arch+=" IF_MPI(-fc=${I_MPI_FC}|)"
  fi
  CFLAGS="${G_CFLAGS} -std=c11 -Wall \$(DFLAGS)"
  CXXFLAGS="${G_CFLAGS} -std=c++14 -Wall \$(DFLAGS)"
  FCFLAGS="${FCFLAGS} -diag-disable=8291 -diag-disable=8293 -fpp -fpscomp logicals -free"
  # Suppress warnings and add include path to omp_lib.mod explicitly.
  # No clue why the Intel oneAPI setup script does not include that path (bug?)
  FCFLAGS="${FCFLAGS} -diag-disable=10448 -I/opt/intel/oneapi/2024.1/opt/compiler/include/intel64"
fi

# Linker flags
# About --whole-archive see: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=52590
STATIC_FLAGS="-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive"
# Get unfortunately ignored: -static-libgcc -static-libstdc++ -static-libgfortran
LDFLAGS="IF_STATIC(${STATIC_FLAGS}|) \$(FCFLAGS) ${CP_LDFLAGS}"

# Library flags
# add standard libs
LIBS="${CP_LIBS} -lstdc++"

if [ "${with_intel}" == "__DONTUSE__" ]; then
  CXXFLAGS+=" --std=c++14 \$(DFLAGS) -Wno-deprecated-declarations"
else
  CXXFLAGS+=" --std=c++14 \$(DFLAGS)"
fi
# CUDA handling
if [ "${ENABLE_CUDA}" = __TRUE__ ] && [ "${GPUVER}" != no ]; then
  CUDA_LIBS="-lcudart -lnvrtc -lcuda -lcufft -lcublasLt -lcublas -lrt IF_DEBUG(-lnvToolsExt|)"
  CUDA_DFLAGS="-D__OFFLOAD_CUDA -D__DBCSR_ACC IF_DEBUG(-D__OFFLOAD_PROFILING|)"
  if [ "${with_cusolvermp}" != "__DONTUSE__" ]; then
    CUDA_LIBS+=" -lcusolverMp -lcusolver -lcal -lnvidia-ml"
    CUDA_DFLAGS+=" -D__CUSOLVERMP"
  fi
  LIBS="${LIBS} IF_CUDA(${CUDA_LIBS}|)"
  DFLAGS="IF_CUDA(${CUDA_DFLAGS}|) ${DFLAGS}"
  NVFLAGS="-g -arch sm_${ARCH_NUM} -O3 -allow-unsupported-compiler -Xcompiler='-fopenmp -Wall -Wextra -Werror' --std=c++11 \$(DFLAGS)"
  check_command nvcc "cuda"
  check_lib -lcudart "cuda"
  check_lib -lnvrtc "cuda"
  check_lib -lcuda "cuda"
  check_lib -lcufft "cuda"
  check_lib -lcublasLt "cuda"
  check_lib -lcublas "cuda"

  # Set include flags
  CUDA_FLAGS=""
  add_include_from_paths CUDA_FLAGS "cuda.h" $INCLUDE_PATHS
  NVFLAGS+=" ${CUDA_FLAGS}"
  NVCC_TOPDIR="$(dirname $(command -v nvcc))/.."
  CUDA_PATH="${CUDA_PATH:-${CUDA_HOME:-${NVCC_TOPDIR:-/CUDA_HOME-notset}}}"
  CFLAGS+=" IF_CUDA(-I${CUDA_PATH}/include|)"
  CXXFLAGS+=" IF_CUDA(-I${CUDA_PATH}/include|)"

  # Set LD-flags
  CUDA_LDFLAGS=""
  add_lib_from_paths CUDA_LDFLAGS "libcudart.*" $LIB_PATHS
  add_lib_from_paths CUDA_LDFLAGS "libnvrtc.*" $LIB_PATHS
  add_lib_from_paths CUDA_LDFLAGS "libcuda.*" $LIB_PATHS
  add_lib_from_paths CUDA_LDFLAGS "libcufft.*" $LIB_PATHS
  add_lib_from_paths CUDA_LDFLAGS "libcublas.*" $LIB_PATHS
  export CUDA_LDFLAGS="${CUDA_LDFLAGS}"
  LDFLAGS+=" IF_CUDA(${CUDA_LDFLAGS}|)"
fi

# HIP handling
if [ "${ENABLE_HIP}" = __TRUE__ ] && [ "${GPUVER}" != no ]; then
  check_command hipcc "hip"
  check_lib -lhipblas "hip"
  add_lib_from_paths HIP_LDFLAGS "libhipblas.*" $LIB_PATHS
  check_lib -lhipfft "hip"
  add_lib_from_paths HIP_LDFLAGS "libhipfft.*" $LIB_PATHS

  HIP_INCLUDES="-I${ROCM_PATH}/include"
  case "${GPUVER}" in
    Mi50)
      check_lib -lamdhip64 "hip"
      add_lib_from_paths HIP_LDFLAGS "libamdhip64.*" $LIB_PATHS
      check_lib -lhipfft "hip"
      add_lib_from_paths HIP_LDFLAGS "libhipfft.*" $LIB_PATHS
      check_lib -lrocblas "hip"
      add_lib_from_paths HIP_LDFLAGS "librocblas.*" $LIB_PATHS
      check_lib -lroctx64 "hip"
      add_lib_from_paths HIP_LDFLAGS "libroctx64.*" $LIB_PATHS
      check_lib -lroctracer64 "hip"
      add_lib_from_paths HIP_LDFLAGS "libroctracer64.*" $LIB_PATHS
      HIP_FLAGS+="-fPIE -D__HIP_PLATFORM_AMD__ -g --offload-arch=gfx906 -O3 --std=c++11 -Wall -Wextra -Werror \$(DFLAGS)"
      LIBS+=" IF_HIP(-lamdhip64 -lhipfft -lhipblas -lrocblas IF_DEBUG(-lroctx64 -lroctracer64|)|)"
      DFLAGS+=" IF_HIP(-D__HIP_PLATFORM_AMD__ -D__OFFLOAD_HIP IF_DEBUG(-D__OFFLOAD_PROFILING|)|)  -D__DBCSR_ACC"
      CXXFLAGS+=" -fopenmp -Wall -Wextra -Werror"
      ;;
    Mi100)
      check_lib -lamdhip64 "hip"
      add_lib_from_paths HIP_LDFLAGS "libamdhip64.*" $LIB_PATHS
      check_lib -lhipfft "hip"
      add_lib_from_paths HIP_LDFLAGS "libhipfft.*" $LIB_PATHS
      check_lib -lrocblas "hip"
      add_lib_from_paths HIP_LDFLAGS "librocblas.*" $LIB_PATHS
      check_lib -lroctx64 "hip"
      add_lib_from_paths HIP_LDFLAGS "libroctx64.*" $LIB_PATHS
      check_lib -lroctracer64 "hip"
      add_lib_from_paths HIP_LDFLAGS "libroctracer64.*" $LIB_PATHS
      HIP_FLAGS+="-fPIE -D__HIP_PLATFORM_AMD__ -g --offload-arch=gfx908 -O3 --std=c++11 -Wall -Wextra -Werror \$(DFLAGS)"
      LIBS+=" IF_HIP(-lamdhip64 -lhipfft -lhipblas -lrocblas IF_DEBUG(-lroctx64 -lroctracer64|)|)"
      DFLAGS+=" IF_HIP(-D__HIP_PLATFORM_AMD__ -D__OFFLOAD_HIP IF_DEBUG(-D__OFFLOAD_PROFILING|)|) -D__DBCSR_ACC"
      CXXFLAGS+=" -fopenmp -Wall -Wextra -Werror"
      ;;
    Mi250)
      check_lib -lamdhip64 "hip"
      add_lib_from_paths HIP_LDFLAGS "libamdhip64.*" $LIB_PATHS
      check_lib -lhipfft "hip"
      add_lib_from_paths HIP_LDFLAGS "libhipfft.*" $LIB_PATHS
      check_lib -lrocblas "hip"
      add_lib_from_paths HIP_LDFLAGS "librocblas.*" $LIB_PATHS
      check_lib -lroctx64 "hip"
      add_lib_from_paths HIP_LDFLAGS "libroctx64.*" $LIB_PATHS
      check_lib -lroctracer64 "hip"
      add_lib_from_paths HIP_LDFLAGS "libroctracer64.*" $LIB_PATHS
      HIP_FLAGS+="-fPIE -D__HIP_PLATFORM_AMD__ -g --offload-arch=gfx90a -munsafe-fp-atomics -O3 --std=c++11 -Wall -Wextra -Werror \$(DFLAGS)"
      LIBS+=" IF_HIP(-lamdhip64 -lhipfft -lhipblas -lrocblas IF_DEBUG(-lroctx64 -lroctracer64|)|)"
      DFLAGS+=" IF_HIP(-D__HIP_PLATFORM_AMD__ -D__OFFLOAD_HIP IF_DEBUG(-D__OFFLOAD_PROFILING|)|) -D__DBCSR_ACC"
      CXXFLAGS+=" -fopenmp -Wall -Wextra -Werror"
      ;;
    *)
      check_command nvcc "cuda"
      check_lib -lcudart "cuda"
      check_lib -lnvrtc "cuda"
      check_lib -lcuda "cuda"
      check_lib -lcufft "cuda"
      check_lib -lcublas "cuda"
      DFLAGS+=" IF_HIP(-D__HIP_PLATFORM_NVIDIA__ -D__HIP_PLATFORM_NVCC__ -D__OFFLOAD_HIP |) -D__DBCSR_ACC"
      HIP_FLAGS+=" -g -arch sm_${ARCH_NUM} -O3 -Xcompiler='-fopenmp -Wall -Wextra -Werror' --std=c++11 \$(DFLAGS)"
      add_include_from_paths CUDA_CFLAGS "cuda.h" $INCLUDE_PATHS
      HIP_INCLUDES+=" -I${CUDA_PATH:-${CUDA_HOME:-/CUDA_HOME-notset}}/include"
      # GCC issues lots of warnings for hip/nvidia_detail/hip_runtime_api.h
      CFLAGS+=" -Wno-error ${CUDA_CFLAGS}"
      CXXFLAGS+=" -Wno-error ${CUDA_CFLAGS}"
      # Set LD-flags
      # Multiple definition because of hip/include/hip/nvidia_detail/nvidia_hiprtc.h
      LDFLAGS+=" -Wl,--allow-multiple-definition"
      LIBS+=" -lhipfft -lhipblas -lhipfft -lnvrtc -lcudart -lcufft -lcublas -lcuda"
      add_lib_from_paths HIP_LDFLAGS "libcudart.*" $LIB_PATHS
      add_lib_from_paths HIP_LDFLAGS "libnvrtc.*" $LIB_PATHS
      add_lib_from_paths HIP_LDFLAGS "libcuda.*" $LIB_PATHS
      add_lib_from_paths HIP_LDFLAGS "libcufft.*" $LIB_PATHS
      add_lib_from_paths HIP_LDFLAGS "libcublas.*" $LIB_PATHS
      ;;
  esac

  LDFLAGS+=" ${HIP_LDFLAGS}"
  CFLAGS+=" ${HIP_INCLUDES}"
  CXXFLAGS+=" ${HIP_INCLUDES}"
fi

# OpenCL handling (GPUVER is not a prerequisite)
if [ "${ENABLE_OPENCL}" = __TRUE__ ]; then
  OPENCL_DFLAGS="-D__DBCSR_ACC"
  # avoid duplicating FLAGS
  if [[ "${GPUVER}" == no || ("${ENABLE_CUDA}" != __TRUE__ && "${ENABLE_HIP}" != __TRUE__) ]]; then
    OPENCL_FLAGS="${CFLAGS} ${OPENCL_DFLAGS} ${DFLAGS}"
    DFLAGS="IF_OPENCL(${OPENCL_DFLAGS} ${DFLAGS}|)"
    # Set include flags
    OPENCL_INCLUDES=""
    add_include_from_paths -p OPENCL_INCLUDES "CL" $INCLUDE_PATHS
    if [ -e "${OPENCL_INCLUDES}/CL/cl.h" ]; then
      OPENCL_FLAGS+=" ${OPENCL_INCLUDES}"
    fi
  fi
  # Append OpenCL library to LIBS
  LIBOPENCL=$(ldconfig -p 2> /dev/null | grep -m1 OpenCL | rev | cut -d' ' -f1 | rev)
  if [ -e "${LIBOPENCL}" ]; then
    echo "Found library ${LIBOPENCL}"
    LIBS+=" IF_OPENCL(${LIBOPENCL}|)"
  else
    LIBS+=" IF_OPENCL(-lOpenCL|)"
  fi
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
  local __full_flag_list="MPI DEBUG CUDA WARNALL COVERAGE"
  local __flag=""
  for __flag in $__full_flag_list; do
    eval "local __${__flag}=off"
  done
  for __flag in $__flags; do
    eval "__${__flag}=on"
  done
  # generate initial arch file
  cat $ARCH_FILE_TEMPLATE > $__filename
  # add additional parts
  if [ "$__CUDA" = "on" ]; then
    cat << EOF >> $__filename
#
GPUVER        = \${GPUVER}
OFFLOAD_CC    = \${NVCC}
OFFLOAD_FLAGS = \${NVFLAGS}
OFFLOAD_TARGET = cuda
EOF
  fi

  if [ "$__HIP" = "on" ]; then
    cat << EOF >> $__filename
#
GPUVER        = \${GPUVER}
OFFLOAD_CC    = \${ROCM_PATH}/hip/bin/hipcc
OFFLOAD_FLAGS = \${HIP_FLAGS} \${HIP_INCLUDES}
OFFLOAD_TARGET = hip
EOF
  fi

  if [ "$__OPENCL" = "on" ]; then
    cat << EOF >> $__filename
#
override DBCSR_USE_ACCEL = opencl
EOF
    if [ "${OPENCL_FLAGS}" ]; then
      cat << EOF >> $__filename
OFFLOAD_FLAGS = \${OPENCL_FLAGS}
EOF
    fi
  fi

  if [ "$__WARNALL" = "on" ]; then
    cat << EOF >> $__filename
#
SHELL       := bash
FC          := set -o pipefail && \\\${FC}
CC          := set -o pipefail && \\\${CC}
CXX         := set -o pipefail && \\\${CXX}
LD          := set -o pipefail && \\\${LD}
FCLOGPIPE   = 2>&1 | tee \\\$(notdir \\\$<).warn
EOF
  fi
  if [ "$with_gcc" != "__DONTUSE__" ]; then
    cat << EOF >> $__filename
#
FYPPFLAGS   = -n --line-marker-format=gfortran5
EOF
  fi
  if [ "${with_intel}" != "__DONTUSE__" ]; then
    cat << EOF >> $__filename
#
# Required due to memory leak that occurs if high optimisations are used
mp2_optimize_ri_basis.o: mp2_optimize_ri_basis.F
	\\\$(FC) -c \\\$(subst -O2,-O0,\\\$(FCFLAGS)) \\\$<
# Required due to SEGFAULTS occurring for higher optimisation levels
paw_basis_types.o: paw_basis_types.F
	\\\$(FC) -c \\\$(subst -O2,-O1,\\\$(FCFLAGS)) \\\$<
# Reduce compilation time
hfx_contraction_methods.o: hfx_contraction_methods.F
	\\\$(FC) -c \\\$(subst -O2,-O1,\\\$(FCFLAGS)) \\\$<
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
if [ "${with_intel}" != "__DONTUSE__" ]; then
  gen_arch_file "local.ssmp"
  gen_arch_file "local.sdbg" DEBUG
else
  gen_arch_file "local.ssmp"
  gen_arch_file "local_static.ssmp" STATIC
  gen_arch_file "local.sdbg" DEBUG
  gen_arch_file "local_asan.ssmp" ASAN
  gen_arch_file "local_coverage.sdbg" COVERAGE
fi
arch_vers="ssmp sdbg"

if [ "$MPI_MODE" != no ]; then
  if [ "${with_intel}" != "__DONTUSE__" ]; then
    gen_arch_file "local.psmp" MPI
    gen_arch_file "local.pdbg" MPI DEBUG
  else
    gen_arch_file "local.psmp" MPI
    gen_arch_file "local.pdbg" MPI DEBUG
    gen_arch_file "local_asan.psmp" MPI ASAN
    gen_arch_file "local_static.psmp" MPI STATIC
    gen_arch_file "local_warn.psmp" MPI WARNALL
    gen_arch_file "local_coverage.pdbg" MPI COVERAGE
  fi
  arch_vers="${arch_vers} psmp pdbg"
fi

# opencl enabled arch files
if [ "$ENABLE_OPENCL" = __TRUE__ ]; then
  gen_arch_file "local_opencl.ssmp" OPENCL
  gen_arch_file "local_opencl.sdbg" OPENCL DEBUG
  if [ "$MPI_MODE" != no ]; then
    gen_arch_file "local_opencl.psmp" OPENCL MPI
    gen_arch_file "local_opencl.pdbg" OPENCL MPI DEBUG
    gen_arch_file "local_opencl_warn.psmp" OPENCL MPI WARNALL
    gen_arch_file "local_coverage_opencl.pdbg" OPENCL MPI COVERAGE
  fi
  DBCSR_OPENCL=OPENCL
fi

# cuda enabled arch files
if [ "$ENABLE_CUDA" = __TRUE__ ]; then
  gen_arch_file "local_cuda.ssmp" CUDA ${DBCSR_OPENCL}
  gen_arch_file "local_cuda.sdbg" CUDA ${DBCSR_OPENCL} DEBUG
  if [ "$MPI_MODE" != no ]; then
    gen_arch_file "local_cuda.psmp" CUDA ${DBCSR_OPENCL} MPI
    gen_arch_file "local_cuda.pdbg" CUDA ${DBCSR_OPENCL} MPI DEBUG
    gen_arch_file "local_cuda_warn.psmp" CUDA ${DBCSR_OPENCL} MPI WARNALL
    gen_arch_file "local_coverage_cuda.pdbg" CUDA ${DBCSR_OPENCL} MPI COVERAGE
  fi
fi

# hip enabled arch files
if [ "$ENABLE_HIP" = __TRUE__ ]; then
  gen_arch_file "local_hip.ssmp" HIP
  gen_arch_file "local_hip.sdbg" HIP DEBUG
  if [ "$MPI_MODE" != no ]; then
    gen_arch_file "local_hip.psmp" HIP MPI
    gen_arch_file "local_hip.pdbg" HIP MPI DEBUG
    gen_arch_file "local_hip_warn.psmp" HIP MPI WARNALL
    gen_arch_file "local_coverage_hip.pdbg" HIP MPI COVERAGE
  fi
fi

cd "${ROOTDIR}"

# -------------------------
# print out user instructions
# -------------------------

cat << EOF
========================== usage =========================
Done!
Now copy:
  cp ${INSTALLDIR}/arch/* to the cp2k/arch/ directory
To use the installed tools and libraries and cp2k version
compiled with it you will first need to execute at the prompt:
  source ${SETUPFILE}
To build CP2K you should change directory:
  cd cp2k/
  make -j $(get_nprocs) ARCH=local VERSION="${arch_vers}"

arch files for GPU enabled CUDA versions are named "local_cuda.*"
arch files for GPU enabled HIP versions are named "local_hip.*"
arch files for OpenCL (GPU) versions are named "local_opencl.*"
arch files for coverage versions are named "local_coverage.*"

Note that these pre-built arch files are for the GNU compiler, users have to adapt them for other compilers.
It is possible to use the provided CP2K arch files as guidance.
EOF

#EOF
