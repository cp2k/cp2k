#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_mkl" ] && rm "${BUILDDIR}/setup_mkl"

MKL_CFLAGS=''
MKL_LDFLAGS=''
MKL_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_mkl" in
  __INSTALL__)
    echo "==================== Installing MKL ===================="
    report_error ${LINENO} "To install MKL you should contact your system administrator."
    exit 1
    ;;
  __SYSTEM__)
    echo "==================== Finding MKL from system paths ===================="
    if ! [ -z "${MKLROOT}" ]; then
      echo "MKLROOT is found to be ${MKLROOT}"
    else
      report_error ${LINENO} "Cannot find env variable MKLROOT, the script relies on it being set. Please check in MKL installation and use --with-mkl=<location> to pass the path to MKL root directory to this script."
      exit 1
    fi
    check_lib -lm
    check_lib -ldl
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking MKL to user paths ===================="
    check_dir "$with_mkl"
    MKLROOT="$with_mkl"
    ;;
esac
if [ "$with_mkl" != "__DONTUSE__" ]; then
  case $OPENBLAS_ARCH in
    x86_64)
      mkl_arch_dir=intel64
      MKL_CFLAGS="-m64"
      ;;
    i386)
      mkl_arch_dir=ia32
      MKL_CFLAGS="-m32"
      ;;
    *)
      report_error $LINENO "MKL only supports intel64 (x86_64) and ia32 (i386) at the moment, and your arch obtained from OpenBLAS prebuild is $OPENBLAS_ARCH"
      exit 1
      ;;
  esac
  mkl_lib_dir="${MKLROOT}/lib/${mkl_arch_dir}"
  # check we have required libraries
  mkl_required_libs="libmkl_gf_lp64.so libmkl_sequential.so libmkl_core.so"
  for ii in $mkl_required_libs; do
    if [ ! -f "$mkl_lib_dir/${ii}" ]; then
      report_error $LINENO "missing MKL library ${ii}"
      exit 1
    fi
  done

  case $MPI_MODE in
    intelmpi|mpich)
      mkl_scalapack_lib="-lmkl_scalapack_lp64"
      mkl_blacs_lib="-lmkl_blacs_intelmpi_lp64"
      ;;
    openmpi)
      mkl_scalapack_lib="-lmkl_scalapack_lp64"
      mkl_blacs_lib="-lmkl_blacs_openmpi_lp64"
      ;;
    *)
      echo "Not using MKL provided ScaLAPACK and BLACS"
      mkl_scalapack_lib=""
      mkl_blacs_lib=""
      ;;
  esac

  # set the correct lib flags from MLK link adviser
  MKL_LIBS="-L${mkl_lib_dir} -Wl,-rpath=${mkl_lib_dir} ${mkl_scalapack_lib}"
  MKL_LIBS+=" -Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core"
  MKL_LIBS+=" ${mkl_blacs_lib} -Wl,--end-group -lpthread -lm -ldl"
  # setup_mkl disables using separate FFTW library (see below)
  MKL_CFLAGS="${MKL_CFLAGS} -I${MKLROOT}/include -I${MKLROOT}/include/fftw"

  # write setup files
  cat << EOF > "${BUILDDIR}/setup_mkl"
export MKLROOT="${MKLROOT}"
EOF
  cat "${BUILDDIR}/setup_mkl" >> ${SETUPFILE}
  cat << EOF >> "${BUILDDIR}/setup_mkl"
export MKL_CFLAGS="${MKL_CFLAGS}"
export MKL_LIBS="${MKL_LIBS}"
export FAST_MATH_CFLAGS="\${FAST_MATH_CFLAGS} ${MKL_CFLAGS}"
export FAST_MATH_LIBS="\${FAST_MATH_LIBS} ${MKL_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__MKL -D__FFTW3 IF_COVERAGE(IF_MPI(|-U__FFTW3)|)"
export with_fftw="__DONTUSE__"
EOF
  if [ -n "${mkl_scalapack_lib}" ]; then
    cat << EOF >> "${BUILDDIR}/setup_mkl"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__SCALAPACK|)"
export with_scalapack="__DONTUSE__"
EOF
  fi
fi

load "${BUILDDIR}/setup_mkl"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "mkl"
