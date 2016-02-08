#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh

with_mkl=${1:-__INSTALL__}

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
        if ! [ -z "MKLROOT" ] ; then
            echo "MKLROOT is found to be $MKLROOT"
        else
            report_error ${LINENO} "Cannot find env variable MKLROOT, the script relies on it being set. Please check in MKL installation and use --with-mkl=<location> to pass the path to MKL root directory to this script."
            exit 1
        fi
        check_lib -lm
        check_lib -ldl
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking MKL to user paths ===================="
        check_dir "$with_mkl"
        MKLROOT="$with_mkl"
        ;;
esac
if [ "$with_mkl" != "__DONTUSE__" ] ; then
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
    mkl_required_libs="libmkl_gf_lp64.a libmkl_core.a libmkl_sequential.a"
    for ii in $mkl_required_libs ; do
        if [ ! -f "$mkl_lib_dir/${ii}" ] ; then
            report_error $LINENO "missing MKL library ${ii}"
            exit 1
        fi
    done
    # set the correct lib flags from  MLK link adviser
    MKL_LIBS="-Wl,--start-group ${mkl_lib_dir}/libmkl_gf_lp64.a ${mkl_lib_dir}/libmkl_core.a ${mkl_lib_dir}/libmkl_sequential.a"
    # check optional libraries
    if [ $MPI_MODE != no ] ; then
        enable_mkl_scalapack="__TRUE__"
        mkl_optional_libs="libmkl_scalapack_lp64.a"
        case $MPI_MODE in
            mpich)
                mkl_optional_libs="$mkl_optional_libs libmkl_blacs_lp64.a"
                mkl_blacs_lib="libmkl_blacs_lp64.a"
                ;;
            openmpi)
                mkl_optional_libs="$mkl_optional_libs libmkl_blacs_openmpi_lp64.a"
                mkl_blacs_lib="libmkl_blacs_openmpi_lp64.a"
                ;;
            *)
                enable_mkl_scalapack="__FALSE__"
                ;;
        esac
        for ii in $mkl_optional_libs ; do
            if ! [ -f "${mkl_lib_dir}/${ii}" ] ; then
                enable_mkl_scalapack="__FALSE__"
            fi
        done
        if [ $enable_mkl_scalapack = "__TRUE__" ] ; then
            echo "Using MKL provided ScaLAPACK and BLACS"
            MKL_LIBS="${mkl_lib_dir}/libmkl_scalapack_lp64.a ${MKL_LIBS} ${mkl_lib_dir}/${mkl_blacs_lib}"
        fi
    else
        enable_mkl_scalapack="__FALSE__"
    fi
    MKL_LIBS="${MKL_LIBS} -Wl,--end-group -lpthread -lm -ldl"
    MKL_CFLAGS="${MKL_CFLAGS} -I${MKLROOT}/include"

    # write setup files
    cat <<EOF > "${BUILDDIR}/setup_mkl"
export MKLROOT="${MKLROOT}"
EOF
    cat "${BUILDDIR}/setup_mkl" >> ${SETUPFILE}
    cat <<EOF >> "${BUILDDIR}/setup_mkl"
export MKL_CFLAGS="${MKL_CFLAGS}"
export MKL_LIBS="${MKL_LIBS}"
export FAST_MATH_CFLAGS="\${FAST_MATH_CFLAGS} ${MKL_CFLAGS}"
export FAST_MATH_LIBS="\${FAST_MATH_LIBS} ${MKL_LIBS}"
EOF
    if [ $enable_mkl_scalapack = "__TRUE__" ] ; then
        cat <<EOF >> "${BUILDDIR}/setup_mkl"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__SCALAPACK|)"
with_scalapack="__DONTUSE__"
EOF
    fi
fi
cd "${ROOTDIR}"
