#!/bin/sh -l
# author: Tiziano Müller
# SPDX-License-Identifier: MIT

set -o errexit
set -o nounset

PYTHON=$(command -v python3 || true)
if [ -z "${PYTHON}" ] ; then
    echo "ERROR: the python3 executable could not be found, but is required for the complete build process"
    exit 1
fi

CURL="$(command -v curl || true)"
if [ -z "${CURL}" ] ; then
    echo "ERROR: curl is missing, but required by Spack"
    exit 1
fi

GIT="$(command -v git || true)"
if [ -z "${GIT}" ] ; then
    echo "ERROR: the git command is missing, but required by Spack"
    exit 1
fi

SCRIPTDIR=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
CP2KDIR="${SCRIPTDIR}/../.."
SPACKDIR="${SCRIPTDIR}/spack"
SPACKENVDIR="${SCRIPTDIR}/envs"

AVAILABLE_VARIANTS="cosma elpa libint libxc pexsi plumed sirius"
VARIANTS_REQUIRING_MPI="cosma elpa pexsi sirius"

with_cosma="__INSTALL__"
with_elpa="__INSTALL__"
with_libint="__INSTALL__"
with_libxc="__INSTALL__"
with_libint="__INSTALL__"
with_pexsi="__DONTUSE__"
with_plumed="__DONTUSE__"
with_sirius="__INSTALL__"

mpi_provider="openmpi"
lapack_provider="openblas"
fftw_provider="fftw"
scalapack_provider="netlib-scalapack"

smm_driver="libxsmm"

libint_lmax=5

nprocs=$(nproc 2>/dev/null || echo 1)

no_check_certificates=0

parse_with_arg() {
    # shellcheck disable=SC2039
    # ... non-posix but almost all shells support it
    local input_var="${1#--with*=}"
    case "${input_var}" in
        "$1"|install)
            # lone --with-<kw> treated as --with-<kw>=install
            echo "__INSTALL__"
            ;;
        system)
            echo "__SYSTEM__"
            ;;
        no)
            echo "__DONTUSE__"
            ;;
        *)
            echo "${input_var}"
            ;;
    esac
}

usage() {
    cat << eof
install_cp2k_toolchain.sh: [args]

    Install Spack environments with CP2K dependencies and generate corresponding CP2K arch files

    -h, --help:                              Show this help screen
    -j NPROC, --jobs NPROC                   Specifies the number of parallel build jobs.
    --no-check-certificate                   Tell Spack to not validate HTTPS connection certificates.
    --install-all                            Install all optional dependencies for CP2K, with mpich as MPI provider.
    --mpi-mode=<mpich|openmpi|intelmpi|no>   Select the MPI provider, default: openmpi.
    --math-mode<cray|mkl|openblas|reflapack> Select math providers, default: openblas (with fftw, netlib-scalapack)
    --with-mpich/openmpi/intelmpi            Specify where to find an MPI provider (also switches the MPI provider)
    --verbose                                Enable verbose logging
eof
}

usage_error() {
    echo "ERROR: $1"
    usage
    exit 1
}

configure_spack_env() {
    # shellcheck disable=SC2039
    local env="$1"
    shift

    # add all the other settings
    while [ "$#" -ge 1 ]; do
        "${SPACK}" --env-dir "${SPACKENVDIR}/${env}" config add "$1"
        shift
    done
}

while [ "$#" -ge 1 ] ; do
    case "$1" in
        -j|--jobs)
            nprocs=$1
            shift
            ;;
        --no-check-certificate)
            no_check_certificates=1
            ;;
        --install-all)
            # enable (almost) all variants
            for variant in ${AVAILABLE_VARIANTS} ; do
                eval "with_${variant}=__INSTALL__"
            done
            mpi_provider="mpich"
            ;;
        --mpi-mode=*)
            case "${1#*=}" in
                mpich)
                    mpi_provider="mpich"
                    ;;
                openmpi)
                    mpi_provider="openmpi"
                    ;;
                intelmpi)
                    mpi_provider="intel-mpi"
                    ;;
                no)
                    mpi_provider="none"
                    ;;
                *)
                    usage_error "'$1' is not supported"
                    ;;
            esac
            break
            ;;
        --math-mode=*)
            case "${1#*=}" in
                cray)
                    lapack_provider="cray-libsci"
                    fftw_provider="cray-libsci"
                    scalapack_provider="cray-libsci"
                    ;;
                mkl)
                    lapack_provider="intel-mkl"
                    fftw_provider="intel-mkl"
                    scalapack_provider="intel-mkl"
                    ;;
                openblas)
                    lapack_provider="openblas"
                    ;;
                reflapack)
                    lapack_provider="netlib-lapack"
                    ;;
                *)
                    usage_error "'$1' is not supported"
            esac
            break
            ;;
        #--gpu-ver=*)
        #    case "${1#*=}" in
        #        K20X)
        #            gpuver="K20X"
        #            ;;
        #        K40)
        #            gpuver="K40"
        #            ;;
        #        K80)
        #            gpuver="K80"
        #            ;;
        #        P100)
        #            gpuver="P100"
        #            ;;
        #        V100)
        #            gpuver="V100"
        #            ;;
        #        no)
        #            gpuver="no"
        #            ;;
        #        *)
        #            usage_error "'$1' is not supported"
        #    esac
        #    ;;
        --libint-lmax=*)
            libint_lmax="${1#*=}"
            ;;
        --with-mpich*)
            with_mpich=$(parse_with_arg "$1")
            if [ "${with_mpich}" != __DONTUSE__ ] ; then
                mpi_provider="mpich"
            fi
            ;;
        --with-openmpi*)
            with_openmpi=$(parse_with_arg "$1")
            if [ "${with_openmpi}" = __INSTALL__ ] ; then
                mpi_provider="openmpi"
            elif [ "${with_openmpi}" = __DONTUSE__ ] ; then
                if [ "${mpi_provider}" = none ] || [ "${mpi_provider}" = openmpi ] ; then
                    echo "ERROR: an MPI library is required to build CP2K since you didn't pass --mpi-mode=no"
                    mpi_provider="openmpi"
                fi
            else
                echo "ERROR: explicitly requesting a system library or setting a path is not yet supported (for $1)"
                exit 1
            fi
            ;;
        --with-intelmpi*)
            with_intelmpi=$(parse_with_arg "$1")
            if [ "${with_intelmpi}" = __INSTALL__ ] ; then
                mpi_provider="intel-mpi"
            elif [ "${with_intelmpi}" = __DONTUSE__ ] ; then
                if [ "${mpi_provider}" = "intel-mpi" ] ; then
                    mpi_provider="openmpi"
                fi
            else
                echo "ERROR: explicitly requesting a system library or setting a path is not yet supported (for $1)"
                exit 1
            fi
            ;;
        --with-fftw*)
            with_fftw=$(parse_with_arg "$1")
            if [ "${with_fftw}" = __INSTALL__ ] ; then
                fftw_provider="fftw"
            elif [ "${with_fftw}" = __DONTUSE__ ] ; then
                if [ "${fftw_provider}" = "fftw" ] ; then
                    echo "ERROR: an FFTW library is required to build CP2K, please use --math-mode to switch the provider."
                    exit 1
                fi
            else
                echo "ERROR: explicitly requesting a system library or setting a path is not yet supported (for $1)"
                exit 1
            fi
            ;;
        --with-reflapack*)
            with_reflapack=$(parse_with_arg "$1")
            if [ "${with_reflapack}" = __INSTALL__ ] ; then
                lapack_provider="netlib-lapack"
            elif [ "${with_reflapack}" = __DONTUSE__ ] ; then
                if [ "${lapack_provider}" = "netlib-lapack" ] ; then
                    echo "ERROR: a BLAS/LAPACK implementation is required to build CP2K, please use --math-mode"
                    echo "       or one of the --with-mkl/openblas options to switch the provider."
                    exit 1
                fi
            else
                echo "ERROR: explicitly requesting a system library or setting a path is not yet supported (for $1)"
                exit 1
            fi
            ;;
        --with-mkl*)
            with_mkl=$(parse_with_arg "$1")
            if [ "${with_mkl}" = __INSTALL__ ] ; then
                lapack_provider="intel-mkl"
                scalapack_provider="intel-mkl"
                fftw_provider="intel-mkl"
            elif [ "${with_mkl}" = __DONTUSE__ ] ; then
                if [ "${lapack_provider}" = "intel-mkl" ] || [ "${fftw_provider}" = "intel-mkl" ] || [ "${scalapack_provider}" = "intel-mkl" ] ; then
                    echo "ERROR: please use --math-mode to switch to a different math provider (MKL is specified for at least one of LAPACK, FFTW, ScaLAPACK)."
                    exit 1
                fi
            else
                echo "ERROR: explicitly requesting a system library or setting a path is not yet supported (for $1)"
                exit 1
            fi
            ;;
        --with-openblas*)
            with_openblas=$(parse_with_arg "$1")
            if [ "${with_openblas}" = __INSTALL__ ] ; then
                lapack_provider="openblas"
            elif [ "${with_openblas}" = __DONTUSE__ ] ; then
                if [ "${lapack_provider}" = "openblas" ] ; then
                    echo "ERROR: a BLAS/LAPACK implementation is required to build CP2K, please use --math-mode"
                    echo "       or one of the --with-mkl/reflapack options to switch the provider."
                    exit 1
                fi
            else
                echo "ERROR: explicitly requesting a system library or setting a path is not yet supported (for $1)"
                exit 1
            fi
            ;;
        --with-scalapack*)
            with_scalapack=$(parse_with_arg "$1")
            if [ "${with_scalapack}" = __INSTALL__ ] ; then
                scalapack_provider="netlib-scalapack"
            elif [ "${with_scalapack}" = __DONTUSE__ ] ; then
                if [ "${scalapack_provider}" = "netlib-scalapack" ] && [ "${mpi_provider}" != "none" ] ; then
                    echo "ERROR: a ScaLAPACK implementation is required to build CP2K with MPI, use --mpi-mode=no"
                    echo "       to disable MPI or the --with-mkl option to switch the provider."
                    exit 1
                fi
            else
                echo "ERROR: explicitly requesting a system library or setting a path is not yet supported (for $1)"
                exit 1
            fi
            ;;
        --with-libxsmm*)
            with_libxsmm=$(parse_with_arg "$1")
            if [ "${with_libxsmm}" = "__DONTUSE__" ] ; then
                smm_driver="libxsmm"
            elif [ "${with_libxsmm}" = "__INSTALL__" ] ; then
                smm_driver="blas"
            else
                echo "ERROR: explicitly requesting a system library or setting a path is not yet supported (for $1)"
                exit 1
            fi
            ;;
        --with-libint*)
            # shellcheck disable=SC2034
            with_libint=$(parse_with_arg "$1")
            ;;
        --with-libxc*)
            # shellcheck disable=SC2034
            with_libxc=$(parse_with_arg "$1")
            ;;
        --with-elpa*)
            # shellcheck disable=SC2034
            with_elpa=$(parse_with_arg "$1")
            ;;
        --with-pexsi*)
            # shellcheck disable=SC2034
            with_pexsi=$(parse_with_arg "$1")
            ;;
        --with-plumed*)
            # shellcheck disable=SC2034
            with_plumed=$(parse_with_arg "$1")
            ;;
        --with-sirius*)
            # shellcheck disable=SC2034
            with_sirius=$(parse_with_arg "$1")
            ;;
        --with-cosma*)
            # shellcheck disable=SC2034
            with_cosma=$(parse_with_arg "$1")
            ;;
        --verbose)
            set -x
            ;;
        -h|--help|-*)
            usage
            exit 0
            ;;
        *)
            break
            ;;
        # --with-gcc*)   # TODO: bootstrapping the compiler (and possibly also the binutils) needs more work
        #--with-cmake*)  # TODO: implement ignoring installed CMake and manually setting up an external
        #                  TODO: separate tools need additional installation steps:
        # --with-valgrind*)
        #--enable-cuda*)
        #--enable-cray*)
        #--with-spglib*)  # not a variant for CP2K in spack
        #--with-quip*)  # not available in Spack
    esac
    shift
done

cp2k_variants="lmax=${libint_lmax} smm=${smm_driver}"
cp2k_variants_mpi_only=""

for variant in ${AVAILABLE_VARIANTS} ; do
    val=$(eval "echo \${with_${variant}}")

    if [ "${val}" = "__INSTALL__" ] ; then
        prefix="+"
    elif [ "${val}" = "__DONTUSE__" ] ; then
        prefix="~"
    else
        echo "ERROR: explicitly requesting a system library or setting a path is not yet supported (for --with-${variant})"
        exit 1
    fi

    if [ "${VARIANTS_REQUIRING_MPI#*${variant}}" = "${VARIANTS_REQUIRING_MPI}" ] ; then
        cp2k_variants="${cp2k_variants} ${prefix}${variant}"
    else
        cp2k_variants_mpi_only="${cp2k_variants_mpi_only} ${prefix}${variant}"
    fi
done

if [ "${lapack_provider}" = "netlib-lapack" ] ; then
    blas_provider="netlib-xblas"
else
    blas_provider="${lapack_provider}"
fi

if [ -d "${SPACKDIR}" ] ; then
    echo "Found a Spack installation, updating..."
    # -C was introduced in git 1.8.5 (Q4 2013)
    git -C "${SPACKDIR}" pull
else
    git clone --depth=1 https://github.com/spack/spack.git "${SPACKDIR}"
fi

SPACK=$(PATH="${SPACKDIR}/bin:${PATH}" command -v spack)
SPACK_ARGS=""

if [ ${no_check_certificates} -eq 1 ] ; then
    SPACK_ARGS="${SPACK_ARGS} --insecure"
fi

echo "Let 'spack' detect system-packages to be used..."
# do it on the site-scope to avoid touching a users Spack config,
# also, if there is a user/system-provided Spack, do not run this
"${SPACK}" external find --scope site

 # workaround Spack issue https://github.com/spack/spack/issues/4635
"${SPACK}" config --scope site rm packages:cmake

# TODO: if there is already a Spack environment/command avvailable, we should setup a chained installation instead

echo "Installing non-MPI environment..."
if [ ! -d "${SPACKENVDIR}/ssmp" ] ; then
    "${SPACK}" env create -d "${SPACKENVDIR}/ssmp"
else
    echo "Existing environment found, reconfiguring..."
fi
configure_spack_env "ssmp" \
    "config:build_language:C.UTF-8" \
    "packages:all:providers:fftw-api:[${fftw_provider}]" \
    "packages:all:providers:blas:[${blas_provider}]" \
    "packages:all:providers:lapack:[${lapack_provider}]" \
    "packages:cp2k:variants:${cp2k_variants} ~mpi" \
    "packages:fftw:variants:+openmp" \
    "packages:openblas:variants:threads=openmp"

if [ "${mpi_provider}" = none ] ; then
    # avoid installing MPI also in the non-MPI environment if fftw is pulled-in
    # note: fftw builds 2 sets of libs with +mpi
    configure_spack_env "ssmp" \
        "packages:fftw:variants:+openmp ~mpi"
fi

echo "Installing packages..."
"${SPACK}" --env-dir "${SPACKENVDIR}/ssmp" dev-build -j "${nprocs}" -d "${CP2KDIR}" --quiet --until edit cp2k@master
echo "Generating environment setup file..."
"${SPACK}" --env-dir "${SPACKENVDIR}/ssmp" build-env --dump "${SPACKENVDIR}/ssmp/setup" cp2k@master --

if [ "${mpi_provider}" != none ] ; then
    echo "Installing MPI environment..."
    if [ ! -d "${SPACKENVDIR}/psmp" ] ; then
        "${SPACK}" env create -d "${SPACKENVDIR}/psmp"
    else
        echo "Existing environment found, reconfiguring..."
    fi
    configure_spack_env "psmp" \
        "config:build_language:C.UTF-8" \
        "packages:all:providers:fftw-api:[${fftw_provider}]" \
        "packages:all:providers:blas:[${lapack_provider}]" \
        "packages:all:providers:lapack:[${blas_provider}]" \
        "packages:all:providers:mpi:[${mpi_provider}]" \
        "packages:all:providers:scalapack:[${scalapack_provider}]" \
        "packages:cp2k:variants:${cp2k_variants} ${cp2k_variants_mpi_only ${cp2k_variants_mpi_only}}" \
        "packages:fftw:variants:+openmp" \
        "packages:openblas:variants:threads=openmp"

    echo "Installing packages..."
    "${SPACK}" --env-dir "${SPACKENVDIR}/psmp" dev-build -j "${nprocs}" -d "${CP2KDIR}" --quiet --until edit cp2k@master
    echo "Generating environment setup file..."
    "${SPACK}" --env-dir "${SPACKENVDIR}/psmp" build-env --dump "${SPACKENVDIR}/psmp/setup" cp2k@master --
fi

# shellcheck disable=SC2016
ARCH=$("${SPACK}" arch)-$("${SPACK}" --env-dir "${SPACKENVDIR}/ssmp" build-env cp2k@master -- sh -c 'echo ${SPACK_COMPILER_SPEC%%@*}')

echo "Building environments using Spack finished successfully..."
echo "To build the single-node implementation of CP2K, run:"
echo ""
echo "    source \"${SPACKENVDIR}/ssmp/setup\""
echo "    make -j${nprocs} ARCH=${ARCH} VERSION=ssmp"
if [ "${mpi_provider}" != none ] ; then
    echo ""
    echo "To build the multi-node (MPI) implementation of CP2K, run:"
    echo ""
    echo "    source \"${SPACKENVDIR}/psmp/setup\""
    echo "    make -j${nprocs} ARCH=${ARCH} VERSION=psmp"
fi
