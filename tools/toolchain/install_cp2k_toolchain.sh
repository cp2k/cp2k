#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

# +-----------------------------------------------------------------------+
# |  CP2K: A general program to perform molecular dynamics simulations    |
# |  Copyright (C) 2000 - 2016  CP2K developers group                     |
# +-----------------------------------------------------------------------+
#
# *****************************************************************************
#> \brief    This script will compile and install or link existing tools and
#>           libraries CP2K depends on and generate a set of ARCH files which
#>           can be used to compile CP2K
#> \history  Created on Friday, 2016/02/05
#> \author   Lianheng Tong (ltong) lianheng.tong@kcl.ac.uk
# *****************************************************************************

# ------------------------------------------------------------------------
# Work directories and used files
# ------------------------------------------------------------------------
export ROOTDIR="${PWD}"
export SCRIPTDIR="${ROOTDIR}/scripts"
export BUILDDIR="${ROOTDIR}/build"
export INSTALLDIR="${ROOTDIR}/install"
export SETUPFILE="${INSTALLDIR}/setup"
export SHA256_CHECKSUM="${SCRIPTDIR}/checksums.sha256"
export ARCH_FILE_TEMPLATE="${SCRIPTDIR}/arch_base.tmpl"

# ------------------------------------------------------------------------
# Load common variables and tools
# ------------------------------------------------------------------------
source "${SCRIPTDIR}"/common_vars.sh
source "${SCRIPTDIR}"/package_versions.sh
source "${SCRIPTDIR}"/tool_kit.sh

# ------------------------------------------------------------------------
# Documentation
# ------------------------------------------------------------------------
show_help() {
    cat<<EOF
This script will help you compile and install, or link libraries
CP2K depends on and setup a set of ARCH files that you can use
to compile CP2K.

USAGE:

$(basename $SCRIPT_NAME) [options]

OPTIONS:

-h, --help                Show this message.
-j <n>                    Number of processors to use for compilation, if
                          this option is not present, then the script
                          automatically tries to determine the number of
                          processors you have and try to use all of the
                          processors.
--no-check-certificate    If you encounter "certificate verification" errors
                          from wget or ones saying that "common name doesn't
                          match requested host name" while at tarball downloading
                          stage, then the recommended solution is to install
                          the newest wget release.  Alternatively, you can use
                          this option to bypass the verification and proceed with
                          the download. Security wise this should still be okay
                          as the installation script will check file checksums
                          after every tarball download. Nevertheless use this
                          option at your own risk.
--install-all             This option will set value of all --with-PKG
                          options to "install". You can selectively set
                          --with-PKG to another value again by having the
                          --with-PKG option placed AFTER this option on the
                          command line.
--mpi-mode                Selects which MPI flavour to use. Available options
                          are: mpich, openmpi and no.  By selecting no, you will
                          be disabling MPI support.  By default the script
                          will try to determine the flavour based on the MPI library
                          currently avaliable in your system path. For CRAY (CLE)
                          systems, the default flavour is mpich. Note that explicitly
                          setting --with-mpich or --with-openmpi options to values
                          other than no will also switch --mpi-mode to the respective
                          mode.
--math-mode               Selects which core math library to use. Available options
                          are: acml, cray, mkl, openblas and reflapack. cray
                          corresponds to cray libsci, and is the default for CRAY
                          (CLE) systems. For non-CRAY systems, if env variable MKLROOT
                          exists then mkl will be default, otherwise openblas is the
                          default option. Note that reflapack corresponds to the
                          reference LAPACK library, and is not really recommended for
                          CP2K binaries used for real simulations. Explicitly setting
                          --with-acml, --with-mkl or --with-openblas options will
                          switch --math-mode to the respective modes, BUT
                          --with-reflapack option do not have this effect.

The --enable-FEATURE options follow the rules:
  --enable-FEATURE=yes    Enable this particular feature
  --enable-FEATHRE=no     Disable this particular feature
  --enable-FEATURE        The option keyword alone is equivalent to
                          --enable-FEATURE=yes

  --enable-tsan           If you are installing GCC using this script
                          this option enables thread sanitizer support.
                          This is only relevant for debugging purposes.
                          Default = no
  --enable-gcc-master     If you are installing GCC using this script
                          this option forces the master development version
                          to be installed.
                          Default = no
  --enable-libxsmm-master If you are installing libxsmm using this script
                          this option forces the master development version
                          to be installed.
                          Default = no
  --enable-omp            Turn on OpenMP (threading) support.
                          Default = yes
  --enable-cuda           Turn on GPU (CUDA) support.
                          Default = no
  --enable-cray           Turn on or off support for CRAY Linux Environment
                          (CLE) manually. By default the script will automatically
                          detect if your system is CLE, and provide support
                          accordingly.

The --with-PKG options follow the rules:
  --with-PKG=install      Will download the package in \$PWD/build and
                          install the library package in \$PWD/install.
  --with-PKG=system       The script will then try to find the required
                          libraries of the package from the system path
                          variables such as PATH, LD_LIBRARY_PATH and
                          CPATH etc.
  --with-PKG=no           Do not use the package.
  --with-PKG=<path>       The package will be assumed to be installed in
                          the given <path>, and be linked accordingly.
  --with-PKG              The option keyword alone will be equivalent to
                          --with-PKG=install

  --with-gcc              The GCC compiler to use to compile CP2K
                          Default = system
  --with-binutils         GNU binutils
                          Default = system
  --with-make             GNU make
                          Default = system
  --with-cmake            Cmake utilities, required for building ParMETIS
                          Default = no
  --with-valgrind         Valgrind memory debugging tool, only used for
                          debugging purposes.
                          Default = no
  --with-lcov             LCOV code coverage utility, mainly used by CP2K developers
                          Default = no
  --with-openmpi          OpenMPI, important if you want parallel version
                          of CP2K.
                          Default = system
  --with-mpich            MPICH, MPI library like OpenMPI. one should
                          only use EITHER OpenMPI or MPICH and not both.
                          Default = system
  --with-libxc            libxc, exchange-correlation library. Needed for
                          QuickStep DFT and hybrid calculations.
                          Default = install
  --with-libint           libint, library for evaluation of two-body molecular
                          integrals, needed for hybrid functional calculations
                          Default = install
  --with-fftw             FFTW3, library for fast fourier transform
                          Default = install
  --with-reflapack        Reference (vanilla) LAPACK and BLAS linear algebra libraries.
                          One should use only ONE linear algebra library. This
                          one is really mostly used for debugging purposes as it is
                          non-optimised.
                          Default = no
  --with-acml             AMD core maths library, which provides LAPACK and BLAS
                          Default = system
  --with-mkl              Intel maths kernel library, which provides LAPACK and BLAS,
                          and depending on your system, may also provide ScaLAPACK.
                          If the MKL version of ScaLAPACK is found, then it will replace
                          the one specified by --with-scalapack option.
                          Default = system
  --with-openblas         OpenBLAS is a free high performance LAPACK and BLAS library,
                          the sucessor to GotoBLAS.
                          Default = install
  --with-scalapack        Parallel linear algebra library, needed for parallel
                          calculations.
                          Default = install
  --with-libsmm           CP2K's own small matrix multiplication library. An optimised
                          libsmm should increase the code performance. If you set
                          --with-libsmm=install, then instead of actually compiling
                          the library (which may take a long time), the script will
                          try to download a preexisting version from the CP2K website
                          that is compatable with your system.
                          Default = no
  --with-libxsmm          Small matrix multiplication library for x86_64 systems. If
                          your system arch is x86_64, then you can use libxsmm
                          instead of libsmm.
                          Default = install
  --with-elpa             Eigenvalue SoLvers for Petaflop-Applications library.
                          Fast library for large parallel jobs.
                          Default = no
  --with-ptscotch         PT-SCOTCH, only used if PEXSI is used
                          Default = no
  --with-parmetis         ParMETIS, and if --with-parmetis=install will also install
                          METIS, only used if PEXSI is used
                          Default = no
  --with-metis            METIS, --with-metis=install actuall does nothing, because
                          METIS is installed together with ParMETIS.  This option
                          is used to specify the METIS library if it is pre-installed
                          else-where. Only used if PEXSI is used
                          Default = no
  --with-superlu          SuperLU DIST, used only if PEXSI is used
                          Default = no
  --with-pexsi            Enable interface to PEXSI library
                          Default = no
  --with-quip             Enable interface to QUIP library
                          Default = no

FURTHER INSTRUCTIONS

All packages to be installed locally will be downloaded and build inside
./build, and then installed into package specific directories inside
./install.

Both ./build and ./install are safe to delete, as they contain
only the files and directories that are generated by this script. However,
once all the packages are installed, and you compile CP2K using the arch
files provided by this script, then you must keep ./install in exactly
the same location as it was first created, as it contains tools and libraries
your version of CP2K binary will depend on.

It should be safe to terminate running of this script in middle of a
build process.  The script will know if a package has been successfully
installed, and will just carry on and recompile and install the last
package it is working on. This is true even if you lose the content of
the entire ./build directory.

  +----------------------------------------------------------------+
  |  YOU SHOULD ALWAYS SOURCE ./install/setup BEFORE YOU RUN CP2K  |
  |  COMPILED WITH THIS TOOLCHAIN                                  |
  +----------------------------------------------------------------+

EOF
}

# ------------------------------------------------------------------------
# PACKAGE LIST: register all new dependent tools and libs here. Order
# is important, the first in the list gets installed first
# ------------------------------------------------------------------------
tool_list="binutils lcov valgrind make cmake gcc"
mpi_list="mpich openmpi"
math_list="mkl acml openblas reflapack"
lib_list="fftw libint libxc libsmm libxsmm scalapack elpa \
          ptscotch parmetis metis superlu pexsi quip"
package_list="$tool_list $mpi_list $math_list $lib_list"
# ------------------------------------------------------------------------

# first set everything to __DONTUSE__
for ii in $package_list ; do
    eval with_${ii}=__DONTUSE__
done

# ------------------------------------------------------------------------
# Work out default settings
# ------------------------------------------------------------------------

# tools to turn on by default:
with_binutils=__SYSTEM__
with_gcc=__SYSTEM__
with_make=__SYSTEM__

# libs to turn on by default, the math and mpi libraries are chosen by there respective modes:
with_fftw=__INSTALL__
with_libint=__INSTALL__
with_libxsmm=__INSTALL__
with_libxc=__INSTALL__
with_scalapack=__INSTALL__

# default math library settings, FAST_MATH_MODE picks the math library
# to use, and with_* defines the default method of installation if it
# is picked. For non-CRAY systems defaults to mkl if $MKLROOT is
# avaliable, otherwise defaults to openblas
if [ "$MKLROOT" ] ; then
    export FAST_MATH_MODE=mkl
else
    export FAST_MATH_MODE=openblas
fi
with_acml=__SYSTEM__
with_mkl=__SYSTEM__
with_openblas=__INSTALL__
with_reflapack=__INSTALL__

# for MPI, we try to detect system MPI variant
with_openmpi=__SYSTEM__
with_mpich=__SYSTEM__
if (command -v mpirun >&- 2>&-) ; then
    # check if we are dealing with openmpi or mpich
    if (mpirun --version 2>&1 | grep -s -q "HYDRA") ; then
        echo "MPI is detected and it appears to be MPICH"
        export MPI_MODE=mpich
    elif (mpirun --version 2>&1 | grep -s -q "Open MPI") ; then
        echo "MPI is detected and it appears to be OpenMPI"
        export MPI_MODE=openmpi
    else
        # default to mpich
        export MPI_MODE=mpich
    fi
else
    report_warning $LINENO "No MPI installation detected on you system. Ignore this message if you are using Cray Linux Environment"
    MPI_MODE=no
fi

# number of processors to use
export NPROCS=$(get_nprocs)

# default enable options
enable_tsan=__FALSE__
enable_gcc_master=__FALSE__
enable_libxsmm_master=__FALSE__
enable_omp=__TRUE__
if (command -v nvcc >&- 2>&-) ; then
   echo "nvcc found, enabling CUDA by default"
   enable_cuda=__TRUE__
else
   echo "nvcc not found, disabling CUDA by default"
   enable_cuda=__FALSE__
fi

# defaults for CRAY Linux Environment
if [ "$CRAY_LD_LIBRARY_PATH" ] ; then
    enable_cray=__TRUE__
    export FAST_MATH_MODE=cray
    # Default MPI used by CLE is assumed to be MPICH, in any case
    # don't use the installers for the MPI libraries
    with_mpich="__DONTUSE__"
    with_openmpi="__DONTUSE__"
    export MPI_MODE=mpich
    # set default value for some installers appropriate for CLE
    with_gcc="__DONTUSE__"
    with_fftw="__SYSTEM__"
    with_scalapack="__DONTUSE__"
else
    enable_cray=__FALSE__
fi

# ------------------------------------------------------------------------
# parse user options
# ------------------------------------------------------------------------
while [ $# -ge 1 ] ; do
    case $1 in
        -j)
            shift
            export NPROCS=$1
            ;;
        --no-check-certificate)
            export DOWNLOADER_FLAGS="-n"
            ;;
        --install-all)
            # set all package to __INSTALL__ status
            for ii in $package_list ; do
                eval with_${ii}=__INSTALL__
            done
            # default mpi-mode to MPICH
            MPI_MODE=mpich
            ;;
        --mpi-mode=*)
            user_input="${1#*=}"
            case "$user_input" in
                mpich)
                    export MPI_MODE=mpich
                    ;;
                openmpi)
                    export MPI_MODE=openmpi
                    ;;
                no)
                    export MPI_MODE=no
                    ;;
                *)
                    report_error ${LINENO} \
                                 "--mpi-mode currently only supports openmpi, mpich and no as options"
                    exit 1
                    ;;
            esac
            ;;
        --math-mode=*)
            user_input="${1#*=}"
            case "$user_input" in
                cray)
                    export FAST_MATH_MODE=cray
                    ;;
                mkl)
                    export FAST_MATH_MODE=mkl
                    ;;
                acml)
                    export FAST_MATH_MODE=acml
                    ;;
                openblas)
                    export FAST_MATH_MODE=openblas
                    ;;
                reflapack)
                    export FAST_MATH_MODE=reflapack
                    ;;
                *)
                    report_error ${LINENO} \
                    "--math-mode currently only supports mkl, acml, openblas and reflapack as options"
            esac
            ;;
        --enable-tsan*)
            enable_tsan=$(read_enable $1)
            if [ $enable_tsan = "__INVALID__" ] ; then
                report_error "invalid value for --enable-tsan, please use yes or no"
                exit 1
            fi
            ;;
        --enable-gcc-master*)
            enable_gcc_master=$(read_enable $1)
            if [ $enable_gcc_master = "__INVALID__" ] ; then
                report_error "invalid value for --enable-gcc-master, please use yes or no"
                exit 1
            fi
            ;;
        --enable-libxsmm-master*)
            enable_libxsmm_master=$(read_enable $1)
            if [ $enable_libxsmm_master = "__INVALID__" ] ; then
                report_error "invalid value for --enable-libxsmm-master, please use yes or no"
                exit 1
            fi
            ;;
        --enable-omp*)
            enable_omp=$(read_enable $1)
            if [ $enable_omp = "__INVALID__" ] ; then
                report_error "invalid value for --enable-omp, please use yes or no"
                exit 1
            fi
            ;;
        --enable-cuda*)
            enable_cuda=$(read_enable $1)
            if [ $enable_cuda = "__INVALID__" ] ; then
                report_error "invalid value for --enable-cuda, please use yes or no"
                exit 1
            fi
            ;;
        --enable-cray*)
            enable_cray=$(read_enable $1)
            if [ $enable_cray = "__INVALID__" ] ; then
                report_error "invalid value for --enable-cray, please use yes or no"
                exit 1
            fi
            ;;
        --with-gcc*)
            with_gcc=$(read_with $1)
            ;;
        --with-binutils*)
            with_binutils=$(read_with $1)
            ;;
        --with-make*)
            with_make=$(read_with $1)
            ;;
        --with-cmake*)
            with_cmake=$(read_with $1)
            ;;
        --with-lcov*)
            with_lcov=$(read_with $1)
            ;;
        --with-valgrind*)
            with_valgrind=$(read_with $1)
            ;;
        --with-mpich*)
            with_mpich=$(read_with $1)
            if [ $with_mpich != __DONTUSE__ ] ; then
                export MPI_MODE=mpich
            fi
            ;;
        --with-openmpi*)
            with_openmpi=$(read_with $1)
            if [ $with_openmpi != __DONTUSE__ ] ; then
                export MPI_MODE=openmpi
            fi
            ;;
        --with-libint*)
            with_libint=$(read_with $1)
            ;;
        --with-libxc*)
            with_libxc=$(read_with $1)
            ;;
        --with-fftw*)
            with_fftw=$(read_with $1)
            ;;
        --with-reflapack*)
            with_reflapack=$(read_with $1)
            ;;
        --with-mkl*)
            with_mkl=$(read_with $1)
            if [ $with_mkl != __DONTUSE__ ] ; then
                export FAST_MATH_MODE=mkl
            fi
            ;;
        --with-acml*)
            with_acml=$(read_with $1)
            if [ $with_acml != __DONTUSE__ ] ; then
                export FAST_MATH_MODE=acml
            fi
            ;;
        --with-openblas*)
            with_openblas=$(read_with $1)
            if [ $with_openblas != __DONTUSE__ ] ; then
                export FAST_MATH_MODE=openblas
            fi
            ;;
        --with-scalapack*)
            with_scalapack=$(read_with $1)
            ;;
        --with-libsmm*)
            with_libsmm=$(read_with $1)
            ;;
        --with-libxsmm*)
            with_libxsmm=$(read_with $1)
            ;;
        --with-elpa*)
            with_elpa=$(read_with $1)
            ;;
        --with-ptscotch*)
            with_ptscotch=$(read_with $1)
            ;;
        --with-parmetis*)
            with_parmetis=$(read_with $1)
            ;;
        --with-metis*)
            with_metis=$(read_with $1)
            ;;
        --with-superlu*)
            with_superlu=$(read_with $1)
            ;;
        --with-pexsi*)
            with_pexsi=$(read_with $1)
            ;;
        --with-quip*)
            with_quip=$(read_with $1)
            ;;
        *)
            show_help
            exit 0
            ;;
    esac
    shift
done

# consolidate settings after user input
export ENABLE_TSAN=$enable_tsan
export ENABLE_OMP=$enable_omp
export ENABLE_CUDA=$enable_cuda
export ENABLE_CRAY=$enable_cray
[ "$enable_gcc_master" = "__TRUE__" ] && export gcc_ver=master
[ "$enable_libxsmm_master" = "__TRUE__" ] && export libxsmm_ver=master
[ "$with_valgrind" != "__DONTUSE__"  ] && export ENABLE_VALGRIND="__TRUE__"
[ "$with_lcov" != "__DONTUSE__" ] && export ENABLE_COVERAGE="__TRUE__"

# ------------------------------------------------------------------------
# Check and solve known conflicts before installations proceed
# ------------------------------------------------------------------------

# GCC thread sanitizer conflicts
if [ $ENABLE_TSAN = "__TRUE__" ] ; then
    if [ "$with_openblas" != "__DONTUSE__" ] ; then
        echo "TSAN is enabled, canoot use openblas, we will use reflapack instead"
        [ "$with_reflapack" = "__DONTUSE__" ] && with_reflapack="__INSTALL__"
        export FAST_MATH_MODE=reflapack
    fi
    echo "TSAN is enabled, canoot use libsmm"
    with_libsmm="__DONTUSE__"
fi

# valgrind conflicts
if [ "$ENABLE_VALGRIND" = "__TRUE__" ] ; then
    if [ "$with_reflapack" = "__DONTUSE__" ] ; then
        echo "reflapack is automatically installed when valgrind is enabled"
        with_reflapack="__INSTALL__"
    fi
fi

# mpi library conflicts
if [ $MPI_MODE = no ] ; then
    if [ "$with_scalapack" != "__DONTUSE__"  ] ; then
        echo "Not using MPI, so scalapack is disabled."
        with_scalapack="__DONTUSE__"
    fi
    if [ "$with_elpa" != "__DONTUSE__" ] ; then
        echo "Not using MPI, so ELPA is disabled."
        with_elpa="__DONTUSE__"
    fi
    if [ "$with_pexi" != "__DONTUSE__" ] ; then
        echo "Not using MPI, so PEXSI is disabled."
        with_pexsi="__DONTUSE__"
    fi
else
    # if gcc is installed, then mpi needs to be installed too
    if [ "$with_gcc" = "__INSTALL__" ] ; then
        echo "You have chosen to install GCC, therefore MPI libraries will have to be installed too"
        with_openmpi="__INSTALL__"
        with_mpich="__INSTALL__"
    fi
fi

# PESXI and its dependencies
if [ "$with_pexsi" = "__DONTUSE__" ] ; then
    if [ "$with_ptscotch" != "__DONTUSE__" ] ; then
        echo "Not using PEXSI, so PT-Scotch is disabled."
        with_ptscotch="__DONTUSE__"
    fi
    if [ "$with_parmetis" != "__DONTUSE__" ] ; then
        echo "Not using PEXSI, so ParMETIS is disabled."
        with_parmetis="__DONTUSE__"
    fi
    if [ "$with_metis" != "__DONTUSE__" ] ; then
        echo "Not using PEXSI, so METIS is disabled."
        with_metis="__DONTUSE__"
    fi
    if [ "$with_superlu" != "__DONTUSE__" ] ; then
        echo "Not using PEXSI, so SuperLU-DIST is disabled."
        with_superlu="__DONTUSE__"
    fi
elif [ "$with_pexsi" = "__INSTALL__" ] ; then
    [ "$with_ptscotch" = "__DONTUSE__" ] && with_ptscotch="__INSTALL__"
    [ "$with_parmetis" = "__DONTUSE__" ] && with_parmetis="__INSTALL__"
    [ "$with_superlu" = "__DONTUSE__" ] && with_superlu="__INSTALL__"
else
    if [ "$with_ptscotch" = "__DONTUSE__" ] ; then
        report_error "For PEXSI to work you need a working PT-Scotch library use --with-ptscotch option to specify if you wish to install the library or specify its location."
        exit 1
    fi
    if [ "$with_parmetis" = "__DONTUSE__" ] ; then
        report_error "For PEXSI to work you need a working PARMETIS library use --with-parmetis option to specify if you wish to install the library or specify its location."
        exit 1
    fi
    if [ "$with_metis" = "__DONTUSE__" ] ; then
        report_error "For PEXSI to work you need a working METIS library use --with-metis option to specify if you wish to install the library or specify its location."
        exit 1
    fi
    if [ "$with_superlu" = "__DONTUSE__" ] ; then
        report_error "For PEXSI to work you need a working SuperLU-DIST library use --with-superlu option to specify if you wish to install the library or specify its location."
        exit 1
    fi
fi
# ParMETIS requires cmake, it also installs METIS if it is chosen
# __INSTALL__ option
if [ "$with_parmetis" = "__INSTALL__" ] ; then
    [ "$with_cmake" = "__DONTUSE__" ] && with_cmake="__INSTALL__"
    with_metis="__INSTALL__"
fi

# ------------------------------------------------------------------------
# Preliminaries
# ------------------------------------------------------------------------

mkdir -p "$INSTALLDIR"

# variables used for generating cp2k ARCH file
export CP_DFLAGS=''
export CP_LIBS=''
export CP_CFLAGS='IF_OMP(-fopenmp|)'
export CP_LDFLAGS="-Wl,--enable-new-dtags"

# ------------------------------------------------------------------------
# Start writing setup file
# ------------------------------------------------------------------------
cat <<EOF > "$SETUPFILE"
#!/bin/bash
source "${SCRIPTDIR}/tool_kit.sh"
EOF

# ------------------------------------------------------------------------
# Special settings for CRAY Linux Environment (CLE)
# ------------------------------------------------------------------------
if [ "$ENABLE_CRAY" = "__TRUE__" ] ; then
    echo "------------------------------------------------------------------------"
    echo "CRAY Linux Environment (CLE) is detected"
    echo "------------------------------------------------------------------------"
    # add cray paths to system search path
    export LIB_PATHS="CRAY_LD_LIBRARY_PATH ${LIB_PATHS}"
    # set compilers to CLE wrappers
    check_command cc
    check_command ftn
    check_command CC
    export CC=cc
    export FC=ftn
    export F77=ftn
    export F90=ftn
    export CXX=CC
    export MPICC=cc
    export MPIFC=ftn
    export MPIF77=ftn
    export MPIF90=ftn
    export MPICXX=CC
    # CRAY libsci should contains core math libraries, scalapack
    # doesn't need LDFLAGS or CFLAGS, nor do the one need to
    # explicitly link the math and scalapack libraries, as all is
    # taken care of by the cray compiler wrappers.
    if [ "$with_scalapack" = "__DONTUSE__" ] ; then
        export CP_DFLAGS="${CP_DFLAGS} IF_MPI(-D__SCALAPACK|)"
    fi
    case $MPI_MODE in
        mpich)
            if [ "$MPICH_DIR" ] ; then
                cray_mpich_include_path="$MPICH_DIR/include"
                cray_mpich_lib_path="$MPICH_DIR/lib"
                export INCLUDE_PATHS="$INCLUDE_PATHS cray_mpich_include_path"
                export LIB_PATHS="$LIB_PATHS cray_mpich_lib_path"
            fi
            if [ "$with_mpich" = "__DONTUSE__" ] ; then
                add_include_from_paths MPI_CFLAGS "mpi.h" $INCLUDE_PATHS
                add_include_from_paths MPI_LDFLAGS "libmpi.*" $LIB_PATHS
                export MPI_CFLAGS
                export MPI_LDFLAGS
                export MPI_LIBS=" "
                export CP_DFLAGS="${CP_DFLAGS} IF_MPI(-D__parallel -D__MPI_VERSION=3|)"
            fi
            ;;
        openmpi)
            if [ "$with_openmpi" = "__DONTUSE__" ] ; then
                add_include_from_paths MPI_CFLAGS "mpi.h" $INCLUDE_PATHS
                add_include_from_paths MPI_LDFLAGS "libmpi.*" $LIB_PATHS
                export MPI_CFLAGS
                export MPI_LDFLAGS
                export MPI_LIBS="-lmpi -lmpi_cxx"
                export CP_DFLAGS="${CP_DFLAGS} IF_MPI(-D__parallel -D__MPI_VERSION=3|)"
            fi
            ;;
    esac
    check_lib -lz
    check_lib -ldl
    export CRAY_EXTRA_LIBS="-lz -ldl"
    # the space is intentional, so that the variable is non-empty and
    # can pass require_env checks
    export SCALAPACK_LDFLAGS=" "
    export SCALAPACK_LIBS=" "
fi

# ------------------------------------------------------------------------
# Installing tools required for building CP2K and associated libraries
# ------------------------------------------------------------------------

echo "Compiling with $NPROCS processes."

# set environment for compiling compilers and tools required for CP2K
# and libraries it depends on
export CFLAGS=${CFLAGS:-"-O2 -g -Wno-error"}
export FFLAGS=${FFLAGS:-"-O2 -g -Wno-error"}
export FCLAGS=${FCLAGS:-"-O2 -g -Wno-error"}
export F90FLAGS=${F90FLAGS:-"-O2 -g -Wno-error"}
export F77FLAGS=${F77FLAGS:-"-O2 -g -Wno-error"}
export CXXFLAGS=${CXXFLAGS:-"-O2 -g -Wno-error"}

# need to setup tools after all of the tools are built. We should use
# consistent pairs of gcc and binutils etc for make. So we use system
# tool sets to compile the tool sets used to compile CP2K
for ii in $tool_list ; do
    install_mode="$(eval echo \${with_${ii}})"
    "${SCRIPTDIR}"/install_${ii}.sh "$install_mode"
done
for ii in $tool_list ; do
    load "${BUILDDIR}/setup_${ii}"
done

# ------------------------------------------------------------------------
# Install or compile packages using newly installed tools
# ------------------------------------------------------------------------

# setup compiler flags, leading to nice stack traces on crashes but
# still optimised
CFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"
FFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"
F77FLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"
F90FLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"
FCFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"
CXXFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"

export CFLAGS=$(allowed_gcc_flags $CFLAGS)
export FFLAGS=$(allowed_gfortran_flags $FFLAGS)
export F77FLAGS=$(allowed_gfortran_flags $F77FLAGS)
export F90FLAGS=$(allowed_gfortran_flags $F90FLAGS)
export FCFLAGS=$(allowed_gfortran_flags $FCFLAGS)
export CXXFLAGS=$(allowed_gxx_flags $CXXFLAGS)

export LDFLAGS="$TSANFLAGS"

# get system arch information using OpenBLAS prebuild
"${SCRIPTDIR}"/get_openblas_arch.sh; load "${BUILDDIR}/openblas_arch"

# MPI libraries
case "$MPI_MODE" in
    mpich)
        "${SCRIPTDIR}"/install_mpich.sh "${with_mpich}"; load "${BUILDDIR}/setup_mpich"
        ;;
    openmpi)
        "${SCRIPTDIR}"/install_openmpi.sh "${with_openmpi}"; load "${BUILDDIR}/setup_openmpi"
        ;;
esac

# math core libraries, need to use reflapack for valgrind builds, as
# many fast libraries are not necesarily valgrind clean
export REF_MATH_CFLAGS=''
export REF_MATH_LDFLAGS=''
export REF_MATH_LIBS=''
export FAST_MATH_CFLAGS=''
export FAST_MATH_LDFLAGS=''
export FAST_MATH_LIBS=''

"${SCRIPTDIR}"/install_reflapack.sh "${with_reflapack}"; load "${BUILDDIR}/setup_reflapack"
case "$FAST_MATH_MODE" in
    mkl)
        "${SCRIPTDIR}"/install_mkl.sh "${with_mkl}"; load "${BUILDDIR}/setup_mkl"
        ;;
    acml)
        "${SCRIPTDIR}"/install_acml.sh "${with_acml}"; load "${BUILDDIR}/setup_acml"
        ;;
    openblas)
        "${SCRIPTDIR}"/install_openblas.sh "${with_openblas}"; load "${BUILDDIR}/setup_openblas"
        ;;
    cray)
        # note the space is intentional so that the variable is
        # non-empty and can pass require_env checks
        export FAST_MATH_LDFLAGS="${FAST_MATH_LDFLAGS} "
        export FAST_MATH_LIBS="${FAST_MATH_LIBS} ${CRAY_EXTRA_LIBS}"
        ;;
esac

if [ $ENABLE_VALGRIND = "__TRUE__" ] ; then
    export MATH_CFLAGS="${REF_MATH_CFLAGS}"
    export MATH_LDFLAGS="${REF_MATH_LDFLAGS}"
    export MATH_LIBS="${REF_MATH_LIBS}"
else
    export MATH_CFLAGS="${FAST_MATH_CFLAGS}"
    export MATH_LDFLAGS="${FAST_MATH_LDFLAGS}"
    export MATH_LIBS="${FAST_MATH_LIBS}"
fi

export CP_CFLAGS="${CP_CFLAGS} IF_DEBUG(${REF_MATH_CFLAGS}|IF_VALGRIND(${REF_MATH_CFLAGS}|${FAST_MATH_CFLAGS}))"
export CP_LDFLAGS="${CP_LDFLAGS} IF_DEBUG(${REF_MATH_LDFLAGS}|IF_VALGRIND(${REF_MATH_LDFLAGS}|${FAST_MATH_LDFLAGS}))"
export CP_LIBS="${CP_LIBS} IF_DEBUG(${REF_MATH_LIBS}|IF_VALGRIND(${REF_MATH_LIBS}|${FAST_MATH_LIBS}))"

# other libraries
for ii in $lib_list ; do
    install_mode="$(eval echo \${with_${ii}})"
    "${SCRIPTDIR}"/install_${ii}.sh "$install_mode"
    load "${BUILDDIR}/setup_${ii}"
done

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
BASEFLAGS="-march=native -fno-omit-frame-pointer -g ${TSANFLAGS}"
OPT_FLAGS="-O3 -funroll-loops -ffast-math"
NOOPT_FLAGS="-O1"

# those flags that do not influence code generation are used always, the others if debug
FCDEB_FLAGS="-ffree-form -std=f2003 -fimplicit-none"
FCDEB_FLAGS_DEBUG="-fsanitize=leak -fcheck=bounds,do,recursion,pointer -ffpe-trap=invalid,zero,overflow -finit-real=snan -fno-fast-math"

# code coverage generation flags
COVERAGE_FLAGS="-O1 -coverage -fkeep-static-functions"
COVERAGE_DFLAGS="-D__NO_ABORT"

# profile based optimization, see https://www.cp2k.org/howto:pgo
PROFOPT_FLAGS="\$(PROFOPT)"

# special flags for gfortran
# https://gcc.gnu.org/onlinedocs/gfortran/Error-and-Warning-Options.html
# we error out for these warnings (-Werror=uninitialized -Wno-maybe-uninitialized -> error on variables that must be used uninitialized)
WFLAGS_ERROR="-Werror=aliasing -Werror=ampersand -Werror=c-binding-type -Werror=intrinsic-shadow -Werror=intrinsics-std -Werror=line-truncation -Werror=tabs -Werror=realloc-lhs-all -Werror=target-lifetime -Werror=underflow -Werror=unused-but-set-variable -Werror=unused-variable -Werror=unused-dummy-argument -Werror=conversion -Werror=zerotrip -Werror=uninitialized -Wno-maybe-uninitialized"
# we just warn for those (that eventually might be promoted to WFLAGSERROR). It is useless to put something here with 100s of warnings.
WFLAGS_WARN="-Wuse-without-only"
# while here we collect all other warnings, some we'll ignore
WFLAGS_WARNALL="-pedantic -Wall -Wextra -Wsurprising -Wunused-parameter -Warray-temporaries -Wcharacter-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wreal-q-constant -Wunused-parameter -Walign-commons -Wfunction-elimination -Wrealloc-lhs -Wcompare-reals -Wzerotrip"

# IEEE_EXCEPTIONS dependency
IEEE_EXCEPTIONS_DFLAGS="-D__HAS_IEEE_EXCEPTIONS"

# check all of the above flags, filter out incompatable flags for the
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

# check if ieee_exeptions module is avaliable for the current version
# of gfortran being used
if ! (check_gfortran_module ieee_exceptions) ; then
    IEEE_EXCEPTIONS_DFLAGS=""
fi

# contagnate the above flags into WFLAGS, FCDEBFLAGS, DFLAGS and
# finally into FCFLAGS and CFLAGS
WFLAGS="$WFLAGS_ERROR $WFLAGS_WARN IF_WARNALL(${WFLAGS_WARNALL}|)"
FCDEBFLAGS="$FCDEB_FLAGS IF_DEBUG($FCDEB_FLAGS_DEBUG|)"
DFLAGS="${CP_DFLAGS} IF_DEBUG($IEEE_EXCEPTIONS_DFLAGS|) IF_COVERAGE($COVERAGE_DFLAGS|)"
# language independent flags
# valgrind with avx can lead to spurious out-of-bound results
G_CFLAGS="$BASEFLAGS IF_VALGRIND(-mno-avx -mno-avx2|)"
G_CFLAGS="$G_CFLAGS IF_COVERAGE($COVERAGE_FLAGS|IF_DEBUG($NOOPT_FLAGS|$OPT_FLAGS))"
G_CFLAGS="$G_CFLAGS IF_DEBUG(|$PROFOPT_FLAGS)"
G_CFLAGS="$G_CFLAGS $CP_CFLAGS"
# FCFLAGS, for gfortran
FCFLAGS="$G_CFLAGS \$(FCDEBFLAGS) \$(WFLAGS) \$(DFLAGS)"
# CFLAGS, spcial flags for gcc (currently none)
CFLAGS="$G_CFLAGS \$(DFLAGS)"

# Linker flags
LDFLAGS="\$(FCFLAGS) ${CP_LDFLAGS}"

# Library flags
# add standard libs
LIBS="${CP_LIBS} -lstdc++"

# CUDA stuff
CUDA_LIBS="-lcudart -lcufft -lcublas -lrt IF_DEBUG(-lnvToolsExt|)"
CUDA_DFLAGS="-D__ACC -D__DBCSR_ACC -D__PW_CUDA IF_DEBUG(-D__CUDA_PROFILING|)"
if [ "$ENABLE_CUDA" = __TRUE__ ] ; then
    LIBS="${LIBS} IF_CUDA(${CUDA_LIBS}|)"
    DFLAGS="IF_CUDA(${CUDA_DFLAGS}|) ${DFLAGS}"
    NVFLAGS="-arch sm_35 \$(DFLAGS)"
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
    local __full_flag_list="MPI OMP DEBUG CUDA WARNALL VALGRIND COVERAGE"
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
NVCC        = \${NVCC} -D__GNUC__=4 -D__GNUC_MINOR__=9 -Xcompiler=--std=gnu++98
NVFLAGS     = \${NVFLAGS}
EOF
    fi
    if [ "$__WARNALL" = "on" ] ; then
        cat <<EOF >> $__filename
#
FCLOGPIPE   =  2> \\\$(notdir \\\$<).warn
export LC_ALL=C
EOF
    fi
    # replace variable values in output file using eval
    local __TMPL=$(cat $__filename)
    eval "printf \"${__TMPL}\n\"" > $__filename
    # pass this to parsers to replace all of the IF_XYZ statements
    python ${SCRIPTDIR}/parse_if.py $__filename $__flags
    echo "Wrote ${INSTALLDIR}/arch/$__filename"
}

rm -f ${INSTALLDIR}/arch/local*
# normal production arch files
    { gen_arch_file "local.sopt" ;          arch_vers="sopt"; }
    { gen_arch_file "local.sdbg" DEBUG;     arch_vers="${arch_vers} sdbg"; }
[ "$ENABLE_OMP" = __TRUE__ ] && \
    { gen_arch_file "local.ssmp" OMP;       arch_vers="${arch_vers} ssmp"; }
[ "$MPI_MODE" != no ] && \
    { gen_arch_file "local.popt" MPI;       arch_vers="${arch_vers} popt"; }
[ "$MPI_MODE" != no ] && \
    { gen_arch_file "local.pdbg" MPI DEBUG; arch_vers="${arch_vers} pdbg"; }
[ "$MPI_MODE" != no ] && \
[ "$ENABLE_OMP" = __TRUE__ ] && \
    { gen_arch_file "local.psmp" MPI OMP;   arch_vers="${arch_vers} psmp"; }
[ "$MPI_MODE" != no ] && \
[ "$ENABLE_OMP" = __TRUE__ ] && \
    gen_arch_file "local_warn.psmp" MPI OMP WARNALL
# cuda enabled arch files
if [ "$ENABLE_CUDA" = __TRUE__ ] ; then
    [ "$ENABLE_OMP" = __TRUE__ ] && \
      gen_arch_file "local_cuda.ssmp"          CUDA OMP
    [ "$MPI_MODE" != no ] && \
    [ "$ENABLE_OMP" = __TRUE__ ] && \
      gen_arch_file "local_cuda.psmp"          CUDA OMP MPI
    [ "$ENABLE_OMP" = __TRUE__ ] && \
      gen_arch_file "local_cuda.sdbg"          CUDA DEBUG OMP
    [ "$MPI_MODE" != no ] && \
    [ "$ENABLE_OMP" = __TRUE__ ] && \
      gen_arch_file "local_cuda.pdbg"          CUDA DEBUG OMP MPI
    [ "$MPI_MODE" != no ] && \
    [ "$ENABLE_OMP" = __TRUE__ ] && \
      gen_arch_file "local_cuda_warn.psmp"     CUDA MPI OMP WARNALL
fi
# valgrind enabled arch files
if [ "$ENABLE_VALGRIND" = __TRUE__ ] ; then
      gen_arch_file "local_valgrind.sdbg"      VALGRIND
    [ "$MPI_MODE" != no ] && \
      gen_arch_file "local_valgrind.pdbg"      VALGRIND MPI
fi
# coverage enabled arch files
if [ "$ENABLE_COVERAGE" = __TRUE__ ]; then
      gen_arch_file "local_coverage.sdbg"      COVERAGE
    [ "$MPI_MODE" != no ] && \
      gen_arch_file "local_coverage.pdbg"      COVERAGE MPI
    [ "$ENABLE_CUDA" = __TRUE__ ] && \
      gen_arch_file "local_coverage_cuda.pdbg" COVERAGE MPI CUDA
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
  cd cp2k/makefiles/
  make -j ${NPROCS} ARCH=local VERSION="${arch_vers}"

arch files for GPU enabled CUDA versions are named "local_cuda.*"
arch files for valgrind versions are named "local_valgrind.*"
arch files for coverage versions are named "local_coverage.*"
EOF

#EOF
