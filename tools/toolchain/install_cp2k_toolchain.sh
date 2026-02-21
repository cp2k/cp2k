#!/bin/bash -e

# Disabled shellcheck items: SC1091 for external scripts, SC2034 for unused
# variables, SC2124 for concatenating toolchain options with $@
# shellcheck disable=SC1091,SC2034,SC2124

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

# +---------------------------------------------------------------------------+
# |   CP2K: A general program to perform molecular dynamics simulations       |
# |   Copyright 2000-2022 CP2K developers group <https://cp2k.org>            |
# |                                                                           |
# |   SPDX-License-Identifier: GPL-2.0-or-later                               |
# +---------------------------------------------------------------------------+
#
#
# *****************************************************************************
#> \brief    This script will compile and install or link existing tools and
#>           libraries CP2K depends on and generate a set of ARCH files which
#>           can be used to compile CP2K
#> \history  Created on Friday, 2016/02/05
#            Update for Intel (17.01.2022, MK)
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

# ------------------------------------------------------------------------
# Make a copy of all options for $SETUPFILE
# ------------------------------------------------------------------------
TOOLCHAIN_OPTIONS="$@"

# ------------------------------------------------------------------------
# Exit this script if it is not called from the ./tools/toolchain directory
# ------------------------------------------------------------------------
if [ "${ROOTDIR}" != "${SCRIPT_DIR}" ]; then
  cat << EOF
ERROR: (${SCRIPT_NAME}, line ${LINENO}) Incorrect execution location.
The absolute path of the main toolchain script is at:
  ${SCRIPT_DIR}
Actual working directory where it is currently called:
  ${ROOTDIR}
Please enter the absolute path above before executing the main toolchain script
so that subsequent scripts can be found and files can be placed correctly.
EOF
  exit 1
fi

# ------------------------------------------------------------------------
# Load common variables and tools
# ------------------------------------------------------------------------
source "${SCRIPTDIR}"/common_vars.sh
source "${SCRIPTDIR}"/tool_kit.sh

# ------------------------------------------------------------------------
# Documentation
# ------------------------------------------------------------------------
show_help() {
  cat << EOF
This script will help you prepare the toolchain for compiling and using CP2K.
There are a number of dependency packages for CP2K, which may be downloaded
from internet and freshly installed, or detected from system path and linked.
Once these dependencies are ready, the CMake options for compiling CP2K and
further instructions will be printed at the end of this script's execution.
See README.md under the toolchain directory for more information.

USAGE:

$(basename "$SCRIPT_NAME") [options]

OPTIONS:

-h, --help                Show this message and exit.
-j <n>                    Number of processors for parallel compiling.
                          If omitted, the script will automatically try to
                          determine the number of available processors and use
                          all of them by default. Set "-j 1" explicitly if this
                          is not desired.
--no-check-certificate    Bypass verification of server's certificate while
                          downloading anything from internet via wget command.
                          In case wget errors about "certificate verification"
                          or "common name doesn't match requested host name"
                          occur while downloading tarballs, the recommended
                          solution is to install the latest wget release;
                          alternatively, this script can be rerun with this
                          option in command line.
                          Security wise this should still be okay as sha256
                          checksums are checked after every tarball download.
                          Nevertheless, use this option at your own risk.
--install-all             Set value of all --with-PKG options to "install"
                          (except --with-intel, --with-intelmpi, --with-amd).
                          By default, GNU compiler and MPICH are installed.
                          AFTER this option on the command line, you can set
                          specific --with-PKG options to another value again
                          for selective fine controls (see below).
--mpi-mode                Select which MPI flavour to use. Available options
                          are: mpich, openmpi, intelmpi, and no.
                          If omitted, the script will automatically try to
                          determine the flavour based on the MPI library
                          available in system path; alternatively, if any of
                          --with-mpich, --with-openmpi or --with-intelmpi
                          options (see below) is explicitly set to values other
                          than no, the script will follow.
                          For CRAY (CLE) systems, the default flavour is mpich.
                          By selecting "no", MPI is not supported and disabled.
--math-mode               Select which core math library to use. Available
                          options are: acml, cray, mkl, and openblas.
                          If omitted, for non-CRAY systems, mkl is the default
                          when environment variable \$MKLROOT exists, otherwise
                          openblas is the drefault; for CRAY (CLE) systems,
                          cray is the default, corresponding to cray libsci.
                          Alternatively, if any of --with-acml, --with-mkl,
                          or --with-openblas options (see below) is explicitly
                          set to values other than no, the script will follow.
--target-cpu              Select the target CPU architecture for compiling.
                          This option determines the value of -mtune flag in
                          CFLAGS for compilers. If omitted or set to "native",
                          compiling will be tuned/optimized for the native host
                          system, and some instruction sets will be detected
                          including AVX, AVX2, AVX512, SSE4, etc. Alternatively
                          this option can be set depending on actual target
                          CPU microarchitecture, e.g. "haswell", "skylake", or
                          just "generic".
                          Default = native
--gpu-ver                 Select the target GPU architecture for compiling.
                          Available options are: K20X, K40, K80, P100, V100,
                          A100, H100, A40, Mi50, Mi100, Mi250, and no.
                          This option determines the value of nvcc -arch flag.
                          Default = no
--libint-lmax             Maximum supported angular momentum by libint if the
                          --with-libint option is set to "install" (see below).
                          Available options are: 4, 5, 6, 7. Higher values will
                          increase build time and library size.
                          Default = 5
--log-lines               In case a package build returns a non-zero exit code,
                          the number of final lines from the log file in the
                          corresponding directory that will be dumped to the
                          command line. Due to limited length, this snippet may
                          not contain the very first warning or error message.
                          Default = 200
--dry-run                 After writing toolchain config and env files to the
                          install directory, just print a list of effective
                          settings after resolving known conflicts, then exit.
                          Like --help, do not actually download tarballs or
                          build packages.

The --enable-FEATURE options follow the rules:
  --enable-FEATURE=yes    Enable this particular feature.
  --enable-FEATURE=no     Disable this particular feature.
  --enable-FEATURE        The option keyword alone is equivalent to
                          --enable-FEATURE=yes
Specific options:
  --enable-tsan           Turn on Thread Sanitizer (TSAN) for GNU compiler.
                          Default = no
  --enable-cuda           Turn on GPU (CUDA) support.
                          Can be combined with --enable-opencl.
                          Default = no
  --enable-hip            Turn on GPU (HIP) support.
                          Default = no
  --enable-opencl         Turn on OpenCL (GPU) support. Requires the OpenCL
                          development packages and runtime. If combined with
                          --enable-cuda, OpenCL alongside of CUDA is used.
                          Default = no
  --enable-cray           Turn on CRAY Linux Environment (CLE) support.
                          If omitted, the script will automatically try to
                          determine if your system is CLE by detecting
                          environment variable \$CRAY_LD_LIBRARY_PATH and
                          provide support accordingly.

The --with-PKG options follow the rules:
  --with-PKG=no           Do not use the package.
  --with-PKG=install      Download and decompress the package in \$PWD/build,
                          and install the library package in \$PWD/install.
  --with-PKG=system       Find the required libraries of the package from the
                          system path variables such as PATH, LD_LIBRARY_PATH
                          and CPATH etc.
  --with-PKG=<path>       The package will be assumed to be installed in
                          the given <path>, and be linked accordingly.
  --with-PKG              The option keyword alone will be equivalent to
                          --with-PKG=install
Specific options:
  --with-gcc              Use the GNU compilers (gcc, g++, gfortran) to build
                          newly installed dependencies and CP2K.
                          Only one of --with-gcc, --with-intel and --with-amd
                          should be used.
                          Default = system
  --with-intel            Use the Intel compilers (icx, icpx, ifort) to build
                          newly installed dependencies and CP2K.
                          Default = system
  --with-ifx              If yes, use the new Intel Fortran compiler ifx
                          instead of ifort to compile CP2K.
                          Default = no
  --with-amd              Use the AMD compilers (clang, clang++, flang) to
                          build newly installed dependencies and CP2K.
                          Default = system
  --with-cmake            CMake utilities.
                          Default = install
  --with-ninja            Ninja utilities.
                          Default = install
  --with-openmpi          Use OpenMPI library for parallel versions of
                          newly installed dependencies and CP2K.
                          Only one of --with-openmpi, --with-mpich and
                          --with-intelmpi should be used, and --mpi-mode option
                          (see above) should be consistent.
                          Default = system
  --with-mpich            Use MPICH library for parallel versions of
                          newly installed dependencies and CP2K.
                          Default = system
  --with-mpich-device     Select the MPICH device, implying use of MPICH.
                          Default = ch4
  --with-intelmpi         Use Intel MPI library for parallel versions of
                          newly installed dependencies and CP2K.
                          Default = system
  --with-openblas         Use OpenBLAS, which provides LAPACK and BLAS library,
                          and is also used to get arch information.
                          --math-mode option (see above) should be consistent.
                          Default = install
  --with-mkl              Use Intel Math Kernel Library (MKL), which provides
                          LAPACK and BLAS library.
                          If MKL's FFTW3 interface is suitable (no FFTW-MPI
                          support), it replaces the FFTW library.
                          If the ScaLAPACK component is found, it replaces the
                          one specified by --with-scalapack.
                          Default = system
  --with-acml             Use AMD core maths library (ACML), which provides
                          LAPACK and BLAS library.
                          Default = system
  --with-gmp              Enable GMP library, optional dependency of GreenX.
                          Default = no
  --with-fftw             Enable FFTW3 library for fast fourier transform.
                          Default = install
  --with-libxc            Enable libxc for exchange-correlation in QuickStep
                          DFT (pure and hybrid functionals) calculations.
                          Default = install
  --with-libint           Enable libint for two-body molecular integrals in
                          Hartree-Fock and hybrid functional calculations.
                          Default = install
  --with-greenx           Enable GreenX library for Minimax grids and Padé
                          analytic continuation in RT-BSE.
                          This package requires CMake, BLAS and LAPACK.
                          Default = no
  --with-cosma            Enable COSMA as a replacement for ScaLAPACK in matrix
                          multiplication. If set to "install", COSTA and TileMM
                          will also be installed; if CUDA and/or HIP support is
                          enabled too, respective versions will all be built.
                          Default = install
  --with-libxsmm          Enable libxsmm as a small matrix multiplication
                          library. Installing is only supported on arch
                          x86_64 or arm64.
                          Default = install
  --with-scalapack        Enable ScaLAPACK for parallel linear algebra
                          calculations.
                          Default = install
  --with-elpa             Enable ELPA library as eigenvalue solver for large
                          parallel jobs.
                          Default = install
  --with-ace              Enable interface to ML-pace library.
                          Default = no
  --with-deepmd           Enable interface to DeePMD-kit library.
                          This does not include other DeepModeling utilities
                          like DP-GEN or dpdata.
                          Default = no
  --with-gsl              Enable the GNU scientific library (GSL).
                          This package is required for PLUMED and SIRIUS.
                          Default = install
  --with-libtorch         Enable libtorch as a machine learning framework.
                          This package is required for NequIP and Allegro, and
                          also for installing DeePMD-kit.
                          Default = no
  --with-plumed           Enable interface to the PLUMED library for enhanced
                          sampling methods.
                          This package requires MPI, GSL and FFTW.
                          Default = no
  --with-hdf5             Enable the hdf5 library for file format support.
                          This package is used by sirius and trexio.
                          Default = install
  --with-libsmeagol       Enable interface to SMEAGOL NEGF library.
                          This package requires MPI.
                          Default = no
  --with-libvdwxc         Enable libvdwxc library for support of Van der Waals
                          interactions in SIRIUS.
                          This library is optinal to SIRIUS.
                          Default = no
  --with-libvori          Enable libvori for the Voronoi integration and the
                          BQB compressed trajectory format.
                          Default = install
  --with-spglib           Enable the spg library for symmetry groups detection.
                          This package depends on CMake.
                          Default = install
  --with-dftd4            Enable the standalone DFTD4 package by Grimme for the
                          DFT-D4 dispersion correction method.
                          This package requires CMake.
                          Default = no
  --with-tblite           Enable the tblite package by Grimme for GFN-xtb and
                          DFT-D4 methods, bundled with multicharge, mctc-lib,
                          mstore, s-dftd3, and toml-f libraries as backends.
                          If tblite is used, standalone DFTD4 package specified
                          by --with-dftd4 will not be used.
                          This package requires CMake.
                          Default = no
  --with-sirius           Enable interface to the plane wave SIRIUS library.
                          This package requires GSL, libspg, ELPA, ScaLAPACK,
                          HDF5, Libxc and pugixml, and libvdwxc is optinal.
                          Default = install
  --with-pugixml          Enable pugixml library for XML parsing.
                          This library is required by SIRIUS.
                          Default = no
  --with-spla             Enable the Specialized Parallel Linear Algebra (SPLA)
                          library.
                          This library is required by SIRIUS and is optional
                          for GPU support.
                          Default = no
  --with-spfft            Enable SpFFT for sparse Fourier Transform.
                          This library is required by SIRIUS.
                          Default = no
  --with-trexio           Enable the trexio library for TREXIO file format.
                          Default = no
  --with-mcl              Install MCL library for MiMiC with toolchain.
                          Default = no
  --with-dbcsr            Install DBCSR library with toolchain.
                          Default = install
  --with-cusolvermp       NVIDIA cusolverMp: CUDA library for distributed dense
                          linear algebra.
                          Default = no

FURTHER INSTRUCTIONS

An alternative way to prepare the dependencies and environment for CP2K is to
use Spack along with podman, which is handled by make_cp2k.sh, another script
that has nothing to do with the current toolchain script.

This toolchain script does not use system-dependent package managers (e.g. apt,
yum, dnf) for finding resource and configuring or installing packages. However,
the script assumes that several prerequisites including (but not limited to)
wget, bzip2 and make are present, which should be ready from package managers.
Prior to executing the toolchain script, the install_requirements.sh script in
the toolchain directory can help collecting them.

All packages to be installed locally will be downloaded and built inside
./build, and then installed into package specific directories inside ./install.

Both ./build and ./install are safe to delete, as they contain only the files
and directories that are generated by this script. However, once all the
packages are installed and you compiled CP2K then you must keep ./install in
exactly the same location as it was first created, as it contains tools and
libraries your version of CP2K binary will depend on.

It should be safe to terminate running of this script in the middle of a
build process. The script will know if a package has been successfully
installed, and will just carry on and recompile and install the last
package it is working on. This is true even if you lose the content of
the entire ./build directory.

For HPC users who wish to install toolchain dependencies and CP2K on public
supercomputer clusters for oneself: as this is a complicated process, it is
strongly advised to contact local system administrators or managers for timely,
specific assistance. Nevertheless, some hints and observations may be useful:
(1) Generally root or sudo power is not necessary, and a convenient directory
with read and write permission as well as sufficient disk space should be okay
when installing toolchain and CP2K for a single user.
(2) The server is very likely to have multiple compilers, MPI libraries, math
libraries and other packages that are managed by module systems, such as LMod
and Environment Modules. Users can load or unload modules to control active
environment variables and paths in runtime without conflicts. It is recommended
to check for available modules (with "module avail", "module show" or similar
commands) beforehand, and activate desired compatible packages when running the
toolchain script with "--with-PKG=system" options to avoid repeated labour.
(3) If no internet connection is available for downloading packages from public
resources on the server, an offline installation of toolchain and CP2K may be
carried out by downloading all packages elsewhere, transferring them to server
and placing them under the ./build directory. The toolchain script will not
attempt to download packages if they are already present in the build directory
with filenames reflecting the correct versions.
(4) An important common feature of clusters is the distinction of node types:
"login node", where users log in and perform tasks with low workload; and
"compute node", where resource-intensive computation jobs are carried out.
They may be hosted on separate machines, and their hardware specifications
(CPU, RAM, disk space, etc.) may be similar or different. Therefore, care must
be taken especially for the latter case; for instance, running toolchain script
on login node with the option of "--target-cpu=native" (which is default too if
omitted) and executing CP2K on compute node afterwards may result in unknown
behaviors due to discrepancies in CPU architectures and supported instruction
sets. The option with best portability for compiling programs at the cost of
reduced (non-optimal) performance is "--target-cpu=generic".
(5) Again, be careful about the environment if CP2K is to be executed with job
submission scripts to the job queue system handling resource allocation. Active
environment variables and paths on the login node (by loading modules, sourcing
scripts, editing ~/.bashrc or /etc/profile files, etc.) may NOT be active on
the compute node where CP2K actually runs, unless all appropriate commands are
explicitly written in the job submission script.

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
tool_list="gcc intel amd cmake ninja"
mpi_list="mpich openmpi intelmpi"
math_list="mkl acml openblas"
lib_list="fftw libint libxc libxsmm cosma scalapack elpa dbcsr
          cusolvermp plumed spfft spla gsl spglib hdf5 libvdwxc sirius
          libvori libtorch deepmd ace dftd4 tblite pugixml libsmeagol
          trexio greenx gmp mcl"
package_list="${tool_list} ${mpi_list} ${math_list} ${lib_list}"
# ------------------------------------------------------------------------

# first set everything to __DONTUSE__
for ii in ${package_list}; do
  eval "with_${ii}=__DONTUSE__"
done

# ------------------------------------------------------------------------
# Work out default settings
# ------------------------------------------------------------------------

# tools to turn on by default:
with_gcc="__SYSTEM__"

# libs to turn on by default:
with_dbcsr="__INSTALL__"
with_fftw="__INSTALL__"
with_libint="__INSTALL__"
with_libxsmm="__INSTALL__"
with_libxc="__INSTALL__"
with_scalapack="__INSTALL__"
with_sirius="__INSTALL__"
with_gsl="__DONTUSE__"
with_spglib="__INSTALL__"
with_hdf5="__DONTUSE__"
with_trexio="__DONTUSE__"
with_elpa="__INSTALL__"
with_cusolvermp="__DONTUSE__"
with_libvdwxc="__DONTUSE__"
with_spfft="__DONTUSE__"
with_spla="__DONTUSE__"
with_cosma="__INSTALL__"
with_libvori="__INSTALL__"
with_libtorch="__DONTUSE__"
with_ninja="__DONTUSE__"
with_dftd4="__DONTUSE__"
with_tblite="__DONTUSE__"
with_libsmeagol="__DONTUSE__"
with_mcl="__DONTUSE__"

# the math and mpi libraries are chosen by their respective modes.
# default math library settings, MATH_MODE picks the math library
# to use, and with_* defines the default method of installation if it
# is picked. For non-CRAY systems defaults to mkl if $MKLROOT is
# available, otherwise defaults to openblas
if [ -n "${MKLROOT}" ]; then
  export MATH_MODE="mkl"
  with_mkl="__SYSTEM__"
else
  export MATH_MODE="openblas"
fi
with_acml="__SYSTEM__"
with_openblas="__INSTALL__"

# for MPI, we try to detect system MPI variant
if (command -v mpiexec > /dev/null 2>&1); then
  # check if we are dealing with openmpi, mpich or intelmpi
  if (mpiexec --version 2>&1 | grep -s -q "HYDRA"); then
    echo "MPI is detected and it appears to be MPICH"
    export MPI_MODE="mpich"
    with_mpich="__SYSTEM__"
  elif (mpiexec --version 2>&1 | grep -s -q "OpenRTE"); then
    echo "MPI is detected and it appears to be OpenMPI 4 (or older)"
    export MPI_MODE="openmpi"
    with_openmpi="__SYSTEM__"
  elif (mpiexec --version 2>&1 | grep -s -q "Open MPI"); then
    echo "MPI is detected and it appears to be OpenMPI 5"
    export MPI_MODE="openmpi"
    with_openmpi="__SYSTEM__"
  elif (mpiexec --version 2>&1 | grep -s -q "Intel"); then
    echo "MPI is detected and it appears to be Intel MPI"
    with_gcc="__DONTUSE__"
    with_amd="__DONTUSE__"
    with_intel="__SYSTEM__"
    with_intelmpi="__SYSTEM__"
    export MPI_MODE="intelmpi"
  else # default to mpich
    echo "MPI is detected and defaults to MPICH"
    export MPI_MODE="mpich"
    with_mpich="__SYSTEM__"
  fi
else
  if [ -n "${CRAY_LD_LIBRARY_PATH}" ]; then
    echo "Cray Linux Environment (CLE) is detected with no MPI"
  else
    echo "No MPI installation detected. (Ignore this message if
    a fresh MPI installation is requested.)"
  fi
  export MPI_MODE="no"
fi

# default enable options
dry_run="__FALSE__"
enable_tsan="__FALSE__"
enable_opencl="__FALSE__"
enable_cuda="__FALSE__"
enable_hip="__FALSE__"
export with_ifx="no"
export GPUVER="no"
export MPICH_DEVICE="ch4"
export TARGET_CPU="native"
NPROCS_OVERWRITE="$(get_nprocs)"
export NPROCS_OVERWRITE

# default for libint
export LIBINT_LMAX="5"

# default for log file dump size
export LOG_LINES="200"

# defaults for CRAY Linux Environment
if [ -n "${CRAY_LD_LIBRARY_PATH}" ]; then
  enable_cray="__TRUE__"
  export MATH_MODE="cray"
  # Default MPI used by CLE is assumed to be MPICH, in any case
  # do not use the installers for the MPI libraries
  with_mpich="__DONTUSE__"
  with_openmpi="__DONTUSE__"
  with_intelmpi="__DONTUSE__"
  export MPI_MODE="mpich"
  # set default value for some installers appropriate for CLE
  with_gcc="__DONTUSE__"
  with_amd="__DONTUSE__"
  with_intel="__DONTUSE__"
  with_fftw="__SYSTEM__"
  with_scalapack="__DONTUSE__"
else
  enable_cray="__FALSE__"
fi

# ------------------------------------------------------------------------
# Parse user options
# ------------------------------------------------------------------------
echo "Toolchain script received the following options:"
echo "  ${TOOLCHAIN_OPTIONS}"
echo "Parsing options and resolving conflicts..."
while [ $# -ge 1 ]; do
  case ${1} in
    -j)
      case "${2}" in
        -*)
          NPROCS_OVERWRITE="$(get_nprocs)"
          export NPROCS_OVERWRITE
          ;;
        [0-9]*)
          shift
          export NPROCS_OVERWRITE="${1}"
          ;;
        *)
          report_error ${LINENO} "Non-integer argument ${2} for -j flag found."
          exit 1
          ;;
      esac
      ;;
    -j[0-9]*)
      export NPROCS_OVERWRITE="${1#-j}"
      ;;
    --no-check-certificate)
      export DOWNLOADER_FLAGS="--no-check-certificate"
      ;;
    --install-all)
      # set all package to the default installation status
      for ii in ${package_list}; do
        if [ "${ii}" != "intel" ] &&
          [ "${ii}" != "intelmpi" ] &&
          [ "${ii}" != "amd" ]; then
          eval "with_${ii}=__INSTALL__"
        fi
      done
      # Use MPICH as default
      export MPI_MODE="mpich"
      ;;
    --mpi-mode=*)
      user_input="${1#*=}"
      case "$user_input" in
        mpich)
          export MPI_MODE="mpich"
          ;;
        openmpi)
          export MPI_MODE="openmpi"
          ;;
        intelmpi)
          export MPI_MODE="intelmpi"
          ;;
        no)
          export MPI_MODE="no"
          ;;
        *)
          report_error ${LINENO} "Invalid value for --mpi-mode found."
          echo "Currently only one of the following options is supported:
            openmpi, mpich, intelmpi.
Otherwise use option no."
          exit 1
          ;;
      esac
      ;;
    --math-mode=*)
      user_input="${1#*=}"
      case "$user_input" in
        cray)
          export MATH_MODE="cray"
          ;;
        mkl)
          export MATH_MODE="mkl"
          ;;
        acml)
          export MATH_MODE="acml"
          ;;
        openblas)
          export MATH_MODE="openblas"
          ;;
        *)
          report_error ${LINENO} "Invalid value for --math-mode found."
          echo "Currently only one of the following options is supported:
            mkl, acml, openblas, cray."
          exit 1
          ;;
      esac
      ;;
    --gpu-ver=*)
      user_input="${1#*=}"
      case "${user_input}" in
        K20X | K40 | K80 | P100 | V100 | A100 | H100 | A40 | Mi50 | Mi100 | Mi250 | no)
          export GPUVER="${user_input}"
          ;;
        *)
          report_error ${LINENO} "Invalid value for --gpu-ver found."
          echo "Currently only one of the following options is supported:
            K20X, K40, K80, P100, V100, A100, H100, A40, Mi50, Mi100, Mi250.
Otherwise use option no."
          exit 1
          ;;
      esac
      ;;
    --target-cpu=*)
      user_input="${1#*=}"
      export TARGET_CPU="${user_input}"
      ;;
    --log-lines=*)
      user_input="${1#*=}"
      case "${user_input}" in
        [0-9]*)
          export LOG_LINES="${user_input}"
          ;;
        *)
          report_error ${LINENO} "Non-integer ${user_input} for --log-lines."
          exit 1
          ;;
      esac
      ;;
    --libint-lmax=*)
      user_input="${1#*=}"
      export LIBINT_LMAX="${user_input}"
      ;;
    --dry-run)
      dry_run="__TRUE__"
      ;;
    --enable-tsan*)
      enable_tsan=$(read_enable "${1}")
      if [ "${enable_tsan}" = "__INVALID__" ]; then
        report_error "invalid value for --enable-tsan, please use yes or no"
        exit 1
      fi
      ;;
    --enable-cuda*)
      enable_cuda=$(read_enable "${1}")
      if [ "${enable_cuda}" = "__INVALID__" ]; then
        report_error "invalid value for --enable-cuda, please use yes or no"
        exit 1
      fi
      ;;
    --enable-hip*)
      enable_hip=$(read_enable "${1}")
      if [ "${enable_hip}" = "__INVALID__" ]; then
        report_error "invalid value for --enable-hip, please use yes or no"
        exit 1
      fi
      ;;
    --enable-opencl*)
      enable_opencl=$(read_enable "${1}")
      if [ "${enable_opencl}" = "__INVALID__" ]; then
        report_error "invalid value for --enable-opencl, please use yes or no"
        exit 1
      fi
      ;;
    --enable-cray*)
      enable_cray=$(read_enable "${1}")
      if [ "${enable_cray}" = "__INVALID__" ]; then
        report_error "invalid value for --enable-cray, please use yes or no"
        exit 1
      fi
      ;;
    --with-gcc*)
      with_gcc=$(read_with "${1}")
      ;;
    --with-cmake*)
      with_cmake=$(read_with "${1}")
      ;;
    --with-ninja*)
      with_ninja=$(read_with "${1}")
      ;;
    --with-mpich-device=*)
      user_input="${1#*=}"
      export MPICH_DEVICE="${user_input}"
      export MPI_MODE="mpich"
      ;;
    --with-mpich*)
      with_mpich=$(read_with "${1}")
      if [ "${with_mpich}" != "__DONTUSE__" ]; then
        export MPI_MODE="mpich"
      fi
      ;;
    --with-openmpi*)
      with_openmpi=$(read_with "${1}")
      if [ "${with_openmpi}" != "__DONTUSE__" ]; then
        export MPI_MODE="openmpi"
      fi
      ;;
    --with-intelmpi*)
      with_intelmpi=$(read_with "${1}" "__SYSTEM__")
      if [ "${with_intelmpi}" != "__DONTUSE__" ]; then
        export MPI_MODE="intelmpi"
      fi
      ;;
    --with-amd*)
      with_amd=$(read_with "${1}" "__SYSTEM__")
      ;;
    --with-ifx*)
      with_ifx=$(read_with "${1}" "yes")
      ;;
    --with-intel*)
      with_intel=$(read_with "${1}" "__SYSTEM__")
      ;;
    --with-libint*)
      with_libint=$(read_with "${1}")
      ;;
    --with-libxc*)
      with_libxc=$(read_with "${1}")
      ;;
    --with-fftw*)
      with_fftw=$(read_with "${1}")
      ;;
    --with-mkl*)
      with_mkl=$(read_with "${1}" "__SYSTEM__")
      if [ "${with_mkl}" != "__DONTUSE__" ]; then
        export MATH_MODE="mkl"
      fi
      ;;
    --with-acml*)
      with_acml=$(read_with "${1}")
      if [ "${with_acml}" != "__DONTUSE__" ]; then
        export MATH_MODE="acml"
      fi
      ;;
    --with-openblas*)
      with_openblas=$(read_with "${1}")
      if [ "${with_openblas}" != "__DONTUSE__" ]; then
        export MATH_MODE="openblas"
      fi
      ;;
    --with-scalapack*)
      with_scalapack=$(read_with "${1}")
      ;;
    --with-libxsmm*)
      with_libxsmm=$(read_with "${1}")
      ;;
    --with-elpa*)
      with_elpa=$(read_with "${1}")
      ;;
    --with-cusolvermp*)
      with_cusolvermp=$(read_with "${1}")
      ;;
    --with-deepmd*)
      with_deepmd=$(read_with "${1}")
      ;;
    --with-ace*)
      with_ace=$(read_with "${1}")
      ;;
    --with-plumed*)
      with_plumed=$(read_with "${1}")
      ;;
    --with-sirius*)
      with_sirius=$(read_with "${1}")
      ;;
    --with-pugixml*)
      with_pugixml=$(read_with "${1}")
      ;;
    --with-gsl*)
      with_gsl=$(read_with "${1}")
      ;;
    --with-spglib*)
      with_spglib=$(read_with "${1}")
      ;;
    --with-hdf5*)
      with_hdf5=$(read_with "${1}")
      ;;
    --with-libvdwxc*)
      with_libvdwxc=$(read_with "${1}")
      ;;
    --with-spfft*)
      with_spfft=$(read_with "${1}")
      ;;
    --with-cosma*)
      with_cosma=$(read_with "${1}")
      ;;
    --with-libvori*)
      with_libvori=$(read_with "${1}")
      ;;
    --with-libtorch*)
      with_libtorch=$(read_with "${1}")
      ;;
    --with-spla*)
      with_spla=$(read_with "${1}")
      ;;
    --with-dftd4*)
      with_dftd4=$(read_with "${1}")
      ;;
    --with-tblite*)
      with_tblite=$(read_with "${1}")
      ;;
    --with-libsmeagol*)
      with_libsmeagol=$(read_with "${1}")
      ;;
    --with-trexio*)
      with_trexio=$(read_with "${1}")
      ;;
    --with-greenx*)
      with_greenx=$(read_with "${1}")
      ;;
    --with-gmp*)
      with_gmp=$(read_with "${1}")
      ;;
    --with-dbcsr*)
      with_dbcsr=$(read_with "${1}")
      ;;
    --with-mcl*)
      with_mcl=$(read_with "${1}")
      ;;
    --help*)
      show_help
      exit 0
      ;;
    -h*)
      show_help
      exit 0
      ;;
    *)
      report_error ${LINENO} "Unknown flag: ${1}
See help message of this script produced by --help option for supported ones."
      exit 1
      ;;
  esac
  shift
done

# consolidate settings after user input
export ENABLE_TSAN="${enable_tsan}"
export ENABLE_CUDA="${enable_cuda}"
export ENABLE_HIP="${enable_hip}"
export ENABLE_OPENCL="${enable_opencl}"
export ENABLE_CRAY="${enable_cray}"

# ------------------------------------------------------------------------
# Check and solve known conflicts before installations proceed
# ------------------------------------------------------------------------
# Compiler conflicts
if [ "${with_gcc}" = "__INSTALL__" ]; then
  if [ "${with_intel}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "With Intel compiler, do not install GNU compiler."
    with_gcc="__SYSTEM__"
  elif [ "${with_amd}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "With AMD compiler, do not install GNU compiler."
    with_gcc="__SYSTEM__"
  fi
fi
if [ "${with_amd}" != "__DONTUSE__" ]; then
  if [ "${with_intel}" != "__DONTUSE__" ]; then
    report_error ${LINENO} "The AMD and Intel compilers can't be used together."
    exit 1
  fi
fi
# MPI library conflicts
if [ "${MPI_MODE}" = "no" ]; then
  if [ "${with_scalapack}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "Not using MPI, so ScaLAPACK is disabled."
    with_scalapack="__DONTUSE__"
  fi
  if [ "${with_elpa}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "Not using MPI, so ELPA is disabled."
    with_elpa="__DONTUSE__"
  fi
  if [ "${with_plumed}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "Not using MPI, so PLUMED is disabled."
    with_plumed="__DONTUSE__"
  fi
  if [ "${with_libsmeagol}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "Not using MPI, so libsmeagol is disabled."
    with_libsmeagol="__DONTUSE__"
  fi
  if [ "${with_sirius}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "Not using MPI, so SIRIUS is disabled."
    with_sirius="__DONTUSE__"
  fi
  if [ "${with_spfft}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "Not using MPI, so SpFFT is disabled."
    with_spfft="__DONTUSE__"
  fi
  if [ "${with_spla}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "Not using MPI, so SpLA is disabled."
    with_spla="__DONTUSE__"
  fi
  if [ "${with_cosma}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "Not using MPI, so COSMA is disabled."
    with_cosma="__DONTUSE__"
  fi
  if [ "${with_mcl}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "Not using MPI, so MCL is disabled."
    with_mcl="__DONTUSE__"
  fi
else
  # if gcc is installed, then mpi needs to be installed too.
  if [ "${with_gcc}" = "__INSTALL__" ]; then
    report_warning ${LINENO} "When installing GNU compiler, Intel compiler and
      Intel MPI are disabled, and MPI library is freshly installed."
    with_intel="__DONTUSE__"
    with_intelmpi="__DONTUSE__"
    case ${MPI_MODE} in
      mpich)
        with_mpich="__INSTALL__"
        with_openmpi="__DONTUSE__"
        ;;
      openmpi)
        with_mpich="__DONTUSE__"
        with_openmpi="__INSTALL__"
        ;;
      intelmpi)
        report_error ${LINENO} "Incompatible --mpi-mode=intelmpi found."
        exit 1
        ;;
    esac
  fi
  # Enable only one MPI implementation.
  case ${MPI_MODE} in
    mpich)
      with_openmpi="__DONTUSE__"
      with_intelmpi="__DONTUSE__"
      if [ "${with_mpich}" = "__DONTUSE__" ]; then
        with_mpich="__INSTALL__"
      fi
      ;;
    openmpi)
      with_mpich="__DONTUSE__"
      with_intelmpi="__DONTUSE__"
      if [ "${with_openmpi}" = "__DONTUSE__" ]; then
        with_openmpi="__INSTALL__"
      fi
      ;;
    intelmpi)
      with_mpich="__DONTUSE__"
      with_openmpi="__DONTUSE__"
      if [ "${with_intelmpi}" = "__DONTUSE__" ]; then
        report_error ${LINENO} "While --mpi-mode=intelmpi is set, no Intel MPI
could be found or linked in the system, and installation by toolchain is not
supported. Please install manually and check executable path before rerunning."
        exit 1
      fi
      ;;
  esac
fi

# If CUDA or HIP are enabled, make sure the GPU version has been defined.
if [ "${ENABLE_CUDA}" = "__TRUE__" ] || [ "${ENABLE_HIP}" = "__TRUE__" ]; then
  if [ "${GPUVER}" = "no" ]; then
    report_error ${LINENO} "Either CUDA or HIP is enabled, but --gpu-ver is not
set to one of the known architectures. See help message of this script produced
by --help option for supported ones."
    exit 1
  fi
fi

# If OpenCL is enabled, make sure LIBXSMM is enabled as well.
if [ "${ENABLE_OPENCL}" = "__TRUE__" ]; then
  if [ "${with_libxsmm}" = "__DONTUSE__" ]; then
    report_warning ${LINENO} "When enabling OpenCL, libxsmm is needed."
    with_libxsmm="__INSTALL__"
  fi
fi

# Since tblite includes dftd4, a separate dftd4 is not needed.
if [ "${with_tblite}" != "__DONTUSE__" ]; then
  if [ "${with_dftd4}" != "__DONTUSE__" ]; then
    report_warning ${LINENO} "Since tblite includes dft-d4, a standalone dft-d4
package will not be used separately."
    with_dftd4="__DONTUSE__"
  fi
fi

# Several packages require cmake.
if [ "${with_spglib}" = "__INSTALL__" ] ||
  [ "${with_libvori}" = "__INSTALL__" ] ||
  [ "${with_scalapack}" = "__INSTALL__" ] ||
  [ "${with_sirius}" = "__INSTALL__" ] ||
  [ "${with_pugixml}" = "__INSTALL__" ] ||
  [ "${with_cosma}" = "__INSTALL__" ] ||
  [ "${with_spfft}" = "__INSTALL__" ] ||
  [ "${with_spla}" = "__INSTALL__" ] ||
  [ "${with_ninja}" = "__INSTALL__" ] ||
  [ "${with_greenx}" = "__INSTALL__" ] ||
  [ "${with_dftd4}" = "__INSTALL__" ] ||
  [ "${with_mcl}" = "__INSTALL__" ] ||
  [ "${with_tblite}" = "__INSTALL__" ]; then
  if [ "${with_cmake}" = "__DONTUSE__" ]; then
    report_warning ${LINENO} "Installing one of the packages requires CMake but
CMake is not found in system, so a new copy of CMake will be installed first."
    with_cmake="__INSTALL__"
  fi
fi

# SIRIUS dependencies
if [ "${with_sirius}" = "__INSTALL__" ]; then
  [ "${with_spfft}" = "__DONTUSE__" ] && with_spfft="__INSTALL__"
  [ "${with_spla}" = "__DONTUSE__" ] && with_spla="__INSTALL__"
  [ "${with_gsl}" = "__DONTUSE__" ] && with_gsl="__INSTALL__"
  [ "${with_libxc}" = "__DONTUSE__" ] && with_libxc="__INSTALL__"
  [ "${with_fftw}" = "__DONTUSE__" ] && with_fftw="__INSTALL__"
  [ "${with_spglib}" = "__DONTUSE__" ] && with_spglib="__INSTALL__"
  [ "${with_hdf5}" = "__DONTUSE__" ] && with_hdf5="__INSTALL__"
  [ "${with_libvdwxc}" = "__DONTUSE__" ] && with_libvdwxc="__INSTALL__"
  [ "${with_cosma}" = "__DONTUSE__" ] && with_cosma="__INSTALL__"
  [ "${with_pugixml}" = "__DONTUSE__" ] && with_pugixml="__INSTALL__"
elif [ "${with_sirius}" = "__DONTUSE__" ]; then
  with_pugixml="__DONTUSE__"
  with_spfft="__DONTUSE__"
  with_libvdwxc="__DONTUSE__"
  [ "${GPUVER}" = "no" ] && with_spla="__DONTUSE__"
fi

if [ "${with_trexio}" = "__INSTALL__" ]; then
  [ "${with_hdf5}" = "__DONTUSE__" ] && with_hdf5="__INSTALL__"
fi

if [ "${with_plumed}" = "__INSTALL__" ]; then
  [ "${with_gsl}" = "__DONTUSE__" ] && with_gsl="__INSTALL__"
  [ "${with_fftw}" = "__DONTUSE__" ] && with_fftw="__INSTALL__"
fi

if [ "${with_deepmd}" = "__INSTALL__" ]; then
  [ "${with_libtorch}" = "__DONTUSE__" ] && with_libtorch="__INSTALL__"
fi

# ------------------------------------------------------------------------
# Preliminaries
# ------------------------------------------------------------------------

mkdir -p "${INSTALLDIR}"

# Select the correct compute number based on the GPU architecture
case ${GPUVER} in
  K20X)
    export ARCH_NUM="35"
    ;;
  K40)
    export ARCH_NUM="35"
    ;;
  K80)
    export ARCH_NUM="37"
    ;;
  P100)
    export ARCH_NUM="60"
    ;;
  V100)
    export ARCH_NUM="70"
    ;;
  A100)
    export ARCH_NUM="80"
    ;;
  A40)
    export ARCH_NUM="86"
    ;;
  H100)
    export ARCH_NUM="90"
    ;;
  Mi50)
    # TODO: export ARCH_NUM=
    ;;
  Mi100)
    # TODO: export ARCH_NUM=
    ;;
  Mi250)
    # TODO: export ARCH_NUM=
    ;;
  no)
    export ARCH_NUM="no"
    ;;
  *)
    report_error ${LINENO} "Invalid value for --gpu-ver found."
    echo "Currently only one of the following options is supported:
      K20X, K40, K80, P100, V100, A100, H100, A40, Mi50, Mi100, Mi250.
Otherwise use option no."
    exit 1
    ;;
esac

# variables used for generating cp2k ARCH file
export CP_DFLAGS=""
export CP_LIBS=""
export CP_CFLAGS=""
export CP_LDFLAGS="-Wl,--enable-new-dtags"

# ------------------------------------------------------------------------
# Special settings for CRAY Linux Environment (CLE)
# TODO: CLE should be handle like gcc or Intel using a with_cray flag and
#       this section should be moved to a separate file install_cray.
# ------------------------------------------------------------------------
if [ "${ENABLE_CRAY}" = "__TRUE__" ]; then
  echo "------------------------------------------------------------------------"
  echo "CRAY Linux Environment (CLE) is detected"
  echo "------------------------------------------------------------------------"
  # add cray paths to system search path
  export LIB_PATHS="CRAY_LD_LIBRARY_PATH ${LIB_PATHS}"
  # set compilers to CLE wrappers
  check_command cc
  check_command ftn
  check_command CC
  export CC="cc"
  export CXX="CC"
  export FC="ftn"
  export F90="${FC}"
  export F77="${FC}"
  export MPICC="${CC}"
  export MPICXX="${CXX}"
  export MPIFC="${FC}"
  export MPIFORT="${MPIFC}"
  export MPIF77="${MPIFC}"
  case $MPI_MODE in
    mpich)
      if [ -n "$MPICH_DIR" ]; then
        cray_mpich_include_path="$MPICH_DIR/include"
        cray_mpich_lib_path="$MPICH_DIR/lib"
        export INCLUDE_PATHS="$INCLUDE_PATHS $cray_mpich_include_path"
        export LIB_PATHS="$LIB_PATHS $cray_mpich_lib_path"
      fi
      if [ "$with_mpich" = "__DONTUSE__" ]; then
        add_include_from_paths MPI_CFLAGS "mpi.h" "$INCLUDE_PATHS"
        add_include_from_paths MPI_LDFLAGS "libmpi.*" "$LIB_PATHS"
        export MPI_CFLAGS
        export MPI_LDFLAGS
        export MPI_LIBS=" "
        export CP_DFLAGS="${CP_DFLAGS} IF_MPI(-D__parallel|)"
      fi
      ;;
    openmpi)
      if [ "$with_openmpi" = "__DONTUSE__" ]; then
        add_include_from_paths MPI_CFLAGS "mpi.h" "$INCLUDE_PATHS"
        add_include_from_paths MPI_LDFLAGS "libmpi.*" "$LIB_PATHS"
        export MPI_CFLAGS
        export MPI_LDFLAGS
        export MPI_LIBS="-lmpi -lmpi_cxx"
        export CP_DFLAGS="${CP_DFLAGS} IF_MPI(-D__parallel|)"
      fi
      ;;
    intelmpi)
      if [ "$with_intelmpi" = "__DONTUSE__" ]; then
        with_gcc="__DONTUSE__"
        with_intel="__SYSTEM__"
        add_include_from_paths MPI_CFLAGS "mpi.h" "$INCLUDE_PATHS"
        add_include_from_paths MPI_LDFLAGS "libmpi.*" "$LIB_PATHS"
        export MPI_CFLAGS
        export MPI_LDFLAGS
        export MPI_LIBS="-lmpi -lmpi_cxx"
        export CP_DFLAGS="${CP_DFLAGS} IF_MPI(-D__parallel|)"
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

# Write head of setup file
cat << EOF > "$SETUPFILE"
#!/bin/bash
source "${SCRIPTDIR}/tool_kit.sh"
export CP2K_TOOLCHAIN_OPTIONS="${TOOLCHAIN_OPTIONS}"
EOF

# Write toolchain environment
write_toolchain_env "${INSTALLDIR}"

# Write toolchain config
echo "tool_list=\"${tool_list}\"" > "${INSTALLDIR}"/toolchain.conf
echo "dry_run=\"${dry_run}\"" >> "${INSTALLDIR}"/toolchain.conf
for ii in ${package_list}; do
  install_mode=$(eval "echo \${with_${ii}}")
  echo "with_${ii}=\"${install_mode}\"" >> "${INSTALLDIR}"/toolchain.conf
done

# ------------------------------------------------------------------------
# Build packages unless dry-run mode is enabled.
# ------------------------------------------------------------------------
if [ "${dry_run}" = "__TRUE__" ]; then
  printf "With --dry-run option, this script concludes with a report.\n"
  printf "The setup, toolchain env and conf files are written to ./install.\n"
  printf "System specifications:\n"
  printf '   -%-20s = %s\n' "j" "${NPROCS_OVERWRITE}"
  printf '  --%-20s = %s\n' "target-cpu" "${TARGET_CPU}"
  printf '  --%-20s = %s\n' "gpu-ver" "${GPUVER}"
  printf '  --%-20s = %s\n' "mpi-mode" "${MPI_MODE}"
  printf '  --%-20s = %s\n' "math-mode" "${MATH_MODE}"
  printf '  --%-20s = %s\n' "enable-tsan" "${enable_tsan}"
  printf '  --%-20s = %s\n' "enable-cuda" "${enable_cuda}"
  printf '  --%-20s = %s\n' "enable-hip" "${enable_hip}"
  printf '  --%-20s = %s\n' "enable-opencl" "${enable_opencl}"
  printf '  --%-20s = %s\n' "enable-cray" "${enable_cray}"
  printf "List of effective settings after resolving package conflicts:\n"
  for ii in ${package_list}; do
    install_mode=$(eval "echo \${with_${ii}}")
    printf '  --with-%-15s = %s\n' "${ii}" "${install_mode}"
  done
else
  echo "Options have been parsed successfully."
  echo "Compiling with ${NPROCS_OVERWRITE} processes for target ${TARGET_CPU}."
  echo "# Leak suppressions" > "${INSTALLDIR}"/lsan.supp
  ./scripts/stage0/install_stage0.sh
  ./scripts/stage1/install_stage1.sh
  ./scripts/stage2/install_stage2.sh
  ./scripts/stage3/install_stage3.sh
  ./scripts/stage4/install_stage4.sh
  ./scripts/stage5/install_stage5.sh
  ./scripts/stage6/install_stage6.sh
  ./scripts/stage7/install_stage7.sh
  ./scripts/stage8/install_stage8.sh
  ./scripts/stage9/install_stage9.sh
fi

# Generate CMake options
./scripts/generate_cmake_options.sh
