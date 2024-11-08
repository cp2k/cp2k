# libxsmm is designed for x86_64 architectures only, see project's README
%bcond libxsmm %[ "%{_arch}" == "x86_64" ]

# Disable LTO due to https://bugzilla.redhat.com/show_bug.cgi?id=2243158
%global _lto_cflags %nil

Name:          cp2k
Version:       0.0.0
Release:       %autorelease
Summary:       Ab Initio Molecular Dynamics
License:       GPL-2.0-or-later
URL:           https://www.cp2k.org/
Source0:       https://github.com/cp2k/cp2k/releases/download/v%{version}/cp2k-%{version}.tar.bz2

# Drop 32bit architectures
# Flaky MPI issues on s390x, and upstream do not officially support it yet
# https://github.com/cp2k/cp2k/issues/3362
ExcludeArch:   %{ix86} s390x

# Build dependencies
BuildRequires: cmake
BuildRequires: gcc
BuildRequires: gcc-gfortran
BuildRequires: gcc-c++
BuildRequires: ninja-build
BuildRequires: python3-fypp
# Project dependencies
BuildRequires: flexiblas-devel
BuildRequires: cmake(DBCSR)
BuildRequires: pkgconfig(fftw3)
%if %{with libxsmm}
BuildRequires: pkgconfig(libxsmm)
%endif
# cmake(Libxc) technically fails
# https://github.com/cp2k/cp2k/issues/3767
BuildRequires: libxc-devel
BuildRequires: cmake(Spglib)
# Test dependencies
BuildRequires: python3

Requires:      %{name}-common = %{version}-%{release}

%global _description %{expand:
CP2K is a freely available (GPL) program, written in Fortran 95, to
perform atomistic and molecular simulations of solid state, liquid,
molecular and biological systems. It provides a general framework for
different methods such as e.g. density functional theory (DFT) using a
mixed Gaussian and plane waves approach (GPW), and classical pair and
many-body potentials.

CP2K does not implement Car-Parinello Molecular Dynamics (CPMD).}

%description
%{_description}

This package contains the non-MPI single process and multi-threaded versions.

%package common
Summary:       Molecular simulations software - common files
BuildArch:     noarch

%description common
%{_description}

This package contains the documentation and the manual.

%package devel
Summary:        Development files for %{name}

Requires:       %{name}%{?_isa} = %{version}-%{release}

%description devel
The %{name}-devel package contains libraries and header files for
developing applications that use %{name}.

%package openmpi
Summary:        Molecular simulations software - openmpi version
BuildRequires:  openmpi-devel
BuildRequires:  blacs-openmpi-devel
BuildRequires:  dbcsr-openmpi-devel
BuildRequires:  scalapack-openmpi-devel

Requires:       %{name}-common = %{version}-%{release}

%description openmpi
%{_description}

This package contains the parallel single- and multi-threaded versions
using OpenMPI.

%package openmpi-devel
Summary:        Development files for %{name}

Requires:       %{name}-openmpi%{?_isa} = %{version}-%{release}

%description openmpi-devel
The %{name}-devel package contains libraries and header files for
developing applications that use %{name}.

%package mpich
Summary:        Molecular simulations software - mpich version
BuildRequires:  mpich-devel
BuildRequires:  blacs-mpich-devel
BuildRequires:  dbcsr-mpich-devel
BuildRequires:  scalapack-mpich-devel

Requires:       %{name}-common = %{version}-%{release}

%description mpich
%{_description}

This package contains the parallel single- and multi-threaded versions
using mpich.

%package mpich-devel
Summary:        Development files for %{name}

Requires:       %{name}-mpich%{?_isa} = %{version}-%{release}

%description mpich-devel
The %{name}-devel package contains libraries and header files for
developing applications that use %{name}.


%prep
%autosetup -p1
rm tools/build_utils/fypp
rm -r exts/dbcsr

# $MPI_SUFFIX will be evaluated in the loops below, set by mpi modules
%global _vpath_builddir %{_vendor}-%{_target_os}-build${MPI_SUFFIX:-_serial}


%build
cmake_common_args=(
  "-G Ninja"
  "-DCP2K_DEBUG_MODE:BOOL=OFF"
  "-DCP2K_BLAS_VENDOR:STRING=FlexiBLAS"
  "-DCP2K_USE_STATIC_BLAS:BOOL=OFF"
  # Dependencies equivalent with Default
  "-DCP2K_USE_FFTW3:BOOL=ON"
  "-DCP2K_USE_COSMA:BOOL=OFF"  # Not packaged
  "-DCP2K_USE_LIBXSMM:BOOL=%{?with_libxsmm:ON}%{?without_libxsmm:OFF}"
  "-DCP2K_USE_LIBXC:BOOL=ON"
  "-DCP2K_USE_LIBINT2:BOOL=OFF"  # Detection is broken
  "-DCP2K_USE_SPGLIB:BOOL=ON"
)
for mpi in '' mpich openmpi; do
  if [ -n "$mpi" ]; then
    module load mpi/${mpi}-%{_arch}
    cmake_mpi_args=(
      "-DCMAKE_INSTALL_PREFIX:PATH=${MPI_HOME}"
      "-DCMAKE_INSTALL_Fortran_MODULES:PATH=${MPI_FORTRAN_MOD_DIR}/cp2k"
      "-DCMAKE_INSTALL_LIBDIR:PATH=lib"
      "-DCP2K_CMAKE_SUFFIX:STRING=${MPI_SUFFIX}"
      "-DCP2K_DATA_DIR:PATH=%{_datadir}/cp2k/data"
      "-DCP2K_USE_MPI_F08:BOOL=ON"
    )
  else
    cmake_mpi_args=(
      "-DCP2K_USE_MPI:BOOL=OFF"
      "-DCMAKE_INSTALL_Fortran_MODULES:PATH=%{_fmoddir}/cp2k"
    )
  fi

  %cmake \
    ${cmake_common_args[@]} \
    ${cmake_mpi_args[@]}
  %cmake_build

  [ -n "$mpi" ] && module unload mpi/${mpi}-%{_arch}
done

%install
for mpi in '' mpich openmpi; do
  [ -n "$mpi" ] && module load mpi/${mpi}-%{_arch}
  %cmake_install
  [ -n "$mpi" ] && module unload mpi/${mpi}-%{_arch}
done

# TODO: Properly separate the installation of unit tests
rm -f %{_buildrootdir}/**/%{_bindir}/*_unittest.*
rm -f %{_buildrootdir}/**/%{_libdir}/openmpi/bin/*_unittest.*
rm -f %{_buildrootdir}/**/%{_libdir}/mpich/bin/*_unittest.*


%check
export CP2K_DATA_DIR=%{buildroot}%{_datadir}/cp2k/data
# See %%_openmpi_load
export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe
test_common_args=(
  "--skip_regtests"
  "--ompthreads 2"
)
for mpi in '' mpich openmpi ; do
  if [ -n "$mpi" ]; then
    # Another module load is done inside the do_regtest.sh. will use that instead
    module load mpi/${mpi}-%{_arch}
    bindir=${MPI_BIN}
    libdir=${MPI_LIB}
    # Note, final position arguments are also here
    test_mpi_args=(
      "--mpiranks 2"
      "psmp"
    )
  else
    bindir=%{_bindir}
    libdir=%{_libdir}
    test_mpi_args=(
      "ssmp"
    )
  fi
  # Run packaged do_regtest.sh with appropriate buildroot runpaths
  # Note: Running unittests only in the spec file which are not packaged,
  # so the binary folder should point to the build directory
  env PATH=%{buildroot}${bindir}:${PATH} \
    LD_LIBRARY_PATH=%{buildroot}${libdir} \
    tests/do_regtest.py ${test_common_args[@]} %{_vpath_builddir}/bin ${test_mpi_args[@]}
  [ -n "$mpi" ] && module unload mpi/${mpi}-%{_arch}
done

%files common
%license LICENSE
%doc README.md
%{_datadir}/cp2k

%files
%{_bindir}/cp2k.ssmp
%{_bindir}/dbm_miniapp.ssmp
%{_bindir}/dumpdcd.ssmp
%{_bindir}/graph.ssmp
%{_bindir}/grid_miniapp.ssmp
%{_bindir}/xyz2dcd.ssmp
%{_libdir}/libcp2k.so.*

%files devel
%{_fmoddir}/cp2k/
%{_includedir}/cp2k/
%{_libdir}/cmake/cp2k/
%{_libdir}/libcp2k.so
%{_libdir}/pkgconfig/libcp2k.pc

%files openmpi
%{_libdir}/openmpi/bin/cp2k.psmp
%{_libdir}/openmpi/bin/dumpdcd.psmp
%{_libdir}/openmpi/bin/dbm_miniapp.psmp
%{_libdir}/openmpi/bin/graph.psmp
%{_libdir}/openmpi/bin/grid_miniapp.psmp
%{_libdir}/openmpi/bin/xyz2dcd.psmp
%{_libdir}/openmpi/lib/libcp2k.so.*

%files openmpi-devel
%{_fmoddir}/openmpi/cp2k/
%{_libdir}/openmpi/include/cp2k/
%{_libdir}/openmpi/lib/cmake/cp2k/
%{_libdir}/openmpi/lib/libcp2k.so
%{_libdir}/openmpi/lib/pkgconfig/libcp2k.pc

%files mpich
%{_libdir}/mpich/bin/cp2k.psmp
%{_libdir}/mpich/bin/dbm_miniapp.psmp
%{_libdir}/mpich/bin/dumpdcd.psmp
%{_libdir}/mpich/bin/graph.psmp
%{_libdir}/mpich/bin/grid_miniapp.psmp
%{_libdir}/mpich/bin/xyz2dcd.psmp
%{_libdir}/mpich/lib/libcp2k.so.*

%files mpich-devel
%{_fmoddir}/mpich/cp2k/
%{_libdir}/mpich/include/cp2k/
%{_libdir}/mpich/lib/cmake/cp2k/
%{_libdir}/mpich/lib/libcp2k.so
%{_libdir}/mpich/lib/pkgconfig/libcp2k.pc

%changelog
%autochangelog
