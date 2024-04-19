# No openmpi on i668 with openmpi 5 in Fedora 40+
%ifarch %{ix86}
%bcond_with openmpi
%else
%bcond_without openmpi
%endif

%ifarch x86_64
%bcond_without libxsmm
%else
# See https://bugzilla.redhat.com/show_bug.cgi?id=1515404
%bcond_with libxsmm
%endif

# Disable LTO due to https://bugzilla.redhat.com/show_bug.cgi?id=2243158
%global _lto_cflags %nil

Name:          cp2k
Version:       0.0.0
Release:       %autorelease
Summary:       Ab Initio Molecular Dynamics
License:       GPL-2.0-or-later
URL:           https://www.cp2k.org/
Source0:       https://github.com/cp2k/cp2k/releases/download/v%{version}/cp2k-%{version}.tar.bz2

# Remove testing packages that were previously packaged
# Can be removed at the end of F40
# Provides should not be necessary but might as well be thorough
Provides:      cp2k-testing = 2024.1-5
Obsoletes:     cp2k-testing < 2024.1-5
Provides:      cp2k-mpich-testing = 2024.1-5
Obsoletes:     cp2k-mpich-testing < 2024.1-5
Provides:      cp2k-openmpi-testing = 2024.1-5
Obsoletes:     cp2k-openmpi-testing < 2024.1-5

# Flaky MPI issues on s390x, and upstream do not officially support it yet
# https://github.com/cp2k/cp2k/issues/3362
ExcludeArch:   s390x

# Build dependencies
BuildRequires: cmake
BuildRequires: gcc-gfortran
BuildRequires: gcc-c++
BuildRequires: ninja-build
BuildRequires: python3-fypp
# Project dependencies
BuildRequires: flexiblas-devel
BuildRequires: dbcsr-devel >= 2.6.0
BuildRequires: fftw-devel
%if %{with libxsmm}
BuildRequires: libxsmm-devel >= 1.8.1-3
%endif
BuildRequires: libxc-devel >= 5.1.0
BuildRequires: spglib-devel
# Test dependencies
BuildRequires: python3
# For pathfix.py
BuildRequires: python3-devel

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

%if %{with openmpi}
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

%endif

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

# Fix test files
%{__python3} %{_rpmconfigdir}/redhat/pathfix.py -i "%{__python3} -Es" -p $(find . -type f -name *.py)

# $MPI_SUFFIX will be evaluated in the loops below, set by mpi modules
%global _vpath_builddir %{_vendor}-%{_target_os}-build${MPI_SUFFIX:-_serial}
# We are running the module load/unload manually until there is a macro-like way to expand this
. /etc/profile.d/modules.sh


%build
cmake_common_args=(
  "-G Ninja"
  "-DCP2K_DEBUG_MODE:BOOL=OFF"
  "-DCP2K_BLAS_VENDOR:STRING=FlexiBLAS"
  "-DCP2K_USE_STATIC_BLAS:BOOL=OFF"
  # Unit tests are included in REGTESTS
  # Note: Enabling this will write build files in the source folder :/
  "-DCP2K_ENABLE_REGTESTS:BOOL=ON"
  # Dependencies equivalent with Default
  "-DCP2K_USE_FFTW3:BOOL=ON"
  "-DCP2K_USE_COSMA:BOOL=OFF"  # Not packaged
  "-DCP2K_USE_LIBXSMM:BOOL=%{?with_libxsmm:ON}%{?without_libxsmm:OFF}"
  "-DCP2K_USE_LIBXC:BOOL=ON"
  "-DCP2K_USE_LIBINT2:BOOL=OFF"  # Detection is broken
  "-DCP2K_USE_SPGLIB:BOOL=ON"
)
for mpi in '' mpich %{?with_openmpi:openmpi}; do
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
for mpi in '' mpich %{?with_openmpi:openmpi}; do
  [ -n "$mpi" ] && module load mpi/${mpi}-%{_arch}
  %cmake_install
  [ -n "$mpi" ] && module unload mpi/${mpi}-%{_arch}
done

# TODO: Properly separate the installation of unit tests
rm -f %{_buildrootdir}/**/%{_bindir}/*_unittest.*
%if %{with openmpi}
rm -f %{_buildrootdir}/**/%{_libdir}/openmpi/bin/*_unittest.*
%endif
rm -f %{_buildrootdir}/**/%{_libdir}/mpich/bin/*_unittest.*


%check
export CP2K_DATA_DIR=%{buildroot}%{_datadir}/cp2k/data
# See %%_openmpi_load
export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe
test_common_args=(
  "--skip_regtests"
  "--maxtasks 4"  # Hard-coding maxtasks to avoid hanging. oversubscription should not matter
  "--ompthreads 2"
)
for mpi in '' mpich %{?with_openmpi:openmpi} ; do
  if [ -n "$mpi" ]; then
    # Another module load is done inside the do_regtest.sh. will use that instead
    module load mpi/${mpi}-%{_arch}
    bindir=${MPI_BIN}
    libdir=${MPI_LIB}
    export CP2K_STEM=%{buildroot}${MPI_BIN}/cp2k
    # Note, final position arguments are also here
    test_mpi_args=(
      "--mpiranks 2"
      "local_${mpi}"
      "psmp"
    )
  else
    bindir=%{_bindir}
    libdir=%{_libdir}
    export CP2K_STEM=%{buildroot}/usr/bin/cp2k
    test_mpi_args=(
      "local"
      "ssmp"
    )
  fi
  # Run packaged do_regtest.sh with appropriate buildroot runpaths
  env PATH=%{buildroot}${bindir}:${PATH} \
    LD_LIBRARY_PATH=%{buildroot}${libdir} \
    tests/do_regtest.py ${test_common_args[@]} ${test_mpi_args[@]}
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

%if %{with openmpi}
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
%endif

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
