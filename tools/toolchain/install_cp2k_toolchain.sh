#!/bin/bash -e
#
# This script installs a fairly complete up to date toolchain for development and use of cp2k.
# It compiles tools etc such that they are suitable for debugging cp2k.
# Full installation / compilation can take a while.
#
# it does currently does not try to build an efficient blas,
# nor libraries such as libsmm, which might be important for performance
#
# updating version could be easy, just change here:
#

binutils_ver=2.24
binutils_sha=4930b2886309112c00a279483eaef2f0f8e1b1b62010e0239c16b22af7c346d4

valgrind_ver=3.10.0
valgrind_sha=03047f82dfc6985a4c7d9d2700e17bc05f5e1a0ca6ad902e5d6c81aeb720edc9

lcov_ver=1.11
lcov_sha=c282de8d678ecbfda32ce4b5c85fc02f77c2a39a062f068bd8e774d29ddc9bf8

gcc_ver=4.9.2
gcc_sha=3e573826ec8b0d62d47821408fbc58721cd020df3e594cd492508de487a43b5e

mpich_ver=3.1.2
mpich_sha=37c3ba2d3cd3f4ea239497d9d34bd57a663a34e2ea25099c2cbef118c9156587

# no version numbering used.
scalapack_ver=XXXXXX
scalapack_sha=8ecae0896a63c7980a71c22520afed6cfefb52b17c2496223026cc01caf07855

libxc_ver=2.0.1
libxc_sha=c332f08648ec2bc7ccce83e45a84776215aa5dfebc64fae2a23f2ac546d41ea4

libint_ver=1.1.4
libint_sha=f67b13bdf1135ecc93b4cff961c1ff33614d9f8409726ddc8451803776885cff

fftw_ver=3.3.4
fftw_sha=8f0cde90929bc05587c3368d2f15cd0530a60b8a9912a8e2979a72dbe5af0982

# will need hand-checking for non-standard dir name in tar
elpa_ver=2013.11.008
elpa_sha=d4a028fddb64a7c1454f08b930525cce0207893c6c770cb7bf92ab6f5d44bd78

#
#
#

rootdir=${PWD}

mkdir -p build
cd build

INSTALLDIR=${rootdir}/install
echo "All tools will be installed in " ${INSTALLDIR}
mkdir -p ${INSTALLDIR}


#
# number of processes to use for compilation
#
nprocs=`nproc --all`

#
# first get an up-to-date toolchain.
#

export CC=gcc
export FC=gfortran
export F77=gfortran
export CXX=g++
export CFLAGS="-O2 -g -Wno-error"
export FFLAGS="-O2 -g -Wno-error"
export FCFLAGS="-O2 -g -Wno-error"
export CXXFLAGS="-O2 -g -Wno-error"

# F90 and F90FLAGS confuse mpich's configure.
unset F90
unset F90FLAGS

echo "==================== Installing binutils ================="
if [ -f binutils-${binutils_ver}.tar.gz  ]; then
  echo "Installation already started, skipping it."
else
  wget http://ftp.gnu.org/gnu/binutils/binutils-${binutils_ver}.tar.gz
  echo "${binutils_sha} *binutils-${binutils_ver}.tar.gz" | sha256sum  --check
  tar -xzf binutils-${binutils_ver}.tar.gz
  cd binutils-${binutils_ver}
  ./configure --prefix=${INSTALLDIR} --enable-gold --enable-plugins >& config.log
  make -j $nprocs >& make.log
  make -j $nprocs install >& install.log
  cd ..
fi

echo "==================== Installing valgrind ================="
if [ -f valgrind-${valgrind_ver}.tar.bz2 ]; then
  echo "Installation already started, skipping it."
else
  wget http://valgrind.org/downloads/valgrind-${valgrind_ver}.tar.bz2
  echo "${valgrind_sha} *valgrind-${valgrind_ver}.tar.bz2" | sha256sum  --check
  tar -xjf valgrind-${valgrind_ver}.tar.bz2
  cd valgrind-${valgrind_ver}
  ./configure --prefix=${INSTALLDIR} >& config.log
  make -j $nprocs >& make.log
  make -j $nprocs install >& install.log
  cd ..
fi

echo "==================== Installing lcov ====================="
if [ -f lcov-${lcov_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://downloads.sourceforge.net/ltp/lcov-${lcov_ver}.tar.gz
  echo "${lcov_sha} *lcov-${lcov_ver}.tar.gz" | sha256sum  --check
  tar -xzf lcov-${lcov_ver}.tar.gz
  cd lcov-${lcov_ver}
  # note.... this installs in ${INSTALLDIR}/usr/bin
  make PREFIX=${INSTALLDIR} install >& make.log
  cd ..
fi

echo "==================== Installing gcc ======================"
if [ -f gcc-${gcc_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget https://ftp.gnu.org/gnu/gcc/gcc-${gcc_ver}/gcc-${gcc_ver}.tar.gz
  echo "${gcc_sha} *gcc-${gcc_ver}.tar.gz" | sha256sum  --check
  tar -xzf gcc-${gcc_ver}.tar.gz
  cd gcc-${gcc_ver}
  ./contrib/download_prerequisites >& prereq.log
  GCCROOT=${PWD}
  mkdir obj
  cd obj
  ${GCCROOT}/configure --prefix=${INSTALLDIR}  --enable-languages=c,c++,fortran --disable-multilib --disable-bootstrap --enable-lto --enable-plugins >& config.log
  make -j $nprocs >& make.log
  make -j $nprocs install >& install.log
  cd ../..
fi

# lsan suppressions for known leaks are created as well, this might need to be adjusted for the versions of the software employed
cat << EOF > ${INSTALLDIR}/lsan.supp
# known leak either related to mpi or scalapack  (e.g. showing randomly for Fist/regtest-7-2/UO2-2x2x2-genpot_units.inp)
leak:__cp_fm_types_MOD_cp_fm_write_unformatted
EOF

# now we need these tools and compiler to be in the path
cat << EOF > ${INSTALLDIR}/setup
if [ -z "\${LD_LIBRARY_PATH}" ]
then
    LD_LIBRARY_PATH=${INSTALLDIR}/lib64:${INSTALLDIR}/lib; export LD_LIBRARY_PATH
else
    LD_LIBRARY_PATH=${INSTALLDIR}/lib64:${INSTALLDIR}/lib:\${LD_LIBRARY_PATH}; export LD_LIBRARY_PATH
fi
if [ -z "\${PATH}" ]
then
    PATH=${INSTALLDIR}/bin:${INSTALLDIR}/usr/bin; export PATH
else
    PATH=${INSTALLDIR}/bin:${INSTALLDIR}/usr/bin:\$PATH; export PATH
fi 
CP2KINSTALLDIR=${INSTALLDIR} ; export CP2KINSTALLDIR
LSAN_OPTIONS=suppressions=$CP2KINSTALLDIR/lsan.supp ; export LSAN_OPTIONS
EOF
SETUPFILE=${INSTALLDIR}/setup
source ${SETUPFILE}

# set some flags, leading to nice stack traces on crashes, yet, are optimized
export CFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math"
export FFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math"
export FCFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math"
export CXXFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math"

echo "=================== Installing mpich ====================="
if [ -f mpich-${mpich_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.mpich.org/static/downloads/${mpich_ver}/mpich-${mpich_ver}.tar.gz
  echo "${mpich_sha} *mpich-${mpich_ver}.tar.gz" | sha256sum  --check
  tar -xzf mpich-${mpich_ver}.tar.gz
  cd mpich-${mpich_ver}
  ./configure --prefix=${INSTALLDIR} >& config.log
  make -j $nprocs >& make.log
  make -j $nprocs install >& install.log
  cd ..
fi

echo "================= Installing scalapack ==================="
if [ -f scalapack_installer.tgz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.netlib.org/scalapack/scalapack_installer.tgz
  echo "${scalapack_sha} *scalapack_installer.tgz" | sha256sum  --check
  tar -xzf scalapack_installer.tgz
  # we dont know the version
  cd scalapack_installer_*
  SLROOT=${PWD}
  # needs fixing for compile options, we use echo as mpirun command to avoid testing,
  # yet download blas / lapack... whole installer is a bit serial as well (and fails with --make="make -j32"
  ./setup.py --mpirun=echo --downblas --downlapack >& make.log
  # copy libraries where we like them
  cp install/lib/* ${INSTALLDIR}/lib/
  cd ..
fi

echo "==================== Installing libxc ===================="
if [ -f libxc-${libxc_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/libxc-${libxc_ver}.tar.gz
  echo "${libxc_sha} *libxc-${libxc_ver}.tar.gz" | sha256sum  --check
  tar -xzf libxc-${libxc_ver}.tar.gz
  cd libxc-${libxc_ver}
  ./configure  --prefix=${INSTALLDIR} >& config.log
  make -j $nprocs >& make.log
  make install >& install.log
  cd ..
fi

echo "=================== Installing libint ===================="
if [ -f libint-${libint_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/libint-${libint_ver}.tar.gz
  echo "${libint_sha} *libint-${libint_ver}.tar.gz" | sha256sum  --check
  tar -xzf libint-${libint_ver}.tar.gz
  cd libint-${libint_ver}
  ./configure  --prefix=${INSTALLDIR} --with-libint-max-am=5 --with-libderiv-max-am1=4 --with-cc-optflags="$CFLAGS" --with-cxx-optflags="$CXXFLAGS" >& config.log
  make -j $nprocs >&  make.log
  make install >& install.log
  cd ..
fi

echo "==================== Installing FFTW ====================="
if [ -f fftw-${fftw_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/fftw-${fftw_ver}.tar.gz
  echo "${fftw_sha} *fftw-${fftw_ver}.tar.gz" | sha256sum  --check
  tar -xzf fftw-${fftw_ver}.tar.gz
  cd fftw-${fftw_ver}
  ./configure  --prefix=${INSTALLDIR} --enable-openmp >& config.log
  make -j $nprocs >&  make.log
  make install >& install.log
  cd ..
fi

echo "==================== Installing ELPA ====================="
# elpa expect FC to be an mpi fortran compiler that's happy with long lines, and that a bunch of libs can be found
export FC="mpif90 -ffree-line-length-none"
export LDFLAGS="-L${INSTALLDIR}/lib"
export LIBS="-lscalapack -lreflapack -lrefblas"
if [ -f elpa-${elpa_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/elpa-${elpa_ver}.tar.gz
  echo "${elpa_sha} *elpa-${elpa_ver}.tar.gz" | sha256sum  --check
  tar -xzf elpa-${elpa_ver}.tar.gz

  # need both flavors ?
  cp -rp ELPA_2013.11 ELPA_2013.11-mt

  cd ELPA_2013.11-mt
  ./configure  --prefix=${INSTALLDIR} --enable-openmp=yes --with-generic --enable-shared=no >& config.log
  # wrong deps, build serially ?
  make -j 1 >&  make.log
  make install >& install.log
  cd ..

  cd ELPA_2013.11
  ./configure  --prefix=${INSTALLDIR} --enable-openmp=no --with-generic --enable-shared=no >& config.log
  # wrong deps, build serially ?
  make -j 1 >&  make.log
  make install >& install.log
  cd ..
fi

echo "================== generating arch files ================="
echo "arch files can be found in the ${INSTALLDIR}/arch subdirectory"
mkdir -p ${INSTALLDIR}/arch


#
# unfortunately, optimal flags depend on compiler etc.
#
WFLAGS="-Waliasing -Wampersand -Wc-binding-type -Wintrinsic-shadow -Wintrinsics-std -Wline-truncation -Wno-tabs -Wrealloc-lhs-all -Wtarget-lifetime -Wunderflow -Wunused-but-set-variable -Wunused-variable -Wconversion -Werror"
DEBFLAGS="-fcheck=bounds,do,recursion,pointer -fsanitize=leak"
BASEFLAGS="-std=f2003 -fimplicit-none -ffree-form -fno-omit-frame-pointer -g -O1"
OPTFLAGS="-O3 -march=native -ffast-math \$(PROFOPT)"
DFLAGS="-D__LIBINT -D__FFTW3 -D__LIBXC2 -D__LIBINT_MAX_AM=6 -D__LIBDERIV_MAX_AM1=5"
CFLAGS="\$(DFLAGS) -I\$(CP2KINSTALLDIR)/include -fno-omit-frame-pointer -g -O1"


cat << EOF > ${INSTALLDIR}/arch/local.pdbg
CC       = gcc
CPP      =
FC       = mpif90 
LD       = mpif90
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS} -D__parallel -D__SCALAPACK
FCFLAGS  = -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${DEBFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib/ \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint  -lscalapack -lreflapack -lrefblas -lstdc++ -lfftw3
EOF

cat << EOF > ${INSTALLDIR}/arch/local.popt
CC       = gcc
CPP      =
FC       = mpif90 
LD       = mpif90
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS} -D__parallel -D__SCALAPACK
FCFLAGS  = -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${OPTFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint  -lscalapack -lreflapack -lrefblas -lstdc++ -lfftw3 
EOF

cat << EOF > ${INSTALLDIR}/arch/local.psmp
CC       = gcc
CPP      =
FC       = mpif90 
LD       = mpif90
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS} -D__parallel -D__SCALAPACK
FCFLAGS  = -fopenmp -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${OPTFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib/ \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint  -lscalapack -lreflapack -lrefblas -lstdc++ -lfftw3 -lfftw3_omp
EOF

cat << EOF > ${INSTALLDIR}/arch/local.sdbg
CC       = gcc
CPP      =
FC       = gfortran 
LD       = gfortran
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS}
FCFLAGS  = -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${DEBFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint -lreflapack -lrefblas -lstdc++ -lfftw3 
EOF

cat << EOF > ${INSTALLDIR}/arch/local.sopt
CC       = gcc
CPP      =
FC       = gfortran 
LD       = gfortran
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS}
FCFLAGS  = -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${OPTFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib/ \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint -lreflapack -lrefblas -lstdc++ -lfftw3 
EOF

cat << EOF > ${INSTALLDIR}/arch/local.ssmp
CC       = gcc
CPP      =
FC       = gfortran 
LD       = gfortran
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS}
FCFLAGS  = -fopenmp -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${OPTFLAGS}  \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -fopenmp -L\$(CP2KINSTALLDIR)/lib/ \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint  -lreflapack -lrefblas -lstdc++ -lfftw3 -lfftw3_omp
EOF

cat << EOF > ${INSTALLDIR}/arch/local_valgrind.sdbg
CC       = gcc
CPP      =
FC       = gfortran 
LD       = gfortran
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS}
FCFLAGS  = -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} -O3  \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint -lreflapack -lrefblas -lstdc++ -lfftw3 
EOF

cat << EOF > ${INSTALLDIR}/arch/local_valgrind.pdbg
CC       = gcc
CPP      =
FC       = mpif90 
LD       = mpif90
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS}  -D__parallel -D__SCALAPACK
FCFLAGS  = -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} -O3 \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint  -lscalapack -lreflapack -lrefblas -lstdc++ -lfftw3
EOF

cat << EOF > ${INSTALLDIR}/arch/local_cuda.psmp
NVCC     = nvcc -D__GNUC_MINOR__=6
CC       = gcc
CPP      =
FC       = mpif90 
LD       = mpif90
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS} -D__ACC -D__DBCSR_ACC -D__PW_CUDA -D__parallel -D__SCALAPACK
FCFLAGS  = -fopenmp -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${OPTFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib/ -L/usr/local/cuda/lib64 \$(FCFLAGS)
NVFLAGS  = \$(DFLAGS) -g -O2 -arch sm_35
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint  -lscalapack -lreflapack -lrefblas -lstdc++ -lfftw3 -lfftw3_omp -lcudart -lcufft -lcublas -lrt
EOF

cat << EOF > ${INSTALLDIR}/arch/local_cuda.ssmp
NVCC     = nvcc -D__GNUC_MINOR__=6
CC       = gcc
CPP      =
FC       = gfortran 
LD       = gfortran
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS} -D__ACC -D__DBCSR_ACC -D__PW_CUDA
FCFLAGS  = -fopenmp -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${OPTFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib/ -L/usr/local/cuda/lib64 \$(FCFLAGS)
NVFLAGS  = \$(DFLAGS) -g -O2 -arch sm_35
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint  -lreflapack -lrefblas -lstdc++ -lfftw3 -lfftw3_omp -lcudart -lcufft -lcublas -lrt
EOF

echo "========================== usage ========================="
echo "done!"
echo "now copy: cp ${INSTALLDIR}/arch/* to the cp2k/arch/ directory"
echo "to use this toolchain or the cp2k version compiled with it you will first need to execute at the prompt:"
echo "source ${SETUPFILE}"
echo "to build CP2K you should change directory cd cp2k/makefiles/"
echo "make -j${nprocs}" 'ARCH=local  VERSION="sdbg sopt ssmp popt pdbg psmp"'

#EOF
