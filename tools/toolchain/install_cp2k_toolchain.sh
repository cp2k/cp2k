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

binutilsver=2.24

valgrindver=3.10.0

lcovver=1.11

gccver=4.9.2

mpichver=3.1.2

# no version numbering used.
scalapackver=XXXXXX

libxcver=2.0.1
libxcsha=c332f08648ec2bc7ccce83e45a84776215aa5dfebc64fae2a23f2ac546d41ea4

libintver=1.1.4
libintsha=f67b13bdf1135ecc93b4cff961c1ff33614d9f8409726ddc8451803776885cff

fftwver=3.3.4
fftwsha=8f0cde90929bc05587c3368d2f15cd0530a60b8a9912a8e2979a72dbe5af0982

# will need hand-checking for non-standard dir name in tar
elpaver=2013.11.008
elpasha=d4a028fddb64a7c1454f08b930525cce0207893c6c770cb7bf92ab6f5d44bd78

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

echo "==================== Installing binutils ================="
if [ -f binutils-${binutilsver}.tar.gz  ]; then
  echo "Installation already started, skipping it."
else
  wget http://ftp.gnu.org/gnu/binutils/binutils-${binutilsver}.tar.gz
  tar -xzf binutils-${binutilsver}.tar.gz
  cd binutils-${binutilsver}
  ./configure --prefix=${INSTALLDIR} --enable-gold --enable-plugins >& config.log
  make -j $nprocs >& make.log
  make -j $nprocs install >& install.log
  cd ..
fi

echo "==================== Installing valgrind ================="
if [ -f valgrind-${valgrindver}.tar.bz2 ]; then
  echo "Installation already started, skipping it."
else
  wget http://valgrind.org/downloads/valgrind-${valgrindver}.tar.bz2
  tar -xjf valgrind-${valgrindver}.tar.bz2
  cd valgrind-${valgrindver}
  ./configure --prefix=${INSTALLDIR} >& config.log
  make -j $nprocs >& make.log
  make -j $nprocs install >& install.log
  cd ..
fi

echo "==================== Installing lcov ======================"
if [ -f lcov-${lcovver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://downloads.sourceforge.net/ltp/lcov-${lcovver}.tar.gz
  tar -xzf lcov-${lcovver}.tar.gz
  cd lcov-${lcovver}
  # note.... this installs in ${INSTALLDIR}/usr/bin
  make PREFIX=${INSTALLDIR} install >& make.log
  cd ..
fi

echo "==================== Installing gcc ======================"
if [ -f gcc-${gccver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget https://ftp.gnu.org/gnu/gcc/gcc-${gccver}/gcc-${gccver}.tar.gz
  tar -xzf gcc-${gccver}.tar.gz
  cd gcc-${gccver}
  ./contrib/download_prerequisites >& prereq.log
  GCCROOT=${PWD}
  mkdir obj
  cd obj
  ${GCCROOT}/configure --prefix=${INSTALLDIR}  --enable-languages=c,c++,fortran --disable-multilib --disable-bootstrap --enable-lto --enable-plugins >& config.log
  make -j $nprocs >& make.log
  make -j $nprocs install >& install.log
  cd ../..
fi

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
EOF
SETUPFILE=${INSTALLDIR}/setup
source ${SETUPFILE}

# set some flags, leading to nice stack traces on crashes, yet, are optimized
export CFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math"
export FFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math"
export FCFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math"
export CXXFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math"

echo "==================== Installing mpich ======================"
if [ -f mpich-${mpichver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.mpich.org/static/downloads/${mpichver}/mpich-${mpichver}.tar.gz
  tar -xzf mpich-${mpichver}.tar.gz
  cd mpich-${mpichver}
  ./configure --prefix=${INSTALLDIR} >& config.log
  make -j $nprocs >& make.log
  make -j $nprocs install >& install.log
  cd ..
fi

echo "==================== Installing scalapack ======================"
if [ -f scalapack_installer.tgz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.netlib.org/scalapack/scalapack_installer.tgz
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
if [ -f libxc-${libxcver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/libxc-${libxcver}.tar.gz
  echo "${libxcsha} *libxc-${libxcver}.tar.gz" | sha256sum  --check
  tar -xzf libxc-${libxcver}.tar.gz
  cd libxc-${libxcver}
  ./configure  --prefix=${INSTALLDIR} >& config.log
  make -j $nprocs >& make.log
  make install >& install.log
  cd ..
fi

echo "==================== Installing libint ===================="
if [ -f libint-${libintver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/libint-${libintver}.tar.gz
  echo "${libintsha} *libint-${libintver}.tar.gz" | sha256sum  --check
  tar -xzf libint-${libintver}.tar.gz
  cd libint-${libintver}
  ./configure  --prefix=${INSTALLDIR} --with-libint-max-am=5 --with-libderiv-max-am1=4 --with-cc-optflags="$CFLAGS" --with-cxx-optflags="$CXXFLAGS" >& config.log
  make -j $nprocs >&  make.log
  make install >& install.log
  cd ..
fi

echo "==================== Installing FFTW ===================="
if [ -f fftw-${fftwver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/fftw-${fftwver}.tar.gz
  echo "${fftwsha} *fftw-${fftwver}.tar.gz" | sha256sum  --check 
  tar -xzf fftw-${fftwver}.tar.gz
  cd fftw-${fftwver}
  ./configure  --prefix=${INSTALLDIR} --enable-openmp >& config.log
  make -j $nprocs >&  make.log
  make install >& install.log
  cd ..
fi

echo "==================== Installing ELPA ===================="
# elpa expect FC to be an mpi fortran compiler that's happy with long lines, and that a bunch of libs can be found
export FC="mpif90 -ffree-line-length-none"
export LDFLAGS="-L${INSTALLDIR}/lib"
export LIBS="-lscalapack -lreflapack -lrefblas"
if [ -f elpa-${elpaver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/elpa-${elpaver}.tar.gz
  echo "${elpasha} *elpa-${elpaver}.tar.gz" | sha256sum  --check
  tar -xzf elpa-${elpaver}.tar.gz

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

echo "==================== generating arch files ===================="
echo "arch files can be found in the ${INSTALLDIR}/arch subdirectory"
mkdir -p ${INSTALLDIR}/arch


#
# unfortunately, optimal flags depend on compiler etc.
#
WFLAGS="-Waliasing -Wampersand -Wc-binding-type -Wintrinsic-shadow -Wintrinsics-std -Wline-truncation -Wno-tabs -Wrealloc-lhs-all -Wtarget-lifetime -Wunderflow -Wunused-but-set-variable -Wunused-variable -Werror"
DEBFLAGS="-fcheck=bounds,do,recursion,pointer -fsanitize=leak"
BASEFLAGS="-std=f2003 -fimplicit-none -ffree-form -fno-omit-frame-pointer -g -O1"
OPTFLAGS="-O3 -march=native -ffast-math"
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

echo "==================== usage ===================="
echo "done!"
echo "now copy: cp ${INSTALLDIR}/arch/* to the cp2k/arch/ directory"
echo "to use this toolchain or the cp2k version compiled with it you will first need to execute at the prompt:"
echo "source ${SETUPFILE}"
echo "to build CP2K you should change directory cd cp2k/makefiles/"
echo "make -j${nprocs}" 'ARCH=local VERSION="sdbg sopt ssmp popt pdbg psmp"'

#EOF
