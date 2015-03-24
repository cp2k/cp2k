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

cmake_ver=3.1.1
cmake_sha=b58694e545d51cde5756a894f53107e3d9e469360e1d92e7f6892b55ebc0bebf

parmetis_ver=4.0.2
parmetis_sha=5acbb700f457d3bda7d4bb944b559d7f21f075bb6fa4c33f42c261019ef2f0b2

scotch_ver=6.0.0
scotch_sha=e57e16c965bab68c1b03389005ecd8a03745ba20fd9c23081c0bb2336972d879

superlu_ver=3.3
superlu_sha=d2fd8dc847ae63ed7980cff2ad4db8d117640ecdf0234c9711e0f6ee1398cac2

pexsi_ver=0.7.3
pexsi_sha=077c0de466bdbef055992eea91e2fee5c10814552fb1f7e2a354a8881dc7d39f
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
export F90=gfortran
export CXX=g++
export CFLAGS="-O2 -g -Wno-error"
export FFLAGS="-O2 -g -Wno-error"
export FCFLAGS="-O2 -g -Wno-error"
export F90FLAGS="-O2 -g -Wno-error"
export F77FLAGS="-O2 -g -Wno-error"
export CXXFLAGS="-O2 -g -Wno-error"

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
  wget http://www.cp2k.org/static/downloads/valgrind-${valgrind_ver}.tar.bz2
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
  wget http://www.cp2k.org/static/downloads/lcov-${lcov_ver}.tar.gz
  echo "${lcov_sha} *lcov-${lcov_ver}.tar.gz" | sha256sum  --check
  tar -xzf lcov-${lcov_ver}.tar.gz
  cd lcov-${lcov_ver}
  # note.... this installs in ${INSTALLDIR}/usr/bin
  make PREFIX=${INSTALLDIR} install >& make.log
  cd ..
fi

echo "================== Installing CMake ================="
if [ -f cmake-${cmake_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/cmake-${cmake_ver}.tar.gz
  echo "${cmake_sha} *cmake-${cmake_ver}.tar.gz" | sha256sum  --check
  tar -xzf cmake-${cmake_ver}.tar.gz
  cd cmake-${cmake_ver}
  ./bootstrap --prefix=${INSTALLDIR} >& config.log
  make -j $nprocs >&  make.log
  make install >& install.log
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
# leaks related to PEXSI
leak:PPEXSIDFTDriver
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
LSAN_OPTIONS=suppressions=${INSTALLDIR}/lsan.supp ; export LSAN_OPTIONS
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
  # needed to install mpich ??
  unset F90; unset F90FLAGS
  wget http://www.cp2k.org/static/downloads/mpich-${mpich_ver}.tar.gz
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
  wget http://www.cp2k.org/static/downloads/scalapack_installer.tgz
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

echo "================== Installing ParMETIS =================="
if [ -f parmetis-${parmetis_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/parmetis-${parmetis_ver}.tar.gz
  echo "${parmetis_sha} *parmetis-${parmetis_ver}.tar.gz" | sha256sum  --check
  tar -xzf parmetis-${parmetis_ver}.tar.gz

  cd parmetis-${parmetis_ver}
  make config prefix=${INSTALLDIR} >& config.log
  make -j $nprocs >& make.log
  make install >& install.log

  # Have to build METIS again independently due to bug in ParMETIS make install
  cd metis
  make config prefix=${INSTALLDIR} >& config.log
  make -j $nprocs >& make.log
  make install >& install.log
  cd ../..
fi

echo "================== Installing PT-Scotch =================="
if [ -f scotch_${scotch_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget https://gforge.inria.fr/frs/download.php/31831/scotch_6.0.0.tar.gz
  echo "${scotch_sha} *scotch_${scotch_ver}.tar.gz" | sha256sum  --check
  tar -xzf scotch_${scotch_ver}.tar.gz
  cd scotch_${scotch_ver}/src
  cat Make.inc/Makefile.inc.x86-64_pc_linux2 | \
  sed "s|\(^CFLAGS\).*|\1 =  $CFLAGS -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Drestrict=__restrict -DIDXSIZE64|" > Makefile.inc

  make scotch -j $nprocs >& make.log
  make ptscotch -j $nrocs >& make.log
  make install prefix=${INSTALLDIR} >& install.log
  cd ../..
fi

echo "================== Installing SuperLU_DIST =================="

if [ -f superlu_dist_${superlu_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/superlu_dist_${superlu_ver}.tar.gz
  echo "${superlu_sha} *superlu_dist_${superlu_ver}.tar.gz" | sha256sum  --check
  tar -xzf superlu_dist_${superlu_ver}.tar.gz

  cd SuperLU_DIST_${superlu_ver}
  mv make.inc make.inc.orig
  cat <<EOF >> make.inc
PLAT=_x86_64
DSUPERLULIB= ${PWD}/lib/libsuperlu_dist_${superlu_ver}.a
LIBS=\$(DSUPERLULIB) -lparmetis -lmetis -lrefblas
ARCH=ar
ARCHFLAGS=cr
RANLIB=ranlib
CC=mpicc
CFLAGS=${CFLAGS}
NOOPTS=-O0
FORTRAN=mpif90
F90FLAGS=${FFLAGS}
LOADER=\$(CC)
LOADOPTS=${CFLAGS}
CDEFS=-DAdd_
EOF
  make &> make.log #-j $nprocs will crash
  # no make install
  chmod a+r lib/* SRC/*.h
  cp lib/libsuperlu_dist_${superlu_ver}.a ${INSTALLDIR}/lib/
  mkdir -p ${INSTALLDIR}/include/superlu_dist_${superlu_ver}
  cp SRC/*.h ${INSTALLDIR}/include/superlu_dist_${superlu_ver}/
  cd ..
fi

echo "================== Installing PEXSI =================="

if [ -f pexsi_v${pexsi_ver}.tar.gz ]; then
  echo "Installation already started, skipping it."
else
  wget http://www.cp2k.org/static/downloads/pexsi_v${pexsi_ver}.tar.gz
  echo "${pexsi_sha} *pexsi_v${pexsi_ver}.tar.gz" | sha256sum  --check
  mkdir pexsi_v${pexsi_ver}
  tar -xzf pexsi_v${pexsi_ver}.tar.gz -C pexsi_v${pexsi_ver}

  cd pexsi_v${pexsi_ver}

  # fix uninitialised variables
  cat <<EOF >> pexsi.patch 
--- src/get_perm_c_parmetis.c
+++ src/get_perm_c_parmetis.c
@@ -131,0 +132 @@
+ dist_order = 0;
--- src/pdsymbfact.c
+++ src/pdsymbfact.c
@@ -87,0 +88 @@
+ symb_mem_usage.total = 0.;
--- src/pzsymbfact.c
+++ src/pzsymbfact.c
@@ -87,0 +88 @@
+ symb_mem_usage.total = 0.;
EOF
  patch -p0 < pexsi.patch >& patch.log

  cat config/make.inc.linux.gnu | \
  sed 's|\(PAR_ND_LIBRARY *=\).*|\1 parmetis|' |\
  sed 's|\(SEQ_ND_LIBRARY *=\).*|\1 metis|' |\
  sed "s|\(PEXSI_DIR *=\).*|\1 ${PWD}|" |\
  sed "s|\(CPP_LIB *=\).*|\1 -lstdc++|" |\
  sed 's|\(LAPACK_LIB *=\).*|\1 -lreflapack|' |\
  sed 's|\(BLAS_LIB *=\).*|\1 -lrefblas|' |\
  sed 's|\(\bMETIS_LIB *=\).*|\1 -lmetis|' |\
  sed 's|\(PARMETIS_LIB *=\).*|\1 -lparmetis|' |\
  sed "s|\(DSUPERLU_LIB *=\).*|\1 -lsuperlu_dist_${superlu_ver}|" |\
  sed 's|#FLOADOPTS *=.*|FLOADOPTS    = ${LIBS} ${CPP_LIB}|' |\
  sed "s|\(DSUPERLU_INCLUDE *=\).*|\1 -I${INSTALLDIR}/include/superlu_dist_${superlu_ver}|" |\
  sed 's|\(INCLUDES *=\).*|\1 ${DSUPERLU_INCLUDE} ${PEXSI_INCLUDE}|' |\
  sed "s|\(COMPILE_FLAG *=\).*|\1 ${CFLAGS}|" |\
  sed "s|\(SUFFIX *=\).*|\1 linux_v${pexsi_ver}|" |\
  sed 's|\(DSUPERLU_DIR *=\).*|\1|' |\
  sed 's|\(METIS_DIR *=\).*|\1|' |\
  sed 's|\(PARMETIS_DIR *=\).*|\1|' |\
  sed 's|\(PTSCOTCH_DIR *=\).*|\1|' |\
  sed 's|\(LAPACK_DIR *=\).*|\1|' |\
  sed 's|\(BLAS_DIR *=\).*|\1|' |\
  sed 's|\(GFORTRAN_LIB *=\).*|\1|' > make.inc
  cd src
  make -j $nprocs >& make.log
  # no make install
  chmod a+r libpexsi_linux_v${pexsi_ver}.a
  cp libpexsi_linux_v${pexsi_ver}.a ${INSTALLDIR}/lib/

  # make fortran interface
  cd ../fortran
  make >& make.log #-j $nprocs will crash
  chmod a+r f_ppexsi_interface.mod
  cp f_ppexsi_interface.mod ${INSTALLDIR}/include/ 
  cd ..
  # no need to install PEXSI headers
  #mkdir -p ${INSTALLDIR}/include/pexsi_v${pexsi_ver}
  #cp include/* ${INSTALLDIR}/include/pexsi_v${pexsi_ver}/
  cd ..
fi

echo "==================== generating arch files ===================="
echo "arch files can be found in the ${INSTALLDIR}/arch subdirectory"
mkdir -p ${INSTALLDIR}/arch

#
# unfortunately, optimal flags depend on compiler etc.
#
WFLAGS="-Waliasing -Wampersand -Wc-binding-type -Wintrinsic-shadow -Wintrinsics-std -Wline-truncation -Wno-tabs -Wrealloc-lhs-all -Wtarget-lifetime -Wunderflow -Wunused-but-set-variable -Wunused-variable -Wconversion -Werror"
DEBFLAGS="-fcheck=bounds,do,recursion,pointer -fsanitize=leak"
BASEFLAGS="-std=f2003 -fimplicit-none -ffree-form -fno-omit-frame-pointer -g -O1"
PARAFLAGS="-D__parallel -D__SCALAPACK -D__LIBPEXSI -D__MPI_VERSION=3"
CUDAFLAGS="-D__ACC -D__DBCSR_ACC -D__PW_CUDA"
OPTFLAGS="-O3 -march=native -ffast-math \$(PROFOPT)"
DFLAGS="-D__LIBINT -D__FFTW3 -D__LIBXC2 -D__LIBINT_MAX_AM=6 -D__LIBDERIV_MAX_AM1=5"
CFLAGS="\$(DFLAGS) -I\$(CP2KINSTALLDIR)/include -fno-omit-frame-pointer -g -O1"

# Link to SCOTCH
LIB_PEXSI="-lpexsi_linux_v${pexsi_ver} -lsuperlu_dist_${superlu_ver} -lptscotchparmetis -lptscotch -lptscotcherr -lscotchmetis -lscotch -lscotcherr"

# Link to ParMETIS
#LIB_PEXSI="-lpexsi_linux_v${pexsi_ver} -lsuperlu_dist_${superlu_ver} -lparmetis -lmetis"

cat << EOF > ${INSTALLDIR}/arch/local.pdbg
CC       = gcc
CPP      =
FC       = mpif90 
LD       = mpif90
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS} ${PARAFLAGS}
FCFLAGS  = -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${DEBFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib/ \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint $LIB_PEXSI -lscalapack -lreflapack -lrefblas -lstdc++ -lfftw3
EOF

cat << EOF > ${INSTALLDIR}/arch/local.popt
CC       = gcc
CPP      =
FC       = mpif90 
LD       = mpif90
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS} ${PARAFLAGS}
FCFLAGS  = -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${OPTFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint $LIB_PEXSI -lscalapack -lreflapack -lrefblas -lstdc++ -lfftw3 
EOF

cat << EOF > ${INSTALLDIR}/arch/local.psmp
CC       = gcc
CPP      =
FC       = mpif90 
LD       = mpif90
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS} ${PARAFLAGS}
FCFLAGS  = -fopenmp -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${OPTFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib/ \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint $LIB_PEXSI -lscalapack -lreflapack -lrefblas -lstdc++ -lfftw3 -lfftw3_omp
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
LIBS     = -lxc -lderiv -lint -lreflapack -lrefblas -lstdc++ -lfftw3 -lfftw3_omp
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
DFLAGS   = ${DFLAGS} ${PARAFLAGS}
FCFLAGS  = -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} -O3 \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib \$(FCFLAGS)
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint $LIB_PEXSI -lscalapack -lreflapack -lrefblas -lstdc++ -lfftw3
EOF

cat << EOF > ${INSTALLDIR}/arch/local_cuda.psmp
NVCC     = nvcc -D__GNUC_MINOR__=6
CC       = gcc
CPP      =
FC       = mpif90 
LD       = mpif90
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS} ${CUDAFLAGS} ${PARAFLAGS}
FCFLAGS  = -fopenmp -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${OPTFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib/ -L/usr/local/cuda/lib64 \$(FCFLAGS)
NVFLAGS  = \$(DFLAGS) -g -O2 -arch sm_35
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint $LIB_PEXSI -lscalapack -lreflapack -lrefblas -lstdc++ -lfftw3 -lfftw3_omp -lcudart -lcufft -lcublas -lrt
EOF

cat << EOF > ${INSTALLDIR}/arch/local_cuda.ssmp
NVCC     = nvcc -D__GNUC_MINOR__=6
CC       = gcc
CPP      =
FC       = gfortran 
LD       = gfortran
AR       = ar -r
WFLAGS   = ${WFLAGS}
DFLAGS   = ${DFLAGS} ${CUDAFLAGS}
FCFLAGS  = -fopenmp -I\$(CP2KINSTALLDIR)/include ${BASEFLAGS} ${OPTFLAGS} \$(DFLAGS) \$(WFLAGS)
LDFLAGS  = -L\$(CP2KINSTALLDIR)/lib/ -L/usr/local/cuda/lib64 \$(FCFLAGS)
NVFLAGS  = \$(DFLAGS) -g -O2 -arch sm_35
CFLAGS   = ${CFLAGS}
LIBS     = -lxc -lderiv -lint -lreflapack -lrefblas -lstdc++ -lfftw3 -lfftw3_omp -lcudart -lcufft -lcublas -lrt
EOF

echo "========================== usage ========================="
echo "done!"
echo "now copy: cp ${INSTALLDIR}/arch/* to the cp2k/arch/ directory"
echo "to use this toolchain or the cp2k version compiled with it you will first need to execute at the prompt:"
echo "source ${SETUPFILE}"
echo "to build CP2K you should change directory cd cp2k/makefiles/"
echo "make -j${nprocs}" 'ARCH=local VERSION="sdbg sopt ssmp popt pdbg psmp"'

#EOF
