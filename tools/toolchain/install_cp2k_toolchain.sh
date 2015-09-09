#!/bin/bash -e
#
# This script installs a fairly complete up to date toolchain for development and use of cp2k.
# It compiles tools etc such that they are suitable for debugging cp2k.
# Full installation / compilation can take a while.
#
# trap errors

error_handler() {
   echo "Non-zero exit code in this script on line $1"
   echo "Aborting, toolchain incomplete"
   exit 1
}
trap 'error_handler ${LINENO}' ERR


# default settings
enable_tsan=false
mpichoice=mpich
mpich_ver=3.1.2
openmpi_ver=1.8.6
openblas_ver=v0.2.14-0-gd0c51c4
scalapack_ver=XXXXXX # no version numbering used.
libxc_ver=2.2.2
libint_ver=1.1.4
fftw_ver=3.3.4
elpa_ver=2015.05.001
cmake_ver=3.1.1
parmetis_ver=4.0.2
scotch_ver=6.0.0
superlu_ver=3.3
pexsi_ver=0.8.0
plumed_ver=2.2b
quip_ver=336cab5c03
#binutils_ver=2.24
binutils_ver=2.25
#valgrind_ver=3.10.0
valgrind_ver=3.10.1
lcov_ver=1.11
#gcc_ver=4.9.2
#gcc_ver=5.1.0
gcc_ver=5.2.0

# parse options
while [ $# -ge 1 ]; do
   case $1 in
   --mpich)
      mpichoice=mpich;;
   --openmpi)
      mpichoice=openmpi;;
   --enable-tsan)
      enable_tsan=true;;
   --enable-gcc-trunk)
      gcc_ver=master;;
   -help|-h|--help)
      echo "Usage: install_cp2k_toolchain.sh [OPTIONS]"
      echo "Installs a well defined development environment for CP2K"
      echo ""
      echo "Options:"
      echo "  -h, -help, --help         print this help screen."
      echo "  --mpich                   use the MPICH package. Default: on."
      echo "  --openmpi                 use the OpenMPI package. Default: off."
      echo "  --enable-tsan             compile entire toolchain with -fsanitize=thread. Default: off."
      echo "                            This is not for normal (production) use, but suitable for"
      echo "                            finding/testing/debugging threading issues during development."
      echo "  --enable-gcc-trunk        use a non-released, development version of gcc for testing."
      echo ""
      echo "For more information visit: <http://cp2k.org>"
      exit 0;;
   -*)
      echo "ERROR: Invalid command line flag $1 found"
      exit 1;;
   # Default case
   *)
      echo "ERROR: Unknown command line string $1 found"
      exit 1;;
   esac
   shift
done


# helper routine to check integrity of downloaded files
checksum() {
   filename=$1
   if grep $filename ../checksums.sha256 | sha256sum --quiet --check ; then
      echo "Checksum of ${filename} Ok"
   else
      echo "Checksum of ${filename} could not be verified, abort."
      rm -v ${filename}
      exit 255
   fi
}

# preliminaries
rootdir=${PWD}
mkdir -p build
cd build
INSTALLDIR=${rootdir}/install
echo "All tools will be installed in " ${INSTALLDIR}
mkdir -p ${INSTALLDIR}
nprocs=`nproc --all` # number of processes to use for compilation

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
   checksum binutils-${binutils_ver}.tar.gz
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
   checksum valgrind-${valgrind_ver}.tar.bz2
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
   checksum lcov-${lcov_ver}.tar.gz
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
   checksum cmake-${cmake_ver}.tar.gz
   tar -xzf cmake-${cmake_ver}.tar.gz
   cd cmake-${cmake_ver}
   ./bootstrap --prefix=${INSTALLDIR} >& config.log
   make -j $nprocs >&  make.log
   make install >& install.log
   cd ..
fi

echo "==================== Installing gcc ======================"
if [ -f gcc-${gcc_ver}.tar.gz -o -f gcc-${gcc_ver}.zip ]; then
   echo "Installation already started, skipping it."
else
   if [ "${gcc_ver}" == "master" ]; then
      # no check since this follows the gcc trunk svn repo and changes constantly
      wget -O gcc-master.zip https://github.com/gcc-mirror/gcc/archive/master.zip
      unzip -q gcc-master.zip 
   else
      wget https://ftp.gnu.org/gnu/gcc/gcc-${gcc_ver}/gcc-${gcc_ver}.tar.gz
      checksum gcc-${gcc_ver}.tar.gz
      tar -xzf gcc-${gcc_ver}.tar.gz
   fi
   cd gcc-${gcc_ver}
   ./contrib/download_prerequisites >& prereq.log
   GCCROOT=${PWD}
   mkdir obj
   cd obj
   ${GCCROOT}/configure --prefix=${INSTALLDIR}  --enable-languages=c,c++,fortran --disable-multilib --disable-bootstrap --enable-lto --enable-plugins >& config.log
   make -j $nprocs >& make.log
   make -j $nprocs install >& install.log

   if $enable_tsan ; then
      # now the tricky bit... we need to recompile in particular libgomp with -fsanitize=thread.. there is not configure option for this (as far as I know).
      # we need to go in the build tree and recompile / reinstall with proper options...
      # this is likely to break for later version of gcc, tested with 5.1.0
      # based on https://gcc.gnu.org/bugzilla/show_bug.cgi?id=55374#c10
      cd x86_64*/libgfortran
      make clean >& clean.log
      make -j $nprocs CFLAGS="-std=gnu99 -g -O2 -fsanitize=thread "  FCFLAGS="-g -O2 -fsanitize=thread" CXXFLAGS="-std=gnu99 -g -O2 -fsanitize=thread " LDFLAGS="-B`pwd`/../libsanitizer/tsan/.libs/ -Wl,-rpath,`pwd`/../libsanitizer/tsan/.libs/ -fsanitize=thread" >& make.log
      make install >& install.log
      cd ../libgomp
      make clean >& clean.log
      make -j $nprocs CFLAGS="-std=gnu99 -g -O2 -fsanitize=thread "  FCFLAGS="-g -O2 -fsanitize=thread" CXXFLAGS="-std=gnu99 -g -O2 -fsanitize=thread " LDFLAGS="-B`pwd`/../libsanitizer/tsan/.libs/ -Wl,-rpath,`pwd`/../libsanitizer/tsan/.libs/ -fsanitize=thread" >& make.log
      make install >& install.log
      cd $GCCROOT/obj/
   fi
   cd ../..
fi
if $enable_tsan ; then
   TSANFLAGS="-fsanitize=thread"
else
   TSANFLAGS=""
fi

# lsan & tsan suppressions for known leaks are created as well, this might need to be adjusted for the versions of the software employed
cat << EOF > ${INSTALLDIR}/lsan.supp
# known leak either related to mpi or scalapack  (e.g. showing randomly for Fist/regtest-7-2/UO2-2x2x2-genpot_units.inp)
leak:__cp_fm_types_MOD_cp_fm_write_unformatted
# leaks related to PEXSI
leak:PPEXSIDFTDriver
# tsan bugs likely related to gcc
# PR66756
deadlock:_gfortran_st_open
mutex:_gfortran_st_open
# PR66761
race:do_spin
race:gomp_team_end
#PR67303
race:gomp_iter_guided_next
# bugs related to removing/filtering blocks in DBCSR.. to be fixed
race:__dbcsr_block_access_MOD_dbcsr_remove_block
race:__dbcsr_operations_MOD_dbcsr_filter_anytype
race:__dbcsr_transformations_MOD_dbcsr_make_untransposed_blocks
EOF

# now we need these tools and compiler to be in the path
cat << EOF > ${INSTALLDIR}/setup
if [ -z "\${LD_LIBRARY_PATH}" ]
then
    export LD_LIBRARY_PATH=${INSTALLDIR}/lib64:${INSTALLDIR}/lib
else
    export LD_LIBRARY_PATH=${INSTALLDIR}/lib64:${INSTALLDIR}/lib:\${LD_LIBRARY_PATH}
fi
if [ -z "\${PATH}" ]
then
    export PATH=${INSTALLDIR}/bin:${INSTALLDIR}/usr/bin
else
    export PATH=${INSTALLDIR}/bin:${INSTALLDIR}/usr/bin:\$PATH
fi
export CP2KINSTALLDIR=${INSTALLDIR}
export LSAN_OPTIONS=suppressions=${INSTALLDIR}/lsan.supp
export TSAN_OPTIONS=suppressions=${INSTALLDIR}/lsan.supp
export CC=gcc
export CXX=g++
export FC=gfortran
export F77=gfortran
export F90=gfortran
EOF
SETUPFILE=${INSTALLDIR}/setup
source ${SETUPFILE}

# set some flags, leading to nice stack traces on crashes, yet, are optimized
export CFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"
export FFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"
export F77FLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"
export F90FLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"
export FCFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"
export CXXFLAGS="-O2 -ftree-vectorize -g -fno-omit-frame-pointer -march=native -ffast-math $TSANFLAGS"
export LDFLAGS=" $TSANFLAGS"

if [ "$mpichoice" == "openmpi" ]; then
   echo "=================== Installing openmpi ====================="
   if [ -f openmpi-${openmpi_ver}.tar.gz ]; then
      echo "Installation already started, skipping it."
   else
      wget http://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-${openmpi_ver}.tar.gz
      checksum openmpi-${openmpi_ver}.tar.gz
      tar -xzf openmpi-${openmpi_ver}.tar.gz
      cd openmpi-${openmpi_ver}
      ./configure --prefix=${INSTALLDIR} >& config.log
      make -j $nprocs >& make.log
      make -j $nprocs install >& install.log
      cd ..
   fi
   DFLAGS="${FLAGS} IF_MPI(-D__parallel -D__MPI_VERSION=3,)"
   #extra libs needed to link with mpif90 also applications based on C++
   LIBS="IF_MPI(-lmpi_cxx,) ${LIBS}"
fi


if [ "$mpichoice" == "mpich" ]; then
   echo "=================== Installing mpich ====================="
   if [ -f mpich-${mpich_ver}.tar.gz ]; then
      echo "Installation already started, skipping it."
   else
      # needed to install mpich ??
      unset F90; unset F90FLAGS
      wget http://www.cp2k.org/static/downloads/mpich-${mpich_ver}.tar.gz
      checksum mpich-${mpich_ver}.tar.gz
      tar -xzf mpich-${mpich_ver}.tar.gz
      cd mpich-${mpich_ver}
      ./configure --prefix=${INSTALLDIR} >& config.log
      make -j $nprocs >& make.log
      make -j $nprocs install >& install.log
      cd ..
   fi
   DFLAGS="${DFLAGS} IF_MPI(-D__parallel -D__MPI_VERSION=3,)"
fi


echo "================= Installing openblas ==================="
if [ -f xianyi-OpenBLAS-${openblas_ver}.zip ]; then
   echo "Installation already started, skipping it."
else
   wget http://www.cp2k.org/static/downloads/xianyi-OpenBLAS-${openblas_ver}.zip
   checksum xianyi-OpenBLAS-${openblas_ver}.zip
   unzip xianyi-OpenBLAS-${openblas_ver}.zip >& unzip.log
   cd xianyi-OpenBLAS-*
   # we install both the serial and the omp threaded version.
   # Unfortunately, neither is thread-safe (i.e. the CP2K ssmp and psmp version need to link to something else, the omp version is unused)
   make -j $nprocs USE_THREAD=0 LIBNAMESUFFIX=serial PREFIX=${INSTALLDIR} >& make.serial.log
   make -j $nprocs USE_THREAD=0 LIBNAMESUFFIX=serial PREFIX=${INSTALLDIR} install >& install.serial.log
   # make clean >& clean.log
   # make -j $nprocs USE_THREAD=1 USE_OPENMP=1 LIBNAMESUFFIX=omp PREFIX=${INSTALLDIR} >& make.omp.log
   # make -j $nprocs USE_THREAD=1 USE_OPENMP=1 LIBNAMESUFFIX=omp PREFIX=${INSTALLDIR} install >& install.omp.log
   cd ..
fi
# currently openblas is not thread safe (neither serial nor omp version),
LIBS="IF_VALGRIND(-lreflapack -lrefblas, IF_OMP(-lreflapack -lrefblas,-lopenblas_serial)) ${LIBS}"

echo "================= Installing libsmm ==================="
if $enable_tsan ; then
   echo "TSAN build ... not downloading libsmm"
   libsmm=""
else
   # Here we attempt to determine which libsmm to download, and do that if it exists.
   # We use info on the architecture / core from the openblas build.

   # helper to check if libsmm is available (return 0) or not (return 8)
   libsmm_exists() {
      wget --spider http://www.cp2k.org/static/downloads/libsmm/$1 >& /dev/null
      echo $?
   }

   # where is the openblas configuration file, which gives us the core
   openblas_conf=`echo ${rootdir}/build/*OpenBLAS*/Makefile.conf`
   if [ ! -f "$openblas_conf" ]; then
      echo "Could not find OpenBLAS' Makefile.conf: $openblas_conf"
      exit 1
   fi
   openblas_libcore=`grep 'LIBCORE=' $openblas_conf | cut -f2 -d=`
   openblas_arch=`grep 'ARCH=' $openblas_conf | cut -f2 -d=`
   libsmm_libcore=libsmm_dnn_${openblas_libcore}.a
   tst=`libsmm_exists $libsmm_libcore`
   if [ "$tst" == "0" ]; then
      libsmm=$libsmm_libcore
      echo "An optimized libsmm $libsmm is available"
   else
      libsmm_arch=libsmm_dnn_${openblas_arch}.a
      tst=`libsmm_exists $libsmm_arch`
      if [ "$tst" == "0" ]; then
         libsmm=$libsmm_arch
         echo "A generic libsmm $libsmm is available."
         echo "Consider building and contributing to CP2K an optimized libsmm for your $openblas_arch $openblas_libcore"
      else
         echo "No libsmm is available"
         echo "Consider building and contributing to CP2K an optimized libsmm for your $openblas_arch $openblas_libcore"
         libsmm=""
      fi
   fi
fi

# we know what to get, proceed with install
if [ "$libsmm" != "" ]; then
   if [ -f $libsmm ]; then
      echo "Installation already started, skipping it."
   else
      wget http://www.cp2k.org/static/downloads/libsmm/$libsmm
      checksum $libsmm
      cp $libsmm ${INSTALLDIR}/lib/
      ln -s ${INSTALLDIR}/lib/$libsmm ${INSTALLDIR}/lib/libsmm_dnn.a
   fi
   DFLAGS="${DFLAGS} IF_VALGRIND(,-D__HAS_smm_dnn)"
   LIBS="IF_VALGRIND(,-lsmm_dnn) ${LIBS}"
fi


echo "================= Installing scalapack ==================="
if [ -f scalapack_installer.tgz ]; then
   echo "Installation already started, skipping it."
else
   wget http://www.cp2k.org/static/downloads/scalapack_installer.tgz
   checksum scalapack_installer.tgz
   tar -xzf scalapack_installer.tgz
   # we dont know the version
   cd scalapack_installer_*
   SLROOT=${PWD}
   # We use echo as mpirun command to avoid testing scalapack,
   # (lapack is still tested, and the --notesting flag causes lapack/blas not to be build, seemingly.)
   # yet download blas / lapack... whole installer is a bit serial as well (and fails with --make="make -j32")
   # also, doesn't seem to stop if something goes wrong in the build process..
   # finally, we should avoid -ffast-math as this seems to cause problems
   ./setup.py --mpirun=echo --downblas --downlapack --fcflags="$FCFLAGS -fno-fast-math" --ccflags="$CFLAGS -fno-fast-math" --ldflags_c="$LDFLAGS -fno-fast-math" --ldflags_fc="$LDFLAGS -fno-fast-math" >& make.log
   # copy libraries where we like them
   cp install/lib/* ${INSTALLDIR}/lib/
   cd ..
fi
DFLAGS="${DFLAGS} IF_MPI(-D__SCALAPACK,)"
LIBS="IF_MPI(-lscalapack,) ${LIBS}"


echo "==================== Installing libxc ===================="
if [ -f libxc-${libxc_ver}.tar.gz ]; then
   echo "Installation already started, skipping it."
else
   wget http://www.cp2k.org/static/downloads/libxc-${libxc_ver}.tar.gz
   checksum libxc-${libxc_ver}.tar.gz
   tar -xzf libxc-${libxc_ver}.tar.gz
   cd libxc-${libxc_ver}
   # patch buggy configure macro (fails with gcc trunk)
   sed -i 's/ax_cv_f90_modext=`ls | sed/ax_cv_f90_modext=`ls -1 | grep -iv smod | sed/g' configure
   ./configure  --prefix=${INSTALLDIR} >& config.log
   make -j $nprocs >& make.log
   make install >& install.log
   cd ..
fi
DFLAGS="${DFLAGS} -D__LIBXC2"
LIBS="-lxcf90 -lxc ${LIBS}"


echo "=================== Installing libint ===================="
if [ -f libint-${libint_ver}.tar.gz ]; then
   echo "Installation already started, skipping it."
else
   wget http://www.cp2k.org/static/downloads/libint-${libint_ver}.tar.gz
   checksum libint-${libint_ver}.tar.gz
   tar -xzf libint-${libint_ver}.tar.gz
   cd libint-${libint_ver}
   # hack for -with-cc (needed for -fsanitize=thread that also needs to be passed to the linker, but seemingly ldflags is ignored by libint's configure)
   ./configure  --prefix=${INSTALLDIR} --with-libint-max-am=5 --with-libderiv-max-am1=4 --with-cc="gcc $CFLAGS" --with-cc-optflags="$CFLAGS" --with-cxx-optflags="$CXXFLAGS" >& config.log
   make -j $nprocs >&  make.log
   make install >& install.log
   cd ..
fi
DFLAGS="${DFLAGS} -D__LIBINT -D__LIBINT_MAX_AM=6 -D__LIBDERIV_MAX_AM1=5"
LIBS="-lderiv -lint ${LIBS}"


echo "==================== Installing FFTW ====================="
if [ -f fftw-${fftw_ver}.tar.gz ]; then
   echo "Installation already started, skipping it."
else
   wget http://www.cp2k.org/static/downloads/fftw-${fftw_ver}.tar.gz
   checksum fftw-${fftw_ver}.tar.gz
   tar -xzf fftw-${fftw_ver}.tar.gz
   cd fftw-${fftw_ver}
   ./configure  --prefix=${INSTALLDIR} --enable-openmp >& config.log
   make -j $nprocs >&  make.log
   make install >& install.log
   cd ..
fi
DFLAGS="${DFLAGS} -D__FFTW3"
DFLAGS="${DFLAGS} IF_COVERAGE(IF_MPI(,-U__FFTW3),)" # also want to cover FFT_SG
LIBS="-lfftw3 IF_OMP(-lfftw3_omp,) ${LIBS}"


echo "==================== Installing ELPA ====================="
if [ -f elpa-${elpa_ver}.tar.gz ]; then
   echo "Installation already started, skipping it."
else
   wget http://www.cp2k.org/static/downloads/elpa-${elpa_ver}.tar.gz
   checksum elpa-${elpa_ver}.tar.gz
   tar -xzf elpa-${elpa_ver}.tar.gz

   # need both flavors ?
   cp -rp elpa-${elpa_ver} elpa-${elpa_ver}_mt

   # elpa expect FC to be an mpi fortran compiler that's happy with long lines, and that a bunch of libs can be found

   cd elpa-${elpa_ver}_mt
   FC="mpif90 -ffree-line-length-none" LDFLAGS="-L${INSTALLDIR}/lib" LIBS="-lscalapack -lreflapack -lrefblas" ./configure  --prefix=${INSTALLDIR} --enable-openmp=yes --with-generic --enable-shared=no >& config.log
   make -j $nprocs >&  make.log
   make install >& install.log
   cd ..

   cd elpa-${elpa_ver}
   FC="mpif90 -ffree-line-length-none" LDFLAGS="-L${INSTALLDIR}/lib" LIBS="-lscalapack -lreflapack -lrefblas" ./configure  --prefix=${INSTALLDIR} --enable-openmp=no --with-generic --enable-shared=no >& config.log
   make -j $nprocs >&  make.log
   make install >& install.log
   cd ..
fi
# Unfortunately, we need two separate include dirs for ELPA w/wo threading.
P1="-I\$(CP2KINSTALLDIR)/include/elpa_openmp-${elpa_ver}/modules"
P2="-I\$(CP2KINSTALLDIR)/include/elpa-${elpa_ver}/modules"
DFLAGS="${DFLAGS} IF_MPI(-D__ELPA2 IF_OMP(${P1},${P2}),)"
LIBS="IF_MPI(IF_OMP(-lelpa_openmp,-lelpa),) ${LIBS}"


echo "================== Installing PT-Scotch =================="
if [ -f scotch_${scotch_ver}.tar.gz ]; then
   echo "Installation already started, skipping it."
else
   wget  http://www.cp2k.org/static/downloads/scotch_${scotch_ver}.tar.gz
   checksum scotch_${scotch_ver}.tar.gz
   tar -xzf scotch_${scotch_ver}.tar.gz
   cd scotch_${scotch_ver}/src
   cat Make.inc/Makefile.inc.x86-64_pc_linux2 | \
   sed "s|\(^CFLAGS\).*|\1 =  $CFLAGS -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Drestrict=__restrict -DIDXSIZE64|" > Makefile.inc

   make scotch -j $nprocs >& make.log
   make ptscotch -j $nrocs >& make.log
   make install prefix=${INSTALLDIR} >& install.log
   cd ../..
fi
LIBS="IF_MPI(-lptscotch -lptscotcherr -lscotchmetis -lscotch -lscotcherr,) ${LIBS}"


echo "================== Installing ParMETIS =================="
if [ -f parmetis-${parmetis_ver}.tar.gz ]; then
   echo "Installation already started, skipping it."
else
   wget http://www.cp2k.org/static/downloads/parmetis-${parmetis_ver}.tar.gz
   checksum parmetis-${parmetis_ver}.tar.gz
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
LIBS="IF_MPI(-lptscotchparmetis,) ${LIBS}"
LIBSOMP="${LIBSOMP} IF_MPI(-lptscotchparmetis,)"


echo "================== Installing SuperLU_DIST =================="
if [ -f superlu_dist_${superlu_ver}.tar.gz ]; then
   echo "Installation already started, skipping it."
else
   wget http://www.cp2k.org/static/downloads/superlu_dist_${superlu_ver}.tar.gz
   checksum superlu_dist_${superlu_ver}.tar.gz
   tar -xzf superlu_dist_${superlu_ver}.tar.gz

   cd SuperLU_DIST_${superlu_ver}
   mv make.inc make.inc.orig
   cat <<EOF >> make.inc
PLAT=_x86_64
DSUPERLULIB= ${PWD}/lib/libsuperlu_dist_${superlu_ver}.a
LIBS=\$(DSUPERLULIB) -L${INSTALLDIR}/lib -lparmetis -lmetis -lrefblas
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
LIBS="IF_MPI(-lsuperlu_dist_${superlu_ver},) ${LIBS}"


echo "================== Installing PEXSI =================="
if [ -f pexsi_v${pexsi_ver}.tar.gz ]; then
   echo "Installation already started, skipping it."
else
   wget http://www.cp2k.org/static/downloads/pexsi_v${pexsi_ver}.tar.gz
   #wget https://math.berkeley.edu/~linlin/pexsi/download/pexsi_v${pexsi_ver}.tar.gz
   checksum pexsi_v${pexsi_ver}.tar.gz

   tar -xzf pexsi_v${pexsi_ver}.tar.gz

   cd pexsi_v${pexsi_ver}

   cat config/make.inc.linux.gnu | \
   sed 's|\(PAR_ND_LIBRARY *=\).*|\1 parmetis|' |\
   sed 's|\(SEQ_ND_LIBRARY *=\).*|\1 metis|' |\
   sed "s|\(PEXSI_DIR *=\).*|\1 ${PWD}|" |\
   sed "s|\(CPP_LIB *=\).*|\1 -lstdc++ ${mpiextralibs}|" |\
   sed "s|\(LAPACK_LIB *=\).*|\1 -L${INSTALLDIR}/lib -lreflapack|" |\
   sed "s|\(BLAS_LIB *=\).*|\1 -L${INSTALLDIR}/lib -lrefblas|" |\
   sed "s|\(\bMETIS_LIB *=\).*|\1 -L${INSTALLDIR}/lib -lmetis|" |\
   sed "s|\(PARMETIS_LIB *=\).*|\1 -L${INSTALLDIR}/lib -lparmetis|" |\
   sed "s|\(DSUPERLU_LIB *=\).*|\1 -L${INSTALLDIR}/lib -lsuperlu_dist_${superlu_ver}|" |\
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
DFLAGS="${DFLAGS} IF_MPI(-D__LIBPEXSI,)"
LIBS="IF_MPI(-lpexsi_linux_v${pexsi_ver},) ${LIBS}"


#echo "==================== Installing PLUMED ====================="
# Unfortunately plumed 2.x does not compile with gcc 5.x at the moment:
# https://groups.google.com/forum/#!msg/plumed-users/Y4q_7bx31ag/dNYdCa-LXZYJ
#if [ -f plumed-${plumed_ver}.tgz ]; then
#  echo "Installation already started, skipping it."
#else
#  wget http://www.cp2k.org/static/downloads/plumed/plumed-${plumed_ver}.tgz
#  checksum plumed-${plumed_ver}.tgz
#  tar -xzf plumed-${plumed_ver}.tgz
#  cd plumed-${plumed_ver}
#  ./configure  --prefix=${INSTALLDIR} >& config.log
#  make -j $nprocs >&  make.log
#  make install >& install.log
#  cd ..
#fi


echo "==================== Installing QUIP ================="
if $enable_tsan ; then
   echo "TSAN build ... will not use QUIP"
   libsmm=""
else
   if [ -f QUIP-${quip_ver}.zip  ]; then
      echo "Installation already started, skipping it."
   else
      wget http://www.cp2k.org/static/downloads/QUIP-${quip_ver}.zip
      checksum QUIP-${quip_ver}.zip
      unzip QUIP-${quip_ver}.zip >& unzip.log
      mv QUIP-public QUIP-${quip_ver}
      cd QUIP-${quip_ver}
      export QUIP_ARCH=linux_x86_64_gfortran
      # hit enter a few times to accept defaults
      echo -e "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" | make config > config.log
      # make -j does not work :-(
      make >& make.log
      cp build/linux_x86_64_gfortran/quip_unified_wrapper_module.mod  ${INSTALLDIR}/include/
      cp build/linux_x86_64_gfortran/*.a                              ${INSTALLDIR}/lib/
      cp src/FoX-4.0.3/objs.linux_x86_64_gfortran/lib/*.a             ${INSTALLDIR}/lib/
      cd ..
   fi
   LIBS="-lquip_core -latoms -lFoX_sax -lFoX_common -lFoX_utils -lFoX_fsys ${LIBS}"
   DFLAGS="${DFLAGS} -D__QUIP"
fi


echo "==================== generating arch files ===================="
echo "arch files can be found in the ${INSTALLDIR}/arch subdirectory"
mkdir -p ${INSTALLDIR}/arch
cd ${INSTALLDIR}/arch

# standart libs
LIBS="${LIBS} -lstdc++ "

# Flags which both gfortran and gcc understand.
BASEFLAGS="IF_OMP(-fopenmp,)"
BASEFLAGS="${BASEFLAGS} -march=native -fno-omit-frame-pointer -g ${TSANFLAGS}"
BASEFLAGS="${BASEFLAGS} IF_COVERAGE(-O0 -coverage, IF_DEBUG(-O1,-O3 -ffast-math))"
BASEFLAGS="${BASEFLAGS} IF_DEBUG(-fsanitize=leak -ffpe-trap='invalid,zero,overflow' -finit-real=snan -fno-fast-math -D__HAS_IEEE_EXCEPTIONS,)"
BASEFLAGS="${BASEFLAGS} \$(PROFOPT)"

# Special flags for gfortran
# https://gcc.gnu.org/onlinedocs/gfortran/Error-and-Warning-Options.html
# we error out for these warnings
WFLAGSERROR="-Werror=aliasing -Werror=ampersand -Werror=c-binding-type -Werror=intrinsic-shadow -Werror=intrinsics-std -Werror=line-truncation -Werror=tabs -Werror=realloc-lhs-all -Werror=target-lifetime -Werror=underflow -Werror=unused-but-set-variable -Werror=unused-variable -Werror=conversion"
# we just warn for those (that eventually might be promoted to WFLAGSERROR). It is useless to put something here with 100s of warnings.
WFLAGSWARN="-Wuse-without-only -Wzerotrip"
# while here we collect all other warnings, some we'll ignore
WFLAGSWARNALL="-pedantic -Wall -Wextra -Wsurprising -Wunused-dummy-argument -Wunused-parameter -Warray-temporaries -Wcharacter-truncation -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -Wreal-q-constant -Wunused-dummy-argument -Wunused-parameter -Walign-commons -Wfunction-elimination -Wrealloc-lhs -Wcompare-reals -Wzerotrip"
# combine warn/error flags
WFLAGS="$WFLAGSERROR $WFLAGSWARN IF_WARNALL(${WFLAGSWARNALL},)"
FCFLAGS="${BASEFLAGS} -I\$(CP2KINSTALLDIR)/include -std=f2003 -fimplicit-none -ffree-form IF_DEBUG(-fcheck='bounds,do,recursion,pointer',) ${WFLAGS} \$(DFLAGS)"
LDFLAGS="-L\$(CP2KINSTALLDIR)/lib \$(FCFLAGS)"

# Spcial flags for gcc (currently none)
CFLAGS="${BASEFLAGS} -I\$(CP2KINSTALLDIR)/include \$(DFLAGS)"

# CUDA stuff
LIBS="${LIBS} IF_CUDA(-lcudart -lcufft -lcublas -lrt IF_DEBUG(-lnvToolsExt,),)"
DFLAGS="IF_CUDA(-D__ACC -D__DBCSR_ACC -D__PW_CUDA IF_DEBUG(-D__CUDA_PROFILING,),) ${DFLAGS}"
NVFLAGS="-arch sm_35 \$(DFLAGS) "

# helper routine for instantiating the arch.tmpl
gen_arch_file() {
 filename=$1
 flags=$2
 TMPL=`cat ../../arch.tmpl`
 eval "printf \"$TMPL\"" | cpp -traditional-cpp -P ${flags} - > $filename
 echo "Wrote install/arch/"$filename
}

rm -f ${INSTALLDIR}/arch/local*
gen_arch_file "local.sopt"                 ""
gen_arch_file "local.sdbg"                 "-DDEBUG"
gen_arch_file "local.ssmp"                 "-DOMP"
gen_arch_file "local.popt"                 "-DMPI"
gen_arch_file "local.pdbg"                 "-DMPI -DDEBUG"
gen_arch_file "local.psmp"                 "-DMPI -DOMP"
gen_arch_file "local_cuda.ssmp"            "-DCUDA -DOMP"
gen_arch_file "local_cuda.psmp"            "-DCUDA -DOMP -DMPI"
gen_arch_file "local_cuda.sdbg"            "-DCUDA -DDEBUG -DOMP"
gen_arch_file "local_cuda.pdbg"            "-DCUDA -DDEBUG -DOMP -DMPI"
gen_arch_file "local_valgrind.sdbg"        "-DVALGRIND"
gen_arch_file "local_valgrind.pdbg"        "-DVALGRIND -DMPI"
gen_arch_file "local_coverage.sdbg"        "-DCOVERAGE -DDEBUG"
gen_arch_file "local_coverage.pdbg"        "-DCOVERAGE -DDEBUG -DMPI"
gen_arch_file "local_coverage_cuda.pdbg"   "-DCOVERAGE -DDEBUG -DMPI -DCUDA"
gen_arch_file "local_cuda_warn.psmp"       "-DCUDA -DMPI -DOMP -DWARNALL"

echo "========================== usage ========================="
echo "done!"
echo "now copy: cp ${INSTALLDIR}/arch/* to the cp2k/arch/ directory"
echo "to use this toolchain or the cp2k version compiled with it you will first need to execute at the prompt:"
echo "source ${SETUPFILE}"
echo "to build CP2K you should change directory cd cp2k/makefiles/"
echo "make -j${nprocs}" 'ARCH=local VERSION="sdbg sopt ssmp popt pdbg psmp"'

#EOF
