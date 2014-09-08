#!/bin/bash -e

#set -x

TMPDIR=`mktemp -d`
REV=`./build_utils/get_revision_number ../ | sed -e "s/://"`
OUTNAME=libdbcsr_$REV
OUTDIR=$TMPDIR/$OUTNAME

#echo "Create tempdir: " $OUTDIR
mkdir $OUTDIR
cp -a ../COPYRIGHT $OUTDIR

mkdir $OUTDIR/src
rsync -axq \
   --exclude=".*" \
   --exclude="**/preprettify/" \
   --exclude="*.instantiation" \
   --exclude="*.template" \
   --exclude="*.pyc" \
   --exclude="libcusmm.cu" \
   ../src/dbcsr/ $OUTDIR/src/dbcsr

rsync -axq \
   --exclude=".*" \
   --exclude="**/preprettify/" \
   --exclude="*.instantiation" \
   --exclude="*.template" \
   --exclude="*.pyc" \
   --exclude="libcusmm.cu" \
   ../src/acc/ $OUTDIR/src/acc

rsync -axq \
   --exclude=".*" \
   --exclude="**/preprettify/" \
   --exclude="*.instantiation" \
   --exclude="*.template" \
   --exclude="*.pyc" \
   --exclude="libcusmm.cu" \
   ../src/base/ $OUTDIR/src/base

rsync -axq \
   --exclude=".*" \
   --exclude="**/preprettify/" \
   --exclude="*.instantiation" \
   --exclude="*.template" \
   --exclude="*.pyc" \
   --exclude="libcusmm.cu" \
   ../src/mpiwrap/ $OUTDIR/src/mpiwrap


mkdir -p $OUTDIR/tools/build_utils
rsync -axq --exclude=".*" ../tools/build_utils/  $OUTDIR/tools/build_utils/

mkdir $OUTDIR/tools/dbcsr_test/
rsync -axq --exclude=".*" ../tools/dbcsr_test/  $OUTDIR/tools/dbcsr_test/

mkdir $OUTDIR/tools/build_libsmm/
rsync -axq \
   --exclude=".*" \
   --exclude="**/run_*/" \
   --exclude="**/lib/" \
   --exclude="*.out" \
   --exclude="*.o" \
   --exclude="*.mod" \
   ../tools/build_libsmm/  $OUTDIR/tools/build_libsmm/

mkdir $OUTDIR/makefiles
cp -a ../makefiles/Makefile  $OUTDIR/makefiles/

mkdir $OUTDIR/arch/
rsync -axq --exclude=".*" ../arch/  $OUTDIR/arch/

# pack everything into a tar.gz-archive
rm -f $OUTNAME.tgz
tar -cz --directory=$TMPDIR --file=$OUTNAME.tgz $OUTNAME

#clean up
rm -r $TMPDIR

echo "Assembled: " $OUTNAME.tgz
#EOF
