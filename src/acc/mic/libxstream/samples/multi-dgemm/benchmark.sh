#!/bin/bash

#TRY="echo"
FILE="benchmark.txt"
STREAMS=2
SIZE=250
STRIDE=1
FRACT=8

if [ "" != "$1" ] ; then
  SIZE=$1
  shift
fi

if [ "" != "$1" ] ; then
  STRIDE=$1
  shift
fi

BSIZE=$((SIZE / FRACT))
BATCH=${STRIDE}

cat /dev/null > ${FILE}
while [[ ${BATCH} -le ${BSIZE} ]] ; do
  ${TRY} ./multi-dgemm.sh ${SIZE} ${BATCH} $* >> ${FILE}
  BATCH=$((BATCH + STRIDE))
done
