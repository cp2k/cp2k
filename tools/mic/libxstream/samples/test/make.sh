#!/bin/bash

LIBXSTREAM_ROOT="../.."
NAME=$(basename ${PWD})

ICCOPT="-O2 -xHost -ansi-alias"
ICCLNK=""

GCCOPT="-O2 -march=native"
GCCLNK=""

OPT="-Wall -std=c++0x"

if [[ "" = "${CXX}" ]] ; then
  CXX=$(which icpc 2> /dev/null)
  if [[ "" != "${CXX}" ]] ; then
    OPT+=" ${ICCOPT}"
    LNK+=" ${ICCLNK}"
  else
    CXX="g++"
    OPT+=" ${GCCOPT}"
    LNK+=" ${GCCLNK}"
  fi
else
  OPT+=" ${GCCOPT}"
  LNK+=" ${GCCLNK}"
fi

if [ "-g" = "$1" ] ; then
  OPT+=" -O0 -g"
  shift
else
  OPT+=" -DNDEBUG"
fi

if [[ "Windows_NT" = "${OS}" ]] ; then
  OPT+=" -D_REENTRANT"
  LNK+=" -lpthread"
else
  OPT+=" -pthread"
fi

${CXX} ${OPT} $* \
  -I${LIBXSTREAM_ROOT}/include -I${LIBXSTREAM_ROOT}/src -DLIBXSTREAM_EXPORTED \
  ${LIBXSTREAM_ROOT}/src/*.cpp *.cpp \
  ${LNK} -o ${NAME}
