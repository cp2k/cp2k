#!/bin/bash
###############################################################################
# Copyright (c) Intel Corporation - All rights reserved.                      #
# This file is part of the XCONFIGURE project.                                #
#                                                                             #
# For information on the license, see the LICENSE file.                       #
# Further information: https://github.com/hfp/xconfigure/                     #
# SPDX-License-Identifier: BSD-3-Clause                                       #
###############################################################################
# Hans Pabst (Intel Corp.)
###############################################################################

# number of systems (clusters nodes)
TOTALNUMNODES=1
# number of physical cores per node
NCORESPERNODE=1
# number of sockets per system
NPROCSPERNODE=1
# number of threads per core
NTHREADSPERCORE=1
# number of ranks per node
PENALTY_MIN=1
# unbalanced rank-count
PENALTY_ODD=3

GREP=$(command -v grep)
SORT=$(command -v sort)
HEAD=$(command -v head)
SEQ=$(command -v seq)
CUT=$(command -v cut)
TR=$(command -v tr)

if [ "" != "${HOME}" ]; then
  CONFIGFILE=${HOME}/.xconfigure-cp2k-plan
else
  HERE=$(
    cd $(dirname $0)
    pwd -P
  )
  CONFIGFILE=${HERE}/.xconfigure-cp2k-plan
fi

function isqrt {
  s=1073741824
  x=$1
  y=0
  while [ "0" != "$((0 < s))" ]; do
    b=$((y | s))
    y=$((y >> 1))
    if [ "0" != "$((b <= x))" ]; then
      x=$((x - b))
      y=$((y | s))
    fi
    s=$((s >> 2))
  done
  echo "${y}"
}

if [ "" != "${GREP}" ] && [ "" != "${SORT}" ] && [ "" != "${HEAD}" ] &&
  [ "" != "${SEQ}" ] && [ "" != "${CUT}" ] && [ "" != "${TR}" ]; then
  HELP=0
  if [ "--help" = "$1" ] || [ "-help" = "$1" ] || [ "-h" = "$1" ]; then
    HELP=1
    shift
  elif [ "" != "$1" ]; then
    TOTALNUMNODES=$1
    shift
  fi
  if [ -e /proc/cpuinfo ] && [ "" != "$(command -v wc)" ]; then
    NS=$(${GREP} "physical id" /proc/cpuinfo | ${SORT} -u | wc -l | ${TR} -d " ")
    NC=$((NS * $(${GREP} -m1 "cpu cores" /proc/cpuinfo | ${TR} -d " " | ${CUT} -d: -f2)))
    NT=$(${GREP} "core id" /proc/cpuinfo | wc -l | ${TR} -d " ")
  elif [ "Darwin" = "$(uname)" ] && [ "" != "$(command -v sysctl)" ]; then
    NS=$(sysctl hw.packages | ${CUT} -d: -f2 | ${TR} -d " ")
    NC=$(sysctl hw.physicalcpu | ${CUT} -d: -f2 | ${TR} -d " ")
    NT=$(sysctl hw.logicalcpu | ${CUT} -d: -f2 | ${TR} -d " ")
  fi
  if [ "" != "${NC}" ] && [ "" != "${NT}" ]; then
    HT=$((NT / NC))
  fi
  OUTPUT=0
  if [ "" = "$1" ]; then
    if [ -e ${CONFIGFILE} ]; then # remind configuration
      NCORESPERNODE=$(${CUT} -d" " -f1 ${CONFIGFILE})
    elif [ "" != "${NC}" ]; then
      NCORESPERNODE=${NC}
    fi
  else
    NCORESPERNODE=$1
    OUTPUT=1
    shift
  fi
  if [ "" = "$1" ]; then
    if [ -e ${CONFIGFILE} ]; then # remind configuration
      NTHREADSPERCORE=$(${CUT} -d" " -f2 ${CONFIGFILE})
    elif [ "" != "${HT}" ]; then
      NTHREADSPERCORE=${HT}
    fi
  else
    NTHREADSPERCORE=$1
    OUTPUT=1
    shift
  fi
  if [ "" = "$1" ]; then
    if [ -e ${CONFIGFILE} ]; then # remind configuration
      NPROCSPERNODE=$(${CUT} -d" " -f3 ${CONFIGFILE})
    elif [ "" != "${NS}" ]; then
      NPROCSPERNODE=${NS}
    fi
  else
    NPROCSPERNODE=$1
    OUTPUT=1
    shift
  fi
  if [ "0" != "${HELP}" ]; then
    echo "Run: $0 [num-nodes [ncores-per-node [nthreads-per-core [nsockets-per-node]]]]"
    echo
    echo "Defaults: [num-nodes] ${TOTALNUMNODES}"
    echo "    [ncores-per-node] ${NCORESPERNODE}"
    echo "  [nthreads-per-core] ${NTHREADSPERCORE}"
    echo "  [nsockets-per-node] ${NPROCSPERNODE}"
    echo
    echo "Settings (except num-nodes) may be setup once and are persistent (run to run)."
    echo "Default of the first run are adopted from the system running $0."
    exit 0
  fi
  # min. number of ranks per node
  MIN_NRANKS=$((PENALTY_MIN * NPROCSPERNODE))
  # remember system configuration
  if [ "0" != "${OUTPUT}" ]; then
    echo "${NCORESPERNODE} ${NTHREADSPERCORE} ${NPROCSPERNODE}" > ${CONFIGFILE} 2> /dev/null
  fi
  NCORESTOTAL=$((TOTALNUMNODES * NCORESPERNODE))
  NCORESOCKET=$((NCORESPERNODE / NPROCSPERNODE))
  echo "================================================================================"
  echo "${NCORESTOTAL} cores: ${TOTALNUMNODES} node(s) with ${NPROCSPERNODE}x${NCORESOCKET} core(s) per node and ${NTHREADSPERCORE} thread(s) per core"
  echo "================================================================================"
  NRANKSMIN=$((TOTALNUMNODES * NPROCSPERNODE))
  NSQRT_MIN=$(isqrt ${NRANKSMIN})
  NSQRT_MAX=$(isqrt ${NCORESTOTAL})
  for NSQRT in $(${SEQ} ${NSQRT_MIN} ${NSQRT_MAX}); do
    NSQR=$((NSQRT * NSQRT))
    NRANKSPERNODE=$((NSQR / TOTALNUMNODES))
    if [ "${NSQR}" = "$((TOTALNUMNODES * NRANKSPERNODE))" ]; then
      REST=$((NCORESPERNODE % NRANKSPERNODE))
      ODD=$(((NRANKSPERNODE % NPROCSPERNODE) != 0))
      # criterion to add penalty in case of unbalanced load
      if [ "0" != "$((PENALTY_ODD * MIN_NRANKS * REST <= NCORESPERNODE))" ] || [ "0" = "${ODD}" ]; then
        if [ "0" != "$((MIN_NRANKS * REST <= NCORESPERNODE))" ] &&
          [ "0" != "$((MIN_NRANKS <= NRANKSPERNODE))" ]; then
          PENALTY=$(((100 * REST + ODD * PENALTY_ODD + NCORESPERNODE - 1) / NCORESPERNODE))
          RESULTS+="${NRANKSPERNODE};${PENALTY};${NSQRT}\n"
        fi
      fi
    fi
  done
  RESULTS=$(echo -e "${RESULTS}" | ${GREP} -v "^$" | ${SORT} -t";" -u -k2n -k1nr)
  NRANKSPERNODE_TOP=$(echo "${RESULTS}" | ${CUT} -d";" -f1 | ${HEAD} -n1)
  NTHREADSPERNODE=$((NCORESPERNODE * NTHREADSPERCORE))
  NSQR_MAX=$((NSQRT_MAX * NSQRT_MAX))
  PENALTY_NCORES=$((NCORESTOTAL - NSQR_MAX))
  PENALTY_TOP=$(((100 * PENALTY_NCORES + NCORESTOTAL - 1) / NCORESTOTAL))
  NRANKSPERNODE=${NCORESPERNODE}
  OUTPUT_POT=0
  while [ "0" != "$((NRANKSPERNODE_TOP < NRANKSPERNODE))" ]; do
    # criterion to add penalty in case of unbalanced load
    ODD=$(((NRANKSPERNODE % NPROCSPERNODE) != 0))
    if [ "0" != "$((PENALTY_ODD * MIN_NRANKS * PENALTY_NCORES <= NCORESTOTAL))" ] || [ "0" = "${ODD}" ]; then
      NTHREADSPERRANK=$((NTHREADSPERNODE / NRANKSPERNODE))
      if [ "0" != "$((MIN_NRANKS * PENALTY_NCORES <= NCORESTOTAL))" ] &&
        [ "0" != "$((MIN_NRANKS <= NRANKSPERNODE))" ]; then
        NRANKS_CUR=$((TOTALNUMNODES * NRANKSPERNODE))
        SQRTNRANKS=$(isqrt $((NRANKS_CUR)))
        NRANKS_SQR=$((SQRTNRANKS * SQRTNRANKS))
        NRANKS_COM=$((NRANKSPERNODE < NRANKS_SQR ? (NRANKS_SQR / NRANKSPERNODE * NRANKSPERNODE) : NRANKS_SQR))
        PENALTY=$((100 * (NRANKS_CUR + ODD * PENALTY_ODD - NRANKS_COM) / NRANKS_COM))
        echo "[${NRANKSPERNODE}x${NTHREADSPERRANK}]: ${NRANKSPERNODE} ranks per node with ${NTHREADSPERRANK} thread(s) per rank (${PENALTY}% penalty)"
        if [ "0" != "$((PENALTY_TOP < PENALTY))" ]; then PENALTY_TOP=${PENALTY}; fi
        OUTPUT_POT=$((OUTPUT_POT + 1))
      fi
    fi
    NRANKSPERNODE=$((NRANKSPERNODE >> 1))
  done
  if [ "0" != "${OUTPUT_POT}" ]; then
    echo "--------------------------------------------------------------------------------"
  fi
  OUTPUT_SQR=0
  # reorder by decreasing rank-count
  RESULTS=$(echo -e "${RESULTS}" | ${TR} " " "\n" | ${SORT} -t";" -k1nr -k2n)
  for RESULT in ${RESULTS}; do
    NRANKSPERNODE=$(echo "${RESULT}" | ${CUT} -d";" -f1)
    NTHREADSPERRANK=$((NTHREADSPERNODE / NRANKSPERNODE))
    PENALTY=$(echo "${RESULT}" | ${CUT} -d";" -f2)
    if [ "0" != "$((OUTPUT_SQR < OUTPUT_POT))" ] ||
      [ "0" != "$((PENALTY <= PENALTY_TOP))" ]; then
      NSQRT=$(echo "${RESULT}" | ${CUT} -d";" -f3)
      echo "[${NRANKSPERNODE}x${NTHREADSPERRANK}]: ${NRANKSPERNODE} ranks per node with ${NTHREADSPERRANK} thread(s) per rank (${PENALTY}% penalty) -> ${NSQRT}x${NSQRT}"
      OUTPUT_SQR=$((OUTPUT_SQR + 1))
    fi
  done
  if [ "0" != "${OUTPUT_SQR}" ]; then
    echo "--------------------------------------------------------------------------------"
  fi
  for NRANKSPERNODE in $(${SEQ} ${MIN_NRANKS} ${NCORESPERNODE}); do
    REST=$((NCORESPERNODE % NRANKSPERNODE))
    PENALTY=$(((100 * REST + NCORESPERNODE - 1) / NCORESPERNODE))
    # criterion to add penalty in case of unbalanced load
    ODD=$(((NRANKSPERNODE % NPROCSPERNODE) != 0))
    if [ "0" != "$((PENALTY_ODD * MIN_NRANKS * REST <= NCORESPERNODE))" ] || [ "0" = "${ODD}" ]; then
      if [ "0" != "$((MIN_NRANKS * REST <= NCORESPERNODE))" ] &&
        [ "0" != "$((MIN_NRANKS <= NRANKSPERNODE))" ]; then
        SQRT=$(isqrt $((TOTALNUMNODES * NRANKSPERNODE)))
        NUMCORES=$((SQRT * SQRT))
        NUMNODES=$((NUMCORES / NRANKSPERNODE))
        if [ "0" = "$((NUMCORES - NUMNODES * NRANKSPERNODE))" ]; then
          SUGGEST_LO+="${NUMNODES}\n"
        fi
        NUMCORES=$(((SQRT - 1) * (SQRT - 1)))
        NUMNODES=$((NUMCORES / NRANKSPERNODE))
        if [ "0" = "$((NUMCORES - NUMNODES * NRANKSPERNODE))" ]; then
          SUGGEST_LO+="${NUMNODES}\n"
        fi
        NUMCORES=$(((SQRT + 1) * (SQRT + 1)))
        NUMNODES=$((NUMCORES / NRANKSPERNODE))
        if [ "0" = "$((NUMCORES - NUMNODES * NRANKSPERNODE))" ]; then
          SUGGEST_HI+="${NUMNODES}\n"
        fi
        NUMCORES=$(((SQRT + 2) * (SQRT + 2)))
        NUMNODES=$((NUMCORES / NRANKSPERNODE))
        if [ "0" = "$((NUMCORES - NUMNODES * NRANKSPERNODE))" ] && [ "" = "${SUGGEST_HI}" ]; then
          SUGGEST_HI+="${NUMNODES}\n"
        fi
      fi
    fi
  done
  SUGGEST_LO=$(echo -e "${SUGGEST_LO}" |
    ${GREP} -vw "${TOTALNUMNODES}" |
    ${SORT} -nr | ${GREP} -v "^$" -m1)
  SUGGEST_HI=$(echo -e "${SUGGEST_HI}" |
    ${GREP} -vw "${TOTALNUMNODES}" |
    ${SORT} -n | ${GREP} -v "^$" -m1)
  if [ "0" != "${SUGGEST_LO}" ]; then
    echo "Try also the following node counts: ${SUGGEST_LO} ${SUGGEST_HI}"
  else
    echo "Try also the following node counts: ${SUGGEST_HI}"
  fi
else
  echo "Error: missing prerequisites!"
  exit 1
fi
