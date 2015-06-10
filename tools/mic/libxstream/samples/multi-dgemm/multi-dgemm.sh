#!/bin/bash

env \
  KMP_AFFINITY=scatter \
  OFFLOAD_INIT=on_start \
  MIC_USE_2MB_BUFFERS=2m \
  MIC_ENV_PREFIX=MIC \
  MIC_KMP_AFFINITY=balanced,granularity=fine \
./multi-dgemm $*
