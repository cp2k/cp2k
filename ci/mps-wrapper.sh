#!/bin/bash

set -eu

mps_prefix="/tmp/$(id -un)/slurm-${SLURM_JOBID}.${SLURM_STEPID}/nvidia"
num_gpus=4

# Reset CUDA environment variables to default values without MPS
export CUDA_DEVICE_MAX_CONNECTIONS=8
export CUDA_DEVICE_MAX_COPY_CONNECTIONS=8

# Each GPU is attached to the corresponding NUMA node
# Disable HWLOC_KEEP_NVIDIA_GPU_NUMA_NODES to avoid GPU NUMA nodes appearing in the list of CUDA devices
numa_node=$(HWLOC_KEEP_NVIDIA_GPU_NUMA_NODES=0 hwloc-calc --physical --intersect NUMAnode "$(hwloc-bind --get --taskset)")

# We expect exactly one non-negative integer for the NUMA node
if ! [[ "${numa_node}" =~ ^[0-9]+$ ]]; then
  echo "The MPS wrapper script only works when the process mask of the rank is adjacent to exactly one GPU. The CPU mask is $(hwloc-bind --get --taskset) and the list of adjacent numa nodes is ${numa_node} for rank ${SLURM_PROCID}."
  exit 1
fi

function start_daemon {
  export CUDA_MPS_PIPE_DIRECTORY=${mps_prefix}-mps-${1}
  export CUDA_MPS_LOG_DIRECTORY=${mps_prefix}-log-${1}
  mkdir -p "${CUDA_MPS_PIPE_DIRECTORY}"
  mkdir -p "${CUDA_MPS_LOG_DIRECTORY}"
  CUDA_VISIBLE_DEVICES=${1} nvidia-cuda-mps-control -d
}

# Start MPS control daemons from a single rank per node, one for each GPU on the node
if [[ $SLURM_LOCALID -eq 0 ]]; then
  # We attempt to kill previous MPS instances, but if we can't (either none
  # have been started or they are unkillable) we ignore it and attempt to run
  # anyway
  pkill --uid "$(id -un)" '^nvidia-cuda-mps-' || true

  for i in $(seq 0 $((num_gpus - 1))); do
    start_daemon "${i}"
  done
fi

# Set MPS options for this rank. Each rank will access the MPS of the GPU
# corresponding to the NUMA node. CUDA_VISIBLE_DEVICES should not be set. The
# chosen MPS determines which device is visible.
export CUDA_MPS_PIPE_DIRECTORY=${mps_prefix}-mps-${numa_node}
export CUDA_MPS_LOG_DIRECTORY=${mps_prefix}-log-${numa_node}

# Wait until the control daemon for our rank is up. The daemon creates a pid
# file which we can wait for. See
# https://docs.nvidia.com/deploy/mps/appendix-tools-and-interface-reference.html.
# Wait up to mps_pid_file_timeout seconds for the pid file to be created. In
# jobs with a large number of ranks some ranks may take a long time to start. If
# that happens consider increasing the timeout.
mps_pid_file_timeout=120
pid_file="${CUDA_MPS_PIPE_DIRECTORY}/nvidia-cuda-mps-control.pid"
if ! timeout ${mps_pid_file_timeout} bash -c "until [[ -f \"${pid_file}\" ]]; do sleep 1; done"; then
  echo "The MPS wrapper script timed out waiting for MPS pid file ${pid_file} on rank ${SLURM_PROCID}. MPS daemon likely did not start correctly or the rank starting the MPS daemons took too long to start."
  exit 1
fi

# Run the command
# We are using `exec`, because we want e.g. signals to be forwarded directly to the application, and not this wrapper script
exec "$@"
