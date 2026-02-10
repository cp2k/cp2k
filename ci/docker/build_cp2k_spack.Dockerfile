# Dockerfile for CP2K continuous integration (CI) runs
#
# A stand-alone build in this folder can be performed with:
# podman build --build-arg DEPS_IMAGE=<image id> --shm-size=1g -f build_cp2k_spack.Dockerfile ../../
#
# Author: Matthias Krack (MK)
#

ARG BASE_IMAGE="ubuntu:24.04"
ARG DEPS_IMAGE=""

###### Stage 2a: Build CP2K ######

FROM "${DEPS_IMAGE}" AS build_cp2k

# Retrieve the number of available CPU cores
ARG NUM_PROCS
ENV NUM_PROCS=${NUM_PROCS:-32}

# Build CP2K
WORKDIR /opt/cp2k
RUN /bin/bash -o pipefail -c "source ./make_cp2k.sh -cv ${CP2K_VERSION} -dlc -j${NUM_PROCS}"

###### Stage 2: Install CP2K ######

FROM "${BASE_IMAGE}" AS install_cp2k

RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    g++ gcc gfortran \
    hwloc \
    libhwloc-dev \
    python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/cp2k

# Install CP2K dependencies built with spack
COPY --from=build_cp2k /opt/cp2k/spack/spack-1.1.0/opt/spack ./spack/spack-1.1.0/opt/spack

# Install CP2K
COPY --from=build_cp2k /opt/cp2k/install ./install

# Install CP2K regression tests
COPY --from=build_cp2k /opt/cp2k/tests ./tests
COPY --from=build_cp2k /opt/cp2k/src/grid/sample_tasks ./src/grid/sample_tasks

# Install CP2K/Quickstep CI benchmarks
COPY --from=build_cp2k /opt/cp2k/benchmarks/CI ./benchmarks/CI

# Create entrypoint and finalise container build
WORKDIR /mnt
ENTRYPOINT ["/opt/cp2k/install/bin/launch"]
CMD ["cp2k", "--help", "--version"]
