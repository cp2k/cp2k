# Dockerfile for CP2K continuous integration (CI) runs
#
# A stand-alone build in this folder can be performed with:
# podman build --build-arg DEPS_IMAGE=<image id> --shm-size=1g -f build_cp2k_spack.Dockerfile ../../
#
# Stage 2: Build CP2K

ARG BASE_IMAGE=${BASE_IMAGE:-ubuntu:26.04}
ARG DEPS_IMAGE=${DEPS_IMAGE:-}

FROM "${DEPS_IMAGE}" AS build_cp2k

# Setup CUDA environment
ENV CUDA_HOME=/usr/local/cuda
ENV LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"

# Retrieve the number of available CPU cores
ARG NUM_PROCS
ENV NUM_PROCS=${NUM_PROCS:-32}

ARG FEATURE_FLAGS
ENV FEATURE_FLAGS=${FEATURE_FLAGS:-}

# Build CP2K
WORKDIR /opt/cp2k
COPY ./src ./src
COPY ./data ./data
COPY ./tools/build_utils ./tools/build_utils
COPY ./cmake ./cmake
COPY ./CMakeLists.txt .

RUN ./make_cp2k.sh -cray -cv ${CP2K_VERSION} -uc no -j${NUM_PROCS} ${FEATURE_FLAGS}

# Stage 3: Install CP2K

FROM "${BASE_IMAGE}" AS install_cp2k

RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    g++ gcc gfortran \
    python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/cp2k

# Install CP2K dependencies built with spack
COPY --from=build_cp2k /opt/cp2k/spack/spack/opt/spack ./spack/spack/opt/spack

# Install CP2K
COPY --from=build_cp2k /opt/cp2k/install ./install

# Install CP2K regression tests
COPY ./tests ./tests
COPY --from=build_cp2k /opt/cp2k/src/grid/sample_tasks ./src/grid/sample_tasks

# Install CP2K/Quickstep CI benchmarks
COPY ./benchmarks/CI ./benchmarks/CI

# Do not rely only on LD_LIBRARY_PATH because it is fragile
COPY --from=build_cp2k /etc/ld.so.conf.d/cp2k.conf /etc/ld.so.conf.d/cp2k.conf
RUN ldconfig

# Create entrypoint and finalise container build
WORKDIR /mnt
ENTRYPOINT ["/opt/cp2k/install/bin/launch"]
CMD ["cp2k", "--help", "--version"]
