# Dockerfile for CP2K continuous integration (CI) runs
#
# A stand-alone build in this folder can be performed with:
# podman build --build-arg DEPS_IMAGE=<image id> --shm-size=1g -f build_cp2k_spack.Dockerfile ../../
#
# Stage 2: Build CP2K

ARG BASE_IMAGE="docker.io/ubuntu:26.04"
ARG DEPS_IMAGE

FROM "${DEPS_IMAGE}" AS build_cp2k

# Number of parallel build jobs (override with --build-arg NUM_PROCS=...)
ARG NUM_PROCS="32"

# Build CP2K
WORKDIR /opt/cp2k
COPY ./CMakeLists.txt ./CMakePresets.json ./
COPY ./cmake ./cmake
COPY ./data ./data
COPY ./src ./src
COPY ./tools/build_utils ./tools/build_utils

RUN ./make_cp2k.sh -cray -cv ${CP2K_VERSION} -j${NUM_PROCS} -uc no ${FEATURE_FLAGS}

# Stage 3: Install CP2K

FROM "${BASE_IMAGE}" AS install_cp2k

ARG APT_GET="apt-get -o APT::Sandbox::User=root -qq"

RUN ${APT_GET} update && \
    ${APT_GET} install --no-install-recommends \
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
COPY ./ci ./ci

# Do not rely only on LD_LIBRARY_PATH because it is fragile
COPY --from=build_cp2k /etc/ld.so.conf.d/cp2k.conf /etc/ld.so.conf.d/cp2k.conf
RUN ldconfig

# Create entrypoint and finalise container build
WORKDIR /mnt
ENTRYPOINT ["/opt/cp2k/install/bin/launch"]
CMD ["cp2k", "--help", "--version"]
