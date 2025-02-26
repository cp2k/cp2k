# Dockerfile for CP2K continous integration (CI) runs
#
# A stand-alone docker build in this folder can be performed using the command:
# DOCKER_BUILDKIT=0 docker build -f build_cp2k_toolchain_psmp.Dockerfile ../../
#
# Author: Matthias Krack
#
# Stage 1a: Build CP2K toolchain

ARG BASE_IMAGE="ubuntu:24.04"

FROM ${BASE_IMAGE} AS build

# Install packages required for the CP2K toolchain build
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    bzip2 \
    ca-certificates \
    g++ \
    gcc \
    gfortran \
    git \
    libtool \
    libtool-bin \
    make \
    patch \
    pkg-config \
    python3 \
    unzip \
    wget \
    zlib1g-dev

# Retrieve the number of available CPU cores
ARG NUM_PROCS=64

# Install an MPI version ABI-compatible with the host MPI on Cray at CSCS
ARG MPICH_VERSION=3.1.4
RUN wget -q https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
    && tar -xf mpich-${MPICH_VERSION}.tar.gz \
    && cd mpich-${MPICH_VERSION} \
    && ./configure --prefix="/usr/local" \
       FFLAGS="${FFLAGS} -fallow-argument-mismatch" \
       FCFLAGS="${FCFLAGS} -fallow-argument-mismatch" \
    && make -j ${NUM_PROCS} \
    && make install \
    && ldconfig \
    && cd .. \
    && rm -rf mpich-${MPICH_VERSION} \
    && rm mpich-${MPICH_VERSION}.tar.gz

# Build CP2K toolchain
ARG BUILD_TYPE="minimal"
COPY ./tools/toolchain /opt/cp2k/tools/toolchain
WORKDIR /opt/cp2k/tools/toolchain
RUN echo "\nBuild type: ${BUILD_TYPE}\n" && \
    if [ "${BUILD_TYPE}" = "minimal" ]; then \
      ./install_cp2k_toolchain.sh -j ${NUM_PROCS} \
        --dry-run \
        --no-arch-files \
        --target-cpu=native \
        --with-gcc=system \
        --with-mpich=system \
        --with-ninja \
        --with-cosma=no \
        --with-dftd4=no \
        --with-elpa=no \
        --with-libgrpp=no \
        --with-libint=no \
        --with-libvori=no \
        --with-libxc=no \
        --with-sirius=no \
        --with-spglib=no; \
    elif [ "${BUILD_TYPE}" = "toolchain" ]; then \
      ./install_cp2k_toolchain.sh -j ${NUM_PROCS} \
        --dry-run \
        --install-all \
        --no-arch-files \
        --target-cpu=native \
        --with-gcc=system \
        --with-mpich=system; \
    fi

# Perform toolchain build step-wise in stages after its initialization with dry-run
COPY ./tools/toolchain/scripts/stage0/ ./scripts/stage0/
RUN  ./scripts/stage0/install_stage0.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage1/ ./scripts/stage1/
RUN  ./scripts/stage1/install_stage1.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage2/ ./scripts/stage2/
RUN  ./scripts/stage2/install_stage2.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage3/ ./scripts/stage3/
RUN  ./scripts/stage3/install_stage3.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage4/ ./scripts/stage4/
RUN  ./scripts/stage4/install_stage4.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage5/ ./scripts/stage5/
RUN  ./scripts/stage5/install_stage5.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage6/ ./scripts/stage6/
RUN  ./scripts/stage6/install_stage6.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage7/ ./scripts/stage7/
RUN  ./scripts/stage7/install_stage7.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage8/ ./scripts/stage8/
RUN  ./scripts/stage8/install_stage8.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage9/ ./scripts/stage9/
RUN  ./scripts/stage9/install_stage9.sh && rm -rf ./build

# Stage 1b: Build CP2K with CMake
WORKDIR /opt/cp2k
COPY ./CMakeLists.txt .
COPY ./benchmarks/CI ./benchmarks/CI
COPY ./cmake ./cmake
COPY ./data ./data
COPY ./src ./src
COPY ./tests ./tests
COPY ./tools/build_utils ./tools/build_utils

# Run CMake
RUN /bin/bash -c -o pipefail " \
    TOOLCHAIN_DIR=/opt/cp2k/tools/toolchain; \
    source ./cmake/cmake_cp2k.sh "${BUILD_TYPE}" psmp |& tee cmake.log"

# Compile CP2K
WORKDIR /opt/cp2k/build
RUN /bin/bash -c -o pipefail " \
    source /opt/cp2k/tools/toolchain/install/setup; \
    echo -e '\nCompiling CP2K ... \c'; \
    if ninja --verbose | tee ninja.log; then \
      echo -e 'done\n'; \
      echo -e 'Installing CP2K ... \c'; \
      if ninja --verbose install &> install.log; then \
        echo -e 'done\n'; \
      else \
        echo -e 'failed\n'; \
        cat install.log; \
      fi; \
    else \
      echo -e 'failed\n'; \
      cat ninja.log; \
    fi"

# Update environment variables
ENV LD_LIBRARY_PATH=/opt/cp2k/lib:${LD_LIBRARY} \
    PATH=/opt/cp2k/bin:${PATH}

# Collect components for installation
WORKDIR /opt/cp2k/bin
RUN /bin/bash -c -o pipefail " \
    source /opt/cp2k/tools/toolchain/install/setup; \
    mkdir -p /toolchain/install /toolchain/scripts; \
    for libdir in \$(ldd /opt/cp2k/bin/cp2k.psmp | \
                     grep /opt/cp2k/tools/toolchain/install | \
                     awk '{print \$3}' | cut -d/ -f7 | \
                     sort | uniq) setup; do \
      cp -ar /opt/cp2k/tools/toolchain/install/\${libdir} /toolchain/install; \
    done; \
    cp /opt/cp2k/tools/toolchain/scripts/tool_kit.sh /toolchain/scripts"

# Stage 2: Install CP2K
FROM ${BASE_IMAGE} AS install

# Install required packages
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    g++ \
    gcc \
    gfortran \
    python3 && rm -rf /var/lib/apt/lists/*

# Copy MPI installation
COPY --from=build /usr/local/bin /usr/bin
COPY --from=build /usr/local/include /usr/include
COPY --from=build /usr/local/lib /usr/lib
RUN ldconfig

# Install CP2K binaries
WORKDIR /opt/cp2k
COPY --from=build /opt/cp2k/bin ./bin

# Install CP2K libraries
COPY --from=build /opt/cp2k/lib ./lib

# Install CP2K database files
COPY --from=build /opt/cp2k/share ./share

# Install CP2K regression tests
COPY --from=build /opt/cp2k/tests ./tests
COPY --from=build /opt/cp2k/src/grid/sample_tasks ./src/grid/sample_tasks

# Install CP2K/Quickstep CI benchmarks
COPY --from=build /opt/cp2k/benchmarks/CI ./benchmarks/CI

# Create links to CP2K binaries
WORKDIR /opt/cp2k/bin
RUN ln -sf cp2k.psmp cp2k
RUN ln -sf cp2k.psmp cp2k.popt
RUN ln -sf cp2k.psmp cp2k_shell

# Install shared libraries required by the CP2K binaries
COPY --from=build /toolchain /opt/cp2k/tools/toolchain

# Create entrypoint script file
RUN printf "#!/bin/bash\n\
ulimit -c 0 -s unlimited\n\
\
export OMP_STACKSIZE=64M\n\
export LD_LIBRARY_PATH=/opt/cp2k/lib:\${LD_LIBRARY_PATH}\n\
export PATH=/opt/cp2k/bin:\${PATH}\n\
source /opt/cp2k/tools/toolchain/install/setup\n\
\"\$@\"" \
>/opt/cp2k/bin/entrypoint.sh && chmod 755 /opt/cp2k/bin/entrypoint.sh

# Create shortcut for regression test
RUN printf "/opt/cp2k/tests/do_regtest.py \$* /opt/cp2k/bin psmp" \
>/opt/cp2k/bin/run_tests && chmod 755 /opt/cp2k/bin/run_tests

# Define entrypoint
WORKDIR /mnt
ENTRYPOINT ["/opt/cp2k/bin/entrypoint.sh"]
CMD ["cp2k", "--help"]

# Label docker image
LABEL author="CP2K Developers" \
      cp2k_version="master" \
      container_type="production_psmp" \
      dockerfile_generator_version="0.1"

# EOF
