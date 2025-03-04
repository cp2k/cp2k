# Dockerfile for CP2K continous integration (CI) runs
#
# A stand-alone docker build in this folder can be performed using the command:
# docker build -f build_deps_toolchain_psmp.Dockerfile ../../
#
# Author: Matthias Krack
#
# Stage 1a: Create a base image providing the dependencies for a CP2K toolchain_psmp build

ARG BASE_IMAGE="ubuntu:24.04"

FROM ${BASE_IMAGE} AS build_deps

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
ARG NUM_PROCS=32

# Retrieve the maximum number log file lines printed on error
ARG LOG_LINES=100

# Install an MPI version ABI-compatible with the host MPI on Cray at CSCS
ARG MPICH_VERSION=3.4.3
RUN /bin/bash -c -o pipefail " \
    wget -q https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz; \
    tar -xf mpich-${MPICH_VERSION}.tar.gz; \
    cd mpich-${MPICH_VERSION}; \
    ./configure --prefix='/usr/local' --enable-fast=all,O3 --with-device=ch3 \
      FFLAGS='${FFLAGS} -fallow-argument-mismatch' \
      FCFLAGS='${FCFLAGS} -fallow-argument-mismatch' \
      &>configure.log || tail -n ${LOG_LINES} configure.log; \
    make -j ${NUM_PROCS} &>make.log || tail -n ${LOG_LINES} make.log; \
    make install &>install.log || tail -n ${LOG_LINES} install.log; \
    ldconfig; cd ..; \
    rm -rf mpich-${MPICH_VERSION} \
    rm mpich-${MPICH_VERSION}.tar.gz"

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

# Stage 1b: Install CP2K toolchain
FROM ${BASE_IMAGE} AS install_deps

# Install required packages
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    g++ \
    gcc \
    gfortran \
    libtool \
    libtool-bin \
    pkg-config \
    python3 && rm -rf /var/lib/apt/lists/*

# Copy MPI installation
COPY --from=build_deps /usr/local/bin /usr/local/bin
COPY --from=build_deps /usr/local/include /usr/local/include
COPY --from=build_deps /usr/local/lib /usr/local/lib

# Install CP2K dependencies (toolchain)
COPY --from=build_deps /opt/cp2k/tools/toolchain/install /opt/cp2k/tools/toolchain/install
COPY --from=build_deps /opt/cp2k/tools/toolchain/scripts /opt/cp2k/tools/toolchain/scripts

# Label docker image
LABEL author="CP2K Developers" \
      cp2k_version="master" \
      image_type="toolchain_psmp_base_image"

# EOF
