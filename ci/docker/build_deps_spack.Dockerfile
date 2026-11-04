# Dockerfile for CP2K continuous integration (CI) runs
#
# A stand-alone build in this folder can be performed with:
# podman build --shm-size=1g -f build_deps_spack.Dockerfile ../../
#
# Stage 1: Create a base image providing the dependencies for building a CP2K binary

ARG BASE_IMAGE="docker.io/ubuntu:26.04"

FROM "${BASE_IMAGE}" AS build_deps

ARG APT_GET="apt-get -o APT::Sandbox::User=root -qq"

# Install packages required to build the CP2K dependencies with Spack
RUN ${APT_GET} update && \
    ${APT_GET} install --no-install-recommends \
    bzip2 \
    ca-certificates \
    cmake \
    g++ gcc gfortran \
    git \
    gnupg \
    libssh-dev \
    libssl-dev \
    libtool \
    libtool-bin \
    lsb-release \
    make \
    patch \
    pkgconf \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    unzip \
    wget \
    xxd \
    xz-utils \
    zstd \
    && rm -rf /var/lib/apt/lists/*

# Setup CUDA environment
ENV CUDA_HOME="/usr/local/cuda"
ENV LD_LIBRARY_PATH="${CUDA_HOME}/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

# Number of parallel build jobs (override with --build-arg NUM_PROCS=...)
ARG NUM_PROCS="32"

ARG CP2K_VERSION="psmp"
ENV CP2K_VERSION="${CP2K_VERSION}"

ARG FEATURE_FLAGS=""
ENV FEATURE_FLAGS="${FEATURE_FLAGS}"

# Build CP2K dependencies
WORKDIR /opt/cp2k
COPY ./make_cp2k.sh .
COPY ./tools/docker ./tools/docker
COPY ./tools/spack ./tools/spack

RUN ./make_cp2k.sh -bd_only -cray -cv ${CP2K_VERSION} -j${NUM_PROCS} -uc no -ue ${FEATURE_FLAGS}
