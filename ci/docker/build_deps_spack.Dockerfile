# Dockerfile for CP2K continuous integration (CI) runs
#
# A stand-alone docker build in this folder can be performed using the command:
# docker build -f build_deps_spack.Dockerfile ../../
#
# Author: Matthias Krack
#
# Stage 1: Create a base image providing the dependencies for building a CP2K binary

ARG BASE_IMAGE="ubuntu:24.04"

FROM ${BASE_IMAGE} AS build_deps

# Install packages required to build the CP2K dependencies with Spack
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    bzip2 \
    ca-certificates \
    cmake \
    g++ \
    gcc \
    gfortran \
    git \
    gnupg \
    hwloc \
    libhwloc-dev \
    libssh-dev \
    libssl-dev \
    libtool \
    libtool-bin \
    lsb-release \
    make \
    ninja-build \
    patch \
    pkgconf \
    python3 \
    python3-dev \
    unzip \
    wget \
    xxd \
    xz-utils \
    zstd && rm -rf /var/lib/apt/lists/*

# Retrieve the number of available CPU cores
ARG NUM_PROCS
ENV NUM_PROCS=${NUM_PROCS:-32}

# Install a recent Spack version
WORKDIR /root/spack
ARG SPACK_VERSION
ENV SPACK_VERSION=${SPACK_VERSION:-develop-2025-05-18}
RUN git init --quiet && \
    git remote add origin https://github.com/spack/spack.git && \
    git fetch --quiet --depth 1 origin ${SPACK_VERSION} --no-tags && \
    git checkout --quiet FETCH_HEAD
ENV PATH="/root/spack/bin:${PATH}"

# Find all compilers
RUN spack compiler find

# Find all external packages
RUN spack external find --all --not-buildable

# Enable Spack build cache
ARG SPACK_BUILD_CACHE
ENV SPACK_BUILD_CACHE="${SPACK_BUILD_CACHE:-develop-2025-05-18}"
RUN spack mirror add ${SPACK_BUILD_CACHE} https://binaries.spack.io/${SPACK_BUILD_CACHE} && \
    spack buildcache keys --install --trust --force && \
    spack mirror remove ${SPACK_BUILD_CACHE}

# Copy Spack configuration and build recipes
ARG CP2K_VERSION
ENV CP2K_VERSION=${CP2K_VERSION:-ssmp}
ARG CP2K_BUILD_TYPE
ENV CP2K_BUILD_TYPE=${CP2K_BUILD_TYPE:-minimal}
COPY ./tools/spack/repo.yaml ./tools/spack/cp2k_deps_${CP2K_BUILD_TYPE}_${CP2K_VERSION}.yaml ./
COPY ./tools/spack/packages ./packages

# Sarus containers must be dynamically linked to an MPI implementation that is ABI-compatible
# with the MPI on the compute nodes at CSCS like MPICH@3
ARG MPICH_VERSION
ENV MPICH_VERSION=${MPICH_VERSION:-3.4.3}
RUN sed -i -e "s/mpich@[0-9.]*/mpich@${MPICH_VERSION}/" cp2k_deps_${CP2K_BUILD_TYPE}_${CP2K_VERSION}.yaml
RUN spack env create myenv cp2k_deps_${CP2K_BUILD_TYPE}_${CP2K_VERSION}.yaml && \
    spack -e myenv repo list

# Install CP2K dependencies via Spack
RUN spack -e myenv concretize -f
ENV SPACK_ENV_VIEW="/root/spack/var/spack/environments/myenv/spack-env/view"
RUN spack -e myenv env depfile -o spack_makefile && \
    make -j${NUM_PROCS} --file=spack_makefile SPACK_COLOR=never --output-sync=recurse && \
    cp -ar ${SPACK_ENV_VIEW}/bin ${SPACK_ENV_VIEW}/include ${SPACK_ENV_VIEW}/lib /opt/spack

# EOF
