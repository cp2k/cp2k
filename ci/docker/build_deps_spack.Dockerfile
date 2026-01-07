# Dockerfile for CP2K continuous integration (CI) runs
#
# A stand-alone build in this folder can be performed with:
# podman build --shm-size=1g -f build_deps_spack.Dockerfile ../../
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

# Install Spack and Spack packages
WORKDIR /root/spack
ARG SPACK_VERSION
ENV SPACK_VERSION=${SPACK_VERSION:-1.1.0}
ARG SPACK_PACKAGES_VERSION
ENV SPACK_PACKAGES_VERSION=${SPACK_PACKAGES_VERSION:-2025.11.0}
ARG SPACK_REPO=https://github.com/spack/spack
ENV SPACK_ROOT=/opt/spack-${SPACK_VERSION}
ARG SPACK_PACKAGES_REPO=https://github.com/spack/spack-packages
ENV SPACK_PACKAGES_ROOT=/opt/spack-packages-${SPACK_PACKAGES_VERSION}
RUN mkdir -p ${SPACK_ROOT} \
    && wget -q ${SPACK_REPO}/archive/v${SPACK_VERSION}.tar.gz \
    && tar -xzf v${SPACK_VERSION}.tar.gz -C /opt && rm -f v${SPACK_VERSION}.tar.gz \
    && mkdir -p ${SPACK_PACKAGES_ROOT} \
    && wget -q ${SPACK_PACKAGES_REPO}/archive/v${SPACK_PACKAGES_VERSION}.tar.gz \
    && tar -xzf v${SPACK_PACKAGES_VERSION}.tar.gz -C /opt && rm -f v${SPACK_PACKAGES_VERSION}.tar.gz

ENV PATH="${SPACK_ROOT}/bin:${PATH}"

# Add Spack packages builtin repository
RUN spack repo add --scope site ${SPACK_PACKAGES_ROOT}/repos/spack_repo/builtin

# Find all compilers
RUN spack compiler find

# Find all external packages
RUN spack external find --all --not-buildable

# Copy Spack configuration and build recipes
ARG CP2K_VERSION
ENV CP2K_VERSION=${CP2K_VERSION:-psmp}
COPY ./tools/spack/cp2k_deps_${CP2K_VERSION}.yaml ./
RUN sed -i -e "s/~xpmem/+xpmem/" cp2k_deps_${CP2K_VERSION}.yaml
COPY ./tools/spack/spack_repo/cp2k_dev ${SPACK_PACKAGES_ROOT}/repos/spack_repo/cp2k_dev/
RUN spack repo add --scope site ${SPACK_PACKAGES_ROOT}/repos/spack_repo/cp2k_dev/
RUN spack env create myenv cp2k_deps_${CP2K_VERSION}.yaml && \
    spack -e myenv repo list

# Install CP2K dependencies via Spack
RUN spack -e myenv concretize -f
RUN spack -e myenv env depfile -o spack_makefile
RUN make -j${NUM_PROCS} --file=spack_makefile SPACK_COLOR=never --output-sync=recurse

# EOF
