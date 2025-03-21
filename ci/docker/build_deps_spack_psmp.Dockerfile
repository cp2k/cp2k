# Dockerfile for CP2K continous integration (CI) runs
#
# A stand-alone docker build in this folder can be performed using the command:
# docker build -f build_deps_spack_psmp.Dockerfile ../../
#
# Author: Matthias Krack
#
# Stage 1: Create a base image providing the dependencies for a CP2K spack_psmp build

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
    unzip \
    wget \
    xxd \
    xz-utils \
    zstd

# Retrieve the number of available CPU cores
ARG NUM_PROCS
ENV NUM_PROCS=${NUM_PROCS:-32}

# Install the latest Spack development version
WORKDIR /root
RUN git clone -c feature.manyFiles=true --depth=2 https://github.com/spack/spack.git
ENV PATH="/root/spack/bin:${PATH}"

# Find all compilers
RUN spack compiler find

# Find all external packages
RUN spack external find --all --not-buildable

# Enable Spack build cache from the latest development version
ARG SPACK_BUILD_CACHE
ENV SPACK_BUILD_CACHE="${SPACK_BUILD_CACHE:-develop}"
RUN spack mirror add ${SPACK_BUILD_CACHE} https://binaries.spack.io/${SPACK_BUILD_CACHE} && \
    spack buildcache keys --install --trust --force && \
    spack mirror remove ${SPACK_BUILD_CACHE}

# Install CP2K dependencies via Spack
ARG CP2K_BUILD_TYPE
ENV CP2K_BUILD_TYPE=${CP2K_BUILD_TYPE:-minimal}
COPY ./ci/spack/cp2k_deps_${CP2K_BUILD_TYPE}.yaml .
RUN spack env create myenv ./cp2k_deps_${CP2K_BUILD_TYPE}.yaml
RUN spack -e myenv concretize -f
ENV SPACK_ENV_VIEW="/root/spack/var/spack/environments/myenv/spack-env/view"
RUN spack -e myenv env depfile -o spack_makefile && \
    make -j${NUM_PROCS} --file=spack_makefile SPACK_COLOR=never --output-sync=recurse && \
    cp -ar ${SPACK_ENV_VIEW}/bin ${SPACK_ENV_VIEW}/include ${SPACK_ENV_VIEW}/lib /opt/spack

# Label docker image
LABEL author="CP2K Developers" \
      cp2k_version="master" \
      image_type="spack_psmp_base_image"

# EOF
