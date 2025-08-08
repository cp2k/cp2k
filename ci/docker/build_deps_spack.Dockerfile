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

# Create dummy xpmem library for the MPICH build. At runtime the
# container engine will inject the xpmem library from the host system
RUN git clone https://github.com/hpc/xpmem \
    && cd xpmem/lib \
    && gcc -I../include -shared -o libxpmem.so.1 libxpmem.c \
    && ln -s libxpmem.so.1 libxpmem.so \
    && mkdir -p /opt/spack/lib /opt/spack/include \
    && mv libxpmem.so* /opt/spack/lib \
    && cp ../include/xpmem.h /opt/spack/include/ \
    && ldconfig /opt/spack/lib \
    && cd ../../ \
    && rm -rf xpmem

# Retrieve the number of available CPU cores
ARG NUM_PROCS
ENV NUM_PROCS=${NUM_PROCS:-32}

# Install Spack and Spack packages
WORKDIR /root/spack
ARG SPACK_VERSION
ENV SPACK_VERSION=${SPACK_VERSION:-1.0.0}
ARG SPACK_PACKAGES_VERSION
ENV SPACK_PACKAGES_VERSION=${SPACK_PACKAGES_VERSION:-2025.07.0}
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

ARG USE_ALPS_SPACK_REPO
ENV USE_ALPS_SPACK_REPO=${USE_ALPS_SPACK_REPO:-0}
ENV ALPS_SPACK_REPO_SHA="57d5c2f39e4f0cf6acc8779dcaedc1e90264d639"
RUN if [ "${USE_ALPS_SPACK_REPO}" -ne 0 ]; then \
    mkdir -p /root/alps-cluster-config && \
    wget -qO- "https://api.github.com/repos/eth-cscs/alps-cluster-config/tarball/$ALPS_SPACK_REPO_SHA" | \
    tar --strip-components=1 -xz -C /root/alps-cluster-config && \
    spack repo add --scope site /root/alps-cluster-config/site/spack_repo/alps; \
    fi

# Find all compilers
RUN spack compiler find

# Find all external packages
RUN spack external find --all --not-buildable

# Copy Spack configuration and build recipes
ARG CP2K_VERSION
ENV CP2K_VERSION=${CP2K_VERSION:-psmp}
COPY ./tools/spack/cp2k_deps_${CP2K_VERSION}.yaml ./
COPY ./tools/spack/cp2k_dev_repo ${SPACK_PACKAGES_ROOT}/repos/spack_repo/cp2k_dev_repo/
RUN spack repo add --scope site ${SPACK_PACKAGES_ROOT}/repos/spack_repo/cp2k_dev_repo/
RUN spack env create myenv cp2k_deps_${CP2K_VERSION}.yaml && \
    spack -e myenv repo list

# Use Cray-MPICH if the ALPS Spack repository is used and add supported dependencies
# NOTE: This is a workaround for ALPS until the CE provides full MPI replacement
# FIXME: SIRIUS with libvdwxc
RUN if [ "${USE_ALPS_SPACK_REPO}" -ne 0 ]; then \
    spack -e myenv config remove "packages:mpi" && \
    spack -e myenv config add "packages:mpi:require:cray-mpich" && \
    spack -e myenv remove "mpich"; \
    spack -e myenv config remove "packages:sirius" && \
    spack -e myenv config add "packages:sirius:require:+fortran+pugixml~apps~vdwxc" && \
    spack -e myenv add "cray-mpich"; \
    spack -e myenv add "dla-future-fortran"; \
    fi

# Install CP2K dependencies via Spack
RUN spack -e myenv config get
RUN spack -e myenv concretize -f
ENV SPACK_ENV_VIEW="${SPACK_ROOT}/var/spack/environments/myenv/spack-env/view"
RUN spack -e myenv env depfile -o spack_makefile && \
    make -j${NUM_PROCS} --file=spack_makefile SPACK_COLOR=never --output-sync=recurse && \
    cp -ar ${SPACK_ENV_VIEW}/bin ${SPACK_ENV_VIEW}/include ${SPACK_ENV_VIEW}/lib /opt/spack

# EOF
