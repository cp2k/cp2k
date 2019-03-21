FROM ubuntu:18.04
USER root

# install Ubuntu packages
# Include also packages for the CI to make check_runs similar to full toolchain builds.
RUN apt-get update && apt-get install -y --no-install-recommends \
    autoconf \
    autogen \
    automake \
    autotools-dev \
    bison \
    bisonc++ \
    build-essential \
    ca-certificates \
    cmake \
    flex \
    flexc++ \
    gfortran \
    git \
    less \
    libc6-dbg \
    libtool \
    make \
    nano \
    pkg-config \
    python \
    rsync \
    unzip \
    wget \
    zlib1g-dev \
   && rm -rf /var/lib/apt/lists/*

# build toolchain
WORKDIR /opt/cp2k-toolchain
COPY . /opt/cp2k-toolchain
RUN ./install_cp2k_toolchain.sh --install-all --with-binutils=system --with-cmake=system --with-make=system \
    && rm -rf ./build

# configure shell
RUN ln -sf bash /bin/sh && \
    echo "source /opt/cp2k-toolchain/install/setup" >> /etc/bash.bashrc

WORKDIR /opt/cp2k-toolchain

#EOF
