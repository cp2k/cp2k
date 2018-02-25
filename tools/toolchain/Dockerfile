FROM ubuntu:latest
USER root

# install Ubuntu packages
RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    gfortran \
    python \
    wget \
    bison \
    bisonc++ \
    flex \
    flexc++ \
    zlib1g-dev \
    less \
    unzip \
    libc6-dbg \
    --no-install-recommends \
  && rm -rf /var/lib/apt/lists/*

# Shifter at CSCS requires MPICH v3.1.4 for ABI-compatibility
# https://user.cscs.ch/scientific_computing/supported_applications/shifter/advanced_shifter/#native-mpi-support
RUN wget -q http://www.mpich.org/static/downloads/3.1.4/mpich-3.1.4.tar.gz \
    && tar xf mpich-3.1.4.tar.gz \
    && cd mpich-3.1.4 \
    && ./configure --enable-fast=all,O3 --prefix=/usr \
    && make -j \
    && make install \
    && ldconfig \
    && cd .. \
    && rm -rf mpich-3.1.4 \
    && rm mpich-3.1.4.tar.gz

# build toolchain
WORKDIR /opt/cp2k-toolchain/
COPY install_cp2k_toolchain.sh ./
COPY scripts ./scripts/
RUN ./install_cp2k_toolchain.sh --install-all  &&  rm -rf ./build

# configure shell
RUN ln -sf bash /bin/sh && \
    echo "source /opt/cp2k-toolchain/install/setup" >> /etc/bash.bashrc

WORKDIR /opt/cp2k-toolchain/

#EOF
