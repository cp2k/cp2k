FROM ubuntu:zesty
USER root

RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    gfortran \
    python \
    wget \
    --no-install-recommends \
  && rm -rf /var/lib/apt/lists/*

# --install-all would additionally require: bison bisonc++ flex flexc++

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

WORKDIR /opt/cp2k-toolchain/
COPY install_cp2k_toolchain.sh ./
COPY scripts ./scripts/
RUN ./install_cp2k_toolchain.sh && rm -rf ./build

#EOF