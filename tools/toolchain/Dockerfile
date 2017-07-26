FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
    build-essential \
    mpich \
    python \
    wget \
  && rm -rf /var/lib/apt/lists/*

# --install-all would additionally require: bison bisonc++ flex flexc++

WORKDIR /opt/cp2k-toolchain/
COPY install_cp2k_toolchain.sh ./
COPY scripts ./scripts/
RUN ./install_cp2k_toolchain.sh

#-------------------------------------------------------------------------------
# leave layers with build files behind to reduce image size
FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
    build-essential \
    mpich \
    python \
    wget \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/cp2k-toolchain/
COPY --from=0 /opt/cp2k-toolchain/install ./install/

#EOF