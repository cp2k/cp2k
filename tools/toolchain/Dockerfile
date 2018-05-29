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

# build toolchain
WORKDIR /opt/cp2k-toolchain/
COPY install_cp2k_toolchain.sh ./
COPY scripts ./scripts/

# Do not compile make until glibc issue is resolved:
#  http://gnu-make.2324884.n4.nabble.com/undefined-reference-to-alloca-td18308.html
#  http://git.savannah.gnu.org/cgit/make.git/commit/?id=48c8a116a914a325a0497721f5d8b58d5bba34d4
RUN ./install_cp2k_toolchain.sh --install-all --with-make=no && rm -rf ./build

# configure shell
RUN ln -sf bash /bin/sh && \
    echo "source /opt/cp2k-toolchain/install/setup" >> /etc/bash.bashrc

WORKDIR /opt/cp2k-toolchain/

#EOF
