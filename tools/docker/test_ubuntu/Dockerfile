FROM ubuntu:rolling

# author: Ole Schuett

# install Ubuntu packages
RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    gfortran \
    python \
    git \
    wget \
    less \
    unzip \
    fftw3-dev \
    libopenblas-dev \
    liblapack-dev \
    libint-dev \
    rsync \
    --no-install-recommends \
  && rm -rf /var/lib/apt/lists/*

RUN ln -sf bash /bin/sh

# download and compile cp2k snapshot
WORKDIR /opt/
RUN wget -q -O cp2k-master.zip https://github.com/cp2k/cp2k/archive/master.zip && \
    unzip -q cp2k-master.zip  && \
    rm cp2k-master.zip

# build toolchain relying mostly on ubuntu packages
WORKDIR /opt/cp2k-master/cp2k/tools/toolchain
RUN ./install_cp2k_toolchain.sh \
    --mpi-mode=no \
    --with-gcc=system \
    --with-binutils=system \
    --with-make=system \
    --with-fftw=system \
    --with-openblas=system \
    --with-reflapack=system \
    --with-libint=system \
    --with-libxc=install \
    --with-libxsmm=install \
  && rm -rf ./build

WORKDIR /opt/cp2k-master/cp2k/arch
RUN ln -vs ../tools/toolchain/install/arch/local* .

# run regtests which lack fixed reference value
WORKDIR /opt/cp2k-master/cp2k/makefiles
RUN source ../tools/toolchain/install/setup && \
    make -j   VERSION=ssmp && \
    make test VERSION=ssmp TESTOPTS="-restrictdir QS/regtest-almo-md -restrictdir QS/regtest-almo-1 -restrictdir SE/regtest-3-4 -restrictdir QS/regtest-ot-1-vib -restrictdir Fist/regtest-5-vib -restrictdir QS/regtest-optbas -restrictdir TMC/regtest_ana_post_proc" && \
    rm -rf ../lib/ ../exe/ ../regtesting/local/ssmp/TEST-*

WORKDIR /opt/cp2k_test_ubuntu
COPY ./cmd.sh .
CMD ["./cmd.sh"]

#EOF
