# Dockerfile for CP2K continuous integration (CI) runs
#
# A stand-alone docker build in this folder can be performed using the command:
# docker build -f build_deps_toolchain.Dockerfile ../../
#
# Author: Matthias Krack
#
# Stage 1a: Create a base image providing the dependencies for building a CP2K binary

ARG BASE_IMAGE="ubuntu:24.04"

FROM ${BASE_IMAGE} AS build_deps

# Install packages required for the CP2K toolchain build
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    autoconf \
    autogen \
    automake \
    autotools-dev \
    bzip2 \
    ca-certificates \
    g++ \
    gcc \
    gfortran \
    git \
    libtool \
    libtool-bin \
    make \
    patch \
    pkg-config \
    python3 \
    unzip \
    wget \
    zlib1g-dev

# Retrieve the number of available CPU cores
ARG NUM_PROCS
ENV NUM_PROCS=${NUM_PROCS:-32}

# Retrieve the maximum number log file lines printed on error
ARG LOG_LINES
ENV LOG_LINES=${LOG_LINES:-200}

ARG CP2K_VERSION
ENV CP2K_VERSION=${CP2K_VERSION:-psmp}
ARG CP2K_BUILD_TYPE
ENV CP2K_BUILD_TYPE=${CP2K_BUILD_TYPE:-minimal}

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

# Install libfabric version currently used by the host system
ARG LIBFABRIC_VERSION
ENV LIBFABRIC_VERSION=${LIBFABRIC_VERSION:-1.22.0}
RUN wget -q https://github.com/ofiwg/libfabric/archive/v${LIBFABRIC_VERSION}.tar.gz \
    && tar -xf v${LIBFABRIC_VERSION}.tar.gz \
    && cd libfabric-${LIBFABRIC_VERSION} \
    && ./autogen.sh \
    && ./configure --prefix=/opt/spack \
    && make -j ${NUM_PROCS} \
    && make install \
    && ldconfig \
    && cd .. \
    && rm -rf v${LIBFABRIC_VERSION}.tar.gz libfabric-${LIBFABRIC_VERSION}

# Install an MPI version ABI-compatible with the host MPI on Cray at CSCS
ARG MPICH_VERSION
ENV MPICH_VERSION=${MPICH_VERSION:-4.3.0}
RUN /bin/bash -c -o pipefail "[[ "${CP2K_VERSION}" == "psmp" ]] && (\
    wget -q https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz; \
    tar -xf mpich-${MPICH_VERSION}.tar.gz; \
    cd mpich-${MPICH_VERSION}; \
    ./configure --prefix='/usr/local' --enable-fast=all,O3 --with-device=ch4:ofi --with-libfabric=/opt/spack --with-xpmem=/opt/spack \
      FFLAGS='${FFLAGS} -fallow-argument-mismatch' \
      FCFLAGS='${FCFLAGS} -fallow-argument-mismatch' \
      &>configure.log || tail -n ${LOG_LINES} configure.log; \
    make -j ${NUM_PROCS} &>make.log || tail -n ${LOG_LINES} make.log; \
    make install &>install.log || tail -n ${LOG_LINES} install.log; \
    ldconfig; cd ..; \
    rm -rf mpich-${MPICH_VERSION} \
    rm mpich-${MPICH_VERSION}.tar.gz) || true"

# Build CP2K toolchain
COPY ./tools/toolchain /opt/cp2k/tools/toolchain
WORKDIR /opt/cp2k/tools/toolchain
RUN echo "\nCP2K build type: ${CP2K_BUILD_TYPE}\n" && \
    if [ "${CP2K_BUILD_TYPE}" = "minimal" ]; then \
      if [ "${CP2K_VERSION}" = "psmp" ]; then \
        ./install_cp2k_toolchain.sh -j${NUM_PROCS} \
          --dry-run \
          --no-arch-files \
          --target-cpu=native \
          --with-gcc=system \
          --with-mpich=system \
          --with-ninja \
          --with-cosma=no \
          --with-dftd4=no \
          --with-elpa=no \
          --with-fftw=no \
          --with-libgrpp=no \
          --with-libint=no \
          --with-libvori=no \
          --with-libxc=no \
          --with-libxsmm=no \
          --with-sirius=no \
          --with-spglib=no; \
      elif [ "${CP2K_VERSION}" = "ssmp" ]; then \
        ./install_cp2k_toolchain.sh \
          --dry-run \
          --mpi-mode=no \
          --no-arch-files \
          --target-cpu=native \
          --with-dbcsr \
          --with-gcc=system \
          --with-ninja \
          --with-dftd4=no \
          --with-fftw=no \
          --with-libgrpp=no \
          --with-libint=no \
          --with-libvori=no \
          --with-libxc=no \
          --with-libxsmm=no \
          --with-spglib=no; \
      fi; \
    elif [ "${CP2K_BUILD_TYPE}" = "all" ]; then \
      if [ "${CP2K_VERSION}" = "psmp" ]; then \
        ./install_cp2k_toolchain.sh -j${NUM_PROCS} \
          --dry-run \
          --install-all \
          --no-arch-files \
          --target-cpu=native \
          --with-gcc=system \
          --with-mpich=system; \
      elif [ "${CP2K_VERSION}" = "ssmp" ]; then \
        ./install_cp2k_toolchain.sh \
          --dry-run \
          --install-all \
          --mpi-mode=no \
          --no-arch-files \
          --target-cpu=native \
          --with-gcc=system; \
      fi; \
    else \
      echo "ERROR: Unknown CP2K_BUILD_TYPE ${CP2K_BUILD_TYPE} was specified"; \
      exit 1; \
    fi

# Perform toolchain build step-wise in stages after its initialization with dry-run
COPY ./tools/toolchain/scripts/stage0/ ./scripts/stage0/
RUN  ./scripts/stage0/install_stage0.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage1/ ./scripts/stage1/
RUN  ./scripts/stage1/install_stage1.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage2/ ./scripts/stage2/
RUN  ./scripts/stage2/install_stage2.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage3/ ./scripts/stage3/
RUN  ./scripts/stage3/install_stage3.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage4/ ./scripts/stage4/
RUN  ./scripts/stage4/install_stage4.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage5/ ./scripts/stage5/
RUN  ./scripts/stage5/install_stage5.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage6/ ./scripts/stage6/
RUN  ./scripts/stage6/install_stage6.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage7/ ./scripts/stage7/
RUN  ./scripts/stage7/install_stage7.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage8/ ./scripts/stage8/
RUN  ./scripts/stage8/install_stage8.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage9/ ./scripts/stage9/
RUN  ./scripts/stage9/install_stage9.sh && rm -rf ./build

# Stage 1b: Install CP2K toolchain
FROM ${BASE_IMAGE} AS install_deps

# Install required packages
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    g++ \
    gcc \
    gfortran \
    libtool \
    libtool-bin \
    pkg-config \
    python3 && rm -rf /var/lib/apt/lists/*

ARG CP2K_VERSION
ENV CP2K_VERSION=${CP2K_VERSION}
ARG CP2K_BUILD_TYPE
ENV CP2K_BUILD_TYPE=${CP2K_BUILD_TYPE}

# Copy MPI installation
COPY --from=build_deps /usr/local/bin /usr/local/bin
COPY --from=build_deps /usr/local/include /usr/local/include
COPY --from=build_deps /usr/local/lib /usr/local/lib
RUN ldconfig

# Install CP2K dependencies (toolchain)
COPY --from=build_deps /opt/cp2k/tools/toolchain/install /opt/cp2k/tools/toolchain/install
COPY --from=build_deps /opt/cp2k/tools/toolchain/scripts /opt/cp2k/tools/toolchain/scripts

# EOF
