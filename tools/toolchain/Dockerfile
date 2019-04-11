FROM ubuntu:18.04
USER root

# install Ubuntu packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    autoconf \
    autogen \
    automake \
    autotools-dev \
    ca-certificates \
    g++ \
    git \
    less \
    libtool \
    make \
    nano \
    pkg-config \
    python \
    unzip \
    wget \
    zlib1g-dev \
   && rm -rf /var/lib/apt/lists/*

# copy helper scripts
WORKDIR /opt/cp2k-toolchain
RUN mkdir scripts
COPY ./scripts/VERSION \
     ./scripts/parse_if.py \
     ./scripts/tool_kit.sh \
     ./scripts/common_vars.sh \
     ./scripts/signal_trap.sh ./scripts/

COPY ./install_cp2k_toolchain.sh .
RUN ./install_cp2k_toolchain.sh --install-all --dry-run

COPY ./scripts/install_valgrind.sh ./scripts/
RUN ./scripts/install_valgrind.sh && rm -rf ./build

COPY ./scripts/install_cmake.sh ./scripts/
RUN ./scripts/install_cmake.sh && rm -rf ./build

COPY ./scripts/install_gcc.sh ./scripts/
RUN ./scripts/install_gcc.sh && rm -rf ./build

COPY ./scripts/get_openblas_arch.sh ./scripts/
COPY ./scripts/setup_buildtools.sh ./scripts/
RUN ./scripts/setup_buildtools.sh && rm -rf ./build

COPY ./scripts/install_mpich.sh ./scripts/
RUN ./scripts/install_mpich.sh && rm -rf ./build

COPY ./scripts/install_openmpi.sh ./scripts/
RUN ./scripts/install_openmpi.sh && rm -rf ./build

COPY ./scripts/install_reflapack.sh \
     ./scripts/install_mkl.sh \
     ./scripts/install_acml.sh \
     ./scripts/install_openblas.sh \
     ./scripts/install_mathlibs.sh ./scripts/
RUN ./scripts/install_mathlibs.sh && rm -rf ./build

COPY ./scripts/install_fftw.sh ./scripts/
RUN ./scripts/install_fftw.sh && rm -rf ./build

COPY ./scripts/install_libint.sh ./scripts/
RUN ./scripts/install_libint.sh && rm -rf ./build

COPY ./scripts/install_libxc.sh ./scripts/
RUN ./scripts/install_libxc.sh && rm -rf ./build

COPY ./scripts/install_libsmm.sh ./scripts/
RUN ./scripts/install_libsmm.sh && rm -rf ./build

COPY ./scripts/install_libxsmm.sh ./scripts/
RUN ./scripts/install_libxsmm.sh && rm -rf ./build

COPY ./scripts/install_scalapack.sh ./scripts/
RUN ./scripts/install_scalapack.sh && rm -rf ./build

COPY ./scripts/install_elpa.sh ./scripts/
RUN ./scripts/install_elpa.sh && rm -rf ./build

COPY ./scripts/install_ptscotch.sh ./scripts/
RUN ./scripts/install_ptscotch.sh && rm -rf ./build

COPY ./scripts/install_parmetis.sh ./scripts/
RUN ./scripts/install_parmetis.sh && rm -rf ./build

COPY ./scripts/install_metis.sh ./scripts/
RUN ./scripts/install_metis.sh && rm -rf ./build

COPY ./scripts/install_superlu.sh ./scripts/
RUN ./scripts/install_superlu.sh && rm -rf ./build

COPY ./scripts/install_pexsi.sh ./scripts/
RUN ./scripts/install_pexsi.sh && rm -rf ./build

COPY ./scripts/install_quip.sh ./scripts/
RUN ./scripts/install_quip.sh && rm -rf ./build

COPY ./scripts/install_gsl.sh ./scripts/
RUN ./scripts/install_gsl.sh && rm -rf ./build

COPY ./scripts/install_spglib.sh ./scripts/
RUN ./scripts/install_spglib.sh && rm -rf ./build

COPY ./scripts/install_hdf5.sh ./scripts/
RUN ./scripts/install_hdf5.sh && rm -rf ./build

COPY ./scripts/install_libvdwxc.sh ./scripts/
RUN ./scripts/install_libvdwxc.sh && rm -rf ./build

COPY ./scripts/install_sirius.sh ./scripts/
RUN ./scripts/install_sirius.sh && rm -rf ./build

COPY ./scripts/install_json_fortran.sh ./scripts/
RUN ./scripts/install_json_fortran.sh && rm -rf ./build

COPY ./scripts/arch_base.tmpl \
     ./scripts/generate_arch_files.sh ./scripts/
RUN ./scripts/generate_arch_files.sh && rm -rf ./build

#EOF
