ARG BASE_IMAGE=ubuntu:20.04
FROM ${BASE_IMAGE}

ARG MPI_MODE=mpich
ARG LIBINT_LMAX=5
ARG BASE_IMAGE

# author: Ole Schuett

# Installs complete toolchain including packages which are skipped by default.

USER root

# install Linux packages.
COPY ./install_requirements*.sh ./
RUN ./install_requirements.sh ${BASE_IMAGE}

# copy helper scripts
WORKDIR /opt/cp2k-toolchain
RUN mkdir scripts
COPY ./scripts/VERSION \
     ./scripts/parse_if.py \
     ./scripts/tool_kit.sh \
     ./scripts/common_vars.sh \
     ./scripts/signal_trap.sh \
     ./scripts/get_openblas_arch.sh \
     ./scripts/

# Dry-run leaves behind config files for the followup install scripts.
# This breaks up the lengthy installation into smaller docker build steps.
COPY ./install_cp2k_toolchain.sh .
RUN ./install_cp2k_toolchain.sh  \
    --install-all \
    --mpi-mode=${MPI_MODE} \
    --libint-lmax=${LIBINT_LMAX} \
    --dry-run

COPY ./scripts/stage0/ ./scripts/stage0/
RUN  ./scripts/stage0/install_stage0.sh && rm -rf ./build

COPY ./scripts/stage1/ ./scripts/stage1/
RUN  ./scripts/stage1/install_stage1.sh && rm -rf ./build

COPY ./scripts/stage2/ ./scripts/stage2/
RUN  ./scripts/stage2/install_stage2.sh && rm -rf ./build

COPY ./scripts/stage3/ ./scripts/stage3/
RUN  ./scripts/stage3/install_stage3.sh && rm -rf ./build

COPY ./scripts/stage4/ ./scripts/stage4/
RUN  ./scripts/stage4/install_stage4.sh && rm -rf ./build

COPY ./scripts/stage5/ ./scripts/stage5/
RUN  ./scripts/stage5/install_stage5.sh && rm -rf ./build

COPY ./scripts/stage6/ ./scripts/stage6/
RUN  ./scripts/stage6/install_stage6.sh && rm -rf ./build

COPY ./scripts/stage7/ ./scripts/stage7/
RUN  ./scripts/stage7/install_stage7.sh && rm -rf ./build

COPY ./scripts/stage8/ ./scripts/stage8/
RUN  ./scripts/stage8/install_stage8.sh && rm -rf ./build

COPY ./scripts/arch_base.tmpl \
     ./scripts/generate_arch_files.sh \
     ./scripts/
RUN ./scripts/generate_arch_files.sh && rm -rf ./build

#EOF
