ARG TOOLCHAIN=cp2k/toolchain:latest
FROM ${TOOLCHAIN}

# author: Ole Schuett

WORKDIR /workspace

COPY ./scripts/install_basics.sh .
RUN ./install_basics.sh

COPY ./scripts/install_regtest.sh .
RUN ./install_regtest.sh local_valgrind sopt

COPY ./scripts/ci_entrypoint.sh ./scripts/test_regtest.sh ./
CMD ["./ci_entrypoint.sh", "./test_regtest.sh", "local_valgrind", "sopt"]

#EOF
