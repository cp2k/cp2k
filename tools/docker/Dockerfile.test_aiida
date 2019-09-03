ARG TOOLCHAIN=cp2k/toolchain:latest
FROM ${TOOLCHAIN}

# author: Ole Schuett

WORKDIR /workspace

COPY ./scripts/install_basics.sh .
RUN ./install_basics.sh

COPY ./scripts/install_aiida.sh .
RUN ./install_aiida.sh

COPY ./scripts/ci_entrypoint.sh ./scripts/test_aiida.sh ./
CMD ["./ci_entrypoint.sh", "./test_aiida.sh"]

#EOF
