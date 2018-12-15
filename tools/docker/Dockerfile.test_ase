ARG TOOLCHAIN=cp2k/toolchain:latest
FROM ${TOOLCHAIN}

# author: Ole Schuett

WORKDIR /workspace

COPY ./scripts/install_basics.sh .
RUN ./install_basics.sh

COPY ./scripts/install_ase.sh .
RUN ./install_ase.sh

COPY ./scripts/ci_entrypoint.sh ./scripts/test_ase.sh ./
CMD ["./ci_entrypoint.sh", "./test_ase.sh"]

#EOF
