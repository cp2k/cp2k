ARG TOOLCHAIN=cp2k/toolchain:latest
FROM ${TOOLCHAIN}

# author: Ole Schuett

WORKDIR /workspace

COPY ./scripts/install_basics.sh .
RUN ./install_basics.sh

COPY ./scripts/install_manual.sh .
RUN ./install_manual.sh

COPY ./scripts/ci_entrypoint.sh ./scripts/test_manual.sh ./
CMD ["./ci_entrypoint.sh", "./test_manual.sh"]

#EOF
