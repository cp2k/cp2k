ARG TOOLCHAIN=cp2k/toolchain:latest
FROM ${TOOLCHAIN}

# author: Ole Schuett

WORKDIR /workspace

COPY ./scripts/install_basics.sh .
RUN ./install_basics.sh

COPY ./scripts/install_scaling.sh .
RUN ./install_scaling.sh

COPY ./scripts/ci_entrypoint.sh ./scripts/test_scaling.sh ./
CMD ["./ci_entrypoint.sh", "./test_scaling.sh"]

#EOF
