ARG TOOLCHAIN=cp2k/toolchain:latest
FROM ${TOOLCHAIN}

# author: Ole Schuett

WORKDIR /workspace

COPY ./scripts/install_basics.sh .
RUN ./install_basics.sh

COPY ./scripts/install_i-pi.sh .
RUN ./install_i-pi.sh

COPY ./scripts/ci_entrypoint.sh ./scripts/test_i-pi.sh ./
CMD ["./ci_entrypoint.sh", "./test_i-pi.sh"]

#EOF
