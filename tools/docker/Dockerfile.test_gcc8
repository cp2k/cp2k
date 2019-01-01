FROM ubuntu:18.04

# author: Ole Schuett

WORKDIR /workspace

COPY ./scripts/install_basics.sh .
RUN ./install_basics.sh

COPY ./scripts/install_ubuntu_toolchain.sh .
RUN ./install_ubuntu_toolchain.sh 8

COPY ./scripts/install_regtest.sh .
RUN ./install_regtest.sh local ssmp

COPY ./scripts/ci_entrypoint.sh ./scripts/test_regtest.sh ./
CMD ["./ci_entrypoint.sh", "./test_regtest.sh", "local", "ssmp"]

#EOF
