FROM ubuntu:18.04

# author: Ole Schuett

WORKDIR /workspace

COPY ./scripts/install_basics.sh .
RUN ./install_basics.sh

COPY ./scripts/install_formatting.sh .
RUN ./install_formatting.sh

COPY ./scripts/ci_entrypoint.sh ./scripts/test_formatting.sh ./
CMD ["./ci_entrypoint.sh", "./test_formatting.sh"]

#EOF
