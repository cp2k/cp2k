FROM ubuntu:24.04

RUN apt-get update -qq
RUN apt-get install -qq --no-install-recommends \
    ca-certificates \
    cmake \
    gcc \
    g++ \
    libgmp-dev \
    libboost-dev \
    make \
    wget \
    xz-utils

WORKDIR /opt
COPY ./build_bundles.sh .
RUN ./build_bundles.sh

#EOF
