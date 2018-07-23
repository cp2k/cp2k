FROM ubuntu:rolling

# author: Ole Schuett

# TODO: maybe use conda instead?

RUN apt-get update && apt-get install -y --no-install-recommends \
    python                \
    python3               \
    python-pip            \
    python3-pip           \
    python-wheel          \
    python3-wheel         \
    python-setuptools     \
    python3-setuptools    \
    python-dev            \
    python3-dev           \
    build-essential       \
    ca-certificates       \
    wget                  \
&& rm -rf /var/lib/apt/lists/*

RUN pip install numpy matplotlib
RUN pip3 install numpy matplotlib

# install python2.6
WORKDIR /opt/
RUN wget https://www.python.org/ftp/python/2.6.9/Python-2.6.9.tgz  && \
    tar -xzf Python-2.6.9.tgz                                      && \
    cd Python-2.6.9                                                && \
    ./configure                                                    && \
    make -j                                                        && \
    make install                                                   && \
    cd ..                                                          && \
    rm -rf Python-2.6.9*

WORKDIR /opt/cp2k_test_python
COPY cmd.sh .
CMD ["./cmd.sh"]

#EOF