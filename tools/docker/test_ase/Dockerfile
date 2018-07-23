FROM cp2k/toolchain:latest

# author: Ole Schuett

# download and compile cp2k snapshot
WORKDIR /opt/
RUN wget -O cp2k-master.zip https://github.com/cp2k/cp2k/archive/master.zip && \
    unzip cp2k-master.zip  && \
    rm cp2k-master.zip

WORKDIR /opt/cp2k-master/cp2k/arch
RUN ln -vs /opt/cp2k-toolchain/install/arch/local* .

WORKDIR /opt/cp2k-master/cp2k/makefiles
RUN source /opt/cp2k-toolchain/install/setup  && \
    make -j VERSION=pdbg                      && \
    rm -rf ../lib/ ../exe/

# install Debian packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3                                                      \
    python3-dev                                                  \
    python3-pip                                                  \
    python3-wheel                                                \
    python3-setuptools                                           \
    build-essential                                              \
    rsync                                                        \
    git                                                          \
  && rm -rf /var/lib/apt/lists/*

# install python packages
RUN pip3 install numpy scipy matplotlib flask

# clone ase reprository
WORKDIR /opt/
RUN git clone https://gitlab.com/ase/ase.git

WORKDIR /opt/cp2k_test_ase
COPY cmd.sh .
CMD ["./cmd.sh"]

#EOF
