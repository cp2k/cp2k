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
RUN export DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true && \
  apt-get update && apt-get install -y --no-install-recommends \
    build-essential       \
    python-setuptools     \
    python-wheel          \
    python-pip            \
    python-dev            \
    postgresql            \
    less                  \
    nano                  \
    sudo                  \
    ssh                   \
    git                   \
    rsync                 \
  && rm -rf /var/lib/apt/lists/*

# install python packages
RUN pip install flake8 aiida ase

# create ubuntu user with sudo powers
RUN adduser --disabled-password --gecos "" ubuntu               && \
    echo "ubuntu ALL=(ALL) NOPASSWD: ALL" >>  /etc/sudoers

WORKDIR /opt/cp2k_test_aiida
COPY cmd.sh .
CMD ["./cmd.sh"]

#EOF
