FROM cp2k/cp2k:trunk

# install Debian packages
RUN apt-get update && apt-get install -y --no-install-recommends \
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
  && rm -rf /var/lib/apt/lists/*

# install python dependencies early to leverage docker build cache
RUN pip install flake8 aiida ase

# install aiida-cp2k
WORKDIR /opt/aiida-cp2k/
RUN git clone https://github.com/cp2k/aiida-cp2k.git .
RUN pip install .

# create ubuntu user with sudo powers
RUN adduser --disabled-password --gecos "" ubuntu               && \
    echo "ubuntu ALL=(ALL) NOPASSWD: ALL" >>  /etc/sudoers


# configure aiida
RUN ln -s /opt/cp2k/exe/local/cp2k.psmp /usr/bin/cp2k
USER ubuntu
WORKDIR /opt/aiida-cp2k/test/
RUN ./configure_aiida.sh

CMD ["./run_tests.sh"]

#EOF
