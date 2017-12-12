FROM cp2k/cp2k:trunk

RUN apt-get update && apt-get install -y --no-install-recommends \
    git                   \
    python-dev            \
    python-pip            \
    python-wheel          \
    python-setuptools     \
    build-essential       \
&& rm -rf /var/lib/apt/lists/*

WORKDIR /opt/ase
RUN git clone https://gitlab.com/ase/ase.git .
RUN pip install .

#RUN pip install --upgrade git+https://gitlab.com/ase/ase.git@master

COPY run_ase_tests.sh .
CMD ["/opt/ase/run_ase_tests.sh"]

#EOF