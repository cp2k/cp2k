FROM cp2k/toolchain:latest

# author: Ole Schuett

# download and compile cp2k snapshot
WORKDIR /opt/
RUN wget -O cp2k-master.zip https://github.com/cp2k/cp2k/archive/master.zip && \
    unzip cp2k-master.zip  && \
    rm cp2k-master.zip

WORKDIR /opt/cp2k-master/cp2k/arch
RUN ln -vs /opt/cp2k-toolchain/install/arch/local* .

ARG ARCH
ARG VERSION
WORKDIR /opt/cp2k-master/cp2k/makefiles
# run regtests which lack fixed reference value
RUN source /opt/cp2k-toolchain/install/setup     && \
    make -j   ARCH=${ARCH} VERSION=${VERSION}    && \
    make test ARCH=${ARCH} VERSION=${VERSION} TESTOPTS="-restrictdir QS/regtest-almo-md -restrictdir QS/regtest-almo-1 -restrictdir SE/regtest-3-4 -restrictdir QS/regtest-ot-1-vib -restrictdir Fist/regtest-5-vib -restrictdir QS/regtest-optbas -restrictdir TMC/regtest_ana_post_proc" && \
    rm -rf ../lib/ ../exe/ ../regtesting/${ARCH}/${VERSION}/TEST-*

# install Debian packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    rsync                                                        \
  && rm -rf /var/lib/apt/lists/*

ARG TESTNAME
ENV TESTNAME=${TESTNAME}
ENV ARCH=${ARCH}
ENV VERSION=${VERSION}
WORKDIR /opt/cp2k_${TESTNAME}
COPY ./cmd.sh .
CMD ["./cmd.sh"]

#EOF
