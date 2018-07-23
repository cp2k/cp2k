FROM cp2k/toolchain:latest

# author: Ole Schuett

# download and compile cp2k snapshot
WORKDIR /opt/
RUN wget -q -O cp2k-master.zip https://github.com/cp2k/cp2k/archive/master.zip && \
    unzip -q cp2k-master.zip  && \
    rm cp2k-master.zip

WORKDIR /opt/cp2k-master/cp2k/arch
RUN ln -vs /opt/cp2k-toolchain/install/arch/local* .

ARG ARCH
ARG VERSION
WORKDIR /opt/cp2k-master/cp2k/makefiles
# run regtests which lack fixed reference value
# Disable LeakSanitizer during docker build as it requires ptrace capabilities.
RUN source /opt/cp2k-toolchain/install/setup     && \
    make -j   ARCH=${ARCH} VERSION=${VERSION}    && \
    export LSAN_OPTIONS="detect_leaks=0"         && \
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
WORKDIR /opt/cp2k_test_${TESTNAME}
COPY ./cmd.sh .
CMD ["./cmd.sh"]

#EOF
