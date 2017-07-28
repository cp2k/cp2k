FROM cp2k/toolchain:latest

WORKDIR /cp2k/
RUN ln -s /cp2k-toolchain/install/arch/ .
COPY makefiles ./makefiles
COPY tools ./tools
COPY tests ./tests
COPY data ./data
COPY src ./src

WORKDIR /cp2k/makefiles
RUN make -j VERSION="sdbg" && make VERSION="sdbg" clean
RUN make -j VERSION="sopt" && make VERSION="sopt" clean
RUN make -j VERSION="pdbg" && make VERSION="pdbg" clean
RUN make -j VERSION="popt" && make VERSION="popt" clean
RUN make -j VERSION="ssmp" && make VERSION="ssmp" clean
RUN make -j VERSION="psmp" && make VERSION="psmp" clean

WORKDIR /cp2k/exe
#EOF