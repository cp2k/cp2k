FROM cp2k/toolchain:latest

# Has to be build from cp2k root directory:
# docker build -f ./infra/test_conventions/Dockerfile .

WORKDIR /opt/cp2k/
COPY makefiles ./makefiles
COPY tools ./tools
COPY arch ./arch
COPY src ./src
RUN ln -s /opt/cp2k-toolchain/install/arch/* ./arch/

RUN ln -sf bash /bin/sh
WORKDIR ./tools/conventions
CMD source /opt/cp2k-toolchain/install/setup && ./test_conventions.sh

#EOF
