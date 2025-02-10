WORKDIR /opt/cp2k
COPY ./src ./src
COPY ./data ./data
COPY ./tests ./tests
COPY ./tools/build_utils ./tools/build_utils
COPY ./cmake ./cmake
COPY ./CMakeLists.txt .

# Build CP2K with CMake and run regression tests
COPY ./tools/docker/scripts/build_cp2k_cmake.sh ./
RUN /bin/bash -o pipefail -c " \
    ./build_cp2k_cmake.sh toolchain psmp |& tee build_cp2k.log"

#EOF
