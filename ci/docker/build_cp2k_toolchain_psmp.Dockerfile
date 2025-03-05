# Dockerfile for CP2K continous integration (CI) runs
#
# A stand-alone docker build in this folder can be performed using the command:
# docker build -f build_cp2k_toolchain_psmp.Dockerfile ../../
#
# Author: Matthias Krack
#
# Stage 2a: Build CP2K

ARG BASE_IMAGE="ubuntu:24.04"
ARG DEPS_IMAGE="mkrack/cp2k_eiger:toolchain_psmp_base_image"

FROM ${DEPS_IMAGE} AS build_cp2k

# Retrieve the maximum number log file lines printed on error
ARG LOG_LINES=100

# Build CP2K with CMake
WORKDIR /opt/cp2k
COPY ./CMakeLists.txt ./
COPY ./benchmarks/CI ./benchmarks/CI
COPY ./cmake ./cmake
COPY ./data ./data
COPY ./src ./src
COPY ./tests ./tests
COPY ./tools/build_utils ./tools/build_utils

# Run CMake
ARG BUILD_TYPE="minimal"
RUN /bin/bash -c -o pipefail " \
    TOOLCHAIN_DIR=/opt/cp2k/tools/toolchain; \
    source ./cmake/cmake_cp2k.sh "${BUILD_TYPE}" psmp"

# Compile CP2K
WORKDIR /opt/cp2k/build
RUN /bin/bash -c -o pipefail " \
    source /opt/cp2k/tools/toolchain/install/setup; \
    echo -e '\nCompiling CP2K ... \c'; \
    if ninja --verbose &>build.log; then \
      echo -e 'done\n'; \
      echo -e 'Installing CP2K ... \c'; \
      if ninja --verbose install &>>build.log; then \
        echo -e 'done\n'; \
      else \
        echo -e 'failed\n'; \
        tail -n ${LOG_LINES} build.log; \
      fi; \
    else \
      echo -e 'failed\n'; \
      tail -n ${LOG_LINES} build.log; \
    fi; \
    gzip build.log"

# Update environment variables
ENV LD_LIBRARY_PATH=/opt/cp2k/lib:/usr/local/lib:${LD_LIBRARY_PATH} \
    PATH=/opt/cp2k/bin:${PATH}

# Collect components for installation
WORKDIR /opt/cp2k/bin
RUN /bin/bash -c -o pipefail " \
    source /opt/cp2k/tools/toolchain/install/setup; \
    mkdir -p /toolchain/install /toolchain/scripts; \
    for libdir in \$(ldd /opt/cp2k/bin/cp2k.psmp | \
                     grep /opt/cp2k/tools/toolchain/install | \
                     awk '{print \$3}' | cut -d/ -f7 | \
                     sort | uniq) setup; do \
      cp -ar /opt/cp2k/tools/toolchain/install/\${libdir} /toolchain/install; \
    done; \
    cp /opt/cp2k/tools/toolchain/scripts/tool_kit.sh /toolchain/scripts"

# Stage 2b: Install CP2K
FROM ${BASE_IMAGE} AS install_cp2k

# Install required packages
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    g++ \
    gcc \
    gfortran \
    python3 && rm -rf /var/lib/apt/lists/*

# Copy MPI installation
COPY --from=build_cp2k /usr/local/bin /usr/local/bin
COPY --from=build_cp2k /usr/local/include /usr/local/include
COPY --from=build_cp2k /usr/local/lib /usr/local/lib

# Install CP2K binaries
WORKDIR /opt/cp2k
COPY --from=build_cp2k /opt/cp2k/bin ./bin

# Install CP2K libraries
COPY --from=build_cp2k /opt/cp2k/lib ./lib

# Install CP2K database files
COPY --from=build_cp2k /opt/cp2k/share ./share

# Install CP2K regression tests
COPY --from=build_cp2k /opt/cp2k/tests ./tests
COPY --from=build_cp2k /opt/cp2k/src/grid/sample_tasks ./src/grid/sample_tasks

# Install CP2K/Quickstep CI benchmarks
COPY --from=build_cp2k /opt/cp2k/benchmarks/CI ./benchmarks/CI

# Create links to CP2K binaries
WORKDIR /opt/cp2k/bin
RUN ln -sf cp2k.psmp cp2k && \
    ln -sf cp2k.psmp cp2k.popt && \
    ln -sf cp2k.psmp cp2k_shell

# Install shared libraries required by the CP2K binaries
COPY --from=build_cp2k /toolchain /opt/cp2k/tools/toolchain

# Import build log file
COPY --from=build_cp2k /opt/cp2k/build/build.log.gz /opt/cp2k/

# Create entrypoint script file
RUN printf "#!/bin/bash\n\
ulimit -c 0 -s unlimited\n\
\
export OMP_STACKSIZE=64M\n\
export LD_LIBRARY_PATH=/opt/cp2k/lib:/usr/local/lib:\${LD_LIBRARY_PATH}\n\
export PATH=/opt/cp2k/bin:\${PATH}\n\
source /opt/cp2k/tools/toolchain/install/setup\n\
\"\$@\"" \
>/opt/cp2k/bin/entrypoint.sh && chmod 755 /opt/cp2k/bin/entrypoint.sh

# Create shortcut for regression test
RUN printf "/opt/cp2k/tests/do_regtest.py \$* /opt/cp2k/bin psmp" \
>/opt/cp2k/bin/run_tests && chmod 755 /opt/cp2k/bin/run_tests

# Define entrypoint
WORKDIR /mnt
ENTRYPOINT ["/opt/cp2k/bin/entrypoint.sh"]
CMD ["cp2k", "--help"]

# Label docker image
LABEL author="CP2K Developers" \
      cp2k_version="master" \
      image_type="toolchain_psmp"

# EOF
