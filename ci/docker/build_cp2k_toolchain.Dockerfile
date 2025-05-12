# Dockerfile for CP2K continuous integration (CI) runs
#
# A stand-alone docker build in this folder can be performed using the command:
# docker build -f build_cp2k_toolchain.Dockerfile ../../
#
# Author: Matthias Krack
#
# Stage 2a: Build CP2K

ARG BASE_IMAGE="ubuntu:24.04"
ARG DEPS_IMAGE=""

FROM ${DEPS_IMAGE} AS build_cp2k

# Store build arguments from base image needed in next stage
RUN echo "${CP2K_VERSION}" >/CP2K_VERSION

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
RUN /bin/bash -c -o pipefail " \
    TOOLCHAIN_DIR=/opt/cp2k/tools/toolchain; \
    source ./cmake/cmake_cp2k.sh toolchain_${CP2K_BUILD_TYPE} ${CP2K_VERSION}"

# Compile CP2K
ARG LOG_LINES
ENV LOG_LINES=${LOG_LINES:-200}
WORKDIR /opt/cp2k/build
RUN /bin/bash -c -o pipefail " \
    source /opt/cp2k/tools/toolchain/install/setup; \
    echo -e '\nCompiling CP2K ... \c'; \
    if ninja --verbose &>ninja.log; then \
      echo -e 'done\n'; \
      echo -e 'Installing CP2K ... \c'; \
      if ninja --verbose install &>install.log; then \
        echo -e 'done\n'; \
      else \
        echo -e 'failed\n'; \
        tail -n ${LOG_LINES} install.log; \
      fi; \
      cat cmake.log ninja.log install.log | gzip >build_cp2k.log.gz; \
    else \
      echo -e 'failed\n'; \
      tail -n ${LOG_LINES} ninja.log; \
      cat cmake.log ninja.log | gzip >build_cp2k.log.gz; \
    fi"

# Update library search path
RUN echo "/opt/cp2k/lib" >/etc/ld.so.conf.d/cp2k.conf && ldconfig

# Collect components for installation
WORKDIR /opt/cp2k/bin
RUN /bin/bash -c -o pipefail " \
    source /opt/cp2k/tools/toolchain/install/setup; \
    mkdir -p /opt/toolchain/install /opt/toolchain/scripts; \
    for libdir in \$(ldd /opt/cp2k/bin/cp2k.$(cat /CP2K_VERSION) | \
                     grep /opt/cp2k/tools/toolchain/install | \
                     awk '{print \$3}' | cut -d/ -f7 | \
                     sort | uniq) setup; do \
      cp -ar /opt/cp2k/tools/toolchain/install/\${libdir} /opt/toolchain/install; \
    done; \
    cp /opt/cp2k/tools/toolchain/scripts/tool_kit.sh /opt/toolchain/scripts"

# Stage 2b: Install CP2K
FROM ${BASE_IMAGE} AS install_cp2k

# Install required packages
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    g++ \
    gcc \
    gfortran \
    python3 && rm -rf /var/lib/apt/lists/*

# Import build arguments from base image
COPY --from=build_cp2k /CP2K_VERSION /

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

# Install shared libraries required by the CP2K binaries
COPY --from=build_cp2k /opt/toolchain /opt/cp2k/tools/toolchain

# Import compressed build log file
COPY --from=build_cp2k /opt/cp2k/build/build_cp2k.log.gz /opt/cp2k/build/build_cp2k.log.gz

# Create links to CP2K binaries
WORKDIR /opt/cp2k/bin
RUN CP2K_VERSION=$(cat /CP2K_VERSION) && \
    ln -sf cp2k.${CP2K_VERSION} cp2k && \
    ln -sf cp2k.${CP2K_VERSION} cp2k.$(echo ${CP2K_VERSION} | sed "s/smp/opt/") && \
    ln -sf cp2k.${CP2K_VERSION} cp2k_shell

# Update library search path
RUN echo "/opt/cp2k/lib" >/etc/ld.so.conf.d/cp2k.conf && ldconfig

# Create entrypoint script file
RUN printf "#!/bin/bash\n\
ulimit -c 0 -s unlimited\n\
\
export OMP_STACKSIZE=64M\n\
export PATH=/opt/cp2k/bin:\${PATH}\n\
source /opt/cp2k/tools/toolchain/install/setup\n\
\"\$@\"" \
>/opt/cp2k/bin/entrypoint.sh && chmod 755 /opt/cp2k/bin/entrypoint.sh

# Create shortcut for regression test
RUN printf "/opt/cp2k/tests/do_regtest.py \$* /opt/cp2k/bin $(cat /CP2K_VERSION)" \
>/opt/cp2k/bin/run_tests && chmod 755 /opt/cp2k/bin/run_tests

# Define entrypoint
WORKDIR /mnt
ENTRYPOINT ["/opt/cp2k/bin/entrypoint.sh"]
CMD ["cp2k", "--help"]

# EOF
