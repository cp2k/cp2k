# Dockerfile for CP2K continous integration (CI) runs
#
# A stand-alone docker build in this folder can be performed using the command:
# docker build -f build_cp2k_spack_psmp.Dockerfile ../../
#
# Author: Matthias Krack
#
# Stage 2a: Build CP2K

ARG BASE_IMAGE="ubuntu:24.04"
ARG DEPS_IMAGE=""

FROM ${DEPS_IMAGE} AS build_cp2k

# Retrieve the maximum number log file lines printed on error
ARG LOG_LINES
ENV LOG_LINES=${LOG_LINES:-100}

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
ARG CP2K_BUILD_TYPE
ENV CP2K_BUILD_TYPE=${CP2K_BUILD_TYPE:-minimal}
RUN /bin/bash -c -o pipefail "source ./cmake/cmake_cp2k.sh spack_${CP2K_BUILD_TYPE} psmp"

# Compile CP2K
WORKDIR /opt/cp2k/build
RUN /bin/bash -c -o pipefail " \
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
    else \
      echo -e 'failed\n'; \
      tail -n ${LOG_LINES} ninja.log; \
    fi; \
    cat cmake.log ninja.log install.log | gzip >build_cp2k.log.gz"

# Stage 2b: Install CP2K
FROM ${BASE_IMAGE} AS install_cp2k

# Install required packages
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    g++ \
    gcc \
    gfortran \
    hwloc \
    libhwloc-dev \
    python3 && rm -rf /var/lib/apt/lists/*

# Install CP2K dependencies built with Spack
WORKDIR /opt
COPY --from=build_cp2k /opt/spack ./spack

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

# Import compressed build log file
COPY --from=build_cp2k /opt/cp2k/build/build_cp2k.log.gz /opt/cp2k/build/build_cp2k.log.gz

# Create links to CP2K binaries
WORKDIR /opt/cp2k/bin
RUN ln -sf cp2k.psmp cp2k && \
    ln -sf cp2k.psmp cp2k.popt && \
    ln -sf cp2k.psmp cp2k_shell

# Update library search path
RUN echo "/opt/cp2k/lib\n/opt/spack/lib\n$(dirname $(find /opt/spack/lib -name libtorch.so))" >/etc/ld.so.conf.d/cp2k.conf && ldconfig

# Create entrypoint script file
RUN printf "#!/bin/bash\n\
ulimit -c 0 -s unlimited\n\
\
export OMP_STACKSIZE=64M\n\
export PATH=/opt/cp2k/bin:/opt/spack/bin:\${PATH}\n\
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
      image_type="spack_psmp"

# EOF
