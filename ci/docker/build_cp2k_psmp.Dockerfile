# Stage 1b: Build CP2K with CMake
WORKDIR /opt/cp2k
COPY ./CMakeLists.txt .
COPY ./cmake ./cmake
COPY ./data ./data
COPY ./src ./src
COPY ./tests ./tests
COPY ./tools/build_utils ./tools/build_utils

# Run CMake
RUN /bin/bash -c -o pipefail " \
    TOOLCHAIN_DIR=/opt/cp2k/tools/toolchain; \
    source ./cmake/cmake_cp2k.sh toolchain psmp |& tee cmake.log"

# Compile CP2K
WORKDIR /opt/cp2k/build
RUN /bin/bash -c -o pipefail " \
    source /opt/cp2k/tools/toolchain/install/setup; \
    echo -en '\nCompiling CP2K ... '; \
    if ninja --verbose &> ninja.log; then \
      echo -e 'done\n'; \
      echo -en 'Installing CP2K ... '; \
      if ninja --verbose install &> install.log; then \
        echo -e 'done\n'; \
      else \
        echo -e 'failed\n'; \
        cat install.log; \
      fi; \
    else \
      echo -e 'failed\n'; \
      cat ninja.log; \
    fi"

# Update environment variables
ENV LD_LIBRARY_PATH=/opt/cp2k/lib:${LD_LIBRARY} \
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

# Stage 2: Install CP2K
FROM ubuntu:24.04 AS install

# Install required packages
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    g++ \
    gcc \
    gfortran \
    python3 && rm -rf /var/lib/apt/lists/*

# Install CP2K binaries and libraries
WORKDIR /opt/cp2k
COPY --from=build /opt/cp2k/bin ./bin
COPY --from=build /opt/cp2k/lib ./lib

# Install CP2K database files
COPY --from=build /opt/cp2k/share ./share

# Install CP2K regression tests
COPY --from=build /opt/cp2k/tests ./tests
COPY --from=build /opt/cp2k/src/grid/sample_tasks ./src/grid/sample_tasks

# Create links to CP2K binaries
WORKDIR /opt/cp2k/bin
RUN ln -sf cp2k.psmp cp2k
RUN ln -sf cp2k.psmp cp2k.popt
RUN ln -sf cp2k.psmp cp2k_shell

# Install shared libraries required by the CP2K binaries
COPY --from=build /toolchain /opt/cp2k/tools/toolchain

# Create entrypoint script file
RUN printf "#!/bin/bash\n\
ulimit -c 0 -s unlimited\n\
\
export OMP_STACKSIZE=64M\n\
export LD_LIBRARY_PATH=/opt/cp2k/lib:\${LD_LIBRARY_PATH}\n\
export PATH=/opt/cp2k/bin:\${PATH}\n\
source /opt/cp2k/tools/toolchain/install/setup\n\
\"\$@\"" \
>/opt/cp2k/bin/entrypoint.sh && chmod 755 /opt/cp2k/bin/entrypoint.sh

# Create shortcut for regression test
RUN printf "/opt/cp2k/tests/do_regtest.py --maxtasks 12 --workbasedir /mnt \$* /opt/cp2k/bin psmp" \
>/opt/cp2k/bin/run_tests && chmod 755 /opt/cp2k/bin/run_tests

# Define entrypoint
WORKDIR /mnt
ENTRYPOINT ["/opt/cp2k/bin/entrypoint.sh"]
CMD ["cp2k", "--help"]

# Label docker image
LABEL author="CP2K Developers" \
      cp2k_version="master" \
      container_type="production_psmp" \
      dockerfile_generator_version="0.1"

# EOF
