#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

# This script assumes a working environment as follows:
# cp2k                                  <- variable ${CP2K_ROOT}; CMake option -S
# ├── CMakeLists.txt                    <- * file to be parsed in this script
# ├── data                              <- CP2K data directory; CMake option -DCP2K_DATA_DIR
# ├── build                             <- to-be-created; CMake option -B
# ├── install                           <- to-be-created; CMake option -DCMAKE_INSTALL_PREFIX
# ├── src                               <- CP2K source code directory
# └── tools
#     └── toolchain                     <- working directory; variable ${ROOTDIR}
#         ├── install_cp2k_toolchain.sh <- script being executed calling this script
#         ├── scripts                   <- variable ${SCRIPT_DIR}
#         │   ├── common_vars.sh
#         │   ├── tool_kit.sh
#         │   └── generate_cmake_options.sh <- this script
#         └── install                   <- variable ${INSTALLDIR}
#             ├── setup                 <- * file to be parsed in this script
#             ├── toolchain.conf        <- * file to be parsed in this script
#             └── toolchain.env         <- * file to be parsed in this script
#
# First, validate completion of relevant upper-level directory and file
printf "\n========================== %s =========================\n" \
  "Generating CMake options for building CP2K"
CP2K_ROOT=$(cd "${ROOTDIR}/../.." && pwd)
if [ -d "${CP2K_ROOT}/src" ]; then
  cat << EOF
Root directory of CP2K with source code is found as ${CP2K_ROOT}
(path is exported to variable \${CP2K_ROOT}).
Build directory will be \${CP2K_ROOT}/build.
Install directory will be \${CP2K_ROOT}/install.
EOF
  CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=./install"
else
  report_error ${LINENO} "\${CP2K_ROOT} does not have subdirectory src."
  return 1
fi
if [ -f "${CP2K_ROOT}/CMakeLists.txt" ] && [ -r "${CP2K_ROOT}/CMakeLists.txt" ]; then
  echo "\${CP2K_ROOT}/CMakeLists.txt exists; will be parsed for CMake options."
else
  report_error ${LINENO} "\${CP2K_ROOT}/CMakeLists.txt cannot be found or read."
  return 1
fi
if [ -d "${CP2K_ROOT}/data" ]; then
  echo "Data directory ${CP2K_ROOT}/data is found and set as CP2K_DATA_DIR."
  CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_DATA_DIR=${CP2K_ROOT}/data"
else
  report_error ${LINENO} "Data directory \${CP2K_ROOT}/data cannot be found."
fi

# ------------------------------------------------------------------------
# generate cmake options for compiling cp2k
# ------------------------------------------------------------------------
# Build the program in source tree for convenience
if [ -n "$(grep -- "--install-all" ${INSTALLDIR}/setup)" ]; then
  CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_EVERYTHING=ON -DCP2K_USE_DLAF=OFF -DCP2K_USE_PEXSI=OFF"
  # Since "--install-all" can be used together with "--with-PKG=no", an extra safeguard is added here
  for toolchain_option in $(grep -i "dontuse" ${INSTALLDIR}/toolchain.conf | grep -vi "gcc" | cut -d'_' -f2 | cut -d'=' -f1); do
    if [ $(eval echo "$with_"'$toolchain_option') != "__DONTUSE__" ]; then
      ADDED_CMAKE_OPTION=$(sed -n '/option(/,/)/p' ${CP2K_ROOT}/CMakeLists.txt | grep -i $toolchain_option | awk '{print $1}' | cut -d'(' -f2 | head -n 1)
      # Use "if-then" below can avoid generating empty "-D=OFF" options
      if [ -n "${ADDED_CMAKE_OPTION}" ]; then
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -D${ADDED_CMAKE_OPTION}=OFF"
      fi
    fi
  done
else
  # If MPI is used, set "CP2K_USE_MPI" to "ON"
  if [ "${MPI_MODE}" != "no" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_MPI=ON"
    if [ -n "$(grep "MPI_F08" "${INSTALLDIR}"/toolchain.env)" ]; then
      CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_MPI_F08=ON"
    fi
  fi
  # If GPU acceleration is used, add the option about GPU acceleration
  if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ACCEL=CUDA -DCP2K_WITH_GPU=${GPUVER}"
  elif [ "${ENABLE_HIP}" = "__TRUE__" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ACCEL=HIP -DCP2K_WITH_GPU=${GPUVER}"
  elif [ "${ENABLE_OPENCL}" = "__TRUE__" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ACCEL=OPENCL"
  fi
  # Some options that should be specially considered:
  # Intel MKL includes FFTW
  if [ "${with_fftw}" != "__DONTUSE__" ] || [ "${MATH_MODE}" = "mkl" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_FFTW3=ON"
  fi
  # There is only MCL (MiMiC Communication Library) and no MiMiC that can be installed via toolchain, but if one chooses to install MCL then MiMiC should have been installed?
  if [ "${with_mcl}" != "__DONTUSE__" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_MIMIC=ON"
  fi
  # Detect if any other dependencies is used and add the proper cmake option. Since "pugixml" and "gsl" are not mentioned in CMakeLists.txt, they will not be considered.
  for toolchain_option in $(grep -vi "dry\|list\|dontuse\|gcc\|cmake\|fftw\|mkl\|dbcsr" ${INSTALLDIR}/toolchain.conf | cut -d'_' -f2 | cut -d'=' -f1); do
    if [ $(eval echo "$with_"'$toolchain_option') != "__DONTUSE__" ]; then
      ADDED_CMAKE_OPTION=$(sed -n '/option(/,/)/p' ${CP2K_ROOT}/CMakeLists.txt | grep -i $toolchain_option | awk '{print $1}' | cut -d'(' -f2 | head -n 1)
      # Use "if-then" below can avoid generating empty "-D=ON" options
      if [ -n "${ADDED_CMAKE_OPTION}" ]; then
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -D${ADDED_CMAKE_OPTION}=ON"
      fi
    fi
  done
fi

# Export variable for CMake options to setup file
cat << EOF >> "${SETUPFILE}"
# ==================== Setup for CP2K ==================== #
export CP2K_ROOT="${CP2K_ROOT}"
export CP2K_CMAKE_OPTIONS="${CMAKE_OPTIONS}"
prepend_path PATH "${CP2K_ROOT}/install/bin"
prepend_path LD_LIBRARY_PATH "${CP2K_ROOT}/install/lib"
prepend_path PKG_CONFIG_PATH "${CP2K_ROOT}/install/lib/pkgconfig"
EOF
if [ "${dry_run}" = "__TRUE__" ]; then
  cat << EOF

Suggested cmake command if toolchain is built with your options:

  cmake ${CMAKE_OPTIONS}
EOF
else
  cat << EOF

Suggested CMake options are collected in the variable \${CP2K_CMAKE_OPTIONS} that is
exported at the end of setup file ${SETUPFILE}.
EOF
fi

# -------------------------
# print out user instructions
# -------------------------
if [ "${dry_run}" != "__TRUE__" ]; then
  echo
  cat << EOF
========================== Epilogue =========================
Toolchain is now ready for building CP2K!

To use the installed tools and libraries and cp2k version compiled with it you
will first need to execute at the prompt:

  source ${SETUPFILE}

Then it's recommended for you to build and install CP2K with following commands:

  cd ${CP2K_ROOT}
  cmake -S . -B build ${CMAKE_OPTIONS}
  cmake --build build --target install -j $(get_nprocs)

For more information about available build options, see:
https://manual.cp2k.org/trunk/getting-started/build-from-source.html

For detailed explanation of above steps of building CP2K, see:
${INSTALLDIR}/cp2k_installation_guide.md

EOF
  cat << EOF > ${INSTALLDIR}/cp2k_installation_guide.md
## Building CP2K with dependencies installed in toolchain

Here is a detailed instruction of building CP2K with dependencies installed via toolchain.

### Required - commands you must execute or the building fails

1. Source setup file to activate toolchain-configured dependencies:

\`\`\`bash
source ${SETUPFILE}
\`\`\`

This setup file **MUST** also be sourced whenever CP2K built with this toolchain is executed.
If modules have been used to estabilish environment variables and paths, remember to load
these modules prior to sourcing setup file.

2. Go to root directory of CP2K and configure CMake with corresponding options:

\`\`\`bash
cd ${CP2K_ROOT}
cmake -S . -B build ${CMAKE_OPTIONS}
\`\`\`

Other commands from \`${CP2K_ROOT}/CMakeLists.txt\` can also be added. For more information
about available build options, see documentation:
<https://manual.cp2k.org/trunk/getting-started/build-from-source.html>.

Alternative to copy-paste long lines in terminal is to use a variable from setup file for
CMake options, which is not to be quoted so that whitespace delimiters allow it to expand
to command options in shell:

\`\`\`bash
cmake -S . -B build \${CP2K_CMAKE_OPTIONS}
\`\`\`

3. Build and install CP2K with command:

\`\`\`bash
cmake --build build --target install -j $(get_nprocs)
\`\`\`

It may be helpful to also save a copy of command line messages to log files:

\`\`\`bash
cmake --build build --target install -j $(get_nprocs) 2>&1 | tee install.log
\`\`\`

If you want another build of CP2K without changing toolchain configuration, simply change
the building directory set by \`-B <dirname>\` in the cmake command above.

### Optional but recommended

At the ending of the output of step 3, CP2K will give you a command that can be used to
do regtests, which can further ensure if you have built CP2K correctly. You can perform
the regtests with that command. You can run
\`${CP2K_ROOT}/tests/do_regtest.py --help\`
to see available options. For a detailed instruction of how to run regtests (especially
on HPC clusters), see:
<https://www.cp2k.org/dev:regtesting>.

### Other optional behaviors

Both for toolchain buildings and CP2K building: once build is completed, the \`build\`
directory can be safely deleted. However, you **MUST** keep the \`install\` directory
as is.

Especially, if you want to save disk space but at the same time keep the cached CMake
files, you can run \`cmake --build build --target clean\` alternatively after installing
successfully.
EOF
fi

#EOF
