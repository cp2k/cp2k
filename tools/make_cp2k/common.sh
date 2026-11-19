# shellcheck shell=bash
#
# Helper functions shared by make_cp2k.sh and the step scripts in tools/make_cp2k
#
# This file only defines functions and is sourced by each script.
# The scripts exchange their settings only via configuration files with
# KEY=value lines (values quoted with printf %q) and unset KEY lines for unset
# variables, which are written by write_config and read by load_config:
#
#   options.conf : Resolved settings of make_cp2k.sh (step 1), also printed by --print_config
#   deps.conf    : Settings provided by the CP2K dependency build with build_deps.sh (step 2)
#   build.conf   : Settings provided by the CP2K build with build_cp2k.sh (step 3)

# Run a command that downloads data with up to 3 attempts and a backoff
# of 60 s and 120 s to bridge transient network issues
retry_download() {
  local attempt delay=60
  for attempt in 1 2 3; do
    "$@" && return 0
    if ((attempt < 3)); then
      echo "WARNING: Attempt ${attempt} of 3 failed, retrying in ${delay} seconds: $*"
      sleep "${delay}"
      delay=$((delay * 2))
    fi
  done
  return 1
}

# Write the given variables as KEY=value lines to stdout, or as unset KEY lines
# for unset variables, so that loading a configuration reproduces its state
write_config() {
  local name
  for name in "$@"; do
    if [[ -v ${name} ]]; then
      printf '%s=%q\n' "${name}" "${!name}"
    else
      printf 'unset %s\n' "${name}"
    fi
  done
}

# Read and export the variables of the given configuration files in the given order
load_config() {
  local config_file
  if (($# == 0)); then
    echo "ERROR: No configuration file specified"
    echo "Usage: $(basename "${0}") CONFIG_FILE..."
    exit 1
  fi
  for config_file in "$@"; do
    if [[ ! -f "${config_file}" ]]; then
      echo "ERROR: The configuration file ${config_file} was not found"
      exit 1
    fi
    set -a
    # shellcheck source=/dev/null
    source "${config_file}"
    set +a
  done
}

# Activate the Spack environment with the CP2K dependencies
activate_spack_env() {
  # shellcheck source=/dev/null
  source "${SPACK_ROOT}/share/spack/setup-env.sh"
  if ! spack env list | grep -q "${CP2K_ENV}"; then
    echo "ERROR: No Spack environment \"${CP2K_ENV}\" found"
    exit 1
  fi
  eval "$(spack env activate --sh "${CP2K_ENV}")"
}
