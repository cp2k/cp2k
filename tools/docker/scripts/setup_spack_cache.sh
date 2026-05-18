#!/bin/bash -e

# author: Ole Schuett

if [[ -n "${SPACK_CACHE}" ]]; then
  if [[ "${SPACK_CACHE}" == *"http://host.containers.internal:9000"* ]]; then
    if ! wget -q --tries=1 "http://host.containers.internal:9000/spack-cache"; then
      echo ""
      echo "ERROR: Could not connect to local Spack cache."
      echo "       Start the cache by running ./spack_cache_start.sh."
      echo "       Alternatively, disable the cache by passing --build-arg SPACK_CACHE=\"\" to podman."
      echo ""
      exit 1
    fi
  fi
  echo "Adding Spack cache: ${SPACK_CACHE}"
  # shellcheck disable=SC2086
  spack mirror add --autopush --unsigned local-cache ${SPACK_CACHE}
else
  echo "No Spack cache provided."
fi

#EOF
