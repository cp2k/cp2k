#!/bin/bash -e

# author: Ole Schuett

if podman start spack-cache; then
  echo "Re-started existing spack cache."

else
  # Start MinIO server.
  podman run --name spack-cache --detach -p 9000:9000 -p 9001:9001 \
    quay.io/minio/minio server /data --console-address ":9001"

  sleep 3
  podman container logs spack-cache

  # Configure localhost.
  podman exec spack-cache mc alias set local http://localhost:9000 minioadmin minioadmin

  # Create bucket.
  podman exec spack-cache mc mb local/spack-cache

  # Make bucket public.
  podman exec spack-cache mc anonymous set public local/spack-cache
  echo "Started new spack cache."
fi

#EOF
