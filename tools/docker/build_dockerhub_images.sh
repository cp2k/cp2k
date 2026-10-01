#!/bin/bash -e

# Usage: ./build_dockerhub_images.sh <cp2k_version>

if [[ -n $1 ]]; then
  VERSION=$1
else
  SHA=$(git rev-parse HEAD)
  DATE=$(git show -s --format=%cs "${SHA}" | sed s/-//g)
  VERSION="dev${DATE}"
fi

./spack_cache_start.sh

# Build OpenMPI container image
TAG="cp2k/cp2k:${VERSION}_openmpi_cascadelake_psmp"
podman build --shm-size=1g --build-arg "GIT_COMMIT_SHA=${SHA}" -t "${TAG}" -f ./Dockerfile.test_spack_openmpi-psmp ../../
echo -e "\n*** Run \"podman push ${TAG} docker.io/${TAG}\" to upload the container image to Dockerhub\n"

# Build MPICH container image
TAG="cp2k/cp2k:${VERSION}_mpich_cascadelake_psmp"
podman build --shm-size=1g --build-arg "GIT_COMMIT_SHA=${SHA}" -t "${TAG}" -f ./Dockerfile.test_spack_psmp-gcc15 ../../
echo -e "\n*** Run \"podman push ${TAG} docker.io/${TAG}\" to upload the container image to Dockerhub"

echo -e "\n*** Run \"podman push ${TAG} docker.io/cp2k/cp2k:latest\" to tag the latest container image on Dockerhub\n"
