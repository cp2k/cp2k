#!/bin/bash -e

# author: Ole Schuett

SHA=$(git rev-parse HEAD)
DATE=$(git show -s --format=%cs "${SHA}" | sed s/-//g)
TAG="cp2k/cp2k:dev${DATE}"

set -x
docker build --shm-size=1g --build-arg "GIT_COMMIT_SHA=${SHA}" -f Dockerfile.prod_generic_psmp -t "${TAG}" ../../

docker run "${TAG}" cp2k --version

#EOF
