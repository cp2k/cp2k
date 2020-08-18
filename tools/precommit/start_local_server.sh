#!/bin/bash -e

# author: Ole Schuett

set -x

SHORT_SHA=$(git rev-parse --short HEAD)
docker build --build-arg "REVISION=${SHORT_SHA}" -t cp2k-precommit .
docker run --rm -p127.0.0.1:8080:8080 cp2k-precommit

#EOF
