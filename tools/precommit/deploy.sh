#!/bin/bash -e

# author: Ole Schuett

set -x

SHORT_SHA=$(git rev-parse --short HEAD)
docker build --build-arg "REVISION=${SHORT_SHA}" -t cp2k-precommit .

docker tag cp2k-precommit gcr.io/cp2k-org-project/img_cp2kprecommit:${SHORT_SHA}
docker tag cp2k-precommit gcr.io/cp2k-org-project/img_cp2kprecommit:latest
docker push gcr.io/cp2k-org-project/img_cp2kprecommit:${SHORT_SHA}
docker push gcr.io/cp2k-org-project/img_cp2kprecommit:latest

gcloud run deploy cp2k-precommit --platform=managed --region=us-central1 \
    --image=gcr.io/cp2k-org-project/img_cp2kprecommit:${SHORT_SHA}

#EOF
