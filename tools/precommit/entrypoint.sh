#!/bin/bash -e

# author: Ole Schuett

trap "exit" INT

if [ -z "${GIT_REF}" ]; then
    # Server mode.
    gunicorn --bind=:8080 --workers=1 --threads=8 --timeout=0 precommit_server:app &
    wait  # Can be interrupted with ctrl-C.
else
    # CI mode.
    echo -e "\n========== Starting Precommit Server =========="
    gunicorn --bind=:8080 --workers=1 --threads=8 --timeout=0 precommit_server:app &> server.logs &
    sleep 3
    cat server.logs
    echo -e "\n========== Fetching Git Commit =========="
    cd /workspace/cp2k
    git fetch --quiet origin "${GIT_BRANCH}"
    git reset --quiet --hard "${GIT_REF}"
    git submodule update --init --recursive
    git --no-pager log -1 --pretty='%nCommitSHA: %H%nCommitTime: %ci%nCommitAuthor: %an%nCommitSubject: %s%n'
    echo -e "\n========== Running Precommit Checks =========="
    export CP2K_PRECOMMIT_SERVER="http://127.0.0.1:8080"
    ./tools/precommit/precommit.py --no-cache --progressbar-wait=10  || true
fi

#EOF
