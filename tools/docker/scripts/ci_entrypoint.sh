#!/bin/bash -e

# author: Ole Schuett

set -eo pipefail
ulimit -c 0  # Disable core dumps as they can take a very long time to write.

# Rsync with common args.
function rsync_changes {
    rsync --update                 \
          --delete                 \
          --executability          \
          --ignore-times           \
          --verbose                \
          --recursive              \
          --checksum               \
          "$@"
}

# Upload to cloud storage.
function upload_file {
    URL=$1
    FILE=$2
    CONTENT_TYPE=$3
    wget --quiet --output-document=- --method=PUT --header="content-type: ${CONTENT_TYPE}" --header="cache-control: no-cache" --body-file="${FILE}" "${URL}" > /dev/null
}

# Calculate checksums of critical files.
CHECKSUMS=/workspace/checksums.md5
shopt -s nullglob  # ignore missing files
md5sum /workspace/cp2k/Makefile \
       /workspace/cp2k/tools/build_utils/* \
       /workspace/cp2k/arch/local* \
       > $CHECKSUMS
shopt -u nullglob

# Get cp2k sources.
if [ -n "${GIT_REF}" ]; then
    echo -e "\n========== Fetching Git Commit =========="
    cd /workspace/cp2k
    git fetch origin "${GIT_BRANCH}"
    git -c advice.detachedHead=false checkout "${GIT_REF}"
    git submodule update --init --recursive
    git --no-pager log -1 --pretty='%nCommitSHA: %H%nCommitTime: %ci%nCommitAuthor: %an%nCommitSubject: %s%n'

elif [ -d  /mnt/cp2k ]; then
    echo -e "\n========== Copying Changed Files =========="
    rsync_changes --exclude="*~"                         \
                  --exclude=".*/"                        \
                  --exclude="*.py[cod]"                  \
                  --exclude="__pycache__"                \
                  --exclude="/obj/"                      \
                  --exclude="/lib/"                      \
                  --exclude="/exe/"                      \
                  --exclude="/regtesting/"               \
                  --exclude="/arch/local*"               \
                  --exclude="/tools/toolchain/build/"    \
                  --exclude="/tools/toolchain/install/"  \
                  /mnt/cp2k/  /workspace/cp2k/
else
    echo "Neither GIT_REF nor /mnt/cp2k found - aborting."
    exit 255
fi

if ! md5sum --status --check ${CHECKSUMS}; then
    echo -e "\n========== Cleaning Build Cache =========="
    cd /workspace/cp2k
    make distclean
fi

# Run actual test.
echo -e "\n========== Running Test =========="
cd /workspace
if ! "$@" ; then
   echo -e "\nSummary: Test had non-zero exit status.\nStatus: FAILED"
fi

# Upload artifacts.
if [ -n "${ARTIFACTS_UPLOAD_URL}" ] &&  [ -d /workspace/artifacts ]; then
    echo "Uploading artifacts..."
    ARTIFACTS_TGZ="/tmp/test_${TESTNAME}_artifacts.tgz"
    cd /workspace/artifacts
    tar -czf "${ARTIFACTS_TGZ}" -- *
    upload_file "${ARTIFACTS_UPLOAD_URL}" "${ARTIFACTS_TGZ}" "application/gzip"
fi

echo "Done :-)"

#EOF
