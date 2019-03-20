#!/bin/bash -e

# author: Ole Schuett

set -eo pipefail

REPORT=/workspace/report.txt
echo -n "StartDate: " | tee -a $REPORT
date --utc --rfc-3339=seconds | tee -a $REPORT
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
          "$@" |& tee -a $REPORT
}

# Upload to cloud storage.
function upload_file {
    URL=$1
    FILE=$2
    CONTENT_TYPE=$3
    wget --quiet --output-document=- --method=PUT --header="content-type: ${CONTENT_TYPE}" --header="cache-control: no-cache" --body-file="${FILE}" "${URL}" > /dev/null
}

# Handle preemption gracefully.
function sigterm_handler {
    echo -e "\nThis job just got preempted. No worries, it should restart soon.\n" | tee -a $REPORT
    echo -n "EndDate: " | tee -a $REPORT
    date --utc --rfc-3339=seconds | tee -a $REPORT
    upload_file "${REPORT_UPLOAD_URL}" "${REPORT}" "text/plain;charset=utf-8"
    exit 1
}
trap sigterm_handler SIGTERM

# Upload preliminary report every 30s in the background.
(
while true ; do
    sleep 1
    count=$(( (count + 1) % 30 ))
    if (( count == 1 )) && [ -n "${REPORT_UPLOAD_URL}" ]; then
        upload_file "${REPORT_UPLOAD_URL}" "${REPORT}" "text/plain;charset=utf-8"
    fi
done
)&

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
    echo -e "\n========== Fetching Git Commit ==========" | tee -a $REPORT
    cd /workspace/cp2k
    git fetch origin "${GIT_BRANCH}"                       |& tee -a $REPORT
    git -c advice.detachedHead=false checkout "${GIT_REF}" |& tee -a $REPORT
    git submodule update --init --recursive                |& tee -a $REPORT
    git --no-pager log -1 --pretty='%nCommitSHA: %H%nCommitTime: %ci%nCommitAuthor: %an%nCommitSubject: %s%n' |& tee -a $REPORT

elif [ -d  /mnt/cp2k ]; then
    echo -e "\n========== Copying Changed Files ==========" | tee -a $REPORT
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

# Update toolchain, if present.
TOOLCHAIN_OK=true
if [ -d /opt/cp2k-toolchain ]; then
    echo -e "\n========== Updating Toolchain ==========" | tee -a $REPORT
    cd /opt/cp2k-toolchain
    rsync_changes --exclude="/build/"    \
                  --exclude="/install/"  \
                  /workspace/cp2k/tools/toolchain/  /opt/cp2k-toolchain/

    # shellcheck disable=SC1091
    source /opt/cp2k-toolchain/install/setup
    # shellcheck disable=SC2086
    if ! ./install_cp2k_toolchain.sh ${CP2K_TOOLCHAIN_OPTIONS} |& tee -a $REPORT ; then
        echo -e "\nSummary: Toolchain update failed." | tee -a $REPORT
        echo -e "Status: FAILED\n" | tee -a $REPORT
        TOOLCHAIN_OK=false
    fi
fi

if ! md5sum --status --check ${CHECKSUMS}; then
    echo -e "\n========== Cleaning Build Cache ==========" | tee -a $REPORT
    cd /workspace/cp2k
    make distclean |& tee -a $REPORT
fi

# Run actual test.
if $TOOLCHAIN_OK ; then
    echo -e "\n========== Running Test ==========" | tee -a $REPORT
    cd /workspace
    if ! "$@" |& tee -a $REPORT ; then
       echo -e "\nSummary: Test had non-zero exit status." | tee -a $REPORT
       echo -e "Status: FAILED\n" | tee -a $REPORT
    fi
fi

# Wrap up.
echo -n "EndDate: " | tee -a $REPORT
date --utc --rfc-3339=seconds | tee -a $REPORT

# Upload final report.
if [ -n "${REPORT_UPLOAD_URL}" ]; then
    echo "Uploading report..."
    upload_file "${REPORT_UPLOAD_URL}" "${REPORT}" "text/plain;charset=utf-8"
fi

# Upload artifacts.
if [ -n "${ARTIFACTS_UPLOAD_URL}" ] &&  [ -d /workspace/artifacts ]; then
    echo "Uploading artifacts..."
    ARTIFACTS_TGZ="/tmp/test_${TESTNAME}_artifacts.tgz"
    cd /workspace/artifacts
    tar -czf "${ARTIFACTS_TGZ}" *
    upload_file "${ARTIFACTS_UPLOAD_URL}" "${ARTIFACTS_TGZ}" "application/gzip"
fi

echo "Done :-)"

#EOF
