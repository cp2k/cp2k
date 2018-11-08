#!/bin/bash -e

# author: Ole Schuett

REPORT=/workspace/report.txt
echo -n "StartDate: " > $REPORT
date --utc --rfc-3339=seconds >> $REPORT

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
    wget -O- --method=PUT --header="content-type: ${CONTENT_TYPE}" --header="cache-control: no-cache" --body-file="${FILE}" "${URL}"
}

# Get cp2k sources.
if [ -n "${GIT_REF}" ]; then
    echo -e "\n========== Fetching Git Commit =========="
    cd /workspace/cp2k
    git fetch origin "${GIT_BRANCH}"
    git checkout "${GIT_REF}"
    git submodule update --init --recursive
    git --no-pager log -1 --pretty='%nCommitSHA: %H%nCommitTime: %ci%nCommitAuthor: %an%nCommitSubject: %s%n' | tee -a $REPORT

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


# Update toolchain.
if [ -d /opt/cp2k-toolchain ]; then
    echo -e "\n========== Updating Toolchain =========="
    cd /opt/cp2k-toolchain
    rsync_changes --exclude="/build/"    \
                  --exclude="/install/"  \
                  /workspace/cp2k/tools/toolchain/  /opt/cp2k-toolchain/

    # shellcheck disable=SC1091
    source /opt/cp2k-toolchain/install/setup
    # shellcheck disable=SC2086
    ./install_cp2k_toolchain.sh ${CP2K_TOOLCHAIN_OPTIONS}
fi

echo -e "\n========== Running Test =========="
cd /workspace
"$@" 2>&1 | tee -a $REPORT &  # Launch in the background.

# Upload preliminary report every 30s while test is running.
while jobs %1 &> /dev/null ; do
    sleep 1
    count=$(( (count + 1) % 30 ))
    if (( count == 1 )) && [ -n "${REPORT_UPLOAD_URL}" ]; then
        upload_file "${REPORT_UPLOAD_URL}" "${REPORT}" "text/plain;charset=utf-8"
    fi
done

# Test has finished.
echo -n "EndDate: "  >> $REPORT
date --utc --rfc-3339=seconds  >> $REPORT

# Upload final report.
if [ -n "${REPORT_UPLOAD_URL}" ]; then
    echo "Uploading report..."
    upload_file "${REPORT_UPLOAD_URL}" "${REPORT}" "text/plain;charset=utf-8"
fi

# Upload artifacts.
if [ -n "${ARTIFACTS_UPLOAD_URL}" ] &&  [ -d /workspace/artifacts ]; then
    echo "Uploading artifacts..."
    ARTIFACTS_TGZ="/tmp/test_${TESTNAME}_artifacts.tgz"
    tar -czf "${ARTIFACTS_TGZ}" -C /workspace/artifacts .
    upload_file "${ARTIFACTS_UPLOAD_URL}" "${ARTIFACTS_TGZ}" "application/gzip"
fi

echo "Done :-)"

#EOF
