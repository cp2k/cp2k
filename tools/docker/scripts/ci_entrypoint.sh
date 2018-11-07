#!/bin/bash -e

# author: Ole Schuett

REPORT=/workspace/report.txt
echo -n "StartDate: " > $REPORT
date --utc --rfc-3339=seconds >> $REPORT

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
    #TODO: Exclusion rules could also come from .gitigore
    rsync --exclude="*~"                      \
          --exclude=".*/"                     \
          --exclude="*.pyc"                   \
          --exclude=/obj/                     \
          --exclude=/lib/                     \
          --exclude=/exe/                     \
          --exclude=/regtesting/              \
          --exclude=/tools/toolchain/         \
          --executability                     \
          --ignore-times                      \
          --update                            \
          --verbose                           \
          --recursive                         \
          --checksum                          \
          /mnt/cp2k/  /workspace/cp2k/
else
    echo "Neighter GIT_REF nor /mnt/cp2k found - aborting."
    exit 255
fi


# Update toolchain.
if [ -d /opt/cp2k-toolchain ]; then
    echo -e "\n========== Updating Toolchain =========="
    cd /opt/cp2k-toolchain
    rsync --exclude="*~"          \
          --exclude=".*/"         \
          --exclude="*.pyc"       \
          --exclude="/build"      \
          --exclude="/install"    \
          --executability         \
          --ignore-times          \
          --update                \
          --verbose               \
          --recursive             \
          --checksum              \
          /workspace/cp2k/tools/toolchain/  /opt/cp2k-toolchain/

    # shellcheck disable=SC1091
    source /opt/cp2k-toolchain/install/setup
    # shellcheck disable=SC2086
    ./install_cp2k_toolchain.sh ${CP2K_TOOLCHAIN_OPTIONS}
fi

echo -e "\n========== Running Test =========="
cd /workspace
"$@" 2>&1 | tee -a $REPORT &  # Launch in the background.
TEST_PID=$!
trap '[ -z $TEST_PID] || kill $TEST_PID' SIGHUP SIGINT SIGQUIT SIGTERM

# Upload preliminary report every 30s while test is running.
while [ -e /proc/${TEST_PID} ]; do
    sleep 1
    count=$(( (count + 1) % 30 ))
    if (( count == 1 )) && [ -n "${REPORT_UPLOAD_URL}" ]; then
        wget -O- --method=PUT --header="content-type: text/plain;charset=utf-8" --header="cache-control: no-cache" --body-file="${REPORT}" "${REPORT_UPLOAD_URL}"
    fi
done

# Test has finished.
echo -n "EndDate: "  >> $REPORT
date --utc --rfc-3339=seconds  >> $REPORT

# Upload final report.
if [ -n "${REPORT_UPLOAD_URL}" ]; then
    echo "Uploading report..."
    wget -O- --method=PUT --header="content-type: text/plain;charset=utf-8" --header="cache-control: no-cache" --body-file="${REPORT}" "${REPORT_UPLOAD_URL}"
fi

# Upload artifacts.
if [ -n "${ARTIFACTS_UPLOAD_URL}" ] &&  [ -d /workspace/artifacts ]; then
    echo "Uploading artifacts..."
    ARTIFACTS_TGZ="/tmp/test_${TESTNAME}_artifacts.tgz"
    tar -czf "${ARTIFACTS_TGZ}" -C /workspace/artifacts .
    wget -O- --method=PUT --header="content-type: application/gzip" --header="cache-control: no-cache" --body-file="${ARTIFACTS_TGZ}" "${ARTIFACTS_UPLOAD_URL}"
fi

echo "Done :-)"

#EOF
