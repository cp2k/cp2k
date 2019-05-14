#!/bin/bash -e

# author: Ole Schuett

set -eo pipefail
ulimit -c 0  # Disable core dumps as they can take a very long time to write.

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
    # Don't skip by matching timestamp since Git (previous mirror) does not preserve timestamps
    # and we end up with an outdated tree if the contents of /mnt/cp2k have been modified before
    # the Docker image was generated.
    # Also ignore the checksum since it is safe to assume that most files have not changed (which would
    # mean they have the same filesize) which would mean we would do a read+read(+write) on most files
    # anyway at which point it is simpler to do a read+write (e.g. replace the whole tree).
    rsync --ignore-times                         \
          --delete                               \
          --executability                        \
          --verbose                              \
          --recursive                            \
          --exclude="*~"                         \
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

#EOF
