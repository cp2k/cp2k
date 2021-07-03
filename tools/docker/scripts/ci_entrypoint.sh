#!/bin/bash -e

# author: Ole Schuett

set -eo pipefail
ulimit -c 0 # Disable core dumps as they can take a very long time to write.

# Calculate checksums of critical files.
CHECKSUMS=/workspace/checksums.md5
shopt -s nullglob # ignore missing files
md5sum /workspace/cp2k/Makefile \
  /workspace/cp2k/tools/build_utils/* \
  /workspace/cp2k/arch/local* \
  > $CHECKSUMS
shopt -u nullglob

# Get cp2k sources.
if [ -n "${GIT_REF}" ]; then
  echo -e "\n========== Fetching Git Commit =========="
  cd /workspace/cp2k
  git fetch --quiet origin "${GIT_BRANCH}"
  git checkout --quiet "${GIT_REF}"
  git submodule update --init --recursive
  git --no-pager log -1 --pretty='%nCommitSHA: %H%nCommitTime: %ci%nCommitAuthor: %an%nCommitSubject: %s%n'

elif [ -d /mnt/cp2k ]; then
  echo -e "\n========== Copying Changed Files =========="
  # We can't rely on timestamps as they depend on the git checkout time.
  # Hence, we use checksums despite the IO because a full rebuild would be even more costly.
  rsync --checksum \
    --delete \
    --executability \
    --verbose \
    --recursive \
    --exclude="*~" \
    --exclude=".*/" \
    --exclude="*.py[cod]" \
    --exclude="__pycache__" \
    --exclude="/obj/" \
    --exclude="/lib/" \
    --exclude="/exe/" \
    --exclude="/regtesting/" \
    --exclude="/arch/local*" \
    --exclude="/tools/toolchain/build/" \
    --exclude="/tools/toolchain/install/" \
    /mnt/cp2k/ /workspace/cp2k/
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
if ! "$@"; then
  echo -e "\nSummary: Test had non-zero exit status.\nStatus: FAILED"
fi

#EOF
