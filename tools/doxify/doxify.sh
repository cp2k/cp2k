#!/bin/bash

set -eo pipefail

ARGS=${ARGS:-} # Set to "-v" for debugging
SCRIPTDIR=$(
  cd $(dirname "$0")
  pwd
) # Pick up full path to scripts from wherever doxify.sh lives

nwarnings=0

for file in "$@"; do
  "${SCRIPTDIR}/is_fypp.py" "${file}" && continue

  # generate temp-file names
  tmp_file=$(mktemp)

  # * Run the fixcomments.pl script. This adds comment blocks to any subroutine/function
  #   definitions that don't have any and checks that existing comments are complete,
  #   fixing these if required.
  # * After adding comments, remove any double comment header lines
  "${SCRIPTDIR}/fixcomments.pl" ${ARGS} "${file}" |
    "${SCRIPTDIR}/remove_extra_comments.pl" ${ARGS} > "${tmp_file}"

  # Copy the final modified source file on top of the original file
  if ! cmp -s "${file}" "${tmp_file}"; then
    cp "${tmp_file}" "${file}"
  fi

  # Remove temp-file
  rm -f "${tmp_file}"

  if grep -e "UNMATCHED_PROCEDURE_ARGUMENT" \
    -e "UNKNOWN_DOXYGEN_COMMENT" \
    -e "UNKNOWN_COMMENT" \
    "${file}"; then
    echo "Found doxify warnings in ${file}"
    ((nwarnings++))
  fi
done

exit ${nwarnings}
