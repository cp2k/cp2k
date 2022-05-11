# signal trapping

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# shellcheck shell=bash

trap 'error_handler ${LINENO}' ERR
