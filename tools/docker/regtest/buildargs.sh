
# author: Ole Schuett

BUILDARGS="$BUILDARGS --build-arg TESTNAME=${TESTNAME}"

if [[ "${TESTNAME}" == "minimal" ]]; then
  BUILDARGS="$BUILDARGS --build-arg ARCH=minimal --build-arg VERSION=sdbg"

elif [[ "${TESTNAME}" == "farming" ]]; then
  BUILDARGS="$BUILDARGS --build-arg ARCH=local --build-arg VERSION=popt"

elif [[ "${TESTNAME}" == "coverage-sdbg" ]]; then
  BUILDARGS="$BUILDARGS --build-arg ARCH=local_coverage --build-arg VERSION=sdbg"

elif [[ "${TESTNAME}" == "coverage-pdbg" ]]; then
  BUILDARGS="$BUILDARGS --build-arg ARCH=local_coverage --build-arg VERSION=pdbg"

elif [[ "${TESTNAME}" =~ ^(sopt|ssmpt|sdbg|popt|psmp|pdbg)$ ]]; then
  BUILDARGS="$BUILDARGS --build-arg ARCH=local --build-arg VERSION=${TESTNAME}"

else
  echo "Unkown test name: ${TESTNAME}"
  exit 1

fi

#EOF
