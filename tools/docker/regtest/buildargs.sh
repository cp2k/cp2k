
# author: Ole Schuett

BUILDARGS="$BUILDARGS --build-arg TESTNAME=${TESTNAME}"

if [[ "${TESTNAME}" =~ ^(sopt|ssmp|sdbg|popt|psmp|pdbg)$ ]]; then
  BUILDARGS="$BUILDARGS --build-arg ARCH=local --build-arg VERSION=${TESTNAME}"

elif [[ "${TESTNAME}" == "farming" ]]; then
  BUILDARGS="$BUILDARGS --build-arg ARCH=local --build-arg VERSION=popt"

elif [[ "${TESTNAME}" == "minimal" ]]; then
  BUILDARGS="$BUILDARGS --build-arg ARCH=minimal --build-arg VERSION=sdbg"

elif [[ "${TESTNAME}" =~ ^coverage-(sdbg|pdbg)$ ]]; then
  BUILDARGS="$BUILDARGS --build-arg ARCH=local_coverage --build-arg VERSION=${TESTNAME:9}"

elif [[ "${TESTNAME}" =~ ^valgrind-(sdbg|sopt|pdbg|popt)$ ]]; then
  BUILDARGS="$BUILDARGS --build-arg ARCH=local_valgrind --build-arg VERSION=${TESTNAME:9}"

else
  echo "Unkown test name: ${TESTNAME}"
  exit 1

fi

#EOF
