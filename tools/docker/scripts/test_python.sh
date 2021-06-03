#!/bin/bash -e

# author: Ole Schuett

function run_test {
  TEST_COMMAND=("$@")
  echo -en "Running \"${TEST_COMMAND[*]}\"... "
  if "${TEST_COMMAND[@]}" &> test.out; then
    echo "done."
  else
    echo -e "failed.\n\n"
    tail -n 100 test.out
    echo -e "\nSummary: Test \"${TEST_COMMAND[*]}\" failed."
    echo -e "Status: FAILED\n"
    exit 0
  fi
}

#===============================================================================
cd /workspace/cp2k

echo "Using $(python3 --version) and $(mypy --version)."
echo ""

run_test ./tools/prettify/prettify_test.py
run_test mypy --strict ./tools/dashboard/generate_dashboard.py
run_test mypy --strict ./tools/regtesting/do_regtest.py
run_test mypy --strict ./tools/regtesting/optimize_test_dirs.py

# Test generate_dashboard.py. Running it twice to also execute its caching.
mkdir -p /workspace/artifacts/dashboard
for _ in {1..2}; do
  run_test ./tools/dashboard/generate_dashboard.py \
    ./tools/dashboard/dashboard.conf \
    /workspace/artifacts/dashboard/status.pickle \
    /workspace/artifacts/dashboard/
done

echo ""
echo "Summary: Python tests passed"
echo "Status: OK"

#EOF
