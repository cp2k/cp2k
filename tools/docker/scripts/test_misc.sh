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
    mkdir -p /workspace/artifacts/
    cp test.out /workspace/artifacts/
    echo -e "\nSummary: Test \"${TEST_COMMAND[*]}\" failed."
    echo -e "Status: FAILED\n"
    exit 0
  fi
}

#===============================================================================
cd /opt/cp2k

echo "Using $(python3 --version) and the following packages:"
pip3 freeze
echo ""

# prepare inputs for minimax_to_fortran_source.py
unzip -q -d ./tools/minimax_tools/1_xData 1_xData.zip

run_test ./tools/precommit/format_fortran_test.py
run_test ./tools/minimax_tools/minimax_to_fortran_source.py --check
run_test ./tools/docker/generate_dockerfiles.py --check

run_test mypy --strict ./tools/minimax_tools/minimax_to_fortran_source.py
run_test mypy --strict ./tools/dashboard/generate_dashboard.py
run_test mypy --strict ./tools/dashboard/generate_regtest_survey.py
run_test mypy --strict ./tools/regtesting/optimize_test_dirs.py
run_test mypy --strict ./tools/precommit/precommit.py
run_test mypy --strict ./tools/precommit/check_file_properties.py
run_test mypy --strict ./tools/precommit/format_makefile.py
run_test mypy --strict ./tools/precommit/format_input_file.py
run_test mypy --strict ./tools/docker/generate_dockerfiles.py
run_test mypy --strict ./tools/conventions/analyze_gfortran_ast.py
run_test mypy --strict ./tests/do_regtest.py
run_test mypy --strict ./docs/generate_input_reference.py
run_test mypy --strict ./docs/fix_github_links.py

# TODO: Find a way to test generate_dashboard.py without git repository.
#
# # Test generate_dashboard.py. Running it twice to also execute its caching.
# mkdir -p /workspace/artifacts/dashboard
# for _ in {1..2}; do
#   run_test ./tools/dashboard/generate_dashboard.py \
#     ./tools/dashboard/dashboard.conf \
#     /workspace/artifacts/dashboard/status.pickle \
#     /workspace/artifacts/dashboard/
# done

run_test cmake -DCP2K_ENABLE_CONSISTENCY_CHECKS=ON .

echo ""
echo "Summary: Miscellaneous tests passed"
echo "Status: OK"

#EOF
