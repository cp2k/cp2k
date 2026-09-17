"""Run the existing Wannier90 regression inputs with their unchanged matchers."""

import argparse
from pathlib import Path
import shutil
import sys
import tomllib

from validate_cp2k import ROOT, run

sys.path.insert(0, str(ROOT / "tests"))
from matchers import run_matcher


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("workdir", type=Path)
    args = parser.parse_args()
    source = ROOT / "tests/QS/regtest-kp-1"
    directory = args.workdir.resolve()
    shutil.copytree(source, directory)
    specs = tomllib.loads((source / "TEST_FILES.toml").read_text())
    count = 0
    for name, checks in specs.items():
        if "wannier90" not in name:
            continue
        run(args.binary.resolve(), directory, (source / name).read_text())
        output = (directory / "run.log").read_text()
        (directory / "run.log").rename(directory / (name + ".log"))
        for check in checks:
            result = run_matcher(output, **check)
            assert result.status == "OK", (name, check, result)
        count += 1
        print(f"Passed {name}", flush=True)
    assert count > 0
    print(f"All {count} legacy Wannier90 regressions passed")


if __name__ == "__main__":
    main()
