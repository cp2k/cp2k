#!/usr/bin/env python3

# author: Ole Schuett

import argparse
from pathlib import Path
import re


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description="Sorts the test dirs by duration.")
    parser.add_argument("regtest_report")
    args = parser.parse_args()

    report = Path(args.regtest_report).read_text(encoding="utf8")
    match = re.search(">>> (.*?)/UNIT\n", report)
    assert match
    report_basedir = match.group(1)

    timings = {}
    for line in re.findall("<<< (.*?)\n", report):
        parts = line.split()
        name = parts[0][len(report_basedir) + 1 :]
        timings[name] = float(parts[6])

    cp2k_root = Path(__file__).resolve().parent.parent.parent
    test_dirs_fn = cp2k_root / "tests" / "TEST_DIRS"
    header_lines = []
    test_dir_lines = []
    for line in test_dirs_fn.read_text(encoding="utf8").split("\n"):
        if not line.strip():
            pass
        elif line.startswith("#"):
            header_lines.append(line)
        else:
            test_dir_lines.append(line)

    def sort_key(line: str) -> float:
        return timings.get(line.split()[0], float("inf"))

    test_dir_lines.sort(key=sort_key, reverse=True)
    output = "\n".join(header_lines + test_dir_lines) + "\n"
    test_dirs_fn.write_text(output, encoding="utf8")

    print(f"Updated {test_dirs_fn}")


# ======================================================================================
if __name__ == "__main__":
    main()

# EOF
