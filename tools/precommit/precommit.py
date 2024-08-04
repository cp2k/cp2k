#!/usr/bin/env python3

# author: Ole Schuett

import ast
import argparse
import json
from time import time, sleep
import re
import os
import sys
import concurrent.futures
from subprocess import PIPE, STDOUT
import subprocess
import traceback
from urllib.request import urlopen, Request
from http.client import HTTPResponse
from urllib.error import HTTPError
import uuid
from difflib import unified_diff
from pathlib import Path
import shutil
from typing import cast, Iterator

SCRATCH_DIR = Path("./obj/precommit")
CACHE_FILE = SCRATCH_DIR / "cache.json"
SERVER = os.environ.get("CP2K_PRECOMMIT_SERVER", "https://precommit.cp2k.org")


# ======================================================================================
def main() -> None:
    # Parse command line arguments.
    parser = argparse.ArgumentParser(
        description="Check source code for formatting and linter problems."
    )
    parser.add_argument(
        "-a",
        "--no-cache",
        action="store_true",
        help="ignore the cache and check all files",
    )
    parser.add_argument(
        "-m",
        "--allow-modifications",
        action="store_true",
        help="allow the tools to modify files",
    )
    parser.add_argument(
        "-j",
        "--num_workers",
        type=int,
        default=min(16, (os.cpu_count() or 0) + 2),
        help="number of parallel workers",
    )
    parser.add_argument(
        "--progressbar-wait",
        metavar="SECONDS",
        type=int,
        default=1,
        help="number seconds in between progressbar updates",
    )
    parser.add_argument("files", metavar="FILE", nargs="*", help="files to process")
    args = parser.parse_args()

    # Say hello to the server.
    print(
        f"Running precommit checks using {args.num_workers} workers and server: {SERVER}"
    )
    server_hello = urlopen(Request(SERVER + "/"), timeout=10).read().decode("utf8")
    assert server_hello.startswith("cp2k precommit server")

    # Store candidate before changing base directory and creating scratch dir.
    file_list = [os.path.abspath(fn) for fn in args.files]
    base_dir = Path(__file__).resolve().parent.parent.parent
    os.chdir(base_dir)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)

    # Collect candidate files.
    if not file_list:
        sys.stdout.write("Searching for files...\r")
        sys.stdout.flush()
        for root, dirs, files in os.walk("."):
            if root.startswith("./tools/toolchain/build"):
                continue
            if root.startswith("./tools/toolchain/install"):
                continue
            if root.startswith("./tools/prettify/fprettify"):
                continue
            if root.startswith("./tools/precommit/fprettify"):
                continue
            if root.startswith("./tools/build_utils/fypp"):
                continue
            if root.startswith("./tools/autotools"):
                continue
            if root.startswith("./tools/minimax_tools/1_xData"):
                continue
            if root.startswith("./tools/fedora"):
                continue
            if root.startswith("./data/DFTB/scc"):
                continue
            if root.startswith("./arch"):
                continue
            if root.startswith("./doxygen"):
                continue
            if root.startswith("./docs/_build"):
                continue
            if root.startswith("./docs/CP2K_INPUT"):
                continue
            if root.startswith("./exts"):
                continue
            if root.startswith("./obj"):
                continue
            if root.startswith("./lib"):
                continue
            if root.startswith("./exe"):
                continue
            if root.startswith("./regtesting"):
                continue
            if root.startswith("./.git"):
                continue
            if "/.mypy_cache/" in root:
                continue
            file_list += [os.path.join(root, fn) for fn in files]

    # Filter symlinks, backup copies, logs, and hidden files.
    file_list = [fn for fn in file_list if not os.path.islink(fn)]
    file_list = [fn for fn in file_list if not fn[-1] in ("~", "#")]
    file_list = [fn for fn in file_list if not fn.endswith(".log")]
    file_list = [fn for fn in file_list if not os.path.basename(fn).startswith(".")]

    # Sort files by size as larger ones will take longer to process.
    file_list.sort(reverse=True, key=lambda fn: os.path.getsize(fn))

    # Load cache.
    should_load_cache = CACHE_FILE.exists() and not args.no_cache
    cache = json.loads(CACHE_FILE.read_text()) if should_load_cache else {}

    # Launch async processing of files.
    futures = {}
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=args.num_workers)
    for fn in file_list:
        if os.path.getmtime(fn) != cache.get(fn, -1):
            futures[fn] = executor.submit(process_file, fn, args.allow_modifications)
    num_skipped = len(file_list) - len(futures)

    # Continuously update progressbar, save cache file, and print errors.
    failed_files = set()
    while True:
        num_done = num_skipped
        for fn, f in futures.items():
            if f.done():
                num_done += 1
                if not f.exception():
                    cache[fn] = os.path.getmtime(fn)
                elif fn not in failed_files:
                    failed_files.add(fn)
                    print_box(fn, str(f.exception()))
        CACHE_FILE.write_text(json.dumps(cache))
        progressbar = "=" * int(60 * num_done / len(file_list))
        sys.stdout.write(
            f"[{progressbar:60s}] {num_done} / {len(file_list)} files processed\r"
        )
        sys.stdout.flush()
        if num_done == len(file_list) or len(failed_files) >= 10:
            executor.shutdown(wait=False)
            break
        sleep(args.progressbar_wait)

    # Print final message.
    print(
        f"Summary: Found {len(file_list)}, "
        f"skipped {num_skipped}, "
        f"checked {num_done - num_skipped}, "
        f"and failed {len(failed_files)} files." + (" " * 50)
    )
    print("Status: " + ("FAILED" if failed_files else "OK"))
    sys.exit(len(failed_files))


# ======================================================================================
def print_box(fn: str, message: str) -> None:
    print("+" + "-" * 160 + "+")
    print(f"| {fn:^158s} |")
    print("+" + "-" * 160 + "+")
    for line in message.strip().split("\n"):
        print(f"| {line:<158s} |")
    print("+" + "-" * 160 + "+\n\n")


# ======================================================================================
def process_file(fn: str, allow_modifications: bool) -> None:
    # Make a backup copy.
    orig_content = Path(fn).read_bytes()
    bak_fn = SCRATCH_DIR / f"{Path(fn).name}_{time()}.bak"
    shutil.copy2(fn, bak_fn)

    if re.match(r".*\.(F|fypp)$", fn):
        run_local_tool("./tools/doxify/doxify.sh", fn)
        run_format_fortran(fn)

    if re.match(r".*\.(c|cu|cl|h)$", fn):
        run_remote_tool("clangformat", fn)

    if re.match(r".*\.(cc|cpp|cxx|hcc|hpp|hxx)$", fn):
        if fn.endswith("/torch_c_api.cpp"):
            # Begrudgingly tolerated because PyTorch has no C API.
            run_remote_tool("clangformat", fn)
        else:
            raise Exception(f"C++ is not supported.")

    if re.match(r"(.*/PACKAGE)|(.*\.py)$", fn):
        ast.parse(orig_content, filename=fn)
        run_remote_tool("black", fn)

    if re.match(r".*\.sh$", fn):
        run_remote_tool("shfmt", fn)
        run_remote_tool("shellcheck", fn)

    if re.match(r".*\.md$", fn):
        run_remote_tool("mdformat", fn)

    if re.match(r"(.*/CMakeLists.txt)|(.*\.cmake)$", fn):
        run_remote_tool("cmakeformat", fn)

    if re.match(r"./data/.*POTENTIALS?$", fn):
        check_data_files()

    if re.match(r".*/Makefile", fn):
        run_local_tool("./tools/precommit/format_makefile.py", fn)

    if re.match(r".*\.inp$", fn):
        run_local_tool("./tools/precommit/format_input_file.py", fn)

    run_check_file_properties(fn)

    new_content = Path(fn).read_bytes()
    if new_content == orig_content:
        bak_fn.unlink()  # remove backup

    elif not allow_modifications:
        bak_fn.replace(fn)  # restore origin content
        diff: Iterator[str]
        try:
            orig_lines = orig_content.decode("utf8").split("\n")
            new_lines = new_content.decode("utf8").split("\n")
            diff = unified_diff(orig_lines, new_lines, "before", "after", lineterm="")
        except Exception:
            diff = iter([])
        raise Exception(f"File modified:\n" + "\n".join(diff))


# ======================================================================================
def run_format_fortran(fn: str) -> None:
    if fn in ("./src/base/base_uses.f90", "./src/common/util.F"):
        return  # Skipping because of prettify bugs.

    # The prettify tool processes only about 1k lines of code per second.
    # Hence, setting a generous timeout as our largest file has 100k lines.
    run_local_tool("./tools/precommit/format_fortran.py", fn, timeout=600)


# ======================================================================================
def run_check_file_properties(fn: str) -> None:
    run_local_tool("./tools/precommit/check_file_properties.py", fn)


# ======================================================================================
def check_data_files() -> None:
    potential_files = [f"./data/{x}_POTENTIALS" for x in ("GTH", "HF", "NLCC", "ALL")]
    expected_content = "".join([Path(fn).read_text() for fn in potential_files])
    if Path("./data/POTENTIAL").read_text() != expected_content:
        raise Exception("The data files are out of sync - please run `make data`.")


# ======================================================================================
def run_local_tool(*cmd: str, timeout: int = 20) -> None:
    p = subprocess.run(cmd, timeout=timeout, stdout=PIPE, stderr=STDOUT)
    if p.returncode != 0:
        raise Exception(p.stdout.decode("utf8"))


# ======================================================================================
def run_remote_tool(tool: str, fn: str) -> None:
    url = f"{SERVER}/{tool}"
    r = http_post(url, fn)
    if r.status == 304:
        pass  # file not modified
    elif r.status == 200:
        Path(fn).write_bytes(r.read())
    else:
        raise Exception(r.read().decode("utf8"))  # something went wrong


# ======================================================================================
def http_post(url: str, fn: str) -> HTTPResponse:
    # This would be so much easier with the Requests library where it'd be a one-liner:
    #    return requests.post(url, files={Path(fn).name: Path(fn).read_bytes()})

    boundary = uuid.uuid1().hex
    name = Path(fn).name
    data = b"".join(
        [
            f"--{boundary}\r\nContent-Disposition: ".encode("utf8"),
            f'form-data; name="{name}"; filename="{name}"\r\n\r\n'.encode("utf8"),
            Path(fn).read_bytes(),
            f"\r\n--{boundary}--\r\n".encode("utf8"),
        ]
    )
    headers = {
        "Content-Length": f"{len(data)}",
        "Content-Type": f"multipart/form-data; boundary={boundary}",
    }
    try:
        response = urlopen(Request(url, data=data, headers=headers), timeout=60)
        return cast(HTTPResponse, response)
    except HTTPError as err:
        # HTTPError has .status and .read() just like a normal HTTPResponse.
        return cast(HTTPResponse, err)


# ======================================================================================
if __name__ == "__main__":
    main()

# EOF
