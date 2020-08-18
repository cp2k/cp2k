#!/usr/bin/python3

# author: Ole Schuett

import ast
import argparse
import json
from time import time, sleep
import re
import os
from os import path
import sys
import concurrent.futures
from subprocess import PIPE, STDOUT
import subprocess
import traceback
from urllib.request import urlopen, Request
from urllib.error import HTTPError
import uuid
from difflib import unified_diff

SCRATCH_DIR = "./obj/precommit"
CACHE_FILE = path.join(SCRATCH_DIR, "cache.json")
SERVER = os.environ.get("CP2K_PRECOMMIT_SERVER", "https://precommit.cp2k.org")


# ======================================================================================
def main():
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
        default=min(32, os.cpu_count() + 4),  # copied from ThreadPoolExecutor
        help="number of parallel workers",
    )
    parser.add_argument(
        "--progressbar-wait",
        type=int,
        default=1,
        help="number seconds in between progressbar updates",
    )
    args = parser.parse_args()

    # Say hello to the server.
    print(
        f"Running precommit checks using {args.num_workers} workers and server: {SERVER}"
    )
    server_hello = urlopen(Request(SERVER + "/"), timeout=10).read().decode("utf8")
    assert server_hello.startswith("cp2k precommit server")

    # Change into base dir and create scratch dir.
    base_dir = path.realpath(path.join(path.dirname(__file__), "../../"))
    os.chdir(base_dir)
    os.makedirs(SCRATCH_DIR, exist_ok=True)

    # Collect eligible files.
    sys.stdout.write("Searching for files...\r")
    sys.stdout.flush()
    filename_pattern = re.compile(r".*\.(F|fypp|c|cu|h|py)$")
    file_list = []
    for root, dirs, files in os.walk("."):
        if root.startswith("./tools/toolchain/build"):
            continue
        if root.startswith("./tools/toolchain/install"):
            continue
        if root.startswith("./tools/prettify/fprettify"):
            continue
        if root.startswith("./tools/build_utils/fypp"):
            continue
        if root.startswith("./tools/autotune_grid"):
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
        file_list += [path.join(root, fn) for fn in files if filename_pattern.match(fn)]

    # Load cache.
    should_load_cache = path.exists(CACHE_FILE) and not args.no_cache
    cache = json.load(open(CACHE_FILE)) if should_load_cache else {}

    # Launch async processing of files.
    futures = {}
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=args.num_workers)
    for fn in file_list:
        if path.getmtime(fn) != cache.get(fn, -1):
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
                    cache[fn] = path.getmtime(fn)
                elif fn not in failed_files:
                    failed_files.add(fn)
                    print_box(fn, str(f.exception()))
        json.dump(cache, open(CACHE_FILE, "w"))
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
def print_box(fn, message):
    print("+" + "-" * 160 + "+")
    print(f"| {fn:^158s} |")
    print("+" + "-" * 160 + "+")
    for line in message.strip().split("\n"):
        print(f"| {line:<158s} |")
    print("+" + "-" * 160 + "+\n\n")


# ======================================================================================
def process_file(fn, allow_modifications):
    orig_content = open(fn, "rb").read()

    if re.match(r".*\.(F|fypp)$", fn):
        run_local_tool("./tools/doxify/doxify.sh", fn)
        run_prettify(fn)
        run_analyze_src(fn)

    elif re.match(r".*\.(c|cu|cpp|h|hpp)$", fn):
        run_remote_tool("clangformat", fn)
        run_analyze_src(fn)

    elif re.match(r".*\.py$", fn):
        ast.parse(orig_content, filename=fn)
        run_remote_tool("black", fn)
        run_analyze_src(fn)

    elif re.match(r".*\.sh$", fn):
        run_remote_tool("shellcheck", fn)

    elif re.match(r".*\.md$", fn):
        run_remote_tool("markdownlint", fn)

    else:
        raise Exception("Unknown file extension: " + fn)

    new_content = open(fn, "rb").read()
    if new_content == orig_content:
        return

    # Deal with a modified file.
    if allow_modifications:
        bak_fn = path.join(SCRATCH_DIR, f"{path.basename(fn)}_{int(time())}.bak")
        open(bak_fn, "wb").write(orig_content)  # make a backup copy
    else:
        open(fn, "wb").write(orig_content)  # restore origin content
        try:
            orig_lines = orig_content.decode("utf8").split("\n")
            new_lines = new_content.decode("utf8").split("\n")
            diff = unified_diff(orig_lines, new_lines, "before", "after", lineterm="")
        except:
            diff = []  #
        raise Exception(f"File modified:\n" + "\n".join(diff))


# ======================================================================================
def run_prettify(fn):
    if fn in ("./src/base/base_uses.f90", "./src/common/util.F"):
        return  # Skipping because of prettify bugs.

    # The prettify tool processes only about 1k lines of code per second.
    # Hence, setting a generous timeout as our largest file has 100k lines.
    run_local_tool(
        "./tools/prettify/prettify.py", "--no-report-errors", fn, timeout=300
    )


# ======================================================================================
def run_analyze_src(fn):
    run_local_tool(
        "./tools/conventions/analyze_src.py",
        "--fail",
        "--suppressions",
        "./tools/conventions/conventions.supp",
        fn,
    )


# ======================================================================================
def run_local_tool(*cmd, timeout=10):
    p = subprocess.run(cmd, timeout=timeout, stdout=PIPE, stderr=STDOUT)
    if p.returncode != 0:
        raise Exception(p.stdout.decode("utf8"))


# ======================================================================================
def run_remote_tool(tool, fn):
    url = f"{SERVER}/{tool}"
    r = http_post(url, fn)
    if r.status == 304:
        pass  # file not modified
    elif r.status == 200:
        open(fn, "wb").write(r.read())
    else:
        raise Exception(r.read().decode("utf8"))  # something went wrong


# ======================================================================================
def http_post(url, fn):
    # This would be so much easier with the Requests library where it'd be a one-liner:
    #    return requests.post(url, files={path.basename(fn): open(fn, "rb").read()})

    boundary = uuid.uuid1().hex
    name = path.basename(fn)
    data = b"".join(
        [
            f"--{boundary}\r\nContent-Disposition: ".encode("utf8"),
            f'form-data; name="{name}"; filename="{name}"\r\n\r\n'.encode("utf8"),
            open(fn, "rb").read(),
            f"\r\n--{boundary}--\r\n".encode("utf8"),
        ]
    )
    headers = {
        "Content-Length": len(data),
        "Content-Type": f"multipart/form-data; boundary={boundary}",
    }
    try:
        return urlopen(Request(url, data=data, headers=headers), timeout=60)
    except HTTPError as err:
        return err  # HTTPError has .status and .read() just like a normal HTTPResponse.


# ======================================================================================
if __name__ == "__main__":
    main()

# EOF
