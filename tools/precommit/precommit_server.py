#!/usr/bin/env python3

# author: Ole Schuett

import os
import logging
from os import path
from time import time
import tempfile
import subprocess
from subprocess import PIPE, STDOUT
from flask import Flask, request, abort

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 1024 * 1024  # 1MB
app.logger.setLevel(logging.INFO)
app.logger.info("CP2K Precommit Server is up and running :-)")


# ======================================================================================
@app.route("/")
def hello():
    return "cp2k precommit server revision: " + os.environ["REVISION"]


# ======================================================================================
@app.route("/black", methods=["POST"])
def black():
    return run_tool(["black"])


# ======================================================================================
@app.route("/spackformat", methods=["POST"])
def spackformat():
    # Run black like https://github.com/spack/spack/blob/develop/pyproject.toml
    return run_tool(["black", "--line-length=99", "--skip-magic-trailing-comma"])


# ======================================================================================
@app.route("/shfmt", methods=["POST"])
def shfmt():
    return run_tool(["shfmt", "-i=2", "-ci", "-sr", "-w"])


# ======================================================================================
@app.route("/shellcheck", methods=["POST"])
def shellcheck():
    return run_tool(["shellcheck"])


# ======================================================================================
@app.route("/mdformat", methods=["POST"])
@app.route("/markdownlint", methods=["POST"])  # for backwards compatibility
def mdformat():
    return run_tool(["mdformat", "--wrap=100"])


# ======================================================================================
@app.route("/clangformat", methods=["POST"])
def clangformat():
    return run_tool(["clang_format_wrapper.sh"])


# ======================================================================================
@app.route("/cmakeformat", methods=["POST"])
def cmakeformat():
    return run_tool(["cmake-format", "-i"])


# ======================================================================================
def run_tool(cmd, timeout=30):
    assert len(request.files) == 1
    orig_fn = list(request.files.keys())[0]
    data_before = request.files[orig_fn].read()
    data_kb = len(data_before) / 1024.0
    fn = path.basename(orig_fn)
    workdir = tempfile.TemporaryDirectory()
    abs_fn = path.join(workdir.name, fn)
    open(abs_fn, "wb").write(data_before)

    t1 = time()
    try:
        p = subprocess.run(
            cmd + [fn], cwd=workdir.name, timeout=timeout, stdout=PIPE, stderr=STDOUT
        )
    except subprocess.TimeoutExpired:
        app.logger.info(f"Timeout of {cmd[0]} on {data_kb:.1f}KB after {timeout}s.")
        return f"Timeout while running {cmd[0]} - please try again.", 504
    t2 = time()
    app.logger.info(f"Ran {cmd[0]} on {data_kb:.1f}KB in {t2-t1:.1f}s.")

    if p.returncode != 0:
        return p.stdout, 422  # Unprocessable Entity
    data_after = open(abs_fn, "rb").read()
    if data_after == data_before:
        return "Not Modified", 304
    return data_after, 200


# ======================================================================================
if __name__ == "__main__":
    app.run()

# EOF
