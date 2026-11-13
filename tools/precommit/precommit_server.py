#!/usr/bin/env python3

# author: Ole Schuett

import os
import gzip
import asyncio
import logging
import argparse
from os import path
from time import time
from typing import List, Optional
import tempfile
import shutil
import subprocess
from subprocess import PIPE, STDOUT

import aiohttp
from aiohttp import web


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8080)
    args = parser.parse_args()

    app = web.Application(client_max_size=1024 * 1024)  # 1 MB

    # Setup public routes.
    app.router.add_get("/", handle_index)
    app.router.add_post("/black", handle_black)
    app.router.add_post("/spackformat", handle_spackformat)
    app.router.add_post("/shfmt", handle_shfmt)
    app.router.add_post("/shellcheck", handle_shellcheck)
    app.router.add_post("/mdformat", handle_mdformat)
    app.router.add_post("/markdownlint", handle_mdformat)  # for backwards compatibility
    app.router.add_post("/clangformat", handle_clangformat)
    app.router.add_post("/tombi", handle_tombi)
    app.router.add_post("/cmakeformat", handle_cmakeformat)
    app.router.add_post("/cmakelint", handle_cmakelint)
    app.router.add_post("/fortitude", handle_fortitude)
    app.router.add_post("/fortitude_nofix", handle_fortitude_nofix)

    # Start listening for requests.
    print("CP2K Precommit Server is up and running :-)")
    web.run_app(app, port=args.port)


# ======================================================================================
async def handle_index(request: web.Request) -> web.Response:
    text = "cp2k precommit server revision: " + os.environ.get("REVISION", "n/a")
    return web.Response(text=text)


# ======================================================================================
async def handle_black(request: web.Request) -> web.Response:
    return await run_tool(request, ["black"])


# ======================================================================================
async def handle_spackformat(request: web.Request) -> web.Response:
    # Run black like https://github.com/spack/spack/blob/develop/pyproject.toml
    return await run_tool(
        request, ["black", "--line-length=99", "--skip-magic-trailing-comma"]
    )


# ======================================================================================
async def handle_shfmt(request: web.Request) -> web.Response:
    return await run_tool(request, ["shfmt", "-i=2", "-ci", "-sr", "-w"])


# ======================================================================================
async def handle_shellcheck(request: web.Request) -> web.Response:
    return await run_tool(request, ["shellcheck"], timeout=30)


# ======================================================================================
async def handle_mdformat(request: web.Request) -> web.Response:
    return await run_tool(request, ["mdformat", "--wrap=100"])


# ======================================================================================
async def handle_clangformat(request: web.Request) -> web.Response:
    return await run_tool(request, ["clang_format_wrapper.sh"])


# ======================================================================================
async def handle_tombi(request: web.Request) -> web.Response:
    config = f"{os.getcwd()}/tombi.toml"
    return await run_tool(request, ["tombi", "format", "--offline"], config=config)


# ======================================================================================
async def handle_cmakeformat(request: web.Request) -> web.Response:
    return await run_tool(request, ["cmake-format", "-i"], timeout=30)


# ======================================================================================
async def handle_cmakelint(request: web.Request) -> web.Response:
    # Note: the config path must be passed as a single "--flag=value" token.
    # cmake-lint's "-c"/"--config-files" option is nargs="+", so passing it as
    # two separate argv tokens ("-c", path) would greedily swallow the target
    # filename that run_tool() appends as another config file, and silently
    # scan zero files instead.
    config = f"--config-files={os.getcwd()}/cmake-lint.yaml"
    return await run_tool(request, ["cmake-lint", config], timeout=30)


# ======================================================================================
async def handle_fortitude(request: web.Request) -> web.Response:
    config = f"--config-file={os.getcwd()}/fortitude.toml"
    return await run_tool(
        request, ["fortitude", config, "check", "--fix", "--unsafe-fixes"], timeout=5
    )


# ======================================================================================
async def handle_fortitude_nofix(request: web.Request) -> web.Response:
    config = f"--config-file={os.getcwd()}/fortitude.toml"
    # Source formatting happens before expansion; generated lines may be longer.
    return await run_tool(
        request,
        ["fortitude", config, "check", "--extend-per-file-ignores=*:line-too-long"],
        timeout=5,
    )


# ======================================================================================
async def run_tool(
    request: web.Request, cmd: List[str], timeout: int = 3, config: Optional[str] = None
) -> web.Response:

    data = await request.post()
    assert len(data) == 1
    field = list(data.values())[0]
    assert isinstance(field, web.FileField)

    data_before = field.file.read()
    data_kb = len(data_before) / 1024.0
    fn = path.basename(field.filename)
    workdir = tempfile.TemporaryDirectory()
    if config is not None:
        shutil.copyfile(config, path.join(workdir.name, path.basename(config)))
    abs_fn = path.join(workdir.name, fn)
    with open(abs_fn, "wb") as f:
        f.write(data_before)

    t1 = time()
    try:
        p = await asyncio.create_subprocess_exec(
            cmd[0], *cmd[1:], fn, cwd=workdir.name, stdout=PIPE, stderr=STDOUT
        )
        await asyncio.wait_for(p.wait(), timeout=timeout)
    except TimeoutError:
        # Note that GCP Cloud Run returns status 504 when its timeout is reached.
        print(f"Timeout of {cmd[0]} on {data_kb:.1f}KB after {timeout}s.")
        msg = f"Timeout while running {cmd[0]} - please try again."
        return web.Response(text=msg, status=408)  # Request Timeout
    t2 = time()
    print(f"Ran {cmd[0]} on {data_kb:.1f}KB in {t2-t1:.1f}s.")
    if p.returncode != 0:
        return web.Response(body=p.stdout, status=422)  # Unprocessable Entity
    with open(abs_fn, "rb") as f:
        data_after = f.read()
    if data_after == data_before:
        return web.Response(status=304)  # Not Modified
    r = web.Response(body=data_after)
    r.enable_compression()
    return r


# ======================================================================================
if __name__ == "__main__":
    main()

# EOF
