# CP2K Precommit

The precommit system consists of the following tools that analyze and format the
source code:

- [doxify](../doxify/)
    to add Doxygen templates.
- [prettify](../prettify/)
    to format Fortran files.
- [analyze_src](../conventions/analyze_src.py)
    to check copyright banners and a few other things.
- [ast.parse](https://docs.python.org/3/library/ast.html)
    to check Python syntax.
- [clang-format](https://clang.llvm.org/docs/ClangFormat.html)
    to format C and Cuda files.
- [black](https://github.com/psf/black)
    to format Python scripts.
- [shellcheck](https://github.com/koalaman/shellcheck)
    to analyze Shell scripts.
- [markdownlint](https://github.com/markdownlint/markdownlint)
    to analyze Markdown files.

In contrast to the [CP2K-CI](https://github.com/cp2k/cp2k-ci) these tools
process each file individually, which makes them much more lightweight.

## Install Git Hook

The [precommit.py](./precommit.py) script can be readily installed as git hook:

```shell
ln -fs ../../tools/precommit/precommit.py .git/hooks/pre-commit
```

## Server

Many of the tools listed above require more than just Python. To avoid their
tedious installation a remote server is used by default. It is hosted at
<https://precommit.cp2k.org> via [Cloud Run](https://cloud.google.com/run).

The same server can also be started locally when Docker is available:

```shell
./start_local_server.sh
```

The server can also be installed without Docker by following the steps in the
[Dockerfile](./Dockerfile).

Once the server is up and running it can be used like this:

```shell
$ export CP2K_PRECOMMIT_SERVER="http://127.0.0.1:8080"
$ ./precommit.py
Running precommit checks using 8 workers and server: http://127.0.0.1:8080
Searching for files...
```

## Backups and Cache

Before precommit run an external tool on a file it create a backup copy in
`obj/precommit`. If the tool leaves the file unmodified then the backup copy is
remove afterwards. If the tool modifies the file and precommit was invoked
without `--allow-modifications` then the backup copy is used to restore the
file's original content.

After a successful tool run the file's timestamp is recorded in
`obj/precommit/cache.json` and precommit will skip it in future unless it's
invoked with `--no-cache`.
