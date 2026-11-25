#!/usr/bin/env python3
"""Reserve local physical cores and exec Open MPI with a native rankfile."""

import argparse
from contextlib import contextmanager, ExitStack
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import stat
import subprocess
import sys


def cpu_list(text):
    cpus = set()
    for part in text.strip().split(","):
        bounds = [int(i) for i in part.split("-")]
        if len(bounds) == 1:
            bounds *= 2
        if len(bounds) != 2 or not 0 <= bounds[0] <= bounds[1]:
            raise ValueError(f"Invalid CPU list: {text!r}")
        cpus.update(range(bounds[0], bounds[1] + 1))
    return cpus


def topology(root=Path("/sys/devices/system/cpu")):
    online = cpu_list((root / "online").read_text())
    groups = {}
    for cpu in online:
        groups[cpu] = tuple(
            sorted(
                cpu_list((root / f"cpu{cpu}/topology/thread_siblings_list").read_text())
            )
        )
    for cpu, siblings in groups.items():
        if cpu not in siblings or any(
            groups[s] != siblings for s in set(siblings) & online
        ):
            raise ValueError(f"Inconsistent sysfs sibling group at CPU {cpu}")
    return online, sorted(set(groups.values()))


def open_lock(path):
    fd = os.open(path, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    info = os.fstat(fd)
    if (
        not stat.S_ISREG(info.st_mode)
        or info.st_uid != os.getuid()
        or info.st_nlink != 1
        or stat.S_IMODE(info.st_mode) != 0o600
    ):
        os.close(fd)
        raise ValueError(f"Unsafe pool lock: {path}")
    return os.fdopen(fd, "r+b")


@contextmanager
def reserve(directory, groups, allowed, count):
    directory = Path(directory)
    directory.mkdir(mode=0o700, exist_ok=True)
    info = directory.lstat()
    if (
        not stat.S_ISDIR(info.st_mode)
        or info.st_uid != os.getuid()
        or stat.S_IMODE(info.st_mode) != 0o700
    ):
        raise ValueError(f"Pool must be owned by this user with mode 0700: {directory}")
    free, locked = [], []
    # The inner stack rolls back partial claims before releasing the mutex.
    with open_lock(directory / "allocation.lock") as mutex, ExitStack() as claims:
        fcntl.flock(mutex, fcntl.LOCK_EX)
        for group in [g for g in groups if set(g) & allowed]:
            path = directory / ("core-" + "-".join(map(str, group)) + ".lock")
            lock = claims.enter_context(open_lock(path))
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                lock.close()
                continue
            locked.append(lock)
            free.append(group)
        if len(free) < count:
            raise ValueError(
                f"Need {count} physical cores, only {len(free)} available after reservations"
            )
        for lock in locked[count:]:
            lock.close()
        held = claims.pop_all()
    # Persistent lock inodes, inherited by mpirun; never unlink or LOCK_UN them.
    with held:
        yield (
            free,
            [min(set(g) & allowed) for g in free[:count]],
            locked[:count],
        )


def positive_integer(value):
    if not re.fullmatch(r"[0-9]+", value) or int(value) == 0:
        raise argparse.ArgumentTypeError(
            f"expected one positive integer, got {value!r}; thread lists are unsupported"
        )
    return int(value)


def argument_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    for short, name, variable in (
        ("-n", "ranks", "MPI_RANKS"),
        ("-t", "threads", "OMP_NUM_THREADS"),
    ):
        # String defaults are converted only if the CLI did not override them.
        parser.add_argument(
            short,
            f"--{name}",
            type=positive_integer,
            default=os.environ.get(variable, "1"),
            metavar="N",
            help=f"{name}: CLI > {variable} > 1 (positive integer)",
        )
    parser.add_argument("--mpirun", default="mpirun", help="Open MPI 5 executable")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print placement, then release the reservation",
    )
    parser.add_argument(
        "command", nargs=argparse.REMAINDER, help="-- cp2k.psmp -i input.inp"
    )
    return parser


def make_rankfile(layout, hwloc_calc):
    lines = []
    for rank, cpus in enumerate(layout):
        # Rankfiles use hwloc logical PUs, never Linux OS CPU numbers.
        slots = subprocess.check_output(
            [
                hwloc_calc,
                "--physical-input",
                "--logical-output",
                "--intersect",
                "PU",
                *[f"pu:{cpu}" for cpu in cpus],
            ],
            text=True,
        ).strip()
        if not re.fullmatch(r"[0-9]+(,[0-9]+)*", slots) or len(
            set(slots.split(","))
        ) != len(cpus):
            raise ValueError(f"hwloc could not map OS CPUs {cpus}: {slots!r}")
        lines.append(f"rank {rank}=localhost slot={slots}\n")
    return "".join(lines)


def main():
    parser = argument_parser()
    args = parser.parse_args()
    if args.command[:1] == ["--"]:
        args.command.pop(0)
    if not args.command:
        parser.error("an application is required")
    if sys.platform != "linux" or os.getuid() == 0:
        parser.error("run as an ordinary user on Linux")
    if any(
        k in os.environ
        for k in ("SLURM_JOB_ID", "PBS_JOBID", "LSB_JOBID", "OMPI_COMM_WORLD_RANK")
    ):
        parser.error("use this launcher outside schedulers and MPI jobs")
    mpirun = shutil.which(args.mpirun)
    if not mpirun or not shutil.which(args.command[0]):
        parser.error("mpirun and the application must be executable")
    version = subprocess.check_output([mpirun, "--version"], text=True)
    if not re.search(r"\(Open MPI\) 5\.", version):
        parser.error("this rankfile launcher requires Open MPI 5")
    calc = shutil.which(str(Path(mpirun).with_name("hwloc-calc"))) or shutil.which(
        "hwloc-calc"
    )
    if not calc:
        parser.error("hwloc-calc from the MPI installation must be available")
    online, groups = topology()
    allowed = online & os.sched_getaffinity(0)
    with reserve(
        f"/tmp/cp2k-local-{os.getuid()}",
        groups,
        allowed,
        args.ranks * args.threads,
    ) as (free, selected, locks):
        layout = [
            selected[i : i + args.threads]
            for i in range(0, len(selected), args.threads)
        ]
        contents = make_rankfile(layout, calc)
        # An inherited anonymous file needs no persistent rankfile or cleanup daemon.
        with os.fdopen(os.memfd_create("cp2k-rankfile", 0), "w+") as rankfile:
            rankfile.write(contents)
            rankfile.flush()
            command = [
                mpirun,
                "--host",
                f"localhost:{args.ranks}",
                "-n",
                str(args.ranks),
                "--map-by",
                f"rankfile:file=/proc/self/fd/{rankfile.fileno()}:hwtcpus",
                "--bind-to",
                "hwthread",
                *args.command,
            ]
            omp = dict(
                OMP_NUM_THREADS=str(args.threads),
                OMP_PLACES="threads",
                OMP_PROC_BIND="close",
                OMP_DYNAMIC="FALSE",
            )
            plan = dict(
                online=sorted(online),
                allowed=sorted(allowed),
                physical_cores=groups,
                available_cores=free,
                rank_cpus=layout,
            )
            print("cp2k-local: plan " + json.dumps(plan), file=sys.stderr)
            print("cp2k-local: rankfile\n" + contents, end="", file=sys.stderr)
            print(
                "cp2k-local: command "
                + shlex.join(["env", *[f"{k}={v}" for k, v in omp.items()], *command]),
                file=sys.stderr,
                flush=True,
            )
            if not args.dry_run:
                for lock in locks:
                    os.set_inheritable(lock.fileno(), True)
                os.execvpe(mpirun, command, dict(os.environ, **omp))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f"cp2k-local: {exc}", file=sys.stderr)
        sys.exit(2)
