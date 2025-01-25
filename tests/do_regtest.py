#!/usr/bin/env python3

# author: Ole Schuett

from asyncio import Semaphore, Task
from asyncio.subprocess import DEVNULL, PIPE, STDOUT, Process
from datetime import datetime
from pathlib import Path
from typing import Any, Coroutine, Dict, List, Optional, TextIO, Tuple, Union
from statistics import mean, stdev
import argparse
import asyncio
import math
import os
import re
import shutil
import subprocess
import sys
import time

try:
    from typing import Literal  # not available before Python 3.8

    TestStatus = Literal[
        "OK", "WRONG RESULT", "RUNTIME FAIL", "TIMED OUT", "HUGE OUTPUT"
    ]
except ImportError:
    TestStatus = str  # type: ignore

# Some tests do not work with --keepalive (which is generally considered a bug).
KEEPALIVE_SKIP_DIRS = [
    "TMC/regtest",
    "TMC/regtest_ana_on_the_fly",
    "TMC/regtest_ana_post_proc",
    "LIBTEST/libvori",
    "LIBTEST/libbqb",
    "QS/regtest-ec",
    "QS/regtest-gapw",
    "QS/regtest-gapw_xc",
    "QS/regtest-gpw-1",
    "QS/regtest-mp2-grad-1",
    "QS/regtest-mp2-grad-2",
]


# ======================================================================================
async def main() -> None:
    parser = argparse.ArgumentParser(description="Runs CP2K regression test suite.")
    parser.add_argument("--mpiranks", type=int, default=2)
    parser.add_argument("--ompthreads", type=int)
    parser.add_argument("--maxtasks", type=int, default=os.cpu_count())
    parser.add_argument("--num_gpus", type=int, default=0)
    parser.add_argument("--timeout", type=int, default=400)
    parser.add_argument("--maxerrors", type=int, default=50)
    help = "Template for launching MPI jobs, {N} is replaced by number of processors."
    parser.add_argument("--mpiexec", default="mpiexec -n {N} --bind-to none", help=help)
    help = "Runs only the first test of each directory."
    parser.add_argument("--smoketest", dest="smoketest", action="store_true", help=help)
    help = "Runs tests under Valgrind memcheck. Best used together with --keepalive."
    parser.add_argument("--valgrind", action="store_true", help=help)
    help = "Use a persistent cp2k-shell process to reduce startup time."
    parser.add_argument("--keepalive", dest="keepalive", action="store_true", help=help)
    help = "Flag slow tests in the final summary and status report."
    parser.add_argument("--flagslow", dest="flagslow", action="store_true", help=help)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--restrictdir", action="append")
    parser.add_argument("--skipdir", action="append")
    parser.add_argument("--workbasedir", type=Path)
    parser.add_argument("--skip_unittests", action="store_true")
    parser.add_argument("--skip_regtests", action="store_true")
    parser.add_argument("binary_dir", type=Path)
    parser.add_argument("version")
    cfg = Config(parser.parse_args())

    print("*************************** Testing started ****************************")
    start_time = time.perf_counter()

    # Query CP2K binary for feature flags.
    version_bytes, _ = await (await cfg.launch_exe("cp2k", "--version")).communicate()
    version_output = version_bytes.decode("utf8", errors="replace")
    flags_line = re.search(r" cp2kflags:(.*)\n", version_output)
    if not flags_line:
        print(version_output + "\nCould not parse feature flags.")
        sys.exit(1)
    else:
        flags = flags_line.group(1).split()

    print("\n----------------------------- Settings ---------------------------------")
    print(f"MPI ranks:      {cfg.mpiranks}")
    print(f"OpenMP threads: {cfg.ompthreads}")
    print(f"GPU devices:    {cfg.num_gpus}")
    print(f"Workers:        {cfg.num_workers}")
    print(f"Timeout [s]:    {cfg.timeout}")
    print(f"Work base dir:  {cfg.work_base_dir}")
    print(f"MPI exec:       {cfg.mpiexec}")
    print(f"Smoke test:     {cfg.smoketest}")
    print(f"Valgrind:       {cfg.valgrind}")
    print(f"Keepalive:      {cfg.keepalive}")
    print(f"Flag slow:      {cfg.flag_slow}")
    print(f"Debug:          {cfg.debug}")
    print(f"Binary dir:     {cfg.binary_dir}")
    print(f"VERSION:        {cfg.version}")
    print(f"Flags:          " + ",".join(flags))

    # Have to copy everything upfront because the test dirs are not self-contained.
    print("------------------------------------------------------------------------")
    print("Copying test files ...", end="")
    shutil.copytree(cfg.cp2k_root / "tests", cfg.work_base_dir)
    print(" done")

    batches: List[Batch] = []

    # Read UNIT_TESTS.
    unit_tests_fn = cfg.cp2k_root / "tests" / "UNIT_TESTS"
    for line in unit_tests_fn.read_text(encoding="utf8").split("\n"):
        line = line.split("#", 1)[0].strip()
        if line:
            batch = Batch(f"UNIT/{line}", cfg)
            batch.workdir.mkdir(parents=True)
            batch.unittests.append(Unittest(line.split()[0], batch.workdir))
            batches.append(batch)

    # Read TEST_TYPES.
    test_types_fn = cfg.cp2k_root / "tests" / "TEST_TYPES"
    test_types: List[Optional[TestType]] = [None]  # test type zero
    lines = test_types_fn.read_text(encoding="utf8").split("\n")
    test_types += [TestType(l) for l in lines[1 : int(lines[0]) + 1]]

    # Read TEST_DIRS.
    test_dirs_fn = cfg.cp2k_root / "tests" / "TEST_DIRS"
    for line in test_dirs_fn.read_text(encoding="utf8").split("\n"):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        batch = Batch(line, cfg)

        # Read TEST_FILES.
        test_files_fn = Path(batch.src_dir / "TEST_FILES")
        for line in test_files_fn.read_text(encoding="utf8").split("\n"):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            batch.regtests.append(Regtest(line, test_types, batch.workdir))
            if cfg.smoketest:
                break  # run only one test per directory
        batches.append(batch)

    # Create async tasks.
    tasks: List[Task[BatchResult]] = []
    num_restrictdirs = num_skipdirs = 0
    for batch in batches:
        if not batch.requirements_satisfied(flags, cfg.mpiranks):
            print(f"Skipping {batch.name} because its requirements are not satisfied.")
        elif not any(re.fullmatch(p, batch.name) for p in cfg.restrictdirs):
            num_restrictdirs += 1
        elif any(re.fullmatch(p, batch.name) for p in cfg.skipdirs):
            num_skipdirs += 1
        else:
            # Need to keep a reference of the task to protect it from garbage collection.
            # asyncio.create_task not available before Python 3.7
            tasks.append(asyncio.get_event_loop().create_task(run_batch(batch, cfg)))

    if num_restrictdirs:
        print(f"Skipping {num_restrictdirs} test directories because of --restrictdir.")
    if num_skipdirs:
        print(f"Skipping {num_skipdirs} test directories because of --skipdir.")
    if not tasks:
        print("\nNo test directories selected, check --restrictdir filter.")
        sys.exit(1)
    print(f"Launched {len(tasks)} test directories and {cfg.num_workers} worker...\n")
    sys.stdout.flush()

    # Wait for tasks to finish and print their results.
    all_results: List[TestResult] = []
    with open(cfg.error_summary, "wt", encoding="utf8", errors="replace") as err_fh:
        for num_done, task in enumerate(asyncio.as_completed(tasks)):
            batch_result = await task
            all_results += batch_result.results
            print(f">>> {batch_result.batch.workdir}")
            print("\n".join(str(r) for r in batch_result.results))
            print(f"<<< {batch_result.batch.workdir} ({num_done + 1}", end="")
            print(f" of {len(tasks)}) done in {batch_result.duration:.2f} sec")
            sys.stdout.flush()
            err_fh.write("\n".join(r.error for r in batch_result.results if r.error))
            err_fh.flush()
            if sum(r.status != "OK" for r in all_results) > cfg.max_errors:
                print(f"\nGot more than {cfg.max_errors} errors, aborting...")
                break

    print("------------------------------- Errors ---------------------------------")
    print("\n".join(r.error for r in all_results if r.error))

    print("\n------------------------------- Timings --------------------------------")
    timings = sorted(r.duration for r in all_results)
    print('Plot: name="timings", title="Timing Distribution", ylabel="time [s]"')
    for p in (100, 99, 98, 95, 90, 80):
        y = percentile(timings, p / 100.0)
        print(f'PlotPoint: name="{p}th_percentile", plot="timings", ', end="")
        print(f'label="{p}th %ile", y={y:.2f}, yerr=0.0')

    if cfg.flag_slow:
        print("\n" + "-" * 15 + "--------------- Slow Tests ---------------" + "-" * 15)
        threshold = 2 * percentile(timings, 0.95)
        outliers = [r for r in all_results if r.duration > threshold]
        maybe_slow = [r for r in outliers if r.fullname not in cfg.slow_suppressions]
        num_suppressed = len(outliers) - len(maybe_slow)
        rerun_tasks: List[Task[BatchResult]] = []
        for b in {r.batch for r in maybe_slow}:
            print(f"Re-running {b.name} to avoid false positives.")
            rerun_tasks.append(asyncio.get_event_loop().create_task(run_batch(b, cfg)))
        rerun_times: Dict[str, float] = {}
        for t in await asyncio.gather(*rerun_tasks):
            rerun_times.update({r.fullname: r.duration for r in t.results})
        stats = {r.fullname: [r.duration, rerun_times[r.fullname]] for r in maybe_slow}
        slow_tests = {k: v for k, v in stats.items() if mean(v) > threshold}
        print(f"Duration threshold (2x 95th %ile): {threshold:.2f} sec")
        print(f"Found {len(slow_tests)} slow tests ({num_suppressed} suppressed):")
        for k, v in slow_tests.items():
            print(f"    {k :<80s} ( {mean(v):6.2f} ±{stdev(v):4.2f} sec)")

    print("\n------------------------------- Summary --------------------------------")
    total_duration = time.perf_counter() - start_time
    num_tests = len(all_results)
    failure_modes = ["RUNTIME FAIL", "TIMED OUT", "HUGE OUTPUT"]
    num_failed = sum(r.status in failure_modes for r in all_results)
    num_wrong = sum(r.status == "WRONG RESULT" for r in all_results)
    num_ok = sum(r.status == "OK" for r in all_results)
    status_ok = (num_ok == num_tests) and (not cfg.flag_slow or not slow_tests)
    print(f"Number of FAILED  tests {num_failed}")
    print(f"Number of WRONG   tests {num_wrong}")
    print(f"Number of CORRECT tests {num_ok}")
    print(f"Total number of   tests {num_tests}")
    summary = f"\nSummary: correct: {num_ok} / {num_tests}"
    summary += f"; wrong: {num_wrong}" if num_wrong > 0 else ""
    summary += f"; failed: {num_failed}" if num_failed > 0 else ""
    summary += f"; slow: {len(slow_tests)}" if cfg.flag_slow and slow_tests else ""
    summary += f"; {total_duration/60.0:.0f}min"
    print(summary)
    print("Status: " + ("OK" if status_ok else "FAILED") + "\n")

    print("*************************** Testing ended ******************************")
    sys.exit(0 if status_ok else 1)


# ======================================================================================
class Config:
    def __init__(self, args: argparse.Namespace):
        self.timeout = args.timeout
        self.use_mpi = args.version.startswith("p")
        default_ompthreads = 2 if "smp" in args.version else 1
        self.ompthreads = args.ompthreads if args.ompthreads else default_ompthreads
        self.mpiranks = args.mpiranks if self.use_mpi else 1
        self.num_workers = int(args.maxtasks / self.ompthreads / self.mpiranks) or 1
        self.workers = Semaphore(self.num_workers)
        self.cp2k_root = Path(__file__).resolve().parent.parent
        self.mpiexec = args.mpiexec
        if "{N}" not in self.mpiexec:  # backwards compatibility
            self.mpiexec = f"{self.mpiexec} ".replace(" ", " -n {N} ", 1).strip()
        self.smoketest = args.smoketest
        self.valgrind = args.valgrind
        self.keepalive = args.keepalive
        self.flag_slow = args.flagslow
        self.binary_dir = args.binary_dir.resolve()
        self.version = args.version
        self.debug = args.debug
        self.max_errors = args.maxerrors
        self.restrictdirs = args.restrictdir if args.restrictdir else [".*"]
        self.skipdirs = args.skipdir if args.skipdir else []
        self.skip_unittests = args.skip_unittests
        self.skip_regtests = args.skip_regtests
        datestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        leaf_dir = f"TEST-{datestamp}"
        self.work_base_dir = (args.workbasedir or args.binary_dir).resolve() / leaf_dir
        self.error_summary = self.work_base_dir / "error_summary"

        # Parse suppression files.
        slow_supps_fn = self.cp2k_root / "tests" / "SLOW_TESTS_SUPPRESSIONS"
        self.slow_suppressions = slow_supps_fn.read_text(encoding="utf8").split("\n")
        huge_supps_fn = self.cp2k_root / "tests" / "HUGE_TESTS_SUPPRESSIONS"
        self.huge_suppressions = huge_supps_fn.read_text(encoding="utf8").split("\n")

        def run_with_capture_stdout(cmd: str) -> bytes:
            # capture_output argument not available before Python 3.7
            return subprocess.run(cmd, shell=True, stdout=PIPE, stderr=DEVNULL).stdout

        # Detect number of GPU devices, if not specified by the user
        if args.num_gpus > 0:
            self.num_gpus = args.num_gpus
        else:
            nv_cmd = "nvidia-smi --query-gpu=gpu_name --format=csv,noheader | wc -l"
            nv_gpus = int(run_with_capture_stdout(nv_cmd))
            amd_cmd = "rocm-smi --showid --csv | grep card | wc -l"
            amd_gpus = int(run_with_capture_stdout(amd_cmd))
            self.num_gpus = nv_gpus + amd_gpus
        self.next_gpu = 0  # Used to assign devices round robin to processes.

    def launch_exe(
        self, exe_stem: str, *args: str, cwd: Optional[Path] = None
    ) -> Coroutine[Any, Any, Process]:
        env = os.environ.copy()
        if self.num_gpus > self.mpiranks:
            visible_gpu_devices: List[str] = []
            for _ in range(self.mpiranks):  # Utilize all available GPU devices.
                self.next_gpu = (self.next_gpu + 1) % self.num_gpus
                visible_gpu_devices.append(str(self.next_gpu))
            env["CUDA_VISIBLE_DEVICES"] = ",".join(visible_gpu_devices)
            env["HIP_VISIBLE_DEVICES"] = ",".join(visible_gpu_devices)
        env["OMP_NUM_THREADS"] = str(self.ompthreads)
        env["PIKA_COMMANDLINE_OPTIONS"] = (
            f"--pika:bind=none --pika:threads={self.ompthreads}"
        )
        exe_name = f"{exe_stem}.{self.version}"
        cmd = [str(self.binary_dir / exe_name)]
        if self.valgrind:
            cmd = ["valgrind", "--error-exitcode=42", "--exit-on-first-error=yes"] + cmd
        if self.use_mpi:
            cmd = self.mpiexec.format(N=self.mpiranks).split() + cmd
        if self.debug:
            print(f"Creating subprocess: {cmd} {args}")
        return asyncio.create_subprocess_exec(
            *cmd, *args, env=env, cwd=cwd, stdin=PIPE, stdout=PIPE, stderr=STDOUT
        )


# ======================================================================================
class TestType:
    """Search pattern for result values, ie. a line in the TEST_TYPES file."""

    def __init__(self, line: str):
        parts = line.rsplit("!", 1)
        self.pattern = parts[0]
        self.pattern = self.pattern.replace("[[:space:]]", r"\s")
        for c in r"|()+":
            self.pattern = self.pattern.replace(c, f"\\{c}")  # escape special chars
        try:
            self.regex = re.compile(self.pattern)
        except Exception:
            print("Bad regex: " + self.pattern)
            raise
        self.column = int(parts[1])

    def grep(self, output: str) -> Optional[str]:
        for line in reversed(output.split("\n")):
            match = self.regex.search(line)
            if match:
                # print("pattern:" + self.pattern)
                # print("line:"+line)
                return line.split()[self.column - 1]
        return None


# ======================================================================================
class Unittest:
    """A unit test, ie. a standalone binary that matches '*_unittest.{cfg.version}'."""

    def __init__(self, name: str, workdir: Path):
        self.name = name
        self.out_path = workdir / (self.name + ".out")


# ======================================================================================
class Regtest:
    """A single input file to test, ie. a line in a TEST_FILES file."""

    def __init__(self, line: str, test_types: List[Optional[TestType]], workdir: Path):
        parts = line.split()
        self.inp_fn = parts[0]
        self.name = self.inp_fn
        self.out_path = workdir / (self.name + ".out")

        self.test_type = test_types[int(parts[1])]
        if self.test_type:
            self.tolerance = float(parts[2])
            self.ref_value_txt = parts[3]
            self.ref_value = float(self.ref_value_txt)


# ======================================================================================
class Batch:
    """A directory of tests, ie. a line in the TEST_DIRS file."""

    def __init__(self, line: str, cfg: Config):
        parts = line.split()
        self.name = parts[0]
        self.requirements = parts[1:]
        self.unittests: List[Unittest] = []
        self.regtests: List[Regtest] = []
        self.src_dir = cfg.cp2k_root / "tests" / self.name
        self.workdir = cfg.work_base_dir / self.name
        self.huge_suppressions = cfg.huge_suppressions

    def requirements_satisfied(self, flags: List[str], mpiranks: int) -> bool:
        result = True
        for r in self.requirements:
            if "mpiranks" in r:
                result &= eval(r.replace("||", " or "))
            elif r.startswith("!"):
                result &= r[1:] not in flags
            else:
                result &= r in flags
        return result


# ======================================================================================
class TestResult:
    def __init__(
        self,
        batch: Batch,
        test: Union[Regtest, Unittest],
        duration: float,
        status: TestStatus,
        value: Optional[float] = None,
        error: Optional[str] = None,
    ):
        self.batch = batch
        self.test = test
        self.duration = duration
        self.status = status
        self.value = value
        self.error = error
        self.fullname = f"{batch.name}/{test.name}"

    def __str__(self) -> str:
        value = f"{self.value:.10g}" if self.value else "-"
        return f"    {self.test.name :<80s} {value :>17} {self.status :>12s} ( {self.duration:6.2f} sec)"


# ======================================================================================
class BatchResult:
    def __init__(self, batch: Batch, results: List[TestResult]):
        self.batch = batch
        self.results = results
        self.duration = sum(float(r.duration) for r in results)


# ======================================================================================
class Cp2kShell:
    def __init__(self, cfg: Config, workdir: Path):
        self.cfg = cfg
        self.workdir = workdir
        self._child: Optional[Process] = None

    async def stop(self) -> None:
        assert self._child
        try:
            self._child.terminate()  # Give mpiexec a chance to shutdown
        except ProcessLookupError:
            pass
        await self._child.communicate()  # Read output to prevent a zombie process.
        self._child = None

    async def start(self) -> None:
        assert self._child is None
        self._child = await self.cfg.launch_exe("cp2k", "--shell", cwd=self.workdir)
        await self.ready()
        await self.sendline("HARSH")  # With harsh mode any error leads to an abort.
        await self.ready()

    async def sendline(self, line: str) -> None:
        if self.cfg.debug:
            print("Sending: " + line)
        assert len(line) < 1024  # input buffer size
        assert self._child and self._child.stdin
        data = (line + "\n").encode("utf8")
        self._child.stdin.write(data)
        try:
            await asyncio.wait_for(self._child.stdin.drain(), timeout=self.cfg.timeout)
        except ConnectionResetError:
            pass

    async def readline(self) -> str:
        assert self._child and self._child.stdout
        data_future = self._child.stdout.readline()
        data = await asyncio.wait_for(data_future, timeout=self.cfg.timeout)
        line = data.decode("utf8", errors="replace")
        if self.cfg.debug:
            print("Received: " + line, end="")
        return line

    async def ready(self, output: Optional[TextIO] = None) -> None:
        while True:
            line = await self.readline()
            if line == "":
                assert self._child
                await self._child.wait()  # child crashed - get return code
                return
            elif line == "* READY\n":
                return  # task finished nominally
            elif output:
                output.write(line)  # task continues

    def returncode(self) -> int:
        assert self._child
        return self._child.returncode if self._child.returncode else 0


# ======================================================================================
async def wait_for_child_process(
    child: Process, timeout: int
) -> Tuple[bytes, int, bool]:
    try:
        output, _ = await asyncio.wait_for(child.communicate(), timeout=timeout)
        timed_out = False
        returncode = child.returncode if child.returncode else 0
    except asyncio.TimeoutError:
        timed_out = True
        returncode = -9
        try:
            child.terminate()  # Give mpiexec a chance to shutdown
        except ProcessLookupError:
            pass
        output, _ = await child.communicate()

    return output, returncode, timed_out


# ======================================================================================
async def run_batch(batch: Batch, cfg: Config) -> BatchResult:
    async with cfg.workers:
        results = []
        if not cfg.skip_unittests:
            results += await run_unittests(batch, cfg)
        if not cfg.skip_regtests:
            results += await run_regtests(batch, cfg)
        return BatchResult(batch, results)


# ======================================================================================
async def run_unittests(batch: Batch, cfg: Config) -> List[TestResult]:
    results: List[TestResult] = []
    for test in batch.unittests:
        start_time = time.perf_counter()
        child = await cfg.launch_exe(test.name, str(cfg.cp2k_root), cwd=batch.workdir)
        output, returncode, timed_out = await wait_for_child_process(child, cfg.timeout)
        duration = time.perf_counter() - start_time
        test.out_path.write_bytes(output)
        output_lines = output.decode("utf8", errors="replace").split("\n")
        output_tail = "\n".join(output_lines[-100:])
        error = "x" * 100 + f"\n{test.out_path}\n{output_tail}\n\n"
        if timed_out:
            error += f"Timed out after {duration} seconds."
            results.append(TestResult(batch, test, duration, "TIMED OUT", error=error))
        elif returncode != 0:
            error += f"Runtime failure with code {returncode}."
            results.append(
                TestResult(batch, test, duration, "RUNTIME FAIL", error=error)
            )
        else:
            results.append(TestResult(batch, test, duration, "OK"))

    return results


# ======================================================================================
async def run_regtests(batch: Batch, cfg: Config) -> List[TestResult]:
    if cfg.keepalive and batch.name not in KEEPALIVE_SKIP_DIRS:
        return await run_regtests_keepalive(batch, cfg)
    else:
        return await run_regtests_classic(batch, cfg)


# ======================================================================================
async def run_regtests_keepalive(batch: Batch, cfg: Config) -> List[TestResult]:
    shell = Cp2kShell(cfg, batch.workdir)
    await shell.start()

    results: List[TestResult] = []
    for test in batch.regtests:
        start_time = time.perf_counter()
        start_dirsize = dirsize(batch.workdir)
        await shell.sendline(f"RUN {test.inp_fn} __STD_OUT__")
        with open(test.out_path, "wt", encoding="utf8", errors="replace") as fh:
            try:
                await asyncio.wait_for(shell.ready(fh), timeout=cfg.timeout)
                timed_out = False
                returncode = shell.returncode()
            except asyncio.TimeoutError:
                timed_out = True
                returncode = -9

        if returncode != 0:
            await shell.stop()
            await shell.start()
        duration = time.perf_counter() - start_time
        output_size = dirsize(batch.workdir) - start_dirsize
        res = eval_regtest(batch, test, duration, output_size, returncode, timed_out)
        results.append(res)

    await shell.stop()
    return results


# ======================================================================================
async def run_regtests_classic(batch: Batch, cfg: Config) -> List[TestResult]:
    results: List[TestResult] = []
    for test in batch.regtests:
        start_time = time.perf_counter()
        start_dirsize = dirsize(batch.workdir)
        child = await cfg.launch_exe("cp2k", test.inp_fn, cwd=batch.workdir)
        output, returncode, timed_out = await wait_for_child_process(child, cfg.timeout)
        test.out_path.write_bytes(output)
        duration = time.perf_counter() - start_time
        output_size = dirsize(batch.workdir) - start_dirsize
        res = eval_regtest(batch, test, duration, output_size, returncode, timed_out)
        results.append(res)

    return results


# ======================================================================================
def dirsize(folder: Path) -> int:
    return sum(f.stat().st_size for f in folder.rglob("*"))


# ======================================================================================
def eval_regtest(
    batch: Batch,
    test: Regtest,
    duration: float,
    output_size: int,
    returncode: int,
    timed_out: bool,
) -> TestResult:
    is_huge_suppressed = f"{batch.name}/{test.name}" in batch.huge_suppressions
    output_bytes = test.out_path.read_bytes() if test.out_path.exists() else b""
    output = output_bytes.decode("utf8", errors="replace")
    output_tail = "\n".join(output.split("\n")[-100:])
    error = "x" * 100 + f"\n{test.out_path}\n"

    # check for timeout
    if timed_out:
        error += f"{output_tail}\n\nTimed out after {duration} seconds."
        return TestResult(batch, test, duration, "TIMED OUT", error=error)

    # check for crash
    if returncode != 0:
        error += f"{output_tail}\n\nRuntime failure with code {returncode}."
        return TestResult(batch, test, duration, "RUNTIME FAIL", error=error)

    # check for huge output
    if output_size > 2 * 1024 * 1024 and not is_huge_suppressed:  # 2MiB limit
        error += f"Test produced {output_size/1024/1024:.2f} MiB of output."
        return TestResult(batch, test, duration, "HUGE OUTPUT", error=error)

    # happy end for test type zero
    if not test.test_type:
        return TestResult(batch, test, duration, "OK")

    # grep result
    value_txt = test.test_type.grep(output)
    if not value_txt:
        error += f"{output_tail}\n\nResult not found: '{test.test_type.pattern}'."
        return TestResult(batch, test, duration, "WRONG RESULT", error=error)

    # parse result
    try:
        value = float(value_txt)
    except:
        error += f"Could not parse result as float: '{value_txt}'."
        return TestResult(batch, test, duration, "WRONG RESULT", error=error)

    # compare result to reference
    diff = value - test.ref_value
    rel_error = abs(diff / test.ref_value if test.ref_value != 0.0 else diff)
    if rel_error > test.tolerance:
        error += f"Difference too large: {rel_error:.2e} > {test.tolerance}, "
        error += f"ref_value: {test.ref_value_txt}, value: {value_txt}."
        return TestResult(batch, test, duration, "WRONG RESULT", value, error)

    # happy end
    return TestResult(batch, test, duration, "OK", value)


# ======================================================================================
def percentile(values: List[float], percent: float) -> float:
    k = (len(values) - 1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return values[int(k)]
    d0 = values[int(f)] * (c - k)
    d1 = values[int(c)] * (k - f)
    return d0 + d1


# ======================================================================================
if __name__ == "__main__":
    # asyncio.run(main()) not available before Python 3.7
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    loop.run_until_complete(main())

# EOF
