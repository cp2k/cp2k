#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Run GauXC/SKALA CP2K inputs on different MPI ranks and write a compact TSV.

The script is intentionally lightweight and tailored for quick force/energy checks
on shared filesystems:

- runs one or more input files for the requested MPI ranks,
- propagates GauXC/SKALA environment toggles,
- records wall-clock runtime and final total FORCE_EVAL energy,
- polls selected GPU metrics while the job is running,
- writes a mergeable `summary.tsv` that is ready for comparison.
"""

from __future__ import annotations

import argparse
import csv
import os
import pathlib
import re
import shlex
import subprocess
import threading
import time
from dataclasses import dataclass
from typing import Dict, List, Optional

ENERGY_VALUE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?"
LOG = re.compile(
    rf"^\s*ENERGY\|\s*Total FORCE_EVAL.*?(?P<energy>{ENERGY_VALUE})\s*$",
    re.MULTILINE,
)


@dataclass
class RunResult:
    label: str
    ranks: int
    rc: int
    wall_s: float
    cp2k_s: float
    energy_ha: Optional[float]
    gpu_peak_mem: Dict[int, int]
    gpu_peak_util: Dict[int, int]


def _to_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "E"))


def parse_last_energy(path: pathlib.Path) -> Optional[float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = None
    for m in LOG.finditer(text):
        match = m
    if match is None:
        return None
    return _to_float(match.group("energy"))


def query_gpu_snapshot(gpu_indices: Optional[List[int]] = None):
    cmd = [
        "nvidia-smi",
        "--query-gpu=index,memory.used,utilization.gpu",
        "--format=csv,noheader,nounits",
    ]
    output = subprocess.check_output(cmd, text=True)
    peaks_mem = {}
    peaks_util = {}

    requested = None if gpu_indices is None else set(gpu_indices)
    for raw_line in output.splitlines():
        parts = [x.strip() for x in raw_line.split(",")]
        if len(parts) < 3:
            continue
        try:
            gpu = int(parts[0])
            mem = int(parts[1])
            util = int(parts[2])
        except ValueError:
            continue
        if requested is not None and gpu not in requested:
            continue
        peaks_mem[gpu] = mem
        peaks_util[gpu] = util
    return peaks_mem, peaks_util


def sample_gpu_metrics(
    gpu_indices: List[int],
    interval_s: float,
    stop: threading.Event,
    peak_mem: Dict[int, int],
    peak_util: Dict[int, int],
):
    while not stop.is_set():
        try:
            mem_now, util_now = query_gpu_snapshot(gpu_indices)
            for gpu in gpu_indices:
                if gpu in mem_now:
                    peak_mem[gpu] = max(peak_mem.get(gpu, 0), mem_now[gpu])
                if gpu in util_now:
                    peak_util[gpu] = max(peak_util.get(gpu, 0), util_now[gpu])
        except (subprocess.CalledProcessError, FileNotFoundError, ValueError):
            # If nvidia-smi is unavailable or busy, skip sampling for this interval.
            pass
        time.sleep(max(interval_s, 0.2))


def launch_run(
    cp2k: pathlib.Path,
    launcher: str,
    input_file: pathlib.Path,
    output_file: pathlib.Path,
    ranks: int,
    launcher_args: List[str],
    env: Dict[str, str],
) -> int:
    cmd = []
    if ranks > 1:
        cmd.extend([launcher, "-np", str(ranks)])
        cmd.extend(launcher_args)
    cmd.extend([str(cp2k), "-i", str(input_file), "-o", str(output_file)])

    stdout_file = output_file.with_suffix(output_file.suffix + ".stdout")
    with stdout_file.open("w", encoding="utf-8") as fhandle:
        proc = subprocess.run(cmd, env=env, stdout=fhandle, stderr=subprocess.STDOUT)
    return proc.returncode


def run_case(
    label: str,
    cp2k: pathlib.Path,
    launcher: str,
    input_file: pathlib.Path,
    ranks: int,
    outdir: pathlib.Path,
    launcher_args: List[str],
    env: Dict[str, str],
    gpu_indices: List[int],
    gpu_interval: float,
) -> RunResult:
    outdir.mkdir(parents=True, exist_ok=True)
    log_file = outdir / f"{label}_r{ranks}.out"

    peak_mem: Dict[int, int] = {gpu: 0 for gpu in gpu_indices}
    peak_util: Dict[int, int] = {gpu: 0 for gpu in gpu_indices}
    stop = threading.Event()

    # Probe once up front so we can include GPU ids even before launch.
    try:
        query_gpu_snapshot(gpu_indices)
        monitor = threading.Thread(
            target=sample_gpu_metrics,
            args=(gpu_indices, gpu_interval, stop, peak_mem, peak_util),
            daemon=True,
        )
        monitor.start()
    except (FileNotFoundError, subprocess.CalledProcessError):
        monitor = None

    t0 = time.perf_counter()
    rc = launch_run(cp2k, launcher, input_file, log_file, ranks, launcher_args, env)
    wall_s = time.perf_counter() - t0

    if monitor is not None:
        stop.set()
        monitor.join(timeout=2.0)

    energy = parse_last_energy(log_file)

    return RunResult(
        label=label,
        ranks=ranks,
        rc=rc,
        wall_s=wall_s,
        cp2k_s=wall_s,
        energy_ha=energy,
        gpu_peak_mem=peak_mem,
        gpu_peak_util=peak_util,
    )


def write_summary(path: pathlib.Path, results: List[RunResult], gpu_indices: List[int]):
    fieldnames = [
        "label",
        "ranks",
        "rc",
        "wall_s",
        "cp2k_s",
        "energy_ha",
    ]
    for gid in gpu_indices:
        fieldnames.extend([f"gpu{gid}_mem_mb", f"gpu{gid}_util_pct"])

    with path.open("w", newline="", encoding="utf-8") as fhandle:
        writer = csv.DictWriter(fhandle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in results:
            row = {
                "label": r.label,
                "ranks": r.ranks,
                "rc": r.rc,
                "wall_s": f"{r.wall_s:.6f}",
                "cp2k_s": f"{r.cp2k_s:.6f}",
                "energy_ha": f"{r.energy_ha:.14f}" if r.energy_ha is not None else "NA",
            }
            for gid in gpu_indices:
                row[f"gpu{gid}_mem_mb"] = r.gpu_peak_mem.get(gid, "NA")
                row[f"gpu{gid}_util_pct"] = r.gpu_peak_util.get(gid, "NA")
            writer.writerow(row)


def parse_comma_ints(value: str) -> List[int]:
    if not value.strip():
        return []
    return [int(part) for part in value.split(",") if part.strip()]


def write_rank_gpu_wrapper(
    cp2k: pathlib.Path, outdir: pathlib.Path, gpu_indices: List[int]
) -> pathlib.Path:
    if not gpu_indices:
        raise ValueError("rank GPU binding requires at least one GPU index")

    wrapper = outdir / "cp2k_rank_gpu.sh"
    devices = " ".join(shlex.quote(str(gpu)) for gpu in gpu_indices)
    script = f"""#!/usr/bin/env bash
set -euo pipefail

rank="${{OMPI_COMM_WORLD_LOCAL_RANK:-}}"
if [[ -z "${{rank}}" ]]; then
  rank="${{PMI_LOCAL_RANK:-}}"
fi
if [[ -z "${{rank}}" ]]; then
  rank="${{SLURM_LOCALID:-0}}"
fi

devices=({devices})
idx=$((rank % {len(gpu_indices)}))
export CUDA_VISIBLE_DEVICES="${{devices[$idx]}}"
exec {shlex.quote(str(cp2k))} "$@"
"""
    wrapper.write_text(script, encoding="utf-8")
    wrapper.chmod(0o755)
    return wrapper


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--cp2k", type=pathlib.Path, required=True, help="cp2k.psmp executable"
    )
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        type=pathlib.Path,
        help="Input file to run. Use multiple times for multiple runs.",
    )
    parser.add_argument(
        "--label",
        action="append",
        default=[],
        help="Per-input label. Defaults to input basename.",
    )
    parser.add_argument("--ranks", nargs="+", default=["1", "2"], help="MPI rank list")
    parser.add_argument("--launcher", default="mpirun", help="MPI launcher")
    parser.add_argument(
        "--outdir", default="skala-gauxc-multigpu-bench", help="Output directory"
    )
    parser.add_argument(
        "--summary", default="summary.tsv", help="Summary file name in outdir"
    )
    parser.add_argument(
        "--gpu-devices", default="0,1", help="CUDA devices visible to the run"
    )
    parser.add_argument("--gpu-sample-interval", type=float, default=1.0)
    parser.add_argument("--fill-fraction", type=float, default=0.95)
    parser.add_argument("--atom-chunk-size", type=int, default=3)
    parser.add_argument("--device-memory-cap", default="12000000000")
    parser.add_argument("--distributed-torch", choices=("0", "1"), default="0")
    parser.add_argument("--gradient-mpi-runtime", choices=("0", "1"), default="0")
    parser.add_argument(
        "--rank-gpu-bind",
        action="store_true",
        help="Bind local MPI ranks round-robin to the requested CUDA devices.",
    )
    parser.add_argument(
        "--launcher-arg", action="append", default=[], help="Extra launcher args"
    )
    args = parser.parse_args()

    cp2k = args.cp2k.expanduser().resolve()
    input_files = [f.expanduser().resolve() for f in args.input]
    labels = args.label or [p.stem for p in input_files]
    if len(labels) != len(input_files):
        raise SystemExit("--label count must match --input count")

    ranks = [int(r) for r in args.ranks]
    gpu_indices = parse_comma_ints(args.gpu_devices)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    base_env = os.environ.copy()
    base_env["OMP_NUM_THREADS"] = "1"
    base_env["GAUXC_DEVICE_RUNTIME_FILL_FRACTION"] = f"{args.fill_fraction}"
    base_env["GAUXC_ONEDFT_ATOM_CHUNK_SIZE"] = str(args.atom_chunk_size)
    base_env["GAUXC_DEVICE_MEMORY_CAP"] = str(args.device_memory_cap)
    base_env["GAUXC_ONEDFT_DISTRIBUTED_TORCH"] = args.distributed_torch
    base_env["GAUXC_GRADIENT_MPI_RUNTIME"] = args.gradient_mpi_runtime
    if gpu_indices:
        base_env["CUDA_VISIBLE_DEVICES"] = ",".join(str(g) for g in gpu_indices)

    cp2k_for_run = cp2k
    if args.rank_gpu_bind:
        cp2k_for_run = write_rank_gpu_wrapper(cp2k, outdir, gpu_indices)

    results: List[RunResult] = []
    for inp, lbl in zip(input_files, labels):
        for nrank in ranks:
            print(f"Running {inp.name} with {nrank} rank(s)...")
            res = run_case(
                lbl,
                cp2k_for_run,
                args.launcher,
                inp,
                nrank,
                outdir,
                args.launcher_arg,
                base_env,
                gpu_indices,
                args.gpu_sample_interval,
            )
            results.append(res)
            if res.energy_ha is not None:
                print(
                    f"  rc={res.rc} wall_s={res.wall_s:.3f} energy={res.energy_ha:.12f}"
                )
            else:
                print(f"  rc={res.rc} wall_s={res.wall_s:.3f} energy=NA")

    summary = outdir / args.summary
    write_summary(summary, results, gpu_indices)
    print(f"Summary written to {summary}")
    return 0 if all(r.rc == 0 for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
