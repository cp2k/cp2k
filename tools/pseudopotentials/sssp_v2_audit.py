#!/usr/bin/env python3
"""Fetch pinned SSSP v2.0 datasets and optionally audit NC import/transformation.

The optional checker is CP2K's atom_upf_unittest executable. Passing this
audit establishes import and finite transformation output, not physical
accuracy, basis convergence, or support for USPP/PAW calculations.
"""

import argparse
from collections import Counter
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile
import urllib.request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--check-executable", type=Path)
    args = parser.parse_args()
    root = args.directory.resolve()
    root.mkdir(parents=True, exist_ok=True)
    config = json.loads(Path(__file__).with_name("sssp_v2_archives.json").read_text())
    manifest = {"sssp_version": config["version"], "archives": []}
    failed = []
    for source in config["archives"]:
        archive = root / source["archive"]
        if not archive.exists():
            with urllib.request.urlopen(source["url"], timeout=120) as response:
                data = response.read()
            if sha256(data).hexdigest() != source["sha256"]:
                raise RuntimeError(f"Checksum mismatch: {source['archive']}")
            archive.write_bytes(data)
        if sha256(archive.read_bytes()).hexdigest() != source["sha256"]:
            raise RuntimeError(f"Checksum mismatch: {archive}")
        dest = root / archive.name.removeprefix("SSSP-lib-").removesuffix("-v2.tar.gz")
        dest.mkdir(exist_ok=True)
        with tarfile.open(archive) as tf:
            tf.extractall(dest, filter="data")
        potentials = []
        for path in sorted(dest.rglob("*.upf")):
            data = path.read_bytes()
            match = re.search(rb"<PP_HEADER\b.*?/>", data, re.DOTALL)
            if match is None:
                raise RuntimeError(f"Missing UPF header: {path}")
            attrs = dict(re.findall(r'(\w+)\s*=\s*"([^"]*)"', match[0].decode()))
            family = attrs["pseudo_type"].strip()
            if family == "USPP":
                family = "US"
            record = {
                "path": str(path.relative_to(root)),
                "sha256": sha256(data).hexdigest(),
                "element": attrs["element"].strip(),
                "family": family,
                "z_valence": float(attrs["z_valence"]),
            }
            if args.check_executable and family == "NC":
                log = path.with_suffix(".check.out")
                env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
                with log.open("w") as stream:
                    try:
                        result = subprocess.run(
                            [str(args.check_executable.resolve()), str(path)],
                            stdout=stream,
                            stderr=subprocess.STDOUT,
                            env=env,
                            timeout=180,
                        )
                        record["check_returncode"] = result.returncode
                    except subprocess.TimeoutExpired:
                        record["check_returncode"] = "timeout"
                if record["check_returncode"] != 0:
                    failed.append(record["path"])
            potentials.append(record)
        counts = dict(Counter(p["family"] for p in potentials))
        expected = {"NC": 20, "US": 40, "PAW": 35}
        if "-prec-" in archive.name:
            expected = {"NC": 21, "US": 37, "PAW": 37}
        if counts != expected:
            raise RuntimeError(f"Unexpected SSSP inventory: {archive.name}: {counts}")
        manifest["archives"].append(
            {**source, "counts": counts, "potentials": potentials}
        )
        print(archive.name, counts, flush=True)
    manifest["failed_checks"] = failed
    (root / "audit.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if failed:
        raise SystemExit(
            f"{len(failed)} NC import/transformation checks failed; see audit.json"
        )


if __name__ == "__main__":
    main()
