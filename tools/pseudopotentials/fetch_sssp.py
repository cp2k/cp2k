#!/usr/bin/env python3
"""Download checksum-pinned SSSP v2.0 PBE/PBEsol pseudopotential libraries.

Usage: python3 tools/pseudopotentials/fetch_sssp.py DIRECTORY [--library pbe-eff]
Without --library, fetch all four Efficiency/Precision collections. Each archive
is verified before extraction into a separate library directory. Existing cached
archives are reused only when their SHA-256 matches the pinned release manifest.
With --install, install unchanged UPFs and cutoff files under
DIRECTORY/{pbe,pbesol}/{efficiency,precision}, without retaining new archives.
DIRECTORY is explicitly chosen by the user. Downloading original files does not
prepare missing partial waves; see README.md for the separate preparation step.
This utility downloads datasets; calculation support depends on the CP2K method.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import tarfile
import tempfile
import urllib.request


def install_library(root, name, source):
    functional, table = name.split("-")
    destination = root / functional / {"eff": "efficiency", "prec": "precision"}[table]
    with tempfile.TemporaryDirectory(prefix="sssp-download-", dir=root) as temporary:
        staging = Path(temporary)
        archive = root / source["archive"]
        if not archive.exists():
            archive = staging / source["archive"]
            with urllib.request.urlopen(source["url"], timeout=120) as response:
                with archive.open("wb") as output:
                    shutil.copyfileobj(response, output)
        if sha256(archive.read_bytes()).hexdigest() != source["sha256"]:
            raise RuntimeError(f"Checksum mismatch: {source['archive']}")
        names = set()
        pending = []
        with tarfile.open(archive) as contents:
            for member in contents:
                if member.isdir():
                    continue
                filename = Path(member.name).name
                if (
                    not member.isfile()
                    or filename in names
                    or not (filename.endswith(".upf") or filename == "cutoffs.json")
                ):
                    raise RuntimeError(f"Unexpected archive member: {member.name}")
                names.add(filename)
                with contents.extractfile(member) as input_file:
                    data = input_file.read()
                target = destination / filename
                if target.exists():
                    if target.read_bytes() != data:
                        raise RuntimeError(f"Existing file differs: {target}")
                    continue
                staged = staging / filename
                staged.write_bytes(data)
                pending.append((staged, target))
        destination.mkdir(parents=True, exist_ok=True)
        for staged, target in pending:
            staged.replace(target)
    print(f"{name}: verified and installed to {destination}", flush=True)


def main():
    manifest = json.loads(Path(__file__).with_name("sssp_v2_archives.json").read_text())
    libraries = {
        item["archive"].removeprefix("SSSP-lib-").removesuffix("-v2.tar.gz"): item
        for item in manifest["archives"]
    }
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory", type=Path, help="Archive cache and extracted libraries"
    )
    parser.add_argument(
        "--install",
        action="store_true",
        help="Install original UPFs and cutoffs under DIRECTORY without retaining new archives",
    )
    parser.add_argument(
        "--library",
        action="append",
        choices=sorted(libraries),
        help="Select a library; repeat to fetch more than one (default: all)",
    )
    args = parser.parse_args()
    root = args.directory.resolve()
    root.mkdir(parents=True, exist_ok=True)
    for name in dict.fromkeys(args.library or libraries):
        source = libraries[name]
        if args.install:
            install_library(root, name, source)
            continue
        archive = root / source["archive"]
        if archive.exists():
            data = archive.read_bytes()
        else:
            with urllib.request.urlopen(source["url"], timeout=120) as response:
                data = response.read()
        if sha256(data).hexdigest() != source["sha256"]:
            raise RuntimeError(f"Checksum mismatch: {source['archive']}")
        if not archive.exists():
            archive.write_bytes(data)
        destination = root / name
        destination.mkdir(exist_ok=True)
        with tarfile.open(archive) as contents:
            contents.extractall(destination, filter="data")
        print(f"{name}: verified and extracted to {destination}", flush=True)


if __name__ == "__main__":
    main()
