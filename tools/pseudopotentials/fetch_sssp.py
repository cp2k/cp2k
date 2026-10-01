#!/usr/bin/env python3
"""Download checksum-pinned SSSP v2.0 PBE/PBEsol pseudopotential libraries.

Usage: python3 tools/pseudopotentials/fetch_sssp.py DIRECTORY [--library pbe-eff]
Without --library, fetch all four Efficiency/Precision collections. Each archive
is verified before extraction into a separate library directory. Existing cached
archives are reused only when their SHA-256 matches the pinned release manifest.
This utility downloads datasets; calculation support depends on the CP2K method.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import tarfile
import urllib.request


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
