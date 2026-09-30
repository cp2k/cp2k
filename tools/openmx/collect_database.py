#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Collect the PAOs linked by the official OpenMX 2019 database, with provenance.

Only explicit links below the database URL are followed. Missing elements are
not invented from the periodic table or from unpublished archive contents.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from html import unescape
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import urljoin, urlsplit, urlunsplit
from urllib.request import urlopen

from convert_basis import ELEMENTS

ROOT = "https://www.openmx-square.org/vps_pao2019/"
GUIDE = "https://www.openmx-square.org/openmx_man3.9/node27.html"
STEM = r"([A-Z][a-z]?)(?:[0-9]+(?:\.[0-9]+)?|\*\.\*)([A-Za-z_0-9]*)"
SELECTION = re.compile(STEM + r"-((?:[spdfg][1-9][0-9]*)+)")


class Links(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.links = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            self.links.extend(v for k, v in attrs if k.lower() == "href" and v)


def database_link(page, link):
    parsed = urlsplit(urljoin(page, link))
    if parsed.hostname not in ("www.openmx-square.org", "openmx-square.org"):
        return None
    normalized = urlunsplit(("https", "www.openmx-square.org", parsed.path, "", ""))
    if not normalized.startswith(ROOT):
        return None
    if normalized.endswith("/"):
        normalized += "index.html"
    return normalized


def family(stem):
    match = re.fullmatch(STEM, stem)
    if not match:
        raise ValueError(f"Unrecognized PAO stem: {stem}")
    return match[1], match[2]


def choices(text, source):
    plain = unescape(re.sub(r"<[^>]*>", "", text))
    found = {}
    for match in SELECTION.finditer(plain):
        key = (match[1], match[2], match[3])
        found[key] = {
            "element": match[1],
            "family_suffix": match[2],
            "shells": match[3],
            "source_page": source,
            "wildcard_cutoff": "*" in match[0],
        }
    return list(found.values())


def cached_fetch(url, path):
    if not path.exists():
        with urlopen(url, timeout=60) as response:
            if response.geturl().split(":", 1)[0] != "https":
                raise ValueError(f"Refusing non-HTTPS redirect for {url}")
            data = response.read()
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as output:
            output.write(data)
    return path.read_bytes()


def collect(cache, workers=4):
    pages, pending, source_pages = {}, {ROOT + "index.html"}, {}
    while pending:
        batch = sorted(pending - pages.keys())
        if not batch:
            break
        pending = set()
        with ThreadPoolExecutor(max_workers=workers) as pool:
            fetched = pool.map(
                lambda url: cached_fetch(url, cache / "html" / url.removeprefix(ROOT)),
                batch,
            )
            for url, raw in zip(batch, fetched):
                html = raw.decode("utf-8", errors="replace")
                pages[url] = {"sha256": hashlib.sha256(raw).hexdigest(), "html": html}
                for link in Links(html).links:
                    target = database_link(url, link)
                    if target is None:
                        continue
                    if target.endswith(".pao"):
                        source_pages.setdefault(target, set()).add(url)
                    elif target.endswith("/index.html") and target not in pages:
                        pending.add(target)
        print(
            f"Discovered {len(pages)} pages, {len(source_pages)} PAO links", flush=True
        )
    if not source_pages:
        raise ValueError("Database contains no PAO links; refusing an empty catalog")
    guide = cached_fetch(GUIDE, cache / "html" / "basis-guideline.html")
    guide_choices = choices(guide.decode("utf-8", errors="replace"), GUIDE)
    sources = []

    def download(url):
        stem = Path(urlsplit(url).path).stem
        element, suffix = family(stem)
        if element not in ELEMENTS[1:]:
            raise ValueError(f"PAO link is not an element: {url}")
        raw = cached_fetch(url, cache / "pao" / (stem + ".pao"))
        available = list(guide_choices)
        for page in sorted(source_pages[url]):
            for choice in choices(pages[page]["html"], page):
                # A leaf page's Nd*.*-s2p1d1 refers to its own _OC PAOs even
                # when the prose omits the suffix. Do not borrow a concrete
                # named basis from another valence family.
                if (
                    choice["wildcard_cutoff"]
                    and choice["element"] == element
                    and not choice["family_suffix"]
                ):
                    choice["family_suffix"] = suffix
                available.append(choice)
        selected = {}
        for choice in available:
            if (choice["element"], choice["family_suffix"]) == (element, suffix):
                selected.setdefault(choice["shells"], set()).add(choice["source_page"])
        # Only used where neither the official table nor the element page gives
        # any explicit shell count. It is marked as a policy choice, never as an
        # OpenMX recommendation. The full PAO remains available for other choices.
        fallback = not selected
        if fallback:
            selected["s3p3d2f1"] = set()
        return {
            "element": element,
            "stem": stem,
            "family_suffix": suffix,
            "url": url,
            "source_pages": sorted(source_pages[url]),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "bytes": len(raw),
            "fallback_shell_policy": fallback,
            "selections": [
                {
                    "label": stem + "-" + shells,
                    "shells": shells,
                    "pattern_source_pages": sorted(origins),
                }
                for shells, origins in sorted(selected.items())
            ],
        }

    urls = sorted(source_pages)
    stems = [Path(urlsplit(url).path).stem for url in urls]
    if len(stems) != len(set(stems)):
        raise ValueError("Distinct source URLs have colliding PAO filenames")
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for source in pool.map(download, urls):
            sources.append(source)
            print(
                f"Cached {source['stem']}: {len(source['selections'])} shell patterns",
                flush=True,
            )
    sources.sort(key=lambda s: (ELEMENTS.index(s["element"]), s["stem"]))
    return {
        "schema": "openmx-2019-source-catalog-v1",
        "database_url": ROOT,
        "scope": "All PAO links reachable through official database index pages; not every combinatorial shell selection",
        "selection_policy": "All published shell-count patterns per element/family, expanded over every linked cutoff; explicit fallback where absent",
        "guide_url": GUIDE,
        "guide_sha256": hashlib.sha256(guide).hexdigest(),
        "pages": [
            {"url": url, "sha256": page["sha256"]}
            for url, page in sorted(pages.items())
        ],
        "elements": sorted({s["element"] for s in sources}, key=ELEMENTS.index),
        "sources": sources,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--overwrite", action="store_true", help="Replace an existing generated catalog"
    )
    parser.add_argument(
        "--fetch-only",
        action="store_true",
        help="Download and verify an existing pinned manifest; do not discover a new catalog",
    )
    args = parser.parse_args()
    if args.manifest.exists() and not (args.fetch_only or args.overwrite):
        parser.error("Manifest already exists; choose a new output")
    if not 1 <= args.workers <= 8:
        parser.error("Use 1..8 download workers")
    if args.fetch_only:
        manifest = json.loads(args.manifest.read_text())
        if manifest.get("schema") != "openmx-2019-source-catalog-v1":
            parser.error("Invalid manifest schema")

        def fetch(source):
            url, stem = source["url"], source["stem"]
            if (
                database_link(ROOT, url) != url
                or not url.startswith("https://")
                or Path(urlsplit(url).path).stem != stem
            ):
                raise ValueError("Invalid PAO URL/stem in manifest")
            family(stem)
            raw = cached_fetch(url, args.cache / "pao" / (stem + ".pao"))
            if hashlib.sha256(raw).hexdigest() != source["sha256"]:
                raise ValueError(
                    f"Checksum mismatch: {stem}; do not use this file for the pinned basis"
                )
            return stem

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for stem in pool.map(fetch, manifest["sources"]):
                print(f"Verified {stem}", flush=True)
        return
    manifest = collect(args.cache, args.workers)
    with args.manifest.open("w" if args.overwrite else "x") as output:
        json.dump(manifest, output, indent=2)
        output.write("\n")
    print(
        f"Catalog: {len(manifest['elements'])} elements, {len(manifest['sources'])} PAOs, "
        f"{sum(len(s['selections']) for s in manifest['sources'])} basis selections"
    )


if __name__ == "__main__":
    main()
