#!/usr/bin/env python3
"""Add verified generator partial waves to a scalar ultrasoft UPF 2.0.1 file.

Usage: add_upf_full_waves.py ORIGINAL GENERATED OUTPUT
Regenerate the original potential with its original recipe and generator,
requesting projector-indexed AE and PS waves (QE ld1: lsave_wfc=.true.).
The mesh, quadrature, potential, projectors and augmentation must reproduce
the original numeric data exactly. Only atomic guess waves/density allow
roundoff differences. Existing original data are preserved byte for byte,
apart from has_wfc; only PP_FULL_WFC is added. No missing core kinetic data
are inferred. OUTPUT must not already exist.
"""

import argparse
from hashlib import sha256
import math
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET
from xml.parsers import expat


def real(text):
    value = float(text.replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError("Nonfinite UPF numeric value")
    return value


def values(node):
    if list(node):
        return []
    result = [real(word) for word in (node.text or "").split()]
    if "size" in node.attrib and int(node.attrib["size"]) != len(result):
        raise ValueError(f"Incorrect array size in {node.tag}")
    return result


def flag(header, name):
    value = header.get(name, "false").strip().lower().strip(".")
    if value not in ("true", "false", "t", "f", "1", "0"):
        raise ValueError(f"Invalid logical attribute {name}")
    return value in ("true", "t", "1")


def normalized(text):
    try:
        value = float(text.replace("D", "E").replace("d", "e"))
    except ValueError:
        return " ".join(text.split()).lower()
    if not math.isfinite(value):
        raise ValueError("Nonfinite UPF numeric attribute")
    return value


def children(node):
    result = {child.tag: child for child in node}
    if len(result) != len(node):
        raise ValueError(f"Repeated child tag in {node.tag}")
    return result


def compare(original, generated):
    if original.tag != generated.tag:
        raise ValueError("Different UPF elements")
    ignored = {"columns", "type"}
    if original.tag == "PP_HEADER":
        ignored.update(("author", "date", "generated", "comment", "has_wfc"))
    a = {
        key: normalized(value)
        for key, value in original.attrib.items()
        if key not in ignored
    }
    b = {
        key: normalized(value)
        for key, value in generated.attrib.items()
        if key not in ignored
    }
    if a != b:
        raise ValueError(f"Different attributes in {original.tag}")
    ca, cb = children(original), children(generated)
    if original.tag == "UPF":
        for name in ("PP_INFO", "PP_FULL_WFC"):
            ca.pop(name, None)
            cb.pop(name, None)
    if ca.keys() != cb.keys():
        raise ValueError(f"Different data sections in {original.tag}")
    if ca:
        for name in ca:
            compare(ca[name], cb[name])
    elif original.tag != "PP_HEADER":
        va, vb = values(original), values(generated)
        if len(va) != len(vb):
            raise ValueError(f"Different array lengths in {original.tag}")
        tolerance = 0.0
        if original.tag.startswith("PP_CHI.") or original.tag == "PP_RHOATOM":
            tolerance = 64 * sys.float_info.epsilon * max(map(abs, va), default=0.0)
        if any(abs(x - y) > tolerance for x, y in zip(va, vb)):
            raise ValueError(f"Different numeric data in {original.tag}")


def top_level_spans(data):
    """Locate actual XML elements, preserving comments and CDATA in the output."""
    parser = expat.ParserCreate()
    depth = 0
    start = 0
    spans = {}
    root_end = None

    def opening(name, attributes):
        nonlocal depth, start
        if depth == 1:
            start = parser.CurrentByteIndex
        depth += 1

    def closing(name):
        nonlocal depth, root_end
        depth -= 1
        end = parser.CurrentByteIndex
        if depth == 1:
            if data[end : end + 2] == b"</":
                end = data.index(b">", end) + 1
            if name in spans:
                raise ValueError(f"Repeated top-level tag {name}")
            spans[name] = (start, end)
        elif depth == 0:
            root_end = end

    parser.StartElementHandler = opening
    parser.EndElementHandler = closing
    parser.Parse(data, True)
    if root_end is None:
        raise ValueError("Missing UPF closing tag")
    return spans, root_end


def enrich(original_data, generated_data):
    original, generated = (
        ET.fromstring(data) for data in (original_data, generated_data)
    )
    for tree in (original, generated):
        if tree.tag != "UPF" or tree.get("version") != "2.0.1":
            raise ValueError("Requires UPF version 2.0.1")
        header = tree.find("PP_HEADER")
        if (
            header is None
            or not flag(header, "is_ultrasoft")
            or flag(header, "is_paw")
            or flag(header, "has_so")
        ):
            raise ValueError(
                "Requires scalar ultrasoft data, without PAW or spin-orbit flags"
            )
    if (
        flag(original.find("PP_HEADER"), "has_wfc")
        or original.find("PP_FULL_WFC") is not None
    ):
        raise ValueError("Original already declares full partial waves")
    if not flag(generated.find("PP_HEADER"), "has_wfc"):
        raise ValueError("Generator file does not declare full partial waves")
    compare(original, generated)
    header = original.find("PP_HEADER")
    mesh, count = int(header.attrib["mesh_size"]), int(header.attrib["number_of_proj"])
    full = generated.find("PP_FULL_WFC")
    if full is None or mesh < 3 or count < 1:
        raise ValueError("Missing full waves, mesh or projectors")
    for name in ("PP_R", "PP_RAB"):
        radial = original.find(f"PP_MESH/{name}")
        if radial is None or len(values(radial)) != mesh:
            raise ValueError(f"Incomplete radial array {name}")
    if int(full.get("number_of_wfc", str(count))) != count:
        raise ValueError("Full waves must be indexed by projectors")
    wave_nodes = children(full)
    expected = {
        f"PP_{kind}WFC.{i}" for kind in ("AE", "PS") for i in range(1, count + 1)
    }
    if wave_nodes.keys() != expected:
        raise ValueError("Incomplete or unexpected full-wave channels")
    for i in range(1, count + 1):
        beta = original.find(f"PP_NONLOCAL/PP_BETA.{i}")
        if beta is None or not 0 < len(values(beta)) <= mesh:
            raise ValueError(f"Missing or incomplete projector {i}")
        angular = int(beta.attrib["angular_momentum"])
        for kind in ("AE", "PS"):
            wave = wave_nodes[f"PP_{kind}WFC.{i}"]
            if (
                int(wave.get("index", str(i))) != i
                or int(wave.attrib["l"]) != angular
                or angular < 0
            ):
                raise ValueError(f"Different full-wave channel for projector {i}")
            if len(values(wave)) != mesh:
                raise ValueError(f"Incomplete full wave for projector {i}")
    spans, root_end = top_level_spans(original_data)
    generated_spans, _ = top_level_spans(generated_data)
    start, end = spans["PP_HEADER"]
    old_header = original_data[start:end]
    if "has_wfc" in header.attrib:
        new_header, replacements = re.subn(
            rb"\bhas_wfc\s*=\s*(['\"])[^'\"]*\1",
            b'has_wfc="true"',
            old_header,
            count=1,
        )
        if replacements != 1:
            raise ValueError("Cannot update has_wfc attribute")
    else:
        opening = re.match(rb"<PP_HEADER\b(?:[^>\"']|\"[^\"]*\"|'[^']*')*>", old_header)
        if opening is None:
            raise ValueError("Cannot locate PP_HEADER opening tag")
        position = opening.end() - (2 if opening.group().endswith(b"/>") else 1)
        new_header = old_header[:position] + b' has_wfc="true"' + old_header[position:]
    full_start, full_end = generated_spans["PP_FULL_WFC"]
    return (
        original_data[:start]
        + new_header
        + original_data[end:root_end]
        + b"\n  "
        + generated_data[full_start:full_end]
        + b"\n"
        + original_data[root_end:]
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path)
    parser.add_argument("generated", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        original, generated = args.original.read_bytes(), args.generated.read_bytes()
        result = enrich(original, generated)
        with args.output.open("xb") as output:
            output.write(result)
    except (OSError, ValueError, KeyError, ET.ParseError, expat.ExpatError) as error:
        parser.error(str(error))
    for name, data in (
        ("Original", original),
        ("Generator", generated),
        ("Output", result),
    ):
        print(f"{name} SHA-256: {sha256(data).hexdigest()}")
    print(f"Verified original data and added full partial waves: {args.output}")


if __name__ == "__main__":
    main()
