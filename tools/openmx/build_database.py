#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Offline, resumable generation of BASIS_OMX from a pinned source catalog.

All linked source families/cutoffs are processed. Numerically invalid entries
are never exported. Stable but inaccurate fits require --include-inaccurate and
carry per-entry warnings; passing thresholds is not an energy/force validation.
"""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

import convert_basis
from convert_basis import BasisSpec, ELEMENTS, Fit, PaoFile, fit_channel, render_basis

DEFAULT_LIMITS = {
    "radial_l2": 0.01,
    "kinetic_absolute_hartree": 0.001,
    "same_center_overlap": 0.001,
    "normalization": 1e-6,
}


def channel_quality(diagnostics, limits):
    observed = {
        "radial_l2": max(diagnostics["normalized_radial_l2_errors"]),
        "kinetic_absolute_hartree": max(diagnostics["kinetic_absolute_errors_hartree"]),
        "same_center_overlap": diagnostics["same_center_overlap_max_error"],
        "normalization": diagnostics["analytic_vs_quadrature_norm_max_error"],
    }
    if not all(math.isfinite(v) for v in observed.values()):
        return False, math.inf, ["nonfinite diagnostics"]
    failed = [key for key, value in observed.items() if value > limits[key]]
    numeric = observed["normalization"] <= limits["normalization"]
    score = max(observed[k] / limits[k] for k in observed if k != "normalization")
    return numeric, score, failed


def fit_payload(fit):
    return {
        "l": fit.l,
        "exponents": fit.exponents.tolist(),
        "coefficients": fit.coefficients.tolist(),
        "diagnostics": fit.diagnostics,
    }


def restore_fit(payload):
    return Fit(
        payload["l"],
        np.array(payload["exponents"]),
        np.array(payload["coefficients"]),
        payload["diagnostics"],
    )


def numerical_settings(settings):
    """Catalog/assembly changes do not invalidate identical radial fits.

    Increment fitting_protocol when changing adaptive selection semantics.
    The converter source hash already invalidates changes to radial fitting.
    """
    return {
        k: settings.get(k)
        for k in ("primitives", "derivative_weight", "limits", "fitter_sha256")
    } | {"fitting_protocol": settings.get("fitting_protocol", 1)}


def process_source(task):
    source, pao_dir, cache_dir, settings = task
    identity = {
        "source_sha256": source["sha256"],
        "selections": source["selections"],
        "settings": settings,
    }
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    cache_path = Path(cache_dir) / (source["stem"] + "-" + key[:20] + ".json")
    pao = PaoFile.read(Path(pao_dir) / (source["stem"] + ".pao"))
    if pao.sha256 != source["sha256"] or pao.element != source["element"]:
        raise ValueError(f"{source['stem']}: PAO does not match pinned catalog")
    if cache_path.exists():
        cached = json.loads(cache_path.read_text())
        if cached["identity"] != identity:
            raise ValueError("Fit cache identity mismatch")
        return cached["result"]
    specs = [BasisSpec.parse(s["label"]) for s in source["selections"]]
    required = sorted(
        {(l, count) for spec in specs for l, count in enumerate(spec.counts) if count}
    )
    channels, errors = {}, {}
    # Reuse individual channels when only catalog selections/provenance changed.
    # A PAO hash, fitting-code hash, thresholds and numerical settings must agree.
    for previous_path in sorted(Path(cache_dir).glob(source["stem"] + "-*.json")):
        previous = json.loads(previous_path.read_text())
        prior_identity = previous["identity"]
        if prior_identity["source_sha256"] == source["sha256"] and numerical_settings(
            prior_identity["settings"]
        ) == numerical_settings(settings):
            channels.update(previous["result"]["channels"])
            errors.update(previous["result"]["errors"])
    wanted = {f"{l}:{count}" for l, count in required}
    channels = {key: value for key, value in channels.items() if key in wanted}
    errors = {key: value for key, value in errors.items() if key in wanted}
    for l, count in required:
        channel_key = f"{l}:{count}"
        if channel_key in channels or channel_key in errors:
            continue
        attempts, candidates = [], []
        for nprimitive in settings["primitives"]:
            try:
                fit = fit_channel(
                    pao,
                    l,
                    count,
                    nprimitive=nprimitive,
                    derivative_weight=settings["derivative_weight"],
                )
                numeric, score, failed = channel_quality(
                    fit.diagnostics, settings["limits"]
                )
                attempts.append(
                    {
                        "nprimitive": nprimitive,
                        "numerically_valid": numeric,
                        "quality_score": score if math.isfinite(score) else None,
                        "failed_criteria": failed,
                    }
                )
                if numeric:
                    candidates.append((score, fit))
                if numeric and not failed:
                    break
            except (ValueError, np.linalg.LinAlgError) as exc:
                attempts.append({"nprimitive": nprimitive, "error": str(exc)})
        if candidates:
            _, selected = min(candidates, key=lambda item: item[0])
            selected.diagnostics["adaptive_attempts"] = attempts
            channels[channel_key] = fit_payload(selected)
        else:
            errors[channel_key] = attempts
    result = {
        "source_file": pao.path.name,
        "source_sha256": pao.sha256,
        "element": pao.element,
        "source_valence": pao.valence,
        "nominal_cutoff_bohr": pao.cutoff,
        "available_radial_functions": {
            str(l): values.shape[1] for l, values in pao.orbitals.items()
        },
        "channels": channels,
        "errors": errors,
    }
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    temporary = cache_path.with_suffix(f".partial-{os.getpid()}")
    temporary.write_text(
        json.dumps({"identity": identity, "result": result}, allow_nan=False)
    )
    temporary.replace(cache_path)
    return result


def assemble(catalog, results, settings, include_inaccurate):
    entries, reports, source_reports = [], [], []
    for source, result in zip(catalog["sources"], results):
        source_reports.append({k: v for k, v in result.items() if k != "channels"})
        pao = PaoFile(
            Path(result["source_file"]),
            source["element"],
            result["nominal_cutoff_bohr"],
            result["source_valence"],
            np.array([]),
            {},
            source["sha256"],
        )
        for selection in source["selections"]:
            spec = BasisSpec.parse(selection["label"])
            required = [f"{l}:{n}" for l, n in enumerate(spec.counts) if n]
            report = {
                **selection,
                "selection": spec.label,
                "cp2k_name": spec.name,
                "element": pao.element,
                "source_file": pao.path.name,
                "source_sha256": pao.sha256,
                "source_url": source["url"],
                "source_valence": pao.valence,
                "nominal_cutoff_bohr": pao.cutoff,
                "nao": spec.nao,
                "cp2k_from_openmx_indices_zero_based": spec.permutation(),
                "fallback_shell_policy": source["fallback_shell_policy"],
            }
            if any(key not in result["channels"] for key in required):
                report.update(
                    status="rejected_numerical",
                    exported=False,
                    failed_channels=[
                        key for key in required if key not in result["channels"]
                    ],
                )
                reports.append(report)
                continue
            fits = [restore_fit(result["channels"][key]) for key in required]
            failed = sorted(
                {
                    criterion
                    for fit in fits
                    for criterion in channel_quality(
                        fit.diagnostics, settings["limits"]
                    )[2]
                }
            )
            passed = not failed
            export = passed or include_inaccurate
            report.update(
                status="passes_fit_thresholds" if passed else "warning_fit_accuracy",
                exported=export,
                failed_criteria=failed,
            )
            if export:
                report["radial_fits"] = [fit.diagnostics for fit in fits]
                quality = (
                    "PASSES NUMERICAL FIT THRESHOLDS; production accuracy unvalidated"
                    if passed
                    else "WARNING: FAILS " + ", ".join(failed)
                )
                entries.append(
                    f"# QUALITY: {quality}\n" + render_basis(spec, pao, fits)
                )
            reports.append(report)
    exported = [r for r in reports if r["exported"]]
    missing_elements = sorted(
        set(catalog["elements"]) - {r["element"] for r in exported}, key=ELEMENTS.index
    )
    summary = {
        "catalog_elements": len(catalog["elements"]),
        "catalog_pao_files": len(catalog["sources"]),
        "catalog_selections": len(reports),
        "exported_basis_sets": len(exported),
        "exported_elements": len({r["element"] for r in exported}),
        "exported_pao_files": len({r["source_file"] for r in exported}),
        "status_counts": dict(Counter(r["status"] for r in reports)),
        "missing_elements": missing_elements,
    }
    policy = (
        "# WARNING entries fail numerical fit thresholds; included by explicit request.\n"
        if include_inaccurate
        else "# Only entries passing all numerical fit thresholds are included.\n"
    )
    header = (
        "# BASIS_OMX: approximate Gaussian fits to the linked OpenMX 2019 PAO database.\n"
        "# Generated by tools/openmx/build_database.py; source catalog and errors in BASIS_OMX.json.\n"
        "# OpenMX database: T. Ozaki and H. Kawai; distributed under GNU GPL, without warranty.\n"
        "# https://www.openmx-square.org/vps_pao2019/\n"
        "# EXPERIMENTAL: not lossless; no pseudopotential conversion; no production energy/force validation.\n"
        + policy
        + "# Default CP2K primitive/contraction normalization; OpenMX DM needs a basis transformation.\n"
        f"# Coverage: {summary['exported_elements']} elements, {summary['exported_pao_files']} PAO files, {len(exported)} selections.\n\n"
    )
    content = header + "\n".join(entries)
    report = {
        "schema": "openmx-cp2k-database-fit-v3",
        "exact_conversion": False,
        "basis_sha256": hashlib.sha256(content.encode()).hexdigest(),
        "source_catalog_sha256": settings["catalog_sha256"],
        "settings": settings,
        "selection_policy": catalog["selection_policy"],
        "include_inaccurate": include_inaccurate,
        "summary": summary,
        "sources": source_reports,
        "basis_sets": reports,
    }
    return content, report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--pao-dir", type=Path, required=True)
    parser.add_argument("--fit-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--primitives", type=int, nargs="+", default=[24, 40, 64])
    parser.add_argument("--include-inaccurate", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    report_path = args.output.with_name(args.output.name + ".json")
    if not args.overwrite and (args.output.exists() or report_path.exists()):
        parser.error("Output exists; use --overwrite explicitly")
    if not 1 <= args.workers <= 16 or any(not 4 <= n <= 80 for n in args.primitives):
        parser.error("Use 1..16 workers and 4..80 primitives")
    catalog = json.loads(args.manifest.read_text())
    if catalog.get("schema") != "openmx-2019-source-catalog-v1" or not catalog.get(
        "sources"
    ):
        parser.error("Invalid or empty source catalog")
    # Prevent manifest-supplied paths from escaping the PAO directory.
    for source in catalog["sources"]:
        stem = source["stem"]
        if Path(stem).name != stem or ".." in stem:
            parser.error("Invalid source stem")
        if any(BasisSpec.parse(s["label"]).stem != stem for s in source["selections"]):
            parser.error("Selection/source mismatch")
    inputs = {args.manifest.resolve()} | {
        (args.pao_dir / (s["stem"] + ".pao")).resolve() for s in catalog["sources"]
    }
    if args.output.resolve() in inputs or report_path.resolve() in inputs:
        parser.error("Output would overwrite a source file or manifest")
    settings = {
        "primitives": args.primitives,
        "derivative_weight": 0.2,
        "limits": DEFAULT_LIMITS,
        "fitter_sha256": hashlib.sha256(
            Path(convert_basis.__file__).read_bytes()
        ).hexdigest(),
        "builder_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "catalog_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
    }
    results = []
    tasks = [
        (s, str(args.pao_dir), str(args.fit_cache), settings)
        for s in catalog["sources"]
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for source, result in zip(catalog["sources"], pool.map(process_source, tasks)):
            results.append(result)
            warnings = sum(
                bool(channel_quality(c["diagnostics"], DEFAULT_LIMITS)[2])
                for c in result["channels"].values()
            )
            print(
                f"{len(results)}/{len(tasks)} {source['stem']}: "
                f"{len(result['channels'])} fitted channels, {warnings} accuracy warnings, {len(result['errors'])} failures",
                flush=True,
            )
    content, report = assemble(catalog, results, settings, args.include_inaccurate)
    args.output.write_text(content)
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report["summary"], indent=2))
    return (
        2
        if report["summary"]["missing_elements"]
        or report["summary"]["status_counts"].get("rejected_numerical")
        else 0
    )


if __name__ == "__main__":
    sys.exit(main())
