# SPDX-License-Identifier: GPL-2.0-or-later
"""Offline catalog, quality-gate and assembly tests."""

import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from build_database import (
    DEFAULT_LIMITS,
    assemble,
    channel_quality,
    fit_payload,
    process_source,
)
from collect_database import ROOT, choices, database_link, family
from convert_basis import Fit, PaoFile, BasisSpec, render_basis
from validate_database import read_basis
from test_convert_basis import fixture


class DatabaseTests(unittest.TestCase):
    def test_links_are_scoped(self):
        page = ROOT + "Fe/index.html"
        self.assertEqual(
            database_link(page, "Fe_Hard/index.html"), ROOT + "Fe/Fe_Hard/index.html"
        )
        self.assertEqual(database_link(page, "../C/C6.0.pao#data"), ROOT + "C/C6.0.pao")
        self.assertIsNone(database_link(page, "https://example.org/file.pao"))
        self.assertIsNone(database_link(page, "../../download.html"))
        self.assertIsNone(
            database_link(page, "https://www.openmx-square.org/elsewhere/C.pao")
        )

    def test_families_and_published_patterns(self):
        self.assertEqual(family("Fe5.5H"), ("Fe", "H"))
        self.assertEqual(family("Nd10.0_OC"), ("Nd", "_OC"))
        patterns = choices(
            "<td>Fe5.5H-s3p2d2f1</td> Fe6.0S-s2p2d1 La*.*-s2p1f1", "source"
        )
        self.assertEqual(
            {(c["element"], c["family_suffix"], c["shells"]) for c in patterns},
            {("Fe", "H", "s3p2d2f1"), ("Fe", "S", "s2p2d1"), ("La", "", "s2p1f1")},
        )

    @staticmethod
    def diagnostics(error=0.001, norm=0.0):
        return {
            "l": 0,
            "n_radial": 1,
            "normalized_radial_l2_errors": [error],
            "kinetic_absolute_errors_hartree": [0.0],
            "same_center_overlap_max_error": 0.0,
            "analytic_vs_quadrature_norm_max_error": norm,
        }

    def test_distinguish_inaccuracy_from_numerical_failure(self):
        valid, score, failed = channel_quality(self.diagnostics(), DEFAULT_LIMITS)
        self.assertTrue(valid)
        self.assertEqual(failed, [])
        valid, score, failed = channel_quality(
            self.diagnostics(error=0.2), DEFAULT_LIMITS
        )
        self.assertTrue(valid)
        self.assertGreater(score, 1)
        self.assertEqual(failed, ["radial_l2"])
        self.assertFalse(channel_quality(self.diagnostics(norm=0.1), DEFAULT_LIMITS)[0])

    def test_assembly_warnings_and_missing_elements(self):
        source = {
            "element": "C",
            "stem": "C6.0",
            "sha256": "abc",
            "url": "source",
            "fallback_shell_policy": False,
            "selections": [{"label": "C6.0-s1", "shells": "s1"}],
        }
        catalog = {"elements": ["C"], "sources": [source], "selection_policy": "test"}
        result = {
            "source_file": "C6.0.pao",
            "source_sha256": "abc",
            "source_valence": 4,
            "nominal_cutoff_bohr": 6,
            "errors": {},
            "channels": {
                "0:1": fit_payload(
                    Fit(
                        0,
                        np.array([1.0]),
                        np.array([[1.0]]),
                        self.diagnostics(error=0.2),
                    )
                )
            },
        }
        settings = {"limits": DEFAULT_LIMITS, "catalog_sha256": "def"}
        text, report = assemble(catalog, [result], settings, False)
        self.assertNotIn("C OMX-FIT", text)
        self.assertIn("Only entries passing all numerical fit thresholds", text)
        self.assertNotIn("radial_fits", report["basis_sets"][0])
        self.assertEqual(report["summary"]["missing_elements"], ["C"])
        text, report = assemble(catalog, [result], settings, True)
        self.assertIn("# QUALITY: WARNING: FAILS radial_l2", text)
        self.assertIn("C OMX-FIT-C6.0-s1", text)
        self.assertEqual(report["summary"]["exported_elements"], 1)
        self.assertIn("radial_fits", report["basis_sets"][0])

    def test_shipped_library_passes_all_thresholds(self):
        root = Path(__file__).resolve().parents[2]
        basis = root / "data" / "BASIS_OMX"
        report = json.loads(basis.with_suffix(".json").read_text())
        catalog = Path(__file__).parent / "openmx2019_sources.json"
        self.assertEqual(
            hashlib.sha256(catalog.read_bytes()).hexdigest(),
            report["source_catalog_sha256"],
        )
        self.assertEqual(
            hashlib.sha256(basis.read_bytes()).hexdigest(), report["basis_sha256"]
        )
        for key, filename in (
            ("fitter_sha256", "convert_basis.py"),
            ("builder_sha256", "build_database.py"),
        ):
            self.assertEqual(
                hashlib.sha256(
                    (Path(__file__).parent / filename).read_bytes()
                ).hexdigest(),
                report["settings"][key],
            )
        self.assertFalse(report["include_inaccurate"])
        entries = read_basis(basis)
        exported = {b["cp2k_name"]: b for b in report["basis_sets"] if b["exported"]}
        self.assertEqual(entries.keys(), exported.keys())
        self.assertEqual(len(entries), 279)
        self.assertEqual(len({b["element"] for b in entries.values()}), 81)
        for name, entry in entries.items():
            self.assertEqual(entry["nao"], exported[name]["nao"])
            self.assertEqual(entry["element"], exported[name]["element"])
            self.assertEqual(exported[name]["status"], "passes_fit_thresholds")
            for channel in exported[name]["radial_fits"]:
                valid, _, failed = channel_quality(channel, DEFAULT_LIMITS)
                self.assertTrue(valid)
                self.assertEqual(failed, [])

    def test_pinned_source_and_cache(self):
        with tempfile.TemporaryDirectory(prefix="basis-catalog-test-") as temp:
            root = Path(temp)
            fixture(root / "C6.0.pao")
            pao = PaoFile.read(root / "C6.0.pao")
            source = {
                "stem": "C6.0",
                "sha256": pao.sha256,
                "element": "C",
                "selections": [{"label": "C6.0-s1"}],
            }
            settings = {
                "limits": DEFAULT_LIMITS,
                "primitives": [24],
                "derivative_weight": 0.2,
            }
            task = source, str(root), str(root / "fits"), settings
            first = process_source(task)
            self.assertEqual(first, process_source(task))
            self.assertEqual(len(list((root / "fits").glob("*.json"))), 1)
            # A catalog/provenance edit must not trigger an expensive refit.
            with patch(
                "build_database.fit_channel",
                side_effect=AssertionError("unexpected refit"),
            ):
                updated = (
                    source,
                    str(root),
                    str(root / "fits"),
                    {**settings, "catalog_sha256": "new catalog"},
                )
                self.assertEqual(first, process_source(updated))
            source["sha256"] = "wrong"
            with self.assertRaises(ValueError):
                process_source(task)

    def test_pinned_catalog(self):
        catalog = json.loads(
            (Path(__file__).parent / "openmx2019_sources.json").read_text()
        )
        self.assertEqual(len(catalog["elements"]), 81)
        self.assertEqual(len(catalog["sources"]), 268)
        names = [
            selection["label"]
            for source in catalog["sources"]
            for selection in source["selections"]
        ]
        self.assertEqual(len(names), 1415)
        self.assertEqual(len(names), len(set(names)))
        self.assertTrue(
            all(not source["fallback_shell_policy"] for source in catalog["sources"])
        )
        self.assertIn("Nd8.0_OC-s2p1d1", names)
        self.assertIn("Cu6.0S-s1d1", names)
        self.assertIn("O7.0-s2p2d1", names)
        self.assertIn("H7.0-s2p1", names)
        for source in catalog["sources"]:
            self.assertEqual(len(source["sha256"]), 64)
            self.assertEqual(database_link(ROOT, source["url"]), source["url"])
            self.assertTrue(
                all(
                    BasisSpec.parse(choice["label"]).stem == source["stem"]
                    for choice in source["selections"]
                )
            )

    def test_independent_basis_reader(self):
        with tempfile.TemporaryDirectory(prefix="basis-reader-unit-") as temp:
            path = Path(temp) / "basis"
            pao = PaoFile(Path("C6.0.pao"), "C", 6.0, 4.0, np.array([]), {}, "abc")
            fit = Fit(0, np.array([1.0]), np.array([[1.0]]), self.diagnostics())
            path.write_text(render_basis(BasisSpec.parse("C6.0-s1"), pao, [fit]))
            self.assertEqual(read_basis(path)["OMX-FIT-C6.0-s1"]["nao"], 1)
            fit.coefficients[0, 0] = 2
            path.write_text(render_basis(BasisSpec.parse("C6.0-s1"), pao, [fit]))
            with self.assertRaises(ValueError):
                read_basis(path)


if __name__ == "__main__":
    unittest.main()
