"""Optional ASE calculator using libcp2k in the Python process."""

# SPDX-License-Identifier: GPL-2.0-or-later

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
import os
from pathlib import Path
from typing import Any, cast

import numpy as np
from ase import Atoms
from ase.calculators.calculator import CalculationFailed, Calculator, all_changes

from ._units import BOHR_TO_ANGSTROM as Bohr, HARTREE_TO_EV as Hartree
from .input import InputMapping, input_to_string
from .library import CP2K, ForceEnvironment, PathLike, SCFConvergenceError


def _uppercase(tree: InputMapping) -> dict[str, Any]:
    if not isinstance(tree, Mapping):
        raise TypeError("inp must be a CP2K input mapping")
    result: dict[str, Any] = {}
    for key, value in tree.items():
        if not isinstance(key, str) or key.upper() in result:
            raise ValueError(f"Invalid or duplicate input key: {key!r}")
        if isinstance(value, Mapping):
            value = _uppercase(value)
        elif isinstance(value, (list, tuple)):
            value = [
                _uppercase(item) if isinstance(item, Mapping) else item
                for item in value
            ]
        result[key.upper()] = deepcopy(value)
    return result


class CP2KCalculator(Calculator):
    """ASE energy/force calculator with no shell or output parsing.

    ``runtime`` is an open :class:`cp2k.CP2K` session, owned by the caller.
    ``inp`` specifies the method, DFT/force-field settings and SUBSYS/KINDs.
    Geometry and periodicity come exclusively from ASE's Atoms. A finite,
    right-handed cell is required even for nonperiodic systems. Positions and
    cells are converted from angstrom, energy/forces back to eV and eV/angstrom.

    Per-atom energies, shell models, topology files, and implicit use
    of ASE initial charges/magnetic moments are not supported. Set charge and
    spin explicitly in inp. Stress uses ASE's tensile-positive eV/angstrom**3
    convention and enables cell optimization/NPT. Close before closing runtime.
    """

    implemented_properties = ["energy", "free_energy", "forces", "stress"]
    default_parameters = {"inp": None}

    def __init__(
        self,
        runtime: CP2K,
        inp: InputMapping,
        *,
        output_file: PathLike | None = None,
        **kwargs: Any,
    ) -> None:
        self._runtime = runtime
        self._env: ForceEnvironment | None = None
        self._closed = False
        self._output_file = output_file
        # ASE's base-class methods do not yet declare their argument types.
        super().__init__(inp=deepcopy(inp), **kwargs)  # type: ignore[no-untyped-call]

    def set(self, **kwargs: Any) -> dict[str, Any]:
        unknown = set(kwargs) - {"inp"}
        if unknown:
            raise TypeError(f"Unknown CP2K calculator parameters: {sorted(unknown)}")
        if "inp" in kwargs:
            # Validate without mutating caller-owned dictionaries.
            kwargs["inp"] = _uppercase(kwargs["inp"])
            input_to_string(kwargs["inp"])
        changed = super().set(**kwargs)  # type: ignore[no-untyped-call]
        if changed:
            self._discard_environment()
            self.reset()  # type: ignore[no-untyped-call]
        return cast(dict[str, Any], changed)

    def _discard_environment(self) -> None:
        if self._env is not None:
            self._env.close()
            self._env = None

    def _make_input(self, atoms: Atoms) -> dict[str, Any]:
        tree: dict[str, Any] = deepcopy(self.parameters["inp"])
        force_eval = tree.setdefault("FORCE_EVAL", {})
        if not isinstance(force_eval, dict):
            raise ValueError("ASE requires exactly one FORCE_EVAL section")
        force_eval.setdefault("STRESS_TENSOR", "ANALYTICAL")
        subsys = force_eval.setdefault("SUBSYS", {})
        if not isinstance(subsys, dict):
            raise ValueError("SUBSYS must be a single section")
        conflicting = {
            "CELL",
            "COORD",
            "TOPOLOGY",
            "VELOCITY",
            "SHELL",
            "CORE",
        } & subsys.keys()
        if conflicting:
            raise ValueError(
                f"ASE owns the geometry; remove SUBSYS {sorted(conflicting)}"
            )
        subsys["CELL"] = {
            axis: row.tolist() for axis, row in zip("ABC", atoms.cell.array)
        }
        subsys["CELL"]["PERIODIC"] = (
            "".join(axis for axis, periodic in zip("XYZ", atoms.pbc) if periodic)
            or "NONE"
        )
        subsys["COORD"] = {
            "_lines": [
                f"{symbol} {x:.17g} {y:.17g} {z:.17g}"
                for symbol, (x, y, z) in zip(
                    atoms.get_chemical_symbols(), atoms.positions
                )
            ]
        }
        glob = tree.setdefault("GLOBAL", {})
        if not isinstance(glob, dict):
            raise ValueError("GLOBAL must be a single section")
        if "PROJECT" not in glob:
            # GLOBAL/PROJECT is limited to 80 bytes in native CP2K. Avoid
            # unnecessarily embedding a potentially long absolute cwd.
            project = os.path.relpath(Path(self.directory) / (self.prefix or "cp2k"))
            if '"' in project or len(os.fsencode(project)) > 80:
                raise ValueError(
                    "Use a shorter calculator directory/label or set GLOBAL/PROJECT "
                    "explicitly; generated PROJECT must fit 80 bytes without quotes"
                )
            glob["PROJECT"] = f'"{project}"'
        glob.setdefault("PRINT_LEVEL", "LOW")
        return tree

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: Sequence[str] = ("energy",),
        system_changes: Sequence[str] = all_changes,
    ) -> None:
        if self._closed:
            raise RuntimeError("This CP2K calculator is closed")
        self._runtime._check()
        super().calculate(atoms, properties, system_changes)  # type: ignore[no-untyped-call]
        self.results = {}
        assert self.atoms is not None
        if len(self.atoms) == 0:
            raise ValueError("CP2K requires at least one atom")
        cell = self.atoms.cell.array
        determinant = np.linalg.det(cell)
        if (
            not np.isfinite(cell).all()
            or not np.isfinite(determinant)
            or determinant <= 0
        ):
            raise ValueError("CP2K requires a finite, right-handed 3D cell")
        if not np.isfinite(self.atoms.positions).all():
            raise ValueError("Atomic positions must be finite")
        if np.any(self.atoms.get_initial_charges()) or np.any(
            self.atoms.get_initial_magnetic_moments()
        ):
            raise ValueError(
                "Specify charge/spin in inp, not ASE initial charges/magmoms"
            )
        if set(system_changes) & {"numbers", "pbc"}:
            self._discard_environment()
        if self._env is None:
            source = self._make_input(self.atoms)
            output = self._output_file or str(
                Path(self.directory) / f"{self.prefix or 'cp2k'}.out"
            )
            self._env = self._runtime.create_force_env(source, output_file=output)
            if self._env.natom != len(self.atoms) or self._env.nparticle != len(
                self.atoms
            ):
                self._discard_environment()
                raise ValueError("ASE requires exactly one CP2K particle per atom")
        self._env.cell = cell / Bohr
        self._env.positions = self.atoms.positions / Bohr
        try:
            result = self._env.calculate(stress="stress" in properties)
        except SCFConvergenceError as error:
            raise CalculationFailed(str(error)) from error
        assert result.forces is not None
        # As in ASE's shell calculator, use CP2K's variational total energy.
        self.results = {
            "energy": result.energy * Hartree,
            "free_energy": result.energy * Hartree,
            "forces": result.forces * (Hartree / Bohr),
        }
        if result.stress is not None:
            tensor = -result.stress * Hartree / Bohr**3
            self.results["stress"] = tensor.flat[[0, 4, 8, 5, 2, 1]]

    def close(self) -> None:
        if not self._closed:
            self._discard_environment()
            self.reset()  # type: ignore[no-untyped-call]
            self._closed = True

    def __enter__(self) -> CP2KCalculator:
        if self._closed:
            raise RuntimeError("This CP2K calculator is closed")
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()
