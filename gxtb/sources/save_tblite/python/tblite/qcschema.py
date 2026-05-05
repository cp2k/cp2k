# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.
"""
Integration with the `QCArchive infrastructure <https://qcarchive.molssi.org>`_.

This module provides a function to run QCSchema input or QCElemental Atomic Input
through the TBLite calculator and return the result as a QCElemental Atomic Result.

The model support the following methods:

- **GFN2-xTB**: 
  Self-consistent extended tight binding Hamiltonian with
  anisotropic second order electrostatic contributions,
  third order on-site contributions and self-consistent D4 dispersion.
  Geometry, frequency and non-covalent interactions parametrisation for
  elements up to Z=86.

- **GFN1-xTB**:
  Self-consistent extended tight binding Hamiltonian with
  isotropic second order electrostatic contributions and
  third order on-site contributions.
  Geometry, frequency and non-covalent interactions parametrisation for
  elements up to Z=86.

- **IPEA1-xTB**:
  Special parametrisation for the GFN1-xTB Hamiltonian to improve the
  description of vertical ionisation potentials and electron affinities.
  Uses additional diffuse s-functions on light main group elements.
  Parametrised up to Z=86.

Supported keywords are:

=================== ==================================== =========================================
 name                description                          default
=================== ==================================== =========================================
 accuracy            Numerical thresholds for SCC         float (1.0)
 guess               Initial guess for wavefunction       integer (0 == SAD)
 max-iter            Maximum number of SCC iterations     integer (250)
 mixer-damping       Parameter for the SCC mixer          float (0.4)
 save-integrals      Keep integral matrices in results    0 (False)
 temperature         Electronic temperature for filling   float (9.500e-4)
 verbosity           Set verbosity of printout            integer (1)
 electric-field      Uniform electric field               Field vector
 spin-polarization   Spin polarization                    Scaling factor
 alpb-solvation      ALPB implicit solvation              Solvent name, solution state (optional)
 gbsa-solvation      GBSA implicit solvation              Solvent name, solution state (optional)
 cpcm-solvation      CPCM implicit solvation              Epsilon
 gbe-solvation       GBε implicit solvation               Epsilon, Born kernel
 gb-solvation        GB implicit solvation                Epsilon, Born kernel
=================== ==================================== =========================================
"""

from io import StringIO
from typing import Any, Dict, Union

import numpy as np
import qcelemental as qcel

from .exceptions import TBLiteRuntimeError, TBLiteTypeError, TBLiteValueError
from .interface import Calculator
from .library import get_version

SUPPORTED_DRIVERS = {
    "energy",
    "gradient",
    # "hessian",
    "properties",
}


def get_provenance() -> qcel.models.Provenance:
    """
    Returns a QCSchema provenance model.

    Returns
    -------
    qcel.models.Provenance
        A QCSchema provenance model.
    """
    return qcel.models.Provenance(
        creator="tblite",
        version=".".join([str(v) for v in get_version()]),
        routine="tblite.qcschema.run_schema",
    )


def get_error(
    input_data: qcel.models.AtomicInput,
    error: Union[Dict[str, Any], qcel.models.ComputeError],
) -> qcel.models.AtomicResult:
    if not isinstance(error, qcel.models.ComputeError):
        error = qcel.models.ComputeError(**error)

    return_data = input_data.dict()
    return_data.update(
        error=error,
        success=False,
        return_result={
            "energy": 0.0,
            "gradient": np.zeros(input_data.molecule.geometry.shape),
            "hessian": np.zeros(
                (
                    input_data.molecule.geometry.size,
                    input_data.molecule.geometry.size,
                )
            ),
            "properties": {},
        }[input_data.driver],
        properties={},
        provenance=get_provenance(),
    )

    return qcel.models.AtomicResult(**return_data)


class _Logger:
    def __init__(self):
        self._buffer = StringIO()

    def __call__(self, message: str) -> None:
        print(message)
        self._buffer.write(message + "\n")

    def __str__(self) -> str:
        return self._buffer.getvalue()


def run_schema(
    input_data: Union[Dict[str, Any], qcel.models.AtomicInput]
) -> qcel.models.AtomicResult:
    """Runs a QCSchema input through the QCEngine stack and returns a QCSchema result.

    Parameters
    ----------
    input_data : qcel.models.AtomicInput
        A QCSchema input dictionary or model.

    Returns
    -------
    qcel.models.AtomicResult
        A QCSchema result model.
    """
    if not isinstance(input_data, qcel.models.AtomicInput):
        input_data = qcel.models.AtomicInput(**input_data)

    if input_data.driver not in SUPPORTED_DRIVERS:
        driver_name = (
            input_data.driver.name
            if hasattr(input_data.driver, "name")
            else str(input_data.driver)
        )
        return get_error(
            input_data,
            qcel.models.ComputeError(
                error_type="input_error",
                error_message=f"Driver '{driver_name}' is not supported by tblite.",
            ),
        )

    if input_data.model.method not in Calculator._loader:
        return get_error(
            input_data,
            qcel.models.ComputeError(
                error_type="input_error",
                error_message=f"Model '{input_data.model.method}' is not supported by tblite.",
            ),
        )

    if input_data.model.basis is not None:
        return get_error(
            input_data,
            qcel.models.ComputeError(
                error_type="input_error",
                error_message="Basis sets are not supported by tblite.",
            ),
        )

    keywords = {
        key: value
        for key, value in input_data.keywords.items()
        if key in Calculator._setter
    }
    interaction = {
        key: value
        for key, value in input_data.keywords.items()
        if key in Calculator._interaction
    }
    unknown_keywords = (
        set(input_data.keywords) - set(keywords) - set(interaction)
    )
    if unknown_keywords:
        return get_error(
            input_data,
            qcel.models.ComputeError(
                error_type="input_error",
                error_message=f"Unknown keywords: {', '.join(unknown_keywords)}.",
            ),
        )

    logger = _Logger()
    try:
        calc = Calculator(
            method=input_data.model.method,
            numbers=input_data.molecule.atomic_numbers,
            positions=input_data.molecule.geometry,
            charge=input_data.molecule.molecular_charge,
            uhf=input_data.molecule.molecular_multiplicity - 1,
            color=False,
            logger=logger,
        )
        for key, value in keywords.items():
            calc.set(key, value)

        for key, value in interaction.items():
            if not isinstance(value, list):
                value = [value]
            calc.add(key, *value)

        result = calc.singlepoint().dict()

        properties = qcel.models.AtomicResultProperties(
            return_energy=result["energy"],
            return_gradient=result["gradient"],
            calcinfo_natom=result["natoms"],
            calcinfo_nbasis=result["norbitals"],
            calcinfo_nmo=result["norbitals"],
            scf_dipole_moment=result["dipole"],
            scf_quadrupole_moment=result["quadrupole"][
                [0, 1, 3, 1, 2, 4, 3, 4, 5]
            ],
            scf_total_energy=result["energy"],
            scf_total_gradient=result["gradient"],
        )

        return_result = {
            "energy": result["energy"],
            "gradient": result["gradient"],
            "properties": {
                "dipole": result["dipole"],
                "mulliken_charges": result["charges"],
                "mayer_indices": result["bond-orders"],
            },
        }[input_data.driver]

    except (TBLiteTypeError, TBLiteValueError) as e:
        return get_error(
            input_data,
            qcel.models.ComputeError(
                error_type="input_error",
                error_message=str(e),
            ),
        )

    except TBLiteRuntimeError as e:
        return get_error(
            input_data,
            qcel.models.ComputeError(
                error_type="execution_error",
                error_message=str(e),
            ),
        )

    return_data = input_data.dict()
    return_data.update(
        success=True,
        return_result=return_result,
        properties=properties,
        provenance=get_provenance(),
        stdout=str(logger),
    )

    return qcel.models.AtomicResult(**return_data)
