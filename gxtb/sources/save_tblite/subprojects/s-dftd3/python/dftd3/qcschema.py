# This file is part of s-dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.
"""
QCSchema Support
================

Integration with the `QCArchive infrastructure <http://docs.qcarchive.molssi.org>`_.

If the QCElemental package is installed the ``dftd3.qcschema`` module becomes
importable and provides the ``run_qcschema`` function supporting QCSchema v1.
If the QCElemental package is >=0.50.0, ``dftd3.qcschema`` supports QCSchema v1
and v2, returning whichever version was submitted. Note that Python 3.14+ only
works with QCSchema v2 due to Pydantic restrictions.

This module provides a way to translate QCSchema or QCElemental Atomic Input
into a format understandable by the ``dftd3`` API which in turn provides the
calculation results in a QCSchema compatible format.

Supported keywords are

======================== =========== ============================================
 Keyword                  Default     Description
======================== =========== ============================================
 level_hint               None        Dispersion correction level
 params_tweaks            None        Optional dict with the damping parameters
 pair_resolved            False       Enable pairwise resolved dispersion energy
======================== =========== ============================================

Allowed level hints are ``"d3bj"``, ``"d3zero"``, ``"d3bjm"``/``"d3mbj"``,
``"d3mzero"``/``"d3zerom"``, ``"d3op"``, and ``"d3cso"``.

The params_tweaks dict contains the damping parameters, at least s8, a1 and a2
must be provided for rational damping, while s8 and rs6 are required in case
of zero damping. For CSO damping, a1 must be provided.

Parameters for (modified) rational damping are:

======================== =========== ============================================
 Tweakable parameter      Default     Description
======================== =========== ============================================
 s6                       1.0         Scaling of the dipole-dipole dispersion
 s8                       None        Scaling of the dipole-quadrupole dispersion
 s9                       1.0         Scaling of the three-body dispersion energy
 a1                       None        Scaling of the critical radii
 a2                       None        Offset of the critical radii
 alp                      14.0        Exponent of the zero damping (ATM only)
======================== =========== ============================================

Parameters for (modified) zero damping are:

======================== =========== ===================================================
 Tweakable parameter      Default     Description
======================== =========== ===================================================
 s6                       1.0         Scaling of the dipole-dipole dispersion
 s8                       None        Scaling of the dipole-quadrupole dispersion
 s9                       1.0         Scaling of the three-body dispersion energy
 rs6                      None        Scaling of the dipole-dipole damping
 rs8                      1.0         Scaling of the dipole-quadrupole damping
 alp                      14.0        Exponent of the zero damping
 bet                      None        Offset for damping radius (modified zero damping)
======================== =========== ===================================================

Parameters for optimized power damping are:

======================== =========== ============================================
 Tweakable parameter      Default     Description
======================== =========== ============================================
 s6                       1.0         Scaling of the dipole-dipole dispersion
 s8                       None        Scaling of the dipole-quadrupole dispersion
 s9                       1.0         Scaling of the three-body dispersion energy
 a1                       None        Scaling of the critical radii
 a2                       None        Offset of the critical radii
 alp                      14.0        Exponent of the zero damping (ATM only)
 bet                      None        Power for the zero-damping component
======================== =========== ============================================

Parameters for CSO (C6-scaled only) damping are:

======================== =========== ============================================
 Tweakable parameter      Default     Description
======================== =========== ============================================
 s6                       1.0         Scaling of the dipole-dipole dispersion
 s9                       1.0         Scaling of the three-body dispersion energy
 a1                       None        Sigmoid amplitude parameter
 a2                       2.5         Sigmoid reference distance scale
 a3                       0.0         Denominator critical radii scale
 a4                       6.25        Denominator constant offset
 alp                      14.0        Exponent of the zero damping (ATM only)
======================== =========== ============================================

.. note::

    input_data.model.method with a full method name and input_data.keywords["params_tweaks"]
    cannot be provided at the same time. It is an error to provide both options at the
    same time.

Example
-------

>>> from dftd3.qcschema import run_qcschema
>>> import qcelemental as qcel
>>> atomic_input = qcel.models.AtomicInput(
...     molecule = qcel.models.Molecule(
...         symbols = ["O", "H", "H"],
...         geometry = [
...             0.00000000000000,  0.00000000000000, -0.73578586109551,
...             1.44183152868459,  0.00000000000000,  0.36789293054775,
...            -1.44183152868459,  0.00000000000000,  0.36789293054775
...         ],
...     ),
...     driver = "energy",
...     model = {
...         "method": "TPSS-D3(BJ)",
...     },
...     keywords = {},
... )
...
>>> atomic_result = run_qcschema(atomic_input)
>>> atomic_result.return_result
-0.00042042440936212056
"""

import sys
from typing import Union, overload
from .interface import (
    DispersionModel,
    RationalDampingParam,
    ZeroDampingParam,
    ModifiedRationalDampingParam,
    ModifiedZeroDampingParam,
    OptimizedPowerDampingParam,
    CSODampingParam,
)
from .library import get_api_version
import numpy as np

if sys.version_info < (3, 14):
    try:
        import qcelemental.models.v1 as qcel_v1
    except ModuleNotFoundError:
        import qcelemental.models as qcel_v1
else:
    qcel_v1 = None

try:
    import qcelemental.models.v2 as qcel_v2
except ModuleNotFoundError:
    qcel_v2 = None


if qcel_v1 is None and qcel_v2 is None:
    raise ModuleNotFoundError(
        "The qcelemental package is required for qcschema support. "
        "Please install it with 'pip install qcelemental'."
    )


_supported_drivers = [
    "energy",
    "gradient",
]

_available_levels = [
    "d3bj",
    "d3zero",
    "d3bjm",
    "d3mbj",
    "d3zerom",
    "d3mzero",
    "d3op",
    "d3cso",
]

_damping_param = {
    "d3bj": RationalDampingParam,
    "d3zero": ZeroDampingParam,
    "d3bjm": ModifiedRationalDampingParam,
    "d3mbj": ModifiedRationalDampingParam,
    "d3zerom": ModifiedZeroDampingParam,
    "d3mzero": ModifiedZeroDampingParam,
    "d3op": OptimizedPowerDampingParam,
    "d3cso": CSODampingParam,
}

_clean_dashlevel = str.maketrans("", "", "()")


def error_return_result(driver, molecule):
    if driver == "energy":
        return 0.0
    if driver == "gradient":
        natoms = len(molecule.symbols)
        return np.zeros((natoms, 3))
    return None


if qcel_v1 is not None:

    @overload
    def run_qcschema(
        input_data: Union[dict, "qcel_v1.AtomicInput"],
    ) -> "qcel_v1.AtomicResult": ...


if qcel_v2 is not None:

    @overload
    def run_qcschema(
        input_data: Union[dict, "qcel_v2.AtomicInput"],
    ) -> "qcel_v2.AtomicResult": ...


def run_qcschema(input_data):
    """Perform disperson correction based on an atomic inputmodel"""

    if qcel_v2 is not None and isinstance(input_data, qcel_v2.AtomicInput):
        atomic_input = input_data
    elif qcel_v1 is not None and isinstance(input_data, qcel_v1.AtomicInput):
        atomic_input = input_data
    elif qcel_v2 is not None and input_data.get("specification"):
        atomic_input = qcel_v2.AtomicInput(**input_data)
    elif qcel_v1 is not None:
        atomic_input = qcel_v1.AtomicInput(**input_data)
    else:
        raise ValueError(
            "Input data is not a valid QCSchema AtomicInput for either v1 or v2."
        )

    schema_version = atomic_input.schema_version
    if schema_version == 1:
        ret_data = atomic_input.dict()
        input_keywords = atomic_input.keywords
        input_method = atomic_input.model.method
        input_driver = atomic_input.driver
    elif schema_version == 2:
        ret_data = {
            "input_data": atomic_input,
            "extras": {},
            "molecule": atomic_input.molecule,
        }
        input_keywords = atomic_input.specification.keywords
        input_method = atomic_input.specification.model.method
        input_driver = atomic_input.specification.driver
    else:
        raise ValueError(
            f"Unsupported QCSchema version: {schema_version}. Only v1 and v2 are supported."
        )

    provenance = {
        "creator": "s-dftd3",
        "version": get_api_version(),
        "routine": "dftd3.qcschema.run_qcschema",
    }
    success = False
    return_result = 0.0
    properties = {}

    # Since it is a level hint we a forgiving if it is not present,
    # we are much less forgiving if the wrong level is hinted here.
    _level = input_keywords.get("level_hint", "d3bj")
    if _level.lower() not in _available_levels:
        error = dict(
            error_type="input error",
            error_message="Level '{}' is invalid for this dispersion correction".format(
                _level
            ),
        )
        if schema_version == 1:
            ret_data.update(
                provenance=provenance,
                success=success,
                properties=properties,
                return_result=return_result,
                error=qcel_v1.ComputeError(**error),
            )
            return qcel_v1.AtomicResult(**ret_data)
        elif schema_version == 2:
            return qcel_v2.FailedOperation(
                input_data=atomic_input, error=qcel_v2.ComputeError(**error)
            )

    # Check if the method is provided and strip the “dashlevel” from the method
    _method = input_method.split("-")
    if _method[-1].lower().translate(_clean_dashlevel) == _level.lower():
        _method.pop()
    _method = "-".join(_method)
    if len(_method) == 0:
        _method = None

    # Obtain the parameters for the damping function
    _input_param = input_keywords.get("params_tweaks", {"method": _method})

    try:
        param = _damping_param[_level](
            **_input_param,
        )

        disp = DispersionModel(
            atomic_input.molecule.atomic_numbers[atomic_input.molecule.real],
            atomic_input.molecule.geometry[atomic_input.molecule.real],
        )

        driver = input_driver
        res = disp.get_dispersion(
            param=param,
            grad=driver == "gradient",
        )
        extras = {"dftd3": res}

        if driver == "gradient":
            if all(atomic_input.molecule.real):
                fullgrad = res.get("gradient")
            else:
                ireal = np.argwhere(atomic_input.molecule.real).reshape((-1))
                fullgrad = np.zeros_like(atomic_input.molecule.geometry)
                fullgrad[ireal, :] = res.get("gradient")

        properties.update(return_energy=res.get("energy"))

        if input_keywords.get("pair_resolved", False):
            res = disp.get_pairwise_dispersion(param=param)
            extras["dftd3"].update(res)

        success = driver in _supported_drivers
        if driver == "energy":
            return_result = properties["return_energy"]
        elif driver == "gradient":
            return_result = fullgrad
        else:
            ret_data.update(
                error=dict(
                    error_type="input error",
                    error_message="Calculation succeeded but invalid driver request provided",
                ),
            )

        ret_data["extras"].update(extras)

    except (RuntimeError, TypeError) as e:
        (
            ret_data.update(
                error=dict(error_type="input error", error_message=str(e)),
            ),
        )

    ret_data.update(
        provenance=provenance,
        success=success,
        properties=properties,
        return_result=return_result,
    )

    if schema_version == 1:
        return qcel_v1.AtomicResult(**ret_data)
    
    if "error" in ret_data:
        return qcel_v2.FailedOperation(
            input_data=atomic_input, error=ret_data["error"]
        )
    return qcel_v2.AtomicResult(**ret_data)
