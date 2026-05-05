# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.
"""
QCSchema Support
----------------

Integration with the `QCArchive infrastructure <http://docs.qcarchive.molssi.org>`_.

This module provides a way to translate QCSchema or QCElemental Atomic Input
into a format understandable by the ``dftd4`` API which in turn provides the
calculation results in a QCSchema compatible format.

Supported keywords are

======================== =========== ============================================
 Keyword                  Default     Description
======================== =========== ============================================
 level_hint               None        Dispersion correction level ("d4" or "d4s")
 damping_hint             None        Optional dict with the damping functions
 params_tweaks            None        Optional dict with the damping parameters
 pair_resolved            False       Enable pairwise resolved dispersion energy
 property                 False       Evaluate dispersion related properties
======================== =========== ============================================

The damping_hint dict contains the two-body ("2b") and three-body ("3b") damping 
function types. If not provided the default damping functions of the moldel will
be used. If the three-body damping function is set to "none" the ATM contribution
will be disabled. 

The optional params_tweaks dict contains the damping parameters (either all
damping parameters or parameters for the model specific default (e.g., rational 
+ zero-avg for d4 requiring at least s8, a1 and a2). Additionally, parameter for
the dispersion models can be specified. The parameters are: 

======================== =========== ============================================
 Tweakable parameter      Default     Description
======================== =========== ============================================
 s6                       1.0         Scaling of the dipole-dipole dispersion
 s8                       None        Scaling of the dipole-quadrupole dispersion
 s9                       1.0         Scaling of the three-body dispersion energy
 a1                       None        Scaling of the critical radii
 a2                       None        Offset of the critical radii
 a3                       None        (Advanced) Additional damping parameter
 a4                       None        (Advanced) Additional damping parameter
 rs6                      None        (Advanced) Radii scaling
 rs8                      None        (Advanced) Radii scaling
 rs9                      None/1.0    (Advanced) Radii scaling for three-body
 alp                      None/16.0   Exponent of the zero damping (ATM only)
 bet                      None        (Advanced) Additional ATM parameter
 ga                       3.0         Charge scaling limiting value
 gc                       2.0         Charge scaling steepness
 wf                       6.0         Coordination number weighting
======================== =========== ============================================

Disabling the three-body dispersion (s9=0.0 or "3b": "none") changes the internal 
selection rules for damping parameters of a given method and prefers special two-body 
only damping parameters if available!

.. note::

    input_data.model.method with a full method name and input_data.keywords["params_tweaks"]
    cannot be provided at the same time. It is an error to provide both options at the
    same time.

Example
-------

>>> from dftd4.qcschema import run_qcschema
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
...         "method": "TPSS-D4",
...     },
...     keywords = {},
... )
>>> atomic_result = run_qcschema(atomic_input)
>>> atomic_result.return_result
-0.0002667885779142513
"""

from typing import Union

import numpy as np
import qcelemental as qcel

from .interface import DampingFunction, DampingParam, DispersionModel
from .library import get_api_version

_supported_drivers = [
    "energy",
    "gradient",
]

_available_levels = [
    "d4",
    "d4s",
]

_clean_dashlevel = str.maketrans("", "", "()")


def run_qcschema(
    input_data: Union[dict, qcel.models.AtomicInput]
) -> qcel.models.AtomicResult:
    """Perform disperson correction based on an atomic inputmodel"""

    if not isinstance(input_data, qcel.models.AtomicInput):
        atomic_input = qcel.models.AtomicInput(**input_data)
    else:
        atomic_input = input_data
    ret_data = atomic_input.dict()

    provenance = {
        "creator": "dftd4",
        "version": get_api_version(),
        "routine": "dftd4.qcschema.run_qcschema",
    }
    success = False
    return_result = 0.0
    properties = {}

    # Since it is a level hint we are forgiving if it is not present,
    # we are much less forgiving if the wrong level is hinted here.
    _level = atomic_input.keywords.get("level_hint", "d4")
    if _level.lower() not in _available_levels:
        ret_data.update(
            provenance=provenance,
            success=success,
            properties=properties,
            return_result=return_result,
            error=qcel.models.ComputeError(
                error_type="input error",
                error_message="Level '{}' is invalid for this dispersion correction".format(
                    _level
                ),
            ),
        )
        return qcel.models.AtomicResult(**ret_data)

    # Check if the method is provided and strip the “dashlevel” from the method
    _method = atomic_input.model.method.split("-")
    if _method[-1].lower().translate(_clean_dashlevel) == _level.lower():
        _method.pop()
    _method = "-".join(_method)
    if len(_method) == 0:
        _method = None

    # Obtain the damping function input
    _input_damp = atomic_input.keywords.get("damping_hint", {})
    _d2 = _input_damp.get("2b")
    _d3 = _input_damp.get("3b")

    # Obtain the input parameters for the damping function
    _input_param = atomic_input.keywords.get("params_tweaks")
    _has_tweaks = bool(_input_param)
    _input_param = _input_param if _has_tweaks else {}

    # Enforce that method and params_tweaks are mutually exclusive
    if _method is not None and _input_param:
        ret_data.update(
            provenance=provenance,
            success=False,
            properties=properties,
            return_result=return_result,
            error=qcel.models.ComputeError(
                error_type="input error",
                error_message=(
                    "input_data.model.method with a full method name and "
                    "input_data.keywords['params_tweaks'] cannot be provided at the same time."
                ),
            ),
        )
        return qcel.models.AtomicResult(**ret_data)

    # Obtain the parameters for the dispersion model
    if(_level.lower() == "d4s"):
        _model_param = {
            key: _input_param.pop(key, default)
            for key, default in (
                ("ga", 3.0),
                ("gc", 2.0),
            )
        }
    else: 
        _model_param = {
            key: _input_param.pop(key, default)
            for key, default in (
                ("ga", 3.0),
                ("gc", 2.0),
                ("wf", 6.0),
            )
        }

    try:
        
        # Damping function
        if _d2 is not None: 
            if _d3 is None:
                # Explicitly specified damping function without ATM
                damp = DampingFunction(damping_2b=_d2)
            else:
                # Explicitly specified damping function with ATM
                damp = DampingFunction(damping_2b=_d2, damping_3b=_d3)
        else: 
            if _d3 is None: 
                # Use default damping function with unmodified ATM
                damp = DampingFunction(model=_level)
            else:
                # Use default damping function with modified ATM
                damp = DampingFunction(model=_level, damping_3b=_d3)

        # Damping parameter
        if _has_tweaks:
            if all(key in _input_param for key in ("s6", "s8", "s9",
                                                   "a1", "a2", "a3", "a4",
                                                   "rs6", "rs8", "rs9",
                                                   "alp", "bet")):
                # Explicitly specify all possible parameters
                param = DampingParam(**_input_param)
            else:
                # Explicitly specified parameters for model default damping
                param = DampingParam(model=_level, **_input_param)
        else:
            if _method is None:
                raise TypeError("Method name or complete damping parameter set required")
            # Default parameters for a given method and model            
            param_kwargs = {"method": _method, "model": _level}
            if _d2 is not None:
                param_kwargs["damping_2b"] = _d2
            if _d3 is not None:
                param_kwargs["damping_3b"] = _d3
            param = DampingParam(**param_kwargs)

        # Dispersion Model
        disp = DispersionModel(
            atomic_input.molecule.atomic_numbers[atomic_input.molecule.real],
            atomic_input.molecule.geometry[atomic_input.molecule.real],
            atomic_input.molecule.molecular_charge,
            model = _level,
            **_model_param,
        )

        res = disp.get_dispersion(
            damp=damp,
            param=param,
            grad=atomic_input.driver == "gradient",
        )
        if atomic_input.keywords.get("property", False):
            res.update(**disp.get_properties())
        extras = {"dftd4": res}

        if atomic_input.driver == "gradient":
            if all(atomic_input.molecule.real):
                fullgrad = res.get("gradient")
            else:
                ireal = np.argwhere(atomic_input.molecule.real).reshape((-1))
                fullgrad = np.zeros_like(atomic_input.molecule.geometry)
                fullgrad[ireal, :] = res.get("gradient")

        properties.update(return_energy=res.get("energy"))

        if atomic_input.keywords.get("pair_resolved", False):
            res = disp.get_pairwise_dispersion(damp=damp, param=param)
            extras["dftd4"].update(res)

        success = atomic_input.driver in _supported_drivers
        if atomic_input.driver == "energy":
            return_result = properties["return_energy"]
        elif atomic_input.driver == "gradient":
            return_result = fullgrad
        else:
            ret_data.update(
                error=qcel.models.ComputeError(
                    error_type="input error",
                    error_message="Calculation succeeded but invalid driver request provided",
                ),
            )

        ret_data["extras"].update(extras)

    except (RuntimeError, TypeError, ValueError) as e:
        ret_data.update(
            error=qcel.models.ComputeError(
                error_type="input error", error_message=str(e)
            ),
        )

    ret_data.update(
        provenance=provenance,
        success=success,
        properties=properties,
        return_result=return_result,
    )

    return qcel.models.AtomicResult(**ret_data)
