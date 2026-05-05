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
ASE Support
-----------

`ASE calculator <https://wiki.fysik.dtu.dk/ase/>`_ implementation
for the ``dftd4`` program.

This module provides a basic single point calculator implementations
to integrate the ``dftd4`` API into existing ASE workflows.
To use DFTD4 as dispersion correction the ``ase.calculators.mixing``
module can be used to combine DFTD4 with a DFT calculator using
the ``SumCalculator``.

Supported properties by this calculator are:

- energy (free_energy)
- forces
- stress

Supported keywords are

======================== ============ ============================================
 Keyword                  Default      Description
======================== ============ ============================================
 method                   None         Method to calculate dispersion for
 params_tweaks            None         Optional dict with the damping parameters
 damping_hint             None         Optional dict with the damping functions
 cache_api                True         Reuse generate API objects (recommended)
 model                    d4           Used dispersion Model (D4S or D4 (default))
======================== ============ ============================================

The optional damping_hint dict contains the two-body ("2b") and three-body ("3b") 
damping function types. If not provided the default damping functions of the model
will be used. If the three-body damping function is set to "none" the ATM contribution
will be disabled.

The optional params_tweaks dict contains the damping parameters (either all
damping parameters or parameters for the model specific default (e.g., rational 
+ zero-avg for d4 requiring at least s8, a1 and a2). The parameters are:

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
======================== =========== ============================================

Disabling the three-body dispersion (s9=0.0 or "3b": "none") changes the internal 
selection rules for damping parameters of a given method and prefers special two-body 
only damping parameters if available!

Example
-------
>>> from ase.build import molecule
>>> from dftd4.ase import DFTD4
>>> atoms = molecule('H2O')
>>> atoms.calc = DFTD4(method="TPSS")
>>> atoms.get_potential_energy()
-0.007310393443152083
>>> atoms.calc.set(method="PBE")
{'method': 'PBE'}
>>> atoms.get_potential_energy()
-0.005358475432239303
>>> atoms.get_forces()
array([[-0.        , -0.        ,  0.00296845],
       [-0.        ,  0.00119152, -0.00148423],
       [-0.        , -0.00119152, -0.00148423]])
"""

try:
    import ase
except ModuleNotFoundError:
    raise ModuleNotFoundError("This submodule requires ASE installed")

from typing import List, Optional

from ase.atoms import Atoms
from ase.calculators.calculator import (
    CalculationFailed,
    Calculator,
    InputError,
    all_changes,
)
from ase.calculators.mixing import SumCalculator
from ase.units import Bohr, Hartree

from .interface import DampingFunction, DampingParam, DispersionModel


class DFTD4(Calculator):
    """
    ASE calculator for DFT-D4 related methods.
    The DFTD4 class can access all methods exposed by the ``dftd4`` API.

    Example
    -------
    >>> from ase.build import molecule
    >>> from ase.calculators.mixing import SumCalculator
    >>> from ase.calculators.nwchem import NWChem
    >>> from dftd4.ase import DFTD4
    >>> atoms = molecule('H2O')
    >>> atoms.calc = SumCalculator([DFTD4(method="PBE", model="d4"), NWChem(xc="PBE")])
    """

    implemented_properties = [
        "energy",
        "forces",
        "stress",
    ]

    default_parameters = {
        "method": None,
        "params_tweaks": None,
        "damping_hint": None,
        "cache_api": True,
        "model" : "d4"
    }

    _disp = None

    def __init__(
        self,
        atoms: Optional[Atoms] = None,
        **kwargs,
    ):
        """Construct the dftd4 dispersion model object."""

        Calculator.__init__(self, atoms=atoms, **kwargs)

    def add_calculator(self, other: Calculator) -> Calculator:
        """
        Convenience function to allow DFTD4 to combine itself with another calculator
        by returning a SumCalculator:

        Example
        -------
        >>> from ase.build import molecule
        >>> from ase.calculators.emt import EMT
        >>> from dftd4.ase import DFTD4
        >>> atoms = molecule("C60")
        >>> atoms.calc = DFTD4(method="pbe", model="d4").add_calculator(EMT())
        >>> atoms.get_potential_energy()
        6.348142387048062
        >>> [calc.get_potential_energy() for calc in atoms.calc.calcs]
        [-6.015477436263984, 12.363619823312046]
        """
        return SumCalculator([self, other])

    def set(self, **kwargs) -> dict:
        """Set new parameters to dftd4"""

        changed_parameters = Calculator.set(self, **kwargs)

        # Always reset the calculation if parameters change
        if changed_parameters:
            self.reset()

        return changed_parameters

    def reset(self) -> None:
        """Clear all information from old calculation"""
        Calculator.reset(self)

        if not self.parameters.cache_api:
            self._disp = None

    def _check_api_calculator(self, system_changes: List[str]) -> None:
        """Check state of API calculator and reset if necessary"""

        # Changes in positions and cell parameters can use a normal update
        _reset = system_changes.copy()
        if "positions" in _reset:
            _reset.remove("positions")
        if "cell" in _reset:
            _reset.remove("cell")

        # Invalidate cached calculator and results object
        if _reset:
            self._disp = None
        else:
            if system_changes and self._disp is not None:
                try:
                    _cell = self.atoms.cell
                    self._disp.update(
                        self.atoms.positions / Bohr,
                        _cell / Bohr,
                    )
                # An exception in this part means the geometry is bad,
                # still we will give a complete reset a try as well
                except RuntimeError:
                    self._disp = None

    def _create_api_calculator(self) -> DispersionModel:
        """Create a new API calculator object"""

        try:
            _cell = self.atoms.cell
            _periodic = self.atoms.pbc
            _charge = self.atoms.get_initial_charges().sum()

            disp = DispersionModel(
                self.atoms.numbers,
                self.atoms.positions / Bohr,
                _charge,
                _cell / Bohr,
                _periodic,
                model=self.parameters.get("model"),
            )

        except (RuntimeError, ValueError) as e:
            raise InputError("Cannot construct dispersion model for dftd4") from e

        return disp

    def _create_damping_function(self) -> DampingFunction:
        """Create a new DampingFunction object."""

        _model = (self.parameters.get("model") or "d4").lower().replace(" ", "")

        _damping = self.parameters.get("damping_hint", None) or {}
        _d2 = _damping.get("2b")
        _d3 = _damping.get("3b")

        try:
            # Damping function
            if _d2 is not None:
                if _d3 is None:
                   # Explicitly specified damping function without ATM
                   return DampingFunction(damping_2b=_d2)
                else:
                   # Explicitly specified damping function with ATM
                   return DampingFunction(damping_2b=_d2, damping_3b=_d3)
            else: 
                if _d3 is None: 
                   # Use default damping function with unmodified ATM
                   return DampingFunction(model=_model)
                else:
                   # Use default damping function with modified ATM
                   return DampingFunction(model=_model, damping_3b=_d3)

        except (RuntimeError, ValueError, TypeError) as e:
            raise InputError("Cannot construct damping function for dftd4") from e

    def _create_damping_param(self) -> DampingParam:
        """Create a new DampingParam object."""

        _model = (self.parameters.get("model") or "d4").lower().replace(" ", "")

        _method = self.parameters.get("method", None)
        if _method is not None and str(_method).strip() == "":
            _method = None

        _tweaks = self.parameters.get("params_tweaks", None)
        _has_tweaks = bool(_tweaks)
        _tweaks = _tweaks if _has_tweaks else {}

        _damping = self.parameters.get("damping_hint", None) or {}
        _d2 = _damping.get("2b")
        _d3 = _damping.get("3b")

        try:
            if _has_tweaks:
                if all(key in _tweaks for key in ("s6", "s8", "s9",
                                                  "a1", "a2", "a3", "a4",
                                                  "rs6", "rs8", "rs9",
                                                  "alp", "bet")):
                    # Explicitly specify all possible parameters
                    return DampingParam(**_tweaks)
                else:
                    # Explicitly specified parameters for model default damping
                    return DampingParam(model=_model, **_tweaks)
            else:
                if _method is None:
                   raise TypeError("Method name or complete damping parameter set required")
                # Default parameters for a given method and model
                param_kwargs = {"method": _method, "model": _model}
                if _d2 is not None:
                   param_kwargs["damping_2b"] = _d2
                if _d3 is not None:
                   param_kwargs["damping_3b"] = _d3
                return DampingParam(**param_kwargs)

        except (RuntimeError, ValueError, TypeError) as e:
            raise InputError("Cannot construct damping parameter for dftd4") from e

    def calculate(
        self,
        atoms: Optional[Atoms] = None,
        properties: List[str] = None,
        system_changes: List[str] = all_changes,
    ) -> None:
        """Perform actual calculation with by calling the dftd4 API"""

        if not properties:
            properties = ["energy"]
        Calculator.calculate(self, atoms, properties, system_changes)

        self._check_api_calculator(system_changes)

        if self._disp is None:
            self._disp = self._create_api_calculator()

        _damp = self._create_damping_function()
        _dpar = self._create_damping_param()

        try:
            _res = self._disp.get_dispersion(damp=_damp, param=_dpar, grad=True)
        except RuntimeError:
            raise CalculationFailed("dftd4 could not evaluate input")

        # These properties are garanteed to exist for all implemented calculators
        self.results["energy"] = _res.get("energy") * Hartree
        self.results["free_energy"] = self.results["energy"]
        self.results["forces"] = -_res.get("gradient") * Hartree / Bohr
        # stress tensor is only returned for periodic systems
        if self.atoms.pbc.any():
            _stress = _res.get("virial") * Hartree / self.atoms.get_volume()
            self.results["stress"] = _stress.flat[[0, 4, 8, 5, 2, 1]]
