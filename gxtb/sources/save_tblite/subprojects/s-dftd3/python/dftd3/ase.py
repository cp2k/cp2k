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
ASE Support
===========

`ASE calculator <https://wiki.fysik.dtu.dk/ase/>`_ implementation
for the ``s-dftd3`` program.

This module provides a basic single point calculator implementations
to integrate the ``s-dftd3`` API into existing ASE workflows.
To use DFTD3 as dispersion correction the ``ase.calculators.mixing``
module can be used to combine DFTD3 with a DFT calculator using
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
 damping                  None         Damping function to use
 params_tweaks            {}           Optional dict with the damping parameters
 realspace_cutoff         {}           Optional dict to override cutoff values
 cache_api                True         Reuse generate API objects (recommended)
======================== ============ ============================================

The params_tweaks dict contains the damping parameters, at least s8, a1 and a2
must be provided

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

Either method or s8, a1 and a2 must be provided, s9 can be used to overwrite
the ATM scaling if the method is provided in the model.
Disabling the three-body dispersion (s9=0.0) changes the internal selection rules
for damping parameters of a given method and prefers special two-body only
damping parameters if available!

The realspace cutoff parameters allow adjusting the distance values for which
interactions are considered

================== =========== ==========================================
 Realspace cutoff   Default     Description
================== =========== ==========================================
 disp2              60 * Bohr   Pairwise dispersion interactions
 disp3              40 * Bohr   Triple dispersion interactions
 cn                 40 * Bohr   Coordination number counting
================== =========== ==========================================

Values provided in the dict are expected to be in Angstrom. When providing values
in Bohr multiply the inputs by the `ase.units.Bohr` constant.

Example
-------
>>> from ase.build import molecule
>>> from dftd3.ase import DFTD3
>>> atoms = molecule('H2O')
>>> atoms.calc = DFTD3(method="TPSS", damping="d3bj")
>>> atoms.get_potential_energy()
-0.011441640787279111
>>> atoms.calc.set(method="PBE")
{'method': 'PBE'}
>>> atoms.get_potential_energy()
-0.009781920198843976
>>> atoms.get_forces()
array([[-0.        , -0.        ,  0.00009572],
       [-0.        , -0.00004060, -0.00004786],
       [-0.        ,  0.00004060, -0.00004786]])
"""

try:
    import ase
except ModuleNotFoundError:
    raise ModuleNotFoundError("This submodule requires ASE installed")

from typing import List, Optional

from .interface import (
    DispersionModel,
    DampingParam,
    RationalDampingParam,
    ZeroDampingParam,
    ModifiedRationalDampingParam,
    ModifiedZeroDampingParam,
    OptimizedPowerDampingParam,
    CSODampingParam,
)
from ase.calculators.calculator import (
    Calculator,
    InputError,
    CalculationFailed,
    all_changes,
)
from ase.calculators.mixing import SumCalculator
from ase.atoms import Atoms
from ase.units import Hartree, Bohr


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

_inv_bohr = 1.0 / Bohr
_hartree_per_bohr = Hartree / Bohr


class DFTD3(Calculator):
    """
    ASE calculator for DFT-D3 related methods.
    The DFTD3 class can access all methods exposed by the ``dftd3`` API.

    Example
    -------
    >>> from ase.build import molecule
    >>> from ase.calculators.mixing import SumCalculator
    >>> from ase.calculators.nwchem import NWChem
    >>> from dftd3.ase import DFTD3
    >>> atoms = molecule('H2O')
    >>> atoms.calc = SumCalculator([DFTD3(method="PBE", damping="d3bj"), NWChem(xc="PBE")])
    """

    implemented_properties = [
        "energy",
        "forces",
        "stress",
    ]

    default_parameters = {
        "method": None,
        "damping": None,
        "params_tweaks": {},
        "realspace_cutoff": {},
        "cache_api": True,
    }

    _disp = None
    _dpar = None

    def __init__(
        self,
        atoms: Optional[Atoms] = None,
        **kwargs,
    ):
        """Construct the DFTD3 dispersion model object."""

        Calculator.__init__(self, atoms=atoms, **kwargs)

    def add_calculator(self, other: Calculator) -> Calculator:
        """
        Convenience function to allow DFTD3 to combine itself with another calculator
        by returning a SumCalculator:

        Example
        -------
        >>> from ase.build import molecule
        >>> from ase.calculators.emt import EMT
        >>> from dftd3.ase import DFTD3
        >>> atoms = molecule("C60")
        >>> atoms.calc = DFTD3(method="pbe", damping="d3bj").add_calculator(EMT())
        >>> atoms.get_potential_energy()
        7.389424410228332
        >>> [calc.get_potential_energy() for calc in atoms.calc.calcs]
        [-4.974195589771621, 12.363619999999953]
        """
        return SumCalculator([self, other])

    def set(self, **kwargs) -> dict:
        """Set new parameters to dftd3"""

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
        self._dpar = None

    def _check_api_calculator(self, system_changes: List[str]) -> None:
        """Check state of API calculator and reset if necessary"""

        if not system_changes:
            return

        # Changes in positions and cell parameters can use a normal update
        _changes = set(system_changes)
        _changes.discard("positions")
        _changes.discard("cell")

        # Invalidate cached calculator and results object
        if _changes:
            self._disp = None
        elif self._disp is not None:
            try:
                self._disp.update(
                    self.atoms.positions * _inv_bohr,
                    self.atoms.cell[:] * _inv_bohr,
                )
            # An exception in this part means the geometry is bad,
            # still we will give a complete reset a try as well
            except RuntimeError:
                self._disp = None

    def _create_api_calculator(self) -> DispersionModel:
        """Create a new API calculator object"""

        try:
            disp = DispersionModel(
                self.atoms.numbers,
                self.atoms.positions * _inv_bohr,
                self.atoms.cell[:] * _inv_bohr,
                self.atoms.pbc,
            )

        except RuntimeError as e:
            raise InputError("Cannot construct dispersion model for dftd3") from e

        return disp

    def _apply_realspace_cutoff(self, disp: DispersionModel) -> None:
        """Apply realspace cutoff parameters to dispersion model"""

        try:
            if self.parameters.realspace_cutoff:
                disp2 = (
                    self.parameters.realspace_cutoff.get("disp2", 60.0 * Bohr) / Bohr
                )
                disp3 = (
                    self.parameters.realspace_cutoff.get("disp3", 40.0 * Bohr) / Bohr
                )
                cn = self.parameters.realspace_cutoff.get("cn", 40.0 * Bohr) / Bohr

                disp.set_realspace_cutoff(disp2=disp2, disp3=disp3, cn=cn)
        except RuntimeError as e:
            raise InputError("Cannot update realspace cutoff for dftd3") from e

    def _create_damping_param(self) -> DampingParam:
        """Create a new API damping parameter object"""

        try:
            params_tweaks = (
                self.parameters.params_tweaks
                if self.parameters.params_tweaks
                else {"method": self.parameters.get("method")}
            )
            dpar = _damping_param[self.parameters.get("damping")](**params_tweaks)

        except RuntimeError as e:
            raise InputError("Cannot construct damping parameter for dftd3") from e

        return dpar

    def calculate(
        self,
        atoms: Optional[Atoms] = None,
        properties: List[str] = None,
        system_changes: List[str] = all_changes,
    ) -> None:
        """Perform actual calculation with by calling the dftd3 API"""

        if not properties:
            properties = ["energy"]
        Calculator.calculate(self, atoms, properties, system_changes)

        self._check_api_calculator(system_changes)

        if self._disp is None:
            self._disp = self._create_api_calculator()

        # Apply realspace cutoff before evaluation (works with cached calculator)
        self._apply_realspace_cutoff(self._disp)

        if self._dpar is None:
            self._dpar = self._create_damping_param()

        _need_grad = "forces" in properties or "stress" in properties

        try:
            _res = self._disp.get_dispersion(param=self._dpar, grad=_need_grad)
        except RuntimeError:
            raise CalculationFailed("dftd3 could not evaluate input")

        # These properties are garanteed to exist for all implemented calculators
        self.results["energy"] = float(_res.get("energy")) * Hartree
        self.results["free_energy"] = self.results["energy"]
        if _need_grad:
            _gradient = _res.get("gradient")
            _gradient *= -_hartree_per_bohr
            self.results["forces"] = _gradient
            # stress tensor is only returned for periodic systems
            if self.atoms.pbc.any():
                _virial = _res.get("virial")
                _virial *= Hartree / self.atoms.get_volume()
                self.results["stress"] = _virial.flat[[0, 4, 8, 5, 2, 1]]
