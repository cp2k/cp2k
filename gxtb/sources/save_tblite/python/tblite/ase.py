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
The Python API of *tblite* natively support integration with the atomic
simulation environment (`ASE`_).
By constructing a calculator most functionality of ASE is readily available.
For details on building the Python API checkout the
:ref:`installation guide <python-build>`.

.. _ase: https://wiki.fysik.dtu.dk/ase/

Example
-------

The calculator can be used like other ASE calculators by attaching it to an atoms object.

>>> from ase.build import molecule
>>> from tblite.ase import TBLite
>>> atoms = molecule("H2O")
>>> atoms.calc = TBLite(method="GFN2-xTB")
>>> atoms.get_potential_energy()  # in eV
-137.96777594361677

The calculator will perform single point calculations as necessary, if properties for
the same atoms are requested the properties will be returned without additional single point
calculations.

>>> atoms.get_dipole_moment()  # in Debye
array([-7.83154860e-17,  1.10157136e-16, -4.76020974e-01])
>>> atoms.get_forces()  # in eV/A
array([[-3.31859412e-15,  7.83871951e-15, -7.49527168e-01],
       [ 3.08351626e-15, -1.53506496e-01,  3.74763584e-01],
       [ 2.35077857e-16,  1.53506496e-01,  3.74763584e-01]])

Changed properties of the atoms object will lead to an reevalution of the single point,
for example updating the initial charges of the atoms to set a total charge.

>>> atoms.set_initial_charges([0, 0, -1])
>>> atoms.get_potential_energy()  # in eV
-131.50015843891816

Some properties, like the total charge, can also be set via the parameters of the
calculator.

>>> atoms = molecule("H2O")
>>> atoms.calc = TBLite(method="GFN2-xTB", charge=-1)
>>> atoms.get_potential_energy()  # in eV
-131.50015843891822

Finally, the parameters of the calculator can be updated at any time using the *set*
method, which resets the calculation as necessary.

>>> atoms.calc.set(method="GFN1-xTB", charge=0)
>>> atoms.get_potential_energy()  # in eV
-156.96750590276173
>>> atoms.get_dipole_moment()  # in Debye
array([ 3.71500789e-16,  3.08439980e-16, -5.97048616e-01])
>>> atoms.get_forces()  # in eV/A
array([[ 2.88997963e-15,  5.70620869e-15, -7.84252969e-01],
       [-4.90812472e-15, -2.28835628e-01,  3.92126485e-01],
       [ 2.01814509e-15,  2.28835628e-01,  3.92126485e-01]])

The calculator also support atom partitioned properties, including energies and charges.

>>> atoms.get_charges()  # in e
array([-0.66558478,  0.33279239,  0.33279239])
>>> atoms.get_potential_energies()  # in eV
array([-132.37935247,  -12.29407672,  -12.29407672])
"""

try:
    import ase
    import ase.calculators.calculator
    from ase.atoms import Atoms
    from ase.units import Bohr, Hartree, kB
except ModuleNotFoundError as e:
    raise ModuleNotFoundError("This submodule requires ASE installed") from e

from typing import Any, Dict, List, Optional

from .interface import Calculator


class TBLite(ase.calculators.calculator.Calculator):
    r"""
    ASE calculator for using xTB Hamiltonians from the tblite library.
    Supported properties by this calculator are:

    - *energy* (*free_energy*): The total energy of the system (in eV)
    - *energies*: The atom-partitioned energy of the system (in eV)
    - *forces*: The derivative of the energy with respect to the atomic positions (in eV/Ang)
    - *stress*: The derivative of the energy with respect to cell deformations (PBC only, in eV/Ang³)
    - *dipole*: The dipole moment of the system (in Debye)
    - *charges*: Mulliken charges of the system (in e)

    Supported keywords are

    ======================== ================= ====================================================
     Keyword                  Default           Description
    ======================== ================= ====================================================
     method                   "GFN2-xTB"        Underlying method for energy and forces
     charge                   None              Total charge of the system
     multiplicity             None              Total multiplicity of the system
     accuracy                 1.0               Numerical accuracy of the calculation
     electronic_temperature   300.0             Electronic temperatur in Kelvin
     max_iterations           250               Iterations for self-consistent evaluation
     initial_guess            "sad"             Initial guess for wavefunction (sad or eeq)
     mixer_damping            0.4               Damping parameter for self-consistent mixer
     electric_field           None              Uniform electric field vector (in V/A)
     spin_polarization        None              Spin polarization (scaling factor)
     solvation                None              Solvation model to use (see below for details)
     cache_api                True              Reuse generate API objects (recommended)
     verbosity                1                 Set verbosity of printout
    ======================== ================= ====================================================

    Solvation models are supported by passing a tuple to the *solvation* keyword.
    The first element of the tuple is the name of the solvation model, the following
    arguments are passed to initialize the solvation model.

    ================== ====================================================
     Solvation model    Arguments
    ================== ====================================================
     "alpb"             Solvent name (str), solution state (optional)
     "gbsa"             Solvent name (str), solution state (optional)
     "cpcm"             Epsilon (float)
     "gbe"              Epsilon (float), Born kernel
     "gb"               Epsilon (float), Born kernel
    ================== ====================================================

    Example
    -------

    An ASE calculator can be constructed by using the *TBLite* class provided
    by the *tblite.ase* module.
    For example to perform a single point calculation for a CO\ :sub:`2`
    crystal use

    >>> from tblite.ase import TBLite
    >>> from ase.atoms import Atoms
    >>> import numpy as np
    >>> atoms = Atoms(
    ...     symbols="C4O8",
    ...     positions=np.array(
    ...         [
    ...             [0.9441259872, 0.9437851680, 0.9543505632],
    ...             [3.7179966528, 0.9556570368, 3.7316862240],
    ...             [3.7159517376, 3.7149292800, 0.9692330016],
    ...             [0.9529872864, 3.7220864832, 3.7296981120],
    ...             [1.6213905408, 1.6190616096, 1.6313879040],
    ...             [0.2656685664, 0.2694175776, 0.2776540416],
    ...             [4.3914553920, 1.6346256864, 3.0545920000],
    ...             [3.0440834880, 0.2764611744, 4.4080419264],
    ...             [4.3910577696, 3.0416409504, 0.2881058304],
    ...             [3.0399936576, 4.3879335936, 1.6497353376],
    ...             [0.2741322432, 4.4003734944, 3.0573754368],
    ...             [1.6312174944, 3.0434586528, 4.4023048032],
    ...         ]
    ...     ),
    ...     cell=np.array([5.68032, 5.68032, 5.68032]),
    ...     pbc=np.array([True, True, True]),
    ... )
    >>> atoms.calc = TBLite(method="GFN1-xTB", verbosity=0)
    >>> atoms.get_potential_energy()  # result in eV
    -1257.0943962462964

    The resulting calculator can be used like most ASE calculator, *e.g.* for optimizing geometries.
    """

    implemented_properties = [
        "energy",
        "energies",
        "forces",
        "charges",
        "dipole",
        "stress",
    ]

    default_parameters = {
        "method": "GFN2-xTB",
        "charge": None,
        "multiplicity": None,
        "accuracy": 1.0,
        "guess": "sad",
        "max_iterations": 250,
        "mixer_damping": 0.4,
        "electric_field": None,
        "spin_polarization": None,
        "solvation": None,
        "electronic_temperature": 300.0,
        "cache_api": True,
        "verbosity": 1,
    }

    _res = None
    _xtb = None

    def __init__(
        self,
        atoms: Optional[Atoms] = None,
        **kwargs,
    ):
        """
        Construct the TBLite base calculator object.
        """

        ase.calculators.calculator.Calculator.__init__(self, atoms=atoms, **kwargs)

    def set(self, **kwargs) -> dict:
        """
        Set new parameters to TBLite. Will automatically reconstruct the underlying
        model in case critical parameters change.

        Example
        -------
        >>> from ase.build import molecule
        >>> from tblite.ase import TBLite
        >>> atoms = molecule("H2O")
        >>> atoms.calc = TBLite(method="GFN2-xTB")
        >>> atoms.get_potential_energy()
        -137.96777625229421
        >>> atoms.calc.set(method="GFN1-xTB")
        {'method': 'GFN1-xTB'}
        >>> atoms.get_potential_energy()
        -156.9675057724589
        """

        _update_parameters(kwargs)

        changed_parameters = ase.calculators.calculator.Calculator.set(self, **kwargs)

        # Always reset the calculation if parameters change
        if changed_parameters:
            self.reset()

        # If the method is changed, invalidate the cached calculator as well
        if (
            "method" in changed_parameters
            or "electric_field" in changed_parameters
            or "spin_polarization" in changed_parameters
            or "solvation" in changed_parameters
        ):
            self._xtb = None
            self._res = None

        # Minor changes can be updated in the API calculator directly
        if self._xtb is not None:
            if "accuracy" in changed_parameters:
                self._xtb.set("accuracy", self.parameters.accuracy)

            if "electronic_temperature" in changed_parameters:
                self._xtb.set(
                    "temperature",
                    self.parameters.electronic_temperature * kB / Hartree,
                )

            if "max_iterations" in changed_parameters:
                self._xtb.set("max-iter", self.parameters.max_iterations)

            if "initial_guess" in changed_parameters:
                self._xtb.set("guess", {"sad": 0, "eeq": 1}[self.parameters.guess])

            if "mixer_damping" in changed_parameters:
                self._xtb.set("mixer-damping", self.parameters.mixer_damping)

            if (
                "charge" in changed_parameters or "multiplicity" in changed_parameters
            ) and self.atoms is not None:
                self._xtb.update(
                    charge=_get_charge(self.atoms, self.parameters),
                    uhf=_get_uhf(self.atoms, self.parameters),
                )

        return changed_parameters

    def reset(self) -> None:
        """
        Clear all information from old calculation. This will only remove the cached
        API objects in case the `cache_api` is set to False.
        """
        ase.calculators.calculator.Calculator.reset(self)

        if not self.parameters.cache_api:
            self._xtb = None
            self._res = None

    def _check_api_calculator(self, system_changes: List[str]) -> None:
        """Check state of API calculator and reset if necessary"""

        # Changes in positions, cell, charges and magnetic moments can use a normal update
        _reset = system_changes.copy()
        for _change in system_changes:
            if _change in ("positions", "cell", "initial_charges", "initial_magmoms"):
                _reset.remove(_change)

        # Invalidate cached calculator and results object
        if _reset:
            self._xtb = None
            self._res = None
        else:
            if system_changes and self._xtb is not None:
                try:
                    _cell = self.atoms.cell
                    self._xtb.update(
                        positions=self.atoms.positions / Bohr,
                        lattice=_cell / Bohr,
                        charge=_get_charge(self.atoms, self.parameters),
                        uhf=_get_uhf(self.atoms, self.parameters),
                    )
                # An exception in this part means the geometry is bad,
                # still we will give a complete reset a try as well
                except RuntimeError:
                    self._xtb = None
                    self._res = None

    def calculate(
        self,
        atoms: Optional[Atoms] = None,
        properties: Optional[List[str]] = None,
        system_changes: List[str] = ase.calculators.calculator.all_changes,
    ) -> None:
        """
        Perform actual calculation by calling the TBLite API

        Example
        -------

        >>> from ase.build import molecule
        >>> from tblite.ase import TBLite
        >>> calc = TBLite(method="GFN2-xTB")
        >>> calc.calculate(molecule("H2O"))
        >>> calc.get_potential_energy()
        -137.96777625229421
        >>> calc.calculate(molecule("CH4"))
        >>> calc.get_potential_energy()
        -113.60956621093894

        Raises
        ------
        ase.calculators.calculator.InputError
            on invalid input passed to the interface module

        ase.calculators.calculator.CalculationFailed
            in case of an `RuntimeError` in the library
        """

        if not properties:
            properties = ["energy"]
        ase.calculators.calculator.Calculator.calculate(self, atoms, properties, system_changes)

        self._check_api_calculator(system_changes)

        if self._xtb is None:
            self._xtb = _create_api_calculator(self.atoms, self.parameters)

        try:
            self._res = self._xtb.singlepoint(self._res)
        except RuntimeError as e:
            raise ase.calculators.calculator.CalculationFailed(str(e)) from e

        # These properties are garanteed to exist for all implemented calculators
        self.results["energy"] = self._res["energy"] * Hartree
        self.results["energies"] = self._res["energies"] * Hartree
        self.results["free_energy"] = self.results["energy"]
        self.results["forces"] = -self._res["gradient"] * Hartree / Bohr
        self.results["charges"] = self._res["charges"]
        self.results["dipole"] = self._res["dipole"] * Bohr
        self.results["bond-orders"] = self._res["bond-orders"]
        # stress tensor is only returned for periodic systems
        if self.atoms.pbc.any():
            _stress = self._res["virial"] * Hartree / self.atoms.get_volume()
            self.results["stress"] = _stress.flat[[0, 4, 8, 5, 2, 1]]


def _create_api_calculator(
    atoms: ase.Atoms,
    parameters: ase.calculators.calculator.Parameters,
) -> Calculator:
    """Create a new API calculator object"""

    try:
        _cell = atoms.cell
        _periodic = atoms.pbc
        _charge = _get_charge(atoms, parameters)
        _uhf = _get_uhf(atoms, parameters)

        calc = Calculator(
            parameters.method,
            atoms.numbers,
            atoms.positions / Bohr,
            _charge,
            _uhf,
            _cell / Bohr,
            _periodic,
        )
        calc.set("accuracy", parameters.accuracy)
        calc.set(
            "temperature",
            parameters.electronic_temperature * kB / Hartree,
        )
        calc.set("max-iter", parameters.max_iterations)
        calc.set("guess", {"sad": 0, "eeq": 1}[parameters.guess])
        calc.set("mixer-damping", parameters.mixer_damping)
        calc.set("verbosity", parameters.verbosity)
        if parameters.electric_field is not None:
            calc.add(
                "electric-field",
                parameters.electric_field * Bohr / Hartree,
            )
        if parameters.spin_polarization is not None:
            calc.add("spin-polarization", parameters.spin_polarization)
        if parameters.solvation is not None:
            solvation_model, *solvation_args = parameters.solvation
            if isinstance(solvation_args, (tuple, list)):
                calc.add(f"{solvation_model}-solvation", *solvation_args)
            else:
                calc.add(f"{solvation_model}-solvation", solvation_args)

    except RuntimeError as e:
        raise ase.calculators.calculator.InputError(str(e)) from e

    return calc


def _get_charge(atoms: ase.Atoms, parameters: ase.calculators.calculator.Parameters) -> int:
    """
    Get the total charge of the system.
    If no charge is provided, the total charge of the system is calculated
    by summing the initial charges of all atoms.
    """
    return atoms.get_initial_charges().sum() if parameters.charge is None else parameters.charge


def _get_uhf(atoms: ase.Atoms, parameters: ase.calculators.calculator.Parameters) -> int:
    """
    Get the number of unpaired electrons.
    If no multiplicity is provided, the number of unpaired electrons
    is calculated by summing the initial magnetic moments of all atoms.
    """
    return (
        int(atoms.get_initial_magnetic_moments().sum().round())
        if parameters.multiplicity is None
        else parameters.multiplicity - 1
    )


def _update_parameters(kwargs: Dict[str, Any]) -> None:
    """
    Update the parameters of the calculator with the provided keyword arguments.
    This function is used to update the parameters of the calculator
    when new values are provided.
    """

    for key in (
        "alpb",
        "gbsa",
        "cpcm",
        "gbe",
        "gb",
    ):
        if f"{key}_solvation" in kwargs:
            value = kwargs.pop(f"{key}_solvation")
            if isinstance(value, (tuple, list)):
                value = (key, *value)
            else:
                value = (key, value)

            if "solvation" in kwargs:
                raise ase.calculators.calculator.InputError(
                    "Multiple solvation models provided, can only use one"
                )
            kwargs["solvation"] = value


if "tblite" not in ase.calculators.calculator.external_calculators:
    ase.calculators.calculator.register_calculator_class("tblite", TBLite)
