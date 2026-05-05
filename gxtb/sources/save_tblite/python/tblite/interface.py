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
Definition of the basic interface to library for most further integration in
other Python frameworks. The classes defined here allow a more Pythonic usage
of the library in actual workflows than the low-level access provided in the
CFFI generated wrappers.
"""

from typing import Any, List, Optional, Union

import numpy as np
from . import library
from .exceptions import TBLiteRuntimeError, TBLiteValueError


class Structure:
    """
    .. Molecular structure data

    Represents a wrapped structure object in ``tblite``.
    The molecular structure data object has a fixed number of atoms
    and immutable atomic identifiers.

    Example
    -------
    >>> from tblite.interface import Structure
    >>> import numpy as np
    >>> mol = Structure(
    ...     positions=np.array([
    ...         [+0.00000000000000, +0.00000000000000, -0.73578586109551],
    ...         [+1.44183152868459, +0.00000000000000, +0.36789293054775],
    ...         [-1.44183152868459, +0.00000000000000, +0.36789293054775],
    ...     ]),
    ...     numbers = np.array([8, 1, 1]),
    ... )
    ...
    >>> len(mol)
    3

    Raises
    ------
    TBLiteValueError
        on invalid input, like incorrect shape / type of the passed arrays
    """

    _mol = library.ffi.NULL

    def __init__(
        self,
        numbers: Union[np.ndarray, List[int]],
        positions: np.ndarray,
        charge: Optional[float] = None,
        uhf: Optional[int] = None,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
    ):
        """
        Create new molecular structure data from arrays. The returned object has
        immutable atomic species and boundary condition, also the total number of
        atoms cannot be changed.

        Raises
        ------
        TBLiteValueError
            on invalid input, like incorrect shape / type of the passed arrays
        """
        if isinstance(numbers, list):
            numbers = np.asarray(numbers)

        if positions.size % 3 != 0:
            raise TBLiteValueError("Expected tripels of cartesian coordinates")

        if 3 * numbers.size != positions.size:
            raise TBLiteValueError(
                "Dimension missmatch between numbers and positions"
            )

        self._natoms = len(numbers)
        _numbers = np.ascontiguousarray(numbers, dtype="i4")
        _positions = np.ascontiguousarray(positions, dtype=float)

        _charge = _ref("double", charge)
        _uhf = _ref("int", uhf)

        if lattice is not None:
            if lattice.size != 9:
                raise TBLiteValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype="float")
        else:
            _lattice = None

        if periodic is not None:
            if periodic.size != 3:
                raise TBLiteValueError("Invalid periodicity provided")
            _periodic = np.ascontiguousarray(periodic, dtype="bool")
        else:
            _periodic = None

        self._mol = library.new_structure(
            self._natoms,
            _cast("int*", _numbers),
            _cast("double*", _positions),
            _charge,
            _uhf,
            _cast("double*", _lattice),
            _cast("bool*", _periodic),
        )

    def __len__(self):
        return self._natoms

    def update(
        self,
        positions: Optional[np.ndarray] = None,
        lattice: Optional[np.ndarray] = None,
        charge: Optional[float] = None,
        uhf: Optional[int] = None,
    ) -> None:
        """Update coordinates and lattice parameters, both provided in
        atomic units (Bohr).
        The lattice update is optional also for periodic structures.

        Generally, only the cartesian coordinates and the lattice parameters
        can be updated, every other modification, regarding total charge,
        total spin, boundary condition, atomic types or number of atoms
        requires the complete reconstruction of the object.

        Raises
        ------
        TBLiteValueError
            on invalid input, like incorrect shape / type of the passed arrays
        """
        if charge is not None:
            _charge = _ref("double", charge)
            library.update_structure_charge(self._mol, _charge)

        if uhf is not None:
            _uhf = _ref("int", uhf)
            library.update_structure_uhf(self._mol, _uhf)

        if positions is None:
            return

        if 3 * len(self) != positions.size:
            raise TBLiteValueError("Dimension missmatch for positions")
        _positions = np.ascontiguousarray(positions, dtype="float")

        if lattice is not None:
            if lattice.size != 9:
                raise TBLiteValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype="float")
        else:
            _lattice = None

        library.update_structure_geometry(
            self._mol,
            _cast("double*", _positions),
            _cast("double*", _lattice),
        )


class Result:
    """
    .. Calculation result and restart data

    Container for calculation results, can be passed to a single point calculation
    as restart data. Allows to retrieve individual results or export all results as dict.

    Example
    -------
    >>> from tblite.interface import Calculator
    >>> import numpy as np
    >>> calc = Calculator(
    ...     method="GFN2-xTB",
    ...     numbers=np.array([14, 1, 1, 1, 1]),
    ...     positions=np.array([
    ...         [ 0.000000000000, 0.000000000000, 0.000000000000],
    ...         [ 1.617683897558, 1.617683897558,-1.617683897558],
    ...         [-1.617683897558,-1.617683897558,-1.617683897558],
    ...         [ 1.617683897558,-1.617683897558, 1.617683897558],
    ...         [-1.617683897558, 1.617683897558, 1.617683897558],
    ...     ]),
    ... )
    >>> res = calc.singlepoint()
    ------------------------------------------------------------
      cycle        total energy    energy error   density error
    ------------------------------------------------------------
          1     -3.710873811182  -3.7424104E+00   1.5535536E-01
          2     -3.763115011046  -5.2241200E-02   4.8803250E-02
          3     -3.763205560205  -9.0549159E-05   2.0456319E-02
          4     -3.763213200846  -7.6406409E-06   2.7461272E-03
          5     -3.763250843634  -3.7642788E-05   2.4071582E-03
          6     -3.763236758520   1.4085114E-05   2.3635321E-04
          7     -3.763237287151  -5.2863152E-07   3.5812281E-04
          8     -3.763233829618   3.4575334E-06   8.5040576E-06
          9     -3.763233750684   7.8934355E-08   1.1095285E-07
    ------------------------------------------------------------
    >>> res.get("energy")
    -3.7632337506836944
    >>> res = calc.singlepoint(res)
    ------------------------------------------------------------
      cycle        total energy    energy error   density error
    ------------------------------------------------------------
          1     -3.763233752583  -3.7947704E+00   1.3823013E-07
          2     -3.763233751557   1.0256169E-09   1.4215249E-08
    ------------------------------------------------------------

    Raises
    ------
    TBliteValueError
        on invalid input, like incorrect shape / type of the passed arrays

    TBLiteRuntimeError
        if a requested quantity is not available in the container
    """

    _res = library.ffi.NULL
    _natoms = 0
    _getter = {
        "energy": library.get_energy,
        "energies": library.get_energies,
        "gradient": library.get_gradient,
        "virial": library.get_virial,
        "charges": library.get_charges,
        "bond-orders": library.get_bond_orders,
        "dipole": library.get_dipole,
        "quadrupole": library.get_quadrupole,
        "orbital-energies": library.get_orbital_energies,
        "orbital-occupations": library.get_orbital_occupations,
        "orbital-coefficients": library.get_orbital_coefficients,
        "density-matrix": library.get_density_matrix,
        "overlap-matrix": library.get_overlap_matrix,
        "hamiltonian-matrix": library.get_hamiltonian_matrix,
        "post-processing-dict": library.get_post_processing_dict,
        "natoms": library.get_number_of_atoms,
        "norbitals": library.get_number_of_orbitals,
    }
    _setter = {}

    def __init__(self, other=None):
        """
        Instantiate a new Result container, in case an existing container is passed
        all existing results will be copied into the newly created container as well.
        """
        if other is not None:
            self._res = library.copy_result(other._res)
        else:
            self._res = library.new_result()

    def get(self, attribute: str):
        """
        Get a quantity stored instade the result container.
        The following quantities are available

        ====================== ================================= ==============
         property               dimension [spin-polarized case]   unit
        ====================== ================================= ==============
         energy                 scalar                            Hartree
         energies               nat                               Hartree
         gradient               nat, 3                            Hartree/Bohr
         virial                 3, 3                              Hartree
         charges                nat                               e
         bond-orders            nat, nat                          e
         dipole                 3                                 e·Bohr
         quadrupole             6                                 e·Bohr²
         orbital-energies       norb [2, norb]                    Hartree
         orbital-occupations    norb [2, norb]                    e
         orbital-coefficients   norb, norb [2, norb, norb]        unitless
         overlap-matrix         norb, norb                        unitless
         hamiltonian-matrix     norb, norb                        Hartree
         density-matrix         norb, norb [2, norb, norb]        e
         natoms                 scalar                            unitless
         norbitals              scalar                            unitless
         post-processing-dict   dependes on the key               /
        ====================== ================================= ==============

        Notes
        -----
        The Hamiltonian matrix is the core Hamiltonian rather than the
        converged full Hamiltonian after selfconsistency. To reconstruct it
        transform the orbital energies from the MO to the AO basis using
        the orbital coefficients.

        Raises
        ------
        TBLiteValueError
            on invalid input, like incorrect shape / type of the passed arrays

        TBLiteRuntimeError
            if a requested quantity is not available in the container
        """

        if attribute not in self._getter:
            raise TBLiteValueError(
                f"Attribute '{attribute}' is not available in this result"
            )

        return self._getter[attribute](self._res)

    def set(self, attribute: str, value):
        """
        Get a quantity stored instade the result container.
        Currently, no quantities can be set in the result container.

        Raises
        ------
        TBLiteValueError
            on invalid input, like incorrect shape / type of the passed arrays

        TBLiteRuntimeError
            if a requested quantity cannot be set in the container
        """

        if attribute not in self._setter:
            raise TBLiteValueError(
                f"Attribute '{attribute}' cannot be set in this result"
            )

        self._setter[attribute](self._res, value)

    def save(self, filename: str) -> None:
        """
        Save the result container to a file. The format is compatible with the
        npz format, which is used by the numpy library and can be read
        by the numpy.load function.

        Raises
        ------
        TBLiteRuntimeError
            if the saving fails
        """
        library.save_wavefunction(self._res, filename)
    
    def load(self, filename: str) -> None:
        """
        Load a result container from a file. The format is compatible with the
        npz format, which is used by the numpy library and can be written
        by the numpy.savez function.

        Raises
        ------
        TBLiteRuntimeError
            if the loading fails
        """
        library.load_wavefunction(self._res, filename)

    def __getitem__(self, key: str) -> Any:
        return self.get(key)

    def __setitem__(self, key: str, value: Any) -> None:
        self.set(key, value)

    def dict(self) -> dict:
        """
        Return all quantities inside the result container as dict. In case no
        results are present an empty dict is returned.
        """
        res = {}

        for key in self._getter:
            try:
                res[key] = self.get(key)
            except (TBLiteValueError, TBLiteRuntimeError):
                pass

        return res


class Calculator(Structure):
    """
    .. Single point calculator

    Represents a wrapped calculator object in ``tblite`` and the associated structure data.
    The calculator is instantiated for the respective structure and is immutable once
    created. The cartesian coordinates and the lattice parameters of the structure can
    be updated, while changing the boundary conditions, atomic types or number of atoms
    require the complete reconstruction of the calculator instance.

    Available methods to parametrization of the calculator are:

    **GFN2-xTB**:

    Self-consistent extended tight binding Hamiltonian with
    anisotropic second order electrostatic contributions,
    third order on-site contributions and self-consistent D4 dispersion.

    Geometry, frequency and non-covalent interactions parametrisation for
    elements up to Z=86.

    Cite as:

    * C. Bannwarth, S. Ehlert and S. Grimme.,
      *J. Chem. Theory Comput.* (2019), **15**, 1652-1671.
      DOI: `10.1021/acs.jctc.8b01176 <https://dx.doi.org/10.1021/acs.jctc.8b01176>`_

    **GFN1-xTB**:

    Self-consistent extended tight binding Hamiltonian with
    isotropic second order electrostatic contributions and
    third order on-site contributions.

    Geometry, frequency and non-covalent interactions parametrisation for
    elements up to Z=86.

    Cite as:

    * S. Grimme, C. Bannwarth, P. Shushkov,
      *J. Chem. Theory Comput.* (2017), **13**, 1989-2009.
      DOI: `10.1021/acs.jctc.7b00118 <https://dx.doi.org/10.1021/acs.jctc.7b00118>`_

    **IPEA1-xTB**:

    Special parametrisation for the GFN1-xTB Hamiltonian to improve the
    description of vertical ionisation potentials and electron affinities.
    Uses additional diffuse s-functions on light main group elements.
    Parametrised up to Z=86.

    Cite as:

    * V. Asgeirsson, C. Bauer and S. Grimme, *Chem. Sci.* (2017), **8**, 4879.
      DOI: `10.1039/c7sc00601b <https://dx.doi.org/10.1039/c7sc00601b>`_

    Example
    -------
    >>> from tblite.interface import Calculator
    >>> import numpy as np
    >>> numbers = np.array([1, 1, 6, 5, 1, 15, 8, 17, 13, 15, 5, 1, 9, 15, 1, 15])
    >>> positions = np.array([  # Coordinates in Bohr
    ...     [+2.79274810283778, +3.82998228828316, -2.79287054959216],
    ...     [-1.43447454186833, +0.43418729987882, +5.53854345129809],
    ...     [-3.26268343665218, -2.50644032426151, -1.56631149351046],
    ...     [+2.14548759959147, -0.88798018953965, -2.24592534506187],
    ...     [-4.30233097423181, -3.93631518670031, -0.48930754109119],
    ...     [+0.06107643564880, -3.82467931731366, -2.22333344469482],
    ...     [+0.41168550401858, +0.58105573172764, +5.56854609916143],
    ...     [+4.41363836635653, +3.92515871809283, +2.57961724984000],
    ...     [+1.33707758998700, +1.40194471661647, +1.97530004949523],
    ...     [+3.08342709834868, +1.72520024666801, -4.42666116106828],
    ...     [-3.02346932078505, +0.04438199934191, -0.27636197425010],
    ...     [+1.11508390868455, -0.97617412809198, +6.25462847718180],
    ...     [+0.61938955433011, +2.17903547389232, -6.21279842416963],
    ...     [-2.67491681346835, +3.00175899761859, +1.05038813614845],
    ...     [-4.13181080289514, -2.34226739863660, -3.44356159392859],
    ...     [+2.85007173009739, -2.64884892757600, +0.71010806424206],
    ... ])
    >>> calc = Calculator("GFN2-xTB", numbers, positions)
    >>> res = calc.singlepoint()
    >>> res.get("energy")  # Results in atomic units
    -31.716159156026254

    Raises
    ------
    TBLiteValueError
        on invalid input, like incorrect shape / type of the passed arrays
    """

    _ctx = library.ffi.NULL
    _calc = library.ffi.NULL
    _loader = {
        "GFN2-xTB": library.new_gfn2_calculator,
        "GFN1-xTB": library.new_gfn1_calculator,
        "IPEA1-xTB": library.new_ipea1_calculator,
    }
    _setter = {
        "accuracy": library.set_calculator_accuracy,
        "guess": library.set_calculator_guess,
        "max-iter": library.set_calculator_max_iter,
        "mixer-damping": library.set_calculator_mixer_damping,
        "save-integrals": library.set_calculator_save_integrals,
        "temperature": library.set_calculator_temperature,
        "verbosity": library.set_calculator_verbosity,
    }
    _getter = {
        "angular-momenta": library.get_calculator_angular_momenta,
        "orbital-map": library.get_calculator_orbital_map,
        "shell-map": library.get_calculator_shell_map,
    }
    _interaction = {
        "electric-field": library.new_electric_field,
        "spin-polarization": library.new_spin_polarization,
        "alpb-solvation": library.new_alpb_solvation,
        "gbsa-solvation": library.new_gbsa_solvation,
        "cpcm-solvation": library.new_cpcm_solvation,
        "gbe-solvation": library.new_gbe_solvation,
        "gb-solvation": library.new_gb_solvation,
    }
    _post_processing = {
        "bond-orders" : "bond-orders",
        "molecular-multipoles" : "molmom",
        "xtbml" : "xtbml",
        "xtbml_xyz" : "xtbml_xyz"
    }

    def __init__(
        self,
        method: str,
        numbers: Union[np.ndarray, List[int]],
        positions: np.ndarray,
        charge: Optional[float] = None,
        uhf: Optional[int] = None,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
        **context_kwargs,
    ):
        """
        Construct new calculator object for a given structure.

        Raises
        ------
        TBLiteValueError
            on invalid input, like incorrect shape / type of the passed arrays
        """
        Structure.__init__(
            self, numbers, positions, charge, uhf, lattice, periodic
        )

        self._ctx = library.new_context(**context_kwargs)
        if method not in self._loader:
            raise TBLiteValueError(
                f"Method '{method}' is not available for this calculator"
            )
        self._calc = self._loader[method](self._ctx, self._mol)
        self._method = method

    def set(self, attribute: str, value) -> None:
        """
        Set an attribute in the calculator instance. Supported attributes are

        ================= ==================================== =================
         name              description                          default
        ================= ==================================== =================
         accuracy          Numerical thresholds for SCC         1.0
         guess             Initial guess for wavefunction       0 (SAD)
         max-iter          Maximum number of SCC iterations     250
         mixer-damping     Parameter for the SCC mixer          0.4
         save-integrals    Keep integral matrices in results    0 (False)
         temperature       Electronic temperature for filling   9.500e-4
         verbosity         Set verbosity of printout            1
        ================= ==================================== =================

        .. note::

           The electronic temperature is given in Hartree, rather than Kelvin.
           The conversion factor from Kelvin to Hartree is the Boltzmann constant
           in Hartree/Kelvin (3.166808578545117e-6).

        Raises
        ------
        TBLiteValueError
            on invalid input, like incorrect shape / type of the passed arrays
        """

        if attribute not in self._setter:
            raise TBLiteValueError(
                f"Attribute '{attribute}' is not supported in this calculator"
            )
        self._setter[attribute](self._ctx, self._calc, value)

    def add(self, interaction, *args) -> None:
        """
        Add an interaction to the calculator instance. Supported interactions are

        =================== =========================== =========================================
         name                description                 Arguments
        =================== =========================== =========================================
         electric-field      Uniform electric field      Field vector (3,)
         spin-polarization   Spin polarization           Scaling factor
         alpb-solvation      ALPB implicit solvation     Solvent name, solution state (optional)
         gbsa-solvation      GBSA implicit solvation     Solvent name, solution state (optional)
         cpcm-solvation      CPCM implicit solvation     Epsilon
         gbe-solvation       GBε implicit solvation      Epsilon, Born kernel
         gb-solvation        GB implicit solvation       Epsilon, Born kernel
        =================== =========================== =========================================

        .. note::

            For GSBA and ALPB:
            For named solvents, uses parametrized GBSA/ALPB with CDS and empirical shift.
            For unnamed solvents (dielectric constant), uses non-empirical GBSA/ALPB.
            Optional solution state correction: gsolv (default), bar1mol, reference.
        """

        if interaction in self._interaction:
            kwargs = {}
            if interaction in ("alpb-solvation", "gbsa-solvation"):
                kwargs["version"] = {
                    "GFN2-xTB": 2,
                    "IPEA1-xTB": 1,
                    "GFN1-xTB": 1,
                }[self._method]
            cont = self._interaction[interaction](
                self._ctx, self._mol, self._calc, *args, **kwargs
            )
            library.calculator_push_back(self._ctx, self._calc, cont)
        elif interaction in self._post_processing:
            library.post_processing_push_back(
                self._ctx, self._calc, self._post_processing[interaction]
            )
        elif ".toml" in interaction:
            library.post_processing_push_back(
                self._ctx, self._calc, interaction
            )
        else:
            raise TBLiteValueError(
                f"Interaction or post processing '{interaction}' is not supported in this calculator"
            )

    def get(self, attribute: str) -> Any:
        """
        Set an attribute from the calculator instance. Supported attributes are

        ================= ====================================
         name              description
        ================= ====================================
         angular-momenta   Angular momenta of shells
         orbital-map       Mapping from orbitals to shells
         shell-map         Mapping from shells to atoms
        ================= ====================================

        Raises
        ------
        TBLiteValueError
            on invalid attributes
        """

        if attribute not in self._getter:
            raise TBLiteValueError(
                f"Attribute '{attribute}' is not supported in this calculator"
            )
        return self._getter[attribute](self._ctx, self._calc)
        

    def singlepoint(
        self, res: Optional[Result] = None, copy: bool = False
    ) -> Result:
        """
        Perform actual single point calculation in the library backend.
        The output of the library will be forwarded to the standard output.

        The routine returns an object containing the results, which can be used
        to restart the calculation. Unless specified the restart object will be
        updated inplace rather than copied.

        Raises
        ------
        TBLiteRuntimeError
            in case the calculation fails
        """

        _res = Result(res) if copy or res is None else res
        library.get_singlepoint(self._ctx, self._mol, self._calc, _res._res)
        return _res


def _cast(ctype, array):
    """Cast a numpy array to an FFI pointer"""
    return (
        library.ffi.NULL
        if array is None
        else library.ffi.cast(ctype, array.ctypes.data)
    )


def _ref(ctype, value):
    """Create a reference to a value"""
    if value is None:
        return library.ffi.NULL
    ref = library.ffi.new(ctype + "*")
    ref[0] = value
    return ref


# fmt: off
ELEMENT_SYMBOLS = [
    *["H", "He"],
    *["Li", "Be", "B", "C", "N", "O", "F", "Ne"],
    *["Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"],
    *["K", "Ca"],
    *["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"],
    *["Ga", "Ge", "As", "Se", "Br", "Kr"],
    *["Rb", "Sr"],
    *["Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"],
    *["In", "Sn", "Sb", "Te", "I", "Xe"],
    *["Cs", "Ba"],
    *["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"],
    *["Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"],
    *["Tl", "Pb", "Bi", "Po", "At", "Rn"],
    *["Fr", "Ra"],
    *["Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No"],
    *["Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn"],
    *["Nh", "Fl", "Mc", "Lv", "Ts", "Og"],
]
# fmt: on

SYMBOL_TO_NUMBER = {
    symbol: number + 1 for number, symbol in enumerate(ELEMENT_SYMBOLS)
}


def symbols_to_numbers(symbols: List[str]) -> List[int]:
    """Convert a list of atomic symbols to atomic numbers."""
    return [SYMBOL_TO_NUMBER[symbol] for symbol in symbols]
