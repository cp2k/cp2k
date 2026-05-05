# This file is part of s-dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.
"""
Wrapper around the C-API of the s-dftd3 shared library.
It provides the definition the basic interface to the library for most further integration
in other Python frameworks.

The classes defined here allow a more Pythonic usage of the API object provided by the
library in actual workflows than the low-level access provided in the CFFI generated wrappers.
"""

from typing import Optional
import numpy as np

from . import library


class Structure:
    """
    .. Molecular structure data

    Represents a wrapped structure object in ``s-dftd3``.
    The molecular structure data object has a fixed number of atoms
    and immutable atomic identifiers

    Example
    -------
    >>> from dftd3.interface import Structure
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
    ValueError
        on invalid input, like incorrect shape / type of the passed arrays
    """

    _mol = library.StructureHandle.null()

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
    ):
        """
        Create new molecular structure data from arrays. The returned object has
        immutable atomic species and boundary condition, also the total number of
        atoms cannot be changed.

        Raises
        ------
        ValueError
            on invalid input, like incorrect shape / type of the passed arrays
        """

        if positions.size % 3 != 0:
            raise ValueError("Expected tripels of cartesian coordinates")

        if 3 * numbers.size != positions.size:
            raise ValueError("Dimension missmatch between numbers and positions")

        self._natoms = len(numbers)
        _numbers = np.ascontiguousarray(numbers, dtype="i4")
        _positions = np.ascontiguousarray(positions, dtype=float)

        if lattice is not None:
            if lattice.size != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype="float")
        else:
            _lattice = None

        if periodic is not None:
            if periodic.size != 3:
                raise ValueError("Invalid periodicity provided")
            _periodic = np.ascontiguousarray(periodic, dtype="bool")
        else:
            _periodic = None

        self._mol = library.new_structure(
            self._natoms,
            _numbers,
            _positions,
            _lattice,
            _periodic,
        )

    def __len__(self):
        return self._natoms

    def update(
        self,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
    ) -> None:
        """
        Update coordinates and lattice parameters, both provided in
        atomic units (Bohr).
        The lattice update is optional also for periodic structures.

        Generally, only the cartesian coordinates and the lattice parameters
        can be updated, every other modification,
        boundary condition, atomic types or number of atoms
        requires the complete reconstruction of the object.

        Raises
        ------
        ValueError
            on invalid input, like incorrect shape / type of the passed arrays
        """

        if 3 * len(self) != positions.size:
            raise ValueError("Dimension missmatch for positions")
        _positions = np.ascontiguousarray(positions, dtype="float")

        if lattice is not None:
            if lattice.size != 9:
                raise ValueError("Invalid lattice provided")
            _lattice = np.ascontiguousarray(lattice, dtype="float")
        else:
            _lattice = None

        library.update_structure(
            self._mol,
            _positions,
            _lattice,
        )


class DampingParam:
    """
    Abstract base class for damping parameters, representing a parametrization of
    a DFT-D3 method.

    The damping parameters contained in the object are immutable. To change the
    parametrization, a new object must be created. Furthermore, the object is
    opaque to the user and the contained data cannot be accessed directly.

    There are two main ways provided to generate a new damping parameter object:

    1. a method name is passed to the constructor, the library will load the
       required data from the *s-dftd3* shared library.

    2. all required parameters are passed to the constructor and the library will
       generate an object from the given parameters.

    .. note::

       Mixing of the two methods is not allowed to avoid partial initialization
       of any created objects. Users who need full control over the creation
       of the object should use the second method.
    """

    _param = library.ParamHandle.null()

    def __init__(self, **kwargs):
        """Create new damping parameter from method name or explicit data"""

        if "method" in kwargs and kwargs["method"] is None:
            del kwargs["method"]

        if "method" in kwargs:
            self._param = self.load_param(**kwargs)
        else:
            self._param = self.new_param(**kwargs)

    @staticmethod
    def load_param(method, **kwargs) -> library.ParamHandle:
        raise NotImplementedError("Child class has to define parameter loading")

    @staticmethod
    def new_param(**kwargs) -> library.ParamHandle:
        raise NotImplementedError("Child class has to define parameter construction")


class RationalDampingParam(DampingParam):
    r"""
    Rational damping function for DFT-D3.
    The original scheme was proposed by Becke and Johnson\ :footcite:`becke2005,johnson2005,johnson2006`
    and implemented in a slightly adjusted form using only the C8/C6 ratio in the critical
    for DFT-D3.\ :footcite:`grimme2011`
    The rational damping scheme has the advantage of damping the dispersion energy to
    finite value, rather than removing it at short distances.

    .. note::

       The zero damping function is retained for the three-body contributions from the ATM
       term.
    """

    def __init__(self, **kwargs):
        _rename_kwargs(kwargs, "alpha6", "alp")
        DampingParam.__init__(self, **kwargs)

    @staticmethod
    def load_param(method: str, atm: bool = False) -> library.ParamHandle:
        return library.load_rational_damping(
            method,
            atm,
        )

    @staticmethod
    def new_param(
        *,
        s6: float = 1.0,
        s8: float,
        s9: float = 1.0,
        a1: float,
        a2: float,
        alp: float = 14.0,
    ) -> library.ParamHandle:
        return library.new_rational_damping(
            s6,
            s8,
            s9,
            a1,
            a2,
            alp,
        )


class ZeroDampingParam(DampingParam):
    r"""
    Original DFT-D3 damping function,\ :footcite:`grimme2010` based on a variant proposed by
    Chai and Head-Gordon.\ :footcite:`chai2008`
    Since it is damping the dispersion energy to zero at short distances it is usually
    called zero damping scheme for simplicity. However, due to this short-range limit
    of the dispersion energy a repulsive contribution to the gradient can arise, which
    is considered artificial.\ :footcite:`grimme2011`
    """

    def __init__(self, **kwargs):
        _rename_kwargs(kwargs, "sr6", "rs6")
        _rename_kwargs(kwargs, "sr8", "rs8")
        _rename_kwargs(kwargs, "alpha6", "alp")
        DampingParam.__init__(self, **kwargs)

    @staticmethod
    def load_param(method: str, atm: bool = False) -> library.ParamHandle:
        return library.load_zero_damping(
            method,
            atm,
        )

    @staticmethod
    def new_param(
        *,
        s6: float = 1.0,
        s8: float,
        s9: float = 1.0,
        rs6: float,
        rs8: float = 1.0,
        alp: float = 14.0,
    ) -> library.ParamHandle:
        return library.new_zero_damping(
            s6,
            s8,
            s9,
            rs6,
            rs8,
            alp,
        )


class ModifiedRationalDampingParam(DampingParam):
    """
    Modified version of the rational damping parameters. The functional form of the
    damping function is *unmodified* with respect to the original rational damping scheme.
    However, for a number of functionals new parameters were introduced.:footcite:`smith2016`

    This constructor allows to automatically load the reparameterized damping function
    from the library rather than the original one. Providing a full parameter set is
    functionally equivalent to using the `RationalDampingParam` constructor.
    """

    def __init__(self, **kwargs):
        _rename_kwargs(kwargs, "alpha6", "alp")
        DampingParam.__init__(self, **kwargs)

    @staticmethod
    def load_param(method: str, atm: bool = False) -> library.ParamHandle:
        return library.load_mrational_damping(
            method,
            atm,
        )

    @staticmethod
    def new_param(
        *,
        s6: float = 1.0,
        s8: float,
        s9: float = 1.0,
        a1: float,
        a2: float,
        alp: float = 14.0,
    ) -> library.ParamHandle:
        return library.new_mrational_damping(
            s6,
            s8,
            s9,
            a1,
            a2,
            alp,
        )


class ModifiedZeroDampingParam(DampingParam):
    r"""
    Modified zero damping function for DFT-D3.\ :footcite:`smith2016`
    This scheme adds an additional offset parameter to the zero damping scheme
    of the original DFT-D3.

    .. note::

       This damping function is identical to zero damping for ``bet=0.0``.
    """

    def __init__(self, **kwargs):
        _rename_kwargs(kwargs, "sr6", "rs6")
        _rename_kwargs(kwargs, "sr8", "rs8")
        _rename_kwargs(kwargs, "alpha6", "alp")
        _rename_kwargs(kwargs, "beta", "bet")
        DampingParam.__init__(self, **kwargs)

    @staticmethod
    def load_param(method: str, atm: bool = False) -> library.ParamHandle:
        return library.load_mzero_damping(
            method,
            atm,
        )

    @staticmethod
    def new_param(
        *,
        s6: float = 1.0,
        s8: float,
        s9: float = 1.0,
        rs6: float,
        rs8: float = 1.0,
        alp: float = 14.0,
        bet: float,
    ) -> library.ParamHandle:
        return library.new_mzero_damping(
            s6,
            s8,
            s9,
            rs6,
            rs8,
            alp,
            bet,
        )


class OptimizedPowerDampingParam(DampingParam):
    r"""
    Optimized power version of the rational damping parameters.\ :footcite:`witte2017`
    The functional form of the damping function is modified by adding an additional
    zero-damping like power function.

    This constructor allows to automatically load the reparameterized damping function
    from the library rather than the original one. Providing the parameter `bet=0` is
    equivalent to using rational the `RationalDampingParam` constructor.
    """

    def __init__(self, **kwargs):
        _rename_kwargs(kwargs, "alpha6", "alp")
        _rename_kwargs(kwargs, "beta", "bet")
        DampingParam.__init__(self, **kwargs)

    @staticmethod
    def load_param(method: str, atm: bool = False) -> library.ParamHandle:
        return library.load_optimizedpower_damping(
            method,
            atm,
        )

    @staticmethod
    def new_param(
        *,
        s6: float = 1.0,
        s8: float,
        s9: float = 1.0,
        a1: float,
        a2: float,
        alp: float = 14.0,
        bet,
    ) -> library.ParamHandle:
        return library.new_optimizedpower_damping(
            s6,
            s8,
            s9,
            a1,
            a2,
            alp,
            bet,
        )


class CSODampingParam(DampingParam):
    r"""
    CSO (C6-scaled only) damping function.\ :footcite:`schroeder2015`
    Reformulation of the D3(Becke-Johnson) dispersion correction using only C6
    dispersion coefficients with a sigmoidal interpolation function.

    This constructor allows to automatically load the parameterized damping function
    from the library. Only eight functionals are parametrized with this scheme.
    """

    def __init__(self, **kwargs):
        _rename_kwargs(kwargs, "alpha6", "alp")
        DampingParam.__init__(self, **kwargs)

    @staticmethod
    def load_param(method: str, atm: bool = False) -> library.ParamHandle:
        return library.load_cso_damping(
            method,
            atm,
        )

    @staticmethod
    def new_param(
        *,
        s6: float = 1.0,
        s9: float = 1.0,
        a1: float,
        a2: float = 2.5,
        a3: float = 0.0,
        a4: float = 6.25,
        alp: float = 14.0,
    ) -> library.ParamHandle:
        return library.new_cso_damping(
            s6,
            s9,
            a1,
            a2,
            a3,
            a4,
            alp,
        )


class DispersionModel(Structure):
    """
    .. Dispersion model

    Contains the required information to evaluate all dispersion related properties,
    like C6 coefficients. It also manages an instance of the geometry it was constructed
    for to ensure that the dispersion model is always used with the correct structure
    input.
    """

    _disp = library.ModelHandle.null()

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
    ):
        Structure.__init__(self, numbers, positions, lattice, periodic)

        self._disp = library.new_d3_model(self._mol)

    def set_realspace_cutoff(self, disp2: float, disp3: float, cn: float):
        """Set realspace cutoff for evaluation of interactions"""

        library.set_model_realspace_cutoff(self._disp, disp2, disp3, cn)

    def get_dispersion(self, param: DampingParam, grad: bool) -> dict:
        """Perform actual evaluation of the dispersion correction"""

        _energy = np.array(0.0)
        if grad:
            _gradient = np.zeros((len(self), 3))
            _sigma = np.zeros((3, 3))
        else:
            _gradient = None
            _sigma = None

        library.get_dispersion(
            self._mol,
            self._disp,
            param._param,
            _energy,
            _gradient,
            _sigma,
        )

        results = dict(energy=_energy)
        if _gradient is not None:
            results.update(gradient=_gradient)
        if _sigma is not None:
            results.update(virial=_sigma)
        return results

    def get_pairwise_dispersion(self, param: DampingParam) -> dict:
        """Evaluate pairwise representation of the dispersion energy"""

        _pair_disp2 = np.zeros((len(self), len(self)))
        _pair_disp3 = np.zeros((len(self), len(self)))

        library.get_pairwise_dispersion(
            self._mol,
            self._disp,
            param._param,
            _pair_disp2,
            _pair_disp3,
        )

        return {
            "additive pairwise energy": _pair_disp2,
            "non-additive pairwise energy": _pair_disp3,
        }


class GeometricCounterpoise(Structure):
    """
    .. Counterpoise correction parameters

    Contains the required information to evaluate the counterpoise correction
    for a given geometry. The counterpoise correction is a method to correct
    the interaction energy in a supermolecular calculation for the basis set
    superposition error (BSSE).
    """

    _gcp = library.GCPHandle.null()

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
        method: Optional[str] = None,
        basis: Optional[str] = None,
    ):
        Structure.__init__(self, numbers, positions, lattice, periodic)

        self._gcp = library.load_gcp_param(self._mol, method, basis)

    def set_realspace_cutoff(self, bas: float, srb: float):
        """Set realspace cutoff for evaluation of interactions"""

        library.set_gcp_realspace_cutoff(self._gcp, bas, srb)

    def get_counterpoise(self, grad: bool) -> dict:
        """Evaluate the counterpoise corrected interaction energy"""

        _energy = np.array(0.0)
        if grad:
            _gradient = np.zeros((len(self), 3))
            _sigma = np.zeros((3, 3))
        else:
            _gradient = None
            _sigma = None

        library.get_counterpoise(self._mol, self._gcp, _energy, _gradient, _sigma)

        results = dict(energy=_energy)
        if _gradient is not None:
            results.update(gradient=_gradient)
        if _sigma is not None:
            results.update(virial=_sigma)
        return results


def _rename_kwargs(kwargs, old_name, new_name):
    if old_name in kwargs and new_name not in kwargs:
        kwargs[new_name] = kwargs[old_name]
        del kwargs[old_name]
