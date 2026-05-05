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
Wrapper around the C-API of the dftd4 shared library.
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

    Represents a wrapped structure object in ``dftd4``.
    The molecular structure data object has a fixed number of atoms
    and immutable atomic identifiers

    Example
    -------
    >>> from dftd4.interface import Structure
    >>> import numpy as np
    >>> mol = Structure(
    ...     positions=np.array([
    ...         [+0.00000000000000, +0.00000000000000, -0.73578586109551],
    ...         [+1.44183152868459, +0.00000000000000, +0.36789293054775],
    ...         [-1.44183152868459, +0.00000000000000, +0.36789293054775],
    ...     ]),
    ...     numbers = np.array([8, 1, 1]),
    ... )
    >>> len(mol)
    3

    Raises
    ------
    ValueError
        on invalid input, like incorrect shape / type of the passed arrays
    """

    _mol = library.ffi.NULL

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        charge: Optional[float] = None,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
    ):
        """Create new molecular structure data"""
        if positions.size % 3 != 0:
            raise ValueError("Expected triples of cartesian coordinates")

        if 3 * numbers.size != positions.size:
            raise ValueError("Dimension missmatch between numbers and positions")

        self._natoms = len(numbers)
        _numbers = np.ascontiguousarray(numbers, dtype="i4")
        _positions = np.ascontiguousarray(positions, dtype=float)

        _charge = _ref("double", charge)

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
            _cast("int*", _numbers),
            _cast("double*", _positions),
            _charge,
            _cast("double*", _lattice),
            _cast("bool*", _periodic),
        )

    def __len__(self):
        return self._natoms

    def update(
        self,
        positions: np.ndarray,
        lattice: Optional[np.ndarray] = None,
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
            _cast("double*", _positions),
            _cast("double*", _lattice),
        )

class DampingFunction:
    """
    .. Damping function definition

    Represents the mathematical form of the short-range damping of the two-body
    and three-body dispersion interaction. The damping function is immutable
    after contruction. To change the damping function, a new object must be created.

    Raises
    ------
    ValueError
        on invalid damping function or method name
    RuntimeError
        failed to construct damping function object in API
    """

    _damp = library.ffi.NULL

    def __init__(self, **kwargs):
        """
        Create a new damping function from explicit names for the two- and 
        three-body damping or from the default of a specified dispersion model.
        Setting the three-body damping function to 'none' deactivates the ATM term.

        Example
        -------
        >>> from dftd4.interface import DampingFunction
        >>> # From explicit names for two- and three-body damping:
        >>> damp = DampingFunction(damping_2b="rational", damping_3b="zero-avg")
        >>> # From the defaults of a dispersion model:
        >>> damp = DampingFunction(model="d4")
        >>> # From the defaults of a dispersion model and deactivated three-body term:
        >>> damp = DampingFunction(model="d4", damping_3b="none")

        Two-body damping functions: 'rational', 'screened', 'zero', 'mzero',
                                    'optpower', 'cso', 'koide'
        Three-body damping functions: 'none', 'rational', 'screened', 'zero',
                                      'zero-avg'
        """

        if "damping_2b" in kwargs and "model" in kwargs:
            raise ValueError("Cannot provide both two-body damping and a model.")
        
        if "damping_2b" in kwargs:
            # Explicitly specified two-body (and three-body) damping functions
            _d2 = kwargs.get("damping_2b")
            _d3 = kwargs.get("damping_3b", "none")
        elif "model" in kwargs:
            # Default damping functions for the specified model
            # (optional modification of the three-body damping to deactivate the ATM term)
            _model = kwargs.pop("model").lower().replace(" ", "")
            if _model in ("d4", "d4s"):
                _d2 = library.DEFAULT_TWOBODY_DAMPING[_model]
                _d3 = kwargs.get("damping_3b", library.DEFAULT_THREEBODY_DAMPING[_model])
            else:
                raise ValueError(f"Unknown dispersion model '{_model}'.")
        else:
            raise ValueError("Either two-body damping or a model name is required.")

        try:
            self._damp = library.new_damping(_d2, _d3)
        except KeyError as e:
            raise ValueError(f"Invalid damping type: {e.args[0]}") from e


    def check_params(self, param: 'DampingParam') -> None:
          """
          Check if the provided damping parameters are compatible with the damping
          function. Will silently pass if the parameters are compatible,
          or raise an error.
          
          Example
          -------
          >>> from dftd4.interface import DampingFunction, DampingParam
          >>> damp = DampingFunction(damping_2b="rational", damping_3b="zero-avg")
          >>> param = DampingParam(model="d4", s6=0.6400, s8=1.16888646, a1=0.44154604,
                                   a2=4.73114642)
          >>> damp.check_params(param)

          Raises
          ------
          RuntimeError
              incompatible or missing parameters required by the damping function
          """
          library.check_params(self._damp, param._param)


class DampingParam:
    """
    .. Damping parameters for the damping function.

    The damping parameters contained in the object are immutable. To change the
    parametrization, a new object must be created. Furthermore, the object is
    opaque to the user and the contained data cannot be accessed directly.

    There are three ways provided to generate a new damping parameter object:

    1. a method name is passed to the constructor, the library will load the
       required data from the *dftd4* shared library.

    2. a model name and all required parameters for the default damping function
       of this dispersion model are passed. The library will generate 
       a compatible object from the specified parameters.

    3. all possible parameters are passed to the constructor without a model name
       and the library will generate a general parameter object specifying
       all parameters.

    Raises
    ------
    TypeError
        incorrect/missing input values provided to constructor
    ValueError
        invalid model, parameters, or damping types provided
    RuntimeError
        failed to construct damping parameter object in API
    """

    _param = library.ffi.NULL

    def __init__(self, **kwargs):
        """Create new damping parameter object from method name or explicit data"""

        if not kwargs:
            raise TypeError("Method name or explicit damping parameters are required.")

        if "method" in kwargs:
            _method = kwargs.pop("method")

            if any (key in kwargs for key in ("s6", "s8", "s9",
                                              "a1", "a2", "a3", "a4",
                                              "rs6", "rs8", "rs9",
                                              "alp", "bet")):
                raise ValueError("Explicit parameters cannot be mixed with a method name.")

            if "model" not in kwargs:
                raise ValueError("Model name must be provided when using method name.")
            _model = kwargs.get("model").lower().replace(" ", "")

            
            if "damping_2b" in kwargs:
                # Explicitly specified two-body (and three-body) damping functions
                _d2 = kwargs.get("damping_2b", library.DEFAULT_TWOBODY_DAMPING[_model])
                _d3 = kwargs.get("damping_3b", "none")
            elif _model in ("d4", "d4s"):
                # Default damping functions for the specified model 
                # (optionally modify three-body damping to deactivate the ATM term) 
                _d2 = library.DEFAULT_TWOBODY_DAMPING[_model]
                _d3 = kwargs.get("damping_3b", library.DEFAULT_THREEBODY_DAMPING[_model])
            else:
                raise ValueError(f"Unknown dispersion model '{_model}'.")

            self._param = self.load_param(_method, _model, _d2, _d3)

        elif "model" in kwargs:
            _model = kwargs.pop("model").lower().replace(" ", "")
            if _model in ("d4", "d4s"):
                # Specify parameters for the default D4 damping (rational + zero-avg)
                self._param = self.new_d4_default_param(**kwargs)
            else:
                raise ValueError(f"Unknown dispersion model '{_model}'.")

        else:
            # Specify all possible parameters
            self._param = self.new_param(**kwargs)

    @staticmethod
    def load_param(method: str, model: str, damping_2b: str, damping_3b: str):
        """
        Create damping parameter API object from internal library storage for a given
        model and damping function combination by searching for the provided method
        name. The method name is case insensitive and hyphens are ignored. Setting
        the three-body damping function to 'none' deactivates the ATM term. In case
        the method name is unknown an exception is raised.

        Example
        -------
        >>> from dftd4.interface import DampingParam
        >>> # From the method name and defaults of the model:
        >>> param = DampingParam(method="pbe", model="d4")
        >>> # From the method name and defaults of the model without ATM:
        >>> param = DampingParam(method="pbe", model="d4", damping_3b="none")

        Raises
        ------
        RuntimeError
            failed to construct damping parameter object in API
        ValueError
            invalid model or damping type provided
        """

        try:
            return library.load_param(method, model, damping_2b, damping_3b)
        except KeyError as e:
            raise ValueError(f"Invalid model or damping type: {e.args[0]}") from e

    @staticmethod
    def new_d4_default_param(*, s6=1.0, s8, s9=1.0, a1, a2, rs9=1.0, alp=16.0):
        """
        Create D4 default damping parameter API object from user provided
        parameters. This object contains damping parameters for a rational
        two-body and averaged zero three-body damping function and requires
        at least the 's8', 'a1', and 'a2' parameters as input. Additionally,
        the parameters 's6', 's9', 'rs9', and 'alp' can be overwritten.
        The user provided damping parameters will be used unchecked.

        Example
        -------
        >>> from dftd4.interface import DampingParam
        >>> param = DampingParam(model="d4", s6=0.6400, s8=1.16888646, a1=0.44154604,
                                 a2=4.73114642)

        Raises
        ------
        RuntimeError
            failed to construct damping parameter object in API
        """
        return library.new_param(s6, s8, s9, a1, a2, 0.0, 0.0, 0.0, 0.0, rs9, alp, 0.0)
    
    @staticmethod
    def new_param(*, s6=1.0, s8, s9=1.0, a1, a2, a3, a4, rs6, rs8, rs9, alp, bet):
        """
        Create damping parameter API object from full set of user provided
        parameters. The user provided damping parameters will be used unchecked.

        Example
        -------
        >>> from dftd4.interface import DampingParam
        >>> param = DampingParam(s6=0.6400, s8=1.16888646, s9=1.0, a1=0.44154604,
                                 a2=4.73114642, a3=0.0, a4=0.0, rs6=0.0, rs8=0.0,
                                 rs9=1.0, alp=16.0, bet=0.0)

        Raises
        ------
        RuntimeError
            failed to construct damping parameter object in API
        """
        return library.new_param(s6, s8, s9, a1, a2, a3, a4, rs6, rs8, rs9, alp, bet)

class DispersionModel(Structure):
    """
    .. Dispersion model

    Representation of a dispersion model to evaluate the dispersion
    coefficients and energies. The model is coupled to the molecular
    structure it has been created from and cannot be transfered to
    another molecular structure without recreating it.

    Example
    -------
    >>> from dftd4.interface import DispersionModel
    >>> import numpy as np
    >>> disp = DispersionModel(
    ...     positions=np.array([  # Coordinates in Bohr
    ...         [+0.00000000000000, +0.00000000000000, -0.73578586109551],
    ...         [+1.44183152868459, +0.00000000000000, +0.36789293054775],
    ...         [-1.44183152868459, +0.00000000000000, +0.36789293054775],
    ...     ]),
    ...     numbers = np.array([8, 1, 1]),
    ... )
    >>> disp.get_properties()["polarizabilities"]
    array([6.74893641, 1.33914933, 1.33914933])

    Raises
    ------
    ValueError
        on incorrect inputs to for the structure data

    RuntimeError
        in case of an error in library
    """

    _disp = library.ffi.NULL

    def __init__(
        self,
        numbers: np.ndarray,
        positions: np.ndarray,
        charge: Optional[float] = None,
        lattice: Optional[np.ndarray] = None,
        periodic: Optional[np.ndarray] = None,
        model: str = "d4",
        **kwargs,
    ):
        """Create new dispersion model"""

        Structure.__init__(self, numbers, positions, charge, lattice, periodic)
        

        if model.lower().replace(" ", "") == "d4":
            if "ga" in kwargs or "gc" in kwargs or "wf" in kwargs:
                self._disp = library.custom_d4_model(
                    self._mol,
                    kwargs.get("ga", 3.0),
                    kwargs.get("gc", 2.0),
                    kwargs.get("wf", 6.0),
                )
            else:
                self._disp = library.new_d4_model(self._mol)
        elif model.lower().replace(" ", "") == "d4s":
            if "ga" in kwargs or "gc" in kwargs:
                self._disp = library.custom_d4s_model(
                    self._mol,
                    kwargs.get("ga", 3.0),
                    kwargs.get("gc", 2.0),
                )
            else: 
                self._disp = library.new_d4s_model(self._mol)
        else: 
            raise ValueError(f"Unknown dispersion model '{model}'.")

    def get_dispersion(self, damp: DampingFunction, param: DampingParam, grad: bool) -> dict:
        """
        Perform actual evaluation of the dispersion correction.

        Example
        -------
        >>> from dftd4.interface import DampingParam, DampingFunction, DispersionModel
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
        >>> model = DispersionModel(numbers, positions)
        >>> damp = DampingFunction(model="d4")
        >>> param = DampingParam(method="scan", model="d4")
        >>> res = model.get_dispersion(damp, param, grad=False)
        >>> res.get("energy")  # Results in atomic units
        -0.005328888532435093

        Raises
        ------
        RuntimeError
            in case the calculation fails in the library
        """

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
            damp._damp,
            param._param,
            _cast("double*", _energy),
            _cast("double*", _gradient),
            _cast("double*", _sigma),
        )

        results = dict(energy=_energy)
        if _gradient is not None:
            results.update(gradient=_gradient)
        if _sigma is not None:
            results.update(virial=_sigma)
        return results

    def get_properties(self) -> dict:
        """
        Evaluate dispersion related properties, like polarizabilities and C6 coefficients.
        Will also return the coordination numbers and partial charges used to derive
        the polarizabilities. Both dipole-dipole and quadrupole-quadrupole static
        polarizabilities are returned.

        Example
        -------
        >>> from dftd4.interface import DispersionModel
        >>> import numpy as np
        >>> disp = DispersionModel(
        ...     numbers=np.array([16, 16, 16, 16, 16, 16, 16, 16]),
        ...     positions=np.array([
        ...         [-4.15128787379191, +1.71951973863958, -0.93066267097296],
        ...         [-4.15128787379191, -1.71951973863958, +0.93066267097296],
        ...         [-1.71951973863958, -4.15128787379191, -0.93066267097296],
        ...         [+1.71951973863958, -4.15128787379191, +0.93066267097296],
        ...         [+4.15128787379191, -1.71951973863958, -0.93066267097296],
        ...         [+4.15128787379191, +1.71951973863958, +0.93066267097296],
        ...         [+1.71951973863958, +4.15128787379191, -0.93066267097296],
        ...         [-1.71951973863958, +4.15128787379191, +0.93066267097296],
        ...     ]),
        ... )
        >>> res = disp.get_properties()
        >>> res.get("coordination numbers")
        array([1.96273847, 1.96273847, 1.96273847, 1.96273847, 1.96273847,
               1.96273847, 1.96273847, 1.96273847])
        >>> res.get("polarizabilities").sum()
        158.748605606818
        """

        _c6 = np.zeros((len(self), len(self)))
        _cn = np.zeros((len(self)))
        _charges = np.zeros((len(self)))
        _alpha = np.zeros((len(self)))
        _alphaqq = np.zeros((len(self)))

        library.get_properties(
            self._mol,
            self._disp,
            _cast("double*", _cn),
            _cast("double*", _charges),
            _cast("double*", _c6),
            _cast("double*", _alpha),
            _cast("double*", _alphaqq),
        )

        return {
            "coordination numbers": _cn,
            "partial charges": _charges,
            "c6 coefficients": _c6,
            "polarizabilities": _alpha,
            "polarizabilities quad-quad": _alphaqq,
        }

    def get_pairwise_dispersion(self, damp: DampingFunction, param: DampingParam) -> dict:
        """
        Evaluate pairwise representation of the dispersion energy

        >>> from dftd4.interface import DispersionModel, DampingParam
        >>> import numpy as np
        >>> disp = DispersionModel(
        ...     numbers=np.array([7, 7, 1, 1, 1, 1, 1, 1]),
        ...     positions=np.array([
        ...         [-2.983345508575, -0.088082052767, +0.000000000000],
        ...         [+2.983345508575, +0.088082052767, +0.000000000000],
        ...         [-4.079203605652, +0.257751166821, +1.529856562614],
        ...         [-1.605268001556, +1.243804812431, +0.000000000000],
        ...         [-4.079203605652, +0.257751166821, -1.529856562614],
        ...         [+4.079203605652, -0.257751166821, -1.529856562614],
        ...         [+1.605268001556, -1.243804812431, +0.000000000000],
        ...         [+4.079203605652, -0.257751166821, +1.529856562614],
        ...     ]),
        ... )
        >>> damp = DampingFunction(model="d4")
        >>> param = DampingParam(method="tpss", model="d4")
        >>> res = disp.get_pairwise_dispersion(damp, param)
        >>> res["additive pairwise energy"].sum()
        -0.0023605238432524104
        >>> res["non-additive pairwise energy"].sum()
        8.794562567135391e-08
        """

        _pair_disp2 = np.zeros((len(self), len(self)))
        _pair_disp3 = np.zeros((len(self), len(self)))

        library.get_pairwise_dispersion(
            self._mol,
            self._disp,
            damp._damp,
            param._param,
            _cast("double*", _pair_disp2),
            _cast("double*", _pair_disp3),
        )

        return {
            "additive pairwise energy": _pair_disp2,
            "non-additive pairwise energy": _pair_disp3,
        }

    def get_numerical_hessian(self, damp: DampingFunction, param: DampingParam) -> np.ndarray:
        """
        Evaluate the numerical Hessian of the dispersion energy.

        Example
        -------
        >>> from dftd4.interface import DampingParam, DampingFunction, DispersionModel
        >>> import numpy as np
        >>> disp = DispersionModel(
        ...     numbers=np.array([7, 7, 1, 1, 1, 1, 1, 1]),
        ...     positions=np.array([
        ...         [-2.983345508575, -0.088082052767, +0.000000000000],
        ...         [+2.983345508575, +0.088082052767, +0.000000000000],
        ...         [-4.079203605652, +0.257751166821, +1.529856562614],
        ...         [-1.605268001556, +1.243804812431, +0.000000000000],
        ...         [-4.079203605652, +0.257751166821, -1.529856562614],
        ...         [+4.079203605652, -0.257751166821, -1.529856562614],
        ...         [+1.605268001556, -1.243804812431, +0.000000000000],
        ...         [+4.079203605652, -0.257751166821, +1.529856562614],
        ...     ]),
        ... )
        >>> damp = DampingFunction(model="d4")
        >>> param = DampingParam(method="blyp", model="d4")
        >>> res = disp.get_numerical_hessian(damp, param)
        >>> res["hessian"].shape
        (8, 3, 8, 3)
        >>> res["hessian"][0, :, 0, :]
        array([[-1.39880830e-04, -8.26672018e-05,  1.76708940e-15],
               [-8.26672036e-05, -1.68686558e-04,  2.70024839e-16],
               [ 1.08420217e-15,  1.66018458e-15, -9.37126387e-05]])
        """

        _hessian = np.zeros((len(self), 3, len(self), 3))

        library.get_numerical_hessian(
            self._mol,
            self._disp,
            damp._damp,
            param._param,
            _cast("double*", _hessian),
        )

        return {
            "hessian": _hessian
        }


def _cast(ctype, array):
    """Cast a numpy array to a FFI pointer"""
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
