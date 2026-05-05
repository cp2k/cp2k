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
PySCF Support
=============

Compatibility layer for supporting DFT-D3 in `pyscf <https://pyscf.org/>`_.
"""

try:
    from pyscf import gto, lib, mcscf, scf
    from pyscf.grad import rhf as rhf_grad
except ModuleNotFoundError:
    raise ModuleNotFoundError("This submodule requires pyscf installed")

import numpy as np
from typing import Dict, Optional, Tuple

from .interface import (
    DispersionModel,
    RationalDampingParam,
    ZeroDampingParam,
    ModifiedRationalDampingParam,
    ModifiedZeroDampingParam,
    OptimizedPowerDampingParam,
    CSODampingParam,
)

GradientsBase = getattr(rhf_grad, "GradientsBase", rhf_grad.Gradients)

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


class DFTD3Dispersion(lib.StreamObject):
    """
    Implementation of the interface for using DFT-D3 in pyscf.
    The `xc` functional can be provided in the constructor together with the
    `version` of the DFT-D3 damping function to use.
    Possible damping functions are

    ``"d3bj"``: (default)
        For rational damping function
    ``"d3zero"``
        For zero damping function
    ``"d3mbj"``
        Modified damping parameters for the rational damping function
    ``"d3mzero"``
        Modified version of the zero damping function
    ``"d3op"``
        Optimized power damping function
    ``"d3cso"``
        CSO (C6-scaled only) damping function

    Custom parameters can be provided with the `param` dictionary.
    The `param` dict contains the damping parameters, at least s8, a1 and a2
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

    The version of the damping can be changed after constructing the dispersion correction.
    With the `atm` boolean the three-body dispersion energy can be enabled, which is
    generally recommended.

    Examples
    --------
    >>> from pyscf import gto
    >>> import dftd3.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          C   -0.189833176  -0.645396435   0.069807761
    ...          C    1.121636324  -0.354065576   0.439096514
    ...          C    1.486520953   0.962572632   0.712107225
    ...          C    0.549329390   1.989209324   0.617868956
    ...          C   -0.757627135   1.681862630   0.246856908
    ...          C   -1.138190460   0.370551816  -0.028582325
    ...          Br  -2.038462778   3.070459841   0.115165429
    ...          H    1.852935245  -1.146434699   0.514119204
    ...          H    0.825048723   3.012176989   0.829385472
    ...          H    2.502259769   1.196433556   1.000317333
    ...          H   -2.157140187   0.151608161  -0.313181471
    ...          H   -0.480820487  -1.664983631  -0.142918416
    ...          S   -4.157443472   5.729584377  -0.878761129
    ...          H   -4.823791426   4.796089466  -1.563433338
    ...          C   -2.828338520   5.970593053  -2.091189515
    ...          H   -2.167577293   6.722356639  -1.668621815
    ...          H   -2.264954814   5.054835899  -2.240198499
    ...          H   -3.218524904   6.337447714  -3.035087058
    ...          '''
    ... )
    >>> d3 = disp.DFTD3Dispersion(mol, xc="PW6B95", version="d3bj")
    >>> d3.kernel()[0]
    array(-0.01009386)
    >>> d3.version = "d3zero"  # Change to zero damping
    >>> d3.kernel()[0]
    array(-0.00574098)
    >>> d3.atm = True  # Activate three-body dispersion
    >>> d3.kernel()[0]
    array(-0.00574289)
    """

    def __init__(
        self,
        mol: gto.Mole,
        xc: str = "hf",
        version: str = "d3bj",
        atm: bool = False,
        param: Optional[Dict[str, float]] = None,
    ):
        self.mol = mol
        self.verbose = mol.verbose
        self.xc = xc
        self.param = param
        self.atm = atm
        self.version = version

    def dump_flags(self, verbose: Optional[bool] = None):
        """
        Show options used for the DFT-D3 dispersion correction.
        """
        lib.logger.info(self, "** DFTD3 parameter **")
        lib.logger.info(self, "func %s", self.xc)
        lib.logger.info(
            self, "version %s", self.version + "-atm" if self.atm else self.version
        )
        return self

    def kernel(self) -> Tuple[float, np.ndarray]:
        """
        Compute the DFT-D3 dispersion correction.

        The dispersion model as well as the parameters are created locally and
        not part of the state of the instance.

        Returns
        -------
        float, ndarray
            The energy and gradient of the DFT-D3 dispersion correction.

        Examples
        --------
        >>> from pyscf import gto
        >>> import dftd3.pyscf as disp
        >>> mol = gto.M(
        ...     atom='''
        ...          Br    0.000000    0.000000    1.919978
        ...          Br    0.000000    0.000000   -0.367147
        ...          N     0.000000    0.000000   -3.235006
        ...          C     0.000000    0.000000   -4.376626
        ...          H     0.000000    0.000000   -5.444276
        ...          '''
        ... )
        >>> d3 = disp.DFTD3Dispersion(mol, xc="PBE0")
        >>> energy, gradient = d3.kernel()
        >>> energy
        array(-0.00303589)
        >>> gradient
        array([[ 0.00000000e+00,  0.00000000e+00,  9.66197638e-05],
               [ 0.00000000e+00,  0.00000000e+00,  2.36000434e-04],
               [ 0.00000000e+00,  0.00000000e+00, -1.16718302e-04],
               [ 0.00000000e+00,  0.00000000e+00, -1.84332770e-04],
               [ 0.00000000e+00,  0.00000000e+00, -3.15691249e-05]])
        """
        mol = self.mol

        lattice = None
        periodic = None
        if hasattr(mol, "lattice_vectors"):
            lattice = mol.lattice_vectors()
            periodic = np.array([True, True, True], dtype=bool)

        disp = DispersionModel(
            np.array([gto.charge(mol.atom_symbol(ia)) for ia in range(mol.natm)]),
            mol.atom_coords(),
            lattice=lattice,
            periodic=periodic,
        )

        if self.param is not None:
            param = _damping_param[self.version](**self.param)
        else:
            param = _damping_param[self.version](
                method=self.xc,
                atm=self.atm,
            )

        res = disp.get_dispersion(param=param, grad=True)

        return res.get("energy"), res.get("gradient")

    def reset(self, mol: gto.Mole):
        """Reset mol and clean up relevant attributes for scanner mode"""
        self.mol = mol
        return self


class _DFTD3:
    """
    Stub class used to identify instances of the `DFTD3` class
    """

    pass


class _DFTD3Grad:
    """
    Stub class used to identify instances of the `DFTD3Grad` class
    """

    pass


def d3_energy(mf: scf.hf.SCF, **kwargs) -> scf.hf.SCF:
    """
    Apply DFT-D3 corrections to SCF or MCSCF methods by returning an
    instance of a new class built from the original instances class.
    The dispersion correction is stored in the `with_dftd3` attribute of
    the class.

    Parameters
    ----------
    mf: scf.hf.SCF
        The method to which DFT-D3 corrections will be applied.
    **kwargs
        Keyword arguments passed to the `DFTD3Dispersion` class.

    Returns
    -------
    The method with DFT-D3 corrections applied.

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import dftd3.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          N  -1.57871857  -0.04661102   0.00000000
    ...          N   1.57871857   0.04661102   0.00000000
    ...          H  -2.15862174   0.13639605   0.80956529
    ...          H  -0.84947130   0.65819321   0.00000000
    ...          H  -2.15862174   0.13639605  -0.80956529
    ...          H   2.15862174  -0.13639605  -0.80956529
    ...          H   0.84947130  -0.65819321   0.00000000
    ...          H   2.15862174  -0.13639605   0.80956529
    ...          '''
    ... )
    >>> mf = disp.d3_energy(scf.RHF(mol)).run()
    converged SCF energy = -110.932603617026
    >>> mf.kernel()
    -110.93260361702605
    """

    if not isinstance(mf, (scf.hf.SCF, mcscf.casci.CASCI)):
        raise TypeError("mf must be an instance of SCF or CASCI")

    with_dftd3 = DFTD3Dispersion(
        mf.mol,
        xc="hf"
        if isinstance(mf, mcscf.casci.CASCI)
        else getattr(mf, "xc", "HF").upper().replace(" ", ""),
        **kwargs,
    )

    if isinstance(mf, _DFTD3):
        mf.with_dftd3 = with_dftd3
        return mf

    class DFTD3(_DFTD3, mf.__class__):
        def __init__(self, method, with_dftd3):
            self.__dict__.update(method.__dict__)
            self.with_dftd3 = with_dftd3
            self._keys.update(["with_dftd3"])

        def dump_flags(self, verbose=None):
            mf.__class__.dump_flags(self, verbose)
            if self.with_dftd3:
                self.with_dftd3.dump_flags(verbose)
            return self

        def energy_nuc(self):
            enuc = mf.__class__.energy_nuc(self)
            if self.with_dftd3:
                edisp = self.with_dftd3.kernel()[0]
                mf.scf_summary["dispersion"] = edisp
                enuc += edisp
            return enuc

        def reset(self, mol=None):
            self.with_dftd3.reset(mol)
            return mf.__class__.reset(self, mol)

        def nuc_grad_method(self):
            scf_grad = mf.__class__.nuc_grad_method(self)
            return grad(scf_grad)

        Gradients = lib.alias(nuc_grad_method, alias_name="Gradients")

    return DFTD3(mf, with_dftd3)


def d3_grad(scf_grad: GradientsBase, **kwargs):
    """
    Apply DFT-D3 corrections to SCF or MCSCF nuclear gradients methods
    by returning an instance of a new class built from the original class.
    The dispersion correction is stored in the `with_dftd3` attribute of
    the class.

    Parameters
    ----------
    scf_grad: rhf_grad.Gradients
        The method to which DFT-D3 corrections will be applied.
    **kwargs
        Keyword arguments passed to the `DFTD3Dispersion` class.

    Returns
    -------
    The method with DFT-D3 corrections applied.

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import dftd3.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          O  -1.65542061  -0.12330038   0.00000000
    ...          O   1.24621244   0.10268870   0.00000000
    ...          H  -0.70409026   0.03193167   0.00000000
    ...          H  -2.03867273   0.75372294   0.00000000
    ...          H   1.57598558  -0.38252146  -0.75856129
    ...          H   1.57598558  -0.38252146   0.75856129
    ...          '''
    ... )
    >>> grad = disp.d3_energy(scf.RHF(mol)).run().nuc_grad_method()
    converged SCF energy = -149.947191000075
    >>> g = grad.kernel()
    --------------- DFTD3 gradients ---------------
             x                y                z
    0 O     0.0171886976     0.0506606246     0.0000000000
    1 O     0.0383596853    -0.0459057549     0.0000000000
    2 H    -0.0313133974    -0.0125865676    -0.0000000000
    3 H     0.0066705789    -0.0380501872     0.0000000000
    4 H    -0.0154527822     0.0229409425     0.0215141991
    5 H    -0.0154527822     0.0229409425    -0.0215141991
    ----------------------------------------------
    """

    if not isinstance(scf_grad, GradientsBase):
        raise TypeError("scf_grad must be an instance of Gradients")

    # Ensure that the zeroth order results include DFTD3 corrections
    if not getattr(scf_grad.base, "with_dftd3", None):
        scf_grad.base = energy(scf_grad.base, **kwargs)

    class DFTD3Grad(_DFTD3Grad, scf_grad.__class__):
        def grad_nuc(self, mol=None, atmlst=None):
            nuc_g = scf_grad.__class__.grad_nuc(self, mol, atmlst)
            with_dftd3 = getattr(self.base, "with_dftd3", None)
            if with_dftd3:
                disp_g = with_dftd3.kernel()[1]
                if atmlst is not None:
                    disp_g = disp_g[atmlst]
                nuc_g += disp_g
            return nuc_g

    mfgrad = DFTD3Grad.__new__(DFTD3Grad)
    mfgrad.__dict__.update(scf_grad.__dict__)
    return mfgrad


energy = d3_energy
grad = d3_grad
