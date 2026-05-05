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
PySCF Support
-------------

Compatibility layer for supporting DFT-D4 in `pyscf <https://pyscf.org/>`_.
"""

try:
    from pyscf import lib, gto, mcscf, scf
    from pyscf.grad import rhf as rhf_grad
except ModuleNotFoundError:
    raise ModuleNotFoundError("This submodule requires pyscf installed")

from typing import Tuple, Optional, Dict

import numpy as np

from .interface import DampingFunction, DampingParam, DispersionModel

GradientsBase = getattr(rhf_grad, "GradientsBase", rhf_grad.Gradients)


class DFTD4Dispersion(lib.StreamObject):
    """
    Implementation of the interface for using DFT-D4 in pyscf.
    The `xc` functional can be provided in the constructor together with the
    DFT-D `model` and optionally the two-body (`damping_2b`) and three-body
    (`damping_3b`) damping functions to use.
    Possible two-body damping functions are

    ``"rational"``: (default for D4 and D4S)
        For rational (Becke-Johnson) damping function
    ``"screened"``:
        For screened rational damping function
    ``"zero"``
        For zero (Chai-Head-Gordon) damping function
    ``"mzero"``
        Modified version of the zero damping function
    ``"optpower"``
        Optimized power damping function
    ``"cso"``
        CSO (C6-scaled only) damping function
    ``"koide"``
        Koide damping function

    Possible three-body damping functions are

    ``"rational"``: (default for D4 and D4S)
        For rational (Becke-Johnson) damping function
    ``"screened"``:
        For screened rational damping function
    ``"zero"``
        For zero (Chai-Head-Gordon) damping function
    ``"zero-avg"``
        For averaged distance zero damping function

    Custom parameters can be provided with the `param` dictionary.
    The `param` dict contains the damping parameters, at least s8, a1 and a2
    must be provided for the default (rational + zero-avg) damping of the d4
    model. For all other damping functions all parameter must be specified
    and optionally set zero if not used. The parameters are:

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

    Defaults are given for the general case where all parameters are required and for 
    the default damping of the d4 model (second value).
    The version of the damping can be changed after constructing the dispersion correction.
    With the `atm` boolean the three-body dispersion energy can be disabled.    

    Examples
    --------
    >>> from pyscf import gto
    >>> import dftd4.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          C   -0.755422531  -0.796459123  -1.023590391
    ...          C    0.634274834  -0.880017014  -1.075233285
    ...          C    1.406955202   0.199695367  -0.653144334
    ...          C    0.798863737   1.361204515  -0.180597909
    ...          C   -0.593166787   1.434312023  -0.133597923
    ...          C   -1.376239198   0.359205222  -0.553258516
    ...          I   -1.514344238   3.173268101   0.573601106
    ...          H    1.110906949  -1.778801728  -1.440619836
    ...          H    1.399172302   2.197767355   0.147412751
    ...          H    2.486417780   0.142466525  -0.689380574
    ...          H   -2.454252250   0.422581120  -0.512807958
    ...          H   -1.362353593  -1.630564523  -1.348743149
    ...          S   -3.112683203   6.289227834   1.226984439
    ...          H   -4.328789697   5.797771251   0.973373089
    ...          C   -2.689135032   6.703163830  -0.489062886
    ...          H   -1.684433029   7.115457372  -0.460265708
    ...          H   -2.683867206   5.816530502  -1.115183775
    ...          H   -3.365330613   7.451201412  -0.890098894
    ...          '''
    ... )
    >>> d4 = disp.DFTD4Dispersion(mol, xc="r2SCAN")
    >>> d4.kernel()[0]
    array(-0.0050011)
    >>> d4.damping = {"3b":"none"}
    >>> d4.kernel()[0]
    array(-0.00499889)
    >>> d4.atm = False
    >>> d4.kernel()[0]
    array(-0.00499638)
    """

    def __init__(self,
                 mol: gto.Mole,
                 xc: str = "hf",
                 atm: bool = True,
                 model: str = "d4",
                 damping: Optional[Dict[str, str]] = None,
                 param: Optional[Dict[str, float]] = None):
        self.mol = mol
        self.verbose = mol.verbose
        self.xc = xc
        self.atm = atm
        self.model = model
        self.damping = damping
        self.param = param

    def dump_flags(self, verbose=None) -> "DFTD4Dispersion":
        """
        Show options used for the DFT-D4 dispersion correction.
        """
        lib.logger.info(self, "** DFTD4 parameter **")
        lib.logger.info(self, "func %s", self.xc)
        lib.logger.info(self, "model %s", self.model)

        damping = self.damping or {}
        _d2 = damping.get("2b", "<default>")
        _d3 = damping.get("3b", "<default>") if self.atm else "none"
        lib.logger.info(self, "damping_2b %s", _d2)
        lib.logger.info(self, "damping_3b %s", _d3 if self.atm else "none")
        return self

    def kernel(self) -> Tuple[float, np.ndarray]:
        """
        Compute the DFT-D4 dispersion correction.

        The dispersion model as well as the parameters are created locally and
        not part of the state of the instance.

        Returns
        -------
        float, ndarray
            The energy and gradient of the DFT-D4 dispersion correction.
        
        Examples
        --------
        >>> from pyscf import gto
        >>> import dftd4.pyscf as disp
        >>> mol = gto.M(
        ...     atom='''
        ...          Br    0.000000    0.000000    1.919978
        ...          Br    0.000000    0.000000   -0.367147
        ...          N     0.000000    0.000000   -3.235006
        ...          C     0.000000    0.000000   -4.376626
        ...          H     0.000000    0.000000   -5.444276
        ...          '''
        ... )
        >>> d4 = disp.DFTD4Dispersion(mol, xc="PBE0", model="d4")
        >>> energy, gradient = d4.kernel()
        >>> energy
        array(-0.00296818)
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
        if hasattr(mol, 'lattice_vectors'):
            lattice = mol.lattice_vectors()
            periodic = np.array([True, True, True], dtype=bool)

        disp = DispersionModel(
            np.asarray([gto.charge(sym) for sym in mol.elements]),
            mol.atom_coords(),
            mol.charge,
            lattice=lattice,
            periodic=periodic,
            model=self.model,
        )

        # Select if two- and three-body damping functions are specified
        if self.damping is not None:
            _d2 = self.damping.get("2b")
            _d3 = self.damping.get("3b") if self.atm else "none"
        else:
            _d2 = None
            _d3 = None if self.atm else "none"

        if _d2 is not None:
            # Explicitly specified damping functions
            damp = DampingFunction(damping_2b=_d2, damping_3b=_d3)
        else:
            if _d3 is None:
                # Default damping functions with unmodified ATM
                damp = DampingFunction(model=self.model)
            else:
                # Default damping functions with modified ATM
                damp = DampingFunction(model=self.model, damping_3b=_d3)

        # Parameter selection
        if self.param is not None:
            if all (key in self.param for key in ("s6", "s8", "s9",
                                                  "a1", "a2", "a3", "a4",
                                                  "rs6", "rs8", "rs9",
                                                  "alp", "bet")):
                # Explicitly specify all possible parameters
                param = DampingParam(**self.param)
            else:
                # Explicitly specified parameters for model default damping
                param = DampingParam(model=self.model, **self.param)
        elif _d2 is not None:
            # Method specific parameters for the model with modified damping
            param = DampingParam(method=self.xc, model=self.model,
                                 damping_2b=_d2, damping_3b=_d3)
        else:
            if _d3 is None:
                # Method and model specific parameters with default damping
                # and unmodified ATM
                param = DampingParam(method=self.xc, model=self.model)
            else:
                # Method and model specific parameters with default damping
                # and modified ATM
                param = DampingParam(method=self.xc, model=self.model,
                                     damping_3b=_d3)
        
        # Actual dispersion calculation
        res = disp.get_dispersion(damp=damp, param=param, grad=True)

        return res.get("energy"), res.get("gradient")

    def reset(self, mol: gto.Mole) -> "DFTD4Dispersion":
        """
        Reset mol and clean up relevant attributes for scanner mode
        """
        self.mol = mol
        return self


class _DFTD4:
    """
    Stub class used to identify instances of the `DFTD4` class
    """

    pass


class _DFTD4Grad:
    """
    Stub class used to identify instances of the `DFTD4Grad` class
    """

    pass


def energy(mf: scf.hf.SCF, **kwargs) -> scf.hf.SCF:
    """
    Apply DFT-D4 corrections to SCF or MCSCF methods by returning an
    instance of a new class built from the original instances class.
    The dispersion correction is stored in the `with_dftd4` attribute of
    the class.

    Parameters
    ----------
    mf: scf.hf.SCF
        The method to which DFT-D4 corrections will be applied.
    **kwargs
        Keyword arguments passed to the `DFTD4Dispersion` class.

    Returns
    -------
    The method with DFT-D4 corrections applied.

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import dftd4.pyscf as disp
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
    >>> mf = disp.energy(scf.RHF(mol)).run()
    converged SCF energy = -110.917424528592
    >>> mf.kernel()
    -110.917424528592
    """

    if not isinstance(mf, (scf.hf.SCF, mcscf.casci.CASCI)):
        raise TypeError("mf must be an instance of SCF or CASCI")

    with_dftd4 = DFTD4Dispersion(
        mf.mol,
        xc="hf"
        if isinstance(mf, mcscf.casci.CASCI)
        else getattr(mf, "xc", "HF").upper().replace(" ", ""),
        **kwargs,
    )

    if isinstance(mf, _DFTD4):
        mf.with_dftd4 = with_dftd4
        return mf

    class DFTD4(_DFTD4, mf.__class__):
        """
        Patched SCF class including DFT-D4 corrections.
        """

        def __init__(self, method, with_dftd4: DFTD4Dispersion):
            self.__dict__.update(method.__dict__)
            self.with_dftd4 = with_dftd4
            self._keys.update(["with_dftd4"])

        def dump_flags(self, verbose=None) -> "DFTD4":
            mf.__class__.dump_flags(self, verbose)
            if self.with_dftd4:
                self.with_dftd4.dump_flags(verbose)
            return self

        def energy_nuc(self) -> float:
            enuc = mf.__class__.energy_nuc(self)
            if self.with_dftd4:
                edisp = self.with_dftd4.kernel()[0]
                self.scf_summary["dispersion"] = edisp
                enuc += edisp
            return enuc

        def reset(self, mol=None) -> "DFTD4":
            self.with_dftd4.reset(mol)
            return mf.__class__.reset(self, mol)

        def nuc_grad_method(self) -> "DFTD4Grad":
            scf_grad = mf.__class__.nuc_grad_method(self)
            return grad(scf_grad)

        Gradients = lib.alias(nuc_grad_method, alias_name="Gradients")

    return DFTD4(mf, with_dftd4)


def grad(mfgrad: GradientsBase, **kwargs):
    """
    Apply DFT-D4 corrections to SCF or MCSCF nuclear gradients methods
    by returning an instance of a new class built from the original class.
    The dispersion correction is stored in the `with_dftd4` attribute of
    the class.

    Parameters
    ----------
    mfgrad
        The method to which DFT-D4 corrections will be applied.
    **kwargs
        Keyword arguments passed to the `DFTD4Dispersion` class.

    Returns
    -------
    The method with DFT-D4 corrections applied.

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import dftd4.pyscf as disp
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
    >>> grad = disp.energy(scf.RHF(mol)).run().nuc_grad_method()
    converged SCF energy = -149.939098424774
    >>> g = grad.kernel()
    --------------- DFTD4 gradients ---------------
             x                y                z
    0 O     0.0172438133     0.0508406920     0.0000000000
    1 O     0.0380018285    -0.0460223790    -0.0000000000
    2 H    -0.0305058266    -0.0126478132    -0.0000000000
    3 H     0.0069233858    -0.0382898692    -0.0000000000
    4 H    -0.0158316004     0.0230596847     0.0218908543
    5 H    -0.0158316004     0.0230596847    -0.0218908543
    ----------------------------------------------
    """

    if not isinstance(mfgrad, GradientsBase):
        raise TypeError("mfgrad must be an instance of Gradients")

    # Ensure that the zeroth order results include DFTD4 corrections
    if not getattr(mfgrad.base, "with_dftd4", None):
        mfgrad.base = energy(mfgrad.base, **kwargs)

    class DFTD4Grad(_DFTD4Grad, mfgrad.__class__):
        """
        Patched SCF class including DFT-D4 corrections.
        """

        def grad_nuc(self, mol=None, atmlst=None):
            nuc_g = mfgrad.__class__.grad_nuc(self, mol, atmlst)
            with_dftd4 = getattr(self.base, "with_dftd4", None)
            if with_dftd4:
                disp_g = with_dftd4.kernel()[1]
                if atmlst is not None:
                    disp_g = disp_g[atmlst]
                nuc_g += disp_g
            return nuc_g

    dgrad = DFTD4Grad.__new__(DFTD4Grad)
    dgrad.__dict__.update(mfgrad.__dict__)
    return dgrad
