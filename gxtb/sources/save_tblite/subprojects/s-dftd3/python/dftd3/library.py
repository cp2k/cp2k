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
Thin wrapper around the CFFI extension. This module mainly acts as a guard
for importing the libdftd3 extension and also provides some FFI based wrappers
for memory handling.
"""

import functools
from typing import Optional

import numpy as np

try:
    from ._libdftd3 import ffi, lib
except ImportError:
    raise ImportError("DFT-D3 C extension unimportable, cannot use C-API")


def get_api_version() -> str:
    """Return the current API version from s-dftd3.
    For easy usage in C the API version is provided as

    10000 * major + 100 * minor + patch

    For Python we want something that looks like a semantic version again.
    """
    api_version = lib.dftd3_get_version()
    return "{}.{}.{}".format(
        api_version // 10000,
        api_version % 10000 // 100,
        api_version % 100,
    )


class Handle:
    def __init__(self, handle):
        self.handle = handle

    @classmethod
    def with_gc(cls, handle):
        return cls(ffi.gc(handle, cls._delete))

    @classmethod
    def null(cls):
        return cls(ffi.NULL)

    @staticmethod
    def _delete(handle):
        raise NotImplementedError("Delete function not implemented")


class StructureHandle(Handle):
    @staticmethod
    def _delete(handle):
        """Delete a DFT-D3 molecular structure data object"""
        ptr = ffi.new("dftd3_structure *")
        ptr[0] = handle
        lib.dftd3_delete_structure(ptr)


class ModelHandle(Handle):
    @staticmethod
    def _delete(handle):
        """Delete a DFT-D3 dispersion model object"""
        ptr = ffi.new("dftd3_model *")
        ptr[0] = handle
        lib.dftd3_delete_model(ptr)


class ParamHandle(Handle):
    @staticmethod
    def _delete(handle):
        """Delete a DFT-D3 damping parameteter object"""
        ptr = ffi.new("dftd3_param *")
        ptr[0] = handle
        lib.dftd3_delete_param(ptr)


class GCPHandle(Handle):
    @staticmethod
    def _delete(handle):
        """Delete a counter-poise parameter object"""
        ptr = ffi.new("dftd3_gcp *")
        ptr[0] = handle
        lib.dftd3_delete_gcp(ptr)


def _delete_error(mol):
    """Delete a DFT-D3 error handle object"""
    ptr = ffi.new("dftd3_error *")
    ptr[0] = mol
    lib.dftd3_delete_error(ptr)


def new_error():
    """Create new DFT-D3 error handle object"""
    return ffi.gc(lib.dftd3_new_error(), _delete_error)


def error_check(func):
    """Handle errors for library functions that require an error handle"""

    @functools.wraps(func)
    def handle_error(*args, **kwargs):
        """Run function and than compare context"""
        _err = new_error()
        value = func(_err, *args, **kwargs)
        if lib.dftd3_check_error(_err):
            _message = ffi.new("char[]", 512)
            lib.dftd3_get_error(_err, _message, ffi.NULL)
            raise RuntimeError(ffi.string(_message).decode())
        return value

    return handle_error


def new_structure(
    natoms: int,
    numbers: np.ndarray,
    positions: np.ndarray,
    lattice: Optional[np.ndarray],
    periodic: Optional[np.ndarray],
) -> StructureHandle:
    """Create new molecular structure data"""
    return StructureHandle.with_gc(
        error_check(lib.dftd3_new_structure)(
            natoms,
            _cast("int*", numbers),
            _cast("double*", positions),
            _cast("double*", lattice),
            _cast("bool*", periodic),
        )
    )


def new_d3_model(mol: StructureHandle) -> ModelHandle:
    """Create new D3 dispersion model"""
    return ModelHandle.with_gc(error_check(lib.dftd3_new_d3_model)(mol.handle))


def set_model_realspace_cutoff(
    disp: ModelHandle, disp2: float, disp3: float, cn: float
) -> None:
    """Set the realspace cutoff for the dispersion model"""
    return error_check(lib.dftd3_set_model_realspace_cutoff)(
        disp.handle, disp2, disp3, cn
    )


def new_zero_damping(
    s6: float, s8: float, s9: float, rs6: float, rs8: float, alp: float
) -> ParamHandle:
    """Create new zero damping parameters"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_new_zero_damping)(s6, s8, s9, rs6, rs8, alp)
    )


def load_zero_damping(method: str, atm: bool) -> ParamHandle:
    """Load zero damping parameters from internal storage"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_load_zero_damping)(_char(method), atm)
    )


def new_rational_damping(
    s6: float, s8: float, s9: float, a1: float, a2: float, alp: float
) -> ParamHandle:
    """Create new rational damping parameters"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_new_rational_damping)(s6, s8, s9, a1, a2, alp)
    )


def load_rational_damping(method: str, atm: bool) -> ParamHandle:
    """Load rational damping parameters from internal storage"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_load_rational_damping)(_char(method), atm)
    )


def new_mzero_damping(
    s6: float, s8: float, s9: float, rs6: float, rs8: float, alp: float, bet: float
) -> ParamHandle:
    """Create new modified zero damping parameters"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_new_mzero_damping)(s6, s8, s9, rs6, rs8, alp, bet)
    )


def load_mzero_damping(method: str, atm: bool) -> ParamHandle:
    """Load modified zero damping parameters from internal storage"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_load_mzero_damping)(_char(method), atm)
    )


def new_mrational_damping(
    s6: float, s8: float, s9: float, a1: float, a2: float, alp: float
) -> ParamHandle:
    """Create new modified rational damping parameters"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_new_mrational_damping)(s6, s8, s9, a1, a2, alp)
    )


def load_mrational_damping(method: str, atm: bool) -> ParamHandle:
    """Load modified rational damping parameters from internal storage"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_load_mrational_damping)(_char(method), atm)
    )


def new_optimizedpower_damping(
    s6: float, s8: float, s9: float, a1: float, a2: float, alp: float, bet: float
) -> ParamHandle:
    """Create new optimized power damping parameters"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_new_optimizedpower_damping)(s6, s8, s9, a1, a2, alp, bet)
    )


def load_optimizedpower_damping(method: str, atm: bool) -> ParamHandle:
    """Load optimized power damping parameters from internal storage"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_load_optimizedpower_damping)(_char(method), atm)
    )


def new_cso_damping(
    s6: float, s9: float, a1: float, a2: float, a3: float, a4: float, alp: float
) -> ParamHandle:
    """Create new CSO (C6-scaled only) damping parameters"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_new_cso_damping)(s6, s9, a1, a2, a3, a4, alp)
    )


def load_cso_damping(method: str, atm: bool) -> ParamHandle:
    """Load CSO damping parameters from internal storage"""
    return ParamHandle.with_gc(
        error_check(lib.dftd3_load_cso_damping)(_char(method), atm)
    )


def update_structure(
    mol: StructureHandle, positions: np.ndarray, lattice: Optional[np.ndarray]
) -> None:
    """Update the molecular structure data"""
    return error_check(lib.dftd3_update_structure)(
        mol.handle,
        _cast("double*", positions),
        _cast("double*", lattice),
    )


def get_dispersion(
    mol: StructureHandle,
    disp: ModelHandle,
    param: ParamHandle,
    energy: np.ndarray,
    gradient: Optional[np.ndarray],
    sigma: Optional[np.ndarray],
) -> None:
    return error_check(lib.dftd3_get_dispersion)(
        mol.handle,
        disp.handle,
        param.handle,
        _cast("double*", energy),
        _cast("double*", gradient),
        _cast("double*", sigma),
    )


def get_pairwise_dispersion(
    mol: StructureHandle,
    disp: ModelHandle,
    param: ParamHandle,
    energy2: np.ndarray,
    energy3: np.ndarray,
) -> None:
    return error_check(lib.dftd3_get_pairwise_dispersion)(
        mol.handle,
        disp.handle,
        param.handle,
        _cast("double*", energy2),
        _cast("double*", energy3),
    )


def load_gcp_param(
    mol: StructureHandle, method: Optional[str], basis: Optional[str]
) -> GCPHandle:
    """Load GCP parameters from internal storage"""
    return GCPHandle.with_gc(
        error_check(lib.dftd3_load_gcp_param)(mol.handle, _char(method), _char(basis))
    )


def set_gcp_realspace_cutoff(gcp: GCPHandle, bas: float, srb: float) -> None:
    error_check(lib.dftd3_set_gcp_realspace_cutoff)(gcp.handle, bas, srb)


def get_counterpoise(
    mol: StructureHandle,
    gcp: GCPHandle,
    energy: np.ndarray,
    gradient: Optional[np.ndarray],
    sigma: Optional[np.ndarray],
) -> None:
    """Get the counterpoise energy"""
    return error_check(lib.dftd3_get_counterpoise)(
        mol.handle,
        gcp.handle,
        _cast("double*", energy),
        _cast("double*", gradient),
        _cast("double*", sigma),
    )


def _char(value: Optional[str]):
    """Convert a string to a C char array"""
    return ffi.new("char[]", value.encode()) if value is not None else ffi.NULL


def _cast(ctype: str, array: Optional[np.ndarray]):
    """Cast a numpy array to a FFI pointer"""
    return ffi.cast(ctype, array.ctypes.data) if array is not None else ffi.NULL
