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
Thin wrapper around CFFI extension module tblite.

This module mainly acts as a guard for importing the libtblite extension and
also provides some FFI based wappers for memory handling.
"""

import functools
import sys
from typing import Callable, Dict

import numpy as np

try:
    from ._libtblite import ffi, lib
except ImportError as e:
    raise ImportError(
        "tblite C extension unimportable, cannot use C-API"
    ) from e

from .exceptions import TBLiteRuntimeError, TBLiteValueError


def get_version() -> tuple:
    """Return the current API version from tblite.
    For easy usage in C the API version is provided as

    10000 * major + 100 * minor + patch

    For Python we want something that looks like a semantic version again.
    """
    version = lib.tblite_get_version()
    return (
        version // 10000,
        version % 10000 // 100,
        version % 100,
    )


def _delete_error(error) -> None:
    """Delete a tblite error handler object"""
    ptr = ffi.new("tblite_error *")
    ptr[0] = error
    lib.tblite_delete_error(ptr)


def new_error():
    """Create new tblite error handler object"""
    return ffi.gc(lib.tblite_new_error(), _delete_error)


def error_check(func):
    """Handle errors for library functions that require an error handle"""

    @functools.wraps(func)
    def handle_error(*args, **kwargs):
        """Run function and than compare context"""
        _err = new_error()
        value = func(_err, *args, **kwargs)
        if lib.tblite_check_error(_err):
            _message = ffi.new("char[]", 512)
            lib.tblite_get_error(_err, _message, ffi.NULL)
            raise TBLiteRuntimeError(ffi.string(_message).decode())
        return value

    return handle_error


@ffi.def_extern()
def logger_callback(error, message, nchar, data):
    """Custom logger callback to write output in a Python friendly way"""
    try:
        callback = ffi.from_handle(data)
        callback(ffi.unpack(message, nchar).decode())
    except Exception as e:
        error_message = ffi.new("char[]", str(e).encode())
        lib.tblite_set_error(
            error,
            error_message,
            ffi.NULL,
        )


def context_check(func):
    """Handle errors for library functions that require a context handle"""

    @functools.wraps(func)
    def handle_context_error(ctx, *args, **kwargs):
        """Run function and than compare context"""
        if isinstance(ctx, tuple):
            ctx, _ = ctx
        value = func(ctx, *args, **kwargs)
        if lib.tblite_check_context(ctx):
            _message = ffi.new("char[]", 512)
            lib.tblite_get_context_error(ctx, _message, ffi.NULL)
            raise TBLiteRuntimeError(ffi.string(_message).decode())
        return value

    return handle_context_error


def _delete_context(context) -> None:
    """Delete a tblite context handler object"""
    ptr = ffi.new("tblite_context *")
    ptr[0] = context
    lib.tblite_delete_context(ptr)


def new_context(color: bool = True, logger: Callable[[str], None] = print):
    """Create new tblite context handler object"""
    ctx = ffi.gc(lib.tblite_new_context(), _delete_context)
    handle = ffi.new_handle(logger)
    context_check(lib.tblite_set_context_logger)(
        ctx, lib.logger_callback, handle
    )
    if color and sys.stdout.isatty():
        context_check(lib.tblite_set_context_color)(ctx, 1)
    return ctx, handle


def _delete_structure(mol) -> None:
    """Delete molecular structure data"""
    ptr = ffi.new("tblite_structure *")
    ptr[0] = mol
    lib.tblite_delete_structure(ptr)


def new_structure(natoms, numbers, positions, charge, uhf, lattice, periodic):
    """Create new molecular structure data"""
    return ffi.gc(
        error_check(lib.tblite_new_structure)(
            natoms,
            numbers,
            positions,
            charge,
            uhf,
            lattice,
            periodic,
        ),
        _delete_structure,
    )


update_structure_geometry = error_check(lib.tblite_update_structure_geometry)
update_structure_charge = error_check(lib.tblite_update_structure_charge)
update_structure_uhf = error_check(lib.tblite_update_structure_uhf)


def _delete_table(table) -> None:
    """Delete a tblite data table object"""
    ptr = ffi.new("tblite_table *")
    ptr[0] = table
    lib.tblite_delete_table(ptr)


def new_table():
    """Create a tblite data table object"""
    return ffi.gc(lib.tblite_new_table(ffi.NULL), _delete_table)


table_set_double = error_check(lib.tblite_table_set_double)
table_set_int64_t = error_check(lib.tblite_table_set_int64_t)
table_set_bool = error_check(lib.tblite_table_set_bool)
table_set_char = error_check(lib.tblite_table_set_char)
table_add_table = error_check(lib.tblite_table_add_table)


def _delete_param(param) -> None:
    """Delete a tblite parametrization record object"""
    ptr = ffi.new("tblite_param *")
    ptr[0] = param
    lib.tblite_delete_param(ptr)


def new_param():
    """Create a tblite data table object"""
    return ffi.gc(lib.tblite_new_param(), _delete_param)


load_param = error_check(lib.tblite_load_param)
dump_param = error_check(lib.tblite_dump_param)
export_gfn2_param = error_check(lib.tblite_export_gfn2_param)
export_gfn1_param = error_check(lib.tblite_export_gfn1_param)
export_ipea1_param = error_check(lib.tblite_export_ipea1_param)


def _delete_result(result) -> None:
    """Delete a tblite result container object"""
    ptr = ffi.new("tblite_result *")
    ptr[0] = result
    lib.tblite_delete_result(ptr)


def new_result():
    """Create new tblite result container object"""
    return ffi.gc(lib.tblite_new_result(), _delete_result)


def copy_result(res):
    """Create new tblite result container object as copy of an existing result"""
    return ffi.gc(lib.tblite_copy_result(res), _delete_result)


def get_number_of_atoms(res) -> int:
    """Retrieve number of atoms from result container"""
    _natoms = ffi.new("int *")
    error_check(lib.tblite_get_result_number_of_atoms)(res, _natoms)
    return _natoms[0]


def get_number_of_orbitals(res) -> int:
    """Retrieve number of orbitals from result container"""
    _norb = ffi.new("int *")
    error_check(lib.tblite_get_result_number_of_orbitals)(res, _norb)
    return _norb[0]


def get_number_of_spins(res) -> int:
    """Retrieve number of spins from result container"""
    _nspin = ffi.new("int *")
    error_check(lib.tblite_get_result_number_of_spins)(res, _nspin)
    return _nspin[0]


def get_energy(res) -> np.ndarray:
    """Retrieve energy from result container"""
    _energy = np.array(0.0)
    error_check(lib.tblite_get_result_energy)(
        res, ffi.cast("double*", _energy.ctypes.data)
    )
    return _energy


def get_energies(res) -> np.ndarray:
    """Retrieve atom-resolved energies from result container"""
    _energies = np.zeros((get_number_of_atoms(res),))
    error_check(lib.tblite_get_result_energies)(
        res, ffi.cast("double*", _energies.ctypes.data)
    )
    return _energies


def get_gradient(res) -> np.ndarray:
    """Retrieve gradient from result container"""
    _gradient = np.zeros((get_number_of_atoms(res), 3))
    error_check(lib.tblite_get_result_gradient)(
        res, ffi.cast("double*", _gradient.ctypes.data)
    )
    return _gradient


def get_virial(res) -> np.ndarray:
    """Retrieve virial from result container"""
    _virial = np.zeros((3, 3))
    error_check(lib.tblite_get_result_virial)(
        res, ffi.cast("double*", _virial.ctypes.data)
    )
    return _virial


def get_charges(res) -> np.ndarray:
    """Retrieve atomic charges from result container"""
    _charges = np.zeros((get_number_of_atoms(res),))
    error_check(lib.tblite_get_result_charges)(
        res, ffi.cast("double*", _charges.ctypes.data)
    )
    return _charges


def get_bond_orders(res) -> np.ndarray:
    """Retrieve Wiberg / Mayer bond orders from result container"""
    _dict = get_post_processing_dict(res=res)
    if "bond-orders" not in _dict:
        raise TBLiteValueError(
            "Bond-orders were not calculated. By default they are computed."
        )
    _bond_orders = _dict["bond-orders"]
    return _bond_orders


def get_dipole(res) -> np.ndarray:
    """Retrieve dipole moment from post processing dictionary"""
    _dict = get_post_processing_dict(res=res)
    if "molecular-dipole" not in _dict:
        raise TBLiteValueError(
            "Molecular dipole was not calculated. By default it is computed."
        )
    _dipole = _dict["molecular-dipole"]
    return _dipole


def get_quadrupole(res) -> np.ndarray:
    """Retrieve quadrupole moment from post processing dictionary"""
    _dict = get_post_processing_dict(res=res)
    if "molecular-quadrupole" not in _dict:
        raise TBLiteValueError(
            "Molecular quadrupole was not calculated. By default it is computed."
        )
    _quadrupole = _dict["molecular-quadrupole"]
    return _quadrupole


def get_orbital_energies(res) -> np.ndarray:
    """Retrieve orbital energies from result container"""
    _norb = get_number_of_orbitals(res)
    _nspin = get_number_of_spins(res)
    _emo = np.zeros((_nspin, _norb))
    error_check(lib.tblite_get_result_orbital_energies)(
        res, ffi.cast("double*", _emo.ctypes.data)
    )
    if _nspin == 1:
        return np.squeeze(_emo, axis=0)
    return _emo


def get_orbital_occupations(res) -> np.ndarray:
    """Retrieve orbital occupations from result container"""
    _norb = get_number_of_orbitals(res)
    _nspin = get_number_of_spins(res)
    _occ = np.zeros((2, _norb))
    error_check(lib.tblite_get_result_orbital_occupations)(
        res, ffi.cast("double*", _occ.ctypes.data)
    )
    if _nspin == 1:
        return np.sum(_occ, axis=0)
    return _occ


def load_wavefunction(res, filename: str) -> None:
    _filename = ffi.new("char[]", filename.encode("ascii"))
    error_check(lib.tblite_load_result_wavefunction)(res, _filename)


def save_wavefunction(res, filename: str) -> None:
    _filename = ffi.new("char[]", filename.encode("ascii"))
    error_check(lib.tblite_save_result_wavefunction)(res, _filename)


def _get_ao_matrix(getter, is_spin_dependent: bool):
    """Correctly set allocation for matrix objects before querying the getter"""

    @functools.wraps(getter)
    def with_allocation(res):
        """Get a matrix property from the results object"""
        _norb = get_number_of_orbitals(res)
        _nspin = get_number_of_spins(res) if is_spin_dependent else 1

        # (_norb, _norb, _nspin) in col-major -> (_nspin, _norb, _norb) in row-major
        # this will allow us to extract alpha- and beta matrices as mat[0] and mat[1]
        _mat = np.zeros((_nspin, _norb, _norb))
        error_check(getter)(res, ffi.cast("double*", _mat.ctypes.data))

        # Transpose actual matrix from col-major to row-major
        # -> important for orbital coefficients
        _mat = np.swapaxes(_mat, 1, 2)

        if _nspin == 1:
            return np.squeeze(_mat, axis=0)
        return _mat

    return with_allocation


get_orbital_coefficients = _get_ao_matrix(
    lib.tblite_get_result_orbital_coefficients, True
)
get_density_matrix = _get_ao_matrix(lib.tblite_get_result_density_matrix, True)
get_overlap_matrix = _get_ao_matrix(lib.tblite_get_result_overlap_matrix, False)
get_hamiltonian_matrix = _get_ao_matrix(
    lib.tblite_get_result_hamiltonian_matrix, False
)


def _delete_calculator(calc) -> None:
    """Delete a tblite calculator object"""
    ptr = ffi.new("tblite_calculator *")
    ptr[0] = calc
    lib.tblite_delete_calculator(ptr)


@context_check
def new_gfn2_calculator(ctx, mol):
    """Create new tblite calculator loaded with GFN2-xTB parametrization data"""
    return ffi.gc(lib.tblite_new_gfn2_calculator(ctx, mol), _delete_calculator)


@context_check
def new_gfn1_calculator(ctx, mol):
    """Create new tblite calculator loaded with GFN1-xTB parametrization data"""
    return ffi.gc(lib.tblite_new_gfn1_calculator(ctx, mol), _delete_calculator)


@context_check
def new_ipea1_calculator(ctx, mol):
    """Create new tblite calculator loaded with IPEA1-xTB parametrization data"""
    return ffi.gc(lib.tblite_new_ipea1_calculator(ctx, mol), _delete_calculator)


@context_check
def new_xtb_calculator(ctx, mol, param):
    """Create new tblite calculator from parametrization records"""
    return ffi.gc(
        lib.tblite_new_xtb_calculator(ctx, mol, param), _delete_calculator
    )


def get_calculator_shell_map(ctx, calc) -> np.ndarray:
    """Retrieve index mapping from shells to atomic centers"""
    _nsh = ffi.new("int *")
    context_check(lib.tblite_get_calculator_shell_count)(ctx, calc, _nsh)
    _map = np.zeros((_nsh[0],), dtype=np.int32)
    context_check(lib.tblite_get_calculator_shell_map)(
        ctx, calc, ffi.cast("int*", _map.ctypes.data)
    )
    return _map


def _delete_double_dictionary(calc) -> None:
    """Delete a tblite double dictionary object"""
    ptr = ffi.new("tblite_double_dictionary *")
    ptr[0] = calc
    lib.tblite_delete_double_dictionary(ptr)


def get_post_processing_dict(res) -> Dict[str, np.ndarray]:
    """Retrieve the dictionary containing all post processing results"""
    _dict = ffi.gc(
        error_check(lib.tblite_get_post_processing_dict)(res),
        _delete_double_dictionary,
    )
    _nentries = error_check(lib.tblite_get_n_entries_dict)(_dict)
    _dict_py = {}
    for i in range(1, _nentries + 1):
        _index = ffi.new("const int*", i)

        _dim1 = ffi.new("int*")
        _dim2 = ffi.new("int*")
        _dim3 = ffi.new("int*")
        error_check(lib.tblite_get_array_size_index)(
            _dict, _index, _dim1, _dim2, _dim3
        )
        if _dim3[0] == 0:
            if _dim2[0] == 0:
                _array = np.zeros((_dim1[0],))
            else:
                _array = np.zeros((_dim1[0], _dim2[0]))
        else:
            _array = np.zeros((_dim1[0], _dim2[0], _dim3[0]))

        error_check(lib.tblite_get_array_entry_index)(
            _dict, _index, ffi.cast("double*", _array.ctypes.data)
        )
        _message = ffi.new("char[]", 512)
        error_check(lib.tblite_get_label_entry_index)(
            _dict, _index, _message, ffi.NULL
        )
        label = ffi.string(_message).decode()
        _dict_py[label] = _array
    return _dict_py


def get_calculator_angular_momenta(ctx, calc) -> np.ndarray:
    """Retrieve angular momenta of shells"""
    _nsh = ffi.new("int *")
    context_check(lib.tblite_get_calculator_shell_count)(ctx, calc, _nsh)
    _am = np.zeros((_nsh[0],), dtype=np.int32)
    context_check(lib.tblite_get_calculator_angular_momenta)(
        ctx, calc, ffi.cast("int*", _am.ctypes.data)
    )
    return _am


def get_calculator_orbital_map(ctx, calc) -> np.ndarray:
    """Retrieve index mapping from atomic orbitals to shells"""
    _nao = ffi.new("int *")
    context_check(lib.tblite_get_calculator_orbital_count)(ctx, calc, _nao)
    _map = np.zeros((_nao[0],), dtype=np.int32)
    context_check(lib.tblite_get_calculator_orbital_map)(
        ctx, calc, ffi.cast("int*", _map.ctypes.data)
    )
    return _map


set_calculator_max_iter = context_check(lib.tblite_set_calculator_max_iter)
set_calculator_accuracy = context_check(lib.tblite_set_calculator_accuracy)
set_calculator_mixer_damping = context_check(
    lib.tblite_set_calculator_mixer_damping
)
set_calculator_guess = context_check(lib.tblite_set_calculator_guess)
set_calculator_temperature = context_check(
    lib.tblite_set_calculator_temperature
)
set_calculator_save_integrals = context_check(
    lib.tblite_set_calculator_save_integrals
)


@context_check
def post_processing_push_back(ctx, calc, s):
    _string = ffi.new("char[]", s.encode("ascii"))
    lib.tblite_push_back_post_processing_str(ctx, calc, _string)


@context_check
def set_calculator_verbosity(ctx, calc, verbosity: int):
    """Set verbosity in context associated with calculator"""
    lib.tblite_set_context_verbosity(ctx, verbosity)


get_singlepoint = context_check(lib.tblite_get_singlepoint)


def _delete_container(cont) -> None:
    """Delete a tblite container object"""
    ptr = ffi.new("tblite_container *")
    ptr[0] = cont
    lib.tblite_delete_container(ptr)


def new_electric_field(ctx, mol, calc, efield):
    """Create new tblite electric field object"""
    return lib.tblite_new_electric_field(efield)


_state_enum = {
    "gsolv": 1,
    "bar1mol": 2,
    "reference": 3,
}

_born_enum = {
    "still": 1,
    "p16": 2,
}


def new_alpb_solvation(
    ctx, mol, calc, solvent: str, state: str = "gsolv", *, version: int
):
    """Create new tblite ALPB solvation object"""
    _solvent = ffi.new("char[]", solvent.encode("ascii"))
    _version = 10 + version
    return error_check(lib.tblite_new_alpb_solvation_solvent)(
        mol, _solvent, _version, _state_enum[state]
    )


def new_gbsa_solvation(
    ctx, mol, calc, solvent: str, state: str = "gsolv", *, version: int
):
    """Create new tblite GBSA solvation object"""
    _solvent = ffi.new("char[]", solvent.encode("ascii"))
    _version = 20 + version
    return error_check(lib.tblite_new_alpb_solvation_solvent)(
        mol, _solvent, _version, _state_enum[state]
    )


def new_gbe_solvation(ctx, mol, calc, epsilon: float, born: str):
    """Create new tblite GBE solvation object"""
    return error_check(lib.tblite_new_gb_solvation_epsilon)(
        mol, float(epsilon), 10, _born_enum[born]
    )


def new_gb_solvation(ctx, mol, calc, epsilon: float, born: str):
    """Create new tblite GB solvation object"""
    return error_check(lib.tblite_new_gb_solvation_epsilon)(
        mol, float(epsilon), 20, _born_enum[born]
    )


def new_cpcm_solvation(ctx, mol, calc, epsilon: float):
    """Create new tblite CPCM solvation object"""
    return error_check(lib.tblite_new_cpcm_solvation_epsilon)(
        mol, float(epsilon)
    )


@context_check
def new_spin_polarization(ctx, mol, calc, wscale: float = 1.0):
    """Create new tblite spin polarization object"""
    return lib.tblite_new_spin_polarization(ctx, mol, calc, wscale)


@context_check
def calculator_push_back(ctx, calc, cont) -> None:
    """Add container to calculator object."""
    ptr = ffi.new("tblite_container *")
    ptr[0] = cont
    lib.tblite_calculator_push_back(ctx, calc, ptr)
