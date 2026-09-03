"""ctypes access to the CP2K AICS process symbols and their array views."""

import ctypes
import math
import os
from typing import Iterable, Tuple

import numpy as np

AICS_ARRAY_SYMBOLS = (
    "cp2k_aics_get_electrostatic_potential",
    "cp2k_aics_get_grid_metadata",
    "cp2k_aics_get_relative_permittivity_parallel",
    "cp2k_aics_get_relative_permittivity_perpendicular",
    "cp2k_aics_get_solute_charge_density",
)

AICS_FD_SENSITIVITY_SYMBOLS = (
    "cp2k_aics_get_fd_field_squared_sensitivity_parallel",
    "cp2k_aics_get_fd_field_squared_sensitivity_perpendicular",
)


_GRID_INFO_ARGUMENTS = (
    *([ctypes.POINTER(ctypes.c_size_t)] * 12),
    *([ctypes.POINTER(ctypes.c_double)] * 3),
)
_SYMBOL_ABI = {
    "cp2k_aics_get_grid_metadata": (_GRID_INFO_ARGUMENTS, None),
    "cp2k_aics_get_solute_charge_density": ((), ctypes.c_void_p),
    "cp2k_aics_get_relative_permittivity_parallel": ((), ctypes.c_void_p),
    "cp2k_aics_get_relative_permittivity_perpendicular": ((), ctypes.c_void_p),
    "cp2k_aics_get_electrostatic_potential": ((), ctypes.c_void_p),
    "cp2k_aics_get_fd_field_squared_sensitivity_parallel": (
        (),
        ctypes.c_void_p,
    ),
    "cp2k_aics_get_fd_field_squared_sensitivity_perpendicular": (
        (),
        ctypes.c_void_p,
    ),
}


def load_aics_library(library_path, required_symbols: Iterable[str]):
    """Bind AICS symbols from the exact already-loaded CP2K library."""
    if not isinstance(library_path, (str, bytes, os.PathLike)):
        raise RuntimeError("Invalid CP2K AICS bridge-library pathname")
    try:
        library = ctypes.CDLL(
            library_path,
            mode=os.RTLD_NOW | os.RTLD_LOCAL | os.RTLD_NOLOAD,
        )
    except OSError as exc:
        raise RuntimeError(
            "CP2K AICS bridge library is not the expected loaded object"
        ) from exc
    for name in required_symbols:
        try:
            symbol = getattr(library, name)
        except AttributeError as exc:
            raise RuntimeError(
                f"Missing required CP2K AICS bridge symbol: {name}"
            ) from exc
        try:
            argtypes, restype = _SYMBOL_ABI[name]
        except KeyError as exc:
            raise RuntimeError(
                f"No ctypes ABI declaration for AICS symbol: {name}"
            ) from exc
        symbol.argtypes = list(argtypes)
        symbol.restype = restype
    return library


def configure_grid_info(library) -> Tuple[ctypes._SimpleCData, ...]:
    """Read and validate CP2K global/local grid metadata."""
    values = (
        *[ctypes.c_size_t() for _ in range(12)],
        *[ctypes.c_double() for _ in range(3)],
    )
    library.cp2k_aics_get_grid_metadata(*(ctypes.byref(value) for value in values))
    global_shape = tuple(value.value for value in values[0:3])
    local_shape = tuple(value.value for value in values[3:6])
    local_start = tuple(value.value for value in values[6:9])
    local_end = tuple(value.value for value in values[9:12])
    if any(extent == 0 for extent in global_shape):
        raise RuntimeError(
            f"CP2K returned an invalid zero global grid extent: {global_shape}"
        )
    for axis, (global_extent, local_extent, start, end) in enumerate(
        zip(
            global_shape,
            local_shape,
            local_start,
            local_end,
        )
    ):
        if end != start + local_extent or end > global_extent:
            raise RuntimeError(
                "CP2K returned inconsistent local grid metadata on axis "
                f"{axis}: global extent {global_extent}, start {start}, "
                f"extent {local_extent}, exclusive end {end}"
            )
    _validated_local_size(local_shape)
    return values


def _validated_local_size(local_shape: Tuple[int, int, int]) -> int:
    """Return a host-indexable shape product without fixed-width overflow."""
    if len(local_shape) != 3 or any(extent < 0 for extent in local_shape):
        raise RuntimeError(f"Invalid CP2K local array shape: {local_shape}")
    local_size = math.prod(local_shape)
    if local_size > np.iinfo(np.intp).max:
        raise RuntimeError(
            f"CP2K local array shape is too large for NumPy: {local_shape}"
        )
    return local_size


def _pointer_view(pointer, name: str, local_shape: Tuple[int, int, int]):
    """Build one synchronous CP2K view using the canonical empty contract."""
    local_size = _validated_local_size(local_shape)
    if local_size == 0:
        if pointer:
            raise RuntimeError(
                f"CP2K returned a non-null pointer for zero-sized {name}"
            )
        return np.empty(local_shape, dtype=np.float64, order="F")
    if not pointer:
        raise RuntimeError(f"CP2K returned a null pointer for nonempty {name}")
    raw = ctypes.cast(pointer, ctypes.POINTER(ctypes.c_double))
    return np.ctypeslib.as_array(raw, shape=(local_size,)).reshape(
        local_shape,
        order="F",
    )


def pointer_views(library, local_shape: Tuple[int, int, int]):
    """Return synchronous Fortran-order views of all CP2K coupling arrays."""
    getters = (
        "cp2k_aics_get_solute_charge_density",
        "cp2k_aics_get_relative_permittivity_parallel",
        "cp2k_aics_get_relative_permittivity_perpendicular",
        "cp2k_aics_get_electrostatic_potential",
    )
    names = (
        "solute_charge_density",
        "relative_permittivity_parallel",
        "relative_permittivity_perpendicular",
        "electrostatic_potential",
    )
    views = []
    for name, getter_name in zip(names, getters):
        pointer = getattr(library, getter_name)()
        views.append(_pointer_view(pointer, name, local_shape))
    return tuple(views)


def fd_sensitivity_pointer_views(library, local_shape: Tuple[int, int, int]):
    """Return zero-copy views of CP2K-owned FD sensitivity output arrays."""
    getters = (
        "cp2k_aics_get_fd_field_squared_sensitivity_parallel",
        "cp2k_aics_get_fd_field_squared_sensitivity_perpendicular",
    )
    names = (
        "dielectric_field_squared_sensitivity_parallel",
        "dielectric_field_squared_sensitivity_perpendicular",
    )
    views = []
    for name, getter_name in zip(names, getters):
        pointer = getattr(library, getter_name)()
        views.append(_pointer_view(pointer, name, local_shape))
    return tuple(views)
