"""Stable public lifecycle API for the CP2K AICS coupling."""

from pathlib import Path
from typing import Union

from .context import (
    finalize_context,
    get_context,
    initialize_prepared_context,
    is_context_initialized,
    prepare_context,
)


def initialize(
    bridge_library_path: Union[str, Path],
    coupling_dir: Union[str, Path],
    comm_f: int,
    poisson_solver: str = "finite_difference",
    normal_axis: int = 2,
    relative_tolerance: float = 1.0e-10,
    absolute_tolerance: float = 1.0e-12,
    max_iterations: int = 100000,
):
    """Initialize the current run's AICS context exactly once."""
    request = prepare_initialize(
        bridge_library_path,
        coupling_dir,
        comm_f,
        poisson_solver,
        normal_axis,
        relative_tolerance,
        absolute_tolerance,
        max_iterations,
    )
    return initialize_prepared(request)


def prepare_initialize(
    bridge_library_path: Union[str, Path],
    coupling_dir: Union[str, Path],
    comm_f: int,
    poisson_solver: str = "finite_difference",
    normal_axis: int = 2,
    relative_tolerance: float = 1.0e-10,
    absolute_tolerance: float = 1.0e-12,
    max_iterations: int = 100000,
):
    """Perform rank-local preparation before CP2K's status agreement."""
    return prepare_context(
        bridge_library_path,
        coupling_dir,
        comm_f,
        poisson_solver,
        normal_axis,
        relative_tolerance,
        absolute_tolerance,
        max_iterations,
    )


def initialize_prepared(request):
    """Enter communicator duplication after CP2K's status agreement."""
    return initialize_prepared_context(request)


def solve() -> None:
    """Solve using the currently initialized run context."""
    get_context().solve()


def finalize() -> None:
    """Release the current run's AICS context, if present."""
    finalize_context()


def is_initialized() -> bool:
    """Return whether an AICS context is initialized for the current run."""
    return is_context_initialized()
