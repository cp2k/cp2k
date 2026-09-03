"""Run-scoped context lifecycle for the CP2K AICS coupling."""

import time
from pathlib import Path
from typing import ClassVar, Optional, Union

import mpi4py

mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

_MAX_PHASE_CAUSE_LENGTH = 768


def _exception_cause(exception):
    """Return a bounded, single-line exception description."""
    try:
        exception_type = type(exception).__name__
    except BaseException:
        exception_type = "Exception"
    try:
        detail = " ".join(str(exception).split())
    except BaseException:
        detail = "exception message unavailable"
    cause = f"{exception_type}: {detail}" if detail else exception_type
    return cause[:_MAX_PHASE_CAUSE_LENGTH]


def run_local_phase(comm, phase_name, operation):
    """Run local work and raise the same rank-attributed error on every rank."""
    result = None
    local_exception = None
    local_cause = None
    try:
        result = operation()
    except BaseException as exception:
        local_exception = exception
        local_cause = _exception_cause(exception)

    failing_rank = comm.allreduce(
        comm.rank if local_exception is not None else comm.size,
        op=MPI.MIN,
    )
    if failing_rank == comm.size:
        return result

    cause = comm.bcast(
        local_cause if comm.rank == failing_rank else None,
        root=failing_rank,
    )
    error = RuntimeError(f"{phase_name} failed on rank {failing_rank}: {cause}")
    raise error from None


class AICSContext:
    """Own all library, grid, Poisson solver, and persistent state."""

    _current: ClassVar[Optional["AICSContext"]] = None

    def __init__(
        self,
        bridge_library_path: Path,
        coupling_dir: Path,
        comm,
        poisson_solver: str,
        normal_axis: int,
        relative_tolerance: float,
        absolute_tolerance: float,
        max_iterations: int,
    ) -> None:
        self._closed = False
        self.bridge_library_path = bridge_library_path
        self.mpi_comm = comm
        self.coupling_dir = coupling_dir
        self.poisson_solver_name = poisson_solver

        def validate_context_options():
            if normal_axis not in (0, 1, 2):
                raise ValueError(f"Invalid AICS Cartesian normal axis: {normal_axis}")
            if poisson_solver not in {"finite_difference", "finite_element"}:
                raise ValueError(f"Unknown AICS Poisson solver: {poisson_solver}")
            return normal_axis

        self.normal_axis = run_local_phase(
            self.mpi_comm,
            "AICS context option validation",
            validate_context_options,
        )
        if poisson_solver == "finite_element":

            def import_finite_element_backend():
                try:
                    from .finite_element_backend import FiniteElementPoissonSolver
                except ImportError as exc:
                    raise RuntimeError(
                        "AICS POISSON_SOLVER FINITE_ELEMENT was requested, but its "
                        "dependencies could not be imported."
                    ) from exc
                return FiniteElementPoissonSolver

            FiniteElementPoissonSolver = run_local_phase(
                self.mpi_comm,
                "AICS FINITE_ELEMENT backend import",
                import_finite_element_backend,
            )
            self.poisson_solver = FiniteElementPoissonSolver(
                bridge_library_path,
                coupling_dir,
                comm,
                relative_tolerance,
                absolute_tolerance,
                max_iterations,
                normal_axis=normal_axis,
            )
            return

        def prepare_finite_difference_backend():
            import numpy as np
            from .interface import (
                AICS_ARRAY_SYMBOLS,
                AICS_FD_SENSITIVITY_SYMBOLS,
                configure_grid_info,
                fd_sensitivity_pointer_views,
                load_aics_library,
                pointer_views,
            )
            from .finite_difference_backend import FiniteDifferencePoissonSolver

            self._np = np
            self._MPI = MPI
            self._pointer_views = pointer_views
            self._fd_sensitivity_pointer_views = fd_sensitivity_pointer_views
            self.lib = load_aics_library(
                self.bridge_library_path,
                AICS_ARRAY_SYMBOLS + AICS_FD_SENSITIVITY_SYMBOLS,
            )
            values = configure_grid_info(self.lib)
            (
                self.n1,
                self.n2,
                self.n3,
                self.ln1,
                self.ln2,
                self.ln3,
                self.l1,
                self.l2,
                self.l3,
                self.e1,
                self.e2,
                self.e3,
                self.dr1,
                self.dr2,
                self.dr3,
            ) = values
            self.nx, self.ny, self.nz = (self.n1.value, self.n2.value, self.n3.value)
            self.grid_spacing = np.array(
                [self.dr1.value, self.dr2.value, self.dr3.value],
                dtype=np.float64,
            )
            return FiniteDifferencePoissonSolver

        FiniteDifferencePoissonSolver = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_DIFFERENCE interface setup",
            prepare_finite_difference_backend,
        )
        self.poisson_solver = FiniteDifferencePoissonSolver(
            self.mpi_comm,
            (self.nx, self.ny, self.nz),
            self.grid_spacing,
            relative_tolerance,
            absolute_tolerance,
            max_iterations,
            normal_axis=normal_axis,
        )

    def _assemble_global(self, phase_name, local_array):
        def prepare_global_array():
            global_array = self._np.zeros(
                (self.nx, self.ny, self.nz),
                dtype=self._np.float64,
            )
            global_array[
                self.l1.value : self.e1.value,
                self.l2.value : self.e2.value,
                self.l3.value : self.e3.value,
            ] = local_array
            return global_array

        global_array = run_local_phase(
            self.mpi_comm,
            phase_name,
            prepare_global_array,
        )
        self.mpi_comm.Allreduce(
            self._MPI.IN_PLACE,
            global_array,
            op=self._MPI.SUM,
        )
        return global_array

    def solve(self) -> None:
        """Read CP2K arrays, solve with the selected discretization, and write phi."""
        if self.poisson_solver_name == "finite_element":
            self.poisson_solver.solve()
            return

        def prepare_pointer_views():
            local_shape = (self.ln1.value, self.ln2.value, self.ln3.value)
            arrays = self._pointer_views(self.lib, local_shape)
            sensitivities = self._fd_sensitivity_pointer_views(
                self.lib,
                local_shape,
            )
            return arrays + sensitivities

        (
            solute_charge_density,
            relative_permittivity_parallel,
            relative_permittivity_perpendicular,
            electrostatic_potential,
            dielectric_field_squared_sensitivity_parallel,
            dielectric_field_squared_sensitivity_perpendicular,
        ) = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_DIFFERENCE CP2K array view setup",
            prepare_pointer_views,
        )
        t0 = time.time()
        solute_charge_density_global = self._assemble_global(
            "AICS FINITE_DIFFERENCE charge-density assembly",
            solute_charge_density,
        )
        relative_permittivity_parallel_global = self._assemble_global(
            "AICS FINITE_DIFFERENCE parallel-permittivity assembly",
            relative_permittivity_parallel,
        )
        relative_permittivity_perpendicular_global = self._assemble_global(
            "AICS FINITE_DIFFERENCE perpendicular-permittivity assembly",
            relative_permittivity_perpendicular,
        )
        assembly_time = time.time() - t0
        self.last_assembly_time = assembly_time
        (
            electrostatic_potential_global,
            dielectric_field_squared_sensitivity_parallel_global,
            dielectric_field_squared_sensitivity_perpendicular_global,
        ) = self.poisson_solver.solve(
            solute_charge_density=solute_charge_density_global,
            relative_permittivity_parallel=relative_permittivity_parallel_global,
            relative_permittivity_perpendicular=relative_permittivity_perpendicular_global,
        )

        def store_local_results():
            local_slice = (
                slice(self.l1.value, self.e1.value),
                slice(self.l2.value, self.e2.value),
                slice(self.l3.value, self.e3.value),
            )
            electrostatic_potential[:, :, :] = electrostatic_potential_global[
                local_slice
            ]
            dielectric_field_squared_sensitivity_parallel[:, :, :] = (
                dielectric_field_squared_sensitivity_parallel_global[local_slice]
            )
            dielectric_field_squared_sensitivity_perpendicular[:, :, :] = (
                dielectric_field_squared_sensitivity_perpendicular_global[local_slice]
            )

        run_local_phase(
            self.mpi_comm,
            "AICS FINITE_DIFFERENCE CP2K result storage",
            store_local_results,
        )

    def close(self) -> None:
        """Release solver-owned resources and this context's references."""
        if self._closed:
            return
        self._closed = True
        first_error = None
        poisson_solver = getattr(self, "poisson_solver", None)
        try:
            comm = getattr(self, "mpi_comm", None)
            if comm is not None and comm != MPI.COMM_NULL:

                def close_poisson_solver():
                    close = getattr(poisson_solver, "close", None)
                    if callable(close):
                        close()

                run_local_phase(
                    comm,
                    "AICS Poisson solver cleanup",
                    close_poisson_solver,
                )
        except BaseException as exc:
            first_error = exc
        finally:
            self.poisson_solver = None
        try:
            comm = getattr(self, "mpi_comm", None)
            if comm is not None and comm != MPI.COMM_NULL:
                comm.Free()
        except BaseException as exc:
            if first_error is None:
                first_error = exc
        finally:
            self.mpi_comm = MPI.COMM_NULL
            if hasattr(self, "lib"):
                self.lib = None
        if first_error is not None:
            raise first_error


def prepare_context(
    bridge_library_path: Union[str, Path],
    coupling_dir: Union[str, Path],
    comm_f: int,
    poisson_solver: str = "finite_difference",
    normal_axis: int = 2,
    relative_tolerance: float = 1.0e-10,
    absolute_tolerance: float = 1.0e-12,
    max_iterations: int = 100000,
) -> tuple:
    """Validate and retain no state for a rank-local pre-Dup request."""
    if AICSContext._current is not None:
        raise RuntimeError(
            "AICS is already initialized; call finalize() before initializing again."
        )
    if comm_f is None:
        raise ValueError("AICS initialization requires a CP2K MPI communicator handle.")
    parent_comm = MPI.Comm.f2py(comm_f)
    if parent_comm == MPI.COMM_NULL:
        raise ValueError("AICS initialization received MPI.COMM_NULL.")
    if parent_comm.Is_inter():
        raise ValueError("AICS initialization requires an MPI intracommunicator.")
    if not isinstance(poisson_solver, str):
        raise ValueError("AICS Poisson solver must be a string.")
    selected_poisson_solver = poisson_solver.strip().lower()
    if selected_poisson_solver not in {"finite_difference", "finite_element"}:
        raise ValueError(
            f"Invalid AICS Poisson solver {selected_poisson_solver!r}; "
            "expected 'finite_difference' or 'finite_element'."
        )
    return (
        Path(bridge_library_path),
        Path(coupling_dir),
        parent_comm,
        selected_poisson_solver,
        normal_axis,
        relative_tolerance,
        absolute_tolerance,
        max_iterations,
    )


def initialize_prepared_context(request: tuple) -> AICSContext:
    """Collectively create a context from a unanimously prepared request."""
    (
        bridge_path,
        coupling_path,
        parent_comm,
        selected_poisson_solver,
        normal_axis,
        relative_tolerance,
        absolute_tolerance,
        max_iterations,
    ) = request
    comm = parent_comm.Dup()
    try:
        selected_poisson_solvers = comm.allgather(selected_poisson_solver)
        if any(
            selection != selected_poisson_solvers[0]
            for selection in selected_poisson_solvers[1:]
        ):
            raise RuntimeError(
                "AICS initialization received inconsistent Poisson solvers "
                "across active ranks."
            )
        context = AICSContext(
            bridge_path,
            coupling_path,
            comm,
            selected_poisson_solver,
            normal_axis,
            relative_tolerance,
            absolute_tolerance,
            max_iterations,
        )
    except BaseException:
        if comm != MPI.COMM_NULL:
            comm.Free()
        raise
    AICSContext._current = context
    return context


def get_context() -> AICSContext:
    """Return the initialized context or fail with an actionable error."""
    if AICSContext._current is None:
        raise RuntimeError(
            "AICS solve called before initialization; call "
            "cp2k_aics.initialize() first."
        )
    return AICSContext._current


def finalize_context() -> None:
    """Close and release the process-wide context, if one exists."""
    context = AICSContext._current
    if context is None:
        return
    try:
        context.close()
    finally:
        AICSContext._current = None


def is_context_initialized() -> bool:
    """Return whether a process-wide context is active."""
    return AICSContext._current is not None
