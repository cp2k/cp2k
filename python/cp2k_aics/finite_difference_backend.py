"""Distributed periodic finite-difference PETSc Poisson solver for CP2K AICS.

The solver discretizes ``-div(epsilon grad(phi)) = 4*pi*rho`` on a periodic,
uniform Cartesian grid. The interface-normal axis uses ``epsilon_perpendicular``
and the other axes use ``epsilon_parallel``. Arithmetic face averaging gives a
seven-point stencil. The constant nullspace is attached, the right-hand side is
projected for compatibility, and the result is shifted to a deterministic
zero-mean gauge. Backend-exact dielectric sensitivities are returned for the
density-dependent KS dielectric chain and direct interface force.
"""

import time

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from .context import run_local_phase


def _petsc_ksp_reason_name(reason):
    """Return PETSc's symbolic KSP reason for concise failure reporting."""
    reason_value = int(reason)
    # petsc4py 3.25.3 retains an older MAX_IT alias, while PETSc itself names
    # KSP reason -3 KSP_DIVERGED_ITS.
    if reason_value == -3:
        return "DIVERGED_ITS"
    for name in dir(PETSc.KSP.ConvergedReason):
        if not name.isupper() or name == "DIVERGED_MAX_IT":
            continue
        try:
            if int(getattr(PETSc.KSP.ConvergedReason, name)) == reason_value:
                return name
        except (TypeError, ValueError):
            continue
    return f"reason {reason_value}"


def _ksp_failure_message(
    backend, ksp, relative_tolerance, absolute_tolerance, max_iterations
):
    """Build a bounded-state diagnostic without masking the PETSc exception."""
    reason_name = None
    iterations = None
    residual = None
    try:
        reason_name = _petsc_ksp_reason_name(ksp.getConvergedReason())
    except Exception:
        pass
    try:
        iterations = int(ksp.getIterationNumber())
    except Exception:
        pass
    try:
        residual = float(ksp.getResidualNorm())
    except Exception:
        pass

    message = f"AICS {backend} Poisson solve failed: PETSc KSP"
    if reason_name is not None:
        message += f" {reason_name}"
    else:
        message += " did not converge"
    if iterations is not None:
        suffix = "iteration" if iterations == 1 else "iterations"
        message += f" after {iterations} {suffix}"
    if residual is not None:
        message += f"; residual norm {residual:.6E}"
    message += (
        f"; configured EPS_REL {relative_tolerance:.6E}, "
        f"EPS_ABS {absolute_tolerance:.6E}, MAX_ITER {max_iterations}."
    )
    return message


def dielectric_field_squared_sensitivities(potential, grid_spacing, normal_axis):
    """Return exact FD field-squared sensitivities for the two dielectric DOFs.

    Both returned cell-centred arrays have units of potential squared per length
    squared.  They exclude the electrostatic energy prefactor ``-dV/(8*pi)``.
    The parallel result is the sum for the two axes transverse to
    ``normal_axis``; the perpendicular result is the normal-axis sensitivity.
    """
    potential = np.asarray(potential, dtype=np.float64)
    grid_spacing = np.asarray(grid_spacing, dtype=np.float64)
    if potential.ndim != 3:
        raise ValueError(f"potential must be three-dimensional, got {potential.shape}")
    if grid_spacing.shape != (3,) or np.any(grid_spacing <= 0.0):
        raise ValueError(
            f"grid_spacing must contain three positive values, got {grid_spacing}"
        )
    if normal_axis not in (0, 1, 2):
        raise ValueError(f"normal_axis must be 0, 1, or 2, got {normal_axis}")
    cartesian = tuple(
        (
            (potential - np.roll(potential, -1, axis=axis)) ** 2
            + (potential - np.roll(potential, +1, axis=axis)) ** 2
        )
        / (2.0 * grid_spacing[axis] ** 2)
        for axis in range(3)
    )
    transverse_axes = tuple(axis for axis in range(3) if axis != normal_axis)
    sensitivity_parallel = cartesian[transverse_axes[0]] + cartesian[transverse_axes[1]]
    sensitivity_perpendicular = cartesian[normal_axis]
    return sensitivity_parallel, sensitivity_perpendicular


class FiniteDifferencePoissonSolver:
    """Periodic variable-coefficient seven-point Poisson solver."""

    def __init__(
        self,
        comm,
        grid_shape,
        grid_spacing,
        relative_tolerance=1.0e-10,
        absolute_tolerance=1.0e-12,
        max_iterations=100000,
        normal_axis=2,
    ) -> None:
        self.comm = comm
        self.relative_tolerance = relative_tolerance
        self.absolute_tolerance = absolute_tolerance
        self.max_iterations = max_iterations
        self.normal_axis = normal_axis
        self._has_solution = False

        def prepare_grid_metadata():
            self.grid_shape = tuple(int(v) for v in grid_shape)
            self.grid_spacing = np.asarray(grid_spacing, dtype=np.float64)
            self.n_unknowns = int(np.prod(self.grid_shape, dtype=np.int64))
            if any(v < 2 for v in self.grid_shape):
                raise ValueError(
                    f"Each grid dimension must be at least 2, got {self.grid_shape}"
                )
            if np.any(self.grid_spacing <= 0.0):
                raise ValueError(
                    f"Grid spacings must be positive, got {self.grid_spacing}"
                )

        run_local_phase(
            self.comm,
            "AICS FINITE_DIFFERENCE grid metadata setup",
            prepare_grid_metadata,
        )
        t0 = time.time()
        self._create_matrix_structure()
        self._create_solver()
        self.communicator_sizes = tuple(
            obj.getComm().getSize()
            for obj in (self.operator, self.rhs, self.potential, self.ksp)
        )
        self.n_coo_stencil_entries = 7 * self.n_unknowns
        self.initialization_time = time.time() - t0

    def _create_matrix_structure(self) -> None:
        def prepare_stencil_indices():
            index = np.arange(self.n_unknowns, dtype=PETSc.IntType).reshape(
                self.grid_shape, order="F"
            )
            center = index.ravel(order="F")
            neighbours = (
                center,
                np.roll(index, -1, axis=0).ravel(order="F"),
                np.roll(index, +1, axis=0).ravel(order="F"),
                np.roll(index, -1, axis=1).ravel(order="F"),
                np.roll(index, +1, axis=1).ravel(order="F"),
                np.roll(index, -1, axis=2).ravel(order="F"),
                np.roll(index, +1, axis=2).ravel(order="F"),
            )
            return center, neighbours

        center, neighbours = run_local_phase(
            self.comm,
            "AICS FINITE_DIFFERENCE stencil-index setup",
            prepare_stencil_indices,
        )
        self.operator = PETSc.Mat().create(comm=self.comm)
        self.operator.setSizes((self.n_unknowns, self.n_unknowns))
        self.operator.setType(PETSc.Mat.Type.AIJ)
        self.row_start, self.row_end = self.operator.getOwnershipRange()

        def prepare_local_matrix_storage():
            owned = slice(self.row_start, self.row_end)
            owned_rows = center[owned]
            coo_rows = np.concatenate((owned_rows,) * 7)
            coo_cols = np.concatenate(tuple(columns[owned] for columns in neighbours))
            self.n_local_unknowns = self.row_end - self.row_start
            self.coo_values = np.empty(
                7 * self.n_local_unknowns,
                dtype=PETSc.ScalarType,
            )
            return coo_rows, coo_cols

        coo_rows, coo_cols = run_local_phase(
            self.comm,
            "AICS FINITE_DIFFERENCE local matrix storage setup",
            prepare_local_matrix_storage,
        )
        self.operator.setPreallocationCOO(coo_rows, coo_cols)
        self.operator.setOption(PETSc.Mat.Option.SYMMETRIC, True)
        self.potential, self.rhs = self.operator.createVecs()

        def validate_local_ownership():
            if self.rhs.getOwnershipRange() != (self.row_start, self.row_end):
                raise RuntimeError("PETSc matrix and vector ownership ranges differ")

        run_local_phase(
            self.comm,
            "AICS FINITE_DIFFERENCE PETSc ownership validation",
            validate_local_ownership,
        )
        ownership_ranges = self.comm.allgather((self.row_start, self.row_end))

        def prepare_solution_layout():
            expected_start = 0
            for start, end in ownership_ranges:
                if start != expected_start or end < start:
                    raise RuntimeError(
                        "PETSc ownership ranges do not partition the grid"
                    )
                expected_start = end
            if expected_start != self.n_unknowns:
                raise RuntimeError("PETSc ownership ranges do not cover the grid")
            self._solution_counts = np.asarray(
                [end - start for start, end in ownership_ranges],
                dtype=np.int64,
            )
            self._solution_displacements = np.asarray(
                [start for start, _ in ownership_ranges],
                dtype=np.int64,
            )

        run_local_phase(
            self.comm,
            "AICS FINITE_DIFFERENCE solution-layout setup",
            prepare_solution_layout,
        )
        self.null_space = PETSc.NullSpace().create(constant=True, comm=self.comm)
        self.operator.setNullSpace(self.null_space)
        self.operator.setNearNullSpace(self.null_space)

    def _create_solver(self) -> None:
        self.ksp = PETSc.KSP().create(comm=self.comm)
        self.ksp.setOptionsPrefix("aics_finite_difference_")
        self.ksp.setType(PETSc.KSP.Type.CG)
        self.ksp.getPC().setType(PETSc.PC.Type.GAMG)
        self.ksp.setTolerances(
            rtol=self.relative_tolerance,
            atol=self.absolute_tolerance,
            max_it=self.max_iterations,
        )
        self.ksp.setOperators(self.operator)
        self.ksp.setErrorIfNotConverged(True)

        self.ksp.setFromOptions()

    @staticmethod
    def _check_field(name, array, shape):
        array = np.asarray(array, dtype=np.float64)
        if array.shape != shape:
            raise ValueError(f"{name} has shape {array.shape}, expected {shape}")
        if not np.all(np.isfinite(array)):
            raise FloatingPointError(f"{name} contains NaN or infinity")
        return array

    def _prepare_matrix_values(
        self,
        relative_permittivity_parallel,
        relative_permittivity_perpendicular,
    ) -> None:
        # The parallel coefficient applies to both directions in the interface plane.
        epsilon_parallel = self._check_field(
            "relative_permittivity_parallel",
            relative_permittivity_parallel,
            self.grid_shape,
        )
        epsilon_perpendicular = self._check_field(
            "relative_permittivity_perpendicular",
            relative_permittivity_perpendicular,
            self.grid_shape,
        )
        eps = (epsilon_parallel, epsilon_perpendicular)
        if min(float(np.min(component)) for component in eps) <= 0.0:
            raise ValueError(
                "All dielectric tensor diagonal components must be positive"
            )
        diagonal = np.zeros(self.grid_shape, dtype=np.float64)
        block = self.n_local_unknowns
        inverse_h2 = 1.0 / np.square(self.grid_spacing)
        components = [epsilon_parallel] * 3
        components[self.normal_axis] = epsilon_perpendicular
        directions = tuple(
            (components[axis], axis, shift, inverse_h2[axis], 1 + 2 * axis + side)
            for axis in range(3)
            for side, shift in enumerate((-1, +1))
        )
        for component, axis, shift, scale, block_id in directions:
            face = 0.5 * (component + np.roll(component, shift, axis=axis)) * scale
            diagonal += face
            face_flat = face.ravel(order="F")
            self.coo_values[block_id * block : (block_id + 1) * block] = -face_flat[
                self.row_start : self.row_end
            ]
        diagonal_flat = diagonal.ravel(order="F")
        self.coo_values[0:block] = diagonal_flat[self.row_start : self.row_end]

    def solve(
        self,
        solute_charge_density,
        relative_permittivity_parallel,
        relative_permittivity_perpendicular,
    ):
        t0 = time.time()

        def prepare_operator_update():
            checked_charge_density = self._check_field(
                "solute_charge_density", solute_charge_density, self.grid_shape
            )
            self._prepare_matrix_values(
                relative_permittivity_parallel,
                relative_permittivity_perpendicular,
            )
            return checked_charge_density

        solute_charge_density = run_local_phase(
            self.comm,
            "AICS FINITE_DIFFERENCE field validation and matrix update",
            prepare_operator_update,
        )
        self.operator.setValuesCOO(self.coo_values, addv=PETSc.InsertMode.INSERT_VALUES)
        matrix_time = time.time() - t0

        def prepare_linear_solve():
            poisson_rhs = 4.0 * np.pi * solute_charge_density.ravel(order="F")
            rhs_mean_before = float(np.mean(poisson_rhs))
            rhs_scale = max(float(np.max(np.abs(poisson_rhs))), 1.0)
            self.last_rhs_mean_before_projection = rhs_mean_before
            self.last_rhs_scale = rhs_scale
            self.rhs.getArray()[:] = poisson_rhs[self.row_start : self.row_end]
            self.null_space.remove(self.rhs)
            if not self._has_solution:
                self.potential.set(0.0)
            self.ksp.setInitialGuessNonzero(self._has_solution)
            self.ksp.setOperators(self.operator)

        run_local_phase(
            self.comm,
            "AICS FINITE_DIFFERENCE PETSc solve preparation",
            prepare_linear_solve,
        )
        t1 = time.time()
        try:
            self.ksp.solve(self.rhs, self.potential)
        except PETSc.Error as exc:
            error = RuntimeError(
                _ksp_failure_message(
                    "FINITE_DIFFERENCE",
                    self.ksp,
                    self.relative_tolerance,
                    self.absolute_tolerance,
                    self.max_iterations,
                )
            )
            error.add_note(f"Original PETSc error:\n{exc}")
            raise error from None
        solve_time = time.time() - t1
        self._has_solution = True
        potential_mean = float(self.potential.sum()) / float(self.n_unknowns)
        self.potential.shift(-potential_mean)

        def prepare_potential_gather():
            reason = self.ksp.getConvergedReason()
            iterations = self.ksp.getIterationNumber()
            residual = self.ksp.getResidualNorm()
            self.last_iterations = iterations
            self.last_residual = residual
            self.last_convergence_reason = reason
            self.last_matrix_update_time = matrix_time
            self.last_linear_solve_time = solve_time
            return (
                np.empty(self.n_unknowns, dtype=np.float64),
                self.potential.getArray(readonly=True),
            )

        potential_global, potential_local = run_local_phase(
            self.comm,
            "AICS FINITE_DIFFERENCE potential-gather setup",
            prepare_potential_gather,
        )
        self.comm.Allgatherv(
            [potential_local, MPI.DOUBLE],
            [
                potential_global,
                self._solution_counts,
                self._solution_displacements,
                MPI.DOUBLE,
            ],
        )

        def prepare_return_fields():
            potential = np.asarray(potential_global, dtype=np.float64).reshape(
                self.grid_shape,
                order="F",
            )
            sensitivity_parallel, sensitivity_perpendicular = (
                dielectric_field_squared_sensitivities(
                    potential,
                    self.grid_spacing,
                    self.normal_axis,
                )
            )
            return potential, sensitivity_parallel, sensitivity_perpendicular

        return run_local_phase(
            self.comm,
            "AICS FINITE_DIFFERENCE sensitivity and return-field setup",
            prepare_return_fields,
        )
