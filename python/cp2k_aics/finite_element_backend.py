"""Persistent distributed periodic CG1 finite-element solver for CP2K AICS.

For total solute charge density ``rho`` and relative permittivities, the weak
form assigns ``epsilon_perpendicular`` to the selected interface-normal axis and
``epsilon_parallel`` to the other two Cartesian axes. The DOLFINx mesh and CG1
space are periodic, and the PETSc operator has a constant nullspace. With the
current zero-initial-guess behavior, PETSc returns the effective compatible,
nullspace-orthogonal solution; this backend does not explicitly remove the RHS
null component or apply a separate zero-mean shift.
"""

import os
import tempfile
import time
from pathlib import Path

import dolfinx
import dolfinx.fem.petsc
import dolfinx.graph
import numpy as np
import ufl
from basix.ufl import element
from dolfinx import fem, mesh
from dolfinx.geometry import (
    bb_tree,
    compute_colliding_cells,
    compute_collisions_points,
)
from dolfinx.io import XDMFFile
from dolfinx.mesh import create_mesh
from mpi4py import MPI
from petsc4py import PETSc
from ufl import TestFunction, TrialFunction, dx, grad

from .context import run_local_phase
from .interface import (
    AICS_ARRAY_SYMBOLS,
    configure_grid_info,
    load_aics_library,
    pointer_views,
)


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


class FiniteElementPoissonSolver:
    """Own the mesh, coefficient functions, and persistent finite-element solver state."""

    def __init__(
        self,
        bridge_library_path: Path,
        coupling_dir: Path,
        comm,
        relative_tolerance=1.0e-10,
        absolute_tolerance=1.0e-12,
        max_iterations=100000,
        normal_axis=2,
    ) -> None:
        self._initialized = False
        self.bridge_library_path = bridge_library_path
        self.coupling_dir = coupling_dir
        self.mpi_comm = comm
        self.relative_tolerance = relative_tolerance
        self.absolute_tolerance = absolute_tolerance
        self.max_iterations = max_iterations
        self.normal_axis = normal_axis
        self.problem = None
        self.null_space = None
        self.phi_solution = None
        self.initialize()

    def initialize(self) -> None:
        """Perform the one-time finite-element setup."""
        self._initialized = False
        self.problem = None
        self.null_space = None
        self.phi_solution = None

        def prepare_grid_metadata():
            self.lib = load_aics_library(
                self.bridge_library_path,
                AICS_ARRAY_SYMBOLS,
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

            self.nx = self.n1.value
            self.ny = self.n2.value
            self.nz = self.n3.value
            self.grid_spacing = np.array(
                [self.dr1.value, self.dr2.value, self.dr3.value]
            )
            self.Lx, self.Ly, self.Lz = self.grid_spacing * np.array(
                [self.nx, self.ny, self.nz]
            )

        run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT grid metadata setup",
            prepare_grid_metadata,
        )
        cartesian_mesh = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT Cartesian mesh creation",
            lambda: mesh.create_box(
                self.mpi_comm,
                [[0.0, 0.0, 0.0], [self.Lx, self.Ly, self.Lz]],
                [self.nx, self.ny, self.nz],
                mesh.CellType.hexahedron,
            ),
        )

        def create_temporary_directory():
            temporary_parent = self.coupling_dir.resolve()
            if self.mpi_comm.rank == 0:
                directory = tempfile.TemporaryDirectory(
                    prefix="cp2k-aics-fe-",
                    dir=temporary_parent,
                )
                return directory, str(Path(directory.name).resolve())
            return None, None

        temporary_directory, root_temporary_path = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT temporary-directory setup",
            create_temporary_directory,
        )
        temporary_path = self.mpi_comm.bcast(
            root_temporary_path,
            root=0,
        )
        mesh_phase_error = None
        try:

            def prepare_mesh_path():
                mesh_path = str(Path(temporary_path) / "mesh.xdmf")
                temporary_directory_visible = Path(
                    temporary_path
                ).is_dir() and os.access(
                    temporary_path,
                    os.R_OK | os.W_OK | os.X_OK,
                )
                return mesh_path, temporary_directory_visible

            mesh_path, temporary_directory_visible = run_local_phase(
                self.mpi_comm,
                "AICS FINITE_ELEMENT temporary-directory visibility check",
                prepare_mesh_path,
            )
            temporary_directory_visible_everywhere = self.mpi_comm.allreduce(
                temporary_directory_visible,
                op=MPI.LAND,
            )
            if not temporary_directory_visible_everywhere:
                raise RuntimeError(
                    "The FEM temporary directory is not visible and accessible to "
                    "every active AICS rank; it must reside on a shared filesystem."
                )

            def write_mesh():
                with XDMFFile(self.mpi_comm, mesh_path, "w") as file:
                    file.write_mesh(cartesian_mesh)

            run_local_phase(
                self.mpi_comm,
                "AICS FINITE_ELEMENT collective mesh write",
                write_mesh,
            )

            # Read all geometry data on all processes.
            def read_global_geometry():
                with XDMFFile(MPI.COMM_SELF, mesh_path, "r") as file:
                    return file.read_geometry_data()

            global_coordinates = run_local_phase(
                self.mpi_comm,
                "AICS FINITE_ELEMENT replicated geometry read",
                read_global_geometry,
            )

            # Read topology data.
            def read_partitioned_mesh():
                with XDMFFile(self.mpi_comm, mesh_path, "r") as file:
                    return (
                        file.read_cell_type(),
                        file.read_geometry_data(),
                        file.read_topology_data(),
                    )

            (cell_shape, cell_degree), local_coordinates, cell_topology = (
                run_local_phase(
                    self.mpi_comm,
                    "AICS FINITE_ELEMENT collective mesh read",
                    read_partitioned_mesh,
                )
            )
        except BaseException as exc:
            mesh_phase_error = exc

        try:

            def cleanup_temporary_directory():
                if temporary_directory is not None:
                    temporary_directory.cleanup()

            run_local_phase(
                self.mpi_comm,
                "AICS FINITE_ELEMENT temporary-directory cleanup",
                cleanup_temporary_directory,
            )
        except BaseException as exc:
            if mesh_phase_error is None:
                mesh_phase_error = exc
        if mesh_phase_error is not None:
            raise mesh_phase_error

        num_local_coor = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT local geometry-size setup",
            lambda: local_coordinates.shape[0],
        )
        all_sizes = self.mpi_comm.allgather(num_local_coor)

        def prepare_repartitioning():
            all_sizes.insert(0, 0)
            all_ranges = np.cumsum(all_sizes)

            # Coordinates must be read contiguously in rank-local chunks.
            self.rank = self.mpi_comm.rank
            if not np.all(
                global_coordinates[all_ranges[self.rank] : all_ranges[self.rank + 1]]
                == local_coordinates
            ):
                raise RuntimeError(
                    "XDMF geometry is not partitioned into contiguous rank-local "
                    "chunks."
                )

            coordinate_domain = ufl.Mesh(
                element("Lagrange", cell_shape.name, cell_degree, shape=(3,))
            )
            midpoints = np.mean(global_coordinates[cell_topology], axis=1)
            dest = np.zeros(len(midpoints), dtype=np.int32)

            max_x, min_x = (
                float(self.Lx) - 2.0 * float(self.Lx) / float(self.nx),
                2.0 * float(self.Lx) / float(self.nx),
            )
            max_y, min_y = (
                float(self.Ly) - 2.0 * float(self.Ly) / float(self.ny),
                2.0 * float(self.Ly) / float(self.ny),
            )
            max_z, min_z = (
                float(self.Lz) - 2.0 * float(self.Lz) / float(self.nz),
                2.0 * float(self.Lz) / float(self.nz),
            )

            boundary_indices = np.where(
                ((midpoints[:, 0] - max_x) > 0.0)
                | ((midpoints[:, 0] - min_x) < 0.0)
                | ((midpoints[:, 1] - max_y) > 0.0)
                | ((midpoints[:, 1] - min_y) < 0.0)
                | ((midpoints[:, 2] - max_z) > 0.0)
                | ((midpoints[:, 2] - min_z) < 0.0)
            )[0]
            dest[boundary_indices] = 0

            if self.mpi_comm.size > 1:
                to_distribute = np.setdiff1d(
                    np.arange(len(midpoints)), boundary_indices
                )
                remaining_ranks = np.arange(1, self.mpi_comm.size)
                num_remaining_points = len(to_distribute)
                num_ranks = len(remaining_ranks)
                repeated_ranks = np.repeat(
                    remaining_ranks, num_remaining_points // num_ranks
                )
                final_dest = np.zeros(len(to_distribute), dtype=np.int32)
                final_dest[: len(repeated_ranks)] = repeated_ranks
                dest[to_distribute] = final_dest
            adjacency = dolfinx.cpp.graph.AdjacencyList_int32(dest)
            return coordinate_domain, adjacency

        coordinate_domain, partition_adjacency = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT repartitioning metadata setup",
            prepare_repartitioning,
        )

        # Partition mesh in layers, capture geometrical and topological data.
        def partitioner(*args):
            return partition_adjacency

        partitioned_mesh = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT distributed mesh repartitioning",
            lambda: create_mesh(
                self.mpi_comm,
                cell_topology,
                coordinate_domain,
                local_coordinates,
                partitioner,
            ),
        )

        def validate_repartitioned_mesh():
            tdim = partitioned_mesh.topology.dim
            if not (
                cartesian_mesh.topology.index_map(tdim).size_global
                == partitioned_mesh.topology.index_map(tdim).size_global
            ):
                raise RuntimeError("Mesh repartitioning changed the global cell count.")

        run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT repartitioned-mesh validation",
            validate_repartitioned_mesh,
        )
        self.periodic_mesh = self._create_periodic_mesh(partitioned_mesh)
        self.function_space = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT function-space construction",
            lambda: dolfinx.fem.functionspace(self.periodic_mesh, ("CG", 1)),
        )

        def prepare_reconstruction_metadata():
            self._reconstruction_points = self.periodic_mesh.geometry.x[
                : self.periodic_mesh.geometry.index_map().size_local
            ]
            tree = bb_tree(self.periodic_mesh, self.periodic_mesh.topology.dim)
            cell_candidates = compute_collisions_points(
                tree, self._reconstruction_points
            )
            colliding_cells = compute_colliding_cells(
                self.periodic_mesh, cell_candidates, self._reconstruction_points
            )
            cells = []
            for i in range(len(self._reconstruction_points)):
                candidates = colliding_cells.links(i)
                if len(candidates) == 0:
                    raise ValueError(
                        f"Point {self._reconstruction_points[i]} is outside the mesh."
                    )
                cells.append(candidates[0])
            self._reconstruction_cells = np.array(cells, dtype=np.int32)
            return self._reconstruction_points.size

        self._reconstruction_local_coord_size = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT reconstruction metadata setup",
            prepare_reconstruction_metadata,
        )
        self._reconstruction_recv_counts_coords = self.mpi_comm.gather(
            self._reconstruction_local_coord_size, root=0
        )

        def prepare_coordinate_gather():
            if self.rank == 0:
                global_size = self.periodic_mesh.geometry.index_map().size_global
                self._reconstruction_global_coords = np.zeros(
                    (global_size, self._reconstruction_points.shape[1]),
                    dtype=self._reconstruction_points.dtype,
                )
                self._reconstruction_displs_coords = np.cumsum(
                    [0] + self._reconstruction_recv_counts_coords[:-1]
                )
            else:
                self._reconstruction_global_coords = None
                self._reconstruction_displs_coords = None

        run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT coordinate-gather setup",
            prepare_coordinate_gather,
        )
        self.mpi_comm.Gatherv(
            sendbuf=self._reconstruction_points,
            recvbuf=(
                self._reconstruction_global_coords,
                self._reconstruction_recv_counts_coords,
                self._reconstruction_displs_coords,
                MPI.DOUBLE,
            ),
            root=0,
        )

        def prepare_coordinate_sort():
            if self.rank == 0:
                self._reconstruction_sorted_indices = np.lexsort(
                    (
                        self._reconstruction_global_coords[:, 2],
                        self._reconstruction_global_coords[:, 1],
                        self._reconstruction_global_coords[:, 0],
                    )
                )
            else:
                self._reconstruction_sorted_indices = None

        run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT coordinate-sort setup",
            prepare_coordinate_sort,
        )
        self.solute_charge_density = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT charge function construction",
            lambda: fem.Function(self.function_space),
        )
        self.relative_permittivity_parallel = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT parallel-permittivity function construction",
            lambda: fem.Function(self.function_space),
        )
        self.relative_permittivity_perpendicular = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT perpendicular-permittivity function construction",
            lambda: fem.Function(self.function_space),
        )

        self.petsc_options = {
            "ksp_type": "cg",
            "pc_type": "gamg",
            "ksp_rtol": self.relative_tolerance,
            "ksp_atol": self.absolute_tolerance,
            "ksp_max_it": self.max_iterations,
            "ksp_error_if_not_converged": True,
        }
        self._initialized = True

    def _indicator(self, x):
        tol = 1e-10
        is_corner = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[1, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
        )
        is_z_axis = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[1, :], 0.0, atol=tol)
            & ~np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_y_axis = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[1, :], 0.0, atol=tol)
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
        )
        is_x_axis = (
            np.isclose(x[1, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[0, :], 0.0, atol=tol)
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
        )
        is_yz_face = (
            np.isclose(x[0, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_zx_face = (
            np.isclose(x[1, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_xy_face = (
            np.isclose(x[2, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
        )
        return (
            is_corner
            | is_x_axis
            | is_y_axis
            | is_z_axis
            | is_xy_face
            | is_yz_face
            | is_zx_face
        )

    def _indicator_dual(self, x):
        tol = 1e-10
        is_corner = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[1, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
        )
        is_z_axis = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[1, :], 0.0, atol=tol)
            & ~np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_y_axis = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[1, :], 0.0, atol=tol)
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
        )
        is_x_axis = (
            np.isclose(x[1, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[0, :], 0.0, atol=tol)
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
        )
        is_yz_face = (
            np.isclose(x[0, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_zx_face = (
            np.isclose(x[1, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_xy_face = (
            np.isclose(x[2, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
        )
        x_condition = np.isclose(x[0, :], 0.0, atol=tol) | np.isclose(
            x[0, :], self.Lx, atol=tol
        )
        y_condition = np.isclose(x[1, :], 0.0, atol=tol) | np.isclose(
            x[1, :], self.Ly, atol=tol
        )
        z_condition = np.isclose(x[2, :], 0.0, atol=tol) | np.isclose(
            x[2, :], self.Lz, atol=tol
        )
        combined_condition = x_condition | y_condition | z_condition
        return (
            ~(
                is_corner
                | is_x_axis
                | is_y_axis
                | is_z_axis
                | is_xy_face
                | is_yz_face
                | is_zx_face
            )
            & combined_condition
        )

    def _mapping_function(self, x):
        tol = 1e-10
        is_corner = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[1, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
        )
        is_z_axis = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[1, :], 0.0, atol=tol)
            & ~np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_y_axis = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[1, :], 0.0, atol=tol)
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
        )
        is_x_axis = (
            np.isclose(x[1, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[0, :], 0.0, atol=tol)
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
        )
        is_yz_face = (
            np.isclose(x[0, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_zx_face = (
            np.isclose(x[1, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_xy_face = (
            np.isclose(x[2, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
        )

        yz_face_x_image = x[:, is_yz_face]
        yz_face_x_image[0, :] += self.Lx
        zx_face_y_image = x[:, is_zx_face]
        zx_face_y_image[1, :] += self.Ly
        xy_face_z_image = x[:, is_xy_face]
        xy_face_z_image[2, :] += self.Lz
        x_axis_y_image = x[:, is_x_axis]
        x_axis_y_image[1, :] += self.Ly
        x_axis_z_image = x[:, is_x_axis]
        x_axis_z_image[2, :] += self.Lz
        x_axis_yz_image = x[:, is_x_axis]
        x_axis_yz_image[2, :] += self.Lz
        x_axis_yz_image[1, :] += self.Ly
        y_axis_x_image = x[:, is_y_axis]
        y_axis_x_image[0, :] += self.Lx
        y_axis_z_image = x[:, is_y_axis]
        y_axis_z_image[2, :] += self.Lz
        y_axis_xz_image = x[:, is_y_axis]
        y_axis_xz_image[0, :] += self.Lx
        y_axis_xz_image[2, :] += self.Lz
        z_axis_y_image = x[:, is_z_axis]
        z_axis_y_image[1, :] += self.Ly
        z_axis_x_image = x[:, is_z_axis]
        z_axis_x_image[0, :] += self.Lx
        z_axis_xy_image = x[:, is_z_axis]
        z_axis_xy_image[0, :] += self.Lx
        z_axis_xy_image[1, :] += self.Ly
        if np.any(is_corner):
            corner_images = np.array(
                [
                    [self.Lx, 0.0, 0.0, self.Lx, 0.0, self.Lx, self.Lx],
                    [0.0, self.Ly, 0.0, self.Ly, self.Ly, 0.0, self.Ly],
                    [0.0, 0.0, self.Lz, 0.0, self.Lz, self.Lz, self.Lz],
                ]
            )
        else:
            corner_images = np.zeros((3, 0), dtype=x.dtype)
        periodic_images = np.hstack(
            (
                yz_face_x_image,
                zx_face_y_image,
                xy_face_z_image,
                x_axis_y_image,
                x_axis_z_image,
                x_axis_yz_image,
                y_axis_x_image,
                y_axis_z_image,
                y_axis_xz_image,
                z_axis_y_image,
                z_axis_x_image,
                z_axis_xy_image,
                corner_images,
            )
        )
        return periodic_images

    def _mapping_function_inv(self, x, left_vertices):
        tol = 1e-10
        is_corner = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[1, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
        )
        is_z_axis = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[1, :], 0.0, atol=tol)
            & ~np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_y_axis = (
            np.isclose(x[0, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[1, :], 0.0, atol=tol)
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
        )
        is_x_axis = (
            np.isclose(x[1, :], 0.0, atol=tol)
            & np.isclose(x[2, :], 0.0, atol=tol)
            & ~np.isclose(x[0, :], 0.0, atol=tol)
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
        )
        is_yz_face = (
            np.isclose(x[0, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_zx_face = (
            np.isclose(x[1, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
            & ~np.isclose(x[2, :], self.Lz, atol=tol)
        )
        is_xy_face = (
            np.isclose(x[2, :], 0, atol=tol)
            & ~is_corner
            & ~is_x_axis
            & ~is_y_axis
            & ~is_z_axis
            & ~np.isclose(x[0, :], self.Lx, atol=tol)
            & ~np.isclose(x[1, :], self.Ly, atol=tol)
        )
        yz_face_sources = left_vertices[is_yz_face]
        zx_face_sources = left_vertices[is_zx_face]
        xy_face_sources = left_vertices[is_xy_face]
        x_axis_y_sources = left_vertices[is_x_axis]
        x_axis_z_sources = left_vertices[is_x_axis]
        x_axis_yz_sources = left_vertices[is_x_axis]
        y_axis_x_sources = left_vertices[is_y_axis]
        y_axis_z_sources = left_vertices[is_y_axis]
        y_axis_xz_sources = left_vertices[is_y_axis]
        z_axis_y_sources = left_vertices[is_z_axis]
        z_axis_x_sources = left_vertices[is_z_axis]
        z_axis_xy_sources = left_vertices[is_z_axis]
        corner_column = left_vertices[is_corner]
        corner_sources = np.repeat(corner_column, 7)
        replace_map = np.concatenate(
            (
                yz_face_sources,
                zx_face_sources,
                xy_face_sources,
                x_axis_y_sources,
                x_axis_z_sources,
                x_axis_yz_sources,
                y_axis_x_sources,
                y_axis_z_sources,
                y_axis_xz_sources,
                z_axis_y_sources,
                z_axis_x_sources,
                z_axis_xy_sources,
                corner_sources,
            )
        )
        return replace_map

    def _create_periodic_mesh(self, source_mesh):
        tdim = source_mesh.topology.dim
        run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT cell-facet connectivity construction",
            lambda: source_mesh.topology.create_connectivity(tdim, tdim - 1),
        )
        run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT facet-cell connectivity construction",
            lambda: source_mesh.topology.create_connectivity(tdim - 1, tdim),
        )
        left_vertices = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT periodic source-boundary location",
            lambda: dolfinx.mesh.locate_entities_boundary(
                source_mesh,
                0,
                self._indicator,
            ),
        )

        def prepare_periodic_target_points():
            left_midpoints = dolfinx.mesh.compute_midpoints(
                source_mesh,
                0,
                left_vertices,
            )
            right_midpoints = self._mapping_function(left_midpoints.T).T
            return left_midpoints, right_midpoints

        left_midpoints, right_midpoints = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT periodic target-point setup",
            prepare_periodic_target_points,
        )
        right_vertices = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT periodic target-boundary location",
            lambda: dolfinx.mesh.locate_entities_boundary(
                source_mesh,
                0,
                self._indicator_dual,
            ),
        )

        def prepare_periodic_vertex_map():
            local_right_vertices = np.asarray(right_vertices, dtype=np.int32)
            bounding_box_tree = dolfinx.geometry.bb_tree(
                source_mesh,
                0,
                entities=local_right_vertices,
            )
            mid_tree = dolfinx.geometry.create_midpoint_tree(
                source_mesh,
                0,
                local_right_vertices,
            )
            closest_vertex = dolfinx.geometry.compute_closest_entity(
                bounding_box_tree,
                mid_tree,
                source_mesh,
                right_midpoints,
            )
            index_map = source_mesh.topology.index_map(0)
            num_vertices_local = index_map.size_local + index_map.num_ghosts
            keep_vertices = np.ones(num_vertices_local, dtype=np.bool_)
            keep_vertices[local_right_vertices] = False
            return closest_vertex, np.flatnonzero(keep_vertices)

        closest_vertex, new_vertices = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT periodic vertex-map setup",
            prepare_periodic_vertex_map,
        )
        new_vertex_map, sub_to_parent = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT periodic sub-index-map construction",
            lambda: dolfinx.cpp.common.create_sub_index_map(
                source_mesh.topology.index_map(0),
                new_vertices,
                allow_owner_change=True,
            ),
        )

        def prepare_periodic_topology():
            index_map = source_mesh.topology.index_map(0)
            num_vertices_local = index_map.size_local + index_map.num_ghosts
            parent_to_sub = np.full(num_vertices_local, -1, dtype=np.int32)
            parent_to_sub[sub_to_parent] = np.arange(
                sub_to_parent.size,
                dtype=np.int32,
            )
            replace_map = np.arange(num_vertices_local, dtype=np.int32)
            replace_map[closest_vertex] = self._mapping_function_inv(
                left_midpoints.T,
                left_vertices,
            )
            c_to_v = source_mesh.topology.connectivity(tdim, 0)
            new_c = parent_to_sub[replace_map[c_to_v.array]]
            new_o = c_to_v.offsets.copy()
            new_c_to_v = dolfinx.graph.adjacencylist(new_c, new_o)
            return dolfinx.cpp.mesh.Topology(
                source_mesh.topology.cell_type,
                new_vertex_map,
                source_mesh.topology.index_map(tdim),
                new_c_to_v._cpp_object,
                np.ascontiguousarray(
                    source_mesh.topology.original_cell_index,
                    dtype=np.int64,
                ),
            )

        topology = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT periodic topology construction",
            prepare_periodic_topology,
        )

        def prepare_periodic_geometry():
            ufl_domain = source_mesh.ufl_domain()
            if ufl_domain is None:
                raise RuntimeError("The source mesh has no UFL domain.")
            coordinate_element = dolfinx.fem.coordinate_element(
                ufl_domain.ufl_coordinate_element().basix_element
            )
            geometry = dolfinx.mesh.create_geometry(
                source_mesh.geometry.index_map(),
                np.ascontiguousarray(
                    source_mesh.geometry.dofmap,
                    dtype=np.int32,
                ),
                coordinate_element._cpp_object,
                np.ascontiguousarray(
                    source_mesh.geometry.x[:, : source_mesh.geometry.dim],
                    dtype=source_mesh.geometry.x.dtype,
                ),
                np.ascontiguousarray(
                    source_mesh.geometry.input_global_indices,
                    dtype=np.int64,
                ),
            )
            return ufl_domain, geometry

        ufl_domain, geometry = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT periodic geometry construction",
            prepare_periodic_geometry,
        )

        def create_cpp_periodic_mesh():
            if source_mesh.geometry.x.dtype == np.float64:
                return dolfinx.cpp.mesh.Mesh_float64(
                    source_mesh.comm,
                    topology,
                    geometry._cpp_object,
                )
            if source_mesh.geometry.x.dtype == np.float32:
                return dolfinx.cpp.mesh.Mesh_float32(
                    source_mesh.comm,
                    topology,
                    geometry._cpp_object,
                )
            raise RuntimeError(
                f"Unsupported dtype for mesh {source_mesh.geometry.x.dtype}"
            )

        cpp_mesh = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT periodic C++ mesh construction",
            create_cpp_periodic_mesh,
        )
        return run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT periodic Python mesh construction",
            lambda: dolfinx.mesh.Mesh(
                cpp_mesh,
                domain=ufl.Mesh(ufl_domain.ufl_coordinate_element()),
            ),
        )

    @staticmethod
    def _interpolate_to_mesh(data, grid_spacing, nx, ny, nz):
        def interpolator(x):
            return data[
                (np.round(x[0] / grid_spacing[0]).astype(int)) % nx,
                (np.round(x[1] / grid_spacing[1]).astype(int)) % ny,
                (np.round(x[2] / grid_spacing[2]).astype(int)) % nz,
            ]

        return interpolator

    def _assemble_global_field(self, phase_name, local_field):
        def prepare_global_field():
            global_field = np.zeros((self.nx, self.ny, self.nz), dtype=np.float64)
            global_field[
                self.l1.value : self.e1.value,
                self.l2.value : self.e2.value,
                self.l3.value : self.e3.value,
            ] = local_field
            return global_field

        global_field = run_local_phase(
            self.mpi_comm,
            phase_name,
            prepare_global_field,
        )
        self.mpi_comm.Allreduce(MPI.IN_PLACE, global_field, op=MPI.SUM)
        return global_field

    def solve(self) -> None:
        """Perform one finite-element solve using the persistent mesh and solver state."""

        def prepare_cp2k_array_views():
            if not self._initialized:
                raise RuntimeError(
                    "Finite-element solve called before initialization completed."
                )

            return pointer_views(
                self.lib,
                (self.ln1.value, self.ln2.value, self.ln3.value),
            )

        (
            solute_charge_density_local,
            relative_permittivity_parallel_local,
            relative_permittivity_perpendicular_local,
            electrostatic_potential_local,
        ) = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT CP2K array view setup",
            prepare_cp2k_array_views,
        )
        interpolation_start = time.time()
        solute_charge_density_global = self._assemble_global_field(
            "AICS FINITE_ELEMENT charge-density assembly",
            solute_charge_density_local,
        )
        run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT charge-density interpolation",
            lambda: self.solute_charge_density.interpolate(
                self._interpolate_to_mesh(
                    solute_charge_density_global,
                    self.grid_spacing,
                    self.nx,
                    self.ny,
                    self.nz,
                )
            ),
        )
        del solute_charge_density_global

        relative_permittivity_parallel_global = self._assemble_global_field(
            "AICS FINITE_ELEMENT parallel-permittivity assembly",
            relative_permittivity_parallel_local,
        )
        run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT parallel-permittivity interpolation",
            lambda: self.relative_permittivity_parallel.interpolate(
                self._interpolate_to_mesh(
                    relative_permittivity_parallel_global,
                    self.grid_spacing,
                    self.nx,
                    self.ny,
                    self.nz,
                )
            ),
        )
        del relative_permittivity_parallel_global

        relative_permittivity_perpendicular_global = self._assemble_global_field(
            "AICS FINITE_ELEMENT perpendicular-permittivity assembly",
            relative_permittivity_perpendicular_local,
        )
        run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT perpendicular-permittivity interpolation",
            lambda: self.relative_permittivity_perpendicular.interpolate(
                self._interpolate_to_mesh(
                    relative_permittivity_perpendicular_global,
                    self.grid_spacing,
                    self.nx,
                    self.ny,
                    self.nz,
                )
            ),
        )
        del relative_permittivity_perpendicular_global

        interpolation_end = time.time()
        self.last_interpolation_time = interpolation_end - interpolation_start

        if self.problem is None:

            def prepare_weak_form():
                phi = TrialFunction(self.function_space)
                v = TestFunction(self.function_space)
                # The parallel coefficient applies to both directions in the interface plane.
                coefficients = [self.relative_permittivity_parallel] * 3
                coefficients[self.normal_axis] = (
                    self.relative_permittivity_perpendicular
                )
                a = (
                    sum(
                        coefficients[axis] * grad(phi)[axis] * grad(v)[axis]
                        for axis in range(3)
                    )
                    * dx
                )
                L = 4.0 * np.pi * self.solute_charge_density * v * dx
                return a, L

            a, L = run_local_phase(
                self.mpi_comm,
                "AICS FINITE_ELEMENT weak-form setup",
                prepare_weak_form,
            )
            problem = run_local_phase(
                self.mpi_comm,
                "AICS FINITE_ELEMENT linear-problem construction",
                lambda: dolfinx.fem.petsc.LinearProblem(
                    a,
                    L,
                    petsc_options_prefix="aics_poisson_",
                    petsc_options=self.petsc_options,
                ),
            )
            null_space = run_local_phase(
                self.mpi_comm,
                "AICS FINITE_ELEMENT null-space construction",
                lambda: PETSc.NullSpace().create(constant=True),
            )
            run_local_phase(
                self.mpi_comm,
                "AICS FINITE_ELEMENT null-space attachment",
                lambda: problem.A.setNullSpace(null_space),
            )
            self.problem = problem
            self.null_space = null_space

        solve_start = time.time()
        try:
            phi_solution = self.problem.solve()
        except PETSc.Error as exc:
            error = RuntimeError(
                _ksp_failure_message(
                    "FINITE_ELEMENT",
                    self.problem.solver,
                    self.relative_tolerance,
                    self.absolute_tolerance,
                    self.max_iterations,
                )
            )
            error.add_note(f"Original PETSc error:\n{exc}")
            raise error from None
        solve_end = time.time()

        def prepare_solution_reconstruction():
            self.phi_solution = phi_solution
            self.last_iterations = self.problem.solver.getIterationNumber()
            self.last_residual = self.problem.solver.getResidualNorm()
            self.last_convergence_reason = self.problem.solver.getConvergedReason()
            self.last_linear_solve_time = solve_end - solve_start
            values_1d = phi_solution.eval(
                self._reconstruction_points,
                self._reconstruction_cells,
            )
            return values_1d, values_1d.size

        values_1d, local_value_size = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT local solution reconstruction",
            prepare_solution_reconstruction,
        )

        recv_counts_values = self.mpi_comm.gather(local_value_size, root=0)

        def prepare_solution_gather():
            if self.rank == 0:
                global_size = self.periodic_mesh.geometry.index_map().size_global
                global_values_1d = np.zeros(global_size, dtype=values_1d.dtype)
                displs_values = np.cumsum([0] + recv_counts_values[:-1])
            else:
                global_values_1d = None
                displs_values = None
            return global_values_1d, displs_values

        global_values_1d, displs_values = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT solution-gather setup",
            prepare_solution_gather,
        )

        self.mpi_comm.Gatherv(
            sendbuf=values_1d,
            recvbuf=(global_values_1d, recv_counts_values, displs_values, MPI.DOUBLE),
            root=0,
        )

        def prepare_nodal_reduction():
            values_3d = np.zeros(
                (self.nx + 1, self.ny + 1, self.nz + 1), dtype=np.float64
            )
            if self.rank == 0:
                sorted_values = global_values_1d[self._reconstruction_sorted_indices]
                values_3d = sorted_values.reshape(
                    (self.nx + 1, self.ny + 1, self.nz + 1), order="C"
                )
            return values_3d

        values_3d = run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT nodal-reduction setup",
            prepare_nodal_reduction,
        )
        self.mpi_comm.Allreduce(MPI.IN_PLACE, values_3d, op=MPI.SUM)

        def store_local_potential():
            electrostatic_potential_local[
                0 : self.ln1.value, 0 : self.ln2.value, 0 : self.ln3.value
            ] = values_3d[
                self.l1.value : self.e1.value,
                self.l2.value : self.e2.value,
                self.l3.value : self.e3.value,
            ]

        run_local_phase(
            self.mpi_comm,
            "AICS FINITE_ELEMENT CP2K potential storage",
            store_local_potential,
        )
        del values_3d
