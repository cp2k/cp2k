"""LAMMPS fix external callback using a persistent, caller-owned libcp2k."""

# SPDX-License-Identifier: GPL-2.0-or-later

import re

import numpy as np

from ._units import BOHR_TO_ANGSTROM, HARTREE_TO_EV, HARTREE_TO_KCALMOL


class ExternalForce:
    """Attach the whole CP2K potential to LAMMPS via ``fix external``.

    Requires identical atoms and periodicity in both inputs. ``atom_ids`` maps
    CP2K's atom order to LAMMPS tags; default 1..N. Atom migration and sorting are
    supported. Units must be ``metal`` or ``real``; cells must be 3D restricted
    triclinic/orthogonal. ``stress=True`` requires CP2K STRESS_TENSOR and enables
    cell changes/NPT. No per-atom energy/stress or automatic QM/MM partitioning
    is supplied. Existing interaction forces are additive: avoid double counting.

    For MPI, initialize mpi4py first and pass the same communicator to LAMMPS
    and CP2K, using the same MPI implementation. All ranks execute every call.
    Use command()/run() below (or call check() after EVERY native command).
    ctypes cannot propagate callback exceptions: a failed callback requests a
    LAMMPS timeout, and check() raises the saved error. Discard any output from
    a failed run. Native CP2K aborts still terminate the process.

    close() removes only the owned fix. Close before LAMMPS/environment/runtime;
    do not clear LAMMPS or unfix externally while attached. Native restart files
    do not store this callback: recreate the adapter after reading a restart.
    """

    def __init__(
        self, lammps, environment, *, atom_ids=None, stress=True, fix_id="cp2k"
    ):
        environment._check()
        if not re.fullmatch(r"[A-Za-z0-9_]+", fix_id):
            raise ValueError(
                "fix_id must contain only letters, digits, and underscores"
            )
        if lammps.has_id("fix", fix_id):
            raise ValueError(f"LAMMPS fix {fix_id!r} already exists")
        if lammps.extract_setting("dimension") != 3:
            raise ValueError("CP2K requires a three-dimensional LAMMPS box")
        self.lammps, self.environment = lammps, environment
        self.fix_id, self.stress = fix_id, bool(stress)
        units = lammps.extract_global("units")
        if units not in ("real", "metal"):
            raise ValueError("The CP2K callback supports LAMMPS real or metal units")
        self.units = units
        self._energy_factor = HARTREE_TO_EV if units == "metal" else HARTREE_TO_KCALMOL
        self._comm = lammps.get_mpi_comm()
        cp_comm = environment._runtime._comm
        if self._comm is not None:
            from mpi4py import MPI

            if cp_comm is None or MPI.Comm.Compare(self._comm, cp_comm) not in (
                MPI.IDENT,
                MPI.CONGRUENT,
            ):
                raise ValueError(
                    "Initialize CP2K and LAMMPS with congruent mpi4py communicators"
                )
        elif lammps.extract_setting("world_size") != 1 or (
            cp_comm is not None and cp_comm.size != 1
        ):
            raise ValueError(
                "Parallel LAMMPS callbacks require mpi4py and matching communicators"
            )
        n = environment.nparticle
        if n != environment.natom or n != int(lammps.get_natoms()):
            raise ValueError(
                "LAMMPS and CP2K must have the same atoms, without shell particles"
            )
        ids = np.arange(1, n + 1) if atom_ids is None else np.asarray(atom_ids)
        if (
            ids.shape != (n,)
            or ids.dtype.kind not in "iu"
            or np.any(ids <= 0)
            or len(np.unique(ids)) != n
        ):
            raise ValueError("atom_ids must contain N distinct positive integer tags")
        self._ids = ids.copy()
        self._indices = {int(tag): index for index, tag in enumerate(ids)}
        self._initial_cell, self._periodic, _ = self._box()
        self._error = None
        self._closed = False
        lammps.command(f"fix {fix_id} all external pf/callback 1 1")
        try:
            lammps.set_fix_external_callback(fix_id, self._callback)
            lammps.command(f"fix_modify {fix_id} energy yes virial yes")
        except Exception:
            lammps.command(f"unfix {fix_id}")
            raise

    def _box(self):
        lo, hi, xy, yz, xz, periodic, _ = self.lammps.extract_box()
        lengths = np.asarray(hi) - lo
        cell = np.array([[lengths[0], 0, 0], [xy, lengths[1], 0], [xz, yz, lengths[2]]])
        return cell, tuple(periodic), np.asarray(lo)

    def _gather(self, value):
        return [value] if self._comm is None else self._comm.allgather(value)

    def _callback(self, caller, step, nlocal, tags, positions, forces):
        # Never let an exception escape through ctypes (it would be swallowed).
        forces[:] = 0
        try:
            self._collective_check()
            pieces = self._gather(
                (np.array(tags, copy=True).reshape(-1), np.array(positions, copy=True))
            )
            all_tags = np.concatenate([piece[0] for piece in pieces])
            all_positions = np.concatenate([piece[1] for piece in pieces])
            if not np.array_equal(np.sort(all_tags), np.sort(self._ids)):
                raise ValueError("LAMMPS atom tags changed or do not match atom_ids")
            cell, periodic, origin = self._box()
            if (
                periodic != self._periodic
                or self.lammps.extract_global("units") != self.units
            ):
                raise ValueError(
                    "LAMMPS periodicity/units changed while CP2K is attached"
                )
            if not self.stress and not np.allclose(
                cell, self._initial_cell, rtol=0, atol=1e-12
            ):
                raise ValueError(
                    "Cell changes require stress=True and CP2K STRESS_TENSOR"
                )
            ordered = np.empty_like(all_positions)
            ordered[[self._indices[int(tag)] for tag in all_tags]] = (
                all_positions - origin
            )
            env = self.environment
            env.cell = cell / BOHR_TO_ANGSTROM
            env.positions = ordered / BOHR_TO_ANGSTROM
            result = env.calculate(stress=self.stress)
            indices = [self._indices[int(tag)] for tag in tags.reshape(-1)]
            forces[:] = result.forces[indices] * self._energy_factor / BOHR_TO_ANGSTROM
            self.lammps.fix_external_set_energy_global(
                self.fix_id, result.energy * self._energy_factor
            )
            if self.stress:
                # LAMMPS expects extensive, pressure-positive [xx yy zz xy xz yz].
                # Its setter performs MPI normalization; do NOT divide here.
                self.lammps.fix_external_set_virial_global(
                    self.fix_id,
                    result.virial.flat[[0, 4, 8, 1, 2, 5]] * self._energy_factor,
                )
        except Exception as error:
            self._error = error
            forces[:] = 0
            self.lammps.force_timeout()

    def check(self):
        if self._closed:
            raise RuntimeError("The CP2K LAMMPS callback is closed")
        if self._error is not None:
            raise RuntimeError(
                "CP2K LAMMPS callback failed; discard this run"
            ) from self._error
        self.environment._check()

    def command(self, command):
        """Execute one LAMMPS command and propagate saved callback errors."""
        self._collective_check()
        self.lammps.command(command)
        self._collective_check()

    def _collective_check(self):
        local_error = None
        try:
            self.check()
        except Exception as error:
            local_error = str(error)
        errors = self._gather(local_error)
        if any(error is not None for error in errors):
            raise RuntimeError(f"Collective callback validation failed: {errors}")

    def run(self, steps):
        if not isinstance(steps, int) or steps < 0:
            raise ValueError("steps must be a nonnegative integer")
        self.command(f"run {steps}")

    def close(self):
        if not self._closed:
            self.lammps.command(f"unfix {self.fix_id}")
            self.lammps.callback.pop(self.fix_id, None)
            self._closed = True

    def __enter__(self):
        self.check()
        return self

    def __exit__(self, *exc):
        self.close()
