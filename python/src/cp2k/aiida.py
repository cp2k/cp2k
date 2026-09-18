"""Optional bridge to the existing AiiDA common CP2K relaxation workflow."""

# SPDX-License-Identifier: GPL-2.0-or-later

from copy import deepcopy
import math


def build_relax_builder(
    structure,
    *,
    code,
    options,
    electronic_type,
    protocol="moderate",
    relax_type="positions",
    spin_type="none",
    magnetization_per_site=None,
    threshold_forces=None,
    threshold_stress=None,
    reference_workchain=None,
):
    """Build, but do not submit, an AiiDA common CP2K relaxation workflow.

    Load an AiiDA profile first. ``structure`` is a StructureData node, ``code``
    a configured CP2K Code, and ``options`` its scheduler options. Electronic
    type (``metal`` or ``insulator``) is an explicit physical choice. Relaxation
    types are ``none``, ``positions`` and ``positions_cell``. Protocols and all
    CP2K input generation belong to aiida-common-workflows, not this interface.
    Thresholds are in eV/angstrom and eV/angstrom**3. This function neither
    stores nor submits the builder, and never loads libcp2k in the AiiDA daemon.
    The supplied common protocols describe neutral, fully periodic systems.
    Use a custom cp2k.base builder for different boundary conditions/charge.
    """
    from aiida import orm
    from aiida.plugins import WorkflowFactory
    from aiida_common_workflows.common import ElectronicType, RelaxType, SpinType

    if not isinstance(structure, orm.StructureData):
        raise TypeError("structure must be an AiiDA StructureData node")
    if not all(structure.pbc):
        raise ValueError(
            "The common CP2K protocols require fully periodic structures; for molecules "
            "or slabs, use cp2k.base with explicit cell/Poisson boundary conditions"
        )
    if not isinstance(options, dict):
        raise TypeError("options must be a scheduler-options dictionary")
    relax_type = RelaxType(relax_type)
    electronic_type = ElectronicType(electronic_type)
    spin_type = SpinType(spin_type)
    if magnetization_per_site is not None and spin_type != SpinType.COLLINEAR:
        raise ValueError("magnetization_per_site requires collinear spin")
    if magnetization_per_site is not None:
        if len(magnetization_per_site) != len(structure.sites) or not all(
            math.isfinite(value) for value in magnetization_per_site
        ):
            raise ValueError(
                "magnetization_per_site must contain one finite value per site"
            )
    if threshold_forces is not None and relax_type == RelaxType.NONE:
        raise ValueError("threshold_forces requires a relaxation")
    if threshold_stress is not None and relax_type != RelaxType.POSITIONS_CELL:
        raise ValueError("threshold_stress requires positions_cell relaxation")
    for name, value in (
        ("threshold_forces", threshold_forces),
        ("threshold_stress", threshold_stress),
    ):
        if value is not None and (not math.isfinite(value) or value <= 0):
            raise ValueError(f"{name} must be finite and positive")
    kwargs = {
        "structure": structure,
        "engines": {"relax": {"code": code, "options": deepcopy(options)}},
        "protocol": protocol,
        "relax_type": relax_type,
        "electronic_type": electronic_type,
        "spin_type": spin_type,
    }
    for key, value in (
        ("magnetization_per_site", magnetization_per_site),
        ("threshold_forces", threshold_forces),
        ("threshold_stress", threshold_stress),
        ("reference_workchain", reference_workchain),
    ):
        if value is not None:
            kwargs[key] = value
    workflow = WorkflowFactory("common_workflows.relax.cp2k")
    builder = workflow.get_input_generator().get_builder(**kwargs)
    parameters = builder.cp2k.parameters.get_dict()
    cell_opt = parameters.get("MOTION", {}).get("CELL_OPT", {})
    if cell_opt.get("TYPE") == "DIRECT_CELL_OPT":
        # Current CP2K always uses DIRECT_CELL_OPT and rejects the old TYPE
        # keyword, including in the unused CELL_OPT section of a GEO_OPT run.
        # Omitting this former default also preserves older CP2K semantics.
        del cell_opt["TYPE"]
    # aiida-cp2k's advanced parser reads XYZ, not the common protocol's DCD.
    # Keep cells at the same cadence for a meaningful variable-cell trajectory.
    printing = parameters.setdefault("MOTION", {}).setdefault("PRINT", {})
    trajectory = printing.setdefault("TRAJECTORY", {})
    trajectory["FORMAT"] = "XYZ"
    trajectory["ADD_LAST"] = "NUMERIC"
    printing["CELL"] = {
        "_": "ON",
        "ADD_LAST": "NUMERIC",
        "EACH": deepcopy(trajectory.get("EACH", {})),
    }
    builder.cp2k.parameters = orm.Dict(dict=parameters)
    return builder
