# SPDX-License-Identifier: GPL-2.0-or-later

from copy import deepcopy
import os
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("aiida")
pytest.importorskip("aiida_cp2k")
pytest.importorskip("aiida_common_workflows")
Atoms = pytest.importorskip("ase").Atoms

from aiida import engine, orm
from aiida.plugins import WorkflowFactory
from aiida.tools.pytest_fixtures import aiida_config_factory, aiida_profile_factory

from cp2k.aiida import build_relax_builder


@pytest.fixture(scope="module")
def workflow_profile(tmp_path_factory, aiida_config_factory, aiida_profile_factory):
    # Isolated SQLite test storage, no daemon/broker, no user-profile mutation.
    with aiida_config_factory(tmp_path_factory.mktemp("workflow-profile")) as config:
        with aiida_profile_factory(config, name="cp2k-workflow-tests") as profile:
            yield profile


@pytest.fixture
def workflow_inputs(workflow_profile, tmp_path):
    computer = orm.Computer(
        label=f"cp2k-{tmp_path.name}",
        hostname="localhost",
        transport_type="core.local",
        scheduler_type="core.direct",
        workdir=str(tmp_path / "jobs"),
    ).store()
    computer.configure()
    computer.set_minimum_job_poll_interval(0.1)
    computer.set_default_mpiprocs_per_machine(1)
    code = orm.InstalledCode(
        computer=computer,
        filepath_executable=os.environ.get("CP2K_TEST_EXECUTABLE", "/bin/true"),
        default_calc_job_plugin="cp2k",
        label="cp2k-tests",
    ).store()
    # A periodic large box keeps the common protocol at Gamma. This is a
    # periodic H2 smoke test, not an isolated-molecule convergence benchmark.
    structure = orm.StructureData(
        ase=Atoms("H2", positions=[[4.6, 5, 5], [5.4, 5, 5]], cell=[10] * 3, pbc=True)
    )
    return {
        "structure": structure,
        "code": code,
        "electronic_type": "insulator",
        "protocol": "fast",
        "options": {
            "resources": {"num_machines": 1, "num_mpiprocs_per_machine": 1},
            "withmpi": False,
            "max_wallclock_seconds": 600,
            "environment_variables": {
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
            },
        },
    }


@pytest.mark.parametrize(
    "relax_type,run_type",
    [
        ("none", "ENERGY_FORCE"),
        ("positions", "GEO_OPT"),
        ("positions_cell", "CELL_OPT"),
    ],
)
def test_builder_uses_existing_workflow(workflow_inputs, relax_type, run_type):
    saved = deepcopy(workflow_inputs["options"])
    nodes_before = orm.QueryBuilder().append(orm.Node).count()
    builder = build_relax_builder(**workflow_inputs, relax_type=relax_type)
    assert builder._process_class is WorkflowFactory("common_workflows.relax.cp2k")
    assert builder.cp2k.parameters["GLOBAL"]["RUN_TYPE"] == run_type
    assert "TYPE" not in builder.cp2k.parameters["MOTION"]["CELL_OPT"]
    assert builder.cp2k.parameters["MOTION"]["GEO_OPT"]["TYPE"] == "MINIMIZATION"
    printing = builder.cp2k.parameters["MOTION"]["PRINT"]
    assert printing["TRAJECTORY"]["FORMAT"] == "XYZ"
    assert printing["CELL"]["_"] == "ON"
    assert printing["CELL"]["EACH"] == printing["TRAJECTORY"]["EACH"]
    assert printing["CELL"]["ADD_LAST"] == printing["TRAJECTORY"]["ADD_LAST"]
    assert builder.cp2k.structure is workflow_inputs["structure"]
    assert builder.cp2k.code is workflow_inputs["code"]
    assert not builder.cp2k.parameters.is_stored
    assert not builder.cp2k.structure.is_stored
    assert builder.cp2k.file["potential"].filename == "GTH_POTENTIALS"
    assert "OT" in builder.cp2k.parameters["FORCE_EVAL"]["DFT"]["SCF"]
    assert builder.cp2k.metadata.options.parser_name == "cp2k_advanced_parser"
    assert workflow_inputs["options"] == saved
    assert orm.QueryBuilder().append(orm.Node).count() == nodes_before


def test_metal_spin_and_thresholds(workflow_inputs):
    workflow_inputs["electronic_type"] = "metal"
    builder = build_relax_builder(
        **workflow_inputs,
        relax_type="positions_cell",
        spin_type="collinear",
        magnetization_per_site=[1.0, -1.0],
        threshold_forces=0.03,
        threshold_stress=0.001,
    )
    params = builder.cp2k.parameters.get_dict()
    dft = params["FORCE_EVAL"]["DFT"]
    assert dft["UKS"]
    assert dft["MULTIPLICITY"] == 1
    assert "SMEAR" in dft["SCF"]
    assert params["MOTION"]["CELL_OPT"]["MAX_FORCE"] == "[eV/angstrom] 0.03"
    assert params["MOTION"]["CELL_OPT"]["PRESSURE_TOLERANCE"] == "[GPa] 0.16021766208"
    assert not workflow_inputs["structure"].is_stored


@pytest.mark.parametrize(
    "kwargs,message",
    [
        ({"magnetization_per_site": [1, -1]}, "collinear"),
        ({"magnetization_per_site": [1], "spin_type": "collinear"}, "one finite value"),
        (
            {"magnetization_per_site": [1, float("nan")], "spin_type": "collinear"},
            "one finite value",
        ),
        ({"threshold_forces": 0.03, "relax_type": "none"}, "requires a relaxation"),
        ({"threshold_stress": 0.01}, "positions_cell"),
        ({"threshold_forces": -1.0}, "finite and positive"),
        ({"threshold_forces": float("nan")}, "finite and positive"),
        ({"relax_type": "positions_volume"}, "not a valid choice"),
        ({"spin_type": "non_collinear"}, "not a valid choice"),
    ],
)
def test_invalid_physical_options(workflow_inputs, kwargs, message):
    with pytest.raises((ValueError, TypeError), match=message):
        build_relax_builder(**workflow_inputs, **kwargs)


def test_nonperiodic_structure_is_not_silently_made_periodic(workflow_inputs):
    workflow_inputs["structure"].pbc = False
    with pytest.raises(ValueError, match="fully periodic"):
        build_relax_builder(**workflow_inputs)


@pytest.mark.integration
@pytest.mark.parametrize("relax_type", ["none", "positions", "positions_cell"])
def test_real_common_relaxation(workflow_inputs, relax_type):
    executable = os.environ.get("CP2K_TEST_EXECUTABLE")
    if not os.environ.get("CP2K_TEST_AIIDA") or not executable:
        pytest.skip(
            "Set CP2K_TEST_AIIDA=1 and CP2K_TEST_EXECUTABLE for a real AiiDA CP2K job"
        )
    assert Path(executable).is_file()
    builder = build_relax_builder(**workflow_inputs, relax_type=relax_type)
    results, node = engine.run_get_node(builder)
    assert node.is_finished_ok, (node.exit_status, node.exit_message)
    assert node.called and all(child.is_finished_ok for child in node.called)
    assert -40 < results["total_energy"].value < -20  # eV, not hartree
    forces = results["forces"].get_array("forces")
    assert forces.shape == (2, 3) and np.isfinite(forces).all()
    calcjobs = [
        child for child in node.called_descendants if isinstance(child, orm.CalcJobNode)
    ]
    assert calcjobs and all(child.is_finished_ok for child in calcjobs)
    assert calcjobs[-1].inputs.file.potential.is_stored
    assert calcjobs[-1].outputs.retrieved.list_object_names()
    if relax_type == "none":
        return
    assert results["relaxed_structure"].is_stored
    assert len(results["relaxed_structure"].sites) == 2
    assert np.max(np.abs(forces)) < 0.05  # common output: eV/angstrom
    trajectory = calcjobs[-1].outputs.output_trajectory
    assert trajectory.numsteps > 1
    # The .cell trajectory stores ten decimal places, the restart stores more.
    np.testing.assert_allclose(
        trajectory.get_cells()[-1], results["relaxed_structure"].cell, rtol=0, atol=1e-9
    )
    np.testing.assert_allclose(
        trajectory.get_positions()[-1],
        results["relaxed_structure"].get_ase().positions,
        atol=1e-6,
    )
