# SPDX-License-Identifier: GPL-2.0-or-later

import pickle
from types import SimpleNamespace

import pytest

from cp2k.openmm import _Computation
from cp2k.lammps import ExternalForce


def test_openmm_ownership(runtime, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with runtime.create_force_env({"GLOBAL": {}}) as env:
        callback = _Computation(env, True)
        with pytest.raises(TypeError, match="cannot be serialized"):
            pickle.dumps(callback)
        with monkeypatch.context() as patch:
            patch.setattr(
                runtime,
                "_mpi",
                SimpleNamespace(COMM_NULL=None, Is_finalized=lambda: False),
            )
            patch.setattr(runtime, "_comm", SimpleNamespace(size=2))
            with pytest.raises(ValueError, match="single-rank"):
                _Computation(env, True)
    with pytest.raises(RuntimeError, match="closed"):
        callback(None)


@pytest.mark.parametrize("fix_id", ["", "a b", "x;run"])
def test_lammps_fix_validation(runtime, tmp_path, monkeypatch, fix_id):
    monkeypatch.chdir(tmp_path)
    with runtime.create_force_env({"GLOBAL": {}}) as env:
        with pytest.raises(ValueError, match="fix_id"):
            ExternalForce(None, env, fix_id=fix_id)
