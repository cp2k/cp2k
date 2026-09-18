# SPDX-License-Identifier: GPL-2.0-or-later

import numpy as np
import pytest


def test_stress_lifecycle(runtime, fake_library, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with runtime.create_force_env({"GLOBAL": {}}) as env:
        with pytest.raises(RuntimeError, match="stress=True"):
            _ = env.stress
        result = env.calculate(stress=True)
        np.testing.assert_array_equal(result.stress, fake_library.stress)
        np.testing.assert_allclose(
            result.virial, result.stress * np.linalg.det(env.cell)
        )
        env.positions = env.positions
        with pytest.raises(RuntimeError, match="stress=True"):
            _ = env.virial
        assert env.calculate(forces=False, stress=True).forces is None
        env.calculate(forces=False)
        with pytest.raises(RuntimeError, match="stress=True"):
            _ = env.stress
        fake_library.stress_available = 0
        with pytest.raises(RuntimeError, match="STRESS_TENSOR"):
            env.calculate(stress=True)
        with pytest.raises(RuntimeError):
            _ = env.potential_energy


def test_older_library(runtime, fake_library, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    del fake_library.cp2k_get_stress_tensor
    with runtime.create_force_env({"GLOBAL": {}}) as env:
        assert env.calculate().energy == -1
        with pytest.raises(RuntimeError, match="lacks"):
            env.calculate(stress=True)
