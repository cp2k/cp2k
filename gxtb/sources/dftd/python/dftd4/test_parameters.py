# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.


from dftd4.parameters import get_all_damping_params, get_damping_param
from pytest import approx


def get_data_file_name() -> str:
    """Make sure we can still test without installing"""
    from os.path import dirname, exists, join

    from dftd4.parameters import get_data_file_name as _get_data_file_name

    data_file = join(dirname(__file__), "..", "..", "assets", "parameters.toml")
    return data_file if exists(data_file) else _get_data_file_name()


def test_get_b3lyp() -> None:
    expected = {
        "s6": 1.0,
        "s9": 1.0,
        "alp": 16.0,
        "s8": 2.02929367,
        "a1": 0.40868035,
        "a2": 4.53807137,
    }
    actual = get_damping_param("b3lyp", data_file=get_data_file_name())

    for key in expected.keys():
        assert approx(actual[key]) == expected[key]


def test_get_b2plyp() -> None:
    expected = {
        "s6": 0.64,
        "s9": 1.0,
        "alp": 16.0,
        "s8": 1.16888646,
        "a1": 0.44154604,
        "a2": 4.73114642,
    }
    actual = get_damping_param("b2plyp", data_file=get_data_file_name())

    for key in expected.keys():
        assert approx(actual[key]) == expected[key]


def test_get_pw6b95() -> None:
    expected = {
        "s6": 1.0,
        "s9": 1.0,
        "alp": 16.0,
        "s8": -0.31629935,
        "a1": 0.03999357,
        "a2": 5.83690254,
    }
    actual = get_damping_param(
        "pw6b95",
        data_file=get_data_file_name(),
        defaults=["bj-eeq-two", "bj-eeq-mbd"],
    )

    for key in expected.keys():
        assert approx(actual[key]) == expected[key]


def test_all_parameters() -> None:
    params = get_all_damping_params()

    assert "b3lyp" in params
    assert "b2plyp" in params
    assert "pw6b95" in params
