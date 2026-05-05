# This file is part of s-dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.


from pytest import approx
from dftd3.parameters import get_damping_param, get_all_damping_params


def get_data_file_name() -> str:
    """Make sure we can still test without installing"""
    from os.path import join, dirname, exists
    from dftd3.parameters import get_data_file_name as _get_data_file_name

    data_file = join(dirname(__file__), "..", "..", "assets", "parameters.toml")
    return data_file if exists(data_file) else _get_data_file_name()


def test_get_b3lyp():

    expected = {
        "s6": 1.0,
        "s9": 1.0,
        "alp": 14.0,
        "a1": 0.3981,
        "s8": 1.9889,
        "a2": 4.4211,
    }
    actual = get_damping_param("b3lyp", data_file=get_data_file_name())

    for key in expected.keys():
        assert approx(actual[key]) == expected[key]


def test_get_m11l():

    expected = {
        "s6": 1.0,
        "s9": 1.0,
        "alp": 14.0,
        "rs8": 1.0,
        "s8": 1.1129,
        "rs6": 2.3933,
    }
    actual = get_damping_param("m11l", data_file=get_data_file_name())

    for key in expected.keys():
        assert approx(actual[key]) == expected[key]


def test_get_pw6b95():

    expected = {
        "s6": 1.0,
        "s9": 1.0,
        "alp": 14.0,
        "a1": 0.2076,
        "s8": 0.7257,
        "a2": 6.3750,
    }
    actual = get_damping_param(
        "pw6b95", data_file=get_data_file_name(), defaults=["bj"]
    )

    for key in expected.keys():
        assert approx(actual[key]) == expected[key]


def test_all_parameters():

    params = get_all_damping_params()

    assert "b3lyp" in params
    assert "b2plyp" in params
    assert "pw6b95" in params
