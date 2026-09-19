# SPDX-License-Identifier: GPL-2.0-or-later

import pytest

from cp2k import input_to_string


def test_sections_keywords_and_repetition():
    text = input_to_string(
        {
            "force_eval": {
                "dft": {"uks": True, "basis_set_file_name": [["BASIS1"], ["BASIS2"]]},
                "subsys": {
                    "kind": [{"_": "H", "basis_set": "DZVP"}, {"_": "O"}],
                    "coord": {"_lines": ["H 0 0 0"]},
                },
            }
        }
    )
    assert "    UKS .TRUE.\n" in text
    assert text.count("&KIND") == 2
    assert "&KIND O\n" in text
    assert "BASIS_SET_FILE_NAME BASIS1" in text
    assert "BASIS_SET_FILE_NAME BASIS2" in text
    assert (
        input_to_string({"GLOBAL": {}, "TEST": None}) == "&GLOBAL\n&END GLOBAL\nTEST\n"
    )


@pytest.mark.parametrize(
    "tree",
    [
        {"GLOBAL": {}, "global": {}},
        {"bad key": 1},
        {"TEST": float("nan")},
        {"TEST": "x\n&END"},
        {"TEST": "x\0"},
        {"TEST": []},
        {"_": "oops"},
        {"_lines": "text"},
        {"_lines": ["&END"]},
        {"TEST": [{}, 3]},
        {"TEST": [[1], 3]},
    ],
)
def test_invalid_input(tree):
    with pytest.raises((TypeError, ValueError)):
        input_to_string(tree)


def test_parameters_arrays_and_booleans():
    text = input_to_string(
        {"CELL": {"ABC": (1.0, 2.0, 3.0), "PERIODIC": "NONE"}, "OT": {"_": False}}
    )
    assert "ABC 1.0 2.0 3.0" in text
    assert "&OT .FALSE." in text
