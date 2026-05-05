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

from typing import Optional

import numpy as np
import pytest
from pytest import approx

try:
    from dftd3.qcschema import run_qcschema, qcel_v1, qcel_v2
except ModuleNotFoundError:
    qcel_v1 = None
    qcel_v2 = None

v1_available = pytest.mark.skipif(
    qcel_v1 is None, reason="QCSchema v1 not available for py314+"
)
v2_available = pytest.mark.skipif(
    qcel_v2 is None, reason="QCSchema v2 not available in current QCElemental"
)


def get_molecule(name: str) -> dict:
    if name == "molecule1":
        return {
            "symbols": "C C C C C C I H H H H H S H C H H H".split(" "),
            "geometry": [
                [-1.42754169820131, -1.50508961850828, -1.93430551124333],
                [+1.19860572924150, -1.66299114873979, -2.03189643761298],
                [+2.65876001301880, +0.37736955363609, -1.23426391650599],
                [+1.50963368042358, +2.57230374419743, -0.34128058818180],
                [-1.12092277855371, +2.71045691257517, -0.25246348639234],
                [-2.60071517756218, +0.67879949508239, -1.04550707592673],
                [-2.86169588073340, +5.99660765711210, +1.08394899986031],
                [+2.09930989272956, -3.36144811062374, -2.72237695164263],
                [+2.64405246349916, +4.15317840474646, +0.27856972788526],
                [+4.69864865613751, +0.26922271535391, -1.30274048619151],
                [-4.63786461351839, +0.79856258572808, -0.96906659938432],
                [-2.57447518692275, -3.08132039046931, -2.54875517521577],
                [-5.88211879210329, 11.88491819358157, +2.31866455902233],
                [-8.18022701418703, 10.95619984550779, +1.83940856333092],
                [-5.08172874482867, 12.66714386256482, -0.92419491629867],
                [-3.18311711399702, 13.44626574330220, -0.86977613647871],
                [-5.07177399637298, 10.99164969235585, -2.10739192258756],
                [-6.35955320518616, 14.08073002965080, -1.68204314084441],
            ],
        }
    if name == "molecule2":
        return {
            "symbols": "C C C C N C S H H H H H".split(),
            "geometry": [
                [-2.56745685564671, -0.02509985979910, 0.00000000000000],
                [-1.39177582455797, +2.27696188880014, 0.00000000000000],
                [+1.27784995624894, +2.45107479759386, 0.00000000000000],
                [+2.62801937615793, +0.25927727028120, 0.00000000000000],
                [+1.41097033661123, -1.99890996077412, 0.00000000000000],
                [-1.17186102298849, -2.34220576284180, 0.00000000000000],
                [-2.39505990368378, -5.22635838332362, 0.00000000000000],
                [+2.41961980455457, -3.62158019253045, 0.00000000000000],
                [-2.51744374846065, +3.98181713686746, 0.00000000000000],
                [+2.24269048384775, +4.24389473203647, 0.00000000000000],
                [+4.66488984573956, +0.17907568006409, 0.00000000000000],
                [-4.60044244782237, -0.17794734637413, 0.00000000000000],
            ],
        }
    if name == "counterpoise":
        return {
            "symbols": "Pb H H H H Bi H H H".split(),
            "geometry": [
                [-0.00000020988889, -4.98043478877778, +0.00000000000000],
                [+3.06964045311111, -6.06324400177778, +0.00000000000000],
                [-1.53482054188889, -6.06324400177778, -2.65838526500000],
                [-1.53482054188889, -6.06324400177778, +2.65838526500000],
                [-0.00000020988889, -1.72196703577778, +0.00000000000000],
                [-0.00000020988889, +4.77334244722222, +0.00000000000000],
                [+1.35700257511111, +6.70626379422222, -2.35039772300000],
                [-2.71400388988889, +6.70626379422222, +0.00000000000000],
                [+1.35700257511111, +6.70626379422222, +2.35039772300000],
            ],
            "real": [True] * 5 + [False] * 4,
        }

    raise ValueError(f"Unknown molecule name: {name}")

def get_atomic_input(
    version: int,
    molecule: dict,
    driver: str,
    method: str,
    params_tweaks: Optional[dict] = None,
    level_hint: Optional[str] = None,
    qcel_object: bool = False,
) -> dict:
    if version == 1:
        input_data = {
            "molecule": molecule,
            "driver": driver,
            "model": {
                "method": method,
            },
            "keywords": {},
        }
        if params_tweaks is not None:
            input_data["keywords"]["params_tweaks"] = params_tweaks

        if level_hint is not None:
            input_data["keywords"]["level_hint"] = level_hint

        if qcel_object:
            return qcel_v1.AtomicInput(**input_data)
        return input_data

    if version == 2:
        input_data = {
            "molecule": molecule,
            "specification": {
                "driver": driver,
                "model": {
                    "method": method,
                },
                "keywords": {},
            },
        }
        if params_tweaks is not None:
            input_data["specification"]["keywords"]["params_tweaks"] = params_tweaks
        if level_hint is not None:
            input_data["specification"]["keywords"]["level_hint"] = level_hint

        if qcel_object:
            return qcel_v2.AtomicInput(**input_data)
        return input_data

    raise ValueError(f"Unsupported version: {version}")


@pytest.fixture(params=[pytest.param(1, marks=v1_available), pytest.param(2, marks=v2_available)])
def qcsk_version(request):
    return request.param

@pytest.fixture(params=["D3", "D3ATM"])
def atm(request):
    return request.param == "D3ATM"


@pytest.mark.skipif(
    qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental models"
)
def test_energy_r2scan_d3bj(atm: bool, qcsk_version: int) -> None:
    thr = 1e-9

    atomic_input = get_atomic_input(
        version=qcsk_version,
        molecule=get_molecule("molecule1"),
        driver="energy",
        method="",
        params_tweaks={
            "method": "r2scan",
            "atm": atm,
        },
        level_hint="d3bj",
        qcel_object=True,
    )
    ref = -0.005790963570050724 if atm else -0.005784012374055654

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == ref


@pytest.mark.skipif(
    qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental models"
)
def test_energy_bp_d3zero(atm: bool, qcsk_version: int) -> None:
    thr = 1e-9

    atomic_input = get_atomic_input(
        version=qcsk_version,
        molecule=get_molecule("molecule1"),
        driver="energy",
        method="",
        params_tweaks={
            "s8": 1.683,
            "rs6": 1.139,
            "s9": 1.0 if atm else 0.0,
        },
        level_hint="d3zero",
        qcel_object=True,
    )
    ref = -0.01410721853585842 if atm else -0.014100267345314462

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == ref


@pytest.mark.skipif(
    qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental models"
)
def test_gradient_b97d_d3bj(atm: bool, qcsk_version: int) -> None:
    thr = 1e-9

    atomic_input = get_atomic_input(
        version=qcsk_version,
        molecule=get_molecule("molecule1"),
        driver="gradient",
        method="b97d-d3(bj)",
        params_tweaks={
            "method": "b97d",
            "atm": atm,
        },
        qcel_object=True,
    )
    if atm:
        gradient = np.array(
            [
                [-2.2443259092095252e-4, -5.9115746657000033e-4, -2.3329260776706518e-4],
                [+1.5024984250557624e-4, -2.8252993039799747e-4, -1.1514673689182528e-4],
                [+6.4022428072682536e-4, -1.9311563977342141e-4, -8.5242235769601618e-5],
                [-3.8718976343890316e-4, +3.3773870712256794e-4, +1.4637097756315060e-4],
                [-1.5079749995767366e-4, +2.7785809307656525e-4, +1.3482408063774428e-4],
                [-5.4334015754870673e-5, +4.8740178552571852e-4, +2.0162061388786579e-4],
                [+4.2524309694851151e-4, -1.0353034250430949e-3, +5.4187321009509747e-4],
                [+2.5238880483243874e-4, -4.7570815576830443e-4, -1.9406096040939865e-4],
                [+5.0428650201581710e-4, +2.5654485601078660e-4, +9.9141715771808528e-5],
                [+5.7124533738733385e-4, -3.7922753617903088e-5, -2.2970835136004695e-5],
                [-5.2629396330476706e-4, -2.1821084773584023e-4, -7.7771997956860579e-5],
                [-3.1447177050715551e-4, -4.4596903523484569e-4, -1.7470534972738137e-4],
                [-4.7241980717575030e-4, +5.7789501566979232e-4, +6.3391340634065063e-4],
                [-4.4295668215519951e-4, +4.4288999267499747e-5, +2.9913326488103452e-4],
                [-5.5532764141529982e-5, +5.0128482842532645e-4, -2.3745756890447453e-4],
                [+1.6636630001049617e-4, +2.9846431594745745e-4, -1.7829054133393723e-4],
                [-4.6582228517614413e-5, +2.5574634046668398e-4, -4.7454953681660685e-4],
                [-3.4993078552581617e-5, +2.4269431262900884e-4, -2.6338889846419609e-4],
            ]
        )
    else:
        gradient = np.array(
            [
                [-2.2562967976217643e-04, -5.8711103078133340e-04, -2.3061961042068857e-04],
                [+1.5075490579531149e-04, -2.8381505981386033e-04, -1.1485932659166404e-04],
                [+6.3724554206790829e-04, -1.8986067913445720e-04, -8.2853482732288822e-05],
                [-3.8832449621583123e-04, +3.3282825489958215e-04, +1.4481440153935071e-04],
                [-1.5045846438610337e-04, +2.7794484580385330e-04, +1.3314407954925452e-04],
                [-4.9560016743308596e-05, +4.8509857651717612e-04, +2.0112762761380736e-04],
                [+4.2275101080035900e-04, -1.0281087399068071e-03, +5.3322479935030958e-04],
                [+2.5240228222254599e-04, -4.7591532248547261e-04, -1.9359997368309395e-04],
                [+5.0291675303562318e-04, +2.5554376917832927e-04, +9.9031452660504099e-05],
                [+5.7163177819497945e-04, -4.0110024194130451e-05, -2.3133113964066147e-05],
                [-5.2462065805044188e-04, -2.1790929154673013e-04, -7.7260813677796259e-05],
                [-3.1304532783227779e-04, -4.4767282257326171e-04, -1.7470716818989833e-04],
                [-4.7193012650844176e-04, +5.7610764389933421e-04, +6.3567846119653526e-04],
                [-4.4301897860924111e-04, +4.4549998176073525e-05, +2.9945109087021988e-04],
                [-5.8493469864357109e-05, +5.0096942382815647e-04, -2.4093109056512457e-04],
                [+1.6568457442979896e-04, +2.9913845360379754e-04, -1.7797313385544664e-04],
                [-4.4555263684831344e-05, +2.5438217773443301e-04, -4.6925288218547617e-04],
                [-3.3750364889515759e-05, +2.4393982679531790e-04, -2.6128131691443809e-04],
            ]
        )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == gradient



@pytest.mark.skipif(
    qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental models"
)
def test_gradient_tpss_d3zero(qcsk_version):
    thr = 1.0e-9

    molecule = {
        "symbols": "O C C F O F H".split(),
        "geometry": [
            [+4.877023733, -3.909030492, +1.796260143],
            [+6.112318716, -2.778558610, +0.091330457],
            [+7.360520527, -4.445334728, -1.932830640],
            [+7.978801077, -6.767751279, -1.031771494],
            [+6.374499300, -0.460299457, -0.213142194],
            [+5.637581753, -4.819746139, -3.831249370],
            [+9.040657008, -3.585225944, -2.750722946],
        ],
        "molecular_charge": -1,
    }

    atomic_input = get_atomic_input(
        version=qcsk_version,
        molecule=molecule,
        driver="gradient",
        method="",
        params_tweaks={
            "sr6": 1.166,
            "s8": 1.105,
            "alpha6": 14.0,
        },
        level_hint="d3zero",
        qcel_object=True,
    )

    gradient = np.array(
        [
            [+8.5996134689276694e-5, +1.3305341130186383e-4, -4.9354141710140030e-5],
            [+1.6831863787900489e-4, -3.6665329921172781e-4, -3.6501894583420186e-4],
            [-1.5843206014794366e-4, +1.9764667697917904e-4, +2.4585929437573800e-4],
            [-1.2283819613327093e-4, +4.7671549768466833e-4, -2.2634211331736250e-4],
            [-6.5711504962469335e-5, -7.6481921946330167e-5, +1.3704153057115390e-4],
            [+4.1919819225194306e-4, -7.9186754243182217e-6, +2.4311111443822667e-4],
            [-3.2653120357654074e-4, -3.5636168938333495e-4, +1.4703261476585888e-5],
        ]
    )

    atomic_result = run_qcschema(atomic_input)

    print(atomic_result.return_result)
    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == gradient
    assert "energy" in atomic_result.extras["dftd3"]
    assert "gradient" in atomic_result.extras["dftd3"]
    assert "virial" in atomic_result.extras["dftd3"]


@pytest.mark.skipif(
    qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental models"
)
def test_error_noargs(qcsk_version):
    molecule = get_molecule("molecule2")

    atomic_input = get_atomic_input(
        version=qcsk_version,
        molecule=molecule,
        driver="energy",
        method="",
    )

    atomic_result = run_qcschema(atomic_input)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input error"


@pytest.mark.skipif(
    qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental models"
)
def test_error_nomethod(qcsk_version):
    molecule = get_molecule("molecule2")

    atomic_input = get_atomic_input(
        version=qcsk_version,
        molecule=molecule,
        driver="energy",
        method="this-method-does-not-exist",
        level_hint="d3bj",
    )

    if qcsk_version == 1:
        from qcelemental.models import ComputeError
    elif qcsk_version == 2:
        from qcelemental.models.v2 import ComputeError
    else:
        raise RuntimeError(f"QCSchema v{qcsk_version} NYI")

    error = ComputeError(
        error_type="input error",
        error_message="No entry for 'this-method-does-not-exist' present",
    )

    atomic_result = run_qcschema(atomic_input)

    assert not atomic_result.success
    assert atomic_result.error == error


@pytest.mark.skipif(
    qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental models"
)
def test_error_level(qcsk_version):
    molecule = get_molecule("molecule2")

    atomic_input = get_atomic_input(
        version=qcsk_version,
        molecule=molecule,
        driver="energy",
        method="SCAN",
        level_hint="D42",
    )

    if qcsk_version == 1:
        from qcelemental.models import ComputeError
    elif qcsk_version == 2:
        from qcelemental.models.v2 import ComputeError
    else:
        raise RuntimeError(f"QCSchema v{qcsk_version} NYI")

    error = ComputeError(
        error_type="input error",
        error_message="Level 'D42' is invalid for this dispersion correction",
    )

    atomic_result = run_qcschema(atomic_input)

    assert not atomic_result.success
    assert atomic_result.error == error


@pytest.mark.skipif(
    qcel_v1 is None and qcel_v2 is None, reason="requires qcelemental models"
)
def test_ghost_pbe_d3bj(atm: bool, qcsk_version: int) -> None:
    thr = 1e-9

    molecule = get_molecule("counterpoise")

    atomic_input = get_atomic_input(
        version=qcsk_version,
        molecule=molecule,
        driver="gradient",
        method="pbe",
        params_tweaks={
            "method": "pbe",
            "atm": atm,
        },
    )

    if atm:
        gradient = np.array(
            [
                [+0.00000000e-0, +1.15093229e-7, +0.00000000e-0],
                [+5.37509663e-5, -1.90067439e-5, +0.00000000e-0],
                [-2.68754977e-5, -1.90067527e-5, -4.65496968e-5],
                [-2.68754977e-5, -1.90067527e-5, +4.65496968e-5],
                [+0.00000000e-0, +5.69051561e-5, +0.00000000e-0],
                [+0.00000000e-0, +0.00000000e-0, +0.00000000e-0],
                [+0.00000000e-0, +0.00000000e-0, +0.00000000e-0],
                [+0.00000000e-0, +0.00000000e-0, +0.00000000e-0],
                [+0.00000000e-0, +0.00000000e-0, +0.00000000e-0],
            ]
        )
    else:
        gradient = np.array(
            [
                [+3.6009409394390283e-11, +1.1526637296579583e-07, +0.0000000000000000e00],
                [+5.3425212496378062e-05, -1.8891583140214365e-05, +0.0000000000000000e00],
                [-2.6712623930671679e-05, -1.8891594924889618e-05, -4.6267592357667043e-05],
                [-2.6712623930671679e-05, -1.8891594924889618e-05, +4.6267592357667043e-05],
                [-6.4444409450829659e-13, +5.6559506617027801e-05, +0.0000000000000000e00],
                [+0.0000000000000000e00, +0.0000000000000000e00, +0.0000000000000000e00],
                [+0.0000000000000000e00, +0.0000000000000000e00, +0.0000000000000000e00],
                [+0.0000000000000000e00, +0.0000000000000000e00, +0.0000000000000000e00],
                [+0.0000000000000000e00, +0.0000000000000000e00, +0.0000000000000000e00],
            ]
        )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == gradient