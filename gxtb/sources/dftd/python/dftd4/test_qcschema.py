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

# pyright: reportPossiblyUnboundVariable=false
# pylint: disable=possibly-used-before-assignment

from typing import Any, Dict

import numpy as np
from pytest import approx, mark

try:
    import qcelemental as qcel
    from dftd4.qcschema import run_qcschema

    has_qcschema = True
except ModuleNotFoundError:
    has_qcschema = False

pytestmark = mark.skipif(not has_qcschema, reason="requires qcelemental")


def get_example_molecule() -> Dict[str, Any]:
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


def test_energy_r2scan_d4() -> None:
    thr = 1e-9

    atomic_input = qcel.models.AtomicInput(
        molecule=get_example_molecule(),
        driver="energy",
        model={
            "method": "r2scan",
        },
    )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == -0.005001101011286166


def test_energy_explicit_damping_hint() -> None:
    thr = 1e-9

    atomic_input = qcel.models.AtomicInput(
        molecule=get_example_molecule(),
        driver="energy",
        model={
            "method": "r2scan",
        },
        keywords={
            "damping_hint": {
                "2b": "rational",
                "3b": "zero-avg",
            }
        },
    )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == -0.005001101011286166


def test_energy_r2scan_d4s() -> None:
    thr = 1e-9

    atomic_input = qcel.models.AtomicInput(
        molecule=get_example_molecule(),
        driver="energy",
        model={
            "method": "r2scan",
        },
        keywords={
            "level_hint": "d4s",
        },
    )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == -0.00509785822000568


def test_energy_r2scan_3c() -> None:
    thr = 1e-9

    atomic_input = qcel.models.AtomicInput(
        molecule=get_example_molecule(),
        driver="energy",
        model={"method": "d4"},
        keywords={
            "params_tweaks": {
                "s6": 1.00,
                "s8": 0.00,
                "s9": 2.00,
                "a1": 0.42,
                "a2": 5.65,
                "ga": 2.00,
                "gc": 1.00,
            },
        },
    )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == -6.0533536923248e-03


def test_energy_lh20t_d4() -> None:
    thr = 1e-9

    atomic_input = qcel.models.AtomicInput(
        molecule=get_example_molecule(),
        driver="energy",
        model={"method": ""},
        keywords={
            "params_tweaks": {
                "s8": 0.113,
                "a1": 0.479,
                "a2": 4.635,
            },
        },
    )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == -0.010064263146257654


def test_energy_lh20t_d4_noatm() -> None:
    thr = 1e-9

    atomic_input = qcel.models.AtomicInput(
        molecule=get_example_molecule(),
        driver="energy",
        model={"method": ""},
        keywords={
            "params_tweaks": {
                "s8": 0.113,
                "a1": 0.479,
                "a2": 4.635,
            },
            "damping_hint": {
                "3b": "none"
            }, 
        },
    )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == -0.010048292243378118


def test_energy_lh20t_d4s() -> None:
    thr = 1e-9

    atomic_input = qcel.models.AtomicInput(
        molecule=get_example_molecule(),
        driver="energy",
        model={"method": ""},
        keywords={
            "params_tweaks": {
                "s8": 0.113,
                "a1": 0.479,
                "a2": 4.635,
            },
            "level_hint": "d4s",
        },
    )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == -0.010252088837042048


def test_energy_m06l_d4() -> None:
    thr = 1e-6

    atomic_input = qcel.models.AtomicInput(
        molecule={
            "symbols": "Li Cl F H H Na B H C H H F C H H H".split(" "),
            "geometry": [
                [+2.06521084486823, +0.08218432748393, -3.31794862285397],
                [-2.32042477766402, -5.46684392277772, -4.04262137940086],
                [+2.26054697581965, -3.03018226193694, -4.64150015772052],
                [-0.14839777550969, -0.43671669092912, +5.46349590128611],
                [+1.25506764958846, +3.72255296450239, +3.12461655723367],
                [-3.70871338035337, -1.70938913801338, -1.44451499871032],
                [-0.50789889989427, -0.08663524018430, +3.25334078665520],
                [+1.95889668026174, -0.69562876271128, +1.49780663184774],
                [-0.13609934796341, +2.63877103476555, +2.01392332556491],
                [-2.03651637484513, -1.53932035918944, +2.35235290216748],
                [+1.01732693290435, +4.50598465234050, -1.41360090365614],
                [-1.21371363315079, +1.45147299379122, -2.46002878669847],
                [+0.99959759042027, +2.58862220349908, -0.61148816085282],
                [+1.82755250562899, -1.68931361321449, +2.62013289657120],
                [-1.90490242357524, +3.74722938063725, +2.11670736907658],
                [+0.59246743346417, -4.08278756806328, -4.51067336050986],
            ],
        },
        driver="energy",
        model={
            "method": "m06l",
        },
        keywords={
            "property": True,
        },
    )

    charges = [
        +6.4130068073040725e-1,
        -2.8383857750967989e-1,
        -4.4179022311738259e-1,
        +7.1454451873716449e-2,
        +9.1415362551287238e-2,
        +6.9802661015384404e-1,
        -3.9007049547150685e-1,
        +2.9957268041964854e-2,
        -1.7345919709167879e-1,
        +2.4649274982777081e-2,
        +8.2070988323221294e-2,
        -3.3909452402278339e-1,
        -2.6535793149085934e-1,
        +3.6883214212240473e-2,
        +8.3750568494545596e-2,
        +1.3410252933988615e-1,
    ]

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert (
        approx(atomic_result.return_result, abs=thr) == -0.0013314656225517764
    )
    assert (
        approx(atomic_result.extras["dftd4"]["partial charges"], abs=thr)
        == charges
    )


def test_gradient_b97m_d4() -> None:
    thr = 1e-9

    atomic_input = qcel.models.AtomicInput(
        molecule=get_example_molecule(),
        driver="gradient",
        model={
            "method": "b97m-D4",
        },
        keywords={},
    )
    gradient = np.array(
        [
            [-2.44611355e-04, -5.62668284e-04, -2.23387025e-04],
            [+2.55729172e-04, -4.81570822e-04, -1.97100251e-04],
            [+6.28561864e-04, -1.60531557e-04, -7.38793802e-05],
            [+5.07019846e-04, +1.84224291e-04, +6.97725810e-05],
            [-1.48274928e-04, +2.74024495e-04, +1.28583321e-04],
            [-4.66509227e-04, -2.63846497e-04, -9.72992491e-05],
            [-1.45316573e-04, +1.14524656e-04, +5.68088686e-04],
            [+8.73699092e-05, -1.64562363e-04, -6.71815237e-05],
            [+1.99967459e-04, +1.12215186e-04, +4.54083544e-05],
            [+2.04308124e-04, -1.82947927e-05, -1.01832111e-05],
            [-2.17520641e-04, -8.08940986e-05, -2.67701158e-05],
            [-1.08257557e-04, -1.62398104e-04, -6.37544043e-05],
            [-3.95087002e-04, +5.42941863e-04, +2.98100565e-04],
            [-8.67554935e-05, +4.46950178e-05, +9.10061596e-05],
            [-4.95209003e-05, +3.52211029e-04, -2.37113109e-04],
            [+2.26441621e-05, +9.26206419e-05, -3.41527716e-05],
            [-3.29465221e-05, +1.14378254e-04, -1.11215443e-04],
            [-1.08003368e-05, +6.29310823e-05, -5.89231846e-05],
        ]
    )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == gradient


def test_gradient_tpss_d4s() -> None:
    thr = 1e-9

    atomic_input = qcel.models.AtomicInput(
        molecule=get_example_molecule(),
        driver="gradient",
        model={
            "method": "TPSS",
        },
        keywords={
            "level_hint": "d4s",
        },
    )
    gradient = np.array(
        [
            [-1.36067766e-04, -4.72731818e-04, -1.88160742e-04],
            [+2.00334329e-04, -3.77344166e-04, -1.54443106e-04],
            [+4.85748705e-04, -1.86465318e-04, -8.16713264e-05],
            [+3.40248245e-04, +1.04514066e-04, +4.31391229e-05],
            [-7.10734227e-05, +1.27160703e-04, +7.31294101e-05],
            [-2.97997359e-04, -1.89214279e-04, -6.66535024e-05],
            [-1.16456272e-04, +4.84627746e-05, +5.51782210e-04],
            [+6.64898493e-05, -1.25194763e-04, -5.12043373e-05],
            [+1.50558513e-04, +7.37719529e-05, +3.08377180e-05],
            [+1.62329419e-04, -2.32012807e-05, -1.15835288e-05],
            [-1.54748174e-04, -6.78861151e-05, -2.21532014e-05],
            [-7.86055402e-05, -1.34246120e-04, -5.28699703e-05],
            [-3.39623373e-04, +5.01841681e-04, +1.79953770e-04],
            [-7.38516701e-05, +7.66383868e-05, +6.65527953e-05],
            [-7.64722239e-05, +3.68207143e-04, -2.27167892e-04],
            [-2.65824875e-05, +9.70334513e-05, -2.37042257e-05],
            [-3.95478770e-05, +1.43973339e-04, -4.59543086e-05],
            [+5.31710551e-06, +3.46803626e-05, -1.98288866e-05],
        ]
    )

    atomic_result = run_qcschema(atomic_input)
    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == gradient


def test_gradient_tpss_d4() -> None:
    thr = 1.0e-9

    atomic_input = qcel.models.AtomicInput(
        molecule={
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
        },
        driver="gradient",
        model={
            "method": "",
        },
        keywords={
            "params_tweaks": {
                "s8": 1.76596355,
                "a1": 0.42822303,
                "a2": 4.54257102,
            },
            "pair_resolved": True,
        },
    )
    gradient = np.array(
        [
            [-1.47959449e-04, +2.95411758e-05, +2.69548700e-04],
            [-9.35127114e-06, +2.37170201e-05, +2.01896552e-05],
            [+4.90329914e-05, -9.08636769e-06, -4.75377709e-05],
            [+1.28728759e-04, -1.98165729e-04, -5.63375526e-05],
            [-1.02103452e-05, +3.04524759e-04, +9.14084679e-05],
            [-4.27301731e-06, -1.00071311e-04, -2.42775990e-04],
            [-5.96766712e-06, -5.04595462e-05, -3.44955102e-05],
        ]
    )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == gradient
    assert "energy" in atomic_result.extras["dftd4"]
    assert "gradient" in atomic_result.extras["dftd4"]
    assert "virial" in atomic_result.extras["dftd4"]
    assert "additive pairwise energy" in atomic_result.extras["dftd4"]
    assert "non-additive pairwise energy" in atomic_result.extras["dftd4"]
    assert (
        approx(atomic_result.extras["dftd4"]["energy"])
        == atomic_result.extras["dftd4"]["additive pairwise energy"].sum()
        + atomic_result.extras["dftd4"]["non-additive pairwise energy"].sum()
    )


def test_error_noargs() -> None:
    atomic_input = qcel.models.AtomicInput(
        molecule={
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
        },
        driver="energy",
        model={"method": ""},
        keywords={},
    )
    error = qcel.models.ComputeError(
        error_type="input error",
        error_message="Method name or complete damping parameter set required",
    )

    atomic_result = run_qcschema(atomic_input)

    assert not atomic_result.success
    assert atomic_result.error == error


def test_error_nomethod() -> None:
    atomic_input = qcel.models.AtomicInput(
        molecule={
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
        },
        driver="energy",
        model={
            "method": "this-method-does-not-exist",
        },
        keywords={
            "level_hint": "D4",
        },
    )
    error = qcel.models.ComputeError(
        error_type="input error",
        error_message="No D4 damping parameters available for this functional.",
    )

    atomic_result = run_qcschema(atomic_input)

    assert not atomic_result.success
    assert atomic_result.error == error


def test_error_level() -> None:
    atomic_input = qcel.models.AtomicInput(
        molecule={
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
        },
        driver="energy",
        model={
            "method": "SCAN",
        },
        keywords={
            "level_hint": "D42",
        },
    )
    error = qcel.models.ComputeError(
        error_type="input error",
        error_message="Level 'D42' is invalid for this dispersion correction",
    )

    atomic_result = run_qcschema(atomic_input)

    assert not atomic_result.success
    assert atomic_result.error == error


def test_ghost_pbe_d4() -> None:
    thr = 1e-9

    atomic_input = qcel.models.AtomicInput(
        molecule={
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
        },
        driver="gradient",
        model={
            "method": "pbe",
        },
    )
    gradient = np.array(
        [
            [+0.00000000e-0, +2.76351835e-7, +0.00000000e-0],
            [+6.38066609e-5, -2.25449242e-5, +0.00000000e-0],
            [-3.19033757e-5, -2.25449510e-5, -5.52582163e-5],
            [-3.19033757e-5, -2.25449510e-5, +5.52582163e-5],
            [+0.00000000e-0, +6.73584743e-5, +0.00000000e-0],
            [+0.00000000e-0, +0.00000000e-0, +0.00000000e-0],
            [+0.00000000e-0, +0.00000000e-0, +0.00000000e-0],
            [+0.00000000e-0, +0.00000000e-0, +0.00000000e-0],
            [+0.00000000e-0, +0.00000000e-0, +0.00000000e-0],
        ]
    )

    atomic_result = run_qcschema(atomic_input)

    assert atomic_result.success
    assert approx(atomic_result.return_result, abs=thr) == gradient
