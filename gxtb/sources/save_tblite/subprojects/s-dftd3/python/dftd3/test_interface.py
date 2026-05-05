# This file is part of s-dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

from dftd3.interface import (
    Structure,
    RationalDampingParam,
    ZeroDampingParam,
    ModifiedZeroDampingParam,
    ModifiedRationalDampingParam,
    OptimizedPowerDampingParam,
    DispersionModel,
    GeometricCounterpoise,
)
import pytest
from pytest import approx, raises
import numpy as np


@pytest.fixture
def numbers() -> np.ndarray:
    return np.array(
        [6, 7, 6, 7, 6, 6, 6, 8, 7, 6, 8, 7, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    )


@pytest.fixture
def positions() -> np.ndarray:
    return np.array(
        [
            [+2.02799738646442, +0.09231312124713, -0.14310895950963],
            [+4.75011007621000, +0.02373496014051, -0.14324124033844],
            [+6.33434307654413, +2.07098865582721, -0.14235306905930],
            [+8.72860718071825, +1.38002919517619, -0.14265542523943],
            [+8.65318821103610, -1.19324866489847, -0.14231527453678],
            [+6.23857175648671, -2.08353643730276, -0.14218299370797],
            [+5.63266886875962, -4.69950321056008, -0.13940509630299],
            [+3.44931709749015, -5.48092386085491, -0.14318454855466],
            [+7.77508917214346, -6.24427872938674, -0.13107140408805],
            [10.30229550927022, -5.39739796609292, -0.13672168520430],
            [12.07410272485492, -6.91573621641911, -0.13666499342053],
            [10.70038521493902, -2.79078533715849, -0.14148379504141],
            [13.24597858727017, -1.76969072232377, -0.14218299370797],
            [+7.40891694074004, -8.95905928176407, -0.11636933482904],
            [+1.38702118184179, +2.05575746325296, -0.14178615122154],
            [+1.34622199478497, -0.86356704498496, +1.55590600570783],
            [+1.34624089204623, -0.86133716815647, -1.84340893849267],
            [+5.65596919189118, +4.00172183859480, -0.14131371969009],
            [14.67430918222276, -3.26230980007732, -0.14344911021228],
            [13.50897177220290, -0.60815166181684, +1.54898960808727],
            [13.50780014200488, -0.60614855212345, -1.83214617078268],
            [+5.41408424778406, -9.49239668625902, -0.11022772492007],
            [+8.31919801555568, -9.74947502841788, +1.56539243085954],
            [+8.31511620712388, -9.76854236502758, -1.79108242206824],
        ]
    )


@pytest.fixture(params=[True, False])
def atm(request) -> bool:
    return request.param


@pytest.fixture
def model() -> DispersionModel:
    numbers = np.array([1, 1, 6, 5, 1, 15, 8, 17, 13, 15, 5, 1, 9, 15, 1, 15])
    positions = np.array(
        [  # Coordinates in Bohr
            [+2.79274810283778, +3.82998228828316, -2.79287054959216],
            [-1.43447454186833, +0.43418729987882, +5.53854345129809],
            [-3.26268343665218, -2.50644032426151, -1.56631149351046],
            [+2.14548759959147, -0.88798018953965, -2.24592534506187],
            [-4.30233097423181, -3.93631518670031, -0.48930754109119],
            [+0.06107643564880, -3.82467931731366, -2.22333344469482],
            [+0.41168550401858, +0.58105573172764, +5.56854609916143],
            [+4.41363836635653, +3.92515871809283, +2.57961724984000],
            [+1.33707758998700, +1.40194471661647, +1.97530004949523],
            [+3.08342709834868, +1.72520024666801, -4.42666116106828],
            [-3.02346932078505, +0.04438199934191, -0.27636197425010],
            [+1.11508390868455, -0.97617412809198, +6.25462847718180],
            [+0.61938955433011, +2.17903547389232, -6.21279842416963],
            [-2.67491681346835, +3.00175899761859, +1.05038813614845],
            [-4.13181080289514, -2.34226739863660, -3.44356159392859],
            [+2.85007173009739, -2.64884892757600, +0.71010806424206],
        ]
    )
    return DispersionModel(numbers, positions)


def test_rational_damping_noargs() -> None:
    """Check constructor of damping parameters for insufficient arguments"""

    with raises(TypeError):
        RationalDampingParam()

    with raises(TypeError, match="s8"):
        RationalDampingParam(a1=0.4, a2=5.0)

    with raises(TypeError, match="a1"):
        RationalDampingParam(s8=1.0, a2=5.0)

    with raises(TypeError, match="a2"):
        RationalDampingParam(s8=1.0, a1=0.4)


def test_zero_damping_noargs() -> None:
    """Check constructor of damping parameters for insufficient arguments"""

    with raises(TypeError):
        ZeroDampingParam()

    with raises(TypeError, match="s8"):
        ZeroDampingParam(rs6=1.2)

    with raises(TypeError, match="rs6"):
        ZeroDampingParam(s8=1.0)


def test_modified_rational_damping_noargs() -> None:
    """Check constructor of damping parameters for insufficient arguments"""

    with raises(TypeError):
        ModifiedRationalDampingParam()

    with raises(TypeError, match="s8"):
        ModifiedRationalDampingParam(a1=0.4, a2=5.0)

    with raises(TypeError, match="a1"):
        ModifiedRationalDampingParam(s8=1.0, a2=5.0)

    with raises(TypeError, match="a2"):
        ModifiedRationalDampingParam(s8=1.0, a1=0.4)


def test_modified_zero_damping_noargs() -> None:
    """Check constructor of damping parameters for insufficient arguments"""

    with raises(TypeError):
        ModifiedZeroDampingParam()

    with raises(TypeError, match="s8"):
        ModifiedZeroDampingParam(rs6=1.2, bet=1.0)

    with raises(TypeError, match="rs6"):
        ModifiedZeroDampingParam(s8=1.0, bet=1.0)

    with raises(TypeError, match="bet"):
        ModifiedZeroDampingParam(s8=1.0, rs6=1.2)

    with raises(TypeError):
        ModifiedZeroDampingParam(s8=1.0, rs6=1.2, bet=1.0, method="abc")


def test_optimized_power_damping_noargs() -> None:
    """Check constructor of damping parameters for insufficient arguments"""

    with raises(TypeError):
        OptimizedPowerDampingParam()

    with raises(TypeError, match="s8"):
        OptimizedPowerDampingParam(a1=0.3, a2=4.2, bet=1.0)

    with raises(TypeError, match="a1"):
        OptimizedPowerDampingParam(s8=1.0, a2=4.2, bet=1.0)

    with raises(TypeError, match="a2"):
        OptimizedPowerDampingParam(s8=1.0, a1=0.3, bet=1.0)

    with raises(TypeError, match="bet"):
        OptimizedPowerDampingParam(s8=1.0, a1=0.3, a2=4.2)

    with raises(TypeError):
        OptimizedPowerDampingParam(s8=1.0, a1=0.3, a2=4.2, bet=1.0, method="abc")


@pytest.mark.parametrize("cls", [Structure, DispersionModel, GeometricCounterpoise])
def test_structure(cls, numbers: np.ndarray, positions: np.ndarray) -> None:
    """check if the molecular structure data is working as expected."""

    rng = np.random.default_rng()

    # Constructor should raise an error for nuclear fusion input
    with raises(RuntimeError, match="Too close interatomic distances found"):
        cls(numbers, np.zeros((24, 3)))

    # The Python class should protect from garbage input like this
    with raises(ValueError, match="Dimension missmatch"):
        cls(np.array([1, 1, 1]), positions)

    # Also check for sane coordinate input
    with raises(ValueError, match="Expected tripels"):
        cls(numbers, rng.random(7))

    # Construct real molecule
    inst = cls(numbers, positions)

    # Try to update a structure with missmatched coordinates
    with raises(ValueError, match="Dimension missmatch for positions"):
        inst.update(rng.random(7))

    # Try to add a missmatched lattice
    with raises(ValueError, match="Invalid lattice provided"):
        inst.update(positions, rng.random(7))

    # Try to update a structure with nuclear fusion coordinates
    with raises(RuntimeError, match="Too close interatomic distances found"):
        inst.update(np.zeros((24, 3)))


def test_pbe0_d3_bj(atm: bool, model: DispersionModel) -> None:
    ref = -0.029489232932494884 if atm else -0.029589132634178342
    res = model.get_dispersion(RationalDampingParam(method="pbe0", atm=atm), grad=False)
    assert approx(res.get("energy")) == ref


def test_b3lyp_d3_zero(atm: bool, model: DispersionModel) -> None:
    ref = -0.022714272555175656 if atm else -0.022814172019166058
    res = model.get_dispersion(ZeroDampingParam(method="b3lyp", atm=atm), grad=False)
    assert approx(res.get("energy")) == ref


def test_pbe_d3_bjm(atm: bool, model: DispersionModel) -> None:
    ref = -0.06327406660942464 if atm else -0.06337396631110809
    res = model.get_dispersion(
        ModifiedRationalDampingParam(method="pbe", atm=atm), grad=False
    )
    assert approx(res.get("energy")) == ref


def test_bp_d3_zerom(atm: bool, model: DispersionModel) -> None:
    ref = -0.026013316869036292 if atm else -0.026113216333026695
    res = model.get_dispersion(
        ModifiedZeroDampingParam(method="bp", atm=atm), grad=False
    )
    assert approx(res.get("energy")) == ref


def test_b97d_d3_op(atm: bool, model: DispersionModel) -> None:
    ref = -0.07681029606751344 if atm else -0.07691018779028679
    res = model.get_dispersion(
        OptimizedPowerDampingParam(method="b97d", atm=atm), grad=False
    )
    assert approx(res.get("energy")) == ref


def test_gcp_empty(numbers: np.ndarray, positions: np.ndarray) -> None:
    gcp = GeometricCounterpoise(
        numbers,
        positions,
    )

    res = gcp.get_counterpoise(grad=False)
    assert approx(res.get("energy")) == 0.0


@pytest.mark.parametrize("method", ["b973c"])
def test_gcp_3c(numbers: np.ndarray, positions: np.ndarray, method: str) -> None:
    gcp = GeometricCounterpoise(
        numbers,
        positions,
        method=method,
    )
    ref = {
        "b973c": -0.07653225860427701,
        "pbeh3c": 0.04977771585466725,
    }[method]

    res = gcp.get_counterpoise(grad=False)
    assert approx(res.get("energy")) == ref, res.get("energy")


def test_pair_resolved() -> None:
    """Calculate pairwise resolved dispersion energy for a molecule"""
    thr = 1.0e-8

    numbers = np.array([16, 16, 16, 16, 16, 16, 16, 16])
    positions = np.array(
        [
            [-4.15128787379191, +1.71951973863958, -0.93066267097296],
            [-4.15128787379191, -1.71951973863958, +0.93066267097296],
            [-1.71951973863958, -4.15128787379191, -0.93066267097296],
            [+1.71951973863958, -4.15128787379191, +0.93066267097296],
            [+4.15128787379191, -1.71951973863958, -0.93066267097296],
            [+4.15128787379191, +1.71951973863958, +0.93066267097296],
            [+1.71951973863958, +4.15128787379191, -0.93066267097296],
            [-1.71951973863958, +4.15128787379191, +0.93066267097296],
        ]
    )
    pair_disp2 = np.array(
        [
            [
                -0.00000000,
                -0.00153111,
                -0.00108052,
                -0.00032865,
                -0.00023796,
                -0.00032865,
                -0.00108052,
                -0.00153111,
            ],
            [
                -0.00153111,
                -0.00000000,
                -0.00153111,
                -0.00108052,
                -0.00032865,
                -0.00023796,
                -0.00032865,
                -0.00108052,
            ],
            [
                -0.00108052,
                -0.00153111,
                -0.00000000,
                -0.00153111,
                -0.00108052,
                -0.00032865,
                -0.00023796,
                -0.00032865,
            ],
            [
                -0.00032865,
                -0.00108052,
                -0.00153111,
                -0.00000000,
                -0.00153111,
                -0.00108052,
                -0.00032865,
                -0.00023796,
            ],
            [
                -0.00023796,
                -0.00032865,
                -0.00108052,
                -0.00153111,
                -0.00000000,
                -0.00153111,
                -0.00108052,
                -0.00032865,
            ],
            [
                -0.00032865,
                -0.00023796,
                -0.00032865,
                -0.00108052,
                -0.00153111,
                -0.00000000,
                -0.00153111,
                -0.00108052,
            ],
            [
                -0.00108052,
                -0.00032865,
                -0.00023796,
                -0.00032865,
                -0.00108052,
                -0.00153111,
                -0.00000000,
                -0.00153111,
            ],
            [
                -0.00153111,
                -0.00108052,
                -0.00032865,
                -0.00023796,
                -0.00032865,
                -0.00108052,
                -0.00153111,
                -0.00000000,
            ],
        ]
    )
    pair_disp3 = np.array(
        [
            [
                0.00000000e-00,
                1.08616452e-07,
                2.91526483e-07,
                3.95872130e-07,
                3.18133443e-07,
                3.95872130e-07,
                2.91526483e-07,
                1.08616452e-07,
            ],
            [
                1.08616452e-07,
                0.00000000e-00,
                1.08616452e-07,
                2.91526483e-07,
                3.95872130e-07,
                3.18133443e-07,
                3.95872130e-07,
                2.91526483e-07,
            ],
            [
                2.91526483e-07,
                1.08616452e-07,
                0.00000000e-00,
                1.08616452e-07,
                2.91526483e-07,
                3.95872130e-07,
                3.18133443e-07,
                3.95872130e-07,
            ],
            [
                3.95872130e-07,
                2.91526483e-07,
                1.08616452e-07,
                0.00000000e-00,
                1.08616452e-07,
                2.91526483e-07,
                3.95872130e-07,
                3.18133443e-07,
            ],
            [
                3.18133443e-07,
                3.95872130e-07,
                2.91526483e-07,
                1.08616452e-07,
                0.00000000e-00,
                1.08616452e-07,
                2.91526483e-07,
                3.95872130e-07,
            ],
            [
                3.95872130e-07,
                3.18133443e-07,
                3.95872130e-07,
                2.91526483e-07,
                1.08616452e-07,
                0.00000000e-00,
                1.08616452e-07,
                2.91526483e-07,
            ],
            [
                2.91526483e-07,
                3.95872130e-07,
                3.18133443e-07,
                3.95872130e-07,
                2.91526483e-07,
                1.08616452e-07,
                0.00000000e-00,
                1.08616452e-07,
            ],
            [
                1.08616452e-07,
                2.91526483e-07,
                3.95872130e-07,
                3.18133443e-07,
                3.95872130e-07,
                2.91526483e-07,
                1.08616452e-07,
                0.00000000e-00,
            ],
        ]
    )

    model = DispersionModel(numbers, positions)

    res = model.get_pairwise_dispersion(
        RationalDampingParam(a1=0.5545, s8=2.2609, a2=3.2297)
    )

    assert approx(res.get("additive pairwise energy"), abs=thr) == pair_disp2
    assert approx(res.get("non-additive pairwise energy"), abs=thr) == pair_disp3
