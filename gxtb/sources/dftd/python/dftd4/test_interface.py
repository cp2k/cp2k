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

import numpy as np
from dftd4.interface import DampingFunction, DampingParam, DispersionModel, Structure
from pytest import approx, raises


def test_param_noargs() -> None:
    """Check constructor of damping parameters for insufficient arguments"""

    with raises(TypeError):
        DampingParam()


    with raises(ValueError, match="Model name must be provided"):
        DampingParam(method="blyp")

    with raises(ValueError, match="Explicit parameters cannot be mixed"):
        DampingParam(method="blyp", model="d4", s8=1.0, a1=0.4, a2=5.0)

    with raises(ValueError, match="Unknown dispersion model"):
        DampingParam(method="blyp", model="D42")

    with raises(ValueError, match="Invalid model or damping type"):
        DampingParam(method="blyp", model="d4", damping_2b="none", damping_3b="none")

    with raises(RuntimeError, match="No D4 damping parameters available"):
        DampingParam(method="abc", model="d4")

    DampingParam(method="blyp", model="d4")
    DampingParam(method="blyp", model="d4", damping_2b="rational", damping_3b="zero-avg")
    DampingParam(method="blyp", model="d4", damping_2b="rational", damping_3b="none")


    with raises(TypeError, match="s8"):
        DampingParam(model="d4", a1=0.4, a2=5.0)

    with raises(TypeError, match="a1"):
        DampingParam(model="d4", s8=1.0, a2=5.0)

    with raises(TypeError, match="a2"):
        DampingParam(model="d4", s8=1.0, a1=0.4)

    with raises(ValueError, match="Unknown dispersion model"):
        DampingParam(model="D42", s8=1.0, a1=0.4, a2=5.0)

    with raises(TypeError, match="a3|unexpected keyword"):
        DampingParam(model="d4s", s8=1.0, a1=0.4, a2=5.0, a3=0.0)

    DampingParam(model="d4", s8=1.0, a1=0.4, a2=5.0)
    DampingParam(model="d4s", s6=0.9, s8=1.0, s9=1.1, a1=0.4, a2=5.0,
                 rs9=1.0, alp=16.0)


    with raises(TypeError, match="a1"):
        DampingParam(s8=1.0, a2=5.0, a3=0.0, a4=0.0, rs6=0.0, rs8=0.0,
                     rs9=1.0, alp=16.0, bet=0.0)

    with raises(TypeError, match="a2"):
        DampingParam(s8=1.0, a1=0.4, a3=0.0, a4=0.0, rs6=0.0, rs8=0.0,
                     rs9=1.0, alp=16.0, bet=0.0)

    with raises(TypeError, match="a3"):
        DampingParam(s8=1.0, a1=0.4, a2=5.0, a4=0.0, rs6=0.0, rs8=0.0,
                     rs9=1.0, alp=16.0, bet=0.0)

    with raises(TypeError, match="a4"):
        DampingParam(s8=1.0, a1=0.4, a2=5.0, a3=0.0, rs6=0.0, rs8=0.0,
                     rs9=1.0, alp=16.0, bet=0.0)

    with raises(TypeError, match="rs6"):
        DampingParam(s8=1.0, a1=0.4, a2=5.0, a3=0.0, a4=0.0, rs8=0.0,
                     rs9=1.0, alp=16.0, bet=0.0)

    with raises(TypeError, match="rs8"):
        DampingParam(s8=1.0, a1=0.4, a2=5.0, a3=0.0, a4=0.0, rs6=0.0,
                     rs9=1.0, alp=16.0, bet=0.0)

    with raises(TypeError, match="rs9"):
        DampingParam(s8=1.0, a1=0.4, a2=5.0, a3=0.0, a4=0.0, rs6=0.0,
                     rs8=1.0, alp=16.0, bet=0.0)

    with raises(TypeError, match="alp"):
        DampingParam(s8=1.0, a1=0.4, a2=5.0, a3=0.0, a4=0.0, rs6=0.0,
                     rs8=1.0, rs9=1.0, bet=0.0)

    with raises(TypeError, match="bet"):
        DampingParam(s8=1.0, a1=0.4, a2=5.0, a3=0.0, a4=0.0, rs6=0.0,
                     rs8=0.0, rs9=1.0, alp=16.0)

    DampingParam(
        s6=1.0, s8=1.0, s9=1.0,
        a1=0.4, a2=5.0, a3=0.0, a4=0.0,
        rs6=0.0, rs8=0.0, rs9=1.0,
        alp=16.0, bet=0.0
    )


def test_damping_noargs() -> None:
    """Check constructor of damping function for insufficient/invalid arguments."""

    with raises(ValueError, match="Either two-body damping"):
        DampingFunction()

    with raises(ValueError, match="Cannot provide both"):
        DampingFunction(model="d4", damping_2b="rational")

    with raises(ValueError, match="Cannot provide both"):
        DampingFunction(model="d4", damping_2b="rational", damping_3b="zero-avg")

    with raises(ValueError, match="Either two-body damping"):
        DampingFunction(damping_3b="zero-avg")

    with raises(ValueError, match="Invalid damping type:"):
        DampingFunction(damping_2b="none")

    with raises(ValueError, match="Invalid damping type:"):
        DampingFunction(damping_2b="rational", damping_3b="abc")

    DampingFunction(damping_2b="rational") 
    DampingFunction(damping_2b="rational", damping_3b="zero-avg")
    DampingFunction(damping_2b="screened", damping_3b="none")

    with raises(ValueError, match="Unknown dispersion model"):
        DampingFunction(model="D42")

    DampingFunction(model="d4")
    DampingFunction(model=" D4S ")
    DampingFunction(model="d4", damping_3b="none")


def test_damping_check_params() -> None:
    """Check compatibility of damping function and parameters."""

    damp = DampingFunction(model="d4")
    param = DampingParam(method="blyp", model="d4")
    damp.check_params(param)

    param = DampingParam(
        s6=1.0, s8=1.0, s9=1.0,
        a1=0.4, a2=5.0, a3=0.0, a4=0.0,
        rs6=0.0, rs8=0.0, rs9=1.0,
        alp=16.0, bet=0.0
    )
    damp.check_params(param)


def test_structure() -> None:
    """check if the molecular structure data is working as expected."""

    numbers = np.array(
        [6, 7, 6, 7, 6, 6, 6, 8, 7, 6, 8, 7, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    )
    positions = np.array(
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

    # Constructor should raise an error for nuclear fusion input
    with raises(RuntimeError, match="Too close interatomic distances found"):
        Structure(numbers, np.zeros((24, 3)))

    # The Python class should protect from garbage input like this
    with raises(ValueError, match="Dimension missmatch"):
        Structure(np.array([1, 1, 1]), positions)

    # Also check for sane coordinate input
    with raises(ValueError, match="Expected triples"):
        Structure(numbers, np.random.rand(7))

    # Construct real molecule
    mol = Structure(numbers, positions)

    # Try to update a structure with missmatched coordinates
    with raises(ValueError, match="Dimension missmatch for positions"):
        mol.update(np.random.rand(7))

    # Try to add a missmatched lattice
    with raises(ValueError, match="Invalid lattice provided"):
        mol.update(positions, np.random.rand(7))

    # Try to update a structure with nuclear fusion coordinates
    with raises(RuntimeError, match="Too close interatomic distances found"):
        mol.update(np.zeros((24, 3)))


def test_blypd4() -> None:
    """Use BLYP-D4 for a mindless molecule"""
    thr = 1.0e-7

    numbers = np.array(
        [
            1,
            1,
            6,
            5,
            1,
            15,
            8,
            17,
            13,
            15,
            5,
            1,
            9,
            15,
            1,
            15,
        ]
    )
    positions = np.array(
        [
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

    model = DispersionModel(numbers, positions)

    res = model.get_dispersion(DampingFunction(model="d4"),
                               DampingParam(method="blyp", model="d4"),
                               grad=False)

    assert approx(res.get("energy"), abs=thr) == -0.06991716314879085

    res = model.get_dispersion(DampingFunction(model="d4"),
                               DampingParam(method="blyp", model="d4"),
                               grad=True)

    assert approx(res.get("energy"), abs=thr) == -0.06991716314879085


def test_tpssd4s() -> None:
    """Use TPSS-D4S for a mindless molecule"""
    thr = 1.0e-7

    numbers = np.array(
        [
            1,
            1,
            6,
            5,
            1,
            15,
            8,
            17,
            13,
            15,
            5,
            1,
            9,
            15,
            1,
            15,
        ]
    )
    positions = np.array(
        [
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

    model = DispersionModel(numbers, positions, model="d4s")

    res = model.get_dispersion(DampingFunction(model="d4s"),
                               DampingParam(method="tpss", model="d4s"),
                               grad=False)

    assert approx(res.get("energy"), abs=thr) == -0.046233140236052253

    res = model.get_dispersion(DampingFunction(model="d4s"),
                               DampingParam(method="tpss", model="d4s"),
                               grad=True)

    assert approx(res.get("energy"), abs=thr) == -0.046233140236052253


def test_pbed4() -> None:
    """Use PBE-D4 for a mindless molecule"""
    thr = 1.0e-7

    numbers = np.array(
        [
            1,
            9,
            15,
            13,
            1,
            1,
            13,
            5,
            3,
            15,
            8,
            1,
            1,
            5,
            16,
            1,
        ]
    )
    positions = np.array(
        [
            [-2.14132037405479, -1.34402701877044, -2.32492500904728],
            [+4.46671289205392, -2.04800110524830, +0.44422406067087],
            [-4.92212517643478, -1.73734240529793, +0.96890323821450],
            [-1.30966093045696, -0.52977363497805, +3.44453452239668],
            [-4.34208759006189, -4.30470270977329, +0.39887431726215],
            [+0.61788392767516, +2.62484136683297, -3.28228926932647],
            [+4.23562873444840, -1.68839322682951, -3.53824299552792],
            [+2.23130060612446, +1.93579813100155, -1.80384647554323],
            [-2.32285463652832, +2.90603947535842, -1.39684847191937],
            [+2.34557941578250, +2.86074312333371, +1.82827238641666],
            [-3.66431367659153, -0.42910188232667, -1.81957402856634],
            [-0.34927881505446, -1.75988134003940, +5.98017466326572],
            [+0.29500802281217, -2.00226104143537, +0.53023447931897],
            [+2.10449364205058, -0.56741404446633, +0.30975625014335],
            [-1.59355304432499, +3.69176153150419, +2.87878226787916],
            [+4.34858700256050, +2.39171478113440, -2.61802993563738],
        ]
    )

    model = DispersionModel(numbers, positions)

    res = model.get_dispersion(DampingFunction(model="d4"),
                               DampingParam(method="pbe", model="d4"),
                               grad=False)

    assert approx(res.get("energy"), abs=thr) == -0.028415184156428127

    res = model.get_dispersion(DampingFunction(model="d4"),
                               DampingParam(method="pbe", model="d4"),
                               grad=True)

    assert approx(res.get("energy"), abs=thr) == -0.028415184156428127


def test_r2scan3c() -> None:
    """Use r2SCAN-3c for a mindless molecule"""
    thr = 1.0e-8

    numbers = np.array(
        [
            1,
            9,
            15,
            13,
            1,
            1,
            13,
            5,
            3,
            15,
            8,
            1,
            1,
            5,
            16,
            1,
        ]
    )
    positions = np.array(
        [
            [-2.14132037405479, -1.34402701877044, -2.32492500904728],
            [+4.46671289205392, -2.04800110524830, +0.44422406067087],
            [-4.92212517643478, -1.73734240529793, +0.96890323821450],
            [-1.30966093045696, -0.52977363497805, +3.44453452239668],
            [-4.34208759006189, -4.30470270977329, +0.39887431726215],
            [+0.61788392767516, +2.62484136683297, -3.28228926932647],
            [+4.23562873444840, -1.68839322682951, -3.53824299552792],
            [+2.23130060612446, +1.93579813100155, -1.80384647554323],
            [-2.32285463652832, +2.90603947535842, -1.39684847191937],
            [+2.34557941578250, +2.86074312333371, +1.82827238641666],
            [-3.66431367659153, -0.42910188232667, -1.81957402856634],
            [-0.34927881505446, -1.75988134003940, +5.98017466326572],
            [+0.29500802281217, -2.00226104143537, +0.53023447931897],
            [+2.10449364205058, -0.56741404446633, +0.30975625014335],
            [-1.59355304432499, +3.69176153150419, +2.87878226787916],
            [+4.34858700256050, +2.39171478113440, -2.61802993563738],
        ]
    )

    model = DispersionModel(numbers, positions, ga=2.0, gc=1.0)

    res = model.get_dispersion(DampingFunction(model="d4"),
        DampingParam(model="d4", s8=0.0, a1=0.42, a2=5.65, s9=2.0), grad=False
    )

    assert approx(res.get("energy"), abs=thr) == -0.008016697276824889

    res = model.get_dispersion(DampingFunction(model="d4"),
        DampingParam(model="d4", s8=0.0, a1=0.42, a2=5.65, s9=2.0), grad=True
    )

    assert approx(res.get("energy"), abs=thr) == -0.008016697276824889


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
                -0.00000000e-0,
                -5.80599854e-4,
                -4.74689854e-4,
                -2.11149449e-4,
                -1.63163128e-4,
                -2.11149449e-4,
                -4.74689854e-4,
                -5.80599854e-4,
            ],
            [
                -5.80599854e-4,
                -0.00000000e-0,
                -5.80599854e-4,
                -4.74689854e-4,
                -2.11149449e-4,
                -1.63163128e-4,
                -2.11149449e-4,
                -4.74689854e-4,
            ],
            [
                -4.74689854e-4,
                -5.80599854e-4,
                -0.00000000e-0,
                -5.80599854e-4,
                -4.74689854e-4,
                -2.11149449e-4,
                -1.63163128e-4,
                -2.11149449e-4,
            ],
            [
                -2.11149449e-4,
                -4.74689854e-4,
                -5.80599854e-4,
                -0.00000000e-0,
                -5.80599854e-4,
                -4.74689854e-4,
                -2.11149449e-4,
                -1.63163128e-4,
            ],
            [
                -1.63163128e-4,
                -2.11149449e-4,
                -4.74689854e-4,
                -5.80599854e-4,
                -0.00000000e-0,
                -5.80599854e-4,
                -4.74689854e-4,
                -2.11149449e-4,
            ],
            [
                -2.11149449e-4,
                -1.63163128e-4,
                -2.11149449e-4,
                -4.74689854e-4,
                -5.80599854e-4,
                -0.00000000e-0,
                -5.80599854e-4,
                -4.74689854e-4,
            ],
            [
                -4.74689854e-4,
                -2.11149449e-4,
                -1.63163128e-4,
                -2.11149449e-4,
                -4.74689854e-4,
                -5.80599854e-4,
                -0.00000000e-0,
                -5.80599854e-4,
            ],
            [
                -5.80599854e-4,
                -4.74689854e-4,
                -2.11149449e-4,
                -1.63163128e-4,
                -2.11149449e-4,
                -4.74689854e-4,
                -5.80599854e-4,
                -0.00000000e-0,
            ],
        ]
    )
    pair_disp3 = np.array(
        [
            [
                0.00000000e-0,
                3.39353850e-7,
                8.74462839e-7,
                1.17634100e-6,
                9.86937310e-7,
                1.17634100e-6,
                8.74462839e-7,
                3.39353850e-7,
            ],
            [
                3.39353850e-7,
                0.00000000e-0,
                3.39353850e-7,
                8.74462839e-7,
                1.17634100e-6,
                9.86937310e-7,
                1.17634100e-6,
                8.74462839e-7,
            ],
            [
                8.74462839e-7,
                3.39353850e-7,
                0.00000000e-0,
                3.39353850e-7,
                8.74462839e-7,
                1.17634100e-6,
                9.86937310e-7,
                1.17634100e-6,
            ],
            [
                1.17634100e-6,
                8.74462839e-7,
                3.39353850e-7,
                0.00000000e-0,
                3.39353850e-7,
                8.74462839e-7,
                1.17634100e-6,
                9.86937310e-7,
            ],
            [
                9.86937310e-7,
                1.17634100e-6,
                8.74462839e-7,
                3.39353850e-7,
                0.00000000e-0,
                3.39353850e-7,
                8.74462839e-7,
                1.17634100e-6,
            ],
            [
                1.17634100e-6,
                9.86937310e-7,
                1.17634100e-6,
                8.74462839e-7,
                3.39353850e-7,
                0.00000000e-0,
                3.39353850e-7,
                8.74462839e-7,
            ],
            [
                8.74462839e-7,
                1.17634100e-6,
                9.86937310e-7,
                1.17634100e-6,
                8.74462839e-7,
                3.39353850e-7,
                0.00000000e-0,
                3.39353850e-7,
            ],
            [
                3.39353850e-7,
                8.74462839e-7,
                1.17634100e-6,
                9.86937310e-7,
                1.17634100e-6,
                8.74462839e-7,
                3.39353850e-7,
                0.00000000e-0,
            ],
        ]
    )

    model = DispersionModel(numbers, positions)

    res = model.get_pairwise_dispersion(
        DampingFunction(model="d4"),
        DampingParam(model="d4", s8=1.20065498, a1=0.40085597, a2=5.02928789)
    )

    assert approx(res.get("additive pairwise energy"), abs=thr) == pair_disp2
    assert (
        approx(res.get("non-additive pairwise energy"), abs=thr) == pair_disp3
    )


def test_properties() -> None:
    """Calculate dispersion related properties for a molecule"""
    thr = 1.0e-7

    numbers = np.array(4 * [7] + 12 * [1])
    positions = np.array(
        [
            [1.97420621099560, 1.97415497783241, 1.97424596974304],
            [6.82182427659395, 2.87346383480995, 7.72099517560089],
            [7.72104957181201, 6.82177051521773, 2.87336561318016],
            [2.87343220660781, 7.72108897828386, 6.82187093171878],
            [3.51863272100286, 2.63865333484548, 1.00652979981286],
            [2.63877594964754, 1.00647313885594, 3.51882748086447],
            [1.00639728563189, 3.51850454450845, 2.63869202592387],
            [8.36624975982697, 2.20896711017229, 8.68870955681018],
            [7.48639684558259, 3.84114715917956, 6.17640982573725],
            [5.85401675167715, 1.32911569888797, 7.05654606696031],
            [7.05646299938990, 5.85409590282274, 1.32879923864813],
            [8.68882633853582, 8.36611541129785, 2.20894120662207],
            [3.84121223226912, 6.17673669892998, 7.48629723649480],
            [1.32897854262127, 7.05658604099926, 5.85414031368096],
            [2.20884896069885, 8.68875820985799, 8.36643568423387],
            [6.17659142004652, 7.48627051643848, 3.84109594690835],
        ]
    )
    lattice = np.array(
        [
            [9.69523775911749, 0.00000000000000, 0.00000000000000],
            [0.00000000000000, 9.69523775911749, 0.00000000000000],
            [0.00000000000000, 0.00000000000000, 9.69523775911749],
        ]
    )
    cn = np.array(
        [
            2.5775156287150218,
            2.5775155620078536,
            2.5775157938667150,
            2.5775157704485387,
            0.8591731475439074,
            0.8591680526841657,
            0.8591744284869478,
            0.8591732359038715,
            0.8591678667769283,
            0.8591744593270527,
            0.8591684383407867,
            0.8591754863625011,
            0.8591751690053771,
            0.8591719636497469,
            0.8591686377934131,
            0.8591718691634261,
        ]
    )
    charges = np.array(
        [
            -0.86974543285813199,
            -0.86974501326130316,
            -0.86974658178316333,
            -0.86974643466661739,
            +0.28992010784471012,
            +0.28989847135893759,
            +0.28992738637932841,
            +0.28992042841052362,
            +0.28989746859382759,
            +0.28992753735979965,
            +0.28990032127437559,
            +0.28993196014690176,
            +0.28993044356403513,
            +0.28991421277784613,
            +0.28990119160907674,
            +0.28991393324985348,
        ]
    )
    alpha = np.array(
        [
            9.9853045768095097,
            9.9853025181287016,
            9.9853102162086920,
            9.9853094944097140,
            1.3513315505023076,
            1.3513939255952379,
            1.3513106323935196,
            1.3513306244831862,
            1.3513968091133417,
            1.3513101978580895,
            1.3513885997538553,
            1.3512974505020905,
            1.3513018164322617,
            1.3513485145877058,
            1.3513860914395939,
            1.3513493246534483,
        ]
    )

    model = DispersionModel(numbers, positions, lattice=lattice)

    res = model.get_properties()

    assert approx(res.get("coordination numbers"), abs=thr) == cn
    assert approx(res.get("partial charges"), abs=thr) == charges
    assert approx(res.get("polarizabilities"), abs=thr) == alpha


def test_hessian_pbed4() -> None:
    """Sanity test for numerical Hessian on a reasonable water structure."""
    thr = 1.0e-6  # numerical FD noise tolerance for symmetry

    numbers = np.array([8, 1, 1], dtype=int)
    positions = np.array(
        [
            [+0.00000000000000, +0.00000000000000, -0.73578586109551],
            [+1.44183152868459, +0.00000000000000, +0.36789293054775],
            [-1.44183152868459, +0.00000000000000, +0.36789293054775],
        ],
        dtype=float,
    )

    hess_ref = np.array(
        [
            [
                [
                    [ 1.43060192e-05,  0.00000000e+00, -4.92337901e-17],
                    [-7.15300960e-06,  0.00000000e+00, -2.56424952e-05],
                    [-7.15300960e-06,  0.00000000e+00,  2.56424952e-05],
                ],
                [
                    [ 0.00000000e+00, -5.26920204e-05,  0.00000000e+00],
                    [ 0.00000000e+00,  2.63460102e-05,  0.00000000e+00],
                    [ 0.00000000e+00,  2.63460102e-05,  0.00000000e+00],
                ],
                [
                    [ 6.77626358e-17,  0.00000000e+00,  1.63815257e-05],
                    [-3.66802752e-05,  0.00000000e+00, -8.19076293e-06],
                    [ 3.66802752e-05,  0.00000000e+00, -8.19076293e-06],
                ],
            ],
            [
                [
                    [-7.15300962e-06,  0.00000000e+00, -3.66802753e-05],
                    [ 1.86530626e-06,  0.00000000e+00,  3.11613853e-05],
                    [ 5.28770334e-06,  0.00000000e+00,  5.51889003e-06],
                ],
                [
                    [ 0.00000000e+00,  2.63460102e-05,  0.00000000e+00],
                    [ 0.00000000e+00, -2.14985217e-05,  0.00000000e+00],
                    [ 0.00000000e+00, -4.84748855e-06,  0.00000000e+00],
                ],
                [
                    [-2.56424952e-05,  0.00000000e+00, -8.19076283e-06],
                    [ 3.11613852e-05,  0.00000000e+00,  5.58414533e-06],
                    [-5.51889001e-06,  0.00000000e+00,  2.60661760e-06],
                ],
            ],
            [
                [
                    [-7.15300962e-06,  0.00000000e+00,  3.66802753e-05],
                    [ 5.28770334e-06,  0.00000000e+00, -5.51889003e-06],
                    [ 1.86530626e-06,  0.00000000e+00, -3.11613853e-05],
                ],
                [
                    [ 0.00000000e+00,  2.63460102e-05,  0.00000000e+00],
                    [ 0.00000000e+00, -4.84748855e-06,  0.00000000e+00],
                    [ 0.00000000e+00, -2.14985217e-05,  0.00000000e+00],
                ],
                [
                    [ 2.56424952e-05,  0.00000000e+00, -8.19076284e-06],
                    [ 5.51889001e-06,  0.00000000e+00,  2.60661760e-06],
                    [-3.11613852e-05,  0.00000000e+00,  5.58414533e-06],
                ],
            ],
        ], 
        dtype=float,
    )


    model = DispersionModel(numbers, positions, model="d4")
    damp = DampingFunction(model="d4")
    param = DampingParam(method="pbe", model="d4")

    res = model.get_numerical_hessian(damp, param)
    hess = res["hessian"]

    assert hess.shape == (3, 3, 3, 3)

    assert np.isfinite(hess).all()

    # symmetric under (i,a) <-> (j,b)
    assert approx(hess, abs=thr) == np.swapaxes(np.swapaxes(hess, 0, 2), 1, 3)

    assert approx(hess, abs=thr) == hess_ref


def test_error_model() -> None:
    """Test the error for unknown dispersion model"""
    numbers = np.array(
        [
            1,
            1,
            8,
        ]
    )
    positions = np.array(
        [
            [-0.02298820517725, 0.00000000000000, -1.76188954246096],
            [1.65369502723146, 0.00000000000000, 0.60848805100320],
            [-0.10273226709885, 0.00000000000000, 0.07266269355725],
        ]
    )

    with raises(ValueError) as exc:
        DispersionModel(numbers, positions, model="D42")

    assert "Unknown dispersion model" in str(exc)
