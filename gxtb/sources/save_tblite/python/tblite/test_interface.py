# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.

from logging import Logger
from tempfile import NamedTemporaryFile

import numpy as np
from pytest import approx, raises

from tblite.exceptions import TBLiteRuntimeError, TBLiteValueError
from tblite.interface import Calculator, Result, symbols_to_numbers

THR = 1.0e-9


def get_ala(conf):
    """Retrieve a conformer of a small peptide"""
    numbers = np.array(
        [7, 1, 6, 1, 6, 6, 1, 1, 1, 8, 6, 6, 1, 1, 1, 8, 7, 1, 6, 1, 1, 1]
    )

    positions = {
        "xab": np.array(
            [
                [+2.65893135608838, -2.39249423371715, -3.66065400053935],
                [+3.49612941769371, -0.88484673975624, -2.85194146578362],
                [-0.06354076626069, -2.63180732150005, -3.28819116275323],
                [-1.07444177498884, -1.92306930149582, -4.93716401361053],
                [-0.83329925447427, -5.37320588052218, -2.81379718546920],
                [-0.90691285352090, -1.04371377845950, -1.04918016247507],
                [-2.86418317801214, -5.46484901686185, -2.49961410229771],
                [-0.34235262692151, -6.52310417728877, -4.43935278498325],
                [+0.13208660968384, -6.10946566962768, -1.15032982743173],
                [-2.96060093623907, +0.01043357425890, -0.99937552379387],
                [+3.76519127865000, -3.27106236675729, -5.83678272799149],
                [+6.47957316843231, -2.46911747464509, -6.21176914665408],
                [+7.32688324906998, -1.67889171278096, -4.51496113512671],
                [+6.54881843238363, -1.06760660462911, -7.71597456720663],
                [+7.56369260941896, -4.10015651865148, -6.82588105651977],
                [+2.64916867837331, -4.60764575400925, -7.35167957128511],
                [+0.77231592220237, -0.92788783332000, +0.90692539619101],
                [+2.18437036702702, -2.20200039553542, +0.92105755612696],
                [+0.01367202674183, +0.22095199845428, +3.27728206652909],
                [+1.67849497305706, +0.53855308534857, +4.43416031916610],
                [-0.89254709011762, +2.01704896333243, +2.87780123699499],
                [-1.32658751691561, -0.95404596601807, +4.30967630773603],
            ]
        ),
        "xac": np.array(
            [
                [+0.89812994422336, -3.49785066920537, -4.04176475935607],
                [+0.22096411844816, -4.86653879303702, -5.16965309356488],
                [+1.50858972981055, -4.22297129418255, -1.44296740345170],
                [+1.26381263207546, -6.26173130725955, -1.36929929792062],
                [+4.22584479532188, -3.58094126128602, -0.72765628187693],
                [-0.36535595783691, -3.23105046121391, +0.53224473210756],
                [+4.59700287054853, -4.26587358691560, +1.17638157393339],
                [+5.53245639246579, -4.49962127256050, -2.02509339376622],
                [+4.54878065243139, -1.55600067058432, -0.78659473595735],
                [-1.40338188617737, -4.70093312894944, +1.99692732352687],
                [+1.19184541389054, -1.17965860262351, -5.08564964981169],
                [+0.66646223547991, -1.02784609220985, -7.88469520699739],
                [+2.46126071643208, -0.84993899853549, -8.87377461956087],
                [-0.42509529147581, +0.67053350705954, -8.25313731918547],
                [-0.33603421332938, -2.66534515804911, -8.61567815458241],
                [+1.87087395689416, +0.69377329947433, -3.88238849685894],
                [-0.66896549096126, -0.70041907517856, +0.57158466492815],
                [+0.24390855421877, +0.32329261951026, -0.76688726802653],
                [-2.34873581523733, +0.49815011834652, +2.36364902826355],
                [-2.51669513442034, -0.71675917183750, +4.00880640863199],
                [-1.57466603805426, +2.31572781563081, +2.92500014923521],
                [-4.23036142386086, +0.78917833057040, +1.57683879578790],
            ]
        ),
        "xag": np.array(
            [
                [-1.98293112719749, -0.68900404398734, +0.37628864997560],
                [-1.54035272113532, -2.33144052962209, +1.23491020893064],
                [+0.09692928120200, +0.99422262169598, -0.11989378229633],
                [-0.09901909722504, +1.72588002674889, -2.03875555767431],
                [+0.14141827166397, +3.20957080884314, +1.73528063332959],
                [+2.50745780392154, -0.55594193008926, +0.12427542626420],
                [+1.69093368632231, +4.49899142526548, +1.31794453504710],
                [-1.63338894984043, +4.23428690562497, +1.59126592199206],
                [+0.36779723029067, +2.50163220083067, +3.65420914638345],
                [+2.56229717302803, -2.53711155537531, +1.32871903135090],
                [-4.36159407476856, -0.12962702640436, -0.39562772950648],
                [-6.34100574546165, -2.09075933777830, +0.21330645338918],
                [-7.23610050442420, -2.68187533006900, -1.53914540310063],
                [-7.78623541548118, -1.21956781217005, +1.38592307383011],
                [-5.57984815064465, -3.73504338108486, +1.18089285147183],
                [-4.86717200167912, +1.83918213313466, -1.51162714188889],
                [+4.55965923505744, +0.48183329726772, -0.99463516518498],
                [+4.30889971335617, +2.04332235450700, -2.04508383499641],
                [+7.03944933546151, -0.66526542357055, -0.85849489668778],
                [+8.44432563879401, +0.71932844719313, -0.28192138996626],
                [+7.59387580143664, -1.47860319740057, -2.66533359746007],
                [+6.94508476368302, -2.16228363458999, +0.53955258865532],
            ]
        ),
    }[conf]

    return (numbers, positions)


def get_crcp2():
    """Get structure for CrCP2"""
    numbers = np.array(
        [24, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 6, 6, 6, 1, 6, 1, 6, 1, 1, 1]
    )
    positions = np.array(
        [
            [+0.00000000000000, +0.00000000000000, -0.06044684528305],
            [+0.00000000000000, +3.19613712523833, +2.30877824528580],
            [+2.18828801115897, +3.32943780995850, +0.70249948585735],
            [+1.33235791539260, +3.55640652898451, -1.83908673090077],
            [-1.33235791539260, +3.55640652898451, -1.83908673090077],
            [-2.18828801115897, +3.32943780995850, +0.70249948585735],
            [+0.00000000000000, +3.10509505378016, +4.34935395653655],
            [+4.13810718850644, +3.28428734944129, +1.31235006648465],
            [+2.52190264478215, +3.60569548880831, -3.50208900904436],
            [-2.52190264478215, +3.60569548880831, -3.50208900904436],
            [-4.13810718850644, +3.28428734944129, +1.31235006648465],
            [+2.18828801115897, -3.32943780995850, +0.70249948585735],
            [+0.00000000000000, -3.19613712523833, +2.30877824528580],
            [+1.33235791539260, -3.55640652898451, -1.83908673090077],
            [+4.13810718850644, -3.28428734944129, +1.31235006648465],
            [-2.18828801115897, -3.32943780995850, +0.70249948585735],
            [+0.00000000000000, -3.10509505378016, +4.34935395653655],
            [-1.33235791539260, -3.55640652898451, -1.83908673090077],
            [+2.52190264478215, -3.60569548880831, -3.50208900904436],
            [-4.13810718850644, -3.28428734944129, +1.31235006648465],
            [-2.52190264478215, -3.60569548880831, -3.50208900904436],
        ]
    )
    return (numbers, positions)


def test_gfn1():
    """Basic test for GFN1-xTB method"""
    gradient = np.array(
        [
            [-2.30890577e-5, -8.23027141e-3, +3.83526103e-3],
            [+8.93064693e-5, +3.76305314e-3, -3.16225182e-4],
            [+2.61812066e-3, +1.43734260e-3, +1.69221930e-3],
            [-1.30022454e-3, +7.36897334e-3, +1.44359655e-3],
            [-1.18095082e-3, -5.07048042e-4, +7.16557787e-4],
            [+1.03129859e-2, -3.22454163e-3, -1.17868424e-2],
            [+2.25696958e-4, -1.52527805e-4, -6.35931277e-4],
            [-2.48556212e-4, +4.95584663e-5, +8.26541624e-4],
            [+1.08580370e-4, -2.32855279e-3, -1.64881140e-4],
            [-7.27516746e-3, -1.08743381e-2, +1.02939725e-2],
            [-6.39918988e-3, -8.43381632e-3, -1.48319990e-2],
            [+8.78977536e-4, +8.30495556e-4, +6.52170578e-4],
            [-1.14377340e-3, -1.13193865e-4, -1.94865181e-4],
            [+8.78368248e-4, -1.62570375e-3, +2.95802372e-4],
            [+4.23367378e-4, +1.05983147e-3, -9.63718801e-4],
            [+3.88369600e-3, +1.27236421e-2, +8.62989820e-3],
            [-3.97039649e-3, +1.20670709e-2, +3.77235583e-3],
            [+1.86128358e-3, -3.42859740e-3, -1.71401975e-3],
            [+5.09298024e-4, -6.34874352e-4, -1.43338011e-3],
            [+1.43706589e-5, +2.39757682e-3, -2.91623626e-3],
            [-1.71237422e-3, -2.17516013e-3, -4.94411233e-4],
            [+1.44967025e-3, +3.10811194e-5, +3.29413462e-3],
        ]
    )

    numbers, positions = get_ala("xab")

    calc = Calculator("GFN1-xTB", numbers, positions)
    calc.set("accuracy", 1.0)

    res = calc.singlepoint()

    assert res.get("energy") == approx(-34.980794815805446, abs=THR)

    numbers, positions = get_ala("xac")
    calc.update(positions)
    calc.set("save-integrals", 1)

    res = calc.singlepoint(res)

    assert res.get("energy") == approx(-34.987786081514066, abs=THR)
    assert res.get("gradient") == approx(gradient, abs=THR)

    res = res.dict()
    assert "density-matrix" in res
    assert "overlap-matrix" in res
    assert "hamiltonian-matrix" in res


def test_gfn2():
    """Basic test for GFN2-xTB method"""
    gradient = np.array(
        [
            [+1.17642910e-2, -5.29400903e-3, +3.72680358e-3],
            [-1.36717977e-3, +8.73017982e-4, -1.58097215e-3],
            [-2.28602657e-3, +3.05314211e-3, +1.29503272e-3],
            [+1.17878429e-3, -1.62260139e-3, +3.54659145e-3],
            [-1.40782520e-3, -1.62625462e-3, -3.06664234e-3],
            [-8.96440228e-3, +1.19118693e-2, -8.33051410e-3],
            [+9.33591498e-4, +1.20346699e-3, +1.89444832e-4],
            [-2.18294922e-4, +3.93504056e-4, -1.53811277e-4],
            [+6.42406310e-4, -7.64866201e-4, +3.03844581e-4],
            [-3.36497732e-4, -1.45985780e-2, +9.86968946e-3],
            [-8.97934174e-4, -1.15903933e-2, +5.90534278e-3],
            [+4.32454422e-4, +1.55768098e-3, -6.84108710e-4],
            [+4.09031189e-4, +1.48996321e-4, +2.50460899e-4],
            [+4.41135937e-4, -1.20304491e-4, -1.06933538e-4],
            [-1.11124347e-3, -1.23817267e-3, +4.80262746e-4],
            [-5.50236606e-3, +1.46035617e-2, -8.54857472e-3],
            [+5.06284695e-3, +7.13309791e-3, -5.28051609e-3],
            [-2.35509860e-3, -1.78589715e-3, +1.35728142e-3],
            [+3.24958615e-3, -2.27290642e-3, +5.08621002e-4],
            [-5.74390261e-4, -1.73044948e-3, -4.21845336e-4],
            [+1.59689965e-4, +8.11151135e-5, +2.26546673e-3],
            [+7.47441343e-4, +1.68498046e-3, -1.52492395e-3],
        ]
    )
    # fmt: off
    shell_map = np.array(
        [
            *[0, 0, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 8, 9, 9],
            *[10, 10, 11, 11, 12, 13, 14, 15, 15, 16, 16, 17, 18, 18, 19, 20, 21],
        ]
    )
    angular_momenta = np.array(
        [
            *[0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1],
            *[0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0],
        ]
    )
    orbital_map = np.array(
        [
            *[0, 1, 1, 1, 2, 3, 4, 4, 4, 5, 6, 7, 7, 7, 8, 9, 9, 9],
            *[10, 11, 12, 13, 14, 14, 14, 15, 16, 16, 16, 17, 18, 18, 18, 19],
            *[20, 21, 22, 23, 23, 23, 24, 25, 25, 25, 26, 27, 28, 28, 28, 29],
            *[30, 31],
        ]
    )
    # fmt: on

    numbers, positions = get_ala("xac")

    calc = Calculator("GFN2-xTB", numbers, positions)

    assert all(shell_map == calc.get("shell-map"))
    assert all(angular_momenta == calc.get("angular-momenta"))
    assert all(orbital_map == calc.get("orbital-map"))

    calc.set("accuracy", 1.0)
    res = Result()

    res = calc.singlepoint(res, copy=True)

    assert res.get("energy") == approx(-32.97134042392298, abs=THR)

    numbers, positions = get_ala("xag")
    calc.update(positions)

    calc.singlepoint(res)

    assert res.get("energy") == approx(-32.97132127543436, abs=THR)
    assert res.get("energy") == approx(sum(res.get("energies")), abs=THR)
    assert res.get("gradient") == approx(gradient, abs=THR)


def test_ipea1():
    """Basic test for IPEA1-xTB method"""
    gradient = np.array(
        [
            [+8.66895375e-3, -8.77412044e-3, +2.97823240e-3],
            [+9.52587611e-4, +3.15940527e-5, -3.12816500e-4],
            [-2.57987382e-3, +3.84961559e-3, +2.77298284e-3],
            [-4.81635429e-4, -2.88802818e-4, +8.40482914e-4],
            [+1.24799170e-3, +3.74900057e-3, +1.76205872e-4],
            [-9.62505451e-3, +9.63955941e-3, -6.95804474e-3],
            [-1.33281533e-3, +8.85183104e-4, +2.69619028e-3],
            [+1.47586309e-3, -1.46453763e-3, -4.78558561e-4],
            [+3.23497148e-4, -2.23069073e-5, -1.99573842e-3],
            [+8.76375745e-4, -1.50759679e-2, +1.07645733e-2],
            [-8.61807851e-3, -9.13562754e-3, +3.40925515e-3],
            [+1.15849566e-3, -6.86610262e-4, +1.00471922e-3],
            [+3.76866861e-4, +1.53381059e-3, +2.32381437e-3],
            [+1.33910593e-3, -1.05118440e-3, -2.32355779e-3],
            [-2.15237913e-3, +3.54288137e-4, -4.07100491e-4],
            [-4.44678667e-3, +1.32270519e-2, -8.43304298e-3],
            [+9.77389138e-3, +9.50568602e-3, -8.40237784e-3],
            [-2.32540898e-3, -4.14061089e-3, +2.43721523e-3],
            [+7.22194104e-3, -1.14254387e-3, -1.31678735e-3],
            [-1.07808668e-3, -2.79837461e-3, -3.51620344e-4],
            [+4.75694235e-4, -1.55429929e-4, +3.34209540e-3],
            [-1.25114510e-3, +1.96032776e-3, -1.76612193e-3],
        ]
    )

    numbers, positions = get_ala("xab")

    calc = Calculator("IPEA1-xTB", numbers, positions)
    for key, value in {"accuracy": 1.0, "verbosity": 2}.items():
        calc.set(key, value)

    res = calc.singlepoint()

    assert res.get("energy") == approx(-38.40436019312474, abs=THR)
    assert res.get("energy") == res["energy"]

    numbers, positions = get_ala("xag")
    calc.update(positions)

    res = calc.singlepoint(res).dict()

    assert res.get("energy") == approx(-38.40675517188822, abs=THR)
    assert res.get("gradient") == approx(gradient, abs=THR)


def test_gfn2_mindless():
    """Test GFN2-xTB for a mindless molecule"""
    numbers = np.array([1, 1, 6, 5, 1, 15, 8, 17, 13, 15, 5, 1, 9, 15, 1, 15])
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
    calc = Calculator("GFN2-xTB", numbers, positions)
    res = calc.singlepoint()

    assert res.get("energy") == approx(-31.716158891203904, abs=THR)

    res2 = Result()
    with NamedTemporaryFile(suffix=".npz") as tmp:
        res.save(tmp.name)
        del res
        res2.load(tmp.name)

    calc.set("max-iter", 3)
    res = calc.singlepoint(res2)

    assert res.get("energy") == approx(-31.716158891203904, abs=THR)


def test_ipea1_charge():
    """Test charged systems in an IP calculation"""
    numbers = np.array([14, 1, 1, 1, 1])
    positions = np.array(
        [
            [+0.00000000000000, +0.00000000000000, +0.00000000000000],
            [+1.61972522566005, -1.61972522566005, +1.61972522566005],
            [-1.61972522566005, +1.61972522566005, +1.61972522566005],
            [-1.61972522566005, -1.61972522566005, -1.61972522566005],
            [+1.61972522566005, +1.61972522566005, -1.61972522566005],
        ]
    )
    calc = Calculator("IPEA1-xTB", numbers, positions)

    res = calc.singlepoint()
    assert res.get("energy") == approx(-4.1365706907641115, abs=THR)

    positions = np.array(
        [
            [-0.00115627677826, -1.28225218579877, -0.61598185229325],
            [+0.00700035129222, +2.20006589600569, -1.15896195249942],
            [-0.00053638554336, +2.26215432948848, +0.32516560218069],
            [+2.43459795927086, -1.59621225620704, +0.72661653583618],
            [-2.43990564824148, -1.58375578348835, +0.72316166677581],
        ]
    )
    calc = Calculator("IPEA1-xTB", numbers, positions, charge=1.0, uhf=1)

    res = calc.singlepoint(res)
    assert res.get("energy") == approx(-3.5168511494225148, abs=THR)


def test_spgfn1():
    """Test SP GFN1-xTB"""
    numbers, positions = get_crcp2()
    calc = Calculator("GFN1-xTB", numbers, positions)

    ls_energy = calc.singlepoint().get("energy")

    calc.add("spin-polarization")
    ls_energy_sp = calc.singlepoint().get("energy")
    assert ls_energy_sp == approx(ls_energy)
    assert ls_energy_sp == approx(-28.349613833732931)

    calc.update(uhf=2)
    hs_energy_sp = calc.singlepoint().get("energy")
    assert hs_energy_sp == approx(-28.37648585236444)


def test_spgfn1_densities():
    """Test density matrices for SP GFN1-xTB"""
    numbers = np.array([1, 8])
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ]
    )
    calc = Calculator("GFN1-xTB", numbers, positions)
    calc.add("spin-polarization")
    calc.set("save-integrals", 1)

    res = calc.singlepoint()

    s = res.get("overlap-matrix")
    pa, pb = res.get("density-matrix")

    assert np.sum(pa * s) == approx(4.0)
    assert np.sum(pb * s) == approx(3.0)


def test_spgfn1_orbital_energies():
    """Test orbital energies for SP GFN1-xTB"""
    numbers, positions = get_crcp2()
    calc = Calculator("GFN1-xTB", numbers, positions)

    orben = calc.singlepoint().get("orbital-energies")

    calc.add("spin-polarization")
    orben_a, orben_b = calc.singlepoint().get("orbital-energies")
    assert orben_a == approx(orben)
    assert orben_b == approx(orben)


def test_spgfn1_orbital_occupations():
    """Test orbital occupations for SP GFN1-xTB"""
    numbers, positions = get_crcp2()
    calc = Calculator("GFN1-xTB", numbers, positions)

    occs = calc.singlepoint().get("orbital-occupations")

    calc.add("spin-polarization")
    occs_a, occs_b = calc.singlepoint().get("orbital-occupations")
    assert occs_a == approx(0.5 * occs, abs=THR)
    assert occs_b == approx(0.5 * occs, abs=THR)


def test_spgfn1_orbital_occupations_and_coefficients():
    """Test orbital occupations and coefficients for SP GFN1-xTB"""
    numbers, positions = get_crcp2()
    calc = Calculator("GFN1-xTB", numbers, positions)
    calc.add("spin-polarization")

    res = calc.singlepoint()

    occs_a, occs_b = res.get("orbital-occupations")
    ca, cb = res.get("orbital-coefficients")

    pa, pb = res.get("density-matrix")

    pa_reconstruct = np.einsum("k,ik,jk->ij", occs_a, ca, ca)
    pb_reconstruct = np.einsum("k,ik,jk->ij", occs_b, cb, cb)

    assert pa_reconstruct == approx(pa, abs=THR)
    assert pb_reconstruct == approx(pb, abs=THR)


def test_post_processing_api():
    """Test post-processing API"""
    numbers, positions = get_crcp2()
    calc = Calculator("GFN1-xTB", numbers, positions)
    calc.add("bond-orders")
    res = calc.singlepoint()
    with raises(
        TBLiteValueError,
        match="Molecular dipole was not calculated. By default it is computed.",
    ):
        res.get("dipole")

    with raises(
        TBLiteValueError,
        match="Molecular quadrupole was not calculated. By default it is computed.",
    ):
        res.get("quadrupole")

    res.dict()

    calc = Calculator("GFN1-xTB", numbers, positions)
    calc.add("molecular-multipoles")
    res = calc.singlepoint()
    with raises(
        TBLiteValueError,
        match="Bond-orders were not calculated. By default they are computed.",
    ):
        res.get("bond-orders")

    calc = Calculator("GFN1-xTB", numbers, positions)
    calc.add("spin-polarization")

    wbo_sp = calc.singlepoint().get("bond-orders")
    assert wbo_sp.ndim == 3

def test_xtbml_api():
    numbers, positions = get_crcp2()
    calc = Calculator("GFN1-xTB", numbers, positions)
    calc.add("xtbml")
    res = calc.singlepoint()
    dict_ = res.dict()
    cn = dict_.get("post-processing-dict")
    cn = cn.get("CN_A")
    assert len(cn) == numbers.size

    dict_xtbml = res.dict()
    dict_xtbml = dict_xtbml.get("post-processing-dict")
    assert len(dict_xtbml.keys()) == 38
    
    toml_inp = ["[post-processing.xtbml] \n",
                "geometry = false \n",
                "density = true \n",
                "orbital = true \n",
                "energy = false \n", 
                "convolution = true\n",
                "a = [1.0, 1.2]\n",
                "tensorial-output = true"]
    
    with open("xtbml.toml", "w") as f:
        f.writelines(toml_inp)
        
    calc = Calculator("GFN1-xTB", numbers, positions)
    calc.add("xtbml.toml")
    
    res = calc.singlepoint()
    dict_ = res.get("post-processing-dict")
    
    assert dict_.get("CN_A") is None
    
    assert len(dict_.keys()) == 129
    

def test_solvation_gfn2_cpcm():
    """Test CPCM solvation with GFN2-xTB"""
    numbers, positions = get_crcp2()

    calc = Calculator("GFN2-xTB", numbers, positions)
    calc.set("accuracy", 1.0)
    calc.add("cpcm-solvation", 7.0)

    energy = calc.singlepoint().get("energy")
    assert energy == approx(-28.43287176929, abs=THR)


def test_solvation_gfn2_alpb():
    """Test ALPB solvation with GFN2-xTB"""
    numbers, positions = get_crcp2()

    calc = Calculator("GFN2-xTB", numbers, positions)
    calc.set("accuracy", 1.0)
    calc.add("alpb-solvation", "ethanol")

    energy = calc.singlepoint().get("energy")
    assert energy == approx(-28.448543424860, abs=THR)


def test_solvation_gfn2_alpb_bar1mol():
    """Test ALPB solvation with GFN2-xTB with state"""
    numbers, positions = get_crcp2()

    calc = Calculator("GFN2-xTB", numbers, positions)
    calc.set("accuracy", 1.0)
    calc.add("alpb-solvation", "ethanol", "bar1mol")

    energy = calc.singlepoint().get("energy")
    assert energy == approx(-28.445512192034, abs=THR)


def test_solvation_gfn2_alpb_reference():
    """Test ALPB solvation with GFN2-xTB"""
    numbers, positions = get_crcp2()

    calc = Calculator("GFN2-xTB", numbers, positions)
    calc.set("accuracy", 1.0)
    calc.add("alpb-solvation", "ethanol", "reference")

    energy = calc.singlepoint().get("energy")
    assert energy == approx(-28.442829777462, abs=THR)


def test_solvation_gfn1_gbe():
    """Test GBE solvation with GFN1-xTB"""
    numbers, positions = get_crcp2()

    calc = Calculator("GFN1-xTB", numbers, positions)
    calc.set("accuracy", 1.0)
    calc.add("gbe-solvation", 7.0, "p16")

    energy = calc.singlepoint().get("energy")
    assert energy == approx(-28.35002525960, abs=THR)


def test_solvation_gfn2_gbsa():
    """Test GBSA solvation with GFN2-xTB"""
    numbers, positions = get_crcp2()

    calc = Calculator("GFN2-xTB", numbers, positions)
    calc.set("accuracy", 1.0)
    calc.add("gbsa-solvation", "water")

    energy = calc.singlepoint().get("energy")
    assert energy == approx(-28.439916768534, abs=THR)


def test_solvation_gfn2_gb():
    """Test GB solvation with GFN2-xTB"""
    numbers, positions = get_crcp2()

    calc = Calculator("GFN2-xTB", numbers, positions)
    calc.set("accuracy", 1.0)
    calc.add("gb-solvation", 7.0, "still")

    energy = calc.singlepoint().get("energy")
    assert energy == approx(-28.43676829542, abs=THR) 

def test_solvation_gfn1_alpb():
    """Test ALPB solvation with GFN1-xTB"""
    numbers, positions = get_crcp2()

    calc = Calculator("GFN1-xTB", numbers, positions)
    calc.set("accuracy", 1.0)
    calc.add("alpb-solvation", "ethanol", "reference")

    energy = calc.singlepoint().get("energy")
    assert energy == approx(-28.3545598426876, abs=THR)


def test_result_getter():
    """Check error handling in result container getter"""

    res = Result()

    with raises(RuntimeError, match="Result does not contain"):
        res.get("energy")

    with raises(RuntimeError, match="Result does not contain"):
        res.get("gradient")

    with raises(RuntimeError, match="Result does not contain"):
        res.get("virial")

    with raises(ValueError, match="Attribute 'unknown' is not available"):
        res.get("unknown")
    


def test_result_setter():
    """Check error handling in result container setter"""

    res = Result()

    with raises(ValueError, match="Attribute 'unknown' cannot be set"):
        res.set("unknown", 1.0)

    with raises(ValueError, match="Attribute 'unknown' cannot be set"):
        res["unknown"] = 1.0


def test_unknown_method():
    """Check handling of non-existing methods"""

    numbers, positions = get_ala("xab")
    with raises(ValueError, match="Method 'GFN-xTB' is not available"):
        Calculator("GFN-xTB", numbers, positions)


def test_unknown_attribute():
    """Check handling of setting unknown attributes"""

    numbers, positions = get_ala("xab")
    calc = Calculator("GFN2-xTB", numbers, positions)

    with raises(ValueError, match="Attribute 'unknown' is not supported"):
        calc.set("unknown", 1.0)


def test_gfn1_logging():
    """Basic test for GFN1-xTB method"""
    numbers, positions = get_ala("xab")

    logger = Logger("test")

    calc = Calculator(
        "GFN1-xTB", numbers, positions, color=False, logger=logger.info
    )
    res = calc.singlepoint()

    assert res.get("energy") == approx(-34.980794815805446, abs=THR)

    def broken_logger(message: str) -> None:
        raise NotImplementedError("This logger is broken")

    calc = Calculator(
        "GFN1-xTB", numbers, positions, color=False, logger=broken_logger
    )
    with raises(TBLiteRuntimeError):
        calc.singlepoint()


def test_symbols():
    """Check initialization with element symbols"""
    symbols = ["Si", "H", "H", "H", "H"]
    positions = np.array(
        [
            [+0.00000000000000, +0.00000000000000, +0.00000000000000],
            [+1.61972522566005, -1.61972522566005, +1.61972522566005],
            [-1.61972522566005, +1.61972522566005, +1.61972522566005],
            [-1.61972522566005, -1.61972522566005, -1.61972522566005],
            [+1.61972522566005, +1.61972522566005, -1.61972522566005],
        ]
    )

    calc = Calculator("GFN2-xTB", symbols_to_numbers(symbols), positions)
    res = calc.singlepoint()

    assert res.get("energy") == approx(-3.763120637211, abs=THR)


def test_numbers():
    """Check initialization with atomic numbers"""
    numbers = [14, 1, 1, 1, 1]
    positions = np.array(
        [
            [+0.00000000000000, +0.00000000000000, +0.00000000000000],
            [+1.61972522566005, -1.61972522566005, +1.61972522566005],
            [-1.61972522566005, +1.61972522566005, +1.61972522566005],
            [-1.61972522566005, -1.61972522566005, -1.61972522566005],
            [+1.61972522566005, +1.61972522566005, -1.61972522566005],
        ]
    )

    calc = Calculator("GFN2-xTB", numbers, positions)
    res = calc.singlepoint()

    assert res.get("energy") == approx(-3.763120637211, abs=THR)
