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
"""Tests for the qcelemental interface."""
from typing import Any

import numpy as np
import pytest

try:
    import qcelemental as qcel
    from qcelemental.models import AtomicInput, Molecule

    from tblite.qcschema import run_schema
except ModuleNotFoundError:
    qcel = None
    AtomicInput = dict
    Molecule = dict


@pytest.fixture
def multiplicity(request) -> int:
    """Multiplicity fixture."""
    return getattr(request, "param", 1)


@pytest.fixture(params=["ala-xab"])
def molecule(request, multiplicity: int) -> Molecule:
    """Get a molecule for testing."""
    if request.param == "ala-xab":
        return Molecule(
            symbols=list("NHCHCCHHHOCCHHHONHCHHH"),
            geometry=np.array(
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
            molecular_multiplicity = multiplicity,
        )

    raise ValueError(f"Unknown molecule: {request.param}")


@pytest.fixture(params=["energy", "gradient"])
def driver(request) -> str:
    """Driver fixture."""
    return request.param


@pytest.fixture(params=["GFN1-xTB", "GFN2-xTB"])
def method(request) -> str:
    """Method fixture."""
    return request.param


@pytest.fixture()
def atomic_input(molecule: Molecule, driver: str, method: str) -> AtomicInput:
    """AtomicInput fixture."""
    return AtomicInput(
        molecule=molecule,
        driver=driver,
        model={"method": method},
        keywords={"spin-polarization": 1.0},
    )


@pytest.fixture()
def return_result(molecule: Molecule, driver: str, method: str) -> Any:
    """Return result fixture."""
    if qcel is None:
        return None

    # fmt: off
    return {
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "energy",
            "GFN1-xTB",
        ): -34.98079481580359,
        (
            "65a3bab7309579268e383249bad8d33095c32db5", 
            "energy",
            "GFN1-xTB",
        ): -34.79969349420318,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "energy",
            "GFN2-xTB",
        ): -32.96247200752864,
        (
            "65a3bab7309579268e383249bad8d33095c32db5", 
            "energy",
            "GFN2-xTB",
        ): -32.79550096196652,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gradient",
            "GFN1-xTB",
        ): np.array(
            [
                [-7.8113833723432776e-3, 1.1927359496178194e-3, 4.5384293534289468e-3],
                [4.6431896826996466e-4, -1.1893457514353986e-3, -2.0196769165192010e-3],
                [1.4651707521730923e-3, -2.1095273552157868e-3, -1.9374661833396272e-3],
                [1.1712050542529663e-3, -6.2696965846751425e-4, 6.1050802435153369e-3],
                [-1.1270171972379136e-3, 9.4432659541408457e-4, -2.0020105757037662e-3],
                [1.1646562626373442e-2, -7.8838363893641971e-3, -9.4098748734209019e-3],
                [4.1886781466182168e-4, -2.0573790837090546e-4, 4.0632563116798048e-4],
                [-1.8053991102095574e-4, 1.1681264331977676e-3, 1.6828801638121846e-3],
                [-5.9900194915516749e-4, 3.4776846711688584e-4, -3.7491749453125091e-4],
                [-1.4520319022330924e-2, 6.0131467543009937e-3, 7.5687331646021375e-4],
                [1.5602405715421146e-2, 7.9473817513023640e-3, 8.0623888147124366e-3],
                [-3.2563959623975588e-4, -2.1680106787522012e-4, -8.8683653969162549e-4],
                [7.5180811527270801e-4, -1.9128778304917517e-4, -1.0174498970762392e-3],
                [7.4132920268697234e-4, -6.3819759962106609e-4, 5.1853393972177886e-4],
                [-7.5646645647444864e-4, 1.5490223231011606e-3, 6.5407919650053525e-4],
                [-9.0634683701016835e-3, -8.8383472527482337e-3, -1.2123112366846918e-2],
                [3.1541524835559096e-3, 2.7233221491356533e-3, 1.0127030544629243e-2],
                [-1.4266036482687263e-3, 1.2132816331079002e-3, -2.7113843360055362e-3],
                [2.0255265706568555e-4, -1.3134487584798992e-3, -7.9291928858605555e-4],
                [-2.3317056624669709e-3, -2.1233032658492240e-4, -1.1563077745307797e-3],
                [1.6004440124424719e-3, -3.2444847372739881e-3, 1.6000422202203360e-3],
                [9.2332778346354430e-4, 3.5712025321916925e-3, -1.9707177917017700e-5],
            ]
        ),
        (
            "65a3bab7309579268e383249bad8d33095c32db5",
            "gradient",
            "GFN1-xTB",
        ): np.array(
            [
                [5.8524614928954187e-3, -2.8273860862564008e-3, -2.5050614815354293e-2],
                [1.1339401892754533e-3, -3.8043823659953784e-3, 4.6415780745915643e-4],
                [6.1161622297322229e-4, 4.2685677262145104e-3, 1.4570851927090144e-3],
                [3.8586597902057181e-3, -2.5177483893864799e-3, 4.5515645848229879e-3],
                [-2.8918206405656940e-3, 2.8455866696309368e-3, -5.9668370428631338e-4],
                [-3.1705370345825074e-2, 2.4083262454045112e-2, 2.7384582155083294e-2],
                [-5.1707111952240632e-4, -2.5643418411340111e-4, 1.2581018565013965e-3],
                [-1.2470352203053496e-3, 2.5614875177244435e-3, 4.3056402494263653e-3],
                [-6.6844836248905481e-4, 2.8089761221454008e-4, 7.8809940986204145e-4],
                [6.2706334682231246e-2, -3.4003886955687332e-2, 8.0586129937534315e-4],
                [-9.4087011001570357e-3, -3.4501995852940378e-3, -1.7483270978818424e-4],
                [-9.2207770819250441e-4, -3.5421985555735337e-3, 9.5662473713666832e-4],
                [1.3078320300956743e-3, 4.1263243281132945e-5, -1.3500611258690431e-3],
                [-9.0178431942333195e-4, -2.0509954569158484e-3, 1.5334520471585486e-3],
                [-6.9184533697939866e-4, 2.1987490376191906e-3, 4.6577555215977903e-4],
                [9.1755875101464449e-3, 1.5842405697885673e-2, 1.4515316734675897e-2],
                [-3.8904615489211355e-2, -1.1230078292625741e-2, -3.0192692178494903e-2],
                [3.2109705087129807e-4, 5.0278073973231176e-3, -4.7478466044375703e-3],
                [1.6316664176902350e-3, 5.0025666147582287e-3, 3.8013881083669714e-3],
                [-2.8338123665745895e-3, 1.1958861407827702e-3, -6.8946026182785858e-4],
                [2.3970530203739990e-3, -3.0498463357009459e-3, 2.4297750760109019e-3],
                [1.6963336024870001e-3, 3.3846760960695773e-3, -1.9152334106901079e-3],
            ]
        ),
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gradient",
            "GFN2-xTB",
        ): np.array(
            [
                [-2.8793438861582673e-3, 3.6459441373180514e-3, 1.0760508724820901e-2],
                [4.4680466445768047e-4, -4.1388790826959765e-4, -1.6588285274993161e-3],
                [-5.2514242258322602e-3, -3.5790063091174704e-3, -3.7701182921630887e-3],
                [1.6245991142079565e-3, 3.3332711451623131e-4, 5.2523471739937032e-3],
                [2.7989793716838601e-4, 1.4306660544144156e-3, -2.4860639014568867e-3],
                [1.3079884332656911e-2, -6.5675701092584695e-3, -9.0219328076389001e-3],
                [-1.1472254803163923e-3, 5.7914197784694460e-5, 5.9834483223732073e-4],
                [-1.2270686170182344e-4, -3.1974725339932176e-4, 5.0166525086703872e-4],
                [1.4366367775558102e-4, 3.3058185471152853e-4, 7.6582251521247196e-4],
                [-1.5793967958079323e-2, 5.5332101481823806e-3, 1.7671361847792130e-4],
                [1.2223318558489288e-2, 6.4330353061192334e-3, 7.2760234115440180e-3],
                [-9.0557546481437291e-4, 2.4489454582212173e-5, -1.5354430992966309e-3],
                [1.0468034102157908e-3, 3.9865637021498047e-4, -1.2501442516385129e-4],
                [2.8779384356944992e-4, 8.3452297808966547e-5, 2.8882016450978579e-5],
                [-6.5407357546850360e-4, 1.9387903892084199e-4, 6.3097965388358664e-4],
                [-8.1387238298003383e-3, -8.7644052230312439e-3, -1.2729476394648730e-2],
                [1.0141709819820969e-2, -1.1848758166900799e-3, 5.2112574163926751e-3],
                [-9.1260711253248795e-4, 5.8340973460203481e-4, -2.0431325551924654e-3],
                [-4.4226928744326823e-3, 1.5005187195653046e-3, 1.5456820632938995e-3],
                [-1.1435369494648283e-3, -2.1402016535110218e-4, -6.6270629007887910e-4],
                [1.3141089852759119e-3, -1.3944191312867211e-3, 1.6665514558433831e-3],
                [7.8329387498344975e-4, 1.8888474876633428e-3, -3.8206183987922015e-4],
            ],
        ),
        (
            "65a3bab7309579268e383249bad8d33095c32db5", 
            "gradient",
            "GFN2-xTB",
        ): np.array(
            [
                [8.9058486746116817e-3, 9.5556561755216762e-4, -2.0209616087179871e-2],
                [1.2421526076339792e-3, -4.1923523693987094e-3, 5.7215715478311808e-4],
                [-2.2576776919733671e-3, 1.6558809431972221e-3, 8.0305415409530898e-4],
                [5.2249531879547508e-3, -1.1939559147596842e-3, 4.8043050324124694e-3],
                [-1.9781590191098415e-3, 1.3443336930798158e-3, -2.2004022007698592e-3],
                [-3.9736135444477710e-2, 2.1195956373900836e-2, 1.5245514578905904e-2],
                [-1.2242649160416944e-3, 2.2701509826147880e-4, 1.5636184382913569e-3],
                [-1.5356203639655174e-3, 1.1772789471972567e-3, 3.9239604716428842e-3],
                [-2.5686637138638515e-4, 6.5675170858745836e-4, 1.4060269271443855e-3],
                [5.7180507952153473e-2, -2.9814915578665419e-2, 3.5461024622819976e-3],
                [-1.6887863561331151e-2, -1.2557533741880246e-2, -1.2840551003415746e-2],
                [-1.2314264675434756e-3, -3.0608713297804904e-3, 1.5889359575349199e-3],
                [2.3018942346149881e-3, 9.7168844510480475e-4, -4.5113837420252472e-4],
                [-1.5507743234806266e-3, -1.9370017908039191e-3, 1.7872142236431562e-3],
                [-5.8814211436887976e-4, 1.1369180735641311e-3, 3.9542469547193223e-4],
                [1.4837352901349667e-2, 2.3174783926639412e-2, 2.4020007638463521e-2],
                [-2.3688848469882187e-2, -1.1167007508474183e-2, -2.4169418811976456e-2],
                [1.7136568601475088e-3, 3.5315178247889707e-3, -3.7553590160028993e-3],
                [-2.4670089463600148e-3, 6.9132702082702405e-3, 4.8206414352870413e-3],
                [-1.4685625067840512e-3, 7.8990871132563574e-4, -1.8312857979123655e-4],
                [1.5860108848151121e-3, -1.6803337279241715e-3, 2.2717947493028079e-3],
                [1.8789728934238306e-3, 1.8731023902173079e-3, -2.9391438459223070e-3],      
            ],
        ),
    }[(molecule.get_hash(), driver, method)]
    # fmt: on


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
@pytest.mark.parametrize("multiplicity", [1, 3], indirect=True)
def test_qcschema(atomic_input: AtomicInput, return_result: Any) -> None:
    """Test qcschema interface."""
    atomic_result = run_schema(atomic_input)

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result) == return_result


@pytest.fixture(params=[
      {"cpcm-solvation": 7.0}, 
      {"alpb-solvation": ["water", "bar1mol"]},
      {"gbsa-solvation": ["methanol", "reference"]},
      {"gbe-solvation": [7.0, "p16"]},
      {"gb-solvation": [7.0, "still"]},
   ])
def solvation(request) -> dict:
    """Solvation fixture."""
    return request.param


@pytest.fixture()
def atomic_input_solvation(molecule: Molecule, method: str, solvation: dict) -> AtomicInput:
    """AtomicInput fixture."""
    return AtomicInput(
        molecule=molecule,
        driver="energy",
        model={"method": method},
        keywords=solvation
    )


@pytest.fixture()
def return_result_solvation(molecule: Molecule, method: str, solvation: dict) -> Any:
    """Return result fixture."""
    if qcel is None:
        return None

    # fmt: off
    return {
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "cpcm-solvation",
            "GFN1-xTB",
        ): -34.9875514437647,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "cpcm-solvation",
            "GFN2-xTB",
        ): -32.9714541886084,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "alpb-solvation",
            "GFN1-xTB",
        ): -34.9966872968042,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "alpb-solvation",
            "GFN2-xTB",
        ): -32.9777888674693,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gbsa-solvation",
            "GFN1-xTB",
        ): -34.9989058302890,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gbsa-solvation",
            "GFN2-xTB",
        ): -32.9763255038126,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gbe-solvation",
            "GFN1-xTB",
        ): -34.9859026499920,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gbe-solvation",
            "GFN2-xTB",
        ): -32.9667327170214,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gb-solvation",
            "GFN1-xTB",
        ): -34.9859968024425,
        (
            "142dbe2f7f02c899c660c08ba85c086a366fbdec",
            "gb-solvation",
            "GFN2-xTB",
        ): -32.9668850575386,
    }[(molecule.get_hash(), list(solvation.keys())[0], method)]
    # fmt: on


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_qcschema_solvation(atomic_input_solvation: AtomicInput, return_result_solvation: Any) -> None:
    """Test qcschema interface."""
    atomic_result = run_schema(atomic_input_solvation)

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result) == return_result_solvation


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_unsupported_driver(molecule: Molecule):
    """Test unsupported driver name."""
    atomic_inp = AtomicInput(
        molecule=molecule,
        driver="hessian",
        model={"method": "GFN1-xTB"},
    )

    atomic_result = run_schema(atomic_inp)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert (
        "Driver 'hessian' is not supported by tblite."
        in atomic_result.error.error_message
    )


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_unsupported_method(molecule: Molecule):
    """Test unsupported method name."""
    atomic_inp = AtomicInput(
        molecule=molecule,
        driver="energy",
        model={"method": "GFN-xTB"},
    )

    atomic_result = run_schema(atomic_inp)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert (
        "Model 'GFN-xTB' is not supported by tblite."
        in atomic_result.error.error_message
    )


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_unsupported_basis(molecule: Molecule):
    """Test unsupported basis set."""
    atomic_inp = AtomicInput(
        molecule=molecule,
        driver="energy",
        model={"method": "GFN1-xTB", "basis": "def2-SVP"},
    )

    atomic_result = run_schema(atomic_inp)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert (
        "Basis sets are not supported by tblite."
        in atomic_result.error.error_message
    )


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_unsupported_keywords(molecule: Molecule):
    """Test unsupported keywords."""
    atomic_inp = AtomicInput(
        molecule=molecule,
        driver="gradient",
        model={"method": "GFN1-xTB"},
        keywords={"unsupported": True},
    )

    atomic_result = run_schema(atomic_inp)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "input_error"
    assert "Unknown keywords: unsupported" in atomic_result.error.error_message


@pytest.mark.skipif(qcel is None, reason="requires qcelemental")
def test_scf_not_converged(molecule: Molecule):
    """Test unconverged SCF."""
    atomic_inp = AtomicInput(
        molecule=molecule,
        driver="gradient",
        model={"method": "GFN1-xTB"},
        keywords={"max-iter": 3},
    )

    atomic_result = run_schema(atomic_inp)

    assert not atomic_result.success
    assert atomic_result.error.error_type == "execution_error"
    assert "SCF not converged in 3 cycles" in atomic_result.error.error_message
