import re
import sys
import traceback
from dataclasses import dataclass
from typing import Any, Dict, Tuple, Optional

if sys.version_info >= (3, 8):
    from typing import Literal, Protocol
else:
    from typing_extensions import Literal, Protocol


# ======================================================================================
@dataclass
class MatchResult:
    status: Literal["OK", "WRONG RESULT", "N/A"]
    error: Optional[str]
    value: Optional[float]


# ======================================================================================
def run_matcher(output: str, **spec: Any) -> MatchResult:
    try:
        kwargs = spec.copy()
        key = kwargs.pop("matcher")
        return registry[key].run(output, **kwargs)
    except:
        return MatchResult("N/A", error=traceback.format_exc(), value=None)


# ======================================================================================
class Matcher(Protocol):
    def run(self, output: str, **kwargs: Any) -> MatchResult: ...


# ======================================================================================
class GenericMatcher(Matcher):
    def __init__(self, pattern: str, col: int):
        self.pattern = pattern
        for c in r"[]()|+*?":
            pattern = pattern.replace(c, f"\\{c}")  # escape special chars
        pattern = pattern.replace("  ", r"\s\s+")
        self.regex = re.compile(pattern)
        self.col = col

    def run(self, output: str, **kwargs: Any) -> MatchResult:
        tol, ref = kwargs["tol"], kwargs["ref"]
        assert isinstance(tol, float) or isinstance(ref, int)
        assert isinstance(ref, float) or isinstance(ref, int)
        # grep result
        for line in reversed(output.split("\n")):
            if self.regex.search(line):
                value_str = line.split()[self.col - 1]
                break
        else:
            error = f"Result not found: '{self.pattern}'.\n"
            return MatchResult("WRONG RESULT", error, value=None)

        # parse result
        try:
            value = float(value_str)
        except:
            error = f"Could not parse result as float: '{value_str}'.\n"
            return MatchResult("WRONG RESULT", error, value=None)

        # compare result to reference
        diff = value - ref
        rel_error = abs(diff / ref if ref != 0.0 else diff)
        if rel_error > tol:
            error = f"Difference too large: {rel_error:.2e} > {tol}, value: {value}.\n"
            return MatchResult("WRONG RESULT", error, value)

        return MatchResult("OK", error=None, value=value)  # passed


# ======================================================================================
class MatcherRegistry(Dict[str, Matcher]):
    def __setitem__(self, key: str, value: Matcher) -> None:
        assert key not in self  # check for name collisions
        super().__setitem__(key, value)


# ======================================================================================
registry = MatcherRegistry()

# Total energy in Hartree
registry["E_total"] = GenericMatcher(r"Total energy:", col=3)

registry["M002"] = GenericMatcher(r"MD| Potential energy", col=5)
registry["M003"] = GenericMatcher(r"Total energy [eV]:", col=4)
registry["M004"] = GenericMatcher(r"Ideal and single determinant", col=8)
registry["M005"] = GenericMatcher(r"BSSE-free interaction energy:", col=5)
registry["M006"] = GenericMatcher(r"Average Energy", col=4)
registry["M007"] = GenericMatcher(r"OPT| Total energy [hartree]", col=5)
registry["M008"] = GenericMatcher(r"VIB|Frequency", col=3)
registry["M009"] = GenericMatcher(r"PINT| Total energy =", col=5)
registry["M010"] = GenericMatcher(r"BAND TOTAL ENERGY [au]", col=6)
registry["M011"] = GenericMatcher(r"ENERGY| Total FORCE_EVAL", col=9)
registry["M012"] = GenericMatcher(r"B2(T) =", col=4)
registry["M013"] = GenericMatcher(r"sparseness function f2 =", col=5)
registry["M014"] = GenericMatcher(r"CheckSum Shifts =", col=4)
registry["M015"] = GenericMatcher(r"CheckSum NICS =", col=4)
registry["M016"] = GenericMatcher(r"CheckSum Chi =", col=4)
registry["M017"] = GenericMatcher(r"Total=", col=8)
registry["M018"] = GenericMatcher(r"MS| TRACKED FREQUENCY", col=6)
registry["M019"] = GenericMatcher(r"CheckSum splines =", col=4)
registry["M020"] = GenericMatcher(r"epr|TOT:checksum", col=2)
registry["M021"] = GenericMatcher(r"VIB|", col=3)
registry["M022"] = GenericMatcher(r"PW exchange energy", col=5)
registry["M023"] = GenericMatcher(r"Total Spread", col=5)
registry["M024"] = GenericMatcher(r"DFT+U energy:", col=3)
registry["M025"] = GenericMatcher(r"OPT| Local curvature", col=4)
registry["M026"] = GenericMatcher(r"Total electronic charge (G-space):", col=6)
registry["M027"] = GenericMatcher(r"Total electronic charge (R-space):", col=5)
registry["M028"] = GenericMatcher(r"Total spin density (R-space)           :", col=6)
registry["M029"] = GenericMatcher(r"Heat of formation [kcal/mol]:", col=5)

# HOMO-LUMO gap / bandgap at Gamma-point of a DFT calculation
registry["E_gap_DFT"] = GenericMatcher(r"HOMO - LUMO gap [eV]", col=7)
registry["E_gap_DFT_2"] = GenericMatcher(r"Band gap:", col=4)  # TODO merge with prev.

registry["M031"] = GenericMatcher(r"STRESS| 1/3 Trace", col=4)
registry["M032"] = GenericMatcher(r"MD| Potential energy", col=6)
registry["M033"] = GenericMatcher(r"Dispersion energy:", col=3)
registry["M034"] = GenericMatcher(r"ISSC| all operator: CheckSum =", col=6)
registry["M035"] = GenericMatcher(r"Energy components", col=7)
registry["M036"] = GenericMatcher(r"ISSC| response: CheckSum =", col=5)
registry["M037"] = GenericMatcher(r"TDDFPT : CheckSum  =", col=5)
registry["M038"] = GenericMatcher(r"ISSC| CheckSum K =", col=5)
registry["M039"] = GenericMatcher(r"X=", col=2)
registry["M040"] = GenericMatcher(r"HELIUM| Total energy =", col=5)
registry["M041"] = GenericMatcher(r"Total charge and spin", col=9)
registry["M042"] = GenericMatcher(r"DEBUG| Sum of differences", col=5)
registry["M043"] = GenericMatcher(r"FORCES| Total core particle force", col=6)
registry["M044"] = GenericMatcher(r"FORCES| Total shell particle force", col=6)
registry["M045"] = GenericMatcher(r"PRM01", col=2)
registry["M046"] = GenericMatcher(r"Final value", col=6)
registry["M047"] = GenericMatcher(r"SUMMARY:: Number of molecule kinds found:", col=7)
registry["M048"] = GenericMatcher(r"E(Fermi):", col=3)
registry["M049"] = GenericMatcher(r"Total Multiplication", col=9)
registry["M050"] = GenericMatcher(r"Direct MP2 Canonical Energy =", col=6)
registry["M051"] = GenericMatcher(r"RESP       1", col=4)
registry["M053"] = GenericMatcher(r"Exchange-correlation energy:", col=3)
registry["M054"] = GenericMatcher(r"FORCES| Grand total force", col=6)
registry["M055"] = GenericMatcher(r"BASOPT| Total residuum value:", col=5)
registry["M056"] = GenericMatcher(r"Final value of function", col=6)
registry["M057"] = GenericMatcher(r"Emp2-RI =", col=3)
registry["M058"] = GenericMatcher(r"Final checksums", col=3)
registry["M059"] = GenericMatcher(r"GLBOPT| Lowest reported potential energy", col=7)
registry["M060"] = GenericMatcher(r"POLAR| xx,yy,zz", col=3)
registry["M061"] = GenericMatcher(r"POWELL| Final value of function", col=6)
registry["M062"] = GenericMatcher(r"Gibbs energy correction", col=6)
registry["M063"] = GenericMatcher(r"^  C", col=2)
registry["M064"] = GenericMatcher(r"^  C", col=3)
registry["M065"] = GenericMatcher(r"^  C", col=4)
registry["M066"] = GenericMatcher(r"HF Etotal", col=3)
registry["M067"] = GenericMatcher(r"Energy Level:", col=9)
registry["M068"] = GenericMatcher(r"TDDFPT|  1", col=3)
registry["M069"] = GenericMatcher(r"Log(1-CN):", col=10)
registry["M070"] = GenericMatcher(r"MD| Temperature [K]", col=4)
registry["M071"] = GenericMatcher(r"Current value of constraint", col=6)
registry["M072"] = GenericMatcher(r"FORCES| Total atomic force", col=5)
registry["M073"] = GenericMatcher(r"Diabatic electronic coupling (rotation", col=6)
registry["M074"] = GenericMatcher(r"Diabatic electronic coupling (wfn", col=7)
registry["M075"] = GenericMatcher(r"Charge transfer energy", col=6)
registry["M076"] = GenericMatcher(r"Diabatic electronic coupling (Lowdin", col=6)
registry["M077"] = GenericMatcher(r"Ground state energy", col=4)

# G0W0 HOMO-LUMO gap of molecule
registry["E_G0W0_gap"] = GenericMatcher(r"G0W0 HOMO-LUMO gap (eV)", col=5)

# G0W0 HOMO-LUMO gap of second spin channel of molecule
registry["E_G0W0_gap_beta"] = GenericMatcher(r"Beta GW HOMO-LUMO gap (eV)", col=6)

registry["IC_gap"] = GenericMatcher(r"IC HOMO-LUMO gap (eV)", col=5)

registry["M081"] = GenericMatcher(r"HOMO SCF Cycle:     4", col=9)
registry["M082"] = GenericMatcher(r"DEBUG| Sum of differences:", col=5)
registry["M083"] = GenericMatcher(r"1[   1] - 2[   1]", col=7)
registry["M084"] = GenericMatcher(r"Ionization potential of the excited atom:", col=7)
registry["M085"] = GenericMatcher(r"Total FORCE_EVAL ( SIRIUS ) energy", col=9)
registry["M086"] = GenericMatcher(r"DIPOLE : CheckSum  =", col=5)
registry["M087"] = GenericMatcher(r"POLAR : CheckSum  =", col=5)
registry["M088"] = GenericMatcher(r"XAS excitation energy (eV):", col=7)
registry["M089"] = GenericMatcher(r"Electronic density on regular grids:", col=7)
registry["M090"] = GenericMatcher(r"Final localization:", col=3)
registry["M091"] = GenericMatcher(r"Ionization potentials for XPS", col=8)
registry["M092"] = GenericMatcher(r"FCIDUMP| Checksum:", col=3)
registry["M093"] = GenericMatcher(r"SPGR| SPACE GROUP NUMBER:", col=5)
registry["M094"] = GenericMatcher(r"KS CSR write|", col=4)
registry["M095"] = GenericMatcher(r"Fermi energy:", col=3)
registry["M096"] = GenericMatcher(r"APT |   1  2", col=7)
registry["M097"] = GenericMatcher(r"r(1)", col=3)

# G0W0 bandgap of solid (old low-scaling GW implementation)
registry["E_G0W0_gap_old"] = GenericMatcher(r"GW bandgap (eV)", col=4)

# G0W0 bandgap of majority spin channel of solid (old low-scaling GW implementation)
registry["E_G0W0_alpha_gap_old"] = GenericMatcher(r"Alpha GW direct gap", col=9)

# G0W0 bandgap of minority spin channel of solid (old low-scaling GW implementation)
registry["E_G0W0_beta_gap_old"] = GenericMatcher(r"Beta GW direct gap", col=9)

# G0W0 + perturbative SOC bandgap of solid (old low-scaling GW implementation)
registry["E_G0W0_SOC_gap_old"] = GenericMatcher(r"GW+SOC bandgap (eV)", col=4)

registry["M099"] = GenericMatcher(r"Total Spread (Berry) :", col=6)
registry["M100"] = GenericMatcher(r"NVP |   1  2", col=8)
registry["M104"] = GenericMatcher(r"Ground state stabilisation:", col=4)
registry["M105"] = GenericMatcher(r"TDDFT+SOC", col=5)

# G0W0 bandgap of solid
registry["E_G0W0_direct_gap"] = GenericMatcher(r"G0W0 direct band gap", col=6)

# DFT + perturbative SOC bandgap of solid
registry["E_SCF_SOC_gap"] = GenericMatcher(r"SCF+SOC direct band gap", col=6)

# G0W0 + perturbative SOC bandgap of solid
registry["E_G0W0_SOC_gap"] = GenericMatcher(r"G0W0+SOC direct band gap", col=6)

# Lowest BSE excitation energy in Tamm-Dancoff approximation (TDA)
registry["BSE_1st_excit_ener_TDA"] = GenericMatcher(r"BSE|  1  Singlet  -TDA-", col=5)

# Lowest BSE excitation energy
registry["BSE_1st_excit_ener"] = GenericMatcher(r"BSE|  1  Singlet  -ABBA-", col=5)

# evGW0 HOMO-LUMO gap of molecule
registry["E_evGW0_gap"] = GenericMatcher(
    r"HOMO-LUMO gap in evGW0 iteration  3 (eV)", col=8
)

# evGW HOMO-LUMO gap of molecule
registry["E_evGW_gap"] = GenericMatcher(
    r"HOMO-LUMO gap in evGW iteration  3 (eV)", col=8
)

registry["M113"] = GenericMatcher(r"BSE|  1  -TDA-", col=7)
registry["M114"] = GenericMatcher(r"BSE|  1  -ABBA-", col=7)
registry["M115"] = GenericMatcher(r"MOMENTS_TRACE_RE|     0.10000000E+000", col=3)
registry["M116"] = GenericMatcher(r"MOMENTS_TRACE_IM|     0.10000000E+000", col=3)
registry["M117"] = GenericMatcher(r"MOMENTS_TRACE_RE|     0.20000000E+000", col=3)
registry["M118"] = GenericMatcher(r"POLARIZABILITY|     0.32436607E+002", col=4)
registry["M119"] = GenericMatcher(r"DEBUG:: Total Force", col=6)
registry["M120"] = GenericMatcher(r"SMEAGOL| Number of electrons:", col=5)
registry["M121"] = GenericMatcher(r"BSE|  3  -TDA-  2", col=5)
registry["M122"] = GenericMatcher(r"BSE|  3  -ABBA-  2", col=5)
registry["M123"] = GenericMatcher(r"Checksum exciton descriptors", col=4)
registry["M124"] = GenericMatcher(
    r"BSE|DEBUG| Averaged dynamical dipole polarizability at 8.2 eV:", col=9
)
registry["M125"] = GenericMatcher(
    r"BSE|DEBUG| Averaged photoabsorption cross section at 8.2 eV:", col=9
)
# Checking maximum polarizability reported by GX-AC@RTBSE
registry["RTBSE_GXAC_H2_pol"] = GenericMatcher(
    r"POLARIZABILITY_PADE|     0.30450000E+002", col=4
)

registry["M126"] = GenericMatcher(r" # Total charge ", col=5)

registry["M127"] = GenericMatcher(r"Checksum (Acoustic Sum Rule):", col=5)
# EOF
