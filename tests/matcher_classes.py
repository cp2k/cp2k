import math
import re
import sys
from dataclasses import dataclass
from typing import Any, Optional

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
class Matcher(Protocol):
    def run(self, output: str, **kwargs: Any) -> MatchResult: ...


# ======================================================================================
class GenericMatcher(Matcher):
    def __init__(
        self,
        pattern: str,
        col: int,
        regex: bool = False,
        abs_value: bool = False,
        first: bool = False,
    ):
        self.pattern = pattern
        self.regex_mode = regex
        self.abs_value = abs_value
        self.first = first
        if not regex:
            for c in r"[]()|+*?":
                pattern = pattern.replace(c, f"\\{c}")
        self.regex = re.compile(pattern)
        self.col = col

    def run(self, output: str, **kwargs: Any) -> MatchResult:
        tol, ref = kwargs["tol"], kwargs["ref"]
        assert isinstance(tol, float) or isinstance(ref, int)
        assert isinstance(ref, float) or isinstance(ref, int)
        # grep result
        lines = output.split("\n")
        for line in lines if self.first else reversed(lines):
            match = self.regex.search(line)
            if match:
                if self.regex_mode and match.groups():
                    value_str = match.group(1)
                else:
                    value_str = line.split()[self.col - 1]
                break
        else:
            error = f"Result not found: '{self.pattern}'.\n"
            return MatchResult("WRONG RESULT", error, value=None)

        # parse result
        try:
            value = float(value_str.replace("D", "E"))
            if self.abs_value:
                value = abs(value)
        except:
            error = f"Could not parse result as float: '{value_str}'.\n"
            return MatchResult("WRONG RESULT", error, value=None)

        if not math.isfinite(value):
            error = f"Result is not finite: '{value_str}'.\n"
            return MatchResult("WRONG RESULT", error, value)

        # compare result to reference
        diff = value - ref
        rel_error = abs(diff / ref if ref != 0.0 else diff)
        if rel_error > tol:
            error = f"Difference too large: {rel_error:.2e} > {tol}, value: {value}.\n"
            return MatchResult("WRONG RESULT", error, value)

        return MatchResult("OK", error=None, value=value)  # passed


# ======================================================================================
class TextPresenceMatcher(Matcher):
    def __init__(self, text: str):
        self.text = text

    def run(self, output: str, **kwargs: Any) -> MatchResult:
        if self.text not in output:
            return MatchResult(
                "WRONG RESULT", f"Text not found: '{self.text}'.\n", value=None
            )
        return MatchResult("OK", error=None, value=None)


# ======================================================================================
class TextAbsenceMatcher(Matcher):
    def __init__(self, text: Optional[str] = None):
        self.text = text

    def run(self, output: str, **kwargs: Any) -> MatchResult:
        text = kwargs.get("text", self.text)
        assert isinstance(text, str) and text
        count = output.count(text)
        if count > 0:
            first_context = ""
            for line_nr, line in enumerate(output.splitlines(), start=1):
                if text in line:
                    first_context = (
                        f"First occurrence at output line {line_nr}: {line}\n"
                    )
                    break
            return MatchResult(
                "WRONG RESULT",
                f"Forbidden text found {count} time(s): '{text}'.\n{first_context}",
                value=None,
            )
        return MatchResult("OK", error=None, value=None)


# ======================================================================================
# TODO remove
class GWOrbitalSelfEnergyMatcher(Matcher):
    """Check the shared real-frequency orbital file and one reference value."""

    def run(self, output: str, **kwargs: Any) -> MatchResult:
        method = kwargs["method"]
        lines = output.splitlines()
        if f"# Method: {method}" not in lines:
            return MatchResult(
                "WRONG RESULT", "Incorrect GW self-energy method.\n", None
            )
        header = "# Real_frequency_eV Re_Sigma_c_eV Im_Sigma_c_eV A_n_1_per_eV"
        if header not in lines:
            return MatchResult(
                "WRONG RESULT", "Missing GW self-energy columns.\n", None
            )
        try:
            rows = [
                tuple(map(float, line.split()))
                for line in lines
                if line and not line.startswith("#")
            ]
        except ValueError:
            return MatchResult("WRONG RESULT", "Invalid GW self-energy row.\n", None)
        if len(rows) != kwargs["points"] or any(
            len(row) != 4 or not all(map(math.isfinite, row)) for row in rows
        ):
            return MatchResult(
                "WRONG RESULT", "Nonfinite or incomplete GW self-energy rows.\n", None
            )
        if any(len(line) > 80 for line in lines):
            return MatchResult(
                "WRONG RESULT", "GW self-energy line exceeds 80 columns.\n", None
            )
        if any(row[3] < -1e-12 for row in rows):
            return MatchResult(
                "WRONG RESULT", "Negative orbital spectral function.\n", None
            )
        if method == "EXACT_RPA_DIAG" and any(row[2] > 1e-12 for row in rows):
            return MatchResult(
                "WRONG RESULT", "Positive retarded RPA imaginary self-energy.\n", None
            )
        if any(rows[i + 1][0] <= rows[i][0] for i in range(len(rows) - 1)):
            return MatchResult(
                "WRONG RESULT", "Nonmonotonic real-frequency grid.\n", None
            )
        metadata = dict(
            line[2:].split(": ", 1)
            for line in lines
            if line.startswith("# ") and ": " in line
        )
        epsilon = float(metadata["SCF_energy_eV"])
        exchange = float(metadata["Sigma_x_minus_vxc_eV"])
        eta = float(metadata["Green_broadening_eV"])
        for omega, real, imag, spectral in rows:
            gamma = eta - imag if method == "EXACT_RPA_DIAG" else eta + abs(imag)
            expected = gamma / (
                math.pi * ((omega - epsilon - exchange - real) ** 2 + gamma**2)
            )
            if not math.isclose(spectral, expected, rel_tol=2e-8, abs_tol=1e-12):
                return MatchResult(
                    "WRONG RESULT",
                    "Spectral function disagrees with documented Dyson equation.\n",
                    None,
                )
        value = min(rows, key=lambda row: abs(row[0]))[1]
        ref, tol = kwargs["ref"], kwargs["tol"]
        if abs(value - ref) > tol * max(abs(ref), 1.0):
            return MatchResult(
                "WRONG RESULT", f"Real self-energy at zero differs: {value}.\n", value
            )
        return MatchResult("OK", None, value)


# EOF
