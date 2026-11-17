import math
import re
import sys
from dataclasses import dataclass
from typing import Any, Callable, Optional

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
def _compare_value(value: float, tol: float, ref: float) -> MatchResult:
    diff = value - ref
    error = abs(diff / ref if ref != 0.0 else diff)
    if error > tol:
        return MatchResult(
            "WRONG RESULT",
            f"Difference too large: {error:.2e} > {tol}, value: {value}.\n",
            value,
        )
    return MatchResult("OK", error=None, value=value)


# ======================================================================================
class ComputedValueMatcher(Matcher):
    """Compare a value derived from the output with a reference."""

    def __init__(self, compute: Callable[[str], float]):
        self.compute = compute

    def run(self, output: str, **kwargs: Any) -> MatchResult:
        try:
            value = self.compute(output)
        except ValueError as error:
            return MatchResult("WRONG RESULT", f"{error}\n", None)
        if not math.isfinite(value):
            return MatchResult(
                "WRONG RESULT", f"Result is not finite: '{value}'.\n", value
            )
        return _compare_value(value, kwargs["tol"], kwargs.get("ref", 0.0))


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

        return _compare_value(value, tol, ref)


# ======================================================================================
def restricted_fcidump_energy_offset(output: str) -> float:
    """Reconstruct the restricted reference determinant energy independently."""
    headers = list(
        re.finditer(
            r"^&FCI\s+NORB=\s*(\d+),NELEC=\s*(\d+),MS2=\s*(\d+),[^/]*?/",
            output,
            re.M,
        )
    )
    if not headers or re.search(r"UHF=\s*1", headers[-1].group()):
        raise ValueError("No restricted FCIDUMP found.")
    norb, nelec, ms2 = map(int, headers[-1].groups())
    na, nb = (nelec + ms2) // 2, (nelec - ms2) // 2
    if not 0 <= nb <= na <= norb or na + nb != nelec:
        raise ValueError("Invalid FCIDUMP occupations.")

    # CP2K writes one record per integral permutation class.
    energy = 0.0
    for line in output[headers[-1].end() :].strip().splitlines():
        value, *indices = line.split()
        i, j, k, l = map(int, indices)
        integral = float(value.replace("D", "E"))
        if i == j == k == l == 0:
            energy += integral
            break
        if k == l == 0 and i == j:
            energy += ((i <= na) + (i <= nb)) * integral
        elif i == j and k == l:
            if i == k:
                energy += (i <= nb) * integral
            else:
                energy += ((i <= na) + (i <= nb)) * ((k <= na) + (k <= nb)) * integral
        elif (i, j) == (k, l):
            energy -= ((max(i, j) <= na) + (max(i, j) <= nb)) * integral
    else:
        raise ValueError("Missing FCIDUMP constant.")

    scf = GenericMatcher("ENERGY| Total FORCE_EVAL", col=9).run(
        output, tol=math.inf, ref=0.0
    )
    if scf.status != "OK" or scf.value is None:
        raise ValueError((scf.error or "No SCF energy found.").rstrip())
    return energy - scf.value


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


# EOF
