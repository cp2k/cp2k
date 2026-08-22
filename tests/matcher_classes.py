import math
import re
import sys
from dataclasses import dataclass
from typing import Any, List, Optional

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
class CP2KMatrixChecksumMatcher(Matcher):
    def run(self, output: str, **kwargs: Any) -> MatchResult:
        tol, ref = kwargs["tol"], kwargs["ref"]
        title, dimension = kwargs["title"], kwargs.get("dimension", 16)
        if (
            not isinstance(tol, (float, int))
            or not isinstance(ref, list)
            or not ref
            or not all(isinstance(value, (float, int)) for value in ref)
            or not isinstance(title, str)
            or not isinstance(dimension, int)
            or dimension < 1
        ):
            return MatchResult(
                "N/A", "Invalid CP2K matrix matcher arguments.\n", value=None
            )

        sections = re.split(rf"(?m)^\s+{re.escape(title)}\s*$", output)[1:]
        if len(sections) != len(ref):
            error = f"Expected {len(ref)} {title} matrices, found {len(sections)}.\n"
            return MatchResult("WRONG RESULT", error, value=None)

        values = []
        for section in sections:
            rows: List[List[float]] = [[] for _ in range(dimension)]
            for line in section.splitlines():
                fields = line.split()
                if (
                    len(fields) < 5
                    or not fields[0].isdigit()
                    or not fields[1].isdigit()
                ):
                    continue
                try:
                    row = int(fields[0])
                    row_values = [
                        float(value.replace("D", "E")) for value in fields[4:]
                    ]
                except ValueError:
                    continue
                if 1 <= row <= dimension and row_values:
                    rows[row - 1].extend(row_values)
                if all(len(row_values) == dimension for row_values in rows):
                    break

            if not all(len(row_values) == dimension for row_values in rows):
                error = f"Could not parse a complete {dimension}x{dimension} {title} matrix.\n"
                return MatchResult("WRONG RESULT", error, value=None)
            values.append(
                sum(
                    (row + 1) * (col + 1) * value
                    for row, row_values in enumerate(rows)
                    for col, value in enumerate(row_values)
                )
            )

        if not all(math.isfinite(value) for value in values):
            return MatchResult("WRONG RESULT", "Checksum is not finite.\n", value=None)

        errors = [
            abs(value - reference) / abs(reference) if reference != 0.0 else abs(value)
            for value, reference in zip(values, ref)
        ]
        if any(error > tol for error in errors):
            checksums = ", ".join(
                f"matrix {i + 1}: {value:.14g}" for i, value in enumerate(values)
            )
            return MatchResult(
                "WRONG RESULT",
                f"{title} checksum mismatch (relative errors: {errors}); values: {checksums}.\n",
                value=None,
            )

        return MatchResult("OK", error=None, value=None)


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


# EOF
