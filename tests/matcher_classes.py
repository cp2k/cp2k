import math
import re
import sys
from dataclasses import dataclass
from typing import Any, Callable, Optional, Union

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
    """Compare a column or a value extracted from the full output with a reference."""

    def __init__(
        self,
        pattern: Union[str, Callable[[str], float]],
        col: Optional[int] = None,
        regex: bool = False,
        abs_value: bool = False,
        first: bool = False,
    ):
        self.abs_value = abs_value
        self.extract: Callable[[str], float]
        if callable(pattern):
            if col is not None or regex or first:
                raise ValueError("Column and pattern options require a string pattern.")
            self.extract = pattern
            return

        if col is None:
            raise ValueError("A column is required for a string pattern.")
        self.pattern = pattern
        self.regex_mode = regex
        self.first = first
        if not regex:
            for c in r"[]()|+*?":
                pattern = pattern.replace(c, f"\\{c}")
        self.regex = re.compile(pattern)
        self.col = col
        self.extract = self._extract_column

    def _extract_column(self, output: str) -> float:
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
            raise ValueError(f"Result not found: '{self.pattern}'.")

        try:
            return float(value_str.replace("D", "E"))
        except ValueError:
            raise ValueError(
                f"Could not parse result as float: '{value_str}'."
            ) from None

    def run(self, output: str, **kwargs: Any) -> MatchResult:
        tol, ref = kwargs["tol"], kwargs["ref"]
        assert isinstance(tol, float) or isinstance(ref, int)
        assert isinstance(ref, float) or isinstance(ref, int)
        try:
            value = self.extract(output)
            if self.abs_value:
                value = abs(value)
        except ValueError as exc:
            return MatchResult("WRONG RESULT", f"{exc}\n", value=None)

        if not math.isfinite(value):
            error = f"Result is not finite: '{value}'.\n"
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


# EOF
