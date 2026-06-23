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
        self, pattern: str, col: int, regex: bool = False, abs_value: bool = False
    ):
        self.pattern = pattern
        self.regex_mode = regex
        self.abs_value = abs_value
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
        for line in reversed(output.split("\n")):
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


# EOF
