#!/usr/bin/env python3

"""Format CP2K test lists while preserving strings, numbers, comments, and order.

Only quoted filenames and matcher tables with string, decimal, or boolean values
are supported.
"""

import argparse
import re
import sys
from pathlib import Path

STRING = r'"(?:\\[^\n]|[^"\\\n])*"' + r"|'[^'\n]*'"
DIGITS = r"[0-9](?:_?[0-9])*"
NUMBER = (
    rf"[+-]?(?:(?:0|[1-9](?:_?[0-9])*)(?:\.{DIGITS})?(?:[eE][+-]?{DIGITS})?|inf|nan)"
)
KEY = rf"(?:{STRING}|[A-Za-z0-9_-]+)"
VALUE = rf"(?:{STRING}|{NUMBER}|true|false)"
TOKEN = re.compile(rf"""\s+|{STRING}|\#[^\n]*|[=\[\]{{}},]|[^\s"'=\[\]{{}},#]+""")


# ======================================================================================
def format_gap(gap: str, indent: str, newline: bool = False) -> str:
    """Preserve comments and blank lines between array elements."""
    if not gap.strip():
        count = gap.count("\n")
        if count > 1:
            return "\n" * count + indent
        return "\n" + indent if newline else ""

    lines = gap.split("\n")
    for i, line in enumerate(lines):
        if line.strip():
            lines[i] = ("  " if i == 0 else indent) + line.strip()
        else:
            lines[i] = indent if i == len(lines) - 1 else ""
    return "\n".join(lines)


# ======================================================================================
def format_test_files(source: str, line_length: int = 120) -> str:
    if line_length <= 0:
        raise ValueError("line length must be positive")
    lexemes = list(TOKEN.finditer(source))
    if sum(len(t[0]) for t in lexemes) != len(source):
        raise ValueError("unsupported or unterminated string")
    tokens = [t for t in lexemes if not t[0].isspace() and not t[0].startswith("#")]
    pos = 0

    def peek() -> str:
        return tokens[pos][0] if pos < len(tokens) else ""

    def take(pattern: str):
        nonlocal pos
        if pos == len(tokens):
            raise ValueError("unexpected end of file")
        token = tokens[pos]
        if not re.fullmatch(pattern, token[0]):
            line = source.count("\n", 0, token.start()) + 1
            raise ValueError(f"line {line}: unexpected token {token[0]!r}")
        pos += 1
        return token

    entries = []
    while peek():
        key = take(STRING)
        start = source.rfind("\n", 0, key.start()) + 1
        if source[start : key.start()].strip():
            raise ValueError("expected one test entry per line")
        take("=")
        opening = take(r"\[")
        if "\n" in source[key.start() : opening.end()]:
            raise ValueError("expected filename, equals sign, and '[' on one line")
        matchers = []
        trailing_comma = False
        gap_start = opening.end()
        gap_prefix = ""
        while peek() != "]":
            table = take(r"\{")
            gap = gap_prefix + source[gap_start : table.start()]
            fields = []
            if peek() != "}":
                while True:
                    name = take(KEY)[0]
                    take("=")
                    value = take(VALUE)[0]
                    fields.append(f"{name}={value}")
                    if peek() == "}":
                        break
                    take(",")
            end = take(r"\}")
            if "\n" in source[table.start() : end.end()]:
                raise ValueError("multiline inline tables are not supported")
            matchers.append((gap, "{" + ", ".join(fields) + "}"))
            gap_start = end.end()
            gap_prefix = ""
            if peek() == "]":
                break
            comma = take(",")
            gap_prefix = source[gap_start : comma.start()]
            gap_start = comma.end()
            trailing_comma = peek() == "]"
        closing = take(r"\]")
        tail = gap_prefix + source[gap_start : closing.start()]
        entries.append((start, closing.end(), key[0], matchers, tail, trailing_comma))

    width = max((len(entry[2]) for entry in entries), default=0)
    width = min(width, max(0, line_length - 4))
    output = []
    end = 0
    for start, stop, key, matchers, tail, trailing_comma in entries:
        prefix = f"{key:<{width}} = ["
        single_line = prefix + ", ".join(table for gap, table in matchers)
        single_line += ("," if trailing_comma else "") + "]"
        suffix = source[stop:].split("\n", 1)[0].rstrip()
        gaps = [gap for gap, table in matchers] + [tail]
        multiline = len(single_line) > line_length or any(
            gap.strip() or gap.count("\n") > 1 for gap in gaps
        )

        def render(block: bool = False, block_indent: str = "  ") -> str:
            text = prefix
            for i, (gap, token) in enumerate(matchers + [(tail, "]")]):
                closing = i == len(matchers)
                indent = block_indent if block else " " * len(prefix)
                if closing:
                    indent = "" if block else indent[:-1]
                    if trailing_comma:
                        text += ","
                elif i:
                    text += "," if multiline or block else ", "
                if block and i == 0 and not gap.lstrip(" \t\r").startswith("\n"):
                    gap = "\n" + gap
                gap = format_gap(
                    gap,
                    indent,
                    newline=block or (bool(i) and not closing and multiline),
                )
                # Move an inline comment to its own line if it would overflow.
                if gap.strip() and not gap.startswith("\n"):
                    column = len(text.rsplit("\n", 1)[-1])
                    if column + len(gap.split("\n", 1)[0]) > line_length:
                        gap = "\n" + indent + gap.lstrip(" ")
                text += gap + token
            if suffix:
                if len(text.rsplit("\n", 1)[-1] + suffix) > line_length:
                    text += "\n" + suffix.lstrip()
                else:
                    text += suffix
            return text

        rendered = render()
        # Hanging indentation must also fit; otherwise put the array below '['.
        if any(len(line) > line_length for line in rendered.splitlines()):
            rendered = render(block=True)
            if any(len(line) > line_length for line in rendered.splitlines()):
                rendered = render(block=True, block_indent="")
        output.append(source[end:start] + rendered)
        end = stop + len(suffix)
    output.append(source[end:])
    text = "".join(output)
    return "\n".join(line.rstrip() for line in text.splitlines()) + (
        "\n" if text else ""
    )


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--line-length", type=int, default=120, help="target line length (default: 120)"
    )
    parser.add_argument("file", type=Path)
    args = parser.parse_args()
    path = args.file
    original = path.read_bytes()
    try:
        formatted = format_test_files(original.decode("utf8"), args.line_length).encode(
            "utf8"
        )
    except ValueError as err:
        sys.exit(f"{path}: {err}")
    if formatted != original:
        path.write_bytes(formatted)


if __name__ == "__main__":
    main()

# EOF
