"""Small, schema-independent serializer for native CP2K input dictionaries."""

# SPDX-License-Identifier: GPL-2.0-or-later

from collections.abc import Mapping
import math
from numbers import Real
import re


def _value(value):
    if isinstance(value, bool):
        return ".TRUE." if value else ".FALSE."
    if isinstance(value, str):
        if any(c in value for c in "\n\r\0"):
            raise ValueError("Input values must be single lines without NUL characters")
        return value
    if isinstance(value, Real) and math.isfinite(value):
        return str(value)
    raise TypeError(f"Unsupported CP2K input value: {value!r}")


def input_to_string(tree):
    """Serialize a mapping to CP2K input (without selecting physical defaults).

    A mapping is a section; a list of mappings repeats a section. ``_`` gives
    section parameters (e.g. ``KIND: {_: H, ...}``), ``_lines`` holds unparsed
    single lines (e.g. coordinates). A list/tuple of scalars gives keyword
    arguments; a list of lists repeats a keyword. ``None`` is a bare keyword.
    Keys are case insensitive. Semantic validation remains CP2K's job.
    """
    if not isinstance(tree, Mapping):
        raise TypeError("CP2K input must be a mapping")

    def contents(section, depth):
        lines = []
        indent = "  " * depth
        seen = set()
        for key, value in section.items():
            if not isinstance(key, str):
                raise TypeError("Input keys must be strings")
            name = key.upper()
            if name in seen:
                raise ValueError(f"Duplicate input key: {key}")
            seen.add(name)
            if name == "_":
                if depth == 0:
                    raise ValueError("Section parameters are not valid at the root")
                continue
            if name == "_LINES":
                if not isinstance(value, (list, tuple)):
                    raise TypeError("_lines must be a list of single lines")
                for line in value:
                    text = _value(line)
                    if text.lstrip().startswith(("&", "@")):
                        raise ValueError("Use mappings for sections, not _lines")
                    lines.append(indent + text)
                continue
            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_+\-]*", key):
                raise ValueError(f"Invalid CP2K input key: {key!r}")
            sections = [value] if isinstance(value, Mapping) else value
            if (
                isinstance(sections, (list, tuple))
                and sections
                and all(isinstance(item, Mapping) for item in sections)
            ):
                for child in sections:
                    parameter = child.get("_")
                    suffix = " " + _value(parameter) if parameter is not None else ""
                    lines.append(f"{indent}&{name}{suffix}")
                    lines.extend(contents(child, depth + 1))
                    lines.append(f"{indent}&END {name}")
            elif value is None:
                lines.append(indent + name)
            elif isinstance(value, (list, tuple)):
                if not value:
                    raise ValueError(f"Empty value list for {name}")
                entries = value if isinstance(value[0], (list, tuple)) else [value]
                for entry in entries:
                    if not isinstance(entry, (list, tuple)) or not entry:
                        raise ValueError(f"Invalid repeated keyword {name}")
                    lines.append(indent + name + " " + " ".join(map(_value, entry)))
            else:
                lines.append(indent + name + " " + _value(value))
        return lines

    return "\n".join(contents(tree, 0)) + "\n"
