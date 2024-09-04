#!/usr/bin/env python3

# author: Ole Schuett

import re
import sys
from pathlib import Path
from typing import List, Tuple, Iterator

VERBATIM_SECTIONS = ["COORD", "TRAINING_SET"]


# ======================================================================================
def separate_comment(line: str) -> Tuple[str, str]:
    m = re.match(r"([^#!]*)([#!].*)?", line)
    assert m
    body = m.group(1).strip() if m.group(1) else ""
    comment = m.group(2).strip() if m.group(2) else ""
    return body, comment


# ======================================================================================
def format_line(line: str) -> str:
    body, comment = separate_comment(line)
    tokens = body.split()
    tokens[0] = tokens[0].upper()
    padded_comment = f"  {comment}" if comment else ""
    return " ".join(tokens) + padded_comment


# ======================================================================================
def indent(lines: List[str]) -> List[str]:
    output = []
    for line in lines:
        output.append(f"  {line}")
    return output


# ======================================================================================
class Child:
    def __init__(self) -> None:
        self.sortkey = ""  # An empty sortkey means it's not sortable.

    def render(self, verbatim: bool) -> List[str]:
        return []


# ======================================================================================
class Section(Child):
    def __init__(self, preamble: List[str], line: str, children: List[Child]) -> None:
        self.preamble = preamble
        self.line = line
        self.children = children
        self.name = self.line.split()[0].upper()[1:]
        self.sortkey = f"2__{self.name}"  # Sections come after keywords.

    def render(self, verbatim: bool) -> List[str]:
        sortable = all(c.sortkey for c in self.children)
        verbatim = self.name in VERBATIM_SECTIONS
        if sortable and not verbatim:
            self.children.sort(key=lambda c: c.sortkey)

        output = self.preamble + [format_line(self.line)]
        for c in self.children:
            output += indent(c.render(verbatim))
        output.append(f"&END {self.name}")
        return output


# ======================================================================================
class Keyword(Child):
    def __init__(self, preamble: List[str], line: str) -> None:
        self.preamble = preamble
        self.line = line
        self.name = self.line.split()[0].upper()
        self.sortable = self.name[0].isalpha()  # Do not sort numeric keywords.
        self.sortkey = f"1__{self.name}" if self.sortable else ""  # Keywords come first

    def render(self, verbatim: bool) -> List[str]:
        if verbatim or not self.sortable:
            return self.preamble + [self.line]
        else:
            return self.preamble + [format_line(self.line)]


# ======================================================================================
class Preprocessor(Child):
    def __init__(self, preamble: List[str], line: str) -> None:
        self.preamble = preamble
        self.line = line
        self.sortkey = ""  # Pre-processor lines are not sortable.

    def render(self, verbatim: bool) -> List[str]:
        tokens = self.line.split(" ", 1)
        tokens[0] = tokens[0].upper()
        return self.preamble + [" ".join(tokens)]


# ======================================================================================
class Epilog(Child):
    def __init__(self, lines: List[str]) -> None:
        self.lines = lines
        self.sortkey = "9__"  # Epilogs come last.

    def render(self, verbatim: bool) -> List[str]:
        return self.lines


# ======================================================================================
def parse_children(lines_iter: Iterator[str]) -> List[Child]:
    children: List[Child] = []
    preamble: List[str] = []

    while True:
        try:
            line = next(lines_iter)
        except StopIteration:
            if preamble:
                children.append(Epilog(preamble))  # left-over preamble
            return children

        # Strip prior indentation and trailing spaces.
        line = line.strip()

        # Remove empty lines.
        if not line:
            continue

        # Split off the comment part.
        body, comment = separate_comment(line)

        # Found comment line.
        if not body:
            preamble.append(comment)

        # Found pre-processor line.
        elif body.startswith("@"):
            children.append(Preprocessor(preamble, line))
            preamble = []

        # Found section end.
        elif body.upper().startswith("&END"):
            if preamble:
                children.append(Epilog(preamble))  # left-over preamble
            return children

        # Found section begining.
        elif body.startswith("&"):
            sub_children = parse_children(lines_iter)
            children.append(Section(preamble, line, sub_children))
            preamble = []

        # Found keyword.
        else:
            children.append(Keyword(preamble, line))
            preamble = []


# ======================================================================================
def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: format_input_file.py <file>")
        sys.exit(1)

    # Parse input file.
    input_file = Path(sys.argv[1])
    lines_iter = iter(input_file.read_text(encoding="utf8").split("\n"))
    children = parse_children(lines_iter)

    # Sort the top sections, but always put &GLOBAL first and &FORCE_EVAL last.
    if not any(isinstance(c, Preprocessor) for c in children):
        sortkey_overwrite = {"2__GLOBAL": "0", "2__FORCE_EVAL": "8"}
        children.sort(key=lambda c: sortkey_overwrite.get(c.sortkey, c.sortkey))

    # Render top sections.
    output: List[str] = []
    for c in children:
        output += c.render(verbatim=False)
        if isinstance(c, Section):
            output.append("")  # Insert empty line after each top section.

    # Write output back to file
    input_file.write_text("\n".join(output), encoding="utf8")


main()

# EOF
