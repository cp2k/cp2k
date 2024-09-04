#!/usr/bin/env python3

# author: Ole Schuett

from typing import Dict, Tuple, List, Set, Optional
import xml.etree.ElementTree as ET
from pathlib import Path
import re
import sys
from datetime import datetime
from collections import defaultdict
from functools import cache


SectionPath = Tuple[str, ...]


# ======================================================================================
def main() -> None:
    if len(sys.argv) != 2:
        print("generate_input_reference.py <cp2k_input.xml>")
        sys.exit(1)

    cp2k_input_xml_fn = sys.argv[1]
    output_dir = Path(__file__).resolve().parent
    root = ET.parse(cp2k_input_xml_fn).getroot()
    build_bibliography(root, output_dir)
    build_units_reference(root, output_dir)
    build_input_reference(root, output_dir)


# ======================================================================================
def get_reference_sortkey(ref: ET.Element) -> str:
    year = int(get_text(ref.find("YEAR")))
    key = ref.attrib.get("key")
    inverse_year = 10000 - year
    return f"{inverse_year:05d}_{key}"


# ======================================================================================
def build_bibliography(root: ET.Element, output_dir: Path) -> None:
    references = root.findall("REFERENCE")
    references.sort(key=get_reference_sortkey)

    output = []
    output += ["%", "% This file was created by generate_input_reference.py", "%"]
    output += ["# Bibliography", ""]

    for r in references:
        key = r.attrib.get("key")
        authors = ", ".join([get_text(a) for a in r.findall("AUTHOR")])
        doi = get_text(r.find("DOI"))
        source = get_text(r.find("SOURCE"))
        volume = get_text(r.find("VOLUME"))
        pages = get_text(r.find("PAGES"))
        year = get_text(r.find("YEAR"))
        title = get_text(r.find("TITLE"))
        ref = source
        ref += f" **{volume}**" if volume else ""
        ref += f", {pages}"
        ref += f" ({year})" if year else ""

        output += [f"({key})=", f"## {key}", ""]
        if doi:
            output += [f"{authors}, _{title},_ [{ref}](https://doi.org/{doi}).", ""]
        else:
            output += [f"{authors}, _{title},_ {ref}.", ""]

    write_file(output_dir / "bibliography.md", "\n".join(output))


# ======================================================================================
def build_units_reference(root: ET.Element, output_dir: Path) -> None:
    output = []
    output += ["%", "% This file was created by generate_input_reference.py", "%"]
    output += ["# Units", "", "Units of measurement available in CP2K's input.", ""]

    for unit_kind in root.findall("UNIT_KIND"):
        kind = unit_kind.attrib.get("name")
        assert kind
        article = "an" if kind[0] in "aeiou" else "a"
        output += [f"## {kind.capitalize()}", ""]
        output += [f"Possible units of measurement for {article} {kind}."]
        output += [f"The `[{kind}]` entry acts like a dummy flag (assumes the"]
        output += [f"unit of measurement of {kind} is in internal units),"]
        output += [f"useful for dimensional analysis.", ""]
        for unit in unit_kind.findall("UNIT"):
            output += [f"- `{get_text(unit)}`"]
        output += [""]

    write_file(output_dir / "units.md", "\n".join(output))


# ======================================================================================
def build_input_reference(root: ET.Element, output_dir: Path) -> None:
    num_files_written = process_section(root, ("CP2K_INPUT",), False, output_dir)

    # Build landing page.
    cp2k_version = get_text(root.find("CP2K_VERSION"))
    compile_revision = get_text(root.find("COMPILE_REVISION"))
    # cp2k_year = get_text(root.find("CP2K_YEAR"))
    # compile_date = get_text(root.find("COMPILE_DATE"))

    output = []
    output += ["%", "% This file was created by generate_input_reference.py", "%"]
    output += ["(CP2K_INPUT)="]
    output += ["# Input Reference", ""]

    assert compile_revision.startswith("git:")
    github_url = f"https://github.com/cp2k/cp2k/tree/{compile_revision[4:]}"
    output += [f"Based on {cp2k_version} ([{compile_revision}]({github_url}))", ""]

    output += ["```{toctree}"]
    output += [":maxdepth: 1"]
    output += [":titlesonly:"]
    output += [":caption: Top-level sections"]
    output += [":glob:", ""]
    output += ["CP2K_INPUT/*", ""]

    write_file(output_dir / "CP2K_INPUT.md", "\n".join(output))  # Replace generic file.
    print(f"Generated markdown files for {num_files_written} input sections.")


# ======================================================================================
def process_section(
    section: ET.Element,
    section_path: SectionPath,
    has_name_collision: bool,
    output_dir: Path,
) -> int:
    # Render section header.
    output, section_name, section_xref = render_section_header(
        section, section_path, has_name_collision
    )

    # Render TOC for subsections
    subsections = sorted(section.findall("SECTION"), key=get_name)
    if subsections:
        output += ["```{toctree}"]
        output += [":maxdepth: 1"]
        output += [":titlesonly:"]
        output += [":caption: Subsections"]
        for subsection in subsections:
            output += [f"{section_name}/{get_name(subsection)}"]
        output += ["```", ""]

    # Collect and sort keywords
    keywords = (
        section.findall("SECTION_PARAMETERS")
        + section.findall("DEFAULT_KEYWORD")
        + sorted(section.findall("KEYWORD"), key=get_name)
    )

    # Filter out removed keywords.
    keywords = [k for k in keywords if k.attrib.get("removed", "no") == "no"]

    # Render keywords
    if keywords:
        # Render TOC for keywords
        output += ["## Keywords", ""]
        for keyword in keywords:
            keyword_name = get_name(keyword)
            keyword_xref = f"{section_xref}.{sanitize_name(keyword_name)}"
            em = "**" if lookup_mentions(keyword_xref) else ""  # emphasize if mentioned
            output += [f"* {em}[{escape_markdown(keyword_name)}](#{keyword_xref}){em}"]
        output += [""]
        # Render keywords
        output += ["## Keyword descriptions", ""]
        for keyword in keywords:
            output += render_keyword(keyword, section_xref)

    # Write output
    section_dir = output_dir / "/".join(section_path[:-1])
    section_dir.mkdir(exist_ok=True)
    write_file(section_dir / f"{section_name}.md", "\n".join(output))
    num_files_written = 1

    # Process subsections
    keyword_names = {get_name(keyword) for keyword in keywords}
    for subsection in subsections:
        subsection_name = get_name(subsection)
        subsection_path = (*section_path, subsection_name)
        has_collision = subsection_name in keyword_names
        if subsection_name == "XC_FUNCTIONAL":
            n = process_xc_functional_section(subsection, subsection_path, output_dir)
        else:
            n = process_section(subsection, subsection_path, has_collision, output_dir)
        num_files_written += n

    return num_files_written


# ======================================================================================
def process_xc_functional_section(
    section: ET.Element,
    section_path: SectionPath,
    output_dir: Path,
) -> int:
    # Render section header and keywords.
    assert section_path[-1] == "XC_FUNCTIONAL"
    output, section_name, section_xref = render_section_header(section, section_path)
    for keyword in section.findall("SECTION_PARAMETERS") + section.findall("KEYWORD"):
        output += render_keyword(keyword, section_xref)

    # Render note.
    output += ["```{note}"]
    output += ["Thanks to [Libxc](https://libxc.gitlab.io/) there are 600+"]
    output += ["functionals available. For ease of browsing their documentation has"]
    output += ["been inlined into this page. Each of the functionals has a"]
    output += ["corresponding subsection that, if present, enables the functional."]
    output += ["```", ""]

    # Index functionals by prefix.
    subsections = sorted(section.findall("SECTION"), key=get_name)
    prefixes = ["", "LDA_X_", "LDA_C_", "LDA_XC_", "LDA_K_", "HYB_LDA_XC_", "GGA_X_"]
    prefixes += ["GGA_C_", "GGA_XC_", "GGA_K_", "HYB_GGA_X_", "HYB_GGA_XC_", "MGGA_X_"]
    prefixes += ["MGGA_C_", "MGGA_XC_", "MGGA_K_", "HYB_MGGA_X_", "HYB_MGGA_XC_"]
    subsections_by_prefix: Dict[str, List[str]] = {prefix: [] for prefix in prefixes}
    for subsection in subsections:
        subsection_name = get_name(subsection)
        if subsection_name == "LDA_X":
            continue  # Special case where libxc deviates from its naming convention.
        for prefix in reversed(prefixes):  # reverse order because "" always matches
            if subsection_name.startswith(prefix):
                subsections_by_prefix[prefix].append(subsection_name)
                break

    # Render list of functionals for each prefix.
    for prefix, section_names in subsections_by_prefix.items():
        category = f"Libxc {prefix[:-1]}" if prefix else "Built-in"
        items = [f"[{s[len(prefix):]}](#{section_xref}.{s})" for s in section_names]
        if prefix == "LDA_X_":
            items.insert(0, f"[LDA_X](#{section_xref}.LDA_X)")  # special case
        output += [f"## {category} Functionals", "", ",\n".join(items), ""]

    # Render inline subsections
    for subsection in subsections:
        subsection_name = get_name(subsection)
        output += [f"({section_xref}.{subsection_name})=", f"## {subsection_name}", ""]
        output += [escape_markdown(get_text(subsection.find("DESCRIPTION"))), ""]
        for keyword in sorted(subsection.findall("KEYWORD"), key=get_name):
            output += render_keyword(keyword, section_xref=None, github=False)

    # Write output
    section_dir = output_dir / "/".join(section_path[:-1])
    write_file(section_dir / "XC_FUNCTIONAL.md", "\n".join(output))
    return 1


# ======================================================================================
def render_section_header(
    section: ET.Element,
    section_path: SectionPath,
    has_name_collision: bool = False,
) -> Tuple[List[str], str, str]:
    # Collect information from section fields.
    repeats = "repeats" in section.attrib and section.attrib["repeats"] == "yes"
    description = get_text(section.find("DESCRIPTION"))
    deprecation_notice = get_text(section.find("DEPRECATION_NOTICE"))
    location = get_text(section.find("LOCATION"))
    section_name = section_path[-1]  # section.find("NAME") doesn't work for root
    section_xref = ".".join(section_path)  # used for cross-referencing
    references = [get_name(ref) for ref in section.findall("REFERENCE")]

    # Render header.
    output = []
    output += ["%", "% This file was created by generate_input_reference.py", "%"]
    # There are a few collisions between cross references for sections and keywords,
    # for example CP2K_INPUT.FORCE_EVAL.SUBSYS.KIND.POTENTIAL
    collision_resolution_suffix = "_SECTION" if has_name_collision else ""
    output += [f"({section_xref}{collision_resolution_suffix})="]
    output += [f"# {section_name}", ""]
    if deprecation_notice:
        output += ["```{warning}"]
        output += ["This section is deprecated and may be removed in a future version."]
        output += ["", escape_markdown(deprecation_notice), "", "```", ""]
    if repeats:
        output += ["**Section can be repeated.**", ""]
    if references:
        citations = ", ".join([f"{{ref}}`{r}`" for r in references])
        output += [f"**References:** {citations}", ""]
    output += [escape_markdown(description), github_link(location), ""]
    return output, section_name, section_xref


# ======================================================================================
def render_keyword(
    keyword: ET.Element,
    section_xref: Optional[str],
    github: bool = True,
) -> List[str]:
    # Find keyword names.
    keyword_names: List[str]
    if keyword.tag == "SECTION_PARAMETERS":
        keyword_names = ["SECTION_PARAMETERS"]
    elif keyword.tag == "DEFAULT_KEYWORD":
        keyword_names = ["DEFAULT_KEYWORD"]
    else:
        keyword_names = [get_text(name) for name in keyword.findall("NAME")]
    assert keyword_names
    canonical_name = sanitize_name(keyword_names[0])
    keyword_xref = f"{section_xref}.{canonical_name}" if section_xref else None
    mentions = lookup_mentions(keyword_xref)

    # Find more keyword fields.
    default_value = get_text(keyword.find("DEFAULT_VALUE"))
    default_unit = get_text(keyword.find("DEFAULT_UNIT"))
    usage = get_text(keyword.find("USAGE"))
    description = get_text(keyword.find("DESCRIPTION"))
    deprecation_notice = get_text(keyword.find("DEPRECATION_NOTICE"))
    location = get_text(keyword.find("LOCATION"))
    lone_keyword_value = get_text(keyword.find("LONE_KEYWORD_VALUE"))

    # Find keyword data type.
    data_type_element = keyword.find("DATA_TYPE")
    assert data_type_element is not None
    data_type = data_type_element.attrib["kind"]
    if data_type == "word":
        data_type = "string"
    if data_type == "keyword":
        data_type = "enum"

    # Need to distiguish between multiple values (n_var) and repeating keyword.
    repeats = keyword.attrib["repeats"] == "yes"
    n_var = int(get_text(data_type_element.find("N_VAR")))

    # Find keyword references.
    references = [get_name(ref) for ref in keyword.findall("REFERENCE")]

    output: List[str] = []

    # Include HTML anchors to preserve old links.
    output += [f"<a id='list_{keyword_names[0]}'></a>"]
    output += [f"<a id='desc_{keyword_names[0]}'></a>"]
    output += [f"<a id='{keyword_names[0]}'></a>", ""]

    # Use Sphinx's py:data directive to document keywords.
    output += [f"```{{py:data}}  {canonical_name}"]
    n_var_brackets = f"[{n_var}]" if n_var > 1 else ""
    if section_xref:
        output += [f":module: {section_xref}"]
    else:
        output += [":noindex:"]
    output += [f":type: '{data_type}{n_var_brackets}'"]
    if default_value or default_unit:
        default_unit_bracketed = f"[{default_unit}]" if default_unit else ""
        output += [f":value: '{default_value} {default_unit_bracketed}'"]
    output += [""]
    if repeats:
        output += ["**Keyword can be repeated.**", ""]
    if len(keyword_names) > 1:
        aliases = " ,".join(keyword_names[1:])
        output += [f"**Aliases:** {escape_markdown(aliases)}", ""]
    if lone_keyword_value:
        output += [f"**Lone keyword:** `{escape_markdown(lone_keyword_value)}`", ""]
    if usage:
        output += [f"**Usage:** _{escape_markdown(usage)}_", ""]
    if data_type == "enum":
        output += ["**Valid values:**"]
        for item in keyword.findall("DATA_TYPE/ENUMERATION/ITEM"):
            item_description = get_text(item.find("DESCRIPTION"))
            output += [f"* `{get_name(item)}`"]
            output += [indent(escape_markdown(item_description))]
        output += [""]
    if references:
        citations = ", ".join([f"{{ref}}`{r}`" for r in references])
        output += [f"**References:** {citations}", ""]
    if mentions:
        mentions_list = ", ".join([f"⭐[](project:{m})" for m in mentions])
        output += [f"**Mentions:** {mentions_list}", ""]
    output += [escape_markdown(description)]
    if github:
        output += [github_link(location)]
    output += ["", "```", ""]  # Close py:data directive.

    if deprecation_notice:
        output += ["```{warning}", f"The keyword [{canonical_name}](#{keyword_xref})"]
        output += ["is deprecated and may be removed in a future version.", ""]
        output += [escape_markdown(deprecation_notice), "", "```", ""]

    return output


# ======================================================================================
@cache
def find_all_mentions() -> Dict[str, Set[Path]]:
    root_dir = Path(__file__).resolve().parent
    mentions = defaultdict(set)
    for fn in (root_dir / "methods").glob("**/*.md"):
        for xref in re.findall(r"\(#(CP2K_INPUT\..*)\)", fn.read_text()):
            mentions[xref].add(fn.relative_to(root_dir))
    return mentions


# ======================================================================================
def lookup_mentions(xref: Optional[str]) -> List[str]:
    if not xref:
        return []
    mentions = find_all_mentions()
    n = xref.count(".") - 1
    return [("../" * n) + str(path) for path in sorted(mentions[xref])]


# ======================================================================================
def get_name(element: ET.Element) -> str:
    return get_text(element.find("NAME"))


# ======================================================================================
def get_text(element: Optional[ET.Element]) -> str:
    if element is not None:
        if element.text is not None:
            return element.text
    return ""


# ======================================================================================
def sanitize_name(name: str) -> str:
    name = name.replace("-", "_")
    name = name.replace("+", "_")
    name = name.replace("[", "_")
    name = name.replace("]", "")
    return name


# ======================================================================================
def escape_markdown(text: str) -> str:
    # Code blocks without a language get mistaken for the end of the py:data directive.
    text = text.replace("\n\n```\n", "\n\n```none\n")

    # Underscores are very common in our docs. Luckily asterisks also work for emphasis.
    text = text.replace(r"__", r"\_\_")

    # Headings mess up the page structure. Please use paragraphs and bold text instead.
    text = text.replace(r"#", r"\#")

    return text


# ======================================================================================
def indent(text: str) -> str:
    return "\n".join(f"  {line}" for line in text.split("\n"))


# ======================================================================================
def github_link(location: str) -> str:
    if not location:
        return ""
    location_url = location.replace(":", "#L")
    github_url = f"https://github.com/cp2k/cp2k/blob/master/src/{location_url}"
    return f"<small>[[Edit on GitHub]({github_url})]</small>"


# ======================================================================================
def write_file(filename: Path, content: str) -> None:
    old_content = filename.read_text() if filename.exists() else None
    if old_content != content:
        filename.write_text(content)
        print(f"Wrote {filename}")


# ======================================================================================
main()

# EOF
