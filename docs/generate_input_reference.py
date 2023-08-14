#!/usr/bin/env python3

# author: Ole Schuett

from typing import Tuple, List, Optional
import lxml.etree as ET
import lxml
from pathlib import Path
import re
import sys

SectionPath = Tuple[str, ...]


# ======================================================================================
def main() -> None:
    if len(sys.argv) != 3:
        print("generate_input_reference.py <cp2k_input.xml> <references.html>")
        sys.exit(1)

    cp2k_input_xml_fn, references_html_fn = sys.argv[1:]
    output_dir = Path(__file__).resolve().parent

    build_bibliography(references_html_fn, output_dir)
    build_input_reference(cp2k_input_xml_fn, output_dir)


# ======================================================================================
def build_bibliography(references_html_fn: str, output_dir: Path) -> None:
    content = Path(references_html_fn).read_text()
    entries = re.findall("<TR>.*?</TR>", content, re.DOTALL)

    output = []
    output += ["%", "% This file was created by generate_input_reference.py", "%"]
    output += ["# Bibliography", ""]

    for entry in entries:
        pattern = (
            '<TR><TD>\[(.*?)\]</TD><TD>\n <A NAME="reference_\d+">(.*?)</A><br>'
            "(.*?)</TD></TR>"
        )
        parts = re.search(pattern, entry, re.DOTALL)
        assert parts
        key = parts.group(1)
        title = parts.group(3).strip()
        if "<br>" in parts.group(2):
            m = re.match("(.*?)<br>(.*)", parts.group(2), re.DOTALL)
            assert m
            authors, mix = m.groups()
        else:
            authors, mix = "", parts.group(2)

        if "https://doi.org" in mix:
            m = re.match('\s*<A HREF="(.*?)">(.*?)</A>', mix, re.DOTALL)
            assert m
            doi, ref = m.groups()
        else:
            doi, ref = "", mix.strip()

        output += [f"({key})=", f"## {key}", ""]
        if doi:
            output += [f"{authors} **{title}** _[{ref}]({doi})_", ""]
        else:
            output += [f"{authors} **{title}** _{ref}_", ""]

    # Write output
    filename = output_dir / "bibliography.md"
    filename.write_text("\n".join(output))
    print(f"Wrote {filename}")


# ======================================================================================
def build_input_reference(cp2k_input_xml_fn: str, output_dir: Path) -> None:
    tree = ET.parse(cp2k_input_xml_fn)
    root = tree.getroot()
    num_files_written = process_section(root, ("CP2K_INPUT",), output_dir)

    # Build landing page.
    cp2k_version = get_text(root.find("CP2K_VERSION"))
    compile_revision = get_text(root.find("COMPILE_REVISION"))
    # cp2k_year = get_text(root.find("CP2K_YEAR"))
    # compile_date = get_text(root.find("COMPILE_DATE"))

    output = []
    output += ["%", "% This file was created by generate_input_reference.py", "%"]
    output += ["# Input reference", ""]

    assert compile_revision.startswith("git:")
    github_url = f"https://github.com/cp2k/cp2k/tree/{compile_revision[4:]}"
    output += [f"Based on {cp2k_version} ([{compile_revision}]({github_url}))", ""]

    output += ["```{toctree}"]
    output += [":maxdepth: 1"]
    output += [":titlesonly:"]
    output += [":caption: Top-level sections"]
    output += [":glob:", ""]
    output += ["CP2K_INPUT/*", ""]

    # Write output
    filename = output_dir / "CP2K_INPUT.md"  # Overwrite generic file.
    filename.write_text("\n".join(output))
    print(f"Wrote markdown files for {num_files_written} input sections.")


# ======================================================================================
def process_section(
    section: lxml.etree._Element, section_path: SectionPath, output_dir: Path
) -> int:
    # Find more section fields.
    repeats = "repeats" in section.attrib and section.attrib["repeats"] == "yes"
    description = get_text(section.find("DESCRIPTION"))
    location = get_text(section.find("LOCATION"))
    section_name = section_path[-1]  # section.find("NAME") doesn't work for root
    section_xref = ".".join(section_path[1:])  # used for cross-referencing

    # Find section references.
    references = [get_name(ref) for ref in section.findall("REFERENCE")]

    output = []
    output += ["%", "% This file was created by generate_input_reference.py", "%"]
    output += [f"# {section_name}", ""]
    if repeats:
        output += ["**Section can be repeated.**", ""]
    if references:
        citations = ", ".join([f"{{ref}}`{r}`" for r in references])
        output += [f"**References:** {citations}", ""]
    output += [f"{escape_markdown(description)} {github_link(location)}", ""]

    # Render TOC for subsections
    if section.findall("SECTION"):
        output += ["```{toctree}"]
        output += [":maxdepth: 1"]
        output += [":titlesonly:"]
        output += [":caption: Subsections"]
        output += [":glob:", ""]
        output += [f"{section_name}/*"]  # TODO maybe list subsection explicitly.
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
            output += [f"* [{escape_markdown(keyword_name)}]({keyword_xref})"]
        output += [""]
        # Render keywords
        output += ["## Keyword descriptions", ""]
        for keyword in keywords:
            output += render_keyword(keyword, section_xref)

    # Write output
    section_dir = output_dir / "/".join(section_path[:-1])
    section_dir.mkdir(exist_ok=True)
    filename = section_dir / f"{section_name}.md"
    filename.write_text("\n".join(output))
    num_files_written = 1

    # Process subsections
    for subsection in section.findall("SECTION"):
        subsection_path = (*section_path, get_name(subsection))
        num_files_written += process_section(subsection, subsection_path, output_dir)

    return num_files_written


# ======================================================================================
def render_keyword(keyword: lxml.etree._Element, section_xref: str) -> List[str]:
    # Find keyword names.
    keyword_names: List[str]
    if keyword.tag == "SECTION_PARAMETERS":
        keyword_names = ["SECTION_PARAMETERS"]
    elif keyword.tag == "DEFAULT_KEYWORD":
        keyword_names = ["DEFAULT_KEYWORD"]
    else:
        keyword_names = [get_text(name) for name in keyword.findall("NAME")]
    assert keyword_names

    # Find more keyword fields.
    default_value = get_text(keyword.find("DEFAULT_VALUE"))
    usage = get_text(keyword.find("USAGE"))
    description = get_text(keyword.find("DESCRIPTION"))
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
    output += [f"```{{py:data}}  {sanitize_name(keyword_names[0])}"]
    n_var_brackets = f"[{n_var}]" if n_var > 1 else ""
    output += [f":module: {section_xref}"]
    output += [f":type: '{data_type}{n_var_brackets}'"]
    if default_value:
        output += [f":value: '{default_value}'"]
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
            output += [f"* `{get_name(item)}` {escape_markdown(item_description)}"]
        output += [""]
    if references:
        citations = ", ".join([f"{{ref}}`{r}`" for r in references])
        output += [f"**References:** {citations}", ""]
    output += [f"{escape_markdown(description)} {github_link(location)}", ""]

    output += ["```", ""]  # Close py:data directive.

    return output


# ======================================================================================
def get_name(element: lxml.etree._Element) -> str:
    return get_text(element.find("NAME"))


# ======================================================================================
def get_text(element: Optional[lxml.etree._Element]) -> str:
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
    text = text.replace("__", "\_\_")
    text = text.replace("#", "\#")
    return text


# ======================================================================================
def github_link(location: str) -> str:
    if not location:
        return ""
    location_url = location.replace(":", "#L")
    github_url = f"https://github.com/cp2k/cp2k/blob/master/src/{location_url}"
    return f"<small>[[Edit on GitHub]({github_url})]</small>"


# ======================================================================================

main()

# EOF
