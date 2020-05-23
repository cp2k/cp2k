#!/usr/bin/env python3

# author: Ole Schuett

import re, sys
from os import path
from os.path import dirname, basename, normpath
import os

# pre-compiled regular expressions
re_module = re.compile(r"(?:^|\n)\s*module\s+(\w+)\s.*\n\s*end\s*module", re.DOTALL)
re_useonly = re.compile(r"\n\s*use\s+(\w+)\s*,\s*only\s*:(.*)(?=\n)")
re_use = re.compile(r"\n\s*use\s+(\w+)\s*(?=\n)")
re_pub = re.compile(r"\n\s*public\s*::\s*(.*)(?=\n)")


# =============================================================================
def main():
    src_dir = "../src"
    parsed_files = {}
    for root, _, files in os.walk(src_dir):
        if "preprettify" in root:
            continue
        for fn in files:
            if not fn.endswith(".F"):
                continue
            absfn = path.join(root, fn)

            parsed_files[absfn] = parse_file(absfn)

    all_used_symboles = set()
    for p in parsed_files.values():
        for m, syms in p["use"]:
            for s in syms:
                all_used_symboles.add(m + "@" + s)

    n = 0
    for fn, p in parsed_files.items():
        if len(p["mod"]) != 1:
            continue
        m = p["mod"][0]
        if m + "@*" in all_used_symboles:
            continue
        unused = []
        for s in p["pub"]:
            if m + "@" + s not in all_used_symboles:
                unused.append(s)
                n += 1
        if len(unused) > 0:
            # print("%s USElessly declares PUBLIC: "%fn+ ", ".join(unused)+"\n")
            clean_publics(fn, unused)

    # print("Found %d unUSEd PUBLIC symbols."%n)


# =============================================================================
def clean_publics(fn, unused):
    content = open(fn).read()
    new_content = ""
    active = False
    protected = False
    new_public = []
    for line in content.split("\n"):
        if line.strip().startswith("!API"):
            protected = True
        if re.match(r"\s*public\s*::.*", line, re.IGNORECASE):
            if protected:
                protected = False
            else:
                active = True

        if not active:
            new_content += line + "\n"
            continue

        prefix, symbols, comment = re.match("^(.*::)?([^!]*)(!.*)?$", line).groups()
        old_symbols = symbols.strip(" &,").split(",")
        new_symbols = []
        for s in old_symbols:
            s = s.strip()
            if s.lower() not in unused:
                new_symbols.append(s)

        if len(new_symbols) > 0:
            new_public.append((", ".join(new_symbols), comment))

        without_comment = re.sub("!.*", "", line).strip()
        if len(without_comment) > 0 and without_comment[-1] != "&":

            # flush new_public
            for i, entry in enumerate(new_public):
                if i == 0:
                    new_content += "  PUBLIC :: "
                else:
                    new_content += "            "
                new_content += entry[0]
                if i < len(new_public) - 1:
                    new_content += ",&"
                if entry[1]:
                    new_content += "  " + entry[1]
                new_content += "\n"

            active = False
            new_public = []

    new_content = new_content[:-1]  # remove last line break

    if new_content != content:
        print("Fixed: ", fn)
        f = open(fn, "w")
        f.write(new_content)


# =============================================================================
def parse_file(fn):
    # print("parsing "+fn)
    content = open(fn).read()
    # re.IGNORECASE is horribly expensive. Converting to lower-case upfront
    content = content.lower()
    content = re.sub("!.*\n", "\n", content)
    content = re.sub("&\s*\n", "", content)
    content = re.sub("&", " ", content)

    mods = re_module.findall(content)

    uses = []
    matches = re_use.findall(content)
    for m in matches:
        uses.append((m.strip(), ("*",)))
        if m.strip() not in ("iso_c_binding", "f77_blas",):
            print("missing ONLY-clause: ", fn, m)

    matches = re_useonly.findall(content)
    for m in matches:
        syms = [p.split("=>")[-1].strip() for p in m[1].split(",")]
        uses.append((m[0].strip(), syms))

    publics = []
    matches = re_pub.findall(content)
    for m in matches:
        publics += [p.strip() for p in m.split(",")]

    return {"mod": mods, "use": uses, "pub": publics}


# ===============================================================================
main()
