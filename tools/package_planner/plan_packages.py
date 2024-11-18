#!/usr/bin/env python3

# author: Ole Schuett

# A tool to help planning CP2K packages, listing currently violated dependencies if any

import os
import re
import sys
from os import path
from os.path import dirname, basename, normpath, abspath

# pre-compiled regular expressions
re_module = re.compile(r"(?:^|\n)\s*module\s+(\w+)\s.*\n\s*end\s*module", re.DOTALL)
re_program = re.compile(r"\n\s*end\s*program")
re_use = re.compile(r"\n\s*use\s+(\w+)")
re_incl1 = re.compile(r"\n\s*include\s+['\"](.+)['\"]")
re_incl2 = re.compile(r"\n#include\s+['\"](.+)['\"]")


# =============================================================================
def main():
    if len(sys.argv) != 2:
        print("Usage: plan_packages.py <plan_file>")
        sys.exit(1)

    planned_pkgs = eval(open(sys.argv[1], encoding="utf8").read())

    srcdir = "../../src"
    abs_srcdir = abspath(srcdir)

    src_files = []
    for root, _, files in os.walk(abs_srcdir):
        if root.endswith("/preprettify"):
            continue
        for fn in files:
            if fn[-2:] in (".F", ".c", ".cu"):
                src_files.append(path.join(root, fn))

    src_basenames = [basename(fn).rsplit(".", 1)[0] for fn in src_files]
    for bfn in src_basenames:
        if src_basenames.count(bfn) > 1:
            error("Multiple source files with the same basename: " + bfn)

    # parse source files
    parsed_files = dict()
    for fn in src_files:
        parse_file(parsed_files, fn)  # parses also included files
    print("Parsed %d source files" % len(parsed_files))

    # create table mapping fortran module-names to file-name
    mod2fn = dict()
    for fn in src_files:
        for m in parsed_files[fn]["module"]:
            if m in mod2fn:
                error('Multiple declarations of module "%s"' % m)
            mod2fn[m] = fn
    print("Created mod2fn table, found %d modules." % len(mod2fn))

    # check "one module per file"-convention
    for m, fn in mod2fn.items():
        if basename(fn) != m + ".F":
            error("Names of module and file do not match: " + fn)

    # check for circular dependencies
    for fn in parsed_files.keys():
        find_cycles(parsed_files, mod2fn, fn)

    # read existsing package manifests
    packages = dict()
    fn2pkg = dict()
    for fn in parsed_files.keys():
        p = normpath(dirname(fn))
        fn2pkg[basename(fn)] = p
        read_manifest(packages, p)

    # update with manifest with planned packages
    for pp in planned_pkgs:
        p = abspath(path.join(srcdir, pp["dirname"]))
        if p not in packages:
            packages[p] = {"problems": []}
        packages[p].update(pp)
        if "files" in pp:
            for fn in pp["files"]:
                fn2pkg[fn] = p
        if "requires+" in pp:
            packages[p]["requires"] += pp["requires+"]
        if "requires-" in pp:
            for i in pp["requires-"]:
                while i in packages[p]["requires"]:
                    packages[p]["requires"].remove(i)

    # process the manifests
    for p in packages.keys():
        process_manifest(packages, p)
        find_pkg_cycles(packages, p)
    print("Read %d package manifests" % len(packages))

    # check dependencies against package manifests
    n_deps = 0
    for fn in src_files:
        p = fn2pkg[basename(fn)]
        deps = collect_include_deps(parsed_files, fn)
        deps += [mod2fn[m] for m in collect_use_deps(parsed_files, fn) if m in mod2fn]
        n_deps += len(deps)
        for d in deps:
            dp = fn2pkg[basename(d)]
            msg = "%32s  ->  %-32s  %-15s " % (
                basename(fn),
                basename(d),
                "[src" + dp[len(abs_srcdir) :] + "]",
            )
            if dp not in packages[p]["allowed_deps"]:
                packages[p]["problems"].append(msg + "(requirement not listed)")
            if dp != p and "public" in packages[dp]:
                if basename(d) not in packages[dp]["public"]:
                    packages[p]["problems"].append(msg + "(file not public)")

    print("Checked %d dependencies\n" % n_deps)

    for p in sorted(packages.keys()):
        if len(packages[p]["problems"]) == 0:
            continue
        assert p.startswith(abs_srcdir)
        print(
            " " * 10
            + "====== Forbidden Dependencies in Package src%s ====="
            % p[len(abs_srcdir) :]
        )
        print("\n".join(packages[p]["problems"]) + "\n")


# =============================================================================
def parse_file(parsed_files, fn):
    if fn in parsed_files:
        return

    content = open(fn, encoding="utf8").read()

    # re.IGNORECASE is horribly expensive. Converting to lower-case upfront
    content_lower = content.lower()

    mods = re_module.findall(content_lower)
    prog = re_program.search(content_lower) is not None
    uses = re_use.findall(content_lower)
    incl1_iter = re_incl1.finditer(content_lower)  # fortran includes
    incls = [content[m.start(1) : m.end(1)] for m in incl1_iter]
    incls += re_incl2.findall(content)  # cpp includes, those are case-sensitiv

    # exclude included files from outside the source tree
    def incl_fn(i):
        return normpath(path.join(dirname(fn), i))

    existing_incl = [i for i in incls if path.exists(incl_fn(i))]

    # store everything in parsed_files cache
    parsed_files[fn] = {
        "module": mods,
        "program": prog,
        "use": uses,
        "include": existing_incl,
    }

    # parse included files
    for i in existing_incl:
        parse_file(parsed_files, incl_fn(i))


# =============================================================================
def read_manifest(packages, p):
    if p in packages:
        return

    fn = p + "/PACKAGE"
    if not path.exists(fn):
        error("Could not open PACKAGE manifest: " + fn)

    content = open(fn, encoding="utf8").read()
    packages[p] = eval(content)
    packages[p]["problems"] = []

    for r in packages[p]["requires"]:
        rp = normpath(path.join(p, r))
        read_manifest(packages, rp)


# =============================================================================
def process_manifest(packages, p):
    if "archive" not in packages[p]:
        packages[p]["archive"] = "libcp2k" + basename(p)
    packages[p]["allowed_deps"] = [normpath(p)]
    packages[p]["allowed_deps"] += [
        normpath(path.join(p, r)) for r in packages[p]["requires"]
    ]

    for r in packages[p]["requires"]:
        rp = normpath(path.join(p, r))
        if rp not in packages:
            error(
                "Unexpected package requirement: "
                + r
                + " for dir "
                + packages[p]["dirname"]
            )


# =============================================================================
def mod2modfile(m, mod_format):
    if mod_format == "no":
        return ""
    if mod_format == "lower":
        return m.lower() + ".mod"
    if mod_format == "upper":
        return m.upper() + ".mod"
    assert False  # modeps unknown


# =============================================================================
def collect_include_deps(parsed_files, fn):
    pf = parsed_files[fn]
    incs = []

    for i in pf["include"]:
        fn_inc = normpath(path.join(dirname(fn), i))
        if fn_inc in parsed_files:
            incs.append(fn_inc)
            incs += collect_include_deps(parsed_files, fn_inc)

    return sorted(set(incs))


# =============================================================================
def collect_use_deps(parsed_files, fn):
    pf = parsed_files[fn]
    uses = pf["use"]

    for i in pf["include"]:
        fn_inc = normpath(path.join(dirname(fn), i))
        if fn_inc in parsed_files:
            uses += collect_use_deps(parsed_files, fn_inc)

    return sorted(set(uses))


# =============================================================================
def find_cycles(parsed_files, mod2fn, fn, S=None):
    pf = parsed_files[fn]
    if "visited" in pf:
        return

    if S is None:
        S = []

    for m in pf["module"]:
        if m in S:
            i = S.index(m)
            error("Circular dependency: " + " -> ".join(S[i:] + [m]))
        S.append(m)

    for m in collect_use_deps(parsed_files, fn):
        if m in mod2fn:
            find_cycles(parsed_files, mod2fn, mod2fn[m], S)

    for m in pf["module"]:
        S.pop()

    pf["visited"] = True


# =============================================================================
def find_pkg_cycles(packages, p, S=None):
    if S is None:
        S = []

    if p in S:
        i = S.index(p)
        error("Circular package dependency: " + " -> ".join(S[i:] + [p]))
    S.append(p)

    for r in packages[p]["requires"]:
        d = normpath(path.join(p, r))
        find_pkg_cycles(packages, d, S)

    S.pop()


# =============================================================================
def error(msg):
    sys.stderr.write("error: %s\n" % msg)
    sys.exit(1)


# =============================================================================
main()
