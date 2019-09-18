#!/usr/bin/env python

from __future__ import print_function

import sys
import re
import tempfile
import os
from os import path
import logging
import argparse
import errno
import traceback

try:
    from hashlib import md5
except ImportError:
    from md5 import new as md5

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from prettify_cp2k import normalizeFortranFile
from prettify_cp2k import replacer

sys.path.append(path.join(path.dirname(path.abspath(__file__)), "fprettify"))
from fprettify import reformat_ffile, fparse_utils, log_exception


TO_UPCASE_RE = re.compile(
    r"""
(?P<toUpcase>
    \.(?:and|eq|eqv|false|ge|gt|le|lt|ne|neqv|not|or|true)\.
    |
    (?<![\w%#])  # do not match stmts/intrinsics midword or called as type-bound procedures
    (?<!%\ )  # do not match stmts/intrinsics when being called as type-bound procedure with a space in between
    (?<!subroutine\ )(?<!function\ )(?<!::\ )  # do not match stmts/intrinsics when used as procedure names
    (?:
        (?: # statements:
            a(?:llocat(?:able|e)|ssign(?:|ment))
            |c(?:a(?:ll|se)|haracter|lose|o(?:m(?:mon|plex)|nt(?:ains|inue))|ycle)
            |d(?:ata|eallocate|imension|o(?:|uble))
            |e(?:lse(?:|if|where)|n(?:d(?:|do|file|if)|try)|quivalence|x(?:it|ternal))
            |f(?:or(?:all|mat)|unction)
            |goto
            |i(?:f|mplicit|n(?:clude|quire|t(?:e(?:ger|nt|rface)|rinsic)))
            |logical
            |module
            |n(?:amelist|one|ullify)
            |o(?:nly|p(?:en|erator|tional))
            |p(?:a(?:rameter|use)|ointer|r(?:ecision|i(?:nt|vate)|o(?:cedure|gram))|ublic)
            |re(?:a[dl]|cursive|sult|turn|wind)
            |s(?:ave|e(?:lect|quence)|top|ubroutine)
            |t(?:arget|hen|ype)
            |use
            |w(?:h(?:ere|ile)|rite)
        )
        | (?: # intrinsic functions:
            a(?:bs|c(?:har|os)|djust[lr]|i(?:mag|nt)|ll(?:|ocated)|n(?:int|y)|s(?:in|sociated)|tan2?)
            |b(?:it_size|test)
            |c(?:eiling|har|mplx|o(?:njg|sh?|unt)|shift)
            |d(?:ate_and_time|ble|i(?:gits|m)|ot_product|prod)
            |e(?:oshift|psilon|xp(?:|onent))
            |f(?:loor|raction)
            |huge
            |i(?:a(?:char|nd)|b(?:clr|its|set)|char|eor|n(?:dex|t)|or|shftc?)
            |kind
            |l(?:bound|en(?:|_trim)|g[et]|l[et]|og(?:|10|ical))
            |m(?:a(?:tmul|x(?:|exponent|loc|val))|erge|in(?:|exponent|loc|val)|od(?:|ulo)|vbits)
            |n(?:earest|int|ot)
            |p(?:ack|r(?:e(?:cision|sent)|oduct))
            |r(?:a(?:dix|n(?:dom_(?:number|seed)|ge))|e(?:peat|shape)|rspacing)
            |s(?:ca(?:le|n)|e(?:lected_(?:int_kind|real_kind)|t_exponent)|hape|i(?:gn|nh?|ze)|p(?:acing|read)|qrt|um|ystem_clock)
            |t(?:anh?|iny|r(?:ans(?:fer|pose)|im))
            |u(?:bound|npack)
            |verify
        ) (?=\ *\()
    )
    (?![\w%])
)
""",
    flags=re.IGNORECASE | re.VERBOSE,
)

TO_UPCASE_OMP_RE = re.compile(
    r"""
(?<![\w%#])
(?P<toUpcase>
    (?:
        atomic|barrier|c(?:apture|ritical)|do|end|flush|if|master|num_threads|ordered|parallel|read
        |s(?:ection(?:|s)|ingle)|t(?:ask(?:|wait|yield)|hreadprivate)|update|w(?:orkshare|rite)|!\$omp
    )
    | (?:
        a|co(?:llapse|py(?:in|private))|default|fi(?:nal|rstprivate)|i(?:and|eor|or)|lastprivate
        |m(?:ax|ergeable|in)|n(?:one|owait)|ordered|private|reduction|shared|untied|\.(?:and|eqv|neqv|or)\.
    )
    | omp_(?:dynamic|max_active_levels|n(?:ested|um_threads)|proc_bind|s(?:tacksize|chedule)|thread_limit|wait_policy)
)
(?![\w%])
""",
    flags=re.IGNORECASE | re.VERBOSE,
)

LINE_PARTS_RE = re.compile(
    r"""
    (?P<commands>[^\"'!]*)
    (?P<comment>!.*)?
    (?P<string>
        (?P<qchar>[\"'])
        .*?
        (?P=qchar))?
""",
    re.VERBOSE,
)


def upcaseStringKeywords(line):
    """Upcases the fortran keywords, operators and intrinsic routines
    in line"""
    res = ""
    start = 0
    while start < len(line):
        m = LINE_PARTS_RE.match(line[start:])
        if not m:
            raise SyntaxError("Syntax error, open string")
        res = res + TO_UPCASE_RE.sub(
            lambda match: match.group("toUpcase").upper(), m.group("commands")
        )
        if m.group("comment"):
            res = res + m.group("comment")
        if m.group("string"):
            res = res + m.group("string")
        start = start + m.end()
    return res


def upcaseKeywords(infile, outfile, upcase_omp):
    """Writes infile to outfile with all the fortran keywords upcased"""

    for line in infile:
        line = upcaseStringKeywords(line)

        if upcase_omp and normalizeFortranFile.OMP_DIR_RE.match(line):
            line = TO_UPCASE_OMP_RE.sub(
                lambda match: match.group("toUpcase").upper(), line
            )

        outfile.write(line)


def prettifyFile(
    infile,
    filename,
    normalize_use,
    decl_linelength,
    decl_offset,
    reformat,
    indent,
    whitespace,
    upcase_keywords,
    upcase_omp,
    replace,
):
    """prettifyes the fortran source in infile into a temporary file that is
    returned. It can be the same as infile.
    if normalize_use normalizes the use statements (defaults to true)
    if upcase_keywords upcases the keywords (defaults to true)
    if replace does the replacements contained in replacer.py (defaults
    to false)

    does not close the input file"""
    max_pretty_iter = 5

    logger = logging.getLogger("fprettify-logger")

    if is_fypp(infile):
        logger.warning(
            "{}: fypp directives not fully supported, running only fprettify".format(
                filename
            )
        )
        replace = False
        normalize_use = False
        upcase_keywords = False

    # create a temporary file first as a copy of the input file
    inbuf = StringIO(infile.read())

    hash_prev = md5(inbuf.getvalue().encode("utf8"))

    for _ in range(max_pretty_iter):
        try:
            if replace:
                outbuf = StringIO()
                replacer.replaceWords(inbuf, outbuf)
                outbuf.seek(0)
                inbuf.close()
                inbuf = outbuf

            if reformat:  # reformat needs to be done first
                outbuf = StringIO()
                try:
                    reformat_ffile(
                        inbuf,
                        outbuf,
                        indent_size=indent,
                        whitespace=whitespace,
                        orig_filename=filename,
                    )
                except fparse_utils.FprettifyParseException as e:
                    log_exception(
                        e, "fprettify could not parse file, file is not prettified"
                    )
                    outbuf.close()
                    inbuf.seek(0)
                else:
                    outbuf.seek(0)
                    inbuf.close()
                    inbuf = outbuf

            if normalize_use:
                outbuf = StringIO()
                normalizeFortranFile.rewriteFortranFile(
                    inbuf,
                    outbuf,
                    indent,
                    decl_linelength,
                    decl_offset,
                    orig_filename=filename,
                )
                outbuf.seek(0)
                inbuf.close()
                inbuf = outbuf

            if upcase_keywords:
                outbuf = StringIO()
                upcaseKeywords(inbuf, outbuf, upcase_omp)
                outbuf.seek(0)
                inbuf.close()
                inbuf = outbuf

            hash_new = md5(inbuf.getvalue().encode("utf8"))

            if hash_prev.digest() == hash_new.digest():
                return inbuf

            hash_prev = hash_new

        except:
            logger.critical("error processing file '{}'".format(filename))
            raise

    else:
        raise RuntimeError(
            "Prettify did not converge in {} steps.".format(max_pretty_iter)
        )


def prettifyInplace(filename, backupdir=None, stdout=False, **kwargs):
    """Same as prettify, but inplace, replaces only if needed"""

    if filename == "stdin":
        infile = tempfile.TemporaryFile(mode="r+")
        infile.write(sys.stdin.read())
        infile.seek(0)
    else:
        infile = open(filename, "r")

    outfile = prettifyFile(infile=infile, filename=filename, **kwargs)

    if stdout:
        outfile.seek(0)
        sys.stdout.write(outfile.read())
        outfile.close()
        return

    if infile == outfile:
        infile.close()
        return

    infile.seek(0)
    outfile.seek(0)
    changed = True
    for line1, line2 in zip(infile, outfile):
        if line1 != line2:
            break
    else:
        changed = False

    if changed:
        if backupdir:
            bkName = path.join(backupdir, path.basename(filename))

            with open(bkName, "w") as fhandle:
                infile.seek(0)
                fhandle.write(infile.read())

        infile.close()  # close it here since we're going to overwrite it

        with open(filename, "w") as fhandle:
            outfile.seek(0)
            fhandle.write(outfile.read())

    else:
        infile.close()

    outfile.close()


def is_fypp(infile):
    FYPP_SYMBOLS = r"(#|\$|@)"
    FYPP_LINE = r"^\s*" + FYPP_SYMBOLS + r":"
    FYPP_INLINE = r"(" + FYPP_SYMBOLS + r"{|}" + FYPP_SYMBOLS + r")"
    FYPP_RE = re.compile(r"(" + FYPP_LINE + r"|" + FYPP_INLINE + r")")

    infile.seek(0)
    for line in infile.readlines():
        if FYPP_RE.search(line):
            return True

    infile.seek(0)
    return False


# based on https://stackoverflow.com/a/31347222
def argparse_add_bool_arg(parser, name, default, helptxt):
    dname = name.replace("-", "_")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument(
        "--{}".format(name), dest=dname, action="store_true", help=helptxt
    )
    group.add_argument("--no-{}".format(name), dest=dname, action="store_false")
    parser.set_defaults(**{dname: default})


# from https://stackoverflow.com/a/600612
def mkdir_p(p):
    try:
        os.makedirs(p)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and path.isdir(p):
            pass
        else:
            raise


# from https://stackoverflow.com/a/14981125
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def abspath(p):
    return path.abspath(path.expanduser(p))


def main(argv):
    parser = argparse.ArgumentParser(
        description="Auto-format F90 source files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""\
            If no files are given, stdin is used.
            Note: for editor integration, use options --no-normalize-use --no-report-errors""",
    )

    parser.add_argument("--indent", type=int, default=3)
    parser.add_argument("--whitespace", type=int, default=2, choices=range(0, 5))
    parser.add_argument("--decl-linelength", type=int, default=100)
    parser.add_argument("--decl-offset", type=int, default=50)
    parser.add_argument("--backup-dir", type=abspath, default=abspath("preprettify"))

    argparse_add_bool_arg(parser, "upcase", True, "Upcasing fortran keywords."),
    argparse_add_bool_arg(
        parser,
        "normalize-use",
        True,
        """\
        Sorting and alignment of variable declarations and USE statements, removal of unused list entries.
        The line length of declarations is controlled by --decl-linelength=n, the offset of the variable list
        is controlled by --decl-offset=n.""",
    )
    argparse_add_bool_arg(parser, "omp-upcase", True, "Upcasing OMP directives.")
    argparse_add_bool_arg(
        parser,
        "reformat",
        True,
        """\
        Auto-indentation, auto-alignment and whitespace formatting.
        Amount of whitespace controlled by --whitespace = 0, 1, 2.
        For indenting with a relative width of n columns specify --indent=n.
        For manual formatting of specific lines:
        - disable auto-alignment by starting line continuation with an ampersand '&'.
        - completely disable reformatting by adding a comment '!&'.
        For manual formatting of a code block, use:
        - start a manually formatted block with a '!&<' comment and close it with a '!&>' comment.""",
    )
    argparse_add_bool_arg(
        parser,
        "replace",
        True,
        "If requested the replacements performed by the replacer.py script are also performed. Note: these replacements are specific to CP2K.",
    )
    argparse_add_bool_arg(parser, "stdout", False, "write output to stdout")
    argparse_add_bool_arg(
        parser,
        "do-backup",
        False,
        "store backups of original files in backup-dir (--backup-dir option)",
    )
    argparse_add_bool_arg(parser, "report-errors", True, "report warnings and errors")
    argparse_add_bool_arg(parser, "debug", False, "increase log level to debug")

    parser.add_argument("files", metavar="file", type=str, nargs="*", default=["stdin"])

    args = parser.parse_args(argv)

    if args.do_backup and not (args.stdout or args.files == ["stdin"]):
        mkdir_p(args.backup_dir)

    failure = 0

    for filename in args.files:
        if not path.isfile(filename) and not filename == "stdin":
            eprint("file '{}' does not exist!".format(filename))
            failure += 1
            continue

        level = logging.CRITICAL

        if args.report_errors:
            if args.debug:
                level = logging.DEBUG
            else:
                level = logging.INFO

        logger = logging.getLogger("fprettify-logger")
        logger.setLevel(level)
        sh = logging.StreamHandler()
        sh.setLevel(level)
        formatter = logging.Formatter("%(levelname)s - %(message)s")
        sh.setFormatter(formatter)
        logger.addHandler(sh)

        try:
            prettifyInplace(
                filename,
                backupdir=args.backup_dir if args.do_backup else None,
                stdout=args.stdout or filename == "stdin",
                normalize_use=args.normalize_use,
                decl_linelength=args.decl_linelength,
                decl_offset=args.decl_offset,
                reformat=args.reformat,
                indent=args.indent,
                whitespace=args.whitespace,
                upcase_keywords=args.upcase,
                upcase_omp=args.omp_upcase,
                replace=args.replace,
            )
        except:
            eprint("-" * 60)
            traceback.print_exc(file=sys.stderr)
            eprint("-" * 60)
            eprint("Processing file '{}'".format(filename))
            failure += 1

    return failure > 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
