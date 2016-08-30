#!/usr/bin/env python

import sys
import re
import tempfile
import os
import os.path
import logging

try:
    from hashlib import md5
except ImportError:
    from md5 import new as md5

from formatting import normalizeFortranFile
from formatting import replacer
from formatting import reformatFortranFile
from formatting import selftest


OPERATORS_STR = r"\.(?:and|eqv?|false|g[et]|l[et]|n(?:e(?:|qv)|ot)|or|true)\."

KEYWORDS_STR = "(?:a(?:llocat(?:able|e)|ssign(?:|ment))|c(?:a(?:ll|se)|haracter|lose|o(?:m(?:mon|plex)|nt(?:ains|inue))|ycle)|d(?:ata|eallocate|imension|o(?:|uble))|e(?:lse(?:|if|where)|n(?:d(?:|do|file|if)|try)|quivalence|x(?:it|ternal))|f(?:or(?:all|mat)|unction)|goto|i(?:f|mplicit|n(?:clude|quire|t(?:e(?:ger|nt|rface)|rinsic)))|logical|module|n(?:amelist|one|ullify)|o(?:nly|p(?:en|erator|tional))|p(?:a(?:rameter|use)|ointer|r(?:ecision|i(?:nt|vate)|o(?:cedure|gram))|ublic)|re(?:a[dl]|cursive|sult|turn|wind)|s(?:ave|e(?:lect|quence)|top|ubroutine)|t(?:arget|hen|ype)|use|w(?:h(?:ere|ile)|rite))"

INTRINSIC_PROCSTR = r"(?:a(?:bs|c(?:har|os)|djust[lr]|i(?:mag|nt)|ll(?:|ocated)|n(?:int|y)|s(?:in|sociated)|tan2?)|b(?:it_size|test)|c(?:eiling|har|mplx|o(?:njg|sh?|unt)|shift)|d(?:ate_and_time|ble|i(?:gits|m)|ot_product|prod)|e(?:oshift|psilon|xp(?:|onent))|f(?:loor|raction)|huge|i(?:a(?:char|nd)|b(?:clr|its|set)|char|eor|n(?:dex|t)|or|shftc?)|kind|l(?:bound|en(?:|_trim)|g[et]|l[et]|og(?:|10|ical))|m(?:a(?:tmul|x(?:|exponent|loc|val))|erge|in(?:|exponent|loc|val)|od(?:|ulo)|vbits)|n(?:earest|int|ot)|p(?:ack|r(?:e(?:cision|sent)|oduct))|r(?:a(?:dix|n(?:dom_(?:number|seed)|ge))|e(?:peat|shape)|rspacing)|s(?:ca(?:le|n)|e(?:lected_(?:int_kind|real_kind)|t_exponent)|hape|i(?:gn|nh?|ze)|p(?:acing|read)|qrt|um|ystem_clock)|t(?:anh?|iny|r(?:ans(?:fer|pose)|im))|u(?:bound|npack)|verify)(?= *\()"

OMP_DIR = r"(?:atomic|barrier|c(?:apture|ritical)|do|end|flush|if|master|num_threads|ordered|parallel|read|s(?:ection(?:|s)|ingle)|t(?:ask(?:|wait|yield)|hreadprivate)|update|w(?:orkshare|rite)|!\$omp)"

OMP_CLAUSE = r"(?:a|co(?:llapse|py(?:in|private))|default|fi(?:nal|rstprivate)|i(?:and|eor|or)|lastprivate|m(?:ax|ergeable|in)|n(?:one|owait)|ordered|private|reduction|shared|untied|\.(?:and|eqv|neqv|or)\.)"

OMP_ENV = r"omp_(?:dynamic|max_active_levels|n(?:ested|um_threads)|proc_bind|s(?:tacksize|chedule)|thread_limit|wait_policy)"

# FIXME: does not correctly match operator '.op.' if it is not separated
# by whitespaces.
TO_UPCASE_RE = re.compile("(?<![A-Za-z0-9_%#])(?<!% )(?P<toUpcase>" + OPERATORS_STR +
                        "|" + KEYWORDS_STR + "|" + INTRINSIC_PROCSTR +
                        ")(?![A-Za-z0-9_%])", flags=re.IGNORECASE)
TO_UPCASE_OMP_RE = re.compile("(?<![A-Za-z0-9_%#])(?P<toUpcase>"
                           + OMP_DIR + "|" + OMP_CLAUSE + "|" + OMP_ENV +
                           ")(?![A-Za-z0-9_%])", flags=re.IGNORECASE)
LINE_PARTS_RE = re.compile("(?P<commands>[^\"'!]*)(?P<comment>!.*)?" +
                         "(?P<string>(?P<qchar>[\"']).*?(?P=qchar))?")


def upcaseStringKeywords(line):
    """Upcases the fortran keywords, operators and intrinsic routines
    in line"""
    res = ""
    start = 0
    while start < len(line):
        m = LINE_PARTS_RE.match(line[start:])
        if not m:
            raise SyntaxError("Syntax error, open string")
        res = res + TO_UPCASE_RE.sub(lambda match: match.group("toUpcase").upper(),
                                   m.group("commands"))
        if m.group("comment"):
            res = res + m.group("comment")
        if m.group("string"):
            res = res + m.group("string")
        start = start + m.end()
    return res


def upcaseOMP(line):
    """Upcases OpenMP stuff."""
    return TO_UPCASE_OMP_RE.sub(lambda match: match.group("toUpcase").upper(), line)


def upcaseKeywords(infile, outfile, upcase_omp):
    """Writes infile to outfile with all the fortran keywords upcased"""
    while 1:
        line = infile.readline()
        if not line:
            break
        line = upcaseStringKeywords(line)
        if upcase_omp:
            if normalizeFortranFile.OMP_DIR_RE.match(line):
                line = upcaseOMP(line)
        outfile.write(line)


def prettifyFile(infile, filename, normalize_use, decl_linelength, decl_offset,
                 reformat, indent, whitespace, upcase_keywords,
                 upcase_omp, replace):
    """prettifyes the fortran source in infile into a temporary file that is
    returned. It can be the same as infile.
    if normalize_use normalizes the use statements (defaults to true)
    if upcase_keywords upcases the keywords (defaults to true)
    if replace does the replacements contained in replacer.py (defaults
    to false)

    does not close the input file"""
    ifile = infile
    orig_filename = filename
    tmpfile = None
    max_pretty_iter = 5
    n_pretty_iter = 0

    while True:
        n_pretty_iter += 1
        hash_prev = md5()
        hash_prev.update(ifile.read().encode("utf8"))
        ifile.seek(0)
        try:
            if replace:
                tmpfile2 = tempfile.TemporaryFile(mode="w+")
                replacer.replaceWords(ifile, tmpfile2)
                tmpfile2.seek(0)
                if tmpfile:
                    tmpfile.close()
                tmpfile = tmpfile2
                ifile = tmpfile
            if reformat:  # reformat needs to be done first
                tmpfile2 = tempfile.TemporaryFile(mode="w+")
                reformatFortranFile.reformat_ffile(ifile, tmpfile2,
                                                   indent_size=indent, whitespace=whitespace,
                                                   orig_filename=orig_filename)
                tmpfile2.seek(0)
                if tmpfile:
                    tmpfile.close()
                tmpfile = tmpfile2
                ifile = tmpfile
            if normalize_use:
                tmpfile2 = tempfile.TemporaryFile(mode="w+")
                normalizeFortranFile.rewriteFortranFile(ifile, tmpfile2, indent,
                                                        decl_linelength, decl_offset,
                                                        orig_filename=orig_filename)
                tmpfile2.seek(0)
                if tmpfile:
                    tmpfile.close()
                tmpfile = tmpfile2
                ifile = tmpfile
            if upcase_keywords:
                tmpfile2 = tempfile.TemporaryFile(mode="w+")
                upcaseKeywords(ifile, tmpfile2, upcase_omp)
                tmpfile2.seek(0)
                if tmpfile:
                    tmpfile.close()
                tmpfile = tmpfile2
                ifile = tmpfile
            hash_next = md5()
            hash_next.update(ifile.read().encode("utf8"))
            ifile.seek(0)
            if hash_prev.digest() == hash_next.digest():
                return ifile
            elif n_pretty_iter >= max_pretty_iter:
                raise RuntimeError(
                    "Prettify did not converge in", max_pretty_iter, "steps.")
        except:
            logging.critical("error processing file '" + infile.name + "'\n")
            raise


def prettfyInplace(fileName, bkDir=None, stdout=False, **kwargs):
    """Same as prettify, but inplace, replaces only if needed"""

    if fileName == 'stdin':
        infile = os.tmpfile()
        infile.write(sys.stdin.read())
    else:
        infile = open(fileName, 'r')

    if stdout:
        outfile = prettifyFile(infile=infile, filename=fileName, **kwargs)
        outfile.seek(0)
        sys.stdout.write(outfile.read())
        outfile.close()
        return

    if bkDir and not os.path.exists(bkDir):
        os.mkdir(bkDir)
    if bkDir and not os.path.isdir(bkDir):
        raise Error("bk-dir must be a directory, was " + bkDir)

    outfile = prettifyFile(infile=infile, filename=fileName, **kwargs)
    if (infile == outfile):
        return
    infile.seek(0)
    outfile.seek(0)
    same = 1

    while 1:
        l1 = outfile.readline()
        l2 = infile.readline()
        if (l1 != l2):
            same = 0
            break
        if not l1:
            break
    if (not same):
        bkFile = None
        if bkDir:
            bkName = os.path.join(bkDir, os.path.basename(fileName))
            bName = bkName
            i = 0
            while os.path.exists(bkName):
                i += 1
                bkName = bName + "." + str(i)
            bkFile = file(bkName, "w")
        infile.seek(0)
        if bkFile:
            bkFile.write(infile.read())
            bkFile.close()
        outfile.seek(0)
        newFile = file(fileName, 'w')
        newFile.write(outfile.read())
        newFile.close()
    infile.close()
    outfile.close()


def main(argv=None):
    if argv is None:
        argv = sys.argv
    defaultsDict = {'upcase': 1, 'normalize-use': 1, 'omp-upcase': 1,
                    'decl-linelength': 100, 'decl-offset': 50,
                    'reformat': 1, 'indent': 3, 'whitespace': 1,
                    'replace': 1,
                    'stdout': 0,
                    'do-backup': 0,
                    'backup-dir': 'preprettify',
                    'report-errors': 1}

    usageDesc = ("usage:\nfprettify" +"""
    [--[no-]upcase] [--[no-]normalize-use] [--[no-]omp-upcase] [--[no-]replace]
    [--[no-]reformat] [--indent=3] [--whitespace=1] [--help]
    [--[no-]stdout] [--[no-]do-backup] [--backup-dir=bk_dir] [--[no-]report-errors] file1 [file2 ...]

    Auto-format F90 source file1, file2, ...:
    If no files are given, stdin is used.
    --normalize-use
             Sorting and alignment of variable declarations and USE statements, removal of unused list entries.
             The line length of declarations is controlled by --decl-linelength=n, the offset of the variable list
             is controlled by --decl-offset=n.
    --reformat
             Auto-indentation, auto-alignment and whitespace formatting.
             Amount of whitespace controlled by --whitespace = 0, 1, 2.
             For indenting with a relative width of n columns specify --indent=n.
             For manual formatting of specific lines:
             - disable auto-alignment by starting line continuation with an ampersand '&'.
             - completely disable reformatting by adding a comment '!&'.
             For manual formatting of a code block, use:
             - start a manually formatted block with a '!&<' comment and close it with a '!&>' comment.
    --upcase
             Upcasing fortran keywords.
    --omp-upcase
             Upcasing OMP directives.
    --replace
             If requested the replacements performed by the replacer.py script are also performed. Note: these replacements are specific to CP2K.
    --stdout
             write output to stdout
    --[no-]do-backup
             store backups of original files in backup-dir (--backup-dir option)
    --[no-]report-errors
             report warnings and errors

    Note: for editor integration, use options --no-normalize-use --no-report-errors

    Defaults:
    """ + str(defaultsDict))

    replace = None

    debug=False
    if "--help" in argv:
        sys.stderr.write(usageDesc + '\n')
        return(0)
    args = []
    for arg in argv[1:]:
        m = re.match(
            r"--(no-)?(normalize-use|upcase|omp-upcase|replace|reformat|stdout|do-backup|report-errors)", arg)
        if m:
            defaultsDict[m.groups()[1]] = not m.groups()[0]
        else:
            m = re.match(
                r"--(indent|whitespace|decl-linelength|decl-offset)=(.*)", arg)
            if m:
                defaultsDict[m.groups()[0]] = int(m.groups()[1])
            else:
                m = re.match(r"--(backup-dir)=(.*)", arg)
                if m:
                    path = os.path.abspath(os.path.expanduser(m.groups()[1]))
                    defaultsDict[m.groups()[0]] = path
                else:
                    if arg.startswith('--'):
                        sys.stderr.write('unknown option ' + arg + '\n')
                    else:
                        args.append(arg)
    bkDir = ''
    if defaultsDict['do-backup']:
        bkDir = defaultsDict['backup-dir']
    if bkDir and not os.path.exists(bkDir):
        # Another parallel running instance might just have created the
        # dir.
        try:
            os.mkdir(bkDir)
        except:
            assert(os.path.exists(bkDir))
    if bkDir and not os.path.isdir(bkDir):
        sys.stderr.write("bk-dir must be a directory" + '\n')
        sys.stderr.write(usageDesc + '\n')
    else:
        failure = 0
        if not args:
            args = ['stdin']
        for fileName in args:
            if not os.path.isfile(fileName) and not fileName == 'stdin':
                sys.stderr.write("file " + fileName + " does not exists!\n")
            else:
                stdout = defaultsDict['stdout'] or fileName == 'stdin'

                if defaultsDict['report-errors']:
                    if debug:
                        level=logging.DEBUG
                    else:
                        level=logging.INFO

                else:
                    level=logging.CRITICAL

                logging.basicConfig(stream=sys.stderr, level=level)

                try:
                    prettfyInplace(fileName, bkDir=bkDir,
                                   stdout=stdout,
                                   normalize_use=defaultsDict[
                                       'normalize-use'],
                                   decl_linelength=defaultsDict[
                                       'decl-linelength'],
                                   decl_offset=defaultsDict[
                                       'decl-offset'],
                                   reformat=defaultsDict['reformat'],
                                   indent=defaultsDict['indent'],
                                   whitespace=defaultsDict[
                                       'whitespace'],
                                   upcase_keywords=defaultsDict[
                                       'upcase'],
                                   upcase_omp=defaultsDict[
                                       'omp-upcase'],
                                   replace=defaultsDict['replace'])
                except:
                    failure += 1
                    import traceback
                    sys.stderr.write('-' * 60 + "\n")
                    traceback.print_exc(file=sys.stderr)
                    sys.stderr.write('-' * 60 + "\n")
                    sys.stderr.write(
                        "Processing file '" + fileName + "'\n")
        return(failure > 0)

#=========================================================================


def run_selftest():
    # create temporary file with example code
    fn = os.path.join(tempfile.gettempdir(), "prettify_selftest.F")
    ref = selftest.content
    f = open(fn, "w")
    f.write(ref)
    f.close()

    # call prettify
    rtn = main([sys.argv[0], fn])
    assert(rtn == 0)

    # check if file was altered
    result = open(fn).read()
    for i, (l1, l2) in enumerate(zip(result.split("\n"), ref.split("\n"))):
        if(l1 != l2):
            print("Error: Line %d is not invariant." % i)
            print("before: " + l1)
            print("after : " + l2)
            os.remove(fn)
            return(1)

    os.remove(fn)
    print("Prettify selftest passed.")
    return(0)

#=========================================================================
if(__name__ == '__main__'):
    if(len(sys.argv) == 2 and sys.argv[-1] == "--selftest"):
        rtn = run_selftest()
    else:
        rtn = main()

    sys.exit(rtn)
# EOF
