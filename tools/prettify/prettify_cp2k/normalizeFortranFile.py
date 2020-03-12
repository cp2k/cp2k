import sys
import re
import logging
from collections import deque
from io import StringIO

# Declare various RE snippets for building larger REs

# Valid fortran argument/variable name
VALID_NAME = r"[a-zA-Z][a-zA-Z_0-9]*"

# Set up valid type declaration
DECL_KIND_KW = rf"\( (?:kind = )? (?:{VALID_NAME}|[0-9]+) \)".replace(" ", r"\s*")
DECL_KIND_STAR = r"\s*\*\s*[0-9]+"
CHAR_KIND_SPEC = rf"(?:(?:len|kind) = )? (?:{VALID_NAME}|[0-9]+|\*|:)"
CHAR_KIND_SPEC = rf"\( {CHAR_KIND_SPEC}(?: , {CHAR_KIND_SPEC}) \)".replace(" ", r"\s*")
DECL_KIND_SPEC = rf"\s*(?:{DECL_KIND_STAR}|{DECL_KIND_KW})"
DECL_TYPE_DEF = rf"""(?:
(?:integer|real|complex)(?:{DECL_KIND_SPEC})?|
character(?:{DECL_KIND_STAR}|{CHAR_KIND_SPEC})?|
logical(?:{DECL_KIND_KW})?|type|class)"""

# String styles
DQUOTED = r'"[^"]*"'
SQUOTED = r"'[^']*'"
QUOTED = rf"(?:{DQUOTED}|{SQUOTED})"

# Stop at punctuation
NOT_PUNC = r"[^()\"'\[\]]"
NOT_PUNC_COMMA = r"[^()\"',\[\]]"
NOT_BRAC = r"[^()]"

# Enable new-style fortran arrays (square brackets)
BRAC_OPEN = r"[(\[]"
BRAC_CLOSE = r"[)\]]"
NOT_PAREN = r"[^()]"

# Nested parens for up to 3-level nesting (1st level's parens are outside this variable)
# Should more than 3 levels be necessary, need to reduplicate this deeper.
NESTED_PAREN_1TO3_CONTENTS = (
    rf"(?:{NOT_PAREN}+|\((?:{NOT_PAREN}+|\({NOT_PAREN}*\))*\))*"
)

# Variable statement
VAR_RE = re.compile(
    rf"""
\s*(?P<var>{VALID_NAME})
\s*(?P<rest>(:?\((?P<param>{NESTED_PAREN_1TO3_CONTENTS})\))?
\s*(?:=\s*(?P<value>(?:
{NOT_PUNC_COMMA}+|
{BRAC_OPEN}(?:{NOT_PUNC}|{BRAC_OPEN}{NOT_PUNC}*{BRAC_CLOSE}|{QUOTED})*{BRAC_CLOSE}|
{QUOTED})+))?)?
\s*(?:(?P<continue>,)|\n?)
\s*""",
    re.IGNORECASE | re.VERBOSE,
)

# Use statement
USE_PARSE_RE = re.compile(
    rf"""\s*
use(\s+|(?P<intrinsic>\s*,\s*Intrinsic\s+::\s*))
(?P<module>{VALID_NAME})
(?P<only>\s*,\s*only\s*:)?\s*
(?P<imports>.*)$""",
    flags=re.IGNORECASE | re.VERBOSE,
)

COMMON_USES_RE = re.compile('^#include\s*"([^"]*(cp_common_uses.f90|base_uses.f90))"')
LOCAL_NAME_RE = re.compile(
    rf"\s*(?P<localName>{VALID_NAME})(?:\s*=>\s*{VALID_NAME})?\s*$", re.VERBOSE
)
VAR_DECL_RE = re.compile(
    rf"""
\s*(?P<type>{DECL_TYPE_DEF})
\s*(?P<parameters>\({NESTED_PAREN_1TO3_CONTENTS}\))?
\s*(?P<attributes>(?:\s*,\s*[a-zA-Z_0-9]+(?:\s*\({NESTED_PAREN_1TO3_CONTENTS}\))?)+)?\s*
(?P<dpnt>::)?
(?P<vars>[^\n]+)\n?
""",
    re.IGNORECASE | re.VERBOSE,
)
INDENT_SIZE = 2
DECL_LINELENGTH = 100
DECL_OFFSET = 50

OMP_DIR_RE = re.compile(r"^\s*(!\$omp)", re.IGNORECASE)
OMP_RE = re.compile(r"^\s*(!\$)", re.IGNORECASE)

FCT_RE = re.compile(
    r"^([^\"'!]* )?FUNCTION\s+\w+\s*(\(.*\))?(\s*RESULT\s*\(\w+\))?\s*;?\s*$",
    re.IGNORECASE,
)

SUBR_RE = re.compile(
    r"^([^\"'!]* )?SUBROUTINE\s+\w+\s*(\(.*\))?\s*;?\s*$", re.IGNORECASE
)

END_RE = re.compile(r" *end\s*(?:subroutine|function)", re.IGNORECASE)
START_ROUTINE_RE = re.compile(
    r"^([^\"'!]* )?(?P<kind>subroutine|function) +(?P<name>[a-zA-Z_][a-zA-Z_0-9]*) *(?:\((?P<arguments>[^()]*)\))? *(?:result *\( *(?P<result>[a-zA-Z_][a-zA-Z_0-9]*) *\))? *(?:bind *\([^()]+\))? *\n?",
    re.IGNORECASE,
)  # $
typeBeginRe = re.compile(
    r" *(?P<type>integer(?: *\* *[0-9]+)?|logical|character(?: *\* *[0-9]+)?|real(?: *\* *[0-9]+)?|complex(?: *\* *[0-9]+)?|type)[,( ]",
    re.IGNORECASE,
)
attributeRe = re.compile(
    r" *, *(?P<attribute>[a-zA-Z_0-9]+) *(?:\( *(?P<param>(?:[^()]+|\((?:[^()]+|\([^()]*\))*\))*)\))? *",
    re.IGNORECASE,
)
IGNORE_RE = re.compile(r" *(?:|implicit +none *)$", re.IGNORECASE)
INTERFACE_START_RE = re.compile(r" *interface *$", re.IGNORECASE)
INTERFACE_END_RE = re.compile(r" *end +interface *$", re.IGNORECASE)

INCLUDE_RE = re.compile(r"#? *include +[\"'](?P<file>.+)[\"'] *$", re.IGNORECASE)

CONTAINS_RE = re.compile(r" *contains *$", re.IGNORECASE)

NON_WORD_RE = re.compile(r"(\(/|/\)|[^-+a-zA-Z0-9_.])")

STR_RE = re.compile(r"('[^'\n]*'|\"[^\"\n]*\")")

MODULE_N_RE = re.compile(
    r".*:: *moduleN *= *(['\"])[a-zA-Z_0-9]+\1", flags=re.IGNORECASE
)

MODULE_RE = re.compile(
    r" *(?:module|program) +(?P<moduleName>[a-zA-Z_][a-zA-Z_0-9]*) *(?:!.*)?$",
    flags=re.IGNORECASE,
)


class InputStreamError(Exception):
    pass


class CharFilter(object):
    """
    An iterator to wrap the iterator returned by `enumerate`
    and ignore comments and characters inside strings
    """

    def __init__(self, it):
        self._it = it
        self._instring = ""

    def __iter__(self):
        return self

    def __next__(self):
        """ python 3 version"""
        pos, char = next(self._it)
        if not self._instring and char == "!":
            raise StopIteration

        # detect start/end of a string
        if char == '"' or char == "'":
            if self._instring == char:
                self._instring = ""
            elif not self._instring:
                self._instring = char

        if self._instring:
            return self.__next__()

        return (pos, char)

    def next(self):
        """ python 2 version"""
        pos, char = self._it.next()
        if not self._instring and char == "!":
            raise StopIteration

        # detect start/end of a string
        if char == '"' or char == "'":
            if self._instring == char:
                self._instring = ""
            elif not self._instring:
                self._instring = char

        if self._instring:
            return self.next()

        return (pos, char)


class InputStream(object):
    """
    Class to read logical Fortran lines from a Fortran file.
    """

    def __init__(self, infile):
        self.line_buffer = deque([])
        self.infile = infile
        self.line_nr = 0

    def nextFortranLine(self):
        """Reads a group of connected lines (connected with &, separated by newline or semicolon)
        returns a touple with the joined line, and a list with the original lines.
        Doesn't support multiline character constants!
        """
        lineRe = re.compile(
            # $
            r"(?:(?P<preprocessor>#.*\n?)| *(&)?(?P<core>(?:!\$|[^&!\"']+|\"[^\"]*\"|'[^']*')*)(?P<continue>&)? *(?P<comment>!.*)?\n?)",
            re.IGNORECASE,
        )
        joinedLine = ""
        comments = []
        lines = []
        continuation = 0

        while True:
            if not self.line_buffer:
                line = self.infile.readline().replace("\t", 8 * " ")
                self.line_nr += 1
                # convert OMP-conditional fortran statements into normal fortran statements
                # but remember to convert them back
                is_omp_conditional = False
                omp_indent = 0
                if OMP_RE.match(line):
                    omp_indent = len(line) - len(line.lstrip(" "))
                    line = OMP_RE.sub("", line, count=1)
                    is_omp_conditional = True
                line_start = 0
                for pos, char in CharFilter(enumerate(line)):
                    if char == ";" or pos + 1 == len(line):
                        self.line_buffer.append(
                            omp_indent * " "
                            + "!$" * is_omp_conditional
                            + line[line_start : pos + 1]
                        )
                        omp_indent = 0
                        is_omp_conditional = False
                        line_start = pos + 1
                if line_start < len(line):
                    # line + comment
                    self.line_buffer.append(
                        "!$" * is_omp_conditional + line[line_start:]
                    )

            if self.line_buffer:
                line = self.line_buffer.popleft()

            if not line:
                break

            lines.append(line)
            m = lineRe.match(line)
            if not m or m.span()[1] != len(line):
                # FIXME: does not handle line continuation of
                # omp conditional fortran statements
                # starting with an ampersand.
                raise InputStreamError("unexpected line format:" + repr(line))
            if m.group("preprocessor"):
                if len(lines) > 1:
                    raise InputStreamError(
                        "continuation to a preprocessor line not supported "
                        + repr(line)
                    )
                comments.append(line)
                break
            coreAtt = m.group("core")
            if OMP_RE.match(coreAtt) and joinedLine.strip():
                # remove omp '!$' for line continuation
                coreAtt = OMP_RE.sub("", coreAtt, count=1).lstrip()
            joinedLine = joinedLine.rstrip("\n") + coreAtt
            if coreAtt and not coreAtt.isspace():
                continuation = 0
            if m.group("continue"):
                continuation = 1
            if line.lstrip().startswith("!") and not OMP_RE.search(line):
                comments.append(line.rstrip("\n"))
            elif m.group("comment"):
                comments.append(m.group("comment"))
            else:
                comments.append("")
            if not continuation:
                break
        return (joinedLine, comments, lines)


def parseRoutine(inFile, logger):
    """Parses a routine"""

    routine = {
        "preRoutine": [],
        "core": [],
        "strippedCore": [],
        "begin": [],
        "end": [],
        "preDeclComments": [],
        "declarations": [],
        "declComments": [],
        "postDeclComments": [],
        "parsedDeclarations": [],
        "postRoutine": [],
        "kind": None,
        "name": None,
        "arguments": None,
        "result": None,
        "interfaceCount": 0,
        "use": [],
    }
    stream = InputStream(inFile)
    while True:
        (jline, _, lines) = stream.nextFortranLine()
        if len(lines) == 0:
            break
        if FCT_RE.match(jline) or SUBR_RE.match(jline):
            break
        routine["preRoutine"].extend(lines)
        m = INCLUDE_RE.match(lines[0])
        if m:
            try:
                subF = open(m.group("file"), "r")
                subStream = InputStream(subF)
                while True:
                    (subjline, _, sublines) = subStream.nextFortranLine()
                    if not sublines:
                        break
                    routine["strippedCore"].append(subjline)
                subF.close()
            except:
                import traceback

                logger.debug(
                    "error trying to follow include '{}', this might lead to the removal of used variables".format(
                        m.group("file")
                    )
                )

                if logger.isEnabledFor(logging.DEBUG):
                    traceback.print_exc()

    if jline:
        routine["begin"] = lines
        m = START_ROUTINE_RE.match(jline)
        if not m or m.span()[1] != len(jline):
            raise SyntaxError("unexpected subroutine start format:" + repr(lines))
        routine["name"] = m.group("name")
        routine["kind"] = m.group("kind")
        if m.group("arguments") and m.group("arguments").strip():
            routine["arguments"] = list(
                x.strip() for x in m.group("arguments").split(",")
            )
        if m.group("result"):
            routine["result"] = m.group("result")
        if (not routine["result"]) and (routine["kind"].lower() == "function"):
            routine["result"] = routine["name"]

    while True:
        (jline, comment_list, lines) = stream.nextFortranLine()
        comments = "\n".join(_ for _ in comment_list)
        if len(lines) == 0:
            break
        if lines[0].lower().startswith("#include"):
            break
        if not IGNORE_RE.match(jline):
            if typeBeginRe.match(jline):
                if routine["postDeclComments"]:
                    routine["declComments"].extend(routine["postDeclComments"])
                routine["postDeclComments"] = []

            if typeBeginRe.match(jline):
                m = VAR_DECL_RE.match(jline)
                if m.group("type").lower() == "type" and not m.group("parameters"):
                    break
                if not m or m.span()[1] != len(jline):
                    raise SyntaxError("unexpected type format:" + repr(jline))
                decl = {
                    "type": m.group("type"),
                    "parameters": None,
                    "attributes": [],
                    "vars": [],
                }
                if m.group("parameters"):
                    decl["parameters"] = (
                        m.group("parameters").replace(" ", "").replace(",", ", ")
                    )
                str = m.group("attributes")
                while str:
                    m2 = attributeRe.match(str)
                    if not m2:
                        raise SyntaxError(
                            "unexpected attribute format "
                            + repr(str)
                            + " in "
                            + repr(lines)
                        )
                    decl["attributes"].append(
                        m2.group().replace(" ", "").replace(",", ", ")[2:]
                    )
                    str = str[m2.span()[1] :]
                str = m.group("vars")
                while True:
                    m2 = VAR_RE.match(str)
                    if not m2:
                        raise SyntaxError(
                            "unexpected var format " + repr(str) + " in " + repr(lines)
                        )
                    var = m2.group("var")
                    if m2.group("param"):
                        var += "(" + m2.group("param") + ")"
                    if m2.group("value"):
                        var += " = "
                        var += m2.group("value")
                    decl["vars"].append(var)
                    str = str[m2.span()[1] :]
                    if not m2.group("continue"):
                        if str:
                            raise SyntaxError(
                                "error parsing vars (leftover="
                                + repr(str)
                                + ") in "
                                + repr(lines)
                            )
                        break
                routine["parsedDeclarations"].append(decl)
            elif INTERFACE_START_RE.match(jline):
                istart = lines
                interfaceDeclFile = StringIO()
                while True:
                    (jline, _, lines) = stream.nextFortranLine()
                    if INTERFACE_END_RE.match(jline):
                        iend = lines
                        break
                    interfaceDeclFile.writelines(lines)
                interfaceDeclFile = StringIO(interfaceDeclFile.getvalue())
                iroutines = []
                while True:
                    iroutine = parseRoutine(interfaceDeclFile, logger)
                    if not iroutine["kind"]:
                        if len(iroutines) == 0:
                            interfaceDeclFile.seek(0)
                            raise SyntaxError(
                                "error parsing interface:"
                                + repr(interfaceDeclFile.read())
                            )
                        iroutines[-1]["postRoutine"].extend(iroutine["preRoutine"])
                        break
                    iroutines.append(iroutine)
                for iroutine in iroutines:
                    routine["interfaceCount"] += 1
                    decl = {
                        "type": "z_interface%02d" % (routine["interfaceCount"]),
                        "parameters": None,
                        "attributes": [],
                        "vars": [iroutine["name"]],
                        "iroutine": iroutine,
                        "istart": istart,
                        "iend": iend,
                    }
                    routine["parsedDeclarations"].append(decl)
            elif USE_PARSE_RE.match(jline):
                routine["use"].append("".join(lines))
            else:
                break
        routine["declarations"].append("".join(lines))
        if (
            len(routine["parsedDeclarations"]) == 0
            and len(routine["use"]) == 0
            and not re.match(" *implicit +none *$", jline, re.IGNORECASE)
        ):
            routine["preDeclComments"].append("".join(lines))
        else:
            routine["postDeclComments"].append(comments)

    while len(lines) > 0:
        if END_RE.match(jline):
            routine["end"] = lines
            break
        routine["strippedCore"].append(jline)
        routine["core"].append("".join(lines))
        if CONTAINS_RE.match(lines[0]):
            break
        m = INCLUDE_RE.match(lines[0])
        if m:
            try:
                subF = open(m.group("file"), "r")
                subStream = InputStream(subF)
                while True:
                    (subjline, _, sublines) = subStream.nextFortranLine()
                    if not sublines:
                        break
                    routine["strippedCore"].append(subjline)
                subF.close()
            except:
                import traceback

                logger.debug(
                    "error trying to follow include '{}', this might lead to the removal of used variables".format(
                        m.group("file")
                    )
                )

                if logger.isEnabledFor(logging.DEBUG):
                    traceback.print_exc()

        (jline, _, lines) = stream.nextFortranLine()
    return routine


def findWord(word, text, options=re.IGNORECASE):
    """Returns the position of word in text or -1 if not found.
    A match is valid only if it is a whole word (i.e. findWord('try','retry')
    returns false)"""
    wordRe = re.compile(
        "(?<![a-zA-Z_0-9%])"
        + word
        + "(?![a-zA-Z_0-9])|(?<=[0-9.]_)"
        + word
        + "(?![a-zA-Z_0-9])",
        options,
    )
    m = wordRe.search(text)
    if m:
        pos = m.span()[0]
    else:
        pos = -1
    return pos


def enforceDeclDependecies(declarations):
    """enforces the dependencies between the vars
    and compacts the declarations, returns the variables needed by other variables"""
    idecl = 0
    ii = 0
    while idecl < len(declarations):
        typeParam = "".join(declarations[idecl]["attributes"])
        if declarations[idecl]["parameters"]:
            typeParam += " " + declarations[idecl]["parameters"]
        typeParam = typeParam.lower()

        ivar = 0
        while ivar < len(declarations[idecl]["vars"]):
            moved = 0
            m = VAR_RE.match(declarations[idecl]["vars"][ivar])
            if not m:
                raise SyntaxError(
                    "could not match var " + repr(declarations[idecl]["vars"][ivar])
                )
            rest = m.group("rest")
            rest = rest.lower()
            if rest:
                for ivar2 in range(ivar + 1, len(declarations[idecl]["vars"])):
                    m = VAR_RE.match(declarations[idecl]["vars"][ivar2])
                    if findWord(m.group("var").lower(), rest) != -1:
                        moved = ivar2 + 1
            if moved:
                declarations[idecl]["vars"][moved:moved] = [
                    declarations[idecl]["vars"][ivar]
                ]
                del declarations[idecl]["vars"][ivar]
            else:
                for idecl2 in range(idecl + 1, len(declarations)):
                    for ivar2 in range(len(declarations[idecl2]["vars"])):
                        ii += 1
                        if ii > 100000:
                            raise Error("could not enforce all constraints")
                        m = VAR_RE.match(declarations[idecl2]["vars"][ivar2])
                        if (
                            ivar == 0
                            and findWord(m.group("var").lower(), typeParam) != -1
                        ):
                            declarations.insert(idecl2 + 1, declarations[idecl])
                            del declarations[idecl]
                            ivar = 0
                            moved = 1
                            break
                        if rest and findWord(m.group("var").lower(), rest) != -1:
                            if len(declarations[idecl]["vars"]) > 1:
                                newDecl = {}
                                newDecl.update(declarations[idecl])
                                newDecl["vars"] = [declarations[idecl]["vars"][ivar]]
                                declarations.insert(idecl2 + 1, newDecl)
                                del declarations[idecl]["vars"][ivar]
                            else:
                                declarations.insert(idecl2 + 1, declarations[idecl])
                                del declarations[idecl]
                                ivar = 0
                            moved = 1
                            break
                    if moved:
                        break
            if not moved:
                ivar += 1
        idecl += 1

    for i in range(len(declarations) - 1, 0, -1):
        if (
            declarations[i]["normalizedType"].lower()
            == declarations[i - 1]["normalizedType"].lower()
        ):
            declarations[i - 1]["vars"].extend(declarations[i]["vars"])
            del declarations[i]


def sortDeclarations(declarations):
    """sorts, compacts declarations and respects dependencies
    normalizedType has to be defined for the declarations"""

    declarations.sort(key=lambda x: x["normalizedType"].lower())

    for i in range(len(declarations) - 1, 0, -1):
        if (
            declarations[i]["normalizedType"].lower()
            == declarations[i - 1]["normalizedType"].lower()
        ):
            declarations[i - 1]["vars"].extend(declarations[i]["vars"])
            del declarations[i]

    for decl in declarations:
        decl["vars"].sort(key=lambda x: x.lower())
    enforceDeclDependecies(declarations)


def writeRoutine(routine, outFile):
    """writes the given routine to outFile"""
    outFile.writelines(routine["preRoutine"])
    outFile.writelines(routine["begin"])
    outFile.writelines(routine["declarations"])
    outFile.writelines(routine["core"])
    outFile.writelines(routine["end"])
    outFile.writelines(routine["postRoutine"])


def writeInCols(dLine, indentCol, maxCol, indentAtt, file):
    """writes out the strings (trying not to cut them) in dLine up to maxCol
    indenting each newline with indentCol.
    The '&' of the continuation line is at maxCol.
    indentAtt is the actual intent, and the new indent is returned"""

    maxSize = maxCol - indentCol - 1
    tol = min(maxSize / 6, 6) + indentCol
    for fragment in dLine:
        if indentAtt + len(fragment) < maxCol:
            file.write(fragment)
            indentAtt += len(fragment)
        elif len(fragment.lstrip()) <= maxSize:
            file.write("&\n" + (" " * indentCol))
            file.write(fragment.lstrip())
            indentAtt = indentCol + len(fragment.lstrip())
        else:
            sPieces = STR_RE.split(fragment)
            for sPiece in sPieces:
                if sPiece and (not (sPiece[0] == '"' or sPiece[0] == "'")):
                    subPieces = NON_WORD_RE.split(sPiece)
                else:
                    subPieces = [sPiece]
                for subPiece in subPieces:
                    if indentAtt == indentCol:
                        file.write(subPiece.lstrip())
                        indentAtt += len(subPiece.lstrip())
                    elif indentAtt < tol or indentAtt + len(subPiece) < maxCol:
                        file.write(subPiece)
                        indentAtt += len(subPiece)
                    else:
                        file.write("&\n" + (" " * indentCol))
                        file.write(subPiece.lstrip())
                        indentAtt = indentCol + len(subPiece.lstrip())
    return indentAtt


def writeCompactDeclaration(declaration, file):
    """Writes a declaration in a compact way"""
    d = declaration
    if "iroutine" in d.keys():
        file.writelines(d["istart"])
        writeRoutine(d["iroutine"], file)
        file.writelines(d["iend"])
    else:
        if len(d["vars"]) > 0:
            decl = " " * INDENT_SIZE * 2 + d["type"]
            if d["parameters"]:  # do not drop empty parameter lists?
                decl += d["parameters"]
            if d["attributes"]:
                for a in d["attributes"]:
                    decl += ", " + a
            decl += " :: "

            dLine = [decl]
            for var in d["vars"]:
                cur_len = sum([len(l) for l in dLine])
                if len(dLine) > 1 and cur_len + len(var) > 600:
                    writeInCols(dLine, 3 * INDENT_SIZE, DECL_LINELENGTH, 0, file)
                    file.write("\n")
                    dLine = [decl]
                if len(dLine) > 1:
                    dLine[-1] += ", "
                dLine.append(var)
            writeInCols(dLine, 3 * INDENT_SIZE, DECL_LINELENGTH, 0, file)
            file.write("\n")


def writeExtendedDeclaration(declaration, file):
    """Writes a declaration in a nicer way (using more space)"""
    d = declaration
    if len(d["vars"]) == 0:
        return
    if "iroutine" in d.keys():
        file.writelines(d["istart"])
        writeRoutine(d["iroutine"], file)
        file.writelines(d["iend"])
    else:
        dLine = []
        dLine.append(" " * INDENT_SIZE * 2 + d["type"])
        if d["parameters"]:  # do not drop empty parameter lists?
            dLine.append(d["parameters"])
        if d["attributes"]:
            for a in d["attributes"]:
                dLine[-1:] = [dLine[-1] + ", "]
                dLine.append(a)

        indentAtt = writeInCols(
            dLine, 3 * INDENT_SIZE, DECL_OFFSET + 1 + 2 * INDENT_SIZE, 0, file
        )
        file.write(" " * (DECL_OFFSET + 2 * INDENT_SIZE - indentAtt))
        file.write(" :: ")
        indentAtt = DECL_OFFSET + 8

        dLine = []
        for var in d["vars"][:-1]:
            dLine.append(var + ", ")
        dLine.append(d["vars"][-1])

        writeInCols(
            dLine, DECL_OFFSET + 4 + 2 * INDENT_SIZE, DECL_LINELENGTH, indentAtt, file
        )
        file.write("\n")


def writeDeclarations(parsedDeclarations, file):
    """Writes the declarations to the given file"""
    for d in parsedDeclarations:
        maxLenVar = 0
        totalLen = 0
        for v in d["vars"]:
            maxLenVar = max(maxLenVar, len(v))
            totalLen += len(v)
        if maxLenVar > 30 or totalLen > DECL_LINELENGTH - 4:
            writeCompactDeclaration(d, file)
        else:
            writeExtendedDeclaration(d, file)


def cleanDeclarations(routine, logger):
    """cleans up the declaration part of the given parsed routine
    removes unused variables"""

    if routine["core"]:
        if CONTAINS_RE.match(routine["core"][-1]):
            logger.debug(
                "routine %s contains other routines, declarations not cleaned"
                % (routine["name"])
            )
            return
    nullifyRe = re.compile(
        r" *nullify *\(([^()]+)\) *\n?", re.IGNORECASE | re.MULTILINE
    )

    if not routine["kind"]:
        return
    if routine["core"]:
        if re.match(" *type *[a-zA-Z_]+ *$", routine["core"][0], re.IGNORECASE):
            logger.debug(
                "routine %s contains local types, not fully cleaned" % (routine["name"])
            )
        if re.match(" *import+ *$", routine["core"][0], re.IGNORECASE):
            logger.debug(
                "routine %s contains import, not fully cleaned" % (routine["name"])
            )
    if re.search("^#", "".join(routine["declarations"]), re.MULTILINE):
        logger.debug(
            "routine %s declarations contain preprocessor directives, declarations not cleaned"
            % (routine["name"])
        )
        return
    try:
        rest = "".join(routine["strippedCore"]).lower()
        nullifys = ",".join(nullifyRe.findall(rest))
        rest = nullifyRe.sub("", rest)
        paramDecl = []
        decls = []
        for d in routine["parsedDeclarations"]:
            d["normalizedType"] = d["type"]
            if d["parameters"]:
                d["normalizedType"] += d["parameters"]
            if d["attributes"]:
                d["attributes"].sort(key=lambda x: x.lower())
                d["normalizedType"] += ", "
                d["normalizedType"] += ", ".join(d["attributes"])
            if any(a.lower() == "parameter" for a in d["attributes"]):
                paramDecl.append(d)
            else:
                decls.append(d)

        sortDeclarations(paramDecl)
        sortDeclarations(decls)
        has_routinen = 0
        pos_routinep = -1
        for d in paramDecl:
            for i in range(len(d["vars"])):
                v = d["vars"][i]
                m = VAR_RE.match(v)
                lowerV = m.group("var").lower()
                if lowerV == "routinen":
                    has_routinen = 1
                    d["vars"][i] = "routineN = '" + routine["name"] + "'"
                elif lowerV == "routinep":
                    pos_routinep = i
                    d["vars"][i] = "routineP = moduleN//':'//routineN"
            if not has_routinen and pos_routinep >= 0:
                d["vars"].insert(pos_routinep, "routineN = '" + routine["name"] + "'")

        if routine["arguments"]:
            routine["lowercaseArguments"] = list(
                x.lower() for x in routine["arguments"]
            )
        else:
            routine["lowercaseArguments"] = []
        if routine["result"]:
            routine["lowercaseArguments"].append(routine["result"].lower())
        argDeclDict = {}
        localDecl = []
        for d in decls:
            localD = {}
            localD.update(d)
            localD["vars"] = []
            argD = None
            for v in d["vars"]:
                m = VAR_RE.match(v)
                lowerV = m.group("var").lower()
                if lowerV in routine["lowercaseArguments"]:
                    argD = {}
                    argD.update(d)
                    argD["vars"] = [v]
                    if lowerV in argDeclDict.keys():
                        raise SyntaxError(
                            "multiple declarations not supported. var="
                            + v
                            + " declaration="
                            + str(d)
                            + "routine="
                            + routine["name"]
                        )
                    argDeclDict[lowerV] = argD
                else:
                    pos = findWord(lowerV, rest)
                    if pos != -1:
                        localD["vars"].append(v)
                    else:
                        if findWord(lowerV, nullifys) != -1:
                            if not rmNullify(lowerV, routine["core"]):
                                raise SyntaxError(
                                    "could not remove nullify of "
                                    + lowerV
                                    + " as expected, routine="
                                    + routine["name"]
                                )
                        logger.info(
                            "removed var %s in routine %s\n" % (lowerV, routine["name"])
                        )
            if len(localD["vars"]):
                localDecl.append(localD)
        argDecl = []
        for arg in routine["lowercaseArguments"]:
            if arg in argDeclDict.keys():
                argDecl.append(argDeclDict[arg])
            else:
                logger.debug(
                    "warning, implicitly typed argument '{}' in routine '{}'".format(
                        arg, routine["name"]
                    )
                )

        if routine["kind"].lower() == "function":
            aDecl = argDecl[:-1]
        else:
            aDecl = argDecl

        # try to have arg/param/local, but checks for dependencies arg/param
        # and param/local
        argDecl.extend(paramDecl)
        enforceDeclDependecies(argDecl)
        splitPos = 0
        for i in range(len(argDecl) - 1, -1, -1):
            if not any(a.lower() == "parameter" for a in argDecl[i]["attributes"]):
                splitPos = i + 1
                break
        paramDecl = argDecl[splitPos:]
        argDecl = argDecl[:splitPos]
        paramDecl.extend(localDecl)
        enforceDeclDependecies(paramDecl)
        splitPos = 0
        for i in range(len(paramDecl) - 1, -1, -1):
            if any(a.lower() == "parameter" for a in paramDecl[i]["attributes"]):
                splitPos = i + 1
                break
        localDecl = paramDecl[splitPos:]
        paramDecl = paramDecl[:splitPos]

        newDecl = StringIO()
        for comment in routine["preDeclComments"]:
            newDecl.write(comment)
        newDecl.writelines(routine["use"])
        writeDeclarations(argDecl, newDecl)
        if argDecl and paramDecl:
            newDecl.write("\n")
        writeDeclarations(paramDecl, newDecl)
        if (argDecl or paramDecl) and localDecl:
            newDecl.write("\n")
        writeDeclarations(localDecl, newDecl)
        if argDecl or paramDecl or localDecl:
            newDecl.write("\n")
        wrote = 0
        for comment in routine["declComments"]:
            if comment.strip():
                newDecl.write(comment.strip())
                newDecl.write("\n")
                wrote = 1
        if wrote:
            newDecl.write("\n")
        routine["declarations"] = [newDecl.getvalue()]
    except:
        if "name" in routine.keys():
            logger.critical("exception cleaning routine " + routine["name"])
        logger.critical("parsedDeclartions={}".format(routine["parsedDeclarations"]))
        raise

    newDecl = StringIO()
    if routine["postDeclComments"]:
        comment_start = 0
        for comment in routine["postDeclComments"]:
            if comment.strip():
                break
            else:
                comment_start += 1

        for comment in routine["postDeclComments"][comment_start:]:
            newDecl.write(comment)
            newDecl.write("\n")
        routine["declarations"][0] += newDecl.getvalue()


def rmNullify(var, strings):
    removed = 0
    var = var.lower()
    nullifyRe = re.compile(r" *nullify *\(", re.IGNORECASE)
    nullify2Re = re.compile(
        r"(?P<nullif> *nullify *\()(?P<vars>[^()!&]+)\)", re.IGNORECASE
    )

    for i in range(len(strings) - 1, -1, -1):
        line = strings[i]
        comments = []
        if nullifyRe.match(line) and findWord(var, line) != -1:
            core = ""
            comments = []
            for l in line.splitlines():
                pos = l.find("&")
                pos2 = l.find("!")
                if pos == -1:
                    if pos2 == -1:
                        core += l
                    else:
                        core += l[:pos2]
                        comments.append(l[pos2:] + "\n")
                else:
                    core += l[:pos]
                    if pos2 != -1:
                        comments.append(l[pos2:] + "\n")
            m = nullify2Re.match(core)
            if not m:
                raise SyntaxError(
                    "could not match nullify to " + repr(core) + "in" + repr(line)
                )
            allVars = []
            vars = m.group("vars")
            v = list(s.strip() for s in vars.split(","))
            removedNow = 0
            for j in range(len(v) - 1, -1, -1):
                if findWord(var, v[j].lower()) != -1:
                    del v[j]
                    removedNow = 1
            if removedNow:
                if len(v) == 0:
                    if not comments:
                        del strings[i]
                    else:
                        strings[i] = "".join(comments)
                else:
                    for j in range(len(v) - 1):
                        v[j] += ", "
                    v[-1] += ")"
                    newS = StringIO()
                    v.insert(0, m.group("nullif"))
                    writeInCols(v, len(v[0]) - len(v[0].lstrip()) + 5, 77, 0, newS)
                    newS.write("\n")
                    if comments:
                        for c in comments:
                            newS.write(c)
                    strings[i] = newS.getvalue()
                removed += 1
    return removed


def parseUse(inFile):
    """Parses the use statements in inFile
    The parsing stops at the first non use statement.
    Returns something like:
    ([{'module':'module1','only':['el1','el2=>el3']},...],
     '! comment1\\n!comment2...\\n',
     'last line (the line that stopped the parsing)')
    """
    lineNr = 0
    preComments = []
    modules = []
    origLines = []
    commonUses = ""
    stream = InputStream(inFile)
    while True:
        (jline, comment_list, lines) = stream.nextFortranLine()
        comments = "\n".join(_ for _ in comment_list if _)
        lineNr = lineNr + len(lines)
        if not lines:
            break
        origLines.append("".join(lines))
        # parse use
        m = USE_PARSE_RE.match(jline)
        if m:
            useAtt = {"module": m.group("module"), "comments": []}

            if m.group("only"):
                useAtt["only"] = list(s.strip() for s in m.group("imports").split(","))
            else:
                useAtt["renames"] = list(
                    s.strip() for s in m.group("imports").split(",")
                )
                if useAtt["renames"] == [""]:
                    del useAtt["renames"]
            if comments:
                useAtt["comments"].append(comments)
            # add use to modules
            modules.append(useAtt)
        elif jline and not jline.isspace():
            break
        else:
            if comments and COMMON_USES_RE.match(comments):
                commonUses += "".join(lines)
            elif len(modules) == 0:
                preComments.append(("".join(lines)))
            elif comments:
                modules[-1]["comments"].append(comments)

    return {
        "modules": modules,
        "preComments": preComments,
        "commonUses": commonUses,
        "postLine": "".join(lines),
        "origLines": origLines[:-1],
    }


def normalizeModules(modules):
    """Sorts the modules and their export and removes duplicates.
    renames aren't sorted correctly"""
    # orders modules
    modules.sort(key=lambda x: x["module"])
    for i in range(len(modules) - 1, 0, -1):
        if modules[i]["module"].lower() == modules[i - 1]["module"].lower():
            if not ("only" in modules[i - 1].keys() and "only" in modules[i].keys()):
                raise SyntaxError(
                    "rejoining of module "
                    + str(modules[i]["module"])
                    + " failed as at least one of the use is not a use ...,only:"
                )
            modules[i - 1]["only"].extend(modules[i]["only"])
            del modules[i]
    # orders imports
    for m in modules:
        if "only" in m.keys():
            m["only"].sort()
            for i in range(len(m["only"]) - 1, 0, -1):
                if m["only"][i - 1].lower() == m["only"][i].lower():
                    del m["only"][i]


def writeUses(modules, outFile):
    """Writes the use declaration using a long or short form depending on how
    many only statements there are"""
    for m in modules:
        if "only" in m.keys() and len(m["only"]) > 8:
            writeUseShort(m, outFile)
        else:
            writeUseLong(m, outFile)


def writeUseLong(m, outFile):
    """Writes a use declaration in a nicer, but longer way"""
    if "only" in m.keys():
        outFile.write(
            INDENT_SIZE * " "
            + "USE "
            + m["module"]
            + ","
            + "ONLY: ".rjust(38 - len(m["module"]))
        )
        if m["only"]:
            outFile.write(m["only"][0])
        for i in range(1, len(m["only"])):
            outFile.write(",&\n" + "".ljust(43 + INDENT_SIZE) + m["only"][i])
    else:
        outFile.write(INDENT_SIZE * " " + "USE " + m["module"])
        if "renames" in m.keys() and m["renames"]:
            outFile.write("," + "".ljust(38) + m["renames"][0])
            for i in range(1, len(m["renames"])):
                outFile.write(",&\n" + "".ljust(43 + INDENT_SIZE) + m["renames"][i])
    if m["comments"]:
        outFile.write("\n")
        outFile.write("\n".join(m["comments"]))
    outFile.write("\n")


def writeUseShort(m, file):
    """Writes a use declaration in a compact way"""
    uLine = []
    if "only" in m.keys():
        file.write(
            INDENT_SIZE * " "
            + "USE "
            + m["module"]
            + ","
            + "ONLY: &\n".rjust(40 - len(m["module"]))
        )
        for k in m["only"][:-1]:
            uLine.append(k + ", ")
        uLine.append(m["only"][-1])
        uLine[0] = " " * (5 + INDENT_SIZE) + uLine[0]
    elif "renames" in m.keys() and m["renames"]:
        uLine.append(INDENT_SIZE * " " + "USE " + m["module"] + ", ")
        for k in m["renames"][:-1]:
            uLine.append(k + ", ")
        uLine.append(m["renames"][-1])
    else:
        uLine.append(INDENT_SIZE * " " + "USE " + m["module"])
    writeInCols(uLine, 5 + INDENT_SIZE, DECL_LINELENGTH, 0, file)
    if m["comments"]:
        file.write("\n")
        file.write("\n".join(m["comments"]))
    file.write("\n")


def prepareImplicitUses(modules):
    """Transforms a modulesDict into an implictUses (dictionary of module names
    each containing a dictionary with the only, and the special key '_WHOLE_'
    wich is true if the whole mosule is implicitly present"""
    mods = {}
    for m in modules:
        m_name = m["module"].lower()
        if m_name not in mods.keys():
            mods[m["module"]] = {"_WHOLE_": 0}
        m_att = mods[m_name]
        if "only" in m.keys():
            for k in m["only"]:
                m = LOCAL_NAME_RE.match(k)
                if not m:
                    raise SyntaxError("could not parse use only:" + repr(k))
                impAtt = m.group("localName").lower()
                m_att[impAtt] = 1
        else:
            m_att["_WHOLE_"] = 1
    return mods


def cleanUse(modulesDict, rest, implicitUses, logger):
    """Removes the unneded modules (the ones that are not used in rest)"""

    exceptions = {}
    modules = modulesDict["modules"]
    rest = rest.lower()

    for i in range(len(modules) - 1, -1, -1):
        m_att = {}
        m_name = modules[i]["module"].lower()
        if implicitUses and m_name in implicitUses.keys():
            m_att = implicitUses[m_name]
        if "_WHOLE_" in m_att.keys() and m_att["_WHOLE_"]:
            logger.info("removed USE of module " + m_name)
            del modules[i]
        elif "only" in modules[i].keys():
            els = modules[i]["only"]
            for j in range(len(els) - 1, -1, -1):
                m = LOCAL_NAME_RE.match(els[j])
                if not m:
                    raise SyntaxError("could not parse use only:" + repr(els[j]))
                impAtt = m.group("localName").lower()
                if impAtt in m_att.keys():
                    logger.info("removed USE " + m_name + ", only: " + repr(els[j]))
                    del els[j]
                elif impAtt not in exceptions.keys():
                    if findWord(impAtt, rest) == -1:
                        logger.info("removed USE " + m_name + ", only: " + repr(els[j]))
                        del els[j]
            if len(modules[i]["only"]) == 0:
                if modules[i]["comments"]:
                    modulesDict["preComments"] += [
                        x + "\n" for x in modules[i]["comments"]
                    ]
                del modules[i]


def resetModuleN(moduleName, lines):
    "resets the moduleN variable to the module name in the lines lines"
    for i in range(len(lines)):
        lines[i] = MODULE_N_RE.sub(
            " " * INDENT_SIZE
            + "CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = '"
            + moduleName
            + "'",
            lines[i],
        )


def rewriteFortranFile(
    inFile, outFile, indent, decl_linelength, decl_offset, orig_filename=None
):
    """rewrites the use statements and declarations of inFile to outFile.
    It sorts them and removes the repetitions."""
    import os.path

    global INDENT_SIZE
    global DECL_OFFSET
    global DECL_LINELENGTH
    INDENT_SIZE = indent
    DECL_OFFSET = decl_offset
    DECL_LINELENGTH = decl_linelength

    coreLines = []
    while True:
        line = inFile.readline()
        if not line:
            break
        if line[0] == "#":
            coreLines.append(line)
        outFile.write(line)
        m = MODULE_RE.match(line)
        if m:
            if not orig_filename:
                orig_filename = inFile.name
            fn = os.path.basename(orig_filename).rsplit(".", 1)[0]
            break

    logger = logging.LoggerAdapter(
        logging.getLogger("prettify-logger"), {"ffilename": orig_filename}
    )

    modulesDict = parseUse(inFile)
    routines = []
    coreLines.append(modulesDict["postLine"])
    routine = parseRoutine(inFile, logger)
    coreLines.extend(routine["preRoutine"])

    if m:
        resetModuleN(m.group("moduleName"), routine["preRoutine"])

    routines.append(routine)

    while routine["kind"]:
        routine = parseRoutine(inFile, logger)
        routines.append(routine)

    for routine in routines:
        cleanDeclarations(routine, logger)  # in-place modification of 'routine'
        coreLines.extend(routine["declarations"])
        coreLines.extend(routine["strippedCore"])

    rest = "".join(coreLines)
    nonStPrep = 0

    for line in modulesDict["origLines"]:
        if re.search("^#", line) and not COMMON_USES_RE.match(line):
            logger.debug("noMatch " + repr(line) + "\n")  # what does it mean?
            nonStPrep = 1

    if nonStPrep:
        logger.debug("use statements contains preprocessor directives, not cleaning")
        outFile.writelines(modulesDict["origLines"])
    else:
        implicitUses = None
        if modulesDict["commonUses"]:
            inc_fn = COMMON_USES_RE.match(modulesDict["commonUses"]).group(1)
            inc_absfn = os.path.join(os.path.dirname(orig_filename), inc_fn)
            try:
                with open(inc_absfn, "r") as fhandle:
                    implicitUsesRaw = parseUse(fhandle)
                implicitUses = prepareImplicitUses(implicitUsesRaw["modules"])
            except:
                logger.critical(
                    "failed to parse use statements contained in common uses precompiler file {}".format(
                        inc_absfn
                    )
                )
                raise

        cleanUse(modulesDict, rest, implicitUses, logger)
        normalizeModules(modulesDict["modules"])
        outFile.writelines(modulesDict["preComments"])
        writeUses(modulesDict["modules"], outFile)
        outFile.write(modulesDict["commonUses"])
        if modulesDict["modules"]:
            outFile.write("\n")

    outFile.write(modulesDict["postLine"])

    for routine in routines:
        writeRoutine(routine, outFile)
