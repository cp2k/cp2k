"""
  Impose white space conventions and indentation based on scopes / subunits

  normalization of white spaces supported for following operators:
  - relational operators:
    .EQ. .NE. .LT. .LE. .GT. .GE.
    ==   /=   <    <=    >   >=
  - logical operators:
    .AND. .OR. .EQV. .NEQV.
    .NOT.
  - bracket delimiters
  - commas and semicolons:
  - arithmetic operators:
    *  /  **  +  -
  - other operators:
    %  - (sign)  = (function argument)
    = (assignment)  => (pointer assignment)

  supported criteria for alignment / indentation:
   Fortran lines:
   - if, else, endif
   - do, enddo
   - select case, case, end select
   - subroutine, end subroutine
   - function, end function
   - module, end module
   - program, end program
   - interface, end interface
   - type, end type
   Actual lines (parts of Fortran lines separated by linebreaks):
   - bracket delimiters (.), (/./), and [.]
   - assignments by value = and pointer =>.

  LIMITATIONS
  - assumes that all subunits are explicitly ended within same file, no treatment of #include statements
  - can not deal with f77 constructs (files are ignored)
"""

import re
import sys
from formatting.normalizeFortranFile import useParseRe, typeRe, InputStream, CharFilter, ompRe, ompDirRe

#=========================================================================
# constants, mostly regular expressions

RE_FLAGS = re.IGNORECASE  # all regex should be case insensitive

FORTRAN_DEFAULT_ERROR_MESSAGE = " Syntax error - this script can not handle invalid Fortran files."
FORMATTER_ERROR_MESSAGE = " Wrong usage of formatting-specific directives '&', '!&', '!&<' or '!&>'."

EOL_STR = r"\s*;?\s*$"  # end of fortran line
EOL_SC = r"\s*;\s*$"  # whether line is ended with semicolon
SOL_STR = r"^\s*"  # start of fortran line

F77_STYLE = re.compile(r"^\s*\d", RE_FLAGS)

# regular expressions for parsing statements that start, continue or end a
# subunit:
IF_RE = re.compile(
    SOL_STR + r"(\w+\s*:)?\s*IF\s*\(.+\)\s*THEN" + EOL_STR, RE_FLAGS)
ELSE_RE = re.compile(
    SOL_STR + r"ELSE(\s*IF\s*\(.+\)\s*THEN)?" + EOL_STR, RE_FLAGS)
ENDIF_RE = re.compile(SOL_STR + r"END\s*IF(\s+\w+)?" + EOL_STR, RE_FLAGS)

DO_RE = re.compile(SOL_STR + r"(\w+\s*:)?\s*DO(" + EOL_STR + r"|\s)", RE_FLAGS)
ENDDO_RE = re.compile(SOL_STR + r"END\s*DO(\s+\w+)?" + EOL_STR, RE_FLAGS)

SELCASE_RE = re.compile(
    SOL_STR + r"SELECT\s*CASE\s*\(.+\)" + EOL_STR, RE_FLAGS)
CASE_RE = re.compile(SOL_STR + r"CASE\s*(\(.+\)|DEFAULT)" + EOL_STR, RE_FLAGS)
ENDSEL_RE = re.compile(SOL_STR + r"END\s*SELECT" + EOL_STR, RE_FLAGS)

SUBR_RE = re.compile(
    r"^[^\"'!]*SUBROUTINE\s+\w+\s*(\(.*\))?" + EOL_STR, RE_FLAGS)
ENDSUBR_RE = re.compile(
    SOL_STR + r"END\s*SUBROUTINE(\s+\w+)?" + EOL_STR, RE_FLAGS)

FCT_RE = re.compile(
    r"^[^\"'!]*FUNCTION\s+\w+\s*(\(.*\))?(\s*RESULT\s*\(\w+\))?" + EOL_STR, RE_FLAGS)
ENDFCT_RE = re.compile(
    SOL_STR + r"END\s*FUNCTION(\s+\w+)?" + EOL_STR, RE_FLAGS)

MOD_RE = re.compile(SOL_STR + r"MODULE\s+\w+" + EOL_STR, RE_FLAGS)
ENDMOD_RE = re.compile(SOL_STR + r"END\s*MODULE(\s+\w+)?" + EOL_STR, RE_FLAGS)

TYPE_RE = re.compile(
    SOL_STR + r"TYPE(\s*,\s*BIND\s*\(\s*C\s*\))?(\s*::\s*|\s+)\w+" + EOL_STR, RE_FLAGS)
ENDTYPE_RE = re.compile(SOL_STR + r"END\s*TYPE(\s+\w+)?" + EOL_STR, RE_FLAGS)

PROG_RE = re.compile(SOL_STR + r"PROGRAM\s+\w+" + EOL_STR, RE_FLAGS)
ENDPROG_RE = re.compile(
    SOL_STR + r"END\s*PROGRAM(\s+\w+)?" + EOL_STR, RE_FLAGS)

INTERFACE_RE = re.compile(SOL_STR + r"INTERFACE(\s+\w+)?" + EOL_STR, RE_FLAGS)
ENDINTERFACE_RE = re.compile(
    SOL_STR + r"END\s*INTERFACE(\s+\w+)?" + EOL_STR, RE_FLAGS)

CONTAINS_RE = re.compile(SOL_STR + r"CONTAINS" + EOL_STR, RE_FLAGS)

PUBLIC_RE = re.compile(SOL_STR + r"PUBLIC\s*::")

# intrinsic statements with parenthesis notation that are not functions
INTR_STMTS_PAR = "(ALLOCATE|DEALLOCATE|REWIND|BACKSPACE|INQUIRE|OPEN|CLOSE|WRITE|READ|FORALL|WHERE|NULLIFY)"

# regular expressions for parsing linebreaks
LINEBREAK_STR = r"(&)[\s]*(?:!.*)?$"

# regular expressions for parsing operators
# Note: +/- in real literals and sign operator is ignored
PLUSMINUS_RE = re.compile(
    r"(?<=[\w\)\]])(?<![\d\.]\w)\s*(\+|-)\s*", RE_FLAGS)
REL_OP_RE = re.compile(
    r"\s*(\.(?:EQ|NE|LT|LE|GT|GE)\.|(?:==|\/=|<(?!=)|<=|(?<!=)>(?!=)|>=))\s*", RE_FLAGS)
LOG_OP_RE = re.compile(r"\s*(\.(?:AND|OR|EQV|NEQV)\.)\s*", RE_FLAGS)

# regular expressions for parsing delimiters
DEL_OPEN_STR = r"(\(\/?|\[)"
DEL_OPEN_RE = re.compile(DEL_OPEN_STR, RE_FLAGS)
DEL_CLOSE_STR = r"(\/?\)|\])"
DEL_CLOSE_RE = re.compile(DEL_CLOSE_STR, RE_FLAGS)

# empty line regex, treats omp statements as exception
EMPTY_RE = re.compile(SOL_STR + r"(![^\$].*)?$", RE_FLAGS)

# two-sided operators
LR_OPS_RE = [REL_OP_RE, LOG_OP_RE, PLUSMINUS_RE]

# markups to deactivate formatter
NO_ALIGN_RE = re.compile(SOL_STR + r"&\s*[^\s*]+")

# combine regex that define subunits
NEW_SCOPE_RE = [IF_RE, DO_RE, SELCASE_RE, SUBR_RE,
                FCT_RE, MOD_RE, PROG_RE, INTERFACE_RE, TYPE_RE]
CONTINUE_SCOPE_RE = [ELSE_RE, None, CASE_RE, CONTAINS_RE,
                     CONTAINS_RE, CONTAINS_RE, CONTAINS_RE, None, None]
END_SCOPE_RE = [ENDIF_RE, ENDDO_RE, ENDSEL_RE, ENDSUBR_RE,
                ENDFCT_RE, ENDMOD_RE, ENDPROG_RE, ENDINTERFACE_RE, ENDTYPE_RE]

#=========================================================================


class F90Indenter(object):
    """
    Parses encapsulation of subunits / scopes line by line and updates the indentation.
    """

    def __init__(self, filename):
        self._scope_storage = []
        self._indent_storage = [0]
        self._filename = filename
        self._line_indents = []
        self._aligner = F90Aligner(filename)

    def process_lines_of_fline(self, f_line, lines, rel_ind, rel_ind_con, line_nr, manual_lines_indent=None):
        """
        Process all lines that belong to a Fortran line `f_line`, impose a relative indent of `rel_ind` for
        current Fortran line, and `rel_ind_con` for line continuation. By default line continuations are
        auto-aligned by F90Aligner - manual offsets can be set by manual_lines_indents.
        """

        self._line_indents = [0] * len(lines)
        br_indent_list = [0] * len(lines)
        line_indents = self._line_indents
        scopes = self._scope_storage
        indents = self._indent_storage
        filename = self._filename

        # check statements that start new scope
        is_new = False
        valid_new = False

        debug = False

        for new_n, newre in enumerate(NEW_SCOPE_RE):
            if newre.search(f_line) and not END_SCOPE_RE[new_n].search(f_line):
                what_new = new_n
                is_new = True
                valid_new = True
                scopes.append(what_new)
                if debug:
                    print(f_line)
                break

        # check statements that continue scope
        is_con = False
        valid_con = False
        for con_n, conre in enumerate(CONTINUE_SCOPE_RE):
            if conre is not None and conre.search(f_line):
                what_con = con_n
                is_con = True
                if len(scopes) > 0:
                    what = scopes[-1]
                    if what == what_con:
                        valid_con = True
                        if debug:
                            print(f_line)
                        break

        # check statements that end scope
        is_end = False
        valid_end = False
        for end_n, endre in enumerate(END_SCOPE_RE):
            if endre.search(f_line):
                what_end = end_n
                is_end = True
                if len(scopes) > 0:
                    what = scopes.pop()
                    if what == what_end:
                        valid_end = True
                        if debug:
                            print(f_line)
                        break

        # deal with line breaks
        if not manual_lines_indent:
            self._aligner.process_lines_of_fline(
                f_line, lines, rel_ind_con, line_nr)
            br_indent_list = self._aligner.get_lines_indent()
        else:
            br_indent_list = manual_lines_indent

        for pos in range(0, len(lines) - 1):
            line_indents[pos + 1] = br_indent_list[pos + 1]

        if is_new:
            if not valid_new:
                raise SyntaxError(filename + ':' + str(line_nr) +
                                  ':' + FORTRAN_DEFAULT_ERROR_MESSAGE)
            else:
                line_indents = [ind + indents[-1] for ind in line_indents]
                old_ind = indents[-1]

                rel_ind += old_ind  # prevent originally unindented do / if blocks from being indented
                indents.append(rel_ind)

        elif is_con:
            if not valid_con:
                raise SyntaxError(filename + ':' + str(line_nr) +
                                  ':' + FORTRAN_DEFAULT_ERROR_MESSAGE)
            else:
                line_indents = [ind + indents[-2] for ind in line_indents]

        elif is_end:
            if not valid_end:
                raise SyntaxError(filename + ':' + str(line_nr) +
                                  ':' + FORTRAN_DEFAULT_ERROR_MESSAGE)
            else:
                line_indents = [ind + indents[-2] for ind in line_indents]
                indents.pop()
        else:
            line_indents = [ind + indents[-1] for ind in line_indents]

        self._line_indents = line_indents

    def get_fline_indent(self):
        """
        after processing, retrieve the indentation of the full Fortran line.
        """
        return self._indent_storage[-1]

    def get_lines_indent(self):
        """
        after processing, retrieve the indents of all line parts.
        """
        return self._line_indents

#=========================================================================


class F90Aligner(object):
    """
    Alignment of continuations of a broken line, based on the following heuristics
    if line break in brackets:     We are parsing the level of nesting and align to most inner bracket delimiter.
    else if line is an assignment: alignment to '=' or '=>'.
                                   note: assignment operator recognized as any '=' that is not
                                   part of another operator and that is not enclosed in bracket
    else if line is a declaration: alignment to '::'
    else                           default indent
    """

    def __init__(self, filename):
        self._filename = filename
        self.__init_line(0)

    def __init_line(self, line_nr):
        self._line_nr = line_nr
        self._line_indents = [0]
        self._level = 0
        self._br_indent_list = [0]

    def process_lines_of_fline(self, f_line, lines, rel_ind, line_nr):
        """
        process all lines that belong to a Fortran line `f_line`, `rel_ind` is the relative indentation size.
        """

        self.__init_line(line_nr)

        is_decl = typeRe.match(f_line) or PUBLIC_RE.match(f_line)
        for pos, line in enumerate(lines):
            self.__align_line_continuations(
                line, is_decl, rel_ind, self._line_nr + pos)
            if pos + 1 < len(lines):
                self._line_indents.append(self._br_indent_list[-1])

        if len(self._br_indent_list) > 2 or self._level:
            raise SyntaxError(self._filename + ':' + str(self._line_nr) +
                              ':' + FORTRAN_DEFAULT_ERROR_MESSAGE)

    def get_lines_indent(self):
        """
        after processing, retrieve the indents of all line parts.
        """
        return self._line_indents

    def __align_line_continuations(self, line, is_decl, indent_size, line_nr):

        indent_list = self._br_indent_list
        level = self._level
        filename = self._filename

        pos_eq = 0
        pos_ldelim = []
        pos_rdelim = []
        ldelim = []
        rdelim = []

        # find delimiters that are not ended on this line.
        # find proper alignment to most inner delimiter
        # or alignment to assignment operator
        rel_ind = indent_list[-1]  # indentation of prev. line

        instring = ''
        end_of_delim = -1

        for pos, char in CharFilter(enumerate(line)):

            what_del_open = None
            what_del_close = None
            if not pos == end_of_delim:
                what_del_open = DEL_OPEN_RE.match(line[pos:pos + 2])
                what_del_close = DEL_CLOSE_RE.match(line[pos:pos + 2])

            if not instring and what_del_open:
                what_del_open = what_del_open.group()
                end_of_delim = pos + len(what_del_open) - 1
                level += 1
                indent_list.append(pos + len(what_del_open) + rel_ind)
                pos_ldelim.append(pos)
                ldelim.append(what_del_open)
            if not instring and what_del_close:
                what_del_close = what_del_close.group()
                end_of_delim = pos + len(what_del_close) - 1
                level += -1
                indent_list.pop()
                if level < 0:
                    raise SyntaxError(filename + ':' + str(line_nr) +
                                      ':' + FORTRAN_DEFAULT_ERROR_MESSAGE)
                if pos_ldelim:
                    pos_ldelim.pop()
                    what_del_open = ldelim.pop()
                    valid = False
                    if what_del_open == r"(":
                        valid = what_del_close == r")"
                    if what_del_open == r"(/":
                        valid = what_del_close == r"/)"
                    if what_del_open == r"[":
                        valid = what_del_close == r"]"
                    if not valid:
                        raise SyntaxError(
                            filename + ':' + str(line_nr) + ':' + FORTRAN_DEFAULT_ERROR_MESSAGE)
                else:
                    pos_rdelim.append(pos)
                    rdelim.append(what_del_close)
            if not instring and not level:
                if not is_decl and char == '=':
                    if not REL_OP_RE.match(line[max(0, pos - 1):min(pos + 2, len(line))]):
                        if pos_eq > 0:
                            raise SyntaxError(
                                filename + ':' + str(line_nr) + ':' + FORTRAN_DEFAULT_ERROR_MESSAGE)
                        is_pointer = line[pos + 1] == '>'
                        pos_eq = pos + 1
                        # don't align if assignment operator directly before
                        # line break
                        if not re.search(r"=>?\s*" + LINEBREAK_STR, line, RE_FLAGS):
                            indent_list.append(
                                pos_eq + 1 + is_pointer + indent_list[-1])
                elif is_decl and line[pos:pos + 2] == '::':
                    if not re.search(r"::\s*" + LINEBREAK_STR, line, RE_FLAGS):
                        indent_list.append(pos + 3 + indent_list[-1])

        # Don't align if delimiter opening directly before line break
        if level and re.search(DEL_OPEN_STR + r"\s*" + LINEBREAK_STR, line, RE_FLAGS):
            if len(indent_list) > 1:
                indent_list[-1] = indent_list[-2]
            else:
                indent_list[-1] = 0

        if not indent_list[-1]:
            indent_list[-1] = indent_size

        self._level = level

#=========================================================================


def inspect_ffile_format(infile, indent_size):
    """
    Determine indentation by inspecting original Fortran file (mainly for finding
    aligned blocks of DO/IF statements). Also check if
    it has f77 constructs.
    """

    adopt = indent_size <= 0

    is_f90 = True
    indents = []
    stream = InputStream(infile)
    prev_offset = 0
    while 1:
        f_line, _, lines = stream.nextFortranLine()
        if not lines:
            break

        offset = len(lines[0]) - len(lines[0].lstrip(' '))
        indents.append(offset - prev_offset)
        if not adopt:  # do not adopt indentations but impose fixed rel. ind.
            # but don't impose indentation for blocked do/if constructs
            if prev_offset != offset or (not IF_RE.search(f_line) and not DO_RE.search(f_line)):
                indents[-1] = indent_size
        prev_offset = offset

        # can not handle f77 style constructs
        if F77_STYLE.search(f_line):
            is_f90 = False
    return indents, is_f90

#=========================================================================


def format_single_fline(f_line, whitespace, linebreak_pos, ampersand_sep, filename, line_nr, auto_format=True):
    """
    format a single Fortran line - imposes white space formatting
    and inserts linebreaks.
    Takes a logical Fortran line `f_line` as input as well as the positions 
    of the linebreaks (`linebreak_pos`), and the number of separating whitespace 
    characters before ampersand (`ampersand_sep`). 
    `filename` and `line_nr` just for error messages.
    The higher `whitespace`, the more white space characters inserted - 
    whitespace = 0, 1, 2 are currently supported.
    auto formatting can be turned off by setting `auto_format` to False.
    """

    # define whether to put whitespaces around operators:
    # 0: comma, semicolon
    # 1: assignment operators
    # 2: relational operators
    # 3: logical operators
    # 4: arithm. operators plus and minus
    if whitespace == 0:
        spacey = [0, 0, 0, 0, 0]
    elif whitespace == 1:
        spacey = [1, 1, 1, 1, 0]
    elif whitespace == 2:
        spacey = [1, 1, 1, 1, 1]
    else:
        raise NotImplementedError("unknown value for whitespace")

    line = f_line
    line_orig = line

    # rm extraneous whitespace chars, except for declarations
    line_ftd = ''
    pos_prev = -1
    for pos, char in CharFilter(enumerate(line)):
        is_decl = line[pos:].lstrip().startswith('::') or line[
            :pos].rstrip().endswith('::')
        if char == ' ':
            # remove double spaces
            if line_ftd and (re.search(r'[\w"]', line_ftd[-1]) or is_decl):
                line_ftd = line_ftd + char
        else:
            if line_ftd and line_ftd[-1] == ' ' and (not re.search(r'[\w"]', char) and not is_decl):
                line_ftd = line_ftd[:-1]  # remove spaces except between words
            line_ftd = line_ftd + line[pos_prev + 1:pos + 1]
        pos_prev = pos
    line = line_ftd

    pos_eq = []
    end_of_delim = -1
    level = 0
    for pos, char in CharFilter(enumerate(line)):
        # offset w.r.t. unformatted line
        offset = len(line_ftd) - len(line)

        # format delimiters
        what_del_open = None
        what_del_close = None
        if pos > end_of_delim:
            what_del_open = DEL_OPEN_RE.match(
                line[pos:pos + 2])  # opening delimiter token
            what_del_close = DEL_CLOSE_RE.match(
                line[pos:pos + 2])  # closing delimiter token

        if what_del_open or what_del_close:
            sep1 = 0
            sep2 = 0

            if what_del_open:
                delim = what_del_open.group()
            else:
                delim = what_del_close.group()

            lhs = line_ftd[:pos + offset]
            rhs = line_ftd[pos + len(delim) + offset:]

            # format opening delimiters
            if what_del_open:
                level += 1  # new scope
                # add separating whitespace before opening delimiter
                # with some exceptions:
                if (not re.search(r"(" + DEL_OPEN_STR + r"|[\w\*/=\+\-:])\s*$", line[:pos], RE_FLAGS) and
                        not EMPTY_RE.search(line[:pos])) or \
                        re.search(SOL_STR + r"(\w+\s*:)?(ELSE)?\s*IF\s*$", line[:pos], RE_FLAGS) or \
                        re.search(SOL_STR + r"(\w+\s*:)?\s*DO\s+WHILE\s*$", line[:pos], RE_FLAGS) or \
                        re.search(SOL_STR + r"(SELECT)?\s*CASE\s*", line[:pos], RE_FLAGS) or \
                        re.search(r"\b" + INTR_STMTS_PAR + r"\s*$", line[:pos], RE_FLAGS):
                    sep1 = 1

            # format closing delimiters
            else:
                level += -1  # close scope
                # add separating whitespace after closing delimiter
                # with some exceptions:
                if not re.search(r"^\s*(" + DEL_CLOSE_STR + r"|[,%:/\*])", line[pos + 1:], RE_FLAGS):
                    sep2 = 1
                elif re.search(r"^\s*::", line[pos + 1:], RE_FLAGS):
                    sep2 = len(rhs) - len(rhs.lstrip(' '))

            # where delimiter token ends
            end_of_delim = pos + len(delim) - 1

            line_ftd = lhs.rstrip(' ') + ' ' * sep1 + \
                delim + ' ' * sep2 + rhs.lstrip(' ')

        # format commas and semicolons
        if char == ',' or char == ';':
            lhs = line_ftd[:pos + offset]
            rhs = line_ftd[pos + 1 + offset:]
            line_ftd = lhs.rstrip(' ') + char + ' ' * \
                spacey[0] + rhs.lstrip(' ')
            line_ftd = line_ftd.rstrip(' ')

        # format .NOT.
        if re.match(r"\.NOT\.", line[pos:pos + 5], RE_FLAGS):
            lhs = line_ftd[:pos + offset]
            rhs = line_ftd[pos + 5 + offset:]
            line_ftd = lhs.rstrip(
                ' ') + line[pos:pos + 5] + ' ' * spacey[3] + rhs.lstrip(' ')

        # strip whitespaces from '=' and prepare assignment operator
        # formatting
        if char == '=':
            if not REL_OP_RE.search(line[pos - 1:pos + 2]):
                lhs = line_ftd[:pos + offset]
                rhs = line_ftd[pos + 1 + offset:]
                line_ftd = lhs.rstrip(' ') + '=' + rhs.lstrip(' ')
                if not level:  # remember position of assignment operator
                    pos_eq.append(len(lhs.rstrip(' ')))

    line = line_ftd

    # format assignments
    for pos in pos_eq:
        offset = len(line_ftd) - len(line)
        is_pointer = line[pos + 1] == '>'
        lhs = line_ftd[:pos + offset]
        rhs = line_ftd[pos + 1 + is_pointer + offset:]
        if is_pointer:
            assign_op = '=>'  # pointer assignment
        else:
            assign_op = '='  # assignment
        line_ftd = lhs.rstrip(
            ' ') + ' ' * spacey[1] + assign_op + ' ' * spacey[1] + rhs.lstrip(' ')
        # offset w.r.t. unformatted line

    line = line_ftd

    # for more advanced replacements we separate comments and strings
    line_parts = []
    str_end = -1
    instring = ''
    for pos, char in enumerate(line):
        if char == '"' or char == "'":  # skip string
            if not instring:
                str_start = pos
                line_parts.append(line[str_end + 1:str_start])
                instring = char
            elif instring == char:
                str_end = pos
                line_parts.append(line[str_start:str_end + 1])
                instring = ''
        if pos == len(line) - 1:
            line_parts.append(line[str_end + 1:])

    # Two-sided operators
    for n_op, lr_re in enumerate(LR_OPS_RE):
        for pos, part in enumerate(line_parts):
            if not re.match(r"['\"!]", part, RE_FLAGS):  # exclude comments, strings
                partsplit = lr_re.split(part)
                line_parts[pos] = (' ' * spacey[n_op + 2]).join(partsplit)

    line = ''.join(line_parts)

    # format ':' for labels
    for newre in NEW_SCOPE_RE[0:2]:
        if newre.search(line) and re.search(SOL_STR + r"\w+\s*:", line):
            line = ': '.join(_.strip() for _ in line.split(':', 1))

    if not auto_format:
        line = line_orig

    # Now it gets messy - we need to shift line break positions from original
    # to reformatted line
    pos_new = 0
    pos_old = 0
    linebreak_pos.sort(reverse=True)
    linebreak_pos_ftd = []
    while 1:
        if pos_new == len(line) or pos_old == len(line_orig):
            break
        assert line[pos_new] is line_orig[pos_old]
        if linebreak_pos and pos_old > linebreak_pos[-1]:
            linebreak_pos.pop()
            linebreak_pos_ftd.append(pos_new)
            continue
        if line[pos_new] is line_orig[pos_old]:
            pos_new += 1
            while pos_new < len(line) and line[pos_new] is ' ':
                pos_new += 1
            pos_old += 1
            while pos_old < len(line_orig) and line_orig[pos_old] is ' ':
                pos_old += 1
        elif line[pos_new] is ' ':
            pos_new += 1
        elif line_orig[pos_old] is ' ':
            pos_old += 1
        else:
            assert False

    linebreak_pos_ftd.insert(0, 0)

    # We do not insert ampersands in empty lines and comments lines
    lines_out = [line[l:r].rstrip(' ') + ' ' * ampersand_sep[pos] + '&' * min(1, r - l)
                 for pos, (l, r) in enumerate(zip(linebreak_pos_ftd[0:-1], linebreak_pos_ftd[1:]))]

    lines_out.append(line[linebreak_pos_ftd[-1]:])

    if level != 0:
        raise SyntaxError(filename + ':' + str(line_nr) +
                          ':' + FORTRAN_DEFAULT_ERROR_MESSAGE)

    return lines_out

#=========================================================================


def reformat_ffile(infile, outfile, logFile=sys.stdout, indent_size=2, whitespace=2, orig_filename=None):
    """
    main method to be invoked for formatting a Fortran file.
    """
    debug = False

    # don't change original indentation if rel-indents set to 0
    adopt_indents = indent_size <= 0

    if not orig_filename:
        orig_filename = infile.name

    indenter = F90Indenter(orig_filename)

    req_indents, is_f90 = inspect_ffile_format(infile, indent_size)
    infile.seek(0)

    if not is_f90:
        logFile.write("*** " + orig_filename +
                      ": formatter can not handle f77 constructs. ***\n")
        outfile.write(infile.read())  # does not handle f77 constructs
        return

    nfl = 0  # fortran line counter

    do_indent = True
    stream = InputStream(infile)
    skip_blank = False
    in_manual_block = False

    while 1:
        f_line, comments, lines = stream.nextFortranLine()
        if not lines:
            break

        comment_lines = []
        for line, comment in zip(lines, comments):
            has_comment = bool(comment.strip())
            sep = has_comment and not comment.strip() == line.strip()
            if line.strip():  # empty lines between linebreaks are ignored
                comment_lines.append(' ' * sep + comment.strip())

        orig_lines = lines
        nfl += 1

        auto_align = not any(NO_ALIGN_RE.search(_) for _ in lines)
        auto_format = not (in_manual_block or any(
            _.lstrip().startswith('!&') for _ in comment_lines))
        if not auto_format:
            auto_align = False
        if (len(lines)) == 1:
            valid_directive = True
            if lines[0].strip().startswith('!&<'):
                if in_manual_block:
                    valid_directive = False
                else:
                    in_manual_block = True
            if lines[0].strip().startswith('!&>'):
                if not in_manual_block:
                    valid_directive = False
                else:
                    in_manual_block = False
            if not valid_directive:
                raise SyntaxError(orig_filename + ':' + str(stream.line_nr) +
                                  ':' + FORMATTER_ERROR_MESSAGE)

        indent = [0] * len(lines)

        is_omp_conditional = False

        if ompRe.match(f_line) and not ompDirRe.match(f_line):
            # convert OMP-conditional fortran statements into normal fortran statements
            # but remember to convert them back
            f_line = ompRe.sub('  ', f_line, count=1)
            lines = [ompRe.sub('  ', l, count=1) for l in lines]
            is_omp_conditional = True

        is_empty = EMPTY_RE.search(f_line)  # blank line or comment only line

        if useParseRe.match(f_line):
            pass  # do not touch use statements cleaned up by existing prettify
        elif ompDirRe.match(f_line):
            # move '!$OMP' to line start, otherwise don't format omp directives
            lines = ['!$OMP' + (len(l) - len(l.lstrip())) *
                     ' ' + ompDirRe.sub('', l, count=1) for l in lines]
            pass
        elif lines[0].startswith('#'):  # preprocessor macros
            assert len(lines) == 1
        elif EMPTY_RE.search(f_line):  # empty lines including comment lines
            assert len(lines) == 1
            if any(comments):
                if lines[0].startswith('!'):
                    pass  # don't indent unindented comment lines
                else:
                    indent[0] = indenter.get_fline_indent()
            elif skip_blank:
                continue

            lines = [l.strip(' ') for l in lines]
        else:

            manual_lines_indent = []
            if not auto_align:
                manual_lines_indent = [
                    len(l) - len(l.lstrip(' ').lstrip('&')) for l in lines]
                manual_lines_indent = [ind - manual_lines_indent[0]
                                       for ind in manual_lines_indent]

            # ampersands at line starts are remembered (pre_ampersand) and recovered later;
            # define the desired number of separating whitespaces before ampersand at line end (ampersand_sep):
            # - insert one whitespace character before ampersand as default formatting
            # - don't do this if next line starts with an ampersand but remember the original formatting
            # this "special rule" is necessary since ampersands starting a line can be used to break literals,
            # so inserting a whitespace in this case leads to invalid syntax.

            pre_ampersand = []
            ampersand_sep = []
            sep_next = None
            for pos, line in enumerate(lines):
                m = re.search(SOL_STR + r'(&\s*)', line)
                if m:
                    pre_ampersand.append(m.group(1))
                    sep = len(
                        re.search(r'(\s*)&[\s]*(?:!.*)?$', lines[pos - 1]).group(1))
                    ampersand_sep.append(sep)
                else:
                    pre_ampersand.append('')
                    if pos > 0:
                        ampersand_sep.append(1)

            lines = [l.strip(' ').strip('&') for l in lines]
            f_line = f_line.strip(' ')

            # find linebreak positions
            linebreak_pos = []
            for pos, line in enumerate(lines):
                found = None
                for char_pos, char in CharFilter(enumerate(line)):
                    if char == "&":
                        found = char_pos
                if found:
                    linebreak_pos.append(found)
                elif line.lstrip(' ').startswith('!'):
                    linebreak_pos.append(0)

            linebreak_pos = [sum(linebreak_pos[0:_ + 1]) -
                             1 for _ in range(0, len(linebreak_pos))]

            lines = format_single_fline(
                f_line, whitespace, linebreak_pos, ampersand_sep, orig_filename, stream.line_nr, auto_format)

            # we need to insert comments in formatted lines
            for pos, (line, comment) in enumerate(zip(lines, comment_lines)):
                if pos < len(lines) - 1:
                    has_nl = True
                else:
                    has_nl = not re.search(EOL_SC, line)
                lines[pos] = lines[pos].rstrip(' ') + comment + '\n' * has_nl

            try:
                rel_indent = req_indents[nfl]
            except IndexError:
                rel_indent = 0

            indenter.process_lines_of_fline(
                f_line, lines, rel_indent, indent_size, stream.line_nr, manual_lines_indent)
            indent = indenter.get_lines_indent()

            # recover ampersands at line start
            for pos, line in enumerate(lines):
                amp_insert = pre_ampersand[pos]
                if amp_insert:
                    indent[pos] += -1
                    lines[pos] = amp_insert + line

        lines = [re.sub(r"\s+$", '\n', l, RE_FLAGS)
                 for l in lines]  # deleting trailing whitespaces

        for ind, line, orig_line in zip(indent, lines, orig_lines):
            # get actual line length excluding comment:
            line_length = 0
            for line_length, _ in CharFilter(enumerate(line)):
                pass
            line_length += 1

            if do_indent:
                ind_use = ind
            else:
                ind_use = 1
            if ind_use + line_length <= 133:  # 132 plus 1 newline char
                outfile.write('!$' * is_omp_conditional + ' ' *
                              (ind_use - 2 * is_omp_conditional +
                               len(line) - len(line.lstrip(' '))) + line.lstrip(' '))
            elif line_length <= 133:
                outfile.write('!$' * is_omp_conditional + ' ' *
                              (133 - 2 * is_omp_conditional -
                               len(line.lstrip(' '))) + line.lstrip(' '))
                if not typeRe.match(f_line):
                    logFile.write("*** " + orig_filename + ":" + str(stream.line_nr) +
                                  ": auto indentation failed due to 132 chars limit, line should be splitted. ***\n")
            else:
                outfile.write(orig_line)
                logFile.write("*** " + orig_filename + ":" + str(stream.line_nr) +
                              (": auto indentation and whitespace formatting failed due to 132 chars limit, line should be splitted. ***\n"))
            if debug:
                print(' ' * ind_use + line)
        # no indentation of blank lines
        if re.search(r";\s*$", f_line, RE_FLAGS):
            do_indent = False
        else:
            do_indent = True

        # rm subsequent blank lines
        skip_blank = is_empty and not any(comments)

try:
    any
except NameError:
    def any(iterable):
        for element in iterable:
            if element:
                return True
        return False

# EOF
