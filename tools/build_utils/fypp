#!/usr/bin/env python3
####################################################################################################
#
# fypp -- Python powered Fortran preprocessor
#
# Copyright (c) 2016-2026 Bálint Aradi, Universität Bremen
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions
# and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of
# conditions and the following disclaimer in the documentation and/or other materials provided with
# the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
####################################################################################################

"""For using the functionality of the Fypp preprocessor from within Python, one usually interacts
with the following two classes:

* `Fypp`_: The actual Fypp preprocessor. It returns for a given input the preprocessed output.

* `FyppOptions`_: Contains customizable settings controlling the behaviour of `Fypp`_.
  Alternatively, the function `get_option_parser()`_ can be used to obtain an option parser, which
  can create settings based on command line arguments.

If processing stops prematurely, an instance of one of the following subclasses of `FyppError`_ is
raised:

* FyppFatalError: Unexpected error (e.g. bad input, missing files, etc.)

* FyppStopRequest: Stop was triggered by an explicit request in the input (by a stop- or an assert-
  directive).
"""

from __future__ import annotations

import builtins
import contextlib
import dataclasses
import inspect
import itertools
import optparse
import os
import pathlib
import platform
import re
import sys
import time
import types
from collections.abc import Callable, Generator, Sequence
from typing import Any, ClassVar, NamedTuple, NoReturn, Protocol, TextIO, TypeAlias

MIN_PYTHON_VERSION = (3, 10)
if sys.version_info < MIN_PYTHON_VERSION:
    sys.exit(
        f"Fypp requires Python version {MIN_PYTHON_VERSION[0]}.{MIN_PYTHON_VERSION[1]} or later.\n"
    )

# Prevent cluttering user directory with Python bytecode
sys.dont_write_bytecode = True

# Fypp version
VERSION = '3.3'

# String used as filename if input is read from standard input
STDIN_INPUT_NAME = "<stdin>"

# String used as filename if input is read from a file object (via API)
FOBJ_INPUT_NAME = "<fileobj>"

# String used as filename if input was read from a string (via API)
STRING_INPUT_NAME = "<string>"

# Exit code when error happened
ERROR_EXIT_CODE = 1

# Exit code when user input (stop, assert) explicitly triggered a stop
USER_ERROR_EXIT_CODE = 2

_INLINE_DIRECTIVE_PATTERN = r"(?P<idirtype>[$\#@])\{[ \t]*(?P<idir>.+?)?[ \t]*\}(?P=idirtype)"

_ALL_DIRECTIVES_PATTERN = r"""
# comment block
(?:^[ \t]*\#!.*(?:\n|\Z))+
|
# line directive (with optional continuation lines)
^[ \t]*(?P<ldirtype>[\#\$@]):[ \t]*
(?P<ldir>.+?(?:&[ \t]*\n(?:[ \t]*&)?.*?)*)?[ \t]*(?:\n|\Z)
|
# inline directive
""" + _INLINE_DIRECTIVE_PATTERN

_ALL_DIRECTIVES_REGEXP = re.compile(_ALL_DIRECTIVES_PATTERN, re.VERBOSE | re.MULTILINE)

_CONTROL_DIR_REGEXP = re.compile(r"(?P<dir>[a-zA-Z_]\w*)[ \t]*(?:[ \t]+(?P<param>[^ \t].*))?$")

_DIRECT_CALL_REGEXP = re.compile(r"(?P<callname>[a-zA-Z_][\w.]*)[ \t]*\((?P<callparams>.+?)?\)$")

_DIRECT_CALL_ARG_REGEXP = re.compile(_INLINE_DIRECTIVE_PATTERN)

_DIRECT_CALL_KWARG_REGEXP = re.compile(r"(?:(?P<kwname>[a-zA-Z_]\w*)\s*=(?=[^=]|$))?")

_DEF_PARAM_REGEXP = re.compile(r"^(?P<name>[a-zA-Z_]\w*)[ \t]*\(\s*(?P<args>.+)?\s*\)$")

_SIMPLE_CALLABLE_REGEXP = re.compile(
    r"^(?P<name>[a-zA-Z_][\w.]*)[ \t]*(?:\([ \t]*(?P<args>.*)[ \t]*\))?$"
)

_IDENTIFIER_NAME_REGEXP = re.compile(r"^(?P<name>[a-zA-Z_]\w*)$")

_PREFIXED_IDENTIFIER_NAME_REGEXP = re.compile(r"^(?P<name>[a-zA-Z_][\w.]*)$")

_SET_PARAM_REGEXP = re.compile(
    r"^(?P<name>(?:[(]\s*)?[a-zA-Z_]\w*(?:\s*,\s*[a-zA-Z_]\w*)*(?:\s*[)])?)\s*"
    r"(?:=\s*(?P<expr>.*))?$"
)

_DEL_PARAM_REGEXP = re.compile(r"^(?:[(]\s*)?[a-zA-Z_]\w*(?:\s*,\s*[a-zA-Z_]\w*)*(?:\s*[)])?$")

_FOR_PARAM_REGEXP = re.compile(
    r"^(?P<loopexpr>[a-zA-Z_]\w*(\s*,\s*[a-zA-Z_]\w*)*)\s+in\s+(?P<iter>.+)$"
)

_INCLUDE_PARAM_REGEXP = re.compile(r'^(\'|")(?P<fname>.*?)\1$')

_COMMENTLINE_REGEXP = re.compile(r"^[ \t]*!.*$")

_CONTLINE_REGEXP = re.compile(r"&[ \t]*\n(?:[ \t]*&)?")

_UNESCAPE_TEXT_REGEXP1 = re.compile(r"([$#@])\\(\\*)([{:])")

_UNESCAPE_TEXT_REGEXP2 = re.compile(r"#\\(\\*)([!])")

_UNESCAPE_TEXT_REGEXP3 = re.compile(r"(\})\\(\\*)([$#@])")

_INLINE_EVAL_REGION_REGEXP = re.compile(r"\${.*?}\$")

_RESERVED_PREFIX = "__"

_RESERVED_NAMES = frozenset(
    {
        "defined",
        "setvar",
        "getvar",
        "delvar",
        "globalvar",
        "_LINE_",
        "_FILE_",
        "_THIS_FILE_",
        "_THIS_LINE_",
        "_TIME_",
        "_DATE_",
        "_SYSTEM_",
        "_MACHINE_",
    }
)

_LINENUM_NEW_FILE = 1

_LINENUM_RETURN_TO_FILE = 2

_QUOTES_FORTRAN = "'\""

_OPENING_BRACKETS_FORTRAN = "{(["

_CLOSING_BRACKETS_FORTRAN = "})]"

_ARGUMENT_SPLIT_CHAR_FORTRAN = ","


class Span(NamedTuple):
    """Beginning and end line of a region in a source file.

    Line numbers start from zero. For directives, which do not consume the end of the line, start
    and end are identical.
    """

    start: int
    end: int


class FyppError(Exception):
    """Signalizes error occurring during preprocessing.

    Args:
        msg: Error message.
        fname: File name. None (default) if file name is not available.
        span: Beginning and end line of the region where error occurred or None if not available. If
            fname was not None, span must not be None.

    Attributes:
        msg: Error message.
        fname: File name or None if not available.
        span: Beginning and end line of the region where error occurred or None if not available.
            Line numbers start from zero. For directives, which do not consume end of the line,
            start and end lines are identical.
    """

    def __init__(self, msg: str, fname: str | None = None, span: Span | None = None):
        assert fname is None or span is not None
        super().__init__()
        self.msg: str = msg
        self.fname: str | None = fname
        self.span: Span | None = span

    def __str__(self) -> str:
        msg = [self.__class__.__name__, ": "]
        if self.fname is not None:
            # Guaranteed by the class contract: fname is only ever set together with span.
            assert self.span is not None
            msg.append("file '" + self.fname + "'")
            if self.span.end > self.span.start + 1:
                msg.append(f", lines {self.span.start + 1}-{self.span.end}")
            else:
                msg.append(f", line {self.span.start + 1}")
            msg.append("\n")
        if self.msg:
            msg.append(self.msg)
        if self.__cause__ is not None:
            msg.append("\n" + str(self.__cause__))
        return "".join(msg)


class FyppFatalError(FyppError):
    """Signalizes an unexpected error during processing."""


class FyppStopRequest(FyppError):
    """Signalizes an explicitly triggered stop (e.g. via stop directive)"""


class _ParserHandlers(Protocol):
    """Structural shape of the callback object receiving events from a Parser.

    Satisfied by HandlerLogger (Parser's default handlers, printing each event) and by Builder
    (turning parser events into a tree representation).
    """

    def handle_include(self, span: Span | None, fname: str) -> None: ...
    def handle_endinclude(self, span: Span | None, fname: str) -> None: ...
    def handle_set(self, span: Span, name: str, expr: str) -> None: ...
    def handle_def(self, span: Span, name: str, argexpr: str | None) -> None: ...
    def handle_enddef(self, span: Span, name: str | None) -> None: ...
    def handle_del(self, span: Span, name: str) -> None: ...
    def handle_if(self, span: Span, cond: str) -> None: ...
    def handle_elif(self, span: Span, cond: str) -> None: ...
    def handle_else(self, span: Span) -> None: ...
    def handle_endif(self, span: Span) -> None: ...
    def handle_for(self, span: Span, loopvar: list[str], iterator: str) -> None: ...
    def handle_endfor(self, span: Span) -> None: ...
    def handle_call(self, span: Span, name: str, argexpr: str | None, blockcall: bool) -> None: ...
    def handle_nextarg(self, span: Span, name: str | None, blockcall: bool) -> None: ...
    def handle_endcall(self, span: Span, name: str | None, blockcall: bool) -> None: ...
    def handle_eval(self, span: Span, expr: str) -> None: ...
    def handle_global(self, span: Span, name: str) -> None: ...
    def handle_text(self, span: Span, txt: str) -> None: ...
    def handle_comment(self, span: Span) -> None: ...
    def handle_mute(self, span: Span) -> None: ...
    def handle_endmute(self, span: Span) -> None: ...
    def handle_stop(self, span: Span, msg: str) -> None: ...
    def handle_assert(self, span: Span, cond: str) -> None: ...


class HandlerLogger:
    """Handlers implementation which logs each event by printing it.

    Used as Parser's default handlers object, so that a Parser can be run stand-alone (without a
    Builder) to inspect the event stream it generates, e.g. for debugging.
    """

    def handle_include(self, span: Span | None, fname: str) -> None:
        """Called when parser starts to process a new file.

        Args:
            span: Start and end line of the include directive or None if called the first time for
                the main input.
            fname: Name of the file.
        """
        self._log_event("include", span, filename=fname)

    def handle_endinclude(self, span: Span | None, fname: str) -> None:
        """Called when parser finished processing a file.

        Args:
            span: Start and end line of the include directive or None if called the first time for
                the main input.
            fname: Name of the file.
        """
        self._log_event("endinclude", span, filename=fname)

    def handle_set(self, span: Span, name: str, expr: str) -> None:
        """Called when parser encounters a set directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the variable.
            expr: String representation of the expression to be assigned to the variable.
        """
        self._log_event("set", span, name=name, expression=expr)

    def handle_def(self, span: Span, name: str, argexpr: str | None) -> None:
        """Called when parser encounters a def directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the macro to be defined.
            argexpr: String with argument definition (or None)
        """
        self._log_event("def", span, name=name, arguments=argexpr)

    def handle_enddef(self, span: Span, name: str | None) -> None:
        """Called when parser encounters an enddef directive.

        Args:
            span: Start and end line of the directive.
            name: Name found after the enddef directive.
        """
        self._log_event("enddef", span, name=name)

    def handle_del(self, span: Span, name: str) -> None:
        """Called when parser encounters a del directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the variable to delete.
        """
        self._log_event("del", span, name=name)

    def handle_if(self, span: Span, cond: str) -> None:
        """Called when parser encounters an if directive.

        Args:
            span: Start and end line of the directive.
            cond: String representation of the branching condition.
        """
        self._log_event("if", span, condition=cond)

    def handle_elif(self, span: Span, cond: str) -> None:
        """Called when parser encounters an elif directive.

        Args:
            span: Start and end line of the directive.
            cond: String representation of the branching condition.
        """
        self._log_event("elif", span, condition=cond)

    def handle_else(self, span: Span) -> None:
        """Called when parser encounters an else directive.

        Args:
            span: Start and end line of the directive.
        """
        self._log_event("else", span)

    def handle_endif(self, span: Span) -> None:
        """Called when parser encounters an endif directive.

        Args:
            span: Start and end line of the directive.
        """
        self._log_event("endif", span)

    def handle_for(self, span: Span, loopvar: list[str], iterator: str) -> None:
        """Called when parser encounters a for directive.

        Args:
            span: Start and end line of the directive.
            loopvar: Names of the loop variable(s).
            iterator: String representation of the iterable.
        """
        self._log_event("for", span, variable=loopvar, iterable=iterator)

    def handle_endfor(self, span: Span) -> None:
        """Called when parser encounters an endfor directive.

        Args:
            span: Start and end line of the directive.
        """
        self._log_event("endfor", span)

    def handle_call(self, span: Span, name: str, argexpr: str | None, blockcall: bool) -> None:
        """Called when parser encounters a call directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the callable to call
            argexpr: Argument expression containing additional arguments for the call.
            blockcall: Whether the alternative "block / contains / endblock" calling directive has
                been used.
        """
        self._log_event("call", span, name=name, argexpr=argexpr, blockcall=blockcall)

    def handle_nextarg(self, span: Span, name: str | None, blockcall: bool) -> None:
        """Called when parser encounters a nextarg directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the argument following next or None if it should be the next positional
                argument.
            blockcall: Whether the alternative "block / contains / endblock" calling directive has
                been used.
        """
        self._log_event("nextarg", span, name=name, blockcall=blockcall)

    def handle_endcall(self, span: Span, name: str | None, blockcall: bool) -> None:
        """Called when parser encounters an endcall directive.

        Args:
            span: Start and end line of the directive.
            name: Name found after the endcall directive.
            blockcall: Whether the alternative "block / contains / endblock" calling directive has
                been used.
        """
        self._log_event("endcall", span, name=name, blockcall=blockcall)

    def handle_eval(self, span: Span, expr: str) -> None:
        """Called when parser encounters an eval directive.

        Args:
            span: Start and end line of the directive.
            expr: String representation of the Python expression to be evaluated.
        """
        self._log_event("eval", span, expression=expr)

    def handle_global(self, span: Span, name: str) -> None:
        """Called when parser encounters a global directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the variable which should be made global.
        """
        self._log_event("global", span, name=name)

    def handle_text(self, span: Span, txt: str) -> None:
        """Called when parser finds text which must left unaltered.

        Args:
            span: Start and end line of the directive.
            txt: Text.
        """
        self._log_event("text", span, content=txt)

    def handle_comment(self, span: Span) -> None:
        """Called when parser finds a preprocessor comment.

        Args:
            span: Start and end line of the directive.
        """
        self._log_event("comment", span)

    def handle_mute(self, span: Span) -> None:
        """Called when parser finds a mute directive.

        Args:
            span: Start and end line of the directive.
        """
        self._log_event("mute", span)

    def handle_endmute(self, span: Span) -> None:
        """Called when parser finds an endmute directive.

        Args:
            span: Start and end line of the directive.
        """
        self._log_event("endmute", span)

    def handle_stop(self, span: Span, msg: str) -> None:
        """Called when parser finds an stop directive.

        Args:
            span: Start and end line of the directive.
            msg: Stop message.
        """
        self._log_event("stop", span, msg=msg)

    def handle_assert(self, span: Span, cond: str) -> None:
        """Called when parser finds an assert directive.

        Args:
            span: Start and end line of the directive.
            cond: String representation of the condition for the assert directive.
        """
        self._log_event("assert", span, condition=cond)

    @staticmethod
    def _log_event(event: str, span: Span | None, **params: Any) -> None:
        start, end = (-1, -1) if span is None else span
        print(f"{event}: {start} --> {end}")
        for parname, parvalue in params.items():
            print(f"  {parname}: ->|{parvalue}|<-")
        print()


class Parser:
    """Parses a text and generates events when encountering Fypp constructs.

    Args:
        includedirs: List of directories, in which include files should be searched for, when they
            are not found at the default location.

        encoding: Encoding to use when reading the file (default: utf-8)

        handlers: Object whose handle_* methods are invoked for the corresponding Fypp constructs.
            If None (default), a HandlerLogger instance is used, which prints each event (useful
            for inspecting the event stream stand-alone, without a Builder).

    Attributes:
        handlers: Object whose handle_* methods are invoked for the corresponding Fypp constructs.
            May be reassigned later (e.g. to a Builder instance) to dispatch events elsewhere.
    """

    def __init__(
        self,
        includedirs: list[str] | None = None,
        encoding: str = "utf-8",
        handlers: _ParserHandlers | None = None,
    ):

        # Directories to search for include files
        if includedirs is None:
            self._includedirs: list[str] = []
        else:
            self._includedirs = includedirs

        # Encoding
        self._encoding: str = encoding

        # Name of current file
        self._file: str | None = None

        # Directory of current file
        self._curdir: str | None = None

        # List of files included via include directive
        self._included_files: list[str] = []

        # Handlers to be called for parsing events
        self.handlers: _ParserHandlers = HandlerLogger() if handlers is None else handlers

    def parsefile(self, fobj: str | pathlib.Path | TextIO) -> None:
        """Parses file or a file like object.

        Args:
            fobj: Name of a file or a file like object.
        """
        self._included_files = []
        if isinstance(fobj, pathlib.Path):
            fobj = str(fobj)
        if isinstance(fobj, str):
            if fobj == STDIN_INPUT_NAME:
                self._includefile(None, sys.stdin, STDIN_INPUT_NAME, os.getcwd())
            else:
                with _open_input_file(fobj, self._encoding) as inpfp:
                    self._includefile(None, inpfp, fobj, os.path.dirname(fobj))
        else:
            self._includefile(None, fobj, FOBJ_INPUT_NAME, os.getcwd())

    def parse(self, txt: str) -> None:
        """Parses string.

        Args:
            txt: Text to parse.
        """
        self._file = STRING_INPUT_NAME
        self._curdir = ""
        self._included_files = []
        self._parse_txt(None, self._file, txt)

    def get_included_files(self) -> list[str]:
        """Returns the list of the files included during parsing"""
        return list(self._included_files)

    def _includefile(self, span: Span | None, fobj: TextIO, fname: str, curdir: str) -> None:
        oldfile = self._file
        olddir = self._curdir
        self._file = fname
        self._curdir = curdir
        try:
            self._parse_txt(span, fname, fobj.read())
        finally:
            self._file = oldfile
            self._curdir = olddir

    def _parse_txt(self, includespan: Span | None, fname: str, txt: str) -> None:
        self.handlers.handle_include(includespan, fname)
        self._parse(txt)
        self.handlers.handle_endinclude(includespan, fname)

    def _parse(self, txt: str, linenr: int = 0, directcall: bool = False) -> None:
        pos = 0
        if directcall:
            matches = _DIRECT_CALL_ARG_REGEXP.finditer(txt)
        else:
            matches = _ALL_DIRECTIVES_REGEXP.finditer(txt)
        for match in matches:
            groups = match.groupdict()
            idirtype, idir = groups.get("idirtype"), groups.get("idir")
            ldirtype, ldir = groups.get("ldirtype"), groups.get("ldir")
            start, end = match.span()
            if start > pos:
                endlinenr = linenr + txt.count("\n", pos, start)
                self._process_text(txt[pos:start], Span(linenr, endlinenr))
                linenr = endlinenr
            endlinenr = linenr + txt.count("\n", start, end)
            if idirtype is None and not txt.endswith("\n", start, end):
                # Line directive or comment terminated by the end of the text instead of a
                # newline: treat it as if the missing newline had been present.
                endlinenr += 1
            span = Span(linenr, endlinenr)
            if directcall and idirtype != "$":
                msg = "only inline eval directives allowed in direct calls"
                raise FyppFatalError(msg, self._file, span)
            elif idirtype is not None:
                if idir is None:
                    msg = "missing inline directive content"
                    raise FyppFatalError(msg, self._file, span)
                dirtype = idirtype
                content = idir
            elif ldirtype is not None:
                if ldir is None:
                    msg = "missing line directive content"
                    raise FyppFatalError(msg, self._file, span)
                dirtype = ldirtype
                content = _CONTLINE_REGEXP.sub("", ldir)
            else:
                # Comment directive
                dirtype = ""
                content = ""
            match dirtype:
                case "$":
                    self.handlers.handle_eval(span, content)
                case "#":
                    self._process_control_dir(content, span)
                case "@":
                    self._process_direct_call(content, span)
                case _:
                    self.handlers.handle_comment(span)
            pos = end
            linenr = endlinenr
        if pos < len(txt):
            endlinenr = linenr + txt.count("\n", pos)
            self._process_text(txt[pos:], Span(linenr, endlinenr))

    def _process_text(self, txt: str, span: Span) -> None:
        unescaped_txt = self._unescape(txt)
        self.handlers.handle_text(span, unescaped_txt)

    def _process_control_dir(self, content: str, span: Span) -> None:
        match = _CONTROL_DIR_REGEXP.match(content)
        if not match:
            msg = f"invalid control directive content '{content}'"
            raise FyppFatalError(msg, self._file, span)
        directive, param = match.groups()
        match directive:
            case "if":
                self._check_param_presence(True, "if", param, span)
                self.handlers.handle_if(span, param)
            case "else":
                self._check_param_presence(False, "else", param, span)
                self.handlers.handle_else(span)
            case "elif":
                self._check_param_presence(True, "elif", param, span)
                self.handlers.handle_elif(span, param)
            case "endif":
                self._check_param_presence(False, "endif", param, span)
                self.handlers.handle_endif(span)
            case "def":
                self._check_param_presence(True, "def", param, span)
                self._check_not_inline_directive("def", span)
                self._process_def(param, span)
            case "enddef":
                self._process_enddef(param, span)
            case "set":
                self._check_param_presence(True, "set", param, span)
                self._process_set(param, span)
            case "del":
                self._check_param_presence(True, "del", param, span)
                self._process_del(param, span)
            case "for":
                self._check_param_presence(True, "for", param, span)
                self._process_for(param, span)
            case "endfor":
                self._check_param_presence(False, "endfor", param, span)
                self.handlers.handle_endfor(span)
            case "call" | "block":
                self._check_param_presence(True, directive, param, span)
                self._process_call(param, span, directive == "block")
            case "nextarg" | "contains":
                self._process_nextarg(param, span, directive == "contains")
            case "endcall" | "endblock":
                self._process_endcall(param, span, directive == "endblock")
            case "include":
                self._check_param_presence(True, "include", param, span)
                self._check_not_inline_directive("include", span)
                self._process_include(param, span)
            case "mute":
                self._check_param_presence(False, "mute", param, span)
                self._check_not_inline_directive("mute", span)
                self._process_mute(span)
            case "endmute":
                self._check_param_presence(False, "endmute", param, span)
                self._check_not_inline_directive("endmute", span)
                self._process_endmute(span)
            case "stop":
                self._check_param_presence(True, "stop", param, span)
                self._check_not_inline_directive("stop", span)
                self.handlers.handle_stop(span, param)
            case "assert":
                self._check_param_presence(True, "assert", param, span)
                self._check_not_inline_directive("assert", span)
                self.handlers.handle_assert(span, param)
            case "global":
                self._check_param_presence(True, "global", param, span)
                self._process_global(param, span)
            case _:
                msg = f"unknown directive '{directive}'"
                raise FyppFatalError(msg, self._file, span)

    def _process_direct_call(self, callexpr: str, span: Span) -> None:
        match = _DIRECT_CALL_REGEXP.match(callexpr)
        if not match:
            msg = "invalid direct call expression"
            raise FyppFatalError(msg, self._file, span)
        callname = match.group("callname")
        self.handlers.handle_call(span, callname, None, False)
        callparams = match.group("callparams")
        if callparams is None or not callparams.strip():
            args = []
        else:
            try:
                args = [arg.strip() for arg in _argsplit_fortran(callparams)]
            except Exception as exc:
                msg = "unable to parse direct call argument"
                raise FyppFatalError(msg, self._file, span) from exc
        for arg in args:
            match = _DIRECT_CALL_KWARG_REGEXP.match(arg)
            assert match is not None
            argval = arg[match.end() :].strip()
            # Remove enclosing braces if present
            if argval.startswith("{"):
                argval = argval[1:-1]
            keyword = match.group("kwname")
            self.handlers.handle_nextarg(span, keyword, False)
            self._parse(argval, linenr=span.start, directcall=True)
        self.handlers.handle_endcall(span, callname, False)

    def _process_def(self, param: str, span: Span) -> None:
        match = _DEF_PARAM_REGEXP.match(param)
        if not match:
            msg = f"invalid macro definition '{param}'"
            raise FyppFatalError(msg, self._file, span)
        name = match.group("name")
        argexpr = match.group("args")
        self.handlers.handle_def(span, name, argexpr)

    def _process_enddef(self, param: str | None, span: Span) -> None:
        if param is not None:
            match = _IDENTIFIER_NAME_REGEXP.match(param)
            if not match:
                msg = f"invalid enddef parameter '{param}'"
                raise FyppFatalError(msg, self._file, span)
            param = match.group("name")
        self.handlers.handle_enddef(span, param)

    def _process_set(self, param: str, span: Span) -> None:
        match = _SET_PARAM_REGEXP.match(param)
        if not match:
            msg = f"invalid variable assignment '{param}'"
            raise FyppFatalError(msg, self._file, span)
        self.handlers.handle_set(span, match.group("name"), match.group("expr"))

    def _process_global(self, param: str, span: Span) -> None:
        match = _DEL_PARAM_REGEXP.match(param)
        if not match:
            msg = f"invalid variable specification '{param}'"
            raise FyppFatalError(msg, self._file, span)
        self.handlers.handle_global(span, param)

    def _process_del(self, param: str, span: Span) -> None:
        match = _DEL_PARAM_REGEXP.match(param)
        if not match:
            msg = f"invalid variable specification '{param}'"
            raise FyppFatalError(msg, self._file, span)
        self.handlers.handle_del(span, param)

    def _process_for(self, param: str, span: Span) -> None:
        match = _FOR_PARAM_REGEXP.match(param)
        if not match:
            msg = f"invalid for loop declaration '{param}'"
            raise FyppFatalError(msg, self._file, span)
        loopexpr = match.group("loopexpr")
        loopvars = [s.strip() for s in loopexpr.split(",")]
        self.handlers.handle_for(span, loopvars, match.group("iter"))

    def _process_call(self, param: str, span: Span, blockcall: bool) -> None:
        match = _SIMPLE_CALLABLE_REGEXP.match(param)
        if not match:
            msg = f"invalid callable expression '{param}'"
            raise FyppFatalError(msg, self._file, span)
        name, args = match.groups()
        self.handlers.handle_call(span, name, args, blockcall)

    def _process_nextarg(self, param: str | None, span: Span, blockcall: bool) -> None:
        if param is not None:
            match = _IDENTIFIER_NAME_REGEXP.match(param)
            if not match:
                msg = f"invalid nextarg parameter '{param}'"
                raise FyppFatalError(msg, self._file, span)
            param = match.group("name")
        self.handlers.handle_nextarg(span, param, blockcall)

    def _process_endcall(self, param: str | None, span: Span, blockcall: bool) -> None:
        if param is not None:
            match = _PREFIXED_IDENTIFIER_NAME_REGEXP.match(param)
            if not match:
                msg = f"invalid endcall parameter '{param}'"
                raise FyppFatalError(msg, self._file, span)
            param = match.group("name")
        self.handlers.handle_endcall(span, param, blockcall)

    def _process_include(self, param: str, span: Span) -> None:
        match = _INCLUDE_PARAM_REGEXP.match(param)
        if not match:
            msg = f"invalid include file declaration '{param}'"
            raise FyppFatalError(msg, self._file, span)
        fname = match.group("fname")
        # Always set while a file/string is being parsed, i.e. whenever this routine is called.
        assert self._curdir is not None
        for incdir in [self._curdir] + self._includedirs:
            fpath = os.path.join(incdir, fname)
            if os.path.exists(fpath):
                break
        else:
            msg = f"include file '{fname}' not found"
            raise FyppFatalError(msg, self._file, span)
        self._included_files.append(os.path.normpath(fpath))
        with _open_input_file(fpath, self._encoding) as inpfp:
            self._includefile(span, inpfp, fpath, os.path.dirname(fpath))

    def _process_mute(self, span: Span) -> None:
        if span.start == span.end:
            msg = "Inline form of mute directive not allowed"
            raise FyppFatalError(msg, self._file, span)
        self.handlers.handle_mute(span)

    def _process_endmute(self, span: Span) -> None:
        if span.start == span.end:
            msg = "Inline form of endmute directive not allowed"
            raise FyppFatalError(msg, self._file, span)
        self.handlers.handle_endmute(span)

    def _check_param_presence(
        self, presence: bool, directive: str, param: str | None, span: Span
    ) -> None:
        if (param is not None) != presence:
            if presence:
                msg = f"missing data in {directive} directive"
            else:
                msg = f"forbidden data in {directive} directive"
            raise FyppFatalError(msg, self._file, span)

    def _check_not_inline_directive(self, directive: str, span: Span) -> None:
        if span.start == span.end:
            msg = f"Inline form of {directive} directive not allowed"
            raise FyppFatalError(msg, self._file, span)

    @staticmethod
    def _unescape(txt: str) -> str:
        txt = _UNESCAPE_TEXT_REGEXP1.sub(r"\1\2\3", txt)
        txt = _UNESCAPE_TEXT_REGEXP2.sub(r"#\1\2", txt)
        txt = _UNESCAPE_TEXT_REGEXP3.sub(r"\1\2\3", txt)
        return txt


#
# Data classes used by the builder and the renderer
#


@dataclasses.dataclass(slots=True)
class _RawText:
    """Literal text passed through to the output unaltered."""

    fname: str | None
    span: Span
    text: str


@dataclasses.dataclass(slots=True)
class _EvalDirective:
    """A '$:'/'@:'/inline '${...}$' eval directive."""

    fname: str | None
    span: Span
    expr: str


@dataclasses.dataclass(slots=True)
class _SetDirective:
    """A '#:set' directive."""

    fname: str | None
    span: Span
    name: str
    expr: str


@dataclasses.dataclass(slots=True)
class _DelDirective:
    """A '#:del' directive."""

    fname: str | None
    span: Span
    name: str


@dataclasses.dataclass(slots=True)
class _GlobalDirective:
    """A '#:global' directive."""

    fname: str | None
    span: Span
    name: str


@dataclasses.dataclass(slots=True)
class _CommentDirective:
    """A comment-only (#!) directive line, kept only for line number bookkeeping."""

    fname: str | None
    span: Span


@dataclasses.dataclass(slots=True)
class _StopDirective:
    """A '#:stop' directive."""

    fname: str | None
    span: Span
    msg: str


@dataclasses.dataclass(slots=True)
class _AssertDirective:
    """An '#:assert' directive."""

    fname: str | None
    span: Span
    cond: str


@dataclasses.dataclass(slots=True)
class _IfBlock:
    """An '#:if'/'#:elif'/'#:else'/'#:endif' construct."""

    directive: ClassVar[str] = "if"
    fname: str | None
    spans: list[Span]
    conditions: list[str]
    trees: list[_Tree]


@dataclasses.dataclass(slots=True)
class _ForBlock:
    """A '#:for'/'#:endfor' construct."""

    directive: ClassVar[str] = "for"
    fname: str | None
    spans: list[Span]
    loopvars: list[str]
    loopiter: str
    tree: _Tree | None = None


@dataclasses.dataclass(slots=True)
class _DefBlock:
    """A '#:def'/'#:enddef' construct."""

    directive: ClassVar[str] = "def"
    fname: str | None
    spans: list[Span]
    name: str
    argexpr: str | None
    tree: _Tree | None = None


@dataclasses.dataclass(slots=True)
class _CallBlock:
    """A '#:call'/'#:nextarg'/'#:endcall' (or block/contains/endblock) construct.

    Field `directive` can be either `call` or `block` depending which on the construct in the input.
    """

    directive: str
    fname: str | None
    spans: list[Span]
    name: str
    argexpr: str | None
    trees: list[_Tree] = dataclasses.field(default_factory=list)
    argnames: list[str] = dataclasses.field(default_factory=list)


@dataclasses.dataclass(slots=True)
class _IncludeBlock:
    """Open '#:include'd file being processed."""

    directive: ClassVar[str] = "include"

    fname: str | None
    # spans[0] and spans[1] are both None for the top-level file
    # (include event triggered without a real existing "include" directive)
    spans: list[Span | None]
    includefname: str
    tree: _Tree | None = None


@dataclasses.dataclass(slots=True)
class _MuteBlock:
    """Open '#:mute'/'#:endmute' construct."""

    directive: ClassVar[str] = "mute"

    fname: str | None
    spans: list[Span]
    tree: _Tree | None = None


# Structure of a parsed document:
#
# A tree (`_Tree`) is a list of nodes. Each node (`_Node`) is either a leaf — a single directive or
# raw text (`_SetDirective`, `_EvalDirective, `_RawText`, ...) — or an open/closeable construct
# (`_IfBlock`, `_ForBlock`, `_DefBlock`, `_CallBlock`, `_IncludeBlock, `_MuteBlock`) that holds its
# own content as further tree(s). A block contains either a single tree (for/def/include/mute) or a
# list of trees (one per if/elif/else branch, or one per call argument). Since a blocks body itself
# is a tree, trees nest arbitrarily deep and represent the AST of the input.

# TODO: use `type` once minimal Python version had been bumped to >= 3.12
_Block: TypeAlias = _IfBlock | _ForBlock | _DefBlock | _CallBlock | _IncludeBlock | _MuteBlock
"""A tree node representing an open/closeable Fypp construct, as tracked by Builder while parsing
is still in progress (see Builder._open_blocks)."""

_Node: TypeAlias = (
    _RawText
    | _EvalDirective
    | _SetDirective
    | _DelDirective
    | _GlobalDirective
    | _CommentDirective
    | _StopDirective
    | _AssertDirective
    | _Block
)
"""A single entry of a fypp tree, as produced by Builder and consumed by Renderer."""

_Tree: TypeAlias = list[_Node]
"""A fypp tree: the sequence of nodes making up a piece of source (or a macro/block body)."""


#
# Builder
#


class Builder:
    """Builds a tree representing a text with preprocessor directives."""

    def __init__(self) -> None:
        # The tree representing the text with the preprocessor directives
        self.tree: _Tree = []

        # List of all open constructs.
        self._open_blocks: list[_Block] = []

        # Stack storing the number of the open-blocks whenever a new file is opened, contains one
        # entry per open file.
        self._open_blocks_per_file: list[int] = []

        # Tree to which new nodes are currently appended
        self._curtree: _Tree = self.tree

        # Stack of trees, so that the previous tree can be restored as _curtree once the
        # corresponding entry in _open_blocks is closed again (one entry per open block).
        self._parent_trees: list[_Tree] = []

        # Current file
        self._file: str | None = None

    def reset(self) -> None:
        """Resets the builder so that it starts to build a new tree."""
        self.tree = []
        self._open_blocks = []
        self._parent_trees = []
        self._open_blocks_per_file = []
        self._curtree = self.tree
        self._file = None

    def handle_include(self, span: Span | None, fname: str) -> None:
        """Should be called to signalize change to new file.

        Args:
            span: Start and end line of the include directive or None if called the first time for
                the main input.
            fname: Name of the file to be included.
        """
        self._parent_trees.append(self._curtree)
        self._open_blocks.append(_IncludeBlock(self._file, [span], fname))
        self._file = fname
        self._curtree = []
        self._open_blocks_per_file.append(len(self._open_blocks))

    def handle_endinclude(self, span: Span | None, fname: str) -> None:
        """Should be called when processing of a file finished.

        Args:
            span: Start and end line of the include directive or None if called the first time for
                the main input.
            fname: Name of the file which has been included.
        """
        block_depth = self._open_blocks_per_file.pop(-1)
        if len(self._open_blocks) > block_depth:
            unclosed = self._open_blocks[-1]
            msg = f"{unclosed.directive} directive still unclosed when reaching end of file"
            raise FyppFatalError(msg, self._file, unclosed.spans[0])
        block = self._open_blocks.pop(-1)
        if not isinstance(block, _IncludeBlock):
            msg = f"internal error: last open block is not 'include' when closing file '{fname}'"
            raise FyppFatalError(msg)
        if span != block.spans[0]:
            msg = (
                "internal error: span for include and endinclude differ "
                f"({span} vs {block.spans[0]})"
            )
            raise FyppFatalError(msg)
        if fname != block.includefname:
            msg = (
                "internal error: mismatching file name in close_file event "
                f"(expected: '{block.includefname}', got: '{fname}')"
            )
            raise FyppFatalError(msg)
        block.tree = self._curtree
        self._curtree = self._parent_trees.pop(-1)
        self._curtree.append(block)
        self._file = block.fname

    def handle_if(self, span: Span, cond: str) -> None:
        """Should be called to signalize an if directive.

        Args:
            span: Start and end line of the directive.
            cond: String representation of the branching condition.
        """
        self._parent_trees.append(self._curtree)
        self._open_blocks.append(_IfBlock(self._file, [span], [cond], []))
        self._curtree = []

    def handle_elif(self, span: Span, cond: str) -> None:
        """Should be called to signalize an elif directive.

        Args:
            span: Start and end line of the directive.
            cond: String representation of the branching condition.
        """
        self._check_current_file_has_open_block(span, "elif")
        block = self._open_blocks[-1]
        self._check_last_directive_matches(block.directive, "if", block.spans[-1], span, "elif")
        assert isinstance(block, _IfBlock)
        block.trees.append(self._curtree)
        block.conditions.append(cond)
        block.spans.append(span)
        self._curtree = []

    def handle_else(self, span: Span) -> None:
        """Should be called to signalize an else directive.

        Args:
            span: Start and end line of the directive.
        """
        self._check_current_file_has_open_block(span, "else")
        block = self._open_blocks[-1]
        self._check_last_directive_matches(block.directive, "if", block.spans[-1], span, "else")
        assert isinstance(block, _IfBlock)
        block.trees.append(self._curtree)
        block.conditions.append("True")
        block.spans.append(span)
        self._curtree = []

    def handle_endif(self, span: Span) -> None:
        """Should be called to signalize an endif directive.

        Args:
            span: Start and end line of the directive.
        """
        self._check_current_file_has_open_block(span, "endif")
        block = self._open_blocks.pop(-1)
        self._check_last_directive_matches(block.directive, "if", block.spans[-1], span, "endif")
        assert isinstance(block, _IfBlock)
        block.trees.append(self._curtree)
        block.spans.append(span)
        self._curtree = self._parent_trees.pop(-1)
        self._curtree.append(block)

    def handle_for(self, span: Span, loopvar: list[str], iterator: str) -> None:
        """Should be called to signalize a for directive.

        Args:
            span: Start and end line of the directive.
            loopvar: String representation of the loop variable expression.
            iterator: String representation of the iterable.
        """
        self._parent_trees.append(self._curtree)
        self._open_blocks.append(_ForBlock(self._file, [span], loopvar, iterator))
        self._curtree = []

    def handle_endfor(self, span: Span) -> None:
        """Should be called to signalize an endfor directive.

        Args:
            span: Start and end line of the directive.
        """
        self._check_current_file_has_open_block(span, "endfor")
        block = self._open_blocks.pop(-1)
        self._check_last_directive_matches(block.directive, "for", block.spans[-1], span, "endfor")
        assert isinstance(block, _ForBlock)
        block.tree = self._curtree
        block.spans.append(span)
        self._curtree = self._parent_trees.pop(-1)
        self._curtree.append(block)

    def handle_def(self, span: Span, name: str, argexpr: str | None) -> None:
        """Should be called to signalize a def directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the macro to be defined.
            argexpr: Macro argument definition or None
        """
        self._parent_trees.append(self._curtree)
        self._open_blocks.append(_DefBlock(self._file, [span], name, argexpr))
        self._curtree = []

    def handle_enddef(self, span: Span, name: str | None) -> None:
        """Should be called to signalize an enddef directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the enddef statement. Could be None, if enddef was specified without name.
        """
        self._check_current_file_has_open_block(span, "enddef")
        block = self._open_blocks.pop(-1)
        self._check_last_directive_matches(block.directive, "def", block.spans[-1], span, "enddef")
        assert isinstance(block, _DefBlock)
        if name is not None and name != block.name:
            msg = f"wrong name in enddef directive (expected '{block.name}', got '{name}')"
            raise FyppFatalError(msg, block.fname, span)
        block.tree = self._curtree
        block.spans.append(span)
        self._curtree = self._parent_trees.pop(-1)
        self._curtree.append(block)

    def handle_call(self, span: Span, name: str, argexpr: str | None, blockcall: bool) -> None:
        """Should be called to signalize a call directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the callable to call
            argexpr: Argument expression containing additional arguments for the call.
            blockcall: Whether the alternative "block / contains / endblock" calling directive has
                been used.
        """
        directive = "block" if blockcall else "call"
        self._parent_trees.append(self._curtree)
        self._open_blocks.append(_CallBlock(directive, self._file, [span, span], name, argexpr))
        self._curtree = []

    def handle_nextarg(self, span: Span, name: str | None, blockcall: bool) -> None:
        """Should be called to signalize a nextarg directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the argument following next or None if it should be the next positional
                argument.
            blockcall: Whether the alternative "block / contains / endblock" calling directive has
                been used.
        """
        self._check_current_file_has_open_block(span, "nextarg")
        block = self._open_blocks[-1]
        if blockcall:
            opened, current = "block", "contains"
        else:
            opened, current = "call", "nextarg"
        self._check_last_directive_matches(block.directive, opened, block.spans[-1], span, current)
        assert isinstance(block, _CallBlock)
        block.trees.append(self._curtree)
        block.spans.append(span)
        if name is not None:
            block.argnames.append(name)
        elif block.argnames:
            msg = "non-keyword argument following keyword argument"
            raise FyppFatalError(msg, block.fname, span)
        self._curtree = []

    def handle_endcall(self, span: Span, name: str | None, blockcall: bool) -> None:
        """Should be called to signalize an endcall directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the endcall statement. Could be None, if endcall was specified without
                name.
            blockcall: Whether the alternative "block / contains / endblock" calling directive has
                been used.
        """
        self._check_current_file_has_open_block(span, "endcall")
        block = self._open_blocks.pop(-1)
        if blockcall:
            opened, current = "block", "endblock"
        else:
            opened, current = "call", "endcall"
        self._check_last_directive_matches(block.directive, opened, block.spans[0], span, current)
        assert isinstance(block, _CallBlock)
        if name is not None and name != block.name:
            msg = f"wrong name in {current} directive (expected '{block.name}', got '{name}')"
            raise FyppFatalError(msg, block.fname, span)
        block.trees.append(self._curtree)
        # If nextarg or endcall immediately followed call, then first argument
        # is empty and should be removed (to allow for calls without arguments
        # and named first argument in calls)
        if block.trees and not block.trees[0]:
            if len(block.argnames) == len(block.trees):
                del block.argnames[0]
            del block.trees[0]
            del block.spans[1]
        block.spans.append(span)
        self._curtree = self._parent_trees.pop(-1)
        self._curtree.append(block)

    def handle_set(self, span: Span, name: str, expr: str) -> None:
        """Should be called to signalize a set directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the variable.
            expr: String representation of the expression to be assigned to the variable.
        """
        self._curtree.append(_SetDirective(self._file, span, name, expr))

    def handle_global(self, span: Span, name: str) -> None:
        """Should be called to signalize a global directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the variable(s) to make global.
        """
        self._curtree.append(_GlobalDirective(self._file, span, name))

    def handle_del(self, span: Span, name: str) -> None:
        """Should be called to signalize a del directive.

        Args:
            span: Start and end line of the directive.
            name: Name of the variable(s) to delete.
        """
        self._curtree.append(_DelDirective(self._file, span, name))

    def handle_eval(self, span: Span, expr: str) -> None:
        """Should be called to signalize an eval directive.

        Args:
            span: Start and end line of the directive.
            expr: String representation of the Python expression to be evaluated.
        """
        self._curtree.append(_EvalDirective(self._file, span, expr))

    def handle_comment(self, span: Span) -> None:
        """Should be called to signalize a comment directive.

        The content of the comment is not needed by the builder, but it needs
        the span of the comment to generate proper line numbers if needed.

        Args:
            span: Start and end line of the directive.
        """
        self._curtree.append(_CommentDirective(self._file, span))

    def handle_text(self, span: Span, txt: str) -> None:
        """Should be called to pass text which goes to output unaltered.

        Args:
            span: Start and end line of the text.
            txt: Text.
        """
        self._curtree.append(_RawText(self._file, span, txt))

    def handle_mute(self, span: Span) -> None:
        """Should be called to signalize a mute directive.

        Args:
            span: Start and end line of the directive.
        """
        self._parent_trees.append(self._curtree)
        self._open_blocks.append(_MuteBlock(self._file, [span]))
        self._curtree = []

    def handle_endmute(self, span: Span) -> None:
        """Should be called to signalize an endmute directive.

        Args:
            span: Start and end line of the directive.
        """
        self._check_current_file_has_open_block(span, "endmute")
        block = self._open_blocks.pop(-1)
        self._check_last_directive_matches(
            block.directive, "mute", block.spans[-1], span, "endmute"
        )
        assert isinstance(block, _MuteBlock)
        block.tree = self._curtree
        block.spans.append(span)
        self._curtree = self._parent_trees.pop(-1)
        self._curtree.append(block)

    def handle_stop(self, span: Span, msg: str) -> None:
        """Should be called to signalize a stop directive.

        Args:
            span: Start and end line of the directive.
            msg: Message to print before stopping.
        """
        self._curtree.append(_StopDirective(self._file, span, msg))

    def handle_assert(self, span: Span, cond: str) -> None:
        """Should be called to signalize an assert directive.

        Args:
            span: Start and end line of the directive.
            cond: Condition to assert
        """
        self._curtree.append(_AssertDirective(self._file, span, cond))

    def _check_current_file_has_open_block(self, span: Span, directive: str) -> None:
        if len(self._open_blocks) <= self._open_blocks_per_file[-1]:
            msg = f"unexpected {directive} directive"
            raise FyppFatalError(msg, self._file, span)

    def _check_last_directive_matches(
        self, lastdir: str, curdir: str, lastspan: Span | None, curspan: Span, directive: str
    ) -> None:
        if curdir != lastdir:
            msg = f"mismatching '{directive}' directive (last block opened was '{lastdir}')"
            raise FyppFatalError(msg, self._file, curspan)
        # Only the top-level file's _IncludeBlock ever has a None span entry, and this method is
        # only called for the other (if/for/def/call/mute) block types.
        assert lastspan is not None
        inline_last = lastspan.start == lastspan.end
        inline_cur = curspan.start == curspan.end
        if inline_last != inline_cur:
            if inline_cur:
                msg = f"expecting line form of directive {directive}"
            else:
                msg = f"expecting inline form of directive {directive}"
            raise FyppFatalError(msg, self._file, curspan)
        if inline_cur and curspan.start != lastspan.start:
            msg = "inline directives of the same construct must be in the same line"
            raise FyppFatalError(msg, self._file, curspan)


#
# Data structures used by the renderer
#


class _RenderResult(NamedTuple):
    """Result of rendering a (sub)tree.

    Fields:
        fragments: Rendered output, one fragment per list entry.
        eval_slots: Indices into ``fragments`` of entries which resulted from evaluating an eval
            directive.
        eval_sources: For each entry in ``eval_slots``, the span and file name it originated from.
    """

    fragments: list[str]
    eval_slots: Sequence[int] = ()
    eval_sources: Sequence[tuple[Span, str | None]] = ()


class _CallableArgSpec(NamedTuple):
    """Argument spec of a callable, as needed for defining a macro from it."""

    args: list[str]
    defaults: dict[str, Any]
    varpos: str | None
    varkw: str | None


#
# Renderer
#


class Renderer:
    """Renders a tree.

    Args:
        evaluator: Evaluator to use when rendering eval directives. If None (default), Evaluator()
            is used.
        linenums: Whether linenums should be generated, defaults to False.
        contlinenums: Whether linenums for continuation should be generated, defaults to False.
        linenumformat: 'std', 'cpp' or 'gfortran5' depending what kind of line directives should be
            created. Default: 'cpp'. Format 'std' emits #line pragmas, 'cpp' resembles GNU cpps
            special format, and 'gfortran5' adds to cpp a workaround for a bug introduced in
            GFortran 5.
        linefolder: Callable to use when folding a line.
        filevarroot: render _FILE_ and _THIS_FILE_ as paths relative to this root directory
            (default: paths are not converted explicitly to relative paths)
    """

    def __init__(
        self,
        evaluator: Evaluator | None = None,
        linenums: bool = False,
        contlinenums: bool = False,
        linenumformat: str | None = None,
        linefolder: Callable[[str], list[str]] | None = None,
        filevarroot: str | None = None,
    ):
        # Evaluator to use for Python expressions
        self._evaluator: Evaluator = Evaluator() if evaluator is None else evaluator
        self._evaluator.updateglobals(_SYSTEM_=platform.system(), _MACHINE_=platform.machine())

        # Whether rendered output is diverted and will be processed
        # further before output (if True: no line numbering and post processing)
        self._diverted: bool = False

        # Whether file name and line numbers should be kept fixed and
        # not updated (typically when rendering macro content)
        self._fixedposition: bool = False

        # Whether line numbering directives should be emitted
        self._linenums: bool = linenums

        # Whether line numbering directives in continuation lines are needed.
        self._contlinenums: bool = contlinenums

        # Line number formatter function and whether gfortran5 fix is needed
        self._linenumdir: Callable[..., str]
        if linenumformat is None or linenumformat in ("cpp", "gfortran5"):
            self._linenumdir = linenumdir_cpp
            self._linenum_gfortran5: bool = linenumformat == "gfortran5"
        else:
            self._linenumdir = linenumdir_std
            self._linenum_gfortran5 = False

        # Callable to be used for folding lines
        self._linefolder: Callable[[str], list[str]]
        if linefolder is None:
            self._linefolder = lambda line: [line]
        else:
            self._linefolder = linefolder

        self._convert_file_path: Callable[[str], str | pathlib.Path]
        if filevarroot is None:
            self._convert_file_path = lambda path: path
        else:
            self._convert_file_path = lambda path: pathlib.Path(path).relative_to(filevarroot)

    def render(self, tree: _Tree, divert: bool = False, fixposition: bool = False) -> str:
        """Renders a tree.

        Args:
            tree: Tree to render.
            divert: Whether output will be diverted and sent for further processing, so that no
                line numbering directives and postprocessing are needed at this stage.
                (Default: False)
            fixposition: Whether file name and line position (variables _FILE_ and _LINE_) should
                be kept at their current values or should be updated continuously. (Default: False).

        Returns:
            Rendered string.
        """
        diverted = self._diverted
        fixedposition_old = self._fixedposition
        self._diverted = divert
        self._fixedposition = self._fixedposition or fixposition
        try:
            output, eval_slots, eval_sources = self._render(tree)
            if not self._diverted and eval_slots:
                self._postprocess_eval_lines(output, eval_slots, eval_sources)
        finally:
            self._diverted = diverted
            self._fixedposition = fixedposition_old
        txt = "".join(output)
        return txt

    def _render(self, tree: _Tree) -> _RenderResult:
        fragments: list[str] = []
        eval_slots: list[int] = []
        eval_sources: list[tuple[Span, str | None]] = []
        for node in tree:
            if isinstance(node, _RawText):
                result = _RenderResult([node.text])
            elif isinstance(node, _EvalDirective):
                result = self._get_eval(node)
            elif isinstance(node, _IfBlock):
                result = self._get_conditional_content(node)
            elif isinstance(node, _DefBlock):
                result = self._define_macro(node)
            elif isinstance(node, _SetDirective):
                result = self._define_variable(node)
            elif isinstance(node, _DelDirective):
                result = self._delete_variable(node)
            elif isinstance(node, _ForBlock):
                result = self._get_iterated_content(node)
            elif isinstance(node, _CallBlock):
                result = self._get_called_content(node)
            elif isinstance(node, _IncludeBlock):
                result = self._get_included_content(node)
            elif isinstance(node, _CommentDirective):
                result = self._get_comment(node)
            elif isinstance(node, _MuteBlock):
                result = self._get_muted_content(node)
            elif isinstance(node, _StopDirective):
                self._handle_stop(node)
            elif isinstance(node, _AssertDirective):
                result = self._handle_assert(node)
            elif isinstance(node, _GlobalDirective):
                result = self._add_global(node)
            else:
                msg = f"internal error: unknown command '{type(node).__name__}'"
                raise FyppFatalError(msg)
            eval_slots += _shiftinds(result.eval_slots, len(fragments))
            eval_sources += result.eval_sources
            fragments += result.fragments
        return _RenderResult(fragments, eval_slots, eval_sources)

    def _get_eval(self, node: _EvalDirective) -> _RenderResult:
        fname = node.fname
        span = node.span
        try:
            result = self._evaluate(node.expr, fname, span.start)
        except Exception as exc:
            msg = f"exception occurred when evaluating '{node.expr}'"
            raise FyppFatalError(msg, fname, span) from exc
        fragments = []
        eval_slots = []
        eval_sources = []
        if result is not None:
            fragments.append(str(result))
            if not self._diverted:
                eval_slots.append(0)
                eval_sources.append((span, fname))
        if span.start != span.end:
            fragments.append("\n")
        return _RenderResult(fragments, eval_slots, eval_sources)

    def _get_conditional_content(self, node: _IfBlock) -> _RenderResult:
        fragments = []
        eval_slots = []
        eval_sources: list[tuple[Span, str | None]] = []
        fname = node.fname
        spans = node.spans
        multiline = spans[0].start != spans[-1].end
        for condition, tree, span in zip(node.conditions, node.trees, spans):
            try:
                cond = bool(self._evaluate(condition, fname, span.start))
            except Exception as exc:
                msg = f"exception occurred when evaluating '{condition}'"
                raise FyppFatalError(msg, fname, span) from exc
            if cond:
                if self._linenums and not self._diverted and multiline:
                    fragments.append(self._linenumdir(span.end, fname))
                result = self._render(tree)
                eval_slots += _shiftinds(result.eval_slots, len(fragments))
                eval_sources += result.eval_sources
                fragments += result.fragments
                break
        if self._linenums and not self._diverted and multiline:
            fragments.append(self._linenumdir(spans[-1].end, fname))
        return _RenderResult(fragments, eval_slots, eval_sources)

    def _get_iterated_content(self, node: _ForBlock) -> _RenderResult:
        fragments = []
        eval_slots = []
        eval_sources: list[tuple[Span, str | None]] = []
        fname = node.fname
        spans = node.spans
        try:
            iterobj = iter(self._evaluate(node.loopiter, fname, spans[0].start))
        except Exception as exc:
            msg = f"exception occurred when evaluating '{node.loopiter}'"
            raise FyppFatalError(msg, fname, spans[0]) from exc
        multiline = spans[0].start != spans[-1].end
        # Set by handle_endfor() once the for/endfor pair is closed; _get_iterated_content() is
        # only ever called on a closed block (i.e. an already-finished tree node).
        assert node.tree is not None
        for var in iterobj:
            if len(node.loopvars) == 1:
                self._define(node.loopvars[0], var)
            else:
                for varname, value in zip(node.loopvars, var):
                    self._define(varname, value)
            if self._linenums and not self._diverted and multiline:
                fragments.append(self._linenumdir(spans[0].end, fname))
            result = self._render(node.tree)
            eval_slots += _shiftinds(result.eval_slots, len(fragments))
            eval_sources += result.eval_sources
            fragments += result.fragments
        if self._linenums and not self._diverted and multiline:
            fragments.append(self._linenumdir(spans[-1].end, fname))
        return _RenderResult(fragments, eval_slots, eval_sources)

    def _get_called_content(self, node: _CallBlock) -> _RenderResult:
        fname = node.fname
        spans = node.spans
        posargs, kwargs = self._get_call_arguments(node)
        try:
            callobj = self._evaluate(node.name, fname, spans[0].start)
            result = callobj(*posargs, **kwargs)
        except Exception as exc:
            msg = f"exception occurred when calling '{node.name}'"
            raise FyppFatalError(msg, fname, spans[0]) from exc
        self._update_predef_globals(fname, spans[0].start)
        span = Span(spans[0].start, spans[-1].end)
        fragments = []
        eval_slots = []
        eval_sources = []
        if result is not None:
            fragments = [str(result)]
            if not self._diverted:
                eval_slots = [0]
                eval_sources = [(span, fname)]
        if span.start != span.end:
            fragments.append("\n")
        return _RenderResult(fragments, eval_slots, eval_sources)

    def _get_call_arguments(self, node: _CallBlock) -> tuple[list[Any], dict[str, Any]]:
        fname = node.fname
        spans = node.spans
        argexpr = node.argexpr
        if argexpr is None:
            posargs = []
            kwargs = {}
        else:
            # Parse and evaluate arguments passed in call header
            with self._evaluator.newscope():
                try:
                    posargs, kwargs = self._evaluate(
                        "__getargvalues(" + argexpr + ")", fname, spans[0].start
                    )
                except Exception as exc:
                    msg = f"unable to parse argument expression '{argexpr}'"
                    raise FyppFatalError(msg, fname, spans[0]) from exc

        # Render arguments passed in call body
        args = []
        for tree in node.trees:
            with self._evaluator.newscope():
                rendered = self.render(tree, divert=True)
            if rendered.endswith("\n"):
                rendered = rendered[:-1]
            args.append(rendered)

        # Separate arguments in call body into positional and keyword ones:
        argnames = node.argnames
        if argnames:
            posargs += args[: len(args) - len(argnames)]
            offset = len(args) - len(argnames)
            for iargname, argname in enumerate(argnames):
                ind = offset + iargname
                if argname in kwargs:
                    msg = f"keyword argument '{argname}' already defined"
                    raise FyppFatalError(msg, fname, spans[ind + 1])
                kwargs[argname] = args[ind]
        else:
            posargs += args

        return posargs, kwargs

    def _get_included_content(self, node: _IncludeBlock) -> _RenderResult:
        fname = node.fname
        spans = node.spans
        includefname = node.includefname
        firstspan = spans[0]
        includefile = firstspan is not None
        fragments: list[str] = []
        if self._linenums and not self._diverted:
            if includefile or self._linenum_gfortran5:
                fragments.append(self._linenumdir(0, includefname, _LINENUM_NEW_FILE))
            else:
                fragments.append(self._linenumdir(0, includefname))
        # Set by handle_endinclude() once the include is closed; _get_included_content() is only
        # ever called on a closed block (i.e. an already-finished tree node).
        assert node.tree is not None
        result = self._render(node.tree)
        eval_slots = _shiftinds(result.eval_slots, len(fragments))
        fragments += result.fragments
        if self._linenums and not self._diverted and includefile:
            assert firstspan is not None
            fragments.append(self._linenumdir(firstspan.end, fname, _LINENUM_RETURN_TO_FILE))
        return _RenderResult(fragments, eval_slots, result.eval_sources)

    def _define_macro(self, node: _DefBlock) -> _RenderResult:
        fname = node.fname
        spans = node.spans
        name = node.name
        argexpr = node.argexpr
        if argexpr is None:
            argspec = _CallableArgSpec([], {}, None, None)
        else:
            # Try to create a lambda function with the argument expression
            lambdaexpr = "lambda " + argexpr + ": None"
            with self._evaluator.newscope():
                try:
                    func = self._evaluate(lambdaexpr, fname, spans[0].start)
                except Exception as exc:
                    msg = f"exception occurred when evaluating argument expression '{argexpr}'"
                    raise FyppFatalError(msg, fname, spans[0]) from exc
            try:
                argspec = _get_callable_argspec(func)
            except Exception as exc:
                msg = f"invalid argument expression '{argexpr}'"
                raise FyppFatalError(msg, fname, spans[0]) from exc
            named_args = argspec.args if argspec.varpos is None else argspec.args + [argspec.varpos]
            named_args = named_args if argspec.varkw is None else named_args + [argspec.varkw]
            for arg in named_args:
                if arg in _RESERVED_NAMES or arg.startswith(_RESERVED_PREFIX):
                    msg = f"invalid argument name '{arg}'"
                    raise FyppFatalError(msg, fname, spans[0])
        try:
            # Set by handle_enddef() once the def/enddef pair is closed; _define_macro() is only
            # ever called on a closed block (i.e. an already-finished tree node).
            assert node.tree is not None
            macro = _Macro(
                name,
                fname,
                spans,
                argspec,
                node.tree,
                self,
                self._evaluator,
                self._evaluator.localscope,
            )
            self._define(name, macro)
        except Exception as exc:
            msg = f"exception occurred when defining macro '{name}'"
            raise FyppFatalError(msg, fname, spans[0]) from exc
        fragments = []
        if self._linenums and not self._diverted:
            fragments.append(self._linenumdir(spans[1].end, fname))
        return _RenderResult(fragments)

    def _define_variable(self, node: _SetDirective) -> _RenderResult:
        fname = node.fname
        span = node.span
        name = node.name
        valstr = node.expr
        try:
            if valstr is None:
                expr = None
            else:
                expr = self._evaluate(valstr, fname, span.start)
            self._define(name, expr)
        except Exception as exc:
            msg = f"exception occurred when setting variable(s) '{name}' to '{valstr}'"
            raise FyppFatalError(msg, fname, span) from exc
        multiline = span.start != span.end
        fragments = []
        if self._linenums and not self._diverted and multiline:
            fragments.append(self._linenumdir(span.end, fname))
        return _RenderResult(fragments)

    def _delete_variable(self, node: _DelDirective) -> _RenderResult:
        fname = node.fname
        span = node.span
        name = node.name
        try:
            self._evaluator.undefine(name)
        except Exception as exc:
            msg = f"exception occurred when deleting variable(s) '{name}'"
            raise FyppFatalError(msg, fname, span) from exc
        multiline = span.start != span.end
        fragments = []
        if self._linenums and not self._diverted and multiline:
            fragments.append(self._linenumdir(span.end, fname))
        return _RenderResult(fragments)

    def _add_global(self, node: _GlobalDirective) -> _RenderResult:
        fname = node.fname
        span = node.span
        name = node.name
        try:
            self._evaluator.addglobal(name)
        except Exception as exc:
            msg = f"exception occurred when making variable(s) '{name}' global"
            raise FyppFatalError(msg, fname, span) from exc
        multiline = span.start != span.end
        fragments = []
        if self._linenums and not self._diverted and multiline:
            fragments.append(self._linenumdir(span.end, fname))
        return _RenderResult(fragments)

    def _get_comment(self, node: _CommentDirective) -> _RenderResult:
        fragments = []
        if self._linenums and not self._diverted:
            fragments.append(self._linenumdir(node.span.end, node.fname))
        return _RenderResult(fragments)

    def _get_muted_content(self, node: _MuteBlock) -> _RenderResult:
        # Set by handle_endmute() once the mute/endmute pair is closed; _get_muted_content() is
        # only ever called on a closed block (i.e. an already-finished tree node).
        assert node.tree is not None
        self._render(node.tree)
        fragments = []
        if self._linenums and not self._diverted:
            fragments.append(self._linenumdir(node.spans[-1].end, node.fname))
        return _RenderResult(fragments)

    def _handle_stop(self, node: _StopDirective) -> NoReturn:
        fname = node.fname
        span = node.span
        msgstr = node.msg
        try:
            msg = str(self._evaluate(msgstr, fname, span.start))
        except Exception as exc:
            msg = f"exception occurred when evaluating stop message '{msgstr}'"
            raise FyppFatalError(msg, fname, span) from exc
        raise FyppStopRequest(msg, fname, span)

    def _handle_assert(self, node: _AssertDirective) -> _RenderResult:
        fname = node.fname
        span = node.span
        expr = node.cond
        try:
            cond = bool(self._evaluate(expr, fname, span.start))
        except Exception as exc:
            msg = f"exception occurred when evaluating assert condition '{expr}'"
            raise FyppFatalError(msg, fname, span) from exc
        if not cond:
            msg = f"Assertion failed ('{expr}')"
            raise FyppStopRequest(msg, fname, span)
        fragments = []
        if self._linenums and not self._diverted:
            fragments.append(self._linenumdir(span.end, fname))
        return _RenderResult(fragments)

    def _evaluate(self, expr: str, fname: str | None, linenr: int) -> Any:
        self._update_predef_globals(fname, linenr)
        result = self._evaluator.evaluate(expr)
        self._update_predef_globals(fname, linenr)
        return result

    def _update_predef_globals(self, fname, linenr):
        # fname is only None for the (nodeless) top-level include event; every node that can
        # reach an eval directive (and hence this method) always carries a real file name.
        assert fname is not None
        fname = _normalized_path(str(self._convert_file_path(fname)))
        self._evaluator.updatelocals(
            _DATE_=time.strftime("%Y-%m-%d"),
            _TIME_=time.strftime("%H:%M:%S"),
            _THIS_FILE_=fname,
            _THIS_LINE_=linenr + 1,
        )
        if not self._fixedposition:
            self._evaluator.updateglobals(_FILE_=fname, _LINE_=linenr + 1)

    def _define(self, var: str, value: Any) -> None:
        self._evaluator.define(var, value)

    def _postprocess_eval_lines(
        self,
        output: list[str],
        eval_slots: Sequence[int],
        eval_sources: Sequence[tuple[Span, str | None]],
    ) -> None:
        ilastproc = -1
        for ieval, ind in enumerate(eval_slots):
            span, fname = eval_sources[ieval]
            if ind <= ilastproc:
                continue
            iprev, eolprev = self._find_last_eol(output, ind)
            inext, eolnext = self._find_next_eol(output, ind)
            curline = self._glue_line(output, ind, iprev, eolprev, inext, eolnext)
            output[iprev + 1 : inext] = [""] * (inext - iprev - 1)
            output[ind] = self._postprocess_eval_line(curline, fname, span)
            ilastproc = inext

    @staticmethod
    def _find_last_eol(output: list[str], ind: int) -> tuple[int, int]:
        "Find last newline before current position."
        iprev = ind - 1
        while iprev >= 0:
            eolprev = output[iprev].rfind("\n")
            if eolprev != -1:
                break
            iprev -= 1
        else:
            iprev = 0
            eolprev = -1
        return iprev, eolprev

    @staticmethod
    def _find_next_eol(output: list[str], ind: int) -> tuple[int, int]:
        "Find next newline after current position."
        # find first eol after expr. evaluation
        inext = ind + 1
        while inext < len(output):
            eolnext = output[inext].find("\n")
            if eolnext != -1:
                break
            inext += 1
        else:
            inext = len(output) - 1
            eolnext = len(output[-1]) - 1
        return inext, eolnext

    @staticmethod
    def _glue_line(
        output: list[str], ind: int, iprev: int, eolprev: int, inext: int, eolnext: int
    ) -> str:
        "Create line from parts between specified boundaries."
        curline_parts = []
        if iprev != ind:
            curline_parts = [output[iprev][eolprev + 1 :]]
            output[iprev] = output[iprev][: eolprev + 1]
        curline_parts.extend(output[iprev + 1 : ind])
        curline_parts.extend(output[ind])
        curline_parts.extend(output[ind + 1 : inext])
        if inext != ind:
            curline_parts.append(output[inext][: eolnext + 1])
            output[inext] = output[inext][eolnext + 1 :]
        return "".join(curline_parts)

    def _postprocess_eval_line(self, evalline: str, fname: str | None, span: Span) -> str:
        lines = evalline.split("\n")
        # If line ended on '\n', last element is ''. We remove it and
        # add the trailing newline later manually.
        trailing_newline = lines[-1] == ""
        if trailing_newline:
            del lines[-1]
        lnum = self._linenumdir(span.start, fname) if self._linenums else ""
        clnum = lnum if self._contlinenums else ""
        linenumsep = "\n" + lnum
        clinenumsep = "\n" + clnum
        foldedlines = [self._foldline(line) for line in lines]
        outlines = [clinenumsep.join(lines) for lines in foldedlines]
        result = linenumsep.join(outlines)
        # Add missing trailing newline
        if trailing_newline:
            trailing = "\n"
            if self._linenums:
                # For inline eval directives span.start == span.end
                # -> next line is span.start + 1 and not span.end as for line eval directives
                nextline = max(span.end, span.start + 1)
                trailing += self._linenumdir(nextline, fname)
        else:
            trailing = ""
        return result + trailing

    def _foldline(self, line: str) -> list[str]:
        if _COMMENTLINE_REGEXP.match(line) is None:
            return self._linefolder(line)
        return [line]


#
# Evaluator
#


class Evaluator:
    """Provides an isolated environment for evaluating Python expressions.

    It restricts the builtins which can be used within this environment to a
    (hopefully safe) subset. Additionally it defines the functions which are
    provided by the preprocessor for the eval directives.

    Args:
        env: Initial definitions for the environment, defaults to None.
    """

    # Restricted builtins working in all supported Python versions. Version
    # specific ones might be added dynamically in _get_restricted_builtins().
    _RESTRICTED_BUILTINS = {
        "abs": builtins.abs,
        "all": builtins.all,
        "any": builtins.any,
        "bin": builtins.bin,
        "bool": builtins.bool,
        "bytearray": builtins.bytearray,
        "bytes": builtins.bytes,
        "chr": builtins.chr,
        "classmethod": builtins.classmethod,
        "complex": builtins.complex,
        "delattr": builtins.delattr,
        "dict": builtins.dict,
        "dir": builtins.dir,
        "divmod": builtins.divmod,
        "enumerate": builtins.enumerate,
        "filter": builtins.filter,
        "float": builtins.float,
        "format": builtins.format,
        "frozenset": builtins.frozenset,
        "getattr": builtins.getattr,
        "globals": builtins.globals,
        "hasattr": builtins.hasattr,
        "hash": builtins.hash,
        "hex": builtins.hex,
        "id": builtins.id,
        "int": builtins.int,
        "isinstance": builtins.isinstance,
        "issubclass": builtins.issubclass,
        "iter": builtins.iter,
        "len": builtins.len,
        "list": builtins.list,
        "locals": builtins.locals,
        "map": builtins.map,
        "max": builtins.max,
        "min": builtins.min,
        "next": builtins.next,
        "object": builtins.object,
        "oct": builtins.oct,
        "ord": builtins.ord,
        "pow": builtins.pow,
        "property": builtins.property,
        "range": builtins.range,
        "repr": builtins.repr,
        "reversed": builtins.reversed,
        "round": builtins.round,
        "set": builtins.set,
        "setattr": builtins.setattr,
        "slice": builtins.slice,
        "sorted": builtins.sorted,
        "staticmethod": builtins.staticmethod,
        "str": builtins.str,
        "sum": builtins.sum,
        "super": builtins.super,
        "tuple": builtins.tuple,
        "type": builtins.type,
        "vars": builtins.vars,
        "zip": builtins.zip,
    }

    def __init__(self, env: dict | None = None):

        # Global scope
        self._globals: dict = env if env is not None else {}

        # Local scope(s)
        self._locals: dict | None = None
        self._locals_stack: list[dict | None] = []

        # Variables which are references to entries in global scope
        self._globalrefs: set[str] | None = None
        self._globalrefs_stack: list[set[str] | None] = []

        # Current scope (globals + locals in all embedding and in current scope)
        self._scope: dict = self._globals

        # Turn on restricted mode
        self._restrict_builtins()

    def evaluate(self, expr: str) -> Any:
        """Evaluate a Python expression using the `eval()` builtin.

        Args:
            expr: String represantion of the expression.

        Return:
            Python object: Result of the expression evaluation.
        """
        result = eval(expr, self._scope)
        return result

    def import_module(self, module: str) -> None:
        """Import a module into the evaluator.

        Note: Import only trustworthy modules! Module imports are global,
        therefore, importing a malicious module which manipulates other global
        modules could affect code behaviour outside of the Evaluator as well.

        Args:
            module: Python module to import.

        Raises:
            FyppFatalError: If module could not be imported.

        """
        rootmod = module.split(".", 1)[0]
        try:
            imported = __import__(module, self._scope)
            self.define(rootmod, imported)
        except Exception as exc:
            msg = f"failed to import module '{module}'"
            raise FyppFatalError(msg) from exc

    def define(self, name: str, value: Any) -> None:
        """Define a Python entity.

        Args:
            name: Name of the entity.
            value: Value of the entity.

        Raises:
            FyppFatalError: If name starts with the reserved prefix or if it is
                a reserved name.
        """
        varnames = self._get_variable_names(name)
        if len(varnames) == 1:
            value = (value,)
        elif len(varnames) != len(value):
            msg = "value for tuple assignment has incompatible length"
            raise FyppFatalError(msg)
        for varname, varvalue in zip(varnames, value):
            self._check_variable_name(varname)
            if self._locals is None:
                self._globals[varname] = varvalue
            else:
                assert self._globalrefs is not None
                if varname in self._globalrefs:
                    self._globals[varname] = varvalue
                else:
                    self._locals[varname] = varvalue
                self._scope[varname] = varvalue

    def undefine(self, name: str) -> None:
        """Undefine a Python entity.

        Args:
            name: Name of the entity to undefine.

        Raises:
            FyppFatalError: If name starts with the reserved prefix or if it is
                a reserved name.
        """
        varnames = self._get_variable_names(name)
        for varname in varnames:
            self._check_variable_name(varname)
            deleted = False
            if self._locals is None:
                if varname in self._globals:
                    del self._globals[varname]
                    deleted = True
            else:
                assert self._globalrefs is not None
                if varname in self._locals:
                    del self._locals[varname]
                    del self._scope[varname]
                    deleted = True
                elif varname in self._globalrefs and varname in self._globals:
                    del self._globals[varname]
                    del self._scope[varname]
                    deleted = True
            if not deleted:
                msg = f"lookup for an erasable instance of '{varname}' failed"
                raise FyppFatalError(msg)

    def addglobal(self, name: str) -> None:
        """Define a given entity as global.

        Args:
            name: Name of the entity to make global.

        Raises:
            FyppFatalError: If entity name is invalid or if the current scope is
                 a local scope and entity is already defined in it.
        """
        varnames = self._get_variable_names(name)
        for varname in varnames:
            self._check_variable_name(varname)
            if self._locals is not None:
                if varname in self._locals:
                    msg = f"variable '{varname}' already defined in local scope"
                    raise FyppFatalError(msg)
                assert self._globalrefs is not None
                self._globalrefs.add(varname)

    def updateglobals(self, **vardict: Any) -> None:
        """Update variables in the global scope.

        This is a shortcut function to inject protected variables in the global
        scope without extensive checks (as in define()). Vardict must not
        contain any global entries which can be shadowed in local scopes
        (e.g. should only contain variables with forbidden prefix).

        Args:
            **vardict: variable definitions.

        """
        self._scope.update(vardict)
        if self._locals is not None:
            self._globals.update(vardict)

    def updatelocals(self, **vardict: Any) -> None:
        """Update variables in the local scope.

        This is a shortcut function to inject variables in the local scope
        without extensive checks (as in define()). Vardict must not contain any
        entries which have been made global via addglobal() before. In order to
        ensure this, updatelocals() should be called immediately after
        openscope(), or with variable names, which are warrantedly not globals
        (e.g variables starting with forbidden prefix)

        Args:
            **vardict: variable definitions.
        """
        self._scope.update(vardict)
        if self._locals is not None:
            self._locals.update(vardict)

    def openscope(self, customlocals: dict | None = None) -> None:
        """Opens a new (embedded) scope.

        Note: consider to use newscope() to handle exceptions between openscope() and closescope().

        Args:
            customlocals: By default, the locals of the embedding scope are visible in the new one.
                When this is not the desired behaviour a dictionary of customized locals
                can be passed, and those locals will become the only visible ones.
        """
        self._locals_stack.append(self._locals)
        self._globalrefs_stack.append(self._globalrefs)
        if customlocals is not None:
            newlocals = customlocals.copy()
        elif self._locals is not None:
            newlocals = self._locals.copy()
        else:
            newlocals = {}
        self._locals = newlocals
        self._globalrefs = set()
        self._scope = self._globals.copy()
        self._scope.update(newlocals)

    def closescope(self) -> None:
        """Close scope and restore embedding scope.

        Note: consider to use newscope() to handle exceptions between openscope() and closescope().
        """
        self._locals = self._locals_stack.pop(-1)
        self._globalrefs = self._globalrefs_stack.pop(-1)
        if self._locals is not None:
            self._scope = self._globals.copy()
            self._scope.update(self._locals)
        else:
            self._scope = self._globals

    @contextlib.contextmanager
    def newscope(self, customlocals: dict | None = None) -> Generator[None]:
        """Context manager opening a new scope and closing it on exit (even if exception occurred).

        Args:
            customlocals: Custom local scope, see openscope().
        """
        self.openscope(customlocals)
        try:
            yield
        finally:
            self.closescope()

    @property
    def globalscope(self) -> dict:
        "Dictionary of the global scope."
        return self._globals

    @property
    def localscope(self) -> dict | None:
        "Dictionary of the current local scope."
        return self._locals

    def _restrict_builtins(self) -> None:
        builtindict = self._get_restricted_builtins()
        builtindict["__import__"] = self._func_import
        builtindict["defined"] = self._func_defined
        builtindict["setvar"] = self._func_setvar
        builtindict["getvar"] = self._func_getvar
        builtindict["delvar"] = self._func_delvar
        builtindict["globalvar"] = self._func_globalvar
        builtindict["__getargvalues"] = self._func_getargvalues
        self._globals["__builtins__"] = builtindict

    @classmethod
    def _get_restricted_builtins(cls) -> dict[str, Any]:
        bidict = dict(cls._RESTRICTED_BUILTINS)
        return bidict

    @staticmethod
    def _get_variable_names(varexpr: str) -> list[str]:
        lpar = varexpr.startswith("(")
        rpar = varexpr.endswith(")")
        if lpar != rpar:
            msg = f"unbalanced parenthesis around variable varexpr(s) in '{varexpr}'"
            raise FyppFatalError(msg, None, None)
        if lpar:
            varexpr = varexpr[1:-1]
        varnames = [s.strip() for s in varexpr.split(",")]
        return varnames

    @staticmethod
    def _check_variable_name(varname: str) -> None:
        if varname.startswith(_RESERVED_PREFIX):
            msg = f"Name '{varname}' starts with reserved prefix '{_RESERVED_PREFIX}'"
            raise FyppFatalError(msg, None, None)
        if varname in _RESERVED_NAMES:
            msg = f"Name '{varname}' is reserved and can not be redefined"
            raise FyppFatalError(msg, None, None)

    def _func_defined(self, var: str) -> bool:
        defined = var in self._scope
        return defined

    def _func_import(self, name: str, *_: Any, **__: Any) -> types.ModuleType:
        module = self._scope.get(name, None)
        if module is not None and isinstance(module, types.ModuleType):
            return module
        msg = f"Import of module '{name}' via '__import__' not allowed"
        raise ImportError(msg)

    def _func_setvar(self, *namesvalues: Any) -> None:
        if len(namesvalues) % 2:
            msg = "setvar function needs an even number of arguments"
            raise FyppFatalError(msg)
        for ind in range(0, len(namesvalues), 2):
            self.define(namesvalues[ind], namesvalues[ind + 1])

    def _func_getvar(self, name: str, defvalue: Any = None) -> Any:
        if name in self._scope:
            return self._scope[name]
        return defvalue

    def _func_delvar(self, *names: str) -> None:
        for name in names:
            self.undefine(name)

    def _func_globalvar(self, *names: str) -> None:
        for name in names:
            self.addglobal(name)

    @staticmethod
    def _func_getargvalues(*args: Any, **kwargs: Any) -> tuple[list, dict]:
        return list(args), kwargs


#
# Macro
#


@dataclasses.dataclass(frozen=True, slots=True)
class _Macro:
    """Represents a user defined macro.

    This object should only be initialized by a Renderer instance, as it needs access to Renderers
    internal variables and methods.

    Args:
        name: Name of the macro.
        fname: The file where the macro was defined.
        spans: Line spans of macro definition.
        argspec: Argument specification of the macro.
        tree: Body of the macro.
        renderer: Renderer to use for evaluating macro content.
        localscope: Dictionary with local variables, which should be used the local scope, when the
            macro is called. Default: None (empty local scope).
    """

    _name: str
    _fname: str | None
    _spans: list[Span]
    _argspec: _CallableArgSpec
    _tree: _Tree
    _renderer: Renderer
    _evaluator: Evaluator
    _localscope: dict | None = None

    def __post_init__(self) -> None:
        if self._localscope is None:
            object.__setattr__(self, "_localscope", {})

    def __call__(self, *args: Any, **keywords: Any) -> str:
        argdict = self._process_arguments(args, keywords)
        with self._evaluator.newscope(customlocals=self._localscope):
            self._evaluator.updatelocals(**argdict)
            output = self._renderer.render(self._tree, divert=True, fixposition=True)
        if output.endswith("\n"):
            return output[:-1]
        return output

    def _process_arguments(self, args: tuple[Any, ...], keywords: dict[str, Any]) -> dict[str, Any]:
        kwdict = dict(keywords)
        argdict: dict[str, Any] = {}
        nargs = min(len(args), len(self._argspec.args))
        for iarg in range(nargs):
            argdict[self._argspec.args[iarg]] = args[iarg]
        if nargs < len(args):
            if self._argspec.varpos is None:
                msg = (
                    f"macro '{self._name}' called with too many positional arguments "
                    f"(expected: {len(self._argspec.args)}, received: {len(args)})"
                )
                raise FyppFatalError(msg, self._fname, self._spans[0])
            argdict[self._argspec.varpos] = list(args[nargs:])
        elif self._argspec.varpos is not None:
            argdict[self._argspec.varpos] = []
        for argname in self._argspec.args[:nargs]:
            if argname in kwdict:
                msg = f"got multiple values for argument '{argname}'"
                raise FyppFatalError(msg, self._fname, self._spans[0])
        if nargs < len(self._argspec.args):
            for argname in self._argspec.args[nargs:]:
                if argname in kwdict:
                    argdict[argname] = kwdict.pop(argname)
                elif argname in self._argspec.defaults:
                    argdict[argname] = self._argspec.defaults[argname]
                else:
                    msg = (
                        f"macro '{self._name}' called without mandatory positional argument "
                        f"'{argname}'"
                    )
                    raise FyppFatalError(msg, self._fname, self._spans[0])
        if kwdict and self._argspec.varkw is None:
            kwstr = "', '".join(kwdict.keys())
            msg = f"macro '{self._name}' called with unknown keyword argument(s) '{kwstr}'"
            raise FyppFatalError(msg, self._fname, self._spans[0])
        if self._argspec.varkw is not None:
            argdict[self._argspec.varkw] = kwdict
        return argdict


#
# Processor
#


class Processor:
    """Connects various objects with each other to create a processor.

    Args:
        parser: Parser to use for parsing text. If None (default), `Parser()` is used.
        builder: Builder to use for building the tree representation of the text. If None (default),
            `Builder()` is used.
        renderer: Renderer to use for rendering the output. If None (default), `Renderer()` is used
            with a default Evaluator().
        evaluator: Evaluator to use for evaluating Python expressions. If None (default),
            `Evaluator()` is used.
    """

    def __init__(
        self,
        parser: Parser | None = None,
        builder: Builder | None = None,
        renderer: Renderer | None = None,
        evaluator: Evaluator | None = None,
    ):
        self._builder: Builder = Builder() if builder is None else builder
        if parser is None:
            self._parser: Parser = Parser(handlers=self._builder)
        else:
            self._parser = parser
            self._parser.handlers = self._builder

        self._renderer: Renderer
        if renderer is None:
            evaluator = Evaluator() if evaluator is None else evaluator
            self._renderer = Renderer(evaluator)
        else:
            self._renderer = renderer

    def process_file(self, fname: str) -> str:
        """Processes a file.

        Args:
            fname: Name of the file to process.

        Returns:
            Processed content.
        """
        self._builder.reset()
        self._parser.parsefile(fname)
        return self._render()

    def process_text(self, txt: str) -> str:
        """Processes a string.

        Args:
            txt: Text to process.

        Returns:
            Processed content.
        """
        self._builder.reset()
        self._parser.parse(txt)
        return self._render()

    def get_included_files(self) -> list[str]:
        """Returns the list of included files."""
        return self._parser.get_included_files()

    def _render(self) -> str:
        output = self._renderer.render(self._builder.tree)
        return output


#
# Fypp
#


class _FyppOptionsLike(Protocol):
    """Structural shape of the options object accepted by Fypp.__init__.

    Satisfied by a FyppOptions instance, an optparse.Values as returned by
    get_option_parser(), or any other object exposing the same attributes.
    """

    encoding: str
    modules: list[str]
    moduledirs: list[str]
    define_mode: str
    defines: list[str]
    defines_str: list[str]
    defines_eval: list[str]
    includes: list[str]
    fixed_format: bool
    no_folding: bool
    folding_mode: str
    line_length: int
    indentation: int
    line_numbering: bool
    line_numbering_mode: str
    create_parent_folder: bool
    line_marker_format: str
    file_var_root: str | None
    depfile: str | None


class Fypp:
    """Fypp preprocessor.

    You can invoke it like ::

        tool = fypp.Fypp()
        tool.process_file('file.in', 'file.out')

    to initialize Fypp with default options, process `file.in` and write the result to `file.out`.
    If the input should be read from a string, the ``process_text()`` method can be used::

        tool = fypp.Fypp()
        output = tool.process_text('#:if DEBUG > 0\\nprint *, "DEBUG"\\n#:endif\\n')

    If you want to fine tune Fypps behaviour, pass a customized `FyppOptions`_ instance at
    initialization::

        options = fypp.FyppOptions()
        options.fixed_format = True
        tool = fypp.Fypp(options)

    Alternatively, you can use the command line parser ``optparse.OptionParser`` to set options for
    Fypp. The function ``get_option_parser()`` returns you a default option parser. You can then use
    its ``parse_args()`` method to obtain settings by reading the command line arguments::

        optparser = fypp.get_option_parser()
        options, leftover = optparser.parse_args()
        tool = fypp.Fypp(options)

    The command line options can also be passed directly as a list when calling ``parse_args()``::

        args = ['-DDEBUG=0', 'input.fpp', 'output.f90']
        optparser = fypp.get_option_parser()
        options, leftover = optparser.parse_args(args=args)
        tool = fypp.Fypp(options)

    For even more fine-grained control over how Fypp works, you can pass in custom factory methods
    that handle construction of the evaluator, parser, builder and renderer components. These
    factory methods must have the same signature as the corresponding component's constructor. As an
    example of using a builder that's customized by subclassing::

        class MyBuilder(fypp.Builder):

            def __init__(self):
               super().__init__()
               ...additional initialization...

        tool = fypp.Fypp(options, builder_factory=MyBuilder)


    Args:
        options: Object containing the settings for Fypp. You typically would pass a customized
            `FyppOptions`_ instance or an ``optparse.Values`` object as returned by the option
            parser. If not present, the default settings in `FyppOptions`_ are used.
        evaluator_factory: Factory function that returns an Evaluator object. Its call signature
            must match that of the Evaluator constructor. If not present, ``Evaluator`` is used.
        parser_factory: Factory function that returns a Parser object.  Its call signature must
            match that of the Parser constructor. If not present, ``Parser`` is used.
        builder_factory: Factory function that returns a Builder object.  Its call signature must
            match that of the Builder constructor. If not present, ``Builder`` is used.
        renderer_factory: Factory function that returns a Renderer object. Its call signature must
            match that of the Renderer constructor.  If not present, ``Renderer`` is used.
    """

    def __init__(
        self,
        options: _FyppOptionsLike | None = None,
        evaluator_factory: Callable[..., Evaluator] = Evaluator,
        parser_factory: Callable[..., Parser] = Parser,
        builder_factory: Callable[..., Builder] = Builder,
        renderer_factory: Callable[..., Renderer] = Renderer,
    ):
        syspath = self._get_syspath_without_scriptdir()
        self._adjust_syspath(syspath)
        if options is None:
            options = FyppOptions()
        if inspect.signature(evaluator_factory) == inspect.signature(Evaluator):
            evaluator = evaluator_factory()
        else:
            raise FyppFatalError("evaluator_factory has incorrect signature")
        self._encoding: str = options.encoding
        if options.modules:
            self._import_modules(options.modules, evaluator, syspath, options.moduledirs)
        evaluate = options.define_mode == "eval"
        if options.defines:
            self._apply_definitions(options.defines, evaluator, evaluate)
        if options.defines_str:
            self._apply_definitions(options.defines_str, evaluator, False)
        if options.defines_eval:
            self._apply_definitions(options.defines_eval, evaluator, True)
        if inspect.signature(parser_factory) == inspect.signature(Parser):
            parser = parser_factory(includedirs=options.includes, encoding=self._encoding)
        else:
            raise FyppFatalError("parser_factory has incorrect signature")
        if inspect.signature(builder_factory) == inspect.signature(Builder):
            builder = builder_factory()
        else:
            raise FyppFatalError("builder_factory has incorrect signature")

        fixed_format = options.fixed_format
        linefolding = not options.no_folding
        linefolder: Callable[[str], list[str]]
        if linefolding:
            folding = "brute" if fixed_format else options.folding_mode
            linelength = 72 if fixed_format else options.line_length
            indentation = 5 if fixed_format else options.indentation
            prefix = "&"
            suffix = "" if fixed_format else "&"
            linefolder = FortranLineFolder(linelength, indentation, folding, prefix, suffix)
        else:
            linefolder = PlaceholderLineFolder()
        linenums = options.line_numbering
        contlinenums = options.line_numbering_mode != "nocontlines"
        self._create_parent_folder: bool = options.create_parent_folder
        if inspect.signature(renderer_factory) == inspect.signature(Renderer):
            renderer = renderer_factory(
                evaluator,
                linenums=linenums,
                contlinenums=contlinenums,
                linenumformat=options.line_marker_format,
                linefolder=linefolder,
                filevarroot=options.file_var_root,
            )
        else:
            raise FyppFatalError("renderer_factory has incorrect signature")
        self._preprocessor = Processor(parser, builder, renderer)
        self._depfile = options.depfile

    def process_file(self, infile: str, outfile: str | None = None) -> str | None:
        """Processes input file and writes result to output file.

        Args:
            infile: Name of the file to read and process. If its value is '-', input is read from
                stdin.
            outfile: Name of the file to write the result to. If its value is '-', result is written
                to stdout. If not present, result will be returned as string.

        Returns:
            Result of processed input, if no outfile was specified.
        """
        if self._depfile and (outfile == "-" or outfile is None):
            raise FyppFatalError(
                "dependency file writing requires an explicit output file for the processed input"
            )
        infile = STDIN_INPUT_NAME if infile == "-" else infile
        output = self._preprocessor.process_file(infile)
        if outfile is None:
            return output
        if outfile == "-":
            sys.stdout.write(output)
        else:
            with _open_output_file(outfile, self._encoding, self._create_parent_folder) as outfp:
                outfp.write(output)
            if self._depfile:
                included_files_uniq = list(dict.fromkeys(self.get_included_files()))
                _write_dependency_file(
                    self._depfile,
                    outfile,
                    included_files_uniq,
                    self._encoding,
                    self._create_parent_folder,
                )
        return None

    def process_text(self, txt: str) -> str:
        """Processes a string.

        Args:
            txt: String to process.

        Returns:
            Processed content.
        """
        return self._preprocessor.process_text(txt)

    def get_included_files(self) -> list[str]:
        """Returns the list of included files

        Returns:
            List of included file paths.
        """
        return self._preprocessor.get_included_files()

    @staticmethod
    def _apply_definitions(defines: list[str], evaluator: Evaluator, evaluate: bool) -> None:
        for define in defines:
            words = define.split("=", 1)
            name = words[0]
            if evaluate:
                value = None
                if len(words) > 1:
                    try:
                        value = evaluator.evaluate(words[1])
                    except Exception as exc:
                        msg = f"exception at evaluating '{words[1]}' in definition for '{name}'"
                        raise FyppFatalError(msg) from exc
            else:
                value = str(words[1]) if len(words) > 1 else ""
            evaluator.define(name, value)

    def _import_modules(
        self,
        modules: list[str],
        evaluator: Evaluator,
        syspath: list[str],
        moduledirs: list[str] | None,
    ) -> None:
        lookuppath = []
        if moduledirs is not None:
            lookuppath += [os.path.abspath(moddir) for moddir in moduledirs]
        lookuppath.append(os.path.abspath("."))
        lookuppath += syspath
        self._adjust_syspath(lookuppath)
        try:
            for module in modules:
                evaluator.import_module(module)
        finally:
            self._adjust_syspath(syspath)

    @staticmethod
    def _get_syspath_without_scriptdir() -> list[str]:
        """Remove the folder of the fypp binary from the search path"""
        syspath = list(sys.path)
        scriptdir = os.path.abspath(os.path.dirname(sys.argv[0]))
        if os.path.abspath(syspath[0]) == scriptdir:
            del syspath[0]
        return syspath

    @staticmethod
    def _adjust_syspath(syspath: list[str]) -> None:
        sys.path = syspath


class FyppOptions(optparse.Values):
    """Container for Fypp options with default values.

    Attributes:
        defines: List of variable definitions in the form of 'VARNAME=VALUE'. Default: []
        includes: List of paths to search when looking for include files. Default: []
        line_numbering: Whether line numbering directives should appear in the output.
            Default: False
        line_numbering_mode: Line numbering mode 'full' or 'nocontlines'. Default: 'full'.
        line_marker_format: Line marker format. Currently 'std', 'cpp' and 'gfortran5' are
            supported, where 'std' emits ``#line`` pragmas similar to standard tools, 'cpp'
            produces line directives as emitted by GNU cpp, and 'gfortran5' cpp line directives
            with a workaround for a bug introduced in GFortran 5. Default: 'cpp'.
        line_length: Length of output lines. Default: 132.
        folding_mode: Folding mode 'smart', 'simple' or 'brute'. Default: 'smart'.
        no_folding: Whether folding should be suppressed. Default: False.
        indentation: Indentation in continuation lines. Default: 4.
        modules: Modules to import at initialization. Default: [].
        moduledirs: Module lookup directories for importing user specified modules. The specified
            paths are looked up *before* the standard module locations in sys.path.
        fixed_format: Whether input file is in fixed format. Default: False.
        encoding: Character encoding for reading/writing files. Allowed values are Pythons codec
            identifiers, e.g. 'ascii', 'utf-8', etc. Default: 'utf-8'. Reading from stdin and
            writing to stdout is always encoded according to the current locale and is not
            affected by this setting.
        create_parent_folder: Whether the parent folder for the output file should be created if it
            does not exist. Default: False.
        file_var_root: Use relative file names in the file related Fypp variables relative to
            this folder. (If None, absolute file names are used).
        depfile: Name of the file to write dependency information into (if None,
            no dependency information is written.)
    """

    def __init__(self) -> None:
        optparse.Values.__init__(self)
        self.defines: list[str] = []
        self.define_mode: str = "eval"
        self.defines_str: list[str] = []
        self.defines_eval: list[str] = []
        self.includes: list[str] = []
        self.line_numbering: bool = False
        self.line_numbering_mode: str = "full"
        self.line_marker_format: str = "cpp"
        self.line_length: int = 132
        self.folding_mode: str = "smart"
        self.no_folding: bool = False
        self.indentation: int = 4
        self.modules: list[str] = []
        self.moduledirs: list[str] = []
        self.fixed_format: bool = False
        self.encoding: str = "utf-8"
        self.create_parent_folder: bool = False
        self.file_var_root: str | None = None
        self.depfile: str | None = None


class FortranLineFolder:
    """Implements line folding with Fortran continuation lines.

    Args:
        maxlen: Maximal line length (default: 132).
        indent: Indentation for continuation lines (default: 4).
        method: Folding method with following options:

            * ``brute``: folding with maximal length of continuation lines,
            * ``simple``: indents with respect of indentation of first line,
            * ``smart``: like ``simple``, but tries to fold at whitespaces.

        prefix: String to use at the beginning of a continuation line (default: '&').
        suffix: String to use at the end of the line preceding a continuation line (default: '&')
    """

    def __init__(
        self,
        maxlen: int = 132,
        indent: int = 4,
        method: str = "smart",
        prefix: str = "&",
        suffix: str = "&",
    ):
        # Line length should be long enough that contintuation lines can host at
        # east one character apart of indentation and two continuation signs
        minmaxlen = indent + len(prefix) + len(suffix) + 1
        if maxlen < minmaxlen:
            msg = f"Maximal line length less than {minmaxlen} when using an indentation of {indent}"
            raise FyppFatalError(msg)
        self._maxlen: int = maxlen
        self._indent: int = indent
        self._prefix: str = " " * self._indent + prefix
        self._suffix: str = suffix

        self._inherit_indent: bool
        self._fold_position_finder: Callable[[str, int, int], int]
        match method:
            case "brute":
                self._inherit_indent = False
                self._fold_position_finder = self._get_maximal_fold_pos
            case "simple":
                self._inherit_indent = True
                self._fold_position_finder = self._get_maximal_fold_pos
            case "smart":
                self._inherit_indent = True
                self._fold_position_finder = self._get_smart_fold_pos
            case _:
                raise FyppFatalError("invalid folding type")

    def __call__(self, line: str) -> list[str]:
        """Folds a line.

        Can be directly called to return the list of folded lines::

            linefolder = FortranLineFolder(maxlen=10)
            linefolder('  print *, "some Fortran line"')

        Args:
            line: Line to fold.

        Returns:
            Components of folded line. They should be assembled via ``\\n.join()`` to obtain the
            string representation.
        """
        if self._maxlen < 0 or len(line) <= self._maxlen:
            return [line]
        if self._inherit_indent:
            indent = len(line) - len(line.lstrip())
            prefix = " " * indent + self._prefix
        else:
            indent = 0
            prefix = self._prefix
        suffix = self._suffix
        return self._split_line(line, self._maxlen, prefix, suffix, self._fold_position_finder)

    @staticmethod
    def _split_line(
        line: str,
        maxlen: int,
        prefix: str,
        suffix: str,
        fold_position_finder: Callable[[str, int, int], int],
    ) -> list[str]:
        # length of continuation lines with 1 or two continuation chars.
        maxlen1 = maxlen - len(prefix)
        maxlen2 = maxlen1 - len(suffix)
        start = 0
        end = fold_position_finder(line, start, maxlen - len(suffix))
        result = [line[start:end] + suffix]
        while end < len(line) - maxlen1:
            start = end
            end = fold_position_finder(line, start, start + maxlen2)
            result.append(prefix + line[start:end] + suffix)
        result.append(prefix + line[end:])
        return result

    @staticmethod
    def _get_maximal_fold_pos(_: str, __: int, end: int) -> int:
        return end

    @staticmethod
    def _get_smart_fold_pos(line: str, start: int, end: int) -> int:
        linelen = end - start
        ispace = line.rfind(" ", start, end)
        # The space we waste for smart folding should be max. 1/3rd of the line
        if ispace != -1 and ispace >= start + (2 * linelen) // 3:
            return ispace
        return end


class PlaceholderLineFolder:
    """Implements a placeholder line folder returning the line unaltered."""

    def __call__(self, line: str) -> list[str]:
        """Returns the entire line without any folding.

        Returns:
            Components of folded line. They should be assembled via ``\\n.join()`` to obtain the
            string representation.
        """
        return [line]


def get_option_parser() -> optparse.OptionParser:
    """Returns an option parser for the Fypp command line tool.

    Returns:
        Parser which can create an optparse.Values object with Fypp settings based on command line
        arguments.
    """
    defs = FyppOptions()
    fypp_name = "fypp"
    fypp_desc = (
        "Preprocesses source code with Fypp directives. The input is "
        "read from INFILE (default: '-', stdin) and written to "
        "OUTFILE (default: '-', stdout)."
    )
    fypp_version = fypp_name + " " + VERSION
    usage = "%prog [options] [INFILE] [OUTFILE]"
    parser = optparse.OptionParser(
        prog=fypp_name, description=fypp_desc, version=fypp_version, usage=usage
    )

    msg = (
        "define a variable; value interpreted depending on the "
        "--define-mode option, either evaluated as a Python expression (e.g "
        "'-DDEBUG=1' sets DEBUG to the integer 1) with None as default "
        "value, or as string literal with the empty string as default value"
    )
    parser.add_option(
        "-D",
        "--define",
        action="append",
        dest="defines",
        metavar="VAR[=VALUE]",
        default=defs.defines,
        help=msg,
    )

    msg = (
        "value evaluation mode in -D options, 'eval': values are evaluated as Python "
        "expressions (default, requires quotation for string values), "
        "'str': values are string literals (no quotation needed)"
    )
    parser.add_option(
        "--define-mode",
        metavar="MODE",
        choices=["eval", "str"],
        default=defs.define_mode,
        dest="define_mode",
        help=msg,
    )

    msg = (
        "define a variable with a string literal value; similar to the -D option, but the value "
        "is always interpreted as a string literal with the empty string as default value"
    )
    parser.add_option(
        "-S",
        "--define-str",
        action="append",
        dest="defines_str",
        metavar="VAR[=VALUE]",
        default=defs.defines_str,
        help=msg,
    )

    msg = (
        "define a variable with an evaluated expression; similar to the -D option, but the "
        "value is always evaluated as a Python expression with None as default value"
    )
    parser.add_option(
        "-E",
        "--define-eval",
        action="append",
        dest="defines_eval",
        metavar="VAR[=VALUE]",
        default=defs.defines_eval,
        help=msg,
    )

    msg = "add directory to the search paths for include files"
    parser.add_option(
        "-I",
        "--include",
        action="append",
        dest="includes",
        metavar="INCDIR",
        default=defs.includes,
        help=msg,
    )

    msg = (
        "import a python module at startup (import only trustworthy modules "
        "as they have access to an **unrestricted** Python environment!)"
    )
    parser.add_option(
        "-m",
        "--module",
        action="append",
        dest="modules",
        metavar="MOD",
        default=defs.modules,
        help=msg,
    )

    msg = (
        "directory to be searched for user imported modules before "
        "looking up standard locations in sys.path"
    )
    parser.add_option(
        "-M",
        "--module-dir",
        action="append",
        dest="moduledirs",
        metavar="MODDIR",
        default=defs.moduledirs,
        help=msg,
    )

    msg = "emit line numbering markers"
    parser.add_option(
        "-n",
        "--line-numbering",
        action="store_true",
        dest="line_numbering",
        default=defs.line_numbering,
        help=msg,
    )

    msg = (
        "line numbering mode, 'full' (default): line numbering "
        "markers generated whenever source and output lines are out "
        "of sync, 'nocontlines': line numbering markers omitted "
        "for continuation lines"
    )
    parser.add_option(
        "-N",
        "--line-numbering-mode",
        metavar="MODE",
        choices=["full", "nocontlines"],
        default=defs.line_numbering_mode,
        dest="line_numbering_mode",
        help=msg,
    )

    msg = (
        "line numbering marker format,  currently 'std', 'cpp' and "
        "'gfortran5' are supported, where 'std' emits #line pragmas "
        "similar to standard tools, 'cpp' produces line directives as "
        "emitted by GNU cpp, and 'gfortran5' cpp line directives with a "
        "workaround for a bug introduced in GFortran 5. Default: 'cpp'."
    )
    parser.add_option(
        "--line-marker-format",
        metavar="FMT",
        choices=["cpp", "gfortran5", "std"],
        dest="line_marker_format",
        default=defs.line_marker_format,
        help=msg,
    )

    msg = (
        "maximal line length (default: 132), lines modified by the "
        "preprocessor are folded if becoming longer"
    )
    parser.add_option(
        "-l",
        "--line-length",
        type=int,
        metavar="LEN",
        dest="line_length",
        default=defs.line_length,
        help=msg,
    )

    msg = (
        "line folding mode, 'smart' (default): indentation context "
        "and whitespace aware, 'simple': indentation context aware, "
        "'brute': mechanical folding"
    )
    parser.add_option(
        "-f",
        "--folding-mode",
        metavar="MODE",
        choices=["smart", "simple", "brute"],
        dest="folding_mode",
        default=defs.folding_mode,
        help=msg,
    )

    msg = "suppress line folding"
    parser.add_option(
        "-F",
        "--no-folding",
        action="store_true",
        dest="no_folding",
        default=defs.no_folding,
        help=msg,
    )

    msg = "indentation to use for continuation lines (default 4)"
    parser.add_option(
        "--indentation",
        type=int,
        metavar="IND",
        dest="indentation",
        default=defs.indentation,
        help=msg,
    )

    msg = (
        "produce fixed format output (any settings for options "
        "--line-length, --folding-method and --indentation are ignored)"
    )
    parser.add_option(
        "--fixed-format",
        action="store_true",
        dest="fixed_format",
        default=defs.fixed_format,
        help=msg,
    )

    msg = (
        "character encoding for reading/writing files. Default: 'utf-8'. "
        "Note: reading from stdin and writing to stdout is encoded "
        "according to the current locale and is not affected by this setting."
    )
    parser.add_option("--encoding", metavar="ENC", default=defs.encoding, help=msg)

    msg = "create parent folders of the output file if they do not exist"
    parser.add_option(
        "-p",
        "--create-parents",
        action="store_true",
        dest="create_parent_folder",
        default=defs.create_parent_folder,
        help=msg,
    )

    msg = (
        "in variables _FILE_ and _THIS_FILE_, use relative paths with DIR "
        "as root directory. Note: the input file and all included files "
        "must be in DIR or in a directory below."
    )
    parser.add_option(
        "--file-var-root", metavar="DIR", dest="file_var_root", default=defs.file_var_root, help=msg
    )

    msg = "write a make-compatible dependency file to this location"
    parser.add_option("--depfile", metavar="DEPFILE", default=defs.depfile, help=msg)

    return parser


def run_fypp() -> None:
    """Run the Fypp command line tool."""
    options = FyppOptions()
    optparser = get_option_parser()
    # Note: parse_args returns first the same object which was passed in as values
    _, leftover = optparser.parse_args(values=options)
    infile = leftover[0] if len(leftover) > 0 else "-"
    outfile = leftover[1] if len(leftover) > 1 else "-"
    try:
        tool = Fypp(options)
        tool.process_file(infile, outfile)
    except FyppStopRequest as exc:
        sys.stderr.write(_formatted_exception(exc))
        sys.exit(USER_ERROR_EXIT_CODE)
    except FyppFatalError as exc:
        sys.stderr.write(_formatted_exception(exc))
        sys.exit(ERROR_EXIT_CODE)


def linenumdir_cpp(linenr: int, fname: str, flag: int | None = None) -> str:
    """Returns a GNU cpp style line directive.

    Args:
        linenr: Line nr (starting with zero).
        fname: File name.
        flag: Optional flag to print after the directive

    Returns:
        Line number directive as string.
    """
    if flag is None:
        return f'# {linenr + 1} "{_linemarker_path(fname)}"\n'
    return f'# {linenr + 1} "{_linemarker_path(fname)}" {flag}\n'


def linenumdir_std(linenr: int, fname: str, flag: int | None = None) -> str:
    """Returns standard #line pragma styled line directive.

    Args:
        linenr: Line nr (starting with zero).
        fname: File name.
        flag: Optional flag to print after the directive. Note, this option is only there to be API
            compatible with linenumdir_cpp(), but is ignored otherwise, since #line pragmas do not
            allow for extra file opening/closing flags.

    Returns:
        Line number directive as string.
    """
    return f'#line {linenr + 1} "{_linemarker_path(fname)}"\n'


def _normalized_path(path: str) -> str:
    """Returns path with forward slashes as separators on all platforms"""
    if os.altsep:
        path = path.replace(os.sep, os.altsep)
    return path


def _make_escaped_path(path: str) -> str:
    """Returns a path escaped for the use in make dependency rules (e.g. in dependency files).
    
    Currently only handles '$', ' ' and '#'.
    """
    return path.replace("$", "$$").replace(" ", "\\ ").replace("#", "\\#")


def _linemarker_path(path: str) -> str:
    """Returns normalized path escaped for a string literal in a line marker"""
    return _normalized_path(path).replace("\\", "\\\\").replace('"', '\\"')


def _shiftinds(inds: Sequence[int], shift: int) -> list[int]:
    return [ind + shift for ind in inds]


def _open_input_file(inpfile: str, encoding: str | None = None) -> TextIO:
    try:
        inpfp = open(inpfile, "r", encoding=encoding)
    except OSError as exc:
        msg = f"Failed to open file '{inpfile}' for read"
        raise FyppFatalError(msg) from exc
    return inpfp


def _open_output_file(
    outfile: str, encoding: str | None = None, create_parents: bool = False
) -> TextIO:
    if create_parents:
        parentdir = os.path.abspath(os.path.dirname(outfile))
        try:
            os.makedirs(parentdir, exist_ok=True)
        except OSError as exc:
            msg = f"Folder '{parentdir}' can not be created"
            raise FyppFatalError(msg) from exc
    try:
        outfp = open(outfile, "w", encoding=encoding)
    except OSError as exc:
        msg = f"Failed to open file '{outfile}' for write"
        raise FyppFatalError(msg) from exc
    return outfp


def _get_callable_argspec(func: Callable) -> _CallableArgSpec:
    sig = inspect.signature(func)
    args = []
    defaults = {}
    varpos = None
    varkw = None
    for param in sig.parameters.values():
        match param.kind:
            case param.POSITIONAL_OR_KEYWORD:
                args.append(param.name)
                if param.default is not param.empty:
                    defaults[param.name] = param.default
            case param.VAR_POSITIONAL:
                varpos = param.name
            case param.VAR_KEYWORD:
                varkw = param.name
            case _:
                msg = f"argument '{param.name}' has invalid argument type"
                raise FyppFatalError(msg)
    return _CallableArgSpec(args, defaults, varpos, varkw)


def _blank_match(match: re.Match[str]) -> str:
    size = match.end() - match.start()
    return " " * size


def _argsplit_fortran(argtxt: str) -> list[str]:
    txt = _INLINE_EVAL_REGION_REGEXP.sub(_blank_match, argtxt)
    splitpos = [-1]
    quote = None
    closing_brace_stack: list[str | None] = []
    closing_brace = None
    for ind, char in enumerate(txt):
        if quote:
            if char == quote:
                quote = None
            continue
        if char in _QUOTES_FORTRAN:
            quote = char
            continue
        if char in _OPENING_BRACKETS_FORTRAN:
            closing_brace_stack.append(closing_brace)
            ibracket = _OPENING_BRACKETS_FORTRAN.index(char)
            closing_brace = _CLOSING_BRACKETS_FORTRAN[ibracket]
            continue
        if char in _CLOSING_BRACKETS_FORTRAN:
            if char == closing_brace:
                closing_brace = closing_brace_stack.pop(-1)
                continue
            msg = (
                f"unexpected closing delimiter '{char}' in expression '{argtxt}' at position "
                f"{ind + 1}"
            )
            raise FyppFatalError(msg)
        if not closing_brace and char == _ARGUMENT_SPLIT_CHAR_FORTRAN:
            splitpos.append(ind)
    if quote or closing_brace:
        msg = f"open quotes or brackets in expression '{argtxt}'"
        raise FyppFatalError(msg)
    splitpos.append(len(txt))
    fragments = [argtxt[start + 1 : end] for start, end in itertools.pairwise(splitpos)]
    return fragments


def _formatted_exception(exc: BaseException) -> str:
    error_header_formstr = "{file}:{line}: "
    error_body_formstr = "error: {errormsg} [{errorclass}]"
    if not isinstance(exc, FyppError):
        return error_body_formstr.format(errormsg=str(exc), errorclass=exc.__class__.__name__)
    out = []
    if exc.fname is not None:
        # Guaranteed by the class contract: fname is only ever set together with span.
        assert exc.span is not None
        if exc.span.end > exc.span.start + 1:
            line = f"{exc.span.start + 1}-{exc.span.end}"
        else:
            line = f"{exc.span.start + 1}"
        out.append(error_header_formstr.format(file=exc.fname, line=line))
    out.append(error_body_formstr.format(errormsg=exc.msg, errorclass=exc.__class__.__name__))
    if exc.__cause__ is not None:
        out.append("\n" + _formatted_exception(exc.__cause__))
    out.append("\n")
    return "".join(out)


def _write_dependency_file(
    depfile: str, target: str, dependencies: Sequence[str], encoding: str, create_parents: bool
) -> None:
    """Writes dependency information into a file."""
    dependencies = [_make_escaped_path(_normalized_path(d)) for d in dependencies]
    depstr = " ".join(dependencies)
    target = _make_escaped_path(_normalized_path(target))
    with _open_output_file(depfile, encoding=encoding, create_parents=create_parents) as fobj:
        fobj.write(f"{target}: {depstr}\n")


if __name__ == "__main__":
    run_fypp()
