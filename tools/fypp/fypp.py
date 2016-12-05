#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
#
# fypp -- Python powered Fortran preprocessor
#
# Copyright (c) 2016, Bálint Aradi
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
################################################################################

'''The functionality of the Fypp preprocessor is mainly realized by
using following classes:

* `Parser`_: Reads a source file, does basic syntax checking and generates
  events.

* `Builder`_: Builds a tree representation of the source file by
  receiving events. Does additional syntax checking.

* `Renderer`_: Renders a tree built by the builder.

* `Evaluator`_: Evaluates Python expressions in a designated environment. It is
  used by `Renderer`_ when rendering eval directives.

* `Processor`_: Connects `Parser`_, `Builder`_, `Renderer`_ and `Evaluator`_
  with each other and returns for a given input the preprocessed output.

* `Fypp`_: The actual Fypp preprocessor. It initializes and drives
  `Processor`_.

* `FyppOptions`_: Contains customizable settings controling the behaviour of
  `Fypp`_. Alternatively, the function `get_option_parser()`_ can be used to
  obtain an argument parser, which can create settings based on command line
  arguments.

If any expected error occurs during processing, `FyppError`_ is raised.

Additional to the ones above, following class is used for fine tuning:

* `FortranLineFolder`_: Folds overlong lines by using Fortran continuation
  lines.

'''

from __future__ import print_function
import sys
if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins
import types
import re
import os
import errno
import time
from argparse import ArgumentParser


VERSION = '1.1'

STDIN = '<stdin>'

FILEOBJ = '<fileobj>'

STRING = '<string>'


_LINE_DIRECTIVES_PATTERN = r'''
# comment block
(?P<comment>(?:^[ \t]*\#!.*\n)+)
|
# line control directive (with optional continuation lines)
^[ \t]*\#:[ \t]*(?P<ldirective>\w+)[ \t]*
(?P<lparam>.*?(?:&[ \t]*\n[ \t]*&.*?)*)?[ \t]*\n
|
# line eval directive (with optional continuation lines)
^[ \t]*\$:[ \t]*(?P<lexpr>.*?(?:&[ \t]*\n(?:[ \t]*&)?.*?)*)[ \t]*\n
|
# direct call directive (with optional continuation lines)
^[ \t]*@:[ \t]*(?P<macro>\w+)[ \t]*
(?P<macroparams>.*?(?:&[ \t]*\n[ \t]*&.*?)*)?[ \t]*\n
'''

_INLINE_DIRECTIVES_PATTERN = r'''
# inline control directive
\#\{[ \t]*(?P<idirective>\w+)[ \t]*(?P<param>.*?)?[ \t]*\}\#
|
# inline eval directive
\$\{[ \t]*(?P<iexpr>.*?)[ \t]*\}\$
'''

_INLINE_DIRECTIVES_REGEXP = re.compile(
    _INLINE_DIRECTIVES_PATTERN, re.VERBOSE | re.MULTILINE)

_ALL_DIRECTIVES_REGEXP = re.compile(
    _LINE_DIRECTIVES_PATTERN + '|' + _INLINE_DIRECTIVES_PATTERN,
    re.VERBOSE | re.MULTILINE)

_DEF_PARAM_REGEXP = re.compile(
    r'^(?P<name>\w+)\(\s*(?P<args>(?:(?:\w+\s*,\s*)*(?:\w+)))?\s*\)$')

_CALL_PARAM_REGEXP = re.compile(r'^(?P<name>\w+)\s*$')

_SETVAR_PARAM_REGEXP = re.compile(r'^(?P<name>\w+)(?:\s+(?P<expr>.*))?\s*$')

_FOR_PARAM_REGEXP = re.compile(
    r'^(?P<loopexpr>\w+(\s*,\s*\w+)*)\s+in\s+(?P<iter>.+)\s*$')

_INCLUDE_PARAM_REGEXP = re.compile(r'^(\'|")(?P<fname>.*?)\1\s*$')

_COMMENTLINE_REGEXP = re.compile(r'^[ \t]*!.*$')

_CONTLINE_REGEXP = re.compile(r'&[ \t]*\n(?:[ \t]*&)?')

_UNESCAPE_REGEXP1 = re.compile(r'([$#])\\(\\*)([{:])')

_UNESCAPE_REGEXP2 = re.compile(r'(\})\\(\\*)([$#])')

_UNESCAPE_REGEXP3 = re.compile(r'(@)\\(\\*)([:@])')

_RESERVED_PREFIX = '__'

_RESERVED_NAMES = ('defined', 'setvar', 'getvar', '_LINE_', '_FILE_',
                   '_TIME_', '_DATE_')


class FyppError(Exception):
    '''Signalizes error occuring during preprocessing.

    Args:
        msg (str): Error message.
        fname (str): File name. None (default) if file name is not available.
        span (tuple of int): Beginning and end line of the region where error
            occured or None if not available.

    Attributes:
        msg (str): Error message.
        fname (str or None): File name or None if not available.
        span (tuple of int or None): Beginning and end line of the region
            where error occured or None if not available. Line numbers start
            from zero. For directives, which do not consume end of the line,
            start and end lines are identical.
    '''

    def __init__(self, msg, fname=None, span=None):
        super(FyppError, self).__init__()
        self.msg = msg
        self.fname = fname
        self.span = span


    def __str__(self):
        msg = [self.__class__.__name__, ':\n']
        if self.fname is not None:
            msg.append("file '" + self.fname + "'")
            if self.span is not None:
                if self.span[1] > self.span[0] + 1:
                    msg.append(', lines {}-{}'.format(
                        self.span[0] + 1, self.span[1]))
                else:
                    msg.append(', line {}'.format(self.span[0] + 1))
            msg.append('\n')
        if self.msg:
            msg.append(self.msg)
        return ''.join(msg)


class Parser:
    '''Parses a text and generates events when encountering Fypp constructs.

    Args:
        includedirs (list): List of directories, in which include files should
            be searched for, when they are not found at the default location.
    '''

    def __init__(self, includedirs=None):

        # Directories to search for include files
        if includedirs is None:
            self._includedirs = []
        else:
            self._includedirs = includedirs

        # Name of current file
        self._curfile = None

        # Directory of current file
        self._curdir = None


    def parsefile(self, fobj):
        '''Parses file or a file like object.

        Args:
            fobj (str or file): Name of a file or a file like object.
        '''
        if isinstance(fobj, str):
            if fobj == STDIN:
                self._includefile(None, sys.stdin, STDIN, os.getcwd())
            else:
                inpfp = _open_input_file(fobj)
                self._includefile(None, inpfp, fobj, os.path.dirname(fobj))
                inpfp.close()
        else:
            self._includefile(None, fobj, FILEOBJ, os.getcwd())


    def _includefile(self, span, fobj, fname, curdir):
        oldfile = self._curfile
        olddir = self._curdir
        self._curfile = fname
        self._curdir = curdir
        self.handle_include(span, fname)
        self._parse(fobj.read())
        self.handle_endinclude(span, fname)
        self._curfile = oldfile
        self._curdir = olddir


    def parse(self, txt):
        '''Parses string.

        Args:
            txt (str): Text to parse.
        '''
        self._curfile = STRING
        self._curdir = os.getcwd()
        self.handle_include(None, self._curfile)
        self._parse(txt)
        self.handle_endinclude(None, self._curfile)


    def handle_include(self, span, fname):
        '''Called when parser starts to process a new file.

        It is a dummy methond and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the include directive
                or None if called the first time for the main input.
            fname (str): Name of the file.
        '''
        self._log_event('include', span, filename=fname)


    def handle_endinclude(self, span, fname):
        '''Called when parser finished processing a file.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the include directive
                or None if called the first time for the main input.
            fname (str): Name of the file.
        '''
        self._log_event('endinclude', span, filename=fname)


    def handle_setvar(self, span, name, expr):
        '''Called when parser encounters a setvar directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
            name (str): Name of the variable.
            expr (str): String representation of the expression to be assigned
                to the variable.
        '''
        self._log_event('setvar', span, name=name, expression=expr)


    def handle_def(self, span, name, args):
        '''Called when parser encounters a def directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
            name (str): Name of the macro to be defined.
            args (list of str): Name of the macro arguments.
        '''
        self._log_event('def', span, name=name, arguments=args)


    def handle_enddef(self, span):
        '''Called when parser encounters an enddef directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._log_event('enddef', span)


    def handle_if(self, span, cond):
        '''Called when parser encounters an if directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
            cond (str): String representation of the branching condition.
        '''
        self._log_event('if', span, condition=cond)


    def handle_elif(self, span, cond):
        '''Called when parser encounters an elif directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
            cond (str): String representation of the branching condition.
        '''
        self._log_event('elif', span, condition=cond)


    def handle_else(self, span):
        '''Called when parser encounters an else directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._log_event('else', span)


    def handle_endif(self, span):
        '''Called when parser encounters an endif directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._log_event('endif', span)


    def handle_for(self, span, varexpr, iterator):
        '''Called when parser encounters a for directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
            varexpr (str): String representation of the loop variable
                expression.
            iterator (str): String representation of the iterable.
        '''
        self._log_event('for', span, variable=varexpr, iterable=iterator)


    def handle_endfor(self, span):
        '''Called when parser encounters an endfor directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._log_event('endfor', span)


    def handle_call(self, span, macro):
        '''Called when parser encounters a call directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
            macro (str): Macro to call.
        '''
        self._log_event('call', span, macro=macro)


    def handle_nextarg(self, span):
        '''Called when parser encounters a nextarg directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._log_event('nextarg', span)


    def handle_endcall(self, span):
        '''Called when parser encounters an endcall directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._log_event('endcall', span)


    def handle_eval(self, span, expr):
        '''Called when parser encounters an eval directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
            expr (str): String representation of the Python expression to
                be evaluated.
        '''
        self._log_event('eval', span, expression=expr)


    def handle_text(self, span, txt):
        '''Called when parser finds text which must left unaltered.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
            txt (str): Text.
        '''
        self._log_event('text', span, content=txt)


    def handle_comment(self, span):
        '''Called when parser finds a preprocessor comment.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._log_event('comment', span)


    def handle_mute(self, span):
        '''Called when parser finds a mute directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._log_event('mute', span)


    def handle_endmute(self, span):
        '''Called when parser finds an endmute directive.

        It is a dummy method and should be overriden for actual use.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._log_event('endmute', span)


    @staticmethod
    def _log_event(event, span=(-1, -1), **params):
        print('{}: {} --> {}'.format(event, span[0], span[1]))
        for parname, parvalue in params.items():
            print('  {}: ->|{}|<-'.format(parname, parvalue))
        print()


    def _parse(self, txt, linenr=0, linedirs=True):
        pos = 0
        if linedirs:
            regexp =  _ALL_DIRECTIVES_REGEXP
        else:
            regexp = _INLINE_DIRECTIVES_REGEXP
        for match in regexp.finditer(txt):
            groups = match.groupdict()
            start, end = match.span()
            if start > pos:
                endlinenr = linenr + txt.count('\n', pos, start)
                self._process_text(txt[pos:start], (linenr, endlinenr))
                linenr = endlinenr
            endlinenr = linenr + txt.count('\n', start, end)
            if linedirs and groups['ldirective'] is not None:
                self._process_directive(
                    groups['ldirective'], groups['lparam'],
                    (linenr, endlinenr))
            elif linedirs and groups['lexpr'] is not None:
                self._process_lexpreval(groups['lexpr'], (linenr, endlinenr))
            elif groups['idirective'] is not None:
                self._process_directive(groups['idirective'], groups['param'],
                                        (linenr, endlinenr))
            elif groups['iexpr'] is not None:
                self._process_iexpreval(groups['iexpr'], (linenr, endlinenr))
            elif linedirs and groups['comment'] is not None:
                self.handle_comment((linenr, endlinenr))
            elif linedirs and groups['macro'] is not None:
                self._process_directcall(
                    groups['macro'], groups['macroparams'], (linenr, endlinenr))
            else:
                msg = 'internal error: unknown matching pattern'
                raise FyppError(msg, self._curfile, (linenr, endlinenr))
            pos = end
            linenr = endlinenr
        if pos < len(txt):
            endlinenr = linenr + txt.count('\n', pos)
            self._process_text(txt[pos:], (linenr, endlinenr))


    def _process_text(self, txt, span):
        escaped_txt = self._unescape(txt)
        self.handle_text(span, escaped_txt)


    def _process_directive(self, directive, param, span):
        param = _CONTLINE_REGEXP.sub('', param)
        if directive == 'if':
            self.handle_if(span, param)
        elif directive == 'else':
            self._check_empty_param('else', param, span)
            self.handle_else(span)
        elif directive == 'elif':
            self.handle_elif(span, param)
        elif directive == 'endif':
            self._check_empty_param('endif', param, span)
            self.handle_endif(span)
        elif directive == 'def':
            self._process_def(param, span)
        elif directive == 'enddef':
            self._check_empty_param('enddef', param, span)
            self.handle_enddef(span)
        elif directive == 'setvar':
            self._process_setvar(param, span)
        elif directive == 'for':
            self._process_for(param, span)
        elif directive == 'endfor':
            self._check_empty_param('endfor', param, span)
            self.handle_endfor(span)
        elif directive == 'call':
            self._process_call(param, span)
        elif directive == 'nextarg':
            self._check_empty_param('nextcall', param, span)
            self.handle_nextarg(span)
        elif directive == 'endcall':
            self._check_empty_param('endcall', param, span)
            self.handle_endcall(span)
        elif directive == 'include':
            self._check_not_inline_directive('include', span)
            self._process_include(param, span)
        elif directive == 'mute':
            self._check_empty_param('mute', param, span)
            self._check_not_inline_directive('mute', span)
            self.handle_mute(span)
        elif directive == 'endmute':
            self._check_empty_param('endmute', param, span)
            self._check_not_inline_directive('endmute', span)
            self.handle_endmute(span)
        else:
            msg = "unknown directive '{}'".format(directive)
            raise FyppError(msg, self._curfile, span)


    def _process_lexpreval(self, expr, span):
        expr = _CONTLINE_REGEXP.sub('', expr)
        self.handle_eval(span, expr)


    def _process_iexpreval(self, expr, span):
        self.handle_eval(span, expr)


    def _process_directcall(self, macroname, macroparams, span):
        macroparams = _CONTLINE_REGEXP.sub('', macroparams)
        self._process_call(macroname, span)
        args = [arg.strip() for arg in macroparams.split('@@')]
        if len(args):
            linenr = span[0]
            self._parse(args[0], linenr=linenr, linedirs=False)
            for arg in args[1:]:
                self.handle_nextarg(span)
                self._parse(arg, linenr=linenr, linedirs=False)
        self.handle_endcall(span)


    def _process_def(self, param, span):
        match = _DEF_PARAM_REGEXP.search(param)
        if not match:
            msg = "invalid macro definition '{}'".format(param)
            raise FyppError(msg, self._curfile, span)
        name = match.group('name')
        argstr = match.group('args')
        if argstr is None:
            args = []
        else:
            args = [s.strip() for s in argstr.split(',')]
        self.handle_def(span, name, args)


    def _process_setvar(self, param, span):
        match = _SETVAR_PARAM_REGEXP.search(param)
        if not match:
            msg = "invalid variable assignment '{}'".format(param)
            raise FyppError(msg, self._curfile, span)
        self.handle_setvar(span, match.group('name'), match.group('expr'))


    def _process_for(self, param, span):
        match = _FOR_PARAM_REGEXP.search(param)
        if not match:
            msg = "invalid for loop declaration '{}'".format(param)
            raise FyppError(msg, self._curfile, span)
        loopexpr = match.group('loopexpr')
        loopvars = [s.strip() for s in loopexpr.split(',')]
        self.handle_for(span, loopvars, match.group('iter'))


    def _process_call(self, param, span):
        match = _CALL_PARAM_REGEXP.search(param)
        if not match:
            msg = "invalid macro call '{}'".format(param)
            raise FyppError(msg, self._curfile, span)
        name = match.group('name')
        self.handle_call(span, name)


    def _process_include(self, param, span):
        match = _INCLUDE_PARAM_REGEXP.search(param)
        if not match:
            msg = "invalid include file declaration '{}'".format(param)
            raise FyppError(msg, self._curfile, span)
        fname = match.group('fname')
        for incdir in [self._curdir] + self._includedirs:
            fpath = os.path.join(incdir, fname)
            if os.path.exists(fpath):
                break
        else:
            msg = "include file '{}' not found".format(fname)
            raise FyppError(msg, self._curfile, span)
        inpfp = _open_input_file(fpath)
        self._includefile(span, inpfp, fpath, os.path.dirname(fpath))
        inpfp.close()


    def _process_mute(self, span):
        if span[0] == span[1]:
            msg = 'Inline form of mute directive not allowed'
            raise FyppError(msg, self._curfile, span)
        self.handle_mute(span)


    def _process_endmute(self, span):
        if span[0] == span[1]:
            msg = 'Inline form of endmute directive not allowed'
            raise FyppError(msg, self._curfile, span)
        self.handle_endmute(span)


    def _check_empty_param(self, directive, param, span):
        if param.strip():
            msg = 'superfluous data in {} directive'.format(directive)
            raise FyppError(msg, self._curfile, span)

    def _check_not_inline_directive(self, directive, span):
        if span[0] == span[1]:
            msg = 'Inline form of {} directive not allowed'.format(directive)
            raise FyppError(msg, self._curfile, span)


    @staticmethod
    def _unescape(txt):
        txt = _UNESCAPE_REGEXP1.sub(r'\1\2\3', txt)
        txt = _UNESCAPE_REGEXP2.sub(r'\1\2\3', txt)
        txt = _UNESCAPE_REGEXP3.sub(r'\1\2\3', txt)
        return txt


class Builder:
    '''Builds a tree representing a text with preprocessor directives.
    '''

    def __init__(self):
        # The tree, which should be built.
        self._tree = []

        # List of all open constructs
        self._open_blocks = []

        # Nodes to which the open blocks have to be appended when closed
        self._path = []

        # Nr. of open blocks when file was opened. Used for checking whether all
        # blocks have been closed, when file processing finishes.
        self._nr_prev_blocks = []

        # Current node, to which content should be added
        self._curnode = self._tree

        # Current file
        self._curfile = None


    def reset(self):
        '''Resets the builder so that it starts to build a new tree.'''
        self._tree = []
        self._open_blocks = []
        self._path = []
        self._nr_prev_blocks = []
        self._curnode = self._tree
        self._curfile = None


    def handle_include(self, span, fname):
        '''Should be called to signalize change to new file.

        Args:
            span (tuple of int): Start and end line of the include directive
                or None if called the first time for the main input.
            fname (str): Name of the file.
        '''
        self._curfile = fname
        self._path.append(self._curnode)
        self._curnode = []
        self._open_blocks.append(('include', [span], fname, None))
        self._nr_prev_blocks.append(len(self._open_blocks))


    def handle_endinclude(self, span, fname):
        '''Should be called when processing of a file finished.

        Args:
            span (tuple of int): Start and end line of the include directive
                or None if called the first time for the main input.
            fname (str): Name of the file.
        '''
        nprev_blocks = self._nr_prev_blocks.pop(-1)
        if len(self._open_blocks) > nprev_blocks:
            directive, spans = self._open_blocks[-1][0:2]
            msg = '{} directive in line {} still unclosed when reaching end '\
                  'of file'.format(directive, spans[0][0] + 1)
            raise FyppError(msg, self._curfile)
        block = self._open_blocks.pop(-1)
        directive, spans = block[0:2]
        if directive != 'include':
            msg = 'internal error: last open block is not \'include\' when '\
                  'closing file \'{}\''.format(fname)
            raise FyppError(msg)
        if span != spans[0]:
            msg = 'internal error: span for include and endinclude differ ('\
                  '{} vs {}'.format(span, spans[0])
            raise FyppError(msg)
        oldfname, _ = block[2:4]
        if fname != oldfname:
            msg = 'internal error: mismatching file name in close_file event'\
                  " (expected: '{}', got '{}')".format(oldfname, fname)
            raise FyppError(msg, fname)
        block = directive, spans, fname, self._curnode
        self._curnode = self._path.pop(-1)
        self._curnode.append(block)


    def handle_if(self, span, cond):
        '''Should be called to signalize an if directive.

        Args:
            span (tuple of int): Start and end line of the directive.
            param (str): String representation of the branching condition.
        '''
        self._path.append(self._curnode)
        self._curnode = []
        self._open_blocks.append(('if', [span], [cond], []))


    def handle_elif(self, span, cond):
        '''Should be called to signalize an elif directive.

        Args:
            span (tuple of int): Start and end line of the directive.
            cond (str): String representation of the branching condition.
        '''
        self._check_for_open_block(span, 'elif')
        block = self._open_blocks[-1]
        directive, spans = block[0:2]
        self._check_if_matches_last(directive, 'if', spans[-1], span, 'elif')
        conds, contents = block[2:4]
        conds.append(cond)
        contents.append(self._curnode)
        spans.append(span)
        self._curnode = []


    def handle_else(self, span):
        '''Should be called to signalize an else directive.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._check_for_open_block(span, 'else')
        block = self._open_blocks[-1]
        directive, spans = block[0:2]
        self._check_if_matches_last(directive, 'if', spans[-1], span, 'else')
        conds, contents = block[2:4]
        conds.append('True')
        contents.append(self._curnode)
        spans.append(span)
        self._curnode = []


    def handle_endif(self, span):
        '''Should be called to signalize an endif directive.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._check_for_open_block(span, 'endif')
        block = self._open_blocks.pop(-1)
        directive, spans = block[0:2]
        self._check_if_matches_last(directive, 'if', spans[-1], span, 'endif')
        _, contents = block[2:4]
        contents.append(self._curnode)
        spans.append(span)
        self._curnode = self._path.pop(-1)
        self._curnode.append(block)


    def handle_for(self, span, loopvar, iterator):
        '''Should be called to signalize a for directive.

        Args:
            span (tuple of int): Start and end line of the directive.
            varexpr (str): String representation of the loop variable
                expression.
            iterator (str): String representation of the iterable.
        '''
        self._path.append(self._curnode)
        self._curnode = []
        self._open_blocks.append(('for', [span], loopvar, iterator, None))


    def handle_endfor(self, span):
        '''Should be called to signalize an endfor directive.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._check_for_open_block(span, 'endfor')
        block = self._open_blocks.pop(-1)
        directive, spans = block[0:2]
        self._check_if_matches_last(directive, 'for', spans[-1], span, 'endfor')
        loopvar, iterator, dummy = block[2:5]
        spans.append(span)
        block = (directive, spans, loopvar, iterator, self._curnode)
        self._curnode = self._path.pop(-1)
        self._curnode.append(block)


    def handle_def(self, span, name, args):
        '''Should be called to signalize a def directive.

        Args:
            span (tuple of int): Start and end line of the directive.
            name (str): Name of the macro to be defined.
            args (list of str): Name of the macro arguments.
        '''
        self._path.append(self._curnode)
        self._curnode = []
        self._open_blocks.append(('def', [span], name, args, None))


    def handle_enddef(self, span):
        '''Should be called to signalize an enddef directive.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._check_for_open_block(span, 'enddef')
        block = self._open_blocks.pop(-1)
        directive, spans = block[0:2]
        self._check_if_matches_last(directive, 'def', spans[-1], span, 'enddef')
        name, args, dummy = block[2:5]
        spans.append(span)
        block = (directive, spans, name, args, self._curnode)
        self._curnode = self._path.pop(-1)
        self._curnode.append(block)


    def handle_call(self, span, name):
        '''Should be called to signalize a call directive.

        Args:
            span (tuple of int): Start and end line of the directive.
            name (str): Name of the macro to call.
        '''
        self._path.append(self._curnode)
        self._curnode = []
        self._open_blocks.append(('call', [span], name, []))


    def handle_nextarg(self, span):
        '''Should be called to signalize a nextarg directive.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._check_for_open_block(span, 'nextarg')
        block = self._open_blocks[-1]
        directive, spans = block[0:2]
        self._check_if_matches_last(directive, 'call', spans[-1], span,
                                    'endcall')
        _, contents = block[2:4]
        contents.append(self._curnode)
        spans.append(span)
        self._curnode = []


    def handle_endcall(self, span):
        '''Should be called to signalize an endcall directive.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._check_for_open_block(span, 'endcall')
        block = self._open_blocks.pop(-1)
        directive, spans = block[0:2]
        self._check_if_matches_last(directive, 'call', spans[-1], span,
                                    'endcall')
        _, contents = block[2:4]
        contents.append(self._curnode)
        spans.append(span)
        self._curnode = self._path.pop(-1)
        self._curnode.append(block)


    def handle_setvar(self, span, name, expr):
        '''Should be called to signalize a setvar directive.

        Args:
            span (tuple of int): Start and end line of the directive.
            name (str): Name of the variable.
            expr (str): String representation of the expression to be assigned
                to the variable.
        '''
        self._curnode.append(('setvar', span, name, expr))


    def handle_eval(self, span, expr):
        '''Should be called to signalize an eval directive.

        Args:
            span (tuple of int): Start and end line of the directive.
            expr (str): String representation of the Python expression to
                be evaluated.
        '''
        self._curnode.append(('eval', span, expr))


    def handle_comment(self, span):
        '''Should be called to signalize a comment directive.

        The content of the comment is not needed by the builder, but it needs
        the span of the comment to generate proper line numbers if needed.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._curnode.append(('comment', span))


    def handle_text(self, span, txt):
        '''Should be called to pass text which goes to output unaltered.

        Args:
            span (tuple of int): Start and end line of the text.
            txt (str): Text.
        '''
        self._curnode.append(('txt', span, txt))


    def handle_mute(self, span):
        '''Should be called to signalize a mute directive.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._path.append(self._curnode)
        self._curnode = []
        self._open_blocks.append(('mute', [span], None))


    def handle_endmute(self, span):
        '''Should be called to signalize an endmute directive.

        Args:
            span (tuple of int): Start and end line of the directive.
        '''
        self._check_for_open_block(span, 'endmute')
        block = self._open_blocks.pop(-1)
        directive, spans = block[0:2]
        self._check_if_matches_last(directive, 'mute', spans[-1], span,
                                    'endmute')
        spans.append(span)
        block = (directive, spans, self._curnode)
        self._curnode = self._path.pop(-1)
        self._curnode.append(block)


    @property
    def tree(self):
        '''Returns the tree built by the Builder.'''
        return self._tree


    def _check_for_open_block(self, span, directive):
        if not len(self._open_blocks) > self._nr_prev_blocks[-1]:
            msg = 'unexpected {} directive'.format(directive)
            raise FyppError(msg, self._curfile, span)


    def _check_if_matches_last(self, lastdir, curdir, lastspan, curspan,
                               directive):
        if curdir != lastdir:
            msg = 'mismatching {} directive'.format(directive)
            raise FyppError(msg, self._curfile, curspan)
        inline_last = lastspan[0] == lastspan[1]
        inline_cur = curspan[0] == curspan[1]
        if inline_last != inline_cur:
            if inline_cur:
                msg = 'expecting line form of directive {}'.format(directive)
            else:
                msg = 'expecting inline form of directive {}'.format(directive)
            raise FyppError(msg, self._curfile, curspan)
        elif inline_cur and curspan[0] != lastspan[0]:
            msg = 'inline directives of the same construct must be in the '\
                  'same row'
            raise FyppError(msg, self._curfile, curspan)


class Renderer:

    ''''Renders a tree.

    Args:
        evaluator (Evaluator, optional): Evaluator to use when rendering eval
            directives. If None (default), Evaluator() is used.
        linenums (bool, optional): Whether linenums should be generated,
            defaults to False.
        contlinenums (bool, optional): Whether linenums for continuation
            should be generated, defaults to False.
        linefolder (callable): Callable to use when folding a line.
    '''

    def __init__(self, evaluator=None, linenums=False, contlinenums=False,
                 linefolder=None):
        # Evaluator to use for Python expressions
        self._evaluator = Evaluator() if evaluator is None else evaluator
        self._evaluator.updateenv(_DATE_=time.strftime('%Y-%m-%d'),
                                  _TIME_=time.strftime('%H:%M:%S'))

        # Name of current file being processed
        self._curfile = None

        # Number of diversions, when > 0 we are within a macro call
        self._diversions = 0

        # Whether line numbering directives should be emitted
        self._linenums = linenums

        # Whether line numbering directives in continuation lines are needed.
        self._contlinenums = contlinenums

        # Callable to be used for folding lines
        if linefolder is None:
            self._linefolder = lambda line: [line]
        else:
            self._linefolder = linefolder


    def render(self, tree, env=None):
        '''Renders a tree.

        Args:
            tree (fypp-tree): Tree to render.
            env (dict, optional): Dictionary containing additional definitions
                for the evaluator. The definitions are removed from the
                from the evaluator, once the rendering finished.

        Returns:
            str: Rendered string.
        '''
        output, eval_inds, eval_pos = self._render(tree, env)
        if eval_inds:
            self._postprocess_eval_lines(output, eval_inds, eval_pos)
        txt = ''.join(output)
        return txt


    def _render(self, tree, env=None):
        newscope = env is not None
        if newscope:
            self._evaluator.pushenv(env)
        output = []
        eval_inds = []
        eval_pos = []
        for node in tree:
            cmd = node[0]
            if cmd == 'txt':
                output.append(node[2])
            elif cmd == 'if':
                out, ieval, peval = self._get_conditional_content(*node[1:4])
                eval_inds += _shiftinds(ieval, len(output))
                eval_pos += peval
                output += out
            elif cmd == 'eval':
                out, ieval, peval = self._get_eval(*node[1:3])
                eval_inds += _shiftinds(ieval, len(output))
                eval_pos += peval
                output += out
            elif cmd == 'def':
                result = self._define_macro(*node[1:5])
                output.append(result)
            elif cmd == 'setvar':
                result = self._define_variable(*node[1:4])
                output.append(result)
            elif cmd == 'for':
                out, ieval, peval = self._get_iterated_content(*node[1:5])
                eval_inds += _shiftinds(ieval, len(output))
                eval_pos += peval
                output += out
            elif cmd == 'call':
                out, ieval, peval = self._get_called_content(*node[1:4])
                eval_inds += _shiftinds(ieval, len(output))
                eval_pos += peval
                output += out
            elif cmd == 'include':
                out, ieval, peval = self._get_included_content(*node[1:4])
                eval_inds += _shiftinds(ieval, len(output))
                eval_pos += peval
                output += out
            elif cmd == 'comment':
                output.append(self._get_comment(*node[1:2]))
            elif cmd == 'mute':
                output.append(self._get_muted_content(*node[1:3]))
            else:
                msg = "internal error: unknown command '{}'".format(cmd)
                raise FyppError(msg, self._curfile)
        if newscope:
            self._evaluator.popenv()
        return output, eval_inds, eval_pos


    def _get_eval(self, span, expr):
        self._update_linenr(span[0])
        try:
            result = self._evaluator.evaluate(expr)
        except Exception as exc:
            msg = "exception occured when evaluating '{}'\n{}".format(expr, exc)
            raise FyppError(msg, self._curfile, span)
        out = []
        ieval = []
        peval = []
        if result is not None:
            out.append(str(result))
            if not self._diversions:
                ieval.append(0)
                peval.append((span, self._curfile))
        if span[0] != span[1]:
            out.append('\n')
        return out, ieval, peval


    def _get_conditional_content(self, spans, conditions, contents):
        out = []
        ieval = []
        peval = []
        multiline = (spans[0][0] != spans[-1][1])
        for condition, content, span in zip(conditions, contents, spans):
            self._update_linenr(span[1])
            try:
                cond = bool(self._evaluator.evaluate(condition))
            except Exception as exc:
                msg = "exception occured when evaluating '{}'\n{}".format(
                    condition, exc)
                raise FyppError(msg, self._curfile, span)
            if cond:
                if self._linenums and not self._diversions and multiline:
                    out.append(linenumdir(span[1], self._curfile))
                outcont, ievalcont, pevalcont = self._render(content)
                ieval += _shiftinds(ievalcont, len(out))
                peval += pevalcont
                out += outcont
                break
        if self._linenums and not self._diversions and multiline:
            out.append(linenumdir(spans[-1][1], self._curfile))
        return out, ieval, peval


    def _get_iterated_content(self, spans, loopvars, loopiter, content):
        out = []
        ieval = []
        peval = []
        self._update_linenr(spans[0][1])
        try:
            iterobj = iter(self._evaluator.evaluate(loopiter))
        except Exception as exc:
            msg = "exception occured when evaluating '{}'\n{}"\
                .format(loopiter, exc)
            raise FyppError(msg, self._curfile, spans[0])
        multiline = (spans[0][0] != spans[-1][1])
        for var in iterobj:
            if len(loopvars) == 1:
                loopscope = {loopvars[0]: var}
            else:
                loopscope = {varname: value
                             for varname, value in zip(loopvars, var)}
            if self._linenums and not self._diversions and multiline:
                out.append(linenumdir(spans[0][1], self._curfile))
            outcont, ievalcont, pevalcont = self._render(content, loopscope)
            ieval += _shiftinds(ievalcont, len(out))
            peval += pevalcont
            out += outcont
        if self._linenums and not self._diversions and multiline:
            out.append(linenumdir(spans[1][1], self._curfile))
        return out, ieval, peval


    def _get_called_content(self, spans, name, contents):
        args = []
        self._divert()
        for content in contents:
            out = self.render(content, {})
            if len(out) and out[-1] == '\n':
                out = out[:-1]
            out_escaped = out.replace('\\', '\\\\')
            out_escaped = out_escaped.replace('"', r'\"')
            args.append('"""' + out_escaped + '"""')
        self._undivert()
        expr = "{}({})".format(name, ','.join(args))
        out, ieval, peval = self._get_eval((spans[0][0], spans[-1][1]), expr)
        return out, ieval, peval


    def _get_included_content(self, spans, fname, content):
        out = []
        oldfile = self._curfile
        self._curfile = fname
        self._evaluator.updateenv(_FILE_=self._curfile)
        if self._linenums and not self._diversions:
            out += linenumdir(0, self._curfile)
        outcont, ieval, peval = self._render(content)
        ieval = _shiftinds(ieval, len(out))
        out += outcont
        self._curfile = oldfile
        self._evaluator.updateenv(_FILE_=self._curfile)
        if self._linenums and not self._diversions and spans[0] is not None:
            out += linenumdir(spans[0][1], self._curfile)
        return out, ieval, peval


    def _define_macro(self, spans, name, args, content):
        result = ''
        try:
            macro = _Macro(name, args, content, self.render, self._divert,
                           self._undivert)
            self._evaluator.define(name, macro)
        except Exception as exc:
            msg = "exception occured when defining macro '{}'\n{}"\
                .format(name, exc)
            raise FyppError(msg, self._curfile, spans[0])
        if self._linenums and not self._diversions:
            result = linenumdir(spans[1][1], self._curfile)
        return result


    def _define_variable(self, span, name, valstr):
        result = ''
        self._update_linenr(span[0])
        try:
            self._evaluator.define(name, self._evaluator.evaluate(valstr))
        except Exception as exc:
            msg = "exception occured when setting variable {} to {}\n{}"\
                .format(name, valstr, exc)
            raise FyppError(msg, self._curfile, span)
        multiline = (span[0] != span[1])
        if self._linenums and not self._diversions and multiline:
            result = linenumdir(span[1], self._curfile)
        return result


    def _get_comment(self, span):
        if self._linenums and not self._diversions:
            return linenumdir(span[1], self._curfile)
        else:
            return ''


    def _get_muted_content(self, spans, content):
        self._render(content)
        if self._linenums and not self._diversions:
            return linenumdir(spans[-1][1], self._curfile)
        else:
            return ''


    def _divert(self):
        self._diversions += 1


    def _undivert(self):
        self._diversions -= 1
        if self._diversions < 0:
            msg = "Internal error: undivert without matching divert"
            raise FyppError(msg, self._curfile)


    def _update_linenr(self, linenr):
        if not self._diversions:
            self._evaluator.updateenv(_LINE_=linenr + 1)


    def _postprocess_eval_lines(self, output, eval_inds, eval_pos):
        ilastproc = -1
        for ieval, ind in enumerate(eval_inds):
            span, fname = eval_pos[ieval]
            if ind <= ilastproc:
                continue
            iprev, eolprev = self._find_last_eol(output, ind)
            inext, eolnext = self._find_next_eol(output, ind)
            curline = self._glue_line(output, ind, iprev, eolprev, inext,
                                      eolnext)
            output[iprev + 1:inext] = [''] * (inext - iprev - 1)
            output[ind] = self._postprocess_eval_line(curline, fname, span)
            ilastproc = inext


    @staticmethod
    def _find_last_eol(output, ind):
        'Find last newline before current position.'
        iprev = ind - 1
        while iprev >= 0:
            eolprev = output[iprev].rfind('\n')
            if eolprev != -1:
                break
            iprev -= 1
        else:
            iprev = 0
            eolprev = -1
        return iprev, eolprev


    @staticmethod
    def _find_next_eol(output, ind):
        'Find last newline before current position.'
        # find first eol after expr. evaluation
        inext = ind + 1
        while inext < len(output):
            eolnext = output[inext].find('\n')
            if eolnext != -1:
                break
            inext += 1
        else:
            inext = len(output) - 1
            eolnext = len(output[-1]) - 1
        return inext, eolnext


    @staticmethod
    def _glue_line(output, ind, iprev, eolprev, inext, eolnext):
        'Create line from parts between specified boundaries.'
        curline_parts = []
        if iprev != ind:
            curline_parts = [output[iprev][eolprev + 1:]]
            output[iprev] = output[iprev][:eolprev + 1]
        curline_parts.extend(output[iprev + 1:ind])
        curline_parts.extend(output[ind])
        curline_parts.extend(output[ind + 1:inext])
        if inext != ind:
            curline_parts.append(output[inext][:eolnext + 1])
            output[inext] = output[inext][eolnext + 1:]
        return ''.join(curline_parts)


    def _postprocess_eval_line(self, evalline, fname, span):
        lines = evalline.split('\n')
        # If line ended on '\n', last element is ''. We remove it and
        # add the trailing newline later manually.
        trailing_newline = (lines[-1] == '')
        if trailing_newline:
            del lines[-1]
        lnum = linenumdir(span[0], fname) if self._linenums else ''
        clnum = lnum if self._contlinenums else ''
        linenumsep = '\n' + lnum
        clinenumsep = '\n' + clnum
        foldedlines = [self._foldline(line) for line in lines]
        outlines = [clinenumsep.join(lines) for lines in foldedlines]
        result = linenumsep.join(outlines)
        # Add missing trailing newline
        if trailing_newline:
            trailing = '\n'
            if self._linenums:
                # Last line was folded, but no linenums were generated for
                # the continuation lines -> current line position is not
                # in sync with the one calculated from the last line number
                unsync = (
                    len(foldedlines) and len(foldedlines[-1]) > 1
                    and not self._contlinenums)
                # Eval directive in source consists of more than one line
                multiline = span[1] - span[0] > 1
                if unsync or multiline:
                    # For inline eval directives span[0] == span[1]
                    # -> next line is span[0] + 1 and not span[1] as for
                    # line eval directives
                    nextline = max(span[1], span[0] + 1)
                    trailing += linenumdir(nextline, fname)
        else:
            trailing = ''
        return result + trailing


    def _foldline(self, line):
        if _COMMENTLINE_REGEXP.match(line) is None:
            return self._linefolder(line)
        else:
            return [line]


class Evaluator:

    '''Provides an isolated environment for evaluating Python expressions.

    It can restrict the builtins which can be used within this environment
    to a (hopefully safe) subset. Additionally it defines the functions
    which are provided by the preprocessor for the eval directives.

    Note, that the restricted environment does not allow importing Python
    modules. If you need a restricted environment with modules loaded,
    launch a non-restricted one, load the modules, export its environment
    and launch a restricted one using that environment.

    Args:
        env (dict, optional): Initial definitions for the environment, defaults
            to None.
        restricted (bool, optional): Whether the restricted builtins should
            be used. Otherwise all Python builtins are accessible. Defaults to
            `True` (restricted environment.
    '''

    RESTRICTED_BUILTINS = {
        'abs': builtins.abs,
        'all': builtins.all,
        'any': builtins.any,
        'bin': builtins.bin,
        'bool': builtins.bool,
        'bytearray': builtins.bytearray,
        'bytes': builtins.bytes,
        'callable': builtins.callable,
        'chr': builtins.chr,
        'classmethod': builtins.classmethod,
        'complex': builtins.complex,
        'delattr': builtins.delattr,
        'dict': builtins.dict,
        'dir': builtins.dir,
        'divmod': builtins.divmod,
        'enumerate': builtins.enumerate,
        'filter': builtins.filter,
        'float': builtins.float,
        'format': builtins.format,
        'frozenset': builtins.frozenset,
        'getattr': builtins.getattr,
        'globals': builtins.globals,
        'hasattr': builtins.hasattr,
        'hash': builtins.hash,
        'hex': builtins.hex,
        'id': builtins.id,
        'int': builtins.int,
        'isinstance': builtins.isinstance,
        'issubclass': builtins.issubclass,
        'iter': builtins.iter,
        'len': builtins.len,
        'list': builtins.list,
        'locals': builtins.locals,
        'map': builtins.map,
        'max': builtins.max,
        'memoryview': builtins.memoryview,
        'min': builtins.min,
        'next': builtins.next,
        'object': builtins.object,
        'oct': builtins.oct,
        'ord': builtins.ord,
        'pow': builtins.pow,
        'property': builtins.property,
        'range': builtins.range,
        'repr': builtins.repr,
        'reversed': builtins.reversed,
        'round': builtins.round,
        'set': builtins.set,
        'setattr': builtins.setattr,
        'slice': builtins.slice,
        'sorted': builtins.sorted,
        'staticmethod': builtins.staticmethod,
        'str': builtins.str,
        'sum': builtins.sum,
        'super': builtins.super,
        'tuple': builtins.tuple,
        'type': builtins.type,
        'vars': builtins.vars,
        'zip': builtins.zip,
        # For Python2 True/False must be explicitely added
        'True': True,
        'False': False,
    }

    def __init__(self, env=None, restricted=True):
        # Definitions (environment) to use when evaluating expressions
        self._env = env.copy() if env is not None else {}

        # Stack for environments to implement nested scopes
        self._envstack = []

        if restricted:
            builtindict = {}
            builtindict.update(self.RESTRICTED_BUILTINS)
            builtindict['__import__'] = self._func_import
        else:
            builtindict = vars(builtins)
        builtindict['defined'] = self._func_defined
        builtindict['setvar'] = self._func_setvar
        builtindict['getvar'] = self._func_getvar

        # Permitted builtins when evaluating expressions
        self._builtins = {'__builtins__': builtindict}


    def evaluate(self, expr):
        '''Evaluate a Python expression using the `eval()` builtin.

        Args:
            expr (str): String represantion of the expression.

        Return:
            Python object: Result of the expression evaluation.
        '''
        result = eval(expr, self._builtins, self._env)
        return result


    def execute(self, code):
        '''Run Python code using the `exec()` builtin.

        Args:
            code (str): Python code to run.
        '''
        exec(code, self._builtins, self._env)


    def define(self, name, value):
        '''Define a Python entity.

        Args:
            name (str): Name of the entity.
            value (Python object): Value of the entity.

        Raises:
            FyppError: If name starts with the reserved prefix or if it is a
                reserved name.
        '''
        if name.startswith(_RESERVED_PREFIX):
            msg = "Name '{}' starts with reserved prefix '{}'"\
                .format(name, _RESERVED_PREFIX)
            raise FyppError(msg, None, None)
        if name in _RESERVED_NAMES:
            msg = "Name '{}' is reserved and can not be redefined"\
                .format(name)
            raise FyppError(msg, None, None)
        self._env[name] = value


    def updateenv(self, **vardict):
        '''Add variables to the environment.

        Args:
            **vardict: variable defintions.
        '''
        self._env.update(vardict)


    def pushenv(self, vardict):
        '''Push current environment to stack, and use its copy with additional
        new defintions instead.

        Args:
            vardict (dict): New variables.
        '''
        self._envstack.append(self._env)
        self._env = self._env.copy()
        self._env.update(vardict)


    def popenv(self):
        '''Replace current environment with pop last one from stack.'''
        self._env = self._envstack.pop(-1)


    @property
    def env(self):
        '''Return current environment.'''
        return self._env


    def _func_defined(self, var):
        return var in self._env


    def _func_import(self, name, *_, **__):
        module = self._env.get(name, None)
        if module is not None and isinstance(module, types.ModuleType):
            return module
        else:
            msg = "Import of module '{}' via '__import__' not allowed" \
                  .format(name)
            raise ImportError(msg)


    def _func_setvar(self, name, value):
        self.define(name, value)
        return ''


    def _func_getvar(self, name, defvalue):
        if name in self._env:
            return self._env[name]
        else:
            return defvalue


class _Macro:

    '''Represents a user defined macro.

    Args:
        name (str): Name of the macro.
        argnames (list of str): Macro dummy arguments.
        content (list): Content of the macro as tree.
        renderfunc (function): Function to call when content should be rendered.
            This is typically the corresponding render routine of the Builder.
        divert (function): Function to call when macro rendering started, in
            order to suppress its output. Typically the corresponding routine
            of the Builder.
        undivert (function): Function to call when macro rendering finished.
            Typically the corresponding routine of the Builder.
    '''

    def __init__(self, name, argnames, content, renderfunc, divert, undivert):
        self._name = name
        self._argnames = argnames
        self._content = content
        self._renderfunc = renderfunc
        self._divert = divert
        self._undivert = undivert


    def __call__(self, *args, **keywords):
        self._divert()
        if len(args) != len(self._argnames):
            msg = "Incorrect nr. of positional arguments for macro '{}' " \
                  "(provided: {}, needed: {})".format(
                      self._name, len(args), len(self._argnames))
            raise FyppError(msg)
        argdict = {}
        for argname, arg in zip(self._argnames, args):
            argdict[argname] = arg
        argdict.update(keywords)
        output = self._renderfunc(self._content, argdict)
        self._undivert()
        if output.endswith('\n'):
            return output[:-1]
        else:
            return output


class Processor:

    '''Connects various objects with each other to create a processor.

    Args:
        parser (Parser, optional): Parser to use for parsing text. If None
            (default), `Parser()` is used.
        builder (Builder, optional): Builder to use for building the tree
            representation of the text. If None (default), `Builder()` is used.
        renderer (Renderer, optional): Renderer to use for rendering the
            output. If None (default), `Renderer()` is used with a default
            Evaluator().
        evaluator (Evaluator, optional): Evaluator to use for evaluating Python
            expressions. If None (default), `Evaluator()` is used.
    '''

    def __init__(self, parser=None, builder=None, renderer=None,
                 evaluator=None):
        self._parser = Parser() if parser is None else parser
        self._builder = Builder() if builder is None else builder
        if renderer is None:
            evaluator = Evaluator() if evaluator is None else evaluator
            self._renderer = Renderer(evaluator)
        else:
            self._renderer = renderer

        self._parser.handle_include = self._builder.handle_include
        self._parser.handle_endinclude = self._builder.handle_endinclude
        self._parser.handle_if = self._builder.handle_if
        self._parser.handle_else = self._builder.handle_else
        self._parser.handle_elif = self._builder.handle_elif
        self._parser.handle_endif = self._builder.handle_endif
        self._parser.handle_eval = self._builder.handle_eval
        self._parser.handle_text = self._builder.handle_text
        self._parser.handle_def = self._builder.handle_def
        self._parser.handle_enddef = self._builder.handle_enddef
        self._parser.handle_setvar = self._builder.handle_setvar
        self._parser.handle_for = self._builder.handle_for
        self._parser.handle_endfor = self._builder.handle_endfor
        self._parser.handle_call = self._builder.handle_call
        self._parser.handle_nextarg = self._builder.handle_nextarg
        self._parser.handle_endcall = self._builder.handle_endcall
        self._parser.handle_comment = self._builder.handle_comment
        self._parser.handle_mute = self._builder.handle_mute
        self._parser.handle_endmute = self._builder.handle_endmute


    def process_file(self, fname, env=None):
        '''Processeses a file.

        Args:
            fname (str): Name of the file to process.
            env (dict): Additional definitons for the evaluator.

        Returns:
            str: Processed content.
        '''
        self._parser.parsefile(fname)
        return self._render(env)


    def process_text(self, txt, env=None):
        '''Processes a string.

        Args:
            txt (str): Text to process.
            env (dict): Additional definitons for the evaluator.

        Returns:
            str: Processed content.
        '''
        self._parser.parse(txt)
        return self._render(env)


    def _render(self, env):
        env = {} if env is None else env
        output = self._renderer.render(self._builder.tree, env)
        self._builder.reset()
        return ''.join(output)


def linenumdir(linenr, fname):
    '''Returns a line numbering directive.

    Args:
        linenr (int): Line nr (starting with 0).
        fname (str): File name.
    '''
    return '# {} "{}"\n'.format(linenr + 1, fname)


class Fypp:

    '''Fypp preprocessor.

    You can invoke it like ::

        tool = Fypp()
        tool.process_file('file.in', 'file.out')

    to initialize Fypp with default options, process `file.in` and write the
    result to `file.out`. If the input should be read from a string, the
    ``process_text()`` method can be used::

        tool = Fypp()
        output = tool.process_text('#:if DEBUG > 0\\nprint *, "DEBUG"\\n#:endif\\n')

    If you want to fine tune Fypps behaviour, pass a customized `FyppOptions`_
    instance at initialization::

        options = FyppOptions()
        options.fixed_format = True
        tool = Fypp(options)

    Alternatively, you can use the command line parser
    ``argparse.ArgumentParser`` to set options for Fypp. The function
    ``get_option_parser()`` returns you a default argument parser. You can then
    use its ``parse_args()`` method to obtain settings by reading the command
    line arguments::

        options = FyppOptions()
        argparser = get_option_parser()
        options = argparser.parse_args(namespace=options)
        tool = fypp.Fypp(options)

    The command line arguments can also be passed directly as a list when
    calling ``parse_args()``::

        options = FyppOptions()
        args = ['-DDEBUG=0', 'input.fpp', 'output.f90']
        argparser = get_option_parser()
        options = argparser.parse_args(args=args, namespace=options)
        tool = fypp.Fypp(options)


    Args:
        options (object): Object containing the settings for Fypp. You typically
            would pass a customized `FyppOptions`_ instance or a ``Namespace``
            object as returned by an argument parser. If not present, the
            default settings in `FyppOptions`_ are used.

    '''

    def __init__(self, options=None):
        if options is None:
            options = FyppOptions()
        inieval = Evaluator(restricted=False)
        if options.modules:
            self._import_modules(options.modules, inieval)
        if options.inifiles:
            self._exec_inifiles(options.inifiles, inieval)
        evaluator = Evaluator(env=inieval.env, restricted=True)
        if options.defines:
            self._apply_definitions(options.defines, evaluator)
        parser = Parser(options.includes)
        builder = Builder()

        fixed_format = options.fixed_format
        linefolding = not options.no_folding
        if linefolding:
            folding = 'brute' if fixed_format else options.folding_mode
            linelength = 72 if fixed_format else options.line_length
            indentation = 5 if fixed_format else options.indentation
            prefix = '&'
            suffix = '' if fixed_format else '&'
            linefolder = FortranLineFolder(linelength, indentation, folding,
                                           prefix, suffix)
        else:
            linefolder = DummyLineFolder()
        linenums = options.line_numbering
        contlinenums = (options.line_numbering_mode != 'nocontlines')
        self._create_parent_folder = options.create_parent_folder
        renderer = Renderer(
            evaluator, linenums=linenums, contlinenums=contlinenums,
            linefolder=linefolder)
        self._preprocessor = Processor(parser, builder, renderer)


    def process_file(self, infile, outfile=None, env=None):
        '''Processes input file and writes result to output file.

        Args:
            infile (str): Name of the file to read and process. If its value is
                '-', input is read from stdin.
            outfile (str, optional): Name of the file to write the result to.
                If its value is '-', result is written to stdout. If not
                present, result will be returned as string.
            env (dict, optional): Additional definitions for the evaluator.

        Returns:
            str: Result of processed input, if no outfile was specified.
        '''
        infile = STDIN if infile == '-' else infile
        output = self._preprocessor.process_file(infile, env)
        if outfile is None:
            return output
        else:
            if outfile == '-':
                outfile = sys.stdout
            else:
                outfile = _open_output_file(outfile, self._create_parent_folder)
            outfile.write(output)
            if outfile != sys.stdout:
                outfile.close()


    def process_text(self, txt, env=None):
        '''Processes a string.

        Args:
            txt (str): String to process.
            env (dict, optional): Additional definitions for the evaluator.

        Returns:
            str: Processed content.
        '''
        return self._preprocessor.process_text(txt, env)


    @staticmethod
    def _apply_definitions(defines, evaluator):
        for define in defines:
            words = define.split('=', 2)
            name = words[0]
            value = None
            if len(words) > 1:
                try:
                    value = evaluator.evaluate(words[1])
                except Exception as exc:
                    msg = "exception at evaluating '{}' in definition for " \
                          "'{}'\n{}".format(words[1], name, exc)
                    raise FyppError(msg)
            evaluator.define(name, value)


    @staticmethod
    def _import_modules(modules, evaluator):
        for module in modules:
            try:
                evaluator.execute('import ' + module)
            except Exception as ex:
                msg = "exception occured during import of module '{}'\n{}"\
                      .format(module, ex)
                raise FyppError(msg)


    @staticmethod
    def _exec_inifiles(inifiles, evaluator):
        for inifile in inifiles:
            try:
                inifp = open(inifile, 'r')
                source = inifp.read()
                inifp.close()
            except IOError as ex:
                msg = "IO error occured at reading file '{}'\n{}"\
                    .format(inifile, ex)
                raise FyppError(msg)
            try:
                code = compile(source, inifile, 'exec', dont_inherit=-1)
                evaluator.execute(code)
            except Exception as ex:
                msg = "exception occured when executing ini-file '{}'\n{}"\
                      .format(inifile, ex)
                raise FyppError(msg)


class FyppOptions:

    '''Container for Fypp options with default values.

    Attributes:
        defines (list of str): List of variable definitions in the form of
            'VARNAME=VALUE'. Default: []
        includes (list of str): List of paths to search when looking for include
            files. Default: []
        line_numbering (bool): Whether line numbering directives should appear
            in the output. Default: False
        line_numbering_mode (str): Line numbering mode 'full' or 'nocontlines'.
            Default: 'full'.
        line_length (int): Length of output lines. Default: 132.
        folding_mode (str): Folding mode 'smart', 'simple' or 'brute'. Default:
            'smart'.
        no_folding (bool): Whether folding should be suppresed. Default: False.
        indentation (int): Indentation in continuation lines. Default: 4.
        modules (list of str): Modules to import at initialization. Default: [].
        inifiles (list of str): Python files to execute at initialization.
            Default: []
        fixed_format (bool): Whether input file is in fixed format.
            Default: False.
        create_parent_folder (bool): Whether the parent folder for the output
            file should be created if it does not exist. Default: False.
    '''

    def __init__(self):
        self.defines = []
        self.includes = []
        self.line_numbering = False
        self.line_numbering_mode = 'full'
        self.line_length = 132
        self.folding_mode = 'smart'
        self.no_folding = False
        self.indentation = 4
        self.modules = []
        self.inifiles = []
        self.fixed_format = False
        self.create_parent_folder = False


class FortranLineFolder:

    '''Implements line folding with Fortran continuation lines.

    Args:
        maxlen (int, optional): Maximal line length (default: 132).
        indent (int, optional): Indentation for continuation lines (default: 4).
        method (str, optional): Folding method with following options:

            * ``brute``: folding with maximal length of continuation lines,
            * ``simple``: indents with respect of indentation of first line,
            * ``smart``: like ``simple``, but tries to fold at whitespaces.

        prefix (str, optional): String to use at the beginning of a continuation
            line (default: '&').
        suffix (str, optional): String to use at the end of the line preceeding
            a continuation line (default: '&')
    '''

    def __init__(self, maxlen=132, indent=4, method='smart', prefix='&',
                 suffix='&'):
        # Line length should be long enough that contintuation lines can host at
        # east one character apart of indentation and two continuation signs
        minmaxlen = indent + len(prefix) + len(suffix) + 1
        if maxlen < minmaxlen:
            msg = 'Maximal line length less than {} when using an indentation' \
                  'of {}'.format(minmaxlen, indent)
            raise FyppError(msg)
        self._maxlen = maxlen
        self._indent = indent
        self._prefix = ' ' * self._indent + prefix
        self._suffix = suffix
        if method not in ['brute', 'smart', 'simple']:
            raise FyppError('invalid folding type')
        if method == 'brute':
            self._inherit_indent = False
            self._fold_position_finder = self._get_maximal_fold_pos
        elif method == 'simple':
            self._inherit_indent = True
            self._fold_position_finder = self._get_maximal_fold_pos
        elif method == 'smart':
            self._inherit_indent = True
            self._fold_position_finder = self._get_smart_fold_pos


    def __call__(self, line):
        '''Folds a line.

        Can be directly called to return the list of folded lines::

            linefolder = FortranLineFolder(maxlen=10)
            linefolder('  print *, "some Fortran line"')

        Args:
            line (str): Line to fold.

        Returns:
            list of str: Components of folded line. They should be
                assembled via ``\\n.join()`` to obtain the string
                representation.
        '''
        if self._maxlen < 0 or len(line) <= self._maxlen:
            return [line]
        if self._inherit_indent:
            indent = len(line) - len(line.lstrip())
            prefix = ' ' * indent + self._prefix
        else:
            indent = 0
            prefix = self._prefix
        suffix = self._suffix
        return self._split_line(line, self._maxlen, prefix, suffix,
                                self._fold_position_finder)


    @staticmethod
    def _split_line(line, maxlen, prefix, suffix, fold_position_finder):
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
    def _get_maximal_fold_pos(_, __, end):
        return end


    @staticmethod
    def _get_smart_fold_pos(line, start, end):
        linelen = end - start
        ispace = line.rfind(' ', start, end)
        # The space we waste for smart folding should be max. 1/3rd of the line
        if ispace != -1 and ispace >= start + (2 * linelen) // 3:
            return ispace
        else:
            return end


class DummyLineFolder:

    '''Implements a dummy line folder returning the line unaltered.'''

    def __call__(self, line):
        '''Returns the entire line without any folding.

        Returns:
            list of str: Components of folded line. They should be
                assembled via ``\\n.join()`` to obtain the string
                representation.
        '''
        return [line]


def get_option_parser():
    '''Returns an option parser for the Fypp command line tool.

    Returns:
        ArgumentParser: Parser which can create a namespace object with
        Fypp settings based on command line arguments.
    '''
    fypp_name = 'fypp'
    fypp_desc = 'Preprocess source files with Fypp directives.'
    parser = ArgumentParser(prog=fypp_name, description=fypp_desc)
    msg = 'define variable, value is interpreted as ' \
          'Python expression (e.g \'-DDEBUG=1\' sets DEBUG to the ' \
          'integer 1) or set to None if ommitted'
    parser.add_argument('-D', '--define', action='append', dest='defines',
                        metavar='VAR[=VALUE]', help=msg)
    msg = 'add directory to the search paths for include files'
    parser.add_argument('-I', '--include', action='append', dest='includes',
                        metavar='INCDIR', help=msg)
    msg = 'put line numbering directives to the output'
    parser.add_argument('-n', '--line-numbering', action='store_true',
                        default=False, help=msg)
    msg = 'line numbering mode, \'full\' (default): line numbering '\
          'directives generated whenever source and output lines are out '\
          'of sync, \'nocontlines\': line numbering directives omitted '\
          'for continuation lines'
    parser.add_argument('-N', '--line-numbering-mode', metavar='MODE',
                        choices=['full', 'nocontlines'], default='full',
                        help=msg)
    msg = 'maximal line length (default: 132), lines modified by the '\
          'preprocessor are folded if becoming longer'
    parser.add_argument('-l', '--line-length', type=int, default=132,
                        metavar='LEN', help=msg)
    msg = 'line folding mode, \'smart\' (default): indentation context '\
          'and whitespace aware, \'simple\': indentation context aware, '\
          '\'brute\': mechnical folding'
    parser.add_argument('-f', '--folding-mode', metavar='MODE',
                        choices=['smart', 'simple', 'brute'],
                        default='smart', help=msg)
    msg = 'suppress line folding'
    parser.add_argument('-F', '--no-folding', action='store_true',
                        dest='no_folding', default=False, help=msg)
    msg = 'indentation to use for continuation lines (default 4)'
    parser.add_argument('--indentation', type=int, metavar='IND',
                        default=4, help=msg)
    msg = 'import python module before starting the processing'
    parser.add_argument('-m', '--module', action='append', dest='modules',
                        metavar='MOD', help=msg)
    msg = 'execute python initialization script before starting processing'
    parser.add_argument('-i', '--ini-file', action='append',
                        dest='inifiles', metavar='INI', help=msg)
    msg = 'produce fixed format output (any settings for options '\
          '--line-length, --folding-method and --indentation are ignored)'
    parser.add_argument('--fixed-format', action='store_true',
                        default=False, help=msg)
    msg = 'create parent folders of the output file if they do not exist'
    parser.add_argument('-p', '--create-parents', action='store_true',
                        default=False, dest='create_parent_folder', help=msg)
    versionstr = '%(prog)s ' + VERSION
    parser.add_argument('-v', '--version', action='version',
                        version=versionstr)
    return parser


def _add_io_arguments(parser):
    msg = "input file to be processed (default: '-', stdin)"
    parser.add_argument('infile', nargs='?', default='-', help=msg)
    msg = "output file where processed content will be written (default: " \
          "'-', stdout)"
    parser.add_argument('outfile', nargs='?', default='-', help=msg)


def run_fypp():
    '''Run the Fypp command line tool.'''
    options = FyppOptions()
    argparser = get_option_parser()
    _add_io_arguments(argparser)
    args = argparser.parse_args(namespace=options)
    try:
        tool = Fypp(args)
        tool.process_file(args.infile, args.outfile)
    except FyppError as exc:
        sys.stderr.write(str(exc))
        sys.stderr.write('\n')
        sys.exit(1)


def _shiftinds(inds, shift):
    return [ind + shift for ind in inds]


def _open_input_file(inpfile):
    try:
        inpfp = open(inpfile, 'r')
    except IOError as ex:
        msg = "Failed to open file '{}' for read\n{}".format(inpfile, ex)
        raise FyppError(msg)
    return inpfp


def _open_output_file(outfile, create_parents=False):
    if create_parents:
        parentdir = os.path.abspath(os.path.dirname(outfile))
        if not os.path.exists(parentdir):
            try:
                os.makedirs(parentdir)
            except OSError as ex:
                if ex.errno != errno.EEXIST:
                    msg = "Folder '{}' can not be created\n{}"\
                        .format(parentdir, ex)
                    raise FyppError(msg)
    try:
        outfp = open(outfile, 'w')
    except IOError as ex:
        msg = "Failed to open file '{}' for write\n{}".format(outfile, ex)
        raise FyppError(msg)
    return outfp


if __name__ == '__main__':
    run_fypp()
