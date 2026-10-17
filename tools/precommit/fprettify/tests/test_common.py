#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
#    This file is part of fprettify.
#    Copyright (C) 2016-2019 Patrick Seewald, CP2K developers group
#
#    fprettify is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    fprettify is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with fprettify. If not, see <http://www.gnu.org/licenses/>.
###############################################################################
import inspect
import io
import os
import sys
import unittest

import fprettify


def joinpath(path1, path2):
    return os.path.normpath(os.path.join(path1, path2))


_MYPATH = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

# path to fprettify
RUNSCRIPT = joinpath(_MYPATH, r"../../fprettify.py")


class FprettifyTestCase(unittest.TestCase):
    """
    test class to be recognized by unittest, specialized for fprettify tests.
    """
