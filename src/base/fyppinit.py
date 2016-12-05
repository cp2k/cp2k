#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Fypp initialization: Python definitions for Fypp generated code.
    Definitions can be used in all CP2K packages (in fyppinit.py and as macros in *.fypp files).
"""

def fypp_header(filename):
    """ Print fypp info. Use as ${fypp_header(_FILE_)}$ in ALL *.fypp files.
    """
    import os.path
    fname = os.path.basename(filename)
    return "!--------------------------------------------------------------------------------------------------!\n"\
        + ("! Generated from " + fname + " using Fypp.").ljust(99)                                          + "!\n"\
        + ("! **DO NOT** modify this file, edit " + fname + " instead.").ljust(99)                          + "!\n"\
        +  "!--------------------------------------------------------------------------------------------------!"

fortran_max_ndim = 7 #: maximum number of dimensions of fortran arrays
