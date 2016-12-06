=============================================
Fypp - Python powered Fortran metaprogramming
=============================================

Fypp is a Python powered preprocessor. It can be used for any programming
languages but its primary aim is to offer a Fortran preprocessor, which helps to
extend Fortran with condititional compiling and template metaprogramming
capabilities. Instead of introducing its own expression syntax, it uses Python
expressions in its preprocessor directives, offering the consistency and
versatility of Python when formulating metaprogramming tasks. It puts strong
emphasis on robustness and on neat integration into developing toolchains.

The project is `hosted on bitbucket <http://bitbucket.org/aradi/fypp>`_.

`Detailed DOCUMENTATION <http://fypp.readthedocs.org>`_ is available on
`readthedocs.org <http://fypp.readthedocs.org>`_. 

Fypp is released under the *BSD 2-clause license*.


Main features
=============

* Definition and evaluation of preprocessor variables::

    #:if DEBUG > 0
      print *, "Some debug information"
    #:endif

    #:set LOGLEVEL 2

* Macro defintions and macro calls (apart of minor syntax differences similar to
  scoped intelligent Fortran macros, which probably will once become part of the
  Fortran standard)::

    #:def assertTrue(cond)
    #:if DEBUG > 0
    if (.not. ${cond}$) then
      print *, "Assert failed in file ${_FILE_}$, line ${_LINE_}$"
      error stop
    end if
    #:endif
    #:enddef

    ! Invoked via direct call (needs no quotation)
    @:assertTrue size(myArray) > 0

    ! Invoked as Python expression (needs quotation)
    $:assertTrue('size(myArray) > 0')
    

* Conditional output::
  
    program test
    #:if defined('WITH_MPI')
      use mpi
    #:elif defined('WITH_OPENMP')
      use openmp
    #:else
      use serial
    #:endif

* Iterated output (e.g. for generating Fortran templates)::

    interface myfunc
    #:for dtype in [ 'real', 'dreal', 'complex', 'dcomplex' ]
      module procedure myfunc_${dtype}$
    #:endfor
    end interface myfunc

* Inline directives::

    logical, parameter :: hasMpi = #{if defined('MPI')}#.true.#{else}#.false.#{endif}#

* Insertion of arbitrary Python expressions::

    character(*), parameter :: comp_date = "${time.strftime('%Y-%m-%d')}$"

* Inclusion of files during preprocessing::

    #:include "macrodefs.fypp"

* Using Fortran-style continutation lines in preprocessor directives::

    #:if var1 > var2 &
        & or var2 > var4
      print *, "Doing something here"
    #:endif

* Passing multiline string arguments to macros::

    #:def debug_code(code)
      #:if DEBUG > 0
        $:code
      #:endif
    #:enddef
    
    #:call debug_code
      if (size(array) > 100) then
        print *, "DEBUG: spuriously large array"
      end if
    #:endcall

* Preprocessor comments::

    #! This will not show up in the output
    #! Also the newline characters at the end of the lines will be suppressed

* Suppressing the preprocessor output in selected regions::

    #! Definitions are read, but no output (e.g. newlines) will be produced
    #:mute
    #:include "macrodefs.fypp"
    #:endmute


Installing
==========

Fypp needs a Python interpreter of version 2.7, 3.2 or above.

Automatic install
-----------------

Use Pythons command line installer ``pip`` in order to download the stable
release from the `Fypp page on PyPI <http://pypi.python.org/pypi/fypp>`_ and
install it on your system::

  pip install fypp

This installs both, the command line tool ``fypp`` and the Python module
``fypp.py``. Latter you can import if you want to access the functionality of
Fypp directly from within your Python scripts.


Manual install
--------------

For a manual install, you can download the source code from the `Fypp project
website <http://bitbucket.org/aradi/fypp>`_ ::

  git clone https://aradi@bitbucket.org/aradi/fypp.git

The project follows `Vincent Driessens git workflow
<http://nvie.com/posts/a-successful-git-branching-model/>`_, so in order to
obtain

* the latest **stable** version, check out the `master` branch::

    cd fypp
    git co master

* the latest **development** snapshot, check out the `develop` branch::

    cd fypp
    git co develop


The command line tool is a single stand-alone script. You can run it directly
from the source folder ::
  
  FYPP_SOURCE_FOLDER/bin/fypp

or after copying it from the `bin` folder to any location listed in your `PATH`
environment variable, by just issuing ::

  fypp

The python module ``fypp.py`` can be found in ``FYP_SOURCE_FOLDER/src``.


Running
=======

The Fypp command line tool reads a file, preprocesses it and writes it to
another file, so you would typically invoke it like::

  fypp source.fpp source.f90

which would process `source.fpp` and write the result to `source.f90`.  If
input and output files are not specified, information is read from stdin and
written to stdout.

The behavior of Fypp can be influenced with various command line options. A
summary of all command line options can be obtained by::

  fypp -h
