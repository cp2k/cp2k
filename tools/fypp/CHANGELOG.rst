==========
Change Log
==========


1.2
===

Added
-----

* Allow (and promote) usage of set directive instead of setvar.

* Implement stop request via stop directive.

* Assignment to variable tuples.

* Hierarchial exception testing.


Fixed
-----

* Wrong file name in error report, when exception occurs in a macro defined in
  an included file.


1.1
===

Added
-----

* Allow inline eval and control directives in direct macro call arguments.

* Add waf integration modules.

* Examples and build system intergration chapters in user guide.

* Change log file.


1.0
===

Added
-----

* Optional suppression of line numbering in continuation lines.

* Optional creation of parent folders for output file.


Changed
-------

* Class Fypp independent of ArgumentParser.


Fixed
-----

* Fix false error, when include was within a directive.

* Wrong line number offset in eval directives.


0.12
====

Added
-----

* Implement direct call.


Changed
-------

* Remove paranthesis from direct call.


0.11
====

Added
-----

* Implement call directive.

* More precise error messages.

* Folding prevention for comment lines.

* Smart line folding, fixed format line folding.

* Python 2.7 compatibility.


Changed
-------

* Control directive prefix changed from ``@`` to ``#``.

* Rename function `default()` into `getvar()`.


Fixed
-----

* Superfluous trailing newlines in macro calls.


0.9
===

Added
-----

* Basic functionality.
