Fortran API
===========

The *dftd4* library seamlessly integrates with other Fortran projects via module interfaces,

.. note::

   Generally, all quantities used in the library are stored in `atomic units <https://en.wikipedia.org/wiki/Hartree_atomic_units>`_.

.. toctree::
   :hidden:

   Full Reference <https://dftd4.github.io/dftd4>


Handling of geometries and structure
------------------------------------

The basic infrastructure to handle molecular and periodic structures is provided by the `modular computation tool chain library <https://github.com/grimme-lab/mctc-lib>`_.
The library provides a structure type which is used to represent all geometry related informations in *dftd4*.
A structure type can be constructed from arrays or read from a file.

The constructor is provided with the generic interface ``new`` and takes an array of atomic numbers (``integer``) or element symbols (``character(len=*)``) as well as the cartesian coordinates in Bohr.
Additionally, the molecular charge and the number of unpaired electrons can be provided the ``charge`` and ``uhf`` keyword, respectively.
To create a periodic structure the lattice parameters can be passed as 3 by 3 matrix with the ``lattice`` keyword.

An example for using the constructor is given here

.. code-block:: fortran

   subroutine example
      use mctc_env, only : wp
      use mctc_io, only : structure_type, new
      implicit none
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: num(:)

      num = [6, 1, 1, 1, 1]
      xyz = reshape([ &
        &  0.00000000000000_wp, -0.00000000000000_wp,  0.00000000000000_wp, &
        & -1.19220800552211_wp,  1.19220800552211_wp,  1.19220800552211_wp, &
        &  1.19220800552211_wp, -1.19220800552211_wp,  1.19220800552211_wp, &
        & -1.19220800552211_wp, -1.19220800552211_wp, -1.19220800552211_wp, &
        &  1.19220800552211_wp,  1.19220800552211_wp, -1.19220800552211_wp],&
        & [3, size(num)])

      call new(mol, num, xyz, charge=0.0_wp, uhf=0)

      ! ...
   end subroutine example


To interact with common input file formats for structures the ``read_structure`` procedure is available.
The file type is inferred from the name of the file automatically or if a file type hint is provided directly from the enumerator of available file types.
The ``read_structure`` routine can also use an already opened unit, but in this case the file type hint is mandatory to select the correct format to read from.

.. code-block:: fortran

   subroutine example
      use mctc_env, only : error_type
      use mctc_io, only : structure_type, read_structure, file_type
      implicit none
      type(structure_type) :: mol
      type(error_type), allocatable :: error
      character(len=:), allocatable :: input

      input = "struc.xyz"

      call read_structure(mol, input, error, file_type%xyz)
      if (allocated(error)) then
         print '(a)', error%message
         stop 1
      end if

      ! ...
   end subroutine example


The structure type as well as the error type contain only allocatable members and can therefore be used without requiring explicit deconstruction.

Certain members of the structure type should be considered immutable, like the number of atoms (``nat``), the identifiers for unique atoms (``id``) and the boundary conditions (``periodic``).
To change those specific structure parameters the structure type and all dependent objects should be reconstructed to ensure a consistent setup.
Other properties, like the geometry (``xyz``), molecular charge (``charge``), number of unpaired electrons (``uhf``) and lattice parameters (``lattice``) can be changed without requiring to reconstruct dependent objects like calculators or restart data.


Error handling
--------------

The basic error handler is an allocatable derived type, available from ``mctc_env`` as ``error_type``, which signals an error by its allocation status.

.. code-block:: fortran

   use mctc_env, only : error_type, fatal_error
   implicit none
   type(error_type), allocatable :: error

   call always_ok(error)
   if (allocated(error)) then
      print '(a)', "Unexpected failure:", error%message
   end if

   call always_failed(error)
   if (allocated(error)) then
      print '(a)', "Error:", error%message
   end if

   contains
      subroutine always_ok(error)
         type(error_type), allocatable, intent(out) :: error
      end subroutine always_ok

      subroutine always_failed(error)
         type(error_type), allocatable, intent(out) :: error

         call fatal_error(error, "Message associated with this error")
      end subroutine always_failed
   end

An unhandled error might get dropped by the next procedure call.


Performing calculations
-----------------------

An example for performing a calculation with DFT-D4 is shown below.

.. tab-set::

   .. tab-item:: ≤ 3.7.0

      .. literalinclude:: ../../assets/examples/api-minimal-3.7.0/app/main.f90
         :language: fortran
         :caption: app/main.f90
         :lines: 3-7
         :dedent: 3

      .. literalinclude:: ../../assets/examples/api-minimal-3.7.0/app/main.f90
         :language: fortran
         :linenos:
         :lines: 44-75

   .. tab-item:: latest

      .. literalinclude:: ../../assets/examples/api-minimal-latest/app/main.f90
         :language: fortran
         :caption: app/main.f90
         :lines: 3-7
         :dedent: 3
         :emphasize-lines: 5

      .. literalinclude:: ../../assets/examples/api-minimal-latest/app/main.f90
         :language: fortran
         :linenos:
         :lines: 44-79
         :emphasize-lines: 1,10,11,18,28


Complete Example
----------------

A minimal program using the snippets from above could look like this:

.. tab-set::

   .. tab-item:: ≤ 3.7.0

      .. literalinclude:: ../../assets/examples/api-minimal-3.7.0/app/main.f90
         :language: fortran
         :caption: app/main.f90
         :linenos:

   .. tab-item:: latest

      .. literalinclude:: ../../assets/examples/api-minimal-latest/app/main.f90
         :language: fortran
         :caption: app/main.f90
         :linenos:
         :emphasize-lines: 7,33,34,45,54,55,62,72


The program can be compiled using the following minimal ``fpm.toml``.

.. code-block:: toml

   name = "api-minimal-3_7_0"

   [dependencies]
   dftd4.git = "https://github.com/dftd4/dftd4"
   dftd4.tag = "v3.7.0"
   multicharge.git = "https://github.com/grimme-lab/multicharge.git"
   multicharge.tag = "v0.3.0"

The examples can also be found in the `assets/examples <https://github.com/dftd4/dftd4/tree/main/assets/examples>`_ directory of the repository.
