Fortran API
===========

The *tblite* library seamlessly integrates with other Fortran projects via module interfaces,

.. note::

   Generally, all quantities used in the library are stored in `atomic units <https://en.wikipedia.org/wiki/Hartree_atomic_units>`_.

.. toctree::

   Full reference <https://tblite.github.io/tblite>


Handling of geometries and structure
------------------------------------

The basic infrastructure to handle molecular and periodic structures is provided by the `modular computation tool chain library <https://github.com/grimme-lab/mctc-lib>`_.
The library provides a structure type which is used to represent all geometry related informations in *tblite*.
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


The structure type as well as the error type are using only allocatable members and can therefore be used without requiring explicit deconstruction.

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


Calculation context
-------------------

The calculation context is available with the ``context_type`` from the ``tblite_context`` module.
The context stores error messages generated while running which can be queried using the type bound function ``failed``.
To access the actual errors the messages can be removed using the type bound subroutine ``get_error``.

An output verbosity is available in the context as the member verbosity, all procedures with access to the context will default to the verbosity of the context unless the verbosity level is overwritten by an argument.
To cutomize the output the ``context_logger`` abstract base class is available.
It must implement a type bound ``message`` procedure, which is used by the context to create output.
This type can be used to create callbacks for customizing or redirecting the output of the library.


High-level interface
--------------------

The high-level interface is defined by the calculation context, the calculator instance and its restart data.
The calculation context is defined with the ``context_type``, which stores general settings regarding the overall method independent setup of the calculation.
The actual parametrisation data is stored in the ``xtb_calculator`` type.
An instance of the calculator can be used in a thread-safe way to perform calculations for a specific structure (defined by its number of atoms, unique elements and boundary conditions).
Changing the specific structure parameters requires to reconstruct the calculator.
Finally the specific persient data for a geometry is stored in a ``wavefunction_type``, which allows to restart calculations based on previous results.
