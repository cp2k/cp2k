.. _c-api:

C API
=====

The C API bindings are provided by using the ``iso_c_binding`` intrinsic module.
Generally, objects are exported as opaque pointers and can only be manipulated within the library.
The API user is required delete all objects created in the library by using the provided deconstructor functions to avoid mamory leaks.

Overall five classes of objects are provided by the library

- error handlers (:c:type:`dftd3_error`),
  used to communicate exceptional conditions and errors from the library to the user
- structure containers (:c:type:`dftd3_structure`),
  used to represent the system specific information and geometry data,
  only the latter are mutable for the user
- dispersion model objects (:c:type:`dftd3_model`),
  general model for calculating dispersion releated properties
- damping function objects (:c:type:`dftd3_param`)
  polymorphic objects to represent the actual method parametrisation
- counter-poise parameter objects (:c:type:`dftd3_gcp`),
  short range correction for compensating basis set superposition errors

.. note::

   Generally, all quantities provided to the library are assumed to be in `atomic units <https://en.wikipedia.org/wiki/Hartree_atomic_units>`_.


Error handling
--------------

.. c:type:: struct _dftd3_error* dftd3_error;

   Error handle class

The library provides a light error handle type (:c:type:`dftd3_error`) for storing error information
The error handle requires only small overhead to construct and can only contain a single error.

The handler is represented by an opaque pointer and can only be manipulated by call from the library.
The user of those objects is required to delete the handlers again using the library provided deconstructors to avoid memory leaks.

.. c:function:: dftd3_error dftd3_new_error();

   :returns: New allocation for error handle

   Create new error handle object

.. c:function:: int dftd3_check_error(dftd3_error error);

   :param error: Error handle
   :returns: Current status of error handle, non-zero in case of error

   Check error handle status

.. c:function:: void dftd3_get_error(dftd3_error error, char* buffer, const int* buffersize);

   :param error: Error handle
   :param buffer: Allocation to store error message in
   :param buffersize: Maximum length of the buffer (optional)

   Get error message from error handle

.. c:function:: void dftd3_delete_error(dftd3_error* error);

   :param error: Error handle

   Delete error handle object. The handle is set to NULL after deletion.


Structure data
--------------

.. c:type:: struct _dftd3_structure* dftd3_structure;

   Molecular structure data class

The structure data is used to represent the system of interest in the library.
It contains immutable system specific information like the number of atoms, the unique atom groups and the boundary conditions as well as mutable geometry data like cartesian coordinates and lattice parameters.

.. c:function:: dftd3_structure dftd3_new_structure(dftd3_error error, const int natoms, const int* numbers, const double* positions, const double* lattice, const bool* periodic);

   :param natoms: Number of atoms in the system
   :param numbers: Atomic numbers of all atoms [natoms]
   :param positions: Cartesian coordinates in Bohr [natoms, 3]
   :param lattice: Lattice parameters in Bohr [3, 3] (optional)
   :param periodic: Periodic dimension of the system [3] (optional)
   :returns: New molecular structure data handle

   Create new molecular structure data (quantities in Bohr)

.. c:function:: void dftd3_update_structure(dftd3_error error, dftd3_structure mol, const double* positions, const double* lattice);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param positions: Cartesian coordinates in Bohr [natoms, 3]
   :param lattice: Lattice parameters in Bohr [3, 3] (optional)

   Update coordinates and lattice parameters (quantities in Bohr)

.. c:function:: void dftd3_delete_structure(dftd3_structure* mol);

   :param mol: Molecular structure data handle

   Delete molecular structure data. The handle is set to NULL after deletion.


Dispersion model
----------------

.. c:type:: struct _dftd3_model* dftd3_model;

   Dispersion model class

Instantiated for a given molecular structure type, it carries no information on the geometry but relies on the atomic species of the structure object.
Recreating a structure object requires to recreate the dispersion model as well.

.. c:function:: dftd3_model dftd3_new_d3_model(dftd3_error error, dftd3_structure mol);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :returns: New dispersion model handle

   Create new D3 dispersion model

.. c:function:: void dftd3_set_model_realspace_cutoff(dftd3_error error, dftd3_model model, double disp2, double disp3, double cn);

   :param error: Error handle
   :param model: Dispersion model handle
   :param disp2: Cutoff for two-body dispersion
   :param disp3: Cutoff for three-body dispersion
   :param cn: Cutoff for coordination number calculation

   Set realspace cutoffs for usage in the dispersion calculation

.. c:function:: void dftd3_delete_model(dftd3_model* disp);

   :param disp: Dispersion model handle

   Delete dispersion model. The handle is set to NULL after deletion.


Damping parameters
------------------

.. c:type:: struct _dftd3_param* dftd3_param;

   Damping parameter class

The damping parameter object determining the short-range behaviour of the dispersion correction.
Standard damping parameters like the rational damping are independent of the molecular structure and can easily be reused for several structures or easily exchanged.

.. c:function:: dftd3_param dftd3_new_zero_damping(dftd3_error error, double s6, double s8, double s9, double rs6, double rs8, double alp);

   :param error: Error handle
   :param s6: Scaling of induced dipole-dipole dispersion energy
   :param s8: Scaling of induced dipole-quadrupole dispersion energy
   :param s9: Scaling of induced triple-dipole dispersion energy
   :param rs6: Range-separation parameter for induced dipole-dipole dispersion energy
   :param rs8: Range-separation parameter for induced dipole-quadrupole dispersion energy
   :param alp: Exponent for the zero damping function
   :returns: New damping parameter handle

   Create new zero damping parameters

.. c:function:: dftd3_param dftd3_load_zero_damping(dftd3_error error, char* method, bool atm);

   :param error: Error handle
   :param method: Name of the method to load parameters for
   :param atm: Use three-body dispersion
   :returns: New damping parameter handle

   Load zero damping parameters from internal storage

.. c:function:: dftd3_param dftd3_new_rational_damping(dftd3_error error, double s6, double s8, double s9, double a1, double a2, double alp);

   :param error: Error handle
   :param s6: Scaling of induced dipole-dipole dispersion energy
   :param s8: Scaling of induced dipole-quadrupole dispersion energy
   :param s9: Scaling of induced triple-dipole dispersion energy
   :param a1: Scaling of atom specific critical radii
   :param a2: Constant offset of critical radii
   :param alp: Exponent for the zero damping function (used for induced triple-dipole dispersion energy)
   :returns: New damping parameter handle

   Create new rational damping parameters

.. c:function:: dftd3_param dftd3_load_rational_damping(dftd3_error error, char* method, bool atm);

   :param error: Error handle
   :param method: Name of the method to load parameters for
   :param atm: Use three-body dispersion
   :returns: New damping parameter handle

   Load rational damping parameters from internal storage

.. c:function:: dftd3_param dftd3_new_mzero_damping(dftd3_error error, double s6, double s8, double s9, double rs6, double rs8, double alp, double bet);

   :param error: Error handle
   :param s6: Scaling of induced dipole-dipole dispersion energy
   :param s8: Scaling of induced dipole-quadrupole dispersion energy
   :param s9: Scaling of induced triple-dipole dispersion energy
   :param rs6: Range-separation parameter for induced dipole-dipole dispersion energy
   :param rs8: Range-separation parameter for induced dipole-quadrupole dispersion energy
   :param alp: Exponent for the zero damping function
   :returns: New damping parameter handle

   Create new modified zero damping parameters

.. c:function:: dftd3_param dftd3_load_mzero_damping(dftd3_error error, char* method, bool atm);

   :param error: Error handle
   :param method: Name of the method to load parameters for
   :param atm: Use three-body dispersion
   :returns: New damping parameter handle

   Load modified zero damping parameters from internal storage

.. c:function:: dftd3_param dftd3_new_mrational_damping(dftd3_error error, double s6, double s8, double s9, double a1, double a2, double alp);

   :param error: Error handle
   :param s6: Scaling of induced dipole-dipole dispersion energy
   :param s8: Scaling of induced dipole-quadrupole dispersion energy
   :param s9: Scaling of induced triple-dipole dispersion energy
   :param a1: Scaling of atom specific critical radii
   :param a2: Constant offset of critical radii
   :param alp: Exponent for the zero damping function (used for induced triple-dipole dispersion energy)
   :returns: New damping parameter handle

   Create new modified rational damping parameters

.. c:function:: dftd3_param dftd3_load_mrational_damping(dftd3_error error, char* method, bool atm);

   :param error: Error handle
   :param method: Name of the method to load parameters for
   :param atm: Use three-body dispersion
   :returns: New damping parameter handle

   Load modified rational damping parameters from internal storage

.. c:function:: dftd3_param dftd3_new_optimizedpower_damping(dftd3_error error, double s6, double s8, double s9, double a1, double a2, double alp, double bet);

   :param error: Error handle
   :param s6: Scaling of induced dipole-dipole dispersion energy
   :param s8: Scaling of induced dipole-quadrupole dispersion energy
   :param s9: Scaling of induced triple-dipole dispersion energy
   :param a1: Scaling of atom specific critical radii
   :param a2: Constant offset of critical radii
   :param alp: Exponent for the zero damping function
   :param bet: Exponent for the rational damping function (used for induced triple-dipole dispersion energy)
   :returns: New damping parameter handle

   Create new optimized power damping parameters

.. c:function:: dftd3_param dftd3_load_optimizedpower_damping(dftd3_error error, char* method, bool atm);

   :param error: Error handle
   :param method: Name of the method to load parameters for
   :param atm: Use three-body dispersion
   :returns: New damping parameter handle

   Load optimized power damping parameters from internal storage

.. c:function:: dftd3_param dftd3_new_cso_damping(dftd3_error error, double s6, double s9, double a1, double a2, double a3, double a4, double alp);

   :param error: Error handle
   :param s6: Scaling of induced dipole-dipole dispersion energy
   :param s9: Scaling of induced triple-dipole dispersion energy
   :param a1: Sigmoid amplitude parameter
   :param a2: Sigmoid reference distance scale
   :param a3: Denominator critical radii scale
   :param a4: Denominator constant offset
   :param alp: Exponent for the zero damping function (used for induced triple-dipole dispersion energy)
   :returns: New damping parameter handle

   Create new CSO (C6-scaled only) damping parameters

.. c:function:: dftd3_param dftd3_load_cso_damping(dftd3_error error, char* method, bool atm);

   :param error: Error handle
   :param method: Name of the method to load parameters for
   :param atm: Use three-body dispersion
   :returns: New damping parameter handle

   Load CSO damping parameters from internal storage

.. c:function:: void dftd3_delete_param(dftd3_param* param);

   :param param: Dispersion parameter handle

   Delete damping parameters. The handle is set to NULL after deletion.


Geometrical counter-poise correction
------------------------------------

.. c:type:: struct _dftd3_gcp* dftd3_gcp;

   Counter-poise parameter class

The counter-poise parameter object provides an additional short ranged correction to account for basis set superposition error in small basis sets.

.. c:function:: dftd3_gcp dftd3_load_gcp_param(dftd3_error error, dftd3_structure mol, char* method, char* basis);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param method: Name of the method to load parameters for
   :param basis: Name of the basis to load parameters for
   :returns: New counter-poise parameter handle

   Load geometrical counter-poise parameters from internal storage

.. c:function:: void dftd3_set_gcp_realspace_cutoff(dftd3_error error, dftd3_gcp gcp, double bas, double srb);

   :param error: Error handle
   :param model: Dispersion model handle
   :param bas: Cutoff for basis set superposition correction
   :param srb: Cutoff for short-range bond correction

   Set realspace cutoffs for usage in the counter-poise calculation

.. c:function:: void dftd3_delete_gcp(dftd3_gcp* gcp);

   :param param: Counter-poise parameter handle

   Delete counter-poise parameters. The handle is set to NULL after deletion.


Calculation entrypoints
-----------------------

To evaluate dispersion energies or related properties the :c:func:`dftd3_get_dispersion` procedure and similar can be used.

.. c:function:: void dftd3_get_dispersion(dftd3_error error, dftd3_structure mol, dftd3_model disp, dftd3_param param, double* energy, double* gradient, double* sigma);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param disp: Dispersion model handle
   :param param: Damping function parameter handle
   :param energy: Dispersion energy
   :param gradient: Dispersion gradient [natoms, 3] (optional)
   :param sigma: Dispersion strain derivatives [3, 3] (optional)

   Evaluate the dispersion energy and its derivatives.

.. c:function:: void dftd3_get_pairwise_dispersion(dftd3_error error, dftd3_structure mol, dftd3_model disp, dftd3_param param, double* pair_energy2, double* pair_energy3);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param disp: Dispersion model handle
   :param param: Damping function parameter handle
   :param energy2: Pairwise additive dispersion energies
   :param energy3: Pairwise non-addititive dispersion energies

   Evaluate the pairwise representation of the dispersion energy

.. c:function:: void dftd3_get_counterpoise(dftd3_error error, dftd3_structure mol, dftd3_gcp gcp, double* energy, double* gradient, double* sigma);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param gcp: Counter-poise parameter handle
   :param energy: Dispersion energy
   :param gradient: Dispersion gradient [natoms, 3] (optional)
   :param sigma: Dispersion strain derivatives [3, 3] (optional)

   Evaluate the counter-poise energy and its derivatives.


Memory management
-----------------

For each object type, a deconstructor function is provided to free the memory allocated by the library.
A type-generic macro is provided to select the correct deconstructor based on the object type.
Note that NULL pointers are allowed and will be ignored by the deconstructor.


.. c:macro:: dftd3_delete(ptr)

   :param ptr: Object handle (e.g., `dftd3_error`, `dftd3_model`, etc.)

   Macro to delete objects created by the library. The handle is set to NULL after deletion. The macro is type-generic and selects the correct deconstructor based on the object type:

   - :c:func:`dftd3_delete_error`
   - :c:func:`dftd3_delete_structure`
   - :c:func:`dftd3_delete_model`
   - :c:func:`dftd3_delete_param`
   - :c:func:`dftd3_delete_gcp`


Performing calculations
-----------------------

An example wrapper to perform a DFT-D3(BJ)-ATM calculation is shown below.


.. code-block:: c

   #include <stdbool.h>
   #include <stdio.h>
   #include <stdlib.h>

   #include "dftd3.h"

   static const buffersize = 512;

   int
   calc_dftd3(int natoms, int* numbers, double* positions,
              double* lattice, bool* periodic, char* method,
              double* energy, double* gradient, double* sigma)
   {
     // Local API objects from the s-dftd3 library
     dftd3_error error = dftd3_new_error();
     dftd3_structure mol = NULL;
     dftd3_model disp = NULL;
     dftd3_param param = NULL;
     int stat = EXIT_SUCCESS;

     // Create a new geometry for the library to work with
     mol = dftd3_new_structure(error, natoms, numbers, positions, lattice, periodic);
     stat = dftd3_check_error(error);

     if (stat) {
       // Initialize the D3 dispersion model for the given structure,
       // this step depends on the atomic numbers, but not on the actual geometry
       disp = dftd3_new_d3_model(error, mol);
       stat = dftd3_check_error(error);
     }

     if (stat) {
       // Load D3(BJ)-ATM parameters for the given method from internal storage,
       // this step depends on the atomic numbers, but not on the actual geometry
       param = dftd3_load_rational_damping(error, mol, method, true);
       stat = dftd3_check_error(error);
     }

     if (stat) {
       // Evaluate the dispersion energy, gradient and virial,
       // the gradient and virial are optional and can be replaced by NULL
       dftd3_get_dispersion(error, mol, disp, param, &energy, gradient, sigma);
       stat = dftd3_check_error(error);
     }

     if (!stat) {
       char buffer[buffersize];
       dftd3_get_error(error, buffer, &buffersize);
       printf("[Error] %s\n", buffer);
     }

     // Always free the used memory
     dftd3_delete(error);
     dftd3_delete(mol);
     dftd3_delete(disp);
     dftd3_delete(param);

     return stat;
   }
