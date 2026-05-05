C API
=====

The C API bindings are provided by using the ``iso_c_binding`` intrinsic module.
Generally, objects are exported as opaque pointers and can only be manipulated within the library.
The API user is required to delete all objects created in the library by using the provided destructor functions to avoid memory leaks.

Overall four classes of objects are provided by the library

- error handlers (``dftd4_error``),
  used to communicate exceptional conditions and errors from the library to the user
- structure containers (``dftd4_structure``),
  used to represent the system specific information and geometry data,
  only the latter are mutable for the user
- dispersion model objects (``dftd4_model``),
  general model for calculating dispersion releated properties
- damping function objects (``dftd4_damp``)
  object containing the two- and three-body damping functions
- damping function objects (``dftd4_param``)
  object to represent the actual method parametrisation

.. note::

   Generally, all quantities provided to the library are assumed to be in `atomic units <https://en.wikipedia.org/wiki/Hartree_atomic_units>`_.


Error handling
--------------

.. c:type:: struct _dftd4_error* dftd4_error;

   Error handle class

The library provides a light error handle type (``dftd4_error``) for storing error information
The error handle requires only small overhead to construct and can only contain a single error.

The handler is represented by an opaque pointer and can only be manipulated by call from the library.
The user of those objects is required to delete the handlers again using the library provided deconstructors to avoid memory leaks.

.. c:function:: dftd4_error dftd4_new_error();

   :returns: New allocation for error handle

   Create new error handle object

.. c:function:: int dftd4_check_error(dftd4_error error);

   :param error: Error handle
   :returns: Current status of error handle, non-zero in case of error

   Check error handle status

.. c:function:: void dftd4_get_error(dftd4_error error, char* buffer, const int* buffersize);

   :param error: Error handle
   :param buffer: Allocation to store error message in
   :param buffersize: Maximum length of the buffer (optional)

   Get error message from error handle

.. c:function:: void dftd4_delete_error(dftd4_error* error);

   :param error: Error handle

   Delete error handle object


Structure data
--------------

.. c:type:: struct _dftd4_structure* dftd4_structure;

   Molecular structure data class

The structure data is used to represent the system of interest in the library.
It contains immutable system specific information like the number of atoms, the unique atom groups and the boundary conditions as well as mutable geometry data like cartesian coordinates and lattice parameters.

.. c:function:: dftd4_structure dftd4_new_structure(dftd4_error error, const int natoms, const int* numbers, const double* positions, const double* charge, const double* lattice, const bool* periodic);

   :param natoms: Number of atoms in the system
   :param numbers: Atomic numbers of all atoms [natoms]
   :param positions: Cartesian coordinates in Bohr [natoms, 3]
   :param charge: Total molecular charge (optional)
   :param lattice: Lattice parameters in Bohr [3, 3] (optional)
   :param periodic: Periodic dimension of the system [3] (optional)
   :returns: New molecular structure data handle

   Create new molecular structure data (quantities in Bohr)

.. c:function:: void dftd4_delete_structure(dftd4_structure* mol);

   :param mol: Molecular structure data handle

   Delete molecular structure data

.. c:function:: void dftd4_update_structure(dftd4_error error, dftd4_structure mol, const double* positions, const double* lattice);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param positions: Cartesian coordinates in Bohr [natoms, 3]
   :param lattice: Lattice parameters in Bohr [3, 3] (optional)

   Update coordinates and lattice parameters (quantities in Bohr)


Dispersion model
----------------

.. c:type:: struct _dftd4_model* dftd4_model;

   Dispersion model class

Instantiated for a given molecular structure type, it carries no information on the geometry but relies on the atomic species of the structure object.
Recreating a structure object requires to recreate the dispersion model as well.

.. c:function:: dftd4_model dftd4_new_d4_model(dftd4_error error, dftd4_structure mol);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :returns: New dispersion model handle

   Create new D4 dispersion model

.. c:function:: dftd4_model dftd4_custom_d4_model(dftd4_error error, dftd4_structure mol, double ga, double gc, double wf);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param ga: Charge scaling height
   :param gc: Charge scaling steepness
   :param wf: Weighting factor for coordination number interpolation
   :returns: New dispersion model handle

   Create new D4 dispersion model with custom parameters

.. c:function:: dftd4_model dftd4_new_d4s_model(dftd4_error error, dftd4_structure mol);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :returns: New dispersion model handle

   Create new D4S dispersion model

.. c:function:: dftd4_model dftd4_custom_d4s_model(dftd4_error error, dftd4_structure mol, double ga, double gc);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param ga: Charge scaling height
   :param gc: Charge scaling steepness
   :returns: New dispersion model handle

   Create new D4S dispersion model with custom parameters   

.. c:function:: void dftd4_delete_model(dftd4_model* disp);

   :param disp: Dispersion model handle

   Delete dispersion model

.. c:enum:: dftd4_dispersion_model

   .. c:enumerator:: dftd4_model_d4
      Standard D4 model with charge-dependent dispersion coefficients.

   .. c:enumerator:: dftd4_model_d4s
      D4S model with charge-dependent dispersion coefficients.

   Specifies the underlying dispersion model.


Damping functions
-----------------

.. c:type:: struct _dftd4_damping* dftd4_damping;

   Damping function class

The damping function object determines the short-range behaviour of the dispersion correction.
A damping function contains a two-body and (optional) three-body damping functions, as well as a damping parameter object which contains the paramterization of both functions. 
The damping functions can be assembled in any combination and can be reused for several calculations with different structures or dispersion models.
Dispersion models have a default damping function, but can be used also with any other damping function after providing the necessary damping parameters.

.. c:function:: dftd4_damping dftd4_new_damping(dftd4_error error, int damping_2b_id, int damping_3b_id);

   :param error: Error handle
   :param damping_2b_id: ID of the two-body damping function to use
   :param damping_3b_id: ID of the three-body damping function to use (negative values deactivate three-body dispersion)
   :returns: New damping function handle

   Create new damping function with custom selection for the two- and three-body damping function.

.. c:function:: dftd4_damping dftd4_new_default_damping(dftd4_error error, dftd4_model disp);

   :param error: Error handle
   :param disp: Dispersion model handle
   :returns: New damping function handle

   Create default damping function for a specific dispersion model.

.. c:function:: void dftd4_check_params(dftd4_error error, dftd4_damping damp, dftd4_param param);

   :param error: Error handle
   :param damp: Damping function handle
   :param param: Damping parameter handle

   Check the availability of the damping parameters required for a specific damping function.

.. c:function:: void dftd4_delete_damping(dftd4_damping* damp);

   :param damp: Damping function handle

   Delete damping function

.. c:enum:: dftd4_damping_twobody

   .. c:enumerator:: dftd4_damping_twobody_rational
      Rational (Beck-Johnson) two-body damping function.

   .. c:enumerator:: dftd4_damping_twobody_screened
      Screened rational (Becke-Johnson) two-body damping function.

   .. c:enumerator:: dftd4_damping_twobody_zero
      Zero (Chai-Head-Gordon) two-body damping function.

   .. c:enumerator:: dftd4_damping_twobody_mzero
      Modified zero (Sherrill) two-body damping function.

   .. c:enumerator:: dftd4_damping_twobody_optpower
      Optimized power damping function.

   .. c:enumerator:: dftd4_damping_twobody_cso
      C-Six-Only damping function.

   .. c:enumerator:: dftd4_damping_twobody_koide
      Sphereical wave (Koide) damping function.

   Available two-body damping functions.

.. c:enum:: dftd4_damping_threebody

   .. c:enumerator:: dftd4_damping_threebody_none
      Deactivate three-body dispersion contributions.

   .. c:enumerator:: dftd4_damping_threebody_rational
      Rational (Becke-Johnson) three-body damping function.

   .. c:enumerator:: dftd4_damping_threebody_screened
      Screened rational (Becke-Johnson) three-body damping function.

   .. c:enumerator:: dftd4_damping_threebody_zero
      Zero (Chai-Head-Gordon) three-body damping function.

   .. c:enumerator:: dftd4_damping_threebody_zero_avg
      Average distance zero (Chai-Head-Gordon) three-body damping function.

   Available three-body (ATM) damping functions.


Damping parameters
------------------

.. c:type:: struct _dftd4_param* dftd4_param;

   Damping parameter class

The damping parameter object parametrizes the short-range damping of the dispersion correction.
Damping parameters are required for setting up a damping function and can be reused for several damping functions provided the appropriate parameters are set (three-body damping parameters set to zero deactivate the three-body dispersion).
Standard damping parameters are independent of the molecular structure and can easily be reused for several structures or easily exchanged.

.. c:function:: dftd4_param dftd4_new_param(double s6, double s8, double s9, double a1, double a2, double a3, double a4, double rs6, double rs8, double rs9, double alp, double bet);

   :param s6: Scaling factor for C6 contribution
   :param s8: Scaling factor for C8 contribution
   :param s9: Scaling factor for C9 contribution
   :param a1: Linear scaling factor for critical radii
   :param a2: Offset distance in Bohr for critical radii
   :param a3: Screening function scaling factor for critical radius or short-range s6 modification (used in: screened, cso)
   :param a4: Screening function exponent or short-range s6 modification (used in: screened, cso)
   :param rs6: Scaling factor for distance-radius-fraction in C6 contribution (used in: zero, mzero, koide)
   :param rs8: Scaling factor for distance-radius-fraction in C8 contribution (used in: zero, mzero, koide)
   :param rs9: Scaling factor for distance-radius-fraction in C9 contribution (used in: zero)
   :param alp: Zero-damping exponent (used in: zero, mzero)
   :param bet: Linear damping radius dependence or exponent of optimized power damping function (used in: mzero, optpower)
   :returns: New damping parameter handle

   Create new damping parameters with custom values. All parameters set to zero are deactivated.

.. c:function:: dftd4_param dftd4_load_param(dftd4_error error, char* method, int model_id, int damping_2b_id, int damping_3b_id);

   :param error: Error handle
   :param method: Name of the method to load parameters for
   :param model_id: ID of the dispersion model to load parameters for
   :param damping_2b_id: ID of the two-body damping function to load parameters for
   :param damping_3b_id: ID of the three-body damping function to load parameters for (negative values deactivate three-body dispersion)
   :returns: New damping parameter handle

   Load damping parameters a specific method/model/damping combination from internal storage

.. c:function:: dftd4_param dftd4_load_default_param(dftd4_error error, char* method, dftd4_model disp);

   :param error: Error handle
   :param method: Name of the method to load parameters for
   :param disp: Dispersion model handle
   :returns: New damping function parameter handle

   Load default damping parameters for a method with a specific model from internal storage

.. c:function:: void dftd4_delete_param(dftd4_param* param);

   :param param: Damping parameter handle

   Delete damping parameters


Calculation entrypoints
-----------------------

To evaluate dispersion energies or related properties the `dftd4_get_dispersion` procedure and similar can be used.

.. c:function:: void dftd4_get_properties(dftd4_error error, dftd4_structure mol, dftd4_model disp, double* cn, double* charges, double* c6, double* alpha);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param disp: Dispersion model handle
   :param cn: Coordination number for all atoms [natoms]
   :param charges: Partial charges for all atoms [natoms]
   :param c6: C6 coefficients for all atom pairs [natoms, natoms]
   :param alpha: Static polarizabilities for all atoms [natoms]

   Evaluate properties related to the dispersion model

.. c:function:: void dftd4_get_dispersion(dftd4_error error, dftd4_structure mol, dftd4_model disp, dftd4_damping damp, dftd4_param param, double* energy, double* gradient, double* sigma);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param disp: Dispersion model handle
   :param damp: Damping function handle
   :param param: Damping parameter handle
   :param energy: Dispersion energy
   :param gradient: Dispersion gradient [natoms, 3] (optional)
   :param sigma: Dispersion strain derivatives [3, 3] (optional)

   Evaluate the dispersion energy and its derivatives

.. c:function:: void dftd4_get_numerical_hessian(dftd4_error error, dftd4_structure mol, dftd4_model disp, dftd4_damping damp, dftd4_param param, double* hess);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param disp: Dispersion model handle
   :param damp: Damping function handle
   :param param: Damping parameter handle
   :param hess: Dispersion hessian [natoms, 3, natoms, 3]

   Evaluate the numerical hessian of the dispersion energy

.. c:function:: void dftd4_get_pairwise_dispersion(dftd4_error error, dftd4_structure mol, dftd4_model disp, dftd4_damping damp, dftd4_param param, double* pair_energy2, double* pair_energy3);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param disp: Dispersion model handle
   :param damp: Damping function handle
   :param param: Damping parameter handle
   :param energy: Pairwise additive dispersion energies
   :param energy: Pairwise non-addititive dispersion energies

   Evaluate the pairwise representation of the dispersion energy
