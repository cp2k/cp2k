.. _parameter:

Parameter file
==============

The tblite parameter files are written in `TOML <https://toml.io>`_.
Each major contribution to the energy is specified in its own table.
The presence of the table enables the contribution in the method.
The parameters are splitted in global and element specific parameters, global parameters are specified in the respective method table, while element specific parameters are collected together in the element records.
The exception are the optional element pairwise parameters for the scaling of the Hamiltonian elements which are specified in the Hamiltonian section.


TOML quickstart
---------------

.. note::

   This guide will only cover the aspects of TOML relevant for representing the parameter file in this project.
   For an extensive guide on TOML visit its `homepage <https://toml.io>`_.

We choose TOML to represent our parameter files since it produces both human and machine readable files and has good support over all major progamming languages.
TOML is a configuration file format, which allows to represent data in tables, arrays and values.
Values can have data types of *integer*, *boolean*, *real*, and *string*.

.. note::

   Integers do not cast implicitly to real numbers in TOML, make sure to actually specify a real number in the parametrization file if a real is required.
   Some TOML parsers offer an implicit conversion but this behaviour is not generally available in standard-compliant TOML parser.

Tables are specified by headers enclosed in brackets, nesting is represented by dots in the header name.

.. code-block:: toml

   [hamiltonian]
   # ...
   [hamiltonian.xtb]
   # ...
   [hamiltonian.xtb.kpair]
   # ...

The toplevel header is optional as long as the table only contains other tables, the *hamiltonian* header can be left away without changing the data structure.
Small tables can be represented by inline tables with curly braces:

.. code-block:: toml

   [dispersion]
   d4 = {sc=true, s8=2.70, a1=0.52, a2=5.00, s9=5.00}

The above structure is equivalent to

.. code-block:: toml

   [dispersion.d4]
   s6 = 1.00
   s8 = 2.70
   a1 = 0.52
   a2 = 5.00
   s9 = 5.00
   sc = true

Arrays can be created by using brackets after a key

.. code-block:: toml

   [element.C]
   shells = [ "2s", "2p" ]
   levels = [ -13.970922, -10.063292 ]
   slater = [ 2.096432, 1.80 ]
   ngauss = [ 4, 4 ]
   refocc = [ 1.0, 3.0 ]
   # ...

.. note::

   In TOML 0.5.0 arrays were constrained to equivalent data only, this restriction was lifted in TOML 1.0.0.
   In this format no mixed data arrays are required, which allows the usage of TOML 0.5.0 compliant parsers as well, if no updated library is available yet.

For libraries to handle TOML documents visit the `TOML wiki <https://github.com/toml-lang/toml/wiki>`_.



Meta section
------------

The meta section allows to specify data describing the parametrization, like the name of the method, the version of the parametrization and the format used to specify it as well as relevant publications for this parametrizations.
Only the data type of the entries is checked and whether the format version is supported by the library.
The current format version is 1, specifying a higher version will result in an error.

.. code-block:: toml
   :caption: meta section example for GFN2-xTB

   [meta]
   format = 1
   name = "GFN2-xTB"
   version = 1
   reference = "DOI: 10.1021/acs.jctc.8b01176"


Allowed entries:

=========== ========================================== =========
 Keyword     Description                                Type
=========== ========================================== =========
 format      version of the parameter file format       integer
 name        name of the parametrization                string
 version     version of the parametrization data        integer
 reference   relevant publications for this parameters  string
=========== ========================================== =========


Hamiltonian section
-------------------

The Hamiltonian section is used to declare the details of the used Hamiltonian.
Currently, only xTB type Hamiltonians can be declared in the *xtb* subtable.

.. code-block:: toml
   :caption: Hamiltonian section example for GFN2-xTB

   [hamiltonian.xtb]
   wexp = 0.5
   enscale = 2.0e-2
   cn = "gfn"
   shell = {ss=1.85, pp=2.23, dd=2.23, sd=2.0, pd=2.0}

Pair parameters can be specified for all elements of the element records.
Specifying an *X-Y* entry implies the *Y-X* entry and only one will be read and used since the pair parameters cannot be asymetric.
The default value is one for each pair.

.. code-block:: toml
   :caption: Start of the GFN1-xTB kpair table

   [hamiltonian.xtb.kpair]
   H-H = 0.96
   H-B = 0.95
   H-N = 1.04
   N-Si = 1.01
   B-P = 0.97
   Sc-Sc = 1.10
   Sc-Ti = 1.10
   # ...

.. note::

   The *kpair* section can be used to tune the Hamiltonian to better describe certain bonding situations.
   The suitable range parameters is close to one and large deviations from unity will create an instable Hamiltonian.

Allowed entries:

========= ================================ =================== =====================
 Keyword   Description                      Type                Unit
========= ================================ =================== =====================
 wexp      scaling from slater exponents    real                dimensionless
 enscale   scaling from EN difference       real                dimensionless
 cn        CN type for selfenergy shift     string
 shell     shell-specific scaling           table of reals      dimensionless
 kpair     pairwise scaling for elements    table of reals      dimensionless
========= ================================ =================== =====================


Dispersion section
------------------

The dispersion section supports the *d4* subtable for DFT-D4 type dispersion corrections and the *d3* subtable for DFT-D3 type dispersion corrections.
The rational (Becke–Johnson) damping scheme is always used.

.. code-block:: toml
   :caption: Self-consistent D4 dispersion section

   [dispersion]
   d4 = {sc=true, s8=2.70, a1=0.52, a2=5.00, s9=5.00}

Allowed entries:

========= ================================= =================== =====================
 Keyword   Description                       Type                Unit
========= ================================= =================== =====================
 s6        scaling for C6 dispersion terms   real                dimensionless
 s8        scaling for C8 dispersion terms   real                dimensionless
 a1        scaling of critical radius        real                dimensionless
 a2        offset for critical radius        real                Bohr
 s9        scaling for triple-dipole terms   real                dimensionless
 sc        use self-consistent dispersion    logical
========= ================================= =================== =====================


Repulsion section
-----------------
The xTB repulsion term can be specified in either the *gfn* or *gxtb* subtable. Only one repulsion term can be active at a time. 


GFN repulsion
~~~~~~~~~~~~~~~~~~~~~~~~~~

GFN repulsion is based on a damped Coulomb-type repulsion (*rexp*=1). The *zeff* parameter in the element records is used as the effective nuclear charge, while *arep* is used after geometric averaging in the exponent in the exponential damping function. 

.. code-block:: toml
   :caption: Repulsion for GFN1-xTB

   [repulsion]
   gfn = {kexp=1.5}

Allowed entries:

========= ================================= =================== =====================
 Keyword   Description                       Type                Unit
========= ================================= =================== =====================
 kexp      exponent for repulsion damping    real                dimensionless
 klight    exponent for light atom pairs     real                dimensionless
 rexp      exponent for inverse polynomial   real                dimensionless
========= ================================= =================== =====================


g-xTB repulsion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

g-xTB repulsion is based on a damped inverse square repulsion (*rexp*=2). The *zeff* parameter in the element records is used as the base effective nuclear charge, which is modified linearly with the CN (*rep_cn*) and the atomic CEH charge (*rep_q*). The covalent radius in the CN is defined by *cn_rcov* in the element records and the exponent of the CN term is defined by *exp_cn* in the repulsion record. The error-function-based damping uses geometrically averaged exponents in *arep* and the damping radii in *rcov_rep*. 

.. code-block:: toml
   :caption: Repulsion for g-xTB

   [repulsion]
   gxtb = {kexp=1.5}

Allowed entries:

========= ================================= =================== =====================
 Keyword   Description                       Type                Unit
========= ================================= =================== =====================
 exp_cn    exponent for repulsion CN         real                dimensionless
 rexp      exponent for inverse polynomial   real                dimensionless
========= ================================= =================== =====================

Halogen section
---------------

The GFN1-xTB specific halogen bonding correction can be specified in the *classical* subtable.

.. code-block:: toml
   :caption: Halogen bonding correction

   [halogen]
   classical = {rscale=1.3, damping=0.44}

Allowed entries:

========= ================================= =================== =====================
 Keyword   Description                       Type                Unit
========= ================================= =================== =====================
 rscale    scaling parameter for radii       real                dimensionless
 damping   damping parameter                 real                dimensionless
========= ================================= =================== =====================


Charge section
--------------

The Klopman–Ohno parametrized electrostatic model is available with the *effective* subtable, while the DFTB γ-functional electrostatic can be enabled with the *gamma* subtable.
Only one electrostatic model can be active at a time.

Klopman–Ohno electrostatic
~~~~~~~~~~~~~~~~~~~~~~~~~~

The *gam* parameter in the element records is used as atomic Hubbard parameters and scaled with the *lgam* parameter to obtain shell-resolved values.
The exponent of the kernel can be modified as well as the averaging scheme for the Hubbard parameters.
Available averaging schemes are *arithmetic* (GFN2-xTB), *harmonic* (GFN1-xTB) and *geometric*.
The electrostatic is always constructed shell-resolved.

.. code-block:: toml
   :caption: Isotropic electrostatic for GFN2-xTB

   [charge]
   effective = {gexp=2.0, average="arithmetic"}

Allowed entries:

========= ================================= =================== =====================
 Keyword   Description                       Type                Unit
========= ================================= =================== =====================
 gexp      exponent of kernel                real                dimensionless
 average   averaging scheme                  string
========= ================================= =================== =====================


DFTB γ-functional electrostatic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The *gam* parameter in the element records is used as Hubbard parameter and scaled with the *lgam* parameter to obtain shell-resolved values.
The electrostatic is always constructed shell-resolved.

.. code-block:: toml
   :caption: Isotropic electrostatic for DFTB

   [charge]
   gamma = {}


Thirdorder section
------------------

An on-site thirdorder charge is supported, to use atomic Hubbard derivatives the *shell* keyword can be set to *false*, while for shell-resolved Hubbard derivatives the scaling parameters for the respective shells have to specified.
The highest specified angular momentum is implicitly used for all but absent higher momenta.

.. code-block:: toml
   :caption: On-site shell-resolved third-order contribution

   [thirdorder]
   shell = {s=1.00, p=0.50, d=0.25}

.. important::

   While the following setup uses the atomic Hubbard derivative for all shells

   .. code:: toml

      [thirdorder]
      shell.s = 1.0

   it is fundamentally different from using an atom-resolved third-order model.


Multipole section
-----------------

The anisotropic electrostatic of GFN2-xTB can be enabled using the *damped* subtable.
It requires five parameters to setup the damping function to reduce the short-range contributions from the multipole electrostatics.

.. code-block:: toml
   :caption: Damped multipole electrostatic for GFN2-xTB

   [multipole]
   damped = {dmp3=3.0, dmp5=4.0, kexp=4.0, shift=1.2, rmax=5.0}

Allowed entries:

========= ================================= =================== =====================
 Keyword   Description                       Type                Unit
========= ================================= =================== =====================
 dmp3      damping for quadratic terms       real                dimensionless
 dmp5      damping for cubic terms           real                dimensionless
 kexp      exponent for multipole radii      real                dimensionless
 shift     shift for valence CN value        real                dimensionless
 rmax      maximum multipole radius          real                dimensionless
========= ================================= =================== =====================


Element records
---------------

The main body of the parameter file contains of element records.
The parameters here are used to initialize contributions from the tables other tables, but are collected in the element records for easy usage.
Most keywords require entries, even if the respective contribution is not used in the method.

.. note::

   Each record is identified by its symbol, which allows to have multiple parameter sets for the same element.
   Input elements which do not match any symbol, will use the parametrization of the first element record with the same atomic number.
   To ensure that the right element is used as fallback an ordered dictionary is recommended to represent the element records.

.. code-block:: toml
   :caption: Hydrogen and carbon records for GFN2-xTB

   [element.H]
   shells = [ "1s" ]
   levels = [ -10.707211 ]
   slater = [ 1.23 ]
   ngauss = [ 3 ]
   refocc = [ 1.0 ]
   kcn = [ -5.0e-2 ]
   gam = 0.405771
   lgam = [ 1.0 ]
   gam3 = 0.08
   zeff = 1.105388
   arep = 2.213717
   en = 2.20
   dkernel = 5.563889e-2
   qkernel = 2.7431e-4
   mprad = 1.4
   mpvcn = 1.0

   [element.C]
   shells = [ "2s", "2p" ]
   levels = [ -13.970922, -10.063292 ]
   slater = [ 2.096432, 1.80 ]
   ngauss = [ 4, 4 ]
   refocc = [ 1.0, 3.0 ]
   kcn = [ -1.02144e-2, 1.61657e-2 ]
   gam = 5.38015e-1
   lgam = [ 1.0, 1.1056358 ]
   gam3 = 1.50e-1
   zeff = 4.231078
   arep = 1.247655
   en = 2.55
   dkernel = -4.11674e-3
   qkernel = 2.13583e-3
   mprad = 3.0
   mpvcn = 3.0


Allowed entries:

========= ================================ =================== =====================
 Keyword   Description                      Type                Unit
========= ================================ =================== =====================
 shells    included valence shells          array of strings    dimensionless
 levels    atomic self-energies             array of reals      eV
 slater    exponents of basis functions     array of reals      1/Bohr²
 ngauss    number of STO-NG primitives      array of integers   dimensionless
 refocc    reference occupation of atom     array of reals      Unitcharge
 kcn       CN dependent self-energy shift   array of reals      eV
 shpoly    polynomial enhancement factor    array of reals      dimensionless
 gam       atomic Hubbard parameter         real                Hartree/Unitcharge²
 lgam      relative shell hardness          array of reals      dimensionless
 gam3      atomic Hubbard derivative        real                Hartree/Unitcharge³
 zeff      effective nuclear charge         real                Unitcharge
 arep      repulsion exponent               real                dimensionless
 rep_cn    linear CN dependence             real                dimensionless
 rep_q     linear CEH charge dependence     real                1/Unitcharge
 rcov_rep  damping radius repulsion         real                Bohr
 cn_rcov   radius repulsion CN              real                Bohr
 dkernel   on-site dipole kernel            real                Hartree
 qkernel   on-site quadrupole kernel        real                Hartree
 mprad     critical multipole radius        real                Bohr
 mpvcn     multipole valence CN             real                dimensionless
 xbond     halogen bonding strength         real                Hartree
 en        atomic electronegativity         real                dimensionless
========= ================================ =================== =====================
