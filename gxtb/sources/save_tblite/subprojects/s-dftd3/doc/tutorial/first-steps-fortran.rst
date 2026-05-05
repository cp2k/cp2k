Testing D3 damping parameters
=============================

In this tutorial we are implementing a command line tool to compute the dispersion energy for a reaction and test multiple damping functions for D3.
For this we will use the Fortran API via the ``dftd3`` module.


Using the Fortran API
---------------------

To make the ``dftd3`` available in our project we will use the Fortran package manager (`fpm <https://fpm.fortran-lang.org>`_) to manage the necessary dependencies.
We create a minimal package manifest with the following content:

.. code-block:: toml
   :caption: fpm.toml

   name = "param-scanner"
   version = "1.0.0"

   [dependencies]
   s-dftd3.git = "https://github.com/dftd3/simple-dftd3"

Our library will be called *param-scanner* because it will be scanning through the available damping parameters for D3.
We also need to declare a version, here we go with *1.0.0* but feel free to choose any version which feels appropriate.
Finally, we declare our dependencies via the URL of the git repository.


Computing dispersion for reactions
----------------------------------

For our parameter scanner to give meaningful results, we would like to target relative energies or reaction energies.
We start by defining a type to hold the reaction definition, combining an array of ``structure_type`` values and an array of stochiometric coefficients.
The ``structure_type`` value is not defined by the ``dftd3`` module but the ``mctc_io`` module, which we already include as a transient dependency.

.. note::

   For more details on ``mctc_io`` checkout its `module documentation <https://grimme-lab.github.io/mctc-lib/module/mctc_io.html>`__.

We setup our ``reaction_type`` as part of our main module ``d3_param_scan``:

.. literalinclude:: first-steps-fortran/src/scan.f90
   :language: fortran
   :caption: src/scan.f90
   :lines: 1-7,9-13,162-163

.. note::

   Here we do not define the precision ourselves but use the one defined by ``mctc_env``.

Now can define a function to evaluate the dispersion energy for a reaction.
For this we need to use the ``get_dispersion`` subroutine provided in the ``dftd3`` module.
We start with checking the signature of this function:

.. code-block:: fortran

   interface
      !> Calculate scalar dispersion energy.
      subroutine get_dispersion(mol, disp, param, cutoff, energy, gradient, sigma)
         use mctc_env, only : wp
         use mctc_io, only : structure_type
         use dftd3, only : d3_model, damping_param, realspace_cutoff
         !> Molecular structure data
         class(structure_type), intent(in) :: mol
         !> Dispersion model
         class(d3_model), intent(in) :: disp
         !> Damping parameters
         class(damping_param), intent(in) :: param
         !> Realspace cutoffs
         type(realspace_cutoff), intent(in) :: cutoff
         !> Dispersion energy
         real(wp), intent(out) :: energy
         !> Dispersion gradient
         real(wp), intent(out), contiguous, optional :: gradient(:, :)
         !> Dispersion virial
         real(wp), intent(out), contiguous, optional :: sigma(:, :)
      end subroutine get_dispersion
   end interface

The main inputs we need to construct are the ``d3_model``, the ``damping_param`` and the ``realspace_cutoff`` objects.
We check the documentation and find that the ``d3_model`` can be constructed from a ``structure_type``.
Similarly, for the ``realspace_cutoff`` we find that the derived type is rather simple and can be easily constructed.
For the ``damping_param`` we find that there are several implementations available.
Our goal is to scan through the damping functions, therefore we will have a closer look to ``damping_param`` soon.
First, we use ``get_dispersion`` to define our own dispersion energy evaluator for our reaction values.

.. literalinclude:: first-steps-fortran/src/scan.f90
   :language: fortran
   :caption: src/scan.f90
   :lines: 17-35

We keep the interface for our ``get_dispersion_for_reaction`` function simple, we focus on the two inputs, the reaction data and the damping parameters, and obtain a single output, the dispersion energy.
For the ``energy`` dummy argument we choose *intent(inout)* in contrast to the ``get_dispersion`` interface where it is *intent(out)*.


Creating a parameter scanner
----------------------------

Next we will look into actually iterating through all the available damping schemes available for D3.
Before we start we define the interface for our procedure:

.. code-block:: fortran

   interface
      !> Scan all available damping parameters for a given method by computing
      !> and comparing the dispersion energy of a reaction.
      subroutine scan_param_for_reaction(error, reaction, method, dft_energy)
         import :: reaction_type
         use mctc_env, only : wp, error_type
         !> Error handling
         type(error_type), allocatable :: error
         !> Reaction data for relative energy computation
         type(reaction_type), intent(in) :: reaction
         !> Method name to query for damping parameters
         character(*), intent(in) :: method
         !> Optional DFT energy of the reaction
         real(wp), intent(in), optional :: dft_energy
      end subroutine scan_param_for_reaction
   end interface

For our interface we will be using the ``error_type`` provided by the ``mctc_env`` module.
In a larger code base we would probably have our own error handling strategy, for this tutorial however we will reuse the same as used in ``dftd3`` to keep the code simple.
Also, we will not return any output back to the calling routine here but rather directly print them out.

With this design decision made we can now have a closer look into the damping schemes and parameters again.
We will start with the zero damping scheme D3 was originally published with, overall there are five different damping schemes currently supported in ``dftd3``.
For this we need to create an instance of the ``zero_damping_param``, which is an implementation of ``damping_param``.
Based on the just defined interface we implement our first block to obtain zero damping parameters for the input method:

.. literalinclude:: first-steps-fortran/src/scan.f90
   :language: fortran
   :caption: src/scan.f90
   :lines: 37-49,52-71,161

For handling the availability of the parameters we use that the ``get_zero_damping`` function fails with populating the error handler.
In case we find the error handler is allocated we leave the block, but we clear the error by deallocating it, the success of loading the parameters will be tracked by the ``has_param`` array.
The dispersion energy computed by our ``get_dispersion_for_reaction`` function will be stored in ``disp_energies``.
For ``disp_energies`` we initialize with the DFT energy if present or zero otherwise and use the fact that the ``energy`` dummy argument will be updated rather than overridden.

The damping parameters are retrieved in a simple derived type ``d3_param``, which we can use to initialize the actually ``damping_param`` implementation.
For ``zero_damping_param`` this is the ``new_zero_damping`` constructor.
The ``d3_param`` type is transparent, we can change parameters if we want, for example to evaluate the dispersion energy onces with two-body contributions only and once with non-additive terms added we can override the ``s9`` value and reconstruct the object.


The command line driver
-----------------------

Now that we have the block for the zero damping finished, we can look into creating our main driver and test our program for the first time.
For our command line driver we will go with positional arguments only for a start.
For example to evaluate the PBE0-D3 energies for the forth system of S66 from the previous tutorial we would use:

.. code-block:: text

  param-scanner PBE0 1 dimer.xyz -1 monomer-a.xyz -1 monomer-b.xyz -0.012069608770

.. dropdown:: Geometries for example command

   .. code-block:: text
      :caption: dimer.xyz

      15
      water-peptide dimer, PBE0/def2-QZVP energy: -324.751193159385
      O    -3.2939688    0.4402024    0.1621802 
      H    -3.8134112    1.2387332    0.2637577 
      H    -2.3770466    0.7564365    0.1766203 
      C    -0.6611637   -1.4159110   -0.1449409 
      H    -0.0112009   -2.2770229   -0.2778563 
      H    -1.3421397   -1.3384389   -0.9888061 
      H    -1.2741806   -1.5547070    0.7420675 
      C     0.0935684   -0.1178981   -0.0123474 
      O    -0.4831471    0.9573968    0.1442414 
      N     1.4442015   -0.2154008   -0.0769653 
      H     1.8451531   -1.1259348   -0.2064804 
      C     2.3124436    0.9365697    0.0324778 
      H     1.6759495    1.8048701    0.1672624 
      H     2.9780331    0.8451145    0.8885706 
      H     2.9069093    1.0659902   -0.8697814 

   .. code-block:: text
      :caption: monomer-a.xyz

      3
      water monomer, PBE0/def2-QZVP energy: -76.386381675761
      O    -3.2939688    0.4402024    0.1621802 
      H    -3.8134112    1.2387332    0.2637577 
      H    -2.3770466    0.7564365    0.1766203 

   .. code-block:: text
      :caption: monomer-b.xyz

      12
      peptide monomer, PBE0/def2-QZVP energy: -248.352741874853
      C    -0.6611637   -1.4159110   -0.1449409 
      H    -0.0112009   -2.2770229   -0.2778563 
      H    -1.3421397   -1.3384389   -0.9888061 
      H    -1.2741806   -1.5547070    0.7420675 
      C     0.0935684   -0.1178981   -0.0123474 
      O    -0.4831471    0.9573968    0.1442414 
      N     1.4442015   -0.2154008   -0.0769653 
      H     1.8451531   -1.1259348   -0.2064804 
      C     2.3124436    0.9365697    0.0324778 
      H     1.6759495    1.8048701    0.1672624 
      H     2.9780331    0.8451145    0.8885706 
      H     2.9069093    1.0659902   -0.8697814 

Again we make use of some of the features conveniently provided in ``mctc_env`` and ``mctc_io`` this time we use the possibility to retrieve command line arguments and the reader for geometry inputs.
Also, we will be using the ``error_type`` handler to report errors.

With this we can create our full main program as

.. literalinclude:: first-steps-fortran/app/main.f90
   :language: fortran
   :caption: app/main.f90
   :lines: 1-25,38-41,76-77

This should already be possible to run, while it does not do anything, we can already verify our error handling.
Let's test our program with no arguments:

.. code-block:: text

   ❯ fpm run
   Usage: param-scanner <method> <coeff1> <mol1> ... [dft energy]
   STOP 1

As a next step we define some helper functions, we want to read a geometry from a command line argument defined by its index.
We define a small subroutine which we will include as a contained procedure in our main program:

.. literalinclude:: first-steps-fortran/app/main.f90
   :language: fortran
   :caption: app/main.f90
   :lines: 51-60

Similarly, we define a procedure to read floating point values from a command line argument index.

.. literalinclude:: first-steps-fortran/app/main.f90
   :language: fortran
   :caption: app/main.f90
   :lines: 62-75

With this we can add a loop over all our arguments to populate the reaction:

.. literalinclude:: first-steps-fortran/app/main.f90
   :language: fortran
   :caption: app/main.f90
   :lines: 26-33

And even already call our scanner function

.. literalinclude:: first-steps-fortran/app/main.f90
   :language: fortran
   :caption: app/main.f90
   :lines: 43-47

.. dropdown:: Full main code

   .. literalinclude:: first-steps-fortran/app/main.f90
      :language: fortran
      :caption: app/main.f90
      :lines: 1-33,38-77


Adding output for the results
-----------------------------

Running the example command line should work now, but not produce any output yet.
Next we are going to look into reporting our results.
We add a block to printout our dispersion energies after the zero damping block:

.. literalinclude:: first-steps-fortran/src/scan.f90
   :language: fortran
   :caption: src/scan.f90
   :lines: 143-159

.. note::

   We choose kJ/mol as unit here, but any other unit you prefer can be used here as well.

Also, we need to define the ``label`` parameter for making the printout more pretty:

.. literalinclude:: first-steps-fortran/src/scan.f90
   :language: fortran
   :caption: src/scan.f90
   :lines: 50-51

.. dropdown:: Full library code

   .. literalinclude:: first-steps-fortran/src/scan.f90
      :language: fortran
      :caption: src/scan.f90
      :lines: 1-72,143-161

Running our example command should give the following output:

.. code-block:: text

   ❯ fpm run -- PBE0 1 dimer.xyz -1 monomer-a.xyz -1 monomer-b.xyz
   Energies in kJ/mol
   ------------------------------------------------------------------
    method                       E(2)         E(2+3)           %E(3)
   ------------------------------------------------------------------
    PBE0-D3(0)                 -4.591         -4.607         -0.350% 
   ------------------------------------------------------------------

Currently, we don't process the DFT energy for the reaction.
Let's add this feature to our main driver by including

.. literalinclude:: first-steps-fortran/app/main.f90
   :language: fortran
   :caption: app/main.f90
   :lines: 34-37

Now we can try our main program again with the DFT relative energy:

.. code-block:: text

   ❯ fpm run -- PBE0 1 dimer.xyz -1 monomer-a.xyz -1 monomer-b.xyz -0.012069608770
   Energies in kJ/mol
   ------------------------------------------------------------------
    method                       E(2)         E(2+3)           %E(3)
   ------------------------------------------------------------------
    PBE0-D3(0)                -36.280        -36.296         -0.044%
   ------------------------------------------------------------------

.. dropdown:: Full main code

   .. literalinclude:: first-steps-fortran/app/main.f90
      :language: fortran
      :caption: app/main.f90
      :lines: 1-77


Supporting all damping functions
--------------------------------

Finally, we fill in the remaining blocks for the other damping parameters.
The procedures needed for this are very similar to the ones used for the zero-damping and therefore we leave this as an exercise.
The output with all damping functions should look like

.. code-block:: text

   ❯ fpm run -- PBE0 1 dimer.xyz -1 monomer-a.xyz -1 monomer-b.xyz -0.012069608770
   Energies in kJ/mol
   ------------------------------------------------------------------
    method                       E(2)         E(2+3)           %E(3)
   ------------------------------------------------------------------
    PBE0-D3(0)                -36.280        -36.296         -0.044%
    PBE0-D3(BJ)               -35.463        -35.479         -0.045%
    PBE0-D3M(BJ)              -35.677        -35.693         -0.045%
    PBE0-D3(op)               -35.616        -35.632         -0.045%
   ------------------------------------------------------------------

.. note::

   We have two main categories of damping functions supported with D3.
   First, rational-type damping which makes the dispersion energy go to a constant value at short distances, in this case constant contributions cancel in reactions if the same pairs are present.
   Second, zero-type damping where the dispersion energy is fully removed at short distances.
   Finally, there are versions like the optimized power damping scheme which can be either of rational-type or zero-type damping depending on the power parameters.

We can compare our numbers with the PBE0-D3(BJ) calculations from the last tutorial to confirm our implementation is indeed correct.
Overall, we can see a bit of spread in the results for PBE0 with D3 based on this small example case, however the dispersion energies are rougly all in the same range.
For deciding on which damping function to choose for a density functional, it is always best to check the performance on the specific systems of interest.

.. dropdown:: Full library code

   .. literalinclude:: first-steps-fortran/src/scan.f90
      :language: fortran
      :caption: src/scan.f90
      :lines: 1-161


Summary
-------

This concludes our tutorial on writing a simple command line tool to scan through damping parameters for D3.

.. important::

   In this tutorial we learned
   
   - how to use the ``dftd3`` module and Fortran API to compute dispersion energies
   - create damping parameters for different damping schemes
   - how to handle geometry and error types used in ``dftd3``