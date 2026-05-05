Parameter optimization
======================

Optimization of parameters is available with the `tblite-fit`_ subcommand.

.. _tblite-fit: https://github.com/awvwgk/tblite/blob/main/man/tblite-fit.1.adoc
.. _tblite-param: https://github.com/awvwgk/tblite/blob/main/man/tblite-param.1.adoc


Setting up a fit
----------------

In this example we try to create a parametrization for titanium oxides based on GFN2-xTB

.. important::

   We will assume a data set to fit against is already existing and prepared here.

First, we create our base parameter file using the `tblite-param`_ command.

.. code-block:: none

   tblite param --method gfn2 --output gfn2-xtb.toml

Now, we can inspect the parameter file, we will find six main sections named *hamiltonian*, *dispersion*, *repulsion*, *charge*, *thirdorder*, and *multipole* as well as the *element* sections.

In the next step we create our input file for the parameter optimization, the input file is written in TOML:

.. code-block:: toml
   :caption: tbfit.toml

   script = "exec ./run.sh"

   [mask]
   hamiltonian = {}
   dispersion = {}
   repulsion = {}
   charge = {}
   thirdorder = {}
   multipole = {}
   [mask.element]
   O = {lgam=[false, true]}
   Ti = {lgam=[true, false, true]}

The most important part is the *mask* section as it selects the parameters to optimize.
Here we include all sections from above to enable them in the element records automatically.
To actually select the elements we add them in *mask.element*, at this place we can also overwrite the automatic defaults.
We will only disable the shell hardnesses for the s shells here, since those should use the atomic hardnesses unmodified.

With our exported parameter file and the fit input we start a test run to check if everything is working correctly

.. code-block:: none

   tblite fit --dry-run gfn2-xtb.toml tbfit.toml

If *tblite* returns an error, check the input file again for typos.

The next step is to create the runner, usually a shell script is sufficient for this purpose.
The script will get its input via environment variables.
The fit driver will set the name of the *tblite* executable in ``TBLITE_EXE``, we can use this variable to use the same version of *tblite* to read the modified parameter file again and perform single points.
Also, we get the name of the created parameter file in ``TBLITE_PAR`` and the name of the data output we are supposed to generate in the command in ``TBLITE_OUT``.

We use the script shown below

.. code-block:: sh
   :caption: run.sh

   #!/usr/bin/env bash

   # Uncomment for debugging
   #set -ex

   # Input from fit runner with defaults for standalone use
   fitpar=$(realpath ${TBLITE_PAR:-./fitpar.toml})
   output=$(realpath ${TBLITE_OUT:-.data})
   tblite=$(which ${TBLITE_EXE:-tblite})

   # Balance to get nproc == njob * nthreads
   nthreads=1
   njob=$(nproc | awk "{print int(\$1/$nthreads)}")

   # Temporary data file
   data=.data

   # Temporary wrapper
   wrapper=./wrapped_runner

   # Arguments for tblite runner
   tblite_args="run --param \"$fitpar\" --grad results.tag coord"

   # Ad hoc error in case the Hamiltonian does not work
   # (SCC does not converge or similar)
   penalty="1.0e3"

   # Create our wrapper script
   cat > "$wrapper" <<-EOF
   #!/usr/bin/env bash
   if [ -d \$1 ]; then
     pushd "\$1" > /dev/null 2>&1
     test -f "$data" && rm "$data"
     OMP_NUM_THREADS=1 "$tblite" $tblite_args  > tblite.out 2> tblite.err \
       || echo "0.0 $penalty  # run: \$1" > "$data"
     "$tblite" tagdiff --fit results.tag reference.tag >> "$data" \
       || echo "0.0 $penalty  # diff: \$1" >> "$data"
   fi
   EOF
   chmod +x "$wrapper"

   # Create the actual multiprocessing queue
   printf "%s\0" data/*/ | xargs -n 1 -P $njob -0 "$wrapper"

   # Collect the data
   cat data/*/$data > "$output"

   # Cleanup
   rm data/*/$data
   rm "$wrapper"

In short this script will process multiple structures in parallel and create the final data output.
The structures we want to evaluate with this script are expected to be in the subdirectories of the ``data`` directory stored as Turbomole ``coord`` files.
Always run the command you want to use before start the fit driver in production mode, for this purpose also uncomment the debugging line in the script and run it with

.. code-block:: none

   TBLITE_PAR=./gfn2-xtb.toml TBLITE_OUT=data.txt time ./run.sh

The most important step is to check the output generated in ``data.txt``.
Failed runs will be marked with the directory name, inspect the output of *tblite* to check whether there is a setup error that has to be fixed first before starting into production.
We also use the ``time`` command here to determine the approximate runtime of our script.
A good target is a maximum runtime of one second.

.. note::

   While long evaluation times are possible, it makes the fit more difficult in practice.
   A good rule of thumb is that everything above 10 seconds will become problematic.

   If necessary the data set has to be cut down by a bit or smaller input have to be used.
   Alternatively, a more powerful machine should be used if available.
   Especially for large numbers of data more cores can help to reduce the wall time.

Finally, if you are happy with the setup start the actual fit in verbose mode.

.. code-block:: none

   tblite fit -v gfn2-xtb.toml tbfit.toml

Using the verbose printout will show the objective function in every step as well as the relative change in all parameters.
Initially, the verbose printout is useful to track the stability of the fit.

.. tip::

   To create long-running fits detach the fit driver from the shell using

   .. code-block:: none

      nohup tblite fit -v gfn2-xtb.toml tbfit.toml > fitlog.txt 2>&1 &

   You can check the output in ``fitlog.txt`` while the fit driver is running in the background.

Once the fit finishes the optimized parameters are found in ``fitpar.toml``.
You should update the *meta* section to specify how those parameters were fitted.

.. code-block:: toml
   :caption: fitpar.toml

   [meta]
   version = 1
   name = "...-xTB"
   reference = "..."

.. important::

   The method we are fitting is an extended tight-binding (xTB) Hamiltonian.
   A parametrization like the second *geometry*, *frequency*, *non-covalent interaction* parametrization (GFN2) results in the final method name GFN2-xTB.
   Parametrizations are always given as prefix to xTB, it is highly discouraged to add something to the xTB part itself, *i.e.* there is no xTB2.
