Running ``tblite`` in parallel
==============================

The ``tblite`` program uses shared memory OpenMP parallelization.
To calculate larger systems, an appropriate OMP stacksize must be provided. Choose a reasonably large number with

.. code:: sh
  
   export OMP_STACKSIZE=4G
  
.. note::

   Note that the memory requirement will increase with the system size *and* the number
   of requested threads.

To distribute the number of threads reasonably within the OpenMP section,
it is recommended to use

.. code:: sh
  
   export OMP_NUM_THREADS=<ncores>,1

You might want to deactivate nested OMP constructs by

.. code:: sh

   export OMP_MAX_ACTIVE_LEVELS=1

.. tip::

   Most OpenMP regions allow customization of the scheduling by setting the ``OMP_SCHEDULE`` environment variable.
   For many threads, the ``dynamic`` schedule has proven to provide good load balance across all threads.

Depending on the linear algebra backend used when compiling ``tblite``, different OpenMP-threaded versions are available.
Usually, those backends repect the settings entered by ``OMP_NUM_THREADS``.
However, you can still adjust the parallelization indivdually for the linear algebra backend.
For Intel's Math Kernel Library, the environment variable is ``MKL_NUM_THREADS``.
For the OpenBLAS backend, use ``OPENBLAS_NUM_THREADS`` instead.
It is then exported for the current session as follows:

.. code:: sh
  
   export MKL_NUM_THREADS=<ncores>

or respectively:

.. code:: sh
  
   export OPENBLAS_NUM_THREADS=<ncores>

When computing large systems, the limit of memory allocated for variables saved on the stack should be adjusted, as exceeding this limit can lead to segmentation faults.
This adjustment can be made on UNIX systems (Linux and macOS) using the ``ulimit`` command, as follows:

.. code:: sh

   ulimit -s unlimited

Parallelisation using the python API
-------------------------------------

When using ``tblite``'s python API, the parallelization behavior is also controlled via the aforementioned environment variables.
These variables can be set in the terminal before launching the python code containing the ``tblite`` calculations.
Another possibility is to set the varaibles from within the python code.
This can be achieved by the ``os.environ`` object, for details consider their `documentation <https://docs.python.org/3/library/os.html#os.environ>`__ for details.

To set up OpenMP in a manner analogous to the above:

.. code:: python

   import os
   import psutil
   os.environ['OMP_STACKSIZE'] = '3G'
   os.environ['OMP_NUM_THREADS'] = f'{len(psutil.Process().cpu_affinity())},1'
   os.environ['OMP_MAX_ACTIVE_LEVELS'] = '1'

The maximum stack size can also set from within python.
We tested this using the `resource <https://docs.python.org/3/library/resource.html#resource-limits>`__ module.

To set the stack size to unlimited, the following code snippet can be used:

.. code:: python

   import resource
   resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
