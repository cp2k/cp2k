#! /usr/bin/env python

###############################################################################
# dbcsr_testing.py 
# Example:
#     ./dbcsr_testing.py 
#     will run the tests with the default number of nodes, threads and mpi call
#
#     ./dbcsr_testing.py --nnodes="1 2 1" --nthreads="2 2 1" --mpirun="mpirun -np "
#     will run the tests with the 1 and 2 nodes, 2 threads and with a call 
#     mpi through mpirun -np
# 
#\author  V. Weber
#\date    2010
#\version 1.1
#\changelog
# - 2011-11-18 [UB] Updates code, removes dependency on Numeric
#
###############################################################################

import os, sys, getopt, time



# Define a function to open the executable
def local_popen( cmdline ):
   pipe = os.popen( cmdline )
   return pipe


# default parameters
nnodes_min = 1
nnodes_max = 1
nnodes_stride = 1
nthreads_min = 1
nthreads_max = 1
nthreads_stride = 1
mpirun = "mpiexec -np "
openmp = "export OMP_NUM_THREADS="
exe = "../dbcsr_test_driver"
use_mpi = 0
use_omp = 0

# parsing
try:
    (optlist, args) = getopt.getopt(sys.argv[1:],"n:t:m:e:",["nnodes=","nthreads=","mpirun=","exe="])
    for o, a in optlist:
        if o in ("-n", "--nnodes"):
           nodes=tuple((int(d) for d in a.split()))
           if len(nodes) != 3:
               raise Exception("The nnodes (min,max,stride) argument must have an associated string of 3 integers (you have %d)." % len(nodes))
           nnodes_min = nodes[0]
           nnodes_max = nodes[1]
           nnodes_stride = nodes[2]
           use_mpi = 1
        elif o in ("-t", "--nthreads"):
           threads=tuple((int(d) for d in a.split()))
           if len(threads) != 3:
                raise Exception("The nthreads (min,max,stride) argument must have an associated string of 3 integers (you have %d)." % len(threads))
           nthreads_min = threads[0]
           nthreads_max = threads[1]
           nthreads_stride = threads[2]
           use_omp = 1
        elif o in ("-m", "--mpirun"):
           mpirun = a
           use_mpi = 1
        elif o in ("-e", "--exe"):
           exe = a
        else:
            raise Exception("Invalid arguments")

except Exception, e:
    print "Error parsing command line arguments: " + str(e)
    sys.exit(1)


# If filename cannot be opened, send output to sys.stderr
filename = "testing_results.txt"
try:
     f = open(filename, 'w')
except IOError:
     f = sys.stdout

executable = " -e \"" + str(exe) + "\""
nodes_triplett = " -n \"" + str(nnodes_min) + " " + str(nnodes_max) + " " + str(nnodes_stride) + "\""
threads_triplett = " -t \"" + str(nthreads_min) + " " + str(nthreads_max) + " " + str(nthreads_stride) + "\""

# Let's go !
print " "
print "  ----------------------- Testing DBCSR Routines ----------------------"
print " "
if use_mpi == 1 and use_omp == 1:
   print "  -- Effective testing command: python dbcsr_testing.py" + nodes_triplett + threads_triplett + executable
elif use_mpi == 1 and use_omp == 0:
   print "  -- Effective testing command: python dbcsr_testing.py" + nodes_triplett + executable
elif use_mpi == 0 and use_omp == 1:
   print "  -- Effective testing command: python dbcsr_testing.py" + threads_triplett + executable
else:
   print "  -- Effective testing command: python dbcsr_testing.py" + executable
print " "
print "  -- Min Number of nodes =", nnodes_min
print "  -- Max Number of nodes =", nnodes_max
print "  -- Nodes stride =", nnodes_stride
print "  -- Min Number of threads =", nthreads_min
print "  -- Max Number of threads =", nthreads_max
print "  -- Thread stride =", nthreads_stride
print "  -- Detailed results are stored in", filename


ntests = {}
nfails = {}
for nnodes in range(nnodes_min, nnodes_max+1, nnodes_stride):
   for nthreads in range(nthreads_min, nthreads_max+1, nthreads_stride):
      ntests[(nthreads,nnodes)] = 0
      nfails[(nthreads,nnodes)] = 0
 

t_start = time.time()
for nnodes in range(nnodes_min, nnodes_max+1, nnodes_stride):
   for nthreads in range(nthreads_min, nthreads_max+1, nthreads_stride):
      #
      # run over parameter files
      for root, dirs, all_files in os.walk('./'):
         par_files = filter(lambda x: x.endswith('.par'), all_files)
         for file in par_files: 
            #
            #
            print " "
            print "  ----------------------- nnodes", nnodes, " nthreads", nthreads,"----------------------"
            print "  -- Parameter file: ", file, " in directory: ",root
            sys.stdout.flush()
            #
            #
            if use_mpi == 1:
               test1 = local_popen( openmp + str(nthreads) + ";" + mpirun + str(nnodes) + " " + exe + " < " + root + "/" + file )
            else:
               test1 = local_popen( openmp + str(nthreads) + ";" + " " + exe + " < " + root + "/" + file )
            for line in test1.readlines():
               f.write(str(line))
               if "TESTING" in line : 
                  print line,
                  ntests[(nthreads,nnodes)] += 1
                  if "FAILED" in line : nfails[(nthreads, nnodes)] += 1
               #
            #
            f.flush()
            sys.stdout.flush()
            #
            #
         #
      #
   #
#
#

t_end = time.time()

def dict_val_sum(d):
   return sum(d.values())

print " "
print "  ------------------------------- Report ------------------------------"
print " "
print "  -- Number of tests =", dict_val_sum(ntests)
print "  -- Number of failures =", dict_val_sum(nfails)
print "  -- Completed in ", t_end-t_start, " seconds"
print " "
print " "
print "  -- Decomposition per node/thread"
#
for nnodes in range(nnodes_min, nnodes_max+1, nnodes_stride):
   for nthreads in range(nthreads_min, nthreads_max+1, nthreads_stride):
      print " "
      print "    -- Number of nodes =", nnodes
      print "    -- Number of threads =", nthreads
      print "    -- Number of tests =", ntests[(nthreads, nnodes)]
      print "    -- Number of failures =", nfails[(nthreads, nnodes)]
      print " "
   #
#
print "  ---------------------------------------------------------------------"
f.flush()
sys.stdout.flush()

# This may close the sys.stdout stream, so make it the last statement
f.close()
