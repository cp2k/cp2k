#
# the build script can generate optimized routines packed in a library for
#
# 1) 'nn' => C=C+MATMUL(A,B)
# 2) 'tn' => C=C+MATMUL(TRANSPOSE(A),B)
# 3) 'nt' => C=C+MATMUL(A,TRANSPOSE(B))
# 4) 'tt' => C=C+MATMUL(TRANPOSE(A),TRANPOSE(B))
#
# select a tranpose_flavor from the list 1 2 3 4
#
transpose_flavor=1

# 1) d => double precision real
# 2) s => single precision real
# 3) z => double precision complex
# 4) c => single precision complex 
#
# select a data_type from the list 1 2 3 4
#
data_type=1

#
# target compiler... this are the options used for building the library.
# They should be aggessive enough to e.g. perform vectorization for the specific CPU (e.g. -ftree-vectorize -march=native),
# and allow some flexibility in reordering floating point expressions (-ffast-math).
# Higher level optimisation (in particular loop nest optimization) should not be used.
#
FC ?= gfortran
FCFLAGS ?= -O2 -funroll-loops -ffast-math -ftree-vectorize -march=native -fno-inline-functions
target_compile= $(FC) $(FCFLAGS)

#
# SIMD registers size (in bytes)
# Set to 32 (AVX) or 64 (Xeon Phi) to generate an optimized vector version
# Any other value will not generate the vector version. 
#
SIMD_size=0

#
# target dgemm link options... these are the options needed to link blas (e.g. -lblas)
# blas is used as a fall back option for sizes not included in the library or in those cases where it is faster
# the same blas library should thus also be used when libsmm is linked.
#
OMP_NUM_THREADS=1
#blas_linking=-L/opt/intel/parallel_studio/composerxe-2011.4.191/mkl/lib/intel64 -Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group
blas_linking=-Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group

#
# matrix dimensions for which optimized routines will be generated. 
# since all combinations of M,N,K are being generated the size of the library becomes very large
# if too many sizes are being optimized for. Numbers have to be ascending.
#
dims_small=1 4 5 6 9 13 16 17 22 23

#
# tiny dimensions used are used as primivitves and generated in an 'exhaustive' search.
# They should be a sequence from 1 to N,
# where N is a number that is large enough to have good in cache performance (e.g. for modern SSE cpus 8 to 12)
# Too large (>12?) is not beneficial, but increases the time needed to build the library
# Too small (<8)   will lead to a slow library, but the build might proceed quickly
# The minimum number for a successful build is 4
#
dims_tiny=1 2 3 4 5 6 7 8 9 10 11 12

#
# host compiler... this is used only to compile a few tools needed to build the library. The library itself is not compiled this way.
# This compiler needs to be able to deal with some Fortran2003 constructs.
#
host_compile=$(FC)

#
# number of processes to use in parallel for compiling / building and benchmarking the library.
# Should *not* be more than the physical (available) number of cores of the machine
#
tasks?=16

