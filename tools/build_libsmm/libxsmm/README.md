libxsmm
=======
Library for small matrix-matrix multiplications targeting Intel Architecture (x86). This initial version of the library is targeting the Intel Xeon Phi coprocessor (an instance of the Intel Many Integrated Core Architecture "MIC") particularly using KNC intrinsic functions.

The code can be compiled to native code which is also usable in an offloaded code section (via a FORTRAN directive or a C/C++ pragma). The prerequisite for offloading the code is to compile it to position-independent code (PIC) even when building a static library.

INSTRUCTIONS
============
To compile the library just use:
```
make
```
The static library is produced inside the directory lib. To remove generated files use:
```
make clean
```
or to remove generated files and the built library use:
```
make realclean
```
The interface to the library is:
```
void xsmm_dnn(int M, int N, int K, const double* a, const double* b, double* c)
```
where C(M,N) = C(M,N) + A(M,K) * B(K,N).

The library can be configured to accept row-major (default) or column-major order matrices;
change the variable ROW_MAJOR inside the Makefile file (1 for row-major,
and column-major otherwise); or one can run:
```
make ROW_MAJOR=0
```
The values of the matrix sizes (M,N,K values) can be set by changing the 
variables inside the Makefile file:
```
INDICES_M
INDICES_N
INDICES_K
```
For example:
```
INDICES_M := 2 4
INDICES_N := 1
INDICES_K := 2 4 5
```
it generates the (M,N,K) values:
```
(2,1,2), (2,1,4), (2,1,5)
(4,1,2), (4,1,4), (4,1,5)
```
The fallback for the library is DGEMM if it is called for values other than specified
by INDICES_M, INDICES_N, or INDICES_K.
