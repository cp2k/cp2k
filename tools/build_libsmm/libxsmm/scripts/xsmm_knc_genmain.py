###############################################################################
## Copyright (c) 2013-2014, Intel Corporation                                ##
## All rights reserved.                                                      ##
##                                                                           ##
## Redistribution and use in source and binary forms, with or without        ##
## modification, are permitted provided that the following conditions        ##
## are met:                                                                  ##
## 1. Redistributions of source code must retain the above copyright         ##
##    notice, this list of conditions and the following disclaimer.          ##
## 2. Redistributions in binary form must reproduce the above copyright      ##
##    notice, this list of conditions and the following disclaimer in the    ##
##    documentation and/or other materials provided with the distribution.   ##
## 3. Neither the name of the copyright holder nor the names of its          ##
##    contributors may be used to endorse or promote products derived        ##
##    from this software without specific prior written permission.          ##
##                                                                           ##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       ##
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         ##
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     ##
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      ##
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    ##
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  ##
## TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    ##
## PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    ##
## LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      ##
## NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        ##
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              ##
###############################################################################
## Christopher Dahnken (Intel Corp.), Hans Pabst (Intel Corp.),
## Alfio Lazzaro (CRAY Inc.), and Gilles Fourestey (CSCS)
###############################################################################
import math
import sys

def create_symmetric_interface(dimsM,dimsN,dimsK,RowMajor):
    print "#include <xsmm_knc.h>"
    print "#include <mkl.h>"
    print "#include <stdio.h>"
    print "void xsmm_dnn(int M, int N, int K, const double* a, const double* b, double* c) {"
    print "  int v=((M-1)<<10)+((N-1)<<5)+(K-1);"
    print "  switch(v) {"
    for m in dimsM:
        for n in dimsN:
           for k in dimsK:
                print "    case "+str(((m-1)<<10)+((n-1)<<5)+(k-1))+":"
                print "         dc_smm_dnn_"+str(m)+"_"+str(n)+"_"+str(k)+"(a,b,c);"
                print "         break;"  
    print "    default:"
    if RowMajor==1:
        print "         cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, a, K, b, N, 1.0, c, N);"
    else:
        print "         cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, a, K, b, N, 1.0, c, N);"
    print "  }"
    print "}"

def load_dims(dims):
    dims=map(int, dims) ; dims.sort()
    return list(set(dims))

dimsM = load_dims(sys.argv[4:4+int(sys.argv[1])])
dimsN = load_dims(sys.argv[4+int(sys.argv[1]):4+int(sys.argv[1])+int(sys.argv[2])])
dimsK = load_dims(sys.argv[4+int(sys.argv[1])+int(sys.argv[2]):])

create_symmetric_interface(dimsM,dimsN,dimsK,int(sys.argv[3]))
