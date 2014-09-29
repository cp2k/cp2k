#!/bin/env python

import math
import sys

def createigemm(M,N,K):
    iparts=int(math.floor(M/8))
    fparts=M%8
    if fparts==0:
        mparts=iparts
    else:
        mparts=iparts+1
    print "#include <immintrin.h>"
    print "#include <micsmmmisc.h>"
    print " "

    print "void dc_smm_dnn_"+str(M)+"_"+str(N)+"_"+str(K)+"(double* a,double* b,double* c) {"
    print "#ifdef __MIC__"
    print "  int i;"
    for k in range(0,K):
        print "  __m512d xa"+str(k)+";"
    for k in range(0,K):
        print "  __m512d xb"+str(k)+";"
    print "  __m512d xc0;"

    for m in range(0,8*mparts,8):
        mm=min(m+7,M-1)
        maskval=(1<<(mm-m+1))-1
        for k in range(0,K):
            print "  xa"+str(k)+" = _MM512_MASK_LOADU_PD(&a["+str(M*k)+"+"+str(m)+"]," +str(maskval)+");"
        print "  for(i=0;i<"+str(N)+";++i) {"
        print "    xc0 = _MM512_MASK_LOADU_PD(&c[i*"+str(M)+"+"+str(m)+"]," +str(maskval)+");"
        for k in range(0,K):
            print "    xb"+str(k)+"=_mm512_set1_pd(b[i*"+str(K)+"+"+str(k)+"]);"
        for k in range(0,K):
            print "    xc0=_mm512_mask3_fmadd_pd(xa"+str(k)+",xb"+str(k)+",xc0," +str(maskval)+");"
        print "    _MM512_MASK_STOREU_PD(&c[i*"+str(M)+"+"+str(m)+"],xc0," +str(maskval)+");"
        print "  }"
    print "#else"
    print "  int m, n, k;"
    print "  for (m=0; m<"+str(M)+"; m++) {"
    print "    for (n=0; n<"+str(N)+"; n++) {"
    print "      for (k=0; k<"+str(K)+"; k++) {"
    print "        c[n*"+str(M)+"+m]+=a[k*"+str(M)+"+m]*b[n*"+str(K)+"+k];"
    print "      }"
    print "    }"
    print "  }"
    print "#endif"
    print "}"
    print " "


createigemm(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]))
