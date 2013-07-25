#!/usr/bin/python
# -*- coding: utf-8 -*-

from subprocess import check_call
import itertools

#===============================================================================
def writefile(fn, content):
    f = open(fn, "w")
    f.write(content)
    f.close()
    if(fn.endswith(".h") or fn.endswith(".cu")):
        check_call(["indent",fn])
    print("Wrote: "+fn)

#===============================================================================
def main():
    all_kernels = []
    #sizes = (5,8,13,23,26,)
    sizes = (5,8,13,16,23,26,)
    for (m,n,k) in itertools.product(sizes, sizes, sizes):
        kernel_name = "stack_mm_mnk_kepler_NxNd_%d_%d_%d_2"%(m,n,k)
        all_kernels.append(kernel_name)
        fn = "./dbcsr_kernels/"+kernel_name

        signature  = "int launch_"+kernel_name+"(int *param_stack, int stack_size, "
        signature += "cudaStream_t stream, int m_max, int n_max, int k_max, "
        signature += "double *a_data, double *b_data, double *c_data)"
        writefile(fn+".h", content=signature+";\n")

        output = "// DBCSR_KERNEL datatype=dbcsr_type_real_8, homogeneous_only=True, m=%d, n=%d, k=%d\n"%(m,n,k)
        output += '#include "dbcsr_kernel.h"\n'
        output += '#include "dbcsr_generic_kernel.h"\n\n'
        output += signature+"{\n"
        output += "int shared_size = 0;\n"
        output += "int careful = (stack_size / GROUPING);\n"
        output += "int nruns = stack_size - careful * GROUPING;\n"
        output += "stack_mm_mnk_kepler_NxNd<%d,%d,%d,2> "%(m,n,k)
        output += "<<< ((stack_size + GROUPING - 1) / GROUPING), 192, shared_size, stream >>>\n"
        output += "(param_stack, careful, nruns, m_max, n_max, k_max, \n"
        output += "a_data, b_data, c_data);\n"
        output += "return(0);\n"
        output += "}\n"
        writefile(fn+".cu", output)

    output = "\n".join(["LIBNV_OBJECTS += %s.o"%k for k in all_kernels])
    writefile("DBCSR_KERNELS_OBJDEFS", output)

    output = "\n".join(["%s.o : dbcsr_kernel.h dbcsr_generic_kernel.h %s.cu  %s.h"%(k,k,k) for k in all_kernels])
    writefile("DBCSR_KERNELS_DEPENDENCIES", output)
#===============================================================================

main()

#EOF
