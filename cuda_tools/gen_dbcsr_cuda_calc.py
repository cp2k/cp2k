#!/usr/bin/python
# -*- coding: utf-8 -*-

from subprocess import check_call
import re

#===============================================================================

class Datatype(object):
    def __init__(self, name, ctype):
        self.name = name
        self.ctype = ctype

    def __str__(self):
        return(self.name)


dbcsr_type_real_4 = Datatype("dbcsr_type_real_4", "float")
dbcsr_type_real_8 = Datatype("dbcsr_type_real_8", "double")
dbcsr_type_complex_4 = Datatype("dbcsr_type_complex_4", "float")
dbcsr_type_complex_8 = Datatype("dbcsr_type_complex_8", "double")
ALL_DATATYPES = (dbcsr_type_real_4, dbcsr_type_real_8, dbcsr_type_complex_4, dbcsr_type_complex_8)

#===============================================================================
class Kernel(object):
    def __init__(self, name, datatype, homogeneous_only, m=None, n=None, k=None):
        self.name=name
        self.datatype=datatype
        self.m=m
        self.n=n
        self.k=k
        self.homogeneous_only=homogeneous_only

    def __repr__(self):
        return("<Kernel "+self.name+">")

    def condition(self):
        if(self.m and self.n and self.k):
            return("(m_max==%s && n_max==%s && k_max==%s)"%(self.m, self.n, self.k))
        return(None)

    def launch(self):
        output  = "stat = launch_"+self.name
        output += "(param_stack, stack_size, *stream, m_max, n_max, k_max,\n"
        output += "("+self.datatype.ctype+" *) a_data,\n"
        output += "("+self.datatype.ctype+" *) b_data,\n"
        output += "("+self.datatype.ctype+" *) c_data);\n"
        return(output)


#===============================================================================
def gen_body(all_kernels):
    output = "switch (datatype){ \n"
    for dt in ALL_DATATYPES:
        output += "case %s:\n"%dt

        hom_kernels = [k for k in all_kernels if k.homogeneous_only and k.datatype==dt]
        inhom_kernels = [k for k in all_kernels if not k.homogeneous_only and k.datatype==dt]

        if(len(hom_kernels) > 0):
            output += "if(def_mnk) {\n"

            spez_kernels = [k for k in hom_kernels if k.condition()!=None]
            univ_kernels = [k for k in hom_kernels if k.condition()==None]

            for k in spez_kernels:
                output += "if "+k.condition() +"{\n"
                output += k.launch() +"\n break; \n"
                output += "}\n"

            assert(len(univ_kernels)<=1)
            if(len(univ_kernels)==1):
                output += univ_kernels[0].launch() +"\n break; \n"
            output += "}\n"


        assert(len(inhom_kernels)<=1) #there musst be at most one inhom kernel
        if(len(inhom_kernels)==1):
            assert(inhom_kernels[0].condition()==None)
            output += inhom_kernels[0].launch() + "break;\n"
        else:
            output += "return(6); //not universal kernel for inhm. stacks avail.\n"


    output += "} \n"
    return(output)


#===============================================================================
def gen_header(all_kernels):
    output  = "/******************************************************************************\n"
    output += "*  CP2K: A general program to perform molecular dynamics simulations\n"
    output += "*  Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group\n"
    output += "*****************************************************************************/\n"
    output += "#include <cuda_runtime.h>\n"
    output += "#include <stdio.h>\n"
    output += "#include <sm_11_atomic_functions.h>\n"
    output += '#include "dbcsr_cuda.h"\n'
    for k in all_kernels:
        output += '#include "dbcsr_kernels/'+k.name+'.h"\n'


    output += r"""

#define dbcsr_type_real_4     1
#define dbcsr_type_real_8     3
#define dbcsr_type_complex_4  5
#define dbcsr_type_complex_8  7

//    static const int verbose_print = 0;


/**
 * \brief Bridge routine to call appropriate CUDA kernel.
 */
extern "C" int
dc_do_stack_cu(int *param_stack, int stack_size, int nparams, int datatype,
        void *a_data, void *b_data, void *c_data,
        int m_max, int n_max, int k_max, int def_mnk, cudaStream_t* stream)
{
    int stat;

/*  if (verbose_print)
    {
      printf ("A data %p.\n", a_data);
      printf ("B data %p.\n", b_data);
      printf ("C data %p.\n", c_data);
      printf ("params %p.\n", param_stack);
    }
*/
   // printf("Got m,n,k: %d %d %d; %d.\n", m_max, n_max, k_max, stack_size);
    """
    return(output)

#===============================================================================
def gen_footer():
    output  = "if(stat!=0) return(stat);\n"
    output += "if (cuda_error_check (cudaGetLastError ())) return 1;\n"
    output += "return 0;\n"
    output += "};\n"
    output += "//EOF\n"
    return(output)

#===============================================================================
def main():
    all_kernel_names  = ["stack_mm_r", "stack_mm_d", "stack_mm_c", "stack_mm_z","stack_mm_mnk_d"]
    #all_kernel_names += ["stack_mm_mnk_sq8_d"]
    #unused: stack_mm_mnk_sq5_d, stack_mm_mnk_sq13_d, stack_mm_mnk_sq13_d_slow

    for line in open("DBCSR_KERNELS_OBJDEFS").readlines():
        kernel_name = line.split("+=")[-1].strip()[:-2]
        all_kernel_names.append(kernel_name)

    all_kernels = []
    for kn in all_kernel_names:
        fn = "./dbcsr_kernels/"+kn+".cu"
        print("Reading: "+fn)
        content = open(fn).read()
        descr_line = re.search("DBCSR_KERNEL(.*)\n", content).group(1)
        #print "found: "+descr_line
        descr = eval("dict("+descr_line+")")
        kernel = Kernel(name=kn, **descr)
        all_kernels.append(kernel)

    print("Found %d dbcsr cuda kernels."%len(all_kernels))

    output  = gen_header(all_kernels)
    output += gen_body(all_kernels)
    output += gen_footer()

    fn = "dbcsr_cuda_calc.cu"
    f = open(fn, "w")
    f.write(output)
    f.close()

    check_call(["indent",fn])
    print("Wrote: "+fn)

#===============================================================================

main()

#EOF
