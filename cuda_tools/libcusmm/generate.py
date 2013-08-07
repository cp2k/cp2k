#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
from os import path
from glob import glob
from itertools import product, chain
from optparse import OptionParser

from lib.kernel_default import DefaultKernel
from lib.kernel_medium import MediumKernel
from lib.kernel_mediumDB import MediumDBKernel
from lib.kernel_small import SmallKernel
from lib.kernel_tiny import TinyKernel


#===============================================================================
def main():
    default_sizes = (5,8,12,13,23,26)

    usage = "Generator of LibCuSMM. The Library for Cuda Small Matrix Multiplications."
    parser = OptionParser(usage)
    parser.add_option("-p", "--params", metavar="filename.txt",
        default="params_default.txt",
        help="Default: %default")
    parser.add_option("-s", "--sizes", metavar="<comma-separated-list>",
        default=",".join([str(i) for i in default_sizes]),
        help="Default: %default")

    (options, args) = parser.parse_args()
    assert(len(args)==0)
    sizes = eval("("+options.sizes+")")
    param_fn = options.params

    plan = make_plan(sizes, param_fn)
    prep_build_dir()
    gen_dispatch_code(plan)
    gen_makefile(plan)

    d = len([x for x in plan.values() if isinstance(x, DefaultKernel)])
    print("Libcusmm: Generated %d kernels (%d defaults)."%(len(plan), d))

#===============================================================================
def make_plan(sizes, param_fn):
    all_kernels = eval(open(param_fn).read())


    plan = {}
    for (m, n, k) in product(sizes, sizes, sizes):
        possible_kernels = [kern for kern in all_kernels if kern.can_handle(m,n,k)]
        if(len(possible_kernels) == 1):
            plan[(m,n,k)] = possible_kernels[0]
        elif(len(possible_kernels) > 1):
            raise(Exception("found more than one kernel for %dx%dx%d"%(m,n,k)))
        else:
            plan[(m,n,k)] = DefaultKernel(m=m, n=n, k=k)

    return(plan)


#===============================================================================
def prep_build_dir():
    if(not path.exists("build")):
        os.mkdir("build")

    for fn in glob("./src/*"):
       bn = path.basename(fn)
       if(not path.exists("./build/"+bn)):
           os.symlink("../src/"+bn, "./build/"+bn)


#===============================================================================
def gen_dispatch_code(plan):
    output  = "/******************************************************************************\n"
    output += "*  CP2K: A general program to perform molecular dynamics simulations\n"
    output += "*  Copyright (C) 2000 - 2013 the CP2K developers group\n"
    output += "*****************************************************************************/\n"

    for i in get_includes(plan):
        output += '#include "%s"\n'%i
    output += "\n\n"

    for kern in plan.values():
        output += kern.launcher_code() + "\n\n"

    output += "int libcusmm_process_d(int *param_stack, int stack_size,"
    output += "cudaStream_t stream, int m, int n, int k, "
    output += "double * a_data, double * b_data, double * c_data){\n"


    #generate jump table -------------------------------------------------------
    m_vals = list(set([m for (m,n,k) in plan.keys()]))
    n_vals = list(set([n for (m,n,k) in plan.keys()]))
    k_vals = list(set([k for (m,n,k) in plan.keys()]))
    assert(len(m_vals) * len(n_vals) * len(k_vals) < pow(2,16))

    output += "int idx = 0;\n"
    output += "bool missing = false;\n\n"
    output += "switch(m){\n"
    for i, m in enumerate(m_vals):
        output += "case %d: idx = %d; break;\n"%(m, i)
    output += "default: missing = true;\n"
    output += "}\n\n"

    output += "idx *= %d;\n"%len(n_vals)
    output += "switch(n){\n"
    for i, n in enumerate(n_vals):
        output += "case %d: idx += %d; break;\n"%(n, i)
    output += "default: missing = true;\n"
    output += "}\n\n"

    output += "idx *= %d;\n"%len(k_vals)
    output += "switch(k){\n"
    for i, k in enumerate(k_vals):
        output += "case %d: idx += %d; break;\n"%(k, i)
    output += "default: missing = true;\n"
    output += "}\n\n"

    # Fall back to non-templated universal kernel for homogenous stacks.
    output += "if(missing) // fallback\n"
    output += "return launch_cusmm_kernel_fallback"
    output += "(param_stack,stack_size,stream, m,n,k,a_data,b_data,c_data);\n\n"

    idx_map = dict()
    for (m,n,k) in plan.keys():
        idx = (m_vals.index(m)*len(n_vals) + n_vals.index(n))*len(k_vals) + k_vals.index(k)
        idx_map[idx] = (m,n,k)

    output += "switch(idx){\n"
    for idx in sorted(idx_map.keys()):
        mnk = idx_map[idx]
        output += "case %d:\n"%idx
        output += "// m=%d, n=%d, k=%d\n"%mnk
        output += "return launch_"+plan[mnk].name
        output += "(param_stack, stack_size, stream, %d, %d, %d, "%mnk
        output += "a_data, b_data, c_data);\n\n"
    output += "}\n\n"

    output += "return -1; // should never happen\n"
    output += "}\n"
    output += "//EOF\n"

    writefile("./build/libcusmm.cu", output)


#===============================================================================
def gen_makefile(plan):
    output  = "libcusmm.a : libcusmm.o cusmm_kernel_fallback.o\n"
    output += "\t$(AR) libcusmm.a libcusmm.o cusmm_kernel_fallback.o\n"

    output += "%.o: %.cu\n"
    output += "\t$(NVCC) -c $(NVFLAGS) $<\n"

    output += "libcusmm.o : libcusmm.cu " +(" ".join(get_includes(plan))) + "\n"
    output += "cusmm_kernel_fallback.o : cusmm_kernel_fallback.cu\n"

    writefile("./build/Makefile", output)

#===============================================================================
def get_includes(plan):
    includes  = list(set([kern.include() for kern in plan.values() ]))
    includes += ["cusmm_kernel_fallback.h",]
    return(includes)


#===============================================================================
def writefile(fn, content):
    if(path.exists(fn)):
        f = open(fn, "r")
        old_content = f.read()
        f.close()
        if(old_content == content):
            return

    f = open(fn, "w")
    f.write(content)
    f.close()

    #if(fn.endswith(".h") or fn.endswith(".cu")):
    #    check_call(["indent",fn])
    #print("Wrote: "+fn)


#===============================================================================

main()

#EOF

