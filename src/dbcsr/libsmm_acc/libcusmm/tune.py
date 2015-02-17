#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from os import path
from os.path import basename
from glob import glob
from itertools import product, chain
from optparse import OptionParser

from kernels.cusmm_dnt_largeDB  import Kernel_dnt_largeDB
from kernels.cusmm_dnt_largeDB2 import Kernel_dnt_largeDB2
from kernels.cusmm_dnt_medium   import Kernel_dnt_medium
from kernels.cusmm_dnt_small    import Kernel_dnt_small
from kernels.cusmm_dnt_tiny     import Kernel_dnt_tiny

ALL_KERNELS = (Kernel_dnt_tiny, Kernel_dnt_small, Kernel_dnt_medium, Kernel_dnt_largeDB, Kernel_dnt_largeDB2,)

#===============================================================================
def main():
    if(len(sys.argv) < 2):
        print("Usage: tune.py <blocksize 1> ... <blocksize N>")
        sys.exit(1)

    blocksizes = [int(i) for i in sys.argv[1:]]
    assert(len(set(blocksizes)) == len(blocksizes))
    blocksizes.sort()

    triples  = combinations(*blocksizes)

    for (m, n, k)  in triples:
        outdir = "tune_%dx%dx%d/"%(m,n,k)
        if(path.exists(outdir)):
            print("Directory %s exists already, skipping."%outdir)
            continue
        os.mkdir(outdir)
        gen_benchmark(outdir, m, n, k)
        gen_jobfile(outdir, m, n, k)
        gen_makefile(outdir)

    #gen_collect(outdir, triples)


#===============================================================================
def format_params(params):
    output = []
    order = ['m','n','k','tile_m', 'tile_n', 'w', 'v', 'split_thread', 'threads', 'blockdim', 'grouping']
    for k in order:
        if(params.has_key(k)):
            output.append("%s=%d"%(k, params[k]))

    for k in params.keys():
        if(k not in order):
            output.append("%s=%d"%(k, params[k]))

    return("(" + ", ".join(output) +")")



#===============================================================================
def gen_benchmark(outdir, m, n, k):
    includes = []
    launcher_codes = []
    launchers = []
    kernel_descr = []

    for kernclass in ALL_KERNELS:
        params = kernclass.promising_parameters(m, n, k)
        if(params == 0):
            continue

        for p in params:
            kern = kernclass(**p)
            includes.append("../kernels/"+kern.include())
            launcher_codes.append(kern.launcher_code())
            launchers.append("launch_"+kern.name)
            kernel_descr.append(kernclass.__name__ + format_params(p))

    print("Found %d parameter sets for %dx%dx%d"%(len(launchers), m, n, k))
    if(len(launchers)==0): return
    #assert(len(launchers)> 0)

    incl_output = '#include "../kernels/cusmm_common.h"\n'
    for i in set(includes):
        incl_output += '#include "%s"\n'%i
    incl_output += "\n\n"

    MAX_LAUNCHERS_PER_EXE = 10000
    LAUNCHERS_PER_OBJ = 100

    n_exe_files = len(launcher_codes)/MAX_LAUNCHERS_PER_EXE + 1
    launchers_per_exe = len(launcher_codes) / n_exe_files + 1

    for i in range(n_exe_files):
        A =  i * launchers_per_exe
        B = min((i+1)*launchers_per_exe, len(launcher_codes))
        for j in range((B-A)/LAUNCHERS_PER_OBJ + 1):
            output = incl_output
            a = A + j*LAUNCHERS_PER_OBJ
            b = min(A + (j+1)*LAUNCHERS_PER_OBJ, B)
            output += "\n\n".join(launcher_codes[a:b])
            fn = outdir+"/tune_%dx%dx%d_exe%d_part%d.cu"%(m, n, k, i, j)
            writefile(fn, output)

        output = '#include "../libcusmm_benchmark.h"\n\n'
        for l in launchers:
            output += "int "+ l + "(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data);\n"

        output += "\n"
        output += "int main(int argc, char** argv){\n"
        output += "KernelLauncher launchers[%d];\n"%(B-A)
        output += "char *kernel_descr[%d];\n"%(B-A)

        for j in range(B-A):
            output += "launchers[%d]    = %s;\n"%(j, launchers[A+j])
            output += 'kernel_descr[%d] = "%s";\n'%(j, kernel_descr[A+j])
        output += "return libcusmm_benchmark(%d, %d, %d, %d, launchers, kernel_descr, true);\n"%(m, n, k, B-A)
        output += "}\n"

        fn = outdir+"/tune_%dx%dx%d_exe%d_main.cu"%(m, n, k, i)
        writefile(fn, output)

#===============================================================================
def gen_jobfile(outdir, m, n, k):
    t = "/tune_%dx%dx%d"%(m,n,k)
    all_exe_src = [basename(fn) for fn in glob(outdir+t+"_*_main.cu")]
    all_exe = sorted([fn.replace("_main.cu", "") for fn in all_exe_src])

    output  = "#!/bin/bash -l\n"
    output += "#SBATCH --nodes=%d\n"%len(all_exe)
    output += "#SBATCH --time=0:30:00\n"
    output += "#SBATCH --account=ch5\n"
    output += "\n"
    output += "source ${MODULESHOME}/init/sh;\n"
    output += "module unload PrgEnv-cray\n"
    output += "module load cudatoolkit PrgEnv-gnu\n"
    output += "module list\n"
    output += "cd $SLURM_SUBMIT_DIR \n"
    output += "\n"
    output += "date\n"
    for exe in all_exe:
        output += "aprun -b -n 1 -N 1 -d 8 make -j 16 %s &\n"%exe
    output += "wait\n"
    output += "date\n"
    output += "\n"
    for exe in all_exe:
        output += "aprun -b -n 1 -N 1 -d 1 ./"+exe+" >"+exe+".log 2>&1 & \n"
    output += "wait\n"
    output += "date\n"
    output += "\n"
    output += "echo Over all winner:\n"
    output += 'grep WINNER .'+t+'_exe*.log  |  sort -n --field-separator="#" -k 2 | tail -n 1\n'
    output += "\n"
    output += "#EOF\n"

    fn = outdir+t+".job"
    writefile(fn, output)



#===============================================================================
def gen_makefile(outdir):
    output  = ".SECONDARY:\n"
    output += "vpath %.cu ../\n\n"
    all_exe_src = sorted([basename(fn) for fn in glob(outdir+"/tune_*_main.cu")])
    build_targets = [fn.replace("_main.cu", "") for fn in all_exe_src]

    output += ".PHONY: do_nothing build_all \n\n"
    output += "do_nothing:\n\n"
    output += "build_all: " +  " ".join(build_targets) + "\n\n"

    output += "libcusmm_benchmark.o : libcusmm_benchmark.cu\n\n"

    headers = " ".join( ["."+fn for fn in glob("./kernels/*.h")] )
    output += "%.o : %.cu "+headers+"\n"
    output += "\tnvcc -arch=sm_35 -c $<\n\n"

    for exe_src in all_exe_src:
        absparts = sorted(glob(outdir+"/"+exe_src.replace("_main.cu", "_part*")))
        parts = [basename(fn) for fn in absparts]
        deps = [exe_src, "libcusmm_benchmark.cu"] + parts
        deps_obj = " ".join([fn.replace(".cu", ".o") for fn in deps])
        exe = exe_src.replace("_main.cu", "")
        output += exe + " : " + deps_obj +"\n"
        output += "\tnvcc -arch=sm_35 -o $@ $^\n\n"

    writefile(outdir+"/Makefile", output)


#===============================================================================
def gen_collect(outdir, triples):
    output  = "#!/bin/bash\n"
    for (m, n, k)  in triples:
        t = "/tune_%dx%dx%d"%(m,n,k)
        output += 'grep WINNER .'+t+'_exe*.log  |  sort -n --field-separator="#" -k 2 | tail -n 1\n'
    output += "#EOF\n"
    fn = outdir+"/collect_winners.sh"
    writefile(fn, output)
    os.system("chmod +x "+fn)


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

    #print("Wrote: "+fn)

#===============================================================================
def combinations(*sizes):
     return(list(product(sizes, sizes, sizes)))

#===============================================================================

main()

#EOF

