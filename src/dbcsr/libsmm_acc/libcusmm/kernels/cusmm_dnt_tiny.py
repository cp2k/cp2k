# -*- coding: utf-8 -*-

class Kernel_dnt_tiny(object):
    def __init__(self, **params):
        self.__dict__.update(params)
        self.name  = "cusmm_dnt_tiny_"
        self.name += "_".join([str(params[k]) for k in sorted(params.keys())])
        assert(self.m * self.n <= self.threads)

    def __repr__(self):
        return("<%s>"%self.name)

    def can_handle(self, m, n, k):
        return(self.m==m and self.n==n and self.k==k)

    def include(self):
        return("cusmm_dnt_tiny.h")

    def launcher_code(self):
       output  = "int launch_"+self.name+"(int *param_stack, int stack_size, "
       output += "cudaStream_t stream, int m_max, int n_max, int k_max, "
       output += "double *a_data, double *b_data, double *c_data){\n"
       output += "int shared_size = 0;\n"
       output += "//%s\n"%str(self.__dict__)
       output += "typedef void (*kernel)(const int*, int, const double*, const double*, double*);\n"
       output += "static kernel kern_func = cusmm_dnt_tiny<%(m)d,%(n)d,%(k)d,%(threads)d,%(grouping)d,%(minblocks)d>;\n"%self.__dict__
       output += "static bool configured = false;\n"
       output += "if(configured == false){\n"
       output += "  cudaError_t err = cudaFuncSetSharedMemConfig(kern_func, cudaSharedMemBankSizeEightByte);\n"
       output += "  if(err != cudaSuccess) return(-1);\n"
       output += "  configured = true;\n"
       output += "}\n"
       output += "kern_func<<< ((stack_size + %(grouping)d - 1) / %(grouping)d), %(threads)d, shared_size, stream >>>\n"%self.__dict__
       output += "(param_stack, stack_size, \n"
       output += "a_data, b_data, c_data);\n"
       output += "return(0);\n"
       output += "}\n"
       return(output)

    @staticmethod
    def promising_parameters(m, n, k):
        params = []
        for minblocks in (7, 8, 14, 28):              # heuristic: kernel dependent optimum
            for grouping in range(1, 33, 1):          # soft: divide stack work in chunks of grouping + the rest
                for threads in (16, 32, 64):          # heuristic: not more than 2 warps per SM (sm_60)
                    if (m * n > threads):
                       continue                       # hard: not enough threads to cover result matrix
                
                    buf_sz = k * (m + n)
                    sizeof_int = 4; sizeof_double = 8
                    smem_tot = buf_sz * sizeof_double + 3 * grouping * sizeof_int
                    if (smem_tot * minblocks > 48 * 1024): # hard: see cudaFuncSetCacheConfig() docu
                       continue                       # hard: uses too much shared memory
                
                    params.append({'m':m, 'n':n, 'k':k,
                                   'threads':threads,
                                   'grouping':grouping,
                                   'minblocks':minblocks})
        return(params)

#EOF
