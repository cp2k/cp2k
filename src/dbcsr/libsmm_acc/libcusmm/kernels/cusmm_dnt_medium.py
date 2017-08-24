# -*- coding: utf-8 -*-

from math import ceil

class Kernel_dnt_medium(object):
    def __init__(self, **params):
        self.__dict__.update(params)
        self.name  = "cusmm_dnt_medium_"
        self.name += "_".join([str(params[k]) for k in sorted(params.keys())])
        assert(self.threads * self.minblocks <= 2048)
        min_threads = ((self.m + self.tile_m - 1) // self.tile_m) * ((self.n + self.tile_n - 1) // self.tile_n)
        assert(min_threads <= self.threads)

    def __repr__(self):
        return("<%s>"%self.name)

    def can_handle(self, m, n, k):
        return(self.m==m and self.n==n and self.k==k)

    def include(self):
        return("cusmm_dnt_medium.h")

    def launcher_code(self):
       output  = "int launch_"+self.name+"(int *param_stack, int stack_size, "
       output += "cudaStream_t stream, int m_max, int n_max, int k_max, "
       output += "double *a_data, double *b_data, double *c_data){\n"
       output += "int shared_size = 0;\n"
       output += "//%s\n"%str(self.__dict__)
       output += "typedef void (*kernel)(const int*, int, const double*, const double*, double*);\n"
       output += "static kernel kern_func = cusmm_dnt_medium<%(m)d,%(n)d,%(k)d,%(tile_m)d,%(tile_n)d,%(threads)d,%(grouping)d,%(minblocks)d>;\n"%self.__dict__
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
        for minblocks in range(1, 28, 1):
            for grouping in range(1, 33, 1):
                for threads in range (32, 257, 32):
                    if(threads * minblocks > 2048): # hard: too much concurrent threads per SM
                        continue
                    for tm in range(1, 7):
                        for tn in range(1, 7):
                            if (tm * tn > 16):
                                continue #heuristic:
                            min_threads = ((m + tm - 1) // tm) * ((n + tn - 1) // tn)
                            if (min_threads > threads):
                                continue #hard: not enough threads to cover result matrix
    
                            if (threads > 4 * min_threads):
                                continue #heuristic: too many threads unused during calculation
    
                            cmax = ((n + tn - 1) // tn)
                            rmax = ((m + tm - 1) // tm)
                            buf_sz = max(m * n, m * k + k * tn * cmax, tm * rmax * k + 1)
                            sizeof_int = 4; sizeof_double = 8
                            smem_tot = buf_sz * sizeof_double + 3 * grouping * sizeof_int
                            if(smem_tot * minblocks > 48 * 1024): # hard: see cudaFuncSetCacheConfig() docu
                                continue #hard: uses too much shared memory
    
                            params.append({'m':m, 'n':n, 'k':k,
                                           'tile_m':tm, 'tile_n':tn,
                                           'threads':threads,
                                           'grouping':grouping,
                                           'minblocks':minblocks})
        return(params)


#EOF
