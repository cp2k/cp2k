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
       output += "int careful = (stack_size / %(grouping)d);\n"%self.__dict__
       output += "int nruns = stack_size - careful * %(grouping)d;\n"%self.__dict__
       output += "cusmm_dnt_tiny<%(m)d,%(n)d,%(k)d,%(split_thread)d,%(threads)d,%(grouping)d,%(minblocks)d> "%self.__dict__
       output += "<<< ((stack_size + %(grouping)d - 1) / %(grouping)d), %(threads)d, shared_size, stream >>>\n"%self.__dict__
       output += "(param_stack, careful, nruns, \n"
       output += "a_data, b_data, c_data);\n"
       output += "return(0);\n"
       output += "}\n"
       return(output)

    @staticmethod
    def promising_parameters(m, n, k):
        params = []
        grouping = 16
        minblocks = 1
        for threads in (64, 96, 128):
            if(m*n > threads):
                continue

            buf_sz = 2*(m*k + k*n)
            sizeof_int = 4; sizeof_double = 8
            smem_tot = buf_sz*sizeof_double + 4*grouping*sizeof_int
            if(smem_tot*minblocks > 48*1024):
                continue # uses too much shared memory

            params.append({'m':m, 'n':n, 'k':k,
                           'threads':threads,
                           'split_thread':32,
                           'grouping':grouping,
                           'minblocks':minblocks})
        return(params)

#EOF
