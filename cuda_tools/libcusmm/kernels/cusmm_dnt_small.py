# -*- coding: utf-8 -*-

class Kernel_dnt_small(object):
    def __init__(self, **params):
        self.__dict__.update(params)
        self.name  = "cusmm_dnt_small_"
        self.name += "_".join([str(params[k]) for k in sorted(params.keys())])

    def __repr__(self):
        return("<%s>"%self.name)

    def can_handle(self, m, n, k):
        return(self.m==m and self.n==n and self.k==k)

    def include(self):
        return("cusmm_dnt_small.h")

    def launcher_code(self):
       output  = "int launch_"+self.name+"(int *param_stack, int stack_size, "
       output += "cudaStream_t stream, int m_max, int n_max, int k_max, "
       output += "double *a_data, double *b_data, double *c_data){\n"
       output += "int shared_size = 0;\n"
       output += "//%s\n"%str(self.__dict__)
       output += "int careful = (stack_size / %(grouping)d);\n"%self.__dict__
       output += "int nruns = stack_size - careful * %(grouping)d;\n"%self.__dict__
       output += "cusmm_dnt_small<%(m)d,%(n)d,%(k)d,%(tile_m)d,%(tile_n)d,%(grouping)d> "%self.__dict__
       output += "<<< ((stack_size + %(grouping)d - 1) / %(grouping)d), %(threads)d, shared_size, stream >>>\n"%self.__dict__
       output += "(param_stack, careful, nruns, \n"
       output += "a_data, b_data, c_data);\n"
       output += "return(0);\n"
       output += "}\n"
       return(output)

#EOF
