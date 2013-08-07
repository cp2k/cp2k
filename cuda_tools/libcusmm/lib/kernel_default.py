# -*- coding: utf-8 -*-

class DefaultKernel(object):
    def __init__(self, **args):
        self.__dict__.update(args)
        self.name = "cusmm_kernel_default_%d_%d_%d"%(self.m,self.n,self.k)

    def __repr__(self):
        return("<%s>"%self.name)

    def can_handle(self, m, n, k):
        return(self.m==m and self.n==n and self.k==k)

    def launcher_code(self):
        output  = "int launch_"+self.name+"(int *param_stack, int stack_size, "
        output += "cudaStream_t stream, int m_max, int n_max, int k_max, "
        output += "double *a_data, double *b_data, double *c_data){\n"
        output += "int shared_size = 0;\n"
        output += "int careful = (stack_size / GROUPING);\n"
        output += "int nruns = stack_size - careful * GROUPING;\n"
        output += "cusmm_kernel_default<%d,%d,%d,2> "%(self.m,self.n,self.k)
        output += "<<< ((stack_size + GROUPING - 1) / GROUPING), 192, shared_size, stream >>>\n"
        output += "(param_stack, careful, nruns, m_max, n_max, k_max, \n"
        output += "a_data, b_data, c_data);\n"
        output += "return(0);\n"
        output += "}\n"
        return(output)

    def include(self):
        return("cusmm_kernel_default.h")

#EOF
