# -*- coding: utf-8 -*-

from math import ceil

class Kernel_dnt_largeDB(object):
    def __init__(self, **params):
        self.__dict__.update(params)
        self.name  = "cusmm_dnt_largeDB_"
        self.name += "_".join([str(params[k]) for k in sorted(params.keys())])
        assert(self.threads * self.minblocks <= 2048)
        min_threads = ((self.m+self.tile_m-1)/self.tile_m) * ((self.n+self.tile_n-1)/self.tile_n)
        assert(min_threads <= self.threads)
        assert(self.tile_m <= self.v)
        assert(self.tile_n <= self.w)

    def __repr__(self):
        return("<%s>"%self.name)

    def can_handle(self, m, n, k):
        return(self.m==m and self.n==n and self.k==k)

    def include(self):
        return("cusmm_dnt_largeDB.h")

    def launcher_code(self):
       output  = "int launch_"+self.name+"(int *param_stack, int stack_size, "
       output += "cudaStream_t stream, int m_max, int n_max, int k_max, "
       output += "double *a_data, double *b_data, double *c_data){\n"
       output += "int shared_size = 0;\n"
       output += "//%s\n"%str(self.__dict__)
       output += "int careful = (stack_size / %(grouping)d);\n"%self.__dict__
       output += "int nruns = stack_size - careful * %(grouping)d;\n"%self.__dict__
       output += "cusmm_dnt_largeDB<%(m)d,%(n)d,%(k)d,%(tile_m)d,%(tile_n)d,%(w)d,%(v)d,%(threads)d,%(grouping)d,%(minblocks)d> "%self.__dict__
       output += "<<< ((stack_size + %(grouping)d - 1) / %(grouping)d), %(threads)d, shared_size, stream >>>\n"%self.__dict__
       output += "(param_stack, careful, nruns, a_data, b_data, c_data);\n"
       output += "return(0);\n"
       output += "}\n"
       return(output)


    @staticmethod
    def promising_parameters(m, n, k):
        params = []
        grouping = 16

        for minblocks in (1,4,8,12):
            for threads in range (96, 513, 32):

                # invalid: too many threads per SM
                if(threads * minblocks > 2048): continue

                for tm in range(1,7):
                    for tn in range(1,7):

                        # invalid: not enough threads to cover result matrix
                        min_threads = ((m+tm-1)/tm) * ((n+tn-1)/tn)
                        if(min_threads > threads): continue

                        #heuristic: too many threads unused during calculation
                        if(threads > 4*min_threads): continue

                        # heuristic: only even numbers
                        for w in range(2, k+1, 2):

                            # invalid: input slap too small
                            if (w < tn): continue

                            # heuristic: do at least one double-buffering step
                            if(2*w > k): continue

                            # heuristic: only even numbers
                            for v in range(2, n+1, 2):

                                # invalid: output slap too small
                                if (v < tm): continue

                                #heuristic: too many registers used
                                n_regs = tm*tn + (w*m+threads-1)/threads + (w*n+threads-1)/threads
                                if(n_regs*threads*minblocks > 15000): continue

                                # invalid: uses too much shared memory
                                buf_sz = max(m*w + w*n, v*m)
                                sizeof_int = 4; sizeof_double = 8
                                smem_tot = buf_sz*sizeof_double + 4*grouping*sizeof_int
                                if(smem_tot*minblocks > 48*1024): continue

                                params.append({'m':m, 'n':n, 'k':k,
                                               'tile_m':tm, 'tile_n':tn,
                                               'w':w, 'v':v,
                                               'threads':threads,
                                               'grouping':grouping,
                                               'minblocks':minblocks})
        return(params)


#EOF
