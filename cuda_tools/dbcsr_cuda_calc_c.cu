/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group
 *****************************************************************************/


__global__ void stack_mm_c
                   (const int* __restrict__ param_stack,
		    int stack_size, int nparams,
		    const float* __restrict__ a_data,
		    const float* __restrict__ b_data,
		    float* __restrict__ c_data,
		    int* __restrict__ c_locks) {
  /**
   *  \var sp        which stack member this thread block is processing
   *  \var our_params  cache for this thread block's multiplication parameters
   *  \var m, n, k   dimensions of the blocks (C is m*n, A is m*k, B is k*n)
   *  \var mn, mk, kn  product of the block dimensions
   *  \var l         multiplication loop index
   *  \var cr        C matrix element (fortran-style) of this thread
   *  \var c, r      C matrix row, column of this thread
   *  \var myc       C matrix accumulator
   *  \var buff      cache for A and B data
   */ 
  int sp, lock_owner, c_id, sp_one;
  int r, c, l;
  int m, n, k;
  int mn, cr, mk, kn;
  float myc_r, myc_i;
  __shared__ int our_params[7];
  float *buff;


  /* Setup shared memory. */
  buff = (float *) cache;

  /* Determine who I am. */
  sp = blockIdx.x;
  cr = threadIdx.x;

  /* Load in the parameters. */
  if (blockDim.x >= 7) {
    if (cr < 7) {
      our_params[cr] = param_stack[cr+7*sp];
    }
  } else if (cr == 0) {
    for (l = 0; l < 7; l++) {
      our_params[l] = param_stack[l+7*sp];
    }
  }
  syncthreads();
  m = our_params[0];
  n = our_params[1];
  k = our_params[2];

  /* Load in the buffers.  The first mk threads load in A while the
     last kn threads load in B. */
  mk = m*k;
  kn = k*n;
  if (cr < mk) {
    buff[2*cr  ] = a_data[2*(our_params[3]-1+cr)  ];
    buff[2*cr+1] = a_data[2*(our_params[3]-1+cr)+1];
  }
  if (cr >= blockDim.x - kn) {
    buff[2*(mk+(cr-(blockDim.x-kn)))  ] = b_data[2*(our_params[4]-1+(cr-(blockDim.x-kn)))];
    buff[2*(mk+(cr-(blockDim.x-kn)))+1] = b_data[2*(our_params[4]-1+(cr-(blockDim.x-kn)))+1];
  }

  /* Calculate who I am. */
  syncthreads();

  mn = m*n;

  /* Do actual multiplication. */
  if (cr < mn) {
    r = cr % m;
    c = cr / m;
    myc_r = 0.0f;
    myc_i = 0.0f;

    for (l = 0; l < k; l++) {
      myc_r = myc_r +
	buff[2*(   l*m+r)  ] *
	buff[2*(mk+c*k+l)  ] -
	buff[2*(   l*m+r)+1] *
	buff[2*(mk+c*k+l)+1];
      myc_i = myc_i +
	buff[2*(   l*m+r)  ] *
	buff[2*(mk+c*k+l)+1] +
	buff[2*(   l*m+r)+1] *
	buff[2*(mk+c*k+l)  ];
    }
  }

  /* Lock the C block. */
  syncthreads();
  if (cr == 0) {
    sp_one = sp + 1;
    c_id = our_params[6]-1;
    lock_owner = 0;
    while ((lock_owner != sp_one))
      lock_owner = atomicCAS (&(c_locks[c_id]), 0, sp_one);
  }

  /* Add our results to the C block. */
  syncthreads();
  if (cr < mn) {
    c_data[2*(our_params[5]-1+cr)  ] += myc_r;
    c_data[2*(our_params[5]-1+cr)+1] += myc_i;
  }

  /* Release the lock on the C block. */
  syncthreads();
  if (cr == 0)
    c_locks[c_id] = 0;

};
