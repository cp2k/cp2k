/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#include "dbcsr_cuda.h"

extern __shared__ double cache[];


__global__ void stack_mm_mnk_d (
	const int *__restrict__ param_stack,
	const int careful, const int nruns,
	const int m, const int n, const int k,
	//const int mn, const int mk, const int kn, const int maxb,
	const int liter,
	const double *__restrict__ a_data,
	const double *__restrict__ b_data,
	double *__restrict__ c_data,
	int *__restrict__ c_locks) {
	
	/**
	 *  \var sp        which stack member this thread block is processing
	 (= CUDA thread block)
	 *  \var psp       pointer to first element of parameters
	 *  \var c_loc     pointer to C data
	 *  \var run       run number
         *  \var nrun      number of runs
	 *  \var my_id     my ID for locking
	 *  \var tn        thread number (of CUDA thread block)
	 *  \var mn        product of the block dimensions
	 *  \var l         multiplication loop index
	 *  \var c, r      C matrix row, column of this thread
	 *  \var myc       C matrix accumulator
	 *  \var buff_l    cache for A data
	 *  \var buff_r    cache for B data
	 *  \var c_id      translated C block number (used in locking)
	 *  \var lock_owner  current C block owner (used in locking)
	 */ 

	int lock_owner, c_id, my_id;
	const int mn = m * n;
	const int mk = m * k;
	const int kn = n * k;
	const int r = threadIdx.x % m;
	const int c = threadIdx.x / m;
	int l, i;
	double myc, tmp;
	const double * __restrict__ buff_l, * __restrict__ buff_r;

	int psp, c_loc;

	int run, nrun;

	__shared__ int our_params[7];
	double *buff;

	buff = (double *) cache;
	buff_l = buff;
	buff_r = &(buff[mk]);

	nrun = GROUPING;
	if (blockIdx.x == careful)
		nrun = nruns;

	for (run = 0; run < nrun; run ++) {
		psp = 7*(blockIdx.x*GROUPING + run);

//		for (l = 0; l <= (mk-1) / blockDim.x; l++) {
//			i = threadIdx.x+blockDim.x*l;
//			if (i < mk)
//				buff[i] = a_data[our_params[3]-1+i];
//		}
//		for (l = 0; l <= (kn-1) / blockDim.x; l++) {
//			i = threadIdx.x+blockDim.x*l;
//			if (i < kn)
//				buff[mk+i] = b_data[our_params[4]-1+i];
//		}
		for (l = 0; l <= liter; l++) {
			i = threadIdx.x+blockDim.x*l;
			if (i < mk)
				buff[i] = a_data[param_stack[psp+3]-1+i];
			if (i < kn)
				buff[mk+i] = b_data[param_stack[psp+4]-1+i];
		}

		syncthreads();

		/* Do actual multiplication. */
		if (threadIdx.x < mn) {
			myc = 0.0l;

			for (l = 0; l < k; l++) {
				myc = myc +
					buff_l[   l*m  + r] *
					buff_r[   c*k + l];
			}
		}

		/* Lock the C block. */
		c_loc = param_stack[psp+5]-1;
		//c_loc = our_params[5]-1;
		syncthreads();
		c_id = param_stack[psp+6]-1;
		//c_id = our_params[6]-1;

		if (threadIdx.x == 0) {
			my_id = blockIdx.x+1;
			lock_owner = 0;
			while ((lock_owner != my_id))
				lock_owner = atomicCAS (&(c_locks[c_id]), 0, my_id);
		} else if (threadIdx.x==1) {
			tmp = c_data[c_loc];
		}
			

		/* Add our results to the C block. */
		syncthreads();
		if (threadIdx.x < mn) {
			c_data[c_loc+threadIdx.x] += myc;
		}

		/* Release the lock on the C block. */
		syncthreads();
		if (threadIdx.x == 0) {
			c_locks[c_id] = 0;
		}
	}


};


__global__ void stack_mm_d
                   (const int *__restrict__ param_stack,
		    int stack_size, int nparams,
		    const double *__restrict__ a_data,
		    const double *__restrict__ b_data,
		    double *__restrict__ c_data,
		    int *__restrict__ c_locks) {

  /**
   *  \var sp        which stack member this thread block is processing
                     (= CUDA thread block)
   *  \var sp_one    translated stack (=sp+1)
   *  \var tn        thread number (of CUDA thread block)
   *  \var nt        number of threads (size of CUDA thread block)
   *  \var m, n, k   dimensions of the blocks (C is m*n, A is m*k, B is k*n)
   *  \var mn, mk, kn  product of the block dimensions
   *  \var l         multiplication loop index
   *  \var c, r      C matrix row, column of this thread
   *  \var myc       C matrix accumulator
   *  \var buff      cache for A and B data
   *  \var c_id      translated C block number (used in locking)
   *  \var lock_owner  current C block owner (used in locking)
   */ 

  int sp, lock_owner, c_id, sp_one;
  int tn;
  int r, c, l;
  int m, n, k;
  int mn;
  double myc;
  const double *buff_l, *buff_r;

  int psp, c_loc;


  /* Setup shared memory. */
  //buff = (double *) cache;

  /* Determine who I am. */
  sp = blockIdx.x;
  tn = threadIdx.x;

  psp = 7*sp;
  m = param_stack[psp];
  n = param_stack[psp+1];
  k = param_stack[psp+2];

  buff_l = &(a_data[param_stack[psp+3]-1]);
  buff_r = &(b_data[param_stack[psp+4]-1]);

  /* Calculate who I am. */

  mn = m*n;

  /* Do actual multiplication. */
  if (tn < mn) {
    r = tn % m;
    c = tn / m;
    myc = 0.0l;

    for (l = 0; l < k; l++) {
      myc = myc +
	buff_l[   l*m+r] *
	buff_r[   c*k+l];
    }
  }

  /* Lock the C block. */
  c_id = param_stack[psp+6]-1;
  c_loc = param_stack[psp+5]-1;
  syncthreads();
  if (tn == 0) {
    sp_one = sp + 1;
    lock_owner = 0;
    while ((lock_owner != sp_one))
      lock_owner = atomicCAS (&(c_locks[c_id]), 0, sp_one);
  }

  /* Add our results to the C block. */
  syncthreads();
  if (tn < mn) {
    c_data[c_loc+tn] += myc;
  }

  /* Release the lock on the C block. */
  syncthreads();
  if (tn == 0) {
    c_locks[c_id] = 0;
    //threadfence();
  }

};


__global__ void stack_mm_mnk_d_direct (
	const int *__restrict__ param_stack,
	const int careful, const int nruns,
	const int m, const int n, const int k, const int mn,
	const double *__restrict__ a_data,
	const double *__restrict__ b_data,
	double *__restrict__ c_data,
	int *__restrict__ c_locks) {

	/**
	 *  \var sp        which stack member this thread block is processing
	 (= CUDA thread block)
	 *  \var psp       pointer to first element of parameters
	 *  \var c_loc     pointer to C data
	 *  \var run       run number
         *  \var nrun      number of runs
	 *  \var my_id     my ID for locking
	 *  \var tn        thread number (of CUDA thread block)
	 *  \var mn        product of the block dimensions
	 *  \var l         multiplication loop index
	 *  \var c, r      C matrix row, column of this thread
	 *  \var myc       C matrix accumulator
	 *  \var buff_l    cache for A data
	 *  \var buff_r    cache for B data
	 *  \var c_id      translated C block number (used in locking)
	 *  \var lock_owner  current C block owner (used in locking)
	 */ 

	int lock_owner, c_id, my_id;
	int l;
	const int r = threadIdx.x % m;
	const int c = threadIdx.x / m;
	double myc, tmp;
	const double *buff_l, *buff_r;

	int psp, c_loc;

	int run, nrun;

	nrun = GROUPING;
	if (blockIdx.x == careful)
		nrun = nruns;

	for (run = 0; run < nrun; run ++) {
		psp = 7*(blockIdx.x*GROUPING + run);

		buff_l = &(a_data[param_stack[psp+3]-1]);
		buff_r = &(b_data[param_stack[psp+4]-1]);
		/* Do actual multiplication. */
		if (threadIdx.x < mn) {
			myc = 0.0l;

			for (l = 0; l < k; l++) {
				myc = myc +
					buff_l[   l*m+r] *
					buff_r[   c*k+l];
			}
		}

		/* Lock the C block. */
		c_loc = param_stack[psp+5]-1;
		syncthreads();
		c_id = param_stack[psp+6]-1;

		if (threadIdx.x == 0) {
			my_id = blockIdx.x+1;
			lock_owner = 0;
			while ((lock_owner != my_id))
				lock_owner = atomicCAS (&(c_locks[c_id]), 0, my_id);
		} else if (threadIdx.x==1) {
			tmp = c_data[c_loc];
		}
			

		/* Add our results to the C block. */
		syncthreads();
		if (threadIdx.x < mn) {
			c_data[c_loc+threadIdx.x] += myc;
		}

		/* Release the lock on the C block. */
		syncthreads();
		if (threadIdx.x == 0) {
			c_locks[c_id] = 0;
			//threadfence();
		}
	}


};


__global__ void stack_mm_mnk_vec_d (
	const int *__restrict__ param_stack,
	const int stack_size, const int nmat,
	const int m, const int n, const int k, const int mn,
	const double *__restrict__ a_data,
	const double *__restrict__ b_data,
	double *__restrict__ c_data,
	int *__restrict__ c_locks) {
	
	/**
	 *  \var sp        which stack member this thread block is processing
	 (= CUDA thread block)
	 *  \var psp       pointer to first element of parameters
	 *  \var c_loc     pointer to C data
	 *  \var run       run number
         *  \var nrun      number of runs
	 *  \var my_id    translated stack (=sp+1)
	 *  \var tn        thread number (of CUDA thread block)
	 *  \var mn        product of the block dimensions
	 *  \var l         multiplication loop index
	 *  \var c, r      C matrix row, column of this thread
	 *  \var myc       C matrix accumulator
	 *  \var buff_l    cache for A data
	 *  \var buff_r    cache for B data
	 *  \var c_id      translated C block number (used in locking)
	 *  \var lock_owner  current C block owner (used in locking)
	 */ 

	int lock_owner, c_id, my_id;
	const int tn = threadIdx.x;
	int nmat_used;
	int nt;
	const int r = threadIdx.x % m;
	int c, l;
	double myc[32];
	double mya[32];
	__shared__ int our_b[32];
	const double *buff_l, *buff_r;

	int psp, c_loc;
	int run, nrun;
	const int my_mat_num = threadIdx.x / m;
	int imat;

	//nrun = GROUPING;
	//if ((blockIdx.x+1) * GROUPING > stack_size)
	//	nrun = stack_size - (blockIdx.x)*GROUPING;

	nmat_used = nmat;
	if ((blockIdx.x+1)*nmat > stack_size)
		nmat_used = stack_size - (blockIdx.x)*nmat;
	nt = m * nmat_used;

	//for (run = 0; run < nrun; run ++) {
	//sp = blockIdx.x*GROUPING + run;

	psp = 7*(blockIdx.x*nmat + my_mat_num);

	buff_l = &(a_data[param_stack[psp+3]-1]);
	buff_r = &(b_data[param_stack[psp+4]-1]);

	/* Do actual multiplication. */
	if (tn < nt) {
		for (l = 0; l < k; l++) {
			mya[l] = buff_l[ l*m + r ];
		}
		for (c = 0; c < n; c++) {
			if (tn < k)
				our_b[l] = buff_r[c*k+tn];
			syncthreads();
			myc[c] = 0.0l;
		
			for (l = 0; l < k; l++) {
				myc[c] = myc[c] +
					mya   [   l    ] *
					our_b [   l    ];
				//buff_r[   c*k+l];
			}
		}
	}
	/* Lock the C block. */
	c_id = param_stack[psp+6]-1;
	syncthreads();
	c_loc = param_stack[psp+5]-1;
	my_id = blockIdx.x + 1;
	for (imat = 0; imat < nmat_used; imat++) {
		if (r == 0 && imat == my_mat_num) {
			lock_owner = 0;
			while ((lock_owner != my_id))
				lock_owner = atomicCAS (&(c_locks[c_id]), 0, my_id);
		}

		/* Add our results to the C block. */
		syncthreads();
		if (tn < nt && imat == my_mat_num) {
			for (c = 0; c < n; c++) {
				c_data[c_loc+r+c*m] += myc[c];
			}
		}

		/* Release the lock on the C block. */
		syncthreads();
		if (r == 0 && imat == my_mat_num) {
			c_locks[c_id] = 0;
		}
	}
};
