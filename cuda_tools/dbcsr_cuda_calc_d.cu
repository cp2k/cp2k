/******************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations
 *  Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group
 *****************************************************************************/

#include "dbcsr_cuda.h"

extern __shared__ double cache[];

/* The following are defined in "dbcsr_cuda.h" */
/* SQUARESIZE, BLOCK, TDIM, NUMTHREADS23SQ, BLOCKSPLIT23SQ, TOTALTHREADS23SQ */

/******************************************************************************/
/* The stack_mm_mnk_sq23_d routine was added by Neil Stringfellow */
/* This is a specialised routine for square matrices of size 23 (m==n==k==23) */
/******************************************************************************/

/* The following is only valid for square matrices of size 23 (m==n==k==23) */
/* Should be called with 64 threads only (blockDim.x==64) */
/* N.B. This is checked for in the calling routine, not down here */

/* Why 64 threads I here you ask ?  What about occupancy !!!!! */
/* !!! Update - using BLOCKSPLIT23SQ we can have 128 threads !!! */
/* I'm glad you asked ... */
/* With 64 threads we can provide the device with lots of work, and with a small
   number of threads we have a better chance of having multiple blocks being
   executed concurrently so that they will overlap core floating point computation
   and global memory traffic, hiding the memory bandwidth pressure that we saw in
   the original implemenation.
   The number of thread-blocks that can be executed concurrently is probably more
   driven by the shared memory requirements per thread-block which is a little
   over 8Kbytes, and therefore we should be able to get up to 5 blocks working at
   the same time per SM.
   If indeed we can have a maximum of 5 thread-blocks each of 64 threads working
   concurrently then we have at most 320 threads, and these threads then can use
   the full complement of 63 32-bit registers available as a maximum per thread.
   This then means that we can hold up to 31 double precision numbers in registers
   and allows the inner kernel (see the code) to operate C=Matmul(A,B) on 3x3
   square matrices entirely out of registers (3 arrays each of 3x3 = 27 registers
   required for the inner kernel */

/* Needs 24x24*8*2 bytes of shared memory when being called */
__global__ void stack_mm_mnk_sq23_d (
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

	/* We have received the data from a Fortran code, so we should try to keep
	   Fortran ordering when considering rows and columns */

	/* Each thread has a row and column beginning number.
	   In each case we have an "r" value which is defined as how many rows we move down the matrix to get our
	   starting point, and this is then multiplied by blocks of (BLOCKS/TDIM) for the next entry.
	   This has changed from the original implementation in case it helps with bank conflicts.
	   For the starting column we have to move along according to TDIM, and here "c" will be
	   later multiplied by TDIM to find the starting point.
	*/
	const int r = (threadIdx.x % NUMTHREADS23SQ) % (BLOCK/TDIM);
	const int c = (threadIdx.x % NUMTHREADS23SQ) / (BLOCK/TDIM);
	const int whichhalf=threadIdx.x/(NUMTHREADS23SQ);
	const int blockbegin=((BLOCK/BLOCKSPLIT23SQ)*whichhalf)/(TDIM);
	const int a_offset_base=r;
	const int b_offset_base=TDIM*c*BLOCK;
	const int a_offset=a_offset_base + blockbegin*TDIM*BLOCK;
	const int b_offset=b_offset_base + blockbegin*TDIM;

	int l, i;

	/* We convert myc into mycarr, a TDIMxTDIM array (i.e. mycarr = myc-array) */
	double mycarr[TDIM][TDIM];

	const double * __restrict__ buff_a, * __restrict__ buff_b;
       
	int psp, c_loc;

	int run, nrun;

	double *buff;

	int j, ki;

	int myarrayindex;
	int buffaddr;

	buff = (double *) cache;

	/* We're going to use a 24x24 block of shared memory to store our 23x23 A and B matrices */
	/* This will allow us to give an independent 3x3 block to each of 64 threads in a 8x8 configuration */
	/* The edge threads will do some redundant calculations, but this is less important than
	   having a balanced workload and clean code. The actual flops are about 9% more than the flops we
	   want, but the threads would be scheduled to do something anyway, and we would need to sprinkle
	   if statements throughout the code. */
	/* So we have buff_a as an alias for the first 24x24 (x8byte double precision) block of shared
	   memory and buff_b as the second part */
	buff_a = buff;
	buff_b = &(buff[BLOCK*BLOCK]);

	nrun = GROUPING;
	if (blockIdx.x == careful)
		nrun = nruns;

	/* First let's zero the 24th row and column of the A and B shared memories */
	if(threadIdx.x<BLOCK){
	  buff[BLOCK*(BLOCK-1)+threadIdx.x]=0.0l;
	  buff[BLOCK*BLOCK + BLOCK*(BLOCK-1)+threadIdx.x]=0.0l;
	  buff[BLOCK*threadIdx.x+(BLOCK-1)]=0.0l;
	  buff[BLOCK*BLOCK + BLOCK*threadIdx.x+(BLOCK-1)]=0.0l;
	}

	/* Set the partial sums to zero. This used to be done in the inner loop, but now we will typically
	   carry the mycarr values over "run" loops as we only update the C matrix when we have run through
	   all iterations in this thread block that would be working on the same piece of C. */
	for(i=0;i<TDIM;i++){
	  for(j=0;j<TDIM;j++){
	    mycarr[i][j] = 0.0l;
	  }
	}

	for (run = 0; run < nrun; run ++) {
		psp = 7*(blockIdx.x*GROUPING + run);

		/* Load from main memory  and store into A and B */
		/* We won't get coallesced accesses on 128-byte boundaries, but we should be close to the
		   optimal performance since we are contiguous in memory and therefore the remainder of
		   each 128-byte boundary should have been loaded in L2 cache. */
		/* We need to map from the 23x23 arrays in a_data and b_data to a 24x24 block that 
		 we will use in shared memory, with the 24th row and column being zero. */

		{
		  const int a_arr_base_index=param_stack[psp+3]-1;
		  const int b_arr_base_index=param_stack[psp+4]-1;

		  for (l=0;l<(SQUARESIZE*SQUARESIZE)/TOTALTHREADS23SQ+1;l++){
		    int myarrayindex=l*TOTALTHREADS23SQ+threadIdx.x;
		    int buffaddr=myarrayindex + (myarrayindex)/SQUARESIZE;
		    if(myarrayindex<(SQUARESIZE*SQUARESIZE)){
		      /* Copy A array */
		      buff[buffaddr]=a_data[a_arr_base_index+myarrayindex];
		      /* Copy B array */
		      buff[(BLOCK*BLOCK)+buffaddr]=b_data[b_arr_base_index+myarrayindex];
		    }
		  }
		  
		}

		syncthreads();

		/* Do multiplication in 3x3 blocks over the appropriate rows/columns of A and B */
		for (l = 0; l < ((BLOCK)/((TDIM)*(BLOCKSPLIT23SQ))); l++) {
		  /* Declare two arrays for a and b that we hope the compiler will place in registers */
		  /* Actually we don't need arrays strictly in the algorithm so let's declare them a
		     little differently - looking at the assembler output with unrolling you get the same
		     effect with either configuration anyway. */
		  double a_reg_scalar;
		  double b_regs1D[TDIM];

		  /* Load the A and B values into registers and compute using those registers */
		  /* This is a rearrangement from the original just to show that the register load and
		     compute can be done together */
#pragma unroll
		  for(ki=0;ki<TDIM;ki++){
#pragma unroll
		    for (i=0;i<TDIM;i++){
		      a_reg_scalar=buff_a[   ((l*TDIM)+ki)*BLOCK  + a_offset + i*(BLOCK/TDIM) ];
#pragma unroll
		      for (j=0;j<TDIM;j++){
			if(i==0){
			  b_regs1D[j]=buff_b[   b_offset + j*BLOCK + (l*TDIM) + ki ];
			}
			  mycarr[i][j] += a_reg_scalar * b_regs1D[j];
		      }
		    }
		  }
		}

		/* Only update c_data if we are in the last iteration, or if the next C-block
		   will be different to this C-block */
		/* param_stack[psp+6] is the current C-block ID, so adding 7 means that param_stack[psp+6+7]
		   should be the next C-block ID */
		if(run==nrun-1 || param_stack[psp+6]-1 != param_stack[psp+6+7]-1) {
		  c_loc = param_stack[psp+5]-1;
		  c_id = param_stack[psp+6]-1;
		  
		  if (threadIdx.x == 0) {
		    my_id = blockIdx.x+1;
		    lock_owner = 0;
		    while ((lock_owner != my_id))
		      lock_owner = atomicCAS (&(c_locks[c_id]), 0, my_id);
		  } 
		  
		  /* Here we need to treat the threads differently depending upon whether they are 
		     in the lower half or the upper half of the thread block in the case of a 128-thread
		     thread block. All threads in the same warp will have the same value for "whichhalf" so
		     hopefully this should work smoothly. */
		  for(ki=0;ki<BLOCKSPLIT23SQ;ki++){
		  
		  /* Need to have finished with A in order to reuse shared memory used for A for a
		     temporary store, so we need a sync */
		    syncthreads();
		    
		    if(ki==whichhalf){
		      
		      /* Add our results into a temporary storage in buff that is normally used for A.
			 This forms our C block in the same 24x24 form as for A and B in shared memory.
			 As we need to refresh A on every run through the "run" loop, we can reuse buff
			 for storing this copy of C. */
		      if(ki==0){
#pragma unroll
			for (i=0;i<TDIM;i++){
#pragma unroll
			  for (j=0;j<TDIM;j++){
			    buff[a_offset_base + b_offset_base + j*BLOCK + i*(BLOCK/TDIM) ]=mycarr[i][j];
			  }
			}
			
		      }else{
			
#pragma unroll
			for (i=0;i<TDIM;i++){
#pragma unroll
			  for (j=0;j<TDIM;j++){
			    buff[a_offset_base + b_offset_base + j*BLOCK + i*(BLOCK/TDIM) ]+=mycarr[i][j];
			  }
			}
			
		      }
		    }
		  }

		  /* We need to ensure that we have a coherent copy of C in the buffer, so that
		     means we have another sync */
		  syncthreads();

		  /* Need to reverse engineer back where these go into memory */
		  /* Need to effect coallesced accesses so as when loading A and B earlier, we use
		   a mapping of the 23x23 C block onto our 24x24 buff block where row 24 and column
		   24 are zeroes so they can be ignored. */
		  for (l=0;l<=(SQUARESIZE*SQUARESIZE)/TOTALTHREADS23SQ+1;l++){
		    myarrayindex=l*TOTALTHREADS23SQ+threadIdx.x;
		    if(myarrayindex<(SQUARESIZE*SQUARESIZE)){
		      buffaddr=myarrayindex + (myarrayindex)/SQUARESIZE;
		      c_data[c_loc+myarrayindex]+=buff[buffaddr];
		    }
		  }

		  /* Release the lock on the C block. */
		  syncthreads();
		  if (threadIdx.x == 0) {
		    c_locks[c_id] = 0;
		  }
		  /* If we have another C-block then we need to reset our partial sum to zero for the new C-block */
		  if(run!=nrun-1){
		    for(i=0;i<TDIM;i++){
		      for(j=0;j<TDIM;j++){
			mycarr[i][j] = 0.0l;
		      }
		    }
		  }

		}

		syncthreads();

	}

};

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
	double myc;
	const double * __restrict__ buff_l, * __restrict__ buff_r;
       
	int psp, c_loc;

	int run, nrun;

	double *buff;

	buff = (double *) cache;
	buff_l = buff;
	buff_r = &(buff[mk]);

	nrun = GROUPING;
	if (blockIdx.x == careful)
		nrun = nruns;

	/* Set the partial sum to zero (this used to be done in the inner loop, but now we might carry it over loops) */
	myc = 0.0l;

	for (run = 0; run < nrun; run ++) {
		psp = 7*(blockIdx.x*GROUPING + run);

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

			for (l = 0; l < k; l++) {
				myc = myc +
				  buff_l[   l*m  + r] *
				  buff_r[   c*k + l ];
			}

		}

		/* Only update c_date if we are in the last iteration, or if the next C-block will be
		   different to this C-block */
		/* param_stack[psp+6] is the current C-block ID, so adding 7 means that param_stack[psp+6+7]
		   should be the next C-block ID */
		if(run==nrun-1 || param_stack[psp+6]-1 != param_stack[psp+6+7]-1) {
		  c_loc = param_stack[psp+5]-1;
		  c_id = param_stack[psp+6]-1;
		  
		  if (threadIdx.x == 0) {
		    my_id = blockIdx.x+1;
		    lock_owner = 0;
		    while ((lock_owner != my_id))
		      lock_owner = atomicCAS (&(c_locks[c_id]), 0, my_id);
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
		  /* If we have another C-block then we need to reset our partial sum to zero for the new C-block */
		  myc = 0.0l;
		}

		syncthreads();

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
	double myc ;
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
//	int run, nrun;
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
