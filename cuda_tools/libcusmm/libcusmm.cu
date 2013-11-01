/******************************************************************************
*  CP2K: A general program to perform molecular dynamics simulations
*  Copyright (C) 2000 - 2013 the CP2K developers group
*****************************************************************************/
#include "./kernels/cusmm_dnt_tiny.h"
#include "./kernels/cusmm_dnt_medium.h"
#include "./kernels/cusmm_dnt_largeDB.h"
#include "./kernels/cusmm_dnt_small.h"
#include "./kernels/cusmm_common.h"
#include "./kernels/cusmm_transpose.h"


int launch_cusmm_dnt_medium_16_13_13_26_96_0_96_2_2(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'name': 'cusmm_dnt_medium_16_13_13_26_96_0_96_2_2', 'panel_in': 96, 'tile_n': 2, 'tile_m': 2, 'm': 13, 'n': 26, 'threads': 96, 'k': 13, 'panel_out': 0, 'grouping': 16}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_medium<13,26,13,2,2,96,0,16 > <<< ((stack_size + 16 - 1) / 16), 96, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_largeDB_96_16_23_23_23_2_3_13_4(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'name': 'cusmm_dnt_largeDB_96_16_23_23_23_2_3_13_4', 'blockdim': 96, 'tile_n': 3, 'tile_m': 2, 'm': 23, 'n': 23, 'w': 4, 'v': 13, 'k': 23, 'grouping': 16}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_largeDB<23,23,23,2,3,4,13,96,16> <<< ((stack_size + 16 - 1) / 16), 96, shared_size, stream >>>
(param_stack, careful, nruns, a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_medium_16_13_13_13_96_16_96_2_2(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'name': 'cusmm_dnt_medium_16_13_13_13_96_16_96_2_2', 'panel_in': 96, 'tile_n': 2, 'tile_m': 2, 'm': 13, 'n': 13, 'threads': 96, 'k': 13, 'panel_out': 16, 'grouping': 16}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_medium<13,13,13,2,2,96,16,16 > <<< ((stack_size + 16 - 1) / 16), 96, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_tiny_16_13_5_5_32_128(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'split_thread': 32, 'threads': 128, 'name': 'cusmm_dnt_tiny_16_13_5_5_32_128', 'k': 13, 'grouping': 16, 'm': 5, 'n': 5}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_tiny<5,5,13,32,16> <<< ((stack_size + 16 - 1) / 16), 128, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_largeDB_128_16_26_26_26_2_3_10_4(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'name': 'cusmm_dnt_largeDB_128_16_26_26_26_2_3_10_4', 'blockdim': 128, 'tile_n': 3, 'tile_m': 2, 'm': 26, 'n': 26, 'w': 4, 'v': 10, 'k': 26, 'grouping': 16}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_largeDB<26,26,26,2,3,4,10,128,16> <<< ((stack_size + 16 - 1) / 16), 128, shared_size, stream >>>
(param_stack, careful, nruns, a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_medium_16_13_26_26_128_0_128_2_3(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'name': 'cusmm_dnt_medium_16_13_26_26_128_0_128_2_3', 'panel_in': 128, 'tile_n': 3, 'tile_m': 2, 'm': 26, 'n': 26, 'threads': 128, 'k': 13, 'panel_out': 0, 'grouping': 16}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_medium<26,26,13,2,3,128,0,16 > <<< ((stack_size + 16 - 1) / 16), 128, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_small_16_5_13_5_96_1_1(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'threads': 96, 'name': 'cusmm_dnt_small_16_5_13_5_96_1_1', 'tile_n': 1, 'tile_m': 1, 'grouping': 16, 'm': 13, 'k': 5, 'n': 5}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_small<13,5,5,1,1,16> <<< ((stack_size + 16 - 1) / 16), 96, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_small_16_5_5_13_96_1_1(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'threads': 96, 'name': 'cusmm_dnt_small_16_5_5_13_96_1_1', 'tile_n': 1, 'tile_m': 1, 'grouping': 16, 'm': 5, 'k': 5, 'n': 13}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_small<5,13,5,1,1,16> <<< ((stack_size + 16 - 1) / 16), 96, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_largeDB_96_16_26_26_13_2_2_13_3(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'name': 'cusmm_dnt_largeDB_96_16_26_26_13_2_2_13_3', 'blockdim': 96, 'tile_n': 2, 'tile_m': 2, 'm': 26, 'n': 13, 'w': 3, 'v': 13, 'k': 26, 'grouping': 16}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_largeDB<26,13,26,2,2,3,13,96,16> <<< ((stack_size + 16 - 1) / 16), 96, shared_size, stream >>>
(param_stack, careful, nruns, a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_small_16_13_5_13_96_1_1(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'threads': 96, 'name': 'cusmm_dnt_small_16_13_5_13_96_1_1', 'tile_n': 1, 'tile_m': 1, 'grouping': 16, 'm': 5, 'k': 13, 'n': 13}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_small<5,13,13,1,1,16> <<< ((stack_size + 16 - 1) / 16), 96, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_small_16_5_13_13_96_1_2(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'threads': 96, 'name': 'cusmm_dnt_small_16_5_13_13_96_1_2', 'tile_n': 2, 'tile_m': 1, 'grouping': 16, 'm': 13, 'k': 5, 'n': 13}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_small<13,13,5,1,2,16> <<< ((stack_size + 16 - 1) / 16), 96, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_tiny_16_5_5_5_32_64(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'split_thread': 32, 'threads': 64, 'name': 'cusmm_dnt_tiny_16_5_5_5_32_64', 'k': 5, 'grouping': 16, 'm': 5, 'n': 5}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_tiny<5,5,5,32,16> <<< ((stack_size + 16 - 1) / 16), 64, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_medium_16_13_26_13_128_0_128_2_2(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'name': 'cusmm_dnt_medium_16_13_26_13_128_0_128_2_2', 'panel_in': 128, 'tile_n': 2, 'tile_m': 2, 'm': 26, 'n': 13, 'threads': 128, 'k': 13, 'panel_out': 0, 'grouping': 16}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_medium<26,13,13,2,2,128,0,16 > <<< ((stack_size + 16 - 1) / 16), 128, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_largeDB_128_16_26_13_26_1_3_26_4(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'name': 'cusmm_dnt_largeDB_128_16_26_13_26_1_3_26_4', 'blockdim': 128, 'tile_n': 3, 'tile_m': 1, 'm': 13, 'n': 26, 'w': 4, 'v': 26, 'k': 26, 'grouping': 16}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_largeDB<13,26,26,1,3,4,26,128,16> <<< ((stack_size + 16 - 1) / 16), 128, shared_size, stream >>>
(param_stack, careful, nruns, a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_medium_16_26_13_13_128_0_128_2_2(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'name': 'cusmm_dnt_medium_16_26_13_13_128_0_128_2_2', 'panel_in': 128, 'tile_n': 2, 'tile_m': 2, 'm': 13, 'n': 13, 'threads': 128, 'k': 26, 'panel_out': 0, 'grouping': 16}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_medium<13,13,26,2,2,128,0,16 > <<< ((stack_size + 16 - 1) / 16), 128, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_small_16_13_13_5_96_1_1(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'threads': 96, 'name': 'cusmm_dnt_small_16_13_13_5_96_1_1', 'tile_n': 1, 'tile_m': 1, 'grouping': 16, 'm': 13, 'k': 13, 'n': 5}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_small<13,5,13,1,1,16> <<< ((stack_size + 16 - 1) / 16), 96, shared_size, stream >>>
(param_stack, careful, nruns, 
a_data, b_data, c_data);
return(0);
}


int launch_cusmm_dnt_largeDB_256_16_54_54_54_2_6_27_10(int *param_stack, int stack_size, cudaStream_t stream, int m_max, int n_max, int k_max, double *a_data, double *b_data, double *c_data){
int shared_size = 0;
//{'name': 'cusmm_dnt_largeDB_256_16_54_54_54_2_6_27_10', 'blockdim': 256, 'tile_n': 6, 'tile_m': 2, 'm': 54, 'n': 54, 'w': 10, 'v': 27, 'k': 54, 'grouping': 16}
int careful = (stack_size / 16);
int nruns = stack_size - careful * 16;
cusmm_dnt_largeDB<54,54,54,2,6,10,27,256,16> <<< ((stack_size + 16 - 1) / 16), 256, shared_size, stream >>>
(param_stack, careful, nruns, a_data, b_data, c_data);
return(0);
}


int libcusmm_process_d(int *param_stack, int stack_size,cudaStream_t stream, int m, int n, int k, double *a_data, double *b_data, double *c_data){
int idx = 0;
bool missing = false;

switch(m){
case 5: idx = 0; break;
case 13: idx = 1; break;
case 23: idx = 2; break;
case 26: idx = 3; break;
case 54: idx = 4; break;
default: missing = true;
}

idx *= 5;
switch(n){
case 5: idx += 0; break;
case 13: idx += 1; break;
case 23: idx += 2; break;
case 26: idx += 3; break;
case 54: idx += 4; break;
default: missing = true;
}

idx *= 5;
switch(k){
case 5: idx += 0; break;
case 13: idx += 1; break;
case 23: idx += 2; break;
case 26: idx += 3; break;
case 54: idx += 4; break;
default: missing = true;
}

if(missing) return -1;
switch(idx){
case 0:
// m=5, n=5, k=5
return launch_cusmm_dnt_tiny_16_5_5_5_32_64(param_stack, stack_size, stream, 5, 5, 5, a_data, b_data, c_data);

case 1:
// m=5, n=5, k=13
return launch_cusmm_dnt_tiny_16_13_5_5_32_128(param_stack, stack_size, stream, 5, 5, 13, a_data, b_data, c_data);

case 5:
// m=5, n=13, k=5
return launch_cusmm_dnt_small_16_5_5_13_96_1_1(param_stack, stack_size, stream, 5, 13, 5, a_data, b_data, c_data);

case 6:
// m=5, n=13, k=13
return launch_cusmm_dnt_small_16_13_5_13_96_1_1(param_stack, stack_size, stream, 5, 13, 13, a_data, b_data, c_data);

case 25:
// m=13, n=5, k=5
return launch_cusmm_dnt_small_16_5_13_5_96_1_1(param_stack, stack_size, stream, 13, 5, 5, a_data, b_data, c_data);

case 26:
// m=13, n=5, k=13
return launch_cusmm_dnt_small_16_13_13_5_96_1_1(param_stack, stack_size, stream, 13, 5, 13, a_data, b_data, c_data);

case 30:
// m=13, n=13, k=5
return launch_cusmm_dnt_small_16_5_13_13_96_1_2(param_stack, stack_size, stream, 13, 13, 5, a_data, b_data, c_data);

case 31:
// m=13, n=13, k=13
return launch_cusmm_dnt_medium_16_13_13_13_96_16_96_2_2(param_stack, stack_size, stream, 13, 13, 13, a_data, b_data, c_data);

case 33:
// m=13, n=13, k=26
return launch_cusmm_dnt_medium_16_26_13_13_128_0_128_2_2(param_stack, stack_size, stream, 13, 13, 26, a_data, b_data, c_data);

case 41:
// m=13, n=26, k=13
return launch_cusmm_dnt_medium_16_13_13_26_96_0_96_2_2(param_stack, stack_size, stream, 13, 26, 13, a_data, b_data, c_data);

case 43:
// m=13, n=26, k=26
return launch_cusmm_dnt_largeDB_128_16_26_13_26_1_3_26_4(param_stack, stack_size, stream, 13, 26, 26, a_data, b_data, c_data);

case 62:
// m=23, n=23, k=23
return launch_cusmm_dnt_largeDB_96_16_23_23_23_2_3_13_4(param_stack, stack_size, stream, 23, 23, 23, a_data, b_data, c_data);

case 81:
// m=26, n=13, k=13
return launch_cusmm_dnt_medium_16_13_26_13_128_0_128_2_2(param_stack, stack_size, stream, 26, 13, 13, a_data, b_data, c_data);

case 83:
// m=26, n=13, k=26
return launch_cusmm_dnt_largeDB_96_16_26_26_13_2_2_13_3(param_stack, stack_size, stream, 26, 13, 26, a_data, b_data, c_data);

case 91:
// m=26, n=26, k=13
return launch_cusmm_dnt_medium_16_13_26_26_128_0_128_2_3(param_stack, stack_size, stream, 26, 26, 13, a_data, b_data, c_data);

case 93:
// m=26, n=26, k=26
return launch_cusmm_dnt_largeDB_128_16_26_26_26_2_3_10_4(param_stack, stack_size, stream, 26, 26, 26, a_data, b_data, c_data);

case 124:
// m=54, n=54, k=54
return launch_cusmm_dnt_largeDB_256_16_54_54_54_2_6_27_10(param_stack, stack_size, stream, 54, 54, 54, a_data, b_data, c_data);

}

return -1; // should never happen
}


int libcusmm_transpose_d(int *trs_stack, int offset, int nblks,
double *buffer, int m, int n, cudaStream_t * stream) {
int idx = 0;
bool missing = false;

switch(m){
case 5: idx = 0; break;
case 13: idx = 1; break;
case 23: idx = 2; break;
case 26: idx = 3; break;
case 54: idx = 4; break;
default: missing = true;
}

idx *= 5;
switch(n){
case 5: idx += 0; break;
case 13: idx += 1; break;
case 23: idx += 2; break;
case 26: idx += 3; break;
case 54: idx += 4; break;
default: missing = true;
}

// If there is no kernel for these blocks, we don't need to transpose them.
if(missing) return 0;

switch(idx){
case 0:
// m=5, n=5
transpose_d<5,5> <<<nblks, 128, 0, *stream>>>(trs_stack+offset, nblks, buffer);
break;
case 1:
// m=5, n=13
transpose_d<5,13> <<<nblks, 128, 0, *stream>>>(trs_stack+offset, nblks, buffer);
break;
case 5:
// m=13, n=5
transpose_d<13,5> <<<nblks, 128, 0, *stream>>>(trs_stack+offset, nblks, buffer);
break;
case 6:
// m=13, n=13
transpose_d<13,13> <<<nblks, 128, 0, *stream>>>(trs_stack+offset, nblks, buffer);
break;
case 8:
// m=13, n=26
transpose_d<13,26> <<<nblks, 128, 0, *stream>>>(trs_stack+offset, nblks, buffer);
break;
case 12:
// m=23, n=23
transpose_d<23,23> <<<nblks, 128, 0, *stream>>>(trs_stack+offset, nblks, buffer);
break;
case 16:
// m=26, n=13
transpose_d<26,13> <<<nblks, 128, 0, *stream>>>(trs_stack+offset, nblks, buffer);
break;
case 18:
// m=26, n=26
transpose_d<26,26> <<<nblks, 128, 0, *stream>>>(trs_stack+offset, nblks, buffer);
break;
case 24:
// m=54, n=54
transpose_d<54,54> <<<nblks, 128, 0, *stream>>>(trs_stack+offset, nblks, buffer);
break;
// If there is no kernel for these blocks, we don't need to transpose them.
default: return(0);
}

return(cudaGetLastError());
}
//EOF
