/****************************************************************************
*  CP2K: A general program to perform molecular dynamics simulations        *
*  Copyright (C) 2000 - 2013 the CP2K developers group                      *
*****************************************************************************/

#ifdef __cplusplus
 extern "C" {
#endif

int libsmm_acc_process (int *param_stack, int stack_size,
    int nparams, int datatype, void *a_data, void *b_data, void *c_data,
    int m_max, int n_max, int k_max, int def_mnk, void* stream);

int libsmm_acc_transpose (int *trs_stack, int offset, int nblks,
    void *buffer,int datatype, int m, int n, void* stream);

#ifdef __cplusplus
 }
#endif

//EOF
