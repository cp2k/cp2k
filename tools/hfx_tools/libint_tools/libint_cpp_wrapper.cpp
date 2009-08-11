#include <libint/libint.h>
#include <libderiv/libderiv.h>
extern "C" {

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Wrappers for libint library
//
//  DESCRIPTION:
//      wrap_to_build_eri_        : Used by Fortran Compilers partially supporting ISO_C_BINDING (NAG, INTEL)        
//      wrap_to_build_eri_        : Used by Fortran Compilers partially supporting ISO_C_BINDING (NAG, INTEL)
//      wrapper_init_lib_         : Used by Fortran Compilers not supporting ISO_C_BINDING (PGI)
//      wrapper_init_deriv_       : Used by Fortran Compilers not supporting ISO_C_BINDING (PGI)
//      wrapper_build_eri_        : Used by Fortran Compilers not supporting ISO_C_BINDING (PGI)
//      wrapper_build_eri_cray_   : Used by Fortran Compilers not supporting ISO_C_BINDING (PGI) but CRAY POINTERS
//      wrapper_build_deriv1_eri_ : Used by Fortran Compilers not supporting ISO_C_BINDING (PGI)
//      wrapper_free_libint_      : Used by Fortran Compilers not supporting ISO_C_BINDING (PGI)    
//      wrapper_
    

  extern REALTYPE* wrap_to_build_eri(int *n_d, int *n_c, int *n_b, int *n_a, Libint_t *lib)
  {
      return build_eri[*n_d][*n_c][*n_b][*n_a] ( lib, 1 );
  }
  extern void wrap_to_build_deriv1_eri(int *n_d, int *n_c, int *n_b, int *n_a, Libderiv_t *deriv)
  {
     build_deriv1_eri[*n_d][*n_c][*n_b][*n_a]( deriv, 1 );
  }
  
  void wrapper_init_lib_(Libint_t *lib, int* max_am,int* max_prim, int* lib_storage)
  {
    init_libint_base();
    *lib_storage = init_libint(lib, *max_am, *max_prim);
    return;
  }

  void wrapper_init_deriv_(Libderiv_t* deriv, int* max_am, int* max_prim, int* max_classes, int* lib_deriv_storage)
  {
    init_libderiv_base(); 
    *lib_deriv_storage = init_libderiv1(deriv, *max_am, *max_prim, *max_classes);
   }
  
  void wrapper_build_eri_(int *n_d, int *n_c, int *n_b, int *n_a, Libint_t *lib, int* size, REALTYPE* array, prim_data* prim)
  {
    REALTYPE* pRetVal;
    int i;
    lib->PrimQuartet = prim; 
    pRetVal = build_eri[*n_d][*n_c][*n_b][*n_a] ( lib, 1 );

    for(i=0;i<*size;i++)
    {
      array[i] = *pRetVal;
      pRetVal++;
    }
  }
  void wrapper_build_eri_cray_(int *n_d, int *n_c, int *n_b, int *n_a, Libint_t *lib, prim_data* prim, REALTYPE** ptr)
  {
    lib->PrimQuartet = prim;
    *ptr = build_eri[*n_d][*n_c][*n_b][*n_a] ( lib, 1 );
  }

  void wrapper_build_deriv1_eri_(int *n_d, int *n_c, int *n_b, int *n_a, Libderiv_t *deriv, int* size, REALTYPE* array, prim_data* prim)
  {
    REALTYPE* pDeriv;
    deriv->PrimQuartet = prim;
    build_deriv1_eri[*n_d][*n_c][*n_b][*n_a]( deriv, 1 );
    for(int i=0;i<12;i++)
    {
      if(i!=3 && i!= 4 && i!=5) 
      {
       pDeriv = deriv->ABCD[i];
       for(int k=0;k<*size;k++)
       {
         array[i*(*size)+k] = *pDeriv;
         pDeriv++;
       } 
      }
    }
  }
  void wrapper_free_libint_(Libint_t* lib)
  {
    free_libint(lib);
  }
  
  void wrapper_free_libderiv_(Libderiv_t* deriv)
  {
    free_libderiv(deriv);
  }

}

