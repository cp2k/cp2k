!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Debug routine for long-range energy (debug value of EWALD vs VALUE KS)
!> \author Teodoro Laino [tlaino] - 04.2009
! *****************************************************************************
INTERFACE

   SUBROUTINE dbg_energy_coulomb_lr(energy, ks_matrix, nspins, qs_env, matrix_p,&
        calculate_forces, para_env, error)
     USE cp_para_types,                   ONLY: cp_para_env_type
     USE f77_blas
     USE qs_energy_types,                 ONLY: init_qs_energy,&
                                                qs_energy_type
     USE qs_environment_types,            ONLY: qs_environment_type
     USE sparse_matrix_types,             ONLY: cp_sm_sm_trace,&
                                                real_matrix_p_type,&
                                                set_matrix
#include "cp_common_uses.h"
   
       TYPE(qs_energy_type), POINTER            :: energy
       TYPE(real_matrix_p_type), DIMENSION(:), &
         POINTER                                :: ks_matrix
       INTEGER, INTENT(IN)                      :: nspins
       TYPE(qs_environment_type), POINTER       :: qs_env
       TYPE(real_matrix_p_type), DIMENSION(:), &
         POINTER                                :: matrix_p
       LOGICAL, INTENT(IN)                      :: calculate_forces
       TYPE(cp_para_env_type), POINTER          :: para_env
       TYPE(cp_error_type), INTENT(inout)       :: error
   
   END SUBROUTINE dbg_energy_coulomb_lr


END INTERFACE
