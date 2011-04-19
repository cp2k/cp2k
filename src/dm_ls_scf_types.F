!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types needed for a linear scaling quickstep SCF run based on the density
!>        matrix
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE dm_ls_scf_types
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: ls_cluster_atomic,&
                                             ls_cluster_molecular
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf_types'

  PUBLIC :: ls_scf_env_type, ls_mstruct_type, ls_cluster_atomic, ls_cluster_molecular

  TYPE ls_mstruct_type
     INTEGER :: cluster_type       
     LOGICAL :: single_precision
     INTEGER, DIMENSION(:), ALLOCATABLE :: atom_to_molecule
  END TYPE

  TYPE ls_scf_env_type
     INTEGER               :: nspins, natoms
     INTEGER               :: nelectron_total
     INTEGER, DIMENSION(2) :: nelectron_spin
     REAL(KIND=dp), DIMENSION(2) ::  mu_spin


     LOGICAL :: has_unit_metric

     INTEGER :: s_preconditioner_type
     INTEGER :: s_inversion_type
     LOGICAL               :: needs_s_inv, has_s_preconditioner, fixed_mu
     REAL(KIND=dp)         :: eps_filter, eps_scf

     REAL(KIND=dp)         :: mixing_fraction
     INTEGER               :: max_scf

     TYPE(cp_dbcsr_type)   :: matrix_s_inv
     TYPE(cp_dbcsr_type)   :: matrix_s
     TYPE(cp_dbcsr_type)   :: matrix_bs_sqrt, matrix_bs_sqrt_inv
     TYPE(cp_dbcsr_type)   :: matrix_s_sqrt, matrix_s_sqrt_inv
     TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_ks
     TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_p

     LOGICAL  :: report_all_sparsities, perform_mu_scan, use_s_sqrt

     TYPE(ls_mstruct_type) :: ls_mstruct

     TYPE(cp_para_env_type), POINTER :: para_env
  END TYPE ls_scf_env_type

END MODULE dm_ls_scf_types


