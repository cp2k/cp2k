!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utility functions that are needed for RTP/EMD in combination with 
!>        HF or hybrid functionals (needs to deal with imaginary KS and P
!> \par History
!>      2014 created [fschiff]
!> \author Florina Schiffmann
! *****************************************************************************
MODULE rt_hfx_utils
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_col_block_sizes,&
                                             cp_dbcsr_create,&
                                             cp_dbcsr_distribution,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_row_block_sizes,&
                                             cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_deallocate_matrix_set
  USE dbcsr_api,                       ONLY: dbcsr_type_antisymmetric
  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE rt_propagation_types,            ONLY: rt_prop_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_hfx_utils'

  PUBLIC :: rtp_hfx_rebuild

!***
CONTAINS
! *****************************************************************************
!> \brief rebuilds the structures of P and KS (imaginary) in case S changed
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE rtp_hfx_rebuild(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rtp_hfx_rebuild', &
      routineP = moduleN//':'//routineN

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,rtp=rtp,error=error)
    IF(dft_control%do_admm)THEN
       CALL get_qs_env(qs_env=qs_env,matrix_s_aux_fit=matrix_s,&
                       sab_aux_fit=sab_orb,error=error)
       CALL rebuild_matrices(rtp%rtp_hfx%matrix_p_admm_im,rtp%rtp_hfx%matrix_ks_admm_im,&
                             sab_orb,matrix_s,dft_control%nspins,error)
    END IF
    CALL get_qs_env(qs_env=qs_env,matrix_s=matrix_s,&
                    sab_orb=sab_orb,error=error)

    CALL rebuild_matrices(rtp%rtp_hfx%matrix_p_im,rtp%rtp_hfx%matrix_ks_im,&
                          sab_orb,matrix_s,dft_control%nspins,error)

  END SUBROUTINE rtp_hfx_rebuild

! *****************************************************************************
!> \brief does the actual rebuilding of P and KS (imaginary) in case S changed
!> \author Florian Schiffmann
! *****************************************************************************

  SUBROUTINE rebuild_matrices(matrix_p,matrix_ks,sab_orb,matrix_s,nspins,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_ks
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    INTEGER                                  :: nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rebuild_matrices', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: headline
    INTEGER                                  :: i

    IF (ASSOCIATED(matrix_p)) THEN
       CALL cp_dbcsr_deallocate_matrix_set(matrix_p,error=error)
    END IF
    ! Create a new density matrix set
    CALL cp_dbcsr_allocate_matrix_set(matrix_p,nspins,error=error)
    DO i=1,nspins
       headline = "Imaginary density matrix"
       ALLOCATE(matrix_p(i)%matrix)
       CALL cp_dbcsr_init(matrix_p(i)%matrix, error=error)
       CALL cp_dbcsr_create(matrix=matrix_p(i)%matrix,&
            name=TRIM(headline),&
            dist=cp_dbcsr_distribution(matrix_s(1)%matrix), matrix_type=dbcsr_type_antisymmetric,&
            row_blk_size=cp_dbcsr_row_block_sizes(matrix_s(1)%matrix),&
            col_blk_size=cp_dbcsr_col_block_sizes(matrix_s(1)%matrix),&
            nblks=0, nze=0, error=error)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_p(i)%matrix,sab_orb,error=error)
       CALL cp_dbcsr_set(matrix_p(i)%matrix,0.0_dp,error=error)
    END DO

    IF (ASSOCIATED(matrix_ks)) THEN
       CALL cp_dbcsr_deallocate_matrix_set(matrix_ks,error=error)
    END IF
    ! Create a new density matrix set
    CALL cp_dbcsr_allocate_matrix_set(matrix_ks,nspins,error=error)
    DO i=1,nspins
       headline = "Imaginary Kohn-Sham matrix"
       ALLOCATE(matrix_ks(i)%matrix)
       CALL cp_dbcsr_init(matrix_ks(i)%matrix, error=error)
       CALL cp_dbcsr_create(matrix=matrix_ks(i)%matrix,&
            name=TRIM(headline),&
            dist=cp_dbcsr_distribution(matrix_s(1)%matrix), matrix_type=dbcsr_type_antisymmetric,&
            row_blk_size=cp_dbcsr_row_block_sizes(matrix_s(1)%matrix),&
            col_blk_size=cp_dbcsr_col_block_sizes(matrix_s(1)%matrix),&
            nblks=0, nze=0, error=error)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks(i)%matrix,sab_orb,error=error)
       CALL cp_dbcsr_set(matrix_ks(i)%matrix,0.0_dp,error=error)
    END DO



   END SUBROUTINE rebuild_matrices

END MODULE rt_hfx_utils
