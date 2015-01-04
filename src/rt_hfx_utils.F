!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_col_block_sizes, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix_set, &
       cp_dbcsr_distribution, cp_dbcsr_init, cp_dbcsr_p_type, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_set, dbcsr_type_antisymmetric
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             set_ks_env
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_set,&
                                             qs_rho_type
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
!> \param qs_env ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE rtp_hfx_rebuild(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rtp_hfx_rebuild', &
      routineP = moduleN//':'//routineN

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_im, &
                                                matrix_ks_im, matrix_s, &
                                                matrix_s_aux, rho_ao_im, &
                                                rho_aux_ao_im
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_aux, sab_orb
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho, rho_aux
    TYPE(rt_prop_type), POINTER              :: rtp

    NULLIFY(ks_env,matrix_s,matrix_ks_im,dft_control,sab_orb,rtp,rho,rho_ao_im)
    NULLIFY(sab_aux, rho_aux, rho_aux_ao_im, matrix_ks_aux_im, matrix_s_aux)

    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    dft_control=dft_control,&
                    rtp=rtp,&
                    matrix_s=matrix_s,&
                    rho=rho,&
                    matrix_ks_im=matrix_ks_im,&
                    sab_orb=sab_orb,&
                    error=error)

    CALL qs_rho_get(rho, rho_ao_im=rho_ao_im, error=error)
    CALL rebuild_matrices(rho_ao_im, matrix_ks_im, sab_orb, matrix_s,&
                          dft_control%nspins, error)
    CALL set_ks_env(ks_env, matrix_ks_im=matrix_ks_im, error=error)
    CALL qs_rho_set(rho, rho_ao_im=rho_ao_im, error=error)

    IF(dft_control%do_admm)THEN
       CALL get_qs_env(qs_env,&
                       matrix_s_aux_fit=matrix_s_aux,&
                       sab_aux_fit=sab_aux,&
                       rho_aux_fit=rho_aux,&
                       matrix_ks_aux_fit_im=matrix_ks_aux_im,&
                       error=error)
       CALL qs_rho_get(rho_aux, rho_ao_im=rho_aux_ao_im, error=error)
       CALL rebuild_matrices(rho_aux_ao_im, matrix_ks_aux_im, sab_aux, matrix_s_aux,&
                             dft_control%nspins,error)
       CALL set_ks_env(ks_env, matrix_ks_aux_fit_im=matrix_ks_aux_im, error=error)
       CALL qs_rho_set(rho_aux, rho_ao_im=rho_aux_ao_im, error=error)
    END IF

  END SUBROUTINE rtp_hfx_rebuild

! *****************************************************************************
!> \brief does the actual rebuilding of P and KS (imaginary) in case S changed
!> \param matrix_p ...
!> \param matrix_ks ...
!> \param sab_orb ...
!> \param matrix_s ...
!> \param nspins ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************

  SUBROUTINE rebuild_matrices(matrix_p,matrix_ks,sab_orb,matrix_s,nspins,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_ks
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    INTEGER, INTENT(in)                      :: nspins
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
