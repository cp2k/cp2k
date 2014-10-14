!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE qs_tddfpt_types
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm,&
                                             fm_pool_give_back_fm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE termination,                     ONLY: stop_memory
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

! *****************************************************************************
  TYPE tddfpt_env_type
     REAL(KIND=dp), DIMENSION(:), POINTER               :: evals  ! eigenvalues
     TYPE(cp_fm_p_type), DIMENSION(:,:), &
          POINTER                                       :: evecs  ! eigenvectors
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: invS   ! the inverse of the metric
     TYPE(cp_fm_pool_p_type), DIMENSION(:), &
          POINTER                                       :: ao_mo_fm_pools
  END TYPE tddfpt_env_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt_types'

  PRIVATE

  PUBLIC :: tddfpt_env_type

  PUBLIC :: tddfpt_env_allocate, &
            tddfpt_env_deallocate

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param t_env ...
!> \param p_env ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE tddfpt_env_allocate(t_env, p_env, qs_env, error)

    TYPE(tddfpt_env_type), INTENT(inout)     :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_env_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n_ev, n_spins, spin, stat
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(dft_control_type), POINTER          :: dft_control

    NULLIFY(dft_control)

    CALL get_qs_env(qs_env, dft_control=dft_control, error=error)
    n_ev    = dft_control%tddfpt_control%n_ev
    n_spins = dft_control%nspins
    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=t_env%ao_mo_fm_pools,error=error)

    ALLOCATE (t_env%evals(n_ev),STAT=stat)
    IF (stat /=0 ) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "t_env%evals",dp_size*n_ev)

    ALLOCATE (t_env%evecs(n_ev,n_spins),STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "t_env%evecs",dp_size*n_ev*n_spins)
    DO spin=1,n_spins
       DO i=1, n_ev
          CALL fm_pool_create_fm(t_env%ao_mo_fm_pools(spin)%pool, &
                                 t_env%evecs(i,spin)%matrix,error=error)
       END DO
    END DO

    !-------------------------------------------------!
    ! allocate memory for the inverse of the S matrix !
    !-------------------------------------------------!
    ALLOCATE (t_env%invS(n_spins),STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "t_env%invS",int_size*n_spins)
    DO spin=1, n_spins
       NULLIFY(fm_struct)
       CALL cp_fm_struct_create(fm_struct, qs_env%para_env, qs_env%blacs_env, &
                                p_env%n_ao(spin), p_env%n_ao(spin),error=error)
       CALL cp_fm_create(t_env%invS(spin)%matrix, fm_struct, routineP//"invS",error=error)
       CALL cp_fm_struct_release(fm_struct,error=error)
    END DO

  END SUBROUTINE tddfpt_env_allocate

! *****************************************************************************
!> \brief ...
!> \param t_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE tddfpt_env_deallocate(t_env,error)

    TYPE(tddfpt_env_type), INTENT(inout)     :: t_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_env_deallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, spin, stat

    DO spin=1, SIZE(t_env%evecs, 2)
       DO i=1, SIZE(t_env%evecs, 1)
          CALL fm_pool_give_back_fm(t_env%ao_mo_fm_pools(spin)%pool, &
                                    t_env%evecs(i,spin)%matrix,error=error)
       END DO
    END DO

    DO spin=1, SIZE(t_env%invS)
       IF (ASSOCIATED(t_env%invS(spin)%matrix)) &
            CALL cp_fm_release(t_env%invS(spin)%matrix,error=error)
    END DO
    DEALLOCATE (t_env%invS,t_env%evecs,t_env%evals,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "t_env%invS,t_env%evecs,t_env%evals")

    NULLIFY (t_env%invS)

  END SUBROUTINE tddfpt_env_deallocate

END MODULE qs_tddfpt_types
