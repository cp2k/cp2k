!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_rho_methods [1.0] *
!!
!!   NAME
!!     qs_rho_methods
!!
!!   FUNCTION
!!     methods of the rho structure (defined in qs_rho_types)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_rho_methods
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_give_back_coeffs,&
                                             pw_pools_init_coeffs
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_rho_methods'
  INTEGER, PRIVATE, SAVE :: last_rho_id_nr=0

  PUBLIC :: qs_rho_update_rho, qs_rho_rebuild

!***
!****************************************************************************
CONTAINS

!!****f* qs_rho_methods/qs_rho_rebuild [1.0] *
!!
!!   NAME
!!     qs_rho_rebuild
!!
!!   SYNOPSIS
!!     Subroutine qs_rho_rebuild(rho, qs_env, &
!!         rebuild_ao, rebuild_grids, error)
!!       Type(qs_rho_type), Pointer:: rho
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Logical, Intent (IN), Optional::  rebuild_ao,&
!!         rebuild_grids
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_rho_rebuild
!!
!!   FUNCTION
!!     rebuilds rho (if necessary allocating and initializing it)
!!
!!   NOTES
!!     needs updated  pw pools, s, s_mstruct and h in qs_env.
!!     The use of p to keep the structure of h (needed for the forces)
!!     is ugly and should be removed.
!!     Change so that it does not allocate a subcomponent if it is not 
!!     associated and not requested?
!!
!!   INPUTS
!!     - rho: the rho type to rebuild (defaults to qs_env%rho)
!!     - qs_env: the environment to which rho belongs
!!     - rebuild_ao: if it is necessary to rebuild rho_ao. Defaults to true.
!!     - rebuild_grids: if it in necessary to rebuild rho_r and rho_g.
!!       Defaults to false.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created replacing qs_rho_create and qs_env_rebuild_rho[fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE qs_rho_rebuild(rho, qs_env, rebuild_ao, rebuild_grids, error)
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: rebuild_ao, rebuild_grids
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_rebuild', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nspins, stat
    LOGICAL                                  :: failure, my_rebuild_ao, &
                                                my_rebuild_grids
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: s

  NULLIFY(pw_env,auxbas_pw_pool,s,dft_control)
  CALL timeset(routineN,'I',' ',handle)
  failure=.FALSE.
  my_rebuild_ao=.TRUE.
  my_rebuild_grids=.TRUE.
  IF (PRESENT(rebuild_ao)) my_rebuild_ao=rebuild_ao
  IF (PRESENT(rebuild_grids)) my_rebuild_grids=rebuild_grids

  CALL get_qs_env(qs_env, pw_env=pw_env,dft_control=dft_control,&
       s=s,error=error)
  CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
  nspins=dft_control%nspins

  IF (.NOT.ASSOCIATED(rho)) THEN
     my_rebuild_grids=.TRUE.
     my_rebuild_ao=.TRUE.
     ALLOCATE(rho,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (.NOT.failure) THEN
        rho%ref_count=1
        rho%rho_g_valid=.FALSE.
        rho%rho_r_valid=.FALSE.
        rho%tau_g_valid=.FALSE.
        rho%tau_r_valid=.FALSE.
        last_rho_id_nr=last_rho_id_nr+1
        rho%id_nr=last_rho_id_nr
        rho%rebuild_each=5
        NULLIFY(rho%rho_r, rho%rho_g, rho%rho_ao, rho%tot_rho_r, rho%tot_rho_g, rho%tau_r, rho%tau_g)
        ALLOCATE(rho%tot_rho_r(nspins),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        IF (.not.failure) rho%tot_rho_r=0.0_dp
     END IF
  END IF

  ! rho_ao
  IF (my_rebuild_ao.OR..NOT.ASSOCIATED(rho%rho_ao)) THEN
     IF (.NOT.failure) THEN
        IF (ASSOCIATED(rho%rho_ao)) THEN
           CALL deallocate_matrix_set(rho%rho_ao)
        END IF
        CALL allocate_matrix_set(rho%rho_ao,nspins)
        DO i=1,nspins
           NULLIFY (rho%rho_ao(i)%matrix)
           CALL replicate_matrix_structure(s(1)%matrix,&
                rho%rho_ao(i)%matrix,&
                "DENSITY_MATRIX_FOR_SPIN_"//TRIM(ADJUSTL(cp_to_string(i))))
        END DO
     END IF
  END IF

  ! rho_r
  IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%rho_r)) THEN
     IF (.NOT.failure) THEN
        IF (ASSOCIATED(rho%rho_r)) THEN
           DO i=1,SIZE(rho%rho_r)
              CALL coeff_deallocate(rho%rho_r(i))
           END DO
           DEALLOCATE(rho%rho_r,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        ALLOCATE(rho%rho_r(nspins),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF (.NOT.failure) THEN
        DO i=1,nspins
           CALL pw_pool_init_coeff(auxbas_pw_pool,rho%rho_r(i),&
                use_data=REALDATA3D,in_space=REALSPACE,error=error)
        END DO
     END IF
  END IF

  ! rho_g
  IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%rho_g)) THEN
     IF (.NOT.failure) THEN
        IF (ASSOCIATED(rho%rho_g)) THEN
           DO i=1,SIZE(rho%rho_g)
              CALL coeff_deallocate(rho%rho_g(i))
           END DO
           DEALLOCATE(rho%rho_g,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        ALLOCATE(rho%rho_g(nspins),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF (.NOT.failure) THEN
        DO i=1,nspins
           CALL pw_pool_init_coeff(auxbas_pw_pool,rho%rho_g(i),&
                use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
        END DO
     END IF
  END IF

  ! allocate tau_r and tau_g if use_kinetic_energy_density 
  IF (dft_control%use_kinetic_energy_density) THEN
     ! tau_r
     IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%tau_r)) THEN
        IF (.NOT.failure) THEN
           IF (ASSOCIATED(rho%tau_r)) THEN
              DO i=1,SIZE(rho%tau_r)
                 CALL coeff_deallocate(rho%tau_r(i))
              END DO
              DEALLOCATE(rho%tau_r,stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
           ALLOCATE(rho%tau_r(nspins),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF (.NOT.failure) THEN
           DO i=1,nspins
              CALL pw_pool_init_coeff(auxbas_pw_pool,rho%tau_r(i),&
                   use_data=REALDATA3D,in_space=REALSPACE,error=error)
           END DO
        END IF
     END IF

     ! tau_g
     IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%tau_g)) THEN
        IF (.NOT.failure) THEN
           IF (ASSOCIATED(rho%tau_g)) THEN
              DO i=1,SIZE(rho%tau_g)
                 CALL coeff_deallocate(rho%tau_g(i))
              END DO
              DEALLOCATE(rho%tau_g,stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
           ALLOCATE(rho%tau_g(nspins),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF (.NOT.failure) THEN
           DO i=1,nspins
              CALL pw_pool_init_coeff(auxbas_pw_pool,rho%tau_g(i),&
                   use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
           END DO
        END IF
     END IF
  END IF  ! use_kinetic_energy_density
  CALL timestop(0.0_dp,handle) 
END SUBROUTINE qs_rho_rebuild
!***************************************************************************

!!****f* qs_rho_methods/qs_rho_update_rho [1.0] *
!!
!!   NAME
!!     qs_rho_update_rho
!!
!!   SYNOPSIS
!!     Subroutine qs_rho_update_rho(rho_struct,  qs_env, &
!!          error)
!!       Type(qs_rho_type), Pointer:: rho_struct
!!       Type(qs_environment_type), Intent (IN):: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_rho_update_rho
!!
!!   FUNCTION
!!     updates rho_r and rho_g to the rho%rho_ao.
!!     if use_kinetic_energy_density also computes tau_r and tau_g
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho_struct: the rho structure that should be updated
!!     - qs_env: the qs_env rho_struct refers to
!!       the integrated charge in r space
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_rho_update_rho(rho_struct, qs_env, error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(qs_environment_type), INTENT(in)    :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_update_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, nspins, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dum
    TYPE(coeff_type), DIMENSION(:), POINTER  :: mgrid_gspace, mgrid_rspace
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

  failure=.FALSE.
  NULLIFY(dft_control, mgrid_rspace, mgrid_gspace, pw_env, auxbas_pw_pool,&
       pw_pools, para_env )

  CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN

     CALL get_qs_env(qs_env, dft_control=dft_control, &
          pw_env=pw_env, para_env=para_env)
     nspins=dft_control%nspins
     ! alloc multigrids
     CALL pw_env_get(pw_env, pw_pools=pw_pools, auxbas_pw_pool=auxbas_pw_pool,&
          error=error)
     ALLOCATE(mgrid_gspace(SIZE(pw_pools)), mgrid_rspace(SIZE(pw_pools)),&
          stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     ! get tmp grids
     CALL pw_pools_init_coeffs(pw_pools, mgrid_rspace,&
          use_data=REALDATA3D, in_space=REALSPACE, error=error)
     CALL pw_pools_init_coeffs(pw_pools, mgrid_gspace,&
          use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, error=error)
    
     DO ispin=1,nspins
        CALL calculate_rho_elec(p=rho_struct%rho_ao(ispin)%matrix,&
             rho=rho_struct%rho_r(ispin),&
             rho_gspace=rho_struct%rho_g(ispin),&
             mgrid_rspace=mgrid_rspace,&
             mgrid_gspace=mgrid_gspace,&
             gridlevel_info=pw_env%gridlevel_info,&
             cube_info=pw_env%cube_info,&
             l_info=pw_env%l_info,&
             total_rho=rho_struct%tot_rho_r(ispin),&
             qs_env=qs_env,error=error)
     END DO
     rho_struct%rho_r_valid=.TRUE.
     rho_struct%rho_g_valid=.TRUE.

     ! if needed compute also the kinetic energy density
     IF (dft_control%use_kinetic_energy_density) THEN
         DO ispin=1,nspins
            CALL calculate_rho_elec(p=rho_struct%rho_ao(ispin)%matrix,&
                 rho=rho_struct%tau_r(ispin),&
                 rho_gspace=rho_struct%tau_g(ispin),&
                 mgrid_rspace=mgrid_rspace,&
                 mgrid_gspace=mgrid_gspace,&
                 gridlevel_info=pw_env%gridlevel_info,&
                 cube_info=pw_env%cube_info,&
                 l_info=pw_env%l_info,&
                 total_rho=dum, &  ! presumably not meaningful
                 qs_env=qs_env, compute_tau=.TRUE., error=error)
         END DO
         rho_struct%tau_r_valid=.TRUE.
         rho_struct%tau_g_valid=.TRUE.
     ENDIF

     ! put back grids
     CALL pw_pools_give_back_coeffs(pw_pools, mgrid_rspace, error=error)
     CALL pw_pools_give_back_coeffs(pw_pools, mgrid_gspace, error=error)
     DEALLOCATE(mgrid_gspace, mgrid_rspace,stat=stat)
  END IF
END SUBROUTINE qs_rho_update_rho
!***************************************************************************

END MODULE qs_rho_methods
