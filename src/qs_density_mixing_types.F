!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module that contains the definitions of the scf types
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE qs_density_mixing_types

  USE ao_util,                         ONLY: exp_radius
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE qs_rho_atom_types,               ONLY: rho_atom_coeff
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_density_mixing_types'

  INTEGER, PARAMETER, PUBLIC :: no_mixing_nr=0, direct_mixing_nr=1,&
                                 gspace_mixing_nr=2, pulay_mixing_nr=3,&
                                 broyden_mixing_nr=4, &
                                 broyden_mixing_new_nr=5,&
                                 multisecant_mixing_nr=6
  PUBLIC :: cp_1d_z_p_type, mixing_storage_create, mixing_storage_type, mixing_storage_release

  TYPE cp_1d_z_p_type
     COMPLEX(dp), DIMENSION(:), POINTER :: cc
  END TYPE cp_1d_z_p_type


  TYPE mixing_storage_type
     INTEGER :: ref_count, ig_max, ncall, nbuffer, n_simple_mix, nskip_mixing, p_metric_method
     INTEGER, POINTER, DIMENSION(:) :: ig_global_index
     LOGICAL, POINTER, DIMENSION(:) :: paw
     CHARACTER (len=15) :: iter_method
     REAL(dp) :: alpha, bconst, beta, broy_w0, max_g2, max_gvec_exp, pulay_alpha, pulay_beta, r_step,&
                 reg_par, sigma_max, wc, wmax
     REAL(dp), DIMENSION(:), POINTER :: p_metric
     REAL(KIND=dp), DIMENSION(:), POINTER :: kerker_factor
     REAL(KIND=dp), DIMENSION(:), POINTER :: special_metric
     REAL(dp), DIMENSION(:,:), POINTER :: weight
     REAL(KIND = dp) , DIMENSION(:,:), POINTER :: norm_res_buffer,pulay_matrix
     REAL(dp), DIMENSION(:,:,:), POINTER :: fmat, gmat, smat
     TYPE(cp_1d_z_p_type), DIMENSION(:), POINTER :: last_res, rhoin, rhoin_old
     TYPE(cp_1d_z_p_type), DIMENSION(:,:), POINTER :: delta_res, u_vec, z_vec
     TYPE(cp_1d_z_p_type), DIMENSION(:,:), POINTER ::drho_buffer, rhoin_buffer, &
      res_buffer
     TYPE(rho_atom_coeff), DIMENSION(:,:), POINTER :: cpc_h_lastres, cpc_s_lastres
     TYPE(rho_atom_coeff), DIMENSION(:,:), POINTER :: cpc_h_in, cpc_s_in
     TYPE(rho_atom_coeff), DIMENSION(:,:), POINTER :: cpc_h_old, cpc_s_old
     TYPE(rho_atom_coeff), DIMENSION(:,:,:), POINTER :: dcpc_h_in, dcpc_s_in
  END TYPE mixing_storage_type

CONTAINS

! *****************************************************************************
!> \brief creates a mixing_storage
!> \par History
!>      05.2009 created [MI]
!> \author [MI]
! *****************************************************************************
  SUBROUTINE mixing_storage_create(mixing_store,  mixing_section, mixing_method, ecut, error)
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(section_vals_type), POINTER         :: mixing_section
    INTEGER, INTENT(IN)                      :: mixing_method
    REAL(dp), INTENT(IN)                     :: ecut
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mixing_storage_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    REAL(dp)                                 :: alpha, eps, gcut

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(mixing_store),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
      ALLOCATE(mixing_store, stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      IF (.NOT.failure) THEN

        mixing_store%ref_count=1
        mixing_store%nbuffer=0
        mixing_store%n_simple_mix=0
        mixing_store%ncall=0
        mixing_store%alpha=1.0_dp
        mixing_store%pulay_beta=1.0_dp
        mixing_store%beta=1.0_dp
        mixing_store%iter_method = "NoMix"
        mixing_store%max_g2 = 2._dp*ecut

        NULLIFY(mixing_store%p_metric)
        NULLIFY(mixing_store%kerker_factor)
        NULLIFY(mixing_store%special_metric)
        NULLIFY(mixing_store%pulay_matrix)
        NULLIFY(mixing_store%weight)
        NULLIFY(mixing_store%fmat)
        NULLIFY(mixing_store%gmat)
        NULLIFY(mixing_store%smat)
        NULLIFY(mixing_store%last_res)
        NULLIFY(mixing_store%rhoin)
        NULLIFY(mixing_store%rhoin_old)
        NULLIFY(mixing_store%delta_res)
        NULLIFY(mixing_store%u_vec)
        NULLIFY(mixing_store%z_vec)
        NULLIFY(mixing_store%drho_buffer)
        NULLIFY(mixing_store%rhoin_buffer)
        NULLIFY(mixing_store%res_buffer)
        NULLIFY(mixing_store%norm_res_buffer)
        NULLIFY(mixing_store%ig_global_index)
        NULLIFY(mixing_store%paw)
        NULLIFY(mixing_store%cpc_h_in)
        NULLIFY(mixing_store%cpc_s_in)
        NULLIFY(mixing_store%cpc_h_old)
        NULLIFY(mixing_store%cpc_s_old)
        NULLIFY(mixing_store%dcpc_h_in)
        NULLIFY(mixing_store%dcpc_s_in)
        NULLIFY(mixing_store%cpc_h_lastres)
        NULLIFY(mixing_store%cpc_s_lastres)

        CALL section_vals_val_get(mixing_section,"ALPHA",r_val=mixing_store%alpha,error=error)
        CALL section_vals_val_get(mixing_section,"BETA",r_val=mixing_store%beta,error=error)
        CALL section_vals_val_get(mixing_section,"N_SIMPLE_MIX",i_val=mixing_store%n_simple_mix,error=error)
        CALL section_vals_val_get(mixing_section,"NBUFFER",i_val=mixing_store%nbuffer,error=error)
        CALL section_vals_val_get(mixing_section,"NSKIP",i_val=mixing_store%nskip_mixing,error=error)
        CALL section_vals_val_get(mixing_section,"MAX_GVEC_EXP",r_val=mixing_store%max_gvec_exp,error=error)

        IF(mixing_store%max_gvec_exp > 0._dp) THEN
          alpha = 0.25_dp/mixing_store%max_gvec_exp
          eps = 1.e-4_dp
          gcut = exp_radius(3,alpha,eps,1.0_dp)
          mixing_store%max_g2 = gcut*gcut
        END IF

        SELECT CASE(mixing_method)
        CASE(gspace_mixing_nr)
           mixing_store%nbuffer = 1
        CASE(pulay_mixing_nr)
           CALL section_vals_val_get(mixing_section,"PULAY_ALPHA",r_val=mixing_store%pulay_alpha,error=error)
           CALL section_vals_val_get(mixing_section,"PULAY_BETA",r_val=mixing_store%pulay_beta,error=error)
        CASE(broyden_mixing_nr)
           CALL section_vals_val_get(mixing_section,"BROY_W0",r_val=mixing_store%broy_w0,error=error)
           mixing_store%bconst = 20.0_dp
        CASE(broyden_mixing_new_nr)
           CALL section_vals_val_get(mixing_section,"BROY_WREF",r_val=mixing_store%wc,error=error)
           CALL section_vals_val_get(mixing_section,"BROY_WMAX",r_val=mixing_store%wmax,error=error)
           mixing_store%bconst = 20.0_dp
           mixing_store%p_metric_method=1
        CASE(multisecant_mixing_nr)
           CALL section_vals_val_get(mixing_section,"REGULARIZATION",r_val=mixing_store%reg_par,error=error)
           CALL section_vals_val_get(mixing_section,"MAX_STEP",r_val=mixing_store%sigma_max,error=error)
           CALL section_vals_val_get(mixing_section,"R_FACTOR",r_val=mixing_store%r_step,error=error)
        END SELECT
      END IF
    END IF

 END SUBROUTINE mixing_storage_create

! *****************************************************************************
!> \brief releases a mixing_storage
!> \par History
!>      05.2009 created [MI]
!> \author [MI]
! *****************************************************************************
 SUBROUTINE mixing_storage_release(mixing_store,error)
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mixing_storage_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(mixing_store)) THEN
     CPPreconditionNoFail(mixing_store%ref_count>0,cp_failure_level,routineP,error)
     mixing_store%ref_count=mixing_store%ref_count-1
     IF (mixing_store%ref_count==0) THEN

        IF (ASSOCIATED(mixing_store%kerker_factor)) THEN
           DEALLOCATE(mixing_store%kerker_factor,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        IF (ASSOCIATED(mixing_store%special_metric)) THEN
           DEALLOCATE(mixing_store%special_metric,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        IF (ASSOCIATED(mixing_store%pulay_matrix)) THEN
           DEALLOCATE(mixing_store%pulay_matrix,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        IF (ASSOCIATED(mixing_store%rhoin_buffer)) THEN
           DO i = 1,SIZE(mixing_store%rhoin_buffer,2)
             DO j = 1,SIZE(mixing_store%rhoin_buffer,1)
                DEALLOCATE(mixing_store%rhoin_buffer(j,i)%cc,STAT=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END DO
           END DO
           DEALLOCATE(mixing_store%rhoin_buffer,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        IF (ASSOCIATED(mixing_store%paw)) THEN
          DEALLOCATE(mixing_store%paw,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(mixing_store%cpc_h_in)) THEN
          DO j = 1,SIZE(mixing_store%cpc_h_in,2)
          DO k = 1,SIZE(mixing_store%cpc_h_in,1)
            IF(ASSOCIATED(mixing_store%cpc_h_in(k,j)%r_coef)) THEN
              DEALLOCATE(mixing_store%cpc_h_in(k,j)%r_coef,STAT=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
              DEALLOCATE(mixing_store%cpc_s_in(k,j)%r_coef,STAT=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
            END IF
          END DO
          END DO
          DEALLOCATE(mixing_store%cpc_h_in,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(mixing_store%cpc_s_in,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(mixing_store%cpc_h_old)) THEN
          DO j = 1,SIZE(mixing_store%cpc_h_old,2)
          DO k = 1,SIZE(mixing_store%cpc_h_old,1)
            IF(ASSOCIATED(mixing_store%cpc_h_old(k,j)%r_coef)) THEN
              DEALLOCATE(mixing_store%cpc_h_old(k,j)%r_coef,STAT=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
              DEALLOCATE(mixing_store%cpc_s_old(k,j)%r_coef,STAT=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
            END IF
          END DO
          END DO
          DEALLOCATE(mixing_store%cpc_h_old,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(mixing_store%cpc_s_old,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(mixing_store%dcpc_h_in)) THEN
          DO i = 1,SIZE(mixing_store%dcpc_h_in,3)
          DO j = 1,SIZE(mixing_store%dcpc_h_in,2)
          DO k = 1,SIZE(mixing_store%dcpc_h_in,1)
            IF(ASSOCIATED(mixing_store%dcpc_h_in(k,j,i)%r_coef)) THEN
              DEALLOCATE(mixing_store%dcpc_h_in(k,j,i)%r_coef,STAT=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
              DEALLOCATE(mixing_store%dcpc_s_in(k,j,i)%r_coef,STAT=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
            END IF
          END DO
          END DO
          END DO
          DEALLOCATE(mixing_store%dcpc_h_in,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(mixing_store%dcpc_s_in,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(mixing_store%cpc_h_lastres)) THEN
          DO j = 1,SIZE(mixing_store%cpc_h_lastres,2)
          DO k = 1,SIZE(mixing_store%cpc_h_lastres,1)
            IF(ASSOCIATED(mixing_store%cpc_h_lastres(k,j)%r_coef)) THEN
              DEALLOCATE(mixing_store%cpc_h_lastres(k,j)%r_coef,STAT=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
              DEALLOCATE(mixing_store%cpc_s_lastres(k,j)%r_coef,STAT=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
          END DO
          END DO
          DEALLOCATE(mixing_store%cpc_h_lastres,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(mixing_store%cpc_s_lastres,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        IF (ASSOCIATED(mixing_store%res_buffer)) THEN
           DO i = 1,SIZE(mixing_store%res_buffer,2)
             DO j = 1,SIZE(mixing_store%res_buffer,1)
                DEALLOCATE(mixing_store%res_buffer(j,i)%cc,STAT=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END DO
           END DO
           DEALLOCATE(mixing_store%res_buffer,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        IF (ASSOCIATED(mixing_store%norm_res_buffer)) THEN
           DEALLOCATE(mixing_store%norm_res_buffer,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        IF (ASSOCIATED(mixing_store%ig_global_index)) THEN
           DEALLOCATE(mixing_store%ig_global_index,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        IF (ASSOCIATED(mixing_store%drho_buffer)) THEN
           DO i = 1,SIZE(mixing_store%drho_buffer,2)
             DO j = 1,SIZE(mixing_store%drho_buffer,1)
                DEALLOCATE(mixing_store%drho_buffer(j,i)%cc,STAT=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END DO
           END DO
           DEALLOCATE(mixing_store%drho_buffer,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF

        IF (ASSOCIATED(mixing_store%last_res)) THEN
           DO i = 1,SIZE(mixing_store%last_res)
             DEALLOCATE(mixing_store%last_res(i)%cc,STAT=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END DO
           DEALLOCATE(mixing_store%last_res,STAT=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        IF (ASSOCIATED(mixing_store%rhoin)) THEN
           DO i = 1,SIZE(mixing_store%rhoin)
             DEALLOCATE(mixing_store%rhoin(i)%cc,STAT=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END DO
           DEALLOCATE(mixing_store%rhoin,STAT=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        IF (ASSOCIATED(mixing_store%rhoin_old)) THEN
           DO i = 1,SIZE(mixing_store%rhoin_old)
             DEALLOCATE(mixing_store%rhoin_old(i)%cc,STAT=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END DO
           DEALLOCATE(mixing_store%rhoin_old,STAT=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        IF (ASSOCIATED(mixing_store%p_metric))THEN
           DEALLOCATE(mixing_store%p_metric, STAT=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        IF (ASSOCIATED(mixing_store%weight))THEN
           DEALLOCATE(mixing_store%weight, STAT=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        IF (ASSOCIATED(mixing_store%fmat))THEN
           DEALLOCATE(mixing_store%fmat, STAT=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        IF (ASSOCIATED(mixing_store%delta_res))THEN
           DO i = 1,SIZE(mixing_store%delta_res,2)
             DO j = 1,SIZE(mixing_store%delta_res,1)
               DEALLOCATE(mixing_store%delta_res(j,i)%cc, STAT=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END DO
           END DO
           DEALLOCATE(mixing_store%delta_res,STAT=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        IF (ASSOCIATED(mixing_store%u_vec))THEN
           DO i = 1,SIZE(mixing_store%u_vec,2)
             DO j = 1,SIZE(mixing_store%u_vec,1)
               DEALLOCATE(mixing_store%u_vec(j,i)%cc, STAT=stat)
               CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END DO
           END DO
           DEALLOCATE(mixing_store%u_vec,STAT=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        IF (ASSOCIATED(mixing_store%z_vec))THEN
           DO i = 1,SIZE(mixing_store%z_vec,2)
             DO j = 1,SIZE(mixing_store%z_vec,1)
               DEALLOCATE(mixing_store%z_vec(j,i)%cc, STAT=stat)
               CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END DO
           END DO
           DEALLOCATE(mixing_store%z_vec,STAT=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        DEALLOCATE(mixing_store, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
     NULLIFY(mixing_store)
  END IF

 END SUBROUTINE mixing_storage_release

END MODULE qs_density_mixing_types
