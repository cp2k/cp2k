!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/qs_linres_nmr_utils *
!!
!!   NAME
!!     qs_linres_nmr_utils
!!
!!   FUNCTION
!!     Chemical shift calculation by dfpt
!!     Initialization of the nmr_env, creation of the special neighbor lists
!!     Perturbation Hamiltonians by application of the p and rxp oprtators to  psi0
!!     Write output
!!     Deallocate everything 
!!
!!   NOTE
!!     The psi0 should be localized
!!     the Sebastiani method works within the assumption that the orbitals are 
!!     completely contained in the simulation box 
!!  
!!   AUTHOR
!!      MI 
!!
!!   MODIFICATION HISTORY
!!      created 07-2005 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_nmr_utils

  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm,&
                                             fm_pool_get_el_struct,&
                                             fm_pool_give_back_fm
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_methods,               ONLY: linres_localize
  USE qs_linres_types,                 ONLY: linres_control_type,&
                                             nmr_env_create, nmr_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE simulation_cell,                 ONLY: cell_type, pbc
  USE sparse_matrix_types,             ONLY: cp_sm_scale_and_add,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             transfer_matrix

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: nmr_env_cleanup, nmr_env_init

  CHARACTER(LEN=*), PARAMETER :: moduleN='qs_linres_nmr_utils'

CONTAINS


! *****************************************************************************
!!****f* cp2k/qs_linres_nmr_utils/nmr_env_init
!!
!!   NAME
!!     nmr_env_init
!!
!!   FUNCTION
!!     Initialize the nmr environment
!!
!!   NOTES
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE nmr_env_init(nmr_env,qs_env,globenv,error)

    TYPE(nmr_env_type)                          :: nmr_env
    TYPE(qs_environment_type), POINTER          :: qs_env
    TYPE(global_environment_type), &
      POINTER                                   :: globenv
    TYPE(cp_error_type), & 
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_env_init', &
                                   routineP = moduleN//'/'//routineN
 
    INTEGER                                     :: homo, idir,ispin,istat,n_mo(2),nao,nmoloc,output_unit
    LOGICAL                                     :: ionode, failure
    REAL(dp), DIMENSION(:,:), POINTER           :: centers
    TYPE(cp_logger_type), POINTER               :: logger
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(linres_control_type), POINTER          :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                   :: mos
    TYPE(qs_matrix_pools_type), POINTER         :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                   :: matrix_s
    TYPE(section_vals_type), POINTER            :: nmr_section

    failure = .FALSE.

    NULLIFY(centers,dft_control,linres_control,logger,matrix_s,mos,mpools,nmr_section)
    n_mo(1:2) = 0
    nao = 0
    nmoloc = 0

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit = cp_logger_get_default_unit_nr(logger)

    IF(nmr_env%ref_count /= 0) THEN
      CALL nmr_env_cleanup(nmr_env,error=error)
    END IF

    CALL nmr_env_create(nmr_env,error=error)

    ! If surrent_density or full_nmr different allocations are required
    nmr_section => section_vals_get_subs_vals(globenv%input_file, &
                                                 "FORCE_EVAL%LINRES%NMR")
    CALL section_vals_val_get(nmr_section,"CURRENT_DENSITY", &
         l_val=nmr_env%current_density,error=error)
    CALL section_vals_val_get(nmr_section,"FULL",l_val=nmr_env%full_nmr,error=error)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    linres_control=linres_control,&
                    matrix_s=matrix_s,&
                    mos=mos,&
                    mpools=mpools)

    ! check that the psi0 are localized and you have all the centers
     CPPrecondition(linres_control%localized_psi0,cp_warning_level,routineP,error,failure)
     IF(failure .AND. ionode) THEN
       WRITE (UNIT=output_unit,'(A)') &
          ' To get NMR parameters within PBC you need localized zero order orbitals '
     END IF
     DO ispin = 1,dft_control%nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=n_mo(ispin),homo=homo)
       centers => linres_control%localized_wfn_control%centers_set(ispin)%array
       nmoloc = size(centers,2)
       CPPrecondition(nmoloc==homo,cp_warning_level,routineP,error,failure)
     END DO  

    ! point to the psi0 centers
     nmr_env%centers_set =  linres_control%localized_wfn_control%centers_set

     CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=nmr_env%ao_mo_fm_pools, error=error)
    ! for the chemical shift we need 6 psi1, i.e. 6 optimization procedures
    ! They become 9 if full nmr is calculated, i.e. with the correction term too
    ! All of them are required at the end of the optimization procedure
    ! if the current density and the induced fiels have to be calculated 
    ! If instead only the shift is needed, only one psi1 should be enough, providing
    ! that after every optimization the corresponding shift contribution is calculated
    ! prepare the psi1
     IF(nmr_env%current_density) THEN

       ALLOCATE(nmr_env%psi1_p(dft_control%nspins,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(nmr_env%psi1_rxp(dft_control%nspins,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ispin = 1,dft_control%nspins
         DO idir = 1,3
           NULLIFY(nmr_env%psi1_p(ispin,idir)%matrix)
           CALL fm_pool_create_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,nmr_env%psi1_p(ispin,idir)%matrix)
           NULLIFY(nmr_env%psi1_rxp(ispin,idir)%matrix)
           CALL fm_pool_create_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,nmr_env%psi1_rxp(ispin,idir)%matrix)
         END DO 
       END DO 

       IF(nmr_env%full_nmr) THEN
         ALLOCATE(nmr_env%psi1_D(dft_control%nspins,3),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DO ispin = 1,dft_control%nspins
           DO idir = 1,3
             NULLIFY(nmr_env%psi1_D(ispin,idir)%matrix)
             CALL fm_pool_create_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,nmr_env%psi1_D(ispin,idir)%matrix)
           END DO
         END DO
       END IF
     END IF

     ! px py pz operator acting on the atomic orbitals
     ALLOCATE(nmr_env%op_p_ao(3),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO idir = 1,3 
       NULLIFY(nmr_env%op_p_ao(idir)%matrix)
       CALL replicate_matrix_structure(matrix_s(1)%matrix, &
            nmr_env%op_p_ao(idir)%matrix,"nmr_env%op_p_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))))
       CALL set_matrix(nmr_env%op_p_ao(idir)%matrix,0.0_dp)
     END DO

     ! for the rxp we cannot calculate it a priori because it is in facts (r-dk)xp 
     ! where dk is the center of the orbital it is applied to. We would need nstate operators
     ! What we can store is (r-dk)xp|psi0k> for each k, which is a full matrix only
     ! Therefore we prepare here the full matrix p_psi0 and rxp_psi0
     ! We also need a temporary sparse matrix where to store the integrals during the calculation
     ALLOCATE(nmr_env%p_psi0(dft_control%nspins,3),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(nmr_env%rxp_psi0(dft_control%nspins,3),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO ispin = 1,dft_control%nspins
       DO idir = 1,3
         NULLIFY(nmr_env%p_psi0(ispin,idir)%matrix)
         CALL fm_pool_create_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,nmr_env%p_psi0(ispin,idir)%matrix)
         NULLIFY(nmr_env%rxp_psi0(ispin,idir)%matrix)
         CALL fm_pool_create_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,nmr_env%rxp_psi0(ispin,idir)%matrix)
       END DO 
     END DO 

     ALLOCATE(nmr_env%op_rmd_ao(3),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO idir = 1,3 
       NULLIFY(nmr_env%op_rmd_ao(idir)%matrix)
       CALL replicate_matrix_structure(matrix_s(1)%matrix, &
            nmr_env%op_rmd_ao(idir)%matrix,"nmr_env%op_rmd_ao"//&
            "-"//TRIM(ADJUSTL(cp_to_string(idir))))
       CALL set_matrix(nmr_env%op_rmd_ao(idir)%matrix,0.0_dp)
     END DO


    ! Current density:
    ! good question, it should have the same structur of the density
    ! but how we define it depends also on how we will use it
    ! If Gapw I suppose we should use the same partitioning used for the density
   


  END SUBROUTINE nmr_env_init

! *****************************************************************************
!!****f* cp2k/qs_linres_nmr_utils/nmr_env_cleanup
!!
!!   NAME
!!     nmr_env_cleanup
!!
!!   FUNCTION
!!     Deallocate the nmr environment
!!
!!   NOTES
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE nmr_env_cleanup(nmr_env,error)

    TYPE(nmr_env_type)                          :: nmr_env
    TYPE(cp_error_type), & 
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_env_cleanup', &
                                   routineP = moduleN//'/'//routineN
 
    INTEGER                                     :: idir,ispin,istat
    LOGICAL                                     :: failure

    IF(.NOT. failure) THEN
      nmr_env%ref_count = nmr_env%ref_count - 1
      IF(nmr_env%ref_count == 0 ) THEN
        !psi1_p
        IF(ASSOCIATED(nmr_env%psi1_p)) THEN
          DO idir = 1,SIZE(nmr_env%psi1_p,2)
            DO ispin = 1,SIZE(nmr_env%psi1_p,1)
              CALL fm_pool_give_back_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,&
                   nmr_env%psi1_p(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(nmr_env%psi1_p, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !psi1_rxp
        IF(ASSOCIATED(nmr_env%psi1_rxp)) THEN
          DO idir = 1,SIZE(nmr_env%psi1_rxp,2)
            DO ispin = 1,SIZE(nmr_env%psi1_rxp,1)
              CALL fm_pool_give_back_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,&
                   nmr_env%psi1_rxp(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(nmr_env%psi1_rxp, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !psi1_D
        IF(ASSOCIATED(nmr_env%psi1_D)) THEN
          DO idir = 1,SIZE(nmr_env%psi1_D,2)
            DO ispin = 1,SIZE(nmr_env%psi1_D,1)
              CALL fm_pool_give_back_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,&
                   nmr_env%psi1_D(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(nmr_env%psi1_D, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !p_psi0
        IF(ASSOCIATED(nmr_env%p_psi0)) THEN
          DO idir = 1,SIZE(nmr_env%p_psi0,2)
            DO ispin = 1,SIZE(nmr_env%p_psi0,1)
              CALL fm_pool_give_back_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,&
                   nmr_env%p_psi0(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(nmr_env%p_psi0, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        !rxp_psi0
        IF(ASSOCIATED(nmr_env%rxp_psi0)) THEN
          DO idir = 1,SIZE(nmr_env%rxp_psi0,2)
            DO ispin = 1,SIZE(nmr_env%rxp_psi0,1)
              CALL fm_pool_give_back_fm(nmr_env%ao_mo_fm_pools(ispin)%pool,&
                   nmr_env%rxp_psi0(ispin,idir)%matrix,error=error)
            END DO
          END DO
          DEALLOCATE(nmr_env%rxp_psi0, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        NULLIFY(nmr_env%ao_mo_fm_pools)
        NULLIFY(nmr_env%centers_set)
        !op_p_ao
        IF(ASSOCIATED(nmr_env%op_p_ao)) THEN
          CALL deallocate_matrix_set(nmr_env%op_p_ao)
        END IF
        !op_p_ao
        IF(ASSOCIATED(nmr_env%op_rmd_ao)) THEN
          CALL deallocate_matrix_set(nmr_env%op_rmd_ao)
        END IF
      END IF
    END IF
    
  END SUBROUTINE nmr_env_cleanup


END MODULE  qs_linres_nmr_utils
