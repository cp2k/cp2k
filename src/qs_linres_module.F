!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****m* cp2k/qs_linres_module *
!!
!!   NAME
!!     qs_linres_module
!!
!!   FUNCTION
!!     Contains the setup for  the calculation of properties by linear response
!!     by the application of second order density functional perturbation theory.
!!     The knowledge of the ground state energy, density and wavefunctions is assumed.
!!     Uses the self consistent approach.
!!     Properties that can be calculated : none
!!
!!   AUTHOR
!!      MI 
!!
!!   MODIFICATION HISTORY
!!      created 06-2005 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_module
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       lr_cg, lr_chemshift, lr_diis, lr_ls_2pnt, lr_ls_3pnt, lr_ls_gold, &
       lr_ls_none, lr_none, lr_precond_full_all, lr_precond_none, &
       lr_precond_s_inverse, lr_precond_single, lr_sd
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_linres_methods,               ONLY: linres_localize
  USE qs_linres_nmr_op,                ONLY: nmr_operators
  USE qs_linres_nmr_utils,             ONLY: nmr_do_nmr,&
                                             nmr_env_cleanup,&
                                             nmr_env_init
  USE qs_linres_types,                 ONLY: linres_control_create,&
                                             linres_control_release,&
                                             linres_control_type,&
                                             nmr_env_type
  USE qs_loc_control,                  ONLY: read_loc_section
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_p_env_methods,                ONLY: p_env_create,&
                                             p_env_psi0_changed
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: linres_calculation

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_module'

CONTAINS

! *****************************************************************************
!!****f* cp2k/qs_linres_module/linres_calculation
!!
!!   NAME
!!     linres_calculation
!!
!!   FUNCTION
!!     Driver for the linear response calculatios
!!     Initialize the perturbation environment
!!     Define which properties is to be calculated
!!     Start up the optimization of the response density and wfn
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     06.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE linres_calculation(qs_env, globenv, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_calculation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
    LOGICAL                                  :: failure, ionode, &
                                                lr_calculation, was_present
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(section_vals_type), POINTER         :: lr_section, prop_section

    CALL timeset(routineN,"I"," ",handle)    
    failure = .FALSE. 
    lr_calculation = .FALSE. 
    was_present = .FALSE.

    NULLIFY(dft_control,p_env,linres_control,logger,lr_section,prop_section)

    logger => cp_error_get_logger(error)
    ionode = (logger%para_env%mepos==logger%para_env%source)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    ! The transition potential method to calculate XAS needs LSD
    lr_section =>  section_vals_get_subs_vals(qs_env%input, &
                                          "PROPERTIES%LINRES")
    CALL section_vals_get(lr_section,explicit=lr_calculation,error=error)
    CALL linres_init(lr_section,p_env,qs_env,globenv,error=error)

    IF (lr_calculation) THEN
       output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
            extension=".linresLog",error=error)    
       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, &
            linres_control=linres_control)

       !The type of perturbation has not been defined yet
       linres_control%property = lr_none

       prop_section => section_vals_get_subs_vals(lr_section,"NMR",error=error)
       CALL section_vals_get(prop_section,explicit=was_present,error=error)
       IF(was_present)THEN
         linres_control%property = lr_chemshift
         linres_control%do_kernel = .FALSE.

         IF(.NOT.linres_control%localized_psi0) THEN
            CALL stop_program(routineP,&
             "Are you sure that you want to calculate the chemical shift without localized psi0??")

           CALL linres_localize(linres_control%localized_wfn_control,qs_env,&
                 dft_control%nspins,globenv,centers_only=.TRUE.,error=error)
         END IF

         !Initialize the nmr environment
         nmr_env%ref_count=0
         CALL nmr_env_init(nmr_env,qs_env,globenv,error=error)
         CALL nmr_operators(nmr_env,qs_env,error=error)
         CALL nmr_do_nmr(nmr_env,p_env,qs_env,globenv,error=error)
     WRITE(*,*) 'donmr'
     STOP 'lr_calculation'
         CALL nmr_env_cleanup(nmr_env,error=error)

       END IF

      ! Other possible LR calculations can be introduced here
      
       CALL p_env_release(p_env)

       IF(output_unit>0) THEN
         WRITE (UNIT=output_unit,FMT="(/,T3,A,/,T25,A,/,T3,A,/)")&
            REPEAT("=",77),&
            "END LINRES CALCULATION",&
            REPEAT("=",77)
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
            "PRINT%PROGRAM_RUN_INFO",error=error)
    ELSE
       IF(ionode) THEN 
          output_unit = cp_logger_get_default_unit_nr(logger)
          WRITE (output_unit, "(2X,A)") "",&
               "-----------------------------------------------------------------------------",&
               "-              No LR calculation has been specified in the input            -",&
               "                        cp2k is going to stop, bye bye                       ",&
               "-----------------------------------------------------------------------------",&
               ""    
       END IF
    END IF 

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE linres_calculation

! *****************************************************************************
!!****f* cp2k/qs_linres_module/linres_init
!!
!!   NAME
!!     linres_init
!!
!!   FUNCTION
!!     Initialize some general settings like the p_env
!!     Localize the psi0 if required
!!
!!   NOTES
!!     - The localization should probably be always for all the occupied states
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     06.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE linres_init(lr_section,p_env,qs_env,globenv,error)

    TYPE(section_vals_type), POINTER         :: lr_section
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ispin, istat, j, nspins, &
                                                output_unit
    LOGICAL                                  :: do_it, failure
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(section_vals_type), POINTER         :: loc_section

    failure = .FALSE.

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)    
    NULLIFY(dft_control, linres_control, loc_section, rho, mos, matrix_ks, &
         matrix_s)

    CPPrecondition(.NOT.ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN


       CALL linres_control_create(linres_control,error=error)
       CALL set_qs_env(qs_env=qs_env, linres_control=linres_control)
       CALL linres_control_release(linres_control,error=error)
       CALL get_qs_env(qs_env=qs_env, linres_control=linres_control,&
                       dft_control=dft_control,matrix_ks=matrix_ks,mos=mos,rho=rho)
  
       ! Localized Psi0 are required when the position operator has to be defined (nmr)
       CALL section_vals_val_get(lr_section,"LOCALIZE_PSI0",&
           l_val=linres_control%localized_psi0,error=error)
 
       IF(linres_control%localized_psi0) THEN
         IF(output_unit>0) THEN
           WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
            "Localization of the ground state orbitals before starting the Linres calculation"
         END IF

         loc_section =>section_vals_get_subs_vals(lr_section,"LOCALIZE",error=error)
         CALL read_loc_section(linres_control%localized_wfn_control,loc_section,error=error)
         CALL linres_localize(linres_control%localized_wfn_control,qs_env,&
               dft_control%nspins,globenv,error=error)

         DO ispin=1,dft_control%nspins
            CALL calculate_density_matrix(mos(ispin)%mo_set,&
                     rho%rho_ao(ispin)%matrix,error=error)
         ENDDO
! ** update qs_env%rho
         CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw=dft_control%qs_control%gapw, error=error)
         IF(dft_control%qs_control%gapw ) THEN
           CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
         ENDIF
       END IF

       CALL section_vals_val_get(lr_section,"OPT_METHOD",i_val=linres_control%opt_method,error=error)
       CALL section_vals_val_get(lr_section,"MAX_SCF",i_val=linres_control%max_scf,error=error)
       CALL section_vals_val_get(lr_section,"EPS_SCF",r_val=linres_control%eps_scf,error=error)
       CALL section_vals_val_get(lr_section,"LINESEARCH",i_val=linres_control%ls_method)
       CALL section_vals_val_get(lr_section,"STEPSIZE",r_val=linres_control%ds_min)

       CALL section_vals_val_get(lr_section,"PRECONDITIONER",i_val=linres_control%preconditioner_type)
       CALL section_vals_val_get(lr_section,"ENERGY_GAP",r_val=linres_control%energy_gap)
       CALL section_vals_val_get(lr_section,"CHECK_ORTHO",&
            l_val=linres_control%check_ortho,error=error)
       
       IF(output_unit>0) THEN
         WRITE (UNIT=output_unit,FMT="(/,T3,A,/,T25,A,/,T3,A,/)")&
            REPEAT("=",77),&
            "START LINRES CALCULATION",&
            REPEAT("=",77)

         WRITE (UNIT=output_unit,FMT="(/,T10,A,/)")&
            "Properties to be Calulated:" 
         CALL section_vals_val_get(lr_section,"NMR%_SECTION_PARAMETERS_",&
            l_val=do_it,error=error)
         IF(do_it) WRITE (UNIT=output_unit,FMT="(T45,A)") &
             "NMR Chemical Shift"

         IF(linres_control%localized_psi0) WRITE (UNIT=output_unit,FMT="(T2,A,T65,A)")&
            "LINRES|","LOCALIZED PSI0" 
         IF(linres_control%opt_method==lr_cg) WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "LINRES|     Optimization algorithm","Conjugate Gradients"
         IF(linres_control%opt_method==lr_sd) WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "LINRES|     Optimization algorithm","    Steepnd Descent"
         IF(linres_control%opt_method==lr_diis) WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "LINRES|     Optimization algorithm","               DIIS"
         IF(linres_control%ls_method==lr_ls_none) &
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "LINRES|                Line Search","               NONE"
         IF(linres_control%ls_method==lr_ls_2pnt) &
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "LINRES|                Line Search","               2pnt"
         IF(linres_control%ls_method==lr_ls_3pnt) &
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "LINRES|                Line Search","               3pnt"
         IF(linres_control%ls_method==lr_ls_gold) &
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "LINRES|                Line Search","               gold"
         IF(linres_control%preconditioner_type==lr_precond_none) &
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "LINRES|             Preconditioner","               NONE"
         IF(linres_control%preconditioner_type==lr_precond_single) &
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "LINRES|             Preconditioner","        FULL SINGLE"
         IF(linres_control%preconditioner_type==lr_precond_s_inverse) &
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "LINRES|             Preconditioner","     FULL S_INVERSE"
         IF(linres_control%preconditioner_type==lr_precond_full_all) &
            WRITE (UNIT=output_unit,FMT="(T2,A,T60,A)")&
            "LINRES|             Preconditioner","           FULL ALL"
         WRITE (UNIT=output_unit,FMT="(/,T2,A,T71,ES8.1)")&
            "LINRES|                    EPS_SCF",linres_control%eps_scf
         WRITE (UNIT=output_unit,FMT="(/,T2,A,T71,I8)")&
            "LINRES|    Max number of iterations per SCF cycle",linres_control%max_scf
       END IF

       
      !------------------!
      ! create the p_env !
      !------------------!
       CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE.,linres_control=linres_control,error=error)

       ! update the m_epsilon matrix
       CALL p_env_psi0_changed(p_env,qs_env,error=error)

       ! calculate eigenvectros and eigenvalues of K
       p_env%os_valid = .FALSE.
       IF(linres_control%preconditioner_type/= lr_precond_none) THEN
         p_env%new_preconditioner = .TRUE.
         CALL get_qs_env(qs_env=qs_env,matrix_ks=matrix_ks,&
              matrix_s=matrix_s,mos=mos)
         nspins = SIZE(p_env%preconditioner,1)
         DO ispin = 1,nspins
           CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
           p_env%ev_h0(ispin)%matrix => mo_coeff
!           CALL cp_fm_to_fm(mo_coeff,p_env%ev_h0(ispin)%matrix,p_env%n_mo(ispin),1,1) 
           ALLOCATE(p_env%preconditioner(ispin)%occ_evals(p_env%n_mo(ispin)),STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
           CALL calculate_subspace_eigenvalues(p_env%ev_h0(ispin)%matrix,&
                matrix_ks(ispin)%matrix,evals_arg=p_env%preconditioner(ispin)%occ_evals,&
                para_env=mo_coeff%matrix_struct%para_env, do_rotation = .FALSE.)

           IF(output_unit>0) THEN
             DO i = 1,p_env%n_mo(ispin),4
               j=MIN(3,p_env%n_mo(ispin)-i)
               SELECT CASE (j)
               CASE(3)
                 WRITE(output_unit,'(1X,4F16.8)') p_env%preconditioner(ispin)%occ_evals(i:i+j)
               CASE(2)
                 WRITE(output_unit,'(1X,3F16.8)') p_env%preconditioner(ispin)%occ_evals(i:i+j)
               CASE(1)
                 WRITE(output_unit,'(1X,2F16.8)') p_env%preconditioner(ispin)%occ_evals(i:i+j)
               CASE(0)
                 WRITE(output_unit,'(1X,1F16.8)') p_env%preconditioner(ispin)%occ_evals(i:i+j)
               END SELECT
             END DO
           END IF
         END DO 
       END IF
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)    
  END SUBROUTINE linres_init


! *****************************************************************************
END MODULE qs_linres_module


