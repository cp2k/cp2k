!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Does all kind of post scf calculations for SCPTB
!> \par History
!>      Started as a copy from the relevant part of qs_scf_post_gpw
!> \author Ralph Koitz (09.2013)
! *****************************************************************************
MODULE qs_scf_post_scptb

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE population_analyses,             ONLY: mulliken_population_analysis
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                             qs_subsys_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  ! Global parameters
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_post_scptb'
  PUBLIC :: scf_post_calculation_scptb

CONTAINS

! *****************************************************************************
!> \brief collects possible post - scf calculations and prints info / computes properties.
!> \param dft_section ...
!> \param scf_env the scf_env whose info should be written out
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>        09.2013  Adapted from qs_scf_post_gpw.F
!> \author  Ralph Koitz
!> \note
! *****************************************************************************
  SUBROUTINE scf_post_calculation_scptb(dft_section, scf_env,qs_env, error)

    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_post_calculation_scptb', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL write_scptb_mo_free_results(qs_env,error)  

    CALL timestop(handle)
  END SUBROUTINE scf_post_calculation_scptb

! *****************************************************************************
!> \brief Write QS results always available in SCPTB calculation
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History 09.2013 adapted from qs_scf_post_gpw.F. [RK] 
! *****************************************************************************
  SUBROUTINE write_scptb_mo_free_results(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_scptb_mo_free_results', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, iw, &
                                                output_unit, print_level, &
                                                unit_nr
    LOGICAL                                  :: failure, print_it, write_ks, &
                                                write_xc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: dft_section, input, print_key

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(dft_control, rho, ks_rmpv, input)
    para_env=>qs_env%para_env
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,dft_control=dft_control, &
            atomic_kind_set=atomic_kind_set,particle_set=particle_set,&
            rho=rho,matrix_ks=ks_rmpv,&
            matrix_s=matrix_s,input=input,cell=cell,subsys=subsys,error=error)
       dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
       CALL qs_subsys_get(subsys,particles=particles,error=error)

       ! Print the total density (electronic + core charge)
       ! Not implemented so far.
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%TOT_DENSITY_CUBE", error=error),cp_p_file)) THEN
           CALL cp_unimplemented_error(fromWhere=routineP, &
                message="Printing of cube files from SCPTB calculation currently " //&
                "not supported", error=error, error_level=cp_warning_level)
       END IF

       ! Write cube file with electron density
       ! Not implemented so far.
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%E_DENSITY_CUBE",error=error),cp_p_file)) THEN
           CALL cp_unimplemented_error(fromWhere=routineP, &
                message="Printing of cube files from SCPTB calculation currently " //&
                "not supported", error=error, error_level=cp_warning_level)
       END IF ! print key


       ! Print the hartree potential
       ! Not implemented so far.
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%V_HARTREE_CUBE",error=error),cp_p_file)) THEN
           CALL cp_unimplemented_error(fromWhere=routineP, &
                message="Printing of cube files from SCPTB calculation currently " //&
                "not supported", error=error, error_level=cp_warning_level)
       ENDIF

       ! Print the Electrical Field Components
       ! Not implemented so far.
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%EFIELD_CUBE",error=error),cp_p_file)) THEN

           CALL cp_unimplemented_error(fromWhere=routineP, &
                message="Printing of cube files from SCPTB calculation currently " //&
                "not supported", error=error, error_level=cp_warning_level)
       END IF

       ! Write the density matrices
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%AO_MATRICES/DENSITY",error=error),cp_p_file)) THEN
          iw = cp_print_key_unit_nr(logger,input,"DFT%PRINT%AO_MATRICES/DENSITY",&
                                    extension=".Log",error=error)
          DO ispin=1,dft_control%nspins
             CALL cp_dbcsr_write_sparse_matrix(rho%rho_ao(ispin)%matrix,4,6,qs_env,&
                                            para_env,output_unit=iw,error=error)
          END DO
          CALL cp_print_key_finished_output(iw,logger,input,&
                                            "DFT%PRINT%AO_MATRICES/DENSITY",&
                                            error=error)
       END IF

       ! Write the Kohn-Sham matrices
       write_ks=BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",error=error),cp_p_file)
       write_xc=BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%AO_MATRICES/MATRIX_VXC",error=error),cp_p_file)
       ! we need to update stuff before writing, potentially computing the matrix_vxc
       IF (write_ks .OR. write_xc) THEN
          IF (write_xc) qs_env%requires_matrix_vxc=.TRUE.
          CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,error=error)
          CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,error=error,&
                                    calculate_forces=.FALSE.,just_energy=.FALSE.)
          IF (write_xc) qs_env%requires_matrix_vxc=.FALSE.
       END IF

       ! Write the Kohn-Sham matrices
       IF (write_ks) THEN
          iw = cp_print_key_unit_nr(logger,input,"DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",&
                                    extension=".Log",error=error)
          DO ispin=1,dft_control%nspins
             CALL cp_dbcsr_write_sparse_matrix(ks_rmpv(ispin)%matrix,4,6,qs_env,&
                                            para_env,output_unit=iw,error=error)
          END DO
          CALL cp_print_key_finished_output(iw,logger,input,&
                                            "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",&
                                            error=error)
       END IF

       ! If user asked for XC Matrix, fail.
       IF (write_xc) THEN
           CALL cp_unimplemented_error(fromWhere=routineP, &
                message="Printing of XC Matrix not possible in SCPTB calculation.", &
                error=error, error_level=cp_warning_level)
       END IF

       ! Compute the Mulliken charges
       print_key => section_vals_get_subs_vals(input,"DFT%PRINT%MULLIKEN", error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MULLIKEN",extension=".mulliken",&
                                       middle_name="",log_filename=.FALSE.,error=error)
          print_level = 1
          CALL section_vals_val_get(print_key,"PRINT_GOP",l_val=print_it,error=error)
          IF (print_it) print_level = 2
          CALL section_vals_val_get(print_key,"PRINT_ALL",l_val=print_it,error=error)
          IF (print_it) print_level = 3
          CALL mulliken_population_analysis(qs_env,unit_nr,print_level,error)
          CALL cp_print_key_finished_output(unit_nr, logger,input,"DFT%PRINT%MULLIKEN",error=error)
       END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE write_scptb_mo_free_results

END MODULE qs_scf_post_scptb
