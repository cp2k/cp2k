!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Does all kind of post scf calculations for semi-empirical
!> \par History
!>      Started printing preliminary stuff for MO_CUBES and MO requires some
!>      more work to complete all other functionalities
!> \author Teodoro Laino (07.2008)
! *****************************************************************************
MODULE qs_scf_post_se
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_nfullrows_total,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_sum
  USE mulliken,                        ONLY: multipoles
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_mo_io,                        ONLY: write_mo_set
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                             qs_subsys_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  ! Global parameters
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_post_se'
  PUBLIC :: scf_post_calculation_se

CONTAINS

! *****************************************************************************
!> \brief collects possible post - scf calculations and prints info / computes properties.
!>        specific for Semi-empirical calculations
!> \param dft_section ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2008 created [tlaino] - Splitted from qs_scf_post (general)
!> \author tlaino
!> \note
!>      this function changes mo_eigenvectors and mo_eigenvalues, depending on the print keys.
!>      In particular, MO_CUBES causes the MOs to be rotated to make them eigenstates of the KS
!>      matrix, and mo_eigenvalues is updated accordingly. This can, for unconverged wavefunctions,
!>      change afterwards slightly the forces (hence small numerical differences between MD
!>      with and without the debug print level). Ideally this should not happen...
! *****************************************************************************
  SUBROUTINE scf_post_calculation_se(dft_section,qs_env,error)

    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_post_calculation_se', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
    LOGICAL                                  :: explicit, failure, &
                                                my_localized_wfn
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: input, print_key, &
                                                wfn_mix_section

    CALL timeset(routineN,handle)

    ! Writes the data that is already available in qs_env
    CALL write_available_results(qs_env,error)

    failure=.FALSE.
    my_localized_wfn = .FALSE.
    NULLIFY(dft_control, mos, rho, matrix_s, &
         subsys, particles, input, print_key, para_env)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    ! Here we start with data that needs a postprocessing...
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,&
                       dft_control=dft_control,&
                       mos=mos,&
                       rho=rho,&
                       matrix_s=matrix_s,&
                       input=input,&
                       subsys=subsys,&
                       para_env=para_env,&
                       error=error)
       CALL qs_subsys_get(subsys,particles=particles,error=error)

       ! Compute Atomic Charges
       CALL qs_scf_post_charges(input, logger, qs_env, rho, matrix_s, mos, para_env, error)

       ! Moments of charge distribution
       CALL qs_scf_post_moments(input, logger, qs_env, output_unit, error)

       ! MO_CUBES
       print_key => section_vals_get_subs_vals(section_vals=input,&
                      subsection_name="DFT%PRINT%MO_CUBES",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
         CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Printing of MO cube files not implemented for Semi-Empirical method.", &
            error=error, error_level=cp_warning_level)
       END IF

       ! STM
       print_key => section_vals_get_subs_vals(section_vals=input,&
                      subsection_name="DFT%PRINT%STM",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
         CALL cp_unimplemented_error(fromWhere=routineP, &
            message="STM not implemented for Semi-Empirical method.", &
            error=error, error_level=cp_warning_level)
       END IF

       ! DFT+U
       print_key => section_vals_get_subs_vals(section_vals=input,&
                      subsection_name="DFT%PRINT%PLUS_U",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
         CALL cp_unimplemented_error(fromWhere=routineP, &
            message="DFT+U not available for Semi-Empirical method.", &
            error=error, error_level=cp_warning_level)
       END IF

       ! Kinetic Energy
       print_key => section_vals_get_subs_vals(section_vals=input,&
                      subsection_name="DFT%PRINT%KINETIC_ENERGY",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
         CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Kinetic energy not available for Semi-Empirical method.", &
            error=error, error_level=cp_warning_level)
       END IF

       ! Wavefunction mixing
       wfn_mix_section => section_vals_get_subs_vals(input,"DFT%PRINT%WFN_MIX",error=error)
       CALL section_vals_get(wfn_mix_section,explicit=explicit,error=error)
       IF(explicit.AND..NOT.qs_env%run_rtp) THEN
         CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Wavefunction mixing not implemented for Semi-Empirical  method.", &
            error=error, error_level=cp_warning_level)
       END IF

       ! Print coherent X-ray diffraction spectrum
       print_key => section_vals_get_subs_vals(section_vals=input,&
                    subsection_name="DFT%PRINT%XRAY_DIFFRACTION_SPECTRUM",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
            message="XRAY_DIFFRACTION_SPECTRUM  not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
       END IF

       ! Calculation of Electric Field Gradients
       print_key => section_vals_get_subs_vals(section_vals=input,&
                    subsection_name="DFT%PRINT%ELECTRIC_FIELD_GRADIENT",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
         CALL cp_unimplemented_error(fromWhere=routineP, &
            message="ELECTRIC_FIELD_GRADIENT not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
       END IF

       ! Calculation of EPR Hyperfine Coupling Tensors
       print_key => section_vals_get_subs_vals(section_vals=input,&
                    subsection_name="DFT%PRINT%HYPERFINE_COUPLING_TENSOR",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),&
              cp_p_file)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
            message="HYPERFINE_COUPLING_TENSOR  not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
       END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE scf_post_calculation_se

! *****************************************************************************
!> \brief Computes and prints electric moments
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param output_unit ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_moments(input, logger, qs_env, output_unit, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_moments', &
      routineP = moduleN//':'//routineN

    TYPE(section_vals_type), POINTER         :: print_key

    print_key => section_vals_get_subs_vals(section_vals=input,&
         subsection_name="DFT%PRINT%MOMENTS",error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="ELECTRIC MOMENTS not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
    END IF

  END SUBROUTINE qs_scf_post_moments

! *****************************************************************************
!> \brief Computes and Prints Atomic Charges with several methods
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param rho ...
!> \param matrix_s ...
!> \param mos ...
!> \param para_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_charges(input, logger, qs_env, rho, matrix_s, mos, &
                                 para_env, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_charges', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: ana
    CHARACTER(LEN=default_string_length)     :: aname
    INTEGER                                  :: i, iat, iatom, ikind, istat, &
                                                j, nao, nat, natom, natorb, &
                                                nkind, nspin, unit_nr
    LOGICAL                                  :: failure, found, scp_nddo
    REAL(KIND=dp)                            :: zeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: mcharge
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: charges
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pblock
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p
    TYPE(cp_dbcsr_type), POINTER             :: pscp
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: print_key
    TYPE(semi_empirical_type), POINTER       :: se_kind

    NULLIFY(particle_set)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    natom=natom,&
                    qs_kind_set=qs_kind_set,&
                    local_particles=local_particles,&
                    particle_set=particle_set,&
                    scp_env=scp_env,&
                    dft_control=dft_control,&
                    error=error)

    scp_nddo = dft_control%qs_control%se_control%scp
    IF (scp_nddo ) CALL get_scp_env(scp_env=scp_env, pscp=pscp, error=error)

    ! Compute the mulliken charges
    print_key => section_vals_get_subs_vals(input,"DFT%PRINT%MULLIKEN", error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MULLIKEN",extension=".mulliken",&
                                    middle_name="",log_filename=.FALSE.,error=error)
       matrix_p => rho%rho_ao
       nao = cp_dbcsr_nfullrows_total(matrix_p(1)%matrix)
       nspin = SIZE(matrix_p)
       ALLOCATE (charges(natom,nspin),mcharge(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       charges = 0.0_dp
       mcharge = 0.0_dp
       ! calculate atomic charges
       nkind = SIZE(atomic_kind_set)
       DO ikind=1,nkind
          CALL get_atomic_kind(atomic_kind_set(ikind),natom=nat)
          CALL get_qs_kind(qs_kind_set(ikind),se_parameter=se_kind)
          CALL get_se_param(se_kind,zeff=zeff,natorb=natorb)
          DO iatom=1,nat
             iat = atomic_kind_set(ikind)%atom_list(iatom)
             DO i=1,nspin
                CALL cp_dbcsr_get_block_p(matrix=matrix_p(i)%matrix,row=iat,col=iat,&
                  block=pblock,found=found)
                IF(found) THEN
                   DO j=1,natorb
                      charges(iat,i) = charges(iat,i) + pblock(j,j)
                   END DO
                END IF
             END DO
             mcharge(iat) = zeff - SUM(charges(iat,1:nspin))
          END DO
       END DO
       !
       CALL mp_sum(charges,para_env%group)
       CALL mp_sum(mcharge,para_env%group)
       !
       IF (unit_nr > 0) THEN
          WRITE (UNIT=unit_nr,FMT="(/,/,T2,A)") "POPULATION ANALYSIS"
          IF (nspin == 1) THEN
             WRITE (UNIT=unit_nr,FMT="(/,T2,A,T70,A)")&
                    " # Atom   Element   Kind        Atomic population"," Net charge"
             DO ikind=1,nkind
                CALL get_atomic_kind(atomic_kind_set(ikind),natom=nat)
                CALL get_qs_kind(qs_kind_set(ikind),se_parameter=se_kind)
                CALL get_se_param(se_kind,name=aname)
                ana = ADJUSTR(TRIM(ADJUSTL(aname)))
                DO iatom=1,nat
                   iat = atomic_kind_set(ikind)%atom_list(iatom)
                   WRITE (UNIT=unit_nr,&
                          FMT="(T2,I7,6X,A2,3X,I6,T39,F12.6,T69,F12.6)")&
                          iat,ana,ikind,charges(iat,1),mcharge(iat)
                END DO
             END DO
             WRITE (UNIT=unit_nr,&
                    FMT="(T2,A,T39,F12.6,T69,F12.6,/)")&
                        "# Total charge",SUM(charges(:,1)),SUM(mcharge(:))
          ELSE
             WRITE (UNIT=unit_nr,FMT="(/,T2,A)")&
                    "# Atom  Element  Kind  Atomic population (alpha,beta)   Net charge  Spin moment"
             DO ikind=1,nkind
                CALL get_atomic_kind(atomic_kind_set(ikind),natom=nat)
                CALL get_qs_kind(qs_kind_set(ikind),se_parameter=se_kind)
                CALL get_se_param(se_kind,name=aname)
                ana = ADJUSTR(TRIM(ADJUSTL(aname)))
                DO iatom=1,nat
                   iat = atomic_kind_set(ikind)%atom_list(iatom)
                   WRITE (UNIT=unit_nr,&
                          FMT="(T2,I6,5X,A2,2X,I6,T29,4(1X,F12.6))")&
                          iat,ana,ikind,charges(iat,1:2),mcharge(iat),charges(iat,1)-charges(iat,2)
                END DO
             END DO
             WRITE (UNIT=unit_nr,&
                    FMT="(T2,A,T29,4(1X,F12.6),/)")&
                        "# Total charge and spin",SUM(charges(:,1)),SUM(charges(:,2)),SUM(mcharge(:))
          END IF
       END IF
       IF ( scp_nddo ) THEN
          CALL multipoles(atomic_kind_set,qs_kind_set,particle_set,&
            rho%rho_ao,pscp,matrix_s(1)%matrix,para_env,unit_nr,"DFT%PRINT%MULLIKEN",error=error)
       END IF

       CALL cp_print_key_finished_output(unit_nr, logger,input,"DFT%PRINT%MULLIKEN",error=error)

       DEALLOCATE (charges,mcharge,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    ! Compute the Lowdin charges
    print_key => section_vals_get_subs_vals(input,"DFT%PRINT%LOWDIN", error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP,&
            message="Lowdin charges not available for semi-empirical calculations!", &
            error=error, error_level=cp_warning_level)
    END IF

  END SUBROUTINE qs_scf_post_charges

! *****************************************************************************
!> \brief Write QS results always available (if switched on through the print_keys)
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE write_available_results(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_available_results', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, iw, output_unit
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: dft_section, input

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(cell, dft_control, mos, atomic_kind_set, particle_set, rho, &
         ks_rmpv, dft_section, input, &
         particles, subsys, para_env)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,&
                       dft_control=dft_control,&
                       mos=mos,&
                       atomic_kind_set=atomic_kind_set,&
                       qs_kind_set=qs_kind_set,&
                       particle_set=particle_set,&
                       rho=rho,&
                       matrix_ks=ks_rmpv,&
                       input=input,&
                       cell=cell,&
                       subsys=subsys,&
                       para_env=para_env,&
                       error=error)
       CALL qs_subsys_get(subsys,particles=particles,error=error)

       ! *** if the dft_section tells you to do so, write last wavefunction to screen
       dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
       IF (dft_control%nspins == 2) THEN
         CALL write_mo_set(mos(1)%mo_set,atomic_kind_set,qs_kind_set,particle_set,4,&
                           dft_section,spin="ALPHA",last=.TRUE.,error=error)
         CALL write_mo_set(mos(2)%mo_set,atomic_kind_set,qs_kind_set,particle_set,4,&
                           dft_section,spin="BETA",last=.TRUE.,error=error)
       ELSE
         CALL write_mo_set(mos(1)%mo_set,atomic_kind_set,qs_kind_set,particle_set,4,&
                           dft_section,last=.TRUE.,error=error)
       END IF

       ! *** at the end of scf print out the projected dos per kind
       IF (BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%PDOS",&
           error=error),cp_p_file) ) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="PDOS not implemented for Semi-Empirical calculations!!", &
               error=error, error_level=cp_warning_level)
       ENDIF

       ! Print the total density (electronic + core charge)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%TOT_DENSITY_CUBE", error=error),cp_p_file)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="TOT_DENSITY_CUBE  not implemented for Semi-Empirical calculations!!", &
               error=error, error_level=cp_warning_level)
       END IF

       ! Write cube file with electron density
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%E_DENSITY_CUBE",error=error),cp_p_file)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="E_DENSITY_CUBE not implemented for Semi-Empirical calculations!!", &
               error=error, error_level=cp_warning_level)
       END IF ! print key

       ! Write cube file with EFIELD 
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%EFIELD_CUBE",error=error),cp_p_file)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="EFIELD_CUBE not implemented for Semi-Empirical calculations!!", &
               error=error, error_level=cp_warning_level)
       END IF ! print key

       ! Write cube file with ELF
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%ELF_CUBE",error=error),cp_p_file)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="ELF function not implemented for Semi-Empirical calculations!!", &
               error=error, error_level=cp_warning_level)
       END IF ! print key

       ! Print the hartree potential
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%V_HARTREE_CUBE",error=error),cp_p_file)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="V_HARTREE_CUBE not implemented for Semi-Empirical calculations!!", &
               error=error, error_level=cp_warning_level)
       ENDIF

       ! Print the XC potential
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%V_XC_CUBE",error=error),cp_p_file)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="V_XC_CUBE not available for Semi-Empirical calculations!!", &
               error=error, error_level=cp_warning_level)
       ENDIF

       !    *** write the density matrix ***
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%AO_MATRICES/DENSITY",error=error),cp_p_file)) THEN
          iw = cp_print_key_unit_nr(logger,input,"DFT%PRINT%AO_MATRICES/DENSITY",&
               extension=".Log",error=error)
          DO ispin=1,dft_control%nspins
             CALL cp_dbcsr_write_sparse_matrix(rho%rho_ao(ispin)%matrix,4,6,qs_env,&
                  para_env,output_unit=iw,error=error)
          END DO
          CALL cp_print_key_finished_output(iw,logger,input,&
               "DFT%PRINT%AO_MATRICES/DENSITY", error=error)
       END IF

       !    **** the Kohn-Sham matrix itself
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",error=error),cp_p_file)) THEN
          CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.FALSE., error=error)
          CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,error=error)
          iw = cp_print_key_unit_nr(logger,input,"DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",&
               extension=".Log",error=error)
          CALL cp_dbcsr_write_sparse_matrix(ks_rmpv(1)%matrix,4,6,Qs_env,para_env,output_unit=iw,error=error)
          CALL cp_print_key_finished_output(iw,logger,input,&
               "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX", error=error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE write_available_results

END MODULE qs_scf_post_se
