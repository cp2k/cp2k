!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief g tensor calculation by dfpt
!>      Initialization of the epr_env, creation of the special neighbor lists
!>      Perturbation Hamiltonians by application of the p and rxp oprtators to  psi0
!>      Write output
!>      Deallocate everything
!> \note
!>      The psi0 should be localized
!>      the Sebastiani method works within the assumption that the orbitals are
!>      completely contained in the simulation box
!> \par History
!>       created 07-2005 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_epr_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             twopi
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: a_fine,&
                                             e_gfactor
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_release
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_types,                 ONLY: deallocate_nablavks_atom_set,&
                                             epr_env_create,&
                                             epr_env_type,&
                                             init_nablavks_atom_set,&
                                             linres_control_type,&
                                             nablavks_atom_type,&
                                             set_epr_env
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE qs_rho_atom_types,               ONLY: deallocate_rho_atom_set
  USE scf_control_types,               ONLY: scf_control_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: epr_env_cleanup, epr_env_init

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_epr_utils'

CONTAINS

! *****************************************************************************
!> \brief Initialize the epr environment
!> \param epr_env ...
!> \param qs_env ...
!> \param error ...
!> \par History
!>      07.2006 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE epr_env_init(epr_env,qs_env,error)
    !
    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_B, idir, ispin, &
                                                istat, n_mo(2), nao, natom, &
                                                nmoloc, nspins, output_unit
    LOGICAL                                  :: failure, gapw
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(nablavks_atom_type), DIMENSION(:), &
      POINTER                                :: nablavks_atom_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: lr_section

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(atomic_kind_set, cell, dft_control, linres_control, scf_control)
    NULLIFY(logger, mos, mpools, particle_set)
    NULLIFY(auxbas_pw_pool, pw_env)
    NULLIFY(nablavks_atom_set)

    n_mo(1:2) = 0
    nao = 0
    nmoloc = 0

    logger => cp_error_get_logger(error)
    !ionode = logger%para_env%mepos==logger%para_env%source
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
                                       extension=".linresLog",error=error)

    IF(epr_env%ref_count /= 0) THEN
      CALL epr_env_cleanup(epr_env,error=error)
    END IF

    IF(output_unit>0) THEN
      WRITE(output_unit,"(/,T20,A,/)") "*** Start EPR g tensor calculation ***"
      WRITE(output_unit,"(T10,A,/)") "Initialization of the EPR environment"
    ENDIF

    CALL epr_env_create(epr_env,error=error)


    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    linres_control=linres_control,&
                    mos=mos,&
                    mpools=mpools,&
                    particle_set=particle_set,&
                    pw_env=pw_env,&
                    scf_control=scf_control,error=error)
    !
    ! Check if restat also psi0 should be restarted
    !IF(epr_env%restart_epr .AND. scf_control%density_guess/=restart_guess)THEN
    !   CALL stop_program(routineN,moduleN,__LINE__,"restart_epr requires density_guess=restart")
    !ENDIF
    !
    ! check that the psi0 are localized and you have all the centers
    CPPrecondition(linres_control%localized_psi0,cp_warning_level,routineP,error,failure)
    IF(failure .AND. (output_unit>0)) THEN
       WRITE(output_unit,'(A)') &
            ' To get EPR parameters within PBC you need localized zero order orbitals '
    ENDIF
    gapw = dft_control%qs_control%gapw
    nspins = dft_control%nspins
    natom = SIZE(particle_set,1)
    !
    ! Conversion factors
    ! Magical constant twopi/cell%deth just like in NMR shift (basically undo scale_fac in qs_linres_nmr_current.F)
    epr_env%g_free_factor = -1.0_dp * e_gfactor
    epr_env%g_zke_factor = e_gfactor * ( a_fine )**2
    epr_env%g_so_factor = ( a_fine )**2 * ( -1.0_dp * e_gfactor - 1.0_dp ) / 2.0_dp * twopi / cell%deth
    epr_env%g_so_factor_gapw =  ( a_fine )**2 * ( -1.0_dp * e_gfactor - 1.0_dp ) / 2.0_dp
    ! * 2 because B_ind = 2 * B_beta
    epr_env%g_soo_factor = 2.0_dp * fourpi * ( a_fine )**2 * twopi / cell%deth
    ! 2 * 2 * 1/4 * e^2 / m * a_0^2 * 2/3 * mu_0 / (omega * 1e-30 )
    epr_env%g_soo_chicorr_factor = 2.0 / 3.0_dp * fourpi * ( a_fine )**2 / cell%deth
    !
    ! If the current density on the grid needs to be stored
    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
    !
    ! Initialize local current density if GAPW calculation
    IF(gapw) THEN
       CALL init_nablavks_atom_set(nablavks_atom_set,atomic_kind_set,&
                   nspins,error=error)
       CALL set_epr_env(epr_env=epr_env,&
                        nablavks_atom_set=nablavks_atom_set,&
                        error=error)
    ENDIF
    !
    ! Bind
    ALLOCATE(epr_env%bind_set(3,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i_B = 1,3
      DO idir = 1,3
        NULLIFY(epr_env%bind_set(idir,i_B)%rho)
        ALLOCATE(epr_env%bind_set(idir,i_B)%rho,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        epr_env%bind_set(idir,i_B)%rho%ref_count = 1
        NULLIFY( epr_env%bind_set(idir,i_B)%rho%rho_r)
        NULLIFY( epr_env%bind_set(idir,i_B)%rho%rho_g)

        ALLOCATE(epr_env%bind_set(idir,i_B)%rho%rho_r(1),stat=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(epr_env%bind_set(idir,i_B)%rho%rho_g(1),stat=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
             epr_env%bind_set(idir,i_B)%rho%rho_r(1)%pw,&
             use_data=REALDATA3D,in_space=REALSPACE,error=error)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
             epr_env%bind_set(idir,i_B)%rho%rho_g(1)%pw,&
             use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
      END DO
    END DO

    ! Nabla_V_ks
    ALLOCATE(epr_env%nablavks_set(3,dft_control%nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir = 1,3
      DO ispin = 1,nspins
        NULLIFY(epr_env%nablavks_set(idir,ispin)%rho)
        ALLOCATE(epr_env%nablavks_set(idir,ispin)%rho,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        epr_env%nablavks_set(idir,ispin)%rho%ref_count = 1
        NULLIFY( epr_env%nablavks_set(idir,ispin)%rho%rho_r)
        NULLIFY( epr_env%nablavks_set(idir,ispin)%rho%rho_g)
        ALLOCATE(epr_env%nablavks_set(idir,ispin)%rho%rho_r(1),stat=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(epr_env%nablavks_set(idir,ispin)%rho%rho_g(1),stat=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
             epr_env%nablavks_set(idir,ispin)%rho%rho_r(1)%pw,&
             use_data=REALDATA3D,in_space=REALSPACE,error=error)
        CALL pw_pool_create_pw(auxbas_pw_pool,&
             epr_env%nablavks_set(idir,ispin)%rho%rho_g(1)%pw,&
             use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
      END DO
    END DO

    ! Initialize the g tensor components
    ALLOCATE(epr_env%g_total(3,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(epr_env%g_so(3,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(epr_env%g_soo(3,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    epr_env%g_total = 0.0_dp
    epr_env%g_zke = 0.0_dp
    epr_env%g_so = 0.0_dp
    epr_env%g_soo = 0.0_dp

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         &                            "PRINT%PROGRAM_RUN_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE epr_env_init

! *****************************************************************************
!> \brief Deallocate the epr environment
!> \param epr_env ...
!> \param error ...
!> \par History
!>      07.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE epr_env_cleanup(epr_env,error)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_env_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i_B, idir, ispin, istat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF(.NOT. failure) THEN
      epr_env%ref_count = epr_env%ref_count - 1
      IF(epr_env%ref_count == 0 ) THEN
        ! nablavks_set
        IF(ASSOCIATED(epr_env%nablavks_set)) THEN
           DO ispin = 1,SIZE(epr_env%nablavks_set,2)
           DO idir = 1,SIZE(epr_env%nablavks_set,1)
              DO i=1,SIZE(epr_env%nablavks_set(idir,ispin)%rho%rho_r)
                 CALL pw_release(epr_env%nablavks_set(idir,ispin)%rho%rho_r(i)%pw,error=error)
              END DO
              DO i=1,SIZE(epr_env%nablavks_set(idir,ispin)%rho%rho_g)
                 CALL pw_release(epr_env%nablavks_set(idir,ispin)%rho%rho_g(i)%pw,error=error)
              END DO
              DEALLOCATE(epr_env%nablavks_set(idir,ispin)%rho%rho_g,stat=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
              DEALLOCATE(epr_env%nablavks_set(idir,ispin)%rho%rho_r,stat=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
              DEALLOCATE(epr_env%nablavks_set(idir,ispin)%rho,STAT=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
           ENDDO
           ENDDO
           DEALLOCATE(epr_env%nablavks_set,STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! nablavks_atom_set
        IF(ASSOCIATED(epr_env%nablavks_atom_set)) THEN
          CALL deallocate_nablavks_atom_set(epr_env%nablavks_atom_set,error=error)
        END IF
        ! vks_atom_set
        IF(ASSOCIATED(epr_env%vks_atom_set)) THEN
          CALL deallocate_rho_atom_set(epr_env%vks_atom_set)
        END IF
        ! bind_set
        IF(ASSOCIATED(epr_env%bind_set)) THEN
           DO i_B = 1,SIZE(epr_env%bind_set,2)
           DO idir = 1,SIZE(epr_env%bind_set,1)
              DO i=1,SIZE(epr_env%bind_set(idir,i_B)%rho%rho_r)
                 CALL pw_release(epr_env%bind_set(idir,i_B)%rho%rho_r(i)%pw,error=error)
              END DO
              DO i=1,SIZE(epr_env%bind_set(idir,i_B)%rho%rho_g)
                 CALL pw_release(epr_env%bind_set(idir,i_B)%rho%rho_g(i)%pw,error=error)
              END DO
              DEALLOCATE(epr_env%bind_set(idir,i_B)%rho%rho_r,STAT=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
              DEALLOCATE(epr_env%bind_set(idir,i_B)%rho%rho_g,STAT=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
              DEALLOCATE(epr_env%bind_set(idir,i_B)%rho,STAT=istat)
              CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
           ENDDO
           ENDDO
          DEALLOCATE(epr_env%bind_set,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! bind_atom_set
        IF(ASSOCIATED(epr_env%bind_atom_set)) THEN
          DEALLOCATE(epr_env%bind_atom_set,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! g_total
        IF(ASSOCIATED(epr_env%g_total)) THEN
          DEALLOCATE(epr_env%g_total,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! g_so
        IF(ASSOCIATED(epr_env%g_so)) THEN
          DEALLOCATE(epr_env%g_so,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
        ! g_soo
        IF(ASSOCIATED(epr_env%g_soo)) THEN
          DEALLOCATE(epr_env%g_soo,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF
      END IF  ! ref count
    END IF ! failure

  END SUBROUTINE epr_env_cleanup

END MODULE qs_linres_epr_utils
