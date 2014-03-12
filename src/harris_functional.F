!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Performs the Harris functional energy correction
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
MODULE harris_functional
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_trace
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE harris_energy_types,             ONLY: harris_energy_type
  USE harris_env_types,                ONLY: harris_env_get,&
                                             harris_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_integral_ab
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_methods,                  ONLY: duplicate_rho_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE qs_vxc,                          ONLY: qs_vxc_create
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'harris_functional'

  ! *** Public subroutines ***
  PUBLIC :: harris_energy_correction, &
            harris_eigenvalue_summation, &
            harris_eigenvalue_calculation, &
            harris_eigenvalue_trace_KS_Pmix, &
            harris_postprocessing

!***

CONTAINS

! *****************************************************************************
!> \brief Performs the harris functional energy correction
!> \param rho The frozen input density, in general the sum of
!>             the atomic densities
!> \param qs_env The QS environment, on which the Harris
!>                functional energy correction should be performed.
!> \param harris_env The harris environment in which all results are stored
!> \param EII_necessary Logical flag which indicates if EII has to be
!>                       computes or not
!> \param fast Logical which prevents some minor double calculations if the
!>              Harris functional is driven on top of the QS SCF-cycle.
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_energy_correction(qs_env, harris_env, para_env, EII_necessary, fast, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN), OPTIONAL            :: EII_necessary, fast
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_energy_correction', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins, stat
    LOGICAL                                  :: failure, fast_flag, &
                                                my_EII_necessary
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: v_hartree_gspace
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_rspace_new, v_tau_rspace
    TYPE(pw_p_type), POINTER                 :: rho_core
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: input, xc_section

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(rho, pw_env, auxbas_pw_pool, v_rspace_new, v_tau_rspace, &
            rho_core, rho, harris_energy, poisson_env, dft_control, input, xc_section)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)

    my_EII_necessary = .TRUE.
    IF (PRESENT(EII_necessary)) my_EII_necessary=EII_necessary

    IF (PRESENT(fast)) THEN
      fast_flag = fast
    ELSE
      fast_flag = .FALSE.
    END IF

    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_energy=harris_energy, error=error)
      CALL get_qs_env(qs_env=qs_env, input=input, pw_env=pw_env, rho=rho, rho_core=rho_core, &
                      dft_control=dft_control, error=error)

      xc_section => section_vals_get_subs_vals(input, "DFT%XC", error=error)

      nspins = dft_control%nspins

      CALL duplicate_rho_type(rho_input=rho, rho_output=harris_env%rho, &
                              qs_env=qs_env, error=error)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env, &
                      error=error)

      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                              use_data=COMPLEXDATA1D, &
                              in_space=RECIPROCALSPACE, error=error)

      DO ispin = 1,nspins
        CALL pw_poisson_solve(poisson_env, rho%rho_g(ispin)%pw, &
                              ehartree=harris_energy%Ehartree_elec, &
                              vhartree=v_hartree_gspace%pw,error=error)
      END DO

      IF (my_EII_necessary) THEN
         CALL pw_poisson_solve(poisson_env, rho_core%pw, ehartree=harris_energy%Ehartree_core, &
                               vhartree=v_hartree_gspace%pw,error=error)
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                                   error=error)

      CALL qs_vxc_create(qs_env=qs_env, rho_struct=harris_env%rho, xc_section=xc_section, &
                         vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=harris_energy%Exc, &
                         just_energy=.FALSE., error=error)

      DO ispin = 1,nspins
        harris_energy%integral_vxc = pw_integral_ab(v_rspace_new(ispin)%pw, &
                                                    rho%rho_r(ispin)%pw,error=error)
      END DO

      IF (ASSOCIATED(v_rspace_new)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_new(ispin)%pw, error=error)
        END DO
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw, error=error)
        END DO
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF

      IF (my_EII_necessary) THEN
        IF (.NOT. fast_flag) THEN
          qs_env%energy%core_overlap = 0.0_dp
          qs_env%energy%core_self = 0.0_dp

          CALL calculate_ecore_overlap(qs_env=qs_env, para_env=para_env, &
                                       calculate_forces=.FALSE., &
                                       E_overlap_core=harris_energy%Ecore_overlap,&
                                       error=error)
          CALL calculate_ecore_self(qs_env, E_self_core=harris_energy%Ecore_self,error=error)
        ELSE
          harris_energy%Ecore_overlap = qs_env%energy%core_overlap
          harris_energy%Ecore_self = qs_env%energy%core_self
        END IF

        harris_energy%EII = harris_energy%Ecore_overlap + harris_energy%Ecore_self + &
                            harris_energy%Ehartree_core
        !IF (PRESENT(EII_necessary)) EII_necessary=.FALSE.
      END IF

      harris_energy%Eharris_correction = - harris_energy%Ehartree_elec &
                                         - harris_energy%integral_vxc &
                                         + harris_energy%Exc + harris_energy%EII
    END IF

    CALL timestop(handle)
  END SUBROUTINE harris_energy_correction

! *****************************************************************************
!> \brief Calculates the sum of eigenvalues needed for the Harris energy
!>      functional out of the eigenvalues extracted from the MOs.
!> \param qs_env The QS environment, on which the Harris
!>                functional energy correction should be performed.
!> \param harris_env The harris environment in which all results are stored
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_eigenvalue_summation(qs_env, harris_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_eigenvalue_summation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, ispin, &
                                                iterator, nspins
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: sum_of_eigenvalues
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues, &
                                                occupation_numbers
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array

!   ------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(eigenvalues, occupation_numbers, mo_array, harris_energy, dft_control)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)

    sum_of_eigenvalues = 0.0_dp

    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_energy=harris_energy, error=error)

      CALL get_qs_env(qs_env, dft_control=dft_control, mos=mo_array, error=error)

      nspins = dft_control%nspins
      DO ispin = 1,nspins
        CALL get_mo_set(mo_array(ispin)%mo_set, homo=homo, eigenvalues=eigenvalues, &
                        occupation_numbers = occupation_numbers)
      END DO

      DO iterator=1,homo
         sum_of_eigenvalues = sum_of_eigenvalues + occupation_numbers(iterator) &
                              * eigenvalues(iterator)
      END DO

      ! Write the sum of eigenvalues back to the harris_energy_type
      harris_energy%sum_of_eigenvalues = sum_of_eigenvalues

      harris_energy%Eharris = harris_energy%sum_of_eigenvalues + &
                              harris_energy%Eharris_correction

    END IF

    CALL timestop(handle)

  END SUBROUTINE harris_eigenvalue_summation

! *****************************************************************************
!> \brief Calculates the sum of eigenvalues needed for the Harris energy
!>      functional as the trace of sparse matrix product K*P, thus:
!>      sum_of_eigenvalues = tr(K*P)
!> \param qs_env The QS environment, on which the Harris
!>                functional energy correction should be performed.
!> \param harris_env The harris environment in which all results are stored
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_eigenvalue_calculation(qs_env, harris_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'harris_eigenvalue_calculation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, rho_ao
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(qs_rho_type), POINTER               :: rho

!   ------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(rho, rho_ao, matrix_ks, harris_energy, dft_control)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_energy=harris_energy, error=error)

      CALL get_qs_env(qs_env, rho=rho, dft_control=dft_control, matrix_ks=matrix_ks, &
                      para_env=para_env, error=error)
      CALL qs_rho_get(rho_struct=rho, rho_ao=rho_ao, error=error)

      nspins = dft_control%nspins
      DO ispin = 1,nspins
        CALL cp_dbcsr_trace(matrix_ks(ispin)%matrix, rho_ao(ispin)%matrix, &
                         trace=harris_energy%sum_of_eigenvalues, error=error)
      END DO

      harris_energy%Eharris = harris_energy%sum_of_eigenvalues &
                            + harris_energy%Eharris_correction
    END IF

    CALL timestop(handle)

  END SUBROUTINE harris_eigenvalue_calculation

! *****************************************************************************
!> \brief Calculates the sum of eigenvalues needed for the Harris energy
!>      functional as the trace of sparse matrix product KS*Pmix_new, thus:
!>      sum_of_eigenvalues = tr(KS*Pmix_new)
!> \param scf_env The SCF environment in which the Pmix matrix is stored
!> \param qs_env The QS environment, on which the Harris
!>                functional energy correction should be performed.
!> \param harris_env The harris environment in which all results are stored
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_eigenvalue_trace_KS_Pmix(scf_env, qs_env, harris_env, error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'harris_eigenvalue_trace_KS_Pmix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_energy_type), POINTER        :: harris_energy

!   ------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(matrix_ks, harris_energy, dft_control)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level, routineP, error, failure)
    CPPrecondition(scf_env%ref_count>0, CP_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_energy=harris_energy, error=error)

      CALL get_qs_env(qs_env, dft_control=dft_control, matrix_ks=matrix_ks, &
                      para_env=para_env, error=error)

      nspins = dft_control%nspins
      DO ispin = 1,nspins
        CALL cp_dbcsr_trace(matrix_ks(ispin)%matrix, scf_env%p_mix_new(ispin)%matrix, &
                            trace=harris_energy%sum_of_eigenvalues, &
                            error=error)
      END DO

      harris_energy%Eharris = harris_energy%sum_of_eigenvalues + &
                              harris_energy%Eharris_correction
    END IF

    CALL timestop(handle)

  END SUBROUTINE harris_eigenvalue_trace_KS_Pmix

! *****************************************************************************
!> \brief Performs the output and diverse analyses
!> \param harris_env The harris environment in which all results are stored
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_postprocessing(harris_env, error)

    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_postprocessing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(harris_energy_type), POINTER        :: harris_energy

!   ------------------------------------------------------------------------

    failure = .FALSE.
    NULLIFY(harris_energy)

    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    logger => cp_error_get_logger(error)

    IF (.NOT. failure) THEN
       CALL harris_env_get(harris_env=harris_env, harris_energy=harris_energy, error=error)

       ! Output
       unit_nr=cp_logger_get_default_io_unit(logger)
       IF (unit_nr>0) THEN
          WRITE (unit_nr,*) ""
          WRITE (unit_nr,*) "The Harris functional energy correction is performed!"
          WRITE (unit_nr,*) ""
          WRITE (unit_nr,*) "Ehartree n_elec            =", harris_energy%Ehartree_elec
          WRITE (unit_nr,*) "Ehartree n_core            =", harris_energy%Ehartree_core

          WRITE (unit_nr,*) "Exc                        =", harris_energy%Exc

          WRITE (unit_nr,*) "The XC potential integral  =", harris_energy%integral_vxc

          WRITE (unit_nr,*) "Ecore_overlap              =", harris_energy%Ecore_overlap
          WRITE (unit_nr,*) "Ecore_self                 =", harris_energy%Ecore_self

          WRITE (unit_nr,*) "EII                        =", harris_energy%EII

          WRITE (unit_nr,*) "Eharris correction energy  =", harris_energy%Eharris_correction

          WRITE (unit_nr,*) "The sum of the eigenvalues =", harris_energy%sum_of_eigenvalues
          WRITE (unit_nr,*) "Eharris                    =", harris_energy%Eharris
          WRITE (unit_nr,*) ""
       END IF
    END IF

  END SUBROUTINE harris_postprocessing

END MODULE harris_functional
