!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2005  CP2K developers group                            !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/harris_functional [1.0] *
!!
!!   NAME
!!     harris_functional
!!
!!   FUNCTION
!!     Performs the Harris functional energy correction
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE harris_functional
  USE pws,                             ONLY: cleanup_pw_poisson_solver, &
                                             init_pw_poisson_solver
  USE qs_rho_types,                    ONLY: qs_rho_type, &
                                             qs_rho_get
  USE qs_environment_types,            ONLY: qs_environment_type, &
                                             get_qs_env
  USE cp_log_handling,                 ONLY: cp_failure_level, &
                                             cp_note_level, &
                                             cp_to_string, &
                                             cp_warning_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type
  USE cp_error_handling,               ONLY: cp_a_l, &
                                             cp_assert, &
                                             cp_assertion_failed, &
                                             cp_debug, &
                                             cp_error_get_logger, &
                                             cp_error_message, &
                                             cp_error_type
  USE kinds,                           ONLY: dp
  USE qs_mo_types,                     ONLY: get_mo_set, &
                                             mo_set_p_type
  USE hartree,                         ONLY: calculate_hartree
  USE qs_ks_methods,                   ONLY: qs_vxc_create
                                             ! calc_rho_tot_gspace
  USE green_types,                     ONLY: greens_function_type
  USE pw_pool_types,                   ONLY: pw_pool_type, &
                                             pw_pool_init_coeff, &
                                             pw_pool_give_back_coeff
  USE pw_env_types,                    ONLY: pw_env_type, &
                                             pw_env_get
  USE pw_types,                        ONLY: COMPLEXDATA1D, &
                                             RECIPROCALSPACE
  USE coefficient_types,               ONLY: coeff_type, &
                                             coeff_add, &
                                             coeff_sumup
  !USE qs_collocate_density,            ONLY: calculate_total_rho
  !USE qs_charges_types,                ONLY: qs_charges_type
  USE pws,                             ONLY: pw_integral_ab
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE timings,                         ONLY: timeset, timestop
  !USE qs_core_energies,                ONLY: calculate_ecore, &
  !                                           calculate_ecore_overlap, &
  !                                           calculate_ecore_self
  USE harris_env_types,                ONLY: harris_env_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type, &
                                             cp_sm_sm_trace
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: moduleN = 'harris_functional'

  ! *** Public subroutines ***
  PUBLIC :: harris_energy_correction, &
            harris_eigenvalue_summation, &
            harris_eigenvalue_calculation, &
            harris_eigenvalue_trace_KS_Pmix, &
            harris_postprocessing

!***
!****************************************************************************

CONTAINS

!!****f* harris_functional/harris_energy_correction [1.0] *
!!
!!   NAME
!!     harris_energy_correction
!!
!!   FUNCTION
!!     Performs the harris functional energy correction
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rho: The frozen input density, in general the sum of 
!!            the atomic densities
!!     - qs_env: The QS environment, on which the Harris 
!!               functional energy correction should be performed.
!!     - harris_env: The harris environment in which all results are stored
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_energy_correction(qs_env, harris_env, EII_necessary, &
                                      error)

    TYPE(qs_environment_type), POINTER              :: qs_env
    TYPE(harris_env_type), POINTER                  :: harris_env
    LOGICAL, INTENT(INOUT), OPTIONAL                :: EII_necessary
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    TYPE(qs_rho_type), POINTER                      :: rho
    INTEGER                                         :: handle
    TYPE(greens_function_type), POINTER             :: green
    TYPE(pw_env_type), POINTER                      :: pw_env
    TYPE(pw_pool_type), POINTER                     :: auxbas_pw_pool
    TYPE(coeff_type)                                :: rho_tot_gspace, &
                                                       v_hartree_gspace
    INTEGER                                         :: ispin
    TYPE(coeff_type), POINTER                       :: rho_core
    !TYPE(qs_charges_type), POINTER                  :: qs_charges
    TYPE(coeff_type), DIMENSION(:), POINTER         :: v_rspace_new, &
                                                       v_tau_rspace
    !INTEGER, DIMENSION(:,:), POINTER                :: bo
    !INTEGER, DIMENSION(2,3)                         :: bo
    TYPE(qs_energy_type), POINTER                   :: qs_energy
    LOGICAL                                         :: failure, my_EII_necessary

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_energy_correction', &
      routineP = moduleN//':'//routineN
    TYPE(cp_logger_type), POINTER :: logger
    INTEGER :: unit_nr

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
    NULLIFY(rho, green, pw_env, auxbas_pw_pool, v_rspace_new, v_tau_rspace, &
            rho_core, qs_energy)
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    my_EII_necessary = .TRUE.
    IF (PRESENT(EII_necessary)) my_EII_necessary=EII_necessary

    IF (logger%para_env%source==logger%para_env%mepos) THEN
          unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.)
          WRITE (unit_nr,*) ""
          WRITE (unit_nr,*) "The Harris functional energy correction is performed!"
    END IF

    IF (.NOT. failure) THEN

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, error=error)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)
      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_gspace, &
                              use_data=COMPLEXDATA1D, &
                              in_space=RECIPROCALSPACE, error=error)
      CALL pw_pool_init_coeff(auxbas_pw_pool, rho_tot_gspace, &
                              use_data=COMPLEXDATA1D, &
                              in_space=RECIPROCALSPACE, error=error)

      ! ** CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho) ** !
      !CALL get_qs_env(qs_env=qs_env, rho_core=rho_core,qs_charges=qs_charges)
      CALL get_qs_env(qs_env=qs_env, rho_core=rho_core, error=error)
      CALL get_qs_env(qs_env, rho=rho, error=error)
      CALL coeff_add(rho_core, rho%rho_g(1), rho_tot_gspace)
      DO ispin=2, qs_env%dft_control%nspins
         CALL coeff_sumup(rho%rho_g(ispin), rho_tot_gspace)
      END DO
      !qs_charges%total_rho_gspace = calculate_total_rho(rho_tot_gspace)

      CALL get_qs_env(qs_env=qs_env, green=green, error=error)
      CALL init_pw_poisson_solver(green, error=error)
      !CALL calculate_hartree(rho_tot_gspace, ehartree=harris_env%Ehartree_tot, &
      !                       vhartree=v_hartree_gspace)
      CALL calculate_hartree(rho%rho_g(1), ehartree=harris_env%Ehartree_elec, &
                             vhartree=v_hartree_gspace)
      IF (my_EII_necessary) THEN
        CALL calculate_hartree(rho_core, ehartree=harris_env%Ehartree_core, &
                               vhartree=v_hartree_gspace)
      END IF
      CALL cleanup_pw_poisson_solver(green, error=error)

      CALL pw_pool_give_back_coeff(auxbas_pw_pool, rho_tot_gspace, &
                                   error=error)
      CALL pw_pool_give_back_coeff(auxbas_pw_pool, v_hartree_gspace, &
                                   error=error)

      !harris_env%Ehartree_elec_core = harris_env%Ehartree_tot &
      !                              - harris_env%Ehartree_elec &
      !                              - harris_env%Ehartree_core

      CALL qs_vxc_create(qs_env=qs_env, vxc_rho=v_rspace_new, &
                         vxc_tau = v_tau_rspace, exc=harris_env%Exc, &
                         just_energy=.FALSE., error=error)

      ! ** Evaluate integral_vxc ** !
      !sum=0.0_dp
      !DO ispin = 1, qs_env%dft_control%nspins
      !bo=v_rspace_new(ispin)%pw%pw_grid%bounds_local
      !DO k=bo(1,3),bo(2,3)
      !   DO j=bo(1,2),bo(2,2)
      !      DO i=bo(1,1),bo(2,1)
      !        sum = sum + (v_rspace_new(ispin)%pw%cr3d(i,j,k) * &
      !                     rho%rho_r(ispin)%pw%cr3d(i,j,k))
      !      END DO
      !   END DO
      !END DO
      !END DO
      !CALL mp_sum(sum)
      !harris_env%integral_vxc = sum * v_rspace_new(1)%pw%pw_grid%dvol

      ispin = 1
      harris_env%integral_vxc = pw_integral_ab(v_rspace_new(ispin)%pw, &
                                               rho%rho_r(ispin)%pw)

      IF (my_EII_necessary) THEN
        CALL get_qs_env(qs_env=qs_env, energy=qs_energy)
        harris_env%Ecore = qs_energy%core
        harris_env%Ecore_overlap = qs_energy%core_overlap
        harris_env%Ecore_self = qs_energy%core_self
        !CALL get_qs_env(qs_env, para_env=para_env, error=error)
        !CALL calculate_ecore(h, p, ecore, para_env)      
        !CALL calculate_ecore_overlap(qs_env, globenv, .FALSE.)
        !CALL calculate_ecore_self(qs_env)

        harris_env%EII = harris_env%Ecore_overlap + harris_env%Ecore_self + &
                         harris_env%Ehartree_core
        IF (PRESENT(EII_necessary)) EII_necessary=.FALSE.
        WRITE (cp_logger_get_default_unit_nr(logger,local=.TRUE.),*) "EII is calculated"
      ELSE
        WRITE (cp_logger_get_default_unit_nr(logger,local=.TRUE.),*) "EII is not calculated"
      END IF

      harris_env%Eharris_correction = - harris_env%Ehartree_elec &
                                      - harris_env%integral_vxc &
                                      + harris_env%Exc + harris_env%EII

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_energy_correction


!!****f* harris_functional/harris_eigenvalue_summation [1.0] *
!!
!!   NAME
!!     harris_eigenvalue_summation
!!
!!   FUNCTION
!!     Calculates the sum of eigenvalues needed for the Harris energy 
!!     functional out of the eigenvalues extracted from the MOs.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: The QS environment, on which the Harris 
!!               functional energy correction should be performed.
!!     - harris_env: The harris environment in which all results are stored
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_eigenvalue_summation(qs_env, harris_env, error)

    TYPE(qs_environment_type), POINTER              :: qs_env
    TYPE(harris_env_type), POINTER                  :: harris_env
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mo_array
    INTEGER                                         :: homo, &
                                                       iterator, &
                                                       handle
    REAL(KIND=dp), DIMENSION(:), POINTER            :: eigenvalues, &
                                                       occupation_numbers
    REAL(KIND=dp)                                   :: sum_of_eigenvalues
    LOGICAL                                         :: failure

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_eigenvalue_summation', &
      routineP = moduleN//':'//routineN

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
    NULLIFY(eigenvalues, occupation_numbers, mo_array)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    
    sum_of_eigenvalues = 0.0_dp
 
    IF (.NOT. failure) THEN

      CALL get_qs_env(qs_env, mos=mo_array, error=error)
      CALL get_mo_set(mo_array(1)%mo_set, homo=homo, eigenvalues=eigenvalues, &
                      occupation_numbers = occupation_numbers)
    
      DO iterator=1,homo
         sum_of_eigenvalues = sum_of_eigenvalues + occupation_numbers(iterator) &
                              * eigenvalues(iterator)
      END DO
    
      ! Write the sum of eigenvalues back to the harris_env_type
      harris_env%sum_of_eigenvalues = sum_of_eigenvalues

      harris_env%Eharris = harris_env%sum_of_eigenvalues + &
                           harris_env%Eharris_correction

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_eigenvalue_summation


!!****f* harris_functional/harris_eigenvalue_calculation [1.0] *
!!
!!   NAME
!!     harris_eigenvalue_calculation
!!
!!   FUNCTION
!!     Calculates the sum of eigenvalues needed for the Harris energy 
!!     functional as the trace of sparse matrix product K*P, thus:
!!     sum_of_eigenvalues = tr(K*P)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: The QS environment, on which the Harris 
!!               functional energy correction should be performed.
!!     - harris_env: The harris environment in which all results are stored
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_eigenvalue_calculation(qs_env, harris_env, error)

    TYPE(qs_environment_type), POINTER              :: qs_env
    TYPE(harris_env_type), POINTER                  :: harris_env
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    INTEGER                                         :: handle
    LOGICAL                                         :: failure
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: matrix_ks, &
                                                       rho_ao
    TYPE(qs_rho_type), POINTER                      :: rho
    TYPE(cp_para_env_type), POINTER                 :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_eigenvalue_calculation', &
      routineP = moduleN//':'//routineN

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
    NULLIFY(rho, rho_ao, matrix_ks)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    
    IF (.NOT. failure) THEN

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, error=error)
      CALL get_qs_env(qs_env, rho=rho, error=error)
      CALL qs_rho_get(rho_struct=rho, rho_ao=rho_ao, error=error)
           
      CALL get_qs_env(qs_env, para_env=para_env, error=error)
      CALL cp_sm_sm_trace(matrix_ks(1)%matrix, rho_ao(1)%matrix, &
                          harris_env%sum_of_eigenvalues, &
                          para_env=para_env, error=error)

      harris_env%Eharris = harris_env%sum_of_eigenvalues + &
                           harris_env%Eharris_correction

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_eigenvalue_calculation


!!****f* harris_functional/harris_eigenvalue_trace_KS_Pmix [1.0] *
!!
!!   NAME
!!     harris_eigenvalue_trace_KS_Pmix
!!
!!   FUNCTION
!!     Calculates the sum of eigenvalues needed for the Harris energy 
!!     functional as the trace of sparse matrix product KS*Pmix_new, thus:
!!     sum_of_eigenvalues = tr(KS*Pmix_new)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: The QS environment, on which the Harris 
!!               functional energy correction should be performed.
!!     - harris_env: The harris environment in which all results are stored
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_eigenvalue_trace_KS_Pmix(scf_env, qs_env, harris_env, error)

    TYPE(qs_scf_env_type), POINTER                  :: scf_env
    TYPE(qs_environment_type), POINTER              :: qs_env
    TYPE(harris_env_type), POINTER                  :: harris_env
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    INTEGER                                         :: handle
    LOGICAL                                         :: failure
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: matrix_ks
    TYPE(cp_para_env_type), POINTER                 :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_eigenvalue_trace_KS_Pmix', &
      routineP = moduleN//':'//routineN

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
    NULLIFY(matrix_ks)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level, routineP, error, failure)
    CPPrecondition(scf_env%ref_count>0, CP_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    
    IF (.NOT. failure) THEN

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, error=error)
                      
      CALL get_qs_env(qs_env, para_env=para_env, error=error)
      CALL cp_sm_sm_trace(matrix_ks(1)%matrix, scf_env%p_mix_new(1)%matrix, &
                          harris_env%sum_of_eigenvalues, &
                          para_env=para_env, error=error)

      harris_env%Eharris = harris_env%sum_of_eigenvalues + &
                           harris_env%Eharris_correction

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_eigenvalue_trace_KS_Pmix


!!****f* harris_functional/harris_postprocessing [1.0] *
!!
!!   NAME
!!     harris_postprocessing
!!
!!   FUNCTION
!!     Performs the output and diverse analyses
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - harris_env: The harris environment in which all results are stored
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_postprocessing(harris_env, error)

    TYPE(harris_env_type), POINTER                  :: harris_env
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    LOGICAL                                         :: failure

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_postprocessing', &
         routineP = moduleN//':'//routineN
    TYPE(cp_logger_type), POINTER :: logger
    INTEGER :: unit_nr

    !   ------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    logger => cp_error_get_logger(error)

    IF (.NOT. failure) THEN

       ! Output
       IF (logger%para_env%source==logger%para_env%mepos) THEN
          unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.)
          WRITE (unit_nr,*) ""
          WRITE (unit_nr,*) "Ehartree n_elec            =", harris_env%Ehartree_elec
          WRITE (unit_nr,*) "Ehartree n_core            =", harris_env%Ehartree_core
          !write (unit_nr,*) "Ehartree n_tot             =", harris_env%Ehartree_tot
          !write (unit_nr,*) "Ehartree_elec_core         =", harris_env%Ehartree_elec_core

          write (unit_nr,*) "Exc                        =", harris_env%Exc

          write (unit_nr,*) "The XC potential integral  =", harris_env%integral_vxc

          write (unit_nr,*) "Ecore                      =", harris_env%Ecore
          write (unit_nr,*) "Ecore_overlap              =", harris_env%Ecore_overlap
          write (unit_nr,*) "Ecore_self                 =", harris_env%Ecore_self

          write (unit_nr,*) "EII                        =", harris_env%EII

          write (unit_nr,*) "Eharris correction energy  =", harris_env%Eharris_correction

          write (unit_nr,*) "The sum of the eigenvalues =", harris_env%sum_of_eigenvalues
          write (unit_nr,*) "Eharris                    =", harris_env%Eharris
       END IF
    END IF

  END SUBROUTINE harris_postprocessing

END MODULE
