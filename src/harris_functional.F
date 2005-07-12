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
  USE coefficient_types,               ONLY: coeff_type
                                             !coeff_add, &
                                             !coeff_sumup, &
                                             !coeff_allocate, &
                                             !coeff_copy
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l, &
                                             cp_assert, &
                                             cp_assertion_failed, &
                                             cp_debug, &
                                             cp_error_get_logger, &
                                             cp_error_message, &
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level, &
                                             cp_note_level, &
                                             cp_to_string, &
                                             cp_warning_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE green_types,                     ONLY: greens_function_type
  USE harris_energy_types,             ONLY: harris_energy_type
  USE harris_env_types,                ONLY: harris_env_type, &
                                             harris_env_get
  USE hartree,                         ONLY: calculate_hartree
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_type, &
                                             pw_env_get
  USE pw_pool_types,                   ONLY: pw_pool_type, &
                                             pw_pool_init_coeff, &
                                             pw_pool_give_back_coeff
  USE pws,                             ONLY: cleanup_pw_poisson_solver, &
                                             init_pw_poisson_solver, &
                                             pw_integral_ab
  USE pw_types,                        ONLY: COMPLEXDATA1D, &
                                             RECIPROCALSPACE
  !USE qs_charges_types,                ONLY: qs_charges_type
  !USE qs_collocate_density,            ONLY: calculate_total_rho
  USE qs_core_energies,                ONLY: calculate_ecore_overlap, &
                                             calculate_ecore_self
                                             !calculate_ecore
  !USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: qs_environment_type, &
                                             get_qs_env
  USE qs_ks_methods,                   ONLY: qs_vxc_create
                                             !calc_rho_tot_gspace
  USE qs_mo_types,                     ONLY: get_mo_set, &
                                             mo_set_p_type
  USE qs_rho_methods,                  ONLY: duplicate_rho_type
  USE qs_rho_types,                    ONLY: qs_rho_type, &
                                             qs_rho_get
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type, &
                                             !replicate_matrix, &
                                             cp_sm_sm_trace
  USE timings,                         ONLY: timeset, timestop
  
 
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
!!     - EII_necessary: Logical flag which indicates if EII has to be 
!!                      computes or not
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
                                      globenv, error)

    TYPE(qs_environment_type), POINTER              :: qs_env
    TYPE(harris_energy_type), POINTER               :: harris_energy
    TYPE(harris_env_type), POINTER                  :: harris_env
    LOGICAL, INTENT(INOUT), OPTIONAL                :: EII_necessary
    TYPE(global_environment_type), POINTER          :: globenv
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    TYPE(qs_rho_type), POINTER                      :: rho
    INTEGER                                         :: handle, &
                                                       stat
    TYPE(greens_function_type), POINTER             :: green
    TYPE(pw_env_type), POINTER                      :: pw_env
    TYPE(pw_pool_type), POINTER                     :: auxbas_pw_pool
    TYPE(coeff_type)                                :: v_hartree_gspace
                                                       !rho_tot_gspace, &
    !INTEGER                                         :: ispin
    TYPE(coeff_type), POINTER                       :: rho_core
    !TYPE(qs_charges_type), POINTER                  :: qs_charges
    TYPE(coeff_type), DIMENSION(:), POINTER         :: v_rspace_new, &
                                                       v_tau_rspace
    !TYPE(qs_energy_type), POINTER                   :: qs_energy
    TYPE(cp_para_env_type), POINTER                 :: para_env
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: matrix_ks
    TYPE(dft_control_type), POINTER                 :: dft_control
    LOGICAL                                         :: failure, my_EII_necessary

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_energy_correction', &
      routineP = moduleN//':'//routineN
    !TYPE(cp_logger_type), POINTER :: logger
    !INTEGER :: unit_nr

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
    NULLIFY(rho, green, pw_env, auxbas_pw_pool, v_rspace_new, v_tau_rspace, &
            rho_core, para_env, matrix_ks, rho, harris_energy)
    !logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)

    my_EII_necessary = .TRUE.
    IF (PRESENT(EII_necessary)) my_EII_necessary=EII_necessary

    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_energy=harris_energy, error=error)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, rho=rho, rho_core=rho_core, &
                      green=green, error=error)

      !CALL coeff_allocate(c_out=harris_env%rho_core, c_model=rho_core, &
      !                    use_basis=rho_core%in_use)
      !CALL coeff_copy(rho_core, harris_env%rho_core)
      CALL duplicate_rho_type(rho_input=rho, rho_output=harris_env%rho, &
                              qs_env=qs_env, error=error)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)
      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_gspace, &
                              use_data=COMPLEXDATA1D, &
                              in_space=RECIPROCALSPACE, error=error)
      !CALL pw_pool_init_coeff(auxbas_pw_pool, rho_tot_gspace, &
      !                        use_data=COMPLEXDATA1D, &
      !                        in_space=RECIPROCALSPACE, error=error)

      ! ** CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho) ** !
      !CALL get_qs_env(qs_env=qs_env,rho=rho,rho_core=rho_core,qs_charges=qs_charges,error=error)

      !CALL coeff_add(rho_core, rho%rho_g(1), rho_tot_gspace)
      !DO ispin=2, qs_env%dft_control%nspins
      !   CALL coeff_sumup(rho%rho_g(ispin), rho_tot_gspace)
      !END DO
      !qs_charges%total_rho_gspace = calculate_total_rho(rho_tot_gspace)

      CALL init_pw_poisson_solver(green, error=error)
      !CALL calculate_hartree(rho_tot_gspace, ehartree=harris_energy%Ehartree_tot, &
      !                       vhartree=v_hartree_gspace)
      CALL calculate_hartree(rho%rho_g(1), ehartree=harris_energy%Ehartree_elec, &
                             vhartree=v_hartree_gspace)
      IF (my_EII_necessary) THEN
        CALL calculate_hartree(rho_core, ehartree=harris_energy%Ehartree_core, &
                               vhartree=v_hartree_gspace)
      END IF
      CALL cleanup_pw_poisson_solver(green, error=error)

      !CALL pw_pool_give_back_coeff(auxbas_pw_pool, rho_tot_gspace, &
      !                             error=error)
      CALL pw_pool_give_back_coeff(auxbas_pw_pool, v_hartree_gspace, &
                                   error=error)

      !harris_energy%Ehartree_elec_core = harris_energy%Ehartree_tot &
      !                                 - harris_energy%Ehartree_elec &
      !                                 - harris_energy%Ehartree_core

      CALL qs_vxc_create(qs_env=qs_env, vxc_rho=v_rspace_new, &
                         vxc_tau = v_tau_rspace, exc=harris_energy%Exc, &
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
      !harris_energy%integral_vxc = sum * v_rspace_new(1)%pw%pw_grid%dvol

      harris_energy%integral_vxc = pw_integral_ab(v_rspace_new(1)%pw, &
                                                  rho%rho_r(1)%pw)

      IF (ASSOCIATED(v_rspace_new)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(1),error=error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_tau_rspace(1),error=error)
      END IF

      IF (my_EII_necessary) THEN
        !CALL get_qs_env(qs_env=qs_env, energy=qs_energy, error=error)
        !harris_energy%Ecore = qs_energy%core
        !harris_energy%Ecore_overlap = qs_energy%core_overlap
        !harris_energy%Ecore_self = qs_energy%core_self
        CALL get_qs_env(qs_env, matrix_ks=matrix_ks, rho=rho, &
                        para_env=para_env, dft_control=dft_control, error=error)
        !CALL calculate_ecore(h=matrix_ks(1)%matrix, p=rho%rho_ao(1)%matrix, &
        !                     ecore=harris_energy%Ecore, para_env=para_env)
        CALL calculate_ecore_overlap(qs_env=qs_env, globenv=globenv, &
                                     calculate_forces=dft_control%forces, &
                                     E_overlap_core=harris_energy%Ecore_overlap)
        !harris_energy%Ecore_overlap = qs_energy%core_overlap
        CALL calculate_ecore_self(qs_env, E_self_core=harris_energy%Ecore_self)
        !harris_energy%Ecore_self = qs_energy%core_self

        harris_energy%EII = harris_energy%Ecore_overlap + harris_energy%Ecore_self + &
                            harris_energy%Ehartree_core
        IF (PRESENT(EII_necessary)) EII_necessary=.FALSE.
      END IF

      harris_energy%Eharris_correction = - harris_energy%Ehartree_elec &
                                         - harris_energy%integral_vxc &
                                         + harris_energy%Exc + harris_energy%EII

      IF (ASSOCIATED(v_rspace_new)) THEN
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      
      IF (ASSOCIATED(v_tau_rspace)) THEN
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_energy_correction
!***************************************************************************

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
    TYPE(harris_energy_type), POINTER               :: harris_energy
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
    NULLIFY(eigenvalues, occupation_numbers, mo_array, harris_energy)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    
    sum_of_eigenvalues = 0.0_dp
 
    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_energy=harris_energy, error=error)

      CALL get_qs_env(qs_env, mos=mo_array, error=error)
      CALL get_mo_set(mo_array(1)%mo_set, homo=homo, eigenvalues=eigenvalues, &
                      occupation_numbers = occupation_numbers)
    
      DO iterator=1,homo
         sum_of_eigenvalues = sum_of_eigenvalues + occupation_numbers(iterator) &
                              * eigenvalues(iterator)
      END DO
    
      ! Write the sum of eigenvalues back to the harris_energy_type
      harris_energy%sum_of_eigenvalues = sum_of_eigenvalues

      harris_energy%Eharris = harris_energy%sum_of_eigenvalues + &
                              harris_energy%Eharris_correction

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_eigenvalue_summation
!***************************************************************************

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
    TYPE(harris_energy_type), POINTER               :: harris_energy
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
    NULLIFY(rho, rho_ao, matrix_ks, harris_energy)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    
    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_energy=harris_energy, error=error)

      CALL get_qs_env(qs_env, rho=rho, matrix_ks=matrix_ks, para_env=para_env,  error=error)
      CALL qs_rho_get(rho_struct=rho, rho_ao=rho_ao, error=error)

      CALL cp_sm_sm_trace(matrix_ks(1)%matrix, rho_ao(1)%matrix, &
                          harris_energy%sum_of_eigenvalues, &
                          para_env=para_env, error=error)

      harris_energy%Eharris = harris_energy%sum_of_eigenvalues + &
                              harris_energy%Eharris_correction

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_eigenvalue_calculation
!***************************************************************************

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
!!     - scf_env: The SCF environment in which the Pmix matrix is stored
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
    TYPE(harris_energy_type), POINTER               :: harris_energy
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
    NULLIFY(matrix_ks, harris_energy)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level, routineP, error, failure)
    CPPrecondition(scf_env%ref_count>0, CP_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    
    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_energy=harris_energy, error=error)

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, para_env=para_env, error=error)
      !CALL replicate_matrix(scf_env%p_mix_new(1)%matrix, &
      !                      harris_env%p_mix_new(1)%matrix, &
      !                      target_name= 'p_mix_new_matrix',&
      !                      allocate_blocks=.TRUE.)
      !harris_energy%p_mix_new => scf_env%p_mix_new

      CALL cp_sm_sm_trace(matrix_ks(1)%matrix, scf_env%p_mix_new(1)%matrix, &
                          harris_energy%sum_of_eigenvalues, &
                          para_env=para_env, error=error)

      harris_energy%Eharris = harris_energy%sum_of_eigenvalues + &
                              harris_energy%Eharris_correction

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_eigenvalue_trace_KS_Pmix
!***************************************************************************

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

    TYPE(harris_energy_type), POINTER               :: harris_energy
    TYPE(harris_env_type), POINTER                  :: harris_env
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    LOGICAL                                         :: failure

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_postprocessing', &
         routineP = moduleN//':'//routineN
    TYPE(cp_logger_type), POINTER :: logger
    INTEGER :: unit_nr

    !   ------------------------------------------------------------------------

    failure = .FALSE.
    NULLIFY(harris_energy)

    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    logger => cp_error_get_logger(error)

    IF (.NOT. failure) THEN
       CALL harris_env_get(harris_env=harris_env, harris_energy=harris_energy, error=error)

       ! Output
       IF (logger%para_env%source==logger%para_env%mepos) THEN
          unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.)
          WRITE (unit_nr,*) ""
          WRITE (unit_nr,*) "The Harris functional energy correction is performed!"
          WRITE (unit_nr,*) ""
          WRITE (unit_nr,*) "Ehartree n_elec            =", harris_energy%Ehartree_elec
          WRITE (unit_nr,*) "Ehartree n_core            =", harris_energy%Ehartree_core
          !WRITE (unit_nr,*) "Ehartree n_tot             =", harris_energy%Ehartree_tot
          !WRITE (unit_nr,*) "Ehartree_elec_core         =", harris_energy%Ehartree_elec_core

          WRITE (unit_nr,*) "Exc                        =", harris_energy%Exc

          WRITE (unit_nr,*) "The XC potential integral  =", harris_energy%integral_vxc

          !WRITE (unit_nr,*) "Ecore                      =", harris_energy%Ecore
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
!***************************************************************************

END MODULE harris_functional
