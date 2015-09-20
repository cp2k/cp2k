!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief basic functionality for using ot in the scf routines.
!> \par History
!>      01.2003 : Joost VandeVondele : adapted for LSD
!> \author Joost VandeVondele (25.08.2002)
! *****************************************************************************
MODULE qs_ot_scf
  USE cp_array_utils_r,                ONLY: cp_1d_r_p_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_copy, cp_dbcsr_get_diag, cp_dbcsr_get_info, cp_dbcsr_init, &
       cp_dbcsr_init_p, cp_dbcsr_multiply, cp_dbcsr_p_type, cp_dbcsr_release, &
       cp_dbcsr_scale_by_vector, cp_dbcsr_set, cp_dbcsr_set_diag, &
       cp_dbcsr_trace, cp_dbcsr_type, dbcsr_type_no_symmetry
  USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                             cp_dbcsr_m_by_n_from_row_template
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE qs_mo_occupation,                ONLY: set_mo_occupation
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_restrict
  USE qs_ot,                           ONLY: qs_ot_get_orbitals,&
                                             qs_ot_get_orbitals_ref,&
                                             qs_ot_get_p
  USE qs_ot_minimizer,                 ONLY: ot_mini
  USE qs_ot_types,                     ONLY: ot_readwrite_input,&
                                             qs_ot_allocate,&
                                             qs_ot_destroy,&
                                             qs_ot_init,&
                                             qs_ot_settings_init,&
                                             qs_ot_type
  USE scf_control_types,               ONLY: smear_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_scf'
  ! *** Public subroutines ***

  PUBLIC :: ot_scf_init
  PUBLIC :: ot_scf_mini
  PUBLIC :: ot_scf_destroy
  PUBLIC :: ot_scf_read_input

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param scf_section ...
! *****************************************************************************
  SUBROUTINE ot_scf_read_input(qs_ot_env,scf_section)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(section_vals_type), POINTER         :: scf_section

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_scf_read_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin, &
                                                output_unit
    LOGICAL                                  :: explicit, failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: ot_section

    failure=.FALSE.

    CALL timeset(routineN,handle)

    logger => cp_get_default_logger()
    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".log")

    ! decide default settings
    CALL qs_ot_settings_init(qs_ot_env(1)%settings)

    ! use ot input new style
    ot_section=>section_vals_get_subs_vals(scf_section,"OT")
    CALL section_vals_get(ot_section,explicit=explicit)

    CALL ot_readwrite_input(qs_ot_env(1)%settings,ot_section,output_unit)

    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
         "PRINT%PROGRAM_RUN_INFO")

    ! copy the ot settings type so it is identical
    nspin=SIZE(qs_ot_env)
    DO ispin=2,nspin
       qs_ot_env(ispin)%settings=qs_ot_env(1)%settings
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE ot_scf_read_input
! *****************************************************************************
  !
  ! performs the actual minimisation, needs only limited info
  ! updated for restricted calculations
  ! matrix_dedc is the derivative of the energy with respect to the orbitals (except for a factor 2*fi)
  ! a null pointer for matrix_s implies that matrix_s is the unit matrix
  !
  !
! *****************************************************************************
!> \brief ...
!> \param mo_array ...
!> \param matrix_dedc ...
!> \param smear ...
!> \param matrix_s ...
!> \param energy ...
!> \param energy_only ...
!> \param delta ...
!> \param qs_ot_env ...
!> \param input ...
! *****************************************************************************
  SUBROUTINE ot_scf_mini(mo_array,matrix_dedc,smear,matrix_s,energy, &
       energy_only, delta, qs_ot_env, input)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_dedc
    TYPE(smear_type), POINTER                :: smear
    TYPE(cp_dbcsr_type), POINTER             :: matrix_s
    REAL(KIND=dp)                            :: energy
    LOGICAL, INTENT(INOUT)                   :: energy_only
    REAL(KIND=dp)                            :: delta
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(section_vals_type), POINTER         :: input

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_scf_mini', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, k, n, nspin, &
                                                output_unit, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ener_nondiag, trace
    TYPE(cp_1d_r_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: expectation_values, &
                                                occupation_numbers, &
                                                scaling_factor
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_dedc_scaled
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(logger)
    logger => cp_get_default_logger()

    nspin=SIZE(mo_array)

    ALLOCATE(occupation_numbers(nspin),stat=stat)
    CPASSERT(stat==0)
    ALLOCATE(scaling_factor(nspin),stat=stat)
    CPASSERT(stat==0)

    IF (qs_ot_env(1)%settings%do_ener) THEN
      ALLOCATE(expectation_values(nspin),stat=stat)
      CPASSERT(stat==0)
    ENDIF

    DO ispin=1,nspin
       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,occupation_numbers=occupation_numbers(ispin)%array)
       ALLOCATE(scaling_factor(ispin)%array(SIZE(occupation_numbers(ispin)%array)),stat=stat)
       CPASSERT(stat==0)
       scaling_factor(ispin)%array=2.0_dp*occupation_numbers(ispin)%array
       IF (qs_ot_env(1)%settings%do_ener) THEN
         ALLOCATE(expectation_values(ispin)%array(SIZE(occupation_numbers(ispin)%array)),stat=stat)
         CPASSERT(stat==0)
       ENDIF
    ENDDO

    ! optimizing orbital energies somehow implies non-equivalent orbitals
    IF (qs_ot_env(1)%settings%do_ener) THEN
       CPASSERT(qs_ot_env(1)%settings%do_rotation)
    ENDIF
    ! add_nondiag_energy requires do_ener
    IF (qs_ot_env(1)%settings%add_nondiag_energy) THEN
       CPASSERT(qs_ot_env(1)%settings%do_ener)
    ENDIF

    ! get a rotational force
    IF ( .NOT. energy_only) THEN
       IF (qs_ot_env(1)%settings%do_rotation) THEN
          DO ispin=1,SIZE(qs_ot_env)
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff_b=mo_coeff)
             CALL cp_dbcsr_get_info(mo_coeff,nfullrows_total=n,nfullcols_total=k)
             CALL cp_dbcsr_multiply('T','N',1.0_dp,mo_coeff,matrix_dedc(ispin)%matrix, &
                  0.0_dp,qs_ot_env(ispin)%rot_mat_chc)
             CALL cp_dbcsr_copy(qs_ot_env(ispin)%matrix_buf1,qs_ot_env(ispin)%rot_mat_chc)

             CALL cp_dbcsr_scale_by_vector(qs_ot_env(ispin)%matrix_buf1,alpha=scaling_factor(ispin)%array,side='right')
             ! create the derivative of the energy wrt to rot_mat_u
             CALL cp_dbcsr_multiply('N','N',1.0_dp,qs_ot_env(ispin)%rot_mat_u,qs_ot_env(ispin)%matrix_buf1, &
                  0.0_dp,qs_ot_env(ispin)%rot_mat_dedu)
          ENDDO

          ! here we construct the derivative of the free energy with respect to the evals
          ! (note that this requires the diagonal elements of chc)
          ! the mo occupations should in principle remain unaltered
          IF (qs_ot_env(1)%settings%do_ener) THEN
            DO ispin=1,SIZE(mo_array)
               CALL cp_dbcsr_get_diag(qs_ot_env(ispin)%rot_mat_chc,expectation_values(ispin)%array)
               qs_ot_env(ispin)%ener_gx=expectation_values(ispin)%array
               CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                       smear=smear, eval_deriv=qs_ot_env(ispin)%ener_gx)
            ENDDO
          ENDIF

          ! chc only needs to be stored in u independent form if we require add_nondiag_energy,
          ! which will use it in non-selfconsistent form for e.g. the linesearch
          ! transform C^T H C -> U C^T H C U ^ T
          IF (qs_ot_env(1)%settings%add_nondiag_energy) THEN
             DO ispin=1,SIZE(qs_ot_env)
                CALL cp_dbcsr_get_info(qs_ot_env(ispin)%rot_mat_u,nfullcols_total=k)
                CALL cp_dbcsr_multiply('N','N',1.0_dp,qs_ot_env(ispin)%rot_mat_u,qs_ot_env(ispin)%rot_mat_chc, &
                  0.0_dp,qs_ot_env(ispin)%matrix_buf1)
                CALL cp_dbcsr_multiply('N','T',1.0_dp,qs_ot_env(ispin)%matrix_buf1,qs_ot_env(ispin)%rot_mat_u, &
                  0.0_dp,qs_ot_env(ispin)%rot_mat_chc)
             ENDDO
          ENDIF
       ENDIF
    ENDIF

    ! evaluate non-diagonal energy contribution
    ener_nondiag=0.0_dp
    IF (qs_ot_env(1)%settings%add_nondiag_energy) THEN
       DO ispin=1,SIZE(qs_ot_env)
          ! transform \tilde H to the current basis of C (assuming non-selfconsistent H)
          CALL cp_dbcsr_get_info(qs_ot_env(ispin)%rot_mat_u,nfullcols_total=k)
          CALL cp_dbcsr_multiply('T','N',1.0_dp,qs_ot_env(ispin)%rot_mat_u,qs_ot_env(ispin)%rot_mat_chc, &
            0.0_dp,qs_ot_env(ispin)%matrix_buf1)
          CALL cp_dbcsr_multiply('N','N',1.0_dp,qs_ot_env(ispin)%matrix_buf1,qs_ot_env(ispin)%rot_mat_u, &
            0.0_dp,qs_ot_env(ispin)%matrix_buf2)

          ! subtract the current ener_x from the diagonal
          CALL cp_dbcsr_get_diag(qs_ot_env(ispin)%matrix_buf2,expectation_values(ispin)%array)
          expectation_values(ispin)%array=expectation_values(ispin)%array-qs_ot_env(ispin)%ener_x
          CALL cp_dbcsr_set_diag(qs_ot_env(ispin)%matrix_buf2,expectation_values(ispin)%array)

          ! get nondiag energy trace (D^T D)
          CALL cp_dbcsr_trace(qs_ot_env(ispin)%matrix_buf2,qs_ot_env(ispin)%matrix_buf2,trace)
          ener_nondiag=ener_nondiag+0.5_dp*qs_ot_env(1)%settings%nondiag_energy_strength*trace

          ! get gradient (again ignoring dependencies of H)
          IF (.NOT. energy_only) THEN
             ! first for the ener_x (-2*(diag(C^T H C)-ener_x))
             qs_ot_env(ispin)%ener_gx = qs_ot_env(ispin)%ener_gx - &
                    qs_ot_env(1)%settings%nondiag_energy_strength * expectation_values(ispin)%array

             ! next for the rot_mat_u derivative (2 * k * \tilde H U D)
             CALL cp_dbcsr_multiply('N','N',1.0_dp,qs_ot_env(ispin)%rot_mat_chc,qs_ot_env(ispin)%rot_mat_u,&
                   0.0_dp,qs_ot_env(ispin)%matrix_buf1)
             CALL cp_dbcsr_multiply('N','N',2.0_dp * qs_ot_env(1)%settings%nondiag_energy_strength, &
                   qs_ot_env(ispin)%matrix_buf1,qs_ot_env(ispin)%matrix_buf2,&
                   1.0_dp,qs_ot_env(ispin)%rot_mat_dedu)
          ENDIF
       ENDDO
    ENDIF

    ! this is kind of a hack so far (costly memory wise), we locally recreate the scaled matrix_hc, and
    ! use it in the following, eventually, as occupations numbers get more integrated, it should become possible
    ! to remove this.
    ALLOCATE(matrix_dedc_scaled(SIZE(matrix_dedc)),stat=stat)
    CPASSERT(stat==0)
    DO ispin=1,SIZE(matrix_dedc)
       ALLOCATE(matrix_dedc_scaled(ispin)%matrix)
       CALL cp_dbcsr_init(matrix_dedc_scaled(ispin)%matrix)
       CALL cp_dbcsr_copy(matrix_dedc_scaled(ispin)%matrix,matrix_dedc(ispin)%matrix)

       ! as a preconditioner, one might want to scale only with a constant, not with f(i)
       ! for the convergence criterium, maybe take it back out
       IF (qs_ot_env(1)%settings%occupation_preconditioner) THEN
          scaling_factor(ispin)%array=2.0_dp
       ENDIF
       CALL cp_dbcsr_scale_by_vector(matrix_dedc_scaled(ispin)%matrix,alpha=scaling_factor(ispin)%array,side='right')
    ENDDO

    ! notice we use qs_ot_env(1) for driving all output and the minimization in case of LSD
    qs_ot_env(1)%etotal=energy+ener_nondiag

    output_unit = cp_print_key_unit_nr(logger,input,"DFT%SCF%PRINT%PROGRAM_RUN_INFO",&
                                       extension=".scfLog")
    CALL ot_mini(qs_ot_env,matrix_dedc_scaled,output_unit)

    delta           =qs_ot_env(1)%delta
    energy_only     =qs_ot_env(1)%energy_only

    ! generate the orbitals using the new matrix_x
    DO ispin=1,SIZE(qs_ot_env)
       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff_b=mo_coeff)
       CALL cp_dbcsr_get_info(mo_coeff,nfullrows_total=n,nfullcols_total=k)
       SELECT CASE(qs_ot_env(1)%settings%ot_algorithm)
       CASE("TOD")
          IF (ASSOCIATED(matrix_s)) THEN
             CALL cp_dbcsr_multiply('N','N',1.0_dp,matrix_s,qs_ot_env(ispin)%matrix_x, &
                  0.0_dp,qs_ot_env(ispin)%matrix_sx)
          ELSE
             CALL cp_dbcsr_copy(qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin)%matrix_x)
          ENDIF
          CALL qs_ot_get_p(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin))
          CALL qs_ot_get_orbitals(mo_coeff,qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin))
       CASE("REF")
          CALL qs_ot_get_orbitals_ref(mo_coeff,matrix_s,qs_ot_env(ispin)%matrix_x,&
                           qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin)%matrix_gx_old,&
                           qs_ot_env(ispin)%matrix_dx,qs_ot_env(ispin),qs_ot_env(1),&
                           output_unit)
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,"Algorithm not yet implemented")
       END SELECT
    ENDDO

    IF (qs_ot_env(1)%restricted) THEN
       CALL mo_set_restrict(mo_array,convert_dbcsr=.TRUE.)
    ENDIF
    !
    ! obtain the new set of OT eigenvalues and set the occupations accordingly
    !
    IF (qs_ot_env(1)%settings%do_ener) THEN
      DO ispin=1,SIZE(mo_array)
         mo_array(ispin)%mo_set%eigenvalues=qs_ot_env(ispin)%ener_x
         CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                 smear=smear)
      ENDDO
    ENDIF

    ! cleanup
    CALL cp_print_key_finished_output(output_unit,logger,input,&
         "DFT%SCF%PRINT%PROGRAM_RUN_INFO")
    DO ispin=1,SIZE(scaling_factor)
      DEALLOCATE(scaling_factor(ispin)%array,stat=stat)
      CPASSERT(stat==0)
    ENDDO
    DEALLOCATE(scaling_factor,stat=stat)
    CPASSERT(stat==0)
    IF (qs_ot_env(1)%settings%do_ener) THEN
       DO ispin=1,SIZE(expectation_values)
         DEALLOCATE(expectation_values(ispin)%array,stat=stat)
         CPASSERT(stat==0)
       ENDDO
       DEALLOCATE(expectation_values,stat=stat)
       CPASSERT(stat==0)
    ENDIF
    DEALLOCATE(occupation_numbers,stat=stat)
    CPASSERT(stat==0)
    DO ispin=1,SIZE(matrix_dedc_scaled)
       CALL cp_dbcsr_release(matrix_dedc_scaled(ispin)%matrix)
       DEALLOCATE(matrix_dedc_scaled(ispin)%matrix)
    ENDDO
    DEALLOCATE(matrix_dedc_scaled,stat=stat)
    CPASSERT(stat==0)

    CALL timestop(handle)

  END SUBROUTINE ot_scf_mini
  !
  ! initialises qs_ot_env so that mo_coeff is the current point
  ! and that the mimizization can be started.
  !
! *****************************************************************************
!> \brief ...
!> \param mo_array ...
!> \param matrix_s ...
!> \param qs_ot_env ...
!> \param matrix_ks ...
!> \param broyden_adaptive_sigma ...
! *****************************************************************************
  SUBROUTINE ot_scf_init ( mo_array, matrix_s, qs_ot_env, matrix_ks, broyden_adaptive_sigma)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_dbcsr_type), POINTER             :: matrix_s
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_dbcsr_type), POINTER             :: matrix_ks
    REAL(KIND=dp)                            :: broyden_adaptive_sigma

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_scf_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, k, n, nspin
    LOGICAL                                  :: is_equal
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff
    TYPE(cp_fm_type), POINTER                :: mo_coeff_fm

    CALL timeset(routineN,handle)


    DO ispin=1,SIZE(mo_array)
       IF(.NOT.ASSOCIATED(mo_array(ispin)%mo_set%mo_coeff_b)) THEN
          CALL stop_program(routineN,moduleN,__LINE__,"Shouldn't get there")
          ! we do ot then copy fm to dbcsr
          ! allocate that somewhere else ! fm -> dbcsr
          CALL cp_dbcsr_init_p(mo_array(ispin)%mo_set%mo_coeff_b)
          CALL cp_dbcsr_m_by_n_from_row_template(mo_array(ispin)%mo_set%mo_coeff_b,template=matrix_ks,&
                                                 n=mo_array(ispin)%mo_set%nmo,&
                                                 sym=dbcsr_type_no_symmetry)
       ENDIF
    ENDDO

    ! *** set a history for broyden
    DO ispin=1,SIZE(qs_ot_env)
       qs_ot_env(ispin)%broyden_adaptive_sigma=broyden_adaptive_sigma
    ENDDO

    ! **** SCP
    ! **** SCP
    ! adapted for work with the restricted keyword
    nspin=SIZE(qs_ot_env)

    DO ispin=1,nspin

       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff_b=mo_coeff,mo_coeff=mo_coeff_fm)
       CALL copy_fm_to_dbcsr(mo_coeff_fm,mo_coeff)!fm -> dbcsr

       CALL cp_dbcsr_get_info(mo_coeff, nfullrows_total=n, nfullcols_total=k)

       ! allocate
       CALL qs_ot_allocate(qs_ot_env(ispin),matrix_ks,mo_coeff_fm%matrix_struct)

       ! set c0,sc0
       CALL cp_dbcsr_copy(qs_ot_env(ispin)%matrix_c0,mo_coeff)
       IF (ASSOCIATED(matrix_s)) THEN
          CALL cp_dbcsr_multiply('N','N',1.0_dp,matrix_s,qs_ot_env(ispin)%matrix_c0, &
               0.0_dp,qs_ot_env(ispin)%matrix_sc0)
       ELSE
          CALL cp_dbcsr_copy(qs_ot_env ( ispin ) % matrix_sc0, qs_ot_env ( ispin ) % matrix_c0)
       ENDIF

       ! init
       CALL qs_ot_init(qs_ot_env(ispin))

       ! set x
       CALL cp_dbcsr_set(qs_ot_env(ispin)%matrix_x,0.0_dp)
       CALL cp_dbcsr_set(qs_ot_env(ispin)%matrix_sx,0.0_dp)

       IF (qs_ot_env(ispin)%settings%do_rotation) THEN
          CALL cp_dbcsr_set(qs_ot_env(ispin)%rot_mat_x,0.0_dp)
       ENDIF

       IF (qs_ot_env(ispin)%settings%do_ener) THEN
          is_equal=SIZE(qs_ot_env ( ispin ) % ener_x)==SIZE(mo_array(ispin)%mo_set%eigenvalues)
          CPASSERT(is_equal)
          qs_ot_env ( ispin ) % ener_x = mo_array(ispin)%mo_set%eigenvalues
       ENDIF

       SELECT CASE(qs_ot_env(1)%settings%ot_algorithm)
       CASE("TOD")
          ! get c
          CALL qs_ot_get_p(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin))
       CASE("REF")
          CALL cp_dbcsr_copy(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_c0)
          CALL cp_dbcsr_copy(qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin)%matrix_sc0)
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,"Algorithm not yet implemented")
       END SELECT

    END DO
    CALL timestop(handle)
  END SUBROUTINE ot_scf_init

! *****************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! *****************************************************************************
  SUBROUTINE ot_scf_destroy(qs_ot_env)

    TYPE(qs_ot_type)                         :: qs_ot_env

    CALL qs_ot_destroy(qs_ot_env)

  END SUBROUTINE ot_scf_destroy

END MODULE qs_ot_scf

