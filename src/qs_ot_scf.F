!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief basic functionality for using ot in the scf routines.
!> \par History
!>      01.2003 : Joost VandeVondele : adapted for LSD
!> \author Joost VandeVondele (25.08.2002)
! *****************************************************************************
MODULE qs_ot_scf
  USE cp_array_r_utils,                ONLY: cp_1d_r_p_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_trace
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_p_type, &
       cp_fm_release, cp_fm_set_all, cp_fm_set_diag, cp_fm_to_fm, cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_restrict,&
                                             set_mo_occupation
  USE qs_ot,                           ONLY: qs_ot_get_orbitals,&
                                             qs_ot_get_orbitals_ref,&
                                             qs_ot_get_p,&
                                             qs_ot_get_scp_coeffs
  USE qs_ot_minimizer,                 ONLY: ot_mini
  USE qs_ot_types,                     ONLY: ot_readwrite_input,&
                                             qs_ot_allocate,&
                                             qs_ot_destroy,&
                                             qs_ot_init,&
                                             qs_ot_settings_init,&
                                             qs_ot_type
  USE scf_control_types,               ONLY: smear_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type,&
                                             aux_coeff_type,&
                                             get_aux_coeff
  USE sparse_matrix_types,             ONLY: real_matrix_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_scf'
  ! *** Public subroutines ***

  PUBLIC :: ot_scf_init
  PUBLIC :: ot_scf_mini
  PUBLIC :: ot_scf_destroy
  PUBLIC :: ot_scf_read_input

CONTAINS

! *****************************************************************************
  SUBROUTINE ot_scf_read_input(qs_ot_env,para_env,scf_section,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_scf_read_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin, &
                                                output_unit
    LOGICAL                                  :: explicit, failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: ot_section

    failure=.FALSE.

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".log",error=error)

    ! decide default settings
    CALL qs_ot_settings_init(qs_ot_env(1)%settings)

    ! use ot input new style
    ot_section=>section_vals_get_subs_vals(scf_section,"OT",error=error)
    CALL section_vals_get(ot_section,explicit=explicit,error=error)

    CALL ot_readwrite_input(qs_ot_env(1)%settings,ot_section,output_unit,error)

    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)

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
  SUBROUTINE ot_scf_mini(mo_array,matrix_dedc,smear,matrix_s,energy, &
       aux_coeff_set, energy_only, delta, qs_ot_env,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrix_dedc
    TYPE(smear_type)                         :: smear
    TYPE(real_matrix_type), POINTER          :: matrix_s
    REAL(KIND=dp)                            :: energy
    TYPE(aux_coeff_set_type), OPTIONAL, &
      POINTER                                :: aux_coeff_set
    LOGICAL, INTENT(INOUT)                   :: energy_only
    REAL(KIND=dp)                            :: delta
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_scf_mini', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, k, n, nspin, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ener_nondiag, trace
    TYPE(cp_1d_r_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: expectation_values, &
                                                occupation_numbers, &
                                                scaling_factor
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrix_dedc_scaled
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset(routineN,handle)

    failure=.FALSE.

    nspin=SIZE(mo_array)

    ALLOCATE(occupation_numbers(nspin),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(scaling_factor(nspin),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (qs_ot_env(1)%settings%do_ener) THEN
      ALLOCATE(expectation_values(nspin),stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    DO ispin=1,nspin
       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,occupation_numbers=occupation_numbers(ispin)%array)
       ALLOCATE(scaling_factor(ispin)%array(SIZE(occupation_numbers(ispin)%array)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       scaling_factor(ispin)%array=2.0_dp*occupation_numbers(ispin)%array
       IF (qs_ot_env(1)%settings%do_ener) THEN
         ALLOCATE(expectation_values(ispin)%array(SIZE(occupation_numbers(ispin)%array)),stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF
    ENDDO

    ! optimizing orbital energies somehow implies non-equivalent orbitals
    IF (qs_ot_env(1)%settings%do_ener) THEN
       CPPrecondition(qs_ot_env(1)%settings%do_rotation,cp_failure_level,routineP,error,failure)
    ENDIF
    ! add_nondiag_energy requires do_ener
    IF (qs_ot_env(1)%settings%add_nondiag_energy) THEN
       CPPrecondition(qs_ot_env(1)%settings%do_ener,cp_failure_level,routineP,error,failure)   
    ENDIF

    ! get a rotational force
    IF ( .NOT. energy_only) THEN
       IF (qs_ot_env(1)%settings%do_rotation) THEN
          DO ispin=1,SIZE(qs_ot_env)
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)
             CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k,error=error)
             CALL cp_fm_gemm('T','N',k,k,n,1.0_dp,mo_coeff,matrix_dedc(ispin)%matrix, &
                  0.0_dp,qs_ot_env(ispin)%rot_mat_chc,error=error)
             CALL cp_fm_to_fm(qs_ot_env(ispin)%rot_mat_chc,qs_ot_env(ispin)%matrix_buf1,error=error)
             CALL cp_fm_column_scale(qs_ot_env(ispin)%matrix_buf1,scaling_factor(ispin)%array)
             ! create the derivative of the energy wrt to rot_mat_u
             CALL cp_fm_gemm('N','N',k,k,k,1.0_dp,qs_ot_env(ispin)%rot_mat_u,qs_ot_env(ispin)%matrix_buf1, &
                  0.0_dp,qs_ot_env(ispin)%rot_mat_dedu,error=error)
          ENDDO

          ! here we construct the derivative of the free energy with respect to the evals
          ! (note that this requires the diagonal elements of chc)
          ! the mo occupations should in principle remain unaltered
          IF (qs_ot_env(1)%settings%do_ener) THEN
            DO ispin=1,SIZE(mo_array)
               CALL cp_fm_get_diag(qs_ot_env(ispin)%rot_mat_chc,expectation_values(ispin)%array,error)
               qs_ot_env(ispin)%ener_gx=expectation_values(ispin)%array
               CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                       smear=smear, eval_deriv=qs_ot_env(ispin)%ener_gx, error=error)
            ENDDO
          ENDIF

          ! chc only needs to be stored in u independent form if we require add_nondiag_energy, 
          ! which will use it in non-selfconsistent form for e.g. the linesearch
          ! transform C^T H C -> U C^T H C U ^ T
          IF (qs_ot_env(1)%settings%add_nondiag_energy) THEN
             DO ispin=1,SIZE(qs_ot_env)
                CALL cp_fm_get_info(qs_ot_env(ispin)%rot_mat_u,ncol_global=k,error=error)
                CALL cp_fm_gemm('N','N',k,k,k,1.0_dp,qs_ot_env(ispin)%rot_mat_u,qs_ot_env(ispin)%rot_mat_chc, &
                  0.0_dp,qs_ot_env(ispin)%matrix_buf1,error=error)
                CALL cp_fm_gemm('N','T',k,k,k,1.0_dp,qs_ot_env(ispin)%matrix_buf1,qs_ot_env(ispin)%rot_mat_u, &
                  0.0_dp,qs_ot_env(ispin)%rot_mat_chc,error=error)
             ENDDO
          ENDIF
       ENDIF
    ENDIF

    ! evaluate non-diagonal energy contribution
    ener_nondiag=0.0_dp
    IF (qs_ot_env(1)%settings%add_nondiag_energy) THEN
       DO ispin=1,SIZE(qs_ot_env)
          ! transform \tilde H to the current basis of C (assuming non-selfconsistent H)
          CALL cp_fm_get_info(qs_ot_env(ispin)%rot_mat_u,ncol_global=k,error=error)
          CALL cp_fm_gemm('T','N',k,k,k,1.0_dp,qs_ot_env(ispin)%rot_mat_u,qs_ot_env(ispin)%rot_mat_chc, &
            0.0_dp,qs_ot_env(ispin)%matrix_buf1,error=error)
          CALL cp_fm_gemm('N','N',k,k,k,1.0_dp,qs_ot_env(ispin)%matrix_buf1,qs_ot_env(ispin)%rot_mat_u, &
            0.0_dp,qs_ot_env(ispin)%matrix_buf2,error=error)

          ! subtract the current ener_x from the diagonal
          CALL cp_fm_get_diag(qs_ot_env(ispin)%matrix_buf2,expectation_values(ispin)%array,error)
          expectation_values(ispin)%array=expectation_values(ispin)%array-qs_ot_env(ispin)%ener_x
          CALL cp_fm_set_diag(qs_ot_env(ispin)%matrix_buf2,expectation_values(ispin)%array,error)
 
          ! get nondiag energy trace (D^T D)
          CALL cp_fm_trace(qs_ot_env(ispin)%matrix_buf2,qs_ot_env(ispin)%matrix_buf2,trace,error)
          ener_nondiag=ener_nondiag+0.5_dp*qs_ot_env(1)%settings%nondiag_energy_strength*trace

          ! get gradient (again ignoring dependencies of H)
          IF (.NOT. energy_only) THEN
             ! first for the ener_x (-2*(diag(C^T H C)-ener_x))
             qs_ot_env(ispin)%ener_gx = qs_ot_env(ispin)%ener_gx - &
                    qs_ot_env(1)%settings%nondiag_energy_strength * expectation_values(ispin)%array

             ! next for the rot_mat_u derivative (2 * k * \tilde H U D)
             CALL cp_fm_gemm('N','N',k,k,k,1.0_dp,qs_ot_env(ispin)%rot_mat_chc,qs_ot_env(ispin)%rot_mat_u,&
                   0.0_dp,qs_ot_env(ispin)%matrix_buf1,error=error)
             CALL cp_fm_gemm('N','N',k,k,k,2.0_dp * qs_ot_env(1)%settings%nondiag_energy_strength, &
                   qs_ot_env(ispin)%matrix_buf1,qs_ot_env(ispin)%matrix_buf2,&
                   1.0_dp,qs_ot_env(ispin)%rot_mat_dedu,error=error)
          ENDIF
       ENDDO
    ENDIF

    ! this is kind of a hack so far (costly memory wise), we locally recreate the scaled matrix_hc, and
    ! use it in the following, eventually, as occupations numbers get more integrated, it should become possible
    ! to remove this.
    ALLOCATE(matrix_dedc_scaled(SIZE(matrix_dedc)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,SIZE(matrix_dedc)
       CALL cp_fm_create(matrix_dedc_scaled(ispin)%matrix,matrix_dedc(ispin)%matrix%matrix_struct,error=error)
       CALL cp_fm_to_fm(matrix_dedc(ispin)%matrix,matrix_dedc_scaled(ispin)%matrix,error=error)

       ! as a preconditioner, one might want to scale only with a constant, not with f(i)
       ! for the convergence criterium, maybe take it back out
       IF (qs_ot_env(1)%settings%occupation_preconditioner) THEN 
          scaling_factor(ispin)%array=2.0_dp
       ENDIF
       CALL cp_fm_column_scale(matrix_dedc_scaled(ispin)%matrix,scaling_factor(ispin)%array)
    ENDDO

    ! notice we use qs_ot_env(1) for driving all output and the minimization in case of LSD
    qs_ot_env(1)%etotal=energy+ener_nondiag
    CALL ot_mini(qs_ot_env,matrix_dedc_scaled,aux_coeff_set, error=error)
    delta           =qs_ot_env(1)%delta
    energy_only     =qs_ot_env(1)%energy_only

    ! generate the orbitals using the new matrix_x
    DO ispin=1,SIZE(qs_ot_env)
       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)
       CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k,error=error)
       SELECT CASE(qs_ot_env(1)%settings%ot_algorithm)
       CASE("TOD")
          IF (ASSOCIATED(matrix_s)) THEN
             CALL cp_sm_fm_multiply(matrix_s,qs_ot_env(ispin)%matrix_x, &
                  qs_ot_env(ispin)%matrix_sx, k,error=error)
          ELSE
             CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,error=error)
          ENDIF
          CALL qs_ot_get_p(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin),error=error)
          CALL qs_ot_get_orbitals(mo_coeff,qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin),error=error)
       CASE("REF")
          CALL qs_ot_get_orbitals_ref(mo_coeff,matrix_s,qs_ot_env(ispin)%matrix_x,&
               &           qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin)%matrix_gx_old,&
               &           qs_ot_env(ispin)%matrix_dx,qs_ot_env(ispin),qs_ot_env(1),&
               &           error=error)
       CASE DEFAULT
          CALL stop_program("ot_scf_mini","ALGORITHM NYI")
       END SELECT
    ENDDO
    IF (qs_ot_env(1)%restricted) THEN
       CALL mo_set_restrict(mo_array,error=error)
    ENDIF
    ! **** SCP
    ! generate the SCP aux coeffs using the new x
    IF (qs_ot_env(1)%settings%scp) CALL qs_ot_get_scp_coeffs ( qs_ot_env ( 1 ), aux_coeff_set, error )
    ! **** SCP

    !
    ! obtain the new set of OT eigenvalues and set the occupations accordingly
    !
    IF (qs_ot_env(1)%settings%do_ener) THEN
      DO ispin=1,SIZE(mo_array)
         mo_array(ispin)%mo_set%eigenvalues=qs_ot_env(ispin)%ener_x
         CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set,&
                 smear=smear, error=error)
      ENDDO
    ENDIF

    ! cleanup
    DO ispin=1,SIZE(scaling_factor)
      DEALLOCATE(scaling_factor(ispin)%array,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    DEALLOCATE(scaling_factor,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (qs_ot_env(1)%settings%do_ener) THEN
       DO ispin=1,SIZE(expectation_values)
         DEALLOCATE(expectation_values(ispin)%array,stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDDO
       DEALLOCATE(expectation_values,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    DEALLOCATE(occupation_numbers,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,SIZE(matrix_dedc_scaled)
       CALL cp_fm_release(matrix_dedc_scaled(ispin)%matrix,error=error)
    ENDDO
    DEALLOCATE(matrix_dedc_scaled,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    CALL timestop(handle)

  END SUBROUTINE ot_scf_mini
  !
  ! initialises qs_ot_env so that mo_coeff is the current point
  ! and that the mimizization can be started.
  !
! *****************************************************************************
  SUBROUTINE ot_scf_init ( mo_array, matrix_s, aux_coeff_set, qs_ot_env, error )

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(real_matrix_type), POINTER          :: matrix_s
    TYPE(aux_coeff_set_type), OPTIONAL, &
      POINTER                                :: aux_coeff_set
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_scf_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, icoef, icoef_atom, &
                                                ikind, ispin, k, n, n_els, &
                                                ncoef_atom, nkind, nspin
    LOGICAL                                  :: failure, is_equal
    REAL(dp), DIMENSION(:, :), POINTER       :: c, fc
    TYPE(aux_coeff_type), POINTER            :: local_coeffs
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset(routineN,handle)

    ! **** SCP
    NULLIFY ( local_coeffs, c, fc )
    ! **** SCP
    ! adapted for work with the restricted keyword
    nspin=SIZE(qs_ot_env)

    DO ispin=1,nspin
       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)

       CALL cp_fm_get_info(mo_coeff, nrow_global=n, ncol_global=k,error=error)

       ! **** SCP
       ! allocate
       CALL qs_ot_allocate(qs_ot_env(ispin),mo_coeff%matrix_struct,  &
            aux_coeff_set = aux_coeff_set,    &
            error=error)
       ! **** SCP

       ! set c0,sc0
       CALL cp_fm_to_fm(mo_coeff,qs_ot_env(ispin)%matrix_c0,error=error)
       IF (ASSOCIATED(matrix_s)) THEN
          CALL cp_sm_fm_multiply(matrix_s,qs_ot_env(ispin)%matrix_c0, &
               qs_ot_env(ispin)%matrix_sc0,k,error=error)
       ELSE
          CALL cp_fm_to_fm ( qs_ot_env ( ispin ) % matrix_c0,  &
               qs_ot_env ( ispin ) % matrix_sc0, error = error )
       ENDIF

       ! init
       CALL qs_ot_init(qs_ot_env(ispin),error=error)

       ! set x
       CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_x,0.0_dp,error=error)
       CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_sx,0.0_dp,error=error)

       IF (qs_ot_env(ispin)%settings%do_rotation) THEN
          CALL cp_fm_set_all(qs_ot_env(ispin)%rot_mat_x,0.0_dp,error=error)
       ENDIF
       IF (qs_ot_env(ispin)%settings%scp) THEN
          qs_ot_env ( ispin ) % x ( : ) = 0.0_dp
       ENDIF

       IF (qs_ot_env(ispin)%settings%do_ener) THEN
          is_equal=SIZE(qs_ot_env ( ispin ) % ener_x)==SIZE(mo_array(ispin)%mo_set%eigenvalues)   
          CPPostcondition(is_equal,cp_failure_level,routineP,error,failure)
          qs_ot_env ( ispin ) % ener_x = mo_array(ispin)%mo_set%eigenvalues
       ENDIF

       SELECT CASE(qs_ot_env(1)%settings%ot_algorithm)
       CASE("TOD")
          ! get c
          CALL qs_ot_get_p(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin),error=error)
       CASE("REF")
          CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_c0,qs_ot_env(ispin)%matrix_x,&
               &           error=error)
          CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_sc0,qs_ot_env(ispin)%matrix_sx,&
               &           error=error)
       CASE DEFAULT
          CALL stop_program("ot_scf_init","ALGORITHM NYI")
       END SELECT
       ! **** SCP
       ! set c = qs_ot_env % x and -fc = qs_ot_env % gx  for SCP
       IF (qs_ot_env(1)%settings%scp .AND. ispin == 1) THEN
          icoef = 0
          nkind = SIZE ( aux_coeff_set % coeffs_of_kind )
          DO ikind = 1, nkind
             local_coeffs => aux_coeff_set % coeffs_of_kind ( ikind ) % coeffs
             IF ( ASSOCIATED ( local_coeffs ) ) THEN
                CALL get_aux_coeff ( coeffs = local_coeffs, c = c, fc = fc,  & 
                     n_els = n_els, ncoef_atom = ncoef_atom, &
                     error = error  )
                DO i = 1, n_els
                   DO icoef_atom = 1, ncoef_atom
                      icoef = icoef + 1
                      qs_ot_env ( 1 ) % x ( icoef ) = c ( i, icoef_atom )
                      qs_ot_env ( 1 ) % gx ( icoef ) = -fc ( i, icoef_atom )
                   END DO
                END DO
             END IF
          END DO
       END IF
       ! **** SCP

    END DO
    CALL timestop(handle)
  END SUBROUTINE ot_scf_init

! *****************************************************************************
  SUBROUTINE ot_scf_destroy(qs_ot_env,error)

    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL qs_ot_destroy(qs_ot_env,error=error)

  END SUBROUTINE ot_scf_destroy

END MODULE qs_ot_scf

