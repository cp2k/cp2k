!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Different diagonalization schemes that can be used 
!>        for the iterative solution of the eigenvalue problem
!> \par History
!>      started from routines previously located in the qs_scf module
!>      05.2009 
! *****************************************************************************
MODULE qs_scf_diagonalization

  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_symm,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_add_to_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_sm_fm_interactions,           ONLY: copy_fm_to_sm,&
                                             copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE input_constants,                 ONLY: cholesky_inverse,&
                                             cholesky_off,&
                                             cholesky_reduce,&
                                             cholesky_restore,&
                                             core_guess,&
                                             general_roks,&
                                             high_spin_roks,&
                                             restart_guess
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE qs_diis,                         ONLY: qs_diis_b_step
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_occupation
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_scf_lanczos,                  ONLY: lanczos_refinement,&
                                             lanczos_refinement_2v
  USE qs_scf_methods,                  ONLY: combine_ks_matrices,&
                                             eigensolver,&
                                             simple_eigensolver
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_diagonalization'

  PUBLIC :: do_general_diag, do_roks_diag, &
            do_special_diag, do_ot_diag, &
            do_block_krylov_diag
            

CONTAINS

! *****************************************************************************
! *****************************************************************************
!> \brief the inner loop of scf, specific to diagonalization with S matrix
!>       basically, in goes the ks matrix out goes a new p matrix
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************

  SUBROUTINE do_general_diag(scf_env,mos,matrix_ks,&
                                         matrix_s,scf_control,scf_section,&
                                         diis_step,use_jacobi,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step, use_jacobi
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ispin
    LOGICAL                                  :: do_level_shift
    REAL(kind=dp)                            :: diis_error

    DO ispin=1,SIZE(matrix_ks)
      CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,&
                         scf_env%scf_work1(ispin)%matrix,&
                         error=error)
    END DO

    IF (scf_env%iter_count > 1 .AND. .NOT. scf_env%skip_diis) THEN
      CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_env%scf_work1,&
                          scf_env%scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                          scf_control%eps_diis,scf_control%nmixing,&
                          s_matrix=matrix_s,&
                          scf_section=scf_section,error=error)
    ELSE
      diis_step = .FALSE.
    END IF

    do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
                      ((scf_control%density_guess == core_guess).OR.&
                       (scf_env%iter_count > 1)))

    IF ((scf_env%iter_count > 1).AND.(scf_env%iter_delta < scf_control%diagonalization%eps_jacobi)) THEN
      use_jacobi = .TRUE.
    ELSE
      use_jacobi = .FALSE.
    END IF

    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      IF (use_jacobi) THEN
        scf_env%iter_method = "DIIS/Jacobi"
      ELSE
        scf_env%iter_method = "DIIS/Diag."
      END IF
    ELSE
      IF(scf_env%mixing_method ==1) THEN
        scf_env%iter_param = scf_env%p_mix_alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = "P_Mix/Jacobi"
        ELSE
          scf_env%iter_method = "P_Mix/Diag."
        END IF
      ELSEIF(scf_env%mixing_method>1) THEN
        scf_env%iter_param = scf_env%mixing_store%alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Jacobi"
        ELSE
          scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Diag."
        END IF
      END IF
    END IF

    scf_env%iter_delta = 0.0_dp

    DO ispin=1,SIZE(matrix_ks)
      CALL eigensolver(matrix_ks=scf_env%scf_work1(ispin)%matrix,&
                       mo_set=mos(ispin)%mo_set,&
                       ortho=scf_env%ortho,&
                       work=scf_env%scf_work2,&
                       do_level_shift=do_level_shift,&
                       level_shift=scf_control%level_shift,&
                       cholesky_method=scf_env%cholesky_method,&
                       use_jacobi=use_jacobi,&
                       jacobi_threshold=scf_control%diagonalization%jacobi_threshold,&
                       error=error)

      CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                             smear=scf_control%smear,&
                             error=error)

      CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                    scf_env%p_mix_new(ispin)%matrix,&
                                    error=error)

    END DO ! ispin

  END SUBROUTINE do_general_diag

! *****************************************************************************
!> \brief the inner loop of scf, specific to diagonalization without S matrix
!>       basically, in goes the ks matrix out goes a new p matrix
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE do_special_diag(scf_env,mos,matrix_ks,scf_control,&
                                         scf_section,diis_step,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: ispin
    LOGICAL                                  :: do_level_shift, use_jacobi
    REAL(kind=dp)                            :: diis_error

    DO ispin=1,SIZE(matrix_ks)
      CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,scf_env%scf_work1(ispin)%matrix,error=error)
    END DO
    IF (scf_env%iter_count > 1 .AND. .NOT. scf_env%skip_diis) THEN
      CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_env%scf_work1,&
                          scf_env%scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                          scf_control%eps_diis,scf_control%nmixing,&
                          scf_section=scf_section,&
                          error=error)
    ELSE
      diis_step = .FALSE.
    END IF

    IF ((scf_env%iter_count > 1).AND.(scf_env%iter_delta < scf_control%diagonalization%eps_jacobi)) THEN
      use_jacobi = .TRUE.
    ELSE
      use_jacobi = .FALSE.
    END IF

    do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
                      ((scf_control%density_guess == core_guess).OR.(scf_env%iter_count > 1)))
    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      IF (use_jacobi) THEN
        scf_env%iter_method = "DIIS/Jacobi"
      ELSE
        scf_env%iter_method = "DIIS/Diag."
      END IF
    ELSE
      IF(scf_env%mixing_method ==1) THEN
        scf_env%iter_param = scf_env%p_mix_alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = "P_Mix/Jacobi"
        ELSE
          scf_env%iter_method = "P_Mix/Diag."
        END IF
      ELSEIF(scf_env%mixing_method>1) THEN
        scf_env%iter_param = scf_env%mixing_store%alpha
        IF (use_jacobi) THEN
          scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Jacobi"
        ELSE
          scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Diag."
        END IF
      END IF
    END IF
    scf_env%iter_delta = 0.0_dp

    DO ispin=1,SIZE(matrix_ks)
      CALL simple_eigensolver(matrix_ks=scf_env%scf_work1(ispin)%matrix,&
                              mo_set=mos(ispin)%mo_set,&
                              work=scf_env%scf_work2,&
                              do_level_shift=do_level_shift,&
                              level_shift=scf_control%level_shift,&
                              use_jacobi=use_jacobi,&
                              jacobi_threshold=scf_control%diagonalization%jacobi_threshold,&
                              error=error)
      CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                             smear=scf_control%smear,&
                             error=error)
      CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                    scf_env%p_mix_new(ispin)%matrix,&
                                    error=error)

    END DO ! ispin

  END SUBROUTINE do_special_diag

! *****************************************************************************
!> \brief the inner loop of scf, specific to iterative diagonalization using OT
!>        with S matrix; basically, in goes the ks matrix out goes a new p matrix
!> \par History
!>      10.2008 created [JGH]
! *****************************************************************************
  SUBROUTINE do_ot_diag(scf_env,mos,matrix_ks,matrix_s,&
               scf_control,scf_section,diis_step,error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: homo, ispin, nmo
    LOGICAL                                  :: do_level_shift
    REAL(kind=dp)                            :: diis_error, eps_iter
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    DO ispin=1,SIZE(matrix_ks)
      CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,scf_env%scf_work1(ispin)%matrix,error=error)
    END DO

    IF (scf_env%iter_count > 1 .AND. .NOT. scf_env%skip_diis) THEN
      CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_env%scf_work1,&
                          scf_env%scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                          scf_control%eps_diis,scf_control%nmixing,&
                          s_matrix=matrix_s,&
                          scf_section=scf_section,error=error)
    ELSE
      diis_step = .FALSE.
    END IF

    do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
                      ((scf_control%density_guess == core_guess).OR.&
                       (scf_env%iter_count > 1)))

    eps_iter = scf_control%diagonalization%eps_iter
    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      scf_env%iter_method = "DIIS/OTdiag"
      DO ispin=1,SIZE(matrix_ks)
        CALL copy_fm_to_sm(scf_env%scf_work1(ispin)%matrix,matrix_ks(ispin)%matrix)
      END DO
      eps_iter = MAX(eps_iter,scf_control%diagonalization%eps_adapt*diis_error)
    ELSE
      IF(scf_env%mixing_method ==1) THEN
        scf_env%iter_param = scf_env%p_mix_alpha
        scf_env%iter_method = "P_Mix/OTdiag."
      ELSEIF(scf_env%mixing_method>1) THEN
        scf_env%iter_param = scf_env%mixing_store%alpha
        scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/OTdiag."
      END IF

    END IF

    scf_env%iter_delta = 0.0_dp

    DO ispin=1,SIZE(matrix_ks)
       CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,nmo=nmo,homo=homo)
       CALL ot_eigensolver(matrix_h=matrix_ks(ispin)%matrix,&
                           matrix_s=matrix_s(1)%matrix,&
                           matrix_c=mo_coeff,&
                           preconditioner=scf_env%ot_preconditioner(1)%preconditioner,&
                           eps_gradient=eps_iter,&
                           iter_max=scf_control%diagonalization%max_iter,&
                           silent=.TRUE.,&
                           ot_settings=scf_control%diagonalization%ot_settings,&
                           error=error)

       CALL calculate_subspace_eigenvalues(mo_coeff,matrix_ks(ispin)%matrix,&
              para_env=mo_coeff%matrix_struct%para_env, do_rotation = .TRUE., error=error)

       CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,smear=scf_control%smear,error=error)
       CALL calculate_density_matrix(mos(ispin)%mo_set,scf_env%p_mix_new(ispin)%matrix,error=error)
    END DO

  END SUBROUTINE do_ot_diag


! *****************************************************************************
!> \brief Solve a set restricted open Kohn-Sham (ROKS) equations based on the
!>         alpha and beta Kohn-Sham matrices from unrestricted Kohn-Sham.
!> \note  
!>         this is only a high-spin ROKS.
!> \par History
!>      04.2006 created [MK]
!>      Revised (01.05.06,MK)
! *****************************************************************************
  SUBROUTINE do_roks_diag(scf_env,mos,matrix_ks,matrix_s,&
                                      scf_control,scf_section,diis_step,&
                                      orthogonal_basis,error)

    ! Literature: - C. C. J. Roothaan, Rev. Mod. Phys. 32, 179 (1960)
    !             - M. F. Guest and V. R. Saunders, Mol. Phys. 28(3), 819 (1974)
    !             - M. Filatov and S. Shaik, Chem. Phys. Lett. 288, 689 (1998)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(INOUT)                   :: diis_step
    LOGICAL, INTENT(IN)                      :: orthogonal_basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_roks_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homoa, homob, imo, &
                                                nalpha, nao, nbeta, nmo
    REAL(KIND=dp)                            :: diis_error, level_shift_loc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eiga, eigb, occa, occb
    TYPE(cp_fm_type), POINTER                :: ksa, ksb, mo2ao, moa, mob, &
                                                ortho, work

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    ortho => scf_env%ortho
    work  => scf_env%scf_work2

    ksa => scf_env%scf_work1(1)%matrix
    ksb => scf_env%scf_work1(2)%matrix

    CALL copy_sm_to_fm(matrix_ks(1)%matrix,ksa,error=error)
    CALL copy_sm_to_fm(matrix_ks(2)%matrix,ksb,error=error)

    ! Get MO information

    CALL get_mo_set(mo_set=mos(1)%mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    nelectron=nalpha,&
                    homo=homoa,&
                    eigenvalues=eiga,&
                    occupation_numbers=occa,&
                    mo_coeff=moa)

    CALL get_mo_set(mo_set=mos(2)%mo_set,&
                    nelectron=nbeta,&
                    homo=homob,&
                    eigenvalues=eigb,&
                    occupation_numbers=occb,&
                    mo_coeff=mob)

    ! Define the amount of level-shifting

    IF ((scf_control%level_shift /= 0.0_dp).AND.&
        ((scf_control%density_guess == core_guess).OR.&
         (scf_control%density_guess == restart_guess).OR.&
         (scf_env%iter_count > 1))) THEN
      level_shift_loc = scf_control%level_shift
    ELSE
      level_shift_loc = 0.0_dp
    END IF

    IF ((scf_env%iter_count > 1).OR.&
        (scf_control%density_guess == core_guess).OR.&
        (scf_control%density_guess == restart_guess)) THEN

      ! Transform the spin unrestricted alpha and beta Kohn-Sham matrices
      ! from AO basis to MO basis: K(MO) = C(T)*K(AO)*C

      CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksa,moa,0.0_dp,work,error=error)
      CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,moa,work,0.0_dp,ksa,error=error)

      CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksb,moa,0.0_dp,work,error=error)
      CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,moa,work,0.0_dp,ksb,error=error)

      ! Combine the spin unrestricted alpha and beta Kohn-Sham matrices
      ! in the MO basis

      IF (scf_control%roks_scheme == general_roks) THEN
        CALL combine_ks_matrices(ksa,ksb,occa,occb,scf_control%roks_f,&
                                 nalpha,nbeta,error=error)
      ELSE IF (scf_control%roks_scheme == high_spin_roks) THEN
        CALL combine_ks_matrices(ksa,ksb,occa,occb,scf_control%roks_parameter,&
                                 error=error)
      ELSE
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Unknown ROKS scheme requested")
      END IF

      ! Back-transform the restricted open Kohn-Sham matrix from MO basis
      ! to AO basis

      IF (orthogonal_basis) THEN
        ! Q = C
        mo2ao => moa
      ELSE
        ! Q = S*C
        mo2ao => mob
!MK     CALL copy_sm_to_fm(matrix_s(1)%matrix,work)
!MK     CALL cp_fm_symm("L","U",nao,nao,1.0_dp,work,moa,0.0_dp,mo2ao)
        CALL cp_sm_fm_multiply(matrix_s(1)%matrix,moa,mo2ao,nao,error=error)
      END IF

      ! K(AO) = Q*K(MO)*Q(T)

      CALL cp_fm_gemm("N","T",nao,nao,nao,1.0_dp,ksa,mo2ao,0.0_dp,work,error=error)
      CALL cp_fm_gemm("N","N",nao,nao,nao,1.0_dp,mo2ao,work,0.0_dp,ksa,error=error)

    ELSE

      ! No transformation matrix available, yet. The closed shell part,
      ! i.e. the beta Kohn-Sham matrix in AO basis, is taken.
      ! There might be better choices, anyhow.

      CALL cp_fm_to_fm(ksb,ksa,error=error)

    END IF

    ! Update DIIS buffer and possibly perform DIIS extrapolation step

    IF (scf_env%iter_count > 1) THEN
      IF (orthogonal_basis) THEN
        CALL qs_diis_b_step(diis_buffer=scf_env%scf_diis_buffer,&
                            mo_array=mos,&
                            kc=scf_env%scf_work1,&
                            sc=work,&
                            delta=scf_env%iter_delta,&
                            error_max=diis_error,&
                            diis_step=diis_step,&
                            eps_diis=scf_control%eps_diis,&
                            scf_section=scf_section,&
                            roks=.TRUE.,&
                            error=error)
      ELSE
        CALL qs_diis_b_step(diis_buffer=scf_env%scf_diis_buffer,&
                            mo_array=mos,&
                            kc=scf_env%scf_work1,&
                            sc=work,&
                            delta=scf_env%iter_delta,&
                            error_max=diis_error,&
                            diis_step=diis_step,&
                            eps_diis=scf_control%eps_diis,&
                            scf_section=scf_section,&
                            s_matrix=matrix_s,&
                            roks=.TRUE.,&
                            error=error)
      END IF
    END IF

    IF (diis_step) THEN
      scf_env%iter_param = diis_error
      scf_env%iter_method = "DIIS/Diag."
    ELSE
      IF(scf_env%mixing_method ==1) THEN
        scf_env%iter_param = scf_env%p_mix_alpha
        scf_env%iter_method = "P_Mix/Diag."
      ELSEIF(scf_env%mixing_method>1) THEN
        scf_env%iter_param = scf_env%mixing_store%alpha
        scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Diag."
      END IF
    END IF

    scf_env%iter_delta = 0.0_dp

    IF (level_shift_loc /= 0.0_dp) THEN

      ! Transform the current Kohn-Sham matrix from AO to MO basis
      ! for level-shifting using the current MO set

      CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksa,moa,0.0_dp,work,error=error)
      CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,moa,work,0.0_dp,ksa,error=error)

      ! Apply level-shifting using 50:50 split of the shift (could be relaxed)

      DO imo=homob+1,homoa
        CALL cp_fm_add_to_element(ksa,imo,imo,0.5_dp*level_shift_loc,error)
      END DO
      DO imo=homoa+1,nmo
        CALL cp_fm_add_to_element(ksa,imo,imo,level_shift_loc,error)
      END DO

    ELSE IF (.NOT.orthogonal_basis) THEN

      ! Transform the current Kohn-Sham matrix to an orthogonal basis
      SELECT CASE(scf_env%cholesky_method)
      CASE(cholesky_reduce)
        CALL cp_fm_cholesky_reduce(ksa,ortho,error=error)
      CASE(cholesky_restore)
         CALL cp_fm_upper_to_full(ksa,work,error=error)
         CALL cp_fm_cholesky_restore(ksa,nao,ortho,work,&
                "SOLVE",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,ksa,&
               "SOLVE",pos="LEFT",transa="T",error=error)
      CASE(cholesky_inverse)
         CALL cp_fm_upper_to_full(ksa,work,error=error)
         CALL cp_fm_cholesky_restore(ksa,nao,ortho,work,&
                "MULTIPLY",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,ksa,&
               "MULTIPLY",pos="LEFT",transa="T",error=error)
      CASE(cholesky_off)
        CALL cp_fm_symm("L","U",nao,nao,1.0_dp,ksa,ortho,0.0_dp,work,error=error)
        CALL cp_fm_gemm("N","N",nao,nao,nao,1.0_dp,ortho,work,0.0_dp,ksa,error=error)
      END SELECT

    END IF

    ! Diagonalization of the ROKS operator matrix

    CALL cp_fm_syevd(ksa,work,eiga,error=error)

    ! Back-transformation of the orthonormal eigenvectors if needed

    IF (level_shift_loc /= 0.0_dp) THEN
      ! Use old MO set for back-transformation if level-shifting was applied
      CALL cp_fm_to_fm(moa,ortho,error=error)
      CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,moa,error=error)
    ELSE
      IF (orthogonal_basis) THEN
        CALL cp_fm_to_fm(work,moa,error=error)
      ELSE
        SELECT CASE(scf_env%cholesky_method)
        CASE(cholesky_reduce,cholesky_restore)
          CALL cp_fm_cholesky_restore(work,nmo,ortho,moa,"SOLVE",error=error)
        CASE(cholesky_inverse)
          CALL cp_fm_cholesky_restore(work,nmo,ortho,moa,"MULTIPLY",error=error)
        CASE(cholesky_off) 
          CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,moa,error=error)
        END SELECT
      END IF
    END IF

    ! Correct MO eigenvalues, if level-shifting was applied

    IF (level_shift_loc /= 0.0_dp) THEN
      DO imo=homob+1,homoa
        eiga(imo) = eiga(imo) - 0.5_dp*level_shift_loc
      END DO
      DO imo=homoa+1,nmo
        eiga(imo) = eiga(imo) - level_shift_loc
      END DO
    END IF

    ! Update also the beta MO set

    eigb(:) = eiga(:)
    CALL cp_fm_to_fm(moa,mob,error=error)

    ! Calculate the new alpha and beta density matrix

    CALL calculate_density_matrix(mos(1)%mo_set,scf_env%p_mix_new(1)%matrix,&
                                  error=error)
    CALL calculate_density_matrix(mos(2)%mo_set,scf_env%p_mix_new(2)%matrix,&
                                  error=error)

    CALL timestop(handle)

  END SUBROUTINE do_roks_diag

! *****************************************************************************
!> \brief iterative diagonalization using the block Krylov-space approach
!> \param 
!> \par History
!>      05.2009 created [MI]
! *****************************************************************************

  SUBROUTINE do_block_krylov_diag(scf_env,mos,matrix_ks, matrix_s,&
                     scf_control, scf_section, check_moconv_only, error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(IN), OPTIONAL            :: check_moconv_only
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'do_block_krylov_diag', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, homo, ispin, iter, &
                                                nao, nmo, output_unit
    LOGICAL                                  :: converged, failure, ionode, &
                                                my_check_moconv_only
    REAL(dp)                                 :: eps_iter, t1, t2
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(cp_fm_type), POINTER                :: c0, c1, chc, evec, ks, &
                                                mo_coeff, ortho, work
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
        output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
        output_unit = -1
    END IF

    CALL timeset(routineN,handle)

    output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%LANCZOS",&
            extension=".scfLog",error=error)


    my_check_moconv_only = .FALSE.
    IF(PRESENT(check_moconv_only)) my_check_moconv_only=check_moconv_only

    NULLIFY (mo_coeff,ortho,work, ks)
    NULLIFY (mo_eigenvalues)
    NULLIFY (c0, c1)

    ortho => scf_env%ortho
    work  => scf_env%scf_work2

    DO ispin=1,SIZE(matrix_ks)
      CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,&
                         scf_env%scf_work1(ispin)%matrix,&
                         error=error)
    END DO

    IF(scf_env%mixing_method ==1) THEN
        scf_env%iter_param = scf_env%p_mix_alpha
        scf_env%iter_method = "P_Mix/Lanczos"
    ELSE
        scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/Lanc."
    END IF

    DO ispin=1,SIZE(matrix_ks)
 
      ks => scf_env%scf_work1(ispin)%matrix
      CALL cp_fm_upper_to_full(ks,work,error=error)

      CALL get_mo_set(mo_set=mos(ispin)%mo_set,&
                      nao=nao,&
                      nmo=nmo,&
                      homo=homo,&
                      eigenvalues=mo_eigenvalues,&
                      mo_coeff=mo_coeff)

      NULLIFY(c0,c1)
      c0 => scf_env%krylov_space%mo_conv(ispin)%matrix
      c1 => scf_env%krylov_space%mo_refine(ispin)%matrix
      SELECT CASE(scf_env%cholesky_method)
      CASE(cholesky_reduce)
         CALL cp_fm_cholesky_reduce(ks,ortho,error=error)
         CALL cp_fm_upper_to_full(ks,work,error=error)
         CALL cp_fm_cholesky_restore(mo_coeff,nmo,ortho,c0,"MULTIPLY",error=error)
      CASE(cholesky_restore)
         CALL cp_fm_cholesky_restore(ks,nao,ortho,work,&
                "SOLVE",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,ks,&
               "SOLVE",pos="LEFT",transa="T",error=error)
         CALL cp_fm_cholesky_restore(mo_coeff,nmo,ortho,c0,"MULTIPLY",error=error)
      CASE(cholesky_inverse)
         CALL cp_fm_cholesky_restore(ks,nao,ortho,work,&
                "MULTIPLY",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(work,nao,ortho,ks,&
               "MULTIPLY",pos="LEFT",transa="T",error=error)
         CALL cp_fm_cholesky_restore(mo_coeff,nmo,ortho,c0,"SOLVE",error=error)
      END SELECT

      scf_env%krylov_space%nmo_nc  = nmo
      scf_env%krylov_space%nmo_conv  = 0

      t1 = m_walltime()
      IF (output_unit > 0) THEN
          WRITE(output_unit,"(/T15,A)") '<<<<<<<<<   LANCZOS REFINEMENT    <<<<<<<<<<'
          WRITE(output_unit,"(T8,A,T15,A,T23,A,T36,A,T49,A,T60,A,/,T8,A)")   &
               " Spin ", " Cycle ", &
               " conv. MOS ", " B2MAX ", " B2MIN ", " Time",  REPEAT("-",60)
      END IF
      eps_iter = MAX(scf_env%krylov_space%eps_conv, scf_env%krylov_space%eps_adapt*scf_env%iter_delta)
      iter = 0
      converged = .FALSE.
      !Check convergence of MOS
      IF(my_check_moconv_only) THEN

         CALL lanczos_refinement(scf_env%krylov_space, ks, c0, c1, mo_eigenvalues,& 
              nao,  eps_iter, check_moconv_only=my_check_moconv_only, error=error)
         t2 = m_walltime()
         IF(output_unit > 0)&
               WRITE(output_unit,'(T8,I3,T16,I5,T24,I6,T33,E12.4,2x,E12.4,T60,F8.3)')  &
               ispin,  iter, scf_env%krylov_space%nmo_conv, &
               scf_env%krylov_space%max_res_norm, scf_env%krylov_space%min_res_norm, t2-t1

         CYCLE
      ELSE
         !Block Lanczos refinement
         DO iter = 1,scf_env%krylov_space%max_iter
           CALL lanczos_refinement_2v(scf_env%krylov_space, ks, c0, c1, mo_eigenvalues,& 
                 nao, homo, eps_iter,  error=error)
            t2 = m_walltime()
            IF (output_unit > 0) THEN
               WRITE(output_unit,'(T8,I3,T16,I5,T24,I6,T33,E12.4,2x,E12.4,T60,F8.3)')&
                      ispin, iter, scf_env%krylov_space%nmo_conv, &
                      scf_env%krylov_space%max_res_norm, scf_env%krylov_space%min_res_norm, t2-t1
            END IF
            t1=m_walltime()
            IF(scf_env%krylov_space%max_res_norm < eps_iter) THEN
              converged = .TRUE.
              IF (output_unit > 0) WRITE(output_unit,*)&
                   " Reached convergence in ",iter," iterations "
              EXIT
            END IF
         END DO
    
         IF(.NOT. converged .AND. output_unit > 0) THEN
           WRITE(output_unit,"(T4, A)") " WARNING Lanczos refinement could "//&
                "not converge all the mos:"
           WRITE(output_unit,"(T40,A,T70,I10)")    " number of not converged mos ",&
                 scf_env%krylov_space%nmo_nc 
           WRITE(output_unit,"(T40,A,T70,E10.2)")  " max norm of the residual " , &
                 scf_env%krylov_space%max_res_norm
    
         END IF

         ! For the moment skip the re-orthogonalization
         IF(.FALSE.) THEN
           !Re-orthogonalization
           NULLIFY( chc, evec)
           chc => scf_env%krylov_space%chc_mat
           evec => scf_env%krylov_space%c_vec
           CALL cp_fm_gemm('N','N',nao,nmo,nao,rone,ks,c0,rzero,work,error=error)
           CALL cp_fm_gemm('T','N',nmo,nmo,nao,rone,c0,work,rzero,chc,error=error)
           !Diagonalize  (C^t)HC 
           CALL cp_fm_syevd(chc,evec,mo_eigenvalues,error=error)
           !Rotate the C vectors
           CALL cp_fm_gemm('N','N',nao,nmo,nmo,rone,c0,evec,rzero,c1,error=error)
           c0 => scf_env%krylov_space%mo_refine(ispin)%matrix
         END IF
    
         IF(scf_env%cholesky_method==cholesky_inverse) THEN
           CALL cp_fm_cholesky_restore(c0,nmo,ortho,mo_coeff,"MULTIPLY",error=error)
         ELSE
           CALL cp_fm_cholesky_restore(c0,nmo,ortho,mo_coeff,"SOLVE",error=error)
         END IF
    
         CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,&
                                smear=scf_control%smear,&
                                error=error)

         CALL calculate_density_matrix(mos(ispin)%mo_set,&
                                       scf_env%p_mix_new(ispin)%matrix,&
                                       error=error)
      END  IF
    END DO ! ispin

    IF (output_unit > 0) THEN
          WRITE(output_unit,"(T15,A/)") '<<<<<<<<< END LANCZOS REFINEMENT  <<<<<<<<<<'
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%LANCZOS", error=error)

    CALL timestop(handle)

  END SUBROUTINE do_block_krylov_diag

END MODULE qs_scf_diagonalization
