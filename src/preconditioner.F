!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief computes preconditioners, and implements methods to apply them
!>      currently used in qs_ot
!> \par History
!>      - [UB] 2009-05-13 Adding stable approximate inverse (full and sparse)
!> \author Joost VandeVondele (09.2002)
! **************************************************************************************************
MODULE preconditioner
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_release_p,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_uplo_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_constants,                 ONLY: &
        cholesky_reduce, ot_low_rank_base_lattice_fft, ot_low_rank_base_overlap, &
        ot_precond_fermi_low_rank, ot_precond_full_all, ot_precond_full_all_covariant, &
        ot_precond_full_kinetic, ot_precond_full_single, ot_precond_full_single_inverse, &
        ot_precond_none, ot_precond_s_inverse, ot_precond_solver_default, &
        ot_precond_solver_inv_chol, ot_precond_solver_update
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE lattice_low_rank_preconditioner, ONLY: build_lattice_inverse
   USE low_rank_preconditioner_model,   ONLY: low_rank_inverse_weight,&
                                              low_rank_select_rank
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE preconditioner_apply,            ONLY: apply_preconditioner_cfm_complex,&
                                              apply_preconditioner_dbcsr,&
                                              apply_preconditioner_dbcsr_complex,&
                                              apply_preconditioner_fm
   USE preconditioner_makes,            ONLY: make_complex_full_all,&
                                              make_complex_full_kinetic,&
                                              make_complex_full_s_inverse,&
                                              make_complex_full_single,&
                                              make_complex_full_single_inverse,&
                                              make_preconditioner_matrix
   USE preconditioner_solvers,          ONLY: solve_preconditioner,&
                                              transfer_dbcsr_to_fm,&
                                              transfer_fm_to_dbcsr
   USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                              init_preconditioner,&
                                              preconditioner_p_type,&
                                              preconditioner_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              get_mo_set,&
                                              mo_set_type,&
                                              set_mo_set
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner'
   REAL(KIND=dp), PARAMETER, PRIVATE    :: fermi_low_rank_gap_factor = 4.0_dp, &
                                           fermi_low_rank_degeneracy_tolerance = 1.0E-8_dp, &
                                           fermi_low_rank_min_window = 1.0_dp

   PUBLIC :: make_preconditioner, make_preconditioner_complex_full_all, &
             make_preconditioner_complex_full_kinetic, &
             make_preconditioner_complex_full_s_inverse, &
             make_preconditioner_complex_full_single, &
             make_preconditioner_complex_full_single_inverse, restart_preconditioner, &
             dbcsr_pair_to_cfm
   PUBLIC :: apply_preconditioner, prepare_preconditioner

! The public interface for apply preconditioner, the routines can be found in preconditioner_apply.F
   INTERFACE apply_preconditioner
      MODULE PROCEDURE apply_preconditioner_cfm_complex
      MODULE PROCEDURE apply_preconditioner_dbcsr
      MODULE PROCEDURE apply_preconditioner_dbcsr_complex
      MODULE PROCEDURE apply_preconditioner_fm
   END INTERFACE

! **************************************************************************************************

CONTAINS

! **************************************************************************************************

! **************************************************************************************************
!> \brief Build the overlap-inverse base of the bounded low-rank spectral preconditioner.
!> \param preconditioner_env ...
!> \param matrix_s ...
! **************************************************************************************************
   SUBROUTINE prepare_low_rank_base(preconditioner_env, matrix_s)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_s

      INTEGER                                            :: i, info, j, n, ncol_local, nrow_local, &
                                                            output_unit
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: used_fft
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlap
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_fm_type)                                   :: overlap_fm, work
      TYPE(cp_logger_type), POINTER                      :: logger

      CPASSERT(ASSOCIATED(preconditioner_env%fm))
      CPASSERT(.NOT. ASSOCIATED(preconditioner_env%base_fm))

      SELECT CASE (preconditioner_env%low_rank_base)
      CASE (ot_low_rank_base_overlap)
         ALLOCATE (preconditioner_env%base_fm)
         CALL cp_fm_create(preconditioner_env%base_fm, preconditioner_env%fm%matrix_struct, &
                           name="FERMI_LOW_RANK overlap inverse")

         IF (PRESENT(matrix_s)) THEN
            CPASSERT(ASSOCIATED(matrix_s))
            CALL copy_dbcsr_to_fm(matrix_s, preconditioner_env%base_fm)
            CALL cp_fm_cholesky_decompose(preconditioner_env%base_fm, info_out=info)
            IF (info /= 0) CPABORT("FERMI_LOW_RANK requires a positive-definite overlap matrix")
            CALL cp_fm_cholesky_invert(preconditioner_env%base_fm)
            CALL cp_fm_create(work, preconditioner_env%base_fm%matrix_struct, name="low-rank base work")
            CALL cp_fm_uplo_to_full(preconditioner_env%base_fm, work)
            CALL cp_fm_release(work)
         ELSE
            CALL cp_fm_set_all(preconditioner_env%base_fm, alpha=0.0_dp, beta=1.0_dp)
         END IF

      CASE (ot_low_rank_base_lattice_fft)
         IF (.NOT. PRESENT(matrix_s)) THEN
            CPABORT("FERMI_LOW_RANK_BASE LATTICE_FFT requires a non-unit overlap metric")
         END IF
         CPASSERT(ASSOCIATED(matrix_s))
         CALL cp_fm_get_info(preconditioner_env%fm, nrow_global=n)
         CALL cp_fm_create(overlap_fm, preconditioner_env%fm%matrix_struct, &
                           name="low-rank lattice overlap")
         CALL cp_fm_set_all(overlap_fm, alpha=0.0_dp, beta=0.0_dp)
         CALL copy_dbcsr_to_fm(matrix_s, overlap_fm)
         CALL cp_fm_create(work, overlap_fm%matrix_struct, name="lattice overlap work")
         CALL cp_fm_uplo_to_full(overlap_fm, work)
         CALL cp_fm_release(work)
         CALL cp_fm_get_info(overlap_fm, nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices, local_data=local_data)
         ALLOCATE (overlap(n, n), source=0.0_dp)
         DO j = 1, ncol_local
            DO i = 1, nrow_local
               IF (row_indices(i) <= col_indices(j)) THEN
                  overlap(row_indices(i), col_indices(j)) = local_data(i, j)
               END IF
            END DO
         END DO
         CALL preconditioner_env%para_env%sum(overlap)
         DO j = 1, n
            DO i = j + 1, n
               overlap(i, j) = overlap(j, i)
            END DO
         END DO
         CALL build_lattice_inverse(overlap, preconditioner_env%lattice_dims, &
                                    preconditioner_env%lattice_inverse_k, &
                                    preconditioner_env%lattice_projection_error, used_fft, info)
         CALL cp_fm_release(overlap_fm)
         DEALLOCATE (overlap)
         IF (info /= 0) THEN
            CPABORT("FERMI_LOW_RANK_BASE LATTICE_FFT found a non-positive Fourier overlap block")
         END IF
         IF (.NOT. used_fft) THEN
            CPABORT("FERMI_LOW_RANK_BASE LATTICE_FFT could not initialize its lattice FFT")
         END IF
         preconditioner_env%lattice_block_size = &
            n/PRODUCT(preconditioner_env%lattice_dims)
         IF (preconditioner_env%para_env%is_source()) THEN
            logger => cp_get_default_logger()
            output_unit = cp_logger_get_default_io_unit(logger)
            IF (output_unit > 0) THEN
               WRITE (output_unit, '(T2,A,3(I0,1X),A,I0)') &
                  "OT| FERMI_LOW_RANK lattice FFT cells: ", preconditioner_env%lattice_dims, &
                  "AO block: ", preconditioner_env%lattice_block_size
               WRITE (output_unit, '(T2,A,ES12.4)') &
                  "OT| Block-circulant overlap projection error: ", &
                  preconditioner_env%lattice_projection_error
            END IF
         END IF

      CASE DEFAULT
         CPABORT("Unknown FERMI_LOW_RANK base operator")
      END SELECT

   END SUBROUTINE prepare_low_rank_base

! **************************************************************************************************

! creates a preconditioner for the system (H-energy_homo S)
! this preconditioner is (must be) symmetric positive definite.
! currently uses a atom-block-diagonal form
! each block will be  ....
! might overwrite matrix_h, matrix_t

! **************************************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param precon_type ...
!> \param solver_type ...
!> \param matrix_h ...
!> \param matrix_s ...
!> \param matrix_t ...
!> \param mo_set ...
!> \param energy_gap ...
!> \param convert_precond_to_dbcsr ...
!> \param chol_type ...
!> \par History
!>      09.2014 removed some unused or unfinished methods
!>              removed sparse preconditioners and the
!>              sparse approximate inverse at rev 14341 [Florian Schiffmann]
! **************************************************************************************************
   SUBROUTINE make_preconditioner(preconditioner_env, precon_type, solver_type, matrix_h, matrix_s, &
                                  matrix_t, mo_set, energy_gap, convert_precond_to_dbcsr, chol_type)

      TYPE(preconditioner_type)                          :: preconditioner_env
      INTEGER, INTENT(IN)                                :: precon_type, solver_type
      TYPE(dbcsr_type), POINTER                          :: matrix_h
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_s, matrix_t
      TYPE(mo_set_type), INTENT(IN)                      :: mo_set
      REAL(KIND=dp)                                      :: energy_gap
      LOGICAL, INTENT(IN), OPTIONAL                      :: convert_precond_to_dbcsr
      INTEGER, INTENT(IN), OPTIONAL                      :: chol_type

      CHARACTER(len=*), PARAMETER :: routineN = 'make_preconditioner'

      INTEGER                                            :: handle, k, max_rank, my_solver_type, &
                                                            nao, nhomo, output_unit
      LOGICAL                                            :: my_convert_precond_to_dbcsr, &
                                                            needs_full_spectrum, needs_homo, &
                                                            use_mo_coeff_b
      REAL(KIND=dp)                                      :: energy_homo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues_ot
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: mo_occ
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_b

      CALL timeset(routineN, handle)

      CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b, homo=nhomo)
      use_mo_coeff_b = mo_set%use_mo_coeff_b
      CALL cp_fm_get_info(mo_coeff, ncol_global=k, nrow_global=nao)

      ! Starting some matrix mess, check where to store the result in preconditioner_env, fm or dbcsr_matrix
      my_convert_precond_to_dbcsr = .FALSE.
      IF (PRESENT(convert_precond_to_dbcsr)) my_convert_precond_to_dbcsr = convert_precond_to_dbcsr
      IF (precon_type == ot_precond_fermi_low_rank .AND. &
          preconditioner_env%low_rank_base == ot_low_rank_base_lattice_fft) THEN
         ! The FFT base is an operator rather than a stored DBCSR matrix. Keep the
         ! complementary vectors in their distributed full-matrix representation.
         my_convert_precond_to_dbcsr = .FALSE.
      END IF

      ! Thanks to the mess with the matrices we need to make sure in this case that the
      ! Previous inverse is properly stored as a sparse matrix, fm gets deallocated here
      ! if it wasn't anyway
      IF (preconditioner_env%solver == ot_precond_solver_update) THEN
         CALL transfer_fm_to_dbcsr(preconditioner_env%fm, preconditioner_env%dbcsr_matrix, matrix_h)
      END IF

      needs_full_spectrum = .FALSE.
      needs_homo = .FALSE.

      SELECT CASE (precon_type)
      CASE (ot_precond_full_all, ot_precond_fermi_low_rank, ot_precond_full_all_covariant)
         needs_full_spectrum = .TRUE.
         ! both of them need the coefficients as fm's, more matrix mess
         IF (use_mo_coeff_b) THEN
            CALL copy_dbcsr_to_fm(mo_coeff_b, mo_coeff)
         END IF
      CASE (ot_precond_full_single)
         needs_homo = .TRUE.
         ! XXXX to be removed if homo estimate only is implemented
         needs_full_spectrum = .TRUE.
      CASE (ot_precond_full_kinetic, ot_precond_s_inverse, ot_precond_full_single_inverse)
         ! these should be happy without an estimate for the homo energy
         ! preconditioning can  not depend on an absolute eigenvalue, only on eigenvalue differences
      CASE DEFAULT
         CPABORT("The preconditioner is unknown ...")
      END SELECT

      ALLOCATE (eigenvalues_ot(k))
      energy_homo = 0.0_dp
      IF (needs_full_spectrum) THEN
         ! XXXXXXXXXXXXXXXX do not touch the initial MOs, could be harmful for either
         !                  the case of non-equivalent MOs but also for the derivate
         ! we could already have all eigenvalues e.g. full_all and we could skip this
         ! to be optimised later.
         ! one flaw is that not all SCF methods (i.e. that go over mo_derivs directly)
         ! have a 'valid' matrix_h... (we even don't know what evals are in that case)
         IF (use_mo_coeff_b) THEN
            CALL calculate_subspace_eigenvalues(mo_coeff_b, matrix_h, &
                                                eigenvalues_ot, do_rotation=.FALSE., &
                                                para_env=mo_coeff%matrix_struct%para_env, &
                                                blacs_env=mo_coeff%matrix_struct%context)
         ELSE
            CALL calculate_subspace_eigenvalues(mo_coeff, matrix_h, &
                                                eigenvalues_ot, do_rotation=.FALSE.)
         END IF
         IF (k > 0) THEN
            CPASSERT(nhomo > 0 .AND. nhomo <= k)
            energy_homo = eigenvalues_ot(nhomo)
         END IF
      ELSE
         IF (needs_homo) THEN
            CPABORT("Not yet implemented")
         END IF
      END IF

      ! After all bits and pieces of checking and initialization, here comes the
      ! part where the preconditioner matrix gets created and solved.
      ! This will give the matrices for later use
      my_solver_type = solver_type
      preconditioner_env%in_use = precon_type
      preconditioner_env%cholesky_use = cholesky_reduce
      IF (PRESENT(chol_type)) preconditioner_env%cholesky_use = chol_type
      preconditioner_env%in_use = precon_type
      IF (nhomo == k) THEN
         CALL make_preconditioner_matrix(preconditioner_env, matrix_h, matrix_s, matrix_t, mo_coeff, &
                                         energy_homo, eigenvalues_ot, energy_gap, my_solver_type)
      ELSE
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nhomo, &
                                  context=preconditioner_env%ctxt, &
                                  para_env=preconditioner_env%para_env)
         CALL cp_fm_create(mo_occ, fm_struct)
         CALL cp_fm_to_fm(mo_coeff, mo_occ, nhomo)
         CALL cp_fm_struct_release(fm_struct)
         !
         CALL make_preconditioner_matrix(preconditioner_env, matrix_h, matrix_s, matrix_t, mo_occ, &
                                         energy_homo, eigenvalues_ot(1:nhomo), energy_gap, my_solver_type)
         !
         CALL cp_fm_release(mo_occ)
      END IF

      CALL solve_preconditioner(my_solver_type, preconditioner_env, matrix_s, matrix_h)

      IF (precon_type == ot_precond_fermi_low_rank) THEN
         preconditioner_env%spectral_reference = energy_homo
         preconditioner_env%spectral_window = &
            MAX(fermi_low_rank_min_window, fermi_low_rank_gap_factor*preconditioner_env%energy_gap)
         preconditioner_env%spectral_base_scale = 1.0_dp/preconditioner_env%spectral_window
         max_rank = preconditioner_env%spectral_max_rank
         IF (max_rank == -1) max_rank = SIZE(preconditioner_env%full_evals) - nhomo
         preconditioner_env%spectral_rank = low_rank_select_rank( &
                                            preconditioner_env%full_evals, nhomo, max_rank, &
                                            fermi_low_rank_degeneracy_tolerance)
         IF (preconditioner_env%low_rank_base == ot_low_rank_base_lattice_fft) THEN
            ! With an approximate base, retain only non-negative spectral updates.
            ! This preserves positive definiteness independently of the projection error.
            DO WHILE (preconditioner_env%spectral_rank > 0)
               IF (low_rank_inverse_weight( &
                   preconditioner_env%full_evals(nhomo + preconditioner_env%spectral_rank), &
                   preconditioner_env%spectral_reference, preconditioner_env%energy_gap) >= &
                   preconditioner_env%spectral_base_scale) EXIT
               preconditioner_env%spectral_rank = preconditioner_env%spectral_rank - 1
            END DO
         END IF
         IF (preconditioner_env%para_env%is_source()) THEN
            logger => cp_get_default_logger()
            output_unit = cp_logger_get_default_io_unit(logger)
            IF (output_unit > 0) THEN
               WRITE (output_unit, '(T2,A,F14.8,A,I0)') &
                  "OT| FERMI_LOW_RANK common reference [hartree]: ", &
                  preconditioner_env%spectral_reference, " retained rank: ", &
                  preconditioner_env%spectral_rank
            END IF
         END IF
         CALL prepare_low_rank_base(preconditioner_env, matrix_s)
      END IF

      ! Here comes more matrix mess, make sure to output the correct matrix format,
      ! A bit pointless to convert the cholesky factorized version as it doesn't work in
      ! dbcsr form and will crash later,...
      IF (my_convert_precond_to_dbcsr) THEN
         CALL transfer_fm_to_dbcsr(preconditioner_env%fm, preconditioner_env%dbcsr_matrix, matrix_h)
         IF (precon_type == ot_precond_fermi_low_rank .AND. &
             preconditioner_env%low_rank_base == ot_low_rank_base_overlap) THEN
            CALL transfer_fm_to_dbcsr(preconditioner_env%base_fm, &
                                      preconditioner_env%base_dbcsr_matrix, matrix_h)
         END IF
      ELSE
         CALL transfer_dbcsr_to_fm(preconditioner_env%dbcsr_matrix, preconditioner_env%fm, &
                                   preconditioner_env%para_env, preconditioner_env%ctxt)
      END IF

      DEALLOCATE (eigenvalues_ot)

      CALL timestop(handle)

   END SUBROUTINE make_preconditioner

! **************************************************************************************************
!> \brief Construct FULL_ALL directly from one complex H(k), S(k), and C(k) channel.
!> \param preconditioner_env preconditioner storage
!> \param matrix_c_re real part of the active OT reference orbitals
!> \param matrix_c_im imaginary part of the active OT reference orbitals
!> \param matrix_h_re real part of H(k)
!> \param matrix_h_im imaginary part of H(k)
!> \param matrix_s_re real part of S(k)
!> \param matrix_s_im imaginary part of S(k)
!> \param mo_set reference-orbital energy labels
!> \param energy_gap denominator floor
!> \param solver_type inverse application selected in the OT input
! **************************************************************************************************
   SUBROUTINE make_preconditioner_complex_full_all(preconditioner_env, &
                                                   matrix_c_re, matrix_c_im, &
                                                   matrix_h_re, matrix_h_im, &
                                                   matrix_s_re, matrix_s_im, &
                                                   mo_set, energy_gap, solver_type)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type), POINTER                          :: matrix_c_re, matrix_c_im, matrix_h_re, &
                                                            matrix_h_im, matrix_s_re, matrix_s_im
      TYPE(mo_set_type), INTENT(IN)                      :: mo_set
      REAL(KIND=dp), INTENT(IN)                          :: energy_gap
      INTEGER, INTENT(IN)                                :: solver_type

      CHARACTER(len=*), PARAMETER :: routineN = 'make_preconditioner_complex_full_all'

      INTEGER                                            :: handle, k
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cp_cfm_type)                                  :: matrix_c, matrix_h, matrix_s

      CALL timeset(routineN, handle)

      NULLIFY (eigenvalues)
      CPASSERT(ASSOCIATED(matrix_c_re))
      CPASSERT(ASSOCIATED(matrix_c_im))
      CPASSERT(ASSOCIATED(matrix_h_re))
      CPASSERT(ASSOCIATED(matrix_h_im))
      CPASSERT(ASSOCIATED(matrix_s_re))
      CPASSERT(ASSOCIATED(matrix_s_im))
      IF (solver_type /= ot_precond_solver_default) THEN
         CPABORT('Complex FULL_ALL supports only PRECOND_SOLVER DEFAULT')
      END IF
      CALL dbcsr_get_info(matrix_c_re, nfullcols_total=k)
      CPASSERT(k > 0)
      CALL get_mo_set(mo_set, eigenvalues=eigenvalues)
      CPASSERT(ASSOCIATED(eigenvalues))
      CPASSERT(SIZE(eigenvalues) >= k)

      CALL dbcsr_pair_to_cfm(matrix_c_re, matrix_c_im, preconditioner_env, &
                             'complex FULL_ALL C', matrix_c)
      CALL dbcsr_pair_to_cfm(matrix_h_re, matrix_h_im, preconditioner_env, &
                             'complex FULL_ALL H', matrix_h)
      CALL dbcsr_pair_to_cfm(matrix_s_re, matrix_s_im, preconditioner_env, &
                             'complex FULL_ALL S', matrix_s)
      CALL make_complex_full_all(preconditioner_env, matrix_c, matrix_h, matrix_s, &
                                 eigenvalues(1:k), energy_gap)

      CALL cp_cfm_release(matrix_c)
      CALL cp_cfm_release(matrix_s)
      CALL cp_cfm_release(matrix_h)

      CALL timestop(handle)

   END SUBROUTINE make_preconditioner_complex_full_all

! **************************************************************************************************
!> \brief Construct a complex FULL_SINGLE preconditioner from H(k) and S(k).
!> \param preconditioner_env preconditioner storage
!> \param matrix_h_re real part of H(k)
!> \param matrix_h_im imaginary part of H(k)
!> \param matrix_s_re real part of S(k)
!> \param matrix_s_im imaginary part of S(k)
!> \param mo_set orbital energy labels defining the occupied edge
!> \param energy_gap denominator floor
!> \param solver_type requested inverse solver
! **************************************************************************************************
   SUBROUTINE make_preconditioner_complex_full_single(preconditioner_env, &
                                                      matrix_h_re, matrix_h_im, &
                                                      matrix_s_re, matrix_s_im, &
                                                      mo_set, energy_gap, solver_type)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type), POINTER                          :: matrix_h_re, matrix_h_im, matrix_s_re, &
                                                            matrix_s_im
      TYPE(mo_set_type), INTENT(IN)                      :: mo_set
      REAL(KIND=dp), INTENT(IN)                          :: energy_gap
      INTEGER, INTENT(IN)                                :: solver_type

      INTEGER                                            :: homo
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cp_cfm_type)                                  :: matrix_h, matrix_s

      NULLIFY (eigenvalues)
      CPASSERT(ASSOCIATED(matrix_h_re))
      CPASSERT(ASSOCIATED(matrix_h_im))
      CPASSERT(ASSOCIATED(matrix_s_re))
      CPASSERT(ASSOCIATED(matrix_s_im))
      IF (solver_type /= ot_precond_solver_default) THEN
         CPABORT('Complex FULL_SINGLE supports only PRECOND_SOLVER DEFAULT')
      END IF
      CALL get_mo_set(mo_set, homo=homo, eigenvalues=eigenvalues)
      CPASSERT(ASSOCIATED(eigenvalues))
      CPASSERT(homo > 0 .AND. homo <= SIZE(eigenvalues))

      CALL dbcsr_pair_to_cfm(matrix_h_re, matrix_h_im, preconditioner_env, &
                             'complex FULL_SINGLE H', matrix_h)
      CALL dbcsr_pair_to_cfm(matrix_s_re, matrix_s_im, preconditioner_env, &
                             'complex FULL_SINGLE S', matrix_s)
      CALL make_complex_full_single(preconditioner_env, matrix_h, matrix_s, &
                                    eigenvalues(homo), energy_gap)
      CALL cp_cfm_release(matrix_s)
      CALL cp_cfm_release(matrix_h)

   END SUBROUTINE make_preconditioner_complex_full_single

! **************************************************************************************************
!> \brief Construct a complex FULL_SINGLE_INVERSE preconditioner without discarding Im(H,S,C).
!> \param preconditioner_env preconditioner storage
!> \param matrix_c_re real part of the occupied reference orbitals
!> \param matrix_c_im imaginary part of the occupied reference orbitals
!> \param matrix_h_re real part of H(k)
!> \param matrix_h_im imaginary part of H(k)
!> \param matrix_s_re real part of S(k)
!> \param matrix_s_im imaginary part of S(k)
!> \param energy_gap lower spectral bound
!> \param solver_type requested inverse solver
! **************************************************************************************************
   SUBROUTINE make_preconditioner_complex_full_single_inverse(preconditioner_env, &
                                                              matrix_c_re, matrix_c_im, &
                                                              matrix_h_re, matrix_h_im, &
                                                              matrix_s_re, matrix_s_im, &
                                                              energy_gap, solver_type)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type), POINTER                          :: matrix_c_re, matrix_c_im, matrix_h_re, &
                                                            matrix_h_im, matrix_s_re, matrix_s_im
      REAL(KIND=dp), INTENT(IN)                          :: energy_gap
      INTEGER, INTENT(IN)                                :: solver_type

      INTEGER                                            :: k, n
      TYPE(cp_cfm_type)                                  :: matrix_c, matrix_h, matrix_s
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_mo, fm_struct_square
      TYPE(cp_fm_type)                                   :: matrix_c_im_fm, matrix_c_re_fm, &
                                                            matrix_h_im_fm, matrix_h_re_fm, &
                                                            matrix_s_im_fm, matrix_s_re_fm

      NULLIFY (fm_struct_mo, fm_struct_square)
      CPASSERT(ASSOCIATED(matrix_c_re))
      CPASSERT(ASSOCIATED(matrix_c_im))
      CPASSERT(ASSOCIATED(matrix_h_re))
      CPASSERT(ASSOCIATED(matrix_h_im))
      CPASSERT(ASSOCIATED(matrix_s_re))
      CPASSERT(ASSOCIATED(matrix_s_im))
      IF (solver_type /= ot_precond_solver_default .AND. &
          solver_type /= ot_precond_solver_inv_chol) THEN
         CPABORT("Complex FULL_SINGLE_INVERSE supports only PRECOND_SOLVER DEFAULT/INVERSE_CHOLESKY")
      END IF
      CALL dbcsr_get_info(matrix_c_re, nfullrows_total=n, nfullcols_total=k)
      CPASSERT(n > 0 .AND. k > 0 .AND. k <= n)

      CALL cp_fm_struct_create(fm_struct_mo, nrow_global=n, ncol_global=k, &
                               context=preconditioner_env%ctxt, &
                               para_env=preconditioner_env%para_env)
      CALL cp_fm_struct_create(fm_struct_square, nrow_global=n, ncol_global=n, &
                               context=preconditioner_env%ctxt, &
                               para_env=preconditioner_env%para_env)
      CALL cp_fm_create(matrix_c_re_fm, fm_struct_mo, name='complex preconditioner C real')
      CALL cp_fm_create(matrix_c_im_fm, fm_struct_mo, name='complex preconditioner C imaginary')
      CALL cp_fm_create(matrix_h_re_fm, fm_struct_square, name='complex preconditioner H real')
      CALL cp_fm_create(matrix_h_im_fm, fm_struct_square, name='complex preconditioner H imaginary')
      CALL cp_fm_create(matrix_s_re_fm, fm_struct_square, name='complex preconditioner S real')
      CALL cp_fm_create(matrix_s_im_fm, fm_struct_square, name='complex preconditioner S imaginary')
      CALL cp_fm_struct_release(fm_struct_mo)
      CALL cp_fm_struct_release(fm_struct_square)
      CALL copy_dbcsr_to_fm(matrix_c_re, matrix_c_re_fm)
      CALL copy_dbcsr_to_fm(matrix_c_im, matrix_c_im_fm)
      CALL copy_dbcsr_to_fm(matrix_h_re, matrix_h_re_fm)
      CALL copy_dbcsr_to_fm(matrix_h_im, matrix_h_im_fm)
      CALL copy_dbcsr_to_fm(matrix_s_re, matrix_s_re_fm)
      CALL copy_dbcsr_to_fm(matrix_s_im, matrix_s_im_fm)

      CALL cp_cfm_create(matrix_c, matrix_c_re_fm%matrix_struct, name='complex preconditioner C')
      CALL cp_cfm_create(matrix_h, matrix_h_re_fm%matrix_struct, name='complex preconditioner H')
      CALL cp_cfm_create(matrix_s, matrix_s_re_fm%matrix_struct, name='complex preconditioner S')
      CALL cp_fm_to_cfm(matrix_c_re_fm, matrix_c_im_fm, matrix_c)
      CALL cp_fm_to_cfm(matrix_h_re_fm, matrix_h_im_fm, matrix_h)
      CALL cp_fm_to_cfm(matrix_s_re_fm, matrix_s_im_fm, matrix_s)
      CALL make_complex_full_single_inverse(preconditioner_env, matrix_c, matrix_h, matrix_s, &
                                            energy_gap)
      preconditioner_env%solver = solver_type

      CALL cp_cfm_release(matrix_s)
      CALL cp_cfm_release(matrix_h)
      CALL cp_cfm_release(matrix_c)
      CALL cp_fm_release(matrix_s_im_fm)
      CALL cp_fm_release(matrix_s_re_fm)
      CALL cp_fm_release(matrix_h_im_fm)
      CALL cp_fm_release(matrix_h_re_fm)
      CALL cp_fm_release(matrix_c_im_fm)
      CALL cp_fm_release(matrix_c_re_fm)

   END SUBROUTINE make_preconditioner_complex_full_single_inverse

! **************************************************************************************************
!> \brief Construct a complex FULL_S_INVERSE preconditioner.
!> \param preconditioner_env preconditioner storage
!> \param matrix_s_re real part of S(k)
!> \param matrix_s_im imaginary part of S(k)
!> \param solver_type requested inverse solver
! **************************************************************************************************
   SUBROUTINE make_preconditioner_complex_full_s_inverse(preconditioner_env, &
                                                         matrix_s_re, matrix_s_im, solver_type)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type), POINTER                          :: matrix_s_re, matrix_s_im
      INTEGER, INTENT(IN)                                :: solver_type

      TYPE(cp_cfm_type)                                  :: matrix_s

      IF (solver_type /= ot_precond_solver_default .AND. &
          solver_type /= ot_precond_solver_inv_chol) THEN
         CPABORT("Complex FULL_S_INVERSE supports only PRECOND_SOLVER DEFAULT/INVERSE_CHOLESKY")
      END IF
      CALL dbcsr_pair_to_cfm(matrix_s_re, matrix_s_im, preconditioner_env, &
                             'complex preconditioner S', matrix_s)
      CALL make_complex_full_s_inverse(preconditioner_env, matrix_s)
      preconditioner_env%solver = solver_type
      CALL cp_cfm_release(matrix_s)

   END SUBROUTINE make_preconditioner_complex_full_s_inverse

! **************************************************************************************************
!> \brief Construct a complex FULL_KINETIC preconditioner.
!> \param preconditioner_env preconditioner storage
!> \param matrix_t_re real part of T(k)
!> \param matrix_t_im imaginary part of T(k)
!> \param matrix_s_re real part of S(k)
!> \param matrix_s_im imaginary part of S(k)
!> \param energy_gap non-negative overlap shift
!> \param solver_type requested inverse solver
! **************************************************************************************************
   SUBROUTINE make_preconditioner_complex_full_kinetic(preconditioner_env, &
                                                       matrix_t_re, matrix_t_im, &
                                                       matrix_s_re, matrix_s_im, &
                                                       energy_gap, solver_type)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type), POINTER                          :: matrix_t_re, matrix_t_im, matrix_s_re, &
                                                            matrix_s_im
      REAL(KIND=dp), INTENT(IN)                          :: energy_gap
      INTEGER, INTENT(IN)                                :: solver_type

      TYPE(cp_cfm_type)                                  :: matrix_s, matrix_t

      IF (solver_type /= ot_precond_solver_default .AND. &
          solver_type /= ot_precond_solver_inv_chol) THEN
         CPABORT("Complex FULL_KINETIC supports only PRECOND_SOLVER DEFAULT/INVERSE_CHOLESKY")
      END IF
      CALL dbcsr_pair_to_cfm(matrix_t_re, matrix_t_im, preconditioner_env, &
                             'complex preconditioner T', matrix_t)
      CALL dbcsr_pair_to_cfm(matrix_s_re, matrix_s_im, preconditioner_env, &
                             'complex preconditioner S', matrix_s)
      CALL make_complex_full_kinetic(preconditioner_env, matrix_t, matrix_s, energy_gap)
      preconditioner_env%solver = solver_type
      CALL cp_cfm_release(matrix_s)
      CALL cp_cfm_release(matrix_t)

   END SUBROUTINE make_preconditioner_complex_full_kinetic

! **************************************************************************************************
!> \brief Copy a real/imaginary DBCSR pair to one distributed complex full matrix.
!> \param matrix_re real matrix component
!> \param matrix_im imaginary matrix component
!> \param preconditioner_env source of the distribution context
!> \param matrix_name matrix label
!> \param matrix complex output matrix
! **************************************************************************************************
   SUBROUTINE dbcsr_pair_to_cfm(matrix_re, matrix_im, preconditioner_env, matrix_name, matrix)

      TYPE(dbcsr_type), POINTER                          :: matrix_re, matrix_im
      TYPE(preconditioner_type)                          :: preconditioner_env
      CHARACTER(LEN=*), INTENT(IN)                       :: matrix_name
      TYPE(cp_cfm_type), INTENT(OUT)                     :: matrix

      INTEGER                                            :: ncol, nrow
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: matrix_im_fm, matrix_re_fm

      NULLIFY (fm_struct)
      CPASSERT(ASSOCIATED(matrix_re))
      CPASSERT(ASSOCIATED(matrix_im))
      CALL dbcsr_get_info(matrix_re, nfullrows_total=nrow, nfullcols_total=ncol)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nrow, ncol_global=ncol, &
                               context=preconditioner_env%ctxt, &
                               para_env=preconditioner_env%para_env)
      CALL cp_fm_create(matrix_re_fm, fm_struct, name=TRIM(matrix_name)//' real')
      CALL cp_fm_create(matrix_im_fm, fm_struct, name=TRIM(matrix_name)//' imaginary')
      CALL cp_fm_struct_release(fm_struct)
      CALL copy_dbcsr_to_fm(matrix_re, matrix_re_fm)
      CALL copy_dbcsr_to_fm(matrix_im, matrix_im_fm)
      CALL cp_cfm_create(matrix, matrix_re_fm%matrix_struct, name=matrix_name)
      CALL cp_fm_to_cfm(matrix_re_fm, matrix_im_fm, matrix)
      CALL cp_fm_release(matrix_im_fm)
      CALL cp_fm_release(matrix_re_fm)

   END SUBROUTINE dbcsr_pair_to_cfm

! **************************************************************************************************
!> \brief Allows for a restart of the preconditioner
!>        depending on the method it purges all arrays or keeps them
!> \param qs_env ...
!> \param preconditioner ...
!> \param prec_type ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE restart_preconditioner(qs_env, preconditioner, prec_type, nspins)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(preconditioner_p_type), DIMENSION(:), POINTER :: preconditioner
      INTEGER, INTENT(IN)                                :: prec_type, nspins

      INTEGER                                            :: ispin
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (para_env, blacs_env)
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      IF (ASSOCIATED(preconditioner)) THEN
         SELECT CASE (prec_type)
         CASE (ot_precond_full_all, ot_precond_fermi_low_rank, &
               ot_precond_full_all_covariant, &
               ot_precond_full_single) ! these depend on the ks matrix
            DO ispin = 1, SIZE(preconditioner)
               CALL destroy_preconditioner(preconditioner(ispin)%preconditioner)
               DEALLOCATE (preconditioner(ispin)%preconditioner)
            END DO
            DEALLOCATE (preconditioner)
         CASE (ot_precond_none, ot_precond_full_kinetic, ot_precond_s_inverse, &
               ot_precond_full_single_inverse) ! these are 'independent'
            ! do nothing
         CASE DEFAULT
            CPABORT("Unknown preconditioner type")
         END SELECT
      END IF

      ! add an OT preconditioner if none is present
      IF (.NOT. ASSOCIATED(preconditioner)) THEN
         SELECT CASE (prec_type)
         CASE (ot_precond_full_all, ot_precond_fermi_low_rank, ot_precond_full_all_covariant, &
               ot_precond_full_single_inverse)
            ALLOCATE (preconditioner(nspins))
         CASE DEFAULT
            ALLOCATE (preconditioner(1))
         END SELECT
         DO ispin = 1, SIZE(preconditioner)
            ALLOCATE (preconditioner(ispin)%preconditioner)
            CALL init_preconditioner(preconditioner(ispin)%preconditioner, &
                                     para_env=para_env, &
                                     blacs_env=blacs_env)
         END DO
      END IF

   END SUBROUTINE restart_preconditioner

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mos ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param ot_preconditioner ...
!> \param prec_type ...
!> \param solver_type ...
!> \param energy_gap ...
!> \param nspins ...
!> \param has_unit_metric ...
!> \param convert_to_dbcsr ...
!> \param chol_type ...
!> \param full_mo_set ...
!> \param chebyshev_degree ...
!> \param low_rank_base ...
!> \param fermi_low_rank_max_rank ...
! **************************************************************************************************
   SUBROUTINE prepare_preconditioner(qs_env, mos, matrix_ks, matrix_s, &
                                     ot_preconditioner, prec_type, solver_type, &
                                     energy_gap, nspins, has_unit_metric, &
                                     convert_to_dbcsr, chol_type, full_mo_set, chebyshev_degree, &
                                     low_rank_base, fermi_low_rank_max_rank)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_type), DIMENSION(:), INTENT(INOUT)     :: mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(preconditioner_p_type), DIMENSION(:), POINTER :: ot_preconditioner
      INTEGER, INTENT(IN)                                :: prec_type, solver_type
      REAL(dp), INTENT(IN)                               :: energy_gap
      INTEGER, INTENT(IN)                                :: nspins
      LOGICAL, INTENT(IN), OPTIONAL                      :: has_unit_metric, convert_to_dbcsr
      INTEGER, INTENT(IN), OPTIONAL                      :: chol_type
      LOGICAL, INTENT(IN), OPTIONAL                      :: full_mo_set
      INTEGER, INTENT(IN), OPTIONAL                      :: chebyshev_degree, low_rank_base, &
                                                            fermi_low_rank_max_rank

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_preconditioner'

      CHARACTER(LEN=default_string_length)               :: msg
      INTEGER                                            :: handle, icall, ispin, &
                                                            my_chebyshev_degree, my_low_rank_base, &
                                                            my_fermi_low_rank_max_rank, n_loops
      INTEGER, DIMENSION(5)                              :: nocc, norb
      INTEGER, DIMENSION(:), POINTER                     :: multiple_unit_cell
      LOGICAL :: do_co_rotate, do_kpoints, my_convert_to_dbcsr, my_full_mo_set, &
         my_has_unit_metric, use_covariant_full_all, use_mo_coeff_b
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: kinetic
      TYPE(dbcsr_type), POINTER                          :: matrix_t, mo_coeff_b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: canonical_mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: input, subsys_section

      CALL timeset(routineN, handle)
      NULLIFY (matrix_t, mo_coeff_b, mo_coeff, kinetic, dft_control, para_env, blacs_env, &
               input, subsys_section, multiple_unit_cell)
      my_has_unit_metric = .FALSE.
      IF (PRESENT(has_unit_metric)) my_has_unit_metric = has_unit_metric
      my_convert_to_dbcsr = .TRUE.
      IF (PRESENT(convert_to_dbcsr)) my_convert_to_dbcsr = convert_to_dbcsr
      my_full_mo_set = .FALSE.
      IF (PRESENT(full_mo_set)) my_full_mo_set = full_mo_set
      my_chebyshev_degree = 8
      IF (PRESENT(chebyshev_degree)) my_chebyshev_degree = chebyshev_degree
      my_low_rank_base = ot_low_rank_base_overlap
      IF (PRESENT(low_rank_base)) my_low_rank_base = low_rank_base
      my_fermi_low_rank_max_rank = 48
      IF (PRESENT(fermi_low_rank_max_rank)) my_fermi_low_rank_max_rank = fermi_low_rank_max_rank
      do_kpoints = .FALSE.
      use_covariant_full_all = prec_type == ot_precond_full_all_covariant

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      input=input, &
                      para_env=para_env, &
                      blacs_env=blacs_env)

      IF (my_low_rank_base == ot_low_rank_base_lattice_fft) THEN
         IF (prec_type /= ot_precond_fermi_low_rank) THEN
            CPABORT("FERMI_LOW_RANK_BASE LATTICE_FFT requires PRECONDITIONER FERMI_LOW_RANK")
         END IF
         IF (my_has_unit_metric) THEN
            CPABORT("FERMI_LOW_RANK_BASE LATTICE_FFT is unavailable with a unit overlap metric")
         END IF
         IF (do_kpoints) THEN
            CPABORT("FERMI_LOW_RANK_BASE LATTICE_FFT currently supports only real Gamma-point OT")
         END IF
         IF (.NOT. ASSOCIATED(input)) THEN
            CPABORT("FERMI_LOW_RANK_BASE LATTICE_FFT could not access the input tree")
         END IF
         subsys_section => section_vals_get_subs_vals(input, "SUBSYS")
         IF (.NOT. ASSOCIATED(subsys_section)) THEN
            CPABORT("FERMI_LOW_RANK_BASE LATTICE_FFT requires a SUBSYS section")
         END IF
         CALL section_vals_val_get(subsys_section, "CELL%MULTIPLE_UNIT_CELL", &
                                   i_vals=multiple_unit_cell)
         IF (.NOT. ASSOCIATED(multiple_unit_cell)) THEN
            CPABORT("FERMI_LOW_RANK_BASE LATTICE_FFT could not read MULTIPLE_UNIT_CELL")
         END IF
         IF (SIZE(multiple_unit_cell) /= 3 .OR. ANY(multiple_unit_cell <= 0) .OR. &
             PRODUCT(multiple_unit_cell) <= 1) THEN
            CPABORT("FERMI_LOW_RANK_BASE LATTICE_FFT requires an explicit replicated multiple unit cell")
         END IF
         my_convert_to_dbcsr = .FALSE.
      ELSE IF (my_low_rank_base /= ot_low_rank_base_overlap) THEN
         CPABORT("Unknown FERMI_LOW_RANK_BASE setting")
      END IF

      IF (dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb .OR. &
          dft_control%qs_control%xtb) THEN
         IF (prec_type == ot_precond_full_kinetic) THEN
            msg = "Full_kinetic not available for semi-empirical methods"
            CPABORT(TRIM(msg))
         END IF
         matrix_t => matrix_s(1)%matrix
      ELSE
         CPASSERT(.NOT. my_has_unit_metric)
         CALL get_qs_env(qs_env, kinetic=kinetic)
         matrix_t => kinetic(1)%matrix
      END IF

      ! use full set of MOs or just occupied MOs
      nocc = 0
      norb = 0
      IF (my_full_mo_set) THEN
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin), homo=nocc(ispin), nmo=norb(ispin))
            CALL set_mo_set(mo_set=mos(ispin), homo=norb(ispin))
         END DO
      END IF
      !determines how often make preconditioner is called, spin dependent methods have to be called twice
      n_loops = 1
      IF (prec_type == ot_precond_full_single_inverse) n_loops = nspins
      ! check whether we need the ev and rotate the MOs
      SELECT CASE (prec_type)
      CASE (ot_precond_full_all)
         ! FULL_ALL needs canonical occupied eigenstates and one preconditioner per spin.
         n_loops = nspins

         do_co_rotate = ASSOCIATED(qs_env%mo_derivs)
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin), mo_coeff_b=mo_coeff_b, mo_coeff=mo_coeff)
            use_mo_coeff_b = mos(ispin)%use_mo_coeff_b
            IF (use_mo_coeff_b .AND. do_co_rotate) THEN
               CALL calculate_subspace_eigenvalues(mo_coeff_b, matrix_ks(ispin)%matrix, &
                                                   do_rotation=.TRUE., &
                                                   co_rotate=qs_env%mo_derivs(ispin)%matrix, &
                                                   para_env=para_env, &
                                                   blacs_env=blacs_env)
            ELSE IF (use_mo_coeff_b) THEN
               CALL calculate_subspace_eigenvalues(mo_coeff_b, matrix_ks(ispin)%matrix, &
                                                   do_rotation=.TRUE., &
                                                   para_env=para_env, &
                                                   blacs_env=blacs_env)
            ELSE
               CALL calculate_subspace_eigenvalues(mo_coeff, matrix_ks(ispin)%matrix, &
                                                   do_rotation=.TRUE.)
            END IF
         END DO
      CASE (ot_precond_fermi_low_rank)
         ! The common-reference low-rank operator depends only on the occupied projector.
         ! Keep the input gauge unchanged so that OT ROTATION remains well defined.
         n_loops = nspins
      CASE (ot_precond_full_all_covariant)
         ! Diagonalize a private occupied copy. The physical MOs and their
         ! derivatives remain in the user-selected gauge.
         n_loops = nspins
         ALLOCATE (canonical_mos(nspins))
         DO ispin = 1, nspins
            CALL duplicate_mo_set(canonical_mos(ispin), mos(ispin))
            CALL get_mo_set(mo_set=canonical_mos(ispin), mo_coeff_b=mo_coeff_b, mo_coeff=mo_coeff)
            IF (canonical_mos(ispin)%use_mo_coeff_b) THEN
               CALL copy_dbcsr_to_fm(mo_coeff_b, mo_coeff)
               canonical_mos(ispin)%use_mo_coeff_b = .FALSE.
            END IF
            CALL canonicalize_full_all_copy(ot_preconditioner(ispin)%preconditioner, &
                                            mo_coeff, matrix_ks(ispin)%matrix, &
                                            canonical_mos(ispin)%homo)
         END DO
      CASE DEFAULT
         ! No need to rotate the MOs
      END SELECT

      ! check whether we have a preconditioner
      SELECT CASE (prec_type)
      CASE (ot_precond_none)
         DO ispin = 1, SIZE(ot_preconditioner)
            ot_preconditioner(ispin)%preconditioner%in_use = 0
         END DO
      CASE DEFAULT
         DO icall = 1, n_loops
            ot_preconditioner(icall)%preconditioner%polynomial_degree = my_chebyshev_degree
            ot_preconditioner(icall)%preconditioner%low_rank_base = my_low_rank_base
            IF (my_low_rank_base == ot_low_rank_base_lattice_fft) THEN
               ot_preconditioner(icall)%preconditioner%lattice_dims = multiple_unit_cell
            END IF
            ot_preconditioner(icall)%preconditioner%spectral_max_rank = my_fermi_low_rank_max_rank
            IF (my_has_unit_metric) THEN
               IF (use_covariant_full_all) THEN
                  CALL make_preconditioner(ot_preconditioner(icall)%preconditioner, &
                                           prec_type, solver_type, &
                                           matrix_h=matrix_ks(icall)%matrix, &
                                           mo_set=canonical_mos(icall), &
                                           energy_gap=energy_gap, &
                                           convert_precond_to_dbcsr=my_convert_to_dbcsr)
               ELSE
                  CALL make_preconditioner(ot_preconditioner(icall)%preconditioner, &
                                           prec_type, solver_type, &
                                           matrix_h=matrix_ks(icall)%matrix, &
                                           mo_set=mos(icall), &
                                           energy_gap=energy_gap, &
                                           convert_precond_to_dbcsr=my_convert_to_dbcsr)
               END IF
            ELSE
               IF (use_covariant_full_all) THEN
                  CALL make_preconditioner(ot_preconditioner(icall)%preconditioner, &
                                           prec_type, solver_type, &
                                           matrix_h=matrix_ks(icall)%matrix, &
                                           matrix_s=matrix_s(1)%matrix, matrix_t=matrix_t, &
                                           mo_set=canonical_mos(icall), &
                                           energy_gap=energy_gap, &
                                           convert_precond_to_dbcsr=my_convert_to_dbcsr, chol_type=chol_type)
               ELSE
                  CALL make_preconditioner(ot_preconditioner(icall)%preconditioner, &
                                           prec_type, solver_type, &
                                           matrix_h=matrix_ks(icall)%matrix, &
                                           matrix_s=matrix_s(1)%matrix, matrix_t=matrix_t, &
                                           mo_set=mos(icall), energy_gap=energy_gap, &
                                           convert_precond_to_dbcsr=my_convert_to_dbcsr, chol_type=chol_type)
               END IF
            END IF
         END DO
      END SELECT

      IF (ALLOCATED(canonical_mos)) THEN
         DO ispin = 1, SIZE(canonical_mos)
            CALL deallocate_mo_set(canonical_mos(ispin))
         END DO
         DEALLOCATE (canonical_mos)
      END IF

      ! reset homo values
      IF (my_full_mo_set) THEN
         DO ispin = 1, nspins
            CALL set_mo_set(mo_set=mos(ispin), homo=nocc(ispin))
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE prepare_preconditioner

! **************************************************************************************************
!> \brief Canonicalize an occupied copy and retain the occupied-space rotation.
!> \param preconditioner_env ...
!> \param mo_coeff ...
!> \param matrix_h ...
!> \param nocc ...
! **************************************************************************************************
   SUBROUTINE canonicalize_full_all_copy(preconditioner_env, mo_coeff, matrix_h, nocc)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_type), POINTER                          :: matrix_h
      INTEGER, INTENT(IN)                                :: nocc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'canonicalize_full_all_copy'

      INTEGER                                            :: handle, nao
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: occupied_evals
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: canonical_coeff, h_block, weighted_coeff

      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(mo_coeff, nrow_global=nao)
      CPASSERT(nocc > 0)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nocc, &
                               context=preconditioner_env%ctxt, &
                               para_env=preconditioner_env%para_env)
      CALL cp_fm_create(canonical_coeff, fm_struct, name=routineN)
      CALL cp_fm_create(weighted_coeff, fm_struct, name=routineN)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_to_fm(mo_coeff, canonical_coeff, nocc)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nocc, ncol_global=nocc, &
                               context=preconditioner_env%ctxt, &
                               para_env=preconditioner_env%para_env)
      CALL cp_fm_create(h_block, fm_struct, name=routineN)
      IF (ASSOCIATED(preconditioner_env%occ_rotation_dbcsr)) THEN
         CALL dbcsr_release_p(preconditioner_env%occ_rotation_dbcsr)
      END IF
      IF (ASSOCIATED(preconditioner_env%occ_rotation)) THEN
         CALL cp_fm_release(preconditioner_env%occ_rotation)
         DEALLOCATE (preconditioner_env%occ_rotation)
      END IF
      ALLOCATE (preconditioner_env%occ_rotation)
      CALL cp_fm_create(preconditioner_env%occ_rotation, fm_struct, name="occupied rotation")
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_dbcsr_sm_fm_multiply(matrix_h, canonical_coeff, weighted_coeff, nocc)
      CALL parallel_gemm('T', 'N', nocc, nocc, nao, 1.0_dp, canonical_coeff, &
                         weighted_coeff, 0.0_dp, h_block)
      ALLOCATE (occupied_evals(nocc))
      CALL choose_eigv_solver(h_block, preconditioner_env%occ_rotation, occupied_evals)

      CALL parallel_gemm('N', 'N', nao, nocc, nocc, 1.0_dp, canonical_coeff, &
                         preconditioner_env%occ_rotation, 0.0_dp, weighted_coeff)
      CALL cp_fm_to_fm(weighted_coeff, mo_coeff, nocc, 1, 1)

      DEALLOCATE (occupied_evals)
      CALL cp_fm_release(canonical_coeff)
      CALL cp_fm_release(h_block)
      CALL cp_fm_release(weighted_coeff)
      CALL timestop(handle)

   END SUBROUTINE canonicalize_full_all_copy

END MODULE preconditioner
