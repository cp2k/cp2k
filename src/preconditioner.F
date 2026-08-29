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
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE input_constants,                 ONLY: &
        cholesky_reduce, ot_precond_full_all, ot_precond_full_kinetic, ot_precond_full_single, &
        ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, &
        ot_precond_solver_default, ot_precond_solver_inv_chol, ot_precond_solver_update
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE preconditioner_apply,            ONLY: apply_preconditioner_cfm_complex,&
                                              apply_preconditioner_dbcsr,&
                                              apply_preconditioner_dbcsr_complex,&
                                              apply_preconditioner_fm
   USE preconditioner_makes,            ONLY: make_complex_full_all,&
                                              make_complex_full_kinetic,&
                                              make_complex_full_s_inverse,&
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
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type,&
                                              set_mo_set
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner'

   PUBLIC :: make_preconditioner, make_preconditioner_complex_full_all, &
             make_preconditioner_complex_full_kinetic, &
             make_preconditioner_complex_full_s_inverse, &
             make_preconditioner_complex_full_single_inverse, restart_preconditioner
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

      INTEGER                                            :: handle, k, my_solver_type, nao, nhomo
      LOGICAL                                            :: my_convert_precond_to_dbcsr, &
                                                            needs_full_spectrum, needs_homo, &
                                                            use_mo_coeff_b
      REAL(KIND=dp)                                      :: energy_homo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues_ot
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: mo_occ
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_b

      CALL timeset(routineN, handle)

      CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b, homo=nhomo)
      use_mo_coeff_b = mo_set%use_mo_coeff_b
      CALL cp_fm_get_info(mo_coeff, ncol_global=k, nrow_global=nao)

      ! Starting some matrix mess, check where to store the result in preconditioner_env, fm or dbcsr_matrix
      my_convert_precond_to_dbcsr = .FALSE.
      IF (PRESENT(convert_precond_to_dbcsr)) my_convert_precond_to_dbcsr = convert_precond_to_dbcsr

      ! Thanks to the mess with the matrices we need to make sure in this case that the
      ! Previous inverse is properly stored as a sparse matrix, fm gets deallocated here
      ! if it wasn't anyway
      IF (preconditioner_env%solver == ot_precond_solver_update) THEN
         CALL transfer_fm_to_dbcsr(preconditioner_env%fm, preconditioner_env%dbcsr_matrix, matrix_h)
      END IF

      needs_full_spectrum = .FALSE.
      needs_homo = .FALSE.

      SELECT CASE (precon_type)
      CASE (ot_precond_full_all)
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

      ! Here comes more matrix mess, make sure to output the correct matrix format,
      ! A bit pointless to convert the cholesky factorized version as it doesn't work in
      ! dbcsr form and will crash later,...
      IF (my_convert_precond_to_dbcsr) THEN
         CALL transfer_fm_to_dbcsr(preconditioner_env%fm, preconditioner_env%dbcsr_matrix, matrix_h)
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
!> \param matrix_h_re real part of H(k)
!> \param matrix_h_im imaginary part of H(k)
!> \param matrix_s_re real part of S(k)
!> \param matrix_s_im imaginary part of S(k)
!> \param mo_set_re real reference orbitals and their energy labels
!> \param mo_set_im imaginary reference orbitals
!> \param energy_gap denominator floor
!> \param solver_type inverse application selected in the OT input
! **************************************************************************************************
   SUBROUTINE make_preconditioner_complex_full_all(preconditioner_env, &
                                                   matrix_h_re, matrix_h_im, &
                                                   matrix_s_re, matrix_s_im, &
                                                   mo_set_re, mo_set_im, energy_gap, solver_type)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type), POINTER                          :: matrix_h_re, matrix_h_im, matrix_s_re, &
                                                            matrix_s_im
      TYPE(mo_set_type), INTENT(IN)                      :: mo_set_re, mo_set_im
      REAL(KIND=dp), INTENT(IN)                          :: energy_gap
      INTEGER, INTENT(IN)                                :: solver_type

      CHARACTER(len=*), PARAMETER :: routineN = 'make_preconditioner_complex_full_all'

      INTEGER                                            :: handle, nao, nmo
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cp_cfm_type)                                  :: matrix_c, matrix_h, matrix_s
      TYPE(cp_fm_type)                                   :: matrix_h_im_fm, matrix_h_re_fm, &
                                                            matrix_s_im_fm, matrix_s_re_fm
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_im, mo_coeff_re

      CALL timeset(routineN, handle)

      NULLIFY (eigenvalues, mo_coeff_im, mo_coeff_re)
      CPASSERT(ASSOCIATED(matrix_h_re))
      CPASSERT(ASSOCIATED(matrix_h_im))
      CPASSERT(ASSOCIATED(matrix_s_re))
      CPASSERT(ASSOCIATED(matrix_s_im))
      IF (solver_type /= ot_precond_solver_default) THEN
         CPABORT('Complex FULL_ALL supports only PRECOND_SOLVER DEFAULT')
      END IF
      CALL get_mo_set(mo_set_re, mo_coeff=mo_coeff_re, eigenvalues=eigenvalues, &
                      nao=nao, nmo=nmo)
      CALL get_mo_set(mo_set_im, mo_coeff=mo_coeff_im)
      CPASSERT(ASSOCIATED(mo_coeff_re))
      CPASSERT(ASSOCIATED(mo_coeff_im))
      CPASSERT(ASSOCIATED(eigenvalues))
      CPASSERT(SIZE(eigenvalues) >= nmo)

      CALL cp_fm_create(matrix_h_re_fm, mo_coeff_re%matrix_struct, nrow=nao, ncol=nao, &
                        name='complex FULL_ALL H real')
      CALL cp_fm_create(matrix_h_im_fm, mo_coeff_re%matrix_struct, nrow=nao, ncol=nao, &
                        name='complex FULL_ALL H imaginary')
      CALL cp_fm_create(matrix_s_re_fm, mo_coeff_re%matrix_struct, nrow=nao, ncol=nao, &
                        name='complex FULL_ALL S real')
      CALL cp_fm_create(matrix_s_im_fm, mo_coeff_re%matrix_struct, nrow=nao, ncol=nao, &
                        name='complex FULL_ALL S imaginary')
      CALL copy_dbcsr_to_fm(matrix_h_re, matrix_h_re_fm)
      CALL copy_dbcsr_to_fm(matrix_h_im, matrix_h_im_fm)
      CALL copy_dbcsr_to_fm(matrix_s_re, matrix_s_re_fm)
      CALL copy_dbcsr_to_fm(matrix_s_im, matrix_s_im_fm)

      CALL cp_cfm_create(matrix_h, matrix_h_re_fm%matrix_struct, name='complex FULL_ALL H')
      CALL cp_cfm_create(matrix_s, matrix_s_re_fm%matrix_struct, name='complex FULL_ALL S')
      CALL cp_cfm_create(matrix_c, mo_coeff_re%matrix_struct, name='complex FULL_ALL C')
      CALL cp_fm_to_cfm(matrix_h_re_fm, matrix_h_im_fm, matrix_h)
      CALL cp_fm_to_cfm(matrix_s_re_fm, matrix_s_im_fm, matrix_s)
      CALL cp_fm_to_cfm(mo_coeff_re, mo_coeff_im, matrix_c)
      CALL make_complex_full_all(preconditioner_env, matrix_c, matrix_h, matrix_s, &
                                 eigenvalues(1:nmo), energy_gap)

      CALL cp_cfm_release(matrix_c)
      CALL cp_cfm_release(matrix_s)
      CALL cp_cfm_release(matrix_h)
      CALL cp_fm_release(matrix_s_im_fm)
      CALL cp_fm_release(matrix_s_re_fm)
      CALL cp_fm_release(matrix_h_im_fm)
      CALL cp_fm_release(matrix_h_re_fm)

      CALL timestop(handle)

   END SUBROUTINE make_preconditioner_complex_full_all

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
         CASE (ot_precond_full_all, ot_precond_full_single) ! these depend on the ks matrix
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
         CASE (ot_precond_full_all, ot_precond_full_single_inverse)
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
! **************************************************************************************************
   SUBROUTINE prepare_preconditioner(qs_env, mos, matrix_ks, matrix_s, &
                                     ot_preconditioner, prec_type, solver_type, &
                                     energy_gap, nspins, has_unit_metric, &
                                     convert_to_dbcsr, chol_type, full_mo_set, chebyshev_degree)

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
      INTEGER, INTENT(IN), OPTIONAL                      :: chebyshev_degree

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_preconditioner'

      CHARACTER(LEN=default_string_length)               :: msg
      INTEGER                                            :: handle, icall, ispin, &
                                                            my_chebyshev_degree, n_loops
      INTEGER, DIMENSION(5)                              :: nocc, norb
      LOGICAL                                            :: do_co_rotate, my_convert_to_dbcsr, &
                                                            my_full_mo_set, my_has_unit_metric, &
                                                            use_mo_coeff_b
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: kinetic
      TYPE(dbcsr_type), POINTER                          :: matrix_t, mo_coeff_b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      NULLIFY (matrix_t, mo_coeff_b, mo_coeff, kinetic, dft_control, para_env, blacs_env)
      my_has_unit_metric = .FALSE.
      IF (PRESENT(has_unit_metric)) my_has_unit_metric = has_unit_metric
      my_convert_to_dbcsr = .TRUE.
      IF (PRESENT(convert_to_dbcsr)) my_convert_to_dbcsr = convert_to_dbcsr
      my_full_mo_set = .FALSE.
      IF (PRESENT(full_mo_set)) my_full_mo_set = full_mo_set
      my_chebyshev_degree = 8
      IF (PRESENT(chebyshev_degree)) my_chebyshev_degree = chebyshev_degree

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      blacs_env=blacs_env)

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
         ! if one of these preconditioners is used every spin needs to call make_preconditioner
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
            IF (my_has_unit_metric) THEN
               CALL make_preconditioner(ot_preconditioner(icall)%preconditioner, &
                                        prec_type, &
                                        solver_type, &
                                        matrix_h=matrix_ks(icall)%matrix, &
                                        mo_set=mos(icall), &
                                        energy_gap=energy_gap, &
                                        convert_precond_to_dbcsr=my_convert_to_dbcsr)
            ELSE
               CALL make_preconditioner(ot_preconditioner(icall)%preconditioner, &
                                        prec_type, &
                                        solver_type, &
                                        matrix_h=matrix_ks(icall)%matrix, &
                                        matrix_s=matrix_s(1)%matrix, &
                                        matrix_t=matrix_t, &
                                        mo_set=mos(icall), &
                                        energy_gap=energy_gap, &
                                        convert_precond_to_dbcsr=my_convert_to_dbcsr, chol_type=chol_type)
            END IF
         END DO
      END SELECT

      ! reset homo values
      IF (my_full_mo_set) THEN
         DO ispin = 1, nspins
            CALL set_mo_set(mo_set=mos(ispin), homo=nocc(ispin))
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE prepare_preconditioner

END MODULE preconditioner
