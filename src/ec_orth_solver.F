!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief AO-based conjugate-graadient response solver routines
!>
!>
!> \date 09.2019
!> \author Fabian Belleflamme
! **************************************************************************************************
MODULE ec_orth_solver
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_external_control,             ONLY: external_control
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_checksum, dbcsr_copy, dbcsr_create, &
        dbcsr_desymmetrize, dbcsr_dot, dbcsr_filter, dbcsr_finalize, dbcsr_get_info, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_scale, dbcsr_transposed, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE ec_methods,                      ONLY: create_kernel
   USE input_constants,                 ONLY: kg_tnadd_embed,&
                                              kg_tnadd_embed_ri,&
                                              ls_s_sqrt_ns,&
                                              ls_s_sqrt_proot,&
                                              precond_mlp
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz,&
                                              matrix_sqrt_proot
   USE kg_correction,                   ONLY: kg_ekin_subset
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathlib,                         ONLY: abnormal_value
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_linres_types,                 ONLY: linres_control_type
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild,&
                                              qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_orth_solver'

   ! Public subroutines

   PUBLIC :: ec_response_ao

CONTAINS

! **************************************************************************************************
!> \brief      Preconditioning of the AO-based CG linear response solver
!>             M * z_0 = r_0
!>             M(X) = [F,B], with B = [X,P]
!>             for M we need F and P in ortho basis
!>             Returns z_0, the preconditioned residual in orthonormal basis
!>
!>             All matrices are in orthonormal Lowdin basis
!>
!> \param qs_env ...
!> \param matrix_ks Ground-state Kohn-Sham matrix
!> \param matrix_p  Ground-state Density matrix
!> \param matrix_rhs Unpreconditioned residual of linear response CG
!> \param matrix_cg_z Preconditioned residual
!> \param eps_filter ...
!> \param iounit ...
!>
!> \date    01.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE preconditioner(qs_env, matrix_ks, matrix_p, matrix_rhs, &
                             matrix_cg_z, eps_filter, iounit)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_ks, matrix_p, matrix_rhs
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: matrix_cg_z
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER, INTENT(IN)                                :: iounit

      CHARACTER(len=*), PARAMETER :: routineN = 'preconditioner', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ispin, max_iter, nao, nspins
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: norm_res, t1, t2
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha, beta, new_norm, norm_cA, norm_rr
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_Ax, matrix_b, matrix_cg, &
                                                            matrix_res
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_p))
      CPASSERT(ASSOCIATED(matrix_rhs))
      CPASSERT(ASSOCIATED(matrix_cg_z))

      NULLIFY (dft_control, linres_control)

      t1 = m_walltime()

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      linres_control=linres_control)
      nspins = dft_control%nspins
      CALL dbcsr_get_info(matrix_ks(1)%matrix, nfullrows_total=nao)

      ALLOCATE (alpha(nspins), beta(nspins), new_norm(nspins), norm_cA(nspins), norm_rr(nspins))

      !----------------------------------------
      ! Create non-symmetric matrices: Ax, B, cg, res
      !----------------------------------------

      NULLIFY (matrix_Ax, matrix_b, matrix_cg, matrix_res)
      CALL dbcsr_allocate_matrix_set(matrix_Ax, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_b, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_cg, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_res, nspins)

      DO ispin = 1, nspins
         ALLOCATE (matrix_Ax(ispin)%matrix)
         ALLOCATE (matrix_b(ispin)%matrix)
         ALLOCATE (matrix_cg(ispin)%matrix)
         ALLOCATE (matrix_res(ispin)%matrix)
         CALL dbcsr_create(matrix_Ax(ispin)%matrix, name="linop MATRIX", &
                           template=matrix_ks(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_b(ispin)%matrix, name="MATRIX B", &
                           template=matrix_ks(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_cg(ispin)%matrix, name="TRIAL MATRIX", &
                           template=matrix_ks(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_res(ispin)%matrix, name="RESIDUE", &
                           template=matrix_ks(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
      END DO

      !----------------------------------------
      ! Get righ-hand-side operators
      !----------------------------------------

      ! Initial guess z_0
      DO ispin = 1, nspins
         CALL dbcsr_copy(matrix_cg_z(ispin)%matrix, matrix_rhs(ispin)%matrix)

         ! r_0 = b
         CALL dbcsr_copy(matrix_res(ispin)%matrix, matrix_rhs(ispin)%matrix)
      END DO

      ! Projector on trial matrix
      ! Projector does not need to be applied here,
      ! as matrix_rhs already had this done before entering preconditioner
      !CALL projector(qs_env, matrix_p, matrix_cg_z, eps_filter)

      ! Mz_0
      CALL hessian_op1(matrix_ks, matrix_p, matrix_cg_z, matrix_b, matrix_Ax, eps_filter)

      ! r_0 = b - Ax_0
      DO ispin = 1, nspins
         CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -1.0_dp)
      END DO

      ! Matrix projector T
      CALL projector(qs_env, matrix_p, matrix_res, eps_filter)

      DO ispin = 1, nspins
         ! cg = p_0 = z_0
         CALL dbcsr_copy(matrix_cg(ispin)%matrix, matrix_res(ispin)%matrix)
      END DO

      ! header
      IF (iounit > 0) THEN
         WRITE (iounit, "(/,T10,A)") "Preconditioning of search direction"
         WRITE (iounit, "(/,T10,A,T25,A,T42,A,T62,A,/,T10,A)") &
            "Iteration", "Stepsize", "Convergence", "Time", &
            REPEAT("-", 58)
      END IF

      alpha(:) = 0.0_dp
      max_iter = 200
      converged = .FALSE.
      norm_res = 0.0_dp

      ! start iteration
      iteration: DO i = 1, max_iter

         ! Hessian Ax = [F,B] is updated preconditioner
         CALL hessian_op1(matrix_ks, matrix_p, matrix_cg, matrix_b, matrix_Ax, eps_filter)

         ! Matrix projector
         CALL projector(qs_env, matrix_p, matrix_Ax, eps_filter)

         DO ispin = 1, nspins

            ! Tr(r_0 * r_0)
            CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_res(ispin)%matrix, norm_rr(ispin))
            IF (abnormal_value(norm_rr(ispin))) &
               CPABORT("Preconditioner: Tr[r_j*r_j] is an abnormal value (NaN/Inf)")

            IF (norm_rr(ispin) .LT. 0.0_dp) CPABORT("norm_rr < 0")
            norm_res = MAX(norm_res, ABS(norm_rr(ispin)/REAL(nao, dp)))

            ! norm_cA = tr(Ap_j * p_j)
            CALL dbcsr_dot(matrix_cg(ispin)%matrix, matrix_Ax(ispin)%matrix, norm_cA(ispin))

            ! Determine step-size
            IF (norm_cA(ispin) .LT. linres_control%eps) THEN
               alpha(ispin) = 1.0_dp
            ELSE
               alpha(ispin) = norm_rr(ispin)/norm_cA(ispin)
            END IF

            ! x_j+1 = x_j + alpha*p_j
            ! save contribution of this iteration
            CALL dbcsr_add(matrix_cg_z(ispin)%matrix, matrix_cg(ispin)%matrix, 1.0_dp, alpha(ispin))

            ! r_j+1 = r_j - alpha * Ap_j
            CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -alpha(ispin))

         END DO

         norm_res = 0.0_dp

         DO ispin = 1, nspins
            ! Tr[r_j+1*z_j+1]
            CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_res(ispin)%matrix, new_norm(ispin))
            IF (new_norm(ispin) .LT. 0.0_dp) CPABORT("tr(r_j+1*z_j+1) < 0")
            IF (abnormal_value(new_norm(ispin))) &
               CPABORT("Preconditioner: Tr[r_j+1*z_j+1] is an abnormal value (NaN/Inf)")
            norm_res = MAX(norm_res, new_norm(ispin)/REAL(nao, dp))

            IF (norm_rr(ispin) .LT. linres_control%eps*0.001_dp &
                .OR. new_norm(ispin) .LT. linres_control%eps*0.001_dp) THEN
               beta(ispin) = 0.0_dp
               converged = .TRUE.
            ELSE
               beta(ispin) = new_norm(ispin)/norm_rr(ispin)
            END IF

            ! update new search vector (matrix cg)
            ! cg_j+1 = z_j+1 + beta*cg_j
            CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_res(ispin)%matrix, beta(ispin), 1.0_dp)
            CALL dbcsr_filter(matrix_cg(ispin)%matrix, eps_filter)

            norm_rr(ispin) = new_norm(ispin)
         END DO

         ! Convergence criteria
         IF (norm_res .LT. linres_control%eps) THEN
            converged = .TRUE.
         END IF

         t2 = m_walltime()
         IF (i .EQ. 1 .OR. MOD(i, 1) .EQ. 0 .OR. converged) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(T10,I5,T25,1E8.2,T33,F25.14,T58,F8.2)") &
                  i, MAXVAL(alpha), norm_res, t2 - t1
               ! Convergence in scientific notation
               !WRITE (iounit, "(T10,I5,T25,1E8.2,T42,1E14.8,T58,F8.2)") &
               !   i, MAXVAL(alpha), norm_res, t2 - t1
               CALL m_flush(iounit)
            END IF
         END IF
         IF (converged) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T10,A,I4,A,/)") "The precon solver converged in ", i, " iterations."
               CALL m_flush(iounit)
            END IF
            EXIT iteration
         END IF

         ! Max number of iteration reached
         IF (i == max_iter) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T10,A/)") &
                  "The precon solver didnt converge! Maximum number of iterations reached."
               CALL m_flush(iounit)
            END IF
            converged = .FALSE.
         END IF

      END DO iteration

      ! Matrix projector
      CALL projector(qs_env, matrix_p, matrix_cg_z, eps_filter)

      ! Release matrices
      CALL dbcsr_deallocate_matrix_set(matrix_Ax)
      CALL dbcsr_deallocate_matrix_set(matrix_b)
      CALL dbcsr_deallocate_matrix_set(matrix_res)
      CALL dbcsr_deallocate_matrix_set(matrix_cg)

      DEALLOCATE (alpha, beta, new_norm, norm_cA, norm_rr)

      CALL timestop(handle)

   END SUBROUTINE preconditioner

! **************************************************************************************************
!> \brief AO-based conjugate radient linear response solver.
!>        In goes the right hand side B of the equation AZ=B, and the linear transformation of the
!>        Hessian matrix A on trial matries is iteratively solved. Result are
!>        the response density matrix_pz, and the energy-weighted response density matrix_wz.
!>
!> \param qs_env ...
!> \param matrix_hz Right hand-side of linear response equation
!> \param matrix_pz Response density
!> \param matrix_wz Energy-weighted response density matrix
!> \param iounit ...
!> \param should_stop ...
!>
!> \date    01.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_response_ao(qs_env, matrix_hz, matrix_pz, matrix_wz, iounit, should_stop)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_hz
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: matrix_pz, matrix_wz
      INTEGER, INTENT(IN)                                :: iounit
      LOGICAL, INTENT(OUT)                               :: should_stop

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_response_ao', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ispin, max_iter_lanczos, nao, &
                                                            nspins, s_sqrt_method, s_sqrt_order
      LOGICAL                                            :: restart
      REAL(KIND=dp)                                      :: eps_filter, eps_lanczos, focc, &
                                                            min_shift, norm_res, old_conv, shift, &
                                                            t1, t2
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha, beta, new_norm, norm_cA, norm_rr, &
                                                            tr_rz00
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: ksmat, matrix_Ax, matrix_cg, matrix_cg_z, &
         matrix_ks, matrix_nsc, matrix_p, matrix_res, matrix_s, matrix_z, matrix_z0, rho_ao
      TYPE(dbcsr_type)                                   :: matrix_s_sqrt, matrix_s_sqrt_inv, &
                                                            matrix_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: solver_section

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(matrix_hz))
      CPASSERT(ASSOCIATED(matrix_pz))
      CPASSERT(ASSOCIATED(matrix_wz))

      NULLIFY (dft_control, ksmat, matrix_s, linres_control, rho)

      t1 = m_walltime()

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      linres_control=linres_control, &
                      matrix_ks=ksmat, &
                      matrix_s=matrix_s, &
                      rho=rho)
      nspins = dft_control%nspins

      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

      solver_section => section_vals_get_subs_vals(qs_env%input, "DFT%ENERGY_CORRECTION%RESPONSE_SOLVER")
      CALL section_vals_val_get(solver_section, "S_SQRT_METHOD", i_val=s_sqrt_method)
      CALL section_vals_val_get(solver_section, "S_SQRT_ORDER", i_val=s_sqrt_order)
      CALL section_vals_val_get(solver_section, "EPS_LANCZOS", r_val=eps_lanczos)
      CALL section_vals_val_get(solver_section, "MAX_ITER_LANCZOS", i_val=max_iter_lanczos)

      eps_filter = linres_control%eps_filter

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ALLOCATE (alpha(nspins), beta(nspins), new_norm(nspins), norm_cA(nspins), norm_rr(nspins))
      ALLOCATE (tr_rz00(nspins))

      ! local matrix P, KS, and NSC
      ! to bring into orthogonal basis
      NULLIFY (matrix_p, matrix_ks, matrix_nsc)
      CALL dbcsr_allocate_matrix_set(matrix_p, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_ks, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_nsc, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_p(ispin)%matrix)
         ALLOCATE (matrix_ks(ispin)%matrix)
         ALLOCATE (matrix_nsc(ispin)%matrix)
         CALL dbcsr_create(matrix_p(ispin)%matrix, name="P_IN ORTHO", &
                           template=ksmat(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_ks(ispin)%matrix, name="KS_IN ORTHO", &
                           template=ksmat(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_nsc(ispin)%matrix, name="NSC IN ORTHO", &
                           template=ksmat(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_desymmetrize(rho_ao(ispin)%matrix, matrix_p(ispin)%matrix)
         CALL dbcsr_desymmetrize(ksmat(ispin)%matrix, matrix_ks(ispin)%matrix)
         CALL dbcsr_desymmetrize(matrix_hz(ispin)%matrix, matrix_nsc(ispin)%matrix)
      END DO

      ! Scale matrix_p by factor 1/2 in closed-shell
      IF (nspins == 1) CALL dbcsr_scale(matrix_p(1)%matrix, 0.5_dp)

      ! Transform P, KS, and Harris kernel matrix into Orthonormal basis
      CALL dbcsr_create(matrix_s_sqrt, template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_s_sqrt_inv, template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      SELECT CASE (s_sqrt_method)
      CASE (ls_s_sqrt_proot)
         CALL matrix_sqrt_proot(matrix_s_sqrt, matrix_s_sqrt_inv, &
                                matrix_s(1)%matrix, eps_filter, &
                                s_sqrt_order, eps_lanczos, max_iter_lanczos, symmetrize=.TRUE.)
      CASE (ls_s_sqrt_ns)
         CALL matrix_sqrt_Newton_Schulz(matrix_s_sqrt, matrix_s_sqrt_inv, &
                                        matrix_s(1)%matrix, eps_filter, &
                                        s_sqrt_order, eps_lanczos, max_iter_lanczos)
      CASE DEFAULT
         CPABORT("Unknown sqrt method.")
      END SELECT

      ! Transform into orthonormal Lowdin basis
      DO ispin = 1, nspins
         CALL transform_m_orth(matrix_p(ispin)%matrix, matrix_s_sqrt, eps_filter)
         CALL transform_m_orth(matrix_ks(ispin)%matrix, matrix_s_sqrt_inv, eps_filter)
         CALL transform_m_orth(matrix_nsc(ispin)%matrix, matrix_s_sqrt_inv, eps_filter)
      END DO

      !----------------------------------------
      ! Create non-symmetric work matrices: Ax, cg, res
      ! Content of Ax, cg, cg_z, res, z0 anti-symmetric
      ! Content of z symmetric
      !----------------------------------------

      CALL dbcsr_create(matrix_tmp, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)

      NULLIFY (matrix_Ax, matrix_cg, matrix_cg_z, matrix_res, matrix_z, matrix_z0)
      CALL dbcsr_allocate_matrix_set(matrix_Ax, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_cg, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_cg_z, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_res, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_z, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_z0, nspins)

      DO ispin = 1, nspins
         ALLOCATE (matrix_Ax(ispin)%matrix)
         ALLOCATE (matrix_cg(ispin)%matrix)
         ALLOCATE (matrix_cg_z(ispin)%matrix)
         ALLOCATE (matrix_res(ispin)%matrix)
         ALLOCATE (matrix_z(ispin)%matrix)
         ALLOCATE (matrix_z0(ispin)%matrix)
         CALL dbcsr_create(matrix_Ax(ispin)%matrix, name="linop MATRIX", &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_cg(ispin)%matrix, name="TRIAL MATRIX", &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_cg_z(ispin)%matrix, name="MATRIX CG-Z", &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_res(ispin)%matrix, name="RESIDUE", &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_z(ispin)%matrix, name="Z-Matrix", &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_z0(ispin)%matrix, name="p after precondi-Matrix", &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
      END DO

      !----------------------------------------
      ! Get righ-hand-side operators
      !----------------------------------------

      ! Spin factor
      focc = -2.0_dp
      IF (nspins == 1) focc = -4.0_dp

      ! E^[1]_Harris = -4*G[\delta P]*Pin - Pin*G[\delta P] = -4*[G[\delta P], Pin]
      CALL commutator(matrix_nsc, matrix_p, matrix_res, eps_filter, .FALSE., alpha=focc)

      ! Initial guess cg_Z
      DO ispin = 1, nspins
         CALL dbcsr_copy(matrix_cg_z(ispin)%matrix, matrix_res(ispin)%matrix)
      END DO

      ! Projector on trial matrix
      CALL projector(qs_env, matrix_p, matrix_cg_z, eps_filter)

      ! Ax0
      CALL build_hessian_op(qs_env=qs_env, &
                            matrix_ks=matrix_ks, &
                            matrix_p=matrix_p, &   ! p
                            matrix_s_sqrt_inv=matrix_s_sqrt_inv, &
                            matrix_cg=matrix_cg_z, & ! cg
                            matrix_Ax=matrix_Ax, &
                            eps_filter=eps_filter)

      ! r_0 = b - Ax0
      DO ispin = 1, nspins
         CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -1.0_dp)
      END DO

      ! Matrix projector T
      CALL projector(qs_env, matrix_p, matrix_res, eps_filter)

      ! Preconditioner
      linres_control%flag = ""
      IF (linres_control%preconditioner_type == precond_mlp) THEN
         ! M * z_0 = r_0
         ! Conjugate gradient returns z_0
         CALL preconditioner(qs_env=qs_env, &
                             matrix_ks=matrix_ks, &
                             matrix_p=matrix_p, &
                             matrix_rhs=matrix_res, &
                             matrix_cg_z=matrix_z0, &
                             eps_filter=eps_filter, &
                             iounit=iounit)
         linres_control%flag = "PCG-AO"
      ELSE
         ! z_0 = r_0
         DO ispin = 1, nspins
            CALL dbcsr_copy(matrix_z0(ispin)%matrix, matrix_res(ispin)%matrix)
            linres_control%flag = "CG-AO"
         END DO
      END IF

      norm_res = 0.0_dp

      DO ispin = 1, nspins
         ! cg = p_0 = z_0
         CALL dbcsr_copy(matrix_cg(ispin)%matrix, matrix_z0(ispin)%matrix)

         ! Tr(r_0 * z_0)
         CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_cg(ispin)%matrix, norm_rr(ispin))

         IF (norm_rr(ispin) .LT. 0.0_dp) CPABORT("norm_rr < 0")
         norm_res = MAX(norm_res, ABS(norm_rr(ispin)/REAL(nao, dp)))
      END DO

      ! eigenvalue shifting
      min_shift = 0.0_dp
      old_conv = norm_rr(1)
      shift = MIN(10.0_dp, MAX(min_shift, 0.05_dp*old_conv))
      old_conv = 100.0_dp

      ! header
      IF (iounit > 0) THEN
         WRITE (iounit, "(/,T3,A,T16,A,T25,A,T38,A,T52,A,/,T3,A)") &
            "Iteration", "Method", "Stepsize", "Convergence", "Time", &
            REPEAT("-", 80)
      END IF

      alpha(:) = 0.0_dp
      restart = .FALSE.
      should_stop = .FALSE.
      linres_control%converged = .FALSE.

      ! start iteration
      iteration: DO i = 1, linres_control%max_iter

         ! Convergence criteria
         ! default for eps 10E-6 in MO_linres
         IF (norm_res .LT. linres_control%eps) THEN
            linres_control%converged = .TRUE.
         END IF

         t2 = m_walltime()
         IF (i .EQ. 1 .OR. MOD(i, 1) .EQ. 0 .OR. linres_control%converged &
             .OR. restart .OR. should_stop) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(T5,I5,T18,A3,T28,L1,T38,1E8.2,T48,F16.10,T68,F8.2)") &
                  i, linres_control%flag, restart, MAXVAL(alpha), norm_res, t2 - t1
               CALL m_flush(iounit)
            END IF
         END IF
         IF (linres_control%converged) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T2,A,I4,A,/)") "The linear solver converged in ", i, " iterations."
               CALL m_flush(iounit)
            END IF
            EXIT iteration
         ELSE IF (should_stop) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T2,A,I4,A,/)") "The linear solver did NOT converge! External stop"
               CALL m_flush(iounit)
            END IF
            EXIT iteration
         END IF

         ! Max number of iteration reached
         IF (i == linres_control%max_iter) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, "(/,T2,A/)") &
                  "The linear solver didnt converge! Maximum number of iterations reached."
               CALL m_flush(iounit)
            END IF
            linres_control%converged = .FALSE.
         END IF

         ! Hessian Ax = [F,B] + [G(B),P]
         CALL build_hessian_op(qs_env=qs_env, &
                               matrix_ks=matrix_ks, &
                               matrix_p=matrix_p, &   ! p
                               matrix_s_sqrt_inv=matrix_s_sqrt_inv, &
                               matrix_cg=matrix_cg, & ! cg
                               matrix_Ax=matrix_Ax, &
                               eps_filter=eps_filter)

         ! Matrix projector T
         CALL projector(qs_env, matrix_p, matrix_Ax, eps_filter)

         DO ispin = 1, nspins

            CALL dbcsr_filter(matrix_Ax(ispin)%matrix, eps_filter)
            ! norm_cA = tr(Ap_j * p_j)
            CALL dbcsr_dot(matrix_cg(ispin)%matrix, matrix_Ax(ispin)%matrix, norm_cA(ispin))

            IF (norm_cA(ispin) .LT. 0.0_dp) THEN

               ! Recalculate w/o preconditioner
               IF (i > 1) THEN
                  ! p = -z + beta*p
                  CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_z0(ispin)%matrix, &
                                 beta(ispin), -1.0_dp)
                  CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_res(ispin)%matrix, new_norm(ispin))
                  beta(ispin) = new_norm(ispin)/tr_rz00(ispin)
                  CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_res(ispin)%matrix, &
                                 beta(ispin), 1.0_dp)
                  norm_rr(ispin) = new_norm(ispin)
               ELSE
                  CALL dbcsr_copy(matrix_res(ispin)%matrix, matrix_cg(ispin)%matrix)
                  CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_res(ispin)%matrix, norm_rr(ispin))
               END IF

               CALL build_hessian_op(qs_env=qs_env, &
                                     matrix_ks=matrix_ks, &
                                     matrix_p=matrix_p, &   ! p
                                     matrix_s_sqrt_inv=matrix_s_sqrt_inv, &
                                     matrix_cg=matrix_cg, & ! cg
                                     matrix_Ax=matrix_Ax, &
                                     eps_filter=eps_filter)

               ! Matrix projector T
               CALL projector(qs_env, matrix_p, matrix_Ax, eps_filter)

               CALL dbcsr_dot(matrix_cg(ispin)%matrix, matrix_Ax(ispin)%matrix, norm_cA(ispin))

               CPABORT("tr(Ap_j*p_j) < 0")
               IF (abnormal_value(norm_cA(ispin))) &
                  CPABORT("Preconditioner: Tr[Ap_j*p_j] is an abnormal value (NaN/Inf)")

            END IF

         END DO

         DO ispin = 1, nspins
            ! Determine step-size
            IF (norm_cA(ispin) .LT. linres_control%eps) THEN
               alpha(ispin) = 1.0_dp
            ELSE
               alpha(ispin) = norm_rr(ispin)/norm_cA(ispin)
            END IF

            ! x_j+1 = x_j + alpha*p_j
            ! save response-denisty of this iteration
            CALL dbcsr_add(matrix_cg_z(ispin)%matrix, matrix_cg(ispin)%matrix, 1.0_dp, alpha(ispin))
         END DO

         ! need to recompute the residue
         restart = .FALSE.
         IF (MOD(i, linres_control%restart_every) .EQ. 0) THEN
            !
            ! r_j+1 = b - A * x_j+1
            CALL build_hessian_op(qs_env=qs_env, &
                                  matrix_ks=matrix_ks, &
                                  matrix_p=matrix_p, &
                                  matrix_s_sqrt_inv=matrix_s_sqrt_inv, &
                                  matrix_cg=matrix_cg_z, & ! cg
                                  matrix_Ax=matrix_Ax, &
                                  eps_filter=eps_filter)
            ! b
            CALL commutator(matrix_nsc, matrix_p, matrix_res, eps_filter, .FALSE., alpha=focc)

            DO ispin = 1, nspins
               CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -1.0_dp)
            END DO

            CALL projector(qs_env, matrix_p, matrix_res, eps_filter)
            !
            restart = .TRUE.
         ELSE
            ! proj Ap onto the virtual subspace
            CALL projector(qs_env, matrix_p, matrix_Ax, eps_filter)
            !
            ! r_j+1 = r_j - alpha * Ap_j
            DO ispin = 1, nspins
               CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -alpha(ispin))
            END DO
            restart = .FALSE.
         END IF

         ! Preconditioner
         linres_control%flag = ""
         IF (linres_control%preconditioner_type == precond_mlp) THEN
            ! M * z_j+1 = r_j+1
            ! Conjugate gradient returns z_j+1
            CALL preconditioner(qs_env=qs_env, &
                                matrix_ks=matrix_ks, &
                                matrix_p=matrix_p, &
                                matrix_rhs=matrix_res, &
                                matrix_cg_z=matrix_z0, &
                                eps_filter=eps_filter, &
                                iounit=iounit)
            linres_control%flag = "PCG-AO"
         ELSE
            DO ispin = 1, nspins
               CALL dbcsr_copy(matrix_z0(ispin)%matrix, matrix_res(ispin)%matrix)
            END DO
            linres_control%flag = "CG-AO"
         END IF

         norm_res = 0.0_dp

         DO ispin = 1, nspins
            ! Tr[r_j+1*z_j+1]
            CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_z0(ispin)%matrix, new_norm(ispin))
            IF (new_norm(ispin) .LT. 0.0_dp) CPABORT("tr(r_j+1*z_j+1) < 0")
            IF (abnormal_value(new_norm(ispin))) &
               CPABORT("Preconditioner: Tr[r_j+1*z_j+1] is an abnormal value (NaN/Inf)")
            norm_res = MAX(norm_res, new_norm(ispin)/REAL(nao, dp))

            IF (norm_rr(ispin) .LT. linres_control%eps .OR. new_norm(ispin) .LT. linres_control%eps) THEN
               beta(ispin) = 0.0_dp
               linres_control%converged = .TRUE.
            ELSE
               beta(ispin) = new_norm(ispin)/norm_rr(ispin)
            END IF

            ! update new search vector (matrix cg)
            ! Here: cg_j+1 = z_j+1 + beta*cg_j
            CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_z0(ispin)%matrix, beta(ispin), 1.0_dp)
            CALL dbcsr_filter(matrix_cg(ispin)%matrix, eps_filter)

            tr_rz00(ispin) = norm_rr(ispin)
            norm_rr(ispin) = new_norm(ispin)
         END DO

         ! Can we exit the loop?
         CALL external_control(should_stop, "LS_SOLVER", target_time=qs_env%target_time, &
                               start_time=qs_env%start_time)

      END DO iteration

      ! Matrix projector
      CALL projector(qs_env, matrix_p, matrix_Ax, eps_filter)

      ! Z = [cg_z,P]
      CALL commutator(matrix_cg_z, matrix_p, matrix_z, eps_filter, .TRUE., alpha=0.5_dp)

      DO ispin = 1, nspins
         ! Transform Z-matrix back into non-orthogonal basis
         CALL transform_m_orth(matrix_z(ispin)%matrix, matrix_s_sqrt_inv, eps_filter)

         ! Export Z-Matrix
         CALL dbcsr_copy(matrix_pz(ispin)%matrix, matrix_z(ispin)%matrix, keep_sparsity=.TRUE.)
      END DO

      ! Calculate energy-weighted response density matrix
      ! AO: Wz = 0.5*(Z*KS*P + P*KS*Z)
      CALL ec_wz_matrix(qs_env, matrix_pz, matrix_wz, eps_filter)

      ! Release matrices
      CALL dbcsr_release(matrix_tmp)

      CALL dbcsr_release(matrix_s_sqrt)
      CALL dbcsr_release(matrix_s_sqrt_inv)

      CALL dbcsr_deallocate_matrix_set(matrix_p)
      CALL dbcsr_deallocate_matrix_set(matrix_ks)
      CALL dbcsr_deallocate_matrix_set(matrix_nsc)
      CALL dbcsr_deallocate_matrix_set(matrix_z)
      CALL dbcsr_deallocate_matrix_set(matrix_Ax)
      CALL dbcsr_deallocate_matrix_set(matrix_res)
      CALL dbcsr_deallocate_matrix_set(matrix_cg)
      CALL dbcsr_deallocate_matrix_set(matrix_cg_z)
      CALL dbcsr_deallocate_matrix_set(matrix_z0)

      DEALLOCATE (alpha, beta, new_norm, norm_cA, norm_rr)
      DEALLOCATE (tr_rz00)

      CALL timestop(handle)

   END SUBROUTINE ec_response_ao

! **************************************************************************************************
!> \brief Compute matrix_wz as needed for the forces
!>        Wz = 0.5*(Z*KS*P + P*KS*Z) (closed-shell)
!> \param qs_env ...
!> \param matrix_z The response density we just calculated
!> \param matrix_wz The energy weighted response-density matrix
!> \param eps_filter ...
!> \par History
!>       2020.2 created [Fabian Belleflamme]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_wz_matrix(qs_env, matrix_z, matrix_wz, eps_filter)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_z
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: matrix_wz
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_wz_matrix', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iounit, ispin, nspins
      REAL(KIND=dp)                                      :: scaling
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_p, matrix_s
      TYPE(dbcsr_type)                                   :: matrix_tmp, matrix_tmp2
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(matrix_z))
      CPASSERT(ASSOCIATED(matrix_wz))

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      matrix_ks=matrix_ks, &
                      matrix_s=matrix_s, &
                      rho=rho)
      nspins = dft_control%nspins
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_unit_nr(logger)

      CALL qs_rho_get(rho, rho_ao=matrix_p)

      ! Init temp matrices
      CALL dbcsr_create(matrix_tmp, template=matrix_z(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp2, template=matrix_z(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      ! Scale matrix_p by factor 1/2 in closed-shell
      scaling = 1.0_dp
      IF (nspins == 1) scaling = 0.5_dp

      ! Whz = ZFP + PFZ = Z(FP) + (Z(FP))^T
      DO ispin = 1, nspins

         ! tmp = FP
         CALL dbcsr_multiply("N", "N", scaling, matrix_ks(ispin)%matrix, matrix_p(ispin)%matrix, &
                             0.0_dp, matrix_tmp, filter_eps=eps_filter, retain_sparsity=.FALSE.)

         ! tmp2 = ZFP
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_z(ispin)%matrix, matrix_tmp, &
                             0.0_dp, matrix_tmp2, filter_eps=eps_filter, retain_sparsity=.FALSE.)

         ! tmp = (ZFP)^T
         CALL dbcsr_transposed(matrix_tmp, matrix_tmp2)

         ! tmp = ZFP + (ZFP)^T
         CALL dbcsr_add(matrix_tmp, matrix_tmp2, 1.0_dp, 1.0_dp)

         CALL dbcsr_filter(matrix_tmp, eps_filter)

         ! Whz = ZFP + PFZ
         CALL dbcsr_copy(matrix_wz(ispin)%matrix, matrix_tmp, keep_sparsity=.TRUE.)

      END DO

      ! Release matrices
      CALL dbcsr_release(matrix_tmp)
      CALL dbcsr_release(matrix_tmp2)

      CALL timestop(handle)

   END SUBROUTINE ec_wz_matrix

! **************************************************************************************************
!> \brief  Calculate first term of electronic Hessian  M = [F, B]
!>         acting as liner transformation on trial matrix (matrix_cg)
!>         with intermediate response density  B = [cg,P] = cg*P - P*cg = cg*P + (cg*P)^T
!>
!>         All matrices are in orthonormal basis
!>
!> \param matrix_ks Ground-state Kohn-Sham matrix
!> \param matrix_p  Ground-state Density matrix
!> \param matrix_cg Trial matrix
!> \param matrix_b  Intermediate response density
!> \param matrix_Ax First term of electronic Hessian applied on trial matrix (matrix_cg)
!>
!> \param eps_filter ...
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE hessian_op1(matrix_ks, matrix_p, matrix_cg, matrix_b, matrix_Ax, eps_filter)

      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_ks, matrix_p, matrix_cg
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: matrix_b, matrix_Ax
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(len=*), PARAMETER :: routineN = 'hessian_op1', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_p))
      CPASSERT(ASSOCIATED(matrix_cg))
      CPASSERT(ASSOCIATED(matrix_b))
      CPASSERT(ASSOCIATED(matrix_Ax))

      ! Build intermediate density matrix
      ! B = [cg, P] = cg*P - P*cg = cg*P + (cg*P)^T
      CALL commutator(matrix_cg, matrix_p, matrix_b, eps_filter, .TRUE.)

      ! Build first part of operator
      ! Ax = [F,[cg,P]] = [F,B]
      CALL commutator(matrix_ks, matrix_b, matrix_Ax, eps_filter, .FALSE.)

      CALL timestop(handle)

   END SUBROUTINE hessian_op1

! **************************************************************************************************
!> \brief  calculate lin transformation of Hessian matrix on a trial vector (matrix) matrix_cg
!>         which is stored in response density B = [cg,P] = cg*P - P*cg = cg*P + (cg*P)^T
!>         Ax = [F, B] + [G(B), Pin] in orthonormal basis
!>
!> \param qs_env ...
!> \param matrix_ks Ground-state Kohn-Sham matrix
!> \param matrix_p  Ground-state Density matrix
!> \param matrix_s_sqrt_inv S^(-1/2) needed for transformation to/from orthonormal basis
!> \param matrix_cg Trial matrix
!> \param matrix_Ax Electronic Hessian applied on trial matrix (matrix_cg)
!> \param eps_filter ...
!>
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE build_hessian_op(qs_env, matrix_ks, matrix_p, matrix_s_sqrt_inv, &
                               matrix_cg, matrix_Ax, eps_filter)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_ks, matrix_p
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s_sqrt_inv
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_cg
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: matrix_Ax
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(len=*), PARAMETER :: routineN = 'build_hessian_op', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: chksum
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_b, rho1_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_rho_type), POINTER                         :: rho, rho1

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_p))
      CPASSERT(ASSOCIATED(matrix_cg))
      CPASSERT(ASSOCIATED(matrix_Ax))

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      rho=rho)
      nspins = dft_control%nspins

      NULLIFY (matrix_b)
      CALL dbcsr_allocate_matrix_set(matrix_b, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_b(ispin)%matrix)
         CALL dbcsr_create(matrix_b(ispin)%matrix, name="[X,P] RSP DNSTY", &
                           template=matrix_p(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
      END DO

      ! Build uncoupled term of Hessian linear transformation
      CALL hessian_op1(matrix_ks, matrix_p, matrix_cg, matrix_b, matrix_Ax, eps_filter)

      ! Avoid the buildup of noisy blocks
      DO ispin = 1, nspins
         CALL dbcsr_filter(matrix_b(ispin)%matrix, eps_filter)
      END DO

      chksum = 0.0_dp
      DO ispin = 1, nspins
         chksum = chksum + dbcsr_checksum(matrix_b(ispin)%matrix)
      END DO

      ! skip the kernel if the DM is very small
      IF (chksum .GT. 1.0E-14_dp) THEN

         CALL get_qs_env(qs_env=qs_env, rho=rho)
         ! Bring B as density on grid
         NULLIFY (rho1)
         CALL qs_rho_create(rho1)
         CALL qs_rho_rebuild(rho1, qs_env, rebuild_ao=.TRUE., rebuild_grids=.TRUE.)
         ! Get response density matrix
         CALL qs_rho_get(rho1, rho_ao=rho1_ao)

         DO ispin = 1, nspins
            ! Transform B into NON-ortho basis for collocation
            CALL transform_m_orth(matrix_b(ispin)%matrix, matrix_s_sqrt_inv, eps_filter)
            ! Filter
            CALL dbcsr_filter(matrix_b(ispin)%matrix, eps_filter)
            ! Keep symmetry of density matrix
            CALL dbcsr_copy(rho1_ao(ispin)%matrix, matrix_b(ispin)%matrix, keep_sparsity=.TRUE.)
         END DO

         ! Updates densities on grid wrt density matrix
         CALL qs_rho_update_rho(rho1, qs_env)

         ! Calculate kernel
         ! Ax = F*B - B*F + G(B)*P - P*G(B)
         !                        IN/OUT     IN        IN             IN
         CALL hessian_op2(qs_env, matrix_Ax, matrix_p, matrix_s_sqrt_inv, rho1, eps_filter)

         ! Add contribution from ADMM to kernel
!       IF (dft_control%do_admm) THEN
!          CALL apply_xc_admm(qs_env, p_env, matrix_Ax, matrix_p, p_env%rho1, eps_filter)
!       END IF

         CALL qs_rho_release(rho1)

      END IF

      CALL dbcsr_deallocate_matrix_set(matrix_b)

      CALL timestop(handle)

   END SUBROUTINE build_hessian_op

! **************************************************************************************************
!> \brief  Calculate lin transformation of Hessian matrix on a trial vector (matrix) matrix_cg
!>         which is stored in response density B = [cg,P] = cg*P - P*cg = cg*P + (cg*P)^T
!>         Ax = [F, B] + [G(B), Pin] in orthonormal basis
!>
!> \param qs_env ...
!> \param matrix_Ax contains first part of Hessian linear transformation, kernel contribution
!>                  is calculated and added in this routine
!> \param matrix_p Density matrix in orthogonal basis
!> \param matrix_s_sqrt_inv contains matrix S^(-1/2) for switching to orthonormal Lowdin basis
!> \param rho1 trial density environment
!> \param eps_filter ...
!>
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE hessian_op2(qs_env, matrix_Ax, matrix_p, matrix_s_sqrt_inv, rho1, eps_filter)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: matrix_Ax
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_p
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s_sqrt_inv
      TYPE(qs_rho_type), INTENT(IN), POINTER             :: rho1
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(len=*), PARAMETER :: routineN = 'hessian_op2', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: ekin_mol
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_G, matrix_s, rho1_ao, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_g, rho1_r, v_xc
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, input, matrix_s, para_env, rho, rho1_g, rho1_r)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      input=input, &
                      matrix_s=matrix_s, &
                      para_env=para_env, &
                      rho=rho)
      nspins = dft_control%nspins

      ! Get non-ortho input density matrix on grid
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      ! Get non-ortho trial density
      CALL qs_rho_get(rho1, rho_g=rho1_g, rho_r=rho1_r)

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))

      NULLIFY (auxbas_pw_pool, poisson_env, pw_pools)
      ! gets the tmp grids
      CALL pw_env_get(pw_env=pw_env, &
                      auxbas_pw_pool=auxbas_pw_pool, &
                      pw_pools=pw_pools, &
                      poisson_env=poisson_env)

      ! Calculate the NSC Hartree potential
      CALL pw_pool_create_pw(pool=auxbas_pw_pool, &
                             pw=v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(pool=auxbas_pw_pool, &
                             pw=rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(pool=auxbas_pw_pool, &
                             pw=v_hartree_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      ! XC-Kernel
      NULLIFY (v_xc, xc_section)
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")

      ! add xc-kernel
      CALL create_kernel(qs_env, &
                         vxc=v_xc, &
                         rho=rho, &
                         rho1_r=rho1_r, &
                         rho1_g=rho1_g, &
                         xc_section=xc_section)

      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
      END DO

      ! take trial density to build G^{H}[B]
      CALL pw_zero(rho_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL pw_axpy(rho1_g(ispin)%pw, rho_tot_gspace%pw)
      END DO

      ! get Hartree potential from rho_tot_gspace
      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, &
                            vhartree=v_hartree_gspace%pw)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

      ! Add v_xc + v_H
      DO ispin = 1, nspins
         CALL pw_axpy(v_hartree_rspace%pw, v_xc(ispin)%pw)
      END DO
      IF (nspins == 1) CALL pw_scale(v_xc(1)%pw, 2.0_dp)

      ! Init response kernel matrix
      ! matrix G(B)
      NULLIFY (matrix_G)
      CALL dbcsr_allocate_matrix_set(matrix_G, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_G(ispin)%matrix)
         CALL dbcsr_copy(matrix_G(ispin)%matrix, matrix_s(1)%matrix, &
                         name="MATRIX Kernel")

         ! Integrate with ground-state density matrix, in non-orthogonal basis
         CALL integrate_v_rspace(v_rspace=v_xc(ispin), &
                                 pmat=rho_ao(ispin), &
                                 hmat=matrix_G(ispin), &
                                 qs_env=qs_env, &
                                 calculate_forces=.FALSE., &
                                 basis_type="ORB")
      END DO

      ! Calculate KG correction to kernel
      IF (dft_control%qs_control%do_kg) THEN
         IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed .OR. &
             qs_env%kg_env%tnadd_method == kg_tnadd_embed_ri) THEN

            CPASSERT(dft_control%nimages == 1)
            ekin_mol = 0.0_dp
            CALL qs_rho_get(rho1, rho_ao=rho1_ao)
            CALL kg_ekin_subset(qs_env=qs_env, &
                                ks_matrix=matrix_G, &
                                ekin_mol=ekin_mol, &
                                calc_force=.FALSE., &
                                do_kernel=.TRUE., &
                                pmat_ext=rho1_ao)
         END IF
      END IF

      ! Transforming G(B) into orthonormal basis
      DO ispin = 1, nspins
         CALL transform_m_orth(matrix_G(ispin)%matrix, matrix_s_sqrt_inv, eps_filter)
         CALL dbcsr_filter(matrix_G(ispin)%matrix, eps_filter)
      END DO

      ! Hessian already contains  Ax = [F,B] (ORTHO), now adding
      ! Ax = Ax + G(B)P - (G(B)P)^T
      CALL commutator(matrix_G, matrix_p, matrix_Ax, eps_filter, .FALSE., 1.0_dp, 1.0_dp)

      ! release pw grids
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (v_xc)

      CALL dbcsr_deallocate_matrix_set(matrix_G)

      CALL timestop(handle)

   END SUBROUTINE hessian_op2

! **************************************************************************************************
!> \brief computes (anti-)commutator exploiting (anti-)symmetry:
!>        A symmetric : RES = beta*RES + k*[A,B] = k*(AB-(AB)^T)
!>        A anti-sym  : RES = beta*RES + k*{A,B} = k*(AB+(AB)^T)
!>
!> \param a          Matrix A
!> \param b          Matrix B
!> \param res        Commutator result
!> \param eps_filter filtering threshold for sparse matrices
!> \param anticomm   Calculate anticommutator
!> \param alpha      Scaling of anti-/commutator
!> \param beta       Scaling of inital content of result matrix
!>
!> \par History
!>       2020.07 Fabian Belleflamme  (based on commutator_symm)
! **************************************************************************************************
   SUBROUTINE commutator(a, b, res, eps_filter, anticomm, alpha, beta)

      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: a, b, res
      REAL(KIND=dp)                                      :: eps_filter
      LOGICAL                                            :: anticomm
      REAL(KIND=dp), OPTIONAL                            :: alpha, beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'commutator', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: facc, myalpha, mybeta
      TYPE(dbcsr_type)                                   :: work, work2

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(a))
      CPASSERT(ASSOCIATED(b))
      CPASSERT(ASSOCIATED(res))

      CALL dbcsr_create(work, template=a(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(work2, template=a(1)%matrix, matrix_type=dbcsr_type_no_symmetry)

      ! Scaling of anti-/commutator
      myalpha = 1.0_dp
      IF (PRESENT(alpha)) myalpha = alpha
      ! Scaling of result matrix
      mybeta = 0.0_dp
      IF (PRESENT(beta)) mybeta = beta
      ! Add/subtract second term when calculating anti-/commutator
      facc = -1.0_dp
      IF (anticomm) facc = 1.0_dp

      DO ispin = 1, SIZE(a)

         CALL dbcsr_multiply("N", "N", myalpha, a(ispin)%matrix, b(ispin)%matrix, &
                             0.0_dp, work, filter_eps=eps_filter)
         CALL dbcsr_transposed(work2, work)

         ! RES= beta*RES + alpha*{A,B} = beta*RES + alpha*[AB+(AB)T]
         ! RES= beta*RES + alpha*[A,B] = beta*RES + alpha*[AB-(AB)T]
         CALL dbcsr_add(work, work2, 1.0_dp, facc)

         CALL dbcsr_add(res(ispin)%matrix, work, mybeta, 1.0_dp)

      END DO

      CALL dbcsr_release(work)
      CALL dbcsr_release(work2)

      CALL timestop(handle)

   END SUBROUTINE commutator

! **************************************************************************************************
!> \brief Projector P(M) = P*M*Q^T + Q*M*P^T
!>        with P = D
!>        with Q = (1-D)
!>
!> \param qs_env ...
!> \param matrix_p  Ground-state density in orthonormal basis
!> \param matrix_io Matrix to which projector is applied.
!>
!> \param eps_filter ...
!> \date    06.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE projector(qs_env, matrix_p, matrix_io, eps_filter)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: matrix_p
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: matrix_io
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(len=*), PARAMETER :: routineN = 'projector', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type)                                   :: matrix_q, matrix_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env)
      nspins = dft_control%nspins

      CALL dbcsr_create(matrix_q, template=matrix_p(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp, template=matrix_p(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      ! Q = (1 - P)
      CALL dbcsr_copy(matrix_q, matrix_p(1)%matrix)
      CALL dbcsr_scale(matrix_q, -1.0_dp)
      CALL dbcsr_add_on_diag(matrix_q, 1.0_dp)
      CALL dbcsr_finalize(matrix_q)

      ! Proj(M) = P*M*Q + Q*M*P
      ! with P = D = CC^T
      ! and  Q = (1 - P)
      DO ispin = 1, nspins

         ! tmp1 = P*M
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin)%matrix, matrix_io(ispin)%matrix, &
                             0.0_dp, matrix_tmp, filter_eps=eps_filter)
         ! m_io = P*M*Q
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp, matrix_q, &
                             0.0_dp, matrix_io(ispin)%matrix, filter_eps=eps_filter)

         ! tmp = (P^T*M^T*Q^T)^T = -(P*M*Q)^T
         CALL dbcsr_transposed(matrix_tmp, matrix_io(ispin)%matrix)
         CALL dbcsr_add(matrix_io(ispin)%matrix, matrix_tmp, 1.0_dp, -1.0_dp)

      END DO

      CALL dbcsr_release(matrix_tmp)
      CALL dbcsr_release(matrix_q)

      CALL timestop(handle)

   END SUBROUTINE projector

! **************************************************************************************************
!> \brief performs a tranformation of a matrix back to/into orthonormal basis
!>        in case of P a scaling of 0.5 has to be applied for closed shell case
!> \param matrix       matrix to be transformed
!> \param matrix_trafo transformation matrix
!> \param eps_filter   filtering threshold for sparse matrices
!> \par History
!>       2012.05 created [Florian Schiffmann]
!> \author Florian Schiffmann
!>
! **************************************************************************************************

   SUBROUTINE transform_m_orth(matrix, matrix_trafo, eps_filter)
      TYPE(dbcsr_type)                                   :: matrix, matrix_trafo
      REAL(KIND=dp)                                      :: eps_filter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'transform_m_orth', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: matrix_tmp, matrix_work

      CALL timeset(routineN, handle)

      CALL dbcsr_create(matrix_work, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp, template=matrix, matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix, matrix_trafo, &
                          0.0_dp, matrix_work, filter_eps=eps_filter)
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_trafo, matrix_work, &
                          0.0_dp, matrix_tmp, filter_eps=eps_filter)
      ! symmetrize results (this is again needed to make sure everything is stable)
      CALL dbcsr_transposed(matrix_work, matrix_tmp)
      CALL dbcsr_add(matrix_tmp, matrix_work, 0.5_dp, 0.5_dp)
      CALL dbcsr_copy(matrix, matrix_tmp)

      ! Avoid the buildup of noisy blocks
      CALL dbcsr_filter(matrix, eps_filter)

      CALL dbcsr_release(matrix_tmp)
      CALL dbcsr_release(matrix_work)
      CALL timestop(handle)

   END SUBROUTINE transform_m_orth

END MODULE ec_orth_solver
