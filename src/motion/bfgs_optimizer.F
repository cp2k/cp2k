!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for Geometry optimization using BFGS algorithm
!> \par History
!>      Module modified by Pierre-André Cazade [pcazade] 01.2020 - University of Limerick.
!>      Modifications enable Space Group Symmetry.
! **************************************************************************************************
MODULE bfgs_optimizer

   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE atomic_kind_types,               ONLY: get_atomic_kind,&
                                              get_atomic_kind_set
   USE bibliography,                    ONLY: Lindh1995,&
                                              cite_reference
   USE cell_opt_utils,                  ONLY: apply_cell_constraints
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE constraint_fxd,                  ONLY: fix_atom_control
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_external_control,             ONLY: external_control
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_matvec,&
                                              cp_fm_transpose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_read_unformatted,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_write_unformatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_iterate,&
                                              cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE global_types,                    ONLY: global_environment_type
   USE gopt_f_methods,                  ONLY: cp_eval_at,&
                                              gopt_f_ii,&
                                              gopt_f_io,&
                                              gopt_f_io_finalize,&
                                              gopt_f_io_init,&
                                              print_geo_opt_header,&
                                              print_geo_opt_nc
   USE gopt_f_types,                    ONLY: gopt_f_type
   USE gopt_param_types,                ONLY: gopt_param_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: bfgs_direct,&
                                              bfgs_inverse,&
                                              default_cell_method_id,&
                                              default_ts_method_id
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE iso_fortran_env,                 ONLY: iostat_end
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_list_types,             ONLY: particle_list_type
   USE space_groups,                    ONLY: identify_space_group,&
                                              print_spgr,&
                                              spgr_apply_rotations_coord,&
                                              spgr_apply_rotations_force
   USE space_groups_types,              ONLY: spgr_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bfgs_optimizer'
   LOGICAL, PARAMETER                   :: debug_this_module = .TRUE.

   ! Both matrix representations share one restart format and the existing file name.
   CHARACTER(LEN=16), PARAMETER         :: bfgs_restart_magic = "CP2K BFGS STATE "
   INTEGER, PARAMETER                   :: bfgs_restart_version = 1
   REAL(KIND=dp), PARAMETER             :: min_hessian_eigenvalue = 1.0e-4_dp

   PUBLIC :: bfgs_read_restart_header, geoopt_bfgs

CONTAINS

! **************************************************************************************************
!> \brief Main driver for BFGS geometry optimizations
!> \param force_env ...
!> \param gopt_param ...
!> \param globenv ...
!> \param geo_section ...
!> \param gopt_env ...
!> \param x0 ...
!> \par History
!>      01.2020 modified to perform Space Group Symmetry [pcazade]
! **************************************************************************************************
   RECURSIVE SUBROUTINE geoopt_bfgs(force_env, gopt_param, globenv, geo_section, gopt_env, x0)

      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x0

      CHARACTER(len=*), PARAMETER                        :: routineN = 'geoopt_bfgs'
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      CHARACTER(LEN=5)                                   :: wildcard
      INTEGER                                            :: handle, info, iter_nr, its, &
                                                            matrix_method, maxiter, ndf, nfree, &
                                                            output_unit
      LOGICAL                                            :: conv, hesrest, hess_reset, inverse, &
                                                            shell_present, should_stop, &
                                                            use_mod_hes, use_rfo
      REAL(KIND=dp)                                      :: ediff, emin, eold, etot, pred, rad, rat, &
                                                            step, t_diff, t_now, t_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dg, dr, dx, eigval, gold, work, xold
      REAL(KIND=dp), DIMENSION(:), POINTER               :: g
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_hes
      TYPE(cp_fm_type)                                   :: eigvec_mat, hess_mat, hess_tmp
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: print_key, root_section
      TYPE(spgr_type), POINTER                           :: spgr

      NULLIFY (logger, g, blacs_env, spgr)
      logger => cp_get_default_logger()
      para_env => force_env%para_env
      root_section => force_env%root_section
      spgr => gopt_env%spgr
      t_old = m_walltime()

      CALL timeset(routineN, handle)
      CALL section_vals_val_get(geo_section, "BFGS%TRUST_RADIUS", r_val=rad)
      print_key => section_vals_get_subs_vals(geo_section, "BFGS%RESTART")
      CALL section_vals_val_get(geo_section, "BFGS%METHOD", i_val=matrix_method)
      inverse = matrix_method == bfgs_inverse
      maxiter = gopt_param%max_iter
      conv = .FALSE.
      rat = 0.0_dp
      wildcard = " BFGS"

      ! Stop if not yet implemented
      SELECT CASE (gopt_env%type_id)
      CASE (default_ts_method_id)
         CPABORT("BFGS method not yet working with DIMER")
      END SELECT

      CALL section_vals_val_get(geo_section, "BFGS%USE_RAT_FUN_OPT", l_val=use_rfo)
      CALL section_vals_val_get(geo_section, "BFGS%USE_MODEL_HESSIAN", l_val=use_mod_hes)
      CALL section_vals_val_get(geo_section, "BFGS%RESTART_HESSIAN", l_val=hesrest)
      IF (inverse .AND. use_rfo) THEN
         CALL cp_abort(__LOCATION__, &
                       "BFGS METHOD INVERSE does not support USE_RAT_FUN_OPT")
      END IF
      IF (inverse) THEN
         IF (.NOT. ieee_is_finite(rad) .OR. rad <= 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "BFGS TRUST_RADIUS must be finite and positive")
         END IF
      END IF
      output_unit = cp_print_key_unit_nr(logger, geo_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".geoLog")
      IF (output_unit > 0) THEN
         IF (inverse) THEN
            WRITE (UNIT=output_unit, FMT="(/,T2,A)") "BFGS| Matrix update method: INVERSE"
         ELSE
            WRITE (UNIT=output_unit, FMT="(/,T2,A)") "BFGS| Matrix update method: DIRECT"
         END IF
         IF (use_rfo) THEN
            WRITE (UNIT=output_unit, FMT="(/,T2,A,T78,A3)") &
               "BFGS| Use rational function optimization for step estimation: ", "YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(/,T2,A,T78,A3)") &
               "BFGS| Use rational function optimization for step estimation: ", " NO"
         END IF
         IF (use_mod_hes) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "BFGS| Use model Hessian for initial guess: ", "YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "BFGS| Use model Hessian for initial guess: ", " NO"
         END IF
         IF (hesrest) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "BFGS| Restart Hessian: ", "YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "BFGS| Restart Hessian: ", " NO"
         END IF
         WRITE (UNIT=output_unit, FMT="(T2,A,T61,F20.3)") &
            "BFGS| Trust radius: ", rad
      END IF

      ndf = SIZE(x0)
      nfree = gopt_env%nfree
      IF (ndf > 3000 .AND. .NOT. inverse) THEN
         CALL cp_warn(__LOCATION__, &
                      "The dimension of the Hessian matrix ("// &
                      TRIM(ADJUSTL(cp_to_string(ndf)))//") is greater than 3000. "// &
                      "The diagonalisation of the full Hessian  matrix needed for BFGS "// &
                      "is computationally expensive. You should consider to use the linear "// &
                      "scaling variant L-BFGS instead.")
      END IF

      ! Initialize hessian (hes = unitary matrix or model hessian )
      CALL cp_blacs_env_create(blacs_env, para_env, globenv%blacs_grid_layout, &
                               globenv%blacs_repeatable)
      CALL cp_fm_struct_create(fm_struct_hes, para_env=para_env, context=blacs_env, &
                               nrow_global=ndf, ncol_global=ndf)
      CALL cp_fm_create(hess_mat, fm_struct_hes, name="hess_mat")

      CALL force_env_get(force_env=force_env, subsys=subsys)
      CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els, shell_present=shell_present)
      IF (use_mod_hes) THEN
         IF (shell_present) THEN
            CALL cp_warn(__LOCATION__, &
                         "No model Hessian is available for core-shell models. "// &
                         "A unit matrix is used as the initial Hessian.")
            use_mod_hes = .FALSE.
         END IF
         IF (gopt_env%type_id == default_cell_method_id) THEN
            CALL cp_warn(__LOCATION__, &
                         "No model Hessian is available for cell optimizations. "// &
                         "A unit matrix is used as the initial Hessian.")
            use_mod_hes = .FALSE.
         END IF
      END IF

      ! hess_mat stores B for DIRECT and H = B^{-1} for INVERSE.
      ! DIRECT keeps the spectral work arrays for step construction; INVERSE
      ! only needs them while building a model-Hessian initial guess.
      IF (.NOT. inverse .OR. (use_mod_hes .AND. .NOT. hesrest)) THEN
         CALL cp_fm_create(hess_tmp, fm_struct_hes, name="hess_tmp")
         CALL cp_fm_create(eigvec_mat, fm_struct_hes, name="eigvec_mat")
         ALLOCATE (eigval(ndf))
         eigval(:) = zero
      END IF

      IF (inverse) THEN
         IF (use_mod_hes .AND. .NOT. hesrest) THEN
            CALL build_initial_hess(hess_mat, hess_tmp, eigvec_mat, eigval, force_env, &
                                    use_mod_hes, .TRUE.)
         ELSE
            CALL cp_fm_set_all(hess_mat, alpha=zero, beta=one)
         END IF
      ELSE
         CALL build_initial_hess(hess_mat, hess_tmp, eigvec_mat, eigval, force_env, &
                                 use_mod_hes, .FALSE.)
      END IF

      IF (inverse .AND. ALLOCATED(eigval)) THEN
         CALL cp_fm_release(eigvec_mat)
         CALL cp_fm_release(hess_tmp)
         DEALLOCATE (eigval)
      END IF

      ! Restore the reference coordinates before symmetry setup and the first evaluation.
      IF (hesrest) CALL read_bfgs_hessian(geo_section, hess_mat, logger, gopt_env, x0, matrix_method)

      ALLOCATE (xold(ndf))
      xold(:) = x0(:)

      ALLOCATE (g(ndf))
      g(:) = 0.0_dp

      ALLOCATE (gold(ndf))
      gold(:) = 0.0_dp

      ALLOCATE (dx(ndf))
      dx(:) = 0.0_dp

      ALLOCATE (dg(ndf))
      dg(:) = 0.0_dp

      ALLOCATE (work(ndf))
      work(:) = 0.0_dp

      ALLOCATE (dr(ndf))
      dr(:) = 0.0_dp

      ! find space_group
      IF (spgr%keep_space_group) THEN
         CALL identify_space_group(subsys, geo_section, gopt_env, output_unit)
         CALL spgr_apply_rotations_coord(spgr, x0)
         CALL print_spgr(spgr)
      END IF

      ! Geometry optimization starts now
      CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=iter_nr)
      CALL print_geo_opt_header(gopt_env, output_unit, wildcard)

      ! Calculate Energy & Gradients
      CALL cp_eval_at(gopt_env, x0, etot, g, master=gopt_env%force_env%para_env%mepos, &
                      para_env=gopt_env%force_env%para_env)

      ! Symmetrize coordinates and forces
      IF (spgr%keep_space_group) THEN
         CALL spgr_apply_rotations_coord(spgr, x0)
         CALL spgr_apply_rotations_force(spgr, g)
      END IF

      ! Print info at time 0
      emin = etot
      t_now = m_walltime()
      t_diff = t_now - t_old
      t_old = t_now
      CALL gopt_f_io_init(gopt_env, output_unit, etot, wildcard=wildcard, its=iter_nr, used_time=t_diff)
      DO its = iter_nr + 1, maxiter
         CALL cp_iterate(logger%iter_info, last=(its == maxiter))
         CALL section_vals_val_set(geo_section, "STEP_START_VAL", i_val=its)
         CALL gopt_f_ii(its, output_unit)

         IF ((its - iter_nr) > 1) THEN
            ! Symmetrize old coordinates and old forces.
            IF (spgr%keep_space_group) THEN
               CALL spgr_apply_rotations_coord(spgr, xold)
               CALL spgr_apply_rotations_force(spgr, gold)
            END IF

            ! cp_eval_at may project the proposed geometry: use the actual displacement.
            dx(:) = x0(:) - xold(:)
            dg(:) = g(:) - gold(:)
            IF (inverse) THEN
               CALL bfgs_inverse_update(ndf, dx, dg, hess_mat, work)
            ELSE
               CALL bfgs(ndf, dx, dg, hess_mat, work, para_env, output_unit, hess_reset)
               IF (hess_reset) THEN
                  CALL build_initial_hess(hess_mat, hess_tmp, eigvec_mat, eigval, force_env, &
                                          use_mod_hes, .FALSE.)
               END IF
            END IF

            IF (spgr%keep_space_group) THEN
               CALL spgr_apply_rotations_force(spgr, dx)
               CALL spgr_apply_rotations_force(spgr, dg)
            END IF

            IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
               CALL write_bfgs_hessian(geo_section, hess_mat, logger, gopt_env, matrix_method)
            END IF
         END IF

         ! Symmetrize coordinates and forces
         IF (spgr%keep_space_group) THEN
            CALL spgr_apply_rotations_coord(spgr, x0)
            CALL spgr_apply_rotations_force(spgr, g)
         END IF

         ! Setting the present positions & gradients as old
         xold(:) = x0
         gold(:) = g

         IF (inverse) THEN
            CALL bfgs_inverse_step(hess_mat, g, dr)
         ELSE
            ! Copying hessian hes to (ndf x ndf) matrix hes_mat for diagonalization
            CALL cp_fm_to_fm(hess_mat, hess_tmp)

            CALL choose_eigv_solver(hess_tmp, eigvec_mat, eigval, info=info)

            ! In rare cases the diagonalization of hess_mat fails (bug in scalapack?)
            IF (info /= 0) THEN
               IF (output_unit > 0) THEN
                  WRITE (output_unit, *) &
                     "BFGS: Matrix diagonalization failed, resetting Hessian to unity."
               END IF
               CALL cp_fm_set_all(hess_mat, alpha=zero, beta=one)
               CALL cp_fm_to_fm(hess_mat, hess_tmp)
               CALL choose_eigv_solver(hess_tmp, eigvec_mat, eigval)
            END IF

            IF (use_rfo) THEN
               CALL set_hes_eig(ndf, eigval, work)
               dx(:) = eigval
               CALL rat_fun_opt(ndf, dg, eigval, work, eigvec_mat, g, para_env)
            END IF
            CALL geoopt_get_step(ndf, eigval, eigvec_mat, hess_tmp, dr, g, para_env, use_rfo)

         END IF
         IF (gopt_env%type_id == default_cell_method_id) THEN
            CALL bfgs_project_cell_step(gopt_env, dr)
         END IF

         ! Symmetrize dr
         IF (spgr%keep_space_group) THEN
            CALL spgr_apply_rotations_force(spgr, dr)
         END IF

         IF (inverse) THEN
            ! Avoid a zero-step division and initialize the reported scale on every step.
            rat = 1.0_dp
            step = MAXVAL(ABS(dr))
            IF (step > rad) THEN
               rat = rad/step
               dr(:) = rat*dr(:)
               step = rad
            END IF
         ELSE
            CALL trust_radius(ndf, step, rad, rat, dr, output_unit)
         END IF

         ! Update the atomic positions
         x0 = x0 + dr

         ! Symmetrize coordinates
         IF (spgr%keep_space_group) THEN
            CALL spgr_apply_rotations_coord(spgr, x0)
         END IF

         eold = etot

         ! Energy & Gradients at new step
         CALL cp_eval_at(gopt_env, x0, etot, g, master=gopt_env%force_env%para_env%mepos, &
                         para_env=gopt_env%force_env%para_env)

         ! Use the displacement actually accepted by constraints and cell projection.
         dr(:) = x0(:) - xold(:)
         step = MAXVAL(ABS(dr))

         ! This diagnostic requires B, not H. It is optional in gopt_f_io.
         IF (.NOT. inverse) CALL energy_predict(ndf, work, hess_mat, dr, gold, conv, pred, para_env)

         ediff = etot - eold

         ! Symmetrize forces
         IF (spgr%keep_space_group) THEN
            CALL spgr_apply_rotations_force(spgr, g)
         END IF

         ! check for an external exit command
         CALL external_control(should_stop, "GEO", globenv=globenv)
         IF (should_stop) EXIT

         ! Some IO and Convergence check
         t_now = m_walltime()
         t_diff = t_now - t_old
         t_old = t_now
         IF (inverse) THEN
            CALL gopt_f_io(gopt_env, force_env, root_section, its, etot, output_unit, &
                           eold, emin, wildcard, gopt_param, ndf, dr, g, conv, &
                           rat=rat, step=step, rad=rad, used_time=t_diff)
         ELSE
            CALL gopt_f_io(gopt_env, force_env, root_section, its, etot, output_unit, &
                           eold, emin, wildcard, gopt_param, ndf, dr, g, conv, pred, rat, &
                           step, rad, used_time=t_diff)
         END IF

         IF (conv .OR. (its == maxiter)) EXIT
         IF (etot < emin) emin = etot
         IF (use_rfo) CALL update_trust_rad(rat, rad, step, ediff)
      END DO

      IF (its == maxiter .AND. (.NOT. conv)) THEN
         CALL print_geo_opt_nc(gopt_env, output_unit)
      END IF

      ! show space_group
      IF (spgr%show_space_group) THEN
         CALL identify_space_group(subsys, geo_section, gopt_env, output_unit)
         CALL print_spgr(spgr)
      END IF

      ! Write final  information, if converged
      CALL cp_iterate(logger%iter_info, last=.TRUE., increment=0)
      CALL write_bfgs_hessian(geo_section, hess_mat, logger, gopt_env, matrix_method)
      CALL gopt_f_io_finalize(gopt_env, force_env, x0, conv, its, root_section, &
                              gopt_env%force_env%para_env, gopt_env%force_env%para_env%mepos, output_unit)

      CALL cp_fm_struct_release(fm_struct_hes)
      CALL cp_fm_release(hess_mat)
      IF (.NOT. inverse) THEN
         CALL cp_fm_release(eigvec_mat)
         CALL cp_fm_release(hess_tmp)
      END IF

      CALL cp_blacs_env_release(blacs_env)
      DEALLOCATE (xold)
      DEALLOCATE (g)
      DEALLOCATE (gold)
      DEALLOCATE (dx)
      DEALLOCATE (dg)
      IF (ALLOCATED(eigval)) DEALLOCATE (eigval)
      DEALLOCATE (work)
      DEALLOCATE (dr)

      CALL cp_print_key_finished_output(output_unit, logger, geo_section, &
                                        "PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE geoopt_bfgs

! **************************************************************************************************
!> \brief Damped inverse-BFGS update, using only a matvec and local rank-two updates.
!> \param ndf Number of optimization variables
!> \param dx Actual change of optimization variables (not modified by damping)
!> \param dg Change of gradients in the same coordinate system
!> \param hess_mat Inverse Hessian, updated in place
!> \param work Replicated workspace for H*dg
! **************************************************************************************************
   SUBROUTINE bfgs_inverse_update(ndf, dx, dg, hess_mat, work)
      INTEGER, INTENT(IN)                                :: ndf
      REAL(KIND=dp), INTENT(IN)                          :: dx(ndf), dg(ndf)
      TYPE(cp_fm_type), INTENT(IN)                       :: hess_mat
      REAL(KIND=dp), INTENT(OUT)                         :: work(ndf)

      REAL(KIND=dp), PARAMETER                           :: min_curvature = 0.2_dp

      INTEGER                                            :: handle, i, ig, invalid, j, jg, nc, nr
      INTEGER, DIMENSION(:), POINTER                     :: cols, rows
      REAL(KIND=dp)                                      :: coeff, rho, sy, theta, yhy
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: data
      REAL(KIND=dp), DIMENSION(ndf)                      :: secant

      CALL timeset("bfgs_inverse_update", handle)
      IF (.NOT. ALL(ieee_is_finite(dx)) .OR. .NOT. ALL(ieee_is_finite(dg))) THEN
         CPABORT("Non-finite displacement or gradient difference in inverse BFGS")
      END IF
      work = 0.0_dp
      IF (NORM2(dx) == 0.0_dp .OR. NORM2(dg) == 0.0_dp) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL cp_fm_matvec(hess_mat, dg, work)
      sy = DOT_PRODUCT(dx, dg)
      yhy = DOT_PRODUCT(dg, work)
      IF (.NOT. ieee_is_finite(sy) .OR. .NOT. ieee_is_finite(yhy)) THEN
         CPABORT("Non-finite curvature in inverse BFGS")
      END IF
      IF (yhy <= 0.0_dp) THEN
         CALL cp_warn(__LOCATION__, &
                      "BFGS: non-positive inverse curvature; resetting to unity.")
         CALL cp_fm_set_all(hess_mat, alpha=0.0_dp, beta=1.0_dp)
         CALL timestop(handle)
         RETURN
      END IF

      secant = dx
      IF (sy < min_curvature*yhy) THEN
         ! Dual-space damping: y^T*s_tilde >= min_curvature*y^T*H*y.
         ! Unlike damping y with B*s, this does not need the direct Hessian.
         theta = ((1.0_dp - min_curvature)*yhy)/(yhy - sy)
         secant = theta*dx + (1.0_dp - theta)*work
         sy = DOT_PRODUCT(secant, dg)
      END IF
      IF (sy <= TINY(1.0_dp)/EPSILON(1.0_dp)) THEN
         CALL timestop(handle)
         RETURN
      END IF
      rho = 1.0_dp/sy
      coeff = (1.0_dp + yhy*rho)*rho
      IF (.NOT. ieee_is_finite(coeff)) THEN
         CALL cp_warn(__LOCATION__, &
                      "BFGS: unrepresentable inverse update; retaining the previous matrix.")
         CALL timestop(handle)
         RETURN
      END IF

      CALL cp_fm_get_info(hess_mat, row_indices=rows, col_indices=cols, &
                          local_data=DATA, nrow_local=nr, ncol_local=nc)
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(nc,nr,rows,cols,data,coeff,rho,secant,work) PRIVATE(i,j,ig,jg)
      DO j = 1, nc
         jg = cols(j)
         DO i = 1, nr
            ig = rows(i)
            DATA(i, j) = DATA(i, j) + coeff*(secant(ig)*secant(jg)) - &
                         rho*(secant(ig)*work(jg) + work(ig)*secant(jg))
         END DO
      END DO
      !$OMP END PARALLEL DO
      invalid = 0
      IF (.NOT. ALL(ieee_is_finite(DATA(1:nr, 1:nc)))) invalid = 1
      CALL hess_mat%matrix_struct%para_env%max(invalid)
      IF (invalid /= 0) THEN
         CALL cp_warn(__LOCATION__, &
                      "BFGS: non-finite inverse matrix; resetting to unity.")
         CALL cp_fm_set_all(hess_mat, alpha=0.0_dp, beta=1.0_dp)
      END IF
      CALL timestop(handle)
   END SUBROUTINE bfgs_inverse_update

! **************************************************************************************************
!> \brief Form -H*g and recover a descent direction if roundoff has corrupted H.
!> \param hess_mat Inverse Hessian
!> \param g Gradient
!> \param dr Search direction
! **************************************************************************************************
   SUBROUTINE bfgs_inverse_step(hess_mat, g, dr)
      TYPE(cp_fm_type), INTENT(IN)                       :: hess_mat
      REAL(KIND=dp), INTENT(IN)                          :: g(:)
      REAL(KIND=dp), INTENT(OUT)                         :: dr(:)

      REAL(KIND=dp)                                      :: slope

      IF (.NOT. ALL(ieee_is_finite(g))) THEN
         CALL cp_abort(__LOCATION__, "Non-finite gradient in inverse BFGS")
      END IF
      dr = 0.0_dp
      IF (MAXVAL(ABS(g)) == 0.0_dp) RETURN
      CALL cp_fm_matvec(hess_mat, g, dr, alpha=-1.0_dp)
      slope = DOT_PRODUCT(g, dr)
      IF (.NOT. ieee_is_finite(slope) .OR. slope >= 0.0_dp) THEN
         CALL cp_warn(__LOCATION__, &
                      "BFGS: inverse matrix is not a descent model; resetting to unity.")
         CALL cp_fm_set_all(hess_mat, alpha=0.0_dp, beta=1.0_dp)
         dr = -g
      END IF
   END SUBROUTINE bfgs_inverse_step

! **************************************************************************************************
!> \brief Apply existing cell restrictions to a proposed BFGS direction.
!> \param gopt_env Optimization environment
!> \param dr Search direction
! **************************************************************************************************
   SUBROUTINE bfgs_project_cell_step(gopt_env, dr)
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT), TARGET :: dr

      REAL(KIND=dp), DIMENSION(:), POINTER               :: cell_step

      cell_step => dr(SIZE(dr) - 5:SIZE(dr))
      CALL apply_cell_constraints(cell_step, gopt_env%cell_env%opt_cell, &
                                  gopt_env%cell_env%keep_angles, gopt_env%cell_env%keep_symmetry, &
                                  gopt_env%cell_env%constraint_id)
   END SUBROUTINE bfgs_project_cell_step

! **************************************************************************************************
!> \brief Convert a restart matrix once, not during normal inverse-BFGS iterations.
!> \param matrix Input matrix, replaced by its (regularized) inverse
!> \param from_inverse Whether the input already represents the inverse Hessian
! **************************************************************************************************
   SUBROUTINE bfgs_invert_matrix(matrix, from_inverse)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix
      LOGICAL, INTENT(IN)                                :: from_inverse

      INTEGER                                            :: handle, info, ndf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type)                                   :: eigenvectors, scratch

      CALL timeset("bfgs_invert_matrix", handle)
      CALL cp_fm_get_info(matrix, matrix_struct=matrix_struct, nrow_global=ndf)
      CALL cp_fm_create(scratch, matrix_struct, name="bfgs_restart_scratch")
      CALL cp_fm_create(eigenvectors, matrix_struct, name="bfgs_restart_eigenvectors")
      ALLOCATE (eigenvalues(ndf))
      CALL cp_fm_to_fm(matrix, scratch)
      CALL choose_eigv_solver(scratch, eigenvectors, eigenvalues, info=info)
      IF (info /= 0) THEN
         CALL cp_abort(__LOCATION__, "BFGS restart matrix diagonalization failed")
      END IF
      IF (.NOT. ALL(ieee_is_finite(eigenvalues))) THEN
         CALL cp_abort(__LOCATION__, "Non-finite BFGS restart spectrum")
      END IF
      IF (from_inverse) THEN
         ! Do not apply the direct-Hessian cutoff to inverse-Hessian eigenvalues.
         IF (MINVAL(eigenvalues) <= 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "BFGS inverse restart matrix is not positive definite")
         END IF
      ELSE
         IF (MINVAL(eigenvalues) < min_hessian_eigenvalue) THEN
            CALL cp_warn(__LOCATION__, &
                         "BFGS: regularizing the direct restart before conversion to INVERSE.")
         END IF
         eigenvalues(:) = MAX(min_hessian_eigenvalue, eigenvalues(:))
      END IF
      eigenvalues(:) = 1.0_dp/eigenvalues(:)
      IF (.NOT. ALL(ieee_is_finite(eigenvalues))) THEN
         CALL cp_abort(__LOCATION__, &
                       "Unrepresentable inverse of BFGS restart matrix")
      END IF
      CALL cp_fm_to_fm(eigenvectors, scratch)
      CALL cp_fm_column_scale(scratch, eigenvalues)
      CALL parallel_gemm("N", "T", ndf, ndf, ndf, 1.0_dp, scratch, eigenvectors, 0.0_dp, matrix)
      CALL cp_fm_release(scratch)
      CALL cp_fm_release(eigenvectors)
      DEALLOCATE (eigenvalues)
      CALL timestop(handle)
   END SUBROUTINE bfgs_invert_matrix

! **************************************************************************************************
!> \brief Read the metadata of a tagged BFGS restart, or identify a legacy raw Hessian.
!> \param unit File unit, returned positioned at the first matrix column
!> \param stored_method Matrix representation in the file
!> \param stored_ndf Matrix dimension stored in a tagged file; zero for legacy files
!> \param stored_type Optimization type stored in a tagged file; zero for legacy files
!> \param tagged Whether the file contains coordinate and representation metadata
!> \param reference_cell Stored CELL_OPT reference cell
!> \param cell_frame Stored canonical-to-input cell frame
! **************************************************************************************************
   SUBROUTINE bfgs_read_restart_header(unit, stored_method, stored_ndf, stored_type, tagged, &
                                       reference_cell, cell_frame)
      INTEGER, INTENT(IN)                                :: unit
      INTEGER, INTENT(OUT)                               :: stored_method, stored_ndf, stored_type
      LOGICAL, INTENT(OUT)                               :: tagged
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: reference_cell, cell_frame

      CHARACTER(LEN=16)                                  :: magic
      INTEGER                                            :: ios, version

      stored_method = bfgs_direct
      stored_ndf = 0
      stored_type = 0
      reference_cell = 0.0_dp
      cell_frame = 0.0_dp
      magic = ""
      READ (unit, IOSTAT=ios) magic
      tagged = ios == 0 .AND. magic == bfgs_restart_magic
      REWIND (unit)
      IF (tagged) THEN
         READ (unit, IOSTAT=ios) magic, version, stored_method, stored_ndf, stored_type, &
            reference_cell, cell_frame
         IF (ios /= 0) CPABORT("Invalid BFGS restart header")
         IF (version /= bfgs_restart_version) THEN
            CALL cp_abort(__LOCATION__, "Unsupported BFGS restart version")
         END IF
         IF (stored_method /= bfgs_direct .AND. stored_method /= bfgs_inverse) THEN
            CPABORT("Unknown matrix representation in BFGS restart")
         END IF
         IF (stored_ndf <= 0) CPABORT("Invalid matrix dimension in BFGS restart")
         IF (.NOT. ALL(ieee_is_finite(reference_cell)) .OR. .NOT. ALL(ieee_is_finite(cell_frame))) THEN
            CPABORT("Non-finite coordinate frame in BFGS restart")
         END IF
      END IF
   END SUBROUTINE bfgs_read_restart_header

! **************************************************************************************************
!> \brief Read and validate a tagged or legacy BFGS restart on the I/O rank.
!> \param unit File unit, returned positioned at the first matrix column
!> \param ndf Expected dimension
!> \param opt_type Expected optimization type
!> \param stored_method Matrix representation in the file
!> \param tagged Whether the file contains coordinate and representation metadata
!> \param reference_cell Stored CELL_OPT reference cell
!> \param cell_frame Stored canonical-to-input cell frame
! **************************************************************************************************
   SUBROUTINE bfgs_read_header(unit, ndf, opt_type, stored_method, tagged, reference_cell, cell_frame)
      INTEGER, INTENT(IN)                                :: unit, ndf, opt_type
      INTEGER, INTENT(OUT)                               :: stored_method
      LOGICAL, INTENT(OUT)                               :: tagged
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: reference_cell, cell_frame

      INTEGER                                            :: ios, j, stored_ndf, stored_type
      REAL(KIND=dp)                                      :: extra
      REAL(KIND=dp), DIMENSION(ndf)                      :: column

      CALL bfgs_read_restart_header(unit, stored_method, stored_ndf, stored_type, tagged, &
                                    reference_cell, cell_frame)
      IF (tagged) THEN
         IF (stored_ndf /= ndf) THEN
            CALL cp_abort(__LOCATION__, &
                          "BFGS restart dimension does not match the optimization variables")
         END IF
         IF (stored_type /= opt_type) THEN
            CALL cp_abort(__LOCATION__, "BFGS restart optimization type does not match")
         END IF
      END IF

      ! Preflight before any collective read: reject truncated/non-finite payloads and
      ! incompatible legacy dimensions with a CP2K error, not a rank-local Fortran I/O failure.
      DO j = 1, ndf
         READ (unit, IOSTAT=ios) column
         IF (ios /= 0) THEN
            CALL cp_abort(__LOCATION__, &
                          "Truncated BFGS restart or incompatible legacy matrix dimension")
         END IF
         IF (.NOT. ALL(ieee_is_finite(column))) THEN
            CALL cp_abort(__LOCATION__, "Non-finite matrix in BFGS restart")
         END IF
      END DO
      READ (unit, IOSTAT=ios) extra
      IF (ios /= iostat_end) THEN
         CALL cp_abort(__LOCATION__, &
                       "Extra matrix columns or invalid trailing data in BFGS restart")
      END IF
      REWIND (unit)
      IF (tagged) READ (unit)
   END SUBROUTINE bfgs_read_header

! **************************************************************************************************
!> \brief Read a common BFGS restart, restoring its coordinates and converting if needed.
!> \param geo_section BFGS input section
!> \param hess_mat Matrix to initialize
!> \param logger Logger
!> \param gopt_env Optimization environment
!> \param x0 Current coordinates, re-expressed in the stored reference cell for CELL_OPT
!> \param matrix_method Requested in-memory representation
! **************************************************************************************************
   SUBROUTINE read_bfgs_hessian(geo_section, hess_mat, logger, gopt_env, x0, matrix_method)
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(cp_fm_type), INTENT(INOUT)                    :: hess_mat
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: x0
      INTEGER, INTENT(IN)                                :: matrix_method

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: ndf, stored_method, unit
      LOGICAL                                            :: tagged
      REAL(KIND=dp), DIMENSION(3, 3)                     :: cell_frame, reference_cell
      TYPE(mp_para_env_type), POINTER                    :: para_env

      para_env => gopt_env%force_env%para_env
      ndf = SIZE(x0)
      unit = -1
      IF (para_env%is_source()) THEN
         CALL section_vals_val_get(geo_section, "BFGS%RESTART_FILE_NAME", c_val=filename)
         IF (LEN_TRIM(filename) == 0) filename = TRIM(logger%iter_info%project_name)//"-BFGS.Hessian"
         CALL open_file(file_name=TRIM(filename), file_status="OLD", &
                        file_form="UNFORMATTED", file_action="READ", unit_number=unit)
         CALL bfgs_read_header(unit, ndf, gopt_env%type_id, stored_method, tagged, reference_cell, cell_frame)
      END IF
      CALL para_env%bcast(stored_method, 0)
      CALL para_env%bcast(tagged, 0)
      CALL para_env%bcast(reference_cell, 0)
      CALL para_env%bcast(cell_frame, 0)
      CALL cp_fm_read_unformatted(hess_mat, unit)
      IF (para_env%is_source()) CALL close_file(unit_number=unit)

      IF (gopt_env%type_id == default_cell_method_id) THEN
         IF (tagged) THEN
            CALL bfgs_restore_cell_frame(gopt_env, x0, reference_cell, cell_frame)
         ELSE
            CALL cp_warn(__LOCATION__, "Legacy BFGS restart has no reference cell; "// &
                         "assuming the current CELL_OPT coordinate frame. Exact frame recovery is not possible.")
         END IF
      END IF
      IF (stored_method /= matrix_method) THEN
         CALL bfgs_invert_matrix(hess_mat, stored_method == bfgs_inverse)
      END IF
   END SUBROUTINE read_bfgs_hessian

! **************************************************************************************************
!> \brief Re-express coordinates in the saved frame without changing physical positions.
!> \param gopt_env Optimization environment
!> \param x0 Optimization coordinates
!> \param reference_cell Reference cell that defines the stored Hessian coordinates
!> \param cell_frame Canonical-to-input cell orientation used when writing the matrix
! **************************************************************************************************
   SUBROUTINE bfgs_restore_cell_frame(gopt_env, x0, reference_cell, cell_frame)
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: x0
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: reference_cell, cell_frame

      INTEGER                                            :: i
      REAL(KIND=dp), DIMENSION(3, 3)                     :: transform

      ! A normal continuation preserves this orientation even when lengths and angles change.
      ! Do not silently attach a matrix to a separately rotated or reordered structure.
      IF (MAXVAL(ABS(cell_frame - gopt_env%cell_env%opt_to_input)) > 1.0E-8_dp) THEN
         CALL cp_abort(__LOCATION__, &
                       "BFGS restart cell orientation differs; use a matching structure or RESTART_HESSIAN F")
      END IF
      IF (ABS(det_3x3(reference_cell)) <= EPSILON(1.0_dp)*MAXVAL(ABS(reference_cell))**3) THEN
         CPABORT("Singular reference cell in BFGS restart")
      END IF
      transform = MATMUL(reference_cell, inv_3x3(gopt_env%h_ref))
      DO i = 1, SIZE(x0) - 6, 3
         x0(i:i + 2) = MATMUL(transform, x0(i:i + 2))
      END DO
      gopt_env%h_ref = reference_cell
   END SUBROUTINE bfgs_restore_cell_frame

! **************************************************************************************************
!> \brief ...
!> \param ndf ...
!> \param dg ...
!> \param eigval ...
!> \param work ...
!> \param eigvec_mat ...
!> \param g ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE rat_fun_opt(ndf, dg, eigval, work, eigvec_mat, g, para_env)

      INTEGER, INTENT(IN)                                :: ndf
      REAL(KIND=dp), INTENT(INOUT)                       :: dg(ndf), eigval(ndf), work(ndf)
      TYPE(cp_fm_type), INTENT(IN)                       :: eigvec_mat
      REAL(KIND=dp), INTENT(INOUT)                       :: g(ndf)
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rat_fun_opt'
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp

      INTEGER                                            :: handle, i, indf, iref, iter, j, k, l, &
                                                            maxit, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: bisec, conv, fail, set
      REAL(KIND=dp)                                      :: fun, fun1, fun2, fun3, fung, lam1, lam2, &
                                                            ln, lp, ssize, step, stol
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data

      CALL timeset(routineN, handle)

      stol = 1.0E-8_dp
      ssize = 0.2_dp
      maxit = 999
      fail = .FALSE.
      bisec = .FALSE.

      dg = 0._dp

      CALL cp_fm_get_info(eigvec_mat, row_indices=row_indices, col_indices=col_indices, &
                          local_data=local_data, nrow_local=nrow_local, ncol_local=ncol_local)

      DO i = 1, nrow_local
         j = row_indices(i)
         DO k = 1, ncol_local
            l = col_indices(k)
            dg(l) = dg(l) + local_data(i, k)*g(j)
         END DO
      END DO
      CALL para_env%sum(dg)

      set = .FALSE.

      DO

!   calculating Lambda

         lp = 0.0_dp
         iref = 1
         ln = 0.0_dp
         IF (eigval(iref) < 0.0_dp) ln = eigval(iref) - 0.01_dp

         conv = .FALSE.
         iter = 0
         DO
            iter = iter + 1
            fun = 0.0_dp
            fung = 0.0_dp
            DO indf = 1, ndf
               fun = fun + dg(indf)**2/(ln - eigval(indf))
               fung = fung - dg(indf)**2/((ln - eigval(indf))**2)
            END DO
            fun = fun - ln
            fung = fung - one
            step = fun/fung
            ln = ln - step
            IF (ABS(step) < stol) THEN
               conv = .TRUE.
               EXIT
            END IF
            IF (iter >= maxit) EXIT
         END DO
         outer: DO
            IF (.NOT. conv) THEN
               conv = .FALSE.
               bisec = .TRUE.
               iter = 0
               maxit = 9999
               lam1 = 0.0_dp
               IF (eigval(iref) < 0.0_dp) lam1 = eigval(iref) - 0.01_dp
               fun1 = 0.0_dp
               DO indf = 1, ndf
                  fun1 = fun1 + dg(indf)**2/(lam1 - eigval(indf))
               END DO
               fun1 = fun1 - lam1
               step = ABS(lam1)/1000.0_dp
               IF (step < ssize) step = ssize
               inner: DO
                  iter = iter + 1
                  IF (iter > maxit) THEN
                     ln = 0.0_dp
                     lp = 0.0_dp
                     fail = .TRUE.
                     EXIT outer
                  END IF
                  fun2 = 0.0_dp
                  lam2 = lam1 - iter*step
                  DO indf = 1, ndf
                     fun2 = fun2 + dg(indf)**2/(lam2 - eigval(indf))
                  END DO
                  fun2 = fun2 - lam2
                  IF (fun2*fun1 < 0.0_dp) THEN
                     iter = 0
                     DO
                        iter = iter + 1
                        IF (iter > maxit) THEN
                           ln = 0.0_dp
                           lp = 0.0_dp
                           fail = .TRUE.
                           EXIT outer
                        END IF
                        step = (lam1 + lam2)/2
                        fun3 = 0.0_dp
                        DO indf = 1, ndf
                           fun3 = fun3 + dg(indf)**2/(step - eigval(indf))
                        END DO
                        fun3 = fun3 - step

                        IF (ABS(step - lam2) < stol) THEN
                           ln = step
                           EXIT inner
                        END IF

                        IF (fun3*fun1 < stol) THEN
                           lam2 = step
                        ELSE
                           lam1 = step
                        END IF
                     END DO
                  END IF
               END DO inner
            END IF
            IF ((ln > eigval(iref)) .OR. ((ln > 0.0_dp) .AND. &
                                          (eigval(iref) > 0.0_dp))) THEN

               IF (.NOT. bisec) THEN
                  conv = .FALSE.
                  CYCLE outer
               END IF
               ln = 0.0_dp
               lp = 0.0_dp
               fail = .TRUE.
            END IF
            EXIT outer
         END DO outer

         IF (fail .AND. .NOT. set) THEN
            set = .TRUE.
            DO indf = 1, ndf
               eigval(indf) = eigval(indf)*work(indf)
            END DO
            CYCLE
         END IF

         IF (.NOT. set) THEN
            work(1:ndf) = one
         END IF

         DO indf = 1, ndf
            eigval(indf) = eigval(indf) - ln
         END DO
         EXIT
      END DO

      CALL timestop(handle)

   END SUBROUTINE rat_fun_opt

! **************************************************************************************************
!> \brief Update the direct Hessian with the damped BFGS formula
!> \param ndf Number of optimization variables
!> \param dx Change of the optimization variables
!> \param dg Change of the gradients in the same coordinate system
!> \param hess_mat Direct Hessian, updated in place
!> \param work Replicated workspace, holds B*dx on exit
!> \param para_env ...
!> \param output_unit Geometry log unit, activates at the default MEDIUM print level
!> \param hess_reset Set if the update had to be refused and the caller must
!>                  collectively rebuild the initial Hessian
! **************************************************************************************************
   SUBROUTINE bfgs(ndf, dx, dg, hess_mat, work, para_env, output_unit, hess_reset)
      INTEGER, INTENT(IN)                                :: ndf
      REAL(KIND=dp), INTENT(INOUT)                       :: dx(ndf), dg(ndf)
      TYPE(cp_fm_type), INTENT(IN)                       :: hess_mat
      REAL(KIND=dp), INTENT(INOUT)                       :: work(ndf)
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env
      INTEGER, INTENT(IN)                                :: output_unit
      LOGICAL, INTENT(OUT)                               :: hess_reset

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'bfgs'
      INTEGER, PARAMETER                                 :: exit_reset = 2, exit_skip = 1, &
                                                            exit_update = 0
      REAL(KIND=dp), PARAMETER                           :: min_curvature = 0.2_dp, one = 1.0_dp, &
                                                            min_sbs = 10.0_dp*TINY(one), &
                                                            zero = 0.0_dp

      INTEGER                                            :: exit_flag, handle, i, invalid, j, k, l, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: dxw, eta_i, eta_j, gdx, sbs, seta, sy, &
                                                            theta
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_hes

      CALL timeset(routineN, handle)
      hess_reset = .FALSE.

      CALL cp_fm_get_info(hess_mat, row_indices=row_indices, col_indices=col_indices, &
                          local_data=local_hes, nrow_local=nrow_local, ncol_local=ncol_local)

      ! Refuse corrupted input instead of poisoning the Hessian with Inf/NaN.
      IF (.NOT. ALL(ieee_is_finite(dx(1:ndf))) .OR. .NOT. ALL(ieee_is_finite(dg(1:ndf)))) THEN
         CPABORT("Non-finite displacement or gradient difference")
      END IF

      work = zero
      DO i = 1, nrow_local
         j = row_indices(i)
         DO k = 1, ncol_local
            l = col_indices(k)
            work(j) = work(j) + local_hes(i, k)*dx(l)
         END DO
      END DO

      IF (PRESENT(para_env)) CALL para_env%sum(work)

      sy = DOT_PRODUCT(dg(1:ndf), dx(1:ndf))
      sbs = DOT_PRODUCT(dx(1:ndf), work(1:ndf))

      ! One collective decision covers every early exit, so that all ranks
      ! take the same path. sbs and sy share the units of the energy, so the
      ! curvature floor is relative to sy. min_sbs only fences off values for
      ! which the divisions in the update below would overflow: with a
      ! smaller sbs, 0.2*sbs (the least possible seta) falls below the
      ! smallest normal number. The floor is never negative, so a non-positive
      ! sbs (lost positive definiteness) is caught as well.
      exit_flag = exit_update
      IF (DOT_PRODUCT(dx(1:ndf), dx(1:ndf)) == zero .OR. &
          DOT_PRODUCT(dg(1:ndf), dg(1:ndf)) == zero) THEN
         exit_flag = exit_skip
      ELSE IF (.NOT. ieee_is_finite(sy) .OR. .NOT. ieee_is_finite(sbs) .OR. &
               sbs <= MAX(EPSILON(one)*ABS(sy), min_sbs)) THEN
         exit_flag = exit_reset
      END IF
      ! The max folds the votes into the most conservative action all ranks
      ! can agree on. The exit codes are ordered by severity, so a skip vote
      ! on one rank and a reset vote on another reset both. Every rank takes
      ! the same branch, so the rebuild stays collective.
      IF (PRESENT(para_env)) CALL para_env%max(exit_flag)

      IF (exit_flag == exit_skip) THEN
         ! A zero displacement or gradient difference carries no curvature
         ! information, so there is nothing to learn from this step.
      ELSE IF (exit_flag == exit_reset) THEN
         ! The Hessian has lost positive definiteness along dx or the predicted
         ! curvature is not representable, so the secant pair cannot be damped:
         ! the caller restarts from the initial Hessian and skips this update.
         CALL cp_warn(__LOCATION__, &
                      "BFGS: vanishing or non-positive predicted curvature, resetting the Hessian")
         hess_reset = .TRUE.
      ELSE
         ! Damped update after M.J.D. Powell, Lecture Notes in Mathematics 630 (1978),
         ! eq. 3.7: whenever the measured curvature sy drops below 20% of the predicted
         ! curvature s.B.s, replace y by eta = theta*y + (1-theta)*B*s. The blended
         ! curvature satisfies eta.s = 0.2*s.B.s exactly, so the update stays positive
         ! definite in exact arithmetic, also for secant pairs that no line search
         ! checked. Because sbs passed the floor above, seta >= 0.2*sbs stays a
         ! normal number, so one/seta cannot overflow and no further floor on seta
         ! is needed.
         theta = one
         IF (sy < min_curvature*sbs) THEN
            theta = ((one - min_curvature)*sbs)/(sbs - sy)
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A,F10.4,A,F3.1,A)") &
                  "BFGS| Damped update, curvature ratio y.s/(s.B.s): ", sy/sbs, &
                  " (threshold ", min_curvature, ")"
               CALL m_flush(output_unit)
            END IF
         END IF
         seta = theta*sy + (one - theta)*sbs
         gdx = one/seta
         dxw = one/sbs

         DO i = 1, nrow_local
            j = row_indices(i)
            eta_i = theta*dg(j) + (one - theta)*work(j)
            DO k = 1, ncol_local
               l = col_indices(k)
               eta_j = theta*dg(l) + (one - theta)*work(l)
               local_hes(i, k) = local_hes(i, k) + gdx*eta_i*eta_j - &
                                 dxw*work(j)*work(l)
            END DO
         END DO

         ! Collective verdict: if any local block is corrupted, every rank resets.
         invalid = 0
         IF (.NOT. ALL(ieee_is_finite(local_hes(1:nrow_local, 1:ncol_local)))) invalid = 1
         IF (PRESENT(para_env)) CALL para_env%max(invalid)
         IF (invalid /= 0) THEN
            CALL cp_warn(__LOCATION__, "BFGS: non-finite Hessian after the update, resetting")
            hess_reset = .TRUE.
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE bfgs

! **************************************************************************************************
!> \brief Build the initial direct or inverse Hessian: sanitized model Hessian or unity
!> \param hess_mat BFGS matrix (direct or inverse Hessian), replaced in place
!> \param hess_tmp Eigensolver scratch matrix, same structure as hess_mat
!> \param eigvec Eigensolver eigenvector matrix, same structure as hess_mat
!> \param eigval Eigensolver workspace of length ndf
!> \param force_env Environment providing the model Hessian at the current geometry
!> \param use_model Build the sanitized model Hessian; otherwise fall back to unity
!> \param inverse Store the inverse of the sanitized model Hessian
!> \note Must be called by all ranks together: construct_initial_hess reduces internally.
! **************************************************************************************************
   SUBROUTINE build_initial_hess(hess_mat, hess_tmp, eigvec, eigval, force_env, use_model, inverse)
      TYPE(cp_fm_type), INTENT(IN)                       :: hess_mat, hess_tmp, eigvec
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigval
      TYPE(force_env_type), POINTER                      :: force_env
      LOGICAL, INTENT(IN)                                :: use_model, inverse

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_initial_hess'
      REAL(KIND=dp), PARAMETER                           :: min_eig = 0.1_dp, one = 1.0_dp, &
                                                            zero = 0.0_dp

      INTEGER                                            :: handle, indf, info, ndf

      CALL timeset(routineN, handle)

      IF (use_model) THEN
         CPASSERT(ASSOCIATED(force_env))
         ndf = SIZE(eigval)
         ! The raw model Hessian is not guaranteed positive definite, so
         ! eigenfactorize, floor the eigenvalues, and rebuild.
         CALL cp_fm_set_all(hess_mat, alpha=zero)
         CALL construct_initial_hess(force_env, hess_mat)
         CALL cp_fm_to_fm(hess_mat, hess_tmp)
         CALL choose_eigv_solver(hess_tmp, eigvec, eigval, info=info)
         IF (info /= 0) THEN
            ! In rare cases the diagonalization fails (bug in scalapack?)
            CALL cp_warn(__LOCATION__, "BFGS: model Hessian diagonalization failed, using unity")
            CALL cp_fm_set_all(hess_mat, alpha=zero, beta=one)
            eigval = one
         ELSE
            DO indf = 1, ndf
               IF (eigval(indf) < min_eig) eigval(indf) = min_eig
               IF (inverse) eigval(indf) = one/eigval(indf)
            END DO
            CALL cp_fm_to_fm(eigvec, hess_tmp)
            CALL cp_fm_column_scale(eigvec, eigval)
            CALL parallel_gemm("N", "T", ndf, ndf, ndf, one, hess_tmp, eigvec, zero, hess_mat)
         END IF
      ELSE
         ! No sanitized prior exists (user opt-out, cell optimization, or
         ! core-shell model): fall back to unity. eigval mirrors the unity
         ! spectrum, keeping the INTENT(OUT) argument defined on every path.
         CALL cp_fm_set_all(hess_mat, alpha=zero, beta=one)
         eigval = one
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_initial_hess

! **************************************************************************************************
!> \brief ...
!> \param ndf ...
!> \param eigval ...
!> \param work ...
! **************************************************************************************************
   SUBROUTINE set_hes_eig(ndf, eigval, work)
      INTEGER, INTENT(IN)                                :: ndf
      REAL(KIND=dp), INTENT(INOUT)                       :: eigval(ndf), work(ndf)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'set_hes_eig'
      REAL(KIND=dp), PARAMETER                           :: max_neg = -0.5_dp, max_pos = 5.0_dp, &
                                                            min_eig = 0.005_dp, one = 1.0_dp

      INTEGER                                            :: handle, indf
      LOGICAL                                            :: neg

      CALL timeset(routineN, handle)

      DO indf = 1, ndf
         IF (eigval(indf) < 0.0_dp) neg = .TRUE.
         IF (eigval(indf) > 1000.0_dp) eigval(indf) = 1000.0_dp
      END DO
      DO indf = 1, ndf
         IF (eigval(indf) < 0.0_dp) THEN
            IF (eigval(indf) < max_neg) THEN
               eigval(indf) = max_neg
            ELSE IF (eigval(indf) > -min_eig) THEN
               eigval(indf) = -min_eig
            END IF
         ELSE IF (eigval(indf) < 1000.0_dp) THEN
            IF (eigval(indf) < min_eig) THEN
               eigval(indf) = min_eig
            ELSE IF (eigval(indf) > max_pos) THEN
               eigval(indf) = max_pos
            END IF
         END IF
      END DO

      DO indf = 1, ndf
         IF (eigval(indf) < 0.0_dp) THEN
            work(indf) = -one
         ELSE
            work(indf) = one
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE set_hes_eig

! **************************************************************************************************
!> \brief ...
!> \param ndf ...
!> \param eigval ...
!> \param eigvec_mat ...
!> \param hess_tmp ...
!> \param dr ...
!> \param g ...
!> \param para_env ...
!> \param use_rfo ...
! **************************************************************************************************
   SUBROUTINE geoopt_get_step(ndf, eigval, eigvec_mat, hess_tmp, dr, g, para_env, use_rfo)

      INTEGER, INTENT(IN)                                :: ndf
      REAL(KIND=dp), INTENT(INOUT)                       :: eigval(ndf)
      TYPE(cp_fm_type), INTENT(IN)                       :: eigvec_mat, hess_tmp
      REAL(KIND=dp), INTENT(INOUT)                       :: dr(ndf), g(ndf)
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env
      LOGICAL                                            :: use_rfo

      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      INTEGER                                            :: i, indf, j, k, l, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type)                                   :: tmp

      CALL cp_fm_to_fm(eigvec_mat, hess_tmp)
      IF (use_rfo) THEN
         DO indf = 1, ndf
            eigval(indf) = one/eigval(indf)
         END DO
      ELSE
         DO indf = 1, ndf
            eigval(indf) = one/MAX(min_hessian_eigenvalue, eigval(indf))
         END DO
      END IF

      CALL cp_fm_column_scale(hess_tmp, eigval)
      CALL cp_fm_get_info(eigvec_mat, matrix_struct=matrix_struct)
      CALL cp_fm_create(tmp, matrix_struct, name="tmp")
      CALL cp_fm_set_all(tmp, alpha=zero)

      CALL parallel_gemm("N", "T", ndf, ndf, ndf, one, hess_tmp, eigvec_mat, zero, tmp)

      CALL cp_fm_transpose(tmp, hess_tmp)
      CALL cp_fm_release(tmp)

      ! New step

      CALL cp_fm_get_info(hess_tmp, row_indices=row_indices, col_indices=col_indices, &
                          local_data=local_data, nrow_local=nrow_local, ncol_local=ncol_local)

      dr = 0.0_dp
      DO i = 1, nrow_local
         j = row_indices(i)
         DO k = 1, ncol_local
            l = col_indices(k)
            dr(j) = dr(j) - local_data(i, k)*g(l)
         END DO
      END DO

      CALL para_env%sum(dr)

   END SUBROUTINE geoopt_get_step

! **************************************************************************************************
!> \brief ...
!> \param ndf ...
!> \param step ...
!> \param rad ...
!> \param rat ...
!> \param dr ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE trust_radius(ndf, step, rad, rat, dr, output_unit)
      INTEGER, INTENT(IN)                                :: ndf
      REAL(KIND=dp), INTENT(INOUT)                       :: step, rad, rat, dr(ndf)
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'trust_radius'
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: scal

      CALL timeset(routineN, handle)

      step = MAXVAL(ABS(dr))
      scal = MAX(one, rad/step)

      IF (step > rad) THEN
         rat = rad/step
         CALL DSCAL(ndf, rat, dr, 1)
         step = rad
         IF (output_unit > 0) THEN
            WRITE (unit=output_unit, FMT="(/,T2,A,F8.5)") &
               " Step is scaled; Scaling factor = ", rat
            CALL m_flush(output_unit)
         END IF
      END IF
      CALL timestop(handle)

   END SUBROUTINE trust_radius

! **************************************************************************************************
!> \brief ...
!> \param ndf ...
!> \param work ...
!> \param hess_mat ...
!> \param dr ...
!> \param g ...
!> \param conv ...
!> \param pred ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE energy_predict(ndf, work, hess_mat, dr, g, conv, pred, para_env)

      INTEGER, INTENT(IN)                                :: ndf
      REAL(KIND=dp), INTENT(INOUT)                       :: work(ndf)
      TYPE(cp_fm_type), INTENT(IN)                       :: hess_mat
      REAL(KIND=dp), INTENT(INOUT)                       :: dr(ndf), g(ndf)
      LOGICAL, INTENT(INOUT)                             :: conv
      REAL(KIND=dp), INTENT(INOUT)                       :: pred
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'energy_predict'
      REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

      INTEGER                                            :: handle, i, j, k, l, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: ener1, ener2
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data

      CALL timeset(routineN, handle)

      ener1 = DOT_PRODUCT(g(1:ndf), dr(1:ndf))

      CALL cp_fm_get_info(hess_mat, row_indices=row_indices, col_indices=col_indices, &
                          local_data=local_data, nrow_local=nrow_local, ncol_local=ncol_local)

      work = zero
      DO i = 1, nrow_local
         j = row_indices(i)
         DO k = 1, ncol_local
            l = col_indices(k)
            work(j) = work(j) + local_data(i, k)*dr(l)
         END DO
      END DO

      CALL para_env%sum(work)
      ener2 = DOT_PRODUCT(dr(1:ndf), work(1:ndf))
      pred = ener1 + 0.5_dp*ener2
      conv = .FALSE.
      CALL timestop(handle)

   END SUBROUTINE energy_predict

! **************************************************************************************************
!> \brief ...
!> \param rat ...
!> \param rad ...
!> \param step ...
!> \param ediff ...
! **************************************************************************************************
   SUBROUTINE update_trust_rad(rat, rad, step, ediff)

      REAL(KIND=dp), INTENT(INOUT)                       :: rat, rad, step, ediff

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'update_trust_rad'
      REAL(KIND=dp), PARAMETER                           :: max_trust = 1.0_dp, min_trust = 0.1_dp

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (rat > 4.0_dp) THEN
         IF (ediff < 0.0_dp) THEN
            rad = step*0.5_dp
         ELSE
            rad = step*0.25_dp
         END IF
      ELSE IF (rat > 2.0_dp) THEN
         IF (ediff < 0.0_dp) THEN
            rad = step*0.75_dp
         ELSE
            rad = step*0.5_dp
         END IF
      ELSE IF (rat > 4.0_dp/3.0_dp) THEN
         IF (ediff < 0.0_dp) THEN
            rad = step
         ELSE
            rad = step*0.75_dp
         END IF
      ELSE IF (rat > 10.0_dp/9.0_dp) THEN
         IF (ediff < 0.0_dp) THEN
            rad = step*1.25_dp
         ELSE
            rad = step
         END IF
      ELSE IF (rat > 0.9_dp) THEN
         IF (ediff < 0.0_dp) THEN
            rad = step*1.5_dp
         ELSE
            rad = step*1.25_dp
         END IF
      ELSE IF (rat > 0.75_dp) THEN
         IF (ediff < 0.0_dp) THEN
            rad = step*1.25_dp
         ELSE
            rad = step
         END IF
      ELSE IF (rat > 0.5_dp) THEN
         IF (ediff < 0.0_dp) THEN
            rad = step
         ELSE
            rad = step*0.75_dp
         END IF
      ELSE IF (rat > 0.25_dp) THEN
         IF (ediff < 0.0_dp) THEN
            rad = step*0.75_dp
         ELSE
            rad = step*0.5_dp
         END IF
      ELSE IF (ediff < 0.0_dp) THEN
         rad = step*0.5_dp
      ELSE
         rad = step*0.25_dp
      END IF

      rad = MAX(rad, min_trust)
      rad = MIN(rad, max_trust)
      CALL timestop(handle)

   END SUBROUTINE update_trust_rad

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param geo_section ...
!> \param hess_mat ...
!> \param logger ...
!> \param gopt_env Optimization coordinate system
!> \param matrix_method Stored matrix representation
! **************************************************************************************************
   SUBROUTINE write_bfgs_hessian(geo_section, hess_mat, logger, gopt_env, matrix_method)

      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(cp_fm_type), INTENT(IN)                       :: hess_mat
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      INTEGER, INTENT(IN)                                :: matrix_method

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_bfgs_hessian'

      INTEGER                                            :: handle, hesunit, ndf
      REAL(KIND=dp), DIMENSION(3, 3)                     :: cell_frame

      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(hess_mat, nrow_global=ndf)
      cell_frame = 0.0_dp
      IF (gopt_env%type_id == default_cell_method_id) cell_frame = gopt_env%cell_env%opt_to_input
      hesunit = cp_print_key_unit_nr(logger, geo_section, "BFGS%RESTART", &
                                     extension=".Hessian", file_form="UNFORMATTED", file_action="WRITE", &
                                     file_position="REWIND")
      IF (hesunit > 0) THEN
         WRITE (hesunit) bfgs_restart_magic, bfgs_restart_version, matrix_method, ndf, &
            gopt_env%type_id, gopt_env%h_ref, cell_frame
      END IF
      ! Store the native representation. Conversions occur only when reading a different METHOD.
      CALL cp_fm_write_unformatted(hess_mat, hesunit)
      CALL cp_print_key_finished_output(hesunit, logger, geo_section, "BFGS%RESTART")
      CALL timestop(handle)

   END SUBROUTINE write_bfgs_hessian

! **************************************************************************************************
!> \brief Constructs model Hessian as described in https://doi.org/10.1016/0009-2614(95)00646-L.
!> \param force_env ...
!> \param hess_mat ...
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE construct_initial_hess(force_env, hess_mat)

      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(cp_fm_type), INTENT(IN)                       :: hess_mat

      INTEGER                                            :: i, iat_col, iat_row, iglobal, iind, j, &
                                                            jat_row, jglobal, jind, k, natom, &
                                                            ncol_local, nrow_local, z
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: at_row
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_ij, rho_ij
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: r_ij
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: fixed, local_data
      REAL(KIND=dp), DIMENSION(3, 3)                     :: alpha, r0
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(particle_list_type), POINTER                  :: particles

      CALL cite_reference(Lindh1995)

      CALL force_env_get(force_env=force_env, subsys=subsys, cell=cell)
      CALL cp_subsys_get(subsys, &
                         particles=particles)

      alpha(1, :) = [1._dp, 0.3949_dp, 0.3949_dp]
      alpha(2, :) = [0.3494_dp, 0.2800_dp, 0.2800_dp]
      alpha(3, :) = [0.3494_dp, 0.2800_dp, 0.1800_dp]

      r0(1, :) = [1.35_dp, 2.10_dp, 2.53_dp]
      r0(2, :) = [2.10_dp, 2.87_dp, 3.40_dp]
      r0(3, :) = [2.53_dp, 3.40_dp, 3.40_dp]

      CALL cp_fm_get_info(hess_mat, row_indices=row_indices, col_indices=col_indices, &
                          local_data=local_data, nrow_local=nrow_local, ncol_local=ncol_local)
      natom = particles%n_els
      ALLOCATE (at_row(natom))
      ALLOCATE (rho_ij(natom, natom))
      ALLOCATE (d_ij(natom, natom))
      ALLOCATE (r_ij(natom, natom, 3))
      ALLOCATE (fixed(3, natom))
      fixed = 1.0_dp
      CALL fix_atom_control(force_env, fixed)
      DO i = 1, 3
         CALL hess_mat%matrix_struct%para_env%min(fixed(i, :))
      END DO
      rho_ij = 0
      !XXXX insert proper rows !XXX
      at_row = 3
      DO i = 1, natom
         CALL get_atomic_kind(atomic_kind=particles%els(i)%atomic_kind, z=z)
         IF (z <= 10) at_row(i) = 2
         IF (z <= 2) at_row(i) = 1
      END DO
      DO i = 2, natom
         iat_row = at_row(i)
         DO j = 1, i - 1
            jat_row = at_row(j)
            !pbc for a distance vector
            r_ij(j, i, :) = pbc(particles%els(i)%r, particles%els(j)%r, cell)
            r_ij(i, j, :) = -r_ij(j, i, :)
            d_ij(j, i) = NORM2(r_ij(j, i, :))
            d_ij(i, j) = d_ij(j, i)
            rho_ij(j, i) = EXP(alpha(jat_row, iat_row)*(r0(jat_row, iat_row)**2 - d_ij(j, i)**2))
            rho_ij(i, j) = rho_ij(j, i)
         END DO
      END DO
      DO i = 1, ncol_local
         iglobal = col_indices(i)
         iind = MOD(iglobal - 1, 3) + 1
         iat_col = (iglobal + 2)/3
         IF (iat_col > natom) CYCLE
         DO j = 1, nrow_local
            jglobal = row_indices(j)
            jind = MOD(jglobal - 1, 3) + 1
            iat_row = (jglobal + 2)/3
            IF (iat_row > natom) CYCLE
            IF (iat_row /= iat_col) THEN
               IF (d_ij(iat_row, iat_col) < 6.0_dp) THEN
                  local_data(j, i) = local_data(j, i) + &
                                     angle_second_deriv(r_ij, d_ij, rho_ij, iind, jind, iat_col, iat_row, natom)
               END IF
            ELSE
               local_data(j, i) = local_data(j, i) + &
                                  angle_second_deriv(r_ij, d_ij, rho_ij, iind, jind, iat_col, iat_row, natom)
            END IF
            IF (iat_col /= iat_row) THEN
               IF (d_ij(iat_row, iat_col) < 6.0_dp) THEN
                  local_data(j, i) = local_data(j, i) - &
                                     dist_second_deriv(r_ij(iat_col, iat_row, :), &
                                                       iind, jind, d_ij(iat_row, iat_col), rho_ij(iat_row, iat_col))
               END IF
            ELSE
               DO k = 1, natom
                  IF (k == iat_col) CYCLE
                  IF (d_ij(iat_row, k) < 6.0_dp) THEN
                     local_data(j, i) = local_data(j, i) + &
                                        dist_second_deriv(r_ij(iat_col, k, :), &
                                                          iind, jind, d_ij(iat_row, k), rho_ij(iat_row, k))
                  END IF
               END DO
            END IF
            IF (fixed(jind, iat_row) < 0.5_dp .OR. fixed(iind, iat_col) < 0.5_dp) THEN
               local_data(j, i) = 0.0_dp
               IF (jind == iind .AND. iat_row == iat_col) local_data(j, i) = 1.0_dp
            END IF
         END DO
      END DO
      DEALLOCATE (fixed)
      DEALLOCATE (rho_ij)
      DEALLOCATE (d_ij)
      DEALLOCATE (r_ij)
      DEALLOCATE (at_row)

   END SUBROUTINE construct_initial_hess

! **************************************************************************************************
!> \brief ...
!> \param r1 ...
!> \param i ...
!> \param j ...
!> \param d ...
!> \param rho ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION dist_second_deriv(r1, i, j, d, rho) RESULT(deriv)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: r1
      INTEGER, INTENT(IN)                                :: i, j
      REAL(KIND=dp), INTENT(IN)                          :: d, rho
      REAL(KIND=dp)                                      :: deriv

      deriv = 0.45_dp*rho*(r1(i)*r1(j))/d**2
   END FUNCTION dist_second_deriv

! **************************************************************************************************
!> \brief ...
!> \param r_ij ...
!> \param d_ij ...
!> \param rho_ij ...
!> \param idir ...
!> \param jdir ...
!> \param iat_der ...
!> \param jat_der ...
!> \param natom ...
!> \return ...
! **************************************************************************************************
   FUNCTION angle_second_deriv(r_ij, d_ij, rho_ij, idir, jdir, iat_der, jat_der, natom) RESULT(deriv)
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: r_ij
      REAL(KIND=dp), DIMENSION(:, :)                     :: d_ij, rho_ij
      INTEGER                                            :: idir, jdir, iat_der, jat_der, natom
      REAL(KIND=dp)                                      :: deriv

      REAL(KIND=dp), PARAMETER                           :: rho_thre = 1.0E-5_dp

      INTEGER                                            :: i, iat, idr, j, jat, jdr
      REAL(KIND=dp)                                      :: d12, d23, d31, D_mat(3, 2), denom1, &
                                                            denom2, denom3, ka1, ka2, ka3, rho12, &
                                                            rho23, rho31, rsst1, rsst2, rsst3
      REAL(KIND=dp), DIMENSION(3)                        :: r12, r23, r31

      deriv = 0._dp
      IF (iat_der == jat_der) THEN
         DO i = 1, natom - 1
            IF (rho_ij(iat_der, i) < rho_thre) CYCLE
            DO j = i + 1, natom
               IF (rho_ij(iat_der, j) < rho_thre) CYCLE
               IF (i == iat_der .OR. j == iat_der) CYCLE
               IF (iat_der < i .OR. iat_der > j) THEN
                  r12 = r_ij(iat_der, i, :); r23 = r_ij(i, j, :); r31 = r_ij(j, iat_der, :)
                  d12 = d_ij(iat_der, i); d23 = d_ij(i, j); d31 = d_ij(j, iat_der)
                  rho12 = rho_ij(iat_der, i); rho23 = rho_ij(i, j); rho31 = rho_ij(j, iat_der)
               ELSE
                  r12 = r_ij(iat_der, j, :); r23 = r_ij(j, i, :); r31 = r_ij(i, iat_der, :)
                  d12 = d_ij(iat_der, j); d23 = d_ij(j, i); d31 = d_ij(i, iat_der)
                  rho12 = rho_ij(iat_der, j); rho23 = rho_ij(j, i); rho31 = rho_ij(i, iat_der)
               END IF
               ka1 = 0.15_dp*rho12*rho23; ka2 = 0.15_dp*rho23*rho31; ka3 = 0.15_dp*rho31*rho12
               rsst1 = DOT_PRODUCT(r12, r23); rsst2 = DOT_PRODUCT(r23, r31); rsst3 = DOT_PRODUCT(r31, r12)
               denom1 = 1.0_dp - rsst1**2/(d12**2*d23**2); denom2 = 1.0_dp - rsst2**2/(d23**2*d31**2)
               denom3 = 1.0_dp - rsst3**2/(d31**2*d12**2)
               denom1 = SIGN(1.0_dp, denom1)*MAX(ABS(denom1), 0.01_dp)
               denom2 = SIGN(1.0_dp, denom2)*MAX(ABS(denom2), 0.01_dp)
               denom3 = SIGN(1.0_dp, denom3)*MAX(ABS(denom3), 0.01_dp)
               D_mat(1, 1) = r23(idir)/(d12*d23) - rsst1*r12(idir)/(d12**3*d23)
               D_mat(1, 2) = r23(jdir)/(d12*d23) - rsst1*r12(jdir)/(d12**3*d23)
               D_mat(2, 1) = -r23(idir)/(d23*d31) + rsst2*r31(idir)/(d23*d31**3)
               D_mat(2, 2) = -r23(jdir)/(d23*d31) + rsst2*r31(jdir)/(d23*d31**3)
               D_mat(3, 1) = (r31(idir) - r12(idir))/(d31*d12) + rsst3*r31(idir)/(d31**3*d12) - &
                             rsst3*r12(idir)/(d31*d12**3)
               D_mat(3, 2) = (r31(jdir) - r12(jdir))/(d31*d12) + rsst3*r31(jdir)/(d31**3*d12) - &
                             rsst3*r12(jdir)/(d31*d12**3)
               IF (ABS(denom1) <= 0.011_dp) D_mat(1, 1) = 0.0_dp
               IF (ABS(denom2) <= 0.011_dp) D_mat(2, 1) = 0.0_dp
               IF (ABS(denom3) <= 0.011_dp) D_mat(3, 1) = 0.0_dp
               deriv = deriv + ka1*D_mat(1, 1)*D_mat(1, 2)/denom1 + &
                       ka2*D_mat(2, 1)*D_mat(2, 2)/denom2 + &
                       ka3*D_mat(3, 1)*D_mat(3, 2)/denom3

            END DO
         END DO
      ELSE
         DO i = 1, natom
            IF (i == iat_der .OR. i == jat_der) CYCLE
            IF (jat_der < iat_der) THEN
               iat = jat_der; jat = iat_der; idr = jdir; jdr = idir
            ELSE
               iat = iat_der; jat = jat_der; idr = idir; jdr = jdir
            END IF
            IF (jat < i .OR. iat > i) THEN
               r12 = r_ij(iat, jat, :); r23 = r_ij(jat, i, :); r31 = r_ij(i, iat, :)
               d12 = d_ij(iat, jat); d23 = d_ij(jat, i); d31 = d_ij(i, iat)
               rho12 = rho_ij(iat, jat); rho23 = rho_ij(jat, i); rho31 = rho_ij(i, iat)
            ELSE
               r12 = r_ij(iat, i, :); r23 = r_ij(i, jat, :); r31 = r_ij(jat, iat, :)
               d12 = d_ij(iat, i); d23 = d_ij(i, jat); d31 = d_ij(jat, iat)
               rho12 = rho_ij(iat, i); rho23 = rho_ij(i, jat); rho31 = rho_ij(jat, iat)
            END IF
            ka1 = 0.15_dp*rho12*rho23; ka2 = 0.15_dp*rho23*rho31; ka3 = 0.15_dp*rho31*rho12
            rsst1 = DOT_PRODUCT(r12, r23); rsst2 = DOT_PRODUCT(r23, r31); rsst3 = DOT_PRODUCT(r31, r12)
            denom1 = 1.0_dp - rsst1**2/(d12**2*d23**2); denom2 = 1.0_dp - rsst2**2/(d23**2*d31**2)
            denom3 = 1.0_dp - rsst3**2/(d31**2*d12**2)
            denom1 = SIGN(1.0_dp, denom1)*MAX(ABS(denom1), 0.01_dp)
            denom2 = SIGN(1.0_dp, denom2)*MAX(ABS(denom2), 0.01_dp)
            denom3 = SIGN(1.0_dp, denom3)*MAX(ABS(denom3), 0.01_dp)
            D_mat(1, 1) = r23(idr)/(d12*d23) - rsst1*r12(idr)/(d12**3*d23)
            D_mat(2, 1) = -r23(idr)/(d23*d31) + rsst2*r31(idr)/(d23*d31**3)
            D_mat(3, 1) = (r31(idr) - r12(idr))/(d31*d12) + rsst3*r31(idr)/(d31**3*d12) - &
                          rsst3*r12(idr)/(d31*d12**3)
            IF (jat < i .OR. iat > i) THEN
               D_mat(1, 2) = (r12(jdr) - r23(jdr))/(d12*d23) + rsst1*r12(jdr)/(d12**3*d23) - &
                             rsst1*r23(jdr)/(d12*d23**3)
               D_mat(2, 2) = r31(jdr)/(d23*d31) - rsst2*r23(jdr)/(d23**3*d31)
               D_mat(3, 2) = -r31(jdr)/(d31*d12) + rsst3*r12(jdr)/(d31*d12**3)
            ELSE
               D_mat(1, 2) = -r12(jdr)/(d12*d23) + rsst1*r23(jdr)/(d12*d23**3)
               D_mat(2, 2) = (r23(jdr) - r31(jdr))/(d23*d31) + rsst2*r23(jdr)/(d23**3*d31) - &
                             rsst2*r31(jdr)/(d23*d31**3)
               D_mat(3, 2) = r12(jdr)/(d31*d12) - rsst3*r31(jdr)/(d31**3*d12)
            END IF
            IF (ABS(denom1) <= 0.011_dp) D_mat(1, 1) = 0.0_dp
            IF (ABS(denom2) <= 0.011_dp) D_mat(2, 1) = 0.0_dp
            IF (ABS(denom3) <= 0.011_dp) D_mat(3, 1) = 0.0_dp

            deriv = deriv + ka1*D_mat(1, 1)*D_mat(1, 2)/denom1 + &
                    ka2*D_mat(2, 1)*D_mat(2, 2)/denom2 + &
                    ka3*D_mat(3, 1)*D_mat(3, 2)/denom3
         END DO
      END IF
      deriv = 0.25_dp*deriv

   END FUNCTION angle_second_deriv

END MODULE bfgs_optimizer
