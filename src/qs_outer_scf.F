!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for performing an outer scf loop
!> \par History
!>      Created [2006.03]
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE qs_outer_scf
   USE cp_control_types,                ONLY: becke_restraint_type,&
                                              cdft_control_type,&
                                              ddapc_restraint_type,&
                                              dft_control_type,&
                                              s2_restraint_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE input_constants,                 ONLY: &
        broyden_type_1, broyden_type_1_explicit, broyden_type_1_explicit_ls, broyden_type_1_ls, &
        broyden_type_2, broyden_type_2_explicit, broyden_type_2_explicit_ls, broyden_type_2_ls, &
        cdft2ot, cdft_combined_constraint, cdft_density_constraint, cdft_magnetization_constraint, &
        do_ddapc_constraint, do_s2_constraint, ot2cdft, outer_scf_basis_center_opt, &
        outer_scf_becke_constraint, outer_scf_cdft_constraint, outer_scf_ddapc_constraint, &
        outer_scf_none, outer_scf_optimizer_bisect, outer_scf_optimizer_broyden, &
        outer_scf_optimizer_diis, outer_scf_optimizer_newton, outer_scf_optimizer_newton_ls, &
        outer_scf_optimizer_none, outer_scf_optimizer_sd, outer_scf_optimizer_secant, &
        outer_scf_s2_constraint
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE qs_basis_gradient,               ONLY: qs_basis_center_gradient,&
                                              qs_update_basis_center_pos,&
                                              return_basis_center_gradient_norm
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_outer_scf'

! *** Public subroutines ***

   PUBLIC :: outer_loop_gradient, outer_loop_optimize, outer_loop_update_qs_env, &
             outer_loop_variables_count, outer_loop_extrapolate, &
             outer_loop_switch, outer_loop_purge_history

CONTAINS

! **************************************************************************************************
!> \brief returns the number of variables that is employed in the outer loop. with a CDFT constraint
!>        this value is returned by the cdft_control type
!> \param scf_control the outer loop control type
!> \param cdft_control the cdft loop control type
!> \retval res the number of variables
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   FUNCTION outer_loop_variables_count(scf_control, cdft_control) RESULT(res)
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(cdft_control_type), INTENT(IN), OPTIONAL, &
         POINTER                                         :: cdft_control
      INTEGER                                            :: res

      SELECT CASE (scf_control%outer_scf%type)
      CASE (outer_scf_ddapc_constraint)
         res = 1
      CASE (outer_scf_s2_constraint)
         res = 1
      CASE (outer_scf_becke_constraint, outer_scf_cdft_constraint)
         IF (PRESENT(cdft_control)) THEN
            SELECT CASE (cdft_control%constraint_type)
            CASE (cdft_density_constraint, cdft_magnetization_constraint)
               res = 1
            CASE (cdft_combined_constraint)
               res = 2
            END SELECT
         ELSE
            res = 1
         END IF
      CASE (outer_scf_basis_center_opt)
         res = 1
      CASE (outer_scf_none) ! just needed to communicate the gradient criterium
         res = 1
      CASE DEFAULT
         res = 0
      END SELECT

   END FUNCTION outer_loop_variables_count

! **************************************************************************************************
!> \brief computes the gradient wrt to the outer loop variables
!> \param qs_env ...
!> \param scf_env ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE outer_loop_gradient(qs_env, scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_gradient', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ihistory, ivar, n
      LOGICAL                                            :: is_constraint
      TYPE(becke_restraint_type), POINTER                :: becke_control
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(ddapc_restraint_type), POINTER                :: ddapc_restraint_control
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(s2_restraint_type), POINTER                   :: s2_restraint_control
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, scf_control=scf_control, &
                      dft_control=dft_control, energy=energy)
      CPASSERT(scf_control%outer_scf%have_scf)

      ihistory = scf_env%outer_scf%iter_count
      CPASSERT(ihistory <= SIZE(scf_env%outer_scf%energy, 1))

      scf_env%outer_scf%energy(ihistory) = energy%total

      SELECT CASE (scf_control%outer_scf%type)
      CASE (outer_scf_none)
         ! just pass the inner loop scf criterium to the outer loop one
         scf_env%outer_scf%variables(1, ihistory) = scf_env%iter_delta
         scf_env%outer_scf%gradient(1, ihistory) = scf_env%iter_delta
      CASE (outer_scf_ddapc_constraint)
         CPASSERT(dft_control%qs_control%ddapc_restraint)
         DO n = 1, SIZE(dft_control%qs_control%ddapc_restraint_control)
            NULLIFY (ddapc_restraint_control)
            ddapc_restraint_control => dft_control%qs_control%ddapc_restraint_control(n)%ddapc_restraint_control
            is_constraint = (ddapc_restraint_control%functional_form == do_ddapc_constraint)
            IF (is_constraint) EXIT
         END DO
         CPASSERT(is_constraint)

         scf_env%outer_scf%variables(:, ihistory) = ddapc_restraint_control%strength
         scf_env%outer_scf%gradient(:, ihistory) = ddapc_restraint_control%ddapc_order_p- &
                                                   ddapc_restraint_control%target
      CASE (outer_scf_s2_constraint)
         CPASSERT(dft_control%qs_control%s2_restraint)
         s2_restraint_control => dft_control%qs_control%s2_restraint_control
         is_constraint = (s2_restraint_control%functional_form == do_s2_constraint)
         CPASSERT(is_constraint)

         scf_env%outer_scf%variables(:, ihistory) = s2_restraint_control%strength
         scf_env%outer_scf%gradient(:, ihistory) = s2_restraint_control%s2_order_p- &
                                                   s2_restraint_control%target
      CASE (outer_scf_becke_constraint)
         CPASSERT(dft_control%qs_control%becke_restraint)
         becke_control => dft_control%qs_control%becke_control
         DO ivar = 1, SIZE(scf_env%outer_scf%gradient, 1)
            scf_env%outer_scf%variables(ivar, ihistory) = becke_control%strength(ivar)
            scf_env%outer_scf%gradient(ivar, ihistory) = becke_control%becke_order_p(ivar)- &
                                                         becke_control%target(ivar)
         END DO
      CASE (outer_scf_cdft_constraint)
         CPASSERT(dft_control%qs_control%cdft)
         cdft_control => dft_control%qs_control%cdft_control
         DO ivar = 1, SIZE(scf_env%outer_scf%gradient, 1)
            scf_env%outer_scf%variables(ivar, ihistory) = cdft_control%strength(ivar)
            scf_env%outer_scf%gradient(ivar, ihistory) = cdft_control%value(ivar)- &
                                                         cdft_control%target(ivar)
         END DO
      CASE (outer_scf_basis_center_opt)
         CALL qs_basis_center_gradient(qs_env)
         scf_env%outer_scf%gradient(:, ihistory) = return_basis_center_gradient_norm(qs_env)

      CASE DEFAULT
         CPABORT("")

      END SELECT

      CALL timestop(handle)

   END SUBROUTINE outer_loop_gradient

! **************************************************************************************************
!> \brief optimizes the parameters of the outer_scf
!> \param scf_env the scf_env where to optimize the parameters
!> \param scf_control control parameters for the optimization
!> \par History
!>      03.2006 created [Joost VandeVondele]
!>      01.2017 added Broyden and Newton optimizers [Nico Holmberg]
!> \note
!>       ought to be general, and independent of the actual kind of variables
! **************************************************************************************************
   SUBROUTINE outer_loop_optimize(scf_env, scf_control)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_optimize', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ibuf, ihigh, ihistory, ilow, &
                                                            j, jbuf, nb, nvar, optimizer_type
      REAL(KIND=dp)                                      :: interval, scale, tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ev
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: a, b, f, x
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: inv_jacobian

      CALL timeset(routineN, handle)

      ihistory = scf_env%outer_scf%iter_count
      optimizer_type = scf_control%outer_scf%optimizer
      NULLIFY (inv_jacobian)

      IF (scf_control%outer_scf%type == outer_scf_basis_center_opt) THEN
         scf_env%outer_scf%variables(:, ihistory+1) = scf_env%outer_scf%variables(:, ihistory)
      ELSE
         DO WHILE (.TRUE.) ! if we need a different run type we'll restart here

            SELECT CASE (optimizer_type)
            CASE (outer_scf_optimizer_bisect) ! bisection on the gradient, needs to be 1D
               CPASSERT(SIZE(scf_env%outer_scf%gradient(:, 1)) == 1)
               ! find the pair of points that bracket a zero of the gradient, with the smallest interval possible
               ilow = -1
               ihigh = -1
               interval = HUGE(interval)
               DO i = 1, ihistory
                  DO j = i+1, ihistory
                     ! distrust often used points
                     IF (scf_env%outer_scf%count(i) .GT. scf_control%outer_scf%bisect_trust_count) CYCLE
                     IF (scf_env%outer_scf%count(j) .GT. scf_control%outer_scf%bisect_trust_count) CYCLE

                     ! if they bracket a zero use them
                     IF (scf_env%outer_scf%gradient(1, i)* &
                         scf_env%outer_scf%gradient(1, j) < 0.0_dp) THEN
                        tmp = ABS(scf_env%outer_scf%variables(1, i)-scf_env%outer_scf%variables(1, j))
                        IF (tmp < interval) THEN
                           ilow = i
                           ihigh = j
                           interval = tmp
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
               IF (ilow == -1) THEN ! we didn't bracket a minimum yet, try something else
                  optimizer_type = outer_scf_optimizer_diis
                  CYCLE
               ENDIF
               scf_env%outer_scf%count(ilow) = scf_env%outer_scf%count(ilow)+1
               scf_env%outer_scf%count(ihigh) = scf_env%outer_scf%count(ihigh)+1
               scf_env%outer_scf%variables(:, ihistory+1) = 0.5_dp*(scf_env%outer_scf%variables(:, ilow)+ &
                                                                    scf_env%outer_scf%variables(:, ihigh))
            CASE (outer_scf_optimizer_none)
               scf_env%outer_scf%variables(:, ihistory+1) = scf_env%outer_scf%variables(:, ihistory)
            CASE (outer_scf_optimizer_sd)
               ! Notice that we are just trying to find a stationary point
               ! e.g. the ddpac_constraint, one maximizes the function, so the stepsize might have
               ! to be negative
               scf_env%outer_scf%variables(:, ihistory+1) = scf_env%outer_scf%variables(:, ihistory)- &
                                                            scf_control%outer_scf%step_size*scf_env%outer_scf%gradient(:, ihistory)
            CASE (outer_scf_optimizer_diis)
               CPASSERT(scf_control%outer_scf%diis_buffer_length > 0)
               ! set up DIIS matrix
               nb = MIN(ihistory, scf_control%outer_scf%diis_buffer_length)
               IF (nb < 2) THEN
                  optimizer_type = outer_scf_optimizer_sd
                  CYCLE
               ELSE
                  ALLOCATE (b(nb+1, nb+1), a(nb+1, nb+1), ev(nb+1))
                  DO I = 1, nb
                     DO J = I, nb
                        ibuf = ihistory-nb+i
                        jbuf = ihistory-nb+j
                        b(I, J) = DOT_PRODUCT(scf_env%outer_scf%gradient(:, ibuf), &
                                              scf_env%outer_scf%gradient(:, jbuf))
                        b(J, I) = b(I, J)
                     ENDDO
                  ENDDO
                  b(nb+1, :) = -1.0_dp
                  b(:, nb+1) = -1.0_dp
                  b(nb+1, nb+1) = 0.0_dp

                  CALL diamat_all(b, ev)
                  a(:, :) = b
                  DO I = 1, nb+1
                     IF (ABS(ev(I)) .LT. 1.0E-12_dp) THEN
                        a(:, I) = 0.0_dp
                     ELSE
                        a(:, I) = a(:, I)/ev(I)
                     ENDIF
                  END DO
                  ev(:) = -MATMUL(a, b(nb+1, :))

                  scf_env%outer_scf%variables(:, ihistory+1) = 0.0_dp
                  DO i = 1, nb
                     ibuf = ihistory-nb+i
                     scf_env%outer_scf%variables(:, ihistory+1) = scf_env%outer_scf%variables(:, ihistory+1)+ &
                                                                  ev(i)*scf_env%outer_scf%variables(:, ibuf)
                  ENDDO
                  DEALLOCATE (a, b, ev)
               ENDIF
            CASE (outer_scf_optimizer_secant)
               CPASSERT(SIZE(scf_env%outer_scf%gradient, 2) >= 3)
               CPASSERT(SIZE(scf_env%outer_scf%gradient, 1) == 1)
               nvar = SIZE(scf_env%outer_scf%gradient, 1)
               IF (ihistory < 2) THEN
                  ! Need two history values to use secant, switch to sd
                  optimizer_type = outer_scf_optimizer_sd
                  CYCLE
               END IF
               ! secant update
               scf_env%outer_scf%variables(1, ihistory+1) = scf_env%outer_scf%variables(1, ihistory)- &
                                                            (scf_env%outer_scf%variables(1, ihistory)- &
                                                             scf_env%outer_scf%variables(1, ihistory-1))/ &
                                                            (scf_env%outer_scf%gradient(1, ihistory)- &
                                                             scf_env%outer_scf%gradient(1, ihistory-1))* &
                                                            scf_env%outer_scf%gradient(1, ihistory)
            CASE (outer_scf_optimizer_broyden)
               IF (.NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
                  ! Inverse Jacobian not yet built, switch to sd
                  optimizer_type = outer_scf_optimizer_sd
                  CYCLE
               END IF
               inv_jacobian => scf_env%outer_scf%inv_jacobian
               IF (ihistory < 2) THEN
                  ! Cannot perform a Broyden update without enough SCF history on this energy evaluation
                  scf_control%outer_scf%broyden_update = .FALSE.
               END IF
               IF (scf_control%outer_scf%broyden_update) THEN
                  ! Perform a Broyden update of the inverse Jacobian J^(-1)
                  IF (SIZE(scf_env%outer_scf%gradient, 2) .LT. 3) &
                     CALL cp_abort(__LOCATION__, &
                                   "Keyword EXTRAPOLATION_ORDER in section OUTER_SCF "// &
                                   "must be greater than or equal to 3 for Broyden optimizers.")
                  nvar = SIZE(scf_env%outer_scf%gradient, 1)
                  ALLOCATE (f(nvar, 1), x(nvar, 1))
                  DO i = 1, nvar
                     f(i, 1) = scf_env%outer_scf%gradient(i, ihistory)-scf_env%outer_scf%gradient(i, ihistory-1)
                     x(i, 1) = scf_env%outer_scf%variables(i, ihistory)-scf_env%outer_scf%variables(i, ihistory-1)
                  END DO
                  SELECT CASE (scf_control%outer_scf%broyden_type)
                  CASE (broyden_type_1, broyden_type_1_explicit, broyden_type_1_ls, broyden_type_1_explicit_ls)
                     ! Broyden's 1st method
                     ! Denote: dx_n = \delta x_n; df_n = \delta f_n
                     ! J_(n+1)^(-1) = J_n^(-1) + (dx_n - J_n^(-1)*df_n)*(dx_n^T * J_n^(-1))/(dx_n^T * J_n^(-1) * df_n)
                     scale = SUM(MATMUL(TRANSPOSE(x), MATMUL(inv_jacobian, f)))
                     scale = 1.0_dp/scale
                     IF (scale < 1.0E-12_dp) scale = 1.0E-12_dp
                     inv_jacobian = inv_jacobian+scale*MATMUL((x-MATMUL(inv_jacobian, f)), &
                                                              MATMUL(TRANSPOSE(x), inv_jacobian))
                  CASE (broyden_type_2, broyden_type_2_explicit, broyden_type_2_ls, broyden_type_2_explicit_ls)
                     ! Broyden's 2nd method
                     ! J_(n+1)^(-1) = J_n^(-1) + (dx_n - J_n^(-1)*df_n)*(df_n^T)/(||df_n||^2)
                     scale = SUM(MATMUL(TRANSPOSE(f), f))
                     scale = 1.0_dp/scale
                     IF (scale < 1.0E-12_dp) scale = 1.0E-12_dp
                     inv_jacobian = inv_jacobian+scale*MATMUL((x-MATMUL(inv_jacobian, f)), TRANSPOSE(inv_jacobian))
                  CASE DEFAULT
                     CALL cp_abort(__LOCATION__, &
                                   "Unknown Broyden type: "// &
                                   cp_to_string(scf_control%outer_scf%broyden_type))
                  END SELECT
                  ! Clean up
                  DEALLOCATE (f, x)
               END IF
               ! Update variables x_(n+1) = x_n - J^(-1)*f(x_n)
               scf_env%outer_scf%variables(:, ihistory+1) = scf_env%outer_scf%variables(:, ihistory)- &
                                                            scf_control%outer_scf%newton_step* &
                                                            MATMUL(inv_jacobian, scf_env%outer_scf%gradient(:, ihistory))
               scf_control%outer_scf%broyden_update = .TRUE.
            CASE (outer_scf_optimizer_newton, outer_scf_optimizer_newton_ls)
               CPASSERT(ASSOCIATED(scf_env%outer_scf%inv_jacobian))
               inv_jacobian => scf_env%outer_scf%inv_jacobian
               scf_env%outer_scf%variables(:, ihistory+1) = scf_env%outer_scf%variables(:, ihistory)- &
                                                            scf_control%outer_scf%newton_step* &
                                                            MATMUL(inv_jacobian, scf_env%outer_scf%gradient(:, ihistory))
            CASE DEFAULT
               CPABORT("")
            END SELECT
            EXIT
         ENDDO
      END IF

      CALL timestop(handle)

   END SUBROUTINE outer_loop_optimize

! **************************************************************************************************
!> \brief propagates the updated variables to wherever they need to be set in
!>       qs_env
!> \param qs_env ...
!> \param scf_env ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE outer_loop_update_qs_env(qs_env, scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_update_qs_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ihistory, n
      LOGICAL                                            :: is_constraint
      TYPE(becke_restraint_type), POINTER                :: becke_control
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(ddapc_restraint_type), POINTER                :: ddapc_restraint_control
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(s2_restraint_type), POINTER                   :: s2_restraint_control
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, scf_control=scf_control, dft_control=dft_control)
      ihistory = scf_env%outer_scf%iter_count

      SELECT CASE (scf_control%outer_scf%type)
      CASE (outer_scf_none)
         ! do nothing
      CASE (outer_scf_ddapc_constraint)
         DO n = 1, SIZE(dft_control%qs_control%ddapc_restraint_control)
            NULLIFY (ddapc_restraint_control)
            ddapc_restraint_control => dft_control%qs_control%ddapc_restraint_control(n)%ddapc_restraint_control
            is_constraint = (ddapc_restraint_control%functional_form == do_ddapc_constraint)
            IF (is_constraint) EXIT
         END DO
         ddapc_restraint_control%strength = scf_env%outer_scf%variables(1, ihistory+1)
      CASE (outer_scf_s2_constraint)
         s2_restraint_control => dft_control%qs_control%s2_restraint_control
         s2_restraint_control%strength = scf_env%outer_scf%variables(1, ihistory+1)
      CASE (outer_scf_becke_constraint)
         becke_control => dft_control%qs_control%becke_control
         becke_control%strength(:) = scf_env%outer_scf%variables(:, ihistory+1)
         IF (dft_control%qs_control%cdft) &
            dft_control%qs_control%cdft_control%strength(:) = becke_control%strength(:)
      CASE (outer_scf_cdft_constraint)
         cdft_control => dft_control%qs_control%cdft_control
         cdft_control%strength(:) = scf_env%outer_scf%variables(:, ihistory+1)
      CASE (outer_scf_basis_center_opt)
         CALL qs_update_basis_center_pos(qs_env)
      CASE DEFAULT
         CPABORT("")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE outer_loop_update_qs_env

! **************************************************************************************************
!> \brief uses the outer_scf_history to extrapolate new values for the variables
!>       and updates their value in qs_env accordingly
!> \param qs_env the qs_environment_type where to update the variables
!> \par History
!>      03.2006 created [Joost VandeVondele]
!> \note
!>       it assumes that the current value of qs_env still needs to be added to the history
!>       simple multilinear extrapolation is employed
! **************************************************************************************************
   SUBROUTINE outer_loop_extrapolate(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_extrapolate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ihis, ivec, n, nhistory, &
                                                            nvariables, nvec, outer_scf_ihistory
      LOGICAL                                            :: is_constraint
      REAL(kind=dp)                                      :: alpha
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: extrapolation
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: outer_scf_history
      TYPE(ddapc_restraint_type), POINTER                :: ddapc_restraint_control
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, outer_scf_history=outer_scf_history, &
                      outer_scf_ihistory=outer_scf_ihistory, &
                      scf_control=scf_control, dft_control=dft_control)

      nvariables = SIZE(outer_scf_history, 1)
      nhistory = SIZE(outer_scf_history, 2)
      ALLOCATE (extrapolation(nvariables))
      CPASSERT(nhistory > 0)

      ! add the current version of qs_env to the history
      outer_scf_ihistory = outer_scf_ihistory+1
      ivec = 1+MODULO(outer_scf_ihistory-1, nhistory)
      SELECT CASE (scf_control%outer_scf%type)
      CASE (outer_scf_none)
         outer_scf_history(1, ivec) = 0.0_dp
      CASE (outer_scf_ddapc_constraint)
         DO n = 1, SIZE(dft_control%qs_control%ddapc_restraint_control)
            NULLIFY (ddapc_restraint_control)
            ddapc_restraint_control => dft_control%qs_control%ddapc_restraint_control(n)%ddapc_restraint_control
            is_constraint = (ddapc_restraint_control%functional_form == do_ddapc_constraint)
            IF (is_constraint) EXIT
         END DO
         outer_scf_history(1, ivec) = &
            ddapc_restraint_control%strength
      CASE (outer_scf_s2_constraint)
         outer_scf_history(1, ivec) = &
            dft_control%qs_control%s2_restraint_control%strength
      CASE (outer_scf_becke_constraint)
         outer_scf_history(:, ivec) = &
            dft_control%qs_control%becke_control%strength(:)
      CASE (outer_scf_cdft_constraint)
         outer_scf_history(:, ivec) = &
            dft_control%qs_control%cdft_control%strength(:)
      CASE (outer_scf_basis_center_opt)
         outer_scf_history(1, ivec) = 0.0_dp
      CASE DEFAULT
         CPABORT("")
      END SELECT
      CALL set_qs_env(qs_env, outer_scf_ihistory=outer_scf_ihistory)
      ! multilinear extrapolation
      nvec = MIN(nhistory, outer_scf_ihistory)
      alpha = nvec
      ivec = 1+MODULO(outer_scf_ihistory-1, nhistory)
      extrapolation(:) = alpha*outer_scf_history(:, ivec)
      DO ihis = 2, nvec
         alpha = -1.0_dp*alpha*REAL(nvec-ihis+1, dp)/REAL(ihis, dp)
         ivec = 1+MODULO(outer_scf_ihistory-ihis, nhistory)
         extrapolation(:) = extrapolation+alpha*outer_scf_history(:, ivec)
      ENDDO

      ! update qs_env to use this extrapolation
      SELECT CASE (scf_control%outer_scf%type)
      CASE (outer_scf_none)
         ! nothing
      CASE (outer_scf_ddapc_constraint)
         ddapc_restraint_control%strength = extrapolation(1)
      CASE (outer_scf_s2_constraint)
         dft_control%qs_control%s2_restraint_control%strength = extrapolation(1)
      CASE (outer_scf_becke_constraint)
         dft_control%qs_control%becke_control%strength(:) = extrapolation(:)
         IF (dft_control%qs_control%cdft) &
            dft_control%qs_control%cdft_control%strength(:) = dft_control%qs_control%becke_control%strength(:)
      CASE (outer_scf_cdft_constraint)
         dft_control%qs_control%cdft_control%strength(:) = extrapolation(:)
      CASE (outer_scf_basis_center_opt)
         ! nothing to do
      CASE DEFAULT
         CPABORT("")
      END SELECT

      DEALLOCATE (extrapolation)

      CALL timestop(handle)

   END SUBROUTINE outer_loop_extrapolate

! **************************************************************************************************
!> \brief switch between two outer_scf envs stored in cdft_control
!> \param scf_env the scf_env where values need to be updated using cdft_control
!> \param scf_control the scf_control where values need to be updated using cdft_control
!> \param cdft_control container for the second outer_scf env
!> \param dir determines what switching operation to perform
!> \par History
!>      12.2015 created [Nico Holmberg]
! **************************************************************************************************

   SUBROUTINE outer_loop_switch(scf_env, scf_control, cdft_control, dir)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      INTEGER, INTENT(IN)                                :: dir

      CHARACTER(len=*), PARAMETER :: routineN = 'outer_loop_switch', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: nvariables

      SELECT CASE (dir)
      CASE (cdft2ot)
         ! Constraint -> OT
         ! Switch data in scf_control: first save values that might have changed
         cdft_control%constraint_control%broyden_update = scf_control%outer_scf%broyden_update
         cdft_control%constraint_control%newton_step = scf_control%outer_scf%newton_step
         cdft_control%constraint_control%ijacobian = scf_control%outer_scf%ijacobian
         ! Now switch
         scf_control%outer_scf%have_scf = cdft_control%ot_control%have_scf
         scf_control%outer_scf%max_scf = cdft_control%ot_control%max_scf
         scf_control%outer_scf%eps_scf = cdft_control%ot_control%eps_scf
         scf_control%outer_scf%step_size = cdft_control%ot_control%step_size
         scf_control%outer_scf%type = cdft_control%ot_control%type
         scf_control%outer_scf%optimizer = cdft_control%ot_control%optimizer
         scf_control%outer_scf%diis_buffer_length = cdft_control%ot_control%diis_buffer_length
         scf_control%outer_scf%bisect_trust_count = cdft_control%ot_control%bisect_trust_count
         scf_control%outer_scf%jacobian_type = cdft_control%ot_control%jacobian_type
         scf_control%outer_scf%jacobian_step = cdft_control%ot_control%jacobian_step
         scf_control%outer_scf%jacobian_freq = cdft_control%ot_control%jacobian_freq
         scf_control%outer_scf%newton_step = cdft_control%ot_control%newton_step
         scf_control%outer_scf%ijacobian = cdft_control%ot_control%ijacobian
         scf_control%outer_scf%build_jacobian = cdft_control%ot_control%build_jacobian
         scf_control%outer_scf%broyden_type = cdft_control%ot_control%broyden_type
         scf_control%outer_scf%broyden_update = cdft_control%ot_control%broyden_update
         ! Switch data in scf_env: first save current values for constraint
         cdft_control%constraint%iter_count = scf_env%outer_scf%iter_count
         cdft_control%constraint%energy = scf_env%outer_scf%energy
         cdft_control%constraint%variables = scf_env%outer_scf%variables
         cdft_control%constraint%gradient = scf_env%outer_scf%gradient
         cdft_control%constraint%count = scf_env%outer_scf%count
         cdft_control%constraint%deallocate_jacobian = scf_env%outer_scf%deallocate_jacobian
         IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
            nvariables = SIZE(scf_env%outer_scf%inv_jacobian, 1)
            IF (.NOT. ASSOCIATED(cdft_control%constraint%inv_jacobian)) &
               ALLOCATE (cdft_control%constraint%inv_jacobian(nvariables, nvariables))
            cdft_control%constraint%inv_jacobian = scf_env%outer_scf%inv_jacobian
         END IF
         ! Now switch
         IF (ASSOCIATED(scf_env%outer_scf%energy)) &
            DEALLOCATE (scf_env%outer_scf%energy)
         ALLOCATE (scf_env%outer_scf%energy(scf_control%outer_scf%max_scf+1))
         scf_env%outer_scf%energy = 0.0_dp
         IF (ASSOCIATED(scf_env%outer_scf%variables)) &
            DEALLOCATE (scf_env%outer_scf%variables)
         ALLOCATE (scf_env%outer_scf%variables(1, scf_control%outer_scf%max_scf+1))
         scf_env%outer_scf%variables = 0.0_dp
         IF (ASSOCIATED(scf_env%outer_scf%gradient)) &
            DEALLOCATE (scf_env%outer_scf%gradient)
         ALLOCATE (scf_env%outer_scf%gradient(1, scf_control%outer_scf%max_scf+1))
         scf_env%outer_scf%gradient = 0.0_dp
         IF (ASSOCIATED(scf_env%outer_scf%count)) &
            DEALLOCATE (scf_env%outer_scf%count)
         ALLOCATE (scf_env%outer_scf%count(scf_control%outer_scf%max_scf+1))
         scf_env%outer_scf%count = 0
         ! OT SCF does not need Jacobian
         scf_env%outer_scf%deallocate_jacobian = .TRUE.
         IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) &
            DEALLOCATE (scf_env%outer_scf%inv_jacobian)
      CASE (ot2cdft)
         ! OT -> constraint
         scf_control%outer_scf%have_scf = cdft_control%constraint_control%have_scf
         scf_control%outer_scf%max_scf = cdft_control%constraint_control%max_scf
         scf_control%outer_scf%eps_scf = cdft_control%constraint_control%eps_scf
         scf_control%outer_scf%step_size = cdft_control%constraint_control%step_size
         scf_control%outer_scf%type = cdft_control%constraint_control%type
         scf_control%outer_scf%optimizer = cdft_control%constraint_control%optimizer
         scf_control%outer_scf%diis_buffer_length = cdft_control%constraint_control%diis_buffer_length
         scf_control%outer_scf%bisect_trust_count = cdft_control%constraint_control%bisect_trust_count
         scf_control%outer_scf%jacobian_type = cdft_control%constraint_control%jacobian_type
         scf_control%outer_scf%jacobian_step = cdft_control%constraint_control%jacobian_step
         scf_control%outer_scf%jacobian_freq = cdft_control%constraint_control%jacobian_freq
         scf_control%outer_scf%ijacobian = cdft_control%constraint_control%ijacobian
         scf_control%outer_scf%build_jacobian = cdft_control%constraint_control%build_jacobian
         scf_control%outer_scf%broyden_type = cdft_control%constraint_control%broyden_type
         scf_control%outer_scf%broyden_update = cdft_control%constraint_control%broyden_update
         scf_control%outer_scf%newton_step = cdft_control%constraint_control%newton_step
         nvariables = SIZE(cdft_control%constraint%variables, 1)
         IF (ASSOCIATED(scf_env%outer_scf%energy)) &
            DEALLOCATE (scf_env%outer_scf%energy)
         ALLOCATE (scf_env%outer_scf%energy(scf_control%outer_scf%max_scf+1))
         scf_env%outer_scf%energy = cdft_control%constraint%energy
         IF (ASSOCIATED(scf_env%outer_scf%variables)) &
            DEALLOCATE (scf_env%outer_scf%variables)
         ALLOCATE (scf_env%outer_scf%variables(nvariables, scf_control%outer_scf%max_scf+1))
         scf_env%outer_scf%variables = cdft_control%constraint%variables
         IF (ASSOCIATED(scf_env%outer_scf%gradient)) &
            DEALLOCATE (scf_env%outer_scf%gradient)
         ALLOCATE (scf_env%outer_scf%gradient(nvariables, scf_control%outer_scf%max_scf+1))
         scf_env%outer_scf%gradient = cdft_control%constraint%gradient
         IF (ASSOCIATED(scf_env%outer_scf%count)) &
            DEALLOCATE (scf_env%outer_scf%count)
         ALLOCATE (scf_env%outer_scf%count(scf_control%outer_scf%max_scf+1))
         scf_env%outer_scf%count = cdft_control%constraint%count
         scf_env%outer_scf%iter_count = cdft_control%constraint%iter_count
         scf_env%outer_scf%deallocate_jacobian = cdft_control%constraint%deallocate_jacobian
         IF (ASSOCIATED(cdft_control%constraint%inv_jacobian)) THEN
            IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) &
               DEALLOCATE (scf_env%outer_scf%inv_jacobian)
            ALLOCATE (scf_env%outer_scf%inv_jacobian(nvariables, nvariables))
            scf_env%outer_scf%inv_jacobian = cdft_control%constraint%inv_jacobian
         END IF
      CASE DEFAULT
         CPABORT("")
      END SELECT

   END SUBROUTINE outer_loop_switch

! **************************************************************************************************
!> \brief purges outer_scf_history zeroing everything except
!>        the latest value of the outer_scf variable stored in qs_control
!> \param qs_env the qs_environment_type where to purge
!> \par History
!>      05.2016 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE outer_loop_purge_history(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_purge_history', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, outer_scf_ihistory
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: gradient_history, outer_scf_history, &
                                                            variable_history

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, outer_scf_history=outer_scf_history, &
                      outer_scf_ihistory=outer_scf_ihistory, &
                      gradient_history=gradient_history, &
                      variable_history=variable_history)
      CPASSERT(SIZE(outer_scf_history, 2) > 0)
      outer_scf_ihistory = 0
      outer_scf_history = 0.0_dp
      gradient_history = 0.0_dp
      variable_history = 0.0_dp
      CALL set_qs_env(qs_env, outer_scf_ihistory=outer_scf_ihistory)

      CALL timestop(handle)

   END SUBROUTINE outer_loop_purge_history

END MODULE qs_outer_scf
