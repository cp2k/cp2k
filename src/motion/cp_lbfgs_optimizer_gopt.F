!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines that optimize a functional using the limited memory bfgs
!>      quasi-newton method.
!>      The process set up so that a master runs the real optimizer and the
!>      others help then to calculate the objective function.
!>      The arguments for the objective function are physically present in
!>      every processor (nedeed in the actual implementation of pao).
!>      In the future tha arguments themselves could be distributed.
!> \par History
!>      09.2003 globenv->para_env, retain/release, better parallel behaviour
!>      01.2020 Space Group Symmetry introduced by Pierre-André Cazade [pcazade]
!> \author Fawzi Mohamed
!>      @version 2.2002
! **************************************************************************************************
MODULE cp_lbfgs_optimizer_gopt
   USE cp_lbfgs, ONLY: setulb
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_type, &
                              cp_to_string
   USE cp_output_handling, ONLY: cp_print_key_finished_output, &
                                 cp_print_key_unit_nr
   USE message_passing, ONLY: mp_para_env_release
   USE message_passing, ONLY: mp_para_env_type
   USE cp_subsys_types, ONLY: cp_subsys_type
   USE force_env_types, ONLY: force_env_get, &
                              force_env_type
   USE gopt_f_methods, ONLY: gopt_f_io
   USE gopt_f_types, ONLY: gopt_f_release, &
                           gopt_f_retain, &
                           gopt_f_type
   USE gopt_param_types, ONLY: gopt_param_type
   USE input_section_types, ONLY: section_vals_type
   USE kinds, ONLY: dp
   USE machine, ONLY: m_walltime
   USE space_groups, ONLY: spgr_apply_rotations_coord, &
                           spgr_apply_rotations_force
   USE space_groups_types, ONLY: spgr_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   #:include "gopt_f77_methods.fypp"

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_lbfgs_optimizer_gopt'

   ! types
   PUBLIC :: cp_lbfgs_opt_gopt_type

   ! core methods

   ! special methos

   ! underlying functions
   PUBLIC :: cp_opt_gopt_create, cp_opt_gopt_release, &
             cp_opt_gopt_next, &
             cp_opt_gopt_stop

! **************************************************************************************************
!> \brief info for the optimizer (see the description of this module)
!> \param task the actual task of the optimizer (in the master it is up to
!>        date, in case of error also the minions one get updated.
!> \param csave internal character string used by the lbfgs optimizer,
!>        meaningful only in the master
!> \param lsave logical array used by the lbfgs optimizer, updated only
!>        in the master
!>        On exit with task = 'NEW_X', the following information is
!>        available:
!>           lsave(1) = .true.  the initial x did not satisfy the bounds;
!>           lsave(2) = .true.  the problem contains bounds;
!>           lsave(3) = .true.  each variable has upper and lower bounds.
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param m the dimension of the subspace used to approximate the second
!>        derivative
!> \param print_every every how many iterations output should be written.
!>        if 0 only at end, if print_every<0 never
!> \param master the pid of the master processor
!> \param max_f_per_iter the maximum number of function evaluations per
!>        iteration
!> \param status 0: just initialized, 1: f g calculation,
!>        2: begin new iteration, 3: ended iteration,
!>        4: normal (converged) exit, 5: abnormal (error) exit,
!>        6: daellocated
!> \param n_iter the actual iteration number
!> \param kind_of_bound an array with 0 (no bound), 1 (lower bound),
!>        2 (both bounds), 3 (upper bound), to describe the bounds
!>        of every variable
!> \param i_work_array an integer workarray of dimension 3*n, present only
!>        in the master
!> \param isave is an INTEGER working array of dimension 44.
!>        On exit with task = 'NEW_X', it contains information that
!>        the user may want to access:
!> \param isave (30) = the current iteration number;
!> \param isave (34) = the total number of function and gradient
!>           evaluations;
!> \param isave (36) = the number of function value or gradient
!>           evaluations in the current iteration;
!> \param isave (38) = the number of free variables in the current
!>           iteration;
!> \param isave (39) = the number of active constraints at the current
!>           iteration;
!> \param f the actual best value of the object function
!> \param wanted_relative_f_delta the wanted relative error on f
!>        (to be multiplied by epsilon), 0.0 -> no check
!> \param wanted_projected_gradient the wanted error on the projected
!>        gradient (hessian times the gradient), 0.0 -> no check
!> \param last_f the value of f in the last iteration
!> \param projected_gradient the value of the sup norm of the projected
!>        gradient
!> \param x the actual evaluation point (best one if converged or stopped)
!> \param lower_bound the lower bounds
!> \param upper_bound the upper bounds
!> \param gradient the actual gradient
!> \param dsave info date for lbfgs (master only)
!> \param work_array a work array for lbfgs (master only)
!> \param para_env the parallel environment for this optimizer
!> \param obj_funct the objective function to be optimized
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! **************************************************************************************************
   TYPE cp_lbfgs_opt_gopt_type
      CHARACTER(len=60) :: task = ""
      CHARACTER(len=60) :: csave = ""
      LOGICAL :: lsave(4) = .FALSE.
      INTEGER :: m = 0, print_every = 0, master = 0, max_f_per_iter = 0, status = 0, n_iter = 0
      INTEGER, DIMENSION(:), POINTER :: kind_of_bound => NULL(), i_work_array => NULL(), isave => NULL()
      REAL(kind=dp) :: f = 0.0_dp, wanted_relative_f_delta = 0.0_dp, wanted_projected_gradient = 0.0_dp, &
                       last_f = 0.0_dp, projected_gradient = 0.0_dp, eold = 0.0_dp, emin = 0.0_dp, trust_radius = 0.0_dp
      REAL(kind=dp), DIMENSION(:), POINTER :: x => NULL(), lower_bound => NULL(), upper_bound => NULL(), &
                                              gradient => NULL(), dsave => NULL(), work_array => NULL()
      TYPE(mp_para_env_type), POINTER :: para_env => NULL()
      TYPE(gopt_f_type), POINTER :: obj_funct => NULL()
   END TYPE cp_lbfgs_opt_gopt_type

CONTAINS

! **************************************************************************************************
!> \brief initializes the optimizer
!> \param optimizer ...
!> \param para_env ...
!> \param obj_funct ...
!> \param x0 ...
!> \param m ...
!> \param print_every ...
!> \param wanted_relative_f_delta ...
!> \param wanted_projected_gradient ...
!> \param lower_bound ...
!> \param upper_bound ...
!> \param kind_of_bound ...
!> \param master ...
!> \param max_f_per_iter ...
!> \param trust_radius ...
!> \par History
!>      02.2002 created [fawzi]
!>      09.2003 refactored (retain/release,para_env) [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      redirects the lbfgs output the the default unit
! **************************************************************************************************
   SUBROUTINE cp_opt_gopt_create(optimizer, para_env, obj_funct, x0, m, print_every, &
                                 wanted_relative_f_delta, wanted_projected_gradient, lower_bound, upper_bound, &
                                 kind_of_bound, master, max_f_per_iter, trust_radius)
      TYPE(cp_lbfgs_opt_gopt_type), INTENT(OUT)          :: optimizer
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(gopt_f_type), POINTER                         :: obj_funct
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: x0
      INTEGER, INTENT(in), OPTIONAL                      :: m, print_every
      REAL(kind=dp), INTENT(in), OPTIONAL                :: wanted_relative_f_delta, &
                                                            wanted_projected_gradient
      REAL(kind=dp), DIMENSION(SIZE(x0)), INTENT(in), &
         OPTIONAL                                        :: lower_bound, upper_bound
      INTEGER, DIMENSION(SIZE(x0)), INTENT(in), OPTIONAL :: kind_of_bound
      INTEGER, INTENT(in), OPTIONAL                      :: master, max_f_per_iter
      REAL(kind=dp), INTENT(in), OPTIONAL                :: trust_radius

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_opt_gopt_create'

      INTEGER                                            :: handle, lenwa, n

      CALL timeset(routineN, handle)

      NULLIFY (optimizer%kind_of_bound, &
               optimizer%i_work_array, &
               optimizer%isave, &
               optimizer%x, &
               optimizer%lower_bound, &
               optimizer%upper_bound, &
               optimizer%gradient, &
               optimizer%dsave, &
               optimizer%work_array, &
               optimizer%para_env, &
               optimizer%obj_funct)
      n = SIZE(x0)
      optimizer%m = 4
      IF (PRESENT(m)) optimizer%m = m
      optimizer%master = para_env%source
      optimizer%para_env => para_env
      CALL para_env%retain()
      optimizer%obj_funct => obj_funct
      CALL gopt_f_retain(obj_funct)
      optimizer%max_f_per_iter = 20
      IF (PRESENT(max_f_per_iter)) optimizer%max_f_per_iter = max_f_per_iter
      optimizer%print_every = -1
      optimizer%n_iter = 0
      optimizer%f = -1.0_dp
      optimizer%last_f = -1.0_dp
      optimizer%projected_gradient = -1.0_dp
      IF (PRESENT(print_every)) optimizer%print_every = print_every
      IF (PRESENT(master)) optimizer%master = master
      IF (optimizer%master == optimizer%para_env%mepos) THEN
         !MK This has to be adapted for a new L-BFGS version possibly
         lenwa = 2*optimizer%m*n + 5*n + 11*optimizer%m*optimizer%m + 8*optimizer%m
         ALLOCATE (optimizer%kind_of_bound(n), optimizer%i_work_array(3*n), &
                   optimizer%isave(44))
         ALLOCATE (optimizer%x(n), optimizer%lower_bound(n), &
                   optimizer%upper_bound(n), optimizer%gradient(n), &
                   optimizer%dsave(29), optimizer%work_array(lenwa))
         optimizer%x = x0
         optimizer%task = 'START'
         optimizer%wanted_relative_f_delta = wanted_relative_f_delta
         optimizer%wanted_projected_gradient = wanted_projected_gradient
         optimizer%kind_of_bound = 0
         IF (PRESENT(kind_of_bound)) optimizer%kind_of_bound = kind_of_bound
         IF (PRESENT(lower_bound)) optimizer%lower_bound = lower_bound
         IF (PRESENT(upper_bound)) optimizer%upper_bound = upper_bound
         IF (PRESENT(trust_radius)) optimizer%trust_radius = trust_radius

         CALL setulb(SIZE(optimizer%x), optimizer%m, optimizer%x, &
                     optimizer%lower_bound, optimizer%upper_bound, &
                     optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
                     optimizer%wanted_relative_f_delta, &
                     optimizer%wanted_projected_gradient, optimizer%work_array, &
                     optimizer%i_work_array, optimizer%task, optimizer%print_every, &
                     optimizer%csave, optimizer%lsave, optimizer%isave, &
                     optimizer%dsave, optimizer%trust_radius)
      ELSE
         NULLIFY ( &
            optimizer%kind_of_bound, optimizer%i_work_array, optimizer%isave, &
            optimizer%lower_bound, optimizer%upper_bound, optimizer%gradient, &
            optimizer%dsave, optimizer%work_array)
         ALLOCATE (optimizer%x(n))
         optimizer%x(:) = 0.0_dp
         ALLOCATE (optimizer%gradient(n))
         optimizer%gradient(:) = 0.0_dp
      END IF
      CALL optimizer%para_env%bcast(optimizer%x, optimizer%master)
      optimizer%status = 0

      CALL timestop(handle)

   END SUBROUTINE cp_opt_gopt_create

! **************************************************************************************************
!> \brief releases the optimizer (see doc/ReferenceCounting.html)
!> \param optimizer the object that should be released
!> \par History
!>      02.2002 created [fawzi]
!>      09.2003 dealloc_ref->release [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cp_opt_gopt_release(optimizer)
      TYPE(cp_lbfgs_opt_gopt_type), INTENT(INOUT)        :: optimizer

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_opt_gopt_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(optimizer%kind_of_bound)) THEN
         DEALLOCATE (optimizer%kind_of_bound)
      END IF
      IF (ASSOCIATED(optimizer%i_work_array)) THEN
         DEALLOCATE (optimizer%i_work_array)
      END IF
      IF (ASSOCIATED(optimizer%isave)) THEN
         DEALLOCATE (optimizer%isave)
      END IF
      IF (ASSOCIATED(optimizer%x)) THEN
         DEALLOCATE (optimizer%x)
      END IF
      IF (ASSOCIATED(optimizer%lower_bound)) THEN
         DEALLOCATE (optimizer%lower_bound)
      END IF
      IF (ASSOCIATED(optimizer%upper_bound)) THEN
         DEALLOCATE (optimizer%upper_bound)
      END IF
      IF (ASSOCIATED(optimizer%gradient)) THEN
         DEALLOCATE (optimizer%gradient)
      END IF
      IF (ASSOCIATED(optimizer%dsave)) THEN
         DEALLOCATE (optimizer%dsave)
      END IF
      IF (ASSOCIATED(optimizer%work_array)) THEN
         DEALLOCATE (optimizer%work_array)
      END IF
      CALL mp_para_env_release(optimizer%para_env)
      CALL gopt_f_release(optimizer%obj_funct)

      CALL timestop(handle)
   END SUBROUTINE cp_opt_gopt_release

! **************************************************************************************************
!> \brief takes different valuse from the optimizer
!> \param optimizer ...
!> \param para_env ...
!> \param obj_funct ...
!> \param m ...
!> \param print_every ...
!> \param wanted_relative_f_delta ...
!> \param wanted_projected_gradient ...
!> \param x ...
!> \param lower_bound ...
!> \param upper_bound ...
!> \param kind_of_bound ...
!> \param master ...
!> \param actual_projected_gradient ...
!> \param n_var ...
!> \param n_iter ...
!> \param status ...
!> \param max_f_per_iter ...
!> \param at_end ...
!> \param is_master ...
!> \param last_f ...
!> \param f ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! **************************************************************************************************
   SUBROUTINE cp_opt_gopt_get(optimizer, para_env, &
                              obj_funct, m, print_every, &
                              wanted_relative_f_delta, wanted_projected_gradient, &
                              x, lower_bound, upper_bound, kind_of_bound, master, &
                              actual_projected_gradient, &
                              n_var, n_iter, status, max_f_per_iter, at_end, &
                              is_master, last_f, f)
      TYPE(cp_lbfgs_opt_gopt_type), INTENT(IN)           :: optimizer
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env
      TYPE(gopt_f_type), OPTIONAL, POINTER               :: obj_funct
      INTEGER, INTENT(out), OPTIONAL                     :: m, print_every
      REAL(kind=dp), INTENT(out), OPTIONAL               :: wanted_relative_f_delta, &
                                                            wanted_projected_gradient
      REAL(kind=dp), DIMENSION(:), OPTIONAL, POINTER     :: x, lower_bound, upper_bound
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: kind_of_bound
      INTEGER, INTENT(out), OPTIONAL                     :: master
      REAL(kind=dp), INTENT(out), OPTIONAL               :: actual_projected_gradient
      INTEGER, INTENT(out), OPTIONAL                     :: n_var, n_iter, status, max_f_per_iter
      LOGICAL, INTENT(out), OPTIONAL                     :: at_end, is_master
      REAL(kind=dp), INTENT(out), OPTIONAL               :: last_f, f

      IF (PRESENT(is_master)) is_master = optimizer%master == optimizer%para_env%mepos
      IF (PRESENT(master)) master = optimizer%master
      IF (PRESENT(status)) status = optimizer%status
      IF (PRESENT(para_env)) para_env => optimizer%para_env
      IF (PRESENT(obj_funct)) obj_funct = optimizer%obj_funct
      IF (PRESENT(m)) m = optimizer%m
      IF (PRESENT(max_f_per_iter)) max_f_per_iter = optimizer%max_f_per_iter
      IF (PRESENT(wanted_projected_gradient)) &
         wanted_projected_gradient = optimizer%wanted_projected_gradient
      IF (PRESENT(wanted_relative_f_delta)) &
         wanted_relative_f_delta = optimizer%wanted_relative_f_delta
      IF (PRESENT(print_every)) print_every = optimizer%print_every
      IF (PRESENT(x)) x => optimizer%x
      IF (PRESENT(n_var)) n_var = SIZE(x)
      IF (PRESENT(lower_bound)) lower_bound => optimizer%lower_bound
      IF (PRESENT(upper_bound)) upper_bound => optimizer%upper_bound
      IF (PRESENT(kind_of_bound)) kind_of_bound => optimizer%kind_of_bound
      IF (PRESENT(n_iter)) n_iter = optimizer%n_iter
      IF (PRESENT(last_f)) last_f = optimizer%last_f
      IF (PRESENT(f)) f = optimizer%f
      IF (PRESENT(at_end)) at_end = optimizer%status > 3
      IF (PRESENT(actual_projected_gradient)) &
         actual_projected_gradient = optimizer%projected_gradient
      IF (optimizer%master == optimizer%para_env%mepos) THEN
         IF (optimizer%isave(30) > 1 .AND. (optimizer%task(1:5) == "NEW_X" .OR. &
                                            optimizer%task(1:4) == "STOP" .AND. optimizer%task(7:9) == "CPU")) THEN
            ! nr iterations >1 .and. dsave contains the wanted data
            IF (PRESENT(last_f)) last_f = optimizer%dsave(2)
            IF (PRESENT(actual_projected_gradient)) &
               actual_projected_gradient = optimizer%dsave(13)
         ELSE
            CPASSERT(.NOT. PRESENT(last_f))
            CPASSERT(.NOT. PRESENT(actual_projected_gradient))
         END IF
      ELSE IF (PRESENT(lower_bound) .OR. PRESENT(upper_bound) .OR. PRESENT(kind_of_bound)) THEN
         CPWARN("asked undefined types")
      END IF

   END SUBROUTINE cp_opt_gopt_get

! **************************************************************************************************
!> \brief does one optimization step
!> \param optimizer ...
!> \param n_iter ...
!> \param f ...
!> \param last_f ...
!> \param projected_gradient ...
!> \param converged ...
!> \param geo_section ...
!> \param force_env ...
!> \param gopt_param ...
!> \param spgr ...
!> \par History
!>      01.2020 modified [pcazade]
!> \author Fawzi Mohamed
!>      @version 2.2002
!> \note
!>      use directly mainlb in place of setulb ??
! **************************************************************************************************
   SUBROUTINE cp_opt_gopt_step(optimizer, n_iter, f, last_f, &
                               projected_gradient, converged, geo_section, force_env, &
                               gopt_param, spgr)
      TYPE(cp_lbfgs_opt_gopt_type), INTENT(INOUT)        :: optimizer
      INTEGER, INTENT(out), OPTIONAL                     :: n_iter
      REAL(kind=dp), INTENT(out), OPTIONAL               :: f, last_f, projected_gradient
      LOGICAL, INTENT(out), OPTIONAL                     :: converged
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(spgr_type), OPTIONAL, POINTER                 :: spgr

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_opt_gopt_step'

      CHARACTER(LEN=5)                                   :: wildcard
      INTEGER                                            :: dataunit, handle, its
      LOGICAL                                            :: conv, is_master, justEntred, &
                                                            keep_space_group
      REAL(KIND=dp)                                      :: t_diff, t_now, t_old
      REAL(KIND=dp), DIMENSION(:), POINTER               :: xold
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys

      NULLIFY (logger, xold)
      logger => cp_get_default_logger()
      CALL timeset(routineN, handle)
      justEntred = .TRUE.
      is_master = optimizer%master == optimizer%para_env%mepos
      IF (PRESENT(converged)) converged = optimizer%status == 4
      ALLOCATE (xold(SIZE(optimizer%x)))

      ! collecting subsys
      CALL force_env_get(force_env, subsys=subsys)

      keep_space_group = .FALSE.
      IF (PRESENT(spgr)) THEN
         IF (ASSOCIATED(spgr)) keep_space_group = spgr%keep_space_group
      END IF

      ! applies rotation matrices to coordinates
      IF (keep_space_group) THEN
         CALL spgr_apply_rotations_coord(spgr, optimizer%x)
      END IF

      xold = optimizer%x
      t_old = m_walltime()

      IF (optimizer%status >= 4) THEN
         CPWARN("status>=4, trying to restart")
         optimizer%status = 0
         IF (is_master) THEN
            optimizer%task = 'START'
            CALL setulb(SIZE(optimizer%x), optimizer%m, optimizer%x, &
                        optimizer%lower_bound, optimizer%upper_bound, &
                        optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
                        optimizer%wanted_relative_f_delta, &
                        optimizer%wanted_projected_gradient, optimizer%work_array, &
                        optimizer%i_work_array, optimizer%task, optimizer%print_every, &
                        optimizer%csave, optimizer%lsave, optimizer%isave, &
                        optimizer%dsave, optimizer%trust_radius, spgr=spgr)
         END IF
      END IF

      DO
         ifMaster: IF (is_master) THEN
            IF (optimizer%task(1:7) == 'RESTART') THEN
               ! restart the optimizer
               optimizer%status = 0
               optimizer%task = 'START'
               ! applies rotation matrices to coordinates and forces
               IF (keep_space_group) THEN
                  CALL spgr_apply_rotations_coord(spgr, optimizer%x)
                  CALL spgr_apply_rotations_force(spgr, optimizer%gradient)
               END IF
               CALL setulb(SIZE(optimizer%x), optimizer%m, optimizer%x, &
                           optimizer%lower_bound, optimizer%upper_bound, &
                           optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
                           optimizer%wanted_relative_f_delta, &
                           optimizer%wanted_projected_gradient, optimizer%work_array, &
                           optimizer%i_work_array, optimizer%task, optimizer%print_every, &
                           optimizer%csave, optimizer%lsave, optimizer%isave, &
                           optimizer%dsave, optimizer%trust_radius, spgr=spgr)
               IF (keep_space_group) THEN
                  CALL spgr_apply_rotations_coord(spgr, optimizer%x)
                  CALL spgr_apply_rotations_force(spgr, optimizer%gradient)
               END IF
            END IF
            IF (optimizer%task(1:2) == 'FG') THEN
               IF (optimizer%isave(36) > optimizer%max_f_per_iter) THEN
                  optimizer%task = 'STOP: CPU, hit max f eval in iter'
                  optimizer%status = 5 ! anormal exit
                  CALL setulb(SIZE(optimizer%x), optimizer%m, optimizer%x, &
                              optimizer%lower_bound, optimizer%upper_bound, &
                              optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
                              optimizer%wanted_relative_f_delta, &
                              optimizer%wanted_projected_gradient, optimizer%work_array, &
                              optimizer%i_work_array, optimizer%task, optimizer%print_every, &
                              optimizer%csave, optimizer%lsave, optimizer%isave, &
                              optimizer%dsave, optimizer%trust_radius, spgr=spgr)
               ELSE
                  optimizer%status = 1
               END IF
            ELSE IF (optimizer%task(1:5) == 'NEW_X') THEN
               IF (justEntred) THEN
                  optimizer%status = 2
                  ! applies rotation matrices to coordinates and forces
                  IF (keep_space_group) THEN
                     CALL spgr_apply_rotations_coord(spgr, optimizer%x)
                     CALL spgr_apply_rotations_force(spgr, optimizer%gradient)
                  END IF
                  CALL setulb(SIZE(optimizer%x), optimizer%m, optimizer%x, &
                              optimizer%lower_bound, optimizer%upper_bound, &
                              optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
                              optimizer%wanted_relative_f_delta, &
                              optimizer%wanted_projected_gradient, optimizer%work_array, &
                              optimizer%i_work_array, optimizer%task, optimizer%print_every, &
                              optimizer%csave, optimizer%lsave, optimizer%isave, &
                              optimizer%dsave, optimizer%trust_radius, spgr=spgr)
                  IF (keep_space_group) THEN
                     CALL spgr_apply_rotations_coord(spgr, optimizer%x)
                     CALL spgr_apply_rotations_force(spgr, optimizer%gradient)
                  END IF
               ELSE
                  ! applies rotation matrices to coordinates and forces
                  IF (keep_space_group) THEN
                     CALL spgr_apply_rotations_coord(spgr, optimizer%x)
                     CALL spgr_apply_rotations_force(spgr, optimizer%gradient)
                  END IF
                  optimizer%status = 3
               END IF
            ELSE IF (optimizer%task(1:4) == 'CONV') THEN
               optimizer%status = 4
            ELSE IF (optimizer%task(1:4) == 'STOP') THEN
               optimizer%status = 5
               CPWARN("task became stop in an unknown way")
            ELSE IF (optimizer%task(1:5) == 'ERROR') THEN
               optimizer%status = 5
            ELSE
               CPWARN("unknown task '"//optimizer%task//"'")
            END IF
         END IF ifMaster
         CALL optimizer%para_env%bcast(optimizer%status, optimizer%master)
         ! Dump info
         IF (optimizer%status == 3) THEN
            its = 0
            IF (is_master) THEN
               ! Iteration level is taken into account in the optimizer external loop
               its = optimizer%isave(30)
            END IF
         END IF
         !
         SELECT CASE (optimizer%status)
         CASE (1)
            !op=1 evaluate f and g
            CALL cp_eval_at(optimizer%obj_funct, x=optimizer%x, &
                            f=optimizer%f, &
                            gradient=optimizer%gradient, &
                            final_evaluation=.FALSE., &
                            master=optimizer%master, para_env=optimizer%para_env)
            ! do not use keywords?
            IF (is_master) THEN
               ! applies rotation matrices to coordinates and forces
               IF (keep_space_group) THEN
                  CALL spgr_apply_rotations_coord(spgr, optimizer%x)
                  CALL spgr_apply_rotations_force(spgr, optimizer%gradient)
               END IF
               CALL setulb(SIZE(optimizer%x), optimizer%m, optimizer%x, &
                           optimizer%lower_bound, optimizer%upper_bound, &
                           optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
                           optimizer%wanted_relative_f_delta, &
                           optimizer%wanted_projected_gradient, optimizer%work_array, &
                           optimizer%i_work_array, optimizer%task, optimizer%print_every, &
                           optimizer%csave, optimizer%lsave, optimizer%isave, &
                           optimizer%dsave, optimizer%trust_radius, spgr=spgr)
               IF (keep_space_group) THEN
                  CALL spgr_apply_rotations_coord(spgr, optimizer%x)
                  CALL spgr_apply_rotations_force(spgr, optimizer%gradient)
               END IF
            END IF
            CALL optimizer%para_env%bcast(optimizer%x, optimizer%master)
         CASE (2)
            !op=2 begin new iter
            CALL optimizer%para_env%bcast(optimizer%x, optimizer%master)
            t_old = m_walltime()
         CASE (3)
            !op=3 ended iter
            wildcard = "LBFGS"
            dataunit = cp_print_key_unit_nr(logger, geo_section, &
                                            "PRINT%PROGRAM_RUN_INFO", extension=".geoLog")
            IF (is_master) its = optimizer%isave(30)
            CALL optimizer%para_env%bcast(its, optimizer%master)

            ! Some IO and Convergence check
            t_now = m_walltime()
            t_diff = t_now - t_old
            t_old = t_now
            CALL gopt_f_io(optimizer%obj_funct, force_env, force_env%root_section, &
                           its, optimizer%f, dataunit, optimizer%eold, optimizer%emin, wildcard, gopt_param, &
                           SIZE(optimizer%x), optimizer%x - xold, optimizer%gradient, conv, used_time=t_diff)
            CALL optimizer%para_env%bcast(conv, optimizer%master)
            CALL cp_print_key_finished_output(dataunit, logger, geo_section, &
                                              "PRINT%PROGRAM_RUN_INFO")
            optimizer%eold = optimizer%f
            optimizer%emin = MIN(optimizer%emin, optimizer%eold)
            xold = optimizer%x
            IF (PRESENT(converged)) converged = conv
            EXIT
         CASE (4)
            !op=4 (convergence - normal exit)
            ! Specific L-BFGS convergence criteria.. overrides the convergence criteria on
            ! stepsize and gradients
            dataunit = cp_print_key_unit_nr(logger, geo_section, &
                                            "PRINT%PROGRAM_RUN_INFO", extension=".geoLog")
            IF (dataunit > 0) THEN
               WRITE (dataunit, '(T2,A)') ""
               WRITE (dataunit, '(T2,A)') "***********************************************"
               WRITE (dataunit, '(T2,A)') "* Specific L-BFGS convergence criteria         "
               WRITE (dataunit, '(T2,A)') "* WANTED_PROJ_GRADIENT and WANTED_REL_F_ERROR  "
               WRITE (dataunit, '(T2,A)') "* satisfied .... run CONVERGED!                "
               WRITE (dataunit, '(T2,A)') "***********************************************"
               WRITE (dataunit, '(T2,A)') ""
            END IF
            CALL cp_print_key_finished_output(dataunit, logger, geo_section, &
                                              "PRINT%PROGRAM_RUN_INFO")
            IF (PRESENT(converged)) converged = .TRUE.
            EXIT
         CASE (5)
            ! op=5 abnormal exit ()
            CALL optimizer%para_env%bcast(optimizer%task, optimizer%master)
         CASE (6)
            ! deallocated
            CPABORT("step on a deallocated opt structure ")
         CASE default
            CALL cp_abort(__LOCATION__, &
                          "unknown status "//cp_to_string(optimizer%status))
            optimizer%status = 5
            EXIT
         END SELECT
         IF (optimizer%status == 1 .AND. justEntred) THEN
            optimizer%eold = optimizer%f
            optimizer%emin = optimizer%eold
         END IF
         justEntred = .FALSE.
      END DO

      CALL optimizer%para_env%bcast(optimizer%x, optimizer%master)
      CALL cp_opt_gopt_bcast_res(optimizer, &
                                 n_iter=optimizer%n_iter, &
                                 f=optimizer%f, last_f=optimizer%last_f, &
                                 projected_gradient=optimizer%projected_gradient)

      DEALLOCATE (xold)
      IF (PRESENT(f)) f = optimizer%f
      IF (PRESENT(last_f)) last_f = optimizer%last_f
      IF (PRESENT(projected_gradient)) &
         projected_gradient = optimizer%projected_gradient
      IF (PRESENT(n_iter)) n_iter = optimizer%n_iter
      CALL timestop(handle)

   END SUBROUTINE cp_opt_gopt_step

! **************************************************************************************************
!> \brief returns the results (and broadcasts them)
!> \param optimizer the optimizer object the info is taken from
!> \param n_iter the number of iterations
!> \param f the actual value of the objective function (f)
!> \param last_f the last value of f
!> \param projected_gradient the infinity norm of the projected gradient
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
!> \note
!>      private routine
! **************************************************************************************************
   SUBROUTINE cp_opt_gopt_bcast_res(optimizer, n_iter, f, last_f, &
                                    projected_gradient)
      TYPE(cp_lbfgs_opt_gopt_type), INTENT(IN)           :: optimizer
      INTEGER, INTENT(out), OPTIONAL                     :: n_iter
      REAL(kind=dp), INTENT(inout), OPTIONAL             :: f, last_f, projected_gradient

      REAL(kind=dp), DIMENSION(4)                        :: results

      IF (optimizer%master == optimizer%para_env%mepos) THEN
         results = (/REAL(optimizer%isave(30), kind=dp), &
                     optimizer%f, optimizer%dsave(2), optimizer%dsave(13)/)
      END IF
      CALL optimizer%para_env%bcast(results, optimizer%master)
      IF (PRESENT(n_iter)) n_iter = NINT(results(1))
      IF (PRESENT(f)) f = results(2)
      IF (PRESENT(last_f)) last_f = results(3)
      IF (PRESENT(projected_gradient)) &
         projected_gradient = results(4)

   END SUBROUTINE cp_opt_gopt_bcast_res

! **************************************************************************************************
!> \brief goes to the next optimal point (after an optimizer iteration)
!>      returns true if converged
!> \param optimizer the optimizer that goes to the next point
!> \param n_iter ...
!> \param f ...
!> \param last_f ...
!> \param projected_gradient ...
!> \param converged ...
!> \param geo_section ...
!> \param force_env ...
!> \param gopt_param ...
!> \param spgr ...
!> \return ...
!> \par History
!>      01.2020 modified [pcazade]
!> \author Fawzi Mohamed
!>      @version 2.2002
!> \note
!>      if you deactivate convergence control it returns never false
! **************************************************************************************************
   FUNCTION cp_opt_gopt_next(optimizer, n_iter, f, last_f, &
                             projected_gradient, converged, geo_section, force_env, &
                             gopt_param, spgr) RESULT(res)
      TYPE(cp_lbfgs_opt_gopt_type), INTENT(INOUT)        :: optimizer
      INTEGER, INTENT(out), OPTIONAL                     :: n_iter
      REAL(kind=dp), INTENT(out), OPTIONAL               :: f, last_f, projected_gradient
      LOGICAL, INTENT(out)                               :: converged
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(spgr_type), OPTIONAL, POINTER                 :: spgr
      LOGICAL                                            :: res

      ! passes spgr structure if present
      CALL cp_opt_gopt_step(optimizer, n_iter=n_iter, f=f, &
                            last_f=last_f, projected_gradient=projected_gradient, &
                            converged=converged, geo_section=geo_section, &
                            force_env=force_env, gopt_param=gopt_param, spgr=spgr)
      res = (optimizer%status < 40) .AND. .NOT. converged

   END FUNCTION cp_opt_gopt_next

! **************************************************************************************************
!> \brief stops the optimization
!> \param optimizer ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 2.2002
! **************************************************************************************************
   SUBROUTINE cp_opt_gopt_stop(optimizer)
      TYPE(cp_lbfgs_opt_gopt_type), INTENT(INOUT)       :: optimizer

      optimizer%task = 'STOPPED on user request'
      optimizer%status = 4 ! normal exit
      IF (optimizer%master == optimizer%para_env%mepos) THEN
         CALL setulb(SIZE(optimizer%x), optimizer%m, optimizer%x, &
                     optimizer%lower_bound, optimizer%upper_bound, &
                     optimizer%kind_of_bound, optimizer%f, optimizer%gradient, &
                     optimizer%wanted_relative_f_delta, &
                     optimizer%wanted_projected_gradient, optimizer%work_array, &
                     optimizer%i_work_array, optimizer%task, optimizer%print_every, &
                     optimizer%csave, optimizer%lsave, optimizer%isave, &
                     optimizer%dsave, optimizer%trust_radius)
      END IF

   END SUBROUTINE cp_opt_gopt_stop

END MODULE cp_lbfgs_optimizer_gopt
