!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The methods which allow to analyze and manipulate the arnoldi procedure
!>        The main routine and this should eb the only public access point for the method
!> \par History
!>       2014.09 created [Florian Schiffmann]
!>       2023.12 Removed support for single-precision [Ole Schuett]
!>       2024.12 Removed support for complex input matrices [Ole Schuett]
!> \author Florian Schiffmann
! **************************************************************************************************
MODULE arnoldi_data_methods
   USE arnoldi_types,                   ONLY: &
        arnoldi_control_type, arnoldi_data_type, arnoldi_env_type, get_control, get_data, &
        get_evals, get_sel_ind, set_control, set_data
   USE arnoldi_vector,                  ONLY: create_col_vec_from_matrix
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_get_data_p, dbcsr_get_info, &
        dbcsr_get_matrix_type, dbcsr_mp_grid_setup, dbcsr_p_type, dbcsr_release, dbcsr_type, &
        dbcsr_type_real_8, dbcsr_type_symmetric
   USE kinds,                           ONLY: dp
   USE util,                            ONLY: sort
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: select_evals, get_selected_ritz_val, arnoldi_is_converged, &
             arnoldi_env_type, get_nrestart, set_arnoldi_initial_vector, &
             setup_arnoldi_env, deallocate_arnoldi_env, get_selected_ritz_vector

CONTAINS

! **************************************************************************************************
!> \brief This routine sets the environment for the arnoldi iteration and
!>        the krylov subspace creation. All simulation parameters have to be given
!>        at this stage so the rest can run fully automated
!>        In addition, this routine allocates the data necessary for
!> \param arnoldi_env this type which gets filled with information and on output contains all
!>                    information necessary to extract whatever the user desires
!> \param matrix vector of matrices, only the first gets used to get some dimensions
!>        and parallel information needed later on
!> \param max_iter maximum dimension of the krylov subspace
!> \param threshold convergence threshold, this is used for both subspace and eigenval
!> \param selection_crit integer defining according to which criterion the
!>        eigenvalues are selected for the subspace
!> \param nval_request for some sel_crit useful, how many eV to select
!> \param nrestarts ...
!> \param generalized_ev ...
!> \param iram ...
! **************************************************************************************************
   SUBROUTINE setup_arnoldi_env(arnoldi_env, matrix, max_iter, threshold, selection_crit, &
                                nval_request, nrestarts, generalized_ev, iram)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      INTEGER                                            :: max_iter
      REAL(dp)                                           :: threshold
      INTEGER                                            :: selection_crit, nval_request, nrestarts
      LOGICAL                                            :: generalized_ev, iram

      CALL setup_arnoldi_control(arnoldi_env, matrix, max_iter, threshold, selection_crit, &
                                 nval_request, nrestarts, generalized_ev, iram)

      CALL setup_arnoldi_data(arnoldi_env, matrix, max_iter)

   END SUBROUTINE setup_arnoldi_env

! **************************************************************************************************
!> \brief Creates the data type for arnoldi, see above for details
!> \param arnoldi_env ...
!> \param matrix ...
!> \param max_iter ...
! **************************************************************************************************
   SUBROUTINE setup_arnoldi_data(arnoldi_env, matrix, max_iter)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      INTEGER                                            :: max_iter

      INTEGER                                            :: data_type, nrow_local
      TYPE(arnoldi_data_type), POINTER                   :: ar_data

      ALLOCATE (ar_data)
      CALL dbcsr_get_info(matrix=matrix(1)%matrix, nfullrows_local=nrow_local, data_type=data_type)
      CPASSERT(data_type == dbcsr_type_real_8)
      ALLOCATE (ar_data%f_vec(nrow_local))
      ALLOCATE (ar_data%x_vec(nrow_local))
      ALLOCATE (ar_data%Hessenberg(max_iter + 1, max_iter))
      ALLOCATE (ar_data%local_history(nrow_local, max_iter))

      ALLOCATE (ar_data%evals(max_iter))
      ALLOCATE (ar_data%revec(max_iter, max_iter))

      CALL set_data(arnoldi_env, ar_data)

   END SUBROUTINE setup_arnoldi_data

! **************************************************************************************************
!> \brief Creates the control type for arnoldi, see above for details
!> \param arnoldi_env ...
!> \param matrix ...
!> \param max_iter ...
!> \param threshold ...
!> \param selection_crit ...
!> \param nval_request ...
!> \param nrestarts ...
!> \param generalized_ev ...
!> \param iram ...
! **************************************************************************************************
   SUBROUTINE setup_arnoldi_control(arnoldi_env, matrix, max_iter, threshold, selection_crit, &
                                    nval_request, nrestarts, generalized_ev, iram)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      INTEGER                                            :: max_iter
      REAL(dp)                                           :: threshold
      INTEGER                                            :: selection_crit, nval_request, nrestarts
      LOGICAL                                            :: generalized_ev, iram

      INTEGER                                            :: group_handle, pcol_handle
      LOGICAL                                            :: subgroups_defined
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(dbcsr_distribution_type)                      :: distri

      ALLOCATE (control)
      ! Fill the information which will later control the arnoldi method and allow synchronization.
      CALL dbcsr_get_info(matrix=matrix(1)%matrix, distribution=distri)
      CALL dbcsr_mp_grid_setup(distri)
      CALL dbcsr_distribution_get(distri, &
                                  group=group_handle, &
                                  mynode=control%myproc, &
                                  subgroups_defined=subgroups_defined, &
                                  pcol_group=pcol_handle)

      CALL control%mp_group%set_handle(group_handle)
      CALL control%pcol_group%set_handle(pcol_handle)

      IF (.NOT. subgroups_defined) &
         CPABORT("arnoldi only with subgroups")

      control%symmetric = .FALSE.
      ! Will need a fix for complex because there it has to be hermitian
      IF (SIZE(matrix) == 1) &
         control%symmetric = dbcsr_get_matrix_type(matrix(1)%matrix) == dbcsr_type_symmetric

      ! Set the control parameters
      control%max_iter = max_iter
      control%current_step = 0
      control%selection_crit = selection_crit
      control%nval_req = nval_request
      control%threshold = threshold
      control%converged = .FALSE.
      control%has_initial_vector = .FALSE.
      control%iram = iram
      control%nrestart = nrestarts
      control%generalized_ev = generalized_ev

      IF (control%nval_req > 1 .AND. control%nrestart > 0 .AND. .NOT. control%iram) &
         CALL cp_abort(__LOCATION__, 'with more than one eigenvalue requested '// &
                       'internal restarting with a previous EVEC is a bad idea, set IRAM or nrestsart=0')

      ! some checks for the generalized EV mode
      IF (control%generalized_ev .AND. selection_crit == 1) &
         CALL cp_abort(__LOCATION__, &
                       'generalized ev can only highest OR lowest EV')
      IF (control%generalized_ev .AND. nval_request .NE. 1) &
         CALL cp_abort(__LOCATION__, &
                       'generalized ev can only compute one EV at the time')
      IF (control%generalized_ev .AND. control%nrestart == 0) &
         CALL cp_abort(__LOCATION__, &
                       'outer loops are mandatory for generalized EV, set nrestart appropriatly')
      IF (SIZE(matrix) .NE. 2 .AND. control%generalized_ev) &
         CALL cp_abort(__LOCATION__, &
                       'generalized ev needs exactly two matrices as input (2nd is the metric)')

      ALLOCATE (control%selected_ind(max_iter))
      CALL set_control(arnoldi_env, control)

   END SUBROUTINE setup_arnoldi_control

! **************************************************************************************************
!> \brief ...
!> \param arnoldi_env ...
!> \param ind ...
!> \param matrix ...
!> \param vector ...
! **************************************************************************************************
   SUBROUTINE get_selected_ritz_vector(arnoldi_env, ind, matrix, vector)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      INTEGER                                            :: ind
      TYPE(dbcsr_type)                                   :: matrix, vector

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:)             :: ritz_v
      INTEGER                                            :: i, myind, sspace_size, vsize
      INTEGER, DIMENSION(:), POINTER                     :: selected_ind
      REAL(kind=dp), DIMENSION(:), POINTER               :: data_vec
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type), POINTER                   :: ar_data

      control => get_control(arnoldi_env)
      selected_ind => get_sel_ind(arnoldi_env)
      ar_data => get_data(arnoldi_env)
      sspace_size = get_subsp_size(arnoldi_env)
      vsize = SIZE(ar_data%f_vec)
      myind = selected_ind(ind)
      ALLOCATE (ritz_v(vsize))
      ritz_v = CMPLX(0.0, 0.0, dp)

      CALL dbcsr_release(vector)
      CALL create_col_vec_from_matrix(vector, matrix, 1)
      IF (control%local_comp) THEN
         DO i = 1, sspace_size
            ritz_v(:) = ritz_v(:) + ar_data%local_history(:, i)*ar_data%revec(i, myind)
         END DO
         data_vec => dbcsr_get_data_p(vector, select_data_type=0.0_dp)
         ! is a bit odd but ritz_v is always complex and matrix type determines where it goes
         ! again I hope the user knows what is required
         data_vec(1:vsize) = REAL(ritz_v(1:vsize), KIND=dp)
      END IF

      DEALLOCATE (ritz_v)

   END SUBROUTINE get_selected_ritz_vector

! **************************************************************************************************
!> \brief Deallocate the data in arnoldi_env
!> \param arnoldi_env ...
! **************************************************************************************************
   SUBROUTINE deallocate_arnoldi_env(arnoldi_env)
      TYPE(arnoldi_env_type)                             :: arnoldi_env

      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type), POINTER                   :: ar_data

      ar_data => get_data(arnoldi_env)
      IF (ASSOCIATED(ar_data%f_vec)) DEALLOCATE (ar_data%f_vec)
      IF (ASSOCIATED(ar_data%x_vec)) DEALLOCATE (ar_data%x_vec)
      IF (ASSOCIATED(ar_data%Hessenberg)) DEALLOCATE (ar_data%Hessenberg)
      IF (ASSOCIATED(ar_data%local_history)) DEALLOCATE (ar_data%local_history)
      IF (ASSOCIATED(ar_data%evals)) DEALLOCATE (ar_data%evals)
      IF (ASSOCIATED(ar_data%revec)) DEALLOCATE (ar_data%revec)
      DEALLOCATE (ar_data)

      control => get_control(arnoldi_env)
      DEALLOCATE (control%selected_ind)
      DEALLOCATE (control)

   END SUBROUTINE deallocate_arnoldi_env

! **************************************************************************************************
!> \brief perform the selection of eigenvalues, fills the selected_ind array
!> \param arnoldi_env ...
! **************************************************************************************************
   SUBROUTINE select_evals(arnoldi_env)
      TYPE(arnoldi_env_type)                             :: arnoldi_env

      INTEGER                                            :: i, last_el, my_crit, my_ind
      REAL(dp)                                           :: convergence
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type), POINTER                   :: ar_data

      control => get_control(arnoldi_env)
      ar_data => get_data(arnoldi_env)

      last_el = control%current_step
      convergence = REAL(0.0, dp)
      my_crit = control%selection_crit
      control%nval_out = MIN(control%nval_req, control%current_step)
      SELECT CASE (my_crit)
         ! minimum and maximum real eval
      CASE (1)
         CALL index_min_max_real_eval(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
         ! n maximum real eval
      CASE (2)
         CALL index_nmax_real_eval(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
         ! n minimum real eval
      CASE (3)
         CALL index_nmin_real_eval(ar_data%evals, control%current_step, control%selected_ind, control%nval_out)
      CASE DEFAULT
         CPABORT("unknown selection index")
      END SELECT
      ! test whether we are converged
      DO i = 1, control%nval_out
         my_ind = control%selected_ind(i)
         convergence = MAX(convergence, &
                           ABS(ar_data%revec(last_el, my_ind)*ar_data%Hessenberg(last_el + 1, last_el)))
      END DO
      control%converged = convergence .LT. control%threshold

   END SUBROUTINE select_evals

! **************************************************************************************************
!> \brief set a new selection type, if you notice you didn't like the initial one
!> \param arnoldi_env ...
!> \param itype ...
! **************************************************************************************************
   SUBROUTINE set_eval_selection(arnoldi_env, itype)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      INTEGER                                            :: itype

      TYPE(arnoldi_control_type), POINTER                :: control

      control => get_control(arnoldi_env)
      control%selection_crit = itype
   END SUBROUTINE set_eval_selection

! **************************************************************************************************
!> \brief returns the number of restarts allowed for arnoldi
!> \param arnoldi_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_nrestart(arnoldi_env) RESULT(nrestart)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      INTEGER                                            :: nrestart

      TYPE(arnoldi_control_type), POINTER                :: control

      control => get_control(arnoldi_env)
      nrestart = control%nrestart

   END FUNCTION get_nrestart

! **************************************************************************************************
!> \brief get the number of eigenvalues matching the search criterion
!> \param arnoldi_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_nval_out(arnoldi_env) RESULT(nval_out)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      INTEGER                                            :: nval_out

      TYPE(arnoldi_control_type), POINTER                :: control

      control => get_control(arnoldi_env)
      nval_out = control%nval_out

   END FUNCTION get_nval_out

! **************************************************************************************************
!> \brief get dimension of the krylov space. Can be less than max_iter if subspace converged early
!> \param arnoldi_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_subsp_size(arnoldi_env) RESULT(current_step)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      INTEGER                                            :: current_step

      TYPE(arnoldi_control_type), POINTER                :: control

      control => get_control(arnoldi_env)
      current_step = control%current_step

   END FUNCTION get_subsp_size

! **************************************************************************************************
!> \brief Find out whether the method with the current search criterion is converged
!> \param arnoldi_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION arnoldi_is_converged(arnoldi_env) RESULT(converged)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      LOGICAL                                            :: converged

      TYPE(arnoldi_control_type), POINTER                :: control

      control => get_control(arnoldi_env)
      converged = control%converged

   END FUNCTION arnoldi_is_converged

! **************************************************************************************************
!> \brief get a single specific Ritz value from the set of selected
!> \param arnoldi_env ...
!> \param ind ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_selected_ritz_val(arnoldi_env, ind) RESULT(eval_out)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      INTEGER                                            :: ind
      COMPLEX(dp)                                        :: eval_out

      COMPLEX(dp), DIMENSION(:), POINTER                 :: evals
      INTEGER                                            :: ev_ind
      INTEGER, DIMENSION(:), POINTER                     :: selected_ind

      IF (ind .GT. get_nval_out(arnoldi_env)) &
         CPABORT('outside range of indexed evals')

      selected_ind => get_sel_ind(arnoldi_env)
      ev_ind = selected_ind(ind)
      evals => get_evals(arnoldi_env)
      eval_out = evals(ev_ind)

   END FUNCTION get_selected_ritz_val

! **************************************************************************************************
!> \brief Get all Ritz values of the selected set. eval_out has to be allocated
!>        at least the size of get_neval_out()
!> \param arnoldi_env ...
!> \param eval_out ...
! **************************************************************************************************
   SUBROUTINE get_all_selected_ritz_val(arnoldi_env, eval_out)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      COMPLEX(dp), DIMENSION(:)                          :: eval_out

      COMPLEX(dp), DIMENSION(:), POINTER                 :: evals
      INTEGER                                            :: ev_ind, ind
      INTEGER, DIMENSION(:), POINTER                     :: selected_ind

      NULLIFY (evals)
      IF (SIZE(eval_out) .LT. get_nval_out(arnoldi_env)) &
         CPABORT('array for eval output too small')
      selected_ind => get_sel_ind(arnoldi_env)

      evals => get_evals(arnoldi_env)

      DO ind = 1, get_nval_out(arnoldi_env)
         ev_ind = selected_ind(ind)
         eval_out(ind) = evals(ev_ind)
      END DO

   END SUBROUTINE get_all_selected_ritz_val

! **************************************************************************************************
!> \brief ...
!> \param arnoldi_env ...
!> \param vector ...
! **************************************************************************************************
   SUBROUTINE set_arnoldi_initial_vector(arnoldi_env, vector)
      TYPE(arnoldi_env_type)                             :: arnoldi_env
      TYPE(dbcsr_type)                                   :: vector

      INTEGER                                            :: ncol_local, nrow_local
      REAL(kind=dp), DIMENSION(:), POINTER               :: data_vec
      TYPE(arnoldi_control_type), POINTER                :: control
      TYPE(arnoldi_data_type), POINTER                   :: ar_data

      control => get_control(arnoldi_env)
      control%has_initial_vector = .TRUE.
      ar_data => get_data(arnoldi_env)

      CALL dbcsr_get_info(matrix=vector, nfullrows_local=nrow_local, nfullcols_local=ncol_local)
      data_vec => dbcsr_get_data_p(vector, select_data_type=0.0_dp)
      IF (nrow_local*ncol_local > 0) ar_data%f_vec(1:nrow_local) = data_vec(1:nrow_local)

   END SUBROUTINE set_arnoldi_initial_vector

!!! Here come the methods handling the selection of eigenvalues and eigenvectors !!!
!!! If you want a personal method, simply created a Subroutine returning the index
!!! array selected ind which contains as the first nval_out entries the index of the evals

! **************************************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! **************************************************************************************************
   SUBROUTINE index_min_max_real_eval(evals, current_step, selected_ind, neval)
      COMPLEX(dp), DIMENSION(:)                          :: evals
      INTEGER, INTENT(IN)                                :: current_step
      INTEGER, DIMENSION(:)                              :: selected_ind
      INTEGER                                            :: neval

      INTEGER                                            :: i
      INTEGER, DIMENSION(current_step)                   :: indexing
      REAL(dp), DIMENSION(current_step)                  :: tmp_array

      neval = 0
      selected_ind = 0
      tmp_array(1:current_step) = REAL(evals(1:current_step), dp)
      CALL sort(tmp_array, current_step, indexing)
      DO i = 1, current_step
         IF (ABS(AIMAG(evals(indexing(i)))) < EPSILON(0.0_dp)) THEN
            selected_ind(1) = indexing(i)
            neval = neval + 1
            EXIT
         END IF
      END DO
      DO i = current_step, 1, -1
         IF (ABS(AIMAG(evals(indexing(i)))) < EPSILON(0.0_dp)) THEN
            selected_ind(2) = indexing(i)
            neval = neval + 1
            EXIT
         END IF
      END DO

   END SUBROUTINE index_min_max_real_eval

! **************************************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! **************************************************************************************************
   SUBROUTINE index_nmax_real_eval(evals, current_step, selected_ind, neval)
      COMPLEX(dp), DIMENSION(:)                          :: evals
      INTEGER, INTENT(IN)                                :: current_step
      INTEGER, DIMENSION(:)                              :: selected_ind
      INTEGER                                            :: neval

      INTEGER                                            :: i, nlimit
      INTEGER, DIMENSION(current_step)                   :: indexing
      REAL(dp), DIMENSION(current_step)                  :: tmp_array

      nlimit = neval; neval = 0
      selected_ind = 0
      tmp_array(1:current_step) = REAL(evals(1:current_step), dp)
      CALL sort(tmp_array, current_step, indexing)
      DO i = 1, current_step
         IF (ABS(AIMAG(evals(indexing(current_step + 1 - i)))) < EPSILON(0.0_dp)) THEN
            selected_ind(i) = indexing(current_step + 1 - i)
            neval = neval + 1
            IF (neval == nlimit) EXIT
         END IF
      END DO

   END SUBROUTINE index_nmax_real_eval

! **************************************************************************************************
!> \brief ...
!> \param evals ...
!> \param current_step ...
!> \param selected_ind ...
!> \param neval ...
! **************************************************************************************************
   SUBROUTINE index_nmin_real_eval(evals, current_step, selected_ind, neval)
      COMPLEX(dp), DIMENSION(:)                          :: evals
      INTEGER, INTENT(IN)                                :: current_step
      INTEGER, DIMENSION(:)                              :: selected_ind
      INTEGER                                            :: neval

      INTEGER                                            :: i, nlimit
      INTEGER, DIMENSION(current_step)                   :: indexing
      REAL(dp), DIMENSION(current_step)                  :: tmp_array

      nlimit = neval; neval = 0
      selected_ind = 0
      tmp_array(1:current_step) = REAL(evals(1:current_step), dp)
      CALL sort(tmp_array, current_step, indexing)
      DO i = 1, current_step
         IF (ABS(AIMAG(evals(indexing(i)))) < EPSILON(0.0_dp)) THEN
            selected_ind(i) = indexing(i)
            neval = neval + 1
            IF (neval == nlimit) EXIT
         END IF
      END DO

   END SUBROUTINE index_nmin_real_eval

END MODULE arnoldi_data_methods
