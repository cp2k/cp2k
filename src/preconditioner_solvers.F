!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief solves the preconditioner, contains to utility function for
!>        fm<->dbcsr transfers, should be moved soon
!> \par History
!>      - [UB] 2009-05-13 Adding stable approximate inverse (full and sparse)
!> \author Joost VandeVondele (09.2002)
! **************************************************************************************************
MODULE preconditioner_solvers
   USE arnoldi_api,                     ONLY: arnoldi_data_type,&
                                              arnoldi_ev,&
                                              deallocate_arnoldi_data,&
                                              get_selected_ritz_val,&
                                              setup_arnoldi_data
   USE bibliography,                    ONLY: Schiffmann2015,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_col_block_sizes, dbcsr_create, dbcsr_distribution, dbcsr_filter, dbcsr_get_info, &
        dbcsr_get_occupation, dbcsr_init_p, dbcsr_p_type, dbcsr_release, dbcsr_row_block_sizes, &
        dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_real_8
   USE input_constants,                 ONLY: ot_precond_solver_default,&
                                              ot_precond_solver_direct,&
                                              ot_precond_solver_inv_chol,&
                                              ot_precond_solver_update
   USE iterate_matrix,                  ONLY: invert_Hotelling
   USE kinds,                           ONLY: dp
   USE preconditioner_types,            ONLY: preconditioner_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_solvers'

   PUBLIC :: solve_preconditioner, transfer_fm_to_dbcsr, transfer_dbcsr_to_fm

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param my_solver_type ...
!> \param preconditioner_env ...
!> \param matrix_s ...
!> \param matrix_h ...
! **************************************************************************************************
   SUBROUTINE solve_preconditioner(my_solver_type, preconditioner_env, matrix_s, &
                                   matrix_h)
      INTEGER                                            :: my_solver_type
      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: matrix_h

      CHARACTER(len=*), PARAMETER :: routineN = 'solve_preconditioner', &
         routineP = moduleN//':'//routineN

      REAL(dp)                                           :: occ_matrix

! here comes the solver

      SELECT CASE (my_solver_type)
      CASE (ot_precond_solver_inv_chol)
         !
         ! compute the full inverse
         preconditioner_env%solver = ot_precond_solver_inv_chol
         CALL make_full_inverse_cholesky(preconditioner_env, matrix_s)
      CASE (ot_precond_solver_direct)
         !
         ! prepare for the direct solver
         preconditioner_env%solver = ot_precond_solver_direct
         CALL make_full_fact_cholesky(preconditioner_env, matrix_s)
      CASE (ot_precond_solver_update)
         !
         ! uses an update of the full inverse (needs to be computed the first time)
         ! make sure preconditioner_env is not destroyed in between
         occ_matrix = 1.0_dp
         IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
            IF (preconditioner_env%condition_num < 0.0_dp) &
               CALL estimate_cond_num(preconditioner_env%sparse_matrix, preconditioner_env%condition_num)
            CALL dbcsr_filter(preconditioner_env%sparse_matrix, &
                              1.0_dp/preconditioner_env%condition_num*0.01_dp)
            occ_matrix = dbcsr_get_occupation(preconditioner_env%sparse_matrix)
         END IF
         ! check whether we are in the first step and if it is a good idea to use cholesky (matrix sparsity)
         IF (preconditioner_env%solver .NE. ot_precond_solver_update .AND. occ_matrix > 0.5_dp) THEN
            preconditioner_env%solver = ot_precond_solver_update
            CALL make_full_inverse_cholesky(preconditioner_env, matrix_s)
         ELSE
            preconditioner_env%solver = ot_precond_solver_update
            CALL make_inverse_update(preconditioner_env, matrix_h)
         END IF
      CASE (ot_precond_solver_default)
         preconditioner_env%solver = ot_precond_solver_default
      CASE DEFAULT
         !
         CPABORT("Doesn't know this type of solver")
      END SELECT

   END SUBROUTINE solve_preconditioner

! **************************************************************************************************
!> \brief Compute the inverse using cholseky factorization
!> \param preconditioner_env ...
!> \param matrix_s ...
! **************************************************************************************************
   SUBROUTINE make_full_inverse_cholesky(preconditioner_env, matrix_s)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_s

      CHARACTER(len=*), PARAMETER :: routineN = 'make_full_inverse_cholesky', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, info
      TYPE(cp_fm_type), POINTER                          :: fm, fm_work

      CALL timeset(routineN, handle)

      ! Maybe we will get a sparse Cholesky at a given point then this can go,
      ! if stuff was stored in fm anyway this simple returns
      CALL transfer_dbcsr_to_fm(preconditioner_env%sparse_matrix, preconditioner_env%fm, &
                                preconditioner_env%para_env, preconditioner_env%ctxt)
      fm => preconditioner_env%fm

      NULLIFY (fm_work)

      CALL cp_fm_create(fm_work, fm%matrix_struct, name="fm_work")
      !
      ! compute the inverse of SPD matrix fm using the Cholesky factorization
      CALL cp_fm_cholesky_decompose(fm, info_out=info)

      !
      ! if fm not SPD we go with the overlap matrix
      IF (info /= 0) THEN
         !
         ! just the overlap matrix
         IF (PRESENT(matrix_s)) THEN
            CALL copy_dbcsr_to_fm(matrix_s, fm)
            CALL cp_fm_cholesky_decompose(fm)
         ELSE
            CALL cp_fm_set_all(fm, alpha=0._dp, beta=1._dp)
         ENDIF
      ENDIF
      CALL cp_fm_cholesky_invert(fm)

      CALL cp_fm_upper_to_full(fm, fm_work)
      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE make_full_inverse_cholesky

! **************************************************************************************************
!> \brief Only perform the factorization, can be used later to solve the linear
!>        system on the fly
!> \param preconditioner_env ...
!> \param matrix_s ...
! **************************************************************************************************
   SUBROUTINE make_full_fact_cholesky(preconditioner_env, matrix_s)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_s

      CHARACTER(len=*), PARAMETER :: routineN = 'make_full_fact_cholesky', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, info_out
      TYPE(cp_fm_type), POINTER                          :: fm

      CALL timeset(routineN, handle)

      ! Maybe we will get a sparse Cholesky at a given point then this can go,
      ! if stuff was stored in fm anyway this simple returns
      CALL transfer_dbcsr_to_fm(preconditioner_env%sparse_matrix, preconditioner_env%fm, &
                                preconditioner_env%para_env, preconditioner_env%ctxt)

      fm => preconditioner_env%fm
      !
      ! compute the inverse of SPD matrix fm using the Cholesky factorization
      CALL cp_fm_cholesky_decompose(fm, info_out=info_out)
      !
      ! if fm not SPD we go with the overlap matrix
      IF (info_out .NE. 0) THEN
         !
         ! just the overlap matrix
         IF (PRESENT(matrix_s)) THEN
            CALL copy_dbcsr_to_fm(matrix_s, fm)
            CALL cp_fm_cholesky_decompose(fm)
         ELSE
            CALL cp_fm_set_all(fm, alpha=0._dp, beta=1._dp)
         ENDIF
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE make_full_fact_cholesky

! **************************************************************************************************
!> \brief computes an approximate inverse using Hotelling iterations
!> \param preconditioner_env ...
!> \param matrix_h as S is not always present this is a safe template for the transfer
! **************************************************************************************************
   SUBROUTINE make_inverse_update(preconditioner_env, matrix_h)
      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(dbcsr_type), POINTER                          :: matrix_h

      CHARACTER(len=*), PARAMETER :: routineN = 'make_inverse_update', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: use_guess
      REAL(KIND=dp)                                      :: filter_eps

      CALL timeset(routineN, handle)
      use_guess = .TRUE.
      !
      ! uses an update of the full inverse (needs to be computed the first time)
      ! make sure preconditioner_env is not destroyed in between

      CALL cite_reference(Schiffmann2015)

      ! Maybe I gonna add a fm Hotelling, ... for now the same as above make sure we are dbcsr
      CALL transfer_fm_to_dbcsr(preconditioner_env%fm, preconditioner_env%sparse_matrix, matrix_h)
      IF (.NOT. ASSOCIATED(preconditioner_env%dbcsr_matrix)) THEN
         use_guess = .FALSE.
         CALL dbcsr_init_p(preconditioner_env%dbcsr_matrix)
         CALL dbcsr_create(preconditioner_env%dbcsr_matrix, "prec_dbcsr", &
                           template=matrix_h, matrix_type=dbcsr_type_no_symmetry)
      END IF

      ! Try to get a reasonbale guess for the filtering threshold
      filter_eps = 1.0_dp/preconditioner_env%condition_num*0.1_dp

      ! Agressive filtering on the initial guess is needed to avoid fill ins and retain sparsity
      CALL dbcsr_filter(preconditioner_env%dbcsr_matrix, filter_eps*100.0_dp)
      ! We don't need a high accuracy for the inverse so 0.4 is reasonable for convergence
      CALL invert_Hotelling(preconditioner_env%dbcsr_matrix, preconditioner_env%sparse_matrix, filter_eps*10.0_dp, &
                            use_inv_as_guess=use_guess, norm_convergence=0.4_dp, filter_eps=filter_eps)

      CALL timestop(handle)

   END SUBROUTINE make_inverse_update

! **************************************************************************************************
!> \brief computes an approximation to the condition number of a matrix using
!>        arnoldi iterations
!> \param matrix ...
!> \param cond_num ...
! **************************************************************************************************
   SUBROUTINE estimate_cond_num(matrix, cond_num)
      TYPE(dbcsr_type), POINTER                          :: matrix
      REAL(KIND=dp)                                      :: cond_num

      CHARACTER(len=*), PARAMETER :: routineN = 'estimate_cond_num', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: max_ev, min_ev
      TYPE(arnoldi_data_type)                            :: my_arnoldi
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrices

      CALL timeset(routineN, handle)

      ! its better to do 2 calculations as the maximum should quickly converge and the minimum won't need iram
      ALLOCATE (matrices(1))
      matrices(1)%matrix => matrix
      ! compute the minimum ev
      CALL setup_arnoldi_data(my_arnoldi, matrices, max_iter=20, threshold=5.0E-4_dp, selection_crit=2, &
                              nval_request=1, nrestarts=15, generalized_ev=.FALSE., iram=.FALSE.)
      CALL arnoldi_ev(matrices, my_arnoldi)
      max_ev = REAL(get_selected_ritz_val(my_arnoldi, 1), dp)
      CALL deallocate_arnoldi_data(my_arnoldi)

      CALL setup_arnoldi_data(my_arnoldi, matrices, max_iter=20, threshold=5.0E-4_dp, selection_crit=3, &
                              nval_request=1, nrestarts=15, generalized_ev=.FALSE., iram=.FALSE.)
      CALL arnoldi_ev(matrices, my_arnoldi)
      min_ev = REAL(get_selected_ritz_val(my_arnoldi, 1), dp)
      CALL deallocate_arnoldi_data(my_arnoldi)

      cond_num = max_ev/min_ev
      DEALLOCATE (matrices)

      CALL timestop(handle)
   END SUBROUTINE estimate_cond_num

! **************************************************************************************************
!> \brief transfers a full matrix to a dbcsr
!> \param fm_matrix a full matrix gets deallocated in the end
!> \param dbcsr_matrix a dbcsr matrix, gets create from a template
!> \param template_mat the template which is used for the structure
! **************************************************************************************************
   SUBROUTINE transfer_fm_to_dbcsr(fm_matrix, dbcsr_matrix, template_mat)

      TYPE(cp_fm_type), POINTER                          :: fm_matrix
      TYPE(dbcsr_type), POINTER                          :: dbcsr_matrix, template_mat

      CHARACTER(len=*), PARAMETER :: routineN = 'transfer_fm_to_dbcsr', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      IF (ASSOCIATED(fm_matrix)) THEN
         IF (.NOT. ASSOCIATED(dbcsr_matrix)) THEN
            CALL dbcsr_init_p(dbcsr_matrix)
            CALL dbcsr_create(dbcsr_matrix, "preconditioner_env%dbcsr_matrix", &
                              dbcsr_distribution(template_mat), dbcsr_type_no_symmetry, &
                              dbcsr_row_block_sizes(template_mat), dbcsr_col_block_sizes(template_mat), &
                              nze=0, data_type=dbcsr_type_real_8)
         ENDIF
!       CALL cp_fm_create(fm_tmp, matrix_struct=fm_matrix%matrix_struct)
!       CALL cp_fm_upper_to_full(fm_matrix, fm_tmp)
         CALL copy_fm_to_dbcsr(fm_matrix, dbcsr_matrix)
!       CALL cp_fm_release(fm_tmp)
         CALL cp_fm_release(fm_matrix)
      END IF

      CALL timestop(handle)

   END SUBROUTINE transfer_fm_to_dbcsr

! **************************************************************************************************
!> \brief transfers a dbcsr to a full matrix
!> \param dbcsr_matrix a dbcsr matrix, gets deallocated at the end
!> \param fm_matrix a full matrix gets created if not yet done
!> \param para_env the para_env
!> \param context the blacs context
! **************************************************************************************************
   SUBROUTINE transfer_dbcsr_to_fm(dbcsr_matrix, fm_matrix, para_env, context)

      TYPE(dbcsr_type), POINTER                          :: dbcsr_matrix
      TYPE(cp_fm_type), POINTER                          :: fm_matrix
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: context

      CHARACTER(len=*), PARAMETER :: routineN = 'transfer_dbcsr_to_fm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, n
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp

      CALL timeset(routineN, handle)
      IF (ASSOCIATED(dbcsr_matrix)) THEN
         NULLIFY (fm_struct_tmp)

         IF (ASSOCIATED(fm_matrix)) THEN
            CALL cp_fm_release(fm_matrix)
         ENDIF

         CALL dbcsr_get_info(dbcsr_matrix, nfullrows_total=n)
         CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=n, ncol_global=n, &
                                  context=context, para_env=para_env)
         CALL cp_fm_create(fm_matrix, fm_struct_tmp)
         CALL cp_fm_struct_release(fm_struct_tmp)

         CALL copy_dbcsr_to_fm(dbcsr_matrix, fm_matrix)
         CALL dbcsr_release(dbcsr_matrix)
         DEALLOCATE (dbcsr_matrix)
      END IF

      CALL timestop(handle)

   END SUBROUTINE transfer_dbcsr_to_fm

END MODULE preconditioner_solvers
