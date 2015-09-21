!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Cayley transformation methods
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE ct_methods
  USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                             cp_dbcsr_cholesky_invert
  USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_desymmetrize, cp_dbcsr_distribution, cp_dbcsr_filter, &
       cp_dbcsr_finalize, cp_dbcsr_frobenius_norm, &
       cp_dbcsr_function_of_elements, cp_dbcsr_get_diag, cp_dbcsr_get_info, &
       cp_dbcsr_get_stored_coordinates, cp_dbcsr_hadamard_product, &
       cp_dbcsr_init, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_nblkcols_total, &
       cp_dbcsr_nblkrows_total, cp_dbcsr_norm, cp_dbcsr_release, &
       cp_dbcsr_reserve_block2d, cp_dbcsr_scale, cp_dbcsr_set, &
       cp_dbcsr_set_diag, cp_dbcsr_trace, cp_dbcsr_transposed, cp_dbcsr_type, &
       cp_dbcsr_work_create, dbcsr_distribution_mp, dbcsr_func_inverse, &
       dbcsr_mp_mynode, dbcsr_norm_maxabsnorm, dbcsr_type_no_symmetry
  USE ct_types,                        ONLY: ct_step_env_type
  USE input_constants,                 ONLY: &
       cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, cg_hager_zhang, &
       cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, cg_zero, &
       tensor_orthogonal, tensor_up_down
  USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ct_methods'

  ! Public subroutines
  PUBLIC :: ct_step_execute, analytic_line_search, diagonalize_diagonal_blocks

CONTAINS

! *****************************************************************************
!> \brief Performs Cayley transformation
!> \param cts_env ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE ct_step_execute(cts_env)
  
    TYPE(ct_step_env_type)                   :: cts_env

    CHARACTER(len=*), PARAMETER :: routineN = 'ct_step_execute', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n, &
                                                preconditioner_type, unit_nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: gap_estimate, safety_margin
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: evals
    TYPE(cp_dbcsr_type) :: matrix_pp, matrix_pq, matrix_qp, matrix_qp_save, &
      matrix_qq, oo1, oo1_sqrt, oo1_sqrt_inv, t_corr, tmp1, u_pp, u_qq
    TYPE(cp_logger_type), POINTER            :: logger

!TYPE(cp_dbcsr_type)                :: rst_x1, rst_x2 
!REAL(KIND=dp)                      :: ener_tmp 
!TYPE(cp_dbcsr_iterator)            :: iter
!INTEGER                            :: iblock_row,iblock_col,&
!                                      iblock_row_size,iblock_col_size
!REAL(KIND=dp), DIMENSION(:,:), POINTER :: data_p

    CALL timeset(routineN,handle)

    failure=.FALSE.
    logger => cp_get_default_logger()
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! check if all input is in place and flags are consistent
    IF (cts_env%update_q.AND.(.NOT.cts_env%update_p)) THEN
       CPABORT("q-update is possible only with p-update")
    ENDIF

    IF (cts_env%tensor_type.eq.tensor_up_down) THEN
       CPABORT("riccati is not implemented for biorthogonal basis")
    ENDIF

    IF (.NOT.ASSOCIATED(cts_env%matrix_ks)) THEN
       CPABORT("KS matrix is not associated")
    ENDIF

    IF (cts_env%use_virt_orbs.AND.(.NOT.cts_env%use_occ_orbs)) THEN
       CPABORT("virtual orbs can be used only with occupied orbs")
    ENDIF

    IF (cts_env%use_occ_orbs) THEN
     IF (.NOT.ASSOCIATED(cts_env%matrix_t)) THEN
       CPABORT("T matrix is not associated")
     ENDIF
     IF (.NOT.ASSOCIATED(cts_env%matrix_qp_template)) THEN
       CPABORT("QP template is not associated")
     ENDIF
     IF (.NOT.ASSOCIATED(cts_env%matrix_pq_template)) THEN
       CPABORT("PQ template is not associated")
     ENDIF
    ENDIF

    IF (cts_env%use_virt_orbs) THEN
     IF (.NOT.ASSOCIATED(cts_env%matrix_v)) THEN
       CPABORT("V matrix is not associated")
     ENDIF
    ELSE
     IF (.NOT.ASSOCIATED(cts_env%matrix_p)) THEN
       CPABORT("P matrix is not associated")
     ENDIF
    ENDIF

    IF (cts_env%tensor_type.NE.tensor_up_down.AND.&
       cts_env%tensor_type.NE.tensor_orthogonal) THEN
       CPABORT("illegal tensor flag")
    ENDIF

    ! start real calculations
    IF (cts_env%use_occ_orbs) THEN

     ! create matrices for various ks blocks
     CALL cp_dbcsr_init(matrix_pp)
     CALL cp_dbcsr_create(matrix_pp,&
            template=cts_env%p_index_up,&
            matrix_type=dbcsr_type_no_symmetry)
     CALL cp_dbcsr_init(matrix_qp)
     CALL cp_dbcsr_create(matrix_qp,&
            template=cts_env%matrix_qp_template,&
            matrix_type=dbcsr_type_no_symmetry)
     CALL cp_dbcsr_init(matrix_qq)
     CALL cp_dbcsr_create(matrix_qq,&
            template=cts_env%q_index_up,&
            matrix_type=dbcsr_type_no_symmetry)
     CALL cp_dbcsr_init(matrix_pq)
     CALL cp_dbcsr_create(matrix_pq,&
            template=cts_env%matrix_pq_template,&
            matrix_type=dbcsr_type_no_symmetry)

     ! create the residue matrix
     CALL cp_dbcsr_init(cts_env%matrix_res)
     CALL cp_dbcsr_create(cts_env%matrix_res,&
            template=cts_env%matrix_qp_template)
     
     CALL assemble_ks_qp_blocks(cts_env%matrix_ks,&
                                cts_env%matrix_p,&
                                cts_env%matrix_t,&
                                cts_env%matrix_v,&
                                cts_env%q_index_down,&
                                cts_env%p_index_up,&
                                cts_env%q_index_up,&
                                matrix_pp,&
                                matrix_qq,&
                                matrix_qp,&
                                matrix_pq,&
                                cts_env%tensor_type,&
                                cts_env%use_virt_orbs,&
                                cts_env%eps_filter)
   
     ! create a matrix of single-excitation amplitudes
     CALL cp_dbcsr_init(cts_env%matrix_x)
     CALL cp_dbcsr_create(cts_env%matrix_x,&
            template=cts_env%matrix_qp_template)
     IF (ASSOCIATED(cts_env%matrix_x_guess)) THEN
        CALL cp_dbcsr_copy(cts_env%matrix_x,&
                cts_env%matrix_x_guess)
        IF (cts_env%tensor_type.eq.tensor_orthogonal) THEN
           ! bring x from contravariant-covariant representation
           ! to the orthogonal/cholesky representation
           ! use res as temporary storage
           CALL cp_dbcsr_multiply("N","N",1.0_dp,cts_env%q_index_down,&
                                  cts_env%matrix_x,0.0_dp,cts_env%matrix_res,&
                                  filter_eps=cts_env%eps_filter)
           CALL cp_dbcsr_multiply("N","N",1.0_dp,cts_env%matrix_res,&
                                  cts_env%p_index_up,0.0_dp,&
                                  cts_env%matrix_x,&
                                  filter_eps=cts_env%eps_filter)
        ENDIF
     ELSE
        ! set amplitudes to zero
        CALL cp_dbcsr_set(cts_env%matrix_x,0.0_dp)
     ENDIF

     !SELECT CASE (cts_env%preconditioner_type)
     !CASE (prec_eigenvector_blocks,prec_eigenvector_full)
     preconditioner_type=1
     safety_margin=2.0_dp
     gap_estimate=0.0001_dp
     SELECT CASE (preconditioner_type)
     CASE (1,2)
!RZK-warning diagonalization works only with orthogonal tensor!!!
        ! find a better basis by diagonalizing diagonal blocks
        ! first pp
        CALL cp_dbcsr_init(u_pp)
        CALL cp_dbcsr_create(u_pp,template=matrix_pp,&
               matrix_type=dbcsr_type_no_symmetry)
        !IF (cts_env%preconditioner_type.eq.prec_eigenvector_full) THEN
        IF (.TRUE.) THEN
           CALL cp_dbcsr_get_info(matrix_pp, nfullrows_total=n )
           ALLOCATE(evals(n))
           CALL cp_dbcsr_syevd(matrix_pp,u_pp,evals,&
                   cts_env%para_env,cts_env%blacs_env)
           DEALLOCATE(evals)
        ELSE
           CALL diagonalize_diagonal_blocks(matrix_pp,u_pp)
        ENDIF
        ! and now qq
        CALL cp_dbcsr_init(u_qq)
        CALL cp_dbcsr_create(u_qq,template=matrix_qq,&
               matrix_type=dbcsr_type_no_symmetry)
        !IF (cts_env%preconditioner_type.eq.prec_eigenvector_full) THEN
        IF (.TRUE.) THEN
           CALL cp_dbcsr_get_info(matrix_qq, nfullrows_total=n )
           ALLOCATE(evals(n))
           CALL cp_dbcsr_syevd(matrix_qq,u_qq,evals,&
                   cts_env%para_env,cts_env%blacs_env)
           DEALLOCATE(evals)
        ELSE
           CALL diagonalize_diagonal_blocks(matrix_qq,u_qq)
        ENDIF
   
        ! apply the transformation to all matrices
        CALL matrix_forward_transform(matrix_pp,u_pp,u_pp,&
                cts_env%eps_filter)
        CALL matrix_forward_transform(matrix_qq,u_qq,u_qq,&
                cts_env%eps_filter)
        CALL matrix_forward_transform(matrix_qp,u_qq,u_pp,&
                cts_env%eps_filter)
        CALL matrix_forward_transform(matrix_pq,u_pp,u_qq,&
                cts_env%eps_filter)
        CALL matrix_forward_transform(cts_env%matrix_x,u_qq,u_pp,&
                cts_env%eps_filter)
   
        IF ( cts_env%max_iter.ge.0 ) THEN
       
           CALL solve_riccati_equation(&
                   pp=matrix_pp,&
                   qq=matrix_qq,&
                   qp=matrix_qp,&
                   pq=matrix_pq,&
                   x=cts_env%matrix_x,&
                   res=cts_env%matrix_res,&
                   neglect_quadratic_term=cts_env%neglect_quadratic_term,&
                   conjugator=cts_env%conjugator,&
                   max_iter=cts_env%max_iter,&
                   eps_convergence=cts_env%eps_convergence,&
                   eps_filter=cts_env%eps_filter,&
                   converged=cts_env%converged)
       
           IF (cts_env%converged) THEN
             !IF (unit_nr>0) THEN
             !   WRITE(unit_nr,*)
             !   WRITE(unit_nr,'(T6,A)') &
             !         "RICCATI equations solved"
             !   CALL m_flush(unit_nr)
             !ENDIF
           ELSE
             CPABORT("RICCATI: CG algorithm has NOT converged")
           ENDIF
   
        ENDIF
   
        IF (cts_env%calculate_energy_corr) THEN
   
          CALL cp_dbcsr_trace(matrix_qp,cts_env%matrix_x,&
             cts_env%energy_correction,"T","N")
   
        ENDIF
   
        CALL cp_dbcsr_release(matrix_pp)
        CALL cp_dbcsr_release(matrix_qp)
        CALL cp_dbcsr_release(matrix_qq)
        CALL cp_dbcsr_release(matrix_pq)
   
        ! back-transform to the original basis
        CALL matrix_backward_transform(cts_env%matrix_x,u_qq,&
                u_pp,cts_env%eps_filter)
   
        CALL cp_dbcsr_release(u_qq)
        CALL cp_dbcsr_release(u_pp)
   
     !CASE (prec_cholesky_inverse)
     CASE (3)

! RZK-warning implemented only for orthogonal tensors!!!
! generalization to up_down should be easy
        CALL cp_dbcsr_init(u_pp)
        CALL cp_dbcsr_create(u_pp,template=matrix_pp,&
               matrix_type=dbcsr_type_no_symmetry)
        CALL cp_dbcsr_copy(u_pp,matrix_pp)
        CALL cp_dbcsr_scale(u_pp,-1.0_dp)
        CALL cp_dbcsr_add_on_diag(u_pp,&
                ABS(safety_margin*gap_estimate))
        CALL cp_dbcsr_cholesky_decompose(u_pp,&
                para_env=cts_env%para_env,&
                blacs_env=cts_env%blacs_env)
        CALL cp_dbcsr_cholesky_invert(u_pp,&
                para_env=cts_env%para_env,&
                blacs_env=cts_env%blacs_env,&
                upper_to_full=.TRUE.)
        !CALL cp_dbcsr_scale(u_pp,-1.0_dp)
        
        CALL cp_dbcsr_init(u_qq)
        CALL cp_dbcsr_create(u_qq,template=matrix_qq,&
               matrix_type=dbcsr_type_no_symmetry)
        CALL cp_dbcsr_copy(u_qq,matrix_qq)
        CALL cp_dbcsr_add_on_diag(u_qq,&
                ABS(safety_margin*gap_estimate))
        CALL cp_dbcsr_cholesky_decompose(u_qq,&
                para_env=cts_env%para_env,&
                blacs_env=cts_env%blacs_env)
        CALL cp_dbcsr_cholesky_invert(u_qq,&
                para_env=cts_env%para_env,&
                blacs_env=cts_env%blacs_env,&
                upper_to_full=.TRUE.)
        
        ! transform all riccati matrices (left-right preconditioner)
        CALL cp_dbcsr_init(tmp1)
        CALL cp_dbcsr_create(tmp1,template=matrix_qq,&
               matrix_type=dbcsr_type_no_symmetry)
        CALL cp_dbcsr_multiply("N","N",1.0_dp,u_qq,&
                matrix_qq,0.0_dp,tmp1,&
                filter_eps=cts_env%eps_filter)
        CALL cp_dbcsr_copy(matrix_qq,tmp1)
        CALL cp_dbcsr_release(tmp1)
        
        CALL cp_dbcsr_init(tmp1)
        CALL cp_dbcsr_create(tmp1,template=matrix_pp,&
               matrix_type=dbcsr_type_no_symmetry)
        CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_pp,&
                u_pp,0.0_dp,tmp1,&
                filter_eps=cts_env%eps_filter)
        CALL cp_dbcsr_copy(matrix_pp,tmp1)
        CALL cp_dbcsr_release(tmp1)

        CALL cp_dbcsr_init(matrix_qp_save)
        CALL cp_dbcsr_create(matrix_qp_save,template=matrix_qp,&
               matrix_type=dbcsr_type_no_symmetry)
        CALL cp_dbcsr_copy(matrix_qp_save,matrix_qp)
        
        CALL cp_dbcsr_init(tmp1)
        CALL cp_dbcsr_create(tmp1,template=matrix_qp,&
               matrix_type=dbcsr_type_no_symmetry)
        CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_qp,&
                u_pp,0.0_dp,tmp1,&
                filter_eps=cts_env%eps_filter)
        CALL cp_dbcsr_multiply("N","N",1.0_dp,u_qq,tmp1,&
                0.0_dp,matrix_qp,&
                filter_eps=cts_env%eps_filter)
        CALL cp_dbcsr_release(tmp1)
!CALL cp_dbcsr_print(matrix_qq)
!CALL cp_dbcsr_print(matrix_qp)
!CALL cp_dbcsr_print(matrix_pp)

        IF ( cts_env%max_iter.ge.0 ) THEN
       
           CALL solve_riccati_equation(&
                   pp=matrix_pp,&
                   qq=matrix_qq,&
                   qp=matrix_qp,&
                   pq=matrix_pq,&
                   oo=u_pp,&
                   vv=u_qq,&
                   x=cts_env%matrix_x,&
                   res=cts_env%matrix_res,&
                   neglect_quadratic_term=cts_env%neglect_quadratic_term,&
                   conjugator=cts_env%conjugator,&
                   max_iter=cts_env%max_iter,&
                   eps_convergence=cts_env%eps_convergence,&
                   eps_filter=cts_env%eps_filter,&
                   converged=cts_env%converged)
       
           IF (cts_env%converged) THEN
             !IF (unit_nr>0) THEN
             !   WRITE(unit_nr,*)
             !   WRITE(unit_nr,'(T6,A)') &
             !         "RICCATI equations solved"
             !   CALL m_flush(unit_nr)
             !ENDIF
           ELSE
             CPABORT("RICCATI: CG algorithm has NOT converged")
           ENDIF
   
        ENDIF
   
        IF (cts_env%calculate_energy_corr) THEN
   
          CALL cp_dbcsr_trace(matrix_qp_save,cts_env%matrix_x,&
             cts_env%energy_correction,"T","N")
   
        ENDIF
        CALL cp_dbcsr_release(matrix_qp_save)
   
        CALL cp_dbcsr_release(matrix_pp)
        CALL cp_dbcsr_release(matrix_qp)
        CALL cp_dbcsr_release(matrix_qq)
        CALL cp_dbcsr_release(matrix_pq)
   
        CALL cp_dbcsr_release(u_qq)
        CALL cp_dbcsr_release(u_pp)
   

     CASE DEFAULT
        CPABORT("illegal preconditioner type")
     END SELECT ! preconditioner type

     IF (cts_env%update_p) THEN
       
        IF (cts_env%tensor_type.eq.tensor_up_down) THEN
           CPABORT("orbital update is NYI for this tensor type")
        ENDIF

        ! transform occupied orbitals
        ! in a way that preserves the overlap metric
        CALL cp_dbcsr_init(oo1)
        CALL cp_dbcsr_create(oo1,&
            template=cts_env%p_index_up,&
            matrix_type=dbcsr_type_no_symmetry)
        CALL cp_dbcsr_init(oo1_sqrt_inv)
        CALL cp_dbcsr_create(oo1_sqrt_inv,&
            template=oo1)
        CALL cp_dbcsr_init(oo1_sqrt)
        CALL cp_dbcsr_create(oo1_sqrt,&
            template=oo1)
       
        ! Compute (1+tr(X).X)^(-1/2)_up_down
        CALL cp_dbcsr_multiply("T","N",1.0_dp,cts_env%matrix_x,&
                               cts_env%matrix_x,0.0_dp,oo1,&
                               filter_eps=cts_env%eps_filter)
        CALL cp_dbcsr_add_on_diag(oo1,1.0_dp)
        CALL matrix_sqrt_Newton_Schulz(oo1_sqrt,&
                                       oo1_sqrt_inv,&
                                       oo1,&
                                       !if cholesky is used then sqrt 
                                       !guess cannot be provided 
                                       !matrix_sqrt_inv_guess=cts_env%p_index_up,&
                                       !matrix_sqrt_guess=cts_env%p_index_down,&
                                       threshold=cts_env%eps_filter,&
                                       order=cts_env%order_lanczos,&
                                       eps_lanczos=cts_env%eps_lancsoz,&
                                       max_iter_lanczos=cts_env%max_iter_lanczos)
        CALL cp_dbcsr_multiply("N","N",1.0_dp,cts_env%p_index_up,&
                               oo1_sqrt_inv,0.0_dp,oo1,&
                               filter_eps=cts_env%eps_filter)
        CALL cp_dbcsr_multiply("N","N",1.0_dp,oo1,&
                               cts_env%p_index_down,0.0_dp,oo1_sqrt,&
                               filter_eps=cts_env%eps_filter)
        CALL cp_dbcsr_release(oo1)
        CALL cp_dbcsr_release(oo1_sqrt_inv)
        
        ! bring x to contravariant-covariant representation now
        CALL cp_dbcsr_init(matrix_qp)
        CALL cp_dbcsr_create(matrix_qp,&
                template=cts_env%matrix_qp_template,&
                matrix_type=dbcsr_type_no_symmetry)
        CALL cp_dbcsr_multiply("N","N",1.0_dp,cts_env%q_index_up,&
                               cts_env%matrix_x,0.0_dp,matrix_qp,&
                               filter_eps=cts_env%eps_filter)
        CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_qp,&
                               cts_env%p_index_down,0.0_dp,&
                               cts_env%matrix_x,&
                               filter_eps=cts_env%eps_filter)
        CALL cp_dbcsr_release(matrix_qp)

        ! update T=T+X or T=T+V.X (whichever is appropriate)
        CALL cp_dbcsr_init(t_corr)
        CALL cp_dbcsr_create(t_corr,template=cts_env%matrix_t)
        IF (cts_env%use_virt_orbs) THEN
           CALL cp_dbcsr_multiply("N","N",1.0_dp,cts_env%matrix_v,&
                   cts_env%matrix_x,0.0_dp,t_corr,&
                   filter_eps=cts_env%eps_filter)
           CALL cp_dbcsr_add(cts_env%matrix_t,t_corr,&
                   1.0_dp,1.0_dp)
        ELSE
           CALL cp_dbcsr_add(cts_env%matrix_t,cts_env%matrix_x,&
                   1.0_dp,1.0_dp)
        ENDIF
        ! adjust T so the metric is preserved: T=(T+X).(1+tr(X).X)^(-1/2)
        CALL cp_dbcsr_multiply("N","N",1.0_dp,cts_env%matrix_t,oo1_sqrt,&
                0.0_dp,t_corr,filter_eps=cts_env%eps_filter)
        CALL cp_dbcsr_copy(cts_env%matrix_t,t_corr)
        
        CALL cp_dbcsr_release(t_corr)
        CALL cp_dbcsr_release(oo1_sqrt)

     ELSE ! do not update p
        
        IF (cts_env%tensor_type.eq.tensor_orthogonal) THEN
           ! bring x to contravariant-covariant representation
           CALL cp_dbcsr_init(matrix_qp)
           CALL cp_dbcsr_create(matrix_qp,&
                   template=cts_env%matrix_qp_template,&
                   matrix_type=dbcsr_type_no_symmetry)
           CALL cp_dbcsr_multiply("N","N",1.0_dp,cts_env%q_index_up,&
                                  cts_env%matrix_x,0.0_dp,matrix_qp,&
                                  filter_eps=cts_env%eps_filter)
           CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_qp,&
                                  cts_env%p_index_down,0.0_dp,&
                                  cts_env%matrix_x,&
                                  filter_eps=cts_env%eps_filter)
           CALL cp_dbcsr_release(matrix_qp)
        ENDIF

     ENDIF
     
    ELSE
     CPABORT("illegal occ option")
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE ct_step_execute

! *****************************************************************************
!> \brief computes oo, ov, vo, and vv blocks of the ks matrix
!> \param ks ...
!> \param p ...
!> \param t ...
!> \param v ...
!> \param q_index_down ...
!> \param p_index_up ...
!> \param q_index_up ...
!> \param pp ...
!> \param qq ...
!> \param qp ...
!> \param pq ...
!> \param tensor_type ...
!> \param use_virt_orbs ...
!> \param eps_filter ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE assemble_ks_qp_blocks(ks,p,t,v,q_index_down,&
    p_index_up,q_index_up,pp,qq,qp,pq,tensor_type,use_virt_orbs,eps_filter)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: ks, p, t, v, q_index_down, &
                                                p_index_up, q_index_up
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: pp, qq, qp, pq
    INTEGER, INTENT(IN)                      :: tensor_type
    LOGICAL, INTENT(IN)                      :: use_virt_orbs
    REAL(KIND=dp), INTENT(IN)                :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'assemble_ks_qp_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: library_fixed
    TYPE(cp_dbcsr_type)                      :: kst, ksv, no, on, oo, &
                                                q_index_up_nosym, sp, spf, &
                                                t_or, v_or

    CALL timeset(routineN,handle)

    IF (use_virt_orbs) THEN

       ! orthogonalize the orbitals
       CALL cp_dbcsr_init(t_or)
       CALL cp_dbcsr_create(t_or,template=t)
       CALL cp_dbcsr_init(v_or)
       CALL cp_dbcsr_create(v_or,template=v)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,t,p_index_up,&
               0.0_dp,t_or,filter_eps=eps_filter)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,v,q_index_up,&
               0.0_dp,v_or,filter_eps=eps_filter)
      
       ! KS.T
       CALL cp_dbcsr_init(kst)
       CALL cp_dbcsr_create(kst,template=t)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,ks,t_or,&
               0.0_dp,kst,filter_eps=eps_filter)
       ! pp=tr(T)*KS.T
       CALL cp_dbcsr_multiply("T","N",1.0_dp,t_or,kst,&
               0.0_dp,pp,filter_eps=eps_filter)
       ! qp=tr(V)*KS.T
       CALL cp_dbcsr_multiply("T","N",1.0_dp,v_or,kst,&
               0.0_dp,qp,filter_eps=eps_filter)
       CALL cp_dbcsr_release(kst)
       
       ! KS.V
       CALL cp_dbcsr_init(ksv)
       CALL cp_dbcsr_create(ksv,template=v)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,ks,v_or,&
               0.0_dp,ksv,filter_eps=eps_filter)
       ! tr(T)*KS.V
       CALL cp_dbcsr_multiply("T","N",1.0_dp,t_or,ksv,&
               0.0_dp,pq,filter_eps=eps_filter)
       ! tr(V)*KS.V
       CALL cp_dbcsr_multiply("T","N",1.0_dp,v_or,ksv,&
               0.0_dp,qq,filter_eps=eps_filter)
       CALL cp_dbcsr_release(ksv)
       
       CALL cp_dbcsr_release(t_or)
       CALL cp_dbcsr_release(v_or)

    ELSE ! no virtuals, use projected AOs

! THIS PROCEDURE HAS NOT BEEN UPDATED FOR CHOLESKY p/q_index_up/down
       CALL cp_dbcsr_init(sp)
       CALL cp_dbcsr_create(sp,template=q_index_down,&
             matrix_type=dbcsr_type_no_symmetry)
       CALL cp_dbcsr_init(spf)
       CALL cp_dbcsr_create(spf,template=q_index_down,&
             matrix_type=dbcsr_type_no_symmetry)
   
       ! qp=KS*T
       CALL cp_dbcsr_multiply("N","N",1.0_dp,ks,t,0.0_dp,qp,&
                                    filter_eps=eps_filter)
       ! pp=tr(T)*KS.T
       CALL cp_dbcsr_multiply("T","N",1.0_dp,t,qp,0.0_dp,pp,&
                                    filter_eps=eps_filter)
       ! sp=-S_*P
       CALL cp_dbcsr_multiply("N","N",-1.0_dp,q_index_down,p,0.0_dp,sp,&
                                    filter_eps=eps_filter)

       ! sp=1/S^-S_.P
       SELECT CASE (tensor_type)
       CASE (tensor_up_down)
         CALL cp_dbcsr_add_on_diag(sp,1.0_dp)
       CASE (tensor_orthogonal)
         CALL cp_dbcsr_init(q_index_up_nosym)
         CALL cp_dbcsr_create(q_index_up_nosym,template=q_index_up,&
                matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_desymmetrize(q_index_up,q_index_up_nosym)
         CALL cp_dbcsr_add(sp,q_index_up_nosym,1.0_dp,1.0_dp)
         CALL cp_dbcsr_release(q_index_up_nosym)
       END SELECT

       ! spf=(1/S^-S_.P)*KS
       CALL cp_dbcsr_multiply("N","N",1.0_dp,sp,ks,0.0_dp,spf,&
                                    filter_eps=eps_filter)
       
       ! qp=spf*T
       CALL cp_dbcsr_multiply("N","N",1.0_dp,spf,t,0.0_dp,qp,&
                                    filter_eps=eps_filter)

       SELECT CASE (tensor_type)
       CASE (tensor_up_down)
         ! pq=tr(qp)
         CALL cp_dbcsr_transposed(pq,qp,transpose_distribution=.FALSE.)
       CASE (tensor_orthogonal)
         ! pq=sig^.tr(qp)
         CALL cp_dbcsr_multiply("N","T",1.0_dp,p_index_up,qp,0.0_dp,pq,&
                                  filter_eps=eps_filter)
         library_fixed=.FALSE.
         IF (library_fixed) THEN
            CALL cp_dbcsr_transposed(qp,pq,transpose_distribution=.FALSE.)
         ELSE
            CALL cp_dbcsr_init(no)
            CALL cp_dbcsr_create(no,template=qp,&
             matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_multiply("N","N",1.0_dp,qp,p_index_up,0.0_dp,no,&
                                     filter_eps=eps_filter)
            CALL cp_dbcsr_copy(qp,no)
            CALL cp_dbcsr_release(no)
         ENDIF
       END SELECT

       ! qq=spf*tr(sp)
       CALL cp_dbcsr_multiply("N","T",1.0_dp,spf,sp,0.0_dp,qq,&
                                    filter_eps=eps_filter)

       SELECT CASE (tensor_type)
       CASE (tensor_up_down)
  
         CALL cp_dbcsr_init(oo)
         CALL cp_dbcsr_create(oo,template=pp,&
             matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_init(no)
         CALL cp_dbcsr_create(no,template=qp,&
             matrix_type=dbcsr_type_no_symmetry)

         ! first index up 
         CALL cp_dbcsr_multiply("N","N",1.0_dp,q_index_up,qq,0.0_dp,spf,&
                                    filter_eps=eps_filter)
         CALL cp_dbcsr_copy(qq,spf)
         CALL cp_dbcsr_multiply("N","N",1.0_dp,q_index_up,qp,0.0_dp,no,&
                                    filter_eps=eps_filter)
         CALL cp_dbcsr_copy(qp,no)
         CALL cp_dbcsr_multiply("N","N",1.0_dp,p_index_up,pp,0.0_dp,oo,&
                                    filter_eps=eps_filter)
         CALL cp_dbcsr_copy(pp,oo)
         CALL cp_dbcsr_multiply("N","N",1.0_dp,p_index_up,pq,0.0_dp,on,&
                                    filter_eps=eps_filter)
         CALL cp_dbcsr_copy(pq,on)

         CALL cp_dbcsr_release(no)
         CALL cp_dbcsr_release(oo)
   
       CASE (tensor_orthogonal)
         
         CALL cp_dbcsr_init(oo)
         CALL cp_dbcsr_create(oo,template=pp,&
             matrix_type=dbcsr_type_no_symmetry)
         
         ! both indeces up in the pp block
         CALL cp_dbcsr_multiply("N","N",1.0_dp,p_index_up,pp,0.0_dp,oo,&
                                    filter_eps=eps_filter)
         CALL cp_dbcsr_multiply("N","N",1.0_dp,oo,p_index_up,0.0_dp,pp,&
                                    filter_eps=eps_filter)
         
         CALL cp_dbcsr_release(oo)

       END SELECT

       CALL cp_dbcsr_release(sp)
       CALL cp_dbcsr_release(spf)

    ENDIF

    CALL timestop(handle)

  END SUBROUTINE assemble_ks_qp_blocks

! *****************************************************************************
!> \brief Solves the generalized Riccati or Sylvester eqation
!>        using the preconditioned conjugate gradient algorithm
!>          qp + qq.x.oo - vv.x.pp - vv.x.pq.x.oo = 0 [oo and vv are optional]
!>          qp + qq.x - x.pp - x.pq.x = 0
!> \param pp ...
!> \param qq ...
!> \param qp ...
!> \param pq ...
!> \param oo ...
!> \param vv ...
!> \param x ...
!> \param res ...
!> \param neglect_quadratic_term ...
!> \param conjugator ...
!> \param max_iter ...
!> \param eps_convergence ...
!> \param eps_filter ...
!> \param converged ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!>       2011.11 generalized [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
   RECURSIVE SUBROUTINE solve_riccati_equation(pp,qq,qp,pq,oo,vv,x,res,&
               neglect_quadratic_term,&
               conjugator,max_iter,eps_convergence,eps_filter,&
               converged)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: pp, qq
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: qp
    TYPE(cp_dbcsr_type), INTENT(IN)          :: pq
    TYPE(cp_dbcsr_type), INTENT(IN), &
      OPTIONAL                               :: oo, vv
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: x
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: res
    LOGICAL, INTENT(IN)                      :: neglect_quadratic_term
    INTEGER, INTENT(IN)                      :: conjugator, max_iter
    REAL(KIND=dp), INTENT(IN)                :: eps_convergence, eps_filter
    LOGICAL, INTENT(OUT)                     :: converged

    CHARACTER(len=*), PARAMETER :: routineN = 'solve_riccati_equation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istep, iteration, &
                                                nsteps, unit_nr, &
                                                update_prec_freq
    LOGICAL                                  :: prepare_to_exit, present_oo, &
                                                present_vv, quadratic_term, &
                                                restart_conjugator
    REAL(KIND=dp)                            :: best_norm, best_step_size, &
                                                beta, c0, c1, c2, c3, denom, &
                                                kappa, numer, obj_function, &
                                                t1, t2, tau
    REAL(KIND=dp), DIMENSION(3)              :: step_size
    TYPE(cp_dbcsr_type)                      :: aux1, aux2, grad, m, n, oo1, &
                                                oo2, prec, res_trial, step, &
                                                step_oo, vv_step
    TYPE(cp_logger_type), POINTER            :: logger

!TYPE(cp_dbcsr_type)                      :: qqqq, pppp, zero_pq, zero_qp

    CALL timeset(routineN,handle)

    logger => cp_get_default_logger()
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    t1 = m_walltime()
   
!IF (level.gt.5) THEN
!  CPErrorMessage(cp_failure_level,routineP,"recursion level is too high")
!  CPPrecondition(.FALSE.,cp_failure_level,routineP,failure)
!ENDIF
!IF (unit_nr>0) THEN
!   WRITE(unit_nr,*) &
!      "========== LEVEL ",level,"=========="
!ENDIF
!CALL cp_dbcsr_print(qq)
!CALL cp_dbcsr_print(pp)
!CALL cp_dbcsr_print(qp)
!!CALL cp_dbcsr_print(pq)
!IF (unit_nr>0) THEN
!   WRITE(unit_nr,*) &
!      "====== END LEVEL ",level,"=========="
!ENDIF

    quadratic_term=.NOT.neglect_quadratic_term
    present_oo=PRESENT(oo)
    present_vv=PRESENT(vv)

    ! create aux1 matrix and init
    CALL cp_dbcsr_init(aux1)
    CALL cp_dbcsr_create(aux1,template=pp)
    CALL cp_dbcsr_copy(aux1,pp)
    CALL cp_dbcsr_scale(aux1,-1.0_dp)

    ! create aux2 matrix and init
    CALL cp_dbcsr_init(aux2)
    CALL cp_dbcsr_create(aux2,template=qq)
    CALL cp_dbcsr_copy(aux2,qq)

    ! create the gradient matrix and init
    CALL cp_dbcsr_init(grad)
    CALL cp_dbcsr_create(grad,template=x)
    CALL cp_dbcsr_set(grad,0.0_dp)

    ! create a preconditioner
    ! RZK-warning how to apply it to up_down tensor?
    CALL cp_dbcsr_init(prec)
    CALL cp_dbcsr_create(prec,template=x)
    !CALL create_preconditioner(prec,aux1,aux2,qp,res,tensor_type,eps_filter)
    !CALL cp_dbcsr_set(prec,1.0_dp)

    ! create the step matrix and init
    CALL cp_dbcsr_init(step)
    CALL cp_dbcsr_create(step,template=x)
    !CALL cp_dbcsr_hadamard_product(prec,grad,step)
    !CALL cp_dbcsr_scale(step,-1.0_dp)

    CALL cp_dbcsr_init(n)
    CALL cp_dbcsr_create(n,template=x)
    CALL cp_dbcsr_init(m)
    CALL cp_dbcsr_create(m,template=x)
    CALL cp_dbcsr_init(oo1)
    CALL cp_dbcsr_create(oo1,template=pp)
    CALL cp_dbcsr_init(oo2)
    CALL cp_dbcsr_create(oo2,template=pp)
    CALL cp_dbcsr_init(res_trial)
    CALL cp_dbcsr_create(res_trial,template=res)
    CALL cp_dbcsr_init(vv_step)
    CALL cp_dbcsr_create(vv_step,template=res)
    CALL cp_dbcsr_init(step_oo)
    CALL cp_dbcsr_create(step_oo,template=res)
    
    ! start conjugate gradient iterations
    iteration=0
    converged=.FALSE.
    prepare_to_exit=.FALSE.
    beta=0.0_dp
    best_step_size=0.0_dp
    best_norm=1.0E+100_dp
    !ecorr=0.0_dp
    !change_ecorr=0.0_dp
    restart_conjugator=.FALSE.
    update_prec_freq=20
    DO

       ! (re)-compute the residuals
       IF (iteration.eq.0) THEN
          CALL cp_dbcsr_copy(res,qp)
          IF (present_oo) THEN
             CALL cp_dbcsr_multiply("N","N",+1.0_dp,qq,x,0.0_dp,res_trial,&
                     filter_eps=eps_filter)
             CALL cp_dbcsr_multiply("N","N",+1.0_dp,res_trial,oo,1.0_dp,res,&
                     filter_eps=eps_filter)
          ELSE
             CALL cp_dbcsr_multiply("N","N",+1.0_dp,qq,x,1.0_dp,res,&
                  filter_eps=eps_filter)
          ENDIF
          IF (present_vv) THEN
             CALL cp_dbcsr_multiply("N","N",-1.0_dp,x,pp,0.0_dp,res_trial,&
                     filter_eps=eps_filter)
             CALL cp_dbcsr_multiply("N","N",+1.0_dp,vv,res_trial,1.0_dp,res,&
                     filter_eps=eps_filter)
          ELSE
             CALL cp_dbcsr_multiply("N","N",-1.0_dp,x,pp,1.0_dp,res,&
                     filter_eps=eps_filter)
          ENDIF
          IF (quadratic_term) THEN
             IF (present_oo) THEN
                CALL cp_dbcsr_multiply("N","N",+1.0_dp,pq,x,0.0_dp,oo1,&
                        filter_eps=eps_filter)
                CALL cp_dbcsr_multiply("N","N",+1.0_dp,oo1,oo,0.0_dp,oo2,&
                        filter_eps=eps_filter)
             ELSE
                CALL cp_dbcsr_multiply("N","N",+1.0_dp,pq,x,0.0_dp,oo2,&
                        filter_eps=eps_filter)
             ENDIF
             IF (present_vv) THEN
                CALL cp_dbcsr_multiply("N","N",-1.0_dp,x,oo2,0.0_dp,res_trial,&
                        filter_eps=eps_filter)
                CALL cp_dbcsr_multiply("N","N",+1.0_dp,vv,res_trial,1.0_dp,res,&
                        filter_eps=eps_filter)
             ELSE
                CALL cp_dbcsr_multiply("N","N",-1.0_dp,x,oo2,1.0_dp,res,&
                     filter_eps=eps_filter)
             ENDIF
          ENDIF
          CALL cp_dbcsr_norm(res, dbcsr_norm_maxabsnorm, norm_scalar=best_norm)
       ELSE
          CALL cp_dbcsr_add(res,m,1.0_dp,best_step_size)
          CALL cp_dbcsr_add(res,n,1.0_dp,-best_step_size*best_step_size)
          CALL cp_dbcsr_filter(res,eps_filter)
       ENDIF

       ! check convergence and other exit criteria
       converged=(best_norm.lt.eps_convergence)
       IF (converged.OR.(iteration.ge.max_iter)) THEN
          prepare_to_exit=.TRUE.
       ENDIF

       IF (.NOT.prepare_to_exit) THEN
          
          ! update aux1=-pp-pq.x.oo and aux2=qq-vv.x.pq
          IF (quadratic_term) THEN
             IF (iteration.eq.0) THEN
                IF (present_oo) THEN
                   CALL cp_dbcsr_multiply("N","N",-1.0_dp,pq,x,0.0_dp,oo1,&
                           filter_eps=eps_filter)
                   CALL cp_dbcsr_multiply("N","N",+1.0_dp,oo1,oo,1.0_dp,aux1,&
                           filter_eps=eps_filter)
                ELSE
                   CALL cp_dbcsr_multiply("N","N",-1.0_dp,pq,x,1.0_dp,aux1,&
                           filter_eps=eps_filter)
                ENDIF
                IF (present_vv) THEN
                   CALL cp_dbcsr_multiply("N","N",-1.0_dp,vv,x,0.0_dp,res_trial,&
                           filter_eps=eps_filter)
                   CALL cp_dbcsr_multiply("N","N",+1.0_dp,res_trial,pq,1.0_dp,aux2,&
                           filter_eps=eps_filter)
                ELSE
                   CALL cp_dbcsr_multiply("N","N",-1.0_dp,x,pq,1.0_dp,aux2,&
                           filter_eps=eps_filter)
                ENDIF
             ELSE
                IF (present_oo) THEN
                   CALL cp_dbcsr_multiply("N","N",-best_step_size,pq,step_oo,1.0_dp,aux1,&
                           filter_eps=eps_filter)
                ELSE
                   CALL cp_dbcsr_multiply("N","N",-best_step_size,pq,step,1.0_dp,aux1,&
                           filter_eps=eps_filter)
                ENDIF
                IF (present_vv) THEN
                   CALL cp_dbcsr_multiply("N","N",-best_step_size,vv_step,pq,1.0_dp,aux2,&
                           filter_eps=eps_filter)
                ELSE
                   CALL cp_dbcsr_multiply("N","N",-best_step_size,step,pq,1.0_dp,aux2,&
                           filter_eps=eps_filter)
                ENDIF
             ENDIF
          ENDIF

          ! recompute the gradient, do not update it yet
          ! use m matrix as a temporary storage
          ! grad=t(vv).res.t(aux1)+t(aux2).res.t(oo)
          IF (present_vv) THEN
             CALL cp_dbcsr_multiply("N","T",1.0_dp,res,aux1,0.0_dp,res_trial,&
                     filter_eps=eps_filter)
             CALL cp_dbcsr_multiply("T","N",1.0_dp,vv,res_trial,0.0_dp,m,&
                     filter_eps=eps_filter)
          ELSE
             CALL cp_dbcsr_multiply("N","T",1.0_dp,res,aux1,0.0_dp,m,&
                     filter_eps=eps_filter)
          ENDIF
          IF (present_oo) THEN
             CALL cp_dbcsr_multiply("T","N",1.0_dp,aux1,res,0.0_dp,res_trial,&
                     filter_eps=eps_filter)
             CALL cp_dbcsr_multiply("N","T",1.0_dp,res_trial,oo,1.0_dp,m,&
                     filter_eps=eps_filter)
          ELSE
             CALL cp_dbcsr_multiply("T","N",1.0_dp,aux2,res,1.0_dp,m,&
                     filter_eps=eps_filter)
          ENDIF

          ! compute preconditioner
          !IF (iteration.eq.0.OR.(mod(iteration,update_prec_freq).eq.0)) THEN
          IF (iteration.eq.0) THEN
             CALL create_preconditioner(prec,aux1,aux2,eps_filter)
             !restart_conjugator=.TRUE.
!CALL cp_dbcsr_set(prec,1.0_dp)
!CALL cp_dbcsr_print(prec)
          ENDIF

          ! compute the conjugation coefficient - beta
          IF ((iteration.eq.0).OR.restart_conjugator) THEN
             beta=0.0_dp
          ELSE
             restart_conjugator=.FALSE.
             SELECT CASE (conjugator)
             CASE (cg_hestenes_stiefel)
                CALL cp_dbcsr_add(grad,m,-1.0_dp,1.0_dp)
                CALL cp_dbcsr_hadamard_product(prec,grad,n)
                CALL cp_dbcsr_trace(n,m,numer,"T","N")
                CALL cp_dbcsr_trace(grad,step,denom,"T","N")
                beta=numer/denom
             CASE (cg_fletcher_reeves)
                CALL cp_dbcsr_hadamard_product(prec,grad,n)
                CALL cp_dbcsr_trace(grad,n,denom,"T","N")
                CALL cp_dbcsr_hadamard_product(prec,m,n)
                CALL cp_dbcsr_trace(m,n,numer,"T","N")
                beta=numer/denom
             CASE (cg_polak_ribiere)
                CALL cp_dbcsr_hadamard_product(prec,grad,n)
                CALL cp_dbcsr_trace(grad,n,denom,"T","N")
                CALL cp_dbcsr_add(grad,m,-1.0_dp,1.0_dp)
                CALL cp_dbcsr_hadamard_product(prec,grad,n)
                CALL cp_dbcsr_trace(n,m,numer,"T","N")
                beta=numer/denom
             CASE (cg_fletcher)
                CALL cp_dbcsr_hadamard_product(prec,m,n)
                CALL cp_dbcsr_trace(m,n,numer,"T","N")
                CALL cp_dbcsr_trace(grad,step,denom,"T","N")
                beta=-1.0_dp*numer/denom
             CASE (cg_liu_storey)
                CALL cp_dbcsr_trace(grad,step,denom,"T","N")
                CALL cp_dbcsr_add(grad,m,-1.0_dp,1.0_dp)
                CALL cp_dbcsr_hadamard_product(prec,grad,n)
                CALL cp_dbcsr_trace(n,m,numer,"T","N")
                beta=-1.0_dp*numer/denom
             CASE (cg_dai_yuan)
                CALL cp_dbcsr_hadamard_product(prec,m,n)
                CALL cp_dbcsr_trace(m,n,numer,"T","N")
                CALL cp_dbcsr_add(grad,m,-1.0_dp,1.0_dp)
                CALL cp_dbcsr_trace(grad,step,denom,"T","N")
                beta=numer/denom
             CASE (cg_hager_zhang)
                CALL cp_dbcsr_add(grad,m,-1.0_dp,1.0_dp)
                CALL cp_dbcsr_trace(grad,step,denom,"T","N")
                CALL cp_dbcsr_hadamard_product(prec,grad,n)
                CALL cp_dbcsr_trace(n,grad,numer,"T","N")
                kappa=2.0_dp*numer/denom
                CALL cp_dbcsr_trace(n,m,numer,"T","N")
                tau=numer/denom
                CALL cp_dbcsr_trace(step,m,numer,"T","N")
                beta=tau-kappa*numer/denom
             CASE (cg_zero)
                beta=0.0_dp
             CASE DEFAULT
                CPABORT("illegal conjugator")
             END SELECT
          ENDIF ! iteration.eq.0
          
          ! move the current gradient to its storage
          CALL cp_dbcsr_copy(grad,m)
          
          ! precondition new gradient (use m as tmp storage)
          CALL cp_dbcsr_hadamard_product(prec,grad,m)
          CALL cp_dbcsr_filter(m,eps_filter)
          
          ! recompute the step direction
          CALL cp_dbcsr_add(step,m,beta,-1.0_dp)
          CALL cp_dbcsr_filter(step,eps_filter)

!! ALTERNATIVE METHOD TO OBTAIN THE STEP FROM THE GRADIENT
!CALL cp_dbcsr_init(qqqq)
!CALL cp_dbcsr_create(qqqq,template=qq)
!CALL cp_dbcsr_init(pppp)
!CALL cp_dbcsr_create(pppp,template=pp)
!CALL cp_dbcsr_init(zero_pq)
!CALL cp_dbcsr_create(zero_pq,template=pq)
!CALL cp_dbcsr_init(zero_qp)
!CALL cp_dbcsr_create(zero_qp,template=qp)
!CALL cp_dbcsr_multiply("T","N",1.0_dp,aux2,aux2,0.0_dp,qqqq,&
!        filter_eps=eps_filter)
!CALL cp_dbcsr_multiply("N","T",-1.0_dp,aux1,aux1,0.0_dp,pppp,&
!        filter_eps=eps_filter)
!CALL cp_dbcsr_set(zero_qp,0.0_dp)
!CALL cp_dbcsr_set(zero_pq,0.0_dp)
!CALL solve_riccati_equation(pppp,qqqq,grad,zero_pq,zero_qp,zero_qp,&
!               .TRUE.,tensor_type,&
!               conjugator,max_iter,eps_convergence,eps_filter,&
!               converged,level+1)
!CALL cp_dbcsr_release(qqqq)
!CALL cp_dbcsr_release(pppp)
!CALL cp_dbcsr_release(zero_qp)
!CALL cp_dbcsr_release(zero_pq)
    
          ! calculate the optimal step size
          ! m=step.aux1+aux2.step
          IF (present_vv) THEN
             CALL cp_dbcsr_multiply("N","N",1.0_dp,vv,step,0.0_dp,vv_step,&
                     filter_eps=eps_filter)
             CALL cp_dbcsr_multiply("N","N",1.0_dp,vv_step,aux1,0.0_dp,m,&
                     filter_eps=eps_filter)
          ELSE
             CALL cp_dbcsr_multiply("N","N",1.0_dp,step,aux1,0.0_dp,m,&
                     filter_eps=eps_filter)
          ENDIF
          IF (present_oo) THEN
             CALL cp_dbcsr_multiply("N","N",1.0_dp,step,oo,0.0_dp,step_oo,&
                     filter_eps=eps_filter)
             CALL cp_dbcsr_multiply("N","N",1.0_dp,aux2,step_oo,1.0_dp,m,&
                     filter_eps=eps_filter)
          ELSE
             CALL cp_dbcsr_multiply("N","N",1.0_dp,aux2,step,1.0_dp,m,&
                     filter_eps=eps_filter)
          ENDIF

          IF (quadratic_term) THEN
             ! n=step.pq.step
             IF (present_oo) THEN
                CALL cp_dbcsr_multiply("N","N",1.0_dp,pq,step,0.0_dp,oo1,&
                        filter_eps=eps_filter)
                CALL cp_dbcsr_multiply("N","N",1.0_dp,oo1,oo,0.0_dp,oo2,&
                        filter_eps=eps_filter)
             ELSE
                CALL cp_dbcsr_multiply("N","N",1.0_dp,pq,step,0.0_dp,oo2,&
                        filter_eps=eps_filter)
             ENDIF
             IF (present_vv) THEN
                CALL cp_dbcsr_multiply("N","N",1.0_dp,step,oo2,0.0_dp,res_trial,&
                        filter_eps=eps_filter)
                CALL cp_dbcsr_multiply("N","N",1.0_dp,vv,res_trial,0.0_dp,n,&
                        filter_eps=eps_filter)
             ELSE
                CALL cp_dbcsr_multiply("N","N",1.0_dp,step,oo2,0.0_dp,n,&
                        filter_eps=eps_filter)
             ENDIF

          ELSE
             CALL cp_dbcsr_set(n,0.0_dp)
          ENDIF

          ! calculate coefficients of the cubic eq for alpha - step size
          c0 = 2.0_dp*(cp_dbcsr_frobenius_norm(n))**2
          
          CALL cp_dbcsr_trace(m,n,c1,"T","N")
          c1 = -3.0_dp*c1

          CALL cp_dbcsr_trace(res,n,c2,"T","N")
          c2 = -2.0_dp*c2+(cp_dbcsr_frobenius_norm(m))**2
          
          CALL cp_dbcsr_trace(res,m,c3,"T","N")

          ! find step size
          CALL analytic_line_search(c0,c1,c2,c3,step_size,nsteps)
          
          IF (nsteps.eq.0) THEN
             CPABORT("no step sizes!")
          ENDIF
          ! if we have several possible step sizes
          ! choose one with the lowest objective function
          best_norm = 1.0E+100_dp 
          best_step_size = 0.0_dp
          DO istep = 1,nsteps
             ! recompute the residues
             CALL cp_dbcsr_copy(res_trial,res)
             CALL cp_dbcsr_add(res_trial,m,1.0_dp,step_size(istep))
             CALL cp_dbcsr_add(res_trial,n,1.0_dp,-step_size(istep)*step_size(istep))
             CALL cp_dbcsr_filter(res_trial,eps_filter)
             ! RZK-warning objective function might be different in the case of
             ! tensor_up_down
             !obj_function=0.5_dp*(cp_dbcsr_frobenius_norm(res_trial))**2
             CALL cp_dbcsr_norm(res_trial, dbcsr_norm_maxabsnorm, norm_scalar=obj_function)
             IF (obj_function.lt.best_norm) THEN
                best_norm=obj_function
                best_step_size=step_size(istep)
             ENDIF
          ENDDO

       ENDIF
       
       ! update X along the line 
       CALL cp_dbcsr_add(x,step,1.0_dp,best_step_size)
       CALL cp_dbcsr_filter(x,eps_filter)
       
       ! evaluate current energy correction
       !change_ecorr=ecorr
       !CALL cp_dbcsr_trace(qp,x,ecorr,"T","N")
       !change_ecorr=ecorr-change_ecorr

       ! check convergence and other exit criteria
       converged=(best_norm.lt.eps_convergence)
       IF (converged.OR.(iteration.ge.max_iter)) THEN
          prepare_to_exit=.TRUE.
       ENDIF

       t2 = m_walltime()
       
       IF (unit_nr>0) THEN
          WRITE(unit_nr,'(T6,A,1X,I4,1X,E12.3,F8.3)') &
             "RICCATI iter ",iteration,best_norm,t2-t1
          !WRITE(unit_nr,'(T6,A,1X,I4,1X,F15.9,F15.9,E12.3,F8.3)') &
          !   "RICCATI iter ",iteration,ecorr,change_ecorr,best_norm,t2-t1
       ENDIF

       t1 = m_walltime()

       iteration=iteration+1

       IF (prepare_to_exit) EXIT
    
    ENDDO

    CALL cp_dbcsr_release(aux1)
    CALL cp_dbcsr_release(aux2)
    CALL cp_dbcsr_release(grad)
    CALL cp_dbcsr_release(step)
    CALL cp_dbcsr_release(n)
    CALL cp_dbcsr_release(m)
    CALL cp_dbcsr_release(oo1)
    CALL cp_dbcsr_release(oo2)
    CALL cp_dbcsr_release(res_trial)
    CALL cp_dbcsr_release(vv_step)
    CALL cp_dbcsr_release(step_oo)
    
    CALL timestop(handle)

  END SUBROUTINE solve_riccati_equation

! *****************************************************************************
!> \brief Computes a preconditioner from diagonal elements of ~f_oo, ~f_vv
!>        The preconditioner is approximately equal to
!>        prec_ai ~ (e_a - e_i)^(-2)
!>        However, the real expression is more complex
!> \param prec ...
!> \param pp ...
!> \param qq ...
!> \param eps_filter ...
!> \par History
!>       2011.07 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE create_preconditioner(prec,pp,qq,eps_filter)

    TYPE(cp_dbcsr_type), INTENT(OUT)         :: prec
    TYPE(cp_dbcsr_type), INTENT(IN)          :: pp, qq
    REAL(KIND=dp), INTENT(IN)                :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'create_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER :: col, handle, hold, iblock_col, iblock_row, mynode, &
      nblkcols_tot, nblkrows_tot, p_nrows, q_nrows, row
    LOGICAL                                  :: tr
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: p_diagonal, q_diagonal
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block
    TYPE(cp_dbcsr_type)                      :: pp_diag, qq_diag, t1, t2, tmp

!LOGICAL, INTENT(IN)                      :: use_virt_orbs

    CALL timeset(routineN,handle)

!    ! copy diagonal elements
!    CALL cp_dbcsr_get_info(pp,nfullrows_total=nrows)
!    CALL cp_dbcsr_init(pp_diag)
!    CALL cp_dbcsr_create(pp_diag,template=pp)
!    ALLOCATE(diagonal(nrows))
!    CALL cp_dbcsr_get_diag(pp,diagonal)
!    CALL cp_dbcsr_add_on_diag(pp_diag,1.0_dp)
!    CALL cp_dbcsr_set_diag(pp_diag,diagonal)
!    DEALLOCATE(diagonal)
!
    ! initialize a matrix to 1.0
    CALL cp_dbcsr_init(tmp)
    CALL cp_dbcsr_create(tmp,template=prec)
    ! use an ugly hack to set all elements of tmp to 1
    ! because cp_dbcsr_set does not do it (despite its name)
    !CALL cp_dbcsr_set(tmp,1.0_dp)
    mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(cp_dbcsr_distribution(tmp)))
    CALL cp_dbcsr_work_create(tmp, work_mutable=.TRUE.)
    nblkrows_tot = cp_dbcsr_nblkrows_total(tmp)
    nblkcols_tot = cp_dbcsr_nblkcols_total(tmp)
    DO row = 1, nblkrows_tot
       DO col = 1, nblkcols_tot
          tr = .FALSE.
          iblock_row = row
          iblock_col = col
          CALL cp_dbcsr_get_stored_coordinates(tmp,iblock_row,iblock_col,tr,hold)
          IF (hold.EQ.mynode) THEN
             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(tmp,iblock_row,iblock_col,p_new_block)
             CPASSERT(ASSOCIATED(p_new_block))
             p_new_block(:,:) = 1.0_dp
          ENDIF ! mynode
       ENDDO
    ENDDO
    CALL cp_dbcsr_finalize(tmp)

    ! copy diagonal elements of pp into cols of a matrix
    CALL cp_dbcsr_get_info(pp,nfullrows_total=p_nrows)
    CALL cp_dbcsr_init(pp_diag)
    CALL cp_dbcsr_create(pp_diag,template=pp)
    ALLOCATE(p_diagonal(p_nrows))
    CALL cp_dbcsr_get_diag(pp,p_diagonal)
    CALL cp_dbcsr_add_on_diag(pp_diag,1.0_dp)
    CALL cp_dbcsr_set_diag(pp_diag,p_diagonal)
    ! RZK-warning is it possible to use cp_dbcsr_scale_by_vector?
    ! or even insert elements directly in the prev cycles
    CALL cp_dbcsr_init(t2)
    CALL cp_dbcsr_create(t2,template=prec)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,tmp,pp_diag,&
             0.0_dp,t2,filter_eps=eps_filter)

    ! copy diagonal elements qq into rows of a matrix
    CALL cp_dbcsr_get_info(qq, nfullrows_total=q_nrows )
    CALL cp_dbcsr_init(qq_diag)
    CALL cp_dbcsr_create(qq_diag,template=qq)
    ALLOCATE(q_diagonal(q_nrows))
    CALL cp_dbcsr_get_diag(qq,q_diagonal)
    CALL cp_dbcsr_add_on_diag(qq_diag,1.0_dp)
    CALL cp_dbcsr_set_diag(qq_diag,q_diagonal)
    CALL cp_dbcsr_set(tmp,1.0_dp)
    CALL cp_dbcsr_init(t1)
    CALL cp_dbcsr_create(t1,template=prec)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,qq_diag,tmp,&
             0.0_dp,t1,filter_eps=eps_filter)

    CALL cp_dbcsr_hadamard_product(t1,t2,prec)
    CALL cp_dbcsr_release(t1)
    CALL cp_dbcsr_scale(prec,2.0_dp)

    ! Get the diagonal of tr(qq).qq 
    CALL cp_dbcsr_multiply("T","N",1.0_dp,qq,qq,&
             0.0_dp,qq_diag,retain_sparsity=.TRUE.,&
             filter_eps=eps_filter)
    CALL cp_dbcsr_get_diag(qq_diag,q_diagonal)
    CALL cp_dbcsr_set(qq_diag,0.0_dp)
    CALL cp_dbcsr_add_on_diag(qq_diag,1.0_dp)
    CALL cp_dbcsr_set_diag(qq_diag,q_diagonal)
    DEALLOCATE(q_diagonal)
    CALL cp_dbcsr_set(tmp,1.0_dp)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,qq_diag,tmp,&
             0.0_dp,t2,filter_eps=eps_filter)
    CALL cp_dbcsr_release(qq_diag)
    CALL cp_dbcsr_add(prec,t2,1.0_dp,1.0_dp)

    ! Get the diagonal of pp.tr(pp)
    CALL cp_dbcsr_multiply("N","T",1.0_dp,pp,pp,&
             0.0_dp,pp_diag,retain_sparsity=.TRUE.,&
             filter_eps=eps_filter)
    CALL cp_dbcsr_get_diag(pp_diag,p_diagonal)
    CALL cp_dbcsr_set(pp_diag,0.0_dp)
    CALL cp_dbcsr_add_on_diag(pp_diag,1.0_dp)
    CALL cp_dbcsr_set_diag(pp_diag,p_diagonal)
    DEALLOCATE(p_diagonal)
    CALL cp_dbcsr_set(tmp,1.0_dp)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,tmp,pp_diag,&
             0.0_dp,t2,filter_eps=eps_filter)
    CALL cp_dbcsr_release(tmp)
    CALL cp_dbcsr_release(pp_diag)
    CALL cp_dbcsr_add(prec,t2,1.0_dp,1.0_dp)

    ! now add the residual component
    !CALL cp_dbcsr_hadamard_product(res,qp,t2)
    !CALL cp_dbcsr_add(prec,t2,1.0_dp,-2.0_dp)
    CALL cp_dbcsr_release(t2)
    CALL cp_dbcsr_function_of_elements(prec,func=dbcsr_func_inverse)
    CALL cp_dbcsr_filter(prec,eps_filter)

    CALL timestop(handle)

  END SUBROUTINE create_preconditioner

! *****************************************************************************
!> \brief Finds real roots of a cubic equation
!>    >        a*x**3 + b*x**2 + c*x + d = 0
!>        and returns only those roots for which the derivative is positive
!>
!>   Step 0: Check the true order of the equation. Cubic, quadratic, linear?
!>   Step 1: Calculate p and q
!>           p = ( 3*c/a - (b/a)**2 ) / 3
!>           q = ( 2*(b/a)**3 - 9*b*c/a/a + 27*d/a ) / 27
!>   Step 2: Calculate discriminant D
!>           D = (p/3)**3 + (q/2)**2
!>   Step 3: Depending on the sign of D, we follow different strategy.
!>           If D<0, three distinct real roots.
!>           If D=0, three real roots of which at least two are equal.
!>           If D>0, one real and two complex roots.
!>   Step 3a: For D>0 and D=0,
!>           Calculate u and v
!>           u = cubic_root(-q/2 + sqrt(D))
!>           v = cubic_root(-q/2 - sqrt(D))
!>           Find the three transformed roots
!>           y1 = u + v
!>           y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
!>           y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
!>   Step 3b Alternately, for D<0, a trigonometric formulation is more convenient
!>           y1 =  2 * sqrt(|p|/3) * cos(phi/3)
!>           y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
!>           y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
!>           where phi = acos(-q/2/sqrt(|p|**3/27))
!>                 pi  = 3.141592654...
!>   Step 4  Find the real roots
!>           x = y - b/a/3
!>   Step 5  Check the derivative and return only those real roots
!>           for which the derivative is positive
!>
!> \param a ...
!> \param b ...
!> \param c ...
!> \param d ...
!> \param minima ...
!> \param nmins ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE analytic_line_search(a,b,c,d,minima,nmins)
  
    REAL(KIND=dp), INTENT(IN)                :: a, b, c, d
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: minima
    INTEGER, INTENT(OUT)                     :: nmins

    CHARACTER(len=*), PARAMETER :: routineN = 'analytic_line_search', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nroots
    REAL(KIND=dp)                            :: DD, der, p, phi, pi, q, &
                                                temp1, temp2, u, v, y1, y2, &
                                                y2i, y2r, y3
    REAL(KIND=dp), DIMENSION(3)              :: x

!    CALL timeset(routineN,handle)

    pi = ACOS(-1.0_dp)

    ! Step 0: Check coefficients and find the true order of the eq
    IF (a .EQ. 0.0_dp) THEN
     IF (b .EQ. 0.0_dp) THEN
      IF (c .EQ. 0.0_dp) THEN
        ! Non-equation, no valid solutions
        nroots = 0
      ELSE
        ! Linear equation with one root.
        nroots = 1
        x(1) = -d/c
      ENDIF
     ELSE
      ! Quadratic equation with max two roots.
      DD = c*c-4.0_dp*b*d
      IF (DD .GT. 0.0_dp) THEN
        nroots = 2
        x(1) = (-c+SQRT(DD))/2.0_dp/b
        x(2) = (-c-SQRT(DD))/2.0_dp/b
      ELSE IF (DD .LT. 0.0_dp) THEN
        nroots = 0
      ELSE
        nroots=1
        x(1) = -c/2.0_dp/b
      ENDIF
     ENDIF
    ELSE
     ! Cubic equation with max three roots
     ! Calculate p and q
     p  = c/a - b*b/a/a/3.0_dp
     q  = (2.0_dp*b*b*b/a/a/a - 9.0_dp*b*c/a/a + 27.0_dp*d/a) / 27.0_dp

     ! Calculate DD
     DD = p*p*p/27.0_dp + q*q/4.0_dp

     IF (DD .LT. 0.0_dp) THEN
      ! three real unequal roots -- use the trigonometric formulation
      phi = ACOS(-q/2.0_dp/SQRT(ABS(p*p*p)/27.0_dp))
      temp1=2.0_dp*SQRT(ABS(p)/3.0_dp)
      y1 =  temp1*COS(phi/3.0_dp)
      y2 = -temp1*COS((phi+pi)/3.0_dp)
      y3 = -temp1*COS((phi-pi)/3.0_dp)
     ELSE
      ! 1 real & 2 conjugate complex roots OR 3 real roots (some are equal)
      temp1 = -q/2.0_dp + SQRT(DD)
      temp2 = -q/2.0_dp - SQRT(DD)
      u = ABS(temp1)**(1.0_dp/3.0_dp)
      v = ABS(temp2)**(1.0_dp/3.0_dp)
      IF (temp1 .LT. 0.0_dp) u=-u
      IF (temp2 .LT. 0.0_dp) v=-v
      y1  = u + v
      y2r = -(u+v)/2.0_dp
      y2i =  (u-v)*SQRT(3.0_dp)/2.0_dp
     ENDIF

     ! Final transformation 
     temp1 = b/a/3.0_dp
     y1 = y1-temp1
     y2 = y2-temp1
     y3 = y3-temp1
     y2r=y2r-temp1

     ! Assign answers
     IF (DD .LT. 0.0_dp) THEN
      nroots = 3
      x(1) = y1
      x(2) = y2
      x(3) = y3
     ELSE IF (DD .EQ. 0.0_dp) THEN
      nroots = 2
      x(1) = y1
      x(2) = y2r
      !x(3) = cmplx(y2r,  0.)
     ELSE
      nroots = 1
      x(1) = y1
      !x(2) = cmplx(y2r, y2i)
      !x(3) = cmplx(y2r,-y2i)
     ENDIF

    ENDIF

!write(*,'(i2,a)') nroots, ' real root(s)'
    nmins=0
    DO i=1, nroots
      ! maximum or minimum? use the derivative
      ! 3*a*x**2+2*b*x+c
      der = 3.0_dp*a*x(i)*x(i)+2.0_dp*b*x(i)+c
      IF (der .GT. 0.0_dp) THEN
       nmins=nmins+1
       minima(nmins)=x(i)
!write(*,'(a,i2,a,f10.5)') 'Minimum ', i, ', value: ', x(i)
      ENDIF
    ENDDO

!    CALL timestop(handle)

  END SUBROUTINE analytic_line_search

! *****************************************************************************
!> \brief Diagonalizes diagonal blocks of a symmetric dbcsr matrix
!>        and returs its eigenvectors
!> \param matrix ...
!> \param c ...
!> \param e ...
!> \par History
!>       2011.07 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE diagonalize_diagonal_blocks(matrix,c,e)
  
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: c
    TYPE(cp_dbcsr_type), INTENT(OUT), &
      OPTIONAL                               :: e

    CHARACTER(len=*), PARAMETER :: routineN = 'diagonalize_diagonal_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, &
                                                iblock_row, iblock_size, &
                                                info, lwork, orbital
    LOGICAL                                  :: block_needed, do_eigenvalues
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, work
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: data_copy
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: data_p, p_new_block
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN,handle)

    IF (PRESENT(e)) THEN
       do_eigenvalues=.TRUE.
    ELSE
       do_eigenvalues=.FALSE.
    ENDIF

    ! create a matrix for eigenvectors
    CALL cp_dbcsr_work_create(c,work_mutable=.TRUE.)
    IF (do_eigenvalues) &
       CALL cp_dbcsr_work_create(e,work_mutable=.TRUE.)

    CALL cp_dbcsr_iterator_start(iter,matrix)

    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

       CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,row_size=iblock_size)

       block_needed=.FALSE.
       IF (iblock_row==iblock_col) block_needed=.TRUE.

       IF (block_needed) THEN

          ! Prepare data
          ALLOCATE(eigenvalues(iblock_size))
          ALLOCATE(data_copy(iblock_size,iblock_size))
          data_copy(:,:)=data_p(:,:)

          ! Query the optimal workspace for dsyev
          LWORK = -1
          ALLOCATE(WORK(MAX(1,LWORK)))
          CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
          LWORK = INT(WORK( 1 ))
          DEALLOCATE(WORK)

          ! Allocate the workspace and solve the eigenproblem
          ALLOCATE(WORK(MAX(1,LWORK)))
          CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
          IF( INFO.NE.0 ) THEN
             CPABORT("DSYEV failed")
          END IF

          ! copy eigenvectors into a cp_dbcsr matrix
          NULLIFY (p_new_block)
          CALL cp_dbcsr_reserve_block2d(c,iblock_row,iblock_col,p_new_block)
          CPASSERT(ASSOCIATED(p_new_block))
          p_new_block(:,:) = data_copy(:,:)

          ! if requested copy eigenvalues into a cp_dbcsr matrix
          IF (do_eigenvalues) THEN
             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(e,iblock_row,iblock_col,p_new_block)
             CPASSERT(ASSOCIATED(p_new_block))
             p_new_block(:,:) = 0.0_dp
             DO orbital=1,iblock_size
                p_new_block(orbital,orbital)=eigenvalues(orbital)
             ENDDO
          ENDIF

          DEALLOCATE(WORK)
          DEALLOCATE(data_copy)
          DEALLOCATE(eigenvalues)

       ENDIF

    ENDDO

    CALL cp_dbcsr_iterator_stop(iter)
    
    CALL cp_dbcsr_finalize(c)
    IF (do_eigenvalues) CALL cp_dbcsr_finalize(e)

    CALL timestop(handle)

  END SUBROUTINE diagonalize_diagonal_blocks

! *****************************************************************************
!> \brief Transforms a matrix M_out = tr(U1) * M_in * U2
!> \param matrix ...
!> \param u1 ...
!> \param u2 ...
!> \param eps_filter ...
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE matrix_forward_transform(matrix,u1,u2,eps_filter)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_dbcsr_type), INTENT(IN)          :: u1, u2
    REAL(KIND=dp), INTENT(IN)                :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_forward_transform', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_dbcsr_type)                      :: tmp

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_init(tmp)
    CALL cp_dbcsr_create(tmp,template=matrix,&
            matrix_type=dbcsr_type_no_symmetry)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix,u2,0.0_dp,tmp,&
               filter_eps=eps_filter)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,u1,tmp,0.0_dp,matrix,&
               filter_eps=eps_filter)
    CALL cp_dbcsr_release(tmp)

    CALL timestop(handle)

  END SUBROUTINE matrix_forward_transform

! *****************************************************************************
!> \brief Transforms a matrix M_out = U1 * M_in * tr(U2)
!> \param matrix ...
!> \param u1 ...
!> \param u2 ...
!> \param eps_filter ...
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE matrix_backward_transform(matrix,u1,u2,eps_filter)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_dbcsr_type), INTENT(IN)          :: u1, u2
    REAL(KIND=dp), INTENT(IN)                :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_backward_transform', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_dbcsr_type)                      :: tmp

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_init(tmp)
    CALL cp_dbcsr_create(tmp,template=matrix,&
            matrix_type=dbcsr_type_no_symmetry)
    CALL cp_dbcsr_multiply("N","T",1.0_dp,matrix,u2,0.0_dp,tmp,&
               filter_eps=eps_filter)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,u1,tmp,0.0_dp,matrix,&
               filter_eps=eps_filter)
    CALL cp_dbcsr_release(tmp)

    CALL timestop(handle)

  END SUBROUTINE matrix_backward_transform

!! *****************************************************************************
!!> \brief Transforms to a representation in which diagonal blocks
!!>        of qq and pp matrices are diagonal. This can improve convergence
!!>        of PCG
!!> \par History
!!>       2011.07 created [Rustam Z Khaliullin]
!!> \author Rustam Z Khaliullin
!! *****************************************************************************
!  SUBROUTINE transform_matrices_to_blk_diag(matrix_pp,matrix_qq,matrix_qp,&
!    matrix_pq,eps_filter)
!
!    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_pp, matrix_qq,&
!                                                matrix_qp, matrix_pq
!    REAL(KIND=dp), INTENT(IN)                :: eps_filter
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'transform_matrices_to_blk_diag',&
!      routineP = moduleN//':'//routineN
!
!    TYPE(cp_dbcsr_type)                      :: tmp_pp, tmp_qq,&
!                                                tmp_qp, tmp_pq,&
!                                                blk, blk2
!    INTEGER                                  :: handle
!
!    CALL timeset(routineN,handle)
!
!    ! find a better basis by diagonalizing diagonal blocks
!    ! first pp
!    CALL cp_dbcsr_init(blk)
!    CALL cp_dbcsr_create(blk,template=matrix_pp)
!    CALL diagonalize_diagonal_blocks(matrix_pp,blk)
!    
!    ! convert matrices to the new basis
!    CALL cp_dbcsr_init(tmp_pp)
!    CALL cp_dbcsr_create(tmp_pp,template=matrix_pp)
!    CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_pp,blk,0.0_dp,tmp_pp,&
!               filter_eps=eps_filter)
!    CALL cp_dbcsr_multiply("T","N",1.0_dp,blk,tmp_pp,0.0_dp,matrix_pp,&
!               filter_eps=eps_filter)
!    CALL cp_dbcsr_release(tmp_pp)
!
!    ! now qq
!    CALL cp_dbcsr_init(blk2)
!    CALL cp_dbcsr_create(blk2,template=matrix_qq)
!    CALL diagonalize_diagonal_blocks(matrix_qq,blk2)
!
!    CALL cp_dbcsr_init(tmp_qq)
!    CALL cp_dbcsr_create(tmp_qq,template=matrix_qq)
!    CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_qq,blk2,0.0_dp,tmp_qq,&
!               filter_eps=eps_filter)
!    CALL cp_dbcsr_multiply("T","N",1.0_dp,blk2,tmp_qq,0.0_dp,matrix_qq,&
!               filter_eps=eps_filter)
!    CALL cp_dbcsr_release(tmp_qq)
!
!    ! transform pq
!    CALL cp_dbcsr_init(tmp_pq)
!    CALL cp_dbcsr_create(tmp_pq,template=matrix_pq)
!    CALL cp_dbcsr_multiply("T","N",1.0_dp,blk,matrix_pq,0.0_dp,tmp_pq,&
!               filter_eps=eps_filter)
!    CALL cp_dbcsr_multiply("N","N",1.0_dp,tmp_pq,blk2,0.0_dp,matrix_pq,&
!               filter_eps=eps_filter)
!    CALL cp_dbcsr_release(tmp_pq)
!   
!    ! transform qp
!    CALL cp_dbcsr_init(tmp_qp)
!    CALL cp_dbcsr_create(tmp_qp,template=matrix_qp)
!    CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_qp,blk,0.0_dp,tmp_qp,&
!               filter_eps=eps_filter)
!    CALL cp_dbcsr_multiply("T","N",1.0_dp,blk2,tmp_qp,0.0_dp,matrix_qp,&
!               filter_eps=eps_filter)
!    CALL cp_dbcsr_release(tmp_qp)
!
!    CALL cp_dbcsr_release(blk2)
!    CALL cp_dbcsr_release(blk)
!
!    CALL timestop(handle)
!
!  END SUBROUTINE transform_matrices_to_blk_diag

! *****************************************************************************
!> \brief computes oo, ov, vo, and vv blocks of the ks matrix
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
!  SUBROUTINE ct_step_env_execute(env)
!  
!    TYPE(ct_step_env_type)                      :: env
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'ct_step_env_execute', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: handle
!
!    CALL timeset(routineN,handle)
!
!
!    CALL timestop(handle)
!
!  END SUBROUTINE ct_step_env_execute

END MODULE ct_methods


