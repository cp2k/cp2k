!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief an eigen-space solver for the generalised symmetric eigenvalue problem
!>      for sparse matrices, needing only multiplications
!> \author Joost VandeVondele (25.08.2002)
! *****************************************************************************
MODULE qs_ot_eigensolver
  USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                             cp_dbcsr_cholesky_invert
  USE cp_dbcsr_interface,              ONLY: &
       array_i1d_obj, array_release, cp_create_bl_distribution, &
       cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_distribution, &
       cp_dbcsr_distribution_release, cp_dbcsr_init_p, cp_dbcsr_multiply, &
       cp_dbcsr_p_type, cp_dbcsr_release_p, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_trace, cp_dbcsr_type, &
       dbcsr_distribution_mp, dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_row_dist, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_type_no_symmetry, dbcsr_type_real_default
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_copy_columns_hack
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE preconditioner_types,            ONLY: preconditioner_in_use,&
                                             preconditioner_type
  USE qs_mo_methods,                   ONLY: make_basis_sv
  USE qs_ot,                           ONLY: qs_ot_get_orbitals,&
                                             qs_ot_get_p,&
                                             qs_ot_new_preconditioner
  USE qs_ot_minimizer,                 ONLY: ot_mini
  USE qs_ot_types,                     ONLY: qs_ot_allocate,&
                                             qs_ot_destroy,&
                                             qs_ot_init,&
                                             qs_ot_settings_init,&
                                             qs_ot_settings_type,&
                                             qs_ot_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_eigensolver'

! *** Public subroutines ***

  PUBLIC :: ot_eigensolver

CONTAINS

! on input c contains the initial guess (should not be zero !)
! on output c spans the subspace
! *****************************************************************************
!> \brief ...
!> \param matrix_h ...
!> \param matrix_s ...
!> \param matrix_orthogonal_space_fm ...
!> \param matrix_c_fm ...
!> \param preconditioner ...
!> \param eps_gradient ...
!> \param iter_max ...
!> \param size_ortho_space ...
!> \param silent ...
!> \param ot_settings ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE ot_eigensolver(matrix_h,matrix_s,matrix_orthogonal_space_fm, &
                            matrix_c_fm,preconditioner,eps_gradient, &
                            iter_max,size_ortho_space,silent,ot_settings,error)

    TYPE(cp_dbcsr_type), POINTER             :: matrix_h, matrix_s
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: matrix_orthogonal_space_fm
    TYPE(cp_fm_type), POINTER                :: matrix_c_fm
    TYPE(preconditioner_type), OPTIONAL, &
      POINTER                                :: preconditioner
    REAL(KIND=dp)                            :: eps_gradient
    INTEGER, INTENT(IN)                      :: iter_max
    INTEGER, INTENT(IN), OPTIONAL            :: size_ortho_space
    LOGICAL, INTENT(IN), OPTIONAL            :: silent
    TYPE(qs_ot_settings_type), INTENT(IN), &
      OPTIONAL                               :: ot_settings
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_eigensolver', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: max_iter_inner_loop = 40
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, ieigensolver, &
                                                iter_total, k, n, ortho_k, &
                                                ortho_space_k, output_unit
    LOGICAL                                  :: energy_only, my_silent, ortho
    REAL(KIND=dp)                            :: delta, energy
    TYPE(array_i1d_obj) :: col_blk_size, col_blk_size_o, col_dist, &
      col_dist_o, row_blk_size, row_blk_size_o, row_dist, row_dist_o
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_hc
    TYPE(cp_dbcsr_type), POINTER :: matrix_buf1_ortho, matrix_buf2_ortho, &
      matrix_c, matrix_orthogonal_space, matrix_os_ortho, matrix_s_ortho
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_distribution_obj)             :: dist, dist_o_o
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    IF ( PRESENT(silent) ) THEN
      my_silent = silent
    ELSE
      my_silent = .FALSE.
    END IF


    NULLIFY(matrix_c)! fm->dbcsr
    CALL cp_fm_get_info(matrix_c_fm, nrow_global=n, ncol_global=k,error=error)! fm->dbcsr
    CALL cp_create_bl_distribution (col_dist, col_blk_size, k, &! fm->dbcsr
         dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_h))))! fm->dbcsr
    CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp (cp_dbcsr_distribution(matrix_h)),&! fm->dbcsr
         dbcsr_distribution_row_dist(cp_dbcsr_distribution(matrix_h)), col_dist)! fm->dbcsr
    CALL array_release (col_dist)! fm->dbcsr

    CALL cp_dbcsr_init_p(matrix_c, error=error)! fm->dbcsr
    CALL cp_dbcsr_create(matrix_c, "matrix_c", dist, dbcsr_type_no_symmetry,&! fm->dbcsr
         cp_dbcsr_row_block_sizes(matrix_h), col_blk_size, 0, 0, dbcsr_type_real_default, error=error)! fm->dbcsr

    CALL copy_fm_to_dbcsr(matrix_c_fm,matrix_c, error=error)! fm->dbcsr

    CALL cp_dbcsr_distribution_release (dist)! fm->dbcsr! fm->dbcsr
    CALL array_release (col_blk_size)

    iter_total=0

    outer_scf: DO

      NULLIFY(qs_ot_env)

      NULLIFY(matrix_s_ortho)
      NULLIFY(matrix_os_ortho)
      NULLIFY(matrix_buf1_ortho)
      NULLIFY(matrix_buf2_ortho)
      NULLIFY(matrix_orthogonal_space)

      ALLOCATE(qs_ot_env(1))
      ALLOCATE(matrix_hc(1))
      NULLIFY(matrix_hc(1)%matrix)
      CALL cp_dbcsr_init_p(matrix_hc(1)%matrix, error=error)
      CALL cp_dbcsr_copy(matrix_hc(1)%matrix,matrix_c,'matrix_hc',error=error)

      ortho=.FALSE.
      IF (PRESENT(matrix_orthogonal_space_fm)) ortho=.TRUE.

      ! decide settings
      IF ( PRESENT(ot_settings) ) THEN
        qs_ot_env(1)%settings = ot_settings
      ELSE
        CALL qs_ot_settings_init(qs_ot_env(1)%settings)
        ! overwrite defaults
        qs_ot_env(1)%settings%ds_min=0.10_dp
      END IF

      IF (ortho) THEN
         CALL cp_fm_get_info(matrix_orthogonal_space_fm,ncol_global=ortho_space_k,error=error)

         CALL cp_create_bl_distribution (col_dist, col_blk_size, ortho_space_k, &
              dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_h))))
         CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp (cp_dbcsr_distribution(matrix_h)),&
              dbcsr_distribution_row_dist(cp_dbcsr_distribution(matrix_h)), col_dist)
         CALL array_release (col_dist)

         CALL cp_dbcsr_init_p(matrix_orthogonal_space, error=error)
         CALL cp_dbcsr_create(matrix_orthogonal_space, "matrix_orthogonal_space", dist, dbcsr_type_no_symmetry,&
              cp_dbcsr_row_block_sizes(matrix_h), col_blk_size, 0, 0, dbcsr_type_real_default, error=error)

         CALL copy_fm_to_dbcsr(matrix_orthogonal_space_fm,matrix_orthogonal_space, error=error)

         CALL cp_dbcsr_distribution_release (dist)
         CALL array_release (col_blk_size)


         IF (PRESENT(size_ortho_space)) ortho_space_k=size_ortho_space
         ortho_k=ortho_space_k+k
      ELSE
         ortho_k=k
      ENDIF

      ! allocate
      CALL qs_ot_allocate(qs_ot_env(1),matrix_s,matrix_c_fm%matrix_struct,ortho_k=ortho_k,error=error)

      IF (ortho) THEN
          ! construct an initial guess that is orthogonal to matrix_orthogonal_space

         CALL cp_dbcsr_init_p(matrix_s_ortho, error=error)
         CALL cp_dbcsr_copy(matrix_s_ortho,matrix_orthogonal_space,name="matrix_s_ortho",error=error)

         CALL cp_create_bl_distribution (col_dist_o, col_blk_size_o, ortho_space_k, &
              dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_h))))
         CALL cp_create_bl_distribution (row_dist_o, row_blk_size_o, ortho_space_k, &
              dbcsr_mp_nprows(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_h))))
         CALL cp_create_bl_distribution (row_dist, row_blk_size, k, &
              dbcsr_mp_nprows(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_h))))
         CALL cp_create_bl_distribution (col_dist, col_blk_size, k, &
              dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_h))))
         CALL dbcsr_distribution_new (dist_o_o, dbcsr_distribution_mp (cp_dbcsr_distribution(matrix_h)),&
              row_dist_o, col_dist_o)
         CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp (cp_dbcsr_distribution(matrix_h)),&
              row_dist_o, col_dist)
         CALL array_release(row_dist_o); CALL array_release(col_dist_o)
         CALL array_release(row_dist); CALL array_release(col_dist)

         CALL cp_dbcsr_init_p(matrix_os_ortho, error=error)
         CALL cp_dbcsr_create(matrix_os_ortho, "matrix_os_ortho", dist_o_o, dbcsr_type_no_symmetry,&
              row_blk_size_o, col_blk_size_o, 0, 0, dbcsr_type_real_default,&
              error=error)

         CALL cp_dbcsr_init_p(matrix_buf1_ortho, error=error)
         CALL cp_dbcsr_create(matrix_buf1_ortho, "matrix_buf1_ortho", dist, dbcsr_type_no_symmetry,&
              row_blk_size_o, col_blk_size, 0, 0, dbcsr_type_real_default,&
              error=error)

         CALL cp_dbcsr_init_p(matrix_buf2_ortho, error=error)
         CALL cp_dbcsr_create(matrix_buf2_ortho, "matrix_buf2_ortho", dist, dbcsr_type_no_symmetry,&
              row_blk_size_o, col_blk_size, 0, 0, dbcsr_type_real_default,&
              error=error)

         CALL cp_dbcsr_distribution_release (dist_o_o)
         CALL cp_dbcsr_distribution_release (dist)
         CALL array_release(col_blk_size_o); CALL array_release(row_blk_size_o)
         CALL array_release(row_blk_size); CALL array_release(col_blk_size)

          CALL cp_dbcsr_multiply('N','N',1.0_dp,matrix_s,matrix_orthogonal_space, &
                             0.0_dp,matrix_s_ortho,error=error)
          CALL cp_dbcsr_multiply('T','N',rone,matrix_s_ortho, matrix_s_ortho, &
              rzero,matrix_os_ortho,error=error)

          CALL cp_dbcsr_cholesky_decompose(matrix_os_ortho,&
               para_env=qs_ot_env(1)%para_env,blacs_env=qs_ot_env(1)%blacs_env,error=error)
          CALL cp_dbcsr_cholesky_invert(matrix_os_ortho,&
               para_env=qs_ot_env(1)%para_env,blacs_env=qs_ot_env(1)%blacs_env,&
               upper_to_full=.TRUE.,error=error)

          CALL cp_dbcsr_multiply('T','N',rone,matrix_s_ortho,matrix_c, &
               rzero,matrix_buf1_ortho, error=error)
          CALL cp_dbcsr_multiply('N','N',rone,matrix_os_ortho,matrix_buf1_ortho, &
               rzero,matrix_buf2_ortho,error=error)
          CALL cp_dbcsr_multiply('N','N',-rone,matrix_s_ortho, matrix_buf2_ortho, &
               rone,matrix_c, error=error)

          ! make matrix_c0 an orthogonal basis, matrix_c contains sc0
          CALL cp_dbcsr_copy(qs_ot_env(1)%matrix_c0,matrix_c,error=error)
          CALL cp_dbcsr_multiply('N','N',1.0_dp,matrix_s,qs_ot_env(1)%matrix_c0, &
               0.0_dp,matrix_c,error=error)

          CALL make_basis_sv(qs_ot_env(1)%matrix_c0,k, matrix_c, &
               qs_ot_env(1)%para_env, qs_ot_env(1)%blacs_env, &
               error=error)

          ! copy sc0 and matrix_s_ortho in qs_ot_env(1)%matrix_sc0
          !CALL cp_dbcsr_copy_columns(qs_ot_env(1)%matrix_sc0,matrix_s_ortho,ortho_space_k,1,1,error=error)
          CALL cp_dbcsr_copy_columns_hack(qs_ot_env(1)%matrix_sc0,matrix_s_ortho,ortho_space_k,1,1,&
               para_env=qs_ot_env(1)%para_env,blacs_env=qs_ot_env(1)%blacs_env,error=error)
          !CALL cp_dbcsr_copy_columns(qs_ot_env(1)%matrix_sc0,matrix_c,k,1,ortho_space_k+1,error=error)
          CALL cp_dbcsr_copy_columns_hack(qs_ot_env(1)%matrix_sc0,matrix_c,k,1,ortho_space_k+1,&
               para_env=qs_ot_env(1)%para_env,blacs_env=qs_ot_env(1)%blacs_env,error=error)

          CALL cp_dbcsr_release_p(matrix_buf1_ortho, error=error)
          CALL cp_dbcsr_release_p(matrix_buf2_ortho, error=error)
          CALL cp_dbcsr_release_p(matrix_os_ortho, error=error)
          CALL cp_dbcsr_release_p(matrix_s_ortho, error=error)

      ELSE

         ! set c0,sc0
         CALL cp_dbcsr_copy(qs_ot_env(1)%matrix_c0,matrix_c,error=error)
         CALL cp_dbcsr_multiply('N','N',1.0_dp,matrix_s,qs_ot_env(1)%matrix_c0, &
           0.0_dp,qs_ot_env(1)%matrix_sc0,error=error)

         CALL make_basis_sv(qs_ot_env(1)%matrix_c0,k, qs_ot_env(1)%matrix_sc0, &
               qs_ot_env(1)%para_env, qs_ot_env(1)%blacs_env, &
               error=error)
      ENDIF

      ! init
      CALL qs_ot_init(qs_ot_env(1),error=error)
      energy_only=qs_ot_env(1)%energy_only

      ! set x
      CALL cp_dbcsr_set(qs_ot_env(1)%matrix_x,0.0_dp,error=error)
      CALL cp_dbcsr_set(qs_ot_env(1)%matrix_sx,0.0_dp,error=error)

      ! get c
      CALL qs_ot_get_p(qs_ot_env(1)%matrix_x,qs_ot_env(1)%matrix_sx,qs_ot_env(1),error=error)
      CALL qs_ot_get_orbitals(matrix_c,qs_ot_env(1)%matrix_x,qs_ot_env(1),error=error)

      ! if present preconditioner, use it

      IF (PRESENT(preconditioner) ) THEN
         IF (ASSOCIATED(preconditioner)) THEN
           IF (preconditioner_in_use(preconditioner)) THEN
               CALL qs_ot_new_preconditioner(qs_ot_env(1),preconditioner,error=error)
           ELSE
               ! we should presumably make one
           END IF
         ENDIF
      ENDIF

!     *** Eigensolver loop ***
      ieigensolver = 0
      eigensolver_loop: DO

        ieigensolver = ieigensolver + 1
        iter_total=iter_total+1

        ! the energy is cHc, the gradient is 2*H*c
        CALL cp_dbcsr_multiply('N','N',1.0_dp,matrix_h,matrix_c, &
            0.0_dp,matrix_hc(1)%matrix,error=error)
        CALL cp_dbcsr_trace(matrix_c,matrix_hc(1)%matrix,energy,error=error)
        IF (.NOT. energy_only) THEN
           CALL cp_dbcsr_scale(matrix_hc(1)%matrix,2.0_dp,error=error)
        ENDIF

        qs_ot_env(1)%etotal=energy
        CALL ot_mini(qs_ot_env,matrix_hc,output_unit=0,error=error)
        delta           =qs_ot_env(1)%delta
        energy_only     =qs_ot_env(1)%energy_only

        CALL cp_dbcsr_multiply('N','N',1.0_dp,matrix_s,qs_ot_env(1)%matrix_x, &
             0.0_dp, qs_ot_env(1)%matrix_sx, error=error)

        CALL qs_ot_get_p(qs_ot_env(1)%matrix_x,qs_ot_env(1)%matrix_sx,qs_ot_env(1),error=error)
        CALL qs_ot_get_orbitals(matrix_c,qs_ot_env(1)%matrix_x,qs_ot_env(1),error=error)

        ! exit on convergence or if maximum of inner loop  cycles is reached
        IF ( delta < eps_gradient .OR. ieigensolver >= max_iter_inner_loop ) EXIT eigensolver_loop
        ! exit if total number of steps is reached, but not during a line search step
        IF ( iter_total >= iter_max .AND. qs_ot_env(1)%OT_METHOD_FULL /= "OT LS" ) EXIT eigensolver_loop

      END DO eigensolver_loop

      CALL qs_ot_destroy(qs_ot_env(1),error=error)
      DEALLOCATE(qs_ot_env)
      CALL cp_dbcsr_release_p(matrix_hc(1)%matrix, error=error)
      DEALLOCATE(matrix_hc)
      CALL cp_dbcsr_release_p(matrix_orthogonal_space, error=error)

      IF (delta < eps_gradient) THEN
          IF ((output_unit>0) .AND. .NOT.my_silent) THEN
             WRITE(output_unit,*) " Reached convergence in ",iter_total," iterations "
          ENDIF
          EXIT outer_scf
      ENDIF
      IF (iter_total >= iter_max) THEN
         IF (output_unit>0) THEN
           IF ( my_silent ) THEN
             WRITE(output_unit,"(A,T60,E20.10)") "  WARNING OT eigensolver did not converge: current gradient",delta
           ELSE
             WRITE(output_unit,*) "WARNING : did not converge in ot_eigensolver"
             WRITE(output_unit,*) "number of iterations ",iter_total," exceeded maximum"
             WRITE(output_unit,*) "current gradient / target gradient",delta," / ",eps_gradient
           END IF
         ENDIF
         EXIT outer_scf
      ENDIF

    ENDDO outer_scf


    CALL copy_dbcsr_to_fm(matrix_c,matrix_c_fm,error=error)! fm->dbcsr
    CALL cp_dbcsr_release_p(matrix_c, error=error)! fm->dbcsr

    CALL timestop(handle)

  END SUBROUTINE ot_eigensolver

END MODULE qs_ot_eigensolver
