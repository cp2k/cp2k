!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_ot_eigensolver [1.0] *
!!
!!   NAME
!!     qs_ot_eigensolver
!!
!!   FUNCTION
!!     an eigen-space solver for the generalised symmetric eigenvalue problem
!!     for sparse matrices, needing only multiplications
!!
!!   AUTHOR
!!     Joost VandeVondele (25.08.2002)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_ot_eigensolver
! *****************************************************************************
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             wp => dp
  USE preconditioner,                  ONLY: preconditioner_in_use,&
                                             preconditioner_type
  USE qs_blacs,                        ONLY: blacs_column_copy,&
                                             blacs_make_basis,&
                                             copy_blacs_to_blacs_matrix,&
                                             cp_fm_add,&
                                             cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_create2,&
                                             cp_fm_gemm,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_symm,&
                                             cp_fm_trace,&
                                             cp_full_matrix_type,&
                                             cp_sm_fm_multiply,&
                                             write_blacs_matrix
  USE qs_ot,                           ONLY: qs_ot_allocate,&
                                             qs_ot_destroy,&
                                             qs_ot_get_derivative,&
                                             qs_ot_get_orbitals,&
                                             qs_ot_get_p,&
                                             qs_ot_init,&
                                             qs_ot_new_preconditioner,&
                                             qs_ot_type
  USE qs_ot_minimizer,                 ONLY: ot_mini
  USE sparse_matrix_types,             ONLY: real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE


! *** Public subroutines ***

  PUBLIC :: ot_eigensolver

! *****************************************************************************

CONTAINS

! *****************************************************************************
! on input c contains the initial guess (should not be zero !)
! on output c spans the subspace
! *****************************************************************************
  SUBROUTINE ot_eigensolver(matrix_h,matrix_s,matrix_orthogonal_space, &
                                 matrix_c,preconditioner,eps_gradient, &
                                                       iter_max,globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(cp_full_matrix_type), POINTER          :: matrix_c
    TYPE(cp_full_matrix_type), POINTER, OPTIONAL:: matrix_orthogonal_space
    REAL(wp)                                  :: eps_gradient
    TYPE(real_matrix_type), POINTER           :: matrix_h,matrix_s
    INTEGER                                   :: iter_max
    TYPE(preconditioner_type), TARGET , OPTIONAL  :: preconditioner

!   *** Local variables ***

    TYPE(qs_ot_type)                          :: qs_ot_env
    INTEGER                                   :: n,k,nrow_block,ncol_block
    INTEGER                                   :: ieigensolver,handle
    REAL(wp)                                  :: energy,delta
    LOGICAL                                   :: energy_only
    INTEGER                                   :: ortho_k,ortho_space_k
    LOGICAL                                   :: ortho
    REAL(wp), PARAMETER :: rzero = 0.0_wp
    REAL(wp), PARAMETER :: rone  = 1.0_wp
    TYPE(cp_full_matrix_type), POINTER          :: matrix_s_ortho, &
                                                 matrix_os_ortho, &
                                                 matrix_buf1_ortho, &
                                                 matrix_buf2_ortho

    
!   ---------------------------------------------------------------------------

    CALL timeset("ot_eigensolver","I","",handle)

    CALL cp_fm_get_info(matrix_c,nrow_block=nrow_block, &
                                               ncol_block=ncol_block, &
                                               nrow_global=n, &
                                               ncol_global=k)

    ortho=.false.
    IF (present(matrix_orthogonal_space)) ortho=.true.

    ! decide settings
    qs_ot_env%settings%ot_method="CG"
    qs_ot_env%settings%line_search_method="2PNT"
    qs_ot_env%settings%ds_min=0.10_wp
    qs_ot_env%settings%diis_m=7 
    qs_ot_env%settings%n=n 
    qs_ot_env%settings%k=k 
    if (ortho) then
       CALL cp_fm_get_info(matrix_orthogonal_space,&
                                  ncol_global=ortho_space_k)
       qs_ot_env%settings%ortho_k=ortho_space_k+k
    else
       qs_ot_env%settings%ortho_k=k
    endif
    qs_ot_env%settings%nrow_block=nrow_block
    qs_ot_env%settings%ncol_block=ncol_block

    ! allocate
    CALL qs_ot_allocate(qs_ot_env,globenv)

    if (ortho) then
        ! construct an initial guess that is orthogonal to matrix_orthogonal_space
        CALL cp_fm_create2(matrix_s_ortho, &
                                   nrow_block=nrow_block, &
                                   ncol_block=ncol_block, &
                                   nrow_global=n, &
                                   ncol_global=ortho_space_k, &
                                   name="matrix_s_ortho", &
                                   globenv=globenv)

        CALL cp_fm_create2(matrix_os_ortho, &
                                   nrow_block=nrow_block, &
                                   ncol_block=ncol_block, &
                                   nrow_global=ortho_space_k, &
                                   ncol_global=ortho_space_k, &
                                   name="matrix_os_ortho", &
                                   globenv=globenv)

        CALL cp_fm_create2(matrix_buf1_ortho, &
                                   nrow_block=nrow_block, &
                                   ncol_block=ncol_block, &
                                   nrow_global=ortho_space_k, &
                                   ncol_global=k, &
                                   name="matrix_buf1_ortho", &
                                   globenv=globenv)

        CALL cp_fm_create2(matrix_buf2_ortho, &
                                   nrow_block=nrow_block, &
                                   ncol_block=ncol_block, &
                                   nrow_global=ortho_space_k, &
                                   ncol_global=k, &
                                   name="matrix_buf2_ortho", &
                                   globenv=globenv)
        CALL cp_sm_fm_multiply(matrix_s,matrix_orthogonal_space, &
                                matrix_s_ortho,ortho_space_k,para_env=globenv%para_env)

        CALL cp_fm_gemm('T','N',ortho_space_k,ortho_space_k,n,rone,matrix_s_ortho, &
                        matrix_s_ortho, rzero,matrix_os_ortho)
        CALL cp_fm_cholesky_decompose(matrix_os_ortho)
        CALL cp_fm_cholesky_invert(matrix_os_ortho)
         
        CALL cp_fm_gemm('T','N',ortho_space_k,k,n,rone,matrix_s_ortho,matrix_c, &
                        rzero,matrix_buf1_ortho)
        CALL cp_fm_symm('L','U',ortho_space_k,k,rone,matrix_os_ortho, &
                        matrix_buf1_ortho, &
                        rzero,matrix_buf2_ortho)
        CALL cp_fm_gemm('N','N',n,k,ortho_space_k,-rone,matrix_s_ortho, &
                        matrix_buf2_ortho, &
                        rone,matrix_c)

        ! make matrix_c0 an orthogonal basis, matrix_c contains sc0
        CALL copy_blacs_to_blacs_matrix(matrix_c,qs_ot_env%matrix_c0)
        CALL cp_sm_fm_multiply(matrix_s,qs_ot_env%matrix_c0, &
                                matrix_c,k, &
                                para_env=globenv%para_env)
        CALL blacs_make_basis(qs_ot_env%matrix_c0,k, &
                              ortho=matrix_c,otype="SV")

        ! copy sc0 and matrix_s_ortho in qs_ot_env%matrix_sc0
        CALL blacs_column_copy(matrix_s_ortho,qs_ot_env%matrix_sc0,ortho_space_k,1,1)
        CALL blacs_column_copy(matrix_c,qs_ot_env%matrix_sc0,k,1,ortho_space_k+1)

        CALL cp_fm_release(matrix_buf1_ortho)
        CALL cp_fm_release(matrix_buf2_ortho)
        CALL cp_fm_release(matrix_os_ortho)
        CALL cp_fm_release(matrix_s_ortho)

    else

    ! set c0,sc0
       CALL copy_blacs_to_blacs_matrix(matrix_c,qs_ot_env%matrix_c0)
       CALL cp_sm_fm_multiply(matrix_s,qs_ot_env%matrix_c0, &
                               qs_ot_env%matrix_sc0,qs_ot_env%settings%k, &
                               para_env=globenv%para_env)
       CALL blacs_make_basis(qs_ot_env%matrix_c0,k, &
                             ortho=qs_ot_env%matrix_sc0,otype="SV")
    endif
    ! init
    CALL qs_ot_init(qs_ot_env,globenv)
    energy_only=qs_ot_env%energy_only

    ! set x
    CALL cp_fm_set_all(qs_ot_env%matrix_x,0.0_wp)
    CALL cp_fm_set_all(qs_ot_env%matrix_sx,0.0_wp)

    ! get c
    CALL qs_ot_get_p(qs_ot_env%matrix_x,qs_ot_env%matrix_sx,qs_ot_env,globenv)
    CALL qs_ot_get_orbitals(matrix_c,qs_ot_env%matrix_x,qs_ot_env,globenv)

    ! if present preconditioner, use it

    IF (present(preconditioner)) THEN
       IF (preconditioner_in_use(preconditioner)) THEN
          CALL qs_ot_new_preconditioner(qs_ot_env,preconditioner) 
       ELSE
          ! we should presumably make one
       END IF
    ENDIF

!   *** Eigensolver loop ***
    ieigensolver = 0 
    eigensolver_loop: DO

      ieigensolver = ieigensolver + 1

      ! the energy is cHc, the gradient is 2*H*c
      CALL cp_sm_fm_multiply(matrix_h,matrix_c, &
                                     qs_ot_env%matrix_hc,  &
                                     qs_ot_env%settings%k,&
                                     para_env=globenv%para_env)
      CALL cp_fm_trace(matrix_c,qs_ot_env%matrix_hc,energy)
      if (.not. energy_only) then
       CALL cp_fm_add(2.0_wp,qs_ot_env%matrix_hc,0.0_wp,qs_ot_env%matrix_hc)
      endif
 
      qs_ot_env%etotal=energy
      CALL ot_mini(qs_ot_env,globenv)
      delta           =qs_ot_env%delta
      energy_only     =qs_ot_env%energy_only

      CALL cp_sm_fm_multiply(matrix_s,qs_ot_env%matrix_x, &
                                         qs_ot_env%matrix_sx, &
                                         qs_ot_env%settings%k,&
                                     para_env=globenv%para_env)

      CALL qs_ot_get_p(qs_ot_env%matrix_x,qs_ot_env%matrix_sx,qs_ot_env,globenv)
      CALL qs_ot_get_orbitals(matrix_c,qs_ot_env%matrix_x,qs_ot_env,globenv)

      IF (delta < eps_gradient) THEN
          if (globenv%ionode) THEN
             write(globenv%scr,*) " Reached convergence in ",ieigensolver, &
                                  " iterations "
          endif
          EXIT eigensolver_loop
      ENDIF
      IF (ieigensolver .ge. iter_max) THEN
         if (globenv%ionode) then 
           write(globenv%scr,*) "WARNING : did not converge in ot_eigensolver"
           write(globenv%scr,*) "number of iterations ",ieigensolver
           write(globenv%scr,*) "current gradient / target gradient",delta," / ",eps_gradient
         endif
         EXIT eigensolver_loop
      ENDIF

    END DO eigensolver_loop

    call qs_ot_destroy(qs_ot_env,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE ot_eigensolver

END MODULE qs_ot_eigensolver
