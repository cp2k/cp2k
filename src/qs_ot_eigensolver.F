!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
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

  USE kinds, ONLY: int_size,&
                   wp => dp
  USE global_types,             ONLY: global_environment_type

  USE qs_blacs,                 ONLY: allocate_blacs_matrix,&
                                      blacs_add, &
                                      blacs_trace, &
                                      blacs_gemm,&
                                      blacs_matrix_type,&
                                      blacs_symm,&
                                      blacs_make_basis, &
                                      blacs_column_copy, &
                                      get_blacs_matrix_info,&
                                      blacs_cholesky_decompose, &
                                      blacs_cholesky_invert, &
                                      deallocate_blacs_matrix, &
                                      blacs_set_all, &
                                      copy_blacs_to_blacs_matrix, &
                                      sparse_times_blacs, &
                                      write_blacs_matrix
  USE sparse_matrix_types,      ONLY: real_matrix_type
  USE termination,              ONLY: stop_memory,&
                                      stop_program
  USE timings,                  ONLY: timeset,&
                                      timestop
  USE timesl,                   ONLY: cputime
  USE qs_ot,                    ONLY: qs_ot_type,qs_ot_allocate, &
                                      qs_ot_destroy,qs_ot_init,  &
                                      qs_ot_get_orbitals,        & 
                                      qs_ot_get_derivative,      &
                                      qs_ot_get_p,               &
                                      qs_ot_mini
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
                                 matrix_c,eps_gradient,iter_max,globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(blacs_matrix_type), POINTER          :: matrix_c
    TYPE(blacs_matrix_type), POINTER, OPTIONAL:: matrix_orthogonal_space
    REAL(wp)                                  :: eps_gradient
    TYPE(real_matrix_type), POINTER           :: matrix_h,matrix_s
    INTEGER                                   :: iter_max

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
    TYPE(blacs_matrix_type), POINTER          :: matrix_s_ortho, &
                                                 matrix_os_ortho, &
                                                 matrix_buf1_ortho, &
                                                 matrix_buf2_ortho

    
!   ---------------------------------------------------------------------------

    CALL timeset("ot_eigensolver","I","",handle)

    CALL get_blacs_matrix_info(matrix_c,nrow_block=nrow_block, &
                                               ncol_block=ncol_block, &
                                               nrow_global=n, &
                                               ncol_global=k)

    ortho=.false.
    IF (present(matrix_orthogonal_space)) ortho=.true.

    ! decide settings
    qs_ot_env%settings%ot_method="CG"
    qs_ot_env%settings%ot_linesearch=.true.
    qs_ot_env%settings%ot_diis=.false.
    qs_ot_env%settings%ds_min=0.10_wp
    if (qs_ot_env%settings%ot_method.eq."QN") then
       qs_ot_env%settings%lbfgs_m=10
    else
       qs_ot_env%settings%lbfgs_m=0
    endif
    if (qs_ot_env%settings%ot_diis) then
       qs_ot_env%settings%diis_m=10 
    else
       qs_ot_env%settings%diis_m=0
    endif
    qs_ot_env%settings%n=n 
    qs_ot_env%settings%k=k 
    if (ortho) then
       CALL get_blacs_matrix_info(matrix_orthogonal_space,&
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
        CALL allocate_blacs_matrix(matrix_s_ortho, &
                                   nrow_block=nrow_block, &
                                   ncol_block=ncol_block, &
                                   nrow_global=n, &
                                   ncol_global=ortho_space_k, &
                                   name="matrix_s_ortho", &
                                   globenv=globenv)

        CALL allocate_blacs_matrix(matrix_os_ortho, &
                                   nrow_block=nrow_block, &
                                   ncol_block=ncol_block, &
                                   nrow_global=ortho_space_k, &
                                   ncol_global=ortho_space_k, &
                                   name="matrix_os_ortho", &
                                   globenv=globenv)

        CALL allocate_blacs_matrix(matrix_buf1_ortho, &
                                   nrow_block=nrow_block, &
                                   ncol_block=ncol_block, &
                                   nrow_global=ortho_space_k, &
                                   ncol_global=k, &
                                   name="matrix_buf1_ortho", &
                                   globenv=globenv)

        CALL allocate_blacs_matrix(matrix_buf2_ortho, &
                                   nrow_block=nrow_block, &
                                   ncol_block=ncol_block, &
                                   nrow_global=ortho_space_k, &
                                   ncol_global=k, &
                                   name="matrix_buf2_ortho", &
                                   globenv=globenv)
        CALL sparse_times_blacs(matrix_s,matrix_orthogonal_space, &
                                matrix_s_ortho,ortho_space_k,para_env=globenv%para_env)

        CALL blacs_gemm('T','N',ortho_space_k,ortho_space_k,n,rone,matrix_s_ortho, &
                        matrix_s_ortho, rzero,matrix_os_ortho)
        CALL blacs_cholesky_decompose(matrix_os_ortho)
        CALL blacs_cholesky_invert(matrix_os_ortho)
         
        CALL blacs_gemm('T','N',ortho_space_k,k,n,rone,matrix_s_ortho,matrix_c, &
                        rzero,matrix_buf1_ortho)
        CALL blacs_symm('L','U',ortho_space_k,k,rone,matrix_os_ortho, &
                        matrix_buf1_ortho, &
                        rzero,matrix_buf2_ortho)
        CALL blacs_gemm('N','N',n,k,ortho_space_k,-rone,matrix_s_ortho, &
                        matrix_buf2_ortho, &
                        rone,matrix_c)

        ! make matrix_c0 an orthogonal basis, matrix_c contains sc0
        CALL copy_blacs_to_blacs_matrix(matrix_c,qs_ot_env%matrix_c0)
        CALL sparse_times_blacs(matrix_s,qs_ot_env%matrix_c0, &
                                matrix_c,k, &
                                para_env=globenv%para_env)
        CALL blacs_make_basis(qs_ot_env%matrix_c0,k, &
                              ortho=matrix_c,otype="SV")

        ! copy sc0 and matrix_s_ortho in qs_ot_env%matrix_sc0
        CALL blacs_column_copy(matrix_s_ortho,qs_ot_env%matrix_sc0,ortho_space_k,1,1)
        CALL blacs_column_copy(matrix_c,qs_ot_env%matrix_sc0,k,1,ortho_space_k+1)

        CALL deallocate_blacs_matrix(matrix_buf1_ortho)
        CALL deallocate_blacs_matrix(matrix_buf2_ortho)
        CALL deallocate_blacs_matrix(matrix_os_ortho)
        CALL deallocate_blacs_matrix(matrix_s_ortho)

    else

    ! set c0,sc0
       CALL copy_blacs_to_blacs_matrix(matrix_c,qs_ot_env%matrix_c0)
       CALL sparse_times_blacs(matrix_s,qs_ot_env%matrix_c0, &
                               qs_ot_env%matrix_sc0,qs_ot_env%settings%k, &
                               para_env=globenv%para_env)
       CALL blacs_make_basis(qs_ot_env%matrix_c0,k, &
                             ortho=qs_ot_env%matrix_sc0,otype="SV")
    endif
    ! init
    CALL qs_ot_init(qs_ot_env,globenv)
    energy_only=qs_ot_env%energy_only

    ! set x
    CALL blacs_set_all(qs_ot_env%matrix_x,0.0_wp)
    CALL blacs_set_all(qs_ot_env%matrix_sx,0.0_wp)

    ! get c
    CALL qs_ot_get_p(qs_ot_env%matrix_x,qs_ot_env%matrix_sx,qs_ot_env,globenv)
    CALL qs_ot_get_orbitals(matrix_c,qs_ot_env%matrix_x,qs_ot_env,globenv)

!   *** Eigensolver loop ***
    ieigensolver = 0 
    eigensolver_loop: DO

      ieigensolver = ieigensolver + 1

      ! the energy is cHc, the gradient is 2*H*c
      CALL sparse_times_blacs(matrix_h,matrix_c, &
                                     qs_ot_env%matrix_hc,  &
                                     qs_ot_env%settings%k,&
                                     para_env=globenv%para_env)
      CALL blacs_trace(matrix_c,qs_ot_env%matrix_hc,energy)
      if (.not. energy_only) then
       CALL blacs_add(2.0_wp,qs_ot_env%matrix_hc,0.0_wp,qs_ot_env%matrix_hc)
      endif
 
      qs_ot_env%etotal=energy
      CALL qs_ot_mini(qs_ot_env,globenv)
      delta           =qs_ot_env%delta
      energy_only     =qs_ot_env%energy_only

      CALL sparse_times_blacs(matrix_s,qs_ot_env%matrix_x, &
                                         qs_ot_env%matrix_sx, &
                                         qs_ot_env%settings%k,&
                                     para_env=globenv%para_env)

      CALL qs_ot_get_p(qs_ot_env%matrix_x,qs_ot_env%matrix_sx,qs_ot_env,globenv)
      CALL qs_ot_get_orbitals(matrix_c,qs_ot_env%matrix_x,qs_ot_env,globenv)

      IF (delta < eps_gradient) EXIT eigensolver_loop
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
