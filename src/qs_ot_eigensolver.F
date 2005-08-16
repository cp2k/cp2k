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
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_symm,&
                                             cp_fm_trace
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE preconditioner,                  ONLY: preconditioner_in_use,&
                                             preconditioner_type
  USE qs_mo_methods,                   ONLY: make_basis
  USE qs_ot,                           ONLY: qs_ot_allocate,&
                                             qs_ot_destroy,&
                                             qs_ot_get_orbitals,&
                                             qs_ot_get_p,&
                                             qs_ot_init,&
                                             qs_ot_new_preconditioner,&
                                             qs_ot_settings_init,&
                                             qs_ot_type
  USE qs_ot_minimizer,                 ONLY: ot_mini
  USE sparse_matrix_types,             ONLY: real_matrix_type
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
                            iter_max,globenv,size_ortho_space)

    TYPE(real_matrix_type), POINTER          :: matrix_h, matrix_s
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: matrix_orthogonal_space
    TYPE(cp_fm_type), POINTER                :: matrix_c
    TYPE(preconditioner_type), OPTIONAL, &
      POINTER                                :: preconditioner
    REAL(KIND=dp)                            :: eps_gradient
    INTEGER                                  :: iter_max
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    INTEGER, OPTIONAL                        :: size_ortho_space

    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, ieigensolver, k, n, &
                                                ortho_k, ortho_space_k
    LOGICAL                                  :: energy_only, ortho
    REAL(KIND=dp)                            :: delta, energy
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrix_hc
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_buf1_ortho, &
                                                matrix_buf2_ortho, &
                                                matrix_os_ortho, &
                                                matrix_s_ortho
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env

!   ---------------------------------------------------------------------------

    CALL timeset("ot_eigensolver","I","",handle)
    NULLIFY(qs_ot_env,fm_struct_tmp,matrix_s_ortho,matrix_os_ortho, &
                      matrix_buf1_ortho,matrix_buf2_ortho)
    ALLOCATE(qs_ot_env(1))

    ALLOCATE(matrix_hc(1))
    CALL cp_fm_create(matrix_hc(1)%matrix,matrix_c%matrix_struct)

    CALL cp_fm_get_info(matrix_c, nrow_global=n, ncol_global=k)

    ortho=.FALSE.
    IF (PRESENT(matrix_orthogonal_space)) ortho=.TRUE.

    ! decide settings
    CALL qs_ot_settings_init(qs_ot_env(1)%settings)
    ! overwrite defaults 
    qs_ot_env(1)%settings%ds_min=0.10_dp

    IF (ortho) THEN
       CALL cp_fm_get_info(matrix_orthogonal_space,&
                                  ncol_global=ortho_space_k)
       IF (PRESENT(size_ortho_space)) ortho_space_k=size_ortho_space
       ortho_k=ortho_space_k+k
    ELSE
       ortho_k=k
    ENDIF

    ! allocate
    CALL qs_ot_allocate(qs_ot_env(1),matrix_c%matrix_struct,qs_ot_env(1)%settings%nprec,ortho_k)

    IF (ortho) THEN
        ! construct an initial guess that is orthogonal to matrix_orthogonal_space
        CALL cp_fm_create(matrix_s_ortho,matrix_orthogonal_space%matrix_struct,"matrix_s_ortho")

        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ortho_space_k, ncol_global=ortho_space_k, &
                                 para_env=matrix_orthogonal_space%matrix_struct%para_env, &
                                 context=matrix_orthogonal_space%matrix_struct%context)
        CALL cp_fm_create(matrix_os_ortho, fm_struct_tmp,"matrix_os_ortho") 
        CALL cp_fm_struct_release(fm_struct_tmp)

        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ortho_space_k, ncol_global=k, &
                                 para_env=matrix_orthogonal_space%matrix_struct%para_env, &
                                 context=matrix_orthogonal_space%matrix_struct%context)
        CALL cp_fm_create(matrix_buf1_ortho, fm_struct_tmp,"matrix_buf1_ortho")
        CALL cp_fm_create(matrix_buf2_ortho, fm_struct_tmp,"matrix_buf2_ortho")
        CALL cp_fm_struct_release(fm_struct_tmp)

        CALL cp_sm_fm_multiply(matrix_s,matrix_orthogonal_space, &
                                matrix_s_ortho,ortho_space_k)

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
        CALL cp_fm_to_fm(matrix_c,qs_ot_env(1)%matrix_c0)
        CALL cp_sm_fm_multiply(matrix_s,qs_ot_env(1)%matrix_c0, &
                                matrix_c,k)
        CALL make_basis(qs_ot_env(1)%matrix_c0,k, ortho=matrix_c,otype="SV")

        ! copy sc0 and matrix_s_ortho in qs_ot_env(1)%matrix_sc0
        CALL cp_fm_to_fm(matrix_s_ortho,qs_ot_env(1)%matrix_sc0,ortho_space_k,1,1)
        CALL cp_fm_to_fm(matrix_c,qs_ot_env(1)%matrix_sc0,k,1,ortho_space_k+1)

        CALL cp_fm_release(matrix_buf1_ortho)
        CALL cp_fm_release(matrix_buf2_ortho)
        CALL cp_fm_release(matrix_os_ortho)
        CALL cp_fm_release(matrix_s_ortho)

    ELSE

    ! set c0,sc0
       CALL cp_fm_get_info(matrix_c,ncol_global=k)
       CALL cp_fm_to_fm(matrix_c,qs_ot_env(1)%matrix_c0)
       CALL cp_sm_fm_multiply(matrix_s,qs_ot_env(1)%matrix_c0, &
                               qs_ot_env(1)%matrix_sc0,k)
       CALL make_basis(qs_ot_env(1)%matrix_c0,k, ortho=qs_ot_env(1)%matrix_sc0,otype="SV")
    ENDIF
    ! init
    CALL qs_ot_init(qs_ot_env(1))
    energy_only=qs_ot_env(1)%energy_only

    ! set x
    CALL cp_fm_set_all(qs_ot_env(1)%matrix_x,0.0_dp)
    CALL cp_fm_set_all(qs_ot_env(1)%matrix_sx,0.0_dp)

    ! get c
    CALL qs_ot_get_p(qs_ot_env(1)%matrix_x,qs_ot_env(1)%matrix_sx,qs_ot_env(1))
    CALL qs_ot_get_orbitals(matrix_c,qs_ot_env(1)%matrix_x,qs_ot_env(1))

    ! if present preconditioner, use it

    IF (PRESENT(preconditioner) ) THEN
       IF (ASSOCIATED(preconditioner)) THEN
         IF (preconditioner_in_use(preconditioner)) THEN
             CALL qs_ot_new_preconditioner(qs_ot_env(1),preconditioner,qs_ot_env(1)%settings%nprec) 
         ELSE
             ! we should presumably make one
         END IF
       ENDIF
    ENDIF

!   *** Eigensolver loop ***
    CALL cp_fm_get_info(matrix_c,ncol_global=k)
    ieigensolver = 0 
    eigensolver_loop: DO

      ieigensolver = ieigensolver + 1

      ! the energy is cHc, the gradient is 2*H*c
      CALL cp_sm_fm_multiply(matrix_h,matrix_c, &
                                     matrix_hc(1)%matrix,k)
      CALL cp_fm_trace(matrix_c,matrix_hc(1)%matrix,energy)
      IF (.NOT. energy_only) THEN
       CALL cp_fm_scale_and_add(2.0_dp,matrix_hc(1)%matrix)
      ENDIF
 
      qs_ot_env(1)%etotal=energy
      CALL ot_mini(qs_ot_env,matrix_hc)
      delta           =qs_ot_env(1)%delta
      energy_only     =qs_ot_env(1)%energy_only

      CALL cp_sm_fm_multiply(matrix_s,qs_ot_env(1)%matrix_x, &
                                         qs_ot_env(1)%matrix_sx, k)

      CALL qs_ot_get_p(qs_ot_env(1)%matrix_x,qs_ot_env(1)%matrix_sx,qs_ot_env(1))
      CALL qs_ot_get_orbitals(matrix_c,qs_ot_env(1)%matrix_x,qs_ot_env(1))

      IF (delta < eps_gradient) THEN
          IF (globenv%ionode) THEN
             WRITE(globenv%scr,*) " Reached convergence in ",ieigensolver, &
                                  " iterations "
          ENDIF
          EXIT eigensolver_loop
      ENDIF
      IF (ieigensolver .GE. iter_max) THEN
         IF (globenv%ionode) THEN 
           WRITE(globenv%scr,*) "WARNING : did not converge in ot_eigensolver"
           WRITE(globenv%scr,*) "number of iterations ",ieigensolver
           WRITE(globenv%scr,*) "current gradient / target gradient",delta," / ",eps_gradient
         ENDIF
         EXIT eigensolver_loop
      ENDIF

    END DO eigensolver_loop

    CALL qs_ot_destroy(qs_ot_env(1))
    DEALLOCATE(qs_ot_env)
    CALL cp_fm_release(matrix_hc(1)%matrix)
    DEALLOCATE(matrix_hc)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE ot_eigensolver

END MODULE qs_ot_eigensolver
!!*****
