!******************************************************************************
!!
!!   NAME
!!     qs_ot
!!
!!   FUNCTION
!!     orbital transformations
!!
!!   AUTHOR
!!     Joost VandeVondele (06.2002)
!!
!!   MODIFICATION HISTORY
!!     None
!!
!!   SOURCE
!******************************************************************************

MODULE qs_ot

! *****************************************************************************
  USE kinds, ONLY: wp => dp, wp_size => dp_size
  USE termination,         ONLY: stop_memory, stop_program
  USE timings,             ONLY: timeset, timestop
  USE timesl,              ONLY: cputime
  USE qs_blacs,            ONLY: blacs_matrix_type,get_blacs_matrix_info, &
                                 blacs_set_element, &
                                 blacs_set_all, blacs_syrk, &
                                 allocate_blacs_matrix, &
                                 deallocate_blacs_matrix, &
                                 blacs_symm,blacs_gemm, blacs_set_element, &
                                 blacs_make_basis,copy_blacs_to_blacs_matrix, &
                                 blacs_init_random,write_blacs_matrix, &
                                 blacs_trace, blacs_transpose, blacs_add, &
                                 blacs_schur_product, &
                                 copy_blacs_to_blacs_matrix, &
                                 blacs_cholesky_decompose, &
                                 blacs_cholesky_invert, &
                                 blacs_cholesky_reduce, &
                                 blacs_syevx, blacs_syevd, &
                                 write_blacs_matrix, &
                                 blacs_scale_matrix, &
                                 blacs_set_local_block

  USE global_types,        ONLY: global_environment_type

  IMPLICIT NONE

  PRIVATE :: qs_ot_get_orbitals
  PUBLIC  :: qs_ot_test_it
  PRIVATE :: qs_ot_p2m
  PRIVATE :: qs_ot_sinc
  PRIVATE :: qs_ot_init
  PRIVATE :: qs_ot_destroy
  PRIVATE :: qs_ot_type
  PRIVATE :: qs_ot_get_p
  PRIVATE :: qs_ot_get_derivative
  PRIVATE :: qs_ot_ortho_to_sc0

  TYPE qs_ot_type
     ! these will change during iterations
     TYPE(blacs_matrix_type), POINTER :: matrix_p
     TYPE(blacs_matrix_type), POINTER :: matrix_r
     TYPE(blacs_matrix_type), POINTER :: matrix_sinp
     TYPE(blacs_matrix_type), POINTER :: matrix_cosp
     TYPE(blacs_matrix_type), POINTER :: matrix_sinp_b
     TYPE(blacs_matrix_type), POINTER :: matrix_cosp_b
     TYPE(blacs_matrix_type), POINTER :: matrix_buf1
     TYPE(blacs_matrix_type), POINTER :: matrix_buf2
     TYPE(blacs_matrix_type), POINTER :: matrix_buf3

     REAL(wp),  DIMENSION(:), POINTER :: evals
     REAL(wp),  DIMENSION(:), POINTER :: dum

     ! for efficient/parallel writing to the blacs_matrix
     INTEGER :: nrow_local
     INTEGER :: ncol_local
     REAL(wp),  DIMENSION(:,:), POINTER :: local_block
     INTEGER, DIMENSION(:), POINTER :: row_indices
     INTEGER, DIMENSION(:), POINTER :: col_indices

     ! this is the constraining set
     TYPE(blacs_matrix_type), POINTER :: matrix_c0
     TYPE(blacs_matrix_type), POINTER :: matrix_sc0
     TYPE(blacs_matrix_type), POINTER :: matrix_os

     ! some useful info about the size
     INTEGER                          :: n,k
  END TYPE qs_ot_type

CONTAINS

!*******************************************************************************
! init qs_ot, assumes c0*S*c0=1
! allocates data, keeps a pointer to c0 and sc0. 
! Notice that these are not copied
!*******************************************************************************
SUBROUTINE qs_ot_init(qs_ot_env,matrix_c0,matrix_sc0,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env 
   TYPE(blacs_matrix_type), POINTER          :: matrix_c0,matrix_sc0
   TYPE(global_environment_type), INTENT(IN) :: globenv

   INTEGER :: nrow_block,ncol_block,k,n
   REAL(wp), PARAMETER :: rone=1.0_wp, rzero=0.0_wp

   CALL get_blacs_matrix_info(matrix_c0,nrow_global=n,ncol_global=k,&
                                nrow_block=nrow_block,ncol_block=ncol_block)
   qs_ot_env%n=n
   qs_ot_env%k=k
   qs_ot_env%matrix_c0 =>matrix_c0
   qs_ot_env%matrix_sc0=>matrix_sc0

   allocate(qs_ot_env%evals(k))
   allocate(qs_ot_env%dum(k))

   CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_os,&
                              nrow_global=k,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_os",&
                              globenv=globenv)

   CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_p,&
                              nrow_global=k,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_p",&
                              globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_r,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_r",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_sinp,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_sinp",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_cosp,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_cosp",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_sinp_b,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_sinp_b",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_cosp_b,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_cosp_b",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_buf1,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_buf1",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_buf2,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_buf2",&
                             globenv=globenv)

  CALL allocate_blacs_matrix(new_matrix=qs_ot_env%matrix_buf3,&
                             nrow_global=k,&
                             ncol_global=k,&
                             nrow_block=nrow_block,&
                             ncol_block=ncol_block,&
                             name="matrix_buf3",&
                             globenv=globenv)

   ! get the matrix_os needed for ortho against sc0
   CALL blacs_gemm('T','N',k,k,n,rone,matrix_sc0,matrix_sc0,rzero,&
                             qs_ot_env%matrix_os,globenv)
   CALL blacs_cholesky_decompose(qs_ot_env%matrix_os,globenv)
   CALL blacs_cholesky_invert(qs_ot_env%matrix_os,globenv)

   !CALL blacs_gemm('T','N',k,k,n,rone,matrix_c0,matrix_sc0,rzero,&
   !                          qs_ot_env%matrix_buf3,globenv)
   !CALL write_blacs_matrix(qs_ot_env%matrix_buf3,6,globenv,.true.)
 
   CALL get_blacs_matrix_info(qs_ot_env%matrix_cosp_b, &
                              nrow_local=qs_ot_env%nrow_local, &
                              ncol_local=qs_ot_env%ncol_local, &
                              row_indices=qs_ot_env%row_indices, &
                              col_indices=qs_ot_env%col_indices, &
                              globenv=globenv)
   allocate(qs_ot_env%local_block(qs_ot_env%nrow_local,qs_ot_env%ncol_local))


END SUBROUTINE qs_ot_init

!******************************************************************************
! deallocates data
!******************************************************************************
SUBROUTINE qs_ot_destroy(qs_ot_env,globenv)
   TYPE(qs_ot_type)                          :: qs_ot_env 
   TYPE(global_environment_type), INTENT(IN) :: globenv

   qs_ot_env%n=0
   qs_ot_env%k=0

   deallocate(qs_ot_env%evals)
   deallocate(qs_ot_env%dum)
   deallocate(qs_ot_env%local_block)
   deallocate(qs_ot_env%row_indices)
   deallocate(qs_ot_env%col_indices)
   call deallocate_blacs_matrix(qs_ot_env%matrix_os)
   call deallocate_blacs_matrix(qs_ot_env%matrix_p)
   call deallocate_blacs_matrix(qs_ot_env%matrix_cosp)
   call deallocate_blacs_matrix(qs_ot_env%matrix_sinp)
   call deallocate_blacs_matrix(qs_ot_env%matrix_r)
   call deallocate_blacs_matrix(qs_ot_env%matrix_cosp_b)
   call deallocate_blacs_matrix(qs_ot_env%matrix_sinp_b)
   call deallocate_blacs_matrix(qs_ot_env%matrix_buf1)
   call deallocate_blacs_matrix(qs_ot_env%matrix_buf2)
   call deallocate_blacs_matrix(qs_ot_env%matrix_buf3)
   NULLIFY(qs_ot_env%matrix_c0)
   NULLIFY(qs_ot_env%matrix_sc0)

END SUBROUTINE qs_ot_destroy

SUBROUTINE qs_ot_test_it(globenv)

   TYPE(global_environment_type), INTENT(IN) :: globenv
   !INTEGER, PARAMETER :: n=6,k=3,nrow_block=1,ncol_block=1
   INTEGER, PARAMETER :: n=1472,k=256,nrow_block=32,ncol_block=32
   !INTEGER, PARAMETER :: n=2944,k=512,nrow_block=32,ncol_block=32
   !INTEGER, PARAMETER :: n=5888,k=1024,nrow_block=32,ncol_block=32

   TYPE(blacs_matrix_type), POINTER :: matrix_h,matrix_s
   TYPE(blacs_matrix_type), POINTER :: matrix_c0,matrix_x,matrix_xp,matrix_xm,matrix_dx 
   TYPE(blacs_matrix_type), POINTER :: matrix_sc0,matrix_sx,matrix_c,matrix_hc
   TYPE(qs_ot_type)                 :: qs_ot_env

   REAL(wp) :: alpha,energy,de,ds,energy_old,timestart,timedum,dt,kin_energy
   REAL(wp), PARAMETER :: rone=1.0_wp, rzero=0.0_wp
   REAL(wp), DIMENSION(:), POINTER :: evals
   INTEGER  :: i

   allocate(evals(n))

   CALL allocate_blacs_matrix(new_matrix=matrix_h,&
                              nrow_global=n,&
                              ncol_global=n,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_h",&
                              globenv=globenv)
   CALL allocate_blacs_matrix(new_matrix=matrix_s,&
                              nrow_global=n,&
                              ncol_global=n,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_s",&
                              globenv=globenv)
   CALL allocate_blacs_matrix(new_matrix=matrix_c,&
                              nrow_global=n,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_c",&
                              globenv=globenv)
   CALL allocate_blacs_matrix(new_matrix=matrix_hc,&
                              nrow_global=n,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_hc",&
                              globenv=globenv)
   CALL allocate_blacs_matrix(new_matrix=matrix_c0,&
                              nrow_global=n,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_c0",&
                              globenv=globenv)
   CALL allocate_blacs_matrix(new_matrix=matrix_x,&
                              nrow_global=n,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_x",&
                              globenv=globenv)
   CALL allocate_blacs_matrix(new_matrix=matrix_xp,&
                              nrow_global=n,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_x",&
                              globenv=globenv)
   CALL allocate_blacs_matrix(new_matrix=matrix_xm,&
                              nrow_global=n,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_x",&
                              globenv=globenv)
   CALL allocate_blacs_matrix(new_matrix=matrix_sc0,&
                              nrow_global=n,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_sc0",&
                              globenv=globenv)
   CALL allocate_blacs_matrix(new_matrix=matrix_sx,&
                              nrow_global=n,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_sx",&
                              globenv=globenv)
   CALL allocate_blacs_matrix(new_matrix=matrix_dx,&
                              nrow_global=n,&
                              ncol_global=k,&
                              nrow_block=nrow_block,&
                              ncol_block=ncol_block,&
                              name="matrix_dx",&
                              globenv=globenv)

   CALL blacs_set_all(matrix_h,0.0_wp,globenv)
   CALL blacs_set_all(matrix_s,0.0_wp,globenv)
   DO i=1,n
      alpha=1+real(i,wp)
      CALL blacs_set_element(matrix_h,i,i,alpha,globenv)
      alpha=1.0_wp
      alpha=1+sqrt(real(i,wp))
      CALL blacs_set_element(matrix_s,i,i,alpha,globenv)
   ENDDO
   DO i=1,n-1
      alpha=sqrt(real(n-i,wp))/n
      CALL blacs_set_element(matrix_h,i,i+1,alpha,globenv)
      CALL blacs_set_element(matrix_h,i+1,i,alpha,globenv)
      alpha=(1.0_wp)/sqrt(real(n-i,wp))
      CALL blacs_set_element(matrix_s,i,i+1,alpha,globenv)
      CALL blacs_set_element(matrix_s,i+1,i,alpha,globenv)
   ENDDO

   ! init c0
   CALL blacs_init_random(matrix_c0,k,globenv)
   CALL blacs_gemm('N','N',n,k,n,rone,matrix_s,matrix_c0,rzero,matrix_sc0,globenv)
   CALL blacs_make_basis(matrix_c0,k,globenv,ortho=matrix_sc0,otype="SV")

   CALL qs_ot_init(qs_ot_env,matrix_c0,matrix_sc0,globenv)

   ! init x
   CALL blacs_init_random(matrix_x,k,globenv)
   CALL qs_ot_ortho_to_sc0(matrix_x,qs_ot_env,globenv)
   CALL copy_blacs_to_blacs_matrix(matrix_x,matrix_xm)
  

   ! start loop
   energy_old=1E+10_wp
   ds=-0.000001_wp
   ds=-0.01_wp
   dt= 0.1_wp
   alpha = 0.97_wp
   DO i=1,2000 
      CALL blacs_gemm('N','N',n,k,n,rone,matrix_s,matrix_x,rzero,matrix_sx,globenv)
      CALL qs_ot_get_p(matrix_x,matrix_sx,qs_ot_env,globenv)
      CALL qs_ot_get_orbitals(matrix_c,matrix_x,qs_ot_env,globenv)
      CALL blacs_gemm('N','N',n,k,n,rone,matrix_h,matrix_c,rzero,matrix_hc,globenv)

      CALL blacs_trace(matrix_c,matrix_hc,energy,globenv)
      energy=energy*0.5_wp
      CALL qs_ot_get_derivative(matrix_hc,matrix_x,matrix_sx, &
                                        matrix_dx,qs_ot_env,globenv)
      CALL copy_blacs_to_blacs_matrix(matrix_dx,matrix_xp)
      CALL blacs_add(-dt*dt,matrix_xp,(1.0_wp+alpha),matrix_x,globenv) 
      CALL blacs_add(1.0_wp,matrix_xp,-alpha,matrix_xm,globenv)
      CALL blacs_add(-1.0_wp,matrix_xm,1.0_wp,matrix_xp,globenv)
      CALL blacs_trace(matrix_xm,matrix_xm,kin_energy,globenv)
      kin_energy=0.5_wp*kin_energy/((2.0_wp*dt)**2)
      if (globenv%ionode) write(6,*) i,energy,kin_energy,energy+kin_energy
 
      
      CALL copy_blacs_to_blacs_matrix(matrix_x,matrix_xm)
      CALL copy_blacs_to_blacs_matrix(matrix_xp,matrix_x)
      !CALL blacs_add(1.0_wp,matrix_x,ds,matrix_dx,globenv)
      !CALL blacs_trace(matrix_dx,matrix_dx,de,globenv)
      !write(6,*) "ds,dg*dg,energy,de" &
      !       ,ds,de,energy,-(energy_old-energy)/ds,i
      !energy_old=energy
      !IF (de.lt.1E-18_wp) EXIT
   ENDDO

   CALL blacs_cholesky_decompose(matrix_s,globenv)
   CALL blacs_cholesky_reduce(matrix_h,matrix_s,globenv)
   CALL blacs_syevx(matrix_h,matrix_s,evals,n,1.0_wp,globenv)
   if (globenv%ionode) write(6,*) "evals",evals
   do i=2,k
      evals(1)=evals(1)+evals(i)
   enddo
   if (globenv%ionode) write(6,*) "diag  ",evals(1)/2.0_wp

   deallocate(evals)
   CALL deallocate_blacs_matrix(matrix_h)
   CALL deallocate_blacs_matrix(matrix_s)
   CALL deallocate_blacs_matrix(matrix_c)
   CALL deallocate_blacs_matrix(matrix_hc)
   CALL deallocate_blacs_matrix(matrix_c0)
   CALL deallocate_blacs_matrix(matrix_x)
   CALL deallocate_blacs_matrix(matrix_sc0)
   CALL deallocate_blacs_matrix(matrix_sx)
   CALL deallocate_blacs_matrix(matrix_dx)

   CALL qs_ot_destroy(qs_ot_env,globenv)

   STOP
END SUBROUTINE qs_ot_test_it

!*******************************************************************************
! computes p=x*S*x
!*******************************************************************************
SUBROUTINE qs_ot_get_p(matrix_x,matrix_sx,qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN) :: globenv
  TYPE(blacs_matrix_type), POINTER          :: matrix_sx,matrix_x
  TYPE(qs_ot_type)                          :: qs_ot_env



  REAL(wp) , PARAMETER            :: rone=1.0_wp , rzero=0.0_wp
  INTEGER                         :: n,k,handle
  REAL(wp)                        :: timedum

  CALL timeset("qs_ot_get_p","I","",handle)

  n=qs_ot_env%n
  k=qs_ot_env%k

! get the overlap
  CALL blacs_gemm('T','N',k,k,n,rone,matrix_x,matrix_sx,rzero,& 
                           qs_ot_env%matrix_p, globenv)

  CALL qs_ot_p2m(qs_ot_env,globenv)

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_ot_get_p

!*******************************************************************************
! makes a set of vectors orthogonal to sc0 
!*******************************************************************************
SUBROUTINE qs_ot_ortho_to_sc0(matrix_x,qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN) :: globenv
  TYPE(blacs_matrix_type), POINTER          :: matrix_x
  TYPE(qs_ot_type)                          :: qs_ot_env

  REAL(wp) , PARAMETER            :: rone=1.0_wp , rzero=0.0_wp
  INTEGER                         :: n,k,handle

  CALL timeset("qs_ot_ortho","I","",handle)

  n=qs_ot_env%n
  k=qs_ot_env%k

  CALL blacs_gemm('T','N',k,k,n,rone,qs_ot_env%matrix_sc0,matrix_x, &
                      rzero,qs_ot_env%matrix_buf1,globenv)

  CALL blacs_symm('L','U',k,k,rone,qs_ot_env%matrix_os,qs_ot_env%matrix_buf1, &
                      rzero,qs_ot_env%matrix_buf2,globenv)

  CALL blacs_gemm('N','N',n,k,k,-rone,qs_ot_env%matrix_sc0,qs_ot_env%matrix_buf2, &
                      rone,matrix_x,globenv)

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_ot_ortho_to_sc0


!*******************************************************************************
! c=c0*cos(p^0.5)+x*sin(p^0.5)*p^(-0.5)
! this assumes that x is already ortho to S*C0, and that p is x*p*x
!*******************************************************************************
SUBROUTINE qs_ot_get_orbitals(matrix_c, matrix_x, qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN) :: globenv
  TYPE(qs_ot_type)                          :: qs_ot_env
  TYPE(blacs_matrix_type), POINTER          :: matrix_c,matrix_x

! locals 

  REAL(wp) , PARAMETER            :: rone=1.0_wp , rzero=0.0_wp
  INTEGER                         :: n,k,handle

  CALL timeset("qs_ot_get_orbitals","I","",handle)

  n=qs_ot_env%n
  k=qs_ot_env%k

  CALL blacs_gemm('N','N',n,k,k,rone,qs_ot_env%matrix_c0,qs_ot_env%matrix_cosp, &
                                               rzero,matrix_c, globenv)
  CALL blacs_gemm('N','N',n,k,k,rone,matrix_x,qs_ot_env%matrix_sinp, &
                                               rone ,matrix_c, globenv)

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_ot_get_orbitals

!*******************************************************************************
! this routines computes dE/dx=dx, with dx ortho to sc0
! needs dE/dC=hc,C0,X,SX,p
!*******************************************************************************
SUBROUTINE qs_ot_get_derivative(matrix_hc,matrix_x,matrix_sx,matrix_dx, &
                                               qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN) :: globenv
  TYPE(blacs_matrix_type), POINTER          :: matrix_hc,matrix_x
  TYPE(blacs_matrix_type), POINTER          :: matrix_sx,matrix_dx
  TYPE(qs_ot_type)                          :: qs_ot_env

! locals 

  REAL(wp) , PARAMETER            :: rone=1.0_wp , rzero=0.0_wp
  INTEGER                         :: n,k,handle

  CALL timeset("qs_ot_get_derivative","I","",handle)

  n=qs_ot_env%n
  k=qs_ot_env%k

  ! this de/dc*(dX/dx)*sinp
  CALL blacs_gemm('N','N',n,k,k,rone,matrix_hc,qs_ot_env%matrix_sinp,rzero,matrix_dx, &
                                                          globenv)
  ! overlap hc*x
  CALL blacs_gemm('T','N',k,k,n,rone,matrix_hc,matrix_x,rzero,qs_ot_env%matrix_buf2, &
                                                        globenv)
  ! get it in the basis of the eigenvectors
  CALL blacs_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_buf2,qs_ot_env%matrix_r, &
                                      rzero,qs_ot_env%matrix_buf1, globenv)
  CALL blacs_gemm('T','N',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                      rzero,qs_ot_env%matrix_buf2, globenv)

  ! get the schur product of O_uv*B_uv
  CALL blacs_schur_product(qs_ot_env%matrix_buf2,qs_ot_env%matrix_sinp_b, &
                                   qs_ot_env%matrix_buf3,globenv)


  ! overlap hc*c0
  CALL blacs_gemm('T','N',k,k,n,rone,matrix_hc,qs_ot_env%matrix_c0,rzero, &
                                    qs_ot_env%matrix_buf2, globenv)
  ! get it in the basis of the eigenvectors
  CALL blacs_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_buf2,qs_ot_env%matrix_r, &
                                    rzero,qs_ot_env%matrix_buf1,globenv)
  CALL blacs_gemm('T','N',k,k,k,rone,qs_ot_env%matrix_r, qs_ot_env%matrix_buf1, &
                                    rzero,qs_ot_env%matrix_buf2, globenv)
  ! get the schur product of O_uv*B_uv
  CALL blacs_schur_product(qs_ot_env%matrix_buf2,qs_ot_env%matrix_cosp_b, &
                                        qs_ot_env%matrix_buf2, globenv)

  ! add the two bs and compute b+b^T

  CALL blacs_add(1.0_wp,qs_ot_env%matrix_buf3,1.0_wp,qs_ot_env%matrix_buf2, &
                                                                 globenv)

  ! get the b in the eigenvector basis
  CALL blacs_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_buf3,qs_ot_env%matrix_r, &
                                     rzero,qs_ot_env%matrix_buf1, globenv)
  CALL blacs_gemm('N','N',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                 rzero,qs_ot_env%matrix_buf3, globenv)


  CALL blacs_transpose(qs_ot_env%matrix_buf3,qs_ot_env%matrix_buf1,globenv)
  CALL blacs_add(1.0_wp,qs_ot_env%matrix_buf3,1.0_wp,qs_ot_env%matrix_buf1, &
                                                               globenv)

  ! and add to the derivative

  CALL blacs_gemm('N','N',n,k,k,rone,matrix_sx,qs_ot_env%matrix_buf3, &
                                          rone,matrix_dx, globenv)
 
  ! and make it orthogonal
  CALL qs_ot_ortho_to_sc0(matrix_dx,qs_ot_env,globenv)

  CALL timestop(0.0_wp,handle)

END SUBROUTINE qs_ot_get_derivative

!*******************************************************************************
! given p, computes  - eigenstuff (matrix_r,evals)
!                    - cos(p^0.5),p^(-0.5)*sin(p^0.5) 
!                    - the real b matrices, needed for the derivatives of these guys
!                    cosp_b_ij=(1/(2pii) * int(cos(z^1/2)/((z-eval(i))*(z-eval(j))))
!                    sinp_b_ij=(1/(2pii) * int(z^(-1/2)*sin(z^1/2)/((z-eval(i))*(z-eval(j))))
!*******************************************************************************
SUBROUTINE qs_ot_p2m(qs_ot_env,globenv)

  TYPE(global_environment_type), INTENT(IN)  :: globenv
  TYPE(qs_ot_type)                           :: qs_ot_env
  ! locals
  INTEGER  :: k,i,j
  REAL(wp) :: dum,a,b,timedum
  REAL(wp), PARAMETER :: rone=1.0_wp, rzero=0.0_wp


  k=qs_ot_env%k 
  CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_p,qs_ot_env%matrix_buf1)
  CALL blacs_syevd(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,qs_ot_env%evals, &
                                                     globenv)
  !CALL blacs_syevx(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,qs_ot_env%evals, &
  !                                                   k,1.0_wp,globenv)
  DO i=1,k
        qs_ot_env%evals(i)=MAX(0.0_wp,qs_ot_env%evals(i))
  ENDDO

  do i=1,k
     qs_ot_env%dum(i)=cos(sqrt(qs_ot_env%evals(i)))
  enddo
  CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_r,qs_ot_env%matrix_buf1)
  CALL blacs_scale_matrix(qs_ot_env%matrix_buf1,qs_ot_env%dum,globenv)
  CALL blacs_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                     rzero,qs_ot_env%matrix_cosp,globenv)

  do i=1,k
     qs_ot_env%dum(i)=qs_ot_sinc(sqrt(qs_ot_env%evals(i)))
  enddo
  CALL copy_blacs_to_blacs_matrix(qs_ot_env%matrix_r,qs_ot_env%matrix_buf1)
  CALL blacs_scale_matrix(qs_ot_env%matrix_buf1,qs_ot_env%dum,globenv)
  CALL blacs_gemm('N','T',k,k,k,rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
                                     rzero,qs_ot_env%matrix_sinp,globenv)

  do j=1,qs_ot_env%ncol_local
    do i=1,qs_ot_env%nrow_local
       a=(sqrt(qs_ot_env%evals(qs_ot_env%row_indices(i))) & 
                  -sqrt(qs_ot_env%evals(qs_ot_env%col_indices(j))))/2.0_wp
       b=(sqrt(qs_ot_env%evals(qs_ot_env%row_indices(i))) &
                  +sqrt(qs_ot_env%evals(qs_ot_env%col_indices(j))))/2.0_wp
       qs_ot_env%local_block(i,j)=-0.5_wp*qs_ot_sinc(a)*qs_ot_sinc(b)
    enddo
  enddo
  CALL blacs_set_local_block(qs_ot_env%matrix_cosp_b,qs_ot_env%local_block, &
                                          globenv)

  do j=1,qs_ot_env%ncol_local
    do i=1,qs_ot_env%nrow_local
       a=sqrt(qs_ot_env%evals(qs_ot_env%row_indices(i))) 
       b=sqrt(qs_ot_env%evals(qs_ot_env%col_indices(j)))
       qs_ot_env%local_block(i,j)=qs_ot_sincf(a,b)
    enddo
  enddo
  CALL blacs_set_local_block(qs_ot_env%matrix_sinp_b,qs_ot_env%local_block, &
                                          globenv)


END SUBROUTINE qs_ot_p2m
!*******************************************************************************


!*******************************************************************************
! computes sin(x)/x for all values of the argument
!*******************************************************************************
FUNCTION qs_ot_sinc(x)

	REAL(wp), INTENT(IN)               :: x
        REAL(wp)                           :: qs_ot_sinc
! locals
        REAL(wp)                           :: y
        REAL(wp), PARAMETER :: q1 =1.0_wp
        REAL(wp), PARAMETER :: q2 =-q1/(2.0_wp *3.0_wp)
        REAL(wp), PARAMETER :: q3 =-q2/(4.0_wp *5.0_wp)
        REAL(wp), PARAMETER :: q4 =-q3/(6.0_wp *7.0_wp)
        REAL(wp), PARAMETER :: q5 =-q4/(8.0_wp *9.0_wp)
        REAL(wp), PARAMETER :: q6 =-q5/(10.0_wp*11.0_wp)
        REAL(wp), PARAMETER :: q7 =-q6/(12.0_wp*13.0_wp)
        REAL(wp), PARAMETER :: q8 =-q7/(14.0_wp*15.0_wp)
        REAL(wp), PARAMETER :: q9 =-q8/(16.0_wp*17.0_wp)
        REAL(wp), PARAMETER :: q10=-q9/(18.0_wp*19.0_wp)

IF (abs(x)>0.5_wp) THEN
   qs_ot_sinc=sin(x)/x
ELSE
   y=x*x
   qs_ot_sinc=q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*(q9+y*(q10)))))))))
ENDIF
END FUNCTION qs_ot_sinc
!*******************************************************************************
! computes (1/(x^2-y^2))*(sinc(x)-sinc(y)) for all positive values of the arguments
!*******************************************************************************
FUNCTION qs_ot_sincf(xa,ya)

	REAL(wp), INTENT(IN)               :: xa,ya
        REAL(wp)                           :: qs_ot_sincf
 
! locals
        REAL(wp)                           :: dum,a,b,ybx,ybxs,sf,rs,xs,x,y
        INTEGER                            :: i

! this is currently a limit of the routine, could be removed rather easily
if (xa.lt.0) call stop_program("qs_ot_sincf","x is negative")
if (ya.lt.0) call stop_program("qs_ot_sincf","y is negative")

if (xa.lt.ya) then 
   x=ya
   y=xa
else
   x=xa
   y=ya
endif

if ( x .lt. 0.5_wp ) then ! use series, keeping in mind that x,y,x+y,x-y can all be zero

  qs_ot_sincf=0.0_wp
  if (x .gt. 0.0_wp) then
     ybx=y/x  
  else ! should be irrelevant  !?
     ybx=0.0_wp
  endif

  sf=-1.0_wp/((1.0_wp+ybx)*6.0_wp)
  rs=1.0_wp
  ybxs=ybx
  xs=1.0_wp

  do i=1,10
     qs_ot_sincf=qs_ot_sincf+sf*rs*xs*(1.0_wp+ybxs)
     sf=-sf/(real((2*i+2),wp)*real((2*i+3),wp))
     rs=rs+ybxs
     ybxs=ybxs*ybx
     xs=xs*x*x
  enddo

else ! no series expansion
  if ( x-y .gt. 0.1_wp ) then  ! safe to use the normal form
     qs_ot_sincf=(qs_ot_sinc(x)-qs_ot_sinc(y))/((x+y)*(x-y))
  else
     a=(x+y)/2.0_wp
     b=(x-y)/2.0_wp ! might be close to zero
                    ! y (=(a-b)) can not be close to zero since it is close to x>0.5
     qs_ot_sincf=(qs_ot_sinc(b)*cos(a)-qs_ot_sinc(a)*cos(b))/(2*x*y)
  endif
endif

END FUNCTION qs_ot_sincf

END MODULE qs_ot

