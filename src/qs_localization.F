!   CP2K: A general program to perform molecular dynamics simulations         !

!!****** cp2k/qs_localization [1.0] *
!!
!!   NAME
!!     qs_localization
!!
!!   FUNCTION
!! 
!!     should be able to find properties related to the position operator in PBC
!!     e.g. location and spread of periodic orbitals
!!          Warnier functions
!!          ...
!!
!!   AUTHORS
!!     Leonardo Bernasconi (01.2003)
!!     Joost VandeVondele (01.2003)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_localization
! *****************************************************************************
  USE coefficient_types,               ONLY: coeff_type, &
                                             coeff_zero
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_full_matrix_type, &
                                             cp_full_matrix_p_type, &
                                             cp_fm_create, &
                                             cp_fm_release, &
                                             cp_fm_get_element
  USE kinds,                           ONLY: wp => dp
  USE cp_para_env,                     ONLY: cp_para_env_retain, &
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_fm_basic_linalg,              ONLY: cp_sm_fm_multiply, &
                                             cp_fm_gemm
  USE mathconstants,                   ONLY: twopi
  USE physcon,                         ONLY: angstrom
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE sparse_matrix_types,             ONLY: real_matrix_type, &
                                             real_matrix_p_type, &
                                             replicate_matrix_structure, &
                                             deallocate_matrix, &
                                             set_matrix
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_types,                        ONLY: REALDATA3D, &
                                             REALSPACE
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE qs_blacs,                        ONLY: write_blacs_matrix
  IMPLICIT NONE

! *** Public ***
  PUBLIC :: qs_loc_init
  PUBLIC :: qs_loc_env_type
  PUBLIC :: qs_loc_compute_set
  PUBLIC :: qs_loc_destroy

! *** Localization Data Type ***
  TYPE qs_loc_env_type
       INTEGER :: is_initialized
       TYPE(real_matrix_p_type), DIMENSION(2,3) :: pos_op
       TYPE(cp_para_env_type), POINTER :: para_env
       REAL(wp) :: box_L(3)
  END TYPE qs_loc_env_type

! *****************************************************************************

CONTAINS

!******************************************************************************
!
! allocates the data, and initializes the operators
!
!******************************************************************************
SUBROUTINE qs_loc_init(qs_loc_env,qs_env)

  TYPE(qs_loc_env_type), POINTER :: qs_loc_env
  TYPE(qs_environment_type) :: qs_env
  TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: s
  INTEGER :: i,j

  TYPE(coeff_type) :: op_rspace
  TYPE(pw_env_type), POINTER :: pw_env
  TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
  INTEGER :: handle
  INTEGER ix,iy,iz
  REAL(wp) :: x,y,z,thearg,theval

  CALL timeset("qs_loc_init","I","",handle)
  NULLIFY(pw_env,auxbas_pw_pool)
  NULLIFY(qs_loc_env%para_env)
  !
  ! start replicating the s structure, we do have 6 operators (x,y,z)*(real,imag)
  ! notice that the pos_op should have the same structure as S
  ! borrow para env from s
  !
  CALL get_qs_env(qs_env=qs_env,s=s,para_env=qs_loc_env%para_env)
  CALL cp_para_env_retain(qs_loc_env%para_env) 
  DO i=1,3
     DO j=1,2
       NULLIFY(qs_loc_env%pos_op(j,i)%matrix)
       CALL replicate_matrix_structure(s(1)%matrix, &
                      qs_loc_env%pos_op(j,i)%matrix,"pos op")
       CALL set_matrix(qs_loc_env%pos_op(j,i)%matrix,0.0_wp)
     ENDDO
  END DO
  !
  ! OK, now we construct these operators the easy way, put them on the pw grid and integrate
  ! analytical construction should be faster.
  !

  ! get a pw grid
  CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
  CALL pw_env_get(pw_env=pw_env,auxbas_pw_pool=auxbas_pw_pool)
  CALL pw_pool_init_coeff(auxbas_pw_pool,op_rspace,&
                          use_data=REALDATA3D,in_space=REALSPACE)
  
  DO i=1,3
   qs_loc_env%box_L(i)=op_rspace%pw%pw_grid%npts(i)*op_rspace%pw%pw_grid%dr(i)
   DO j=1,2
      ! CALL coeff_zero(op_rspace) 
      ! the definition of the position on the PW grid has to be 
      ! compatible with the one used in the integration procedure / collocation procedure
      ! integrate_v_rspace calculate_rho_elec
      do iz=op_rspace%pw%pw_grid%bounds_local(1,3),op_rspace%pw%pw_grid%bounds_local(2,3)
      do iy=op_rspace%pw%pw_grid%bounds_local(1,2),op_rspace%pw%pw_grid%bounds_local(2,2)
      do ix=op_rspace%pw%pw_grid%bounds_local(1,1),op_rspace%pw%pw_grid%bounds_local(2,1)
         x=dble(ix-op_rspace%pw%pw_grid%bounds(1,1))*twopi/dble(op_rspace%pw%pw_grid%npts(1))
         y=dble(iy-op_rspace%pw%pw_grid%bounds(1,2))*twopi/dble(op_rspace%pw%pw_grid%npts(2))
         z=dble(iz-op_rspace%pw%pw_grid%bounds(1,3))*twopi/dble(op_rspace%pw%pw_grid%npts(3)) 
         SELECT CASE (i)
         CASE (1)
            thearg=x
         CASE (2)
            thearg=y
         CASE (3)
            thearg=z
         END SELECT
         SELECT CASE (j)
         CASE (1)
            theval=cos(thearg)
         CASE (2)
            theval=sin(thearg)
         END SELECT
         ! notice the volume element, slightly ugly I think
         op_rspace%pw%cr3d(ix,iy,iz)=theval*op_rspace%pw%pw_grid%dvol 
      enddo
      enddo
      enddo
      ! do the actual integration
      CALL integrate_v_rspace(op_rspace,qs_loc_env%pos_op(j,i),qs_env)
   ENDDO
  ENDDO

  ! done
  CALL pw_pool_give_back_coeff(auxbas_pw_pool, op_rspace)
  qs_loc_env%is_initialized = 777
  CALL timestop(0.0_wp,handle) 
END SUBROUTINE qs_loc_init


!******************************************************************************
!
! applies the operators to obtain the necessary data
!
!******************************************************************************
SUBROUTINE qs_loc_compute_set(qs_loc_env,vectors,ionode)

  TYPE(qs_loc_env_type), POINTER :: qs_loc_env
  TYPE(cp_full_matrix_type), POINTER :: vectors
  LOGICAL :: ionode

  TYPE(cp_full_matrix_type), POINTER :: opvec
  TYPE(cp_full_matrix_p_type), DIMENSION(2,3) :: result
  TYPE(cp_fm_struct_type), POINTER :: tmp_fm_struct

  INTEGER :: n,k
  INTEGER :: i,j,l
  INTEGER :: handle

  COMPLEX(kind=wp) :: z
  REAL(wp) :: realpart,imagpart,pos(3),width(3)

  CALL timeset("qs_loc_compute_set","I","",handle)
  NULLIFY(tmp_fm_struct)

  ! get rows and cols of the input
  CALL cp_fm_get_info(vectors,nrow_global=n,ncol_global=k)
  ! replicate the input kind of matrix
  CALL cp_fm_create(opvec,vectors%matrix_struct)
  ! create a new matrix structure k x k in size, borrowing most of info
  CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k, &
            ncol_global=k,para_env=vectors%matrix_struct%para_env, &
            context=vectors%matrix_struct%context)
  write(6,*) "nrow_local",tmp_fm_struct%nrow_local,"ncol_local",tmp_fm_struct%ncol_local,&
             "nrow_block",tmp_fm_struct%nrow_block,"ncol_block",tmp_fm_struct%ncol_block
  ! make the actual matrices
  DO i=1,3
    DO j=1,2
        NULLIFY(result(j,i)%matrix)
        CALL cp_fm_create(result(j,i)%matrix,tmp_fm_struct)
    ENDDO
  ENDDO
  ! release one instance of the fm_struct
  CALL cp_fm_struct_release(tmp_fm_struct) 

  ! Compute results here
  DO i=1,3
    DO j=1,2
       CALL cp_sm_fm_multiply(qs_loc_env%pos_op(j,i)%matrix, &
                              vectors,opvec,ncol=k,para_env=qs_loc_env%para_env)
       CALL cp_fm_gemm("T","N",k,k,n,1.0_wp,vectors,opvec,0.0_wp,result(j,i)%matrix)
       ! CALL write_blacs_matrix(result(j,i)%matrix,6,formatted=.true.)
    ENDDO
  ENDDO

  ! OK print out the center of these friends
  IF (ionode) OPEN(UNIT=17,FILE="WAVEFUNCTION_CENTERS.xyz")
  DO l=1,k
    DO i=1,3
       CALL cp_fm_get_element(result(1,i)%matrix,l,l,realpart)
       CALL cp_fm_get_element(result(2,i)%matrix,l,l,imagpart)
       z=cmplx(realpart,imagpart)
       pos(i)=(angstrom*qs_loc_env%box_L(i)/(twopi))*IMAG(log(z))
       if (pos(i).lt.0.0_wp) pos(i)=pos(i)+qs_loc_env%box_L(i)*angstrom        ! :-(
       width(i)=(angstrom*qs_loc_env%box_L(i)/(twopi))*sqrt(-2.0_wp*log(ABS(Z)))    ! sqrt of eq. 7 
    ENDDO
    IF (ionode) write(17,'(A,6F12.8)') "He",pos(1:3),width(1:3)
  ENDDO
  IF (ionode) CLOSE(17)

  ! release results
  DO i=1,3
    DO j=1,2
       CALL cp_fm_release(result(j,i)%matrix)
    ENDDO
  ENDDO
  CALL timestop(0.0_wp,handle) 

END SUBROUTINE qs_loc_compute_set
!******************************************************************************
!
! deallocates all data 
!
!******************************************************************************
SUBROUTINE qs_loc_destroy(qs_loc_env)

  TYPE(qs_loc_env_type), POINTER :: qs_loc_env

  INTEGER :: i , j

  if (qs_loc_env%is_initialized.ne.777) call stop_program("qs_localization","problems")

  DO i=1,3
   DO j=1,2
     CALL deallocate_matrix(qs_loc_env%pos_op(j,i)%matrix)
   ENDDO
  END DO
  CALL cp_para_env_release(qs_loc_env%para_env) 

END SUBROUTINE qs_loc_destroy

END MODULE qs_localization

