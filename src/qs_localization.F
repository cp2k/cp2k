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
!!          Wannier functions
!!
!!   AUTHORS
!!     Leonardo Bernasconi (01.2003)
!!     Joost VandeVondele (01.2003)
!!     CJM (04.2003)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_localization
! *****************************************************************************
  USE coefficient_types,               ONLY: coeff_type,&
                                             coeff_zero
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_sm_fm_multiply, &
                                             cp_fm_lu_decompose
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_lu_decompose
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_complex_full_matrix_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: wp => dp
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_sum
  USE physcon,                         ONLY: angstrom
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE
  USE qs_blacs,                        ONLY: write_blacs_matrix
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_localization_methods,         ONLY: jacobi_rotations, initialize_weights
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

! *** Public ***
  PUBLIC :: qs_loc_init
  PUBLIC :: qs_loc_env_type, PIPEK_MENZY, POSITION_OPERATOR
  PUBLIC :: qs_loc_compute_properties
  PUBLIC :: qs_loc_destroy

! *** Localization Data Type ***
  TYPE qs_loc_env_type
       INTEGER :: operator_type
       LOGICAL :: localization
       INTEGER :: is_initialized
       TYPE(real_matrix_p_type),POINTER :: operator(:,:)
       TYPE(cp_para_env_type), POINTER :: para_env
       TYPE ( cell_type ), POINTER :: cell
       REAL ( wp ) :: weights ( 3 )
  END TYPE qs_loc_env_type
  
  INTEGER, PARAMETER :: PIPEK_MENZY=301, POSITION_OPERATOR=302

! *****************************************************************************

CONTAINS

!******************************************************************************
!
! allocates the data, and initializes the operators
!
!******************************************************************************
SUBROUTINE qs_loc_init(qs_loc_env,qs_env)
  IMPLICIT NONE

  TYPE(qs_loc_env_type), POINTER :: qs_loc_env
  TYPE(qs_environment_type) :: qs_env
! locals
  TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: s
  INTEGER :: i, j, natoms, isos, handle

  CALL timeset("qs_loc_init","I","",handle)

  NULLIFY(qs_loc_env%para_env)
  !
  ! start replicating the s structure, we do have
  ! 6 operators (x,y,z)*(real,imag) for POSITION_OPERATOR
  ! and we have NATOMS operators for PIPEK_MENZY
  ! notice that the operator should have the same structure as S
  ! borrow para env from S
  !
  CALL get_qs_env(qs_env=qs_env,s=s,cell=qs_loc_env%cell,para_env=qs_loc_env%para_env)
  CALL cp_para_env_retain(qs_loc_env%para_env) 

! dbg
  qs_loc_env % operator_type = POSITION_OPERATOR
  qs_loc_env % localization = .TRUE.
! dbg
  IF ( qs_loc_env % operator_type == POSITION_OPERATOR )  THEN
    ALLOCATE ( qs_loc_env % operator ( 2, 3 ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'operator' )
  ELSEIF ( qs_loc_env % operator_type == PIPEK_MENZY ) THEN
    natoms = SIZE ( qs_env % particle_set )
    ALLOCATE ( qs_loc_env % operator ( 1, natoms ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'operator' )
  ENDIF

  DO i=1, SIZE ( qs_loc_env % operator, 2 )
     DO j=1, SIZE ( qs_loc_env % operator, 1 )
       NULLIFY(qs_loc_env%operator(j,i)%matrix)
       CALL replicate_matrix_structure(s(1)%matrix, &
            qs_loc_env%operator(j,i)%matrix,"operator")
       CALL set_matrix(qs_loc_env%operator(j,i)%matrix,0.0_wp)
     ENDDO
  END DO


  IF ( qs_loc_env % operator_type == POSITION_OPERATOR ) THEN 
    CALL initialize_weights ( qs_loc_env % cell, qs_loc_env % weights ) 
    CALL get_boys_foster_operator ( qs_loc_env, qs_env )
  ENDIF

  IF ( qs_loc_env % operator_type == PIPEK_MENZY ) & 
  CALL get_pipek_menzy_operator ( qs_loc_env, qs_env )
  

! done
  qs_loc_env%is_initialized = 777
  CALL timestop(0.0_wp,handle) 
END SUBROUTINE qs_loc_init

!******************************************************************************
!
! Computes the Boys-Foster operator for periodic systems
!
!******************************************************************************
SUBROUTINE get_boys_foster_operator ( qs_loc_env, qs_env ) 

! globals
  TYPE(qs_loc_env_type), POINTER :: qs_loc_env
  TYPE(qs_environment_type) :: qs_env
! locals
  TYPE(coeff_type) :: op_rspace
  TYPE(pw_env_type), POINTER :: pw_env
  TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
  INTEGER ix, iy, iz, i, j
  REAL(wp) :: x, y, z, thearg, theval

  NULLIFY ( pw_env, auxbas_pw_pool )
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
   DO j=1,2
      ! CALL coeff_zero(op_rspace) 
      ! the definition of the position on the PW grid has to be 
      ! compatible with the one used in the integration procedure / collocation procedure
      ! integrate_v_rspace calculate_rho_elec
      DO iz=op_rspace%pw%pw_grid%bounds_local(1,3),op_rspace%pw%pw_grid%bounds_local(2,3)
      DO iy=op_rspace%pw%pw_grid%bounds_local(1,2),op_rspace%pw%pw_grid%bounds_local(2,2)
      DO ix=op_rspace%pw%pw_grid%bounds_local(1,1),op_rspace%pw%pw_grid%bounds_local(2,1)
         x=REAL(ix-op_rspace%pw%pw_grid%bounds(1,1),wp)*twopi/REAL(op_rspace%pw%pw_grid%npts(1),wp)
         y=REAL(iy-op_rspace%pw%pw_grid%bounds(1,2),wp)*twopi/REAL(op_rspace%pw%pw_grid%npts(2),wp)
         z=REAL(iz-op_rspace%pw%pw_grid%bounds(1,3),wp)*twopi/REAL(op_rspace%pw%pw_grid%npts(3),wp) 
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
            theval=COS(thearg)
         CASE (2)
            theval=SIN(thearg)
         END SELECT
         ! notice the volume element, slightly ugly I think
         op_rspace%pw%cr3d(ix,iy,iz)=theval*op_rspace%pw%pw_grid%dvol 
      ENDDO
      ENDDO
      ENDDO
! do the actual integration
      CALL integrate_v_rspace(op_rspace,qs_loc_env%operator(j,i),qs_env)
   ENDDO
  ENDDO

  CALL pw_pool_give_back_coeff ( auxbas_pw_pool, op_rspace )

END SUBROUTINE get_boys_foster_operator

!******************************************************************************
!
! Computes the Pipek-Menzy operator for periodic systems
!
!******************************************************************************
SUBROUTINE get_pipek_menzy_operator ( qs_loc_env, qs_env )
  IMPLICIT NONE
! globals
  TYPE(qs_loc_env_type), POINTER :: qs_loc_env
  TYPE(qs_environment_type) :: qs_env
  CALL stop_program("qs_localization","PIPEK MENZY not implemented")
END SUBROUTINE get_pipek_menzy_operator

!******************************************************************************
!
! applies the operators to obtain the necessary data
!
!******************************************************************************
SUBROUTINE qs_loc_compute_properties(qs_loc_env,vectors,ionode,iunit)

  TYPE(qs_loc_env_type), POINTER :: qs_loc_env
  TYPE(cp_full_matrix_type), POINTER :: vectors
  LOGICAL :: ionode
  INTEGER :: iunit
! locals
  TYPE(cp_full_matrix_type), POINTER :: opvec
  TYPE(cp_full_matrix_p_type), ALLOCATABLE :: ZIJ ( :, : )
  TYPE(cp_fm_struct_type), POINTER :: tmp_fm_struct

  INTEGER :: n,k,natoms
  INTEGER :: i,j,l,isos
  INTEGER :: handle

  CALL timeset("qs_loc_compute_properties","I","",handle)
  NULLIFY(tmp_fm_struct,opvec)

  ! get rows and cols of the input
  CALL cp_fm_get_info(vectors,nrow_global=n,ncol_global=k)
  ! replicate the input kind of matrix
  CALL cp_fm_create(opvec,vectors%matrix_struct)
  ! create a new matrix structure k x k in size, borrowing most of info
  CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k, &
            ncol_global=k,para_env=vectors%matrix_struct%para_env, &
            context=vectors%matrix_struct%context)

! Allocating ZIJ
  IF ( qs_loc_env % operator_type == POSITION_OPERATOR ) THEN
    ALLOCATE ( ZIJ ( 2, 3 ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'ZIJ' )
  ELSEIF ( qs_loc_env % operator_type == PIPEK_MENZY ) THEN
    natoms = SIZE ( qs_loc_env % operator, 2 ) 
    ALLOCATE ( ZIJ ( 1, natoms ), stat = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'ZIJ' )
  ENDIF
! make the actual matrices
  DO i=1,SIZE ( ZIJ, 2 )
    DO j=1, SIZE ( ZIJ, 1 )
        NULLIFY ( ZIJ ( j, i ) % matrix )
        CALL cp_fm_create ( ZIJ ( j, i ) % matrix, tmp_fm_struct )
    ENDDO
  ENDDO

! release one instance of the fm_struct
  CALL cp_fm_struct_release ( tmp_fm_struct ) 

! Compute zij here
  DO i=1,SIZE ( ZIJ, 2 )
    DO j=1,SIZE ( ZIJ, 1 )
       CALL cp_sm_fm_multiply(qs_loc_env%operator(j,i)%matrix, &
                            vectors,opvec,ncol=k,para_env=qs_loc_env%para_env)
       CALL cp_fm_gemm("T","N",k,k,n,1.0_wp,vectors,opvec,0.0_wp,ZIJ(j,i)%matrix)
!       CALL write_blacs_matrix(ZIJ(j,i)%matrix,6,formatted=.true.)
    ENDDO
  ENDDO

  IF ( qs_loc_env % localization ) CALL jacobi_rotations ( qs_loc_env % weights, zij, vectors ) 

  IF ( qs_loc_env % operator_type == POSITION_OPERATOR )  & 
  CALL centers_and_spreads ( ionode, iunit, qs_loc_env, zij ) 


! release results

  DO i=1,SIZE ( ZIJ, 2 )
    DO j=1, SIZE ( ZIJ, 1 )
       CALL cp_fm_release(ZIJ(j,i)%matrix)
    ENDDO
  ENDDO

  CALL timestop(0.0_wp,handle) 

END SUBROUTINE qs_loc_compute_properties

!******************************************************************************
!
! computes spread and centers
!
!******************************************************************************

SUBROUTINE centers_and_spreads ( ionode, iunit, qs_loc_env, zij )
! globals
  INTEGER, INTENT ( IN ) :: iunit
  LOGICAL, INTENT ( IN ) :: ionode
  TYPE ( cp_full_matrix_p_type ), INTENT ( INOUT ) :: zij ( :, : )
  TYPE(qs_loc_env_type), POINTER :: qs_loc_env

! locals
  TYPE ( cell_type ), POINTER :: cell
  TYPE(cp_complex_full_matrix_type), POINTER :: complexresult
  COMPLEX ( KIND = wp ) :: almost_determinant
  COMPLEX ( KIND = wp ) :: z
  REAL ( KIND = wp ) :: realpart, imagpart, center ( 3 ), spread ( 3 )
  INTEGER :: l, i, nstates

! get nstates
  CALL cp_fm_get_info ( zij ( 1, 1 ) % matrix, nrow_global = nstates )

! initialize complexresult for boys-foster
  NULLIFY ( complexresult )
  CALL cp_cfm_create ( complexresult, zij ( 1, 1 ) % matrix % matrix_struct )
!
! OK print out the center of these friends
  cell => qs_loc_env % cell
  DO l=1,nstates
    DO i=1,3
! only the cpu that has the data can return it. others give 0
       CALL cp_fm_get_element(zij(1,i)%matrix,l,l,realpart)
       CALL cp_fm_get_element(zij(2,i)%matrix,l,l,imagpart)
       z=CMPLX(realpart,imagpart)
       center ( i )=( cell % hmat ( i, i ) / twopi ) * AIMAG ( LOG ( z ) )
       spread ( i )=SQRT( -qs_loc_env % weights ( i ) * 2.0_wp * LOG ( ABS ( Z ) ) ) / twopi
    ENDDO
    IF (ionode) WRITE(iunit,'(A,6F12.8)') "He", angstrom*center( 1 : 3 ),  &
                                                angstrom*spread( 1 : 3 )
  ENDDO
!
  ! now, find the spread of the orbitals
  ZIJ(1,1)%matrix%local_data(1,1)=0.0_wp
  DO i=1,3
     complexresult%local_data ( :, : ) = CMPLX (ZIJ(1,i)%matrix%local_data(:,:), &
                                         ZIJ(2,i)%matrix%local_data(:,:))
     CALL cp_cfm_lu_decompose ( complexresult, almost_determinant )
     IF (ionode) WRITE(6,*) "D Eq. 11 ",i,SQRT(-2.0*nstates*LOG(ABS(almost_determinant)))
  ENDDO

  CALL cp_cfm_release(complexresult)

END SUBROUTINE centers_and_spreads

!******************************************************************************
!
! deallocates all data 
!
!******************************************************************************
SUBROUTINE qs_loc_destroy ( qs_loc_env )

  TYPE(qs_loc_env_type), POINTER :: qs_loc_env

  INTEGER :: i , j

  IF (qs_loc_env%is_initialized.ne.777) CALL stop_program("qs_localization","problems")

  DO i=1,SIZE ( qs_loc_env % operator, 2 )
   DO j=1,SIZE ( qs_loc_env % operator, 1 )
     CALL deallocate_matrix(qs_loc_env%operator(j,i)%matrix)
   ENDDO
  END DO
  CALL cp_para_env_release(qs_loc_env%para_env) 

END SUBROUTINE qs_loc_destroy

END MODULE qs_localization

