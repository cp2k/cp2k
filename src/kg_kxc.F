!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_force [1.0] *
!!
!!   NAME
!!     kg_kxc 
!!
!!   FUNCTION
!!     calculation of the Kinetic + Exchange correlation
!!
!!   AUTHOR
!!     gloria,30.09.2002 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE kg_kxc
  USE coefficient_types,               ONLY: coeff_transform_space,&
                                             coeff_type
  USE dft_types,                       ONLY: dft_control_type
  USE functionals
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: pi,&
                                             zero
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_max,&
                                             mp_range,&
                                             mp_sum,&
                                             mp_sync
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE pw_types,                        ONLY: pw_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_calculate_pw1
  USE xc_derivative_types,             ONLY: xc_derivative_set_type, &
                                             xc_derivative_set_create, &
                                             xc_derivative_set_release
  USE xc_derivatives,                  ONLY: xc_calc_derivatives

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calculate_kxc_potential
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='kg_kxc'

  INTERFACE calculate_kxc_potential
    MODULE PROCEDURE calc_kxc_pot_kg,calc_kxc_pot_kgpol 
  END INTERFACE

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** kgkxc/calc_kxc_pot_kg [1.1] *
!!
!!   NAME
!!     calc_kxc_pot_kg
!!
!!   FUNCTION
!!     Calculates the Kinetic-XC energy and potential        
!!
!!   AUTHOR
!!     gt         
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

 SUBROUTINE calc_kxc_pot_kg (pw_pool, rho_r, rho_g, v_r, e_xc, dft_control )

  IMPLICIT NONE

! Arguments
  TYPE ( pw_pool_type ), POINTER :: pw_pool
  TYPE ( coeff_type ), INTENT (in) :: rho_r
  TYPE ( coeff_type ), INTENT (in) :: rho_g
  TYPE ( coeff_type ), INTENT (out) :: v_r
  REAL ( dbl ), INTENT (out) :: e_xc
  TYPE ( dft_control_type ), INTENT ( IN ), TARGET :: dft_control

! Local parameters
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calc_kxc_pot_kg (MODULE kg_kxc)"
! Locals
  REAL (dbl), DIMENSION(:,:,:), POINTER :: rho
  REAL (dbl), DIMENSION(:,:,:), POINTER :: vxc
  COMPLEX (dbl), DIMENSION(:), POINTER :: rhog
!--------------------------------------------------------

  e_xc = 0.0_dbl
  rho => rho_r%pw%cr3d(:,:,:)
  vxc => v_r%pw%cr3d(:,:,:)
  rhog => rho_g%pw%cc(:)

  CALL xc_calculate_pw1 (rho_r%pw,rho_g%pw,v_r%pw, &
       e_xc,dft_control,pw_pool)
 
  vxc (:,:,:) =  vxc(:,:,:) *rho_r%pw%pw_grid%dvol

 END SUBROUTINE calc_kxc_pot_kg

!******************************************************************************
!!****** kgkxc/calc_kxc_pot_kgpol [1.1] *
!!
!!   NAME
!!     calc_kxc_pot_kgpol
!!
!!   FUNCTION
!!     Calculates the K-XC energy  and potential, and potentials associated
!!     to  second and third derivatives of the KXC kernel 
!!     V1_r = dvol* rhop_r * (partial^2  (E_kxc) / partial^2 rho0)
!!     V2_r = dvol**2 * rhop_r**2 * (partial ^3 (E_kxc) /partial^3 rho0 )
!!     Full gradient functional for Exc and Vxc, only lda functional
!!     for higher order derivatives.
!!
!!   AUTHOR
!!     gt         
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

 SUBROUTINE calc_kxc_pot_kgpol ( pw_pool,rho_r, rho_g, rhop_r, v_r, v1_r, &
                           v2_r, e_xc, e1_xc, e2_xc, dft_control )

  IMPLICIT NONE

! Arguments
  TYPE (pw_pool_type), POINTER :: pw_pool
  TYPE ( coeff_type ), INTENT (in) :: rho_r
  TYPE ( coeff_type ), INTENT (in) :: rho_g
  TYPE ( coeff_type ), INTENT (in) :: rhop_r
  TYPE ( coeff_type ), INTENT (out) :: v_r
  TYPE ( coeff_type ), INTENT (out) :: v1_r
  TYPE ( coeff_type ), INTENT (out) :: v2_r
  REAL ( dbl ), INTENT (out) :: e_xc
  REAL ( dbl ), INTENT (out) :: e1_xc
  REAL ( dbl ), INTENT (out) :: e2_xc
  TYPE ( dft_control_type ), INTENT ( IN ), TARGET :: dft_control

! Local parameters
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calc_kxc_pot_kgpol (MODULE kg_kxc)"
! Locals
  REAL (dbl), DIMENSION(:,:,:), POINTER :: rho
  REAL (dbl), DIMENSION(:,:,:), POINTER :: vxc
  REAL (dbl), DIMENSION(:,:,:), POINTER :: v1xc
  REAL (dbl), DIMENSION(:,:,:), POINTER :: v2xc
  COMPLEX (dbl), DIMENSION(:), POINTER :: rhog
!--------------------------------------------------------

  e_xc = 0.0_dbl
  vxc => v_r%pw%cr3d(:,:,:)
  v1xc => v1_r%pw%cr3d(:,:,:)
  v2xc => v2_r%pw%cr3d(:,:,:)

  CALL xc_calculate_pw1 (rho_r%pw,rho_g%pw,v_r%pw, &
       e_xc,dft_control,pw_pool)

  CALL xc_calculate_derivatives (rho_r%pw,rhop_r%pw,v_r%pw,v1_r%pw, &
                             v2_r%pw, e1_xc, e2_xc, dft_control)
 
  e1_xc =  rho_r%pw%pw_grid%dvol*e1_xc      
  e2_xc =  0.5_dbl*rho_r%pw%pw_grid%dvol*e2_xc      
  vxc (:,:,:) =   vxc(:,:,:) *rho_r%pw%pw_grid%dvol
  v1xc (:,:,:) =  v1xc(:,:,:) *rho_r%pw%pw_grid%dvol
  v2xc (:,:,:) =  0.5_dbl * v2xc(:,:,:) *rho_r%pw%pw_grid%dvol

 END SUBROUTINE calc_kxc_pot_kgpol
!******************************************************************************
!!****** xc/xc_calculate_derivatives [1.0] *
!!
!!   NAME
!!     xc_calculate_derivatives
!!
!!   FUNCTION
!!     Exchange and Correlation functional derivatices
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE xc_calculate_derivatives ( rho_r,rhop_r,v,v1,v2,e1,e2,dft_control )

  IMPLICIT NONE

! Arguments
  TYPE(pw_type), INTENT(IN) :: rho_r
  TYPE(pw_type), INTENT(IN) :: rhop_r
  TYPE(pw_type), INTENT(IN) :: v
  TYPE(pw_type), INTENT(INOUT) :: v1
  TYPE(pw_type), INTENT(INOUT) :: v2
  REAL (dbl), INTENT (out) :: e1
  REAL (dbl), INTENT (out) :: e2
  TYPE(dft_control_type), INTENT(IN) :: dft_control

! Local
  REAL ( dbl ) , DIMENSION (:), ALLOCATABLE :: rgrid
  LOGICAL, DIMENSION  (3) :: gradient_functionals, crossterms
  INTEGER :: stat, ntot, n(3), i, j, k, ii, order, ik, npot
  INTEGER, DIMENSION (3,3) :: nd
  INTEGER, DIMENSION (2,3) :: bo
  LOGICAL :: failure, spezgrid, do_lsd
  TYPE(xc_derivative_set_type), pointer :: derivative_set

! Parameters
  CHARACTER(len=*), PARAMETER :: routineN='xc_calculate_derivatives',&
                                 routineP=moduleN//':'//routineN
!------------------------------------------------------------------------------
  failure=.FALSE.
  nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))
  e1 = 0.0_dbl
  e2 = 0.0_dbl
  ! is a special grid needed
  spezgrid = .NOT. dft_control%integration_grid%pw_default
  ! transform to higher cutoff
  IF ( spezgrid ) THEN
     CALL stop_program ( routineP, "NYI" )
  ELSE
     ! allocate and copy to grids
     n(1:3) = rho_r%pw_grid%npts_local (1:3)
     ntot = PRODUCT ( n )
     ALLOCATE ( rgrid ( ntot ), STAT = stat )
     IF ( stat /= 0 ) CALL stop_memory ( routineP, "rgrid", ntot )
     bo = rho_r%pw_grid%bounds_local
     ii = 0
     DO k = bo(1,3), bo(2,3)
        DO j = bo(1,2), bo(2,2)
           DO i = bo(1,1), bo(2,1)
              ii = ii + 1
              rgrid(ii) = rho_r%cr3d(i,j,k)
           END DO
        END DO
     END DO
  END IF

  order = 3

! no gradient corrections are implemented for derivatives
  gradient_functionals = dft_control%xc_control%gradient_functionals
  do_lsd = (dft_control%nspins > 1)
  crossterms = (dft_control%xc_control%crossterms.AND.do_lsd)
  IF (ANY(gradient_functionals).OR.ANY(crossterms).OR.do_lsd) THEN
     CALL stop_program(routineP, "not yet implemented")
  END IF

  CALL xc_derivative_set_create(derivative_set, dft_control%xc_control, ntot)
  CALL xc_calc_derivatives(derivative_set, order, rgrid)
  ! copy back to structure
  IF ( spezgrid ) THEN
     ! transform to lower cutoff
     CALL stop_program ( routineP, "NYI" )
  ELSE
     ii = 0
     DO k = bo(1,3), bo(2,3)
        DO j = bo(1,2), bo(2,2)
           DO i = bo(1,1), bo(2,1)
              ii = ii + 1
              v1%cr3d(i,j,k) = derivative_set%data(ii,2)*rhop_r%cr3d(i,j,k)
              v2%cr3d(i,j,k) = derivative_set%data(ii,3)*rhop_r%cr3d(i,j,k)
           END DO
        END DO
     END DO
     DO k = bo(1,3), bo(2,3)
        DO j = bo(1,2), bo(2,2)
           DO i = bo(1,1), bo(2,1)
              v2%cr3d(i,j,k) =v2%cr3d(i,j,k)*rhop_r%cr3d(i,j,k)
              e1= e1+ v%cr3d(i,j,k)*rhop_r%cr3d(i,j,k)
              e2= e2+ v1%cr3d(i,j,k)*rhop_r%cr3d(i,j,k)
           END DO
        END DO
     END DO
     DEALLOCATE ( rgrid, STAT = stat )
     IF ( stat /= 0 ) CALL stop_memory ( routineP, "rgrid" )
  END IF
  CALL xc_derivative_set_release(derivative_set)

END SUBROUTINE xc_calculate_derivatives

!******************************************************************************

END MODULE kg_kxc  

!******************************************************************************
