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
  USE coefficient_types,               ONLY: coeff_type
  USE cp_control_types,                ONLY: dft_control_type
  USE functionals
  USE kinds,                           ONLY: dbl
  USE message_passing,                 ONLY: mp_sum
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE pw_types,                        ONLY: pw_type
  USE pw_grid_types,                   ONLY: PW_MODE_DISTRIBUTED
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_release
  USE xc_derivatives,                  ONLY: xc_calc_derivatives

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calculate_kxc_derivatives
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='kg_kxc'

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** kgkxc/calculate_kxc_derivatives [1.1] *
!!
!!   NAME
!!     calculate_kxc_derivatives
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

 SUBROUTINE calculate_kxc_derivatives ( pw_pool,rho_r, rho_g, rhop_r, v_r, v1_r, &
                           v2_r, e1_xc, e2_xc, dft_control,calculate_forces )


    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(coeff_type), INTENT(in)             :: rho_r, rho_g, rhop_r, v_r
    TYPE(coeff_type), INTENT(inout)          :: v1_r, v2_r
    REAL(dbl), INTENT(out)                   :: e1_xc, e2_xc
    TYPE(dft_control_type), INTENT(IN), &
      TARGET                                 :: dft_control
    LOGICAL, INTENT(in)                      :: calculate_forces

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calc_kxc_derivatives (MODULE kg_kxc)"

    INTEGER                                  :: handle, order
    REAL(dbl), DIMENSION(:, :, :), POINTER   :: v1xc, v2xc

!--------------------------------------------------------
  CALL timeset ( 'routine','I',' ',handle )

  v1xc => v1_r%pw%cr3d(:,:,:)
  v2xc => v2_r%pw%cr3d(:,:,:)

  IF(calculate_forces) THEN
    order=3
  ELSE
    order=2
  END IF

  CALL xc_calculate_derivatives (rho_r%pw,rhop_r%pw,v_r%pw,v1_r%pw, &
                             v2_r%pw, e1_xc, e2_xc, order, dft_control, pw_pool)
 
  IF ( rho_r%pw%pw_grid%para%mode == PW_MODE_DISTRIBUTED ) THEN
    CALL mp_sum ( e1_xc, rho_r%pw%pw_grid%para%group )
    CALL mp_sum ( e2_xc, rho_r%pw%pw_grid%para%group )
  END IF

  e2_xc =  0.5_dbl*rho_r%pw%pw_grid%dvol*e2_xc      

  v1xc (:,:,:) =  v1xc(:,:,:) *rho_r%pw%pw_grid%dvol
  IF (calculate_forces)THEN
    v2xc (:,:,:) =  0.5_dbl * v2xc(:,:,:) *rho_r%pw%pw_grid%dvol
  END IF
  CALL timestop(0.0_dbl,handle)
 END SUBROUTINE calculate_kxc_derivatives
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
SUBROUTINE xc_calculate_derivatives ( rho_r,rhop_r,v,v1,v2,e1,e2,order,&
                                      dft_control, pw_pool )


    TYPE(pw_type), INTENT(IN)                :: rho_r, rhop_r, v
    TYPE(pw_type), INTENT(INOUT)             :: v1, v2
    REAL(dbl), INTENT(out)                   :: e1, e2
    INTEGER, intent(IN)                      :: order
    TYPE(dft_control_type), INTENT(IN)       :: dft_control
    TYPE(pw_pool_type), POINTER              :: pw_pool

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_calculate_derivatives', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, j, k, n(3), ntot, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3, 3)                 :: nd
    LOGICAL                                  :: do_lsd, failure, spezgrid
    LOGICAL, DIMENSION(3)                    :: crossterms, &
                                                gradient_functionals
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: rgrid
    TYPE(xc_derivative_set_type), POINTER    :: derivative_set

!------------------------------------------------------------------------------

  failure=.FALSE.

  NULLIFY(derivative_set)
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

! no gradient corrections are implemented for derivatives
  gradient_functionals = dft_control%xc_control%gradient_functionals
  do_lsd = (dft_control%nspins > 1)
  crossterms = (dft_control%xc_control%crossterms.AND.do_lsd)
  IF (ANY(gradient_functionals).OR.ANY(crossterms).OR.do_lsd) THEN
     CALL stop_program(routineP, "not yet implemented")
  END IF

  CALL xc_dset_create(derivative_set, pw_pool)
  CALL xc_calc_derivatives(derivative_set, dft_control%xc_control, order, rgrid)
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
              IF(order==3)THEN
                 v2%cr3d(i,j,k) = derivative_set%data(ii,3)*rhop_r%cr3d(i,j,k)
              END IF
           END DO
        END DO
     END DO
     DO k = bo(1,3), bo(2,3)
        DO j = bo(1,2), bo(2,2)
           DO i = bo(1,1), bo(2,1)
              IF(order==3)THEN
                 v2%cr3d(i,j,k) =v2%cr3d(i,j,k)*rhop_r%cr3d(i,j,k)
              END IF
              e1= e1+ v%cr3d(i,j,k)*rhop_r%cr3d(i,j,k)
              e2= e2+ v1%cr3d(i,j,k)*rhop_r%cr3d(i,j,k)
           END DO
        END DO
     END DO
     DEALLOCATE ( rgrid, STAT = stat )
     IF ( stat /= 0 ) CALL stop_memory ( routineP, "rgrid" )
  END IF
  CALL xc_dset_release(derivative_set)

END SUBROUTINE xc_calculate_derivatives

!******************************************************************************

END MODULE kg_kxc

!******************************************************************************
