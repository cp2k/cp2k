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
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE pw_grid_types,                   ONLY: PW_MODE_DISTRIBUTED
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE pw_types,                        ONLY: pw_type, pw_p_type,&
       pw_zero
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc, only: xc_rho_set_and_dset_create
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative,&
                                             xc_dset_release
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_rho_set_types, ONLY: xc_rho_set_type, xc_rho_set_release
  USE cp_subsystem_types, only: cp_subsystem_p_type
  USE simulation_cell, only:cell_type
  USE cp_linked_list_xc_deriv,         ONLY: cp_sll_xc_deriv_next,&
                                             cp_sll_xc_deriv_type
  USE xc_derivative_desc,              ONLY: MAX_DERIVATIVE_DESC_LENGTH,&
                                             MAX_LABEL_LENGTH

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
    REAL(KIND=dp), INTENT(out)                   :: e1_xc, e2_xc
    TYPE(dft_control_type), POINTER          :: dft_control
    LOGICAL, INTENT(in)                      :: calculate_forces

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "calc_kxc_derivatives (MODULE kg_kxc)"

    INTEGER                                  :: handle, order
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER   :: v1xc, v2xc

!--------------------------------------------------------

  CALL timeset (routine,'I',' ',handle )

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

  e2_xc =  0.5_dp*rho_r%pw%pw_grid%dvol*e2_xc      

  v1xc (:,:,:) =  v1xc(:,:,:) *rho_r%pw%pw_grid%dvol
  IF (calculate_forces)THEN
    v2xc (:,:,:) =  0.5_dp * v2xc(:,:,:) *rho_r%pw%pw_grid%dvol
  END IF
  CALL timestop(0.0_dp,handle)
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


   TYPE(pw_type), INTENT(IN), target        :: rho_r, rhop_r, v
   TYPE(pw_type), INTENT(INOUT)             :: v1, v2
   REAL(KIND=dp), INTENT(out)                   :: e1, e2
   INTEGER, intent(IN)                      :: order
   TYPE(dft_control_type), POINTER          :: dft_control
   TYPE(pw_pool_type), POINTER              :: pw_pool
   TYPE(cell_type), pointer                 :: cell

   CHARACTER(len=*), PARAMETER :: routineN = 'xc_calculate_derivatives', &
        routineP = moduleN//':'//routineN

   INTEGER                                  :: i, ii, j, k, n(3), ntot, stat,&
        deriv_order
   INTEGER, DIMENSION(2, 3)                 :: bo
   INTEGER, DIMENSION(3, 3)                 :: nd
   LOGICAL                                  :: do_lsd, failure, gradient_f
   LOGICAL, DIMENSION(3)                    :: crossterms, &
        gradient_functionals
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)     :: rgrid
   TYPE(xc_derivative_set_type), POINTER    :: derivative_set
   TYPE(xc_rho_set_type), pointer           :: rho_set
   TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_rho_g,my_tau
   TYPE(pw_p_type), DIMENSION(1), target    :: my_rho_r
   TYPE(cp_sll_xc_deriv_type), POINTER      :: pos
    CHARACTER&
      (len=MAX_DERIVATIVE_DESC_LENGTH)       :: desc
    TYPE(xc_derivative_type), POINTER        :: deriv
    REAL(KIND=dp), DIMENSION(:,:,:), pointer :: pot
   !------------------------------------------------------------------------------

   failure=.FALSE.

   NULLIFY(derivative_set,rho_set)
   nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))
   e1 = 0.0_dp
   e2 = 0.0_dp

   ! no gradient corrections are implemented for derivatives
   gradient_functionals = dft_control%xc_control%gradient_functionals
   do_lsd = (dft_control%nspins > 1)

   NULLIFY(rho_set,derivative_set,my_rho_g,my_tau)
   my_rho_r(1)%pw => rho_r
   bo=rho_r%pw_grid%bounds_local
   CALL xc_rho_set_and_dset_create(rho_set=rho_set,&
        deriv_set=derivative_set,&
        deriv_order=order,&
        rho_r=my_rho_r,rho_g=my_rho_g,tau=my_tau,dft_control=dft_control,&
        xc_control=dft_control%xc_control,cell=cell,pw_pool=pw_pool,&
        gradient_f=gradient_f,&
        needs_basic_components=.FALSE.)

   ! check for unsupported derivatives
   pos => derivative_set%derivs
   DO WHILE (cp_sll_xc_deriv_next(pos,el_att=deriv))
      CALL xc_derivative_get(deriv,order=deriv_order,&
           desc=desc)
      IF (deriv_order==2.AND.desc/="(rho)(rho)".OR.&
           order==3.and.deriv_order==3.AND.desc/="(rho)(rho)(rho)") THEN
         call stop_program(routine=routineP,&
              error_message="unsupported xc derivative: "//&
              trim(desc))
      END IF
   END DO

   call xc_rho_set_release(rho_set)
   deriv => xc_dset_get_derivative(derivative_set,"(rho)(rho)")
   IF (ASSOCIATED(deriv)) THEN
      CALL xc_derivative_get(deriv,deriv_data=pot)
      DO k = bo(1,3), bo(2,3)
         DO j = bo(1,2), bo(2,2)
            DO i = bo(1,1), bo(2,1)
               v1%cr3d(i,j,k) = pot(i,j,k)*rhop_r%cr3d(i,j,k)
            END DO
         END DO
      END DO
   ELSE
      CALL pw_zero(v1)
   END IF
   IF(order==3)THEN
      deriv => xc_dset_get_derivative(derivative_set,"(rho)(rho)(rho)")
      IF (ASSOCIATED(deriv)) THEN
         CALL xc_derivative_get(deriv,deriv_data=pot)
         DO k = bo(1,3), bo(2,3)
            DO j = bo(1,2), bo(2,2)
               DO i = bo(1,1), bo(2,1)
                  v2%cr3d(i,j,k) = pot(i,j,k)*rhop_r%cr3d(i,j,k)
               END DO
            END DO
         END DO
      ELSE
         call pw_zero(v2)
      END IF
   END IF
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
   CALL xc_dset_release(derivative_set)

 END SUBROUTINE xc_calculate_derivatives

!******************************************************************************

END MODULE kg_kxc

!******************************************************************************
