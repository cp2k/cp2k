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

  USE coefficient_types, ONLY: coeff_type
  USE functionals
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi, zero
  USE message_passing, ONLY : mp_sum, mp_max, mp_sync,&
                  mp_range, mp_allgather
  USE termination, ONLY : stop_program, stop_memory
  USE timings, ONLY : timeset, timestop
  USE xc_thomas_fermi, only: thomas_fermi

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calculate_tf_potential,  calculate_xc_potential

!******************************************************************************

CONTAINS

!******************************************************************************
!!****** kgkxc/calculate_kxc_potential [1.1] *
!!
!!   NAME
!!     calculate_kxc_potential
!!
!!   FUNCTION
!!     Calculates the total potential energy, total force, and the
!!     total pressure tensor from the potentials
!!
!!   AUTHOR
!!     gt         
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE calculate_tf_potential ( rho_elec, v_r, e_tf, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( coeff_type ), INTENT (in) :: rho_elec
  TYPE ( coeff_type ), INTENT (out) :: v_r
  REAL ( dbl ), INTENT (out) :: e_tf
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Local parameters
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_tf_potential (MODULE kg_kxc)"
! Locals
  INTEGER :: n1, n2, n3, npoints, i, j, k, ip, handle, isos
  INTEGER :: order
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: pot
  REAL (dbl), DIMENSION(:,:,:), POINTER :: rhoe
  REAL (dbl), DIMENSION(:,:,:), POINTER :: vxc
  REAL(DBL), DIMENSION(:), ALLOCATABLE :: rho
!------------------------------------------------------------------------------

  CALL timeset ( 'calculate_tfxc_potential','I',' ',handle )

  IF(0==0) THEN
    rhoe =>rho_elec%pw%cr3d(:,:,:)
    vxc =>v_r%pw%cr3d(:,:,:)
    n1 =  SIZE(rhoe,1)
    n2 =  SIZE(rhoe,2)
    n3 =  SIZE(rhoe,3)

! allocates auxiliary pot array                    
    npoints = n1*n2*n3  
    ALLOCATE (pot(npoints,2),STAT=isos)
    IF (isos /= 0) CALL stop_memory(routine,"pot",npoints*2)
    ALLOCATE (rho(npoints),STAT=isos)
    IF (isos /= 0) CALL stop_memory(routine,"rho",npoints)
    pot(:,:)=0.0_dbl
    e_tf = 0.0_dbl

    ip = 0
    DO i = 1, n1
      DO j = 1, n2
        DO k = 1, n3
          ip = ip + 1
          rho(ip) =  -rhoe(i,j,k)
        END DO
      END DO
    END DO

! calculate tf energy 
    order = 1
    CALL thomas_fermi ( rho, pot, order)                    
    DO ip = 1, npoints
      e_tf = e_tf + pot(ip,1)
    END DO
    e_tf = rho_elec%pw%pw_grid%dvol * e_tf
  
! get pot in the correct form 
    ip = 0
    DO i = 1, n1
      DO j = 1, n2
        DO k = 1, n3
          ip = ip + 1
          vxc(i,j,k) = -pot (ip, 2) * rho_elec%pw%pw_grid%dvol
        END DO
      END DO
    END DO

    DEALLOCATE (pot,STAT=isos)
    DEALLOCATE (rho,STAT=isos)

  ELSE
!
!   call xc_calculate ()
!
  ENDIF


 END SUBROUTINE calculate_tf_potential

!******************************************************************************

 SUBROUTINE calculate_xc_potential ( rho_elec, v_r, e_xc, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( coeff_type ), INTENT (in) :: rho_elec
  TYPE ( coeff_type ), INTENT (out) :: v_r
  REAL ( dbl ), INTENT (out) :: e_xc
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Local parameters
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_xc_potential (MODULE kg_kxc)"
! Locals
  INTEGER :: ORDER
  REAL (dbl), DIMENSION(:,:,:), POINTER :: rho
  REAL (dbl), DIMENSION(:,:,:), POINTER :: vxc
!--------------------------------------------------------

  IF(0==0) THEN
    e_xc = 0.0_dbl
    rho => rho_elec%pw%cr3d(:,:,:)
    vxc => v_r%pw%cr3d(:,:,:)
    rho (:,:,:) = -1.0_dbl * rho(:,:,:)
    call pade(rho(:,:,:),e_xc,vxc(:,:,:))
    e_xc = rho_elec%pw%pw_grid%dvol*e_xc   
    rho (:,:,:) = -1.0_dbl * rho(:,:,:)
    vxc (:,:,:) = -1.0_dbl * vxc(:,:,:)* rho_elec%pw%pw_grid%dvol
  ELSE
!
!   call xc_calculate ()
!
  ENDIF

!!!call timestop(0.0_dbl,handle)

 END SUBROUTINE calculate_xc_potential

!******************************************************************************

END MODULE kg_kxc  

!******************************************************************************
