!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!
!!****** cp2k/energy_optimize [1.0] *
!!
!!   NAME
!!     energy_optimize
!!
!!   FUNCTION
!!   energy optimization routine
!!
!!   AUTHOR
!!     gloria   
!!
!!   SOURCE
!******************************************************************************
!------------------------------------------------------------------------------!
    MODULE energy_optimize
!------------------------------------------------------------------------------!
!
  USE brillouin, ONLY : kpoint_type
  USE dft_types, ONLY : dft_control_type
  USE energy_input, ONLY : energy_min_type, min_info_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE force_control, ONLY : force
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type, thermodynamic_type              
  USE pol_diis, ONLY : diis
  USE structure_types, ONLY : structure_type, interaction_type
  USE termination, ONLY : stop_program, stop_memory
  IMPLICIT NONE
!
  PUBLIC :: optimize                                          
  PRIVATE
!      
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!
 SUBROUTINE optimize ( struc, inter, thermo, simpar, ewald_param, box_change, &
                        globenv, dft_control, kp, energy_min, iter)        

  IMPLICIT NONE
! Arguments
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  LOGICAL, INTENT ( IN ) :: box_change
  TYPE ( interaction_type ), INTENT ( IN ) :: inter
  TYPE ( global_environment_type), INTENT ( IN ) :: globenv
  TYPE ( dft_control_type ), INTENT ( IN ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  TYPE ( energy_min_type ), INTENT (INOUT) :: energy_min
  INTEGER, intent(IN) :: iter      

! locals
  INTEGER :: i, itimes, ncoef, isos, nhist
  INTEGER, SAVE :: ihist
  REAL(dbl) :: eold,cons,fconv,tol, deltae
  REAL (dbl), SAVE :: step
  REAL(dbl), DIMENSION(:), ALLOCATABLE :: cgrad              
  REAL(dbl), DIMENSION(:), ALLOCATABLE, SAVE ::  hess
  REAL(dbl), DIMENSION(:,:), ALLOCATABLE, SAVE :: c_hist,e_hist
!
  ncoef = size(struc % coef_pos (1) % ao % cr)
  nhist = energy_min % ndiis

  IF ( .NOT. ALLOCATED(cgrad)) THEN 
     ALLOCATE ( cgrad (ncoef), STAT=isos )
     IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'cgrad', ncoef )
  END IF

  SELECT CASE ( energy_min % mintype )
  CASE ('SD')
    CALL sd ( struc, inter, thermo, simpar, ewald_param, box_change, &
                      globenv, dft_control, kp, energy_min, iter, step)
    DO i =1,ncoef
       cgrad (i) = struc % coef_force (1) % ao % cr (i)
    END DO
    tol = maxval(abs(cgrad(:)))

! diis optimization
  CASE ('DIIS')

    IF(iter==1)THEN
      IF( .NOT. ALLOCATED (hess)) THEN
         ALLOCATE ( hess (ncoef), STAT = isos)
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'hess', ncoef )
      END IF
      IF( .NOT. ALLOCATED(c_hist)) THEN
         ALLOCATE ( c_hist (nhist,ncoef), STAT = isos )  
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'c_hist', ncoef*nhist)
      END IF
      IF( .NOT. ALLOCATED (e_hist)) THEN
          ALLOCATE ( e_hist (nhist,ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'e_hist', ncoef*nhist )
      END IF
      ihist=0
      DO i = 1,ncoef 
        hess(i)=1.0_dbl                                               
      END DO  
    ENDIF
    CALL diis( struc % coef_pos (1) %ao %cr, struc % coef_force (1) %ao % cr, &
               energy_min % sd_delta, &
               fconv, hess, e_hist, c_hist, nhist, ihist)
    CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
       globenv, dft_control, kp )
    DO i =1,ncoef
       cgrad (i) = struc % coef_force (1) % ao % cr (i)
    END DO
    tol = maxval(abs(cgrad(:)))
  END SELECT                                                                     

  eold = energy_min % min_info %e
  energy_min % min_info %e = thermo % pot
  energy_min % min_info %e_class = thermo % e_class      
  energy_min % min_info %e_over = thermo % e_overlap
  energy_min % min_info %e_electro = thermo % e_electro
  energy_min % min_info %e_thk = thermo % e_thk
  energy_min % min_info %de = thermo % pot - eold
  energy_min % min_info %tol = tol                     
  energy_min % min_info %err = fconv                     

  END SUBROUTINE optimize
!------------------------------------------------------------------------------!
  SUBROUTINE sd ( struc, inter, thermo, simpar, ewald_param, box_change, &
                        globenv, dft_control, kp, energy_min, iter, step)          

  IMPLICIT NONE
! Arguments
  TYPE ( structure_type ), INTENT ( INOUT ) :: struc
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  LOGICAL, INTENT ( IN ) :: box_change
  TYPE ( interaction_type ), INTENT ( IN ) :: inter
  TYPE ( global_environment_type), INTENT ( IN ) :: globenv
  TYPE ( dft_control_type ), INTENT ( IN ) :: dft_control
  TYPE ( kpoint_type ), INTENT ( IN ) :: kp
  TYPE ( energy_min_type ), INTENT (INOUT) :: energy_min
  INTEGER, intent (in) :: iter
  REAL ( dbl ), intent (inout) :: step

! locals
  REAL( dbl ), DIMENSION(:), ALLOCATABLE, SAVE ::  c_old, f_old
  INteger :: ios, ncoef, iw, i
  REAL ( dbl ) :: ezero, deltae

  ezero = 1.0E-50_dbl
  ncoef = size ( struc % coef_pos (1) % ao % cr )
  iw = globenv % scr

  IF ( .NOT. allocated(c_old)) THEN 
     ALLOCATE ( c_old (ncoef), STAT=ios)
     IF ( ios /= 0 )  CALL stop_memory ( 'sd', 'c_old', ncoef )
  END IF
  IF ( .NOT. allocated(f_old)) THEN
     ALLOCATE (f_old(ncoef),STAT=ios)
     IF ( ios /= 0 )  CALL stop_memory ( 'sd', 'f_old', ncoef)
  END IF

!set coefficient velocities to 0 and move coefficients
  IF (iter==1) THEN
     step = energy_min % sd_delta
     DO i = 1,ncoef
        c_old (i) = struc % coef_pos (1) % ao % cr (i)
        f_old (i) = struc % coef_force (1) % ao % cr (i)
     END DO
  END IF
  struc % coef_vel (1) % ao % cr (:) = 0.0_dbl
  DO i = 1,ncoef
     struc % coef_pos (1) % ao % cr (i) =  struc % coef_pos (1) % ao % cr (i) - &
                     step * struc % coef_force (1) % ao % cr (i)
  END DO

! calculates new forces on coefficients
  struc % coef_force (1) % ao % cr (:) = 0.0_dbl
  CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
              globenv, dft_control, kp )

! change time step if needed
  deltae = thermo%pot - energy_min % min_info %e
  IF (deltae >= ezero) THEN
     DO i = 1, ncoef
        struc % coef_pos (1) % ao % cr (i) = c_old (i)
        struc % coef_force (1) % ao % cr (i) = f_old (i)
     END DO
     step = step*0.5_dbl
     WRITE (iw,'(A)') ' WARNING: RESET COEFFICIENTS, DELTAE > 0 '
     WRITE (iw,'(A,T66,E15.7)') ' NEW SD STEP      = ', STEP
  ELSE
     DO i = 1, ncoef
        c_old(i) = struc % coef_pos (1) % ao % cr (i)
        f_old(i) = struc % coef_force (1) % ao % cr (i)
     END DO
  END IF
  END SUBROUTINE sd
!      
!------------------------------------------------------------------------------!
    END MODULE energy_optimize
!--------------------------------------------------------------------------------!
