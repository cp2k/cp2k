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
!!   SOURCE
!******************************************************************************
!------------------------------------------------------------------------------!
    MODULE kg_energy_optimize
!------------------------------------------------------------------------------!
!
  USE coefficient_types,               ONLY: coeff_type
  USE global_types,                    ONLY: global_environment_type
  USE kg_energy_utils,                 ONLY: diis,&
                                             energy_min_type,&
                                             min_info_type
  USE kg_environment_types,            ONLY: kg_environment_type,&
                                             kg_thermodynamic_type
  USE kg_force,                        ONLY: kg_force_control
  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE
!
  PUBLIC :: optimize_ao                                                       

  PRIVATE
!      
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!
 SUBROUTINE optimize_ao( kg_env, globenv, energy_min, iter )

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min
    INTEGER, intent(IN)                      :: iter

    INTEGER                                  :: i, isos, ncoef, nhist
    INTEGER, SAVE                            :: ihist = 0
    REAL(dbl)                                :: fconv, tol
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: cgrad
    REAL(dbl), ALLOCATABLE, DIMENSION(:), &
      SAVE                                   :: hess
    REAL(dbl), ALLOCATABLE, &
      DIMENSION(:, :), SAVE                  :: c_hist, e_hist
    REAL(dbl), SAVE                          :: step = 0.0_dbl
    TYPE(coeff_type), POINTER                :: coef_force, coef_pos, coef_vel
    TYPE(kg_thermodynamic_type), POINTER     :: thermo

!

  thermo => kg_env % thermo
  coef_pos => kg_env % coef_pos (1)
  coef_vel => kg_env % coef_vel (1)
  coef_force => kg_env % coef_force (1)

  ncoef = size(coef_pos % ao % cr)
  nhist = energy_min % ndiis

  IF ( .NOT. ALLOCATED(cgrad)) THEN 
     ALLOCATE ( cgrad (ncoef), STAT=isos )
     IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'cgrad', ncoef )
  END IF

  SELECT CASE ( energy_min % mintype )
  CASE ('SD')
    CALL sd_ao ( kg_env, globenv, energy_min, iter, step)     
    DO i =1,ncoef
       cgrad (i) =  coef_force % ao % cr (i)
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
    CALL diis( coef_pos %ao %cr, coef_force %ao % cr, &
               energy_min % sd_delta, &
               fconv, hess, e_hist, c_hist, nhist, ihist)
    coef_vel % ao % cr (:) = 0.0_dbl
    CALL kg_force_control ( kg_env, globenv )
    DO i =1,ncoef
       cgrad (i) = coef_force % ao % cr (i)
    END DO
    tol = maxval(abs(cgrad(:)))
  END SELECT                                                                     

  CALL copy_energy_info ( energy_min % min_info, thermo, tol, fconv )

  IF (ALLOCATED (cgrad)) DEALLOCATE (cgrad, STAT = isos)

  END SUBROUTINE optimize_ao
!------------------------------------------------------------------------------!
!
  SUBROUTINE sd_ao ( kg_env, globenv, energy_min, iter, step)

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min
    INTEGER, intent(in)                      :: iter
    REAL(dbl), intent(inout)                 :: step

    INteger                                  :: i, ios, iw, ncoef
    REAL(dbl)                                :: deltae, ezero
    REAL(dbl), ALLOCATABLE, DIMENSION(:), &
      SAVE                                   :: c_old, f_old
    TYPE(coeff_type), POINTER                :: coef_force, coef_pos, coef_vel
    TYPE(kg_thermodynamic_type), POINTER     :: thermo

!

  thermo => kg_env % thermo
  coef_pos => kg_env % coef_pos (1)
  coef_vel => kg_env % coef_vel (1)
  coef_force => kg_env % coef_force (1)

  ezero = 1.0E-50_dbl
  ncoef = size ( coef_pos % ao % cr )
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
        c_old (i) = coef_pos % ao % cr (i)
        f_old (i) = coef_force % ao % cr (i)
     END DO
  END IF
  coef_vel % ao % cr (:) = 0.0_dbl
  DO i = 1,ncoef
     coef_pos % ao % cr (i) =  coef_pos % ao % cr (i) - &
                     step * coef_force% ao % cr (i)
  END DO

! calculates new forces on coefficients
  coef_force % ao % cr (:) = 0.0_dbl
  CALL kg_force_control ( kg_env, globenv )

! change time step if needed
  deltae = thermo%pot - energy_min % min_info %e
  IF (deltae >= ezero) THEN
     DO i = 1, ncoef
        coef_pos % ao % cr (i) = c_old (i)
        coef_force % ao % cr (i) = f_old (i)
     END DO
     step = step*0.5_dbl
     WRITE (iw,'(A)') ' WARNING: RESET COEFFICIENTS, DELTAE > 0 '
     WRITE (iw,'(A,T66,E15.7)') ' NEW SD STEP      = ', STEP
  ELSE
     DO i = 1, ncoef
        c_old(i) = coef_pos % ao % cr (i)
        f_old(i) = coef_force % ao % cr (i)
     END DO
  END IF
  END SUBROUTINE sd_ao
!      
!------------------------------------------------------------------------------!
  SUBROUTINE copy_energy_info ( min_info, thermo, tol, fconv )
    TYPE(min_info_type), INTENT(INOUT)       :: min_info
    TYPE(kg_thermodynamic_type), INTENT(IN)  :: thermo
    REAL(dbl), intent(in)                    :: tol, fconv

    REAL(dbl)                                :: eold

!

  eold = min_info %e
  min_info %e = thermo % pot
  min_info %e_xc = thermo % e_xc
  min_info %e_hartree = thermo % e_hartree
  min_info %e_ppl = thermo % e_ppl
  min_info %e_ppl2 = thermo % e_ppl2
  min_info %de = thermo % pot - eold
  min_info %tol = tol
  min_info %err = fconv
 END  SUBROUTINE copy_energy_info 
!------------------------------------------------------------------------------!
    END MODULE kg_energy_optimize
!--------------------------------------------------------------------------------!
