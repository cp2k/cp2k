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
  USE brillouin,                       ONLY: kpoint_type
  USE cp_control_types,                ONLY: dft_control_type
  USE energy_input,                    ONLY: energy_min_type,&
                                             min_info_type
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE fit_force_control,               ONLY: fit_force
  USE fit_parameters_types,            ONLY: fit_parameter_type
  USE force_env_types,                 ONLY: force
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE md,                              ONLY: simulation_parameters_type,&
                                             thermodynamic_type
  USE pol_diis,                        ONLY: diis
  USE structure_types,                 ONLY: interaction_type,&
                                             structure_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE
!
  PUBLIC :: optimize_ao, optimize_pw                                          
  INTERFACE optimize_ao
    MODULE PROCEDURE optimize_ao_pol, optimize_ao_fit
  END INTERFACE
  INTERFACE optimize_pw
    MODULE PROCEDURE optimize_pw_pol, optimize_pw_fit
  END INTERFACE

  PRIVATE
!      
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!
 SUBROUTINE optimize_ao_pol ( struc, inter, thermo, simpar, ewald_param, box_change, &
                        globenv, dft_control, kp, energy_min, iter)        

    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(interaction_type), INTENT(IN)       :: inter
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    LOGICAL, INTENT(IN)                      :: box_change
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(dft_control_type), INTENT(IN)       :: dft_control
    TYPE(kpoint_type), INTENT(IN)            :: kp
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min
    INTEGER, intent(IN)                      :: iter

    INTEGER                                  :: i, isos, ncoef, nhist
    INTEGER, SAVE                            :: ihist
    REAL(dbl)                                :: fconv, tol
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: cgrad
    REAL(dbl), ALLOCATABLE, DIMENSION(:), &
      SAVE                                   :: hess
    REAL(dbl), ALLOCATABLE, &
      DIMENSION(:, :), SAVE                  :: c_hist, e_hist
    REAL(dbl), SAVE                          :: step

!

  ncoef = size(struc % coef_pos (1) % ao % cr)
  nhist = energy_min % ndiis

  IF ( .NOT. ALLOCATED(cgrad)) THEN 
     ALLOCATE ( cgrad (ncoef), STAT=isos )
     IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'cgrad', ncoef )
  END IF

  SELECT CASE ( energy_min % mintype )
  CASE ('SD')
    CALL sd_ao ( struc, inter, thermo, simpar, ewald_param, box_change, &
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

  CALL copy_energy_info ( energy_min % min_info, thermo, tol, fconv )

  IF (ALLOCATED (cgrad)) DEALLOCATE (cgrad, STAT = isos)

  END SUBROUTINE optimize_ao_pol
!------------------------------------------------------------------------------!
!
 SUBROUTINE optimize_pw_pol ( struc, inter, thermo, simpar, ewald_param, box_change, &
                        globenv, dft_control, kp, energy_min, iter)        

    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(interaction_type), INTENT(IN)       :: inter
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    LOGICAL, INTENT(IN)                      :: box_change
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(dft_control_type), INTENT(IN)       :: dft_control
    TYPE(kpoint_type), INTENT(IN)            :: kp
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min
    INTEGER, intent(IN)                      :: iter

    INTEGER                                  :: i, isos, ncoef, nhist
    INTEGER, SAVE                            :: ihist, jhist
    REAL(dbl)                                :: fconv, ivol, tol, tol1, tol2
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: agrad, bgrad
    REAL(dbl), ALLOCATABLE, DIMENSION(:), &
      SAVE                                   :: ci, cr, fi, fr, hess_i, hess_r
    REAL(dbl), ALLOCATABLE, &
      DIMENSION(:, :), SAVE                  :: ci_hist, cr_hist, ei_hist, &
                                                er_hist
    REAL(dbl), SAVE                          :: step

!

  ncoef = size(struc % coef_pos (1) % pw % cc)
  nhist = energy_min % ndiis
  ivol = 1._dbl / struc % box % deth

  IF ( .NOT. ALLOCATED(agrad)) THEN 
     ALLOCATE ( agrad (ncoef), STAT=isos )
     IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'agrad', ncoef )
  END IF

  IF ( .NOT. ALLOCATED(bgrad)) THEN 
     ALLOCATE ( bgrad (ncoef), STAT=isos )
     IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'bgrad', ncoef )
  END IF

  SELECT CASE ( energy_min % mintype )
  CASE ('SD')
    CALL sd_pw ( struc, inter, thermo, simpar, ewald_param, box_change, &
                      globenv, dft_control, kp, energy_min, iter, step)
    DO i =1,ncoef
       agrad (i) =  REAL (struc % coef_force (1) % pw % cc (i), dbl)
       bgrad (i) =  AIMAG (struc % coef_force (1) % pw % cc (i)) 
    END DO
    tol1 = maxval(abs(agrad(:)))
    tol2 = maxval(abs(bgrad(:)))
    tol = max(tol1,tol2)

! diis optimization
  CASE ('DIIS')

    IF(iter==1)THEN
      IF( .NOT. ALLOCATED (hess_r)) THEN
         ALLOCATE ( hess_r (ncoef), STAT = isos)
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'hess', ncoef )
      END IF
      IF( .NOT. ALLOCATED(cr_hist)) THEN
         ALLOCATE ( cr_hist (nhist,ncoef), STAT = isos )  
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'cr_hist', ncoef*nhist)
      END IF
      IF( .NOT. ALLOCATED (er_hist)) THEN
          ALLOCATE ( er_hist (nhist,ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'er_hist', ncoef*nhist )
      END IF
      IF( .NOT. ALLOCATED (hess_i)) THEN
         ALLOCATE ( hess_i (ncoef), STAT = isos)
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'hess', ncoef )
      END IF
      IF( .NOT. ALLOCATED(ci_hist)) THEN
         ALLOCATE ( ci_hist (nhist,ncoef), STAT = isos )  
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'ci_hist', ncoef*nhist)
      END IF
      IF( .NOT. ALLOCATED (ei_hist)) THEN
          ALLOCATE ( ei_hist (nhist,ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'ei_hist', ncoef*nhist )
      END IF
      IF( .NOT. ALLOCATED (fr)) THEN
          ALLOCATE ( fr (ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'fr', ncoef )
      END IF
      IF( .NOT. ALLOCATED (fi)) THEN
          ALLOCATE ( fi (ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'fi', ncoef )
      END IF
      IF( .NOT. ALLOCATED (cr)) THEN
          ALLOCATE ( cr (ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'cr', ncoef )
      END IF
      IF( .NOT. ALLOCATED (ci)) THEN
          ALLOCATE ( ci (ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'ci', ncoef )
      END IF
      ihist=0
      jhist=0
      DO i = 2,ncoef 
        hess_r(i)=1.0_dbl/struc % coef_pos (1) % pw % pw_grid % gsq ( i )   
        hess_i(i)=1.0_dbl /struc % coef_pos (1) % pw % pw_grid % gsq ( i )   
      END DO  
      hess_r(1) = 1.0_dbl
      hess_i(1) = 1.0_dbl
    ENDIF
    DO i = 1, ncoef
      cr (i) = REAL ( struc % coef_pos (1) %pw %cc (i))
      ci (i) = AIMAG ( struc % coef_pos (1) %pw %cc(i))
      fr (i) = REAL ( struc % coef_force (1) %pw %cc (i))
      fi (i) = AIMAG ( struc % coef_force (1) %pw %cc(i))
    END DO
    CALL diis( cr, fr, energy_min % sd_delta, &
              fconv, hess_r, er_hist, cr_hist, nhist, ihist)
    CALL diis( ci, fi, energy_min % sd_delta, &
              fconv, hess_i, ei_hist, ci_hist, nhist, jhist)
    DO I = 1, ncoef
      struc % coef_pos (1) %pw %cc (i) =  CMPLX ( (cr ( i )), (ci (i)), dbl )
      struc % coef_force (1) %pw %cc (i) = CMPLX ( (fr ( i )), (fi (i)), dbl )
    END DO
    CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
       globenv, dft_control, kp )
    DO i =1,ncoef
       agrad (i) =  REAL (struc % coef_force (1) % pw % cc (i), dbl)
       bgrad (i) =  AIMAG (struc % coef_force (1) % pw % cc (i)) 
    END DO
    tol1 = maxval(abs(agrad(:)))
    tol2 = maxval(abs(bgrad(:)))
    tol = max(tol1,tol2)
  END SELECT                                                                     

  CALL copy_energy_info ( energy_min % min_info, thermo, tol, fconv )

  IF (ALLOCATED (agrad)) DEALLOCATE (agrad, STAT = isos)
  IF (ALLOCATED (bgrad)) DEALLOCATE (bgrad, STAT = isos)

  END SUBROUTINE optimize_pw_pol
!------------------------------------------------------------------------------!
!
 SUBROUTINE optimize_ao_fit ( struc, inter, thermo, simpar, ewald_param, box_change, &
                        globenv,fit_param, energy_min, iter)

    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(interaction_type), INTENT(INOUT)    :: inter
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    LOGICAL, INTENT(IN)                      :: box_change
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(fit_parameter_type), INTENT(INOUT)  :: fit_param
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min
    INTEGER, intent(IN)                      :: iter

    INTEGER                                  :: i, isos, ncoef, nhist
    INTEGER, SAVE                            :: ihist
    REAL(dbl)                                :: fconv, tol
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: cgrad
    REAL(dbl), ALLOCATABLE, DIMENSION(:), &
      SAVE                                   :: hess
    REAL(dbl), ALLOCATABLE, &
      DIMENSION(:, :), SAVE                  :: c_hist, e_hist

!

  ncoef = size(struc % coef_pos (1) % ao % cr)
  nhist = energy_min % ndiis

  IF ( .NOT. ALLOCATED(cgrad)) THEN
     ALLOCATE ( cgrad (ncoef), STAT=isos )
     IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'cgrad', ncoef )
  END IF

  SELECT CASE ( energy_min % mintype )
!  CASE ('SD')
!    CALL sd ( struc, inter, thermo, simpar, ewald_param, box_change, &
!                      globenv, dft_control, kp, energy_min, iter, step)
!    DO i =1,ncoef
!       cgrad (i) = struc % coef_force (1) % ao % cr (i)
!    END DO
!     tol = maxval(abs(cgrad(:)))

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
    CALL fit_force ( struc, inter, thermo, simpar, ewald_param, box_change, fit_param, &
       globenv )
    DO i =1,ncoef
       cgrad (i) = struc % coef_force (1) % ao % cr (i)
    END DO
    tol = maxval(abs(cgrad(:)))
  END SELECT

  CALL copy_energy_info ( energy_min % min_info, thermo, tol, fconv ) 

  IF (ALLOCATED (cgrad)) DEALLOCATE (cgrad, STAT = isos)

  END SUBROUTINE optimize_ao_fit
!
!------------------------------------------------------------------------------!
!
 SUBROUTINE optimize_pw_fit ( struc, inter, thermo, simpar, ewald_param, box_change, &
                        globenv,fit_param, energy_min, iter)

    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(interaction_type), INTENT(INOUT)    :: inter
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    LOGICAL, INTENT(IN)                      :: box_change
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(fit_parameter_type), INTENT(INOUT)  :: fit_param
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min
    INTEGER, intent(IN)                      :: iter

    INTEGER                                  :: i, isos, ncoef, nhist
    INTEGER, SAVE                            :: ihist, jhist
    REAL(dbl)                                :: fconv, ivol, tol, tol1, tol2
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: agrad, bgrad
    REAL(dbl), ALLOCATABLE, DIMENSION(:), &
      SAVE                                   :: ci, cr, fi, fr, hess_i, hess_r
    REAL(dbl), ALLOCATABLE, &
      DIMENSION(:, :), SAVE                  :: ci_hist, cr_hist, ei_hist, &
                                                er_hist

!

  ncoef = size(struc % coef_pos (1) % pw % cc)
  nhist = energy_min % ndiis
  ivol = 1._dbl / struc % box % deth

  IF ( .NOT. ALLOCATED(agrad)) THEN
     ALLOCATE ( agrad (ncoef), STAT=isos )
     IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'agrad', ncoef )
  END IF

  IF ( .NOT. ALLOCATED(bgrad)) THEN 
     ALLOCATE ( bgrad (ncoef), STAT=isos )
     IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'bgrad', ncoef )
  END IF

  SELECT CASE ( energy_min % mintype )
!   CASE ('SD')
!   CALL sd_pw ( struc, inter, thermo, simpar, ewald_param, box_change, &
!                     globenv, dft_control, kp, energy_min, iter, step)
!   DO i =1,ncoef
!      agrad (i) =  REAL (struc % coef_force (1) % pw % cc (i), dbl)
!      bgrad (i) =  AIMAG (struc % coef_force (1) % pw % cc (i)) 
!   END DO
!   tol1 = maxval(abs(agrad(:)))
!   tol2 = maxval(abs(bgrad(:)))
!   tol = max(tol1,tol2)

! diis optimization
  CASE ('DIIS')

    IF(iter==1)THEN
      IF( .NOT. ALLOCATED (hess_r)) THEN
         ALLOCATE ( hess_r (ncoef), STAT = isos)
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'hess', ncoef )
      END IF
      IF( .NOT. ALLOCATED(cr_hist)) THEN
         ALLOCATE ( cr_hist (nhist,ncoef), STAT = isos )  
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'cr_hist', ncoef*nhist)
      END IF
      IF( .NOT. ALLOCATED (er_hist)) THEN
          ALLOCATE ( er_hist (nhist,ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'er_hist', ncoef*nhist )
      END IF
      IF( .NOT. ALLOCATED (hess_i)) THEN
         ALLOCATE ( hess_i (ncoef), STAT = isos)
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'hess', ncoef )
      END IF
      IF( .NOT. ALLOCATED(ci_hist)) THEN
         ALLOCATE ( ci_hist (nhist,ncoef), STAT = isos )  
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'ci_hist', ncoef*nhist)
      END IF
      IF( .NOT. ALLOCATED (ei_hist)) THEN
          ALLOCATE ( ei_hist (nhist,ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'ei_hist', ncoef*nhist )
      END IF
      IF( .NOT. ALLOCATED (fr)) THEN
          ALLOCATE ( fr (ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'fr', ncoef )
      END IF
      IF( .NOT. ALLOCATED (fi)) THEN
          ALLOCATE ( fi (ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'fi', ncoef )
      END IF
      IF( .NOT. ALLOCATED (cr)) THEN
          ALLOCATE ( cr (ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'cr', ncoef )
      END IF
      IF( .NOT. ALLOCATED (ci)) THEN
          ALLOCATE ( ci (ncoef), STAT = isos )
         IF ( isos /= 0 )  CALL stop_memory ( 'minimizer', 'ci', ncoef )
      END IF
      ihist=0
      jhist=0
      DO i = 2,ncoef 
        hess_r(i)=1.0_dbl/struc % coef_pos (1) % pw % pw_grid % gsq ( i )   
        hess_i(i)=1.0_dbl /struc % coef_pos (1) % pw % pw_grid % gsq ( i )   
      END DO  
      hess_r(1) = 1.0_dbl
      hess_i(1) = 1.0_dbl
    ENDIF
    DO i = 1, ncoef
      cr (i) = REAL ( struc % coef_pos (1) %pw %cc (i))
      ci (i) = AIMAG ( struc % coef_pos (1) %pw %cc(i))
      fr (i) = REAL ( struc % coef_force (1) %pw %cc (i))
      fi (i) = AIMAG ( struc % coef_force (1) %pw %cc(i))
    END DO
    CALL diis( cr, fr, energy_min % sd_delta, &
              fconv, hess_r, er_hist, cr_hist, nhist, ihist)
    CALL diis( ci, fi, energy_min % sd_delta, &
              fconv, hess_i, ei_hist, ci_hist, nhist, jhist)
    DO I = 1, ncoef
      struc % coef_pos (1) %pw %cc (i) =  CMPLX ( (cr ( i )), (ci (i)), dbl )
      struc % coef_force (1) %pw %cc (i) = CMPLX ( (fr ( i )), (fi (i)), dbl )
    END DO
   CALL fit_force ( struc, inter, thermo, simpar, ewald_param, box_change, fit_param, &
       globenv )
    DO i =1,ncoef
       agrad (i) =  REAL (struc % coef_force (1) % pw % cc (i), dbl)
       bgrad (i) =  AIMAG (struc % coef_force (1) % pw % cc (i)) 
    END DO
    tol1 = maxval(abs(agrad(:)))
    tol2 = maxval(abs(bgrad(:)))
    tol = max(tol1,tol2)
  END SELECT                                                                     

  CALL copy_energy_info ( energy_min % min_info, thermo, tol, fconv ) 

  IF (ALLOCATED (agrad)) DEALLOCATE (agrad, STAT = isos)
  IF (ALLOCATED (bgrad)) DEALLOCATE (bgrad, STAT = isos)

  END SUBROUTINE optimize_pw_fit
!
!------------------------------------------------------------------------------!
  SUBROUTINE sd_ao ( struc, inter, thermo, simpar, ewald_param, box_change, &
                        globenv, dft_control, kp, energy_min, iter, step)          

    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(interaction_type), INTENT(IN)       :: inter
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    LOGICAL, INTENT(IN)                      :: box_change
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(dft_control_type), INTENT(IN)       :: dft_control
    TYPE(kpoint_type), INTENT(IN)            :: kp
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min
    INTEGER, intent(in)                      :: iter
    REAL(dbl), intent(inout)                 :: step

    INteger                                  :: i, ios, iw, ncoef
    REAL(dbl)                                :: deltae, ezero
    REAL(dbl), ALLOCATABLE, DIMENSION(:), &
      SAVE                                   :: c_old, f_old

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
  END SUBROUTINE sd_ao
!      
!------------------------------------------------------------------------------!
!      
  SUBROUTINE sd_pw ( struc, inter, thermo, simpar, ewald_param, box_change, &
                        globenv, dft_control, kp, energy_min, iter, step)          

    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(interaction_type), INTENT(IN)       :: inter
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    LOGICAL, INTENT(IN)                      :: box_change
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(dft_control_type), INTENT(IN)       :: dft_control
    TYPE(kpoint_type), INTENT(IN)            :: kp
    TYPE(energy_min_type), INTENT(INOUT)     :: energy_min
    INTEGER, intent(in)                      :: iter
    REAL(dbl), intent(inout)                 :: step

    COMPLEX(dbl), ALLOCATABLE, &
      DIMENSION(:), SAVE                     :: c_old, f_old
    INteger                                  :: i, ios, iw, ncoef
    REAL(dbl)                                :: deltae, ezero, ivol

  ezero = 1.0E-50_dbl
  ivol = 1._dbl / struc % box % deth
  ncoef = size ( struc % coef_pos (1) % pw % cc )
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
        c_old (i) = struc % coef_pos (1) % pw % cc (i)
        f_old (i) = struc % coef_force (1) % pw % cc (i)
     END DO
  END IF
  struc % coef_vel (1) % pw % cc (:) = CMPLX ( 0.0_dbl, 0.0_dbl, dbl )
  DO i = 1,ncoef
     struc % coef_pos (1) % pw % cc (i) =  (struc % coef_pos (1) % pw % cc (i) - &
                     step * struc % coef_force (1) % pw % cc (i))
  END DO

! calculates new forces on coefficients
  struc % coef_force (1) % pw % cc (:) = CMPLX ( 0.0_dbl, 0.0_dbl, dbl )
  CALL force ( struc, inter, thermo, simpar, ewald_param, box_change, &
              globenv, dft_control, kp )

! change time step if needed
  deltae = thermo%pot - energy_min % min_info %e
  IF (deltae >= ezero) THEN
     DO i = 1, ncoef
        struc % coef_pos (1) % pw % cc (i) = c_old (i)
        struc % coef_force (1) % pw % cc (i) = f_old (i)
     END DO
     step = step*0.5_dbl
     WRITE (iw,'(A)') ' WARNING: RESET COEFFICIENTS, DELTAE > 0 '
     WRITE (iw,'(A,T66,E15.7)') ' NEW SD STEP      = ', STEP
  ELSE
     DO i = 1, ncoef
        c_old(i) = struc % coef_pos (1) % pw % cc (i)
        f_old(i) = struc % coef_force (1) % pw % cc (i)
     END DO
  END IF
  END SUBROUTINE sd_pw
!------------------------------------------------------------------------------!
  SUBROUTINE copy_energy_info ( min_info, thermo, tol, fconv )
    TYPE(min_info_type), INTENT(INOUT)       :: min_info
    TYPE(thermodynamic_type), INTENT(IN)     :: thermo
    REAL(dbl), intent(in)                    :: tol, fconv

    REAL(dbl)                                :: eold

!

  eold = min_info %e
  min_info %e = thermo % pot
  min_info %e_class = thermo % e_class
  min_info %e_over = thermo % e_overlap
  min_info %e_over_ii = thermo % e_overlap_ii
  min_info %e_over_ij = thermo % e_overlap_ij
  min_info %e_electro = thermo % e_electro
  min_info %e_electro_0 = thermo % e_electro_0
  min_info %e_electro_1 = thermo % e_electro_1
  min_info %e_electro_2 = thermo % e_electro_2
  min_info %e_thk = thermo % e_thk
  min_info %e_pol = thermo % e_pol
  min_info %de = thermo % pot - eold
  min_info %tol = tol
  min_info %err = fconv
 END  SUBROUTINE copy_energy_info 
!------------------------------------------------------------------------------!
    END MODULE energy_optimize
!--------------------------------------------------------------------------------!
