!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart               !
!------------------------------------------------------------------------------!
    MODULE fit_pair_potential
!-----------------------------------------------------------------------------!
!!****** cp2k/fit_pair_potential [1.0] *
!!
!!   NAME
!!     fit_pair_potential
!!
!!   FUNCTION
!!
!*****************************************************************
!------------------------------------------------------------------------------!
!
  USE amoeba,                          ONLY: amoeba_evaluate
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE fit_force_control,               ONLY: fit_force
  USE fit_parameters_types,            ONLY: fit_averages_type,&
                                             fit_info_type,&
                                             fit_parameter_type,&
                                             linklist_parameters_type
  USE fit_types,                       ONLY: copy_params,&
                                             initialize_params
  USE fit_utilities,                   ONLY: fit_print_info,&
                                             fit_read_restart,&
                                             fit_write_restart,&
                                             get_lj_params,&
                                             get_will_params
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE md,                              ONLY: simulation_parameters_type,&
                                             thermodynamic_type
  USE pair_potential,                  ONLY: potentialparm_type
  USE particle_types,                  ONLY: particle_prop_type
  USE simulation_cell,                 ONLY: cell_type
  USE structure_types,                 ONLY: interaction_type,&
                                             structure_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE
!
  PUBLIC :: fit_potparm, local
  PRIVATE
!
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!
!!>----------------------------------------------------------------------------!
!! SECTION: &fit ... &end                                                      !
!!                                                                             !
!!  type:    [LJ, WILL]                  type of fit                     !
!!  restart: [INITIAL, ALL]              restart option                        !
!!  data type: [POS_FORCE, POS_VEL]      type of data to read                  !
!!  out unit: [ATOMIC, KELVIN]           output unit                           !
!!  tol:      max_tol                    maximum tolerance (fractional)        !
!!  iter:     max_iter                   maximum number of iterations          !
!!  runtype: [AMOEBA_RI, FORCES_POL          
!!  nsteps:   ntime                      number of trajectory steps            !
!!  iskip:    iskip                      analyse every iskip step              !
!!  idump:    idump                      write to dump file every idump steps  !
!!                                                                             !
!!                                                                             !
!!<----------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!
  SUBROUTINE fit_potparm(struc,inter,simpar,fit_param, pstat, &
                         ewald_param,box_ref,finfo, atom_names,globenv)
    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(interaction_type), intent(inout)    :: inter
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(fit_parameter_type), INTENT(INOUT)  :: fit_param
    TYPE(particle_prop_type), DIMENSION(:), &
      POINTER                                :: pstat
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    TYPE(cell_type), DIMENSION(:), POINTER   :: box_ref
    TYPE(fit_info_type), DIMENSION(:, :), &
      POINTER                                :: finfo
    CHARACTER(LEN=*), DIMENSION(:), POINTER  :: atom_names
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: isos
    LOGICAL                                  :: conv_flag
    TYPE(fit_averages_type)                  :: fit_ave
    TYPE(linklist_parameters_type), &
      DIMENSION(:, :), POINTER               :: param
    TYPE(potentialparm_type), &
      DIMENSION(:, :), POINTER               :: potparm
    TYPE(thermodynamic_type)                 :: thermo

!    TYPE (potentialparm_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: potparm 
!

    potparm => inter % potparm
!
! copy all relevent variables to local, global variables

    ALLOCATE (param(size(potparm,1),size(potparm,2)),STAT=isos)

    IF (isos/=0) CALL stop_memory('fit_potparm','param')

! allocating the arrays for the analysis

    SELECT CASE (fit_param%pot_type)
      
    CASE ('LJ')
      ALLOCATE (fit_ave % ave_sig(size(potparm,1),size(potparm,2)),STAT=isos)
      IF ( isos /= 0 ) CALL stop_memory ( 'fit_potparm', 'ave_sig')

      ALLOCATE (fit_ave % ave_eps(size(potparm,1),size(potparm,2)),STAT=isos)
      IF (isos/=0) CALL stop_memory('fit_potparm', 'ave_eps')
   
    CASE ('WILL')

      ALLOCATE (fit_ave % ave_a(size(potparm,1),size(potparm,2)),STAT=isos)
      IF ( isos /= 0 ) CALL stop_memory ( 'fit_potparm', 'ave_a')

      ALLOCATE (fit_ave % ave_b(size(potparm,1),size(potparm,2)),STAT=isos)
      IF ( isos /= 0 ) CALL stop_memory ( 'fit_potparm', 'ave_b')
     
      ALLOCATE (fit_ave % ave_c(size(potparm,1),size(potparm,2)),STAT=isos)
      IF ( isos /= 0 ) CALL stop_memory ( 'fit_potparm', 'ave_c')
      
    END SELECT

    conv_flag = .FALSE.

    CALL do_fit(struc,inter,thermo,simpar,fit_param,param,pstat,fit_ave, &
              ewald_param,atom_names,box_ref,finfo,globenv,conv_flag)

    IF (conv_flag) THEN

      WRITE (*,*)
      WRITE (*,'(T23,A)') '******************************'
      WRITE (*,'(T23,A)') '*                            *'
      WRITE (*,'(T23,A)') '*  FIT|CONVERGENCE ACHIEVED  *'
      WRITE (*,'(T23,A)') '*                            *'
      WRITE (*,'(T23,A)') '******************************'
      WRITE (*,*)

    ELSE

      WRITE (*,*)
      WRITE (*,'(T23,A)') '******************************'
      WRITE (*,'(T23,A)') '*                            *'
      WRITE (*,'(T23,A)') '* FIT|WARNING:NO CONVERGENCE *'
      WRITE (*,'(T23,A)') '*                            *'
      WRITE (*,'(T23,A)') '******************************'
      WRITE (*,*)

    END IF

! deallocating the arrays for the analysis

    SELECT CASE (fit_param%pot_type)

    CASE ('LJ')

      IF ( ASSOCIATED ( fit_ave % ave_sig ) ) THEN   
        DEALLOCATE (fit_ave % ave_sig , STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'fit_potparm', 'deallocate ave_sig')
      END IF   

      IF ( ASSOCIATED ( fit_ave % ave_eps ) ) THEN   
        DEALLOCATE (fit_ave % ave_eps , STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'fit_potparm', 'deallocate ave_eps')
      END IF   
   
    CASE ('WILL')

      IF ( ASSOCIATED ( fit_ave % ave_a ) ) THEN   
        DEALLOCATE (fit_ave % ave_a , STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'fit_potparm', 'deallocate ave_a' )
      END IF  

      IF ( ASSOCIATED ( fit_ave % ave_b ) ) THEN   
        DEALLOCATE (fit_ave % ave_b , STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'fit_potparm', 'deallocate ave_b' )
      END IF  

      IF ( ASSOCIATED ( fit_ave % ave_c ) ) THEN   
        DEALLOCATE (fit_ave % ave_c , STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'fit_potparm', 'deallocate ave_c' )
      END IF  
 
    END SELECT

    IF ( ASSOCIATED ( param ) ) THEN   
      DEALLOCATE ( param , STAT=isos)
      IF ( isos /= 0 ) CALL stop_memory ( 'fit_potparm', 'deallocate param' )
    END IF  

  END SUBROUTINE fit_potparm
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE do_fit(struc,inter,thermo,simpar,fit_param,param,pstat,fit_ave,&
                 ewald_param,atom_names,box_ref,finfo,globenv,conv_flag)

! Given the reference force, and the derived force, optimize
! paramters with a least squares fit.


    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(interaction_type), INTENT(INOUT)    :: inter
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(fit_parameter_type), INTENT(INOUT)  :: fit_param
    TYPE(linklist_parameters_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: param
    TYPE(particle_prop_type), DIMENSION(:), &
      POINTER                                :: pstat
    TYPE(fit_averages_type)                  :: fit_ave
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    CHARACTER(LEN=*), DIMENSION(:), POINTER  :: atom_names
    TYPE(cell_type), DIMENSION(:), POINTER   :: box_ref
    TYPE(fit_info_type), DIMENSION(:, :), &
      POINTER                                :: finfo
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    LOGICAL, INTENT(INOUT)                   :: conv_flag

    INTEGER                                  :: i, iter
    REAL(dbl)                                :: tol
    REAL(dbl), DIMENSION(:), POINTER         :: y
    REAL(dbl), DIMENSION(:, :), POINTER      :: p
  INTERFACE
    
    FUNCTION func(x,struc,inter,thermo,simpar,fit_param,box_ref, &
                  finfo,param,pstat,ewald_param, globenv)

      USE kinds, ONLY : dbl
      USE structure_types, ONLY : structure_type, interaction_type
      USE particle_types, ONLY : particle_prop_type 
      USE md, ONLY : simulation_parameters_type, thermodynamic_type
      USE simulation_cell, ONLY : cell_type
      USE fit_parameters_types, ONLY : fit_parameter_type, linklist_parameters_type, &
                                       fit_info_type
      USE ewald_parameters_types, ONLY : ewald_parameters_type
      USE global_types, ONLY : global_environment_type
      IMPLICIT NONE
      TYPE (structure_type), INTENT (INOUT) :: struc
      TYPE (interaction_type), INTENT (INOUT) :: inter
      TYPE (thermodynamic_type), INTENT (INOUT) :: thermo
      TYPE (simulation_parameters_type), INTENT (IN) :: simpar
      REAL (dbl), INTENT (IN), DIMENSION (:) :: x
      TYPE (fit_parameter_type), INTENT(IN) :: fit_param
      TYPE (cell_type) , POINTER, DIMENSION (:) :: box_ref
      TYPE (fit_info_type), DIMENSION (:,:), POINTER :: finfo
      TYPE (linklist_parameters_type ), INTENT(INOUT), DIMENSION ( :, : ) :: param 
      TYPE (particle_prop_type ), DIMENSION ( : ), POINTER :: pstat 
      TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param  
      TYPE ( global_environment_type), INTENT ( IN ) :: globenv

      REAL (dbl) :: func

    END FUNCTION func

  END INTERFACE

!!  TYPE (potentialparm_type), DIMENSION (:,:), INTENT (INOUT) :: potparm
! Allocate and intialize arrays for amoeba

  SELECT CASE (fit_param%fit_type)

  CASE ('AMOEBA_RI')

    CALL initialize_params(fit_param, inter % potparm, param, y, p )

  CASE ('AMOEBA_POL')

    CALL stop_program("fit", "not implemanted")

  END SELECT  

! initialize y
  SELECT CASE ( fit_param % restart )
  CASE ( 'ALL' )

     CALL fit_read_restart ( y, p )

  CASE ( 'INITIAL' )

     DO i = 1, SIZE ( y )
        y ( i ) = func ( p ( i, : ), struc, inter, thermo, simpar, &
          fit_param, box_ref, finfo, param, pstat, ewald_param, globenv )
     END DO

  CASE DEFAULT
     CALL stop_program ( "fit", "no valid restart option" )

  END SELECT

! optimize

  DO iter = 1, fit_param%max_iter, fit_param%iwrite
  
    CALL amoeba_evaluate(p,y,tol,struc, thermo, inter, simpar, fit_param, &
            box_ref, finfo, param, pstat, ewald_param, globenv, func )

    SELECT CASE (fit_param%fit_type)

    CASE ('AMOEBA_RI')
      SELECT CASE (fit_param%pot_type)
      CASE ('LJ')
        CALL get_lj_params(inter % potparm, y, param, fit_ave%ave_sig, fit_ave%ave_eps)
      CASE ('WILL')
        CALL get_will_params(inter % potparm, y, param, fit_ave% ave_a, &
                              fit_ave%ave_b, fit_ave%ave_c)
      END SELECT
      CALL fit_print_info(tol,inter % potparm,fit_param,y,atom_names, &
                           fit_ave,globenv,iter)
    CASE ('AMOEBA_POL')
    END SELECT
    
! write restart file
     IF ( MOD ( iter, fit_param % idump ) == 0 ) &
          CALL fit_write_restart ( y, p )

! convergence acheived
    IF (tol<=fit_param%max_tol) THEN
    
      conv_flag = .TRUE.
        CALL fit_write_restart ( y, p )
      
      EXIT
      
    END IF
  END DO
  END SUBROUTINE do_fit
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE forcefit(struc, inter, thermo, simpar,box_ref,finfo, &
                      fit_param, ewald_param, globenv)

    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(interaction_type), INTENT(INOUT)    :: inter
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(cell_type), DIMENSION(:), POINTER   :: box_ref
    TYPE(fit_info_type), DIMENSION(:, :), &
      POINTER                                :: finfo
    TYPE(fit_parameter_type), INTENT(INOUT)  :: fit_param
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: ii, ipart, itime, jj, npart
    LOGICAL                                  :: box_change
    REAL(dbl)                                :: boxtol, deltabox

! From the reference trajectory file obtain the fit force
! using a CP2K force field.

  deltabox = 0.0_dbl
  npart = size(struc%pnode)
  box_change = .FALSE.
  boxtol = 1.E-6_dbl

  DO itime = 1, fit_param%nsteps

    box_change = .FALSE.
    DO ipart = 1, npart
      struc%pnode(ipart) %p%r(:) = finfo(itime,ipart) %r(:)
    END DO

!check if the box changes
out: DO ii = 1,3
inn:   DO JJ = 1,3
          deltabox = box_ref (itime) % hmat(ii,jj) - struc % box % hmat (ii,jj)
           IF (ABS(deltabox) >= boxtol ) THEN
              box_change = .TRUE.
              EXIT out
           END IF
        END DO inn
     END DO out

!pass the box  
    IF ( box_change ) struc %box=box_ref(itime)

    CALL fit_force ( struc, inter, thermo, simpar, &
                           ewald_param, box_change, fit_param, globenv )
    
    DO ipart = 1, npart
      finfo(itime,ipart) %f_fit(:) = struc%pnode(ipart) %p%f(:)
    END DO
  END DO
  END SUBROUTINE forcefit
!
!------------------------------------------------------------------------------!
!
  FUNCTION local(x,struc,inter,thermo,simpar,fit_param,box_ref,finfo, &
                 param, pstat, ewald_param, globenv )

    REAL(dbl), DIMENSION(:), INTENT(IN)      :: x
    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(interaction_type), INTENT(INOUT)    :: inter
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(fit_parameter_type), INTENT(INOUT)  :: fit_param
    TYPE(cell_type), DIMENSION(:), POINTER   :: box_ref
    TYPE(fit_info_type), DIMENSION(:, :), &
      POINTER                                :: finfo
    TYPE(linklist_parameters_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: param
    TYPE(particle_prop_type), DIMENSION(:), &
      POINTER                                :: pstat
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    REAL(dbl)                                :: local

    INTEGER                                  :: ipart, itime, npart
    REAL(dbl), DIMENSION(3)                  :: vec

  npart = size(struc%pnode)
  
  SELECT CASE (fit_param%fit_type) 
  
  CASE ('AMOEBA_RI')
  
    CALL copy_params(fit_param, inter % potparm , pstat, param, ewald_param, x)

    CALL forcefit(struc,inter,thermo, simpar, box_ref,finfo, &
                   fit_param, ewald_param, globenv)

  CASE default    
  
    CALL stop_program('local', 'not implemented')

  END SELECT

  local =  0._dbl
  
  DO itime = 1, fit_param%nsteps
    DO ipart = 1, npart
      vec(:) = finfo(itime,ipart) %f_fit(:)  - finfo(itime,ipart) %f_ref(:)
      local = local + dot_product(vec,vec)
    END DO
  END DO
  
! The addition of the 1 is a safety mechanism so
! that one can define the FRACTIONAL tolerance for
! systems whos func value may be converged at zero.
! It does not effect other converged func values.
! dbg
!        local = local/(3._dbl*npart*fitpar%ntime) + 1._dbl
! dbg

  local=local/(3._dbl*npart*fit_param%nsteps)

  END FUNCTION local
!
!------------------------------------------------------------------------------!
!
    END MODULE fit_pair_potential
!------------------------------------------------------------------------------!
