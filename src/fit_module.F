!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart               !
!------------------------------------------------------------------------------!
    MODULE fit_module
!-----------------------------------------------------------------------------!
!!****** cp2k/fit_module [1.0] *
!!
!!   NAME
!!     fit_module
!!
!!   FUNCTION
!!
!!   enables to do different types of fit (rigid_ion, polarizable ..)
!!   using the amoeba procedire
!!
!*****************************************************************
!------------------------------------------------------------------------------!
!
  USE amoeba, ONLY : amoeba_evaluate
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type, thermodynamic_type
  USE pair_potential, ONLY : potentialparm_type, spline_nonbond_control
  USE particle_types, ONLY : particle_prop_type 
  USE structure_types, ONLY : structure_type
  USE termination, ONLY : stop_program, stop_memory
  USE fit_parameters_types, ONLY : fit_parameter_type,  fit_averages_type, &
      linklist_parameters_type, fit_info_type
  USE fit_utilities, ONLY : force_ref, fit_print_info, get_will_params, &
                            get_lj_params, fit_read_restart, fit_write_restart
  USE fit_types, ONLY : copy_params, initialize_params
  USE simulation_cell, ONLY :  get_hinv, cell_type      
  USE fit_force_control, ONLY : fit_force
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  
!gt this is to fix to do the polarizable fit  
!      USE force_control, ONLY : force
!      USE basis_kind_info, ONLY : kind_info_type
!      USE minimize, ONLY : fit_min_control, dump_coef  => dump      
!      USE pol_lists, ONLY : pol_es_lists
!      USE pol_electrostatic_force, ONLY : force_electrostatic
!      USE pol_overlap_force, ONLY : force_overlap
!      USE pol_classic_force, ONLY : force_classic

  IMPLICIT NONE
!
  PUBLIC :: control_fit, local
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
  SUBROUTINE control_fit(struc,potparm,simpar,fit_param, pstat, &
                         ewald_param, atom_names,globenv)
    IMPLICIT NONE
    TYPE (structure_type), INTENT (INOUT) :: struc
    TYPE (potentialparm_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: potparm 
    TYPE (simulation_parameters_type), INTENT ( IN ) :: simpar      
    TYPE (particle_prop_type ), DIMENSION ( : ), POINTER :: pstat 
    TYPE (fit_parameter_type), INTENT(INOUT) :: fit_param
    TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param  
    TYPE ( global_environment_type ), INTENT ( IN ) :: globenv  
    CHARACTER ( LEN = * ), DIMENSION ( : ), POINTER :: atom_names 

!
    INTEGER :: isos
! locals
    TYPE ( thermodynamic_type ) :: thermo
    TYPE ( fit_averages_type ) :: fit_ave
    TYPE ( linklist_parameters_type ), POINTER, DIMENSION ( :, : ) :: param   
    LOGICAL :: conv_flag
!
    SELECT CASE (fit_param%fit_type)
!
! before fitting short-range parameters to the cp forces, correct them 
! by subtracting the polarization contribution (and Coulomb if present)
!
    CASE ('FORCES_POL')

! gt to be implemanted
!         CALL calculate_polar_forces(struc,thermo,simpar,potparm)
! gt to be implemanted
!       
!       amoeba fit on the cpmd forces , or on the forces obtained
!       by subtracting the polarizable forces to the cpmd ones.
!       This works both for rigid ions model and polarizable force fields
!       (no variation of the coefficients, only amoeba minimizer)
!
    CASE ('AMOEBA_RI')

! copy all relevent variables to local, global variables

      ALLOCATE (param(size(potparm,1),size(potparm,2)),STAT=isos)

      IF (isos/=0) CALL stop_memory('fit_control','param')

! allocating the arrays for the analysis

      SELECT CASE (fit_param%pot_type)
      
      CASE ('LJ')
        ALLOCATE (fit_ave % ave_sig(size(potparm,1),size(potparm,2)),STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'control_fit', 'ave_sig')

        ALLOCATE (fit_ave % ave_eps(size(potparm,1),size(potparm,2)),STAT=isos)
        IF (isos/=0) CALL stop_memory('control_fit', 'ave_eps')
   
      CASE ('WILL')

        ALLOCATE (fit_ave % ave_a(size(potparm,1),size(potparm,2)),STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'control_fit', 'ave_a')

        ALLOCATE (fit_ave % ave_b(size(potparm,1),size(potparm,2)),STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'control_fit', 'ave_b')
     
        ALLOCATE (fit_ave % ave_c(size(potparm,1),size(potparm,2)),STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'control_fit', 'ave_c')
      
      END SELECT

      conv_flag = .FALSE.

      CALL fit(struc,potparm,thermo,simpar,fit_param,param,pstat,fit_ave, &
                ewald_param,atom_names,globenv,conv_flag)

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
          IF ( isos /= 0 ) CALL stop_memory ( 'control_fit', 'deallocate ave_sig')
        END IF   

        IF ( ASSOCIATED ( fit_ave % ave_eps ) ) THEN   
          DEALLOCATE (fit_ave % ave_eps , STAT=isos)
          IF ( isos /= 0 ) CALL stop_memory ( 'control_fit', 'deallocate ave_eps')
        END IF   
   
      CASE ('WILL')

        IF ( ASSOCIATED ( fit_ave % ave_a ) ) THEN   
          DEALLOCATE (fit_ave % ave_a , STAT=isos)
          IF ( isos /= 0 ) CALL stop_memory ( 'control_fit', 'deallocate ave_a' )
        END IF  

        IF ( ASSOCIATED ( fit_ave % ave_b ) ) THEN   
          DEALLOCATE (fit_ave % ave_b , STAT=isos)
          IF ( isos /= 0 ) CALL stop_memory ( 'control_fit', 'deallocate ave_b' )
        END IF  

        IF ( ASSOCIATED ( fit_ave % ave_c ) ) THEN   
          DEALLOCATE (fit_ave % ave_c , STAT=isos)
          IF ( isos /= 0 ) CALL stop_memory ( 'control_fit', 'deallocate ave_c' )
        END IF  
   
      END SELECT

      IF ( ASSOCIATED ( param ) ) THEN   
        DEALLOCATE ( param , STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'control_fit', 'deallocate param' )
      END IF  

! to perform amoeba fit in the case of polarizable force fields
! in which the amoeba parameters can also depend on the coefficents
! TO BE IMPLEMENTED!!!!
!        CASE ('AMOEBA_POL')

    END SELECT

  END SUBROUTINE control_fit
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE fit(struc,potparm,thermo,simpar,fit_param,param,pstat,fit_ave,&
                 ewald_param,atom_names,globenv,conv_flag)

! Given the reference force, and the derived force, optimize
! paramters with a least squares fit.

  IMPLICIT NONE

  TYPE (structure_type), INTENT (INOUT) :: struc
  TYPE (potentialparm_type), DIMENSION (:,:), INTENT (INOUT) :: potparm
  TYPE (thermodynamic_type), INTENT (INOUT) :: thermo
  TYPE (simulation_parameters_type), INTENT (IN) :: simpar
  TYPE (fit_parameter_type), INTENT(INOUT) :: fit_param
  TYPE (linklist_parameters_type ), INTENT(INOUT), DIMENSION ( :, : ) :: param   
  TYPE ( particle_prop_type ), DIMENSION ( : ), POINTER :: pstat 
  TYPE (fit_averages_type) :: fit_ave
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param  
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv  
  CHARACTER ( LEN = * ), DIMENSION ( : ), POINTER :: atom_names 
  LOGICAL, INTENT (INOUT) :: conv_flag

! locals

  INTEGER :: i, j, iter
  REAL (dbl) :: tol
  REAL ( dbl ), DIMENSION ( : ), POINTER :: y
  REAL ( dbl ), DIMENSION ( :, : ), POINTER :: p
  TYPE (fit_info_type), DIMENSION (:,:), POINTER :: finfo
  TYPE (cell_type) , POINTER, DIMENSION (:) :: box_ref


  INTERFACE
    
    FUNCTION func(x,struc,potparm,thermo,simpar,fit_param,box_ref, &
                  finfo,param,pstat,ewald_param, globenv)

      USE kinds, ONLY : dbl
      USE structure_types, ONLY : structure_type
      USE pair_potential, ONLY : potentialparm_type
      USE particle_types, ONLY : particle_prop_type 
      USE md, ONLY : simulation_parameters_type, thermodynamic_type
      USE simulation_cell, ONLY : cell_type
      USE fit_parameters_types, ONLY : fit_parameter_type, linklist_parameters_type, &
                                       fit_info_type
      USE ewald_parameters_types, ONLY : ewald_parameters_type
      USE global_types, ONLY : global_environment_type
      IMPLICIT NONE
      TYPE (structure_type), INTENT (INOUT) :: struc
      TYPE (potentialparm_type), DIMENSION (:,:), INTENT (INOUT) :: potparm
      TYPE (thermodynamic_type), INTENT (INOUT) :: thermo
      TYPE (simulation_parameters_type), INTENT (IN) :: simpar
      REAL (dbl), INTENT (IN), DIMENSION (:) :: x
      TYPE (fit_parameter_type), INTENT(INOUT) :: fit_param
      TYPE (cell_type) , POINTER, DIMENSION (:) :: box_ref
      TYPE (fit_info_type), DIMENSION (:,:), POINTER :: finfo
      TYPE (linklist_parameters_type ), INTENT(INOUT), DIMENSION ( :, : ) :: param 
      TYPE (particle_prop_type ), DIMENSION ( : ), POINTER :: pstat 
      TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param  
      TYPE ( global_environment_type), INTENT ( IN ) :: globenv

      REAL (dbl) :: func

    END FUNCTION func

  END INTERFACE

! Allocate and intialize arrays for amoeba

  SELECT CASE (fit_param%fit_type)

  CASE ('AMOEBA_RI')

    CALL initialize_params(fit_param, potparm, param, y, p )

  CASE ('AMOEBA_POL')

    CALL stop_program("fit", "not implemanted")

  END SELECT  

! Get reference force

  CALL force_ref(struc,potparm,simpar,fit_param,finfo,box_ref)

! initialize y
  SELECT CASE ( fit_param % restart )
  CASE ( 'ALL' )

     CALL fit_read_restart ( y, p )

  CASE ( 'INITIAL' )

     DO i = 1, SIZE ( y )
        y ( i ) = func ( p ( i, : ), struc, potparm, thermo, simpar, &
          fit_param, box_ref, finfo, param, pstat, ewald_param, globenv )
     END DO

  CASE DEFAULT
     CALL stop_program ( "fit", "no valid restart option" )

  END SELECT

! optimize

  DO iter = 1, fit_param%max_iter, fit_param%iwrite
  
    CALL amoeba_evaluate(p,y,tol,struc, thermo, potparm, simpar, fit_param, &
            box_ref, finfo, param, pstat, ewald_param, globenv, func )

    SELECT CASE (fit_param%fit_type)

    CASE ('AMOEBA_RI')
      SELECT CASE (fit_param%pot_type)
      CASE ('LJ')
        CALL get_lj_params(potparm, y, param, fit_ave%ave_sig, fit_ave%ave_eps)
      CASE ('WILL')
        CALL get_will_params(potparm, y, param, fit_ave% ave_a, &
                              fit_ave%ave_b, fit_ave%ave_c)
      END SELECT
      CALL fit_print_info(tol,potparm,fit_param,y,atom_names, &
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
  END SUBROUTINE fit
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE forcefit(struc, potparm, thermo,box_ref,finfo, &
                      fit_param, ewald_param, globenv)

  IMPLICIT NONE
  TYPE (structure_type), INTENT (INOUT) :: struc
  TYPE (potentialparm_type), DIMENSION (:,:), INTENT (INOUT) :: potparm
  TYPE (thermodynamic_type), INTENT (INOUT) :: thermo
  TYPE (cell_type) , POINTER, DIMENSION (:) :: box_ref
  TYPE (fit_info_type), DIMENSION (:,:), POINTER :: finfo
  TYPE (fit_parameter_type), INTENT(INOUT) :: fit_param
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param  
  TYPE ( global_environment_type), INTENT ( IN ) :: globenv

!locals

  INTEGER :: itime, ipart, npart
  
! From the reference trajectory file obtain the fit force
! using a CP2K force field.

  npart = size(struc%pnode)

  DO itime = 1, fit_param%nsteps

    DO ipart = 1, npart
      struc%pnode(ipart) %p%r(:) = finfo(itime,ipart) %r(:)
    END DO

!pass the box  

    struc%box=box_ref(itime)

    CALL fit_force ( struc, potparm, thermo, &
                           ewald_param, fit_param, globenv )
    
    DO ipart = 1, npart
      finfo(itime,ipart) %f_fit(:) = struc%pnode(ipart) %p%f(:)
    END DO
  END DO
  END SUBROUTINE forcefit
!
!------------------------------------------------------------------------------!
!
  FUNCTION local(x,struc,potparm,thermo,simpar,fit_param,box_ref,finfo, &
                 param, pstat, ewald_param, globenv )

  IMPLICIT NONE
  REAL (dbl), DIMENSION (:), INTENT (IN) :: x
  TYPE (structure_type), INTENT (INOUT) :: struc
  TYPE (thermodynamic_type), INTENT (INOUT) :: thermo
  TYPE (simulation_parameters_type), INTENT (IN) :: simpar
  TYPE (potentialparm_type), DIMENSION (:,:), INTENT (INOUT) :: potparm
  TYPE (fit_parameter_type), INTENT(INOUT) :: fit_param
  TYPE (cell_type) , POINTER, DIMENSION (:) :: box_ref
  TYPE (fit_info_type), DIMENSION (:,:), POINTER :: finfo
  TYPE (linklist_parameters_type ), INTENT(INOUT), DIMENSION ( :, : ) :: param 
  TYPE (particle_prop_type ), DIMENSION ( : ), POINTER :: pstat 
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param  
  TYPE ( global_environment_type), INTENT ( IN ) :: globenv


! locals
  REAL (dbl), DIMENSION (3) :: vec
  REAL (dbl) :: local
  INTEGER :: itime, ipart, npart

  npart = size(struc%pnode)
  
  SELECT CASE (fit_param%fit_type) 
  
  CASE ('AMOEBA_RI')
  
    CALL copy_params(fit_param, potparm , pstat, param, ewald_param, x)

    CALL forcefit(struc,potparm,thermo,box_ref,finfo, &
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
    END MODULE fit_module
!------------------------------------------------------------------------------!
