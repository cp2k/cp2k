!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart               !
!------------------------------------------------------------------------------!
    MODULE fit_types
!-----------------------------------------------------------------------------!
!!****** cp2k/fit_types [1.0] *
!!
!!   NAME
!!     fit_types
!!
!!   FUNCTION
!!
!!   prepares the array for the amoeba for different types of potentials
!!
!*****************************************************************
!------------------------------------------------------------------------------!
!
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE fit_parameters_types,            ONLY: fit_parameter_type,&
                                             linklist_parameters_type,&
                                             linklist_real_type
  USE kinds,                           ONLY: dp
  USE pair_potential,                  ONLY: potentialparm_type,&
                                             spline_nonbond_control
  USE particle_types,                  ONLY: particle_prop_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE
!
  PUBLIC :: copy_params, initialize_params
  PRIVATE
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE copy_params(fit_param, potparm, pstat, param, ewald_param, x)

    TYPE(fit_parameter_type), INTENT(IN)     :: fit_param
    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: potparm
    TYPE(particle_prop_type), DIMENSION(:), &
      POINTER                                :: pstat
    TYPE(linklist_parameters_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: param
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: x

!------------------------------------------------------------------------------!

  SELECT CASE (fit_param%pot_type)

  CASE ('LJ')

     CALL lj_fit(potparm,pstat,param,ewald_param,x)

  CASE ('WILL')

    CALL will_fit(potparm,pstat,param,ewald_param,x)

  CASE DEFAULT 
  
    CALL stop_program ( "copy_params", "no valid potential type" )
    
  END SELECT

  END SUBROUTINE copy_params
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE initialize_params(fit_param, potparm, param, y, p)

    TYPE(fit_parameter_type), INTENT(IN)     :: fit_param
    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: potparm
    TYPE(linklist_parameters_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: param
    REAL(KIND=dp), DIMENSION(:), POINTER         :: y
    REAL(KIND=dp), DIMENSION(:, :), POINTER      :: p

!------------------------------------------------------------------------------!

  SELECT CASE (fit_param%pot_type)

  CASE ('LJ')

     CALL lj_initialize(potparm,param,y,p)

  CASE ('WILL')

    CALL will_initialize(potparm,param,y,p)

  CASE DEFAULT 
  
    CALL stop_program ( "copy_params", "no valid potential type" )
    
  END SELECT

  END SUBROUTINE initialize_params
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE lj_fit(potparm,pstat,param,ewald_param,x)


    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: potparm
    TYPE(particle_prop_type), DIMENSION(:), &
      POINTER                                :: pstat
    TYPE(linklist_parameters_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: param
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: x

    INTEGER                                  :: iatom, icount, jatom

    icount = 0
    
    DO iatom = 1, size(potparm,1)
    
      DO jatom = iatom, size(potparm,2)
      
        icount = icount + 1

        potparm(iatom,jatom) %lj%epsilon = x(icount)*x(icount)
! symmetrizing
        potparm(jatom,iatom) %lj%epsilon = potparm(iatom,jatom) &
          %lj%epsilon
        icount = icount + 1
        potparm(iatom,jatom) %lj%sigma6 = x(icount)*x(icount)
        potparm(iatom,jatom) %lj%sigma12 = x(icount)*x(icount) &
                                                    *x(icount)*x(icount)
! symmetrizing
        potparm(jatom,iatom) %lj%sigma6 = potparm(iatom,jatom) &
          %lj%sigma6
        potparm(jatom,iatom) %lj%sigma12 = potparm(iatom,jatom) &
          %lj%sigma12
      END DO

    END DO

! initialize splines

    CALL spline_nonbond_control(potparm, pstat, 2500, ewald_param)

  END SUBROUTINE lj_fit
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE will_fit(potparm,pstat,param,ewald_param,x)

    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: potparm
    TYPE(particle_prop_type), DIMENSION(:), &
      POINTER                                :: pstat
    TYPE(linklist_parameters_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: param
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: x

    INTEGER                                  :: iatom, icount, jatom

    icount = 0

    DO iatom = 1, size(potparm,1)
      DO jatom = iatom, size(potparm,2)
        icount = icount + 1
        potparm(iatom,jatom) %willis%a = x(icount)*x(icount)
! symmetrizing
        potparm(jatom,iatom) %willis%a = potparm(iatom,jatom) &
          %willis%a
        icount = icount + 1
        potparm(iatom,jatom) %willis%b = x(icount)*x(icount)
! symmetrizing
        potparm(jatom,iatom) %willis%b = potparm(iatom,jatom) &
          %willis%b
        icount = icount + 1
        potparm(iatom,jatom) %willis%c = x(icount)*x(icount)
! symmetrizing
        potparm(jatom,iatom) %willis%c = potparm(iatom,jatom) &
          %willis%c
      END DO
    END DO
    
! initialize splines

    CALL spline_nonbond_control(potparm,pstat,2500,ewald_param)

  END SUBROUTINE will_fit
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE lj_initialize(potparm,param,y,p)


    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: potparm
    TYPE(linklist_parameters_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: param
    REAL(KIND=dp), DIMENSION(:), POINTER         :: y
    REAL(KIND=dp), DIMENSION(:, :), POINTER      :: p

    INTEGER                                  :: i, iatom, icount, isos, &
                                                jatom, numpar
    REAL(KIND=dp)                                :: eps, sig
    TYPE(linklist_real_type), POINTER        :: lleps, llsig

! number of potential parameters to be optimized is 2*number 
! of interactions (remember: we only need the upper triangle)
! e.g. for LJ each interaction has 2 parameters, sigma and epsilon

    numpar = 2*size(potparm,1)*(size(potparm,1)+1)/2

    ALLOCATE (y(numpar+1),STAT=isos)
    IF (isos/=0) CALL stop_memory('lj_initialize','alloc y')
    
    ALLOCATE (p(numpar+1,numpar),STAT=isos)

    IF (isos/=0) CALL stop_memory('lj_initialize','alloc p')

! i=0,numpar is the numpar+1 vector in p(i,j).  It is also
! used to generate the initial guesses for the verticies
! of the amoeba

    ALLOCATE (lleps,STAT=isos)
    IF (isos/=0) CALL stop_memory('lj_initialize', 'lleps')

    ALLOCATE (llsig,STAT=isos)
    IF (isos/=0) CALL stop_memory('lj_initialize', 'llsig')
     
    icount = 0

    DO iatom = 1, size(potparm,1)

      DO jatom = iatom, size(potparm,2)

        param(iatom,jatom) %lleps => lleps
        param(iatom,jatom) %llsig => llsig

        icount = icount + 1
! loop for eps

        DO i = 1, numpar + 1
          eps=potparm(iatom,jatom) %lj%epsilon
          p(i,icount) = SQRT(eps)
          lleps%real_value => p(i,icount)
          ALLOCATE (lleps%next,STAT=isos)
          IF (isos/=0) CALL stop_memory('lj_initialize', 'lleps%next' )
          lleps => lleps%next
        END DO

        icount = icount + 1

! loop for sig

        DO i = 1, numpar + 1
          sig=potparm(iatom,jatom) %lj%sigma6
          p(i,icount) = SQRT(sig)
          llsig%real_value => p(i,icount)
          ALLOCATE (llsig%next,STAT=isos)
          IF (isos/=0) CALL stop_memory('lj_initialize', 'llsig%next')
          llsig => llsig%next
        END DO
      END DO
    END DO

! scaling the appropriate interaction parameter by 1.1
! to construct the amoeba

    DO i = 0, numpar
      IF (i>0) p(i+1,i) = p(i+1,i)*1.1_dp
    END DO

  END SUBROUTINE lj_initialize
!
!------------------------------------------------------------------------------!
!
  SUBROUTINE will_initialize(potparm,param,y,p)

    TYPE(potentialparm_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: potparm
    TYPE(linklist_parameters_type), &
      DIMENSION(:, :), INTENT(INOUT)         :: param
    REAL(KIND=dp), DIMENSION(:), POINTER         :: y
    REAL(KIND=dp), DIMENSION(:, :), POINTER      :: p

    INTEGER                                  :: i, iatom, icount, isos, &
                                                jatom, numpar
    REAL(KIND=dp)                                :: a, b, c
    TYPE(linklist_real_type), POINTER        :: lla, llb, llc

! of interactions (remember: we only need the upper triangle)
! e.g. for WILL each interaction has 3 parameters, a, b, and c.

    numpar =  3_dp*size(potparm,1)*(size(potparm,1)+1)/2
    
    ALLOCATE (y(numpar+1),STAT=isos)
    IF (isos/=0) CALL stop_memory('will_initialize','alloc y')

    ALLOCATE (p(numpar+1,numpar),STAT=isos)
    IF (isos/=0) CALL stop_memory('will_initialize','alloc p')
    
! i=0,numpar is the numpar+1 vector in p(i,j).  It is also
! used to generate the initial guesses for the verticies
! of the amoeba

    ALLOCATE (lla,STAT=isos)
    IF (isos/=0) CALL stop_memory('will_initialize','alloc lla')

    ALLOCATE (llb,STAT=isos)
    IF (isos/=0) CALL stop_memory('will_initialize','alloc llb')
    
    ALLOCATE (llc,STAT=isos)
    IF (isos/=0) CALL stop_memory('will_initialize','alloc llc')

    icount = 0
    
    DO iatom = 1, size(potparm,1)

      DO jatom = iatom, size(potparm,2)

        param(iatom,jatom) %lla => lla
        param(iatom,jatom) %llb => llb
        param(iatom,jatom) %llc => llc

        icount = icount + 1

! loop for a

        DO i = 1, numpar + 1
          a=potparm(iatom,jatom) %willis%a
          p(i,icount) = SQRT(a)
          lla%real_value => p(i,icount)
          ALLOCATE (lla%next,STAT=isos)
          IF (isos/=0) CALL stop_memory('will_initialize', 'alloc lla%next')
          lla => lla%next
        END DO

        icount = icount + 1

! loop for b

        DO i = 1, numpar + 1
          b=potparm(iatom,jatom) %willis%b
          p(i,icount) = SQRT(b)
          llb%real_value => p(i,icount)
          ALLOCATE (llb%next,STAT=isos)
          IF (isos/=0) CALL stop_memory('will_initialize', 'alloc llb%next')
          llb => llb%next
        END DO

        icount = icount + 1

! loop for c

        DO i = 1, numpar + 1
          c=potparm(iatom,jatom) %willis%c
          p(i,icount) = SQRT(c)
          llc%real_value => p(i,icount)
          ALLOCATE (llc%next,STAT=isos)
          IF (isos/=0) CALL stop_memory('will_initialize', 'alloc llc%next')
          llc => llc%next
        END DO

      END DO
      
    END DO

! scaling the appropriate interaction parameter by 1.1
! to construct the amoeba

    DO i = 0, numpar
      IF (i>0) p(i+1,i) = p(i+1,i)*1.1_dp
    END DO

  END SUBROUTINE will_initialize
!
!------------------------------------------------------------------------------!
    END MODULE fit_types
!-----------------------------------------------------------------------------!
