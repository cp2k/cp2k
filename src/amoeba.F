!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/amoeba [1.0] *
!!
!!   NAME
!!     amoeba, NUMERICAL RECIPIES SUBROUTINE AMOEBA ADAPTED FOR USE IN CP2K
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

 MODULE amoeba

  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE fit_parameters_types,            ONLY: fit_info_type,&
                                             fit_parameter_type,&
                                             linklist_parameters_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type,&
                                             thermodynamic_type
  USE nrutil,                          ONLY: assert_eq,&
                                             imaxloc,&
                                             iminloc,&
                                             swap
  USE particle_types,                  ONLY: particle_prop_type
  USE simulation_cell,                 ONLY: cell_type
  USE structure_types,                 ONLY: interaction_type,&
                                             structure_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: amoeba_evaluate

!!*****

CONTAINS

!!****** amoeba/amoeba_evaluate [1.0] *
!!
!!   NAME
!!     amoeba_evaluate
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
!-----------------------------------------------------------------------------!
!   NUMERICAL RECIPIES SUBROUTINE AMOEBA ADAPTED FOR USE IN CP2K              !
!-----------------------------------------------------------------------------!

  SUBROUTINE amoeba_evaluate ( p, y, rtol, struc, thermo, inter, simpar, &
           fit_param, box_ref, finfo, param, pstat, ewald_param, globenv, func )


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: p
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)   :: y
    REAL(KIND=dp), INTENT(OUT)                   :: rtol
    TYPE(structure_type), INTENT(INOUT)      :: struc
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(interaction_type), INTENT(INOUT)    :: inter
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
  INTERFACE
     FUNCTION func ( x ,struc,inter,thermo,simpar,fit_param,box_ref,&
       finfo,param,pstat,ewald_param, globenv)
       
      USE kinds, ONLY: dp
      USE structure_types, ONLY : structure_type, interaction_type
      USE particle_types, ONLY : particle_prop_type
      USE md, ONLY : simulation_parameters_type, thermodynamic_type
      USE fit_parameters_types, ONLY : fit_parameter_type, linklist_parameters_type, &
                                       fit_info_type
      USE ewald_parameters_types, ONLY : ewald_parameters_type
      USE global_types, ONLY : global_environment_type
      USE simulation_cell, ONLY : cell_type

  
      IMPLICIT NONE
      TYPE (structure_type), INTENT (INOUT) :: struc
      TYPE (interaction_type), INTENT (INOUT) :: inter
      TYPE (thermodynamic_type), INTENT (INOUT) :: thermo
      TYPE (simulation_parameters_type), INTENT (IN) :: simpar
      REAL(KIND=dp), INTENT (IN), DIMENSION (:) :: x
      TYPE (fit_parameter_type), INTENT(INOUT) :: fit_param
      TYPE (cell_type) , POINTER, DIMENSION (:) :: box_ref
      TYPE (fit_info_type), DIMENSION (:,:), POINTER :: finfo
      TYPE (linklist_parameters_type ), INTENT(INOUT), DIMENSION ( :, : ) :: param 
      TYPE (particle_prop_type ), DIMENSION ( : ), POINTER :: pstat 
      TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param  
      TYPE ( global_environment_type), INTENT ( IN ) :: globenv
      REAL (KIND=dp) :: func
     END FUNCTION func

  END INTERFACE

    INTEGER                                  :: ihi, itmax, ndim
    REAL(KIND=dp)                                :: ftol
    REAL(KIND=dp), DIMENSION(size(p, 2))         :: psum

  ftol = fit_param % max_tol
  itmax = fit_param % iwrite

  CALL amoeba_private()

!------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE amoeba_private
    IMPLICIT NONE
    INTEGER :: i, ilo, inhi
    INTEGER :: iter
    REAL(KIND=dp) :: ysave, ytry, ytmp

    ndim = assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
    iter = 0
    psum(:) = sum(p(:,:),dim=1)
    DO
       ilo = iminloc(y(:))
       ihi = imaxloc(y(:))
       ytmp = y(ihi)
       y(ihi) = y(ilo)
       inhi = imaxloc(y(:))
       y(ihi) = ytmp
       rtol = 2.0_dp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
       IF (rtol<ftol) THEN
          CALL swap(y(1),y(ilo))
          CALL swap(p(1,:),p(ilo,:))
          RETURN
       END IF
       IF (iter>=itmax) RETURN
       ytry = amotry(-1.0_dp)
       iter = iter + 1
       IF (ytry<=y(ilo)) THEN
          ytry = amotry(2.0_dp)
          iter = iter + 1
       ELSE IF (ytry>=y(inhi)) THEN
          ysave = y(ihi)
          ytry = amotry(0.5_dp)
          iter = iter + 1
          IF (ytry>=ysave) THEN
             p(:,:) = 0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
             DO i = 1, ndim + 1
                IF (i/=ilo) y(i) = func(p(i,:),struc,inter,thermo,simpar,fit_param,&
                 box_ref,finfo,param, pstat, ewald_param, globenv)
             END DO
             iter = iter + ndim
             psum(:) = sum(p(:,:),dim=1)
          END IF
       END IF
    END DO
  END SUBROUTINE amoeba_private

!BL
  FUNCTION amotry(fac)
    REAL(KIND=dp), INTENT(IN)                    :: fac
    REAL(KIND=dp)                                :: amotry

    REAL(KIND=dp)                                :: fac1, fac2, ytry
    REAL(KIND=dp), DIMENSION(size(p, 2))         :: ptry

    fac1 = (1.0_dp-fac)/ndim
    fac2 = fac1 - fac
    ptry(:) = psum(:)*fac1 - p(ihi,:)*fac2
    ytry = func(ptry,struc,inter,thermo,simpar,fit_param,&
               box_ref,finfo,param, pstat, ewald_param, globenv)
    IF (ytry<y(ihi)) THEN
       y(ihi) = ytry
       psum(:) = psum(:) - p(ihi,:) + ptry(:)
       p(ihi,:) = ptry(:)
    END IF
    amotry = ytry
  END FUNCTION amotry

END SUBROUTINE amoeba_evaluate

!!*****
!******************************************************************************

END MODULE amoeba
