!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002,2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****h* cp2k/colvar_methods [1.0] *
!!
!!   NAME
!!     colvar_methods
!!
!!   FUNCTION
!!     defines collective variables s({R}) and the derivative of this variable wrt R
!!     these can then be used in constraints, restraints and metadynamics ...
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Alessandro Laio,Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     04.2004 created
!!     01.2006 Refactored [Joost VandeVondele]
!!
!!   SOURCE
!****************************************************************************
MODULE colvar_methods
  USE cell_types,                      ONLY: cell_type
  USE global_types,                    ONLY: global_environment_type
  USE colvar_types,                    ONLY: colvar_type,&
                                             dist_colvar_id,&
                                             coord_colvar_id,&
                                             torsion_colvar_id,&
                                             plane_distance_colvar_id,&
                                             angle_colvar_id
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_list_types,             ONLY: particle_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE mathconstants,                   ONLY: pi
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colvar_methods'

  PUBLIC :: colvar_read, colvar_eval_f

!***
!****************************************************************************
CONTAINS
!!****f* colvar_methods/colvar_read *
!!
!!   NAME
!!     colvar_read
!!
!!   FUNCTION
!!     reads a colvar from the input
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - colvar: the place where to store what will be read
!!     - extended_lagrange: if the extended lagrangian formalismus should 
!!       be used
!!     - icol: number of the current colvar (repetition in colvar_section)
!!     - colvar_section: the colvar section
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!     04.2004 created [alessandro laio and fawzi mohamed]
!!
!!*** **********************************************************************
  SUBROUTINE colvar_read(colvar, icol, colvar_section, globenv, error)
    TYPE(colvar_type), INTENT(inout)         :: colvar
    INTEGER, INTENT(IN)                      :: icol
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: colvar_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ii, stat, iw, kk, iatm
    INTEGER, DIMENSION(:), POINTER           :: iatms
    LOGICAL                                  :: failure
    LOGICAL, DIMENSION(5)                    :: my_subsection
    TYPE(section_vals_type), POINTER         :: angle_section, &
                                                coordination_section, &
                                                distance_section, &
                                                torsion_section,&
                                                plane_dist_section

    CALL timeset(routineN,"I","",handle)

    NULLIFY(iatms)
    failure=.FALSE.

    distance_section     => section_vals_get_subs_vals(colvar_section,"DISTANCE",i_rep_section=icol,error=error)
    angle_section        => section_vals_get_subs_vals(colvar_section,"ANGLE",i_rep_section=icol,error=error)
    torsion_section      => section_vals_get_subs_vals(colvar_section,"TORSION",i_rep_section=icol,error=error)
    coordination_section => section_vals_get_subs_vals(colvar_section,"COORDINATION",i_rep_section=icol,error=error)
    plane_dist_section   => section_vals_get_subs_vals(colvar_section,"DISTANCE_POINT_PLANE",i_rep_section=icol,error=error)

    CALL section_vals_get(distance_section,explicit=my_subsection(1),error=error)
    CALL section_vals_get(angle_section,explicit=my_subsection(2),error=error)
    CALL section_vals_get(torsion_section,explicit=my_subsection(3),error=error)
    CALL section_vals_get(coordination_section,explicit=my_subsection(4),error=error)
    CALL section_vals_get(plane_dist_section,explicit=my_subsection(5),error=error)

    ! only one colvar can be present 
    CPPostcondition(COUNT(my_subsection) == 1,cp_failure_level,routinep,error,failure)

    IF      (my_subsection(1)) THEN
       ! Distance
       colvar%type_id=dist_colvar_id
       CALL section_vals_val_get(distance_section,"ATOMS",i_vals=iatms,error=error)
       colvar%dist_param%i_at = iatms(1)
       colvar%dist_param%j_at = iatms(2)
    ELSE IF (my_subsection(2)) THEN
       ! Angle
       colvar%type_id=angle_colvar_id
       CALL section_vals_val_get(angle_section,"ATOMS",i_vals=iatms,error=error)
       colvar%angle_param%i_at_angle = iatms
    ELSE IF (my_subsection(3)) THEN
       ! Torsion
       colvar%type_id=torsion_colvar_id
       CALL section_vals_val_get(torsion_section,"ATOMS",i_vals=iatms,error=error)
       colvar%torsion_param%i_at_tors = iatms
    ELSE IF (my_subsection(4)) THEN
       ! Coordination
       colvar%type_id=coord_colvar_id
       CALL section_vals_val_get(coordination_section,"ATOMS_FROM",i_vals=iatms,error=error)
       colvar%coord_param%n_atoms_from = SIZE(iatms)
       ALLOCATE(colvar%coord_param%i_at_from(SIZE(iatms)),stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       colvar%coord_param%i_at_from = iatms
       CALL section_vals_val_get(coordination_section,"ATOMS_TO",i_vals=iatms,error=error)
       colvar%coord_param%n_atoms_to = SIZE(iatms)
       ALLOCATE(colvar%coord_param%i_at_to(SIZE(iatms)),stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       colvar%coord_param%i_at_to = iatms       
       CALL section_vals_val_get(coordination_section,"R0",r_val=colvar%coord_param%r_0,error=error)
       CALL section_vals_val_get(coordination_section,"NN",i_val=colvar%coord_param%nncrd,error=error)
       CALL section_vals_val_get(coordination_section,"ND",i_val=colvar%coord_param%ndcrd,error=error)
    ELSE IF (my_subsection(5)) THEN
       ! Distance point from plane
       colvar%type_id=plane_distance_colvar_id
       CALL section_vals_val_get(plane_dist_section,"ATOMS_PLANE",i_vals=iatms,error=error)
       CPPostcondition(SIZE(iatms) == 3,cp_failure_level,routinep,error,failure)
       colvar%plane_distance_param%plane = iatms
       CALL section_vals_val_get(plane_dist_section,"ATOM_POINT",i_val=iatm,error=error)
       colvar%plane_distance_param%point = iatm
    END IF
    
    SELECT CASE(colvar%type_id)
    CASE(dist_colvar_id)
       colvar%n_atom_s=2
       ALLOCATE(colvar%dsdr(3,colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(colvar%i_atom(colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       colvar%i_atom(1)=colvar%dist_param%i_at
       colvar%i_atom(2)=colvar%dist_param%j_at
    CASE(coord_colvar_id)
       colvar%n_atom_s=colvar%coord_param%n_atoms_from+colvar%coord_param%n_atoms_to
       ALLOCATE(colvar%dsdr(3,colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(colvar%i_atom(colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DO ii=1,colvar%coord_param%n_atoms_from
          colvar%i_atom(ii)=colvar%coord_param%i_at_from(ii)
       ENDDO
       DO ii=1,colvar%coord_param%n_atoms_to
          colvar%i_atom(colvar%coord_param%n_atoms_from+ii)=colvar%coord_param%i_at_to(ii)
       ENDDO
    CASE(angle_colvar_id)
       colvar%n_atom_s=3
       ALLOCATE(colvar%dsdr(3,colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(colvar%i_atom(colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DO ii=1,3
          colvar%i_atom(ii)=colvar%angle_param%i_at_angle(ii)
       ENDDO
    CASE(torsion_colvar_id)
       colvar%n_atom_s=4
       ALLOCATE(colvar%dsdr(3,colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(colvar%i_atom(colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DO ii=1,4
          colvar%i_atom(ii)=colvar%torsion_param%i_at_tors(ii)
       ENDDO
    CASE(plane_distance_colvar_id)
       colvar%n_atom_s=4
       ALLOCATE(colvar%dsdr(3,colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(colvar%i_atom(colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DO ii=1,3
          colvar%i_atom(ii)=colvar%plane_distance_param%plane(ii)
       ENDDO
       colvar%i_atom(4)=colvar%plane_distance_param%point
    END SELECT

    IF (globenv % ionode .AND. globenv % print_level>=0) THEN
       iw = globenv % scr
       WRITE ( iw, '( A )')'          ----------------------------------------------------------------------'
       SELECT CASE(colvar%type_id)
       CASE(angle_colvar_id)
          WRITE ( iw, '( A,T71,3I5)' )   ' COLVARS| ANGLE          >>> ATOMS:',colvar%angle_param%i_at_angle
       CASE(plane_distance_colvar_id)
          WRITE ( iw, '( A,T61,3I5)' )   ' COLVARS| PLANE DISTANCE    - PLANE  >>> ATOMS:',colvar%plane_distance_param%plane
          WRITE ( iw, '( A,T61,1I5)' )   ' COLVARS| PLANE DISTANCE    - POINT  >>>  ATOM:',colvar%plane_distance_param%point
       CASE(torsion_colvar_id)
          WRITE ( iw, '( A,T61,4I5)' )   ' COLVARS| TORSION       >>> ATOMS:',colvar%torsion_param%i_at_tors
       CASE(dist_colvar_id)
          WRITE ( iw, '( A,T71,2I5)' )   ' COLVARS| BOND          >>> ATOMS:',colvar%dist_param%i_at,colvar%dist_param%j_at
       CASE(coord_colvar_id)
          WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COORDINATION  >>> FROM ATOM:', colvar%coord_param%i_at_from(kk),&
               kk=1,SIZE(colvar%coord_param%i_at_from))
          WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COORDINATION  >>>   TO ATOM:', colvar%coord_param%i_at_to(kk),&
               kk=1,SIZE(colvar%coord_param%i_at_to))
          WRITE ( iw, '( A,T71,F10.5)' ) ' COLVARS| R0',colvar%coord_param%r_0
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| NN',colvar%coord_param%nncrd
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| ND',colvar%coord_param%ndcrd
       END SELECT
       WRITE ( iw, '( A )')'          ----------------------------------------------------------------------'
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE colvar_read

!!****f* colvar/colvar_eval_f [1.0] *
!!
!!   NAME
!!     colvar_eval_f
!!
!!   FUNCTION
!!     evaluates the force (dsdr) given and due to the given colvar
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - colvar: the collective variable to evaluate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Alessandro Laio and fawzi
!!
!!*** **********************************************************************
SUBROUTINE colvar_eval_f(isubsys,icolvar,force_env,error)
    INTEGER                                  :: isubsys, icolvar
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_eval_f', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: colvar_ok, failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys

  failure=.FALSE.
  NULLIFY(subsys,cell)
  CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
  colvar_ok=ASSOCIATED(subsys(isubsys)%subsys%colvar)
  CPAssert(colvar_ok,cp_failure_level,routineP,error,failure)
  
  IF (.NOT. failure) THEN
     SELECT CASE(subsys(isubsys)%subsys%colvar(icolvar)%type_id)
     CASE (dist_colvar_id)
        CALL dist_colvar(icolvar,subsys(isubsys)%subsys,cell,error=error)
     CASE (coord_colvar_id)
        CALL coord_colvar(icolvar,subsys(isubsys)%subsys,cell, error=error)
     CASE (torsion_colvar_id)
        CALL torsion_colvar(icolvar,subsys(isubsys)%subsys,cell, error=error)
     CASE (angle_colvar_id)
        CALL angle_colvar(icolvar,subsys(isubsys)%subsys,cell, error=error)
     CASE (plane_distance_colvar_id)
        CALL plane_distance_colvar(icolvar,subsys(isubsys)%subsys,cell, error=error)
     CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     END SELECT
  END IF
END SUBROUTINE colvar_eval_f
!***************************************************************************

!!****f* colvar/plane_distance_colvar [1.0] *
!!
!!   NAME
!!     plane_distance_colvar
!!
!!   FUNCTION
!!     evaluates the force due (and on) the distance from the plane collective variable
!!
!!   NOTES
!!     -
!!   AUTHOR
!!     Teodoro Laino 02.2006 [created]
!!
!!*** **********************************************************************
  SUBROUTINE plane_distance_colvar(icolvar,subsys,cell,error)
    INTEGER                                  :: icolvar
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'plane_distance_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, l
    LOGICAL                                  :: failure
    REAL(dp)                                 :: xpij(3), a, b, &
                                                xpkj(3),&
                                                dxpndxi(3,3),&
                                                dxpndxj(3,3),&
                                                dxpndxk(3,3),&
                                                xpn(3), xpl(3), r12, dsdxpn(3)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(particle_list_type), POINTER        :: particles_i

    failure=.FALSE.
    NULLIFY(particles_i)
    colvar=>subsys%colvar(icolvar)

    CPPrecondition(colvar%type_id==plane_distance_colvar_id,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys,particles=particles_i,error=error)

    i=colvar%plane_distance_param%plane(1)
    j=colvar%plane_distance_param%plane(2)
    k=colvar%plane_distance_param%plane(3)
    l=colvar%plane_distance_param%point
    xpij=particles_i%els(i)%r(1:3)-particles_i%els(j)%r(1:3)
    xpkj=particles_i%els(k)%r(1:3)-particles_i%els(j)%r(1:3)
    xpl =particles_i%els(l)%r(1:3)-(particles_i%els(i)%r(1:3)+&
                                    particles_i%els(j)%r(1:3)+&
                                    particles_i%els(k)%r(1:3))/3.0_dp
    xpn(1)  = xpij(2)*xpkj(3)-xpij(3)*xpkj(2)
    xpn(2)  = xpij(3)*xpkj(1)-xpij(1)*xpkj(3)
    xpn(3)  = xpij(1)*xpkj(2)-xpij(2)*xpkj(1)
    a = DOT_PRODUCT(xpn,xpn)
    b = DOT_PRODUCT(xpl,xpn)
    r12=SQRT(a)
    colvar%ss=b/r12
    dsdxpn(1) = xpl(1)/r12 - b*xpn(1)/(r12*a)
    dsdxpn(2) = xpl(2)/r12 - b*xpn(2)/(r12*a)
    dsdxpn(3) = xpl(3)/r12 - b*xpn(3)/(r12*a)
    !
    dxpndxi(1,1)=  0.0_dp
    dxpndxi(1,2)=  1.0_dp * xpkj(3)
    dxpndxi(1,3)= -1.0_dp * xpkj(2)
    dxpndxi(2,1)= -1.0_dp * xpkj(3)
    dxpndxi(2,2)=  0.0_dp
    dxpndxi(2,3)=  1.0_dp * xpkj(1)
    dxpndxi(3,1)=  1.0_dp * xpkj(2)
    dxpndxi(3,2)= -1.0_dp * xpkj(1)
    dxpndxi(3,3)=  0.0_dp
    !
    dxpndxj(1,1)=  0.0_dp
    dxpndxj(1,2)= -1.0_dp * xpkj(3) + xpij(3)
    dxpndxj(1,3)= -1.0_dp * xpij(2) + xpkj(2)
    dxpndxj(2,1)= -1.0_dp * xpij(3) + xpkj(3)
    dxpndxj(2,2)=  0.0_dp
    dxpndxj(2,3)= -1.0_dp * xpkj(1) + xpij(1)
    dxpndxj(3,1)= -1.0_dp * xpkj(2) + xpij(2)
    dxpndxj(3,2)= -1.0_dp * xpij(1) + xpkj(1)
    dxpndxj(3,3)=  0.0_dp
    !
    dxpndxk(1,1)=  0.0_dp
    dxpndxk(1,2)= -1.0_dp * xpij(3)
    dxpndxk(1,3)=  1.0_dp * xpij(2)
    dxpndxk(2,1)=  1.0_dp * xpij(3)
    dxpndxk(2,2)=  0.0_dp
    dxpndxk(2,3)= -1.0_dp * xpij(1)
    dxpndxk(3,1)= -1.0_dp * xpij(2)
    dxpndxk(3,2)=  1.0_dp * xpij(1)
    dxpndxk(3,3)=  0.0_dp
    !
    colvar%dsdr(:,1)=MATMUL(dsdxpn,dxpndxi)-xpn/(3.0_dp*r12)
    colvar%dsdr(:,2)=MATMUL(dsdxpn,dxpndxj)-xpn/(3.0_dp*r12)
    colvar%dsdr(:,3)=MATMUL(dsdxpn,dxpndxk)-xpn/(3.0_dp*r12)
    colvar%dsdr(:,4)=xpn/r12

  END SUBROUTINE plane_distance_colvar

!!****f* colvar/angle_colvar [1.0] *
!!
!!   NAME
!!     angle_colvar
!!
!!   FUNCTION
!!     evaluates the force due (and on) the distance from the plane collective variable
!!
!!   NOTES
!!     -
!!   AUTHOR
!!     Teodoro Laino 02.2006 [created]
!!
!!*** **********************************************************************
  SUBROUTINE angle_colvar(icolvar,subsys,cell,error)
    INTEGER                                  :: icolvar
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'angle_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k
    LOGICAL                                  :: failure
    REAL(dp)                                 :: xij(3), xkj(3), a, b, t0, t1, t2, t3, fmod
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(particle_list_type), POINTER        :: particles_i

    failure=.FALSE.
    NULLIFY(particles_i)
    colvar=>subsys%colvar(icolvar)

    CPPrecondition(colvar%type_id==angle_colvar_id,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys,particles=particles_i,error=error)

    i=colvar%angle_param%i_at_angle(1)
    j=colvar%angle_param%i_at_angle(2)
    k=colvar%angle_param%i_at_angle(3)
    xij=particles_i%els(i)%r(1:3)-particles_i%els(j)%r(1:3)
    xkj=particles_i%els(k)%r(1:3)-particles_i%els(j)%r(1:3)
    a = SQRT(DOT_PRODUCT(xij,xij))
    b = SQRT(DOT_PRODUCT(xkj,xkj))
    t0 = 1.0_dp/(a*b)
    t1 = 1.0_dp/(a**3.0_dp*b)
    t2 = 1.0_dp/(a*b**3.0_dp)
    t3 = DOT_PRODUCT(xij,xkj)
    colvar%ss = ACOS(t3*t0)
    IF ((ABS(colvar%ss).LT.0.00001D0).OR.(ABS(colvar%ss-pi).LT.0.00001D0)) THEN
       fmod =   0.D0
    ELSE
       fmod = - 1.D0 / SIN(colvar%ss)
    ENDIF
    colvar%dsdr(:,1)= xkj(:)*t0-xij(:)*t1*t3
    colvar%dsdr(:,2)=-xkj(:)*t0+xij(:)*t1*t3-xij(:)*t0+xkj(:)*t2*t3
    colvar%dsdr(:,3)= xij(:)*t0-xkj(:)*t2*t3
    colvar%dsdr = colvar%dsdr * fmod
    
  END SUBROUTINE angle_colvar


!!****f* colvar/dist_colvar [1.0] *
!!
!!   NAME
!!     dist_colvar
!!
!!   FUNCTION
!!     evaluates the force due (and on) the distance collective variable
!!
!!   NOTES
!!     -
!!   AUTHOR
!!     Alessandro Laio, Fawzi Mohamed
!!
!!*** **********************************************************************
  SUBROUTINE dist_colvar(icolvar,subsys,cell,error)
    INTEGER                                  :: icolvar
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dist_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure
    REAL(dp)                                 :: r12, ss(3), xij(3), xpi(3), &
                                                xpj(3)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(particle_list_type), POINTER        :: particles_i

    failure=.FALSE.
    NULLIFY(particles_i)
    colvar=>subsys%colvar(icolvar)

    CPPrecondition(colvar%type_id==dist_colvar_id,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys,particles=particles_i,error=error)

    i=colvar%dist_param%i_at
    j=colvar%dist_param%j_at
    xpi=particles_i%els(i)%r(1:3)
    xpj=particles_i%els(j)%r(1:3)
    ss=MATMUL(cell%h_inv,xpi-xpj)
    ss=ss-NINT(ss)
    xij=MATMUL(cell%hmat,ss)
    r12=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
    colvar%ss=r12
    colvar%dsdr(:,1)=xij/r12
    colvar%dsdr(:,2)=-xij/r12

  END SUBROUTINE dist_colvar
!***************************************************************************


!!****f* colvar/coord_colvar [1.0] *
!!
!!   NAME
!!     coord_colvar
!!
!!   FUNCTION
!!     evaluates the force due (and on) the coordination collective variable
!!
!!   NOTES
!!
!!   AUTHOR
!!     Alessandro Laio, Fawzi Mohamed
!!
!!*** **********************************************************************
  SUBROUTINE coord_colvar(icolvar,subsys,cell,error)
    INTEGER                                  :: icolvar
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coord_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, jj, n_atoms_from, &
                                                n_atoms_to, ndcrd, nncrd
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dfunc, func, invden, ncoord, &
                                                num, r12, r_0, rdist, ss(3), &
                                                xij(3)
    REAL(dp), DIMENSION(3)                   :: xpi, xpj
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(particle_list_type), POINTER        :: particles_i

!

    failure=.FALSE.
    colvar=>subsys%colvar(icolvar)
    n_atoms_to=colvar%coord_param%n_atoms_to
    n_atoms_from=colvar%coord_param%n_atoms_from
    nncrd=colvar%coord_param%nncrd
    ndcrd=colvar%coord_param%ndcrd
    r_0=colvar%coord_param%r_0
    NULLIFY(particles_i)
    CPPrecondition(colvar%type_id==coord_colvar_id,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys,particles=particles_i,error=error)

    ncoord=0.0_dp
    colvar%dsdr=0.0_dp
    DO ii=1,n_atoms_from
      i=colvar%coord_param%i_at_from(ii)
      xpi(:)=particles_i%els(i)%r(1:3)
      DO jj=1,n_atoms_to
        i=colvar%coord_param%i_at_to(jj)
        xpj(:)=particles_i%els(i)%r(1:3)
        ss=MATMUL(cell%h_inv,xpi(:)-xpj(:))
        ss=ss-NINT(ss)
        xij=MATMUL(cell%hmat,ss)
        r12=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
        IF(r12 < 1.0e-8_dp)CYCLE
        rdist = r12/r_0
        num=(1.0_dp-rdist**nncrd)
        invden=1.0_dp/(1.0_dp-rdist**ndcrd)
        func=num*invden
        dfunc= (- nncrd*rdist **(nncrd-1) *invden  &
          + num*(invden)**2 * ndcrd *rdist **(ndcrd-1))/(r12*r_0)
                                                                               
        ncoord=ncoord+func
        colvar%dsdr(:,ii)=colvar%dsdr(:,ii)+dfunc*xij
        colvar%dsdr(:,n_atoms_from+jj)=colvar%dsdr(:,n_atoms_from+jj)-dfunc*xij
      ENDDO
    ENDDO
    colvar%ss=ncoord

  END SUBROUTINE coord_colvar
!***************************************************************************

!!****f* colvar/torsion_colvar [1.0] *
!!
!!   NAME
!!     torsion_colvar
!!
!!   FUNCTION
!!     evaluates the force due to the torsion collective variable
!!
!!   NOTES
!!
!!   AUTHOR
!!     Alessandro Laio, Fawzi Mohamed
!!
!!*** **********************************************************************
  SUBROUTINE torsion_colvar(icolvar,subsys,cell,error)
    INTEGER                                  :: icolvar
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'torsion_colvar', &
         routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii
    LOGICAL                                  :: failure
    REAL(dp) :: angle, cosine, dedphi, dedxia, dedxib, dedxic, dedxid, dedxt, &
         dedxu, dedyia, dedyib, dedyic, dedyid, dedyt, dedyu, dedzia, dedzib, &
         dedzic, dedzid, dedzt, dedzu, dt, e, radian, rcb, rt2, rtru, ru2, sine, &
         xba, xca, xcb, xdb, xdc, xia, xib, xic, xid, xt, xtu, xu, yba, yca, &
         ycb, ydb, ydc, yia, yib, yic, yid, yt, ytu, yu, zba, zca, zcb, zdb, &
         zdc, zia, zib, zic, zid, zt, ztu, zu
    REAL(dp), DIMENSION(3, 4)                :: rr
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(particle_list_type), POINTER        :: particles_i

    ! locals (tinker)
    !

    failure=.FALSE.
    colvar=>subsys%colvar(icolvar)
    NULLIFY(particles_i)
    CPPrecondition(colvar%type_id==torsion_colvar_id,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys,particles=particles_i,error=error)

    DO ii=1,4
       i=colvar%torsion_param%i_at_tors(ii)
       rr(:,ii)=particles_i%els(i)%r(1:3)
    ENDDO

    ! temporary: from tinker improper torsion routines
    radian=180.0_dp/pi
    xia = rr(1,1)
    yia = rr(2,1)
    zia = rr(3,1)
    xib = rr(1,2)
    yib = rr(2,2)
    zib = rr(3,2)
    xic = rr(1,3)
    yic = rr(2,3)
    zic = rr(3,3)
    xid = rr(1,4)
    yid = rr(2,4)
    zid = rr(3,4)
    xba = xib - xia
    yba = yib - yia
    zba = zib - zia
    xcb = xic - xib
    ycb = yic - yib
    zcb = zic - zib
    xdc = xid - xic
    ydc = yid - yic
    zdc = zid - zic
    xt = yba*zcb - ycb*zba
    yt = zba*xcb - zcb*xba
    zt = xba*ycb - xcb*yba
    xu = ycb*zdc - ydc*zcb
    yu = zcb*xdc - zdc*xcb
    zu = xcb*ydc - xdc*ycb
    xtu = yt*zu - yu*zt
    ytu = zt*xu - zu*xt
    ztu = xt*yu - xu*yt
    rt2 = xt*xt + yt*yt + zt*zt
    ru2 = xu*xu + yu*yu + zu*zu
    rtru = SQRT(rt2 * ru2)
    IF (rtru .NE. 0.0e0_dp) THEN
       rcb = SQRT(xcb*xcb + ycb*ycb + zcb*zcb)
       cosine = (xt*xu + yt*yu + zt*zu) / rtru
       sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
       cosine = MIN(1.0e0_dp,MAX(-1.0e0_dp,cosine))
       angle = radian * ACOS(cosine)
       IF (sine .LT. 0.0e0_dp)  angle = -angle
       dt = angle 
       DO ii=1,100
          IF(dt<180._dp) EXIT
          dt = dt - 360.0e0_dp
       ENDDO
       DO ii=1,100
          IF(dt>-180._dp) EXIT
          dt = dt + 360.0e0_dp
       ENDDO
       dt = dt / radian
       !
       !     calculate improper energy and master chain rule term
       !
       e = dt
       dedphi = 1.0_dp
       !
       !     chain rule terms for first derivative components
       !
       xca = xic - xia
       yca = yic - yia
       zca = zic - zia
       xdb = xid - xib
       ydb = yid - yib
       zdb = zid - zib
       dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
       dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
       dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
       dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb)
       dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb)
       dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb)
       !
       !     compute first derivative components for this angle
       !
       dedxia = zcb*dedyt - ycb*dedzt
       dedyia = xcb*dedzt - zcb*dedxt
       dedzia = ycb*dedxt - xcb*dedyt
       dedzia = ycb*dedxt - xcb*dedyt
       dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu
       dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu
       dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu
       dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu
       dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu
       dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu
       dedxid = zcb*dedyu - ycb*dedzu
       dedyid = xcb*dedzu - zcb*dedxu
       dedzid = ycb*dedxu - xcb*dedyu
    ENDIF
    !
    colvar%ss=e
    colvar%dsdr(1,1)=dedxia
    colvar%dsdr(2,1)=dedyia
    colvar%dsdr(3,1)=dedzia
    colvar%dsdr(1,2)=dedxib
    colvar%dsdr(2,2)=dedyib
    colvar%dsdr(3,2)=dedzib
    colvar%dsdr(1,3)=dedxic
    colvar%dsdr(2,3)=dedyic
    colvar%dsdr(3,3)=dedzic
    colvar%dsdr(1,4)=dedxid
    colvar%dsdr(2,4)=dedyid
    colvar%dsdr(3,4)=dedzid

  END SUBROUTINE torsion_colvar
!***************************************************************************

END MODULE colvar_methods
