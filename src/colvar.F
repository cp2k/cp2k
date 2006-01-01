MODULE colvar
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/colvar [1.0] *
!!
!!   NAME
!!     colvar
!!
!!   FUNCTION
!!     collective variables for metadynamics
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Alessandro Laio,Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     04.2004 created
!!
!!   SOURCE
!****************************************************************************
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE metadynamics_types,              ONLY: colvar_type,&
                                             coord_colvar_id,&
                                             dist_colvar_id,&
                                             no_colvar_id,&
                                             torsion_colvar_id
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: massunit
  USE simulation_cell,                 ONLY: cell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_get
#include "cp_common_uses.h"
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='colvar'

  PUBLIC :: colvar_read
  PUBLIC :: colvar_eval_f

!***
!****************************************************************************
CONTAINS

!!****f* metadynamics/colvar_read *
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
  SUBROUTINE colvar_read(colvar, extended_lagrange, icol, colvar_section,  error)
    IMPLICIT NONE
    ! Arguments
    TYPE(colvar_type), INTENT(inout)         :: colvar
    LOGICAL                                  :: extended_lagrange
    TYPE(section_vals_type), POINTER, OPTIONAL :: colvar_section    
    INTEGER, INTENT(IN)                        :: icol
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Variables
    CHARACTER(len=*), PARAMETER :: routinen = 'colvar_read', &
         routinep = modulen//':'//routinen
    CHARACTER(len=3)                         :: test_result
    CHARACTER(len=default_string_length)     :: colvar_type_name, keyword
    INTEGER                                  :: i, ii, n_opt, n_opt_cv, &
                                                n_opt_var, ntemp, stat, n_rep
    LOGICAL                                  :: failure, RECOGNIZED_KEYWORD, &
                                                section_found
    LOGICAL, DIMENSION(4)                    :: my_subsection
    TYPE(section_vals_type), POINTER         :: distance_section, angle_section,&
                                                torsion_section, coordination_section    
    INTEGER, DIMENSION(:), POINTER           :: iatms
    

    NULLIFY(iatms)
    failure=.FALSE.
    n_opt=0
    n_opt_cv=0
    n_opt_var=0
    colvar%do_wall_p=.FALSE.
    colvar%do_wall_m=.FALSE.

    distance_section     => section_vals_get_subs_vals(colvar_section,"DISTANCE",i_rep_section=icol,error=error)
    angle_section        => section_vals_get_subs_vals(colvar_section,"ANGLE",i_rep_section=icol,error=error)
    torsion_section      => section_vals_get_subs_vals(colvar_section,"TORSION",i_rep_section=icol,error=error)
    coordination_section => section_vals_get_subs_vals(colvar_section,"COORDINATION",i_rep_section=icol,error=error)
    
    CALL section_vals_get(distance_section,explicit=my_subsection(1),error=error)
    CALL section_vals_get(angle_section,explicit=my_subsection(2),error=error)
    CALL section_vals_get(torsion_section,explicit=my_subsection(3),error=error)
    CALL section_vals_get(coordination_section,explicit=my_subsection(4),error=error)
    
    CPPostcondition(COUNT(my_subsection) == 1,cp_failure_level,routinep,error,failure)
    CALL section_vals_val_get(colvar_section,"SUBSYS",i_rep_section=icol,i_val=colvar%i_subsys,error=error)
    CALL section_vals_val_get(colvar_section,"MASS",i_rep_section=icol,r_val=colvar%mass,error=error)
    CALL section_vals_val_get(colvar_section,"SCALE",i_rep_section=icol,r_val=colvar%delta_s,error=error)
    CALL section_vals_val_get(colvar_section,"WALL_PLUS",i_rep_section=icol,n_rep_val=n_rep,error=error)
    IF (n_rep /= 0) THEN 
       colvar%do_wall_p=.TRUE.
       CALL section_vals_val_get(colvar_section,"WALL_PLUS",i_rep_section=icol,r_val=colvar%wall_p,error=error)
    END IF
    CALL section_vals_val_get(colvar_section,"WALL_MINUS",i_rep_section=icol,n_rep_val=n_rep,error=error)
    IF (n_rep /= 0) THEN 
       colvar%do_wall_m=.TRUE.
       CALL section_vals_val_get(colvar_section,"WALL_MINUS",i_rep_section=icol,r_val=colvar%wall_m,error=error)
    END IF
    CALL section_vals_val_get(colvar_section,"LAMBDA",i_rep_section=icol,r_val=colvar%lambda,error=error)
    
    IF      (my_subsection(1)) THEN
       ! Distance
       colvar%type_id=dist_colvar_id
       CALL section_vals_val_get(distance_section,"ATOMS",i_vals=iatms,error=error)
       colvar%dist_param%i_at = iatms(1)
       colvar%dist_param%j_at = iatms(2)
    ELSE IF (my_subsection(2)) THEN
       ! Angle
       WRITE(*,*)"Not implemented Yet!"
       STOP
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
    CASE(torsion_colvar_id)
       colvar%n_atom_s=4
       ALLOCATE(colvar%dsdr(3,colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(colvar%i_atom(colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DO ii=1,4
          colvar%i_atom(ii)=colvar%torsion_param%i_at_tors(ii)
       ENDDO
    END SELECT
    !
    ! Check Input Inconsistencies..
    !

    IF(extended_lagrange)THEN
       CALL section_vals_val_get(colvar_section,"MASS",i_rep_section=icol,  n_rep_val=n_rep, error=error)
       CPPostcondition(n_rep /= 0,cp_failure_level,routinep,error,failure)
       CALL section_vals_val_get(colvar_section,"SCALE",i_rep_section=icol, n_rep_val=n_rep, error=error)
       CPPostcondition(n_rep /= 0,cp_failure_level,routinep,error,failure)
       CALL section_vals_val_get(colvar_section,"LAMBDA",i_rep_section=icol, n_rep_val=n_rep, error=error)
       CPPostcondition(n_rep /= 0,cp_failure_level,routinep,error,failure)
    ENDIF

  END SUBROUTINE colvar_read


!!****f* colvar/colvar_eval_f [1.0] *
!!
!!   NAME
!!     colvar_eval_f
!!
!!   FUNCTION
!!     evaluates the force given and due to the given colvar
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
SUBROUTINE colvar_eval_f(colvar,force_env,error)
    TYPE(colvar_type), INTENT(inout)         :: colvar
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_eval_f', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys

  failure=.FALSE.
  NULLIFY(subsys,cell)
  CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
  IF (.NOT. failure) THEN
     SELECT CASE(colvar%type_id)
     CASE (dist_colvar_id)
        CALL dist_colvar(colvar,subsys,cell,error=error)
     CASE (coord_colvar_id)
        CALL coord_colvar(colvar,subsys,cell, error=error)
     CASE (torsion_colvar_id)
        CALL torsion_colvar(colvar,subsys,cell, error=error)
     CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     END SELECT
  END IF
END SUBROUTINE colvar_eval_f
!***************************************************************************

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
  SUBROUTINE dist_colvar(colvar,subsys,cell,error)
    TYPE(colvar_type)                        :: colvar
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dist_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure
    REAL(dp)                                 :: r12, ss(3), xij(3), xpi(3), &
                                                xpj(3)
    TYPE(particle_list_type), POINTER        :: particles_i

    failure=.FALSE.
    NULLIFY(particles_i)
    CPPrecondition(colvar%type_id==dist_colvar_id,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys(colvar%i_subsys)%subsys,&
         particles=particles_i,error=error)

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
  SUBROUTINE coord_colvar(colvar,subsys,cell,error)
    TYPE(colvar_type)                        :: colvar
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
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
    TYPE(particle_list_type), POINTER        :: particles_i

!

    failure=.FALSE.
    n_atoms_to=colvar%coord_param%n_atoms_to
    n_atoms_from=colvar%coord_param%n_atoms_from
    nncrd=colvar%coord_param%nncrd
    ndcrd=colvar%coord_param%ndcrd
    r_0=colvar%coord_param%r_0
    NULLIFY(particles_i)
    CPPrecondition(colvar%type_id==coord_colvar_id,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys(colvar%i_subsys)%subsys,&
         particles=particles_i,error=error)

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
  SUBROUTINE torsion_colvar(colvar,subsys,cell,error)
    TYPE(colvar_type)                        :: colvar
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
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
    TYPE(particle_list_type), POINTER        :: particles_i

! locals (tinker)
!

    failure=.FALSE.
    NULLIFY(particles_i)
    CPPrecondition(colvar%type_id==torsion_colvar_id,cp_failure_level,routineP,error,failure)
    CALL cp_subsys_get(subsys(colvar%i_subsys)%subsys,&
         particles=particles_i,error=error)

    DO ii=1,4
      i=colvar%torsion_param%i_at_tors(ii)
      rr(:,ii)=particles_i%els(i)%r(1:3)
    ENDDO


! temporary: from tinker improper torsion routines
            radian=57.29577951308232088_dp
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
                 IF(dt<180._dp) exit
                 dt = dt - 360.0e0_dp
               ENDDO
               DO ii=1,100
                 IF(dt>-180._dp) exit
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

END MODULE colvar
