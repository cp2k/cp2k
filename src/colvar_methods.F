!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
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
  USE colvar_types,                    ONLY: &
       angle_colvar_id, colvar_check_points, colvar_create, colvar_setup, &
       colvar_type, coord_colvar_id, ddiff_colvar_id, dist_colvar_id, &
       eval_point_der, eval_point_pos, hydronium_colvar_id, &
       plane_distance_colvar_id, qparm_colvar_id, rotation_colvar_id, &
       torsion_colvar_id
  USE constraint_fxd,                  ONLY: check_fixed_atom_cns_colv
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: fac,&
                                             maxfac,&
                                             pi
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: fixd_constraint_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE spherical_harmonics,             ONLY: dlegendre,&
                                             legendre
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colvar_methods'

  PUBLIC :: colvar_read,&
            colvar_eval_glob_f,&
            colvar_eval_mol_f

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
  SUBROUTINE colvar_read(colvar, icol, colvar_section, error)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER, INTENT(IN)                      :: icol
    TYPE(section_vals_type), POINTER         :: colvar_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_read', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=7)                         :: tag
    CHARACTER(LEN=default_string_length)     :: tmpStr, tmpStr2
    INTEGER                                  :: handle, iatm, iw, k, kk, &
                                                n_var, ndim, stat
    INTEGER, DIMENSION(:), POINTER           :: iatms
    LOGICAL                                  :: failure
    LOGICAL, DIMENSION(9)                    :: my_subsection
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER :: angle_section, coordination_section, &
      ddiff_section, distance_section, hydronium_section, plane_dist_section, &
      point_section, qparm_section, rotation_section, torsion_section, &
      wrk_section

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,"I","",handle)

    NULLIFY(iatms)
    failure=.FALSE.

    distance_section     => section_vals_get_subs_vals(colvar_section,"DISTANCE",i_rep_section=icol,error=error)
    ddiff_section        => section_vals_get_subs_vals(colvar_section,"DISTANCE_DIFFERENCE",&
         i_rep_section=icol,error=error)
    angle_section        => section_vals_get_subs_vals(colvar_section,"ANGLE",i_rep_section=icol,error=error)
    torsion_section      => section_vals_get_subs_vals(colvar_section,"TORSION",i_rep_section=icol,error=error)
    coordination_section => section_vals_get_subs_vals(colvar_section,"COORDINATION",i_rep_section=icol,error=error)
    plane_dist_section   => section_vals_get_subs_vals(colvar_section,"DISTANCE_POINT_PLANE",i_rep_section=icol,&
                            error=error)
    rotation_section     => section_vals_get_subs_vals(colvar_section,"BOND_ROTATION",i_rep_section=icol,error=error)
    qparm_section        => section_vals_get_subs_vals(colvar_section,"QPARM",i_rep_section=icol,error=error)
    hydronium_section    => section_vals_get_subs_vals(colvar_section,"HYDRONIUM",i_rep_section=icol,error=error)

    CALL section_vals_get(distance_section,     explicit=my_subsection(1), error=error)
    CALL section_vals_get(angle_section,        explicit=my_subsection(2), error=error)
    CALL section_vals_get(torsion_section,      explicit=my_subsection(3), error=error)
    CALL section_vals_get(coordination_section, explicit=my_subsection(4), error=error)
    CALL section_vals_get(plane_dist_section,   explicit=my_subsection(5), error=error)
    CALL section_vals_get(rotation_section,     explicit=my_subsection(6), error=error)
    CALL section_vals_get(ddiff_section,        explicit=my_subsection(7), error=error)
    CALL section_vals_get(qparm_section,        explicit=my_subsection(8), error=error)
    CALL section_vals_get(hydronium_section,    explicit=my_subsection(9),error=error)
    ! Only one colvar can be present
    CPPostcondition(COUNT(my_subsection) == 1,cp_failure_level,routinep,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(colvar),cp_failure_level,routinep,error,failure)

    IF      (my_subsection(1)) THEN
       ! Distance
       wrk_section => distance_section
       CALL colvar_create(colvar, dist_colvar_id, error)
       CALL colvar_check_points(colvar, distance_section, error)
       CALL section_vals_val_get(distance_section,"ATOMS",i_vals=iatms,error=error)
       colvar%dist_param%i_at = iatms(1)
       colvar%dist_param%j_at = iatms(2)
    ELSE IF (my_subsection(2)) THEN
       ! Angle
       wrk_section => angle_section
       CALL colvar_create(colvar, angle_colvar_id, error)
       CALL colvar_check_points(colvar, angle_section, error)
       CALL section_vals_val_get(angle_section,"ATOMS",i_vals=iatms,error=error)
       colvar%angle_param%i_at_angle = iatms
    ELSE IF (my_subsection(3)) THEN
       ! Torsion
       wrk_section => torsion_section
       CALL colvar_create(colvar, torsion_colvar_id, error)
       CALL colvar_check_points(colvar, torsion_section, error)
       CALL section_vals_val_get(torsion_section,"ATOMS",i_vals=iatms,error=error)
       colvar%torsion_param%i_at_tors = iatms
       colvar%torsion_param%o0        = 0.0_dp
    ELSE IF (my_subsection(4)) THEN
       ! Coordination
       wrk_section => coordination_section
       CALL colvar_create(colvar, coord_colvar_id, error)
       CALL colvar_check_points(colvar, coordination_section, error)
       ! This section can be repeated
       CALL section_vals_val_get(coordination_section,"ATOMS_FROM",n_rep_val=n_var,error=error)
       ndim = 0
       NULLIFY(colvar%coord_param%i_at_from)
       NULLIFY(colvar%coord_param%i_at_to)
       DO k = 1, n_var
          CALL section_vals_val_get(coordination_section,"ATOMS_FROM",i_rep_val=k,i_vals=iatms,error=error)
          CALL reallocate(colvar%coord_param%i_at_from,1, ndim+SIZE(iatms))
          colvar%coord_param%i_at_from(ndim+1:ndim+SIZE(iatms)) = iatms 
          ndim = ndim + SIZE(iatms)
       END DO
       colvar%coord_param%n_atoms_from = ndim
       ! This section can be repeated
       CALL section_vals_val_get(coordination_section,"ATOMS_TO",n_rep_val=n_var,error=error)
       ndim = 0
       DO k = 1, n_var
          CALL section_vals_val_get(coordination_section,"ATOMS_TO",i_rep_val=k,i_vals=iatms,error=error)
          CALL reallocate(colvar%coord_param%i_at_to,1, ndim+SIZE(iatms))
          colvar%coord_param%i_at_to(ndim+1:ndim+SIZE(iatms)) = iatms 
          ndim = ndim + SIZE(iatms)
       END DO
       colvar%coord_param%n_atoms_to = ndim
       ! Let's finish reading the other parameters
       CALL section_vals_val_get(coordination_section,"R0",r_val=colvar%coord_param%r_0,error=error)
       CALL section_vals_val_get(coordination_section,"NN",i_val=colvar%coord_param%nncrd,error=error)
       CALL section_vals_val_get(coordination_section,"ND",i_val=colvar%coord_param%ndcrd,error=error)
    ELSE IF (my_subsection(5)) THEN
       ! Distance point from plane
       wrk_section => plane_dist_section
       CALL colvar_create(colvar, plane_distance_colvar_id, error)
       CALL colvar_check_points(colvar, plane_dist_section, error)
       CALL section_vals_val_get(plane_dist_section,"ATOMS_PLANE",i_vals=iatms,error=error)
       CPPostcondition(SIZE(iatms) == 3,cp_failure_level,routinep,error,failure)
       colvar%plane_distance_param%plane = iatms
       CALL section_vals_val_get(plane_dist_section,"ATOM_POINT",i_val=iatm,error=error)
       colvar%plane_distance_param%point = iatm
    ELSE IF (my_subsection(6)) THEN
       ! Rotation colvar of a segment w.r.t. another segment
       wrk_section => rotation_section
       CALL colvar_create(colvar, rotation_colvar_id, error)
       CALL colvar_check_points(colvar, rotation_section, error)
       CALL section_vals_val_get(rotation_section,"P1_BOND1",i_val=colvar%rotation_param%i_at1_bond1,error=error)
       CALL section_vals_val_get(rotation_section,"P2_BOND1",i_val=colvar%rotation_param%i_at2_bond1,error=error)
       CALL section_vals_val_get(rotation_section,"P1_BOND2",i_val=colvar%rotation_param%i_at1_bond2,error=error)
       CALL section_vals_val_get(rotation_section,"P2_BOND2",i_val=colvar%rotation_param%i_at2_bond2,error=error)
    ELSE IF (my_subsection(7)) THEN
       ! Difference of two distances
       wrk_section => ddiff_section
       CALL colvar_create(colvar, ddiff_colvar_id, error)
       CALL colvar_check_points(colvar, ddiff_section, error)
       CALL section_vals_val_get(ddiff_section,"ATOMS",i_vals=iatms,error=error)
       colvar%ddiff_param%i_at_ddiff = iatms
    ELSE IF (my_subsection(8)) THEN
       ! Q Parameter
       wrk_section => qparm_section
       CALL colvar_create(colvar, qparm_colvar_id, error)
       CALL colvar_check_points(colvar, qparm_section, error)
       CALL section_vals_val_get(qparm_section,"RCUT",r_val=colvar%qparm_param%rcut,error=error)
       CALL section_vals_val_get(qparm_section,"ALPHA",r_val=colvar%qparm_param%alpha,error=error)
       CALL section_vals_val_get(qparm_section,"L",i_val=colvar%qparm_param%l,error=error)
       NULLIFY(colvar%qparm_param%i_at_from)
       NULLIFY(colvar%qparm_param%i_at_to)
       CALL section_vals_val_get(qparm_section,"ATOMS_FROM",n_rep_val=n_var,error=error)
       ndim = 0
       DO k = 1, n_var
          CALL section_vals_val_get(qparm_section,"ATOMS_FROM",i_rep_val=k,i_vals=iatms,error=error)
          CALL reallocate(colvar%qparm_param%i_at_from,1, ndim+SIZE(iatms))
          colvar%qparm_param%i_at_from(ndim+1:ndim+SIZE(iatms)) = iatms 
          ndim = ndim + SIZE(iatms)
       END DO
       colvar%qparm_param%n_atoms_from = ndim
       ! This section can be repeated
       CALL section_vals_val_get(qparm_section,"ATOMS_TO",n_rep_val=n_var,error=error)
       ndim = 0
       DO k = 1, n_var
          CALL section_vals_val_get(qparm_section,"ATOMS_TO",i_rep_val=k,i_vals=iatms,error=error)
          CALL reallocate(colvar%qparm_param%i_at_to,1, ndim+SIZE(iatms))
          colvar%qparm_param%i_at_to(ndim+1:ndim+SIZE(iatms)) = iatms 
          ndim = ndim + SIZE(iatms)
       END DO
       colvar%qparm_param%n_atoms_to = ndim
    ELSE IF (my_subsection(9)) THEN
       CALL colvar_create(colvar,hydronium_colvar_id, error)
       CALL section_vals_val_get(hydronium_section,"OXYGENS",i_vals=iatms,error=error)
       colvar%hydronium_param%n_oxygens = SIZE(iatms)
       ALLOCATE(colvar%hydronium_param%i_oxygens(SIZE(iatms)),stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       colvar%hydronium_param%i_oxygens = iatms
       CALL section_vals_val_get(hydronium_section,"HYDROGENS",i_vals=iatms,error=error)
       colvar%hydronium_param%n_hydrogens = SIZE(iatms)
       ALLOCATE(colvar%hydronium_param%i_hydrogens(SIZE(iatms)),stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       colvar%hydronium_param%i_hydrogens = iatms       
       CALL section_vals_val_get(hydronium_section,"ROO",r_val=colvar%hydronium_param%r_OO,error=error)
       CALL section_vals_val_get(hydronium_section,"ROH",r_val=colvar%hydronium_param%r_OH,error=error)
       CALL section_vals_val_get(hydronium_section,"pNH",i_val=colvar%hydronium_param%pnh,error=error)
       CALL section_vals_val_get(hydronium_section,"qNH",i_val=colvar%hydronium_param%qnh,error=error)
       CALL section_vals_val_get(hydronium_section,"pNO",i_val=colvar%hydronium_param%pno,error=error)
       CALL section_vals_val_get(hydronium_section,"qNO",i_val=colvar%hydronium_param%qno,error=error)
       CALL section_vals_val_get(hydronium_section,"p",i_val=colvar%hydronium_param%p,error=error)
       CALL section_vals_val_get(hydronium_section,"q",i_val=colvar%hydronium_param%q,error=error)
       CALL section_vals_val_get(hydronium_section,"NH",r_val=colvar%hydronium_param%nh,error=error)
       CALL section_vals_val_get(hydronium_section,"LAMBDA",r_val=colvar%hydronium_param%lambda,error=error)
    END IF

    CALL colvar_setup(colvar, error)

    iw = cp_print_key_unit_nr(logger,colvar_section,&
         "PRINT%PROGRAM_RUN_INFO",extension=".colvarLog",error=error)
    IF (iw>0) THEN
       tag = "ATOMS: "
       IF (colvar%use_points) tag = "POINTS:"
       WRITE ( iw, '( A )')'          ----------------------------------------------------------------------'
       SELECT CASE(colvar%type_id)
       CASE(angle_colvar_id)
          WRITE ( iw, '( A,T66,3I5)' )   ' COLVARS| ANGLE          >>> '//tag,&
               colvar%angle_param%i_at_angle
       CASE(ddiff_colvar_id)
          WRITE ( iw, '( A,T61,4I5)' )   ' COLVARS| DIFFERENCE BETWEEN TWO DISTANCES   >>> '//tag,&
               colvar%ddiff_param%i_at_ddiff
       CASE(plane_distance_colvar_id)
          WRITE ( iw, '( A,T66,3I5)' )   ' COLVARS| PLANE DISTANCE    - PLANE  >>> '//tag,&
               colvar%plane_distance_param%plane
          WRITE ( iw, '( A,T76,1I5)' )   ' COLVARS| PLANE DISTANCE    - POINT  >>> '//tag,&
               colvar%plane_distance_param%point
       CASE(torsion_colvar_id)
          WRITE ( iw, '( A,T61,4I5)' )   ' COLVARS| TORSION       >>> '//tag,&
               colvar%torsion_param%i_at_tors
       CASE(dist_colvar_id)
          WRITE ( iw, '( A,T71,2I5)' )   ' COLVARS| BOND          >>> '//tag,&
               colvar%dist_param%i_at,colvar%dist_param%j_at
       CASE(coord_colvar_id)
          WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COORDINATION  >>> FROM '//tag,&
               colvar%coord_param%i_at_from(kk),&
               kk=1,SIZE(colvar%coord_param%i_at_from))
          WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COORDINATION  >>>   TO '//tag,&
               colvar%coord_param%i_at_to(kk),&
               kk=1,SIZE(colvar%coord_param%i_at_to))
          WRITE ( iw, '( A,T71,F10.5)' ) ' COLVARS| R0',colvar%coord_param%r_0
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| NN',colvar%coord_param%nncrd
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| ND',colvar%coord_param%ndcrd
       CASE(rotation_colvar_id)
           WRITE ( iw, '( A,T71,I10)' )  ' COLVARS| BOND_ROTATION   - POINT 1 LINE 1  >>> '//tag,&
               colvar%rotation_param%i_at1_bond1
           WRITE ( iw, '( A,T71,I10)' )  ' COLVARS| BOND_ROTATION   - POINT 2 LINE 1  >>> '//tag,&
               colvar%rotation_param%i_at2_bond1         
           WRITE ( iw, '( A,T71,I10)' )  ' COLVARS| BOND_ROTATION   - POINT 1 LINE 2  >>> '//tag,&
               colvar%rotation_param%i_at1_bond2
           WRITE ( iw, '( A,T71,I10)' )  ' COLVARS| BOND_ROTATION   - POINT 2 LINE 2  >>> '//tag,&
               colvar%rotation_param%i_at2_bond2
       CASE(qparm_colvar_id)
          WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| Q-PARM  >>> FROM '//tag,&
               colvar%qparm_param%i_at_from(kk),&
               kk=1,SIZE(colvar%qparm_param%i_at_from))
          WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| Q-PARM  >>>   TO '//tag,&
               colvar%qparm_param%i_at_to(kk),&
               kk=1,SIZE(colvar%qparm_param%i_at_to))
          WRITE ( iw, '( A,T71,F10.5)' ) ' COLVARS| RCUT',colvar%qparm_param%rcut
          WRITE ( iw, '( A,T71,F10.5)' ) ' COLVARS| ALPHA',colvar%qparm_param%alpha
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| L',colvar%qparm_param%l
       END SELECT
       IF (colvar%use_points) THEN
          WRITE ( iw, '( A)') ' COLVARS| INFORMATION ON DEFINED GEOMETRICAL POINTS'
          DO kk = 1, SIZE(colvar%points)
             point_section => section_vals_get_subs_vals(wrk_section,"POINT",error=error)
             CALL section_vals_val_get(point_section,"TYPE",i_rep_section=kk,c_val=tmpStr,error=error)
             tmpStr2 = cp_to_string(kk)
             WRITE ( iw, '( A)') ' COLVARS|  POINT Nr.'//TRIM(tmpStr2)//' OF TYPE: '//TRIM(tmpStr)
             IF (ASSOCIATED(colvar%points(kk)%atoms)) THEN
                WRITE ( iw, '( A)') ' COLVARS|   ATOMS BUILDING THE GEOMETRICAL POINT'
                WRITE ( iw, '( A, I10)')(' COLVARS|   ATOM:',colvar%points(kk)%atoms(k),k=1,SIZE(colvar%points(kk)%atoms))
             ELSE
              WRITE ( iw, '( A,4X,3F12.6)') ' COLVARS|   XYZ POSITION OF FIXED POINT:', colvar%points(kk)%r
             END IF
          END DO
       END IF
       WRITE ( iw, '( A )')'          ----------------------------------------------------------------------'
    END IF
    CALL cp_print_key_finished_output(iw,logger,colvar_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)    

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE colvar_read

  !!****f* colvar/colvar_eval_mol_f [1.0] *
  !!
  !!   NAME
  !!     colvar_eval_mol_f
  !!
  !!   FUNCTION
  !!     evaluates the force (dsdr) given and due to the given colvar
  !!     variables in a molecular environment
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
  !!     Teodoro Laino
  !!
  !!*** **********************************************************************
  SUBROUTINE colvar_eval_mol_f(colvar, cell, particles, pos, fixd_list, error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: pos
    TYPE(fixd_constraint_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: fixd_list
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_eval_mol_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: colvar_ok, failure

    failure=.FALSE.

    colvar_ok=ASSOCIATED(colvar)
    CPAssert(colvar_ok,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(pos)) THEN
          DO i = 1, SIZE(colvar%i_atom)
             j = colvar%i_atom(i)
             particles(j)%r=pos(:,j)
          END DO
       END IF
       SELECT CASE(colvar%type_id)
       CASE (dist_colvar_id)
          CALL dist_colvar(colvar,cell,particles=particles,error=error)
       CASE (coord_colvar_id)
          CALL coord_colvar(colvar,cell,particles=particles, error=error)
       CASE (torsion_colvar_id)
          CALL torsion_colvar(colvar, cell, particles=particles, error=error)
       CASE (angle_colvar_id)
          CALL angle_colvar(colvar, cell, particles=particles, error=error)
       CASE (ddiff_colvar_id)
          CALL ddiff_colvar(colvar, cell, particles=particles, error=error)
       CASE (plane_distance_colvar_id)
          CALL plane_distance_colvar(colvar,cell, particles=particles, error=error)
       CASE (rotation_colvar_id)
          CALL rotation_colvar(colvar,cell, particles=particles, error=error)
       CASE (qparm_colvar_id)
          CALL qparm_colvar(colvar,cell, particles=particles, error=error)
       CASE (hydronium_colvar_id)
          CALL hydronium_colvar(colvar,cell,particles=particles, error=error)
       CASE default
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       ! Check for fixed atom constraints
       IF (PRESENT(fixd_list)) CALL check_fixed_atom_cns_colv(fixd_list, colvar)
    END IF

  END SUBROUTINE colvar_eval_mol_f
  !***************************************************************************

  !!****f* colvar/colvar_eval_glob_f [1.0] *
  !!
  !!   NAME
  !!     colvar_eval_glob_f
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
  SUBROUTINE colvar_eval_glob_f(isubsys,icolvar,force_env,error)
    INTEGER                                  :: isubsys, icolvar
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_eval_glob_f', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: colvar_ok, failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys

    failure=.FALSE.
    NULLIFY(subsys,cell,colvar)
    CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
    colvar_ok=ASSOCIATED(subsys(isubsys)%subsys%colvar_p)
    CPAssert(colvar_ok,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       colvar => subsys(isubsys)%subsys%colvar_p(icolvar)%colvar
       SELECT CASE(colvar%type_id)
       CASE (dist_colvar_id)
          CALL dist_colvar(colvar,cell,subsys=subsys(isubsys)%subsys,error=error)
       CASE (coord_colvar_id)
          CALL coord_colvar(colvar,cell,subsys=subsys(isubsys)%subsys, error=error)
       CASE (torsion_colvar_id)
          CALL torsion_colvar(colvar,cell,subsys=subsys(isubsys)%subsys, error=error)
       CASE (angle_colvar_id)
          CALL angle_colvar(colvar,cell,subsys=subsys(isubsys)%subsys, error=error)
       CASE (ddiff_colvar_id)
          CALL ddiff_colvar(colvar,cell,subsys=subsys(isubsys)%subsys, error=error)
       CASE (plane_distance_colvar_id)
          CALL plane_distance_colvar(colvar,cell,subsys=subsys(isubsys)%subsys, error=error)
       CASE (rotation_colvar_id)
          CALL rotation_colvar(colvar,cell,subsys=subsys(isubsys)%subsys, error=error)
       CASE (qparm_colvar_id)
          CALL qparm_colvar(colvar,cell,subsys=subsys(isubsys)%subsys, error=error)
       CASE (hydronium_colvar_id)
          CALL hydronium_colvar(colvar,cell,subsys=subsys(isubsys)%subsys, error=error)
       CASE default
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       ! Check for fixed atom constraints
       CALL check_fixed_atom_cns_colv(subsys(isubsys)%subsys%gci%fixd_list,colvar)
    END IF
  END SUBROUTINE colvar_eval_glob_f
  !***************************************************************************

  !!****f* colvar/get_coordinates [1.0] *
  !!
  !!   NAME
  !!     get_coordinates
  !!
  !!   FUNCTION
  !!     Get coordinates of atoms or of geoemtrical points
  !!
  !!   NOTES
  !!     -
  !!   AUTHOR
  !!     Teodoro Laino 03.2007 [created]
  !!
  !!*** **********************************************************************
  SUBROUTINE get_coordinates(colvar, i, ri, my_particles)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER, INTENT(IN)                      :: i
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: ri
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    IF (colvar%use_points) THEN
       CALL eval_point_pos(colvar%points(i),&
                           my_particles,&
                           ri)
    ELSE
       ri(:) = my_particles(i)%r(:)
    END IF
    
  END SUBROUTINE get_coordinates

  !!****f* colvar/put_derivative [1.0] *
  !!
  !!   NAME
  !!     put
  !!
  !!   FUNCTION
  !!     Transfer derivatives to ds/dr
  !!
  !!   NOTES
  !!     -
  !!   AUTHOR
  !!     Teodoro Laino 03.2007 [created]
  !!
  !!*** **********************************************************************
  SUBROUTINE put_derivative(colvar, i, fi)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER, INTENT(IN)                      :: i
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: fi

    IF (colvar%use_points) THEN
       CALL eval_point_der(colvar%points,&
                           i,&
                           colvar%dsdr,&
                           fi)
    ELSE
       colvar%dsdr(:,i) = colvar%dsdr(:,i) + fi
    END IF
    
  END SUBROUTINE put_derivative

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
  SUBROUTINE plane_distance_colvar(colvar,cell,subsys,particles,error)

    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'plane_distance_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, l
    LOGICAL                                  :: failure
    REAL(dp) :: a, b, dsdxpn(3), dxpndxi(3,3), dxpndxj(3,3), dxpndxk(3,3), &
      fi(3), fj(3), fk(3), fl(3), r12, ri(3), rj(3), rk(3), rl(3), ss(3), &
      xpij(3), xpkj(3), xpl(3), xpn(3)
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    NULLIFY(particles_i)

    CPPrecondition(colvar%type_id==plane_distance_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    colvar%dsdr=0.0_dp
    i=colvar%plane_distance_param%plane(1)
    j=colvar%plane_distance_param%plane(2)
    k=colvar%plane_distance_param%plane(3)
    l=colvar%plane_distance_param%point
    ! Get coordinates of atoms or points
    CALL get_coordinates(colvar, i, ri, my_particles)
    CALL get_coordinates(colvar, j, rj, my_particles)
    CALL get_coordinates(colvar, k, rk, my_particles)
    CALL get_coordinates(colvar, l, rl, my_particles)
    ! xpij
    ss=MATMUL(cell%h_inv,ri-rj)
    ss=ss-NINT(ss)
    xpij=MATMUL(cell%hmat,ss)
    ! xpkj
    ss=MATMUL(cell%h_inv,rk-rj)
    ss=ss-NINT(ss)
    xpkj=MATMUL(cell%hmat,ss)
    ! xpl
    ss=MATMUL(cell%h_inv,rl-(ri+rj+rk)/3.0_dp)
    ss=ss-NINT(ss)
    xpl=MATMUL(cell%hmat,ss)
    ! xpn
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
    fi(:)=MATMUL(dsdxpn,dxpndxi)-xpn/(3.0_dp*r12)
    fj(:)=MATMUL(dsdxpn,dxpndxj)-xpn/(3.0_dp*r12)
    fk(:)=MATMUL(dsdxpn,dxpndxk)-xpn/(3.0_dp*r12)
    fl(:)=xpn/r12
    ! Transfer derivatives on atoms
    CALL put_derivative(colvar, 1, fi)
    CALL put_derivative(colvar, 2, fj)
    CALL put_derivative(colvar, 3, fk)
    CALL put_derivative(colvar, 4, fl)

  END SUBROUTINE plane_distance_colvar

  !!****f* colvar/rotation_colvar [1.0] *
  !!
  !!   NAME
  !!     rotation_colvar
  !!
  !!   FUNCTION
  !!     Evaluates the value of the rotation angle between two bonds 
  !!
  !!   NOTES
  !!     -
  !!   AUTHOR
  !!     Teodoro Laino 02.2006 [created]
  !!
  !!*** **********************************************************************
  SUBROUTINE rotation_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotation_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, idum
    LOGICAL                                  :: failure
    REAL(dp)                                 :: a, b, fmod, t0, t1, t2, t3, &
                                                xdum(3), xij(3), xkj(3)
    REAL(KIND=dp)                            :: dp1b1(3), dp1b2(3), dp2b1(3), &
                                                dp2b2(3), ss(3), xp1b1(3), &
                                                xp1b2(3), xp2b1(3), xp2b2(3)
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    NULLIFY(particles_i)

    CPPrecondition(colvar%type_id==rotation_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF

    colvar%dsdr=0.0_dp
    i = colvar%rotation_param%i_at1_bond1
    CALL get_coordinates(colvar, i, xp1b1, my_particles)
    i = colvar%rotation_param%i_at2_bond1
    CALL get_coordinates(colvar, i, xp2b1, my_particles)
    i = colvar%rotation_param%i_at1_bond2
    CALL get_coordinates(colvar, i, xp1b2, my_particles)
    i = colvar%rotation_param%i_at2_bond2
    CALL get_coordinates(colvar, i, xp2b2, my_particles)
    ! xij
    ss=MATMUL(cell%h_inv,xp1b1-xp2b1)
    ss=ss-NINT(ss)
    xij=MATMUL(cell%hmat,ss)
    ! xkj
    ss=MATMUL(cell%h_inv,xp1b2-xp2b2)
    ss=ss-NINT(ss)
    xkj=MATMUL(cell%hmat,ss)
    ! evaluation of the angle..
    a = SQRT(DOT_PRODUCT(xij,xij))
    b = SQRT(DOT_PRODUCT(xkj,xkj))
    t0 = 1.0_dp/(a*b)
    t1 = 1.0_dp/(a**3.0_dp*b)
    t2 = 1.0_dp/(a*b**3.0_dp)
    t3 = DOT_PRODUCT(xij,xkj)
    colvar%ss = ACOS(t3*t0)
    IF ((ABS(colvar%ss).LT.0.00001_dp).OR.(ABS(colvar%ss-pi).LT.0.00001_dp)) THEN
       fmod =   0.0_dp
    ELSE
       fmod = - 1.0_dp / SIN(colvar%ss)
    ENDIF
    dp1b1 =  xkj(:)*t0-xij(:)*t1*t3
    dp2b1 = -xkj(:)*t0+xij(:)*t1*t3
    dp1b2 =  xij(:)*t0-xkj(:)*t2*t3
    dp2b2 = -xij(:)*t0+xkj(:)*t2*t3

    xdum = dp1b1 * fmod
    idum = colvar%rotation_param%i_at1_bond1
    CALL put_derivative(colvar, idum, xdum)
    xdum = dp2b1 * fmod
    idum = colvar%rotation_param%i_at2_bond1
    CALL put_derivative(colvar, idum, xdum)
    xdum = dp1b2 * fmod
    idum = colvar%rotation_param%i_at1_bond2
    CALL put_derivative(colvar, idum, xdum)
    xdum = dp2b2 * fmod
    idum = colvar%rotation_param%i_at2_bond2
    CALL put_derivative(colvar, idum, xdum)

  END SUBROUTINE rotation_colvar

  !!****f* colvar/ddiff_colvar [1.0] *
  !!
  !!   NAME
  !!     ddiff_colvar
  !!
  !!   FUNCTION
  !!     evaluates the force due to the difference of two distances
  !!
  !!   NOTES
  !!     -
  !!   AUTHOR
  !!     Teodoro Laino 02.2006 [created]
  !!
  !!*** **********************************************************************
  SUBROUTINE ddiff_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ddiff_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, l
    LOGICAL                                  :: failure
    REAL(dp)                                 :: fi(3), fj(3), fk(3), fl(3), &
                                                r12, r34, ss(3), xij(3), &
                                                xkl(3), xpi(3), xpj(3), &
                                                xpk(3), xpl(3)
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    NULLIFY(particles_i)

    CPPrecondition(colvar%type_id==ddiff_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    colvar%dsdr=0.0_dp
    i=colvar%ddiff_param%i_at_ddiff(1)
    j=colvar%ddiff_param%i_at_ddiff(2)
    ! First bond
    CALL get_coordinates(colvar, i, xpi, my_particles)
    CALL get_coordinates(colvar, j, xpj, my_particles)
    ss=MATMUL(cell%h_inv,xpi-xpj)
    ss=ss-NINT(ss)
    xij=MATMUL(cell%hmat,ss)
    r12=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
    ! Second bond
    k=colvar%ddiff_param%i_at_ddiff(3)
    l=colvar%ddiff_param%i_at_ddiff(4)
    CALL get_coordinates(colvar, k, xpk, my_particles)
    CALL get_coordinates(colvar, l, xpl, my_particles)
    ss=MATMUL(cell%h_inv,xpk-xpl)
    ss=ss-NINT(ss)
    xkl=MATMUL(cell%hmat,ss)
    r34=SQRT(xkl(1)**2+xkl(2)**2+xkl(3)**2)
    !
    colvar%ss=r12-r34
    fi(:)= xij/r12
    fj(:)=-xij/r12
    fk(:)=-xkl/r34
    fl(:)= xkl/r34
    CALL put_derivative(colvar, 1, fi)
    CALL put_derivative(colvar, 2, fj)
    CALL put_derivative(colvar, 3, fk)
    CALL put_derivative(colvar, 4, fl)
    
  END SUBROUTINE ddiff_colvar

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
  SUBROUTINE angle_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'angle_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k
    LOGICAL                                  :: failure
    REAL(dp)                                 :: a, b, fi(3), fj(3), fk(3), &
                                                fmod, ri(3), rj(3), rk(3), &
                                                ss(3), t0, t1, t2, t3, &
                                                xij(3), xkj(3)
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    NULLIFY(particles_i)

    CPPrecondition(colvar%type_id==angle_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    colvar%dsdr=0.0_dp
    i=colvar%angle_param%i_at_angle(1)
    j=colvar%angle_param%i_at_angle(2)
    k=colvar%angle_param%i_at_angle(3)
    CALL get_coordinates(colvar, i, ri, my_particles)
    CALL get_coordinates(colvar, j, rj, my_particles)
    CALL get_coordinates(colvar, k, rk, my_particles)
    ! xij
    ss=MATMUL(cell%h_inv,ri-rj)
    ss=ss-NINT(ss)
    xij=MATMUL(cell%hmat,ss)
    ! xkj
    ss=MATMUL(cell%h_inv,rk-rj)
    ss=ss-NINT(ss)
    xkj=MATMUL(cell%hmat,ss)
    ! Evaluation of the angle..
    a = SQRT(DOT_PRODUCT(xij,xij))
    b = SQRT(DOT_PRODUCT(xkj,xkj))
    t0 = 1.0_dp/(a*b)
    t1 = 1.0_dp/(a**3.0_dp*b)
    t2 = 1.0_dp/(a*b**3.0_dp)
    t3 = DOT_PRODUCT(xij,xkj)
    colvar%ss = ACOS(t3*t0)
    IF ((ABS(colvar%ss).LT.0.00001_dp).OR.(ABS(colvar%ss-pi).LT.0.00001_dp)) THEN
       fmod =   0.0_dp
    ELSE
       fmod = - 1.0_dp / SIN(colvar%ss)
    ENDIF
    fi(:) = xkj(:)*t0-xij(:)*t1*t3
    fj(:) =-xkj(:)*t0+xij(:)*t1*t3-xij(:)*t0+xkj(:)*t2*t3
    fk(:) = xij(:)*t0-xkj(:)*t2*t3
    fi = fi * fmod
    fj = fj * fmod
    fk = fk * fmod
    CALL put_derivative(colvar, 1, fi)
    CALL put_derivative(colvar, 2, fj)
    CALL put_derivative(colvar, 3, fk)

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
  SUBROUTINE dist_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dist_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure
    REAL(dp)                                 :: fi(3), fj(3), r12, ss(3), &
                                                xij(3), xpi(3), xpj(3)
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    NULLIFY(particles_i)

    CPPrecondition(colvar%type_id==dist_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    colvar%dsdr=0.0_dp
    i=colvar%dist_param%i_at
    j=colvar%dist_param%j_at
    CALL get_coordinates(colvar, i, xpi, my_particles)
    CALL get_coordinates(colvar, j, xpj, my_particles)
    ss=MATMUL(cell%h_inv,xpi-xpj)
    ss=ss-NINT(ss)
    xij=MATMUL(cell%hmat,ss)
    r12=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
    colvar%ss=r12
    fi(:)=xij/r12
    fj(:)=-xij/r12
    CALL put_derivative(colvar, 1, fi)
    CALL put_derivative(colvar, 2, fj)
    
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
  SUBROUTINE coord_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coord_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, jj, n_atoms_from, &
                                                n_atoms_to, ndcrd, nncrd
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dfunc, ftmp(3), func, &
                                                inv_n_atoms_from, invden, &
                                                ncoord, num, r12, r_0, rdist, &
                                                ss(3), xij(3)
    REAL(dp), DIMENSION(3)                   :: xpi, xpj
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    n_atoms_to=colvar%coord_param%n_atoms_to
    n_atoms_from=colvar%coord_param%n_atoms_from
    nncrd=colvar%coord_param%nncrd
    ndcrd=colvar%coord_param%ndcrd
    r_0=colvar%coord_param%r_0
    NULLIFY(particles_i)
    CPPrecondition(colvar%type_id==coord_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF

    ncoord=0.0_dp
    colvar%dsdr=0.0_dp
    inv_n_atoms_from=1.0_dp/REAL(n_atoms_from)
    DO ii=1,n_atoms_from
       i=colvar%coord_param%i_at_from(ii)
       CALL get_coordinates(colvar, i, xpi, my_particles)
       DO jj=1,n_atoms_to
          i=colvar%coord_param%i_at_to(jj)
          CALL get_coordinates(colvar, i, xpj, my_particles)       
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
          ftmp = dfunc*xij
          CALL put_derivative(colvar, ii, ftmp)
          ftmp = -dfunc*xij
          CALL put_derivative(colvar, n_atoms_from+jj, ftmp)
       ENDDO
    ENDDO
    colvar%ss=ncoord*inv_n_atoms_from
    colvar%dsdr(:,:) = colvar%dsdr(:,:)*inv_n_atoms_from

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
  SUBROUTINE torsion_colvar(colvar,cell,subsys,particles,error)

    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'torsion_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii
    LOGICAL                                  :: failure
    REAL(dp) :: angle, cosine, dedphi, dedxia, dedxib, dedxic, dedxid, dedxt, &
      dedxu, dedyia, dedyib, dedyic, dedyid, dedyt, dedyu, dedzia, dedzib, &
      dedzic, dedzid, dedzt, dedzu, dt, e, ftmp(3), o0, rcb, rt2, rtmp(3), &
      rtru, ru2, sine, ss(3), xba, xca, xcb, xdb, xdc, xt, xtu, xu, yba, yca, &
      ycb, ydb, ydc, yt, ytu, yu, zba, zca, zcb, zdb, zdc, zt, ztu, zu
    REAL(dp), DIMENSION(3, 4)                :: rr
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    NULLIFY(particles_i)
    CPPrecondition(colvar%type_id==torsion_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    colvar%dsdr=0.0_dp
    DO ii=1,4
       i=colvar%torsion_param%i_at_tors(ii)
       CALL get_coordinates(colvar, i, rtmp, my_particles)
       rr(:,ii)=rtmp(1:3)
    ENDDO
    o0 = colvar%torsion_param%o0
    ! ba
    ss=MATMUL(cell%h_inv,rr(:,2)-rr(:,1))
    ss=ss-NINT(ss)
    ss=MATMUL(cell%hmat,ss)
    xba = ss(1)
    yba = ss(2)
    zba = ss(3)
    ! cb
    ss=MATMUL(cell%h_inv,rr(:,3)-rr(:,2))
    ss=ss-NINT(ss)
    ss=MATMUL(cell%hmat,ss)
    xcb = ss(1)
    ycb = ss(2)
    zcb = ss(3)
    ! dc
    ss=MATMUL(cell%h_inv,rr(:,4)-rr(:,3))
    ss=ss-NINT(ss)
    ss=MATMUL(cell%hmat,ss)    
    xdc = ss(1)
    ydc = ss(2)
    zdc = ss(3)
    !
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
    IF (rtru .NE. 0.0_dp) THEN
       rcb = SQRT(xcb*xcb + ycb*ycb + zcb*zcb)
       cosine = (xt*xu + yt*yu + zt*zu) / rtru
       sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
       cosine = MIN(1.0_dp,MAX(-1.0_dp,cosine))
       angle  = ACOS(cosine)
       IF (sine .LT. 0.0_dp)  angle = -angle
       ! 
       dt = angle ! [rad]
       dt=MOD(2.0E4_dp*pi+dt-o0,2.0_dp*pi)
       IF ( dt .GT. pi ) dt = dt - 2.0_dp*pi
       dt=o0+dt
       colvar%torsion_param%o0 = dt
       !
       !     calculate improper energy and master chain rule term
       !
       e = dt
       dedphi = 1.0_dp
       !
       !     chain rule terms for first derivative components
       !
       ! ca
       ss=MATMUL(cell%h_inv,rr(:,3)-rr(:,1))
       ss=ss-NINT(ss)
       ss=MATMUL(cell%hmat,ss)
       xca = ss(1)
       yca = ss(2)
       zca = ss(3)   
       ! db
       ss=MATMUL(cell%h_inv,rr(:,4)-rr(:,2))
       ss=ss-NINT(ss)
       ss=MATMUL(cell%hmat,ss)
       xdb = ss(1)
       ydb = ss(2)
       zdb = ss(3)
       !
       dedxt =  dedphi * (yt*zcb - ycb*zt) / (rt2*rcb)
       dedyt =  dedphi * (zt*xcb - zcb*xt) / (rt2*rcb)
       dedzt =  dedphi * (xt*ycb - xcb*yt) / (rt2*rcb)
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
    ftmp(1)=dedxia
    ftmp(2)=dedyia
    ftmp(3)=dedzia
    CALL put_derivative(colvar, 1, ftmp)
    ftmp(1)=dedxib
    ftmp(2)=dedyib
    ftmp(3)=dedzib
    CALL put_derivative(colvar, 2, ftmp)
    ftmp(1)=dedxic
    ftmp(2)=dedyic
    ftmp(3)=dedzic
    CALL put_derivative(colvar, 3, ftmp)
    ftmp(1)=dedxid
    ftmp(2)=dedyid
    ftmp(3)=dedzid
    CALL put_derivative(colvar, 4, ftmp)
  END SUBROUTINE torsion_colvar

  !***************************************************************************

  !!****f* colvar/qparm_colvar [1.0] *
  !!
  !!   NAME
  !!     qparm_colvar
  !!
  !!   FUNCTION
  !!     evaluates the force due (and on) the Q PARM collective variable
  !!
  !!   NOTES
  !!
  !!   AUTHOR
  !!
  !!*** **********************************************************************
  SUBROUTINE qparm_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qparm_colvar', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: ylm, ylm_calc, ylm_calc2
    INTEGER                                  :: i, ii, j, jj, l, lmm, lpm, m, &
                                                n_atoms_from, n_atoms_to
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: alpha, bond, costheta, dplm, dylm, exp_fac, fact, fi, &
      ftmp(3), inv_n_atoms_from, nbond, plm, pre_fac, qlm, qlm2, qparm, r12, &
      r_tmp, rcut, ss(3), theta, x_tmp, xij(3), y_tmp, ymag, ymag1, z_tmp
    REAL(KIND=dp), DIMENSION(3)              :: grad, grad_nb, grad_nb_tot, &
                                                grad_tot, xpi, xpj
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

   failure=.FALSE.
   n_atoms_to=colvar%qparm_param%n_atoms_to
   n_atoms_from=colvar%qparm_param%n_atoms_from
   rcut=colvar%qparm_param%rcut
   l=colvar%qparm_param%l
   alpha=colvar%qparm_param%alpha
   NULLIFY(particles_i)
   CPPrecondition(colvar%type_id==qparm_colvar_id,cp_failure_level,routineP,error,failure)
   IF (PRESENT(particles)) THEN
      my_particles => particles
   ELSE
      CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
      CALL cp_subsys_get(subsys,particles=particles_i,error=error)
      my_particles => particles_i%els
   END IF

   qparm=0.0_dp
   colvar%dsdr=0.0_dp
   inv_n_atoms_from=1.0_dp/REAL(n_atoms_from)
   DO ii=1,n_atoms_from
      i=colvar%qparm_param%i_at_from(ii)
      CALL get_coordinates(colvar, i, xpi, my_particles)
      qlm = 0.0_dp
      qlm2 = 0.0_dp
      grad_tot(:) = 0.0_dp
      grad_nb_tot(:) = 0.0_dp
      DO m = 0, l 
        ylm = 0.0_dp
        grad(:) = 0.0_dp
        grad_nb(:) = 0.0_dp
        nbond = 0
        DO jj=1,n_atoms_to
           j=colvar%qparm_param%i_at_to(jj)
           IF (i==j) CYCLE
           CALL get_coordinates(colvar, j, xpj, my_particles)       
           ss=MATMUL(cell%h_inv,xpj(:)-xpi(:))
           ss=ss-NINT(ss)
           xij=MATMUL(cell%hmat,ss)
           r12=SQRT(DOT_PRODUCT(xij,xij))
           IF(r12 < rcut) nbond = nbond + 1.0_dp
           x_tmp = xij(1)
           y_tmp = xij(2)
           z_tmp = xij(3)
           r_tmp = r12 
           IF (ABS(x_tmp) .GT. 1.0E-8_dp) THEN
             fi = ATAN(y_tmp/x_tmp)
           ELSE
             fi = 0.5_dp*pi
           ENDIF
           costheta = z_tmp/r_tmp
           IF(costheta > 1.0_dp) costheta = 1.0_dp
           IF(costheta < -1.0_dp) costheta = -1.0_dp
           theta = ACOS(costheta)
           plm = legendre (costheta, l, m)
           dplm = dlegendre (costheta, l, m)
           IF((l+m) .GT. maxfac) THEN
             CALL stop_program("qparm_colvar","l-m GT maxfac")
           ENDIF
           lmm = fac(l-m)
           lpm = fac(l+m)
           bond  = 1.0_dp/(1.0_dp + EXP(alpha*(r_tmp - rcut)))
           IF (bond > 1.0_dp) THEN 
             CALL stop_program("qparm_colvar","bond GT 1.0_dp")
           ENDIF
           pre_fac = bond*SQRT(((2*l+1)*lmm)/(4*pi*lpm))
           IF (plm < 0.0_dp) THEN
             dylm = -pre_fac*dplm !*plm/dabs(plm)
           ELSE
             dylm = pre_fac*dplm
           ENDIF  
           ylm_calc = pre_fac*plm*(CMPLX(COS(m*fi),SIN(m*fi)))
           ylm = ylm + ylm_calc
           ylm_calc2 = pre_fac/bond*plm*(CMPLX(COS(m*fi),SIN(m*fi)))
           ymag = (SQRT(ylm_calc2*CONJG(ylm_calc2)))
           ! Fi = -dUi/dr, but r_tmp = rj - ri, so have double negative
           grad(1) = grad(1) - dylm*(z_tmp*x_tmp)/(r_tmp*r_tmp*r_tmp)
           grad(2) = grad(2) - dylm*(z_tmp*y_tmp)/(r_tmp*r_tmp*r_tmp)
           grad(3) = grad(3) + dylm*((1.0_dp/r_tmp) - (z_tmp**2/r_tmp**3))
           exp_fac = alpha*EXP(alpha*(r_tmp-rcut))/((1.0_dp+EXP(alpha*(r_tmp-rcut)))**2)
           grad_nb(1) = grad_nb(1) - ymag*exp_fac*(x_tmp/r_tmp)
           grad_nb(2) = grad_nb(2) - ymag*exp_fac*(y_tmp/r_tmp)
           grad_nb(3) = grad_nb(3) - ymag*exp_fac*(z_tmp/r_tmp)
        ENDDO
        IF (m .GT. 0) THEN
          fact = 2.0_dp
        ELSE
          fact = 1.0_dp
        ENDIF
        IF (nbond < 0.0001_dp) THEN
          CALL stop_program("qparm_colvar","NBOND = 0.0")
        ENDIF
        ylm = ylm/nbond
        grad(:) = grad(:)/nbond
        grad_nb(:) = grad_nb(:)/nbond
        ymag1 = ylm*CONJG(ylm)
        qlm = qlm + fact*ymag1
        grad_tot(:) = grad_tot(:) + fact*grad(:)
        grad_nb_tot(:) = grad_nb_tot(:) + fact*grad_nb(:)
      ENDDO
      pre_fac = (4.0_dp*pi)/(2.0_dp*l + 1) 
      qlm = SQRT(pre_fac*qlm)
      qlm2 = SQRT(pre_fac*qlm2)
      qparm=qparm+qlm
      ftmp(:) = SQRT(pre_fac)*(grad_tot(:) + grad_nb_tot(:))
      CALL put_derivative(colvar, ii, ftmp)
   ENDDO
   colvar%ss=qparm*inv_n_atoms_from
   colvar%dsdr(:,:) = colvar%dsdr(:,:)*inv_n_atoms_from

  END SUBROUTINE qparm_colvar

  !!****f* colvar/hydronium_colvar [1.0] *
  !!
  !!   NAME
  !!     hydronium_colvar
  !!
  !!   FUNCTION
  !!     evaluates the force due (and on) the hydronium collective variable
  !!
  !!   NOTES
  !!
  !!   AUTHOR
  !!     Marcel Baer
  !!
  !!*** **********************************************************************
  SUBROUTINE hydronium_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hydronium_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, jj, n_hydrogens, &
                                                n_oxygens, p, pnh, pno, q, &
                                                qnh, qno, stat
    LOGICAL                                  :: failure
    REAL(dp)                                 :: fscalar, invden, lambda, nh, &
                                                num, qtot, r12, r_OH, r_OO, &
                                                rdist, ss(3), xij(3)
    REAL(dp), DIMENSION(3)                   :: xpi, xpj
    REAL(dp), DIMENSION(:)                   :: M, nhcoord, no, qloc

! DEBUG
! DEBUG

    ALLOCATABLE                              :: nhcoord,M,&
                                                qloc,no
    REAL(dp), DIMENSION(:,:,:)               :: dfunc_nh,dfunc_no,dM
    ALLOCATABLE                              :: dfunc_nh,dfunc_no,dM
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    n_oxygens=colvar%hydronium_param%n_oxygens
    n_hydrogens=colvar%hydronium_param%n_hydrogens
    nh=colvar%hydronium_param%nh
    pnh=colvar%hydronium_param%pnh
    qnh=colvar%hydronium_param%qnh
    pno=colvar%hydronium_param%pno
    qno=colvar%hydronium_param%qno
    r_OO=colvar%hydronium_param%r_OO
    r_OH=colvar%hydronium_param%r_OH
    lambda=colvar%hydronium_param%lambda
    p=colvar%hydronium_param%p
    q=colvar%hydronium_param%q


    NULLIFY(particles_i)
    CPPrecondition(colvar%type_id==hydronium_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF


    ALLOCATE ( dfunc_nh ( 3, n_oxygens, n_hydrogens ), stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE ( nhcoord  ( n_oxygens), stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE ( M (  n_oxygens ), stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE ( dM ( 3, n_oxygens,n_hydrogens ), stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)


    ALLOCATE ( dfunc_no ( 3, n_oxygens, n_oxygens ), stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE ( no       (  n_oxygens), stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    
    ALLOCATE ( qloc (  n_oxygens), stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

! Zero Arrays:
    colvar%dsdr=0.0_dp
    dfunc_nh = 0._dp 
    dfunc_no = 0._dp 
    M = 0._dp 
    dM = 0._dp 
    no = 0._dp 
    qloc = 0._dp 
    nhcoord = 0._dp
    DO ii=1,n_oxygens
       i=colvar%hydronium_param%i_oxygens(ii)
       xpi(:)=my_particles(i)%r(1:3)
! Computing M( n ( ii ) ) 
       DO jj=1,n_hydrogens
          i=colvar%hydronium_param%i_hydrogens(jj)
          xpj(:)=my_particles(i)%r(1:3)
          ss=MATMUL(cell%h_inv,xpi(:)-xpj(:))
          ss=ss-NINT(ss)
          xij=MATMUL(cell%hmat,ss)
          r12=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
          rdist = r12/r_OH
          num=(1.0_dp-rdist**pnh)
          invden=1.0_dp/(1.0_dp-rdist**qnh)
          IF ( ABS ( invden )  < 1.e-10_dp )  invden = 1.e-10_dp
          fscalar = ((-pnh*(rdist**(pnh-1))*invden) &
                     + num*(invden)**2*qnh*(rdist**(qnh-1)))/(r12*r_OH) 

          dfunc_nh(1,ii,jj)= xij (1) * fscalar
          dfunc_nh(2,ii,jj)= xij (2) * fscalar
          dfunc_nh(3,ii,jj)= xij (3) * fscalar

          nhcoord(ii)=nhcoord(ii) + num*invden 
       END DO 
       M(ii)  = 1.0_dp-(1.0_dp - (nhcoord(ii)/nh)**p)/&
                       (1.0_dp - (nhcoord(ii)/nh)**q) 

! Computing no ( ii )
       DO jj=1,n_oxygens
          i=colvar%hydronium_param%i_oxygens(jj)
          xpj(:)=my_particles(i)%r(1:3)
          ss=MATMUL(cell%h_inv,xpi(:)-xpj(:))
          ss=ss-NINT(ss)
          xij=MATMUL(cell%hmat,ss)
          r12=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
          IF ( r12 < 1.e-3_dp ) CYCLE
          rdist = r12/r_OO
          num=(1.0_dp-rdist**pno)
          invden=1.0_dp/(1.0_dp-rdist**qno)
          IF ( ABS ( invden )  < 1.e-10_dp )  invden = 1.e-10_dp
          fscalar = ((-pno*(rdist**(pno-1))*invden) &
                     + num*(invden)**2*qno*(rdist**(qno-1)))/(r12*r_OO) 
          
          dfunc_no(1,ii,jj)= xij(1)*fscalar
          dfunc_no(2,ii,jj)= xij(2)*fscalar
          dfunc_no(3,ii,jj)= xij(3)*fscalar

          no(ii)=no(ii)+ num*invden 
       END DO 
    END DO 
! no(ii) done

! computing qloc and Q
    qtot = 0._dp
    DO ii = 1, n_oxygens
     qloc ( ii ) = EXP ( lambda * M (ii) * no ( ii ) )
     qtot = qtot + qloc ( ii ) 
    END DO 
! compute forces
    DO ii = 1, n_oxygens
! Computing f_OH
       DO jj=1,n_hydrogens
          dM(1,ii,jj)=(p*((nhcoord(ii)/nh)**(p-1))*dfunc_nh(1,ii,jj))/nh/ &
                      (1.0_dp - (nhcoord(ii)/nh)**q) - &
                      (1.0_dp - (nhcoord(ii)/nh)**p)/ &
                      ((1.0_dp - (nhcoord(ii)/nh)**q)**2) *&
                      q*dfunc_nh(1,ii,jj)*(nhcoord(ii)/nh)**(q-1)/nh
          dM(2,ii,jj)=(p*((nhcoord(ii)/nh)**(p-1))*dfunc_nh(2,ii,jj))/nh/ &
                      (1.0_dp - (nhcoord(ii)/nh)**q) - &
                      (1.0_dp - (nhcoord(ii)/nh)**p)/ &
                      ((1.0_dp - (nhcoord(ii)/nh)**q)**2) *&
                      q*dfunc_nh(2,ii,jj)*(nhcoord(ii)/nh)**(q-1)/nh
          dM(3,ii,jj)=(p*((nhcoord(ii)/nh)**(p-1))*dfunc_nh(3,ii,jj))/nh/ &
                      (1.0_dp - (nhcoord(ii)/nh)**q) - &
                      (1.0_dp - (nhcoord(ii)/nh)**p)/ &
                      ((1.0_dp - (nhcoord(ii)/nh)**q)**2) *&
                      q*dfunc_nh(3,ii,jj)*(nhcoord(ii)/nh)**(q-1)/nh

          colvar%dsdr(1,ii)=colvar%dsdr(1,ii)+qloc(ii)*dM(1,ii,jj)*no(ii)/qtot
          colvar%dsdr(1,n_oxygens+jj)=colvar%dsdr(1,n_oxygens+jj) &
                                              -qloc(ii)*dM(1,ii,jj)*no(ii)/qtot
          colvar%dsdr(2,ii)=colvar%dsdr(2,ii)+qloc(ii)*dM(2,ii,jj)*no(ii)/qtot
          colvar%dsdr(2,n_oxygens+jj)=colvar%dsdr(2,n_oxygens+jj) &
                                              -qloc(ii)*dM(2,ii,jj)*no(ii)/qtot
          colvar%dsdr(3,ii)=colvar%dsdr(3,ii)+qloc(ii)*dM(3,ii,jj)*no(ii)/qtot
          colvar%dsdr(3,n_oxygens+jj)=colvar%dsdr(3,n_oxygens+jj) &
                                              -qloc(ii)*dM(3,ii,jj)*no(ii)/qtot
       END DO 
! Computing f_OO
       DO jj=1,n_oxygens
          colvar%dsdr(1,ii)=colvar%dsdr(1,ii)+qloc(ii)*M(ii)*dfunc_no(1,ii,jj)/qtot
          colvar%dsdr(1,jj)=colvar%dsdr(1,jj) &
                                              -qloc(ii)*M(ii)*dfunc_no(1,ii,jj)/qtot
          colvar%dsdr(2,ii)=colvar%dsdr(2,ii)+qloc(ii)*M(ii)*dfunc_no(2,ii,jj)/qtot
          colvar%dsdr(2,jj)=colvar%dsdr(2,jj) &
                                              -qloc(ii)*M(ii)*dfunc_no(2,ii,jj)/qtot
          colvar%dsdr(3,ii)=colvar%dsdr(3,ii)+qloc(ii)*M(ii)*dfunc_no(3,ii,jj)/qtot
          colvar%dsdr(3,jj)=colvar%dsdr(3,jj) &
                                              -qloc(ii)*M(ii)*dfunc_no(3,ii,jj)/qtot
       END DO 
    END DO 
    
    colvar%ss=LOG(qtot)/lambda
    DEALLOCATE ( dfunc_nh, stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE ( nhcoord, stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE ( M, stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE ( dM, stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE ( dfunc_no, stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE ( no, stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE ( qloc, stat=stat )
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
  END SUBROUTINE hydronium_colvar

END MODULE colvar_methods
