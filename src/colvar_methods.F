!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief defines collective variables s({R}) and the derivative of this variable wrt R
!>      these can then be used in constraints, restraints and metadynamics ...
!> \par History
!>      04.2004 created
!>      01.2006 Refactored [Joost VandeVondele]
!> \author Alessandro Laio,Fawzi Mohamed
! *****************************************************************************
MODULE colvar_methods
  USE cell_types,                      ONLY: cell_type
  USE colvar_types,                    ONLY: &
       angle_colvar_id, colvar_check_points, colvar_create, colvar_setup, &
       colvar_type, combine_colvar_id, coord_colvar_id, dfunct_colvar_id, &
       dist_colvar_id, eval_point_der, eval_point_pos, hydronium_colvar_id, &
       plane_distance_colvar_id, qparm_colvar_id, reaction_path_colvar_id, &
       rotation_colvar_id, torsion_colvar_id
  USE constraint_fxd,                  ONLY: check_fixed_atom_cns_colv
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE fparser,                         ONLY: evalf,&
                                             evalfd,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE input_constants,                 ONLY: do_clv_x,&
                                             do_clv_xy,&
                                             do_clv_xz,&
                                             do_clv_y,&
                                             do_clv_yz,&
                                             do_clv_z
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
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
  USE string_utilities,                ONLY: compress,&
                                             uppercase
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

CONTAINS

! *****************************************************************************
!> \brief reads a colvar from the input
!> \param colvar the place where to store what will be read
!> \param extended_lagrange if the extended lagrangian formalismus should 
!>        be used
!> \param icol number of the current colvar (repetition in colvar_section)
!> \param colvar_section the colvar section
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      04.2004 created [alessandro laio and fawzi mohamed]
!> \author teo
! *****************************************************************************
  RECURSIVE SUBROUTINE colvar_read(colvar, icol, colvar_section, error)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER, INTENT(IN)                      :: icol
    TYPE(section_vals_type), POINTER         :: colvar_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_read', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=7)                         :: tag
    CHARACTER(LEN=default_path_length)       :: path_function
    CHARACTER(LEN=default_string_length)     :: tmpStr, tmpStr2
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: c_kinds, my_par
    INTEGER                                  :: handle, i, iatm, iend, isize, &
                                                istart, iw, iw1, j, k, kk, &
                                                n_var, n_var_k, ncol, ndim, &
                                                stat, v_count
    INTEGER, DIMENSION(:), POINTER           :: iatms
    INTEGER, DIMENSION(:, :), POINTER        :: p_bounds
    LOGICAL                                  :: check, failure
    LOGICAL, DIMENSION(11)                   :: my_subsection
    REAL(dp), DIMENSION(:), POINTER          :: s1
    REAL(dp), DIMENSION(:, :), POINTER       :: p_range, s1v
    REAL(KIND=dp)                            :: coeff
    REAL(KIND=dp), DIMENSION(1)              :: my_val
    REAL(KIND=dp), DIMENSION(:), POINTER     :: g_range, grid_point, grid_sp, &
                                                my_vals, range
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER :: angle_section, colvar_subsection, &
      combine_section, coordination_section, dfunct_section, &
      distance_section, hydronium_section, plane_dist_section, point_section, &
      qparm_section, reaction_path_section, rotation_section, &
      torsion_section, wrk_section

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(logger, c_kinds, iatms)
    logger => cp_error_get_logger(error)
    my_subsection        = .FALSE.
    distance_section     => section_vals_get_subs_vals(colvar_section,"DISTANCE",i_rep_section=icol,error=error)
    dfunct_section       => section_vals_get_subs_vals(colvar_section,"DISTANCE_FUNCTION",&
                            i_rep_section=icol,error=error)
    angle_section        => section_vals_get_subs_vals(colvar_section,"ANGLE",i_rep_section=icol,error=error)
    torsion_section      => section_vals_get_subs_vals(colvar_section,"TORSION",i_rep_section=icol,error=error)
    coordination_section => section_vals_get_subs_vals(colvar_section,"COORDINATION",i_rep_section=icol,error=error)
    plane_dist_section   => section_vals_get_subs_vals(colvar_section,"DISTANCE_POINT_PLANE",i_rep_section=icol,&
                            error=error)
    rotation_section     => section_vals_get_subs_vals(colvar_section,"BOND_ROTATION",i_rep_section=icol,error=error)
    qparm_section        => section_vals_get_subs_vals(colvar_section,"QPARM",i_rep_section=icol,error=error)
    hydronium_section    => section_vals_get_subs_vals(colvar_section,"HYDRONIUM",i_rep_section=icol,error=error)
    reaction_path_section=> section_vals_get_subs_vals(colvar_section,"REACTION_PATH",i_rep_section=icol,&
                            can_return_null=.TRUE.,error=error)
    combine_section      => section_vals_get_subs_vals(colvar_section,"COMBINE_COLVAR",i_rep_section=icol,&
                            can_return_null=.TRUE.,error=error)

    CALL section_vals_get(distance_section,     explicit=my_subsection( 1), error=error)
    CALL section_vals_get(angle_section,        explicit=my_subsection( 2), error=error)
    CALL section_vals_get(torsion_section,      explicit=my_subsection( 3), error=error)
    CALL section_vals_get(coordination_section, explicit=my_subsection( 4), error=error)
    CALL section_vals_get(plane_dist_section,   explicit=my_subsection( 5), error=error)
    CALL section_vals_get(rotation_section,     explicit=my_subsection( 6), error=error)
    CALL section_vals_get(dfunct_section,       explicit=my_subsection( 7), error=error)
    CALL section_vals_get(qparm_section,        explicit=my_subsection( 8), error=error)
    CALL section_vals_get(hydronium_section,    explicit=my_subsection( 9), error=error)
    ! These are just special cases since they are not present in their own defition of COLVARS
    IF (ASSOCIATED(reaction_path_section)) THEN
       CALL section_vals_get(reaction_path_section,explicit=my_subsection(10), error=error)
    END IF
    IF (ASSOCIATED(combine_section)) THEN
       CALL section_vals_get(combine_section,explicit=my_subsection(11), error=error)
    END IF

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
       CALL section_vals_val_get(distance_section,"AXIS",i_val=colvar%dist_param%axis_id,error=error)
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
       NULLIFY(colvar%coord_param%i_at_from, colvar%coord_param%c_kinds_from)
       NULLIFY(colvar%coord_param%i_at_to, colvar%coord_param%c_kinds_to)
       NULLIFY(colvar%coord_param%i_at_to_b, colvar%coord_param%c_kinds_to_b)
       ! This section can be repeated
       CALL section_vals_val_get(coordination_section,"ATOMS_FROM",n_rep_val=n_var,error=error)
       ndim = 0
       IF (n_var /= 0) THEN
          ! INDEX LIST
          DO k = 1, n_var
             CALL section_vals_val_get(coordination_section,"ATOMS_FROM",i_rep_val=k,i_vals=iatms,error=error)
             CALL reallocate(colvar%coord_param%i_at_from,1, ndim+SIZE(iatms))
             colvar%coord_param%i_at_from(ndim+1:ndim+SIZE(iatms)) = iatms 
             ndim = ndim + SIZE(iatms)
          END DO
          colvar%coord_param%n_atoms_from = ndim
          colvar%coord_param%use_kinds_from = .FALSE.
       ELSE
          ! KINDS
          CALL section_vals_val_get(coordination_section,"KINDS_FROM",n_rep_val=n_var,error=error)
          CPPostcondition(n_var>0,cp_failure_level,routinep,error,failure)
          DO k = 1, n_var
             CALL section_vals_val_get(coordination_section,"KINDS_FROM",i_rep_val=k,c_vals=c_kinds,error=error)
             CALL reallocate(colvar%coord_param%c_kinds_from,1, ndim+SIZE(c_kinds))
             colvar%coord_param%c_kinds_from(ndim+1:ndim+SIZE(c_kinds)) = c_kinds
             ndim = ndim + SIZE(c_kinds)
          END DO
          colvar%coord_param%n_atoms_from = 0
          colvar%coord_param%use_kinds_from = .TRUE.
          ! Uppercase the label
          DO k = 1, ndim
             CALL uppercase(colvar%coord_param%c_kinds_from(k))
          END DO
       END IF
       ! This section can be repeated
       CALL section_vals_val_get(coordination_section,"ATOMS_TO",n_rep_val=n_var,error=error)
       ndim = 0
       IF (n_var /= 0) THEN
          ! INDEX LIST
          DO k = 1, n_var
             CALL section_vals_val_get(coordination_section,"ATOMS_TO",i_rep_val=k,i_vals=iatms,error=error)
             CALL reallocate(colvar%coord_param%i_at_to,1, ndim+SIZE(iatms))
             colvar%coord_param%i_at_to(ndim+1:ndim+SIZE(iatms)) = iatms 
             ndim = ndim + SIZE(iatms)
          END DO
          colvar%coord_param%n_atoms_to = ndim
          colvar%coord_param%use_kinds_to = .FALSE.
       ELSE
          ! KINDS
          CALL section_vals_val_get(coordination_section,"KINDS_TO",n_rep_val=n_var,error=error)
          CPPostcondition(n_var>0,cp_failure_level,routinep,error,failure)
          DO k = 1, n_var
             CALL section_vals_val_get(coordination_section,"KINDS_TO",i_rep_val=k,c_vals=c_kinds,error=error)
             CALL reallocate(colvar%coord_param%c_kinds_to,1, ndim+SIZE(c_kinds))
             colvar%coord_param%c_kinds_to(ndim+1:ndim+SIZE(c_kinds)) = c_kinds
             ndim = ndim + SIZE(c_kinds)
          END DO
          colvar%coord_param%n_atoms_to = 0
          colvar%coord_param%use_kinds_to = .TRUE.
          ! Uppercase the label
          DO k = 1, ndim
             CALL uppercase(colvar%coord_param%c_kinds_to(k))
          END DO
       END IF
       ! Let's finish reading the other parameters
       CALL section_vals_val_get(coordination_section,"R0",r_val=colvar%coord_param%r_0,error=error)
       CALL section_vals_val_get(coordination_section,"NN",i_val=colvar%coord_param%nncrd,error=error)
       CALL section_vals_val_get(coordination_section,"ND",i_val=colvar%coord_param%ndcrd,error=error)
       ! This section can be repeated
       CALL section_vals_val_get(coordination_section,"ATOMS_TO_B",n_rep_val=n_var,error=error)
       CALL section_vals_val_get(coordination_section,"KINDS_TO_B",n_rep_val=n_var_k,error=error)
       ndim = 0
       IF (n_var /= 0 .OR. n_var_k /= 0) THEN
          colvar%coord_param%do_chain = .TRUE.
          IF (n_var /= 0) THEN
             ! INDEX LIST
             DO k = 1, n_var
                CALL section_vals_val_get(coordination_section,"ATOMS_TO_B",i_rep_val=k,i_vals=iatms,error=error)
                CALL reallocate(colvar%coord_param%i_at_to_b,1, ndim+SIZE(iatms))
                colvar%coord_param%i_at_to_b(ndim+1:ndim+SIZE(iatms)) = iatms 
                ndim = ndim + SIZE(iatms)
             END DO
             colvar%coord_param%n_atoms_to_b = ndim
             colvar%coord_param%use_kinds_to_b = .FALSE.
          ELSE
             ! KINDS
             CALL section_vals_val_get(coordination_section,"KINDS_TO_B",n_rep_val=n_var_k,error=error)
             CPPostcondition(n_var_k>0,cp_failure_level,routinep,error,failure)
             DO k = 1, n_var_k
                CALL section_vals_val_get(coordination_section,"KINDS_TO_B",i_rep_val=k,c_vals=c_kinds,error=error)
                CALL reallocate(colvar%coord_param%c_kinds_to_b,1, ndim+SIZE(c_kinds))
                colvar%coord_param%c_kinds_to_b(ndim+1:ndim+SIZE(c_kinds)) = c_kinds
                ndim = ndim + SIZE(c_kinds)
             END DO
             colvar%coord_param%n_atoms_to_b = 0
             colvar%coord_param%use_kinds_to_b = .TRUE.
             ! Uppercase the label
             DO k = 1, ndim
                CALL uppercase(colvar%coord_param%c_kinds_to_b(k))
             END DO
          END IF
          ! Let's finish reading the other parameters
          CALL section_vals_val_get(coordination_section,"R0_B",r_val=colvar%coord_param%r_0_b,error=error)
          CALL section_vals_val_get(coordination_section,"NN_B",i_val=colvar%coord_param%nncrd_b,error=error)
          CALL section_vals_val_get(coordination_section,"ND_B",i_val=colvar%coord_param%ndcrd_b,error=error)
       ELSE
          colvar%coord_param%do_chain = .FALSE.
          colvar%coord_param%n_atoms_to_b = 0
          colvar%coord_param%use_kinds_to_b = .FALSE.
          NULLIFY(colvar%coord_param%i_at_to_b)
          NULLIFY(colvar%coord_param%c_kinds_to_b)
          colvar%coord_param%nncrd_b = 0
          colvar%coord_param%ndcrd_b = 0
          colvar%coord_param%r_0_b = 0._dp
       END IF

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
       wrk_section => dfunct_section
       CALL colvar_create(colvar, dfunct_colvar_id, error)
       CALL colvar_check_points(colvar, dfunct_section, error)
       CALL section_vals_val_get(dfunct_section,"ATOMS",i_vals=iatms,error=error)
       colvar%dfunct_param%i_at_dfunct = iatms
       CALL section_vals_val_get(dfunct_section,"COEFFICIENT",r_val=colvar%dfunct_param%coeff,error=error)
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
    ELSE IF(my_subsection(10))THEN
       !reaction path
       CALL colvar_create(colvar,reaction_path_colvar_id, error)
       colvar%use_points=.FALSE.
       colvar_subsection => section_vals_get_subs_vals(reaction_path_section,"COLVAR",error=error)
       CALL section_vals_get(colvar_subsection,n_repetition=ncol,error=error)
       ALLOCATE(colvar%reaction_path_param%colvar_p(ncol),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.NOT. failure) THEN
          DO i= 1, ncol
             NULLIFY(colvar%reaction_path_param%colvar_p(i)%colvar)
             CALL colvar_read(colvar%reaction_path_param%colvar_p(i)%colvar,i,colvar_subsection, error=error)
          ENDDO
       ENDIF
       colvar%reaction_path_param%n_components=ncol
       NULLIFY(range)
       CALL section_vals_val_get(reaction_path_section,"RANGE",r_vals=range,error=error)
       CALL section_vals_val_get(reaction_path_section,"STEP_SIZE",r_val=colvar%reaction_path_param%step_size,error=error)     
       CALL section_vals_val_get(reaction_path_section,"LAMBDA",r_val=colvar%reaction_path_param%lambda,error=error)
       iend=CEILING(MAX(RANGE(1),RANGE(2))/colvar%reaction_path_param%step_size)
       istart=FLOOR(MIN(RANGE(1),RANGE(2))/colvar%reaction_path_param%step_size)
       colvar%reaction_path_param%function_bounds(1)=istart
       colvar%reaction_path_param%function_bounds(2)=iend
       ALLOCATE(colvar%reaction_path_param%f_vals(ncol,istart:iend),stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       CALL section_vals_val_get(reaction_path_section,"VARIABLE",c_vals=my_par,i_rep_val=1,error=error)
       CALL section_vals_val_get(reaction_path_section,"FUNCTION",n_rep_val=ncol,error=error)
       check = (ncol==SIZE(colvar%reaction_path_param%colvar_p))
       CPPostcondition(check,cp_failure_level,routinep,error,failure)     
       CALL initf(ncol)
       DO i=1,ncol
          CALL section_vals_val_get(reaction_path_section,"FUNCTION",c_val=path_function,i_rep_val=i,error=error)
          CALL compress(path_function, full=.TRUE.)
          CALL parsef(i,TRIM(path_function),my_par)
          DO j=istart,iend
             my_val=REAL(j,kind=dp)*colvar%reaction_path_param%step_size
             colvar%reaction_path_param%f_vals(i,j)=evalf(i,my_val)
          END DO
       END DO
       CALL finalizef()

       iw1= cp_print_key_unit_nr(logger,reaction_path_section,&
            "MAP",extension=".dat",file_status="REPLACE",error=error)
       IF(iw1>0)THEN
          CALL section_vals_val_get(reaction_path_section,"MAP%GRID_SPACING",n_rep_val=ncol,error=error)
          ALLOCATE(grid_sp(ncol),stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          DO i=1,ncol
             CALL section_vals_val_get(reaction_path_section,"MAP%GRID_SPACING",r_val=grid_sp(i),error=error)
          END DO
          CALL section_vals_val_get(reaction_path_section,"MAP%RANGE",n_rep_val=ncol,error=error)
          CPPostcondition(ncol ==SIZE(grid_sp) ,cp_failure_level,routinep,error,failure)
          ALLOCATE(p_range(2,ncol),stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          ALLOCATE(p_bounds(2,ncol),stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          DO i=1,ncol
             CALL section_vals_val_get(reaction_path_section,"MAP%RANGE",r_vals=g_range,error=error)
             p_range(:,i)=g_range(:)
             p_bounds(2,i)=CEILING( MAX(p_range(1,i),p_range(2,i))/grid_sp(i))
             p_bounds(1,i)=FLOOR(MIN(p_range(1,i) , p_range(2,i))/grid_sp(i))
          END DO
          ALLOCATE(s1v(2,istart:iend),stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          ALLOCATE(s1(2),stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          ALLOCATE(grid_point(ncol),stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          v_count=0
          kk=rec_eval_grid(iw1,ncol,colvar%reaction_path_param%f_vals,v_count,&
               grid_point,grid_sp,colvar%reaction_path_param%step_size,istart,&
               iend,s1v,s1,p_bounds,colvar%reaction_path_param%lambda)
          DEALLOCATE(grid_sp,stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          DEALLOCATE(p_range,stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          DEALLOCATE(p_bounds,stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          DEALLOCATE(s1v,stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          DEALLOCATE(s1,stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          DEALLOCATE(grid_point,stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       END IF
       CALL cp_print_key_finished_output(iw1,logger,reaction_path_section,&
            "MAP", error=error)
    ELSE IF(my_subsection(11))THEN
       ! combine colvar
       CALL colvar_create(colvar,combine_colvar_id, error)
       colvar%use_points=.FALSE.
       colvar_subsection => section_vals_get_subs_vals(combine_section,"COLVAR",error=error)
       CALL section_vals_get(colvar_subsection,n_repetition=ncol,error=error)
       ALLOCATE(colvar%combine_cvs_param%colvar_p(ncol),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       ! In case we need to print some information..
       iw = cp_print_key_unit_nr(logger,colvar_section,&
            "PRINT%PROGRAM_RUN_INFO",extension=".colvarLog",error=error)
       IF (iw>0) THEN
          WRITE ( iw, '( A )')'          '//&
               '**********************************************************************'
          WRITE ( iw, '( A,I8)' )   ' COLVARS| COLVAR INPUT INDEX: ',icol
          WRITE ( iw, '( A,T49,4I8)' )   ' COLVARS| COMBINATION OF THE FOLOWING COLVARS:'
       END IF
       CALL cp_print_key_finished_output(iw,logger,colvar_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)
       ! Parsing the real COLVARs
       DO i= 1, ncol
          NULLIFY(colvar%combine_cvs_param%colvar_p(i)%colvar)
          CALL colvar_read(colvar%combine_cvs_param%colvar_p(i)%colvar,i,colvar_subsection, error=error)
       END DO
       ! Function definition
       CALL section_vals_val_get(combine_section,"FUNCTION",c_val=colvar%combine_cvs_param%function,error=error)
       CALL compress(colvar%combine_cvs_param%function, full=.TRUE.)
       ! Variables
       CALL section_vals_val_get(combine_section,"VARIABLES",c_vals=my_par,error=error)
       ALLOCATE(colvar%combine_cvs_param%variables(SIZE(my_par)),stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       colvar%combine_cvs_param%variables = my_par
       ! Check that the number of COLVAR provided is equal to the number of variables..
       CALL cp_assert(SIZE(my_par)==ncol,cp_fatal_level,cp_assertion_failed,routineP,&
            "Number of defined COLVAR for COMBINE_COLVAR is different from the "//&
            "number of variables! It is not possible to define COLVARs in a COMBINE_COLVAR "//&
            "and avoid their usage in the combininig function!"//&
CPSourceFileRef)
       ! Parameters
       ALLOCATE(colvar%combine_cvs_param%c_parameters(0),stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       CALL section_vals_val_get(combine_section,"PARAMETERS",n_rep_val=ncol,error=error)
       DO i = 1,ncol
          isize = SIZE(colvar%combine_cvs_param%c_parameters)
          CALL section_vals_val_get(combine_section,"PARAMETERS",c_vals=my_par,i_rep_val=i,error=error)
          CALL reallocate(colvar%combine_cvs_param%c_parameters,1,isize+SIZE(my_par))
          colvar%combine_cvs_param%c_parameters(isize+1:isize+SIZE(my_par)) = my_par
       END DO
       ALLOCATE(colvar%combine_cvs_param%v_parameters(0),stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       CALL section_vals_val_get(combine_section,"VALUES",n_rep_val=ncol,error=error)
       DO i = 1,ncol
          isize = SIZE(colvar%combine_cvs_param%v_parameters)
          CALL section_vals_val_get(combine_section,"VALUES",r_vals=my_vals,i_rep_val=i,error=error)
          CALL reallocate(colvar%combine_cvs_param%v_parameters,1,isize+SIZE(my_vals))
          colvar%combine_cvs_param%v_parameters(isize+1:isize+SIZE(my_vals)) = my_vals
       END DO
       ! Info on derivative evaluation
       CALL section_vals_val_get(combine_section,"DX",r_val=colvar%combine_cvs_param%dx,error=error)
       CALL section_vals_val_get(combine_section,"ERROR_LIMIT",r_val=colvar%combine_cvs_param%lerr,&
            error=error)
    END IF
    CALL colvar_setup(colvar, error)

    iw = cp_print_key_unit_nr(logger,colvar_section,&
         "PRINT%PROGRAM_RUN_INFO",extension=".colvarLog",error=error)
    IF (iw>0) THEN
       tag = "ATOMS: "
       IF (colvar%use_points) tag = "POINTS:"
       ! Description header
       IF (colvar%type_id/=combine_colvar_id) THEN
          WRITE ( iw, '( A )')'          '//&
               '----------------------------------------------------------------------'
          WRITE ( iw, '( A,I8)' )   ' COLVARS| COLVAR INPUT INDEX: ',icol
       END IF
       ! Colvar Description
       SELECT CASE(colvar%type_id)
       CASE(angle_colvar_id)
          WRITE ( iw, '( A,T57,3I8)' )   ' COLVARS| ANGLE          >>> '//tag,&
               colvar%angle_param%i_at_angle
       CASE(dfunct_colvar_id)
          WRITE ( iw, '( A,T49,4I8)' )   ' COLVARS| DISTANCE DIFFERENCE  >>> '//tag,&
               colvar%dfunct_param%i_at_dfunct
       CASE(plane_distance_colvar_id)
          WRITE ( iw, '( A,T57,3I8)' )   ' COLVARS| PLANE DISTANCE - PLANE  >>> '//tag,&
               colvar%plane_distance_param%plane
          WRITE ( iw, '( A,T73,1I8)' )   ' COLVARS| PLANE DISTANCE - POINT  >>> '//tag,&
               colvar%plane_distance_param%point
       CASE(torsion_colvar_id)
          WRITE ( iw, '( A,T49,4I8)' )   ' COLVARS| TORSION       >>> '//tag,&
               colvar%torsion_param%i_at_tors
       CASE(dist_colvar_id)
          WRITE ( iw, '( A,T65,2I8)' )   ' COLVARS| BOND          >>> '//tag,&
               colvar%dist_param%i_at,colvar%dist_param%j_at
       CASE(coord_colvar_id)
          IF(colvar%coord_param%do_chain) THEN
             WRITE ( iw, '( A)' ) ' COLVARS| COORDINATION CHAIN FC(from->to)*FC(to->to_B)>> '
          END IF
          IF (colvar%coord_param%use_kinds_from) THEN
             WRITE ( iw, '( A,T71,A10)' )  (' COLVARS| COORDINATION  >>> FROM KINDS',&
                  ADJUSTR(colvar%coord_param%c_kinds_from(kk)(1:10)),&
                  kk=1,SIZE(colvar%coord_param%c_kinds_from))
          ELSE
             WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COORDINATION  >>> FROM '//tag,&
                  colvar%coord_param%i_at_from(kk),&
                  kk=1,SIZE(colvar%coord_param%i_at_from))
          END IF
          IF (colvar%coord_param%use_kinds_to) THEN
             WRITE ( iw, '( A,T71,A10)' )  (' COLVARS| COORDINATION  >>>   TO KINDS',&
                  ADJUSTR(colvar%coord_param%c_kinds_to(kk)(1:10)),&
                  kk=1,SIZE(colvar%coord_param%c_kinds_to))
          ELSE
             WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COORDINATION  >>>   TO '//tag,&
                  colvar%coord_param%i_at_to(kk),&
                  kk=1,SIZE(colvar%coord_param%i_at_to))
          END IF
          WRITE ( iw, '( A,T71,F10.5)' ) ' COLVARS| R0',colvar%coord_param%r_0
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| NN',colvar%coord_param%nncrd
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| ND',colvar%coord_param%ndcrd
          IF(colvar%coord_param%do_chain) THEN
             IF (colvar%coord_param%use_kinds_to_b) THEN
                WRITE ( iw, '( A,T71,A10)' )  (' COLVARS| COORDINATION  >>>   TO KINDS B',&
                     ADJUSTR(colvar%coord_param%c_kinds_to_b(kk)(1:10)),&
                     kk=1,SIZE(colvar%coord_param%c_kinds_to_b))
             ELSE
                WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COORDINATION  >>>   TO '//tag//' B',&
                     colvar%coord_param%i_at_to_b(kk),&
                     kk=1,SIZE(colvar%coord_param%i_at_to_b))
             END IF
             WRITE ( iw, '( A,T71,F10.5)' ) ' COLVARS| R0 B',colvar%coord_param%r_0_b
             WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| NN B',colvar%coord_param%nncrd_b
             WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| ND B',colvar%coord_param%ndcrd_b
          END IF
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
       CASE(combine_colvar_id)
          WRITE ( iw, '( A)' ) ' COLVARS| COMBINING FUNCTION : '//&
               TRIM(colvar%combine_cvs_param%function)
          WRITE ( iw, '( A)', ADVANCE="NO") ' COLVARS| VARIABLES : '
          DO i = 1, SIZE(colvar%combine_cvs_param%variables)
             WRITE ( iw, '( A)', ADVANCE="NO")&
                  TRIM(colvar%combine_cvs_param%variables(i))//" "
          END DO
          WRITE ( iw, '(/)')
          WRITE ( iw, '( A)' ) ' COLVARS| DEFINED PARAMETERS [label]  [value]:'
          DO i = 1, SIZE(colvar%combine_cvs_param%c_parameters)
             WRITE ( iw, '( A,A7,F9.3)' ) '                            ',&
                  TRIM(colvar%combine_cvs_param%c_parameters(i)),colvar%combine_cvs_param%v_parameters(i)
          END DO
          WRITE ( iw, '( A,T71,G10.5)' ) ' COLVARS| ERROR ON DERIVATIVE EVALUATION',&
               colvar%combine_cvs_param%lerr
          WRITE ( iw, '( A,T71,G10.5)' ) ' COLVARS| DX',&
               colvar%combine_cvs_param%dx
       CASE(reaction_path_colvar_id)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Description header for REACTION_PATH COLVAR missing!!", &
               error=error, error_level=cp_warning_level)
       CASE(hydronium_colvar_id)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Description header for HYDRONIUM COLVAR missing!!", &
               error=error, error_level=cp_warning_level)
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
       ! Close the description layer
       IF (colvar%type_id/=combine_colvar_id) THEN
          WRITE ( iw, '( A )')'          '//&
               '----------------------------------------------------------------------'
       ELSE
          WRITE ( iw, '( A )')'          '//&
               '**********************************************************************'
       END IF
    END IF
    CALL cp_print_key_finished_output(iw,logger,colvar_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)

    CALL timestop(handle)
  END SUBROUTINE colvar_read

! *****************************************************************************
!> \brief evaluates the derivatives (dsdr) given and due to the given colvar
!>      variables in a molecular environment
!> \param colvar the collective variable to evaluate
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author Teodoro Laino
! *****************************************************************************
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
       ! Initialize the content of the derivative
       colvar%dsdr=0.0_dp
       SELECT CASE(colvar%type_id)
       CASE (dist_colvar_id)
          CALL dist_colvar(colvar,cell,particles=particles,error=error)
       CASE (coord_colvar_id)
          CALL coord_colvar(colvar,cell,particles=particles, error=error)
       CASE (torsion_colvar_id)
          CALL torsion_colvar(colvar, cell, particles=particles, error=error)
       CASE (angle_colvar_id)
          CALL angle_colvar(colvar, cell, particles=particles, error=error)
       CASE (dfunct_colvar_id)
          CALL dfunct_colvar(colvar, cell, particles=particles, error=error)
       CASE (plane_distance_colvar_id)
          CALL plane_distance_colvar(colvar,cell, particles=particles, error=error)
       CASE (rotation_colvar_id)
          CALL rotation_colvar(colvar,cell, particles=particles, error=error)
       CASE (qparm_colvar_id)
          CALL qparm_colvar(colvar,cell, particles=particles, error=error)
       CASE (hydronium_colvar_id)
          CALL hydronium_colvar(colvar,cell,particles=particles, error=error)
       CASE (reaction_path_colvar_id)
          CALL reaction_path_colvar(colvar,cell,particles=particles, error=error)
       CASE (combine_colvar_id)
          CALL combine_colvar(colvar,cell,particles=particles, error=error)
       CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       ! Check for fixed atom constraints
       IF (PRESENT(fixd_list)) CALL check_fixed_atom_cns_colv(fixd_list, colvar)
    END IF

  END SUBROUTINE colvar_eval_mol_f

! *****************************************************************************
!> \brief evaluates the derivatives (dsdr) given and due to the given colvar
!> \param colvar the collective variable to evaluate
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \note
!>      The torsion that generally is defined without the continuity problem 
!>      here (for free energy calculations) is defined only for (-pi,pi]
!> \author Alessandro Laio and fawzi
! *****************************************************************************
  SUBROUTINE colvar_eval_glob_f(icolvar,force_env,error)
    INTEGER                                  :: icolvar
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_eval_glob_f', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: colvar_ok, failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cp_subsys_type), POINTER            :: subsys

    failure=.FALSE.
    NULLIFY(subsys,cell,colvar)
    CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
    colvar_ok=ASSOCIATED(subsys%colvar_p)
    CPAssert(colvar_ok,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       colvar => subsys%colvar_p(icolvar)%colvar
       ! Initialize the content of the derivative
       colvar%dsdr=0.0_dp
       SELECT CASE(colvar%type_id)
       CASE (dist_colvar_id)
          CALL dist_colvar(colvar,cell,subsys=subsys,error=error)
       CASE (coord_colvar_id)
          CALL coord_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (torsion_colvar_id)
          CALL torsion_colvar(colvar,cell,subsys=subsys, no_riemann_sheet_op=.TRUE.,error=error)
       CASE (angle_colvar_id)
          CALL angle_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (dfunct_colvar_id)
          CALL dfunct_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (plane_distance_colvar_id)
          CALL plane_distance_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (rotation_colvar_id)
          CALL rotation_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (qparm_colvar_id)
          CALL qparm_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (hydronium_colvar_id)
          CALL hydronium_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (reaction_path_colvar_id)
          CALL reaction_path_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (combine_colvar_id)
          CALL combine_colvar(colvar,cell,subsys=subsys, error=error)
       CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       ! Check for fixed atom constraints
       CALL check_fixed_atom_cns_colv(subsys%gci%fixd_list,colvar)
    END IF
  END SUBROUTINE colvar_eval_glob_f

! *****************************************************************************
!> \brief evaluates the derivatives (dsdr) given and due to the given colvar
!>        for the specification of a recursive colvar type
!> \param colvar the collective variable to evaluate
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author sfchiff
! *****************************************************************************
  SUBROUTINE colvar_recursive_eval(colvar,cell,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_recursive_eval', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT. failure) THEN
       ! Initialize the content of the derivative
       colvar%dsdr=0.0_dp
       SELECT CASE(colvar%type_id)
       CASE (dist_colvar_id)
          CALL dist_colvar(colvar,cell,particles=particles,error=error)
       CASE (coord_colvar_id)
          CALL coord_colvar(colvar,cell,particles=particles, error=error)
       CASE (torsion_colvar_id)
          CALL torsion_colvar(colvar, cell, particles=particles, error=error)
       CASE (angle_colvar_id)
          CALL angle_colvar(colvar, cell, particles=particles, error=error)
       CASE (dfunct_colvar_id)
          CALL dfunct_colvar(colvar, cell, particles=particles, error=error)
       CASE (plane_distance_colvar_id)
          CALL plane_distance_colvar(colvar,cell, particles=particles, error=error)
       CASE (rotation_colvar_id)
          CALL rotation_colvar(colvar,cell, particles=particles, error=error)
       CASE (qparm_colvar_id)
          CALL qparm_colvar(colvar,cell, particles=particles, error=error)
       CASE (hydronium_colvar_id)
          CALL hydronium_colvar(colvar,cell,particles=particles, error=error)
       CASE (reaction_path_colvar_id)
          CALL reaction_path_colvar(colvar,cell,particles=particles, error=error)
       CASE (combine_colvar_id)
          CALL combine_colvar(colvar,cell,particles=particles, error=error)
       CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
    END IF
  END SUBROUTINE colvar_recursive_eval

! *****************************************************************************
!> \brief Get coordinates of atoms or of geoemtrical points
!> \author Teodoro Laino 03.2007 [created]
! *****************************************************************************
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

! *****************************************************************************
!> \brief Transfer derivatives to ds/dr
!> \author Teodoro Laino 03.2007 [created]
! *****************************************************************************
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

! *****************************************************************************
!> \brief evaluates the force due (and on) the distance from the plane collective variable
!> \author Teodoro Laino 02.2006 [created]
! *****************************************************************************
  SUBROUTINE plane_distance_colvar(colvar,cell,subsys,particles,error)

    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
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

! *****************************************************************************
!> \brief Evaluates the value of the rotation angle between two bonds 
!> \author Teodoro Laino 02.2006 [created]
! *****************************************************************************
  SUBROUTINE rotation_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
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

! *****************************************************************************
!> \brief evaluates the force due to the function of two distances
!> \author Teodoro Laino 02.2006 [created]
!        modified Florian Schiffmann 08.2008
! *****************************************************************************
  SUBROUTINE dfunct_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dfunct_colvar', &
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

    CPPrecondition(colvar%type_id==dfunct_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    i=colvar%dfunct_param%i_at_dfunct(1)
    j=colvar%dfunct_param%i_at_dfunct(2)
    ! First bond
    CALL get_coordinates(colvar, i, xpi, my_particles)
    CALL get_coordinates(colvar, j, xpj, my_particles)
    ss=MATMUL(cell%h_inv,xpi-xpj)
    ss=ss-NINT(ss)
    xij=MATMUL(cell%hmat,ss)
    r12=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
    ! Second bond
    k=colvar%dfunct_param%i_at_dfunct(3)
    l=colvar%dfunct_param%i_at_dfunct(4)
    CALL get_coordinates(colvar, k, xpk, my_particles)
    CALL get_coordinates(colvar, l, xpl, my_particles)
    ss=MATMUL(cell%h_inv,xpk-xpl)
    ss=ss-NINT(ss)
    xkl=MATMUL(cell%hmat,ss)
    r34=SQRT(xkl(1)**2+xkl(2)**2+xkl(3)**2)
    !
    colvar%ss=r12+colvar%dfunct_param%coeff*r34
    fi(:)= xij/r12
    fj(:)=-xij/r12
    fk(:)=colvar%dfunct_param%coeff*xkl/r34
    fl(:)=-colvar%dfunct_param%coeff*xkl/r34
    CALL put_derivative(colvar, 1, fi)
    CALL put_derivative(colvar, 2, fj)
    CALL put_derivative(colvar, 3, fk)
    CALL put_derivative(colvar, 4, fl)

  END SUBROUTINE dfunct_colvar

! *****************************************************************************
!> \brief evaluates the force due (and on) the distance from the plane collective variable
!> \author Teodoro Laino 02.2006 [created]
! *****************************************************************************
  SUBROUTINE angle_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
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

! *****************************************************************************
!> \brief evaluates the force due (and on) the distance collective variable
!> \author Alessandro Laio, Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE dist_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
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
    i=colvar%dist_param%i_at
    j=colvar%dist_param%j_at
    CALL get_coordinates(colvar, i, xpi, my_particles)
    CALL get_coordinates(colvar, j, xpj, my_particles)
    ss=MATMUL(cell%h_inv,xpi-xpj)
    ss=ss-NINT(ss)
    xij=MATMUL(cell%hmat,ss)
    SELECT CASE (colvar%dist_param%axis_id)
    CASE(do_clv_x)
       xij(2)=0.0_dp
       xij(3)=0.0_dp
    CASE(do_clv_y) 
       xij(1)=0.0_dp
       xij(3)=0.0_dp
    CASE(do_clv_z) 
       xij(1)=0.0_dp
       xij(2)=0.0_dp
    CASE(do_clv_xy) 
       xij(3)=0.0_dp
    CASE(do_clv_xz) 
       xij(2)=0.0_dp
    CASE(do_clv_yz) 
       xij(1)=0.0_dp
    CASE DEFAULT
       !do_clv_xyz
    END SELECT
    r12=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)

    colvar%ss=r12
    fi(:)=xij/r12
    fj(:)=-xij/r12

    CALL put_derivative(colvar, 1, fi)
    CALL put_derivative(colvar, 2, fj)

  END SUBROUTINE dist_colvar

! *****************************************************************************
!> \brief evaluates the force due to the torsion collective variable
!> \author Alessandro Laio, Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE torsion_colvar(colvar,cell,subsys,particles,no_riemann_sheet_op,error)

    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    LOGICAL, INTENT(IN), OPTIONAL            :: no_riemann_sheet_op
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'torsion_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii
    LOGICAL                                  :: failure, no_riemann_sheet
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
    no_riemann_sheet = .FALSE.
    IF (PRESENT(no_riemann_sheet_op)) no_riemann_sheet = no_riemann_sheet_op
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
    IF (no_riemann_sheet) colvar%ss = ATAN2(SIN(e),COS(e))
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

! *****************************************************************************
!> \brief evaluates the force due (and on) the Q PARM collective variable
! *****************************************************************************
  SUBROUTINE qparm_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
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

! *****************************************************************************
!> \brief evaluates the force due (and on) the hydronium collective variable
!> \author Marcel Baer
!> \note This function needs to be extended to the POINT structure!! 
!>       non-standard conform.. it's a breach in the colvar module.
! *****************************************************************************
  SUBROUTINE hydronium_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
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

! *****************************************************************************
!> \brief evaluates the force due (and on) the coordination-chain collective variable
!> \author MI
!> \note When the third set of atoms is not defined, this variable is equivalent
!>       to the simple coordination number. 
! *****************************************************************************
  SUBROUTINE coord_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coord_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, jj, k, kk, &
                                                n_atoms_from, n_atoms_to_a, &
                                                n_atoms_to_b, p_a, p_b, q_a, &
                                                q_b
    LOGICAL                                  :: failure
    REAL(dp) :: dfunc_ij, dfunc_jk, func_ij, func_jk, func_k, &
      inv_n_atoms_from, invden_ij, invden_jk, ncoord, num_ij, num_jk, r_0_a, &
      r_0_b, rdist_ij, rdist_jk, rij, rjk
    REAL(dp), DIMENSION(3)                   :: ftmp_i, ftmp_j, ftmp_k, ss, &
                                                xij, xjk, xpi, xpj, xpk
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    ! If we defined the coordination number with KINDS then we have still
    ! to fill few missing informations...
    NULLIFY(particles_i)
    CPPrecondition(colvar%type_id==coord_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    n_atoms_to_a=colvar%coord_param%n_atoms_to
    n_atoms_to_b=colvar%coord_param%n_atoms_to_b
    n_atoms_from=colvar%coord_param%n_atoms_from
    p_a=colvar%coord_param%nncrd
    q_a=colvar%coord_param%ndcrd
    r_0_a=colvar%coord_param%r_0
    p_b=colvar%coord_param%nncrd_b
    q_b=colvar%coord_param%ndcrd_b
    r_0_b=colvar%coord_param%r_0_b

    ncoord=0.0_dp
    inv_n_atoms_from=1.0_dp/REAL(n_atoms_from)
    DO ii=1,n_atoms_from
       i=colvar%coord_param%i_at_from(ii)
       CALL get_coordinates(colvar, i, xpi, my_particles)
       DO jj=1,n_atoms_to_a
          i=colvar%coord_param%i_at_to(jj)
          CALL get_coordinates(colvar, i, xpj, my_particles)       
          ss=MATMUL(cell%h_inv,xpi(:)-xpj(:))
          ss=ss-NINT(ss)
          xij=MATMUL(cell%hmat,ss)
          rij=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
          IF(rij < 1.0e-8_dp)CYCLE
          rdist_ij = rij/r_0_a
          num_ij=(1.0_dp-rdist_ij**p_a)
          invden_ij=1.0_dp/(1.0_dp-rdist_ij**q_a)
          func_ij=num_ij*invden_ij
          dfunc_ij= (- p_a*rdist_ij **(p_a-1) *invden_ij  &
               + num_ij*(invden_ij)**2 * q_a *rdist_ij **(q_a-1))/(rij*r_0_a)
          IF(n_atoms_to_b /=0) THEN
             func_k = 0.0_dp
             DO kk = 1,n_atoms_to_b
                k=colvar%coord_param%i_at_to_b(kk)
                CALL get_coordinates(colvar, k, xpk, my_particles)       
                ss=MATMUL(cell%h_inv,xpj(:)-xpk(:))
                ss=ss-NINT(ss)
                xjk=MATMUL(cell%hmat,ss)
                rjk=SQRT(xjk(1)**2+xjk(2)**2+xjk(3)**2)
                IF(rjk < 1.0e-8_dp)CYCLE
                rdist_jk = rjk/r_0_b
                num_jk=(1.0_dp-rdist_jk**p_b)
                invden_jk=1.0_dp/(1.0_dp-rdist_jk**q_b)
                func_jk=num_jk*invden_jk
                dfunc_jk= (- p_b*rdist_jk **(p_b-1) *invden_jk  &
                     + num_jk*(invden_jk)**2 * q_b *rdist_jk **(q_b-1))/(rjk*r_0_b)
                func_k = func_k + func_jk
                ftmp_k = - func_ij * dfunc_jk * xjk
                CALL put_derivative(colvar, n_atoms_from+n_atoms_to_a+kk, ftmp_k)

                ftmp_j = - dfunc_ij * xij * func_jk + func_ij * dfunc_jk * xjk
                CALL put_derivative(colvar, n_atoms_from+jj, ftmp_j)
             END DO
          ELSE
             func_k = 1.0_dp
             dfunc_jk = 0.0_dp
             ftmp_j = -dfunc_ij*xij
             CALL put_derivative(colvar, n_atoms_from+jj, ftmp_j)
          END IF
          ncoord=ncoord+func_ij*func_k
          ftmp_i = dfunc_ij*xij*func_k
          CALL put_derivative(colvar, ii, ftmp_i)
       ENDDO
    ENDDO
    colvar%ss=ncoord*inv_n_atoms_from
    colvar%dsdr(:,:) = colvar%dsdr(:,:)*inv_n_atoms_from
  END SUBROUTINE coord_colvar

! *****************************************************************************
!> \brief  evaluates function and forces due to a combination of COLVARs
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  SUBROUTINE combine_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'combine_colvar', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: def_error, this_error
    CHARACTER(LEN=default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: my_par
    INTEGER                                  :: i, ii, j, ncolv, ndim, stat
    LOGICAL                                  :: failure
    REAL(dp)                                 :: err
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: dss_vals, my_val, ss_vals
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: fi
    REAL(dp), DIMENSION(3)                   :: ss
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    CPPrecondition(colvar%type_id==combine_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF

    ncolv=SIZE(colvar%combine_cvs_param%colvar_p)
    ALLOCATE(ss_vals(ncolv),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(dss_vals(ncolv),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    
    ! Evaluate the individual COLVARs
    DO i=1,ncolv
       CALL colvar_recursive_eval(colvar%combine_cvs_param%colvar_p(i)%colvar,cell,my_particles,error)
       ss_vals(i)=colvar%combine_cvs_param%colvar_p(i)%colvar%ss
    ENDDO

    ! Evaluate the combination of the COLVARs
    CALL initf(1)
    ndim = SIZE(colvar%combine_cvs_param%c_parameters)+&
           SIZE(colvar%combine_cvs_param%variables)
    ALLOCATE(my_par(ndim),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure) 
    my_par(1:SIZE(colvar%combine_cvs_param%variables)) =colvar%combine_cvs_param%variables
    my_par(SIZE(colvar%combine_cvs_param%variables)+1:)=colvar%combine_cvs_param%c_parameters
    ALLOCATE(my_val(ndim),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure) 
    my_val(1:SIZE(colvar%combine_cvs_param%variables)) =ss_vals
    my_val(SIZE(colvar%combine_cvs_param%variables)+1:)=colvar%combine_cvs_param%v_parameters
    CALL parsef(1,TRIM(colvar%combine_cvs_param%function),my_par)
    colvar%ss=evalf(1,my_val)
    DO i = 1, ncolv
       dss_vals(i) = evalfd(1,i,my_val,colvar%combine_cvs_param%dx, err)
       IF ((ABS(err)>colvar%combine_cvs_param%lerr)) THEN
          WRITE(this_error,"(A,G12.6,A)")"(",err,")"
          WRITE(def_error,"(A,G12.6,A)")"(",colvar%combine_cvs_param%lerr,")"
          CALL compress(this_error,.TRUE.)
          CALL compress(def_error,.TRUE.)
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
               ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
               TRIM(def_error)//' . '//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
       END IF
    END DO
    DEALLOCATE(my_val,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(my_par,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    CALL finalizef()

    ! Evaluate forces
    ALLOCATE(fi(3,colvar%n_atom_s),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ii=0
    DO i=1,ncolv
       DO j=1,colvar%combine_cvs_param%colvar_p(i)%colvar%n_atom_s
          ii=ii+1
          fi(:,ii)=colvar%combine_cvs_param%colvar_p(i)%colvar%dsdr(:,j)*dss_vals(i)
       END DO
    END DO

    DO i=1,colvar%n_atom_s
       CALL put_derivative(colvar,i,fi(:,i))
    END DO

    DEALLOCATE(fi,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ss_vals,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)  
    DEALLOCATE(dss_vals,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)  
  END SUBROUTINE combine_colvar

! *****************************************************************************
!> \brief evaluates the force due (and on) reaction path collective variable
!> \author fschiff
!> \note the system is still able to move in the space spanned by the CV
!>       perpendicular to the path 
! *****************************************************************************
  SUBROUTINE reaction_path_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reaction_path_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iend, ii, istart, j, k, &
                                                ncolv, stat
    LOGICAL                                  :: failure
    REAL(dp)                                 :: lambda, step_size
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: s1, ss_vals
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: ds1, f_vals, fi, s1v
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ds1v
    REAL(dp), DIMENSION(3)                   :: ss
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.

    CPPrecondition(colvar%type_id==reaction_path_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF

    istart=colvar%reaction_path_param%function_bounds(1)
    iend=colvar%reaction_path_param%function_bounds(2)

    step_size=colvar%reaction_path_param%step_size
    ncolv=colvar%reaction_path_param%n_components 
    lambda=colvar%reaction_path_param%lambda
    ALLOCATE(f_vals(ncolv,istart:iend),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)    
    f_vals=colvar%reaction_path_param%f_vals
    ALLOCATE(ss_vals(ncolv),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)    

    DO i=1,ncolv
       CALL colvar_recursive_eval(colvar%reaction_path_param%colvar_p(i)%colvar,cell,my_particles,error)
       ss_vals(i)=colvar%reaction_path_param%colvar_p(i)%colvar%ss
    ENDDO

    ALLOCATE(s1v(2,istart:iend),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1v(ncolv,2,istart:iend),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(s1(2),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1(ncolv,2),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    DO k=istart,iend
       s1v(1,k)=REAL(k,kind=dp)*step_size*EXP(-lambda*DOT_PRODUCT(ss_vals(:)-f_vals(:,k),ss_vals(:)-f_vals(:,k)))
       s1v(2,k)=EXP(-lambda*DOT_PRODUCT(ss_vals(:)-f_vals(:,k),ss_vals(:)-f_vals(:,k)))
       DO j=1,ncolv
          ds1v(j,1,k) =REAL(k,kind=dp)*step_size*f_vals(j,k)*&
               EXP(-lambda*DOT_PRODUCT(ss_vals(:)-f_vals(:,k),ss_vals(:)-f_vals(:,k)))
          ds1v(j,2,k) =f_vals(j,k)*EXP(-lambda*DOT_PRODUCT(ss_vals(:)-f_vals(:,k),ss_vals(:)-f_vals(:,k)))
       END DO
    END DO
    DO i=1,2
       s1(i)=accurate_sum(s1v(i,:))
       DO j=1,ncolv
          ds1(j,i)=accurate_sum(ds1v(j,i,:))
       END DO
    END DO

    colvar%ss=s1(1)/s1(2)
    ALLOCATE(fi(3,colvar%n_atom_s),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    ii=0
    DO i=1,ncolv
       DO j=1,colvar%reaction_path_param%colvar_p(i)%colvar%n_atom_s
          ii=ii+1
          fi(:,ii)=colvar%reaction_path_param%colvar_p(i)%colvar%dsdr(:,j)*lambda*(ds1(i,1) &
               /s1(2)-colvar%ss*ds1(i,2)/s1(2))*2.0_dp
       END DO
    END DO

    DO i=1,colvar%n_atom_s
       CALL put_derivative(colvar,i,fi(:,i))
    END DO

    DEALLOCATE(fi,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(f_vals,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ss_vals,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)  
    DEALLOCATE(s1v,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ds1v,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)       
    DEALLOCATE(s1,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ds1,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
  END SUBROUTINE reaction_path_colvar

! *****************************************************************************
!> \brief used to print reaction_path function values on an arbitrary dimensional grid
!> \author fschiff
! *****************************************************************************
  RECURSIVE FUNCTION rec_eval_grid(iw1,ncol,f_vals,v_count,&
       gp,grid_sp,step_size,istart,iend,s1v,s1,p_bounds,lambda) RESULT(k)
    INTEGER                                  :: iw1, ncol
    REAL(dp), DIMENSION(:, :), POINTER       :: f_vals
    INTEGER                                  :: v_count
    REAL(dp), DIMENSION(:), POINTER          :: gp, grid_sp
    REAL(dp)                                 :: step_size
    INTEGER                                  :: istart, iend
    REAL(dp), DIMENSION(:, :), POINTER       :: s1v
    REAL(dp), DIMENSION(:), POINTER          :: s1
    INTEGER, DIMENSION(:, :), POINTER        :: p_bounds
    REAL(dp)                                 :: lambda
    INTEGER                                  :: k

    INTEGER                                  :: count1, i

    k=1
    IF(v_count.LT.ncol)THEN
       count1=v_count+1
       DO i=p_bounds(1,count1),p_bounds(2,count1)
          gp(count1)=REAL(i,KIND=dp)*grid_sp(count1)
          k=rec_eval_grid(iw1,ncol,f_vals,count1,gp,grid_sp,step_size,istart,iend,s1v,s1,p_bounds,lambda)
       END DO
    ELSE IF(v_count==ncol)THEN
       DO i=istart,iend
          s1v(1,i)=REAL(i,kind=dp)*step_size*EXP(-lambda*DOT_PRODUCT(gp(:)-f_vals(:,i),&
               gp(:)-f_vals(:,i)))
          s1v(2,i)=EXP(-lambda*DOT_PRODUCT(gp(:)-f_vals(:,i),gp(:)-f_vals(:,i)))
       END DO
       DO i=1,2
          s1(i)=accurate_sum(s1v(i,:))

       END DO
       WRITE(iw1, '(5F10.5)')gp(:),s1(1)/s1(2)
    END IF
  END FUNCTION rec_eval_grid

END MODULE colvar_methods
