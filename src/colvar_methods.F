!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE colvar_types,                    ONLY: &
       HBP_colvar_id, Wc_colvar_id, angle_colvar_id, colvar_check_points, &
       colvar_create, colvar_setup, colvar_type, combine_colvar_id, &
       coord_colvar_id, dfunct_colvar_id, dist_colvar_id, &
       distance_from_path_colvar_id, eval_point_der, eval_point_mass, &
       eval_point_pos, gyration_colvar_id, hydronium_colvar_id, &
       mindist_colvar_id, plane_distance_colvar_id, &
       plane_plane_angle_colvar_id, population_colvar_id, qparm_colvar_id, &
       reaction_path_colvar_id, ring_puckering_colvar_id, rmsd_colvar_id, &
       rotation_colvar_id, torsion_colvar_id, u_colvar_id, &
       xyz_diag_colvar_id, xyz_outerdiag_colvar_id
  USE constraint_fxd,                  ONLY: check_fixed_atom_cns_colv
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_p_type,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             use_mixed_force
  USE force_fields_util,               ONLY: get_generic_info
  USE fparser,                         ONLY: EvalErrType,&
                                             evalf,&
                                             evalfd,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE input_constants,                 ONLY: &
       do_clv_x, do_clv_xy, do_clv_xz, do_clv_y, do_clv_yz, do_clv_z, &
       plane_def_atoms, plane_def_vec, rmsd_all, rmsd_list, rmsd_weightlist
  USE input_cp2k_colvar,               ONLY: create_colvar_xyz_d_section,&
                                             create_colvar_xyz_od_section
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: section_get_keyword,&
                                             section_release,&
                                             section_type,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: fac,&
                                             maxfac,&
                                             pi,&
                                             twopi
  USE mathlib,                         ONLY: vector_product
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE mixed_energy_types,              ONLY: mixed_force_type
  USE mixed_environment_utils,         ONLY: get_subsys_map_index
  USE molecule_kind_types,             ONLY: fixd_constraint_type
  USE particle_list_types,             ONLY: particle_list_p_type,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE rmsd,                            ONLY: rmsd3
  USE spherical_harmonics,             ONLY: dlegendre,&
                                             legendre
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE wannier_states_types,            ONLY: wannier_centres_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colvar_methods'
  REAL(KIND=dp), PRIVATE, PARAMETER    :: tolerance_acos = 1.0E-5_dp

  PUBLIC :: colvar_read,&
            colvar_eval_glob_f,&
            colvar_eval_mol_f

CONTAINS

! *****************************************************************************
!> \brief reads a colvar from the input
!> \param colvar the place where to store what will be read
!> \param icol number of the current colvar (repetition in colvar_section)
!> \param colvar_section the colvar section
!> \param para_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param extended_lagrange if the extended lagrangian formalismus should
!>        be used UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      04.2004 created [alessandro laio and fawzi mohamed]
!> \author teo
! *****************************************************************************
  RECURSIVE SUBROUTINE colvar_read(colvar, icol, colvar_section, para_env, error)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER, INTENT(IN)                      :: icol
    TYPE(section_vals_type), POINTER         :: colvar_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_read', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3)                         :: fmid
    CHARACTER(LEN=7)                         :: tag, tag_comp, tag_comp1, &
                                                tag_comp2
    CHARACTER(LEN=default_path_length)       :: path_function
    CHARACTER(LEN=default_string_length)     :: tmpStr, tmpStr2
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: c_kinds, my_par
    INTEGER :: handle, i, iatm, icomponent, iend, ifunc, ii, isize, istart, &
      iw, iw1, j, k, kk, n_var, n_var_k, ncol, ndim, nr_frame, stat, v_count
    INTEGER, DIMENSION(:), POINTER           :: iatms
    INTEGER, DIMENSION(:, :), POINTER        :: p_bounds
    LOGICAL                                  :: check, failure, &
                                                use_mixed_energy
    LOGICAL, DIMENSION(23)                   :: my_subsection
    REAL(dp), DIMENSION(:), POINTER          :: s1, wei, weights
    REAL(dp), DIMENSION(:, :), POINTER       :: p_range, s1v
    REAL(KIND=dp), DIMENSION(1)              :: my_val
    REAL(KIND=dp), DIMENSION(:), POINTER     :: g_range, grid_point, grid_sp, &
                                                my_vals, range
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER :: angle_section, colvar_subsection, &
      combine_section, coordination_section, dfunct_section, &
      distance_from_path_section, distance_section, frame_section, &
      gyration_section, HBP_section, hydronium_section, mindist_section, &
      path_section, plane_dist_section, plane_plane_angle_section, &
      plane_sections, point_section, population_section, qparm_section, &
      reaction_path_section, ring_puckering_section, rmsd_section, &
      rotation_section, torsion_section, u_section, Wc_section, wrk_section, &
      xyz_diag_section, xyz_outerdiag_section

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(logger, c_kinds, iatms)
    logger => cp_error_get_logger(error)
    my_subsection        = .FALSE.
    distance_section      => section_vals_get_subs_vals(colvar_section,"DISTANCE",i_rep_section=icol,error=error)
    dfunct_section        => section_vals_get_subs_vals(colvar_section,"DISTANCE_FUNCTION",&
                             i_rep_section=icol,error=error)
    angle_section         => section_vals_get_subs_vals(colvar_section,"ANGLE",i_rep_section=icol,error=error)
    torsion_section       => section_vals_get_subs_vals(colvar_section,"TORSION",i_rep_section=icol,error=error)
    coordination_section  => section_vals_get_subs_vals(colvar_section,"COORDINATION",i_rep_section=icol,&
                             error=error)
    plane_dist_section    => section_vals_get_subs_vals(colvar_section,"DISTANCE_POINT_PLANE",i_rep_section=icol,&
                             error=error)
    plane_plane_angle_section &
                          => section_vals_get_subs_vals(colvar_section,"ANGLE_PLANE_PLANE",i_rep_section=icol,&
                             error=error)
    rotation_section      => section_vals_get_subs_vals(colvar_section,"BOND_ROTATION",i_rep_section=icol,&
                             error=error)
    qparm_section         => section_vals_get_subs_vals(colvar_section,"QPARM",i_rep_section=icol,error=error)
    hydronium_section     => section_vals_get_subs_vals(colvar_section,"HYDRONIUM",i_rep_section=icol,error=error)
    reaction_path_section => section_vals_get_subs_vals(colvar_section,"REACTION_PATH",i_rep_section=icol,&
                             can_return_null=.TRUE.,error=error)
    distance_from_path_section &
                          => section_vals_get_subs_vals(colvar_section,"DISTANCE_FROM_PATH",&
                             i_rep_section=icol, can_return_null=.TRUE.,error=error)
    combine_section       => section_vals_get_subs_vals(colvar_section,"COMBINE_COLVAR",i_rep_section=icol,&
                             can_return_null=.TRUE.,error=error)
    population_section    => section_vals_get_subs_vals(colvar_section,"POPULATION",i_rep_section=icol,error=error)
    gyration_section      => section_vals_get_subs_vals(colvar_section,"GYRATION_RADIUS",i_rep_section=icol,&
                             error=error)
    rmsd_section          => section_vals_get_subs_vals(colvar_section,"RMSD",i_rep_section=icol,error=error)
    xyz_diag_section      => section_vals_get_subs_vals(colvar_section,"XYZ_DIAG",i_rep_section=icol,error=error)
    xyz_outerdiag_section => section_vals_get_subs_vals(colvar_section,"XYZ_OUTERDIAG",i_rep_section=icol,&
                             error=error)
    u_section             => section_vals_get_subs_vals(colvar_section,"U",i_rep_section=icol,error=error)
    Wc_section            => section_vals_get_subs_vals(colvar_section,"WC",i_rep_section=icol,error=error)
    HBP_section           => section_vals_get_subs_vals(colvar_section,"HBP",i_rep_section=icol,error=error)
    ring_puckering_section&
                          => section_vals_get_subs_vals(colvar_section,"RING_PUCKERING",i_rep_section=icol,error=error)
    mindist_section       => section_vals_get_subs_vals(colvar_section,"CONDITIONED_DISTANCE",i_rep_section=icol,error=error)

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
       CALL section_vals_get(reaction_path_section,&
                                                explicit=my_subsection(10), error=error)
    END IF
    IF (ASSOCIATED(distance_from_path_section)) THEN
       CALL section_vals_get(distance_from_path_section,&
                                                explicit=my_subsection(16), error=error)
    END IF
    IF (ASSOCIATED(combine_section)) THEN
       CALL section_vals_get(combine_section,   explicit=my_subsection(11), error=error)
    END IF
    CALL section_vals_get(population_section,   explicit=my_subsection(12), error=error)
    CALL section_vals_get(plane_plane_angle_section,&
                                                explicit=my_subsection(13), error=error)
    CALL section_vals_get(gyration_section,     explicit=my_subsection(14), error=error)
    CALL section_vals_get(rmsd_section,         explicit=my_subsection(15), error=error)
    CALL section_vals_get(xyz_diag_section,     explicit=my_subsection(17), error=error)
    CALL section_vals_get(xyz_outerdiag_section,explicit=my_subsection(18), error=error)
    CALL section_vals_get(u_section,            explicit=my_subsection(19), error=error)
    CALL section_vals_get(Wc_section,           explicit=my_subsection(20), error=error)
    CALL section_vals_get(HBP_section,          explicit=my_subsection(21), error=error)
    CALL section_vals_get(ring_puckering_section,&
                                                explicit=my_subsection(22), error=error)
    CALL section_vals_get(mindist_section,      explicit=my_subsection(23), error=error)

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
       CALL section_vals_val_get(plane_dist_section,"PBC",l_val=colvar%plane_distance_param%use_pbc,error=error)
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
       CALL section_vals_val_get(dfunct_section,"PBC",l_val=colvar%dfunct_param%use_pbc,error=error)
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
       ! Hydronium
       CALL colvar_create(colvar,hydronium_colvar_id, error)
       NULLIFY(colvar%hydronium_param%i_oxygens)
       NULLIFY(colvar%hydronium_param%i_hydrogens)

       CALL section_vals_val_get(hydronium_section,"OXYGENS",n_rep_val=n_var,error=error)
       ndim = 0
       DO k = 1, n_var
         CALL section_vals_val_get(hydronium_section,"OXYGENS",i_vals=iatms,error=error)
         CALL reallocate(colvar%hydronium_param%i_oxygens,1,ndim+SIZE(iatms))
         colvar%hydronium_param%i_oxygens(ndim+1:ndim+SIZE(iatms)) = iatms
         ndim = ndim + SIZE(iatms)
       END DO
       colvar%hydronium_param%n_oxygens = ndim

       CALL section_vals_val_get(hydronium_section,"HYDROGENS",n_rep_val=n_var,error=error)
       ndim = 0
       DO k = 1, n_var
         CALL section_vals_val_get(hydronium_section,"HYDROGENS",i_vals=iatms,error=error)
         CALL reallocate(colvar%hydronium_param%i_hydrogens,1,ndim+SIZE(iatms))
         colvar%hydronium_param%i_hydrogens(ndim+1:ndim+SIZE(iatms)) = iatms
         ndim = ndim + SIZE(iatms)
       END DO
       colvar%hydronium_param%n_hydrogens = ndim

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
    ELSE IF(my_subsection(10) .OR. my_subsection(16))THEN
       !reaction path or distance from reaction path
       IF (my_subsection(10) ) THEN
          path_section => reaction_path_section
          CALL colvar_create(colvar,reaction_path_colvar_id, error)
          fmid="POS"
          ifunc=1
       ELSE IF (my_subsection(16) ) THEN
          path_section => distance_from_path_section
          CALL colvar_create(colvar,distance_from_path_colvar_id, error)
          fmid="DIS"
          ifunc=2
       END IF
       colvar%use_points=.FALSE.
       CALL section_vals_val_get(path_section,"LAMBDA",r_val=colvar%reaction_path_param%lambda,error=error)
       CALL section_vals_val_get(path_section,"DISTANCES_RMSD",l_val=colvar%reaction_path_param%dist_rmsd,&
            error=error)
       CALL section_vals_val_get(path_section,"RMSD",l_val=colvar%reaction_path_param%rmsd,&
            error=error)
       IF(colvar%reaction_path_param%dist_rmsd .AND. colvar%reaction_path_param%rmsd)THEN
          CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,routineP,&
                     "CV REACTION PATH: only one between DISTANCES_RMSD and RMSD can be used "//&
                     CPSourceFileRef)
       END IF
       IF(colvar%reaction_path_param%dist_rmsd .OR. colvar%reaction_path_param%rmsd) THEN
          NULLIFY(colvar%reaction_path_param%i_rmsd, colvar%reaction_path_param%r_ref)
          frame_section => section_vals_get_subs_vals(path_section,"FRAME",error=error)
          CALL section_vals_get(frame_section,n_repetition=nr_frame,error=error)

          colvar%reaction_path_param%nr_frames=nr_frame
          CALL read_frames(frame_section,para_env,nr_frame,colvar%reaction_path_param%r_ref,&
               colvar%reaction_path_param%n_components,error=error)
          CALL section_vals_val_get(path_section,"SUBSET_TYPE",i_val=colvar%reaction_path_param%subset,&
               error=error)
          IF (colvar%reaction_path_param%subset==rmsd_all) THEN
             ALLOCATE(colvar%reaction_path_param%i_rmsd(colvar%reaction_path_param%n_components),STAT=stat)
             CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
             DO i = 1, colvar%reaction_path_param%n_components
                colvar%reaction_path_param%i_rmsd(i) = i
             END DO
          ELSE IF (colvar%reaction_path_param%subset==rmsd_list) THEN
             ! This section can be repeated
             CALL section_vals_val_get(path_section,"ATOMS",n_rep_val=n_var,error=error)
             ndim = 0
             IF (n_var /= 0) THEN
                ! INDEX LIST
                DO k = 1, n_var
                   CALL section_vals_val_get(path_section,"ATOMS",i_rep_val=k,i_vals=iatms,error=error)
                   CALL reallocate(colvar%reaction_path_param%i_rmsd,1, ndim+SIZE(iatms))
                   colvar%reaction_path_param%i_rmsd(ndim+1:ndim+SIZE(iatms)) = iatms
                   ndim = ndim + SIZE(iatms)
                END DO
                colvar%reaction_path_param%n_components = ndim
             ELSE
                CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,routineP,&
                     "CV REACTION PATH: if SUBSET_TYPE=LIST a list of atoms needs to be provided "//&
                     CPSourceFileRef)
             END IF
          END IF

          CALL section_vals_val_get(path_section,"ALIGN_FRAMES",l_val=colvar%reaction_path_param%align_frames,&
               error=error)
       ELSE
          colvar_subsection => section_vals_get_subs_vals(path_section,"COLVAR",error=error)
          CALL section_vals_get(colvar_subsection,n_repetition=ncol,error=error)
          ALLOCATE(colvar%reaction_path_param%colvar_p(ncol),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (ncol>0) THEN
             DO i= 1, ncol
                NULLIFY(colvar%reaction_path_param%colvar_p(i)%colvar)
                CALL colvar_read(colvar%reaction_path_param%colvar_p(i)%colvar,i,colvar_subsection, para_env,&
                     error=error)
             ENDDO
          ELSE
             CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,routineP,&
                    "CV REACTION PATH: the number of CV to define the path must be >0 "//&
                    CPSourceFileRef)
          ENDIF
          colvar%reaction_path_param%n_components=ncol
          NULLIFY(range)
          CALL section_vals_val_get(path_section,"RANGE",r_vals=range,error=error)
          CALL section_vals_val_get(path_section,"STEP_SIZE",r_val=colvar%reaction_path_param%step_size,error=error)
          iend=CEILING(MAX(RANGE(1),RANGE(2))/colvar%reaction_path_param%step_size)
          istart=FLOOR(MIN(RANGE(1),RANGE(2))/colvar%reaction_path_param%step_size)
          colvar%reaction_path_param%function_bounds(1)=istart
          colvar%reaction_path_param%function_bounds(2)=iend
          colvar%reaction_path_param%nr_frames= 2 !iend - istart + 1
          ALLOCATE(colvar%reaction_path_param%f_vals(ncol,istart:iend),stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          CALL section_vals_val_get(path_section,"VARIABLE",c_vals=my_par,i_rep_val=1,error=error)
          CALL section_vals_val_get(path_section,"FUNCTION",n_rep_val=ncol,error=error)
          check = (ncol==SIZE(colvar%reaction_path_param%colvar_p))
          CPPostcondition(check,cp_failure_level,routinep,error,failure)
          CALL initf(ncol)
          DO i=1,ncol
             CALL section_vals_val_get(path_section,"FUNCTION",c_val=path_function,i_rep_val=i,error=error)
             CALL compress(path_function, full=.TRUE.)
             CALL parsef(i,TRIM(path_function),my_par)
             DO j=istart,iend
                my_val=REAL(j,kind=dp)*colvar%reaction_path_param%step_size
                colvar%reaction_path_param%f_vals(i,j)=evalf(i,my_val)
             END DO
          END DO
          CALL finalizef()

          iw1= cp_print_key_unit_nr(logger,path_section,&
               "MAP",middle_name=fmid,extension=".dat",file_status="REPLACE",error=error)
          IF(iw1>0)THEN
             CALL section_vals_val_get(path_section,"MAP%GRID_SPACING",n_rep_val=ncol,error=error)
             ALLOCATE(grid_sp(ncol),stat=stat)
             CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
             DO i=1,ncol
                CALL section_vals_val_get(path_section,"MAP%GRID_SPACING",r_val=grid_sp(i),error=error)
             END DO
             CALL section_vals_val_get(path_section,"MAP%RANGE",n_rep_val=ncol,error=error)
             CPPostcondition(ncol ==SIZE(grid_sp) ,cp_failure_level,routinep,error,failure)
             ALLOCATE(p_range(2,ncol),stat=stat)
             CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
             ALLOCATE(p_bounds(2,ncol),stat=stat)
             CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
             DO i=1,ncol
                CALL section_vals_val_get(path_section,"MAP%RANGE",r_vals=g_range,error=error)
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
                  iend,s1v,s1,p_bounds,colvar%reaction_path_param%lambda,ifunc=ifunc,&
                  nconf=colvar%reaction_path_param%nr_frames)
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
          CALL cp_print_key_finished_output(iw1,logger,path_section,&
               "MAP", error=error)
       END IF

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
          CALL colvar_read(colvar%combine_cvs_param%colvar_p(i)%colvar,i,colvar_subsection, para_env, error=error)
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
    ELSE IF (my_subsection(12)) THEN
       ! Population
       wrk_section => population_section
       CALL colvar_create(colvar, population_colvar_id, error)
       CALL colvar_check_points(colvar, population_section, error)

       NULLIFY(colvar%population_param%i_at_from, colvar%population_param%c_kinds_from)
       NULLIFY(colvar%population_param%i_at_to, colvar%population_param%c_kinds_to)
       ! This section can be repeated

       CALL section_vals_val_get(population_section,"ATOMS_FROM",n_rep_val=n_var,error=error)
       ndim = 0
       IF (n_var /= 0) THEN
          ! INDEX LIST
          DO k = 1, n_var
             CALL section_vals_val_get(population_section,"ATOMS_FROM",i_rep_val=k,i_vals=iatms,error=error)
             CALL reallocate(colvar%population_param%i_at_from,1, ndim+SIZE(iatms))
             colvar%population_param%i_at_from(ndim+1:ndim+SIZE(iatms)) = iatms
             ndim = ndim + SIZE(iatms)
          END DO
          colvar%population_param%n_atoms_from = ndim
          colvar%population_param%use_kinds_from = .FALSE.
       ELSE
          ! KINDS
          CALL section_vals_val_get(population_section,"KINDS_FROM",n_rep_val=n_var,error=error)
          CPPostcondition(n_var>0,cp_failure_level,routinep,error,failure)
          DO k = 1, n_var
             CALL section_vals_val_get(population_section,"KINDS_FROM",i_rep_val=k,c_vals=c_kinds,error=error)
             CALL reallocate(colvar%population_param%c_kinds_from,1, ndim+SIZE(c_kinds))
             colvar%population_param%c_kinds_from(ndim+1:ndim+SIZE(c_kinds)) = c_kinds
             ndim = ndim + SIZE(c_kinds)
          END DO
          colvar%population_param%n_atoms_from = 0
          colvar%population_param%use_kinds_from = .TRUE.
          ! Uppercase the label
          DO k = 1, ndim
             CALL uppercase(colvar%population_param%c_kinds_from(k))
          END DO
       END IF
       ! This section can be repeated
       CALL section_vals_val_get(population_section,"ATOMS_TO",n_rep_val=n_var,error=error)
       ndim = 0
       IF (n_var /= 0) THEN
          ! INDEX LIST
          DO k = 1, n_var
             CALL section_vals_val_get(population_section,"ATOMS_TO",i_rep_val=k,i_vals=iatms,error=error)
             CALL reallocate(colvar%population_param%i_at_to,1, ndim+SIZE(iatms))
             colvar%population_param%i_at_to(ndim+1:ndim+SIZE(iatms)) = iatms
             ndim = ndim + SIZE(iatms)
          END DO
          colvar%population_param%n_atoms_to = ndim
          colvar%population_param%use_kinds_to = .FALSE.
       ELSE
          ! KINDS
          CALL section_vals_val_get(population_section,"KINDS_TO",n_rep_val=n_var,error=error)
          CPPostcondition(n_var>0,cp_failure_level,routinep,error,failure)
          DO k = 1, n_var
             CALL section_vals_val_get(population_section,"KINDS_TO",i_rep_val=k,c_vals=c_kinds,error=error)
             CALL reallocate(colvar%population_param%c_kinds_to,1, ndim+SIZE(c_kinds))
             colvar%population_param%c_kinds_to(ndim+1:ndim+SIZE(c_kinds)) = c_kinds
             ndim = ndim + SIZE(c_kinds)
          END DO
          colvar%population_param%n_atoms_to = 0
          colvar%population_param%use_kinds_to = .TRUE.
          ! Uppercase the label
          DO k = 1, ndim
             CALL uppercase(colvar%population_param%c_kinds_to(k))
          END DO
       END IF
       ! Let's finish reading the other parameters
       CALL section_vals_val_get(population_section,"R0",r_val=colvar%population_param%r_0,error=error)
       CALL section_vals_val_get(population_section,"NN",i_val=colvar%population_param%nncrd,error=error)
       CALL section_vals_val_get(population_section,"ND",i_val=colvar%population_param%ndcrd,error=error)
       CALL section_vals_val_get(population_section,"N0",i_val=colvar%population_param%n0,error=error)
       CALL section_vals_val_get(population_section,"SIGMA",r_val=colvar%population_param%sigma,error=error)
    ELSE IF (my_subsection(13)) THEN
       ! Angle between two planes
       wrk_section => plane_plane_angle_section
       CALL colvar_create(colvar, plane_plane_angle_colvar_id, error)
       CALL colvar_check_points(colvar, plane_plane_angle_section, error)
       ! Read the specification of the two planes
       plane_sections => section_vals_get_subs_vals(plane_plane_angle_section,"PLANE",error=error)
       CALL section_vals_get(plane_sections, n_repetition=n_var, error=error)
       CALL cp_assert(n_var==2,cp_fatal_level,cp_assertion_failed,routineP,&
               "PLANE_PLANE_ANGLE Colvar section: Two PLANE sections must be provided!"//&
CPSourceFileRef)
       ! Plane 1
       CALL section_vals_val_get(plane_sections,"DEF_TYPE",i_rep_section=1,&
            i_val=colvar%plane_plane_angle_param%plane1%type_of_def,error=error)
       IF (colvar%plane_plane_angle_param%plane1%type_of_def==plane_def_vec) THEN
          CALL section_vals_val_get(plane_sections,"NORMAL_VECTOR",i_rep_section=1,&
               r_vals=s1,error=error)
          colvar%plane_plane_angle_param%plane1%normal_vec = s1
       ELSE
          CALL section_vals_val_get(plane_sections,"ATOMS",i_rep_section=1,&
               i_vals=iatms,error=error)
          colvar%plane_plane_angle_param%plane1%points = iatms
       END IF

       ! Plane 2
       CALL section_vals_val_get(plane_sections,"DEF_TYPE",i_rep_section=2,&
            i_val=colvar%plane_plane_angle_param%plane2%type_of_def,error=error)
       IF (colvar%plane_plane_angle_param%plane2%type_of_def==plane_def_vec) THEN
          CALL section_vals_val_get(plane_sections,"NORMAL_VECTOR",i_rep_section=2,&
               r_vals=s1,error=error)
          colvar%plane_plane_angle_param%plane2%normal_vec = s1
       ELSE
          CALL section_vals_val_get(plane_sections,"ATOMS",i_rep_section=2,&
               i_vals=iatms,error=error)
          colvar%plane_plane_angle_param%plane2%points = iatms
       END IF
    ELSE IF (my_subsection(14)) THEN
       ! Gyration Radius
       wrk_section => gyration_section
       CALL colvar_create(colvar, gyration_colvar_id, error)
       CALL colvar_check_points(colvar, gyration_section, error)

       NULLIFY(colvar%gyration_param%i_at, colvar%gyration_param%c_kinds)

       ! This section can be repeated
       CALL section_vals_val_get(gyration_section,"ATOMS",n_rep_val=n_var,error=error)
       ndim = 0
       IF (n_var /= 0) THEN
          ! INDEX LIST
          DO k = 1, n_var
             CALL section_vals_val_get(gyration_section,"ATOMS",i_rep_val=k,i_vals=iatms,error=error)
             CALL reallocate(colvar%gyration_param%i_at,1, ndim+SIZE(iatms))
             colvar%gyration_param%i_at(ndim+1:ndim+SIZE(iatms)) = iatms
             ndim = ndim + SIZE(iatms)
          END DO
          colvar%gyration_param%n_atoms = ndim
          colvar%gyration_param%use_kinds = .FALSE.
       ELSE
          ! KINDS
          CALL section_vals_val_get(gyration_section,"KINDS",n_rep_val=n_var,error=error)
          CPPostcondition(n_var>0,cp_failure_level,routinep,error,failure)
          DO k = 1, n_var
             CALL section_vals_val_get(gyration_section,"KINDS",i_rep_val=k,c_vals=c_kinds,error=error)
             CALL reallocate(colvar%gyration_param%c_kinds,1, ndim+SIZE(c_kinds))
             colvar%gyration_param%c_kinds(ndim+1:ndim+SIZE(c_kinds)) = c_kinds
             ndim = ndim + SIZE(c_kinds)
          END DO
          colvar%gyration_param%n_atoms = 0
          colvar%gyration_param%use_kinds = .TRUE.
          ! Uppercase the label
          DO k = 1, ndim
             CALL uppercase(colvar%gyration_param%c_kinds(k))
          END DO
       END IF
    ELSE IF (my_subsection(15)) THEN
       ! RMSD_AB
       wrk_section => rmsd_section
       CALL colvar_create(colvar, rmsd_colvar_id, error)

       NULLIFY(colvar%rmsd_param%i_rmsd, colvar%rmsd_param%r_ref, colvar%rmsd_param%weights)

       frame_section => section_vals_get_subs_vals(rmsd_section,"FRAME",error=error)
       CALL section_vals_get(frame_section,n_repetition=nr_frame,error=error)

       colvar%rmsd_param%nr_frames=nr_frame
       ! Calculation is aborted if reference frame are less than 2
       CPPostcondition((nr_frame>=1.AND.nr_frame<=2),cp_failure_level,routineP,error,failure)
       CALL read_frames(frame_section,para_env,nr_frame,colvar%rmsd_param%r_ref,&
              colvar%rmsd_param%n_atoms,error=error)

       ALLOCATE(colvar%rmsd_param%weights(colvar%rmsd_param%n_atoms), STAT=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       colvar%rmsd_param%weights = 0.0_dp

       CALL section_vals_val_get(rmsd_section,"SUBSET_TYPE",i_val=colvar%rmsd_param%subset,error=error)
       IF (colvar%rmsd_param%subset==rmsd_all) THEN
          ALLOCATE(colvar%rmsd_param%i_rmsd(colvar%rmsd_param%n_atoms),STAT=stat)
          CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          DO i = 1, colvar%rmsd_param%n_atoms
             colvar%rmsd_param%i_rmsd(i) = i
          END DO
       ELSE IF (colvar%rmsd_param%subset==rmsd_list) THEN
          ! This section can be repeated
          CALL section_vals_val_get(rmsd_section,"ATOMS",n_rep_val=n_var,error=error)
          ndim = 0
          IF (n_var /= 0) THEN
             ! INDEX LIST
             DO k = 1, n_var
                CALL section_vals_val_get(rmsd_section,"ATOMS",i_rep_val=k,i_vals=iatms,error=error)
                CALL reallocate(colvar%rmsd_param%i_rmsd,1, ndim+SIZE(iatms))
                colvar%rmsd_param%i_rmsd(ndim+1:ndim+SIZE(iatms)) = iatms
                ndim = ndim + SIZE(iatms)
             END DO
             colvar%rmsd_param%n_atoms = ndim
          ELSE
             CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,routineP,&
                  "CV RMSD: if SUBSET_TYPE=LIST a list of atoms needs to be provided "//&
                  CPSourceFileRef)
          END IF
       ELSE IF (colvar%rmsd_param%subset==rmsd_weightlist) THEN
          CALL section_vals_val_get(rmsd_section,"ATOMS",n_rep_val=n_var,error=error)
          ndim = 0
          IF (n_var /= 0) THEN
             ! INDEX LIST
             DO k = 1, n_var
                CALL section_vals_val_get(rmsd_section,"ATOMS",i_rep_val=k,i_vals=iatms,error=error)
                CALL reallocate(colvar%rmsd_param%i_rmsd,1, ndim+SIZE(iatms))
                colvar%rmsd_param%i_rmsd(ndim+1:ndim+SIZE(iatms)) = iatms
                ndim = ndim + SIZE(iatms)
             END DO
             colvar%rmsd_param%n_atoms = ndim
          ELSE
             CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,routineP,&
                  "CV RMSD: if SUBSET_TYPE=WEIGHT_LIST a list of atoms needs to be provided "//&
                  CPSourceFileRef)
          END IF
          CALL section_vals_val_get(rmsd_section,"WEIGHTS",n_rep_val=n_var,error=error)
          ndim = 0
          IF (n_var /= 0) THEN
             ! INDEX LIST
             DO k = 1, n_var
                CALL section_vals_val_get(rmsd_section,"WEIGHTS",i_rep_val=k,r_vals=wei,error=error)
                CALL reallocate(weights,1, ndim+SIZE(wei))
                weights(ndim+1:ndim+SIZE(wei)) = wei
                ndim = ndim + SIZE(wei)
             END DO
             CALL cp_assert(ndim==colvar%rmsd_param%n_atoms,cp_fatal_level,cp_assertion_failed,routineP,&
                  "CV RMSD: list of atoms and list of weights need to contain same number of entries. "//&
                  CPSourceFileRef)
             DO i = 1, ndim
                ii = colvar%rmsd_param%i_rmsd(i)
                colvar%rmsd_param%weights(ii) = weights(i)
             END DO
             DEALLOCATE (weights, STAT=stat)
             CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
          ELSE
             CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,routineP,&
                  "CV RMSD: if SUBSET_TYPE=WEIGHT_LIST a list of weights need to be provided. "//&
                  CPSourceFileRef)
          END IF

       ELSE
          CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,routineP,&
               "CV RMSD: unknown SUBSET_TYPE."//&
               CPSourceFileRef)
       END IF

       CALL section_vals_val_get(rmsd_section,"ALIGN_FRAMES",l_val=colvar%rmsd_param%align_frames,&
               error=error)
    ELSE IF (my_subsection(17)) THEN
       ! Work on XYZ positions of atoms
       wrk_section => xyz_diag_section
       CALL colvar_create(colvar, xyz_diag_colvar_id, error)
       CALL colvar_check_points(colvar, wrk_section, error)
       CALL section_vals_val_get(wrk_section,"ATOM",i_val=iatm,error=error)
       CALL section_vals_val_get(wrk_section,"COMPONENT",i_val=icomponent,error=error)
       CALL section_vals_val_get(wrk_section,"PBC",l_val=colvar%xyz_diag_param%use_pbc,error=error)
       CALL section_vals_val_get(wrk_section,"ABSOLUTE_POSITION",l_val=colvar%xyz_diag_param%use_absolute_position,error=error)
       colvar%xyz_diag_param%i_atom     = iatm
       colvar%xyz_diag_param%component  = icomponent
    ELSE IF (my_subsection(18)) THEN
       ! Work on the outer diagonal (two atoms A,B) XYZ positions
       wrk_section => xyz_outerdiag_section
       CALL colvar_create(colvar, xyz_outerdiag_colvar_id, error)
       CALL colvar_check_points(colvar, wrk_section, error)
       CALL section_vals_val_get(wrk_section,"ATOMS",i_vals=iatms,error=error)
       colvar%xyz_outerdiag_param%i_atoms       = iatms
       CALL section_vals_val_get(wrk_section,"COMPONENT_A",i_val=icomponent,error=error)
       colvar%xyz_outerdiag_param%components(1) = icomponent
       CALL section_vals_val_get(wrk_section,"COMPONENT_B",i_val=icomponent,error=error)
       colvar%xyz_outerdiag_param%components(2) = icomponent
       CALL section_vals_val_get(wrk_section,"PBC",l_val=colvar%xyz_outerdiag_param%use_pbc,error=error)
    ELSE IF (my_subsection(19)) THEN
       ! Energy
       wrk_section => u_section
       CALL colvar_create(colvar, u_colvar_id, error)
       colvar%u_param%mixed_energy_section => section_vals_get_subs_vals(wrk_section,"MIXED",error=error)
       CALL section_vals_get(colvar%u_param%mixed_energy_section,explicit=use_mixed_energy,error=error)
       IF (.NOT.use_mixed_energy) NULLIFY(colvar%u_param%mixed_energy_section)
    ELSE IF (my_subsection(20)) THEN
       ! Wc hydrogen bond
       wrk_section => Wc_section
       CALL colvar_create(colvar, Wc_colvar_id, error)
       CALL colvar_check_points(colvar, Wc_section, error)
       CALL section_vals_val_get(Wc_section,"ATOMS",i_vals=iatms,error=error)
       CALL section_vals_val_get(wrk_section,"RCUT",r_val=my_val(1),error=error)
       colvar%Wc%rcut = cp_unit_to_cp2k(my_val(1),"angstrom",error=error)
       colvar%Wc%ids = iatms
    ELSE IF (my_subsection(21)) THEN
       ! HBP colvar
       wrk_section => HBP_section
       CALL colvar_create(colvar, HBP_colvar_id, error)
       CALL colvar_check_points(colvar, HBP_section, error)
       CALL section_vals_val_get(wrk_section,"NPOINTS",i_val=colvar%HBP%nPoints,error=error)
       CALL section_vals_val_get(wrk_section,"RCUT",r_val=my_val(1),error=error)
       colvar%HBP%rcut = cp_unit_to_cp2k(my_val(1),"angstrom",error=error)
       CALL section_vals_val_get(wrk_section,"RCUT",r_val=colvar%HBP%shift,error=error)

       ALLOCATE(colvar%HBP%ids(colvar%HBP%nPoints,3),STAT=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       ALLOCATE(colvar%HBP%ewc(colvar%HBP%nPoints),STAT=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       DO i = 1, colvar%HBP%nPoints
        CALL section_vals_val_get(wrk_section,"ATOMS",i_rep_val=i,i_vals=iatms,error=error)
        colvar%HBP%ids(i,:) = iatms
       ENDDO
    ELSE IF (my_subsection(22)) THEN
       ! Ring Puckering
       CALL colvar_create(colvar,ring_puckering_colvar_id, error)
       CALL section_vals_val_get(ring_puckering_section,"ATOMS",i_vals=iatms,error=error)
       colvar%ring_puckering_param%nring = SIZE(iatms)
       ALLOCATE(colvar%ring_puckering_param%atoms(SIZE(iatms)),stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
       colvar%ring_puckering_param%atoms = iatms
       CALL section_vals_val_get(ring_puckering_section,"COORDINATE",&
            i_val=colvar%ring_puckering_param%iq,error=error)
       ! test the validity of the parameters
       ndim = colvar%ring_puckering_param%nring
       CALL cp_assert(ndim > 3,cp_fatal_level,cp_assertion_failed,routineP,&
            "CV Ring Puckering: Ring size has to be 4 or larger. "//&
            CPSourceFileRef)
       ii = colvar%ring_puckering_param%iq
       CALL cp_assert( (ABS(ii)/=1 .AND. ii>=-(ndim-1)/2 .AND. ii<=ndim/2) ,&
            cp_fatal_level,cp_assertion_failed,routineP,&
            "CV Ring Puckering: Invalid coordinate number. "//&
            CPSourceFileRef)
    ELSE IF (my_subsection(23)) THEN
       ! Minimum Distance
       wrk_section => mindist_section
       CALL colvar_create(colvar, mindist_colvar_id, error)
       CALL colvar_check_points(colvar, mindist_section, error)
       NULLIFY(colvar%mindist_param%i_dist_from,colvar%mindist_param%i_coord_from,&
               colvar%mindist_param%k_coord_from,colvar%mindist_param%i_coord_to,&
               colvar%mindist_param%k_coord_to)
       CALL section_vals_val_get(mindist_section,"ATOMS_DISTANCE",i_vals=iatms,error=error)
       colvar%mindist_param%n_dist_from = SIZE(iatms)
       ALLOCATE(colvar%mindist_param%i_dist_from(SIZE(iatms)),stat=stat)
       colvar%mindist_param%i_dist_from = iatms
       CALL section_vals_val_get(mindist_section,"ATOMS_FROM",n_rep_val=n_var,error=error)
       ndim = 0 
       IF (n_var /= 0) THEN
         ! INDEX LIST
         DO k = 1, n_var
           CALL section_vals_val_get(mindist_section,"ATOMS_FROM",i_rep_val=k,i_vals=iatms,error=error)
           CALL reallocate(colvar%mindist_param%i_coord_from,1, ndim+SIZE(iatms))
           colvar%mindist_param%i_coord_from(ndim+1:ndim+SIZE(iatms)) = iatms
           ndim = ndim + SIZE(iatms)
         END DO
         colvar%mindist_param%n_coord_from = ndim
         colvar%mindist_param%use_kinds_from = .FALSE.
       ELSE
         !KINDS
          CALL section_vals_val_get(mindist_section,"KINDS_FROM",n_rep_val=n_var,error=error)
          CPPostcondition(n_var>0,cp_failure_level,routinep,error,failure)
          DO k = 1, n_var
             CALL section_vals_val_get(mindist_section,"KINDS_FROM",i_rep_val=k,c_vals=c_kinds,error=error)
             CALL reallocate(colvar%mindist_param%k_coord_from,1, ndim+SIZE(c_kinds))
             colvar%mindist_param%k_coord_from(ndim+1:ndim+SIZE(c_kinds)) = c_kinds
             ndim = ndim + SIZE(c_kinds)
          END DO
          colvar%mindist_param%n_coord_from = 0
          colvar%mindist_param%use_kinds_from = .TRUE.
          ! Uppercase the label
          DO k = 1, ndim
             CALL uppercase(colvar%mindist_param%k_coord_from(k))
          END DO
       END IF


       CALL section_vals_val_get(mindist_section,"ATOMS_TO",n_rep_val=n_var,error=error)
       ndim = 0
       IF (n_var /= 0) THEN
          ! INDEX LIST
         DO k = 1, n_var
           CALL section_vals_val_get(mindist_section,"ATOMS_TO",i_rep_val=k,i_vals=iatms,error=error)
           CALL reallocate(colvar%mindist_param%i_coord_to,1, ndim+SIZE(iatms))
           colvar%mindist_param%i_coord_to(ndim+1:ndim+SIZE(iatms)) = iatms
           ndim = ndim + SIZE(iatms)
         END DO
         colvar%mindist_param%n_coord_to = ndim
         colvar%mindist_param%use_kinds_to = .FALSE.
       ELSE
         !KINDS
          CALL section_vals_val_get(mindist_section,"KINDS_TO",n_rep_val=n_var,error=error)
          CPPostcondition(n_var>0,cp_failure_level,routinep,error,failure)
          DO k = 1, n_var
             CALL section_vals_val_get(mindist_section,"KINDS_TO",i_rep_val=k,c_vals=c_kinds,error=error)
             CALL reallocate(colvar%mindist_param%k_coord_to,1, ndim+SIZE(c_kinds))
             colvar%mindist_param%k_coord_to(ndim+1:ndim+SIZE(c_kinds)) = c_kinds
             ndim = ndim + SIZE(c_kinds)
          END DO
          colvar%mindist_param%n_coord_to = 0
          colvar%mindist_param%use_kinds_to = .TRUE.
          ! Uppercase the label
          DO k = 1, ndim
             CALL uppercase(colvar%mindist_param%k_coord_to(k))
          END DO
       END IF

       CALL section_vals_val_get(mindist_section,"R0",r_val=colvar%mindist_param%r_cut,error=error)
       CALL section_vals_val_get(mindist_section,"NN",i_val=colvar%mindist_param%p_exp,error=error)
       CALL section_vals_val_get(mindist_section,"ND",i_val=colvar%mindist_param%q_exp,error=error)
!       CALL section_vals_val_get(mindist_section,"NC",r_val=colvar%mindist_param%n_cut,error=error)
       CALL section_vals_val_get(mindist_section,"LAMBDA",r_val=colvar%mindist_param%lambda,error=error)
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
       CASE(plane_plane_angle_colvar_id)
          IF (colvar%plane_plane_angle_param%plane1%type_of_def==plane_def_atoms) THEN
             WRITE ( iw, '( A,T57,3I8)' )   ' COLVARS| PLANE-PLANE ANGLE - PLANE 1  (ATOMS) >>> '//tag,&
                  colvar%plane_plane_angle_param%plane1%points
          ELSE
             WRITE ( iw, '( A,T57,3F8.3)' )   ' COLVARS| PLANE-PLANE ANGLE - PLANE 1 (VECTOR) >>> '//tag,&
                  colvar%plane_plane_angle_param%plane1%normal_vec
          END IF

          IF (colvar%plane_plane_angle_param%plane2%type_of_def==plane_def_atoms) THEN
             WRITE ( iw, '( A,T57,3I8)' )   ' COLVARS| PLANE-PLANE ANGLE - PLANE 1  (ATOMS) >>> '//tag,&
                  colvar%plane_plane_angle_param%plane2%points
          ELSE
             WRITE ( iw, '( A,T57,3F8.3)' )   ' COLVARS| PLANE-PLANE ANGLE - PLANE 1 (VECTOR) >>> '//tag,&
                  colvar%plane_plane_angle_param%plane2%normal_vec
          END IF
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
       CASE(population_colvar_id)
          IF (colvar%population_param%use_kinds_from) THEN
             WRITE ( iw, '( A,T71,A10)' )  (' COLVARS| POPULATION based on coordination >>> FROM KINDS',&
                  ADJUSTR(colvar%population_param%c_kinds_from(kk)(1:10)),&
                  kk=1,SIZE(colvar%population_param%c_kinds_from))
          ELSE
             WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| POPULATION based on coordination >>> FROM '//tag,&
                  colvar%population_param%i_at_from(kk),&
                  kk=1,SIZE(colvar%population_param%i_at_from))
          END IF
          IF (colvar%population_param%use_kinds_to) THEN
             WRITE ( iw, '( A,T71,A10)' )  (' COLVARS| POPULATION based on coordination >>>   TO KINDS',&
                  ADJUSTR(colvar%population_param%c_kinds_to(kk)(1:10)),&
                  kk=1,SIZE(colvar%population_param%c_kinds_to))
          ELSE
             WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| POPULATION based on coordination   >>>   TO '//tag,&
                  colvar%population_param%i_at_to(kk),&
                  kk=1,SIZE(colvar%population_param%i_at_to))
          END IF
          WRITE ( iw, '( A,T71,F10.5)' ) ' COLVARS| R0',colvar%population_param%r_0
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| NN',colvar%population_param%nncrd
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| ND',colvar%population_param%ndcrd
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| N0',colvar%population_param%n0
          WRITE ( iw, '( A,T71,F10.5)' )   ' COLVARS| SIGMA',colvar%population_param%sigma
       CASE(gyration_colvar_id)
          IF (colvar%gyration_param%use_kinds) THEN
             WRITE ( iw, '( A,T71,A10)' )  (' COLVARS| Gyration Radius >>> KINDS',&
                  ADJUSTR(colvar%gyration_param%c_kinds(kk)(1:10)),&
                  kk=1,SIZE(colvar%gyration_param%c_kinds))
          ELSE
             WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| Gyration Radius >>> ATOMS '//tag,&
                  colvar%gyration_param%i_at(kk),&
                  kk=1,SIZE(colvar%gyration_param%i_at))
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
       CASE(distance_from_path_colvar_id)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Description header for REACTION_PATH COLVAR missing!!", &
               error=error, error_level=cp_warning_level)
       CASE(hydronium_colvar_id)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Description header for HYDRONIUM COLVAR missing!!", &
               error=error, error_level=cp_warning_level)
       CASE(rmsd_colvar_id)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Description header for RMSD COLVAR missing!!", &
               error=error, error_level=cp_warning_level)
       CASE(xyz_diag_colvar_id)
          NULLIFY(section, keyword, enum)
          CALL create_colvar_xyz_d_section(section,error=error)
          keyword => section_get_keyword(section,"COMPONENT",error=error)
          CALL keyword_get(keyword,enum=enum,error=error)
          tag_comp = enum_i2c(enum,colvar%xyz_diag_param%component,error=error)
          CALL section_release(section,error=error)

          WRITE ( iw, '( A,T57,3I8)' )   ' COLVARS| POSITION ('//TRIM(tag_comp)&
               //')  >>> '//tag,colvar%xyz_diag_param%i_atom
       CASE(xyz_outerdiag_colvar_id)
          NULLIFY(section, keyword, enum)
          CALL create_colvar_xyz_od_section(section,error=error)
          keyword => section_get_keyword(section,"COMPONENT_A",error=error)
          CALL keyword_get(keyword,enum=enum,error=error)
          tag_comp1 = enum_i2c(enum,colvar%xyz_outerdiag_param%components(1),error=error)
          keyword => section_get_keyword(section,"COMPONENT_B",error=error)
          CALL keyword_get(keyword,enum=enum,error=error)
          tag_comp2 = enum_i2c(enum,colvar%xyz_outerdiag_param%components(2),error=error)
          CALL section_release(section,error=error)

          WRITE ( iw, '( A,T57,3I8)' )   ' COLVARS| CROSS TERM POSITION ('//TRIM(tag_comp1)&
               //" * "//TRIM(tag_comp2)//')  >>> '//tag, colvar%xyz_outerdiag_param%i_atoms
       CASE(u_colvar_id)
          WRITE ( iw, '( A,T77,A4)' ) ' COLVARS| ENERGY          >>> '//tag,'all!'
       CASE(Wc_colvar_id)
          WRITE ( iw, '( A,T57,F16.8)' )   ' COLVARS| Wc          >>> RCUT: ',&
               colvar%Wc%rcut
          WRITE ( iw, '( A,T57,3I8)' )   ' COLVARS| Wc          >>> '//tag,&
               colvar%Wc%ids
       CASE(HBP_colvar_id)
          WRITE ( iw, '( A,T57,I8)' )   ' COLVARS| HBP          >>> NPOINTS',&
               colvar%HBP%nPoints
          WRITE ( iw, '( A,T57,F16.8)' )   ' COLVARS| HBP          >>> RCUT',&
               colvar%HBP%rcut
          WRITE ( iw, '( A,T57,F16.8)' )   ' COLVARS| HBP          >>> RCUT',&
               colvar%HBP%shift
          DO i =1, colvar%HBP%nPoints
             WRITE ( iw, '( A,T57,3I8)' )   ' COLVARS| HBP          >>> '//tag,&
               colvar%HBP%ids(i,:)
          ENDDO
       CASE(ring_puckering_colvar_id)
          WRITE ( iw, '( A,T57,I8)' )   ' COLVARS| Ring Puckering      >>> ring size',&
               colvar%ring_puckering_param%nring
          IF (colvar%ring_puckering_param%iq == 0) THEN
             WRITE ( iw, '( A,T40,A)' )   ' COLVARS| Ring Puckering      >>> coordinate',&
                                          ' Total Puckering Amplitude'
          ELSEIF (colvar%ring_puckering_param%iq > 0) THEN
             WRITE ( iw, '( A,T35,A,T57,I8)' )   ' COLVARS| Ring Puckering      >>> coordinate',&
                                          ' Puckering Amplitude',&
               colvar%ring_puckering_param%iq
          ELSE
             WRITE ( iw, '( A,T35,A,T57,I8)' )   ' COLVARS| Ring Puckering      >>> coordinate',&
                                          ' Puckering Angle',&
               colvar%ring_puckering_param%iq
          END IF
       CASE(mindist_colvar_id)
          WRITE ( iw, '( A)' ) ' COLVARS| CONDITIONED DISTANCE>> '
          WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COND.DISTANCE  >>> DISTANCE FROM '//tag,&
                colvar%mindist_param%i_dist_from(kk),&
                kk=1,SIZE(colvar%mindist_param%i_dist_from))
          IF (colvar%mindist_param%use_kinds_from) THEN
             WRITE ( iw, '( A,T71,A10)' )  (' COLVARS| COND.DIST.  >>> COORDINATION FROM KINDS ',&
                  ADJUSTR(colvar%mindist_param%k_coord_from(kk)(1:10)),&
                  kk=1,SIZE(colvar%mindist_param%k_coord_from))
          ELSE
             WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COND.DIST.  >>> COORDINATION FROM '//tag,&
                  colvar%mindist_param%i_coord_from(kk),&
                  kk=1,SIZE(colvar%mindist_param%i_coord_from))
          END IF
          IF (colvar%mindist_param%use_kinds_to) THEN
             WRITE ( iw, '( A,T71,A10)' )  (' COLVARS| COND.DIST.  >>> COORDINATION TO KINDS ',&
                  ADJUSTR(colvar%mindist_param%k_coord_to(kk)(1:10)),&
                  kk=1,SIZE(colvar%mindist_param%k_coord_to))
          ELSE
             WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COND.DIST.  >>> COORDINATION TO '//tag,&
                  colvar%mindist_param%i_coord_to(kk),&
                  kk=1,SIZE(colvar%mindist_param%i_coord_to))
          END IF
          WRITE ( iw, '( A,T71,F10.5)' ) ' COLVARS| R0',colvar%mindist_param%r_cut
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| NN',colvar%mindist_param%p_exp
          WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| ND',colvar%mindist_param%q_exp
          WRITE ( iw, '( A,T71,F10.5)' ) ' COLVARS| LAMBDA',colvar%mindist_param%lambda

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
!> \param cell ...
!> \param particles ...
!> \param pos ...
!> \param fixd_list ...
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
       CASE (population_colvar_id)
          CALL population_colvar(colvar,cell,particles=particles, error=error)
       CASE (gyration_colvar_id)
          CALL gyration_radius_colvar(colvar,cell,particles=particles, error=error)
       CASE (torsion_colvar_id)
          CALL torsion_colvar(colvar, cell, particles=particles, error=error)
       CASE (angle_colvar_id)
          CALL angle_colvar(colvar, cell, particles=particles, error=error)
       CASE (dfunct_colvar_id)
          CALL dfunct_colvar(colvar, cell, particles=particles, error=error)
       CASE (plane_distance_colvar_id)
          CALL plane_distance_colvar(colvar,cell, particles=particles, error=error)
       CASE (plane_plane_angle_colvar_id)
          CALL plane_plane_angle_colvar(colvar,cell, particles=particles, error=error)
       CASE (rotation_colvar_id)
          CALL rotation_colvar(colvar,cell, particles=particles, error=error)
       CASE (qparm_colvar_id)
          CALL qparm_colvar(colvar,cell, particles=particles, error=error)
       CASE (hydronium_colvar_id)
          CALL hydronium_colvar(colvar,cell,particles=particles, error=error)
       CASE(rmsd_colvar_id)
          CALL rmsd_colvar(colvar,particles=particles,error=error)
       CASE (reaction_path_colvar_id)
          CALL reaction_path_colvar(colvar,cell,particles=particles, error=error)
       CASE (distance_from_path_colvar_id)
          CALL distance_from_path_colvar(colvar,cell,particles=particles, error=error)
       CASE (combine_colvar_id)
          CALL combine_colvar(colvar,cell,particles=particles, error=error)
       CASE (xyz_diag_colvar_id)
          CALL xyz_diag_colvar(colvar,cell,particles=particles, error=error)
       CASE (xyz_outerdiag_colvar_id)
          CALL xyz_outerdiag_colvar(colvar,cell,particles=particles, error=error)
       CASE (ring_puckering_colvar_id)
          CALL ring_puckering_colvar(colvar,cell,particles=particles, error=error)
       CASE (mindist_colvar_id)
          CALL mindist_colvar(colvar,cell,particles=particles, error=error)
       CASE (u_colvar_id)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="need force_env!", error=error, error_level=cp_failure_level)
       CASE (Wc_colvar_id)
          !!! FIXME this is rubbish at the moment as we have no force to be computed on this
          CALL Wc_colvar(colvar, cell, particles=particles, error=error)
       CASE (HBP_colvar_id)
          !!! FIXME this is rubbish at the moment as we have no force to be computed on this
          CALL HBP_colvar(colvar, cell, particles=particles, error=error)
       CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       ! Check for fixed atom constraints
       IF (PRESENT(fixd_list)) CALL check_fixed_atom_cns_colv(fixd_list, colvar)
    END IF

  END SUBROUTINE colvar_eval_mol_f

! *****************************************************************************
!> \brief evaluates the derivatives (dsdr) given and due to the given colvar
!> \param icolvar ...
!> \param force_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param colvar the collective variable to evaluate UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author Alessandro Laio and fawzi
!> \note
!>      The torsion that generally is defined without the continuity problem
!>      here (for free energy calculations) is defined only for (-pi,pi]
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
    TYPE(qs_environment_type), POINTER       :: qs_env

    failure=.FALSE.
    NULLIFY(subsys,cell,colvar,qs_env)
    CALL force_env_get(force_env,subsys=subsys,cell=cell,qs_env=qs_env,error=error)
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
       CASE (population_colvar_id)
          CALL population_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (gyration_colvar_id)
          CALL gyration_radius_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (torsion_colvar_id)
          CALL torsion_colvar(colvar,cell,subsys=subsys, no_riemann_sheet_op=.TRUE.,error=error)
       CASE (angle_colvar_id)
          CALL angle_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (dfunct_colvar_id)
          CALL dfunct_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (plane_distance_colvar_id)
          CALL plane_distance_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (plane_plane_angle_colvar_id)
          CALL plane_plane_angle_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (rotation_colvar_id)
          CALL rotation_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (qparm_colvar_id)
          CALL qparm_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (hydronium_colvar_id)
          CALL hydronium_colvar(colvar,cell,subsys=subsys, error=error)
       CASE(rmsd_colvar_id)
          CALL rmsd_colvar(colvar,subsys=subsys,error=error)
       CASE (reaction_path_colvar_id)
          CALL reaction_path_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (distance_from_path_colvar_id)
          CALL distance_from_path_colvar(colvar,cell, subsys=subsys, error=error)
       CASE (combine_colvar_id)
          CALL combine_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (xyz_diag_colvar_id)
          CALL xyz_diag_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (xyz_outerdiag_colvar_id)
          CALL xyz_outerdiag_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (u_colvar_id)
          CALL u_colvar(colvar,force_env=force_env,subsys=subsys,error=error)
       CASE (Wc_colvar_id)
          CALL Wc_colvar(colvar,cell,subsys=subsys, qs_env=qs_env, error=error)
       CASE (HBP_colvar_id)
          CALL HBP_colvar(colvar,cell,subsys=subsys, qs_env=qs_env, error=error)
       CASE (ring_puckering_colvar_id)
          CALL ring_puckering_colvar(colvar,cell,subsys=subsys, error=error)
       CASE (mindist_colvar_id)
          CALL mindist_colvar(colvar,cell,subsys=subsys , error=error)
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
!> \param cell ...
!> \param particles ...
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
       CASE (plane_plane_angle_colvar_id)
          CALL plane_plane_angle_colvar(colvar,cell, particles=particles, error=error)
       CASE (rotation_colvar_id)
          CALL rotation_colvar(colvar,cell, particles=particles, error=error)
       CASE (qparm_colvar_id)
          CALL qparm_colvar(colvar,cell, particles=particles, error=error)
       CASE (hydronium_colvar_id)
          CALL hydronium_colvar(colvar,cell,particles=particles, error=error)
       CASE(rmsd_colvar_id)
          CALL rmsd_colvar(colvar,particles=particles,error=error)
       CASE (reaction_path_colvar_id)
          CALL reaction_path_colvar(colvar,cell,particles=particles, error=error)
       CASE (distance_from_path_colvar_id)
          CALL distance_from_path_colvar(colvar,cell,particles=particles, error=error)
       CASE (combine_colvar_id)
          CALL combine_colvar(colvar,cell,particles=particles, error=error)
       CASE (xyz_diag_colvar_id)
          CALL xyz_diag_colvar(colvar,cell,particles=particles, error=error)
       CASE (xyz_outerdiag_colvar_id)
          CALL xyz_outerdiag_colvar(colvar,cell,particles=particles, error=error)
       CASE (ring_puckering_colvar_id)
          CALL ring_puckering_colvar(colvar,cell,particles=particles, error=error)
       CASE (mindist_colvar_id)
          CALL mindist_colvar(colvar,cell,particles=particles, error=error)
       CASE (u_colvar_id)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="need force_env!", error=error, error_level=cp_failure_level)
       CASE (Wc_colvar_id)
          CALL Wc_colvar(colvar, cell, particles=particles, error=error)
       CASE (HBP_colvar_id)
          CALL HBP_colvar(colvar, cell, particles=particles, error=error)
       CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
    END IF
  END SUBROUTINE colvar_recursive_eval

! *****************************************************************************
!> \brief Get coordinates of atoms or of geometrical points
!> \param colvar ...
!> \param i ...
!> \param ri ...
!> \param my_particles ...
!> \author Teodoro Laino 03.2007 [created]
! *****************************************************************************
  SUBROUTINE get_coordinates(colvar, i, ri, my_particles)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER, INTENT(IN)                      :: i
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: ri
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    IF (colvar%use_points) THEN
       CALL eval_point_pos(colvar%points(i),my_particles,ri)
    ELSE
       ri(:) = my_particles(i)%r(:)
    END IF

  END SUBROUTINE get_coordinates


! *****************************************************************************
!> \brief Get masses of atoms or of geometrical points
!> \param colvar ...
!> \param i ...
!> \param mi ...
!> \param my_particles ...
!> \author Teodoro Laino 03.2007 [created]
! *****************************************************************************
  SUBROUTINE get_mass(colvar, i, mi, my_particles)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER, INTENT(IN)                      :: i
    REAL(KIND=dp), INTENT(OUT)               :: mi
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    IF (colvar%use_points) THEN
       CALL eval_point_mass(colvar%points(i),my_particles,mi)
    ELSE
       mi = my_particles(i)%atomic_kind%mass
    END IF

  END SUBROUTINE get_mass

! *****************************************************************************
!> \brief Transfer derivatives to ds/dr
!> \param colvar ...
!> \param i ...
!> \param fi ...
!> \author Teodoro Laino 03.2007 [created]
! *****************************************************************************
  SUBROUTINE put_derivative(colvar, i, fi)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER, INTENT(IN)                      :: i
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: fi

    IF (colvar%use_points) THEN
       CALL eval_point_der(colvar%points,i,colvar%dsdr,fi)
    ELSE
       colvar%dsdr(:,i) = colvar%dsdr(:,i) + fi
    END IF

  END SUBROUTINE put_derivative

! *****************************************************************************
!> \brief  evaluates the force due to the position colvar
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
!> \author Teodoro Laino 02.2010 [created]
! *****************************************************************************
  SUBROUTINE xyz_diag_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xyz_diag_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(dp)                                 :: fi(3), r, r0(3), ss(3), &
                                                xi(3), xpi(3)
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    NULLIFY(particles_i)

    CPPrecondition(colvar%type_id==xyz_diag_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    i=colvar%xyz_diag_param%i_atom
    ! Atom coordinates
    CALL get_coordinates(colvar, i, xpi, my_particles)
    ! Use the current coordinates as initial coordinates, if no initialization
    ! was performed yet
    IF (.NOT.colvar%xyz_diag_param%use_absolute_position) THEN
       IF (ALL(colvar%xyz_diag_param%r0 == HUGE(0.0_dp))) THEN
          colvar%xyz_diag_param%r0 = xpi
       END IF
       r0 = colvar%xyz_diag_param%r0
    ELSE
       r0 = 0.0_dp
    ENDIF

    IF(colvar%xyz_diag_param%use_pbc)THEN
       ss=MATMUL(cell%h_inv,xpi-r0)
       ss=ss-NINT(ss)
       xi=MATMUL(cell%hmat,ss)
    ELSE
       xi=xpi-r0
    END IF

    IF (.NOT.colvar%xyz_diag_param%use_absolute_position) THEN
       SELECT CASE (colvar%xyz_diag_param%component)
       CASE(do_clv_x)
          xi(2)=0.0_dp
          xi(3)=0.0_dp
       CASE(do_clv_y)
          xi(1)=0.0_dp
          xi(3)=0.0_dp
       CASE(do_clv_z)
          xi(1)=0.0_dp
          xi(2)=0.0_dp
       CASE(do_clv_xy)
          xi(3)=0.0_dp
       CASE(do_clv_xz)
          xi(2)=0.0_dp
       CASE(do_clv_yz)
          xi(1)=0.0_dp
       CASE DEFAULT
          ! do_clv_xyz
       END SELECT

       r=xi(1)**2+xi(2)**2+xi(3)**2
       fi(:)= 2.0_dp*xi
    ELSE
       SELECT CASE (colvar%xyz_diag_param%component)
       CASE(do_clv_x)
          r=xi(1)
          xi(1)=1.0_dp
          xi(2)=0.0_dp
          xi(3)=0.0_dp
       CASE(do_clv_y)
          r=xi(2)
          xi(1)=0.0_dp
          xi(2)=1.0_dp
          xi(3)=0.0_dp
       CASE(do_clv_z)
          r=xi(3)
          xi(1)=0.0_dp
          xi(2)=0.0_dp
          xi(3)=1.0_dp
       CASE DEFAULT
          !Not implemented for anything which is not a single component.
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       fi(:)= xi
    ENDIF

    colvar%ss=r
    CALL put_derivative(colvar, 1, fi)

  END SUBROUTINE xyz_diag_colvar

! *****************************************************************************
!> \brief  evaluates the force due to the position colvar
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
!> \author Teodoro Laino 02.2010 [created]
! *****************************************************************************
  SUBROUTINE xyz_outerdiag_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xyz_outerdiag_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k, l
    LOGICAL                                  :: failure
    REAL(dp)                                 :: fi(3,2), r, r0(3), ss(3), &
                                                xi(3,2), xpi(3)
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    NULLIFY(particles_i)

    CPPrecondition(colvar%type_id==xyz_outerdiag_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    DO k = 1, 2
       i=colvar%xyz_outerdiag_param%i_atoms(k)
       ! Atom coordinates
       CALL get_coordinates(colvar, i, xpi, my_particles)
       r0 = colvar%xyz_outerdiag_param%r0(:,k)
       IF (ALL(colvar%xyz_outerdiag_param%r0(:,k)==HUGE(0.0_dp))) r0=xpi

       IF(colvar%xyz_outerdiag_param%use_pbc)THEN
          ss=MATMUL(cell%h_inv,xpi-r0)
          ss=ss-NINT(ss)
          xi(:,k)=MATMUL(cell%hmat,ss)
       ELSE
          xi(:,k)=xpi-r0
       END IF

       SELECT CASE (colvar%xyz_outerdiag_param%components(k))
       CASE(do_clv_x)
          xi(2,k)=0.0_dp
          xi(3,k)=0.0_dp
       CASE(do_clv_y)
          xi(1,k)=0.0_dp
          xi(3,k)=0.0_dp
       CASE(do_clv_z)
          xi(1,k)=0.0_dp
          xi(2,k)=0.0_dp
       CASE(do_clv_xy)
          xi(3,k)=0.0_dp
       CASE(do_clv_xz)
          xi(2,k)=0.0_dp
       CASE(do_clv_yz)
          xi(1,k)=0.0_dp
       CASE DEFAULT
          ! do_clv_xyz
       END SELECT
    END DO

    r =0.0_dp
    fi=0.0_dp
    DO i =1,3
       DO l =1,3
          IF (xi(l,1)/=0.0_dp) fi(l,1)=fi(l,1)+xi(i,2)
          r = r + xi(l,1)*xi(i,2)
       END DO
       IF (xi(i,2)/=0.0_dp) fi(i,2)=SUM(xi(:,1))
    END DO

    colvar%ss=r
    CALL put_derivative(colvar, 1, fi(:,1))
    CALL put_derivative(colvar, 2, fi(:,2))

  END SUBROUTINE xyz_outerdiag_colvar

! *****************************************************************************
!> \brief evaluates the force due (and on) the energy as collective variable
!> \param colvar ...
!> \param force_env ...
!> \param subsys ...
!> \param error ...
!> \par History Modified to allow functions of energy in a mixed_env environment
!>              Teodoro Laino [tlaino] - 02.2011
!> \author Sebastiano Caravati
! *****************************************************************************
  SUBROUTINE u_colvar(colvar,force_env,subsys,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'u_colvar', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: coupling_function
    CHARACTER(LEN=default_string_length)     :: def_error, this_error
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: parameters
    INTEGER                                  :: iatom, iforce_eval, &
                                                iparticle, jparticle, natom, &
                                                natom_iforce, nforce_eval, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: glob_natoms, map_index
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dedf, dx, err, fi(3), lerr, &
                                                potential_energy
    REAL(KIND=dp), DIMENSION(:), POINTER     :: values
    TYPE(cp_subsys_p_type), DIMENSION(:), &
      POINTER                                :: subsystems
    TYPE(cp_subsys_type), POINTER            :: subsys_main
    TYPE(mixed_force_type), DIMENSION(:), &
      POINTER                                :: global_forces
    TYPE(particle_list_p_type), &
      DIMENSION(:), POINTER                  :: particles
    TYPE(particle_list_type), POINTER        :: particles_main
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                mapping_section, wrk_section

    failure = .FALSE.
    IF (PRESENT(force_env)) THEN
       NULLIFY(particles_main, subsys_main)
       CALL force_env_get(force_env=force_env,subsys=subsys_main, error=error)
       CALL cp_subsys_get(subsys=subsys_main, particles=particles_main,error=error)
       natom = SIZE(particles_main%els)
       colvar%n_atom_s      = natom
       colvar%u_param%natom = natom
       CALL reallocate(colvar%i_atom,  1,natom)
       CALL reallocate(colvar%dsdr,1,3,1,natom)
       DO iatom=1,natom
          colvar%i_atom(iatom)=iatom
       ENDDO

       IF (.NOT.ASSOCIATED(colvar%u_param%mixed_energy_section)) THEN
          CALL force_env_get(force_env,potential_energy=potential_energy,error=error)
          colvar%ss            = potential_energy

          DO iatom=1,natom
             ! store derivative
             fi(:)=-particles_main%els(iatom)%f
             CALL put_derivative(colvar, iatom, fi)
          ENDDO
       ELSE
          CALL cp_assert(force_env%in_use==use_mixed_force,cp_fatal_level,-300,routineP,&
               'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
               ' A combination of mixed force_eval energies has been requested as '//&
               ' collective variable, but the MIXED env is not in use! Aborting.',&
               error=error,only_ionode=.TRUE.)
          CALL force_env_get(force_env, force_env_section=force_env_section, error=error)
          mapping_section => section_vals_get_subs_vals(force_env_section,"MIXED%MAPPING",error=error)
          NULLIFY(values, parameters, subsystems, particles, global_forces, map_index, glob_natoms)
          nforce_eval = SIZE(force_env%sub_force_env)
          ALLOCATE(glob_natoms(nforce_eval), stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(subsystems(nforce_eval), stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(particles(nforce_eval), stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ! Local Info to sync
          ALLOCATE(global_forces(nforce_eval), stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

          glob_natoms=0
          DO iforce_eval = 1, nforce_eval
             NULLIFY(subsystems(iforce_eval)%subsys, particles(iforce_eval)%list)
             IF (.NOT.ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) CYCLE
             ! Get all available subsys
             CALL force_env_get(force_env=force_env%sub_force_env(iforce_eval)%force_env,&
                                subsys=subsystems(iforce_eval)%subsys,error=error)
             ! Get available particles
             CALL cp_subsys_get(subsys=subsystems(iforce_eval)%subsys,&
                                particles=particles(iforce_eval)%list,error=error)

             ! Get Mapping index array
             natom_iforce = SIZE(particles(iforce_eval)%list%els)

             ! Only the rank 0 process collect info for each computation
             IF ( force_env%sub_force_env(iforce_eval)%force_env%para_env%mepos==&
                  force_env%sub_force_env(iforce_eval)%force_env%para_env%source) THEN
                glob_natoms(iforce_eval) = natom_iforce
             END IF
          END DO

          ! Handling Parallel execution
          CALL mp_sync(force_env%para_env%group)
          CALL mp_sum(glob_natoms, force_env%para_env%group)

          ! Transfer forces
          DO iforce_eval = 1, nforce_eval
             ALLOCATE(global_forces(iforce_eval)%forces(3,glob_natoms(iforce_eval)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             global_forces(iforce_eval)%forces = 0.0_dp
             IF (ASSOCIATED(force_env%sub_force_env(iforce_eval)%force_env)) THEN
                IF ( force_env%sub_force_env(iforce_eval)%force_env%para_env%mepos==&
                     force_env%sub_force_env(iforce_eval)%force_env%para_env%source) THEN
                   ! Forces
                   DO iparticle = 1, glob_natoms(iforce_eval)
                      global_forces(iforce_eval)%forces(:,iparticle) = &
                           particles(iforce_eval)%list%els(iparticle)%f
                   END DO
                END IF
             END IF
             CALL mp_sum(global_forces(iforce_eval)%forces, force_env%para_env%group)
          END DO

          wrk_section => colvar%u_param%mixed_energy_section
          ! Support any number of force_eval sections
          CALL get_generic_info(wrk_section, "ENERGY_FUNCTION", coupling_function, parameters,&
               values, force_env%mixed_env%energies, error=error)
          CALL initf(1)
          CALL parsef(1,TRIM(coupling_function),parameters)
          ! Store the value of the COLVAR
          colvar%ss = evalf(1,values)
          CPPrecondition(EvalErrType<=0,cp_failure_level,routineP,error,failure)

          DO iforce_eval = 1, nforce_eval
             CALL section_vals_val_get(wrk_section,"DX",r_val=dx,error=error)
             CALL section_vals_val_get(wrk_section,"ERROR_LIMIT",r_val=lerr,error=error)
             dedf = evalfd(1,iforce_eval,values,dx,err)
             IF (ABS(err)>lerr) THEN
                WRITE(this_error,"(A,G12.6,A)")"(",err,")"
                WRITE(def_error,"(A,G12.6,A)")"(",lerr,")"
                CALL compress(this_error,.TRUE.)
                CALL compress(def_error,.TRUE.)
                CALL cp_assert(.FALSE.,cp_warning_level,-300,routineP,&
                     'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
                     ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
                     TRIM(def_error)//' .',error=error,only_ionode=.TRUE.)
             END IF
             ! General Mapping of forces...
             ! First: Get Mapping index array
             CALL get_subsys_map_index(mapping_section, glob_natoms(iforce_eval), iforce_eval,&
                                       nforce_eval, map_index, error)

             ! Second: store derivatives
             DO iparticle = 1, glob_natoms(iforce_eval)
                jparticle = map_index(iparticle)
                fi= - dedf * global_forces(iforce_eval)%forces(:,iparticle)
                CALL put_derivative(colvar, jparticle, fi)
             END DO
             ! Deallocate map_index array
             IF (ASSOCIATED(map_index)) THEN
                DEALLOCATE(map_index, stat=stat)
                CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END DO
          CALL finalizef()
          DO iforce_eval = 1, nforce_eval
             DEALLOCATE(global_forces(iforce_eval)%forces,stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          END DO
          DEALLOCATE(glob_natoms, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(values, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(parameters, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(global_forces, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(subsystems, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(particles, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    ELSE
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="need force_env!", error=error, error_level=cp_failure_level)
    ENDIF
  END SUBROUTINE u_colvar

! *****************************************************************************
!> \brief evaluates the force due (and on) the distance from the plane collective variable
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
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
    xpij=ri-rj
    xpkj=rk-rj
    xpl=rl-(ri+rj+rk)/3.0_dp
    IF (colvar%plane_distance_param%use_pbc) THEN
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
    END IF
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
!> \brief evaluates the force due (and on) the angle between two planes.
!>        plane-plane angle collective variable
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
!> \author Teodoro Laino 02.2009 [created]
! *****************************************************************************
  SUBROUTINE plane_plane_angle_colvar(colvar,cell,subsys,particles,error)

    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'plane_plane_angle_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, j1, j2, k1, k2, np
    LOGICAL                                  :: check, failure
    REAL(dp) :: a1, a2, d, dnorm_dxpn(3), dprod12_dxpn(3), dsdxpn(3), &
      dt_dxpn(3), dxpndxi(3,3), dxpndxj(3,3), dxpndxk(3,3), fi(3), fj(3), &
      fk(3), fmod, norm1, norm2, prod_12, ri1(3), ri2(3), rj1(3), rj2(3), &
      rk1(3), rk2(3), ss(3), t, xpij1(3), xpij2(3), xpkj1(3), xpkj2(3), &
      xpn1(3), xpn2(3)
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    NULLIFY(particles_i)

    check = colvar%type_id==plane_plane_angle_colvar_id
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF

    ! Plane 1
    IF (colvar%plane_plane_angle_param%plane1%type_of_def==plane_def_atoms) THEN
       i1=colvar%plane_plane_angle_param%plane1%points(1)
       j1=colvar%plane_plane_angle_param%plane1%points(2)
       k1=colvar%plane_plane_angle_param%plane1%points(3)

       ! Get coordinates of atoms or points
       CALL get_coordinates(colvar, i1, ri1, my_particles)
       CALL get_coordinates(colvar, j1, rj1, my_particles)
       CALL get_coordinates(colvar, k1, rk1, my_particles)

       ! xpij
       ss=MATMUL(cell%h_inv,ri1-rj1)
       ss=ss-NINT(ss)
       xpij1=MATMUL(cell%hmat,ss)

       ! xpkj
       ss=MATMUL(cell%h_inv,rk1-rj1)
       ss=ss-NINT(ss)
       xpkj1=MATMUL(cell%hmat,ss)

       ! xpn
       xpn1(1)  = xpij1(2)*xpkj1(3)-xpij1(3)*xpkj1(2)
       xpn1(2)  = xpij1(3)*xpkj1(1)-xpij1(1)*xpkj1(3)
       xpn1(3)  = xpij1(1)*xpkj1(2)-xpij1(2)*xpkj1(1)
    ELSE
       xpn1     = colvar%plane_plane_angle_param%plane1%normal_vec
    END IF
    a1    = DOT_PRODUCT(xpn1,xpn1)
    norm1 = SQRT(a1)
    CPPrecondition(norm1/=0.0_dp,cp_failure_level,routineP,error,failure)

    ! Plane 2
    IF (colvar%plane_plane_angle_param%plane2%type_of_def==plane_def_atoms) THEN
       i2=colvar%plane_plane_angle_param%plane2%points(1)
       j2=colvar%plane_plane_angle_param%plane2%points(2)
       k2=colvar%plane_plane_angle_param%plane2%points(3)

       ! Get coordinates of atoms or points
       CALL get_coordinates(colvar, i2, ri2, my_particles)
       CALL get_coordinates(colvar, j2, rj2, my_particles)
       CALL get_coordinates(colvar, k2, rk2, my_particles)

       ! xpij
       ss=MATMUL(cell%h_inv,ri2-rj2)
       ss=ss-NINT(ss)
       xpij2=MATMUL(cell%hmat,ss)

       ! xpkj
       ss=MATMUL(cell%h_inv,rk2-rj2)
       ss=ss-NINT(ss)
       xpkj2=MATMUL(cell%hmat,ss)

       ! xpn
       xpn2(1)  = xpij2(2)*xpkj2(3)-xpij2(3)*xpkj2(2)
       xpn2(2)  = xpij2(3)*xpkj2(1)-xpij2(1)*xpkj2(3)
       xpn2(3)  = xpij2(1)*xpkj2(2)-xpij2(2)*xpkj2(1)
    ELSE
       xpn2     = colvar%plane_plane_angle_param%plane2%normal_vec
    END IF
    a2    = DOT_PRODUCT(xpn2,xpn2)
    norm2 = SQRT(a2)
    CPPrecondition(norm2/=0.0_dp,cp_failure_level,routineP,error,failure)

    ! The value of the angle is defined only between 0 and Pi
    prod_12  = DOT_PRODUCT(xpn1,xpn2)

    d        = norm1*norm2
    t        = prod_12/d
    t        = MIN(1.0_dp,ABS(t))*SIGN(1.0_dp,t)
    colvar%ss= ACOS(t)

    IF ((ABS(colvar%ss).LT.tolerance_acos).OR.(ABS(colvar%ss-pi).LT.tolerance_acos)) THEN
       fmod =   0.0_dp
    ELSE
       fmod = - 1.0_dp / SIN(colvar%ss)
    ENDIF
    ! Compute derivatives
    np = 0
    ! Plane 1
    IF (colvar%plane_plane_angle_param%plane1%type_of_def==plane_def_atoms) THEN
       dprod12_dxpn = xpn2
       dnorm_dxpn   = 1.0_dp/norm1*xpn1
       dt_dxpn      = (dprod12_dxpn * d - prod_12 * dnorm_dxpn * norm2)/d**2

       dsdxpn(1) = fmod * dt_dxpn(1)
       dsdxpn(2) = fmod * dt_dxpn(2)
       dsdxpn(3) = fmod * dt_dxpn(3)
       !
       dxpndxi(1,1)=  0.0_dp
       dxpndxi(1,2)=  1.0_dp * xpkj1(3)
       dxpndxi(1,3)= -1.0_dp * xpkj1(2)
       dxpndxi(2,1)= -1.0_dp * xpkj1(3)
       dxpndxi(2,2)=  0.0_dp
       dxpndxi(2,3)=  1.0_dp * xpkj1(1)
       dxpndxi(3,1)=  1.0_dp * xpkj1(2)
       dxpndxi(3,2)= -1.0_dp * xpkj1(1)
       dxpndxi(3,3)=  0.0_dp
       !
       dxpndxj(1,1)=  0.0_dp
       dxpndxj(1,2)= -1.0_dp * xpkj1(3) + xpij1(3)
       dxpndxj(1,3)= -1.0_dp * xpij1(2) + xpkj1(2)
       dxpndxj(2,1)= -1.0_dp * xpij1(3) + xpkj1(3)
       dxpndxj(2,2)=  0.0_dp
       dxpndxj(2,3)= -1.0_dp * xpkj1(1) + xpij1(1)
       dxpndxj(3,1)= -1.0_dp * xpkj1(2) + xpij1(2)
       dxpndxj(3,2)= -1.0_dp * xpij1(1) + xpkj1(1)
       dxpndxj(3,3)=  0.0_dp
       !
       dxpndxk(1,1)=  0.0_dp
       dxpndxk(1,2)= -1.0_dp * xpij1(3)
       dxpndxk(1,3)=  1.0_dp * xpij1(2)
       dxpndxk(2,1)=  1.0_dp * xpij1(3)
       dxpndxk(2,2)=  0.0_dp
       dxpndxk(2,3)= -1.0_dp * xpij1(1)
       dxpndxk(3,1)= -1.0_dp * xpij1(2)
       dxpndxk(3,2)=  1.0_dp * xpij1(1)
       dxpndxk(3,3)=  0.0_dp
       !
       fi=MATMUL(dsdxpn,dxpndxi)
       fj=MATMUL(dsdxpn,dxpndxj)
       fk=MATMUL(dsdxpn,dxpndxk)

       ! Transfer derivatives on atoms
       CALL put_derivative(colvar, np+1, fi)
       CALL put_derivative(colvar, np+2, fj)
       CALL put_derivative(colvar, np+3, fk)
       np = 3
    END IF

    ! Plane 2
    IF (colvar%plane_plane_angle_param%plane2%type_of_def==plane_def_atoms) THEN
       dprod12_dxpn = xpn1
       dnorm_dxpn   = 1.0_dp/norm2*xpn2
       dt_dxpn      = (dprod12_dxpn * d - prod_12 * dnorm_dxpn * norm1)/d**2

       dsdxpn(1) = fmod * dt_dxpn(1)
       dsdxpn(2) = fmod * dt_dxpn(2)
       dsdxpn(3) = fmod * dt_dxpn(3)
       !
       dxpndxi(1,1)=  0.0_dp
       dxpndxi(1,2)=  1.0_dp * xpkj1(3)
       dxpndxi(1,3)= -1.0_dp * xpkj1(2)
       dxpndxi(2,1)= -1.0_dp * xpkj1(3)
       dxpndxi(2,2)=  0.0_dp
       dxpndxi(2,3)=  1.0_dp * xpkj1(1)
       dxpndxi(3,1)=  1.0_dp * xpkj1(2)
       dxpndxi(3,2)= -1.0_dp * xpkj1(1)
       dxpndxi(3,3)=  0.0_dp
       !
       dxpndxj(1,1)=  0.0_dp
       dxpndxj(1,2)= -1.0_dp * xpkj1(3) + xpij1(3)
       dxpndxj(1,3)= -1.0_dp * xpij1(2) + xpkj1(2)
       dxpndxj(2,1)= -1.0_dp * xpij1(3) + xpkj1(3)
       dxpndxj(2,2)=  0.0_dp
       dxpndxj(2,3)= -1.0_dp * xpkj1(1) + xpij1(1)
       dxpndxj(3,1)= -1.0_dp * xpkj1(2) + xpij1(2)
       dxpndxj(3,2)= -1.0_dp * xpij1(1) + xpkj1(1)
       dxpndxj(3,3)=  0.0_dp
       !
       dxpndxk(1,1)=  0.0_dp
       dxpndxk(1,2)= -1.0_dp * xpij1(3)
       dxpndxk(1,3)=  1.0_dp * xpij1(2)
       dxpndxk(2,1)=  1.0_dp * xpij1(3)
       dxpndxk(2,2)=  0.0_dp
       dxpndxk(2,3)= -1.0_dp * xpij1(1)
       dxpndxk(3,1)= -1.0_dp * xpij1(2)
       dxpndxk(3,2)=  1.0_dp * xpij1(1)
       dxpndxk(3,3)=  0.0_dp
       !
       fi=MATMUL(dsdxpn,dxpndxi)
       fj=MATMUL(dsdxpn,dxpndxj)
       fk=MATMUL(dsdxpn,dxpndxk)

       ! Transfer derivatives on atoms
       CALL put_derivative(colvar, np+1, fi)
       CALL put_derivative(colvar, np+2, fj)
       CALL put_derivative(colvar, np+3, fk)
    END IF

  END SUBROUTINE plane_plane_angle_colvar

! *****************************************************************************
!> \brief Evaluates the value of the rotation angle between two bonds
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
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
    IF ((ABS(colvar%ss).LT.tolerance_acos).OR.(ABS(colvar%ss-pi).LT.tolerance_acos)) THEN
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
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
!> \author Teodoro Laino 02.2006 [created]
!> \note modified Florian Schiffmann 08.2008
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
    IF(colvar%dfunct_param%use_pbc)THEN
       ss=MATMUL(cell%h_inv,xpi-xpj)
       ss=ss-NINT(ss)
       xij=MATMUL(cell%hmat,ss)
    ELSE
       xij=xpi-xpj
    END IF
    r12=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
    ! Second bond
    k=colvar%dfunct_param%i_at_dfunct(3)
    l=colvar%dfunct_param%i_at_dfunct(4)
    CALL get_coordinates(colvar, k, xpk, my_particles)
    CALL get_coordinates(colvar, l, xpl, my_particles)
    IF(colvar%dfunct_param%use_pbc)THEN
       ss=MATMUL(cell%h_inv,xpk-xpl)
       ss=ss-NINT(ss)
       xkl=MATMUL(cell%hmat,ss)
    ELSE
       xkl=xpk-xpl
    END IF
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
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
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
    IF ((ABS(colvar%ss).LT.tolerance_acos).OR.(ABS(colvar%ss-pi).LT.tolerance_acos)) THEN
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
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
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
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param no_riemann_sheet_op ...
!> \param error ...
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
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
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
             IF ((l+m) > maxfac) THEN
                CALL stop_program(routineN,moduleN,__LINE__,&
                                  "(l-m) > maxfac")
             END IF
             lmm = fac(l-m)
             lpm = fac(l+m)
             bond  = 1.0_dp/(1.0_dp + EXP(alpha*(r_tmp - rcut)))
             IF (bond > 1.0_dp) THEN
                CALL stop_program(routineN,moduleN,__LINE__,&
                                  "bond > 1.0_dp")
             END IF
             pre_fac = bond*SQRT(((2*l+1)*lmm)/(4*pi*lpm))
             IF (plm < 0.0_dp) THEN
                dylm = -pre_fac*dplm
             ELSE
                dylm = pre_fac*dplm
             ENDIF
             ylm_calc = pre_fac*plm*(CMPLX(COS(m*fi),SIN(m*fi),KIND=dp))
             ylm = ylm + ylm_calc
             ylm_calc2 = pre_fac/bond*plm*(CMPLX(COS(m*fi),SIN(m*fi),KIND=dp))
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
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "NBOND = 0.0")
          END IF
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
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
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
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
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

    INTEGER                                  :: i, ii, j, jj, k, kk, &
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
          j=colvar%coord_param%i_at_to(jj)
          CALL get_coordinates(colvar, j, xpj, my_particles)
          ! define coordination of atom A with itself to be 0. also fixes rij==0 for the force calculation
          IF (i.EQ.j) CYCLE
          ss=MATMUL(cell%h_inv,xpi(:)-xpj(:))
          ss=ss-NINT(ss)
          xij=MATMUL(cell%hmat,ss)
          rij=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
          IF(rij < 1.0e-8_dp)CYCLE
          rdist_ij = rij/r_0_a
          IF(ABS(1.0_dp-rdist_ij) > EPSILON(0.0_dp)*1.0E+4_dp) THEN
             num_ij=(1.0_dp-rdist_ij**p_a)
             invden_ij=1.0_dp/(1.0_dp-rdist_ij**q_a)
             func_ij=num_ij*invden_ij
             IF (rij < 1.0E-8_dp) THEN
                ! provide the correct limit of the derivative
                dfunc_ij = 0.0_dp
             ELSE
                dfunc_ij= (- p_a*rdist_ij **(p_a-1) *invden_ij  &
                     + num_ij*(invden_ij)**2 * q_a *rdist_ij **(q_a-1))/(rij*r_0_a)
             END IF
          ELSE
             ! Provide the correct limit for function value and derivative
             func_ij  = REAL(p_a,KIND=dp)/REAL(q_a,KIND=dp)
             dfunc_ij = REAL(p_a,KIND=dp)*REAL((-q_a+p_a),KIND=dp)/(REAL(2*q_a,KIND=dp)*r_0_a)
          END IF
          IF(n_atoms_to_b /=0) THEN
             func_k = 0.0_dp
             DO kk = 1,n_atoms_to_b
                k=colvar%coord_param%i_at_to_b(kk)
                IF (k.EQ.j) CYCLE
                CALL get_coordinates(colvar, k, xpk, my_particles)
                ss=MATMUL(cell%h_inv,xpj(:)-xpk(:))
                ss=ss-NINT(ss)
                xjk=MATMUL(cell%hmat,ss)
                rjk=SQRT(xjk(1)**2+xjk(2)**2+xjk(3)**2)
                IF(rjk < 1.0e-8_dp)CYCLE
                rdist_jk = rjk/r_0_b
                IF(ABS(1.0_dp-rdist_jk) > EPSILON(0.0_dp)*1.0E+4_dp) THEN
                   num_jk=(1.0_dp-rdist_jk**p_b)
                   invden_jk=1.0_dp/(1.0_dp-rdist_jk**q_b)
                   func_jk=num_jk*invden_jk
                   IF (rjk < 1.0E-8_dp) THEN
                      ! provide the correct limit of the derivative
                      dfunc_jk = 0.0_dp
                   ELSE
                      dfunc_jk= (- p_b*rdist_jk **(p_b-1) *invden_jk  &
                           + num_jk*(invden_jk)**2 * q_b *rdist_jk **(q_b-1))/(rjk*r_0_b)
                   END IF
                ELSE
                   ! Provide the correct limit for function value and derivative
                   func_jk  = REAL(p_b,KIND=dp)/REAL(q_b,KIND=dp)
                   dfunc_jk = REAL(p_b,KIND=dp)*REAL((-q_b+p_b),KIND=dp)/(REAL(2*q_b,KIND=dp)*r_0_b)
                ENDIF
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
!> \brief ...
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE mindist_colvar(colvar,cell, subsys, particles, error)

    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mindist_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, j, jj, n_coord_from, &
                                                n_coord_to, n_dist_from, p, &
                                                q, stat
    LOGICAL                                  :: failure
    REAL(dp) :: den_n, den_Q, fscalar, ftemp_i(3), inv_den_n, inv_den_Q, &
      lambda, num_n, num_Q, Qfunc, r12, r_cut, rfact, rij(3), rpi(3), rpj(3)
    REAL(dp), DIMENSION(:), POINTER          :: dqfunc_dnL, expnL, nLcoord, &
                                                sum_rij
    REAL(dp), DIMENSION(:, :, :), POINTER    :: dnLcoord, dqfunc_dr
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    ! If we defined the coordination number with KINDS then we have still
    ! to fill few missing informations...
    NULLIFY(particles_i)
    CPPrecondition(colvar%type_id==mindist_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF

    n_dist_from=colvar%mindist_param%n_dist_from
    n_coord_from=colvar%mindist_param%n_coord_from
    n_coord_to=colvar%mindist_param%n_coord_to
    p=colvar%mindist_param%p_exp
    q=colvar%mindist_param%q_exp
    r_cut=colvar%mindist_param%r_cut
    lambda=colvar%mindist_param%lambda


    NULLIFY(nLcoord,dnLcoord,dqfunc_dr,dqfunc_dnL,expnL,sum_rij)
    ALLOCATE (nLcoord(n_coord_from),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE (dnLcoord(3,n_coord_from,n_coord_to),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE (expnL(n_coord_from),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE (sum_rij(n_coord_from),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE (dqfunc_dr(3,n_dist_from,n_coord_from),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE (dqfunc_dnL(n_coord_from),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    ! coordination numbers
    nLcoord = 0.0_dp
    dnLcoord = 0.0_dp 
    expnL = 0.0_dp
    den_Q = 0.0_dp
    DO i = 1,n_coord_from
      ii = colvar%mindist_param%i_coord_from(i)
      rpi=my_particles(ii)%r(1:3)
      DO j = 1,n_coord_to
         jj = colvar%mindist_param%i_coord_to(j)
         rpj=my_particles(jj)%r(1:3)
         rij = pbc(rpj,rpi,cell)
         r12 = SQRT(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
         rfact=r12/r_cut
         num_n = 1.0_dp-rfact**p
         den_n = 1.0_dp-rfact**q
         inv_den_n = 1.0_dp/den_n
         IF(ABS(inv_den_n)<1.e-10_dp ) THEN
             inv_den_n=1.e-10_dp
             num_n = ABS(num_n)
         END IF

         fscalar = (-p*rfact**(p-1)+num_n*q*rfact**(q-1)*inv_den_n)*inv_den_n/(r_cut*r12)

         dnLcoord(1,i,j) = rij(1)*fscalar
         dnLcoord(2,i,j) = rij(2)*fscalar
         dnLcoord(3,i,j) = rij(3)*fscalar

         nLcoord(i) =  nLcoord(i)  + num_n*inv_den_n
      END DO
      expnL(i) =  EXP(lambda*nLcoord(i)) 
!dbg
! write(*,*) ii,nLcoord(i),expnL(i) 
!dbg
      den_Q = den_Q +  expnL(i)
    END DO
    inv_den_Q = 1.0_dp/den_Q

    qfunc = 0.0_dp
    dqfunc_dr = 0.0_dp
    dqfunc_dnL = 0.0_dp
    num_Q = 0.0_dp
    sum_rij = 0.0_dp
    DO i = 1,n_dist_from
      ii =  colvar%mindist_param%i_dist_from(i)
      rpi=my_particles(ii)%r(1:3)
      DO j = 1,n_coord_from
         jj =  colvar%mindist_param%i_coord_from(j)
         rpj=my_particles(jj)%r(1:3)
         rij = pbc(rpj,rpi,cell)
         r12 = SQRT(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))

!dbg
! write(*,*) ii,jj,rpi(1:3),rpj(1:3),rij(1:3),r12
!dbg
         num_Q = num_Q +  r12 *  expnL(j)

         sum_rij(j) = sum_rij(j) + r12
         dqfunc_dr(1,i,j) = expnL(j)* rij(1)/r12
         dqfunc_dr(2,i,j) = expnL(j)* rij(2)/r12
         dqfunc_dr(3,i,j) = expnL(j)* rij(3)/r12

      END DO

    END DO

    ! Function and derivatives
    qfunc = num_Q * inv_den_Q
    dqfunc_dr = dqfunc_dr  * inv_den_Q
    colvar%ss=qfunc
!dbg
!      write(*,*) ' ss ',  colvar%ss
! stop
!dbg

    DO i = 1,n_coord_from
      dqfunc_dnL(i) = lambda * expnL(i) * inv_den_Q*( sum_rij(i) - num_Q * inv_den_Q)
    END DO 

    !Compute Forces
    DO i = 1,n_dist_from
      DO j = 1,n_coord_from
        ftemp_i(1) =  dqfunc_dr(1,i,j)
        ftemp_i(2) =  dqfunc_dr(2,i,j)
        ftemp_i(3) =  dqfunc_dr(3,i,j)

        CALL put_derivative(colvar,i,ftemp_i)
        CALL put_derivative(colvar,j+n_dist_from,-ftemp_i)

      END DO
    END DO
    DO i = 1,n_coord_from
      DO j = 1,n_coord_to
        ftemp_i(1) =  dqfunc_dnL(i)*dnLcoord(1,i,j)
        ftemp_i(2) =  dqfunc_dnL(i)*dnLcoord(2,i,j)
        ftemp_i(3) =  dqfunc_dnL(i)*dnLcoord(3,i,j)

        CALL put_derivative(colvar,i+n_dist_from,ftemp_i)
        CALL put_derivative(colvar,j+n_dist_from+n_coord_from,-ftemp_i)

      END DO
    END DO




    DEALLOCATE(nLcoord,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(dnLcoord,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE (expnL,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE (dqfunc_dr,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE (sum_rij,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE (dqfunc_dnL,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

  END SUBROUTINE mindist_colvar

! *****************************************************************************
!> \brief  evaluates function and forces due to a combination of COLVARs
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
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
!>             ss(R) = [\sum_i i*dt exp{-\lambda \sum_a(S_a(R)-f_a(i))^2}]/
!>                     [\sum_i exp{-\lambda \sum_a(S_a(R)-f_a(i))^2}]
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
!> \par History
!>      extended MI 01.2010
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

    LOGICAL                                  :: failure
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

    IF(colvar%reaction_path_param%dist_rmsd) THEN
      CALL rpath_dist_rmsd(colvar,my_particles,error=error)
    ELSEIF(colvar%reaction_path_param%rmsd) THEN
      CALL rpath_rmsd(colvar,my_particles,error=error)
    ELSE
       CALL rpath_colvar(colvar,cell,my_particles,error)
    END IF

  END SUBROUTINE reaction_path_colvar

! *****************************************************************************
!> \brief  position along the path calculated using selected colvars
!>         as compared to functions describing the variation of these same colvars
!>         along the path given as reference
!> \param colvar ...
!> \param cell ...
!> \param particles ...
!> \param error ...
!> \author fschiff
! *****************************************************************************
  SUBROUTINE rpath_colvar(colvar,cell,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rpath_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iend, ii, istart, j, k, &
                                                ncolv, nconf, stat
    LOGICAL                                  :: failure
    REAL(dp)                                 :: lambda, step_size
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: s1, ss_vals
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: ds1, f_vals, fi, s1v
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ds1v

    istart=colvar%reaction_path_param%function_bounds(1)
    iend=colvar%reaction_path_param%function_bounds(2)

    nconf=colvar%reaction_path_param%nr_frames
    step_size=colvar%reaction_path_param%step_size
    ncolv=colvar%reaction_path_param%n_components
    lambda=colvar%reaction_path_param%lambda
    ALLOCATE(f_vals(ncolv,istart:iend),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    f_vals(:,:)=colvar%reaction_path_param%f_vals
    ALLOCATE(ss_vals(ncolv),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    DO i=1,ncolv
       CALL colvar_recursive_eval(colvar%reaction_path_param%colvar_p(i)%colvar,cell,particles,error)
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
           ds1v(j,1,k) = f_vals(j,k)* s1v(1,k)
           ds1v(j,2,k) = f_vals(j,k) * s1v(2,k)
       END DO
    END DO
    DO i=1,2
       s1(i)=accurate_sum(s1v(i,:))
       DO j=1,ncolv
          ds1(j,i)=accurate_sum(ds1v(j,i,:))
       END DO
    END DO

    colvar%ss=s1(1)/s1(2)/REAL(nconf-1,dp)

    ALLOCATE(fi(3,colvar%n_atom_s),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    ii=0
    DO i=1,ncolv
       DO j=1,colvar%reaction_path_param%colvar_p(i)%colvar%n_atom_s
          ii=ii+1
          fi(:,ii)=colvar%reaction_path_param%colvar_p(i)%colvar%dsdr(:,j)*lambda*(ds1(i,1) &
               /s1(2)/REAL(nconf-1,dp)-colvar%ss*ds1(i,2)/s1(2))*2.0_dp
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

  END SUBROUTINE rpath_colvar

! *****************************************************************************
!> \brief  position along the path calculated from the positions of a selected list of
!>         atoms as compared to the same positions in reference
!>         configurations belonging to the given path.
!> \param colvar ...
!> \param particles ...
!> \param error ...
!> \date  01.2010
!> \author MI
! *****************************************************************************
  SUBROUTINE rpath_dist_rmsd(colvar,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rpath_dist_rmsd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iat, ii, ik, natom, nconf, &
                                                rmsd_atom, stat
    INTEGER, DIMENSION(:), POINTER           :: iatom
    LOGICAL                                  :: failure
    REAL(dp)                                 :: lambda, my_rmsd, s1(2), &
                                                sum_exp
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: r, r0, vec_dif
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: dvec_dif, fi, riat, s1v
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ds1
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: ds1v
    REAL(dp), DIMENSION(:, :), POINTER       :: path_conf

    failure=.FALSE.

    nconf=colvar%reaction_path_param%nr_frames
    rmsd_atom=colvar%reaction_path_param%n_components
    lambda=colvar%reaction_path_param%lambda
    path_conf => colvar%reaction_path_param%r_ref
    iatom => colvar%reaction_path_param%i_rmsd

    natom = SIZE(particles)

    ALLOCATE(r0(3*natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(r(3*natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(riat(3,rmsd_atom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(vec_dif(rmsd_atom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(dvec_dif(3,rmsd_atom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(s1v(2,nconf),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1v(3,rmsd_atom,2,nconf),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1(3,rmsd_atom,2),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DO i = 1,natom
      ii = (i-1)*3
      r0(ii+1) = particles(i)%r(1)
      r0(ii+2) = particles(i)%r(2)
      r0(ii+3) = particles(i)%r(3)
    END DO

    DO iat = 1,rmsd_atom
      ii = iatom(iat)
      riat(:,iat) = particles(ii)%r
    END DO

    DO ik = 1, nconf
       DO i = 1,natom
        ii = (i-1)*3
        r(ii+1) = path_conf(ii+1,ik)
        r(ii+2) = path_conf(ii+2,ik)
        r(ii+3) = path_conf(ii+3,ik)
       END DO

       CALL rmsd3(particles,r,r0,output_unit=-1,my_val=my_rmsd,rotate=.TRUE.,error=error)

       sum_exp = 0.0_dp
       DO iat  = 1,rmsd_atom
         i = iatom(iat)
         ii = (i-1)*3
         vec_dif(iat) = (riat(1,iat)-r(ii+1))**2+(riat(2,iat)-r(ii+2))**2&
                       +(riat(3,iat)-r(ii+3))**2
         sum_exp =  sum_exp + vec_dif(iat)
       END DO

       s1v(1,ik) = REAL(ik-1,dp)*EXP(-lambda*sum_exp)
       s1v(2,ik) = EXP(-lambda*sum_exp)
       DO  iat  = 1,rmsd_atom
         i = iatom(iat)
         ii = (i-1)*3
         ds1v(1,iat,1,ik) = r(ii+1)*s1v(1,ik)
         ds1v(1,iat,2,ik) = r(ii+1)*s1v(2,ik)
         ds1v(2,iat,1,ik) = r(ii+2)*s1v(1,ik)
         ds1v(2,iat,2,ik) = r(ii+2)*s1v(2,ik)
         ds1v(3,iat,1,ik) = r(ii+3)*s1v(1,ik)
         ds1v(3,iat,2,ik) = r(ii+3)*s1v(2,ik)
       END DO

    END DO
    s1(1) = accurate_sum(s1v(1,:))
    s1(2) = accurate_sum(s1v(2,:))
    DO i = 1,2
      DO iat =  1,rmsd_atom
        ds1(1,iat,i) = accurate_sum(ds1v(1,iat,i,:))
        ds1(2,iat,i) = accurate_sum(ds1v(2,iat,i,:))
        ds1(3,iat,i) = accurate_sum(ds1v(3,iat,i,:))
      END DO
    END DO

    colvar%ss = s1(1)/s1(2)/REAL(nconf-1,dp)

    ALLOCATE(fi(3,rmsd_atom),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    DO iat = 1,rmsd_atom
      fi(1,iat) = 2.0_dp*lambda/s1(2)/REAL(nconf-1,dp)*(ds1(1,iat,1) - ds1(1,iat,2)*s1(1)/s1(2))
      fi(2,iat) = 2.0_dp*lambda/s1(2)/REAL(nconf-1,dp)*(ds1(2,iat,1) - ds1(2,iat,2)*s1(1)/s1(2))
      fi(3,iat) = 2.0_dp*lambda/s1(2)/REAL(nconf-1,dp)*(ds1(3,iat,1) - ds1(3,iat,2)*s1(1)/s1(2))
      CALL put_derivative(colvar,iat,fi(:,iat))
    END DO

    DEALLOCATE(fi,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(r0,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(r,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(riat,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(vec_dif,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(dvec_dif,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(s1v,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ds1v,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ds1,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

  END SUBROUTINE rpath_dist_rmsd


! *****************************************************************************
!> \brief ...
!> \param colvar ...
!> \param particles ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE  rpath_rmsd(colvar,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rpath_rmsd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iat, ii, ik, natom, nconf, &
                                                rmsd_atom, stat
    INTEGER, DIMENSION(:), POINTER           :: iatom
    LOGICAL                                  :: failure
    REAL(dp)                                 :: lambda, my_rmsd, s1(2)
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: r, r0
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: fi, riat, s1v
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ds1
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: ds1v
    REAL(dp), DIMENSION(:, :), POINTER       :: path_conf
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: weight
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: drmsd

    failure=.FALSE.

    nconf=colvar%reaction_path_param%nr_frames
    rmsd_atom=colvar%reaction_path_param%n_components
    lambda=colvar%reaction_path_param%lambda
    path_conf => colvar%reaction_path_param%r_ref
    iatom => colvar%reaction_path_param%i_rmsd

    natom = SIZE(particles)

    ALLOCATE(r0(3*natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(r(3*natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(riat(3,rmsd_atom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(s1v(2,nconf),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1v(3,rmsd_atom,2,nconf),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1(3,rmsd_atom,2),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(drmsd(3,natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    drmsd = 0.0_dp
    ALLOCATE (weight(natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    DO i = 1,natom
      ii = (i-1)*3
      r0(ii+1) = particles(i)%r(1)
      r0(ii+2) = particles(i)%r(2)
      r0(ii+3) = particles(i)%r(3)
    END DO

    DO iat = 1,rmsd_atom
      ii = iatom(iat)
      riat(:,iat) = particles(ii)%r
    END DO

! set weights of atoms in the rmsd list
    weight = 0.0_dp
    DO iat  = 1,rmsd_atom
       i = iatom(iat)
       weight(i) = 1.0_dp
    END DO


    DO ik = 1, nconf
       DO i = 1,natom
        ii = (i-1)*3
        r(ii+1) = path_conf(ii+1,ik)
        r(ii+2) = path_conf(ii+2,ik)
        r(ii+3) = path_conf(ii+3,ik)
       END DO

       CALL rmsd3(particles,r0,r,output_unit=-1,weights=weight,my_val=my_rmsd,&
           rotate=.FALSE.,drmsd3=drmsd,error=error)

       s1v(1,ik) = REAL(ik-1,dp)*EXP(-lambda*my_rmsd)
       s1v(2,ik) = EXP(-lambda*my_rmsd)
       DO  iat  = 1,rmsd_atom
         i = iatom(iat)
         ds1v(1,iat,1,ik) = drmsd(1,i)*s1v(1,ik)
         ds1v(1,iat,2,ik) = drmsd(1,i)*s1v(2,ik)
         ds1v(2,iat,1,ik) = drmsd(2,i)*s1v(1,ik)
         ds1v(2,iat,2,ik) = drmsd(2,i)*s1v(2,ik)
         ds1v(3,iat,1,ik) = drmsd(3,i)*s1v(1,ik)
         ds1v(3,iat,2,ik) = drmsd(3,i)*s1v(2,ik)
       END DO
    END DO  ! ik

    s1(1) = accurate_sum(s1v(1,:))
    s1(2) = accurate_sum(s1v(2,:))
    DO i = 1,2
      DO iat =  1,rmsd_atom
        ds1(1,iat,i) = accurate_sum(ds1v(1,iat,i,:))
        ds1(2,iat,i) = accurate_sum(ds1v(2,iat,i,:))
        ds1(3,iat,i) = accurate_sum(ds1v(3,iat,i,:))
      END DO
    END DO

   colvar%ss = s1(1)/s1(2)/REAL(nconf-1,dp)

    ALLOCATE(fi(3,rmsd_atom),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    DO iat = 1,rmsd_atom
      fi(1,iat) = -lambda/s1(2)/REAL(nconf-1,dp)*(ds1(1,iat,1) - ds1(1,iat,2)*s1(1)/s1(2))
      fi(2,iat) = -lambda/s1(2)/REAL(nconf-1,dp)*(ds1(2,iat,1) - ds1(2,iat,2)*s1(1)/s1(2))
      fi(3,iat) = -lambda/s1(2)/REAL(nconf-1,dp)*(ds1(3,iat,1) - ds1(3,iat,2)*s1(1)/s1(2))
      CALL put_derivative(colvar,iat,fi(:,iat))
    END DO

    DEALLOCATE(fi,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(r0,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(r,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(riat,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(s1v,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ds1v,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ds1,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(drmsd,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(weight,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

  END SUBROUTINE  rpath_rmsd

! *****************************************************************************
!> \brief evaluates the force due (and on) distance from reaction path collective variable
!>             ss(R) = -1/\lambda \log[\sum_i exp{-\lambda \sum_a(S_a(R)-f_a(i))^2}]
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
!> \date 01.2010
!> \author MI
! *****************************************************************************
  SUBROUTINE  distance_from_path_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distance_from_path_colvar', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.

    CPPrecondition(colvar%type_id==distance_from_path_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF

    IF(colvar%reaction_path_param%dist_rmsd) THEN
      CALL dpath_dist_rmsd(colvar,my_particles,error=error)
    ELSEIF(colvar%reaction_path_param%rmsd) THEN
      CALL dpath_rmsd(colvar,my_particles,error=error)
    ELSE
      CALL dpath_colvar(colvar,cell,my_particles,error=error)
    END IF

  END SUBROUTINE  distance_from_path_colvar

! *****************************************************************************
!> \brief  distance from path calculated using selected colvars
!>         as compared to functions describing the variation of these same colvars
!>         along the path given as reference
!> \param colvar ...
!> \param cell ...
!> \param particles ...
!> \param error ...
!> \date  01.2010
!> \author MI
! *****************************************************************************
  SUBROUTINE  dpath_colvar(colvar,cell,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dpath_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iend, ii, istart, j, k, &
                                                ncolv, stat
    LOGICAL                                  :: failure
    REAL(dp)                                 :: lambda, s1
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: ds1, s1v, ss_vals
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: ds1v, f_vals, fi

    istart=colvar%reaction_path_param%function_bounds(1)
    iend=colvar%reaction_path_param%function_bounds(2)

    ncolv=colvar%reaction_path_param%n_components
    lambda=colvar%reaction_path_param%lambda
    ALLOCATE(f_vals(ncolv,istart:iend),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    f_vals(:,:)=colvar%reaction_path_param%f_vals
    ALLOCATE(ss_vals(ncolv),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    DO i=1,ncolv
       CALL colvar_recursive_eval(colvar%reaction_path_param%colvar_p(i)%colvar,cell,particles,error)
       ss_vals(i)=colvar%reaction_path_param%colvar_p(i)%colvar%ss
    ENDDO

    ALLOCATE(s1v(istart:iend),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1v(ncolv,istart:iend),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1(ncolv),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    DO k=istart,iend
       s1v(k)=EXP(-lambda*DOT_PRODUCT(ss_vals(:)-f_vals(:,k),ss_vals(:)-f_vals(:,k)))
       DO j=1,ncolv
         ds1v(j,k)=f_vals(j,k)*s1v(k)
       END DO
    END DO

    s1 = accurate_sum(s1v(:))
    DO j=1,ncolv
      ds1(j) = accurate_sum(ds1v(j,:))
    END DO
    colvar%ss=-1.0_dp/lambda*LOG(s1)

    ALLOCATE(fi(3,colvar%n_atom_s),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    ii=0
    DO i=1,ncolv
       DO j=1,colvar%reaction_path_param%colvar_p(i)%colvar%n_atom_s
          ii=ii+1
          fi(:,ii)=colvar%reaction_path_param%colvar_p(i)%colvar%dsdr(:,j)*&
                   2.0_dp*( ss_vals(i) - ds1(i)/s1)
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
    DEALLOCATE(ds1,stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

  END SUBROUTINE dpath_colvar

! *****************************************************************************
!> \brief  distance from path calculated from the positions of a selected list of
!>         atoms as compared to the same positions in reference
!>         configurations belonging to the given path.
!> \param colvar ...
!> \param particles ...
!> \param error ...
!> \date  01.2010
!> \author MI
! *****************************************************************************
  SUBROUTINE dpath_dist_rmsd(colvar,particles,error)

    TYPE(colvar_type), POINTER               :: colvar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dpath_dist_rmsd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iat, ii, ik, natom, nconf, &
                                                rmsd_atom, stat
    INTEGER, DIMENSION(:), POINTER           :: iatom
    LOGICAL                                  :: failure
    REAL(dp)                                 :: lambda, s1, sum_exp
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: r, r0, s1v, vec_dif
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: ds1, dvec_dif, fi, riat
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ds1v
    REAL(dp), DIMENSION(:, :), POINTER       :: path_conf

    failure=.FALSE.

    nconf=colvar%reaction_path_param%nr_frames
    rmsd_atom=colvar%reaction_path_param%n_components
    lambda=colvar%reaction_path_param%lambda
    path_conf => colvar%reaction_path_param%r_ref
    iatom => colvar%reaction_path_param%i_rmsd

    natom = SIZE(particles)

    ALLOCATE(r0(3*natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(r(3*natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(riat(3,rmsd_atom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(vec_dif(rmsd_atom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(dvec_dif(3,rmsd_atom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(s1v(nconf),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1v(3,rmsd_atom,nconf),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1(3,rmsd_atom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DO i = 1,natom
      ii = (i-1)*3
      r0(ii+1) = particles(i)%r(1)
      r0(ii+2) = particles(i)%r(2)
      r0(ii+3) = particles(i)%r(3)
    END DO

    DO iat = 1,rmsd_atom
      ii = iatom(iat)
      riat(:,iat) = particles(ii)%r
    END DO

    DO ik = 1, nconf
       DO i = 1,natom
        ii = (i-1)*3
        r(ii+1) = path_conf(ii+1,ik)
        r(ii+2) = path_conf(ii+2,ik)
        r(ii+3) = path_conf(ii+3,ik)
       END DO

       CALL rmsd3(particles,r,r0,output_unit=-1,rotate=.TRUE.,error=error)

       sum_exp = 0.0_dp
       DO iat  = 1,rmsd_atom
         i = iatom(iat)
         ii = (i-1)*3
         vec_dif(iat) = (riat(1,iat)-r(ii+1))**2+(riat(2,iat)-r(ii+2))**2+(riat(3,iat)-r(ii+3))**2
         sum_exp =  sum_exp + vec_dif(iat)
         dvec_dif(1,iat) = r(ii+1)
         dvec_dif(2,iat) = r(ii+2)
         dvec_dif(3,iat) = r(ii+3)
       END DO
       s1v(ik)=EXP(-lambda*sum_exp)
       DO iat = 1,rmsd_atom
         ds1v(1,iat,ik) = dvec_dif(1,iat) * s1v(ik)
         ds1v(2,iat,ik) = dvec_dif(2,iat) * s1v(ik)
         ds1v(3,iat,ik) = dvec_dif(3,iat) * s1v(ik)
       END DO
    END DO

    s1 = accurate_sum(s1v(:))
    DO iat = 1,rmsd_atom
      ds1(1,iat) = accurate_sum( ds1v(1,iat,:))
      ds1(2,iat) = accurate_sum( ds1v(2,iat,:))
      ds1(3,iat) = accurate_sum( ds1v(3,iat,:))
    END DO
    colvar%ss=-1.0_dp/lambda*LOG(s1)

    ALLOCATE(fi(3,rmsd_atom),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    DO iat = 1,rmsd_atom
      fi(:,iat) = 2.0_dp*(riat(:,iat)-ds1(:,iat)/s1)
      CALL put_derivative(colvar,iat,fi(:,iat))
    END DO

    DEALLOCATE(fi,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(r0,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(r,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(riat,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(vec_dif,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(dvec_dif,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(s1v,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ds1v,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ds1,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
  END SUBROUTINE  dpath_dist_rmsd

! *****************************************************************************
!> \brief ...
!> \param colvar ...
!> \param particles ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dpath_rmsd(colvar,particles,error)

    TYPE(colvar_type), POINTER               :: colvar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dpath_rmsd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iat, ii, ik, natom, nconf, &
                                                rmsd_atom, stat
    INTEGER, DIMENSION(:), POINTER           :: iatom
    LOGICAL                                  :: failure
    REAL(dp)                                 :: lambda, my_rmsd, s1
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: r, r0, s1v
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: ds1, fi, riat
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ds1v
    REAL(dp), DIMENSION(:, :), POINTER       :: path_conf
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: weight
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: drmsd

    failure=.FALSE.

    nconf=colvar%reaction_path_param%nr_frames
    rmsd_atom=colvar%reaction_path_param%n_components
    lambda=colvar%reaction_path_param%lambda
    path_conf => colvar%reaction_path_param%r_ref
    iatom => colvar%reaction_path_param%i_rmsd

    natom = SIZE(particles)

    ALLOCATE(r0(3*natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(r(3*natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(riat(3,rmsd_atom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(s1v(nconf),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1v(3,rmsd_atom,nconf),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(ds1(3,rmsd_atom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    ALLOCATE(drmsd(3,natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    drmsd = 0.0_dp
    ALLOCATE (weight(natom),STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)


    DO i = 1,natom
      ii = (i-1)*3
      r0(ii+1) = particles(i)%r(1)
      r0(ii+2) = particles(i)%r(2)
      r0(ii+3) = particles(i)%r(3)
    END DO

    DO iat = 1,rmsd_atom
      ii = iatom(iat)
      riat(:,iat) = particles(ii)%r
    END DO

! set weights of atoms in the rmsd list
    weight = 0.0_dp
    DO iat  = 1,rmsd_atom
       i = iatom(iat)
       weight(i) = 1.0_dp
    END DO

    DO ik = 1, nconf
       DO i = 1,natom
        ii = (i-1)*3
        r(ii+1) = path_conf(ii+1,ik)
        r(ii+2) = path_conf(ii+2,ik)
        r(ii+3) = path_conf(ii+3,ik)
       END DO

       CALL rmsd3(particles,r0,r,output_unit=-1,weights=weight,my_val=my_rmsd,&
             rotate=.FALSE.,drmsd3=drmsd,error=error)

       s1v(ik)=EXP(-lambda*my_rmsd)
       DO iat = 1,rmsd_atom
         i = iatom(iat)
         ds1v(1,iat,ik) = drmsd(1,i) * s1v(ik)
         ds1v(2,iat,ik) = drmsd(2,i) * s1v(ik)
         ds1v(3,iat,ik) = drmsd(3,i) * s1v(ik)
       END DO
    END DO

    s1 = accurate_sum(s1v(:))
    DO iat = 1,rmsd_atom
      ds1(1,iat) = accurate_sum( ds1v(1,iat,:))
      ds1(2,iat) = accurate_sum( ds1v(2,iat,:))
      ds1(3,iat) = accurate_sum( ds1v(3,iat,:))
    END DO
    colvar%ss=-1.0_dp/lambda*LOG(s1)

    ALLOCATE(fi(3,rmsd_atom),stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

    DO iat = 1,rmsd_atom
      fi(:,iat) = ds1(:,iat)/s1
      CALL put_derivative(colvar,iat,fi(:,iat))
    END DO

    DEALLOCATE(fi,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(r0,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(r,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(riat,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(s1v,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ds1v,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(ds1,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(drmsd,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
    DEALLOCATE(weight,STAT=stat)
    CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)

  END SUBROUTINE dpath_rmsd

! *****************************************************************************
!> \brief evaluates the force due to population colvar
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
!> \date  01.2009
!> \author fsterpone
! *****************************************************************************
  SUBROUTINE population_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'population_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, jj, n_atoms_from, &
                                                n_atoms_to, ndcrd, nncrd
    LOGICAL                                  :: failure
    REAL(dp) :: dfunc, dfunc_coord, ftmp(3), func, func_coord, &
      inv_n_atoms_from, invden, n_0, ncoord, norm, num, population, r12, r_0, &
      rdist, sigma, ss(3), xij(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: ftmp_coord
    REAL(dp), DIMENSION(3)                   :: xpi, xpj
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.
    ! If we defined the coordination number with KINDS then we have still
    ! to fill few missing informations...
    NULLIFY(particles_i)
    CPPrecondition(colvar%type_id==population_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    n_atoms_to=colvar%population_param%n_atoms_to
    n_atoms_from=colvar%population_param%n_atoms_from
    nncrd=colvar%population_param%nncrd
    ndcrd=colvar%population_param%ndcrd
    r_0=colvar%population_param%r_0
    n_0=colvar%population_param%n0
    sigma=colvar%population_param%sigma

    ALLOCATE(ftmp_coord(3,n_atoms_to))
    ftmp_coord=0.0_dp

    ncoord=0.0_dp
    population=0.0_dp

    colvar%dsdr=0.0_dp
    inv_n_atoms_from=1.0_dp/REAL(n_atoms_from)

    norm=SQRT(pi*2.0_dp)*sigma
    norm=1/norm

    DO ii=1,n_atoms_from
       i=colvar%population_param%i_at_from(ii)
       CALL get_coordinates(colvar, i, xpi, my_particles)
       DO jj=1,n_atoms_to
          i=colvar%population_param%i_at_to(jj)
          CALL get_coordinates(colvar, i, xpj, my_particles)
          ss=MATMUL(cell%h_inv,xpi(:)-xpj(:))
          ss=ss-NINT(ss)
          xij=MATMUL(cell%hmat,ss)
          r12=SQRT(xij(1)**2+xij(2)**2+xij(3)**2)
          IF(r12 < 1.0e-8_dp)CYCLE
          rdist = r12/r_0
          num=(1.0_dp-rdist**nncrd)
          invden=1.0_dp/(1.0_dp-rdist**ndcrd)
          func_coord=num*invden
          dfunc_coord= (- nncrd*rdist **(nncrd-1) *invden  &
               + num*(invden)**2 * ndcrd *rdist **(ndcrd-1))/(r12*r_0)

          ncoord=ncoord+func_coord
          ftmp_coord(1,jj) = dfunc_coord*xij(1)
          ftmp_coord(2,jj) = dfunc_coord*xij(2)
          ftmp_coord(3,jj) = dfunc_coord*xij(3)
       END DO

       func=EXP(-(ncoord-n_0)**2/(2.0_dp*sigma*sigma))
       dfunc=-func*(ncoord-n_0)/(sigma*sigma)

       population=population+norm*func
       DO jj=1,n_atoms_to
          ftmp(1)=ftmp_coord(1,jj)*dfunc
          ftmp(2)=ftmp_coord(2,jj)*dfunc
          ftmp(3)=ftmp_coord(3,jj)*dfunc
          CALL put_derivative(colvar, ii, ftmp)
          ftmp(1)=-ftmp_coord(1,jj)*dfunc
          ftmp(2)=-ftmp_coord(2,jj)*dfunc
          ftmp(3)=-ftmp_coord(3,jj)*dfunc
          CALL put_derivative(colvar, n_atoms_from+jj, ftmp)
       ENDDO
       ncoord=0.0_dp
    ENDDO
    colvar%ss=population
  END SUBROUTINE population_colvar

! *****************************************************************************
!> \brief evaluates the force due to the gyration radius colvar
!>        sum_i (r_i-rcom)^2/N
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
!> \date  03.2009
!> \author MI
! *****************************************************************************
  SUBROUTINE gyration_radius_colvar(colvar,cell,subsys,particles,error)

    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gyration_radius_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, n_atoms
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dri2, func, gyration, inv_n, &
                                                mass_tot, mi
    REAL(dp), DIMENSION(3)                   :: dfunc, dxi, ftmp, ss, xpcom, &
                                                xpi
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.

    NULLIFY(particles_i,my_particles)
    CPPrecondition(colvar%type_id==gyration_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    n_atoms=colvar%gyration_param%n_atoms
    inv_n = 1.0_dp/n_atoms


    !compute COM position
    xpcom = 0.0_dp
    mass_tot = 0.0_dp
    DO ii = 1,n_atoms
       i=colvar%gyration_param%i_at(ii)
       CALL get_coordinates(colvar, i, xpi, my_particles)
       CALL get_mass(colvar, i, mi, my_particles)
        xpcom(:) = xpcom(:) + xpi(:)*mi
       mass_tot = mass_tot + mi
    END DO
    xpcom(:) = xpcom(:)/mass_tot

    func  = 0.0_dp
    ftmp  = 0.0_dp
    dfunc = 0.0_dp
    DO ii = 1,n_atoms
       i=colvar%gyration_param%i_at(ii)
       CALL get_coordinates(colvar, i, xpi, my_particles)
       ss=MATMUL(cell%h_inv,xpi(:)-xpcom(:))
       ss=ss-NINT(ss)
       dxi=MATMUL(cell%hmat,ss)
       dri2=(dxi(1)**2+dxi(2)**2+dxi(3)**2)
       func = func + dri2
       dfunc(:) = dfunc(:) + dxi(:)
    END DO
    gyration = SQRT(inv_n*func)

    DO ii = 1,n_atoms
      i=colvar%gyration_param%i_at(ii)
      CALL get_coordinates(colvar, i, xpi, my_particles)
      CALL get_mass(colvar, i, mi, my_particles)
      ss=MATMUL(cell%h_inv,xpi(:)-xpcom(:))
      ss=ss-NINT(ss)
      dxi=MATMUL(cell%hmat,ss)
      ftmp(1) = dxi(1) - dfunc(1)*mi/mass_tot
      ftmp(2) = dxi(2) - dfunc(2)*mi/mass_tot
      ftmp(3) = dxi(3) - dfunc(3)*mi/mass_tot
      ftmp(:) = ftmp(:)*inv_n/gyration
      CALL put_derivative(colvar, ii, ftmp)
    END DO
    colvar%ss=gyration

  END SUBROUTINE gyration_radius_colvar

! *****************************************************************************
!> \brief evaluates the force due to the rmsd colvar
!> \param colvar ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
!> \date  12.2009
!> \author MI
!> \note  could be extended to be used with more than 2 reference structures
! *****************************************************************************
  SUBROUTINE rmsd_colvar(colvar,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rmsd_colvar', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL rmsd_colvar_low(colvar,subsys,particles,error)
  END SUBROUTINE rmsd_colvar

! *****************************************************************************
!> \brief  evaluates the force due to the rmsd colvar
!>        ss = (RMSDA-RMSDB)/(RMSDA+RMSDB)
!>        RMSD is calculated with respect to two reference structures, A and B,
!>        considering all the atoms of the system or only a subset of them,
!>        as selected by the input keyword LIST
!> \param colvar ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
!> \date  12.2009
!> \par History TL 2012 (generalized to any number of frames) 
!> \author MI
! *****************************************************************************
  SUBROUTINE rmsd_colvar_low(colvar,subsys,particles,error)

    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rmsd_colvar_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, natom, nframes, stat
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: cv_val, f1, ftmp(3)
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: der, r, rmsd
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r0
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: drmsd
    REAL(kind=dp), DIMENSION(:), POINTER     :: weights
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.

    NULLIFY(my_particles,particles_i,weights)
    CPPrecondition(colvar%type_id==rmsd_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF

    natom = SIZE(my_particles)
    nframes=colvar%rmsd_param%nr_frames
    ALLOCATE (drmsd(3,natom,nframes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    drmsd = 0.0_dp

    ALLOCATE (r0(3*natom,nframes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (rmsd(nframes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (der(nframes),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (r(3*natom),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    weights => colvar%rmsd_param%weights
    DO i = 1,natom
      ii = (i-1)*3
      r(ii+1) = my_particles(i)%r(1)
      r(ii+2) = my_particles(i)%r(2)
      r(ii+3) = my_particles(i)%r(3)
    END DO
    r0(:,:) = colvar%rmsd_param%r_ref
    rmsd = 0.0_dp

    CALL rmsd3(my_particles,r,r0(:,1),output_unit=-1,weights=weights,my_val=rmsd(1),rotate=.FALSE.,drmsd3=drmsd(:,:,1),error=error)

    IF (nframes==2) THEN
       CALL rmsd3(my_particles,r,r0(:,2),output_unit=-1,weights=weights,my_val=rmsd(2),rotate=.FALSE.,drmsd3=drmsd(:,:,2),&
                  error=error)

       f1 = 1.0_dp/(rmsd(1)+rmsd(2))
       ! (rmsdA-rmsdB)/(rmsdA+rmsdB)
       cv_val = (rmsd(1)-rmsd(2))*f1
       ! (rmsdA+rmsdB)^-1-(rmsdA-rmsdB)/(rmsdA+rmsdB)^2
       der(1) = f1-cv_val*f1
       ! -(rmsdA+rmsdB)^-1-(rmsdA-rmsdB)/(rmsdA+rmsdB)^2
       der(2) = -f1-cv_val*f1

       DO i = 1,colvar%rmsd_param%n_atoms
         ii = colvar%rmsd_param%i_rmsd(i)
         IF(weights(ii)>0.0_dp) THEN
           ftmp(1) = der(1)*drmsd(1,ii,1)+der(2)*drmsd(1,ii,2)
           ftmp(2) = der(1)*drmsd(2,ii,1)+der(2)*drmsd(2,ii,2)
           ftmp(3) = der(1)*drmsd(3,ii,1)+der(2)*drmsd(3,ii,2)
           CALL put_derivative(colvar, i, ftmp)
         END IF
       END DO
    ELSE IF (nframes==1) THEN
       ! Protect in case of numerical issues (for two identical frames!)
       rmsd(1)=ABS(rmsd(1))
       cv_val = SQRT(rmsd(1))
       f1     = 0.0_dp
       IF (cv_val /= 0.0_dp) f1 = 0.5_dp/cv_val
       DO i = 1,colvar%rmsd_param%n_atoms
         ii = colvar%rmsd_param%i_rmsd(i)
         IF(weights(ii)>0.0_dp) THEN
           ftmp(1) = f1*drmsd(1,ii,1)
           ftmp(2) = f1*drmsd(2,ii,1)
           ftmp(3) = f1*drmsd(3,ii,1)
           CALL put_derivative(colvar, i, ftmp)
         END IF
       END DO
    ELSE
        CALL cp_unimplemented_error(fromWhere=routineP, &
             message="RMSD implemented only for 1 and 2 reference frames!", error=error, error_level=cp_failure_level)
    END IF
    colvar%ss=cv_val

    DEALLOCATE(der,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(r0,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(r,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(drmsd,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(rmsd,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE rmsd_colvar_low

! *****************************************************************************
!> \brief evaluates the force from ring puckering collective variables
!>   Cramer and Pople, JACS 97 1354 (1975)
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param error ...
!> \date 08.2012
!> \author JGH
! *****************************************************************************
  SUBROUTINE  ring_puckering_colvar(colvar,cell,subsys,particles,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ring_puckering_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, istat, j, jj, m, nring
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, at, b, da, db, ds, kr, &
                                                rpxpp, svar
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cosj, sinj, z
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: nforce, zforce
    REAL(KIND=dp), DIMENSION(3)              :: ftmp, nv, r0, rp, rpp, uv
    REAL(KIND=dp), DIMENSION(3, 3)           :: dnvp, dnvpp
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles

    failure=.FALSE.

    CPPrecondition(colvar%type_id==ring_puckering_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF

    nring = colvar%ring_puckering_param%nring 
    ALLOCATE (r(3,nring),z(nring),cosj(nring),sinj(nring),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (nforce(3,3,nring),zforce(nring,nring,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ii = 1,nring
       i=colvar%ring_puckering_param%atoms(ii)
       CALL get_coordinates(colvar, i, r(:,ii), my_particles)
    END DO
    ! get all atoms within PBC distance of atom 1
    r0(:) = r(:,1)
    DO ii = 1,nring
       r(:,ii) = pbc(r(:,ii),r0,cell)
    END DO
    !compute origin position
    r0 = 0.0_dp
    DO ii = 1,nring
       r0(:) = r0(:) + r(:,ii)
    END DO
    kr = 1._dp/REAL(nring,KIND=dp)
    r0(:) = r0(:)*kr
    DO ii = 1,nring
       r(:,ii) = r(:,ii) - r0(:)
    END DO
    ! orientation vectors
    rp = 0._dp
    rpp = 0._dp
    DO ii = 1,nring
       cosj(ii) = COS(twopi*(ii-1)*kr)
       sinj(ii) = SIN(twopi*(ii-1)*kr)
       rp(:) = rp(:) + r(:,ii)*sinj(ii)
       rpp(:) = rpp(:) + r(:,ii)*cosj(ii)
    END DO
    nv = vector_product(rp,rpp)
    nv = nv/SQRT(SUM(nv**2))

    ! derivatives of normal
    uv = vector_product(rp,rpp)
    rpxpp = SQRT(SUM(uv**2))
    DO i=1,3
      uv = 0._dp
      uv(i) = 1._dp
      uv = vector_product(uv,rpp)/rpxpp
      dnvp(:,i) = uv - nv * SUM(uv*nv)
      uv = 0._dp
      uv(i) = 1._dp
      uv = vector_product(rp,uv)/rpxpp
      dnvpp(:,i) = uv - nv * SUM(uv*nv)
    END DO
    DO ii = 1,nring
      nforce(:,:,ii) = dnvp(:,:)*sinj(ii) + dnvpp(:,:)*cosj(ii)
    END DO

    ! molecular z-coordinate
    DO ii = 1,nring
       z(ii) = SUM(r(:,ii)*nv(:))
    END DO
    ! z-force
    DO ii = 1,nring
      DO jj = 1,nring
         IF(ii==jj) THEN
            zforce(ii,jj,:) = nv
         ELSE
            zforce(ii,jj,:) = 0._dp
         END IF
         DO i=1,3
           DO j=1,3
             zforce(ii,jj,i) = zforce(ii,jj,i) + r(j,ii)*nforce(j,i,jj)
           END DO
         END DO
      END DO
    END DO

    IF(colvar%ring_puckering_param%iq == 0) THEN
       ! total puckering amplitude
       svar = SQRT(SUM(z**2))
       DO ii = 1,nring
          ftmp = 0._dp
          DO jj = 1,nring
             ftmp(:) = ftmp(:) + zforce(jj,ii,:)*z(jj)
          END DO
          ftmp = ftmp/svar
          CALL put_derivative(colvar, ii, ftmp)
       END DO
    ELSE
       m = ABS(colvar%ring_puckering_param%iq)
       CPPrecondition(m/=1,cp_failure_level,routineP,error,failure)
       IF(MOD(nring,2)==0 .AND. colvar%ring_puckering_param%iq == nring/2) THEN
          ! single puckering amplitude
          svar = 0._dp
          DO ii = 1,nring
             IF(MOD(ii,2)==0) THEN
                svar = svar - z(ii)
             ELSE
                svar = svar + z(ii)
             END IF
          END DO
          svar = svar*SQRT(kr)
          DO ii = 1,nring
             ftmp = 0._dp
             DO jj=1,nring
                IF(MOD(jj,2)==0) THEN
                  ftmp(:) = ftmp(:) - zforce(jj,ii,:)*SQRT(kr)
                ELSE
                  ftmp(:) = ftmp(:) + zforce(jj,ii,:)*SQRT(kr)
                END IF
             END DO
             CALL put_derivative(colvar, ii, -ftmp)
          END DO
       ELSE
          CPPrecondition(m<=(nring-1)/2,cp_failure_level,routineP,error,failure)
          a = 0._dp
          b = 0._dp
          DO ii = 1,nring
             a = a + z(ii)*COS(twopi*m*(ii-1)*kr)
             b = b - z(ii)*SIN(twopi*m*(ii-1)*kr)
          END DO
          a = a*SQRT(2._dp*kr)
          b = b*SQRT(2._dp*kr)
          IF(colvar%ring_puckering_param%iq > 0) THEN
             ! puckering amplitude
             svar = SQRT(a*a+b*b)
             da = a/svar 
             db = b/svar 
          ELSE 
             ! puckering phase angle
             at = ATAN2(a,b)
             IF(at>pi/2._dp) THEN
               svar = 2.5_dp*pi - at
             ELSE
               svar = 0.5_dp*pi - at
             END IF
             da = -b/(a*a+b*b)
             db = a/(a*a+b*b)
          END IF
          DO jj = 1,nring
             ftmp = 0._dp
             DO ii = 1,nring
                ds = da*COS(twopi*m*(ii-1)*kr)
                ds = ds - db*SIN(twopi*m*(ii-1)*kr)
                ftmp(:) = ftmp(:) + ds*SQRT(2._dp*kr)*zforce(ii,jj,:)
             END DO
             CALL put_derivative(colvar, jj, ftmp)
          END DO
       END IF
    END IF

    colvar%ss=svar

    DEALLOCATE (r,z,cosj,sinj,nforce,zforce,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE  ring_puckering_colvar

! *****************************************************************************
!> \brief used to print reaction_path function values on an arbitrary dimensional grid
!> \param iw1 ...
!> \param ncol ...
!> \param f_vals ...
!> \param v_count ...
!> \param gp ...
!> \param grid_sp ...
!> \param step_size ...
!> \param istart ...
!> \param iend ...
!> \param s1v ...
!> \param s1 ...
!> \param p_bounds ...
!> \param lambda ...
!> \param ifunc ...
!> \param nconf ...
!> \retval k ...
!> \author fschiff
! *****************************************************************************
  RECURSIVE FUNCTION rec_eval_grid(iw1,ncol,f_vals,v_count,&
       gp,grid_sp,step_size,istart,iend,s1v,s1,p_bounds,lambda,ifunc,nconf) RESULT(k)
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
    INTEGER                                  :: ifunc, nconf, k

    INTEGER                                  :: count1, i

    k=1
    IF(v_count.LT.ncol)THEN
       count1=v_count+1
       DO i=p_bounds(1,count1),p_bounds(2,count1)
          gp(count1)=REAL(i,KIND=dp)*grid_sp(count1)
          k=rec_eval_grid(iw1,ncol,f_vals,count1,gp,grid_sp,step_size,istart,iend,s1v,s1,p_bounds,lambda,ifunc,nconf)
       END DO
    ELSE IF(v_count==ncol .AND. ifunc==1)THEN
       DO i=istart,iend
          s1v(1,i)=REAL(i,kind=dp)*step_size*EXP(-lambda*DOT_PRODUCT(gp(:)-f_vals(:,i),&
               gp(:)-f_vals(:,i)))
          s1v(2,i)=EXP(-lambda*DOT_PRODUCT(gp(:)-f_vals(:,i),gp(:)-f_vals(:,i)))
       END DO
       DO i=1,2
          s1(i)=accurate_sum(s1v(i,:))
       END DO
       WRITE(iw1, '(5F10.5)')gp(:),s1(1)/s1(2)/REAL(nconf-1,dp)
    ELSE IF(v_count==ncol .AND. ifunc==2)THEN
       DO i=istart,iend
          s1v(1,i)=EXP(-lambda*DOT_PRODUCT(gp(:)-f_vals(:,i),gp(:)-f_vals(:,i)))
       END DO
       s1(1)=accurate_sum(s1v(1,:))

       WRITE(iw1, '(5F10.5)')gp(:),-lambda*LOG(s1(1))
    END IF
  END FUNCTION rec_eval_grid

! *****************************************************************************
!> \brief  Reads the coordinates of reference configurations given in input
!>         either as xyz files or in &COORD section
!> \param frame_section ...
!> \param para_env ...
!> \param nr_frames ...
!> \param r_ref ...
!> \param n_atoms ...
!> \param error ...
!> \date 01.2010
!> \author MI
! *****************************************************************************
  SUBROUTINE read_frames(frame_section,para_env,nr_frames,r_ref,n_atoms,error)

    TYPE(section_vals_type), POINTER         :: frame_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: nr_frames
    REAL(dp), DIMENSION(:, :), POINTER       :: r_ref
    INTEGER, INTENT(OUT)                     :: n_atoms
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_frames', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: dummy_char
    INTEGER                                  :: i, j, natom, stat
    LOGICAL                                  :: explicit, failure, my_end
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rptr
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(section_vals_type), POINTER         :: coord_section

     failure = .FALSE.
     NULLIFY(rptr)

     DO i = 1,nr_frames
       coord_section => section_vals_get_subs_vals(frame_section,"COORD",i_rep_section=i,error=error)
       CALL section_vals_get(coord_section,explicit=explicit,error=error)
       ! Cartesian Coordinates
       IF (explicit) THEN
          CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
                n_rep_val=natom,error=error)
          IF(i==1) THEN
            ALLOCATE(r_ref(3*natom,nr_frames), STAT=stat)
            CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
            n_atoms = natom
          ELSE
           CPPostcondition(3*natom==SIZE(r_ref,1),cp_failure_level,routineP,error,failure)
          END IF
          DO j = 1, natom
              CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
                     i_rep_val=j,r_vals=rptr,error=error)
               r_ref((j-1)*3+1:(j-1)*3+3,i) = rptr(1:3)
          END DO ! natom
       ELSE
          CALL section_vals_val_get(frame_section,"COORD_FILE_NAME",i_rep_section=i,c_val=filename,error=error)
          CPPostcondition(TRIM(filename)/="",cp_failure_level,routineP,error,failure)
          NULLIFY(parser)
          ALLOCATE(rptr(3))
          CALL parser_create(parser,filename,para_env=para_env,parse_white_lines=.TRUE.,&
                              error=error)
          CALL parser_get_next_line(parser,1,error=error)
          ! Start parser
          CALL parser_get_object(parser,natom,error=error)
          CALL parser_get_next_line(parser,1,error=error)
          IF(i==1) THEN
             ALLOCATE(r_ref(3*natom,nr_frames), STAT=stat)
             CPPostcondition(stat == 0,cp_failure_level,routinep,error,failure)
             n_atoms = natom
          ELSE
           CPPostcondition(3*natom==SIZE(r_ref,1),cp_failure_level,routineP,error,failure)
          END IF
          DO j = 1, natom
            ! Atom coordinates
            CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
            CALL cp_assert(.NOT.my_end,cp_fatal_level,cp_assertion_failed,routineP,&
                 "Number of lines in XYZ format not equal to the number of atoms."//&
                 " Error in XYZ format for COORD_A (CV rmsd). Very probably the"//&
                 " line with title is missing or is empty. Please check the XYZ file and rerun your job!"//&
CPSourceFileRef,&
            only_ionode=.TRUE.)
            READ(parser%input_line,*) dummy_char,rptr(1:3)
            r_ref((j-1)*3+1,i) = cp_unit_to_cp2k(rptr(1),"angstrom",error=error)
            r_ref((j-1)*3+2,i) = cp_unit_to_cp2k(rptr(2),"angstrom",error=error)
            r_ref((j-1)*3+3,i) = cp_unit_to_cp2k(rptr(3),"angstrom",error=error)
          END DO ! natom
          CALL parser_release(parser,error=error)
          DEALLOCATE(rptr)
       END IF
     END DO  ! nr_frames

  END SUBROUTINE  read_frames


! *****************************************************************************
!> \brief evaluates the collective variable associated with a hydrogen bond
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param qs_env ...
!> \param error ...
!> \author alin m elena
! *****************************************************************************
  SUBROUTINE Wc_colvar(colvar,cell,subsys,particles,qs_env,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(qs_environment_type), POINTER, OPTIONAL       :: qs_env ! optional just because I am lazy... but I should get rid of it...
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'Wc_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: Od, H, Oa
    LOGICAL                                  :: failure
    REAL(dp)                                 :: rOd(3), rOa(3), rH(3), &
                                                x,y,s(3),xv(3),dmin,amin
    INTEGER                                  :: idmin, iamin,i,j
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles
    TYPE(wannier_centres_type), DIMENSION(:),POINTER :: wc
    INTEGER, ALLOCATABLE                     :: wcai(:),wcdi(:) ! contains the indeces of the wannier centres closed to the donor and acceptor
    INTEGER                                  :: nwca,nwcd,ierror
    REAL(dp)                                 :: rcut

    failure=.FALSE.
    NULLIFY(particles_i,wc)

    CPPrecondition(colvar%type_id==Wc_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    CALL get_qs_env(qs_env,WannierCentres=wc,error=error)
    rcut =colvar%Wc%rcut ! distances are in bohr as far as I remember
    Od=colvar%Wc%ids(1)
    H=colvar%Wc%ids(2)
    Oa=colvar%Wc%ids(3)
    CALL get_coordinates(colvar, Od, rOd, my_particles)
    CALL get_coordinates(colvar, H, rH, my_particles)
    CALL get_coordinates(colvar, Oa, rOa, my_particles)
    ALLOCATE(wcai(SIZE(wc(1)%WannierHamDiag)),stat=ierror)
    CPPrecondition(ierror==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wcdi(SIZE(wc(1)%WannierHamDiag)),stat=ierror)
    CPPrecondition(ierror==0,cp_failure_level,routineP,error,failure)
    nwca=0
    nwcd=0
    DO j=1,SIZE(wc(1)%WannierHamDiag)
      x=distance(rOd-wc(1)%centres(:,j))
      y=distance(rOa-wc(1)%centres(:,j))
      IF (x<rcut) THEN
        nwcd=nwcd+1
        wcdi(nwcd)=j
        CYCLE
      ENDIF
      IF (y<rcut) THEN
        nwca=nwca+1
        wcai(nwca)=j
      ENDIF
    ENDDO

    dmin=distance(rH-wc(1)%centres(:,wcdi(1)))
    amin=distance(rH-wc(1)%centres(:,wcai(1)))
    idmin=wcdi(1)
    iamin=wcai(1)
    !dmin constains the smallest numer, amin the next smallest
    DO i=2,nwcd
      x=distance(rH-wc(1)%centres(:,wcdi(i)))
      IF (x<dmin) THEN
        dmin=x
        idmin=wcdi(i)
      ENDIF
    ENDDO
    DO i=2,nwca
      x=distance(rH-wc(1)%centres(:,wcai(i)))
      IF (x<amin) THEN
        amin=x
        iamin=wcai(i)
      ENDIF
    ENDDO
!     zero=0.0_dp
!     CALL put_derivative(colvar, 1, zero)
!     CALL put_derivative(colvar, 2,zero)
!     CALL put_derivative(colvar, 3, zero)

!     write(*,'(2(i0,1x),4(f16.8,1x))')idmin,iamin,wc(1)%WannierHamDiag(idmin),wc(1)%WannierHamDiag(iamin),dmin,amin
    colvar%ss = wc(1)%WannierHamDiag(idmin)-wc(1)%WannierHamDiag(iamin)
    DEALLOCATE(wcai,stat=ierror)
    CPPrecondition(ierror==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wcdi,stat=ierror)
    CPPrecondition(ierror==0,cp_failure_level,routineP,error,failure)

    CONTAINS
! *****************************************************************************
!> \brief ...
!> \param rij ...
! *****************************************************************************
    REAL(dp) FUNCTION distance(rij)
    REAL(dp), INTENT(in) :: rij(3)

      s=MATMUL(cell%h_inv,rij)
      s=s-NINT(s)
      xv=MATMUL(cell%hmat,s)
      distance=SQRT(DOT_PRODUCT(xv,xv))
    END FUNCTION distance

  END SUBROUTINE Wc_colvar

! *****************************************************************************
!> \brief evaluates the collective variable associated with a hydrogen bond wire
!> \param colvar ...
!> \param cell ...
!> \param subsys ...
!> \param particles ...
!> \param qs_env ...
!> \param error ...
!> \author alin m elena
! *****************************************************************************
  SUBROUTINE HBP_colvar(colvar,cell,subsys,particles,qs_env,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles
    TYPE(qs_environment_type), POINTER, OPTIONAL       :: qs_env ! optional just because I am lazy... but I should get rid of it...
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'HBP_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: Od, H, Oa
    LOGICAL                                  :: failure
    REAL(dp)                                 :: rOd(3), rOa(3), rH(3), &
                                                x,y,s(3),xv(3),dmin,amin
    INTEGER                                  :: idmin, iamin,i,j,il,output_unit
    TYPE(particle_list_type), POINTER        :: particles_i
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: my_particles
    TYPE(wannier_centres_type),&
                        DIMENSION(:),POINTER :: wc
    INTEGER, ALLOCATABLE                     :: wcai(:),wcdi(:) ! contains the indeces of the wannier centres closed to the donor and acceptor
    INTEGER                                  :: nwca,nwcd,ierror
    REAL(dp)                                 :: rcut
!     real(dp)                                 :: zero(3)
    TYPE(cp_logger_type), POINTER            :: logger



    failure=.FALSE.
    NULLIFY(particles_i,wc,logger)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CPPrecondition(colvar%type_id==HBP_colvar_id,cp_failure_level,routineP,error,failure)
    IF (PRESENT(particles)) THEN
       my_particles => particles
    ELSE
       CPPrecondition(PRESENT(subsys),cp_failure_level,routineP,error,failure)
       CALL cp_subsys_get(subsys,particles=particles_i,error=error)
       my_particles => particles_i%els
    END IF
    CALL get_qs_env(qs_env,WannierCentres=wc,error=error)
    rcut = colvar%HBP%rcut ! distances are in bohr as far as I remember
    ALLOCATE(wcai(SIZE(wc(1)%WannierHamDiag)),stat=ierror)
    CPPrecondition(ierror==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(wcdi(SIZE(wc(1)%WannierHamDiag)),stat=ierror)
    CPPrecondition(ierror==0,cp_failure_level,routineP,error,failure)
    colvar%ss=0.0_dp
    DO il=1,colvar%HBP%nPoints
      Od=colvar%HBP%ids(il,1)
      H=colvar%HBP%ids(il,2)
      Oa=colvar%HBP%ids(il,3)
      CALL get_coordinates(colvar, Od, rOd, my_particles)
      CALL get_coordinates(colvar, H, rH, my_particles)
      CALL get_coordinates(colvar, Oa, rOa, my_particles)
      nwca=0
      nwcd=0
      DO j=1,SIZE(wc(1)%WannierHamDiag)
        x=distance(rOd-wc(1)%centres(:,j))
        y=distance(rOa-wc(1)%centres(:,j))
        IF (x<rcut) THEN
          nwcd=nwcd+1
          wcdi(nwcd)=j
          CYCLE
        ENDIF
        IF (y<rcut) THEN
          nwca=nwca+1
          wcai(nwca)=j
        ENDIF
      ENDDO

      dmin=distance(rH-wc(1)%centres(:,wcdi(1)))
      amin=distance(rH-wc(1)%centres(:,wcai(1)))
      idmin=wcdi(1)
      iamin=wcai(1)
      !dmin constains the smallest numer, amin the next smallest
      DO i=2,nwcd
        x=distance(rH-wc(1)%centres(:,wcdi(i)))
        IF (x<dmin) THEN
          dmin=x
          idmin=wcdi(i)
        ENDIF
      ENDDO
      DO i=2,nwca
        x=distance(rH-wc(1)%centres(:,wcai(i)))
        IF (x<amin) THEN
          amin=x
          iamin=wcai(i)
        ENDIF
      ENDDO
!       if (output_unit>0) then
!          write(output_unit,'(a,2(i0,1x),5(f16.8,1x))')"HBP|",idmin,iamin,wc(1)%WannierHamDiag(idmin),wc(1)%WannierHamDiag(iamin),dmin,amin, wc(1)%WannierHamDiag(idmin)-wc(1)%WannierHamDiag(iamin)
!       endif
      colvar%HBP%ewc(il) = colvar%HBP%shift + wc(1)%WannierHamDiag(idmin)-wc(1)%WannierHamDiag(iamin)
      colvar%ss = colvar%ss + colvar%HBP%shift + wc(1)%WannierHamDiag(idmin)-wc(1)%WannierHamDiag(iamin)
    ENDDO
     IF (output_unit>0) THEN
       DO il=1,colvar%HBP%nPoints
         WRITE(output_unit,'(a,1(f16.8,1x))')"HBP| = ",colvar%HBP%ewc(il)
       ENDDO
       WRITE(output_unit,'(a,1(f16.8,1x))')"HBP|\theta(x) = ",colvar%ss
    ENDIF
    DEALLOCATE(wcai,stat=ierror)
    CPPrecondition(ierror==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(wcdi,stat=ierror)
    CPPrecondition(ierror==0,cp_failure_level,routineP,error,failure)

  CONTAINS
! *****************************************************************************
!> \brief ...
!> \param rij ...
! *****************************************************************************
    REAL(dp) FUNCTION distance(rij)
    REAL(dp), INTENT(in) :: rij(3)

      s=MATMUL(cell%h_inv,rij)
      s=s-NINT(s)
      xv=MATMUL(cell%hmat,s)
      distance=SQRT(DOT_PRODUCT(xv,xv))
    END FUNCTION distance

  END SUBROUTINE HBP_colvar


END MODULE colvar_methods
