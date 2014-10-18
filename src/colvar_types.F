!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initialize the collective variables types
!> \par History
!>      5.2004 created [fawzi and alessandro]
!>      1.2009 Fabio Sterpone : added the population COLVAR
!> \author Teodoro Laino
! *****************************************************************************
MODULE colvar_types
  USE f77_blas
  USE input_constants,                 ONLY: do_clv_fix_point,&
                                             do_clv_geo_center,&
                                             do_clv_xyz,&
                                             plane_def_atoms,&
                                             plane_def_vec
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE particle_types,                  ONLY: particle_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colvar_types'

  PUBLIC :: colvar_type,&
            colvar_p_type,&
            colvar_p_reallocate,&
            colvar_p_release,&
            colvar_create,&
            colvar_clone,&
            colvar_setup,&
            colvar_release,&
            colvar_counters,&
            colvar_check_points,&
            eval_point_der,&
            eval_point_pos,&
            eval_point_mass,&
            diff_colvar

  INTEGER, PARAMETER, PUBLIC :: no_colvar_id                  = -2,&
                                dist_colvar_id                =  1,&
                                coord_colvar_id               =  2,&
                                torsion_colvar_id             =  3,&
                                angle_colvar_id               =  4,&
                                plane_distance_colvar_id      =  5,&
                                rotation_colvar_id            =  6,&
                                dfunct_colvar_id              =  7,&
                                qparm_colvar_id               =  8,&
                                hydronium_colvar_id           =  9,&
                                reaction_path_colvar_id       = 10,&
                                combine_colvar_id             = 11,&
                                population_colvar_id          = 12,&
                                plane_plane_angle_colvar_id   = 13,&
                                gyration_colvar_id            = 14,&
                                rmsd_colvar_id                = 15,&
                                distance_from_path_colvar_id  = 16,&
                                xyz_diag_colvar_id            = 17,&
                                xyz_outerdiag_colvar_id       = 18,&
                                u_colvar_id                   = 19,&
                                Wc_colvar_id                  = 20,&
                                hbp_colvar_id                 = 21,&
                                ring_puckering_colvar_id      = 22,&
                                mindist_colvar_id             = 23


! *****************************************************************************
!> \brief parameters for the distance collective variable
!> \param i_at ,j_at: indexes of the two atoms between which you calculate
!>        the distance
!> \author alessandro laio and fawzi mohamed
! *****************************************************************************
  TYPE dist_colvar_type
     INTEGER       :: i_at, j_at, axis_id
  END TYPE dist_colvar_type

! *****************************************************************************
  TYPE coord_colvar_type
          LOGICAL                    :: do_chain, use_kinds_from, use_kinds_to,&
                                        use_kinds_to_b
     INTEGER                         :: n_atoms_to,&
                                        n_atoms_from,&
                                        nncrd,&
                                        ndcrd,&
                                        n_atoms_to_b,&
                                        nncrd_b,&
                                        ndcrd_b
     INTEGER, POINTER, DIMENSION (:) :: i_at_from,&
                                        i_at_to,&
                                        i_at_to_b
     CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: c_kinds_from,&
                                                                    c_kinds_to,&
                                                                    c_kinds_to_b
     REAL(KIND=dp)                   :: r_0, r_0_b
  END TYPE coord_colvar_type

! *****************************************************************************
 TYPE population_colvar_type
     LOGICAL                         :: use_kinds_from, use_kinds_to
     INTEGER                         :: n_atoms_to,&
                                        n_atoms_from,&
                                        nncrd,&
                                        ndcrd,&
                                        n0
     INTEGER, POINTER, DIMENSION (:) :: i_at_from,&
                                        i_at_to
     CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: c_kinds_from,&
                                                                    c_kinds_to
     REAL(KIND=dp)                   :: r_0,sigma
  END TYPE population_colvar_type

! *****************************************************************************
  TYPE gyration_colvar_type
     LOGICAL                         :: use_kinds
     INTEGER                         :: n_atoms
     INTEGER, POINTER, DIMENSION (:) :: i_at
     CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: c_kinds
  END TYPE gyration_colvar_type

! *****************************************************************************
  TYPE torsion_colvar_type
     REAL(KIND=dp)                   :: o0
     INTEGER, DIMENSION (4)          :: i_at_tors
  END TYPE torsion_colvar_type

! *****************************************************************************
  TYPE plane_distance_colvar_type
     LOGICAL               :: use_pbc
     INTEGER, DIMENSION(3) :: plane
     INTEGER               :: point
  END TYPE plane_distance_colvar_type

! *****************************************************************************
  TYPE plane_def_type
     INTEGER                     :: type_of_def
     INTEGER, DIMENSION(3)       :: points
     REAL(KIND=dp), DIMENSION(3) :: normal_vec
  END TYPE plane_def_type

  TYPE plane_plane_angle_colvar_type
     TYPE(plane_def_type)        :: plane1, plane2
  END TYPE plane_plane_angle_colvar_type

! *****************************************************************************
  TYPE angle_colvar_type
     INTEGER, DIMENSION(3) :: i_at_angle
  END TYPE angle_colvar_type

! *****************************************************************************
  TYPE rotation_colvar_type
     INTEGER :: i_at1_bond1,&
                i_at2_bond1,&
                i_at1_bond2,&
                i_at2_bond2
  END TYPE rotation_colvar_type

! *****************************************************************************
  TYPE dfunct_colvar_type
     INTEGER, DIMENSION(4)         :: i_at_dfunct
     LOGICAL                       :: use_pbc
     REAL(KIND=dp)                 :: coeff
  END TYPE dfunct_colvar_type

! *****************************************************************************
  TYPE qparm_colvar_type
     INTEGER                         :: l
     INTEGER                         :: n_atoms_to,&
                                        n_atoms_from
     INTEGER, POINTER, DIMENSION (:) :: i_at_from,&
                                        i_at_to
     REAL(KIND=dp)                   :: rcut, alpha
  END TYPE qparm_colvar_type

! *****************************************************************************
  TYPE hydronium_colvar_type
     INTEGER                         :: n_oxygens,&
                                        n_hydrogens,&
                                        pnh,qnh,pno,qno,&
                                        p,q
     INTEGER, POINTER, DIMENSION (:) :: i_oxygens, &
                                        i_hydrogens
     REAL(KIND=dp)                   :: r_OO, r_OH, lambda,nh
  END TYPE hydronium_colvar_type

! *****************************************************************************
  TYPE reaction_path_colvar_type
     INTEGER                                  :: type_id
     INTEGER                                  :: n_components, nr_frames, subset
     INTEGER,DIMENSION(2)                     :: function_bounds
     INTEGER, POINTER, DIMENSION (:)          :: i_rmsd
     LOGICAL                                  :: align_frames, dist_rmsd, rmsd
     REAL(KIND=dp), DIMENSION (:,:), POINTER  :: f_vals
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: r_ref
     REAL(KIND=dp)                            :: lambda
     REAL(KIND=dp)                            :: step_size
     TYPE(colvar_p_type),POINTER,DIMENSION(:) :: colvar_p
  END TYPE reaction_path_colvar_type

! *****************************************************************************
  TYPE combine_colvar_type
     INTEGER                                  :: type_id
     TYPE(colvar_p_type),POINTER,DIMENSION(:) :: colvar_p
     REAL(KIND=dp)                            :: lerr, dx
     CHARACTER(LEN=default_path_length)       :: FUNCTION
     CHARACTER(LEN=default_string_length),&
          DIMENSION(:), POINTER               :: c_parameters, variables
     REAL(KIND=dp), DIMENSION(:), POINTER     :: v_parameters
  END TYPE combine_colvar_type
! *****************************************************************************
  TYPE rmsd_colvar_type
     INTEGER                                  :: n_atoms, nr_frames, subset
     INTEGER, POINTER, DIMENSION (:)          :: i_rmsd
     LOGICAL                                  :: align_frames
     REAL(KIND=dp), DIMENSION(:), POINTER     :: weights
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: r_ref
  END TYPE rmsd_colvar_type

! *****************************************************************************
  TYPE point_type
     INTEGER :: type_id
     INTEGER, DIMENSION (:), POINTER          :: atoms
     REAL(KIND=dp), DIMENSION (:), POINTER    :: weights
     REAL(KIND=dp), DIMENSION(3)              :: r
  END TYPE point_type

! *****************************************************************************
  TYPE xyz_diag_colvar_type
     LOGICAL                                  :: use_pbc
     LOGICAL                                  :: use_absolute_position
     INTEGER                                  :: i_atom
     INTEGER                                  :: component
     REAL(KIND=dp), DIMENSION(3)              :: r0
  END TYPE xyz_diag_colvar_type

! *****************************************************************************
  TYPE xyz_outerdiag_colvar_type
     LOGICAL                                  :: use_pbc
     INTEGER, DIMENSION(2)                    :: i_atoms
     INTEGER, DIMENSION(2)                    :: components
     REAL(KIND=dp), DIMENSION(3,2)            :: r0
  END TYPE xyz_outerdiag_colvar_type

! *****************************************************************************
  TYPE u_colvar_type
     TYPE(section_vals_type), POINTER         :: mixed_energy_section
     INTEGER       :: natom
  END TYPE u_colvar_type

! *****************************************************************************
  TYPE Wc_colvar_type
    INTEGER                         :: ids(3) ! first is the Od, second the H, third the Oa
    REAL(KIND=dp)                   :: ewc
    REAL(KIND=dp)                   :: rcut
  END TYPE Wc_colvar_type

! *****************************************************************************
  TYPE HBP_colvar_type
    INTEGER                         :: nPoints ! number of the points in the path
    INTEGER, POINTER                :: ids(:,:) ! first is the Od, second the H, 
                                                ! third the Oa and contains a row for each intermediate point in the path
    REAL(KIND=dp),POINTER           :: ewc(:) ! one for each point in the path
    REAL(KIND=dp)                   :: rcut
    REAL(KIND=dp)                   :: shift  ! shift applied for each term in the collective variable
  END TYPE HBP_colvar_type

! *****************************************************************************
  TYPE ring_puckering_colvar_type
     INTEGER                         :: nring
     INTEGER, POINTER, DIMENSION (:) :: atoms
     INTEGER                         :: iq
  END TYPE ring_puckering_colvar_type

! *****************************************************************************
  TYPE mindist_colvar_type
     LOGICAL                    :: use_kinds_from, use_kinds_to
     INTEGER                         :: n_coord_to,&
                                        n_coord_from,&
                                        n_dist_from,&
                                        p_exp,&
                                        q_exp 
     INTEGER, POINTER, DIMENSION (:) :: i_coord_from,&
                                        i_coord_to,&
                                        i_dist_from
     CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: k_coord_from,&
                                                                    k_coord_to
     REAL(KIND=dp)                   :: lambda, r_cut
  END TYPE mindist_colvar_type

! *****************************************************************************
!> \brief parameters for a collective variable
!> \author alessandro laio and fawzi mohamed
! *****************************************************************************
  TYPE colvar_type
     INTEGER                                     :: type_id
     LOGICAL                                     :: use_points
     REAL(kind=dp)                               :: ss               ! Value of the colvar
     REAL(kind=dp), DIMENSION (:,:), POINTER     :: dsdr             ! Derivative of colvar (3,:)
     INTEGER, DIMENSION (:), POINTER             :: i_atom           ! Mapping of dsdr
     INTEGER                                     :: n_atom_s
     ! Available COLVAR types
     TYPE(dist_colvar_type),POINTER              :: dist_param
     TYPE(coord_colvar_type),POINTER             :: coord_param
     TYPE(population_colvar_type),POINTER        :: population_param
     TYPE(gyration_colvar_type),POINTER          :: gyration_param
     TYPE(torsion_colvar_type),POINTER           :: torsion_param
     TYPE(angle_colvar_type),POINTER             :: angle_param
     TYPE(plane_distance_colvar_type),POINTER    :: plane_distance_param
     TYPE(plane_plane_angle_colvar_type),POINTER :: plane_plane_angle_param
     TYPE(rotation_colvar_type), POINTER         :: rotation_param
     TYPE(dfunct_colvar_type), POINTER           :: dfunct_param
     TYPE(qparm_colvar_type), POINTER            :: qparm_param
     TYPE(hydronium_colvar_type), POINTER        :: hydronium_param
     TYPE(reaction_path_colvar_type),POINTER     :: reaction_path_param
     TYPE(combine_colvar_type), POINTER          :: combine_cvs_param
     TYPE(rmsd_colvar_type), POINTER             :: rmsd_param
     TYPE(xyz_diag_colvar_type), POINTER         :: xyz_diag_param
     TYPE(xyz_outerdiag_colvar_type), POINTER    :: xyz_outerdiag_param
     TYPE(u_colvar_type),POINTER                 :: u_param
     TYPE(point_type), DIMENSION(:), POINTER     :: points
     TYPE(Wc_colvar_type), POINTER               :: Wc
     TYPE(HBP_colvar_type), POINTER              :: HBP
     TYPE(ring_puckering_colvar_type), POINTER   :: ring_puckering_param
     TYPE(mindist_colvar_type), POINTER          :: mindist_param
  END TYPE colvar_type

! *****************************************************************************
  TYPE colvar_p_type
     TYPE(colvar_type), POINTER :: colvar
  END TYPE colvar_p_type

! *****************************************************************************
  TYPE colvar_counters
     INTEGER :: ndist
     INTEGER :: nangle
     INTEGER :: ntorsion
     INTEGER :: ncoord
     INTEGER :: nplane_dist
     INTEGER :: nplane_angle
     INTEGER :: nrot
     INTEGER :: ndfunct
     INTEGER :: nqparm
     INTEGER :: nhydronium
     INTEGER :: nreactionpath
     INTEGER :: ncombinecvs
     INTEGER :: nrestraint
     INTEGER :: npopulation
     INTEGER :: ngyration
     INTEGER :: nrmsd
     INTEGER :: nxyz_diag
     INTEGER :: nxyz_outerdiag
     INTEGER :: ntot
     INTEGER :: nWc
     INTEGER :: nHBP
     INTEGER :: nring_puckering
  END TYPE colvar_counters

CONTAINS

! *****************************************************************************
!> \brief initializes a colvar_param type
!> \param colvar the colvat to initialize
!> \param colvar_id ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author alessandro laio and fawzi mohamed
! *****************************************************************************
  SUBROUTINE colvar_create(colvar,colvar_id,error)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER, INTENT(in)                      :: colvar_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPostcondition(.NOT.ASSOCIATED(colvar),cp_warning_level,routineP,error,failure)
    ALLOCATE(colvar, stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    ! Nullify all available COLVAR types
    NULLIFY(colvar%dist_param)
    NULLIFY(colvar%angle_param)
    NULLIFY(colvar%torsion_param)
    NULLIFY(colvar%plane_distance_param)
    NULLIFY(colvar%plane_plane_angle_param)
    NULLIFY(colvar%rotation_param)
    NULLIFY(colvar%dfunct_param)
    NULLIFY(colvar%coord_param)
    NULLIFY(colvar%population_param)
    NULLIFY(colvar%gyration_param)
    NULLIFY(colvar%qparm_param)
    NULLIFY(colvar%xyz_diag_param)
    NULLIFY(colvar%xyz_outerdiag_param)
    NULLIFY(colvar%u_param)
    NULLIFY(colvar%hydronium_param)
    NULLIFY(colvar%reaction_path_param)
    NULLIFY(colvar%combine_cvs_param)
    NULLIFY(colvar%rmsd_param)
    ! More types to nullify
    NULLIFY(colvar%points)
    NULLIFY(colvar%dsdr)
    NULLIFY(colvar%i_atom)
    NULLIFY(colvar%Wc)
    NULLIFY(colvar%ring_puckering_param)
    NULLIFY(colvar%mindist_param)
    IF (.NOT.failure) THEN
       colvar%type_id    = colvar_id
       colvar%use_points = .FALSE.
       SELECT CASE (colvar_id)
       CASE(dist_colvar_id)
          ALLOCATE(colvar%dist_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar%dist_param%axis_id = do_clv_xyz
       CASE(coord_colvar_id)
          ALLOCATE(colvar%coord_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(population_colvar_id)
          ALLOCATE(colvar%population_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(gyration_colvar_id)
          ALLOCATE(colvar%gyration_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(angle_colvar_id)
          ALLOCATE(colvar%angle_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(torsion_colvar_id)
          ALLOCATE(colvar%torsion_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(plane_distance_colvar_id)
          ALLOCATE(colvar%plane_distance_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(plane_plane_angle_colvar_id)
          ALLOCATE(colvar%plane_plane_angle_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(rotation_colvar_id)
          ALLOCATE(colvar%rotation_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(dfunct_colvar_id)
          ALLOCATE(colvar%dfunct_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(qparm_colvar_id)
          ALLOCATE(colvar%qparm_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(xyz_diag_colvar_id)
          ALLOCATE(colvar%xyz_diag_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          ! Initialize r0 with dummy..
          colvar%xyz_diag_param%r0=HUGE(0.0_dp)
       CASE(xyz_outerdiag_colvar_id)
          ALLOCATE(colvar%xyz_outerdiag_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          ! Initialize r0 with dummy..
          colvar%xyz_outerdiag_param%r0=HUGE(0.0_dp)
       CASE(u_colvar_id)
          ALLOCATE(colvar%u_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          NULLIFY(colvar%u_param%mixed_energy_section)
       CASE(hydronium_colvar_id)
          ALLOCATE(colvar%hydronium_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(reaction_path_colvar_id)
          ALLOCATE(colvar%reaction_path_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(distance_from_path_colvar_id)
          ALLOCATE(colvar%reaction_path_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(combine_colvar_id)
          ALLOCATE(colvar%combine_cvs_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(rmsd_colvar_id)
          ALLOCATE(colvar%rmsd_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(Wc_colvar_id)
          ALLOCATE(colvar%Wc, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(HBP_colvar_id)
          ALLOCATE(colvar%HBP, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(ring_puckering_colvar_id)
          ALLOCATE(colvar%ring_puckering_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(mindist_colvar_id)
          ALLOCATE(colvar%mindist_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(no_colvar_id)
          ! Do nothing
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_warning_level,routineP,error,failure)
       END SELECT
    END IF

  END SUBROUTINE colvar_create

! *****************************************************************************
!> \brief Check and setup about the use of geometrical points instead of atoms
!> \param colvar the colvat to initialize
!> \param section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino, [teo] 03.2007
! *****************************************************************************
  SUBROUTINE  colvar_check_points(colvar, section, error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_check_points', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, irep, natoms, npoints, &
                                                nrep, nweights, stat
    INTEGER, DIMENSION(:), POINTER           :: atoms
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r, weights
    TYPE(section_vals_type), POINTER         :: point_sections

    failure=.FALSE.
    NULLIFY(point_sections)
    NULLIFY(atoms)
    NULLIFY(weights)
    CPPostcondition(ASSOCIATED(colvar),cp_warning_level,routineP,error,failure)
    point_sections => section_vals_get_subs_vals(section,"POINT",error=error)
    CALL section_vals_get(point_sections, explicit=explicit, error=error)
    IF (explicit) THEN
       colvar%use_points = .TRUE.
       CALL section_vals_get(point_sections,n_repetition=npoints,error=error)
       ALLOCATE(colvar%points(npoints),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Read points definition
       DO i = 1, npoints
          natoms = 0
          nweights = 0
          NULLIFY(colvar%points(i)%atoms)
          NULLIFY(colvar%points(i)%weights)
          CALL section_vals_val_get(point_sections,"TYPE",i_rep_section=i,i_val=colvar%points(i)%type_id,error=error)
          SELECT CASE(colvar%points(i)%type_id)
          CASE(do_clv_geo_center)
             ! Define a point through a list of atoms..
             CALL section_vals_val_get(point_sections,"ATOMS",i_rep_section=i,n_rep_val=nrep,i_vals=atoms,error=error)
             DO irep=1,nrep
               CALL section_vals_val_get(point_sections,"ATOMS",i_rep_section=i,i_rep_val=irep,i_vals=atoms,error=error)
               natoms = natoms + SIZE(atoms)
             END DO
             ALLOCATE(colvar%points(i)%atoms(natoms),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             natoms = 0
             DO irep=1,nrep
               CALL section_vals_val_get(point_sections,"ATOMS",i_rep_section=i,i_rep_val=irep,i_vals=atoms,error=error)
               colvar%points(i)%atoms(natoms+1:) = atoms(:)
               natoms = natoms + SIZE(atoms)
             END DO
             ! Define weights of the list
             ALLOCATE(colvar%points(i)%weights(natoms),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             colvar%points(i)%weights = 1.0_dp/REAL(natoms,KIND=dp)
             CALL section_vals_val_get(point_sections,"WEIGHTS",i_rep_section=i,n_rep_val=nrep,error=error)
             IF (nrep/=0) THEN
                DO irep=1,nrep
                   CALL section_vals_val_get(point_sections,"WEIGHTS",i_rep_section=i,i_rep_val=irep,&
                        r_vals=weights,error=error)
                   colvar%points(i)%weights(nweights+1:) = weights(:)
                   nweights = nweights + SIZE(weights)
                END DO
                CPPostcondition(natoms==nweights,cp_failure_level,routineP,error,failure)
             END IF
          CASE(do_clv_fix_point )
             ! Define the point as a fixed point in space..
             CALL section_vals_val_get(point_sections,"XYZ",i_rep_section=i,r_vals=r,error=error)
             colvar%points(i)%r = r
          END SELECT
       END DO
    END IF
  END SUBROUTINE colvar_check_points

! *****************************************************************************
!> \brief Finalize the setup of the collective variable
!> \param colvar the colvar to initialize
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino, [teo] 09.03.2006
! *****************************************************************************
  SUBROUTINE colvar_setup(colvar,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, idum, iend, ii, istart, j, &
                                                np, stat
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: failure

    CPPostcondition(ASSOCIATED(colvar),cp_warning_level,routineP,error,failure)
    failure=.FALSE.
    IF (.NOT.failure) THEN
       SELECT CASE(colvar%type_id)
       CASE(dist_colvar_id)
          np = 2
          i = colvar%dist_param%i_at
          j = colvar%dist_param%j_at
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=COLV_SIZE(colvar,i)+&
                          COLV_SIZE(colvar,j)
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          list(1) = colvar%dist_param%i_at
          list(2) = colvar%dist_param%j_at
       CASE(coord_colvar_id)
          np = colvar%coord_param%n_atoms_from+colvar%coord_param%n_atoms_to&
              + colvar%coord_param%n_atoms_to_b
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=0
          DO ii=1,colvar%coord_param%n_atoms_from
             i = colvar%coord_param%i_at_from(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          DO ii=1,colvar%coord_param%n_atoms_to
             i = colvar%coord_param%i_at_to(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          IF(colvar%coord_param%n_atoms_to_b/=0) THEN
            DO ii=1,colvar%coord_param%n_atoms_to_b
               i = colvar%coord_param%i_at_to_b(ii)
               colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
            ENDDO
          END IF
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1,colvar%coord_param%n_atoms_from
             idum = idum + 1
             i = colvar%coord_param%i_at_from(ii)
             list(idum) = i
          ENDDO
          DO ii=1,colvar%coord_param%n_atoms_to
             idum = idum + 1
             i = colvar%coord_param%i_at_to(ii)
             list(idum) = i
          ENDDO
          IF(colvar%coord_param%n_atoms_to_b/=0) THEN
            DO ii=1,colvar%coord_param%n_atoms_to_b
               idum = idum + 1
               i = colvar%coord_param%i_at_to_b(ii)
               list(idum) = i
            ENDDO
          END IF
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(population_colvar_id)
          np = colvar%population_param%n_atoms_from+colvar%population_param%n_atoms_to
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=0
          DO ii=1,colvar%population_param%n_atoms_from
             i = colvar%population_param%i_at_from(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          DO ii=1,colvar%population_param%n_atoms_to
             i = colvar%population_param%i_at_to(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1,colvar%population_param%n_atoms_from
             idum = idum + 1
             i = colvar%population_param%i_at_from(ii)
             list(idum) = i
          ENDDO
          DO ii=1,colvar%population_param%n_atoms_to
             idum = idum + 1
             i = colvar%population_param%i_at_to(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(gyration_colvar_id)
          np = colvar%gyration_param%n_atoms
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=0
          DO ii=1,colvar%gyration_param%n_atoms
             i = colvar%gyration_param%i_at(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
         ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1,colvar%gyration_param%n_atoms
             idum = idum + 1
             i = colvar%gyration_param%i_at(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(angle_colvar_id)
          np = 3
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          DO ii=1,3
             i = colvar%angle_param%i_at_angle(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1, 3
             idum = idum + 1
             i = colvar%angle_param%i_at_angle(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(torsion_colvar_id)
          np = 4
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          DO ii=1,4
             i = colvar%torsion_param%i_at_tors(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1, 4
             idum = idum + 1
             i = colvar%torsion_param%i_at_tors(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(plane_distance_colvar_id)
          np = 4
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          DO ii=1,3
             i = colvar%plane_distance_param%plane(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          i = colvar%plane_distance_param%point
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1, 3
             idum = idum + 1
             i = colvar%plane_distance_param%plane(ii)
             list(idum) = i
          ENDDO
          i = colvar%plane_distance_param%point
          list(4) = i
          idum = idum + 1
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(plane_plane_angle_colvar_id)
          np = 0
          IF (colvar%plane_plane_angle_param%plane1%type_of_def==plane_def_atoms) np = np +3
          IF (colvar%plane_plane_angle_param%plane2%type_of_def==plane_def_atoms) np = np +3
          ! if np is equal to zero this means that this is not a COLLECTIVE variable..
          CALL cp_assert(np/=0,cp_fatal_level,cp_assertion_failed,routineP,&
               "PLANE_PLANE_ANGLE Colvar defined using two normal vectors! This is "//&
               "not a COLLECTIVE VARIABLE! One of the two planes must be defined "//&
               "using atomic positions. "//&
CPSourceFileRef)

          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          IF (colvar%plane_plane_angle_param%plane1%type_of_def==plane_def_atoms) THEN
             DO ii = 1, 3
                i = colvar%plane_plane_angle_param%plane1%points(ii)
                colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
             ENDDO
          END IF
          IF (colvar%plane_plane_angle_param%plane2%type_of_def==plane_def_atoms) THEN
             DO ii = 1, 3
                i = colvar%plane_plane_angle_param%plane2%points(ii)
                colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
             ENDDO
          END IF

          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          IF (colvar%plane_plane_angle_param%plane1%type_of_def==plane_def_atoms) THEN
             DO ii=1, 3
                idum = idum + 1
                i = colvar%plane_plane_angle_param%plane1%points(ii)
                list(idum) = i
             ENDDO
          END IF
          IF (colvar%plane_plane_angle_param%plane2%type_of_def==plane_def_atoms) THEN
             DO ii=1, 3
                idum = idum + 1
                i = colvar%plane_plane_angle_param%plane2%points(ii)
                list(idum) = i
             ENDDO
          END IF
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(dfunct_colvar_id)
          np = 4
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          DO ii=1,4
             i = colvar%dfunct_param%i_at_dfunct(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1, 4
             idum = idum + 1
             i = colvar%dfunct_param%i_at_dfunct(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(rotation_colvar_id)
          np = 4
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          i = colvar%rotation_param%i_at1_bond1
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          i = colvar%rotation_param%i_at2_bond1
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          i = colvar%rotation_param%i_at1_bond2
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          i = colvar%rotation_param%i_at2_bond2
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          i = colvar%rotation_param%i_at1_bond1
          list(1) = i
          i = colvar%rotation_param%i_at2_bond1
          list(2) = i
          i = colvar%rotation_param%i_at1_bond2
          list(3) = i
          i = colvar%rotation_param%i_at2_bond2
          list(4) = i
       CASE(qparm_colvar_id)
          np = colvar%qparm_param%n_atoms_from+colvar%qparm_param%n_atoms_to
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=0
          DO ii=1,colvar%qparm_param%n_atoms_from
             i = colvar%qparm_param%i_at_from(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          DO ii=1,colvar%qparm_param%n_atoms_to
             i = colvar%qparm_param%i_at_to(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1,colvar%qparm_param%n_atoms_from
             idum = idum + 1
             i = colvar%qparm_param%i_at_from(ii)
             list(idum) = i
          ENDDO
          DO ii=1,colvar%qparm_param%n_atoms_to
             idum = idum + 1
             i = colvar%qparm_param%i_at_to(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(hydronium_colvar_id)
          np = colvar%hydronium_param%n_oxygens+colvar%hydronium_param%n_hydrogens
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=0
          DO ii=1,colvar%hydronium_param%n_oxygens
             i = colvar%hydronium_param%i_oxygens(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          DO ii=1,colvar%hydronium_param%n_hydrogens
             i = colvar%hydronium_param%i_hydrogens(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1,colvar%hydronium_param%n_oxygens
             idum = idum + 1
             i = colvar%hydronium_param%i_oxygens(ii)
             list(idum) = i
          ENDDO
          DO ii=1,colvar%hydronium_param%n_hydrogens
             idum = idum + 1
             i = colvar%hydronium_param%i_hydrogens(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(rmsd_colvar_id)
          np = colvar%rmsd_param%n_atoms
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=0
          DO ii=1,colvar%rmsd_param%n_atoms
             i = colvar%rmsd_param%i_rmsd(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1,colvar%rmsd_param%n_atoms
            idum = idum+1
            i = colvar%rmsd_param%i_rmsd(ii)
            list(idum) = i
          END DO
       CASE(reaction_path_colvar_id, distance_from_path_colvar_id)
          colvar%n_atom_s=0
          IF(colvar%reaction_path_param%dist_rmsd .OR. colvar%reaction_path_param%rmsd) THEN
            colvar%n_atom_s = colvar%reaction_path_param%n_components
          ELSE
            DO ii=1,SIZE(colvar%reaction_path_param%colvar_p)
               colvar%n_atom_s=colvar%n_atom_s+colvar%reaction_path_param%colvar_p(ii)%colvar%n_atom_s
            END DO
          END IF
          ALLOCATE(list(colvar%n_atom_s))
          idum=0
          IF(colvar%reaction_path_param%dist_rmsd .OR. colvar%reaction_path_param%rmsd) THEN
            DO ii=1,SIZE(colvar%reaction_path_param%i_rmsd)
               idum = idum + 1
               i = colvar%reaction_path_param%i_rmsd(ii)
               list(idum) = i
            END DO
          ELSE
            DO ii=1,SIZE(colvar%reaction_path_param%colvar_p)
               DO j=1,colvar%reaction_path_param%colvar_p(ii)%colvar%n_atom_s
                  idum=idum+1
                  list(idum)=colvar%reaction_path_param%colvar_p(ii)%colvar%i_atom(j)
               END DO
            END DO
          END IF
       CASE(xyz_diag_colvar_id)
          np = 1
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          i = colvar%xyz_diag_param%i_atom
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          i = colvar%xyz_diag_param%i_atom
          list(1) = i
       CASE(xyz_outerdiag_colvar_id)
          np = 2
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          i = colvar%xyz_outerdiag_param%i_atoms(1)
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          i = colvar%xyz_outerdiag_param%i_atoms(2)
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          i = colvar%xyz_outerdiag_param%i_atoms(1)
          list(1) = i
          i = colvar%xyz_outerdiag_param%i_atoms(2)
          list(2) = i
       CASE(u_colvar_id)
          np = 1; ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          colvar%n_atom_s=np; list(1)=1
       CASE(Wc_colvar_id)
          np = 3
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          DO ii=1,3
             i = colvar%Wc%ids(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1, 3
             idum = idum + 1
             i = colvar%Wc%ids(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(HBP_colvar_id)
          np = 3*colvar%HBP%nPoints
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          DO j=1,colvar%HBP%nPoints
            DO ii=1,3
              i = colvar%HBP%ids(j,ii)
              colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
            ENDDO
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO j=1,colvar%HBP%nPoints
            DO ii=1, 3
              idum = idum + 1
              i = colvar%HBP%ids(j,ii)
              list(idum) = i
            ENDDO
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(ring_puckering_colvar_id)
          np = colvar%ring_puckering_param%nring
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=0
          DO ii=1,colvar%ring_puckering_param%nring
             i = colvar%ring_puckering_param%atoms(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1,colvar%ring_puckering_param%nring
             idum = idum + 1
             i = colvar%ring_puckering_param%atoms(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(mindist_colvar_id)
          np =colvar%mindist_param%n_dist_from +&
              colvar%mindist_param%n_coord_from+colvar%mindist_param%n_coord_to
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=0
          DO ii=1,colvar%mindist_param%n_dist_from
               i = colvar%mindist_param%i_dist_from(ii)
               colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          DO ii=1,colvar%mindist_param%n_coord_from
             i = colvar%mindist_param%i_coord_from(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          DO ii=1,colvar%mindist_param%n_coord_to
             i = colvar%mindist_param%i_coord_to(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1,colvar%mindist_param%n_dist_from
             idum = idum + 1
             i = colvar%mindist_param%i_dist_from(ii)
             list(idum) = i
          ENDDO
          DO ii=1,colvar%mindist_param%n_coord_from
             idum = idum + 1
             i = colvar%mindist_param%i_coord_from(ii)
             list(idum) = i
          ENDDO
          DO ii=1,colvar%mindist_param%n_coord_to
             idum = idum + 1
             i = colvar%mindist_param%i_coord_to(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(combine_colvar_id)
          colvar%n_atom_s=0
          DO ii=1,SIZE(colvar%combine_cvs_param%colvar_p)
             colvar%n_atom_s=colvar%n_atom_s+colvar%combine_cvs_param%colvar_p(ii)%colvar%n_atom_s
          END DO
          ALLOCATE(list(colvar%n_atom_s))
          idum=0
          DO ii=1,SIZE(colvar%combine_cvs_param%colvar_p)
             DO j=1,colvar%combine_cvs_param%colvar_p(ii)%colvar%n_atom_s
                idum=idum+1
                list(idum)=colvar%combine_cvs_param%colvar_p(ii)%colvar%i_atom(j)
             END DO
          END DO
       END SELECT

       IF (ASSOCIATED(colvar%dsdr)) THEN
          DEALLOCATE(colvar%dsdr, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       IF (ASSOCIATED(colvar%i_atom)) THEN
          DEALLOCATE(colvar%i_atom, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       ALLOCATE(colvar%dsdr(3,colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(colvar%i_atom(colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ! And now map real atoms
       istart = 0
       iend   = 0
       DO i = 1, SIZE(list)
          IF (.NOT.colvar%use_points) THEN
             ! No point centers
             colvar%i_atom(i) = list(i)
             iend = iend + 1
          ELSE
             IF (ASSOCIATED(colvar%points(list(i))%atoms)) THEN
                iend = istart + SIZE(colvar%points(list(i))%atoms)
                colvar%i_atom(istart+1:iend) = colvar%points(list(i))%atoms
                istart = iend
             END IF
          END IF
       END DO
       CPPostcondition(iend==colvar%n_atom_s,cp_failure_level,routinep,error,failure)
       DEALLOCATE(list, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

    END IF
  END SUBROUTINE colvar_setup

! *****************************************************************************
!> \brief Gives back the size of an array of integer. If not associated gives back 1
!> \param colvar ...
!> \param i ...
!> \retval my_size ...
!> \author Teodoro Laino - 03.2007
! *****************************************************************************
  FUNCTION colv_size(colvar,i) RESULT(my_size)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER                                  :: i, my_size

    my_size = 1
    IF (ASSOCIATED(colvar%points)) THEN
       IF (ASSOCIATED(colvar%points(i)%atoms)) THEN
          my_size = SIZE(colvar%points(i)%atoms)
       ELSE
          my_size = 0
       END IF
    END IF
  END FUNCTION colv_size

! *****************************************************************************
!> \brief releases the memory that might have been allocated by the colvar
!> \param colvar the colvar to deallocate
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author alessandro laio and fawzi mohamed
! *****************************************************************************
  RECURSIVE SUBROUTINE colvar_release(colvar,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPostcondition(ASSOCIATED(colvar),cp_warning_level,routineP,error,failure)
    IF (ASSOCIATED(colvar%dsdr)) THEN
       DEALLOCATE(colvar%dsdr, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    END IF
    IF (ASSOCIATED(colvar%i_atom)) THEN
        DEALLOCATE(colvar%i_atom, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    END IF
    IF (ASSOCIATED(colvar%points)) THEN
       DO i = 1, SIZE(colvar%points)
          IF (ASSOCIATED(colvar%points(i)%atoms)) THEN
             DEALLOCATE(colvar%points(i)%atoms, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          END IF
          IF (ASSOCIATED(colvar%points(i)%weights)) THEN
             DEALLOCATE(colvar%points(i)%weights, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          END IF
       END DO
       DEALLOCATE(colvar%points, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    END IF
    SELECT CASE(colvar%type_id)
    CASE (dist_colvar_id)
       DEALLOCATE(colvar%dist_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (coord_colvar_id)
       IF (ASSOCIATED(colvar%coord_param%i_at_from)) THEN
          DEALLOCATE(colvar%coord_param%i_at_from, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       IF (ASSOCIATED(colvar%coord_param%i_at_to)) THEN
          DEALLOCATE(colvar%coord_param%i_at_to,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(colvar%coord_param%c_kinds_from)) THEN
          DEALLOCATE(colvar%coord_param%c_kinds_from, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       IF (ASSOCIATED(colvar%coord_param%c_kinds_to)) THEN
          DEALLOCATE(colvar%coord_param%c_kinds_to,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(colvar%coord_param%i_at_to_b)) THEN
          DEALLOCATE(colvar%coord_param%i_at_to_b,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(colvar%coord_param%c_kinds_to_b)) THEN
          DEALLOCATE(colvar%coord_param%c_kinds_to_b,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       DEALLOCATE(colvar%coord_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (population_colvar_id)
       IF (ASSOCIATED(colvar%population_param%i_at_from)) THEN
          DEALLOCATE(colvar%population_param%i_at_from, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       IF (ASSOCIATED(colvar%population_param%i_at_to)) THEN
          DEALLOCATE(colvar%population_param%i_at_to,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(colvar%population_param%c_kinds_from)) THEN
          DEALLOCATE(colvar%population_param%c_kinds_from, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       IF (ASSOCIATED(colvar%population_param%c_kinds_to)) THEN
          DEALLOCATE(colvar%population_param%c_kinds_to,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       DEALLOCATE(colvar%population_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (gyration_colvar_id)
       IF (ASSOCIATED(colvar%gyration_param%i_at)) THEN
          DEALLOCATE(colvar%gyration_param%i_at, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       IF (ASSOCIATED(colvar%gyration_param%c_kinds)) THEN
          DEALLOCATE(colvar%gyration_param%c_kinds, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       DEALLOCATE(colvar%gyration_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE(angle_colvar_id)
       DEALLOCATE(colvar%angle_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (torsion_colvar_id)
       DEALLOCATE(colvar%torsion_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE(plane_distance_colvar_id)
       DEALLOCATE(colvar%plane_distance_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE(plane_plane_angle_colvar_id)
       DEALLOCATE(colvar%plane_plane_angle_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE(dfunct_colvar_id)
       DEALLOCATE(colvar%dfunct_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (rotation_colvar_id)
       DEALLOCATE(colvar%rotation_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (qparm_colvar_id)
       DEALLOCATE(colvar%qparm_param%i_at_from, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%qparm_param%i_at_to,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       DEALLOCATE(colvar%qparm_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (xyz_diag_colvar_id)
       DEALLOCATE(colvar%xyz_diag_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (xyz_outerdiag_colvar_id)
       DEALLOCATE(colvar%xyz_outerdiag_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (u_colvar_id)
       NULLIFY(colvar%u_param%mixed_energy_section)
       DEALLOCATE(colvar%u_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (hydronium_colvar_id)
       DEALLOCATE(colvar%hydronium_param%i_oxygens, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%hydronium_param%i_hydrogens,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       DEALLOCATE(colvar%hydronium_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (reaction_path_colvar_id, distance_from_path_colvar_id)
       IF(colvar%reaction_path_param%dist_rmsd .OR. colvar%reaction_path_param%rmsd) THEN
         DEALLOCATE(colvar%reaction_path_param%r_ref,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
         DEALLOCATE(colvar%reaction_path_param%i_rmsd,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ELSE
         DO i=1,SIZE(colvar%reaction_path_param%colvar_p)
             CALL colvar_release(colvar%reaction_path_param%colvar_p(i)%colvar,error)
         END DO
         DEALLOCATE(colvar%reaction_path_param%colvar_p, stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
         DEALLOCATE(colvar%reaction_path_param%f_vals, stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       DEALLOCATE(colvar%reaction_path_param,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (combine_colvar_id)
       DO i=1,SIZE(colvar%combine_cvs_param%colvar_p)
          CALL colvar_release(colvar%combine_cvs_param%colvar_p(i)%colvar,error)
       END DO
       DEALLOCATE(colvar%combine_cvs_param%colvar_p, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%combine_cvs_param%c_parameters, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%combine_cvs_param%v_parameters, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%combine_cvs_param%variables, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%combine_cvs_param,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (rmsd_colvar_id)
       DEALLOCATE(colvar%rmsd_param%weights, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%rmsd_param%r_ref, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%rmsd_param%i_rmsd, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%rmsd_param,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE(Wc_colvar_id)
       DEALLOCATE(colvar%Wc, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE(HBP_colvar_id)
       DEALLOCATE(colvar%HBP%ewc, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%HBP%ids, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%HBP, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (ring_puckering_colvar_id)
       DEALLOCATE(colvar%ring_puckering_param%atoms, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%ring_puckering_param, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    CASE (mindist_colvar_id)
       IF (ASSOCIATED(colvar%mindist_param%i_dist_from)) THEN
          DEALLOCATE(colvar%mindist_param%i_dist_from, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       IF (ASSOCIATED(colvar%mindist_param%i_coord_from)) THEN
          DEALLOCATE(colvar%mindist_param%i_coord_from, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       IF (ASSOCIATED(colvar%mindist_param%i_coord_to)) THEN
          DEALLOCATE(colvar%mindist_param%i_coord_to,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(colvar%mindist_param%k_coord_from)) THEN
          DEALLOCATE(colvar%mindist_param%k_coord_from, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       END IF
       IF (ASSOCIATED(colvar%mindist_param%k_coord_to)) THEN
          DEALLOCATE(colvar%mindist_param%k_coord_to,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       DEALLOCATE(colvar%mindist_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE(no_colvar_id)
       ! Do nothing
    CASE default
       CPAssert(.FALSE.,cp_failure_level,routinep,error,failure)
    END SELECT
    DEALLOCATE(colvar, stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)

  END SUBROUTINE colvar_release

! *****************************************************************************
!> \brief Clone a colvar type
!> \param colvar_out ...
!> \param colvar_in the colvar to deallocate
!> \param i_atom_offset ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] 04.2006
! *****************************************************************************
  RECURSIVE SUBROUTINE colvar_clone(colvar_out, colvar_in, i_atom_offset, error)
    TYPE(colvar_type), POINTER               :: colvar_out, colvar_in
    INTEGER, INTENT(IN), OPTIONAL            :: i_atom_offset
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_clone', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, my_offset, ndim, ndim2, &
                                                stat
    LOGICAL                                  :: failure

    my_offset = 0
    failure=.FALSE.
    CPPostcondition(ASSOCIATED(colvar_in),cp_warning_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(colvar_out),cp_warning_level,routineP,error,failure)
    IF (PRESENT(i_atom_offset)) my_offset = i_atom_offset
    CALL colvar_create(colvar_out, colvar_in%type_id, error)
    CALL colvar_clone_points(colvar_out, colvar_in, my_offset, error)
    IF (colvar_in%use_points) my_offset = 0
    SELECT CASE(colvar_out%type_id)
    CASE(dist_colvar_id)
       colvar_out%dist_param%i_at = colvar_in%dist_param%i_at + my_offset
       colvar_out%dist_param%j_at = colvar_in%dist_param%j_at + my_offset
       colvar_out%dist_param%axis_id = colvar_in%dist_param%axis_id
    CASE(coord_colvar_id)
       colvar_out%coord_param%n_atoms_to   = colvar_in%coord_param%n_atoms_to
       colvar_out%coord_param%n_atoms_to_b   = colvar_in%coord_param%n_atoms_to_b
       colvar_out%coord_param%n_atoms_from = colvar_in%coord_param%n_atoms_from
       colvar_out%coord_param%nncrd        = colvar_in%coord_param%nncrd
       colvar_out%coord_param%ndcrd        = colvar_in%coord_param%ndcrd
       colvar_out%coord_param%r_0          = colvar_in%coord_param%r_0
       colvar_out%coord_param%nncrd_b        = colvar_in%coord_param%nncrd_b
       colvar_out%coord_param%ndcrd_b        = colvar_in%coord_param%ndcrd_b
       colvar_out%coord_param%r_0_b          = colvar_in%coord_param%r_0_b
       colvar_out%coord_param%use_kinds_from  = colvar_in%coord_param%use_kinds_from
       colvar_out%coord_param%use_kinds_to    = colvar_in%coord_param%use_kinds_to
       colvar_out%coord_param%use_kinds_to_b    = colvar_in%coord_param%use_kinds_to_b
       NULLIFY(colvar_out%coord_param%c_kinds_from, colvar_out%coord_param%i_at_from)
       NULLIFY(colvar_out%coord_param%c_kinds_to, colvar_out%coord_param%i_at_to)
       NULLIFY(colvar_out%coord_param%c_kinds_to_b, colvar_out%coord_param%i_at_to_b)
       IF(colvar_in%coord_param%use_kinds_from) THEN
          ! KINDS
          ndim = SIZE(colvar_in%coord_param%c_kinds_from)
          ALLOCATE(colvar_out%coord_param%c_kinds_from(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%coord_param%c_kinds_from    = colvar_in%coord_param%c_kinds_from
       ELSE
          ! INDEX
          ndim = SIZE(colvar_in%coord_param%i_at_from)
          ALLOCATE(colvar_out%coord_param%i_at_from(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%coord_param%i_at_from    = colvar_in%coord_param%i_at_from + my_offset
       END IF
       IF(colvar_in%coord_param%use_kinds_to) THEN
          ! KINDS
          ndim = SIZE(colvar_in%coord_param%c_kinds_to)
          ALLOCATE(colvar_out%coord_param%c_kinds_to(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%coord_param%c_kinds_to    = colvar_in%coord_param%c_kinds_to
       ELSE
          ! INDEX
          ndim = SIZE(colvar_in%coord_param%i_at_to)
          ALLOCATE(colvar_out%coord_param%i_at_to(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%coord_param%i_at_to      = colvar_in%coord_param%i_at_to + my_offset
       END IF
       IF(colvar_in%coord_param%use_kinds_to_b) THEN
          ! KINDS
          ndim = SIZE(colvar_in%coord_param%c_kinds_to_b)
          ALLOCATE(colvar_out%coord_param%c_kinds_to_b(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%coord_param%c_kinds_to_b    = colvar_in%coord_param%c_kinds_to_b
       ELSEIF(ASSOCIATED(colvar_in%coord_param%i_at_to_b)) THEN
          ! INDEX
          ndim = SIZE(colvar_in%coord_param%i_at_to_b)
          ALLOCATE(colvar_out%coord_param%i_at_to_b(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%coord_param%i_at_to_b      = colvar_in%coord_param%i_at_to_b + my_offset
       END IF


     CASE(population_colvar_id)
       colvar_out%population_param%n_atoms_to   = colvar_in%population_param%n_atoms_to
       colvar_out%population_param%n_atoms_from = colvar_in%population_param%n_atoms_from
       colvar_out%population_param%nncrd        = colvar_in%population_param%nncrd
       colvar_out%population_param%ndcrd        = colvar_in%population_param%ndcrd
       colvar_out%population_param%r_0          = colvar_in%population_param%r_0
       colvar_out%population_param%use_kinds_from  = colvar_in%population_param%use_kinds_from
       colvar_out%population_param%use_kinds_to    = colvar_in%population_param%use_kinds_to
       NULLIFY(colvar_out%population_param%c_kinds_from, colvar_out%population_param%i_at_from)
       NULLIFY(colvar_out%population_param%c_kinds_to, colvar_out%population_param%i_at_to)
       IF(colvar_in%population_param%use_kinds_from) THEN
          ! KINDS
          ndim = SIZE(colvar_in%population_param%c_kinds_from)
          ALLOCATE(colvar_out%population_param%c_kinds_from(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%population_param%c_kinds_from    = colvar_in%population_param%c_kinds_from
       ELSE
          ! INDEX
          ndim = SIZE(colvar_in%population_param%i_at_from)
          ALLOCATE(colvar_out%population_param%i_at_from(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%population_param%i_at_from    = colvar_in%population_param%i_at_from + my_offset
       END IF
       IF(colvar_in%population_param%use_kinds_to) THEN
          ! KINDS
          ndim = SIZE(colvar_in%population_param%c_kinds_to)
          ALLOCATE(colvar_out%population_param%c_kinds_to(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%population_param%c_kinds_to    = colvar_in%population_param%c_kinds_to
       ELSE
          ! INDEX
          ndim = SIZE(colvar_in%population_param%i_at_to)
          ALLOCATE(colvar_out%population_param%i_at_to(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%population_param%i_at_to      = colvar_in%population_param%i_at_to + my_offset
       END IF

     CASE(gyration_colvar_id)
       colvar_out%gyration_param%n_atoms   = colvar_in%gyration_param%n_atoms
       colvar_out%gyration_param%use_kinds    = colvar_in%gyration_param%use_kinds
       NULLIFY(colvar_out%gyration_param%c_kinds, colvar_out%gyration_param%i_at)
       IF(colvar_in%gyration_param%use_kinds) THEN
          ! KINDS
          ndim = SIZE(colvar_in%gyration_param%c_kinds)
          ALLOCATE(colvar_out%gyration_param%c_kinds(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%gyration_param%c_kinds   = colvar_in%gyration_param%c_kinds
       ELSE
          ! INDEX
          ndim = SIZE(colvar_in%gyration_param%i_at)
          ALLOCATE(colvar_out%gyration_param%i_at(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%gyration_param%i_at    = colvar_in%gyration_param%i_at + my_offset
       END IF
    CASE(angle_colvar_id)
       colvar_out%angle_param%i_at_angle   = colvar_in%angle_param%i_at_angle + my_offset
    CASE(torsion_colvar_id)
       colvar_out%torsion_param%i_at_tors     = colvar_in%torsion_param%i_at_tors + my_offset
       colvar_out%torsion_param%o0            = colvar_in%torsion_param%o0
    CASE(plane_distance_colvar_id)
       colvar_out%plane_distance_param%use_pbc = colvar_in%plane_distance_param%use_pbc
       colvar_out%plane_distance_param%plane = colvar_in%plane_distance_param%plane + my_offset
       colvar_out%plane_distance_param%point = colvar_in%plane_distance_param%point + my_offset
    CASE(plane_plane_angle_colvar_id)
       colvar_out%plane_plane_angle_param%plane1%type_of_def = colvar_in%plane_plane_angle_param%plane1%type_of_def
       IF (colvar_out%plane_plane_angle_param%plane1%type_of_def == plane_def_vec) THEN
          colvar_out%plane_plane_angle_param%plane1%normal_vec  = colvar_in%plane_plane_angle_param%plane1%normal_vec
       ELSE
          colvar_out%plane_plane_angle_param%plane1%points = colvar_in%plane_plane_angle_param%plane1%points + my_offset
       END IF

       colvar_out%plane_plane_angle_param%plane2%type_of_def = colvar_in%plane_plane_angle_param%plane2%type_of_def
       IF (colvar_out%plane_plane_angle_param%plane2%type_of_def == plane_def_vec) THEN
          colvar_out%plane_plane_angle_param%plane2%normal_vec  = colvar_in%plane_plane_angle_param%plane2%normal_vec
       ELSE
          colvar_out%plane_plane_angle_param%plane2%points = colvar_in%plane_plane_angle_param%plane2%points + my_offset
       END IF
    CASE(rotation_colvar_id)
       colvar_out%rotation_param%i_at1_bond1  = colvar_in%rotation_param%i_at1_bond1 + my_offset
       colvar_out%rotation_param%i_at2_bond1  = colvar_in%rotation_param%i_at2_bond1 + my_offset
       colvar_out%rotation_param%i_at1_bond2  = colvar_in%rotation_param%i_at1_bond2 + my_offset
       colvar_out%rotation_param%i_at2_bond2  = colvar_in%rotation_param%i_at2_bond2 + my_offset
    CASE(dfunct_colvar_id)
       colvar_out%dfunct_param%i_at_dfunct   = colvar_in%dfunct_param%i_at_dfunct + my_offset
       colvar_out%dfunct_param%coeff         = colvar_in%dfunct_param%coeff
       colvar_out%dfunct_param%use_pbc       = colvar_in%dfunct_param%use_pbc
    CASE(qparm_colvar_id)
       colvar_out%qparm_param%n_atoms_to   = colvar_in%qparm_param%n_atoms_to
       colvar_out%qparm_param%n_atoms_from = colvar_in%qparm_param%n_atoms_from
       colvar_out%qparm_param%rcut         = colvar_in%qparm_param%rcut
       colvar_out%qparm_param%l            = colvar_in%qparm_param%l
       colvar_out%qparm_param%alpha        = colvar_in%qparm_param%alpha
       ndim = SIZE(colvar_in%qparm_param%i_at_from)
       ALLOCATE(colvar_out%qparm_param%i_at_from(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ndim = SIZE(colvar_in%qparm_param%i_at_to)
       ALLOCATE(colvar_out%qparm_param%i_at_to(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%qparm_param%i_at_from    = colvar_in%qparm_param%i_at_from + my_offset
       colvar_out%qparm_param%i_at_to      = colvar_in%qparm_param%i_at_from + my_offset
     CASE(xyz_diag_colvar_id)
       colvar_out%xyz_diag_param%i_atom       = colvar_in%xyz_diag_param%i_atom+ my_offset
       colvar_out%xyz_diag_param%component    = colvar_in%xyz_diag_param%component
       colvar_out%xyz_diag_param%r0           = colvar_in%xyz_diag_param%r0
       colvar_out%xyz_diag_param%use_pbc      = colvar_in%xyz_diag_param%use_pbc
       colvar_out%xyz_diag_param%use_absolute_position      = colvar_in%xyz_diag_param%use_absolute_position
     CASE(xyz_outerdiag_colvar_id)
       colvar_out%xyz_outerdiag_param%i_atoms    = colvar_in%xyz_outerdiag_param%i_atoms+ my_offset
       colvar_out%xyz_outerdiag_param%components = colvar_in%xyz_outerdiag_param%components
       colvar_out%xyz_outerdiag_param%r0         = colvar_in%xyz_outerdiag_param%r0
       colvar_out%xyz_outerdiag_param%use_pbc    = colvar_in%xyz_outerdiag_param%use_pbc
     CASE(u_colvar_id)
       colvar_out%u_param%natom = colvar_in%u_param%natom
     CASE(hydronium_colvar_id)
       colvar_out%hydronium_param%n_hydrogens = colvar_in%hydronium_param%n_hydrogens
       colvar_out%hydronium_param%n_oxygens   = colvar_in%hydronium_param%n_oxygens
       colvar_out%hydronium_param%nh          = colvar_in%hydronium_param%nh
       colvar_out%hydronium_param%pnh         = colvar_in%hydronium_param%pnh
       colvar_out%hydronium_param%pno         = colvar_in%hydronium_param%pno
       colvar_out%hydronium_param%qnh         = colvar_in%hydronium_param%qnh
       colvar_out%hydronium_param%qno         = colvar_in%hydronium_param%qno
       colvar_out%hydronium_param%p           = colvar_in%hydronium_param%p
       colvar_out%hydronium_param%q           = colvar_in%hydronium_param%q
       colvar_out%hydronium_param%r_OO        = colvar_in%hydronium_param%r_OO
       colvar_out%hydronium_param%r_OH        = colvar_in%hydronium_param%r_OH
       ndim = SIZE(colvar_in%hydronium_param%i_oxygens)
       ALLOCATE(colvar_out%hydronium_param%i_oxygens(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ndim = SIZE(colvar_in%hydronium_param%i_hydrogens)
       ALLOCATE(colvar_out%hydronium_param%i_hydrogens(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%hydronium_param%i_oxygens    = colvar_in%hydronium_param%i_oxygens + my_offset
       colvar_out%hydronium_param%i_hydrogens      = colvar_in%hydronium_param%i_oxygens + my_offset
    CASE(reaction_path_colvar_id, distance_from_path_colvar_id)
       colvar_out%reaction_path_param%dist_rmsd = colvar_in%reaction_path_param%dist_rmsd
       colvar_out%reaction_path_param%rmsd = colvar_in%reaction_path_param%rmsd
       colvar_out%reaction_path_param%nr_frames = colvar_in%reaction_path_param%nr_frames
       IF(colvar_in%reaction_path_param%dist_rmsd .OR. colvar_in%reaction_path_param%rmsd) THEN
         colvar_out%reaction_path_param%align_frames = colvar_in%reaction_path_param%align_frames
         colvar_out%reaction_path_param%subset = colvar_in%reaction_path_param%subset
         NULLIFY(colvar_out%reaction_path_param%i_rmsd)
         NULLIFY(colvar_out%reaction_path_param%r_ref)
         ndim = SIZE(colvar_in%reaction_path_param%i_rmsd)
         ALLOCATE(colvar_out%reaction_path_param%i_rmsd(ndim),stat=stat)
         colvar_out%reaction_path_param%i_rmsd = colvar_in%reaction_path_param%i_rmsd
         ndim = SIZE(colvar_in%reaction_path_param%r_ref,1)
         ndim2 = SIZE(colvar_in%reaction_path_param%r_ref,2)
         ALLOCATE(colvar_out%reaction_path_param%r_ref(ndim,ndim2),stat=stat)
         colvar_out%reaction_path_param%r_ref = colvar_in%reaction_path_param%r_ref
       ELSE
         ndim=SIZE(colvar_in%reaction_path_param%colvar_p)
         ALLOCATE(colvar_out%reaction_path_param%colvar_p(ndim),stat=stat)
         CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
         DO i=1,ndim
            NULLIFY(colvar_out%reaction_path_param%colvar_p(i)%colvar)
            CALL colvar_clone(colvar_out%reaction_path_param%colvar_p(i)%colvar,&
                              colvar_in%reaction_path_param%colvar_p(i)%colvar,&
                              my_offset,error)
         END DO
         colvar_out%reaction_path_param% function_bounds = colvar_in%reaction_path_param%  function_bounds
         NULLIFY(colvar_out%reaction_path_param%f_vals)
         ndim=SIZE(colvar_in%reaction_path_param%f_vals,1)
         ndim2=SIZE(colvar_in%reaction_path_param%f_vals,2)
         ALLOCATE(colvar_out%reaction_path_param%f_vals(ndim,ndim2),stat=stat)
         CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
         colvar_out%reaction_path_param%f_vals           = colvar_in%reaction_path_param%f_vals
       END IF
       colvar_out%reaction_path_param%step_size        = colvar_in%reaction_path_param%step_size
       colvar_out%reaction_path_param%n_components     = colvar_in%reaction_path_param%n_components
       colvar_out%reaction_path_param%lambda           = colvar_in%reaction_path_param%lambda
    CASE(combine_colvar_id)
       ndim=SIZE(colvar_in%combine_cvs_param%colvar_p)
       ALLOCATE(colvar_out%combine_cvs_param%colvar_p(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       DO i=1,ndim
          NULLIFY(colvar_out%combine_cvs_param%colvar_p(i)%colvar)
          CALL colvar_clone(colvar_out%combine_cvs_param%colvar_p(i)%colvar,&
                            colvar_in%combine_cvs_param%colvar_p(i)%colvar,&
                            my_offset,error)
       END DO
       colvar_out%combine_cvs_param%lerr         = colvar_in%combine_cvs_param%lerr
       colvar_out%combine_cvs_param%dx           = colvar_in%combine_cvs_param%dx
       colvar_out%combine_cvs_param%function     = colvar_in%combine_cvs_param%function
       !
       ndim = SIZE(colvar_in%combine_cvs_param%c_parameters)
       ALLOCATE(colvar_out%combine_cvs_param%c_parameters(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%combine_cvs_param%c_parameters = colvar_in%combine_cvs_param%c_parameters
       !
       ndim = SIZE(colvar_in%combine_cvs_param%v_parameters)
       ALLOCATE(colvar_out%combine_cvs_param%v_parameters(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%combine_cvs_param%v_parameters = colvar_in%combine_cvs_param%v_parameters
       !
       ndim = SIZE(colvar_in%combine_cvs_param%variables)
       ALLOCATE(colvar_out%combine_cvs_param%variables(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%combine_cvs_param%variables    = colvar_in%combine_cvs_param%variables
     CASE(rmsd_colvar_id)
       colvar_out%rmsd_param%n_atoms   = colvar_in%rmsd_param%n_atoms
       colvar_out%rmsd_param%align_frames  = colvar_in%rmsd_param%align_frames
       colvar_out%rmsd_param%nr_frames   = colvar_in%rmsd_param%nr_frames
       colvar_out%rmsd_param%subset   = colvar_in%rmsd_param%subset
       NULLIFY( colvar_out%rmsd_param%i_rmsd)
       NULLIFY( colvar_out%rmsd_param%weights)
       NULLIFY( colvar_out%rmsd_param%r_ref)
       ! INDEX
       ndim = SIZE(colvar_in%rmsd_param%i_rmsd)
       ALLOCATE(colvar_out%rmsd_param%i_rmsd(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%rmsd_param%i_rmsd    = colvar_in%rmsd_param%i_rmsd + my_offset
       ! A and Bconfigurations and weights
       ndim = SIZE(colvar_in%rmsd_param%weights)
       ALLOCATE(colvar_out%rmsd_param%weights(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%rmsd_param%weights    = colvar_in%rmsd_param%weights
       ndim = SIZE(colvar_in%rmsd_param%r_ref,1)
       ndim2 = SIZE(colvar_in%rmsd_param%r_ref,2)
       ALLOCATE(colvar_out%rmsd_param%r_ref(ndim,ndim2),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%rmsd_param%r_ref    = colvar_in%rmsd_param%r_ref
    CASE(Wc_colvar_id)
       colvar_out%Wc%ids   = colvar_in%Wc%ids + my_offset
       colvar_out%Wc%rcut   = colvar_in%Wc%rcut
    CASE(HBP_colvar_id)
       ndim=colvar_out%HBP%nPoints
       ALLOCATE(colvar_out%HBP%ids(ndim,3),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ALLOCATE(colvar_out%HBP%ewc(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%HBP%ids   = colvar_in%HBP%ids + my_offset
       colvar_out%HBP%ewc   = colvar_in%HBP%ewc + my_offset
       colvar_out%HBP%nPoints   = colvar_in%HBP%nPoints
       colvar_out%HBP%rcut   = colvar_in%HBP%rcut
       colvar_out%HBP%shift   = colvar_in%HBP%shift
     CASE(ring_puckering_colvar_id)
       ndim = colvar_in%ring_puckering_param%nring
       colvar_out%ring_puckering_param%nring = colvar_in%ring_puckering_param%nring
       colvar_out%ring_puckering_param%iq    = colvar_in%ring_puckering_param%iq
       ALLOCATE(colvar_out%ring_puckering_param%atoms(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%ring_puckering_param%atoms = colvar_in%ring_puckering_param%atoms + my_offset
    CASE(mindist_colvar_id)
       colvar_out%mindist_param%n_dist_from   = colvar_in%mindist_param%n_dist_from
       colvar_out%mindist_param%n_coord_to   = colvar_in%mindist_param%n_coord_to
       colvar_out%mindist_param%n_coord_from = colvar_in%mindist_param%n_coord_from
       colvar_out%mindist_param%p_exp        = colvar_in%mindist_param%p_exp
       colvar_out%mindist_param%q_exp        = colvar_in%mindist_param%q_exp
       colvar_out%mindist_param%r_cut          = colvar_in%mindist_param%r_cut
       colvar_out%mindist_param%lambda          = colvar_in%mindist_param%lambda
       colvar_out%mindist_param%use_kinds_from  = colvar_in%mindist_param%use_kinds_from
       colvar_out%mindist_param%use_kinds_to    = colvar_in%mindist_param%use_kinds_to
       NULLIFY(colvar_out%mindist_param%k_coord_from, colvar_out%mindist_param%i_coord_from)
       NULLIFY(colvar_out%mindist_param%k_coord_to, colvar_out%mindist_param%i_coord_to)
       NULLIFY(colvar_out%mindist_param%i_dist_from)
       ! INDEX
       ndim = SIZE(colvar_in%mindist_param%i_dist_from)
       ALLOCATE(colvar_out%mindist_param%i_dist_from(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%mindist_param%i_dist_from    = colvar_in%mindist_param%i_dist_from + my_offset
       IF(colvar_in%mindist_param%use_kinds_from) THEN
          ! KINDS
          ndim = SIZE(colvar_in%mindist_param%k_coord_from)
          ALLOCATE(colvar_out%mindist_param%k_coord_from(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%mindist_param%k_coord_from    = colvar_in%mindist_param%k_coord_from
       ELSE
          ! INDEX
          ndim = SIZE(colvar_in%mindist_param%i_coord_from)
          ALLOCATE(colvar_out%mindist_param%i_coord_from(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%mindist_param%i_coord_from    = colvar_in%mindist_param%i_coord_from + my_offset
       END IF
       IF(colvar_in%mindist_param%use_kinds_to) THEN
          ! KINDS
          ndim = SIZE(colvar_in%mindist_param%k_coord_to)
          ALLOCATE(colvar_out%mindist_param%k_coord_to(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%mindist_param%k_coord_to    = colvar_in%mindist_param%k_coord_to
       ELSE
          ! INDEX
          ndim = SIZE(colvar_in%mindist_param%i_coord_to)
          ALLOCATE(colvar_out%mindist_param%i_coord_to(ndim),stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          colvar_out%mindist_param%i_coord_to      = colvar_in%mindist_param%i_coord_to + my_offset
       END IF

    END SELECT
    CALL colvar_setup(colvar_out, error)
  END SUBROUTINE colvar_clone

! *****************************************************************************
!> \brief Clone points type of a colvar type
!> \param colvar_out ...
!> \param colvar_in the colvar to deallocate
!> \param offset ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] 03.2007
! *****************************************************************************
  SUBROUTINE colvar_clone_points(colvar_out, colvar_in, offset, error)
    TYPE(colvar_type), POINTER               :: colvar_out, colvar_in
    INTEGER, INTENT(IN)                      :: offset
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_clone_points', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, natoms, npoints, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    colvar_out%use_points = colvar_in%use_points
    IF (colvar_in%use_points) THEN
       CPPostcondition(ASSOCIATED(colvar_in%points),cp_failure_level,routineP,error,failure)
       npoints = SIZE(colvar_in%points)
       ALLOCATE(colvar_out%points(npoints),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, npoints
          IF (ASSOCIATED(colvar_in%points(i)%atoms)) THEN
             natoms = SIZE(colvar_in%points(i)%atoms)
             ALLOCATE(colvar_out%points(i)%atoms(natoms),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             colvar_out%points(i)%atoms = colvar_in%points(i)%atoms+offset
          ELSE
             NULLIFY(colvar_out%points(i)%atoms)
          END IF
          IF (ASSOCIATED(colvar_in%points(i)%weights)) THEN
             natoms = SIZE(colvar_in%points(i)%weights)
             ALLOCATE(colvar_out%points(i)%weights(natoms),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             colvar_out%points(i)%weights = colvar_in%points(i)%weights
          ELSE
             NULLIFY(colvar_out%points(i)%weights)
          END IF
          colvar_out%points(i)%type_id = colvar_in%points(i)%type_id
          colvar_out%points(i)%r       = colvar_in%points(i)%r
       END DO
    ELSE
       NULLIFY(colvar_out%points)
    END IF

  END SUBROUTINE colvar_clone_points

! *****************************************************************************
!> \brief Change the dimension of a colvar_p_type
!> \param colvar_set ...
!> \param lb1_new ...
!> \param ub1_new ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] 04.2006
! *****************************************************************************
  SUBROUTINE colvar_p_reallocate(colvar_set,lb1_new,ub1_new, error)
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: colvar_set
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_p_reallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, lb1, lb1_old, stat, ub1, &
                                                ub1_old
    LOGICAL                                  :: failure
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: work

    failure = .FALSE.
    NULLIFY(work)
    IF (ASSOCIATED(colvar_set)) THEN
       lb1_old = LBOUND(colvar_set,1)
       ub1_old = UBOUND(colvar_set,1)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       ALLOCATE (work(lb1:ub1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO j = lb1, ub1
          NULLIFY(work(j)%colvar)
          CALL colvar_clone(work(j)%colvar, colvar_set(j)%colvar, error=error)
       END DO
       DO j = lb1, ub1
          CALL colvar_release(colvar_set(j)%colvar, error)
          NULLIFY(colvar_set(j)%colvar)
       END DO
       DEALLOCATE(colvar_set, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    ALLOCATE (colvar_set(lb1_new:ub1_new),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO j = lb1_new, ub1_new
       NULLIFY(colvar_set(j)%colvar)
    END DO
    IF (ASSOCIATED(work)) THEN
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       DO j = lb1, ub1
          CALL colvar_clone(colvar_set(j)%colvar, work(j)%colvar, error=error)
       END DO
       DO j = lb1, ub1
          CALL colvar_release(work(j)%colvar, error)
          NULLIFY(work(j)%colvar)
       END DO
       DEALLOCATE (work,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE colvar_p_reallocate


! *****************************************************************************
!> \brief Deallocate a set of colvar_p_type
!> \param colvar_p ...
!> \param error ...
!> \par History
!>      07.2003 created [fawzi]
!>      01.2014 moved from cp_subsys_release() into separate routine.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE colvar_p_release(colvar_p, error)
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: colvar_p
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i

! Colvar info

     IF(ASSOCIATED(colvar_p)) THEN
        DO i=1,SIZE(colvar_p)
           IF (ASSOCIATED(colvar_p(i)%colvar)) &
              CALL colvar_release(colvar_p(i)%colvar,error)
        ENDDO
        DEALLOCATE(colvar_p)
     END IF
  END SUBROUTINE colvar_p_release

! *****************************************************************************
!> \brief Evaluate the position of the geometrical point
!> \param point ...
!> \param particles ...
!> \param r ...
!> \author Teodoro Laino - 03.2007
! *****************************************************************************
  SUBROUTINE eval_point_pos(point, particles, r)
    TYPE(point_type)                         :: point
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: r

    INTEGER                                  :: i

    SELECT CASE(point%type_id)
    CASE(do_clv_geo_center)
       r = 0.0_dp
       DO i = 1,SIZE(point%atoms)
          r = r + particles(point%atoms(i))%r*point%weights(i)
       END DO
    CASE(do_clv_fix_point )
       r = point%r
    END SELECT

  END SUBROUTINE eval_point_pos

! *****************************************************************************
!> \brief ...
!> \param point ...
!> \param particles ...
!> \param m ...
! *****************************************************************************
  SUBROUTINE eval_point_mass(point, particles, m)
    TYPE(point_type)                         :: point
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), INTENT(OUT)               :: m

    INTEGER                                  :: i

    SELECT CASE(point%type_id)
    CASE(do_clv_geo_center)
       m = 0.0_dp
       DO i = 1,SIZE(point%atoms)
          m = m + particles(point%atoms(i))%atomic_kind%mass*point%weights(i)
       END DO
    CASE(do_clv_fix_point )
       m = 0.0_dp
    END SELECT

  END SUBROUTINE eval_point_mass

! *****************************************************************************
!> \brief Evaluate the position of the geometrical point
!> \param points ...
!> \param i ...
!> \param dsdr ...
!> \param f ...
!> \author Teodoro Laino - 03.2007
! *****************************************************************************
  SUBROUTINE eval_point_der(points, i, dsdr, f)
    TYPE(point_type), DIMENSION(:), POINTER  :: points
    INTEGER, INTENT(IN)                      :: i
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dsdr
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: f

    INTEGER                                  :: ind, j
    REAL(KIND=dp)                            :: fac

    SELECT CASE(points(i)%type_id)
    CASE(do_clv_geo_center)
       ind = 0
       DO j = 1, i-1
          IF (ASSOCIATED(points(j)%atoms)) THEN
             ind = ind + SIZE(points(j)%atoms)
          END IF
       END DO
       DO j = 1,SIZE(points(i)%atoms)
          fac = points(i)%weights(j)
          dsdr(:,ind+j) = dsdr(:,ind+j) + f * fac
       END DO
    CASE(do_clv_fix_point )
       ! Do nothing if it's a fixed point in space
    END SELECT

  END SUBROUTINE eval_point_der

! *****************************************************************************
!> \brief subtract b from the ss value of a colvar: general function for handling
!>        periodic/non-periodic colvar
!> \param colvar ...
!> \param b ...
!> \retval diff ...
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  FUNCTION diff_colvar(colvar, b) RESULT(diff)
    TYPE(colvar_type), POINTER               :: colvar
    REAL(KIND=dp), INTENT(IN)                :: b
    REAL(KIND=dp)                            :: diff

    CHARACTER(len=*), PARAMETER :: routineN = 'diff_colvar', &
      routineP = moduleN//':'//routineN

    diff = colvar%ss - b
    IF (colvar%type_id==torsion_colvar_id) THEN
       ! The difference of a periodic COLVAR is always within [-pi,pi]
       diff = SIGN(1.0_dp,ASIN(SIN(diff)))*ACOS(COS(diff))
    END IF
  END FUNCTION diff_colvar

END MODULE colvar_types
