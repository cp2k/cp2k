!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module provides variables for the TMC analysis tool
!> \par History
!>      02.2013 created [Mandes Schönherr]
!> \author Mandes 
! *****************************************************************************

MODULE tmc_analysis_types
  USE cell_types,                      ONLY: cell_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE tmc_tree_types,                  ONLY: tree_type
  USE tmc_types,                       ONLY: tmc_atom_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tmc_analysis_types'

  PUBLIC :: tmc_analysis_env, tmc_ana_list_type
  PUBLIC :: tmc_ana_env_create, tmc_ana_env_release
  PUBLIC :: tmc_ana_density_create
  PUBLIC :: pair_correl_type, tmc_ana_pair_correl_create, &
            search_pair_in_list, atom_pairs_type
  PUBLIC :: dipole_moment_type, tmc_ana_dipole_moment_create
  PUBLIC :: dipole_analysis_type, tmc_ana_dipole_analysis_create
  PUBLIC :: displacement_type, tmc_ana_displacement_create


  CHARACTER(LEN= default_path_length ), PARAMETER, &
    PUBLIC :: tmc_ana_density_file_name            = "tmc_ana_density.dat"
  CHARACTER(LEN= default_path_length ), PARAMETER, &
    PUBLIC :: tmc_ana_pair_correl_file_name        = "tmc_ana_g_r.dat"

  INTEGER, PARAMETER, PUBLIC                      :: ana_type_default = 0
  INTEGER, PARAMETER, PUBLIC                      :: ana_type_ice     = 1
  INTEGER, PARAMETER, PUBLIC                      :: ana_type_sym_xyz = 2


  TYPE tmc_ana_list_type
    TYPE(tmc_analysis_env), POINTER               :: temp => NULL()
  END TYPE tmc_ana_list_type

  TYPE tmc_analysis_env
    INTEGER                                       :: io_unit
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                       :: dirs
    CHARACTER(LEN= default_path_length )          :: out_file_prefix
    INTEGER                                       :: conf_offset
    TYPE(cell_type), POINTER                      :: cell
    TYPE(tmc_atom_type), DIMENSION(:), POINTER    :: atoms
    INTEGER                                       :: dim_per_elem = 3
    INTEGER                                       :: nr_dim = -1
    REAL(KIND=dp)                                 :: temperature
    TYPE(tree_type), POINTER                      :: last_elem
    INTEGER                                       :: from_elem, to_elem
    INTEGER                                       :: id_traj, id_cell, id_frc, id_dip, id_ener
    INTEGER                                       :: lc_traj, lc_cell, lc_frc, lc_dip, lc_ener
    CHARACTER(LEN= default_path_length )          :: costum_pos_file_name 
    CHARACTER(LEN= default_path_length )          :: costum_dip_file_name
    CHARACTER(LEN= default_path_length )          :: costum_cell_file_name
    LOGICAL                                       :: restart, restarted
    LOGICAL                                       :: print_test_output

    TYPE(density_3d_type), POINTER                :: density_3d
    TYPE(pair_correl_type), POINTER               :: pair_correl
    TYPE(dipole_moment_type), POINTER             :: dip_mom
    TYPE(dipole_analysis_type), POINTER           :: dip_ana
    TYPE(displacement_type), POINTER              :: displace
  END TYPE tmc_analysis_env

  TYPE density_3d_type
    INTEGER                                       :: conf_counter
    INTEGER, DIMENSION(3)                         :: nr_bins
    REAL(KIND=dp)                                 :: sum_vol
    REAL(KIND=dp)                                 :: sum_vol2
    REAL(KIND=dp), DIMENSION(3)                   :: sum_box_length
    REAL(KIND=dp), DIMENSION(3)                   :: sum_box_length2
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER      :: sum_density, sum_dens2
    LOGICAL                                       :: print_dens
  END TYPE density_3d_type

  TYPE pair_correl_type
    INTEGER                                       :: conf_counter
    INTEGER                                       :: nr_bins
    REAL(KIND=dp)                                 :: step_lenght
    TYPE(atom_pairs_type), DIMENSION(:), POINTER  :: pairs
    REAL(KIND=dp), DIMENSION(:,:), POINTER        :: g_r
  END TYPE pair_correl_type

  TYPE atom_pairs_type
    CHARACTER(LEN=default_string_length)          :: f_n = ""
    CHARACTER(LEN=default_string_length)          :: s_n = ""
    INTEGER                                       :: pair_count
  END TYPE atom_pairs_type

  TYPE dipole_moment_type
    INTEGER                                       :: conf_counter
    TYPE(tmc_atom_type), DIMENSION(:), POINTER    :: charges_inp
    REAL(KIND=dp), DIMENSION(:), POINTER          :: charges
    REAL(KIND=dp), DIMENSION(:), POINTER          :: last_dip_cl
    LOGICAL                                       :: print_cl_dip = .TRUE.
  END TYPE dipole_moment_type

  TYPE dipole_analysis_type
    REAL(KIND=dp)                                 :: conf_counter
    INTEGER                                       :: ana_type
    LOGICAL                                       :: print_diel_const_traj = .TRUE.
    ! squared dipoles per volume
    REAL(KIND=dp)                                 :: mu2_pv_s
    ! dipole per square root ov volume per direction
    REAL(KIND=dp), DIMENSION(:), POINTER          :: mu_psv, mu_pv, mu2_pv
    ! dipole dipole correlation matrix (per volume)
    REAL(KIND=dp), DIMENSION(:,:), POINTER        :: mu2_pv_mat
    
  END TYPE dipole_analysis_type

  TYPE displacement_type
    INTEGER                                       :: conf_counter
    REAL(KIND=dp)                                 :: disp
    LOGICAL                                       :: print_disp
  END TYPE displacement_type

CONTAINS 

! *****************************************************************************
!> \brief creates a new structure environment for TMC analysis
!> \param tmc_ana structure with parameters for TMC analysis
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_env_create(tmc_ana,error)
    TYPE(tmc_analysis_env), POINTER          :: tmc_ana
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_ana_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(tmc_ana),cp_failure_level,routineP,error,failure)

    ALLOCATE(tmc_ana,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY(tmc_ana%dirs, tmc_ana%cell, tmc_ana%atoms, tmc_ana%density_3d, &
            tmc_ana%pair_correl, tmc_ana%dip_mom, tmc_ana%last_elem, &
            tmc_ana%displace, tmc_ana%dip_ana)

    ! file IO units set 
    tmc_ana%id_traj = -1
    tmc_ana%id_cell = -1
    tmc_ana%id_frc  = -1
    tmc_ana%id_dip  = -1
    tmc_ana%id_ener = -1
    ! line counters
    tmc_ana%lc_traj = 0
    tmc_ana%lc_cell = 0
    tmc_ana%lc_frc  = 0
    tmc_ana%lc_dip  = 0
    tmc_ana%lc_ener = 0
    ! file names
    tmc_ana%out_file_prefix       = ""
    tmc_ana%costum_pos_file_name  = ""
    tmc_ana%costum_dip_file_name  = ""
    tmc_ana%costum_cell_file_name = ""
    tmc_ana%restart               = .TRUE.
    tmc_ana%restarted             = .FALSE.
    tmc_ana%print_test_output     = .FALSE.

    ! configuration control
    tmc_ana%conf_offset = 0
    tmc_ana%from_elem = -1
    tmc_ana%to_elem   = -1
  END SUBROUTINE tmc_ana_env_create

! *****************************************************************************
!> \brief releases the structure environment for TMC analysis
!> \param tmc_ana structure with parameters for TMC analysis
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_env_release(tmc_ana,error)
    TYPE(tmc_analysis_env), POINTER          :: tmc_ana
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_ana_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPostcondition(ASSOCIATED(tmc_ana),cp_failure_level,routineP,error,failure)

    IF(ASSOCIATED(tmc_ana%dirs)) &
      DEALLOCATE(tmc_ana%dirs)

    IF(ASSOCIATED(tmc_ana%density_3d)) &
      CALL tmc_ana_dens_release(tmc_ana%density_3d, error)
    IF(ASSOCIATED(tmc_ana%pair_correl)) &
      CALL tmc_ana_pair_correl_release(tmc_ana%pair_correl, error)

    IF(ASSOCIATED(tmc_ana%dip_mom)) &
      CALL tmc_ana_dipole_moment_release(tmc_ana%dip_mom, error)

    IF(ASSOCIATED(tmc_ana%dip_ana)) &
      CALL tmc_ana_dipole_analysis_release(tmc_ana%dip_ana, error)

    IF(ASSOCIATED(tmc_ana%displace)) &
      CALL tmc_ana_displacement_release(ana_disp=tmc_ana%displace, error=error)

    DEALLOCATE(tmc_ana,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE tmc_ana_env_release

  !============================================================================
  ! density calculations
  !============================================================================

! *****************************************************************************
!> \brief creates a new structure environment for TMC analysis
!> \param ana_dens structure with parameters for TMC density analysis
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_density_create(ana_dens, nr_bins, error)
    TYPE(density_3d_type), POINTER           :: ana_dens
    INTEGER, DIMENSION(3)                    :: nr_bins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_ana_density_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(ana_dens),cp_failure_level,routineP,error,failure)

    ALLOCATE(ana_dens,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY(ana_dens%sum_density, ana_dens%sum_dens2)

    ana_dens%nr_bins(:)   = nr_bins(:)
    ana_dens%conf_counter = 0
    ana_dens%sum_vol      = 0.0_dp
    ana_dens%sum_box_length(:) = 0.0_dp
    ana_dens%sum_vol2     = 0.0_dp
    ana_dens%sum_box_length2(:) = 0.0_dp

    ALLOCATE(ana_dens%sum_density(nr_bins(1), nr_bins(2), nr_bins(3)))
    ALLOCATE(ana_dens%sum_dens2(nr_bins(1), nr_bins(2), nr_bins(3)))
    ana_dens%sum_density(:,:,:) = 0.0_dp
    ana_dens%sum_dens2(:,:,:) = 0.0_dp
    ana_dens%print_dens = .TRUE.
  END SUBROUTINE tmc_ana_density_create

! *****************************************************************************
!> \brief releases the structure environment for TMC analysis
!> \param ana_dens structure with parameters for TMC analysis
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_dens_release(ana_dens,error)
    TYPE(density_3d_type), POINTER           :: ana_dens
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_ana_dens_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPostcondition(ASSOCIATED(ana_dens),cp_failure_level,routineP,error,failure)

    DEALLOCATE(ana_dens%sum_density)
    DEALLOCATE(ana_dens%sum_dens2)
    DEALLOCATE(ana_dens,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE tmc_ana_dens_release

  !============================================================================
  ! radial distribution function
  !============================================================================

! *****************************************************************************
!> \brief creates a new structure environment for TMC analysis
!> \param
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_pair_correl_create(ana_pair_correl, nr_bins, error)
    TYPE(pair_correl_type), POINTER          :: ana_pair_correl
    INTEGER                                  :: nr_bins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_ana_pair_correl_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(ana_pair_correl),cp_failure_level,routineP,error,failure)
    ALLOCATE(ana_pair_correl,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY(ana_pair_correl%g_r, ana_pair_correl%pairs)
    ana_pair_correl%conf_counter = 0
    ana_pair_correl%nr_bins = nr_bins
    ana_pair_correl%step_lenght = -1.0_dp
    ! still the initialization routine has to be called
  END SUBROUTINE tmc_ana_pair_correl_create

! *****************************************************************************
!> \brief releases the structure environment for TMC analysis
!> \param 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_pair_correl_release(ana_pair_correl,error)
    TYPE(pair_correl_type), POINTER          :: ana_pair_correl
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_ana_pair_correl_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPostcondition(ASSOCIATED(ana_pair_correl),cp_failure_level,routineP,error,failure)

    DEALLOCATE(ana_pair_correl%g_r)
    DEALLOCATE(ana_pair_correl%pairs)
    DEALLOCATE(ana_pair_correl, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE tmc_ana_pair_correl_release

! *****************************************************************************
!> \brief search the pair of two atom types in list
!> \param n1,n2 atom names
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  FUNCTION search_pair_in_list(pair_list, n1, n2, list_end, error) RESULT(ind)
    TYPE(atom_pairs_type), DIMENSION(:), &
      POINTER                                :: pair_list
    CHARACTER(LEN=default_string_length)     :: n1, n2
    INTEGER, OPTIONAL                        :: list_end
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: ind

    CHARACTER(LEN=*), PARAMETER :: routineN = 'search_pair_in_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last, list_nr
    LOGICAL                                  :: failure

    CPPrecondition(ASSOCIATED(pair_list),cp_failure_level,routineP,error,failure)
    IF(PRESENT(list_end)) THEN
      CPPrecondition(list_end.LE.SIZE(pair_list),cp_failure_level,routineP,error,failure)
      last = list_end
    ELSE
      last = SIZE(pair_list)
    END IF

    ind = -1

    list_search:DO list_nr=1, last
      IF((pair_list(list_nr)%f_n.EQ.n1 .AND. &
          pair_list(list_nr)%s_n.EQ.n2) .OR. &
         (pair_list(list_nr)%f_n.EQ.n2 .AND. &
          pair_list(list_nr)%s_n.EQ.n1)) THEN
        ind = list_nr
        EXIT list_search
      END IF
    END DO list_search
  END FUNCTION search_pair_in_list

  !============================================================================
  ! classical cell dipole moment
  !============================================================================

! *****************************************************************************
!> \brief creates a new structure environment for TMC analysis
!> \param
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_dipole_moment_create(ana_dip_mom, charge_atm, charge, &
                                          dim_per_elem, error)
    TYPE(dipole_moment_type), POINTER        :: ana_dip_mom
    CHARACTER(LEN=default_string_length), &
      POINTER                                :: charge_atm(:)
    REAL(KIND=dp), POINTER                   :: charge(:)
    INTEGER                                  :: dim_per_elem
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_ana_dipole_moment_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(ana_dip_mom),cp_failure_level,routineP,error,failure)
    ALLOCATE(ana_dip_mom,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY(ana_dip_mom%charges_inp, ana_dip_mom%charges)

    ALLOCATE(ana_dip_mom%charges_inp(SIZE(charge)))
    DO i=1, SIZE(charge)
      ana_dip_mom%charges_inp(i)%name = charge_atm(i)
      ana_dip_mom%charges_inp(i)%mass = charge(i)
    END DO

    ana_dip_mom%conf_counter = 0
    ALLOCATE(ana_dip_mom%last_dip_cl(dim_per_elem))
    ! still the initialization routine has to be called

  END SUBROUTINE tmc_ana_dipole_moment_create

! *****************************************************************************
!> \brief releases the structure environment for TMC analysis
!> \param 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_dipole_moment_release(ana_dip_mom, error)
    TYPE(dipole_moment_type), POINTER        :: ana_dip_mom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'tmc_ana_dipole_moment_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPostcondition(ASSOCIATED(ana_dip_mom),cp_failure_level,routineP,error,failure)

    IF(ASSOCIATED(ana_dip_mom%charges_inp)) DEALLOCATE(ana_dip_mom%charges_inp)
    IF(ASSOCIATED(ana_dip_mom%charges))     DEALLOCATE(ana_dip_mom%charges)
    DEALLOCATE(ana_dip_mom%last_dip_cl)
    DEALLOCATE(ana_dip_mom, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE tmc_ana_dipole_moment_release

! *****************************************************************************
!> \brief creates a new structure environment for TMC analysis
!> \param
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_dipole_analysis_create(ana_dip_ana, error)
    TYPE(dipole_analysis_type), POINTER      :: ana_dip_ana
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'tmc_ana_dipole_analysis_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(ana_dip_ana),cp_failure_level,routineP,error,failure)
    ALLOCATE(ana_dip_ana,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ana_dip_ana%conf_counter = 0
    ana_dip_ana%ana_type     = -1
    ana_dip_ana%mu2_pv_s     = 0.0_dp
    ALLOCATE(ana_dip_ana%mu_psv(3))
    ana_dip_ana%mu_psv(:)    = 0.0_dp
    ALLOCATE(ana_dip_ana%mu_pv(3))
    ana_dip_ana%mu_pv(:)     = 0.0_dp
    ALLOCATE(ana_dip_ana%mu2_pv(3))
    ana_dip_ana%mu2_pv(:)    = 0.0_dp
    ALLOCATE(ana_dip_ana%mu2_pv_mat(3,3))
    ana_dip_ana%mu2_pv_mat(:,:) = 0.0_dp
  END SUBROUTINE tmc_ana_dipole_analysis_create

! *****************************************************************************
!> \brief releases the structure environment for TMC analysis
!> \param 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_dipole_analysis_release(ana_dip_ana, error)
    TYPE(dipole_analysis_type), POINTER      :: ana_dip_ana
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'tmc_ana_dipole_analysis_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPostcondition(ASSOCIATED(ana_dip_ana),cp_failure_level,routineP,error,failure)

    DEALLOCATE(ana_dip_ana%mu_psv)
    DEALLOCATE(ana_dip_ana%mu_pv)
    DEALLOCATE(ana_dip_ana%mu2_pv)
    DEALLOCATE(ana_dip_ana%mu2_pv_mat)

    DEALLOCATE(ana_dip_ana, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE tmc_ana_dipole_analysis_release

  !============================================================================
  ! particle displacement in cell (from one configuration to the next)
  !============================================================================

! *****************************************************************************
!> \brief creates a new structure environment for TMC analysis
!> \param
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_displacement_create(ana_disp, dim_per_elem, error)
    TYPE(displacement_type), POINTER         :: ana_disp
    INTEGER                                  :: dim_per_elem
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_ana_displacement_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(ana_disp),cp_failure_level,routineP,error,failure)
    CPPrecondition(dim_per_elem.GT.0,cp_failure_level,routineP,error,failure)

    ALLOCATE(ana_disp,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ana_disp%disp = 0.0_dp
    ana_disp%conf_counter = 0
    ana_disp%print_disp = .TRUE.
  END SUBROUTINE tmc_ana_displacement_create

! *****************************************************************************
!> \brief releases a structure environment for TMC analysis
!> \param
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes 02.2013
! *****************************************************************************
  SUBROUTINE tmc_ana_displacement_release(ana_disp, error)
    TYPE(displacement_type), POINTER         :: ana_disp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tmc_ana_displacement_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(ana_disp),cp_failure_level,routineP,error,failure)

    DEALLOCATE(ana_disp,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE tmc_ana_displacement_release
END MODULE tmc_analysis_types
