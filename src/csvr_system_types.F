!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Type for the canonical sampling through velocity rescaling
!> \author Teodoro Laino - 09.2007 University of Zurich [tlaino]
! *****************************************************************************
MODULE csvr_system_types
  USE bibliography,                    ONLY: Bussi2007,&
                                             cite_reference
  USE extended_system_types,           ONLY: create_map_info_type,&
                                             map_info_type,&
                                             release_map_info_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_rng_seed,&
                                             rng_stream_type
  USE simpar_types,                    ONLY: simpar_type
  USE string_utilities,                ONLY: compress
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: csvr_system_type,&
            csvr_thermo_type,&
            csvr_init,&
            csvr_dealloc,&
            csvr_thermo_create

! *****************************************************************************
  TYPE csvr_thermo_type
     INTEGER                                 :: degrees_of_freedom
     REAL(KIND=dp)                           :: nkt
     REAL(KIND=dp)                           :: thermostat_energy
     REAL(KIND=dp)                           :: region_kin_energy
     TYPE(rng_stream_type), POINTER          :: gaussian_rng_stream
  END TYPE csvr_thermo_type

! *****************************************************************************
  TYPE csvr_system_type
     INTEGER                                 :: region, glob_num_csvr, loc_num_csvr
     REAL(KIND=dp)                           :: tau_csvr, dt_fact
     TYPE(csvr_thermo_type), POINTER         :: nvt(:)
     TYPE(map_info_type), POINTER            :: map_info
  END TYPE csvr_system_type

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'csvr_system_types'

CONTAINS

! *****************************************************************************
!> \brief Initialize type for Canonical Sampling through Velocity Rescaling (CSVR)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
  SUBROUTINE csvr_init(csvr, simpar, section, error)
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_init', &
      routineP = moduleN//':'//routineN

    NULLIFY(csvr%nvt)
    NULLIFY(csvr%map_info)
    csvr%loc_num_csvr=0
    csvr%glob_num_csvr=0
    csvr%dt_fact=1.0_dp
    CALL cite_reference(Bussi2007)
    CALL section_vals_val_get(section,"TIMECON",r_val=csvr%tau_csvr,error=error)
    ! The CSVR library expects the tau_csv to be in unit of integration timestep
    ! if applied once.. divided by two if the process is applied both to the first
    ! and the second verlet step
    csvr%tau_csvr = csvr%tau_csvr/(0.5_dp*simpar%dt)
    CALL create_map_info_type(csvr%map_info, error)

  END SUBROUTINE csvr_init

! *****************************************************************************
!> \brief Initialize NVT type for CSVR thermostat
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
  SUBROUTINE csvr_thermo_create(csvr, error)
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_thermo_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=40)                        :: name
    INTEGER                                  :: i, ithermo, my_index, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: seed
    REAL(KIND=dp), DIMENSION(3, 2)           :: initial_seed, my_seed

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(csvr),cp_fatal_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(csvr%nvt),cp_fatal_level,routineP,error,failure)

    ALLOCATE ( csvr%nvt(csvr%loc_num_csvr),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, csvr%loc_num_csvr
       csvr%nvt(i)%thermostat_energy = 0.0_dp
       NULLIFY(csvr%nvt(i)%gaussian_rng_stream)
    END DO
    ! Initialize the gaussian stream random number
    ALLOCATE (seed(3,2,csvr%glob_num_csvr),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    initial_seed = next_rng_seed(error=error)

    seed(:,:,1) = initial_seed
    DO ithermo=2,csvr%glob_num_csvr
       seed(:,:,ithermo) = next_rng_seed(seed(:,:,ithermo-1),error=error)
    END DO
    ! Update initial seed
    initial_seed = next_rng_seed(seed(:,:,csvr%glob_num_csvr),error=error)
    DO ithermo = 1, csvr%loc_num_csvr
       my_index = csvr%map_info%index(ithermo)
       my_seed  = seed(:,:,my_index)
       WRITE (UNIT=name,FMT="(A,I8)") "Wiener process for Thermostat #",my_index
       CALL compress(name)
       CALL create_rng_stream(rng_stream=csvr%nvt(ithermo)%gaussian_rng_stream,&
            name=name,distribution_type=GAUSSIAN, extended_precision=.TRUE.,&
            seed=my_seed,error=error)
    END DO
    DEALLOCATE (seed,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE csvr_thermo_create

! *****************************************************************************
!> \brief Deallocate type for CSVR thermostat
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
  SUBROUTINE csvr_dealloc ( csvr, error )
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(csvr)) THEN
       CALL csvr_thermo_dealloc(csvr%nvt, error)
       CALL release_map_info_type(csvr%map_info, error)
       DEALLOCATE (csvr, STAT=stat)
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    ENDIF

  END SUBROUTINE csvr_dealloc

! *****************************************************************************
!> \brief Deallocate NVT type for CSVR thermostat
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] 10.2007- University of Zurich
! *****************************************************************************
  SUBROUTINE csvr_thermo_dealloc ( nvt, error )
    TYPE(csvr_thermo_type), DIMENSION(:), &
      POINTER                                :: nvt
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'csvr_thermo_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(nvt)) THEN
       DO i = 1, SIZE(nvt)
          IF (ASSOCIATED(nvt(i)%gaussian_rng_stream)) THEN
             CALL delete_rng_stream(nvt(i)%gaussian_rng_stream,error=error)
          ENDIF
       END DO
       DEALLOCATE (nvt, STAT=stat)
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    ENDIF
  END SUBROUTINE csvr_thermo_dealloc

END MODULE csvr_system_types

