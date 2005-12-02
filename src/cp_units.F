!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_units *
!!
!!   NAME
!!     cp_units
!!
!!   FUNCTION
!!     unit conversion facility
!!
!!     Units are complex, this module does not try to be very smart, for 
!!     example SI prefixes are not supported automatically, and
!!     which kinds are really basic can change depending on the system of
!!     units choosen, and equivalences are not always catched.
!!
!!     This is thought as a simple conversion facility for the input and output.
!!     If you need something more you are probably better off using the
!!     physcon module directly.
!!
!!   NOTES
!!     One design choice was not to use dynamically allocated elements to
!!     reduce the possibility of leaks.
!!     Needs to be extended (for example charge, dipole,...)
!!     I just added the units and kinds that I needed.
!!     Used by the parser
!!     Should keep an unsorted/uncompressed version for nicer labels?
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     01.2005 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_units
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger, cp_a_l, cp_error_check
  USE kinds, ONLY: dp, default_string_length
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE physcon, ONLY: evolt, kcalmol, joule, wavenum, bohr, a_bohr,&
       femtoseconds, e_mass,massunit,kelvin, seconds, hz, h_bar
  USE mathconstants, only: radians
  USE string_utilities, ONLY: compress, uppercase
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_units'
  INTEGER, SAVE, PRIVATE :: last_unit_id=0, last_unit_set_id=0
  
  INTEGER, PARAMETER, PUBLIC :: cp_ukind_none=0, cp_ukind_energy=1, &
       cp_ukind_length=2, cp_ukind_temperature=3, cp_ukind_angle=4, &
       cp_ukind_pressure=5, cp_ukind_time=6, cp_ukind_mass=7,&
       cp_ukind_max=7

  INTEGER, PARAMETER, PUBLIC ::  cp_units_none=0,cp_units_au=11
  INTEGER, PARAMETER, PUBLIC :: cp_units_m_e=100,cp_units_amu=101, cp_units_kg=102
  INTEGER, PARAMETER, PUBLIC :: cp_units_hartree=120, cp_units_wavenum=121,&
       cp_units_joule=122, cp_units_kcalmol=123, cp_units_Ry=124, cp_units_eV=125
  INTEGER, PARAMETER, PUBLIC :: cp_units_bohr=140, cp_units_angstrom=141,&
       cp_units_m=142
  INTEGER, PARAMETER, PUBLIC ::  cp_units_k=150
  INTEGER, PARAMETER, PUBLIC ::  cp_units_bar=161
  INTEGER, PARAMETER, PUBLIC ::  cp_units_rad=170, cp_units_deg=171
  INTEGER, PARAMETER, PUBLIC ::  cp_units_fs=180, cp_units_s=181, cp_units_wn=182
  
  INTEGER, PARAMETER, PUBLIC :: cp_unit_max_kinds=8, &
       cp_unit_basic_desc_length=15, &
       cp_unit_desc_length=cp_unit_max_kinds*cp_unit_basic_desc_length

  PUBLIC :: cp_unit_type, cp_unit_p_type, cp_unit_set_type
  PUBLIC :: cp_unit_create, cp_unit_retain, cp_unit_release, &
       cp_unit_to_cp2k, cp_unit_from_cp2k, cp_unit_desc,&
       cp_unit_set_create, cp_unit_set_retain, cp_unit_set_release,&
       cp_unit_to_cp2k1, cp_unit_from_cp2k1
  
!!***
!****************************************************************************

  !!****s* cp_units/cp_unit_type *
  !!
  !!   NAME
  !!     cp_unit_type
  !!
  !!   FUNCTION
  !!     stores a unit
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - ref_count: the reference count (see doc/ReferenceCounting.html)
  !!     - id_nr: identification number (unique for each instance)
  !!     - kind: the kind of unit (energy, length,...)
  !!     - unit: the actual unit (Joule, eV,...)
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE cp_unit_type
     INTEGER :: id_nr,ref_count,n_kinds
     INTEGER, DIMENSION(cp_unit_max_kinds):: kind_id, unit_id, power
  END TYPE cp_unit_type
  !!***
  !****************************************************************************

  !!****s* cp_units/cp_unit_p_type *
  !!
  !!   NAME
  !!     cp_unit_p_type
  !!
  !!   FUNCTION
  !!     represent a pointer to a unit (to build arrays of pointers)
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - unit: the pointer to the unit
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE cp_unit_p_type
     TYPE(cp_unit_type), pointer :: unit
  END TYPE cp_unit_p_type
  !!***
  !****************************************************************************

  !!****s* cp_units/cp_unit_set_type *
  !!
  !!   NAME
  !!     cp_unit_set_type
  !!
  !!   FUNCTION
  !!     stores the default units to be used
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE cp_unit_set_type
     INTEGER :: id_nr, ref_count
     TYPE(cp_unit_p_type), dimension(cp_ukind_max) :: units
  END TYPE cp_unit_set_type
  !!***
  !****************************************************************************

contains

!!****f* cp_units/cp_unit_create *
!!
!!   NAME
!!     cp_unit_create
!!
!!   FUNCTION
!!     creates a unit parsing a string
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - unit: the unit to initialize
!!     - string: the string containing the description of the unit
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE cp_unit_create(unit, string, error)
    TYPE(cp_unit_type), pointer :: unit
    CHARACTER(len=*), intent(in) :: string
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    INTEGER, DIMENSION(cp_unit_max_kinds) :: unit_id, kind_id, power
    INTEGER :: i_low, i_high, len_string, i_unit, next_power
    character(len=*), parameter :: routineN='cp_unit_create',&
         routineP=moduleN//':'//routineN
    character(default_string_length) :: desc
    character(LEN=40)                :: formatstr

    failure=.FALSE.
    unit_id=0
    kind_id=0
    power=0
    i_low=1
    i_high=1
    len_string=LEN(string)
    i_unit=0
    next_power=1
    DO WHILE(i_low<len_string)
       IF (string(i_low:i_low)/=' ') EXIT
       i_low=i_low+1
    END DO
    i_high=i_low
    DO WHILE(i_high<=len_string)
       IF ( string(i_high:i_high)==' '.OR.string(i_high:i_high)=='^'.OR.&
            string(i_high:i_high)=='*'.OR.string(i_high:i_high)=='/') EXIT
       i_high=i_high+1
    END DO
    DO WHILE(.NOT.failure)
       IF (i_high<=i_low.OR.i_low>len_string) EXIT
       i_unit=i_unit+1
       IF (i_unit>cp_unit_max_kinds) THEN
          CALL cp_assert(.FALSE.,cp_failure_level,&
               cp_assertion_failed,routineP,&
               "maximum number of combined units exceded",&
               error,failure)
          EXIT
       END IF
       ! read unit
       SELECT CASE(string(i_low:i_high-1))
       CASE("hartree")
          unit_id(i_unit)=cp_units_hartree
          kind_id(i_unit)=cp_ukind_energy
       CASE("wavenumber_e")
          unit_id(i_unit)=cp_units_wavenum
          kind_id(i_unit)=cp_ukind_energy
       CASE("joule")
          unit_id(i_unit)=cp_units_joule
          kind_id(i_unit)=cp_ukind_energy
       CASE("kcalmol")
          unit_id(i_unit)=cp_units_kcalmol
          kind_id(i_unit)=cp_ukind_energy
       CASE("Ry")
          unit_id(i_unit)=cp_units_Ry
          kind_id(i_unit)=cp_ukind_energy
       CASE("eV")
          unit_id(i_unit)=cp_units_eV
          kind_id(i_unit)=cp_ukind_energy
       CASE("K_e")
          unit_id(i_unit)=cp_units_k
          kind_id(i_unit)=cp_ukind_energy
       CASE("energy")
          unit_id(i_unit)=cp_units_none
          kind_id(i_unit)=cp_ukind_energy
       CASE("au_l")
          unit_id(i_unit)=cp_units_au
          kind_id(i_unit)=cp_ukind_length
       CASE("bohr")
          unit_id(i_unit)=cp_units_bohr
          kind_id(i_unit)=cp_ukind_length
       CASE("m")
          unit_id(i_unit)=cp_units_m
          kind_id(i_unit)=cp_ukind_length
       CASE("angstrom")
          unit_id(i_unit)=cp_units_angstrom
          kind_id(i_unit)=cp_ukind_length
       CASE ("length")
          unit_id(i_unit)=cp_units_none
          kind_id(i_unit)=cp_ukind_length
       CASE("K","K_temp")
          unit_id(i_unit)=cp_units_k
          kind_id(i_unit)=cp_ukind_temperature
       CASE("au_temp")
          unit_id(i_unit)=cp_units_au
          kind_id(i_unit)=cp_ukind_temperature
       CASE("temperature")
          unit_id(i_unit)=cp_units_none
          kind_id(i_unit)=cp_ukind_temperature
       CASE("bar")
          unit_id(i_unit)=cp_units_bar
          kind_id(i_unit)=cp_ukind_pressure
       CASE("au_p")
          unit_id(i_unit)=cp_units_au
          kind_id(i_unit)=cp_ukind_pressure
       CASE("pressure")
          unit_id(i_unit)=cp_units_none
          kind_id(i_unit)=cp_ukind_pressure
       CASE("rad")
          unit_id(i_unit)=cp_units_rad
          kind_id(i_unit)=cp_ukind_angle
       CASE("deg")
          unit_id(i_unit)=cp_units_deg
          kind_id(i_unit)=cp_ukind_angle
       CASE("angle")
          unit_id(i_unit)=cp_units_none
          kind_id(i_unit)=cp_ukind_angle
       CASE("s")
          unit_id(i_unit)=cp_units_s
          kind_id(i_unit)=cp_ukind_time
       CASE("fs")
          unit_id(i_unit)=cp_units_fs
          kind_id(i_unit)=cp_ukind_time
       CASE("wavenumber_t")
          unit_id(i_unit)=cp_units_wn
          kind_id(i_unit)=cp_ukind_time          
       CASE("au_t")
          unit_id(i_unit)=cp_units_au
          kind_id(i_unit)=cp_ukind_time
       CASE("time")
          unit_id(i_unit)=cp_units_none
          kind_id(i_unit)=cp_ukind_time
       CASE("kg")
          unit_id(i_unit)=cp_units_kg
          kind_id(i_unit)=cp_ukind_mass
       CASE("amu")
          unit_id(i_unit)=cp_units_amu
          kind_id(i_unit)=cp_ukind_mass
       CASE("m_e")
          unit_id(i_unit)=cp_units_m_e
          kind_id(i_unit)=cp_ukind_mass
       CASE("au_m")
          unit_id(i_unit)=cp_units_au
          kind_id(i_unit)=cp_ukind_mass
       CASE("mass")
          unit_id(i_unit)=cp_units_none
          kind_id(i_unit)=cp_ukind_mass
       CASE("au")
          CALL cp_assert(.false.,cp_failure_level,&
               cp_assertion_failed,routineP,&
               "au unit without specifing its kind not accepted, use (au_e, au_t, au_temp, au_l)",&
               error,failure)
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown unit:"//string(i_low:i_high-1),&
               error,failure)
       END SELECT
       power(i_unit)=next_power
       ! parse op
       i_low=i_high
       DO WHILE(i_low<=len_string)
          IF (string(i_low:i_low)/=' ') EXIT
          i_low=i_low+1
       END DO
       i_high=i_low
       DO WHILE(i_high<=len_string)
          IF ( string(i_high:i_high)==' '.OR.string(i_high:i_high)=='^'.OR.&
               string(i_high:i_high)=='*'.OR.string(i_high:i_high)=='/') EXIT
          i_high=i_high+1
       END DO
       IF (i_high<i_low.OR.i_low>len_string) EXIT
       IF (string(i_low:i_high)=='^') THEN
          i_low=i_high+1
          DO WHILE(i_low<=len_string)
             IF (string(i_low:i_low)/=' ') EXIT
             i_low=i_low+1
          END DO
          i_high=i_low
          DO WHILE(i_high<=len_string)
             SELECT CASE(string(i_high:i_high))
             CASE('+','-','0','1','2','3','4','5','6','7','8','9')
                i_high=i_high+1
             CASE default
                EXIT
             END SELECT
          END DO
          IF (i_high<=i_low.OR.i_low>len_string) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,&
                  cp_assertion_failed,routineP,&
                  "an integer number is expected after a '^'",&
                  error,failure)
             EXIT
          END IF
          formatstr="(i"//cp_to_string(i_high-i_low+1)//")"
          READ (string(i_low:i_high-1),formatstr)&
               next_power
          power(i_unit)=power(i_unit)*next_power
          ! next op
          i_low=i_high
          DO WHILE(i_low<len_string)
             IF (string(i_low:i_low)/=' ') EXIT
             i_low=i_low+1
          END DO
          i_high=i_low
          DO WHILE(i_high<=len_string)
             IF ( string(i_high:i_high)==' '.OR.string(i_high:i_high)=='^'.OR.&
                  string(i_high:i_high)=='*'.OR.string(i_high:i_high)=='/') EXIT
             i_high=i_high+1
          END DO
       END IF
       IF (i_low>len_string) EXIT
       next_power=1
       IF (string(i_low:i_high)=="*".OR.string(i_low:i_high)=='/') THEN
          IF (string(i_low:i_high)=='/') next_power=-1
          i_low=i_high+1
          DO WHILE(i_low<=len_string)
             IF (string(i_low:i_low)/=' ') EXIT
             i_low=i_low+1
          END DO
          i_high=i_low
          DO WHILE(i_high<=len_string)
             IF ( string(i_high:i_high)==' '.OR.string(i_high:i_high)=='^'.OR.&
                  string(i_high:i_high)=='*'.OR.string(i_high:i_high)=='/') EXIT
             i_high=i_high+1
          END DO
       END IF
    END DO
    CALL cp_unit_create2(unit,kind_id=kind_id, unit_id=unit_id, &
         power=power, error=error)
    desc=cp_unit_desc(unit,error=error)
  END SUBROUTINE cp_unit_create
!***************************************************************************

!!****f* cp_units/cp_unit_create2 *
!!
!!   NAME
!!     cp_unit_create2
!!
!!   FUNCTION
!!     creates and initializes the given unit of mesure (performs some error
!!     check)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - unit: the unit descriptor to be initialized
!!     - kind_id: the kind of unit (length,energy,...), use the constants 
!!       cp_ukind_*
!!     - unit_id: the actual unit (use constants cp_units_*)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_unit_create2(unit, kind_id, unit_id, power, error)
  TYPE(cp_unit_type), pointer :: unit
  INTEGER, INTENT(in), DIMENSION(:) :: kind_id, unit_id
  INTEGER, INTENT(in), DIMENSION(:), OPTIONAL :: power
  type(cp_error_type), optional, intent(inout) :: error
  
  LOGICAL :: failure,repeat
  INTEGER :: stat, i, max_kind, max_pos,j
  character(len=*), parameter :: routineN='cp_unit_init',&
        routineP=moduleN//':'//routineN

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(unit),cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(kind_id)<=cp_unit_max_kinds,cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(unit_id)<=cp_unit_max_kinds,cp_failure_level,routineP,error,failure)
  ALLOCATE(unit,stat=stat)
  CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
  unit%ref_count=1
  last_unit_id=last_unit_id+1
  unit%id_nr=last_unit_id
  unit%kind_id(1:SIZE(kind_id))=kind_id
  unit%kind_id(SIZE(kind_id)+1:)=0
  unit%unit_id(1:SIZE(unit_id))=unit_id
  unit%unit_id(SIZE(unit_id):)=0
  IF (PRESENT(power)) THEN
     unit%power(1:size(power))=power
     unit%power(SIZE(power)+1:)=0
     DO i=1,SIZE(unit%power)
        IF (unit%power(i)==0) THEN
           unit%kind_id(i)=0
           unit%unit_id(i)=0
        END IF
     END DO
  ELSE
     DO i=1,SIZE(unit%power)
        IF (unit%unit_id(i)/=0) THEN
           unit%power(i)=1
        ELSE
           unit%power(i)=0
        END IF
     END DO
  END IF
  
  ! remove unnecessary units
  ! reorder & compress
  unit%n_kinds=0
  DO i=1,SIZE(unit%kind_id)
     ! find max and compress in the rest
     DO
        max_kind=unit%kind_id(i)
        max_pos=i
        repeat=.false.
        DO j=i+1,SIZE(unit%kind_id)
           IF (unit%kind_id(j)>=max_kind) THEN
              IF (unit%kind_id(j)/=0.AND.unit%kind_id(j)==max_kind.AND.&
                   unit%unit_id(j)==unit%unit_id(max_pos)) THEN
                 unit%power(max_pos)=unit%power(max_pos)+unit%power(j)
                 unit%kind_id(j)=0
                 unit%unit_id(j)=0
                 unit%power(j)=0
                 IF (unit%power(max_pos)==0) THEN
                    unit%kind_id(max_pos)=0
                    unit%unit_id(max_pos)=0
                    unit%power(max_pos)=0
                    repeat=.TRUE.
                    exit
                 END IF
              ELSE IF (unit%kind_id(j)>max_kind.OR.&
                   (unit%kind_id(j)==max_kind.AND.&
                   unit%unit_id(j)>unit%unit_id(max_pos))) THEN
                 max_kind=unit%kind_id(j)
                 max_pos=j
              END IF
           END IF
        END DO
        if (.not.repeat) exit
     END DO
     IF (max_kind/=0) unit%n_kinds=unit%n_kinds+1
     ! put the max at pos i
     IF (max_pos/=i) THEN
        unit%kind_id(max_pos)=unit%kind_id(i)
        unit%kind_id(i)=max_kind
        max_kind=unit%unit_id(max_pos)
        unit%unit_id(max_pos)=unit%unit_id(i)
        unit%unit_id(i)=max_kind
        max_kind=unit%power(max_pos)
        unit%power(max_pos)=unit%power(i)
        unit%power(i)=max_kind
     END IF
     ! check unit
     failure=failure.OR..NOT.cp_basic_unit_check(basic_kind=unit%kind_id(i),&
          basic_unit=unit%unit_id(i),error=error)
  END DO

END SUBROUTINE cp_unit_create2
!***************************************************************************

!!****f* cp_units/cp_unit_retain *
!!
!!   NAME
!!     cp_unit_retain
!!
!!   FUNCTION
!!     retains the given unit
!!
!!   NOTES
!!     at the moment not needed, there for completeness
!!
!!   INPUTS
!!     - unit: the unit to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_unit_retain(unit,error)
  TYPE(cp_unit_type), pointer :: unit
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_unit_dealloc_ref',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(ASSOCIATED(unit),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(unit%ref_count>0,cp_failure_level,routineP,error)
     unit%ref_count=unit%ref_count+1
  END IF
END SUBROUTINE cp_unit_retain
!***************************************************************************

!!****f* cp_units/cp_unit_release *
!!
!!   NAME
!!     cp_unit_release
!!
!!   FUNCTION
!!     releases the given unit
!!
!!   NOTES
!!     at the moment not needed, there for completeness
!!
!!   INPUTS
!!     - unit: the unit to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_unit_release(unit,error)
  TYPE(cp_unit_type), pointer :: unit
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_unit_dealloc_ref',&
        routineP=moduleN//':'//routineN
  integer :: stat

  failure=.false.
  
  IF (associated(unit)) THEN
     CPPreconditionNoFail(unit%ref_count>0,cp_failure_level,routineP,error)
     unit%ref_count=unit%ref_count-1
     IF (unit%ref_count==0) THEN
        DEALLOCATE(unit,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(unit)
END SUBROUTINE cp_unit_release
!***************************************************************************

!!****f* cp_units/cp_basic_unit_check *
!!
!!   NAME
!!     cp_basic_unit_check
!!
!!   FUNCTION
!!     controls that the kind and contains meaningful information
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - basic_kind: the kind of the unit
!!     - basic_unit: the unit to check
!!     - error_level: error level of the errors due to invalid values.
!!       Defaults to cp_failure_level
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_basic_unit_check(basic_kind,basic_unit,error_level,error)&
     result(res)
  INTEGER, INTENT(in) :: basic_kind, basic_unit
  INTEGER, INTENT(in), OPTIONAL :: error_level
  type(cp_error_type), optional, intent(inout) :: error
  logical :: res

  LOGICAL :: failure
  integer :: my_error_level
  character(len=*), parameter :: routineN='cp_basic_unit_check',&
       routineP=moduleN//':'//routineN

  failure=.false.
  my_error_level=cp_failure_level
  if (present(error_level)) my_error_level=error_level

  IF (.NOT. failure) THEN
     SELECT CASE(basic_kind)
     CASE(cp_ukind_energy)
        SELECT CASE (basic_unit)
        CASE(cp_units_hartree, cp_units_wavenum, cp_units_joule, cp_units_kcalmol,&
             cp_units_Ry, cp_units_eV,cp_units_au,cp_units_k, cp_units_none)
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown energy unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_length)
        SELECT CASE (basic_unit)
        CASE(cp_units_bohr, cp_units_angstrom, cp_units_au,cp_units_none, cp_units_m)
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown length unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_temperature)
        SELECT CASE (basic_unit)
        CASE(cp_units_k,cp_units_au,cp_units_none)
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown temperature unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_pressure)
        SELECT CASE (basic_unit)
        CASE(cp_units_bar, cp_units_au,cp_units_none)
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown pressure unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_angle)
        SELECT CASE (basic_unit)
        CASE(cp_units_rad, cp_units_deg,cp_units_none)
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown angle unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_time)
        SELECT CASE (basic_unit)
        CASE(cp_units_s, cp_units_fs, cp_units_au, cp_units_wn, cp_units_none)
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown time unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_mass)
        SELECT CASE (basic_unit)
        CASE(cp_units_kg, cp_units_amu, cp_units_m_e, cp_units_au, cp_units_none)
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown mass unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_none)
        CALL cp_assert(basic_unit==cp_units_none,my_error_level,&
             cp_assertion_failed, routineP,&
             "if the kind of the unit is none also unit must be undefined,not:"&
             //TRIM(cp_to_string(basic_unit)),&
             error,failure)
     CASE default
        CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
             routineP,"unknown kind of unit:"//TRIM(cp_to_string(basic_kind)),&
             error,failure)
     END SELECT
  END IF
  res=.not.failure
END FUNCTION cp_basic_unit_check
!***************************************************************************

!!****f* cp_units/cp_basic_unit_to_cp2k *
!!
!!   NAME
!!     cp_basic_unit_to_cp2k
!!
!!   FUNCTION
!!     converts a value to the internal cp2k units
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - value: the value to convert
!!     - basic_kind: the kind of the unit of the value
!!     - basic_unit: the unit of the value
!!     - power: the power of the unit (defaults to 1)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_basic_unit_to_cp2k(value,basic_kind,basic_unit,power,error)&
     RESULT(res)
  REAL(kind=dp), INTENT(in) :: value
  INTEGER, INTENT(in) :: basic_kind, basic_unit
  INTEGER, INTENT(in), optional :: power
  type(cp_error_type), optional, intent(inout) :: error
  real(kind=dp) :: res
  
  logical :: failure
  integer :: my_power
  character(len=*), parameter :: routineN='cp_basic_unit_to_cp2k',&
       routineP=moduleN//':'//routineN
  INTEGER, parameter :: my_error_level=cp_failure_level

  failure=.false.
  my_power=1
  IF (PRESENT(power)) my_power=power
  IF (basic_unit==cp_units_none) THEN
     CALL cp_assert(basic_kind==cp_units_none,my_error_level,&
          cp_assertion_failed,routineP,&
          "unit not yet fully specified, unit of kind "//&
          TRIM(cp_to_string(basic_unit)),error,failure)
  END IF
  IF (.NOT.failure) THEN
     SELECT CASE(basic_kind)
     CASE(cp_ukind_energy)
        SELECT CASE (basic_unit)
        CASE(cp_units_hartree,cp_units_au)
           res=value
        CASE(cp_units_wavenum)
           res=wavenum**(-my_power)*value
        CASE(cp_units_joule)
           res=joule**(-my_power)*value
        CASE(cp_units_kcalmol)
           res=kcalmol**(-my_power)*value
        CASE(cp_units_Ry)
           res=0.5_dp**my_power*value
        CASE(cp_units_eV)
           res=evolt**(-my_power)*value
        CASE(cp_units_k)
           res=kelvin**(-my_power)*value
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown energy unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_length)
        SELECT CASE (basic_unit)
        CASE(cp_units_bohr,cp_units_au)
           res=value
        CASE(cp_units_m)
           res=value*(1.e10_dp*bohr)**my_power
        CASE(cp_units_angstrom)
           res=value*bohr**my_power
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown length unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_temperature)
        SELECT CASE (basic_unit)
        CASE(cp_units_k)
           res=kelvin**(-my_power)*value
        CASE(cp_units_au)
           res=value
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown temperature unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_pressure)
        SELECT CASE (basic_unit)
        CASE(cp_units_bar)
           res=(joule/(a_bohr**3*1.0E+05_dp))**(-my_power)*value
        CASE(cp_units_au)
           res=value
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown pressure unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_angle)
        SELECT CASE (basic_unit)
        CASE(cp_units_rad)
           res=value
        CASE(cp_units_deg)
           res=value*(radians)**my_power
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown angle unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_time)
        SELECT CASE (basic_unit)
        CASE(cp_units_s)
           res=value*seconds**(-my_power)
        CASE(cp_units_fs)
           res=value*femtoseconds**(-my_power)
        CASE(cp_units_au)
           res=value
        CASE(cp_units_wn)
           res=(joule * wavenum / hz / h_bar)**(my_power)/value
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown time unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_mass)
        SELECT CASE (basic_unit)
        CASE(cp_units_kg)
           res=e_mass**my_power*value
        CASE(cp_units_amu)
           res=massunit**my_power*value
        CASE(cp_units_m_e,cp_units_au)
           res=value
        CASE default
           CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
                routineP,"unknown mass unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_none)
        CALL cp_assert(.false.,my_error_level,&
             cp_assertion_failed, routineP,&
             "if the kind of the unit is none also unit must be undefined,not:"&
             //TRIM(cp_to_string(basic_unit)),&
             error,failure)        
     CASE default
        CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
             routineP,"unknown kind of unit:"//TRIM(cp_to_string(basic_kind)),&
             error,failure)
     END SELECT
  END IF
END FUNCTION cp_basic_unit_to_cp2k
!***************************************************************************

!!****f* cp_units/cp_basic_unit_desc *
!!
!!   NAME
!!     cp_basic_unit_desc
!!
!!   FUNCTION
!!     returns the label of the current basic unit
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - basic_kind: the kind of the unit of the value
!!     - basic_unit: the unit of the value
!!     - power: the power of the unit (defaults to 1)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_basic_unit_desc(basic_kind,basic_unit,power,&
     accept_undefined,error)&
     RESULT(res)
  INTEGER, INTENT(in) :: basic_kind, basic_unit
  INTEGER, INTENT(in), optional :: power
  LOGICAL, INTENT(in), optional :: accept_undefined
  type(cp_error_type), optional, intent(inout) :: error
  character(len=cp_unit_basic_desc_length) :: res

  LOGICAL :: failure, my_accept_undefined
  INTEGER :: my_power, a
  character(len=*), parameter :: routineN='cp_basic_unit_desc',&
       routineP=moduleN//':'//routineN
  INTEGER, parameter :: my_error_level=cp_failure_level

  failure=.false.
  my_power=1
  res=""
  my_accept_undefined=.false.
  if (accept_undefined) my_accept_undefined=accept_undefined
  IF (PRESENT(power)) my_power=power
  IF (basic_unit==cp_units_none) THEN
     CALL cp_assert(my_accept_undefined.OR.basic_kind/=cp_units_none,&
          my_error_level,cp_assertion_failed,routineP,&
          "unit not yet fully specified, unit of kind "//&
          TRIM(cp_to_string(basic_kind)),error,failure)
  END IF
  IF (.NOT.failure) THEN
     SELECT CASE(basic_kind)
     CASE(cp_ukind_energy)
        SELECT CASE (basic_unit)
        CASE(cp_units_hartree,cp_units_au)
           res="hartree"
        CASE(cp_units_wavenum)
           res="wavenumber_e"
        CASE(cp_units_joule)
           res="joule"
        CASE(cp_units_kcalmol)
           res="kcalmol"
        CASE(cp_units_Ry)
           res="Ry"
        CASE(cp_units_eV)
           res="eV"
        CASE(cp_units_k)
           res="K_e"
        CASE(cp_units_none)
           res="energy"
           CALL cp_assert(my_accept_undefined,my_error_level,&
                cp_assertion_failed,routineP,&
                "unit not yet fully specified, unit of kind "//&
                TRIM(res),error,failure)
        CASE default
           CALL cp_assert(.false.,my_error_level,cp_assertion_failed,&
                routineP,"unknown energy unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_length)
        SELECT CASE (basic_unit)
        CASE(cp_units_bohr,cp_units_au)
           res="bohr"
        CASE(cp_units_m)
           res="m"
        CASE(cp_units_angstrom)
           res="angstrom"
        CASE default
           res="length"
           CALL cp_assert(.false.,my_error_level,cp_assertion_failed,&
                routineP,"unknown length unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_temperature)
        SELECT CASE (basic_unit)
        CASE(cp_units_k)
           res="K"
        CASE(cp_units_au)
           res="au_temp"
        CASE(cp_units_none)
           res="temperature"
           CALL cp_assert(my_accept_undefined,my_error_level,&
                cp_assertion_failed,routineP,&
                "unit not yet fully specified, unit of kind "//&
                TRIM(res),error,failure)
        CASE default
           CALL cp_assert(.false.,my_error_level,cp_assertion_failed,&
                routineP,"unknown temperature unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_pressure)
        SELECT CASE (basic_unit)
        CASE(cp_units_bar)
           res="bar"
        CASE(cp_units_au)
           res="au_p"
        CASE(cp_units_none)
           res="pressure"
           CALL cp_assert(my_accept_undefined,my_error_level,&
                cp_assertion_failed,routineP,&
                "unit not yet fully specified, unit of kind "//&
                TRIM(res),error,failure)
        CASE default
           CALL cp_assert(.false.,my_error_level,cp_assertion_failed,&
                routineP,"unknown pressure unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_angle)
        SELECT CASE (basic_unit)
        CASE(cp_units_rad)
           res="rad"
        CASE(cp_units_deg)
           res="deg"
        CASE(cp_units_none)
           res="angle"
           CALL cp_assert(my_accept_undefined,my_error_level,&
                cp_assertion_failed,routineP,&
                "unit not yet fully specified, unit of kind "//&
                TRIM(res),error,failure)
        CASE default
           CALL cp_assert(.false.,my_error_level,cp_assertion_failed,&
                routineP,"unknown angle unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_time)
        SELECT CASE (basic_unit)
        CASE(cp_units_s)
           res="s"
        CASE(cp_units_fs)
           res="fs"
        CASE(cp_units_au)
           res="au_t"
        CASE(cp_units_wn)
           res="wavenumber_t"
        CASE(cp_units_none)
           res="time"
           CALL cp_assert(my_accept_undefined,my_error_level,&
                cp_assertion_failed,routineP,&
                "unit not yet fully specified, unit of kind "//&
                TRIM(res),error,failure)
        CASE default
           CALL cp_assert(.false.,my_error_level,cp_assertion_failed,&
                routineP,"unknown time unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_mass)
        SELECT CASE (basic_unit)
        CASE(cp_units_kg)
           res="kg"
        CASE(cp_units_amu)
           res="amu"
        CASE(cp_units_m_e,cp_units_au)
           res="m_e"
        CASE(cp_units_none)
           res="mass"
           CALL cp_assert(my_accept_undefined,my_error_level,&
                cp_assertion_failed,routineP,&
                "unit not yet fully specified, unit of kind "//&
                TRIM(res),error,failure)
        CASE default
           CALL cp_assert(.false.,my_error_level,cp_assertion_failed,&
                routineP,"unknown mass unit:"//TRIM(cp_to_string(basic_unit)),&
                error,failure)
        END SELECT
     CASE(cp_ukind_none)
        CALL cp_assert(.false.,my_error_level,&
             cp_assertion_failed, routineP,&
             "if the kind of the unit is none also unit must be undefined,not:"&
             //TRIM(cp_to_string(basic_unit)),&
             error,failure)        
     CASE default
        CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
             routineP,"unknown kind of unit:"//TRIM(cp_to_string(basic_kind)),&
             error,failure)
     END SELECT
     IF (my_power/=1) THEN
        a=LEN_TRIM(res)
        CPPrecondition(LEN(res)-a>=3,cp_failure_level,routineP,error,failure)
        WRITE (res(a+1:),"('^',i3)") my_power
        CALL compress(res,.true.)
     END IF
  END IF
END FUNCTION cp_basic_unit_desc
!***************************************************************************

!!****f* cp_units/cp_unit_desc *
!!
!!   NAME
!!     cp_unit_desc
!!
!!   FUNCTION
!!     returns the "name" of the given unit
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - unit: the unit to describe
!!     - defaults: defaults for the undefined units, optional
!!     - accept_undefined: if defaults is not present or is not associated
!!       whether undefined units should be accepted (defaults to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_unit_desc(unit,defaults,accept_undefined,error)&
     result(res)
  TYPE(cp_unit_type), pointer :: unit
  TYPE(cp_unit_set_type), POINTER, optional :: defaults
  LOGICAL, INTENT(in), optional :: accept_undefined
  type(cp_error_type), optional, intent(inout) :: error
  character(len=cp_unit_desc_length) :: res
  
  LOGICAL :: failure,has_defaults,my_accept_undefined
  INTEGER :: i
  character(len=*), parameter :: routineN='cp_unit_desc',&
        routineP=moduleN//':'//routineN
  INTEGER :: pos, my_unit
  failure=.false.
  
  CPPrecondition(ASSOCIATED(unit),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     res=""
     pos=1
     my_accept_undefined=.false.
     IF (PRESENT(accept_undefined)) my_accept_undefined=accept_undefined
     DO i=1,unit%n_kinds
        CPPrecondition(unit%kind_id(i)/=0,cp_failure_level,routineP,error,failure)
        CPPrecondition(pos<LEN(res),cp_failure_level,routineP,error,failure)
        IF (failure) EXIT
        my_unit=unit%unit_id(i)
        has_defaults=.false.
        IF (PRESENT(defaults)) has_defaults=ASSOCIATED(defaults)
        IF (my_unit==0) THEN
           IF (has_defaults) THEN
              my_unit=defaults%units(unit%kind_id(i))%unit%unit_id(1)
           ELSE
              CPPrecondition(my_accept_undefined,cp_failure_level,routineP,error,failure)
           END IF
        END IF
        res(pos:)=TRIM(cp_basic_unit_desc(basic_kind=unit%kind_id(i),&
             basic_unit=my_unit,accept_undefined=my_accept_undefined,&
             power=unit%power(i),error=error))
        pos=LEN_TRIM(res)+1
     END DO
  END IF
END FUNCTION cp_unit_desc
!***************************************************************************

!!****f* cp_units/cp_unit_to_cp2k *
!!
!!   NAME
!!     cp_unit_to_cp2k
!!
!!   FUNCTION
!!     tranform a value to the internal cp2k units
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - value: the value to convert
!!     - unit: the unit of the result
!!     - defaults: the defaults unit for those that are left free
!!       (cp_units_none)
!!     - power: the power of the unit (defaults to 1)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_unit_to_cp2k1(value,unit,defaults,power,error) RESULT(res)
  REAL(kind=dp), intent(in) :: value
  TYPE(cp_unit_type), pointer :: unit
  TYPE(cp_unit_set_type), POINTER, optional :: defaults
  INTEGER, OPTIONAL, intent(in) :: power
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  REAL(kind=dp) :: res
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_unit_to_cp2k1',&
        routineP=moduleN//':'//routineN
  INTEGER :: my_power,my_basic_unit, i_unit
  
  failure=.false.
  my_power=1
  if (present(power)) my_power=power
  res=value
  CPPrecondition(ASSOCIATED(unit),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i_unit=1,unit%n_kinds
        CPPrecondition(unit%kind_id(i_unit)>0,cp_failure_level,routineP,error,failure)
        my_basic_unit=unit%unit_id(i_unit)
        IF (my_basic_unit==0) THEN
           CPPrecondition(PRESENT(defaults),cp_failure_level,routineP,error,failure)
           IF (failure) EXIT
           CPPrecondition(ASSOCIATED(defaults),cp_failure_level,routineP,error,failure)
           IF (failure) EXIT
           CALL cp_assert(ASSOCIATED(defaults%units(unit%kind_id(i_unit))%unit),&
                cp_failure_level,cp_assertion_failed,routineP,&
                CPSourceFileRef,&
                error,failure)
           IF (failure) EXIT
           my_basic_unit=defaults%units(unit%kind_id(i_unit))%unit%unit_id(1)
        END IF
        res=cp_basic_unit_to_cp2k(value=res,basic_unit=my_basic_unit,&
             basic_kind=unit%kind_id(i_unit),&
             power=my_power*unit%power(i_unit),error=error)
     END DO
  END IF
END FUNCTION cp_unit_to_cp2k1
!***************************************************************************

!!****f* cp_units/cp_unit_from_cp2k *
!!
!!   NAME
!!     cp_unit_from_cp2k
!!
!!   FUNCTION
!!     converts from the internal cp2k units to the given unit
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - value: the value to convert
!!     - unit: the unit of the result
!!     - defaults: the defaults unit for those that are left free
!!       (cp_units_none)
!!     - power: the power of the unit (defaults to 1)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_unit_from_cp2k1(value,unit,defaults,power,error) RESULT(res)
  REAL(kind=dp), intent(in) :: value
  TYPE(cp_unit_type), POINTER :: unit
  TYPE(cp_unit_set_type), POINTER, OPTIONAL :: defaults
  INTEGER, INTENT(in), optional :: power
  type(cp_error_type), optional, intent(inout) :: error
  REAL(kind=dp) :: res

  logical :: failure
  integer :: my_power
  character(len=*), parameter :: routineN='cp_unit_from_cp2k1',&
       routineP=moduleN//':'//routineN

  failure=.false.
  my_power=1
  if (present(power)) my_power=power
  IF (PRESENT(defaults)) THEN
     res=cp_unit_to_cp2k1(value=value,unit=unit,defaults=defaults,&
          power=-my_power,error=error)
  ELSE
     res=cp_unit_to_cp2k1(value=value,unit=unit,power=-my_power,error=error)
  END IF
END FUNCTION cp_unit_from_cp2k1
!***************************************************************************

!!****f* cp_units/cp_unit_to_cp2k *
!!
!!   NAME
!!     cp_unit_to_cp2k
!!
!!   FUNCTION
!!     converts to the internal cp2k units to the given unit
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - value: the value to convert
!!     - unit_str: the unit of the result as string
!!     - defaults: the defaults unit for those that are left free
!!       (cp_units_none)
!!     - power: the power of the unit (defaults to 1)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_unit_to_cp2k(value,unit_str,defaults,power,error) RESULT(res)
  REAL(kind=dp), intent(in) :: value
  CHARACTER(len=*), INTENT(in) :: unit_str
  TYPE(cp_unit_set_type), POINTER, OPTIONAL :: defaults
  INTEGER, INTENT(in), optional :: power
  type(cp_error_type), optional, intent(inout) :: error
  REAL(kind=dp) :: res

  TYPE(cp_unit_type), POINTER :: my_unit
  logical :: failure
  character(len=*), parameter :: routineN='cp_unit_to_cp2k',&
       routineP=moduleN//':'//routineN

  nullify(my_unit)
  CALL cp_unit_create(my_unit,unit_str,error=error)
  IF (PRESENT(defaults)) THEN
     res=cp_unit_to_cp2k1(value=value,unit=my_unit,defaults=defaults,&
          power=power,error=error)
  ELSE
     res=cp_unit_to_cp2k1(value=value,unit=my_unit,power=power,error=error)
  END IF
  CALL cp_unit_release(my_unit,error=error)
  
END FUNCTION cp_unit_to_cp2k
!***************************************************************************

!!****f* cp_units/cp_unit_from_cp2k *
!!
!!   NAME
!!     cp_unit_from_cp2k
!!
!!   FUNCTION
!!     converts from the internal cp2k units to the given unit
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - value: the value to convert
!!     - unit_str: the unit of the result as string
!!     - defaults: the defaults unit for those that are left free
!!       (cp_units_none)
!!     - power: the power of the unit (defaults to 1)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_unit_from_cp2k(value,unit_str,defaults,power,error) RESULT(res)
  REAL(kind=dp), intent(in) :: value
  CHARACTER(len=*), INTENT(in) :: unit_str
  TYPE(cp_unit_set_type), POINTER, OPTIONAL :: defaults
  INTEGER, INTENT(in), optional :: power
  type(cp_error_type), optional, intent(inout) :: error
  REAL(kind=dp) :: res

  TYPE(cp_unit_type), POINTER :: my_unit
  logical :: failure
  character(len=*), parameter :: routineN='cp_unit_from_cp2k',&
       routineP=moduleN//':'//routineN

  nullify(my_unit)
  CALL cp_unit_create(my_unit,unit_str,error=error)
  IF (PRESENT(defaults)) THEN
     res=cp_unit_from_cp2k1(value=value,unit=my_unit,defaults=defaults,&
          power=power,error=error)
  ELSE
     res=cp_unit_from_cp2k1(value=value,unit=my_unit,power=power,error=error)
  END IF
  CALL cp_unit_release(my_unit,error=error)
  
END FUNCTION cp_unit_from_cp2k
!***************************************************************************

!!****f* cp_units/cp_unit_transform *
!!
!!   NAME
!!     cp_unit_transform
!!
!!   FUNCTION
!!     converts from a unit to another
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - value: the value to convert
!!     - unit_from: the unit of the value
!!     - unit_to: the unit of the result
!!     - power: the power of the unit (defaults to 1)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_unit_transform(value,unit_from,unit_to,defaults,power,error) &
     RESULT(res)
  REAL(kind=dp), intent(in) :: value
  TYPE(cp_unit_type), pointer :: unit_from, unit_to
  TYPE(cp_unit_set_type), POINTER, optional :: defaults
  INTEGER, INTENT(in), optional :: power
  type(cp_error_type), optional, intent(inout) :: error
  real(kind=dp) :: res

  logical :: failure
  integer :: my_power
  character(len=*), parameter :: routineN='cp_unit_transform',&
       routineP=moduleN//':'//routineN

  failure=.false.
  my_power=1
  IF (PRESENT(power)) my_power=power
  call cp_assert(cp_unit_compatible(unit_from,unit_to,error=error),&
       cp_failure_level,cp_assertion_failed,routineP,"incompatible units"//&
       TRIM(cp_unit_desc(unit_from,error=error))//", "//&
       cp_unit_desc(unit_to,error=error),error,failure)
  IF (.NOT.failure) THEN
     IF (PRESENT(defaults)) THEN
        res=cp_unit_to_cp2k1(value,unit_from,power=my_power,defaults=defaults,&
             error=error)
     ELSE
        res=cp_unit_from_cp2k1(res,unit_to,power=my_power,error=error)
     END IF
  END IF
END FUNCTION cp_unit_transform
!***************************************************************************

!!****f* cp_units/cp_unit_compatible *
!!
!!   NAME
!!     cp_unit_compatible
!!
!!   FUNCTION
!!     returs true if the two units are compatible
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_unit_compatible(unit1,unit2,error) RESULT(res)
  logical ::res
  TYPE(cp_unit_type), POINTER :: unit1, unit2
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_unit_compatible',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  !FM to do (unit=none, compress, substitute reducing to basic
  !FM        SI kinds of unit, compress, compare)
  res=.TRUE.
END FUNCTION cp_unit_compatible
!***************************************************************************

!!****f* cp_units/cp_unit_set_create *
!!
!!   NAME
!!     cp_unit_set_init
!!
!!   FUNCTION
!!     initializes the given unit set
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - unit_set: the set to initialize
!!     - name: the name of the set, used for the dafault initialization of
!!       the various units
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_unit_set_create(unit_set,name,error)
  TYPE(cp_unit_set_type), pointer :: unit_set
  CHARACTER(len=*), intent(in) :: name
  type(cp_error_type), optional, intent(inout) :: error

  logical :: failure
  character(len=*), parameter :: routineN='cp_unit_set_init',&
       routineP=moduleN//':'//routineN
  character(len=default_string_length) :: my_name
  INTEGER :: stat, i

  failure=.false.
  CPPrecondition(.NOT.ASSOCIATED(unit_set),cp_failure_level,routineP,error,failure)
  ALLOCATE(unit_set,stat=stat)
  CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
  unit_set%ref_count=1
  last_unit_set_id=last_unit_set_id+1
  unit_set%id_nr=last_unit_set_id
  my_name=name
  call uppercase(my_name)
  IF (.NOT.failure) THEN

     DO i=1,cp_ukind_max
        NULLIFY(unit_set%units(i)%unit)
     END DO
     DO i=1,cp_ukind_max
        SELECT CASE(name)
        CASE('ATOM','ATOMIC','INTERNAL','CP2K')
           IF (i==cp_ukind_angle) THEN
              CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/),&
                   unit_id=(/cp_units_rad/), power=(/1/), error=error)
           ELSE
              CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/),&
                   unit_id=(/cp_units_au/), power=(/1/), error=error)
           END IF
        CASE('OUTPUT')
           SELECT CASE(i)
           CASE(cp_ukind_energy)
              CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/), unit_id=(/cp_units_hartree/),&
                   power=(/1/), error=error)
           CASE (cp_ukind_length)
              CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/), unit_id=(/cp_units_angstrom/),&
                   power=(/1/), error=error)
           CASE (cp_ukind_temperature)
              CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/), unit_id=(/cp_units_k/),&
                   power=(/1/), error=error)
           CASE (cp_ukind_angle)
              CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/), unit_id=(/cp_units_deg/),&
                   power=(/1/), error=error)
           CASE (cp_ukind_pressure)
              CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/), unit_id=(/cp_units_bar/),&
                   power=(/1/), error=error)
           CASE (cp_ukind_time)
              CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/), unit_id=(/cp_units_fs/),&
                   power=(/1/), error=error)
           CASE (cp_ukind_mass)
              CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/), unit_id=(/cp_units_amu/),&
                   power=(/1/), error=error)
           CASE default
              CALL cp_assert(.FALSE.,cp_assertion_failed,cp_failure_level,routineP,&
                   "unhandled unit type "//TRIM(cp_to_string(i)),error,failure)
              EXIT
           END SELECT
        CASE default
           CALL cp_assert(.FALSE.,cp_assertion_failed,cp_failure_level,&
                routineP,'unknown parameter set name '//TRIM(name),&
                error,failure)
        END SELECT
     END DO
  END IF
END SUBROUTINE cp_unit_set_create
!***************************************************************************

!!****f* cp_units/cp_unit_set_retain *
!!
!!   NAME
!!     cp_unit_set_retain
!!
!!   FUNCTION
!!     retains the given unit set
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - unit_set: the unit set to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_unit_set_retain(unit_set,error)
  TYPE(cp_unit_set_type), pointer :: unit_set
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_unit_set_retain',&
        routineP=moduleN//':'//routineN
  failure=.false.

  CPPrecondition(ASSOCIATED(unit_set),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     CPPrecondition(unit_set%ref_count>0,cp_failure_level,routineP,error,failure)
     unit_set%ref_count=unit_set%ref_count+1
  end if
END SUBROUTINE cp_unit_set_retain
!***************************************************************************

!!****f* cp_units/cp_unit_set_release *
!!
!!   NAME
!!     cp_unit_set_release
!!
!!   FUNCTION
!!     releases the given unit set
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - unit_set: the unit set to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_unit_set_release(unit_set,error)
  TYPE(cp_unit_set_type), pointer :: unit_set
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  INTEGER :: stat, i
  character(len=*), parameter :: routineN='cp_unit_set_release',&
        routineP=moduleN//':'//routineN
  failure=.false.

  IF (ASSOCIATED(unit_set)) THEN
     CPPrecondition(unit_set%ref_count>0,cp_failure_level,routineP,error,failure)
     unit_set%ref_count=unit_set%ref_count-1
     IF (unit_set%ref_count == 0) THEN
        DO i = 1, SIZE(unit_set%units)
           CALL cp_unit_release(unit_set%units(i)%unit,error)
        END DO
        DEALLOCATE(unit_set,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE cp_unit_set_release
!***************************************************************************

END MODULE cp_units
