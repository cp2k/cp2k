!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief unit conversion facility
!>
!>      Units are complex, this module does not try to be very smart, for
!>      example SI prefixes are not supported automatically, and
!>      which kinds are really basic can change depending on the system of
!>      units choosen, and equivalences are not always catched.
!>
!>      This is thought as a simple conversion facility for the input and output.
!>      If you need something more you are probably better off using the
!>      physcon module directly.
!> \note
!>      One design choice was not to use dynamically allocated elements to
!>      reduce the possibility of leaks.
!>      Needs to be extended (for example charge, dipole,...)
!>      I just added the units and kinds that I needed.
!>      Used by the parser
!>      Should keep an unsorted/uncompressed version for nicer labels?
!> \par History
!>      01.2005 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE cp_units
  
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: radians,&
                                             twopi
  USE physcon,                         ONLY: &
       atm, bar, bohr, e_mass, evolt, femtoseconds, joule, kcalmol, kelvin, &
       kjmol, massunit, newton, pascal, picoseconds, seconds, wavenumbers
  USE string_utilities,                ONLY: compress,&
                                             s2a,&
                                             uppercase
#include "cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER          :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_units'
  INTEGER, SAVE, PRIVATE               :: last_unit_id=0, last_unit_set_id=0

  INTEGER, PARAMETER, PUBLIC :: cp_ukind_none=0,&
                                cp_ukind_energy=1,&
                                cp_ukind_length=2,&
                                cp_ukind_temperature=3,&
                                cp_ukind_angle=4,&
                                cp_ukind_pressure=5,&
                                cp_ukind_time=6,&
                                cp_ukind_mass=7,&
                                cp_ukind_undef=8,&
                                cp_ukind_potential=9,&
                                cp_ukind_force=10,&
                                cp_ukind_max=10

  ! General
  INTEGER, PARAMETER, PUBLIC :: cp_units_none=100,&
                                cp_units_au=101
  ! Mass
  INTEGER, PARAMETER, PUBLIC :: cp_units_m_e=110,&
                                cp_units_amu=111,&
                                cp_units_kg=112
  ! Energy
  INTEGER, PARAMETER, PUBLIC :: cp_units_hartree=130,&
                                cp_units_wavenum=131,&
                                cp_units_joule=132,&
                                cp_units_kcalmol=133,&
                                cp_units_Ry=134,&
                                cp_units_eV=135,&
                                cp_units_kjmol=136,&
                                cp_units_jmol=137,&
                                cp_units_keV=138

  ! Length
  INTEGER, PARAMETER, PUBLIC :: cp_units_bohr=140,&
                                cp_units_angstrom=141,&
                                cp_units_m=142,&
                                cp_units_pm=143,&
                                cp_units_nm=144

  ! Temperature
  INTEGER, PARAMETER, PUBLIC :: cp_units_k=150

  ! Pressure
  INTEGER, PARAMETER, PUBLIC :: cp_units_bar=161
  INTEGER, PARAMETER, PUBLIC :: cp_units_atm=162
  INTEGER, PARAMETER, PUBLIC :: cp_units_kbar=163
  INTEGER, PARAMETER, PUBLIC :: cp_units_Pa=164
  INTEGER, PARAMETER, PUBLIC :: cp_units_MPa=165
  INTEGER, PARAMETER, PUBLIC :: cp_units_GPa=166

  ! Angles
  INTEGER, PARAMETER, PUBLIC :: cp_units_rad=170,&
                                cp_units_deg=171

  ! Time
  INTEGER, PARAMETER, PUBLIC :: cp_units_fs=180,&
                                cp_units_s=181,&
                                cp_units_wn=182,&
                                cp_units_ps=183

  ! Potential
  INTEGER, PARAMETER, PUBLIC :: cp_units_volt=190

  ! Force
  INTEGER, PARAMETER, PUBLIC :: cp_units_Newton=200,&
                                cp_units_mNewton=201

  INTEGER, PARAMETER, PUBLIC :: cp_unit_max_kinds=8, cp_unit_basic_desc_length=15,&
       cp_unit_desc_length=cp_unit_max_kinds*cp_unit_basic_desc_length

  PUBLIC :: cp_unit_type, cp_unit_set_type
  PUBLIC :: cp_unit_create, cp_unit_retain, cp_unit_release,&
            cp_unit_to_cp2k, cp_unit_from_cp2k, cp_unit_desc,&
            cp_unit_set_create, cp_unit_set_release,&
            cp_unit_to_cp2k1, cp_unit_from_cp2k1, cp_unit_compatible, print_all_units

! *****************************************************************************
!> \brief stores a unit
!> \param ref_count the reference count (see doc/ReferenceCounting.html)
!> \param id_nr identification number (unique for each instance)
!> \param kind the kind of unit (energy, length,...)
!> \param unit the actual unit (Joule, eV,...)
!> \author fawzi
! *****************************************************************************
  TYPE cp_unit_type
     INTEGER :: id_nr,ref_count,n_kinds
     INTEGER, DIMENSION(cp_unit_max_kinds):: kind_id, unit_id, power
  END TYPE cp_unit_type

! *****************************************************************************
!> \brief represent a pointer to a unit (to build arrays of pointers)
!> \param unit the pointer to the unit
!> \author fawzi
! *****************************************************************************
  TYPE cp_unit_p_type
     TYPE(cp_unit_type), POINTER :: unit
  END TYPE cp_unit_p_type

! *****************************************************************************
!> \brief stores the default units to be used
!> \author fawzi
! *****************************************************************************
  TYPE cp_unit_set_type
     INTEGER :: id_nr, ref_count
     TYPE(cp_unit_p_type), DIMENSION(cp_ukind_max) :: units
  END TYPE cp_unit_set_type

CONTAINS

! *****************************************************************************
!> \brief creates a unit parsing a string
!> \param unit the unit to initialize
!> \param string the string containing the description of the unit
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE cp_unit_create(unit, string, error)
    TYPE(cp_unit_type), POINTER              :: unit
    CHARACTER(len=*), INTENT(in)             :: string
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_string_length)         :: desc
    CHARACTER(LEN=40)                        :: formatstr
    INTEGER                                  :: i_high, i_low, i_unit, &
                                                len_string, next_power
    INTEGER, DIMENSION(cp_unit_max_kinds)    :: kind_id, power, unit_id
    LOGICAL                                  :: failure

    failure=.FALSE.
    unit_id=cp_units_none
    kind_id=cp_ukind_none
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
               "Maximum number of combined units exceeded",&
               error,failure)
          EXIT
       END IF
       ! read unit
       SELECT CASE(string(i_low:i_high-1))
       CASE("internal_cp2k")
          unit_id(i_unit)=cp_units_none
          kind_id(i_unit)=cp_ukind_undef
       CASE("hartree")
          unit_id(i_unit)=cp_units_hartree
          kind_id(i_unit)=cp_ukind_energy
       CASE("au_e")
          unit_id(i_unit)=cp_units_au
          kind_id(i_unit)=cp_ukind_energy
       CASE("wavenumber_e")
          unit_id(i_unit)=cp_units_wavenum
          kind_id(i_unit)=cp_ukind_energy
       CASE("joule","J")
          unit_id(i_unit)=cp_units_joule
          kind_id(i_unit)=cp_ukind_energy
       CASE("kcalmol")
          unit_id(i_unit)=cp_units_kcalmol
          kind_id(i_unit)=cp_ukind_energy
       CASE("kjmol")
          unit_id(i_unit)=cp_units_kjmol
          kind_id(i_unit)=cp_ukind_energy
       CASE("jmol")
          unit_id(i_unit)=cp_units_jmol
          kind_id(i_unit)=cp_ukind_energy
       CASE("Ry")
          unit_id(i_unit)=cp_units_Ry
          kind_id(i_unit)=cp_ukind_energy
       CASE("eV")
          unit_id(i_unit)=cp_units_eV
          kind_id(i_unit)=cp_ukind_energy
       CASE("keV")
          unit_id(i_unit)=cp_units_keV
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
       CASE("pm")
          unit_id(i_unit)=cp_units_pm
          kind_id(i_unit)=cp_ukind_length
       CASE("nm")
          unit_id(i_unit)=cp_units_nm
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
       CASE("atm")
          unit_id(i_unit)=cp_units_atm
          kind_id(i_unit)=cp_ukind_pressure
       CASE("bar")
          unit_id(i_unit)=cp_units_bar
          kind_id(i_unit)=cp_ukind_pressure
       CASE("kbar")
          unit_id(i_unit)=cp_units_kbar
          kind_id(i_unit)=cp_ukind_pressure
       CASE("Pa")
          unit_id(i_unit)=cp_units_Pa
          kind_id(i_unit)=cp_ukind_pressure
       CASE("MPa")
          unit_id(i_unit)=cp_units_MPa
          kind_id(i_unit)=cp_ukind_pressure
       CASE("GPa")
          unit_id(i_unit)=cp_units_GPa
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
       CASE("ps")
          unit_id(i_unit)=cp_units_ps
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
       CASE("volt")
          unit_id(i_unit)=cp_units_volt
          kind_id(i_unit)=cp_ukind_potential
       CASE("au_pot")
          unit_id(i_unit)=cp_units_au
          kind_id(i_unit)=cp_ukind_potential
       CASE("potential")
          unit_id(i_unit)=cp_units_none
          kind_id(i_unit)=cp_ukind_potential
       CASE("N","Newton")
          unit_id(i_unit)=cp_units_Newton
          kind_id(i_unit)=cp_ukind_force
       CASE("mN","mNewton")
          unit_id(i_unit)=cp_units_mNewton
          kind_id(i_unit)=cp_ukind_force
       CASE("au_f")
          unit_id(i_unit)=cp_units_au
          kind_id(i_unit)=cp_ukind_force
       CASE("force")
          unit_id(i_unit)=cp_units_none
          kind_id(i_unit)=cp_ukind_force
       CASE("au")
          CALL cp_assert(.FALSE.,cp_failure_level,&
               cp_assertion_failed,routineP,&
               "au unit without specifing its kind not accepted, use "//&
               "(au_e, au_f, au_t, au_temp, au_l, au_m, au_p, au_pot)",&
               error,failure)
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"Unknown unit: "//string(i_low:i_high-1),&
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

       IF (i_high<=len_string) THEN
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
       ENDIF
       IF (i_low>len_string) EXIT
       next_power=1
       IF (i_high<=len_string) THEN
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
       ENDIF
    END DO
    CALL cp_unit_create2(unit,kind_id=kind_id, unit_id=unit_id, &
         power=power, error=error)
    desc=cp_unit_desc(unit,error=error)
  END SUBROUTINE cp_unit_create

! *****************************************************************************
!> \brief creates and initializes the given unit of mesure (performs some error
!>      check)
!> \param unit the unit descriptor to be initialized
!> \param kind_id the kind of unit (length,energy,...), use the constants
!>        cp_ukind_*
!> \param unit_id the actual unit (use constants cp_units_*)
!> \param power ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE cp_unit_create2(unit, kind_id, unit_id, power, error)
    TYPE(cp_unit_type), POINTER              :: unit
    INTEGER, DIMENSION(:), INTENT(in)        :: kind_id, unit_id
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: power
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_create2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, max_kind, max_pos, stat
    LOGICAL                                  :: failure, repeat

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(unit),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(kind_id)<=cp_unit_max_kinds,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(unit_id)<=cp_unit_max_kinds,cp_failure_level,routineP,error,failure)
    ALLOCATE(unit,stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    unit%ref_count=1
    last_unit_id=last_unit_id+1
    unit%id_nr=last_unit_id
    unit%kind_id(1:SIZE(kind_id))=kind_id
    unit%kind_id(SIZE(kind_id)+1:)=cp_ukind_none
    unit%unit_id(1:SIZE(unit_id))=unit_id
    unit%unit_id(SIZE(unit_id):)=cp_units_none
    IF (PRESENT(power)) THEN
       unit%power(1:SIZE(power))=power
       unit%power(SIZE(power)+1:)=0
       DO i=1,SIZE(unit%power)
          IF (unit%power(i)==0) THEN
             unit%kind_id(i)=cp_ukind_none
             unit%unit_id(i)=cp_units_none
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
          repeat=.FALSE.
          DO j=i+1,SIZE(unit%kind_id)
             IF (unit%kind_id(j)>=max_kind) THEN
                IF (unit%kind_id(j)/=0.AND.unit%kind_id(j)==max_kind.AND.&
                     unit%unit_id(j)==unit%unit_id(max_pos)) THEN
                   unit%power(max_pos)=unit%power(max_pos)+unit%power(j)
                   unit%kind_id(j)=cp_ukind_none
                   unit%unit_id(j)=cp_units_none
                   unit%power(j)=0
                   IF (unit%power(max_pos)==0) THEN
                      unit%kind_id(max_pos)=cp_ukind_none
                      unit%unit_id(max_pos)=cp_units_none
                      unit%power(max_pos)=0
                      repeat=.TRUE.
                      EXIT
                   END IF
                ELSE IF (unit%kind_id(j)>max_kind.OR.&
                     (unit%kind_id(j)==max_kind.AND.&
                     unit%unit_id(j)>unit%unit_id(max_pos))) THEN
                   max_kind=unit%kind_id(j)
                   max_pos=j
                END IF
             END IF
          END DO
          IF (.not.repeat) EXIT
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

! *****************************************************************************
!> \brief retains the given unit
!> \param unit the unit to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      at the moment not needed, there for completeness
! *****************************************************************************
  SUBROUTINE cp_unit_retain(unit,error)
    TYPE(cp_unit_type), POINTER              :: unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(unit),cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(unit%ref_count>0,cp_failure_level,routineP,error)
    unit%ref_count=unit%ref_count+1
  END SUBROUTINE cp_unit_retain

! *****************************************************************************
!> \brief releases the given unit
!> \param unit the unit to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      at the moment not needed, there for completeness
! *****************************************************************************
  SUBROUTINE cp_unit_release(unit,error)
    TYPE(cp_unit_type), POINTER              :: unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(unit)) THEN
       CPPreconditionNoFail(unit%ref_count>0,cp_failure_level,routineP,error)
       unit%ref_count=unit%ref_count-1
       IF (unit%ref_count==0) THEN
          DEALLOCATE(unit,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
  END SUBROUTINE cp_unit_release

! *****************************************************************************
!> \brief controls that the kind and contains meaningful information
!> \param basic_kind the kind of the unit
!> \param basic_unit the unit to check
!> \param error_level error level of the errors due to invalid values.
!>        Defaults to cp_failure_level
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_basic_unit_check(basic_kind,basic_unit,error_level,error)&
       RESULT(res)
    INTEGER, INTENT(in)                      :: basic_kind, basic_unit
    INTEGER, INTENT(in), OPTIONAL            :: error_level
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_basic_unit_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_error_level
    LOGICAL                                  :: failure

    failure=.FALSE.
    my_error_level=cp_failure_level
    IF (PRESENT(error_level)) my_error_level=error_level

    SELECT CASE(basic_kind)
    CASE(cp_ukind_undef)
       SELECT CASE (basic_unit)
       CASE(cp_units_none)
       CASE default
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown undef unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_energy)
       SELECT CASE (basic_unit)
       CASE(cp_units_hartree, cp_units_wavenum, cp_units_joule, cp_units_kcalmol,&
            cp_units_kjmol, cp_units_Ry, cp_units_eV, cp_units_keV, cp_units_au, cp_units_k,&
            cp_units_jmol, cp_units_none)
       CASE default
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown energy unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_length)
       SELECT CASE (basic_unit)
       CASE(cp_units_bohr, cp_units_angstrom, cp_units_au, cp_units_none, cp_units_m,&
            cp_units_pm, cp_units_nm)
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
       CASE(cp_units_bar,cp_units_atm,cp_units_kbar,cp_units_Pa,cp_units_MPa,cp_units_GPa,cp_units_au,cp_units_none)
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
       CASE(cp_units_s, cp_units_fs, cp_units_ps, cp_units_au, cp_units_wn, cp_units_none)
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
    CASE(cp_ukind_potential)
       SELECT CASE (basic_unit)
       CASE(cp_units_volt, cp_units_au, cp_units_none)
       CASE default
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown potential unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_force)
       SELECT CASE (basic_unit)
       CASE(cp_units_Newton, cp_units_mNewton, cp_units_au, cp_units_none)
       CASE default
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown force unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_none)
       IF(basic_unit/=cp_units_none)&
            CALL cp_assert(.FALSE.,my_error_level,&
            cp_assertion_failed, routineP,&
            "if the kind of the unit is none also unit must be undefined,not:"&
            //TRIM(cp_to_string(basic_unit)),&
            error,failure)
    CASE default
       CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
            routineP,"unknown kind of unit:"//TRIM(cp_to_string(basic_kind)),&
            error,failure)
    END SELECT
    res=.not.failure
  END FUNCTION cp_basic_unit_check

! *****************************************************************************
!> \brief converts a value to the internal cp2k units
!> \param value the value to convert
!> \param basic_kind the kind of the unit of the value
!> \param basic_unit the unit of the value
!> \param power the power of the unit (defaults to 1)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_basic_unit_to_cp2k(value,basic_kind,basic_unit,power,error)  RESULT(res)
    REAL(kind=dp), INTENT(in)                :: value
    INTEGER, INTENT(in)                      :: basic_kind, basic_unit
    INTEGER, INTENT(in), OPTIONAL            :: power
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(kind=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_basic_unit_to_cp2k', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER :: my_error_level = cp_failure_level

    INTEGER                                  :: my_power
    LOGICAL                                  :: failure

    failure=.FALSE.
    my_power=1
    IF (PRESENT(power)) my_power=power
    IF (basic_unit==cp_units_none.AND.basic_kind/=cp_ukind_undef) THEN
       CALL cp_assert(basic_kind==cp_units_none,my_error_level,&
            cp_assertion_failed,routineP,&
            "unit not yet fully specified, unit of kind "//&
            TRIM(cp_to_string(basic_unit)),error,failure)
    END IF
    SELECT CASE(basic_kind)
    CASE(cp_ukind_undef)
       SELECT CASE (basic_unit)
       CASE(cp_units_none)
          res = value
       CASE default
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown energy unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_energy)
       SELECT CASE (basic_unit)
       CASE(cp_units_hartree,cp_units_au)
          res=value
       CASE(cp_units_wavenum)
          res=wavenumbers**(-my_power)*value
       CASE(cp_units_joule)
          res=joule**(-my_power)*value
       CASE(cp_units_kcalmol)
          res=kcalmol**(-my_power)*value
       CASE(cp_units_kjmol)
          res=kjmol**(-my_power)*value
       CASE(cp_units_jmol)
          res=(kjmol*1.0E+3_dp)**(-my_power)*value
       CASE(cp_units_Ry)
          res=0.5_dp**my_power*value
       CASE(cp_units_eV)
          res=evolt**(-my_power)*value
       CASE(cp_units_keV)
          res=(1.0E-3_dp*evolt)**(-my_power)*value
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
          res=value*(1.0E10_dp*bohr)**my_power
       CASE(cp_units_pm)
          res=value*(0.01_dp*bohr)**my_power
       CASE(cp_units_nm)
          res=value*(10.0_dp*bohr)**my_power
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
          res = bar**(-my_power)*value
       CASE(cp_units_atm)
          res = atm**(-my_power)*value
       CASE(cp_units_kbar)
          res = (1.0E-3_dp*bar)**(-my_power)*value
       CASE(cp_units_Pa)
          res = pascal**(-my_power)*value
       CASE(cp_units_MPa)
          res = (1.0E-6_dp*pascal)**(-my_power)*value
       CASE(cp_units_GPa)
          res = (1.0E-9_dp*pascal)**(-my_power)*value
       CASE(cp_units_au)
          res = value
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
       CASE(cp_units_ps)
          res=value*picoseconds**(-my_power)
       CASE(cp_units_au)
          res=value
       CASE(cp_units_wn)
          res=(twopi*wavenumbers)**(my_power)/value
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
    CASE(cp_ukind_potential)
       SELECT CASE (basic_unit)
       CASE(cp_units_volt)
          res=evolt**(-my_power)*value
       CASE(cp_units_au)
          res=value
       CASE default
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown potential unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_force)
       SELECT CASE (basic_unit)
       CASE(cp_units_Newton)
          res=value*newton**(-my_power)
       CASE(cp_units_mNewton)
          res=value*(1.0E+3*newton)**(-my_power)
       CASE(cp_units_au)
          res=value
       CASE default
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown force unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_none)
       CALL cp_assert(.FALSE.,my_error_level,&
            cp_assertion_failed, routineP,&
            "if the kind of the unit is none also unit must be undefined,not:"&
            //TRIM(cp_to_string(basic_unit)),&
            error,failure)
    CASE default
       CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
            routineP,"unknown kind of unit:"//TRIM(cp_to_string(basic_kind)),&
            error,failure)
    END SELECT
  END FUNCTION cp_basic_unit_to_cp2k

! *****************************************************************************
!> \brief returns the label of the current basic unit
!> \param basic_kind the kind of the unit of the value
!> \param basic_unit the unit of the value
!> \param power the power of the unit (defaults to 1)
!> \param accept_undefined ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_basic_unit_desc(basic_kind,basic_unit,power,accept_undefined,error)&
       RESULT(res)
    INTEGER, INTENT(in)                      :: basic_kind, basic_unit
    INTEGER, INTENT(in), OPTIONAL            :: power
    LOGICAL, INTENT(in), OPTIONAL            :: accept_undefined
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(len=cp_unit_basic_desc_length) :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_basic_unit_desc', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER :: my_error_level = cp_failure_level

    INTEGER                                  :: a, my_power
    LOGICAL                                  :: failure, my_accept_undefined

    failure=.FALSE.
    my_power=1
    res=""
    my_accept_undefined=.FALSE.
    IF (accept_undefined) my_accept_undefined=accept_undefined
    IF (PRESENT(power)) my_power=power
    IF (basic_unit==cp_units_none) THEN
       CALL cp_assert(my_accept_undefined.OR.basic_kind/=cp_units_none,&
            my_error_level,cp_assertion_failed,routineP,&
            "unit not yet fully specified, unit of kind "//&
            TRIM(cp_to_string(basic_kind)),error,failure)
    END IF
    SELECT CASE(basic_kind)
    CASE(cp_ukind_undef)
       SELECT CASE (basic_unit)
       CASE(cp_units_none)
          res="internal_cp2k"
       CASE DEFAULT
          CALL cp_assert(.FALSE.,my_error_level,&
               cp_assertion_failed,routineP,&
               "unit not yet fully specified, unit of kind "//&
               TRIM(res),error,failure)
       END SELECT
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
       CASE(cp_units_kjmol)
          res="kjmol"
       CASE(cp_units_jmol)
          res="jmol"
       CASE(cp_units_Ry)
          res="Ry"
       CASE(cp_units_eV)
          res="eV"
       CASE(cp_units_keV)
          res="keV"
       CASE(cp_units_k)
          res="K_e"
       CASE(cp_units_none)
          res="energy"
          CALL cp_assert(my_accept_undefined,my_error_level,&
               cp_assertion_failed,routineP,&
               "unit not yet fully specified, unit of kind "//&
               TRIM(res),error,failure)
       CASE default
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown energy unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_length)
       SELECT CASE (basic_unit)
       CASE(cp_units_bohr,cp_units_au)
          res="bohr"
       CASE(cp_units_m)
          res="m"
       CASE(cp_units_pm)
          res="pm"
       CASE(cp_units_nm)
          res="nm"
       CASE(cp_units_angstrom)
          res="angstrom"
       CASE default
          res="length"
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
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
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown temperature unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_pressure)
       SELECT CASE (basic_unit)
       CASE(cp_units_bar)
          res="bar"
       CASE(cp_units_atm)
          res="atm"
       CASE(cp_units_kbar)
          res="kbar"
       CASE(cp_units_Pa)
          res="Pa"
       CASE(cp_units_MPa)
          res="MPa"
       CASE(cp_units_GPa)
          res="GPa"
       CASE(cp_units_au)
          res="au_p"
       CASE(cp_units_none)
          res="pressure"
          CALL cp_assert(my_accept_undefined,my_error_level,&
               cp_assertion_failed,routineP,&
               "unit not yet fully specified, unit of kind "//&
               TRIM(res),error,failure)
       CASE default
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
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
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown angle unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_time)
       SELECT CASE (basic_unit)
       CASE(cp_units_s)
          res="s"
       CASE(cp_units_fs)
          res="fs"
       CASE(cp_units_ps)
          res="ps"
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
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
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
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown mass unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_potential)
       SELECT CASE (basic_unit)
       CASE(cp_units_volt)
          res="volt"
       CASE(cp_units_au)
          res="au_pot"
       CASE(cp_units_none)
          res="potential"
          CALL cp_assert(my_accept_undefined,my_error_level,&
               cp_assertion_failed,routineP,&
               "unit not yet fully specified, unit of kind "//&
               TRIM(res),error,failure)
       CASE default
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown potential unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_force)
       SELECT CASE (basic_unit)
       CASE(cp_units_Newton)
          res="N"
       CASE(cp_units_mNewton)
          res="mN"
       CASE(cp_units_au)
          res="au_f"
       CASE(cp_units_none)
          res="force"
          CALL cp_assert(my_accept_undefined,my_error_level,&
               cp_assertion_failed,routineP,&
               "unit not yet fully specified, unit of kind "//&
               TRIM(res),error,failure)
       CASE default
          CALL cp_assert(.FALSE.,my_error_level,cp_assertion_failed,&
               routineP,"unknown potential unit:"//TRIM(cp_to_string(basic_unit)),&
               error,failure)
       END SELECT
    CASE(cp_ukind_none)
       CALL cp_assert(.FALSE.,my_error_level,&
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
       CALL compress(res,.TRUE.)
    END IF
  END FUNCTION cp_basic_unit_desc

! *****************************************************************************
!> \brief returns the "name" of the given unit
!> \param unit the unit to describe
!> \param defaults defaults for the undefined units, optional
!> \param accept_undefined if defaults is not present or is not associated
!>        whether undefined units should be accepted (defaults to false)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_unit_desc(unit,defaults,accept_undefined,error)&
       RESULT(res)
    TYPE(cp_unit_type), POINTER              :: unit
    TYPE(cp_unit_set_type), OPTIONAL, &
      POINTER                                :: defaults
    LOGICAL, INTENT(in), OPTIONAL            :: accept_undefined
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(len=cp_unit_desc_length)       :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_desc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, my_unit, pos
    LOGICAL                                  :: check, failure, has_defaults, &
                                                my_accept_undefined

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(unit),cp_failure_level,routineP,error,failure)
    res=""
    pos=1
    my_accept_undefined=.FALSE.
    IF (PRESENT(accept_undefined)) my_accept_undefined=accept_undefined
    DO i=1,unit%n_kinds
       CPPrecondition(unit%kind_id(i)/=0,cp_failure_level,routineP,error,failure)
       CPPrecondition(pos<LEN(res),cp_failure_level,routineP,error,failure)
       IF (failure) EXIT
       my_unit=unit%unit_id(i)
       has_defaults=.FALSE.
       IF (PRESENT(defaults)) has_defaults=ASSOCIATED(defaults)
       IF (my_unit==0) THEN
          IF (has_defaults) THEN
             my_unit=defaults%units(unit%kind_id(i))%unit%unit_id(1)
          ELSE
             check = my_accept_undefined.OR.unit%kind_id(i)/=0
             CPPrecondition(check,cp_failure_level,routineP,error,failure)
          END IF
       END IF
       IF (i > 1) THEN
          res(pos:pos) = "*"
          pos = pos + 1
       END IF
       res(pos:)=TRIM(cp_basic_unit_desc(basic_kind=unit%kind_id(i),&
            basic_unit=my_unit,accept_undefined=my_accept_undefined,&
            power=unit%power(i),error=error))
       pos=LEN_TRIM(res)+1
    END DO

  END FUNCTION cp_unit_desc

! *****************************************************************************
!> \brief tranform a value to the internal cp2k units
!> \param value the value to convert
!> \param unit the unit of the result
!> \param defaults the defaults unit for those that are left free
!>        (cp_units_none)
!> \param power the power of the unit (defaults to 1)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_unit_to_cp2k1(value,unit,defaults,power,error) RESULT(res)
    REAL(kind=dp), INTENT(in)                :: value
    TYPE(cp_unit_type), POINTER              :: unit
    TYPE(cp_unit_set_type), OPTIONAL, &
      POINTER                                :: defaults
    INTEGER, INTENT(in), OPTIONAL            :: power
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(kind=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_to_cp2k1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_unit, my_basic_unit, &
                                                my_power
    LOGICAL                                  :: failure

    failure=.FALSE.
    my_power=1
    IF (PRESENT(power)) my_power=power
    res=value
    CPPrecondition(ASSOCIATED(unit),cp_failure_level,routineP,error,failure)
    DO i_unit=1,unit%n_kinds
       CPPrecondition(unit%kind_id(i_unit)>0,cp_failure_level,routineP,error,failure)
       my_basic_unit=unit%unit_id(i_unit)
       IF (my_basic_unit==0.AND.unit%kind_id(i_unit)/=cp_ukind_undef) THEN
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
  END FUNCTION cp_unit_to_cp2k1

! *****************************************************************************
!> \brief converts from the internal cp2k units to the given unit
!> \param value the value to convert
!> \param unit the unit of the result
!> \param defaults the defaults unit for those that are left free
!>        (cp_units_none)
!> \param power the power of the unit (defaults to 1)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_unit_from_cp2k1(value,unit,defaults,power,error) RESULT(res)
    REAL(kind=dp), INTENT(in)                :: value
    TYPE(cp_unit_type), POINTER              :: unit
    TYPE(cp_unit_set_type), OPTIONAL, &
      POINTER                                :: defaults
    INTEGER, INTENT(in), OPTIONAL            :: power
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(kind=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_from_cp2k1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_power
    LOGICAL                                  :: failure

    failure=.FALSE.
    my_power=1
    IF (PRESENT(power)) my_power=power
    IF (PRESENT(defaults)) THEN
       res=cp_unit_to_cp2k1(value=value,unit=unit,defaults=defaults,&
            power=-my_power,error=error)
    ELSE
       res=cp_unit_to_cp2k1(value=value,unit=unit,power=-my_power,error=error)
    END IF
  END FUNCTION cp_unit_from_cp2k1

! *****************************************************************************
!> \brief converts to the internal cp2k units to the given unit
!> \param value the value to convert
!> \param unit_str the unit of the result as string
!> \param defaults the defaults unit for those that are left free
!>        (cp_units_none)
!> \param power the power of the unit (defaults to 1)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_unit_to_cp2k(value,unit_str,defaults,power,error) RESULT(res)
    REAL(kind=dp), INTENT(in)                :: value
    CHARACTER(len=*), INTENT(in)             :: unit_str
    TYPE(cp_unit_set_type), OPTIONAL, &
      POINTER                                :: defaults
    INTEGER, INTENT(in), OPTIONAL            :: power
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(kind=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_to_cp2k', &
      routineP = moduleN//':'//routineN

    TYPE(cp_unit_type), POINTER              :: my_unit

    NULLIFY(my_unit)
    CALL cp_unit_create(my_unit,unit_str,error=error)
    IF (PRESENT(defaults)) THEN
       res=cp_unit_to_cp2k1(value=value,unit=my_unit,defaults=defaults,&
            power=power,error=error)
    ELSE
       res=cp_unit_to_cp2k1(value=value,unit=my_unit,power=power,error=error)
    END IF
    CALL cp_unit_release(my_unit,error=error)
  END FUNCTION cp_unit_to_cp2k

! *****************************************************************************
!> \brief converts from the internal cp2k units to the given unit
!> \param value the value to convert
!> \param unit_str the unit of the result as string
!> \param defaults the defaults unit for those that are left free
!>        (cp_units_none)
!> \param power the power of the unit (defaults to 1)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_unit_from_cp2k(value,unit_str,defaults,power,error) RESULT(res)
    REAL(kind=dp), INTENT(in)                :: value
    CHARACTER(len=*), INTENT(in)             :: unit_str
    TYPE(cp_unit_set_type), OPTIONAL, &
      POINTER                                :: defaults
    INTEGER, INTENT(in), OPTIONAL            :: power
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(kind=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_from_cp2k', &
      routineP = moduleN//':'//routineN

    TYPE(cp_unit_type), POINTER              :: my_unit

    NULLIFY(my_unit)
    CALL cp_unit_create(my_unit,unit_str,error=error)
    IF (PRESENT(defaults)) THEN
       res=cp_unit_from_cp2k1(value=value,unit=my_unit,defaults=defaults,&
            power=power,error=error)
    ELSE
       res=cp_unit_from_cp2k1(value=value,unit=my_unit,power=power,error=error)
    END IF
    CALL cp_unit_release(my_unit,error=error)
  END FUNCTION cp_unit_from_cp2k

! *****************************************************************************
!> \brief returs true if the two units are compatible
!> \param ref_unit ...
!> \param unit ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author Teodoro Laino [tlaino] - 11.2007 - University of Zurich
! *****************************************************************************
  FUNCTION cp_unit_compatible(ref_unit,unit,error) RESULT(res)
    TYPE(cp_unit_type), POINTER              :: ref_unit, unit
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_compatible', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    res = .TRUE.
    DO i = 1, SIZE(ref_unit%kind_id)
       IF  (ref_unit%kind_id(i)==unit%kind_id(i)) CYCLE
       IF ((ref_unit%kind_id(1)==cp_ukind_undef).AND.(ALL(ref_unit%kind_id(2:)==cp_ukind_none))) CYCLE
       res = .FALSE.
       EXIT
    END DO

  END FUNCTION cp_unit_compatible

! *****************************************************************************
!> \brief initializes the given unit set
!> \param unit_set the set to initialize
!> \param name the name of the set, used for the dafault initialization of
!>        the various units
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE cp_unit_set_create(unit_set,name,error)
    TYPE(cp_unit_set_type), POINTER          :: unit_set
    CHARACTER(len=*), INTENT(in)             :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_set_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: my_name
    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(unit_set),cp_failure_level,routineP,error,failure)
    ALLOCATE(unit_set,stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    unit_set%ref_count=1
    last_unit_set_id=last_unit_set_id+1
    unit_set%id_nr=last_unit_set_id
    my_name=name
    CALL uppercase(my_name)

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
          CASE(cp_ukind_undef)
             CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/), unit_id=(/cp_units_none/),&
                  power=(/1/), error=error)
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
          CASE (cp_ukind_potential)
             CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/), unit_id=(/cp_units_volt/),&
                  power=(/1/), error=error)
          CASE (cp_ukind_force)
             CALL cp_unit_create2(unit_set%units(i)%unit, kind_id=(/i/), unit_id=(/cp_units_newton/),&
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
  END SUBROUTINE cp_unit_set_create

! *****************************************************************************
!> \brief retains the given unit set
!> \param unit_set the unit set to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE cp_unit_set_retain(unit_set,error)
    TYPE(cp_unit_set_type), POINTER          :: unit_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_set_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(unit_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(unit_set%ref_count>0,cp_failure_level,routineP,error,failure)
    unit_set%ref_count=unit_set%ref_count+1
  END SUBROUTINE cp_unit_set_retain

! *****************************************************************************
!> \brief releases the given unit set
!> \param unit_set the unit set to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE cp_unit_set_release(unit_set,error)
    TYPE(cp_unit_set_type), POINTER          :: unit_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unit_set_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

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

! *****************************************************************************
!> \brief Prints info on all available units in CP2K
!> \param unit_nr ...
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE print_all_units(unit_nr)
    INTEGER, INTENT(IN)                      :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'print_all_units', &
      routineP = moduleN//':'//routineN

    CALL print_section_unit(unit_label="Undefined",description="If the default unit "//&
         "of a keyword is explicitly undefined, all possible units of measurement can "//&
         "be used to define a proper value.",&
         units_set=s2a("internal_cp2k"), unit_nr=unit_nr)

    CALL print_section_unit(unit_label="Energy",description="Possible units of measurement "//&
         "for Energies. The [energy] entry acts like a dummy flag (assumes the unit of "//&
         "measurement of energy is in internal units), useful for dimensional analysis.",&
         units_set=s2a("hartree","wavenumber_e","joule","kcalmol","kjmol","Ry",&
         "eV","K_e","energy"), unit_nr=unit_nr)

    CALL print_section_unit(unit_label="Length",description="Possible units of measurement "//&
         "for Lengths. The [length] entry acts like a dummy flag (assumes the unit of "//&
         "measurement of length is in internal units), useful for dimensional analysis.",&
         units_set=s2a("bohr","m","pm","nm","angstrom","length"), unit_nr=unit_nr)

    CALL print_section_unit(unit_label="Temperature",description="Possible units of measurement "//&
         "for Temperature. The [temperature] entry acts like a dummy flag (assumes the unit of "//&
         "measurement of temperature is in internal units), useful for dimensional analysis.",&
         units_set=s2a("K","au_temp","temperature"), unit_nr=unit_nr)

    CALL print_section_unit(unit_label="Pressure",description="Possible units of measurement "//&
         "for Pressure. The [pressure] entry acts like a dummy flag (assumes the unit of "//&
         "measurement of pressure is in internal units), useful for dimensional analysis.",&
         units_set=s2a("bar","atm","kbar","Pa","MPa","GPa","au_p","pressure"),&
         unit_nr=unit_nr)

    CALL print_section_unit(unit_label="Angle",description="Possible units of measurement "//&
         "for Angles. The [angle] entry acts like a dummy flag (assumes the unit of "//&
         "measurement of angle is in internal units), useful for dimensional analysis.",&
         units_set=s2a("rad","deg","angle"),unit_nr=unit_nr)

    CALL print_section_unit(unit_label="Time",description="Possible units of measurement "//&
         "for Time. The [time] entry acts like a dummy flag (assumes the unit of "//&
         "measurement of time is in internal units), useful for dimensional analysis.",&
         units_set=s2a("s","fs","ps","au_t","wavenumber_t","time"),unit_nr=unit_nr)

    CALL print_section_unit(unit_label="Mass",description="Possible units of measurement "//&
         "for Masses. The [mass] entry acts like a dummy flag (assumes the unit of "//&
         "measurement of mass is in internal units), useful for dimensional analysis.",&
         units_set=s2a("kg","amu","m_e","mass"),unit_nr=unit_nr)

    CALL print_section_unit(unit_label="Potential",description="Possible units of measurement "//&
         "for potentials. The [potential] entry acts like a dummy flag (assumes the unit of "//&
         "measurement of potential is in internal units), useful for dimensional analysis.",&
         units_set=s2a("volt","au_pot","potential"),unit_nr=unit_nr)

    CALL print_section_unit(unit_label="Force",description="Possible units of measurement "//&
         "for forces. The [force] entry acts like a dummy flag (assumes the unit of "//&
         "measurement of force is in internal units), useful for dimensional analysis.",&
         units_set=s2a("N","Newton","mN","mNewton","au_f","force"),unit_nr=unit_nr)

  END SUBROUTINE print_all_units

! *****************************************************************************
!> \brief Prints info on all available units in CP2K - Low level
!> \param unit_label ...
!> \param description ...
!> \param units_set ...
!> \param unit_nr ...
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE print_section_unit(unit_label, description, units_set, unit_nr)
    CHARACTER(LEN=*), INTENT(IN)             :: unit_label, description
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(IN)                             :: units_set
    INTEGER, INTENT(IN)                      :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'print_section_unit', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    WRITE(unit_nr,FMT='(A)') "<H2>"//TRIM(unit_label)//"</H2>"
    WRITE(unit_nr,FMT='(A)')description//"<BR><DL>"
    DO i =1, SIZE(units_set)
       WRITE(unit_nr,FMT='(A)')"<DD><B>"//TRIM(units_set(i))//"</B></DD>"
    END DO
    WRITE(unit_nr,FMT='(A)')"</DL><P>"
  END SUBROUTINE print_section_unit

END MODULE cp_units
