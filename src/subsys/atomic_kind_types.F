!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Define the atomic kind types and their sub types
!> \author  Matthias Krack (MK)
!> \date    02.01.2002
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - 01.2002 creation [MK]
!> - 04.2002 added pao [fawzi]
!> - 09.2002 adapted for POL/KG use [GT]
!> - 02.2004 flexible normalization of basis sets [jgh]
!> - 03.2004 attach/detach routines [jgh]
!> - 10.2004 removed pao [fawzi]
!> - 08.2014 moevd qs-related stuff into new qs_kind_types.F [Ole Schuett]
! *****************************************************************************
MODULE atomic_kind_types
  USE damping_dipole_types,            ONLY: damping_p_release,&
                                             damping_p_type
  USE external_potential_types,        ONLY: deallocate_potential,&
                                             fist_potential_type,&
                                             get_potential
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE periodic_table,                  ONLY: get_ptable_info
  USE shell_potential_types,           ONLY: shell_kind_type,&
                                             shell_release,&
                                             shell_retain
  USE termination,                     ONLY: stop_program
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atomic_kind_types'


!> Provides all information about an atomic kind
! *****************************************************************************
  TYPE atomic_kind_type
     TYPE(fist_potential_type), POINTER     :: fist_potential => Null()
     CHARACTER(LEN=default_string_length)   :: name = ""
     CHARACTER(LEN=2)                       :: element_symbol = ""
     REAL(KIND = dp)                        :: mass = 0.0_dp
     INTEGER                                :: kind_number = -1
     INTEGER                                :: natom = -1
     INTEGER, DIMENSION(:), POINTER         :: atom_list => Null()
     LOGICAL                                :: shell_active = .FALSE.
     TYPE(shell_kind_type), POINTER         :: shell => Null()
     TYPE(damping_p_type), POINTER          :: damping => Null()
  END TYPE atomic_kind_type

!> Provides a vector of pointers of type atomic_kind_type
! *****************************************************************************
  TYPE atomic_kind_p_type
     TYPE(atomic_kind_type),DIMENSION(:),&
          POINTER                             :: atomic_kind_set
  END TYPE atomic_kind_p_type

  ! Public subroutines

  PUBLIC :: deallocate_atomic_kind_set,&
            get_atomic_kind,&
            get_atomic_kind_set,&
            set_atomic_kind,&
            is_hydrogen

  ! Public data types
  PUBLIC :: atomic_kind_type

CONTAINS

! *****************************************************************************
!> \brief   Destructor routine for a set of atomic kinds
!> \param atomic_kind_set ...
!> \param error ...
!> \date    02.01.2002
!> \author  Matthias Krack (MK)
!> \version 2.0
! *****************************************************************************
  SUBROUTINE deallocate_atomic_kind_set(atomic_kind_set,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_atomic_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, nkind, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(atomic_kind_set)) THEN

       nkind = SIZE(atomic_kind_set)

       DO ikind=1,nkind
          IF (ASSOCIATED(atomic_kind_set(ikind)%fist_potential)) THEN
             CALL deallocate_potential(atomic_kind_set(ikind)%fist_potential,error)
          END IF
          IF (ASSOCIATED(atomic_kind_set(ikind)%atom_list)) THEN
             DEALLOCATE (atomic_kind_set(ikind)%atom_list,STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          CALL shell_release(atomic_kind_set(ikind)%shell,error)

          CALL damping_p_release(atomic_kind_set(ikind)%damping,error)
       END DO
       DEALLOCATE (atomic_kind_set,STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer atomic_kind_set is not associated and "//&
            "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_atomic_kind_set

! *****************************************************************************
!> \brief Get attributes of an atomic kind.
!> \param atomic_kind ...
!> \param fist_potential ...
!> \param element_symbol ...
!> \param name ...
!> \param mass ...
!> \param kind_number ...
!> \param natom ...
!> \param atom_list ...
!> \param rcov ...
!> \param rvdw ...
!> \param z ...
!> \param qeff ...
!> \param apol ...
!> \param cpol ...
!> \param mm_radius ...
!> \param shell ...
!> \param shell_active ...
!> \param damping ...
! *****************************************************************************
  SUBROUTINE get_atomic_kind(atomic_kind, fist_potential,&
       element_symbol, name, mass, kind_number, natom, atom_list, &
       rcov, rvdw, z, qeff, apol, cpol, mm_radius, &
       shell, shell_active, damping)

    TYPE(atomic_kind_type)                   :: atomic_kind
    TYPE(fist_potential_type), OPTIONAL, &
      POINTER                                :: fist_potential
    CHARACTER(LEN=2), INTENT(OUT), OPTIONAL  :: element_symbol
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: mass
    INTEGER, INTENT(OUT), OPTIONAL           :: kind_number, natom
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: atom_list
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: rcov, rvdw
    INTEGER, INTENT(OUT), OPTIONAL           :: z
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: qeff, apol, cpol, mm_radius
    TYPE(shell_kind_type), OPTIONAL, POINTER :: shell
    LOGICAL, INTENT(OUT), OPTIONAL           :: shell_active
    TYPE(damping_p_type), OPTIONAL, POINTER  :: damping

    CHARACTER(len=*), PARAMETER :: routineN = 'get_atomic_kind', &
      routineP = moduleN//':'//routineN

    IF (PRESENT(fist_potential)) fist_potential => atomic_kind%fist_potential
    IF (PRESENT(element_symbol)) element_symbol = atomic_kind%element_symbol
    IF (PRESENT(name)) name = atomic_kind%name
    IF (PRESENT(mass)) mass = atomic_kind%mass
    IF (PRESENT(kind_number)) kind_number = atomic_kind%kind_number
    IF (PRESENT(natom)) natom = atomic_kind%natom
    IF (PRESENT(atom_list)) atom_list => atomic_kind%atom_list

    IF (PRESENT(z)) THEN
       CALL get_ptable_info(atomic_kind%element_symbol, number=z)
    END IF
    IF (PRESENT(rcov)) THEN
       CALL get_ptable_info(atomic_kind%element_symbol, covalent_radius=rcov)
    END IF
    IF (PRESENT(rvdw)) THEN
       CALL get_ptable_info(atomic_kind%element_symbol, vdw_radius=rvdw)
    END IF
    IF (PRESENT(qeff)) THEN
       IF (ASSOCIATED(atomic_kind%fist_potential)) THEN
          CALL get_potential(potential=atomic_kind%fist_potential,qeff=qeff)
       ELSE
          qeff = -HUGE(0.0_dp)
       END IF
    END IF
    IF (PRESENT(apol)) THEN
       IF (ASSOCIATED(atomic_kind%fist_potential)) THEN
          CALL get_potential(potential=atomic_kind%fist_potential,apol=apol)
       ELSE
          apol = -HUGE(0.0_dp)
       END IF
    END IF
    IF (PRESENT(cpol)) THEN
       IF (ASSOCIATED(atomic_kind%fist_potential)) THEN
          CALL get_potential(potential=atomic_kind%fist_potential,cpol=cpol)
       ELSE
          cpol = -HUGE(0.0_dp)
       END IF
    END IF
    IF (PRESENT(mm_radius)) THEN
       IF (ASSOCIATED(atomic_kind%fist_potential)) THEN
          CALL get_potential(potential=atomic_kind%fist_potential,mm_radius=mm_radius)
       ELSE
          mm_radius = -HUGE(0.0_dp)
       END IF
    END IF
    IF (PRESENT(shell)) shell => atomic_kind%shell
    IF (PRESENT(shell_active)) shell_active = atomic_kind%shell_active
    IF (PRESENT(damping)) damping => atomic_kind%damping

  END SUBROUTINE get_atomic_kind

! *****************************************************************************
!> \brief Get attributes of an atomic kind set.
!> \param atomic_kind_set ...
!> \param atom_of_kind ...
!> \param kind_of ...
!> \param natom_of_kind ...
!> \param maxatom ...
!> \param natom ...
!> \param nshell ...
!> \param fist_potential_present ...
!> \param shell_present ...
!> \param shell_adiabatic ...
!> \param shell_check_distance ...
!> \param damping_present ...
! *****************************************************************************
  SUBROUTINE get_atomic_kind_set(atomic_kind_set,&
       atom_of_kind, kind_of, natom_of_kind,&
       maxatom,&
       natom,&
       nshell,&
       fist_potential_present,&
       shell_present, shell_adiabatic,&
       shell_check_distance,&
       damping_present)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: atom_of_kind, kind_of, &
                                                natom_of_kind
    INTEGER, INTENT(OUT), OPTIONAL           :: maxatom, natom, nshell
    LOGICAL, INTENT(OUT), OPTIONAL :: fist_potential_present, shell_present, &
      shell_adiabatic, shell_check_distance, damping_present

    CHARACTER(len=*), PARAMETER :: routineN = 'get_atomic_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, iatom, ikind, nkind
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(damping_p_type), POINTER            :: damping
    TYPE(fist_potential_type), POINTER       :: fist_potential
    TYPE(shell_kind_type), POINTER           :: shell

    IF (ASSOCIATED(atomic_kind_set)) THEN

       IF (PRESENT(maxatom)) maxatom = 0
       IF (PRESENT(natom)) natom = 0
       IF (PRESENT(nshell)) nshell = 0
       IF (PRESENT(shell_present)) shell_present = .FALSE.
       IF (PRESENT(shell_adiabatic)) shell_adiabatic = .FALSE.
       IF (PRESENT(shell_check_distance)) shell_check_distance = .FALSE.
       IF (PRESENT(damping_present)) damping_present = .FALSE.
       IF (PRESENT(atom_of_kind)) atom_of_kind(:) = 0
       IF (PRESENT(kind_of)) kind_of(:) = 0
       IF (PRESENT(natom_of_kind)) natom_of_kind(:) = 0

       nkind = SIZE(atomic_kind_set)
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               fist_potential=fist_potential,&
               shell=shell,&
               damping=damping)
          IF (PRESENT(maxatom)) THEN
             maxatom = MAX(maxatom,atomic_kind%natom)
          END IF
          IF (PRESENT(natom)) THEN
             natom = natom + atomic_kind_set(ikind)%natom
          END IF
          IF (PRESENT(fist_potential_present)) THEN
             IF (ASSOCIATED(fist_potential)) THEN
                fist_potential_present = .TRUE.
             END IF
          END IF
          IF (PRESENT(shell_present)) THEN
             IF (ASSOCIATED(shell)) THEN
                shell_present = .TRUE.
             END IF
          END IF
          IF (PRESENT(shell_adiabatic) .AND. ASSOCIATED(shell) )THEN
             IF(.NOT. shell_adiabatic) shell_adiabatic = (shell%massfrac /= 0.0_dp)
          END IF
          IF (PRESENT(shell_check_distance) .AND. ASSOCIATED(shell) )THEN
             IF(.NOT. shell_check_distance) shell_check_distance = (shell%max_dist > 0.0_dp)
          END IF
          IF (PRESENT(damping_present))THEN
             IF (ASSOCIATED(damping)) THEN
                damping_present=.TRUE.
             END IF
          END IF
          IF (PRESENT(atom_of_kind)) THEN
             DO iatom=1,atomic_kind%natom
                atom_a = atomic_kind%atom_list(iatom)
                atom_of_kind(atom_a) = iatom
             END DO
          END IF
          IF (PRESENT(kind_of)) THEN
             DO iatom=1,atomic_kind%natom
                atom_a = atomic_kind%atom_list(iatom)
                kind_of(atom_a) = ikind
             END DO
          END IF
          IF (PRESENT(natom_of_kind)) THEN
             natom_of_kind(ikind) = atomic_kind_set(ikind)%natom
          END IF
       END DO
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer atomic_kind_set is not associated")
    END IF

  END SUBROUTINE get_atomic_kind_set

! *****************************************************************************
!> \brief Set the components of an atomic kind data set.
!> \param atomic_kind ...
!> \param element_symbol ...
!> \param name ...
!> \param mass ...
!> \param kind_number ...
!> \param natom ...
!> \param atom_list ...
!> \param fist_potential ...
!> \param shell ...
!> \param shell_active ...
!> \param damping ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE set_atomic_kind(atomic_kind,element_symbol,name,mass,kind_number,&
       natom,atom_list,&
       fist_potential,shell,&
       shell_active,damping,error)

    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: element_symbol, name
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: mass
    INTEGER, INTENT(IN), OPTIONAL            :: kind_number, natom
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: atom_list
    TYPE(fist_potential_type), OPTIONAL, &
      POINTER                                :: fist_potential
    TYPE(shell_kind_type), OPTIONAL, POINTER :: shell
    LOGICAL, INTENT(IN), OPTIONAL            :: shell_active
    TYPE(damping_p_type), OPTIONAL, POINTER  :: damping
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_atomic_kind', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n, stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(atomic_kind)) THEN

       IF (PRESENT(element_symbol)) atomic_kind%element_symbol = element_symbol
       IF (PRESENT(name)) atomic_kind%name = name
       IF (PRESENT(mass)) atomic_kind%mass = mass
       IF (PRESENT(kind_number)) atomic_kind%kind_number = kind_number
       IF (PRESENT(natom)) atomic_kind%natom = natom
       IF (PRESENT(atom_list)) THEN
          n = SIZE(atom_list)
          IF (n > 0) THEN
             IF (ASSOCIATED(atomic_kind%atom_list)) THEN
                DEALLOCATE (atomic_kind%atom_list,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             ALLOCATE (atomic_kind%atom_list(n),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             atomic_kind%atom_list(:) = atom_list(:)
             atomic_kind%natom = n
          ELSE
             CALL stop_program(routineN,moduleN,__LINE__,&
                  "An invalid atom_list was supplied")
          END IF
       END IF
       IF (PRESENT(fist_potential)) atomic_kind%fist_potential => fist_potential
       IF (PRESENT(shell)) THEN
          atomic_kind%shell => shell
          CALL shell_retain(shell,error)
       END IF
       IF (PRESENT(shell_active)) atomic_kind%shell_active=shell_active

       IF (PRESENT(damping)) atomic_kind%damping => damping

    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer atomic_kind is not associated")
    END IF


  END SUBROUTINE set_atomic_kind

! *****************************************************************************
!> \brief Determines if the atomic_kind is HYDROGEN
!> \param atomic_kind ...
!> \retval res ...
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  PURE FUNCTION is_hydrogen(atomic_kind) RESULT(res)
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    LOGICAL                                  :: res

    res = TRIM(atomic_kind%element_symbol)=="H"
  END FUNCTION is_hydrogen

END MODULE atomic_kind_types
