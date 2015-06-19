!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE qs_fb_atomic_halo_types

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_max,&
                                             mp_min,&
                                             mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: locate,&
                                             sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! public types
  PUBLIC :: fb_atomic_halo_obj,&
            fb_atomic_halo_list_obj

! public methods
  PUBLIC :: fb_atomic_halo_release,&
            fb_atomic_halo_nullify,&
            fb_atomic_halo_has_data,&
            fb_atomic_halo_create,&
            fb_atomic_halo_init,&
            fb_atomic_halo_get,&
            fb_atomic_halo_set,&
            fb_atomic_halo_sort,&
            fb_atomic_halo_atom_global2halo,&
            fb_atomic_halo_nelectrons_estimate_Z,&
            fb_atomic_halo_cost,&
            fb_atomic_halo_build_halo_atoms,&
            fb_atomic_halo_list_retain,&
            fb_atomic_halo_list_release,&
            fb_atomic_halo_list_nullify,&
            fb_atomic_halo_list_has_data,&
            fb_atomic_halo_list_associate,&
            fb_atomic_halo_list_create,&
            fb_atomic_halo_list_get,&
            fb_atomic_halo_list_set,&
            fb_atomic_halo_list_write_info,&
            fb_build_pair_radii


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_atomic_halo_types'
  INTEGER, PRIVATE, SAVE :: last_fb_atomic_halo_id = 0
  INTEGER, PRIVATE, SAVE :: last_fb_atomic_halo_list_id = 0


! *****************************************************************************
!> \brief derived type containing the list of atoms in an atomic halo,
!>        used by filtered-basis diagonalisation method
!> \param id_nr            : unique id for the object
!> \param ref_count        : reference count on the object
!> \param owner_atom       : global atomic id of the atom this halo belongs to
!> \param owner_id_in_halo : index of the owner_atom in the halo_atoms array
!> \param natoms           : number of atoms in the halo
!> \param nelectrons       : estimate of total number of electrons in halo
!> \param halo_atoms       : the list of global id of atoms in the halo
!> \param sorted           : whether the halo_atoms list is sorted or not
!> \param cost             : computational cost for the atomic matrix associated
!>                           to this atomic halo
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_atomic_halo_data
     INTEGER :: id_nr, ref_count
     INTEGER :: owner_atom
     INTEGER :: owner_id_in_halo
     INTEGER :: natoms
     INTEGER :: nelectrons
     INTEGER, DIMENSION(:), POINTER :: halo_atoms
     LOGICAL :: sorted
     REAL(KIND=dp) :: cost
  END TYPE fb_atomic_halo_data


! *****************************************************************************
!> \brief defines a fb_atomic_halo object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_atomic_halo_obj
     TYPE(fb_atomic_halo_data), POINTER, PRIVATE :: obj
  END TYPE fb_atomic_halo_obj


! *****************************************************************************
!> \brief derived type describing an atomic halo list used by
!>        filtered-basis diagonalisation method
!> \param id_nr     : unique id for the object
!> \param ref_count : reference count on the object
!> \param nhalos  : number of halos in the list
!> \param max_nhalos : maximum of the number of halos amongst all of the procs
!> \param halos   : halos(ihalo) gives the ihalo-th fb_atomic_halo object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_atomic_halo_list_data
     INTEGER :: id_nr, ref_count
     INTEGER :: nhalos
     INTEGER :: max_nhalos
     TYPE(fb_atomic_halo_obj), DIMENSION(:), POINTER :: halos
  END TYPE fb_atomic_halo_list_data


! *****************************************************************************
!> \brief defines a fb_atomic_halo_list object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_atomic_halo_list_obj
     TYPE(fb_atomic_halo_list_data), POINTER, PRIVATE :: obj
  END TYPE fb_atomic_halo_list_obj


CONTAINS


! *****************************************************************************
!> \brief Retains an fb_atomic_halo object
!> \param atomic_halo the fb_atomic_halo object, its content must
!>                     not be NULL or UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_retain(atomic_halo, error)
    ! note INTENT(IN) is okay because the obj pointer contained in the
    ! obj type will not be changed
    TYPE(fb_atomic_halo_obj), INTENT(IN)     :: atomic_halo
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halo%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(atomic_halo%obj%ref_count>0, cp_failure_level, routineP, error)
       atomic_halo%obj%ref_count = atomic_halo%obj%ref_count + 1
    END IF
  END SUBROUTINE fb_atomic_halo_retain


! *****************************************************************************
!> \brief Releases an fb_atomic_halo object
!> \param atomic_halo the fb_atomic_halo object, its content must
!>                     not be UNDEFINED, and the subroutine does nothing
!>                     if the content points to NULL
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_release(atomic_halo, error)
    TYPE(fb_atomic_halo_obj), INTENT(INOUT)  :: atomic_halo
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(atomic_halo%obj)) THEN
       CPPreconditionNoFail(atomic_halo%obj%ref_count>0, cp_failure_level, routineP, error)
       atomic_halo%obj%ref_count = atomic_halo%obj%ref_count - 1
       IF (atomic_halo%obj%ref_count == 0) THEN
          atomic_halo%obj%ref_count = 1
          IF (ASSOCIATED(atomic_halo%obj%halo_atoms)) THEN
             ! note that if there are other pointers associated to the memory pointed
             ! by atomic_halo%obj%halo_atoms, their behaviour becomes undefined per
             ! FORTRAN standard (and thus becomes compiler dependent and unreliable)
             ! after the following DEALLOCATE
             DEALLOCATE(atomic_halo%obj%halo_atoms, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          atomic_halo%obj%ref_count = 0
          DEALLOCATE(atomic_halo%obj, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    ELSE
       NULLIFY(atomic_halo%obj)
    END IF
  END SUBROUTINE fb_atomic_halo_release


! *****************************************************************************
!> \brief Nullifies a fb_atomic_halo object, note that it does not
!>        release the original object. This procedure is used to nullify
!>        the pointer contained in the object which is used to associate
!>        to the actual object content
!> \param atomic_halo the fb_atomic_halo object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_nullify(atomic_halo)
    TYPE(fb_atomic_halo_obj), INTENT(INOUT)  :: atomic_halo

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(atomic_halo%obj)
  END SUBROUTINE fb_atomic_halo_nullify


! *****************************************************************************
!> \brief Associates one fb_atomic_halo object to another
!> \param a the fb_atomic_halo object to be associated
!> \param b the fb_atomic_halo object that a is to be associated to
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_associate(a, b)
    TYPE(fb_atomic_halo_obj), INTENT(OUT)    :: a
    TYPE(fb_atomic_halo_obj), INTENT(IN)     :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_atomic_halo_associate


! *****************************************************************************
!> \brief Checks if a fb_atomic_halo object is associated with an actual
!>        data content or not
!> \param atomic_halo the fb_atomic_halo object
!> \retval res ...
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  FUNCTION fb_atomic_halo_has_data(atomic_halo) RESULT(res)
    TYPE(fb_atomic_halo_obj), INTENT(IN)     :: atomic_halo
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_has_data', &
      routineP = moduleN//':'//routineN

    res = ASSOCIATED(atomic_halo%obj)
  END FUNCTION fb_atomic_halo_has_data


! *****************************************************************************
!> \brief Creates and initialises an empty fb_atomic_halo object
!> \param atomic_halo the fb_atomic_halo object, its content must
!>                     be NULL and cannot be UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_create(atomic_halo, error)
    TYPE(fb_atomic_halo_obj), INTENT(INOUT)  :: atomic_halo
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(atomic_halo%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(atomic_halo%obj, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       atomic_halo%obj%owner_atom = 0
       atomic_halo%obj%owner_id_in_halo = 0
       atomic_halo%obj%natoms = 0
       atomic_halo%obj%nelectrons = 0
       atomic_halo%obj%sorted = .FALSE.
       atomic_halo%obj%cost = 0.0_dp
       NULLIFY(atomic_halo%obj%halo_atoms)
       atomic_halo%obj%ref_count = 1
       atomic_halo%obj%id_nr = last_fb_atomic_halo_id + 1
       last_fb_atomic_halo_id = atomic_halo%obj%id_nr
    END IF
  END SUBROUTINE fb_atomic_halo_create


! *****************************************************************************
!> \brief Initialises an fb_atomic_halo object, and makes it empty
!> \param atomic_halo the fb_atomic_halo object, its content must
!>                     not be NULL or UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_init(atomic_halo, error)
    TYPE(fb_atomic_halo_obj), INTENT(INOUT)  :: atomic_halo
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halo%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ! if halo_atoms are associated, then deallocate and de-associate
       IF (ASSOCIATED(atomic_halo%obj%halo_atoms)) THEN
          DEALLOCATE(atomic_halo%obj%halo_atoms, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          NULLIFY(atomic_halo%obj%halo_atoms)
       END IF
       atomic_halo%obj%owner_atom = 0
       atomic_halo%obj%owner_id_in_halo = 0
       atomic_halo%obj%natoms = 0
       atomic_halo%obj%nelectrons = 0
       atomic_halo%obj%sorted = .FALSE.
       atomic_halo%obj%cost = 0.0_dp
    END IF
  END SUBROUTINE fb_atomic_halo_init


! *****************************************************************************
!> \brief Gets attributes from a fb_atomic_halo object, one should
!>        only access the data content in a fb_atomic_halo outside
!>        this module via this procedure.
!> \param atomic_halo the fb_atomic_halo object, its content must
!>                     not be NULL or UNDEFINED
!> \param owner_atom [OPTIONAL]: if present, outputs atmic_halo%obj%owner_atom
!> \param owner_id_in_halo ...
!> \param natoms [OPTIONAL]: if present, outputs atomic_halo%obj%natoms
!> \param nelectrons [OPTIONAL]: if present, outputs atomic_halo%obj%nelectrons
!> \param halo_atoms [OPTIONAL]: if present, outputs pointer
!>                               atomic_halo%obj%halo_atoms
!> \param sorted [OPTIONAL]: if present, outputs atomic_halo%obj%sorted
!> \param cost [OPTIONAL]: if present, outputs atomic_halo%obj%cost
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_get(atomic_halo, &
                                owner_atom, &
                                owner_id_in_halo, &
                                natoms, &
                                nelectrons, &
                                halo_atoms, &
                                sorted, &
                                cost, &
                                error)
    TYPE(fb_atomic_halo_obj), INTENT(IN)     :: atomic_halo
    INTEGER, INTENT(OUT), OPTIONAL           :: owner_atom, owner_id_in_halo, &
                                                natoms, nelectrons
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: halo_atoms
    LOGICAL, INTENT(OUT), OPTIONAL           :: sorted
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: cost
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halo%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(owner_atom)) owner_atom = atomic_halo%obj%owner_atom
       IF (PRESENT(owner_id_in_halo)) owner_id_in_halo = atomic_halo%obj%owner_id_in_halo
       IF (PRESENT(natoms)) natoms = atomic_halo%obj%natoms
       IF (PRESENT(nelectrons)) nelectrons = atomic_halo%obj%nelectrons
       IF (PRESENT(halo_atoms)) halo_atoms => atomic_halo%obj%halo_atoms
       IF (PRESENT(sorted)) sorted = atomic_halo%obj%sorted
       IF (PRESENT(cost)) cost = atomic_halo%obj%cost
    END IF
  END SUBROUTINE fb_atomic_halo_get


! *****************************************************************************
!> \brief Sets attributes in a fb_atomic_halo object, one should
!>        only set the data content in a fb_atomic_halo from outside
!>        this module via this procedure.
!> \param atomic_halo the fb_atomic_halo object, its content must
!>                     not be NULL or UNDEFINED
!> \param owner_atom [OPTIONAL]: if present, sets
!>                               atmic_halo%obj%owner_atom = owner_atom
!> \param owner_id_in_halo ...
!> \param natoms [OPTIONAL]: if present, sets atomic_halo%obj%natoms = natoms
!> \param nelectrons [OPTIONAL]: if present, sets atomic_halo%obj%nelectrons = nelectrons
!> \param halo_atoms [OPTIONAL]: if present, reallocates atomic_halo%obj%halo_atoms
!>                               to the size of halo_atoms, and copies
!>                               contents of halo_atoms to atomic_halo%obj%halo_atoms
!> \param sorted [OPTIONAL]: if present, sets atomic_halo%obj%sorted = sorted
!> \param cost [OPTIONAL]: if present, sets atomic_halo%obj%cost = cost
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_set(atomic_halo, &
                                owner_atom, &
                                owner_id_in_halo, &
                                natoms, &
                                nelectrons, &
                                halo_atoms, &
                                sorted, &
                                cost, &
                                error)
    TYPE(fb_atomic_halo_obj), INTENT(INOUT)  :: atomic_halo
    INTEGER, INTENT(IN), OPTIONAL            :: owner_atom, owner_id_in_halo, &
                                                natoms, nelectrons
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: halo_atoms
    LOGICAL, INTENT(IN), OPTIONAL            :: sorted
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: cost
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halo%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(owner_atom)) atomic_halo%obj%owner_atom = owner_atom
       IF (PRESENT(owner_id_in_halo)) atomic_halo%obj%owner_id_in_halo = owner_id_in_halo
       IF (PRESENT(natoms)) atomic_halo%obj%natoms = natoms
       IF (PRESENT(nelectrons)) atomic_halo%obj%nelectrons = nelectrons
       IF (PRESENT(halo_atoms)) THEN
          IF (ASSOCIATED(atomic_halo%obj%halo_atoms)) THEN
             DEALLOCATE(atomic_halo%obj%halo_atoms, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          atomic_halo%obj%halo_atoms => halo_atoms
       END IF
       IF (PRESENT(nelectrons)) atomic_halo%obj%nelectrons = nelectrons
       IF (PRESENT(sorted)) atomic_halo%obj%sorted = sorted
       IF (PRESENT(cost)) atomic_halo%obj%cost = cost
    END IF
  END SUBROUTINE fb_atomic_halo_set


! *****************************************************************************
!> \brief Sort the list of atomic indices in the halo in accending order.
!>        The atomic_halo must not be empty
!> \param atomic_halo  the atomic_halo object
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_sort(atomic_halo, error)
    TYPE(fb_atomic_halo_obj), INTENT(INOUT)  :: atomic_halo
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_sort', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_index
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(SIZE(atomic_halo%obj%halo_atoms) > 0, cp_failure_level, routineP, error, failure)
    ALLOCATE(tmp_index(atomic_halo%obj%natoms), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL sort(atomic_halo%obj%halo_atoms, atomic_halo%obj%natoms, tmp_index)
    DEALLOCATE(tmp_index, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    atomic_halo%obj%sorted = .TRUE.
  END SUBROUTINE fb_atomic_halo_sort


! *****************************************************************************
!> \brief Given a global atomic index, convert it to its index in a
!>        given atomic halo, if found.
!>        The atomic_halo object must already have been sorted
!> \param atomic_halo  the atomic_halo object
!> \param iatom_global  the global atomic index
!> \param iatom_halo  the atomic index inside the halo
!> \param found  returns true if given atom is in the halo, otherwise
!>                false
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_atom_global2halo(atomic_halo, &
                                             iatom_global, &
                                             iatom_halo, &
                                             found, &
                                             error)
    TYPE(fb_atomic_halo_obj), INTENT(IN)     :: atomic_halo
    INTEGER, INTENT(IN)                      :: iatom_global
    INTEGER, INTENT(OUT)                     :: iatom_halo
    LOGICAL, INTENT(OUT)                     :: found
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'fb_atomic_halo_atom_global2halo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN, handle)

    failure = .FALSE.
    CPPrecondition(atomic_halo%obj%sorted, cp_failure_level, routineP, error, failure)
    iatom_halo = locate(atomic_halo%obj%halo_atoms, iatom_global)
    IF (iatom_halo == 0) THEN
       found = .FALSE.
    ELSE
       found = .TRUE.
    END IF

    CALL timestop(handle)

  END SUBROUTINE fb_atomic_halo_atom_global2halo


! *****************************************************************************
!> \brief Estimates the total number of electrons in a halo using atomic
!>        numbers
!> \param atomic_halo the atomic_halo object
!> \param particle_set an array of cp2k particle set objects (this
!>                       gives atomic information)
!> \retval nelectrons estimate of electron number
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  FUNCTION fb_atomic_halo_nelectrons_estimate_Z(atomic_halo, particle_set) RESULT(nelectrons)
    TYPE(fb_atomic_halo_obj), INTENT(IN)     :: atomic_halo
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: particle_set
    INTEGER                                  :: nelectrons

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'fb_atomic_halo_nelectrons_estimate_Z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom_global, iatom_halo, z
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    failure = .FALSE.
    nelectrons = 0
    IF (ASSOCIATED(atomic_halo%obj)) THEN
       DO iatom_halo = 1, atomic_halo%obj%natoms
          iatom_global = atomic_halo%obj%halo_atoms(iatom_halo)
          atomic_kind => particle_set(iatom_global)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind, &
                               z=z)
          nelectrons = nelectrons + z
       END DO
    END IF
  END FUNCTION fb_atomic_halo_nelectrons_estimate_Z


! *****************************************************************************
!> \brief Estimates the computational cost with respect to the filter matrix
!>        calculation associated to an atomic halo. Given the bottle neck
!>        of the filter matrix generation will be the diagoanlisation of the
!>        atomic matrices (each consists of atoms in an atomic halo), the cost
!>        can be estimated by counting the total number of contracted gaussians
!>        in the halo
!> \param atomic_halo  : the atomic_halo object in question
!> \param particle_set : an array of cp2k particle set objects, this
!>                       provides atomic information
!> \param qs_kind_set  : cp2k qs_kind objects, provides information on the
!>                       number of contracted gaussian functions each kind
!>                       has
!> \retval cost        : computation cost w.r.t. the filter matrix
!>                       calculation for this atomic halo
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  FUNCTION fb_atomic_halo_cost(atomic_halo, &
                               particle_set, &
                               qs_kind_set) &
    RESULT(cost)
    TYPE(fb_atomic_halo_obj), INTENT(IN)     :: atomic_halo
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      INTENT(IN)                             :: qs_kind_set
    REAL(KIND=dp)                            :: cost

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_cost', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ii, ikind, ncgf

    cost = 0.0_dp
    DO ii = 1, atomic_halo%obj%natoms
       iatom = atomic_halo%obj%halo_atoms(ii)
       CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                            kind_number=ikind)
       CALL get_qs_kind(qs_kind=qs_kind_set(ikind), &
                        ncgf=ncgf)
       cost = cost + REAL(ncgf, dp)
    END DO
    ! diagonalisation is N**3 process, so cost must reflect that
    cost = cost**3
  END FUNCTION fb_atomic_halo_cost


! *****************************************************************************
!> \brief Builds halo atoms for a given (owner) atom
!> \param owner_atom   : the atom the halo is going to be built for
!> \param particle_set : an array of cp2k particle set objects, this
!>                       provides atomic information
!> \param cell         : cp2k cell object, used for resolving periodic
!>                       boundary conditions
!> \param pair_radii   : 2D array storing interaction radii between two kinds
!> \param halo_atoms   : must be NULL pointer on input, and outputs an
!>                       array of halo atoms corresponding to the owner atom
!> \param nhalo_atoms  : outputs number of halo atoms
!> \param owner_id_in_halo : the index of the owner atom in the halo_atoms list
!> \param error        : CP2K error handler container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_build_halo_atoms(owner_atom, &
                                             particle_set, &
                                             cell, &
                                             pair_radii, &
                                             halo_atoms, &
                                             nhalo_atoms, &
                                             owner_id_in_halo, &
                                             error)
    INTEGER, INTENT(IN)                      :: owner_atom
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: pair_radii
    INTEGER, DIMENSION(:), POINTER           :: halo_atoms
    INTEGER, INTENT(OUT)                     :: nhalo_atoms, owner_id_in_halo
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'fb_atomic_halo_build_halo_atoms', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ikind, jatom, jkind, &
                                                natoms_global, stat
    LOGICAL                                  :: check_ok, failure
    REAL(KIND=dp)                            :: rij
    REAL(KIND=dp), DIMENSION(3)              :: ri, rij_pbc, rj
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    failure = .FALSE.
    check_ok = .NOT. ASSOCIATED(halo_atoms)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)

    NULLIFY(atomic_kind)

    iatom = owner_atom
    atomic_kind => particle_set(iatom)%atomic_kind
    CALL get_atomic_kind(atomic_kind=atomic_kind, &
                         kind_number=ikind)
    natoms_global = SIZE(particle_set)
    ALLOCATE(halo_atoms(natoms_global), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    owner_id_in_halo = 0
    nhalo_atoms = 0
    DO jatom = 1, natoms_global
       atomic_kind => particle_set(jatom)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, &
                            kind_number=jkind)
       ! calculate the minimum distance between iatom and
       ! jatom, taking account of the periodic boundary
       ! conditions
       ri(1:3) = particle_set(iatom)%r(1:3)
       rj(1:3) = particle_set(jatom)%r(1:3)
       rij_pbc = pbc(ri, rj, cell)
       rij = rij_pbc(1)*rij_pbc(1) + &
             rij_pbc(2)*rij_pbc(2) + &
             rij_pbc(3)*rij_pbc(3)
       rij = SQRT(rij)
       IF (rij .LE. pair_radii(ikind, jkind)) THEN
          ! jatom is in iatom's halo
          nhalo_atoms = nhalo_atoms + 1
          halo_atoms(nhalo_atoms) = jatom
          IF (jatom == iatom) owner_id_in_halo = nhalo_atoms
       END IF
    END DO
    CALL reallocate(halo_atoms, 1, nhalo_atoms)
  END SUBROUTINE fb_atomic_halo_build_halo_atoms


! *****************************************************************************
!> \brief Retains an fb_atomic_halo_list object
!> \param atomic_halos the fb_atomic_halo object, its content must
!>                      not be NULL or UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_list_retain(atomic_halos, error)
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(IN)                             :: atomic_halos
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_list_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halos%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(atomic_halos%obj%ref_count>0, cp_failure_level, routineP, error)
       atomic_halos%obj%ref_count = atomic_halos%obj%ref_count + 1
    END IF
  END SUBROUTINE fb_atomic_halo_list_retain


! *****************************************************************************
!> \brief Releases an fb_atomic_halo_list object
!> \param atomic_halos the fb_atomic_halo object, its content must
!>                      not be UNDEFINED, and does nothing if it is NULL
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_list_release(atomic_halos, error)
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(INOUT)                          :: atomic_halos
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_list_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(atomic_halos%obj)) THEN
       CPPreconditionNoFail(atomic_halos%obj%ref_count>0, cp_failure_level, routineP, error)
       atomic_halos%obj%ref_count = atomic_halos%obj%ref_count - 1
       IF (atomic_halos%obj%ref_count == 0) THEN
          atomic_halos%obj%ref_count = 1
          IF (ASSOCIATED(atomic_halos%obj%halos)) THEN
             DO ii = 1, SIZE(atomic_halos%obj%halos)
                CALL fb_atomic_halo_release(atomic_halos%obj%halos(ii), error=error)
             END DO
             DEALLOCATE(atomic_halos%obj%halos, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          atomic_halos%obj%ref_count = 0
          DEALLOCATE(atomic_halos%obj, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    ELSE
       NULLIFY(atomic_halos%obj)
    END IF
  END SUBROUTINE fb_atomic_halo_list_release


! *****************************************************************************
!> \brief Nullifies a fb_atomic_halo_list object, note that it does
!>        not release the original object. This procedure is used to
!>        nullify the pointer contained in the object which is used to
!>        associate to the actual object content
!> \param atomic_halos the fb_atomic_halo_list object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_list_nullify(atomic_halos)
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(INOUT)                          :: atomic_halos

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_list_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(atomic_halos%obj)
  END SUBROUTINE fb_atomic_halo_list_nullify


! *****************************************************************************
!> \brief Checks if a fb_atomic_halo_list object is associated with
!>        an actual data content or not
!> \param atomic_halos the fb_atomic_halo_list object
!> \retval res ...
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  FUNCTION fb_atomic_halo_list_has_data(atomic_halos) RESULT(res)
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(IN)                             :: atomic_halos
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_list_has_data', &
      routineP = moduleN//':'//routineN

    res = ASSOCIATED(atomic_halos%obj)
  END FUNCTION fb_atomic_halo_list_has_data


! *****************************************************************************
!> \brief Associates one fb_atomic_halo_list object to another
!> \param a the fb_atomic_halo_list object to be associated
!> \param b the fb_atomic_halo_list object that a is to be associated to
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_list_associate(a, b)
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(OUT)                            :: a
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(IN)                             :: b

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'fb_atomic_halo_list_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_atomic_halo_list_associate


! *****************************************************************************
!> \brief Creates and initialises an empty fb_atomic_halo_list object
!> \param atomic_halos the fb_atomic_halo object, its content must
!>                      not be NULL or UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_list_create(atomic_halos, error)
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(INOUT)                          :: atomic_halos
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_list_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(atomic_halos%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(atomic_halos%obj, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       atomic_halos%obj%nhalos = 0
       atomic_halos%obj%max_nhalos = 0
       NULLIFY(atomic_halos%obj%halos)
       atomic_halos%obj%ref_count = 1
       atomic_halos%obj%id_nr = last_fb_atomic_halo_list_id + 1
       last_fb_atomic_halo_list_id = atomic_halos%obj%id_nr
    END IF
  END SUBROUTINE fb_atomic_halo_list_create


! *****************************************************************************
!> \brief Initialises an fb_atomic_halo_list object and make it empty
!> \param atomic_halos the fb_atomic_halo object, its content must
!>                      not be NULL or UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_list_init(atomic_halos, error)
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(INOUT)                          :: atomic_halos
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_list_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halos%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ! if the arrays are associated, then deallocate and de-associate
       IF (ASSOCIATED(atomic_halos%obj%halos)) THEN
          DO ii = 1, SIZE(atomic_halos%obj%halos)
             CALL fb_atomic_halo_release(atomic_halos%obj%halos(ii), error)
          END DO
          DEALLOCATE(atomic_halos%obj%halos, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          NULLIFY(atomic_halos%obj%halos)
       END IF
       atomic_halos%obj%nhalos = 0
       atomic_halos%obj%max_nhalos = 0
    END IF
  END SUBROUTINE fb_atomic_halo_list_init


! *****************************************************************************
!> \brief Gets attributes from an fb_atomic_halo_list object, one should
!>        only access the data content in a fb_atomic_halo_list outside
!>        this module via this procedure.
!> \param atomic_halos the fb_atomic_halo object, its content must
!>                      not be NULL or UNDEFINED
!> \param nhalos [OPTIONAL]: if present, gives nhalos = atomic_halos%obj%nhalos
!> \param max_nhalos [OPTIONAL]: if present, gives max_nhalos = atomic_halos%obj%max_nhalos
!> \param halos [OPTIONAL]: if present, gives halos => atomic_halos%obj%halos
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_list_get(atomic_halos, nhalos, max_nhalos, halos, error)
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(IN)                             :: atomic_halos
    INTEGER, INTENT(OUT), OPTIONAL           :: nhalos, max_nhalos
    TYPE(fb_atomic_halo_obj), DIMENSION(:), &
      OPTIONAL, POINTER                      :: halos
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_list_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halos%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(nhalos)) nhalos = atomic_halos%obj%nhalos
       IF (PRESENT(max_nhalos)) max_nhalos = atomic_halos%obj%max_nhalos
       IF (PRESENT(halos)) halos => atomic_halos%obj%halos
    END IF
  END SUBROUTINE fb_atomic_halo_list_get


! *****************************************************************************
!> \brief Sets attributes from an fb_atomic_halo_list object, one should
!>        only set the data content in a fb_atomic_halo_list outside
!>        this module via this procedure.
!> \param atomic_halos the fb_atomic_halo object, its content must
!>                      not be NULL or UNDEFINED
!> \param nhalos [OPTIONAL]: if present, sets atomic_halos%obj%nhalos = nhalos
!> \param max_nhalos [OPTIONAL]: if present, sets atomic_halos%obj%max_nhalos = max_nhalos
!> \param halos [OPTIONAL]: if present, reallocates atomic_halos%obj%halos
!>                          to the size of halos
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_list_set(atomic_halos, nhalos, max_nhalos, halos, error)
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(INOUT)                          :: atomic_halos
    INTEGER, INTENT(IN), OPTIONAL            :: nhalos, max_nhalos
    TYPE(fb_atomic_halo_obj), DIMENSION(:), &
      OPTIONAL, POINTER                      :: halos
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_list_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ihalo, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atomic_halos%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(nhalos)) atomic_halos%obj%nhalos = nhalos
       IF (PRESENT(max_nhalos)) atomic_halos%obj%max_nhalos = max_nhalos
       IF (PRESENT(halos)) THEN
          IF (ASSOCIATED(atomic_halos%obj%halos)) THEN
             DO ihalo = 1, SIZE(atomic_halos%obj%halos)
                CALL fb_atomic_halo_release(atomic_halos%obj%halos(ihalo), error)
             END DO
             DEALLOCATE(atomic_halos%obj%halos, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          atomic_halos%obj%halos => halos
       END IF
    END IF
  END SUBROUTINE fb_atomic_halo_list_set


! *****************************************************************************
!> \brief Writes out the atomic halo list from an fb_atomic_halo_list
!>        object using information
!> \param atomic_halos the fb_atomic_halo object
!> \param para_env pointer to a para_env_type object containing MPI info
!> \param fb_section pointer to the input section to filtered basis method
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_list_write(atomic_halos, para_env, fb_section, error)
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(IN)                             :: atomic_halos
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: fb_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_atomic_halo_list_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: string
    INTEGER                                  :: ihalo, jatom, mype, &
                                                nhalo_atoms, nhalos, &
                                                owner_atom, print_unit
    INTEGER, DIMENSION(:), POINTER           :: halo_atoms
    LOGICAL                                  :: failure, new_file
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fb_atomic_halo_obj), DIMENSION(:), &
      POINTER                                :: halos

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info, fb_section, &
                                         "PRINT%ATOMIC_HALOS", error=error), &
              cp_p_file)) THEN
       print_unit = cp_print_key_unit_nr(logger=logger, &
                                         basis_section=fb_section, &
                                         print_key_path="PRINT%ATOMIC_HALOS", &
                                         extension=".out", &
                                         local=.TRUE., &
                                         log_filename=.FALSE., &
                                         file_position="REWIND", &
                                         file_action="WRITE", &
                                         is_new_file=new_file, &
                                         error=error)
       mype = para_env%mepos
       ! print headline
       string = ""
       WRITE (UNIT=string, FMT="(A,I5,A)") &
             "ATOMIC HALOS IN (PROCESS ", mype, ")"
       CALL compress(string)
       IF (print_unit > 0) THEN
          WRITE (UNIT=print_unit, FMT="(/,/,T2,A)") TRIM(string)
          WRITE (UNIT=print_unit, FMT="(/,T2,A)") &
                "atom : list of atoms in the atomic halo"
       END IF
       ! print content
       CALL fb_atomic_halo_list_get(atomic_halos=atomic_halos, &
                                    nhalos=nhalos, &
                                    halos=halos, &
                                    error=error)
       DO ihalo = 1, nhalos
          CALL fb_atomic_halo_get(halos(ihalo), &
                                  owner_atom=owner_atom, &
                                  natoms=nhalo_atoms, &
                                  halo_atoms=halo_atoms, &
                                  error=error)
          WRITE (UNIT=print_unit, FMT="(2X,I6,A)", ADVANCE="no") &
                owner_atom, " : "
          DO jatom = 1, nhalo_atoms
             WRITE (UNIT=print_unit, FMT="(I6)", ADVANCE="no") &
                   halo_atoms(jatom)
          END DO
          WRITE (UNIT=print_unit) ""
       END DO
       ! finish
       CALL cp_print_key_finished_output(print_unit, logger, fb_section, &
                                         "PRINT%ATOMIC_HALOS", error=error)
    END IF
  END SUBROUTINE fb_atomic_halo_list_write


! *****************************************************************************
!> \brief Writes out the atomic halo list summary, no detailed neighbour lists,
!>        just average, min and max number of halo atoms in the halo list
!> \param atomic_halos : the fb_atomic_halo object
!> \param para_env     : pointer to a para_env_type object containing MPI info
!> \param scf_section  : pointer to the scf input section
!> \param error        : CP2K data container for error container
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atomic_halo_list_write_info(atomic_halos, para_env, scf_section, error)
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(IN)                             :: atomic_halos
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'fb_atomic_halo_list_write_info', &
      routineP = moduleN//':'//routineN

    INTEGER :: ihalo, max_natoms, min_natoms, nhalo_atoms, nhalos, &
      total_n_halo_atoms, total_n_halos, unit_nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ave_natoms
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fb_atomic_halo_obj), DIMENSION(:), &
      POINTER                                :: halos

    failure = .FALSE.
    NULLIFY(logger, halos)
    logger => cp_error_get_logger(error)
    unit_nr = cp_print_key_unit_nr(logger, scf_section, &
                                   "PRINT%FILTER_MATRIX", &
                                   extension="", &
                                   error=error)

    ! obtain data
    CALL fb_atomic_halo_list_get(atomic_halos=atomic_halos, &
         halos=halos, &
         nhalos=nhalos, &
         error=error)
    max_natoms = 0
    min_natoms = HUGE(0)
    total_n_halo_atoms = 0
    total_n_halos = nhalos
    DO ihalo = 1, nhalos
       CALL fb_atomic_halo_get(atomic_halo=halos(ihalo), &
            natoms=nhalo_atoms, &
            error=error)
       total_n_halo_atoms = total_n_halo_atoms + nhalo_atoms
       max_natoms = MAX(max_natoms, nhalo_atoms)
       min_natoms = MIN(min_natoms, nhalo_atoms)
    END DO
    CALL mp_max(max_natoms, para_env%group)
    CALL mp_min(min_natoms, para_env%group)
    CALL mp_sum(total_n_halos, para_env%group)
    CALL mp_sum(total_n_halo_atoms, para_env%group)
    ave_natoms = REAL(total_n_halo_atoms,dp) / REAL(total_n_halos,dp)
    ! write info
    IF (unit_nr > 0) THEN
       WRITE (UNIT=unit_nr, FMT="(/,A)") &
            " FILTER_MAT_DIAG| Atomic matrix neighbor lists information:"
       WRITE (UNIT=unit_nr, FMT="(A,I10)") &
            " FILTER_MAT_DIAG|   Number of atomic matrices: ", &
            total_n_halos
       WRITE (UNIT=unit_nr, &
              FMT="(A,T45,A,T57,A,T69,A,T81,A)") &
            " FILTER_MAT_DIAG| ", "Average", "Max", "Min"
       WRITE (UNIT=unit_nr, &
              FMT="(A,T45,F10.1,T57,I10,T69,I10,T81,I10)") &
            " FILTER_MAT_DIAG|   N neighbors per atom:", &
            ave_natoms, max_natoms, min_natoms
    END IF
    ! finish
    CALL cp_print_key_finished_output(unit_nr, logger, scf_section, &
                                      "PRINT%FILTER_MATRIX", error=error)
  END SUBROUTINE fb_atomic_halo_list_write_info


! *****************************************************************************
!> \brief Builds the required pair_radii array required for building the
!>        halo atoms from a given set of cut off radii
!> \param rcut   : rcut(ikind) is the cutoff radii for determining the halo
!>                 corresponding to atomic kind ikind
!> \param nkinds : total number of atomic kinds in rcut
!> \param pair_radii : output array pair_radii(ikind,jkind) gives the
!>                     corresponding interaction range between a pair of atoms
!>                     of kinds ikind and jkind
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  PURE SUBROUTINE fb_build_pair_radii(rcut, nkinds, pair_radii)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rcut
    INTEGER, INTENT(IN)                      :: nkinds
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: pair_radii

    INTEGER                                  :: ii, jj

    pair_radii = 0.0_dp
    DO ii = 1, nkinds
       DO jj = 1, nkinds
          pair_radii(ii,jj) = rcut(ii) + rcut(jj)
       END DO
    END DO
  END SUBROUTINE fb_build_pair_radii

END MODULE qs_fb_atomic_halo_types
