MODULE qs_fb_com_tasks_types

  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: int_4,&
                                             int_8
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_alltoall
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! public parameters:
  PUBLIC :: TASK_N_RECORDS, &
            TASK_DEST,      &
            TASK_SRC,       &
            TASK_PAIR,      &
            TASK_COST

! public types
  PUBLIC :: fb_com_tasks_obj, &
            fb_com_atom_pairs_obj

! public methods
  PUBLIC :: fb_com_tasks_retain,             &
            fb_com_tasks_release,            &
            fb_com_tasks_nullify,            &
            fb_com_tasks_has_data,           &
            fb_com_tasks_associate,          &
            fb_com_tasks_create,             &
            fb_com_tasks_init,               &
            fb_com_tasks_get,                &
            fb_com_tasks_set,                &
            fb_com_tasks_transpose_dest_src, &
            fb_com_tasks_build_atom_pairs,   &
            fb_com_tasks_encode_pair,        &
            fb_com_tasks_decode_pair,        &
            fb_com_atom_pairs_retain,        &
            fb_com_atom_pairs_release,       &
            fb_com_atom_pairs_nullify,       &
            fb_com_atom_pairs_has_data,      &
            fb_com_atom_pairs_associate,     &
            fb_com_atom_pairs_create,        &
            fb_com_atom_pairs_init,          &
            fb_com_atom_pairs_get,           &
            fb_com_atom_pairs_set,           &
            fb_com_atom_pairs_encode,        &
            fb_com_atom_pairs_decode,        &
            fb_com_atom_pairs_calc_buffer_sizes

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_com_tasks_types'
  INTEGER, PRIVATE, SAVE :: last_fb_com_tasks_id = 0
  INTEGER, PRIVATE, SAVE :: last_fb_com_atom_pairs_id = 0

! **********************************************************************
! explaination on format of task lists (same for tasks_recv and tasks_send):
! tasks_recv has dimension (4, ntasks_recv), and stores information on
! the block to be copied or transfered
! - tasks_recv(TASK_DEST,itask) = destination MPI rank of itask-th task
! - tasks_recv(TASK_SRC,itask) = source MPI rank of itask-th task
! - tasks_recv(TASK_PAIR,itask) = compressed pair indices of the block of itask-th task
! - tasks_recv(TASK_COST,itask) = the cost of itask-th task
!
! number of record slots in each task in the task lists
  INTEGER, PARAMETER :: TASK_N_RECORDS = 4
! the indices for the records (1:TASK_DIM) in a task
  INTEGER, PARAMETER :: TASK_DEST = 1, &
                        TASK_SRC  = 2, &
                        TASK_PAIR = 3, &
                        TASK_COST = 4
! **********************************************************************


! **********************************************************************
!> \brief data content for communication tasks used for send and receive
!>        matrix blocks
!> \param id_nr     : unique id for the object
!> \param ref_count : reference count on the object
!> \param tasks     : the list of communication tasks, which is
!>                    represented by a 2D array, first dim stores
!>                    info for the communication: src and desc procs
!>                    and the atomic pair indexing the matrix block
!>                    to be communicated, etc.
!> \param task_dim  : the size of the first dimension of tasks
!> \param ntasks    : total number of local tasks
!> \param nencode   : the total number of atoms used for encoding
!>                    the block coordinates (iatom, jatom)
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **********************************************************************
  TYPE fb_com_tasks_data
     INTEGER :: id_nr, ref_count
     ! use pure integer arrays to facilitate easier MPI coms
     INTEGER(KIND=int_8), DIMENSION(:,:), POINTER :: tasks
     INTEGER :: task_dim
     INTEGER :: ntasks
     INTEGER :: nencode
  END TYPE fb_com_tasks_data


!**********************************************************************
!> \brief defines a fb_com_tasks object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
!**********************************************************************
  TYPE fb_com_tasks_obj
     TYPE(fb_com_tasks_data), POINTER, PRIVATE :: obj
  END TYPE fb_com_tasks_obj


! **********************************************************************
!> \brief data content for the list of block coordinates with the
!>        associated src/dest proc id for communication. These will be
!>        generated from the fb_com_tasks object
!> \param id_nr         : unique id for the object
!> \param ref_count     : reference count on the object
!> \param pairs         : the list of communication tasks, which is
!>                        represented by a 2D array, first dim stores
!>                        info for the communication: src and desc procs
!>                        and the atomic pair indexing the matrix block
!>                        to be communicated, etc.
!> \param npairs        : number of blks to be communicated in the atom
!>                        pair list
!> \param natoms_encode : the total number of atoms used for encoding
!>                        the proc + block coordinates (pe, iatom, jatom)
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **********************************************************************
  TYPE fb_com_atom_pairs_data
     INTEGER :: id_nr, ref_count
     INTEGER(KIND=int_8), DIMENSION(:), POINTER :: pairs
     INTEGER :: npairs
     INTEGER :: natoms_encode
  END TYPE fb_com_atom_pairs_data


! **********************************************************************
!> \brief defines a fb_com_atom_pairs object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **********************************************************************
  TYPE fb_com_atom_pairs_obj
     TYPE(fb_com_atom_pairs_data), POINTER, PRIVATE :: obj
  END TYPE fb_com_atom_pairs_obj


CONTAINS


! **********************************************************************
!> \brief Retains an fb_com_tasks object
!> \param com_tasks the fb_com_tasks object, its content must not be
!>                   NULL or UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_retain(com_tasks, error)
    TYPE(fb_com_tasks_obj), INTENT(IN)       :: com_tasks
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(com_tasks%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(com_tasks%obj%ref_count>0, cp_failure_level, routineP, error)
       com_tasks%obj%ref_count = com_tasks%obj%ref_count + 1
    END IF
  END SUBROUTINE fb_com_tasks_retain


! **********************************************************************
!> \brief Retains an fb_com_atom_pairs object
!> \param atom_pairs the fb_com_atom_pairs object, its content must not be
!>                   NULL or UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_atom_pairs_retain(atom_pairs, error)
    TYPE(fb_com_atom_pairs_obj), INTENT(IN)  :: atom_pairs
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_atom_pairs_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atom_pairs%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(atom_pairs%obj%ref_count>0, cp_failure_level, routineP, error)
       atom_pairs%obj%ref_count = atom_pairs%obj%ref_count + 1
    END IF
  END SUBROUTINE fb_com_atom_pairs_retain


! **********************************************************************
!> \brief Releases an fb_com_tasks object
!> \param com_tasks the fb_com_tasks object, its content must not be
!>                   UNDEFINED, and the subroutine does nothing if the
!>                   content points to NULL
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_release(com_tasks, error)
    TYPE(fb_com_tasks_obj), INTENT(INOUT)    :: com_tasks
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(com_tasks%obj)) THEN
       CPPreconditionNoFail(com_tasks%obj%ref_count>0, cp_failure_level, routineP, error)
       com_tasks%obj%ref_count = com_tasks%obj%ref_count - 1
       IF (com_tasks%obj%ref_count == 0) THEN
          com_tasks%obj%ref_count = 1
          IF (ASSOCIATED(com_tasks%obj%tasks)) THEN
             DEALLOCATE(com_tasks%obj%tasks, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          com_tasks%obj%ref_count = 0
          DEALLOCATE(com_tasks%obj, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    ELSE
       NULLIFY(com_tasks%obj)
    END IF
  END SUBROUTINE fb_com_tasks_release


! **********************************************************************
!> \brief Releases an fb_com_atom_pairs object
!> \param atom_pairs the fb_com_atom_pairs object, its content must not
!>                    be UNDEFINED, and the subroutine does nothing if
!>                    the content points to NULL
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_atom_pairs_release(atom_pairs, error)
    TYPE(fb_com_atom_pairs_obj), &
      INTENT(INOUT)                          :: atom_pairs
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_atom_pairs_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(atom_pairs%obj)) THEN
       CPPreconditionNoFail(atom_pairs%obj%ref_count>0, cp_failure_level, routineP, error)
       atom_pairs%obj%ref_count = atom_pairs%obj%ref_count - 1
       IF (atom_pairs%obj%ref_count == 0) THEN
          atom_pairs%obj%ref_count = 1
          IF (ASSOCIATED(atom_pairs%obj%pairs)) THEN
             DEALLOCATE(atom_pairs%obj%pairs, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          atom_pairs%obj%ref_count = 0
          DEALLOCATE(atom_pairs%obj, stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
    ELSE
       NULLIFY(atom_pairs%obj)
    END IF
  END SUBROUTINE fb_com_atom_pairs_release


! **********************************************************************
!> \brief Nullifies a fb_com_tasks object, note that it does not release
!>        the original object. This procedure is used to nullify the
!>        pointer contained in the object which is used to associate to
!>        the actual object content
!> \param com_tasks the com_tasks object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_nullify(com_tasks)
    TYPE(fb_com_tasks_obj), INTENT(INOUT)    :: com_tasks

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(com_tasks%obj)
  END SUBROUTINE fb_com_tasks_nullify


! **********************************************************************
!> \brief Nullifies a fb_com_atom_pairs object, note that it does not
!>        release the original object. This procedure is used to nullify
!>        the pointer contained in the object which is used to associate
!>        to the actual object content
!> \param atom_pairs the fb_com_atom_pairs object
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_atom_pairs_nullify(atom_pairs)
    TYPE(fb_com_atom_pairs_obj), &
      INTENT(INOUT)                          :: atom_pairs

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_atom_pairs_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(atom_pairs%obj)
  END SUBROUTINE fb_com_atom_pairs_nullify


! **********************************************************************
!> \brief Associates one fb_com_tasks object to another
!> \param a the fb_com_tasks object to be associated
!> \param b the fb_com_tasks object that a is to be associated to
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_associate(a, b)
    TYPE(fb_com_tasks_obj), INTENT(OUT)      :: a
    TYPE(fb_com_tasks_obj), INTENT(IN)       :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_com_tasks_associate


! **********************************************************************
!> \brief Associates one fb_com_atom_pairs object to another
!> \param a the fb_com_atom_pairs object to be associated
!> \param b the fb_com_atom_pairs object that a is to be associated to
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_atom_pairs_associate(a, b)
    TYPE(fb_com_atom_pairs_obj), INTENT(OUT) :: a
    TYPE(fb_com_atom_pairs_obj), INTENT(IN)  :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_atom_pairs_associate', &
      routineP = moduleN//':'//routineN

    a%obj => b%obj
  END SUBROUTINE fb_com_atom_pairs_associate


! **********************************************************************
!> \brief Checks if a fb_com_tasks object is associated with an actual
!>        data content or not
!> \param com_tasks the fb_com_tasks object
!> \retval res ...
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  FUNCTION fb_com_tasks_has_data(com_tasks) RESULT(res)
    TYPE(fb_com_tasks_obj), INTENT(IN)       :: com_tasks
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_has_data', &
      routineP = moduleN//':'//routineN

    res = ASSOCIATED(com_tasks%obj)
  END FUNCTION fb_com_tasks_has_data


! **********************************************************************
!> \brief Checks if a fb_com_atom_pairs object is associated with an actual
!>        data content or not
!> \param atom_pairs the fb_com_atom_pairs object
!> \retval res ...
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  FUNCTION fb_com_atom_pairs_has_data(atom_pairs) RESULT(res)
    TYPE(fb_com_atom_pairs_obj), INTENT(IN)  :: atom_pairs
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_atom_pairs_has_data', &
      routineP = moduleN//':'//routineN

    res = ASSOCIATED(atom_pairs%obj)
  END FUNCTION fb_com_atom_pairs_has_data


! **********************************************************************
!> \brief Creates and initialises an empty fb_com_tasks object
!> \param com_tasks the fb_com_tasks object, its content must be NULL
!>                   and cannot be UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_create(com_tasks, error)
    TYPE(fb_com_tasks_obj), INTENT(INOUT)    :: com_tasks
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(com_tasks%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(com_tasks%obj, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       com_tasks%obj%task_dim = TASK_N_RECORDS
       com_tasks%obj%ntasks = 0
       com_tasks%obj%nencode = 0
       NULLIFY(com_tasks%obj%tasks)
       com_tasks%obj%ref_count = 1
       com_tasks%obj%id_nr = last_fb_com_tasks_id + 1
       last_fb_com_tasks_id = com_tasks%obj%id_nr
    END IF
  END SUBROUTINE fb_com_tasks_create


! **********************************************************************
!> \brief Creates and initialises an empty fb_com_atom_pairs object
!> \param atom_pairs the fb_com_atom_pairs object, its content must be
!>                    NULL and cannot be UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_atom_pairs_create(atom_pairs, error)
    TYPE(fb_com_atom_pairs_obj), &
      INTENT(INOUT)                          :: atom_pairs
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_atom_pairs_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(atom_pairs%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       ALLOCATE(atom_pairs%obj, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       atom_pairs%obj%npairs = 0
       atom_pairs%obj%natoms_encode = 0
       NULLIFY(atom_pairs%obj%pairs)
       atom_pairs%obj%ref_count = 1
       atom_pairs%obj%id_nr = last_fb_com_atom_pairs_id + 1
       last_fb_com_atom_pairs_id = atom_pairs%obj%id_nr
    END IF
  END SUBROUTINE fb_com_atom_pairs_create


! **********************************************************************
!> \brief Initialises an fb_com_tasks object, and makes it empty
!> \param com_tasks the fb_com_tasks object, its content must not be
!>                   NULL or UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_init(com_tasks, error)
    TYPE(fb_com_tasks_obj), INTENT(INOUT)    :: com_tasks
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(com_tasks%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(com_tasks%obj%tasks)) THEN
          DEALLOCATE(com_tasks%obj%tasks, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          NULLIFY(com_tasks%obj%tasks)
       END IF
       com_tasks%obj%task_dim = TASK_N_RECORDS
       com_tasks%obj%ntasks = 0
       com_tasks%obj%nencode = 0
    END IF
  END SUBROUTINE fb_com_tasks_init


! **********************************************************************
!> \brief Initialises an fb_com_atom_pairs object, and makes it empty
!> \param atom_pairs the fb_com_atom_pairs object, its content must not
!>                    be NULL or UNDEFINED
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_atom_pairs_init(atom_pairs, error)
    TYPE(fb_com_atom_pairs_obj), &
      INTENT(INOUT)                          :: atom_pairs
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_atom_pairs_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atom_pairs%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (ASSOCIATED(atom_pairs%obj%pairs)) THEN
          DEALLOCATE(atom_pairs%obj%pairs, STAT=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          NULLIFY(atom_pairs%obj%pairs)
       END IF
       atom_pairs%obj%npairs = 0
       atom_pairs%obj%natoms_encode = 0
    END IF
  END SUBROUTINE fb_com_atom_pairs_init


! **********************************************************************
!> \brief Gets attributes from a fb_com_tasks object, one should only
!>        access the data content in a fb_com_tasks object outside this
!>        module via this procedure.
!> \param com_tasks the fb_com_tasks object, its content must not be
!>                   NULL or UNDEFINED
!> \param task_dim [OPTIONAL]: if present, outputs com_tasks%obj%task_dim
!> \param ntasks [OPTIONAL]: if present, outputs com_tasks%obj%ntasks
!> \param nencode [OPTIONAL]: if present, outputs com_tasks%obj%nencode
!> \param tasks [OPTIONAL]: if present, outputs pointer com_tasks%obj%tasks
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_get(com_tasks, &
                              task_dim, &
                              ntasks, &
                              nencode, &
                              tasks, &
                              error)
    TYPE(fb_com_tasks_obj), INTENT(IN)       :: com_tasks
    INTEGER, INTENT(OUT), OPTIONAL           :: task_dim, ntasks, nencode
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: tasks
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(com_tasks%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(task_dim)) task_dim = com_tasks%obj%task_dim
       IF (PRESENT(ntasks)) ntasks = com_tasks%obj%ntasks
       IF (PRESENT(nencode)) nencode = com_tasks%obj%nencode
       IF (PRESENT(tasks)) tasks => com_tasks%obj%tasks
    END IF
  END SUBROUTINE fb_com_tasks_get


! **********************************************************************
!> \brief Gets attributes from a fb_com_atom_pairs object, one should
!>        only access the data content in a fb_com_atom_pairs object
!>        outside this module via this procedure.
!> \param atom_pairs the fb_com_atom_pairs object, its content must not
!>                    be NULL or UNDEFINED
!> \param npairs [OPTIONAL]: if present, outputs atom_pairs%obj%npairs
!> \param natoms_encode [OPTIONAL]: if present, outputs atom_pairs%obj%natoms_encode
!> \param pairs [OPTIONAL]: if present, outputs pointer atom_pairs%obj%pairs
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_atom_pairs_get(atom_pairs, &
                                   npairs, &
                                   natoms_encode, &
                                   pairs, &
                                   error)
    TYPE(fb_com_atom_pairs_obj), INTENT(IN)  :: atom_pairs
    INTEGER, INTENT(OUT), OPTIONAL           :: npairs, natoms_encode
    INTEGER(KIND=int_8), DIMENSION(:), &
      OPTIONAL, POINTER                      :: pairs
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_atom_pairs_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atom_pairs%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(npairs)) npairs = atom_pairs%obj%npairs
       IF (PRESENT(natoms_encode)) natoms_encode = atom_pairs%obj%natoms_encode
       IF (PRESENT(pairs)) pairs => atom_pairs%obj%pairs
    END IF
  END SUBROUTINE fb_com_atom_pairs_get


! **********************************************************************
!> \brief Sets attributes in a fb_com_tasks object, one should only
!>        access the data content in a fb_com_tasks object outside this
!>        module via this procedure.
!> \param com_tasks the fb_com_tasks object, its content must not be
!>                   NULL or UNDEFINED
!> \param task_dim [OPTIONAL]: if present, sets com_tasks%obj%task_dim
!> \param ntasks [OPTIONAL]: if present, sets com_tasks%obj%ntasks
!> \param nencode [OPTIONAL]: if present, sets com_tasks%obj%nencode
!> \param tasks [OPTIONAL]: if present, associates pointer com_tasks%obj%tasks
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_set(com_tasks, &
                              task_dim, &
                              ntasks, &
                              nencode, &
                              tasks, &
                              error)
    TYPE(fb_com_tasks_obj), INTENT(INOUT)    :: com_tasks
    INTEGER, INTENT(IN), OPTIONAL            :: task_dim, ntasks, nencode
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: tasks
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_tasks_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(com_tasks%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(task_dim)) com_tasks%obj%task_dim = task_dim
       IF (PRESENT(ntasks)) com_tasks%obj%ntasks = ntasks
       IF (PRESENT(nencode)) com_tasks%obj%nencode = nencode
       IF (PRESENT(tasks)) THEN
          IF (ASSOCIATED(com_tasks%obj%tasks)) THEN
             DEALLOCATE(com_tasks%obj%tasks, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          com_tasks%obj%tasks => tasks
       END IF
    END IF
  END SUBROUTINE fb_com_tasks_set


! **********************************************************************
!> \brief Sets attributes in a fb_com_atom_pairs object, one should only
!>        access the data content in a fb_com_atom_pairs object outside
!>        this module via this procedure.
!> \param atom_pairs the fb_com_atom_pairs object, its content must not
!>                    be NULL or UNDEFINED
!> \param npairs [OPTIONAL]: if present, sets atom_pairs%obj%npairs
!> \param natoms_encode [OPTIONAL]: if present, sets atom_pairs%obj%natoms_encode
!> \param pairs [OPTIONAL]: if present, associates pointer atom_pairs%obj%pairs
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_atom_pairs_set(atom_pairs, &
                                   npairs, &
                                   natoms_encode, &
                                   pairs, &
                                   error)
    TYPE(fb_com_atom_pairs_obj), &
      INTENT(INOUT)                          :: atom_pairs
    INTEGER, INTENT(IN), OPTIONAL            :: npairs, natoms_encode
    INTEGER(KIND=int_8), DIMENSION(:), &
      OPTIONAL, POINTER                      :: pairs
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_com_atom_pairs_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(atom_pairs%obj), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(npairs)) atom_pairs%obj%npairs = npairs
       IF (PRESENT(natoms_encode)) atom_pairs%obj%natoms_encode = natoms_encode
       IF (PRESENT(pairs)) THEN
          IF (ASSOCIATED(atom_pairs%obj%pairs)) THEN
             DEALLOCATE(atom_pairs%obj%pairs, STAT=stat)
             CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          END IF
          atom_pairs%obj%pairs => pairs
       END IF
    END IF
  END SUBROUTINE fb_com_atom_pairs_set


! **********************************************************************
!> \brief Start from a local set of tasks that has desc/src process equal
!>        to the local MPI rank, communicate with other processes so
!>        that a new local set of tasks is constructed with src/desc
!>        process equal to the local MPI rank
!> \param tasks_dest_is_me the local com_task object with all tasks
!>                          having the desc process id equal to my_id
!> \param direction direction of operation:
!>                   ">" means from tasks_dest_is_me construct tasks_src_is_me
!>                   "<" means from tasks_src_is_me construct tasks_dest_is_me
!> \param tasks_src_is_me the local com_task object with all tasks
!>                          having the src process id equal to my_id
!> \param para_env CP2K parallel environment object that stores MPI related
!>                  information of the current run
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_transpose_dest_src(tasks_dest_is_me, &
                                             direction, &
                                             tasks_src_is_me, &
                                             para_env, &
                                             error)
    TYPE(fb_com_tasks_obj), INTENT(INOUT)    :: tasks_dest_is_me
    CHARACTER, INTENT(IN)                    :: direction
    TYPE(fb_com_tasks_obj), INTENT(INOUT)    :: tasks_src_is_me
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'fb_com_tasks_transpose_dest_src', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ii, ind, ipe, itask, jj, LT_UNIT, nencode, ntasks_in, &
      ntasks_out, rank, rank_pos, stat, task_dim
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks_in, tasks_out
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: recv_buf, recv_disps, &
                                                recv_sizes, send_buf, &
                                                send_disps, send_sizes
    LOGICAL                                  :: failure

!BEG:DEBUG:LT:2015/01/29
!END:DEBUG:LT:2015/01/29
!BEG:DEBUG:LT:2015/01/30

    LT_UNIT = 100 + para_env%mepos
!END:DEBUG:LT:2015/01/30

    CALL timeset(routineN, handle)

    failure = .FALSE.

    NULLIFY(tasks_in, tasks_out)

    IF (direction == "<") THEN
       CALL fb_com_tasks_get(com_tasks=tasks_src_is_me, &
                             task_dim=task_dim, &
                             ntasks=ntasks_in, &
                             tasks=tasks_in, &
                             nencode=nencode, &
                             error=error)
       rank_pos = TASK_DEST
    ELSE
       CALL fb_com_tasks_get(com_tasks=tasks_dest_is_me, &
                             task_dim=task_dim, &
                             ntasks=ntasks_in, &
                             tasks=tasks_in, &
                             nencode=nencode, &
                             error=error)
       rank_pos = TASK_SRC
    END IF

    ! allocate local arrays
    ALLOCATE(send_sizes(para_env%num_pe), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(send_disps(para_env%num_pe), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(send_buf(para_env%num_pe), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ALLOCATE(recv_sizes(para_env%num_pe), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(recv_disps(para_env%num_pe), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(recv_buf(para_env%num_pe), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

!BEG:DEBUG:LT:2015/01/20
    ! write (LT_UNIT,*) "size(send_buf) = ", SIZE(send_buf)
    ! write (LT_UNIT,*) "ntasks_in = ", ntasks_in
    ! write (LT_UNIT,*) "size(tasks_in) = ", SIZE(tasks_in)
!END:DEBUG:LT:2015/01/20

    ! first count how many local recv/send tasks need to be sent to
    ! other processes, and share this information with the other
    ! processes.  using send_buf as a temporary array for counting
    send_buf = 0
    ! looping over local task list
    DO itask = 1, ntasks_in
       rank = INT(tasks_in(rank_pos,itask)) + 1

!BEG:DEBUG:LT:2015/01/20
       ! write (LT_UNIT,*) "itask, rank", itask, rank
!END:DEBUG:LT:2015/01/20

       send_buf(rank) = send_buf(rank) + 1
    END DO

!BEG:DEBUG:LT:2015/01/20
    ! write (LT_UNIT,*) "send_buf = ", send_buf
!END:DEBUG:LT:2015/01/20

    CALL mp_alltoall(send_buf, recv_buf, 1, para_env%group)

    ! now that we know how many recv/send tasks to send, pack the
    ! tasks, and send them around, so that the recv/send tasks are
    ! sent to the correct src/dest processes, and these then are
    ! collected into the send/recv tasks list on each of the src/dest
    ! processes

    send_sizes = 0
    send_disps = 0
    recv_sizes = 0
    recv_disps = 0

    ! work out the sizes of send and recv buffers and allocate them
    send_sizes(1) = send_buf(1) * task_dim
    recv_sizes(1) = recv_buf(1) * task_dim
    DO ipe = 2, para_env%num_pe
       send_sizes(ipe) = send_buf(ipe) * task_dim
       send_disps(ipe) = send_disps(ipe-1) + send_sizes(ipe-1)
       recv_sizes(ipe) = recv_buf(ipe) * task_dim
       recv_disps(ipe) = recv_disps(ipe-1) + recv_sizes(ipe-1)
    END DO

    ! reallocate send and recv buffers to the correct sizes for
    ! transferring the actual tasks
    DEALLOCATE(send_buf, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_buf, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(send_buf(SUM(send_sizes)), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(recv_buf(SUM(recv_sizes)), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! now that the send buffer is of correct size, do packing
    ! send_buf and recv_buf may be zero sized
    IF (SIZE(send_buf) > 0) send_buf = 0
    IF (SIZE(recv_buf) > 0) recv_buf = 0
    send_sizes = 0
    DO itask = 1, ntasks_in
       rank = INT(tasks_in(rank_pos,itask)) + 1
       DO ii = 1, task_dim
          ind = send_disps(rank) + send_sizes(rank) + ii
          send_buf(ind) = INT(tasks_in(ii,itask))
       END DO
       send_sizes(rank) = send_sizes(rank) + task_dim
    END DO
    ! do communication
    CALL mp_alltoall(send_buf, send_sizes, send_disps, &
                     recv_buf, recv_sizes, recv_disps, &
                     para_env%group)

    ! deallocate send buffers
    DEALLOCATE(send_buf, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(send_sizes, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(send_disps, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! allocate the output task list
    ntasks_out = SUM(recv_sizes) / task_dim
    ! this will not be deallocated in this subroutine
    ALLOCATE(tasks_out(task_dim,ntasks_out), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! do unpacking
    itask = 0
    DO ipe = 1, para_env%num_pe
       DO ii = 0, recv_sizes(ipe) / task_dim - 1
          itask = itask + 1
          DO jj = 1, task_dim
             ind = recv_disps(ipe) + ii * task_dim + jj
             tasks_out(jj,itask) = recv_buf(ind)
          END DO
       END DO
    END DO

    ! set output tasks
    IF (direction == "<") THEN
       CALL fb_com_tasks_set(com_tasks=tasks_dest_is_me, &
                             task_dim=task_dim, &
                             ntasks=ntasks_out, &
                             tasks=tasks_out, &
                             nencode=nencode, &
                             error=error)
    ELSE
       CALL fb_com_tasks_set(com_tasks=tasks_src_is_me, &
                             task_dim=task_dim, &
                             ntasks=ntasks_out, &
                             tasks=tasks_out, &
                             nencode=nencode, &
                             error=error)
    END IF

    ! deallocate recv buffers
    DEALLOCATE(recv_buf, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_sizes, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_disps, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    CALL timestop(handle)

  END SUBROUTINE fb_com_tasks_transpose_dest_src


! **********************************************************************
!> \brief Generate send or receive atom_pair lists from a com_tasks
!>        object. atom_pair list is used as a condensed index for the
!>        local/remote matrix blocks to be sent/received.
!> \param com_tasks the com_tasks object
!> \param atom_pairs fb_com_atom_pairs_obj containing  list of encoded
!>                    atomic pair indices and the dest/src proc id for
!>                    the matrix block to be sent/received.
!> \param natoms_encode the total number of atoms the atomic pair indices
!>                       corresponds to, and it is used for encode the
!>                       atom_pairs values
!> \param send_or_recv whether the atom_pair to be generated is for
!>                      the local matrix blocks to be sent or the
!>                      remote matrix blocks to be received for this MPI
!>                      process
!> \param symmetric whether the matrix the blocks are to be taken from
!>                   is symmetric
!> \param error CP2K data container for error logging
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_build_atom_pairs(com_tasks, &
                                           atom_pairs, &
                                           natoms_encode, &
                                           send_or_recv, &
                                           symmetric, &
                                           error)
    TYPE(fb_com_tasks_obj), INTENT(IN)       :: com_tasks
    TYPE(fb_com_atom_pairs_obj), &
      INTENT(INOUT)                          :: atom_pairs
    INTEGER, INTENT(IN)                      :: natoms_encode
    CHARACTER(len=*), INTENT(IN)             :: send_or_recv
    LOGICAL, INTENT(IN), OPTIONAL            :: symmetric
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'fb_com_tasks_build_atom_pairs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, iatom, &
                                                ii, itask, jatom, LT_UNIT, &
                                                npairs, ntasks, rank, &
                                                rank_pos, stat
    INTEGER(KIND=int_8)                      :: pair
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: pairs
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_index
    LOGICAL                                  :: check_ok, failure, &
                                                my_symmetric

!BEG:DEBUG:LT:2015/01/30
!END:DEBUG:LT:2015/01/30
!BEG:DEBUG:LT:2015/01/30

    LT_UNIT = 200
!END:DEBUG:LT:2015/01/30

    CALL timeset(routineN, handle)

    failure = .FALSE.

    NULLIFY(pairs, tasks)

    check_ok = fb_com_atom_pairs_has_data(atom_pairs)
    CPPrecondition(check_ok, cp_failure_level, routineP, error, failure)

    ! initialise atom_pairs
    CALL fb_com_atom_pairs_init(atom_pairs, error)

    IF (TRIM(send_or_recv) == "send") THEN
       rank_pos = TASK_DEST
    ELSE
       rank_pos = TASK_SRC
    END IF

    IF (PRESENT(symmetric)) THEN
       my_symmetric = symmetric
    ELSE
       my_symmetric = .FALSE.
    END IF

    CALL fb_com_tasks_get(com_tasks=com_tasks, &
                          ntasks=ntasks, &
                          tasks=tasks, &
                          error=error)

    ALLOCATE(pairs(ntasks), STAT=stat)
    CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
    ! we can have cases where ntasks == 0
    IF (SIZE(pairs) > 0) pairs = 0
    npairs = ntasks

    DO itask = 1, ntasks
       pair = tasks(TASK_PAIR,itask)
       CALL fb_com_tasks_decode_pair(pair, iatom, jatom, natoms_encode)

       IF (my_symmetric) THEN
          IF (iatom .LE. jatom) THEN
             arow = iatom
             acol = jatom
          ELSE
             arow = jatom
             acol = iatom
          END IF
       ELSE
          arow = iatom
          acol = jatom
       END IF
       rank = INT(tasks(rank_pos,itask))
       CALL fb_com_atom_pairs_encode(pairs(itask), &
                                     rank, arow, acol, natoms_encode)

!BEG:DEBUG:LT:2015/01/25
       ! write (LT_UNIT,*) &
       !       "|LT| fb_com_tasks_build_atom_pairs: itask = ", itask
       ! write (LT_UNIT,*) &
       !       "|LT| fb_com_tasks_build_atom_pairs: rank, arow, acol, natoms_encode = ", &
       !       rank, arow, acol, natoms_encode
       ! write (LT_UNIT,*) &
       !       "|LT| fb_com_tasks_build_atom_pairs: pairs(itask) = ", pairs(itask)
!END:DEBUG:LT:2015/01/25

    END DO

!BEG:DEBUG:LT:2015/01/25
    ! write (LT_UNIT,*) "|LT| fb_com_tasks_build_atom_pairs: SIZE(pairs) = ", SIZE(pairs)
    ! write (LT_UNIT,*) "|LT| fb_com_tasks_build_atom_pairs: pairs = ", pairs
!END:DEBUG:LT:2015/01/25

    ! sort atom_pairs so that the pairs are ordered process blocks and
    ! that possible duplicates may be found (we don't want to send or
    ! receive same information to the same destination or source more
    ! than once)
    IF (npairs > 0) THEN
       ALLOCATE(tmp_index(npairs), STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
       ! only sort the actual pairs recorded in the send list
       CALL sort(pairs, npairs, tmp_index)
       DEALLOCATE(tmp_index, STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF

!BEG:DEBUG:LT:2015/01/25
    ! write (LT_UNIT,*) "|LT| fb_com_tasks_build_atom_pairs: After sort"
    ! write (LT_UNIT,*) "|LT| fb_com_tasks_build_atom_pairs: SIZE(pairs) = ", SIZE(pairs)
    ! write (LT_UNIT,*) "|LT| fb_com_tasks_build_atom_pairs: pairs = ", pairs
!END:DEBUG:LT:2015/01/25

    ! remove duplicates
    IF (npairs > 1) THEN
       npairs = 1
       ! first atom pair must be allowed
       DO ii = 2, ntasks
          IF (pairs(ii) > pairs(ii-1)) THEN
             npairs = npairs + 1
             pairs(npairs) = pairs(ii)
          END IF
       END DO
       ! reallocate the pairs list
       CALL reallocate(pairs, 1, npairs)
    END IF

    CALL fb_com_atom_pairs_set(atom_pairs=atom_pairs, &
                               pairs=pairs, &
                               npairs=npairs, &
                               natoms_encode=natoms_encode, &
                               error=error)

!BEG:DEBUG:LT:2015/01/25
    ! write (LT_UNIT,*) &
    !       "|LT| fb_com_tasks_build_atom_pairs: After fb_com_atom_pairs_set"
    ! CALL fb_com_atom_pairs_get(atom_pairs=atom_pairs, &
    !                            pairs=pairs, &
    !                            npairs=npairs, &
    !                            error=error)
    ! write (LT_UNIT,*) "|LT| fb_com_tasks_build_atom_pairs: npairs = ", npairs
    ! write (LT_UNIT,*) "|LT| fb_com_tasks_build_atom_pairs: pairs = ", pairs
!END:DEBUG:LT:2015/01/25

    CALL timestop(handle)

  END SUBROUTINE fb_com_tasks_build_atom_pairs


! **********************************************************************
!> \brief Encodes (iatom, jatom) pair index of a block into a single
!>        integer
!> \param ind encoded integer
!> \param iatom the first index of the (iatom, jatom) block index
!> \param jatom the second index of the (iatom, jatom) block index
!> \param natoms the total number of atoms iatom and jatom indexes
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_encode_pair(ind, iatom, jatom, natoms)
    INTEGER(KIND=int_8), INTENT(OUT)         :: ind
    INTEGER, INTENT(IN)                      :: iatom, jatom, natoms

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_com_tasks_encode_pair', &
      routineP = moduleN//':'//routineN

    INTEGER(KIND=int_8)                      :: iatom8, jatom8, natoms8

    natoms8 = INT(natoms, int_8)
    iatom8 = INT(iatom, int_8)
    jatom8 = INT(jatom, int_8)

    ind = (iatom8 - 1_int_8) * natoms8 + (jatom8 - 1_int_8)
  END SUBROUTINE fb_com_tasks_encode_pair


! **********************************************************************
!> \brief Dncodes a single integer into (iatom, jatom) pair index of
!>        a block into a single
!> \param ind encoded integer
!> \param iatom the first index of the (iatom, jatom) block index
!> \param jatom the second index of the (iatom, jatom) block index
!> \param natoms the total number of atoms iatom and jatom indexes
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_tasks_decode_pair(ind, iatom, jatom, natoms)
    INTEGER(KIND=int_8), INTENT(IN)          :: ind
    INTEGER, INTENT(OUT)                     :: iatom, jatom
    INTEGER, INTENT(IN)                      :: natoms

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_com_tasks_decode_pair', &
      routineP = moduleN//':'//routineN

    INTEGER(KIND=int_8)                      :: iatom8, jatom8, natoms8

    natoms8 = INT(natoms, int_8)
    iatom8 = ind / natoms8 + 1_int_8
    jatom8 = MOD(ind, natoms8) + 1_int_8
    iatom = INT(iatom8, int_4)
    jatom = INT(jatom8, int_4)
  END SUBROUTINE fb_com_tasks_decode_pair


! **********************************************************************
!> \brief Encodes (rank, iatom, jatom) index of a communication task---to
!>         send/receive a block to/from a process---into a single integer
!> \param ind encoded integer
!> \param pe the rank of the process the block to be send to or receive
!>            from
!> \param iatom the first index of the (iatom, jatom) block index
!> \param jatom the second index of the (iatom, jatom) block index
!> \param natoms the total number of atoms iatom and jatom indexes
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_atom_pairs_encode(ind, pe, iatom, jatom, natoms)
    INTEGER(KIND=int_8), INTENT(OUT)         :: ind
    INTEGER, INTENT(IN)                      :: pe, iatom, jatom, natoms

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_com_atom_pairs_encode', &
      routineP = moduleN//':'//routineN

    INTEGER(KIND=int_8)                      :: natoms8, pair

! pe must start count from 0 (i.e same as MPI convension)

    natoms8 = INT(natoms, int_8)
    CALL fb_com_tasks_encode_pair(pair, iatom, jatom, natoms)
    ind = INT(pe, int_8) * natoms8 * natoms8 + pair
  END SUBROUTINE fb_com_atom_pairs_encode


! **********************************************************************
!> \brief Decodes a single integer into the (rank, iatom, jatom) index
!>        of a communication task to send/receive a block to/from a
!>        process
!> \param ind    : encoded integer
!> \param pe     : the rank of the process the block to be send to or receive
!>            from
!> \param iatom  : the first index of the (iatom, jatom) block index
!> \param jatom  : the second index of the (iatom, jatom) block index
!> \param natoms : the total number of atoms iatom and jatom indexes
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_atom_pairs_decode(ind, pe, iatom, jatom, natoms)
    INTEGER(KIND=int_8), INTENT(IN)          :: ind
    INTEGER, INTENT(OUT)                     :: pe, iatom, jatom
    INTEGER, INTENT(IN)                      :: natoms

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_com_atom_pairs_decode', &
      routineP = moduleN//':'//routineN

    INTEGER(KIND=int_8)                      :: natoms8, pair

! pe start count from 0 (i.e same as MPI convension)

    natoms8 = INT(natoms, int_8)
    pe = INT(ind / (natoms8*natoms8), int_4)
    pair =  MOD(ind, natoms8*natoms8)
    CALL fb_com_tasks_decode_pair(pair, iatom, jatom, natoms)
  END SUBROUTINE fb_com_atom_pairs_decode


! **********************************************************************
!> \brief Calculate the MPI send or recv buffer sizes according to the
!>        communication pairs (atom_pairs) and DBCSR matrix data.
!>        Each atom_pair corresponds to one DBCSR matrix block that
!>        needs to be sent or recerived.
!> \param atom_pairs : the communication pair object for either sending
!>                     or receiving
!> \param dbcsr_matrix : the DBCSR matrix the blocks originate from
!>                       or has to be sent to
!> \param nprocs : total number of MPI processes in communicator
!> \param row_blk_sizes : row_blk_sizes(iblkrow) = number of element rows
!>                        in each block in the iblkrow-th block row of
!>                        the DBCSR matrix
!> \param col_blk_sizes : col_blk_sizes(iblkcol) = number of element cols
!>                        in each block in the iblkcol-th block col of
!>                        the DBCSR matrix
!> \param sendrecv_sizes : size required for the send of recv buffer
!>                         for each dest/src process
!> \param sendrecv_disps : sendrecv_disps(ipe) + 1 = starting location
!>                         in send/recv buffer for data destined for
!>                         process ipe
!> \param sendrecv_pair_counts : sendrecv_pair_counts(ipe) = number of
!>                               pairs (blocks) to be sent to or recv
!>                               from process ipe
!> \param sendrecv_pair_disps send_recv_pair_disps(ipe) + 1 = start
!>                               location in atom_pairs array for
!>                               all the pairs to be sent to or recv
!>                               from process ipe
!> \param error : cp2k error container
!> \param row_map : optional blk row map for the DBCSR blocks
!> \param col_map : optional blk col map for the DBCSR blocks
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_com_atom_pairs_calc_buffer_sizes(atom_pairs, &
                                                 dbcsr_matrix, &
                                                 nprocs, &
                                                 row_blk_sizes, &
                                                 col_blk_sizes, &
                                                 sendrecv_sizes, &
                                                 sendrecv_disps, &
                                                 sendrecv_pair_counts, &
                                                 sendrecv_pair_disps, &
                                                 error, &
                                                 row_map, &
                                                 col_map)
    TYPE(fb_com_atom_pairs_obj), INTENT(IN)  :: atom_pairs
    TYPE(cp_dbcsr_type), POINTER             :: dbcsr_matrix
    INTEGER, INTENT(IN)                      :: nprocs
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_blk_sizes, col_blk_sizes
    INTEGER, DIMENSION(:), INTENT(OUT)       :: sendrecv_sizes, &
                                                sendrecv_disps, &
                                                sendrecv_pair_counts, &
                                                sendrecv_pair_disps
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: row_map, col_map

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'fb_com_atom_pairs_calc_buffer_sizes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ipair, ipe, jatom, &
                                                natoms_encode, ncols_blk, &
                                                npairs, nrows_blk, pe
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: pairs
    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.

    NULLIFY(pairs)

    check_ok = SIZE(sendrecv_sizes) == nprocs .AND. &
               SIZE(sendrecv_disps) == nprocs .AND. &
               SIZE(sendrecv_pair_counts) == nprocs .AND. &
               SIZE(sendrecv_pair_disps) == nprocs
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)

    check_ok = fb_com_atom_pairs_has_data(atom_pairs)
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)

    CALL fb_com_atom_pairs_get(atom_pairs=atom_pairs, &
                               pairs=pairs, &
                               npairs=npairs, &
                               natoms_encode=natoms_encode, &
                               error=error)

    sendrecv_sizes = 0
    sendrecv_pair_counts = 0
    DO ipair = 1, npairs
       ! decode processor and (iatom, jatom) information
       CALL fb_com_atom_pairs_decode(pairs(ipair), &
                                     pe, iatom, jatom, natoms_encode)
       pe = pe + 1 ! we need proc to count from 1
       IF (PRESENT(row_map)) iatom = row_map(iatom)
       IF (PRESENT(col_map)) jatom = row_map(jatom)
       nrows_blk = row_blk_sizes(iatom)
       ncols_blk = col_blk_sizes(jatom)
       sendrecv_sizes(pe) = sendrecv_sizes(pe) + nrows_blk * ncols_blk
       sendrecv_pair_counts(pe) = sendrecv_pair_counts(pe) + 1
    END DO
    ! calculate displacements of the data of each destibation pe in
    ! send buffer and in the list of pairs to be sent
    sendrecv_disps = 0
    sendrecv_pair_disps = 0
    DO ipe = 2, nprocs
       sendrecv_disps(ipe) = sendrecv_disps(ipe-1) + sendrecv_sizes(ipe-1)
       sendrecv_pair_disps(ipe) = sendrecv_pair_disps(ipe-1) + sendrecv_pair_counts(ipe-1)
    END DO

  END SUBROUTINE fb_com_atom_pairs_calc_buffer_sizes


END MODULE qs_fb_com_tasks_types
