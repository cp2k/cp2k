!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief utilities to represent orbitals that are localized on few atoms
!> \par History
!>      07.2002 created [fawzi]
!>      11.2002 updated and checked [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_p_sparse_psi
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type
  USE kinds,                           ONLY: dp
  USE particle_methods,                ONLY: get_particle_set
  USE particle_types,                  ONLY: particle_type
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_p_sparse_psi'
  INTEGER, SAVE, PRIVATE :: last_proj_id_nr=0

!***

! *****************************************************************************
!> \brief represent a simple projection that involves orbitals from a small
!>      group of atoms
!> \param atoms array with the indexes of the atoms that form this group
!> \param proj_indexes the mapping between the indexes of the reduced
!>        vectors and the full vectors
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_p_projection_type
     INTEGER, DIMENSION(:), POINTER :: atoms
     INTEGER, DIMENSION(:), POINTER :: proj_indexes
     INTEGER :: ref_count, id_nr
  END TYPE qs_p_projection_type

! *****************************************************************************
!> \brief to be able to build array of pointers to projections
!> \param projection the projection in this object
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_p_projection_p_type
     TYPE(qs_p_projection_type), POINTER :: projection
  END TYPE qs_p_projection_p_type

CONTAINS

! *****************************************************************************
!> \brief allocates and initializes a qs_p_projection_type
!> \param p_proj the projection to initialize
!> \param atoms the atoms that are part of the projection group.
!> \param particle_set the particle set from where to extract the indexes
!>        of the ao
!> \param qs_kind_set ...
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE p_proj_create(p_proj, atoms, particle_set, qs_kind_set)
    TYPE(qs_p_projection_type), POINTER      :: p_proj
    INTEGER, DIMENSION(:), INTENT(in)        :: atoms
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CHARACTER(len=*), PARAMETER :: routineN = 'p_proj_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, nparticle, pos, &
                                                size_proj
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf

    ALLOCATE(p_proj)
    NULLIFY(p_proj%atoms, p_proj%proj_indexes)
    ALLOCATE(p_proj%atoms(SIZE(atoms)))

    p_proj%atoms=atoms
    nparticle = SIZE(particle_set)
    ALLOCATE (first_sgf(nparticle), last_sgf(nparticle))

    CALL get_particle_set(particle_set, qs_kind_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf)
    size_proj = 0
    DO i=1,SIZE(p_proj%atoms)
      size_proj = size_proj + last_sgf(p_proj%atoms(i)) -&
                  first_sgf(p_proj%atoms(i)) + 1
    END DO

    ALLOCATE(p_proj%proj_indexes(size_proj))
    pos = 1
    DO i=1,SIZE(p_proj%atoms)
      DO j=first_sgf(p_proj%atoms(i)),last_sgf(p_proj%atoms(i))
        p_proj%proj_indexes(pos) = j
        pos = pos + 1
      END DO
    END DO

    DEALLOCATE (first_sgf, last_sgf)

    last_proj_id_nr=last_proj_id_nr+1
    p_proj%id_nr=last_proj_id_nr
    p_proj%ref_count=1
END SUBROUTINE p_proj_create

! *****************************************************************************
!> \brief augments the retain count by one, to be called to hold a shared copy
!>      of this object
!> \param p_proj the projection to retain
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE p_proj_retain(p_proj)
    TYPE(qs_p_projection_type), POINTER      :: p_proj

    CHARACTER(len=*), PARAMETER :: routineN = 'p_proj_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(p_proj),cp_failure_level,routineP,failure)
  CPPrecondition(p_proj%ref_count>0,cp_failure_level,routineP,failure)
  p_proj%ref_count=p_proj%ref_count+1
END SUBROUTINE p_proj_retain

! *****************************************************************************
!> \brief decreases the retain count by one, deleting the objecs when it hits 0.
!>      To be called when you no longer need a shared copy you retained
!> \param p_proj the projection to release
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE p_proj_release(p_proj)
    TYPE(qs_p_projection_type), POINTER      :: p_proj

    CHARACTER(len=*), PARAMETER :: routineN = 'p_proj_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(p_proj)) THEN
     CPPrecondition(p_proj%ref_count>0,cp_failure_level,routineP,failure)
     p_proj%ref_count=p_proj%ref_count-1
     IF (p_proj%ref_count<1) THEN
        DEALLOCATE(p_proj%atoms,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,failure)
        DEALLOCATE(p_proj%proj_indexes,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,failure)
        p_proj%ref_count=0
        DEALLOCATE(p_proj, stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,failure)
     END IF
  END IF
END SUBROUTINE p_proj_release

! *****************************************************************************
!> \brief injects the columns of a blacs matrix of the minimal system
!>      into the full basis:
!>      full_m= i min_m
!> \param p_proj the projection that defines the injection
!> \param min_m the full matrix in the reduced basis
!> \param full_m the full matrix in the full basis that gets the new
!>        values
!> \param start_col_full the column of full_m at which the first column of
!>        (i min_m) will be written (defaults to 1) 
!> \param start_col_min ...
!> \param ncol number of columns to set (defaults to all the columns of
!>        min_m
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE p_proj_transfer_blacs_to_f(p_proj, min_m, full_m,&
     start_col_full, start_col_min, ncol)
    TYPE(qs_p_projection_type), POINTER      :: p_proj
    TYPE(cp_fm_type), POINTER                :: min_m, full_m
    INTEGER, INTENT(in), OPTIONAL            :: start_col_full, &
                                                start_col_min, ncol

    CHARACTER(len=*), PARAMETER :: routineN = 'p_proj_transfer_blacs_to_f', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: max_blocksize = 100

    INTEGER :: blocksize, handle, i, icol, my_ncol, my_start_col_full, &
      my_start_col_min, ncol_full, ncol_min, nrow_full, nrow_min, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: tmp_full, tmp_min

  CALL timeset(routineN,handle)
  failure=.FALSE.

  CPPrecondition(ASSOCIATED(p_proj),cp_failure_level,routineP,failure)
  CPPrecondition(p_proj%ref_count>0,cp_failure_level,routineP,failure)
  CALL cp_fm_get_info(min_m, nrow_global=nrow_min, ncol_global=ncol_min)
  CALL cp_fm_get_info(full_m, nrow_global=nrow_full, ncol_global=ncol_full)

  my_start_col_full=1
  IF (PRESENT(start_col_full)) my_start_col_full=start_col_full
  my_start_col_min=1
  IF (PRESENT(start_col_min)) my_start_col_min=start_col_min
  my_ncol=ncol_min-my_start_col_min
  IF (PRESENT(ncol)) my_ncol=ncol
  blocksize=MIN(my_ncol,max_blocksize)

  CPPrecondition(my_ncol>=0,cp_failure_level,routineP,failure)
  CPPrecondition(my_start_col_min+my_ncol<ncol_min+2,cp_failure_level,routineP,failure)
  CPPrecondition(my_start_col_full+my_ncol<ncol_full+2,cp_failure_level,routineP,failure)
  CPPrecondition(nrow_min>=SIZE(p_proj%proj_indexes),cp_failure_level,routineP,failure)

  ALLOCATE(tmp_min(blocksize,nrow_min), tmp_full(blocksize,nrow_full),&
       stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,failure)

  DO icol=0,my_ncol-1,blocksize
     IF (icol+blocksize>my_ncol) blocksize=my_ncol-icol
     CALL cp_fm_get_submatrix(min_m, target_m=tmp_min, &
          start_col=icol+my_start_col_min, n_cols=blocksize, &
          transpose=.TRUE.)

     CALL dcopy(SIZE(tmp_full,1)*SIZE(tmp_full,2),0.0_dp,0,tmp_full(1,1),1)
     DO i=1,SIZE(p_proj%proj_indexes)
        tmp_full(1:blocksize,p_proj%proj_indexes(i))=tmp_min(1:blocksize,i)
     END DO

     CALL cp_fm_set_submatrix(full_m, new_values=tmp_full, &
          start_col=icol+my_start_col_full, n_cols=blocksize, &
          transpose=.TRUE.)
  END DO

  DEALLOCATE(tmp_min, tmp_full, stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP)

  CALL timestop(handle)
END SUBROUTINE p_proj_transfer_blacs_to_f

! *****************************************************************************
!> \brief restrains the requested columns of the full matrix to the space
!>      of the minimal basis, i.e. computes:
!>      full_m(:,start_col:star_col+ncol) = i pi full_m(:,start_col:star_col+ncol)
!> \param p_proj ...
!> \param full_m the matrix in the full basis to be restrained
!> \param start_col the first column to restrain (defaults to 1)
!> \param ncol number of columns to restrain
!> \par History
!>      03.2003 created [fawzi]
!> \author fawzi
!> \note
!>      keep_row could easily reduced to (min(proj_indexes):max(proj_indexes))
! *****************************************************************************
SUBROUTINE p_proj_restrain_f(p_proj, full_m,start_col, ncol)
    TYPE(qs_p_projection_type), POINTER      :: p_proj
    TYPE(cp_fm_type), POINTER                :: full_m
    INTEGER, INTENT(in), OPTIONAL            :: start_col, ncol

    CHARACTER(len=*), PARAMETER :: routineN = 'p_proj_restrain_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, icol, irow, &
                                                max_index, min_index, &
                                                my_end_col, my_start_col, &
                                                row_index, stat
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: keep_row

  CALL timeset(routineN,handle)

  failure=.FALSE.
  NULLIFY(row_indices,col_indices)

  CPPrecondition(ASSOCIATED(p_proj),cp_failure_level,routineP,failure)
  CPPrecondition(p_proj%ref_count>0,cp_failure_level,routineP,failure)
  CPPrecondition(ASSOCIATED(full_m),cp_failure_level,routineP,failure)
  CPPrecondition(full_m%ref_count>0,cp_failure_level,routineP,failure)
  min_index=MINVAL(p_proj%proj_indexes)
  max_index=MAXVAL(p_proj%proj_indexes)
  ALLOCATE(keep_row(min_index:max_index),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,failure)

  my_start_col=1
  IF (PRESENT(start_col)) my_start_col=start_col
  CALL cp_fm_get_info(full_m,ncol_global=my_end_col, &
       row_indices=row_indices, col_indices=col_indices)
  IF (PRESENT(ncol)) my_end_col=start_col+ncol

  keep_row=.FALSE.
  DO irow=1,SIZE(p_proj%proj_indexes)
     keep_row(p_proj%proj_indexes(irow))=.TRUE.
  END DO

  DO icol=1,SIZE(full_m%local_data,2)
     IF (col_indices(icol)>start_col .AND.&
          col_indices(icol)<my_end_col) THEN
        DO irow=1,SIZE(full_m%local_data,1)
           row_index=row_indices(irow)
           IF (row_index<min_index.OR.row_index>max_index) THEN
              full_m%local_data(irow,icol)=0.0_dp
           ELSE IF (.NOT.keep_row(row_index)) THEN
              full_m%local_data(irow,irow)=0.0_dp
           END IF
        END DO
     END IF
  END DO

  DEALLOCATE(keep_row,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP)

  CALL timestop(handle)
END SUBROUTINE p_proj_restrain_f

END MODULE qs_p_sparse_psi
