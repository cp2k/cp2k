!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****s* cp2k/realspace_task_selection [1.0] *
!!
!!   NAME
!!     realspace_task_selection
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (14-May-2007)
!!
!!   MODIFICATION HISTORY
!!     JGH (14-May-2007) : split off from realspace_grid_types
!!
!!   NOTES
!!     Organizes the task selection for parallel/distributed grids
!!
!*****
!******************************************************************************

MODULE realspace_task_selection
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_cart_sub,&
                                             mp_comm_free,&
                                             mp_environ,&
                                             mp_max,&
                                             mp_shift
  USE realspace_grid_types,            ONLY: realspace_grid_type
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             get_block_node,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rs_get_my_tasks,&
            rs_get_loop_vars

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'realspace_task_selection'

!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!***************************************************************************
!!****f* realspace_grid_types/rs_get_my_tasks [1.0] *
!!
!!   NAME
!!     rs_get_my_tasks
!!
!!   FUNCTION
!!     Assembles tasks to be performed on local grid
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rs: the grid
!!     - tasks: the task set generate on this processor
!!     - tasks_local: the task set to be processed localy
!!     - npme: Number of tasks for local processing
!!
!!   AUTHOR
!!     JGH (19.06.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************

  SUBROUTINE rs_get_my_tasks ( rs, tasks, npme, ival, rval, &
       pmat, pmat2, pmat3, pcor ,symmetric, error)

    TYPE(realspace_grid_type), POINTER       :: rs
    INTEGER, DIMENSION(:, :), POINTER        :: tasks
    INTEGER, INTENT(OUT)                     :: npme
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: ival
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: rval
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: pmat, pmat2, pmat3
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: pcor
    LOGICAL, INTENT(IN), OPTIONAL            :: symmetric
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER :: acol, arow, dir, handle, i, ic, icmax, ileft, isend, j, k, n1, &
      npmax, nppp, plength, plength2, plength3, stat, subgroup, subpos, &
      subsize
    INTEGER, DIMENSION(3)                    :: lb, ub
    INTEGER, DIMENSION(:, :), POINTER        :: ilist, plist, pls, pls2, &
                                                pls3, tlist
    LOGICAL                                  :: matrix, my_symmetric
    LOGICAL, DIMENSION(3)                    :: subdim
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ppack, ppack2, ppack3
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block, rlist

    CALL timeset("rs_get_my_tasks",'I',' ',handle)

    IF ( .NOT. ASSOCIATED ( tasks ) ) &
         CALL stop_program ( "get_my_tasks", "tasks not associated" )

    IF ( PRESENT ( pmat ) .AND. PRESENT ( pcor ) ) THEN
       matrix = .TRUE.
    ELSE
       matrix = .FALSE.
    END IF
    IF(PRESENT(symmetric) )THEN
       my_symmetric = symmetric
    ELSE
       my_symmetric = .TRUE.
    END IF

    ! get number of tasks available locally
    npme = SIZE ( tasks, 2 )
    DO i = 1, SIZE ( tasks, 2)
       IF ( tasks ( 1, i ) <= 0 ) THEN
          npme = i - 1
          EXIT
       END IF
    END DO

    IF ( rs%parallel .AND. .NOT. ALL (rs%perd == 1 ) ) THEN

       DO dir = 1,3
          IF ( rs % perd (dir) /= 1) THEN
             ! bounds of local grid
             lb(dir) = rs%lb_local (dir) + rs%border
             ub(dir) = rs%ub_local (dir) - rs%border
          ELSE
             lb(dir) = rs%lb_local (dir) 
             ub(dir) = rs%ub_local (dir) 
          END IF
       END DO

       ! HACK, what is the MPI doing now?
       subdim = .TRUE.

       CALL mp_cart_sub ( rs % group, subdim, subgroup )
       CALL mp_environ ( subsize, subpos, subgroup )


       ! determine maximum number of tasks
       npmax = npme
       CALL mp_max ( npmax, subgroup )
       nppp = npmax

       ! allocate local arrays
       ALLOCATE(tlist(4,npmax),STAT=stat)
       IF (stat/=0) CALL stop_memory("get_my_tasks","tlist",4*npmax)
       tlist(1:4,1:npme)=tasks(1:4,1:npme)
       tlist(1:4,npme+1:npmax)=0
       IF ( PRESENT ( ival ) ) THEN
          n1 = SIZE ( ival, 1 )
          ALLOCATE(ilist(n1,npmax),STAT=stat)
          IF (stat/=0) CALL stop_memory("get_my_tasks","ilist",n1*npmax)
          ilist(1:n1,1:npme)=ival(1:n1,1:npme)
          ilist(1:n1,npme+1:npmax)=0
          ival=0
       END IF
       IF ( PRESENT ( rval ) ) THEN
          n1 = SIZE ( rval, 1 )
          ALLOCATE(rlist(n1,npmax),STAT=stat)
          IF (stat/=0) CALL stop_memory("get_my_tasks","rlist",n1*npmax)
          rlist(1:n1,1:npme)=rval(1:n1,1:npme)
          rlist(1:n1,npme+1:npmax)=0.0_dp
          rval=0.0_dp
       END IF
       IF ( matrix ) THEN
          n1 = SIZE ( pcor, 1 )
          ALLOCATE(plist(n1,npmax),STAT=stat)
          IF (stat/=0) CALL stop_memory("get_my_tasks","plist",n1*npmax)
          plist(1:n1,1:npme)=pcor(1:n1,1:npme)
          plist(1:n1,npme+1:npmax)=0
          ALLOCATE(pls(n1,npmax),STAT=stat)
          IF (stat/=0) CALL stop_memory("get_my_tasks","pls",n1*npmax)
          pls = 0
          ALLOCATE(ppack(npmax),STAT=stat)
          IF (stat/=0) CALL stop_memory("get_my_tasks","ppack",npmax)
          IF(PRESENT(pmat2))THEN
             ALLOCATE(pls2(n1,npmax),STAT=stat)
             IF (stat/=0) CALL stop_memory("get_my_tasks","pls",n1*npmax)
             pls2= 0
             ALLOCATE(ppack2(npmax),STAT=stat)
             IF (stat/=0) CALL stop_memory("get_my_tasks","ppack",npmax)
          END IF
          IF(PRESENT(pmat3))THEN
             ALLOCATE(pls3(n1,npmax),STAT=stat)
             IF (stat/=0) CALL stop_memory("get_my_tasks","pls",n1*npmax)
             pls3= 0
             ALLOCATE(ppack3(npmax),STAT=stat)
             IF (stat/=0) CALL stop_memory("get_my_tasks","ppack",npmax)
          END IF
       END IF

       ! keep tasks to be processed locally
       ! send remaining tasks to the next processor
       npme = 0
       DO isend = 0, subsize - 1
          IF ( isend == 0 ) THEN
             ileft = SIZE ( tlist, 2)
             DO i = 1, SIZE ( tlist, 2)
                IF ( tlist ( 1, i ) <= 0 ) THEN
                   ileft = i - 1
                   EXIT
                END IF
             END DO
             npmax = ileft
             CALL mp_max ( npmax, subgroup )
          ELSE
             ! count left over tasks
             j = 0
             DO i = 1, npmax
                IF ( tlist ( 1, i ) > 0 ) THEN
                   j = j + 1
                   tlist(:,j) = tlist(:,i)
                   IF ( PRESENT ( ival ) ) ilist(:,j) = ilist(:,i)
                   IF ( PRESENT ( rval ) ) rlist(:,j) = rlist(:,i)
                   IF ( matrix ) plist(:,j) = plist(:,i)
                END IF
             END DO
             tlist(:,j+1:npmax) = 0
             ileft = j
             npmax = ileft
             CALL mp_max ( npmax, subgroup )
             IF ( npmax > nppp ) &
                  CALL stop_program ( "get_my_tasks","npmax too large")
             IF ( npmax == 0 ) EXIT
             ! send/receive tasks
             CALL mp_shift ( tlist ( 1:4, 1:npmax ), subgroup )
             IF ( PRESENT ( ival ) ) THEN
                n1 = SIZE ( ilist, 1 )
                CALL mp_shift ( ilist ( 1:n1, 1:npmax ), subgroup )
             END IF
             IF ( PRESENT ( rval ) ) THEN
                n1 = SIZE ( rlist, 1 )
                CALL mp_shift ( rlist ( 1:n1, 1:npmax ), subgroup )
             END IF
             IF ( matrix ) THEN
                n1 = SIZE ( plist, 1 )
                CALL mp_shift ( plist ( 1:n1, 1:npmax ), subgroup )
             END IF
          END IF
          ! look for tasks do be done on this processor
          ic = 0
          DO i = 1, npmax
             IF ( tlist(1,i) > 0 ) THEN
                IF ( tlist(2,i) >= lb(1) .AND. tlist(2,i) <= ub(1) .AND. &
                     tlist(3,i) >= lb(2) .AND. tlist(3,i) <= ub(2) .AND. &
                     tlist(4,i) >= lb(3) .AND. tlist(4,i) <= ub(3) ) THEN
                   ! found new local task
                   npme = npme + 1
                   IF ( npme > SIZE ( tasks, 2 ) ) &
                        CALL reallocate ( tasks, 1, 4, 1, 2*npme )
                   tasks ( :, npme ) = tlist ( :, i )
                   tlist ( :, i ) = 0
                   IF ( PRESENT ( ival ) ) THEN
                      IF ( npme > SIZE ( ival, 2 ) ) THEN
                         n1 = SIZE ( ival, 1 )
                         CALL reallocate ( ival, 1, n1, 1, 2*npme )
                      END IF
                      ival ( :, npme ) = ilist ( :, i )
                   END IF
                   IF ( PRESENT ( rval ) ) THEN
                      IF ( npme > SIZE ( rval, 2 ) ) THEN
                         n1 = SIZE ( rval, 1 )
                         CALL reallocate ( rval, 1, n1, 1, 2*npme )
                      END IF
                      rval ( :, npme ) = rlist ( :, i )
                   END IF
                   IF ( matrix ) THEN
                      ic = ic + 1
                      IF ( ic > SIZE ( pcor, 2 ) ) THEN
                         n1 = SIZE ( pcor, 1 )
                         CALL reallocate ( pcor, 1, n1, 1, 2*ic )
                      END IF
                      pcor ( :, ic ) = plist ( :, i )
                   END IF
                END IF
             END IF
          END DO
          ! now distribute matrix blocks
          IF ( matrix .AND. isend /= 0 ) THEN
             ! eliminate double entries
             IF ( ic > 1 ) THEN
                j=1
                DO i = 2, ic
                   DO k = j, 1, -1
                      IF ( pcor(1,i)==pcor(1,k) .AND. pcor(2,i)==pcor(2,k) ) EXIT
                      IF ( k == 1 ) THEN
                         j = j + 1
                         pcor(:,j)=pcor(:,i)
                      END IF
                   END DO
                END DO
                ic = j
             END IF
             ! eliminate blocks that are already local
             j = 0
             DO i = 1, ic
                IF(my_symmetric) THEN
                   IF ( pcor(1,i) <= pcor(2,i) ) THEN
                      arow = pcor(1,i)
                      acol = pcor(2,i)
                   ELSE
                      arow = pcor(2,i)
                      acol = pcor(1,i)
                   END IF
                ELSE
                   arow = pcor(1,i)
                   acol = pcor(2,i)
                END IF
                CALL get_block_node(matrix=pmat,&
                     block_row=arow,&
                     block_col=acol,&
                     BLOCK=p_block)
                IF ( .NOT. ASSOCIATED ( p_block ) ) THEN
                   j = j + 1
                   pcor(:,j) = pcor(:,i)
                END IF
             END DO
             ic = j
             icmax = ic
             CALL mp_max ( icmax, subgroup )
             ! are there missing blocks?
             IF ( icmax > 0 ) THEN
                ! send block coordinates back to original pe
                pls(:,1:ic) = pcor(:,1:ic)
                pls(:,ic+1:) = 0
                CALL mp_shift ( pls ( 1:2, 1:icmax ), subgroup, -isend )
                ! pack the data to be sent
                CALL pack_matrix ( pmat, pls, ppack, plength, subgroup ,symmetric=my_symmetric)
                CALL mp_shift ( ppack ( 1:plength ), subgroup, isend )
                CALL mp_shift ( pls ( 1:2, 1:icmax ), subgroup, isend )
                CALL unpack_matrix ( pmat, pcor, pls, ic, ppack , symmetric=my_symmetric,error=error)
                IF(PRESENT(pmat2)) THEN
                   pls2(:,1:ic) = pcor(:,1:ic)
                   pls2(:,ic+1:) = 0
                   CALL mp_shift ( pls2 ( 1:2, 1:icmax ), subgroup, -isend )
                   ! pack the data to be sent
                   CALL pack_matrix ( pmat2, pls2, ppack2, plength2, subgroup,symmetric=my_symmetric )
                   CALL mp_shift ( ppack2 ( 1:plength2 ), subgroup, isend )
                   CALL mp_shift ( pls2 ( 1:2, 1:icmax ), subgroup, isend )
                   CALL unpack_matrix ( pmat2, pcor, pls2, ic, ppack2,symmetric=my_symmetric ,error=error)
                END IF
                IF(PRESENT(pmat3)) THEN
                   pls3(:,1:ic) = pcor(:,1:ic)
                   pls3(:,ic+1:) = 0
                   CALL mp_shift ( pls3 ( 1:2, 1:icmax ), subgroup, -isend )
                   ! pack the data to be sent
                   CALL pack_matrix ( pmat3, pls3, ppack3, plength3, subgroup,symmetric=my_symmetric )
                   CALL mp_shift ( ppack3 ( 1:plength3 ), subgroup, isend )
                   CALL mp_shift ( pls3 ( 1:2, 1:icmax ), subgroup, isend )
                   CALL unpack_matrix ( pmat3, pcor, pls3, ic, ppack3,symmetric=my_symmetric ,error=error)
                END IF
             END IF
          END IF
       END DO

       ! Release the communicator
       CALL mp_comm_free ( subgroup )

       ! are all tasks distributed?
       IF ( ANY ( tlist(1,:) > 0 ) ) THEN
          WRITE(6,*) "# of unallocated tasks is ",SIZE ( tlist , 2) 
          CALL stop_program ( "get_my_tasks", "left over tasks" )
       END IF

       DEALLOCATE(tlist,STAT=stat)
       IF (stat/=0) CALL stop_memory("get_my_tasks","tlist")
       IF ( PRESENT ( ival ) ) THEN
          DEALLOCATE(ilist,STAT=stat)
          IF (stat/=0) CALL stop_memory("get_my_tasks","ilist")
       END IF
       IF ( PRESENT ( rval ) ) THEN
          DEALLOCATE(rlist,STAT=stat)
          IF (stat/=0) CALL stop_memory("get_my_tasks","rlist")
       END IF
       IF ( matrix ) THEN
          DEALLOCATE(plist,STAT=stat)
          IF (stat/=0) CALL stop_memory("get_my_tasks","plist")
          DEALLOCATE(pls,STAT=stat)
          IF (stat/=0) CALL stop_memory("get_my_tasks","pls")
          DEALLOCATE(ppack,STAT=stat)
          IF (stat/=0) CALL stop_memory("get_my_tasks","ppack")
          IF(PRESENT(pmat2)) THEN
             DEALLOCATE(pls2,STAT=stat)
             IF (stat/=0) CALL stop_memory("get_my_tasks","pls")
             DEALLOCATE(ppack2,STAT=stat)
             IF (stat/=0) CALL stop_memory("get_my_tasks","ppack")
          END IF
          IF(PRESENT(pmat3)) THEN
             DEALLOCATE(pls3,STAT=stat)
             IF (stat/=0) CALL stop_memory("get_my_tasks","pls")
             DEALLOCATE(ppack3,STAT=stat)
             IF (stat/=0) CALL stop_memory("get_my_tasks","ppack")
          END IF
       END IF

    ELSE

       ! fully replicated grids, each processor can process all its tasks

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE rs_get_my_tasks

!***************************************************************************
SUBROUTINE pack_matrix ( pmat, pls, ppack , plength, group , symmetric)

    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: pmat
    INTEGER, DIMENSION(:, :), POINTER        :: pls
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ppack
    INTEGER, INTENT(OUT)                     :: plength
    INTEGER, INTENT(IN)                      :: group
    LOGICAL, INTENT(IN), OPTIONAL            :: symmetric

    INTEGER                                  :: acol, arow, i, j, k, n, nc, nr
    LOGICAL                                  :: my_symmetric
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block

  my_symmetric = .TRUE.
  IF(PRESENT(symmetric)) my_symmetric = symmetric
  plength = 0
  DO i = 1, SIZE ( pls, 2 )
    IF ( pls ( 1, i ) == 0 ) EXIT
    IF(my_symmetric) THEN
      IF ( pls(1,i) <= pls(2,i) ) THEN
         arow = pls(1,i)
         acol = pls(2,i)
      ELSE
         arow = pls(2,i)
         acol = pls(1,i)
      END IF
    ELSE
       arow = pls(1,i)
       acol = pls(2,i)
    END IF
    CALL get_block_node(matrix=pmat,&
                        block_row=arow,&
                        block_col=acol,&
                        BLOCK=p_block)
    IF ( .NOT. ASSOCIATED ( p_block ) ) &
       CALL stop_program ( "pack_matrix", "Matrix block not found" )
    nr = SIZE ( p_block, 1 )
    nc = SIZE ( p_block, 2 )
    pls ( 1, i ) = nr
    pls ( 2, i ) = nc
    n = nc * nr
    IF ( plength + n > SIZE ( ppack ) ) THEN
      CALL reallocate ( ppack, 1, plength+5*n )
    END IF
    DO j = 1, nc
      DO k = 1, nr
        plength = plength + 1
        ppack(plength) = p_block(k,j)
      END DO
    END DO
  END DO
  CALL mp_max ( plength, group )
  IF ( plength > SIZE ( ppack ) ) THEN
    CALL reallocate ( ppack, 1, plength )
  END IF

END SUBROUTINE pack_matrix

SUBROUTINE unpack_matrix ( pmat, pcor, pls, ic, ppack , symmetric, error)

    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: pmat
    INTEGER, DIMENSION(:, :), POINTER        :: pcor, pls
    INTEGER, INTENT(IN)                      :: ic
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ppack
    LOGICAL, INTENT(IN), OPTIONAL            :: symmetric
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: acol, arow, i, j, k, nc, nr, &
                                                pl
    LOGICAL                                  :: my_symmetric
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block

  my_symmetric = .TRUE.
  IF(PRESENT(symmetric)) my_symmetric = symmetric
  pl = 0
  DO i = 1, ic
    IF(my_symmetric) THEN
      IF ( pcor(1,i) <= pcor(2,i) ) THEN
         arow = pcor(1,i)
         acol = pcor(2,i)
      ELSE
         arow = pcor(2,i)
         acol = pcor(1,i)
      END IF
    ELSE
       arow = pcor(1,i)
       acol = pcor(2,i)
    END IF
    nr = pls(1,i)
    nc = pls(2,i)
    NULLIFY ( p_block )
    CALL add_block_node ( pmat, arow, acol, p_block ,error=error)
    DO j = 1, nc
      DO k = 1, nr
        pl = pl + 1
        p_block(k,j) = ppack(pl)
      END DO
    END DO

  END DO

END SUBROUTINE unpack_matrix

!***************************************************************************

SUBROUTINE rs_get_loop_vars ( npme, ival, natom_pairs, asets, atasks )

    INTEGER, INTENT(IN)                      :: npme
    INTEGER, DIMENSION(6, npme), INTENT(IN)  :: ival
    INTEGER, INTENT(OUT)                     :: natom_pairs
    INTEGER, DIMENSION(:, :), POINTER        :: asets, atasks

    INTEGER                                  :: iatom, iatom_old, iset, &
                                                iset_old, itask, jatom, &
                                                jatom_old, jset, jset_old, &
                                                nset_pairs

     IF(SIZE(asets,2) < npme) CALL reallocate(asets,1,2,1,npme)
     IF(SIZE(atasks,2) < npme) CALL reallocate(atasks,1,2,1,npme)
     natom_pairs = 0
     nset_pairs = 0
     iatom_old = 0
     jatom_old = 0
     DO itask = 1, npme
        iatom = ival(1,itask)
        jatom = ival(2,itask)
        iset = ival(3,itask)
        jset = ival(4,itask)
        IF ( iatom /= iatom_old .OR. jatom /= jatom_old ) THEN
          IF(natom_pairs>0) asets(2,natom_pairs) = nset_pairs
          natom_pairs = natom_pairs + 1
          iatom_old = iatom
          jatom_old = jatom
          IF(nset_pairs>0) atasks(2,nset_pairs) = itask - 1
          nset_pairs = nset_pairs + 1
          asets(1,natom_pairs) = nset_pairs
          atasks(1,nset_pairs) = itask
          iset_old = iset
          jset_old = jset
        ELSE IF ( iset /= iset_old .OR. jset /= jset_old ) THEN
          atasks(2,nset_pairs) = itask - 1
          nset_pairs = nset_pairs + 1
          atasks(1,nset_pairs) = itask
          iset_old = iset
          jset_old = jset
        END IF
     END DO
     IF(natom_pairs>0) asets(2,natom_pairs) = nset_pairs
     IF(nset_pairs>0) atasks(2,nset_pairs) = npme

END SUBROUTINE rs_get_loop_vars

!******************************************************************************

END MODULE realspace_task_selection

!******************************************************************************
