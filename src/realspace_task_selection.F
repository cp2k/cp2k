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
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_alltoall
  USE realspace_grid_types,            ONLY: realspace_grid_type
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             get_block_node,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rs_get_my_tasks,&
            pair2int,&
            int2pair,&
            distribute_matrix

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'realspace_task_selection'

!-----------------------------------------------------------------------------!

CONTAINS

!!****f* realspace_task_selection/pair2int *
!!
!!   NAME
!!     pair2int
!!
!!   FUNCTION
!!     converts a pgf index pair (ipgf,iset,iatom),(jpgf,jset,jatom)
!!     to a unique integer.
!!     a list of integers can be sorted, and will result in a list of pgf pairs
!!     for which all atom pairs are grouped, and for each atom pair all set pairs are grouped
!!     and for each set pair, all pgfs are grouped
!!
!!   NOTES
!!     will hopefully not overflow any more
!!
!!   INPUTS
!!    -
!!    -
!!
!!   MODIFICATION HISTORY
!!     11.2007 created [Joost]
!!
!!*** **********************************************************************
  SUBROUTINE pair2int(res,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)
    INTEGER(kind=int_8)                      :: res
    INTEGER                                  :: iatom, jatom, iset, jset, &
                                                ipgf, jpgf, natom, maxset, &
                                                maxpgf

    INTEGER(kind=int_8)                      :: maxpgf8, maxset8, natom8

    natom8=natom ; maxset8=maxset ; maxpgf8=maxpgf

    res=((iatom-1)*natom8+jatom-1)*maxset8**2*maxpgf8**2 + &
        ((iset-1)*maxset8+jset-1)*maxpgf8**2+ &
         (ipgf-1)*maxpgf8+jpgf-1

  END SUBROUTINE pair2int

  SUBROUTINE int2pair(res,iatom,jatom,iset,jset,ipgf,jpgf,natom,maxset,maxpgf)
    INTEGER(kind=int_8)                      :: res
    INTEGER                                  :: iatom, jatom, iset, jset, &
                                                ipgf, jpgf, natom, maxset, &
                                                maxpgf

    INTEGER(kind=int_8)                      :: iatom8, ijatom, ijpgf, ijset, &
                                                ipgf8, iset8, jatom8, jpgf8, &
                                                jset8, maxpgf8, maxset8, &
                                                natom8, tmp

    natom8=natom ; maxset8=maxset ; maxpgf8=maxpgf

    ijatom=res/(maxpgf8**2*maxset8**2)
    iatom8=ijatom/natom8+1
    jatom8=MOD(ijatom,natom8)+1
    tmp=MOD(res,maxpgf8**2*maxset8**2)
    ijset=tmp/maxpgf8**2
    iset8=ijset/maxset8+1
    jset8=MOD(ijset,maxset8)+1
    ijpgf=MOD(tmp,maxpgf8**2)
    ipgf8=ijpgf/maxpgf8+1
    jpgf8=MOD(ijpgf,maxpgf8)+1

    iatom=iatom8 ; jatom=jatom8; iset=iset8 ; jset=jset8 
    ipgf=ipgf8 ; jpgf=jpgf8

  END SUBROUTINE int2pair

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
!!     - tasks: the task set generated on this processor
!!     - npme: Number of tasks for local processing
!!
!!   AUTHOR
!!     MattW 21/11/2007
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************

  SUBROUTINE rs_get_my_tasks ( rs, ntasks, natoms,&
       maxset,maxpgf,  tasks, rval, atom_pair_send, atom_pair_recv,&
       symmetric, error)

    TYPE(realspace_grid_type), POINTER       :: rs
    INTEGER                                  :: ntasks, natoms, maxset, maxpgf
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rval
    INTEGER, DIMENSION(:), POINTER           :: atom_pair_send, atom_pair_recv
    LOGICAL, INTENT(IN), OPTIONAL            :: symmetric
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: handle, i, j, k, l, &
                                                ntasks_recv, stat, task_dim
    INTEGER(kind=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: recv_buf_i, send_buf_i
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks_recv
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index, recv_disps, &
                                                recv_sizes, send_disps, &
                                                send_sizes, taskid
    LOGICAL                                  :: my_symmetric, scatter, send
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_buf_r, send_buf_r

    CALL timeset("rs_get_my_tasks",'I',' ',handle)

    IF ( .NOT. ASSOCIATED ( tasks ) ) &
         CALL stop_program ( "get_my_tasks", "tasks not associated" )

    IF(PRESENT(symmetric) )THEN
       my_symmetric = symmetric
    ELSE
       my_symmetric = .TRUE.
    END IF

    ! allocate local arrays
    ALLOCATE(send_sizes(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
    ALLOCATE(recv_sizes(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
    ALLOCATE(send_disps(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
    ALLOCATE(recv_disps(rs%group_size),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)

    ! communication step 1 : compute the number of tasks we will send to each cpu
    !                       and send this information around

    ! tell targets how many tasks they're going to receive
    send_sizes = 1
    recv_sizes = 1
    send_disps(1)=0
    recv_disps(1)=0
    DO i = 2, rs % group_size
       send_disps(i) = send_disps(i-1) + send_sizes(i-1)
       recv_disps(i) = recv_disps(i-1) + recv_sizes(i-1)
    ENDDO

    ! allocate buffers
    ALLOCATE(send_buf_i(SUM(send_sizes)),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf",rs%group_size)
    ALLOCATE(recv_buf_i(SUM(recv_sizes)),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf",rs%group_size)

    ! fill send buffer
    send_buf_i=0
    DO i = 1, ntasks
       send_buf_i( tasks(1,i)+1 ) = send_buf_i( tasks(1,i)+1 ) + 1 
    END DO

    CALL mp_alltoall(send_buf_i, send_sizes, send_disps, recv_buf_i, recv_sizes, recv_disps, rs % group)

    ! communication step 2 : pack the tasks, and send them around

    task_dim = SIZE(tasks,1)

    send_sizes = 0
    send_disps = 0
    recv_sizes = 0
    recv_disps = 0

    send_sizes(1) = send_buf_i(1) * task_dim
    recv_sizes(1) = recv_buf_i(1) * task_dim
    DO i = 2,rs%group_size
       send_sizes(i) = send_buf_i(i) * task_dim
       recv_sizes(i) = recv_buf_i(i) * task_dim
       send_disps(i)=send_disps(i-1)+send_sizes(i-1)           
       recv_disps(i)=recv_disps(i-1)+recv_sizes(i-1)
    ENDDO

    ! deallocate old send/recv buffers
    DEALLOCATE(send_buf_i,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf_i",rs%group_size)
    DEALLOCATE(recv_buf_i,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf_i",rs%group_size)

    ! allocate them with new sizes
    ALLOCATE(send_buf_i(SUM(send_sizes)),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf",rs%group_size)
    ALLOCATE(recv_buf_i(SUM(recv_sizes)),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf",rs%group_size)

    ! do packing
    send_buf_i = 0
    send_sizes = 0
    DO j = 1,ntasks
       i = tasks(1,j)+1
       DO k=1,task_dim
          send_buf_i(send_disps(i)+send_sizes(i)+k)=tasks(k,j)
       ENDDO
       send_sizes(i)=send_sizes(i) + task_dim
    ENDDO

    ! do communication
    CALL mp_alltoall(send_buf_i, send_sizes, send_disps,&
         recv_buf_i, recv_sizes, recv_disps, rs % group)

    DEALLOCATE(send_buf_i,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf",rs%group_size)

    ntasks_recv=SUM(recv_sizes)/task_dim
    ALLOCATE(tasks_recv(task_dim,ntasks_recv),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","tasks_recv",rs%group_size)

    ! do unpacking
    l = 0
    DO i = 1,rs % group_size
       DO j = 0,recv_sizes(i)/task_dim-1
          l = l + 1
          DO k=1,task_dim 
             tasks_recv(k,l)=recv_buf_i(recv_disps(i)+j*task_dim+k)
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(recv_buf_i,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf",rs%group_size)

    !***************************************************************************
    ! send rvals (to be removed :-)

    ! take care of the new task_dim in the send_sizes
    send_sizes=(send_sizes/task_dim)*SIZE(rval,1)
    recv_sizes=(recv_sizes/task_dim)*SIZE(rval,1)
    send_disps=(send_disps/task_dim)*SIZE(rval,1)
    recv_disps=(recv_disps/task_dim)*SIZE(rval,1)
    task_dim=SIZE(rval,1)

    ALLOCATE(send_buf_r(SUM(send_sizes)),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf_r",rs%group_size)
    ALLOCATE(recv_buf_r(SUM(recv_sizes)),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf_r",rs%group_size)

    !do packing
    send_sizes = 0
    DO j = 1,ntasks
       i = tasks(1,j)+1
       DO k=1,task_dim
          send_buf_r(send_disps(i)+send_sizes(i)+k)=rval(k,j)
       ENDDO
       send_sizes(i)=send_sizes(i) + task_dim
    ENDDO

    ! do communication
    CALL mp_alltoall(send_buf_r, send_sizes, send_disps,&
         recv_buf_r, recv_sizes, recv_disps, rs % group)

    DEALLOCATE(send_buf_r,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf_r",rs%group_size)
    DEALLOCATE(rval,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","rval",rs%group_size)
    ALLOCATE(rval(task_dim,SUM(recv_sizes)/task_dim),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","rval",rs%group_size)

    ! do unpacking
    l = 0
    DO i = 1,rs % group_size
       DO j = 0,recv_sizes(i)/task_dim-1
          l = l + 1
          DO k=1,task_dim
             rval(k,l)=recv_buf_r(recv_disps(i)+j*task_dim+k)
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(recv_buf_r,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf_r",rs%group_size)

    !******************************************************************************
    ! prepare to send matrices

    ! calculate list of atom pairs I'll be sending
    ALLOCATE(atom_pair_send(ntasks),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","atom_pair_send",rs%group_size)

    ALLOCATE(atom_pair_recv(ntasks_recv),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","atom_pair_recv",rs%group_size)

    send = .TRUE.

    CALL get_atom_pair ( atom_pair_send, tasks, send,&
         my_symmetric, natoms, maxset, maxpgf )

    send = .FALSE.

    CALL get_atom_pair ( atom_pair_recv, tasks_recv, send,&
         my_symmetric, natoms, maxset, maxpgf )

    !********************************************************************************************

    DEALLOCATE(send_sizes,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
    DEALLOCATE(recv_sizes,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
    DEALLOCATE(send_disps,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
    DEALLOCATE(recv_disps,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)

    DEALLOCATE(tasks,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","tasks",rs%group_size)

    ALLOCATE(taskid(ntasks_recv),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","taskid",rs%group_size)
    ALLOCATE(INDEX(ntasks_recv),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","index",rs%group_size)

    taskid=tasks_recv(3,:)
    CALL sort(taskid,SIZE(taskid),index)

    DO k=1,SIZE(tasks_recv,1)
       tasks_recv(k,:)=tasks_recv(k,index)
    ENDDO

    DO k=1,SIZE(rval,1)
       rval(k,:)=rval(k,index)
    ENDDO

    DEALLOCATE(taskid,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","task_id",rs%group_size)
    DEALLOCATE(index,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","index",rs%group_size)

    tasks=>tasks_recv
    ntasks=ntasks_recv

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE rs_get_my_tasks

! ****************************************************************************************

SUBROUTINE get_atom_pair ( atom_pair, my_tasks, send, my_symmetric, natoms, maxset, maxpgf)

    INTEGER, DIMENSION(:), POINTER           :: atom_pair
    INTEGER(kind=int_8), DIMENSION(:, :), &
      POINTER                                :: my_tasks
    LOGICAL                                  :: send, my_symmetric
    INTEGER                                  :: natoms, maxset, maxpgf

    INTEGER                                  :: acol, arow, i, iatom, ipgf, &
                                                iset, j, jatom, jpgf, jset, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: index

! calculate list of atom pairs
! fill pair list taking into acount symmetry

  atom_pair = 0
  DO i = 1,SIZE(atom_pair)
     CALL int2pair(my_tasks(3,i),iatom,jatom,iset,jset,ipgf,jpgf,natoms,maxset,maxpgf)
     IF( my_symmetric ) THEN
        IF(iatom.LE.jatom) THEN
           arow = iatom
           acol = jatom
        ELSE
           arow = jatom
           acol = iatom
        ENDIF
     ELSE
        arow = iatom
        acol = jatom
     ENDIF
     IF ( send ) THEN
        atom_pair(i) = my_tasks(1,i)*natoms*natoms+(arow-1)*natoms + (acol-1)
     ELSE
        atom_pair(i) = my_tasks(2,i)*natoms*natoms+(arow-1)*natoms + (acol-1)
     ENDIF
  ENDDO

  ALLOCATE(INDEX(SIZE(atom_pair)),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_atom_pair","index",SIZE(atom_pair))  

  ! find unique atom pairs that I'm sending
  CALL sort(atom_pair,SIZE(atom_pair),index)

  DEALLOCATE(index,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_atom_pair","index",SIZE(atom_pair))  

  IF (SIZE(atom_pair)>1) THEN
     j=1
     ! first atom pair must be allowed
     DO i = 2,SIZE(atom_pair)
        IF( atom_pair(i) .GT. atom_pair(i-1) ) THEN 
           j = j + 1
           atom_pair(j) = atom_pair(i)
        ENDIF
     ENDDO
     ! reallocate the atom pair list
     CALL  reallocate(atom_pair,1,j)
  ENDIF

END SUBROUTINE get_atom_pair

! **************************************************************************

SUBROUTINE distribute_matrix( rs, pmat, atom_pair_send, atom_pair_recv, natoms, scatter, error, hmat )

    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(real_matrix_type), POINTER          :: pmat
    INTEGER, DIMENSION(:), POINTER           :: atom_pair_send, atom_pair_recv
    INTEGER                                  :: natoms
    LOGICAL                                  :: scatter
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: hmat

    INTEGER                                  :: acol, arow, handle, i, j, k, &
                                                l, ncol, nrow, pair, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: recv_disps, recv_sizes, &
                                                send_disps, send_sizes
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_buf_r, send_buf_r
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block

  CALL timeset("distribute_matrix",'I',' ',handle)

  ! allocate local arrays
  ALLOCATE(send_sizes(rs%group_size),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
  ALLOCATE(recv_sizes(rs%group_size),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
  ALLOCATE(send_disps(rs%group_size),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
  ALLOCATE(recv_disps(rs%group_size),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)

  ! set up send buffer sizes

  send_sizes = 0
  DO i = 1, SIZE(atom_pair_send)

     ! proc we're sending this block to
     k = atom_pair_send(i) / natoms**2 + 1

     pair = MOD(atom_pair_send(i),natoms**2)

     arow = pair / natoms + 1
     acol = MOD(pair, natoms) + 1

     nrow = pmat%last_row(arow) - pmat%first_row(arow) + 1
     ncol = pmat%last_col(acol) - pmat%first_col(acol) + 1       

     send_sizes(k) = send_sizes(k) + nrow * ncol        

  ENDDO

  send_disps = 0
  DO i = 2, rs % group_size
     send_disps(i) = send_disps(i-1) + send_sizes(i-1)
  ENDDO

  ALLOCATE(send_buf_r(SUM(send_sizes)),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf_r",SUM(send_sizes))

  send_buf_r = 0

  ! do packing
  send_sizes=0
  DO i = 1, SIZE(atom_pair_send)

     l = atom_pair_send(i) / natoms**2 + 1

     pair = MOD(atom_pair_send(i),natoms**2)

     arow = pair / natoms + 1
     acol = MOD(pair, natoms) + 1

     nrow = pmat%last_row(arow)-pmat%first_row(arow) + 1
     ncol = pmat%last_col(acol)-pmat%first_col(acol) + 1

     CALL get_block_node(matrix=pmat,&
          block_row=arow,&
          block_col=acol,&
          BLOCK=p_block)
     IF ( .NOT. ASSOCIATED ( p_block ) ) THEN 
        CALL stop_program ( "pack_matrix almost there", "Matrix block not found" )
     ENDIF

     DO k = 1, ncol
        DO j = 1, nrow
           send_buf_r(send_disps(l)+send_sizes(l)+j+(k-1)*nrow) = p_block(j,k)
        ENDDO
     ENDDO
     send_sizes(l)=send_sizes(l)+nrow*ncol
  ENDDO

  ! set up recv buffer

  recv_sizes = 0
  DO i = 1, SIZE(atom_pair_recv)

     ! proc we're receiving this data from
     k = atom_pair_recv(i) / natoms**2 + 1

     pair = MOD(atom_pair_recv(i),natoms**2)

     arow = pair / natoms + 1
     acol = MOD(pair,natoms) + 1

     nrow = pmat%last_row(arow) - pmat%first_row(arow) + 1
     ncol = pmat%last_col(acol) - pmat%first_col(acol) + 1       

     recv_sizes(k) = recv_sizes(k) + nrow * ncol        

  ENDDO

  recv_disps = 0     
  DO i = 2, rs % group_size
     recv_disps(i) = recv_disps(i-1) + recv_sizes(i-1)
  ENDDO

  ALLOCATE(recv_buf_r(SUM(recv_sizes)),STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf_r",SUM(recv_sizes))

  recv_buf_r = 0

  ! do communication
  CALL mp_alltoall(send_buf_r, send_sizes, send_disps,&
       recv_buf_r, recv_sizes, recv_disps, rs % group)

  !do unpacking
  recv_sizes=0
  DO i = 1, SIZE(atom_pair_recv)
     l = atom_pair_recv(i) / natoms**2 + 1
     pair = MOD(atom_pair_recv(i),natoms**2)

     arow = pair / natoms + 1
     acol = MOD(pair, natoms) + 1

     nrow = pmat%last_row(arow)-pmat%first_row(arow) + 1
     ncol = pmat%last_col(acol)-pmat%first_col(acol) + 1

     CALL get_block_node(matrix=pmat,&
          block_row=arow,&
          block_col=acol,&
          BLOCK=p_block)

     IF ( PRESENT ( hmat ) ) THEN
     CALL get_block_node(matrix=hmat,&
          block_row=arow,&
          block_col=acol,&
          BLOCK=h_block)
     ENDIF

     IF ( ASSOCIATED ( p_block ) .AND. l .NE. rs%my_pos+1) THEN
        IF( scatter ) THEN
           WRITE(6,*) arow,acol,l 
           CALL stop_program ( "unpack_matrix here and new", "Matrix block already present" )
        ENDIF
     ELSE
        CALL add_block_node ( pmat, arow, acol, p_block ,error=error)
     ENDIF

     DO k = 1, ncol
        DO j = 1, nrow
           IF ( scatter ) THEN
              p_block(j,k) = recv_buf_r( recv_disps(l) + recv_sizes(l)+ j + (k-1)*nrow )
           ELSE
              h_block(j,k) = h_block(j,k) + recv_buf_r( recv_disps(l) + recv_sizes(l)+ j + (k-1)*nrow )
           ENDIF
        ENDDO
     ENDDO
     recv_sizes(l)=recv_sizes(l)+nrow*ncol
  ENDDO

  DEALLOCATE(send_buf_r,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_buf_r",rs%group_size)
  DEALLOCATE(recv_buf_r,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_buf_r",rs%group_size)

  DEALLOCATE(send_sizes,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_sizes",rs%group_size)
  DEALLOCATE(recv_sizes,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_sizes",rs%group_size)
  DEALLOCATE(send_disps,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","send_disps",rs%group_size)
  DEALLOCATE(recv_disps,STAT=stat)
  IF (stat/=0) CALL stop_memory("get_my_tasks","recv_disps",rs%group_size)

  CALL timestop(handle)
  

END SUBROUTINE distribute_matrix


!******************************************************************************

END MODULE realspace_task_selection

!******************************************************************************
