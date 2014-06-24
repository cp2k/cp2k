!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulationsi        !
!   Copyright (C) 2000 - 2014  CP2K developers group group!
!-----------------------------------------------------------------------------!
!> \brief   Process mapping support for cp2k
!> \author  Christiane Pousa Ribeiro
!> \date    2012-01-17
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2012-01-17
MODULE ma_process_mapping

  USE ma_errors,                       ONLY: ma_error_allocation,&
                                             ma_error_stop
  USE ma_kinds,                        ONLY: dp,&
                                             int_size
  USE ma_topology,                     ONLY: allocated_topology,&
                                             ma_2dgrid_dimensions,&
                                             ma_allocated_topology,&
                                             ma_comm_graph,&
                                             ma_compute_contention,&
                                             ma_get_job_allocation,&
                                             ma_net_topology,&
                                             net_topology
  USE machine_architecture,            ONLY: ma_get_ncores,&
                                             ma_get_nmachines
  USE machine_architecture_types,      ONLY: &
       cannon, cannon_graph, complete_graph, hilbert, hilbert_peano, own, &
       packed, peano, round_robin, snake, switch
  USE message_passing,                 ONLY: mp_environ,&
                                             mp_reordering,&
                                             mp_sum

IMPLICIT NONE

 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ma_process_mapping'

  PUBLIC :: ma_mpi_reordering

CONTAINS


! *****************************************************************************
!> \brief Apply a process mapping to CP2K MPI
!> \param mp_comm [input/output] : handle of the default communicator
!> \param reorder [input] : the strategy to be applied
!> \par History
!>      2.2012 created [ Christiane Pousa ]
! *****************************************************************************
SUBROUTINE ma_mpi_reordering(reorder,mp_comm)
   INTEGER, INTENT(INOUT)    :: mp_comm
   INTEGER, INTENT(IN)       :: reorder
  
#if defined(__parallel)
   INTEGER                   :: new_comm, numtask,taskid,&
                                ncol,nrow

  IF( reorder > hilbert_peano) THEN
          CALL ma_net_topology(mp_comm)
          CALL ma_allocated_topology(mp_comm)
  ENDIF

  SELECT CASE (reorder)
  CASE (packed)
   CALL ma_packed(new_comm, mp_comm)
  CASE (hilbert)
   CALL ma_hilbert_curve(new_comm, mp_comm) 
  CASE (peano)
   CALL ma_peano_curve(new_comm, mp_comm)
  CASE (snake)
   CALL ma_snake_curve(new_comm, mp_comm)
  CASE (round_robin)
   CALL ma_round_robin(new_comm, mp_comm) 
  CASE (hilbert_peano)
    CALL mp_environ(numtask,taskid,mp_comm)
    CALL ma_2dgrid_dimensions(ncol,nrow,numtask)
   IF ((ncol .LT. 4) .OR. (nrow .LT. 4)) THEN   
    CALL ma_packed(new_comm,mp_comm)  
   ELSE
    CALL ma_hilbert_peano(new_comm,mp_comm)
   ENDIF
  CASE (switch)
   CALL ma_switch(new_comm,mp_comm)
  CASE (cannon)
    CALL mp_environ(numtask,taskid,mp_comm)
    CALL ma_2dgrid_dimensions(ncol,nrow,numtask)
   IF ((ncol .LT. 4) .OR. (nrow .LT. 4)) THEN
    CALL ma_switch(new_comm,mp_comm,kernel=cannon_graph)
   ELSE
    CALL ma_cannon(new_comm,mp_comm) 
   ENDIF
  CASE (own)
    CALL ma_designed(new_comm,mp_comm)
  END SELECT 
  mp_comm = new_comm
#endif
END SUBROUTINE

! *****************************************************************************
!> \brief Find the closest square grid for the current number of processes
!> \brief the grid will be 2**N x 2**N
!> \param numtask [input] : total number of MPI tasks
!> \param n [output] : the grid size in one dimension
!> \par History
!>      2.2012 created [ Christiane Pousa ]
! *****************************************************************************
FUNCTION ma_find_min_square_grid(numtask) RESULT (n)
   INTEGER, INTENT(IN)    :: numtask
   INTEGER                :: n


#if !defined(__parallel)
   n = numtask
#else
    INTEGER                                  :: tmp_n, square
    REAL                                     :: ntasks
    LOGICAL                                  :: power_two

    ntasks = numtask
    n = CEILING(SQRT(REAL((ntasks),KIND=dp)))
    power_two = .FALSE.

    tmp_n = 1
    DO WHILE (.NOT. power_two)
      square = (2**tmp_n)*(2**tmp_n)
      IF ( n*n .GT. square) THEN
          tmp_n = tmp_n + 1
      ELSE
        n = square
        power_two = .TRUE. 
      ENDIF
    ENDDO
#endif

END FUNCTION ma_find_min_square_grid

! *****************************************************************************
!> \brief Find the max square grid that holds the current number of processes
!> \brief the grid will be 2**N x 2**N
!> \param numtask [input] : total number of MPI tasks
!> \param n [output] : the grid size in one dimension
!> \par History
!>      6.2012 created [ Christiane Pousa ]
! *****************************************************************************
FUNCTION ma_find_max_square_grid(numtask) RESULT(n)
   INTEGER, INTENT(IN)                      :: numtask
   INTEGER                                  :: n

#if defined(__parallel)
   INTEGER                                  :: tmp_n, square, old

   IF (numtask .NE. 0) THEN

     n = CEILING(SQRT(REAL((numtask),KIND=dp)))
     old = n
     tmp_n = 1
     square = (2**tmp_n)*(2**tmp_n)
     DO WHILE (square .LE. n*n)
       tmp_n = tmp_n + 1
       old = square 
       square = (2**tmp_n)*(2**tmp_n)
     END DO
     n = old
   ELSE
     n = 0
   END IF
#else 
   n = numtask
#endif

END FUNCTION ma_find_max_square_grid


! *****************************************************************************
!> \brief Find the closest 3^level grid for the current number of processes
!> \brief the grid will be 3**N x 3**N
!> \param numtask [input] : the total number of MPI tasks
!> \param n [output] : the grid level in one dimension
!> \par History
!>      4.2012 created [ Christiane Pousa ]
! *****************************************************************************
FUNCTION ma_find_grid_level(numtask) RESULT (n)
   INTEGER, INTENT(IN)   :: numtask
   INTEGER                :: n

#if !defined(__parallel)
   n = numtask
#else 
    INTEGER                                  :: level,grid

    level = 1
    grid = (3**level)*(3**level)
    DO WHILE(numtask .GT. grid)
     level = level + 1
     grid = (3**level)*(3**level) 
    ENDDO

    n = level

#endif
END FUNCTION

SUBROUTINE move(angle,coordX,coordY,pos)
    INTEGER                                  :: angle, coordX, coordY, pos

   DO WHILE (angle .GT. 270) 
    angle = angle - 360  
   ENDDO 
   DO WHILE (angle .LT. 0) 
    angle = angle + 360 
   ENDDO
   IF  (angle .EQ. 0) THEN
       coordX = coordX + 1  
   ELSE IF (angle .EQ. 90) THEN
       coordY = coordY + 1  
   ELSE IF (angle .EQ. 180) THEN 
       coordX = coordX - 1  
   ELSE IF (angle .EQ. 270) THEN
        coordY = coordY - 1  
   ENDIF
   pos = pos + 1

END SUBROUTINE move

RECURSIVE SUBROUTINE peano_curve (tmp, original_grid,grid_level_t, &
                                  ori,angle,coordX,coordY,pos)
    INTEGER, DIMENSION(:)                    :: tmp
    INTEGER, DIMENSION(:, :)                 :: original_grid
    INTEGER, INTENT(IN)                      :: grid_level_t, ori
    INTEGER, INTENT(INOUT)                   :: angle, coordX, coordY
    INTEGER                                  :: pos

    INTEGER                                  :: grid_level, ori_tmp

    IF (grid_level_t .GT. 0) THEN
      grid_level = grid_level_t - 1
      ori_tmp = ori
      CALL peano_curve(tmp,original_grid,grid_level,ori_tmp,angle,coordX,coordY,pos)  
      CALL move(angle,coordX,coordY,pos)
      tmp(pos+1) = original_grid(coordX+1,coordY+1)
      ori_tmp = -(ori)
      CALL peano_curve(tmp,original_grid,grid_level,ori_tmp,angle,coordX,coordY,pos)
      CALL move(angle,coordX,coordY,pos)
      tmp(pos+1) = original_grid(coordX+1,coordY+1)      
      ori_tmp=ori
      CALL peano_curve(tmp,original_grid,grid_level,ori_tmp,angle,coordX,coordY,pos)
      angle = angle - (ori*90)
      CALL move(angle,coordX,coordY,pos)
      tmp(pos+1) = original_grid(coordX+1,coordY+1)      
      angle = angle - (ori*90)
      ori_tmp = -(ori)
      CALL peano_curve(tmp,original_grid,grid_level,ori_tmp,angle,coordX,coordY,pos)
      CALL move(angle,coordX,coordY,pos)
      tmp(pos+1) = original_grid(coordX+1,coordY+1)
      ori_tmp=ori
      CALL peano_curve(tmp,original_grid,grid_level,ori_tmp,angle,coordX,coordY,pos)
      CALL move(angle,coordX,coordY,pos)
      tmp(pos+1) = original_grid(coordX+1,coordY+1)
      ori_tmp = -(ori)
      CALL peano_curve(tmp,original_grid,grid_level,ori_tmp,angle,coordX,coordY,pos)
      angle = angle + (ori*90)
      CALL move(angle,coordX,coordY,pos)
      tmp(pos+1) = original_grid(coordX+1,coordY+1)
      angle = angle + (ori*90)
      ori_tmp=ori
      CALL peano_curve(tmp,original_grid,grid_level,ori_tmp,angle,coordX,coordY,pos)
      CALL move(angle,coordX,coordY,pos)
      tmp(pos+1) = original_grid(coordX+1,coordY+1)
      ori_tmp = -(ori)
      CALL peano_curve(tmp,original_grid,grid_level,ori_tmp,angle,coordX,coordY,pos)
      CALL move(angle,coordX,coordY,pos)
      tmp(pos+1) = original_grid(coordX+1,coordY+1)             
      ori_tmp = (ori)
      CALL peano_curve(tmp,original_grid,grid_level,ori_tmp,angle,coordX,coordY,pos)
    ELSE
      angle = 0 
      coordX = 0
      coordY = 0
    ENDIF

END SUBROUTINE peano_curve

FUNCTION is_power_three(numtask) RESULT(ispower)
    INTEGER, INTENT(IN)                      :: numtask
    LOGICAL                                  :: ispower

    INTEGER                                  :: n

   n = numtask

   DO WHILE(MOD(n,9) .EQ. 0)
     n = n/9
   ENDDO

   IF((n .EQ. 1) .OR. (n .EQ. 3)) THEN
     ispower = .TRUE.
   ELSE
     ispower = .FALSE.
   ENDIF

END FUNCTION is_power_three

! *****************************************************************************
!> \brief Create the peano curve for a two dimensional grid 
!> \param mp_comm [output] : handle of the default communicator
!> \param mp_new_comm [output] : handle of the new  communicator
!> note: Algorithm presented on C Code From Graphics Gems II, Academic Press, Inc
!> \par History
!>      4.2012 created [ Christiane Pousa ]
! *****************************************************************************
SUBROUTINE ma_peano_curve(mp_new_comm,mp_comm)   
    INTEGER, INTENT(IN)                   :: mp_comm
    INTEGER, INTENT(out)                  :: mp_new_comm
#if !defined(__parallel)
    mp_new_comm = mp_comm 
#else  
    INTEGER, DIMENSION(:), POINTER           :: ranks_order, tmp
    INTEGER, DIMENSION(:,:), POINTER         :: original_grid             
    INTEGER                                  :: numtask,grid_level,stat
    INTEGER                                  :: grid_dim,pos,i,j,taskid    
    INTEGER                                  :: angle,coordX,ori,coordY    

    ! Get a new MPI rank ordering, reduce comm costs
    
  CALL mp_environ(numtask,taskid,mp_comm)

IF (is_power_three(numtask) .AND. ( numtask .GE. 3**2)) THEN ! the smallest grid for peano
    grid_level = ma_find_grid_level(numtask)
    grid_dim = 3**grid_level


    ALLOCATE(ranks_order(numtask),stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
    ALLOCATE(tmp(grid_dim*grid_dim),stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

    ! the original grid order on CP2K
    ALLOCATE(original_grid(grid_dim,grid_dim))
    pos = 0
    DO i=1, grid_dim 
      DO j=1, grid_dim
        original_grid(i,j) = pos
        pos = pos + 1
      ENDDO 
    ENDDO
    

    tmp(:) = 0
    ori=-1;angle=0;coordX=0;coordY=0;pos=0
    tmp(1) = 0
    CALL peano_curve(tmp,original_grid,grid_level,ori,angle,&
                     coordX,coordY,pos)
    ! Create peano curve
    ranks_order(:) = 0
    
    pos = 0
    DO i=1,SIZE(tmp)
      IF(tmp(i) .LT. numtask) THEN
        ranks_order(tmp(i)+1) = pos
        pos = pos + 1
      ENDIF
    ENDDO

    ! Create the new communicator
    CALL mp_reordering(mp_comm,mp_new_comm,ranks_order)

    stat = 0
    IF(ASSOCIATED(ranks_order)) DEALLOCATE(ranks_order,stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation) 
    IF(ASSOCIATED(tmp)) DEALLOCATE(tmp,stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation) 
    IF(ASSOCIATED(original_grid)) DEALLOCATE(original_grid,stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation) 
  ELSE
      CALL ma_packed(mp_new_comm,mp_comm)
  ENDIF
  
#endif
END SUBROUTINE ma_peano_curve

! *****************************************************************************
!> \brief Define a packed (sequential) distribution for MPI ranks
!> \param mp_comm [output] : handle of the default communicator
!> \param mp_new_comm [output] : handle of the new  communicator
!> \par History
!>      4.2012 created [ Christiane Pousa ]
! *****************************************************************************
SUBROUTINE ma_packed(mp_new_comm,mp_comm)
    INTEGER, INTENT(IN)                   :: mp_comm
    INTEGER, INTENT(out)                  :: mp_new_comm

#if !defined(__parallel)
    mp_new_comm = mp_comm
#else
    INTEGER :: i
    INTEGER, DIMENSION(:), POINTER           :: ranks_order
    INTEGER                                  :: numtask,taskid,stat

    ! Get a new MPI rank ordering, reduce comm costs
    CALL mp_environ(numtask,taskid,mp_comm)   
    ALLOCATE(ranks_order(numtask),stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation) 

   DO i=1, numtask
     ranks_order(i) = i - 1
   ENDDO 

    ! Create the new communicator
    CALL mp_reordering(mp_comm,mp_new_comm,ranks_order)

    IF(ASSOCIATED(ranks_order)) DEALLOCATE(ranks_order,stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
#endif    
END SUBROUTINE ma_packed


! *****************************************************************************
!> \brief Define a round-robin distribution for MPI ranks, the round is made
!> \brief using the number of allocated machines and cores
!> \param mp_comm [output] : handle of the default communicator
!> \param mp_new_comm [output] : handle of the new  communicator
!> \par History
!>      4.2012 created [ Christiane Pousa ]
! *****************************************************************************
SUBROUTINE ma_round_robin(mp_new_comm,mp_comm)
    INTEGER, INTENT(out)                     :: mp_new_comm
    INTEGER, INTENT(IN)                      :: mp_comm

#if defined(__parallel)
    INTEGER                                  :: unit_nr
#if defined(__HWLOC)
    INTEGER                                  :: i, ind, n, ncores, nnodes, &
                                                numtask, rounds, stat, taskid
    INTEGER, DIMENSION(:), POINTER           :: ranks_order
#endif
#endif

#if ! defined (__parallel)
   mp_new_comm = mp_comm
#else
#if defined __HWLOC
    ! get number of machines and number of cores per machine
    nnodes = ma_get_nmachines(mp_comm)
    ncores = ma_get_ncores()
 
    IF (nnodes .GT. 0) THEN
      ! Get a new MPI rank ordering, reduce comm costs
       CALL mp_environ(numtask,taskid,mp_comm)
       ALLOCATE(ranks_order(numtask),stat=stat)
       IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

     ranks_order(:) = 0
     ind = 1
     rounds = CEILING(REAL((numtask/ncores),KIND=dp))

     IF(rounds .EQ. 0)THEN
        DO i=0,numtask-1
          ranks_order(ind) = i
          ind = ind + 1
        ENDDO
      ELSE
        i = 0
        DO ind=1, numtask
          ranks_order(ind) =  MOD((ind-1),nnodes)*ncores+i
          IF (MOD(ind,nnodes) == 0 ) i = i + 1
        ENDDO
      ENDIF

      ! Create the new communicator
      CALL mp_reordering(mp_comm,mp_new_comm,ranks_order)

      stat = 0
      IF(ASSOCIATED(ranks_order)) DEALLOCATE(ranks_order,stat=stat)
      IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
    ELSE
       unit_nr = 6
       WRITE(unit_nr,'(T2,A,A)') "WARNING: round-robin not", &
            " applied"
       mp_new_comm = mp_comm
    ENDIF
#else
   unit_nr = 6
   WRITE(unit_nr,'(T2,A,A)') "WARNING: You need hwloc library", &
            " to use round-robin strategy"
    mp_new_comm = mp_comm
#endif
#endif

END SUBROUTINE ma_round_robin


! *****************************************************************************
!> \brief Create the hilbert curve for a two dimensional grid
!> \param mp_comm [output] : handle of the default communicator
!> \param mp_new_comm [output] : handle of the new  communicator
!> \note  Based on the algorithm proposed Liu and SCHRACK, 1996
!>        Enconding and Decoding the Hilbert Order
!> \par History
!>      2.2012 created [ Christiane Pousa ]
! *****************************************************************************
SUBROUTINE ma_hilbert_bits(ranks_order,ntask,mp_comm)
    INTEGER, INTENT(IN)                   :: mp_comm
    INTEGER, INTENT(IN), OPTIONAL         :: ntask
    INTEGER, DIMENSION(:)                 :: ranks_order
#if defined(__parallel)
    INTEGER, DIMENSION(:),POINTER          :: tmp_order
    INTEGER                                :: numtask,taskid,stat

    INTEGER                                :: mask,notx,noty,heven,&
                                              hodd,temp,v1,v0,k,i
    INTEGER                                :: grid_dim,pos,x,y,z

    IF(.NOT. PRESENT(ntask)) THEN
     ! Get a new MPI rank ordering, reduce comm costs
     CALL mp_environ(numtask,taskid,mp_comm)
     grid_dim = SQRT(REAL(ma_find_min_square_grid(numtask),KIND=dp))
     ALLOCATE(tmp_order(grid_dim*grid_dim),stat=stat)
     IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
    ELSE
     numtask = ntask
     grid_dim = SQRT(REAL(ma_find_max_square_grid(numtask),KIND=dp))
     ALLOCATE(tmp_order(grid_dim*grid_dim),stat=stat)
     IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
    ENDIF
 
    ranks_order(:) = 0

    ! Create hilbert curve - encoding from (x,y) to a hilbert point
    pos = 0
    DO x=0, grid_dim-1
     DO y=0, grid_dim-1
      mask = ISHFT(1,grid_dim)-1
      heven = IEOR(x,y)
      temp = NOT(x)
      notx = IAND(temp,mask)
      temp = NOT(y)
      noty = IAND(temp,mask)
      temp = IEOR(notx,y)
      v0 = 0
      v1 = 0

      DO k=0, grid_dim-1
       v1 = ISHFT(IOR(IAND(v1,heven),IAND(IEOR(v0,noty),temp)),-1)

       v0 = ISHFT(IOR(IAND(v0,IEOR(v1,notx)),IAND(NOT(v0),IEOR(v1,noty))),-1)
      ENDDO
      hodd = IOR(IAND(NOT(v0),IEOR(v1,x)),IAND(v0,IEOR(v1,noty)))

      z=0
      DO i=0, int_size*8-1
       temp = IOR(ISHFT(IAND(heven,ISHFT(1,i)),i),ISHFT(IAND(hodd,ISHFT(1,i)),i+1))
       z = IOR(z,temp)
      ENDDO
        tmp_order(z+1) = pos
        pos = pos + 1
     ENDDO
    ENDDO

    pos = 0 
    DO k=1,SIZE(tmp_order)
      IF(tmp_order(k) .LT. numtask) THEN
        ranks_order(tmp_order(k)+1) = pos
        pos = pos + 1
      ENDIF 
    ENDDO

    stat = 0
    IF(ASSOCIATED(tmp_order)) DEALLOCATE(tmp_order,stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
#endif
END SUBROUTINE ma_hilbert_bits

! *****************************************************************************
!> \brief Create the hilbert curve for a two dimensional grid
!> \param mp_comm [output] : handle of the default communicator
!> \param mp_new_comm [output] : handle of the new  communicator
!> \note  Based on the algorithm proposed Liu and SCHRACK, 1996
!>        Enconding and Decoding the Hilbert Order
!> \par History
!>      2.2012 created [ Christiane Pousa ]
! *****************************************************************************
SUBROUTINE ma_hilbert_curve(mp_new_comm,mp_comm)
    INTEGER, INTENT(IN)                   :: mp_comm
    INTEGER, INTENT(out)                  :: mp_new_comm

#if !defined(__parallel)
    mp_new_comm = mp_comm 
#else
    INTEGER, DIMENSION(:), POINTER           :: ranks_order
    INTEGER                                  :: numtask,taskid,stat


    ! Get a new MPI rank ordering, reduce comm costs
    CALL mp_environ(numtask,taskid,mp_comm)

  IF(square_grid(numtask) .OR. (IAND(numtask,numtask-1) .EQ. 0)) THEN
  
    ALLOCATE(ranks_order(numtask),stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
 
    ranks_order(:) = 0

    CALL ma_hilbert_bits(ranks_order,mp_comm=mp_comm)

    ! Create the new communicator
    CALL mp_reordering(mp_comm,mp_new_comm,ranks_order)

    IF(ASSOCIATED(ranks_order)) DEALLOCATE(ranks_order,stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
  ELSE
    mp_new_comm = mp_comm
  ENDIF
#endif

END SUBROUTINE ma_hilbert_curve

RECURSIVE SUBROUTINE snake_curve (myrank, pos, order,grid_dim,neighbor,done)
    INTEGER, INTENT(IN)                      :: grid_dim
    INTEGER                                  :: pos, myrank
    INTEGER, DIMENSION(:), POINTER           :: order
    LOGICAL, DIMENSION(:), POINTER           :: neighbor, done

#if defined(__parallel)
    INTEGER                                  :: nup,nright
   
    IF ( pos .LT. SIZE(order) .AND. .NOT.(done(myrank))) THEN
     IF(.NOT.neighbor(pos)) THEN
         order(pos) = myrank
         pos = pos + 1
         neighbor(pos) = .TRUE.      
     ENDIF
     nup = myrank + grid_dim
     IF(nup .LT. SIZE(order) .AND. .NOT.neighbor(nup)) THEN  
         order(pos) = nup
         pos = pos + 1
         neighbor(nup) = .TRUE.        
     ENDIF
     nright = myrank + 1
     IF(nright .LT. SIZE(order) .AND. .NOT.neighbor(nright)) THEN
         order(pos) = nright
         pos = pos + 1
         neighbor(nright) = .TRUE.
     ENDIF

     CALL snake_curve(nup,pos,order,grid_dim,neighbor,done)
     CALL snake_curve(nright,pos,order,grid_dim,neighbor,done)     
    ENDIF
#endif    

END SUBROUTINE snake_curve

! *****************************************************************************
!> \brief Create the interleaved rows order for a two dimensional grid
!> \param mp_comm [output] : handle of the default communicator
!> \param mp_new_comm [output] : handle of the new  communicator
!> \par History
!>      2.2012 created [ Christiane Pousa ]
! *****************************************************************************
SUBROUTINE ma_snake_curve (mp_new_comm,mp_comm)
    INTEGER, INTENT(IN)                   :: mp_comm
    INTEGER, INTENT(out)                  :: mp_new_comm

#if !defined(__parallel)
    mp_new_comm = mp_comm
#else
    INTEGER, DIMENSION(:), POINTER           :: ranks_order,order
    LOGICAL, DIMENSION(:), POINTER           :: neighbor, done
    INTEGER                                  :: numtask,taskid,&
                                                stat, i, pos, grid_dim

    ! Get a new MPI rank ordering, reduce comm costs
    CALL mp_environ(numtask,taskid,mp_comm)
  IF (numtask .GE. 4) THEN
    grid_dim = ma_find_min_square_grid(numtask)
    ALLOCATE(ranks_order(numtask),stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
    ALLOCATE(neighbor(grid_dim*grid_dim),stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
    ALLOCATE(done(grid_dim*grid_dim),stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)    
    ALLOCATE(order(grid_dim*grid_dim),stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)    

    ranks_order(:) = 0
    order(:) = 0
    neighbor(:) = .FALSE.
    done(:) = .FALSE.

    ! Create the interleaved rows order in a recursive way
    ! Starts with rank 0
    CALL snake_curve (0,1,order,grid_dim,neighbor,done)  
    pos = 1
    DO i=1, SIZE(order) - 1
     IF (order(i) .LT. SIZE(ranks_order) ) THEN
      ranks_order(pos) = order(i)
      pos = pos + 1
     ENDIF
    ENDDO
   
    ! Create the new communicator
    CALL mp_reordering(mp_comm,mp_new_comm,ranks_order)

    IF(ASSOCIATED(ranks_order)) DEALLOCATE(ranks_order,stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
    IF(ASSOCIATED(order)) DEALLOCATE(order,stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
    IF(ASSOCIATED(neighbor)) DEALLOCATE(neighbor,stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
    IF(ASSOCIATED(done)) DEALLOCATE(done,stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
  ELSE 
   CALL ma_packed(mp_new_comm,mp_comm) 
  ENDIF

#endif
END SUBROUTINE ma_snake_curve

SUBROUTINE rotate(s,x,y,rx,ry)
    INTEGER                                  :: s, x, y, rx, ry

    INTEGER                                  :: tmp

 IF (ry .EQ. 0) THEN
    IF (rx .EQ. 1) THEN
          x = s-1 - x
          y = s-1 - y
    ENDIF
          tmp = x
          x = y
          y = tmp
 ENDIF


END SUBROUTINE rotate

FUNCTION distance(grid_dim,countx,county) RESULT(dist)
    INTEGER                                  :: grid_dim, countx, county, dist

    INTEGER                                  :: rx, ry, s, tmp, x, y

 dist=0
 s = grid_dim/2
 x = countx
 y = county
 DO WHILE (s .GT. 0)
      IF (IAND(x,s) .GT. 0) THEN
       rx = 1
      ELSE
       rx = 0
      ENDIF
      IF (IAND(y,s) .GT. 0) THEN
       ry = 1
      ELSE
       ry = 0
      ENDIF
      tmp = 3*rx
      tmp = IEOR(tmp,ry)
      dist = dist + (s*s*tmp)
      CALL rotate(s,x,y,rx,ry)
      s = s/2
 ENDDO

END FUNCTION distance

SUBROUTINE ma_hilbert_order(tmp_order,grid_dim)

  INTEGER, DIMENSION(:)                  :: tmp_order
  INTEGER                                :: grid_dim

#if defined(__parallel)
  INTEGER                                :: pos, countx, county

  tmp_order(:) = 0  
  ! Create hilbert curve - encoding from (x,y) to a hilbert point

  pos = 1
  countx = 0
  DO WHILE (countx .LE. grid_dim-1)
   county=0
   DO WHILE (county .LE. grid_dim-1)
     tmp_order(pos) = distance(grid_dim,countx,county)
     pos = pos + 1 
     county = county + 1
  ENDDO
    countx = countx + 1
  ENDDO  

#endif

END SUBROUTINE ma_hilbert_order 

! *****************************************************************************
!> \brief Create a hilbert order for all sub-grids within a grid
!> \param first_grid [input] : grid dimension for the biggest grid
!> \param rest_grid [input] : dim for the grid rest
!> \param numtask [input] : number of tasks
!> return has_2d one if has more 2d gris and zero if not
!>      7.2012 created [ Christiane Pousa ]
! *****************************************************************************
RECURSIVE SUBROUTINE ma_hilbert_recursive(ranks_order,ncol,nrow,numtask,& 
                                          current_numtask,current_nrow, &
                                          current_ncol, start_grid,&
                                          working_col,mp_comm)
 INTEGER, DIMENSION(:)                 :: ranks_order
 INTEGER                               :: ncol,nrow,numtask,current_numtask,&
                                          current_nrow,current_ncol,start_grid,&
                                          mp_comm
 LOGICAL                               :: working_col

#if defined(__parallel)
 INTEGER                               :: grid_dim,count_itens,index_order, &
                                          row_counter,col_counter,col_childs,&
                                          col_rest,i,last_grid,row_child,&
                                          row_rest,ntasks,taskid,stat,&
                                          local_start_grid, local_current_numtask,&
                                          local_numtask,local_current_ncol
 INTEGER, DIMENSION(:), POINTER        :: tmp_order, real_grid

 IF( (current_numtask .GE. 4) .AND. (numtask .GT. 0))THEN 
  IF (working_col) THEN

   grid_dim = SQRT(REAL(ma_find_max_square_grid(current_numtask),KIND=dp))
  
   ALLOCATE(tmp_order(grid_dim*grid_dim),stat=stat)
   IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
   tmp_order(:) = 0
  
   CALL mp_environ(local_numtask,taskid,mp_comm)

  !Select the curve direction
   IF ( local_numtask .EQ. numtask) THEN  
    CALL ma_hilbert_bits(tmp_order,numtask,mp_comm)
   ELSE
    CALL ma_hilbert_order(tmp_order,grid_dim)
   ENDIF

   ALLOCATE(real_grid(grid_dim*grid_dim),stat=stat)
   IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
   real_grid(:) = 0   

   IF ( local_numtask .EQ. numtask) THEN
     !create the real sub-grid
     count_itens = 1
     DO row_counter=0, grid_dim-1
      DO col_counter=1, grid_dim
       real_grid(count_itens) = row_counter*nrow+start_grid+col_counter-1
       count_itens = count_itens + 1
      ENDDO
     ENDDO 
   ELSE ! revert grid ranks
     !create the real sub-grid
     count_itens = 1
     DO col_counter=1, grid_dim
      DO row_counter=0, grid_dim-1
       real_grid(count_itens) = row_counter*nrow+start_grid+col_counter-1
       count_itens = count_itens + 1
      ENDDO
     ENDDO
   ENDIF 

   !include the ordered ranks
   count_itens = 0
   row_counter = 0
   col_counter = 1 
   DO WHILE ( count_itens .LT. SIZE(tmp_order))
     index_order = row_counter*nrow + start_grid + col_counter 
     ranks_order(index_order) = real_grid(tmp_order(count_itens+1)+1)
     count_itens = count_itens + 1 
     col_counter = col_counter + 1
     IF ( col_counter .GT. grid_dim) THEN
        col_counter = 1
        row_counter = row_counter + 1
     ENDIF
   ENDDO

   IF(ASSOCIATED(tmp_order)) DEALLOCATE(tmp_order,stat=stat)
   IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)  
 
   IF(ASSOCIATED(real_grid)) DEALLOCATE(real_grid,stat=stat)
   IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

   last_grid = grid_dim
   !updating current_numtask
   local_current_numtask = (current_ncol - last_grid) * (current_ncol - last_grid)
   
   IF ( local_current_numtask .NE. 0 ) THEN
     grid_dim = INT(CEILING(SQRT(REAL(ma_find_max_square_grid(local_current_numtask),KIND=dp))))
     col_childs = INT(current_nrow/grid_dim)
     col_rest = MOD(nrow,grid_dim) 
     local_numtask = numtask
     
  
     DO i=1, col_childs
     local_start_grid = ((i-1)*nrow*grid_dim)+last_grid
      IF( i .EQ. 1) THEN
         local_numtask = local_numtask-i*last_grid*last_grid
      ELSE
         local_numtask = local_numtask-grid_dim*grid_dim
      ENDIF
      !                        local_start_grid,working_col
      CALL ma_hilbert_recursive(ranks_order,ncol,nrow,numtask=local_numtask,&
                                current_numtask=local_current_numtask, current_nrow=grid_dim, & 
                                current_ncol=grid_dim,start_grid=local_start_grid,&
                                working_col=.TRUE.,mp_comm=mp_comm)
     ENDDO
   ENDIF
  ELSE ! perform rows

   local_current_numtask = current_numtask
   grid_dim = SQRT(REAL(ma_find_max_square_grid(local_current_numtask),KIND=dp))

   ALLOCATE(tmp_order(grid_dim*grid_dim),stat=stat)
   IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
   tmp_order(:) = 0

   CALL mp_environ(local_numtask,taskid,mp_comm)

   IF ( local_numtask .EQ. numtask) THEN
    CALL ma_hilbert_bits(tmp_order,numtask,mp_comm)
   ELSE
    CALL ma_hilbert_order(tmp_order,grid_dim)
   ENDIF

   ALLOCATE(real_grid(grid_dim*grid_dim),stat=stat)
   IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
   real_grid(:) = 0

   IF ( local_numtask .EQ. numtask) THEN
    !create the real sub-grid
     count_itens = 1 
     DO row_counter=0, grid_dim-1
      DO col_counter=1, grid_dim
       real_grid(count_itens) = row_counter*nrow+start_grid+col_counter-1
       count_itens = count_itens + 1
      ENDDO
     ENDDO

   ELSE
    !create the real sub-grid
     count_itens = 1
     DO col_counter=grid_dim,1,-1
      DO row_counter=grid_dim-1,0, -1
       real_grid(count_itens) = row_counter*nrow+start_grid+col_counter-1
       count_itens = count_itens + 1
      ENDDO
     ENDDO
   ENDIF

   !include the ordered ranks
   count_itens = 0
   row_counter = 0
   col_counter = 1
   DO WHILE ( count_itens .LT. SIZE(tmp_order))
     index_order = row_counter*nrow + start_grid + col_counter
     ranks_order(index_order) = real_grid(tmp_order(count_itens+1)+1)
     count_itens = count_itens + 1
     col_counter = col_counter + 1
     IF ( col_counter .GT. grid_dim) THEN
      col_counter = 1
      row_counter = row_counter + 1
     ENDIF
   ENDDO

   IF(ASSOCIATED(tmp_order)) DEALLOCATE(tmp_order,stat=stat)
   IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

   IF(ASSOCIATED(real_grid)) DEALLOCATE(real_grid,stat=stat)
   IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)


   last_grid = grid_dim
   local_current_numtask = (current_nrow - last_grid) * (current_nrow - last_grid)
   
   IF ( local_current_numtask .NE. 0 ) THEN
      grid_dim = SQRT(REAL(ma_find_max_square_grid(local_current_numtask),KIND=dp))
      row_child = INT(last_grid/grid_dim)
      row_rest = MOD(local_current_numtask+last_grid,grid_dim)
      local_current_ncol = current_ncol
      local_numtask = numtask
  
 
      DO i=1, row_child
       local_start_grid = last_grid*local_current_ncol+(i-1)*grid_dim
       IF( i== 1) THEN
          local_numtask = local_numtask-last_grid*last_grid
       ELSE
          local_numtask = local_numtask-grid_dim*grid_dim
       ENDIF

        CALL ma_hilbert_recursive(ranks_order,ncol,nrow,numtask=local_numtask,&
                                   current_numtask=local_current_numtask, current_nrow=grid_dim, &
                                   current_ncol=local_current_ncol,start_grid=local_start_grid,&
                                   working_col=.FALSE.,mp_comm=mp_comm)
      
   
    
     ENDDO
    ENDIF
   ENDIF
 ELSE IF (numtask .GT. 0)THEN
    ! just include numbers in the ordering
   CALL mp_environ(ntasks,taskid,mp_comm)
   i = start_grid + 1
   DO count_itens=1, current_numtask
    IF ( ranks_order(i) .EQ. 0) ranks_order(i) = i - 1
    i = i + 1
   ENDDO
ENDIF
#endif

END SUBROUTINE ma_hilbert_recursive

FUNCTION is_square (grid_dim) RESULT (square)
    INTEGER                                  :: grid_dim
    LOGICAL                                  :: square

    INTEGER                                  :: i, j

 square = .FALSE. 
 
!Assuming the biggestes number of cores in a system os
! 16777216
 j = 2
 i = 2**j
 DO WHILE ( (grid_dim .GE. i) .AND. (i .LT. 16777216))
   IF ( grid_dim .EQ. i) square = .TRUE.
   j = j +1
   i = 2**j
 ENDDO

END FUNCTION is_square

! *****************************************************************************
!> \brief Verify if there are any other hilbert grids
!> \param ncol [input] : number of columns
!> \param nrows [input] : number of rows
!> \param numtask [input] : number of tasks
!> return has_2d one if has more 2d gris and zero if not
!>      7.2012 created [ Christiane Pousa ]
! *****************************************************************************
FUNCTION has_more_2dgrid(ncol,nrow,numtask) RESULT(has_2d)

    INTEGER, INTENT(IN)                      :: ncol, nrow, numtask
    INTEGER                                  :: has_2d

#if defined(__parallel)
    INTEGER                                  :: grid_dim, rest, rest_col
#endif

 has_2d = 0 

#if defined(__parallel)
 grid_dim = SQRT(REAL(ma_find_max_square_grid(numtask),KIND=dp))
 rest =  nrow - grid_dim
 rest_col = ncol - grid_dim
 !Now compute the remaining grids
 IF ( rest .EQ. rest_col ) THEN
   grid_dim = SQRT(REAL(ma_find_max_square_grid(rest*rest),KIND=dp))
   IF (is_square(grid_dim)) has_2d = 1
 ENDIF
#endif

END FUNCTION has_more_2dgrid

! *****************************************************************************
!> \brief Verify if the number of tasks is square
!> \param numtask [input] : number of tasks
!> return square true if it is square
!>      10.2012 created [ Christiane Pousa ]
! *****************************************************************************
FUNCTION square_grid(numtask) RESULT(square)
    INTEGER, INTENT(IN)                      :: numtask
    LOGICAL                                  :: square

    INTEGER                                  :: grid_dim

  grid_dim = INT(SQRT(REAL(numtask,KIND=dp)))
  IF (grid_dim*grid_dim .EQ. numtask) THEN
    square = .TRUE.
  ELSE
    square = .FALSE.
  ENDIF
   
END FUNCTION  square_grid

! *****************************************************************************
!> \brief Create the hilbert curve for a two dimensional grid
!> \param mp_comm [input] : handle of the default communicator
!> \param mp_new_comm [output] : handle of the new  communicator
!> \par History
!>      6.2012 created [ Christiane Pousa ]
! *****************************************************************************
SUBROUTINE ma_hilbert_peano(mp_new_comm,mp_comm)
  INTEGER, INTENT(IN)                   :: mp_comm
  INTEGER, INTENT(OUT)                  :: mp_new_comm
#if !defined(__parallel)
  mp_new_comm = mp_comm
#else
  INTEGER, DIMENSION(:), POINTER           :: ranks_order, tmp_order, &
                                              coord_order, tmp
  INTEGER                                  :: numtask,taskid,stat, &
                                              h_grid_dim, grid_dim
  INTEGER                                  :: ipe,jpe, index_order, &
                                              rest_col,ncol, nrow, &
                                              end_col,start_row,end_row, &
                                              done, rest_row,kpe, has_2dgrid

  ! Get a new MPI rank ordering, reduce comm costs
  CALL mp_environ(numtask,taskid,mp_comm)

  IF(square_grid(numtask) .OR. (IAND(numtask,numtask-1) .EQ. 0)) THEN  

    ALLOCATE(ranks_order(numtask),stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
    ranks_order(:) = 0
  
    CALL ma_2dgrid_dimensions(ncol,nrow,numtask)
    has_2dgrid = has_more_2dgrid(ncol,nrow,numtask)
  
    IF( numtask .GT. 4) THEN  !smallest grid
     h_grid_dim = SQRT(REAL(ma_find_max_square_grid(numtask),KIND=dp))
  
     ! perfect hilbert grid
     IF ( numtask .EQ. (h_grid_dim*h_grid_dim)) THEN 
      CALL ma_hilbert_order(ranks_order,h_grid_dim)
     ! there are other sub grids
     ELSE IF ( has_2dgrid .GT. 0 ) THEN
      CALL ma_hilbert_recursive(ranks_order,ncol=ncol,nrow=nrow,numtask=numtask,&
                                current_numtask=numtask, current_nrow=nrow, &
                                current_ncol=ncol, start_grid=0, &
                                working_col=.TRUE.,mp_comm=mp_comm)
  
      CALL ma_hilbert_recursive(ranks_order,ncol=ncol,nrow=nrow,numtask=numtask,&
                                current_numtask=numtask, current_nrow=nrow, &
                                current_ncol=ncol, start_grid=0, &
                                working_col=.FALSE.,mp_comm=mp_comm)
  
     DO ipe=2,numtask
      IF( ranks_order(ipe) == 0 ) ranks_order(ipe) = ipe-1
     ENDDO
  
     IF(taskid==0) THEN
       ALLOCATE(tmp(numtask),stat=stat)
       IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)   
       DO ipe=1,numtask
        tmp(ranks_order(ipe)+1) = ipe -1
       ENDDO
     ENDIF
     !no perfect hilbert grids, use peano
     ELSE
     ALLOCATE(tmp_order(h_grid_dim*h_grid_dim),stat=stat)
     IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)   
     ALLOCATE(coord_order(h_grid_dim*h_grid_dim),stat=stat)
     IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)   
  
      coord_order(:) = 0
      tmp_order(:) = 0 
   
      CALL ma_hilbert_order(coord_order,h_grid_dim)
      
      done = 1
      DO ipe=1,h_grid_dim
       DO jpe=0,h_grid_dim-1
         tmp_order(jpe*h_grid_dim+ipe) = coord_order(done)
         done = done + 1
       ENDDO
      ENDDO
  
      coord_order(:) = tmp_order(:)
   
      DEALLOCATE(tmp_order,stat=stat)
      IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
      ALLOCATE(tmp_order(numtask),stat=stat)
      IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
      tmp_order(:) = 0    
  
      !Now include the column ranks that were not in the Hilbert grid
      rest_col = ncol-h_grid_dim
      rest_row = nrow-h_grid_dim
      index_order = 1
      grid_dim = SQRT(REAL(numtask,KIND=dp))
  
      DO ipe=0,nrow-1
       DO jpe=ipe*ncol+h_grid_dim, ipe*ncol+h_grid_dim+rest_col-1
        tmp_order(index_order) = jpe
        index_order = index_order + 1
       ENDDO
      ENDDO
  
      end_col = index_order
      start_row = index_order 
  
      !Now include de row ranks that were not in the Hilbert grid
      DO ipe=0,rest_row-1
  !      DO jpe=numtask-rest_col-(ipe*nrow)-1,ncol*h_grid_dim, -1
        DO jpe=numtask-rest_col-(ipe*nrow)-1,(numtask-rest_col-(ipe*nrow))-h_grid_dim, -1
          tmp_order(index_order) = jpe
          index_order = index_order + 1
        ENDDO
      ENDDO
  
      end_row = index_order
  
      !And now copy all them to the rank_order array
      index_order = 1
      jpe = 0 
      kpe = 1    
      DO ipe=1,numtask
       IF ( (ipe .GT. jpe*ncol) .AND. (ipe .LT. jpe*ncol+h_grid_dim+1) .AND. &
             (index_order .LE. SIZE(coord_order))) THEN
           ranks_order(ipe) = coord_order(index_order)
           index_order = index_order + 1
        ENDIF
       IF (ipe == h_grid_dim+jpe*ncol) jpe = jpe + 1
      ENDDO
  
      jpe = 0
      kpe = 1        
      index_order = index_order - 1 
      DO ipe=1,numtask+1
       IF ( (ipe .GE. jpe*ncol+nrow) .AND. (ipe .LT. jpe*ncol+nrow+rest_col)) THEN
           ranks_order(tmp_order(kpe)+1) = index_order 
           kpe = kpe + 1
           index_order = index_order + 1
       ENDIF
       IF (ipe == (jpe+1)*ncol+1) jpe = jpe + 1
      ENDDO
  
      IF ( rest_row == 1) THEN
         DO ipe=kpe,kpe+h_grid_dim-1
          ranks_order(tmp_order(ipe)+1) = index_order
          index_order = index_order + 1
         ENDDO
      ELSE
        
        DO WHILE ( index_order .LT. numtask)
          IF (kpe .LT. SIZE(tmp_order) .AND. tmp_order(kpe) .NE. 0) THEN
            ranks_order(tmp_order(kpe)+1) = index_order
            index_order = index_order + 1           
            tmp_order(kpe) = 0
          ENDIF
          kpe = kpe + h_grid_dim
          IF (kpe .LT. SIZE(tmp_order) .AND. tmp_order(kpe) .NE. 0) THEN
            ranks_order(tmp_order(kpe)+1) = index_order
            index_order = index_order + 1
            tmp_order(kpe) = 0
          ENDIF
          kpe = kpe + 1
          IF (kpe .LT. SIZE(tmp_order) .AND. tmp_order(kpe) .NE. 0) THEN
            ranks_order(tmp_order(kpe)+1) = index_order
            index_order = index_order + 1
            tmp_order(kpe) = 0
          ENDIF
          kpe = kpe - h_grid_dim
          IF (kpe .LT. SIZE(tmp_order) .AND. tmp_order(kpe) .NE. 0) THEN
            ranks_order(tmp_order(kpe)+1) = index_order
            index_order = index_order + 1
            tmp_order(kpe) = 0
          ENDIF
          kpe = kpe + 1      
        ENDDO
      ENDIF
  
     ENDIF
    ! Create the new communicator
     CALL mp_reordering(mp_comm,mp_new_comm,ranks_order)
    ELSE
     CALL ma_packed(mp_new_comm,mp_comm)   
    ENDIF  
  ELSE ! for now we just copy the old ordering
    mp_new_comm = mp_comm
  ENDIF !perfect square and power of two

  IF(ASSOCIATED(tmp_order))DEALLOCATE(tmp_order,stat=stat)
  IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
 
  IF(ASSOCIATED(coord_order))DEALLOCATE(coord_order,stat=stat)
  IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)  

  IF(ASSOCIATED(ranks_order))DEALLOCATE(ranks_order,stat=stat)
  IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
#endif

END SUBROUTINE ma_hilbert_peano


! *****************************************************************************
!> \brief Create the hilbert curve for a two dimensional grid
!> \param mp_comm [input] : handle of the default communicator
!> \param mp_new_comm [output] : handle of the new  communicator
!> \param order [input] : the current MPI ordering
!> \param kernel [input]: for each kernel the switch should be applied 
!> \par History
!>      6.2012 created [ Christiane Pousa ]
! *****************************************************************************
SUBROUTINE ma_switch(mp_new_comm,mp_comm,order,kernel)
  INTEGER, INTENT(OUT)                                   :: mp_new_comm
  INTEGER, INTENT(IN)                                    :: mp_comm
  INTEGER, DIMENSION(:), POINTER, INTENT(IN), OPTIONAL   :: order
  INTEGER, INTENT(IN), OPTIONAL                          :: kernel

#if defined(__parallel) 
#if defined __GEMINI || __SEASTAR || __NET
  INTEGER, DIMENSION(:), POINTER                         :: ordering
#endif
#endif

#if ! defined(__parallel) 
  mp_new_comm = mp_comm
#else
#if defined __GEMINI || __SEASTAR || __NET
   INTEGER                                    :: numtask, taskid,i, &
                                                 ipe,jpe, tmp_order, &
                                                 not_exchange, stat, jump,&
                                                 exchange
   REAL, DIMENSION(3)                         :: contention
   INTEGER, DIMENSION(:), POINTER             :: tmp_topology, tmp_topo
   REAL                                       :: old_contention, new_contention 

   ! Get a new MPI rank ordering, reduce comm costs
   CALL mp_environ(numtask,taskid,mp_comm)

#ifdef __DEBUG_HWTOPO
   IF(taskid==0) WRITE(*,*) 'Starting switch - each MPI rank is exchanged with its 4 neighbors' 
#endif

   not_exchange = 0
   exchange = 0
   jump = numtask-1

   ALLOCATE(ordering(numtask),stat=stat)
   IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
   ordering(:) = 0

   IF ( .NOT. PRESENT(order)) THEN
     DO i=1, numtask
        ordering(i) = i - 1
     ENDDO
   ELSE
    ordering(:) = order(:)

    ALLOCATE(tmp_topo(numtask),stat=stat)
    IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

    DO ipe=1, numtask
      tmp_topo(ordering(ipe)+1) = allocated_topology(ipe)
    ENDDO  
     
    allocated_topology(:) = tmp_topo(:)
 
   ENDIF
 
   !In the case that there is not a communication graph, we assume 
   ! that all tasks communicate with each other. 
   IF ( .NOT. PRESENT(kernel)) THEN
#ifdef __DEBUG_HWTOPO
   IF(taskid==0) WRITE(*,*) 'Complete graph'
#endif
    CALL ma_comm_graph (mp_comm,complete_graph)
   ELSE
#ifdef __DEBUG_HWTOPO
   IF(taskid==0) WRITE(*,*) 'Cannon graph'
#endif
   !Create the comm graph
    CALL ma_comm_graph (mp_comm,kernel)
   ENDIF
 
   !Get the current contention
   CALL ma_compute_contention(contention,mp_comm)
   old_contention = SUM(contention)

   IF (kernel .EQ. cannon_graph) jump = 4 ! 4 neighbors
   

   !Make the switch between each pair of MPI and compute the new contention
   !The computed contention take into account the communication graph 
   !in the ma_topology

IF(old_contention .GT. 0) THEN
   new_contention = 0
   DO ipe=1,numtask
     DO jpe=ipe, ipe+jump
       tmp_order = ordering(jpe)
       ordering(jpe) = ordering(jpe+1)        
       ordering(jpe+1) = tmp_order

       tmp_order = allocated_topology(jpe)
       allocated_topology(jpe) = allocated_topology(jpe+1)
       allocated_topology(jpe+1) = tmp_order    

#ifdef __DEBUG_HWTOPO
   IF(taskid==0) WRITE(*,*) 'Exchanged ',ordering(jpe+1), ' with ', ordering(jpe)
#endif

       CALL ma_compute_contention(contention,mp_comm)
       new_contention = SUM(contention)

#ifdef __DEBUG_HWTOPO
   IF(taskid==0) WRITE(*,*) 'New contention ', new_contention
#endif

       !we just change if the new contention is at least 10% better than 
       IF ( new_contention .LT. old_contention) THEN
          old_contention = new_contention
          exchange = exchange + 1
       ELSE
         tmp_order = ordering(jpe+1)
         ordering(jpe+1) = ordering(jpe)
         ordering(jpe) = tmp_order

         tmp_order = allocated_topology(jpe+1)
         allocated_topology(jpe+1) = allocated_topology(jpe)
         allocated_topology(jpe) = tmp_order

         not_exchange = not_exchange + 1
#ifdef __DEBUG_HWTOPO
   IF(taskid==0) WRITE(*,*) 'Exchanged back ',ordering(jpe), ' with ', ordering(jpe+1)
#endif
      ENDIF
     ENDDO

       IF (( ipe .GE. numtask/2) .OR. (not_exchange .GE. numtask/2) ) EXIT 
   ENDDO
ENDIF

#ifdef __DEBUG_HWTOPO
   IF(taskid==0) WRITE(*,*) 'Finishing switch - have changed ', exchange
#endif

 ! Create the new communicator
   CALL mp_reordering(mp_comm,mp_new_comm,ordering)

   IF(ASSOCIATED(tmp_topo)) DEALLOCATE(tmp_topo,stat=stat)
   IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

   IF(ASSOCIATED(ordering)) DEALLOCATE(ordering,stat=stat)
   IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

#endif
#endif

END SUBROUTINE ma_switch

! *****************************************************************************
!> \brief Create the near-optimal ordering for cannon communication
!> \param mp_comm [input] : handle of the default communicator
!> \param mp_new_comm [output] : handle of the new  communicator
!> \par History
!>      6.2012 created [ Christiane Pousa ]
! *****************************************************************************
SUBROUTINE ma_cannon(mp_new_comm,mp_comm)
  INTEGER, INTENT(IN)                           :: mp_comm
  INTEGER, INTENT(OUT)                          :: mp_new_comm

#if defined(__parallel)
#if defined  __GEMINI || __SEASTAR || __NET
  INTEGER, DIMENSION(:), POINTER           :: ranks_order, tmp_order, &
                                              coord_order, orderX, orderY, &
                                              orderZ
  REAL, DIMENSION(3)                       :: contention
  INTEGER                                  :: numtask,taskid,stat, &
                                              h_grid_dim, grid_dim
  INTEGER                                  :: ipe,jpe, index_order, &
                                              rest_col,ncol, nrow, &
                                              end_col,start_row,end_row, &
                                              done, rest_row,kpe,has_2dgrid, &
                                              processor
  REAL                                     :: contention_X,contention_Y, contention_Z,&
                                              contention_hp, contention_packed  


 ! Get a new MPI rank ordering, reduce comm costs
  CALL mp_environ(numtask,taskid,mp_comm)

#ifdef __DEBUG_HWTOPO
   IF(taskid==0) WRITE(*,*) 'Starting Cannon heuristic'
#endif
  IF(square_grid(numtask) .OR. (IAND(numtask,numtask-1) .EQ. 0)) THEN
   
     !create the cannon graph
     CALL ma_comm_graph (mp_comm,cannon_graph)
   
     ALLOCATE(ranks_order(numtask),stat=stat)
     IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
     ranks_order(:) = 0
   
     CALL ma_2dgrid_dimensions(ncol,nrow,numtask)
     has_2dgrid = has_more_2dgrid(ncol,nrow,numtask)
   
     CALL ma_compute_contention(contention,mp_comm)
     contention_packed = SUM (contention)
     contention_packed = contention_packed+contention_packed/10
    
      h_grid_dim = SQRT(REAL(ma_find_max_square_grid(numtask),KIND=dp))
   
      IF ( numtask .EQ. (h_grid_dim*h_grid_dim)) THEN
       CALL ma_hilbert_order(ranks_order,h_grid_dim)
   
#ifdef __DEBUG_HWTOPO
    IF(taskid==0) WRITE(*,*) 'Computing Hilbert curve - Only one curve is required'
#endif
       
      ! there are other sub grids
      ELSE IF ( has_2dgrid .GT. 0 ) THEN
       CALL ma_hilbert_recursive(ranks_order,ncol=ncol,nrow=nrow,numtask=numtask,&
                                 current_numtask=numtask, current_nrow=nrow, &
                                 current_ncol=ncol, start_grid=0, &
                                 working_col=.TRUE.,mp_comm=mp_comm)
   
       CALL ma_hilbert_recursive(ranks_order,ncol=ncol,nrow=nrow,numtask=numtask,&
                                 current_numtask=numtask, current_nrow=nrow, &
                                 current_ncol=ncol, start_grid=0, &
                                 working_col=.FALSE.,mp_comm=mp_comm)
   
      DO ipe=2,numtask
       IF( ranks_order(ipe) == 0 ) ranks_order(ipe) = ipe-1
      ENDDO
   
#ifdef __DEBUG_HWTOPO
     IF(taskid==0) WRITE(*,*) 'Computing Hilbert curve - More than one curve is required'
#endif
   
     !no perfect hilbert grids, use peano
      ELSE
   
#ifdef __DEBUG_HWTOPO
      IF(taskid==0) WRITE(*,*) 'Computing Hilbert-peano curve'
#endif
   
       ALLOCATE(tmp_order(h_grid_dim*h_grid_dim),stat=stat)
       IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
   
       ALLOCATE(coord_order(h_grid_dim*h_grid_dim),stat=stat)
       IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
   
       coord_order(:) = 0
       tmp_order(:) = 0
   
       CALL ma_hilbert_order(coord_order,h_grid_dim)
   
      done = 1
      DO ipe=1,h_grid_dim
       DO jpe=0,h_grid_dim-1
         tmp_order(jpe*h_grid_dim+ipe) = coord_order(done)
         done = done + 1
       ENDDO
      ENDDO
   
       coord_order(:) = tmp_order(:)
   
       DEALLOCATE(tmp_order,stat=stat)
       IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
       ALLOCATE(tmp_order(numtask),stat=stat)
       IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
       tmp_order(:) = 0
   
       !Now include the column ranks that were not in the Hilbert grid
       CALL ma_2dgrid_dimensions(ncol,nrow,numtask)
       rest_col = ncol-h_grid_dim
       rest_row = nrow-h_grid_dim
       index_order = 1
       grid_dim = SQRT(REAL(numtask,KIND=dp))
   
       DO ipe=0,nrow-1
        DO jpe=ipe*ncol+h_grid_dim, ipe*ncol+h_grid_dim+rest_col-1
         tmp_order(index_order) = jpe
         index_order = index_order + 1
        ENDDO
       ENDDO
   
       end_col = index_order
       start_row = index_order
   
       !Now include de row ranks that were not in the Hilbert grid
       DO ipe=0,rest_row-1
   !      DO jpe=numtask-rest_col-(ipe*nrow)-1,ncol*h_grid_dim, -1
         DO jpe=numtask-rest_col-(ipe*nrow)-1,(numtask-rest_col-(ipe*nrow))-h_grid_dim, -1
           tmp_order(index_order) = jpe
           index_order = index_order + 1
         ENDDO
       ENDDO
   
       end_row = index_order
   
       !And now copy all them to the rank_order array
       index_order = 1
       jpe = 0
       kpe = 1
       DO ipe=1,numtask
        IF ( (ipe .GT. jpe*ncol) .AND. (ipe .LT. jpe*ncol+h_grid_dim+1) .AND. &
              (index_order .LE. SIZE(coord_order))) THEN
            ranks_order(ipe) = coord_order(index_order)
            index_order = index_order + 1
        ENDIF
        IF (ipe == h_grid_dim+jpe*ncol) jpe = jpe + 1
       ENDDO
   
       jpe = 0
       kpe = 1
       index_order = index_order - 1
       DO ipe=1,numtask+1
        IF ( (ipe .GE. jpe*ncol+nrow) .AND. (ipe .LT. jpe*ncol+nrow+rest_col)) THEN
            ranks_order(tmp_order(kpe)+1) = index_order
            kpe = kpe + 1
            index_order = index_order + 1
        ENDIF
        IF (ipe == (jpe+1)*ncol+1) jpe = jpe + 1
       ENDDO
   
       IF ( rest_row == 1) THEN
          DO ipe=kpe,kpe+h_grid_dim-1
           ranks_order(tmp_order(ipe)+1) = index_order
           index_order = index_order + 1
          ENDDO
       ELSE
   
         DO WHILE ( index_order .LT. numtask)
           IF (kpe .LT. SIZE(tmp_order) .AND. tmp_order(kpe) .NE. 0) THEN
             ranks_order(tmp_order(kpe)+1) = index_order
             index_order = index_order + 1
             tmp_order(kpe) = 0
           ENDIF
           kpe = kpe + h_grid_dim
           IF (kpe .LT. SIZE(tmp_order) .AND. tmp_order(kpe) .NE. 0) THEN
             ranks_order(tmp_order(kpe)+1) = index_order
             index_order = index_order + 1
             tmp_order(kpe) = 0
           ENDIF
           kpe = kpe + 1
           IF (kpe .LT. SIZE(tmp_order) .AND. tmp_order(kpe) .NE. 0) THEN
             ranks_order(tmp_order(kpe)+1) = index_order
             index_order = index_order + 1
             tmp_order(kpe) = 0
           ENDIF
           kpe = kpe - h_grid_dim
           IF (kpe .LT. SIZE(tmp_order) .AND. tmp_order(kpe) .NE. 0) THEN
             ranks_order(tmp_order(kpe)+1) = index_order
             index_order = index_order + 1
             tmp_order(kpe) = 0
           ENDIF
           kpe = kpe + 1
         ENDDO
       ENDIF
   
   ENDIF
   
     ALLOCATE(tmp_order(numtask),stat=stat)
     IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
     tmp_order(:) = 0
    
     tmp_order(:) = allocated_topology(:)
   
     DO ipe=1, numtask
       allocated_topology(ranks_order(ipe)+1) = tmp_order(ipe)
     ENDDO
     
     CALL ma_compute_contention(contention,mp_comm)
     contention_hp = SUM (contention)
   
     allocated_topology(:) = tmp_order(:)
   
#if defined __GEMINI || __SEASTAR
#ifdef __DEBUG_HWTOPO
      IF(taskid==0) WRITE(*,*) 'Contention with Hilbert curve ', contention_hp, ' Contention with packed ', contention_packed
#endif  
   
     IF( contention_hp .LT. contention_packed) THEN
   
#ifdef __DEBUG_HWTOPO
      IF(taskid==0) WRITE(*,*) 'Contention is smaller with Hilbert curve'
#endif
   
     ALLOCATE(orderX(numtask),stat=stat)
     IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
     orderX(:) = 0
     
     ALLOCATE(orderY(numtask),stat=stat)
     IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
     orderY(:) = 0
   
     ALLOCATE(orderZ(numtask),stat=stat)
     IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
     orderZ(:) = 0
   
    !Project hilbert-peano in each direction and check for contention
     ipe = 1
     kpe = 0
     DO WHILE (ipe .LE. numtask)
       DO jpe=1, numtask
         processor = allocated_topology(jpe)
         IF (net_topology(processor+1)%coord_x .EQ. kpe) THEN
           orderX(jpe) = ranks_order(ipe)
           ipe = ipe + 1
         ENDIF
       ENDDO
      kpe = kpe + 1
     ENDDO 
   
     tmp_order(:) = allocated_topology(:)
   
     DO ipe=1, numtask
       allocated_topology(orderX(ipe)+1) = tmp_order(ipe)
     ENDDO
   
     CALL ma_compute_contention(contention,mp_comm)
     contention_X = SUM (contention)
   
      allocated_topology(:) = tmp_order(:)
   
     ipe = 1
     kpe = 0
     DO WHILE (ipe .LE. numtask)
       DO jpe=1, numtask
         processor = allocated_topology(jpe)
         IF (net_topology(processor+1)%coord_y .EQ. kpe) THEN
           orderY(jpe) = ranks_order(ipe)
           ipe = ipe + 1
         ENDIF
       ENDDO
       kpe = kpe + 1
     ENDDO
   
     tmp_order(:) = allocated_topology(:)
   
     DO ipe=1, numtask
       allocated_topology(orderY(ipe)+1) = tmp_order(ipe)
     ENDDO
   
     CALL ma_compute_contention(contention,mp_comm)
     contention_Y = SUM (contention)
   
     allocated_topology(:) = tmp_order(:)
   
     ipe = 1
     kpe = 0
     DO WHILE (ipe .LE. numtask)
       DO jpe=1, numtask
         processor = allocated_topology(jpe)
         IF (net_topology(processor+1)%coord_z .EQ. kpe) THEN
           orderZ(jpe) = ranks_order(ipe)
           ipe = ipe + 1
         ENDIF
       ENDDO
       kpe = kpe + 1
     ENDDO
   
     tmp_order(:) = allocated_topology(:)
    
     DO ipe=1, numtask
       allocated_topology(orderZ(ipe)+1) = tmp_order(ipe)
     ENDDO
   
     CALL ma_compute_contention(contention,mp_comm)
     contention_Z = SUM (contention)
   
     allocated_topology(:) = tmp_order(:)
   
#ifdef __DEBUG_HWTOPO
      IF(taskid==0) WRITE(*,*) 'Contention: ', contention_X,contention_Y,contention_Z
#endif
    
     IF ( (contention_X .LE. contention_Y+0.4*contention_Y) .AND. &
          (contention_X .LE. contention_Z+0.2*contention_Z)) THEN
        ranks_order(:) = orderX(:)
        tmp_order(:) = allocated_topology(:)
   
        DO ipe=1, numtask
           allocated_topology(orderX(ipe)+1) = tmp_order(ipe)
        ENDDO
   
#ifdef __DEBUG_HWTOPO
      IF(taskid==0) WRITE(*,*) 'Use hilbert curve in X direction'
#endif
      
     ELSE IF ( (contention_Y .LE. contention_X) .AND. &
          (contention_Y .LE. contention_Z)) THEN
         ranks_order(:) = orderY(:)
         tmp_order(:) = allocated_topology(:)
         DO ipe=1, numtask
           allocated_topology(orderY(ipe)+1) = tmp_order(ipe)
        ENDDO
#ifdef __DEBUG_HWTOPO
     IF(taskid==0) WRITE(*,*) 'Use hilbert curve in Y direction'
#endif
     ELSE 
        ranks_order(:) = orderZ(:)
        tmp_order(:) = allocated_topology(:)
        DO ipe=1, numtask
           allocated_topology(orderZ(ipe)+1) = tmp_order(ipe)
        ENDDO
#ifdef __DEBUG_HWTOPO
      IF(taskid==0) WRITE(*,*) 'Use hilbert curve in Z direction'
#endif
    ENDIF
   ENDIF
#endif
   
     !hilbert does not improve the overall bw
     IF( contention_hp .GT. contention_packed) THEN
#ifdef __DEBUG_HWTOPO
      IF(taskid==0) WRITE(*,*) 'Back to packed! Hilbert is not improving contention'
#endif
       DO ipe=1, numtask
           ranks_order(ipe) = ipe - 1
       ENDDO
#if defined __GEMINI || __SEASTAR
      IF(ma_get_job_allocation() .NE. -1) &
        CALL ma_switch(mp_new_comm,mp_comm,ranks_order,complete_graph)
#else
       CALL ma_switch(mp_new_comm,mp_comm,ranks_order,complete_graph)
#endif
     ELSE
#if defined __GEMINI || __SEASTAR
      IF(ma_get_job_allocation() .NE. -1)&
       CALL ma_switch(mp_new_comm,mp_comm,ranks_order,cannon_graph)
#else
       CALL ma_switch(mp_new_comm,mp_comm,ranks_order,cannon_graph)
#endif
     ENDIF
   
#ifdef __DEBUG_HWTOPO
      IF(taskid==0) WRITE(*,*) 'New communicator generated'
#endif

  IF(ASSOCIATED(coord_order))DEALLOCATE(coord_order,stat=stat)
  IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

 IF(ASSOCIATED(tmp_order))DEALLOCATE(tmp_order,stat=stat)
  IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

  IF(ASSOCIATED(orderX))DEALLOCATE(orderX,stat=stat)
  IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

  IF(ASSOCIATED(orderY))DEALLOCATE(orderY,stat=stat)
  IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

  IF(ASSOCIATED(orderZ))DEALLOCATE(orderZ,stat=stat)
  IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

  IF(ASSOCIATED(ranks_order))DEALLOCATE(ranks_order,stat=stat)
  IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
 ELSE
   mp_new_comm = mp_comm
 ENDIF ! end power of two and square

#endif
#else
 mp_new_comm = mp_comm 
#endif

END SUBROUTINE ma_cannon

! *****************************************************************************
!> \brief Create a new communicator with the ordering designed by the user
!> the file should be named cp2k_ordering
!> \param mp_comm [input] : handle of the default communicator
!> \param mp_new_comm [output] : handle of the new  communicator
!> \par History
!>      7.2012 created [ Christiane Pousa ]
! *****************************************************************************
SUBROUTINE ma_designed (mp_new_comm,mp_comm)
  INTEGER, INTENT(IN)                           :: mp_comm
  INTEGER, INTENT(OUT)                          :: mp_new_comm

#if defined(__parallel)
  INTEGER, DIMENSION(:),POINTER                 :: ordering
  INTEGER                                       :: numtask, taskid, &
                                                   stat, rst_unit, &
                                                   ipe

  CALL mp_environ(numtask,taskid,mp_comm)

  ALLOCATE(ordering(numtask),stat=stat)
  IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)
  ordering(:) = 0
  
  IF (taskid==0) THEN
    rst_unit = 121316
    OPEN(rst_unit,FILE="cp2k_ordering",ACTION="READ",STATUS="OLD",ACCESS="SEQUENTIAL")
    DO ipe=1, numtask  
      READ(rst_unit,*) ordering(ipe)
    ENDDO
    CLOSE(rst_unit)
  ENDIF

  CALL mp_sum(ordering(:),mp_comm)

  ! Create the new communicator
  CALL mp_reordering(mp_comm,mp_new_comm,ordering)   

  IF(ASSOCIATED(ordering))DEALLOCATE(ordering,stat=stat)
  IF ( stat /= 0 ) CALL ma_error_stop(ma_error_allocation)

#else
  mp_new_comm = mp_comm
#endif  
END SUBROUTINE ma_designed

END MODULE ma_process_mapping
