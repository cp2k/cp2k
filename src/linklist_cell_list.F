!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/linklist_cell_list [1.0] *
!!
!!   NAME
!!     linklist_cell_list
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE linklist_cell_list

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE linklist_cell_types,             ONLY: cell_atoms_linklist_type,&
                                             cell_neighborlist_type,&
                                             neighbor_cells_type
  USE molecule_types,                  ONLY: particle_node_type
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell_param
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: cell_list

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

! to construct cells within the box so that it can  be
! used in the construction of all neighbor lists

SUBROUTINE cell_list ( natom_types, n_cell, rcut, tag, pnode, part, &
     box, cell_ll, start_list, neighbor, globenv, list_type )


    INTEGER, INTENT(IN)                      :: natom_types
    INTEGER, DIMENSION(:), INTENT(IN)        :: n_cell
    REAL(dbl), DIMENSION(:, :), INTENT(IN)   :: rcut
    CHARACTER(LEN=*), INTENT(IN)             :: tag
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(INOUT)                          :: pnode
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(cell_neighborlist_type), &
      DIMENSION(:, :, :), POINTER            :: cell_ll
    TYPE(cell_atoms_linklist_type), &
      DIMENSION(:, :, :, :), POINTER         :: start_list
    TYPE(neighbor_cells_type), &
      DIMENSION(:, :), POINTER               :: neighbor
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    INTEGER, INTENT(IN)                      :: list_type

    INTEGER                                  :: handle, i, ii, ineighbor, &
                                                isos, itype, iw, jj, jtype, &
                                                kk, natoms, ncell_tot, &
                                                ncells, nnodes
    INTEGER, DIMENSION(3)                    :: bin, cell_i, irange, next_cell
    INTEGER, POINTER                         :: index_ij( :, : ), &
                                                index_ji( :, : )
    LOGICAL                                  :: first_time
    REAL(dbl)                                :: delta, kint, pl
    REAL(dbl), DIMENSION(3)                  :: cell_ang, cell_lengths, &
                                                list_cell, perd, r, s, &
                                                s_list_cell
    REAL(dbl), DIMENSION(3, 3)               :: h_inv
    TYPE(cell_atoms_linklist_type), POINTER  :: loc_ll

!------------------------------------------------------------------------------

  CALL timeset( 'CELL_INDEX', 'I', 'KIs', handle )
!
! print banner if this is the first call to the routine
!

  iw = globenv % scr

  first_time = .NOT. ASSOCIATED(cell_ll)

  IF (globenv % ionode .AND. first_time .AND. globenv % print_level>0) THEN
     WRITE ( iw, '( A,A36,A )' ) ' **********************', tag, &
          '*********************'
  END IF

! this avoids problems for atoms that are exactly on a cell boundary
  pl = 0.5_dbl - EPSILON ( 0.0_dbl )

  CALL get_cell_param ( box, cell_lengths, cell_ang )

  list_cell ( : ) = cell_lengths ( : )/real(n_cell ( : ))

  h_inv = box%h_inv
  perd = box%perd

  IF (first_time) THEN
     IF (globenv % ionode .AND. globenv % print_level>0 ) THEN
        ncell_tot = n_cell(1)*n_cell(2)*n_cell(3)
        WRITE ( iw, '( A,T71,I10 )' ) &
             ' Total # of cells used for Verlet neighbor list ', ncell_tot
        WRITE ( iw, '( A,T60,i5,a3,i5,a3,i5 )' ) &
             ' # of cells in each direction = ', n_cell(1), ' X ', n_cell(2), &
             ' X ', n_cell(3)
     END IF
  END IF

! find out allocation parameters for neighboring cell calculation
  IF (first_time .AND. ALL(n_cell>1)) THEN
     ALLOCATE (neighbor(natom_types,natom_types),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'cell_list', &
          'neighbor_cell', natom_types ** 2 )
     DO itype = 1, natom_types
        DO jtype = 1, natom_types
           irange ( : ) = INT(rcut(itype,jtype )/list_cell ( : )) + 1
           CALL get_neighbor_cells_est(irange,ncells)
           IF (ncells>=ncell_tot) CALL stop_program( 'cell_list', &
                'NUMBER OF NEIGHBOR CELLS > TOTAL CELLS. REINITIALIZE')
           ALLOCATE (neighbor(itype,jtype ) % index(3,ncells),STAT=isos)
           IF ( isos /= 0 ) CALL stop_memory ( 'cell_list', &
                'neighbor%index', 3 * ncells )
           CALL get_neighbor_cells ( irange, neighbor ( itype, jtype ) )
        END DO
     END DO
  END IF

  nnodes = SIZE(pnode)
  natoms = SIZE(part)

  IF ( .NOT. ASSOCIATED(cell_ll)) THEN
     ALLOCATE (cell_ll(n_cell(1),n_cell(2),n_cell(3)),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'cell_list', 'cell_ll', 0 )
     DO kk = 1, n_cell(3)
        DO jj = 1, n_cell(2)
           DO ii = 1, n_cell(1)
              ALLOCATE (cell_ll(ii,jj,kk) % natoms(natom_types),STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'cell_list', &
                   'cell_ll%natoms', natom_types )
           END DO
        END DO
     END DO
  END IF

  s_list_cell = list_cell / cell_lengths

! linked list that forms a list of all atoms in a particular cell
  DO kk = 1, n_cell(3)
     DO jj = 1, n_cell(2)
        DO ii = 1, n_cell(1)
           cell_ll(ii,jj,kk) % natoms ( : ) = 0
        END DO
     END DO
  END DO

  IF (first_time) THEN
     DO kk = 1, n_cell(3)
        DO jj = 1, n_cell(2)
           DO ii = 1, n_cell(1)
              ALLOCATE (cell_ll(ii,jj,kk) % ptype ( natom_types),STAT=isos)
              IF ( isos /= 0 ) CALL stop_memory ( 'cell_list', &
                   'cell_ll%ptype', natom_types )
              DO itype = 1, natom_types
                 ALLOCATE (start_list(itype,ii,jj,kk),STAT=isos)
                 IF ( isos /= 0 ) CALL stop_memory ( 'cell_list', &
                      'start_list', 0 )
                 NULLIFY (start_list(itype,ii,jj,kk) % next)
              END DO
           END DO
        END DO
     END DO

#if defined ( __PGI ) || defined ( __AIX ) || defined ( __DEC )
!*apsi* Due to some strange "feature" in PGI-3.2
     DO kk = 1, n_cell(3)
        DO jj = 1, n_cell(2)
           DO ii = 1, n_cell(1)
              DO itype = 1, natom_types
                 NULLIFY (start_list(itype,ii,jj,kk) % next)
              END DO
           END DO
        END DO
     END DO
#endif
  END IF

  DO kk = 1, n_cell(3)
     DO jj = 1, n_cell(2)
        DO ii = 1, n_cell(1)
           DO itype = 1, natom_types
              cell_ll(ii,jj,kk) % ptype ( itype ) % ll     &
                   => start_list(itype,ii,jj,kk)

              loc_ll =>  cell_ll(ii,jj,kk) % ptype ( itype ) % ll
              loc_ll % used = .FALSE.

              DO
                 IF ( ASSOCIATED ( loc_ll % next)) THEN
                    loc_ll => loc_ll % next
                    loc_ll % used = .FALSE.
                 ELSE
                    EXIT
                 END IF
              END DO

           END DO
        END DO
     END DO
  END DO

  DO i = 1, natoms
     SELECT CASE (list_type)
     CASE (1)
        itype = part(i) % prop % ptype
     CASE (2)
        itype = part(i) % kind            
        IF (itype <=0 ) CYCLE
     END SELECT
     IF (itype > natom_types) CYCLE
! pbc so that particles outside the central box are also mapped
!    onto proper cells
     r(1) = part(i) % r(1)
     r(2) = part(i) % r(2)
     r(3) = part(i) % r(3)
     s(1) = h_inv(1,1)*r(1) + h_inv(1,2)*r(2) + h_inv(1,3)*r(3)
     s(2) = h_inv(2,1)*r(1) + h_inv(2,2)*r(2) + h_inv(2,3)*r(3)
     s(3) = h_inv(3,1)*r(1) + h_inv(3,2)*r(2) + h_inv(3,3)*r(3)
     IF(ABS(s(1)) > 0.5_dbl) THEN
        s(1) = s(1) - perd(1)*INT(s(1)+SIGN(0.5_dbl,s(1)))
     END IF
     IF(ABS(s(2)) > 0.5_dbl) THEN
        s(2) = s(2) - perd(2)*INT(s(2)+SIGN(0.5_dbl,s(2)))
     END IF
     IF(ABS(s(3)) > 0.5_dbl) THEN
        s(3) = s(3) - perd(3)*INT(s(3)+SIGN(0.5_dbl,s(3)))
     END IF

! compute the bin number
     bin = int((s+pl)/s_list_cell) + 1

! increment the counter
     cell_ll(bin(1),bin(2),bin(3)) % natoms(itype ) = cell_ll(bin(1), &
          bin(2),bin(3)) %natoms(itype ) + 1

! load in the atom number
     loc_ll => cell_ll(bin(1),bin(2),bin(3)) % pTYPE ( itype ) % ll
     DO
        IF ( .NOT. loc_ll % used ) THEN
           loc_ll % used = .TRUE.
           EXIT
        ELSE IF ( ASSOCIATED ( loc_ll % next)) THEN
           loc_ll => loc_ll % next
        ELSE
           ALLOCATE (loc_ll % next, STAT=isos)
           IF ( isos /= 0 ) CALL stop_memory ( 'cell_list', 'll%next', 0 )
           loc_ll => loc_ll % next
           NULLIFY (loc_ll % next)
           loc_ll % used = .FALSE.
        END IF
     END DO
     loc_ll % atom = i
     SELECT CASE (list_type)
     CASE (1)
       part(i) % cell => cell_ll(bin(1),bin(2),bin(3))
     CASE (2)
       part(i) % cell_ol => cell_ll(bin(1),bin(2),bin(3))
     END SELECT
  END DO

! linked list that forms a list of all neighboring cells (done only once)
  IF ( first_time .AND. ALL ( n_cell > 1 ) ) THEN
     delta = 1.0E-5_dbl

! Note:
! floor returns the largest integer < or = to its REAL argument.
! n_cell is the periodicity. Thus, floor(n_cell/n_cell) gives back 1,
! although we want it to be zero. So, we subtract the argument by a
! small number (delta) to give us zero.
     DO kk = 1, n_cell(3)
        DO jj = 1, n_cell(2)
           DO ii = 1, n_cell(1)
              ALLOCATE ( cell_ll ( ii, jj, kk ) &
                   % neighbor ( natom_types, natom_types), STAT = isos )
              IF ( isos /= 0 ) CALL stop_memory ( 'cell_index', 'neighbor', 0 )
              cell_i(1) = ii
              cell_i(2) = jj
              cell_i(3) = kk

! loop over number of types
              DO itype = 1, natom_types
                 DO jtype = itype, natom_types

! Define number of cells to loop over
                    ncells=SIZE(neighbor(itype,jtype ) % index, 2)

                    ALLOCATE( cell_ll (ii,jj,kk) % neighbor(itype,jtype ) &
                         % index ( 3, ncells ) )
                    IF ( isos /= 0 ) CALL stop_memory ( 'cell_index', 'index', 0 )
                    IF (itype/=jtype ) THEN
                       ALLOCATE( cell_ll (ii,jj,kk) % neighbor(jtype,itype ) &
                            % index (3, ncells ) )
                       IF ( isos /= 0 ) CALL stop_memory ( 'cell_index', 'index', 0 )
                    END IF

! POINT index_ij => cell_ll % neighbor % index
                    index_ij => cell_ll(ii,jj,kk) % neighbor(itype,jtype ) % index
                    IF (itype/=jtype ) THEN
                       index_ji => cell_ll(ii,jj,kk) % neighbor(jtype,itype ) % index
                    END IF

! loop over the number of cells to find the cell neighbors of
!    cell_ll(ii,jj,kk)
                    DO ineighbor = 1, ncells
                       next_cell(1) = neighbor(itype,jtype ) % index(1, &
                            ineighbor) + cell_i(1)
                       next_cell(2) = neighbor(itype,jtype ) % index(2, &
                            ineighbor) + cell_i(2)
                       next_cell(3) = neighbor(itype,jtype ) % index(3, &
                            ineighbor) + cell_i(3)

! pbc the cell number
                       next_cell(1) = next_cell(1) - n_cell(1)* &
                            floor((real(next_cell(1)))/ real(n_cell(1))-delta)
                       next_cell(2) = next_cell(2) - n_cell(2)* &
                            floor((real(next_cell(2)))/ real(n_cell(2))-delta)
                       next_cell(3) = next_cell(3) - n_cell(3)* &
                            floor((real(next_cell(3)))/ real(n_cell(3))-delta)

! load in the index
                       index_ij(1,ineighbor) = next_cell(1)
                       index_ij(2,ineighbor) = next_cell(2)
                       index_ij(3,ineighbor) = next_cell(3)
                       IF (itype/=jtype ) THEN
                          index_ji(1,ineighbor) = next_cell(1)
                          index_ji(2,ineighbor) = next_cell(2)
                          index_ji(3,ineighbor) = next_cell(3)
                       END IF

                    END DO
                 END DO

              END DO
           END DO
        END DO
     END DO
  END IF

  IF (globenv % ionode .AND. first_time .AND. globenv % print_level>0) THEN
     WRITE ( iw, '( A,A,A,/ )' ) ' **************************', &
          ' end subroutine cell_list ', '***************************'
  END IF

  kint = PRODUCT ( n_cell ( 1:3 ) ) * 0.001_dbl
  CALL timestop ( kint, handle )

END SUBROUTINE cell_list

!******************************************************************************

! this routine gives back the indices of cell neighbors of a cell (0,0,0).
! only the neighbors on the right side are returned

SUBROUTINE get_neighbor_cells ( n_max, neighbor )


    INTEGER, DIMENSION(:), INTENT(IN)        :: n_max
    TYPE(neighbor_cells_type), INTENT(OUT)   :: neighbor

    INTEGER                                  :: i, j, k, n

!------------------------------------------------------------------------------
! along nz=0, ny=0, nx=1, nx_max

  n = 0
  DO i = 1, n_max(1)
     n = n + 1
     neighbor % index(1,n) = i
     neighbor % index(2,n) = 0
     neighbor % index(3,n) = 0
  END DO

! along nz=0, ny=1, ny_max, nx=-nx_max, nx_max
  DO j = 1, n_max(2)
     DO i = -n_max(1), n_max(1)
        n = n + 1
        neighbor % index(1,n) = i
        neighbor % index(2,n) = j
        neighbor % index(3,n) = 0
     END DO
  END DO

! along nz=1, nz_max, ny=-ny_max, ny_max, nx=-nx_max, nx_max
  DO k = 1, n_max(3)
     DO j = -n_max(2), n_max(2)
        DO i = -n_max(1), n_max(1)
           n = n + 1
           neighbor % index(1,n) = i
           neighbor % index(2,n) = j
           neighbor % index(3,n) = k
        END DO
     END DO
  END DO

END SUBROUTINE get_neighbor_cells

!******************************************************************************

! Returns the number of neighbor cells so we can ALLOCATE neighbor_cells

SUBROUTINE get_neighbor_cells_est ( n_max, n )


    INTEGER, DIMENSION(:), INTENT(IN)        :: n_max
    INTEGER, INTENT(OUT)                     :: n

    INTEGER                                  :: i, j, k

!------------------------------------------------------------------------------

  n = 0
  DO i = 1, n_max(1)
     n = n + 1
  END DO

! along nz=0, ny=1, ny_max, nx=-nx_max, nx_max
  DO j = 1, n_max(2)
     DO i = -n_max(1), n_max(1)
        n = n + 1
     END DO
  END DO

! along nz=1, nz_max, ny=-ny_max, ny_max, nx=-nx_max, nx_max
  DO k = 1, n_max(3)
     DO j = -n_max(2), n_max(2)
        DO i = -n_max(1), n_max(1)
           n = n + 1
        END DO
     END DO
  END DO

END SUBROUTINE get_neighbor_cells_est

!******************************************************************************

END MODULE linklist_cell_list


