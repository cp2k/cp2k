!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Define the neighbor list data types and the corresponding functionality
!> \par History
!>      - cleaned (23.07.2003,MK)
!>      - full refactoring, list iterators (20.10.2010, JGH)
!>      - add get_neighbor_list_set_p, return info for a set of neighborlists
!>                                                             (07.2014,JGH)
!> \author Matthias Krack (21.06.2000)
! *****************************************************************************
MODULE qs_neighbor_list_types

  USE kinds,                           ONLY: dp,&
                                             int_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE util,                            ONLY: locate,&
                                             sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_neighbor_list_types'

! *** Definition of the data types for a linked list of neighbors ***

! *****************************************************************************
  TYPE neighbor_node_type
    PRIVATE
    TYPE(neighbor_node_type), POINTER :: next_neighbor_node
    REAL(dp), DIMENSION(3)            :: r
    INTEGER, DIMENSION(3)             :: cell
    INTEGER                           :: neighbor
  END TYPE neighbor_node_type

! *****************************************************************************
  TYPE neighbor_list_type
    PRIVATE
    TYPE(neighbor_list_type), POINTER :: next_neighbor_list
    TYPE(neighbor_node_type), POINTER :: first_neighbor_node,&
                                         last_neighbor_node
    INTEGER                           :: atom,nnode
  END TYPE neighbor_list_type

! *****************************************************************************
  TYPE neighbor_list_set_type
    PRIVATE
    TYPE(neighbor_list_type), POINTER :: first_neighbor_list,&
                                         last_neighbor_list
    REAL(dp)                          :: r_max
    INTEGER                           :: nlist
    LOGICAL                           :: symmetric
    LOGICAL                           :: mic
    LOGICAL                           :: molecular
  END TYPE neighbor_list_set_type

! *****************************************************************************
  TYPE neighbor_list_p_type
    TYPE(neighbor_list_type), POINTER :: neighbor_list
  END TYPE neighbor_list_p_type

! *****************************************************************************
  TYPE neighbor_list_set_p_type
    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
  END TYPE neighbor_list_set_p_type

! *****************************************************************************
  TYPE list_search_type
    PRIVATE
    INTEGER                               :: nlist
    INTEGER, DIMENSION(:), POINTER        :: atom_list
    INTEGER, DIMENSION(:), POINTER        :: atom_index
    TYPE(neighbor_list_p_type), &
      DIMENSION(:), POINTER               :: neighbor_list
  END TYPE list_search_type

! *****************************************************************************
! Neighbor List Iterator
! *****************************************************************************
  TYPE neighbor_list_iterator_type
    PRIVATE
    INTEGER                               :: ikind, jkind, ilist, inode
    INTEGER                               :: nkind, nlist, nnode
    INTEGER                               :: iatom, jatom
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER               :: nl
    TYPE(neighbor_list_type), POINTER     :: neighbor_list
    TYPE(neighbor_node_type), POINTER     :: neighbor_node
    TYPE(list_search_type), &
      DIMENSION(:), POINTER               :: list_search
  END TYPE neighbor_list_iterator_type

  TYPE neighbor_list_iterator_p_type
    PRIVATE
    TYPE(neighbor_list_iterator_type), POINTER :: neighbor_list_iterator
    INTEGER                                    :: last
  END TYPE neighbor_list_iterator_p_type
! *****************************************************************************

! *** Public data types ***

  PUBLIC :: neighbor_list_p_type,&
            neighbor_list_set_type,&
            neighbor_list_set_p_type

! *** Public subroutines ***

  PUBLIC :: add_neighbor_list,&
            add_neighbor_node,&
            allocate_neighbor_list_set,&
            deallocate_neighbor_list_set,&
            get_neighbor_list_set,&
            get_neighbor_list_set_p

! *** Iterator functions and types ***

  PUBLIC :: neighbor_list_iterator_p_type,&
            neighbor_list_iterator_create,&
            neighbor_list_iterator_release,&
            neighbor_list_iterate,&
            nl_set_sub_iterator,&
            nl_sub_iterate,&
            get_iterator_info

CONTAINS

! *****************************************************************************
!> \brief   Neighbor list iterator functions
!> \param iterator_set ...
!> \param nl ...
!> \param search ...
!> \param nthread ...
!> \date    28.07.2010
!> \author  jhu
!> \version 1.0
! *****************************************************************************
  SUBROUTINE neighbor_list_iterator_create(iterator_set,nl,search,nthread)
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: iterator_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: nl
    LOGICAL, INTENT(IN), OPTIONAL            :: search
    INTEGER, INTENT(IN), OPTIONAL            :: nthread

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'neighbor_list_iterator_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, il, ilist, istat, &
                                                mthread, nlist
    TYPE(list_search_type), DIMENSION(:), &
      POINTER                                :: list_search
    TYPE(neighbor_list_iterator_type), &
      POINTER                                :: iterator
    TYPE(neighbor_list_type), POINTER        :: neighbor_list

    mthread = 1
    IF(PRESENT(nthread)) mthread=nthread

    ALLOCATE(iterator_set(0:mthread-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "iterator_set",int_size*mthread)

    DO il=0,mthread-1
       ALLOCATE (iterator_set(il)%neighbor_list_iterator,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "neighbor_list_iterator",0)

       iterator => iterator_set(il)%neighbor_list_iterator

       iterator%nl => nl

       iterator%ikind = 0
       iterator%jkind = 0
       iterator%nkind = NINT(SQRT(REAL(SIZE(nl),dp)))

       iterator%ilist = 0
       iterator%nlist = 0
       iterator%inode = 0
       iterator%nnode = 0

       iterator%iatom = 0
       iterator%jatom = 0

       NULLIFY(iterator%neighbor_list)
       NULLIFY(iterator%neighbor_node)
       NULLIFY(iterator%list_search)
    END DO

    iterator_set(:)%last = 0

    IF (PRESENT(search)) THEN
      IF (search) THEN
         ALLOCATE(list_search(SIZE(nl)),STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                          "list_search",0)
         DO il=1,SIZE(nl)
           IF (ASSOCIATED(nl(il)%neighbor_list_set)) THEN
              CALL get_neighbor_list_set(neighbor_list_set=nl(il)%neighbor_list_set,nlist=nlist)
              list_search(il)%nlist = nlist
              ALLOCATE (list_search(il)%atom_list(nlist),STAT=istat)
              IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                               "list_search(il)%atom_list",int_size*nlist)
              ALLOCATE (list_search(il)%atom_index(nlist),STAT=istat)
              IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                               "list_search(il)%atom_index",int_size*nlist)
              ALLOCATE (list_search(il)%neighbor_list(nlist),STAT=istat)
              IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                               "list_search(il)%neighbor_list",int_size*nlist)

              NULLIFY ( neighbor_list )
              DO ilist=1,nlist
                 IF ( .NOT. ASSOCIATED(neighbor_list) ) THEN
                    neighbor_list => first_list(nl(il)%neighbor_list_set)
                 ELSE
                    neighbor_list => neighbor_list%next_neighbor_list
                 END IF
                 CALL get_neighbor_list(neighbor_list=neighbor_list,atom=iatom)
                 list_search(il)%atom_list(ilist) = iatom
                 list_search(il)%neighbor_list(ilist)%neighbor_list => neighbor_list
              END DO
              CALL sort (list_search(il)%atom_list,nlist,list_search(il)%atom_index)

           ELSE
              list_search(il)%nlist = -1
              NULLIFY(list_search(il)%atom_list,list_search(il)%atom_index,list_search(il)%neighbor_list)
           END IF
         END DO
         DO il=0,mthread-1
            iterator => iterator_set(il)%neighbor_list_iterator
            iterator%list_search => list_search
         END DO
      END IF
    END IF

  END SUBROUTINE neighbor_list_iterator_create

! *****************************************************************************
!> \brief ...
!> \param iterator_set ...
! *****************************************************************************
  SUBROUTINE neighbor_list_iterator_release(iterator_set)
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: iterator_set

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'neighbor_list_iterator_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: il, istat, mthread
    TYPE(neighbor_list_iterator_type), &
      POINTER                                :: iterator

!all threads have the same search list

    iterator => iterator_set(0)%neighbor_list_iterator
    IF (ASSOCIATED(iterator%list_search)) THEN
       DO il=1,SIZE(iterator%list_search)
          IF (iterator%list_search(il)%nlist >= 0) THEN
             DEALLOCATE (iterator%list_search(il)%atom_list,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "iterator%list_search(il)%atom_list")
             DEALLOCATE (iterator%list_search(il)%atom_index,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "iterator%list_search(il)%atom_index")
             DEALLOCATE (iterator%list_search(il)%neighbor_list,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "iterator%list_search(il)%neighbor_list")
          END IF
       END DO
       DEALLOCATE (iterator%list_search,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"iterator%list_search")
    END IF

    mthread = SIZE(iterator_set)
    DO il=0,mthread-1
       DEALLOCATE (iterator_set(il)%neighbor_list_iterator,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"neighbor_list_iterator")
    END DO
    DEALLOCATE (iterator_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"iterator_set")

  END SUBROUTINE neighbor_list_iterator_release

! *****************************************************************************
!> \brief ...
!> \param iterator_set ...
!> \param ikind ...
!> \param jkind ...
!> \param iatom ...
! *****************************************************************************
  SUBROUTINE nl_set_sub_iterator(iterator_set,ikind,jkind,iatom)
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: iterator_set
    INTEGER, INTENT(IN)                      :: ikind, jkind, iatom

    CHARACTER(len=*), PARAMETER :: routineN = 'nl_set_sub_iterator', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ij, il, ilist, mthread, &
                                                nlist, nnode
    TYPE(list_search_type), POINTER          :: list_search
    TYPE(neighbor_list_iterator_type), &
      POINTER                                :: iterator
    TYPE(neighbor_list_type), POINTER        :: neighbor_list

    iterator => iterator_set(0)%neighbor_list_iterator
    ij = ikind + iterator%nkind*(jkind - 1)
    IF (ASSOCIATED(iterator%list_search)) THEN
       list_search => iterator%list_search(ij)
       nlist = list_search%nlist
       ilist = 0
       NULLIFY(neighbor_list)
       IF (nlist > 0) THEN
         i = locate(list_search%atom_list,iatom)
         i = list_search%atom_index(i)
         IF (i > 0) neighbor_list => list_search%neighbor_list(i)%neighbor_list
         ilist = i
       END IF
       IF (ASSOCIATED(neighbor_list)) THEN
          CALL get_neighbor_list(neighbor_list=neighbor_list,nnode=nnode)
       ELSE
          nnode=0
       END IF
    ELSE
       STOP
    END IF

    mthread = SIZE(iterator_set)
    DO il=0,mthread-1
       iterator => iterator_set(il)%neighbor_list_iterator

       iterator%ikind = ikind
       iterator%jkind = jkind

       iterator%ilist = ilist
       iterator%nlist = nlist
       iterator%inode = 0
       iterator%nnode = nnode

       iterator%iatom = iatom
       iterator%jatom = 0

       iterator%neighbor_list => neighbor_list
       NULLIFY(iterator%neighbor_node)
    END DO

    iterator_set(:)%last = 0

  END SUBROUTINE nl_set_sub_iterator

! *****************************************************************************
!> \brief ...
!> \param iterator_set ...
!> \param mepos ...
!> \retval istat ...
! *****************************************************************************
  FUNCTION neighbor_list_iterate(iterator_set,mepos) RESULT(istat)
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: iterator_set
    INTEGER, OPTIONAL                        :: mepos
    INTEGER                                  :: istat

    CHARACTER(len=*), PARAMETER :: routineN = 'neighbor_list_iterate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iab, last, me
    TYPE(neighbor_list_iterator_type), &
      POINTER                                :: iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: nl

    IF(PRESENT(mepos)) THEN
      me=mepos
    ELSE
      me=0
    END IF

    istat  = 0

!$omp critical(neighbour_list_iterate_critical)
    last = iterator_set(0)%last
    IF(last /= me) THEN
      iterator_set(me)%neighbor_list_iterator = iterator_set(last)%neighbor_list_iterator
    END IF
    iterator => iterator_set(me)%neighbor_list_iterator
    nl => iterator%nl

    IF(iterator%inode < iterator%nnode) THEN
      ! we can be sure that there is another node in this list
      iterator%inode = iterator%inode + 1
      iterator%neighbor_node => iterator%neighbor_node%next_neighbor_node
    ELSE
      iab = MAX(iterator%ikind + iterator%nkind*(iterator%jkind - 1),0)
      kindloop : DO ! look for the next list with nnode /= 0
         listloop : DO
            IF(iterator%ilist >= iterator%nlist) EXIT listloop
            iterator%ilist = iterator%ilist + 1
            IF(ASSOCIATED(iterator%neighbor_list)) THEN
               iterator%neighbor_list => iterator%neighbor_list%next_neighbor_list
            ELSE
               iterator%neighbor_list => first_list(nl(iab)%neighbor_list_set)
            END IF
            CALL get_neighbor_list(neighbor_list=iterator%neighbor_list,atom=iterator%iatom,&
                                 nnode=iterator%nnode)
            IF(iterator%nnode > 0) EXIT kindloop
         END DO listloop
         IF (iab >= iterator%nkind**2) THEN
           istat = 1
           EXIT kindloop
         ELSE
           iab = iab + 1
           iterator%jkind = (iab-1)/iterator%nkind + 1
           iterator%ikind = iab - iterator%nkind*(iterator%jkind - 1)
           iterator%ilist = 0
           IF (.NOT.ASSOCIATED(nl(iab)%neighbor_list_set)) THEN
             iterator%ilist = 0
             iterator%nlist = 0
           ELSE
             CALL get_neighbor_list_set(neighbor_list_set=&
                     nl(iab)%neighbor_list_set,nlist=iterator%nlist)
             iterator%ilist = 0
           END IF
           NULLIFY(iterator%neighbor_list)
         END IF
      END DO kindloop
      IF (istat==0) THEN
        iterator%inode = 1
        iterator%neighbor_node => first_node(iterator%neighbor_list)
      END IF
    END IF
    IF (istat==0) THEN
      CALL get_neighbor_node(neighbor_node=iterator%neighbor_node,neighbor=iterator%jatom)
    END IF

    ! mark the last iterator updated
    iterator_set(:)%last = me
!$omp end critical(neighbour_list_iterate_critical)

  END FUNCTION neighbor_list_iterate

! *****************************************************************************
!> \brief ...
!> \param iterator_set ...
!> \retval istat ...
! *****************************************************************************
  FUNCTION nl_sub_iterate(iterator_set) RESULT(istat)
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: iterator_set
    INTEGER                                  :: istat

    CHARACTER(len=*), PARAMETER :: routineN = 'nl_sub_iterate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last, me
    TYPE(neighbor_list_iterator_type), &
      POINTER                                :: iterator

! Sub-iterator single-threaded, or with private sub-iterators only

    me = 0
    istat = 0

    last = iterator_set(0)%last
    IF(last /= me) THEN
      iterator_set(me)%neighbor_list_iterator = iterator_set(last)%neighbor_list_iterator
    END IF
    iterator => iterator_set(me)%neighbor_list_iterator

    IF (ASSOCIATED(iterator%neighbor_list)) THEN
       IF (iterator%inode >= iterator%nnode) THEN
          ! end of loop
          istat  = 1
       ELSEIF (iterator%inode == 0) THEN
          iterator%inode = 1
          iterator%neighbor_node => first_node(iterator%neighbor_list)
       ELSEIF (iterator%inode > 0) THEN
          ! we can be sure that there is another node in this list
          iterator%inode = iterator%inode + 1
          iterator%neighbor_node => iterator%neighbor_node%next_neighbor_node
       ELSE
          STOP "wrong"
       END IF
    ELSE
       ! no list available
       istat  = 1
    END IF
    IF (istat==0) THEN
      CALL get_neighbor_node(neighbor_node=iterator%neighbor_node,neighbor=iterator%jatom)
    END IF

    ! mark the last iterator updated
    iterator_set(:)%last = me

  END FUNCTION nl_sub_iterate

! *****************************************************************************
!> \brief ...
!> \param iterator_set ...
!> \param mepos ...
!> \param ikind ...
!> \param jkind ...
!> \param nkind ...
!> \param ilist ...
!> \param nlist ...
!> \param inode ...
!> \param nnode ...
!> \param iatom ...
!> \param jatom ...
!> \param r ...
!> \param cell ...
! *****************************************************************************
  SUBROUTINE get_iterator_info(iterator_set,mepos,&
    ikind,jkind,nkind,ilist,nlist,inode,nnode,iatom,jatom,r,cell)
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: iterator_set
    INTEGER, OPTIONAL                        :: mepos, ikind, jkind, nkind, &
                                                ilist, nlist, inode, nnode, &
                                                iatom, jatom
    REAL(dp), DIMENSION(3), OPTIONAL         :: r
    INTEGER, DIMENSION(3), OPTIONAL          :: cell

    CHARACTER(len=*), PARAMETER :: routineN = 'get_iterator_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: me
    TYPE(neighbor_list_iterator_type), &
      POINTER                                :: iterator

    IF(PRESENT(mepos)) THEN
      me=mepos
    ELSE
      me=0
    END IF
    iterator => iterator_set(me)%neighbor_list_iterator

    IF(PRESENT(ikind)) ikind=iterator%ikind
    IF(PRESENT(jkind)) jkind=iterator%jkind
    IF(PRESENT(nkind)) nkind=iterator%nkind
    IF(PRESENT(ilist)) ilist=iterator%ilist
    IF(PRESENT(nlist)) nlist=iterator%nlist
    IF(PRESENT(inode)) inode=iterator%inode
    IF(PRESENT(nnode)) nnode=iterator%nnode
    IF(PRESENT(iatom)) iatom=iterator%iatom
    IF(PRESENT(jatom)) jatom=iterator%jatom
    IF(PRESENT(r)) THEN
      CALL get_neighbor_node(neighbor_node=iterator%neighbor_node,r=r)
    END IF
    IF(PRESENT(cell)) THEN
      CALL get_neighbor_node(neighbor_node=iterator%neighbor_node,cell=cell)
    END IF

  END SUBROUTINE get_iterator_info

! *****************************************************************************
!> \brief   Add a new neighbor list to a neighbor list set.
!> \param neighbor_list_set ...
!> \param atom ...
!> \param neighbor_list ...
!> \date    13.09.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE add_neighbor_list(neighbor_list_set,atom,neighbor_list)

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    INTEGER, INTENT(IN)                      :: atom
    TYPE(neighbor_list_type), POINTER        :: neighbor_list

    CHARACTER(len=*), PARAMETER :: routineN = 'add_neighbor_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    TYPE(neighbor_list_type), POINTER        :: new_neighbor_list

    IF (ASSOCIATED(neighbor_list_set)) THEN

      IF (ASSOCIATED(neighbor_list_set%last_neighbor_list)) THEN

        new_neighbor_list =>&
          neighbor_list_set%last_neighbor_list%next_neighbor_list

        IF (.NOT.ASSOCIATED(new_neighbor_list)) THEN

!         *** Allocate a new neighbor list ***

          ALLOCATE (new_neighbor_list,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routineN,moduleN,__LINE__,&
                             "new_neighbor_list",0)
          END IF

          NULLIFY (new_neighbor_list%next_neighbor_list)
          NULLIFY (new_neighbor_list%first_neighbor_node)

!         *** Link the new neighbor list to the neighbor list set ***

          neighbor_list_set%last_neighbor_list%next_neighbor_list => new_neighbor_list

        END IF

      ELSE

        new_neighbor_list => neighbor_list_set%first_neighbor_list

        IF (.NOT.ASSOCIATED(new_neighbor_list)) THEN

!         *** Allocate a new first neighbor list ***

          ALLOCATE (new_neighbor_list,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routineN,moduleN,__LINE__,&
                             "new_neighbor_list",0)
          END IF

          NULLIFY (new_neighbor_list%next_neighbor_list)
          NULLIFY (new_neighbor_list%first_neighbor_node)

!         *** Link the new first neighbor list to the neighbor list set ***

          neighbor_list_set%first_neighbor_list => new_neighbor_list

        END IF

      END IF

!     *** Store the data set of the new neighbor list ***

      NULLIFY (new_neighbor_list%last_neighbor_node)
      new_neighbor_list%atom = atom
      new_neighbor_list%nnode = 0

!     *** Update the pointer to the last neighbor ***
!     *** list of the neighbor list set           ***

      neighbor_list_set%last_neighbor_list => new_neighbor_list

!     *** Increment the neighbor list counter ***

      neighbor_list_set%nlist = neighbor_list_set%nlist + 1

!     *** Return a pointer to the new neighbor list ***

      neighbor_list => new_neighbor_list

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The requested neighbor list set is not associated")

    END IF

  END SUBROUTINE add_neighbor_list

! *****************************************************************************
!> \brief   Add a new neighbor list node to a neighbor list.
!> \param neighbor_list ...
!> \param neighbor ...
!> \param cell ...
!> \param r ...
!> \param exclusion_list ...
!> \param nkind ...
!> \date    23.06.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE add_neighbor_node(neighbor_list,neighbor,cell,r,exclusion_list,nkind)

    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    INTEGER, INTENT(IN)                      :: neighbor
    INTEGER, DIMENSION(3), INTENT(IN)        :: cell
    REAL(dp), DIMENSION(3), INTENT(IN)       :: r
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: exclusion_list
    INTEGER, INTENT(IN), OPTIONAL            :: nkind

    CHARACTER(len=*), PARAMETER :: routineN = 'add_neighbor_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, istat, my_nkind
    TYPE(neighbor_node_type), POINTER        :: new_neighbor_node

    IF (ASSOCIATED(neighbor_list)) THEN

!     *** Check for exclusions ***

      IF (PRESENT(exclusion_list)) THEN
        IF ( ASSOCIATED ( exclusion_list ) ) THEN
          DO iatom=1,SIZE(exclusion_list)
            IF (exclusion_list(iatom) == 0) EXIT
            IF (exclusion_list(iatom) == neighbor) RETURN
          END DO
        END IF
      END IF

      my_nkind = 0
      IF(PRESENT(nkind)) my_nkind = nkind

      IF (ASSOCIATED(neighbor_list%last_neighbor_node)) THEN

        new_neighbor_node => neighbor_list%last_neighbor_node%next_neighbor_node

        IF (.NOT.ASSOCIATED(new_neighbor_node)) THEN

!         *** Allocate a new neighbor node ***

          ALLOCATE (new_neighbor_node,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routineN,moduleN,__LINE__,&
                             "new_neighbor_node",0)
          END IF

          NULLIFY (new_neighbor_node%next_neighbor_node)

!         *** Link the new neighbor node to the neighbor list ***

          neighbor_list%last_neighbor_node%next_neighbor_node => new_neighbor_node

        END IF

      ELSE

        new_neighbor_node => neighbor_list%first_neighbor_node

        IF (.NOT.ASSOCIATED(new_neighbor_node)) THEN

!         *** Allocate a new first neighbor node ***

          ALLOCATE (new_neighbor_node,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routineN,moduleN,__LINE__,&
                             "new_neighbor_node",0)
          END IF

          NULLIFY (new_neighbor_node%next_neighbor_node)

!         *** Link the new first neighbor node to the neighbor list ***

          neighbor_list%first_neighbor_node => new_neighbor_node

        END IF

      END IF

!     *** Store the data set of the new neighbor ***

      new_neighbor_node%neighbor = neighbor
      new_neighbor_node%cell(:) = cell(:)
      new_neighbor_node%r(:) = r(:)

!     *** Update the pointer to the last neighbor node of the neighbor list ***

      neighbor_list%last_neighbor_node => new_neighbor_node

!     *** Increment the neighbor node counter ***

      neighbor_list%nnode = neighbor_list%nnode + 1

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The requested neighbor list is not associated")

    END IF

  END SUBROUTINE add_neighbor_node

! *****************************************************************************
!> \brief   Allocate and initialize a set of neighbor lists.
!> \param neighbor_list_set ...
!> \param r_max ...
!> \param symmetric ...
!> \param mic ...
!> \param molecular ...
!> \date    23.06.2000
!> \author MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_neighbor_list_set(neighbor_list_set,r_max,symmetric,mic,molecular)

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    REAL(dp), INTENT(IN)                     :: r_max
    LOGICAL, INTENT(IN)                      :: symmetric, mic, molecular

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_neighbor_list_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

!   *** Deallocate the old neighbor list set ***

    IF (ASSOCIATED(neighbor_list_set)) THEN
      CALL deallocate_neighbor_list_set(neighbor_list_set)
    END IF

!   *** Allocate a set of neighbor lists ***

    ALLOCATE (neighbor_list_set,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
                       "neighbor_list_set",0)
    END IF

    NULLIFY (neighbor_list_set%first_neighbor_list)

!   *** Initialize the pointers to the first neighbor list ***

    CALL init_neighbor_list_set(neighbor_list_set,r_max,symmetric,mic,molecular)

  END SUBROUTINE allocate_neighbor_list_set

! *****************************************************************************
!> \brief   Deallocate a neighbor list.
!> \param neighbor_list ...
!> \date    20.09.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_neighbor_list(neighbor_list)

    TYPE(neighbor_list_type), POINTER        :: neighbor_list

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_neighbor_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    TYPE(neighbor_node_type), POINTER        :: neighbor_node, &
                                                next_neighbor_node

    IF (ASSOCIATED(neighbor_list)) THEN

      neighbor_node => neighbor_list%first_neighbor_node

      DO WHILE (ASSOCIATED(neighbor_node))
        next_neighbor_node => neighbor_node%next_neighbor_node
        DEALLOCATE (neighbor_node,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"neighbor_node")
        END IF
        neighbor_node => next_neighbor_node
      END DO

      DEALLOCATE (neighbor_list,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"neighbor_list")
      END IF

    END IF

  END SUBROUTINE deallocate_neighbor_list

! *****************************************************************************
!> \brief   Deallocate a neighbor list set.
!> \param neighbor_list_set ...
!> \date    03.11.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_neighbor_list_set(neighbor_list_set)
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_neighbor_list_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    TYPE(neighbor_list_type), POINTER        :: neighbor_list, &
                                                next_neighbor_list

    IF (ASSOCIATED(neighbor_list_set)) THEN

      neighbor_list => neighbor_list_set%first_neighbor_list

      DO WHILE (ASSOCIATED(neighbor_list))
        next_neighbor_list => neighbor_list%next_neighbor_list
        CALL  deallocate_neighbor_list(neighbor_list)
        neighbor_list => next_neighbor_list
      END DO

      DEALLOCATE (neighbor_list_set,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"neighbor_list_set")
      END IF

    END IF

  END SUBROUTINE deallocate_neighbor_list_set

! *****************************************************************************
!> \brief   Return a pointer to the first neighbor list of a neighbor list set.
!> \param neighbor_list_set ...
!> \retval first_neighbor_list ...
!> \date    13.09.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION first_list(neighbor_list_set) RESULT(first_neighbor_list)

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: first_neighbor_list

    first_neighbor_list => neighbor_list_set%first_neighbor_list

  END FUNCTION first_list

! *****************************************************************************
!> \brief   Return a pointer to the first neighbor node of a neighbor list.
!> \param neighbor_list ...
!> \retval first_neighbor_node ...
!> \date    23.06.2000,
!> \author  MK
!> \version 1.0
! *****************************************************************************
  FUNCTION first_node(neighbor_list) RESULT(first_neighbor_node)

    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: first_neighbor_node

    first_neighbor_node => neighbor_list%first_neighbor_node

  END FUNCTION first_node

! *****************************************************************************
!> \brief   Return the reqested data of a neighbor list.
!> \param neighbor_list ...
!> \param atom ...
!> \param nnode ...
!> \date    13.09.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_neighbor_list(neighbor_list,atom,nnode)

    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    INTEGER, INTENT(OUT), OPTIONAL           :: atom, nnode

    CHARACTER(len=*), PARAMETER :: routineN = 'get_neighbor_list', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(neighbor_list)) THEN

      IF (PRESENT(atom)) atom = neighbor_list%atom
      IF (PRESENT(nnode)) nnode = neighbor_list%nnode

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The requested neighbor list is not associated")

    END IF

  END SUBROUTINE get_neighbor_list

! *****************************************************************************
!> \brief   Return the components of a neighbor list set.
!> \param neighbor_list_set ...
!> \param r_max ...
!> \param nlist ...
!> \param symmetric ...
!> \param mic ...
!> \param molecular ...
!> \date    10.11.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_neighbor_list_set(neighbor_list_set,r_max,nlist,symmetric,mic,molecular)

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    REAL(dp), INTENT(OUT), OPTIONAL          :: r_max
    INTEGER, INTENT(OUT), OPTIONAL           :: nlist
    LOGICAL, INTENT(OUT), OPTIONAL           :: symmetric, mic, molecular

    CHARACTER(len=*), PARAMETER :: routineN = 'get_neighbor_list_set', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(neighbor_list_set)) THEN

      IF (PRESENT(r_max)) r_max = neighbor_list_set%r_max
      IF (PRESENT(nlist)) nlist = neighbor_list_set%nlist
      IF (PRESENT(symmetric)) symmetric = neighbor_list_set%symmetric
      IF (PRESENT(mic)) mic = neighbor_list_set%mic
      IF (PRESENT(molecular)) molecular = neighbor_list_set%molecular

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The requested neighbor list set is not associated")

    END IF

  END SUBROUTINE get_neighbor_list_set

! *****************************************************************************
!> \brief   Return the components of the first neighbor list set.
!> \param neighbor_list_sets ...
!> \param r_max ...
!> \param nlist ...
!> \param symmetric ...
!> \param mic ...
!> \param molecular ...
!> \date    07.2014
!> \author  JGH
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_neighbor_list_set_p(neighbor_list_sets,r_max,nlist,symmetric,mic,molecular)

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: neighbor_list_sets
    REAL(dp), INTENT(OUT), OPTIONAL          :: r_max
    INTEGER, INTENT(OUT), OPTIONAL           :: nlist
    LOGICAL, INTENT(OUT), OPTIONAL           :: symmetric, mic, molecular

    CHARACTER(len=*), PARAMETER :: routineN = 'get_neighbor_list_set_p', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set

    IF (ASSOCIATED(neighbor_list_sets)) THEN

      NULLIFY(neighbor_list_set)
      DO i = 1, SIZE(neighbor_list_sets)
         neighbor_list_set => neighbor_list_sets(i)%neighbor_list_set
         IF (ASSOCIATED(neighbor_list_set)) EXIT
      END DO

      IF (ASSOCIATED(neighbor_list_set)) THEN
         IF (PRESENT(r_max)) r_max = neighbor_list_set%r_max
         IF (PRESENT(nlist)) nlist = neighbor_list_set%nlist
         IF (PRESENT(symmetric)) symmetric = neighbor_list_set%symmetric
         IF (PRESENT(mic)) mic = neighbor_list_set%mic
         IF (PRESENT(molecular)) molecular = neighbor_list_set%molecular
      ELSE
         CALL stop_program(routineN,moduleN,__LINE__,&
                        "No neighbor list set is associated")
      END IF

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The requested neighbor list sets are not associated")

    END IF

  END SUBROUTINE get_neighbor_list_set_p

! *****************************************************************************
!> \brief   Return the reqested data of a neighbor node.
!> \param neighbor_node ...
!> \param neighbor ...
!> \param cell ...
!> \param r ...
!> \date    23.06.2000
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_neighbor_node(neighbor_node,neighbor,cell,r)

    TYPE(neighbor_node_type), POINTER        :: neighbor_node
    INTEGER, INTENT(OUT), OPTIONAL           :: neighbor
    INTEGER, DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: cell
    REAL(dp), DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: r

    CHARACTER(len=*), PARAMETER :: routineN = 'get_neighbor_node', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(neighbor_node)) THEN

      IF (PRESENT(neighbor)) neighbor = neighbor_node%neighbor
      IF (PRESENT(r)) r(:) = neighbor_node%r(:)
      IF (PRESENT(cell)) cell(:) = neighbor_node%cell(:)

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The requested neighbor node is not associated")

    END IF

  END SUBROUTINE get_neighbor_node

! *****************************************************************************
!> \brief Initialize a neighbor list set. Nothing is (de)allocated here.
!>         This routine is also used to prepare a neighbor list set for
!>         overwriting.
!> \param neighbor_list_set ...
!> \param r_max ...
!> \param symmetric ...
!> \param mic ...
!> \param molecular ...
!> \date  20.09.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_neighbor_list_set(neighbor_list_set,r_max,symmetric,mic,molecular)

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    REAL(dp), INTENT(IN)                     :: r_max
    LOGICAL, INTENT(IN)                      :: symmetric, mic, molecular

    CHARACTER(len=*), PARAMETER :: routineN = 'init_neighbor_list_set', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(neighbor_list_set)) THEN

      ! *** Initialize the pointers to the last neighbor list ***
      NULLIFY (neighbor_list_set%last_neighbor_list)

      ! *** Initialize the neighbor list counter ***
      neighbor_list_set%nlist = 0

      ! *** Initialize the maximum interaction radius (optional) ***
      neighbor_list_set%r_max = r_max

      ! *** Initialize the neighbor list build properties
      neighbor_list_set%symmetric = symmetric
      neighbor_list_set%mic       = mic
      neighbor_list_set%molecular = molecular

    ELSE

      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The requested neighbor list set is not associated")

    END IF

  END SUBROUTINE init_neighbor_list_set

END MODULE qs_neighbor_list_types
