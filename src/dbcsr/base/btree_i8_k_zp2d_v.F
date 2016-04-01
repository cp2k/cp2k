!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   B-tree implementation template
!> \author  Urban Borstnik
!> \date    2000-05-12
!> \version 1.0
!> <b>Modification history:</b>
!> - Created
! **************************************************************************************************

MODULE btree_i8_k_zp2d_v

   IMPLICIT NONE
   PRIVATE

   !API
   PUBLIC :: btree, btree_new, btree_add, btree_find, &
             btree_delete, btree_get_entries

   INTEGER, PARAMETER :: keyt = SELECTED_INT_KIND(10)
   INTEGER, PARAMETER :: valt = SELECTED_INT_KIND(5); 
   INTEGER, PARAMETER :: dp = KIND(0.0d0)

   TYPE zp2d
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER :: p
      LOGICAL :: tr
   END TYPE zp2d

   PUBLIC :: zp2d

   TYPE btree_node
      INTEGER id
      INTEGER :: filled
      INTEGER(KIND=keyt), DIMENSION(:), POINTER :: keys
      TYPE(zp2d), DIMENSION(:), POINTER :: values
      TYPE(btree_node_p), DIMENSION(:), POINTER :: subtrees
      TYPE(btree_node), POINTER :: parent
   END TYPE btree_node

   TYPE btree_node_p
      TYPE(btree_node), POINTER :: node
   END TYPE btree_node_p

   TYPE btree_node_structure
      INTEGER :: min_fill, max_fill
      INTEGER :: n
      INTEGER :: lastid
      INTEGER :: refcount
      TYPE(btree_node), POINTER :: root
   END TYPE btree_node_structure

   TYPE btree
      TYPE(btree_node_structure) :: b
   END TYPE btree

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE btree_new(tree, order)
      TYPE(btree), INTENT(OUT)                           :: tree
      INTEGER, INTENT(IN), OPTIONAL                      :: order

      INTEGER                                            :: maxs, mins

!

      IF (PRESENT(order)) THEN
         maxs = order-1
      ELSE
         maxs = 15
      ENDIF
      mins = ISHFT(maxs, -1)
      IF (mins*2 .GT. maxs) maxs = 2*maxs
      IF (mins .LT. 1) mins = 1
      IF (maxs .LT. 3) maxs = 3
      tree%b%min_fill = mins
      tree%b%max_fill = maxs
      tree%b%refcount = 1
      tree%b%n = 0
      NULLIFY (tree%b%root)
      tree%b%lastid = 0
   END SUBROUTINE btree_new

! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \retval num_entries ...
! **************************************************************************************************
   FUNCTION btree_get_entries(tree) RESULT(num_entries)
      TYPE(btree), INTENT(INOUT)                         :: tree
      INTEGER                                            :: num_entries

      num_entries = tree%b%n
   END FUNCTION btree_get_entries

! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \param keys ...
!> \param values ...
! **************************************************************************************************
   SUBROUTINE btree_delete(tree, keys, values)
      TYPE(btree), INTENT(INOUT)                         :: tree
      INTEGER(KIND=keyt), DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: keys
      TYPE(zp2d), DIMENSION(:), INTENT(OUT), OPTIONAL    :: values

      INTEGER                                            :: pos

!

      IF (ASSOCIATED(tree%b%root)) THEN
         pos = 0
         IF (PRESENT(keys) .AND. PRESENT(values)) THEN
            pos = 1
            CALL btree_delete_node(tree%b%root, pos, keys, values)
         ELSE
            CALL btree_delete_node(tree%b%root)
         ENDIF
      ENDIF
      NULLIFY (tree%b%root)
   END SUBROUTINE btree_delete

! **************************************************************************************************
!> \brief ...
!> \param node ...
!> \param pos ...
!> \param keys ...
!> \param values ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE btree_delete_node(node, pos, keys, values)
      TYPE(btree_node), POINTER                          :: node
      INTEGER, INTENT(INOUT), OPTIONAL                   :: pos
      INTEGER(KIND=keyt), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: keys
      TYPE(zp2d), DIMENSION(:), INTENT(INOUT), OPTIONAL  :: values

      INTEGER                                            :: i

!
!

      IF (node%filled .GT. 0 .AND. ASSOCIATED(node%subtrees(1)%node)) THEN
         DO i = 1, node%filled+1
            IF (PRESENT(pos)) THEN
               CALL btree_delete_node(node%subtrees(i)%node, pos, keys, values)
            ELSE
               CALL btree_delete_node(node%subtrees(i)%node)
            ENDIF
            IF (PRESENT(pos) .AND. i .LE. node%filled) THEN
               keys(pos) = node%keys(i)
               values(pos) = node%values(i)
               pos = pos+1
            ENDIF
         ENDDO
      ELSEIF (PRESENT(pos) .AND. node%filled .GT. 0) THEN
         keys(pos:pos+node%filled-1) = node%keys(1:node%filled)
         values(pos:pos+node%filled-1) = node%values(1:node%filled)
         pos = pos+node%filled
      ENDIF
      CALL btree_free_node(node)
   END SUBROUTINE btree_delete_node

   ! Find the key
   ! IF node still has space, insert & update the node
   ! else
   ! 1. select median
   ! 2. split keys into two nodes (one is new)
   ! 3. insert separation key put into parent, and repeat upwards
! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \param key ...
!> \param value ...
!> \param exists ...
!> \param existing_value ...
!> \param replace ...
! **************************************************************************************************
   SUBROUTINE btree_add(tree, key, value, exists, existing_value, replace)
      TYPE(btree), INTENT(INOUT)                         :: tree
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      TYPE(zp2d), INTENT(IN)                             :: value
      LOGICAL, INTENT(OUT), OPTIONAL                     :: exists
      TYPE(zp2d), INTENT(OUT), OPTIONAL                  :: existing_value
      LOGICAL, INTENT(IN), OPTIONAL                      :: replace

      INTEGER                                            :: ge_pos, position
      TYPE(btree_node), POINTER                          :: node

!
!

      IF (PRESENT(exists)) THEN
         CALL btree_find_full(tree, key, node, position, ge_pos, short=.TRUE.)
         IF (position .GT. 0) THEN
            exists = .TRUE.
            existing_value = node%values(position)
            IF (PRESENT(replace)) THEN
               IF (replace) THEN
                  node%values(position) = value
               ENDIF
            ENDIF
            RETURN
         ELSE
            exists = .FALSE.
         ENDIF
      ELSE
         CALL btree_find_leaf(tree, key, node, ge_pos)
      ENDIF
      CALL btree_add_into(tree, node, key, value, before=ge_pos)
      IF (PRESENT(exists)) existing_value = value
      tree%b%n = tree%b%n+1
   END SUBROUTINE btree_add

! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \param node ...
!> \param key ...
!> \param value ...
!> \param before ...
!> \param subtree ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE btree_add_into(tree, node, key, value, before, subtree)
      TYPE(btree), INTENT(INOUT)                         :: tree
      TYPE(btree_node), POINTER                          :: node
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      TYPE(zp2d), INTENT(IN)                             :: value
      INTEGER, INTENT(IN), OPTIONAL                      :: before
      TYPE(btree_node), OPTIONAL, POINTER                :: subtree

      INTEGER                                            :: ge_pos, split_pos
      INTEGER(KIND=keyt)                                 :: upgrade_key
      LOGICAL                                            :: leaf
      TYPE(btree_node), POINTER                          :: new_node
      TYPE(zp2d)                                         :: upgrade_value

!
!
! Root is special

      IF (.NOT. ASSOCIATED(node)) THEN
         CALL btree_new_root(tree, key, value)
         IF (PRESENT(subtree)) THEN
            tree%b%root%subtrees(2)%node => subtree
            subtree%parent => tree%b%root
         ENDIF
         RETURN
      ENDIF
      ! Where the insertion takes place.
      IF (PRESENT(before)) THEN
         ge_pos = before
      ELSE
         CALL btree_node_find_gt_pos(node%keys, key, ge_pos, node%filled)
      ENDIF
      ! Addition is easy if the node has enough space.
      leaf = .NOT. ASSOCIATED(node%subtrees(1)%node)
      IF (node%filled .LT. tree%b%max_fill) THEN
         IF (PRESENT(subtree)) THEN
            CALL btree_simple_insertion(node, key, value, ge_pos, subtree)
         ELSE
            CALL btree_simple_insertion(node, key, value, ge_pos)
         ENDIF
         RETURN
      ELSE
         split_pos = ISHFT(tree%b%max_fill+1, -1)
         ! I assert that split_pos <= SIZE(node%keys)
         CALL btree_new_node(tree, new_node)
         ! The key to be added falls in the left node
         node%filled = split_pos-1
         IF (ge_pos .LE. split_pos) THEN
            IF (ge_pos .EQ. split_pos) THEN
               upgrade_key = key
               upgrade_value = value
            ELSE
               upgrade_key = node%keys(split_pos-1)
               upgrade_value = node%values(split_pos-1)
            ENDIF
            IF (PRESENT(subtree)) THEN
               CALL btree_left_insertion(tree, node, new_node, key, value, &
                                         ge_pos, split_pos, subtree)
               !CALL btree_adopt_subtrees (new_node)
            ELSE
               CALL btree_left_insertion(tree, node, new_node, key, value, &
                                         ge_pos, split_pos)
            ENDIF
            !
         ELSE
            upgrade_key = node%keys(split_pos)
            upgrade_value = node%values(split_pos)
            IF (PRESENT(subtree)) THEN
               CALL btree_right_insertion(tree, node, new_node, key, value, &
                                          ge_pos, split_pos, subtree)
               !CALL btree_adopt_subtrees (new_node)
            ELSE
               CALL btree_right_insertion(tree, node, new_node, key, value, &
                                          ge_pos, split_pos)
            ENDIF
            !
         ENDIF
         !
         new_node%parent => node%parent
         !
         IF (.NOT. leaf) THEN
            CALL btree_adopt_subtrees(new_node)
         ENDIF
         !
         CALL btree_add_into(tree, node%parent, upgrade_key, upgrade_value, &
                             subtree=new_node)
         !
      ENDIF
   END SUBROUTINE btree_add_into

! **************************************************************************************************
!> \brief ...
!> \param node ...
!> \param key ...
!> \param value ...
!> \param before ...
!> \param subtree ...
! **************************************************************************************************
   SUBROUTINE btree_simple_insertion(node, key, value, before, subtree)
      TYPE(btree_node), INTENT(INOUT)                    :: node
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      TYPE(zp2d), INTENT(IN)                             :: value
      INTEGER, INTENT(IN)                                :: before
      TYPE(btree_node), OPTIONAL, POINTER                :: subtree

!
! Shift keys

      node%keys(before+1:node%filled+1) = node%keys(before:node%filled)
      node%keys(before) = key
      ! Shift values
      node%values(before+1:node%filled+1) = node%values(before:node%filled)
      node%values(before) = value
      ! Shift subtree pointers, but only if node is not a leaf ; assume
      ! leaf <=> present(subtree)
      IF (PRESENT(subtree)) THEN
         node%subtrees(before+2:node%filled+2) = &
            node%subtrees(before+1:node%filled+1)
         node%subtrees(before+1)%node => subtree
      ENDIF
      node%filled = node%filled+1
   END SUBROUTINE btree_simple_insertion

! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \param node ...
!> \param new_node ...
!> \param key ...
!> \param value ...
!> \param before ...
!> \param split_pos ...
!> \param subtree ...
! **************************************************************************************************
   SUBROUTINE btree_left_insertion(tree, node, new_node, key, value, before, split_pos, subtree)
      TYPE(btree), INTENT(IN)                            :: tree
      TYPE(btree_node), INTENT(INOUT)                    :: node, new_node
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      TYPE(zp2d), INTENT(IN)                             :: value
      INTEGER, INTENT(IN)                                :: before, split_pos
      TYPE(btree_node), OPTIONAL, POINTER                :: subtree

!

      new_node%filled = (tree%b%max_fill)-(split_pos-1)
      new_node%keys(1:new_node%filled) = &
         node%keys(split_pos:tree%b%max_fill)
      new_node%values(1:new_node%filled) = &
         node%values(split_pos:tree%b%max_fill)
      !IF (ASSOCIATED (node%subtrees(1)%node)) THEN
      IF (PRESENT(subtree)) THEN
         IF (before .EQ. split_pos) THEN
            new_node%subtrees(2:new_node%filled+1) = &
               node%subtrees(split_pos+1:tree%b%max_fill+1)
            new_node%subtrees(1)%node => subtree
         ELSE
            new_node%subtrees(1:new_node%filled+1) = &
               node%subtrees(split_pos:tree%b%max_fill+1)
         ENDIF
      ENDIF
      ! Fill node%{keys,values}(1:node%filled), where node%filled
      ! is split_pos-1, but do insert the new value at ge_pos. The
      ! key/value at split_pos is to be inserted into the
      ! parent.
      ! The new tree is added to the right of the new insertion.
      node%keys(before+1:node%filled) = node%keys(before:node%filled-1)
      node%keys(before) = key
      node%values(before+1:node%filled) = node%values(before:node%filled-1)
      node%values(before) = value
      IF (PRESENT(subtree)) THEN
         node%subtrees(before+2:node%filled+1) = &
            node%subtrees(before+1:node%filled)
         node%subtrees(before+1)%node => subtree
      ELSE
         NULLIFY (node%subtrees(before+1)%node)
      ENDIF
   END SUBROUTINE btree_left_insertion

! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \param node ...
!> \param new_node ...
!> \param key ...
!> \param value ...
!> \param before ...
!> \param split_pos ...
!> \param subtree ...
! **************************************************************************************************
   SUBROUTINE btree_right_insertion(tree, node, new_node, key, value, before, split_pos, subtree)
      TYPE(btree), INTENT(IN)                            :: tree
      TYPE(btree_node), INTENT(INOUT)                    :: node, new_node
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      TYPE(zp2d), INTENT(IN)                             :: value
      INTEGER, INTENT(IN)                                :: before, split_pos
      TYPE(btree_node), OPTIONAL, POINTER                :: subtree

!

      new_node%filled = (tree%b%max_fill+1)-split_pos
      new_node%keys(1:before-split_pos-1) = &
         node%keys(split_pos+1:before-1)
      new_node%keys(before-split_pos) = key
      new_node%keys(before-split_pos+1:new_node%filled) = &
         node%keys(before:tree%b%max_fill)
      new_node%values(1:before-split_pos-1) = &
         node%values(split_pos+1:before-1)
      new_node%values(before-split_pos) = value
      new_node%values(before-split_pos+1:new_node%filled) = &
         node%values(before:tree%b%max_fill)
      IF (PRESENT(subtree)) THEN
         new_node%subtrees(1:before-split_pos) = &
            node%subtrees(split_pos+1:before)
         new_node%subtrees(before-split_pos+1)%node => subtree
         new_node%subtrees(before-split_pos+2:new_node%filled+1) = &
            node%subtrees(before+1:tree%b%max_fill+1)
      ENDIF
   END SUBROUTINE btree_right_insertion

   ! node is a non-leaf node
! **************************************************************************************************
!> \brief ...
!> \param node ...
! **************************************************************************************************
   SUBROUTINE btree_adopt_subtrees(node)
      TYPE(btree_node), POINTER                          :: node

      INTEGER                                            :: i

!
! Assume that node is not a leaf!

      DO i = 1, node%filled+1
         !IF (ASSOCIATED (node%subtrees(i)%node)) THEN
         !IF (.NOT. ASSOCIATED (node%subtrees(i)%node%parent,&
         ! node)) THEN
         node%subtrees(i)%node%parent => node
         !ENDIF
         !ENDIF
      ENDDO
   END SUBROUTINE btree_adopt_subtrees

! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \param key ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE btree_new_root(tree, key, value)
      TYPE(btree), INTENT(INOUT)                         :: tree
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      TYPE(zp2d), INTENT(IN)                             :: value

      TYPE(btree_node), POINTER                          :: new_root, old_root

!

      CALL btree_new_node(tree, new_root)
      new_root%filled = 1
      new_root%keys(1) = key
      new_root%values(1) = value
      IF (ASSOCIATED(tree%b%root)) THEN
         old_root => tree%b%root
         old_root%parent => new_root
         new_root%subtrees(1)%node => old_root
         old_root%parent => new_root
      ENDIF
      tree%b%root => new_root
   END SUBROUTINE btree_new_root

! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \param node ...
! **************************************************************************************************
   SUBROUTINE btree_new_node(tree, node)
      TYPE(btree), INTENT(INOUT)                         :: tree
      TYPE(btree_node), POINTER                          :: node

      INTEGER                                            :: i

!

      ALLOCATE (node)
      ALLOCATE (node%keys(tree%b%max_fill))
      ALLOCATE (node%values(tree%b%max_fill))
      ALLOCATE (node%subtrees(tree%b%max_fill+1))
      DO i = 1, tree%b%max_fill+1
         NULLIFY (node%subtrees(i)%node)
      ENDDO
      node%filled = 0
      NULLIFY (node%parent)
      tree%b%lastid = tree%b%lastid+1
      node%id = tree%b%lastid
   END SUBROUTINE btree_new_node

! **************************************************************************************************
!> \brief ...
!> \param node ...
! **************************************************************************************************
   SUBROUTINE btree_free_node(node)
      TYPE(btree_node), POINTER                          :: node

!

      DEALLOCATE (node%keys)
      DEALLOCATE (node%values)
      DEALLOCATE (node%subtrees)
      DEALLOCATE (node)
   END SUBROUTINE btree_free_node

! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \param key ...
!> \param value ...
!> \param exists ...
! **************************************************************************************************
   SUBROUTINE btree_find(tree, key, value, exists)
      TYPE(btree), INTENT(IN)                            :: tree
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      TYPE(zp2d), INTENT(OUT)                            :: value
      LOGICAL, INTENT(OUT), OPTIONAL                     :: exists

      INTEGER                                            :: position
      TYPE(btree_node), POINTER                          :: node

!
!

      CALL btree_find_full(tree, key, node, position, short=.TRUE.)
      IF (PRESENT(exists)) THEN
         exists = position .GT. 0
      ENDIF
      IF (position .GT. 0) THEN
         value = node%values(position)
      ENDIF
   END SUBROUTINE btree_find

! **************************************************************************************************
!> \brief ...
!> \param keys ...
!> \param key ...
!> \param position ...
!> \param filled ...
! **************************************************************************************************
   SUBROUTINE btree_node_find_ge_pos(keys, key, position, filled)
      INTEGER(KIND=keyt), DIMENSION(:)                   :: keys
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      INTEGER, INTENT(OUT)                               :: position
      INTEGER, INTENT(IN)                                :: filled

      INTEGER                                            :: left, right

!

      IF (keys(1) .GE. key) THEN
         position = 1
         RETURN
      ENDIF
      IF (keys(filled) .LT. key) THEN
         position = filled+1
         RETURN
      ENDIF
      left = 2
      right = filled
      position = MAX(ISHFT(left+right, -1), left)
      DO WHILE (left .LE. right)
         IF (keys(position) .GE. key .AND. keys(position-1) .LT. key) THEN
            RETURN
         ENDIF
         IF (keys(position) .GE. key) right = MIN(position, right-1)
         IF (keys(position) .LT. key) left = MAX(position, left+1)
         position = MAX(ISHFT(left+right, -1), left)
      ENDDO
   END SUBROUTINE btree_node_find_ge_pos

! **************************************************************************************************
!> \brief ...
!> \param keys ...
!> \param key ...
!> \param position ...
!> \param filled ...
! **************************************************************************************************
   SUBROUTINE btree_node_find_gt_pos(keys, key, position, filled)
      INTEGER(KIND=keyt), DIMENSION(:)                   :: keys
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      INTEGER, INTENT(OUT)                               :: position
      INTEGER, INTENT(IN)                                :: filled

      INTEGER                                            :: left, right

!

      IF (keys(1) .GT. key) THEN
         position = 1
         RETURN
      ENDIF
      IF (keys(filled) .LE. key) THEN
         position = filled+1
         RETURN
      ENDIF
      left = 2
      right = filled
      position = MAX(ISHFT(left+right, -1), left)
      DO WHILE (left .LE. right)
         IF (keys(position) .GT. key .AND. keys(position-1) .LE. key) THEN
            RETURN
         ENDIF
         IF (keys(position) .GT. key) right = MIN(position, right-1)
         IF (keys(position) .LE. key) left = MAX(position, left+1)
         position = MAX(ISHFT(left+right, -1), left)
      ENDDO
   END SUBROUTINE btree_node_find_gt_pos
! **************************************************************************************************
!> \brief ...
!> \param keys ...
!> \param key ...
!> \param position ...
!> \param filled ...
!> \param first ...
! **************************************************************************************************
   SUBROUTINE btree_node_find_gte_pos(keys, key, position, filled, first)
      INTEGER(KIND=keyt), DIMENSION(:)                   :: keys
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      INTEGER, INTENT(OUT)                               :: position
      INTEGER, INTENT(IN)                                :: filled
      INTEGER, INTENT(IN), OPTIONAL                      :: first

      INTEGER                                            :: left, one, right

!

      one = 1
      IF (PRESENT(FIRST)) one = first
      IF (one .LE. filled) THEN
         IF (keys(one) .GT. key) THEN
            position = one
            RETURN
         ENDIF
      ENDIF
      IF (keys(filled) .LE. key) THEN
         position = filled+1
         RETURN
      ENDIF
      left = one+1
      right = filled
      position = MAX(ISHFT(left+right, -1), left)
      DO WHILE (left .LE. right)
         IF (keys(position) .GT. key .AND. keys(position-1) .LE. key) THEN
            RETURN
         ENDIF
         IF (keys(position) .GT. key) right = MIN(position, right-1)
         IF (keys(position) .LE. key) left = MAX(position, left+1)
         position = MAX(ISHFT(left+right, -1), left)
      ENDDO
   END SUBROUTINE btree_node_find_gte_pos

   ! node is unassociated and position=0 if not found
   ! Precondition: The key is tree or its subtree.
! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \param key ...
!> \param node ...
!> \param position ...
!> \param ge_position ...
!> \param short ...
! **************************************************************************************************
   SUBROUTINE btree_find_full(tree, key, node, position, ge_position, short)
      TYPE(btree), INTENT(IN)                            :: tree
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      TYPE(btree_node), POINTER                          :: node
      INTEGER, INTENT(OUT)                               :: position
      INTEGER, INTENT(OUT), OPTIONAL                     :: ge_position
      LOGICAL, INTENT(IN), OPTIONAL                      :: short

      INTEGER                                            :: gti
      LOGICAL                                            :: stop_short

! Used mark searches
!

      stop_short = .FALSE.
      IF (PRESENT(short)) stop_short = short
      NULLIFY (node)
      position = 0
      IF (PRESENT(ge_position)) ge_position = 0
      !IF (tree%b%n .EQ. 0) RETURN
      IF (.NOT. ASSOCIATED(tree%b%root)) RETURN
      gti = 1
      ! Try to find the key in the given node. If it's found, then
      ! return the node.
      node => tree%b%root
      descent: DO WHILE (.TRUE.)
         ! Try to find the first element equal to or greater than the
         ! one we're searching for.
         CALL btree_node_find_ge_pos(node%keys, key, position, node%filled)
         ! One of three things is now true about position: it's now
         ! greater than the number of keys (if all keys are smaller), or
         ! it points to the key that is equal to or greater than the one
         ! we are searching for. If it is found and we are just
         ! searching for one equal element (i.e., user search), we can
         ! return.
         IF (stop_short .AND. position .LE. node%filled) THEN
            IF (node%keys(position) .EQ. key) THEN
               IF (PRESENT(ge_position)) ge_position = position
               RETURN
            ENDIF
         ENDIF
         ! If the key is not found, then either return the GE position
         ! if we're in a leaf (case 2 here), otherwise descend into the
         ! subtrees.
         !CALL btree_node_find_gt_pos (node%keys, key, gti, node%filled, position)
         CALL btree_node_find_gte_pos(node%keys, key, gti, node%filled, position)
         IF (ASSOCIATED(node%subtrees(1)%node)) THEN
            node => node%subtrees(gti)%node
         ELSE
            IF (PRESENT(ge_position)) ge_position = gti
            position = 0
            RETURN
         ENDIF
      END DO descent
   END SUBROUTINE btree_find_full

   ! node is unassociated and position=0 if not found
   ! Precondition: The key is tree or its subtree.
! **************************************************************************************************
!> \brief ...
!> \param tree ...
!> \param key ...
!> \param node ...
!> \param gti ...
! **************************************************************************************************
   SUBROUTINE btree_find_leaf(tree, key, node, gti)
      TYPE(btree), INTENT(IN)                            :: tree
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      TYPE(btree_node), POINTER                          :: node
      INTEGER, INTENT(OUT)                               :: gti

!

      NULLIFY (node)
      !IF (tree%b%n .EQ. 0) RETURN
      IF (.NOT. ASSOCIATED(tree%b%root)) RETURN
      gti = 1
      ! Try to find the key in the given node. If it's found, then
      ! return the node.
      node => tree%b%root
      descent: DO WHILE (.TRUE.)
         ! Try to find the first element equal to or greater than the
         ! one we're searching for.
         !CALL btree_node_find_ge_pos (node%keys, key, position, node%filled)
         ! One of three things is now true about position: it's now
         ! greater than the number of keys (if all keys are smaller), or
         ! it points to the key that is equal to or greater than the one
         ! we are searching for. If it is found and we are just
         ! searching for one equal element (i.e., user search), we can
         ! return.
         !
         ! If the key is not found, then either return the GE position
         ! if we're in a leaf (case 2 here), otherwise descend into the
         ! subtrees.
         CALL btree_node_find_gt_pos(node%keys, key, gti, node%filled)
         !CALL btree_node_find_gt2_pos (node%keys, key, i, node%filled)
         !IF (i .NE. gti) WRITE(*,*)'XXXX difference',i,gti
         IF (ASSOCIATED(node%subtrees(1)%node)) THEN
            node => node%subtrees(gti)%node
         ELSE
            RETURN
         ENDIF
      END DO descent
   END SUBROUTINE btree_find_leaf

END MODULE btree_i8_k_zp2d_v
