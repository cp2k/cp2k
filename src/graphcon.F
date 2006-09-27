!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****h* cp2k/graphcon *
!!
!!   NAME
!!     graphcon
!!
!!   FUNCTION
!!     uses a combination of graphs and hashing to determine if two molecules
!!     are topologically equivalent, and if so, finds the one by one mapping
!!
!!   NOTES
!!     the graph isomorphism being solved is a computationally hard one
!!     and can not be solved in polynomial time in the general case
!!     http://mathworld.wolfram.com/IsomorphicGraphs.html
!!     the problem arises if many atoms are topologically equivalent
!!     the current algorithm is able to solve the problem for benzene (C6H6)
!!     but not for a fullerene (C60). Large systems are not really a problem (JAC).
!!     as almost all atoms are topologically unique.
!!
!!   AUTHOR
!!     Joost VandeVondele
!!     
!!
!!   MODIFICATION HISTORY
!!     09.2006 [Joost VandeVondele]
!!
!!   SOURCE
!****************************************************************************
MODULE graphcon
    USE util, ONLY: sort
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: vertex, reorder_graph, hash_molecule

    ! a molecule is an array of vertices, each vertex has a kind
    ! and a list of edges (bonds). 
    ! (the number is the index of the other vertex in the array that builds the molecule)
    TYPE vertex
       INTEGER :: kind
       INTEGER, POINTER, DIMENSION(:) :: bonds
    END TYPE vertex

    TYPE class
       INTEGER, DIMENSION(:), POINTER :: reference
       INTEGER, DIMENSION(:), POINTER :: unordered
       INTEGER :: kind
       INTEGER :: Nele
       LOGICAL :: first
       INTEGER, DIMENSION(:), POINTER :: order
       INTEGER, DIMENSION(:), POINTER :: q
    END TYPE class

    TYPE superclass
       INTEGER :: Nele
       INTEGER, DIMENSION(:), POINTER :: classes
    END TYPE

CONTAINS
!!****f* graphcon/hash_molecule *
!!
!!   NAME
!!     hash_molecule
!!
!!   FUNCTION
!!     hashes a molecule to a number. Molecules that are the (topologically) the same
!!     have the same hash. However, there is a small chance that molecules with the same hash
!!     are different
!!
!!   NOTES
!!     Although relatively fast in general, might be quadratic with molecule size for
!!     some systems (e.g. linear alkanes)
!!
!!   INPUTS
!!     reference : IN  : molecule with atomic kinds and bonds
!!     kind_ref  : OUT : an atomic hash which is the same for topologically equivalent atoms
!!     hash      : OUT : a hash which is the same for topologically equivalent molecules
!!
!!   MODIFICATION HISTORY
!!     09.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
    SUBROUTINE hash_molecule(reference,kind_ref,hash)
        TYPE(vertex), DIMENSION(:),INTENT(IN)  :: reference
        INTEGER, DIMENSION(:),INTENT(OUT)      :: kind_ref
        INTEGER, INTENT(OUT)                   :: hash
        INTEGER            :: Ihash,N,I,old_class,Nclasses,Nclasses_old
        INTEGER, DIMENSION(:), ALLOCATABLE :: kind_new, index
        N=SIZE(kind_ref)
        ALLOCATE(kind_new(N),index(N))
        kind_ref=reference%kind
        Nclasses_old=0
        DO Ihash=1,N
           ! generate a hash based on the the kind of each atom and the kind of its bonded atoms
           DO I=1,N
              kind_new(I)=hash_kind(kind_ref(I),kind_ref(reference(I)%bonds))
           ENDDO
           kind_ref=kind_new
           ! find the number of equivalent atoms
           CALL sort(kind_new,N,index)
           Nclasses=1 
           old_class=kind_new(1)
           DO i=2,N
              IF (kind_new(I).NE.old_class) THEN
                 Nclasses=Nclasses+1
                 old_class=kind_new(I)
              ENDIF
           ENDDO
           ! if we have not generated new classes, we have presumably found all equivalence classes
           IF (Nclasses==Nclasses_old) EXIT
           Nclasses_old=Nclasses
           ! write(6,*) "Classes",Ihash, Nclasses
        ENDDO
        ! hash (sorted) kinds to a molecular hash
        hash=joaat_hash_i(kind_new)
        DEALLOCATE(kind_new,index)
    END SUBROUTINE hash_molecule

!!****f* graphcon/reorder_graph *
!!
!!   NAME
!!     reorder_graph
!!
!!   FUNCTION
!!     If two molecules are topologically the same, finds the ordering that maps 
!!     the unordered one on the ordered one.
!!
!!   NOTES
!!     See not at the top of the file about why this algorithm might consider
!!     molecules with a large number of equivalent atoms as different
!!     despite the fact that an ordering could exist for which they are the same
!!     
!!
!!   INPUTS
!!     reference, unordered : molecular description (see type definition) 
!!     order                : the mapping reference=order(unordred) if matches=.TRUE.
!!                            undefined if matches=.FALSE.
!!     matches              : .TRUE. = the ordering was found
!!
!!   MODIFICATION HISTORY
!!     09.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
    SUBROUTINE reorder_graph(reference, unordered, order, matches)
        TYPE(vertex), DIMENSION(:),INTENT(IN)  :: reference
        TYPE(vertex), DIMENSION(:),INTENT(IN)  :: unordered
        INTEGER, DIMENSION(:),INTENT(OUT)     :: order
        LOGICAL, INTENT(OUT)                  :: matches
        ! allows for worst case matching of two benzenes ... (6!)*(6!)/6=86400
        ! with some margin for other molecules
        ! molecules with no symmetry e.g. JAC need less than 500 tries
        INTEGER, PARAMETER                     :: max_tries=1000000
        LOGICAL :: try_it
        INTEGER :: I,N,Nclasses,Iclass,Nele,Ihash,J,old_class,hash_re,hash_un,isuperclass,iele,itries
        TYPE(class), DIMENSION(:), ALLOCATABLE :: classes
        TYPE(superclass), DIMENSION(:), ALLOCATABLE :: superclasses
        INTEGER, DIMENSION(:), ALLOCATABLE :: kind_ref,kind_un,index_ref,index_un, &
                                              kind_ref_ordered,kind_un_ordered, &
                                              class_of_atom, superclass_of_atom
        REAL :: t2,t1

        ! catch the cases where the molecules are trivially different 
        IF (SIZE(reference).NE.SIZE(unordered)) THEN
           matches=.FALSE.
           return
        ENDIF

        ! catch the case where the molecules are already in the right order
        N=SIZE(order)
        order=(/(i,i=1,N)/)
        IF (matrix_equal(reference,unordered,order)) THEN
           matches=.TRUE.
           return
        ENDIF

        ! determine the kind of each atom, and the hash of the whole molecule
        ALLOCATE(kind_ref(N),kind_un(N),index_ref(N),index_un(N), &
                 kind_ref_ordered(N),kind_un_ordered(N), & 
                 class_of_atom(N),superclass_of_atom(N))
        CALL hash_molecule(reference,kind_ref,hash_re)
        CALL hash_molecule(unordered,kind_un,hash_un)
        IF (hash_re .NE. hash_un) THEN
           matches=.FALSE.
           return
        ENDIF

        ! generate the classes of equivalent atoms, i.e. the groups of atoms of the same topological kind
        kind_ref_ordered=kind_ref
        CALL sort(kind_ref_ordered,N,index_ref)
        kind_un_ordered=kind_un
        CALL sort(kind_un_ordered,N,index_un)
        IF (ANY(kind_ref_ordered.NE.kind_un_ordered)) THEN
           matches=.FALSE.
           return
        ENDIF

        ! count different classes, assign their kinds, and the number of elements
        Nclasses=1 
        old_class=kind_ref_ordered(1)
        DO i=2,N
           IF (kind_ref_ordered(I).NE.old_class) THEN
              Nclasses=Nclasses+1
              old_class=kind_ref_ordered(I)
           ENDIF
        ENDDO
        ALLOCATE(classes(Nclasses))
        classes(1)%kind=kind_ref_ordered(1)
        Nclasses=1
        classes(1)%Nele=1
        DO i=2,N
           IF (kind_ref_ordered(I).NE.classes(Nclasses)%kind) THEN
              Nclasses=Nclasses+1
              classes(Nclasses)%kind=kind_ref_ordered(I)
              classes(Nclasses)%Nele=1
           ELSE 
              classes(Nclasses)%Nele=classes(Nclasses)%Nele+1
           ENDIF
        ENDDO

        ! assign the atoms to their classes
        iele=0
        DO I=1,Nclasses
           Nele=classes(I)%Nele
           ALLOCATE(classes(I)%reference(Nele))
           ALLOCATE(classes(I)%unordered(Nele))
           DO J=1,Nele
              iele=iele+1
              classes(I)%reference(J)=index_ref(iele)
              classes(I)%unordered(J)=index_un(iele)
           ENDDO
           class_of_atom(classes(I)%reference)=I
           ALLOCATE(classes(I)%order(Nele))
           ALLOCATE(classes(I)%q(Nele)) 
           classes(I)%order=(/(J,J=1,Nele)/)
           classes(I)%first=.TRUE.
        ENDDO

        ! find which groups of classes (superclasses) that can be solved independently. 
        ! only classes with more than one element that are connected need to be reordered simultaniously

        ! find these connected components in a recursive way
        superclass_of_atom=-1
        isuperclass=0
        DO I=1,N
           ! this atom belongs to a class with several equivalent atoms, and has not yet been found
           IF (superclass_of_atom(I).EQ.-1 .AND. classes(class_of_atom(I))%Nele>1) THEN
              isuperclass=isuperclass+1
              CALL spread_superclass(I,isuperclass,superclass_of_atom,class_of_atom,classes,reference)
           ENDIF
        ENDDO

        ! put classes into superclasses
        ALLOCATE(superclasses(isuperclass))
        superclasses%Nele=0
        DO I=1,Nclasses 
           J=superclass_of_atom(classes(I)%reference(1))
           IF (J>0) superclasses(J)%Nele=superclasses(J)%Nele+1
        ENDDO
        DO I=1,isuperclass
           ALLOCATE(superclasses(I)%classes(superclasses(I)%Nele))
           superclasses(I)%Nele=0
        ENDDO
        DO I=1,Nclasses
           J=superclass_of_atom(classes(I)%reference(1))
           IF (J>0) THEN
               superclasses(J)%Nele=superclasses(J)%Nele+1
               superclasses(J)%classes(superclasses(J)%Nele)=I
           ENDIF
        ENDDO

        ! write(6,*) "Class generation time",t2-t1
        ! WRITE(6,*) "Nclasses, max size, total-non-1 ",Nclasses,MAXVAL(classes%Nele),COUNT(classes%Nele>1)
        ! write(6,*) "isuperclass ",isuperclass

        ! assign the order array to their initial value
        DO Iclass=1,Nclasses
           order(classes(Iclass)%unordered)=classes(Iclass)%reference(classes(Iclass)%order)
        ENDDO

        ! reorder the atoms superclass after superclass 
        itries=0 
        DO I=1,isuperclass
          DO
            itries=itries+1

            ! assign the current order
            DO iclass=1,superclasses(I)%Nele
               J=superclasses(I)%classes(iclass)
               order(classes(J)%unordered)=classes(J)%reference(classes(J)%order)
            ENDDO

            ! check for matches within this superclass only, be happy if we have a match
            matches =  matrix_superclass_equal(reference,unordered,order,superclasses(I),classes)
            IF (itries>max_tries) THEN
               write(6,*) "Could not find the 1-to-1 mapping to prove graph isomorphism"
               write(6,*) "Reordering failed, assuming these molecules are different"
               EXIT
            ENDIF
            IF (matches) EXIT

            ! generate next permutation within this superclass
            DO iclass=1,superclasses(I)%Nele
               J=superclasses(I)%classes(iclass)
               CALL all_permutations(classes(J)%order,classes(J)%Nele, &
                                     classes(J)%q,classes(J)%first)
               IF (.NOT. classes(J)%first) EXIT
            ENDDO

            ! we are back at the original permutation so we're unable to match this superclass.
            IF (iclass.EQ.superclasses(I)%Nele .AND. &
                classes(superclasses(I)%classes(superclasses(I)%Nele))%first) EXIT
          ENDDO
          ! failed in this superblock, can exit now
          IF (.NOT. matches) EXIT
        ENDDO

        ! the final check, just to be sure
        matches=matrix_equal(reference,unordered,order)

        DO Iclass=1,Nclasses
           DEALLOCATE(classes(Iclass)%reference)
           DEALLOCATE(classes(Iclass)%unordered)
           DEALLOCATE(classes(Iclass)%order)
           DEALLOCATE(classes(Iclass)%q)
        ENDDO
        DEALLOCATE(classes)

    END SUBROUTINE reorder_graph
!!****f* graphcon/spread_superclass *
!!
!!   NAME
!!     spread_superclass
!!
!!   FUNCTION
!!     spreads the superclass over all atoms of this class and all their bonded atoms
!!     provided that the latter belong to a class which contains more than one element
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     09.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
    RECURSIVE SUBROUTINE spread_superclass(I,isuperclass,superclass_of_atom,class_of_atom, &
                                             classes, reference)
        INTEGER, INTENT(IN)  :: i,isuperclass
        INTEGER, DIMENSION(:), INTENT(INOUT)   :: superclass_of_atom
        INTEGER, DIMENSION(:), INTENT(IN)      :: class_of_atom
        TYPE(class), DIMENSION(:), INTENT(IN)  :: classes
        TYPE(vertex), DIMENSION(:),INTENT(IN)  :: reference
        INTEGER :: J
        IF (superclass_of_atom(I).EQ.-1 .AND. classes(class_of_atom(I))%Nele>1) THEN
           superclass_of_atom(I)=isuperclass
           DO J=1, classes(class_of_atom(I))%Nele
              CALL spread_superclass(classes(class_of_atom(I))%reference(J), isuperclass, &
                                     superclass_of_atom,class_of_atom,classes,reference)
           ENDDO
           DO J=1,SIZE(reference(I)%bonds)
              CALL spread_superclass(reference(I)%bonds(J),isuperclass, &
                                     superclass_of_atom,class_of_atom,classes,reference)
           ENDDO
        ENDIF
    END SUBROUTINE spread_superclass
!!****f* graphcon/matrix_superclass_equal *
!!
!!   NAME
!!     matrix_superclass_equal
!!
!!   FUNCTION
!!     determines of the vertices of this superclass have the same edges
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     09.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
    FUNCTION matrix_superclass_equal(reference,unordered,order,super,classes) RESULT(res)
        TYPE(vertex), DIMENSION(:),INTENT(IN)  :: reference
        TYPE(vertex), DIMENSION(:),INTENT(IN)  :: unordered
        TYPE(superclass), INTENT(IN)           :: super
        TYPE(class), DIMENSION(:), INTENT(IN)  :: classes
        INTEGER, DIMENSION(:), INTENT(IN) :: order
        LOGICAL                  :: res
        INTEGER :: I,iclass,iele,J
        ! I is the atom in the unordered set
loop:   DO iclass=1,super%Nele
         DO iele=1,classes(super%classes(iclass))%Nele
           I=classes(super%classes(iclass))%unordered(iele)
           res=( reference(order(I))%kind == unordered(I)%kind .AND. &
               SIZE(reference(order(I))%bonds) == SIZE(unordered(I)%bonds) )
           IF (res) THEN
              DO J=1,SIZE(reference(order(I))%bonds)
                IF (ALL(reference(order(I))%bonds(:).NE.order(unordered(I)%bonds(J)))) THEN
                    res=.FALSE.
                    EXIT loop
                ENDIF
              ENDDO
           ELSE
              EXIT loop
           ENDIF
         ENDDO 
        ENDDO loop
    END FUNCTION matrix_superclass_equal
!!****f* graphcon/matrix_equal *
!!
!!   NAME
!!     matrix_equal
!!
!!   FUNCTION
!!     determines of the vertices of the full set is equal, i.e.
!!     we have the same connectivity graph
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     09.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
    FUNCTION matrix_equal(reference,unordered,order) RESULT(res)
        TYPE(vertex), DIMENSION(:),INTENT(IN)  :: reference
        TYPE(vertex), DIMENSION(:),INTENT(IN)  :: unordered
        INTEGER, DIMENSION(:), INTENT(IN) :: order
        LOGICAL                  :: res
        INTEGER :: I,J
loop:   DO I=1,SIZE(reference)
           res=( reference(order(I))%kind == unordered(I)%kind .AND. &
               SIZE(reference(order(I))%bonds) == SIZE(unordered(I)%bonds) )
           IF (res) THEN
              DO J=1,SIZE(reference(order(I))%bonds)
                IF (ALL(reference(order(I))%bonds(:).NE.order(unordered(I)%bonds(J)))) THEN
                    res=.FALSE.
                    EXIT loop
                ENDIF
              ENDDO
           ELSE
              EXIT loop
           ENDIF
        ENDDO loop
    END FUNCTION matrix_equal
!!****f* graphcon/hash_kind *
!!
!!   NAME
!!      hash_kind
!!
!!   FUNCTION
!!      creates a hash for an atom based on its own kind and on the kinds
!!      of its bonded neighbors
!!
!!   NOTES
!!      bonds are sorted so that the order of neighbors appearing in the bonded list
!!      is not important
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     09.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
    FUNCTION hash_kind(me,bonds) Result(res)
       INTEGER, INTENT(IN) :: me
       INTEGER, DIMENSION(:), INTENT(IN) :: bonds
       INTEGER             :: res
       CHARACTER(LEN=250) :: line
       INTEGER :: N,I
       INTEGER, DIMENSION(:), ALLOCATABLE :: ordered_bonds, index
       N=SIZE(bonds)
       ALLOCATE(ordered_bonds(N+1),index(N))
       DO I=1,N
          ordered_bonds(I)=bonds(I)
       ENDDO
       ordered_bonds(N+1)=me
       ! N: only sort the bonds, not me
       CALL sort(ordered_bonds,N,index)
       res=joaat_hash_i(ordered_bonds)
    END FUNCTION hash_kind

!!****f* string_table/joaat_hash_i *
!!
!!   NAME
!!     joaat_hash_i
!!
!!   FUNCTION
!!     generates the hash of an array of integers and the index in the table
!!
!!   NOTES
!!      http://en.wikipedia.org/wiki/Hash_table
!!      http://www.burtleburtle.net/bob/hash/doobs.html
!!      However, since fortran doesn't have an unsigned 4 byte int
!!      we compute it using an integer with the appropriate range
!!      we return already the index in the table as a final result
!!
!!   INPUTS
!!      key : an integer array of any length
!!
!!   MODIFICATION HISTORY
!!      09.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  FUNCTION joaat_hash_i(key) RESULT(hash_index)
    INTEGER, DIMENSION(:), INTENT(IN)        :: key
    INTEGER                                  :: hash_index

    INTEGER, PARAMETER                       :: k64 = SELECTED_INT_KIND(10)
    INTEGER(KIND=k64), PARAMETER             :: b32 = 2_k64**32-1_k64

    INTEGER                                  :: i
    INTEGER(KIND=k64)                        :: hash

    hash=0_k64
    DO i=1,SIZE(key)
       hash=IAND(hash+IBITS(key(i),0,8)                 ,b32)
       hash=IAND(     hash+IAND(ISHFT(hash,10),b32)     ,b32)
       hash=IAND(IEOR(hash,IAND(ISHFT(hash,-6),b32))    ,b32)
       hash=IAND(hash+IBITS(key(i),8,8)                 ,b32)
       hash=IAND(     hash+IAND(ISHFT(hash,10),b32)     ,b32)
       hash=IAND(IEOR(hash,IAND(ISHFT(hash,-6),b32))    ,b32)
       hash=IAND(hash+IBITS(key(i),16,8)                ,b32)
       hash=IAND(     hash+IAND(ISHFT(hash,10),b32)     ,b32)
       hash=IAND(IEOR(hash,IAND(ISHFT(hash,-6),b32))    ,b32)
       hash=IAND(hash+IBITS(key(i),24,8)                ,b32)
       hash=IAND(     hash+IAND(ISHFT(hash,10),b32)     ,b32)
       hash=IAND(IEOR(hash,IAND(ISHFT(hash,-6),b32))    ,b32)
    ENDDO
    hash=IAND(     hash+IAND(ISHFT(hash,  3),b32)  ,b32)
    hash=IAND(IEOR(hash,IAND(ISHFT(hash,-11),b32)) ,b32)
    hash=IAND(     hash+IAND(ISHFT(hash, 15),b32)  ,b32)
    ! hash is the real 32bit hash value of the string,
    ! hash_index is an index in the hash_table
    hash_index=MOD(hash,INT(HUGE(hash_index),KIND=k64))
  END FUNCTION joaat_hash_i

!===ACM Algorithm 323, Generation of Permutations in Lexicographic
!   Order (G6) by R. J. Ord-Smith, CACM 11 (Feb. 1968):117
!   Original Algorithm modified via Certification by I.M. Leitch,
!   17 March 1969.
! Algol to Fortran 77 by H.D.Knoble <hdkLESS at SPAM psu dot edu>,
!                                          May 1995.
!   x = initial values (/1...n/), first=.TRUE.
!   q = scratch
!   first = .TRUE. if you're back at the original order
        SUBROUTINE all_permutations(x,n,q,first)
        Integer n,k,m,t
        Integer x(n),q(n)
        Logical first
        if (n==1) return
        if (first) then
          first=.false.
          do m=1,n-1
            q(m)=n
          end do
        endif
        if(q(n-1).eq.n) then
          q(n-1)=n-1
          t=x(n)
          x(n)=x(n-1)
          x(n-1)=t
          return
        endif
        do k=n-1,1,-1
          if(q(k).eq.k) then
            q(k)=n
          else
            go to 1
          endif
        end do
        first=.true.
        k=1
        goto 2
1       m=q(k)
        t=x(m)
        x(m)=x(k)
        x(k)=t
        q(k)=m-1
        k=k+1
2       m=n
3       t=x(m)
        x(m)=x(k)
        x(k)=t
        m=m-1
        k=k+1
        if(k.lt.m) goto 3
        END SUBROUTINE
END MODULE graphcon

