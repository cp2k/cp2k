!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for a Kim-Gordon-like partitioning into molecular subunits
!>        unsing a vertex coloring algorithm (DSATUR) to find non-interating 
!>        subsets, such that two molecules within the same subset have
!>        small/zero overlap (in other words: this molecular pair is not present in
!>        the neighborlist sab_orb for the current value of EPS_DEFAULT)
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! *****************************************************************************
MODULE kg_vertex_coloring_methods
  USE bibliography,                    ONLY: Brelaz1979,&
                                             cite_reference
  USE f77_blas
  USE kg_environment_types,            ONLY: kg_environment_type
  USE kinds,                           ONLY: int_4
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  TYPE vertex
    INTEGER :: id
    INTEGER :: color
    INTEGER :: degree ! degree (number of neighbors)
    INTEGER :: dsat   ! degree of saturation
    TYPE (vertex), POINTER :: next => NULL() 
    TYPE (vertex_p_type), DIMENSION(:), POINTER :: neighbors => NULL()
  END TYPE vertex
 
  TYPE vertex_p_type
    TYPE (vertex), POINTER :: vertex
  END TYPE vertex_p_type

  TYPE node_type
    INTEGER                                   :: color, nstart
    TYPE (vertex), POINTER                    :: vertex => NULL()
  END TYPE node_type

  TYPE node_p_type
    TYPE(node_type), POINTER :: node
  END TYPE node_p_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_vertex_coloring_methods'

  PUBLIC :: kg_create_graph, kg_vertex_coloring

CONTAINS

  SUBROUTINE kg_create_graph(kg_env, pairs, graph)
    TYPE(kg_environment_type), POINTER       :: kg_env
    INTEGER(KIND=int_4), ALLOCATABLE, &
      DIMENSION(:, :), INTENT(IN)            :: pairs
    TYPE(vertex_p_type), DIMENSION(:), &
      POINTER                                :: graph

    INTEGER                                  :: i, imol, ineighbor, jmol, &
                                                kmol, maxdegree, nneighbors, &
                                                nnodes

! The array pairs contains all interacting (overlapping) pairs of molecules.
! It is assumed to be ordered in the following way:
! (1,2), (1,3), (1,4), ..., (2,3), (2,4), ...
! There are no entries (i,i)
! get the number of nodes = total number of molecules

    nnodes=SIZE(kg_env%molecule_set)

    NULLIFY(graph)
    ALLOCATE(graph(nnodes))

    ! allocate and initialize all vertices
    DO i=1,nnodes
      ALLOCATE(graph(i)%vertex)
      graph(i)%vertex%id = i   ! id = imol (molecular index)
      graph(i)%vertex%color=0  ! means vertex is not colored yet
      graph(i)%vertex%dsat=0   ! no colored neighbors yet
    END DO

    ! allocate the neighbor lists
    imol = 0

    maxdegree=0

    DO i=1,SIZE(pairs,2)
      jmol = pairs(1,i)
      ! counting loop
      IF (jmol.NE.imol) THEN
        IF (imol.NE.0) THEN
          ALLOCATE(graph(imol)%vertex%neighbors(nneighbors))
          graph(imol)%vertex%degree = nneighbors
          IF (nneighbors.GT.maxdegree) maxdegree=nneighbors
        END IF
        imol=jmol
        nneighbors=0
      END IF
      nneighbors=nneighbors+1
    END DO

    IF (imol.NE.0) THEN
      ALLOCATE(graph(imol)%vertex%neighbors(nneighbors))
      graph(imol)%vertex%degree = nneighbors
      IF (nneighbors.GT.maxdegree) maxdegree=nneighbors
    END IF

    kg_env%maxdegree=maxdegree
    
    ! there can be now some nodes that have no neighbors, thus vertex%neighbors
    ! is NOT ASSOCIATED

    ! now add the neighbors - if there are any
    imol=0
    ineighbor=0

    DO i=1,SIZE(pairs,2)
      jmol = pairs(1,i)
      IF (jmol.NE.imol) THEN
        ineighbor=0
        imol = jmol
      END IF
      ineighbor=ineighbor+1
      kmol = pairs(2,i)
      graph(imol)%vertex%neighbors(ineighbor)%vertex => graph(kmol)%vertex
    END DO

  END SUBROUTINE

  SUBROUTINE color_graph(graph, maxdegree, ncolors)
    TYPE(vertex_p_type), DIMENSION(:), &
      POINTER                                :: graph
    INTEGER, INTENT(in)                      :: maxdegree
    INTEGER, INTENT(out)                     :: ncolors

    INTEGER                                  :: color, i, j, newcolor, &
                                                nneighbors, nnodes
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: color_present

    ncolors=0

    nnodes=SIZE(graph)

    ALLOCATE(color_present(nnodes))

    DO i=1,nnodes
      color_present(:)=.FALSE.
      IF (ASSOCIATED(graph(i)%vertex%neighbors)) THEN
        nneighbors=SIZE(graph(i)%vertex%neighbors)
        DO j=1,nneighbors
          color = graph(i)%vertex%neighbors(j)%vertex%color
          IF (color.NE.0) color_present(color)=.TRUE.
        END DO
      END IF
      DO j=1,maxdegree+1 !nnodes
        IF (color_present(j).EQV..FALSE.) THEN
          newcolor=j
          EXIT
        END IF
      END DO
      IF (newcolor.GT.ncolors) ncolors=newcolor
      graph(i)%vertex%color = newcolor ! smallest possible
    END DO

    DEALLOCATE(color_present)

  END SUBROUTINE

  ! returns the number of items in the todolist 'list'
  INTEGER FUNCTION nitem(list)
    TYPE(vertex), POINTER                    :: list, this
    nitem=0
    this => list
    DO WHILE (ASSOCIATED(this%next))
      this=>this%next
      nitem=nitem+1
    END DO
  END FUNCTION

  ! prints the subset info to the screen - useful for vmd visualization
  ! note that the index starts with '0' and not with '1'
  SUBROUTINE print_subsets(graph, ncolors, unit_nr, error)
    TYPE(vertex_p_type), DIMENSION(:), &
      POINTER                                :: graph
    INTEGER, INTENT(IN)                      :: ncolors, unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: counter, i, j, nnodes

    IF (unit_nr>0) THEN

      WRITE(unit_nr,'(T2,A,T10,A)'), "Color |", "Molecules in this subset (IDs start from 0)"

      nnodes = SIZE(graph)
 
      DO i=1,ncolors
        WRITE(unit_nr,'(T2,I5,1X,A)', ADVANCE='NO'), i,"|"
        counter=0
        DO j=1, nnodes
          IF (graph(j)%vertex%color.EQ.i) THEN
            counter=counter+1
            IF (MOD(counter, 13).EQ.0) THEN
              counter=counter+1
              WRITE(unit_nr,'()') ! line break
              WRITE(unit_nr,'(6X,A)', ADVANCE='NO') " |" ! indent next line
            END IF
            WRITE(unit_nr, '(I5,1X)', ADVANCE='NO') graph(j)%vertex%id-1
          END IF
        END DO
        WRITE(unit_nr,'()')
      END DO

    END IF 

  END SUBROUTINE

 ! allocate and initialize the todolist
  SUBROUTINE init_todo(graph, list)
    TYPE(vertex_p_type), DIMENSION(:), &
      POINTER                                :: graph
    TYPE(vertex), POINTER                    :: list

    INTEGER                                  :: i, maxdeg, nnodes
    TYPE(vertex), POINTER                    :: maxnode, this

    nnodes = SIZE(graph)

    ALLOCATE(list)

    ! find a node of maximum degree
    maxdeg = 0

    NULLIFY(maxnode)

    DO i=1, nnodes
      IF (graph(i)%vertex%degree.GT.maxdeg) THEN
        maxdeg = graph(i)%vertex%degree
        maxnode => graph(i)%vertex
      END IF
    END DO

    IF (ASSOCIATED(maxnode)) THEN

      ! first node in the todolist has max degree
      list%next => maxnode

      this => maxnode
 
      ! add all OTHER nodes in arbitrary order...
      DO i=1, nnodes
        IF (graph(i)%vertex%id.NE.maxnode%id) THEN
          this%next => graph(i)%vertex
          this => graph(i)%vertex
        END IF
      END DO

    ELSE
      ! no node with degree > 0

      this => list

      ! add all nodes in arbitrary order
      DO i=1, nnodes
        this%next => graph(i)%vertex
        this => graph(i)%vertex
      END DO

    END IF

  END SUBROUTINE

  ! calculates the dsat of the given node
  SUBROUTINE calculate_dsat(graph, node, dsat) 
    TYPE(vertex_p_type), DIMENSION(:), &
      POINTER                                :: graph
    TYPE(vertex), POINTER                    :: node
    INTEGER, INTENT(out)                     :: dsat

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_dsat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: color, i, nneighbors, nnodes
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: color_present

    nnodes = SIZE(graph)
    nneighbors = SIZE(node%neighbors)
 
    ALLOCATE(color_present(nnodes))
 
    color_present(:)=.FALSE.

    DO i=1,nneighbors
      color = node%neighbors(i)%vertex%color
      IF (color.NE.0) color_present(color) = .TRUE.
    END DO

    dsat=0

    DO i=1,SIZE(color_present)
      IF (color_present(i)) dsat=dsat+1
    END DO

    DEALLOCATE(color_present)

  END SUBROUTINE

  ! removes a given node from the list,
  ! if present - does nothing otherwise
  SUBROUTINE remove(node, list, stat)
    TYPE(vertex), POINTER                    :: node, list
    LOGICAL, INTENT(out), OPTIONAL           :: stat

    TYPE(vertex), POINTER                    :: this

    stat = .FALSE.

    this => list

    DO WHILE (ASSOCIATED(this%next))

      IF (this%next%id.EQ.node%id) THEN
        stat=.TRUE.
        IF (ASSOCIATED(this%next%next)) THEN
          this%next => this%next%next
        ELSE
          this%next => NULL()
        END IF
        RETURN
      END IF 

      this => this%next

    END DO

  END SUBROUTINE

  ! adds a given node to the list, if not already present
  ! does nothing otherwise
  SUBROUTINE add(list, node, stat)
    TYPE(vertex), POINTER                    :: list, node
    LOGICAL, INTENT(out), OPTIONAL           :: stat

    TYPE(vertex), POINTER                    :: this

    stat = .TRUE.
    this => list

    DO WHILE (ASSOCIATED(this%next))

      IF (this%next%id.EQ.node%id) THEN
        stat=.FALSE. ! node already in list
        RETURN
      END IF 

      this => this%next

    END DO

    ! add the node at the top of the list
    node%next => list%next
    list%next => node

  END SUBROUTINE

  ! update position of a given node in the todolist
  SUBROUTINE update_todo(graph, list, node)
    TYPE(vertex_p_type), DIMENSION(:), &
      POINTER                                :: graph
    TYPE(vertex), POINTER                    :: list, node

    CHARACTER(len=*), PARAMETER :: routineN = 'update_todo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: degree, dsat
    LOGICAL                                  :: stat
    TYPE(vertex), POINTER                    :: this

    degree = node%degree
    dsat = node%dsat

    ! remove node from the list (if already present)
    CALL remove(node, list, stat)

    IF (stat) THEN 
      this => list
  
      DO WHILE (ASSOCIATED(this%next))
  
        IF (this%next%dsat.LE.dsat) THEN
          IF (this%next%dsat.EQ.dsat) THEN
            IF (this%next%degree.LE.degree) EXIT
          ELSE
            EXIT
          END IF
        END IF 
  
        this => this%next
  
      END DO
  
      ! insert node after this
      node%next => this%next
      this%next => node
    END IF

  END SUBROUTINE

  SUBROUTINE kg_dsatur(kg_env, graph, ncolors, error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(vertex_p_type), DIMENSION(:), &
      POINTER                                :: graph
    INTEGER(KIND=int_4), INTENT(OUT)         :: ncolors
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_dsatur', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: color, dsat, handle, i, &
                                                newcolor, nneighbors, nnodes
    LOGICAL                                  :: stat
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: color_present
    TYPE(vertex), POINTER                    :: list, neighbor, this

    CALL timeset(routineN,handle)

    CALL cite_reference(Brelaz1979)

    ncolors=0

    CALL init_todo(graph, list)

    nnodes=SIZE(graph)

    IF (kg_env%maxdegree.EQ.0) THEN
      ! all nodes are disconnected

      ncolors=1

      DO i=1,nnodes
        ! only one color needed to color the graph
        graph(i)%vertex%color=1
      END DO
      
    ELSE
  
      ALLOCATE(color_present(nnodes))
  
      DO WHILE (ASSOCIATED(list%next))
  
        this => list%next
  
        color_present(:)=.FALSE.
  
        IF (ASSOCIATED(this%neighbors)) THEN
          nneighbors=SIZE(this%neighbors)
          DO i=1,nneighbors
            color = this%neighbors(i)%vertex%color
            IF (color.NE.0) color_present(color)=.TRUE.
          END DO
        END IF 
        DO i=1,kg_env%maxdegree+1 !nnodes
          IF (color_present(i).EQV..FALSE.) THEN
            newcolor=i
            EXIT
          END IF
        END DO
        this%color = newcolor ! smallest possible
        IF (newcolor.GT.ncolors) ncolors=newcolor
  
        ! remove current node from the todolist
        CALL remove(this, list, stat)
  
        ! update all neighboring nodes
        DO i=1,nneighbors
  
          neighbor => this%neighbors(i)%vertex
  
          CALL calculate_dsat(graph, neighbor, dsat)
          IF (neighbor%dsat.NE.dsat) THEN
  
            ! update dsat
            neighbor%dsat = dsat
  
            ! update position of this neighbor in the todolist
            CALL update_todo(graph, list, neighbor)
  
          END IF
        END DO
      END DO
  
      DEALLOCATE(color_present)   
      DEALLOCATE(list)

    END IF

    CALL timestop(handle)

  END SUBROUTINE
 
  SUBROUTINE deallocate_graph(graph)
    TYPE(vertex_p_type), DIMENSION(:), &
      POINTER                                :: graph

    INTEGER                                  :: i, nnodes

    nnodes = SIZE(graph)

    DO i=1,nnodes
      IF (ASSOCIATED(graph(i)%vertex%neighbors)) DEALLOCATE(graph(i)%vertex%neighbors)
      DEALLOCATE(graph(i)%vertex)
    END DO
    DEALLOCATE(graph)

  END SUBROUTINE


  SUBROUTINE kg_vertex_coloring(kg_env, pairs, ncolors, color_of_node, error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    INTEGER(KIND=int_4), ALLOCATABLE, &
      DIMENSION(:, :), INTENT(IN)            :: pairs
    INTEGER(KIND=int_4), INTENT(OUT)         :: ncolors
    INTEGER(KIND=int_4), ALLOCATABLE, &
      DIMENSION(:), INTENT(out)              :: color_of_node
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_vertex_coloring', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: color, handle, i, nnodes, &
                                                unit_nr
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(vertex_p_type), DIMENSION(:), &
      POINTER                                :: graph

! get a useful output_unit

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    CALL timeset(routineN, handle)

    CALL kg_create_graph(kg_env, pairs, graph)

    ! CALL color_graph(graph, kg_env%maxdegree, ncolors) ! simple greedy algorithm
    CALL kg_dsatur(kg_env, graph, ncolors, error)

    nnodes = SIZE(kg_env%molecule_set)

    ALLOCATE(color_of_node(nnodes))

    ! gather the subset info in a simple integer array
    DO i=1, nnodes
      color = graph(i)%vertex%color
      color_of_node(i) = color
    END DO

    IF (unit_nr>0) THEN

      WRITE(unit_nr,'()')
      WRITE(unit_nr,'(T2,A,A,A)') REPEAT("-",30)," KG coloring result ",REPEAT("-",29)
      !WRITE(unit_nr,'()')
      !CALL print_subsets(graph, ncolors, unit_nr, error)
      !WRITE(unit_nr,'()')
      WRITE(unit_nr,'(T2, A16,50X,I6,1X,A6)') 'Number of colors', ncolors, 'colors'

      WRITE(unit_nr,'(T2,A)') REPEAT("-",79)
      WRITE(unit_nr,'()')

    END IF

    CALL deallocate_graph(graph)

    CALL timestop(handle)

  END SUBROUTINE

END MODULE
