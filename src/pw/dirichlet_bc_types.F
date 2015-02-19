!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Dirichlet boundary condition data types
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE dirichlet_bc_types

  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: unit_matrix,&
                                             vector_product
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: pw_release,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dirichlet_bc_types'

! main types
  TYPE :: dirichlet_bc_parameters
     LOGICAL                                  :: debug

     INTEGER                                  :: n_aa_planar
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: aa_planar_xxtnt => NULL()
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: aa_planar_yxtnt => NULL()
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: aa_planar_zxtnt => NULL()
     INTEGER, DIMENSION(:,:), POINTER         :: aa_planar_nprtn => NULL()
     REAL(KIND=dp), DIMENSION(:), POINTER     :: aa_planar_vD => NULL()
     REAL(KIND=dp), DIMENSION(:), POINTER     :: aa_planar_zeta => NULL()
     LOGICAL, DIMENSION(:), POINTER           :: aa_planar_smooth => NULL()

     INTEGER                                  :: n_planar
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: planar_Avtx  => NULL()
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: planar_Bvtx  => NULL()
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: planar_Cvtx  => NULL()
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: planar_Dvtx  => NULL()
     INTEGER, DIMENSION(:,:), POINTER         :: planar_nprtn => NULL()
     REAL(KIND=dp), DIMENSION(:), POINTER     :: planar_vD => NULL()
     REAL(KIND=dp), DIMENSION(:), POINTER     :: planar_zeta => NULL()
     LOGICAL, DIMENSION(:), POINTER           :: planar_smooth => NULL()

     INTEGER                                  :: n_xaa_cylindrical
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: xaa_cylindrical_xxtnt  => NULL()
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: xaa_cylindrical_bctr   => NULL()
     REAL(KIND=dp), DIMENSION(:), POINTER     :: xaa_cylindrical_brad   => NULL()
     INTEGER, DIMENSION(:), POINTER           :: xaa_cylindrical_nsides => NULL()
     INTEGER, DIMENSION(:), POINTER           :: xaa_cylindrical_apxtyp => NULL()
     INTEGER, DIMENSION(:,:), POINTER         :: xaa_cylindrical_nprtn  => NULL()
     REAL(KIND=dp), DIMENSION(:), POINTER     :: xaa_cylindrical_vD => NULL()
     REAL(KIND=dp), DIMENSION(:), POINTER     :: xaa_cylindrical_zeta => NULL()
     LOGICAL, DIMENSION(:), POINTER           :: xaa_cylindrical_smooth => NULL()

     INTEGER                                  :: n_aa_cuboidal
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: aa_cuboidal_xxtnt => NULL()
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: aa_cuboidal_yxtnt => NULL()
     REAL(KIND=dp), DIMENSION(:,:), POINTER   :: aa_cuboidal_zxtnt => NULL()
     INTEGER, DIMENSION(:,:), POINTER         :: aa_cuboidal_nprtn => NULL()
     REAL(KIND=dp), DIMENSION(:), POINTER     :: aa_cuboidal_vD => NULL()
     REAL(KIND=dp), DIMENSION(:), POINTER     :: aa_cuboidal_zeta => NULL()
     LOGICAL, DIMENSION(:), POINTER           :: aa_cuboidal_smooth => NULL()
  END TYPE dirichlet_bc_parameters

  TYPE :: dirichlet_bc_type
     REAL(dp)                          :: v_D
     INTEGER                           :: dbc_id
     INTEGER                           :: dbc_geom
     TYPE(cs_rectangle), POINTER       :: rectangle => NULL()
     TYPE(cs_box), POINTER             :: box => NULL()
     TYPE(tile_p_type), &
       DIMENSION(:), POINTER           :: tiles => NULL()
     INTEGER                           :: n_tiles
     LOGICAL                           :: smooth
     REAL(dp)                          :: mollifier_zeta
  END TYPE dirichlet_bc_type

  TYPE :: dirichlet_bc_p_type
     TYPE(dirichlet_bc_type), POINTER  :: dirichlet_bc
     INTEGER                           :: tag
  END TYPE dirichlet_bc_p_type

  TYPE :: tile_type
     REAL(dp)                          :: v_D
     INTEGER                           :: tile_id
     INTEGER                           :: tile_geom
     TYPE(cs_rectangle), POINTER       :: rectangle => NULL()
     TYPE(cs_box), POINTER             :: box => NULL()
     TYPE(pw_type), POINTER            :: tile_pw => NULL()
     INTEGER                           :: npts
     LOGICAL                           :: smooth
     REAL(dp)                          :: mollifier_zeta
  END TYPE tile_type

  TYPE :: tile_p_type
     TYPE(tile_type), POINTER          :: tile 
     INTEGER                           :: tag
  END TYPE tile_p_type

! auxiliary types (3D Cartesian space primitive objects)
  TYPE :: cs_segment
     REAL(dp), DIMENSION(3)            :: A
     REAL(dp), DIMENSION(3)            :: B
     REAL(dp), DIMENSION(3)            :: AB
  END TYPE cs_segment

  TYPE :: cs_edge
     TYPE(cs_segment), POINTER         :: segment => NULL()
     INTEGER                           :: tag
  END TYPE cs_edge

  TYPE :: cs_plane
     REAL(dp), DIMENSION(3)            :: A, B, C
     TYPE(cs_segment), POINTER         :: AB => NULL()
     TYPE(cs_segment), POINTER         :: AC => NULL()
     REAL(dp), DIMENSION(3)            :: normal
     REAL(dp), DIMENSION(3)            :: unit_normal
     REAL(dp)                          :: d
  END TYPE cs_plane

  TYPE :: cs_triangle
     REAL(dp), DIMENSION(3)            :: A, B, C
     REAL(dp), DIMENSION(3,3)          :: vertices
     TYPE(cs_edge), DIMENSION(3)       :: edges
     TYPE(cs_plane), POINTER           :: plane => NULL()
  END TYPE cs_triangle

  TYPE :: cs_rectangle
     REAL(dp), DIMENSION(3)            :: A, B, C, D
     REAL(dp), DIMENSION(3,4)          :: vertices
     TYPE(cs_edge), DIMENSION(4)       :: edges
     TYPE(cs_plane), POINTER           :: plane => NULL()
  END TYPE cs_rectangle

  TYPE :: cs_polygon
     REAL(dp), DIMENSION(:,:), POINTER :: vertices => NULL()
     TYPE(cs_edge), &
       DIMENSION(:), POINTER           :: edges => NULL()
     TYPE(cs_plane), POINTER           :: plane => NULL()
     INTEGER                           :: n_vertices
  END TYPE cs_polygon

  TYPE :: cs_face
     TYPE(cs_rectangle), POINTER       :: rectangle => NULL()
     INTEGER                           :: tag
  END TYPE cs_face

  TYPE :: cs_box
     REAL(dp), DIMENSION(3,8)          :: vertices
     TYPE(cs_edge), DIMENSION(12)      :: edges
     TYPE(cs_face), DIMENSION(6)       :: faces
  END TYPE cs_box

!==
  TYPE :: cs_intersection_obj
     REAL(dp), DIMENSION(3)            :: point
     TYPE(cs_segment), POINTER         :: segment => NULL()
     TYPE(cs_triangle), POINTER        :: triangle => NULL()
     TYPE(cs_rectangle), POINTER       :: rectangle => NULL()
     TYPE(cs_polygon), POINTER         :: polygon => NULL()

     INTEGER                           :: n_intxn
     INTEGER                           :: intxn_type
  END TYPE cs_intersection_obj

!==
  TYPE :: cs_grid_segment
     REAL(dp), DIMENSION(3)            :: orig
     INTEGER, DIMENSION(3)             :: grid_segment_idx
     INTEGER, DIMENSION(2)             :: weight
     INTEGER                           :: alignment
     TYPE(cs_segment), POINTER         :: segment => NULL()
     LOGICAL                           :: is_on
  END TYPE cs_grid_segment

  TYPE :: cs_grid_tile
     REAL(dp), DIMENSION(3)            :: orig
     INTEGER, DIMENSION(3)             :: grid_tile_idx
     INTEGER, DIMENSION(4)             :: weight
     INTEGER                           :: alignment
     TYPE(cs_rectangle), POINTER       :: rectangle => NULL()
     LOGICAL                           :: is_on
  END TYPE cs_grid_tile

  TYPE :: cs_voxel
     REAL(dp), DIMENSION(3)            :: orig
     INTEGER, DIMENSION(3)             :: voxel_idx
     INTEGER, DIMENSION(8)             :: weight
     TYPE(cs_box), POINTER             :: box => NULL()
     LOGICAL                           :: is_on
  END TYPE cs_voxel

! publicly accessible entities
  PUBLIC dirichlet_bc_type, dirichlet_bc_p_type, tile_type, tile_p_type
  PUBLIC dirichlet_bc_parameters

  PUBLIC dirichlet_boundary_region_release,&
         dirichlet_bc_parameters_dealloc

  PUBLIC cs_segment, cs_plane, &
         cs_triangle, cs_rectangle, cs_polygon, cs_box, &
         cs_edge, cs_face, &
         cs_grid_segment, cs_grid_tile, cs_voxel, &
         cs_intersection_obj

  PUBLIC cs_segment_create, cs_plane_create, &
         cs_triangle_create, cs_rectangle_create, cs_polygon_create, cs_box_create, &
         cs_edge_create, cs_face_create, &
         cs_grid_segment_create, cs_grid_tile_create, cs_voxel_create, &
         cs_intersection_obj_create

  PUBLIC cs_segment_release, cs_plane_release, &
         cs_triangle_release, cs_rectangle_release, cs_polygon_release, cs_box_release, &
         cs_edge_release, cs_face_release, &
         cs_grid_segment_release, cs_grid_tile_release, cs_voxel_release, &
         cs_intersection_obj_release

  PUBLIC cast_triangle_to_polygon, cast_rectangle_to_polygon, sort_vertices

! interface declaration
  INTERFACE cs_plane_create
     MODULE PROCEDURE cs_plane_ABC_create, &
                      cs_plane_An_create
  END INTERFACE cs_plane_create

  INTERFACE cs_rectangle_create
     MODULE PROCEDURE cs_rectangle_ABCD_create, &
                      cs_rectangle_ABC_create
  END INTERFACE cs_rectangle_create

! parameters
  INTEGER, PARAMETER, PUBLIC           :: AA_PLANAR   = 100, &
                                          PLANAR      = 200, &
                                          CYLINDRICAL = 300, &
                                          AA_CUBOIDAL = 400

  INTEGER, PARAMETER, PUBLIC           :: rectangular_tile = 1000, &
                                          cuboidal_tile    = 2000

  INTEGER, PARAMETER, PUBLIC           :: CIRCUMSCRIBED = 10, &
                                          INSCRIBED     = 20

  INTEGER, PARAMETER, PUBLIC           :: NO_INTXN   = 0, &
                                          PT_INTXN   = 1, &
                                          SEG_INTXN  = 2, &
                                          TRI_INTXN  = 3, &
                                          RECT_INTXN = 4, &
                                          PY_INTXN   = 5

  INTEGER, PARAMETER, PUBLIC           :: no_grid_point         = 100, &
                                          single_grid_point     = 110, &
                                          xaxis_aligned_segment = 121, &
                                          yaxis_aligned_segment = 122, &
                                          zaxis_aligned_segment = 123, &
                                          xy_aligned_rectangle  = 141, &
                                          xz_aligned_rectangle  = 142, &
                                          yz_aligned_rectangle  = 143, &
                                          axis_aligned_box      = 180

  INTEGER, PARAMETER, PUBLIC           :: FWROT =  1, &
                                          BWROT = -1
! private parameters
  REAL(dp), PARAMETER, PRIVATE         :: small_value = 1.0E-8_dp
  REAL(dp), PARAMETER, PRIVATE         :: huge_value = HUGE(0.0_dp)

CONTAINS

! *****************************************************************************
!> \brief   releases the defined Dirichlet boundary region 
!> \param gates the DBC region to be released  
!> \param pw_pool pool of the plane wave grid
!> \param error cp2k error
!>
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dirichlet_boundary_region_release(gates, pw_pool, error)

    TYPE(dirichlet_bc_p_type), &
      DIMENSION(:), INTENT(INOUT), POINTER   :: gates
    TYPE(pw_pool_type), INTENT(IN), &
      OPTIONAL, POINTER                      :: pw_pool
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dirichlet_boundary_region_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dbc_geom, handle, i, k, &
                                                n_gates, n_tiles, tile_geom

    CALL timeset(routineN,handle)

    n_gates = SIZE(gates)

    IF (ASSOCIATED(gates)) THEN
       IF (PRESENT(pw_pool)) THEN
          DO i = 1, n_gates
             dbc_geom = gates(i)%dirichlet_bc%dbc_geom
             IF (dbc_geom .EQ. AA_CUBOIDAL) THEN
                CALL cs_box_release(gates(i)%dirichlet_bc%box, error=error)
             ELSE
                CALL cs_rectangle_release(gates(i)%dirichlet_bc%rectangle, error=error)
             END IF

             n_tiles = gates(i)%dirichlet_bc%n_tiles
             DO k = 1, n_tiles
                tile_geom = gates(i)%dirichlet_bc%tiles(k)%tile%tile_geom
                CALL pw_pool_give_back_pw(pw_pool, gates(i)%dirichlet_bc%tiles(k)%tile%tile_pw, error=error)

                IF (tile_geom .EQ. rectangular_tile) THEN
                   CALL cs_rectangle_release(gates(i)%dirichlet_bc%tiles(k)%tile%rectangle, error)
                ELSE IF (tile_geom .EQ. cuboidal_tile) THEN
                   CALL cs_box_release(gates(i)%dirichlet_bc%tiles(k)%tile%box, error)
                END IF

                DEALLOCATE(gates(i)%dirichlet_bc%tiles(k)%tile)
             END DO
             DEALLOCATE(gates(i)%dirichlet_bc%tiles)

             DEALLOCATE(gates(i)%dirichlet_bc)
          END DO
       ELSE 
          DO i = 1, n_gates
             dbc_geom = gates(i)%dirichlet_bc%dbc_geom
             IF (dbc_geom .EQ. AA_CUBOIDAL) THEN
                CALL cs_box_release(gates(i)%dirichlet_bc%box, error=error)
             ELSE
                CALL cs_rectangle_release(gates(i)%dirichlet_bc%rectangle, error=error)
             END IF

             n_tiles = gates(i)%dirichlet_bc%n_tiles
             DO k = 1, n_tiles
                tile_geom = gates(i)%dirichlet_bc%tiles(k)%tile%tile_geom
                CALL pw_release(gates(i)%dirichlet_bc%tiles(k)%tile%tile_pw, error=error)

                IF (tile_geom .EQ. rectangular_tile) THEN
                   CALL cs_rectangle_release(gates(i)%dirichlet_bc%tiles(k)%tile%rectangle, error)
                ELSE IF (tile_geom .EQ. cuboidal_tile) THEN
                   CALL cs_box_release(gates(i)%dirichlet_bc%tiles(k)%tile%box, error)
                END IF

                DEALLOCATE(gates(i)%dirichlet_bc%tiles(k)%tile)
             END DO
             DEALLOCATE(gates(i)%dirichlet_bc%tiles)

             DEALLOCATE(gates(i)%dirichlet_bc)
          END DO
       END IF
       DEALLOCATE(gates)
    END IF

    CALL timestop(handle)

  END SUBROUTINE dirichlet_boundary_region_release

! *****************************************************************************
!> \brief   deallocates dirichlet_bc_parameters type
!> \param dbc_params dbc parameters
!> \param error cp2k error 
!>
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dirichlet_bc_parameters_dealloc(dbc_params, error)

    TYPE(dirichlet_bc_parameters), &
      INTENT(INOUT)                          :: dbc_params
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dirichlet_bc_parameters_dealloc', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(dbc_params%aa_planar_xxtnt))  DEALLOCATE(dbc_params%aa_planar_xxtnt)
    IF (ASSOCIATED(dbc_params%aa_planar_yxtnt))  DEALLOCATE(dbc_params%aa_planar_yxtnt)
    IF (ASSOCIATED(dbc_params%aa_planar_zxtnt))  DEALLOCATE(dbc_params%aa_planar_zxtnt)
    IF (ASSOCIATED(dbc_params%aa_planar_nprtn))  DEALLOCATE(dbc_params%aa_planar_nprtn)
    IF (ASSOCIATED(dbc_params%aa_planar_vD))     DEALLOCATE(dbc_params%aa_planar_vD)
    IF (ASSOCIATED(dbc_params%aa_planar_zeta))   DEALLOCATE(dbc_params%aa_planar_zeta)
    IF (ASSOCIATED(dbc_params%aa_planar_smooth)) DEALLOCATE(dbc_params%aa_planar_smooth)

    IF (ASSOCIATED(dbc_params%xaa_cylindrical_xxtnt))  DEALLOCATE(dbc_params%xaa_cylindrical_xxtnt)
    IF (ASSOCIATED(dbc_params%xaa_cylindrical_bctr))   DEALLOCATE(dbc_params%xaa_cylindrical_bctr)
    IF (ASSOCIATED(dbc_params%xaa_cylindrical_brad))   DEALLOCATE(dbc_params%xaa_cylindrical_brad)
    IF (ASSOCIATED(dbc_params%xaa_cylindrical_apxtyp)) DEALLOCATE(dbc_params%xaa_cylindrical_apxtyp)
    IF (ASSOCIATED(dbc_params%xaa_cylindrical_nsides)) DEALLOCATE(dbc_params%xaa_cylindrical_nsides)
    IF (ASSOCIATED(dbc_params%xaa_cylindrical_nprtn))  DEALLOCATE(dbc_params%xaa_cylindrical_nprtn)
    IF (ASSOCIATED(dbc_params%xaa_cylindrical_vD))     DEALLOCATE(dbc_params%xaa_cylindrical_vD)
    IF (ASSOCIATED(dbc_params%xaa_cylindrical_zeta))   DEALLOCATE(dbc_params%xaa_cylindrical_zeta)
    IF (ASSOCIATED(dbc_params%xaa_cylindrical_smooth)) DEALLOCATE(dbc_params%xaa_cylindrical_smooth)

    IF (ASSOCIATED(dbc_params%planar_Avtx))   DEALLOCATE(dbc_params%planar_Avtx)
    IF (ASSOCIATED(dbc_params%planar_Bvtx))   DEALLOCATE(dbc_params%planar_Bvtx)
    IF (ASSOCIATED(dbc_params%planar_Cvtx))   DEALLOCATE(dbc_params%planar_Cvtx)
    IF (ASSOCIATED(dbc_params%planar_Dvtx))   DEALLOCATE(dbc_params%planar_Dvtx)
    IF (ASSOCIATED(dbc_params%planar_nprtn))  DEALLOCATE(dbc_params%planar_nprtn)
    IF (ASSOCIATED(dbc_params%planar_vD))     DEALLOCATE(dbc_params%planar_vD)
    IF (ASSOCIATED(dbc_params%planar_zeta))   DEALLOCATE(dbc_params%planar_zeta)
    IF (ASSOCIATED(dbc_params%planar_smooth)) DEALLOCATE(dbc_params%planar_smooth)

    IF (ASSOCIATED(dbc_params%aa_cuboidal_xxtnt))  DEALLOCATE(dbc_params%aa_cuboidal_xxtnt)
    IF (ASSOCIATED(dbc_params%aa_cuboidal_yxtnt))  DEALLOCATE(dbc_params%aa_cuboidal_yxtnt)
    IF (ASSOCIATED(dbc_params%aa_cuboidal_zxtnt))  DEALLOCATE(dbc_params%aa_cuboidal_zxtnt)
    IF (ASSOCIATED(dbc_params%aa_cuboidal_nprtn))  DEALLOCATE(dbc_params%aa_cuboidal_nprtn)
    IF (ASSOCIATED(dbc_params%aa_cuboidal_vD))     DEALLOCATE(dbc_params%aa_cuboidal_vD)
    IF (ASSOCIATED(dbc_params%aa_cuboidal_zeta))   DEALLOCATE(dbc_params%aa_cuboidal_zeta)
    IF (ASSOCIATED(dbc_params%aa_cuboidal_smooth)) DEALLOCATE(dbc_params%aa_cuboidal_smooth)

  END SUBROUTINE dirichlet_bc_parameters_dealloc

! *****************************************************************************
!> \brief   generates a segment
!> \param segment the segment to be created
!> \param A coordinates of the end point A
!> \param B coordinates of the end point B
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_segment_create(segment, A, B, error)

    TYPE(cs_segment), INTENT(INOUT), POINTER :: segment
    REAL(dp), DIMENSION(3), INTENT(IN)       :: A, B
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_segment_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    segment%A = A
    segment%B = B
    segment%AB = B - A

    CALL timestop(handle)

  END SUBROUTINE cs_segment_create

! *****************************************************************************
!> \brief   generates an edge joining two vertices of an object in the 3D Cartesian space
!> \param edge the edge to be created
!> \param A coordinates of the vertex A
!> \param B coordinates of the vertex B
!> \param tag unique tag to be assigned to the edge
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_edge_create(edge, A, B, tag, error)

    TYPE(cs_edge), INTENT(INOUT), POINTER    :: edge
    REAL(dp), DIMENSION(3), INTENT(IN)       :: A, B
    INTEGER, INTENT(IN)                      :: tag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_edge_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    ALLOCATE(edge%segment)
    CALL cs_segment_create(edge%segment, A, B, error)
    edge%tag = tag

    CALL timestop(handle)

  END SUBROUTINE cs_edge_create

! *****************************************************************************
!> \brief   generates a plane given the coordinates of three points on it
!> \param plane the plane to be created
!> \param A coordinates of the vertex A
!> \param B coordinates of the vertex B
!> \param C coordinates of the vertex C
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_plane_ABC_create(plane, A, B, C, error)

    TYPE(cs_plane), INTENT(INOUT), POINTER   :: plane
    REAL(dp), DIMENSION(3), INTENT(IN)       :: A, B, C
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_plane_ABC_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: AB_norm, AC_norm, BC_norm
    REAL(dp), DIMENSION(3)                   :: AB, AC, BC, xy_normal, &
                                                xz_normal, yz_normal

    CALL timeset(routineN,handle)

    xy_normal = (/0.0_dp, 0.0_dp, 1.0_dp/)
    xz_normal = (/0.0_dp, 1.0_dp, 0.0_dp/)
    yz_normal = (/1.0_dp, 0.0_dp, 0.0_dp/)

    AB = B - A
    AC = C - A
    BC = C - B

    plane%A = A
    plane%B = B
    plane%C = C

    ALLOCATE(plane%AB, plane%AC)
    CALL cs_segment_create(plane%AB, A, B, error)
    CALL cs_segment_create(plane%AC, A, C, error)

    AB_norm = SQRT(SUM(AB**2))
    AC_norm = SQRT(SUM(AC**2))
    BC_norm = SQRT(SUM(BC**2))

    IF ( (AB_norm .LE. small_value) .OR. (AC_norm.LE. small_value) .OR. (BC_norm .LE. small_value) ) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "Not enough points to define a cartesian space plane",error)
    ELSE IF ( ABS(ABS(SUM(AB/AC) - 3)) .LE. small_value ) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "Unable to define a cartesian space plane. reason: A, B, C are colinear.",error)
    ELSE
       plane%normal = vector_product(AB, AC)
       plane%unit_normal = plane%normal / SQRT(SUM(plane%normal**2))
!       plane%d = -1.0_dp * DOT_PRODUCT(plane%unit_normal, A)
       plane%d = DOT_PRODUCT(plane%unit_normal, A)

    END IF

    CALL timestop(handle)

  END SUBROUTINE cs_plane_ABC_create

! *****************************************************************************
!> \brief   generates a plane using a point and the unit normal vector to the plane
!> \param plane the plane to be created
!> \param A coordinates of a point on the plane
!> \param n the unit normal vector to the plane
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_plane_An_create(plane, A, n, error)

    TYPE(cs_plane), INTENT(INOUT), POINTER   :: plane
    REAL(dp), DIMENSION(3), INTENT(IN)       :: A, n
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_plane_An_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp), DIMENSION(3)                   :: AB, AC, B, BC, C, xy_normal, &
                                                xz_normal, yz_normal
    REAL(dp), DIMENSION(3, 1)                :: n_prime
    REAL(dp), DIMENSION(3, 3)                :: H, I

    CALL timeset(routineN,handle)

! create the mirror vector
    IF (ABS(n(1) - 1) .GT. ABS(n(1) + 1)) THEN
       n_prime(1,1) = n(1) - 1
    ELSE
       n_prime(1,1) = n(1) + 1
    END IF
    n_prime(2,1) = n(2)
    n_prime(3,1) = n(3)
! create the Householder matrix
    CALL unit_matrix(I)
    H = I - 2 * MATMUL(n_prime,TRANSPOSE(n_prime))/SUM(n_prime**2)
! find two more points on the plane
    B = A + H(2,:)
    C = A + H(3,:)

    yz_normal = I(1,:);    xz_normal = I(2,:);    xy_normal = I(3,:)

    AB = B - A
    AC = C - A
    BC = C - B

    plane%A = A
    plane%B = B
    plane%C = C

    ALLOCATE(plane%AB, plane%AC)
    CALL cs_segment_create(plane%AB, A, B, error)
    CALL cs_segment_create(plane%AC, A, C, error)
! NB: the Housholder procedure is quite stable no need to check if the points are colinear
! or any two of them coincide
    plane%normal = n
    plane%unit_normal = plane%normal / SQRT(SUM(plane%normal**2))
!    plane%d = -1.0_dp * DOT_PRODUCT(plane%unit_normal, A)
    plane%d = DOT_PRODUCT(plane%unit_normal, A)

    CALL timestop(handle)

  END SUBROUTINE cs_plane_An_create

! *****************************************************************************
!> \brief   generates a triangle given the coordinates of its three vertices
!> \param triangle the triangle to be created
!> \param A coordinates of the vertex A
!> \param B coordinates of the vertex B
!> \param C coordinates of the vertex C
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_triangle_create(triangle, A, B, C, error)

    TYPE(cs_triangle), INTENT(INOUT), &
      POINTER                                :: triangle
    REAL(dp), DIMENSION(3), INTENT(IN)       :: A, B, C
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_triangle_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, tag
    REAL(dp), DIMENSION(3, 0:2)              :: tri_vertices
    TYPE(cs_edge), POINTER                   :: edge

    CALL timeset(routineN,handle)

    triangle%A = A
    triangle%B = B
    triangle%C = C
    ALLOCATE(triangle%plane)
    CALL cs_plane_create(triangle%plane, A, B, C, error)

    triangle%vertices(:,1) = A
    triangle%vertices(:,2) = B
    triangle%vertices(:,3) = C
    CALL sort_vertices(triangle%vertices, A, triangle%plane%unit_normal, error)
    tri_vertices = triangle%vertices

    ALLOCATE(edge)
    DO i = 0, 2
       tag = (i+1)*10 + MOD(i+1, 3)+1
       CALL cs_edge_create(edge, tri_vertices(:,i), tri_vertices(:,MOD(i+1, 3)), tag, error)
       triangle%edges(i+1) = edge
    END DO
    DEALLOCATE(edge)

    CALL timestop(handle)

  END SUBROUTINE cs_triangle_create

! *****************************************************************************
!> \brief   generates a rectangle given the coordinates of its four vertices
!> \param rectangle the rectangle to be created
!> \param A coordinates of the vertex A
!> \param B coordinates of the vertex B
!> \param C coordinates of the vertex C
!> \param D coordinates of the vertex D
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_rectangle_ABCD_create(rectangle, A, B, C, D, error)

    TYPE(cs_rectangle), INTENT(INOUT), &
      POINTER                                :: rectangle
    REAL(dp), DIMENSION(3), INTENT(IN)       :: A, B, C, D
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_rectangle_ABCD_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, tag
    REAL(dp), DIMENSION(3, 0:3)              :: rect_vertices
    TYPE(cs_edge), POINTER                   :: edge

    CALL timeset(routineN,handle)

    rectangle%A = A
    rectangle%B = B
    rectangle%C = C
    rectangle%D = D
    ALLOCATE(rectangle%plane)
    CALL cs_plane_create(rectangle%plane, A, B, C, error)

    rectangle%vertices(:,1) = A
    rectangle%vertices(:,2) = B
    rectangle%vertices(:,3) = C
    rectangle%vertices(:,4) = D
    CALL sort_vertices(rectangle%vertices, A, rectangle%plane%unit_normal, error)
    rect_vertices = rectangle%vertices

    ALLOCATE(edge)
    DO i = 0, 3
       tag = (i+1)*10 + MOD(i+1, 4)+1
       CALL cs_edge_create(edge, rect_vertices(:,i), rect_vertices(:,MOD(i+1, 4)), tag, error)
       rectangle%edges(i+1) = edge
    END DO
    DEALLOCATE(edge)

    CALL timestop(handle)

  END SUBROUTINE cs_rectangle_ABCD_create

! *****************************************************************************
!> \brief    generates a rectangle given the coordinates of three veritces. The
!>           vertices are assumed to be ordered counter-clockwise wrt the outer
!>           unit normal to the rectangle's plane.
!> \param rectangle the rectangle to be created
!> \param A coordinates of the vertex A
!> \param B coordinates of the vertex B
!> \param C coordinates of the vertex C
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_rectangle_ABC_create(rectangle, A, B, C, error)

    TYPE(cs_rectangle), INTENT(INOUT), &
      POINTER                                :: rectangle
    REAL(dp), DIMENSION(3), INTENT(IN)       :: A, B, C
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_rectangle_ABC_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, tag
    REAL(dp), DIMENSION(3)                   :: D
    REAL(dp), DIMENSION(3, 0:3)              :: rect_vertices
    TYPE(cs_edge), POINTER                   :: edge

    CALL timeset(routineN,handle)

    D = A + (C - B)

    rectangle%A = A
    rectangle%B = B
    rectangle%C = C
    rectangle%D = D
    ALLOCATE(rectangle%plane)
    CALL cs_plane_create(rectangle%plane, A, B, C, error)

    rectangle%vertices(:,1) = A
    rectangle%vertices(:,2) = B
    rectangle%vertices(:,3) = C
    rectangle%vertices(:,4) = D
    CALL sort_vertices(rectangle%vertices, A, rectangle%plane%unit_normal, error)
    rect_vertices = rectangle%vertices

    ALLOCATE(edge)
    DO i = 0, 3
       tag = (i+1)*10 + MOD(i+1, 4)+1
       CALL cs_edge_create(edge, rect_vertices(:,i), rect_vertices(:,MOD(i+1, 4)), tag, error)
       rectangle%edges(i+1) = edge
    END DO
    DEALLOCATE(edge)

    CALL timestop(handle)

  END SUBROUTINE cs_rectangle_ABC_create

! *****************************************************************************
!> \brief   generates a polygon give the coordinates of its vertices
!> \param polygon the polygon to be created
!> \param vertices coordinates of the vertices of the polygon
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_polygon_create(polygon, vertices, error)

    TYPE(cs_polygon), INTENT(INOUT), POINTER :: polygon
    REAL(dp), DIMENSION(:, :), INTENT(IN), &
      POINTER                                :: vertices
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_polygon_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, num_vertices, tag
    REAL(dp) :: py_vertices(3, 0:SIZE(vertices)/3-1)
    TYPE(cs_edge), POINTER                   :: edge

    CALL timeset(routineN,handle)

    num_vertices = SIZE(vertices)/3
    ALLOCATE(polygon%vertices(3,num_vertices))
    polygon%vertices = vertices
    ALLOCATE(polygon%plane)
    CALL cs_plane_create(polygon%plane, vertices(:,1), vertices(:,2), vertices(:,3), error)
    polygon%n_vertices = num_vertices

! assume vertices are already sorted
    py_vertices = polygon%vertices
    ALLOCATE(polygon%edges(num_vertices))
    ALLOCATE(edge)
    DO i = 0, num_vertices-1
       tag = (i+1)*10 + MOD(i+1, num_vertices)+1
       CALL cs_edge_create(edge, py_vertices(:,i), py_vertices(:,MOD(i+1, num_vertices)), tag, error)
       polygon%edges(i+1) = edge
    END DO
    DEALLOCATE(edge)

    CALL timestop(handle)

  END SUBROUTINE cs_polygon_create

! *****************************************************************************
!> \brief   generates a planar face as part of the boundary of a solid 3D object
!> \param face the face to be created
!> \param A coordinates of the vertex A
!> \param B coordinates of the vertex B
!> \param C coordinates of the vertex C
!> \param tag unique tag to be assigned to the face
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_face_create(face, A, B, C, tag, error)

    TYPE(cs_face), INTENT(INOUT), POINTER    :: face
    REAL(dp), DIMENSION(3), INTENT(IN)       :: A, B, C
    INTEGER, INTENT(IN)                      :: tag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_face_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    ALLOCATE(face%rectangle)
    CALL cs_rectangle_create(face%rectangle, A, B, C, error)
    face%tag = tag

    CALL timestop(handle)

  END SUBROUTINE cs_face_create

! *****************************************************************************
!> \brief  generates an axis-aligned box in the 3D Cartesian space with vertices
!>         labeled as:        6--------5
!>                           /|       /|
!>                          7-|------8 |
!>                          | 1------|-4
!>                          |/       |/
!>                          2--------3
!> \param box the box to be generated
!> \param vertices vertices of the box
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_box_create(box, vertices, error)

    TYPE(cs_box), INTENT(INOUT), POINTER     :: box
    REAL(dp), DIMENSION(3, 8), INTENT(IN)    :: vertices
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_box_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cs_edge), POINTER                   :: edge
    TYPE(cs_face), POINTER                   :: face

    CALL timeset(routineN,handle)

    ALLOCATE(edge, face)
    box%vertices = vertices

! faces taged in a way that vertices are counter-clockwise with respect to the outer unit normal to the face
    CALL cs_face_create(face, vertices(:,2), vertices(:,1), vertices(:,4), 1432, error);   box%faces(1) = face
    CALL cs_face_create(face, vertices(:,1), vertices(:,6), vertices(:,5), 1654, error);   box%faces(2) = face
    CALL cs_face_create(face, vertices(:,2), vertices(:,7), vertices(:,6), 2761, error);   box%faces(3) = face
    CALL cs_face_create(face, vertices(:,6), vertices(:,7), vertices(:,8), 6785, error);   box%faces(4) = face
    CALL cs_face_create(face, vertices(:,7), vertices(:,2), vertices(:,3), 7238, error);   box%faces(5) = face
    CALL cs_face_create(face, vertices(:,8), vertices(:,3), vertices(:,4), 8345, error);   box%faces(6) = face

    CALL cs_edge_create(edge, vertices(:,1), vertices(:,2), 12, error);   box%edges(1)  = edge
    CALL cs_edge_create(edge, vertices(:,2), vertices(:,3), 23, error);   box%edges(2)  = edge
    CALL cs_edge_create(edge, vertices(:,3), vertices(:,4), 34, error);   box%edges(3)  = edge
    CALL cs_edge_create(edge, vertices(:,4), vertices(:,1), 41, error);   box%edges(4)  = edge
    CALL cs_edge_create(edge, vertices(:,5), vertices(:,6), 56, error);   box%edges(5)  = edge
    CALL cs_edge_create(edge, vertices(:,6), vertices(:,7), 67, error);   box%edges(6)  = edge
    CALL cs_edge_create(edge, vertices(:,7), vertices(:,8), 78, error);   box%edges(7)  = edge
    CALL cs_edge_create(edge, vertices(:,8), vertices(:,5), 85, error);   box%edges(8)  = edge
    CALL cs_edge_create(edge, vertices(:,7), vertices(:,2), 72, error);   box%edges(9)  = edge
    CALL cs_edge_create(edge, vertices(:,8), vertices(:,3), 83, error);   box%edges(10) = edge
    CALL cs_edge_create(edge, vertices(:,5), vertices(:,4), 54, error);   box%edges(11) = edge
    CALL cs_edge_create(edge, vertices(:,6), vertices(:,1), 61, error);   box%edges(12) = edge

    DEALLOCATE(edge, face)

    CALL timestop(handle)

  END SUBROUTINE cs_box_create

! *****************************************************************************
!> \brief  generates a cs_intersection data object
!> \param intersection the intersection object to be created
!> \param intxn_type intersection type
!> \param n_intersections the number of vertices of the intersection object
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_intersection_obj_create(intersection, intxn_type, n_intersections, error)

    TYPE(cs_intersection_obj), POINTER       :: intersection
    INTEGER, INTENT(IN)                      :: intxn_type
    INTEGER, INTENT(IN), OPTIONAL            :: n_intersections
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_intersection_obj_create', &
      routineP = moduleN//':'//routineN

    SELECT CASE (intxn_type)
       CASE(NO_INTXN)
          intersection%n_intxn = 0
          intersection%intxn_type = NO_INTXN
          intersection%point = (/huge_value, huge_value, huge_value/)
       CASE(PT_INTXN)
          intersection%n_intxn = 1
          intersection%intxn_type = PT_INTXN
       CASE(SEG_INTXN)
          intersection%n_intxn = 2
          intersection%intxn_type = SEG_INTXN
          ALLOCATE(intersection%segment)
       CASE(TRI_INTXN)
          intersection%n_intxn = 3
          intersection%intxn_type = TRI_INTXN
          ALLOCATE(intersection%triangle)
       CASE(RECT_INTXN)
          intersection%n_intxn = 4
          intersection%intxn_type = RECT_INTXN
          ALLOCATE(intersection%rectangle)
       CASE(PY_INTXN)
          IF (PRESENT(n_intersections)) THEN
             intersection%n_intxn = n_intersections
          ELSE
             intersection%n_intxn = 5
          END IF
          intersection%intxn_type = PY_INTXN
          ALLOCATE(intersection%polygon)
    END SELECT

  END SUBROUTINE cs_intersection_obj_create

! *****************************************************************************
!> \brief   generates a grid segment (the smallest axis-aligned segment that can be
!>          created on a grid)
!> \param grid_seg the grid segment to be generated
!> \param origin the coordinates of the vertex labeled 'A'
!> \param grid_segment_index the indices of the coordinates of the vertex labeled 'A'
!> \param dr the grid spacing
!> \param alignment the alignment of the segment
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_grid_segment_create(grid_seg, origin, grid_segment_index, dr, alignment, error)

    TYPE(cs_grid_segment), INTENT(INOUT), &
      POINTER                                :: grid_seg
    REAL(dp), DIMENSION(3), INTENT(IN)       :: origin
    INTEGER, DIMENSION(3), INTENT(IN)        :: grid_segment_index
    REAL(dp), DIMENSION(3), INTENT(IN)       :: dr
    INTEGER, INTENT(IN)                      :: alignment
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_grid_segment_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: dx, dy, dz, ox, oy, oz
    REAL(dp), DIMENSION(3)                   :: A, B

    CALL timeset(routineN,handle)

    A = 0.0_dp;    B = 0.0_dp

    grid_seg%orig = origin
    grid_seg%grid_segment_idx = grid_segment_index
    grid_seg%weight = 0
    grid_seg%alignment = alignment

    ox = origin(1);    oy = origin(2);    oz = origin(3)
    dx = dr(1);        dy = dr(2);        dz = dr(3)

    SELECT CASE (alignment)
      CASE (xaxis_aligned_segment)
         A = (/ox   , oy   , oz   /)
         B = (/ox+dx, oy   , oz   /)
      CASE (yaxis_aligned_segment)
         A = (/ox   , oy   , oz   /)
         B = (/ox   , oy+dy, oz   /)
      CASE (zaxis_aligned_segment)
         A = (/ox   , oy   , oz   /)
         B = (/ox   , oy   , oz+dz/)
    END SELECT

    ALLOCATE(grid_seg%segment)
    CALL cs_segment_create(grid_seg%segment, A, B, error)

    grid_seg%is_on = .FALSE.

    CALL timestop(handle)

  END SUBROUTINE cs_grid_segment_create

! *****************************************************************************
!> \brief   generates a grid tile (the smallest axis-aligned rectangle that can be
!>          created on a grid)
!> \param grid_tile the grid tile to be generated
!> \param origin  the coordinates of the vertex labeled 'A'
!> \param grid_tile_index the indices of the coordinates of the vertex labeled 'A'
!> \param dr the grid spacing
!> \param alignment alignment of the tile
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_grid_tile_create(grid_tile, origin, grid_tile_index, dr, alignment, error)

    TYPE(cs_grid_tile), INTENT(INOUT), &
      POINTER                                :: grid_tile
    REAL(dp), DIMENSION(3), INTENT(IN)       :: origin
    INTEGER, DIMENSION(3), INTENT(IN)        :: grid_tile_index
    REAL(dp), DIMENSION(3), INTENT(IN)       :: dr
    INTEGER, INTENT(IN)                      :: alignment
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_grid_tile_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: dx, dy, dz, ox, oy, oz
    REAL(dp), DIMENSION(3)                   :: A, B, C, D

    CALL timeset(routineN,handle)

    A = 0.0_dp;    B = 0.0_dp;    C = 0.0_dp;    D = 0.0_dp

    grid_tile%orig = origin
    grid_tile%grid_tile_idx = grid_tile_index
    grid_tile%weight = 0
    grid_tile%alignment = alignment

    ox = origin(1);    oy = origin(2);    oz = origin(3)
    dx = dr(1);        dy = dr(2);        dz = dr(3)

    SELECT CASE (alignment)
      CASE (xy_aligned_rectangle)
         A = (/ox   , oy   , oz   /)
         B = (/ox+dx, oy   , oz   /)
         C = (/ox+dx, oy+dy, oz   /)
         D = (/ox   , oy+dy, oz   /)
      CASE (yz_aligned_rectangle)
         A = (/ox   , oy   , oz   /)
         B = (/ox   , oy+dy, oz   /)
         C = (/ox   , oy+dy, oz+dz/)
         D = (/ox   , oy   , oz+dz/)
      CASE (xz_aligned_rectangle)
         A = (/ox+dx, oy   , oz   /)
         B = (/ox   , oy   , oz   /)
         C = (/ox   , oy   , oz+dz/)
         D = (/ox+dx, oy   , oz+dz/)
    END SELECT

    ALLOCATE(grid_tile%rectangle)
    CALL cs_rectangle_create(grid_tile%rectangle, A, B, C, D, error)

    grid_tile%is_on = .FALSE.

    CALL timestop(handle)

  END SUBROUTINE cs_grid_tile_create

! *****************************************************************************
!> \brief   generates a voxel (the smallest axis-aligned box that can be created
!>          on a grid)
!> \param voxel the voxel to be generated
!> \param origin the coordinates of the vertex labeled '1'
!> \param voxel_index the indices of the coordinates of the vertex labeled '1'
!> \param dr the grid spacing
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_voxel_create(voxel, origin, voxel_index, dr, error)

    TYPE(cs_voxel), INTENT(INOUT), POINTER   :: voxel
    REAL(dp), DIMENSION(3), INTENT(IN)       :: origin
    INTEGER, DIMENSION(3), INTENT(IN)        :: voxel_index
    REAL(dp), DIMENSION(3), INTENT(IN)       :: dr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_voxel_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: dx, dy, dz, ox, oy, oz
    REAL(dp), DIMENSION(3, 8)                :: vertices

    CALL timeset(routineN,handle)

    vertices = 0.0_dp

    voxel%orig = origin
    voxel%voxel_idx = voxel_index
    voxel%weight = 0

    ox = origin(1);    oy = origin(2);    oz = origin(3)
    dx = dr(1);        dy = dr(2);        dz = dr(3)

    vertices(:,1) = (/ox   , oy   , oz   /)
    vertices(:,2) = (/ox+dx, oy   , oz   /)
    vertices(:,3) = (/ox+dx, oy+dy, oz   /)
    vertices(:,4) = (/ox   , oy+dy, oz   /)
    vertices(:,5) = (/ox   , oy+dy, oz+dz/)
    vertices(:,6) = (/ox   , oy   , oz+dz/)
    vertices(:,7) = (/ox+dx, oy   , oz+dz/)
    vertices(:,8) = (/ox+dx, oy+dy, oz+dz/)

    ALLOCATE(voxel%box)
    CALL cs_box_create(voxel%box, vertices, error)

    voxel%is_on = .FALSE.

    CALL timestop(handle)

  END SUBROUTINE cs_voxel_create

! *****************************************************************************
!> \brief   releases a cs_segment data object
!> \param segment the segment to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_segment_release(segment, error)

    TYPE(cs_segment), POINTER                :: segment
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_segment_release', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(segment)) DEALLOCATE(segment)

  END SUBROUTINE cs_segment_release

! *****************************************************************************
!> \brief   releases a cs_edge data object
!> \param edge the edge to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_edge_release(edge, error)

    TYPE(cs_edge), POINTER                   :: edge
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_edge_release', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(edge)) THEN
       CALL cs_segment_release(edge%segment, error)
       DEALLOCATE(edge)
    END IF

  END SUBROUTINE cs_edge_release

! *****************************************************************************
!> \brief   releases a cs_plane data object
!> \param plane the plane to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_plane_release(plane, error)

    TYPE(cs_plane), POINTER                  :: plane
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_plane_release', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(plane)) THEN
       DEALLOCATE(plane%AB, plane%AC)
       DEALLOCATE(plane)
    END IF

  END SUBROUTINE cs_plane_release

! *****************************************************************************
!> \brief   releases a cs_triangle data object
!> \param triangle the triangle to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_triangle_release(triangle, error)

    TYPE(cs_triangle), POINTER               :: triangle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_triangle_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF(ASSOCIATED(triangle)) THEN
       DO i = 1, 3
          DEALLOCATE(triangle%edges(i)%segment)
       END DO
       CALL cs_plane_release(triangle%plane, error)
       DEALLOCATE(triangle)
    END IF

  END SUBROUTINE cs_triangle_release

! *****************************************************************************
!> \brief   releases a cs_rectangle data object
!> \param rectangle the rectangle to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_rectangle_release(rectangle, error)

    TYPE(cs_rectangle), POINTER              :: rectangle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_rectangle_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF(ASSOCIATED(rectangle)) THEN
       DO i = 1, 4
          DEALLOCATE(rectangle%edges(i)%segment)
       END DO
       CALL cs_plane_release(rectangle%plane, error)
       DEALLOCATE(rectangle)
    END IF

  END SUBROUTINE cs_rectangle_release

! *****************************************************************************
!> \brief   releases a cs_polygon data type
!> \param polygon the polygon to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_polygon_release(polygon, error)

    TYPE(cs_polygon), POINTER                :: polygon
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_polygon_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF(ASSOCIATED(polygon)) THEN
       DEALLOCATE(polygon%vertices)
       DO i = 1, polygon%n_vertices
          DEALLOCATE(polygon%edges(i)%segment)
       END DO
       DEALLOCATE(polygon%edges)
       CALL cs_plane_release(polygon%plane, error)
       DEALLOCATE(polygon)
    END IF

  END SUBROUTINE cs_polygon_release

! *****************************************************************************
!> \brief   releases a cs_face data object
!> \param face the face to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_face_release(face, error)

    TYPE(cs_face), POINTER                   :: face
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_face_release', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(face)) THEN
       CALL cs_rectangle_release(face%rectangle, error)
       DEALLOCATE(face)
    END IF

  END SUBROUTINE cs_face_release

! *****************************************************************************
!> \brief   releases a cs_box data object
!> \param box the box to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_box_release(box, error)

    TYPE(cs_box), POINTER                    :: box
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_box_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF(ASSOCIATED(box)) THEN
       DO i = 1, 12
          DEALLOCATE(box%edges(i)%segment)
       END DO
       DO i = 1, 6
          CALL cs_rectangle_release(box%faces(i)%rectangle, error)
       END DO
       DEALLOCATE(box)
    END IF

  END SUBROUTINE cs_box_release

! *****************************************************************************
!> \brief   releases a cs_intersection data object
!> \param intersection the intersection object to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_intersection_obj_release(intersection, error)

    TYPE(cs_intersection_obj), POINTER       :: intersection
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_intersection_obj_release', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(intersection)) THEN
       SELECT CASE (intersection%intxn_type)
          CASE(NO_INTXN)
             ! do nothing
          CASE(PT_INTXN)
             ! do nothing
          CASE(SEG_INTXN)
             CALL cs_segment_release(intersection%segment, error)
          CASE(TRI_INTXN)
             CALL cs_triangle_release(intersection%triangle, error)
          CASE(RECT_INTXN)
             CALL cs_rectangle_release(intersection%rectangle, error)
          CASE(PY_INTXN)
             CALL cs_polygon_release(intersection%polygon, error)
       END SELECT
       DEALLOCATE(intersection)
    END IF

  END SUBROUTINE cs_intersection_obj_release

! *****************************************************************************
!> \brief   releases a cs_grid_segment data object
!> \param grid_seg the grid_segment to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_grid_segment_release(grid_seg, error)

    TYPE(cs_grid_segment), POINTER           :: grid_seg
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_grid_segment_release', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(grid_seg)) THEN
       CALL cs_segment_release(grid_seg%segment, error)
       DEALLOCATE(grid_seg)
    END IF

  END SUBROUTINE cs_grid_segment_release

! *****************************************************************************
!> \brief   releases a cs_grid_tile data object
!> \param grid_tile the grid_tile to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_grid_tile_release(grid_tile, error)

    TYPE(cs_grid_tile), POINTER              :: grid_tile
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_grid_tile_release', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(grid_tile)) THEN
       CALL cs_rectangle_release(grid_tile%rectangle, error)
       DEALLOCATE(grid_tile)
    END IF

  END SUBROUTINE cs_grid_tile_release

! *****************************************************************************
!> \brief   releases a cs_voxel data object
!> \param voxel the voxel to be released
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cs_voxel_release(voxel, error)

    TYPE(cs_voxel), POINTER                  :: voxel
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cs_voxel_release', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(voxel)) THEN
       CALL cs_box_release(voxel%box, error)
       DEALLOCATE(voxel)
    END IF

  END SUBROUTINE cs_voxel_release

! *****************************************************************************
!> \brief sorts vertices to get a well-defined convex polygon (or rectangle):
!>        since maximum number of vertices to be sorted is 6, no fancy sorting
!>        algorithm is needed and a simple selection sort is OK.
!>        In the sorting algorithm minimum is found according to the following
!>        strict weak ordering:
!>        Vertex A precede vertex B if and only if OAxOB (the cross product of
!>        OA and OB) and the normal vector to the plane point in different di-
!>        rections.
!>        O is a user defined reference vertex
!> \param vertices the vertices to be sorted
!> \param O reference vertex
!> \param unit_normal the unit normal to the polygon's plane
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE sort_vertices(vertices, O, unit_normal, error)

    REAL(dp), DIMENSION(:,:), INTENT(INOUT)    :: vertices
    REAL(dp), DIMENSION(3), INTENT(IN)         :: O, unit_normal
    TYPE(cp_error_type), INTENT(INOUT)         :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sort_vertices', &
      routineP = moduleN//':'//routineN

    INTEGER                       :: handle, i, n_vertices, minimum_index

    CALL timeset(routineN,handle)

    n_vertices = SIZE(vertices)/3
    DO i = 1, n_vertices - 1
       minimum_index = find_minimum(vertices, O, unit_normal, i, n_vertices)
       CALL swap_columns(vertices, i, minimum_index)
    END DO

    CALL timestop(handle)

    CONTAINS
! *****************************************************************************
!> \brief ...
!> \param O ...
!> \param A ...
!> \param B ...
!> \param unit_normal ...
!> \retval a_precedes_b ...
! *****************************************************************************
    PURE FUNCTION a_precedes_b(O, A, B, unit_normal)
    REAL(dp), DIMENSION(3), INTENT(IN)       :: O, A, B, unit_normal
    LOGICAL                                  :: a_precedes_b

    REAL(dp), DIMENSION(3)                   :: OA, OAxOB, OB

        a_precedes_b = .TRUE.
        OA = A - O
        OB = B - O
        OAxOB = vector_product(OA, OB)
        IF (DOT_PRODUCT(OAxOB, unit_normal) .GT. 0) a_precedes_b = .FALSE.
   END FUNCTION
! *****************************************************************************
!> \brief ...
!> \param arr ...
!> \param O ...
!> \param unit_normal ...
!> \param start_ind ...
!> \param end_ind ...
!> \retval find_minimum ...
! *****************************************************************************
    PURE FUNCTION find_minimum(arr, O, unit_normal, start_ind, end_ind)
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: arr
    REAL(dp), DIMENSION(3), INTENT(IN)       :: O, unit_normal
    INTEGER, INTENT(IN)                      :: start_ind, end_ind
    INTEGER                                  :: find_minimum

    INTEGER                                  :: i, minimum_index
    REAL(dp), DIMENSION(3)                   :: minimum

        minimum = arr(:,start_ind)
        minimum_index = start_ind
        DO i = start_ind+1 , end_ind
           IF (.NOT. a_precedes_b(O, arr(:,i), minimum, unit_normal)) THEN
              minimum = arr(:,i)
              minimum_index = i
           END IF
        END DO
        find_minimum = minimum_index
    END FUNCTION  find_minimum
! *****************************************************************************
!> \brief ...
!> \param arr ...
!> \param col_ind1 ...
!> \param col_ind2 ...
! *****************************************************************************
    SUBROUTINE swap_columns(arr, col_ind1, col_ind2)
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: arr
    INTEGER, INTENT(IN)                      :: col_ind1, col_ind2

    REAL(dp), DIMENSION(3)                   :: column

        column = arr(:,col_ind1)
        arr(:,col_ind1) = arr(:,col_ind2)
        arr(:,col_ind2) = column
    END SUBROUTINE swap_columns

  END SUBROUTINE sort_vertices

! *****************************************************************************
!> \brief   casts a triangle to a polygon
!> \param triangle the triangle to be casted
!> \param polygon the output polygon
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cast_triangle_to_polygon(triangle, polygon, error)

    TYPE(cs_triangle), INTENT(IN), POINTER   :: triangle
    TYPE(cs_polygon), INTENT(INOUT), POINTER :: polygon
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cast_triangle_to_polygon', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp), DIMENSION(:, :), POINTER       :: py_vertices

    CALL timeset(routineN,handle)

    py_vertices => triangle%vertices
    CALL cs_polygon_create(polygon, py_vertices, error=error)

    CALL timestop(handle)

  END SUBROUTINE cast_triangle_to_polygon

! *****************************************************************************
!> \brief   casts a rectangle to a polygon
!> \param rectangle the rectangle to be casted
!> \param polygon the output polygon
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE cast_rectangle_to_polygon(rectangle, polygon, error)

    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    TYPE(cs_polygon), INTENT(INOUT), POINTER :: polygon
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cast_rectangle_to_polygon', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp), DIMENSION(:, :), POINTER       :: py_vertices

    CALL timeset(routineN,handle)

    py_vertices => rectangle%vertices
    CALL cs_polygon_create(polygon, py_vertices, error=error)

    CALL timestop(handle)

  END SUBROUTINE cast_rectangle_to_polygon

END MODULE dirichlet_bc_types
