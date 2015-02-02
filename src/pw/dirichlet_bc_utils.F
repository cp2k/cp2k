!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief utilities for perfoming Cartesian geometry operations and setting up 
!>        the Dirichlet regions of the domain
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE dirichlet_bc_utils

  USE dirichlet_bc_types,              ONLY: &
       BWROT, FWROT, NO_INTXN, PT_INTXN, PY_INTXN, SEG_INTXN, TRI_INTXN, &
       axis_aligned_box, cast_rectangle_to_polygon, cast_triangle_to_polygon, &
       cs_box, cs_box_create, cs_edge, cs_grid_segment, &
       cs_grid_segment_create, cs_grid_segment_release, cs_grid_tile, &
       cs_grid_tile_create, cs_grid_tile_release, cs_intersection_obj, &
       cs_intersection_obj_create, cs_intersection_obj_release, cs_plane, &
       cs_plane_create, cs_plane_release, cs_polygon, cs_polygon_create, &
       cs_polygon_release, cs_rectangle, cs_rectangle_create, &
       cs_rectangle_release, cs_segment, cs_segment_create, &
       cs_segment_release, cs_triangle_create, cs_voxel, cs_voxel_create, &
       cs_voxel_release, no_grid_point, single_grid_point, sort_vertices, &
       tile_p_type, xaxis_aligned_segment, xy_aligned_rectangle, &
       xz_aligned_rectangle, yaxis_aligned_segment, yz_aligned_rectangle, &
       zaxis_aligned_segment
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: degree,&
                                             twopi
  USE mathlib,                         ONLY: angle,&
                                             rotate_vector,&
                                             vector_product
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_types,                        ONLY: pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dirichlet_bc_utils'

  PUBLIC setup_grid_axes, &
         linspace,&
         find_fst_ge,&
         find_lst_le,&
         rotate_rectangle,&
         convert_to_aa_rectangle,&
         partition_aa_cuboid_into_tiles,&
         partition_aa_rectangle_into_tiles,&
         voxelize_rectangle,&
         voxelize_aa_cuboid

  REAL(dp), PARAMETER, PRIVATE         :: small_value = 1.0E-8_dp
  REAL(dp), PARAMETER, PRIVATE         :: huge_value  = HUGE(0.0_dp)

CONTAINS

! *****************************************************************************
!> \brief returns the global axes and the portion of the axes that are local to 
!>        the current mpi rank
!> \param pw_grid plane wave grid   
!> \param x_glbl x grid vetor of the simulation box
!> \param y_glbl y grid vetor of the simulation box
!> \param z_glbl z grid vetor of the simulation box
!> \param x_locl x grid vetor of the simulation box local to this process
!> \param y_locl y grid vetor of the simulation box local to this process
!> \param z_locl z grid vetor of the simulation box local to this process
!> \param error cp2k error 
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE setup_grid_axes(pw_grid, x_glbl, y_glbl, z_glbl, x_locl, y_locl, z_locl, error)

    TYPE(pw_grid_type), INTENT(IN), POINTER  :: pw_grid
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(OUT)                            :: x_glbl, y_glbl, z_glbl, &
                                                x_locl, y_locl, z_locl
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_grid_axes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: glb1, glb2, glb3, gub1, gub2, &
                                                gub3, handle, i, lb1, lb2, &
                                                lb3, ub1, ub2, ub3
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: gindx, gindy, gindz, lindx, &
                                                lindy, lindz
    INTEGER, DIMENSION(2, 3)                 :: bounds, bounds_local
    INTEGER, DIMENSION(3)                    :: npts, npts_local
    REAL(dp), DIMENSION(3)                   :: dr

    CALL timeset(routineN,handle)

    dr = pw_grid%dr
    bounds = pw_grid%bounds
    bounds_local = pw_grid%bounds_local
    npts = pw_grid%npts
    npts_local = pw_grid%npts_local

! local and global lower and upper bounds
    glb1 = bounds(1,1);          gub1 = bounds(2,1)
    glb2 = bounds(1,2);          gub2 = bounds(2,2)
    glb3 = bounds(1,3);          gub3 = bounds(2,3)
    lb1  = bounds_local(1,1);    ub1  = bounds_local(2,1)
    lb2  = bounds_local(1,2);    ub2  = bounds_local(2,2)
    lb3  = bounds_local(1,3);    ub3  = bounds_local(2,3)

    ALLOCATE(lindx( lb1 : ub1 ), lindy( lb2 : ub2 ), lindz( lb3 : ub3 ))
    ALLOCATE(gindx(glb1 : gub1), gindy(glb2 : gub2), gindz(glb3 : gub3))
    ALLOCATE(x_locl( lb1 : ub1 ), y_locl( lb2 : ub2 ), z_locl( lb3 : ub3 ))
    ALLOCATE(x_glbl(glb1 : gub1), y_glbl(glb2 : gub2), z_glbl(glb3 : gub3))

    gindx(:) = (/ (i, i=0, npts(1) - 1) /)
    gindy(:) = (/ (i, i=0, npts(2) - 1) /)
    gindz(:) = (/ (i, i=0, npts(3) - 1) /)
    lindx(:) = (/ (i, i=0, npts_local(1) - 1) /)
    lindy(:) = (/ (i, i=0, npts_local(2) - 1) /)
    lindz(:) = (/ (i, i=0, npts_local(3) - 1) /)

    x_glbl(:) = gindx * dr(1)
    y_glbl(:) = gindy * dr(2)
    z_glbl(:) = gindz * dr(3)

! map [0 .. (npts_local-1)] --> [lb .. ub]
    IF (lb1 .EQ. ub1) THEN 
       lindx = lb1 
    ELSE
       lindx = lindx * ((ub1 - lb1) / (npts_local(1) - 1)) + lb1
    END IF 
    IF (lb2 .EQ. ub2) THEN 
       lindy = lb2
    ELSE
       lindy = lindy * ((ub2 - lb2) / (npts_local(2) - 1)) + lb2
    END IF 
    IF (lb3 .EQ. ub3) THEN 
       lindz = lb3 
    ELSE
       lindz = lindz * ((ub3 - lb3) / (npts_local(3) - 1)) + lb3
    END IF 

    x_locl(:) = x_glbl(lindx)
    y_locl(:) = y_glbl(lindy)
    z_locl(:) = z_glbl(lindz)

    CALL timestop(handle)

  END SUBROUTINE setup_grid_axes

! *****************************************************************************
!> \brief Checks if a point is the infinity point
!> \param point the point to be checked 
!> \retval is_infinity_point whether or not the point is the infinity point
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  PURE FUNCTION is_infinity_point(point)

    REAL(dp), DIMENSION(3), INTENT(IN)       :: point
    LOGICAL                                  :: is_infinity_point

    IF ((point(1) .EQ. huge_value) .AND. &
        (point(2) .EQ. huge_value) .AND. &
        (point(3) .EQ. huge_value)) THEN
       is_infinity_point = .TRUE.
    ELSE
       is_infinity_point = .FALSE.
    END IF

  END FUNCTION is_infinity_point

! *****************************************************************************
!> \brief Checks if two points are equal
!> \param point1 the point to be checked 
!> \param point2 the point to be checked 
!> \retval are_equal whether or not the points are equal
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  PURE FUNCTION are_equal(point1, point2)

    REAL(dp), DIMENSION(3), INTENT(IN)       :: point1, point2
    LOGICAL                                  :: are_equal

    IF ((point1(1) .EQ. point2(1)) .AND. &
        (point1(2) .EQ. point2(2)) .AND. &
        (point1(3) .EQ. point2(3))) THEN
       are_equal = .TRUE.
    ELSE
       are_equal = .FALSE.
    END IF

  END FUNCTION are_equal

! *****************************************************************************
!> \brief  Creteas an evenly-spaced 1D array between two values
!> \param a starting value
!> \param b end value
!> \param N number of evenly-spaced points between a and b 
!> \retval arr array containg N evenly-spaced points between a and b
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  FUNCTION linspace(a, b, N) RESULT(arr)

    REAL(dp), INTENT(IN)                     :: a, b
    INTEGER, INTENT(IN)                      :: N
    REAL(dp), DIMENSION(N)                   :: arr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(dp)                                 :: dx

    dx = (b-a)/(N-1)
    arr = (/(a+(i-1)*dx, i=1,N)/)

  END FUNCTION linspace

! *****************************************************************************
!> \brief  Finds the first element in a 1D array that is greater than or equal to a 
!>         given value
!> \param arr the 1D array to be searched
!> \param val the value to be found
!> \retval loc the index of the first element of arr greater than or equal to val
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  FUNCTION find_fst_ge(arr, val) RESULT(loc)

    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: arr
    REAL(dp), INTENT(IN)                     :: val
    INTEGER                                  :: loc

    CHARACTER(LEN=*), PARAMETER :: routineN = 'find_fst_ge', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, lb, ub
    LOGICAL                                  :: found

    CALL timeset(routineN,handle)

    found = .FALSE.
    lb = LBOUND(arr,1);    ub = UBOUND(arr,1)

    DO i = lb, ub
       IF (arr(i) .GE. val) THEN
          loc = i
          found = .TRUE.
          EXIT
       END IF
    END DO
! if val is larger than the last element in arr return the index of the last element
    IF (found .EQV. .FALSE.) loc = ub

    CALL timestop(handle)

  END FUNCTION find_fst_ge

! *****************************************************************************
!> \brief  Finds the last element in a 1D array that is less than or equal to a 
!>         given value
!> \param arr the 1D array to be searched
!> \param val the value to be found
!> \retval loc the index of the last element of arr less than or equal to val
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  FUNCTION find_lst_le(arr, val) RESULT(loc)
     
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: arr
    REAL(dp), INTENT(IN)                     :: val
    INTEGER                                  :: loc

    CHARACTER(LEN=*), PARAMETER :: routineN = 'find_lst_le', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, lb, ub
    LOGICAL                                  :: found

    CALL timeset(routineN,handle)

    found = .FALSE. 
    lb = LBOUND(arr,1);    ub = UBOUND(arr,1)

    DO i = lb, ub
       IF (arr(i) .LE. val) THEN
          loc = i
          found = .TRUE.
       END IF
    END DO
! if val is smaller than the first element in arr return the index of the first element
    IF (found .EQV. .FALSE.) loc = lb

    CALL timestop(handle)

  END FUNCTION find_lst_le

! *****************************************************************************
!> \brief Checks if a point in 3D space is located inside a rectangle
!> \param point a point in 3D Cartesian space
!> \param rectangle rectangle to be tested
!> \retval is_inside_rectangle  whether or not point is inside rectangle
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  PURE FUNCTION is_inside_rectangle(point, rectangle)

    REAL(dp), DIMENSION(3), INTENT(IN)       :: point
    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    LOGICAL                                  :: is_inside_rectangle

    INTEGER                                  :: i
    REAL(dp)                                 :: costheta, prodnorm, sum_angle
    REAL(dp), DIMENSION(3)                   :: v1, v2
    REAL(dp), DIMENSION(3, 0:3)              :: vertices

    vertices = rectangle%vertices

    sum_angle = 0.0_dp
    DO i = 0, 3 
       v1 = vertices(:,i) - point
       v2 = vertices(:,MOD(i+1,4)) - point

       prodnorm = SQRT(SUM(v1**2)) * SQRT(SUM(v2**2))
       IF (prodnorm .LE. small_value) THEN
          sum_angle = twopi
          EXIT
       ELSE
          costheta = DOT_PRODUCT(v1, v2) / prodnorm
       END IF
       sum_angle = sum_angle + ACOS(costheta)
    END DO

    IF (ABS(sum_angle - twopi) .LE. small_value) THEN 
       is_inside_rectangle = .TRUE.
    ELSE
       is_inside_rectangle = .FALSE.
    END IF 

  END FUNCTION is_inside_rectangle

! *****************************************************************************
!> \brief Checks if a point in 3D space is located inside a polygon
!> \param point a point in 3D Cartesian space
!> \param polygon polygon to be tested
!> \retval is_inside_polygon  whether or not point is inside polygon
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  PURE FUNCTION is_inside_polygon(point, polygon)

    REAL(dp), DIMENSION(3), INTENT(IN)       :: point
    TYPE(cs_polygon), INTENT(IN), POINTER    :: polygon
    LOGICAL                                  :: is_inside_polygon

    INTEGER                                  :: i, n_vertices
    REAL(dp)                                 :: costheta, prodnorm, sum_angle
    REAL(dp), DIMENSION(3)                   :: v1, v2
    REAL(dp) :: vertices(3, 0:polygon%n_vertices-1)

    vertices = polygon%vertices
    n_vertices = polygon%n_vertices

    sum_angle = 0.0_dp
    DO i = 0, n_vertices-1
       v1 = vertices(:,i) - point
       v2 = vertices(:,MOD(i+1, n_vertices)) - point

       prodnorm = SQRT(SUM(v1**2)) * SQRT(SUM(v2**2))
       IF (prodnorm .LE. small_value) THEN
          sum_angle = twopi
          EXIT
       ELSE
          costheta = DOT_PRODUCT(v1, v2) / prodnorm
       END IF
       sum_angle = sum_angle + ACOS(costheta)
    END DO

    IF (ABS(sum_angle - twopi) .LE. small_value) THEN 
       is_inside_polygon = .TRUE.
    ELSE
       is_inside_polygon = .FALSE.
    END IF

  END FUNCTION is_inside_polygon

! *****************************************************************************
!> \brief Removes duplicate vertices from an array of vertices
!> \param vertices_in array containing the coordinates of a set of vertices 
!> \param vertices_out array containing unique vertices in vertices_in
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE remove_duplicate_vertices(vertices_in, vertices_out, error)

    REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
      INTENT(IN)                             :: vertices_in
    REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
      INTENT(OUT)                            :: vertices_out
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_duplicate_vertices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, n_cols, n_unique
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: temp_arr

    CALL timeset(routineN,handle)

    n_cols = SIZE(vertices_in)/3
    ALLOCATE(temp_arr(3,n_cols))

    n_unique = 1
    temp_arr(:,1) = vertices_in(:,1)

    outer: DO i = 2, n_cols
       DO j = 1, n_unique
          IF (SQRT(SUM((temp_arr(:,j) - vertices_in(:,i))**2)) .LE. small_value) CYCLE outer
       END DO
       n_unique = n_unique + 1
       temp_arr(:,n_unique) = vertices_in(:,i)
    END DO outer

    ALLOCATE(vertices_out(3,n_unique))
    vertices_out(:,1:n_unique) = temp_arr(:,1:n_unique)

    DEALLOCATE(temp_arr)

    CALL timestop(handle)

  END SUBROUTINE remove_duplicate_vertices

! *****************************************************************************
!> \brief   Finds the extents of the axis-aligned bounding box of a given rectangle
!> \param rectangle input rectangle
!> \param x_glbl x grid vetor of the simulation box
!> \param y_glbl y grid vetor of the simulation box
!> \param z_glbl z grid vetor of the simulation box
!> \param x_locl x grid vetor of the simulation box local to this process
!> \param y_locl y grid vetor of the simulation box local to this process
!> \param z_locl z grid vetor of the simulation box local to this process
!> \param aabb_extents global extents
!> \param aabb_extents_local local extents
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE compute_rectangle_aabb_extents(rectangle, x_glbl, y_glbl, z_glbl, &
                           x_locl, y_locl, z_locl, aabb_extents, aabb_extents_local, error)

    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl, &
                                                x_locl, y_locl, z_locl
    INTEGER, DIMENSION(2, 3), INTENT(OUT)    :: aabb_extents, &
                                                aabb_extents_local
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'compute_rectangle_aabb_extents', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: glb1, glb2, glb3, gub1, gub2, &
                                                gub3, handle, lb1, lb2, lb3, &
                                                ub1, ub2, ub3
    REAL(dp)                                 :: maxx, maxy, maxz, minx, miny, &
                                                minz

    CALL timeset(routineN,handle)

    minx = MINVAL(rectangle%vertices(1,:))
    maxx = MAXVAL(rectangle%vertices(1,:))
    miny = MINVAL(rectangle%vertices(2,:))
    maxy = MAXVAL(rectangle%vertices(2,:))
    minz = MINVAL(rectangle%vertices(3,:))
    maxz = MAXVAL(rectangle%vertices(3,:))

    glb1 = LBOUND(x_glbl,1);    gub1 = UBOUND(x_glbl,1)
    glb2 = LBOUND(y_glbl,1);    gub2 = UBOUND(y_glbl,1)
    glb3 = LBOUND(z_glbl,1);    gub3 = UBOUND(z_glbl,1)
    lb1  = LBOUND(x_locl,1);    ub1  = UBOUND(x_locl,1)
    lb2  = LBOUND(y_locl,1);    ub2  = UBOUND(y_locl,1)
    lb3  = LBOUND(z_locl,1);    ub3  = UBOUND(z_locl,1)

    aabb_extents_local(1,1) = find_lst_le(x_locl, minx)
    aabb_extents_local(2,1) = find_fst_ge(x_locl, maxx)
    aabb_extents_local(1,2) = find_lst_le(y_locl, miny)
    aabb_extents_local(2,2) = find_fst_ge(y_locl, maxy)
    aabb_extents_local(1,3) = find_lst_le(z_locl, minz)
    aabb_extents_local(2,3) = find_fst_ge(z_locl, maxz)

    IF (aabb_extents_local(1,1) .EQ. aabb_extents_local(2,1)) THEN
       IF (aabb_extents_local(1,1) .EQ. ub1) THEN 
          aabb_extents_local(1,1) = aabb_extents_local(1,1) - 1
       ELSE
          aabb_extents_local(2,1) = aabb_extents_local(2,1) + 1
       END IF
    END IF
    IF (aabb_extents_local(1,2) .EQ. aabb_extents_local(2,2)) THEN
       IF (aabb_extents_local(1,2) .EQ. ub2) THEN 
          aabb_extents_local(1,2) = aabb_extents_local(1,2) - 1
       ELSE
          aabb_extents_local(2,2) = aabb_extents_local(2,2) + 1
       END IF
    END IF
    IF (aabb_extents_local(1,3) .EQ. aabb_extents_local(2,3)) THEN
       IF (aabb_extents_local(1,3) .EQ. ub3) THEN 
          aabb_extents_local(1,3) = aabb_extents_local(1,3) - 1
       ELSE
          aabb_extents_local(2,3) = aabb_extents_local(2,3) + 1
       END IF
    END IF

    aabb_extents(1,1) = find_lst_le(x_glbl, minx)
    aabb_extents(2,1) = find_fst_ge(x_glbl, maxx)
    aabb_extents(1,2) = find_lst_le(y_glbl, miny)
    aabb_extents(2,2) = find_fst_ge(y_glbl, maxy)
    aabb_extents(1,3) = find_lst_le(z_glbl, minz)
    aabb_extents(2,3) = find_fst_ge(z_glbl, maxz)

    IF (aabb_extents(1,1) .EQ. aabb_extents(2,1)) THEN
       IF (aabb_extents(1,1) .EQ. gub1) THEN 
          aabb_extents(1,1) = aabb_extents(1,1) - 1
       ELSE
          aabb_extents(2,1) = aabb_extents(2,1) + 1
       END IF
    END IF
    IF (aabb_extents(1,2) .EQ. aabb_extents(2,2)) THEN
       IF (aabb_extents(1,2) .EQ. gub2) THEN 
          aabb_extents(1,2) = aabb_extents(1,2) - 1
       ELSE
          aabb_extents(2,2) = aabb_extents(2,2) + 1
       END IF
    END IF
    IF (aabb_extents(1,3) .EQ. aabb_extents(2,3)) THEN
       IF (aabb_extents(1,3) .EQ. gub3) THEN 
          aabb_extents(1,3) = aabb_extents(1,3) - 1
       ELSE
          aabb_extents(2,3) = aabb_extents(2,3) + 1
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE compute_rectangle_aabb_extents

! *****************************************************************************
!> \brief Rotates a rectangle such that it becomes parallel to a given plane,  
!>        then projects it onto the plane if it is asked to do so.
!> \param rectangle rectangle to be rotated
!> \param plane the given plane
!> \param rectangle_rot the rotated rectangle  
!> \param phi the rotation angle
!> \param rot_axis the rotation axis
!> \param project whether or not to project
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE rotate_rectangle_onto_plane(rectangle, plane, rectangle_rot, phi, rot_axis, project, error)

    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    TYPE(cs_plane), INTENT(IN), POINTER      :: plane
    TYPE(cs_rectangle), INTENT(INOUT), &
      POINTER                                :: rectangle_rot
    REAL(dp), INTENT(OUT), OPTIONAL          :: phi
    REAL(dp), DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: rot_axis
    LOGICAL, INTENT(IN), OPTIONAL            :: project
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rotate_rectangle_onto_plane', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: intern_project
    REAL(dp)                                 :: intern_phi, n1dotn2, norm_n1, &
                                                norm_n1xn2, norm_n2
    REAL(dp), DIMENSION(3)                   :: intern_rot_axis, n1, n1xn2, n2
    TYPE(cs_rectangle), POINTER              :: rectangle_tmp

    CALL timeset(routineN,handle)
    intern_project = .FALSE.
    IF (PRESENT(project)) intern_project = project

    n1 = rectangle%plane%unit_normal
    n2 = plane%unit_normal

    norm_n1 = SQRT(SUM(n1**2))
    norm_n2 = SQRT(SUM(n2**2))
   
    n1dotn2 = DOT_PRODUCT(n1, n2)/(norm_n1 * norm_n2)
! check if the rectangle is parallel to the plane
    IF (ABS(ABS(n1dotn2) - 1) .LE. small_value) THEN
! project it onto the plane if the user wants
       IF (intern_project .EQV. .TRUE.) THEN
          CALL project_rectangle_onto_plane(rectangle, plane, rectangle_rot, error)
       ELSE
          CALL cs_rectangle_create(rectangle_rot, rectangle%A, rectangle%B, &
                                     rectangle%C, rectangle%D, error)
!          rectangle_rot = rectangle
       END IF
       intern_rot_axis = n2
       intern_phi = 0.0_dp
    ELSE
       n1xn2 = vector_product(n1, n2) 
       norm_n1xn2 = SQRT(SUM(n1xn2**2))
       intern_rot_axis = n1xn2 / norm_n1xn2
       intern_phi = angle(n1, n2)

       ALLOCATE(rectangle_tmp)
       CALL rotate_rectangle(rectangle, intern_phi, intern_rot_axis, FWROT, rectangle_tmp, error)

       IF (intern_project .EQV. .TRUE.) THEN
          CALL project_rectangle_onto_plane(rectangle_tmp, plane, rectangle_rot, error)
          CALL cs_rectangle_release(rectangle_tmp, error=error)
       ELSE
          rectangle_rot = rectangle_tmp
          DEALLOCATE(rectangle_tmp)
       END IF
    END IF
   
    IF (PRESENT(phi)) phi = intern_phi
    IF (PRESENT(rot_axis)) rot_axis = intern_rot_axis

    CALL timestop(handle)

  END SUBROUTINE rotate_rectangle_onto_plane 

! *****************************************************************************
!> \brief Given the rotation angle and axis, rotates a rectangle. 
!> \param rectangle rectangle to be rotated
!> \param phi the rotation angle
!> \param rot_axis the rotation axis 
!> \param direction forward or backward rotation
!> \param rectangle_rot the rotated rectangle  
!> \param error cp2k error
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE rotate_rectangle(rectangle, phi, rot_axis, direction, rectangle_rot, error)

    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    REAL(dp), INTENT(IN)                     :: phi
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rot_axis
    INTEGER, INTENT(IN)                      :: direction
    TYPE(cs_rectangle), INTENT(INOUT), &
      POINTER                                :: rectangle_rot
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rotate_rectangle', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: intern_phi
    REAL(dp), DIMENSION(3)                   :: Aprime, Bprime, Cprime, Dprime

    CALL timeset(routineN,handle)

    IF (direction .EQ. FWROT) THEN 
       intern_phi = phi
    ELSE IF (direction .EQ. BWROT) THEN
       intern_phi = - phi
    END IF

    Aprime = rotate_vector(rectangle%A, intern_phi, rot_axis)
    Bprime = rotate_vector(rectangle%B, intern_phi, rot_axis)
    Cprime = rotate_vector(rectangle%C, intern_phi, rot_axis)
    Dprime = rotate_vector(rectangle%D, intern_phi, rot_axis)
    CALL cs_rectangle_create(rectangle_rot, Aprime, Bprime, Cprime, Dprime, error)

    CALL timestop(handle)

  END SUBROUTINE rotate_rectangle 

! *****************************************************************************
!> \brief Converts an arbitrary rectangle to an axis-aligned one in 3 steps: 
!>   (1) choose a plane parallel to one of the coordinate planes (let's call it plane)
!>   (2) rotate rectangle onto plane (let's name the obtained rectangle rectangle_rot)
!>   (3) rotate rectangle_rot on plane to make it axis-aligned
!> \param rectangle rectangle to be converted
!> \param rectangle_aa the axis-aligned rectangle  
!> \param phi1 the first rotation angle 
!> \param rot_axis1 the first rotation axis  
!> \param phi2 the second rotation angle 
!> \param rot_axis2 the second rotation axis 
!> \param error cp2k error
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE convert_to_aa_rectangle(rectangle, rectangle_aa, phi1, rot_axis1, phi2, rot_axis2, error)

    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    TYPE(cs_rectangle), INTENT(INOUT), &
      POINTER                                :: rectangle_aa
    REAL(dp), INTENT(OUT)                    :: phi1
    REAL(dp), DIMENSION(3), INTENT(OUT)      :: rot_axis1
    REAL(dp), INTENT(OUT)                    :: phi2
    REAL(dp), DIMENSION(3), INTENT(OUT)      :: rot_axis2
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_to_aa_rectangle', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, l
    REAL(dp)                                 :: AB_e1_ang, AB_e2_ang, &
                                                AB_e3_ang
    REAL(dp), DIMENSION(3)                   :: AB, e1, e2, e3, O, p_normal
    TYPE(cs_plane), POINTER                  :: plane
    TYPE(cs_rectangle), POINTER              :: rectangle_rot

    CALL timeset(routineN,handle)
    ALLOCATE(plane, rectangle_rot)

! choose a plane parallel to one of the coordinate planes
    O = 0.0_dp
    l = MAXLOC(ABS(rectangle%plane%unit_normal),1)
    p_normal = 0.0_dp
    p_normal(l) = 1.0_dp
    CALL cs_plane_create(plane, O, p_normal, error)
! rotate the rectangle onto the plane
    CALL rotate_rectangle_onto_plane(rectangle, plane, rectangle_rot, phi1, rot_axis1, .FALSE., error)

! test if the rotated rectangle is axis-aligned
    AB = rectangle_rot%B - rectangle_rot%A
    e1 = (/1.0_dp, 0.0_dp, 0.0_dp/)
    e2 = (/0.0_dp, 1.0_dp, 0.0_dp/)
    e3 = (/0.0_dp, 0.0_dp, 1.0_dp/)
    AB_e1_ang = angle(e1, AB)
    AB_e2_ang = angle(e2, AB)
    AB_e3_ang = angle(e3, AB)

    IF ( (ABS((AB_e1_ang + AB_e2_ang)*degree - 180.0_dp) .LE. small_value) .OR. & 
         (ABS((AB_e2_ang + AB_e3_ang)*degree - 180.0_dp) .LE. small_value) .OR. &
         (ABS((AB_e1_ang + AB_e3_ang)*degree - 180.0_dp) .LE. small_value) ) THEN
       rot_axis2 = rot_axis1
       phi2 = 0.0_dp
       CALL cs_rectangle_create(rectangle_aa, rectangle_rot%A, rectangle_rot%B, &
                                  rectangle_rot%C, rectangle_rot%D, error)
    ELSE
       IF (ABS(ABS(AB_e1_ang*degree) - 90.0_dp) .LE. small_value) THEN
          rot_axis2 = e1
          IF (((AB(2) .GT. 0.0_dp) .AND. (AB(3) .GT. 0.0_dp)) .OR. &
              ((AB(2) .LT. 0.0_dp) .AND. (AB(3) .GT. 0.0_dp))) THEN
             phi2 = - AB_e2_ang
          ELSE
             phi2 = AB_e2_ang
          END IF
       ELSE IF (ABS(ABS(AB_e2_ang*degree) - 90.0_dp) .LE. small_value) THEN
          rot_axis2 = e2
          IF (((AB(1) .GT. 0.0_dp) .AND. (AB(3) .GT. 0.0_dp)) .OR. &
              ((AB(1) .LT. 0.0_dp) .AND. (AB(3) .GT. 0.0_dp))) THEN
             phi2 = - AB_e3_ang
          ELSE
             phi2 = AB_e3_ang
          END IF
       ELSE IF (ABS(ABS(AB_e3_ang*degree) - 90.0_dp) .LE. small_value) THEN
          rot_axis2 = e3
          IF (((AB(1) .GT. 0.0_dp) .AND. (AB(2) .GT. 0.0_dp)) .OR. &
              ((AB(1) .LT. 0.0_dp) .AND. (AB(2) .GT. 0.0_dp))) THEN
             phi2 = - AB_e1_ang
          ELSE
             phi2 = AB_e1_ang
          END IF
       END IF
       CALL rotate_rectangle(rectangle_rot, phi2, rot_axis2, FWROT , rectangle_aa, error)
    END IF

    CALL cs_rectangle_release(rectangle_rot, error=error)
    CALL cs_plane_release(plane, error=error)
    CALL timestop(handle)

  END SUBROUTINE convert_to_aa_rectangle 

! *****************************************************************************
!> \brief Projects a rectangle onto a plane 
!> \param rectangle the rectangle to be projected
!> \param plane the given plane
!> \param rectangle_proj the projected rectangle
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE project_rectangle_onto_plane(rectangle, plane, rectangle_proj, error)

    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    TYPE(cs_plane), INTENT(IN), POINTER      :: plane
    TYPE(cs_rectangle), INTENT(INOUT), &
      POINTER                                :: rectangle_proj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'project_rectangle_onto_plane', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: ndotPA, ndotPB, ndotPC, ndotPD
    REAL(dp), DIMENSION(3)                   :: Aprime, Bprime, Cprime, &
                                                Dprime, P

    CALL timeset(routineN,handle)

    P = plane%AB%A

    ndotPA = DOT_PRODUCT(rectangle%A - P, plane%unit_normal)
    ndotPB = DOT_PRODUCT(rectangle%B - P, plane%unit_normal)
    ndotPC = DOT_PRODUCT(rectangle%C - P, plane%unit_normal)
    ndotPD = DOT_PRODUCT(rectangle%D - P, plane%unit_normal)

    Aprime = rectangle%A - ndotPA * plane%unit_normal
    Bprime = rectangle%B - ndotPB * plane%unit_normal
    Cprime = rectangle%C - ndotPC * plane%unit_normal
    Dprime = rectangle%D - ndotPD * plane%unit_normal

    CALL cs_rectangle_create(rectangle_proj, Aprime, Bprime, Cprime, Dprime, error)

    CALL timestop(handle)

  END SUBROUTINE project_rectangle_onto_plane 

! *****************************************************************************
!> \brief    Partitions an axis-aligned rectangle into tiles. The rectangle is either
!>    a gate or a tile itself.
!> \param rectangle axis-aligned rectangle to be partitioned
!> \param x_glbl x grid vetor of the simulation box
!> \param y_glbl y grid vetor of the simulation box
!> \param z_glbl z grid vetor of the simulation box
!> \param n_prtn vetor of size 3 specifying the number of times that the x, y and
!>               z interval (defining the region) should be partitioned into
!> \param tiles the obtained tiles    
!> \param error cp2k error
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE partition_aa_rectangle_into_tiles(rectangle, x_glbl, y_glbl, z_glbl, n_prtn, tiles, error)

    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl
    INTEGER, DIMENSION(3), INTENT(IN)        :: n_prtn
    TYPE(tile_p_type), DIMENSION(:), &
      INTENT(OUT), POINTER                   :: tiles
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'partition_aa_rectangle_into_tiles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, xxtnt_fst, &
                                                xxtnt_ult, yxtnt_fst, &
                                                yxtnt_ult, zxtnt_fst, &
                                                zxtnt_ult
    INTEGER, DIMENSION(n_prtn(1))            :: xprtn_lbind, xprtn_npts, &
                                                xprtn_ubind
    INTEGER, DIMENSION(n_prtn(2))            :: yprtn_lbind, yprtn_npts, &
                                                yprtn_ubind
    INTEGER, DIMENSION(n_prtn(3))            :: zprtn_lbind, zprtn_npts, &
                                                zprtn_ubind
    REAL(dp)                                 :: step
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: xglbl_sym, yglbl_sym, &
                                                zglbl_sym
    REAL(dp), DIMENSION(2)                   :: x_xtnt, y_xtnt, z_xtnt
    REAL(dp), DIMENSION(3)                   :: A, B, C, D
    REAL(dp), DIMENSION(n_prtn(1))           :: xprtn_lb, xprtn_ub
    REAL(dp), DIMENSION(n_prtn(2))           :: yprtn_lb, yprtn_ub
    REAL(dp), DIMENSION(n_prtn(3))           :: zprtn_lb, zprtn_ub

    CALL timeset(routineN,handle)

    ALLOCATE(xglbl_sym(2*SIZE(x_glbl)-1), yglbl_sym(2*SIZE(y_glbl)-1), zglbl_sym(2*SIZE(z_glbl)-1))

! extend the global axes and make them symmetric
    xglbl_sym(:) = (/ (- x_glbl(i), i = UBOUND(x_glbl,1), LBOUND(x_glbl,1)+1, -1), x_glbl /)
    yglbl_sym(:) = (/ (- y_glbl(i), i = UBOUND(y_glbl,1), LBOUND(y_glbl,1)+1, -1), y_glbl /)
    zglbl_sym(:) = (/ (- z_glbl(i), i = UBOUND(z_glbl,1), LBOUND(z_glbl,1)+1, -1), z_glbl /)

! find the extents of the rectangle in all three directions
    x_xtnt(1) = MINVAL(rectangle%vertices(1,:));    x_xtnt(2) = MAXVAL(rectangle%vertices(1,:))
    y_xtnt(1) = MINVAL(rectangle%vertices(2,:));    y_xtnt(2) = MAXVAL(rectangle%vertices(2,:))
    z_xtnt(1) = MINVAL(rectangle%vertices(3,:));    z_xtnt(2) = MAXVAL(rectangle%vertices(3,:))

! find the position of the both ends of xtnt in the symmetric global axes
    xxtnt_fst = find_fst_ge(xglbl_sym, x_xtnt(1))
!    xxtnt_ult = find_fst_ge(xglbl_sym, x_xtnt(2))
    xxtnt_ult = find_lst_le(xglbl_sym, x_xtnt(2))
    yxtnt_fst = find_fst_ge(yglbl_sym, y_xtnt(1))
!    yxtnt_ult = find_fst_ge(yglbl_sym, y_xtnt(2))
    yxtnt_ult = find_lst_le(yglbl_sym, y_xtnt(2))
    zxtnt_fst = find_fst_ge(zglbl_sym, z_xtnt(1))
!    zxtnt_ult = find_fst_ge(zglbl_sym, z_xtnt(2))
    zxtnt_ult = find_lst_le(zglbl_sym, z_xtnt(2))

! devide the integer sequence of indices {xtnt_fst .. xtnt_ult} into n_prtn partitions
    step = REAL((xxtnt_ult-xxtnt_fst),kind=dp)/REAL(n_prtn(1), kind=dp)
    xprtn_lbind = xxtnt_fst + NINT( (/(i, i=0, n_prtn(1)-1)/) * step ) + 1   ! lower bounds
    xprtn_ubind = xxtnt_fst + NINT( ((/(i, i=0, n_prtn(1)-1)/) + 1) * step ) ! upper bounds
    IF (step .EQ. 0.0_dp) xprtn_lbind = xprtn_ubind
! count the number of grid points
    xprtn_npts = xprtn_ubind - xprtn_lbind + 1
! find the corresponding values on the symmetric axes
    xprtn_lb = xglbl_sym(xprtn_lbind)
    xprtn_ub = xglbl_sym(xprtn_ubind)
! make sure that the tiles will cover the whole rectangle
    xprtn_lb(1) = MIN(xprtn_lb(1), x_xtnt(1))
    xprtn_ub(n_prtn(1)) = MIN(xprtn_ub(n_prtn(1)), x_xtnt(2))

    step = REAL((yxtnt_ult-yxtnt_fst),kind=dp)/REAL(n_prtn(2), kind=dp)
    yprtn_lbind = yxtnt_fst + NINT( (/(i, i=0, n_prtn(2)-1)/) * step ) + 1
    yprtn_ubind = yxtnt_fst + NINT( ((/(i, i=0, n_prtn(2)-1)/) + 1) * step )
    IF (step .EQ. 0.0_dp) yprtn_lbind = yprtn_ubind
    yprtn_npts = yprtn_ubind - yprtn_lbind + 1
    yprtn_lb = yglbl_sym(yprtn_lbind)
    yprtn_ub = yglbl_sym(yprtn_ubind)
    yprtn_lb(1) = MIN(yprtn_lb(1), y_xtnt(1))
    yprtn_ub(n_prtn(2)) = MIN(yprtn_ub(n_prtn(2)), y_xtnt(2))

    step = REAL((zxtnt_ult-zxtnt_fst),kind=dp)/REAL(n_prtn(3), kind=dp)
    zprtn_lbind = zxtnt_fst + NINT( (/(i, i=0, n_prtn(3)-1)/) * step ) + 1
    zprtn_ubind = zxtnt_fst + NINT( ((/(i, i=0, n_prtn(3)-1)/) + 1) * step )
    IF (step .EQ. 0.0_dp) zprtn_lbind = zprtn_ubind
    zprtn_npts = zprtn_ubind - zprtn_lbind + 1
    zprtn_lb = zglbl_sym(zprtn_lbind)
    zprtn_ub = zglbl_sym(zprtn_ubind)
    zprtn_lb(1) = MIN(zprtn_lb(1), z_xtnt(1))
    zprtn_ub(n_prtn(3)) = MIN(zprtn_ub(n_prtn(3)), z_xtnt(2))

    IF (ABS(z_xtnt(1)-z_xtnt(2)) .LE. small_value) THEN
! the rectangle is parallel to the xy plane
       ALLOCATE(tiles(n_prtn(1)*n_prtn(2)))
       k = 1
       DO i = 1, n_prtn(1)
          DO j = 1, n_prtn(2)
             A = (/ xprtn_lb(i), yprtn_ub(j), z_xtnt(1) /)
             B = (/ xprtn_lb(i), yprtn_lb(j), z_xtnt(1) /)
             C = (/ xprtn_ub(i), yprtn_lb(j), z_xtnt(1) /)
             D = (/ xprtn_ub(i), yprtn_ub(j), z_xtnt(1) /)

             ALLOCATE(tiles(k)%tile)
             tiles(k)%tile%tile_id = k
             ALLOCATE(tiles(k)%tile%rectangle)
             CALL cs_rectangle_create(tiles(k)%tile%rectangle, A, B, C, D, error)

             tiles(k)%tile%npts = xprtn_npts(i) * yprtn_npts(j)
             k = k + 1
          END DO
       END DO
    ELSE IF (ABS(y_xtnt(1)-y_xtnt(2)) .LE. small_value) THEN
! the rectangle is parallel to the xz plane
       ALLOCATE(tiles(n_prtn(1)*n_prtn(3)))
       k = 1
       DO i = 1, n_prtn(1)
          DO j = 1, n_prtn(3)
             A = (/ xprtn_lb(i), y_xtnt(1), zprtn_ub(j) /)
             B = (/ xprtn_lb(i), y_xtnt(1), zprtn_lb(j) /)
             C = (/ xprtn_ub(i), y_xtnt(1), zprtn_lb(j) /)
             D = (/ xprtn_ub(i), y_xtnt(1), zprtn_ub(j) /)

             ALLOCATE(tiles(k)%tile)
             tiles(k)%tile%tile_id = k
             ALLOCATE(tiles(k)%tile%rectangle)
             CALL cs_rectangle_create(tiles(k)%tile%rectangle, A, B, C, D, error)

             tiles(k)%tile%npts = xprtn_npts(i) * zprtn_npts(j)
             k = k + 1
          END DO
       END DO
    ELSE IF (ABS(x_xtnt(1)-x_xtnt(2)) .LE. small_value) THEN
! the rectangle is parallel to the yz plane
       ALLOCATE(tiles(n_prtn(2)*n_prtn(3)))
       k = 1
       DO i = 1, n_prtn(2)
          DO j = 1, n_prtn(3)
             A = (/ x_xtnt(1), yprtn_lb(i), zprtn_ub(j) /)
             B = (/ x_xtnt(1), yprtn_lb(i), zprtn_lb(j) /)
             C = (/ x_xtnt(1), yprtn_ub(i), zprtn_lb(j) /)
             D = (/ x_xtnt(1), yprtn_ub(i), zprtn_ub(j) /)

             ALLOCATE(tiles(k)%tile)
             tiles(k)%tile%tile_id = k
             ALLOCATE(tiles(k)%tile%rectangle)
             CALL cs_rectangle_create(tiles(k)%tile%rectangle, A, B, C, D, error)

             tiles(k)%tile%npts = yprtn_npts(i) * zprtn_npts(j)
             k = k + 1
          END DO
       END DO
    END IF

    DEALLOCATE(xglbl_sym, yglbl_sym, zglbl_sym)

    CALL timestop(handle)

  END SUBROUTINE partition_aa_rectangle_into_tiles

! *****************************************************************************
!> \brief    Partitions an axis-aligned cuboid into tiles. The cuboid is either
!>    a gate or a tile itself.
!> \param box axis-aligned cuboid to be partitioned
!> \param x_glbl x grid vetor of the simulation box
!> \param y_glbl y grid vetor of the simulation box
!> \param z_glbl z grid vetor of the simulation box
!> \param n_prtn vetor of size 3 specifying the number of times that the x, y and
!>               z interval (defining the region) should be partitioned into
!> \param tiles the obtained tiles    
!> \param error cp2k error
!> \par History
!>       12.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE partition_aa_cuboid_into_tiles(box, x_glbl, y_glbl, z_glbl, n_prtn, tiles, error)

    TYPE(cs_box), INTENT(IN), POINTER        :: box
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl
    INTEGER, DIMENSION(3), INTENT(IN)        :: n_prtn
    TYPE(tile_p_type), DIMENSION(:), &
      INTENT(OUT), POINTER                   :: tiles
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'partition_aa_cuboid_into_tiles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ii, jj, k, kk, &
                                                xxtnt_fst, xxtnt_ult, &
                                                yxtnt_fst, yxtnt_ult, &
                                                zxtnt_fst, zxtnt_ult
    INTEGER, DIMENSION(n_prtn(1))            :: xprtn_lbind, xprtn_npts, &
                                                xprtn_ubind
    INTEGER, DIMENSION(n_prtn(2))            :: yprtn_lbind, yprtn_npts, &
                                                yprtn_ubind
    INTEGER, DIMENSION(n_prtn(3))            :: zprtn_lbind, zprtn_npts, &
                                                zprtn_ubind
    REAL(dp)                                 :: step
    REAL(dp), DIMENSION(2)                   :: x_xtnt, y_xtnt, z_xtnt
    REAL(dp), DIMENSION(3, 8)                :: tile_vtx
    REAL(dp), DIMENSION(n_prtn(1))           :: xprtn_lb, xprtn_ub
    REAL(dp), DIMENSION(n_prtn(2))           :: yprtn_lb, yprtn_ub
    REAL(dp), DIMENSION(n_prtn(3))           :: zprtn_lb, zprtn_ub

    CALL timeset(routineN,handle)

! find the extents of the cuboid in all three directions
    x_xtnt(1) = MINVAL(box%vertices(1,:));    x_xtnt(2) = MAXVAL(box%vertices(1,:))
    y_xtnt(1) = MINVAL(box%vertices(2,:));    y_xtnt(2) = MAXVAL(box%vertices(2,:))
    z_xtnt(1) = MINVAL(box%vertices(3,:));    z_xtnt(2) = MAXVAL(box%vertices(3,:))

! find the position of the both ends of xtnt in the global axes
    xxtnt_fst = find_fst_ge(x_glbl, x_xtnt(1))
    xxtnt_ult = find_lst_le(x_glbl, x_xtnt(2))
    yxtnt_fst = find_fst_ge(y_glbl, y_xtnt(1))
    yxtnt_ult = find_lst_le(y_glbl, y_xtnt(2))
    zxtnt_fst = find_fst_ge(z_glbl, z_xtnt(1))
    zxtnt_ult = find_lst_le(z_glbl, z_xtnt(2))

! devide the integer sequence of indices {xtnt_fst .. xtnt_ult} into n_prtn partitions
    step = REAL((xxtnt_ult-xxtnt_fst),kind=dp)/REAL(n_prtn(1), kind=dp)
    xprtn_lbind = xxtnt_fst + NINT( (/(i, i=0, n_prtn(1)-1)/) * step ) + 1   ! lower bounds
    xprtn_ubind = xxtnt_fst + NINT( ((/(i, i=0, n_prtn(1)-1)/) + 1) * step ) ! upper bounds
    IF (step .EQ. 0.0_dp) xprtn_lbind = xprtn_ubind
! count the number of grid points
    xprtn_npts = xprtn_ubind - xprtn_lbind + 1
! find the corresponding values on the axes
    xprtn_lb = x_glbl(xprtn_lbind)
    xprtn_ub = x_glbl(xprtn_ubind)
! make sure that the tiles will cover the whole box
    xprtn_lb(1) = MIN(xprtn_lb(1), x_xtnt(1))
    xprtn_ub(n_prtn(1)) = MIN(xprtn_ub(n_prtn(1)), x_xtnt(2))

    step = REAL((yxtnt_ult-yxtnt_fst),kind=dp)/REAL(n_prtn(2), kind=dp)
    yprtn_lbind = yxtnt_fst + NINT( (/(i, i=0, n_prtn(2)-1)/) * step ) + 1
    yprtn_ubind = yxtnt_fst + NINT( ((/(i, i=0, n_prtn(2)-1)/) + 1) * step )
    IF (step .EQ. 0.0_dp) yprtn_lbind = yprtn_ubind
    yprtn_npts = yprtn_ubind - yprtn_lbind + 1
    yprtn_lb = y_glbl(yprtn_lbind)
    yprtn_ub = y_glbl(yprtn_ubind)
    yprtn_lb(1) = MIN(yprtn_lb(1), y_xtnt(1))
    yprtn_ub(n_prtn(2)) = MIN(yprtn_ub(n_prtn(2)), y_xtnt(2))

    step = REAL((zxtnt_ult-zxtnt_fst),kind=dp)/REAL(n_prtn(3), kind=dp)
    zprtn_lbind = zxtnt_fst + NINT( (/(i, i=0, n_prtn(3)-1)/) * step ) + 1
    zprtn_ubind = zxtnt_fst + NINT( ((/(i, i=0, n_prtn(3)-1)/) + 1) * step )
    IF (step .EQ. 0.0_dp) zprtn_lbind = zprtn_ubind
    zprtn_npts = zprtn_ubind - zprtn_lbind + 1
    zprtn_lb = z_glbl(zprtn_lbind)
    zprtn_ub = z_glbl(zprtn_ubind)
    zprtn_lb(1) = MIN(zprtn_lb(1), z_xtnt(1))
    zprtn_ub(n_prtn(3)) = MIN(zprtn_ub(n_prtn(3)), z_xtnt(2))

    ALLOCATE(tiles(n_prtn(1)*n_prtn(2)*n_prtn(3)))
    k = 1
    DO kk = 1, n_prtn(3)
       DO jj = 1, n_prtn(2)
          DO ii = 1, n_prtn(1)
             tile_vtx(1:3,1) = (/xprtn_lb(ii), yprtn_lb(jj), zprtn_lb(kk)/)
             tile_vtx(1:3,2) = (/xprtn_ub(ii), yprtn_lb(jj), zprtn_lb(kk)/)
             tile_vtx(1:3,3) = (/xprtn_ub(ii), yprtn_ub(jj), zprtn_lb(kk)/)
             tile_vtx(1:3,4) = (/xprtn_lb(ii), yprtn_ub(jj), zprtn_lb(kk)/)
             tile_vtx(1:3,5) = (/xprtn_lb(ii), yprtn_ub(jj), zprtn_ub(kk)/)
             tile_vtx(1:3,6) = (/xprtn_lb(ii), yprtn_lb(jj), zprtn_ub(kk)/)
             tile_vtx(1:3,7) = (/xprtn_ub(ii), yprtn_lb(jj), zprtn_ub(kk)/)
             tile_vtx(1:3,8) = (/xprtn_ub(ii), yprtn_ub(jj), zprtn_ub(kk)/)

             ALLOCATE(tiles(k)%tile)
             tiles(k)%tile%tile_id = k
             ALLOCATE(tiles(k)%tile%box)
             CALL cs_box_create(tiles(k)%tile%box, tile_vtx, error)

             tiles(k)%tile%npts = xprtn_npts(ii)*yprtn_npts(jj)*zprtn_npts(kk)
             k = k + 1
          END DO
       END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE partition_aa_cuboid_into_tiles

! *****************************************************************************
!> \brief   Returns the intersection of two coplanar segments
!> \param segment1  segment 1
!> \param segment2  segment 2
!> \param plane the plane that includes both segment 1 and segment 2
!> \param intersection the intersection of the two segments
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE coplanar_segment_segment_intersection(segment1, segment2, plane, intersection, error)

    TYPE(cs_segment), POINTER, INTENT(IN)             :: segment1, segment2
    TYPE(cs_plane), POINTER, INTENT(IN)               :: plane
    TYPE(cs_intersection_obj), POINTER, INTENT(INOUT) :: intersection
    TYPE(cp_error_type), INTENT(INOUT)                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'coplanar_segment_segment_intersection', &
      routineP = moduleN//':'//routineN

    LOGICAL                       :: intersect, cond1, cond2
    REAL(dp), DIMENSION(3)        :: A, B, C, D, unit_normal, AB ,AC, CD
    REAL(dp)                      :: tnum, tden
    INTEGER                       :: handle

    CALL timeset(routineN,handle)

    A = segment1%A;    B = segment1%B
    C = segment2%A;    D = segment2%B
    unit_normal = plane%unit_normal

! first check if the line segments intersect (4 counter-clockwise orientation test)
! does NOT handle the case where the line segments are colinear
! here the case is already taken care of using is_inside_polygon and is_inside_rectangle
    cond1 = ordered_counter_clockwise(A, C, D, unit_normal) .NEQV. &
            ordered_counter_clockwise(B, C, D, unit_normal)
    cond2 = ordered_counter_clockwise(A, B, C, unit_normal) .NEQV. &
            ordered_counter_clockwise(A, B, D, unit_normal)
    IF (cond1 .AND. cond2) THEN
       intersect = .TRUE.
    ELSE
       intersect = .FALSE.
    END IF

! if they inersect find the intersection point
    IF (intersect .EQV. .TRUE.) THEN
       AB = segment1%AB
       CD = segment2%AB
       AC = C - A
       tnum = SQRT(SUM(vector_product(AC, CD)**2))
       tden = SQRT(SUM(vector_product(AB, CD)**2))
! I might need to check if tden is non-zero or is less than small_val ????
       CALL cs_intersection_obj_create(intersection, PT_INTXN, error=error)
       intersection%point = A + AB * tnum / tden
    ELSE
       CALL cs_intersection_obj_create(intersection, NO_INTXN, error= error)
    END IF

    CALL timestop(handle)

    CONTAINS
! *****************************************************************************
!> \brief checks if three points A, B, C are ordered counter-clockwise according
!>        to the normal vector to the plane that they are located on
!> \param A ...
!> \param B ...
!> \param C ...
!> \param unit_normal ...
!> \retval ordered_counter_clockwise ...
! *****************************************************************************
    PURE FUNCTION ordered_counter_clockwise(A, B, C, unit_normal)

    REAL(dp), DIMENSION(3), INTENT(IN)       :: A, B, C, unit_normal
    LOGICAL                                  :: ordered_counter_clockwise

    REAL(dp), DIMENSION(3)                   :: AB, ABxBC, BC

       ordered_counter_clockwise = .TRUE.
       AB = B - A
       BC = C - B
       ABxBC = vector_product(AB, BC)
       IF (DOT_PRODUCT(ABxBC, unit_normal) .GE. 0.0_dp) ordered_counter_clockwise = .FALSE.

    END FUNCTION ordered_counter_clockwise

  END SUBROUTINE coplanar_segment_segment_intersection

! *****************************************************************************
!> \brief   Returns the intersection of a plane and a segment
!> \param segment input segment
!> \param plane input plane
!> \param intersection the intersection of the plane and the segment
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE segment_plane_intersection(segment, plane, intersection, error)

    TYPE(cs_segment), INTENT(IN), POINTER    :: segment
    TYPE(cs_plane), INTENT(IN), POINTER      :: plane
    TYPE(cs_intersection_obj), &
      INTENT(INOUT), POINTER                 :: intersection
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'segment_plane_intersection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(dp)                                 :: ndotA, ndotAB, t

    CALL timeset(routineN,handle)

    ndotA = DOT_PRODUCT(plane%unit_normal, segment%A)
    ndotAB = DOT_PRODUCT(plane%unit_normal, segment%AB)

! if AB and the plane are parallel
    IF (ABS(ndotAB) .LE. small_value) THEN
! AB is entirely on the plane ... return A and B
       IF ((plane%d - ndotA) .EQ. 0.0_dp) THEN
          CALL cs_intersection_obj_create(intersection, SEG_INTXN, error=error)
          intersection%segment = segment
! 0 intersection
       ELSE
          CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
       END IF
! if AB and the plane are NOT parallel
    ELSE
       t = (plane%d - ndotA) / ndotAB
! AB intersects the plane
       IF (((0 .LE. t) .AND. (t .LE. 1)) .OR. &
           (ABS(t) .LE. small_value) .OR.     &
           (ABS(t-1) .LE. small_value)) THEN
          CALL cs_intersection_obj_create(intersection, PT_INTXN, error=error)
          intersection%point = segment%A + (t * segment%AB)
! 0 intersection
       ELSE
          CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE segment_plane_intersection

! *****************************************************************************
!> \brief   Returns the intersection of a segment and a rectangle
!> \param segment input segment
!> \param rectangle input rectangle
!> \param intersection the intersection of the segment and the rectangle
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE segment_rectangle_intersection(segment, rectangle, intersection, error)

    TYPE(cs_segment), INTENT(IN), POINTER    :: segment
    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    TYPE(cs_intersection_obj), &
      INTENT(INOUT), POINTER                 :: intersection
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'segment_rectangle_intersection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, &
                                                n_intersections, n_intxn_tot
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: unique_intxn, valid_intxn
    REAL(dp), DIMENSION(3)                   :: A, B, crossing_pt
    REAL(dp), DIMENSION(3, 3)                :: all_intxn
    TYPE(cs_edge), POINTER                   :: rect_edge
    TYPE(cs_intersection_obj), POINTER       :: seg_pln_intxn, seg_seg_intxn
    TYPE(cs_plane), POINTER                  :: plane

    CALL timeset(routineN,handle)

    ALLOCATE(seg_pln_intxn)
    plane => rectangle%plane

! check if the segment crosses the rectangle's plane
    CALL segment_plane_intersection(segment, plane, seg_pln_intxn, error)
    SELECT CASE (seg_pln_intxn%intxn_type)
! the segment does not cross the rectangle's plane
    CASE (NO_INTXN)
       CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
! the segment crosses the rectangle's plane at one point
    CASE (PT_INTXN)
       IF (is_inside_rectangle(seg_pln_intxn%point, rectangle)) THEN
          CALL cs_intersection_obj_create(intersection, PT_INTXN, error=error)
          intersection%point = seg_pln_intxn%point
       ELSE
          CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
       END IF
! the segment and the rectangle are coplanar
    CASE (SEG_INTXN)
       A = seg_pln_intxn%segment%A
       B = seg_pln_intxn%segment%B

       ! the both ends of the segment, A and B, are inside the rectangle
       IF (is_inside_rectangle(A, rectangle) .AND. is_inside_rectangle(B, rectangle)) THEN
          CALL cs_intersection_obj_create(intersection, SEG_INTXN, error=error)
          intersection%segment => seg_pln_intxn%segment

       ! A is inside rectangle and B is not inside
       ELSE IF (is_inside_rectangle(A, rectangle) .AND. .NOT. is_inside_rectangle(B, rectangle)) THEN
          ALLOCATE(seg_seg_intxn)
          crossing_pt = A
          DO i = 1, 4
             rect_edge => rectangle%edges(i)
             CALL coplanar_segment_segment_intersection(rect_edge%segment, segment, &
                                                        plane, seg_seg_intxn, error)
             IF (seg_seg_intxn%intxn_type .NE. NO_INTXN) THEN
                crossing_pt = seg_seg_intxn%point
                EXIT
             END IF
          END DO
          IF (are_equal(crossing_pt, A)) THEN
             CALL cs_intersection_obj_create(intersection, PT_INTXN, error=error)
             intersection%point = A
          ELSE
             CALL cs_intersection_obj_create(intersection, SEG_INTXN, error=error)
             CALL cs_segment_create(intersection%segment, A, crossing_pt, error=error)
          END IF
          CALL cs_intersection_obj_release(seg_seg_intxn, error)

       ! A is not inside rectangle and B is inside
       ELSE IF (.NOT. is_inside_rectangle(A, rectangle) .AND. is_inside_rectangle(B, rectangle)) THEN
          ALLOCATE(seg_seg_intxn)
          crossing_pt = B
          DO i = 1, 4
             rect_edge => rectangle%edges(i)
             CALL coplanar_segment_segment_intersection(rect_edge%segment, segment, &
                                                        plane, seg_seg_intxn, error)
             IF (seg_seg_intxn%intxn_type .NE. NO_INTXN) THEN
                crossing_pt = seg_seg_intxn%point
                EXIT
             END IF
          END DO
          IF (are_equal(crossing_pt, B)) THEN
             CALL cs_intersection_obj_create(intersection, PT_INTXN, error=error)
             intersection%point = B
          ELSE
             CALL cs_intersection_obj_create(intersection, SEG_INTXN, error=error)
             CALL cs_segment_create(intersection%segment, crossing_pt, B, error=error)
          END IF
          CALL cs_intersection_obj_release(seg_seg_intxn, error)

       ! A and B are not inside rectangle
       ELSE
       ! find the intersection between the segment and rectangle edges
          ALLOCATE(seg_seg_intxn)
          n_intxn_tot = 0
          j = 1
          DO i = 1, 4
             rect_edge => rectangle%edges(i)
             CALL coplanar_segment_segment_intersection(rect_edge%segment, segment, &
                                                        plane, seg_seg_intxn, error)
             n_intxn_tot = n_intxn_tot + seg_seg_intxn%n_intxn
             all_intxn(:,j) = seg_seg_intxn%point
             j = j + seg_seg_intxn%n_intxn
          END DO

          ALLOCATE(valid_intxn(3,n_intxn_tot))
          valid_intxn(:,:) = all_intxn(:,1:n_intxn_tot)

          IF (n_intxn_tot .NE. 0) THEN
          ! if there is any intersection remove duplicate points
             CALL remove_duplicate_vertices(valid_intxn, unique_intxn, error)
             n_intersections = SIZE(unique_intxn)/3

             SELECT CASE (n_intersections)
             CASE (1)
                CALL cs_intersection_obj_create(intersection, PT_INTXN, error=error)
                intersection%point = unique_intxn(:,1)
             CASE (2)
                CALL cs_intersection_obj_create(intersection, SEG_INTXN, error=error)
                CALL cs_segment_create(intersection%segment, &
                                       unique_intxn(:,1), unique_intxn(:,2), &
                                       error)
             CASE DEFAULT
             END SELECT
          ELSE
             CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
          END IF

          DEALLOCATE(valid_intxn)
          CALL cs_intersection_obj_release(seg_seg_intxn, error)
       END IF
    END SELECT ! status of the segment wrt the rectangle's plane

    CALL cs_intersection_obj_release(seg_pln_intxn, error)
    CALL timestop(handle)

  END SUBROUTINE segment_rectangle_intersection

! *****************************************************************************
!> \brief Returns the intersection of a plane and a box in 3D space.
!>        The output is either a point, triangle, polygon (including rectangle)
!>        or infinity (no intersection)
!> \param plane input plane
!> \param box input box
!> \param intersection the intersection of the plane and the box
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE plane_box_intersection(plane, box, intersection, error)

    TYPE(cs_plane), INTENT(IN), POINTER      :: plane
    TYPE(cs_box), INTENT(IN), POINTER        :: box
    TYPE(cs_intersection_obj), &
      INTENT(INOUT), POINTER                 :: intersection
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'plane_box_intersection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, &
                                                n_intersections, n_intxn_tot
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: unique_intxn, valid_intxn
    REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
      TARGET                                 :: intersection_points
    REAL(dp), DIMENSION(3)                   :: O, unit_normal
    REAL(dp), DIMENSION(3, 24)               :: all_intxn
    REAL(dp), DIMENSION(:, :), POINTER       :: polygon_vertices
    TYPE(cs_intersection_obj), POINTER       :: seg_pln_intxn

    CALL timeset(routineN,handle)

    ALLOCATE(seg_pln_intxn)
! find the intersection of the edges of the box with the plane
    n_intxn_tot = 0
    j = 1
    DO i = 1, 12
       CALL segment_plane_intersection(box%edges(i)%segment, plane, &
                                       seg_pln_intxn, error)
       n_intxn_tot = n_intxn_tot + seg_pln_intxn%n_intxn
       SELECT CASE (seg_pln_intxn%intxn_type)
       CASE (PT_INTXN)
          all_intxn(:,j) = seg_pln_intxn%point
       CASE (SEG_INTXN)
          all_intxn(:,j) = seg_pln_intxn%segment%A
          all_intxn(:,j+1) = seg_pln_intxn%segment%B
          CALL cs_segment_release(seg_pln_intxn%segment, error)
       END SELECT
       j = j + seg_pln_intxn%n_intxn
    END DO
    CALL cs_intersection_obj_release(seg_pln_intxn, error)

! get actual intersection points i.e. remove infinity points
    ALLOCATE(valid_intxn(3,n_intxn_tot))
    valid_intxn(:,:) = all_intxn(:,1:n_intxn_tot)

    IF (n_intxn_tot .NE. 0) THEN
! if there is any intersection remove duplicate points
       CALL remove_duplicate_vertices(valid_intxn, unique_intxn, error)
       n_intersections = SIZE(unique_intxn)/3

       SELECT CASE (n_intersections)
       CASE (1)
          CALL cs_intersection_obj_create(intersection, PT_INTXN, error=error)
          intersection%point = unique_intxn(:,1)
       CASE (2)
          CALL cs_intersection_obj_create(intersection, SEG_INTXN, error=error)
          CALL cs_segment_create(intersection%segment, &
                                 unique_intxn(:,1), unique_intxn(:,2), &
                                 error)
       CASE (3)
          CALL cs_intersection_obj_create(intersection, TRI_INTXN, error=error)
          CALL cs_triangle_create(intersection%triangle, &
                                 unique_intxn(:,1), unique_intxn(:,2), unique_intxn(:,3), &
                                 error)
     ! CASE (4)
          ! consider rectangle as a regular polygon
       CASE DEFAULT
          CALL cs_intersection_obj_create(intersection, PY_INTXN, n_intersections, error=error)

          ALLOCATE(intersection_points(3,n_intersections))
          intersection_points(:,:) = unique_intxn(:,1:n_intersections)

          ! sort intersection points
          O = intersection_points(:,1)     ! reference vertex
          unit_normal = plane%unit_normal
          IF (n_intersections .GT. 3) CALL sort_vertices(intersection_points, O, unit_normal, error)
          polygon_vertices => intersection_points
          CALL cs_polygon_create(intersection%polygon, polygon_vertices, error=error)

          DEALLOCATE(intersection_points)
       END SELECT

    ELSE
! if there is no intersection
       CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
    END IF

    DEALLOCATE(valid_intxn)

    CALL timestop(handle)

  END SUBROUTINE plane_box_intersection

! *****************************************************************************
!> \brief Returns the overlap of a polygon and a rectangle that are located on
!>        the same plane. The output is either a point, segment, triangle, polygon
!>        (including rectangle) or infinity (no intersection)
!> \param polygon input polygon
!> \param rectangle input rectangle
!> \param intersection the intersection of the polygon and the rectangle
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE coplanar_polygon_rectangle_intersection(polygon, rectangle, intersection, error)
    TYPE(cs_polygon), INTENT(IN), POINTER    :: polygon
    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    TYPE(cs_intersection_obj), &
      INTENT(INOUT), POINTER                 :: intersection
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'coplanar_polygon_rectangle_intersection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, handle, i, j, k, &
                                                ovlp_n_vtx, py_n_vtx
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: all_pts, crossing_pts, &
                                                internal_pts, &
                                                unique_vertices, &
                                                valid_vertices
    REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
      TARGET                                 :: overlap_vertices
    REAL(dp), DIMENSION(3)                   :: O, unit_normal, vertex
    REAL(dp), DIMENSION(:, :), POINTER       :: ovlp_vertices
    TYPE(cs_edge), POINTER                   :: py_edge, rect_edge
    TYPE(cs_intersection_obj), POINTER       :: seg_seg_intxn
    TYPE(cs_plane), POINTER                  :: plane

    CALL timeset(routineN,handle)

    ALLOCATE(seg_seg_intxn)
    plane => polygon%plane
    py_n_vtx = polygon%n_vertices

! find internal vertices
    ALLOCATE(internal_pts(3, 4 + py_n_vtx))
    ALLOCATE(crossing_pts(3, 4 * py_n_vtx))
    internal_pts = huge_value
    crossing_pts = 0.0_dp


    counter = 0              ! count valid vertices
    DO i = 1, 4              ! loop over rectangle's vertices
       vertex = rectangle%vertices(:,i)
       IF (is_inside_polygon(vertex, polygon)) THEN
          internal_pts(:, i) = vertex
          counter = counter + 1
       END IF
    END DO

    DO i = 1, py_n_vtx       ! loop over polygon's vertices
       vertex = polygon%vertices(:,i)
       IF (is_inside_rectangle(vertex, rectangle)) THEN
          internal_pts(:, i+4) = vertex
          counter = counter + 1
       END IF
    END DO

! find crossing points between the rectangle and the polygon edges
    k = 1
    DO i = 1, py_n_vtx       ! loop over polygon's edges
       DO j = 1, 4           ! lopp over rectangle's edges
          py_edge => polygon%edges(i)
          rect_edge => rectangle%edges(j)
          CALL coplanar_segment_segment_intersection(py_edge%segment, rect_edge%segment, &
                                                     plane, seg_seg_intxn, error)
          crossing_pts(:,k) = seg_seg_intxn%point
          counter = counter + seg_seg_intxn%n_intxn
          k = k + 1
       END DO
     END DO

! concatenate internal_pts and crossing_pts
    ALLOCATE(all_pts(3, (4+py_n_vtx) + 4*py_n_vtx))
    all_pts( : , 1:(4+py_n_vtx) ) = internal_pts
    all_pts( : , (4+py_n_vtx)+1:(4+py_n_vtx) + 4*py_n_vtx ) = crossing_pts

! remove infinity points
    ALLOCATE(valid_vertices(3,counter))
    j = 1
    DO i = 1, (4+py_n_vtx) + 4*py_n_vtx
       vertex = all_pts(:,i)
       IF (.NOT. is_infinity_point(vertex)) THEN
          valid_vertices(:,j) = vertex
          j = j + 1
       END IF
    END DO

    IF (counter .NE. 0) THEN
! if there is any intersection remove duplicate points
       CALL remove_duplicate_vertices(valid_vertices, unique_vertices, error)
       ovlp_n_vtx = SIZE(unique_vertices)/3

       SELECT CASE (ovlp_n_vtx)
       CASE (1)
          CALL cs_intersection_obj_create(intersection, PT_INTXN, error=error)
          intersection%point = unique_vertices(:,1)
       CASE (2)
          CALL cs_intersection_obj_create(intersection, SEG_INTXN, error=error)
          CALL cs_segment_create(intersection%segment, &
                                 unique_vertices(:,1), unique_vertices(:,2), &
                                 error)
       CASE (3)
          CALL cs_intersection_obj_create(intersection, TRI_INTXN, error=error)
          CALL cs_triangle_create(intersection%triangle, &
                                 unique_vertices(:,1), unique_vertices(:,2), unique_vertices(:,3), &
                                 error)
     ! CASE (4)
          ! consider rectangle as a regular polygon
       CASE DEFAULT
          CALL cs_intersection_obj_create(intersection, PY_INTXN, ovlp_n_vtx, error=error)

          ALLOCATE(overlap_vertices(3,ovlp_n_vtx))
          overlap_vertices(:,:) = unique_vertices(:,1:ovlp_n_vtx)

          ! sort intersection points
          O = overlap_vertices(:,1)
          unit_normal = plane%unit_normal
          IF (ovlp_n_vtx .GT. 3) CALL sort_vertices(overlap_vertices, O, unit_normal, error)
          ovlp_vertices => overlap_vertices
          CALL cs_polygon_create(intersection%polygon, ovlp_vertices, error=error)

          DEALLOCATE(overlap_vertices)
       END SELECT

    ELSE
! if there is no intersection
       CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
    END IF

    DEALLOCATE(internal_pts, crossing_pts, all_pts, valid_vertices)
    CALL cs_intersection_obj_release(seg_seg_intxn, error)

    CALL timestop(handle)

  END SUBROUTINE coplanar_polygon_rectangle_intersection

! *****************************************************************************
!> \brief   finds the intersection between two rectangles
!> \param rectangle1 rectangle 1
!> \param rectangle2 rectangle 2
!> \param intersection the intersection between rectangle 1 and rectangle 2
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE rectangle_rectangle_intersection(rectangle1, rectangle2, intersection, error)

    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle1, rectangle2
    TYPE(cs_intersection_obj), &
      INTENT(INOUT), POINTER                 :: intersection
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'rectangle_rectangle_intersection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, handle, i, j, &
                                                n_intersections, n_intxn_tot
    LOGICAL                                  :: xranges_overlap, &
                                                yranges_overlap, &
                                                zranges_overlap
    LOGICAL, DIMENSION(3)                    :: mask
    REAL(dp)                                 :: maxx1, maxx2, maxy1, maxy2, &
                                                maxz1, maxz2, minx1, minx2, &
                                                miny1, miny2, minz1, minz2, &
                                                n1dotA1A2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: unique_intxn, valid_intxn
    REAL(dp), DIMENSION(3)                   :: A1A2, n1, n2
    REAL(dp), DIMENSION(3, 4)                :: all_intxn
    TYPE(cs_intersection_obj), POINTER       :: seg_rect_intxn
    TYPE(cs_polygon), POINTER                :: polygon_tmp

    CALL timeset(routineN,handle)

    n1 = rectangle1%plane%unit_normal
    n2 = rectangle2%plane%unit_normal

! rectangles are parallel
    IF (SQRT(SUM(vector_product(n1, n2)**2)) .LE. small_value) THEN
       A1A2 = rectangle1%plane%A - rectangle2%plane%A
       n1dotA1A2 = DOT_PRODUCT(rectangle1%plane%unit_normal, A1A2)
! rectangles are co-planar
       IF (n1dotA1A2 .EQ. 0.0_dp) THEN
          ALLOCATE(polygon_tmp)
          CALL cast_rectangle_to_polygon(rectangle1, polygon_tmp, error)
          CALL coplanar_polygon_rectangle_intersection(polygon_tmp, rectangle2, &
                                                       intersection, error)
          CALL cs_polygon_release(polygon_tmp, error)
! rectangles are parallel but not co-planar
       ELSE
          CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
       END IF
! rectangles are not parallel
    ELSE
       ! check the extents of the rectangles
       minx1 = MINVAL(rectangle1%vertices(1,:));    maxx1 = MAXVAL(rectangle1%vertices(1,:))
       miny1 = MINVAL(rectangle1%vertices(2,:));    maxy1 = MAXVAL(rectangle1%vertices(2,:))
       minz1 = MINVAL(rectangle1%vertices(3,:));    maxz1 = MAXVAL(rectangle1%vertices(3,:))
       minx2 = MINVAL(rectangle2%vertices(1,:));    maxx2 = MAXVAL(rectangle2%vertices(1,:))
       miny2 = MINVAL(rectangle2%vertices(2,:));    maxy2 = MAXVAL(rectangle2%vertices(2,:))
       minz2 = MINVAL(rectangle2%vertices(3,:));    maxz2 = MAXVAL(rectangle2%vertices(3,:))

       xranges_overlap = MAXVAL((/minx1, minx2/)) .LE. MINVAL((/maxx1, maxx2/))
       yranges_overlap = MAXVAL((/miny1, miny2/)) .LE. MINVAL((/maxy1, maxy2/))
       zranges_overlap = MAXVAL((/minz1, minz2/)) .LE. MINVAL((/maxz1, maxz2/))

       mask = .TRUE.
       counter = COUNT((/xranges_overlap, yranges_overlap, zranges_overlap/) .EQV. mask)
! rectangles might intersect
       IF (counter .GT. 2) THEN
       ! find the intersection of the edges of rectangle1 with rectangle2
          all_intxn = huge_value
          ALLOCATE(seg_rect_intxn)
          n_intxn_tot = 0
          j = 1
          DO i = 1, 4
             CALL segment_rectangle_intersection(rectangle1%edges(i)%segment, rectangle2, &
                                                 seg_rect_intxn, error)
             n_intxn_tot = n_intxn_tot + seg_rect_intxn%n_intxn
             SELECT CASE (seg_rect_intxn%intxn_type)
             CASE (PT_INTXN)
                all_intxn(:,j) = seg_rect_intxn%point
             CASE (SEG_INTXN)
                all_intxn(:,j) = seg_rect_intxn%segment%A
                all_intxn(:,j+1) = seg_rect_intxn%segment%B
                CALL cs_segment_release(seg_rect_intxn%segment, error)
             END SELECT
             j = j + seg_rect_intxn%n_intxn
          END DO
          CALL cs_intersection_obj_release(seg_rect_intxn, error)

          ! get actual intersection points i.e. remove infinity points
          ALLOCATE(valid_intxn(3,n_intxn_tot))
          valid_intxn(:,:) = all_intxn(:,1:n_intxn_tot)

          IF (n_intxn_tot .NE. 0) THEN
          ! if there is any intersection remove duplicate points
             CALL remove_duplicate_vertices(valid_intxn, unique_intxn, error)
             n_intersections = SIZE(unique_intxn)/3

             SELECT CASE (n_intersections)
             CASE (1)
                CALL cs_intersection_obj_create(intersection, PT_INTXN, error=error)
                intersection%point = unique_intxn(:,1)
             CASE (2)
                CALL cs_intersection_obj_create(intersection, SEG_INTXN, error=error)
                CALL cs_segment_create(intersection%segment, &
                                       unique_intxn(:,1), unique_intxn(:,2), &
                                       error)
             END SELECT

          ELSE
          ! if there is no intersection
             CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
          END IF
! rectangles are not parallel and do not intersect
       ELSE
          CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE rectangle_rectangle_intersection

! *****************************************************************************
!> \brief   Returns the intersection of a rectangle and a box
!> \param rectangle input rectangle
!> \param box input box
!> \param intersection the intersection of the rectangle and the box
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE rectangle_box_intersection(rectangle, box, intersection, error)

    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    TYPE(cs_box), INTENT(IN), POINTER        :: box
    TYPE(cs_intersection_obj), &
      INTENT(INOUT), POINTER                 :: intersection
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rectangle_box_intersection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cs_intersection_obj), POINTER       :: pln_box_intxn
    TYPE(cs_polygon), POINTER                :: polygon_tmp

    CALL timeset(routineN,handle)

    ALLOCATE(pln_box_intxn)

! first find the intersection of the rectangle's plane with the box
    CALL plane_box_intersection(rectangle%plane, box, pln_box_intxn, error)

! then find the overlap of the intersection with the rectangle
    SELECT CASE (pln_box_intxn%intxn_type)
    CASE (NO_INTXN)
       CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
    CASE (PT_INTXN)
       IF (is_inside_rectangle(pln_box_intxn%point, rectangle)) THEN
          CALL cs_intersection_obj_create(intersection, PT_INTXN, error=error)
          intersection%point = pln_box_intxn%point
       ELSE
          CALL cs_intersection_obj_create(intersection, NO_INTXN, error=error)
       END IF
    CASE (SEG_INTXN)
       CALL segment_rectangle_intersection(pln_box_intxn%segment, rectangle, &
                                           intersection, error)
    CASE (TRI_INTXN)
       ALLOCATE(polygon_tmp)
       CALL cast_triangle_to_polygon(pln_box_intxn%triangle, polygon_tmp, error)
       CALL coplanar_polygon_rectangle_intersection(polygon_tmp, rectangle, &
                                                    intersection, error)
       CALL cs_polygon_release(polygon_tmp, error)
    CASE DEFAULT
       CALL coplanar_polygon_rectangle_intersection(pln_box_intxn%polygon, rectangle, &
                                                    intersection, error)
    END SELECT

    CALL cs_intersection_obj_release(pln_box_intxn, error)

    CALL timestop(handle)

  END SUBROUTINE rectangle_box_intersection

! *****************************************************************************
!> \brief  toggles a grid_segment on if a given rectangle intersects with it
!> \param grid_seg grid segment whose end points to be weighted
!> \param rectangle the given rectangle
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE toggle_grid_segment(grid_seg, rectangle, error)

    TYPE(cs_grid_segment), INTENT(INOUT), &
      POINTER                                :: grid_seg
    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'toggle_grid_segment', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cs_intersection_obj), POINTER       :: seg_rect_intxn

    CALL timeset(routineN,handle)

    ALLOCATE(seg_rect_intxn)
    CALL segment_rectangle_intersection(grid_seg%segment, rectangle, seg_rect_intxn, error)

    IF (seg_rect_intxn%intxn_type .EQ. NO_INTXN) THEN
       grid_seg%is_on  = .FALSE.
       grid_seg%weight = 0
    ELSE
       grid_seg%is_on  = .TRUE.
       grid_seg%weight = 1
    END IF

    CALL cs_intersection_obj_release(seg_rect_intxn, error)
    CALL timestop(handle)

  END SUBROUTINE toggle_grid_segment

! *****************************************************************************
!> \brief  toggles a grid_tile on if a given rectangle intersects with it
!> \param grid_tile grid tile whose vertices to be weighted
!> \param rectangle input rectangle
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE toggle_grid_tile(grid_tile, rectangle, error)

    TYPE(cs_grid_tile), INTENT(INOUT), &
      POINTER                                :: grid_tile
    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'toggle_grid_tile', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cs_intersection_obj), POINTER       :: rect_rect_intxn

    CALL timeset(routineN,handle)

    ALLOCATE(rect_rect_intxn)
    CALL rectangle_rectangle_intersection(rectangle, grid_tile%rectangle, rect_rect_intxn, error)

    IF (rect_rect_intxn%intxn_type .EQ. NO_INTXN) THEN
       grid_tile%is_on  = .FALSE.
       grid_tile%weight = 0
    ELSE
       grid_tile%is_on  = .TRUE.
       grid_tile%weight = 1
    END IF

    CALL cs_intersection_obj_release(rect_rect_intxn, error)

    CALL timestop(handle)

  END SUBROUTINE toggle_grid_tile

! *****************************************************************************
!> \brief  toggles a voxel on if a given rectangle intersects with it
!> \param voxel voxel whose vertices to be weighted
!> \param rectangle the given rectangle
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE toggle_voxel(voxel, rectangle, error)

    TYPE(cs_voxel), INTENT(INOUT), POINTER   :: voxel
    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'toggle_voxel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cs_intersection_obj), POINTER       :: rect_box_intxn

    CALL timeset(routineN,handle)
    ALLOCATE(rect_box_intxn)
    CALL rectangle_box_intersection(rectangle, voxel%box, rect_box_intxn, error)

    IF (rect_box_intxn%intxn_type .EQ. NO_INTXN) THEN
       voxel%is_on  = .FALSE.
       voxel%weight = 0
    ELSE
       voxel%is_on  = .TRUE.
       voxel%weight = 1
    END IF

    CALL cs_intersection_obj_release(rect_box_intxn, error)
    CALL timestop(handle)

  END SUBROUTINE toggle_voxel

! *****************************************************************************
!> \brief Returns the type of the local domain
!> \param bounds_local current process's local bounds
!> \retval ld_type local domain type
!> \par History
!>       12.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  FUNCTION local_domain_type(bounds_local) RESULT(ld_type)

    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bounds_local
    INTEGER                                  :: ld_type

    CHARACTER(LEN=*), PARAMETER :: routineN = 'local_domain_type', &
      routineP = moduleN//':'//routineN
    LOGICAL, DIMENSION(3), PARAMETER         :: test_eq_bounds = .TRUE.

    INTEGER                                  :: lb1, lb2, lb3, n_eq_bounds, &
                                                ub1, ub2, ub3

    lb1 = bounds_local(1,1);    ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);    ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);    ub3 = bounds_local(2,3)

    n_eq_bounds = COUNT((/lb1 .EQ. ub1,lb2 .EQ. ub2,lb3 .EQ. ub3/) .EQV. test_eq_bounds)

    SELECT CASE(n_eq_bounds)
       CASE (0) 
          ld_type = axis_aligned_box
       CASE (1) 
          IF (lb1 .EQ. ub1) THEN       
             ld_type = yz_aligned_rectangle
          ELSE IF (lb2 .EQ. ub2) THEN 
             ld_type = xz_aligned_rectangle
          ELSE IF (lb3 .EQ. ub3) THEN 
             ld_type = xy_aligned_rectangle
          END IF 
       CASE (2) 
          IF ((lb1 .EQ. ub1) .AND. (lb2 .EQ. ub2)) THEN      
             ld_type = zaxis_aligned_segment
          ELSE IF ((lb1 .EQ. ub1) .AND. (lb3 .EQ. ub3)) THEN 
             ld_type = yaxis_aligned_segment
          ELSE IF ((lb2 .EQ. ub2) .AND. (lb3 .EQ. ub3)) THEN 
             ld_type = xaxis_aligned_segment
          END IF 
       CASE (3) 
          ld_type = single_grid_point
       CASE (4) 
          ld_type = no_grid_point
    END SELECT

  END FUNCTION local_domain_type

! *****************************************************************************
!> \brief  maps the weights of vertices of a grid segment on a pw data type
!> \param grid_seg input grid segment
!> \param pw the output plane wave data type
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE convert_grid_segment_to_grid(grid_seg, pw, error)

    TYPE(cs_grid_segment), INTENT(IN), &
      POINTER                                :: grid_seg
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_grid_segment_to_grid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idx_x, idx_y, idx_z
    REAL(dp), DIMENSION(2)                   :: wt

    CALL timeset(routineN,handle)

    IF (grid_seg%is_on) THEN
       idx_x = grid_seg%grid_segment_idx(1)   
       idx_y = grid_seg%grid_segment_idx(2)
       idx_z = grid_seg%grid_segment_idx(3)

       wt = REAL(grid_seg%weight, kind=dp)

       SELECT CASE (grid_seg%alignment)
         CASE (xaxis_aligned_segment)
            pw%cr3d(idx_x  , idx_y  , idx_z  ) = pw%cr3d(idx_x  , idx_y  , idx_z  ) + wt(1)
            pw%cr3d(idx_x+1, idx_y  , idx_z  ) = pw%cr3d(idx_x+1, idx_y  , idx_z  ) + wt(2)
         CASE (yaxis_aligned_segment)
            pw%cr3d(idx_x  , idx_y  , idx_z  ) = pw%cr3d(idx_x  , idx_y  , idx_z  ) + wt(1)
            pw%cr3d(idx_x  , idx_y+1, idx_z  ) = pw%cr3d(idx_x  , idx_y+1, idx_z  ) + wt(2)
         CASE (zaxis_aligned_segment)
            pw%cr3d(idx_x  , idx_y  , idx_z  ) = pw%cr3d(idx_x  , idx_y  , idx_z  ) + wt(1)
            pw%cr3d(idx_x  , idx_y  , idx_z+1) = pw%cr3d(idx_x  , idx_y  , idx_z+1) + wt(2)
       END SELECT
    ELSE
       ! do nothing
    END IF
    
    CALL timestop(handle)

  END SUBROUTINE convert_grid_segment_to_grid

! *****************************************************************************
!> \brief  maps the weights of vertices of a grid tile on a pw data type
!> \param grid_tile input grid tile
!> \param pw the output plane wave data type
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE convert_grid_tile_to_grid(grid_tile, pw, error)

    TYPE(cs_grid_tile), INTENT(IN), POINTER  :: grid_tile
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_grid_tile_to_grid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idx_x, idx_y, idx_z
    REAL(dp), DIMENSION(4)                   :: wt

    CALL timeset(routineN,handle)

    IF (grid_tile%is_on) THEN
       idx_x = grid_tile%grid_tile_idx(1)   
       idx_y = grid_tile%grid_tile_idx(2)
       idx_z = grid_tile%grid_tile_idx(3)

       wt = REAL(grid_tile%weight, kind=dp)

       SELECT CASE (grid_tile%alignment)
         CASE (xy_aligned_rectangle)
            pw%cr3d(idx_x  , idx_y  , idx_z  ) = pw%cr3d(idx_x  , idx_y  , idx_z  ) + wt(1)
            pw%cr3d(idx_x+1, idx_y  , idx_z  ) = pw%cr3d(idx_x+1, idx_y  , idx_z  ) + wt(2)
            pw%cr3d(idx_x+1, idx_y+1, idx_z  ) = pw%cr3d(idx_x+1, idx_y+1, idx_z  ) + wt(3)
            pw%cr3d(idx_x  , idx_y+1, idx_z  ) = pw%cr3d(idx_x  , idx_y+1, idx_z  ) + wt(4)
         CASE (yz_aligned_rectangle)
            pw%cr3d(idx_x  , idx_y  , idx_z  ) = pw%cr3d(idx_x  , idx_y  , idx_z  ) + wt(1)
            pw%cr3d(idx_x  , idx_y+1, idx_z  ) = pw%cr3d(idx_x  , idx_y+1, idx_z  ) + wt(2)
            pw%cr3d(idx_x  , idx_y+1, idx_z+1) = pw%cr3d(idx_x  , idx_y+1, idx_z+1) + wt(3)
            pw%cr3d(idx_x  , idx_y  , idx_z+1) = pw%cr3d(idx_x  , idx_y  , idx_z+1) + wt(4)
         CASE (xz_aligned_rectangle)
            pw%cr3d(idx_x+1, idx_y  , idx_z  ) = pw%cr3d(idx_x+1, idx_y  , idx_z  ) + wt(1)
            pw%cr3d(idx_x  , idx_y  , idx_z  ) = pw%cr3d(idx_x  , idx_y  , idx_z  ) + wt(2)
            pw%cr3d(idx_x  , idx_y  , idx_z+1) = pw%cr3d(idx_x  , idx_y  , idx_z+1) + wt(3)
            pw%cr3d(idx_x+1, idx_y  , idx_z+1) = pw%cr3d(idx_x+1, idx_y  , idx_z+1) + wt(4)
       END SELECT
    ELSE
       ! do nothing
    END IF
    
    CALL timestop(handle)

  END SUBROUTINE convert_grid_tile_to_grid

! *****************************************************************************
!> \brief  maps the weights of vertices of a voxel on a pw data type
!> \param voxel input voxel
!> \param pw the output plane wave data type
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE convert_voxel_to_grid(voxel, pw, error)

    TYPE(cs_voxel), INTENT(IN), POINTER      :: voxel
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_voxel_to_grid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idx_x, idx_y, idx_z
    REAL(dp), DIMENSION(8)                   :: wt

    CALL timeset(routineN,handle)

    IF (voxel%is_on) THEN
       idx_x = voxel%voxel_idx(1)   
       idx_y = voxel%voxel_idx(2)
       idx_z = voxel%voxel_idx(3)

       wt = REAL(voxel%weight, kind=dp)

       pw%cr3d(idx_x  , idx_y  , idx_z  ) = pw%cr3d(idx_x  , idx_y  , idx_z  ) + wt(1)
       pw%cr3d(idx_x+1, idx_y  , idx_z  ) = pw%cr3d(idx_x+1, idx_y  , idx_z  ) + wt(2)
       pw%cr3d(idx_x+1, idx_y+1, idx_z  ) = pw%cr3d(idx_x+1, idx_y+1, idx_z  ) + wt(3)
       pw%cr3d(idx_x  , idx_y+1, idx_z  ) = pw%cr3d(idx_x  , idx_y+1, idx_z  ) + wt(4)
       pw%cr3d(idx_x  , idx_y+1, idx_z+1) = pw%cr3d(idx_x  , idx_y+1, idx_z+1) + wt(5)
       pw%cr3d(idx_x  , idx_y  , idx_z+1) = pw%cr3d(idx_x  , idx_y  , idx_z+1) + wt(6)
       pw%cr3d(idx_x+1, idx_y  , idx_z+1) = pw%cr3d(idx_x+1, idx_y  , idx_z+1) + wt(7)
       pw%cr3d(idx_x+1, idx_y+1, idx_z+1) = pw%cr3d(idx_x+1, idx_y+1, idx_z+1) + wt(8)
    ELSE
       ! do nothing
    END IF
    
    CALL timestop(handle)

  END SUBROUTINE convert_voxel_to_grid

! *****************************************************************************
!> \brief   voxelizes a given rectangle
!> \param rectangle input rectangle
!> \param x_glbl x grid vetor of the simulation box
!> \param y_glbl y grid vetor of the simulation box
!> \param z_glbl z grid vetor of the simulation box
!> \param x_locl x grid vetor of the simulation box local to this process
!> \param y_locl y grid vetor of the simulation box local to this process
!> \param z_locl z grid vetor of the simulation box local to this process
!> \param pw pw grid data to be created
!> \param error cp2k error
!> \par History
!>       09.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE voxelize_rectangle(rectangle, x_glbl, y_glbl, z_glbl, &
                                           x_locl, y_locl, z_locl, pw, error)

    TYPE(cs_rectangle), INTENT(IN), POINTER  :: rectangle
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl, &
                                                x_locl, y_locl, z_locl
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'voxelize_rectangle', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, lb1, lb2, &
                                                lb3, ld_type, ub1, ub2, ub3
    INTEGER, DIMENSION(2, 3)                 :: aabb_extents, &
                                                aabb_extents_local, &
                                                bounds_local
    INTEGER, DIMENSION(3)                    :: grid_segment_index, &
                                                grid_tile_index, voxel_index
    REAL(dp), DIMENSION(3)                   :: dr, origin, P
    TYPE(cs_grid_segment), POINTER           :: grid_seg
    TYPE(cs_grid_tile), POINTER              :: grid_tile
    TYPE(cs_voxel), POINTER                  :: voxel

    CALL timeset(routineN,handle)

    dr = pw%pw_grid%dr
    bounds_local = pw%pw_grid%bounds_local
    lb1 = bounds_local(1,1);    ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);    ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);    ub3 = bounds_local(2,3)

    ld_type = local_domain_type(bounds_local)
    CALL compute_rectangle_aabb_extents(rectangle, x_glbl, y_glbl, z_glbl, &
                  x_locl, y_locl, z_locl, aabb_extents, aabb_extents_local, error)

    SELECT CASE (ld_type)
! if the current process's local domain is a box:
    CASE (axis_aligned_box)
       DO k = aabb_extents_local(1,3), aabb_extents_local(2,3)-1
          DO j = aabb_extents_local(1,2), aabb_extents_local(2,2)-1
             DO i = aabb_extents_local(1,1), aabb_extents_local(2,1)-1
                ALLOCATE(voxel)
                origin = (/x_locl(i),y_locl(j),z_locl(k)/)
                voxel_index = (/i,j,k/)
                CALL cs_voxel_create(voxel, origin, voxel_index, dr, error)
                CALL toggle_voxel(voxel, rectangle, error)
                CALL convert_voxel_to_grid(voxel, pw, error)
                CALL cs_voxel_release(voxel, error)
             END DO
          END DO
       END DO
! if the current process's local domain is a rectangle:
    CASE (yz_aligned_rectangle)
!       DO k = lb3, ub3-1
!          DO j = lb2, ub2-1
       DO k = aabb_extents_local(1,3), aabb_extents_local(2,3)-1
          DO j = aabb_extents_local(1,2), aabb_extents_local(2,2)-1
             ALLOCATE(grid_tile)
             origin = (/x_locl(lb1),y_locl(j),z_locl(k)/)
             grid_tile_index = (/lb1,j,k/)
             CALL cs_grid_tile_create(grid_tile, origin, grid_tile_index, dr, ld_type, error)
             CALL toggle_grid_tile(grid_tile, rectangle, error)
             CALL convert_grid_tile_to_grid(grid_tile, pw, error)
             CALL cs_grid_tile_release(grid_tile, error)
          END DO
       END DO
    CASE (xz_aligned_rectangle)
!       DO k = lb3, ub3-1
!          DO i = lb1, ub1-1
       DO k = aabb_extents_local(1,3), aabb_extents_local(2,3)-1
          DO i = aabb_extents_local(1,1), aabb_extents_local(2,1)-1
             ALLOCATE(grid_tile)
             origin = (/x_locl(i),y_locl(lb2),z_locl(k)/)
             grid_tile_index = (/i,lb2,k/)
             CALL cs_grid_tile_create(grid_tile, origin, grid_tile_index, dr, ld_type, error)
             CALL toggle_grid_tile(grid_tile, rectangle, error)
             CALL convert_grid_tile_to_grid(grid_tile, pw, error)
             CALL cs_grid_tile_release(grid_tile, error)
          END DO
       END DO
    CASE (xy_aligned_rectangle)
!       DO j = lb2, ub2-1
!          DO i = lb1, ub1-1
       DO j = aabb_extents_local(1,2), aabb_extents_local(2,2)-1
          DO i = aabb_extents_local(1,1), aabb_extents_local(2,1)-1
             ALLOCATE(grid_tile)
             origin = (/x_locl(i),y_locl(j),z_locl(lb3)/)
             grid_tile_index = (/i,j,lb3/)
             CALL cs_grid_tile_create(grid_tile, origin, grid_tile_index, dr, ld_type, error)
             CALL toggle_grid_tile(grid_tile, rectangle, error)
             CALL convert_grid_tile_to_grid(grid_tile, pw, error)
             CALL cs_grid_tile_release(grid_tile, error)
          END DO
       END DO
! if the current process's local domain is a segment:
    CASE (zaxis_aligned_segment)
!       DO k = lb3, ub3-1
       DO k = aabb_extents_local(1,3), aabb_extents_local(2,3)-1
          ALLOCATE(grid_seg)
          origin = (/x_locl(lb1),y_locl(lb2),z_locl(k)/)
          grid_segment_index = (/lb1,lb2,k/)
          CALL cs_grid_segment_create(grid_seg, origin, grid_segment_index, dr, ld_type, error)
          CALL toggle_grid_segment(grid_seg, rectangle, error)
          CALL convert_grid_segment_to_grid(grid_seg, pw, error)
          CALL cs_grid_segment_release(grid_seg, error)
       END DO
    CASE (yaxis_aligned_segment)
!       DO j = lb2, ub2-1
       DO j = aabb_extents_local(1,2), aabb_extents_local(2,2)-1
          ALLOCATE(grid_seg)
          origin = (/x_locl(lb1),y_locl(j),z_locl(lb3)/)
          grid_segment_index = (/lb1,j,lb3/)
          CALL cs_grid_segment_create(grid_seg, origin, grid_segment_index, dr, ld_type, error)
          CALL toggle_grid_segment(grid_seg, rectangle, error)
          CALL convert_grid_segment_to_grid(grid_seg, pw, error)
          CALL cs_grid_segment_release(grid_seg, error)
       END DO
    CASE (xaxis_aligned_segment)
!       DO i = lb1, ub1-1
       DO i = aabb_extents_local(1,1), aabb_extents_local(2,1)-1
          ALLOCATE(grid_seg)
          origin = (/x_locl(i),y_locl(lb2),z_locl(lb3)/)
          grid_segment_index = (/i,lb2,lb3/)
          CALL cs_grid_segment_create(grid_seg, origin, grid_segment_index, dr, ld_type, error)
          CALL toggle_grid_segment(grid_seg, rectangle, error)
          CALL convert_grid_segment_to_grid(grid_seg, pw, error)
          CALL cs_grid_segment_release(grid_seg, error)
       END DO
    CASE (single_grid_point)
! if the current process's local domain is a point:
       P = pw%cr3d(lb1,lb2,lb3)
       IF (is_inside_rectangle(P, rectangle)) pw%cr3d = 1.0_dp
    CASE (no_grid_point)
! if the current process does not hold any portion of the grid:
       ! do nothing
    END SELECT

! set all nonzero elements equal to 1 (consider every nonzero grid point only once)
    DO k = lb3, ub3
       DO j = lb2, ub2
          DO i = lb1, ub1
             IF (pw%cr3d(i,j,k) .NE. 0.0_dp) pw%cr3d(i,j,k) = 1.0_dp
          END DO
       END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE voxelize_rectangle

! *****************************************************************************
!> \brief   voxelizes a given axis-aligned cuboid (box)
!> \param cuboid the axis-aligned cuboid (box) to be voxelized
!> \param x_locl x grid vetor of the simulation box local to this process
!> \param y_locl y grid vetor of the simulation box local to this process
!> \param z_locl z grid vetor of the simulation box local to this process
!> \param pw pw grid data to be created
!> \param error cp2k error
!> \par History
!>       09.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE voxelize_aa_cuboid(cuboid, x_locl, y_locl, z_locl, pw, error)

    TYPE(cs_box), INTENT(IN), POINTER        :: cuboid
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_locl, y_locl, z_locl
    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'voxelize_aa_cuboid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, lb1, lb2, &
                                                lb3, ub1, ub2, ub3
    INTEGER, DIMENSION(2, 3)                 :: aac_extents_local, &
                                                bounds_local
    REAL(dp)                                 :: maxx, maxy, maxz, minx, miny, &
                                                minz

    CALL timeset(routineN,handle)

    bounds_local = pw%pw_grid%bounds_local
    lb1 = bounds_local(1,1);    ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);    ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);    ub3 = bounds_local(2,3)

    minx = MINVAL(cuboid%vertices(1,:))
    maxx = MAXVAL(cuboid%vertices(1,:))
    miny = MINVAL(cuboid%vertices(2,:))
    maxy = MAXVAL(cuboid%vertices(2,:))
    minz = MINVAL(cuboid%vertices(3,:))
    maxz = MAXVAL(cuboid%vertices(3,:))

    aac_extents_local(1,1) = find_lst_le(x_locl, minx)
    aac_extents_local(2,1) = find_fst_ge(x_locl, maxx)
    aac_extents_local(1,2) = find_lst_le(y_locl, miny)
    aac_extents_local(2,2) = find_fst_ge(y_locl, maxy)
    aac_extents_local(1,3) = find_lst_le(z_locl, minz)
    aac_extents_local(2,3) = find_fst_ge(z_locl, maxz)

    DO k = aac_extents_local(1,3), aac_extents_local(2,3)-1
       DO j = aac_extents_local(1,2), aac_extents_local(2,2)-1
          DO i = aac_extents_local(1,1), aac_extents_local(2,1)-1
             pw%cr3d(i,j,k) = 1.0_dp
          END DO
       END DO
    END DO

! set all nonzero elements equal to 1 (consider every nonzero grid point only once)
    DO k = lb3, ub3
       DO j = lb2, ub2
          DO i = lb1, ub1
             IF (pw%cr3d(i,j,k) .NE. 0.0_dp) pw%cr3d(i,j,k) = 1.0_dp
          END DO
       END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE voxelize_aa_cuboid

END MODULE dirichlet_bc_utils
