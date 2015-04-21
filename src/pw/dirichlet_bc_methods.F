!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief subroutines for defining and creating the Dirichlet parts of the domain
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE dirichlet_bc_methods

  USE dirichlet_bc_types,              ONLY: &
       AA_CUBOIDAL, AA_PLANAR, BWROT, CIRCUMSCRIBED, CYLINDRICAL, INSCRIBED, &
       PLANAR, cs_box_create, cs_rectangle, cs_rectangle_create, &
       cs_rectangle_release, cuboidal_tile, dirichlet_bc_p_type, &
       dirichlet_bc_parameters_dealloc, dirichlet_bc_type, rectangular_tile, &
       tile_p_type
  USE dirichlet_bc_utils,              ONLY: &
       convert_to_aa_rectangle, find_fst_ge, find_lst_le, linspace, &
       partition_aa_cuboid_into_tiles, partition_aa_rectangle_into_tiles, &
       rotate_rectangle, setup_grid_axes, voxelize_aa_cuboid, &
       voxelize_rectangle
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi,&
                                             twopi
  USE mathlib,                         ONLY: vector_product
  USE message_passing,                 ONLY: mp_sum
  USE physcon,                         ONLY: angstrom
  USE ps_implicit_types,               ONLY: MIXED_BC,&
                                             MIXED_PERIODIC_BC,&
                                             NEUMANN_BC,&
                                             PERIODIC_BC
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_zero
  USE pw_poisson_types,                ONLY: pw_poisson_parameter_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dirichlet_bc_methods'
  LOGICAL, PRIVATE, PARAMETER          :: debug = .FALSE.

  PUBLIC dirichlet_boundary_region_setup

  REAL(dp), PARAMETER, PRIVATE         :: small_value = 1.0E-8_dp

CONTAINS

! *****************************************************************************
!> \brief  Sets up the Dirichlet boundary condition 
!> \param pw_pool pool of plane wave grid 
!> \param poisson_params poisson_env parameters   
!> \param dbcs the DBC region to be created 
!> \param error cp2k error
!> \par History
!>       10.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dirichlet_boundary_region_setup(pw_pool, poisson_params, dbcs, error)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    TYPE(pw_poisson_parameter_type), &
      INTENT(INOUT)                          :: poisson_params
    TYPE(dirichlet_bc_p_type), &
      DIMENSION(:), INTENT(INOUT), POINTER   :: dbcs
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dirichlet_boundary_region_setup', &
      routineP = moduleN//':'//routineN

    INTEGER :: apx_type, dbc_id, handle, ind_end, ind_start, j, l, &
      n_aa_cuboidal, n_aa_planar, n_dbcs, n_planar, n_xaa_cylindrical, unit_nr
    INTEGER, DIMENSION(3)                    :: n_prtn
    INTEGER, DIMENSION(:), POINTER           :: xaa_cylindrical_nsides
    LOGICAL                                  :: smooth
    REAL(dp)                                 :: base_radius, v_D, zeta
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: x_glbl, x_locl, y_glbl, &
                                                y_locl, z_glbl, z_locl
    REAL(dp), DIMENSION(2)                   :: base_center
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pw_grid_type), POINTER              :: pw_grid

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
       unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
    ELSE
       unit_nr = -1
    ENDIF   

    n_aa_planar = poisson_params%dbc_params%n_aa_planar
    n_aa_cuboidal = poisson_params%dbc_params%n_aa_cuboidal
    n_planar = poisson_params%dbc_params%n_planar
    n_xaa_cylindrical = poisson_params%dbc_params%n_xaa_cylindrical
    xaa_cylindrical_nsides => poisson_params%dbc_params%xaa_cylindrical_nsides
    pw_grid => pw_pool%pw_grid
    CALL setup_grid_axes(pw_grid, x_glbl, y_glbl, z_glbl, x_locl, y_locl, z_locl, error)
    n_dbcs = n_aa_planar + n_aa_cuboidal + n_planar + SUM(xaa_cylindrical_nsides)

    SELECT CASE (poisson_params%ps_implicit_params%boundary_condition)
    CASE (MIXED_BC, MIXED_PERIODIC_BC)
      IF (n_dbcs .EQ. 0) THEN 
         CALL cp_assert(.FALSE., cp_failure_level, cp_assertion_failed, routineP, &
                        "Dirichlet Boundary regions are not defined.", error)
      END IF

      ALLOCATE(dbcs(n_dbcs))
      IF (unit_nr .GT. 0) THEN 
         WRITE(unit_nr,'(/,T3,A,A,/,T3,A)') "POISSON| IMPLICIT (GENERALIZED) SOLVER ", REPEAT('-', 39), &
                                            "POISSON| Preparing Dirichlet boundary regions ..."
      END IF

      DO j = 1, n_aa_planar
         ALLOCATE(dbcs(j)%dirichlet_bc)
         n_prtn = poisson_params%dbc_params%aa_planar_nprtn(:,j)
         dbc_id = AA_PLANAR + j
         v_D = poisson_params%dbc_params%aa_planar_vD(j)
         smooth = poisson_params%dbc_params%aa_planar_smooth(j)
         zeta = poisson_params%dbc_params%aa_planar_zeta(j)

         IF (unit_nr .GT. 0) THEN
            WRITE(unit_nr,'(T3,A,I5)') "POISSON| DBC", j
            WRITE(unit_nr,'(T3,A)') "POISSON|     type : AXIS-ALIGNED PLANAR"
            WRITE(unit_nr,'(T3,A,E13.4,2X,A)') "POISSON|     applied potential :", v_D, "[Eh/e]"
         END IF
         CALL create_aa_planar_dbc(x_glbl, y_glbl, z_glbl, &
                                   poisson_params%dbc_params%aa_planar_xxtnt(:,j), &
                                   poisson_params%dbc_params%aa_planar_yxtnt(:,j), &
                                   poisson_params%dbc_params%aa_planar_zxtnt(:,j), &
                                   v_D, smooth, zeta, dbc_id, dbcs(j)%dirichlet_bc, error)

         CALL dirichlet_bc_partition(v_D, smooth, zeta, n_prtn, pw_pool, x_glbl, y_glbl, z_glbl, &
                                     x_locl, y_locl, z_locl, dbcs(j)%dirichlet_bc, error)
      END DO

      l = n_aa_planar
      DO j = l + 1, l + n_aa_cuboidal
         ALLOCATE(dbcs(j)%dirichlet_bc)
         n_prtn = poisson_params%dbc_params%aa_cuboidal_nprtn(:,j - l)
         dbc_id = AA_CUBOIDAL + j - l
         v_D = poisson_params%dbc_params%aa_cuboidal_vD(j - l)
         smooth = poisson_params%dbc_params%aa_cuboidal_smooth(j - l)
         zeta = poisson_params%dbc_params%aa_cuboidal_zeta(j - l)

         IF (unit_nr .GT. 0) THEN
            WRITE(unit_nr,'(T3,A,I5)') "POISSON| DBC", j
            WRITE(unit_nr,'(T3,A)') "POISSON|     type : AXIS-ALIGNED CUBOIDAL"
            WRITE(unit_nr,'(T3,A,E13.4,2X,A)') "POISSON|     applied potential :", v_D, "[Eh/e]"
         END IF
         CALL create_aa_cuboidal_dbc(x_glbl, y_glbl, z_glbl, &
                                   poisson_params%dbc_params%aa_cuboidal_xxtnt(:,j - l), &
                                   poisson_params%dbc_params%aa_cuboidal_yxtnt(:,j - l), &
                                   poisson_params%dbc_params%aa_cuboidal_zxtnt(:,j - l), &
                                   v_D, smooth, zeta, dbc_id, dbcs(j)%dirichlet_bc, error)

         CALL dirichlet_bc_partition(v_D, smooth, zeta, n_prtn, pw_pool, x_glbl, y_glbl, z_glbl, &
                                     x_locl, y_locl, z_locl, dbcs(j)%dirichlet_bc, error)
      END DO

      l = n_aa_planar + n_aa_cuboidal
      DO j = l + 1, l + n_planar
         ALLOCATE(dbcs(j)%dirichlet_bc)
         n_prtn = poisson_params%dbc_params%planar_nprtn(:,j - l)
         dbc_id = PLANAR + j - l
         v_D = poisson_params%dbc_params%planar_vD(j - l)
         smooth = poisson_params%dbc_params%planar_smooth(j - l)
         zeta = poisson_params%dbc_params%planar_zeta(j - l)

         IF (unit_nr .GT. 0) THEN
            WRITE(unit_nr,'(T3,A,I5)') "POISSON| DBC", j
            WRITE(unit_nr,'(T3,A)') "POISSON|     type : PLANAR"
            WRITE(unit_nr,'(T3,A,E13.4,2X,A)') "POISSON|     applied potential :", v_D, "[Eh/e]"
         END IF
         CALL create_arbitrary_planar_dbc(x_glbl, y_glbl, z_glbl, &
                                   poisson_params%dbc_params%planar_Avtx(:,j - l), &
                                   poisson_params%dbc_params%planar_Bvtx(:,j - l), &
                                   poisson_params%dbc_params%planar_Cvtx(:,j - l), &
                                   v_D, smooth, zeta, dbc_id, dbcs(j)%dirichlet_bc, error)

         CALL dirichlet_bc_partition(v_D, smooth, zeta, n_prtn, pw_pool, x_glbl, y_glbl, z_glbl, &
                                     x_locl, y_locl, z_locl, dbcs(j)%dirichlet_bc, error)
      END DO

      l = n_aa_planar + n_aa_cuboidal + n_planar
      DO j = 1, n_xaa_cylindrical
         ind_start = l + 1
         ind_end   = l + xaa_cylindrical_nsides(j)

         n_prtn = poisson_params%dbc_params%xaa_cylindrical_nprtn(:,j)
         base_center = poisson_params%dbc_params%xaa_cylindrical_bctr(:,j)
         base_radius = poisson_params%dbc_params%xaa_cylindrical_brad(j)
         v_D = poisson_params%dbc_params%xaa_cylindrical_vD(j)
         smooth = poisson_params%dbc_params%xaa_cylindrical_smooth(j)
         zeta = poisson_params%dbc_params%xaa_cylindrical_zeta(j)
         apx_type = poisson_params%dbc_params%xaa_cylindrical_apxtyp(j)

         IF (unit_nr .GT. 0) THEN
            WRITE(unit_nr,'(T3,A,I5)') "POISSON| DBC", l + j
            WRITE(unit_nr,'(T3,A)') "POISSON|     type : X-AXIS-ALIGNED CYLINDRICAL"
            WRITE(unit_nr,'(T3,A,E13.4,2X,A)') "POISSON|     applied potential :", v_D, "[Eh/e]"
         END IF
         CALL create_xaa_cylindrical_dbc(pw_pool, x_glbl, y_glbl, z_glbl, &
                                   x_locl, y_locl, z_locl, &
                                   poisson_params%dbc_params%xaa_cylindrical_xxtnt(:,j), &
                                   base_center, base_radius, v_D, smooth, zeta, n_prtn, &
                                   dbcs(ind_start : ind_end), apx_type, error)

         l = l + xaa_cylindrical_nsides(j)
      END DO

    CASE (PERIODIC_BC, NEUMANN_BC)
       ! do nothing
    END SELECT

! we won't need parameters anymore so deallocate them 
    CALL dirichlet_bc_parameters_dealloc(poisson_params%dbc_params, error)

    CALL timestop(handle)

  END SUBROUTINE dirichlet_boundary_region_setup

! *****************************************************************************
!> \brief Partitions a dirichlet_bc_type into tiles  
!> \param v_D potential imposed at the Dirichlet part of the boundary (gate voltage)
!> \param smooth whether or not to smoothen the dirichlet region 
!> \param zeta the mollifier's width
!> \param n_prtn vetor of size 3 specifying the number of times that the x, y and
!>               z interval (defining the region) should be partitioned into  
!> \param pw_pool pool of plane-wave grid
!> \param x_glbl x grid vetor of the simulation box  
!> \param y_glbl y grid vetor of the simulation box  
!> \param z_glbl z grid vetor of the simulation box  
!> \param x_locl x grid vetor of the simulation box local to this process 
!> \param y_locl y grid vetor of the simulation box local to this process 
!> \param z_locl z grid vetor of the simulation box local to this process 
!> \param dirichlet_bc the dirichlet_bc object to be partitioned 
!> \param error cp2k error
!>
!> \par History
!>       10.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dirichlet_bc_partition(v_D, smooth, zeta, n_prtn, pw_pool, x_glbl, y_glbl, z_glbl, &
                                    x_locl, y_locl, z_locl, dirichlet_bc, error)

    REAL(dp), INTENT(IN)                     :: v_D
    LOGICAL, INTENT(IN)                      :: smooth
    REAL(dp), INTENT(IN)                     :: zeta
    INTEGER, DIMENSION(3), INTENT(IN)        :: n_prtn
    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl, &
                                                x_locl, y_locl, z_locl
    TYPE(dirichlet_bc_type), INTENT(INOUT), &
      POINTER                                :: dirichlet_bc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dirichlet_bc_partition', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, k, n_tiles, &
                                                tile_npts, unit_nr
    REAL(dp)                                 :: phi1, phi2
    REAL(dp), DIMENSION(3)                   :: rot_axis1, rot_axis2
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cs_rectangle), POINTER              :: rectangle_aa
    TYPE(tile_p_type), DIMENSION(:), POINTER :: tiles

    CALL timeset(routineN,handle)
    
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
       unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
    ELSE
       unit_nr = -1
    ENDIF   

    SELECT CASE (dirichlet_bc%dbc_geom)
    CASE (AA_PLANAR, PLANAR, CYLINDRICAL)
! In order to partition a rectangular shaped dirichlet_bc_type:
!   (1) Choose a plane parallel to the coordinate planes
!   (2) Rotate the rectangle onto the plane
!   (3) If the result is not axis-aligned rotate it on the plane until it becomes one
!   (4) Partition the final axis_aligned rectangle 
!   (5) Rotate back the rectangle and all the tiles created in the last step 

      ALLOCATE(rectangle_aa)
      CALL convert_to_aa_rectangle(dirichlet_bc%rectangle, rectangle_aa, &
                                   phi1, rot_axis1, phi2, rot_axis2, error)
      CALL partition_aa_rectangle_into_tiles(rectangle_aa, x_glbl, y_glbl, z_glbl, n_prtn, tiles, error)
      CALL cs_rectangle_release(rectangle_aa, error=error)
  
      n_tiles = SIZE(tiles)
      dirichlet_bc%n_tiles = n_tiles
  
      ALLOCATE(dirichlet_bc%tiles(n_tiles))
      IF ((unit_nr .GT. 0) .AND. debug) WRITE(unit_nr,'(T7,A,I5)') "number of partitions : ", n_tiles
  
      DO k = 1, n_tiles
         ALLOCATE(dirichlet_bc%tiles(k)%tile)
         ALLOCATE(dirichlet_bc%tiles(k)%tile%rectangle)
         CALL rotate_rectangle(tiles(k)%tile%rectangle, phi1, rot_axis1, BWROT, &
                               dirichlet_bc%tiles(k)%tile%rectangle, error)
  
         dirichlet_bc%tiles(k)%tile%tile_id = 8000 + k
         dirichlet_bc%tiles(k)%tile%v_D     = v_D
         dirichlet_bc%tiles(k)%tile%smooth  = smooth
         dirichlet_bc%tiles(k)%tile%tile_geom = rectangular_tile
         dirichlet_bc%tiles(k)%tile%mollifier_zeta = zeta
   
         CALL pw_pool_create_pw(pw_pool, dirichlet_bc%tiles(k)%tile%tile_pw, &
                                use_data=REALDATA3D, in_space=REALSPACE, error=error)
         CALL pw_zero(dirichlet_bc%tiles(k)%tile%tile_pw, error=error)
  
         IF ((unit_nr .GT. 0) .AND. debug) THEN 
            WRITE(unit_nr,'(T7,A,I5)') "tile", k
            DO i = 1, 4
               WRITE(unit_nr,'(T10,A,I1,3F10.3)') "   vertex ", i, &
                             angstrom * dirichlet_bc%tiles(k)%tile%rectangle%vertices(:,i)
            END DO
         END IF
  
         CALL voxelize_rectangle(dirichlet_bc%tiles(k)%tile%rectangle, x_glbl, y_glbl, z_glbl, &
                                 x_locl, y_locl, z_locl, dirichlet_bc%tiles(k)%tile%tile_pw, error)
         tile_npts = NINT(SUM(dirichlet_bc%tiles(k)%tile%tile_pw%cr3d),KIND=KIND(tile_npts))
         CALL mp_sum(tile_npts, pw_pool%pw_grid%para%group)
         dirichlet_bc%tiles(k)%tile%npts = tile_npts
      END DO
  
      DO k = 1, n_tiles
         CALL cs_rectangle_release(tiles(k)%tile%rectangle, error)
         DEALLOCATE(tiles(k)%tile)
      END DO
      DEALLOCATE(tiles)

    CASE (AA_CUBOIDAL)

      CALL partition_aa_cuboid_into_tiles(dirichlet_bc%box, x_glbl, y_glbl, z_glbl, n_prtn, tiles, error)
      n_tiles = SIZE(tiles)
      dirichlet_bc%n_tiles = n_tiles
      ALLOCATE(dirichlet_bc%tiles(n_tiles))
      dirichlet_bc%tiles = tiles
      IF ((unit_nr .GT. 0) .AND. debug) WRITE(unit_nr,'(T7,A,I5)') "number of partitions : ", n_tiles
  
      DO k = 1, n_tiles
         ! update id and initialize missing info 
         dirichlet_bc%tiles(k)%tile%tile_id = 8000 + k
         dirichlet_bc%tiles(k)%tile%v_D     = v_D
         dirichlet_bc%tiles(k)%tile%smooth  = smooth
         dirichlet_bc%tiles(k)%tile%tile_geom = cuboidal_tile
         dirichlet_bc%tiles(k)%tile%mollifier_zeta = zeta
   
         CALL pw_pool_create_pw(pw_pool, dirichlet_bc%tiles(k)%tile%tile_pw, &
                                use_data=REALDATA3D, in_space=REALSPACE, error=error)
         CALL pw_zero(dirichlet_bc%tiles(k)%tile%tile_pw, error=error)
  
         IF ((unit_nr .GT. 0) .AND. debug) THEN 
            WRITE(unit_nr,'(T7,A,I5)') "tile", k
            DO i = 1, 8
               WRITE(unit_nr,'(T10,A,I1,3F10.3)') &
                             "   vertex ", i, angstrom * dirichlet_bc%box%vertices(:,i)
            END DO
         END IF
  
         CALL voxelize_aa_cuboid(dirichlet_bc%tiles(k)%tile%box, x_locl, y_locl, z_locl, &
                                 dirichlet_bc%tiles(k)%tile%tile_pw, error)
         tile_npts = NINT(SUM(dirichlet_bc%tiles(k)%tile%tile_pw%cr3d),KIND=KIND(tile_npts))
         CALL mp_sum(tile_npts, pw_pool%pw_grid%para%group)
         dirichlet_bc%tiles(k)%tile%npts = tile_npts
      END DO
  
      DEALLOCATE(tiles)

    END SELECT

    CALL timestop(handle)

  END SUBROUTINE dirichlet_bc_partition

! *****************************************************************************
!> \brief   Creates an x-axis-aligned cylindrical gate region. The cylindrical shape 
!>     is approximated by an n-gonal (circumscribed or inscribed) uniform prism-like 
!>     object. The circular base is approximated by a polygon. A second polygon is 
!>     then created by rotating the first one by delta_alpha radians. Let's name the 
!>     vertices of the first and the second polygon as v1, v2, v3, ... vn-1, vn and 
!>     w1, w2, w3, ... wn-1, wn, respectively. Then w1v2, w2v3, ... wn-1vn, wnv1 form
!>     the edges of the cross-section of the approximating object. This way the dbcs
!>     do not share any edge.   
!> \param pw_pool  pool of plane wave grid 
!> \param x_glbl x grid vetor of the simulation box  
!> \param y_glbl y grid vetor of the simulation box  
!> \param z_glbl z grid vetor of the simulation box  
!> \param x_locl x grid vetor of the simulation box local to this process 
!> \param y_locl y grid vetor of the simulation box local to this process 
!> \param z_locl z grid vetor of the simulation box local to this process 
!> \param x_xtnt the x extent of the cylindrical region 
!> \param base_center the y and z coordinates of the cylinder's base 
!> \param base_radius the radius of the cylinder's base 
!> \param v_D potential imposed at the Dirichlet part of the boundary (gate voltage)
!> \param smooth whether or not to smoothen the dirichlet region 
!> \param zeta the mollifier's width
!> \param n_prtn vetor of size 3 specifying the number of times that the x, y and
!>               z interval (defining the region) should be partitioned into
!> \param dbcs  the x-axis-aligned cylindrical gate region to be created
!> \param apx_type the type of the n-gonal prism approximating the cylinder
!> \param error cp2k error
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE create_xaa_cylindrical_dbc(pw_pool, &
                       x_glbl, y_glbl, z_glbl, x_locl, y_locl, z_locl, &
                       x_xtnt, base_center, base_radius, v_D, smooth, zeta, &
                       n_prtn, dbcs, apx_type, error)

    TYPE(pw_pool_type), POINTER              :: pw_pool
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl, &
                                                x_locl, y_locl, z_locl
    REAL(dp), DIMENSION(2), INTENT(IN)       :: x_xtnt, base_center
    REAL(dp), INTENT(IN)                     :: base_radius, v_D
    LOGICAL, INTENT(IN)                      :: smooth
    REAL(dp), INTENT(IN)                     :: zeta
    INTEGER, DIMENSION(3), INTENT(IN)        :: n_prtn
    TYPE(dirichlet_bc_p_type), &
      DIMENSION(:), INTENT(INOUT)            :: dbcs
    INTEGER, INTENT(IN), OPTIONAL            :: apx_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_xaa_cylindrical_dbc', &
      routineP = moduleN//':'//routineN

    INTEGER :: glb1, glb2, glb3, gub1, gub2, gub3, handle, i, &
      intern_apx_type, j, n_dbcs, unit_nr, xxtnt_fst, xxtnt_ult
    LOGICAL                                  :: forb_xtnt1, forb_xtnt2, &
                                                forb_xtnt3, forb_xtnt4
    REAL(dp)                                 :: delta_alpha, dx, dy, dz, h, &
                                                Lx, Ly, Lz, theta, xlb, xub
    REAL(dp), DIMENSION(3)                   :: A, B, C, D
    TYPE(cp_logger_type), POINTER            :: logger
    REAL(dp), DIMENSION(SIZE(dbcs)+1)        :: alpha, alpha_rotd, ylb, yub, &
                                                zlb, zub

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
       unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
    ELSE
       unit_nr = -1
    ENDIF

    glb1 = LBOUND(x_glbl,1);     gub1 = UBOUND(x_glbl,1)
    glb2 = LBOUND(y_glbl,1);     gub2 = UBOUND(y_glbl,1)
    glb3 = LBOUND(z_glbl,1);     gub3 = UBOUND(z_glbl,1)

    dx = x_glbl(glb1+1) - x_glbl(glb1)
    dy = y_glbl(glb2+1) - y_glbl(glb2)
    dz = z_glbl(glb3+1) - z_glbl(glb3)
    Lx = x_glbl(gub1) - x_glbl(glb1)
    Ly = y_glbl(gub2) - y_glbl(glb2)
    Lz = z_glbl(gub3) - z_glbl(glb3)

    IF( (x_xtnt(1) .LT. x_glbl(glb1)) .OR. (x_xtnt(2) .GT. x_glbl(gub1)+dx) ) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP, &
                      "The length of the cylindrical Dirichlet region is larger than the "//&
                      "x range of the simulation cell.", error)
    END IF
    forb_xtnt1 = base_center(1)-base_radius .LT. x_glbl(glb1)
    forb_xtnt2 = base_center(1)+base_radius .GT. x_glbl(gub1)+dx
    forb_xtnt3 = base_center(2)-base_radius .LT. y_glbl(glb2)
    forb_xtnt4 = base_center(2)+base_radius .GT. y_glbl(gub2)+dy
    IF( forb_xtnt1 .OR. forb_xtnt2 .OR. forb_xtnt3 .OR. forb_xtnt4 ) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP, &
                      "The cylinder does not fit entirely inside the simulation cell.", error)
    END IF

    intern_apx_type = CIRCUMSCRIBED
    IF (PRESENT(apx_type)) intern_apx_type = apx_type

    n_dbcs = SIZE(dbcs)

    xxtnt_fst = find_fst_ge(x_glbl, x_xtnt(1))
    xxtnt_ult = find_lst_le(x_glbl, x_xtnt(2))

    xlb = x_glbl(xxtnt_fst)
    xub = x_glbl(xxtnt_ult)
    
    theta = twopi/n_dbcs

    IF (intern_apx_type .EQ. INSCRIBED) THEN
       h = base_radius                ! inscribed uniform prism
    ELSE IF (intern_apx_type .EQ. CIRCUMSCRIBED) THEN
       h = base_radius/COS(0.5*theta) ! circumscribed uniform prism
    ELSE 
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "Unknown approximation type for cylinder.",error)
    END IF
    IF( h .GT. MINVAL((/Ly, Lz/)) ) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP, &
                      "Reduce the base radius!", error)
    END IF

    delta_alpha = 0.05_dp
    alpha = linspace(0.0_dp, 2*pi, n_dbcs+1)
    alpha_rotd = alpha + delta_alpha;

    DO j = 1, n_dbcs
       ylb(j) = base_center(1) + h*SIN(alpha(j)); 
       zlb(j) = base_center(2) + h*COS(alpha(j));
       yub(j) = base_center(1) + h*SIN(alpha_rotd(j));
       zub(j) = base_center(2) + h*COS(alpha_rotd(j));
    END DO
    ylb(n_dbcs+1) = ylb(1);
    yub(n_dbcs+1) = yub(1);
    zlb(n_dbcs+1) = zlb(1);
    zub(n_dbcs+1) = zub(1);

    DO j = 1, n_dbcs
       ALLOCATE(dbcs(j)%dirichlet_bc)
       dbcs(j)%dirichlet_bc%dbc_id = CYLINDRICAL + j
       dbcs(j)%dirichlet_bc%v_D = v_D
       dbcs(j)%dirichlet_bc%smooth = smooth
       dbcs(j)%dirichlet_bc%mollifier_zeta = zeta
       dbcs(j)%dirichlet_bc%dbc_geom = CYLINDRICAL

       A = (/ xlb  , yub(j)  , zub(j)   /)
       B = (/ xlb  , ylb(j+1), zlb(j+1) /)
       C = (/ xub  , ylb(j+1), zlb(j+1) /)
       D = (/ xub  , yub(j)  , zub(j)   /)
     
       ALLOCATE(dbcs(j)%dirichlet_bc%rectangle)
       CALL cs_rectangle_create(dbcs(j)%dirichlet_bc%rectangle, A, B, C, D, error)

       dbcs(j)%dirichlet_bc%n_tiles = 1

       IF ((unit_nr .GT. 0) .AND. debug) THEN 
          WRITE(unit_nr,'(T7,A,I5,T20,A,I5,A)') "edge", j, "of the", n_dbcs, &
                                                "-gonal prism approximating the cylinder" 
          DO i = 1, 4
             WRITE(unit_nr,'(T10,A,I1,3F10.3)') "   vertex ", i, &
                                      angstrom * dbcs(j)%dirichlet_bc%rectangle%vertices(:,i)
          END DO
       END IF

       CALL dirichlet_bc_partition(v_D, smooth, zeta, n_prtn, pw_pool, x_glbl, y_glbl, z_glbl, &
                                   x_locl, y_locl, z_locl, dbcs(j)%dirichlet_bc, error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE create_xaa_cylindrical_dbc

! *****************************************************************************
!> \brief   Creates an axis-aligned planar gate region.
!> \param x_glbl x grid vetor of the simulation box  
!> \param y_glbl y grid vetor of the simulation box  
!> \param z_glbl z grid vetor of the simulation box  
!> \param x_xtnt the x extent of the planar region 
!> \param y_xtnt the y extent of the planar region 
!> \param z_xtnt the z extent of the planar region 
!> \param v_D potential imposed at the Dirichlet part of the boundary (gate voltage)
!> \param smooth whether or not to smoothen the dirichlet region 
!> \param zeta the mollifier's width
!> \param dbc_id unique ID for the planar Dirichlet region
!> \param dirichlet_bc the dirichlet_bc object to be created 
!> \param error cp2k error
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE create_aa_planar_dbc(x_glbl, y_glbl, z_glbl, x_xtnt, y_xtnt, z_xtnt, &
                       v_D, smooth, zeta, dbc_id, dirichlet_bc, error)

    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl
    REAL(dp), DIMENSION(2), INTENT(IN)       :: x_xtnt, y_xtnt, z_xtnt
    REAL(dp), INTENT(IN)                     :: v_D
    LOGICAL, INTENT(IN)                      :: smooth
    REAL(dp), INTENT(IN)                     :: zeta
    INTEGER, INTENT(IN)                      :: dbc_id
    TYPE(dirichlet_bc_type), INTENT(INOUT), &
      POINTER                                :: dirichlet_bc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_aa_planar_dbc', &
      routineP = moduleN//':'//routineN
    LOGICAL, DIMENSION(6), PARAMETER         :: test_forb_xtnts = .TRUE. 

    INTEGER :: handle, i, n_forb_xtnts, unit_nr, xxtnt_fst, xxtnt_ult, &
      yxtnt_fst, yxtnt_ult, zxtnt_fst, zxtnt_ult
    LOGICAL                                  :: forb_xtnt1, forb_xtnt2, &
                                                forb_xtnt3, forb_xtnt4, &
                                                forb_xtnt5, forb_xtnt6
    REAL(dp)                                 :: dx, dy, dz, xlb, xub, ylb, &
                                                yub, zlb, zub
    REAL(dp), DIMENSION(3)                   :: A, B, C, D
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
       unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
    ELSE
       unit_nr = -1
    ENDIF

    dirichlet_bc%dbc_id = dbc_id
    dirichlet_bc%dbc_geom = AA_PLANAR
    dirichlet_bc%v_D = v_D
    dirichlet_bc%smooth = smooth
    dirichlet_bc%mollifier_zeta = zeta

    dx = x_glbl(LBOUND(x_glbl,1)+1) - x_glbl(LBOUND(x_glbl,1))
    dy = y_glbl(LBOUND(y_glbl,1)+1) - y_glbl(LBOUND(y_glbl,1))
    dz = z_glbl(LBOUND(z_glbl,1)+1) - z_glbl(LBOUND(z_glbl,1))

    forb_xtnt1 = x_xtnt(1) .LT. x_glbl(LBOUND(x_glbl,1))
    forb_xtnt2 = x_xtnt(2) .GT. x_glbl(UBOUND(x_glbl,1))+dx
    forb_xtnt3 = y_xtnt(1) .LT. y_glbl(LBOUND(y_glbl,1))
    forb_xtnt4 = y_xtnt(2) .GT. y_glbl(UBOUND(y_glbl,1))+dy
    forb_xtnt5 = z_xtnt(1) .LT. z_glbl(LBOUND(z_glbl,1))
    forb_xtnt6 = z_xtnt(2) .GT. z_glbl(UBOUND(z_glbl,1))+dz
    n_forb_xtnts = COUNT((/forb_xtnt1, forb_xtnt2, forb_xtnt3, &
                           forb_xtnt4, forb_xtnt5, forb_xtnt6/) .EQV. test_forb_xtnts)
    IF (n_forb_xtnts .GT. 0) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "The given extents for the Dirichlet part of the boundary are outside "//&
                      "the range of the simulation cell.", error)
    END IF 

    xxtnt_fst = find_fst_ge(x_glbl, x_xtnt(1))
    xxtnt_ult = find_lst_le(x_glbl, x_xtnt(2))
    yxtnt_fst = find_fst_ge(y_glbl, y_xtnt(1))
    yxtnt_ult = find_lst_le(y_glbl, y_xtnt(2))
    zxtnt_fst = find_fst_ge(z_glbl, z_xtnt(1))
    zxtnt_ult = find_lst_le(z_glbl, z_xtnt(2))

    xlb = x_glbl(xxtnt_fst)
    xub = x_glbl(xxtnt_ult)
    ylb = y_glbl(yxtnt_fst)
    yub = y_glbl(yxtnt_ult)
    zlb = z_glbl(zxtnt_fst)
    zub = z_glbl(zxtnt_ult)
    
    IF (ABS(z_xtnt(1)-z_xtnt(2)) .LE. small_value) THEN
       A = (/ xlb, yub, zlb /)
       B = (/ xlb, ylb, zlb /)
       C = (/ xub, ylb, zlb /)
       D = (/ xub, yub, zlb /)
       ALLOCATE(dirichlet_bc%rectangle)
       CALL cs_rectangle_create(dirichlet_bc%rectangle, A, B, C, D, error)
       dirichlet_bc%n_tiles = 1
    ELSE IF (ABS(y_xtnt(1)-y_xtnt(2)) .LE. small_value) THEN
       A = (/ xlb, ylb, zub /)
       B = (/ xlb, ylb, zlb /)
       C = (/ xub, ylb, zlb /)
       D = (/ xub, ylb, zub /)
       ALLOCATE(dirichlet_bc%rectangle)
       CALL cs_rectangle_create(dirichlet_bc%rectangle, A, B, C, D, error)
       dirichlet_bc%n_tiles = 1
    ELSE IF (ABS(x_xtnt(1)-x_xtnt(2)) .LE. small_value) THEN
       A = (/ xlb, ylb, zub /)
       B = (/ xlb, ylb, zlb /)
       C = (/ xlb, yub, zlb /)
       D = (/ xlb, yub, zub /)
       ALLOCATE(dirichlet_bc%rectangle)
       CALL cs_rectangle_create(dirichlet_bc%rectangle, A, B, C, D, error)
       dirichlet_bc%n_tiles = 1
    ELSE 
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "According to the given extents, the Dirichlet region is not axis-aligned or a plane. "//&
                      "Try the input section PLANAR or AA_CUBOIDAL instead.",error)
    END IF
      
    IF ((unit_nr .GT. 0) .AND. debug) THEN 
       DO i = 1, 4
          WRITE(unit_nr,'(T10,A,I1,3F10.3)') "   vertex ", i, &
                                   angstrom * dirichlet_bc%rectangle%vertices(:,i)
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE create_aa_planar_dbc

! *****************************************************************************
!> \brief   Creates an arbitrary planar gate given the coordinates of its rectangular shape.
!> \param x_glbl x grid vetor of the simulation box  
!> \param y_glbl y grid vetor of the simulation box  
!> \param z_glbl z grid vetor of the simulation box  
!> \param A coordinates of the vertex A
!> \param B coordinates of the vertex B
!> \param C coordinates of the vertex C
!> \param v_D potential imposed at the Dirichlet part of the boundary (gate voltage)
!> \param smooth whether or not to smoothen the dirichlet region 
!> \param zeta the mollifier's width
!> \param dbc_id unique ID for the planar Dirichlet region
!> \param dirichlet_bc the dirichlet_bc object to be created 
!> \param error cp2k error
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE create_arbitrary_planar_dbc(x_glbl, y_glbl, z_glbl, A, B, C, v_D, &
                       smooth, zeta, dbc_id, dirichlet_bc, error)

    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl
    REAL(dp), DIMENSION(3), INTENT(IN)       :: A, B, C
    REAL(dp), INTENT(IN)                     :: v_D
    LOGICAL, INTENT(IN)                      :: smooth
    REAL(dp), INTENT(IN)                     :: zeta
    INTEGER, INTENT(IN)                      :: dbc_id
    TYPE(dirichlet_bc_type), INTENT(INOUT), &
      POINTER                                :: dirichlet_bc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_arbitrary_planar_dbc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: glb1, glb2, glb3, gub1, gub2, &
                                                gub3, handle, i, unit_nr
    LOGICAL :: A_is_inside, are_coplanar, are_orthogonal, B_is_inside, &
      C_is_inside, D_is_inside, is_rectangle
    REAL(dp)                                 :: dist1, dist2, dist3, dist4, &
                                                dx, dy, dz
    REAL(dp), DIMENSION(3)                   :: AB, AC, AD, BC, cm, D
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
       unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
    ELSE
       unit_nr = -1
    ENDIF

    glb1 = LBOUND(x_glbl,1);     gub1 = UBOUND(x_glbl,1)
    glb2 = LBOUND(y_glbl,1);     gub2 = UBOUND(y_glbl,1)
    glb3 = LBOUND(z_glbl,1);     gub3 = UBOUND(z_glbl,1)
    dx = x_glbl(glb1+1) - x_glbl(glb1)
    dy = y_glbl(glb2+1) - y_glbl(glb2)
    dz = z_glbl(glb3+1) - z_glbl(glb3)

    D = A + (C - B)
    AB = B - A;     AC = C - A;     AD = D - A;     BC = C - B;
    are_orthogonal = ABS( DOT_PRODUCT(AB, BC)/(SQRT(SUM(AB**2)) * SQRT(SUM(BC**2))) ) .LE. small_value
    IF (.NOT. are_orthogonal) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP, &
            "The given vertices for defining a planar Dirichlet region do not form orthogonal edges.", error)
    END IF

    A_is_inside = ((A(1) .GT. x_glbl(glb1)) .AND. (A(1) .LT. x_glbl(gub1)+dx)) .AND. &
                  ((A(2) .GT. y_glbl(glb2)) .AND. (A(2) .LT. y_glbl(gub2)+dy)) .AND. &
                  ((A(3) .GT. z_glbl(glb3)) .AND. (A(3) .LT. z_glbl(gub3)+dz)) 
    B_is_inside = ((B(1) .GT. x_glbl(glb1)) .AND. (B(1) .LT. x_glbl(gub1)+dx)) .AND. &
                  ((B(2) .GT. y_glbl(glb2)) .AND. (B(2) .LT. y_glbl(gub2)+dy)) .AND. &
                  ((B(3) .GT. z_glbl(glb3)) .AND. (B(3) .LT. z_glbl(gub3)+dz)) 
    C_is_inside = ((C(1) .GT. x_glbl(glb1)) .AND. (C(1) .LT. x_glbl(gub1)+dx)) .AND. &
                  ((C(2) .GT. y_glbl(glb2)) .AND. (C(2) .LT. y_glbl(gub2)+dy)) .AND. &
                  ((C(3) .GT. z_glbl(glb3)) .AND. (C(3) .LT. z_glbl(gub3)+dz)) 
    D_is_inside = ((D(1) .GT. x_glbl(glb1)) .AND. (D(1) .LT. x_glbl(gub1)+dx)) .AND. &
                  ((D(2) .GT. y_glbl(glb2)) .AND. (D(2) .LT. y_glbl(gub2)+dy)) .AND. &
                  ((D(3) .GT. z_glbl(glb3)) .AND. (D(3) .LT. z_glbl(gub3)+dz)) 
    IF (.NOT. (A_is_inside .AND. B_is_inside .AND. C_is_inside .AND. D_is_inside)) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP, &
                      "At least one of the given vertices for defining a planar Dirichlet "//&
                      "region is outside the simulation box.", error)
    END IF

    are_coplanar = ABS(DOT_PRODUCT(vector_product(AB, AC), AD)) .LE. small_value
    IF (.NOT. are_coplanar) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP, &
            "The given vertices for defining a planar Dirichlet region are not coplanar.", error)
    END IF

    cm(1) = (A(1) + B(1) + C(1) + D(1))/4.0_dp
    cm(2) = (A(2) + B(2) + C(2) + D(2))/4.0_dp
    cm(3) = (A(3) + B(3) + C(3) + D(3))/4.0_dp
    dist1 = (A(1)-cm(1))**2 + (A(2)-cm(2))**2 + (A(3)-cm(3))**2
    dist2 = (B(1)-cm(1))**2 + (B(2)-cm(2))**2 + (B(3)-cm(3))**2
    dist3 = (C(1)-cm(1))**2 + (C(2)-cm(2))**2 + (C(3)-cm(3))**2
    dist4 = (D(1)-cm(1))**2 + (D(2)-cm(2))**2 + (D(3)-cm(3))**2
    is_rectangle = (ABS(dist1-dist2) .LE. small_value) .AND. &
                   (ABS(dist1-dist3) .LE. small_value) .AND. &
                   (ABS(dist1-dist4) .LE. small_value)  
    IF (.NOT. is_rectangle) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP, &
            "The given vertices for defining a planar Dirichlet region do not form a rectangle.", error)
    END IF

    dirichlet_bc%dbc_id = dbc_id
    dirichlet_bc%dbc_geom = PLANAR
    ALLOCATE(dirichlet_bc%rectangle)
    CALL cs_rectangle_create(dirichlet_bc%rectangle, A, B, C, D, error)
    dirichlet_bc%n_tiles = 1
    dirichlet_bc%v_D = v_D
    dirichlet_bc%smooth = smooth
    dirichlet_bc%mollifier_zeta = zeta

    IF ((unit_nr .GT. 0) .AND. debug) THEN 
       DO i = 1, 4
          WRITE(unit_nr,'(T10,A,I1,3F10.3)') "   vertex ", i, &
                                   angstrom * dirichlet_bc%rectangle%vertices(:,i)
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE create_arbitrary_planar_dbc

! *****************************************************************************
!> \brief   Creates an axis-aligned planar gate region.
!> \param x_glbl x grid vetor of the simulation box  
!> \param y_glbl y grid vetor of the simulation box  
!> \param z_glbl z grid vetor of the simulation box  
!> \param x_xtnt the x extent of the planar region 
!> \param y_xtnt the y extent of the planar region 
!> \param z_xtnt the z extent of the planar region 
!> \param v_D potential imposed at the Dirichlet part of the boundary (gate voltage)
!> \param smooth whether or not to smoothen the dirichlet region 
!> \param zeta the mollifier's width
!> \param dbc_id unique ID for the planar Dirichlet region
!> \param dirichlet_bc the dirichlet_bc object to be created 
!> \param error cp2k error
!> \par History
!>       12.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE create_aa_cuboidal_dbc(x_glbl, y_glbl, z_glbl, x_xtnt, y_xtnt, z_xtnt, &
                       v_D, smooth, zeta, dbc_id, dirichlet_bc, error)

    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl
    REAL(dp), DIMENSION(2), INTENT(IN)       :: x_xtnt, y_xtnt, z_xtnt
    REAL(dp), INTENT(IN)                     :: v_D
    LOGICAL, INTENT(IN)                      :: smooth
    REAL(dp), INTENT(IN)                     :: zeta
    INTEGER, INTENT(IN)                      :: dbc_id
    TYPE(dirichlet_bc_type), INTENT(INOUT), &
      POINTER                                :: dirichlet_bc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_aa_cuboidal_dbc', &
      routineP = moduleN//':'//routineN
    LOGICAL, DIMENSION(6), PARAMETER         :: test_forb_xtnts = .TRUE. 

    INTEGER :: handle, i, n_forb_xtnts, unit_nr, xxtnt_fst, xxtnt_ult, &
      yxtnt_fst, yxtnt_ult, zxtnt_fst, zxtnt_ult
    LOGICAL                                  :: forb_xtnt1, forb_xtnt2, &
                                                forb_xtnt3, forb_xtnt4, &
                                                forb_xtnt5, forb_xtnt6
    REAL(dp)                                 :: dx, dy, dz, xlb, xub, ylb, &
                                                yub, zlb, zub
    REAL(dp), DIMENSION(3, 8)                :: box_vtx
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
       unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
    ELSE
       unit_nr = -1
    ENDIF

    dirichlet_bc%dbc_id = dbc_id
    dirichlet_bc%dbc_geom = AA_CUBOIDAL
    dirichlet_bc%v_D = v_D
    dirichlet_bc%smooth = smooth
    dirichlet_bc%mollifier_zeta = zeta

    dx = x_glbl(LBOUND(x_glbl,1)+1) - x_glbl(LBOUND(x_glbl,1))
    dy = y_glbl(LBOUND(y_glbl,1)+1) - y_glbl(LBOUND(y_glbl,1))
    dz = z_glbl(LBOUND(z_glbl,1)+1) - z_glbl(LBOUND(z_glbl,1))

    forb_xtnt1 = x_xtnt(1) .LT. x_glbl(LBOUND(x_glbl,1))
    forb_xtnt2 = x_xtnt(2) .GT. x_glbl(UBOUND(x_glbl,1))+dx
    forb_xtnt3 = y_xtnt(1) .LT. y_glbl(LBOUND(y_glbl,1))
    forb_xtnt4 = y_xtnt(2) .GT. y_glbl(UBOUND(y_glbl,1))+dy
    forb_xtnt5 = z_xtnt(1) .LT. z_glbl(LBOUND(z_glbl,1))
    forb_xtnt6 = z_xtnt(2) .GT. z_glbl(UBOUND(z_glbl,1))+dz
    n_forb_xtnts = COUNT((/forb_xtnt1, forb_xtnt2, forb_xtnt3, &
                           forb_xtnt4, forb_xtnt5, forb_xtnt6/) .EQV. test_forb_xtnts)
    IF (n_forb_xtnts .GT. 0) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "The given extents for the Dirichlet part of the boundary are outside "//&
                      "the range of the simulation cell.", error)
    END IF 

    xxtnt_fst = find_fst_ge(x_glbl, x_xtnt(1))
    xxtnt_ult = find_lst_le(x_glbl, x_xtnt(2))
    yxtnt_fst = find_fst_ge(y_glbl, y_xtnt(1))
    yxtnt_ult = find_lst_le(y_glbl, y_xtnt(2))
    zxtnt_fst = find_fst_ge(z_glbl, z_xtnt(1))
    zxtnt_ult = find_lst_le(z_glbl, z_xtnt(2))

    xlb = x_glbl(xxtnt_fst)
    xub = x_glbl(xxtnt_ult)
    ylb = y_glbl(yxtnt_fst)
    yub = y_glbl(yxtnt_ult)
    zlb = z_glbl(zxtnt_fst)
    zub = z_glbl(zxtnt_ult)

    IF ((ABS(z_xtnt(1)-z_xtnt(2)) .LE. small_value) .OR. &
        (ABS(y_xtnt(1)-y_xtnt(2)) .LE. small_value) .OR. &
        (ABS(x_xtnt(1)-x_xtnt(2)) .LE. small_value)) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                      "According to the given extents, the gate's geometry is not a cuboid. "//&
                      "Try the input section AA_PLANAR instead.",error)
    END IF

    box_vtx(1:3,1) = (/xlb, ylb, zlb/)
    box_vtx(1:3,2) = (/xub, ylb, zlb/)
    box_vtx(1:3,3) = (/xub, yub, zlb/)
    box_vtx(1:3,4) = (/xlb, yub, zlb/)
    box_vtx(1:3,5) = (/xlb, yub, zub/)
    box_vtx(1:3,6) = (/xlb, ylb, zub/)
    box_vtx(1:3,7) = (/xub, ylb, zub/)
    box_vtx(1:3,8) = (/xub, yub, zub/)
      
    ALLOCATE(dirichlet_bc%box)
    CALL cs_box_create(dirichlet_bc%box, box_vtx, error)
    dirichlet_bc%n_tiles = 1

    IF ((unit_nr .GT. 0) .AND. debug) THEN 
       DO i = 1, 8
          WRITE(unit_nr,'(T10,A,I1,3F10.3)') "   vertex ", i, angstrom * dirichlet_bc%box%vertices(:,i)
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE create_aa_cuboidal_dbc

END MODULE dirichlet_bc_methods
