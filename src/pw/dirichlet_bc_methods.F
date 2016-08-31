!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief subroutines for defining and creating Dirichlet type subdomains
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!>       10.2015 completely revised [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
MODULE dirichlet_bc_methods

   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE dirichlet_bc_types,              ONLY: &
        AA_CUBOIDAL, AA_PLANAR, CIRCUMSCRIBED, CYLINDRICAL, INSCRIBED, PLANAR, &
        dbc_parameters_dealloc, dirichlet_bc_p_type, dirichlet_bc_type, tile_p_type, x_axis, &
        xy_plane, xz_plane, y_axis, yz_plane, z_axis
   USE erf_fn,                          ONLY: erf
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE mathlib,                         ONLY: inv_3x3,&
                                              matvec_3x3,&
                                              vector_product
   USE message_passing,                 ONLY: mp_max,&
                                              mp_sum
   USE physcon,                         ONLY: angstrom
   USE ps_implicit_types,               ONLY: MIXED_BC,&
                                              MIXED_PERIODIC_BC,&
                                              NEUMANN_BC,&
                                              PERIODIC_BC
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_zero
   USE pw_poisson_types,                ONLY: pw_poisson_parameter_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_type
   USE rs_methods,                      ONLY: setup_grid_axes
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dirichlet_bc_methods'

   PUBLIC dirichlet_boundary_region_setup

   REAL(dp), PARAMETER, PRIVATE         :: small_value = 1.0E-8_dp
   INTEGER, PARAMETER, PRIVATE          :: FWROT = 1, &
                                           BWROT = -1

CONTAINS

! **************************************************************************************************
!> \brief  Sets up the Dirichlet boundary condition
!> \param pw_pool pool of plane wave grid
!> \param poisson_params poisson_env parameters
!> \param dbcs the DBC region to be created
!> \par History
!>       10.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE dirichlet_boundary_region_setup(pw_pool, poisson_params, dbcs)

      TYPE(pw_pool_type), INTENT(IN), POINTER            :: pw_pool
      TYPE(pw_poisson_parameter_type), INTENT(INOUT)     :: poisson_params
      TYPE(dirichlet_bc_p_type), DIMENSION(:), &
         INTENT(INOUT), POINTER                          :: dbcs

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dirichlet_boundary_region_setup', &
         routineP = moduleN//':'//routineN

      INTEGER :: apx_type, dbc_id, handle, ind_end, ind_start, j, l, n_aa_cuboidal, &
         n_aa_cylindrical, n_aa_planar, n_dbcs, n_planar, parallel_axis, parallel_plane, unit_nr
      INTEGER, DIMENSION(3)                              :: n_prtn, npts
      INTEGER, DIMENSION(:), POINTER                     :: aa_cylindrical_nsides
      LOGICAL                                            :: is_periodic, verbose
      REAL(dp)                                           :: base_radius, delta_alpha, sigma, &
                                                            thickness, v_D
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: x_glbl, x_locl, y_glbl, y_locl, z_glbl, &
                                                            z_locl
      REAL(dp), DIMENSION(2)                             :: base_center
      REAL(dp), DIMENSION(2, 3)                          :: cell_xtnts
      REAL(dp), DIMENSION(3)                             :: dr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      dr = pw_pool%pw_grid%dr
      npts = pw_pool%pw_grid%npts
      cell_xtnts = 0.0_dp
      cell_xtnts(2, 1) = npts(1)*dr(1)
      cell_xtnts(2, 2) = npts(2)*dr(2)
      cell_xtnts(2, 3) = npts(3)*dr(3)

      verbose = poisson_params%dbc_params%verbose_output

      n_aa_planar = poisson_params%dbc_params%n_aa_planar
      n_aa_cuboidal = poisson_params%dbc_params%n_aa_cuboidal
      n_planar = poisson_params%dbc_params%n_planar
      n_aa_cylindrical = poisson_params%dbc_params%n_aa_cylindrical
      aa_cylindrical_nsides => poisson_params%dbc_params%aa_cylindrical_nsides
      pw_grid => pw_pool%pw_grid
      CALL setup_grid_axes(pw_grid, x_glbl, y_glbl, z_glbl, x_locl, y_locl, z_locl)
      n_dbcs = n_aa_planar+n_aa_cuboidal+n_planar+SUM(aa_cylindrical_nsides)

      SELECT CASE (poisson_params%ps_implicit_params%boundary_condition)
      CASE (MIXED_BC, MIXED_PERIODIC_BC)
         IF (n_dbcs .EQ. 0) THEN
            CPABORT("No Dirichlet region is defined.")
         END IF

         ALLOCATE (dbcs(n_dbcs))
         IF (unit_nr .GT. 0) THEN
            WRITE (unit_nr, '(/,T3,A,A,/,T3,A)') "POISSON| IMPLICIT (GENERALIZED) SOLVER ", REPEAT('-', 39), &
               "POISSON| Preparing Dirichlet regions ..."
         END IF

         DO j = 1, n_aa_planar
            ALLOCATE (dbcs(j)%dirichlet_bc)
            n_prtn = poisson_params%dbc_params%aa_planar_nprtn(:, j)
            dbc_id = AA_PLANAR+j
            v_D = poisson_params%dbc_params%aa_planar_vD(j)
            sigma = poisson_params%dbc_params%aa_planar_sigma(j)
            thickness = poisson_params%dbc_params%aa_planar_thickness(j)
            parallel_plane = poisson_params%dbc_params%aa_planar_pplane(j)
            is_periodic = poisson_params%dbc_params%aa_planar_is_periodic(j)

            IF (unit_nr .GT. 0) THEN
               WRITE (unit_nr, '(T3,A,I5)') "POISSON| Dirichlet region", j
               WRITE (unit_nr, '(T3,A)') "POISSON|     type : axis-aligned planar"
               WRITE (unit_nr, '(T3,A,E13.4,2X,A)') "POISSON|     applied potential :", v_D, "[Eh/e]"
            END IF
            CALL aa_planar_dbc_setup(cell_xtnts, parallel_plane, &
                                     poisson_params%dbc_params%aa_planar_xxtnt(:, j), &
                                     poisson_params%dbc_params%aa_planar_yxtnt(:, j), &
                                     poisson_params%dbc_params%aa_planar_zxtnt(:, j), &
                                     thickness, sigma, v_D, dbc_id, verbose, dbcs(j)%dirichlet_bc)

            CALL dirichlet_bc_partition(sigma, n_prtn, pw_pool, cell_xtnts, &
                                        x_locl, y_locl, z_locl, is_periodic, verbose, dbcs(j)%dirichlet_bc)
         END DO

         l = n_aa_planar
         DO j = l+1, l+n_aa_cuboidal
            ALLOCATE (dbcs(j)%dirichlet_bc)
            n_prtn = poisson_params%dbc_params%aa_cuboidal_nprtn(:, j-l)
            dbc_id = AA_CUBOIDAL+j-l
            v_D = poisson_params%dbc_params%aa_cuboidal_vD(j-l)
            sigma = poisson_params%dbc_params%aa_cuboidal_sigma(j-l)
            is_periodic = poisson_params%dbc_params%aa_cuboidal_is_periodic(j-l)

            IF (unit_nr .GT. 0) THEN
               WRITE (unit_nr, '(T3,A,I5)') "POISSON| Dirichlet region", j
               WRITE (unit_nr, '(T3,A)') "POISSON|     type : axis-aligned cuboidal"
               WRITE (unit_nr, '(T3,A,E13.4,2X,A)') "POISSON|     applied potential :", v_D, "[Eh/e]"
            END IF
            CALL aa_cuboidal_dbc_setup(cell_xtnts, &
                                       poisson_params%dbc_params%aa_cuboidal_xxtnt(:, j-l), &
                                       poisson_params%dbc_params%aa_cuboidal_yxtnt(:, j-l), &
                                       poisson_params%dbc_params%aa_cuboidal_zxtnt(:, j-l), &
                                       sigma, v_D, dbc_id, verbose, dbcs(j)%dirichlet_bc)

            CALL dirichlet_bc_partition(sigma, n_prtn, pw_pool, cell_xtnts, &
                                        x_locl, y_locl, z_locl, is_periodic, verbose, dbcs(j)%dirichlet_bc)
         END DO

         l = n_aa_planar+n_aa_cuboidal
         DO j = l+1, l+n_planar
            ALLOCATE (dbcs(j)%dirichlet_bc)
            n_prtn = 1
            n_prtn(1:2) = poisson_params%dbc_params%planar_nprtn(:, j-l)
            dbc_id = PLANAR+j-l
            v_D = poisson_params%dbc_params%planar_vD(j-l)
            sigma = poisson_params%dbc_params%planar_sigma(j-l)
            thickness = poisson_params%dbc_params%planar_thickness(j-l)
            is_periodic = poisson_params%dbc_params%planar_is_periodic(j-l)

            IF (unit_nr .GT. 0) THEN
               WRITE (unit_nr, '(T3,A,I5)') "POISSON| Dirichlet region", j
               WRITE (unit_nr, '(T3,A)') "POISSON|     type : planar"
               WRITE (unit_nr, '(T3,A,E13.4,2X,A)') "POISSON|     applied potential :", v_D, "[Eh/e]"
            END IF
            CALL arbitrary_planar_dbc_setup(cell_xtnts, &
                                            poisson_params%dbc_params%planar_Avtx(:, j-l), &
                                            poisson_params%dbc_params%planar_Bvtx(:, j-l), &
                                            poisson_params%dbc_params%planar_Cvtx(:, j-l), &
                                            thickness, sigma, v_D, dbc_id, verbose, dbcs(j)%dirichlet_bc)

            CALL dirichlet_bc_partition(sigma, n_prtn, pw_pool, cell_xtnts, &
                                        x_locl, y_locl, z_locl, is_periodic, verbose, dbcs(j)%dirichlet_bc)
         END DO

         l = n_aa_planar+n_aa_cuboidal+n_planar
         DO j = 1, n_aa_cylindrical
            ind_start = l+1
            ind_end = l+aa_cylindrical_nsides(j)

            n_prtn = 1
            n_prtn(1:2) = poisson_params%dbc_params%aa_cylindrical_nprtn(:, j)
            parallel_axis = poisson_params%dbc_params%aa_cylindrical_paxis(j)
            apx_type = poisson_params%dbc_params%aa_cylindrical_apxtyp(j)
            base_center = poisson_params%dbc_params%aa_cylindrical_bctr(:, j)
            base_radius = poisson_params%dbc_params%aa_cylindrical_brad(j)
            thickness = poisson_params%dbc_params%aa_cylindrical_thickness(j)
            delta_alpha = poisson_params%dbc_params%aa_cylindrical_sgap(j)
            sigma = poisson_params%dbc_params%aa_cylindrical_sigma(j)
            v_D = poisson_params%dbc_params%aa_cylindrical_vD(j)
            is_periodic = poisson_params%dbc_params%aa_cylindrical_is_periodic(j)

            IF (unit_nr .GT. 0) THEN
               WRITE (unit_nr, '(T3,A,I5)') "POISSON| Dirichlet region", l+j
               WRITE (unit_nr, '(T3,A)') "POISSON|     type : axis-aligned cylindrical"
               WRITE (unit_nr, '(T3,A,E13.4,2X,A)') "POISSON|     applied potential :", v_D, "[Eh/e]"
            END IF
            CALL aa_cylindrical_dbc_setup(pw_pool, cell_xtnts, x_locl, y_locl, z_locl, apx_type, parallel_axis, &
                                          poisson_params%dbc_params%aa_cylindrical_xtnt(:, j), &
                                          base_center, base_radius, thickness, delta_alpha, sigma, v_D, &
                                          n_prtn, is_periodic, verbose, dbcs(ind_start:ind_end))

            l = l+aa_cylindrical_nsides(j)
         END DO

      CASE (PERIODIC_BC, NEUMANN_BC)
         ! do nothing
      END SELECT

! we won't need parameters anymore so deallocate them
      CALL dbc_parameters_dealloc(poisson_params%dbc_params)

      CALL timestop(handle)

   END SUBROUTINE dirichlet_boundary_region_setup

! **************************************************************************************************
!> \brief Partitions a dirichlet_bc_type into tiles
!> \param sigma the width of the transition region between 0 and 1 at the edges of the tile
!> \param n_prtn vetor of size 3 specifying the number of times that the x, y and
!>               z interval (defining the region) should be partitioned into
!> \param pw_pool pool of plane-wave grid
!> \param cell_xtnts extents of the simulation cell
!> \param x_locl x grid vetor of the simulation box local to this process
!> \param y_locl y grid vetor of the simulation box local to this process
!> \param z_locl z grid vetor of the simulation box local to this process
!> \param is_periodic whether or not the tile is periodic with a periodicity of one unit cell
!> \param verbose whether or not to print out the coordinates of the vertices
!> \param dirichlet_bc the dirichlet_bc object to be partitioned
!> \par History
!>       10.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE dirichlet_bc_partition(sigma, n_prtn, pw_pool, cell_xtnts, &
                                     x_locl, y_locl, z_locl, is_periodic, verbose, dirichlet_bc)

      REAL(dp), INTENT(IN)                               :: sigma
      INTEGER, DIMENSION(3), INTENT(IN)                  :: n_prtn
      TYPE(pw_pool_type), INTENT(IN), POINTER            :: pw_pool
      REAL(dp), DIMENSION(2, 3), INTENT(IN)              :: cell_xtnts
      REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(IN)    :: x_locl, y_locl, z_locl
      LOGICAL, INTENT(IN)                                :: is_periodic, verbose
      TYPE(dirichlet_bc_type), INTENT(INOUT), POINTER    :: dirichlet_bc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dirichlet_bc_partition', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, k, n_tiles, unit_nr
      REAL(dp)                                           :: cyl_maxval, tile_volume
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: tile_maxvals
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(pw_type), POINTER                             :: cylinder_pw
      TYPE(tile_p_type), DIMENSION(:), POINTER           :: tiles

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      SELECT CASE (dirichlet_bc%dbc_geom)
      CASE (PLANAR)

         CALL arbitrary_dbc_partition(dirichlet_bc%vertices, n_prtn(1:2), tiles)
         n_tiles = SIZE(tiles)
         dirichlet_bc%n_tiles = n_tiles
         ALLOCATE (dirichlet_bc%tiles(n_tiles))
         dirichlet_bc%tiles = tiles
         IF ((unit_nr .GT. 0) .AND. verbose) WRITE (unit_nr, '(T7,A,I5)') "number of partitions : ", n_tiles

         ALLOCATE (tile_maxvals(n_tiles))
         tile_maxvals = 0.0_dp
         DO k = 1, n_tiles
            CALL pw_pool_create_pw(pw_pool, dirichlet_bc%tiles(k)%tile%tile_pw, use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(dirichlet_bc%tiles(k)%tile%tile_pw)

            IF ((unit_nr .GT. 0) .AND. verbose) THEN
               WRITE (unit_nr, '(T7,A,I5)') "tile", k
               DO i = 1, 8
                  WRITE (unit_nr, '(T10,A,I1,3F10.3)') &
                     "   vertex ", i, angstrom*dirichlet_bc%tiles(k)%tile%vertices(:, i)
               END DO
            END IF

            CALL arbitrary_tile_pw_compute(cell_xtnts, x_locl, y_locl, z_locl, dirichlet_bc%tiles(k)%tile%vertices, &
                                           sigma, dirichlet_bc%tiles(k)%tile%tile_pw)
            tile_volume = SUM(dirichlet_bc%tiles(k)%tile%tile_pw%cr3d)
            CALL mp_sum(tile_volume, pw_pool%pw_grid%para%group)
            dirichlet_bc%tiles(k)%tile%volume = tile_volume
         END DO

         IF ((unit_nr .GT. 0) .AND. verbose) WRITE (unit_nr, '(T3,A)') REPEAT('=', 78)

         DEALLOCATE (tiles)

      CASE (CYLINDRICAL)

         CALL arbitrary_dbc_partition(dirichlet_bc%vertices, n_prtn(1:2), tiles)
         n_tiles = SIZE(tiles)
         dirichlet_bc%n_tiles = n_tiles
         ALLOCATE (dirichlet_bc%tiles(n_tiles))
         dirichlet_bc%tiles = tiles
         IF ((unit_nr .GT. 0) .AND. verbose) WRITE (unit_nr, '(T7,A,I5)') "number of partitions : ", n_tiles

         CALL pw_pool_create_pw(pw_pool, cylinder_pw, use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(cylinder_pw)

         DO k = 1, n_tiles
            CALL pw_pool_create_pw(pw_pool, dirichlet_bc%tiles(k)%tile%tile_pw, use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(dirichlet_bc%tiles(k)%tile%tile_pw)

            IF ((unit_nr .GT. 0) .AND. verbose) THEN
               WRITE (unit_nr, '(T7,A,I5)') "tile", k
               DO i = 1, 8
                  WRITE (unit_nr, '(T10,A,I1,3F10.3)') &
                     "   vertex ", i, angstrom*dirichlet_bc%tiles(k)%tile%vertices(:, i)
               END DO
            END IF

            CALL arbitrary_tile_pw_compute(cell_xtnts, x_locl, y_locl, z_locl, dirichlet_bc%tiles(k)%tile%vertices, &
                                           sigma, dirichlet_bc%tiles(k)%tile%tile_pw)
            CALL pw_axpy(dirichlet_bc%tiles(k)%tile%tile_pw, cylinder_pw)
            tile_volume = SUM(dirichlet_bc%tiles(k)%tile%tile_pw%cr3d)
            CALL mp_sum(tile_volume, pw_pool%pw_grid%para%group)
            dirichlet_bc%tiles(k)%tile%volume = tile_volume
         END DO

         cyl_maxval = MAXVAL(cylinder_pw%cr3d)
         CALL mp_max(cyl_maxval, pw_pool%pw_grid%para%group)
         IF (cyl_maxval .GT. 1.01_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "Inaccurate approximation of unit step function at cylindrical Dirichlet region. "// &
                          "Increase the number of sides (N_SIDES) and/or the base radius. Alternatively, "// &
                          "one can increase DELTA_ALPHA (see the reference manual).")
         END IF

         CALL pw_pool_give_back_pw(pw_pool, cylinder_pw)

         IF ((unit_nr .GT. 0) .AND. verbose) WRITE (unit_nr, '(T3,A)') REPEAT('=', 78)

         DEALLOCATE (tiles)

      CASE (AA_PLANAR, AA_CUBOIDAL)

         CALL aa_dbc_partition(dirichlet_bc%vertices, n_prtn, tiles)
         n_tiles = SIZE(tiles)
         dirichlet_bc%n_tiles = n_tiles
         ALLOCATE (dirichlet_bc%tiles(n_tiles))
         dirichlet_bc%tiles = tiles
         IF ((unit_nr .GT. 0) .AND. verbose) WRITE (unit_nr, '(T7,A,I5)') "number of partitions : ", n_tiles

         DO k = 1, n_tiles
            CALL pw_pool_create_pw(pw_pool, dirichlet_bc%tiles(k)%tile%tile_pw, use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(dirichlet_bc%tiles(k)%tile%tile_pw)

            IF ((unit_nr .GT. 0) .AND. verbose) THEN
               WRITE (unit_nr, '(T7,A,I5)') "tile", k
               DO i = 1, 8
                  WRITE (unit_nr, '(T10,A,I1,3F10.3)') &
                     "   vertex ", i, angstrom*dirichlet_bc%tiles(k)%tile%vertices(:, i)
               END DO
            END IF

            CALL aa_tile_pw_compute(cell_xtnts, x_locl, y_locl, z_locl, dirichlet_bc%tiles(k)%tile%vertices, &
                                    sigma, is_periodic, dirichlet_bc%tiles(k)%tile%tile_pw)
            tile_volume = SUM(dirichlet_bc%tiles(k)%tile%tile_pw%cr3d)
            CALL mp_sum(tile_volume, pw_pool%pw_grid%para%group)
            dirichlet_bc%tiles(k)%tile%volume = tile_volume
         END DO

         IF ((unit_nr .GT. 0) .AND. verbose) WRITE (unit_nr, '(T3,A)') REPEAT('=', 78)

         DEALLOCATE (tiles)

      END SELECT

      CALL timestop(handle)

   END SUBROUTINE dirichlet_bc_partition

! **************************************************************************************************
!> \brief   Creates an axis-aligned planar Dirichlet region.
!> \param cell_xtnts extents of the simulation cell
!> \param parallel_plane the coordinate plane that the region is parallel to
!> \param x_xtnt the x extent of the planar region
!> \param y_xtnt the y extent of the planar region
!> \param z_xtnt the z extent of the planar region
!> \param thickness the thickness of the region
!> \param sigma the width of the transition region between 0 and 1 at the edges of the tile
!> \param v_D potential imposed at the Dirichlet part of the boundary (gate voltage)
!> \param dbc_id unique ID for the Dirichlet region
!> \param verbose whether or not to print out the coordinates of the vertices
!> \param dirichlet_bc the dirichlet_bc object to be created
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE aa_planar_dbc_setup(cell_xtnts, parallel_plane, x_xtnt, y_xtnt, z_xtnt, &
                                  thickness, sigma, v_D, dbc_id, verbose, dirichlet_bc)

      REAL(dp), DIMENSION(2, 3), INTENT(IN)              :: cell_xtnts
      INTEGER, INTENT(IN)                                :: parallel_plane
      REAL(dp), DIMENSION(2), INTENT(IN)                 :: x_xtnt, y_xtnt, z_xtnt
      REAL(dp), INTENT(IN)                               :: thickness, sigma, v_D
      INTEGER, INTENT(IN)                                :: dbc_id
      LOGICAL, INTENT(IN)                                :: verbose
      TYPE(dirichlet_bc_type), INTENT(INOUT), POINTER    :: dirichlet_bc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'aa_planar_dbc_setup', &
         routineP = moduleN//':'//routineN
      LOGICAL, DIMENSION(6), PARAMETER                   :: test_forb_xtnts = .TRUE.

      INTEGER                                            :: handle, i, n_forb_xtnts, unit_nr
      LOGICAL                                            :: forb_xtnt1, forb_xtnt2, forb_xtnt3, &
                                                            forb_xtnt4, forb_xtnt5, forb_xtnt6
      REAL(dp)                                           :: xlb, xub, ylb, yub, zlb, zub
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      xlb = x_xtnt(1); xub = x_xtnt(2)
      ylb = y_xtnt(1); yub = y_xtnt(2)
      zlb = z_xtnt(1); zub = z_xtnt(2)

      SELECT CASE (parallel_plane)
      CASE (xy_plane)
         zlb = zlb-thickness*0.5_dp
         zub = zub+thickness*0.5_dp
      CASE (xz_plane)
         ylb = ylb-thickness*0.5_dp
         yub = yub+thickness*0.5_dp
      CASE (yz_plane)
         xlb = xlb-thickness*0.5_dp
         xub = xub+thickness*0.5_dp
      END SELECT

      forb_xtnt1 = xlb .LT. cell_xtnts(1, 1)
      forb_xtnt2 = xub .GT. cell_xtnts(2, 1)
      forb_xtnt3 = ylb .LT. cell_xtnts(1, 2)
      forb_xtnt4 = yub .GT. cell_xtnts(2, 2)
      forb_xtnt5 = zlb .LT. cell_xtnts(1, 3)
      forb_xtnt6 = zub .GT. cell_xtnts(2, 3)
      n_forb_xtnts = COUNT((/forb_xtnt1, forb_xtnt2, forb_xtnt3, &
                             forb_xtnt4, forb_xtnt5, forb_xtnt6/) .EQV. test_forb_xtnts)
      IF (n_forb_xtnts .GT. 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "The given extents for the Dirichlet region are outside the range of "// &
                       "the simulation cell.")
      END IF

      dirichlet_bc%v_D = v_D
      dirichlet_bc%dbc_id = dbc_id
      dirichlet_bc%dbc_geom = AA_PLANAR
      dirichlet_bc%smoothing_width = sigma
      dirichlet_bc%n_tiles = 1

      dirichlet_bc%vertices(1:3, 1) = (/xlb, ylb, zlb/)
      dirichlet_bc%vertices(1:3, 2) = (/xub, ylb, zlb/)
      dirichlet_bc%vertices(1:3, 3) = (/xub, yub, zlb/)
      dirichlet_bc%vertices(1:3, 4) = (/xlb, yub, zlb/)
      dirichlet_bc%vertices(1:3, 5) = (/xlb, yub, zub/)
      dirichlet_bc%vertices(1:3, 6) = (/xlb, ylb, zub/)
      dirichlet_bc%vertices(1:3, 7) = (/xub, ylb, zub/)
      dirichlet_bc%vertices(1:3, 8) = (/xub, yub, zub/)

      IF ((unit_nr .GT. 0) .AND. verbose) THEN
         WRITE (unit_nr, '(T3,A,A)') "======== verbose ", REPEAT('=', 61)
         DO i = 1, 8
            WRITE (unit_nr, '(T10,A,I1,3F10.3)') "   vertex ", i, angstrom*dirichlet_bc%vertices(:, i)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE aa_planar_dbc_setup

! **************************************************************************************************
!> \brief   Creates an axis-aligned cuboidal Dirichlet region.
!> \param cell_xtnts extents of the simulation cell
!> \param x_xtnt the x extent of the cuboidal region
!> \param y_xtnt the y extent of the cuboidal region
!> \param z_xtnt the z extent of the cuboidal region
!> \param sigma the width of the transition region between 0 and 1 at the edges of the tile
!> \param v_D potential imposed at the Dirichlet part of the boundary (gate voltage)
!> \param dbc_id unique ID for the Dirichlet region
!> \param verbose whether or not to print out the coordinates of the vertices
!> \param dirichlet_bc the dirichlet_bc object to be created
!> \par History
!>       12.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE aa_cuboidal_dbc_setup(cell_xtnts, x_xtnt, y_xtnt, z_xtnt, &
                                    sigma, v_D, dbc_id, verbose, dirichlet_bc)

      REAL(dp), DIMENSION(2, 3), INTENT(IN)              :: cell_xtnts
      REAL(dp), DIMENSION(2), INTENT(IN)                 :: x_xtnt, y_xtnt, z_xtnt
      REAL(dp), INTENT(IN)                               :: sigma, v_D
      INTEGER, INTENT(IN)                                :: dbc_id
      LOGICAL, INTENT(IN)                                :: verbose
      TYPE(dirichlet_bc_type), INTENT(INOUT), POINTER    :: dirichlet_bc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'aa_cuboidal_dbc_setup', &
         routineP = moduleN//':'//routineN
      LOGICAL, DIMENSION(6), PARAMETER                   :: test_forb_xtnts = .TRUE.

      INTEGER                                            :: handle, i, n_forb_xtnts, unit_nr
      LOGICAL                                            :: forb_xtnt1, forb_xtnt2, forb_xtnt3, &
                                                            forb_xtnt4, forb_xtnt5, forb_xtnt6
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      forb_xtnt1 = x_xtnt(1) .LT. cell_xtnts(1, 1)
      forb_xtnt2 = x_xtnt(2) .GT. cell_xtnts(2, 1)
      forb_xtnt3 = y_xtnt(1) .LT. cell_xtnts(1, 2)
      forb_xtnt4 = y_xtnt(2) .GT. cell_xtnts(2, 2)
      forb_xtnt5 = z_xtnt(1) .LT. cell_xtnts(1, 3)
      forb_xtnt6 = z_xtnt(2) .GT. cell_xtnts(2, 3)
      n_forb_xtnts = COUNT((/forb_xtnt1, forb_xtnt2, forb_xtnt3, &
                             forb_xtnt4, forb_xtnt5, forb_xtnt6/) .EQV. test_forb_xtnts)
      IF (n_forb_xtnts .GT. 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "The given extents for the Dirichlet region are outside the range of "// &
                       "the simulation cell.")
      END IF

      dirichlet_bc%v_D = v_D
      dirichlet_bc%dbc_id = dbc_id
      dirichlet_bc%dbc_geom = AA_CUBOIDAL
      dirichlet_bc%smoothing_width = sigma
      dirichlet_bc%n_tiles = 1

      dirichlet_bc%vertices(1:3, 1) = (/x_xtnt(1), y_xtnt(1), z_xtnt(1)/)
      dirichlet_bc%vertices(1:3, 2) = (/x_xtnt(2), y_xtnt(1), z_xtnt(1)/)
      dirichlet_bc%vertices(1:3, 3) = (/x_xtnt(2), y_xtnt(2), z_xtnt(1)/)
      dirichlet_bc%vertices(1:3, 4) = (/x_xtnt(1), y_xtnt(2), z_xtnt(1)/)
      dirichlet_bc%vertices(1:3, 5) = (/x_xtnt(1), y_xtnt(2), z_xtnt(2)/)
      dirichlet_bc%vertices(1:3, 6) = (/x_xtnt(1), y_xtnt(1), z_xtnt(2)/)
      dirichlet_bc%vertices(1:3, 7) = (/x_xtnt(2), y_xtnt(1), z_xtnt(2)/)
      dirichlet_bc%vertices(1:3, 8) = (/x_xtnt(2), y_xtnt(2), z_xtnt(2)/)

      IF ((unit_nr .GT. 0) .AND. verbose) THEN
         WRITE (unit_nr, '(T3,A,A)') "======== verbose ", REPEAT('=', 61)
         DO i = 1, 8
            WRITE (unit_nr, '(T10,A,I1,3F10.3)') "   vertex ", i, angstrom*dirichlet_bc%vertices(:, i)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE aa_cuboidal_dbc_setup

! **************************************************************************************************
!> \brief   Creates an arbitrary (rectangular-shaped) planar Dirichlet region given
!>     the coordinates of its vertices.
!> \param cell_xtnts extents of the simulation cell
!> \param A coordinates of the vertex A
!> \param B coordinates of the vertex B
!> \param C coordinates of the vertex C
!> \param thickness the thickness of the region
!> \param sigma the width of the transition region between 0 and 1 at the edges of the tile
!> \param v_D potential imposed at the Dirichlet part of the boundary (gate voltage)
!> \param dbc_id unique ID for the Dirichlet region
!> \param verbose whether or not to print out the coordinates of the vertices
!> \param dirichlet_bc the dirichlet_bc object to be created
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE arbitrary_planar_dbc_setup(cell_xtnts, A, B, C, thickness, &
                                         sigma, v_D, dbc_id, verbose, dirichlet_bc)

      REAL(dp), DIMENSION(2, 3), INTENT(IN)              :: cell_xtnts
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: A, B, C
      REAL(dp), INTENT(IN)                               :: thickness, sigma, v_D
      INTEGER, INTENT(IN)                                :: dbc_id
      LOGICAL, INTENT(IN)                                :: verbose
      TYPE(dirichlet_bc_type), INTENT(INOUT), POINTER    :: dirichlet_bc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'arbitrary_planar_dbc_setup', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, unit_nr
      LOGICAL                                            :: A_is_inside, are_coplanar, &
                                                            are_orthogonal, B_is_inside, &
                                                            C_is_inside, D_is_inside, is_rectangle
      REAL(dp)                                           :: dist1, dist2, dist3, dist4
      REAL(dp), DIMENSION(3)                             :: AB, AC, AD, BC, cm, D, normal_vector, &
                                                            unit_normal
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      D = A+(C-B)
      AB = B-A; AC = C-A; AD = D-A; BC = C-B
      are_orthogonal = ABS(DOT_PRODUCT(AB, BC)/(SQRT(SUM(AB**2))*SQRT(SUM(BC**2)))) .LE. small_value
      IF (.NOT. are_orthogonal) THEN
         CALL cp_abort(__LOCATION__, "The given vertices for defining a "// &
                       "planar Dirichlet region do not form orthogonal edges.")
      END IF

      normal_vector = vector_product(AB, AC)
      unit_normal = normal_vector/SQRT(SUM(normal_vector**2))

      A_is_inside = (A(1) .GT. cell_xtnts(1, 1) .AND. A(1) .LT. cell_xtnts(2, 1)) .AND. &
                    (A(2) .GT. cell_xtnts(1, 2) .AND. A(2) .LT. cell_xtnts(2, 2)) .AND. &
                    (A(3) .GT. cell_xtnts(1, 3) .AND. A(3) .LT. cell_xtnts(2, 3))
      B_is_inside = (B(1) .GT. cell_xtnts(1, 1) .AND. B(1) .LT. cell_xtnts(2, 1)) .AND. &
                    (B(2) .GT. cell_xtnts(1, 2) .AND. B(2) .LT. cell_xtnts(2, 2)) .AND. &
                    (B(3) .GT. cell_xtnts(1, 3) .AND. B(3) .LT. cell_xtnts(2, 3))
      C_is_inside = (C(1) .GT. cell_xtnts(1, 1) .AND. C(1) .LT. cell_xtnts(2, 1)) .AND. &
                    (C(2) .GT. cell_xtnts(1, 2) .AND. C(2) .LT. cell_xtnts(2, 2)) .AND. &
                    (C(3) .GT. cell_xtnts(1, 3) .AND. C(3) .LT. cell_xtnts(2, 3))
      D_is_inside = (D(1) .GT. cell_xtnts(1, 1) .AND. D(1) .LT. cell_xtnts(2, 1)) .AND. &
                    (D(2) .GT. cell_xtnts(1, 2) .AND. D(2) .LT. cell_xtnts(2, 2)) .AND. &
                    (D(3) .GT. cell_xtnts(1, 3) .AND. D(3) .LT. cell_xtnts(2, 3))
      IF (.NOT. (A_is_inside .AND. B_is_inside .AND. C_is_inside .AND. D_is_inside)) THEN
         CALL cp_abort(__LOCATION__, &
                       "At least one of the given vertices for defining a planar Dirichlet "// &
                       "region is outside the simulation box.")
      END IF

      are_coplanar = ABS(DOT_PRODUCT(vector_product(AB, AC), AD)) .LE. small_value
      IF (.NOT. are_coplanar) THEN
         CPABORT("The given vertices for defining a planar Dirichlet region are not coplanar.")
      END IF

      cm = (A+B+C+D)/4.0_dp
      dist1 = SUM((A-cm)**2)
      dist2 = SUM((B-cm)**2)
      dist3 = SUM((C-cm)**2)
      dist4 = SUM((D-cm)**2)
      is_rectangle = (ABS(dist1-dist2) .LE. small_value) .AND. &
                     (ABS(dist1-dist3) .LE. small_value) .AND. &
                     (ABS(dist1-dist4) .LE. small_value)
      IF (.NOT. is_rectangle) THEN
         CALL cp_abort(__LOCATION__, "The given vertices for defining "// &
                       "a planar Dirichlet region do not form a rectangle.")
      END IF

      dirichlet_bc%v_D = v_D
      dirichlet_bc%dbc_id = dbc_id
      dirichlet_bc%dbc_geom = PLANAR
      dirichlet_bc%smoothing_width = sigma
      dirichlet_bc%n_tiles = 1

      dirichlet_bc%vertices(1:3, 1) = A-0.5_dp*thickness*unit_normal
      dirichlet_bc%vertices(1:3, 2) = B-0.5_dp*thickness*unit_normal
      dirichlet_bc%vertices(1:3, 3) = C-0.5_dp*thickness*unit_normal
      dirichlet_bc%vertices(1:3, 4) = D-0.5_dp*thickness*unit_normal
      dirichlet_bc%vertices(1:3, 5) = D+0.5_dp*thickness*unit_normal
      dirichlet_bc%vertices(1:3, 6) = A+0.5_dp*thickness*unit_normal
      dirichlet_bc%vertices(1:3, 7) = B+0.5_dp*thickness*unit_normal
      dirichlet_bc%vertices(1:3, 8) = C+0.5_dp*thickness*unit_normal

      IF ((unit_nr .GT. 0) .AND. verbose) THEN
         WRITE (unit_nr, '(T3,A,A)') "======== verbose ", REPEAT('=', 61)
         DO i = 1, 8
            WRITE (unit_nr, '(T10,A,I1,3F10.3)') "   vertex ", i, angstrom*dirichlet_bc%vertices(:, i)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE arbitrary_planar_dbc_setup

! **************************************************************************************************
!> \brief   Creates an x-axis-aligned cylindrical Dirichlet region. The cylindrical
!>     shape is approximated by an n-gonal (circumscribed or inscribed) uniform prism-
!>     like object. The circular base is approximated by a polygon. A second polygon is
!>     then created by rotating the first one by delta_alpha radians. Let's name the
!>     vertices of the first and the second polygon as v1, v2, v3, ... vn-1, vn and
!>     w1, w2, w3, ... wn-1, wn, respectively. Then w1v2, w2v3, ... wn-1vn, wnv1 form
!>     the edges of the cross-section of the approximating object. This way the dbcs
!>     do not share any edge.
!> \param pw_pool  pool of plane wave grid
!> \param cell_xtnts extents of the simulation cell
!> \param x_locl x grid vetor of the simulation box local to this process
!> \param y_locl y grid vetor of the simulation box local to this process
!> \param z_locl z grid vetor of the simulation box local to this process
!> \param apx_type the type of the n-gonal prism approximating the cylinder
!> \param parallel_axis the coordinate axis that the cylindrical region extends along
!> \param xtnt the x extent of the cylindrical region along its central axis
!> \param base_center the y and z coordinates of the cylinder's base
!> \param base_radius the radius of the cylinder's base
!> \param thickness the thickness of the region
!> \param delta_alpha ...
!> \param sigma the width of the transition region between 0 and 1 at the edges of the tile
!> \param v_D potential imposed at the Dirichlet part of the boundary (gate voltage)
!> \param n_prtn number of partitions along the edges of the faces of the prism
!> \param is_periodic whether or not the tile is periodic with a periodicity of one unit cell
!> \param verbose whether or not to print out the coordinates of the vertices
!> \param dbcs  the x-axis-aligned cylindrical gate region to be created
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE aa_cylindrical_dbc_setup(pw_pool, cell_xtnts, x_locl, y_locl, z_locl, apx_type, parallel_axis, &
                                       xtnt, base_center, base_radius, thickness, delta_alpha, &
                                       sigma, v_D, n_prtn, is_periodic, verbose, dbcs)

      TYPE(pw_pool_type), POINTER                        :: pw_pool
      REAL(dp), DIMENSION(2, 3), INTENT(IN)              :: cell_xtnts
      REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(IN)    :: x_locl, y_locl, z_locl
      INTEGER, INTENT(IN), OPTIONAL                      :: apx_type
      INTEGER, INTENT(IN)                                :: parallel_axis
      REAL(dp), DIMENSION(2), INTENT(IN)                 :: xtnt, base_center
      REAL(dp), INTENT(IN)                               :: base_radius, thickness, delta_alpha, &
                                                            sigma, v_D
      INTEGER, DIMENSION(3), INTENT(IN)                  :: n_prtn
      LOGICAL, INTENT(IN)                                :: is_periodic, verbose
      TYPE(dirichlet_bc_p_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: dbcs

      CHARACTER(LEN=*), PARAMETER :: routineN = 'aa_cylindrical_dbc_setup', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, intern_apx_type, j, n_dbcs, &
                                                            unit_nr
      LOGICAL                                            :: forb_xtnt1, forb_xtnt2, forb_xtnt3, &
                                                            forb_xtnt4
      REAL(dp)                                           :: cntrlaxis_lb, cntrlaxis_ub, h, Lx, Ly, &
                                                            Lz, theta
      REAL(dp), DIMENSION(3)                             :: A, B, C, D, normal_vec, unit_normal
      TYPE(cp_logger_type), POINTER                      :: logger
      REAL(dp), DIMENSION(SIZE(dbcs)+1)                  :: alpha, alpha_rotd, xlb, xub, ylb, yub, &
                                                            zlb, zub

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos .EQ. logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      Lx = cell_xtnts(2, 1)-cell_xtnts(1, 1)
      Ly = cell_xtnts(2, 2)-cell_xtnts(1, 2)
      Lz = cell_xtnts(2, 3)-cell_xtnts(1, 3)

      SELECT CASE (parallel_axis)
      CASE (x_axis)
         IF (xtnt(1) .LT. cell_xtnts(1, 1) .OR. xtnt(2) .GT. cell_xtnts(2, 1)) THEN
            CALL cp_abort(__LOCATION__, &
                          "The length of the cylindrical Dirichlet region is larger than the "// &
                          "x range of the simulation cell.")
         END IF
         forb_xtnt1 = base_center(1)-base_radius .LT. cell_xtnts(1, 2)
         forb_xtnt2 = base_center(1)+base_radius .GT. cell_xtnts(2, 2)
         forb_xtnt3 = base_center(2)-base_radius .LT. cell_xtnts(1, 3)
         forb_xtnt4 = base_center(2)+base_radius .GT. cell_xtnts(2, 3)
         IF (forb_xtnt1 .OR. forb_xtnt2 .OR. forb_xtnt3 .OR. forb_xtnt4) THEN
            CPABORT("The cylinder does not fit entirely inside the simulation cell.")
         END IF
      CASE (y_axis)
         IF (xtnt(1) .LT. cell_xtnts(1, 2) .OR. xtnt(2) .GT. cell_xtnts(2, 2)) THEN
            CALL cp_abort(__LOCATION__, &
                          "The length of the cylindrical Dirichlet region is larger than the "// &
                          "y range of the simulation cell.")
         END IF
         forb_xtnt1 = base_center(1)-base_radius .LT. cell_xtnts(1, 1)
         forb_xtnt2 = base_center(1)+base_radius .GT. cell_xtnts(2, 1)
         forb_xtnt3 = base_center(2)-base_radius .LT. cell_xtnts(1, 3)
         forb_xtnt4 = base_center(2)+base_radius .GT. cell_xtnts(2, 3)
         IF (forb_xtnt1 .OR. forb_xtnt2 .OR. forb_xtnt3 .OR. forb_xtnt4) THEN
            CPABORT("The cylinder does not fit entirely inside the simulation cell.")
         END IF
      CASE (z_axis)
         IF (xtnt(1) .LT. cell_xtnts(1, 3) .OR. xtnt(2) .GT. cell_xtnts(2, 3)) THEN
            CALL cp_abort(__LOCATION__, &
                          "The length of the cylindrical Dirichlet region is larger than the "// &
                          "z range of the simulation cell.")
         END IF
         forb_xtnt1 = base_center(1)-base_radius .LT. cell_xtnts(1, 1)
         forb_xtnt2 = base_center(1)+base_radius .GT. cell_xtnts(2, 1)
         forb_xtnt3 = base_center(2)-base_radius .LT. cell_xtnts(1, 2)
         forb_xtnt4 = base_center(2)+base_radius .GT. cell_xtnts(2, 2)
         IF (forb_xtnt1 .OR. forb_xtnt2 .OR. forb_xtnt3 .OR. forb_xtnt4) THEN
            CPABORT("The cylinder does not fit entirely inside the simulation cell.")
         END IF
      END SELECT

      intern_apx_type = CIRCUMSCRIBED
      IF (PRESENT(apx_type)) intern_apx_type = apx_type

      n_dbcs = SIZE(dbcs)

      cntrlaxis_lb = xtnt(1)
      cntrlaxis_ub = xtnt(2)

      theta = twopi/n_dbcs

      IF (intern_apx_type .EQ. INSCRIBED) THEN
         h = base_radius ! inscribed uniform prism
      ELSE IF (intern_apx_type .EQ. CIRCUMSCRIBED) THEN
         h = base_radius/COS(0.5*theta) ! circumscribed uniform prism
      ELSE
         CPABORT("Unknown approximation type for cylinder.")
      END IF
      SELECT CASE (parallel_axis)
      CASE (x_axis)
         IF (h .GT. MINVAL((/Ly, Lz/))) THEN
            CPABORT("Reduce the base radius!")
         END IF
      CASE (y_axis)
         IF (h .GT. MINVAL((/Lx, Lz/))) THEN
            CPABORT("Reduce the base radius!")
         END IF
      CASE (z_axis)
         IF (h .GT. MINVAL((/Lx, Ly/))) THEN
            CPABORT("Reduce the base radius!")
         END IF
      END SELECT

      alpha = linspace(0.0_dp, 2*pi, n_dbcs+1)
      alpha_rotd = alpha+delta_alpha; 
      SELECT CASE (parallel_axis)
      CASE (x_axis)
         DO j = 1, n_dbcs
            ylb(j) = base_center(1)+h*SIN(alpha(j))
            zlb(j) = base_center(2)+h*COS(alpha(j))
            yub(j) = base_center(1)+h*SIN(alpha_rotd(j))
            zub(j) = base_center(2)+h*COS(alpha_rotd(j))
         END DO
         ylb(n_dbcs+1) = ylb(1)
         yub(n_dbcs+1) = yub(1)
         zlb(n_dbcs+1) = zlb(1)
         zub(n_dbcs+1) = zub(1)
         DO j = 1, n_dbcs
            ALLOCATE (dbcs(j)%dirichlet_bc)
            dbcs(j)%dirichlet_bc%dbc_geom = CYLINDRICAL
            dbcs(j)%dirichlet_bc%v_D = v_D
            dbcs(j)%dirichlet_bc%dbc_id = CYLINDRICAL+j
            dbcs(j)%dirichlet_bc%smoothing_width = sigma

            A = (/cntrlaxis_lb, yub(j), zub(j)/)
            B = (/cntrlaxis_lb, ylb(j+1), zlb(j+1)/)
            C = (/cntrlaxis_ub, ylb(j+1), zlb(j+1)/)
            D = (/cntrlaxis_ub, yub(j), zub(j)/)
            normal_vec = vector_product((A-C), (D-B))
            unit_normal = normal_vec/SQRT(SUM(normal_vec**2))

            dbcs(j)%dirichlet_bc%vertices(1:3, 1) = A
            dbcs(j)%dirichlet_bc%vertices(1:3, 2) = B
            dbcs(j)%dirichlet_bc%vertices(1:3, 3) = C
            dbcs(j)%dirichlet_bc%vertices(1:3, 4) = D
            dbcs(j)%dirichlet_bc%vertices(1:3, 5) = D+thickness*unit_normal
            dbcs(j)%dirichlet_bc%vertices(1:3, 6) = A+thickness*unit_normal
            dbcs(j)%dirichlet_bc%vertices(1:3, 7) = B+thickness*unit_normal
            dbcs(j)%dirichlet_bc%vertices(1:3, 8) = C+thickness*unit_normal

            dbcs(j)%dirichlet_bc%n_tiles = 1

            IF ((unit_nr .GT. 0) .AND. verbose) THEN
               WRITE (unit_nr, '(T3,A,A)') "======== verbose ", REPEAT('=', 61)
               WRITE (unit_nr, '(T7,A,I5,T20,A,I5,A)') "edge", j, "of the", n_dbcs, &
                  "-gonal prism approximating the cylinder"
               DO i = 1, 8
                  WRITE (unit_nr, '(T10,A,I1,3F10.3)') "   vertex ", i, &
                     angstrom*dbcs(j)%dirichlet_bc%vertices(:, i)
               END DO
            END IF

            CALL dirichlet_bc_partition(sigma, n_prtn, pw_pool, cell_xtnts, &
                                        x_locl, y_locl, z_locl, is_periodic, verbose, dbcs(j)%dirichlet_bc)
         END DO
      CASE (y_axis)
         DO j = 1, n_dbcs
            xlb(j) = base_center(1)+h*SIN(alpha(j))
            zlb(j) = base_center(2)+h*COS(alpha(j))
            xub(j) = base_center(1)+h*SIN(alpha_rotd(j))
            zub(j) = base_center(2)+h*COS(alpha_rotd(j))
         END DO
         xlb(n_dbcs+1) = xlb(1)
         xub(n_dbcs+1) = xub(1)
         zlb(n_dbcs+1) = zlb(1)
         zub(n_dbcs+1) = zub(1)
         DO j = 1, n_dbcs
            ALLOCATE (dbcs(j)%dirichlet_bc)
            dbcs(j)%dirichlet_bc%dbc_geom = CYLINDRICAL
            dbcs(j)%dirichlet_bc%v_D = v_D
            dbcs(j)%dirichlet_bc%dbc_id = CYLINDRICAL+j
            dbcs(j)%dirichlet_bc%smoothing_width = sigma

            A = (/xub(j), cntrlaxis_lb, zub(j)/)
            B = (/xlb(j+1), cntrlaxis_lb, zlb(j+1)/)
            C = (/xlb(j+1), cntrlaxis_ub, zlb(j+1)/)
            D = (/xub(j), cntrlaxis_ub, zub(j)/)
            normal_vec = vector_product((A-C), (D-B))
            unit_normal = normal_vec/SQRT(SUM(normal_vec**2))

            dbcs(j)%dirichlet_bc%vertices(1:3, 1) = A
            dbcs(j)%dirichlet_bc%vertices(1:3, 2) = B
            dbcs(j)%dirichlet_bc%vertices(1:3, 3) = C
            dbcs(j)%dirichlet_bc%vertices(1:3, 4) = D
            dbcs(j)%dirichlet_bc%vertices(1:3, 5) = D+thickness*unit_normal
            dbcs(j)%dirichlet_bc%vertices(1:3, 6) = A+thickness*unit_normal
            dbcs(j)%dirichlet_bc%vertices(1:3, 7) = B+thickness*unit_normal
            dbcs(j)%dirichlet_bc%vertices(1:3, 8) = C+thickness*unit_normal

            dbcs(j)%dirichlet_bc%n_tiles = 1

            IF ((unit_nr .GT. 0) .AND. verbose) THEN
               WRITE (unit_nr, '(T3,A,A)') "======== verbose ", REPEAT('=', 61)
               WRITE (unit_nr, '(T7,A,I5,T20,A,I5,A)') "edge", j, "of the", n_dbcs, &
                  "-gonal prism approximating the cylinder"
               DO i = 1, 8
                  WRITE (unit_nr, '(T10,A,I1,3F10.3)') "   vertex ", i, &
                     angstrom*dbcs(j)%dirichlet_bc%vertices(:, i)
               END DO
            END IF

            CALL dirichlet_bc_partition(sigma, n_prtn, pw_pool, cell_xtnts, &
                                        x_locl, y_locl, z_locl, is_periodic, verbose, dbcs(j)%dirichlet_bc)
         END DO
      CASE (z_axis)
         DO j = 1, n_dbcs
            xlb(j) = base_center(1)+h*SIN(alpha(j))
            ylb(j) = base_center(2)+h*COS(alpha(j))
            xub(j) = base_center(1)+h*SIN(alpha_rotd(j))
            yub(j) = base_center(2)+h*COS(alpha_rotd(j))
         END DO
         xlb(n_dbcs+1) = xlb(1)
         xub(n_dbcs+1) = xub(1)
         ylb(n_dbcs+1) = ylb(1)
         yub(n_dbcs+1) = yub(1)
         DO j = 1, n_dbcs
            ALLOCATE (dbcs(j)%dirichlet_bc)
            dbcs(j)%dirichlet_bc%dbc_geom = CYLINDRICAL
            dbcs(j)%dirichlet_bc%v_D = v_D
            dbcs(j)%dirichlet_bc%dbc_id = CYLINDRICAL+j
            dbcs(j)%dirichlet_bc%smoothing_width = sigma

            A = (/xub(j), yub(j), cntrlaxis_lb/)
            B = (/xlb(j+1), ylb(j+1), cntrlaxis_lb/)
            C = (/xlb(j+1), ylb(j+1), cntrlaxis_ub/)
            D = (/xub(j), yub(j), cntrlaxis_ub/)
            normal_vec = vector_product((A-C), (D-B))
            unit_normal = normal_vec/SQRT(SUM(normal_vec**2))

            dbcs(j)%dirichlet_bc%vertices(1:3, 1) = A
            dbcs(j)%dirichlet_bc%vertices(1:3, 2) = B
            dbcs(j)%dirichlet_bc%vertices(1:3, 3) = C
            dbcs(j)%dirichlet_bc%vertices(1:3, 4) = D
            dbcs(j)%dirichlet_bc%vertices(1:3, 5) = D+thickness*unit_normal
            dbcs(j)%dirichlet_bc%vertices(1:3, 6) = A+thickness*unit_normal
            dbcs(j)%dirichlet_bc%vertices(1:3, 7) = B+thickness*unit_normal
            dbcs(j)%dirichlet_bc%vertices(1:3, 8) = C+thickness*unit_normal

            dbcs(j)%dirichlet_bc%n_tiles = 1

            IF ((unit_nr .GT. 0) .AND. verbose) THEN
               WRITE (unit_nr, '(T3,A,A)') "======== verbose ", REPEAT('=', 61)
               WRITE (unit_nr, '(T7,A,I5,T20,A,I5,A)') "edge", j, "of the", n_dbcs, &
                  "-gonal prism approximating the cylinder"
               DO i = 1, 8
                  WRITE (unit_nr, '(T10,A,I1,3F10.3)') "   vertex ", i, &
                     angstrom*dbcs(j)%dirichlet_bc%vertices(:, i)
               END DO
            END IF

            CALL dirichlet_bc_partition(sigma, n_prtn, pw_pool, cell_xtnts, &
                                        x_locl, y_locl, z_locl, is_periodic, verbose, dbcs(j)%dirichlet_bc)
         END DO
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE aa_cylindrical_dbc_setup

! **************************************************************************************************
!> \brief  Creteas an evenly-spaced 1D array between two values
!> \param a starting value
!> \param b end value
!> \param N number of evenly-spaced points between a and b
!> \retval arr array containg N evenly-spaced points between a and b
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   FUNCTION linspace(a, b, N) RESULT(arr)

      REAL(dp), INTENT(IN)                               :: a, b
      INTEGER, INTENT(IN)                                :: N
      REAL(dp), DIMENSION(N)                             :: arr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'linspace', routineP = moduleN//':'//routineN

      INTEGER                                            :: i
      REAL(dp)                                           :: dx

      dx = (b-a)/(N-1)
      arr = (/(a+(i-1)*dx, i=1, N)/)

   END FUNCTION linspace

! **************************************************************************************************
!> \brief rotates and translates a cuboid
!> \param cuboid_vtx vertices of the cuboid
!> \param Rmat rotation matrix
!> \param Tpnt translation vector (translates the origin to the center of the cuboid)
!> \param cuboid_transfd_vtx vertices of the transformed cuboid
!> \par History
!>       11.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE rotate_translate_cuboid(cuboid_vtx, Rmat, Tpnt, cuboid_transfd_vtx)

      REAL(dp), DIMENSION(3, 8), INTENT(IN)              :: cuboid_vtx
      REAL(dp), DIMENSION(3, 3), INTENT(OUT)             :: Rmat
      REAL(dp), DIMENSION(3), INTENT(OUT)                :: Tpnt
      REAL(dp), DIMENSION(3, 8), INTENT(OUT)             :: cuboid_transfd_vtx

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rotate_translate_cuboid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(dp), DIMENSION(3)                             :: A, A2, AA2, AB, AD, B, B2, C, C2, D, D2, &
                                                            e1_locl, e2_locl, e3_locl

      CALL timeset(routineN, handle)

      Rmat = 0.0_dp

      A = cuboid_vtx(1:3, 1); A2 = cuboid_vtx(1:3, 6)
      B = cuboid_vtx(1:3, 2); B2 = cuboid_vtx(1:3, 7)
      C = cuboid_vtx(1:3, 3); C2 = cuboid_vtx(1:3, 8)
      D = cuboid_vtx(1:3, 4); D2 = cuboid_vtx(1:3, 5)

      Tpnt = (A+C2)/2.0_dp

      AB = B-A
      AD = D-A
      AA2 = A2-A

! unit vectors generating the local coordinate system
      e1_locl = AB/SQRT(SUM(AB**2))
      e2_locl = AD/SQRT(SUM(AD**2))
      e3_locl = AA2/SQRT(SUM(AA2**2))
! rotation matrix
      Rmat(1:3, 1) = e1_locl
      Rmat(1:3, 2) = e2_locl
      Rmat(1:3, 3) = e3_locl

      cuboid_transfd_vtx(1:3, 1) = rotate_translate_vector(Rmat, Tpnt, BWROT, A)
      cuboid_transfd_vtx(1:3, 2) = rotate_translate_vector(Rmat, Tpnt, BWROT, B)
      cuboid_transfd_vtx(1:3, 3) = rotate_translate_vector(Rmat, Tpnt, BWROT, C)
      cuboid_transfd_vtx(1:3, 4) = rotate_translate_vector(Rmat, Tpnt, BWROT, D)
      cuboid_transfd_vtx(1:3, 5) = rotate_translate_vector(Rmat, Tpnt, BWROT, D2)
      cuboid_transfd_vtx(1:3, 6) = rotate_translate_vector(Rmat, Tpnt, BWROT, A2)
      cuboid_transfd_vtx(1:3, 7) = rotate_translate_vector(Rmat, Tpnt, BWROT, B2)
      cuboid_transfd_vtx(1:3, 8) = rotate_translate_vector(Rmat, Tpnt, BWROT, C2)

      CALL timestop(handle)

   END SUBROUTINE rotate_translate_cuboid

! **************************************************************************************************
!> \brief transforms a position vector according to a given rotation matrix and a
!>        translation vector
!> \param Rmat rotation matrix
!> \param Tp image of the origin under the translation
!> \param direction forward or backward transformation
!> \param vec input vector
!> \retval vec_transfd transformed vector
!> \par History
!>       11.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   FUNCTION rotate_translate_vector(Rmat, Tp, direction, vec) RESULT(vec_transfd)
      REAL(dp), DIMENSION(3, 3), INTENT(IN)              :: Rmat
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: Tp
      INTEGER, INTENT(IN)                                :: direction
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: vec
      REAL(dp), DIMENSION(3)                             :: vec_transfd

      REAL(dp), DIMENSION(3)                             :: Tpoint, vec_tmp
      REAL(dp), DIMENSION(3, 3)                          :: Rmat_inv

      Tpoint = 0.0_dp

      IF (direction .EQ. FWROT) THEN
         CALL matvec_3x3(vec_tmp, Rmat, vec)
         vec_transfd = vec_tmp+Tp
      ELSEIF (direction .EQ. BWROT) THEN
         Rmat_inv = inv_3x3(Rmat)
         Tpoint(1) = Rmat_inv(1, 1)*Tp(1)+Rmat_inv(1, 2)*Tp(2)+Rmat_inv(1, 3)*Tp(3)
         Tpoint(2) = Rmat_inv(2, 1)*Tp(1)+Rmat_inv(2, 2)*Tp(2)+Rmat_inv(2, 3)*Tp(3)
         Tpoint(3) = Rmat_inv(3, 1)*Tp(1)+Rmat_inv(3, 2)*Tp(2)+Rmat_inv(3, 3)*Tp(3)
         CALL matvec_3x3(vec_tmp, Rmat_inv, vec)
         vec_transfd = vec_tmp-Tpoint
      END IF

   END FUNCTION rotate_translate_vector

! **************************************************************************************************
!> \brief  Partitions an axis-aligned cuboid into tiles.
!> \param dbc_vertices vertices of the axis-aligned Dirichlet region to be partitioned
!> \param n_prtn vetor of size 3 specifying the number of times that the x, y and
!>               z interval (defining the region) should be partitioned into
!> \param tiles the obtained tiles
!> \par History
!>       10.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE aa_dbc_partition(dbc_vertices, n_prtn, tiles)

      REAL(dp), DIMENSION(3, 8), INTENT(IN)              :: dbc_vertices
      INTEGER, DIMENSION(3), INTENT(IN)                  :: n_prtn
      TYPE(tile_p_type), DIMENSION(:), INTENT(OUT), &
         POINTER                                         :: tiles

      CHARACTER(LEN=*), PARAMETER :: routineN = 'aa_dbc_partition', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ii, jj, k, kk
      REAL(dp)                                           :: step
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: xprtn_lb, xprtn_ub, yprtn_lb, yprtn_ub, &
                                                            zprtn_lb, zprtn_ub
      REAL(dp), DIMENSION(2)                             :: x_xtnt, y_xtnt, z_xtnt

      CALL timeset(routineN, handle)

! vertices are labeled as:
!               6-------5       A2------D2
!              /|      /|       /|      /|
!             7-|-----8 |     B2-|----C2 |
!             | 1-----|-4  or  | A-----|-D
!             |/      |/       |/      |/
!             2-------3        B-------C

      ALLOCATE (xprtn_lb(n_prtn(1)), xprtn_ub(n_prtn(1)))
      ALLOCATE (yprtn_lb(n_prtn(2)), yprtn_ub(n_prtn(2)))
      ALLOCATE (zprtn_lb(n_prtn(3)), zprtn_ub(n_prtn(3)))

! find the extents of the cuboid in all three directions
      x_xtnt(1) = MINVAL(dbc_vertices(1, :)); x_xtnt(2) = MAXVAL(dbc_vertices(1, :))
      y_xtnt(1) = MINVAL(dbc_vertices(2, :)); y_xtnt(2) = MAXVAL(dbc_vertices(2, :))
      z_xtnt(1) = MINVAL(dbc_vertices(3, :)); z_xtnt(2) = MAXVAL(dbc_vertices(3, :))

! devide the x, y and z extents into n_prtn partitions
      step = (x_xtnt(2)-x_xtnt(1))/REAL(n_prtn(1), kind=dp)
      xprtn_lb(:) = x_xtnt(1)+(/(i, i=0, n_prtn(1)-1)/)*step ! lower bounds
      xprtn_ub(:) = x_xtnt(1)+(/(i, i=1, n_prtn(1))/)*step ! upper bounds

      step = (y_xtnt(2)-y_xtnt(1))/REAL(n_prtn(2), kind=dp)
      yprtn_lb(:) = y_xtnt(1)+(/(i, i=0, n_prtn(2)-1)/)*step
      yprtn_ub(:) = y_xtnt(1)+(/(i, i=1, n_prtn(2))/)*step

      step = (z_xtnt(2)-z_xtnt(1))/REAL(n_prtn(3), kind=dp)
      zprtn_lb(:) = z_xtnt(1)+(/(i, i=0, n_prtn(3)-1)/)*step
      zprtn_ub(:) = z_xtnt(1)+(/(i, i=1, n_prtn(3))/)*step

      ALLOCATE (tiles(n_prtn(1)*n_prtn(2)*n_prtn(3)))
      k = 1
      DO kk = 1, n_prtn(3)
         DO jj = 1, n_prtn(2)
            DO ii = 1, n_prtn(1)
               ALLOCATE (tiles(k)%tile)
               tiles(k)%tile%tile_id = k

               tiles(k)%tile%vertices(1:3, 1) = (/xprtn_lb(ii), yprtn_lb(jj), zprtn_lb(kk)/)
               tiles(k)%tile%vertices(1:3, 2) = (/xprtn_ub(ii), yprtn_lb(jj), zprtn_lb(kk)/)
               tiles(k)%tile%vertices(1:3, 3) = (/xprtn_ub(ii), yprtn_ub(jj), zprtn_lb(kk)/)
               tiles(k)%tile%vertices(1:3, 4) = (/xprtn_lb(ii), yprtn_ub(jj), zprtn_lb(kk)/)
               tiles(k)%tile%vertices(1:3, 5) = (/xprtn_lb(ii), yprtn_ub(jj), zprtn_ub(kk)/)
               tiles(k)%tile%vertices(1:3, 6) = (/xprtn_lb(ii), yprtn_lb(jj), zprtn_ub(kk)/)
               tiles(k)%tile%vertices(1:3, 7) = (/xprtn_ub(ii), yprtn_lb(jj), zprtn_ub(kk)/)
               tiles(k)%tile%vertices(1:3, 8) = (/xprtn_ub(ii), yprtn_ub(jj), zprtn_ub(kk)/)

               tiles(k)%tile%volume = 0.0_dp
               k = k+1
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE aa_dbc_partition

! **************************************************************************************************
!> \brief  Partitions an arbitrary cuboid into tiles.
!> \param dbc_vertices vertices of the axis-aligned Dirichlet region to be partitioned
!> \param n_prtn number of partitions along the edges
!> \param tiles the obtained tiles
!> \par History
!>       10.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE arbitrary_dbc_partition(dbc_vertices, n_prtn, tiles)

      REAL(dp), DIMENSION(3, 8), INTENT(IN)              :: dbc_vertices
      INTEGER, DIMENSION(2), INTENT(IN)                  :: n_prtn
      TYPE(tile_p_type), DIMENSION(:), INTENT(OUT), &
         POINTER                                         :: tiles

      CHARACTER(LEN=*), PARAMETER :: routineN = 'arbitrary_dbc_partition', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, k, l, np1, np2
      REAL(dp)                                           :: ABlength, ADlength, step, thickness
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: end_points, X
      REAL(dp), DIMENSION(3)                             :: A, AB, AC, AD, B, C, D, D2, &
                                                            normal_vector, point1, point2, &
                                                            unit_normal

      CALL timeset(routineN, handle)

      A = dbc_vertices(:, 1)
      B = dbc_vertices(:, 2)
      C = dbc_vertices(:, 3)
      D = dbc_vertices(:, 4)
      D2 = dbc_vertices(:, 5)

      AB = B-A
      AC = C-A
      AD = D-A
      normal_vector = vector_product(AB, AC)
      unit_normal = normal_vector/SQRT(SUM(normal_vector**2))
      thickness = SQRT(SUM((D-D2)**2))

! the larger n_prtn is assigned to the longer edge
      ABlength = SQRT(SUM(AB**2))
      ADlength = SQRT(SUM(AD**2))
      IF (ADlength .GE. ABlength) THEN
         np1 = MAX(n_prtn(1)+1, n_prtn(2)+1)
         np2 = MIN(n_prtn(1)+1, n_prtn(2)+1)
      ELSE
         np1 = MIN(n_prtn(1)+1, n_prtn(2)+1)
         np2 = MAX(n_prtn(1)+1, n_prtn(2)+1)
      END IF

      ALLOCATE (X(np1, np2, 3))
      ALLOCATE (end_points(2, np1, 3))

! partition AD and BC
      DO l = 1, np1
         step = REAL((l-1), kind=dp)/REAL((np1-1), kind=dp)
         end_points(1, l, :) = A*(1.0_dp-step)+D*step
         end_points(2, l, :) = B*(1.0_dp-step)+C*step
      END DO

! partition in the second direction along the line segments with endpoints from
! the previous partitioning step
      DO l = 1, np1
         point1(:) = end_points(1, l, :)
         point2(:) = end_points(2, l, :)
         DO i = 1, np2
            step = REAL((i-1), kind=dp)/REAL((np2-1), kind=dp)
            X(l, i, :) = point1*(1.0_dp-step)+point2*step
         END DO
      END DO

      ALLOCATE (tiles((np1-1)*(np2-1)))
      k = 1
      DO l = 1, np1-1
         DO i = 1, np2-1
            ALLOCATE (tiles(k)%tile)
            tiles(k)%tile%tile_id = k

            tiles(k)%tile%vertices(1:3, 1) = X(l, i, :)
            tiles(k)%tile%vertices(1:3, 2) = X(l, i+1, :)
            tiles(k)%tile%vertices(1:3, 3) = X(l+1, i+1, :)
            tiles(k)%tile%vertices(1:3, 4) = X(l+1, i, :)
            tiles(k)%tile%vertices(1:3, 5) = X(l+1, i, :)+thickness*unit_normal
            tiles(k)%tile%vertices(1:3, 6) = X(l, i, :)+thickness*unit_normal
            tiles(k)%tile%vertices(1:3, 7) = X(l, i+1, :)+thickness*unit_normal
            tiles(k)%tile%vertices(1:3, 8) = X(l+1, i+1, :)+thickness*unit_normal

            tiles(k)%tile%volume = 0.0_dp
            k = k+1
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE arbitrary_dbc_partition

! **************************************************************************************************
!> \brief computes the pw data corresponding to an axis-aligned tile
!> \param cell_xtnts extents of the simulation cell
!> \param x_locl x grid vetor of the simulation box local to this process
!> \param y_locl y grid vetor of the simulation box local to this process
!> \param z_locl z grid vetor of the simulation box local to this process
!> \param tile_vertices coordinates of the vertices of the input tile
!> \param sigma the width of the transition region between 0 and 1 at the edges of the tile
!> \param is_periodic whether or not the tile is periodic with a periodicity of one unit cell
!> \param tile_pw the output pw data
!> \par History
!>       10.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE aa_tile_pw_compute(cell_xtnts, x_locl, y_locl, z_locl, tile_vertices, &
                                 sigma, is_periodic, tile_pw)

      REAL(dp), DIMENSION(2, 3), INTENT(IN)              :: cell_xtnts
      REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(IN)    :: x_locl, y_locl, z_locl
      REAL(dp), DIMENSION(3, 8), INTENT(IN)              :: tile_vertices
      REAL(dp), INTENT(IN)                               :: sigma
      LOGICAL, INTENT(IN)                                :: is_periodic
      TYPE(pw_type), INTENT(INOUT), POINTER              :: tile_pw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'aa_tile_pw_compute', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, k, lb1, lb2, lb3, ub1, &
                                                            ub2, ub3
      INTEGER, DIMENSION(2, 3)                           :: bounds_local
      REAL(dp)                                           :: fx, fy, fz, gx, gy, gz, hx, hy, hz, xi0, &
                                                            xil, xir, yj0, yjl, yjr, zk0, zkl, zkr
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: xlocll, xloclr, ylocll, yloclr, zlocll, &
                                                            zloclr
      REAL(dp), DIMENSION(2)                             :: x_xtnt, y_xtnt, z_xtnt
      REAL(dp), DIMENSION(3)                             :: per
      REAL(dp), DIMENSION(3, 3)                          :: prefactor

      CALL timeset(routineN, handle)

      bounds_local = tile_pw%pw_grid%bounds_local
      lb1 = bounds_local(1, 1); ub1 = bounds_local(2, 1)
      lb2 = bounds_local(1, 2); ub2 = bounds_local(2, 2)
      lb3 = bounds_local(1, 3); ub3 = bounds_local(2, 3)

      ALLOCATE (xlocll(lb1:ub1), ylocll(lb2:ub2), zlocll(lb3:ub3))
      ALLOCATE (xloclr(lb1:ub1), yloclr(lb2:ub2), zloclr(lb3:ub3))

! periodicities
      per = cell_xtnts(2, :)-cell_xtnts(1, :)

      x_xtnt(1) = MINVAL(tile_vertices(1, :)); x_xtnt(2) = MAXVAL(tile_vertices(1, :))
      y_xtnt(1) = MINVAL(tile_vertices(2, :)); y_xtnt(2) = MAXVAL(tile_vertices(2, :))
      z_xtnt(1) = MINVAL(tile_vertices(3, :)); z_xtnt(2) = MAXVAL(tile_vertices(3, :))

! prefactor: columns: left, original, right - rows: x , y, z directions
      prefactor = 0.5_dp

      IF (is_periodic) THEN
         DO k = lb3, ub3
            DO j = lb2, ub2
               DO i = lb1, ub1
                  xi0 = x_locl(i); yj0 = y_locl(j); zk0 = z_locl(k)
                  xil = x_locl(i)-per(1); yjl = y_locl(j)-per(2); zkl = z_locl(k)-per(3)
                  xir = x_locl(i)+per(1); yjr = y_locl(j)+per(2); zkr = z_locl(k)+per(3)

                  ! points from the original cell local to the current processor
                  fx = prefactor(1, 2)*(erf((xi0-x_xtnt(1))/sigma)-erf((xi0-x_xtnt(2))/sigma))
                  fy = prefactor(2, 2)*(erf((yj0-y_xtnt(1))/sigma)-erf((yj0-y_xtnt(2))/sigma))
                  fz = prefactor(3, 2)*(erf((zk0-z_xtnt(1))/sigma)-erf((zk0-z_xtnt(2))/sigma))

                  ! periodically replicated cell on the left, points local to the current processor
                  gx = prefactor(1, 1)*(erf((xil-x_xtnt(1))/sigma)-erf((xil-x_xtnt(2))/sigma))
                  gy = prefactor(2, 1)*(erf((yjl-y_xtnt(1))/sigma)-erf((yjl-y_xtnt(2))/sigma))
                  gz = prefactor(3, 1)*(erf((zkl-z_xtnt(1))/sigma)-erf((zkl-z_xtnt(2))/sigma))

                  ! periodically replicated cell on the right, points local to the current processor
                  hx = prefactor(1, 3)*(erf((xir-x_xtnt(1))/sigma)-erf((xir-x_xtnt(2))/sigma))
                  hy = prefactor(2, 3)*(erf((yjr-y_xtnt(1))/sigma)-erf((yjr-y_xtnt(2))/sigma))
                  hz = prefactor(3, 3)*(erf((zkr-z_xtnt(1))/sigma)-erf((zkr-z_xtnt(2))/sigma))

                  tile_pw%cr3d(i, j, k) = (fx+gx+hx)*(fy+gy+hy)*(fz+gz+hz)
               END DO
            END DO
         END DO
      ELSE
         DO k = lb3, ub3
            DO j = lb2, ub2
               DO i = lb1, ub1
                  xi0 = x_locl(i); yj0 = y_locl(j); zk0 = z_locl(k)

                  fx = 0.5_dp*(erf((xi0-x_xtnt(1))/sigma)-erf((xi0-x_xtnt(2))/sigma))
                  fy = 0.5_dp*(erf((yj0-y_xtnt(1))/sigma)-erf((yj0-y_xtnt(2))/sigma))
                  fz = 0.5_dp*(erf((zk0-z_xtnt(1))/sigma)-erf((zk0-z_xtnt(2))/sigma))

                  tile_pw%cr3d(i, j, k) = fx*fy*fz
               END DO
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE aa_tile_pw_compute

! **************************************************************************************************
!> \brief computes the pw data corresponding to an arbitrary (rectangular-shaped) tile
!> \param cell_xtnts extents of the simulation cell
!> \param x_locl x grid vetor of the simulation box local to this process
!> \param y_locl y grid vetor of the simulation box local to this process
!> \param z_locl z grid vetor of the simulation box local to this process
!> \param tile_vertices coordinates of the vertices of the input tile
!> \param sigma the width of the transition region between 0 and 1 at the edges of the tile
!> \param tile_pw the output pw data
!> \par History
!>       11.2015 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE arbitrary_tile_pw_compute(cell_xtnts, x_locl, y_locl, z_locl, tile_vertices, &
                                        sigma, tile_pw)

      REAL(dp), DIMENSION(2, 3), INTENT(IN)              :: cell_xtnts
      REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(IN)    :: x_locl, y_locl, z_locl
      REAL(dp), DIMENSION(3, 8), INTENT(IN)              :: tile_vertices
      REAL(dp), INTENT(IN)                               :: sigma
      TYPE(pw_type), INTENT(INOUT), POINTER              :: tile_pw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'arbitrary_tile_pw_compute', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, k, lb1, lb2, lb3, ub1, &
                                                            ub2, ub3
      INTEGER, DIMENSION(2, 3)                           :: bounds_local
      REAL(dp)                                           :: cm_x, cm_y, cm_z, fx, fy, fz, xi0, yj0, &
                                                            zk0
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: xlocll, xloclr, ylocll, yloclr, zlocll, &
                                                            zloclr
      REAL(dp), DIMENSION(2)                             :: transfdx_xtnt, transfdy_xtnt, &
                                                            transfdz_xtnt, x_xtnt, y_xtnt, z_xtnt
      REAL(dp), DIMENSION(3)                             :: per, pnt0, Tpnt
      REAL(dp), DIMENSION(3, 3)                          :: prefactor, Rmat
      REAL(dp), DIMENSION(3, 8)                          :: tile_transfd_vertices

      CALL timeset(routineN, handle)

      bounds_local = tile_pw%pw_grid%bounds_local
      lb1 = bounds_local(1, 1); ub1 = bounds_local(2, 1)
      lb2 = bounds_local(1, 2); ub2 = bounds_local(2, 2)
      lb3 = bounds_local(1, 3); ub3 = bounds_local(2, 3)

      ALLOCATE (xlocll(lb1:ub1), ylocll(lb2:ub2), zlocll(lb3:ub3))
      ALLOCATE (xloclr(lb1:ub1), yloclr(lb2:ub2), zloclr(lb3:ub3))

      per = cell_xtnts(2, :)-cell_xtnts(1, :)

      CALL rotate_translate_cuboid(tile_vertices, Rmat, Tpnt, tile_transfd_vertices)

      transfdx_xtnt(1) = MINVAL(tile_transfd_vertices(1, :))
      transfdx_xtnt(2) = MAXVAL(tile_transfd_vertices(1, :))
      transfdy_xtnt(1) = MINVAL(tile_transfd_vertices(2, :))
      transfdy_xtnt(2) = MAXVAL(tile_transfd_vertices(2, :))
      transfdz_xtnt(1) = MINVAL(tile_transfd_vertices(3, :))
      transfdz_xtnt(2) = MAXVAL(tile_transfd_vertices(3, :))

      cm_x = 0.5_dp*(transfdx_xtnt(2)+transfdx_xtnt(1))
      cm_y = 0.5_dp*(transfdy_xtnt(2)+transfdy_xtnt(1))
      cm_z = 0.5_dp*(transfdz_xtnt(2)+transfdz_xtnt(1))

      x_xtnt = transfdx_xtnt-cm_x
      y_xtnt = transfdy_xtnt-cm_y
      z_xtnt = transfdz_xtnt-cm_z

      prefactor = 0.5_dp

      DO k = lb3, ub3
         DO j = lb2, ub2
            DO i = lb1, ub1
               pnt0 = rotate_translate_vector(Rmat, Tpnt, BWROT, (/x_locl(i), y_locl(j), z_locl(k)/))

               xi0 = pnt0(1); yj0 = pnt0(2); zk0 = pnt0(3)

               fx = 0.5_dp*(erf((xi0-x_xtnt(1))/sigma)-erf((xi0-x_xtnt(2))/sigma))
               fy = 0.5_dp*(erf((yj0-y_xtnt(1))/sigma)-erf((yj0-y_xtnt(2))/sigma))
               fz = 0.5_dp*(erf((zk0-z_xtnt(1))/sigma)-erf((zk0-z_xtnt(2))/sigma))

               tile_pw%cr3d(i, j, k) = fx*fy*fz
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE arbitrary_tile_pw_compute

END MODULE dirichlet_bc_methods
