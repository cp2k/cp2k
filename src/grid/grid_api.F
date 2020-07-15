!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
MODULE grid_api
   USE ISO_C_BINDING,                   ONLY: C_BOOL,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_PTR
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cube_utils,                      ONLY: cube_info_type
   USE grid_collocate,                  ONLY: &
        GRID_FUNC_AB, GRID_FUNC_ADBmDAB_X, GRID_FUNC_ADBmDAB_Y, GRID_FUNC_ADBmDAB_Z, &
        GRID_FUNC_ARDBmDARB_XX, GRID_FUNC_ARDBmDARB_XY, GRID_FUNC_ARDBmDARB_XZ, &
        GRID_FUNC_ARDBmDARB_YX, GRID_FUNC_ARDBmDARB_YY, GRID_FUNC_ARDBmDARB_YZ, &
        GRID_FUNC_ARDBmDARB_ZX, GRID_FUNC_ARDBmDARB_ZY, GRID_FUNC_ARDBmDARB_ZZ, &
        GRID_FUNC_DABpADB_X, GRID_FUNC_DABpADB_Y, GRID_FUNC_DABpADB_Z, GRID_FUNC_DADB, &
        GRID_FUNC_DX, GRID_FUNC_DXDX, GRID_FUNC_DXDY, GRID_FUNC_DY, GRID_FUNC_DYDY, &
        GRID_FUNC_DYDZ, GRID_FUNC_DZ, GRID_FUNC_DZDX, GRID_FUNC_DZDZ, &
        collocate_pgf_product_legacy => collocate_pgf_product
   USE grid_integrate,                  ONLY: integrate_pgf_product
   USE kinds,                           ONLY: dp
   USE realspace_grid_types,            ONLY: realspace_grid_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grid_api'

   LOGICAL, PARAMETER, PRIVATE :: validate_collocate = .FALSE.

   PUBLIC :: GRID_FUNC_AB
   PUBLIC :: GRID_FUNC_DADB
   PUBLIC :: GRID_FUNC_ADBmDAB_X
   PUBLIC :: GRID_FUNC_ADBmDAB_Y
   PUBLIC :: GRID_FUNC_ADBmDAB_Z
   PUBLIC :: GRID_FUNC_ARDBmDARB_XX
   PUBLIC :: GRID_FUNC_ARDBmDARB_XY
   PUBLIC :: GRID_FUNC_ARDBmDARB_XZ
   PUBLIC :: GRID_FUNC_ARDBmDARB_YX
   PUBLIC :: GRID_FUNC_ARDBmDARB_YY
   PUBLIC :: GRID_FUNC_ARDBmDARB_YZ
   PUBLIC :: GRID_FUNC_ARDBmDARB_ZX
   PUBLIC :: GRID_FUNC_ARDBmDARB_ZY
   PUBLIC :: GRID_FUNC_ARDBmDARB_ZZ
   PUBLIC :: GRID_FUNC_DABpADB_X
   PUBLIC :: GRID_FUNC_DABpADB_Y
   PUBLIC :: GRID_FUNC_DABpADB_Z
   PUBLIC :: GRID_FUNC_DX
   PUBLIC :: GRID_FUNC_DY
   PUBLIC :: GRID_FUNC_DZ
   PUBLIC :: GRID_FUNC_DXDY
   PUBLIC :: GRID_FUNC_DYDZ
   PUBLIC :: GRID_FUNC_DZDX
   PUBLIC :: GRID_FUNC_DXDX
   PUBLIC :: GRID_FUNC_DYDY
   PUBLIC :: GRID_FUNC_DZDZ

   PUBLIC :: collocate_pgf_product
   PUBLIC :: integrate_pgf_product

CONTAINS

! **************************************************************************************************
!> \brief low level collocation of primitive gaussian functions
!> \param la_max ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param zetb ...
!> \param lb_min ...
!> \param ra ...
!> \param rab ...
!> \param scale ...
!> \param pab ...
!> \param o1 ...
!> \param o2 ...
!> \param rsgrid ...
!> \param cell ...
!> \param cube_info ...
!> \param ga_gb_function ...
!> \param radius ...
!> \param use_subpatch ...
!> \param subpatch_pattern ...
! **************************************************************************************************
   SUBROUTINE collocate_pgf_product(la_max, zeta, la_min, &
                                    lb_max, zetb, lb_min, &
                                    ra, rab, scale, pab, o1, o2, &
                                    rsgrid, cell, cube_info, &
                                    ga_gb_function, radius, &
                                    use_subpatch, subpatch_pattern)

      INTEGER, INTENT(IN)                                :: la_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zetb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), TARGET    :: ra, rab
      REAL(KIND=dp), INTENT(IN)                          :: scale
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      INTEGER, INTENT(IN)                                :: o1, o2
      TYPE(realspace_grid_type)                          :: rsgrid
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type), INTENT(IN)                   :: cube_info
      INTEGER, INTENT(IN)                                :: ga_gb_function
      REAL(KIND=dp), INTENT(IN)                          :: radius
      LOGICAL, OPTIONAL                                  :: use_subpatch
      INTEGER, INTENT(IN), OPTIONAL                      :: subpatch_pattern

      CHARACTER(len=*), PARAMETER :: routineN = 'collocate_pgf_product', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iw, j, k, my_pos, my_subpatch_pattern
      INTEGER, DIMENSION(3), TARGET                      :: npts_global, npts_local, shift_local
      INTEGER, SAVE                                      :: counter = 0
      LOGICAL(KIND=C_BOOL)                               :: my_use_subpatch, orthorhombic
      REAL(KIND=dp)                                      :: diff, ref_value, rel_diff, tolerance
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid, grid_before
      INTERFACE
         SUBROUTINE grid_collocate_pgf_product_cpu(orthorhombic, &
                                                   use_subpatch, subpatch, border, func, &
                                                   la_max, la_min, lb_max, lb_min, &
                                                   zeta, zetb, rscale, dh, dh_inv, ra, rab, &
                                                   npts_global, npts_local, shift_local, &
                                                   radius, o1, o2, n1, n2, pab, &
                                                   grid) &
            BIND(C, name="grid_collocate_pgf_product_cpu")
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_BOOL
            LOGICAL(KIND=C_BOOL), VALUE               :: orthorhombic
            LOGICAL(KIND=C_BOOL), VALUE               :: use_subpatch
            INTEGER(KIND=C_INT), VALUE                :: subpatch
            INTEGER(KIND=C_INT), VALUE                :: border
            INTEGER(KIND=C_INT), VALUE                :: func
            INTEGER(KIND=C_INT), VALUE                :: la_max
            INTEGER(KIND=C_INT), VALUE                :: la_min
            INTEGER(KIND=C_INT), VALUE                :: lb_max
            INTEGER(KIND=C_INT), VALUE                :: lb_min
            REAL(KIND=C_DOUBLE), VALUE                :: zeta
            REAL(KIND=C_DOUBLE), VALUE                :: zetb
            REAL(KIND=C_DOUBLE), VALUE                :: rscale
            TYPE(C_PTR), VALUE                        :: dh
            TYPE(C_PTR), VALUE                        :: dh_inv
            TYPE(C_PTR), VALUE                        :: ra
            TYPE(C_PTR), VALUE                        :: rab
            TYPE(C_PTR), VALUE                        :: npts_global
            TYPE(C_PTR), VALUE                        :: npts_local
            TYPE(C_PTR), VALUE                        :: shift_local
            REAL(KIND=C_DOUBLE), VALUE                :: radius
            INTEGER(KIND=C_INT), VALUE                :: o1
            INTEGER(KIND=C_INT), VALUE                :: o2
            INTEGER(KIND=C_INT), VALUE                :: n1
            INTEGER(KIND=C_INT), VALUE                :: n2
            TYPE(C_PTR), VALUE                        :: pab
            TYPE(C_PTR), VALUE                        :: grid
         END SUBROUTINE grid_collocate_pgf_product_cpu
      END INTERFACE

      !TODO remove unused dummy arguments
      MARK_USED(cell)
      MARK_USED(cube_info)

      IF (PRESENT(use_subpatch)) THEN
         my_use_subpatch = use_subpatch
      ELSE
         my_use_subpatch = .FALSE.
      ENDIF
      CPASSERT(.NOT. my_use_subpatch .OR. PRESENT(subpatch_pattern))

      my_subpatch_pattern = 0
      IF (PRESENT(subpatch_pattern)) THEN
         my_subpatch_pattern = subpatch_pattern
      ENDIF

      IF (2.0_dp*radius < MAXVAL(ABS(rsgrid%desc%dh))) THEN
         RETURN  ! radius too small to be mapped onto grid of given resolution
      END IF

      orthorhombic = rsgrid%desc%orthorhombic

      CPASSERT(LBOUND(pab, 1) == 1)
      CPASSERT(LBOUND(pab, 2) == 1)

      CALL get_rsgrid_properties(rsgrid, npts_global=npts_global, npts_local=npts_local, shift_local=shift_local)

      IF (validate_collocate) THEN
         ALLOCATE (grid(SIZE(rsgrid%r, 1), SIZE(rsgrid%r, 2), SIZE(rsgrid%r, 3)))
         ALLOCATE (grid_before(SIZE(rsgrid%r, 1), SIZE(rsgrid%r, 2), SIZE(rsgrid%r, 3)))
         grid_before = rsgrid%r ! copy
         grid = 0.0_dp
         rsgrid%r = 0.0_dp
      ELSE
         grid(1:, 1:, 1:) => rsgrid%r(:, :, :)  ! pointer assignment
      ENDIF

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(rsgrid%desc%dh))
      CPASSERT(IS_CONTIGUOUS(rsgrid%desc%dh_inv))
      CPASSERT(IS_CONTIGUOUS(ra))
      CPASSERT(IS_CONTIGUOUS(rab))
      CPASSERT(IS_CONTIGUOUS(npts_global))
      CPASSERT(IS_CONTIGUOUS(npts_local))
      CPASSERT(IS_CONTIGUOUS(shift_local))
      CPASSERT(IS_CONTIGUOUS(pab))
      CPASSERT(IS_CONTIGUOUS(grid))
#endif

      CALL grid_collocate_pgf_product_cpu(orthorhombic=orthorhombic, &
                                          use_subpatch=my_use_subpatch, &
                                          subpatch=my_subpatch_pattern, &
                                          border=rsgrid%desc%border, &
                                          func=ga_gb_function, &
                                          la_max=la_max, &
                                          la_min=la_min, &
                                          lb_max=lb_max, &
                                          lb_min=lb_min, &
                                          zeta=zeta, &
                                          zetb=zetb, &
                                          rscale=scale, &
                                          dh=C_LOC(rsgrid%desc%dh(1, 1)), &
                                          dh_inv=C_LOC(rsgrid%desc%dh_inv(1, 1)), &
                                          ra=C_LOC(ra(1)), &
                                          rab=C_LOC(rab(1)), &
                                          npts_global=C_LOC(npts_global(1)), &
                                          npts_local=C_LOC(npts_local(1)), &
                                          shift_local=C_LOC(shift_local(1)), &
                                          radius=radius, &
                                          o1=o1, &
                                          o2=o2, &
                                          n1=SIZE(pab, 1), &
                                          n2=SIZE(pab, 2), &
                                          pab=C_LOC(pab(1, 1)), &
                                          grid=C_LOC(grid(1, 1, 1)))

      IF (validate_collocate) THEN
         CALL collocate_pgf_product_legacy(la_max, zeta, la_min, &
                                           lb_max, zetb, lb_min, &
                                           ra, rab, scale, pab, o1, o2, &
                                           rsgrid, cell, cube_info, &
                                           ga_gb_function, radius, &
                                           use_subpatch, subpatch_pattern)
         IF (orthorhombic .AND. .NOT. my_use_subpatch) THEN
            tolerance = 6e-14
         ELSE
            tolerance = 2e-10
         END IF
         counter = counter + 1
         my_pos = rsgrid%desc%my_pos
         DO i = 1, SIZE(rsgrid%r, 1)
         DO j = 1, SIZE(rsgrid%r, 2)
         DO k = 1, SIZE(rsgrid%r, 3)
            ref_value = rsgrid%r(i - 1 + LBOUND(rsgrid%r, 1), j - 1 + LBOUND(rsgrid%r, 2), k - 1 + LBOUND(rsgrid%r, 3))
            diff = ABS(grid(i, j, k) - ref_value)
            rel_diff = diff/MAX(1.0_dp, ABS(ref_value))
            !IF ((ref_value==0.0_dp) .NEQV. (grid(i, j, k)==0.0_dp)) THEN
            IF (rel_diff > tolerance) THEN
               iw = cp_logger_get_default_io_unit()
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "diff: ", diff
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "rel. diff: ", rel_diff
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "value:", grid(i, j, k), "ref_value:", ref_value
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "ijk:", i - 1, j - 1, k - 1
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "npts_global:", npts_global
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "npts_local:", npts_local
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "shift_local:", shift_local
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "orthorhombic:", orthorhombic
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "ga_gb_function:", ga_gb_function
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "use_subpatch:", my_use_subpatch
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "subpatch_pattern:", my_subpatch_pattern
               WRITE (iw, *) "counter: ", counter, "my_pos: ", my_pos, "border:", rsgrid%desc%border
               CPABORT("Collocate validation failed.")
            END IF
         ENDDO
         ENDDO
         ENDDO
         rsgrid%r = rsgrid%r + grid_before
         DEALLOCATE (grid, grid_before)
      ENDIF

   END SUBROUTINE collocate_pgf_product

! **************************************************************************************************
!> \brief Helper routines for getting rsgrid properties and asserting underlying assumptions.
!> \param rsgrid ...
!> \param npts_global ...
!> \param npts_local ...
!> \param shift_local ...
! **************************************************************************************************
   SUBROUTINE get_rsgrid_properties(rsgrid, npts_global, npts_local, shift_local)
      TYPE(realspace_grid_type), INTENT(IN)              :: rsgrid
      INTEGER, DIMENSION(:)                              :: npts_global, npts_local, shift_local

      INTEGER                                            :: i

      ! See rs_grid_create() in ./src/pw/realspace_grid_types.F.
      CPASSERT(LBOUND(rsgrid%r, 1) == rsgrid%lb_local(1))
      CPASSERT(UBOUND(rsgrid%r, 1) == rsgrid%ub_local(1))
      CPASSERT(LBOUND(rsgrid%r, 2) == rsgrid%lb_local(2))
      CPASSERT(UBOUND(rsgrid%r, 2) == rsgrid%ub_local(2))
      CPASSERT(LBOUND(rsgrid%r, 3) == rsgrid%lb_local(3))
      CPASSERT(UBOUND(rsgrid%r, 3) == rsgrid%ub_local(3))

      ! While the rsgrid code assumes that the grid starts at rsgrid%lb,
      ! the collocate code assumes that the grid starts at (1,1,1) in Fortran, or (0,0,0) in C.
      ! So, a point rp(:) gets the following grid coordinates MODULO(rp(:)/dr(:),npts_global(:))

      ! Number of global grid points in each direction.
      npts_global = rsgrid%desc%ub - rsgrid%desc%lb + 1

      ! Number of local grid points in each direction.
      npts_local = rsgrid%ub_local - rsgrid%lb_local + 1

      ! Number of points the local grid is shifted wrt global grid.
      shift_local = rsgrid%lb_local - rsgrid%desc%lb

      ! Check that rsgrid%desc%perd means what I think it means.
      DO i = 1, 3
         IF (rsgrid%desc%perd(i) == 1) THEN
            ! Periodic meaning the grid in this direction is entriely present on every processor.
            CPASSERT(npts_local(i) == npts_global(i))
            CPASSERT(shift_local(i) == 0)
         ELSE
            ! Not periodic meaning the grid in this direction is distributed among processors.
            CPASSERT(npts_local(i) < npts_global(i))
            ! Check bounds of grid section that is owned by this processor.
            CPASSERT(rsgrid%lb_real(i) == rsgrid%lb_local(i) + rsgrid%desc%border)
            CPASSERT(rsgrid%ub_real(i) == rsgrid%ub_local(i) - rsgrid%desc%border)
         ENDIF
      ENDDO
   END SUBROUTINE get_rsgrid_properties

END MODULE grid_api
