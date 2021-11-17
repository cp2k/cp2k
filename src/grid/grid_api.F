!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: BSD-3-Clause                                                          !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Fortran API for the grid package, which is written in C.
!> \author Ole Schuett
! **************************************************************************************************
MODULE grid_api
   USE ISO_C_BINDING,                   ONLY: &
        C_ASSOCIATED, C_BOOL, C_CHAR, C_DOUBLE, C_FUNLOC, C_FUNPTR, C_INT, C_LOC, C_LONG, &
        C_NULL_PTR, C_PTR
   USE cell_types,                      ONLY: cell_type
   USE cube_utils,                      ONLY: cube_info_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE offload_api,                     ONLY: offload_buffer_type
   USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                              realspace_grid_type
   USE string_utilities,                ONLY: strlcpy_c2f
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grid_api'

   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_AB = 100
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DADB = 200
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ADBmDAB_X = 301
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ADBmDAB_Y = 302
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ADBmDAB_Z = 303
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_XX = 411
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_XY = 412
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_XZ = 413
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_YX = 421
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_YY = 422
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_YZ = 423
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_ZX = 431
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_ZY = 432
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_ZZ = 433
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DABpADB_X = 501
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DABpADB_Y = 502
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DABpADB_Z = 503
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DX = 601
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DY = 602
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DZ = 603
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DXDY = 701
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DYDZ = 702
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DZDX = 703
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DXDX = 801
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DYDY = 802
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DZDZ = 803
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DAB_X = 901
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DAB_Y = 902
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DAB_Z = 903
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ADB_X = 904
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ADB_Y = 905
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ADB_Z = 906

   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_CORE_X = 1001
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_CORE_Y = 1002
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_CORE_Z = 1003

   INTEGER, PARAMETER, PUBLIC :: GRID_BACKEND_AUTO = 10
   INTEGER, PARAMETER, PUBLIC :: GRID_BACKEND_REF = 11
   INTEGER, PARAMETER, PUBLIC :: GRID_BACKEND_CPU = 12
   INTEGER, PARAMETER, PUBLIC :: GRID_BACKEND_GPU = 13

   PUBLIC :: grid_library_init, grid_library_finalize
   PUBLIC :: grid_library_set_config, grid_library_print_stats
   PUBLIC :: collocate_pgf_product, integrate_pgf_product
   PUBLIC :: grid_basis_set_type, grid_create_basis_set, grid_free_basis_set
   PUBLIC :: grid_task_list_type, grid_create_task_list, grid_free_task_list
   PUBLIC :: grid_collocate_task_list, grid_integrate_task_list

   TYPE grid_basis_set_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   END TYPE grid_basis_set_type

   TYPE grid_task_list_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   END TYPE grid_task_list_type

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
!> \author Ole Schuett
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

      INTEGER                                            :: border_mask
      INTEGER, DIMENSION(3), TARGET                      :: border_width, npts_global, npts_local, &
                                                            shift_local
      LOGICAL(KIND=C_BOOL)                               :: orthorhombic
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid
      INTERFACE
         SUBROUTINE grid_ref_collocate_pgf_product_c(orthorhombic, &
                                                     border_mask, func, &
                                                     la_max, la_min, lb_max, lb_min, &
                                                     zeta, zetb, rscale, dh, dh_inv, ra, rab, &
                                                     npts_global, npts_local, shift_local, border_width, &
                                                     radius, o1, o2, n1, n2, pab, &
                                                     grid) &
            BIND(C, name="grid_ref_collocate_pgf_product")
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_BOOL
            LOGICAL(KIND=C_BOOL), VALUE               :: orthorhombic
            INTEGER(KIND=C_INT), VALUE                :: border_mask
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
            TYPE(C_PTR), VALUE                        :: border_width
            REAL(KIND=C_DOUBLE), VALUE                :: radius
            INTEGER(KIND=C_INT), VALUE                :: o1
            INTEGER(KIND=C_INT), VALUE                :: o2
            INTEGER(KIND=C_INT), VALUE                :: n1
            INTEGER(KIND=C_INT), VALUE                :: n2
            TYPE(C_PTR), VALUE                        :: pab
            TYPE(C_PTR), VALUE                        :: grid
         END SUBROUTINE grid_ref_collocate_pgf_product_c
      END INTERFACE

      !TODO remove unused dummy arguments
      MARK_USED(cell)
      MARK_USED(cube_info)

      border_mask = 0
      IF (PRESENT(use_subpatch)) THEN
         IF (use_subpatch) THEN
            CPASSERT(PRESENT(subpatch_pattern))
            border_mask = IAND(63, NOT(subpatch_pattern))  ! invert last 6 bits
         END IF
      END IF

      orthorhombic = LOGICAL(rsgrid%desc%orthorhombic, C_BOOL)

      CPASSERT(LBOUND(pab, 1) == 1)
      CPASSERT(LBOUND(pab, 2) == 1)

      CALL get_rsgrid_properties(rsgrid, npts_global=npts_global, &
                                 npts_local=npts_local, &
                                 shift_local=shift_local, &
                                 border_width=border_width)

      grid(1:, 1:, 1:) => rsgrid%r(:, :, :)  ! pointer assignment

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(rsgrid%desc%dh))
      CPASSERT(IS_CONTIGUOUS(rsgrid%desc%dh_inv))
      CPASSERT(IS_CONTIGUOUS(ra))
      CPASSERT(IS_CONTIGUOUS(rab))
      CPASSERT(IS_CONTIGUOUS(npts_global))
      CPASSERT(IS_CONTIGUOUS(npts_local))
      CPASSERT(IS_CONTIGUOUS(shift_local))
      CPASSERT(IS_CONTIGUOUS(border_width))
      CPASSERT(IS_CONTIGUOUS(pab))
      CPASSERT(IS_CONTIGUOUS(grid))
#endif

      ! For collocating single pgf products there is only the reference implementation.

      CALL grid_ref_collocate_pgf_product_c(orthorhombic=orthorhombic, &
                                            border_mask=border_mask, &
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
                                            border_width=C_LOC(border_width(1)), &
                                            radius=radius, &
                                            o1=o1, &
                                            o2=o2, &
                                            n1=SIZE(pab, 1), &
                                            n2=SIZE(pab, 2), &
                                            pab=C_LOC(pab(1, 1)), &
                                            grid=C_LOC(grid(1, 1, 1)))

   END SUBROUTINE collocate_pgf_product

! **************************************************************************************************
!> \brief low level function to compute matrix elements of primitive gaussian functions
!> \param la_max ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param zetb ...
!> \param lb_min ...
!> \param ra ...
!> \param rab ...
!> \param rsgrid ...
!> \param cell ...
!> \param cube_info ...
!> \param hab ...
!> \param pab ...
!> \param o1 ...
!> \param o2 ...
!> \param radius ...
!> \param calculate_forces ...
!> \param force_a ...
!> \param force_b ...
!> \param compute_tau ...
!> \param use_virial ...
!> \param my_virial_a ...
!> \param my_virial_b ...
!> \param hdab Derivative with respect to the primitive on the left.
!> \param hadb Derivative with respect to the primitive on the right.
!> \param a_hdab ...
!> \param use_subpatch ...
!> \param subpatch_pattern ...
! **************************************************************************************************
   SUBROUTINE integrate_pgf_product(la_max, zeta, la_min, &
                                    lb_max, zetb, lb_min, &
                                    ra, rab, rsgrid, cell, &
                                    cube_info, hab, pab, o1, o2, &
                                    radius, &
                                    calculate_forces, force_a, force_b, &
                                    compute_tau, &
                                    use_virial, my_virial_a, &
                                    my_virial_b, hdab, hadb, a_hdab, use_subpatch, subpatch_pattern)

      INTEGER, INTENT(IN)                                :: la_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zetb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), TARGET    :: ra, rab
      TYPE(realspace_grid_type), POINTER                 :: rsgrid
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type), INTENT(IN)                   :: cube_info
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hab
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: pab
      INTEGER, INTENT(IN)                                :: o1, o2
      REAL(KIND=dp), INTENT(IN)                          :: radius
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT), &
         OPTIONAL                                        :: force_a, force_b
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_tau, use_virial
      REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL           :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: hdab, hadb
      REAL(KIND=dp), DIMENSION(:, :, :, :), OPTIONAL, &
         POINTER                                         :: a_hdab
      LOGICAL, OPTIONAL                                  :: use_subpatch
      INTEGER, INTENT(IN), OPTIONAL                      :: subpatch_pattern

      INTEGER                                            :: border_mask
      INTEGER, DIMENSION(3), TARGET                      :: border_width, npts_global, npts_local, &
                                                            shift_local
      LOGICAL                                            :: my_use_virial
      LOGICAL(KIND=C_BOOL)                               :: my_compute_tau, orthorhombic
      REAL(KIND=dp), DIMENSION(3, 2), TARGET             :: forces
      REAL(KIND=dp), DIMENSION(3, 3, 2), TARGET          :: virials
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid
      TYPE(C_PTR)                                        :: a_hdab_cptr, forces_cptr, hadb_cptr, &
                                                            hdab_cptr, pab_cptr, virials_cptr
      INTERFACE
         SUBROUTINE grid_ref_integrate_pgf_product_c(orthorhombic, compute_tau, &
                                                     border_mask, &
                                                     la_max, la_min, lb_max, lb_min, &
                                                     zeta, zetb, dh, dh_inv, ra, rab, &
                                                     npts_global, npts_local, shift_local, border_width, &
                                                     radius, o1, o2, n1, n2, grid, hab, pab, &
                                                     forces, virials, hdab, hadb, a_hdab) &
            BIND(C, name="grid_ref_integrate_pgf_product")
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_BOOL
            LOGICAL(KIND=C_BOOL), VALUE               :: orthorhombic
            LOGICAL(KIND=C_BOOL), VALUE               :: compute_tau
            INTEGER(KIND=C_INT), VALUE                :: border_mask
            INTEGER(KIND=C_INT), VALUE                :: la_max
            INTEGER(KIND=C_INT), VALUE                :: la_min
            INTEGER(KIND=C_INT), VALUE                :: lb_max
            INTEGER(KIND=C_INT), VALUE                :: lb_min
            REAL(KIND=C_DOUBLE), VALUE                :: zeta
            REAL(KIND=C_DOUBLE), VALUE                :: zetb
            TYPE(C_PTR), VALUE                        :: dh
            TYPE(C_PTR), VALUE                        :: dh_inv
            TYPE(C_PTR), VALUE                        :: ra
            TYPE(C_PTR), VALUE                        :: rab
            TYPE(C_PTR), VALUE                        :: npts_global
            TYPE(C_PTR), VALUE                        :: npts_local
            TYPE(C_PTR), VALUE                        :: shift_local
            TYPE(C_PTR), VALUE                        :: border_width
            REAL(KIND=C_DOUBLE), VALUE                :: radius
            INTEGER(KIND=C_INT), VALUE                :: o1
            INTEGER(KIND=C_INT), VALUE                :: o2
            INTEGER(KIND=C_INT), VALUE                :: n1
            INTEGER(KIND=C_INT), VALUE                :: n2
            TYPE(C_PTR), VALUE                        :: grid
            TYPE(C_PTR), VALUE                        :: hab
            TYPE(C_PTR), VALUE                        :: pab
            TYPE(C_PTR), VALUE                        :: forces
            TYPE(C_PTR), VALUE                        :: virials
            TYPE(C_PTR), VALUE                        :: hdab
            TYPE(C_PTR), VALUE                        :: hadb
            TYPE(C_PTR), VALUE                        :: a_hdab
         END SUBROUTINE grid_ref_integrate_pgf_product_c
      END INTERFACE

      MARK_USED(cube_info)
      MARK_USED(cell)

      IF (radius == 0.0_dp) THEN
         RETURN
      END IF

      border_mask = 0
      IF (PRESENT(use_subpatch)) THEN
         IF (use_subpatch) THEN
            CPASSERT(PRESENT(subpatch_pattern))
            border_mask = IAND(63, NOT(subpatch_pattern))  ! invert last 6 bits
         END IF
      END IF

      ! When true then 0.5 * (nabla x_a).(v(r) nabla x_b) is computed.
      IF (PRESENT(compute_tau)) THEN
         my_compute_tau = LOGICAL(compute_tau, C_BOOL)
      ELSE
         my_compute_tau = .FALSE.
      END IF

      IF (PRESENT(use_virial)) THEN
         my_use_virial = use_virial
      ELSE
         my_use_virial = .FALSE.
      END IF

      IF (calculate_forces) THEN
         CPASSERT(PRESENT(pab))
         pab_cptr = C_LOC(pab(1, 1))
         forces(:, :) = 0.0_dp
         forces_cptr = C_LOC(forces(1, 1))
      ELSE
         pab_cptr = C_NULL_PTR
         forces_cptr = C_NULL_PTR
      END IF

      IF (calculate_forces .AND. my_use_virial) THEN
         virials(:, :, :) = 0.0_dp
         virials_cptr = C_LOC(virials(1, 1, 1))
      ELSE
         virials_cptr = C_NULL_PTR
      END IF

      IF (calculate_forces .AND. PRESENT(hdab)) THEN
         hdab_cptr = C_LOC(hdab(1, 1, 1))
      ELSE
         hdab_cptr = C_NULL_PTR
      END IF

      IF (calculate_forces .AND. PRESENT(hadb)) THEN
         hadb_cptr = C_LOC(hadb(1, 1, 1))
      ELSE
         hadb_cptr = C_NULL_PTR
      END IF

      IF (calculate_forces .AND. my_use_virial .AND. PRESENT(a_hdab)) THEN
         a_hdab_cptr = C_LOC(a_hdab(1, 1, 1, 1))
      ELSE
         a_hdab_cptr = C_NULL_PTR
      END IF

      orthorhombic = LOGICAL(rsgrid%desc%orthorhombic, C_BOOL)

      CALL get_rsgrid_properties(rsgrid, npts_global=npts_global, &
                                 npts_local=npts_local, &
                                 shift_local=shift_local, &
                                 border_width=border_width)

      grid(1:, 1:, 1:) => rsgrid%r(:, :, :)  ! pointer assignment

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(rsgrid%desc%dh))
      CPASSERT(IS_CONTIGUOUS(rsgrid%desc%dh_inv))
      CPASSERT(IS_CONTIGUOUS(ra))
      CPASSERT(IS_CONTIGUOUS(rab))
      CPASSERT(IS_CONTIGUOUS(npts_global))
      CPASSERT(IS_CONTIGUOUS(npts_local))
      CPASSERT(IS_CONTIGUOUS(shift_local))
      CPASSERT(IS_CONTIGUOUS(border_width))
      CPASSERT(IS_CONTIGUOUS(grid))
      CPASSERT(IS_CONTIGUOUS(hab))
      CPASSERT(IS_CONTIGUOUS(forces))
      CPASSERT(IS_CONTIGUOUS(virials))
      IF (PRESENT(pab)) THEN
         CPASSERT(IS_CONTIGUOUS(pab))
      END IF
      IF (PRESENT(hdab)) THEN
         CPASSERT(IS_CONTIGUOUS(hdab))
      END IF
      IF (PRESENT(a_hdab)) THEN
         CPASSERT(IS_CONTIGUOUS(a_hdab))
      END IF
#endif

      CALL grid_ref_integrate_pgf_product_c(orthorhombic=orthorhombic, &
                                            compute_tau=my_compute_tau, &
                                            border_mask=border_mask, &
                                            la_max=la_max, &
                                            la_min=la_min, &
                                            lb_max=lb_max, &
                                            lb_min=lb_min, &
                                            zeta=zeta, &
                                            zetb=zetb, &
                                            dh=C_LOC(rsgrid%desc%dh(1, 1)), &
                                            dh_inv=C_LOC(rsgrid%desc%dh_inv(1, 1)), &
                                            ra=C_LOC(ra(1)), &
                                            rab=C_LOC(rab(1)), &
                                            npts_global=C_LOC(npts_global(1)), &
                                            npts_local=C_LOC(npts_local(1)), &
                                            shift_local=C_LOC(shift_local(1)), &
                                            border_width=C_LOC(border_width(1)), &
                                            radius=radius, &
                                            o1=o1, &
                                            o2=o2, &
                                            n1=SIZE(hab, 1), &
                                            n2=SIZE(hab, 2), &
                                            grid=C_LOC(grid(1, 1, 1)), &
                                            hab=C_LOC(hab(1, 1)), &
                                            pab=pab_cptr, &
                                            forces=forces_cptr, &
                                            virials=virials_cptr, &
                                            hdab=hdab_cptr, &
                                            hadb=hadb_cptr, &
                                            a_hdab=a_hdab_cptr)

      IF (PRESENT(force_a) .AND. C_ASSOCIATED(forces_cptr)) &
         force_a = force_a + forces(:, 1)
      IF (PRESENT(force_b) .AND. C_ASSOCIATED(forces_cptr)) &
         force_b = force_b + forces(:, 2)
      IF (PRESENT(my_virial_a) .AND. C_ASSOCIATED(virials_cptr)) &
         my_virial_a = my_virial_a + virials(:, :, 1)
      IF (PRESENT(my_virial_b) .AND. C_ASSOCIATED(virials_cptr)) &
         my_virial_b = my_virial_b + virials(:, :, 2)

   END SUBROUTINE integrate_pgf_product

! **************************************************************************************************
!> \brief Helper routines for getting rsgrid properties and asserting underlying assumptions.
!> \param rsgrid ...
!> \param npts_global ...
!> \param npts_local ...
!> \param shift_local ...
!> \param border_width ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE get_rsgrid_properties(rsgrid, npts_global, npts_local, shift_local, border_width)
      TYPE(realspace_grid_type), INTENT(IN)              :: rsgrid
      INTEGER, DIMENSION(:)                              :: npts_global, npts_local, shift_local, &
                                                            border_width

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

      ! Convert rsgrid%desc%border and rsgrid%desc%perd into the more convenient border_width array.
      DO i = 1, 3
         IF (rsgrid%desc%perd(i) == 1) THEN
            ! Periodic meaning the grid in this direction is entriely present on every processor.
            CPASSERT(npts_local(i) == npts_global(i))
            CPASSERT(shift_local(i) == 0)
            ! No need for halo regions.
            border_width(i) = 0
         ELSE
            ! Not periodic meaning the grid in this direction is distributed among processors.
            CPASSERT(npts_local(i) <= npts_global(i))
            ! Check bounds of grid section that is owned by this processor.
            CPASSERT(rsgrid%lb_real(i) == rsgrid%lb_local(i) + rsgrid%desc%border)
            CPASSERT(rsgrid%ub_real(i) == rsgrid%ub_local(i) - rsgrid%desc%border)
            ! We have halo regions.
            border_width(i) = rsgrid%desc%border
         END IF
      END DO
   END SUBROUTINE get_rsgrid_properties

! **************************************************************************************************
!> \brief Allocates a basis set which can be passed to grid_create_task_list.
!> \param nset ...
!> \param nsgf ...
!> \param maxco ...
!> \param maxpgf ...
!> \param lmin ...
!> \param lmax ...
!> \param npgf ...
!> \param nsgf_set ...
!> \param first_sgf ...
!> \param sphi ...
!> \param zet ...
!> \param basis_set ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE grid_create_basis_set(nset, nsgf, maxco, maxpgf, &
                                    lmin, lmax, npgf, nsgf_set, first_sgf, sphi, zet, &
                                    basis_set)
      INTEGER, INTENT(IN)                                :: nset, nsgf, maxco, maxpgf
      INTEGER, DIMENSION(:), INTENT(IN), TARGET          :: lmin, lmax, npgf, nsgf_set
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: first_sgf
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), TARGET :: sphi, zet
      TYPE(grid_basis_set_type), INTENT(INOUT)           :: basis_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grid_create_basis_set'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(nset), TARGET                   :: my_first_sgf
      TYPE(C_PTR)                                        :: first_sgf_c, lmax_c, lmin_c, npgf_c, &
                                                            nsgf_set_c, sphi_c, zet_c
      INTERFACE
         SUBROUTINE grid_create_basis_set_c(nset, nsgf, maxco, maxpgf, &
                                            lmin, lmax, npgf, nsgf_set, first_sgf, sphi, zet, &
                                            basis_set) &
            BIND(C, name="grid_create_basis_set")
            IMPORT :: C_PTR, C_INT
            INTEGER(KIND=C_INT), VALUE                :: nset
            INTEGER(KIND=C_INT), VALUE                :: nsgf
            INTEGER(KIND=C_INT), VALUE                :: maxco
            INTEGER(KIND=C_INT), VALUE                :: maxpgf
            TYPE(C_PTR), VALUE                        :: lmin
            TYPE(C_PTR), VALUE                        :: lmax
            TYPE(C_PTR), VALUE                        :: npgf
            TYPE(C_PTR), VALUE                        :: nsgf_set
            TYPE(C_PTR), VALUE                        :: first_sgf
            TYPE(C_PTR), VALUE                        :: sphi
            TYPE(C_PTR), VALUE                        :: zet
            TYPE(C_PTR)                               :: basis_set
         END SUBROUTINE grid_create_basis_set_c
      END INTERFACE

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(lmin) == nset)
      CPASSERT(SIZE(lmin) == nset)
      CPASSERT(SIZE(lmax) == nset)
      CPASSERT(SIZE(npgf) == nset)
      CPASSERT(SIZE(nsgf_set) == nset)
      CPASSERT(SIZE(first_sgf, 2) == nset)
      CPASSERT(SIZE(sphi, 1) == maxco .AND. SIZE(sphi, 2) == nsgf)
      CPASSERT(SIZE(zet, 1) == maxpgf .AND. SIZE(zet, 2) == nset)
      CPASSERT(.NOT. C_ASSOCIATED(basis_set%c_ptr))

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(lmin))
      CPASSERT(IS_CONTIGUOUS(lmax))
      CPASSERT(IS_CONTIGUOUS(npgf))
      CPASSERT(IS_CONTIGUOUS(nsgf_set))
      CPASSERT(IS_CONTIGUOUS(my_first_sgf))
      CPASSERT(IS_CONTIGUOUS(sphi))
      CPASSERT(IS_CONTIGUOUS(zet))
#endif

      lmin_c = C_NULL_PTR
      lmax_c = C_NULL_PTR
      npgf_c = C_NULL_PTR
      nsgf_set_c = C_NULL_PTR
      first_sgf_c = C_NULL_PTR
      sphi_c = C_NULL_PTR
      zet_c = C_NULL_PTR

      ! Basis sets arrays can be empty, need to check before accessing the first element.
      IF (nset > 0) THEN
         lmin_c = C_LOC(lmin(1))
         lmax_c = C_LOC(lmax(1))
         npgf_c = C_LOC(npgf(1))
         nsgf_set_c = C_LOC(nsgf_set(1))
      END IF
      IF (SIZE(first_sgf) > 0) THEN
         my_first_sgf(:) = first_sgf(1, :)  ! make a contiguous copy
         first_sgf_c = C_LOC(my_first_sgf(1))
      END IF
      IF (SIZE(sphi) > 0) THEN
         sphi_c = C_LOC(sphi(1, 1))
      END IF
      IF (SIZE(zet) > 0) THEN
         zet_c = C_LOC(zet(1, 1))
      END IF

      CALL grid_create_basis_set_c(nset=nset, &
                                   nsgf=nsgf, &
                                   maxco=maxco, &
                                   maxpgf=maxpgf, &
                                   lmin=lmin_c, &
                                   lmax=lmax_c, &
                                   npgf=npgf_c, &
                                   nsgf_set=nsgf_set_c, &
                                   first_sgf=first_sgf_c, &
                                   sphi=sphi_c, &
                                   zet=zet_c, &
                                   basis_set=basis_set%c_ptr)
      CPASSERT(C_ASSOCIATED(basis_set%c_ptr))

      CALL timestop(handle)
   END SUBROUTINE grid_create_basis_set

! **************************************************************************************************
!> \brief Deallocates given basis set.
!> \param basis_set ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE grid_free_basis_set(basis_set)
      TYPE(grid_basis_set_type), INTENT(INOUT)           :: basis_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grid_free_basis_set'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE grid_free_basis_set_c(basis_set) &
            BIND(C, name="grid_free_basis_set")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: basis_set
         END SUBROUTINE grid_free_basis_set_c
      END INTERFACE

      CALL timeset(routineN, handle)

      CPASSERT(C_ASSOCIATED(basis_set%c_ptr))

      CALL grid_free_basis_set_c(basis_set%c_ptr)

      basis_set%c_ptr = C_NULL_PTR

      CALL timestop(handle)
   END SUBROUTINE grid_free_basis_set

! **************************************************************************************************
!> \brief Allocates a task list which can be passed to grid_collocate_task_list.
!> \param ntasks ...
!> \param natoms ...
!> \param nkinds ...
!> \param nblocks ...
!> \param block_offsets ...
!> \param atom_positions ...
!> \param atom_kinds ...
!> \param basis_sets ...
!> \param level_list ...
!> \param iatom_list ...
!> \param jatom_list ...
!> \param iset_list ...
!> \param jset_list ...
!> \param ipgf_list ...
!> \param jpgf_list ...
!> \param border_mask_list ...
!> \param block_num_list ...
!> \param radius_list ...
!> \param rab_list ...
!> \param rs_grids ...
!> \param task_list ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE grid_create_task_list(ntasks, natoms, nkinds, nblocks, &
                                    block_offsets, atom_positions, atom_kinds, basis_sets, &
                                    level_list, iatom_list, jatom_list, &
                                    iset_list, jset_list, ipgf_list, jpgf_list, &
                                    border_mask_list, block_num_list, &
                                    radius_list, rab_list, rs_grids, task_list)

      INTEGER, INTENT(IN)                                :: ntasks, natoms, nkinds, nblocks
      INTEGER, DIMENSION(:), INTENT(IN), TARGET          :: block_offsets
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), TARGET :: atom_positions
      INTEGER, DIMENSION(:), INTENT(IN), TARGET          :: atom_kinds
      TYPE(grid_basis_set_type), DIMENSION(:), &
         INTENT(IN), TARGET                              :: basis_sets
      INTEGER, DIMENSION(:), INTENT(IN), TARGET          :: level_list, iatom_list, jatom_list, &
                                                            iset_list, jset_list, ipgf_list, &
                                                            jpgf_list, border_mask_list, &
                                                            block_num_list
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), TARGET    :: radius_list
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), TARGET :: rab_list
      TYPE(realspace_grid_p_type), DIMENSION(:), &
         INTENT(IN)                                      :: rs_grids
      TYPE(grid_task_list_type), INTENT(INOUT)           :: task_list

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grid_create_task_list'

      INTEGER                                            :: handle, ikind, ilevel, nlevels
      INTEGER, ALLOCATABLE, DIMENSION(:, :), TARGET      :: border_width, npts_global, npts_local, &
                                                            shift_local
      LOGICAL(KIND=C_BOOL)                               :: orthorhombic
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         TARGET                                          :: dh, dh_inv
      TYPE(C_PTR) :: block_num_list_c, block_offsets_c, border_mask_list_c, iatom_list_c, &
         ipgf_list_c, iset_list_c, jatom_list_c, jpgf_list_c, jset_list_c, level_list_c, &
         rab_list_c, radius_list_c
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:), TARGET     :: basis_sets_c
      TYPE(realspace_grid_type), POINTER                 :: rsgrid
      INTERFACE
         SUBROUTINE grid_create_task_list_c(orthorhombic, &
                                            ntasks, nlevels, natoms, nkinds, nblocks, &
                                            block_offsets, atom_positions, atom_kinds, basis_sets, &
                                            level_list, iatom_list, jatom_list, &
                                            iset_list, jset_list, ipgf_list, jpgf_list, &
                                            border_mask_list, block_num_list, &
                                            radius_list, rab_list, &
                                            npts_global, npts_local, shift_local, &
                                            border_width, dh, dh_inv, task_list) &
            BIND(C, name="grid_create_task_list")
            IMPORT :: C_PTR, C_INT, C_BOOL
            LOGICAL(KIND=C_BOOL), VALUE               :: orthorhombic
            INTEGER(KIND=C_INT), VALUE                :: ntasks
            INTEGER(KIND=C_INT), VALUE                :: nlevels
            INTEGER(KIND=C_INT), VALUE                :: natoms
            INTEGER(KIND=C_INT), VALUE                :: nkinds
            INTEGER(KIND=C_INT), VALUE                :: nblocks
            TYPE(C_PTR), VALUE                        :: block_offsets
            TYPE(C_PTR), VALUE                        :: atom_positions
            TYPE(C_PTR), VALUE                        :: atom_kinds
            TYPE(C_PTR), VALUE                        :: basis_sets
            TYPE(C_PTR), VALUE                        :: level_list
            TYPE(C_PTR), VALUE                        :: iatom_list
            TYPE(C_PTR), VALUE                        :: jatom_list
            TYPE(C_PTR), VALUE                        :: iset_list
            TYPE(C_PTR), VALUE                        :: jset_list
            TYPE(C_PTR), VALUE                        :: ipgf_list
            TYPE(C_PTR), VALUE                        :: jpgf_list
            TYPE(C_PTR), VALUE                        :: border_mask_list
            TYPE(C_PTR), VALUE                        :: block_num_list
            TYPE(C_PTR), VALUE                        :: radius_list
            TYPE(C_PTR), VALUE                        :: rab_list
            TYPE(C_PTR), VALUE                        :: npts_global
            TYPE(C_PTR), VALUE                        :: npts_local
            TYPE(C_PTR), VALUE                        :: shift_local
            TYPE(C_PTR), VALUE                        :: border_width
            TYPE(C_PTR), VALUE                        :: dh
            TYPE(C_PTR), VALUE                        :: dh_inv
            TYPE(C_PTR)                               :: task_list
         END SUBROUTINE grid_create_task_list_c
      END INTERFACE

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(block_offsets) == nblocks)
      CPASSERT(SIZE(atom_positions, 1) == 3 .AND. SIZE(atom_positions, 2) == natoms)
      CPASSERT(SIZE(atom_kinds) == natoms)
      CPASSERT(SIZE(basis_sets) == nkinds)
      CPASSERT(SIZE(level_list) == ntasks)
      CPASSERT(SIZE(iatom_list) == ntasks)
      CPASSERT(SIZE(jatom_list) == ntasks)
      CPASSERT(SIZE(iset_list) == ntasks)
      CPASSERT(SIZE(jset_list) == ntasks)
      CPASSERT(SIZE(ipgf_list) == ntasks)
      CPASSERT(SIZE(jpgf_list) == ntasks)
      CPASSERT(SIZE(border_mask_list) == ntasks)
      CPASSERT(SIZE(block_num_list) == ntasks)
      CPASSERT(SIZE(radius_list) == ntasks)
      CPASSERT(SIZE(rab_list, 1) == 3 .AND. SIZE(rab_list, 2) == ntasks)

      ALLOCATE (basis_sets_c(nkinds))
      DO ikind = 1, nkinds
         basis_sets_c(ikind) = basis_sets(ikind)%c_ptr
      END DO

      nlevels = SIZE(rs_grids)
      CPASSERT(nlevels > 0)
      orthorhombic = LOGICAL(rs_grids(1)%rs_grid%desc%orthorhombic, C_BOOL)

      ALLOCATE (npts_global(3, nlevels), npts_local(3, nlevels))
      ALLOCATE (shift_local(3, nlevels), border_width(3, nlevels))
      ALLOCATE (dh(3, 3, nlevels), dh_inv(3, 3, nlevels))
      DO ilevel = 1, nlevels
         rsgrid => rs_grids(ilevel)%rs_grid
         CALL get_rsgrid_properties(rsgrid=rsgrid, &
                                    npts_global=npts_global(:, ilevel), &
                                    npts_local=npts_local(:, ilevel), &
                                    shift_local=shift_local(:, ilevel), &
                                    border_width=border_width(:, ilevel))
         CPASSERT(rsgrid%desc%orthorhombic .EQV. orthorhombic)  ! should be the same for all levels
         dh(:, :, ilevel) = rsgrid%desc%dh(:, :)
         dh_inv(:, :, ilevel) = rsgrid%desc%dh_inv(:, :)
      END DO

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(block_offsets))
      CPASSERT(IS_CONTIGUOUS(atom_positions))
      CPASSERT(IS_CONTIGUOUS(atom_kinds))
      CPASSERT(IS_CONTIGUOUS(basis_sets))
      CPASSERT(IS_CONTIGUOUS(level_list))
      CPASSERT(IS_CONTIGUOUS(iatom_list))
      CPASSERT(IS_CONTIGUOUS(jatom_list))
      CPASSERT(IS_CONTIGUOUS(iset_list))
      CPASSERT(IS_CONTIGUOUS(jset_list))
      CPASSERT(IS_CONTIGUOUS(ipgf_list))
      CPASSERT(IS_CONTIGUOUS(jpgf_list))
      CPASSERT(IS_CONTIGUOUS(border_mask_list))
      CPASSERT(IS_CONTIGUOUS(block_num_list))
      CPASSERT(IS_CONTIGUOUS(radius_list))
      CPASSERT(IS_CONTIGUOUS(rab_list))
      CPASSERT(IS_CONTIGUOUS(npts_global))
      CPASSERT(IS_CONTIGUOUS(npts_local))
      CPASSERT(IS_CONTIGUOUS(shift_local))
      CPASSERT(IS_CONTIGUOUS(border_width))
      CPASSERT(IS_CONTIGUOUS(dh))
      CPASSERT(IS_CONTIGUOUS(dh_inv))
#endif

      IF (ntasks > 0) THEN
         block_offsets_c = C_LOC(block_offsets(1))
         level_list_c = C_LOC(level_list(1))
         iatom_list_c = C_LOC(iatom_list(1))
         jatom_list_c = C_LOC(jatom_list(1))
         iset_list_c = C_LOC(iset_list(1))
         jset_list_c = C_LOC(jset_list(1))
         ipgf_list_c = C_LOC(ipgf_list(1))
         jpgf_list_c = C_LOC(jpgf_list(1))
         border_mask_list_c = C_LOC(border_mask_list(1))
         block_num_list_c = C_LOC(block_num_list(1))
         radius_list_c = C_LOC(radius_list(1))
         rab_list_c = C_LOC(rab_list(1, 1))
      ELSE
         ! Without tasks the lists are empty and there is no first element to call C_LOC on.
         block_offsets_c = C_NULL_PTR
         level_list_c = C_NULL_PTR
         iatom_list_c = C_NULL_PTR
         jatom_list_c = C_NULL_PTR
         iset_list_c = C_NULL_PTR
         jset_list_c = C_NULL_PTR
         ipgf_list_c = C_NULL_PTR
         jpgf_list_c = C_NULL_PTR
         border_mask_list_c = C_NULL_PTR
         block_num_list_c = C_NULL_PTR
         radius_list_c = C_NULL_PTR
         rab_list_c = C_NULL_PTR
      END IF

      !If task_list%c_ptr is already allocated, then its memory will be reused or freed.
      CALL grid_create_task_list_c(orthorhombic=orthorhombic, &
                                   ntasks=ntasks, &
                                   nlevels=nlevels, &
                                   natoms=natoms, &
                                   nkinds=nkinds, &
                                   nblocks=nblocks, &
                                   block_offsets=block_offsets_c, &
                                   atom_positions=C_LOC(atom_positions(1, 1)), &
                                   atom_kinds=C_LOC(atom_kinds(1)), &
                                   basis_sets=C_LOC(basis_sets_c(1)), &
                                   level_list=level_list_c, &
                                   iatom_list=iatom_list_c, &
                                   jatom_list=jatom_list_c, &
                                   iset_list=iset_list_c, &
                                   jset_list=jset_list_c, &
                                   ipgf_list=ipgf_list_c, &
                                   jpgf_list=jpgf_list_c, &
                                   border_mask_list=border_mask_list_c, &
                                   block_num_list=block_num_list_c, &
                                   radius_list=radius_list_c, &
                                   rab_list=rab_list_c, &
                                   npts_global=C_LOC(npts_global(1, 1)), &
                                   npts_local=C_LOC(npts_local(1, 1)), &
                                   shift_local=C_LOC(shift_local(1, 1)), &
                                   border_width=C_LOC(border_width(1, 1)), &
                                   dh=C_LOC(dh(1, 1, 1)), &
                                   dh_inv=C_LOC(dh_inv(1, 1, 1)), &
                                   task_list=task_list%c_ptr)

      CPASSERT(C_ASSOCIATED(task_list%c_ptr))

      CALL timestop(handle)
   END SUBROUTINE grid_create_task_list

! **************************************************************************************************
!> \brief Deallocates given task list, basis_sets have to be freed separately.
!> \param task_list ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE grid_free_task_list(task_list)
      TYPE(grid_task_list_type), INTENT(INOUT)           :: task_list

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grid_free_task_list'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE grid_free_task_list_c(task_list) &
            BIND(C, name="grid_free_task_list")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: task_list
         END SUBROUTINE grid_free_task_list_c
      END INTERFACE

      CALL timeset(routineN, handle)

      IF (C_ASSOCIATED(task_list%c_ptr)) THEN
         CALL grid_free_task_list_c(task_list%c_ptr)
      END IF

      task_list%c_ptr = C_NULL_PTR

      CALL timestop(handle)
   END SUBROUTINE grid_free_task_list

! **************************************************************************************************
!> \brief Collocate all tasks of in given list onto given grids.
!> \param task_list ...
!> \param ga_gb_function ...
!> \param pab_blocks ...
!> \param rs_grids ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE grid_collocate_task_list(task_list, ga_gb_function, pab_blocks, rs_grids)
      TYPE(grid_task_list_type), INTENT(IN)              :: task_list
      INTEGER, INTENT(IN)                                :: ga_gb_function
      TYPE(offload_buffer_type), INTENT(IN)              :: pab_blocks
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_grids

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grid_collocate_task_list'

      INTEGER                                            :: handle, ilevel, nlevels
      INTEGER, ALLOCATABLE, DIMENSION(:, :), TARGET      :: npts_local
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:), TARGET     :: grids_c
      TYPE(realspace_grid_type), POINTER                 :: rsgrid
      INTERFACE
         SUBROUTINE grid_collocate_task_list_c(task_list, func, nlevels, &
                                               npts_local, pab_blocks, grids) &
            BIND(C, name="grid_collocate_task_list")
            IMPORT :: C_PTR, C_INT, C_BOOL
            TYPE(C_PTR), VALUE                        :: task_list
            INTEGER(KIND=C_INT), VALUE                :: func
            INTEGER(KIND=C_INT), VALUE                :: nlevels
            TYPE(C_PTR), VALUE                        :: npts_local
            TYPE(C_PTR), VALUE                        :: pab_blocks
            TYPE(C_PTR), VALUE                        :: grids
         END SUBROUTINE grid_collocate_task_list_c
      END INTERFACE

      CALL timeset(routineN, handle)

      nlevels = SIZE(rs_grids)
      CPASSERT(nlevels > 0)

      ALLOCATE (grids_c(nlevels))
      ALLOCATE (npts_local(3, nlevels))
      DO ilevel = 1, nlevels
         rsgrid => rs_grids(ilevel)%rs_grid
         npts_local(:, ilevel) = rsgrid%ub_local - rsgrid%lb_local + 1
         grids_c(ilevel) = rsgrid%buffer%c_ptr
      END DO

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(npts_local))
      CPASSERT(IS_CONTIGUOUS(grids_c))
#endif

      CPASSERT(C_ASSOCIATED(task_list%c_ptr))
      CPASSERT(C_ASSOCIATED(pab_blocks%c_ptr))

      CALL grid_collocate_task_list_c(task_list=task_list%c_ptr, &
                                      func=ga_gb_function, &
                                      nlevels=nlevels, &
                                      npts_local=C_LOC(npts_local(1, 1)), &
                                      pab_blocks=pab_blocks%c_ptr, &
                                      grids=C_LOC(grids_c(1)))

      CALL timestop(handle)
   END SUBROUTINE grid_collocate_task_list

! **************************************************************************************************
!> \brief Integrate all tasks of in given list from given grids.
!> \param task_list ...
!> \param compute_tau ...
!> \param calculate_forces ...
!> \param calculate_virial ...
!> \param pab_blocks ...
!> \param rs_grids ...
!> \param hab_blocks ...
!> \param forces ...
!> \param virial ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE grid_integrate_task_list(task_list, compute_tau, calculate_forces, calculate_virial, &
                                       pab_blocks, rs_grids, hab_blocks, forces, virial)
      TYPE(grid_task_list_type), INTENT(IN)              :: task_list
      LOGICAL, INTENT(IN)                                :: compute_tau, calculate_forces, &
                                                            calculate_virial
      TYPE(offload_buffer_type), INTENT(IN)              :: pab_blocks
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_grids
      TYPE(offload_buffer_type), INTENT(INOUT)           :: hab_blocks
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         TARGET                                          :: forces
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         TARGET                                          :: virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grid_integrate_task_list'

      INTEGER                                            :: handle, ilevel, nlevels
      INTEGER, ALLOCATABLE, DIMENSION(:, :), TARGET      :: npts_local
      TYPE(C_PTR)                                        :: forces_c, virial_c
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:), TARGET     :: grids_c
      TYPE(realspace_grid_type), POINTER                 :: rsgrid
      INTERFACE
         SUBROUTINE grid_integrate_task_list_c(task_list, compute_tau, natoms, &
                                               nlevels, npts_local, &
                                               pab_blocks, grids, hab_blocks, forces, virial) &
            BIND(C, name="grid_integrate_task_list")
            IMPORT :: C_PTR, C_INT, C_BOOL
            TYPE(C_PTR), VALUE                        :: task_list
            LOGICAL(KIND=C_BOOL), VALUE               :: compute_tau
            INTEGER(KIND=C_INT), VALUE                :: natoms
            INTEGER(KIND=C_INT), VALUE                :: nlevels
            TYPE(C_PTR), VALUE                        :: npts_local
            TYPE(C_PTR), VALUE                        :: pab_blocks
            TYPE(C_PTR), VALUE                        :: grids
            TYPE(C_PTR), VALUE                        :: hab_blocks
            TYPE(C_PTR), VALUE                        :: forces
            TYPE(C_PTR), VALUE                        :: virial
         END SUBROUTINE grid_integrate_task_list_c
      END INTERFACE

      CALL timeset(routineN, handle)

      nlevels = SIZE(rs_grids)
      CPASSERT(nlevels > 0)

      ALLOCATE (grids_c(nlevels))
      ALLOCATE (npts_local(3, nlevels))
      DO ilevel = 1, nlevels
         rsgrid => rs_grids(ilevel)%rs_grid
         npts_local(:, ilevel) = rsgrid%ub_local - rsgrid%lb_local + 1
         grids_c(ilevel) = rsgrid%buffer%c_ptr
      END DO

      IF (calculate_forces) THEN
         forces_c = C_LOC(forces(1, 1))
      ELSE
         forces_c = C_NULL_PTR
      END IF

      IF (calculate_virial) THEN
         virial_c = C_LOC(virial(1, 1))
      ELSE
         virial_c = C_NULL_PTR
      END IF

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(npts_local))
      CPASSERT(IS_CONTIGUOUS(grids_c))
      CPASSERT(IS_CONTIGUOUS(forces))
      CPASSERT(IS_CONTIGUOUS(virial))
#endif

      CPASSERT(SIZE(forces, 1) == 3)
      CPASSERT(C_ASSOCIATED(task_list%c_ptr))
      CPASSERT(C_ASSOCIATED(hab_blocks%c_ptr))
      CPASSERT(C_ASSOCIATED(pab_blocks%c_ptr) .OR. .NOT. calculate_forces)
      CPASSERT(C_ASSOCIATED(pab_blocks%c_ptr) .OR. .NOT. calculate_virial)

      CALL grid_integrate_task_list_c(task_list=task_list%c_ptr, &
                                      compute_tau=LOGICAL(compute_tau, C_BOOL), &
                                      natoms=SIZE(forces, 2), &
                                      nlevels=nlevels, &
                                      npts_local=C_LOC(npts_local(1, 1)), &
                                      pab_blocks=pab_blocks%c_ptr, &
                                      grids=C_LOC(grids_c(1)), &
                                      hab_blocks=hab_blocks%c_ptr, &
                                      forces=forces_c, &
                                      virial=virial_c)

      CALL timestop(handle)
   END SUBROUTINE grid_integrate_task_list

! **************************************************************************************************
!> \brief Initialize grid library
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE grid_library_init()
      INTERFACE
         SUBROUTINE grid_library_init_c() BIND(C, name="grid_library_init")
         END SUBROUTINE grid_library_init_c
      END INTERFACE

      CALL grid_library_init_c()

   END SUBROUTINE grid_library_init

! **************************************************************************************************
!> \brief Finalize grid library
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE grid_library_finalize()
      INTERFACE
         SUBROUTINE grid_library_finalize_c() BIND(C, name="grid_library_finalize")
         END SUBROUTINE grid_library_finalize_c
      END INTERFACE

      CALL grid_library_finalize_c()

   END SUBROUTINE grid_library_finalize

! **************************************************************************************************
!> \brief Configures the grid library
!> \param backend : backend to be used for collocate/integrate, possible values are REF, CPU, GPU
!> \param validate : if set to true, compare the results of all backend to the reference backend
!> \param apply_cutoff : apply a spherical cutoff before collocating or integrating. Only relevant for CPU backend
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE grid_library_set_config(backend, validate, apply_cutoff)
      INTEGER, INTENT(IN)                                :: backend
      LOGICAL, INTENT(IN)                                :: validate, apply_cutoff

      INTERFACE
         SUBROUTINE grid_library_set_config_c(backend, validate, apply_cutoff) &
            BIND(C, name="grid_library_set_config")
            IMPORT :: C_INT, C_BOOL
            INTEGER(KIND=C_INT), VALUE                :: backend
            LOGICAL(KIND=C_BOOL), VALUE               :: validate
            LOGICAL(KIND=C_BOOL), VALUE               :: apply_cutoff
         END SUBROUTINE grid_library_set_config_c
      END INTERFACE

      CALL grid_library_set_config_c(backend=backend, &
                                     validate=LOGICAL(validate, C_BOOL), &
                                     apply_cutoff=LOGICAL(apply_cutoff, C_BOOL))

   END SUBROUTINE grid_library_set_config

! **************************************************************************************************
!> \brief Print grid library statistics
!> \param mpi_comm ...
!> \param output_unit ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE grid_library_print_stats(mpi_comm, output_unit)
      INTEGER, INTENT(IN)                                :: mpi_comm, output_unit

      INTERFACE
         SUBROUTINE grid_library_print_stats_c(mpi_sum_func, mpi_comm, print_func, output_unit) &
            BIND(C, name="grid_library_print_stats")
            IMPORT :: C_FUNPTR, C_INT
            TYPE(C_FUNPTR), VALUE                     :: mpi_sum_func
            INTEGER(KIND=C_INT), VALUE                :: mpi_comm
            TYPE(C_FUNPTR), VALUE                     :: print_func
            INTEGER(KIND=C_INT), VALUE                :: output_unit
         END SUBROUTINE grid_library_print_stats_c
      END INTERFACE

      ! Since Fortran units and mpi groups can't be used from C, we pass function pointers instead.
      CALL grid_library_print_stats_c(mpi_sum_func=C_FUNLOC(mpi_sum_func), &
                                      mpi_comm=mpi_comm, &
                                      print_func=C_FUNLOC(print_func), &
                                      output_unit=output_unit)

   END SUBROUTINE grid_library_print_stats

! **************************************************************************************************
!> \brief Callback to run mpi_sum on a Fortran MPI communicator.
!> \param number ...
!> \param mpi_comm ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE mpi_sum_func(number, mpi_comm) BIND(C)
      INTEGER(KIND=C_LONG), INTENT(INOUT)                :: number
      INTEGER(KIND=C_INT), INTENT(IN), VALUE             :: mpi_comm

      CALL mp_sum(number, mpi_comm)
   END SUBROUTINE mpi_sum_func

! **************************************************************************************************
!> \brief Callback to write to a Fortran output unit.
!> \param message ...
!> \param output_unit ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE print_func(message, output_unit) BIND(C)
      CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN)          :: message(*)
      INTEGER(KIND=C_INT), INTENT(IN), VALUE             :: output_unit

      CHARACTER(LEN=1000)                                :: buffer
      INTEGER                                            :: nchars

      IF (output_unit <= 0) &
         RETURN

      ! Convert C char array into Fortran string.
      nchars = strlcpy_c2f(buffer, message)

      ! Print the message.
      WRITE (output_unit, FMT="(A)", ADVANCE="NO") buffer(1:nchars)
   END SUBROUTINE print_func

END MODULE grid_api
