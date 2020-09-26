!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE grid_integrate
   USE ISO_C_BINDING,                   ONLY: C_BOOL,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_PTR
   USE cell_types,                      ONLY: cell_type
   USE cube_utils,                      ONLY: cube_info_type
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
   USE realspace_grid_types,            ONLY: realspace_grid_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grid_integrate'

! *** Don't include this routines directly, use the interface to
! *** qs_integrate_potential

   PUBLIC :: integrate_pgf_product

CONTAINS

! **************************************************************************************************
!> \brief Low level integrate of primitive gaussian functions.
!> \param la_max ...
!> \param la_min ...
!> \param zeta ...
!> \param lb_max ...
!> \param lb_min ...
!> \param zetb ...
!> \param rscale ...
!> \param ra ...
!> \param rab ...
!> \param rsgrid ...
!> \param radius ...
!> \param vab ...
!> \param use_subpatch ...
!> \param subpatch_pattern ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE grid_integrate_vab(la_max, la_min, zeta, lb_max, lb_min, zetb, rscale, ra, rab, &
                                 rsgrid, radius, vab, use_subpatch, subpatch_pattern)

      INTEGER, INTENT(IN)                                :: la_max, la_min
      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, lb_min
      REAL(KIND=dp), INTENT(IN)                          :: zetb, rscale
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), TARGET    :: ra, rab
      TYPE(realspace_grid_type), INTENT(IN)              :: rsgrid
      REAL(KIND=dp), INTENT(IN)                          :: radius
      REAL(KIND=dp), DIMENSION(:, :), TARGET             :: vab
      LOGICAL, OPTIONAL                                  :: use_subpatch
      INTEGER, INTENT(IN), OPTIONAL                      :: subpatch_pattern

      CHARACTER(len=*), PARAMETER :: routineN = 'grid_integrate_vab', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: border_mask
      INTEGER, DIMENSION(3), TARGET                      :: border_width, npts_global, npts_local, &
                                                            shift_local
      LOGICAL(KIND=C_BOOL)                               :: orthorhombic
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid
      INTERFACE
         SUBROUTINE grid_ref_integrate_vab_c(orthorhombic, &
                                             border_mask, &
                                             la_max, la_min, lb_max, lb_min, &
                                             zeta, zetb, rscale, dh, dh_inv, ra, rab, &
                                             npts_global, npts_local, shift_local, border_width, &
                                             radius, grid, vab) &
            BIND(C, name="grid_ref_integrate_vab")
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_BOOL
            LOGICAL(KIND=C_BOOL), VALUE               :: orthorhombic
            INTEGER(KIND=C_INT), VALUE                :: border_mask
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
            TYPE(C_PTR), VALUE                        :: grid
            TYPE(C_PTR), VALUE                        :: vab
         END SUBROUTINE grid_ref_integrate_vab_c
      END INTERFACE

      border_mask = 0
      IF (PRESENT(use_subpatch)) THEN
         IF (use_subpatch) THEN
            CPASSERT(PRESENT(subpatch_pattern))
            border_mask = IAND(63, NOT(subpatch_pattern))  ! invert last 6 bits
         END IF
      END IF

      orthorhombic = LOGICAL(rsgrid%desc%orthorhombic, C_BOOL)

      CALL get_rsgrid_properties(rsgrid, npts_global=npts_global, &
                                 npts_local=npts_local, &
                                 shift_local=shift_local, &
                                 border_width=border_width)

      grid(1:, 1:, 1:) => rsgrid%r(:, :, :)  ! pointer assignment

      CPASSERT(LBOUND(vab, 1) == 1 .AND. LBOUND(vab, 2) == 1)

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
      CPASSERT(IS_CONTIGUOUS(vab))
#endif

      CALL grid_ref_integrate_vab_c(orthorhombic=orthorhombic, &
                                    border_mask=border_mask, &
                                    la_max=la_max, &
                                    la_min=la_min, &
                                    lb_max=lb_max, &
                                    lb_min=lb_min, &
                                    zeta=zeta, &
                                    zetb=zetb, &
                                    rscale=rscale, &
                                    dh=C_LOC(rsgrid%desc%dh(1, 1)), &
                                    dh_inv=C_LOC(rsgrid%desc%dh_inv(1, 1)), &
                                    ra=C_LOC(ra(1)), &
                                    rab=C_LOC(rab(1)), &
                                    npts_global=C_LOC(npts_global(1)), &
                                    npts_local=C_LOC(npts_local(1)), &
                                    shift_local=C_LOC(shift_local(1)), &
                                    border_width=C_LOC(border_width(1)), &
                                    radius=radius, &
                                    grid=C_LOC(grid(1, 1, 1)), &
                                    vab=C_LOC(vab(1, 1)))

   END SUBROUTINE grid_integrate_vab

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
         ENDIF
      ENDDO
   END SUBROUTINE get_rsgrid_properties

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
!> \param hdab ...
!> \param hadb ...
!> \param force_a ...
!> \param force_b ...
!> \param compute_tau ...
!> \param use_virial ...
!> \param my_virial_a ...
!> \param my_virial_b ...
!> \param a_hdab ...
!> \param use_subpatch ...
!> \param subpatch_pattern ...
! **************************************************************************************************
   SUBROUTINE integrate_pgf_product(la_max, zeta, la_min, &
                                    lb_max, zetb, lb_min, &
                                    ra, rab, rsgrid, cell, &
                                    cube_info, hab, pab, o1, o2, &
                                    radius, &
                                    calculate_forces, hdab, hadb, force_a, force_b, &
                                    compute_tau, &
                                    use_virial, my_virial_a, &
                                    my_virial_b, a_hdab, use_subpatch, subpatch_pattern)

      INTEGER, INTENT(IN)                                :: la_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zetb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rab
      TYPE(realspace_grid_type), POINTER                 :: rsgrid
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type), INTENT(IN)                   :: cube_info
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hab
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: pab
      INTEGER, INTENT(IN)                                :: o1, o2
      REAL(KIND=dp), INTENT(IN)                          :: radius
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: hdab, hadb
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT), &
         OPTIONAL                                        :: force_a, force_b
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_tau, use_virial
      REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL           :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:, :, :, :), OPTIONAL, &
         POINTER                                         :: a_hdab
      LOGICAL, OPTIONAL                                  :: use_subpatch
      INTEGER, INTENT(IN), OPTIONAL                      :: subpatch_pattern

      INTEGER                                            :: ax, ay, az, bx, by, bz, ico, j, jco, la, &
                                                            la_max_local, la_min_local, lb, &
                                                            lb_max_local, lb_min_local
      LOGICAL                                            :: my_compute_tau, my_use_virial, &
                                                            subpatch_integrate
      REAL(KIND=dp)                                      :: axpm0, der_a(3), der_b(3), ftza, ftzb, &
                                                            pabval
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: vab

      MARK_USED(cube_info)
      MARK_USED(cell)

      subpatch_integrate = .FALSE.

      IF (PRESENT(use_subpatch)) THEN
         IF (use_subpatch) THEN
            subpatch_integrate = .TRUE.
            CPASSERT(PRESENT(subpatch_pattern))
         ENDIF
      ENDIF

      IF (PRESENT(use_virial)) THEN
         my_use_virial = use_virial
      ELSE
         my_use_virial = .FALSE.
      ENDIF

      ! my_compute_tau defaults to .FALSE.
      ! IF (.true.) it will compute 0.5 * (nabla x_a).(v(r) nabla x_b)
      IF (PRESENT(compute_tau)) THEN
         my_compute_tau = compute_tau
      ELSE
         my_compute_tau = .FALSE.
      ENDIF

      IF (calculate_forces) THEN
         la_max_local = la_max + 1 ! needed for the derivative of the gaussian, unimportant which one
         la_min_local = MAX(la_min - 1, 0) ! just in case the la_min,lb_min is not zero
         lb_min_local = MAX(lb_min - 1, 0)
         lb_max_local = lb_max
         IF (my_use_virial) THEN
            la_max_local = la_max_local + 1
            lb_max_local = lb_max_local + 1
         ENDIF
      ELSE
         la_max_local = la_max
         la_min_local = la_min
         lb_min_local = lb_min
         lb_max_local = lb_max
      END IF

      IF (my_compute_tau) THEN
         la_max_local = la_max_local + 1
         lb_max_local = lb_max_local + 1
         la_min_local = MAX(la_min_local - 1, 0)
         lb_min_local = MAX(lb_min_local - 1, 0)
      ENDIF

      IF (radius == 0.0_dp) THEN
         RETURN
      ENDIF

      ALLOCATE (vab(ncoset(la_max_local), ncoset(lb_max_local)))
      vab = 0.0_dp

      CALL grid_integrate_vab(la_max=la_max_local, &
                              la_min=la_min_local, &
                              zeta=zeta, &
                              lb_max=lb_max_local, &
                              lb_min=lb_min_local, &
                              zetb=zetb, &
                              rscale=1.0_dp, &
                              ra=ra, &
                              rab=rab, &
                              rsgrid=rsgrid, &
                              radius=radius, &
                              vab=vab, &
                              use_subpatch=use_subpatch, &
                              subpatch_pattern=subpatch_pattern)

!   *** vab contains all the information needed to find the elements of hab
!   *** and optionally of derivatives of these elements

      ftza = 2.0_dp*zeta
      ftzb = 2.0_dp*zetb

      DO la = la_min, la_max
         DO ax = 0, la
            DO ay = 0, la - ax
               az = la - ax - ay
               ico = coset(ax, ay, az)
               DO lb = lb_min, lb_max
                  DO bx = 0, lb
                     DO by = 0, lb - bx
                        bz = lb - bx - by
                        jco = coset(bx, by, bz)
                        IF (.NOT. my_compute_tau) THEN
                           axpm0 = vab(ico, jco)
                        ELSE
                           axpm0 = 0.5_dp*(ax*bx*vab(coset(MAX(ax - 1, 0), ay, az), coset(MAX(bx - 1, 0), by, bz)) + &
                                           ay*by*vab(coset(ax, MAX(ay - 1, 0), az), coset(bx, MAX(by - 1, 0), bz)) + &
                                           az*bz*vab(coset(ax, ay, MAX(az - 1, 0)), coset(bx, by, MAX(bz - 1, 0))) &
                                           - ftza*bx*vab(coset(ax + 1, ay, az), coset(MAX(bx - 1, 0), by, bz)) &
                                           - ftza*by*vab(coset(ax, ay + 1, az), coset(bx, MAX(by - 1, 0), bz)) &
                                           - ftza*bz*vab(coset(ax, ay, az + 1), coset(bx, by, MAX(bz - 1, 0))) &
                                           - ax*ftzb*vab(coset(MAX(ax - 1, 0), ay, az), coset(bx + 1, by, bz)) &
                                           - ay*ftzb*vab(coset(ax, MAX(ay - 1, 0), az), coset(bx, by + 1, bz)) &
                                           - az*ftzb*vab(coset(ax, ay, MAX(az - 1, 0)), coset(bx, by, bz + 1)) + &
                                           ftza*ftzb*vab(coset(ax + 1, ay, az), coset(bx + 1, by, bz)) + &
                                           ftza*ftzb*vab(coset(ax, ay + 1, az), coset(bx, by + 1, bz)) + &
                                           ftza*ftzb*vab(coset(ax, ay, az + 1), coset(bx, by, bz + 1)))
                        ENDIF
                        hab(o1 + ico, o2 + jco) = hab(o1 + ico, o2 + jco) + axpm0
                        IF (calculate_forces .AND. PRESENT(force_a)) THEN
                           IF (my_compute_tau) THEN
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*ax*bx
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                MAX(ax - 1, 0), ay, az, MAX(bx - 1, 0), by, bz, vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    MAX(ax - 1, 0), ay, az, MAX(bx - 1, 0), by, bz, vab)
                              ENDIF
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*ay*by
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                ax, MAX(ay - 1, 0), az, bx, MAX(by - 1, 0), bz, vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    ax, MAX(ay - 1, 0), az, bx, MAX(by - 1, 0), bz, vab)
                              ENDIF
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*az*bz
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                ax, ay, MAX(az - 1, 0), bx, by, MAX(bz - 1, 0), vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    ax, ay, MAX(az - 1, 0), bx, by, MAX(bz - 1, 0), vab)
                              ENDIF
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*(-ftza*bx)
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                ax + 1, ay, az, MAX(bx - 1, 0), by, bz, vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    ax + 1, ay, az, MAX(bx - 1, 0), by, bz, vab)
                              ENDIF
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*(-ftza*by)
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                ax, ay + 1, az, bx, MAX(by - 1, 0), bz, vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    ax, ay + 1, az, bx, MAX(by - 1, 0), bz, vab)
                              ENDIF
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*(-ftza*bz)
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                ax, ay, az + 1, bx, by, MAX(bz - 1, 0), vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    ax, ay, az + 1, bx, by, MAX(bz - 1, 0), vab)
                              ENDIF
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*(-ax*ftzb)
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                MAX(ax - 1, 0), ay, az, bx + 1, by, bz, vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    MAX(ax - 1, 0), ay, az, bx + 1, by, bz, vab)
                              ENDIF
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*(-ay*ftzb)
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                ax, MAX(ay - 1, 0), az, bx, by + 1, bz, vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    ax, MAX(ay - 1, 0), az, bx, by + 1, bz, vab)
                              ENDIF
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*(-az*ftzb)
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                ax, ay, MAX(az - 1, 0), bx, by, bz + 1, vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    ax, ay, MAX(az - 1, 0), bx, by, bz + 1, vab)
                              ENDIF
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*(ftza*ftzb)
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                ax + 1, ay, az, bx + 1, by, bz, vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    ax + 1, ay, az, bx + 1, by, bz, vab)
                              ENDIF
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*(ftza*ftzb)
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                ax, ay + 1, az, bx, by + 1, bz, vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    ax, ay + 1, az, bx, by + 1, bz, vab)
                              ENDIF
                              pabval = pab(o1 + ico, o2 + jco)*0.5_dp*(ftza*ftzb)
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                ax, ay, az + 1, bx, by, bz + 1, vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    ax, ay, az + 1, bx, by, bz + 1, vab)
                              ENDIF
                           ELSE
                              pabval = pab(o1 + ico, o2 + jco)
                              CALL force_update(force_a, force_b, rab, pabval, ftza, ftzb, &
                                                ax, ay, az, bx, by, bz, vab)
                              IF (my_use_virial) THEN
                                 CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                    ax, ay, az, bx, by, bz, vab)
                              ENDIF
                           ENDIF
                        END IF
                        IF (calculate_forces .AND. PRESENT(hdab)) THEN
                           der_a(1:3) = 0.0_dp
                           der_b(1:3) = 0.0_dp
                           CALL hab_derivatives(der_a, der_b, rab, ftza, ftzb, &
                                                ax, ay, az, bx, by, bz, vab)
                           hdab(1:3, o1 + ico, o2 + jco) = der_a(1:3)
                           hadb(1:3, o1 + ico, o2 + jco) = der_b(1:3)
                           pabval = 1.0_dp
                           IF (my_use_virial .AND. PRESENT(a_hdab)) THEN
                              my_virial_a = 0.0_dp
                              my_virial_b = 0.0_dp
                              CALL virial_update(my_virial_a, my_virial_b, rab, pabval, ftza, ftzb, &
                                                 ax, ay, az, bx, by, bz, vab)
                              DO j = 1, 3
                                 a_hdab(1:3, j, o1 + ico, o2 + jco) = &
                                    a_hdab(1:3, j, o1 + ico, o2 + jco) + my_virial_a(1:3, j)
                              END DO
                           ENDIF
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (vab)

   END SUBROUTINE integrate_pgf_product

! **************************************************************************************************
!> \brief given a set of matrix elements, perform the correct contraction to obtain the virial
!> \param my_virial_a ...
!> \param my_virial_b ...
!> \param rab ...
!> \param pab ...
!> \param ftza ...
!> \param ftzb ...
!> \param ax ...
!> \param ay ...
!> \param az ...
!> \param bx ...
!> \param by ...
!> \param bz ...
!> \param vab ...
! **************************************************************************************************
   SUBROUTINE virial_update(my_virial_a, my_virial_b, rab, pab, &
                            ftza, ftzb, ax, ay, az, bx, by, bz, vab)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), INTENT(IN)                          :: pab, ftza, ftzb
      INTEGER, INTENT(IN)                                :: ax, ay, az, bx, by, bz
      REAL(KIND=dp)                                      :: vab(:, :)

      my_virial_a(1, 1) = my_virial_a(1, 1) &
                          + pab*ftza*vab(coset(ax + 2, ay, az), coset(bx, by, bz)) &
                          - pab*REAL(ax, dp)*vab(coset(MAX(0, ax - 1) + 1, ay, az), coset(bx, by, bz))
      my_virial_a(1, 2) = my_virial_a(1, 2) &
                          + pab*ftza*vab(coset(ax + 1, ay + 1, az), coset(bx, by, bz)) &
                          - pab*REAL(ax, dp)*vab(coset(MAX(0, ax - 1), ay + 1, az), coset(bx, by, bz))
      my_virial_a(1, 3) = my_virial_a(1, 3) &
                          + pab*ftza*vab(coset(ax + 1, ay, az + 1), coset(bx, by, bz)) &
                          - pab*REAL(ax, dp)*vab(coset(MAX(0, ax - 1), ay, az + 1), coset(bx, by, bz))
      my_virial_a(2, 1) = my_virial_a(2, 1) &
                          + pab*ftza*vab(coset(ax + 1, ay + 1, az), coset(bx, by, bz)) &
                          - pab*REAL(ay, dp)*vab(coset(ax + 1, MAX(0, ay - 1), az), coset(bx, by, bz))
      my_virial_a(2, 2) = my_virial_a(2, 2) &
                          + pab*ftza*vab(coset(ax, ay + 2, az), coset(bx, by, bz)) &
                          - pab*REAL(ay, dp)*vab(coset(ax, MAX(0, ay - 1) + 1, az), coset(bx, by, bz))
      my_virial_a(2, 3) = my_virial_a(2, 3) &
                          + pab*ftza*vab(coset(ax, ay + 1, az + 1), coset(bx, by, bz)) &
                          - pab*REAL(ay, dp)*vab(coset(ax, MAX(0, ay - 1), az + 1), coset(bx, by, bz))
      my_virial_a(3, 1) = my_virial_a(3, 1) &
                          + pab*ftza*vab(coset(ax + 1, ay, az + 1), coset(bx, by, bz)) &
                          - pab*REAL(az, dp)*vab(coset(ax + 1, ay, MAX(0, az - 1)), coset(bx, by, bz))
      my_virial_a(3, 2) = my_virial_a(3, 2) &
                          + pab*ftza*vab(coset(ax, ay + 1, az + 1), coset(bx, by, bz)) &
                          - pab*REAL(az, dp)*vab(coset(ax, ay + 1, MAX(0, az - 1)), coset(bx, by, bz))
      my_virial_a(3, 3) = my_virial_a(3, 3) &
                          + pab*ftza*vab(coset(ax, ay, az + 2), coset(bx, by, bz)) &
                          - pab*REAL(az, dp)*vab(coset(ax, ay, MAX(0, az - 1) + 1), coset(bx, by, bz))

      my_virial_b(1, 1) = my_virial_b(1, 1) + pab*ftzb*( &
                          vab(coset(ax + 2, ay, az), coset(bx, by, bz)) &
                          - vab(coset(ax + 1, ay, az), coset(bx, by, bz))*rab(1) &
                          - vab(coset(ax + 1, ay, az), coset(bx, by, bz))*rab(1) &
                          + vab(coset(ax, ay, az), coset(bx, by, bz))*rab(1)*rab(1)) &
                          - pab*REAL(bx, dp)*vab(coset(ax, ay, az), coset(MAX(0, bx - 1) + 1, by, bz))
      my_virial_b(1, 2) = my_virial_b(1, 2) + pab*ftzb*( &
                          vab(coset(ax + 1, ay + 1, az), coset(bx, by, bz)) &
                          - vab(coset(ax, ay + 1, az), coset(bx, by, bz))*rab(1) &
                          - vab(coset(ax + 1, ay, az), coset(bx, by, bz))*rab(2) &
                          + vab(coset(ax, ay, az), coset(bx, by, bz))*rab(1)*rab(2)) &
                          - pab*REAL(bx, dp)*vab(coset(ax, ay, az), coset(MAX(0, bx - 1), by + 1, bz))
      my_virial_b(1, 3) = my_virial_b(1, 3) + pab*ftzb*( &
                          vab(coset(ax + 1, ay, az + 1), coset(bx, by, bz)) &
                          - vab(coset(ax, ay, az + 1), coset(bx, by, bz))*rab(1) &
                          - vab(coset(ax + 1, ay, az), coset(bx, by, bz))*rab(3) &
                          + vab(coset(ax, ay, az), coset(bx, by, bz))*rab(1)*rab(3)) &
                          - pab*REAL(bx, dp)*vab(coset(ax, ay, az), coset(MAX(0, bx - 1), by, bz + 1))
      my_virial_b(2, 1) = my_virial_b(2, 1) + pab*ftzb*( &
                          vab(coset(ax + 1, ay + 1, az), coset(bx, by, bz)) &
                          - vab(coset(ax + 1, ay, az), coset(bx, by, bz))*rab(2) &
                          - vab(coset(ax, ay + 1, az), coset(bx, by, bz))*rab(1) &
                          + vab(coset(ax, ay, az), coset(bx, by, bz))*rab(2)*rab(1)) &
                          - pab*REAL(by, dp)*vab(coset(ax, ay, az), coset(bx + 1, MAX(0, by - 1), bz))
      my_virial_b(2, 2) = my_virial_b(2, 2) + pab*ftzb*( &
                          vab(coset(ax, ay + 2, az), coset(bx, by, bz)) &
                          - vab(coset(ax, ay + 1, az), coset(bx, by, bz))*rab(2) &
                          - vab(coset(ax, ay + 1, az), coset(bx, by, bz))*rab(2) &
                          + vab(coset(ax, ay, az), coset(bx, by, bz))*rab(2)*rab(2)) &
                          - pab*REAL(by, dp)*vab(coset(ax, ay, az), coset(bx, MAX(0, by - 1) + 1, bz))
      my_virial_b(2, 3) = my_virial_b(2, 3) + pab*ftzb*( &
                          vab(coset(ax, ay + 1, az + 1), coset(bx, by, bz)) &
                          - vab(coset(ax, ay, az + 1), coset(bx, by, bz))*rab(2) &
                          - vab(coset(ax, ay + 1, az), coset(bx, by, bz))*rab(3) &
                          + vab(coset(ax, ay, az), coset(bx, by, bz))*rab(2)*rab(3)) &
                          - pab*REAL(by, dp)*vab(coset(ax, ay, az), coset(bx, MAX(0, by - 1), bz + 1))
      my_virial_b(3, 1) = my_virial_b(3, 1) + pab*ftzb*( &
                          vab(coset(ax + 1, ay, az + 1), coset(bx, by, bz)) &
                          - vab(coset(ax + 1, ay, az), coset(bx, by, bz))*rab(3) &
                          - vab(coset(ax, ay, az + 1), coset(bx, by, bz))*rab(1) &
                          + vab(coset(ax, ay, az), coset(bx, by, bz))*rab(3)*rab(1)) &
                          - pab*REAL(bz, dp)*vab(coset(ax, ay, az), coset(bx + 1, by, MAX(0, bz - 1)))
      my_virial_b(3, 2) = my_virial_b(3, 2) + pab*ftzb*( &
                          vab(coset(ax, ay + 1, az + 1), coset(bx, by, bz)) &
                          - vab(coset(ax, ay + 1, az), coset(bx, by, bz))*rab(3) &
                          - vab(coset(ax, ay, az + 1), coset(bx, by, bz))*rab(2) &
                          + vab(coset(ax, ay, az), coset(bx, by, bz))*rab(3)*rab(2)) &
                          - pab*REAL(bz, dp)*vab(coset(ax, ay, az), coset(bx, by + 1, MAX(0, bz - 1)))
      my_virial_b(3, 3) = my_virial_b(3, 3) + pab*ftzb*( &
                          vab(coset(ax, ay, az + 2), coset(bx, by, bz)) &
                          - vab(coset(ax, ay, az + 1), coset(bx, by, bz))*rab(3) &
                          - vab(coset(ax, ay, az + 1), coset(bx, by, bz))*rab(3) &
                          + vab(coset(ax, ay, az), coset(bx, by, bz))*rab(3)*rab(3)) &
                          - pab*REAL(bz, dp)*vab(coset(ax, ay, az), coset(bx, by, MAX(0, bz - 1) + 1))

   END SUBROUTINE virial_update

! **************************************************************************************************
!> \brief given a bunch of matrix elements, performe the right contractions to obtain the forces
!> \param force_a ...
!> \param force_b ...
!> \param rab ...
!> \param pab ...
!> \param ftza ...
!> \param ftzb ...
!> \param ax ...
!> \param ay ...
!> \param az ...
!> \param bx ...
!> \param by ...
!> \param bz ...
!> \param vab ...
! **************************************************************************************************
   SUBROUTINE force_update(force_a, force_b, rab, pab, ftza, ftzb, ax, ay, az, bx, by, bz, vab)
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT)         :: force_a, force_b
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), INTENT(IN)                          :: pab, ftza, ftzb
      INTEGER, INTENT(IN)                                :: ax, ay, az, bx, by, bz
      REAL(KIND=dp)                                      :: vab(:, :)

      REAL(KIND=dp)                                      :: axm1, axp1, axpm0, aym1, ayp1, azm1, &
                                                            azp1, bxm1, bym1, bzm1

      axpm0 = vab(coset(ax, ay, az), coset(bx, by, bz))
      axp1 = vab(coset(ax + 1, ay, az), coset(bx, by, bz))
      axm1 = vab(coset(MAX(0, ax - 1), ay, az), coset(bx, by, bz))
      ayp1 = vab(coset(ax, ay + 1, az), coset(bx, by, bz))
      aym1 = vab(coset(ax, MAX(0, ay - 1), az), coset(bx, by, bz))
      azp1 = vab(coset(ax, ay, az + 1), coset(bx, by, bz))
      azm1 = vab(coset(ax, ay, MAX(0, az - 1)), coset(bx, by, bz))
      bxm1 = vab(coset(ax, ay, az), coset(MAX(0, bx - 1), by, bz))
      bym1 = vab(coset(ax, ay, az), coset(bx, MAX(0, by - 1), bz))
      bzm1 = vab(coset(ax, ay, az), coset(bx, by, MAX(0, bz - 1)))
      force_a(1) = force_a(1) + pab*(ftza*axp1 - REAL(ax, dp)*axm1)
      force_a(2) = force_a(2) + pab*(ftza*ayp1 - REAL(ay, dp)*aym1)
      force_a(3) = force_a(3) + pab*(ftza*azp1 - REAL(az, dp)*azm1)
      force_b(1) = force_b(1) + pab*(ftzb*(axp1 - rab(1)*axpm0) - REAL(bx, dp)*bxm1)
      force_b(2) = force_b(2) + pab*(ftzb*(ayp1 - rab(2)*axpm0) - REAL(by, dp)*bym1)
      force_b(3) = force_b(3) + pab*(ftzb*(azp1 - rab(3)*axpm0) - REAL(bz, dp)*bzm1)

   END SUBROUTINE force_update

! **************************************************************************************************
!> \brief given a bunch of matrix elements perform the right contractions to obtain the
!>      derivatives of the hab matirx
!> \param der_a ...
!> \param der_b ...
!> \param rab ...
!> \param ftza ...
!> \param ftzb ...
!> \param ax ...
!> \param ay ...
!> \param az ...
!> \param bx ...
!> \param by ...
!> \param bz ...
!> \param vab ...
! **************************************************************************************************
   SUBROUTINE hab_derivatives(der_a, der_b, rab, ftza, ftzb, ax, ay, az, bx, by, bz, vab)
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT)         :: der_a, der_b
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), INTENT(IN)                          :: ftza, ftzb
      INTEGER, INTENT(IN)                                :: ax, ay, az, bx, by, bz
      REAL(KIND=dp)                                      :: vab(:, :)

      REAL(KIND=dp)                                      :: axm1, axp1, axpm0, aym1, ayp1, azm1, &
                                                            azp1, bxm1, bym1, bzm1

      axpm0 = vab(coset(ax, ay, az), coset(bx, by, bz))
      axp1 = vab(coset(ax + 1, ay, az), coset(bx, by, bz))
      axm1 = vab(coset(MAX(0, ax - 1), ay, az), coset(bx, by, bz))
      ayp1 = vab(coset(ax, ay + 1, az), coset(bx, by, bz))
      aym1 = vab(coset(ax, MAX(0, ay - 1), az), coset(bx, by, bz))
      azp1 = vab(coset(ax, ay, az + 1), coset(bx, by, bz))
      azm1 = vab(coset(ax, ay, MAX(0, az - 1)), coset(bx, by, bz))
      bxm1 = vab(coset(ax, ay, az), coset(MAX(0, bx - 1), by, bz))
      bym1 = vab(coset(ax, ay, az), coset(bx, MAX(0, by - 1), bz))
      bzm1 = vab(coset(ax, ay, az), coset(bx, by, MAX(0, bz - 1)))
      der_a(1) = (ftza*axp1 - REAL(ax, dp)*axm1)
      der_a(2) = (ftza*ayp1 - REAL(ay, dp)*aym1)
      der_a(3) = (ftza*azp1 - REAL(az, dp)*azm1)
      der_b(1) = (ftzb*(axp1 - rab(1)*axpm0) - REAL(bx, dp)*bxm1)
      der_b(2) = (ftzb*(ayp1 - rab(2)*axpm0) - REAL(by, dp)*bym1)
      der_b(3) = (ftzb*(azp1 - rab(3)*axpm0) - REAL(bz, dp)*bzm1)

   END SUBROUTINE hab_derivatives

END MODULE grid_integrate
