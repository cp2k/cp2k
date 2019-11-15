!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
MODULE qs_grid_atom

   USE input_constants,                 ONLY: do_gapw_gcs,&
                                              do_gapw_gct,&
                                              do_gapw_log
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: pi
   USE memory_utilities,                ONLY: reallocate
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_grid_atom'

   TYPE grid_batch_type
      INTEGER                                             :: np
      REAL(KIND=dp), DIMENSION(3)                         :: rcenter
      REAL(KIND=dp)                                       :: rad
      REAL(dp), DIMENSION(:, :), ALLOCATABLE              :: rco
      REAL(dp), DIMENSION(:), ALLOCATABLE                 :: weight
   END TYPE grid_batch_type

   TYPE atom_integration_grid_type
      INTEGER                                             :: nr, na
      INTEGER                                             :: np, ntot
      INTEGER                                             :: lebedev_grid
      REAL(dp), DIMENSION(:), ALLOCATABLE                 :: rr
      REAL(dp), DIMENSION(:), ALLOCATABLE                 :: wr, wa
      INTEGER                                             :: nbatch
      TYPE(grid_batch_type), DIMENSION(:), ALLOCATABLE    :: batch
   END TYPE atom_integration_grid_type

   TYPE grid_atom_type
      INTEGER                         :: quadrature
      INTEGER                         :: nr, ng_sphere
      REAL(dp), DIMENSION(:), POINTER :: rad, rad2, &
                                         wr, wa, &
                                         azi, cos_azi, sin_azi, &
                                         pol, cos_pol, sin_pol, usin_azi
      REAL(dp), DIMENSION(:, :), &
         POINTER :: rad2l, oorad2l, weight
   END TYPE grid_atom_type

   PUBLIC :: allocate_grid_atom, create_grid_atom, deallocate_grid_atom
   PUBLIC :: grid_atom_type
   PUBLIC :: initialize_atomic_grid
   PUBLIC :: atom_integration_grid_type, deallocate_atom_int_grid

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief   Initialize components of the grid_atom_type structure
!> \param grid_atom ...
!> \date    03.11.2000
!> \author  MK
!> \author Matthias Krack (MK)
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE allocate_grid_atom(grid_atom)

      TYPE(grid_atom_type), POINTER                      :: grid_atom

      CHARACTER(len=*), PARAMETER :: routineN = 'allocate_grid_atom', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(grid_atom)) CALL deallocate_grid_atom(grid_atom)

      ALLOCATE (grid_atom)

      NULLIFY (grid_atom%rad)
      NULLIFY (grid_atom%rad2)
      NULLIFY (grid_atom%wr)
      NULLIFY (grid_atom%wa)
      NULLIFY (grid_atom%weight)
      NULLIFY (grid_atom%azi)
      NULLIFY (grid_atom%cos_azi)
      NULLIFY (grid_atom%sin_azi)
      NULLIFY (grid_atom%pol)
      NULLIFY (grid_atom%cos_pol)
      NULLIFY (grid_atom%sin_pol)
      NULLIFY (grid_atom%usin_azi)
      NULLIFY (grid_atom%rad2l)
      NULLIFY (grid_atom%oorad2l)

   END SUBROUTINE allocate_grid_atom

! **************************************************************************************************
!> \brief   Deallocate a Gaussian-type orbital (GTO) basis set data set.
!> \param grid_atom ...
!> \date    03.11.2000
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE deallocate_grid_atom(grid_atom)
      TYPE(grid_atom_type), POINTER                      :: grid_atom

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_grid_atom', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(grid_atom)) THEN

         IF (ASSOCIATED(grid_atom%rad)) THEN
            DEALLOCATE (grid_atom%rad)
         END IF

         IF (ASSOCIATED(grid_atom%rad2)) THEN
            DEALLOCATE (grid_atom%rad2)
         END IF

         IF (ASSOCIATED(grid_atom%wr)) THEN
            DEALLOCATE (grid_atom%wr)
         END IF

         IF (ASSOCIATED(grid_atom%wa)) THEN
            DEALLOCATE (grid_atom%wa)
         END IF

         IF (ASSOCIATED(grid_atom%weight)) THEN
            DEALLOCATE (grid_atom%weight)
         END IF

         IF (ASSOCIATED(grid_atom%azi)) THEN
            DEALLOCATE (grid_atom%azi)
         END IF

         IF (ASSOCIATED(grid_atom%cos_azi)) THEN
            DEALLOCATE (grid_atom%cos_azi)
         END IF

         IF (ASSOCIATED(grid_atom%sin_azi)) THEN
            DEALLOCATE (grid_atom%sin_azi)
         END IF

         IF (ASSOCIATED(grid_atom%pol)) THEN
            DEALLOCATE (grid_atom%pol)
         END IF

         IF (ASSOCIATED(grid_atom%cos_pol)) THEN
            DEALLOCATE (grid_atom%cos_pol)
         END IF

         IF (ASSOCIATED(grid_atom%sin_pol)) THEN
            DEALLOCATE (grid_atom%sin_pol)
         END IF

         IF (ASSOCIATED(grid_atom%usin_azi)) THEN
            DEALLOCATE (grid_atom%usin_azi)
         END IF

         IF (ASSOCIATED(grid_atom%rad2l)) THEN
            DEALLOCATE (grid_atom%rad2l)
         END IF

         IF (ASSOCIATED(grid_atom%oorad2l)) THEN
            DEALLOCATE (grid_atom%oorad2l)
         END IF

         DEALLOCATE (grid_atom)
      ELSE
         CALL cp_abort(__LOCATION__, &
                       "The pointer grid_atom is not associated and "// &
                       "cannot be deallocated")
      END IF
   END SUBROUTINE deallocate_grid_atom

! **************************************************************************************************
!> \brief ...
!> \param grid_atom ...
!> \param nr ...
!> \param na ...
!> \param llmax ...
!> \param ll ...
!> \param quadrature ...
! **************************************************************************************************
   SUBROUTINE create_grid_atom(grid_atom, nr, na, llmax, ll, quadrature)

      TYPE(grid_atom_type), POINTER                      :: grid_atom
      INTEGER, INTENT(IN)                                :: nr, na, llmax, ll, quadrature

      CHARACTER(len=*), PARAMETER :: routineN = 'create_grid_atom', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ia, ir, l
      REAL(dp)                                           :: cosia, pol
      REAL(dp), DIMENSION(:), POINTER                    :: rad, rad2, wr

      CALL timeset(routineN, handle)

      NULLIFY (rad, rad2, wr)

      IF (ASSOCIATED(grid_atom)) THEN

         ! Allocate the radial grid arrays
         CALL reallocate(grid_atom%rad, 1, nr)
         CALL reallocate(grid_atom%rad2, 1, nr)
         CALL reallocate(grid_atom%wr, 1, nr)
         CALL reallocate(grid_atom%wa, 1, na)
         CALL reallocate(grid_atom%weight, 1, na, 1, nr)
         CALL reallocate(grid_atom%azi, 1, na)
         CALL reallocate(grid_atom%cos_azi, 1, na)
         CALL reallocate(grid_atom%sin_azi, 1, na)
         CALL reallocate(grid_atom%pol, 1, na)
         CALL reallocate(grid_atom%cos_pol, 1, na)
         CALL reallocate(grid_atom%sin_pol, 1, na)
         CALL reallocate(grid_atom%usin_azi, 1, na)
         CALL reallocate(grid_atom%rad2l, 1, nr, 0, llmax + 1)
         CALL reallocate(grid_atom%oorad2l, 1, nr, 0, llmax + 1)

         ! Calculate the radial grid for this kind
         rad => grid_atom%rad
         rad2 => grid_atom%rad2
         wr => grid_atom%wr

         grid_atom%quadrature = quadrature
         CALL radial_grid(nr, rad, rad2, wr, quadrature)

         grid_atom%rad2l(:, 0) = 1._dp
         grid_atom%oorad2l(:, 0) = 1._dp
         DO l = 1, llmax + 1
            grid_atom%rad2l(:, l) = grid_atom%rad2l(:, l - 1)*rad(:)
            grid_atom%oorad2l(:, l) = grid_atom%oorad2l(:, l - 1)/rad(:)
         ENDDO

         IF (ll > 0) THEN
            grid_atom%wa(1:na) = 4._dp*pi*lebedev_grid(ll)%w(1:na)
            DO ir = 1, nr
               DO ia = 1, na
                  grid_atom%weight(ia, ir) = grid_atom%wr(ir)*grid_atom%wa(ia)
               END DO
            END DO

            DO ia = 1, na
               ! polar angle: pol = acos(r(3))
               cosia = lebedev_grid(ll)%r(3, ia)
               grid_atom%cos_pol(ia) = cosia
               ! azimuthal angle: pol = atan(r(2)/r(1))
               IF (ABS(lebedev_grid(ll)%r(2, ia)) < EPSILON(1.0_dp) .AND. &
                   ABS(lebedev_grid(ll)%r(1, ia)) < EPSILON(1.0_dp)) THEN
                  grid_atom%azi(ia) = 0.0_dp
               ELSE
                  grid_atom%azi(ia) = ATAN2(lebedev_grid(ll)%r(2, ia), lebedev_grid(ll)%r(1, ia))
               END IF
               grid_atom%cos_azi(ia) = COS(grid_atom%azi(ia))
               pol = ACOS(cosia)
               IF (grid_atom%sin_pol(ia) < 0.0_dp) pol = -pol
               grid_atom%pol(ia) = pol
               grid_atom%sin_pol(ia) = SIN(grid_atom%pol(ia))

               grid_atom%sin_azi(ia) = SIN(grid_atom%azi(ia))
               IF (ABS(grid_atom%sin_azi(ia)) > EPSILON(1.0_dp)) THEN
                  grid_atom%usin_azi(ia) = 1.0_dp/grid_atom%sin_azi(ia)
               ELSE
                  grid_atom%usin_azi(ia) = 1.0_dp
               END IF

            ENDDO

         END IF

      ELSE
         CPABORT("The pointer grid_atom is not associated")
      END IF

      CALL timestop(handle)

   END SUBROUTINE create_grid_atom

! **************************************************************************************************
!> \brief   Initialize atomic grid
!> \param   int_grid ...
!> \param nr ...
!> \param na ...
!> \param rmax ...
!> \param quadrature ...
!> \param iunit ...
!> \date    02.2018
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE initialize_atomic_grid(int_grid, nr, na, rmax, quadrature, iunit)
      TYPE(atom_integration_grid_type), POINTER          :: int_grid
      INTEGER, INTENT(IN)                                :: nr, na
      REAL(KIND=dp), INTENT(IN)                          :: rmax
      INTEGER, INTENT(IN), OPTIONAL                      :: quadrature, iunit

      CHARACTER(len=*), PARAMETER :: routineN = 'initialize_atomic_grid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ig, ir, ix, iy, iz, la, ll, my_quad, &
                                                            n1, n2, n3, nbatch, ng, no, np, ntot, &
                                                            nu, nx
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: icell
      REAL(KIND=dp)                                      :: ag, dd, dmax, r1, r2, r3
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rad, rad2, wa, wc, wr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rang, rco
      REAL(KIND=dp), DIMENSION(10)                       :: dco
      REAL(KIND=dp), DIMENSION(3)                        :: rm
      TYPE(atom_integration_grid_type), POINTER          :: igr

      ALLOCATE (igr)

      ! type of quadrature grid
      IF (PRESENT(quadrature)) THEN
         my_quad = quadrature
      ELSE
         my_quad = do_gapw_log
      END IF

      ! radial grid
      CPASSERT(nr > 1)
      ALLOCATE (rad(nr), rad2(nr), wr(nr))
      CALL radial_grid(nr, rad, rad2, wr, my_quad)
      !
      igr%nr = nr
      ALLOCATE (igr%rr(nr))
      ALLOCATE (igr%wr(nr))
      ! store grid points always in ascending order
      IF (rad(1) > rad(nr)) THEN
         DO ir = nr, 1, -1
            igr%rr(nr - ir + 1) = rad(ir)
            igr%wr(nr - ir + 1) = wr(ir)
         END DO
      ELSE
         igr%rr(1:nr) = rad(1:nr)
         igr%wr(1:nr) = wr(1:nr)
      END IF
      ! only include grid points smaller than rmax
      np = 0
      DO ir = 1, nr
         IF (igr%rr(ir) < rmax) THEN
            np = np + 1
            rad(np) = igr%rr(ir)
            wr(np) = igr%wr(ir)
         END IF
      END DO
      igr%np = np
      !
      ! angular grid
      CPASSERT(na > 1)
      ll = get_number_of_lebedev_grid(n=na)
      np = lebedev_grid(ll)%n
      la = lebedev_grid(ll)%l
      ALLOCATE (rang(3, np), wa(np))
      wa(1:na) = 4._dp*pi*lebedev_grid(ll)%w(1:np)
      rang(1:3, 1:np) = lebedev_grid(ll)%r(1:3, 1:np)
      igr%lebedev_grid = ll
      ALLOCATE (igr%wa(np))
      igr%na = np
      igr%wa(1:np) = wa(1:np)
      !
      ! total grid points
      ntot = igr%na*igr%np
      igr%ntot = ntot
      ALLOCATE (rco(3, ntot), wc(ntot))
      ig = 0
      DO ir = 1, igr%np
         DO ia = 1, igr%na
            ig = ig + 1
            rco(1:3, ig) = rang(1:3, ia)*rad(ir)
            wc(ig) = wa(ia)*wr(ir)
         END DO
      END DO
      ! grid for batches, odd number of cells
      ng = NINT((REAL(ntot, dp)/32._dp)**0.33333_dp)
      ng = ng + MOD(ng + 1, 2)
      ! avarage number of points along radial grid
      dco = 0.0_dp
      ag = REAL(igr%np, dp)/ng
      CPASSERT(SIZE(dco) >= (ng + 1)/2)
      DO ig = 1, ng, 2
         ir = MIN(NINT(ag)*ig, igr%np)
         ia = (ig + 1)/2
         dco(ia) = rad(ir)
      END DO
      ! batches
      ALLOCATE (icell(ntot))
      icell = 0
      nx = (ng - 1)/2
      DO ig = 1, ntot
         ix = grid_coord(rco(1, ig), dco, nx + 1) + nx
         iy = grid_coord(rco(2, ig), dco, nx + 1) + nx
         iz = grid_coord(rco(3, ig), dco, nx + 1) + nx
         icell(ig) = iz*ng*ng + iy*ng + ix + 1
      END DO
      !
      igr%nbatch = ng*ng*ng
      ALLOCATE (igr%batch(igr%nbatch))
      igr%batch(:)%np = 0
      DO ig = 1, ntot
         ia = icell(ig)
         igr%batch(ia)%np = igr%batch(ia)%np + 1
      END DO
      DO ig = 1, igr%nbatch
         np = igr%batch(ig)%np
         ALLOCATE (igr%batch(ig)%rco(3, np), igr%batch(ig)%weight(np))
         igr%batch(ig)%np = 0
      END DO
      DO ig = 1, ntot
         ia = icell(ig)
         igr%batch(ia)%np = igr%batch(ia)%np + 1
         np = igr%batch(ia)%np
         igr%batch(ia)%rco(1:3, np) = rco(1:3, ig)
         igr%batch(ia)%weight(np) = wc(ig)
      END DO
      !
      DEALLOCATE (rad, rad2, rang, wr, wa)
      DEALLOCATE (rco, wc, icell)
      !
      IF (ASSOCIATED(int_grid)) CALL deallocate_atom_int_grid(int_grid)
      ALLOCATE (int_grid)
      ALLOCATE (int_grid%rr(igr%nr), int_grid%wr(igr%nr), int_grid%wa(igr%na))
      int_grid%nr = igr%nr
      int_grid%na = igr%na
      int_grid%np = igr%np
      int_grid%ntot = igr%ntot
      int_grid%lebedev_grid = igr%lebedev_grid
      int_grid%rr(:) = igr%rr(:)
      int_grid%wr(:) = igr%wr(:)
      int_grid%wa(:) = igr%wa(:)
      !
      ! count batches
      nbatch = 0
      DO ig = 1, igr%nbatch
         IF (igr%batch(ig)%np == 0) THEN
            ! empty batch
         ELSE IF (igr%batch(ig)%np <= 48) THEN
            ! single batch
            nbatch = nbatch + 1
         ELSE
            ! multiple batches
            nbatch = nbatch + NINT(igr%batch(ig)%np/32._dp)
         END IF
      END DO
      int_grid%nbatch = nbatch
      ALLOCATE (int_grid%batch(nbatch))
      ! fill batches
      n1 = 0
      DO ig = 1, igr%nbatch
         IF (igr%batch(ig)%np == 0) THEN
            ! empty batch
         ELSE IF (igr%batch(ig)%np <= 48) THEN
            ! single batch
            n1 = n1 + 1
            np = igr%batch(ig)%np
            ALLOCATE (int_grid%batch(n1)%rco(3, np), int_grid%batch(n1)%weight(np))
            int_grid%batch(n1)%np = np
            int_grid%batch(n1)%rco(1:3, 1:np) = igr%batch(ig)%rco(1:3, 1:np)
            int_grid%batch(n1)%weight(1:np) = igr%batch(ig)%weight(1:np)
         ELSE
            ! multiple batches
            n2 = NINT(igr%batch(ig)%np/32._dp)
            n3 = igr%batch(ig)%np/n2
            DO ia = n1 + 1, n1 + n2
               nu = (ia - n1 - 1)*n3 + 1
               no = nu + n3 - 1
               IF (ia == n1 + n2) no = igr%batch(ig)%np
               np = no - nu + 1
               ALLOCATE (int_grid%batch(ia)%rco(3, np), int_grid%batch(ia)%weight(np))
               int_grid%batch(ia)%np = np
               int_grid%batch(ia)%rco(1:3, 1:np) = igr%batch(ig)%rco(1:3, nu:no)
               int_grid%batch(ia)%weight(1:np) = igr%batch(ig)%weight(nu:no)
            END DO
            n1 = n1 + n2
         END IF
      END DO
      CPASSERT(nbatch == n1)
      ! batch center and radius
      DO ig = 1, int_grid%nbatch
         np = int_grid%batch(ig)%np
         IF (np > 0) THEN
            rm(1) = SUM(int_grid%batch(ig)%rco(1, 1:np))
            rm(2) = SUM(int_grid%batch(ig)%rco(2, 1:np))
            rm(3) = SUM(int_grid%batch(ig)%rco(3, 1:np))
            rm(1:3) = rm(1:3)/REAL(np, KIND=dp)
         END IF
         int_grid%batch(ig)%rcenter(1:3) = rm(1:3)
         dmax = 0.0_dp
         DO ia = 1, np
            dd = SUM((int_grid%batch(ig)%rco(1:3, ia) - rm(1:3))**2)
            dmax = MAX(dd, dmax)
         END DO
         int_grid%batch(ig)%rad = SQRT(dmax)
      END DO
      !
      CALL deallocate_atom_int_grid(igr)
      !
      IF (PRESENT(iunit)) THEN
         IF (iunit > 0) THEN
            WRITE (iunit, "(/,A)") " Atomic Integration Grid Information"
            WRITE (iunit, "(A,T51,3I10)") "    Number of grid points [radial,angular,total]", &
               int_grid%np, int_grid%na, int_grid%ntot
            WRITE (iunit, "(A,T71,I10)") "    Lebedev grid number", int_grid%lebedev_grid
            WRITE (iunit, "(A,T61,F20.5)") "    Maximum of radial grid [Bohr]", &
               int_grid%rr(int_grid%np)
            nbatch = int_grid%nbatch
            WRITE (iunit, "(A,T71,I10)") "    Total number of gridpoint batches", nbatch
            n1 = int_grid%ntot
            n2 = 0
            n3 = NINT(REAL(int_grid%ntot, dp)/REAL(nbatch, dp))
            DO ig = 1, nbatch
               n1 = MIN(n1, int_grid%batch(ig)%np)
               n2 = MAX(n2, int_grid%batch(ig)%np)
            END DO
            WRITE (iunit, "(A,T51,3I10)") "    Number of grid points/batch [min,max,ave]", n1, n2, n3
            r1 = 1000._dp
            r2 = 0.0_dp
            r3 = 0.0_dp
            DO ig = 1, int_grid%nbatch
               r1 = MIN(r1, int_grid%batch(ig)%rad)
               r2 = MAX(r2, int_grid%batch(ig)%rad)
               r3 = r3 + int_grid%batch(ig)%rad
            END DO
            r3 = r3/REAL(ng*ng*ng, KIND=dp)
            WRITE (iunit, "(A,T51,3F10.2)") "    Batch radius (bohr) [min,max,ave]", r1, r2, r3
         END IF
      END IF

   END SUBROUTINE initialize_atomic_grid

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param dco ...
!> \param ng ...
!> \return ...
!> \retval igrid ...
! **************************************************************************************************
   FUNCTION grid_coord(x, dco, ng) RESULT(igrid)
      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: dco
      INTEGER, INTENT(IN)                                :: ng
      INTEGER                                            :: igrid

      INTEGER                                            :: ig
      REAL(KIND=dp)                                      :: xval

      xval = ABS(x)
      igrid = ng
      DO ig = 1, ng
         IF (xval <= dco(ig)) THEN
            igrid = ig - 1
            EXIT
         END IF
      END DO
      IF (x < 0.0_dp) igrid = -igrid
      CPASSERT(ABS(igrid) < ng)
   END FUNCTION grid_coord

! **************************************************************************************************
!> \brief ...
!> \param int_grid ...
! **************************************************************************************************
   SUBROUTINE deallocate_atom_int_grid(int_grid)
      TYPE(atom_integration_grid_type), POINTER          :: int_grid

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_atom_int_grid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ib

      IF (ASSOCIATED(int_grid)) THEN
         IF (ALLOCATED(int_grid%rr)) DEALLOCATE (int_grid%rr)
         IF (ALLOCATED(int_grid%wr)) DEALLOCATE (int_grid%wr)
         IF (ALLOCATED(int_grid%wa)) DEALLOCATE (int_grid%wa)
         ! batch
         IF (ALLOCATED(int_grid%batch)) THEN
            DO ib = 1, SIZE(int_grid%batch)
               IF (ALLOCATED(int_grid%batch(ib)%rco)) DEALLOCATE (int_grid%batch(ib)%rco)
               IF (ALLOCATED(int_grid%batch(ib)%weight)) DEALLOCATE (int_grid%batch(ib)%weight)
            END DO
            DEALLOCATE (int_grid%batch)
         END IF
         !
         DEALLOCATE (int_grid)
         NULLIFY (int_grid)
      END IF

   END SUBROUTINE deallocate_atom_int_grid
! **************************************************************************************************
!> \brief   Generate a radial grid with n points by a quadrature rule.
!> \param n ...
!> \param r ...
!> \param r2 ...
!> \param wr ...
!> \param radial_quadrature ...
!> \date    20.09.1999
!> \par Literature
!>           - A. D. Becke, J. Chem. Phys. 88, 2547 (1988)
!>           - J. M. Perez-Jorda, A. D. Becke and E. San-Fabian,
!>             J. Chem. Phys. 100, 6520 (1994)
!>           - M. Krack and A. M. Koester, J. Chem. Phys. 108, 3226 (1998)
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE radial_grid(n, r, r2, wr, radial_quadrature)

      INTEGER, INTENT(IN)                                :: n
      REAL(dp), DIMENSION(:), INTENT(INOUT)              :: r, r2, wr
      INTEGER, INTENT(IN)                                :: radial_quadrature

      CHARACTER(len=*), PARAMETER :: routineN = 'radial_grid', routineP = moduleN//':'//routineN

      INTEGER                                            :: i
      REAL(dp)                                           :: cost, f, sint, sint2, t, w, x

      f = pi/REAL(n + 1, dp)

      SELECT CASE (radial_quadrature)

      CASE (do_gapw_gcs)

         !     *** Gauss-Chebyshev quadrature formula of the second kind ***
         !     *** u [-1,+1] -> r [0,infinity] =>  r = (1 + u)/(1 - u)   ***

         DO i = 1, n
            t = REAL(i, dp)*f
            x = COS(t)
            w = f*SIN(t)**2
            r(i) = (1.0_dp + x)/(1.0_dp - x)
            r2(i) = r(i)**2
            wr(i) = w/SQRT(1.0_dp - x**2)
            wr(i) = 2.0_dp*wr(i)*r2(i)/(1.0_dp - x)**2
         END DO

      CASE (do_gapw_gct)

         !     *** Transformed Gauss-Chebyshev quadrature formula of the second kind ***
         !     *** u [-1,+1] -> r [0,infinity] => r = (1 + u)/(1 - u)                ***

         DO i = 1, n
            t = REAL(i, dp)*f
            cost = COS(t)
            sint = SIN(t)
            sint2 = sint**2
            x = REAL(2*i - n - 1, dp)/REAL(n + 1, dp) - &
                2.0_dp*(1.0_dp + 2.0_dp*sint2/3.0_dp)*cost*sint/pi
            w = 16.0_dp*sint2**2/REAL(3*(n + 1), dp)
            r(n + 1 - i) = (1.0_dp + x)/(1.0_dp - x)
            r2(n + 1 - i) = r(n + 1 - i)**2
            wr(n + 1 - i) = 2.0_dp*w*r2(n + 1 - i)/(1.0_dp - x)**2
         END DO

      CASE (do_gapw_log)

         !     *** Transformed Gauss-Chebyshev quadrature formula of the second kind ***
         !     *** u [-1,+1] -> r [0,infinity] => r = ln(2/(1 - u))/ln(2)            ***

         DO i = 1, n
            t = REAL(i, dp)*f
            cost = COS(t)
            sint = SIN(t)
            sint2 = sint**2
            x = REAL(2*i - n - 1, dp)/REAL(n + 1, dp) - &
                2.0_dp*(1.0_dp + 2.0_dp*sint2/3.0_dp)*cost*sint/pi
            w = 16.0_dp*sint2**2/REAL(3*(n + 1), dp)
            r(n + 1 - i) = LOG(2.0_dp/(1.0_dp - x))/LOG(2.0_dp)
            r2(n + 1 - i) = r(n + 1 - i)**2
            wr(n + 1 - i) = w*r2(n + 1 - i)/(LOG(2.0_dp)*(1.0_dp - x))
         END DO

      CASE DEFAULT

         CPABORT("Invalid radial quadrature type specified")

      END SELECT

   END SUBROUTINE radial_grid

END MODULE qs_grid_atom
