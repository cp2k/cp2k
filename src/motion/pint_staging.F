!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Data type and methods dealing with PI calcs in staging coordinates
!> \author fawzi
!> \par    History
!>         2006-02 created
!>         2006-11 modified so it might actually work [hforbert]
!>         2009-04-07 moved from pint_types module to a separate file [lwalewski]
! **************************************************************************************************
MODULE pint_staging
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE pint_types,                      ONLY: staging_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_staging'

   PUBLIC :: staging_env_create, staging_release, &
             staging_init_masses, &
             staging_x2u, staging_u2x, staging_f2uf, &
             staging_calc_uf_h

CONTAINS

! ***************************************************************************
!> \brief creates the data needed for a staging transformation
!> \param staging_env ...
!> \param staging_section ...
!> \param p ...
!> \param kT ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE staging_env_create(staging_env, staging_section, p, kT)
      TYPE(staging_env_type), INTENT(OUT)                :: staging_env
      TYPE(section_vals_type), POINTER                   :: staging_section
      INTEGER, INTENT(in)                                :: p
      REAL(kind=dp), INTENT(in)                          :: kT

      CALL section_vals_val_get(staging_section, "j", i_val=staging_env%j)
      CALL section_vals_val_get(staging_section, "Q_end", i_val=staging_env%j)
      staging_env%p = p
      staging_env%nseg = staging_env%p/staging_env%j

      staging_env%w_p = SQRT(REAL(staging_env%p, dp))*kT
      staging_env%w_j = kT*SQRT(REAL(staging_env%nseg, dp))
      staging_env%Q_stage = kT/staging_env%w_p**2
      IF (staging_env%Q_end <= 0._dp) THEN
         staging_env%Q_end = staging_env%j*staging_env%Q_stage
      END IF
   END SUBROUTINE staging_env_create

! ***************************************************************************
!> \brief releases the staging environment, kept for symmetry reasons with staging_env_create
!> \param staging_env the staging_env to release
!> \author Fawzi Mohamed
! **************************************************************************************************
   ELEMENTAL SUBROUTINE staging_release(staging_env)
      TYPE(staging_env_type), INTENT(IN)                 :: staging_env

      MARK_USED(staging_env)
   END SUBROUTINE staging_release

! ***************************************************************************
!> \brief initializes the masses and fictitious masses compatibly with the
!>      staging information
!> \param staging_env the definition of the staging transformation
!> \param mass *input* the masses of the particles
!> \param mass_beads masses of the beads
!> \param mass_fict the fictitious masses
!> \param Q masses of the nose thermostats
!> \par History
!>      11.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   PURE SUBROUTINE staging_init_masses(staging_env, mass, mass_beads, mass_fict, &
                                       Q)
      TYPE(staging_env_type), INTENT(IN)                 :: staging_env
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: mass
      REAL(kind=dp), DIMENSION(:, :), INTENT(out), &
         OPTIONAL                                        :: mass_beads, mass_fict
      REAL(kind=dp), DIMENSION(:), INTENT(out), OPTIONAL :: Q

      INTEGER                                            :: i, iat, ib, iseg
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: scal

      IF (PRESENT(Q)) THEN
         Q = staging_env%Q_stage
         DO i = 1, staging_env%p, staging_env%j
            Q(i) = staging_env%Q_end
         END DO
      END IF
      IF (PRESENT(mass_beads) .OR. PRESENT(mass_fict)) THEN

         ALLOCATE (scal(staging_env%p))
         DO iseg = 1, staging_env%nseg
            DO i = 1, staging_env%j ! check order!!!
               scal(staging_env%j*(iseg - 1) + i) = REAL(i, dp)/REAL(MAX(1, i - 1), dp)
            END DO
         END DO
         !   scal=zeros(staging_env%j,Float64)
         !   divide(arange(2,staging_env%j+1,typecode=Float64),
         !          arange(1,staging_env%j,typecode=Float64),scal[1:])
         !   scal[0]=1.
         !   scal=outerproduct(ones(staging_env%nseg),scal)

         IF (PRESENT(mass_beads)) THEN
            DO iat = 1, SIZE(mass)
               DO ib = 1, staging_env%p
                  mass_beads(ib, iat) = scal(ib)*mass(iat)
               END DO
            END DO
         END IF
         IF (PRESENT(mass_fict)) THEN
            DO iat = 1, SIZE(mass)
               DO ib = 1, staging_env%p
                  mass_fict(ib, iat) = scal(ib)*mass(iat)
               END DO
            END DO
         END IF
      END IF
   END SUBROUTINE staging_init_masses

! ***************************************************************************
!> \brief Transforms from the x into the u variables using a staging
!>      transformation for the positions
!> \param staging_env the environment for the staging transformation
!> \param ux will contain the u variable
!> \param x the positions to transform
!> \author fawzi
! **************************************************************************************************
   PURE SUBROUTINE staging_x2u(staging_env, ux, x)
      TYPE(staging_env_type), INTENT(IN)                 :: staging_env
      REAL(kind=dp), DIMENSION(:, :), INTENT(out)        :: ux
      REAL(kind=dp), DIMENSION(:, :), INTENT(in)         :: x

      INTEGER                                            :: k, s

      ux = x
      DO s = 0, staging_env%nseg - 1
         DO k = 2, staging_env%j
            ux(staging_env%j*s + k, :) = ux(staging_env%j*s + k, :) &
                                         - ((REAL(k - 1, dp)/REAL(k, dp) &
                                             *x(MODULO((staging_env%j*s + k + 1), staging_env%p), :) + &
                                             x(staging_env%j*s + 1, :)/REAL(k, dp)))
         END DO
      END DO
   END SUBROUTINE staging_x2u

! ***************************************************************************
!> \brief transform from the u variable to the x (back staging transformation
!>      for the positions)
!> \param staging_env the environment for the staging transformation
!> \param ux the u variable (positions to be backtransformed)
!> \param x will contain the positions
!> \author fawzi
! **************************************************************************************************
   PURE SUBROUTINE staging_u2x(staging_env, ux, x)
      TYPE(staging_env_type), INTENT(IN)                 :: staging_env
      REAL(kind=dp), DIMENSION(:, :), INTENT(in)         :: ux
      REAL(kind=dp), DIMENSION(:, :), INTENT(out)        :: x

      INTEGER                                            :: i, ist, j
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iii, jjj
      REAL(kind=dp)                                      :: const, const2

      j = staging_env%j
      const = REAL(j - 1, dp)/REAL(j, dp)
      const2 = 1._dp/REAL(j, dp)
      ALLOCATE (iii(staging_env%nseg), jjj(staging_env%nseg))
      DO i = 1, staging_env%nseg
         iii(i) = staging_env%j*(i - 1) + 1 !first el
      END DO
      DO i = 1, staging_env%nseg - 1
         jjj(i) = iii(i) + j ! next first el (pbc)
      END DO
      jjj(staging_env%nseg) = 1

      x = ux
      DO i = 1, staging_env%nseg
         x(j - 1 + iii(i), :) = x(j - 1 + iii(i), :) + &
                                const*ux(jjj(i), :) + ux(iii(i), :)*const2
      END DO
      DO ist = 1, staging_env%nseg
         DO i = staging_env%j - 2, 2, -1
            x(i + iii(ist), :) = x(i + iii(ist), :) + &
                                 REAL(i - 1, dp)/REAL(i, dp)*x(i + iii(ist) + 1, :) &
                                 + ux(iii(ist), :)/REAL(i, dp)
         END DO
      END DO
   END SUBROUTINE staging_u2x

! ***************************************************************************
!> \brief staging transformation for the forces
!> \param staging_env the environment for the staging transformation
!> \param uf will contain the forces after for the transformed variable
!> \param f the forces to transform
!> \author fawzi
! **************************************************************************************************
   PURE SUBROUTINE staging_f2uf(staging_env, uf, f)
      TYPE(staging_env_type), INTENT(IN)                 :: staging_env
      REAL(kind=dp), DIMENSION(:, :), INTENT(out)        :: uf
      REAL(kind=dp), DIMENSION(:, :), INTENT(in)         :: f

      INTEGER                                            :: i, idim, ij, ist, k
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iii, jjj, kkk
      REAL(kind=dp)                                      :: const, sum_f

      const = REAL(staging_env%j - 1, dp)/REAL(staging_env%j, dp)
      ALLOCATE (iii(staging_env%j), jjj(staging_env%j), &
                kkk(staging_env%j))
      DO ist = 1, staging_env%j - 1
         iii(ist) = (ist - 1)*staging_env%j + 1 ! first el
         jjj(ist) = iii(ist) + staging_env%j - 1 ! last el
         kkk(ist) = iii(ist) - 1 ! prev el
      END DO
      kkk(1) = staging_env%p

      uf = f
      ! staging beads
      DO k = 1, staging_env%nseg
         DO i = 2, staging_env%j
            uf(i + iii(k), :) = uf(i + iii(k), :) &
                                + REAL(i - 1, dp)/REAL(i, dp)*uf(i + iii(k) - 1, :)
         END DO
      END DO
      ! end point beads
      DO idim = 1, SIZE(uf, 2)
         DO k = 1, staging_env%nseg
            sum_f = 0._dp
            DO ij = 2, staging_env%j - 1
               sum_f = sum_f + uf((k - 1)*staging_env%j + ij, idim)
            END DO
            uf(iii(k), idim) = uf(iii(k), idim) + &
                               sum_f - const*(uf(jjj(k), idim) - uf(kkk(k), idim))
         END DO
      END DO
   END SUBROUTINE staging_f2uf

! ***************************************************************************
!> \brief calculates the harmonic force in the staging basis
!> \param staging_env the staging environment
!> \param mass_beads the masses of the beads
!> \param ux the positions of the beads in the staging basis
!> \param uf_h the harmonic forces (not accelerations)
!> \param e_h ...
!> \author fawzi
! **************************************************************************************************
   PURE SUBROUTINE staging_calc_uf_h(staging_env, mass_beads, ux, uf_h, e_h)
      TYPE(staging_env_type), INTENT(IN)                 :: staging_env
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: mass_beads, ux, uf_h
      REAL(KIND=dp), INTENT(OUT)                         :: e_h

      INTEGER                                            :: idim, isg, ist
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iii, jjj, kkk
      REAL(kind=dp)                                      :: d, f

      e_h = 0.0_dp

      ALLOCATE (iii(staging_env%nseg), jjj(staging_env%nseg), &
                kkk(staging_env%nseg))

      DO ist = 1, staging_env%nseg
         iii(ist) = (ist - 1)*staging_env%j + 1 ! first el
         jjj(ist) = iii(ist) + staging_env%j ! next first (pbc)
         kkk(ist) = iii(ist) - staging_env%j ! prev first el (pbc)
      END DO
      jjj(staging_env%nseg) = 1
      kkk(1) = staging_env%p - staging_env%j

      DO idim = 1, SIZE(mass_beads, 2)
         DO ist = 1, staging_env%nseg
            e_h = e_h + 0.5*mass_beads(1, idim)*staging_env%w_j**2* &
                  (ux(iii(ist), idim) - ux(jjj(ist), idim))**2
            uf_h(iii(ist), idim) = mass_beads(1, idim)*staging_env%w_j**2*( &
                                   2._dp*ux(iii(ist), idim) &
                                   - ux(jjj(ist), idim) &
                                   - ux(kkk(ist), idim) &
                                   )
            DO isg = 2, staging_env%j ! use 3 as start?
               d = ux((ist - 1)*staging_env%j + isg, idim)
               f = mass_beads((ist - 1)*staging_env%j + isg, idim)*staging_env%w_j**2*d
               e_h = e_h + 0.5_dp*f*d
               uf_h((ist - 1)*staging_env%j + isg, idim) = f
            END DO
         END DO
      END DO
   END SUBROUTINE staging_calc_uf_h

END MODULE pint_staging
