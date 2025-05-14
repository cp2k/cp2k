!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE hairy_probes

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: hairy_probes_type
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_type
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: nso
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hairy_probes'

   INTEGER, PARAMETER, PRIVATE          :: BISECT_MAX_ITER = 400
   PUBLIC :: probe_occupancy, probe_occupancy_kp, AO_boundaries

CONTAINS

!**************************************************************************************
!> \brief subroutine to calculate occupation number and 'Fermi' level using the
!> \brief HAIR PROBE approach; gamma point calculation.
!> \param occ occupation numbers
!> \param fermi fermi level
!> \param kTS entropic energy contribution
!> \param energies MOs eigenvalues
!> \param coeff MOs coefficient
!> \param maxocc maximum allowed occupation number of an MO (1 or 2)
!> \param probe hairy probe
!> \param N number of electrons
! **************************************************************************************************

   SUBROUTINE probe_occupancy(occ, fermi, kTS, energies, coeff, maxocc, probe, N)

      !i/o variables and arrays
      REAL(KIND=dp), INTENT(out)                         :: occ(:), fermi, kTS
      REAL(KIND=dp), INTENT(IN)                          :: energies(:)
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: coeff
      REAL(KIND=dp), INTENT(IN)                          :: maxocc
      TYPE(hairy_probes_type), INTENT(INOUT)             :: probe(:)
      REAL(KIND=dp), INTENT(IN)                          :: N

      REAL(KIND=dp), PARAMETER                           :: epsocc = 1.0e-12_dp

      INTEGER                                            :: iter, ncol_global, nrow_global
      REAL(KIND=dp)                                      :: de, delta_fermi, fermi_fit, fermi_half, &
                                                            fermi_max, fermi_min, h, N_fit, &
                                                            N_half, N_max, N_min, N_now, y0, y1, y2
      REAL(KIND=dp), ALLOCATABLE                         :: smatrix_squared(:, :)
      REAL(KIND=dp), POINTER                             :: smatrix(:, :)

!subroutine variables and arrays
!smatrix: Spherical MOs matrix
!squared Spherical MOs matrix

      CALL cp_fm_get_info(coeff, &
                          nrow_global=nrow_global, &
                          ncol_global=ncol_global)
      ALLOCATE (smatrix(nrow_global, ncol_global))
      CALL cp_fm_get_submatrix(coeff, smatrix)

      ALLOCATE (smatrix_squared(nrow_global, ncol_global))
      smatrix_squared(:, :) = smatrix(:, :)**2.0d0

!**************************************************************************************
!1)calculate Fermi energy using the hairy probe formula
!**************************************************************************************
      de = probe(1)%T*LOG((1.0_dp - epsocc)/epsocc)
      de = MAX(de, 0.5_dp)
      fermi_max = MAXVAL(probe%mu) + de
      fermi_min = MINVAL(probe%mu) - de

      CALL HP_occupancy(probe=probe, matrix=smatrix_squared, energies=energies, &
                        maxocc=maxocc, fermi=fermi_max, occupancy=occ, kTS=kTS)
      N_max = accurate_sum(occ)

      CALL HP_occupancy(probe=probe, matrix=smatrix_squared, energies=energies, &
                        maxocc=maxocc, fermi=fermi_min, occupancy=occ, kTS=kTS)
      N_min = accurate_sum(occ)

      iter = 0
      DO WHILE (ABS(N_max - N_min) > N*epsocc)
         iter = iter + 1

         fermi_half = (fermi_max + fermi_min)/2.0_dp
         CALL HP_occupancy(probe=probe, matrix=smatrix_squared, energies=energies, &
                           maxocc=maxocc, fermi=fermi_half, occupancy=occ, kTS=kTS)
         N_half = accurate_sum(occ)

         h = fermi_half - fermi_min
         IF (h .GT. N*epsocc*100) THEN
            y0 = N_min - N
            y1 = N_half - N
            y2 = N_max - N

            CALL three_point_zero(y0, y1, y2, h, delta_fermi)
            fermi_fit = fermi_min + delta_fermi

            CALL HP_occupancy(probe=probe, matrix=smatrix_squared, energies=energies, &
                              maxocc=maxocc, fermi=fermi_fit, occupancy=occ, kTS=kTS)
            N_fit = accurate_sum(occ)
         END IF

         !define 1st bracked using fermi_half
         IF (N_half < N) THEN
            fermi_min = fermi_half
            N_min = N_half
         ELSE IF (N_half > N) THEN
            fermi_max = fermi_half
            N_max = N_half
         ELSE
            fermi_min = fermi_half
            N_min = N_half
            fermi_max = fermi_half
            N_max = N_half
            h = 0.0d0
         END IF

         !define 2nd bracker using fermi_fit
         IF (h .GT. N*epsocc*100) THEN
            IF (fermi_fit .GE. fermi_min .AND. fermi_fit .LE. fermi_max) THEN
               IF (N_fit < N) THEN
                  fermi_min = fermi_fit
                  N_min = N_fit
               ELSE IF (N_fit > N) THEN
                  fermi_max = fermi_fit
                  N_max = N_fit
               END IF
            END IF
         END IF

         IF (ABS(N_max - N) < N*epsocc) THEN
            fermi = fermi_max
            EXIT
         ELSE IF (ABS(N_min - N) < N*epsocc) THEN
            fermi = fermi_min
            EXIT
         END IF

         IF (iter > BISECT_MAX_ITER) THEN
            CPWARN("Maximum number of iterations reached while finding the Fermi energy")
            EXIT
         END IF
      END DO

!**************************************************************************************
!2)calculate occupation numbers according to hairy probe formula
!**************************************************************************************
      occ(:) = 0.0_dp
      N_now = 0.0_dp
      CALL HP_occupancy(probe=probe, matrix=smatrix_squared, energies=energies, &
                        maxocc=maxocc, fermi=fermi, occupancy=occ, kTS=kTS)
      N_now = accurate_sum(occ)

      IF (ABS(N_now - N) > N*epsocc) CPWARN("Total number of electrons is not accurate - HP")

      DEALLOCATE (smatrix, smatrix_squared)

   END SUBROUTINE probe_occupancy

!**************************************************************************************
!> \brief subroutine to calculate occupation number and 'Fermi' level using the
!> \brief HAIR PROBE approach; kpoints calculation.
!> \param occ occupation numbers
!> \param fermi fermi level
!> \param kTS entropic energy contribution
!> \param energies eigenvalues
!> \param rcoeff ...
!> \param icoeff ...
!> \param maxocc maximum allowed occupation number of an MO (1 or 2)
!> \param probe hairy probe
!> \param N number of electrons
!> \param wk weight of kpoints
! **************************************************************************************************
   SUBROUTINE probe_occupancy_kp(occ, fermi, kTS, energies, rcoeff, icoeff, maxocc, probe, N, wk)

      REAL(KIND=dp), INTENT(OUT)                         :: occ(:, :, :), fermi, kTS
      REAL(KIND=dp), INTENT(IN)                          :: energies(:, :, :), rcoeff(:, :, :, :), &
                                                            icoeff(:, :, :, :), maxocc
      TYPE(hairy_probes_type), INTENT(IN)                :: probe(:)
      REAL(KIND=dp), INTENT(IN)                          :: N, wk(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'probe_occupancy_kp'
      REAL(KIND=dp), PARAMETER                           :: epsocc = 1.0e-12_dp

      INTEGER                                            :: handle, ikp, ispin, iter, nAO, nkp, nMO, &
                                                            nspin
      REAL(KIND=dp)                                      :: de, delta_fermi, fermi_fit, fermi_half, &
                                                            fermi_max, fermi_min, h, kTS_kp, &
                                                            N_fit, N_half, N_max, N_min, N_now, &
                                                            y0, y1, y2
      REAL(KIND=dp), ALLOCATABLE                         :: coeff_squared(:, :, :, :)

      CALL timeset(routineN, handle)

!**************************************************************************************
!1)calculate Fermi energy using the hairy probe formula
!**************************************************************************************
      nAO = SIZE(rcoeff, 1)
      nMO = SIZE(rcoeff, 2)
      nkp = SIZE(rcoeff, 3)
      nspin = SIZE(rcoeff, 4)

      ALLOCATE (coeff_squared(nAO, nMO, nkp, nspin))
      coeff_squared(:, :, :, :) = rcoeff(:, :, :, :)**2.0d0 + icoeff(:, :, :, :)**2.0d0

      occ(:, :, :) = 0.0_dp

      !define initial brackets
      de = probe(1)%T*LOG((1.0_dp - epsocc)/epsocc)
      de = MAX(de, 0.5_dp)
      fermi_max = MAXVAL(probe%mu) + de
      fermi_min = MINVAL(probe%mu) - de

      N_max = 0.0_dp
      kTS = 0.0_dp
      !***HP loop
      DO ispin = 1, nspin
         DO ikp = 1, nkp
            CALL HP_occupancy(probe, coeff_squared(:, :, ikp, ispin), energies(:, ikp, ispin), &
                              maxocc, fermi_max, occ(:, ikp, ispin), kTS_kp)

            kTS = kTS + kTS_kp*wk(ikp) !entropic contribution
            N_max = N_max + accurate_sum(occ(1:nMo, ikp, ispin))*wk(ikp)
         END DO
      END DO
      !***HP loop

      N_min = 0.0_dp
      kTS = 0.0_dp
      !***HP loop
      DO ispin = 1, nspin
         DO ikp = 1, nkp
            CALL HP_occupancy(probe, coeff_squared(:, :, ikp, ispin), energies(:, ikp, ispin), &
                              maxocc, fermi_min, occ(:, ikp, ispin), kTS_kp)

            kTS = kTS + kTS_kp*wk(ikp) !entropic contribution
            N_min = N_min + accurate_sum(occ(1:nMo, ikp, ispin))*wk(ikp)
         END DO
      END DO
      !***HP loop

      iter = 0
      DO WHILE (ABS(N_max - N_min) > N*epsocc)
         iter = iter + 1
         fermi_half = (fermi_max + fermi_min)/2.0_dp
         N_half = 0.0_dp
         kTS = 0.0_dp
         !***HP loop
         DO ispin = 1, nspin
            DO ikp = 1, nkp
               CALL HP_occupancy(probe, coeff_squared(:, :, ikp, ispin), energies(:, ikp, ispin), &
                                 maxocc, fermi_half, occ(:, ikp, ispin), kTS_kp)

               kTS = kTS + kTS_kp*wk(ikp) !entropic contribution
               N_half = N_half + accurate_sum(occ(1:nMo, ikp, ispin))*wk(ikp)
            END DO
         END DO
         !***HP loop

         h = fermi_half - fermi_min
         IF (h .GT. N*epsocc*100) THEN
            y0 = N_min - N
            y1 = N_half - N
            y2 = N_max - N

            CALL three_point_zero(y0, y1, y2, h, delta_fermi)
            fermi_fit = fermi_min + delta_fermi
            N_fit = 0.0_dp
            kTS = 0.0_dp

            !***HP loop
            DO ispin = 1, nspin
               DO ikp = 1, nkp
                  CALL HP_occupancy(probe, coeff_squared(:, :, ikp, ispin), energies(:, ikp, ispin), &
                                    maxocc, fermi_fit, occ(:, ikp, ispin), kTS_kp)

                  kTS = kTS + kTS_kp*wk(ikp) !entropic contribution
                  N_fit = N_fit + accurate_sum(occ(1:nMo, ikp, ispin))*wk(ikp)
               END DO
            END DO
            !***HP loop

         END IF

         !define 1st bracked using fermi_half
         IF (N_half < N) THEN
            fermi_min = fermi_half
            N_min = N_half
         ELSE IF (N_half > N) THEN
            fermi_max = fermi_half
            N_max = N_half
         ELSE
            fermi_min = fermi_half
            N_min = N_half
            fermi_max = fermi_half
            N_max = N_half
            h = 0.0d0
         END IF

         !define 2nd bracker using fermi_fit
         IF (h .GT. N*epsocc*100) THEN
            IF (fermi_fit .GE. fermi_min .AND. fermi_fit .LE. fermi_max) THEN
               IF (N_fit < N) THEN
                  fermi_min = fermi_fit
                  N_min = N_fit
               ELSE IF (N_fit > N) THEN
                  fermi_max = fermi_fit
                  N_max = N_fit
               END IF
            END IF
         END IF

         IF (ABS(N_max - N) < N*epsocc) THEN
            fermi = fermi_max
            EXIT
         ELSE IF (ABS(N_min - N) < N*epsocc) THEN
            fermi = fermi_min
            EXIT
         END IF

         IF (iter > BISECT_MAX_ITER) THEN
            CPWARN("Maximum number of iterations reached while finding the Fermi energy")
            EXIT
         END IF
      END DO

!**************************************************************************************
!3)calculate occupation numbers using the hairy probe formula
!**************************************************************************************
      N_now = 0.0_dp
      kTS = 0.0_dp

      !***HP loop
      DO ispin = 1, nspin
         DO ikp = 1, nkp
            CALL HP_occupancy(probe, coeff_squared(:, :, ikp, ispin), energies(:, ikp, ispin), &
                              maxocc, fermi, occ(:, ikp, ispin), kTS_kp)

            kTS = kTS + kTS_kp*wk(ikp) !entropic contribution
            N_now = N_now + accurate_sum(occ(1:nMo, ikp, ispin))*wk(ikp)
         END DO
      END DO
      !***HP loop

      DEALLOCATE (coeff_squared)

      CALL timestop(handle)
   END SUBROUTINE probe_occupancy_kp

!**************************************************************************************
!
!
!
!**************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param probe ...
!> \param matrix ...
!> \param energies ...
!> \param maxocc ...
!> \param fermi ...
!> \param occupancy ...
!> \param kTS ...
! **************************************************************************************************
   SUBROUTINE HP_occupancy(probe, matrix, energies, maxocc, fermi, occupancy, kTS)

      TYPE(hairy_probes_type), INTENT(IN)                :: probe(:)
      REAL(KIND=dp), INTENT(IN)                          :: matrix(:, :), energies(:), maxocc, fermi
      REAL(KIND=dp), INTENT(OUT)                         :: occupancy(:), kTS

      INTEGER                                            :: iMO, ip, nMO, np
      REAL(KIND=dp)                                      :: alpha, C, f, fermi_fun, fermi_fun_sol, &
                                                            mu, s, sum_coeff, sum_coeff_sol, &
                                                            sum_fermi_fun, sum_fermi_fun_sol

!squared coefficient matrix

      nMO = SIZE(matrix, 2)
      np = SIZE(probe)
      kTS = 0.0_dp
      alpha = 1.0_dp
      ! kTS is the entropic contribution to the electronic energy

      mos2: DO iMO = 1, nMO

         sum_fermi_fun = 0.0_dp
         sum_coeff = 0.0_dp

         sum_fermi_fun_sol = 0.0_dp
         sum_coeff_sol = 0.0_dp
         fermi_fun_sol = 0.0_dp

         probes2: DO ip = 1, np
            IF (probe(ip)%alpha .LT. 1.0_dp) THEN
               alpha = probe(ip)%alpha
               !fermi distribution, solution probes
               CALL fermi_distribution(energies(iMO), fermi, probe(ip)%T, fermi_fun_sol)
               !sum of coefficients, solution probes
               sum_coeff_sol = sum_coeff_sol + SUM(matrix(probe(ip)%first_ao:probe(ip)%last_ao, iMO))
            ELSE
               C = SUM(matrix(probe(ip)%first_ao:probe(ip)%last_ao, iMO))
               !bias probes
               mu = fermi - probe(ip)%mu
               !fermi distribution, main probes
               CALL fermi_distribution(energies(iMO), mu, probe(ip)%T, fermi_fun)
               !sum fermi distribution * coefficients
               sum_fermi_fun = sum_fermi_fun + (fermi_fun*C)
               !sum cofficients
               sum_coeff = sum_coeff + C
            END IF
         END DO probes2

         sum_fermi_fun_sol = alpha*fermi_fun_sol*sum_coeff_sol
         sum_coeff_sol = alpha*sum_coeff_sol
         f = (sum_fermi_fun_sol + sum_fermi_fun)/(sum_coeff_sol + sum_coeff)
         occupancy(iMO) = f*maxocc

         !entropy kTS= kT*[f ln f + (1-f) ln (1-f)]
         IF (f .EQ. 0.0d0 .OR. f .EQ. 1.0d0) THEN
            s = 0.0d0
         ELSE
            s = f*LOG(f) + (1.0d0 - f)*LOG(1.0d0 - f)
         END IF
         kTS = kTS + probe(np)%T*maxocc*s
      END DO MOs2

   END SUBROUTINE HP_occupancy

!**************************************************************************************
!
!
!
!**************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param probe ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param nAO ...
! **************************************************************************************************
   SUBROUTINE AO_boundaries(probe, atomic_kind_set, qs_kind_set, particle_set, nAO)

      TYPE(hairy_probes_type), INTENT(INOUT)             :: probe
      TYPE(atomic_kind_type), INTENT(IN), POINTER        :: atomic_kind_set(:)
      TYPE(qs_kind_type), INTENT(IN), POINTER            :: qs_kind_set(:)
      TYPE(particle_type), INTENT(IN), POINTER           :: particle_set(:)
      INTEGER, INTENT(IN)                                :: nAO

      INTEGER                                            :: iatom, ii, ikind, iset, isgf, ishell, &
                                                            lshell, natom, nset, nsgf, p_atom
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

!from get_atomic_kind_set
!from get_atomic_kind
!from get_qs_kind_set
!subroutine variables and arrays
!from get_qs_kind
!from get_gto_basis_set
!from get_gto_basis_set
!from get_gto_basis_set

      !get sets from different modules:
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, nsgf=nsgf) !indexes for orbital symbols

      !get iAO boundaries:
      p_atom = SIZE(probe%atom_ids)
      probe%first_ao = nsgf !nAO
      probe%last_ao = 1
      isgf = 0

      atoms: DO iatom = 1, natom
         NULLIFY (orb_basis_set)
         !1) iatom is used to find the correct atom kind and its index (ikind) is the particle set
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)

         !2) ikind is used to find the basis set associate to that atomic kind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)

         !3) orb_basis_set is used to get the gto basis set variables
         IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                   nset=nset, nshell=nshell, l=l) !? ,cgf_symbol=bcgf_symbol )
         END IF

         !4) get iAO boundaries
         sets: DO iset = 1, nset

            shells: DO ishell = 1, nshell(iset)
               lshell = l(ishell, iset)

               isgf = isgf + nso(lshell)

               boundaries: DO ii = 1, p_atom
                  IF (iatom .NE. probe%atom_ids(ii)) THEN
                     CYCLE boundaries
                  ELSE
                     !defines iAO boundaries***********
                     probe%first_ao = MIN(probe%first_ao, isgf)
                     probe%last_ao = MAX(probe%last_ao, isgf)
                  END IF
               END DO boundaries

            END DO shells
         END DO sets
      END DO atoms

      IF (isgf .NE. nAO) CPWARN("row count does not correspond to nAO, number of rows in mo_coeff")

   END SUBROUTINE AO_boundaries

!*******************************************************************************************************
!*******************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param E ...
!> \param pot ...
!> \param temp ...
!> \param f_p ...
! **************************************************************************************************
   SUBROUTINE fermi_distribution(E, pot, temp, f_p)

      REAL(kind=dp), INTENT(IN)                          :: E, pot, temp
      REAL(kind=dp), INTENT(OUT)                         :: f_p

      REAL(KIND=dp)                                      :: arg, exponential, exponential_plus_1, f, &
                                                            one_minus_f

      ! have the result of exp go to zero instead of overflowing
      IF (E > pot) THEN
         arg = -(E - pot)/temp
         ! exponential is smaller than 1
         exponential = EXP(arg)
         exponential_plus_1 = exponential + 1.0_dp

         one_minus_f = exponential/exponential_plus_1
         f = 1.0_dp/exponential_plus_1
         f_p = one_minus_f
      ELSE
         arg = (E - pot)/temp
         ! exponential is smaller than 1
         exponential = EXP(arg)
         exponential_plus_1 = exponential + 1.0_dp

         f = 1.0_dp/exponential_plus_1
         one_minus_f = exponential/exponential_plus_1
         f_p = f
      END IF

   END SUBROUTINE fermi_distribution

!*******************************************************************************************************
!*******************************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param y0 ...
!> \param y1 ...
!> \param y2 ...
!> \param h ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE three_point_zero(y0, y1, y2, h, delta)

      REAL(kind=dp), INTENT(IN)                          :: y0, y1, y2, h
      REAL(kind=dp), INTENT(OUT)                         :: delta

      REAL(kind=dp)                                      :: a, b, c, d, u0

      a = (y2 - 2.0_dp*y1 + y0)/(2.0_dp*h*h)
      b = (4.0_dp*y1 - 3.0_dp*y0 - y2)/(2.0_dp*h)
      c = y0

      IF (ABS(a) .LT. 1.0e-15_dp) THEN
         delta = 0.0_dp
         RETURN
      END IF

      d = SQRT(b*b - 4.0_dp*a*c)

      u0 = (-b + d)/(2.0_dp*a)

      IF (u0 .GE. 0.0_dp .AND. u0 .LE. 2.0_dp*h) THEN
         delta = u0
         !RETURN
      ELSE
         u0 = (-b - d)/(2.0_dp*a)
         IF (u0 .GE. 0.0_dp .AND. u0 .LE. 2.0_dp*h) THEN
            delta = u0
            !RETURN
         ELSE
            IF (y1 .LT. 0.0_dp) delta = 2.0_dp*h
            IF (y1 .GE. 0.0_dp) delta = 0.0_dp
         END IF
      END IF

   END SUBROUTINE three_point_zero

END MODULE hairy_probes
