!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of charge equilibration in xTB
!> \author JGH
! **************************************************************************************************
MODULE xtb_eeq
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_control_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert,&
                                              cp_fm_matvec,&
                                              cp_fm_solve
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_element,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_unit_nr
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE ewald_pw_types,                  ONLY: ewald_pw_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_poisson_types,                ONLY: do_ewald_spme
   USE qs_dispersion_pairpot,           ONLY: dcnum_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE spme,                            ONLY: spme_forces,&
                                              spme_potential
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_eeq'

   PUBLIC :: xtb_eeq_calculation, xtb_eeq_forces

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param cnumbers ...
!> \param eeq_energy ...
!> \param lambda ...
! **************************************************************************************************
   SUBROUTINE xtb_eeq_calculation(qs_env, charges, cnumbers, eeq_energy, lambda)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: charges
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cnumbers
      REAL(KIND=dp), INTENT(INOUT)                       :: eeq_energy, lambda

      CHARACTER(len=*), PARAMETER :: routineN = 'xtb_eeq_calculation'

      INTEGER                                            :: ewald_type, handle, i, ia, iac, iar, &
                                                            iatom, ic, ikind, ir, iunit, jatom, &
                                                            jkind, natom, ncloc, ncvloc, nkind, &
                                                            nrloc, nrvloc, ns
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(3)                              :: periodic
      INTEGER, DIMENSION(:), POINTER                     :: cind, cvind, rind, rvind
      LOGICAL                                            :: defined, do_ewald, do_sparse
      REAL(KIND=dp)                                      :: ad, ala, alb, alk, alpha, betak, deth, &
                                                            dr, edis, gam, grc, grc1, grc2, kappa, &
                                                            papk, rdl, res, rmax, rzk, rzkm, rzkp, &
                                                            scn, totalcharge, xi, xr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: adiag, chia
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gab
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rj
      REAL(KIND=dp), DIMENSION(:), POINTER               :: av0, pv0, rv0, xv0, zv0
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: mat_struct, vec_struct
      TYPE(cp_fm_type)                                   :: eeq_mat, mmat, pmat, rhs_vec
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_tbe
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      iunit = cp_logger_get_default_unit_nr()

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      cell=cell, &
                      atprop=atprop, &
                      dft_control=dft_control)
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)

      xtb_control => dft_control%qs_control%xtb_control

      totalcharge = dft_control%charge

      IF (atprop%energy) THEN
         CALL atprop_array_init(atprop%atecoul, natom)
      END IF

      ! gamma[a,b]
      ALLOCATE (gab(nkind, nkind))
      gab = 0.0_dp
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined)
         IF (.NOT. defined) CYCLE
         CALL get_xtb_atom_param(xtb_atom_a, alpg=ala)
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, defined=defined)
            IF (.NOT. defined) CYCLE
            CALL get_xtb_atom_param(xtb_atom_b, alpg=alb)
            !
            gab(ikind, jkind) = SQRT(1._dp/(ala*ala + alb*alb))
            !
         END DO
      END DO

      ! Chi[a,a]
      ALLOCATE (chia(natom), adiag(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, xi=xi, kappa0=kappa, eta=gam)
         !
         scn = SQRT(cnumbers(iatom)) + 1.0e-14_dp
         chia(iatom) = xi - kappa*scn
         adiag(iatom) = gam + 2.0_dp*gab(ikind, ikind)*oorootpi
         !
      END DO

      do_ewald = xtb_control%do_ewald
      do_sparse = .FALSE.
      IF (natom > 100000 .AND. do_ewald) do_sparse = .TRUE.
      IF (do_ewald) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw)
         CALL get_cell(cell=cell, periodic=periodic, deth=deth)
         CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
         ad = 2.0_dp*alpha*oorootpi
         IF (ewald_type /= do_ewald_spme) THEN
            CALL cp_abort(__LOCATION__, "Only SPME Ewald method available with xTB.")
         END IF
         ! check alpha for MIC
         !
         rmax = MIN(8.0_dp, 0.5_dp*deth**0.333_dp)
         !
         CALL get_qs_env(qs_env=qs_env, sab_tbe=sab_tbe)
         !
         ns = natom
         IF (do_sparse) THEN
            CPABORT("sparse")
         ELSE
            CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
            CALL cp_fm_struct_create(mat_struct, context=blacs_env, para_env=para_env, &
                                     nrow_global=ns, ncol_global=ns)
            CALL cp_fm_create(eeq_mat, mat_struct)
            CALL cp_fm_set_all(eeq_mat, 1.0_dp, 0.0_dp)
            CALL cp_fm_create(pmat, mat_struct)
            CALL cp_fm_set_all(pmat, 1.0_dp, 0.0_dp)
            CALL cp_fm_create(mmat, mat_struct)
         END IF
         ! response matrix
         CALL neighbor_list_iterator_create(nl_iterator, sab_tbe)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                   iatom=iatom, jatom=jatom, r=rij)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
            CALL get_xtb_atom_param(xtb_atom_a, defined=defined)
            IF (.NOT. defined) CYCLE
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, defined=defined)
            IF (.NOT. defined) CYCLE
            dr = SQRT(SUM(rij**2))
            IF (iatom == jatom) THEN
               CPASSERT(dr < 1.e-6_dp)
               grc1 = adiag(iatom) - ad
               grc2 = adiag(iatom)
               IF (do_sparse) THEN
                  CPABORT("sparse")
               ELSE
                  CALL cp_fm_set_element(eeq_mat, iatom, jatom, grc1)
                  CALL cp_fm_set_element(pmat, iatom, jatom, grc2)
               END IF
            ELSE
               grc1 = erf(gab(ikind, jkind)*dr)/dr - erf(alpha*dr)/dr
               IF (dr < rmax) THEN
                  grc2 = erf(gab(ikind, jkind)*dr)/dr
               ELSE
                  grc2 = 0.0_dp
               END IF
               IF (do_sparse) THEN
                  CPABORT("sparse")
               ELSE
                  CALL cp_fm_set_element(eeq_mat, iatom, jatom, grc1)
                  CALL cp_fm_set_element(eeq_mat, jatom, iatom, grc1)
                  CALL cp_fm_set_element(pmat, iatom, jatom, grc2)
                  CALL cp_fm_set_element(pmat, jatom, iatom, grc2)
               END IF
            END IF
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)
         !
         ! preconditioner invers
         CALL cp_fm_invert(pmat, mmat)
         !
         ALLOCATE (av0(natom), xv0(natom), rv0(natom), pv0(natom), zv0(natom))
         ! Initial guess
         CALL cp_fm_matvec(mmat, -chia, xv0)
         grc = (SUM(xv0) - totalcharge)/REAL(natom, KIND=dp)
         xv0 = xv0 - grc
         CALL cp_fm_matvec(eeq_mat, xv0, rv0)
         CALL spme_potential(ewald_env, ewald_pw, cell, particle_set, xv0, particle_set, rv0)
         edis = 0.5_dp*DOT_PRODUCT(xv0, rv0) + DOT_PRODUCT(xv0, chia)
         rv0 = -(chia + rv0)
         grc = SUM(rv0)/REAL(natom, KIND=dp)
         rv0 = rv0 - grc
         res = SQRT(SUM(rv0*rv0))
         IF (iunit > 0) THEN
            WRITE (iunit, '(A)', advance='NO') " EEQ| "
         END IF
         ! PCG
         CALL cp_fm_matvec(mmat, rv0, zv0)
         grc = SUM(zv0)/REAL(natom, KIND=dp)
         zv0 = zv0 - grc
         pv0 = zv0
         DO i = 1, natom
            rzk = DOT_PRODUCT(rv0, zv0)
            CALL cp_fm_matvec(eeq_mat, pv0, av0)
            CALL spme_potential(ewald_env, ewald_pw, cell, particle_set, pv0, particle_set, av0)
            grc = SUM(av0)/REAL(natom, KIND=dp)
            av0 = av0 - grc
            papk = DOT_PRODUCT(pv0, av0)
            alk = rzk/papk
            xv0 = xv0 + alk*pv0
            rv0 = rv0 - alk*av0
            rzkm = DOT_PRODUCT(rv0, zv0)
            res = SQRT(DOT_PRODUCT(rv0, rv0))
            !
            CALL cp_fm_matvec(eeq_mat, xv0, zv0)
            CALL spme_potential(ewald_env, ewald_pw, cell, particle_set, xv0, particle_set, zv0)
            edis = 0.5_dp*DOT_PRODUCT(xv0, zv0) + DOT_PRODUCT(xv0, chia)
            IF (iunit > 0) THEN
               rdl = MAX(-LOG10(res), 0.0_dp)
               WRITE (iunit, '(I1)', advance='NO') MIN(NINT(rdl), 9)
               IF (MOD(i, 75) == 0) THEN
                  WRITE (iunit, *)
                  WRITE (iunit, '(A)', advance='NO') " EEQ| "
               END IF
            END IF
            IF (res < 1.e-7_dp) EXIT
            CALL cp_fm_matvec(mmat, rv0, zv0)
            grc = SUM(zv0)/REAL(natom, KIND=dp)
            zv0 = zv0 - grc
            rzkp = DOT_PRODUCT(rv0, zv0)
            betak = rzkp/rzk - rzkm/rzk
            pv0 = zv0 + betak*pv0
         END DO
         !
         charges = xv0
         CALL cp_fm_matvec(eeq_mat, xv0, rv0)
         CALL spme_potential(ewald_env, ewald_pw, cell, particle_set, xv0, particle_set, rv0)
         rv0 = -(chia + rv0)
         lambda = SUM(rv0)/REAL(natom, KIND=dp)
         rv0 = rv0 - lambda
         res = SQRT(SUM(rv0*rv0))
         ! energy:  -0.5*q^T.X - lambda*totalcharge
         eeq_energy = 0.5*SUM(charges(1:natom)*chia(1:natom)) - lambda*totalcharge
         IF (iunit > 0) WRITE (iunit, *)
         !
         DEALLOCATE (av0, xv0, rv0, pv0, zv0)
         !
      ELSE
         ns = natom + 1
         CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
         CALL cp_fm_struct_create(mat_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=ns, ncol_global=ns)
         CALL cp_fm_create(eeq_mat, mat_struct)
         CALL cp_fm_get_info(eeq_mat, nrow_local=nrloc, ncol_local=ncloc, &
                             row_indices=rind, col_indices=cind)
         CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=ns, ncol_global=1)
         CALL cp_fm_create(rhs_vec, vec_struct)
         CALL cp_fm_get_info(rhs_vec, nrow_local=nrvloc, ncol_local=ncvloc, &
                             row_indices=rvind, col_indices=cvind)
         !
         ! set up matrix
         CALL cp_fm_set_all(eeq_mat, 1.0_dp, 0.0_dp)
         CALL cp_fm_set_all(rhs_vec, 0.0_dp)
         DO ir = 1, nrloc
            iar = rind(ir)
            IF (iar > natom) CYCLE
            ri(1:3) = particle_set(iar)%r(1:3)
            DO ic = 1, ncloc
               iac = cind(ic)
               IF (iac > natom) CYCLE
               rj(1:3) = particle_set(iac)%r(1:3)
               IF (iar == iac) THEN
                  grc = adiag(iar)
               ELSE
                  ikind = kind_of(iar)
                  jkind = kind_of(iac)
                  rij(1:3) = ri(1:3) - rj(1:3)
                  rij = pbc(rij, cell)
                  dr = SQRT(SUM(rij**2))
                  grc = erf(gab(ikind, jkind)*dr)/dr
               END IF
               CALL cp_fm_set_element(eeq_mat, iar, iac, grc)
            END DO
         END DO
         ! set up rhs vector
         DO ir = 1, nrvloc
            iar = rvind(ir)
            DO ic = 1, ncvloc
               iac = cvind(ic)
               ia = MAX(iar, iac)
               IF (ia > natom) THEN
                  xr = totalcharge
               ELSE
                  xr = -chia(ia)
               END IF
               CALL cp_fm_set_element(rhs_vec, iar, iac, xr)
            END DO
         END DO
         !
         ! fm matrices are overwritten!
         CALL cp_fm_solve(eeq_mat, rhs_vec)
         !
         charges = 0.0_dp
         lambda = 0.0_dp
         DO ir = 1, nrvloc
            iar = rvind(ir)
            DO ic = 1, ncvloc
               iac = cvind(ic)
               ia = MAX(iar, iac)
               IF (ia <= natom) THEN
                  xr = rhs_vec%local_data(ir, ic)
                  charges(ia) = xr
               ELSE
                  lambda = rhs_vec%local_data(ir, ic)
               END IF
            END DO
         END DO
         CALL para_env%sum(lambda)
         CALL para_env%sum(charges)
         !
         ! energy:  -0.5*q^T.X - lambda*totalcharge
         eeq_energy = 0.5*SUM(charges(1:natom)*chia(1:natom)) - lambda*totalcharge
         !
      END IF

      IF (do_ewald) THEN
         IF (do_sparse) THEN
            CPABORT("NAN")
         ELSE
            CALL cp_fm_struct_release(mat_struct)
            CALL cp_fm_release(eeq_mat)
            CALL cp_fm_release(pmat)
            CALL cp_fm_release(mmat)
         END IF
      ELSE
         IF (do_sparse) THEN
            CPABORT("NAN")
         ELSE
            CALL cp_fm_struct_release(mat_struct)
            CALL cp_fm_struct_release(vec_struct)
            CALL cp_fm_release(eeq_mat)
            CALL cp_fm_release(rhs_vec)
         END IF
      END IF

      DEALLOCATE (gab, chia, adiag)

      CALL timestop(handle)

   END SUBROUTINE xtb_eeq_calculation

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param dcharges ...
!> \param cnumbers ...
!> \param dcnum ...
! **************************************************************************************************
   SUBROUTINE xtb_eeq_forces(qs_env, charges, dcharges, cnumbers, dcnum)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: charges, dcharges, cnumbers
      TYPE(dcnum_type), DIMENSION(:), INTENT(IN)         :: dcnum

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xtb_eeq_forces'

      INTEGER :: atom_a, atom_b, atom_c, ewald_type, handle, i, ia, iac, iar, iatom, ic, ikind, &
         ir, iunit, jatom, jkind, katom, kkind, natom, ncloc, ncvloc, nkind, nrloc, nrvloc, ns
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(3)                              :: periodic
      INTEGER, DIMENSION(:), POINTER                     :: cind, cvind, rind, rvind
      LOGICAL                                            :: defined, do_ewald, do_sparse, use_virial
      REAL(KIND=dp) :: ad, ala, alb, alk, alpha, betak, ctot, deth, dr, dr2, drk, gam, gam2, grc, &
         grc1, grc2, kappa, papk, qlam, qq, rdl, res, rmax, rzk, rzkm, rzkp, scn, totalcharge, xi, &
         xr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: adiag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: epforce, gab
      REAL(KIND=dp), DIMENSION(3)                        :: fdik, ri, rij, rik, rj
      REAL(KIND=dp), DIMENSION(:), POINTER               :: av0, chrgx, dchia, pv0, qlag, rv0, xv0, &
                                                            zv0
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: mat_struct, vec_struct
      TYPE(cp_fm_type)                                   :: eeq_mat, mmat, pmat, rhs_vec
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_tbe
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      iunit = cp_logger_get_default_unit_nr()

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      atprop=atprop, &
                      force=force, &
                      virial=virial, &
                      cell=cell, &
                      dft_control=dft_control)
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      xtb_control => dft_control%qs_control%xtb_control

      totalcharge = dft_control%charge

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of)

      ! gamma[a,b]
      ALLOCATE (gab(nkind, nkind))
      gab = 0.0_dp
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined)
         IF (.NOT. defined) CYCLE
         CALL get_xtb_atom_param(xtb_atom_a, alpg=ala)
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, defined=defined)
            IF (.NOT. defined) CYCLE
            CALL get_xtb_atom_param(xtb_atom_b, alpg=alb)
            !
            gab(ikind, jkind) = SQRT(1._dp/(ala*ala + alb*alb))
            !
         END DO
      END DO

      ! Chi[a,a]
      ALLOCATE (adiag(natom), qlag(natom))
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, eta=gam)
         !
         adiag(iatom) = gam + 2.0_dp*gab(ikind, ikind)*oorootpi
         !
      END DO

      do_ewald = xtb_control%do_ewald
      do_sparse = .FALSE.
      IF (natom > 100000 .AND. do_ewald) do_sparse = .TRUE.
      IF (do_ewald) THEN
         ! response matrix
         CALL get_qs_env(qs_env=qs_env, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw)
         CALL get_cell(cell=cell, periodic=periodic, deth=deth)
         CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
         ad = 2.0_dp*alpha*oorootpi
         IF (ewald_type /= do_ewald_spme) THEN
            CALL cp_abort(__LOCATION__, "Only SPME Ewald method available with xTB.")
         END IF
         ! check alpha for MIC
         !
         rmax = MIN(8.0_dp, 0.5_dp*deth**0.333_dp)
         !
         CALL get_qs_env(qs_env=qs_env, sab_tbe=sab_tbe)
         !
         ns = natom
         IF (do_sparse) THEN
            CPABORT("sparse")
         ELSE
            CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
            CALL cp_fm_struct_create(mat_struct, context=blacs_env, para_env=para_env, &
                                     nrow_global=ns, ncol_global=ns)
            CALL cp_fm_create(eeq_mat, mat_struct)
            CALL cp_fm_set_all(eeq_mat, 1.0_dp, 0.0_dp)
            CALL cp_fm_create(pmat, mat_struct)
            CALL cp_fm_set_all(pmat, 1.0_dp, 0.0_dp)
            CALL cp_fm_create(mmat, mat_struct)
            CALL cp_fm_set_all(mmat, 0.0_dp, 0.0_dp)
         END IF
         ! response matrix
         CALL neighbor_list_iterator_create(nl_iterator, sab_tbe)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                   iatom=iatom, jatom=jatom, r=rij)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
            CALL get_xtb_atom_param(xtb_atom_a, defined=defined)
            IF (.NOT. defined) CYCLE
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, defined=defined)
            IF (.NOT. defined) CYCLE
            dr = SQRT(SUM(rij**2))
            IF (iatom == jatom) THEN
               CPASSERT(dr < 1.e-6_dp)
               grc1 = adiag(iatom) - ad
               grc2 = adiag(iatom)
               IF (do_sparse) THEN
                  CPABORT("sparse")
               ELSE
                  CALL cp_fm_set_element(eeq_mat, iatom, iatom, grc1)
                  CALL cp_fm_set_element(pmat, iatom, iatom, grc2)
               END IF
            ELSE
               grc1 = erf(gab(ikind, jkind)*dr)/dr - erf(alpha*dr)/dr
               IF (dr < rmax) THEN
                  grc2 = erf(gab(ikind, jkind)*dr)/dr
               ELSE
                  grc2 = 0.0_dp
               END IF
               IF (do_sparse) THEN
                  CPABORT("sparse")
               ELSE
                  CALL cp_fm_set_element(eeq_mat, iatom, jatom, grc1)
                  CALL cp_fm_set_element(eeq_mat, jatom, iatom, grc1)
                  CALL cp_fm_set_element(pmat, iatom, jatom, grc2)
                  CALL cp_fm_set_element(pmat, jatom, iatom, grc2)
               END IF
            END IF
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)
         !
         ! preconditioner invers
         CALL cp_fm_invert(pmat, mmat)
         !
         ALLOCATE (av0(natom), xv0(natom), rv0(natom), pv0(natom), zv0(natom))
         ! Initial guess
         CALL cp_fm_matvec(mmat, -dcharges, xv0)
         grc = SUM(xv0)/REAL(natom, KIND=dp)
         xv0 = xv0 - grc
         CALL cp_fm_matvec(eeq_mat, xv0, rv0)
         CALL spme_potential(ewald_env, ewald_pw, cell, particle_set, xv0, particle_set, rv0)
         rv0 = -(dcharges + rv0)
         grc = SUM(rv0)/REAL(natom, KIND=dp)
         rv0 = rv0 - grc
         IF (iunit > 0) THEN
            WRITE (iunit, '(A)', advance='NO') " EEQ| "
         END IF
         ! PCG
         CALL cp_fm_matvec(mmat, rv0, zv0)
         grc = SUM(zv0)/REAL(natom, KIND=dp)
         zv0 = zv0 - grc
         pv0 = zv0
         DO i = 1, natom
            rzk = DOT_PRODUCT(rv0, zv0)
            CALL cp_fm_matvec(eeq_mat, pv0, av0)
            CALL spme_potential(ewald_env, ewald_pw, cell, particle_set, pv0, particle_set, av0)
            grc = SUM(av0)/REAL(natom, KIND=dp)
            av0 = av0 - grc
            papk = DOT_PRODUCT(pv0, av0)
            alk = rzk/papk
            xv0 = xv0 + alk*pv0
            rv0 = rv0 - alk*av0
            rzkm = DOT_PRODUCT(rv0, zv0)
            res = SQRT(DOT_PRODUCT(rv0, rv0))
            !
            CALL cp_fm_matvec(eeq_mat, xv0, zv0)
            CALL spme_potential(ewald_env, ewald_pw, cell, particle_set, xv0, particle_set, zv0)
            IF (iunit > 0) THEN
               rdl = MAX(-LOG10(res), 0.0_dp)
               WRITE (iunit, '(I1)', advance='NO') MIN(NINT(rdl), 9)
               IF (MOD(i, 75) == 0) THEN
                  WRITE (iunit, *)
                  WRITE (iunit, '(A)', advance='NO') " EEQ| "
               END IF
            END IF
            IF (res < 1.e-7_dp) EXIT
            CALL cp_fm_matvec(mmat, rv0, zv0)
            grc = SUM(zv0)/REAL(natom, KIND=dp)
            zv0 = zv0 - grc
            rzkp = DOT_PRODUCT(rv0, zv0)
            betak = rzkp/rzk - rzkm/rzk
            pv0 = zv0 + betak*pv0
         END DO
         qlag = -xv0
         CALL cp_fm_matvec(eeq_mat, xv0, rv0)
         CALL spme_potential(ewald_env, ewald_pw, cell, particle_set, xv0, particle_set, rv0)
         rv0 = -(dcharges + rv0)
         qlam = -SUM(rv0)/REAL(natom, KIND=dp)
         IF (iunit > 0) WRITE (iunit, *)
         !
         DEALLOCATE (av0, xv0, rv0, pv0, zv0)
         !
         CALL cp_fm_struct_release(mat_struct)
         CALL cp_fm_release(eeq_mat)
         CALL cp_fm_release(pmat)
         CALL cp_fm_release(mmat)

      ELSE
         ns = natom + 1
         CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
         CALL cp_fm_struct_create(mat_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=ns, ncol_global=ns)
         CALL cp_fm_create(eeq_mat, mat_struct)
         CALL cp_fm_get_info(eeq_mat, nrow_local=nrloc, ncol_local=ncloc, &
                             row_indices=rind, col_indices=cind)
         !
         CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=ns, ncol_global=1)
         CALL cp_fm_create(rhs_vec, vec_struct)
         CALL cp_fm_get_info(rhs_vec, nrow_local=nrvloc, ncol_local=ncvloc, &
                             row_indices=rvind, col_indices=cvind)
         ! set up matrix
         CALL cp_fm_set_all(eeq_mat, 1.0_dp, 0.0_dp)
         CALL cp_fm_set_all(rhs_vec, 0.0_dp)
         DO ir = 1, nrloc
            iar = rind(ir)
            IF (iar > natom) CYCLE
            ri(1:3) = particle_set(iar)%r(1:3)
            DO ic = 1, ncloc
               iac = cind(ic)
               IF (iac > natom) CYCLE
               rj(1:3) = particle_set(iac)%r(1:3)
               IF (iar == iac) THEN
                  grc = adiag(iar)
               ELSE
                  ikind = kind_of(iar)
                  jkind = kind_of(iac)
                  rij(1:3) = ri(1:3) - rj(1:3)
                  rij = pbc(rij, cell)
                  dr = SQRT(SUM(rij**2))
                  grc = erf(gab(ikind, jkind)*dr)/dr
               END IF
               CALL cp_fm_set_element(eeq_mat, iar, iac, grc)
            END DO
         END DO
         ! set up rhs vector
         DO ir = 1, nrvloc
            iar = rvind(ir)
            DO ic = 1, ncvloc
               iac = cvind(ic)
               ia = MAX(iar, iac)
               IF (ia > natom) THEN
                  xr = 0.0_dp
               ELSE
                  xr = dcharges(ia)
               END IF
               CALL cp_fm_set_element(rhs_vec, iar, iac, xr)
            END DO
         END DO
         !
         ! fm matrices are overwritten!
         CALL cp_fm_solve(eeq_mat, rhs_vec)
         !
         qlag = 0.0_dp
         qlam = 0.0_dp
         DO ir = 1, nrvloc
            iar = rvind(ir)
            DO ic = 1, ncvloc
               iac = cvind(ic)
               ia = MAX(iar, iac)
               IF (ia <= natom) THEN
                  xr = rhs_vec%local_data(ir, ic)
                  qlag(iar) = xr
               ELSE
                  qlam = rhs_vec%local_data(ir, ic)
               END IF
            END DO
         END DO
         CALL para_env%sum(qlam)
         CALL para_env%sum(qlag)

         CALL cp_fm_struct_release(mat_struct)
         CALL cp_fm_struct_release(vec_struct)
         CALL cp_fm_release(eeq_mat)
         CALL cp_fm_release(rhs_vec)

      END IF

      ALLOCATE (chrgx(natom), dchia(natom))
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, xi=xi, kappa0=kappa, eta=gam)
         !
         scn = SQRT(cnumbers(iatom)) + 1.0e-14_dp
         ctot = 0.5_dp*(charges(iatom) - qlag(iatom))
         dchia(iatom) = -ctot*kappa/scn
      END DO

      ! Forces from q*X
      CALL get_qs_env(qs_env=qs_env, &
                      local_particles=local_particles)
      DO ikind = 1, nkind
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            atom_a = atom_of_kind(iatom)
            DO i = 1, dcnum(iatom)%neighbors
               katom = dcnum(iatom)%nlist(i)
               kkind = kind_of(katom)
               atom_c = atom_of_kind(katom)
               rik = dcnum(iatom)%rik(:, i)
               drk = SQRT(SUM(rik(:)**2))
               IF (drk > 1.e-3_dp) THEN
                  fdik(:) = dchia(iatom)*dcnum(iatom)%dvals(i)*rik(:)/drk
                  force(ikind)%rho_elec(:, atom_a) = force(ikind)%rho_elec(:, atom_a) - fdik(:)
                  force(kkind)%rho_elec(:, atom_c) = force(kkind)%rho_elec(:, atom_c) + fdik(:)
                  IF (use_virial) THEN
                     CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
                     IF (atprop%stress) THEN
                        CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fdik, rik)
                        CALL virial_pair_force(atprop%atstress(:, :, katom), -0.5_dp, fdik, rik)
                     END IF
                  END IF
               END IF
            END DO
         END DO
      END DO
      ! Forces from (0.5*q+l)*dA/dR*q
      DO ikind = 1, nkind
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            atom_a = atom_of_kind(iatom)
            ri(1:3) = particle_set(iatom)%r(1:3)
            DO jatom = 1, natom
               IF (iatom == jatom) CYCLE
               jkind = kind_of(jatom)
               atom_b = atom_of_kind(jatom)
               qq = (0.5_dp*charges(iatom) - qlag(iatom))*charges(jatom)
               rj(1:3) = particle_set(jatom)%r(1:3)
               rij(1:3) = ri(1:3) - rj(1:3)
               rij = pbc(rij, cell)
               dr2 = SUM(rij**2)
               dr = SQRT(dr2)
               gam = gab(ikind, jkind)
               gam2 = gam*gam
               grc = 2._dp*gam*EXP(-gam2*dr2)*oorootpi/dr - erf(gam*dr)/dr2
               IF (do_ewald) THEN
                  grc = grc - 2._dp*alpha*EXP(-alpha**2*dr2)*oorootpi/dr + erf(alpha*dr)/dr2
               END IF
               fdik(:) = qq*grc*rij(:)/dr
               force(ikind)%rho_elec(:, atom_a) = force(ikind)%rho_elec(:, atom_a) + fdik(:)
               force(jkind)%rho_elec(:, atom_b) = force(jkind)%rho_elec(:, atom_b) - fdik(:)
               IF (use_virial) THEN
                  CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
                  IF (atprop%stress) THEN
                     CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fdik, rik)
                     CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp, fdik, rik)
                  END IF
               END IF
            END DO
         END DO
      END DO
      ! Forces from Ewald potential: (q+l)*A*q
      IF (do_ewald) THEN
         ALLOCATE (epforce(3, natom))
         epforce = 0.0_dp
         dchia = -charges + qlag
         chrgx = charges
         CALL spme_forces(ewald_env, ewald_pw, cell, particle_set, chrgx, &
                          particle_set, dchia, epforce)
         dchia = charges
         chrgx = qlag
         CALL spme_forces(ewald_env, ewald_pw, cell, particle_set, chrgx, &
                          particle_set, dchia, epforce)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            i = atom_of_kind(iatom)
            force(ikind)%rho_elec(:, i) = force(ikind)%rho_elec(:, i) + epforce(:, iatom)
!debf       IF (use_virial) THEN
!              CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
!              IF (atprop%stress) THEN
!                 CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fdik, rik)
!                 CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp, fdik, rik)
!              END IF
!debf       END IF
         END DO
         DEALLOCATE (epforce)
      END IF

      DEALLOCATE (gab, chrgx, dchia, adiag, qlag)

      CALL timestop(handle)

   END SUBROUTINE xtb_eeq_forces

END MODULE xtb_eeq
