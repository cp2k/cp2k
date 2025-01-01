!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of dispersion using pair potentials
!> \author Johann Pototschnig
! **************************************************************************************************
MODULE qs_dispersion_d4
   USE atomic_kind_types, ONLY: atomic_kind_type, &
                                get_atomic_kind, &
                                get_atomic_kind_set
   USE distribution_1d_types, ONLY: distribution_1d_type
   USE eeq_method, ONLY: eeq_charges, eeq_forces
   USE machine, ONLY: m_flush, &
                      m_walltime
   USE cell_types, ONLY: cell_type, &
                         plane_distance, &
                         pbc, &
                         get_cell
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_force_types, ONLY: qs_force_type
   USE qs_kind_types, ONLY: get_qs_kind, &
                            qs_kind_type, &
                            set_qs_kind
   USE qs_neighbor_list_types, ONLY: get_iterator_info, &
                                     neighbor_list_iterate, &
                                     neighbor_list_iterator_create, &
                                     neighbor_list_iterator_p_type, &
                                     neighbor_list_iterator_release, &
                                     neighbor_list_set_p_type
   USE virial_methods, ONLY: virial_pair_force
   USE virial_types, ONLY: virial_type
   USE kinds, ONLY: dp
   USE particle_types, ONLY: particle_type
   USE qs_dispersion_types, ONLY: qs_dispersion_type
   USE qs_dispersion_utils, ONLY: cellhash
   USE qs_dispersion_cnum, ONLY: cnumber_init, dcnum_type, cnumber_release
   USE message_passing, ONLY: mp_para_env_type

#if defined(__DFTD4)
!&<
   USE dftd4,                           ONLY: d4_model, &
                                              damping_param, &
                                              get_dispersion, &
                                              get_rational_damping, &
                                              new, &
                                              new_d4_model, &
                                              realspace_cutoff, &
                                              structure_type, &
                                              rational_damping_param, &
                                              get_coordination_number, &
                                              get_lattice_points
   USE dftd4_charge,                    ONLY: get_charges
!&>
#endif
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_d4'

   PUBLIC :: calculate_dispersion_d4_pairpot

! **************************************************************************************************

CONTAINS

#if defined(__DFTD4)
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param evdw ...
!> \param calculate_forces ...
!> \param iw ...
!> \param atomic_energy ...
! **************************************************************************************************
   SUBROUTINE calculate_dispersion_d4_pairpot(qs_env, dispersion_env, evdw, calculate_forces, iw, &
                                              atomic_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_dispersion_type), INTENT(IN), POINTER      :: dispersion_env
      REAL(KIND=dp), INTENT(INOUT)                       :: evdw
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: iw
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: atomic_energy

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dispersion_d4_pairpot'

      INTEGER                                            :: atoma, cnfun, enshift, handle, iatom, &
                                                            ikind, mref, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, atomtype, kind_of
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: debug, grad, use_virial
      LOGICAL, DIMENSION(3)                              :: lperiod
      REAL(KIND=dp)                                      :: ed2, ed3, ev1, ev2, ev3, ev4, pd2, pd3, &
                                                            ta, tb, tc, td, te, ts
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cn, cnd, dEdcn, dEdq, edcn, edq, enerd2, &
                                                            enerd3, energies, energies3, q, qd
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ga, gradient, gwdcn, gwdq, gwvec, tvec, &
                                                            xyz
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: gdeb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: sigma, stress
      REAL(KIND=dp), DIMENSION(3, 3, 4)                  :: sdeb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      CLASS(damping_param), ALLOCATABLE                  :: param
      TYPE(d4_model)                                     :: disp
      TYPE(structure_type)                               :: mol
      TYPE(realspace_cutoff)                             :: cutoff

      CALL timeset(routineN, handle)

      debug = dispersion_env%d4_debug

      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, atomic_kind_set=atomic_kind_set, &
                      cell=cell, force=force, virial=virial, para_env=para_env)
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)

      !get information about particles
      natom = SIZE(particle_set)
      ALLOCATE (xyz(3, natom), atomtype(natom))
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      DO iatom = 1, natom
         xyz(:, iatom) = particle_set(iatom)%r(:)
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), zatom=atomtype(iatom))
      END DO

      !get information about cell / lattice
      CALL get_cell(cell=cell, periodic=periodic)
      lperiod(1) = periodic(1) == 1
      lperiod(2) = periodic(2) == 1
      lperiod(3) = periodic(3) == 1
      ! enforce en shift method 1 (original/molecular)
      ! method 2 from paper on PBC seems not to work
      enshift = 1
      !IF (ALL(periodic == 0)) enshift = 1

      !prepare for the call to the dispersion function
      CALL new(mol, atomtype, xyz, lattice=cell%hmat, periodic=lperiod)
      CALL new_d4_model(disp, mol)

      IF (dispersion_env%ref_functional == "none") THEN
         CALL get_rational_damping("pbe", param, s9=0.0_dp)
         SELECT TYPE (param)
         TYPE is (rational_damping_param)
            param%s6 = dispersion_env%s6
            param%s8 = dispersion_env%s8
            param%a1 = dispersion_env%a1
            param%a2 = dispersion_env%a2
            param%alp = dispersion_env%alp
         END SELECT
      ELSE
         CALL get_rational_damping(dispersion_env%ref_functional, param, s9=dispersion_env%s9)
         SELECT TYPE (param)
         TYPE is (rational_damping_param)
            dispersion_env%s6 = param%s6
            dispersion_env%s8 = param%s8
            dispersion_env%a1 = param%a1
            dispersion_env%a2 = param%a2
            dispersion_env%alp = param%alp
         END SELECT
      END IF

      ! Coordination number cutoff
      cutoff%cn = dispersion_env%rc_cn
      ! Two-body interaction cutoff
      cutoff%disp2 = dispersion_env%rc_d4*2._dp
      ! Three-body interaction cutoff
      cutoff%disp3 = dispersion_env%rc_disp*2._dp
      IF (cutoff%disp3 > cutoff%disp2) THEN
         CPABORT("D4: Three-body cutoff should be smaller than two-body cutoff")
      END IF

      IF (calculate_forces) THEN
         grad = .TRUE.
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      ELSE
         grad = .FALSE.
         use_virial = .FALSE.
      END IF

      IF (dispersion_env%d4_reference_code) THEN

         !> Wrapper to handle the evaluation of dispersion energy and derivatives
         IF (.NOT. dispersion_env%doabc) THEN
            CPWARN("Using D4_REFERENCE_CODE enforces calculation of C9 term.")
         END IF
         IF (grad) THEN
            ALLOCATE (gradient(3, natom))
            CALL get_dispersion(mol, disp, param, cutoff, evdw, gradient, stress)
            IF (calculate_forces) THEN
               IF (use_virial) THEN
                  virial%pv_virial = virial%pv_virial - stress/para_env%num_pe
               END IF
               DO iatom = 1, natom
                  ikind = kind_of(iatom)
                  atoma = atom_of_kind(iatom)
                  force(ikind)%dispersion(:, atoma) = &
                     force(ikind)%dispersion(:, atoma) + gradient(:, iatom)/para_env%num_pe
               END DO
            END IF
            DEALLOCATE (gradient)
         ELSE
            CALL get_dispersion(mol, disp, param, cutoff, evdw)
         END IF
         !dispersion energy is computed by every MPI process
         evdw = evdw/para_env%num_pe
         IF (dispersion_env%ext_charges) dispersion_env%dcharges = 0.0_dp
         IF (PRESENT(atomic_energy)) THEN
            CPWARN("Atomic energies not available for D4 reference code")
            atomic_energy = 0.0_dp
         END IF

      ELSE

         IF (iw > 0) THEN
            WRITE (iw, '(/,T2,A)') '!-----------------------------------------------------------------------------!'
            WRITE (iw, FMT="(T32,A)") "DEBUG D4 DISPERSION"
            WRITE (iw, '(T2,A)') '!-----------------------------------------------------------------------------!'
            WRITE (iw, '(A,T71,A10)') " DEBUG D4| Reference functional   ", TRIM(dispersion_env%ref_functional)
            WRITE (iw, '(A,T71,F10.4)') " DEBUG D4| Scaling parameter (s6) ", dispersion_env%s6
            WRITE (iw, '(A,T71,F10.4)') " DEBUG D4| Scaling parameter (s8) ", dispersion_env%s8
            WRITE (iw, '(A,T71,F10.4)') " DEBUG D4| BJ Damping parameter (a1) ", dispersion_env%a1
            WRITE (iw, '(A,T71,F10.4)') " DEBUG D4| BJ Damping parameter (a2) ", dispersion_env%a2
            WRITE (iw, '(A,T71,E10.4)') " DEBUG D4| Cutoff value coordination numbers ", dispersion_env%eps_cn
            WRITE (iw, '(A,T71,F10.4)') " DEBUG D4| Cutoff radius coordination numbers ", dispersion_env%rc_cn
            WRITE (iw, '(A,T71,I10)') " DEBUG D4| Coordination number function type ", dispersion_env%cnfun
            WRITE (iw, '(A,T71,F10.4)') " DEBUG D4| Cutoff radius 2-body terms [bohr]", 2._dp*dispersion_env%rc_d4
            WRITE (iw, '(A,T71,F10.4)') " DEBUG D4| Cutoff radius 3-body terms [bohr]", 2._dp*dispersion_env%rc_disp
         END IF

         td = 0.0_dp
         IF (debug .AND. iw > 0) THEN
            ts = m_walltime()
            CALL refd4_debug(param, disp, mol, cutoff, grad, dispersion_env%doabc, &
                             enerd2, enerd3, cnd, qd, Edcn, Edq, gdeb, sdeb)
            te = m_walltime()
            td = te - ts
         END IF

         tc = 0.0_dp
         ts = m_walltime()

         mref = MAXVAL(disp%ref)
         ! Coordination numbers
         cnfun = dispersion_env%cnfun
         CALL cnumber_init(qs_env, cn, dcnum, cnfun, grad)
         IF (debug .AND. iw > 0) THEN
            WRITE (iw, '(A,T71,F10.6)') " DEBUG D4| CN differences (max)", MAXVAL(ABS(cn - cnd))
            WRITE (iw, '(A,T71,F10.6)') " DEBUG D4| CN differences (ave)", SUM(ABS(cn - cnd))/natom
         END IF

         ! EEQ charges
         ALLOCATE (q(natom))
         IF (dispersion_env%ext_charges) THEN
            q(1:natom) = dispersion_env%charges(1:natom)
         ELSE
            CALL eeq_charges(qs_env, q, dispersion_env%eeq_sparam, 2, enshift)
         END IF
         IF (debug .AND. iw > 0) THEN
            WRITE (iw, '(A,T71,F10.6)') " DEBUG D4| Charge differences (max)", MAXVAL(ABS(q - qd))
            WRITE (iw, '(A,T71,F10.6)') " DEBUG D4| Charge differences (ave)", SUM(ABS(q - qd))/natom
         END IF
         ! Weights for C6 calculation
         ALLOCATE (gwvec(mref, natom))
         IF (grad) ALLOCATE (gwdcn(mref, natom), gwdq(mref, natom))
         CALL disp%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)

         ALLOCATE (energies(natom))
         energies(:) = 0.0_dp
         IF (grad) THEN
            ALLOCATE (gradient(3, natom), ga(3, natom))
            ALLOCATE (dEdcn(natom), dEdq(natom))
            dEdcn(:) = 0.0_dp; dEdq(:) = 0.0_dp
            ga(:, :) = 0.0_dp
            sigma(:, :) = 0.0_dp
         END IF
         CALL dispersion_2b(dispersion_env, cutoff%disp2, disp%r4r2, &
                            gwvec, gwdcn, gwdq, disp%c6, disp%ref, &
                            energies, dEdcn, dEdq, grad, ga, sigma)
         IF (grad) THEN
            gradient(1:3, 1:natom) = ga(1:3, 1:natom)
            stress = sigma
            IF (debug) THEN
               CALL para_env%sum(ga)
               CALL para_env%sum(sigma)
               IF (iw > 0) THEN
                  CALL gerror(ga, gdeb(:, :, 1), ev1, ev2, ev3, ev4)
                  WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| RMS error Gradient [2B]", ev1, ev2, " %"
                  WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| MAV error Gradient [2B]", ev3, ev4, " %"
                  IF (use_virial) THEN
                     CALL serror(sigma, sdeb(:, :, 1), ev1, ev2)
                     WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| MAV error Stress [2B]", ev1, ev2, " %"
                  END IF
               END IF
            END IF
         END IF
         ! no contribution from dispersion_3b as q=0 (but q is changed!)
         ! so we callculate this here
         IF (grad) THEN
            IF (dispersion_env%ext_charges) THEN
               dispersion_env%dcharges = dEdq
            ELSE
               CALL para_env%sum(dEdq)
               ga(:, :) = 0.0_dp
               sigma = 0.0_dp
               CALL eeq_forces(qs_env, q, dEdq, ga, sigma, dispersion_env%eeq_sparam, &
                               2, enshift, response_only=.TRUE.)
               gradient(1:3, 1:natom) = gradient(1:3, 1:natom) + ga(1:3, 1:natom)
               stress = stress + sigma
               IF (debug) THEN
                  CALL para_env%sum(ga)
                  CALL para_env%sum(sigma)
                  IF (iw > 0) THEN
                     CALL verror(dEdq, Edq, ev1, ev2)
                     WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| MAV error Derivative dEdq", ev1, ev2, " %"
                     CALL gerror(ga, gdeb(:, :, 2), ev1, ev2, ev3, ev4)
                     WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| RMS error Gradient [dEdq]", ev1, ev2, " %"
                     WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| MAV error Gradient [dEdq]", ev3, ev4, " %"
                     IF (use_virial) THEN
                        CALL serror(sigma, sdeb(:, :, 2), ev1, ev2)
                        WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| MAV error Stress [dEdq]", ev1, ev2, " %"
                     END IF
                  END IF
               END IF
            END IF
         END IF

         IF (dispersion_env%doabc) THEN
            ALLOCATE (energies3(natom))
            energies3(:) = 0.0_dp
            q(:) = 0.0_dp
            ! i.e. dc6dq = dEdq = 0
            CALL disp%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)
            !
            IF (grad) THEN
               gwdq = 0.0_dp
               ga(:, :) = 0.0_dp
               sigma = 0.0_dp
            END IF
            CALL get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, tvec)
            CALL dispersion_3b(qs_env, dispersion_env, tvec, cutoff%disp3, disp%r4r2, &
                               gwvec, gwdcn, gwdq, disp%c6, disp%ref, &
                               energies3, dEdcn, dEdq, grad, ga, sigma)
            IF (grad) THEN
               gradient(1:3, 1:natom) = gradient(1:3, 1:natom) + ga(1:3, 1:natom)
               stress = stress + sigma
               IF (debug) THEN
                  CALL para_env%sum(ga)
                  CALL para_env%sum(sigma)
                  IF (iw > 0) THEN
                     CALL gerror(ga, gdeb(:, :, 3), ev1, ev2, ev3, ev4)
                     WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| RMS error Gradient [3B]", ev1, ev2, " %"
                     WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| MAV error Gradient [3B]", ev3, ev4, " %"
                     IF (use_virial) THEN
                        CALL serror(sigma, sdeb(:, :, 3), ev1, ev2)
                        WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| MAV error Stress [3B]", ev1, ev2, " %"
                     END IF
                  END IF
               END IF
            END IF
         END IF

         IF (grad) THEN
            CALL para_env%sum(dEdcn)
            ga(:, :) = 0.0_dp
            sigma = 0.0_dp
            CALL dEdcn_force(qs_env, dEdcn, dcnum, ga, sigma)
            gradient(1:3, 1:natom) = gradient(1:3, 1:natom) + ga(1:3, 1:natom)
            stress = stress + sigma
            IF (debug) THEN
               CALL para_env%sum(ga)
               CALL para_env%sum(sigma)
               IF (iw > 0) THEN
                  CALL verror(dEdcn, Edcn, ev1, ev2)
                  WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| MAV error Derivative dEdcn", ev1, ev2, " %"
                  CALL gerror(ga, gdeb(:, :, 4), ev1, ev2, ev3, ev4)
                  WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| RMS error Gradient [dEdcn]", ev1, ev2, " %"
                  WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| MAV error Gradient [dEdcn]", ev3, ev4, " %"
                  IF (use_virial) THEN
                     CALL serror(sigma, sdeb(:, :, 4), ev1, ev2)
                     WRITE (iw, '(A,T51,F14.10,T69,F10.4,A)') " DEBUG D4| MAV error Stress [dEdcn]", ev1, ev2, " %"
                  END IF
               END IF
            END IF
         END IF
         DEALLOCATE (q)
         CALL cnumber_release(cn, dcnum, grad)
         te = m_walltime()
         tc = tc + te - ts

         IF (debug) THEN
            ta = SUM(energies)
            CALL para_env%sum(ta)
            IF (iw > 0) THEN
               tb = SUM(enerd2)
               ed2 = ta - tb
               pd2 = ABS(ed2)/ABS(tb)*100.
               WRITE (iw, '(A,T51,F14.8,T69,F10.4,A)') " DEBUG D4| Energy error 2-body", ed2, pd2, " %"
            END IF
            IF (dispersion_env%doabc) THEN
               ta = SUM(energies3)
               CALL para_env%sum(ta)
               IF (iw > 0) THEN
                  tb = SUM(enerd3)
                  ed3 = ta - tb
                  pd3 = ABS(ed3)/ABS(tb)*100.
                  WRITE (iw, '(A,T51,F14.8,T69,F10.4,A)') " DEBUG D4| Energy error 3-body", ed3, pd3, " %"
               END IF
            END IF
            IF (iw > 0) THEN
               WRITE (iw, '(A,T67,F14.4)') " DEBUG D4| Time for reference code [s]", td
               WRITE (iw, '(A,T67,F14.4)') " DEBUG D4| Time for production code [s]", tc
            END IF
         END IF

         IF (dispersion_env%doabc) THEN
            energies(:) = energies(:) + energies3(:)
         END IF
         evdw = SUM(energies)
         IF (PRESENT(atomic_energy)) THEN
            atomic_energy(1:natom) = energies(1:natom)
         END IF

         IF (use_virial .AND. calculate_forces) THEN
            virial%pv_virial = virial%pv_virial - stress
         END IF
         IF (calculate_forces) THEN
            DO iatom = 1, natom
               ikind = kind_of(iatom)
               atoma = atom_of_kind(iatom)
               force(ikind)%dispersion(:, atoma) = &
                  force(ikind)%dispersion(:, atoma) + gradient(:, iatom)
            END DO
         END IF

         DEALLOCATE (energies)
         IF (dispersion_env%doabc) DEALLOCATE (energies3)
         IF (grad) THEN
            DEALLOCATE (gradient, ga)
         END IF

      END IF

      DEALLOCATE (xyz, atomtype)

      CALL timestop(handle)

   END SUBROUTINE calculate_dispersion_d4_pairpot

! **************************************************************************************************
!> \brief ...
!> \param param ...
!> \param disp ...
!> \param mol ...
!> \param cutoff ...
!> \param grad ...
!> \param doabc ...
!> \param enerd2 ...
!> \param enerd3 ...
!> \param cnd ...
!> \param qd ...
!> \param dEdcn ...
!> \param dEdq ...
!> \param gradient ...
!> \param stress ...
! **************************************************************************************************
   SUBROUTINE refd4_debug(param, disp, mol, cutoff, grad, doabc, &
                          enerd2, enerd3, cnd, qd, dEdcn, dEdq, gradient, stress)
      CLASS(damping_param)                               :: param
      TYPE(d4_model)                                     :: disp
      TYPE(structure_type)                               :: mol
      TYPE(realspace_cutoff)                             :: cutoff
      LOGICAL, INTENT(IN)                                :: grad, doabc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: enerd2, enerd3, cnd, qd, dEdcn, dEdq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: gradient
      REAL(KIND=dp), DIMENSION(3, 3, 4)                  :: stress

      INTEGER                                            :: mref, natom, i
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: q, qq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: lattr, gwdcn, gwdq, gwvec, &
                                                            c6, dc6dcn, dc6dq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: cndr, cndL, qdr, qdL

      mref = MAXVAL(disp%ref)
      natom = mol%nat

      ! Coordination numbers
      ALLOCATE (cnd(natom))
      IF (grad) ALLOCATE (cndr(3, natom, natom), cndL(3, 3, natom))
      CALL get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
      CALL get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, &
                                   cnd, cndr, cndL)
      ! EEQ charges
      ALLOCATE (qd(natom))
      IF (grad) ALLOCATE (qdr(3, natom, natom), qdL(3, 3, natom))
      CALL get_charges(mol, qd, qdr, qdL)
      ! C6 interpolation
      ALLOCATE (gwvec(mref, natom))
      IF (grad) ALLOCATE (gwdcn(mref, natom), gwdq(mref, natom))
      CALL disp%weight_references(mol, cnd, qd, gwvec, gwdcn, gwdq)
      ALLOCATE (c6(natom, natom))
      IF (grad) ALLOCATE (dc6dcn(natom, natom), dc6dq(natom, natom))
      CALL disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)
      CALL get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
      !
      IF (grad) THEN
         ALLOCATE (gradient(3, natom, 4))
         gradient = 0.0_dp
         stress = 0.0_dp
      END IF
      !
      ALLOCATE (enerd2(natom))
      enerd2(:) = 0.0_dp
      IF (grad) THEN
         ALLOCATE (dEdcn(natom), dEdq(natom))
         dEdcn(:) = 0.0_dp; dEdq(:) = 0.0_dp
      END IF
      CALL param%get_dispersion2(mol, lattr, cutoff%disp2, disp%r4r2, c6, dc6dcn, dc6dq, &
                                 enerd2, dEdcn, dEdq, gradient(:, :, 1), stress(:, :, 1))
      !
      IF (grad) THEN
         DO i = 1, 3
            gradient(i, :, 2) = MATMUL(qdr(i, :, :), dEdq(:))
            stress(i, :, 2) = MATMUL(qdL(i, :, :), dEdq(:))
         END DO
      END IF
      !
      IF (doabc) THEN
         ALLOCATE (q(natom), qq(natom))
         q(:) = 0.0_dp; qq(:) = 0.0_dp
         ALLOCATE (enerd3(natom))
         enerd3(:) = 0.0_dp
         CALL disp%weight_references(mol, cnd, q, gwvec, gwdcn, gwdq)
         CALL disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)
         CALL get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
         CALL param%get_dispersion3(mol, lattr, cutoff%disp3, disp%r4r2, c6, dc6dcn, dc6dq, &
                                    enerd3, dEdcn, qq, gradient(:, :, 3), stress(:, :, 3))
      END IF
      IF (grad) THEN
         DO i = 1, 3
            gradient(i, :, 4) = MATMUL(cndr(i, :, :), dEdcn(:))
            stress(i, :, 4) = MATMUL(cndL(i, :, :), dEdcn(:))
         END DO
      END IF

   END SUBROUTINE refd4_debug

#else

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param evdw ...
!> \param calculate_forces ...
!> \param iw ...
!> \param atomic_energy ...
! **************************************************************************************************
   SUBROUTINE calculate_dispersion_d4_pairpot(qs_env, dispersion_env, evdw, calculate_forces, &
                                              iw, atomic_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_dispersion_type), INTENT(IN), POINTER      :: dispersion_env
      REAL(KIND=dp), INTENT(INOUT)                       :: evdw
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: iw
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: atomic_energy

      MARK_USED(qs_env)
      MARK_USED(dispersion_env)
      MARK_USED(evdw)
      MARK_USED(calculate_forces)
      MARK_USED(iw)
      MARK_USED(atomic_energy)

      CPABORT("CP2K build without DFTD4")

   END SUBROUTINE calculate_dispersion_d4_pairpot

#endif

! **************************************************************************************************
!> \brief ...
!> \param dispersion_env ...
!> \param cutoff ...
!> \param r4r2 ...
!> \param gwvec ...
!> \param gwdcn ...
!> \param gwdq ...
!> \param c6ref ...
!> \param mrefs ...
!> \param energies ...
!> \param dEdcn ...
!> \param dEdq ...
!> \param calculate_forces ...
!> \param gradient ...
!> \param stress ...
! **************************************************************************************************
   SUBROUTINE dispersion_2b(dispersion_env, cutoff, r4r2, &
                            gwvec, gwdcn, gwdq, c6ref, mrefs, &
                            energies, dEdcn, dEdq, &
                            calculate_forces, gradient, stress)
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: r4r2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gwvec, gwdcn, gwdq
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: c6ref
      INTEGER, DIMENSION(:), INTENT(IN)                  :: mrefs
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: energies, dEdcn, dEdq
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: gradient, stress

      INTEGER                                            :: iatom, ikind, jatom, jkind, mepos, num_pe
      REAL(KINd=dp)                                      :: a1, a2, c6ij, cutoff2, d6, d8, dE, dr2, &
                                                            edisp, fac, gdisp, r0ij, rrij, s6, s8, &
                                                            t6, t8
      REAL(KINd=dp), DIMENSION(2)                        :: dcdcn, dcdq
      REAL(KINd=dp), DIMENSION(3)                        :: dG, rij
      REAL(KINd=dp), DIMENSION(3, 3)                     :: dS
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_vdw

      a1 = dispersion_env%a1
      a2 = dispersion_env%a2
      s6 = dispersion_env%s6
      s8 = dispersion_env%s8
      cutoff2 = cutoff*cutoff

      sab_vdw => dispersion_env%sab_vdw

      num_pe = 1
      CALL neighbor_list_iterator_create(nl_iterator, sab_vdw, nthread=num_pe)

      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij)
         ! vdW potential
         dr2 = SUM(rij(:)**2)
         IF (dr2 <= cutoff2 .AND. dr2 > 0.0000001_dp) THEN
            rrij = 3._dp*r4r2(ikind)*r4r2(jkind)
            r0ij = a1*SQRT(rrij) + a2
            IF (calculate_forces) THEN
               CALL get_c6derivs(c6ij, dcdcn, dcdq, iatom, jatom, ikind, jkind, &
                                 gwvec, gwdcn, gwdq, c6ref, mrefs)
            ELSE
               CALL get_c6value(c6ij, iatom, jatom, ikind, jkind, gwvec, c6ref, mrefs)
            END IF
            fac = 1._dp
            IF (iatom == jatom) fac = 0.5_dp
            t6 = 1.0_dp/(dr2**3 + r0ij**6)
            t8 = 1.0_dp/(dr2**4 + r0ij**8)

            edisp = (s6*t6 + s8*rrij*t8)*fac
            dE = -c6ij*edisp
            energies(iatom) = energies(iatom) + dE*0.5_dp
            energies(jatom) = energies(jatom) + dE*0.5_dp

            IF (calculate_forces) THEN
               d6 = -6.0_dp*dr2**2*t6**2
               d8 = -8.0_dp*dr2**3*t8**2
               gdisp = (s6*d6 + s8*rrij*d8)*fac
               dG(:) = -c6ij*gdisp*rij(:)
               gradient(:, iatom) = gradient(:, iatom) - dG
               gradient(:, jatom) = gradient(:, jatom) + dG
               dS(:, :) = SPREAD(dG, 1, 3)*SPREAD(rij, 2, 3)
               stress(:, :) = stress(:, :) + dS(:, :)
               dEdcn(iatom) = dEdcn(iatom) - dcdcn(1)*edisp
               dEdq(iatom) = dEdq(iatom) - dcdq(1)*edisp
               dEdcn(jatom) = dEdcn(jatom) - dcdcn(2)*edisp
               dEdq(jatom) = dEdq(jatom) - dcdq(2)*edisp
            END IF
         END IF
      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

   END SUBROUTINE dispersion_2b

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param tvec ...
!> \param cutoff ...
!> \param r4r2 ...
!> \param gwvec ...
!> \param gwdcn ...
!> \param gwdq ...
!> \param c6ref ...
!> \param mrefs ...
!> \param energies ...
!> \param dEdcn ...
!> \param dEdq ...
!> \param calculate_forces ...
!> \param gradient ...
!> \param stress ...
! **************************************************************************************************
   SUBROUTINE dispersion_3b(qs_env, dispersion_env, tvec, cutoff, r4r2, &
                            gwvec, gwdcn, gwdq, c6ref, mrefs, &
                            energies, dEdcn, dEdq, &
                            calculate_forces, gradient, stress)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: tvec
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: r4r2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gwvec, gwdcn, gwdq
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: c6ref
      INTEGER, DIMENSION(:), INTENT(IN)                  :: mrefs
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: energies, dEdcn, dEdq
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: gradient, stress

      INTEGER                                            :: iatom, ikind, jatom, jkind, katom, &
                                                            kkind, ktr, mepos, natom, num_pe
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(3)                              :: cell_b
      REAL(KINd=dp)                                      :: a1, a2, alp, ang, c6ij, c6ik, c6jk, c9, &
                                                            cutoff2, dang, dE, dfdmp, fac, fdmp, &
                                                            r0, r0ij, r0ik, r0jk, r1, r2, r2ij, &
                                                            r2ik, r2jk, r3, r5, rr, s6, s8, s9
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rcpbc
      REAL(KINd=dp), DIMENSION(2)                        :: dc6dcnij, dc6dcnik, dc6dcnjk, dc6dqij, &
                                                            dc6dqik, dc6dqjk
      REAL(KINd=dp), DIMENSION(3)                        :: dGij, dGik, dGjk, ra, rb, rb0, rij, vij, &
                                                            vik, vjk
      REAL(KINd=dp), DIMENSION(3, 3)                     :: dS
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_vdw
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL get_qs_env(qs_env=qs_env, natom=natom, cell=cell, &
                      atomic_kind_set=atomic_kind_set, particle_set=particle_set)

      ALLOCATE (rcpbc(3, natom))
      DO iatom = 1, natom
         rcpbc(:, iatom) = pbc(particle_set(iatom)%r(:), cell)
      END DO
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      a1 = dispersion_env%a1
      a2 = dispersion_env%a2
      s6 = dispersion_env%s6
      s8 = dispersion_env%s8
      s9 = dispersion_env%s9
      alp = dispersion_env%alp

      cutoff2 = cutoff**2

      sab_vdw => dispersion_env%sab_vdw

      num_pe = 1
      CALL neighbor_list_iterator_create(nl_iterator, sab_vdw, nthread=num_pe)

      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, jatom=jatom, r=rij)

         r2ij = SUM(rij(:)**2)
         IF (calculate_forces) THEN
            CALL get_c6derivs(c6ij, dc6dcnij, dc6dqij, iatom, jatom, ikind, jkind, &
                              gwvec, gwdcn, gwdq, c6ref, mrefs)
         ELSE
            CALL get_c6value(c6ij, iatom, jatom, ikind, jkind, gwvec, c6ref, mrefs)
         END IF
         r0ij = a1*SQRT(3._dp*r4r2(jkind)*r4r2(ikind)) + a2
         IF (r2ij <= cutoff2 .AND. r2ij > EPSILON(1._dp)) THEN
            CALL get_iterator_info(nl_iterator, cell=cell_b)
            rb0(:) = MATMUL(cell%hmat, cell_b)
            ra(:) = rcpbc(:, iatom)
            rb(:) = rcpbc(:, jatom) + rb0
            vij(:) = rb(:) - ra(:)

            DO katom = 1, MIN(iatom, jatom)
               kkind = kind_of(katom)
               IF (calculate_forces) THEN
                  CALL get_c6derivs(c6ik, dc6dcnik, dc6dqik, katom, iatom, kkind, ikind, &
                                    gwvec, gwdcn, gwdq, c6ref, mrefs)
                  CALL get_c6derivs(c6jk, dc6dcnjk, dc6dqjk, katom, jatom, kkind, jkind, &
                                    gwvec, gwdcn, gwdq, c6ref, mrefs)
               ELSE
                  CALL get_c6value(c6ik, katom, iatom, kkind, ikind, gwvec, c6ref, mrefs)
                  CALL get_c6value(c6jk, katom, jatom, kkind, jkind, gwvec, c6ref, mrefs)
               END IF
               c9 = -s9*SQRT(ABS(c6ij*c6ik*c6jk))
               r0ik = a1*SQRT(3._dp*r4r2(kkind)*r4r2(ikind)) + a2
               r0jk = a1*SQRT(3._dp*r4r2(kkind)*r4r2(jkind)) + a2
               r0 = r0ij*r0ik*r0jk
               fac = triple_scale(iatom, jatom, katom)
               DO ktr = 1, SIZE(tvec, 2)
                  vik(:) = rcpbc(:, katom) + tvec(:, ktr) - rcpbc(:, iatom)
                  r2ik = vik(1)*vik(1) + vik(2)*vik(2) + vik(3)*vik(3)
                  IF (r2ik > cutoff2 .OR. r2ik < EPSILON(1.0_dp)) CYCLE
                  vjk(:) = rcpbc(:, katom) + tvec(:, ktr) - rb(:)
                  r2jk = vjk(1)*vjk(1) + vjk(2)*vjk(2) + vjk(3)*vjk(3)
                  IF (r2jk > cutoff2 .OR. r2jk < EPSILON(1.0_dp)) CYCLE
                  r2 = r2ij*r2ik*r2jk
                  r1 = SQRT(r2)
                  r3 = r2*r1
                  r5 = r3*r2

                  fdmp = 1.0_dp/(1.0_dp + 6.0_dp*(r0/r1)**(alp/3.0_dp))
                  ang = 0.375_dp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)* &
                        (-r2ij + r2jk + r2ik)/r5 + 1.0_dp/r3

                  rr = ang*fdmp
                  dE = rr*c9*fac
                  energies(iatom) = energies(iatom) - dE/3._dp
                  energies(jatom) = energies(jatom) - dE/3._dp
                  energies(katom) = energies(katom) - dE/3._dp

                  IF (calculate_forces) THEN

                     dfdmp = -2.0_dp*alp*(r0/r1)**(alp/3.0_dp)*fdmp**2

                     ! d/drij
                     dang = -0.375_dp*(r2ij**3 + r2ij**2*(r2jk + r2ik) &
                                       + r2ij*(3.0_dp*r2jk**2 + 2.0_dp*r2jk*r2ik &
                                               + 3.0_dp*r2ik**2) &
                                       - 5.0_dp*(r2jk - r2ik)**2*(r2jk + r2ik))/r5
                     dGij(:) = c9*(-dang*fdmp + ang*dfdmp)/r2ij*vij

                     ! d/drik
                     dang = -0.375_dp*(r2ik**3 + r2ik**2*(r2jk + r2ij) &
                                       + r2ik*(3.0_dp*r2jk**2 + 2.0_dp*r2jk*r2ij &
                                               + 3.0_dp*r2ij**2) &
                                       - 5.0_dp*(r2jk - r2ij)**2*(r2jk + r2ij))/r5
                     dGik(:) = c9*(-dang*fdmp + ang*dfdmp)/r2ik*vik

                     ! d/drjk
                     dang = -0.375_dp*(r2jk**3 + r2jk**2*(r2ik + r2ij) &
                                       + r2jk*(3.0_dp*r2ik**2 + 2.0_dp*r2ik*r2ij &
                                               + 3.0_dp*r2ij**2) &
                                       - 5.0_dp*(r2ik - r2ij)**2*(r2ik + r2ij))/r5
                     dGjk(:) = c9*(-dang*fdmp + ang*dfdmp)/r2jk*vjk

                     gradient(:, iatom) = gradient(:, iatom) - dGij - dGik
                     gradient(:, jatom) = gradient(:, jatom) + dGij - dGjk
                     gradient(:, katom) = gradient(:, katom) + dGik + dGjk

                     dS(:, :) = SPREAD(dGij, 1, 3)*SPREAD(vij, 2, 3) &
                                + SPREAD(dGik, 1, 3)*SPREAD(vik, 2, 3) &
                                + SPREAD(dGjk, 1, 3)*SPREAD(vjk, 2, 3)

                     stress(:, :) = stress + dS*fac

                     dEdcn(iatom) = dEdcn(iatom) - dE*0.5_dp &
                                    *(dc6dcnij(1)/c6ij + dc6dcnik(2)/c6ik)
                     dEdcn(jatom) = dEdcn(jatom) - dE*0.5_dp &
                                    *(dc6dcnij(2)/c6ij + dc6dcnjk(2)/c6jk)
                     dEdcn(katom) = dEdcn(katom) - dE*0.5_dp &
                                    *(dc6dcnik(1)/c6ik + dc6dcnjk(1)/c6jk)

                     dEdq(iatom) = dEdq(iatom) - dE*0.5_dp &
                                   *(dc6dqij(1)/c6ij + dc6dqik(2)/c6ik)
                     dEdq(jatom) = dEdq(jatom) - dE*0.5_dp &
                                   *(dc6dqij(2)/c6ij + dc6dqjk(2)/c6jk)
                     dEdq(katom) = dEdq(katom) - dE*0.5_dp &
                                   *(dc6dqik(1)/c6ik + dc6dqjk(1)/c6jk)

                  END IF

               END DO
            END DO
         END IF
      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (rcpbc)

   END SUBROUTINE dispersion_3b

! **************************************************************************************************
!> \brief ...
!> \param ii ...
!> \param jj ...
!> \param kk ...
!> \return ...
! **************************************************************************************************
   FUNCTION triple_scale(ii, jj, kk) RESULT(triple)
      INTEGER, INTENT(IN)                                :: ii, jj, kk
      REAL(KIND=dp)                                      :: triple

      IF (ii == jj) THEN
         IF (ii == kk) THEN
            ! ii'i" -> 1/6
            triple = 1.0_dp/6.0_dp
         ELSE
            ! ii'j -> 1/2
            triple = 0.5_dp
         END IF
      ELSE
         IF (ii /= kk .AND. jj /= kk) THEN
            ! ijk -> 1 (full)
            triple = 1.0_dp
         ELSE
            ! ijj' and iji' -> 1/2
            triple = 0.5_dp
         END IF
      END IF

   END FUNCTION triple_scale

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dEdcn ...
!> \param dcnum ...
!> \param gradient ...
!> \param stress ...
! **************************************************************************************************
   SUBROUTINE dEdcn_force(qs_env, dEdcn, dcnum, gradient, stress)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: dEdcn
      TYPE(dcnum_type), DIMENSION(:), INTENT(IN)         :: dcnum
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: gradient
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: stress

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dEdcn_force'

      INTEGER                                            :: handle, i, ia, iatom, ikind, katom, &
                                                            natom, nkind
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: drk
      REAL(KIND=dp), DIMENSION(3)                        :: fdik, rik
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind, natom=natom, &
                      local_particles=local_particles, &
                      virial=virial)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      DO ikind = 1, nkind
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            DO i = 1, dcnum(iatom)%neighbors
               katom = dcnum(iatom)%nlist(i)
               rik = dcnum(iatom)%rik(:, i)
               drk = SQRT(SUM(rik(:)**2))
               fdik(:) = -(dEdcn(iatom) + dEdcn(katom))*dcnum(iatom)%dvals(i)*rik(:)/drk
               gradient(:, iatom) = gradient(:, iatom) + fdik(:)
               IF (use_virial) THEN
                  CALL virial_pair_force(stress, -0.5_dp, fdik, rik)
               END IF
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE dEdcn_force

! **************************************************************************************************
!> \brief ...
!> \param c6ij ...
!> \param ia ...
!> \param ja ...
!> \param ik ...
!> \param jk ...
!> \param gwvec ...
!> \param c6ref ...
!> \param mrefs ...
! **************************************************************************************************
   SUBROUTINE get_c6value(c6ij, ia, ja, ik, jk, gwvec, c6ref, mrefs)
      REAL(KIND=dp), INTENT(OUT)                         :: c6ij
      INTEGER, INTENT(IN)                                :: ia, ja, ik, jk
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gwvec
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: c6ref
      INTEGER, DIMENSION(:), INTENT(IN)                  :: mrefs

      INTEGER                                            :: iref, jref
      REAL(KIND=dp)                                      :: refc6

      c6ij = 0.0_dp
      DO jref = 1, mrefs(jk)
         DO iref = 1, mrefs(ik)
            refc6 = c6ref(iref, jref, ik, jk)
            c6ij = c6ij + gwvec(iref, ia)*gwvec(jref, ja)*refc6
         END DO
      END DO

   END SUBROUTINE get_c6value

! **************************************************************************************************
!> \brief ...
!> \param c6ij ...
!> \param dc6dcn ...
!> \param dc6dq ...
!> \param ia ...
!> \param ja ...
!> \param ik ...
!> \param jk ...
!> \param gwvec ...
!> \param gwdcn ...
!> \param gwdq ...
!> \param c6ref ...
!> \param mrefs ...
! **************************************************************************************************
   SUBROUTINE get_c6derivs(c6ij, dc6dcn, dc6dq, ia, ja, ik, jk, &
                           gwvec, gwdcn, gwdq, c6ref, mrefs)
      REAL(KIND=dp), INTENT(OUT)                         :: c6ij
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT)           :: dc6dcn, dc6dq
      INTEGER, INTENT(IN)                                :: ia, ja, ik, jk
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gwvec, gwdcn, gwdq
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: c6ref
      INTEGER, DIMENSION(:), INTENT(IN)                  :: mrefs

      INTEGER                                            :: iref, jref
      REAL(KIND=dp)                                      :: refc6

      c6ij = 0.0_dp
      dc6dcn = 0.0_dp
      dc6dq = 0.0_dp
      DO jref = 1, mrefs(jk)
         DO iref = 1, mrefs(ik)
            refc6 = c6ref(iref, jref, ik, jk)
            c6ij = c6ij + gwvec(iref, ia)*gwvec(jref, ja)*refc6
            dc6dcn(1) = dc6dcn(1) + gwdcn(iref, ia)*gwvec(jref, ja)*refc6
            dc6dcn(2) = dc6dcn(2) + gwvec(iref, ia)*gwdcn(jref, ja)*refc6
            dc6dq(1) = dc6dq(1) + gwdq(iref, ia)*gwvec(jref, ja)*refc6
            dc6dq(2) = dc6dq(2) + gwvec(iref, ia)*gwdq(jref, ja)*refc6
         END DO
      END DO

   END SUBROUTINE get_c6derivs

! **************************************************************************************************
!> \brief ...
!> \param ga ...
!> \param gd ...
!> \param ev1 ...
!> \param ev2 ...
!> \param ev3 ...
!> \param ev4 ...
! **************************************************************************************************
   SUBROUTINE gerror(ga, gd, ev1, ev2, ev3, ev4)
      REAL(KIND=dp), DIMENSION(:, :)                     :: ga, gd
      REAL(KIND=dp), INTENT(OUT)                         :: ev1, ev2, ev3, ev4

      INTEGER                                            :: na, np(2)

      na = SIZE(ga, 2)
      ev1 = SQRT(SUM((gd - ga)**2)/na)
      ev2 = ev1/SQRT(SUM(gd**2)/na)*100._dp
      np = MAXLOC(ABS(gd - ga))
      ev3 = ABS(gd(np(1), np(2)) - ga(np(1), np(2)))
      ev4 = ABS(gd(np(1), np(2)))
      IF (ev4 > 1.E-6) THEN
         ev4 = ev3/ev4*100._dp
      ELSE
         ev4 = 100.0_dp
      END IF

   END SUBROUTINE gerror

! **************************************************************************************************
!> \brief ...
!> \param sa ...
!> \param sd ...
!> \param ev1 ...
!> \param ev2 ...
! **************************************************************************************************
   SUBROUTINE serror(sa, sd, ev1, ev2)
      REAL(KIND=dp), DIMENSION(3, 3)                     :: sa, sd
      REAL(KIND=dp), INTENT(OUT)                         :: ev1, ev2

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: rel

      ev1 = MAXVAL(ABS(sd - sa))
      ev2 = 0.0_dp
      DO i = 1, 3
         DO j = 1, 3
            IF (ABS(sd(i, j)) > 1.E-6_dp) THEN
               rel = ABS(sd(i, j) - sa(i, j))/ABS(sd(i, j))*100._dp
               ev2 = MAX(ev2, rel)
            END IF
         END DO
      END DO

   END SUBROUTINE serror

! **************************************************************************************************
!> \brief ...
!> \param va ...
!> \param vd ...
!> \param ev1 ...
!> \param ev2 ...
! **************************************************************************************************
   SUBROUTINE verror(va, vd, ev1, ev2)
      REAL(KIND=dp), DIMENSION(:)                        :: va, vd
      REAL(KIND=dp), INTENT(OUT)                         :: ev1, ev2

      INTEGER                                            :: i, na
      REAL(KIND=dp)                                      :: rel

      na = SIZE(va)
      ev1 = MAXVAL(ABS(vd - va))
      ev2 = 0.0_dp
      DO i = 1, na
         IF (ABS(vd(i)) > 1.E-8_dp) THEN
            rel = ABS(vd(i) - va(i))/ABS(vd(i))*100._dp
            ev2 = MAX(ev2, rel)
         END IF
      END DO

   END SUBROUTINE verror

END MODULE qs_dispersion_d4
