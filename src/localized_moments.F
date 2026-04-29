!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the calculation of moments from Wannier functions
!> \author Lukas Schreder
! **************************************************************************************************
MODULE localized_moments

   USE bibliography,                    ONLY: Schreder2021,&
                                              Schreder2024a,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_element,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_init_p, dbcsr_iterator_type, dbcsr_release, dbcsr_set, &
        dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_allocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_iter_string,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_constants,                 ONLY: use_mom_ref_com,&
                                              use_mom_ref_zero
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_types,                  ONLY: get_molecule_set_info,&
                                              molecule_type
   USE moments_utils,                   ONLY: get_reference_point
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              persistent_states_type,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_moments,                      ONLY: build_local_magmom_matrix,&
                                              build_local_moment_matrix,&
                                              build_local_moments_der_matrix,&
                                              calculate_commutator_nl_terms,&
                                              print_moments,&
                                              print_moments_nl
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE wannier_states_types,            ONLY: wannier_centres_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** public subroutines ***

   PUBLIC :: calculate_kg_moments,&
             calculate_localized_moments

! *** body ***

CONTAINS

! **************************************************************************************************
!> \brief Calculates localised electric and magnetic moments from Wannier functions or Kim-Gordon
!>        subdensities
!> \param qs_env ...
!> \param mo_local ...
! **************************************************************************************************
   SUBROUTINE calculate_localized_moments(qs_env, maxmom, mo_local)

      TYPE(qs_environment_type), POINTER  :: qs_env
      INTEGER                             :: maxmom
      TYPE(cp_fm_type), OPTIONAL, POINTER :: mo_local(:)

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_localized_moments', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)  :: name
      INTEGER                               :: handle, iatom, imo, isubset, jatom, nao, &
                                               natoms, nmo, nmoments, nsubsets, nwcs
      INTEGER, ALLOCATABLE                  :: first_atoms(:), first_basis(:), &
                                               last_atoms(:), last_basis(:), &
                                               nbasis(:), nelectrons(:), &
                                               wcs_to_atoms(:)
      INTEGER, POINTER                      :: states(:)
      LOGICAL                               :: do_periodic
      REAL(dp)                              :: dist, r_min, rcc(3), rcc0(3)
      REAL(dp), ALLOCATABLE                 :: moments(:)
      REAL(dp), POINTER                     :: rccs(:, :), ref_point(:)
      TYPE(cell_type), POINTER              :: cell
      TYPE(cp_logger_type), POINTER         :: logger
      TYPE(mp_para_env_type), POINTER       :: para_env
      TYPE(dbcsr_type), POINTER             :: matrix_s(:)
      TYPE(dft_control_type), POINTER       :: dft_control
      TYPE(molecule_type), POINTER          :: molecules(:)
      TYPE(neighbor_list_set_type), POINTER :: sab_orb(:)
      TYPE(particle_type), POINTER          :: particles(:)
      TYPE(qs_kind_type), POINTER           :: qs_kinds(:)
      TYPE(qs_rho_type), POINTER            :: rho
      TYPE(section_vals_type), POINTER      :: input
      TYPE(wannier_centres_type), POINTER   :: wcs(:)

      CALL timeset(routineN, handle)

      CALL cite_reference(Schreder2021)
      CALL cite_reference(Schreder2024a)

      logger => cp_get_default_logger()
      iteration = iteration + 1

      ! environment infos
      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      input=input, &
                      matrix_s=matrix_s, &
                      molecule_set=molecules, &
                      para_env=para_env, &
                      particle_set=particles, &
                      qs_kind_set=qs_kinds, &
                      rho=rho, &
                      sab_orb=sab_orb)

      ! charge(3), eldip(3), elqup(6), mdip(3), mom(3), elqupvel(6), mdipnl(3), momnl(3), elqupvelnl(6)
      !nmoments = 36
      natoms = SIZE(particles)
      nsubsets = SIZE(molecules)
      ALLOCATE (first_atoms(nsubsets), first_basis(nsubsets), last_atoms(nsubsets), &
                last_basis(nsubsets), nbasis(nsubsets), nelectrons(nsubsets), moments(nmoments), &
                rccs(nsubsets, 3), ref_point(3))

      CALL get_molecule_set_info(molecules, &
                                 mol_to_first_atom=first_atoms, &
                                 mol_to_last_atom=last_atoms, &
                                 mol_to_nbasis=nbasis, &
                                 mol_to_nelectrons=nelectrons)

      ! reference points
      DO isubset = 1, nsubsets
         CALL get_reference_point(rccs(isubset, :), &
                                  qs_env=qs_env, &
                                  reference=use_mom_ref_com, &
                                  ifirst=first_atoms(isubset), &
                                  ilast=last_atoms(isubset))
      END DO

      IF (iteration == 1) THEN
         CALL cp_fm_get_info(mo_local(1), ncol_global=nmo, nrow_global=nao)
         CALL get_qs_env(qs_env, wanniercentres=wcs)
         ALLOCATE (persistent_states(nsubsets, 1), wcs_to_atoms(nmo))
         ! states -> subset
         DO imo = 1, nmo
            r_min = HUGE(0._dp)
            jatom = 0
            DO iatom = 1, natoms
               dist = NORM2(pbc(wcs(1)%centres(:, imo) - particles(iatom)%r, cell))
               IF (dist < r_min) THEN
                  jatom = iatom
                  r_min = dist
               END IF
            END DO
            wcs_to_atoms(imo) = jatom
         END DO

         ! subsets -> states map
         nelectrons = 0
         DO isubset = 1, nsubsets
            DO imo = 1, nmo
               IF ((first_atoms(isubset) <= wcs_to_atoms(imo)) .AND. &
                   (last_atoms(isubset) >= wcs_to_atoms(imo))) THEN
                  nelectrons(isubset) = nelectrons(isubset) + 1
               END IF
            END DO
         END DO

         DO isubset = 1, nsubsets
            ALLOCATE (persistent_states(isubset, 1)%states(nelectrons(isubset)))
            nwcs = 1
            DO imo = 1, nmo
               IF ((first_atoms(isubset) <= wcs_to_atoms(imo)) .AND. &
                   (last_atoms(isubset) >= wcs_to_atoms(imo))) THEN
                  persistent_states(isubset, 1)%states(nwcs) = imo
                  nwcs = nwcs + 1
               END IF
            END DO
            PRINT *, "dbg| states(isubset)", isubset, SIZE(persistent_states(isubset, 1)%states), &
               persistent_states(isubset, 1)%states
         END DO
         DEALLOCATE (wcs_to_atoms)

         ! create operators
         CALL create_operators(qs_env, rcc, rcc0, rccs, mo_local)
      END IF

      CALL calculate_loc_moments(qs_env, rcc, mo_local, opcode=-1, tag="ref")

      ! subsystems
      DO isubset = 1, nsubsets
         states => persistent_states(isubset, 1)%states
         name = molecules(isubset)%molecule_kind%name
         CALL calculate_loc_moments(qs_env, rccs(isubset, :), mo_local, isubset, states, &
                                    TRIM(name))
      END DO

      DEALLOCATE (first_atoms, first_basis, last_atoms, last_basis, moments, nelectrons, rccs)

      CALL timestop(handle)

   END SUBROUTINE calculate_localized_moments
! **************************************************************************************************

   SUBROUTINE create_operators(qs_env, rcc, rcc0, rccs, mos)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp)                                           :: rcc(3), rcc0(3)
      REAL(dp), POINTER                                  :: rccs(:, :)
      TYPE(cp_fm_type), POINTER                          :: mos(:)

      CHARACTER(*), PARAMETER                            :: routineN = 'create_operators'

      INTEGER                                            :: handle, i, isubset, j, nao, nmm, nmo, &
                                                            nsubsets
      REAL(dp), POINTER                                  :: myrcc(:)
      TYPE(cp_fm_struct_type), POINTER                   :: o_struct
      TYPE(cp_fm_type), POINTER                          :: or
      TYPE(dbcsr_type), POINTER                          :: matrix_s(:), mop(:), pop(:), rop(:), &
                                                            ropder(:, :)
      TYPE(neighbor_list_set_type), POINTER              :: sab_all(:), sab_orb(:)

      CALL timeset(routinen, handle)

      CALL get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb)
      CALL cp_fm_get_info(mos(1), nrow_global=nao, ncol_global=nmo)

      nsubsets = SIZE(rccs, 1)

      ALLOCATE ( &
         ro(-1:SIZE(rccs, 1), 9), &
         roder(-1:SIZE(rccs, 1), 3, 3), &
         po(-1:SIZE(rccs, 1), 3), &
         mo(-1:SIZE(rccs, 1), 3))

      CALL cp_fm_struct_create(o_struct, nrow_global=nao, ncol_global=nao, &
                               template_fmstruct=mos(1)%matrix_struct)
      CALL cp_fm_create(or, o_struct)  ! tmp

      DO isubset = -1, nsubsets

         PRINT *, "create operator", isubset

         IF (isubset == -1) THEN
            ALLOCATE (myrcc(3))
            myrcc(:) = rcc(:)
         ELSE IF (isubset == 0) THEN
            ALLOCATE (myrcc(3))
            myrcc(:) = rcc0(:)
         ELSE
            myrcc => rccs(isubset, :)
         END IF

         ! r and derivatives
         NULLIFY (rop)
         CALL dbcsr_allocate_matrix_set(rop, nm)
         DO i = 1, nm
            ALLOCATE (rop(i))
            CALL dbcsr_create(rop(i), TEMPLATE=matrix_s(1), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(rop(i), sab_all)
            CALL dbcsr_set(rop(i), 0._dp)
         END DO
         NULLIFY (ropder)
         CALL dbcsr_allocate_matrix_set(ropder, 3, 3)
         DO i = 1, 3
            DO j = 1, 3
               CALL dbcsr_init_p(ropder(i, j))
               CALL dbcsr_create(ropder(i, j), TEMPLATE=matrix_s(1), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL cp_dbcsr_alloc_block_from_nbl(ropder(i, j), sab_all)
               CALL dbcsr_set(ropder(i, j), 0._dp)
            END DO
         END DO
         CALL build_local_moments_der_matrix(qs_env, rop, ropder, 2, 1, myrcc)
         DO i = 1, nm
            CALL cp_cfm_create(ro(isubset, i), o_struct)
            CALL copy_dbcsr_to_fm(rop(i), or)
            CALL dbcsr_release(rop(i))
            DEALLOCATE (rop(i))
            CALL cp_fm_to_cfm(or, mtarget=ro(isubset, i))
         END DO
         DEALLOCATE (rop)
         DO i = 1, 3
            DO j = 1, 3
               CALL cp_cfm_create(roder(isubset, i, j), o_struct)
               CALL copy_dbcsr_to_fm(ropder(i, j), or)
               CALL dbcsr_release(ropder(i, j))
               DEALLOCATE (ropder(i, j))
               CALL cp_fm_to_cfm(or, mtarget=roder(isubset, i, j))
            END DO
         END DO
         DEALLOCATE (ropder)

         ! m
         NULLIFY (mop)
         CALL build_local_magmom_matrix(qs_env, mop, 1, myrcc)
         nmm = SIZE(mop)
         DO i = 1, nmm
            CALL cp_cfm_create(mo(isubset, i), o_struct)
            CALL copy_dbcsr_to_fm(mop(i), or)
            CALL dbcsr_release(mop(i))
            DEALLOCATE (mop(i))
            CALL cp_fm_to_cfm(or, mtarget=mo(isubset, i))
         END DO
         DEALLOCATE (mop)

         ! p
         NULLIFY (pop)
         CALL dbcsr_allocate_matrix_set(pop, 3)
         DO i = 1, 3
            CALL dbcsr_init_p(pop(i))
            CALL dbcsr_create(pop(i), TEMPLATE=matrix_s(1), &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(pop(i), sab_orb)
            CALL dbcsr_set(pop(i), 0._dp)
         END DO
         CALL build_lin_mom_matrix(qs_env, pop)
         DO i = 1, 3
            CALL cp_cfm_create(po(isubset, i), o_struct)
            CALL copy_dbcsr_to_fm(pop(i), or)
            CALL dbcsr_release(pop(i))
            DEALLOCATE (pop(i))
            CALL cp_fm_to_cfm(or, mtarget=po(isubset, i))
         END DO
         DEALLOCATE (pop)

         IF (isubset == -1) THEN
            DEALLOCATE (myrcc)
         ELSE IF (isubset == 0) THEN
            DEALLOCATE (myrcc)
         END IF
      END DO

      CALL cp_fm_release(or)  ! tmp
      CALL cp_fm_struct_release(o_struct)

      CALL timestop(handle)

   END SUBROUTINE create_operators
! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rcc ...
!> \param mos ...
!> \param opcode ...
!> \param statelist ...
!> \param tag ...
! **************************************************************************************************
   SUBROUTINE calculate_loc_moments(qs_env, rcc, mos, opcode, statelist, tag)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp)                                           :: rcc(3)
      TYPE(cp_fm_type), POINTER                          :: mos(:)
      INTEGER                                            :: opcode
      INTEGER, OPTIONAL, POINTER                         :: statelist(:)
      CHARACTER(*), OPTIONAL                             :: tag

      CHARACTER(*), PARAMETER :: routineN = 'calculate_loc_moments'

      COMPLEX(dp)                                        :: di
      INTEGER                                            :: file_unit, handle, i, j, k, nao, nmo, ns
      INTEGER, POINTER                                   :: states(:)
      LOGICAL                                            :: do_rtp
      REAL(dp)                                           :: d, m(3), mu(3), muv(3), q, th(6), &
                                                            thder(9), thv(6), zwfc
      REAL(dp), ALLOCATABLE                              :: rrv(:), rv(:), rxrv(:)
      TYPE(cp_cfm_type), POINTER                         :: o, ov, v, vov
      TYPE(cp_fm_struct_type), POINTER                   :: v_struct, vov_struct
      TYPE(cp_fm_type), POINTER                          :: vi, vr
      TYPE(dbcsr_type), POINTER                          :: matrix_s(:)
      TYPE(neighbor_list_set_type), POINTER              :: sab_all(:), sab_orb(:)

      CALL timeset(routinen, handle)

      CALL get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb)
      CALL cp_fm_get_info(mos(1), nrow_global=nao, ncol_global=nmo)
      zwfc = 2._dp
      do_rtp = (SIZE(mos) == 2)

      IF (PRESENT(statelist)) THEN
         ns = SIZE(statelist)
         states => statelist
      ELSE
         ns = nmo
         ALLOCATE (states(ns))
         DO i = 1, ns
            states(i) = i
         END DO
      END IF

      ! matrices
      CALL cp_fm_struct_create(v_struct, nrow_global=nao, ncol_global=ns, &
                               template_fmstruct=mos(1)%matrix_struct)
      CALL cp_fm_struct_create(vov_struct, nrow_global=ns, ncol_global=ns, &
                               template_fmstruct=mos(1)%matrix_struct)

      CALL cp_cfm_create(ov, v_struct, name="ov")
      CALL cp_cfm_create(vov, vov_struct, name="vov")
      CALL cp_fm_struct_release(vov_struct)

      ! MOs
      CALL cp_fm_create(vr, v_struct, name="vr")
      CALL cp_fm_create(vi, v_struct, name="vi")
      CALL cp_cfm_create(v, v_struct, name="v")
      CALL cp_fm_struct_release(v_struct)
      DO i = 1, ns
         CALL cp_fm_to_fm(mos(1), vr, 1, states(i), i)
         IF (do_rtp) THEN
            CALL cp_fm_to_fm(mos(2), vi, 1, states(i), i)
            CALL cp_fm_to_cfm(vr, vi, v)
         ELSE
            CALL cp_fm_to_cfm(vr, mtarget=v)
         END IF
      END DO
      CALL cp_fm_release(vr)
      CALL cp_fm_release(vi)

      ! charge
      q = -ns*zwfc
      CALL open_file("q_"//tag//".csv", unit_number=file_unit, file_action="write", &
                     file_position="append", file_status="unknown")
      IF (file_unit > 0) THEN
         WRITE (file_unit, "(E14.8)") q
      END IF
      CALL close_file(file_unit)

      ! mu
      mu = 0._dp
      DO i = 1, 3
         o => ro(opcode, i)
         CALL cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o, v, (0._dp, 0._dp), ov)
         CALL cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), v, ov, (0._dp, 0._dp), vov)
         d = 0._dp
         DO j = 1, ns
            CALL cp_cfm_get_element(vov, j, j, di)
            d = d + REAL(di)
         END DO
         mu(i) = -zwfc*d
      END DO
      CALL open_file("mu_"//tag//".csv", unit_number=file_unit, file_action="write", &
                     file_position="append", file_status="unknown")
      IF (file_unit > 0) THEN
         WRITE (file_unit, "(2(E14.8,','),E14.8)") mu
      END IF
      CALL close_file(file_unit)

      ! muv
      muv = 0._dp
      DO i = 1, 3
         o => po(opcode, i)
         CALL cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o, v, (0._dp, 0._dp), ov)
         CALL cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), v, ov, (0._dp, 0._dp), vov)
         d = 0._dp
         DO j = 1, ns
            CALL cp_cfm_get_element(vov, j, j, di)
            d = d + AIMAG(di)
         END DO
         muv(i) = -zwfc*d
      END DO
      CALL open_file("muv_"//tag//".csv", unit_number=file_unit, file_action="write", &
                     file_position="append", file_status="unknown")
      IF (file_unit > 0) THEN
         WRITE (file_unit, "(2(E14.8,','),E14.8)") muv
      END IF
      CALL close_file(file_unit)

      ! m
      m = 0._dp
      DO i = 1, 3
         o => mo(opcode, i)
         CALL cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o, v, (0._dp, 0._dp), ov)
         CALL cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), v, ov, (0._dp, 0._dp), vov)
         d = 0._dp
         DO j = 1, ns
            CALL cp_cfm_get_element(vov, j, j, di)
            d = d + AIMAG(di)
         END DO
         m(i) = -zwfc*d
      END DO
      CALL open_file("m_"//tag//".csv", unit_number=file_unit, file_action="write", &
                     file_position="append", file_status="unknown")
      IF (file_unit > 0) THEN
         WRITE (file_unit, "(2(E14.8,','),E14.8)") m
      END IF
      CALL close_file(file_unit)

      ! th
      th = 0._dp
      DO i = 1, 6
         o => ro(opcode, i + 3)
         CALL cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o, v, (0._dp, 0._dp), ov)
         CALL cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), v, ov, (0._dp, 0._dp), vov)
         d = 0._dp
         DO j = 1, ns
            CALL cp_cfm_get_element(vov, j, j, di)
            d = d + zwfc*REAL(di)
         END DO
         th(i) = -d
      END DO
      CALL open_file("th_"//tag//".csv", unit_number=file_unit, file_action="write", &
                     file_position="append", file_status="unknown")
      IF (file_unit > 0) THEN
         WRITE (file_unit, "(5(E14.8,','),E14.8)") th
      END IF
      CALL close_file(file_unit)

      ! thv
      thv = 0._dp
      DO i = 1, 3
         DO j = 1, 3
            o => roder(opcode, i, j)
            CALL cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o, v, (0._dp, 0._dp), ov)
            CALL cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), v, ov, (0._dp, 0._dp), vov)
            d = 0._dp
            DO k = 1, ns
               CALL cp_cfm_get_element(vov, k, k, di)
               d = d - AIMAG(di)
               IF (i == j) THEN
                  d = d + REAL(di)
               END IF
            END DO
            thder((i - 1)*3 + j) = zwfc*d
         END DO
      END DO
      thv(1) = -2*thder(1) - q
      thv(2) = -thder(2) - thder(4)
      thv(3) = -thder(3) - thder(7)
      thv(4) = -2*thder(5) - q
      thv(5) = -thder(6) - thder(8)
      thv(6) = -2*thder(9) - q
      CALL open_file("thv_"//tag//".csv", unit_number=file_unit, file_action="write", &
                     file_position="append", file_status="unknown")
      IF (file_unit > 0) THEN
         WRITE (file_unit, "(5(E14.8,','),E14.8)") thv
      END IF
      CALL close_file(file_unit)

      IF (.NOT. PRESENT(statelist)) THEN
         DEALLOCATE (states)
      END IF
      CALL cp_cfm_release(v)
      CALL cp_cfm_release(ov)
      CALL cp_cfm_release(vov)

      ! nl terms
      ALLOCATE (rv(3), rrv(6), rxrv(3))
      rv = 0._dp
      rrv = 0._dp
      rxrv = 0._dp
      CALL calculate_commutator_nl_terms(qs_env, rv, rxrv, rrv, rcc)

      CALL open_file("rv_"//tag//".csv", unit_number=file_unit, file_action="write", &
                     file_position="append", file_status="unknown")
      IF (file_unit > 0) THEN
         WRITE (file_unit, "(2(E14.8,','),E14.8)") rv
      END IF
      CALL close_file(file_unit)

      CALL open_file("rrv_"//tag//".csv", unit_number=file_unit, file_action="write", &
                     file_position="append", file_status="unknown")
      IF (file_unit > 0) THEN
         WRITE (file_unit, "(5(E14.8,','),E14.8)") rrv
      END IF
      CALL close_file(file_unit)

      CALL open_file("rxrv_"//tag//".csv", unit_number=file_unit, file_action="write", &
                     file_position="append", file_status="unknown")
      IF (file_unit > 0) THEN
         WRITE (file_unit, "(2(E14.8,','),E14.8)") rxrv
      END IF
      CALL close_file(file_unit)

      CALL open_file("rcc_"//tag//".csv", unit_number=file_unit, file_action="write", &
                     file_position="append", file_status="unknown")
      IF (file_unit > 0) THEN
         WRITE (file_unit, "(2(E14.8,','),E14.8)") rcc
      END IF
      CALL close_file(file_unit)

      CALL timestop(handle)

   END SUBROUTINE calculate_loc_moments

! **************************************************************************************************
!> \brief ...
!> \param moments ...
!> \param qs_env ...
!> \param rcc ...
!> \param mo_local ...
!> \param particles ...
!> \param qs_kinds ...
!> \param cell ...
!> \param first_atom ...
!> \param last_atom ...
!> \param statelist ...
!> \param periodic ...
!> \param ref ...
! **************************************************************************************************
   SUBROUTINE calculate_localized_loc_moments(moments, qs_env, rcc, mo_local, particles, qs_kinds, &
                                              cell, first_atom, last_atom, statelist, periodic, ref)

      REAL(dp)                                           :: moments(:)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp)                                           :: rcc(3)
      TYPE(cp_fm_type), POINTER                          :: mo_local(:)
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kinds(:)
      TYPE(cell_type), POINTER                           :: cell
      INTEGER                                            :: first_atom, last_atom
      INTEGER, OPTIONAL, POINTER                         :: statelist(:)
      LOGICAL, OPTIONAL                                  :: periodic, ref

      COMPLEX(dp)                                        :: cd
      INTEGER                                            :: i, j, nao, nmm, nmo, ns
      INTEGER, POINTER                                   :: states(:)
      LOGICAL                                            :: do_periodic, do_ref, do_rtp
      REAL(dp)                                           :: d, zwfc
      REAL(dp), ALLOCATABLE                              :: nlrrv(:), nlrv(:), nlrxrv(:)
      TYPE(cp_cfm_type), POINTER                         :: mago(:), momo(:), o(:), oder(:, :)
      TYPE(cp_cfm_type), POINTER                         :: m, mom, om
      TYPE(cp_fm_struct_type), POINTER                   :: m_struct, mom_struct, o_struct
      TYPE(cp_fm_type), POINTER                          :: mi, mr, or
      TYPE(dbcsr_type), POINTER                          :: magop(:), matrix_s(:), momop(:), op(:), &
                                                            opder(:, :)
      TYPE(neighbor_list_set_type), POINTER              :: sab_all(:), sab_orb(:)

      !print *, "dbg| loc"
      !print *, "dbg| calc_loc_mom rcc", rcc

      ! info
      IF (PRESENT(periodic)) THEN
         do_periodic = periodic
      ELSE
         do_periodic = .FALSE.
      END IF
      IF (PRESENT(ref)) THEN
         do_ref = ref
      ELSE
         do_ref = .FALSE.
      END IF
      CALL get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb)
      do_rtp = (SIZE(mo_local) == 2)
      CALL cp_fm_get_info(mo_local(1), nrow_global=nao, ncol_global=nmo)
      zwfc = 2._dp
      ALLOCATE (states(nmo))
      IF (PRESENT(statelist)) THEN
         ns = SIZE(statelist)
         DO i = 1, ns
            states(i) = statelist(i)
         END DO
      ELSE
         ns = nmo
         DO i = 1, ns
            states(i) = i
         END DO
      END IF

      ! structs
      CALL cp_fm_struct_create(m_struct, nrow_global=nao, ncol_global=ns, &
                               template_fmstruct=mo_local(1)%matrix_struct)
      CALL cp_fm_struct_create(o_struct, nrow_global=nao, ncol_global=nao, &
                               template_fmstruct=mo_local(1)%matrix_struct)
      CALL cp_fm_struct_create(mom_struct, nrow_global=ns, ncol_global=ns, &
                               template_fmstruct=mo_local(1)%matrix_struct)

      ! get only this subsystem's states
      CALL cp_fm_create(mr, m_struct, name="mr")
      CALL cp_fm_set_all(mr, 0._dp)
      IF (do_rtp) THEN
         CALL cp_fm_create(mi, m_struct, name="mi")
         CALL cp_fm_set_all(mi, 0._dp)
      END IF
      CALL cp_cfm_create(m, m_struct, name="m")
      DO i = 1, ns
         CALL cp_fm_to_fm(mo_local(1), mr, 1, states(i), i)
         IF (do_rtp) THEN
            CALL cp_fm_to_fm(mo_local(2), mi, 1, states(i), i)
            CALL cp_fm_to_cfm(mr, mi, m)
         ELSE
            CALL cp_fm_to_cfm(mr, mtarget=m)
         END IF
      END DO
      CALL cp_fm_release(mr)
      IF (do_rtp) THEN
         CALL cp_fm_release(mi)
      END IF

     !! get the operators
      CALL cp_fm_create(or, o_struct, name="or")
      ! el moments with derivatives
      NULLIFY (op)
      CALL dbcsr_allocate_matrix_set(op, nm)
      DO i = 1, nm
         CALL dbcsr_init_p(op(i))
         CALL dbcsr_create(op(i), TEMPLATE=matrix_s(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_alloc_block_from_nbl(op(i), sab_all)
         CALL dbcsr_set(op(i), 0._dp)
      END DO
      NULLIFY (opder)
      CALL dbcsr_allocate_matrix_set(opder, 3, 3)
      DO i = 1, 3
         DO j = 1, 3
            CALL dbcsr_init_p(opder(i, j))
            CALL dbcsr_create(opder(i, j), TEMPLATE=matrix_s(1), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_alloc_block_from_nbl(opder(i, j), sab_all)
            CALL dbcsr_set(opder(i, j), 0._dp)
         END DO
      END DO
      ALLOCATE (o(9))
      IF (ref) THEN
         CALL build_local_moments_der_matrix(qs_env, op, opder, 2, 1, rcc)
      ELSE
         CALL build_local_moment_matrix(qs_env, op, 2, rcc)
      END IF
      DO i = 1, nm
         CALL cp_cfm_create(o(i), o_struct)
         IF (ref) THEN
            CALL copy_dbcsr_to_fm(op(i), or)
         END IF
         CALL dbcsr_release(op(i))
         DEALLOCATE (op(i))
         IF (ref) THEN
            CALL cp_fm_to_cfm(or, mtarget=o(i))
         ELSE
            CALL cp_cfm_set_all(o(i), (0._dp, 0._dp))
         END IF
      END DO
      DEALLOCATE (op)
      ALLOCATE (oder(3, 3))
      DO i = 1, 3
         DO j = 1, 3
            CALL cp_cfm_create(oder(i, j), o_struct)
            IF (ref) THEN
               CALL copy_dbcsr_to_fm(opder(i, j), or)
            END IF
            CALL dbcsr_release(opder(i, j))
            DEALLOCATE (opder(i, j))
            IF (ref) THEN
               CALL cp_fm_to_cfm(or, mtarget=oder(i, j))
            ELSE
               CALL cp_cfm_set_all(oder(i, j), (0._dp, 0._dp))
            END IF
         END DO
      END DO
      DEALLOCATE (opder)

      ! momentum
      NULLIFY (momop)
      CALL dbcsr_allocate_matrix_set(momop, 3)
      DO i = 1, 3
         CALL dbcsr_init_p(momop(i))
         CALL dbcsr_create(momop(i), TEMPLATE=matrix_s(1), &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(momop(i), sab_orb)
         CALL dbcsr_set(momop(i), 0._dp)
      END DO
      ALLOCATE (momo(3))
      CALL build_lin_mom_matrix(qs_env, momop)
      DO i = 1, 3
         CALL cp_cfm_create(momo(i), o_struct)
         CALL copy_dbcsr_to_fm(momop(i), or)
         CALL dbcsr_release(momop(i))
         DEALLOCATE (momop(i))
         CALL cp_fm_to_cfm(or, mtarget=momo(i))
      END DO
      DEALLOCATE (momop)

      ! mag moments
      NULLIFY (magop)
      CALL dbcsr_allocate_matrix_set(magop, 3)
      DO i = 1, 3
         CALL dbcsr_init_p(magop(i))
         CALL dbcsr_create(magop(i), TEMPLATE=matrix_s(1), &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(magop(i), sab_orb)
         CALL dbcsr_set(magop(i), 0._dp)
      END DO
      NULLIFY (op)
      ALLOCATE (mago(3))
      CALL build_local_magmom_matrix(qs_env, magop, 1, ref_point=rcc)
      nmm = SIZE(magop)
      DO i = 1, nmm
         CALL cp_cfm_create(mago(i), o_struct, name="mago")
         IF (ref) THEN
            CALL copy_dbcsr_to_fm(magop(i), or)
         END IF
         CALL dbcsr_release(magop(i))
         DEALLOCATE (magop(i))
         IF (ref) THEN
            CALL cp_fm_to_cfm(or, mtarget=mago(i))
         ELSE
            CALL cp_cfm_set_all(mago(i), (0._dp, 0._dp))
         END IF
      END DO
      DEALLOCATE (magop)
      CALL cp_fm_release(or)

      ! get the nonlocal terms
      ALLOCATE (nlrv(3), nlrrv(6), nlrxrv(3))
      nlrv = 0._dp
      nlrrv = 0._dp
      nlrxrv = 0._dp
      IF (ref) THEN
         CALL calculate_commutator_nl_terms(qs_env, nlrv, nlrxrv, nlrrv, rcc)
      END IF

      ! nuclear contribution
      CALL calculate_localized_nuc_moments(moments, particles, qs_kinds, cell, first_atom, &
                                           last_atom, rcc)

      ! expectation values
      CALL cp_cfm_create(om, m_struct, name="om")
      CALL cp_cfm_create(mom, mom_struct, name="mom")
      CALL cp_fm_struct_release(m_struct)
      CALL cp_fm_struct_release(o_struct)
      CALL cp_fm_struct_release(mom_struct)

      ! charge
      !call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, m, (0._dp, 0._dp), mom)
      !d = 0.
      !do i = 1, ns
      !   call cp_cfm_get_element(mom, i, i, cd)
      !   moments(2) = moments(2) + real(cd)
      !end do
      !moments(2) = -zwfc * ns! moments(2)
      !moments(3) = moments(1) + moments(2)

      ! edip
      DO i = 1, 3
         CALL cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o(i), m, (0._dp, 0._dp), om)
         CALL cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
         d = 0.
         DO j = 1, ns
            CALL cp_cfm_get_element(mom, j, j, cd)
            d = d + REAL(cd)
         END DO
         moments(i + 3) = moments(i + 3) - zwfc*d
      END DO
      ! vel
      !do i = 1, 3
      !   call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), momo(i), m, (0._dp, 0._dp), om)
      !   call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
      !   d = 0.
      !   do j = 1, ns
      !      call cp_cfm_get_element(mom, j, j, cd)
      !      d = d + aimag(cd)
      !   end do
      !   moments(i+15) = moments(i+15) - zwfc * d
      !end do
     !! nl
      !moments(28:30) = moments(16:18) + nlrv

     !! mdip
      !do i = 1, 3
      !   call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), mago(i), m, (0._dp, 0._dp), om)
      !   call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
      !   d = 0.
      !   do j = 1, ns
      !      call cp_cfm_get_element(mom, j, j, cd)
      !      d = d + aimag(cd)
      !   end do
      !   moments(i+12) = moments(i+12) - zwfc * d
      !end do
     !! nl
      !moments(25:27) = moments(13:15) + nlrxrv

     !! equp
      !do i = 4, 9
      !   call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o(i), m, (0._dp, 0._dp), om)
      !   call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
      !   d = 0.
      !   do j = 1, ns
      !      call cp_cfm_get_element(mom, j, j, cd)
      !      d = d + zwfc * real(cd)
      !   end do
      !   moments(i+3) = moments(i+3) - d
      !end do
     !! vel
      !qupoleder = 0._dp
      !do i = 1, 3
      !   do j = 1, 3
      !      call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), oder(i, j), m, (0._dp, 0._dp), om)
      !      call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
      !      d = 0.
      !      do k = 1, ns
      !         call cp_cfm_get_element(mom, k, k, cd)
      !         d = d - aimag(cd)
      !         if (i == j) then
      !            d = d + real(cd)
      !         end if
      !      end do
      !      qupoleder((i-1)*3+j) = qupoleder((i-1)*3+j) + zwfc * d
      !   end do
      !end do
      !moments(19) = -2*qupoleder(1) - moments(2)
      !moments(20) = -qupoleder(2) - qupoleder(4)
      !moments(21) = -qupoleder(3) - qupoleder(7)
      !moments(22) = -2*qupoleder(5) - moments(2)
      !moments(23) = -qupoleder(6) - qupoleder(8)
      !moments(24) = -2*qupoleder(9) - moments(2)
     !! nl
      !moments(31:36) = moments(19:24) + nlrrv

      ! cleanup
      CALL cp_cfm_release(m)
      CALL cp_cfm_release(om)
      CALL cp_cfm_release(mom)
      DO i = 1, SIZE(o)
         CALL cp_cfm_release(o(i))
      END DO
      DEALLOCATE (o)
      DO i = 1, SIZE(oder, 1)
         DO j = 1, SIZE(oder, 2)
            CALL cp_cfm_release(oder(i, j))
         END DO
      END DO
      DEALLOCATE (oder)
      DO i = 1, SIZE(momo)
         CALL cp_cfm_release(momo(i))
      END DO
      DEALLOCATE (momo)
      DO i = 1, SIZE(mago)
         CALL cp_cfm_release(mago(i))
      END DO
      DEALLOCATE (mago)
      DEALLOCATE (nlrv, nlrrv, nlrxrv, states)

   END SUBROUTINE calculate_localized_loc_moments

! **************************************************************************************************
! **************************************************************************************************
   subroutine calculate_localized_kg_moments(moments, qs_env, rcc, rho, particles, qs_kinds, cell, &
        first_atom, last_atom, do_periodic)

     real(dp)                            :: moments(:)
     type(qs_environment_type), pointer  :: qs_env
     real(dp)                            :: rcc(3)
     type(qs_rho_type), pointer          :: rho
     type(particle_type), pointer        :: particles(:)
     type(qs_kind_type), pointer         :: qs_kinds(:)
     type(cell_type), pointer            :: cell
     integer                             :: first_atom, last_atom
     logical  :: do_periodic

     integer                                  :: i, j
     logical                                  :: found, tr
     real(dp)                                 :: focc, qupoleder(9), strace, trace
     real(dp), allocatable                    :: nlrv(:), nlrrv(:), nlrxrv(:)
     real(dp), pointer                        :: blk_p(:)
     type(dbcsr_iterator_type), pointer       :: iterator
     type(dbcsr_type), pointer              :: o(:), oder(:, :), mago(:), momo(:), matrix_s(:), &
          rho_ao(:), rho_ao_im(:)
     type(dbcsr_type), pointer                :: densmat, densmat_im, tmp
     type(neighbor_list_set_type), pointer  :: sab_all(:), sab_orb(:)

     ! info
     call get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb)
     call qs_rho_get(rho, rho_ao=rho_ao, rho_ao_im=rho_ao_im)
     focc = 1._dp

     ! get real and imaginary part of the subdensity
     allocate(densmat)
     call dbcsr_copy(densmat, rho_ao(1))
     allocate(iterator)
     call dbcsr_iterator_start(iterator, densmat)
     do while (dbcsr_iterator_blocks_left(iterator))
        call dbcsr_iterator_next_block(iterator, i, j, blk_p)
        if ((i < first_atom) .or. (i > last_atom) .or. (j < first_atom) .or. (j > last_atom)) then
           call dbcsr_get_block_p(densmat, i, j, blk_p, tr, found)
           blk_p = 0._dp
        end if
     end do
     call dbcsr_iterator_stop(iterator)
     allocate(densmat_im)
     call dbcsr_copy(densmat_im, rho_ao_im(1))
     call dbcsr_iterator_start(iterator, densmat_im)
     do while (dbcsr_iterator_blocks_left(iterator))
        call dbcsr_iterator_next_block(iterator, i, j, blk_p)
        if ((i < first_atom) .or. (i > last_atom) .or. (j < first_atom) .or. (j > last_atom)) then
           call dbcsr_get_block_p(densmat_im, i, j, blk_p, tr, found)
           blk_p = 0._dp
        end if
     end do
     call dbcsr_iterator_stop(iterator)
     deallocate(iterator)

     ! get the operators
     ! el moments with derivatives
     nullify(o)
     call dbcsr_allocate_matrix_set(o, nm)
     do i = 1, nm
        call dbcsr_init_p(o(i))
        call dbcsr_create(o(i), template=matrix_s(1), &
             matrix_type=dbcsr_type_no_symmetry)
        call cp_dbcsr_alloc_block_from_nbl(o(i), sab_all)
        call dbcsr_set(o(i), 0._dp)
     end do
     nullify(oder)
     call dbcsr_allocate_matrix_set(oder, 3, 3)
     do i = 1, 3
        do j = 1, 3
           call dbcsr_init_p(oder(i, j))
           call dbcsr_create(oder(i, j), template=matrix_s(1), &
                matrix_type=dbcsr_type_no_symmetry)
           call cp_dbcsr_alloc_block_from_nbl(oder(i, j), sab_all)
           call dbcsr_set(oder(i, j), 0._dp)
        end do
     end do
     call build_local_moments_der_matrix(qs_env, o, oder, 2, 1, rcc)

     ! momentum
     nullify(momo)
     call dbcsr_allocate_matrix_set(momo, 3)
     do i = 1, 3
        call dbcsr_init_p(momo(i))
        call dbcsr_create(momo(i), template=matrix_s(1), &
             matrix_type=dbcsr_type_antisymmetric)
        call cp_dbcsr_alloc_block_from_nbl(momo(i), sab_orb)
        call dbcsr_set(momo(i), 0._dp)
     end do
     call build_lin_mom_matrix(qs_env, momo)

     ! mag moments
     nullify(mago)
     call dbcsr_allocate_matrix_set(mago, 3)
     do i = 1, 3
        call dbcsr_init_p(mago(i))
        call dbcsr_create(mago(i), template=matrix_s(1), &
             matrix_type=dbcsr_type_antisymmetric)
        call cp_dbcsr_alloc_block_from_nbl(mago(i), sab_orb)
        call dbcsr_set(mago(i), 0._dp)
     end do
     call build_local_magmom_matrix(qs_env, mago, nmm, ref_point=rcc)

     ! get the nonlocal terms
     allocate(nlrv(3), nlrrv(6), nlrxrv(3))
     nlrv = 0._dp
     nlrrv = 0._dp
     nlrxrv = 0._dp
     call calculate_commutator_nl_terms(qs_env, nlrv, nlrxrv, nlrrv, rcc)

     ! nuclear contribution
     call calculate_localized_nuc_moments(moments, particles, qs_kinds, cell, first_atom, &
          last_atom, rcc)

     ! temp matrix
     allocate(tmp)
     call dbcsr_desymmetrize(matrix_s(1), tmp)

     ! charge
     call dbcsr_dot(densmat, matrix_s(1), trace)
     moments(2) = - focc * trace
     moments(3) = moments(1) + moments(2)

     ! edip
     do i = 1, 3
        call dbcsr_multiply("T", "N", 1._dp, densmat, o(i), 0._dp, tmp)
        call dbcsr_trace(tmp, trace)
        moments(i+3) = moments(i+3) - focc * trace
     end do
     ! vel
     do i = 1, 3
        call dbcsr_multiply("T", "N", 1._dp, densmat_im, momo(i), 0._dp, tmp)
        call dbcsr_trace(tmp, trace)
        moments(i+15) = moments(i+15) + focc * trace
     end do
     ! nl
     moments(28:30) = moments(16:18) + nlrv

     ! mdip
     do i = 1, 3
        call dbcsr_multiply("T", "N", 1._dp, densmat_im, mago(i), 0._dp, tmp)
        call dbcsr_trace(tmp, trace)
        moments(i+12) = moments(i+12) + focc * trace
     end do
     ! nl
     moments(25:27) = moments(13:15) + nlrxrv

     ! equp
     do i = 4, 9
        call dbcsr_multiply("T", "N", 1._dp, densmat, o(i), 0._dp, tmp)
        call dbcsr_trace(tmp, trace)
        moments(i+3) = moments(i+3) - focc * trace
     end do
     ! vel
     qupoleder = 0._dp
     do i = 1, 3
        do j = 1, 3
           call dbcsr_multiply("T", "N", 1._dp, densmat, oder(i, j), 0._dp, tmp)
           call dbcsr_trace(tmp, trace)
           qupoleder((i-1)*3+j) = qupoleder((i-1)*3+j) + focc * trace
           call dbcsr_multiply("T", "N", 1._dp, densmat_im, oder(i, j), 0._dp, tmp)
           call dbcsr_trace(tmp, trace)
           qupoleder((i-1)*3+j) = qupoleder((i-1)*3+j) + focc * trace
        end do
     end do
     moments(19) = -2*qupoleder(1) - moments(2)
     moments(20) = -qupoleder(2) - qupoleder(4)
     moments(21) = -qupoleder(3) - qupoleder(7)
     moments(22) = -2*qupoleder(5) - moments(2)
     moments(23) = -qupoleder(6) - qupoleder(8)
     moments(24) = -2*qupoleder(9) - moments(2)
     ! nl
     moments(31:36) = moments(19:24) + nlrrv

     do i = 1, size(o)
        call dbcsr_release(o(i))
        deallocate(o(i))
     end do
     deallocate(o)
     do i = 1, size(oder, 1)
        do j = 1, size(oder, 2)
           call dbcsr_release(oder(i, j))
           deallocate(oder(i, j))
        end do
     end do
     deallocate(oder)
     do i = 1, size(momo)
        call dbcsr_release(momo(i))
        deallocate(momo(i))
     end do
     deallocate(momo)
     do i = 1, size(mago)
        call dbcsr_release(mago(i))
        deallocate(mago(i))
     end do
     deallocate(mago)
     call dbcsr_release(densmat)
     deallocate(densmat)
     call dbcsr_release(densmat_im)
     deallocate(densmat_im)
     call dbcsr_release(tmp)
     deallocate(tmp)
     deallocate(nlrv, nlrrv, nlrxrv)

   end subroutine calculate_localized_kg_moments

! **************************************************************************************************
!> \brief ...
!> \param logger ...
!> \param input ...
!> \param moments ...
!> \param rcc ...
!> \param molname ...
!> \param suffix ...
! **************************************************************************************************
   SUBROUTINE print_localized_moments(logger, input, moments, rcc, molname, suffix)

      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: input
      REAL(dp)                                           :: moments(:), rcc(:)
      CHARACTER(*)                                       :: molname
      CHARACTER(*), OPTIONAL                             :: suffix

      CHARACTER(len=default_string_length)               :: iter
      INTEGER                                            :: output_unit

      IF (PRESENT(suffix)) THEN
         output_unit = cp_print_key_unit_nr(logger=logger, basis_section=input, &
                                            print_key_path="DFT%PRINT%MOMENTS", extension=".csv", &
                                            middle_name="locmoments_"//TRIM(molname)//"_"//TRIM(suffix), log_filename=.FALSE.)
      ELSE
         output_unit = cp_print_key_unit_nr(logger=logger, basis_section=input, &
                                            print_key_path="DFT%PRINT%MOMENTS", extension=".csv", &
                                            middle_name="locmoments_"//TRIM(molname), log_filename=.FALSE.)
      END IF
      iter = cp_iter_string(logger%iter_info)

      IF (output_unit > 0) THEN
         IF (iteration == 1) THEN
            WRITE (output_unit, "(A)") header
         END IF
         WRITE (output_unit, "(I6,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,"// &
                "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,"// &
                "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,"// &
                "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,"// &
                "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3)") &
            iteration, ", "//TRIM(iter)//", ", &
            moments(1), ", ", moments(2), ", ", moments(3), ", ", &
            moments(4), ", ", moments(5), ", ", moments(6), ", ", moments(7), ", ", &
            moments(8), ", ", moments(9), ", ", moments(10), ", ", moments(11), ", ", &
            moments(12), ", ", moments(13), ", ", moments(14), ", ", moments(15), ", ", &
            moments(16), ", ", moments(17), ", ", moments(18), ", ", moments(19), ", ", &
            moments(20), ", ", moments(21), ", ", moments(22), ", ", moments(23), ", ", &
            moments(24), ", ", moments(25), ", ", moments(26), ", ", moments(27), ", ", &
            moments(28), ", ", moments(29), ", ", moments(30), ", ", moments(31), ", ", &
            moments(32), ", ", moments(33), ", ", moments(34), ", ", moments(35), ", ", &
            moments(36), ", ", rcc(1), ", ", rcc(2), ",", rcc(3)
      END IF

      CALL cp_print_key_finished_output(output_unit, logger, input, "DFT%PRINT%MOMENTS")

   END SUBROUTINE print_localized_moments

END MODULE localized_moments

