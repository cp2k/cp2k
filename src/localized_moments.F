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

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE bibliography,                    ONLY: Schreder2021,&
                                              Schreder2024a,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_element,&
                                              cp_cfm_release,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_get_readonly_block_p, dbcsr_init_p, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_readonly_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_set, dbcsr_type, &
        dbcsr_type_antisymmetric, dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
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
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_constants,                 ONLY: use_mom_ref_com
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                              molecule_kind_type
   USE molecule_types,                  ONLY: molecule_type
   USE moments_utils,                   ONLY: get_reference_point
   USE orbital_pointers,                ONLY: current_maxl,&
                                              indco
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE kg_environment_types,            ONLY: kg_environment_type
   USE qs_loc_types,                    ONLY: qs_loc_env_type
   USE qs_moments,                      ONLY: build_local_magmom_matrix,&
                                              build_local_moment_matrix,&
                                              build_local_moments_der_matrix,&
                                              calculate_commutator_nl_terms,&
                                              print_moments,&
                                              print_moments_nl,&
                                              set_label
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** public subroutines ***

   PUBLIC :: calculate_kg_moments,&
             calculate_localized_moments

! *** body ***

CONTAINS

! **************************************************************************************************
!> \brief Calculates localized multipole moments per molecule in the Wannier MO basis
!> \param qs_env ...
!> \param qs_loc_env ...
!> \param mo_local ...
!> \param max_moment ...
!> \param magnetic ...
!> \param vel_reprs ...
!> \param com_nl ...
!> \param loc_print_section ...
! **************************************************************************************************
   SUBROUTINE calculate_localized_moments(qs_env, qs_loc_env, mo_local, max_moment, magnetic, &
                                          vel_reprs, com_nl, loc_print_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_loc_env_type), INTENT(IN)                  :: qs_loc_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_local
      INTEGER, INTENT(IN)                                :: max_moment
      LOGICAL, INTENT(IN)                                :: magnetic, vel_reprs, com_nl
      TYPE(section_vals_type), POINTER                   :: loc_print_section

      CHARACTER(LEN=*), PARAMETER :: routineN = "calculate_localized_moments"

      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)        :: rlab
      CHARACTER(LEN=default_string_length)               :: mol_name
      INTEGER :: akind, first_atom, handle, i, iatom, idir, imo_im, imo_re, imol, iounit, iproc, &
                 ispin, ix, iy, iz, l, last_atom, nao, natom, nm, nmom, nmols, nspins, ns, order, &
                 rmom_vel_size
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: states
      INTEGER, DIMENSION(2)                              :: nstates
      LOGICAL                                            :: do_rtp, floating, ghost
      REAL(dp)                                           :: charge, dd, im_part, re_part, zwfc
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: mmom, nlcom_rrv, nlcom_rrv_vrr, nlcom_rv, nlcom_rvr, &
                                             nlcom_rxrv, qupole_der, rmom_vel
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rmom
      REAL(dp), DIMENSION(3)                             :: rcc, ria
      REAL(dp), DIMENSION(:), POINTER                    :: ref_point
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_type)                                  :: ov_cfm, v_cfm, vov_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: o_struct, v_struct, vov_struct
      TYPE(cp_fm_type)                                   :: or_fm, vi_fm, vr_fm
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magmom, matrix_s, moments, momentum
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: moments_der
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)

      CALL timeset(routineN, handle)

      CALL cite_reference(Schreder2021)
      IF (max_moment > 1 .OR. magnetic .OR. vel_reprs .OR. com_nl) THEN
         CALL cite_reference(Schreder2024a)
      END IF

      logger => cp_get_default_logger()

      NULLIFY (cell, dft_control, matrix_s, molecule_set, qs_kind_set, sab_orb)
      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      molecule_set=molecule_set, &
                      qs_kind_set=qs_kind_set, &
                      sab_orb=sab_orb)
      particle_set => qs_loc_env%particle_set
      para_env => qs_loc_env%para_env

      nspins = dft_control%nspins
      zwfc = 3.0_dp - REAL(nspins, KIND=dp)
      ! mo_local is sized nspins for SCF and 2*nspins for RTP (real, imag interleaved per spin)
      do_rtp = qs_env%run_rtp .OR. (SIZE(mo_local) == 2*nspins)

      nmom = MIN(max_moment, current_maxl)
      ! number of independent multipole components for orders 1..nmom
      nm = (6 + 11*nmom + 6*nmom**2 + nmom**3)/6 - 1
      nmols = SIZE(molecule_set)

      iounit = cp_print_key_unit_nr(logger, loc_print_section, "LOCALIZED_MOMENTS", &
                                    extension=".LocMom", middle_name="LOCALIZED_MOMENTS")

      ALLOCATE (rmom(nm + 1, 3), rlab(nm + 1))
      IF (magnetic) ALLOCATE (mmom(3))
      IF (vel_reprs) THEN
         IF (com_nl .AND. (nmom >= 2)) THEN
            rmom_vel_size = 21
         ELSE
            rmom_vel_size = MAX(nm, 9)
         END IF
         ALLOCATE (rmom_vel(rmom_vel_size))
      END IF

      DO imol = 1, nmols
         molecule_kind => molecule_set(imol)%molecule_kind
         first_atom = molecule_set(imol)%first_atom
         last_atom = molecule_set(imol)%last_atom
         CALL get_molecule_kind(molecule_kind=molecule_kind, natom=natom, name=mol_name)
         NULLIFY (ref_point)
         CALL get_reference_point(rcc, qs_env=qs_env, reference=use_mom_ref_com, &
                                  ref_point=ref_point, ifirst=first_atom, &
                                  ilast=last_atom)

         rmom = 0.0_dp
         rlab = ""

         NULLIFY (moments)
         CALL dbcsr_allocate_matrix_set(moments, nm)
         DO i = 1, nm
            ALLOCATE (moments(i)%matrix)
            IF (vel_reprs .AND. (nmom >= 2)) THEN
               CALL dbcsr_create(moments(i)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_symmetric)
               CALL cp_dbcsr_alloc_block_from_nbl(moments(i)%matrix, sab_orb)
            ELSE
               CALL dbcsr_copy(moments(i)%matrix, matrix_s(1)%matrix, "Moments")
            END IF
            CALL dbcsr_set(moments(i)%matrix, 0.0_dp)
         END DO
         IF (vel_reprs .AND. (nmom >= 2)) THEN
            NULLIFY (moments_der)
            CALL dbcsr_allocate_matrix_set(moments_der, 3, 3)
            DO i = 1, 3
               DO idir = 1, 3
                  CALL dbcsr_init_p(moments_der(i, idir)%matrix)
                  CALL dbcsr_create(moments_der(i, idir)%matrix, template=matrix_s(1)%matrix, &
                                    matrix_type=dbcsr_type_antisymmetric)
                  CALL cp_dbcsr_alloc_block_from_nbl(moments_der(i, idir)%matrix, sab_orb)
                  CALL dbcsr_set(moments_der(i, idir)%matrix, 0.0_dp)
               END DO
            END DO
            CALL build_local_moments_der_matrix(qs_env, moments_der, 1, 2, ref_point=rcc, &
                                                moments=moments)
         ELSE
            CALL build_local_moment_matrix(qs_env, moments, nmom, ref_point=rcc)
         END IF

         IF (magnetic) THEN
            NULLIFY (magmom)
            CALL dbcsr_allocate_matrix_set(magmom, 3)
            DO i = 1, 3
               CALL dbcsr_init_p(magmom(i)%matrix)
               CALL dbcsr_create(magmom(i)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               CALL cp_dbcsr_alloc_block_from_nbl(magmom(i)%matrix, sab_orb)
               CALL dbcsr_set(magmom(i)%matrix, 0.0_dp)
            END DO
            CALL build_local_magmom_matrix(qs_env, magmom, nmom, ref_point=rcc)
            mmom = 0.0_dp
         END IF

         IF (vel_reprs) THEN
            NULLIFY (momentum)
            CALL dbcsr_allocate_matrix_set(momentum, 3)
            DO i = 1, 3
               CALL dbcsr_init_p(momentum(i)%matrix)
               CALL dbcsr_create(momentum(i)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               CALL cp_dbcsr_alloc_block_from_nbl(momentum(i)%matrix, sab_orb)
               CALL dbcsr_set(momentum(i)%matrix, 0.0_dp)
            END DO
            CALL build_lin_mom_matrix(qs_env, momentum)
            rmom_vel = 0.0_dp
         END IF

         CALL cp_fm_get_info(mo_local(1), nrow_global=nao)
         CALL cp_fm_struct_create(o_struct, nrow_global=nao, ncol_global=nao, &
                                  template_fmstruct=mo_local(1)%matrix_struct)
         CALL cp_fm_create(or_fm, o_struct)

         DO ispin = 1, nspins
            IF (ASSOCIATED(molecule_set(imol)%lmi)) THEN
               nstates(1) = molecule_set(imol)%lmi(ispin)%nstates
            ELSE
               nstates(1) = 0
            END IF
            nstates(2) = para_env%mepos
            CALL para_env%maxloc(nstates)
            IF (nstates(1) == 0) CYCLE
            ns = nstates(1)
            iproc = nstates(2)
            ALLOCATE (states(ns))
            IF (iproc == para_env%mepos) THEN
               states(:) = molecule_set(imol)%lmi(ispin)%states(:)
            ELSE
               states(:) = 0
            END IF
            CALL para_env%bcast(states, iproc)

            IF (do_rtp) THEN
               imo_re = 2*ispin - 1
               imo_im = 2*ispin
            ELSE
               imo_re = ispin
               imo_im = 0
            END IF

            CALL cp_fm_struct_create(v_struct, nrow_global=nao, ncol_global=ns, &
                                     template_fmstruct=mo_local(imo_re)%matrix_struct)
            CALL cp_fm_struct_create(vov_struct, nrow_global=ns, ncol_global=ns, &
                                     template_fmstruct=mo_local(imo_re)%matrix_struct)
            CALL cp_fm_create(vr_fm, v_struct, name="vr")
            CALL cp_fm_set_all(vr_fm, 0.0_dp)
            IF (do_rtp) THEN
               CALL cp_fm_create(vi_fm, v_struct, name="vi")
               CALL cp_fm_set_all(vi_fm, 0.0_dp)
            END IF
            CALL cp_cfm_create(v_cfm, v_struct, name="v")
            CALL cp_cfm_create(ov_cfm, v_struct, name="ov")
            CALL cp_cfm_create(vov_cfm, vov_struct, name="vov")
            CALL cp_fm_struct_release(v_struct)
            CALL cp_fm_struct_release(vov_struct)

            DO i = 1, ns
               CALL cp_fm_to_fm(mo_local(imo_re), vr_fm, 1, states(i), i)
               IF (do_rtp) THEN
                  CALL cp_fm_to_fm(mo_local(imo_im), vi_fm, 1, states(i), i)
               END IF
            END DO
            IF (do_rtp) THEN
               CALL cp_fm_to_cfm(vr_fm, vi_fm, v_cfm)
               CALL cp_fm_release(vi_fm)
            ELSE
               CALL cp_fm_to_cfm(vr_fm, mtarget=v_cfm)
            END IF
            CALL cp_fm_release(vr_fm)

            DO i = 1, nm
               CALL trace_op_in_mo(moments(i)%matrix, or_fm, v_cfm, ov_cfm, vov_cfm, &
                                   o_struct, ns, nao, re_part, im_part)
               rmom(i + 1, 1) = rmom(i + 1, 1) + zwfc*re_part
            END DO
            rmom(1, 1) = rmom(1, 1) + zwfc*REAL(ns, KIND=dp)

            IF (magnetic) THEN
               DO i = 1, 3
                  CALL trace_op_in_mo(magmom(i)%matrix, or_fm, v_cfm, ov_cfm, vov_cfm, &
                                      o_struct, ns, nao, re_part, im_part)
                  mmom(i) = mmom(i) - zwfc*im_part
               END DO
            END IF

            IF (vel_reprs) THEN
               DO order = 1, nmom
                  SELECT CASE (order)
                  CASE (1)
                     DO i = 1, 3
                        CALL trace_op_in_mo(momentum(i)%matrix, or_fm, v_cfm, ov_cfm, vov_cfm, &
                                            o_struct, ns, nao, re_part, im_part)
                        rmom_vel(i) = rmom_vel(i) - zwfc*im_part
                     END DO
                  CASE (2)
                     ALLOCATE (qupole_der(9))
                     qupole_der = 0.0_dp
                     DO i = 1, 3
                        DO idir = 1, 3
                           CALL trace_op_in_mo(moments_der(i, idir)%matrix, or_fm, v_cfm, ov_cfm, &
                                               vov_cfm, o_struct, ns, nao, re_part, im_part)
                           qupole_der((i - 1)*3 + idir) = &
                              qupole_der((i - 1)*3 + idir) + zwfc*(re_part - im_part)
                        END DO
                     END DO
                     rmom_vel(4) = -2.0_dp*qupole_der(1) - rmom(1, 1)
                     rmom_vel(5) = -qupole_der(2) - qupole_der(4)
                     rmom_vel(6) = -qupole_der(3) - qupole_der(7)
                     rmom_vel(7) = -2.0_dp*qupole_der(5) - rmom(1, 1)
                     rmom_vel(8) = -qupole_der(6) - qupole_der(8)
                     rmom_vel(9) = -2.0_dp*qupole_der(9) - rmom(1, 1)
                     DEALLOCATE (qupole_der)
                  END SELECT
               END DO
            END IF

            CALL cp_cfm_release(v_cfm)
            CALL cp_cfm_release(ov_cfm)
            CALL cp_cfm_release(vov_cfm)
            DEALLOCATE (states)
         END DO

         CALL cp_fm_release(or_fm)
         CALL cp_fm_struct_release(o_struct)

         CALL para_env%sum(rmom(:, 1))
         IF (magnetic) CALL para_env%sum(mmom)
         IF (vel_reprs) CALL para_env%sum(rmom_vel)

         DO iatom = first_atom, first_atom + natom - 1
            ria = pbc(particle_set(iatom)%r - rcc, cell) + rcc
            ria = ria - rcc
            atomic_kind => particle_set(iatom)%atomic_kind
            CALL get_atomic_kind(atomic_kind, kind_number=akind)
            CALL get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating, &
                             core_charge=charge)
            IF (ghost .OR. floating) CYCLE
            rmom(1, 2) = rmom(1, 2) - charge
            DO l = 1, nm
               ix = indco(1, l + 1)
               iy = indco(2, l + 1)
               iz = indco(3, l + 1)
               dd = 1.0_dp
               IF (ix > 0) dd = dd*ria(1)**ix
               IF (iy > 0) dd = dd*ria(2)**iy
               IF (iz > 0) dd = dd*ria(3)**iz
               rmom(l + 1, 2) = rmom(l + 1, 2) - charge*dd
               CALL set_label(rlab(l + 1), ix, iy, iz)
            END DO
         END DO
         rmom(:, :) = -rmom(:, :)
         rmom(:, 3) = rmom(:, 1) + rmom(:, 2)

         IF (com_nl) THEN
            IF ((nmom >= 1) .AND. vel_reprs) THEN
               ALLOCATE (nlcom_rv(3))
               nlcom_rv(:) = 0.0_dp
            END IF
            IF ((nmom >= 2) .AND. vel_reprs) THEN
               ALLOCATE (nlcom_rrv(6), nlcom_rvr(6), nlcom_rrv_vrr(6))
               nlcom_rrv(:) = 0.0_dp
               nlcom_rvr(:) = 0.0_dp
               nlcom_rrv_vrr(:) = 0.0_dp
            END IF
            IF (magnetic) THEN
               ALLOCATE (nlcom_rxrv(3))
               nlcom_rxrv(:) = 0.0_dp
            END IF
            CALL calculate_commutator_nl_terms(qs_env, nlcom_rv, nlcom_rxrv, nlcom_rrv, &
                                               nlcom_rvr, nlcom_rrv_vrr, rcc)
         END IF

         IF (iounit > 0) THEN
            WRITE (UNIT=iounit, FMT='(/,T3,A,I6,2X,A)') &
               "# Molecule:", imol, TRIM(mol_name)
            IF (magnetic .AND. vel_reprs) THEN
               CALL print_moments(iounit, nmom, rmom, rlab, rcc, cell, periodic=.FALSE., &
                                  mmom=mmom, rmom_vel=rmom_vel)
            ELSE IF (magnetic) THEN
               CALL print_moments(iounit, nmom, rmom, rlab, rcc, cell, periodic=.FALSE., &
                                  mmom=mmom)
            ELSE IF (vel_reprs) THEN
               CALL print_moments(iounit, nmom, rmom, rlab, rcc, cell, periodic=.FALSE., &
                                  rmom_vel=rmom_vel)
            ELSE
               CALL print_moments(iounit, nmom, rmom, rlab, rcc, cell, periodic=.FALSE.)
            END IF

            IF (com_nl) THEN
               IF (magnetic) mmom(:) = nlcom_rxrv(:)
               IF (vel_reprs .AND. (nmom >= 1)) rmom_vel(1:3) = nlcom_rv
               IF (vel_reprs .AND. (nmom >= 2)) THEN
                  rmom_vel(4:9) = nlcom_rrv
                  rmom_vel(10:15) = nlcom_rvr
                  rmom_vel(16:21) = nlcom_rrv_vrr
               END IF
               IF (magnetic .AND. vel_reprs) THEN
                  CALL print_moments_nl(iounit, nmom, rlab, mmom=mmom, rmom_vel=rmom_vel)
               ELSE IF (magnetic) THEN
                  CALL print_moments_nl(iounit, nmom, rlab, mmom=mmom)
               ELSE IF (vel_reprs) THEN
                  CALL print_moments_nl(iounit, nmom, rlab, rmom_vel=rmom_vel)
               END IF
            END IF
         END IF

         IF (com_nl) THEN
            IF ((nmom >= 1) .AND. vel_reprs) DEALLOCATE (nlcom_rv)
            IF ((nmom >= 2) .AND. vel_reprs) DEALLOCATE (nlcom_rrv, nlcom_rvr, nlcom_rrv_vrr)
            IF (magnetic) DEALLOCATE (nlcom_rxrv)
         END IF
         CALL dbcsr_deallocate_matrix_set(moments)
         IF (vel_reprs .AND. (nmom >= 2)) CALL dbcsr_deallocate_matrix_set(moments_der)
         IF (magnetic) CALL dbcsr_deallocate_matrix_set(magmom)
         IF (vel_reprs) CALL dbcsr_deallocate_matrix_set(momentum)
      END DO

      DEALLOCATE (rmom, rlab)
      IF (magnetic) DEALLOCATE (mmom)
      IF (vel_reprs) DEALLOCATE (rmom_vel)

      CALL cp_print_key_finished_output(iounit, logger, loc_print_section, "LOCALIZED_MOMENTS")

      CALL timestop(handle)

   END SUBROUTINE calculate_localized_moments

! **************************************************************************************************
!> \brief Trace a real dbcsr operator against a complex MO matrix
!> \param op_dbcsr ...
!> \param or_fm ...
!> \param v_cfm ...
!> \param ov_cfm ...
!> \param vov_cfm ...
!> \param o_struct ...
!> \param ns ...
!> \param nao ...
!> \param re_sum ...
!> \param im_sum ...
! **************************************************************************************************
   SUBROUTINE trace_op_in_mo(op_dbcsr, or_fm, v_cfm, ov_cfm, vov_cfm, o_struct, ns, nao, &
                             re_sum, im_sum)
      TYPE(dbcsr_type), POINTER                          :: op_dbcsr
      TYPE(cp_fm_type), INTENT(INOUT)                    :: or_fm
      TYPE(cp_cfm_type), INTENT(IN)                      :: v_cfm
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: ov_cfm, vov_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: o_struct
      INTEGER, INTENT(IN)                                :: ns, nao
      REAL(dp), INTENT(OUT)                              :: re_sum, im_sum

      COMPLEX(dp)                                        :: di
      INTEGER                                            :: j
      TYPE(cp_cfm_type)                                  :: o_cfm

      CALL copy_dbcsr_to_fm(op_dbcsr, or_fm)
      CALL cp_cfm_create(o_cfm, o_struct)
      CALL cp_fm_to_cfm(or_fm, mtarget=o_cfm)
      CALL cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o_cfm, v_cfm, &
                       (0._dp, 0._dp), ov_cfm)
      CALL cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), v_cfm, ov_cfm, &
                       (0._dp, 0._dp), vov_cfm)
      re_sum = 0.0_dp
      im_sum = 0.0_dp
      DO j = 1, ns
         CALL cp_cfm_get_element(vov_cfm, j, j, di)
         re_sum = re_sum + REAL(di, KIND=dp)
         im_sum = im_sum + AIMAG(di)
      END DO
      CALL cp_cfm_release(o_cfm)

   END SUBROUTINE trace_op_in_mo
   
! **************************************************************************************************
!> \brief Calculates multipole moments per molecule from the Kim-Gordon AO density matrix
!> \param qs_env ...
!> \param unit_nr  output unit (the open DFT%PRINT%MOMENTS unit from the caller)
!> \param max_moment ...
!> \param magnetic ...
!> \param vel_reprs ...
!> \param com_nl ...
! **************************************************************************************************
   SUBROUTINE calculate_kg_moments(qs_env, unit_nr, max_moment, magnetic, vel_reprs, com_nl)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr, max_moment
      LOGICAL, INTENT(IN)                                :: magnetic, vel_reprs, com_nl

      CHARACTER(LEN=*), PARAMETER :: routineN = "calculate_kg_moments"

      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)        :: rlab
      CHARACTER(LEN=default_string_length)               :: mol_name
      INTEGER :: akind, first_atom, handle, i, iatom, idir, imol, ispin, ix, iy, iz, &
                 jatom, l, last_atom, natom, nm, nmol, nmom, nspins, rmom_vel_size
      INTEGER, DIMENSION(:), POINTER                     :: atom_to_molecule
      LOGICAL                                            :: do_rtp, floating, found, ghost
      REAL(dp)                                           :: charge, dd, factor
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: mmom, nlcom_rrv, nlcom_rrv_vrr, nlcom_rv, &
                                             nlcom_rvr, nlcom_rxrv, qupole_der, rmom_vel
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rmom
      REAL(dp), DIMENSION(3)                             :: rcc, ria
      REAL(dp), DIMENSION(:), POINTER                    :: ref_point
      REAL(dp), DIMENSION(:, :), POINTER                 :: oblock, pblock, sblock
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magmom, matrix_s, moments, &
                                                            momentum, rho_ao, rho_ao_im
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: moments_der
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      CALL cite_reference(Schreder2021)
      IF (max_moment > 1 .OR. magnetic .OR. vel_reprs .OR. com_nl) THEN
         CALL cite_reference(Schreder2024a)
      END IF

      NULLIFY (cell, dft_control, kg_env, matrix_s, para_env, particle_set, &
               qs_kind_set, rho, sab_orb)
      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      kg_env=kg_env, &
                      matrix_s=matrix_s, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      rho=rho, &
                      sab_orb=sab_orb)

      CPASSERT(ASSOCIATED(kg_env))
      molecule_set => kg_env%molecule_set
      atom_to_molecule => kg_env%atom_to_molecule

      nspins = dft_control%nspins
      do_rtp = qs_env%run_rtp

      nmom = MIN(max_moment, current_maxl)
      nm = (6 + 11*nmom + 6*nmom**2 + nmom**3)/6 - 1
      nmol = SIZE(molecule_set)

      NULLIFY (rho_ao, rho_ao_im)
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      IF (do_rtp) CALL qs_rho_get(rho, rho_ao_im=rho_ao_im)

      ALLOCATE (rmom(nm + 1, 3), rlab(nm + 1))
      IF (magnetic) ALLOCATE (mmom(3))
      IF (vel_reprs) THEN
         IF (com_nl .AND. (nmom >= 2)) THEN
            rmom_vel_size = 21
         ELSE
            rmom_vel_size = MAX(nm, 9)
         END IF
         ALLOCATE (rmom_vel(rmom_vel_size))
      END IF

      DO imol = 1, nmol
         molecule_kind => molecule_set(imol)%molecule_kind
         first_atom = molecule_set(imol)%first_atom
         last_atom = molecule_set(imol)%last_atom
         CALL get_molecule_kind(molecule_kind=molecule_kind, natom=natom, name=mol_name)
         NULLIFY (ref_point)
         CALL get_reference_point(rcc, qs_env=qs_env, reference=use_mom_ref_com, &
                                  ref_point=ref_point, ifirst=first_atom, ilast=last_atom)

         rmom = 0.0_dp
         rlab = ""

         NULLIFY (moments)
         CALL dbcsr_allocate_matrix_set(moments, nm)
         DO i = 1, nm
            ALLOCATE (moments(i)%matrix)
            IF (vel_reprs .AND. (nmom >= 2)) THEN
               CALL dbcsr_create(moments(i)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_symmetric)
               CALL cp_dbcsr_alloc_block_from_nbl(moments(i)%matrix, sab_orb)
            ELSE
               CALL dbcsr_copy(moments(i)%matrix, matrix_s(1)%matrix, "Moments")
            END IF
            CALL dbcsr_set(moments(i)%matrix, 0.0_dp)
         END DO
         IF (vel_reprs .AND. (nmom >= 2)) THEN
            NULLIFY (moments_der)
            CALL dbcsr_allocate_matrix_set(moments_der, 3, 3)
            DO i = 1, 3
               DO idir = 1, 3
                  CALL dbcsr_init_p(moments_der(i, idir)%matrix)
                  CALL dbcsr_create(moments_der(i, idir)%matrix, template=matrix_s(1)%matrix, &
                                    matrix_type=dbcsr_type_antisymmetric)
                  CALL cp_dbcsr_alloc_block_from_nbl(moments_der(i, idir)%matrix, sab_orb)
                  CALL dbcsr_set(moments_der(i, idir)%matrix, 0.0_dp)
               END DO
            END DO
            CALL build_local_moments_der_matrix(qs_env, moments_der, 1, 2, ref_point=rcc, &
                                                moments=moments)
         ELSE
            CALL build_local_moment_matrix(qs_env, moments, nmom, ref_point=rcc)
         END IF

         IF (magnetic) THEN
            NULLIFY (magmom)
            CALL dbcsr_allocate_matrix_set(magmom, 3)
            DO i = 1, 3
               CALL dbcsr_init_p(magmom(i)%matrix)
               CALL dbcsr_create(magmom(i)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               CALL cp_dbcsr_alloc_block_from_nbl(magmom(i)%matrix, sab_orb)
               CALL dbcsr_set(magmom(i)%matrix, 0.0_dp)
            END DO
            CALL build_local_magmom_matrix(qs_env, magmom, nmom, ref_point=rcc)
            mmom = 0.0_dp
         END IF

         IF (vel_reprs) THEN
            NULLIFY (momentum)
            CALL dbcsr_allocate_matrix_set(momentum, 3)
            DO i = 1, 3
               CALL dbcsr_init_p(momentum(i)%matrix)
               CALL dbcsr_create(momentum(i)%matrix, template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               CALL cp_dbcsr_alloc_block_from_nbl(momentum(i)%matrix, sab_orb)
               CALL dbcsr_set(momentum(i)%matrix, 0.0_dp)
            END DO
            CALL build_lin_mom_matrix(qs_env, momentum)
            rmom_vel = 0.0_dp
         END IF

         IF (vel_reprs .AND. (nmom >= 2)) THEN
            ALLOCATE (qupole_der(9))
            qupole_der = 0.0_dp
         END IF

         DO ispin = 1, nspins
            CALL dbcsr_iterator_readonly_start(iter, rho_ao(ispin)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iatom, jatom, pblock)
               IF (SIZE(pblock) == 0) CYCLE
               IF (atom_to_molecule(iatom) /= imol) CYCLE
               IF (atom_to_molecule(jatom) /= imol) CYCLE
               factor = MERGE(2.0_dp, 1.0_dp, iatom /= jatom)
               CALL dbcsr_get_readonly_block_p(matrix_s(1)%matrix, iatom, jatom, sblock, found)
               IF (found) rmom(1, 1) = rmom(1, 1) + factor*SUM(pblock*sblock)
               DO i = 1, nm
                  CALL dbcsr_get_readonly_block_p(moments(i)%matrix, iatom, jatom, oblock, found)
                  IF (found) rmom(i + 1, 1) = rmom(i + 1, 1) + factor*SUM(pblock*oblock)
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)

            IF (do_rtp .AND. (magnetic .OR. vel_reprs)) THEN
               CALL dbcsr_iterator_readonly_start(iter, rho_ao_im(ispin)%matrix)
               DO WHILE (dbcsr_iterator_blocks_left(iter))
                  CALL dbcsr_iterator_next_block(iter, iatom, jatom, pblock)
                  IF (SIZE(pblock) == 0) CYCLE
                  IF (atom_to_molecule(iatom) /= imol) CYCLE
                  IF (atom_to_molecule(jatom) /= imol) CYCLE
                  factor = MERGE(2.0_dp, 1.0_dp, iatom /= jatom)
                  IF (magnetic) THEN
                     DO i = 1, 3
                        CALL dbcsr_get_readonly_block_p(magmom(i)%matrix, iatom, jatom, oblock, found)
                        IF (found) mmom(i) = mmom(i) + factor*SUM(pblock*oblock)
                     END DO
                  END IF
                  IF (vel_reprs) THEN
                     DO i = 1, 3
                        CALL dbcsr_get_readonly_block_p(momentum(i)%matrix, iatom, jatom, oblock, found)
                        IF (found) rmom_vel(i) = rmom_vel(i) + factor*SUM(pblock*oblock)
                     END DO
                     IF (nmom >= 2) THEN
                        DO i = 1, 3
                           DO idir = 1, 3
                              CALL dbcsr_get_readonly_block_p(moments_der(i, idir)%matrix, &
                                                              iatom, jatom, oblock, found)
                              IF (found) &
                                 qupole_der((i - 1)*3 + idir) = &
                                 qupole_der((i - 1)*3 + idir) - factor*SUM(pblock*oblock)
                           END DO
                        END DO
                     END IF
                  END IF
               END DO
               CALL dbcsr_iterator_stop(iter)
            END IF
         END DO

         IF (vel_reprs .AND. (nmom >= 2)) THEN
            rmom_vel(4) = -2.0_dp*qupole_der(1) - rmom(1, 1)
            rmom_vel(5) = -qupole_der(2) - qupole_der(4)
            rmom_vel(6) = -qupole_der(3) - qupole_der(7)
            rmom_vel(7) = -2.0_dp*qupole_der(5) - rmom(1, 1)
            rmom_vel(8) = -qupole_der(6) - qupole_der(8)
            rmom_vel(9) = -2.0_dp*qupole_der(9) - rmom(1, 1)
            DEALLOCATE (qupole_der)
         END IF

         CALL para_env%sum(rmom(:, 1))
         IF (magnetic) CALL para_env%sum(mmom)
         IF (vel_reprs) CALL para_env%sum(rmom_vel)

         DO iatom = first_atom, first_atom + natom - 1
            ria = pbc(particle_set(iatom)%r - rcc, cell) + rcc
            ria = ria - rcc
            atomic_kind => particle_set(iatom)%atomic_kind
            CALL get_atomic_kind(atomic_kind, kind_number=akind)
            CALL get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating, &
                             core_charge=charge)
            IF (ghost .OR. floating) CYCLE
            rmom(1, 2) = rmom(1, 2) - charge
            DO l = 1, nm
               ix = indco(1, l + 1)
               iy = indco(2, l + 1)
               iz = indco(3, l + 1)
               dd = 1.0_dp
               IF (ix > 0) dd = dd*ria(1)**ix
               IF (iy > 0) dd = dd*ria(2)**iy
               IF (iz > 0) dd = dd*ria(3)**iz
               rmom(l + 1, 2) = rmom(l + 1, 2) - charge*dd
               CALL set_label(rlab(l + 1), ix, iy, iz)
            END DO
         END DO
         rmom(:, :) = -rmom(:, :)
         rmom(:, 3) = rmom(:, 1) + rmom(:, 2)

         IF (com_nl) THEN
            IF ((nmom >= 1) .AND. vel_reprs) THEN
               ALLOCATE (nlcom_rv(3))
               nlcom_rv(:) = 0.0_dp
            END IF
            IF ((nmom >= 2) .AND. vel_reprs) THEN
               ALLOCATE (nlcom_rrv(6), nlcom_rvr(6), nlcom_rrv_vrr(6))
               nlcom_rrv(:) = 0.0_dp
               nlcom_rvr(:) = 0.0_dp
               nlcom_rrv_vrr(:) = 0.0_dp
            END IF
            IF (magnetic) THEN
               ALLOCATE (nlcom_rxrv(3))
               nlcom_rxrv(:) = 0.0_dp
            END IF
            CALL calculate_commutator_nl_terms(qs_env, nlcom_rv, nlcom_rxrv, nlcom_rrv, &
                                               nlcom_rvr, nlcom_rrv_vrr, rcc)
         END IF

         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT='(/,T3,A,I6,2X,A)') &
               "# Molecule:", imol, TRIM(mol_name)
            IF (magnetic .AND. vel_reprs) THEN
               CALL print_moments(unit_nr, nmom, rmom, rlab, rcc, cell, periodic=.FALSE., &
                                  mmom=mmom, rmom_vel=rmom_vel)
            ELSE IF (magnetic) THEN
               CALL print_moments(unit_nr, nmom, rmom, rlab, rcc, cell, periodic=.FALSE., &
                                  mmom=mmom)
            ELSE IF (vel_reprs) THEN
               CALL print_moments(unit_nr, nmom, rmom, rlab, rcc, cell, periodic=.FALSE., &
                                  rmom_vel=rmom_vel)
            ELSE
               CALL print_moments(unit_nr, nmom, rmom, rlab, rcc, cell, periodic=.FALSE.)
            END IF

            IF (com_nl) THEN
               IF (magnetic) mmom(:) = nlcom_rxrv(:)
               IF (vel_reprs .AND. (nmom >= 1)) rmom_vel(1:3) = nlcom_rv
               IF (vel_reprs .AND. (nmom >= 2)) THEN
                  rmom_vel(4:9) = nlcom_rrv
                  rmom_vel(10:15) = nlcom_rvr
                  rmom_vel(16:21) = nlcom_rrv_vrr
               END IF
               IF (magnetic .AND. vel_reprs) THEN
                  CALL print_moments_nl(unit_nr, nmom, rlab, mmom=mmom, rmom_vel=rmom_vel)
               ELSE IF (magnetic) THEN
                  CALL print_moments_nl(unit_nr, nmom, rlab, mmom=mmom)
               ELSE IF (vel_reprs) THEN
                  CALL print_moments_nl(unit_nr, nmom, rlab, rmom_vel=rmom_vel)
               END IF
            END IF
         END IF

         IF (com_nl) THEN
            IF ((nmom >= 1) .AND. vel_reprs) DEALLOCATE (nlcom_rv)
            IF ((nmom >= 2) .AND. vel_reprs) DEALLOCATE (nlcom_rrv, nlcom_rvr, nlcom_rrv_vrr)
            IF (magnetic) DEALLOCATE (nlcom_rxrv)
         END IF
         CALL dbcsr_deallocate_matrix_set(moments)
         IF (vel_reprs .AND. (nmom >= 2)) CALL dbcsr_deallocate_matrix_set(moments_der)
         IF (magnetic) CALL dbcsr_deallocate_matrix_set(magmom)
         IF (vel_reprs) CALL dbcsr_deallocate_matrix_set(momentum)
      END DO

      DEALLOCATE (rmom, rlab)
      IF (magnetic) DEALLOCATE (mmom)
      IF (vel_reprs) DEALLOCATE (rmom_vel)

      CALL timestop(handle)

   END SUBROUTINE calculate_kg_moments

END MODULE localized_moments

