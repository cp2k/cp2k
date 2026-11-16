!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate ntrinsic atomic orbitals and analyze wavefunctions
!> \par History
!>      03.2023 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE iao_types
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type,&
                                              cp_fm_write_formatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE input_constants,                 ONLY: do_iaoloc_enone,&
                                              do_iaoloc_pm2
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC ::  iao_env_type, iao_read_input, iao_set_default, &
             organise_printout_for_rose, print_fragment_association

! **************************************************************************************************
   TYPE iao_env_type
      LOGICAL                                :: do_iao = .FALSE.
      !
      REAL(KIND=dp)                          :: eps_svd = 0.0_dp
      REAL(KIND=dp)                          :: eps_occ = 0.0_dp
      ! chages
      LOGICAL                                :: do_charges = .FALSE.
      ! one-center expansion
      LOGICAL                                :: do_oce = .FALSE.
      INTEGER                                :: lmax_oce = 0
      INTEGER                                :: nbas_oce = 0
      LOGICAL                                :: append_oce = .FALSE.
      ! Bond orbitals
      LOGICAL                                :: do_bondorbitals = .FALSE.
      ! Wannier centers
      LOGICAL                                :: do_center = .FALSE.
      LOGICAL                                :: pos_periodic = .FALSE.
      INTEGER                                :: loc_operator = 0
      INTEGER                                :: eloc_function = 0
      REAL(KIND=dp)                          :: eloc_weight = 0.0_dp
      ! Molden
      LOGICAL                                :: molden_iao = .FALSE.
      LOGICAL                                :: molden_ibo = .FALSE.
      ! CUBE files
      LOGICAL                                :: cubes_iao = .FALSE.
      LOGICAL                                :: cubes_ibo = .FALSE.
      ! Input sections
      TYPE(section_vals_type), POINTER       :: iao_cubes_section => NULL(), &
                                                iao_molden_section => NULL(), &
                                                ibo_cubes_section => NULL(), &
                                                ibo_molden_section => NULL(), &
                                                ibo_cc_section => NULL()
      ! Fragments
      LOGICAL                                 :: do_fragments = .FALSE.
      INTEGER                                 :: nfrags = 0
   END TYPE iao_env_type

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param iao_env ...
! **************************************************************************************************
   SUBROUTINE iao_set_default(iao_env)
      TYPE(iao_env_type), INTENT(INOUT)                  :: iao_env

      !iao
      iao_env%do_iao = .FALSE.
      iao_env%eps_svd = 0.0_dp
      iao_env%eps_occ = 0.0_dp
      ! charges
      iao_env%do_charges = .FALSE.
      ! one-center expansion
      iao_env%do_oce = .FALSE.
      iao_env%lmax_oce = 3
      iao_env%nbas_oce = 10
      iao_env%append_oce = .FALSE.
      ! Bond orbitals
      iao_env%do_bondorbitals = .FALSE.
      ! Wannier centers
      iao_env%do_center = .FALSE.
      iao_env%pos_periodic = .FALSE.
      iao_env%loc_operator = do_iaoloc_pm2
      iao_env%eloc_function = do_iaoloc_enone
      iao_env%eloc_weight = 0.0_dp
      ! i/o
      iao_env%molden_iao = .FALSE.
      iao_env%molden_ibo = .FALSE.
      iao_env%cubes_iao = .FALSE.
      iao_env%cubes_ibo = .FALSE.
      ! Input sections
      NULLIFY (iao_env%iao_cubes_section, iao_env%iao_molden_section)
      NULLIFY (iao_env%ibo_cubes_section, iao_env%ibo_molden_section)
      NULLIFY (iao_env%ibo_cc_section)

   END SUBROUTINE iao_set_default

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param iao_env ...
!> \param iao_section ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE iao_read_input(iao_env, iao_section, cell)
      TYPE(iao_env_type), INTENT(INOUT)                  :: iao_env
      TYPE(section_vals_type), POINTER                   :: iao_section
      TYPE(cell_type), OPTIONAL                          :: cell

      LOGICAL                                            :: explicit, iao_explicit
      TYPE(section_vals_type), POINTER                   :: subsection

      CALL iao_set_default(iao_env)

      CALL section_vals_get(iao_section, explicit=iao_explicit)
      IF (iao_explicit) THEN
         iao_env%do_iao = .TRUE.
         ! input options
         CALL section_vals_val_get(iao_section, "EPS_SVD", r_val=iao_env%eps_svd)
         CALL section_vals_val_get(iao_section, "EPS_OCC", r_val=iao_env%eps_occ)
         CALL section_vals_val_get(iao_section, "ATOMIC_CHARGES", l_val=iao_env%do_charges)
         iao_env%iao_molden_section => section_vals_get_subs_vals(iao_section, "IAO_MOLDEN")
         CALL section_vals_get(iao_env%iao_molden_section, explicit=iao_env%molden_iao)
         iao_env%iao_cubes_section => section_vals_get_subs_vals(iao_section, "IAO_CUBES")
         CALL section_vals_get(iao_env%iao_cubes_section, explicit=iao_env%cubes_iao)
         subsection => section_vals_get_subs_vals(iao_section, "ONE_CENTER_EXPANSION")
         CALL section_vals_get(subsection, explicit=iao_env%do_oce)
         IF (iao_env%do_oce) THEN
            subsection => section_vals_get_subs_vals(iao_section, "ONE_CENTER_EXPANSION")
            CALL section_vals_val_get(subsection, "LMAX", i_val=iao_env%lmax_oce)
            CALL section_vals_val_get(subsection, "NBAS", i_val=iao_env%nbas_oce)
            CALL section_vals_val_get(subsection, "APPEND", l_val=iao_env%append_oce)
         END IF
         subsection => section_vals_get_subs_vals(iao_section, "BOND_ORBITALS")
         CALL section_vals_get(subsection, explicit=iao_env%do_bondorbitals)
         IF (iao_env%do_bondorbitals) THEN
            subsection => section_vals_get_subs_vals(iao_section, "BOND_ORBITALS")
            CALL section_vals_val_get(subsection, "LOCALIZATION_OPERATOR", i_val=iao_env%loc_operator)
            CALL section_vals_val_get(subsection, "ENERGY_LOCALIZATION_FUNCTION", i_val=iao_env%eloc_function)
            CALL section_vals_val_get(subsection, "ENERGY_LOCALIZATION_WEIGHT", r_val=iao_env%eloc_weight)
            iao_env%ibo_molden_section => section_vals_get_subs_vals(subsection, "IBO_MOLDEN")
            CALL section_vals_get(iao_env%ibo_molden_section, explicit=iao_env%molden_ibo)
            iao_env%ibo_cubes_section => section_vals_get_subs_vals(subsection, "IBO_CUBES")
            CALL section_vals_get(iao_env%ibo_cubes_section, explicit=iao_env%cubes_ibo)
            iao_env%ibo_cc_section => section_vals_get_subs_vals(subsection, "CHARGE_CENTER")
            CALL section_vals_get(iao_env%ibo_cc_section, explicit=iao_env%do_center)
            IF (iao_env%do_center) THEN
               CALL section_vals_val_get(iao_env%ibo_cc_section, "POSITION_OPERATOR_BERRY", &
                                         l_val=iao_env%pos_periodic, explicit=explicit)
               IF (.NOT. explicit) THEN
                  ! set default according to cell periodicity
                  iao_env%pos_periodic = .TRUE.
                  IF (PRESENT(cell)) THEN
                     IF (ALL(cell%perd == 0)) iao_env%pos_periodic = .FALSE.
                  END IF
               END IF
            END IF
         END IF
      END IF

   END SUBROUTINE iao_read_input
! **************************************************************************************************
!> \brief Prints fragmented overlap matrices for interface with ROSE
!> \param qs_env ...
!> \param orb_basis_set_list ...
! **************************************************************************************************
   SUBROUTINE organise_printout_for_rose(qs_env, orb_basis_set_list)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: orb_basis_set_list

      CHARACTER(len=*), PARAMETER :: routineN = 'organise_printout_for_rose'

      INTEGER                                            :: after, colskip, handle, iatom, ii, &
                                                            ikind, iounit, jj, ncol, nfrags, nrow, &
                                                            rowfrag, rowskip, totalcol, totalrow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      LOGICAL                                            :: lfirstcol, lfirstcoljj
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type)                                   :: fm_sorb, fm_sorb_frags12, fm_sorb_frags22
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: smat
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      !get number of fragments
      NULLIFY (particle_set, atomic_kind_set)
      CALL get_qs_env(qs_env, particle_set=particle_set, atomic_kind_set=atomic_kind_set, &
                      matrix_s_kp=matrix_s)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
      smat => matrix_s(1, 1)%matrix

      nfrags = 0
      DO ii = 1, SIZE(particle_set)
         IF (nfrags < particle_set(ii)%fragment_index) THEN
            nfrags = particle_set(ii)%fragment_index
         END IF
      END DO

      logger => cp_get_default_logger()
      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           qs_env%input, "DFT%PRINT%IAO_ANALYSIS%IAO_OVERLAP"), cp_p_file)) THEN
         iounit = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%IAO_ANALYSIS%IAO_OVERLAP", &
                                       extension=".Log", file_form="FORMATTED", file_action="WRITE", &
                                       file_status="REPLACE")
         CALL section_vals_val_get(qs_env%input, "DFT%PRINT%IAO_ANALYSIS%IAO_OVERLAP%NDIGITS", i_val=after)
         after = MIN(MAX(after, 1), 16)

         DO ii = 1, nfrags
            DO jj = 1, nfrags
               NULLIFY (para_env, blacs_env)
               CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, para_env=para_env)
               totalcol = 0
               totalrow = 0
               rowfrag = 0
               colskip = 0
               rowskip = 0
               lfirstcol = .TRUE.
               lfirstcoljj = .TRUE.

               DO iatom = 1, SIZE(particle_set)
                  ikind = kind_of(iatom)
                  IF (particle_set(iatom)%fragment_index == jj) THEN
                     totalcol = totalcol + orb_basis_set_list(ikind)%gto_basis_set%nsgf
                     lfirstcoljj = .FALSE.
                  ELSE IF (lfirstcoljj) THEN
                     colskip = colskip + orb_basis_set_list(ikind)%gto_basis_set%nsgf
                  END IF
                  IF (particle_set(iatom)%fragment_index == ii) THEN
                     totalrow = totalrow + orb_basis_set_list(ikind)%gto_basis_set%nsgf
                     lfirstcol = .FALSE.
                  ELSE IF (lfirstcol) THEN
                     rowskip = rowskip + orb_basis_set_list(ikind)%gto_basis_set%nsgf
                  END IF
                  rowfrag = rowfrag + orb_basis_set_list(ikind)%gto_basis_set%nsgf
               END DO

               CALL dbcsr_get_info(matrix_s(1, 1)%matrix, nfullrows_total=nrow, nfullcols_total=ncol)
               CALL cp_fm_struct_create(fmstruct=fmstruct, para_env=para_env, context=blacs_env, &
                                        nrow_global=nrow, ncol_global=ncol)
               CALL cp_fm_create(matrix=fm_sorb, matrix_struct=fmstruct, name="Overlap matrix S_B1")
               CALL cp_fm_struct_release(fmstruct=fmstruct)
               CALL copy_dbcsr_to_fm(matrix_s(1, 1)%matrix, fm_sorb)

               CALL cp_fm_struct_create(fmstruct=fmstruct, para_env=para_env, context=blacs_env, &
                                        nrow_global=totalrow, ncol_global=totalcol)
               CALL cp_fm_create(matrix=fm_sorb_frags22, matrix_struct=fmstruct, name="Overlap matrix S_B1^kk' (S22)")
               CALL cp_fm_struct_release(fmstruct=fmstruct)
               CALL cp_fm_to_fm_submat(fm_sorb, fm_sorb_frags22, &
                                       totalrow, totalcol, 1 + rowskip, 1 + colskip, 1, 1)

               CALL cp_fm_struct_create(fmstruct=fmstruct, para_env=para_env, context=blacs_env, &
                                        nrow_global=rowfrag, ncol_global=totalcol)
               CALL cp_fm_create(matrix=fm_sorb_frags12, matrix_struct=fmstruct, name="Overlap matrix S_B1^k' (S12)")
               CALL cp_fm_struct_release(fmstruct=fmstruct)
               CALL cp_fm_to_fm_submat(fm_sorb, fm_sorb_frags12, rowfrag, totalcol, 1, 1 + colskip, 1, 1)

               IF (iounit > 0) WRITE (iounit, *) "FRAGMENTS kk' NR. ", ii, jj
               CALL cp_fm_write_formatted(fm_sorb_frags22, iounit, "Overlap matrix S_B1^kk' (S22)")
               IF (ii /= jj) THEN
                  IF (iounit > 0) WRITE (iounit, *) "FRAGMENTS kk' NR. ", ii, jj
                  CALL cp_fm_write_formatted(fm_sorb_frags12, iounit, "Overlap matrix S_B1^k' (S12)")
               END IF

               CALL cp_fm_release(fm_sorb)
               CALL cp_fm_release(fm_sorb_frags22)
               CALL cp_fm_release(fm_sorb_frags12)

            END DO
         END DO

         CALL cp_print_key_finished_output(iounit, logger, qs_env%input, "DFT%PRINT%IAO_ANALYSIS%IAO_OVERLAP")
      END IF

      CALL timestop(handle)

   END SUBROUTINE organise_printout_for_rose
! **************************************************************************************************
!> \brief Associates IBO centers to fragments
!> \param qs_env ...
!> \param moments ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_fragment_association(qs_env, moments, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: moments
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'print_fragment_association'

      INTEGER                                            :: handle, iatom, is, ispin, jcenter, &
                                                            natom, ncenters, nspin
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: atomcenter
      REAL(KIND=dp)                                      :: conv, dab, mindab
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env=qs_env, cell=cell, particle_set=particle_set, natom=natom)
      conv = cp_unit_from_cp2k(1.0_dp, TRIM("angstrom"))
      nspin = SIZE(moments, 3)
      ncenters = SIZE(moments, 2)
      ALLOCATE (atomcenter(ncenters, nspin))
      atomcenter = 0
      DO ispin = 1, nspin
      DO jcenter = 1, SIZE(moments, 2)
         mindab = 10.0_dp
         DO iatom = 1, natom
            rab(:) = pbc(particle_set(iatom)%r(:), moments(:, jcenter, ispin), cell)
            dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
            IF (dab <= mindab) THEN
               mindab = dab
               atomcenter(jcenter, ispin) = iatom
            END IF
         END DO
      END DO

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, "(/,T2,A,i1)") "Intrinsic Bond Orbitals: Centers associated to fragments for Spin ", &
            ispin
         WRITE (unit_nr, "(T2,A6,T30,A6,T60,A4,T70,A8)") "Center", "Moment", "Atom", "Fragment"
         DO is = 1, ncenters
            WRITE (unit_nr, "(i7,3F15.8,4X,I7,4X,I7)") is, moments(1:3, is, ispin), atomcenter(is, ispin), &
               particle_set(atomcenter(is, ispin))%fragment_index
         END DO
         WRITE (unit_nr, "(/)")
      END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE print_fragment_association

END MODULE iao_types
