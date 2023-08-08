!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate intrinsic atomic orbitals and analyze wavefunctions
!> \par History
!>      03.2023 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE iao_analysis
   USE ai_contraction,                  ONLY: block_add,&
                                              contraction
   USE ai_overlap,                      ONLY: overlap_ab
   USE atomic_charges,                  ONLY: print_atomic_charges
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE auto_basis,                      ONLY: create_oce_basis
   USE basis_set_types,                 ONLY: deallocate_gto_basis_set,&
                                              get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type,&
                                              init_orb_basis_set
   USE bibliography,                    ONLY: Knizia2013,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_geadd,&
                                              cp_fm_invert
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_element,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_distribution_type, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_replicate_all, dbcsr_reserve_diag_blocks, dbcsr_set, dbcsr_trace, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE iao_types,                       ONLY: iao_env_type
   USE input_constants,                 ONLY: do_iaoloc_enone,&
                                              do_iaoloc_pm2
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: twopi
   USE mathlib,                         ONLY: invmat_symm,&
                                              jacobi
   USE message_passing,                 ONLY: mp_comm_type
   USE min_basis_set,                   ONLY: create_minbas_set
   USE molden_utils,                    ONLY: write_mos_molden
   USE orbital_pointers,                ONLY: ncoset
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_loc_methods,                  ONLY: jacobi_rotation_pipek
   USE qs_loc_utils,                    ONLY: compute_berry_operator
   USE qs_localization_methods,         ONLY: initialize_weights
   USE qs_mo_methods,                   ONLY: make_basis_lowdin
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_type,&
                                              set_mo_set
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: setup_neighbor_list
   USE qs_overlap,                      ONLY: build_overlap_matrix_simple
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'iao_analysis'

   PUBLIC ::  iao_wfn_analysis, iao_charges, iao_calculate_dmat

   INTERFACE iao_calculate_dmat
      MODULE PROCEDURE iao_calculate_dmat_diag, & ! diagonal occupation matrix
         iao_calculate_dmat_full    ! full occupation matrix
   END INTERFACE

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param iao_env ...
!> \param unit_nr ...
!> \param c_iao_coef ...
!> \param mos ...
!> \param bond_centers ...
! **************************************************************************************************
   SUBROUTINE iao_wfn_analysis(qs_env, iao_env, unit_nr, c_iao_coef, mos, bond_centers)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(iao_env_type), INTENT(IN)                     :: iao_env
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:), &
         OPTIONAL                                        :: c_iao_coef
      TYPE(mo_set_type), DIMENSION(:), OPTIONAL, POINTER :: mos
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL        :: bond_centers

      CHARACTER(len=*), PARAMETER                        :: routineN = 'iao_wfn_analysis'

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: bname
      INTEGER                                            :: handle, i, iatom, ikind, isgf, ispin, &
                                                            nao, natom, nimages, nkind, no, norb, &
                                                            nref, ns, nsgf, nspin, nvec, nx, sweeps
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, last_sgf, nsgfat
      INTEGER, DIMENSION(2)                              :: nocc
      LOGICAL                                            :: cubes_iao, cubes_ibo, molden_iao, &
                                                            molden_ibo, uniform_occupation
      REAL(KIND=dp)                                      :: trace, zval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mcharge
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: moments
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occupation_numbers
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: smat_kind
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:, :) :: oce_atom
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: p_orb_ref, p_ref_orb, smo
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: c_loc_orb, ciao, cloc, cvec, iao_coef
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: zij_atom
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: sref, sro, sxo
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_type)                                   :: dmat
      TYPE(dbcsr_type), POINTER                          :: smat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: oce_basis_set_list, orb_basis_set_list, &
                                                            ref_basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: ocebasis, orbbasis, refbasis
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: mo_iao, mo_loc
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: my_mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sro_list, srr_list, sxo_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: iao_cubes_section, iao_molden_section, &
                                                            ibo_cc_section, ibo_cubes_section, &
                                                            ibo_molden_section

      ! only do IAO analysis if explicitly requested
      CPASSERT(iao_env%do_iao)

      ! k-points?
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspin = dft_control%nspins
      nimages = dft_control%nimages
      IF (nimages > 1) THEN
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T2,A)") &
               "K-Points: Intrinsic Atomic Orbitals Analysis not available."
         END IF
      END IF
      IF (nimages > 1) RETURN

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/,T2,A)') '!-----------------------------------------------------------------------------!'
         WRITE (UNIT=unit_nr, FMT="(T24,A)") "INTRINSIC ATOMIC ORBITALS ANALYSIS"
         WRITE (UNIT=unit_nr, FMT="(T13,A)") "G. Knizia, J. Chem. Theory Comput. 9, 4834-4843 (2013)"
         WRITE (unit_nr, '(T2,A)') '!-----------------------------------------------------------------------------!'
      END IF
      CALL cite_reference(Knizia2013)

      ! input sections
      iao_molden_section => iao_env%iao_molden_section
      iao_cubes_section => iao_env%iao_cubes_section
      ibo_molden_section => iao_env%ibo_molden_section
      ibo_cubes_section => iao_env%ibo_cubes_section
      ibo_cc_section => iao_env%ibo_cc_section

      ! check for or generate reference basis
      CALL create_minbas_set(qs_env, unit_nr)

      ! overlap matrices
      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s)
      smat => matrix_s(1, 1)%matrix
      !
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)
      ALLOCATE (ref_basis_set_list(nkind), orb_basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         NULLIFY (ref_basis_set_list(ikind)%gto_basis_set)
         NULLIFY (orb_basis_set_list(ikind)%gto_basis_set)
         NULLIFY (refbasis, orbbasis)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=orbbasis, basis_type="ORB")
         IF (ASSOCIATED(orbbasis)) orb_basis_set_list(ikind)%gto_basis_set => orbbasis
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=refbasis, basis_type="MIN")
         IF (ASSOCIATED(refbasis)) ref_basis_set_list(ikind)%gto_basis_set => refbasis
      END DO

      ! neighbor lists
      NULLIFY (srr_list, sro_list)
      CALL setup_neighbor_list(srr_list, ref_basis_set_list, qs_env=qs_env)
      CALL setup_neighbor_list(sro_list, ref_basis_set_list, orb_basis_set_list, qs_env=qs_env)

      ! Srr and Sro overlap matrices
      NULLIFY (sref, sro)
      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL build_overlap_matrix_simple(ks_env, sref, &
                                       ref_basis_set_list, ref_basis_set_list, srr_list)
      CALL build_overlap_matrix_simple(ks_env, sro, &
                                       ref_basis_set_list, orb_basis_set_list, sro_list)
      !
      IF (PRESENT(mos)) THEN
         my_mos => mos
      ELSE
         CALL get_qs_env(qs_env, mos=my_mos)
      END IF
      CALL get_mo_set(my_mos(1), mo_coeff=mo_coeff)
      CALL dbcsr_get_info(sro(1)%matrix, nfullrows_total=nref, nfullcols_total=nao)

      ! Projectors
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=mo_coeff%matrix_struct%context, nrow_global=nref, &
                               ncol_global=nao, para_env=mo_coeff%matrix_struct%para_env)
      CALL cp_fm_create(p_ref_orb, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=mo_coeff%matrix_struct%context, nrow_global=nao, &
                               ncol_global=nref, para_env=mo_coeff%matrix_struct%para_env)
      CALL cp_fm_create(p_orb_ref, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      !
      CALL iao_projectors(smat, sref(1)%matrix, sro(1)%matrix, p_orb_ref, p_ref_orb, iao_env%eps_svd)

      ! occupied orbitals, orbitals considered for IAO generation
      nocc(1:2) = 0
      DO ispin = 1, nspin
         CALL get_mo_set(my_mos(ispin), mo_coeff=mo_coeff, nmo=norb, uniform_occupation=uniform_occupation, &
                         occupation_numbers=occupation_numbers)
         IF (uniform_occupation) THEN
            nvec = norb
         ELSE
            nvec = 0
            DO i = 1, norb
               IF (occupation_numbers(i) > iao_env%eps_occ) THEN
                  nvec = i
               ELSE
                  EXIT
               END IF
            END DO
         END IF
         nocc(ispin) = nvec
      END DO
      ! generate IAOs
      ALLOCATE (iao_coef(nspin), cvec(nspin))
      DO ispin = 1, nspin
         CALL get_mo_set(my_mos(ispin), mo_coeff=mo_coeff)
         nvec = nocc(ispin)
         IF (nvec > 0) THEN
            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, ncol_global=nvec, template_fmstruct=mo_coeff%matrix_struct)
            CALL cp_fm_create(cvec(ispin), fm_struct)
            CALL cp_fm_struct_release(fm_struct)
            CALL cp_fm_to_fm(mo_coeff, cvec(ispin), nvec)
            !
            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, ncol_global=nref, template_fmstruct=mo_coeff%matrix_struct)
            CALL cp_fm_create(iao_coef(ispin), fm_struct)
            CALL cp_fm_struct_release(fm_struct)
            !
            CALL intrinsic_ao_calc(smat, p_orb_ref, p_ref_orb, cvec(ispin), iao_coef(ispin))
         END IF
      END DO
      !
      IF (iao_env%do_charges) THEN
         ! MOs in IAO basis
         ALLOCATE (ciao(nspin))
         DO ispin = 1, nspin
            CALL get_mo_set(my_mos(ispin), mo_coeff=mo_coeff, nao=nao, nmo=norb)
            CALL cp_fm_create(smo, mo_coeff%matrix_struct)
            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nref, template_fmstruct=mo_coeff%matrix_struct)
            CALL cp_fm_create(ciao(ispin), fm_struct)
            CALL cp_fm_struct_release(fm_struct)
            ! CIAO = A(T)*S*C
            CALL cp_dbcsr_sm_fm_multiply(smat, mo_coeff, smo, ncol=norb)
            CALL parallel_gemm('T', 'N', nref, norb, nao, 1.0_dp, iao_coef(ispin), smo, 0.0_dp, ciao(ispin))
            CALL cp_fm_release(smo)
         END DO
         !
         ! population analysis
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(/,T2,A)') 'Intrinsic AO Population Analysis '
         END IF
         CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, natom=natom)
         CALL get_ks_env(ks_env=ks_env, particle_set=particle_set)
         ALLOCATE (mcharge(natom, nspin))
         CALL dbcsr_get_info(sref(1)%matrix, nfullrows_total=nref)
         DO ispin = 1, nspin
            CALL get_mo_set(my_mos(ispin), &
                            uniform_occupation=uniform_occupation, &
                            occupation_numbers=occupation_numbers)
            ! diagonal block matrix
            CALL dbcsr_create(dmat, template=sref(1)%matrix)
            CALL dbcsr_reserve_diag_blocks(dmat)
            !
            CALL iao_calculate_dmat(ciao(ispin), dmat, occupation_numbers, uniform_occupation)
            CALL dbcsr_trace(dmat, trace)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,I2,T66,F15.4)') 'Number of Electrons: Trace(Piao) Spin ', ispin, trace
            END IF
            CALL iao_charges(dmat, mcharge(:, ispin))
            CALL dbcsr_release(dmat)
         END DO
         CALL print_atomic_charges(particle_set, qs_kind_set, unit_nr, "Intrinsic Atomic Orbital Charges", &
                                   electronic_charges=mcharge)
         DEALLOCATE (mcharge)
         CALL cp_fm_release(ciao)
      END IF
      !
      ! Deallocate the neighbor list structure
      CALL release_neighbor_list_sets(srr_list)
      CALL release_neighbor_list_sets(sro_list)
      IF (ASSOCIATED(sref)) CALL dbcsr_deallocate_matrix_set(sref)
      IF (ASSOCIATED(sro)) CALL dbcsr_deallocate_matrix_set(sro)
      CALL cp_fm_release(p_ref_orb)
      CALL cp_fm_release(p_orb_ref)
      CALL cp_fm_release(cvec)
      DEALLOCATE (orb_basis_set_list)
      !
      IF (iao_env%do_oce) THEN
         ! generate OCE basis
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A)') "IAO One-Center Expansion: OCE Basis Set"
         END IF
         CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
         nkind = SIZE(qs_kind_set)
         ALLOCATE (oce_basis_set_list(nkind), orb_basis_set_list(nkind))
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            NULLIFY (orb_basis_set_list(ikind)%gto_basis_set)
            NULLIFY (orbbasis)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=orbbasis, basis_type="ORB")
            IF (ASSOCIATED(orbbasis)) orb_basis_set_list(ikind)%gto_basis_set => orbbasis
         END DO
         DO ikind = 1, nkind
            orbbasis => orb_basis_set_list(ikind)%gto_basis_set
            CPASSERT(ASSOCIATED(orbbasis))
            NULLIFY (ocebasis)
            CALL create_oce_basis(ocebasis, orbbasis, iao_env%lmax_oce, iao_env%nbas_oce)
            CALL init_orb_basis_set(ocebasis)
            CALL init_interaction_radii_orb_basis(ocebasis, dft_control%qs_control%eps_pgf_orb)
            oce_basis_set_list(ikind)%gto_basis_set => ocebasis
            IF (unit_nr > 0) THEN
               qs_kind => qs_kind_set(ikind)
               CALL get_qs_kind(qs_kind, zeff=zval, element_symbol=element_symbol)
               CALL get_gto_basis_set(ocebasis, name=bname, nsgf=nsgf)
               WRITE (unit_nr, '(T2,A,A,T14,A,I4,T40,A,A30)') "Kind: ", element_symbol, "NBasFun: ", nsgf, &
                  "OCE Basis: ", ADJUSTL(TRIM(bname))
            END IF
         END DO
         IF (unit_nr > 0) WRITE (unit_nr, *)
         ! OCE basis overlap
         ALLOCATE (smat_kind(nkind))
         DO ikind = 1, nkind
            ocebasis => oce_basis_set_list(ikind)%gto_basis_set
            CALL get_gto_basis_set(gto_basis_set=ocebasis, nsgf=nsgf)
            ALLOCATE (smat_kind(ikind)%array(nsgf, nsgf))
         END DO
         CALL oce_overlap_matrix(smat_kind, oce_basis_set_list)
         ! overlap between IAO OCE basis and orbital basis
         NULLIFY (sxo, sxo_list)
         CALL setup_neighbor_list(sxo_list, oce_basis_set_list, orb_basis_set_list, qs_env=qs_env)
         CALL get_qs_env(qs_env, ks_env=ks_env)
         CALL build_overlap_matrix_simple(ks_env, sxo, &
                                          oce_basis_set_list, orb_basis_set_list, sxo_list)
         !
         ! One Center Expansion of IAOs
         CALL get_qs_env(qs_env=qs_env, natom=natom)
         ALLOCATE (oce_atom(natom, nspin))
         CALL iao_oce_expansion(qs_env, oce_basis_set_list, smat_kind, sxo, iao_coef, oce_atom, unit_nr)
         !
         DO ikind = 1, nkind
            ocebasis => oce_basis_set_list(ikind)%gto_basis_set
            CALL deallocate_gto_basis_set(ocebasis)
         END DO
         DEALLOCATE (oce_basis_set_list, orb_basis_set_list)
         !
         CALL release_neighbor_list_sets(sxo_list)
         IF (ASSOCIATED(sxo)) CALL dbcsr_deallocate_matrix_set(sxo)
         DO ikind = 1, nkind
            DEALLOCATE (smat_kind(ikind)%array)
         END DO
         DEALLOCATE (smat_kind)
         DO ispin = 1, nspin
            DO iatom = 1, natom
               DEALLOCATE (oce_atom(iatom, ispin)%array)
            END DO
         END DO
         DEALLOCATE (oce_atom)
      END IF
      !
      IF (molden_iao) THEN
         ! Print basis functions: molden file
         CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
         CALL get_ks_env(ks_env=ks_env, particle_set=particle_set)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A)') '  Write IAO in MOLDEN Format'
         END IF
         ALLOCATE (mo_iao(nspin))
         DO ispin = 1, nspin
            CALL get_mo_set(my_mos(ispin), nao=nao)
            CALL allocate_mo_set(mo_iao(ispin), nao, nref, nref, 0.0_dp, 1.0_dp, 0.0_dp)
            CALL init_mo_set(mo_iao(ispin), fm_ref=iao_coef(ispin), name="iao_set")
            CALL cp_fm_to_fm(iao_coef(ispin), mo_iao(ispin)%mo_coeff)
            CALL set_mo_set(mo_iao(ispin), homo=nref)
            mo_iao(ispin)%occupation_numbers = 1.0_dp
         END DO
         ! Print IAO basis functions: molden format
         CALL write_mos_molden(mo_iao, qs_kind_set, particle_set, iao_molden_section)
         DO ispin = 1, nspin
            CALL deallocate_mo_set(mo_iao(ispin))
         END DO
         DEALLOCATE (mo_iao)
      END IF
      IF (cubes_iao) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A)') '  Write IAO as CUBE Files'
         END IF
         ! Print basis functions: cube file
         CALL print_iao_cubes(qs_env, iao_cubes_section, iao_coef, ref_basis_set_list)
      END IF
      !
      ! Intrinsic bond orbitals
      IF (iao_env%do_bondorbitals) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(/,T2,A)') 'Intrinsic Bond Orbital Generation'
         END IF
         ! MOs in IAO basis -> ciao
         ALLOCATE (ciao(nspin))
         DO ispin = 1, nspin
            CALL get_mo_set(my_mos(ispin), mo_coeff=mo_coeff, nao=nao, nmo=norb)
            CALL cp_fm_create(smo, mo_coeff%matrix_struct)
            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nref, template_fmstruct=mo_coeff%matrix_struct)
            CALL cp_fm_create(ciao(ispin), fm_struct)
            CALL cp_fm_struct_release(fm_struct)
            ! CIAO = A(T)*S*C
            CALL cp_dbcsr_sm_fm_multiply(smat, mo_coeff, smo, ncol=norb)
            CALL parallel_gemm('T', 'N', nref, norb, nao, 1.0_dp, iao_coef(ispin), smo, 0.0_dp, ciao(ispin))
            CALL cp_fm_release(smo)
         END DO
         !
         ! localize occupied orbitals nocc(ispin), see IAO generation
         !
         ALLOCATE (cloc(nspin), c_loc_orb(nspin))
         DO ispin = 1, nspin
            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, ncol_global=nocc(ispin), &
                                     template_fmstruct=ciao(ispin)%matrix_struct)
            CALL cp_fm_create(cloc(ispin), fm_struct)
            CALL cp_fm_struct_release(fm_struct)
            CALL cp_fm_to_fm(ciao(ispin), cloc(ispin), ncol=nocc(ispin))
         END DO
         CALL cp_fm_release(ciao)
         CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, natom=natom)
         ALLOCATE (first_sgf(natom), last_sgf(natom), nsgfat(natom))
         CALL get_ks_env(ks_env=ks_env, particle_set=particle_set)
         CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf, &
                               nsgf=nsgfat, basis=ref_basis_set_list)

         IF (iao_env%loc_operator /= do_iaoloc_pm2) THEN
            CPABORT("IAO localization operator NYA")
         END IF
         IF (iao_env%eloc_function /= do_iaoloc_enone) THEN
            CPABORT("IAO energy weight localization NYA")
         END IF

         DO ispin = 1, nspin
            nvec = nocc(ispin)
            IF (nvec > 0) THEN
               ALLOCATE (zij_atom(natom, 1))
               NULLIFY (fm_struct)
               CALL cp_fm_struct_create(fm_struct, nrow_global=nvec, ncol_global=nvec, &
                                        template_fmstruct=cloc(ispin)%matrix_struct)
               DO iatom = 1, natom
                  CALL cp_fm_create(zij_atom(iatom, 1), fm_struct)
                  isgf = first_sgf(iatom)
                  CALL parallel_gemm('T', 'N', nvec, nvec, nsgfat(iatom), 1.0_dp, cloc(ispin), cloc(ispin), &
                                     0.0_dp, zij_atom(iatom, 1), &
                                     a_first_col=1, a_first_row=isgf, b_first_col=1, b_first_row=isgf)
               END DO
               CALL cp_fm_struct_release(fm_struct)
               !
               CALL jacobi_rotation_pipek(zij_atom, cloc(ispin), sweeps)
               !
               CALL cp_fm_release(zij_atom)
               CALL get_mo_set(my_mos(ispin), nao=nao)
               NULLIFY (fm_struct)
               CALL cp_fm_struct_create(fm_struct, nrow_global=nao, &
                                        template_fmstruct=cloc(ispin)%matrix_struct)
               CALL cp_fm_create(c_loc_orb(ispin), fm_struct)
               CALL cp_fm_struct_release(fm_struct)
               CALL parallel_gemm('N', 'N', nao, nvec, nref, 1.0_dp, iao_coef(ispin), cloc(ispin), &
                                  0.0_dp, c_loc_orb(ispin))
            END IF
         END DO
         DEALLOCATE (first_sgf, last_sgf, nsgfat)
         CALL cp_fm_release(cloc)
         !
         IF (iao_env%do_center) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A)') '  Calculate Localized Orbital Centers and Spread'
               IF (iao_env%pos_periodic) THEN
                  WRITE (unit_nr, '(T2,A)') '    Use Berry Phase Position Operator'
               ELSE
                  WRITE (unit_nr, '(T2,A)') '    Use Local Position Operator'
               END IF
            END IF
            nvec = MAXVAL(nocc)
            ! x  y  z  m2(1)  m2(2)
            ALLOCATE (moments(5, nvec, nspin))
            moments = 0.0_dp
            IF (iao_env%pos_periodic) THEN
               CALL center_spread_berry(qs_env, c_loc_orb, moments)
            ELSE
               CALL center_spread_loc(qs_env, c_loc_orb, moments)
            END IF
            !
            IF (ASSOCIATED(ibo_cc_section)) THEN
               CALL print_center_spread(moments, nocc, ibo_cc_section)
            END IF
            !
            IF (PRESENT(bond_centers)) THEN
               nx = SIZE(bond_centers, 1)
               no = SIZE(bond_centers, 2)
               ns = SIZE(bond_centers, 3)
               CPASSERT(no >= nvec)
               CPASSERT(ns == nspin)
               CPASSERT(nx >= 3)
               nx = MIN(nx, 5)
               bond_centers(1:nx, 1:nvec, 1:ns) = moments(1:nx, 1:nvec, 1:ns)
            END IF
            DEALLOCATE (moments)
         END IF
         !
         IF (molden_ibo) THEN
            CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
            CALL get_ks_env(ks_env=ks_env, particle_set=particle_set)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A)') '  Write IBO in MOLDEN Format'
            END IF
            ALLOCATE (mo_loc(nspin))
            DO ispin = 1, nspin
               CALL get_mo_set(my_mos(ispin), nao=nao)
               nvec = nocc(ispin)
               CALL allocate_mo_set(mo_loc(ispin), nao, nvec, nvec, 0.0_dp, 1.0_dp, 0.0_dp)
               CALL init_mo_set(mo_loc(ispin), fm_ref=c_loc_orb(ispin), name="ibo_orb")
               CALL cp_fm_to_fm(c_loc_orb(ispin), mo_loc(ispin)%mo_coeff)
               CALL set_mo_set(mo_loc(ispin), homo=nvec)
               mo_loc(ispin)%occupation_numbers = 1.0_dp
            END DO
            ! Print IBO functions: molden format
            CALL write_mos_molden(mo_loc, qs_kind_set, particle_set, ibo_molden_section)
            DO ispin = 1, nspin
               CALL deallocate_mo_set(mo_loc(ispin))
            END DO
            DEALLOCATE (mo_loc)
         END IF
         !
         IF (cubes_ibo) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A)') '  Write IBO on CUBE files'
            END IF
            ! Print localized orbital functions: cube file
            CALL print_ibo_cubes(qs_env, ibo_cubes_section, c_loc_orb)
         END IF
         !
         IF (PRESENT(mos)) THEN
            ! c_loc_orb -> mos
            DO ispin = 1, nspin
               CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
               nvec = nocc(ispin)
               CALL cp_fm_to_fm(c_loc_orb(ispin), mo_coeff, ncol=nvec)
            END DO
         END IF
         CALL cp_fm_release(c_loc_orb)
      END IF
      DEALLOCATE (ref_basis_set_list)

      IF (PRESENT(c_iao_coef)) THEN
         CALL cp_fm_release(c_iao_coef)
         ALLOCATE (c_iao_coef(nspin))
         DO ispin = 1, nspin
            CALL cp_fm_create(c_iao_coef(ispin), iao_coef(ispin)%matrix_struct)
            CALL cp_fm_to_fm(iao_coef(ispin), c_iao_coef(ispin))
         END DO
      END IF
      CALL cp_fm_release(iao_coef)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') &
            '!----------------------------END OF IAO ANALYSIS------------------------------!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE iao_wfn_analysis

! **************************************************************************************************
!> \brief Computes projector matrices for ref to orb basis and reverse
!> \param smat ...
!> \param sref ...
!> \param s_r_o ...
!> \param p_o_r ...
!> \param p_r_o ...
!> \param eps_svd ...
! **************************************************************************************************
   SUBROUTINE iao_projectors(smat, sref, s_r_o, p_o_r, p_r_o, eps_svd)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: smat, sref, s_r_o
      TYPE(cp_fm_type), INTENT(INOUT)                    :: p_o_r, p_r_o
      REAL(KIND=dp), INTENT(IN)                          :: eps_svd

      CHARACTER(len=*), PARAMETER                        :: routineN = 'iao_projectors'

      INTEGER                                            :: handle, norb, nref
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_inv, fm_s, fm_sro

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(p_r_o, nrow_global=nref, ncol_global=norb)
      CALL cp_fm_create(fm_sro, p_r_o%matrix_struct)
      CALL copy_dbcsr_to_fm(s_r_o, fm_sro)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, nrow_global=norb, ncol_global=norb, &
                               template_fmstruct=p_r_o%matrix_struct)
      CALL cp_fm_create(fm_s, fm_struct, name="smat")
      CALL cp_fm_create(fm_inv, fm_struct, name="sinv")
      CALL cp_fm_struct_release(fm_struct)
      CALL copy_dbcsr_to_fm(smat, fm_s)
      CALL cp_fm_invert(fm_s, fm_inv, eps_svd=eps_svd)
      CALL parallel_gemm('N', 'T', norb, nref, norb, 1.0_dp, fm_inv, fm_sro, 0.0_dp, p_o_r)
      CALL cp_fm_release(fm_s)
      CALL cp_fm_release(fm_inv)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nref, ncol_global=nref, &
                               template_fmstruct=p_r_o%matrix_struct)
      CALL cp_fm_create(fm_s, fm_struct, name="sref")
      CALL cp_fm_create(fm_inv, fm_struct, name="sinv")
      CALL cp_fm_struct_release(fm_struct)
      CALL copy_dbcsr_to_fm(sref, fm_s)
      CALL cp_fm_invert(fm_s, fm_inv, eps_svd=eps_svd)
      CALL parallel_gemm('N', 'N', nref, norb, nref, 1.0_dp, fm_inv, fm_sro, 0.0_dp, p_r_o)
      CALL cp_fm_release(fm_s)
      CALL cp_fm_release(fm_inv)

      CALL cp_fm_release(fm_sro)

      CALL timestop(handle)

   END SUBROUTINE iao_projectors

! **************************************************************************************************
!> \brief Computes intrinsic orbitals for a given MO vector set
!> \param smat ...
!> \param p_o_r ...
!> \param p_r_o ...
!> \param cvec ...
!> \param avec ...
! **************************************************************************************************
   SUBROUTINE intrinsic_ao_calc(smat, p_o_r, p_r_o, cvec, avec)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: smat
      TYPE(cp_fm_type), INTENT(INOUT)                    :: p_o_r, p_r_o, cvec, avec

      CHARACTER(len=*), PARAMETER                        :: routineN = 'intrinsic_ao_calc'

      INTEGER                                            :: handle, nao, norb, nref
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: ctvec, pc, sc, sct, vec1

      CALL timeset(routineN, handle)

      ! number of orbitals
      CALL cp_fm_get_info(cvec, ncol_global=norb)
      ! basis set sizes
      CALL cp_fm_get_info(p_r_o, nrow_global=nref, ncol_global=nao)
      ! temp array
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nref, ncol_global=norb, &
                               template_fmstruct=cvec%matrix_struct)
      CALL cp_fm_create(pc, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      ! CT = orth(Por.Pro.C)
      CALL parallel_gemm('N', 'N', nref, norb, nao, 1.0_dp, p_r_o, cvec, 0.0_dp, pc)
      CALL cp_fm_create(ctvec, cvec%matrix_struct)
      CALL parallel_gemm('N', 'N', nao, norb, nref, 1.0_dp, p_o_r, pc, 0.0_dp, ctvec)
      CALL cp_fm_release(pc)
      CALL make_basis_lowdin(ctvec, norb, smat)
      ! S*C and S*CT
      CALL cp_fm_create(sc, cvec%matrix_struct)
      CALL cp_dbcsr_sm_fm_multiply(smat, cvec, sc, ncol=norb)
      CALL cp_fm_create(sct, cvec%matrix_struct)
      CALL cp_dbcsr_sm_fm_multiply(smat, ctvec, sct, ncol=norb)
      CALL cp_fm_struct_create(fm_struct, nrow_global=norb, ncol_global=nref, &
                               template_fmstruct=cvec%matrix_struct)
      CALL cp_fm_create(pc, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      ! V1 = (CT*SCT(T))Por
      CALL parallel_gemm('T', 'N', norb, nref, nao, 1.0_dp, sct, p_o_r, 0.0_dp, pc)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nref, &
                               template_fmstruct=cvec%matrix_struct)
      CALL cp_fm_create(vec1, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL parallel_gemm('N', 'N', nao, nref, norb, 1.0_dp, ctvec, pc, 0.0_dp, vec1)
      ! A = C*SC(T)*V1
      CALL parallel_gemm('T', 'N', norb, nref, nao, 1.0_dp, sc, vec1, 0.0_dp, pc)
      CALL parallel_gemm('N', 'N', nao, nref, norb, 1.0_dp, cvec, pc, 0.0_dp, avec)
      ! V1 = Por - V1
      CALL cp_fm_geadd(1.0_dp, 'N', p_o_r, -1.0_dp, vec1)
      ! A = A + V1
      CALL cp_fm_geadd(1.0_dp, 'N', vec1, 1.0_dp, avec)
      ! A = A - C*SC(T)*V1
      CALL parallel_gemm('T', 'N', norb, nref, nao, 1.0_dp, sc, vec1, 0.0_dp, pc)
      CALL parallel_gemm('N', 'N', nao, nref, norb, -1.0_dp, cvec, pc, 1.0_dp, avec)
      ! A = orth(A)
      CALL make_basis_lowdin(avec, nref, smat)

      ! clean up
      CALL cp_fm_release(pc)
      CALL cp_fm_release(vec1)
      CALL cp_fm_release(sc)
      CALL cp_fm_release(sct)
      CALL cp_fm_release(ctvec)

      CALL timestop(handle)

   END SUBROUTINE intrinsic_ao_calc

! **************************************************************************************************
!> \brief Calculate the density matrix from fm vectors using occupation numbers
!> \param cvec ...
!> \param density_matrix ...
!> \param occupation ...
!> \param uniform_occupation ...
! **************************************************************************************************
   SUBROUTINE iao_calculate_dmat_diag(cvec, density_matrix, occupation, uniform_occupation)

      TYPE(cp_fm_type), INTENT(IN)                       :: cvec
      TYPE(dbcsr_type)                                   :: density_matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation
      LOGICAL, INTENT(IN)                                :: uniform_occupation

      CHARACTER(len=*), PARAMETER :: routineN = 'iao_calculate_dmat_diag'

      INTEGER                                            :: handle, ncol
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_fm_type)                                   :: fm_tmp

      CALL timeset(routineN, handle)

      CALL dbcsr_set(density_matrix, 0.0_dp)

      CALL cp_fm_get_info(cvec, ncol_global=ncol)
      IF (.NOT. uniform_occupation) THEN
         CALL cp_fm_create(fm_tmp, cvec%matrix_struct)
         CALL cp_fm_to_fm(cvec, fm_tmp)
         CALL cp_fm_column_scale(fm_tmp, occupation(1:ncol))
         alpha = 1.0_dp
         CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix, &
                                    matrix_v=cvec, matrix_g=fm_tmp, &
                                    ncol=ncol, alpha=alpha)
         CALL cp_fm_release(fm_tmp)
      ELSE
         alpha = occupation(1)
         CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix, &
                                    matrix_v=cvec, ncol=ncol, alpha=alpha)
      END IF

      CALL timestop(handle)

   END SUBROUTINE iao_calculate_dmat_diag

! **************************************************************************************************
!> \brief Calculate the density matrix from fm vectors using an occupation matrix
!> \param cvec ...
!> \param density_matrix ...
!> \param weight ...
!> \param occmat ...
! **************************************************************************************************
   SUBROUTINE iao_calculate_dmat_full(cvec, density_matrix, weight, occmat)

      TYPE(cp_fm_type), INTENT(IN)                       :: cvec
      TYPE(dbcsr_type)                                   :: density_matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: weight
      TYPE(cp_fm_type), INTENT(IN)                       :: occmat

      CHARACTER(len=*), PARAMETER :: routineN = 'iao_calculate_dmat_full'

      INTEGER                                            :: handle, ic, jc, ncol
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm1, fm2

      CALL timeset(routineN, handle)

      CALL dbcsr_set(density_matrix, 0.0_dp)

      CALL cp_fm_struct_create(fm_struct, ncol_global=1, &
                               template_fmstruct=cvec%matrix_struct)
      CALL cp_fm_create(fm1, fm_struct)
      CALL cp_fm_create(fm2, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_get_info(cvec, ncol_global=ncol)
      DO ic = 1, ncol
         CALL cp_fm_to_fm(cvec, fm1, ncol=1, source_start=ic, target_start=1)
         DO jc = 1, ncol
            CALL cp_fm_to_fm(cvec, fm2, ncol=1, source_start=jc, target_start=1)
            CALL cp_fm_get_element(occmat, ic, jc, alpha)
            alpha = weight(ic)*alpha
            CALL cp_dbcsr_plus_fm_fm_t(density_matrix, fm1, fm2, ncol=1, &
                                       alpha=alpha, symmetry_mode=1)
         END DO
      END DO
      CALL cp_fm_release(fm1)
      CALL cp_fm_release(fm2)

      CALL timestop(handle)

   END SUBROUTINE iao_calculate_dmat_full

! **************************************************************************************************
!> \brief compute the atomic charges (orthogonal basis)
!> \param p_matrix ...
!> \param charges ...
! **************************************************************************************************
   SUBROUTINE iao_charges(p_matrix, charges)
      TYPE(dbcsr_type)                                   :: p_matrix
      REAL(KIND=dp), DIMENSION(:)                        :: charges

      INTEGER                                            :: blk, group_handle, i, iblock_col, &
                                                            iblock_row
      REAL(kind=dp)                                      :: trace
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: group

      charges = 0.0_dp

      CALL dbcsr_iterator_start(iter, p_matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         NULLIFY (p_block)
         CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, p_block, blk)
         CPASSERT(iblock_row == iblock_col)
         trace = 0.0_dp
         DO i = 1, SIZE(p_block, 1)
            trace = trace + p_block(i, i)
         END DO
         charges(iblock_row) = trace
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_get_info(p_matrix, group=group_handle)
      CALL group%set_handle(group_handle)

      CALL group%sum(charges)

   END SUBROUTINE iao_charges

! **************************************************************************************************
!> \brief Computes and prints the Cube Files for the intrinsic atomic orbitals
!> \param qs_env ...
!> \param print_section ...
!> \param iao_coef ...
!> \param basis_set_list ...
! **************************************************************************************************
   SUBROUTINE print_iao_cubes(qs_env, print_section, iao_coef, basis_set_list)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: print_section
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: iao_coef
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list

      CHARACTER(LEN=default_path_length)                 :: filename, title
      INTEGER                                            :: i, i_rep, ispin, ivec, iw, j, n_rep, &
                                                            nat, natom, norb, nspins, nstart
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, blk_sizes, first_bas, stride
      LOGICAL                                            :: mpi_io
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: wf_g, wf_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys

      logger => cp_get_default_logger()
      stride => section_get_ivals(print_section, "STRIDE")
      CALL section_vals_val_get(print_section, "ATOM_LIST", n_rep_val=n_rep)

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      subsys=subsys, cell=cell, particle_set=particle_set, pw_env=pw_env, dft_control=dft_control)
      CALL qs_subsys_get(subsys, particles=particles)

      CALL get_qs_env(qs_env=qs_env, natom=natom)
      ALLOCATE (blk_sizes(natom), first_bas(0:natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=blk_sizes, basis=basis_set_list)
      first_bas(0) = 0
      DO i = 1, natom
         first_bas(i) = first_bas(i - 1) + blk_sizes(i)
      END DO

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_r, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_g, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      nspins = SIZE(iao_coef)
      nstart = MIN(1, n_rep)

      DO ispin = 1, nspins
         DO i_rep = nstart, n_rep
            CALL cp_fm_get_info(iao_coef(ispin), ncol_global=norb)
            IF (i_rep == 0) THEN
               nat = natom
            ELSE
               CALL section_vals_val_get(print_section, "ATOM_LIST", i_rep_val=i_rep, i_vals=atom_list)
               nat = SIZE(atom_list)
            END IF
            DO i = 1, nat
               IF (i_rep == 0) THEN
                  j = i
               ELSE
                  j = atom_list(i)
               END IF
               CPASSERT(j >= 1 .AND. j <= natom)
               DO ivec = first_bas(j - 1) + 1, first_bas(j)
                  WRITE (filename, '(a4,I5.5,a1,I1.1)') "IAO_", ivec, "_", ispin
                  WRITE (title, *) "Intrinsic Atomic Orbitals ", ivec, " atom ", j, " spin ", ispin
                  mpi_io = .TRUE.
                  iw = cp_print_key_unit_nr(logger, print_section, "", extension=".cube", &
                                            middle_name=TRIM(filename), file_position="REWIND", log_filename=.FALSE., &
                                            mpi_io=mpi_io)
                  CALL calculate_wavefunction(iao_coef(ispin), ivec, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                              cell, dft_control, particle_set, pw_env)
                  CALL cp_pw_to_cube(wf_r, iw, title, particles=particles, stride=stride, mpi_io=mpi_io)
                  CALL cp_print_key_finished_output(iw, logger, print_section, "", mpi_io=mpi_io)
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (blk_sizes, first_bas)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g)

   END SUBROUTINE print_iao_cubes

! **************************************************************************************************
!> \brief Computes and prints the Cube Files for the intrinsic bond orbitals
!> \param qs_env ...
!> \param print_section ...
!> \param ibo_coef ...
! **************************************************************************************************
   SUBROUTINE print_ibo_cubes(qs_env, print_section, ibo_coef)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: print_section
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: ibo_coef

      CHARACTER(LEN=default_path_length)                 :: filename, title
      INTEGER                                            :: i, i_rep, ispin, iw, j, n_rep, norb, &
                                                            nspins, nstart, nstate
      INTEGER, DIMENSION(:), POINTER                     :: state_list, stride
      LOGICAL                                            :: mpi_io
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: wf_g, wf_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CPASSERT(ASSOCIATED(print_section))

      logger => cp_get_default_logger()
      stride => section_get_ivals(print_section, "STRIDE")
      CALL section_vals_val_get(print_section, "STATE_LIST", n_rep_val=n_rep)

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      subsys=subsys, cell=cell, particle_set=particle_set, pw_env=pw_env, dft_control=dft_control)
      CALL qs_subsys_get(subsys, particles=particles)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_r, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_g, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      nspins = SIZE(ibo_coef)
      nstart = MIN(1, n_rep)

      DO ispin = 1, nspins
         DO i_rep = nstart, n_rep
            CALL cp_fm_get_info(ibo_coef(ispin), ncol_global=norb)
            IF (i_rep == 0) THEN
               nstate = norb
            ELSE
               CALL section_vals_val_get(print_section, "STATE_LIST", i_rep_val=i_rep, i_vals=state_list)
               nstate = SIZE(state_list)
            END IF
            DO i = 1, nstate
               IF (i_rep == 0) THEN
                  j = i
               ELSE
                  j = state_list(i)
               END IF
               CPASSERT(j >= 1 .AND. j <= norb)
               WRITE (filename, '(a4,I5.5,a1,I1.1)') "IBO_", j, "_", ispin
               WRITE (title, *) "Intrinsic Atomic Orbitals ", j, " spin ", ispin
               mpi_io = .TRUE.
               iw = cp_print_key_unit_nr(logger, print_section, "", extension=".cube", &
                                         middle_name=TRIM(filename), file_position="REWIND", log_filename=.FALSE., &
                                         mpi_io=mpi_io)
               CALL calculate_wavefunction(ibo_coef(ispin), j, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                           cell, dft_control, particle_set, pw_env)
               CALL cp_pw_to_cube(wf_r, iw, title, particles=particles, stride=stride, mpi_io=mpi_io)
               CALL cp_print_key_finished_output(iw, logger, print_section, "", mpi_io=mpi_io)
            END DO
         END DO
      END DO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g)

   END SUBROUTINE print_ibo_cubes

! **************************************************************************************************
!> \brief prints charge center and spreads for all orbitals
!> \param moments ...
!> \param nocc ...
!> \param print_section ...
! **************************************************************************************************
   SUBROUTINE print_center_spread(moments, nocc, print_section)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: moments
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nocc
      TYPE(section_vals_type), POINTER                   :: print_section

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: is, ispin, iw, nspin
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()
      nspin = SIZE(moments, 3)

      WRITE (filename, '(A18,I1.1)') "IBO_CENTERS_SPREAD"
      iw = cp_print_key_unit_nr(logger, print_section, "", extension=".csp", &
                                middle_name=TRIM(filename), file_position="REWIND", log_filename=.FALSE.)
      IF (iw > 0) THEN
         DO ispin = 1, nspin
            WRITE (iw, "(A,i1)") "Intrinsic Bond Orbitals: Centers and Spread for Spin ", ispin
            WRITE (iw, "(A7,T30,A6,T68,A7)") "State", "Center", "Spreads"
            DO is = 1, nocc(ispin)
               WRITE (iw, "(i7,3F15.8,8X,2F10.5)") is, moments(1:5, is, ispin)
            END DO
         END DO
      END IF
      CALL cp_print_key_finished_output(iw, logger, print_section, "")

   END SUBROUTINE print_center_spread

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param c_loc_orb ...
!> \param moments ...
! **************************************************************************************************
   SUBROUTINE center_spread_loc(qs_env, c_loc_orb, moments)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: c_loc_orb
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: moments

      CHARACTER(len=*), PARAMETER                        :: routineN = 'center_spread_loc'
      INTEGER, DIMENSION(6), PARAMETER                   :: list = (/1, 2, 3, 4, 7, 9/)

      INTEGER                                            :: handle, i, iop, iorb, ispin, nao, norb, &
                                                            nspin
      REAL(KIND=dp)                                      :: xmii
      REAL(KIND=dp), DIMENSION(3)                        :: rpoint
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: ccmat, ocvec
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipmat
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: omat, smat

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrix_s)
      smat => matrix_s(1, 1)%matrix
      rpoint = 0.0_dp
      nspin = SIZE(c_loc_orb)
      moments = 0.0_dp

      ALLOCATE (dipmat(9))
      DO i = 1, 9
         ALLOCATE (dipmat(i)%matrix)
         CALL dbcsr_copy(dipmat(i)%matrix, smat)
         CALL dbcsr_set(dipmat(i)%matrix, 0.0_dp)
      END DO

      CALL build_local_moment_matrix(qs_env, dipmat, 2, rpoint)

      DO ispin = 1, nspin
         CALL cp_fm_get_info(c_loc_orb(ispin), nrow_global=nao, ncol_global=norb)
         CALL cp_fm_create(ocvec, c_loc_orb(ispin)%matrix_struct)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, nrow_global=norb, ncol_global=norb, &
                                  template_fmstruct=c_loc_orb(ispin)%matrix_struct)
         CALL cp_fm_create(ccmat, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
         DO i = 1, 6
            iop = list(i)
            omat => dipmat(iop)%matrix
            CALL cp_dbcsr_sm_fm_multiply(omat, c_loc_orb(ispin), ocvec, ncol=norb)
            CALL parallel_gemm('T', 'N', norb, norb, nao, 1.0_dp, c_loc_orb(ispin), ocvec, &
                               0.0_dp, ccmat)
            DO iorb = 1, norb
               CALL cp_fm_get_element(ccmat, iorb, iorb, xmii)
               IF (iop <= 3) THEN
                  moments(iop, iorb, ispin) = moments(iop, iorb, ispin) + xmii
                  moments(4, iorb, ispin) = moments(4, iorb, ispin) - xmii**2
               ELSE
                  moments(4, iorb, ispin) = moments(4, iorb, ispin) + xmii
               END IF
            END DO
         END DO
         CALL cp_fm_release(ocvec)
         CALL cp_fm_release(ccmat)
      END DO

      DO i = 1, 9
         CALL dbcsr_release(dipmat(i)%matrix)
         DEALLOCATE (dipmat(i)%matrix)
      END DO
      DEALLOCATE (dipmat)

      CALL get_qs_env(qs_env=qs_env, cell=cell)
      DO ispin = 1, nspin
         DO iorb = 1, norb
            moments(1:3, iorb, ispin) = pbc(moments(1:3, iorb, ispin), cell)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE center_spread_loc

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param c_loc_orb ...
!> \param moments ...
! **************************************************************************************************
   SUBROUTINE center_spread_berry(qs_env, c_loc_orb, moments)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: c_loc_orb
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: moments

      CHARACTER(len=*), PARAMETER :: routineN = 'center_spread_berry'

      COMPLEX(KIND=dp)                                   :: z
      INTEGER                                            :: dim_op, handle, i, idir, ispin, istate, &
                                                            j, jdir, nao, norb, nspin
      REAL(dp), DIMENSION(3)                             :: c, cpbc
      REAL(dp), DIMENSION(6)                             :: weights
      REAL(KIND=dp)                                      :: imagpart, realpart, spread_i, spread_ii
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: opvec
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: zij
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: op_sm_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, cell=cell)

      IF (cell%orthorhombic) THEN
         dim_op = 3
      ELSE
         dim_op = 6
      END IF
      ALLOCATE (op_sm_set(2, dim_op))
      DO i = 1, dim_op
         DO j = 1, 2
            NULLIFY (op_sm_set(j, i)%matrix)
            ALLOCATE (op_sm_set(j, i)%matrix)
            CALL dbcsr_copy(op_sm_set(j, i)%matrix, matrix_s(1)%matrix)
            CALL dbcsr_set(op_sm_set(j, i)%matrix, 0.0_dp)
         END DO
      END DO
      CALL initialize_weights(cell, weights)

      CALL compute_berry_operator(qs_env, cell, op_sm_set, dim_op)

      nspin = SIZE(c_loc_orb, 1)
      DO ispin = 1, nspin
         CALL cp_fm_get_info(c_loc_orb(ispin), nrow_global=nao, ncol_global=norb)
         CALL cp_fm_create(opvec, c_loc_orb(ispin)%matrix_struct)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, nrow_global=norb, ncol_global=norb, &
                                  template_fmstruct=c_loc_orb(ispin)%matrix_struct)
         ALLOCATE (zij(2, dim_op))

         ! Compute zij here
         DO i = 1, dim_op
            DO j = 1, 2
               CALL cp_fm_create(zij(j, i), fm_struct)
               CALL cp_fm_set_all(zij(j, i), 0.0_dp)
               CALL cp_dbcsr_sm_fm_multiply(op_sm_set(j, i)%matrix, c_loc_orb(ispin), opvec, ncol=norb)
               CALL parallel_gemm("T", "N", norb, norb, nao, 1.0_dp, c_loc_orb(ispin), opvec, 0.0_dp, zij(j, i))
            END DO
         END DO

         CALL cp_fm_struct_release(fm_struct)
         CALL cp_fm_release(opvec)

         DO istate = 1, norb
            c = 0.0_dp
            spread_i = 0.0_dp
            spread_ii = 0.0_dp
            DO jdir = 1, dim_op
               CALL cp_fm_get_element(zij(1, jdir), istate, istate, realpart)
               CALL cp_fm_get_element(zij(2, jdir), istate, istate, imagpart)
               z = CMPLX(realpart, imagpart, dp)
               spread_i = spread_i - weights(jdir)* &
                          LOG(realpart*realpart + imagpart*imagpart)/twopi/twopi
               spread_ii = spread_ii + weights(jdir)* &
                           (1.0_dp - (realpart*realpart + imagpart*imagpart))/twopi/twopi
               IF (jdir < 4) THEN
                  DO idir = 1, 3
                     c(idir) = c(idir) + &
                               (cell%hmat(idir, jdir)/twopi)*AIMAG(LOG(z))
                  END DO
               END IF
            END DO
            cpbc = pbc(c, cell)
            moments(1:3, istate, ispin) = cpbc(1:3)
            moments(4, istate, ispin) = spread_i
            moments(5, istate, ispin) = spread_ii
         END DO

         CALL cp_fm_release(zij)

      END DO

      DO i = 1, dim_op
         DO j = 1, 2
            CALL dbcsr_release(op_sm_set(j, i)%matrix)
            DEALLOCATE (op_sm_set(j, i)%matrix)
         END DO
      END DO
      DEALLOCATE (op_sm_set)

      CALL timestop(handle)

   END SUBROUTINE center_spread_berry

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param oce_basis_set_list ...
!> \param smat_kind ...
!> \param sxo ...
!> \param iao_coef ...
!> \param oce_atom ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE iao_oce_expansion(qs_env, oce_basis_set_list, smat_kind, sxo, iao_coef, oce_atom, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: oce_basis_set_list
      TYPE(cp_2d_r_p_type), DIMENSION(:)                 :: smat_kind
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: sxo
      TYPE(cp_fm_type), DIMENSION(:)                     :: iao_coef
      TYPE(cp_2d_r_p_type), DIMENSION(:, :)              :: oce_atom
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER                        :: routineN = 'iao_oce_expansion'

      INTEGER                                            :: handle, i, i1, i2, iao, iatom, ikind, &
                                                            iset, ishell, ispin, l, m, maxl, n, &
                                                            natom, nkind, noce, ns, nset, nsgf, &
                                                            nspin, para_group_handle
      INTEGER, DIMENSION(:), POINTER                     :: iao_blk_sizes, nshell, oce_blk_sizes, &
                                                            orb_blk_sizes
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf, last_sgf, lval
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ev, vector
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat, bmat, filter, oce_comp, prol, vec
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ablock
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: sinv_kind
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_type)                                   :: iao_vec, sx_vec
      TYPE(gto_basis_set_type), POINTER                  :: oce_basis
      TYPE(mp_comm_type)                                 :: para_type
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      ! basis sets: block sizes
      CALL dbcsr_get_info(sxo(1)%matrix, row_blk_size=oce_blk_sizes, col_blk_size=orb_blk_sizes, &
                          distribution=dbcsr_dist)
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, natom=natom)
      CALL get_ks_env(ks_env=ks_env, particle_set=particle_set, qs_kind_set=qs_kind_set)
      ALLOCATE (iao_blk_sizes(natom))
      DO iatom = 1, natom
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), nsgf=ns, basis_type="MIN")
         iao_blk_sizes(iatom) = ns
      END DO

      CALL dbcsr_create(iao_vec, name="IAO_VEC", dist=dbcsr_dist, &
                        matrix_type=dbcsr_type_no_symmetry, row_blk_size=orb_blk_sizes, &
                        col_blk_size=iao_blk_sizes, nze=0)
      CALL dbcsr_create(sx_vec, name="SX_VEC", dist=dbcsr_dist, &
                        matrix_type=dbcsr_type_no_symmetry, row_blk_size=oce_blk_sizes, &
                        col_blk_size=iao_blk_sizes, nze=0)
      CALL dbcsr_reserve_diag_blocks(matrix=sx_vec)

      nkind = SIZE(smat_kind)
      ALLOCATE (sinv_kind(nkind))
      DO ikind = 1, nkind
         noce = SIZE(smat_kind(ikind)%array, 1)
         ALLOCATE (sinv_kind(ikind)%array(noce, noce))
         sinv_kind(ikind)%array = smat_kind(ikind)%array
         CALL invmat_symm(sinv_kind(ikind)%array)
      END DO
      CALL dbcsr_get_info(iao_vec, group=para_group_handle)
      CALL para_type%set_handle(para_group_handle)

      nspin = SIZE(iao_coef, 1)
      ALLOCATE (oce_comp(natom, nspin))
      oce_comp = 0.0_dp
      DO ispin = 1, nspin
         CALL copy_fm_to_dbcsr(iao_coef(ispin), iao_vec)
         CALL dbcsr_multiply("N", "N", 1.0_dp, sxo(1)%matrix, iao_vec, 0.0_dp, sx_vec, &
                             retain_sparsity=.TRUE.)
         DO iatom = 1, natom
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL dbcsr_get_block_p(matrix=sx_vec, &
                                   row=iatom, col=iatom, BLOCK=ablock, found=found)
            IF (found) THEN
               n = SIZE(ablock, 1)
               m = SIZE(ablock, 2)
               ablock = MATMUL(sinv_kind(ikind)%array, ablock)
               ALLOCATE (amat(n, m), bmat(m, m), ev(m), vec(m, m))
               amat(1:n, 1:m) = MATMUL(smat_kind(ikind)%array(1:n, 1:n), ablock(1:n, 1:m))
               bmat(1:m, 1:m) = MATMUL(TRANSPOSE(ablock(1:n, 1:m)), amat(1:n, 1:m))
               CALL jacobi(bmat, ev, vec)
               oce_comp(iatom, ispin) = SUM(ev(1:m))/REAL(m, KIND=dp)
               DO i = 1, m
                  ev(i) = 1._dp/SQRT(ev(i))
                  bmat(1:m, i) = vec(1:m, i)*ev(i)
               END DO
               bmat(1:m, 1:m) = MATMUL(bmat(1:m, 1:m), TRANSPOSE(vec(1:m, 1:m)))
               ablock(1:n, 1:m) = MATMUL(ablock(1:n, 1:m), bmat(1:m, 1:m))
               DEALLOCATE (amat, bmat, ev, vec)
            END IF
         END DO
         CALL dbcsr_replicate_all(sx_vec)
         DO iatom = 1, natom
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL dbcsr_get_block_p(matrix=sx_vec, &
                                   row=iatom, col=iatom, BLOCK=ablock, found=found)
            CPASSERT(found)
            n = SIZE(ablock, 1)
            m = SIZE(ablock, 2)
            ALLOCATE (oce_atom(iatom, ispin)%array(n, m))
            oce_atom(iatom, ispin)%array = ablock
         END DO
      END DO

      CALL para_type%sum(oce_comp)
      IF (unit_nr > 0) THEN
         DO iatom = 1, natom
            WRITE (unit_nr, "(T4,A,I6,T30,A,2F12.4)") "Atom:", iatom, " Completeness: ", oce_comp(iatom, 1:nspin)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            oce_basis => oce_basis_set_list(ikind)%gto_basis_set
            CALL get_gto_basis_set(oce_basis, nset=nset, nshell=nshell, nsgf=nsgf, maxl=maxl, &
                                   l=lval, first_sgf=first_sgf, last_sgf=last_sgf)
            ALLOCATE (filter(nsgf, 0:maxl), vector(nsgf), prol(0:maxl, nspin))
            filter = 0.0_dp
            DO iset = 1, nset
               DO ishell = 1, nshell(iset)
                  l = lval(ishell, iset)
                  i1 = first_sgf(ishell, iset)
                  i2 = last_sgf(ishell, iset)
                  filter(i1:i2, l) = 1.0_dp
               END DO
            END DO
            !
            n = SIZE(oce_atom(iatom, 1)%array, 1)
            m = SIZE(oce_atom(iatom, 1)%array, 2)
            CPASSERT(n == nsgf)
            DO iao = 1, m
               prol = 0.0_dp
               DO ispin = 1, nspin
                  DO l = 0, maxl
                     vector(1:n) = oce_atom(iatom, ispin)%array(1:n, iao)*filter(1:n, l)
                     prol(l, ispin) = SUM(MATMUL(smat_kind(ikind)%array(1:n, 1:n), vector(1:n))*vector(1:n))
                  END DO
               END DO
               WRITE (unit_nr, "(T4,I3,T15,A,T39,(6F7.4))") iao, " l-contributions:", (SUM(prol(l, :)), l=0, maxl)
            END DO
            DEALLOCATE (filter, vector, prol)
         END DO
         WRITE (unit_nr, *)
      END IF

      DEALLOCATE (oce_comp)
      DO ikind = 1, nkind
         DEALLOCATE (sinv_kind(ikind)%array)
      END DO
      DEALLOCATE (sinv_kind)
      DEALLOCATE (iao_blk_sizes)
      CALL dbcsr_release(iao_vec)
      CALL dbcsr_release(sx_vec)

      CALL timestop(handle)

   END SUBROUTINE iao_oce_expansion

! **************************************************************************************************
!> \brief ...
!> \param smat_kind ...
!> \param basis_set_list ...
! **************************************************************************************************
   SUBROUTINE oce_overlap_matrix(smat_kind, basis_set_list)
      TYPE(cp_2d_r_p_type), DIMENSION(:)                 :: smat_kind
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list

      CHARACTER(len=*), PARAMETER :: routineN = 'oce_overlap_matrix'

      INTEGER                                            :: handle, ikind, iset, jset, ldsab, m1, &
                                                            m2, n1, n2, ncoa, ncob, nkind, nseta, &
                                                            sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: oint, owork
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, scon_a, smat, zeta
      TYPE(gto_basis_set_type), POINTER                  :: basis_set

      CALL timeset(routineN, handle)

      rab(1:3) = 0.0_dp

      nkind = SIZE(smat_kind)
      ldsab = 0
      DO ikind = 1, nkind
         basis_set => basis_set_list(ikind)%gto_basis_set
         CPASSERT(ASSOCIATED(basis_set))
         CALL get_gto_basis_set(gto_basis_set=basis_set, maxco=m1, nsgf=m2)
         ldsab = MAX(m1, m2, ldsab)
      END DO

      ALLOCATE (oint(ldsab, ldsab), owork(ldsab, ldsab))

      DO ikind = 1, nkind
         basis_set => basis_set_list(ikind)%gto_basis_set
         CPASSERT(ASSOCIATED(basis_set))
         smat => smat_kind(ikind)%array
         smat = 0.0_dp
         ! basis ikind
         first_sgfa => basis_set%first_sgf
         la_max => basis_set%lmax
         la_min => basis_set%lmin
         npgfa => basis_set%npgf
         nseta = basis_set%nset
         nsgfa => basis_set%nsgf_set
         rpgfa => basis_set%pgf_radius
         scon_a => basis_set%scon
         zeta => basis_set%zet

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            n1 = npgfa(iset)*(ncoset(la_max(iset)) - ncoset(la_min(iset) - 1))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nseta

               ncob = npgfa(jset)*ncoset(la_max(jset))
               n2 = npgfa(jset)*(ncoset(la_max(jset)) - ncoset(la_min(jset) - 1))
               sgfb = first_sgfa(1, jset)

               ! calculate integrals and derivatives
               CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                               la_max(jset), la_min(jset), npgfa(jset), rpgfa(:, jset), zeta(:, jset), &
                               rab, sab=oint(:, :))
               ! Contraction
               CALL contraction(oint(:, :), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgfa(iset), &
                                cb=scon_a(:, sgfb:), nb=n2, mb=nsgfa(jset), fscale=1.0_dp, trans=.FALSE.)
               CALL block_add("IN", owork, nsgfa(iset), nsgfa(jset), smat, &
                              sgfa, sgfb, trans=.FALSE.)

            END DO
         END DO

      END DO
      DEALLOCATE (oint, owork)

      CALL timestop(handle)

   END SUBROUTINE oce_overlap_matrix

! **************************************************************************************************

END MODULE iao_analysis
