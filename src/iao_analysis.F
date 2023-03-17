!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate ntrinsic atomic orbitals and analyze wavefunctions
!> \par History
!>      03.2023 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE iao_analysis
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_container_types,       ONLY: add_basis_set_to_container,&
                                              basis_set_container_type,&
                                              get_basis_from_container,&
                                              remove_basis_from_container,&
                                              remove_basis_set_container
   USE basis_set_types,                 ONLY: allocate_sto_basis_set,&
                                              create_gto_from_sto_basis,&
                                              deallocate_sto_basis_set,&
                                              get_gto_basis_set,&
                                              gto_basis_set_type,&
                                              gto_basis_set_p_type,&
                                              init_orb_basis_set,&
                                              set_sto_basis_set,&
                                              srules,&
                                              sto_basis_set_type
   USE bibliography,                    ONLY: Knizia2013,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              copy_fm_to_dbcsr,&
                                              copy_dbcsr_to_fm
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_copy_general, cp_fm_create, cp_fm_get_info, cp_fm_release, cp_fm_set_all, &
        cp_fm_set_element, cp_fm_to_fm, cp_fm_to_fm_submat_general, cp_fm_type, cp_fm_release
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE dbcsr_api,                       ONLY: dbcsr_get_block_p,&
                                              dbcsr_distribution_type,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_multiply,&
                                              dbcsr_p_type,&
                                              dbcsr_create,&
                                              dbcsr_release,&
                                              dbcsr_get_info,&
                                              dbcsr_type
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_comm_type,&
                                              mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: get_ptable_info,&
                                              ptable
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
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_mo_types,                     ONLY: mo_set_type, get_mo_set
   USE qs_mo_methods,                   ONLY: make_basis_lowdin
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

   PUBLIC ::  iao_wfn_analysis

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param input_section ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE iao_wfn_analysis(qs_env, input_section, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_section
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER                        :: routineN = 'iao_wfn_analysis'

      CHARACTER(LEN=default_string_length)               :: bname
      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: handle, ikind, ispin, nimages, nkind, &
                                                            nspin, z, nsgf, nref, norb
      INTEGER, DIMENSION(:), POINTER                     :: econf
      LOGICAL                                            :: do_bondorder, explicit, do_charges, do_printbas
      REAL(KIND=dp)                                      :: eps_filter, eps_svd, zval
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: p_ref_orb, p_orb_ref
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: iao_coef
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dbcsr_type), POINTER                          :: smat
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: sref, sro 
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(gto_basis_set_type), POINTER                  :: refbasis, orbbasis
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: ref_basis_set_list, orb_basis_set_list
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: srr_list, sro_list
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      ! only do IAO analysis if explicitly requested
      CALL section_vals_get(input_section, explicit=explicit)
      IF (.NOT. explicit) RETURN

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

      ! input options
      CALL section_vals_val_get(input_section, "EPS_FILTER", r_val=eps_filter)
      CALL section_vals_val_get(input_section, "EPS_SVD", r_val=eps_svd)
      CALL section_vals_val_get(input_section, "ATOMIC_CHARGES", l_val=do_charges)
      CALL section_vals_val_get(input_section, "BOND_ORDER", l_val=do_bondorder)
      CALL section_vals_val_get(input_section, "PRINT_BASIS", l_val=do_printbas)

      ! check for or generate reference basis
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind, zeff=zval, elec_conf=econf, element_symbol=element_symbol)
         CALL get_ptable_info(element_symbol, ielement=z)
         NULLIFY (refbasis)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=refbasis, basis_type="MIN")
         IF (.NOT. ASSOCIATED(refbasis)) THEN
            ! generate a minimal basis set
            CALL create_min_basis(refbasis, z, econf)
            CALL add_basis_set_to_container(qs_kind%basis_sets, refbasis, "MIN")
         END IF
         IF (unit_nr > 0) THEN
            CALL get_gto_basis_set(refbasis,name=bname,nsgf=nsgf)
            WRITE (unit_nr, '(T2,A,A,T14,A,I4,T40,A,A24)') "Kind: ",element_symbol,"NBasFun: ",nsgf,&
                "Reference Basis: ",ADJUSTL(TRIM(bname))
         END IF
      END DO

      ! overlap matrices
      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s)
      smat => matrix_s(1, 1)%matrix
      !
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
      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mos(1), mo_coeff=mo_coeff)
      CALL dbcsr_get_info(sro(1)%matrix, nfullrows_total=nref, nfullcols_total=norb)
      !
      CALL cp_fm_struct_create(fm_struct, context=mo_coeff%matrix_struct%context, nrow_global=nref, &
                               ncol_global=norb, para_env=mo_coeff%matrix_struct%para_env)
      CALL cp_fm_create(p_ref_orb, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=mo_coeff%matrix_struct%context, nrow_global=norb, &
                               ncol_global=nref, para_env=mo_coeff%matrix_struct%para_env)
      CALL cp_fm_create(p_orb_ref, fm_struct)
      !
      CALL iao_projectors(smat, sref(1)%matrix, sro(1)%matrix, &
                          p_orb_ref, p_ref_orb, eps_svd, eps_filter)

      ! generate IAOs
      ALLOCATE (iao_coef(nspin))
      DO ispin = 1, nspin
!        CALL intrinsic_ao_calc(smat, p_orb_ref, p_ref_orb, cvec, avec)
      END DO
!        ! Mulliken population
!        IF (unit_nr > 0) THEN
!           WRITE (unit_nr, '(/,T2,A)') 'Mulliken Population Analysis '
!        END IF
!        ALLOCATE (mcharge(natom, nspin))
!        DO ispin = 1, nspin
!           CALL dbcsr_dot(pqmat(ispin)%matrix, sqmat(ispin)%matrix, trace)
!           IF (unit_nr > 0) THEN
!              WRITE (unit_nr, '(T2,A,I2,T66,F15.4)') 'Number of Electrons: Trace(PS) Spin ', ispin, trace
!           END IF
!           CALL mulliken_charges(pqmat(ispin)%matrix, sqmat(ispin)%matrix, para_env, mcharge(:, ispin))
!        END DO
!        CALL print_atomic_charges(particle_set, qs_kind_set, unit_nr, "Minimal Basis Mulliken Charges", &
!                                  electronic_charges=mcharge)
!        ! Mayer bond orders
!        IF (do_bondorder) THEN
!           ALLOCATE (border(natom, natom))
!           border = 0.0_dp
!           CALL dbcsr_create(psmat, name="PS", template=sqmat(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
!           CALL dbcsr_create(spmat, name="SP", template=sqmat(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
!           filter_eps = 1.e-6_dp
!           DO ispin = 1, nspin
!              CALL dbcsr_multiply("N", "N", 1.0_dp, pqmat(ispin)%matrix, sqmat(ispin)%matrix, 0.0_dp, psmat, &
!                                  filter_eps=filter_eps)
!              CALL dbcsr_multiply("N", "N", 1.0_dp, sqmat(ispin)%matrix, pqmat(ispin)%matrix, 0.0_dp, spmat, &
!                                  filter_eps=filter_eps)
!              CALL compute_bond_order(psmat, spmat, border)
!           END DO
!           CALL para_env%sum(border)
!           border = border*REAL(nspin, KIND=dp)
!           CALL dbcsr_release(psmat)
!           CALL dbcsr_release(spmat)
!           CALL print_bond_orders(particle_set, unit_nr, border)
!           DEALLOCATE (border)
!        END IF

!        ! for printing purposes we now copy the QUAMBOs into MO format
!        ALLOCATE (mbas(nspin))
!        DO ispin = 1, nspin
!           CALL allocate_mo_set(mbas(ispin), nao, nmao, nmao, 0.0_dp, 1.0_dp, 0.0_dp)
!           CALL set_mo_set(mbas(ispin), homo=nmao)
!           ALLOCATE (mbas(ispin)%eigenvalues(nmao))
!           mbas(ispin)%eigenvalues = 0.0_dp
!           ALLOCATE (mbas(ispin)%occupation_numbers(nmao))
!           mbas(ispin)%occupation_numbers = 1.0_dp
!           CALL cp_fm_create(mbas(ispin)%mo_coeff, fm_struct_a)
!           CALL copy_dbcsr_to_fm(quambo(ispin)%matrix, mbas(ispin)%mo_coeff)
!        END DO

!        ! Print basis functions: cube files
!        DO ispin = 1, nspin
!           CALL get_mo_set(mbas(ispin), mo_coeff=fm_mos)
!           CALL post_minbas_cubes(qs_env, input_section, fm_mos, ispin)
!        END DO
!        ! Print basis functions: molden format
!        molden_section => section_vals_get_subs_vals(input_section, "MINBAS_MOLDEN")
!        CALL write_mos_molden(mbas, qs_kind_set, particle_set, molden_section)
!        DO ispin = 1, nspin
!           CALL deallocate_mo_set(mbas(ispin))
!        END DO
!        DEALLOCATE (mbas)

      ! Deallocate the neighbor list structure
      CALL release_neighbor_list_sets(srr_list)
      CALL release_neighbor_list_sets(sro_list)

      DEALLOCATE (ref_basis_set_list, orb_basis_set_list)

      IF (ASSOCIATED(sref)) CALL dbcsr_deallocate_matrix_set(sref)
      IF (ASSOCIATED(sro)) CALL dbcsr_deallocate_matrix_set(sro)
      CALL cp_fm_release(p_ref_orb)
      CALL cp_fm_release(p_orb_ref)

      DO ispin = 1, nspin
         CALL cp_fm_release(iao_coef(ispin))
      END DO
      DEALLOCATE (iao_coef)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/,T2,A)') &
            '!----------------------------END OF IAO ANALYSIS------------------------------!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE iao_wfn_analysis

! **************************************************************************************************
!> \brief Computes projector matrices for ref to orb basis and reverse
! **************************************************************************************************
   SUBROUTINE iao_projectors(smat, sref, s_r_o, p_o_r, p_r_o, eps_svd, eps_filter)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: smat, sref
      TYPE(dbcsr_type), INTENT(INOUT)                    :: s_r_o
      TYPE(cp_fm_type), INTENT(INOUT)                    :: p_o_r, p_r_o
      REAL(KIND=dp), INTENT(IN)                          :: eps_svd, eps_filter

      CHARACTER(len=*), PARAMETER                        :: routineN = 'iao_projectors'

      INTEGER                                            :: handle, norb, nref
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_sro, fm_s, fm_inv

      MARK_USED(eps_filter)

      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(p_r_o, nrow_global=nref, ncol_global=norb)
      CALL cp_fm_create(fm_sro, p_r_o%matrix_struct)
      CALL copy_dbcsr_to_fm(s_r_o, fm_sro)

      CALL cp_fm_struct_create(fm_struct, nrow_global=norb, ncol_global=norb, &
                               template_fmstruct=p_r_o%matrix_struct)
      CALL cp_fm_create(fm_s, fm_struct, name="smat")
      CALL cp_fm_create(fm_inv, fm_struct, name="sinv")
      CALL cp_fm_struct_release(fm_struct)
      CALL copy_dbcsr_to_fm(smat, fm_s)
      CALL cp_fm_invert(fm_s, fm_inv, eps_svd=eps_svd)
      CALL parallel_gemm('N', 'T', norb, norb, nref, 1.0_dp, fm_inv, fm_sro, 0.0_dp, p_o_r)
      CALL cp_fm_release(fm_s)
      CALL cp_fm_release(fm_inv)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nref, ncol_global=nref, &
                               template_fmstruct=p_r_o%matrix_struct)
      CALL cp_fm_create(fm_s, fm_struct, name="sref")
      CALL cp_fm_create(fm_inv, fm_struct, name="sinv")
      CALL cp_fm_struct_release(fm_struct)
      CALL copy_dbcsr_to_fm(sref, fm_s)
      CALL cp_fm_invert(fm_s, fm_inv, eps_svd=eps_svd)
      CALL parallel_gemm('N', 'N', nref, nref, norb, 1.0_dp, fm_inv, fm_sro, 0.0_dp, p_r_o)
      CALL cp_fm_release(fm_s)
      CALL cp_fm_release(fm_inv)

      CALL cp_fm_release(fm_sro)

      CALL timestop(handle)

   END SUBROUTINE iao_projectors

! **************************************************************************************************
!> \brief Computes intrinsic orbitals for a given MO vector set
! **************************************************************************************************
   SUBROUTINE intrinsic_ao_calc(smat, p_o_r, p_r_o, cvec, avec)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: smat
      TYPE(cp_fm_type), INTENT(INOUT)                    :: p_r_o, p_o_r
      TYPE(cp_fm_type), INTENT(INOUT)                    :: cvec, avec

      CHARACTER(len=*), PARAMETER                        :: routineN = 'intrinsic_ao_calc'

      INTEGER                                            :: handle, norb, nao, nref
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: sc, sct, pc, vec1, ctvec

      CALL timeset(routineN, handle)

      ! number of orbitals
      CALL cp_fm_get_info(cvec, ncol_global=norb)
      ! basis set sizes
      CALL cp_fm_get_info(p_r_o, nrow_global=nref, ncol_global=nao)
      ! temp array
      CALL cp_fm_struct_create(fm_struct, nrow_global=nref, ncol_global=norb, &
                               template_fmstruct=cvec%matrix_struct)
      CALL cp_fm_create(pc, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      ! CT = orth(Por.Pro.C)
      CALL parallel_gemm('N', 'N', nref, nao, norb, 1.0_dp, p_r_o, cvec, 0.0_dp, pc)
      CALL cp_fm_create(ctvec, cvec%matrix_struct)
      CALL parallel_gemm('N', 'N', nao, nref, norb, 1.0_dp, p_o_r, pc, 0.0_dp, ctvec)
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
      CALL parallel_gemm('T', 'N', norb, nao, nref, 1.0_dp, sct, p_r_o, 0.0_dp, pc)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nref, &
                               template_fmstruct=cvec%matrix_struct)
      CALL cp_fm_create(vec1, cvec%matrix_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL parallel_gemm('N', 'N', nao, norb, nref, 1.0_dp, ctvec, pc, 0.0_dp, vec1)
      ! A = C*SC(T)*V1
      CALL parallel_gemm('T', 'N', norb, nao, nref, 1.0_dp, sc, vec1, 0.0_dp, pc)
      CALL parallel_gemm('N', 'N', nao, norb, nref, 1.0_dp, cvec, pc, 0.0_dp, avec)
      ! V1 = Por - V1
      CALL cp_fm_scale_and_add(-1.0_dp, vec1, 1.0_dp, p_r_o)
      ! A = A + V1
      CALL cp_fm_scale_and_add(1.0_dp, avec, 1.0_dp, vec1)
      ! A = A - C*SC(T)*V1
      CALL parallel_gemm('T', 'N', norb, nao, nref, 1.0_dp, sc, vec1, 0.0_dp, pc)
      CALL parallel_gemm('N', 'N', nao, norb, nref, -1.0_dp, cvec, pc, 1.0_dp, avec)

      ! clean up
      CALL cp_fm_release(pc)
      CALL cp_fm_release(vec1)
      CALL cp_fm_release(sc)
      CALL cp_fm_release(sct)

      CALL timestop(handle)

   END SUBROUTINE intrinsic_ao_calc

! **************************************************************************************************
!> \brief Computes and prints the Cube Files for the intrinsic atomic orbitals
!> \param qs_env ...
!> \param print_section ...
!> \param iao_coeff ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE print_iao_cubes(qs_env, print_section, iao_coeff, ispin)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: print_section
      TYPE(cp_fm_type), INTENT(IN)                       :: iao_coeff
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(LEN=default_path_length)                 :: filename, title
      INTEGER                                            :: i, i_rep, ivec, iw, j, n_rep, natom
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes, first_bas, ilist, stride
      LOGICAL                                            :: explicit, mpi_io
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
      TYPE(section_vals_type), POINTER                   :: iao_section

      iao_section => section_vals_get_subs_vals(print_section, "IAO_CUBE")
      CALL section_vals_get(iao_section, explicit=explicit)
      IF (.NOT. explicit) RETURN

      logger => cp_get_default_logger()
      stride => section_get_ivals(print_section, "IAO_CUBE%STRIDE")

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      subsys=subsys, cell=cell, particle_set=particle_set, pw_env=pw_env, dft_control=dft_control)
      CALL qs_subsys_get(subsys, particles=particles)

      CALL get_qs_env(qs_env=qs_env, natom=natom)
      ALLOCATE (blk_sizes(natom), first_bas(0:natom))
      CALL get_particle_set(particle_set, qs_kind_set, nmao=blk_sizes)
      first_bas(0) = 0
      DO i = 1, natom
         first_bas(i) = first_bas(i - 1) + blk_sizes(i)
      END DO

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_r, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_g, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      ! loop over list of atoms
      CALL section_vals_val_get(iao_section, "ATOM_LIST", n_rep_val=n_rep)
      IF (n_rep == 0) THEN
         DO i = 1, natom
            DO ivec = first_bas(i - 1) + 1, first_bas(i)
               WRITE (filename, '(a4,I5.5,a1,I1.1)') "IAO_", ivec, "_", ispin
               WRITE (title, *) "Intrinsic Atomic Orbitals ", ivec, " atom ", i, " spin ", ispin
               mpi_io = .TRUE.
               iw = cp_print_key_unit_nr(logger, print_section, "IAO_CUBE", extension=".cube", &
                                         middle_name=TRIM(filename), file_position="REWIND", log_filename=.FALSE., &
                                         mpi_io=mpi_io)
               CALL calculate_wavefunction(iao_coeff, ivec, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                           cell, dft_control, particle_set, pw_env)
               CALL cp_pw_to_cube(wf_r, iw, title, particles=particles, stride=stride, mpi_io=mpi_io)
               CALL cp_print_key_finished_output(iw, logger, print_section, "IAO_CUBE", mpi_io=mpi_io)
            END DO
         END DO
      ELSE
         DO i_rep = 1, n_rep
            CALL section_vals_val_get(iao_section, "ATOM_LIST", i_rep_val=i_rep, i_vals=ilist)
            DO i = 1, SIZE(ilist, 1)
               j = ilist(i)
               DO ivec = first_bas(j - 1) + 1, first_bas(j)
                  WRITE (filename, '(a4,I5.5,a1,I1.1)') "IAO_", ivec, "_", ispin
                  WRITE (title, *) "Intrinsic Atomic Orbitals ", ivec, " atom ", j, " spin ", ispin
                  mpi_io = .TRUE.
                  iw = cp_print_key_unit_nr(logger, print_section, "IAO_CUBE", extension=".cube", &
                                            middle_name=TRIM(filename), file_position="REWIND", log_filename=.FALSE., &
                                            mpi_io=mpi_io)
                  CALL calculate_wavefunction(iao_coeff, ivec, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                              cell, dft_control, particle_set, pw_env)
                  CALL cp_pw_to_cube(wf_r, iw, title, particles=particles, stride=stride, mpi_io=mpi_io)
                  CALL cp_print_key_finished_output(iw, logger, print_section, "IAO_CUBE", mpi_io=mpi_io)
               END DO
            END DO
         END DO
      END IF
      DEALLOCATE (blk_sizes, first_bas)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g)

   END SUBROUTINE print_iao_cubes

! **************************************************************************************************
!> \brief Creates a minimal basis set based on Slater rules
!> \param min_basis_set ...
!> \param zval ...
!> \param ishell ...
!> \param nshell ...
!> \param lnam ...
!> \par History
!>      03.2023 created [JGH]
! **************************************************************************************************
   SUBROUTINE create_min_basis(min_basis_set, zval, econf)
      TYPE(gto_basis_set_type), POINTER                  :: min_basis_set
      INTEGER, INTENT(IN)                                :: zval
      INTEGER, DIMENSION(0:3)                            :: econf

      CHARACTER(LEN=default_string_length)               :: bname
      CHARACTER(len=6)                                   :: str
      CHARACTER(len=6), DIMENSION(:), POINTER            :: sym
      INTEGER, DIMENSION(:), POINTER                     :: lq, nq
      REAL(KIND=dp), DIMENSION(:), POINTER               :: zet

      INTEGER                                            :: i, l, nss, iss, n, nn
      INTEGER, DIMENSION(0:3)                            :: nae, npe
      INTEGER, DIMENSION(4, 7)                           :: ne
      CHARACTER(len=1), DIMENSION(0:3), PARAMETER        :: lnam = (/"S","P","D","F"/)
      TYPE(sto_basis_set_type), POINTER                  :: sto_basis_set

      CPASSERT(.NOT. ASSOCIATED(min_basis_set))
      NULLIFY (sto_basis_set)

      ! electronic configuration
      ne = 0
      DO l = 1, 4 
         nn = 2*(l - 1) + 1
         DO i = l, 7
            ne(l, i) = ptable(zval)%e_conv(l - 1) - 2*nn*(i - l)
            ne(l, i) = MAX(ne(l, i), 0)
            ne(l, i) = MIN(ne(l, i), 2*nn)
         END DO
      END DO
      ! STO definition
      nae = 0
      npe = 0
      DO l=0,3
         nn = 2*(2*l+1)
         nae(l) = ptable(zval)%e_conv(l)/nn
         IF(MOD(ptable(zval)%e_conv(l),nn) /= 0) nae(l) = nae(l)+1
         npe(l) = econf(l)/nn
         IF(MOD(econf(l),nn) /= 0) npe(l) = npe(l)+1
      END DO
      CPASSERT(ALL(nae-npe >= 0))
      nss = SUM(npe)
      ALLOCATE(sym(nss), lq(nss), nq(nss), zet(nss))
      iss = 0
      DO l=0,3
         DO i=1,npe(l)
            iss = iss + 1
            lq(iss) = l
            n = nae(l) - npe(l) + l
            nq(iss) = n + i
            str = "      "
            WRITE(str(5:5),FMT='(I1)') nq(iss)
            str(6:6) = lnam(l)
            sym(iss) = str
            zet(iss) = srules(zval, ne, nq(iss), lq(iss))
         END DO
      END DO

      bname = ADJUSTR(ptable(zval)%symbol)//"_MBas"
      CALL allocate_sto_basis_set(sto_basis_set)
      CALL set_sto_basis_set(sto_basis_set, name=bname, nshell=nss, nq=nq, &
                             lq=lq, zet=zet, symbol=sym)
      CALL create_gto_from_sto_basis(sto_basis_set, min_basis_set, 6)
      min_basis_set%norm_type = 2
      CALL init_orb_basis_set(min_basis_set)
      CALL deallocate_sto_basis_set(sto_basis_set)

      DEALLOCATE(sym, lq, nq, zet)

   END SUBROUTINE create_min_basis

END MODULE iao_analysis
