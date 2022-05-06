!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE optbas_opt_utils
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert,&
                                              cp_fm_trace
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_distribution_type,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_transposed,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE qs_condnum,                      ONLY: overlap_condnum
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: evaluate_optvals, fit_mo_coeffs, optbas_build_neighborlist

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optbas_opt_utils'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mos ...
!> \param mos_aux_fit ...
!> \param matrix_ks ...
!> \param Q ...
!> \param Snew ...
!> \param S_inv_orb ...
!> \param fval ...
!> \param energy ...
!> \param S_cond_number ...
! **************************************************************************************************
   SUBROUTINE evaluate_optvals(mos, mos_aux_fit, matrix_ks, Q, Snew, S_inv_orb, &
                               fval, energy, S_cond_number)
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos, mos_aux_fit
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dbcsr_type), POINTER                          :: Q, Snew
      TYPE(cp_fm_type), POINTER                          :: S_inv_orb
      REAL(KIND=dp)                                      :: fval, energy, S_cond_number

      CHARACTER(len=*), PARAMETER                        :: routineN = 'evaluate_optvals'

      INTEGER                                            :: handle, ispin, iunit, naux, nmo, norb, &
                                                            nspins
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      REAL(KIND=dp)                                      :: tmp_energy, trace
      REAL(KIND=dp), DIMENSION(2)                        :: condnum
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_aux_fit, tmp1, tmp2
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: smat
      TYPE(dbcsr_type)                                   :: Qt

      CALL timeset(routineN, handle)

      nspins = SIZE(mos)

      NULLIFY (col_blk_sizes, row_blk_sizes)
      CALL dbcsr_get_info(Q, distribution=dbcsr_dist, &
                          nfullrows_total=naux, nfullcols_total=norb, &
                          row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes)
      CALL dbcsr_create(matrix=Qt, name="Qt", &
                        dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=col_blk_sizes, col_blk_size=row_blk_sizes, &
                        nze=0)
      CALL dbcsr_transposed(Qt, Q)
      !
      fval = 0.0_dp
      energy = 0.0_dp
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
         CALL get_mo_set(mos_aux_fit(ispin)%mo_set, mo_coeff=mo_coeff_aux_fit)
         CALL cp_fm_get_info(mo_coeff, ncol_global=nmo)
         CALL cp_fm_struct_create(fm_struct, ncol_global=nmo, &
                                  template_fmstruct=S_inv_orb%matrix_struct)
         CALL cp_fm_create(tmp1, matrix_struct=fm_struct)
         CALL cp_dbcsr_sm_fm_multiply(Qt, mo_coeff_aux_fit, tmp1, nmo)
         CALL cp_fm_trace(tmp1, mo_coeff, trace)
         fval = fval - 2.0_dp*trace + 2.0_dp*nmo
         !
         CALL cp_fm_create(tmp2, matrix_struct=fm_struct)
         CALL cp_fm_struct_release(fm_struct)
         CALL cp_gemm('N', 'N', norb, nmo, norb, 1.0_dp, S_inv_orb, tmp1, 0.0_dp, tmp2)
         CALL cp_dbcsr_sm_fm_multiply(matrix_ks(ispin)%matrix, tmp2, tmp1, nmo)
         CALL cp_fm_trace(tmp2, tmp1, tmp_energy)
         energy = energy + tmp_energy*(3.0_dp - REAL(nspins, KIND=dp))
         CALL cp_fm_release(tmp1)
         CALL cp_fm_release(tmp2)
      END DO
      CALL dbcsr_release(Qt)

      ALLOCATE (smat(1, 1))
      smat(1, 1)%matrix => Snew
      iunit = -1
      CALL cp_fm_get_info(S_inv_orb, context=blacs_env)
      CALL overlap_condnum(smat, condnum, iunit, .FALSE., .TRUE., .FALSE., blacs_env)
      S_cond_number = condnum(2)
      DEALLOCATE (smat)

      CALL timestop(handle)

   END SUBROUTINE evaluate_optvals

! **************************************************************************************************
!> \brief ...
!> \param saux ...
!> \param sauxorb ...
!> \param mos ...
!> \param mosaux ...
! **************************************************************************************************
   SUBROUTINE fit_mo_coeffs(saux, sauxorb, mos, mosaux)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: saux, sauxorb
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos, mosaux

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fit_mo_coeffs'
      REAL(KIND=dp), PARAMETER                           :: threshold = 1.E-12_dp

      INTEGER                                            :: handle, ispin, naux, ndep, nmo, norb, &
                                                            nspins
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_s, fm_sinv, mo_coeff, mo_coeff_aux, &
                                                            tmat, tmp1, tmp2, work

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(saux(1)%matrix, nfullrows_total=naux)
      CALL dbcsr_get_info(sauxorb(1)%matrix, nfullcols_total=norb)
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff)

      CALL cp_fm_struct_create(fm_struct, nrow_global=naux, ncol_global=naux, &
                               context=mo_coeff%matrix_struct%context, &
                               para_env=mo_coeff%matrix_struct%para_env)
      CALL cp_fm_create(fm_s, fm_struct, name="s_aux")
      CALL cp_fm_create(fm_sinv, fm_struct, name="s_aux_inv")
      CALL copy_dbcsr_to_fm(saux(1)%matrix, fm_s)
      CALL cp_fm_invert(fm_s, fm_sinv)
      CALL cp_fm_release(fm_s)
      CALL cp_fm_struct_release(fm_struct)
      nspins = SIZE(mos)
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
         CALL get_mo_set(mosaux(ispin)%mo_set, mo_coeff=mo_coeff_aux)
         CALL cp_fm_get_info(mo_coeff, ncol_global=nmo)
         CALL cp_fm_create(tmp1, matrix_struct=mo_coeff_aux%matrix_struct)
         CALL cp_fm_create(tmp2, matrix_struct=mo_coeff_aux%matrix_struct)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nmo, ncol_global=nmo, &
                                  context=mo_coeff%matrix_struct%context, &
                                  para_env=mo_coeff%matrix_struct%para_env)
         CALL cp_fm_create(tmat, fm_struct, name="tmat")
         CALL cp_fm_create(work, fm_struct, name="work")
         CALL cp_fm_struct_release(fm_struct)
         !
         CALL cp_dbcsr_sm_fm_multiply(sauxorb(1)%matrix, mo_coeff, tmp1, nmo)
         CALL cp_gemm('N', 'N', naux, nmo, naux, 1.0_dp, fm_sinv, tmp1, 0.0_dp, tmp2)
         CALL cp_gemm('T', 'N', nmo, nmo, naux, 1.0_dp, tmp1, tmp2, 0.0_dp, tmat)
         CALL cp_fm_power(tmat, work, -0.5_dp, threshold, ndep)
         CALL cp_gemm('N', 'N', naux, nmo, nmo, 1.0_dp, tmp2, tmat, 0.0_dp, mo_coeff_aux)
         !
         CALL cp_fm_release(work)
         CALL cp_fm_release(tmat)
         CALL cp_fm_release(tmp1)
         CALL cp_fm_release(tmp2)
      END DO
      CALL cp_fm_release(fm_sinv)

      CALL timestop(handle)

   END SUBROUTINE fit_mo_coeffs

! **************************************************************************************************
!> \brief rebuilds neighborlist for absis sets
!> \param qs_env ...
!> \param sab_aux ...
!> \param sab_aux_orb ...
!> \param basis_type ...
!> \par History
!>       adapted from kg_build_neighborlist
! **************************************************************************************************
   SUBROUTINE optbas_build_neighborlist(qs_env, sab_aux, sab_aux_orb, basis_type)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux, sab_aux_orb
      CHARACTER(*)                                       :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'optbas_build_neighborlist'

      INTEGER                                            :: handle, ikind, nkind
      LOGICAL                                            :: mic, molecule_only
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: aux_fit_present, orb_present
      REAL(dp)                                           :: subcells
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: aux_fit_radius, orb_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(gto_basis_set_type), POINTER                  :: aux_fit_basis_set, orb_basis_set
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)
      NULLIFY (para_env)

      ! restrict lists to molecular subgroups
      molecule_only = .FALSE.
      mic = molecule_only

      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      distribution_2d=distribution_2d, &
                      molecule_set=molecule_set, &
                      local_particles=distribution_1d, &
                      particle_set=particle_set, &
                      para_env=para_env)

      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)

      ! Allocate work storage
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (orb_radius(nkind), aux_fit_radius(nkind))
      orb_radius(:) = 0.0_dp
      aux_fit_radius(:) = 0.0_dp
      ALLOCATE (orb_present(nkind), aux_fit_present(nkind))
      ALLOCATE (pair_radius(nkind, nkind))
      ALLOCATE (atom2d(nkind))

      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, molecule_only, particle_set=particle_set)

      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom2d(ikind)%list)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type="ORB")
         IF (ASSOCIATED(orb_basis_set)) THEN
            orb_present(ikind) = .TRUE.
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, kind_radius=orb_radius(ikind))
         ELSE
            orb_present(ikind) = .FALSE.
            orb_radius(ikind) = 0.0_dp
         END IF
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=aux_fit_basis_set, basis_type=basis_type)
         IF (ASSOCIATED(aux_fit_basis_set)) THEN
            aux_fit_present(ikind) = .TRUE.
            CALL get_gto_basis_set(gto_basis_set=aux_fit_basis_set, kind_radius=aux_fit_radius(ikind))
         ELSE
            aux_fit_present(ikind) = .FALSE.
            aux_fit_radius(ikind) = 0.0_dp
         END IF
      END DO
      !
      CALL pair_radius_setup(aux_fit_present, aux_fit_present, aux_fit_radius, aux_fit_radius, pair_radius)
      CALL build_neighbor_lists(sab_aux, particle_set, atom2d, cell, pair_radius, &
                                mic=mic, molecular=molecule_only, subcells=subcells, nlname="sab_aux")
      CALL pair_radius_setup(aux_fit_present, orb_present, aux_fit_radius, orb_radius, pair_radius)
      CALL build_neighbor_lists(sab_aux_orb, particle_set, atom2d, cell, pair_radius, &
                                mic=mic, symmetric=.FALSE., molecular=molecule_only, subcells=subcells, &
                                nlname="sab_aux_orb")

      ! Release work storage
      CALL atom2d_cleanup(atom2d)
      DEALLOCATE (atom2d)
      DEALLOCATE (orb_present, aux_fit_present)
      DEALLOCATE (orb_radius, aux_fit_radius)
      DEALLOCATE (pair_radius)

      CALL timestop(handle)

   END SUBROUTINE optbas_build_neighborlist

END MODULE optbas_opt_utils
