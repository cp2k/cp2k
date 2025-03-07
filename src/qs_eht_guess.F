!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Generate an initial guess (dm and orb) from EHT calculation
! **************************************************************************************************
MODULE qs_eht_guess
   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE input_constants,                 ONLY: do_method_xtb
   USE input_section_types,             ONLY: section_vals_duplicate,&
                                              section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_energy_init,                  ONLY: qs_energies_init
   USE qs_environment,                  ONLY: qs_init
   USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_create,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix_simple
   USE xtb_ks_matrix,                   ONLY: build_xtb_ks_matrix
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_eht_guess'

   PUBLIC ::  calculate_eht_guess

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief EHT MO guess calclulation
!> \param qs_env ...
!> \param mo_array ...
! **************************************************************************************************
   SUBROUTINE calculate_eht_guess(qs_env, mo_array)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mo_array

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_eht_guess'

      INTEGER                                            :: handle, ispin, nao, nbas, neeht, neorb, &
                                                            nkind, nmo, nspins, zero
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigval
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: mstruct_ee, mstruct_oe, mstruct_oo
      TYPE(cp_fm_type)                                   :: fmksmat, fmorb, fmscr, fmsmat, fmvec, &
                                                            fmwork, sfull, sinv
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ksmat, matrix_s, matrix_t, smat
      TYPE(dbcsr_type)                                   :: tempmat, tmat
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list_a, basis_set_list_b
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_environment_type), POINTER                 :: eht_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: dft_section, eht_force_env_section, &
                                                            force_env_section, qs_section, &
                                                            subsys_section, xtb_section

      CALL timeset(routineN, handle)

      NULLIFY (subsys_section)
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      para_env=para_env, &
                      input=force_env_section, &
                      cp_subsys=cp_subsys)
      NULLIFY (eht_force_env_section)
      CALL section_vals_duplicate(force_env_section, eht_force_env_section)
      dft_section => section_vals_get_subs_vals(eht_force_env_section, "DFT")
      qs_section => section_vals_get_subs_vals(dft_section, "QS")
      CALL section_vals_val_set(qs_section, "METHOD", i_val=do_method_xtb)
      xtb_section => section_vals_get_subs_vals(qs_section, "xTB")
      zero = 0
      CALL section_vals_val_set(xtb_section, "GFN_TYPE", i_val=zero)
      !
      ALLOCATE (eht_env)
      CALL qs_env_create(eht_env)
      CALL qs_init(eht_env, para_env, cp_subsys=cp_subsys, &
                   force_env_section=eht_force_env_section, &
                   subsys_section=subsys_section, &
                   use_motion_section=.FALSE., silent=.TRUE.)
      !
      CALL get_qs_env(qs_env, nelectron_total=neorb)
      CALL get_qs_env(eht_env, nelectron_total=neeht)
      IF (neorb /= neeht) THEN
         CPWARN("EHT has different number of electrons than calculation method.")
         CPABORT("EHT Initial Guess")
      END IF
      !
      CALL qs_env_rebuild_pw_env(eht_env)
      CALL qs_energies_init(eht_env, calc_forces=.FALSE.)
      CALL build_xtb_ks_matrix(eht_env, .FALSE., .FALSE.)
      !
      CALL get_qs_env(eht_env, &
                      matrix_s=smat, matrix_ks=ksmat)
      nspins = SIZE(ksmat, 1)
      CALL get_qs_env(eht_env, para_env=para_env, blacs_env=blacs_env)
      CALL dbcsr_get_info(smat(1)%matrix, nfullrows_total=nao)
      CALL cp_fm_struct_create(fmstruct=mstruct_ee, context=blacs_env, &
                               nrow_global=nao, ncol_global=nao)
      CALL cp_fm_create(fmksmat, mstruct_ee)
      CALL cp_fm_create(fmsmat, mstruct_ee)
      CALL cp_fm_create(fmvec, mstruct_ee)
      CALL cp_fm_create(fmwork, mstruct_ee)
      ALLOCATE (eigenvalues(nao))

      ! DBCSR matrix
      CALL dbcsr_create(tempmat, template=smat(1)%matrix, matrix_type=dbcsr_type_no_symmetry)

      ! transfer to FM
      CALL dbcsr_desymmetrize(smat(1)%matrix, tempmat)
      CALL copy_dbcsr_to_fm(tempmat, fmsmat)

      !SINV of origianl basis
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nbas)
      CALL dbcsr_create(tmat, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL cp_fm_struct_create(fmstruct=mstruct_oo, context=blacs_env, &
                               nrow_global=nbas, ncol_global=nbas)
      CALL cp_fm_create(sfull, mstruct_oo)
      CALL cp_fm_create(sinv, mstruct_oo)
      CALL dbcsr_desymmetrize(matrix_s(1)%matrix, tmat)
      CALL copy_dbcsr_to_fm(tmat, sfull)
      CALL cp_fm_invert(sfull, sinv)
      CALL dbcsr_release(tmat)
      CALL cp_fm_release(sfull)
      !TMAT(bas1, bas2)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, sab_all=sab_nl, nkind=nkind)
      IF (.NOT. ASSOCIATED(sab_nl)) THEN
         CPWARN("Full neighborlist not available for this method. EHT initial guess not possible.")
         CPABORT("EHT Initial Guess")
      END IF
      ALLOCATE (basis_set_list_a(nkind), basis_set_list_b(nkind))
      CALL basis_set_list_setup(basis_set_list_a, "ORB", qs_kind_set)
      CALL get_qs_env(eht_env, qs_kind_set=qs_kind_set)
      CALL basis_set_list_setup(basis_set_list_b, "ORB", qs_kind_set)
      !
      NULLIFY (matrix_t)
      CALL build_overlap_matrix_simple(ks_env, matrix_t, &
                                       basis_set_list_a, basis_set_list_b, sab_nl)
      DEALLOCATE (basis_set_list_a, basis_set_list_b)

      ! KS matrix is not spin dependent!
      CALL dbcsr_desymmetrize(ksmat(1)%matrix, tempmat)
      CALL copy_dbcsr_to_fm(tempmat, fmksmat)
      ! diagonalize
      CALL cp_fm_geeig(fmksmat, fmsmat, fmvec, eigenvalues, fmwork)
      ! Sinv*T*d
      CALL cp_fm_struct_create(fmstruct=mstruct_oe, context=blacs_env, &
                               nrow_global=nbas, ncol_global=nao)
      CALL cp_fm_create(fmscr, mstruct_oe)
      CALL cp_fm_create(fmorb, mstruct_oe)
      CALL cp_dbcsr_sm_fm_multiply(matrix_t(1)%matrix, fmvec, fmscr, ncol=nao)
      CALL parallel_gemm('N', 'N', nbas, nao, nbas, 1.0_dp, sinv, fmscr, 0.0_dp, fmorb)
      !
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mo_array(ispin), mo_coeff=mo_coeff, nmo=nmo)
         CALL cp_fm_to_fm(fmorb, mo_coeff, nmo, 1, 1)
         NULLIFY (eigval)
         CALL get_mo_set(mo_set=mo_array(ispin), eigenvalues=eigval)
         IF (ASSOCIATED(eigval)) THEN
            eigval(1:nmo) = eigenvalues(1:nmo)
         END IF
      END DO
      CALL set_mo_occupation(mo_array, smear=qs_env%scf_control%smear)

      DEALLOCATE (eigenvalues)
      CALL dbcsr_release(tempmat)
      CALL dbcsr_deallocate_matrix_set(matrix_t)
      CALL cp_fm_release(fmksmat)
      CALL cp_fm_release(fmsmat)
      CALL cp_fm_release(fmvec)
      CALL cp_fm_release(fmwork)
      CALL cp_fm_release(fmscr)
      CALL cp_fm_release(fmorb)
      CALL cp_fm_release(sinv)
      CALL cp_fm_struct_release(mstruct_ee)
      CALL cp_fm_struct_release(mstruct_oe)
      CALL cp_fm_struct_release(mstruct_oo)
      !
      CALL qs_env_release(eht_env)
      DEALLOCATE (eht_env)
      CALL section_vals_release(eht_force_env_section)

      CALL timestop(handle)

   END SUBROUTINE calculate_eht_guess

END MODULE qs_eht_guess
