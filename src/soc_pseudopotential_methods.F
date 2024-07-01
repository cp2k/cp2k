!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE soc_pseudopotential_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_cfm_types,                    ONLY: cp_cfm_get_info,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_neighbor_list_set_p,&
                                              neighbor_list_set_p_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'soc_pseudopotential_methods'

   PUBLIC :: V_SOC_xyz_from_pseudopotential, &
             remove_soc_outside_energy_window_ao, &
             remove_soc_outside_energy_window_mo

CONTAINS

! **************************************************************************************************
!> \brief V^SOC_µν^(α),R = ħ/2 < ϕ_µ cell O | sum_ℓ ΔV_ℓ^SO(r,r') L^(α) | ϕ_ν cell R>, α = x,y,z
!>        see Hartwigsen, Goedecker, Hutter, Eq.(18), (19) (doi.org/10.1103/PhysRevB.58.3641)
!>        Caution: V^SOC_µν^(α) is purely imaginary and Hermitian; V^SOC_µν^(α) is stored as real
!>                 dbcsr matrix mat_V_SOC_xyz without symmetry; V^SOC_µν^(α) is stored without
!>                 the imaginary unit, i.e. mat_V_SOC_xyz is real and antisymmetric
!> \param qs_env ...
!> \param mat_V_SOC_xyz ...
!> \par History
!>    * 09.2023 created
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE V_SOC_xyz_from_pseudopotential(qs_env, mat_V_SOC_xyz)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_SOC_xyz

      CHARACTER(LEN=*), PARAMETER :: routineN = 'V_SOC_xyz_from_pseudopotential'

      INTEGER                                            :: handle, img, nder, nimages, xyz
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: calculate_forces, do_kp, do_symmetric, &
                                                            use_virial
      REAL(KIND=dp)                                      :: eps_ppnl
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_l, mat_l_nosym, mat_pot_dummy, &
                                                            matrix_dummy, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, dft_control, sab_orb, sap_ppnl, particle_set, atomic_kind_set, &
               cell_to_index)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, dft_control=dft_control, &
                      matrix_s_kp=matrix_s, kpoints=kpoints, atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, sab_orb=sab_orb, sap_ppnl=sap_ppnl)

      eps_ppnl = dft_control%qs_control%eps_ppnl
      nimages = dft_control%nimages
      do_kp = (nimages > 1)
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_orb, symmetric=do_symmetric)
      CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)

      NULLIFY (mat_l, mat_pot_dummy)
      CALL dbcsr_allocate_matrix_set(mat_l, 3, nimages)
      DO xyz = 1, 3
         DO img = 1, nimages
            ALLOCATE (mat_l(xyz, img)%matrix)
            CALL dbcsr_create(mat_l(xyz, img)%matrix, template=matrix_s(1, 1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(mat_l(xyz, img)%matrix, sab_orb)
            CALL dbcsr_set(mat_l(xyz, img)%matrix, 0.0_dp)
         END DO
      END DO

      ! get mat_l; the next CPASSERT fails if the atoms do not have any SOC parameters, i.e.
      ! SOC is zero and one should not activate the SOC section
      CPASSERT(ASSOCIATED(sap_ppnl))
      nder = 0
      use_virial = .FALSE.
      calculate_forces = .FALSE.

      NULLIFY (mat_pot_dummy)
      CALL dbcsr_allocate_matrix_set(mat_pot_dummy, 1, nimages)
      DO img = 1, nimages
         ALLOCATE (mat_pot_dummy(1, img)%matrix)
         CALL dbcsr_create(mat_pot_dummy(1, img)%matrix, template=matrix_s(1, 1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(mat_pot_dummy(1, img)%matrix, sab_orb)
         CALL dbcsr_set(mat_pot_dummy(1, img)%matrix, 0.0_dp)
      END DO

      CALL build_core_ppnl(mat_pot_dummy, matrix_dummy, force, virial, &
                           calculate_forces, use_virial, nder, &
                           qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, &
                           eps_ppnl, nimages=nimages, cell_to_index=cell_to_index, &
                           basis_type="ORB", matrix_l=mat_l)

      NULLIFY (mat_l_nosym)
      CALL dbcsr_allocate_matrix_set(mat_l_nosym, 3, nimages)
      DO xyz = 1, 3
         DO img = 1, nimages

            ALLOCATE (mat_l_nosym(xyz, img)%matrix)
            IF (do_kp) THEN
               CALL dbcsr_create(mat_l_nosym(xyz, img)%matrix, template=matrix_s(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               CALL dbcsr_copy(mat_l_nosym(xyz, img)%matrix, mat_l(xyz, img)%matrix)
            ELSE
               CALL dbcsr_create(mat_l_nosym(xyz, img)%matrix, template=matrix_s(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_desymmetrize(mat_l(xyz, img)%matrix, mat_l_nosym(xyz, img)%matrix)
            END IF

         END DO
      END DO

      NULLIFY (mat_V_SOC_xyz)
      CALL dbcsr_allocate_matrix_set(mat_V_SOC_xyz, 3, nimages)
      DO xyz = 1, 3
         DO img = 1, nimages
            ALLOCATE (mat_V_SOC_xyz(xyz, img)%matrix)
            IF (do_kp) THEN
               ! mat_V_SOC_xyz^R with neighbor cell R actually has no symmetry
               ! mat_V_SOC_xyz^R_µν = mat_V_SOC_xyz^R_νµ* (the actual symmetry is
               ! mat_V_SOC_xyz^R_µν = mat_V_SOC_xyz^-R_νµ* )  but rskp_transform
               ! for mat_V_SOC_xyz^R -> mat_V_SOC_xyz(k) requires symmetry...
               CALL dbcsr_create(mat_V_SOC_xyz(xyz, img)%matrix, template=matrix_s(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
            ELSE
               CALL dbcsr_create(mat_V_SOC_xyz(xyz, img)%matrix, template=matrix_s(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)
            END IF
            CALL cp_dbcsr_alloc_block_from_nbl(mat_V_SOC_xyz(xyz, img)%matrix, sab_orb)
            ! factor 0.5 from ħ/2 prefactor
            CALL dbcsr_add(mat_V_SOC_xyz(xyz, img)%matrix, mat_l_nosym(xyz, img)%matrix, &
                           0.0_dp, 0.5_dp)
         END DO
      END DO

      CALL dbcsr_deallocate_matrix_set(mat_pot_dummy)
      CALL dbcsr_deallocate_matrix_set(mat_l_nosym)
      CALL dbcsr_deallocate_matrix_set(mat_l)

      CALL timestop(handle)

   END SUBROUTINE V_SOC_xyz_from_pseudopotential

! **************************************************************************************************
!> \brief Remove SOC outside of energy window (otherwise, numerical problems arise
!>        because energetically low semicore states and energetically very high
!>        unbound states couple to the states around the Fermi level).
!>        This routine is for mat_V_SOC_xyz being in the atomic-orbital (ao) basis.
!> \param mat_V_SOC_xyz ...
!> \param e_win ...
!> \param fm_mo_coeff ...
!> \param homo ...
!> \param eigenval ...
!> \param fm_s ...
! **************************************************************************************************
   SUBROUTINE remove_soc_outside_energy_window_ao(mat_V_SOC_xyz, e_win, fm_mo_coeff, homo, &
                                                  eigenval, fm_s)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: mat_V_SOC_xyz
      REAL(KIND=dp)                                      :: e_win
      TYPE(cp_fm_type)                                   :: fm_mo_coeff
      INTEGER                                            :: homo
      REAL(KIND=dp), DIMENSION(:)                        :: eigenval
      TYPE(cp_fm_type)                                   :: fm_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_soc_outside_energy_window_ao'

      INTEGER                                            :: handle, i_glob, iiB, j_glob, jjB, nao, &
                                                            ncol_local, nrow_local, xyz
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: E_HOMO, E_i, E_j, E_LUMO
      TYPE(cp_fm_type)                                   :: fm_V_ao, fm_V_mo, fm_work

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_work, fm_s%matrix_struct)
      CALL cp_fm_create(fm_V_ao, fm_s%matrix_struct)
      CALL cp_fm_create(fm_V_mo, fm_s%matrix_struct)

      CALL cp_fm_get_info(matrix=fm_s, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      nao = SIZE(eigenval)

      E_HOMO = eigenval(homo)
      E_LUMO = eigenval(homo + 1)

      DO xyz = 1, 3

         CALL copy_dbcsr_to_fm(mat_V_SOC_xyz(xyz)%matrix, fm_V_ao)

         ! V_MO = C^T*V_AO*C
         CALL parallel_gemm(transa="N", transb="N", m=nao, n=nao, k=nao, alpha=1.0_dp, &
                            matrix_a=fm_V_ao, matrix_b=fm_mo_coeff, beta=0.0_dp, matrix_c=fm_work)

         CALL parallel_gemm(transa="T", transb="N", m=nao, n=nao, k=nao, alpha=1.0_dp, &
                            matrix_a=fm_mo_coeff, matrix_b=fm_work, beta=0.0_dp, matrix_c=fm_V_mo)

         DO jjB = 1, ncol_local
            j_glob = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_glob = row_indices(iiB)

               E_i = eigenval(i_glob)
               E_j = eigenval(j_glob)

               IF (E_i < E_HOMO - 0.5_dp*e_win .OR. E_i > E_LUMO + 0.5_dp*e_win .OR. &
                   E_j < E_HOMO - 0.5_dp*e_win .OR. E_j > E_LUMO + 0.5_dp*e_win) THEN
                  fm_V_mo%local_data(iiB, jjB) = 0.0_dp
               END IF

            END DO
         END DO

         ! V_AO = S*C*V_MO*C^T*S
         CALL parallel_gemm(transa="N", transb="T", m=nao, n=nao, k=nao, alpha=1.0_dp, &
                            matrix_a=fm_V_mo, matrix_b=fm_mo_coeff, beta=0.0_dp, matrix_c=fm_work)

         CALL parallel_gemm(transa="N", transb="N", m=nao, n=nao, k=nao, alpha=1.0_dp, &
                            matrix_a=fm_mo_coeff, matrix_b=fm_work, beta=0.0_dp, matrix_c=fm_V_ao)

         CALL parallel_gemm(transa="N", transb="N", m=nao, n=nao, k=nao, alpha=1.0_dp, &
                            matrix_a=fm_s, matrix_b=fm_V_ao, beta=0.0_dp, matrix_c=fm_work)

         CALL parallel_gemm(transa="N", transb="N", m=nao, n=nao, k=nao, alpha=1.0_dp, &
                            matrix_a=fm_work, matrix_b=fm_s, beta=0.0_dp, matrix_c=fm_V_ao)

         CALL copy_fm_to_dbcsr(fm_V_ao, mat_V_SOC_xyz(xyz)%matrix)

      END DO

      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_V_ao)
      CALL cp_fm_release(fm_V_mo)

      CALL timestop(handle)

   END SUBROUTINE remove_soc_outside_energy_window_ao

! **************************************************************************************************
!> \brief ...
!> \param cfm_ks_spinor ...
!> \param e_win ...
!> \param eigenval ...
!> \param E_HOMO ...
!> \param E_LUMO ...
! **************************************************************************************************
   SUBROUTINE remove_soc_outside_energy_window_mo(cfm_ks_spinor, e_win, eigenval, E_HOMO, E_LUMO)
      TYPE(cp_cfm_type)                                  :: cfm_ks_spinor
      REAL(KIND=dp)                                      :: e_win
      REAL(KIND=dp), DIMENSION(:)                        :: eigenval
      REAL(KIND=dp)                                      :: E_HOMO, E_LUMO

      CHARACTER(LEN=*), PARAMETER :: routineN = 'remove_soc_outside_energy_window_mo'

      INTEGER                                            :: handle, i_glob, iiB, j_glob, jjB, &
                                                            ncol_global, ncol_local, nrow_global, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: E_i, E_j

      ! Remove SOC outside of energy window (otherwise, numerical problems arise
      ! because energetically low semicore states and energetically very high
      ! unbound states couple to the states around the Fermi level).
      ! This routine is for cfm_ks_spinor being in the molecular-orbital (mo) with
      ! corresponding eigenvalues "eigenval".

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(matrix=cfm_ks_spinor, &
                           nrow_global=nrow_global, &
                           ncol_global=ncol_global, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      CPASSERT(nrow_global == SIZE(eigenval))
      CPASSERT(ncol_global == SIZE(eigenval))

      DO jjB = 1, ncol_local
         j_glob = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_glob = row_indices(iiB)

            E_i = eigenval(i_glob)
            E_j = eigenval(j_glob)

            IF (E_i < E_HOMO - 0.5_dp*e_win .OR. E_i > E_LUMO + 0.5_dp*e_win .OR. &
                E_j < E_HOMO - 0.5_dp*e_win .OR. E_j > E_LUMO + 0.5_dp*e_win) THEN
               cfm_ks_spinor%local_data(iiB, jjB) = 0.0_dp
            END IF

         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE remove_soc_outside_energy_window_mo

END MODULE soc_pseudopotential_methods
