!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief define create destroy get and put information
!>      in xas_env to calculate the x-ray absorption spectra
!> \par History
!>      created 05.2005
!> \author MI (05.2005)
! **************************************************************************************************
MODULE xas_env_types

   USE basis_set_types,                 ONLY: deallocate_gto_basis_set,&
                                              gto_basis_set_p_type
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pools_give_back_fm_vect
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE qs_loc_types,                    ONLY: qs_loc_env_release,&
                                              qs_loc_env_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                              scf_env_release
   USE scf_control_types,               ONLY: scf_c_release,&
                                              scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_env_types'

! *** Public data types ***

   PUBLIC :: xas_environment_type

! *** Public subroutines ***

   PUBLIC :: get_xas_env, set_xas_env, xas_env_create, xas_env_release

! **************************************************************************************************
!> \param nao number of atomic orbitals in the basis
!> \param exc_state state that is now excited (this change atom by atom)
!> \param nvirtual number of empy states to take into account for the spectrum
!> \param state_of_atom for each atom the states that have to be excited (global index)
!>        dimension is the number of atoms to be excited by the largest number of included states
!> \param atom_of_state atom to which each state is assigned,
!>        dimension is the number of states occupied that might be excited
!> \param nexc_states number of states to be excited per atom
!>        dimension is the number of atoms to be excited
!> \param type_of_state character of the state (1s,2s,2p...)
!> \param spectrum for each excitation the energy and the oscillator strength
!> \param centers_wfn for each wfn the center of charge (optimized by localization)
!> \param groundstate_coeff temporary storage for the original mos coefficients
!> \param ostrength_sm sin and cos integrals computed for the contracted GTO functions
!> \param dip_fm_set fm for the sin and cos integrals to define the pos operator
!> \param qs_loc_env environment for the localization procedure
!> \par History
!>       created 05-2005
!> \author MI
! **************************************************************************************************
   TYPE xas_environment_type
      INTEGER :: nao = 0, exc_state = 0, xas_estate = 0
      INTEGER :: nexc_search = 0, nexc_atoms = 0
      INTEGER :: spin_channel = 0
      INTEGER :: nvirtual = 0, nvirtual2 = 0
      INTEGER :: unoccupied_max_iter = 0

      INTEGER, DIMENSION(:), POINTER :: atom_of_state => NULL()
      INTEGER, DIMENSION(:), POINTER :: type_of_state => NULL()
      INTEGER, DIMENSION(:), POINTER :: mykind_of_atom => NULL()
      INTEGER, DIMENSION(:), POINTER :: mykind_of_kind => NULL()
      INTEGER, DIMENSION(:), POINTER :: exc_atoms => NULL()
      INTEGER, DIMENSION(:), POINTER :: nexc_states => NULL()
      INTEGER, DIMENSION(:, :), POINTER :: state_of_atom => NULL()

      REAL(dp) :: ip_energy = 0.0_dp, occ_estate = 0.0_dp, unoccupied_eps = 0.0_dp, xas_nelectron = 0.0_dp, homo_occ = 0.0_dp
      REAL(dp), DIMENSION(:), POINTER :: all_evals => NULL()
      REAL(dp), DIMENSION(:), POINTER :: unoccupied_evals => NULL()
      REAL(dp), DIMENSION(:, :), POINTER :: spectrum => NULL()
      REAL(dp), DIMENSION(:, :), POINTER :: centers_wfn => NULL()
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER :: stogto_overlap => NULL()
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER :: my_gto_basis => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER :: groundstate_coeff => NULL()
      TYPE(cp_fm_type), DIMENSION(:, :), POINTER :: dip_fm_set => NULL()
      TYPE(cp_fm_pool_p_type), DIMENSION(:), &
         POINTER                                   :: ao_mo_fm_pools => NULL()
      TYPE(cp_fm_type), POINTER :: excvec_coeff => NULL()
      TYPE(cp_fm_type), POINTER :: excvec_overlap => NULL()
      TYPE(cp_fm_type), POINTER :: unoccupied_orbs => NULL()
      TYPE(cp_fm_type), POINTER :: all_vectors => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: ostrength_sm => NULL()
      TYPE(qs_loc_env_type), POINTER :: qs_loc_env => NULL()
      TYPE(qs_scf_env_type), POINTER                        :: scf_env => NULL()
      TYPE(scf_control_type), POINTER          :: scf_control => NULL()

   END TYPE xas_environment_type

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param xas_env ...
!> \param exc_state ...
!> \param nao ...
!> \param nvirtual ...
!> \param nvirtual2 ...
!> \param centers_wfn ...
!> \param atom_of_state ...
!> \param exc_atoms ...
!> \param nexc_states ...
!> \param type_of_state ...
!> \param mykind_of_atom ...
!> \param mykind_of_kind ...
!> \param state_of_atom ...
!> \param spectrum ...
!> \param groundstate_coeff ...
!> \param ostrength_sm ...
!> \param dip_fm_set ...
!> \param excvec_coeff ...
!> \param excvec_overlap ...
!> \param unoccupied_orbs ...
!> \param unoccupied_evals ...
!> \param unoccupied_max_iter ...
!> \param unoccupied_eps ...
!> \param all_vectors ...
!> \param all_evals ...
!> \param my_gto_basis ...
!> \param qs_loc_env ...
!> \param stogto_overlap ...
!> \param occ_estate ...
!> \param xas_nelectron ...
!> \param xas_estate ...
!> \param nexc_atoms ...
!> \param nexc_search ...
!> \param spin_channel ...
!> \param scf_env ...
!> \param scf_control ...
! **************************************************************************************************
   SUBROUTINE get_xas_env(xas_env, exc_state, nao, nvirtual, nvirtual2, &
                          centers_wfn, atom_of_state, exc_atoms, nexc_states, type_of_state, mykind_of_atom, &
                          mykind_of_kind, state_of_atom, spectrum, groundstate_coeff, ostrength_sm, &
                          dip_fm_set, excvec_coeff, excvec_overlap, &
                          unoccupied_orbs, unoccupied_evals, unoccupied_max_iter, unoccupied_eps, &
                          all_vectors, all_evals, my_gto_basis, qs_loc_env, &
                          stogto_overlap, occ_estate, xas_nelectron, xas_estate, nexc_atoms, nexc_search, spin_channel, &
                          scf_env, scf_control)

      TYPE(xas_environment_type), INTENT(IN)             :: xas_env
      INTEGER, INTENT(OUT), OPTIONAL                     :: exc_state, nao, nvirtual, nvirtual2
      REAL(dp), DIMENSION(:, :), OPTIONAL, POINTER       :: centers_wfn
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: atom_of_state, exc_atoms, nexc_states, &
                                                            type_of_state, mykind_of_atom, &
                                                            mykind_of_kind
      INTEGER, DIMENSION(:, :), OPTIONAL, POINTER        :: state_of_atom
      REAL(dp), DIMENSION(:, :), OPTIONAL, POINTER       :: spectrum
      TYPE(cp_fm_type), DIMENSION(:), OPTIONAL, POINTER  :: groundstate_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: ostrength_sm
      TYPE(cp_fm_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: dip_fm_set
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: excvec_coeff, excvec_overlap, &
                                                            unoccupied_orbs
      REAL(dp), DIMENSION(:), OPTIONAL, POINTER          :: unoccupied_evals
      INTEGER, INTENT(OUT), OPTIONAL                     :: unoccupied_max_iter
      REAL(dp), OPTIONAL                                 :: unoccupied_eps
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: all_vectors
      REAL(dp), DIMENSION(:), OPTIONAL, POINTER          :: all_evals
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: my_gto_basis
      TYPE(qs_loc_env_type), OPTIONAL, POINTER           :: qs_loc_env
      TYPE(cp_2d_r_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: stogto_overlap
      REAL(dp), INTENT(OUT), OPTIONAL                    :: occ_estate, xas_nelectron
      INTEGER, INTENT(OUT), OPTIONAL                     :: xas_estate, nexc_atoms, nexc_search, &
                                                            spin_channel
      TYPE(qs_scf_env_type), OPTIONAL, POINTER           :: scf_env
      TYPE(scf_control_type), OPTIONAL, POINTER          :: scf_control

      IF (PRESENT(exc_state)) exc_state = xas_env%exc_state
      IF (PRESENT(nao)) nao = xas_env%nao
      IF (PRESENT(nvirtual)) nvirtual = xas_env%nvirtual
      IF (PRESENT(nvirtual2)) nvirtual2 = xas_env%nvirtual2
      IF (PRESENT(xas_nelectron)) xas_nelectron = xas_env%xas_nelectron
      IF (PRESENT(occ_estate)) occ_estate = xas_env%occ_estate
      IF (PRESENT(xas_estate)) xas_estate = xas_env%xas_estate
      IF (PRESENT(nexc_search)) nexc_search = xas_env%nexc_search
      IF (PRESENT(nexc_states)) nexc_states => xas_env%nexc_states
      IF (PRESENT(spin_channel)) spin_channel = xas_env%spin_channel
      IF (PRESENT(nexc_atoms)) nexc_atoms = xas_env%nexc_atoms
      IF (PRESENT(unoccupied_eps)) unoccupied_eps = xas_env%unoccupied_eps
      IF (PRESENT(unoccupied_max_iter)) unoccupied_max_iter = xas_env%unoccupied_max_iter
      IF (PRESENT(centers_wfn)) centers_wfn => xas_env%centers_wfn
      IF (PRESENT(atom_of_state)) atom_of_state => xas_env%atom_of_state
      IF (PRESENT(exc_atoms)) exc_atoms => xas_env%exc_atoms
      IF (PRESENT(type_of_state)) type_of_state => xas_env%type_of_state
      IF (PRESENT(state_of_atom)) state_of_atom => xas_env%state_of_atom
      IF (PRESENT(mykind_of_atom)) mykind_of_atom => xas_env%mykind_of_atom
      IF (PRESENT(mykind_of_kind)) mykind_of_kind => xas_env%mykind_of_kind
      IF (PRESENT(unoccupied_evals)) unoccupied_evals => xas_env%unoccupied_evals
      IF (PRESENT(all_evals)) all_evals => xas_env%all_evals
      IF (PRESENT(spectrum)) spectrum => xas_env%spectrum
      IF (PRESENT(groundstate_coeff)) groundstate_coeff => xas_env%groundstate_coeff
      IF (PRESENT(ostrength_sm)) ostrength_sm => xas_env%ostrength_sm
      IF (PRESENT(excvec_overlap)) excvec_overlap => xas_env%excvec_overlap
      IF (PRESENT(unoccupied_orbs)) unoccupied_orbs => xas_env%unoccupied_orbs
      IF (PRESENT(all_vectors)) all_vectors => xas_env%all_vectors
      IF (PRESENT(dip_fm_set)) dip_fm_set => xas_env%dip_fm_set
      IF (PRESENT(qs_loc_env)) qs_loc_env => xas_env%qs_loc_env
      IF (PRESENT(excvec_coeff)) excvec_coeff => xas_env%excvec_coeff
      IF (PRESENT(my_gto_basis)) my_gto_basis => xas_env%my_gto_basis
      IF (PRESENT(stogto_overlap)) stogto_overlap => xas_env%stogto_overlap
      IF (PRESENT(scf_env)) scf_env => xas_env%scf_env
      IF (PRESENT(scf_control)) scf_control => xas_env%scf_control
   END SUBROUTINE get_xas_env

! **************************************************************************************************
!> \brief ...
!> \param xas_env ...
!> \param nexc_search ...
!> \param spin_channel ...
!> \param nexc_atoms ...
!> \param nvirtual ...
!> \param nvirtual2 ...
!> \param ip_energy ...
!> \param occ_estate ...
!> \param qs_loc_env ...
!> \param xas_estate ...
!> \param xas_nelectron ...
!> \param homo_occ ...
!> \param scf_env ...
!> \param scf_control ...
! **************************************************************************************************
   SUBROUTINE set_xas_env(xas_env, nexc_search, spin_channel, nexc_atoms, &
                          nvirtual, nvirtual2, ip_energy, occ_estate, qs_loc_env, &
                          xas_estate, xas_nelectron, homo_occ, scf_env, scf_control)

      TYPE(xas_environment_type), INTENT(INOUT)          :: xas_env
      INTEGER, INTENT(IN), OPTIONAL                      :: nexc_search, spin_channel, nexc_atoms, &
                                                            nvirtual, nvirtual2
      REAL(dp), INTENT(IN), OPTIONAL                     :: ip_energy, occ_estate
      TYPE(qs_loc_env_type), OPTIONAL, POINTER           :: qs_loc_env
      INTEGER, INTENT(IN), OPTIONAL                      :: xas_estate
      REAL(dp), INTENT(IN), OPTIONAL                     :: xas_nelectron, homo_occ
      TYPE(qs_scf_env_type), OPTIONAL, POINTER           :: scf_env
      TYPE(scf_control_type), OPTIONAL, POINTER          :: scf_control

      IF (PRESENT(nexc_search)) xas_env%nexc_search = nexc_search
      IF (PRESENT(spin_channel)) xas_env%spin_channel = spin_channel
      IF (PRESENT(nexc_atoms)) xas_env%nexc_atoms = nexc_atoms
      IF (PRESENT(nvirtual)) xas_env%nvirtual = nvirtual
      IF (PRESENT(nvirtual2)) xas_env%nvirtual2 = nvirtual2
      IF (PRESENT(occ_estate)) xas_env%occ_estate = occ_estate
      IF (PRESENT(xas_nelectron)) xas_env%xas_nelectron = xas_nelectron
      IF (PRESENT(homo_occ)) xas_env%homo_occ = homo_occ
      IF (PRESENT(xas_estate)) xas_env%xas_estate = xas_estate
      IF (PRESENT(ip_energy)) xas_env%ip_energy = ip_energy
      IF (PRESENT(qs_loc_env)) THEN
         IF (ASSOCIATED(xas_env%qs_loc_env)) THEN
            IF (.NOT. ASSOCIATED(xas_env%qs_loc_env, qs_loc_env)) THEN
               CALL qs_loc_env_release(xas_env%qs_loc_env)
               DEALLOCATE (xas_env%qs_loc_env)
            END IF
         END IF
         xas_env%qs_loc_env => qs_loc_env
      END IF
      IF (PRESENT(scf_env)) THEN ! accept also null pointers ?
         IF (ASSOCIATED(xas_env%scf_env)) THEN
         IF (ASSOCIATED(xas_env%scf_env, scf_env)) THEN
            CALL scf_env_release(xas_env%scf_env)
            DEALLOCATE (xas_env%scf_env)
         END IF
         END IF
         xas_env%scf_env => scf_env
      END IF
      IF (PRESENT(scf_control)) THEN ! accept also null pointers?
         IF (ASSOCIATED(xas_env%scf_control)) THEN
            IF (.NOT. ASSOCIATED(xas_env%scf_control, scf_control)) THEN
               CALL scf_c_release(xas_env%scf_control)
               DEALLOCATE (xas_env%scf_control)
            END IF
         END IF
         xas_env%scf_control => scf_control
      END IF

   END SUBROUTINE set_xas_env

! **************************************************************************************************
!> \brief ...
!> \param xas_env ...
! **************************************************************************************************
   SUBROUTINE xas_env_create(xas_env)

      TYPE(xas_environment_type), INTENT(OUT)            :: xas_env

      xas_env%nvirtual = 0
      xas_env%nvirtual2 = 0

      NULLIFY (xas_env%ao_mo_fm_pools)
      NULLIFY (xas_env%my_gto_basis)
      NULLIFY (xas_env%atom_of_state)
      NULLIFY (xas_env%nexc_states)
      NULLIFY (xas_env%state_of_atom)
      NULLIFY (xas_env%exc_atoms)
      NULLIFY (xas_env%excvec_coeff, xas_env%excvec_overlap)
      NULLIFY (xas_env%type_of_state, xas_env%mykind_of_atom)
      NULLIFY (xas_env%type_of_state, xas_env%mykind_of_kind)
      NULLIFY (xas_env%groundstate_coeff, xas_env%dip_fm_set)
      NULLIFY (xas_env%ostrength_sm, xas_env%qs_loc_env, xas_env%spectrum)
      NULLIFY (xas_env%all_evals, xas_env%all_vectors)
      NULLIFY (xas_env%unoccupied_evals, xas_env%unoccupied_orbs)
      NULLIFY (xas_env%stogto_overlap)
      NULLIFY (xas_env%scf_env)
      NULLIFY (xas_env%scf_control)

   END SUBROUTINE xas_env_create

! **************************************************************************************************
!> \brief ...
!> \param xas_env ...
! **************************************************************************************************
   SUBROUTINE xas_env_release(xas_env)

      TYPE(xas_environment_type), INTENT(INOUT)          :: xas_env

      INTEGER                                            :: ik

      DEALLOCATE (xas_env%state_of_atom, xas_env%atom_of_state)
      DEALLOCATE (xas_env%nexc_states)
      DEALLOCATE (xas_env%type_of_state)
      DEALLOCATE (xas_env%mykind_of_atom)
      DEALLOCATE (xas_env%mykind_of_kind)
      DEALLOCATE (xas_env%exc_atoms)
      DEALLOCATE (xas_env%centers_wfn)
      IF (ASSOCIATED(xas_env%all_evals)) THEN
         DEALLOCATE (xas_env%all_evals)
      END IF
      IF (ASSOCIATED(xas_env%unoccupied_evals)) THEN
         DEALLOCATE (xas_env%unoccupied_evals)
      END IF
      CALL fm_pools_give_back_fm_vect(xas_env%ao_mo_fm_pools, &
                                      xas_env%groundstate_coeff)

      CALL cp_fm_release(xas_env%dip_fm_set)

      IF (ASSOCIATED(xas_env%excvec_coeff)) THEN
         CALL cp_fm_release(xas_env%excvec_coeff)
         DEALLOCATE (xas_env%excvec_coeff)
         NULLIFY (xas_env%excvec_coeff)
      END IF
      IF (ASSOCIATED(xas_env%excvec_overlap)) THEN
         CALL cp_fm_release(xas_env%excvec_overlap)
         DEALLOCATE (xas_env%excvec_overlap)
         NULLIFY (xas_env%excvec_overlap)
      END IF
      IF (ASSOCIATED(xas_env%unoccupied_orbs)) THEN
         CALL cp_fm_release(xas_env%unoccupied_orbs)
         DEALLOCATE (xas_env%unoccupied_orbs)
         NULLIFY (xas_env%unoccupied_orbs)
      END IF
      NULLIFY (xas_env%ao_mo_fm_pools)
      IF (ASSOCIATED(xas_env%all_vectors) .AND. xas_env%nvirtual .GT. 0) THEN
         CALL cp_fm_release(xas_env%all_vectors)
         DEALLOCATE (xas_env%all_vectors)
         NULLIFY (xas_env%all_vectors)
      ELSE
         NULLIFY (xas_env%all_vectors)
      END IF

      IF (ASSOCIATED(xas_env%ostrength_sm)) THEN
         CALL dbcsr_deallocate_matrix_set(xas_env%ostrength_sm)
      END IF
      IF (ASSOCIATED(xas_env%qs_loc_env)) THEN
         CALL qs_loc_env_release(xas_env%qs_loc_env)
         DEALLOCATE (xas_env%qs_loc_env)
      END IF

      IF (ASSOCIATED(xas_env%my_gto_basis)) THEN
         DO ik = 1, SIZE(xas_env%my_gto_basis, 1)
            CALL deallocate_gto_basis_set(xas_env%my_gto_basis(ik)%gto_basis_set)
         END DO
         DEALLOCATE (xas_env%my_gto_basis)
      END IF

      IF (ASSOCIATED(xas_env%stogto_overlap)) THEN
         DO ik = 1, SIZE(xas_env%stogto_overlap, 1)
            DEALLOCATE (xas_env%stogto_overlap(ik)%array)
         END DO
         DEALLOCATE (xas_env%stogto_overlap)
      END IF

      IF (ASSOCIATED(xas_env%scf_env)) THEN
         CALL scf_env_release(xas_env%scf_env)
         DEALLOCATE (xas_env%scf_env)
      END IF
      IF (ASSOCIATED(xas_env%scf_control)) THEN
         CALL scf_c_release(xas_env%scf_control)
         DEALLOCATE (xas_env%scf_control)
      END IF

   END SUBROUTINE xas_env_release

END MODULE xas_env_types

