!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for an energy correction on top of a Kohn-Sham calculation
!> \par History
!>       03.2026 created from energy_correction
!> \author JGH
! **************************************************************************************************
MODULE ec_diag_solver
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_filter, dbcsr_get_info, dbcsr_p_type, dbcsr_release, dbcsr_scale, &
        dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_init_random,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_unit_nr
   USE dm_ls_scf,                       ONLY: calculate_w_matrix_ls,&
                                              post_scf_sparsities
   USE dm_ls_scf_methods,               ONLY: apply_matrix_preconditioner,&
                                              density_matrix_sign,&
                                              density_matrix_tc2,&
                                              density_matrix_trs4,&
                                              ls_scf_init_matrix_s
   USE dm_ls_scf_qs,                    ONLY: matrix_ls_create,&
                                              matrix_ls_to_qs,&
                                              matrix_qs_to_ls
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE ec_env_types,                    ONLY: energy_correction_type
   USE ec_methods,                      ONLY: ec_mos_init
   USE input_constants,                 ONLY: ec_matrix_sign,&
                                              ec_matrix_tc2,&
                                              ec_matrix_trs4,&
                                              ec_ot_atomic,&
                                              ec_ot_gs,&
                                              ot_precond_full_single_inverse,&
                                              ot_precond_solver_default
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE preconditioner,                  ONLY: make_preconditioner
   USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                              init_preconditioner,&
                                              preconditioner_type
   USE qs_atomic_block,                 ONLY: calculate_atomic_block_dm
   USE qs_density_matrices,             ONLY: calculate_density_matrix,&
                                              calculate_w_matrix
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues,&
                                              make_basis_sm
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_type
   USE qs_ot_eigensolver,               ONLY: ot_eigensolver
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_diag_solver'

   PUBLIC :: ec_diag_solver_gamma, ec_diag_solver_kp, ec_ot_diag_solver, ec_ls_init, ec_ls_solver

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Solve KS equation using diagonalization
!> \param qs_env ...
!> \param ec_env ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param matrix_p ...
!> \param matrix_w ...
!> \par History
!>      03.2014 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ec_diag_solver_gamma(qs_env, ec_env, matrix_ks, matrix_s, matrix_p, matrix_w)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type)                       :: ec_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s, matrix_p, matrix_w

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_diag_solver_gamma'

      INTEGER                                            :: handle, ispin, nmo(2), nsize, nspins
      REAL(KIND=dp)                                      :: eps_filter, flexible_electron_count, &
                                                            focc(2), n_el_f
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigvals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_k, fm_s, fm_w
      TYPE(cp_fm_type), POINTER                          :: fm_mos
      TYPE(dbcsr_type), POINTER                          :: ref_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: moset
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      eps_filter = dft_control%qs_control%eps_filter_matrix
      nspins = dft_control%nspins

      nmo = 0
      CALL get_qs_env(qs_env=qs_env, nelectron_spin=nmo)
      focc = 1._dp
      IF (nspins == 1) focc = 2._dp

      CALL dbcsr_get_info(matrix_ks(1, 1)%matrix, nfullrows_total=nsize)

      NULLIFY (blacs_env, para_env)
      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, para_env=para_env)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nsize, &
                               ncol_global=nsize, para_env=para_env)
      CALL cp_fm_create(fm_k, fm_struct)
      CALL cp_fm_create(fm_s, fm_struct)
      CALL cp_fm_create(fm_w, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ! MO sets
      flexible_electron_count = 0.0_dp
      IF (ec_env%smear%do_smear) flexible_electron_count = 0.0001_dp

      ALLOCATE (moset(nspins))
      DO ispin = 1, nspins
         n_el_f = nmo(ispin)
         CALL allocate_mo_set(moset(ispin), nsize, nsize, nmo(ispin), n_el_f, &
                              focc(ispin), flexible_electron_count)
         CALL init_mo_set(moset(ispin), fm_ref=fm_w, name="MO")
      END DO

      DO ispin = 1, nspins
         ref_matrix => matrix_s(1, 1)%matrix
         CALL copy_dbcsr_to_fm(ref_matrix, fm_s)
         ref_matrix => matrix_ks(ispin, 1)%matrix
         CALL copy_dbcsr_to_fm(ref_matrix, fm_k)
         CALL get_mo_set(moset(ispin), eigenvalues=eigvals, mo_coeff=fm_mos)
         CALL cp_fm_geeig(fm_k, fm_s, fm_mos, eigvals, fm_w)
      END DO
      !
      CALL set_mo_occupation(moset, ec_env%smear)
      !
      DO ispin = 1, nspins
         ec_env%ekTS = ec_env%ekTS + moset(ispin)%kTS
         CALL calculate_density_matrix(moset(ispin), matrix_p(ispin, 1)%matrix)
         CALL calculate_w_matrix(moset(ispin), matrix_w(ispin, 1)%matrix)
      END DO
      !
      CALL cp_fm_release(fm_k)
      CALL cp_fm_release(fm_s)
      CALL cp_fm_release(fm_w)
      DO ispin = 1, nspins
         CALL deallocate_mo_set(moset(ispin))
      END DO
      DEALLOCATE (moset)

      CALL timestop(handle)

   END SUBROUTINE ec_diag_solver_gamma

! **************************************************************************************************
!> \brief Solve Kpoint-KS equation using diagonalization
!> \param qs_env ...
!> \param ec_env ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param matrix_p ...
!> \param matrix_w ...
!> \par History
!>      03.2026 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ec_diag_solver_kp(qs_env, ec_env, matrix_ks, matrix_s, matrix_p, matrix_w)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type)                       :: ec_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s, matrix_p, matrix_w

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_diag_solver_kp'

      INTEGER                                            :: handle, ispin, nmo(2), nsize, nspins
      REAL(KIND=dp)                                      :: eps_filter, flexible_electron_count, &
                                                            focc(2), n_el_f
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigvals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_k, fm_s, fm_w
      TYPE(cp_fm_type), POINTER                          :: fm_mos
      TYPE(dbcsr_type), POINTER                          :: ref_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), ALLOCATABLE, DIMENSION(:)       :: moset
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      eps_filter = dft_control%qs_control%eps_filter_matrix
      nspins = dft_control%nspins

      nmo = 0
      CALL get_qs_env(qs_env=qs_env, nelectron_spin=nmo)
      focc = 1._dp
      IF (nspins == 1) focc = 2._dp

      CALL dbcsr_get_info(matrix_ks(1, 1)%matrix, nfullrows_total=nsize)

      NULLIFY (blacs_env, para_env)
      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, para_env=para_env)
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nsize, &
                               ncol_global=nsize, para_env=para_env)
      CALL cp_fm_create(fm_k, fm_struct)
      CALL cp_fm_create(fm_s, fm_struct)
      CALL cp_fm_create(fm_w, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ! MO sets
      flexible_electron_count = 0.0_dp
      IF (ec_env%smear%do_smear) flexible_electron_count = 0.0001_dp

      ALLOCATE (moset(nspins))
      DO ispin = 1, nspins
         n_el_f = nmo(ispin)
         CALL allocate_mo_set(moset(ispin), nsize, nsize, nmo(ispin), n_el_f, &
                              focc(ispin), flexible_electron_count)
         CALL init_mo_set(moset(ispin), fm_ref=fm_w, name="MO")
      END DO

      DO ispin = 1, nspins
         ref_matrix => matrix_s(1, 1)%matrix
         CALL copy_dbcsr_to_fm(ref_matrix, fm_s)
         ref_matrix => matrix_ks(ispin, 1)%matrix
         CALL copy_dbcsr_to_fm(ref_matrix, fm_k)
         CALL get_mo_set(moset(ispin), eigenvalues=eigvals, mo_coeff=fm_mos)
         CALL cp_fm_geeig(fm_k, fm_s, fm_mos, eigvals, fm_w)
      END DO
!deb  CALL diag_kp_basic(matrix_ks, matrix_s, kpoints, fmwork)
!deb
!     ! MO occupations
!     CALL kpoint_set_mo_occupation(kpoints, smear)
!     ! density matrices
!     CALL kpoint_density_matrices(kpoints)
!     ! density matrices in real space
!     CALL kpoint_density_transform(kpoints, matrix_p, .FALSE., tempmat, sab_nl, fmwork)

!     energy%kTS = 0.0_dp
!     energy%efermi = 0.0_dp
!     mos => kpoints%kp_env(1)%kpoint_env%mos
!     DO ispin = 1, SIZE(mos, 2)
!        energy%kTS = energy%kTS + mos(1, ispin)%kTS
!        energy%efermi = energy%efermi + mos(1, ispin)%mu
!     END DO
!     energy%efermi = energy%efermi/REAL(SIZE(mos, 2), KIND=dp)

!deb
      !
      CALL set_mo_occupation(moset, ec_env%smear)
      !
      DO ispin = 1, nspins
         ec_env%ekTS = ec_env%ekTS + moset(ispin)%kTS
         CALL calculate_density_matrix(moset(ispin), matrix_p(ispin, 1)%matrix)
         CALL calculate_w_matrix(moset(ispin), matrix_w(ispin, 1)%matrix)
      END DO
      !
      CALL cp_fm_release(fm_k)
      CALL cp_fm_release(fm_s)
      CALL cp_fm_release(fm_w)
      DO ispin = 1, nspins
         CALL deallocate_mo_set(moset(ispin))
      END DO
      DEALLOCATE (moset)

      CALL timestop(handle)

   END SUBROUTINE ec_diag_solver_kp

! **************************************************************************************************
!> \brief Use OT-diagonalziation to obtain density matrix from Harris Kohn-Sham matrix
!>        Initial guess of density matrix is either the atomic block initial guess from SCF
!>        or the ground-state density matrix. The latter only works if the same basis is used
!>
!> \param qs_env ...
!> \param ec_env ...
!> \param matrix_ks Harris Kohn-Sham matrix
!> \param matrix_s Overlap matrix in Harris functional basis
!> \param matrix_p Harris dentiy matrix, calculated here
!> \param matrix_w Harris energy weighted density matrix, calculated here
!>
!> \par History
!>       09.2020 created
!> \author F.Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_ot_diag_solver(qs_env, ec_env, matrix_ks, matrix_s, matrix_p, matrix_w)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN), &
         POINTER                                         :: matrix_ks, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), &
         INTENT(INOUT), POINTER                          :: matrix_p, matrix_w

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ec_ot_diag_solver'

      INTEGER                                            :: handle, homo, ikind, iounit, ispin, &
                                                            max_iter, nao, nkind, nmo, nspins
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type)                                   :: sv
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set, harris_basis
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(preconditioner_type), POINTER                 :: local_preconditioner
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      iounit = cp_logger_get_default_unit_nr(local=.FALSE.)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(ec_env))
      CPASSERT(ASSOCIATED(matrix_ks))
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(matrix_p))
      CPASSERT(ASSOCIATED(matrix_w))

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      blacs_env=blacs_env, &
                      dft_control=dft_control, &
                      nelectron_spin=nelectron_spin, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)
      nspins = dft_control%nspins

      ! Maximum number of OT iterations for diagonalization
      max_iter = 200

      ! If linear scaling, need to allocate and init MO set
      ! set it to qs_env%mos
      IF (dft_control%qs_control%do_ls_scf) THEN
         CALL ec_mos_init(qs_env, matrix_s(1, 1)%matrix)
      END IF

      CALL get_qs_env(qs_env, mos=mos)

      ! Inital guess to use
      NULLIFY (p_rmpv)

      ! Using ether ground-state DM or ATOMIC GUESS requires
      ! Harris functional to use the same basis set
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, nkind=nkind)
      CALL uppercase(ec_env%basis)
      ! Harris basis only differs from ground-state basis if explicitly added
      ! thus only two cases that need to be tested
      ! 1) explicit Harris basis present?
      IF (ec_env%basis == "HARRIS") THEN
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            ! Basis sets of ground-state
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
            ! Basis sets of energy correction
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=harris_basis, basis_type="HARRIS")

            IF (basis_set%name /= harris_basis%name) THEN
               CPABORT("OT-Diag initial guess: Harris and ground state need to use the same basis")
            END IF
         END DO
      END IF
      ! 2) Harris uses MAOs
      IF (ec_env%mao) THEN
         CPABORT("OT-Diag initial guess: not implemented for MAOs")
      END IF

      ! Initital guess obtained for OT Diagonalization
      SELECT CASE (ec_env%ec_initial_guess)
      CASE (ec_ot_atomic)

         p_rmpv => matrix_p(:, 1)

         CALL calculate_atomic_block_dm(p_rmpv, matrix_s(1, 1)%matrix, atomic_kind_set, qs_kind_set, &
                                        nspins, nelectron_spin, iounit, para_env)

      CASE (ec_ot_gs)

         CALL get_qs_env(qs_env, rho=rho)
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao)
         p_rmpv => rho_ao(:, 1)

      CASE DEFAULT
         CPABORT("Unknown inital guess for OT-Diagonalization (Harris functional)")
      END SELECT

      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), &
                         mo_coeff=mo_coeff, &
                         nmo=nmo, &
                         nao=nao, &
                         homo=homo)

         ! Calculate first MOs
         CALL cp_fm_set_all(mo_coeff, 0.0_dp)
         CALL cp_fm_init_random(mo_coeff, nmo)

         CALL cp_fm_create(sv, mo_coeff%matrix_struct, "SV")
         ! multiply times PS
         ! PS*C(:,1:nomo)+C(:,nomo+1:nmo) (nomo=NINT(nelectron/maxocc))
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1, 1)%matrix, mo_coeff, sv, nmo)
         CALL cp_dbcsr_sm_fm_multiply(p_rmpv(ispin)%matrix, sv, mo_coeff, homo)
         CALL cp_fm_release(sv)
         ! and ortho the result
         ! If DFBT or SE, then needs has_unit_metrix option
         CALL make_basis_sm(mo_coeff, nmo, matrix_s(1, 1)%matrix)
      END DO

      ! Preconditioner
      NULLIFY (local_preconditioner)
      ALLOCATE (local_preconditioner)
      CALL init_preconditioner(local_preconditioner, para_env=para_env, &
                               blacs_env=blacs_env)
      DO ispin = 1, nspins
         CALL make_preconditioner(local_preconditioner, &
                                  precon_type=ot_precond_full_single_inverse, &
                                  solver_type=ot_precond_solver_default, &
                                  matrix_h=matrix_ks(ispin, 1)%matrix, &
                                  matrix_s=matrix_s(ispin, 1)%matrix, &
                                  convert_precond_to_dbcsr=.TRUE., &
                                  mo_set=mos(ispin), energy_gap=0.2_dp)

         CALL get_mo_set(mos(ispin), &
                         mo_coeff=mo_coeff, &
                         eigenvalues=eigenvalues, &
                         nmo=nmo, &
                         homo=homo)
         CALL ot_eigensolver(matrix_h=matrix_ks(ispin, 1)%matrix, &
                             matrix_s=matrix_s(1, 1)%matrix, &
                             matrix_c_fm=mo_coeff, &
                             preconditioner=local_preconditioner, &
                             eps_gradient=ec_env%eps_default, &
                             iter_max=max_iter, &
                             silent=.FALSE.)
         CALL calculate_subspace_eigenvalues(mo_coeff, matrix_ks(ispin, 1)%matrix, &
                                             evals_arg=eigenvalues, do_rotation=.TRUE.)

         ! Deallocate preconditioner
         CALL destroy_preconditioner(local_preconditioner)
         DEALLOCATE (local_preconditioner)

         !fm->dbcsr
         CALL copy_fm_to_dbcsr(mos(ispin)%mo_coeff, &
                               mos(ispin)%mo_coeff_b)
      END DO

      ! Calculate density matrix from MOs
      DO ispin = 1, nspins
         CALL calculate_density_matrix(mos(ispin), matrix_p(ispin, 1)%matrix)

         CALL calculate_w_matrix(mos(ispin), matrix_w(ispin, 1)%matrix)
      END DO

      ! Get rid of MO environment again
      IF (dft_control%qs_control%do_ls_scf) THEN
         DO ispin = 1, nspins
            CALL deallocate_mo_set(mos(ispin))
         END DO
         IF (ASSOCIATED(qs_env%mos)) THEN
            DO ispin = 1, SIZE(qs_env%mos)
               CALL deallocate_mo_set(qs_env%mos(ispin))
            END DO
            DEALLOCATE (qs_env%mos)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_ot_diag_solver

! **************************************************************************************************
!> \brief Solve the Harris functional by linear scaling density purification scheme,
!>        instead of the diagonalization performed in ec_diag_solver
!>
!> \param qs_env ...
!> \param matrix_ks Harris Kohn-Sham matrix
!> \param matrix_s Overlap matrix in Harris functional basis
!> \par History
!>       09.2020 created
!> \author F.Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_ls_init(qs_env, matrix_ks, matrix_s)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ec_ls_init'

      INTEGER                                            :: handle, ispin, nspins
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(ls_scf_env_type), POINTER                     :: ls_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      ec_env=ec_env)
      nspins = dft_control%nspins
      ls_env => ec_env%ls_env

      ! create the matrix template for use in the ls procedures
      CALL matrix_ls_create(matrix_ls=ls_env%matrix_s, matrix_qs=matrix_s(1, 1)%matrix, &
                            ls_mstruct=ls_env%ls_mstruct)

      IF (ALLOCATED(ls_env%matrix_p)) THEN
         DO ispin = 1, SIZE(ls_env%matrix_p)
            CALL dbcsr_release(ls_env%matrix_p(ispin))
         END DO
      ELSE
         ALLOCATE (ls_env%matrix_p(nspins))
      END IF

      DO ispin = 1, nspins
         CALL dbcsr_create(ls_env%matrix_p(ispin), template=ls_env%matrix_s, &
                           matrix_type=dbcsr_type_no_symmetry)
      END DO

      ALLOCATE (ls_env%matrix_ks(nspins))
      DO ispin = 1, nspins
         CALL dbcsr_create(ls_env%matrix_ks(ispin), template=ls_env%matrix_s, &
                           matrix_type=dbcsr_type_no_symmetry)
      END DO

      ! Set up S matrix and needed functions of S
      CALL ls_scf_init_matrix_s(matrix_s(1, 1)%matrix, ls_env)

      ! Bring KS matrix from QS to LS form
      ! EC KS-matrix already calculated
      DO ispin = 1, nspins
         CALL matrix_qs_to_ls(matrix_ls=ls_env%matrix_ks(ispin), &
                              matrix_qs=matrix_ks(ispin, 1)%matrix, &
                              ls_mstruct=ls_env%ls_mstruct, &
                              covariant=.TRUE.)
      END DO

      CALL timestop(handle)

   END SUBROUTINE ec_ls_init

! **************************************************************************************************
!> \brief Solve the Harris functional by linear scaling density purification scheme,
!>        instead of the diagonalization performed in ec_diag_solver
!>
!> \param qs_env ...
!> \param matrix_p Harris dentiy matrix, calculated here
!> \param matrix_w Harris energy weighted density matrix, calculated here
!> \param ec_ls_method which purification scheme should be used
!> \par History
!>      12.2019 created [JGH]
!>      08.2020 refactoring [fbelle]
!> \author Fabian Belleflamme
! **************************************************************************************************

   SUBROUTINE ec_ls_solver(qs_env, matrix_p, matrix_w, ec_ls_method)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_w
      INTEGER, INTENT(IN)                                :: ec_ls_method

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_ls_solver'

      INTEGER                                            :: handle, ispin, nelectron_spin_real, &
                                                            nspins
      INTEGER, DIMENSION(2)                              :: nmo
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: wmat
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: matrix_ks_deviation
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(ls_scf_env_type), POINTER                     :: ls_env
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (para_env)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env)
      nspins = dft_control%nspins
      ec_env => qs_env%ec_env
      ls_env => ec_env%ls_env

      nmo = 0
      CALL get_qs_env(qs_env=qs_env, nelectron_spin=nmo)
      IF (nspins == 1) nmo(1) = nmo(1)/2
      ls_env%homo_spin(:) = 0.0_dp
      ls_env%lumo_spin(:) = 0.0_dp

      ALLOCATE (matrix_ks_deviation(nspins))
      DO ispin = 1, nspins
         CALL dbcsr_create(matrix_ks_deviation(ispin), template=ls_env%matrix_ks(ispin))
         CALL dbcsr_set(matrix_ks_deviation(ispin), 0.0_dp)
      END DO

      ! F = S^-1/2 * F * S^-1/2
      IF (ls_env%has_s_preconditioner) THEN
         DO ispin = 1, nspins
            CALL apply_matrix_preconditioner(ls_env%matrix_ks(ispin), "forward", &
                                             ls_env%matrix_bs_sqrt, ls_env%matrix_bs_sqrt_inv)

            CALL dbcsr_filter(ls_env%matrix_ks(ispin), ls_env%eps_filter)
         END DO
      END IF

      DO ispin = 1, nspins
         nelectron_spin_real = ls_env%nelectron_spin(ispin)
         IF (ls_env%nspins == 1) nelectron_spin_real = nelectron_spin_real/2

         SELECT CASE (ec_ls_method)
         CASE (ec_matrix_sign)
            CALL density_matrix_sign(ls_env%matrix_p(ispin), &
                                     ls_env%mu_spin(ispin), &
                                     ls_env%fixed_mu, &
                                     ls_env%sign_method, &
                                     ls_env%sign_order, &
                                     ls_env%matrix_ks(ispin), &
                                     ls_env%matrix_s, &
                                     ls_env%matrix_s_inv, &
                                     nelectron_spin_real, &
                                     ec_env%eps_default)

         CASE (ec_matrix_trs4)
            CALL density_matrix_trs4( &
               ls_env%matrix_p(ispin), &
               ls_env%matrix_ks(ispin), &
               ls_env%matrix_s_sqrt_inv, &
               nelectron_spin_real, &
               ec_env%eps_default, &
               ls_env%homo_spin(ispin), &
               ls_env%lumo_spin(ispin), &
               ls_env%mu_spin(ispin), &
               matrix_ks_deviation=matrix_ks_deviation(ispin), &
               dynamic_threshold=ls_env%dynamic_threshold, &
               eps_lanczos=ls_env%eps_lanczos, &
               max_iter_lanczos=ls_env%max_iter_lanczos)

         CASE (ec_matrix_tc2)
            CALL density_matrix_tc2( &
               ls_env%matrix_p(ispin), &
               ls_env%matrix_ks(ispin), &
               ls_env%matrix_s_sqrt_inv, &
               nelectron_spin_real, &
               ec_env%eps_default, &
               ls_env%homo_spin(ispin), &
               ls_env%lumo_spin(ispin), &
               non_monotonic=ls_env%non_monotonic, &
               eps_lanczos=ls_env%eps_lanczos, &
               max_iter_lanczos=ls_env%max_iter_lanczos)

         END SELECT

      END DO

      ! de-orthonormalize
      IF (ls_env%has_s_preconditioner) THEN
         DO ispin = 1, nspins
            ! P = S^-1/2 * P_tilde * S^-1/2 (forward)
            CALL apply_matrix_preconditioner(ls_env%matrix_p(ispin), "forward", &
                                             ls_env%matrix_bs_sqrt, ls_env%matrix_bs_sqrt_inv)

            CALL dbcsr_filter(ls_env%matrix_p(ispin), ls_env%eps_filter)
         END DO
      END IF

      ! Closed-shell
      IF (nspins == 1) CALL dbcsr_scale(ls_env%matrix_p(1), 2.0_dp)

      IF (ls_env%report_all_sparsities) CALL post_scf_sparsities(ls_env)

      ! ls_scf_dm_to_ks
      ! Density matrix from LS to EC
      DO ispin = 1, nspins
         CALL matrix_ls_to_qs(matrix_qs=matrix_p(ispin, 1)%matrix, &
                              matrix_ls=ls_env%matrix_p(ispin), &
                              ls_mstruct=ls_env%ls_mstruct, &
                              covariant=.FALSE.)
      END DO

      wmat => matrix_w(:, 1)
      CALL calculate_w_matrix_ls(wmat, ec_env%ls_env)

      ! clean up
      CALL dbcsr_release(ls_env%matrix_s)
      IF (ls_env%has_s_preconditioner) THEN
         CALL dbcsr_release(ls_env%matrix_bs_sqrt)
         CALL dbcsr_release(ls_env%matrix_bs_sqrt_inv)
      END IF
      IF (ls_env%needs_s_inv) THEN
         CALL dbcsr_release(ls_env%matrix_s_inv)
      END IF
      IF (ls_env%use_s_sqrt) THEN
         CALL dbcsr_release(ls_env%matrix_s_sqrt)
         CALL dbcsr_release(ls_env%matrix_s_sqrt_inv)
      END IF

      DO ispin = 1, SIZE(ls_env%matrix_ks)
         CALL dbcsr_release(ls_env%matrix_ks(ispin))
      END DO
      DEALLOCATE (ls_env%matrix_ks)

      DO ispin = 1, nspins
         CALL dbcsr_release(matrix_ks_deviation(ispin))
      END DO
      DEALLOCATE (matrix_ks_deviation)

      CALL timestop(handle)

   END SUBROUTINE ec_ls_solver

END MODULE ec_diag_solver

