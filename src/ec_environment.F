!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Energy correction environment setup and handling
!> \par History
!>       2019.09 created
!> \author JGH
! **************************************************************************************************
MODULE ec_environment
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_container_types,       ONLY: add_basis_set_to_container,&
                                              remove_basis_from_container
   USE basis_set_types,                 ONLY: copy_gto_basis_set,&
                                              create_primitive_basis_set,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Niklasson2003,&
                                              Niklasson2014,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_env,                     ONLY: cp_para_env_retain
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE ec_env_types,                    ONLY: energy_correction_type
   USE input_constants,                 ONLY: &
        ec_diagonalization, ec_functional_harris, ec_matrix_sign, ec_matrix_tc2, ec_matrix_trs4, &
        ec_ot_atomic, ec_ot_diag, ec_ot_gs, kg_cholesky, ls_cluster_atomic, ls_cluster_molecular, &
        ls_s_inversion_hotelling, ls_s_inversion_none, ls_s_inversion_sign_sqrt, &
        ls_s_preconditioner_atomic, ls_s_preconditioner_molecular, ls_s_preconditioner_none, &
        ls_s_sqrt_ns, ls_s_sqrt_proot, xc_vdw_fun_nonloc, xc_vdw_fun_pairpot
   USE input_cp2k_check,                ONLY: xc_functionals_expand
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE molecule_types,                  ONLY: molecule_of_atom,&
                                              molecule_type
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE particle_types,                  ONLY: particle_type
   USE qs_dispersion_nonloc,            ONLY: qs_dispersion_nonloc_init
   USE qs_dispersion_pairpot,           ONLY: qs_dispersion_pairpot_init
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_dispersion_utils,             ONLY: qs_dispersion_env_set
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_environment'

   PUBLIC :: ec_env_create
   PUBLIC :: ec_write_input

CONTAINS

! **************************************************************************************************
!> \brief Allocates and intitializes ec_env
!> \param qs_env ...
!> \param ec_env the object to create
!> \param dft_section ...
!> \param input ...
!> \par History
!>       2019.09 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ec_env_create(qs_env, ec_env, dft_section, input)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(section_vals_type), POINTER                   :: dft_section
      TYPE(section_vals_type), OPTIONAL, POINTER         :: input

      CPASSERT(.NOT. ASSOCIATED(ec_env))
      ALLOCATE (ec_env)
      CALL init_ec_env(qs_env, ec_env, dft_section, input)

   END SUBROUTINE ec_env_create

! **************************************************************************************************
!> \brief Initializes ec_env
!> \param qs_env ...
!> \param ec_env ...
!> \param dft_section ...
!> \param input ...
!> \par History
!>       2019.09 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE init_ec_env(qs_env, ec_env, dft_section, input)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(section_vals_type), OPTIONAL, POINTER         :: dft_section, input

      INTEGER                                            :: ikind, maxlgto, nkind, unit_nr
      LOGICAL                                            :: explicit
      REAL(KIND=dp)                                      :: eps_pgf_orb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set, harris_basis
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(section_vals_type), POINTER                   :: ec_section, nl_section, pp_section, &
                                                            section1, section2, xc_section

      NULLIFY (atomic_kind_set, dispersion_env, ec_env%ls_env, para_env)
      NULLIFY (ec_env%sab_orb, ec_env%sac_ppl, ec_env%sap_ppnl)
      NULLIFY (ec_env%matrix_ks, ec_env%matrix_h, ec_env%matrix_s)
      NULLIFY (ec_env%matrix_t, ec_env%matrix_p, ec_env%matrix_w)
      NULLIFY (ec_env%task_list)
      NULLIFY (ec_env%mao_coef)
      NULLIFY (ec_env%dispersion_env)
      NULLIFY (ec_env%xc_section)
      NULLIFY (ec_env%matrix_z)
      NULLIFY (ec_env%matrix_hz)
      NULLIFY (ec_env%matrix_wz)
      NULLIFY (ec_env%z_admm)
      NULLIFY (ec_env%p_env)
      NULLIFY (ec_env%vh_rspace)
      NULLIFY (ec_env%vxc_rspace)
      NULLIFY (ec_env%vtau_rspace)
      NULLIFY (ec_env%vadmm_rspace)
      ec_env%should_update = .TRUE.
      ec_env%mao = .FALSE.

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (qs_env%energy_correction) THEN

         ec_section => section_vals_get_subs_vals(input, "DFT%ENERGY_CORRECTION")
         CALL section_vals_val_get(ec_section, "ALGORITHM", &
                                   i_val=ec_env%ks_solver)
         CALL section_vals_val_get(ec_section, "ENERGY_FUNCTIONAL", &
                                   i_val=ec_env%energy_functional)
         CALL section_vals_val_get(ec_section, "FACTORIZATION", &
                                   i_val=ec_env%factorization)
         CALL section_vals_val_get(ec_section, "OT_INITIAL_GUESS", &
                                   i_val=ec_env%ec_initial_guess)
         CALL section_vals_val_get(ec_section, "EPS_DEFAULT", &
                                   r_val=ec_env%eps_default)
         CALL section_vals_val_get(ec_section, "HARRIS_BASIS", &
                                   c_val=ec_env%basis)
         CALL section_vals_val_get(ec_section, "MAO", &
                                   l_val=ec_env%mao)
         CALL section_vals_val_get(ec_section, "MAO_MAX_ITER", &
                                   i_val=ec_env%mao_max_iter)
         CALL section_vals_val_get(ec_section, "MAO_EPS_GRAD", &
                                   r_val=ec_env%mao_eps_grad)
         ! Skip EC calculation if ground-state calculation did not converge
         CALL section_vals_val_get(ec_section, "SKIP_EC", &
                                   l_val=ec_env%skip_ec)

         ec_env%do_skip = .FALSE.

         ! set basis
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, nkind=nkind)
         CALL uppercase(ec_env%basis)
         SELECT CASE (ec_env%basis)
         CASE ("ORBITAL")
            DO ikind = 1, nkind
               qs_kind => qs_kind_set(ikind)
               CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
               IF (ASSOCIATED(basis_set)) THEN
                  NULLIFY (harris_basis)
                  CALL get_qs_kind(qs_kind=qs_kind, basis_set=harris_basis, basis_type="HARRIS")
                  IF (ASSOCIATED(harris_basis)) THEN
                     CALL remove_basis_from_container(qs_kind%basis_sets, basis_type="HARRIS")
                  END IF
                  NULLIFY (harris_basis)
                  CALL copy_gto_basis_set(basis_set, harris_basis)
                  CALL add_basis_set_to_container(qs_kind%basis_sets, harris_basis, "HARRIS")
               END IF
            END DO
         CASE ("PRIMITIVE")
            DO ikind = 1, nkind
               qs_kind => qs_kind_set(ikind)
               CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
               IF (ASSOCIATED(basis_set)) THEN
                  NULLIFY (harris_basis)
                  CALL get_qs_kind(qs_kind=qs_kind, basis_set=harris_basis, basis_type="HARRIS")
                  IF (ASSOCIATED(harris_basis)) THEN
                     CALL remove_basis_from_container(qs_kind%basis_sets, basis_type="HARRIS")
                  END IF
                  NULLIFY (harris_basis)
                  CALL create_primitive_basis_set(basis_set, harris_basis)
                  CALL get_qs_env(qs_env, dft_control=dft_control)
                  eps_pgf_orb = dft_control%qs_control%eps_pgf_orb
                  CALL init_interaction_radii_orb_basis(harris_basis, eps_pgf_orb)
                  harris_basis%kind_radius = basis_set%kind_radius
                  CALL add_basis_set_to_container(qs_kind%basis_sets, harris_basis, "HARRIS")
               END IF
            END DO
         CASE ("HARRIS")
            DO ikind = 1, nkind
               qs_kind => qs_kind_set(ikind)
               NULLIFY (harris_basis)
               CALL get_qs_kind(qs_kind=qs_kind, basis_set=harris_basis, basis_type="HARRIS")
               IF (.NOT. ASSOCIATED(harris_basis)) THEN
                  CPWARN("Harris Basis not defined for all types of atoms.")
               END IF
            END DO
         CASE DEFAULT
            CPABORT("Unknown basis set for energy correction (Harris functional)")
         END SELECT
         !
         CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, basis_type="HARRIS")
         CALL init_orbital_pointers(maxlgto + 1)
         ! set functional
         SELECT CASE (ec_env%energy_functional)
         CASE (ec_functional_harris)
            ec_env%ec_name = "Harris"
         CASE DEFAULT
            CPABORT("unknown energy correction")
         END SELECT
         ! select the XC section
         NULLIFY (xc_section)
         xc_section => section_vals_get_subs_vals(dft_section, "XC")
         section1 => section_vals_get_subs_vals(ec_section, "XC")
         section2 => section_vals_get_subs_vals(ec_section, "XC%XC_FUNCTIONAL")
         CALL section_vals_get(section2, explicit=explicit)
         IF (explicit) THEN
            CALL xc_functionals_expand(section2, section1)
            ec_env%xc_section => section1
         ELSE
            ec_env%xc_section => xc_section
         END IF
         ! dispersion
         ALLOCATE (dispersion_env)
         NULLIFY (xc_section)
         xc_section => ec_env%xc_section
         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, para_env=para_env)
         CALL qs_dispersion_env_set(dispersion_env, xc_section)
         IF (dispersion_env%type == xc_vdw_fun_pairpot) THEN
            NULLIFY (pp_section)
            pp_section => section_vals_get_subs_vals(xc_section, "VDW_POTENTIAL%PAIR_POTENTIAL")
            CALL qs_dispersion_pairpot_init(atomic_kind_set, qs_kind_set, dispersion_env, pp_section, para_env)
         ELSE IF (dispersion_env%type == xc_vdw_fun_nonloc) THEN
            NULLIFY (nl_section)
            nl_section => section_vals_get_subs_vals(xc_section, "VDW_POTENTIAL%NON_LOCAL")
            CALL qs_dispersion_nonloc_init(dispersion_env, para_env)
         END IF
         ec_env%dispersion_env => dispersion_env

         ! Initialize Harris LS solver environment
         ec_env%use_ls_solver = .FALSE.
         ec_env%use_ls_solver = (ec_env%ks_solver .EQ. ec_matrix_sign) &
                                .OR. (ec_env%ks_solver .EQ. ec_matrix_trs4) &
                                .OR. (ec_env%ks_solver .EQ. ec_matrix_tc2)

         IF (ec_env%use_ls_solver) THEN
            CALL ec_ls_create(qs_env, ec_env)
         END IF

         ! Write input
         IF (unit_nr > 0) THEN
            CALL ec_write_input(ec_env, unit_nr)
         END IF

      END IF

   END SUBROUTINE init_ec_env

! **************************************************************************************************
!> \brief Initializes linear scaling environment for LS based solver of
!>        Harris energy functional and parses input section
!> \param qs_env ...
!> \param ec_env ...
!> \par History
!>       2020.10 created [Fabian Belleflamme]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_ls_create(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_ls_create', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: mu
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ls_scf_env_type), POINTER                     :: ls_env
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: ec_section, input

      CALL timeset(routineN, handle)

      ALLOCATE (ec_env%ls_env)
      ls_env => ec_env%ls_env

      NULLIFY (dft_control, input, ls_env%para_env)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      input=input, &
                      molecule_set=molecule_set, &
                      particle_set=particle_set, &
                      para_env=ls_env%para_env, &
                      nelectron_spin=ls_env%nelectron_spin)

      ! copy some basic stuff
      ls_env%nspins = dft_control%nspins
      ls_env%natoms = SIZE(particle_set, 1)
      CALL cp_para_env_retain(ls_env%para_env)

      ! initialize block to group to defined molecules
      ALLOCATE (ls_env%ls_mstruct%atom_to_molecule(ls_env%natoms))
      CALL molecule_of_atom(molecule_set, atom_to_mol=ls_env%ls_mstruct%atom_to_molecule)

      ls_env%do_transport = .FALSE.
      ls_env%do_pao = .FALSE.
      ls_env%ls_mstruct%do_pao = ls_env%do_pao
      ls_env%do_pexsi = .FALSE.
      ls_env%has_unit_metric = .FALSE.

      ec_section => section_vals_get_subs_vals(input, "DFT%ENERGY_CORRECTION")
      CALL section_vals_val_get(ec_section, "EPS_FILTER", r_val=ls_env%eps_filter)
      CALL section_vals_val_get(ec_section, "MU", r_val=mu)
      CALL section_vals_val_get(ec_section, "FIXED_MU", l_val=ls_env%fixed_mu)
      ls_env%mu_spin = mu
      CALL section_vals_val_get(ec_section, "S_PRECONDITIONER", i_val=ls_env%s_preconditioner_type)
      CALL section_vals_val_get(ec_section, "MATRIX_CLUSTER_TYPE", i_val=ls_env%ls_mstruct%cluster_type)
      CALL section_vals_val_get(ec_section, "SINGLE_PRECISION_MATRICES", l_val=ls_env%ls_mstruct%single_precision)
      CALL section_vals_val_get(ec_section, "S_INVERSION", i_val=ls_env%s_inversion_type)
      CALL section_vals_val_get(ec_section, "CHECK_S_INV", l_val=ls_env%check_s_inv)
      CALL section_vals_val_get(ec_section, "REPORT_ALL_SPARSITIES", l_val=ls_env%report_all_sparsities)
      CALL section_vals_val_get(ec_section, "SIGN_METHOD", i_val=ls_env%sign_method)
      CALL section_vals_val_get(ec_section, "SIGN_ORDER", i_val=ls_env%sign_order)
      CALL section_vals_val_get(ec_section, "DYNAMIC_THRESHOLD", l_val=ls_env%dynamic_threshold)
      CALL section_vals_val_get(ec_section, "NON_MONOTONIC", l_val=ls_env%non_monotonic)
      CALL section_vals_val_get(ec_section, "S_SQRT_METHOD", i_val=ls_env%s_sqrt_method)
      CALL section_vals_val_get(ec_section, "S_SQRT_ORDER", i_val=ls_env%s_sqrt_order)
      CALL section_vals_val_get(ec_section, "EPS_LANCZOS", r_val=ls_env%eps_lanczos)
      CALL section_vals_val_get(ec_section, "MAX_ITER_LANCZOS", i_val=ls_env%max_iter_lanczos)

      SELECT CASE (qs_env%ec_env%ks_solver)
      CASE (ec_matrix_sign)
         ! S inverse required for Sign matrix algorithm,
         ! calculated either by Hotelling or multiplying S matrix sqrt inv
         SELECT CASE (ls_env%s_inversion_type)
         CASE (ls_s_inversion_sign_sqrt)
            ls_env%needs_s_inv = .TRUE.
            ls_env%use_s_sqrt = .TRUE.
         CASE (ls_s_inversion_hotelling)
            ls_env%needs_s_inv = .TRUE.
            ls_env%use_s_sqrt = .FALSE.
         CASE (ls_s_inversion_none)
            ls_env%needs_s_inv = .FALSE.
            ls_env%use_s_sqrt = .FALSE.
         CASE DEFAULT
            CPABORT("")
         END SELECT
      CASE (ec_matrix_trs4, ec_matrix_tc2)
         ls_env%needs_s_inv = .FALSE.
         ls_env%use_s_sqrt = .TRUE.
      CASE DEFAULT
         CPABORT("")
      END SELECT

      SELECT CASE (ls_env%s_preconditioner_type)
      CASE (ls_s_preconditioner_none)
         ls_env%has_s_preconditioner = .FALSE.
      CASE DEFAULT
         ls_env%has_s_preconditioner = .TRUE.
      END SELECT

      ! buffer for the history of matrices, not needed here
      ls_env%extrapolation_order = 0
      ls_env%scf_history%nstore = 0
      ls_env%scf_history%istore = 0
      ALLOCATE (ls_env%scf_history%matrix(ls_env%nspins, ls_env%scf_history%nstore))

      NULLIFY (ls_env%mixing_store)

      CALL timestop(handle)

   END SUBROUTINE ec_ls_create

! **************************************************************************************************
!> \brief Initializes linear scaling environment for LS based solver of
!>        Harris energy functional and parses input section
!> \param ec_env ...
!> \param unit_nr ...
!> \par History
!>       2020.10 created [Fabian Belleflamme]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_write_input(ec_env, unit_nr)
      TYPE(energy_correction_type), POINTER              :: ec_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_write_input', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(ls_scf_env_type), POINTER                     :: ls_env

      CALL timeset(routineN, handle)

      !WRITE (unit_nr, '()')
      WRITE (unit_nr, '(T2,A,A,A,A,A)') "!", REPEAT("-", 29), &
         " Energy Correction ", REPEAT("-", 29), "!"

      ! Algorithm
      SELECT CASE (ec_env%ks_solver)
      CASE (ec_diagonalization)
         WRITE (unit_nr, '(T2,A,T61,A20)') "Algorithm: ", "DIAGONALIZATION"
         WRITE (*, *) "print diagonalization"
      CASE (ec_ot_diag)
         WRITE (unit_nr, '(T2,A,T61,A20)') "Algorithm: ", "OT DIAGONALIZATION"
         WRITE (*, *) "print ot diagonalization"
      CASE (ec_matrix_sign)
         WRITE (unit_nr, '(T2,A,T61,A20)') "Algorithm: ", "MATRIX_SIGN"
      CASE (ec_matrix_trs4)
         WRITE (unit_nr, '(T2,A,T61,A20)') "Algorithm: ", "TRS4"
         CALL cite_reference(Niklasson2003)
      CASE (ec_matrix_tc2)
         WRITE (unit_nr, '(T2,A,T61,A20)') "Algorithm: ", "TC2"
         CALL cite_reference(Niklasson2014)
      END SELECT

      ! Harris functional parameters
      WRITE (unit_nr, '(T2,A,T61,E20.3)') "eps_default:", ec_env%eps_default

      CALL uppercase(ec_env%basis)
      SELECT CASE (ec_env%basis)
      CASE ("ORBITAL")
         WRITE (unit_nr, '(T2,A,T61,A20)') "Harris basis: ", "ORBITAL"
      CASE ("PRIMITIVE")
         WRITE (unit_nr, '(T2,A,T61,A20)') "Harris basis: ", "PRIMITIVE"
      CASE ("HARRIS")
         WRITE (unit_nr, '(T2,A,T61,A20)') "Harris Basis: ", "HARRIS"
      END SELECT
      WRITE (unit_nr, '()')

      ! MAO
      IF (ec_env%mao) THEN
         WRITE (unit_nr, '(T2,A,T61,L20)') "MAO:", ec_env%mao
         WRITE (unit_nr, '(T2,A,T61,I20)') "MAO_MAX_ITER:", ec_env%mao_max_iter
         WRITE (unit_nr, '(T2,A,T61,E20.3)') "MAO_EPS_GRAD:", ec_env%mao_eps_grad
         WRITE (unit_nr, '()')
      END IF

      ! Parameters depending on chosen algorithm to solve Harris functional
      IF (.NOT. ec_env%use_ls_solver) THEN

         WRITE (unit_nr, '(T2,A)') "MO Solver"
         WRITE (unit_nr, '()')

         SELECT CASE (ec_env%ks_solver)
         CASE (ec_diagonalization)

            SELECT CASE (ec_env%factorization)
            CASE (kg_cholesky)
               WRITE (unit_nr, '(T2,A,T61,A20)') "Factorization: ", "CHOLESKY"
            END SELECT

         CASE (ec_ot_diag)

            ! OT Diagonalization
            ! Initial guess : 1) block diagonal initial guess
            !                 2) GS-density matrix (might require trafo if basis diff)

            SELECT CASE (ec_env%ec_initial_guess)
            CASE (ec_ot_atomic)
               WRITE (unit_nr, '(T2,A,T61,A20)') "OT Diag initial guess: ", "ATOMIC"
            CASE (ec_ot_gs)
               WRITE (unit_nr, '(T2,A,T61,A20)') "OT Diag initial guess: ", "GROUND STATE DM"
            END SELECT

         CASE DEFAULT
            CPABORT("Unknown Diagonalization algorithm for Harris functional")
         END SELECT

      ELSE

         WRITE (unit_nr, '(T2,A)') "AO Solver"
         WRITE (unit_nr, '()')

         ls_env => ec_env%ls_env
         WRITE (unit_nr, '(T2,A,T61,E20.3)') "eps_filter:", ls_env%eps_filter
         WRITE (unit_nr, '(T2,A,T61,L20)') "fixed chemical potential (mu)", ls_env%fixed_mu
         WRITE (unit_nr, '(T2,A,T61,L20)') "Computing inv(S):", ls_env%needs_s_inv
         WRITE (unit_nr, '(T2,A,T61,L20)') "Computing sqrt(S):", ls_env%use_s_sqrt
         WRITE (unit_nr, '(T2,A,T61,L20)') "Computing S preconditioner ", ls_env%has_s_preconditioner
         WRITE (unit_nr, '(T2,A,T61,L20)') "has unit metric:", ls_env%has_unit_metric
         WRITE (unit_nr, '(T2,A,T61,L20)') "Use single precision matrices", &
            ls_env%ls_mstruct%single_precision

         IF (ls_env%use_s_sqrt) THEN
            SELECT CASE (ls_env%s_sqrt_method)
            CASE (ls_s_sqrt_ns)
               WRITE (unit_nr, '(T2,A,T61,A20)') "S sqrt method:", "NEWTONSCHULZ"
            CASE (ls_s_sqrt_proot)
               WRITE (unit_nr, '(T2,A,T61,A20)') "S sqrt method:", "PROOT"
            CASE DEFAULT
               CPABORT("Unknown sqrt method.")
            END SELECT
            WRITE (unit_nr, '(T2,A,T61,I20)') "S sqrt order:", ls_env%s_sqrt_order

         END IF

         SELECT CASE (ls_env%s_preconditioner_type)
         CASE (ls_s_preconditioner_none)
            WRITE (unit_nr, '(T2,A,T61,A20)') "S preconditioner type ", "NONE"
         CASE (ls_s_preconditioner_atomic)
            WRITE (unit_nr, '(T2,A,T61,A20)') "S preconditioner type ", "ATOMIC"
         CASE (ls_s_preconditioner_molecular)
            WRITE (unit_nr, '(T2,A,T61,A20)') "S preconditioner type ", "MOLECULAR"
         END SELECT

         SELECT CASE (ls_env%ls_mstruct%cluster_type)
         CASE (ls_cluster_atomic)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Cluster type", ADJUSTR("ATOMIC")
         CASE (ls_cluster_molecular)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Cluster type", ADJUSTR("MOLECULAR")
         CASE DEFAULT
            CPABORT("Unknown cluster type")
         END SELECT

      END IF

      WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
      WRITE (unit_nr, '()')

      CALL timestop(handle)

   END SUBROUTINE ec_write_input

END MODULE ec_environment
