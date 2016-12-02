!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for a linear scaling quickstep SCF run based on the density
!>        matrix
!> \par History
!>       2010.10 created [Joost VandeVondele]
!>       2016.11 created from dm_ls_scf to avoid circular dependencies
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE dm_ls_scf_create
   USE bibliography,                    ONLY: Lin2009,&
                                              Lin2013,&
                                              Niklasson2003,&
                                              Niklasson2014,&
                                              Shao2003,&
                                              VandeVondele2012,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_env,                     ONLY: cp_para_env_retain
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE input_constants,                 ONLY: &
        ls_cluster_atomic, ls_cluster_molecular, ls_s_inversion_hotelling, ls_s_inversion_none, &
        ls_s_inversion_sign_sqrt, ls_s_preconditioner_atomic, ls_s_preconditioner_molecular, &
        ls_s_preconditioner_none, ls_scf_ns, ls_scf_pexsi, ls_scf_tc2, ls_scf_trs4
   USE input_enumeration_types,         ONLY: enum_i2c,&
                                              enumeration_type
   USE input_keyword_types,             ONLY: keyword_get,&
                                              keyword_type
   USE input_section_types,             ONLY: section_get_keyword,&
                                              section_release,&
                                              section_type,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_retain,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE molecule_types_new,              ONLY: molecule_of_atom,&
                                              molecule_type
   USE pao_main,                        ONLY: pao_init
   USE particle_types,                  ONLY: particle_type
   USE pexsi_methods,                   ONLY: pexsi_init_read_input
   USE pexsi_types,                     ONLY: lib_pexsi_init
   USE qs_density_mixing_types,         ONLY: create_mixing_section,&
                                              mixing_storage_create
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf_create'

   PUBLIC :: ls_scf_create

CONTAINS

! **************************************************************************************************
!> \brief Creation and basic initialization of the LS type.
!> \param qs_env ...
!> \par History
!>       2012.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE ls_scf_create(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_create', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, unit_nr
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ls_scf_env_type), POINTER                     :: ls_scf_env
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: input, mixing_section

      CALL timeset(routineN, handle)

      CALL cite_reference(VandeVondele2012)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ALLOCATE (ls_scf_env)

      ! get basic quantities from the qs_env
      CALL get_qs_env(qs_env, nelectron_total=ls_scf_env%nelectron_total, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      molecule_set=molecule_set, &
                      input=input, &
                      has_unit_metric=ls_scf_env%has_unit_metric, &
                      para_env=ls_scf_env%para_env, &
                      do_transport=ls_scf_env%do_transport, &
                      nelectron_spin=ls_scf_env%nelectron_spin)

      ! copy some basic stuff
      ls_scf_env%nspins = dft_control%nspins
      ls_scf_env%natoms = SIZE(particle_set, 1)
      CALL cp_para_env_retain(ls_scf_env%para_env)

      ! initialize block to group to defined molecules
      ALLOCATE (ls_scf_env%ls_mstruct%atom_to_molecule(ls_scf_env%natoms))

      CALL molecule_of_atom(molecule_set, atom_to_mol=ls_scf_env%ls_mstruct%atom_to_molecule)

      ! parse the ls_scf section and set derived quantities
      CALL ls_scf_init_read_write_input(input, ls_scf_env, unit_nr)

      ! set up the buffer for the history of matrices
      ls_scf_env%scf_history%nstore = ls_scf_env%extrapolation_order
      ls_scf_env%scf_history%istore = 0
      ALLOCATE (ls_scf_env%scf_history%matrix(ls_scf_env%nspins, ls_scf_env%scf_history%nstore))

      NULLIFY (ls_scf_env%mixing_store)
      mixing_section => section_vals_get_subs_vals(input, "DFT%LS_SCF%RHO_MIXING")
      CALL mixing_storage_create(ls_scf_env%mixing_store, mixing_section, &
                                 ls_scf_env%density_mixing_method, &
                                 dft_control%qs_control%cutoff)

      ! initialize PEXSI
      IF (ls_scf_env%do_pexsi) THEN
         IF (dft_control%qs_control%eps_filter_matrix .NE. 0.0_dp) &
            CPABORT("EPS_FILTER_MATRIX must be set to 0 for PEXSI.")
         CALL lib_pexsi_init(ls_scf_env%pexsi, ls_scf_env%para_env%group, ls_scf_env%nspins)
      ENDIF

      ! initialize PAO
      CALL pao_init(qs_env, ls_scf_env)

      ! put the ls_scf_env in qs_env
      CALL set_qs_env(qs_env, ls_scf_env=ls_scf_env)

      CALL timestop(handle)

   END SUBROUTINE ls_scf_create

! **************************************************************************************************
!> \brief parse the input section, no need to pass it around
!> \param input ...
!> \param ls_scf_env ...
!> \param unit_nr ...
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE ls_scf_init_read_write_input(input, ls_scf_env, unit_nr)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(ls_scf_env_type), INTENT(INOUT)               :: ls_scf_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_init_read_write_input', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: mu
      TYPE(enumeration_type), POINTER                    :: enum
      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: section
      TYPE(section_vals_type), POINTER                   :: chebyshev_section, curvy_section, &
                                                            ls_scf_section, mixing_section, &
                                                            pao_section, pexsi_section

      CALL timeset(routineN, handle)
      CALL cite_reference(VandeVondele2012)
      ls_scf_section => section_vals_get_subs_vals(input, "DFT%LS_SCF")
      curvy_section => section_vals_get_subs_vals(ls_scf_section, "CURVY_STEPS")

      ! should come from input
      CALL section_vals_val_get(ls_scf_section, "LS_DIIS", l_val=ls_scf_env%ls_diis)
      CALL section_vals_val_get(ls_scf_section, "INI_DIIS", i_val=ls_scf_env%iter_ini_diis)
      CALL section_vals_val_get(ls_scf_section, "MAX_DIIS", i_val=ls_scf_env%max_diis)
      CALL section_vals_val_get(ls_scf_section, "NMIXING", i_val=ls_scf_env%nmixing)
      CALL section_vals_val_get(ls_scf_section, "EPS_DIIS", r_val=ls_scf_env%eps_diis)
      CALL section_vals_val_get(ls_scf_section, "EPS_SCF", r_val=ls_scf_env%eps_scf)
      CALL section_vals_val_get(ls_scf_section, "EPS_FILTER", r_val=ls_scf_env%eps_filter)
      CALL section_vals_val_get(ls_scf_section, "MU", r_val=mu)
      CALL section_vals_val_get(ls_scf_section, "FIXED_MU", l_val=ls_scf_env%fixed_mu)
      ls_scf_env%mu_spin = mu
      CALL section_vals_val_get(ls_scf_section, "MIXING_FRACTION", r_val=ls_scf_env%mixing_fraction)
      CALL section_vals_val_get(ls_scf_section, "MAX_SCF", i_val=ls_scf_env%max_scf)
      CALL section_vals_val_get(ls_scf_section, "S_PRECONDITIONER", i_val=ls_scf_env%s_preconditioner_type)
      CALL section_vals_val_get(ls_scf_section, "MATRIX_CLUSTER_TYPE", i_val=ls_scf_env%ls_mstruct%cluster_type)
      CALL section_vals_val_get(ls_scf_section, "SINGLE_PRECISION_MATRICES", l_val=ls_scf_env%ls_mstruct%single_precision)
      CALL section_vals_val_get(ls_scf_section, "S_INVERSION", i_val=ls_scf_env%s_inversion_type)
      CALL section_vals_val_get(ls_scf_section, "CHECK_S_INV", l_val=ls_scf_env%check_s_inv)
      CALL section_vals_val_get(ls_scf_section, "REPORT_ALL_SPARSITIES", l_val=ls_scf_env%report_all_sparsities)
      CALL section_vals_val_get(ls_scf_section, "PERFORM_MU_SCAN", l_val=ls_scf_env%perform_mu_scan)
      CALL section_vals_val_get(ls_scf_section, "PURIFICATION_METHOD", i_val=ls_scf_env%purification_method)
      CALL section_vals_val_get(ls_scf_section, "DYNAMIC_THRESHOLD", l_val=ls_scf_env%dynamic_threshold)
      CALL section_vals_val_get(ls_scf_section, "NON_MONOTONIC", l_val=ls_scf_env%non_monotonic)
      CALL section_vals_val_get(ls_scf_section, "SIGN_SQRT_ORDER", i_val=ls_scf_env%sign_sqrt_order)
      CALL section_vals_val_get(ls_scf_section, "EXTRAPOLATION_ORDER", i_val=ls_scf_env%extrapolation_order)
      CALL section_vals_val_get(ls_scf_section, "RESTART_READ", l_val=ls_scf_env%restart_read)
      CALL section_vals_val_get(ls_scf_section, "RESTART_WRITE", l_val=ls_scf_env%restart_write)
      CALL section_vals_val_get(ls_scf_section, "EPS_LANCZOS", r_val=ls_scf_env%eps_lanczos)
      CALL section_vals_val_get(ls_scf_section, "MAX_ITER_LANCZOS", i_val=ls_scf_env%max_iter_lanczos)

      CALL section_vals_get(curvy_section, explicit=ls_scf_env%curvy_steps)
      CALL section_vals_val_get(curvy_section, "LINE_SEARCH", i_val=ls_scf_env%curvy_data%line_search_type)
      CALL section_vals_val_get(curvy_section, "N_BCH_HISTORY", i_val=ls_scf_env%curvy_data%n_bch_hist)
      CALL section_vals_val_get(curvy_section, "MIN_HESSIAN_SHIFT", r_val=ls_scf_env%curvy_data%min_shift)
      CALL section_vals_val_get(curvy_section, "FILTER_FACTOR", r_val=ls_scf_env%curvy_data%filter_factor)
      CALL section_vals_val_get(curvy_section, "FILTER_FACTOR_SCALE", r_val=ls_scf_env%curvy_data%scale_filter)
      CALL section_vals_val_get(curvy_section, "MIN_FILTER", r_val=ls_scf_env%curvy_data%min_filter)

      ls_scf_env%extrapolation_order = MAX(0, ls_scf_env%extrapolation_order)

      chebyshev_section => section_vals_get_subs_vals(input, "DFT%LS_SCF%CHEBYSHEV")
      CALL section_vals_get(chebyshev_section, explicit=ls_scf_env%chebyshev%compute_chebyshev)
      IF (ls_scf_env%chebyshev%compute_chebyshev) THEN
         CALL section_vals_val_get(chebyshev_section, "N_CHEBYSHEV", i_val=ls_scf_env%chebyshev%n_chebyshev)
         CALL section_vals_val_get(chebyshev_section, "DOS%N_GRIDPOINTS", i_val=ls_scf_env%chebyshev%n_gridpoint_dos)

         ls_scf_env%chebyshev%print_key_dos => &
            section_vals_get_subs_vals(chebyshev_section, "DOS")
         CALL section_vals_retain(ls_scf_env%chebyshev%print_key_dos)

         ls_scf_env%chebyshev%print_key_cube => &
            section_vals_get_subs_vals(chebyshev_section, "PRINT_SPECIFIC_E_DENSITY_CUBE")
         CALL section_vals_retain(ls_scf_env%chebyshev%print_key_cube)
      ENDIF

      mixing_section => section_vals_get_subs_vals(input, "DFT%LS_SCF%RHO_MIXING")
      CALL section_vals_get(mixing_section, explicit=ls_scf_env%do_rho_mixing)

      CALL section_vals_val_get(mixing_section, "METHOD", i_val=ls_scf_env%density_mixing_method)
      IF (ls_scf_env%ls_diis .AND. ls_scf_env%do_rho_mixing) &
         CPABORT("LS_DIIS and RHO_MIXING are not compatible.")

      pexsi_section => section_vals_get_subs_vals(input, "DFT%LS_SCF%PEXSI")
      CALL section_vals_get(pexsi_section)

      ls_scf_env%do_pexsi = ls_scf_env%purification_method .EQ. ls_scf_pexsi &
                            .AND. .NOT. ls_scf_env%do_transport
      IF (ls_scf_env%do_pexsi) THEN
         CALL pexsi_init_read_input(pexsi_section, ls_scf_env%pexsi)
         ! Turn off S inversion (not used for PEXSI).
         ! Methods such as purification must thus be avoided... which is OK, as the density matrix computed in pexsi is
         ! a finite temperature one, and thus not idempotent
         ls_scf_env%s_inversion_type = ls_s_inversion_none
         ! PEXSI needs the sparsity pattern of S to be fixed by the upper bound ~ Int[|phi_a||phi_b|],
         ! they can not be filtered based on magnitude, as small elements in S (e.g. symmetry) do not necessarily
         ! correspond to small elements in the density matrix, with non-zero contributions to the total density.
         ! the easiest way to make sure S is untouched is to set eps_filter to zero (which should be inconsequential,
         ! as a run based on pexsi should execute exactly zero multiplications)
         ls_scf_env%eps_filter = 0.0_dp
      ENDIF

      ! Turn off S inversion and set eps_filter to zero for transport
      IF (ls_scf_env%do_transport) THEN
         ls_scf_env%s_inversion_type = ls_s_inversion_none
         ls_scf_env%eps_filter = 0.0_dp
      END IF

      SELECT CASE (ls_scf_env%s_inversion_type)
      CASE (ls_s_inversion_sign_sqrt)
         ls_scf_env%needs_s_inv = .TRUE.
         ls_scf_env%use_s_sqrt = .TRUE.
      CASE (ls_s_inversion_hotelling)
         ls_scf_env%needs_s_inv = .TRUE.
         ls_scf_env%use_s_sqrt = .FALSE.
      CASE (ls_s_inversion_none)
         ls_scf_env%needs_s_inv = .FALSE.
         ls_scf_env%use_s_sqrt = .FALSE.
      CASE DEFAULT
         CPABORT("")
      END SELECT

      SELECT CASE (ls_scf_env%s_preconditioner_type)
      CASE (ls_s_preconditioner_none)
         ls_scf_env%has_s_preconditioner = .FALSE.
      CASE DEFAULT
         ls_scf_env%has_s_preconditioner = .TRUE.
      END SELECT

      ! verify some requirements for the curvy steps
      IF (ls_scf_env%curvy_steps .AND. ls_scf_env%do_pexsi) &
         CPABORT("CURVY_STEPS cannot be used together with PEXSI.")
      IF (ls_scf_env%curvy_steps .AND. ls_scf_env%do_transport) &
         CPABORT("CURVY_STEPS cannot be used together with TRANSPORT.")
      IF (ls_scf_env%curvy_steps .AND. ls_scf_env%has_s_preconditioner) &
         CPABORT("S Preconditioning not implemented in combination with CURVY_STEPS.")
      IF (ls_scf_env%curvy_steps .AND. .NOT. ls_scf_env%use_s_sqrt) &
         CPABORT("CURVY_STEPS requires the use of the sqrt inversion.")

      ! an undocumented feature ... allows for just doing the initial guess, no expensive stuff
      IF (ls_scf_env%max_scf < 0) THEN
         ls_scf_env%needs_s_inv = .FALSE.
         ls_scf_env%use_s_sqrt = .FALSE.
         ls_scf_env%has_s_preconditioner = .FALSE.
      ENDIF

      pao_section => section_vals_get_subs_vals(input, "DFT%LS_SCF%PAO")
      CALL section_vals_get(pao_section, explicit=ls_scf_env%do_pao)
      ls_scf_env%ls_mstruct%do_pao = ls_scf_env%do_pao

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 30), " Linear scaling SCF ", REPEAT("-", 29)
         WRITE (unit_nr, '(T2,A,T38,E20.3)') "eps_scf:", ls_scf_env%eps_scf
         WRITE (unit_nr, '(T2,A,T38,E20.3)') "eps_filter:", ls_scf_env%eps_filter
         IF (ls_scf_env%do_rho_mixing) THEN
            IF (ls_scf_env%density_mixing_method > 0) THEN
               NULLIFY (section)
               CALL create_mixing_section(section, ls_scf=.TRUE.)
               keyword => section_get_keyword(section, "METHOD")
               CALL keyword_get(keyword, enum=enum)
               WRITE (unit_nr, "(T2,A,T38,A20)") &
                  "Density mixing in g-space:", ADJUSTR(TRIM(enum_i2c(enum, &
                                                                      ls_scf_env%density_mixing_method)))
               CALL section_release(section)
            END IF
         ELSE
            WRITE (unit_nr, '(T2,A,T38,E20.3)') "mixing_fraction:", ls_scf_env%mixing_fraction
         ENDIF
         WRITE (unit_nr, '(T2,A,T38,I20)') "max_scf:", ls_scf_env%max_scf
         IF (ls_scf_env%ls_diis) THEN
            WRITE (unit_nr, '(T2,A,T38,I20)') "DIIS: max_diis:", ls_scf_env%max_diis
            WRITE (unit_nr, '(T2,A,T38,E20.3)') "DIIS: eps_diis:", ls_scf_env%eps_diis
            WRITE (unit_nr, '(T2,A,T38,I20)') "DIIS: ini_diis:", ls_scf_env%iter_ini_diis
            WRITE (unit_nr, '(T2,A,T38,I20)') "DIIS: nmixing:", ls_scf_env%nmixing
         ENDIF
         WRITE (unit_nr, '(T2,A,T38,L20)') "fixed chemical potential (mu)", ls_scf_env%fixed_mu
         WRITE (unit_nr, '(T2,A,T38,L20)') "has unit metric:", ls_scf_env%has_unit_metric
         WRITE (unit_nr, '(T2,A,T38,L20)') "Computing inv(S):", ls_scf_env%needs_s_inv
         WRITE (unit_nr, '(T2,A,T38,L20)') "Computing sqrt(S):", ls_scf_env%use_s_sqrt
         WRITE (unit_nr, '(T2,A,T38,L20)') "Computing s preconditioner ", ls_scf_env%has_s_preconditioner
         WRITE (unit_nr, '(T2,A,T38,I20)') "sign sqrt order:", ls_scf_env%sign_sqrt_order
         WRITE (unit_nr, '(T2,A,T38,I20)') "Extrapolation order:", ls_scf_env%extrapolation_order
         WRITE (unit_nr, '(T2,A,T38,L20)') "Use single precision matrices", ls_scf_env%ls_mstruct%single_precision

         SELECT CASE (ls_scf_env%s_preconditioner_type)
         CASE (ls_s_preconditioner_none)
            WRITE (unit_nr, '(T2,A,T38,A20)') "S preconditioner type ", "NONE"
         CASE (ls_s_preconditioner_atomic)
            WRITE (unit_nr, '(T2,A,T38,A20)') "S preconditioner type ", "ATOMIC"
         CASE (ls_s_preconditioner_molecular)
            WRITE (unit_nr, '(T2,A,T38,A20)') "S preconditioner type ", "MOLECULAR"
         END SELECT

         WRITE (unit_nr, '(T2,A,T38,L20)') "Polarized Atomic Orbitals (PAO) ", ls_scf_env%do_pao

         IF (ls_scf_env%curvy_steps) THEN
            WRITE (unit_nr, '(T2,A,T38,A30)') "Using curvy steps to optimize the density matrix"
            CALL cite_reference(Shao2003)
         ENDIF

         SELECT CASE (ls_scf_env%purification_method)
         CASE (ls_scf_ns)
            WRITE (unit_nr, '(T2,A,T38,A30)') "Purification method", "Newton-Schulz sign iter"
         CASE (ls_scf_tc2)
            CALL cite_reference(Niklasson2014)
            WRITE (unit_nr, '(T2,A,T38,A30)') "Purification method", "Trace conserving 2nd  order"
         CASE (ls_scf_trs4)
            CALL cite_reference(Niklasson2003)
            WRITE (unit_nr, '(T2,A,T38,A30)') "Purification method", "Trace resetting 4th order"
         CASE (ls_scf_pexsi)
            CALL cite_reference(Lin2009)
            CALL cite_reference(Lin2013)
            WRITE (unit_nr, '(T2,A,T38,A20)') "Purification method", "PEXSI"
         CASE DEFAULT
            CPABORT("")
         END SELECT

         SELECT CASE (ls_scf_env%ls_mstruct%cluster_type)
         CASE (ls_cluster_atomic)
            WRITE (unit_nr, '(T2,A,T38,A20)') "Cluster type", "ATOMIC"
         CASE (ls_cluster_molecular)
            WRITE (unit_nr, '(T2,A,T38,A20)') "Cluster type", "MOLECULAR"
         CASE DEFAULT
            CPABORT("Unkown cluster type")
         END SELECT

         IF (ls_scf_env%chebyshev%compute_chebyshev) THEN
            WRITE (unit_nr, '(T2,A,T38,A20)') "Computing Chebyshev", "TRUE"
            WRITE (unit_nr, '(T2,A,T38,I20)') "N_CHEBYSHEV:", ls_scf_env%chebyshev%n_chebyshev
            WRITE (unit_nr, '(T2,A,T38,I20)') "N_GRIDPOINT_DOS:", ls_scf_env%chebyshev%n_gridpoint_dos
         ELSE
            WRITE (unit_nr, '(T2,A,T38,A20)') "Computing Chebyshev", "FALSE"
         ENDIF

         WRITE (unit_nr, '(T2,A,T38,L20)') "Using PAO", ls_scf_env%do_pao

         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
         WRITE (unit_nr, '()')
         CALL m_flush(unit_nr)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE ls_scf_init_read_write_input

END MODULE dm_ls_scf_create
