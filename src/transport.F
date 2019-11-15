!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines for DFT+NEGF calculations (coupling with the quantum transport code OMEN)
!> \par History
!>       12.2012 created external_scf_method [Hossein Bani-Hashemian]
!>       05.2013 created rotines to work with C-interoperable matrices [Hossein Bani-Hashemian]
!>       07.2013 created transport_env routines [Hossein Bani-Hashemian]
!>       11.2014 switch to CSR matrices [Hossein Bani-Hashemian]
!>       12.2014 merged [Hossein Bani-Hashemian]
!>       04.2016 added imaginary DM and current density cube [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
MODULE transport
   USE ISO_C_BINDING,                   ONLY: C_ASSOCIATED,&
                                              C_DOUBLE,&
                                              C_F_PROCPOINTER,&
                                              C_INT,&
                                              C_LOC,&
                                              C_NULL_PTR,&
                                              C_PTR
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE bibliography,                    ONLY: Bruck2014,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_to_csr_screening
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE dbcsr_api,                       ONLY: &
        dbcsr_convert_csr_to_dbcsr, dbcsr_convert_dbcsr_to_csr, dbcsr_copy, &
        dbcsr_copy_into_existing, dbcsr_create, dbcsr_csr_create, dbcsr_csr_create_from_dbcsr, &
        dbcsr_csr_dbcsr_blkrow_dist, dbcsr_csr_print_sparsity, dbcsr_csr_type, &
        dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_has_symmetry, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry, dbcsr_type_real_8
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_allgather,&
                                              mp_environ
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: boltzmann,&
                                              e_charge,&
                                              evolt,&
                                              h_bar
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_current,               ONLY: calculate_jrho_resp
   USE qs_linres_types,                 ONLY: current_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE transport_env_types,             ONLY: cp2k_csr_interop_type,&
                                              cp2k_transport_parameters,&
                                              csr_interop_matrix_get_info,&
                                              csr_interop_nullify,&
                                              transport_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'transport'

   PUBLIC :: transport_env_create, transport_initialize, external_scf_method
   PUBLIC :: qs_scf_post_transport

!> interface between C/C++ and FORTRAN
   INTERFACE c_func_interface
! **************************************************************************************************
!> \brief C routine that takes the S and H matrices as input and outputs a P matrix
!> \param cp2k_transport_params transport parameters read form a CP2K input file
!> \param s_mat  C-interoperable overlap matrix
!> \param ks_mat C-interoperable Kohn-Sham matrix
!> \param p_mat  C-interoperable density matrix
!> \param imagp_mat C-interoperable imaginary density matrix
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
      SUBROUTINE c_scf_routine(cp2k_transport_params, s_mat, ks_mat, p_mat, imagp_mat) BIND(C)
         IMPORT :: C_INT, C_PTR, cp2k_csr_interop_type, cp2k_transport_parameters
         IMPLICIT NONE
         TYPE(cp2k_transport_parameters), VALUE, INTENT(IN) :: cp2k_transport_params
         TYPE(cp2k_csr_interop_type), VALUE, INTENT(IN)     :: s_mat
         TYPE(cp2k_csr_interop_type), VALUE, INTENT(IN)     :: ks_mat
         TYPE(cp2k_csr_interop_type), INTENT(INOUT)         :: p_mat
         TYPE(cp2k_csr_interop_type), INTENT(INOUT)         :: imagp_mat
      END SUBROUTINE c_scf_routine
   END INTERFACE c_func_interface

CONTAINS

! **************************************************************************************************
!> \brief creates the transport environment
!> \param[inout] qs_env the qs_env containing the transport_env
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE transport_env_create(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'transport_env_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(transport_env_type), POINTER                  :: transport_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      transport_env=transport_env, &
                      dft_control=dft_control, &
                      input=input)

      CPASSERT(.NOT. ASSOCIATED(transport_env))

      ALLOCATE (transport_env)

      CALL transport_init_read_input(input, transport_env)
      CALL transport_set_contact_params(qs_env, transport_env)

      CALL set_qs_env(qs_env, transport_env=transport_env)

      CALL timestop(handle)

   END SUBROUTINE transport_env_create

! **************************************************************************************************
!> \brief intitializes all fields of transport_env using the parameters read from
!>        the corresponding input section
!> \param[inout] input         the input file
!> \param[inout] transport_env the transport_env to be initialized
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE transport_init_read_input(input, transport_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(transport_env_type), INTENT(INOUT)            :: transport_env

      CHARACTER(len=*), PARAMETER :: routineN = 'transport_init_read_input', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: contact_bandwidth, contact_injsign, &
                                                            contact_natoms, contact_start, handle, &
                                                            i, n_contacts, stride_contacts
      INTEGER, DIMENSION(:), POINTER                     :: i_vals
      LOGICAL                                            :: contact_explicit, injecting_contact, &
                                                            obc_equilibrium, one_circle
      TYPE(section_vals_type), POINTER                   :: beyn_section, contact_section, &
                                                            pexsi_section, transport_section

      CALL timeset(routineN, handle)

      transport_section => section_vals_get_subs_vals(input, "DFT%TRANSPORT")
      contact_section => section_vals_get_subs_vals(transport_section, "CONTACT")
      beyn_section => section_vals_get_subs_vals(transport_section, "BEYN")
      pexsi_section => section_vals_get_subs_vals(transport_section, "PEXSI")
      CALL section_vals_get(contact_section, explicit=contact_explicit, n_repetition=n_contacts)

      NULLIFY (i_vals)
! read from input
      CALL section_vals_val_get(transport_section, "TRANSPORT_METHOD", i_val=transport_env%params%method)
      CALL section_vals_val_get(transport_section, "INJECTION_METHOD", i_val=transport_env%params%injection_method)
      CALL section_vals_val_get(transport_section, "REAL_AXIS_INTEGRATION_METHOD", &
                                i_val=transport_env%params%rlaxis_integration_method)
      CALL section_vals_val_get(transport_section, "QT_FORMALISM", i_val=transport_env%params%qt_formalism)
      CALL section_vals_val_get(transport_section, "LINEAR_SOLVER", i_val=transport_env%params%linear_solver)
      CALL section_vals_val_get(transport_section, "MATRIX_INVERSION_METHOD", &
                                i_val=transport_env%params%matrixinv_method)
      CALL section_vals_val_get(transport_section, "CONTACT_FILLING", i_val=transport_env%params%transport_neutral)
      CALL section_vals_val_get(transport_section, "N_KPOINTS", i_val=transport_env%params%n_kpoint)
      CALL section_vals_val_get(transport_section, "NUM_INTERVAL", i_val=transport_env%params%num_interval)
      CALL section_vals_val_get(transport_section, "TASKS_PER_ENERGY_POINT", &
                                i_val=transport_env%params%tasks_per_energy_point)
      CALL section_vals_val_get(transport_section, "TASKS_PER_POLE", i_val=transport_env%params%tasks_per_pole)
      CALL section_vals_val_get(transport_section, "NUM_POLE", i_val=transport_env%params%num_pole)
      CALL section_vals_val_get(transport_section, "GPUS_PER_POINT", i_val=transport_env%params%gpus_per_point)
      CALL section_vals_val_get(transport_section, "N_POINTS_INV", i_val=transport_env%params%n_points_inv)
      CALL section_vals_val_get(transport_section, "COLZERO_THRESHOLD", r_val=transport_env%params%colzero_threshold)
      CALL section_vals_val_get(transport_section, "EPS_LIMIT", r_val=transport_env%params%eps_limit)
      CALL section_vals_val_get(transport_section, "EPS_LIMIT_CC", r_val=transport_env%params%eps_limit_cc)
      CALL section_vals_val_get(transport_section, "EPS_DECAY", r_val=transport_env%params%eps_decay)
      CALL section_vals_val_get(transport_section, "EPS_SINGULARITY_CURVATURES", &
                                r_val=transport_env%params%eps_singularity_curvatures)
      CALL section_vals_val_get(transport_section, "EPS_MU", r_val=transport_env%params%eps_mu)
      CALL section_vals_val_get(transport_section, "EPS_EIGVAL_DEGEN", r_val=transport_env%params%eps_eigval_degen)
      CALL section_vals_val_get(transport_section, "EPS_FERMI", r_val=transport_env%params%eps_fermi)
      CALL section_vals_val_get(transport_section, "ENERGY_INTERVAL", r_val=transport_env%params%energy_interval)
      CALL section_vals_val_get(transport_section, "MIN_INTERVAL", r_val=transport_env%params%min_interval)
      CALL section_vals_val_get(transport_section, "TEMPERATURE", r_val=transport_env%params%temperature)
      CALL section_vals_val_get(transport_section, "DENSITY_MIXING", r_val=transport_env%params%dens_mixing)
      CALL section_vals_val_get(transport_section, "CSR_SCREENING", l_val=transport_env%csr_screening)

      ! logical*1 to logical*4 , l_val is logical*1 and c_bool is equivalent to logical*4
      CALL section_vals_val_get(transport_section, "OBC_EQUILIBRIUM", l_val=obc_equilibrium)
      IF (obc_equilibrium) THEN
         transport_env%params%obc_equilibrium = .TRUE.
      ELSE
         transport_env%params%obc_equilibrium = .FALSE.
      END IF

      CALL section_vals_val_get(transport_section, "CUTOUT", i_vals=i_vals)
      transport_env%params%cutout = i_vals

      CALL section_vals_val_get(beyn_section, "TASKS_PER_INTEGRATION_POINT", &
                                i_val=transport_env%params%tasks_per_integration_point)
      CALL section_vals_val_get(beyn_section, "N_POINTS_BEYN", i_val=transport_env%params%n_points_beyn)
      CALL section_vals_val_get(beyn_section, "N_RAND", r_val=transport_env%params%n_rand_beyn)
      CALL section_vals_val_get(beyn_section, "N_RAND_CC", r_val=transport_env%params%n_rand_cc_beyn)
      CALL section_vals_val_get(beyn_section, "SVD_CUTOFF", r_val=transport_env%params%svd_cutoff)
      CALL section_vals_val_get(beyn_section, "ONE_CIRCLE", l_val=one_circle)
      IF (one_circle) THEN
         transport_env%params%ncrc_beyn = 1
      ELSE
         transport_env%params%ncrc_beyn = 2
      END IF

      CALL section_vals_val_get(pexsi_section, "ORDERING", i_val=transport_env%params%ordering)
      CALL section_vals_val_get(pexsi_section, "ROW_ORDERING", i_val=transport_env%params%row_ordering)
      CALL section_vals_val_get(pexsi_section, "VERBOSITY", i_val=transport_env%params%verbosity)
      CALL section_vals_val_get(pexsi_section, "NP_SYMB_FACT", i_val=transport_env%params%pexsi_np_symb_fact)

      IF (contact_explicit) THEN
         transport_env%params%num_contacts = n_contacts
         stride_contacts = 5
         transport_env%params%stride_contacts = stride_contacts
         ALLOCATE (transport_env%contacts_data(stride_contacts*n_contacts))

         DO i = 1, n_contacts
            CALL section_vals_val_get(contact_section, "BANDWIDTH", i_rep_section=i, i_val=contact_bandwidth)
            CALL section_vals_val_get(contact_section, "START", i_rep_section=i, i_val=contact_start)
            CALL section_vals_val_get(contact_section, "N_ATOMS", i_rep_section=i, i_val=contact_natoms)
            CALL section_vals_val_get(contact_section, "INJECTION_SIGN", i_rep_section=i, i_val=contact_injsign)
            CALL section_vals_val_get(contact_section, "INJECTING_CONTACT", i_rep_section=i, l_val=injecting_contact)

            IF (contact_natoms .LE. 0) CPABORT("Number of atoms in contact region needs to be defined.")

            transport_env%contacts_data((i - 1)*stride_contacts + 1) = contact_bandwidth
            transport_env%contacts_data((i - 1)*stride_contacts + 2) = contact_start - 1 ! C indexing
            transport_env%contacts_data((i - 1)*stride_contacts + 3) = contact_natoms
            transport_env%contacts_data((i - 1)*stride_contacts + 4) = contact_injsign

            IF (injecting_contact) THEN
               transport_env%contacts_data((i - 1)*stride_contacts + 5) = 1
            ELSE
               transport_env%contacts_data((i - 1)*stride_contacts + 5) = 0
            END IF
         END DO
         transport_env%params%contacts_data = C_LOC(transport_env%contacts_data(1))
      ELSE
         CPABORT("No contact region is defined.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE transport_init_read_input

! **************************************************************************************************
!> \brief initializes the transport environment
!> \param ks_env ...
!> \param[inout] transport_env the transport env to be initialized
!> \param[in]    template_matrix   template matrix to keep the sparsity of matrices fixed
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE transport_initialize(ks_env, transport_env, template_matrix)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(transport_env_type), INTENT(INOUT)            :: transport_env
      TYPE(dbcsr_type), INTENT(IN)                       :: template_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'transport_initialize', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, mp_group, mynode, numnodes, &
                                                            unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      CALL cite_reference(Bruck2014)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      mp_group = logger%para_env%group
      CALL mp_environ(numnodes, mynode, mp_group)

      IF (dbcsr_has_symmetry(template_matrix)) THEN
         CALL dbcsr_copy(transport_env%template_matrix_sym, template_matrix)
         CALL dbcsr_desymmetrize(transport_env%template_matrix_sym, transport_env%template_matrix_nosym)
      ELSE
         CALL dbcsr_copy(transport_env%template_matrix_nosym, template_matrix)
         CALL dbcsr_copy(transport_env%template_matrix_sym, template_matrix)
      END IF

      ALLOCATE (transport_env%dm_imag)
      CALL dbcsr_create(transport_env%dm_imag, "imaginary DM", &
                        template=template_matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_set(transport_env%dm_imag, 0.0_dp)

      CALL dbcsr_create(transport_env%csr_sparsity, "CSR sparsity", &
                        template=transport_env%template_matrix_sym, &
                        data_type=dbcsr_type_real_8)
      CALL dbcsr_copy(transport_env%csr_sparsity, transport_env%template_matrix_sym)

      CALL cp_dbcsr_to_csr_screening(ks_env, transport_env%csr_sparsity)

      IF (.NOT. transport_env%csr_screening) CALL dbcsr_set(transport_env%csr_sparsity, 1.0)
      CALL dbcsr_csr_create_from_dbcsr(transport_env%template_matrix_nosym, &
                                       transport_env%s_matrix, &
                                       dbcsr_csr_dbcsr_blkrow_dist, &
                                       csr_sparsity=transport_env%csr_sparsity, &
                                       numnodes=numnodes)

      CALL dbcsr_csr_print_sparsity(transport_env%s_matrix, unit_nr)

      CALL dbcsr_convert_dbcsr_to_csr(transport_env%template_matrix_nosym, transport_env%s_matrix)

      CALL dbcsr_csr_create(transport_env%ks_matrix, transport_env%s_matrix)
      CALL dbcsr_csr_create(transport_env%p_matrix, transport_env%s_matrix)
      CALL dbcsr_csr_create(transport_env%imagp_matrix, transport_env%s_matrix)

      CALL timestop(handle)

   END SUBROUTINE transport_initialize

! **************************************************************************************************
!> \brief SCF calcualtion with an externally evaluated density matrix
!> \param[inout] transport_env  transport environment
!> \param[in]    matrix_s       DBCSR overlap matrix
!> \param[in]    matrix_ks      DBCSR Kohn-Sham matrix
!> \param[inout] matrix_p       DBCSR density matrix
!> \param[in]    nelectron_spin number of electrons
!> \param[in]    natoms         number of atoms
!> \param[in]    energy_diff    scf energy difference
!> \param[in]    iscf           the current scf iteration
!> \param[in]    extra_scf      whether or not an extra scf step will be performed
!> \par History
!>       12.2012 created [Hossein Bani-Hashemian]
!>       12.2014 revised [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE external_scf_method(transport_env, matrix_s, matrix_ks, matrix_p, &
                                  nelectron_spin, natoms, energy_diff, iscf, extra_scf)

      TYPE(transport_env_type), INTENT(INOUT)            :: transport_env
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s, matrix_ks
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_p
      INTEGER, INTENT(IN)                                :: nelectron_spin, natoms
      REAL(dp), INTENT(IN)                               :: energy_diff
      INTEGER, INTENT(IN)                                :: iscf
      LOGICAL, INTENT(IN)                                :: extra_scf

      CHARACTER(len=*), PARAMETER :: routineN = 'external_scf_method', &
         routineP = moduleN//':'//routineN

      TYPE(cp2k_csr_interop_type)                        :: imagp_mat, ks_mat, p_mat, s_mat

      PROCEDURE(c_scf_routine), POINTER        :: c_method
      INTEGER                                  :: handle

      CALL timeset(routineN, handle)

      CALL C_F_PROCPOINTER(transport_env%ext_c_method_ptr, c_method)
      IF (.NOT. C_ASSOCIATED(transport_env%ext_c_method_ptr)) &
         CALL cp_abort(__LOCATION__, &
                       "MISSING C/C++ ROUTINE: The TRANSPORT section is meant to be used together with an external "// &
                       "program, e.g. the quantum transport code OMEN, that provides CP2K with a density matrix.")

      transport_env%params%n_occ = nelectron_spin
      transport_env%params%n_atoms = natoms
      transport_env%params%energy_diff = energy_diff
      transport_env%params%evoltfactor = evolt
      transport_env%params%e_charge = e_charge
      transport_env%params%boltzmann = boltzmann
      transport_env%params%h_bar = h_bar
      transport_env%params%iscf = iscf
      transport_env%params%extra_scf = extra_scf

      CALL csr_interop_nullify(s_mat)
      CALL csr_interop_nullify(ks_mat)
      CALL csr_interop_nullify(p_mat)
      CALL csr_interop_nullify(imagp_mat)

      CALL dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_s)
      CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%s_matrix, s_mat)

      CALL dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_ks)
      CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%ks_matrix, ks_mat)

      CALL dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_p)
      CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%p_matrix, p_mat)

      CALL dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_s)
      CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%imagp_matrix, imagp_mat)

      CALL c_method(transport_env%params, s_mat, ks_mat, p_mat, imagp_mat)

      CALL convert_csr_interop_to_dbcsr(p_mat, transport_env%p_matrix, transport_env%template_matrix_nosym)
      CALL dbcsr_copy(matrix_p, transport_env%template_matrix_nosym)

      CALL convert_csr_interop_to_dbcsr(imagp_mat, transport_env%imagp_matrix, transport_env%template_matrix_nosym)
      CALL dbcsr_copy(transport_env%dm_imag, transport_env%template_matrix_nosym)

      CALL timestop(handle)

   END SUBROUTINE external_scf_method

! **************************************************************************************************
!> \brief converts a DBCSR matrix to a C-interoperable CSR matrix
!> \param[in]    dbcsr_mat  DBCSR matrix to be converted
!> \param[inout] csr_mat    auxiliary CSR matrix
!> \param[inout] csr_interop_mat C-interoperable CSR matrix
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE convert_dbcsr_to_csr_interop(dbcsr_mat, csr_mat, csr_interop_mat)

      TYPE(dbcsr_type), INTENT(IN)          :: dbcsr_mat
      TYPE(dbcsr_csr_type), INTENT(INOUT)            :: csr_mat
      TYPE(cp2k_csr_interop_type), INTENT(INOUT)      :: csr_interop_mat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_dbcsr_to_csr_interop', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle, mepos, num_pe, mp_group
      INTEGER, ALLOCATABLE, DIMENSION(:)       :: nrows_local_all, first_row_all
      INTEGER(C_INT), DIMENSION(:), POINTER    :: colind_local, rowptr_local, nzerow_local
      REAL(C_DOUBLE), DIMENSION(:), POINTER    :: nzvals_local
      TYPE(cp_logger_type), POINTER            :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      mp_group = logger%para_env%group
      mepos = logger%para_env%mepos
      num_pe = logger%para_env%num_pe

! dbcsr to csr
      CALL dbcsr_convert_dbcsr_to_csr(dbcsr_mat, csr_mat)

! csr to csr_interop
      rowptr_local => csr_mat%rowptr_local
      colind_local => csr_mat%colind_local
      nzerow_local => csr_mat%nzerow_local
      nzvals_local => csr_mat%nzval_local%r_dp ! support real double percision for now

      IF (SIZE(rowptr_local) .EQ. 0) THEN
         csr_interop_mat%rowptr_local = C_NULL_PTR
      ELSE
         csr_interop_mat%rowptr_local = C_LOC(rowptr_local(1))
      END IF

      IF (SIZE(colind_local) .EQ. 0) THEN
         csr_interop_mat%colind_local = C_NULL_PTR
      ELSE
         csr_interop_mat%colind_local = C_LOC(colind_local(1))
      END IF

      IF (SIZE(nzerow_local) .EQ. 0) THEN
         csr_interop_mat%nzerow_local = C_NULL_PTR
      ELSE
         csr_interop_mat%nzerow_local = C_LOC(nzerow_local(1))
      END IF

      IF (SIZE(nzvals_local) .EQ. 0) THEN
         csr_interop_mat%nzvals_local = C_NULL_PTR
      ELSE
         csr_interop_mat%nzvals_local = C_LOC(nzvals_local(1))
      END IF

      ALLOCATE (nrows_local_all(0:num_pe - 1), first_row_all(0:num_pe - 1))
      CALL mp_allgather(csr_mat%nrows_local, nrows_local_all, mp_group)
      CALL cumsum_i(nrows_local_all, first_row_all)

      IF (mepos .EQ. 0) THEN
         csr_interop_mat%first_row = 0
      ELSE
         csr_interop_mat%first_row = first_row_all(mepos - 1)
      END IF
      csr_interop_mat%nrows_total = csr_mat%nrows_total
      csr_interop_mat%ncols_total = csr_mat%ncols_total
      csr_interop_mat%nze_local = csr_mat%nze_local
      IF (csr_mat%nze_total > HUGE(csr_interop_mat%nze_total)) THEN
         CPABORT("overflow in nze")
      ENDIF
      csr_interop_mat%nze_total = INT(csr_mat%nze_total, KIND=KIND(csr_interop_mat%nze_total))
      csr_interop_mat%nrows_local = csr_mat%nrows_local
      csr_interop_mat%data_type = csr_mat%nzval_local%data_type

      CALL timestop(handle)

   CONTAINS
! **************************************************************************************************
!> \brief cumulative sum of a 1d array of integers
!> \param[in]  arr    input array
!> \param[out] cumsum cumulative sum of the input array
! **************************************************************************************************
      SUBROUTINE cumsum_i(arr, cumsum)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: arr
      INTEGER, DIMENSION(SIZE(arr)), INTENT(OUT)         :: cumsum

      INTEGER                                            :: i

         cumsum(1) = arr(1)
         DO i = 2, SIZE(arr)
            cumsum(i) = cumsum(i - 1) + arr(i)
         END DO
      END SUBROUTINE cumsum_i

   END SUBROUTINE convert_dbcsr_to_csr_interop

! **************************************************************************************************
!> \brief converts a C-interoperable CSR matrix to a DBCSR matrix
!> \param[in] csr_interop_mat C-interoperable CSR matrix
!> \param[inout] csr_mat         auxiliary CSR matrix
!> \param[inout] dbcsr_mat       DBCSR matrix
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE convert_csr_interop_to_dbcsr(csr_interop_mat, csr_mat, dbcsr_mat)

      TYPE(cp2k_csr_interop_type), INTENT(IN)            :: csr_interop_mat
      TYPE(dbcsr_csr_type), INTENT(INOUT)                      :: csr_mat
      TYPE(dbcsr_type), INTENT(INOUT)                    :: dbcsr_mat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'convert_csr_interop_to_dbcsr', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: data_type, handle, ncols_total, &
                                                            nrows_local, nrows_total, nze_local, &
                                                            nze_total
      INTEGER, DIMENSION(:), POINTER                     :: colind_local, nzerow_local, rowptr_local
      REAL(dp), DIMENSION(:), POINTER                    :: nzvals_local

      CALL timeset(routineN, handle)

! csr_interop to csr
      CALL csr_interop_matrix_get_info(csr_interop_mat, &
                                       nrows_total=nrows_total, ncols_total=ncols_total, nze_local=nze_local, &
                                       nze_total=nze_total, nrows_local=nrows_local, data_type=data_type, &
                                       rowptr_local=rowptr_local, colind_local=colind_local, &
                                       nzerow_local=nzerow_local, nzvals_local=nzvals_local)

      csr_mat%nrows_total = nrows_total
      csr_mat%ncols_total = ncols_total
      csr_mat%nze_local = nze_local
      csr_mat%nze_total = nze_total
      csr_mat%nrows_local = nrows_local
      csr_mat%nzval_local%data_type = data_type

      csr_mat%rowptr_local = rowptr_local
      csr_mat%colind_local = colind_local
      csr_mat%nzerow_local = nzerow_local
      csr_mat%nzval_local%r_dp = nzvals_local

! csr to dbcsr
      CALL dbcsr_convert_csr_to_dbcsr(dbcsr_mat, csr_mat)

      CALL timestop(handle)

   END SUBROUTINE convert_csr_interop_to_dbcsr

! **************************************************************************************************
!> \brief extraxts zeff (effective nuclear charges per atom) and nsgf (the size
!>   of spherical Gaussian basis functions per atom) from qs_env and intializes
!>   the corresponding arrays in transport_env%params
!> \param[in] qs_env qs environment
!> \param[inout] transport_env transport environment
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE transport_set_contact_params(qs_env, transport_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(transport_env_type), INTENT(INOUT)            :: transport_env

      CHARACTER(len=*), PARAMETER :: routineN = 'transport_set_contact_params', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iat, ikind, natom, nkind
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      REAL(KIND=dp)                                      :: zeff
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)
      CALL get_qs_env(qs_env, particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      ALLOCATE (transport_env%nsgf(natom))
      ALLOCATE (transport_env%zeff(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=transport_env%nsgf)

      ! reference charges
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind, atom_list=atom_list)
         DO iat = 1, SIZE(atom_list)
            i = atom_list(iat)
            transport_env%zeff(i) = zeff
         END DO
      END DO

      IF (natom .EQ. 0) THEN
         transport_env%params%nsgf = C_NULL_PTR
         transport_env%params%zeff = C_NULL_PTR
      ELSE
         transport_env%params%nsgf = C_LOC(transport_env%nsgf(1))
         transport_env%params%zeff = C_LOC(transport_env%zeff(1))
      END IF

   END SUBROUTINE transport_set_contact_params

!**************************************************************************************************
!> \brief evaluates current density using the imaginary part of the density matrix and prints it
!>        into cube files
!> \param[in] input the input
!> \param[in] qs_env qs environment
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE transport_current(input, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'transport_current', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=14)                                  :: ext
      CHARACTER(len=2)                                   :: sdir
      INTEGER                                            :: dir, handle, unit_nr
      LOGICAL                                            :: do_current_cube, do_transport, mpi_io
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(current_env_type)                             :: current_env
      TYPE(dbcsr_type), POINTER                          :: zero
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r
      TYPE(pw_p_type), POINTER                           :: gs, rs
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section
      TYPE(transport_env_type), POINTER                  :: transport_env

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      dft_section => section_vals_get_subs_vals(input, "DFT")
      CALL get_qs_env(qs_env=qs_env, &
                      subsys=subsys, &
                      pw_env=pw_env, &
                      transport_env=transport_env, &
                      do_transport=do_transport, &
                      dft_control=dft_control, &
                      rho=rho)
      CALL qs_subsys_get(subsys, particles=particles)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g)

      do_current_cube = BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                         "DFT%TRANSPORT%PRINT%CURRENT"), cp_p_file)

      ! check if transport is active i.e. the imaginary dm has been actually passed to cp2k by the external code
      IF (do_transport) THEN
         IF (C_ASSOCIATED(transport_env%ext_c_method_ptr)) THEN

            ! calculate_jrho_resp uses sab_all which is not associated in DFTB environment
            IF (dft_control%qs_control%dftb) CPABORT("Current is not available for DFTB.")
            IF (dft_control%qs_control%xtb) CPABORT("Current is not available for xTB.")

            ! now do the current
            IF (do_current_cube) THEN
               current_env%gauge = -1
               current_env%gauge_init = .FALSE.

               NULLIFY (rs, gs)
               ALLOCATE (rs, gs)
               NULLIFY (rs%pw, gs%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, rs%pw, use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_pool_create_pw(auxbas_pw_pool, gs%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

               NULLIFY (zero)
               ALLOCATE (zero)
               CALL dbcsr_create(zero, template=transport_env%dm_imag)
               CALL dbcsr_copy(zero, transport_env%dm_imag)
               CALL dbcsr_set(zero, 0.0_dp)

               DO dir = 1, 3
                  CALL pw_zero(rs%pw)
                  CALL pw_zero(gs%pw)
                  CALL calculate_jrho_resp(zero, transport_env%dm_imag, &
                                           zero, zero, dir, dir, rs, gs, qs_env, current_env, &
                                           retain_rsgrid=.TRUE.)
                  SELECT CASE (dir)
                  CASE (1)
                     sdir = "-x"
                  CASE (2)
                     sdir = "-y"
                  CASE (3)
                     sdir = "-z"
                  END SELECT
                  ext = sdir//".cube"
                  mpi_io = .TRUE.
                  unit_nr = cp_print_key_unit_nr(logger, dft_section, "TRANSPORT%PRINT%CURRENT", &
                                                 extension=ext, file_status="REPLACE", file_action="WRITE", &
                                                 log_filename=.FALSE., mpi_io=mpi_io)
                  CALL cp_pw_to_cube(rs%pw, unit_nr, "Transport current", particles=particles, &
                                     stride=section_get_ivals(dft_section, "TRANSPORT%PRINT%CURRENT%STRIDE"), &
                                     mpi_io=mpi_io)
                  CALL cp_print_key_finished_output(unit_nr, logger, dft_section, "TRANSPORT%PRINT%CURRENT", &
                                                    mpi_io=mpi_io)
               END DO

               CALL dbcsr_deallocate_matrix(zero)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rs%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, gs%pw)
               DEALLOCATE (rs, gs)
            END IF
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE transport_current

!**************************************************************************************************
!> \brief post scf calculations for transport
!> \param[in] qs_env qs environment
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE qs_scf_post_transport(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_transport', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      NULLIFY (input)

      CALL get_qs_env(qs_env=qs_env, &
                      input=input)

      CALL transport_current(input, qs_env)

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_transport

END MODULE transport

