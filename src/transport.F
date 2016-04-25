!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines for DFT+NEGF calculations (coupling with the quantum transport code OMEN)
!> \par History
!>       12.2012 created external_scf_method [Hossein Bani-Hashemian]
!>       05.2013 created rotines to work with C-interoperable matrices [Hossein Bani-Hashemian]
!>       07.2013 created transport_env routines [Hossein Bani-Hashemian]
!>       11.2014 switch to CSR matrices [Hossein Bani-Hashemian]
!>       12.2014 merged [Hossein Bani-Hashemian]
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
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Bruck2014,&
                                              cite_reference
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_to_csr_screening
   USE cp_dbcsr_interface,              ONLY: &
        cp_convert_csr_to_dbcsr, cp_convert_dbcsr_to_csr, cp_csr_create, cp_csr_create_from_dbcsr, &
        cp_csr_print_sparsity, cp_dbcsr_copy, cp_dbcsr_copy_into_existing, cp_dbcsr_create, &
        cp_dbcsr_desymmetrize, cp_dbcsr_has_symmetry, cp_dbcsr_init, cp_dbcsr_set, cp_dbcsr_type, &
        csr_eqrow_ceil_dist, csr_type, dbcsr_type_real_8
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_constants,                 ONLY: transport_diagscalapack
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp,&
                                              int_1
   USE message_passing,                 ONLY: mp_allgather,&
                                              mp_bcast,&
                                              mp_environ
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              get_neighbor_list_set_p,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
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

!> interface between C/C++ and FORTRAN
   INTERFACE c_func_interface
! **************************************************************************************************
!> \brief C routine that takes the S and H matrices as input and outputs a P matrix
!> \param cp2k_transport_params transport parameters read form a CP2K input file
!> \param s_mat  C-interoperable overlap matrix
!> \param ks_mat C-interoperable Kohn-Sham matrix
!> \param p_mat  C-interoperable density matrix
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
      SUBROUTINE c_scf_routine(cp2k_transport_params, s_mat, ks_mat, p_mat) BIND(C)
         IMPORT :: C_INT, C_PTR, cp2k_csr_interop_type, cp2k_transport_parameters
         IMPLICIT NONE
         TYPE(cp2k_transport_parameters), VALUE, INTENT(IN) :: cp2k_transport_params
         TYPE(cp2k_csr_interop_type), VALUE, INTENT(IN)     :: s_mat
         TYPE(cp2k_csr_interop_type), VALUE, INTENT(IN)     :: ks_mat
         TYPE(cp2k_csr_interop_type), INTENT(INOUT)         :: p_mat
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
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(transport_env_type), POINTER                  :: transport_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, transport_env=transport_env, input=input)

      CPASSERT(.NOT. ASSOCIATED(transport_env))

      ALLOCATE (transport_env)
      CALL transport_init_read_input(input, transport_env)
      CALL transport_set_contact_params(qs_env, transport_env)
      CALL cp_dbcsr_init(transport_env%template_matrix_sym)
      CALL cp_dbcsr_init(transport_env%template_matrix_nosym)
      CALL cp_dbcsr_init(transport_env%csr_sparsity)
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
                                                            one_circle
      TYPE(section_vals_type), POINTER                   :: beyn_section, contact_section, &
                                                            transport_section

      CALL timeset(routineN, handle)

      transport_section => section_vals_get_subs_vals(input, "DFT%TRANSPORT")
      contact_section => section_vals_get_subs_vals(transport_section, "CONTACT")
      beyn_section => section_vals_get_subs_vals(transport_section, "BEYN")
      CALL section_vals_get(contact_section, explicit=contact_explicit, n_repetition=n_contacts)

      NULLIFY (i_vals)
! read from input
      CALL section_vals_val_get(transport_section, "TRANSPORT_METHOD", i_val=transport_env%params%method)
      CALL section_vals_val_get(transport_section, "INJECTION_METHOD", i_val=transport_env%params%injection_method)
      CALL section_vals_val_get(transport_section, "REAL_AXIS_INTEGRATION_METHOD", &
                                i_val=transport_env%params%rlaxis_integration_method)
      CALL section_vals_val_get(transport_section, "LINEAR_SOLVER", i_val=transport_env%params%linear_solver)
      CALL section_vals_val_get(transport_section, "N_ABSCISSAE", i_val=transport_env%params%n_abscissae)
      CALL section_vals_val_get(transport_section, "N_KPOINTS", i_val=transport_env%params%n_kpoint)
      CALL section_vals_val_get(transport_section, "NUM_INTERVAL", i_val=transport_env%params%num_interval)
      CALL section_vals_val_get(transport_section, "TASKS_PER_POINT", i_val=transport_env%params%tasks_per_point)
      CALL section_vals_val_get(transport_section, "GPUS_PER_POINT", i_val=transport_env%params%gpus_per_point)
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
      CALL section_vals_val_get(transport_section, "CSR_SCREENING", l_val=transport_env%csr_screening)

      CALL section_vals_val_get(transport_section, "CUTOUT", i_vals=i_vals)
      transport_env%params%cutout = i_vals

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

      IF (contact_explicit) THEN
         transport_env%params%num_contacts = n_contacts
         stride_contacts = 6; 
         ALLOCATE (transport_env%contacts_data(stride_contacts*n_contacts))

         DO i = 1, n_contacts
            CALL section_vals_val_get(contact_section, "BANDWIDTH", i_rep_section=i, i_val=contact_bandwidth)
            CALL section_vals_val_get(contact_section, "START", i_rep_section=i, i_val=contact_start)
            CALL section_vals_val_get(contact_section, "N_ATOMS", i_rep_section=i, i_val=contact_natoms)
            CALL section_vals_val_get(contact_section, "INJECTION_SIGN", i_rep_section=i, i_val=contact_injsign)
            CALL section_vals_val_get(contact_section, "INJECTING_CONTACT", i_rep_section=i, l_val=injecting_contact)

            IF (contact_natoms .LE. 0) CPABORT("Number of atoms in contact region needs to be defined.")

            transport_env%contacts_data((i-1)*stride_contacts+1) = contact_bandwidth
            transport_env%contacts_data((i-1)*stride_contacts+2) = contact_start-1 ! C indexing
            transport_env%contacts_data((i-1)*stride_contacts+3) = contact_natoms
            transport_env%contacts_data((i-1)*stride_contacts+4) = contact_injsign

            IF (injecting_contact) THEN
               transport_env%contacts_data((i-1)*stride_contacts+5) = 1
            ELSE
               transport_env%contacts_data((i-1)*stride_contacts+5) = 0
            END IF
            transport_env%contacts_data((i-1)*stride_contacts+6) = contact_bandwidth*contact_natoms
         END DO
         transport_env%params%contacts_data = C_LOC(transport_env%contacts_data(1))
      ELSE
         IF (transport_env%params%method .EQ. transport_diagscalapack) THEN
            transport_env%params%num_contacts = 0
            transport_env%params%contacts_data = C_NULL_PTR
         ELSE
            CPABORT("No contact region is defined.")
         END IF
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
      TYPE(cp_dbcsr_type), INTENT(IN)                    :: template_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'transport_initialize', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, mp_group, mynode, numnodes, &
                                                            unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      CALL cite_reference(Bruck2014)

      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      mp_group = logger%para_env%group
      CALL mp_environ(numnodes, mynode, mp_group)

      ! Use all tasks if user input for TASKS_PER_POINT is 0 or greater than MPI size.
      IF ((transport_env%params%tasks_per_point .GT. numnodes) &
          .OR. (transport_env%params%tasks_per_point .EQ. 0)) THEN
         transport_env%params%tasks_per_point = numnodes
      ENDIF

      ! Determine valid TASKS_PER_POINT from user input
      DO WHILE (MOD(numnodes, transport_env%params%tasks_per_point) .NE. 0)
         transport_env%params%tasks_per_point = transport_env%params%tasks_per_point+1
      ENDDO

      IF (cp_dbcsr_has_symmetry(template_matrix)) THEN
         CALL cp_dbcsr_copy(transport_env%template_matrix_sym, template_matrix)
         CALL cp_dbcsr_desymmetrize(transport_env%template_matrix_sym, transport_env%template_matrix_nosym)
      ELSE
         CALL cp_dbcsr_copy(transport_env%template_matrix_nosym, template_matrix)
         CALL cp_dbcsr_copy(transport_env%template_matrix_sym, template_matrix)
      END IF

      CALL cp_dbcsr_create(transport_env%csr_sparsity, "CSR sparsity", &
                           template=transport_env%template_matrix_sym, &
                           data_type=dbcsr_type_real_8)
      CALL cp_dbcsr_copy(transport_env%csr_sparsity, transport_env%template_matrix_sym)

      CALL cp_dbcsr_to_csr_screening(ks_env, transport_env%csr_sparsity)

      IF (.NOT. transport_env%csr_screening) CALL cp_dbcsr_set(transport_env%csr_sparsity, 1.0)
      CALL cp_csr_create_from_dbcsr(transport_env%template_matrix_nosym, &
                                    transport_env%s_matrix, &
                                    csr_eqrow_ceil_dist, &
                                    csr_sparsity=transport_env%csr_sparsity, &
                                    numnodes=numnodes)

      CALL cp_csr_print_sparsity(transport_env%s_matrix, unit_nr)

      CALL cp_convert_dbcsr_to_csr(transport_env%template_matrix_nosym, transport_env%s_matrix)

      CALL cp_csr_create(transport_env%ks_matrix, transport_env%s_matrix)
      CALL cp_csr_create(transport_env%p_matrix, transport_env%s_matrix)

      CALL transport_eval_blocks(ks_env, transport_env)

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
!> \param[in]    extra_scf      whether or not an extra scf step will be performed
!> \par History
!>       12.2012 created [Hossein Bani-Hashemian]
!>       12.2014 revised [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE external_scf_method(transport_env, matrix_s, matrix_ks, matrix_p, &
                                  nelectron_spin, natoms, energy_diff, extra_scf)

      TYPE(transport_env_type), INTENT(INOUT)            :: transport_env
      TYPE(cp_dbcsr_type), INTENT(IN)                    :: matrix_s, matrix_ks
      TYPE(cp_dbcsr_type), INTENT(INOUT)                 :: matrix_p
      INTEGER, INTENT(IN)                                :: nelectron_spin, natoms
      REAL(dp), INTENT(IN)                               :: energy_diff
      LOGICAL, INTENT(IN)                                :: extra_scf

      CHARACTER(len=*), PARAMETER :: routineN = 'external_scf_method', &
         routineP = moduleN//':'//routineN

      TYPE(cp2k_csr_interop_type)                        :: ks_mat, p_mat, s_mat

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
      transport_env%params%extra_scf = extra_scf

      CALL csr_interop_nullify(s_mat)
      CALL csr_interop_nullify(ks_mat)
      CALL csr_interop_nullify(p_mat)

      CALL cp_dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_s)
      CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%s_matrix, s_mat)

      CALL cp_dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_ks)
      CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%ks_matrix, ks_mat)

      CALL cp_dbcsr_copy_into_existing(transport_env%template_matrix_sym, matrix_p)
      CALL convert_dbcsr_to_csr_interop(transport_env%template_matrix_sym, transport_env%p_matrix, p_mat)

      CALL c_method(transport_env%params, s_mat, ks_mat, p_mat)

      CALL convert_csr_interop_to_dbcsr(p_mat, transport_env%p_matrix, transport_env%template_matrix_nosym)
      CALL cp_dbcsr_copy(matrix_p, transport_env%template_matrix_nosym)

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

      TYPE(cp_dbcsr_type), INTENT(IN)          :: dbcsr_mat
      TYPE(csr_type), INTENT(INOUT)            :: csr_mat
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
      CALL cp_convert_dbcsr_to_csr(dbcsr_mat, csr_mat)

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

      ALLOCATE (nrows_local_all(0:num_pe-1), first_row_all(0:num_pe-1))
      CALL mp_allgather(csr_mat%nrows_local, nrows_local_all, mp_group)
      CALL cumsum_i(nrows_local_all, first_row_all)

      IF (mepos .EQ. 0) THEN
         csr_interop_mat%first_row = 0
      ELSE
         csr_interop_mat%first_row = first_row_all(mepos-1)
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
            cumsum(i) = cumsum(i-1)+arr(i)
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
      TYPE(csr_type), INTENT(INOUT)                      :: csr_mat
      TYPE(cp_dbcsr_type), INTENT(INOUT)                 :: dbcsr_mat

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
      CALL cp_convert_csr_to_dbcsr(dbcsr_mat, csr_mat)

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

      INTEGER                                            :: i, iat, ikind, natom, nkind, &
                                                            stride_contacts
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

      transport_env%params%nsgf = C_LOC(transport_env%nsgf(1))
      transport_env%params%zeff = C_LOC(transport_env%zeff(1))

      stride_contacts = 6
      DO i = 1, transport_env%params%num_contacts
         IF (transport_env%contacts_data((i-1)*stride_contacts+2) .LT. 0) THEN
            IF (transport_env%contacts_data((i-1)*stride_contacts+5) .EQ. 0) THEN
               transport_env%contacts_data((i-1)*stride_contacts+2) = natom/2; 
            ELSE IF (transport_env%contacts_data((i-1)*stride_contacts+4) .EQ. -1) THEN
               transport_env%contacts_data((i-1)*stride_contacts+2) = &
                  natom-transport_env%contacts_data((i-1)*stride_contacts+3)-transport_env%params%cutout(2)
            ELSE
               transport_env%contacts_data((i-1)*stride_contacts+2) = &
                  transport_env%params%cutout(1)
            END IF
         END IF
      END DO

   END SUBROUTINE transport_set_contact_params

! **************************************************************************************************
!> \brief determine the interaction range of the contact unit cells, the size of the self energy
!>        matrices and the adjacency matrix to partition the system into tridiagonal blocks
!> \param[inout] ks_env ks environment
!> \param[inout] transport_env transport environment
!> \author Sascha Brueck
! **************************************************************************************************
   SUBROUTINE transport_eval_blocks(ks_env, transport_env)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(transport_env_type), INTENT(INOUT)            :: transport_env

      CHARACTER(len=*), PARAMETER :: routineN = 'transport_eval_blocks', &
         routineP = moduleN//':'//routineN

      INTEGER :: colend, contact_bandwidth, contact_injsign, contact_natoms, contact_start, &
         handle, i, iatom, icol, icol_end, icol_start, ikind, irow, jatom, jkind, maxcol, mincol, &
         natom, nkind, num_glo_atom_pairs, num_loc_atom_pairs, num_tridiag_blocks, sigmasize, &
         stride_contacts, tridiag_end, tridiag_start
      INTEGER(KIND=int_1), ALLOCATABLE, DIMENSION(:, :)  :: fulladjmat
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, glo_atom_pairs, &
                                                            loc_atom_pairs, loc_atom_pairs_big, &
                                                            rcount, rdispl, tridiag_blocks_start
      LOGICAL                                            :: block_exists, do_symmetric
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list_a, basis_set_list_b
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: neighbour_list
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

! fulladjmat could probably be int_8 but is defined int_1 to save memory, just to be on the safe side

      NULLIFY (atomic_kind_set, qs_kind_set, nl_iterator)
      NULLIFY (basis_set_list_a, basis_set_list_b, basis_set_a, basis_set_b)

      CALL timeset(routineN, handle)

      CALL get_ks_env(ks_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      sab_orb=neighbour_list, &
                      natom=natom, &
                      para_env=para_env)

      ALLOCATE (loc_atom_pairs_big(2*natom*natom))
      num_loc_atom_pairs = 0

      nkind = SIZE(qs_kind_set)
      ALLOCATE (atom_of_kind(natom))
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind)
      CPASSERT(SIZE(neighbour_list) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=neighbour_list, symmetric=do_symmetric)
      CPASSERT(do_symmetric)
      ALLOCATE (basis_set_list_a(nkind), basis_set_list_b(nkind))
      CALL basis_set_list_setup(basis_set_list_a, "ORB", qs_kind_set)
      CALL basis_set_list_setup(basis_set_list_b, "ORB", qs_kind_set)

      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, neighbour_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, &
                                ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom)

         basis_set_a => basis_set_list_a(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list_b(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE

         loc_atom_pairs_big(num_loc_atom_pairs+1) = iatom
         loc_atom_pairs_big(num_loc_atom_pairs+2) = jatom
         num_loc_atom_pairs = num_loc_atom_pairs+2
      ENDDO

      ALLOCATE (loc_atom_pairs(num_loc_atom_pairs))
      DO i = 1, num_loc_atom_pairs
         loc_atom_pairs(i) = loc_atom_pairs_big(i)
      END DO
      DEALLOCATE (loc_atom_pairs_big)
      CALL neighbor_list_iterator_release(nl_iterator)
      DEALLOCATE (basis_set_list_a, basis_set_list_b)

      ALLOCATE (rcount(para_env%num_pe), rdispl(para_env%num_pe))
      CALL mp_allgather(num_loc_atom_pairs, rcount, para_env%group)
      num_glo_atom_pairs = SUM(rcount)
      rdispl(1) = 0
      DO i = 2, para_env%num_pe
         rdispl(i) = rdispl(i-1)+rcount(i-1)
      ENDDO

      ALLOCATE (glo_atom_pairs(num_glo_atom_pairs))
      CALL mp_allgather(loc_atom_pairs, glo_atom_pairs, rcount, rdispl, para_env%group)
      DEALLOCATE (loc_atom_pairs, rcount, rdispl)

      IF (para_env%mepos == para_env%source) THEN
         ALLOCATE (fulladjmat(natom, natom))
         fulladjmat = 0
         DO i = 1, num_glo_atom_pairs, 2
            irow = glo_atom_pairs(i)
            icol = glo_atom_pairs(i+1)
            fulladjmat(irow, icol) = 1
            fulladjmat(icol, irow) = 1
         END DO
         DO i = 1, transport_env%params%num_contacts
            stride_contacts = 6
            contact_bandwidth = transport_env%contacts_data((i-1)*stride_contacts+1)
            contact_start = transport_env%contacts_data((i-1)*stride_contacts+2)
            contact_natoms = transport_env%contacts_data((i-1)*stride_contacts+3)
            contact_injsign = transport_env%contacts_data((i-1)*stride_contacts+4)
            IF (contact_bandwidth .LE. 0) THEN
               contact_bandwidth = 0
               block_exists = .TRUE.
               DO WHILE (block_exists)
                  block_exists = .FALSE.
                  DO irow = contact_start+1, contact_start+contact_natoms
                     icol_start = contact_start+1+contact_injsign*(contact_bandwidth+1)*contact_natoms
                     icol_end = contact_start+contact_natoms+contact_injsign*(contact_bandwidth+1)*contact_natoms
                     DO icol = icol_start, icol_end
                        IF (fulladjmat(irow, icol) .EQ. 1) THEN
                           block_exists = .TRUE.
                           contact_bandwidth = contact_bandwidth+1
                           EXIT
                        END IF
                     END DO
                     IF (block_exists) EXIT
                  END DO
               END DO
               transport_env%contacts_data((i-1)*stride_contacts+1) = contact_bandwidth
            END IF
            tridiag_start = contact_start
            IF (contact_injsign .EQ. -1) tridiag_start = contact_start-(contact_bandwidth-1)*contact_natoms
            mincol = natom
            maxcol = 0
            DO irow = tridiag_start+1, tridiag_start+contact_bandwidth*contact_natoms
               icol_start = tridiag_start+1+contact_injsign*contact_bandwidth*contact_natoms
               icol_end = tridiag_start+contact_bandwidth*contact_natoms+contact_injsign*contact_bandwidth*contact_natoms
               DO icol = icol_start, icol_end
                  IF (fulladjmat(irow, icol) .EQ. 1) THEN
                     IF (icol < mincol) mincol = icol; 
                     IF (icol > maxcol) maxcol = icol; 
                  END IF
               END DO
            END DO
            sigmasize = maxcol-(tridiag_start+contact_bandwidth*contact_natoms)
            IF (contact_injsign .EQ. -1) sigmasize = tridiag_start+1-mincol
            transport_env%contacts_data((i-1)*stride_contacts+6) = sigmasize
         END DO
         ALLOCATE (tridiag_blocks_start(MAX(3, natom)))
         num_tridiag_blocks = 1
         tridiag_blocks_start(1) = transport_env%params%cutout(1)
         tridiag_blocks_start(2) = transport_env%contacts_data(2)+transport_env%contacts_data(6)
         tridiag_end = transport_env%contacts_data(stride_contacts+2)+ &
                       transport_env%contacts_data(stride_contacts+3)- &
                       transport_env%contacts_data(stride_contacts+6)
         maxcol = 0
         DO WHILE (tridiag_blocks_start(num_tridiag_blocks+1) < tridiag_end)
            DO irow = tridiag_blocks_start(num_tridiag_blocks)+1, tridiag_blocks_start(num_tridiag_blocks+1)
               colend = natom-transport_env%params%cutout(2)
               IF (irow .LE. transport_env%contacts_data(1)*transport_env%contacts_data(3)) THEN
                  colend = transport_env%contacts_data(2)+2*transport_env%contacts_data(1)*transport_env%contacts_data(3)
               END IF
               DO icol = irow, colend
                  IF (fulladjmat(irow, icol) .EQ. 1) THEN
                     IF (icol > maxcol) maxcol = icol; 
                  END IF
               END DO
            END DO
            num_tridiag_blocks = num_tridiag_blocks+1
            tridiag_blocks_start(num_tridiag_blocks+1) = maxcol
         END DO
         tridiag_blocks_start(1) = 0
         tridiag_blocks_start(num_tridiag_blocks+1) = natom
         ALLOCATE (transport_env%tridiag_blocks(num_tridiag_blocks))
         DO i = 1, num_tridiag_blocks
            transport_env%tridiag_blocks(i) = tridiag_blocks_start(i+1)-tridiag_blocks_start(i)
         END DO
         DEALLOCATE (tridiag_blocks_start)
         DEALLOCATE (fulladjmat)
      END IF
      DEALLOCATE (glo_atom_pairs)

      CALL mp_bcast(num_tridiag_blocks, para_env%source, para_env%group)
      transport_env%params%n_blocks = num_tridiag_blocks
      IF (para_env%mepos .NE. para_env%source) ALLOCATE (transport_env%tridiag_blocks(num_tridiag_blocks))
      CALL mp_bcast(transport_env%tridiag_blocks, para_env%source, para_env%group)
      transport_env%params%tridiag_blocks = C_LOC(transport_env%tridiag_blocks(1))

      DO i = 1, transport_env%params%num_contacts
         stride_contacts = 6
         CALL mp_bcast(transport_env%contacts_data((i-1)*stride_contacts+1), para_env%source, para_env%group)
         CALL mp_bcast(transport_env%contacts_data((i-1)*stride_contacts+6), para_env%source, para_env%group)
      END DO

      CALL timestop(handle)
   END SUBROUTINE transport_eval_blocks

END MODULE transport

