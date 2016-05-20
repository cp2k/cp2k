!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Determine active space Hamiltonian
!> \par History
!>      04.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE qs_active_space_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: &
        allocate_sto_basis_set, create_gto_from_sto_basis, deallocate_gto_basis_set, &
        deallocate_sto_basis_set, gto_basis_set_type, init_orb_basis_set, set_sto_basis_set, &
        srules, sto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_interface,              ONLY: cp_csr_create,&
                                              cp_dbcsr_allocate_matrix_set,&
                                              cp_dbcsr_copy,&
                                              cp_dbcsr_init,&
                                              cp_dbcsr_p_type,&
                                              cp_dbcsr_type,&
                                              csr_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_element, cp_fm_get_info, cp_fm_p_type, cp_fm_release, &
        cp_fm_set_all, cp_fm_set_element, cp_fm_to_fm, cp_fm_type
   USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: &
        cp_p_file, cp_print_key_finished_output, cp_print_key_should_output, cp_print_key_unit_nr, &
        debug_print_level, high_print_level, low_print_level, medium_print_level, &
        silent_print_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE input_constants,                 ONLY: &
        casci_canonical, eri_method_full_gpw, eri_method_gpw_ht, eri_operator_coulomb, &
        eri_operator_erf, eri_operator_erfc, eri_operator_gaussian, eri_operator_yukawa, &
        mao_projection, wannier_projection
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
   USE mathconstants,                   ONLY: fourpi
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_environ,&
                                              mp_sum,&
                                              mp_sync
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: ptable
   USE pw_env_methods,                  ONLY: pw_env_create,&
                                              pw_env_rebuild
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_release,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_integrate_function,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_rebuild,&
                                              pw_poisson_solve
   USE pw_poisson_types,                ONLY: ANALYTIC0D,&
                                              PERIODIC3D,&
                                              greens_fn_type,&
                                              pw_poisson_analytic,&
                                              pw_poisson_periodic,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_release,&
                                              pw_type
   USE qs_active_space_types,           ONLY: active_space_type,&
                                              create_active_space_type,&
                                              eri_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_loc_methods,                  ONLY: qs_loc_driver
   USE qs_loc_types,                    ONLY: qs_loc_env_create,&
                                              qs_loc_env_new_type,&
                                              qs_loc_env_release
   USE qs_loc_utils,                    ONLY: qs_loc_control_init,&
                                              qs_loc_env_init,&
                                              qs_loc_init
   USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                              calculate_subspace_eigenvalues
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              deallocate_task_list,&
                                              task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_active_space_methods'

   PUBLIC ::  active_space_main

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param input ...
!> \param logger ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE active_space_main(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'active_space_main', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=10)                                  :: cshell, lnam(5)
      INTEGER :: eri_method, eri_operator, eri_print, handle, i, iatom, ishell, isp, ispin, iw, j, &
         jm, m, mselect, n1, n2, nactive, nao, natom, ncol, nel, ninactive, ninactive_alpha, &
         ninactive_beta, nmo, nn1, nn2, nspins, ntot, zval
      INTEGER, DIMENSION(2)                              :: nactive_orb, nepol, ninactive_orb
      INTEGER, DIMENSION(5)                              :: nshell
      INTEGER, DIMENSION(:), POINTER                     :: invals
      LOGICAL                                            :: do_kpoints, ex_operator, ex_perd, &
                                                            explicit, isolated, store_wfn
      REAL(KIND=dp) :: eref, eri_eps_grid, eri_eps_int, eri_gpw_cutoff, eri_op_param, eri_rcut, &
         eri_rel_cutoff, fel, focc, maxocc, nze_percentage
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenvalues
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: ks_matrix, rho_ao
      TYPE(cp_dbcsr_type), POINTER                       :: denmat
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_localized
      TYPE(cp_fm_type), POINTER                          :: fm_ref, fm_target, fmat, mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(csr_type), POINTER                            :: eri_mat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: pro_basis_set
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env_occ, qs_loc_env_virt
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: as_input, loc_print, loc_section, &
                                                            print_orb

      !--------------------------------------------------------------------------------------------!

      as_input => section_vals_get_subs_vals(input, "DFT%PRINT%ACTIVE_SPACE")
      CALL section_vals_get(as_input, explicit=explicit)
      IF (.NOT. explicit) RETURN
      CALL timeset(routineN, handle)

      iw = cp_logger_get_default_io_unit(logger)

      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') &
            '!-----------------------------------------------------------------------------!'
         WRITE (iw, '(T26,A)') "Generate Active Space Hamiltonian"
         WRITE (iw, '(T27,A)') "Interface to CAS-CI and DMRG-CI"
         WRITE (iw, '(T2,A)') &
            '!-----------------------------------------------------------------------------!'
      END IF

      ! k-points?
      CALL get_qs_env(qs_env, do_kpoints=do_kpoints)
      IF (do_kpoints) THEN
         CALL cp_abort(__LOCATION__, "K-points not allowd for Active Space Interface")
      END IF

      NULLIFY (active_space_env)
      CALL create_active_space_type(active_space_env)

      ! input options

      ! isolated (molecular) system?
      CALL section_vals_val_get(as_input, "ISOLATED_SYSTEM", l_val=isolated)
      active_space_env%molecule = isolated
      IF (iw > 0) THEN
         IF (active_space_env%molecule) THEN
            WRITE (iw, '(T3,A)') "System is treated without periodicity"
         END IF
      END IF

      CALL section_vals_val_get(as_input, "ACTIVE_ELECTRONS", i_val=nactive)
      IF (nactive <= 0) CPABORT("Specify a positive number of active electrons.")
      CALL get_qs_env(qs_env, nelectron_total=ntot, nelectron_spin=nepol)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins
      CALL section_vals_val_get(as_input, "INACTIVE_ELECTRONS", explicit=explicit)
      CALL section_vals_val_get(as_input, "INACTIVE_ELECTRONS", i_vals=invals)
      IF (nspins > 1) THEN
         IF (.NOT. explicit) THEN
            CALL cp_abort(__LOCATION__, "Number of Inactive Electrons has to be specified"// &
                          " in spin polarised case.")
         END IF
         ninactive_alpha = invals(1)
         ninactive_beta = invals(2)
         ninactive = ninactive_alpha+ninactive_beta
      ELSE
         IF (explicit) THEN
            ninactive = invals(1)
         ELSE
            ninactive = ntot-nactive
         END IF
         ninactive_alpha = ninactive
         ninactive_beta = 0
      END IF
      CPASSERT(ninactive >= 0)
      CPASSERT(ninactive_alpha >= 0)
      CPASSERT(ninactive_beta >= 0)
      CPASSERT(ntot == ninactive+nactive)
      IF (nspins > 1) THEN
         CPASSERT(nepol(1) >= ninactive_alpha)
         CPASSERT(nepol(2) >= ninactive_beta)
      END IF

      IF (iw > 0) THEN
         IF (nspins < 2) THEN
            WRITE (iw, '(T4,A,T69,I10)') "Total number of electrons", ntot
            WRITE (iw, '(T4,A,T69,I10)') "Number of active electrons", nactive
            WRITE (iw, '(T4,A,T69,I10)') "Number of inactive electrons", ninactive
         ELSE
            WRITE (iw, '(T4,A,T69,I10)') "Total number of electrons", ntot
            WRITE (iw, '(T4,A,T69,I10)') "Number of active electrons", nactive
            WRITE (iw, '(T4,A,T69,I10)') "Number of inactive electrons", ninactive
            WRITE (iw, '(T4,A,T69,I10)') "Number of inactive electrons (alpha)", ninactive_alpha
            WRITE (iw, '(T4,A,T69,I10)') "Number of inactive electrons (beta)", ninactive_beta
         END IF
      END IF
      active_space_env%nactive = nactive
      active_space_env%ninactive = ninactive
      IF (nspins > 1) THEN
         active_space_env%nelectrons(1) = ntot
         active_space_env%nelectrons(2) = 0
      ELSE
         active_space_env%nelectrons(1) = nepol(1)
         active_space_env%nelectrons(2) = nepol(2)
      END IF
      active_space_env%ninspin(1) = ninactive_alpha
      active_space_env%ninspin(2) = ninactive_beta
      active_space_env%multiplicity = dft_control%multiplicity
      active_space_env%nspins = nspins

      ! define the active/inactive space orbitals
      CALL section_vals_val_get(as_input, "ACTIVE_ORBITALS", explicit=explicit)
      CALL section_vals_val_get(as_input, "ACTIVE_ORBITALS", i_vals=invals)
      IF (.NOT. explicit) THEN
         CALL cp_abort(__LOCATION__, "Number of Active Orbitals has to be specified.")
      END IF
      CPASSERT(SIZE(invals) >= nspins)
      DO ispin = 1, nspins
         nactive_orb(ispin) = invals(ispin)
      END DO
      ninactive_orb = 0
      CALL section_vals_val_get(as_input, "INACTIVE_ORBITALS", explicit=explicit)
      CALL section_vals_val_get(as_input, "INACTIVE_ORBITALS", i_vals=invals)
      IF (.NOT. explicit) THEN
         IF (nspins == 1) THEN
            ninactive_orb(1) = active_space_env%ninspin(1)/2
         ELSE
            DO ispin = 1, nspins
               ninactive_orb(ispin) = active_space_env%ninspin(ispin)
            END DO
         END IF
      ELSE
         CPASSERT(SIZE(invals) >= nspins)
         DO ispin = 1, nspins
            ninactive_orb(ispin) = invals(ispin)
         END DO
      END IF
      IF (nspins == 1) THEN
         CPASSERT(MOD(ninactive, 2) == 0)
         CPASSERT(ninactive_orb(1) >= ninactive/2)
      ELSE
         CPASSERT(ninactive_orb(1) >= ninactive_alpha)
         CPASSERT(ninactive_orb(2) >= ninactive_beta)
      END IF

      CALL get_qs_env(qs_env, mos=mos)
      IF (.NOT. ASSOCIATED(mos(1)%mo_set)) THEN
         CALL cp_abort(__LOCATION__, "A reference MO set has to be available when determining"// &
                       " the active space orbitals.")
      END IF
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=fm_ref, nao=nao)
      maxocc = 2.0_dp
      IF (nspins > 1) maxocc = 1.0_dp
      ALLOCATE (active_space_env%mos_active(nspins))
      ALLOCATE (active_space_env%mos_inactive(nspins))
      DO ispin = 1, nspins
         nmo = nactive_orb(ispin)
         NULLIFY (mo_set)
         CALL allocate_mo_set(mo_set, nao, nmo, 0, 0.0_dp, maxocc, 0.0_dp)
         CALL init_mo_set(mo_set, fm_ref=fm_ref, name="Active Space MO")
         active_space_env%mos_active(ispin)%mo_set => mo_set
         mo_set => active_space_env%mos_inactive(ispin)%mo_set
         nmo = ninactive_orb(ispin)
         nel = active_space_env%ninspin(ispin)
         NULLIFY (mo_set)
         CALL allocate_mo_set(mo_set, nao, nmo, nel, REAL(nel, KIND=dp), maxocc, 0.0_dp)
         CALL init_mo_set(mo_set, fm_ref=fm_ref, name="Inactive Space MO")
         active_space_env%mos_inactive(ispin)%mo_set => mo_set
      END DO

      CALL section_vals_val_get(as_input, "ORBITAL_SELECTION", i_val=mselect)
      IF (iw > 0) THEN
         SELECT CASE (mselect)
         CASE DEFAULT
            CPABORT("Unknown orbital selection method")
         CASE (casci_canonical)
            WRITE (iw, '(/,T4,A)') &
               "Active space orbitals selected using energy ordered canonical orbitals"
         CASE (wannier_projection)
            WRITE (iw, '(/,T4,A)') &
               "Active space orbitals selected using projected Wannier orbitals"
         CASE (mao_projection)
            WRITE (iw, '(/,T4,A)') &
               "Active space orbitals selected using modified atomic orbitals (MAO)"
         END SELECT
         IF (nspins < 2) THEN
            WRITE (iw, '(T4,A,T69,I10)') "Number of active orbitals", nactive_orb(1)
            WRITE (iw, '(T4,A,T69,I10)') "Number of inactive orbitals", ninactive_orb(1)
         ELSE
            WRITE (iw, '(T4,A,T69,I10)') "Number of active orbitals", SUM(nactive_orb)
            WRITE (iw, '(T4,A,T69,I10)') "Number of active orbitals (alpha)", nactive_orb(1)
            WRITE (iw, '(T4,A,T69,I10)') "Number of active orbitals (beta)", nactive_orb(2)
            WRITE (iw, '(T4,A,T69,I10)') "Number of inactive orbitals", SUM(ninactive_orb)
            WRITE (iw, '(T4,A,T69,I10)') "Number of inactive orbitals (alpha)", ninactive_orb(1)
            WRITE (iw, '(T4,A,T69,I10)') "Number of inactive orbitals (beta)", ninactive_orb(2)
         END IF
      END IF

      ! get projection spaces
      CALL section_vals_val_get(as_input, "SUBSPACE_ATOM", i_val=iatom, explicit=explicit)
      IF (explicit) THEN
         CALL get_qs_env(qs_env, natom=natom)
         IF (iatom <= 0 .OR. iatom > natom) THEN
            IF (iw > 0) THEN
               WRITE (iw, '(/,T4,A,I3)') "ERROR: SUBSPACE_ATOM number is not valid", iatom
            END IF
            CPABORT("Select a valid SUBSPACE_ATOM")
         END IF
      END IF
      CALL section_vals_val_get(as_input, "SUBSPACE_SHELL", c_val=cshell, explicit=explicit)
      nshell = 0
      lnam = ""
      IF (explicit) THEN
         cshell = ADJUSTL(cshell)
         n1 = 1
         DO i = 1, 5
            ishell = i
            IF (cshell(n1:n1) == " ") THEN
               ishell = ishell-1
               EXIT
            END IF
            READ (cshell(n1:), "(I1,A1)") nshell(i), lnam(i)
            n1 = n1+2
         END DO
      END IF

      ! generate orbitals
      SELECT CASE (mselect)
      CASE DEFAULT
         CPABORT("Unknown orbital selection method")
      CASE (casci_canonical)
         ! create canonical orbitals
         IF (dft_control%restricted) THEN
            CPABORT(" Unclear how we define MOs in the restricted case ... stopping")
         ELSE
            IF (dft_control%do_admm) THEN
               CPABORT("ADMM currently not possible for canonical orbital options")
            END IF
            nmo = 0
            DO ispin = 1, nspins
               CALL get_mo_set(mos(ispin)%mo_set, nmo=m)
               IF (m < ninactive_orb(ispin)+nactive_orb(ispin)) THEN
                  CPABORT("Not enough canonical orbitals available.")
               END IF
               nmo = MAX(m, nmo)
            END DO
            ALLOCATE (eigenvalues(nmo, nspins))
            eigenvalues = 100000._dp
            CALL get_qs_env(qs_env, matrix_ks=ks_matrix)
            IF (mos(1)%mo_set%use_mo_coeff_b) THEN
               CPABORT("Needs MOs in FM format, missing code for DBCSR")
            ELSE
               DO ispin = 1, nspins
                  CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=fm_ref, nmo=nmo)
                  CALL calculate_subspace_eigenvalues(fm_ref, ks_matrix(ispin)%matrix, &
                                                      evals_arg=eigenvalues(:, ispin), do_rotation=.TRUE.)
                  ! copy inactive orbitals
                  mo_set => active_space_env%mos_inactive(ispin)%mo_set
                  CALL get_mo_set(mo_set, mo_coeff=fm_target)
                  ncol = ninactive_orb(ispin)
                  CALL cp_fm_to_fm(fm_ref, fm_target, ncol, 1, 1)
                  ! copy active orbitals
                  mo_set => active_space_env%mos_active(ispin)%mo_set
                  CALL get_mo_set(mo_set, mo_coeff=fm_target)
                  ncol = nactive_orb(ispin)
                  CALL cp_fm_to_fm(fm_ref, fm_target, ncol, ninactive_orb(ispin)+1, 1)
               END DO
            END IF
            IF (iw > 0) THEN
               DO ispin = 1, nspins
                  WRITE (iw, '(/,T4,A,I3,T65,A)') "Canonical Orbital Selection for spin", ispin, &
                     "[atomic units]"
                  DO i = 1, ninactive_orb(ispin), 4
                     jm = MIN(3, ninactive_orb(ispin)-i)
                     WRITE (iw, '(T3,4(F14.6,A5))') (eigenvalues(i+j, 1), " [I]", j=0, jm)
                  END DO
                  DO i = ninactive_orb(ispin)+1, ninactive_orb(ispin)+nactive_orb(ispin), 4
                     jm = MIN(3, ninactive_orb(ispin)+nactive_orb(ispin)-i)
                     WRITE (iw, '(T3,4(F14.6,A5))') (eigenvalues(i+j, 1), " [A]", j=0, jm)
                  END DO
               END DO
            END IF
            DEALLOCATE (eigenvalues)
         END IF

      CASE (wannier_projection)
         NULLIFY (loc_section, loc_print)
         loc_section => section_vals_get_subs_vals(as_input, "LOCALIZE")
         CPASSERT(ASSOCIATED(loc_section))
         loc_print => section_vals_get_subs_vals(as_input, "LOCALIZE%PRINT")
         NULLIFY (qs_loc_env_occ, qs_loc_env_virt)
         CALL qs_loc_env_create(qs_loc_env_occ)
         ALLOCATE (mos_localized(nspins))
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
            CALL cp_fm_create(mos_localized(ispin)%matrix, mo_coeff%matrix_struct)
            CALL cp_fm_to_fm(mo_coeff, mos_localized(ispin)%matrix)
         END DO
         ! occupied states
         NULLIFY (qs_loc_env_occ)
         CALL qs_loc_env_create(qs_loc_env_occ)
         CALL qs_loc_control_init(qs_loc_env_occ, loc_section, do_homo=.TRUE.)
         CALL qs_loc_init(qs_env, qs_loc_env_occ, loc_section, &
                          mos_localized=mos_localized, do_homo=.TRUE.)
         DO ispin = 1, nspins
            CALL qs_loc_driver(qs_env, qs_loc_env_occ, loc_print, myspin=ispin, &
                               ext_mo_coeff=mos_localized(ispin)%matrix)
         END DO
         ! virtual states
         NULLIFY (qs_loc_env_virt)
         CALL qs_loc_env_create(qs_loc_env_virt)
         CALL qs_loc_control_init(qs_loc_env_virt, loc_section, do_homo=.FALSE.)
         CALL qs_loc_init(qs_env, qs_loc_env_virt, loc_section, &
                          mos_localized=mos_localized, do_homo=.FALSE.)
         CALL qs_loc_env_init(qs_loc_env_virt, qs_loc_env_virt%localized_wfn_control, qs_env, &
                              loc_coeff=mos_localized)
         DO ispin = 1, nspins
            CALL qs_loc_driver(qs_env, qs_loc_env_virt, loc_print, myspin=ispin, &
                               ext_mo_coeff=mos_localized(ispin)%matrix)
         END DO
         !
         ! get definition of subspace
         ! iatom, ishell, nshell, lnam
         CPASSERT(iatom > 0)
         CPASSERT(ishell > 0)
         CALL get_qs_env(qs_env, particle_set=particle_set)
         atomic_kind => particle_set(iatom)%atomic_kind
         CALL get_atomic_kind(atomic_kind=atomic_kind, z=zval)
         NULLIFY (pro_basis_set)
         CALL create_pro_basis(pro_basis_set, zval, ishell, nshell, lnam)
         !
         CALL deallocate_gto_basis_set(pro_basis_set)
         !
         CPABORT("not yet available")
         !
         CALL qs_loc_env_release(qs_loc_env_occ)
         CALL qs_loc_env_release(qs_loc_env_virt)
         CALL cp_fm_vect_dealloc(mos_localized)
         !
         CPABORT("not yet available")
         !
      CASE (mao_projection)
         !
         CPABORT("not yet available")
         !
      END SELECT

      ! Print orbitals on Cube files
      print_orb => section_vals_get_subs_vals(as_input, "PRINT_ORBITAL_CUBES")
      CALL section_vals_get(print_orb, explicit=explicit)
      IF (explicit) THEN
         !
         CALL print_orbital_cubes(print_orb, qs_env, active_space_env%mos_active)
         !
      END IF

      ! calculate inactive density matrix
      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CPASSERT(ASSOCIATED(rho_ao))
      CALL cp_dbcsr_allocate_matrix_set(active_space_env%pmat_inactive, nspins)
      DO ispin = 1, nspins
         ALLOCATE (denmat)
         CALL cp_dbcsr_init(denmat)
         CALL cp_dbcsr_copy(denmat, rho_ao(ispin)%matrix)
         mo_set => active_space_env%mos_inactive(ispin)%mo_set
         CALL calculate_density_matrix(mo_set, denmat)
         active_space_env%pmat_inactive(ispin)%matrix => denmat
      END DO

      ! generate integrals
      ! make sure that defaults are set correctly (from basic calculation)
      ! make sure that periodicity is consistent
      CALL section_vals_val_get(as_input, "ERI%METHOD", i_val=eri_method)
      active_space_env%eri%method = eri_method
      CALL section_vals_val_get(as_input, "ERI%OPERATOR", i_val=eri_operator, explicit=ex_operator)
      active_space_env%eri%operator = eri_operator
      CALL section_vals_val_get(as_input, "ERI%OPERATOR_PARAMETER", r_val=eri_op_param)
      active_space_env%eri%operator_parameter = eri_op_param
      CALL section_vals_val_get(as_input, "ERI%CUTOFF_RADIUS", r_val=eri_rcut)
      active_space_env%eri%cutoff_radius = eri_rcut
      CALL section_vals_val_get(as_input, "ERI%PERIODICITY", i_vals=invals, explicit=ex_perd)
      CALL section_vals_val_get(as_input, "ERI%EPS_INTEGRAL", r_val=eri_eps_int)
      active_space_env%eri%eps_integral = eri_eps_int
      IF (active_space_env%molecule) THEN
         ! check that we are in a non-periodic setting
         CALL get_qs_env(qs_env, cell=cell)
         IF (SUM(cell%perd) /= 0) THEN
            CPABORT("Active space option ISOLATED_SYSTEM requires non-periodic setting")
         END IF
         IF (ex_perd) THEN
            IF (SUM(invals) /= 0) THEN
               CPABORT("Active space option ISOLATED_SYSTEM requires non-periodic setting")
            END IF
         END IF
         active_space_env%eri%periodicity(1:3) = 0
      ELSE IF (ex_perd) THEN
         IF (SIZE(invals) == 1) THEN
            active_space_env%eri%periodicity(1:3) = invals(1)
         ELSE
            active_space_env%eri%periodicity(1:3) = invals(1:3)
         END IF
      ELSE
         CALL get_qs_env(qs_env, cell=cell)
         active_space_env%eri%periodicity(1:3) = cell%perd(1:3)
      END IF
      IF (iw > 0) THEN
         WRITE (iw, '(/,T4,A)') "Calculation or Electron Repulsion Intergrals"
         SELECT CASE (eri_method)
         CASE (eri_method_full_gpw)
            WRITE (iw, '(T4,A,T48,A)') "Integration method", "GPW Fourier transform over MOs"
         CASE (eri_method_gpw_ht)
            WRITE (iw, '(T4,A,T43,A)') "Integration method", "Half transformed integrals from GPW"
         CASE DEFAULT
            CPABORT("Unknown ERI method")
         END SELECT
         SELECT CASE (eri_operator)
         CASE (eri_operator_coulomb)
            WRITE (iw, '(T4,A,T40,A,T73,A)') "ERI operator", "Coulomb", "<1/R>"
         CASE (eri_operator_yukawa)
            WRITE (iw, '(T4,A,T40,A,T67,A)') "ERI operator", "Yukawa", "<EXP(-a*R)/R>"
            WRITE (iw, '(T4,A,T63,F14.6)') "ERI operator parameter", eri_op_param
         CASE (eri_operator_erf)
            WRITE (iw, '(T4,A,T40,A,T67,A)') "ERI operator", "Error function", "<ERF(a*R)/R>"
            WRITE (iw, '(T4,A,T63,F14.6)') "ERI operator parameter", eri_op_param
         CASE (eri_operator_erfc)
            WRITE (iw, '(T4,A,T40,A,T66,A)') "ERI operator", "Compl. error function", "<ERFC(a*R)/R>"
            WRITE (iw, '(T4,A,T63,F14.6)') "ERI operator parameter", eri_op_param
         CASE (eri_operator_gaussian)
            WRITE (iw, '(T4,A,T40,A,T66,A)') "ERI operator", "Gaussian attenuated", "<EXP(-a*R^2)/R>"
            WRITE (iw, '(T4,A,T63,F14.6)') "ERI operator parameter", eri_op_param
         CASE DEFAULT
            CPABORT("Unknown ERI operator")
         END SELECT
         WRITE (iw, '(T4,A,T66,E12.4)') "Accuracy of ERI", eri_eps_int
         WRITE (iw, '(T4,A,T69,3I3)') "Periodicity", active_space_env%eri%periodicity(1:3)
         IF (PRODUCT(active_space_env%eri%periodicity(1:3)) == 0) THEN
            IF (eri_rcut > 0.0_dp) WRITE (iw, '(T4,A,T63,F14.6)') "Periodicity (Cutoff)", eri_rcut
         END IF
         IF (nspins < 2) THEN
            WRITE (iw, '(T4,A,T66,I12)') "Total Number of ERI", (nactive_orb(1)**4)/8
         ELSE
            WRITE (iw, '(T4,A,T66,I12)') "Total Number of ERI (aa|aa)", (nactive_orb(1)**4)/8
            WRITE (iw, '(T4,A,T66,I12)') "Total Number of ERI (bb|bb)", (nactive_orb(2)**4)/8
            WRITE (iw, '(T4,A,T66,I12)') "Total Number of ERI (aa|bb)", &
               (nactive_orb(1)**2)*(nactive_orb(2)**2)/4
         END IF
      END IF

      ! allocate container for integrals (CSR matrix)
      CALL get_qs_env(qs_env, para_env=para_env)
      m = (nspins*(nspins+1))/2
      ALLOCATE (active_space_env%eri%eri(m))
      DO i = 1, m
         ALLOCATE (active_space_env%eri%eri(i)%csr_mat)
         eri_mat => active_space_env%eri%eri(i)%csr_mat
         IF (i == 1) THEN
            n1 = nactive_orb(1)
            n2 = nactive_orb(1)
         ELSEIF (i == 2) THEN
            n1 = nactive_orb(1)
            n2 = nactive_orb(2)
         ELSE
            n1 = nactive_orb(2)
            n2 = nactive_orb(2)
         END IF
         nn1 = (n1*(n1+1))/2
         nn2 = (n2*(n2+1))/2
         CALL cp_csr_create(eri_mat, nn1, nn2, 0_int_8, 0, 0, para_env%group)
      END DO

      SELECT CASE (eri_method)
      CASE (eri_method_full_gpw, eri_method_gpw_ht)
         CALL section_vals_val_get(as_input, "ERI_GPW%EPS_GRID", r_val=eri_eps_grid)
         active_space_env%eri%eri_gpw%eps_grid = eri_eps_grid
         CALL section_vals_val_get(as_input, "ERI_GPW%CUTOFF", r_val=eri_gpw_cutoff)
         active_space_env%eri%eri_gpw%cutoff = eri_gpw_cutoff
         CALL section_vals_val_get(as_input, "ERI_GPW%REL_CUTOFF", r_val=eri_rel_cutoff)
         active_space_env%eri%eri_gpw%rel_cutoff = eri_rel_cutoff
         CALL section_vals_val_get(as_input, "ERI_GPW%PRINT_LEVEL", i_val=eri_print)
         active_space_env%eri%eri_gpw%print_level = eri_print
         CALL section_vals_val_get(as_input, "ERI_GPW%STORE_WFN", l_val=store_wfn)
         active_space_env%eri%eri_gpw%store_wfn = store_wfn
         active_space_env%eri%eri_gpw%redo_poisson = (ex_operator .OR. ex_perd)
         IF (iw > 0) THEN
            WRITE (iw, '(T4,A,T68,F10.4)') "ERI_GPW| Energy cutoff [Ry]", eri_gpw_cutoff
            WRITE (iw, '(T4,A,T68,F10.4)') "ERI_GPW| Relative energy cutoff [Ry]", eri_rel_cutoff
         END IF
         !
         CALL calculate_eri_gpw(active_space_env%mos_active, active_space_env%eri, qs_env, iw)
         !
      CASE DEFAULT
         CPABORT("Unknown ERI method")
      END SELECT
      IF (iw > 0) THEN
         DO isp = 1, SIZE(active_space_env%eri%eri)
            eri_mat => active_space_env%eri%eri(isp)%csr_mat
            nze_percentage = 100.0_dp*(REAL(eri_mat%nze_total, KIND=dp) &
                                       /REAL(eri_mat%nrows_total, KIND=dp))/REAL(eri_mat%ncols_total, KIND=dp)
            WRITE (iw, '(T4,A,I2,T30,A,T66,I12)') "ERI_GPW| Spinmatrix:", isp, &
               "Number of  CSR non-zero elements:", eri_mat%nze_total
            WRITE (iw, '(T4,A,I2,T30,A,T66,F12.4)') "ERI_GPW| Spinmatrix:", isp, &
               "Percentage CSR non-zero elements:", nze_percentage
         END DO
      END IF

      ! set the reference active space density matrix
      nspins = active_space_env%nspins
      ALLOCATE (active_space_env%p_ref(nspins))
      DO isp = 1, nspins
         mo_set => active_space_env%mos_active(isp)%mo_set
         CALL get_mo_set(mo_set, mo_coeff=mo_coeff)
         NULLIFY (active_space_env%p_ref(isp)%matrix)
         CALL create_subspace_matrix(mo_coeff, active_space_env%p_ref(isp)%matrix)
      END DO
      SELECT CASE (mselect)
      CASE DEFAULT
         CPABORT("Unknown orbital selection method")
      CASE (casci_canonical)
         focc = 2.0_dp
         IF (nspins == 2) focc = 1.0_dp
         DO isp = 1, nspins
            fmat => active_space_env%p_ref(isp)%matrix
            CALL cp_fm_set_all(fmat, alpha=0.0_dp)
            n1 = active_space_env%nelectrons(isp)
            DO i = 1, nactive_orb(isp)
               fel = MIN(focc, REAL(n1, KIND=dp))
               CALL cp_fm_set_element(fmat, i, i, fel)
               n1 = n1-NINT(fel)
               n1 = MAX(n1, 0)
            END DO
         END DO
      CASE (wannier_projection)
         CPABORT("NOT IMPLEMENTED")
      CASE (mao_projection)
         CPABORT("NOT IMPLEMENTED")
      END SELECT

      ! calculate one-electron operators in the subspace
      CALL calculate_operators(active_space_env%mos_active, qs_env, active_space_env)

      ! Output a FCIDUMP file
      CALL get_qs_env(qs_env, energy=energy)
      eref = energy%total
      CALL fcidump(active_space_env, as_input, eref)

      CALL set_qs_env(qs_env, active_space=active_space_env)

      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') &
            '!-------------------- End of Active Space Interface --------------------------!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE active_space_main

! **************************************************************************************************
!> \brief computes the one-electron operators in the subspace of the provided orbital set
!> \param mos ...
!> \param qs_env ...
!> \param active_space_env ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE calculate_operators(mos, qs_env, active_space_env)

      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(active_space_type), POINTER                   :: active_space_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_operators', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, is, nmo, nspins
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: h_matrix, ks_matrix
      TYPE(cp_fm_type), POINTER                          :: hmat, ksmat, mo_coeff

      CALL timeset(routineN, handle)

      ! Kohn-Sham / Fock operator
      CALL get_qs_env(qs_env, matrix_ks_kp=ks_matrix)
      IF (SIZE(ks_matrix, 2) > 1) THEN
         CPABORT("No k-points allowed at this point")
      END IF
      nspins = active_space_env%nspins
      ALLOCATE (active_space_env%ks_sub(nspins))
      DO is = 1, nspins
         CALL get_mo_set(mo_set=mos(is)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
         NULLIFY (ksmat)
         CALL subspace_operator(mo_coeff, ks_matrix(is, 1)%matrix, ksmat)
         active_space_env%ks_sub(is)%matrix => ksmat
      END DO

      ! Core Hamiltonian
      NULLIFY (h_matrix)
      CALL get_qs_env(qs_env=qs_env, matrix_h_kp=h_matrix)
      IF (SIZE(h_matrix, 2) > 1) THEN
         CPABORT("No k-points allowed at this point")
      END IF
      ALLOCATE (active_space_env%h_sub(nspins))
      DO is = 1, nspins
         CALL get_mo_set(mo_set=mos(is)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
         NULLIFY (hmat)
         CALL subspace_operator(mo_coeff, h_matrix(1, 1)%matrix, hmat)
         active_space_env%h_sub(is)%matrix => hmat
      END DO

      CALL timestop(handle)

   END SUBROUTINE calculate_operators

! **************************************************************************************************
!> \brief computes a one-electron operator in the subspace of the provided orbital set
!> \param orbitals ...
!> \param op_matrix operator matrix in AO basis
!> \param op_sub operator in orbital basis
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE subspace_operator(orbitals, op_matrix, op_sub)

      TYPE(cp_fm_type), POINTER                          :: orbitals
      TYPE(cp_dbcsr_type), POINTER                       :: op_matrix
      TYPE(cp_fm_type), POINTER                          :: op_sub

      CHARACTER(len=*), PARAMETER :: routineN = 'subspace_operator', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ncol, nrow
      TYPE(cp_fm_type), POINTER                          :: vectors

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(op_sub))

      CALL cp_fm_get_info(matrix=orbitals, ncol_global=ncol, nrow_global=nrow)

      IF (ncol .GT. 0) THEN

         CALL cp_fm_create(vectors, orbitals%matrix_struct, "vectors")

         CALL create_subspace_matrix(orbitals, op_sub)

         CALL cp_dbcsr_sm_fm_multiply(op_matrix, orbitals, vectors, ncol)

         CALL cp_gemm('T', 'N', ncol, ncol, nrow, 1.0_dp, orbitals, vectors, 0.0_dp, op_sub)

         CALL cp_fm_release(vectors)

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE subspace_operator

! **************************************************************************************************
!> \brief creates a matrix of subspace size
!> \param orbitals ...
!> \param op_sub operator in orbital basis
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE create_subspace_matrix(orbitals, op_sub)

      TYPE(cp_fm_type), POINTER                          :: orbitals, op_sub

      CHARACTER(len=*), PARAMETER :: routineN = 'create_subspace_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ncol, nrow
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CPASSERT(.NOT. ASSOCIATED(op_sub))

      CALL cp_fm_get_info(matrix=orbitals, ncol_global=ncol, nrow_global=nrow)

      IF (ncol .GT. 0) THEN

         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, nrow_global=ncol, ncol_global=ncol, &
                                  para_env=orbitals%matrix_struct%para_env, &
                                  context=orbitals%matrix_struct%context)
         CALL cp_fm_create(op_sub, fm_struct, name="Subspace operator")
         CALL cp_fm_struct_release(fm_struct)

      ENDIF

   END SUBROUTINE create_subspace_matrix

! **************************************************************************************************
!> \brief computes a electron repulsion integrals using GPW technology
!> \param mos ...
!> \param eri_env ...
!> \param qs_env ...
!> \param iw ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE calculate_eri_gpw(mos, eri_env, qs_env, iw)

      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(eri_type)                                     :: eri_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_eri_gpw', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i_multigrid, icount2, intcount, irange(2), isp, isp1, isp2, ispin, iwa1, &
         iwa12, iwa2, iwb1, iwb2, iwbs, iwbt, iwfn, n_multigrid, nmm, nmo, nmo1, nmo2, nspins, nx
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: eri_index
      LOGICAL                                            :: print1, print2, &
                                                            skip_load_balance_distributed
      REAL(KIND=dp)                                      :: cutoff_old, dvol, erint, pair_int, &
                                                            progression_factor, rc, &
                                                            relative_cutoff_old, rsize
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: e_cutoff_old, eri
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_p_type)                                    :: pot_g, rho_g, rho_r, wfn_r
      TYPE(pw_p_type), DIMENSION(:, :), POINTER          :: wfn_a
      TYPE(pw_p_type), POINTER                           :: wfn1, wfn2, wfn3, wfn4
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(task_list_type), POINTER                      :: task_list_sub

      CALL timeset(routineN, handle)

      ! print levels
      SELECT CASE (eri_env%eri_gpw%print_level)
      CASE (silent_print_level)
         print1 = .FALSE.
         print2 = .FALSE.
      CASE (low_print_level)
         print1 = .FALSE.
         print2 = .FALSE.
      CASE (medium_print_level)
         print1 = .TRUE.
         print2 = .FALSE.
      CASE (high_print_level)
         print1 = .TRUE.
         print2 = .TRUE.
      CASE (debug_print_level)
         print1 = .TRUE.
         print2 = .TRUE.
      CASE DEFAULT
         ! do nothing
      END SELECT

      ! This should be done differently! Copied from MP2 code
      CALL get_qs_env(qs_env, dft_control=dft_control)
      progression_factor = dft_control%qs_control%progression_factor
      n_multigrid = SIZE(dft_control%qs_control%e_cutoff)
      ALLOCATE (e_cutoff_old(n_multigrid))
      e_cutoff_old(:) = dft_control%qs_control%e_cutoff
      cutoff_old = dft_control%qs_control%cutoff

      dft_control%qs_control%cutoff = eri_env%eri_gpw%cutoff*0.5_dp
      dft_control%qs_control%e_cutoff(1) = dft_control%qs_control%cutoff
      DO i_multigrid = 2, n_multigrid
         dft_control%qs_control%e_cutoff(i_multigrid) = dft_control%qs_control%e_cutoff(i_multigrid-1) &
                                                        /progression_factor
      END DO

      relative_cutoff_old = dft_control%qs_control%relative_cutoff
      dft_control%qs_control%relative_cutoff = eri_env%eri_gpw%rel_cutoff*0.5_dp

      ! Generate the appropriate  pw_env
      NULLIFY (pw_env_sub)
      CALL pw_env_create(pw_env_sub)
      CALL pw_env_rebuild(pw_env_sub, qs_env)
      CALL pw_env_get(pw_env_sub, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      IF (eri_env%eri_gpw%redo_poisson) THEN
         IF (.NOT. (poisson_env%parameters%solver == pw_poisson_analytic .OR. &
                    poisson_env%parameters%solver == pw_poisson_periodic)) THEN
            CPABORT("Only use simple analytic Poisson solvers with ERI calculation")
         END IF
         poisson_env%parameters%periodic = eri_env%periodicity
         CALL pw_poisson_rebuild(poisson_env)
         IF (eri_env%cutoff_radius > 0.0_dp) THEN
            poisson_env%green_fft%radius = eri_env%cutoff_radius
         ELSE
            CALL get_qs_env(qs_env, cell=cell)
            rc = cell%hmat(1, 1)
            DO iwa1 = 1, 3
               rc = MIN(rc, 0.5_dp*cell%hmat(iwa1, iwa1))
            END DO
            poisson_env%green_fft%radius = rc
         END IF
         CALL pw_eri_green_create(poisson_env%green_fft, eri_env)
         IF (print1 .AND. iw > 0) THEN
            WRITE (iw, "(T4,'ERI_GPW|',' Redefine Poisson Greens function ')")
         END IF
      END IF

      IF (eri_env%method == eri_method_gpw_ht) THEN
         ! We need a task list
         NULLIFY (task_list_sub)
         skip_load_balance_distributed = dft_control%qs_control%skip_load_balance_distributed
         CALL get_qs_env(qs_env, ks_env=ks_env, sab_orb=sab_orb_sub)
         CALL allocate_task_list(task_list_sub)
         CALL generate_qs_task_list(ks_env, task_list_sub, &
                                    reorder_rs_grid_ranks=.TRUE., soft_valid=.FALSE., &
                                    skip_load_balance_distributed=skip_load_balance_distributed, &
                                    pw_env_external=pw_env_sub, sab_orb_external=sab_orb_sub)
      END IF

      NULLIFY (wfn_r%pw, rho_g%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, wfn_r%pw, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, cell=cell, &
                      particle_set=particle_set, atomic_kind_set=atomic_kind_set)

      ! pre-calculate wavefunctions on reals space grid
      nspins = SIZE(mos)
      IF (eri_env%eri_gpw%store_wfn) THEN
         rsize = 0.0_dp
         nmo = 0
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, nmo=nx)
            nmo = MAX(nmo, nx)
            rsize = REAL(SIZE(wfn_r%pw%cr3d), KIND=dp)*nx
         END DO
         IF (print1 .AND. iw > 0) THEN
            rsize = rsize*8._dp/1000000._dp
            WRITE (iw, "(T4,'ERI_GPW|',' Store active orbitals on real space grid ',T63,F12.3,' MB')") rsize
         END IF
         ALLOCATE (wfn_a(nmo, nspins))
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
            DO iwfn = 1, nmo
               NULLIFY (wfn_a(iwfn, ispin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, wfn_a(iwfn, ispin)%pw, &
                                      use_data=REALDATA3D, &
                                      in_space=REALSPACE)
               CALL calculate_wavefunction(mo_coeff, iwfn, wfn_a(iwfn, ispin), rho_g, atomic_kind_set, &
                                           qs_kind_set, cell, dft_control, particle_set, pw_env_sub)
               IF (print2 .AND. iw > 0) THEN
                  WRITE (iw, "(T4,'ERI_GPW|',' Orbital stored ',I4,'  Spin ',i1)") iwfn, ispin
               END IF
            END DO
         END DO
      END IF

      ! get some of the grids ready
      NULLIFY (rho_r%pw, pot_g%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_r%pw, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, pot_g%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      ! run the FFT once, to set up buffers and to take into account the memory
      rho_r%pw%cr3d = 0.0D0
      CALL pw_transfer(rho_r%pw, rho_g%pw)
      dvol = rho_r%pw%pw_grid%dvol

      ! calculate the integrals
      intcount = 0
      DO isp1 = 1, nspins
         CALL get_mo_set(mo_set=mos(isp1)%mo_set, nmo=nmo1)
         nmm = (nmo1*(nmo1+1))/2
         irange = get_irange_csr(nmm, eri_env%eri(1)%csr_mat%mp_group)
         DO iwa1 = 1, nmo1
            IF (eri_env%eri_gpw%store_wfn) THEN
               wfn1 => wfn_a(iwa1, isp1)
            ELSE
               CPABORT("")
            END IF
            DO iwa2 = iwa1, nmo1
               iwa12 = get_index_ij(iwa1, iwa2, nmo1)
               IF (iwa12 >= irange(1) .AND. iwa12 <= irange(2)) THEN
                  iwa12 = iwa12-irange(1)+1
               ELSE
                  iwa12 = 0
               END IF
               IF (eri_env%eri_gpw%store_wfn) THEN
                  wfn2 => wfn_a(iwa2, isp1)
               ELSE
                  CPABORT("")
               END IF
               ! calculate charge distribution and potential
               rho_r%pw%cr3d = wfn1%pw%cr3d*wfn2%pw%cr3d
               IF (print2 .AND. iw > 0) THEN
                  erint = pw_integrate_function(rho_r%pw)/dvol
                  WRITE (iw, "(T4,'ERI_GPW| Integral rho_ab ',T32,2I4,' [',I1,']',T58,G20.14)") &
                     iwa1, iwa2, isp1, erint
               END IF
               CALL pw_transfer(rho_r%pw, rho_g%pw)
               CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_int, pot_g%pw)
               ! screening using pair_int
               IF (pair_int < eri_env%eps_integral) CYCLE
               CALL pw_transfer(pot_g%pw, rho_r%pw)
               !
               IF (eri_env%method == eri_method_gpw_ht) THEN
                  CPABORT("Not available")
               ELSEIF (eri_env%method == eri_method_full_gpw) THEN
                  DO isp2 = isp1, nspins
                     CALL get_mo_set(mo_set=mos(isp1)%mo_set, nmo=nmo2)
                     nx = (nmo2*(nmo2+1))/2
                     ALLOCATE (eri(nx), eri_index(nx))
                     icount2 = 0
                     iwbs = 1
                     IF (isp1 == isp2) iwbs = iwa1
                     isp = (isp1-1)*isp2-((isp1-1)*(isp1-2))/2+(isp2-isp1+1)
                     DO iwb1 = iwbs, nmo2
                        IF (eri_env%eri_gpw%store_wfn) THEN
                           wfn3 => wfn_a(iwb1, isp2)
                        ELSE
                           CPABORT("")
                        END IF
                        iwbt = iwb1
                        IF (isp1 == isp2 .AND. iwa1 == iwb1) iwbt = iwa2
                        DO iwb2 = iwbt, nmo2
                           IF (eri_env%eri_gpw%store_wfn) THEN
                              wfn4 => wfn_a(iwb2, isp2)
                           ELSE
                              CPABORT("")
                           END IF
                           wfn_r%pw%cr3d = rho_r%pw%cr3d*wfn3%pw%cr3d*wfn4%pw%cr3d
                           erint = pw_integrate_function(wfn_r%pw)
                           IF (erint > eri_env%eps_integral) THEN
                              intcount = intcount+1
                              IF (print2 .AND. iw > 0) THEN
                                 WRITE (iw, "(T4,'ERI_GPW|',T20,2I4,' [',I1,']',2I4,' [',I1,']',T58,G20.14)") &
                                    iwa1, iwa2, isp1, iwb1, iwb2, isp2, erint
                              END IF
                              icount2 = icount2+1
                              eri(icount2) = erint
                              eri_index(icount2) = get_index_ij(iwb1, iwb2, nmo2)
                           END IF
                        END DO
                     END DO
                     !
                     CALL update_csr_matrix(eri_env%eri(isp)%csr_mat, icount2, eri, eri_index, iwa12)
                     !
                     DEALLOCATE (eri, eri_index)
                  END DO
               ELSE
                  CPABORT("Unknown option")
               END IF
            END DO
         END DO
      END DO

      IF (print1 .AND. iw > 0) THEN
         WRITE (iw, "(T4,'ERI_GPW|',' Number of Integrals stored ',T68,I10)") intcount
      END IF

      IF (eri_env%eri_gpw%store_wfn) THEN
         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, nmo=nmo)
            DO iwfn = 1, nmo
               CALL pw_release(wfn_a(iwfn, ispin)%pw)
            END DO
         END DO
         DEALLOCATE (wfn_a)
      END IF
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wfn_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, pot_g%pw)

      IF (eri_env%method == eri_method_gpw_ht) THEN
         CALL deallocate_task_list(task_list_sub)
      END IF
      CALL pw_env_release(pw_env_sub)
      ! restore the initial value of the cutoff
      dft_control%qs_control%e_cutoff = e_cutoff_old
      dft_control%qs_control%cutoff = cutoff_old
      dft_control%qs_control%relative_cutoff = relative_cutoff_old
      DEALLOCATE (e_cutoff_old)

      CALL timestop(handle)

   END SUBROUTINE calculate_eri_gpw

! **************************************************************************************************
!> \brief Sets the Green's function
!> \param green ...
!> \param eri_env ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE pw_eri_green_create(green, eri_env)

      TYPE(greens_fn_type), POINTER                      :: green
      TYPE(eri_type)                                     :: eri_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_eri_green_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ig
      REAL(KIND=dp)                                      :: a, g2, g3d, gg, rlength
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_type), POINTER                             :: gf

      ! initialize influence function
      gf => green%influence_fn
      grid => green%influence_fn%pw_grid
      SELECT CASE (green%method)
      CASE (PERIODIC3D)

         SELECT CASE (eri_env%operator)
         CASE (eri_operator_coulomb)
            DO ig = grid%first_gne0, grid%ngpts_cut_local
               g2 = grid%gsq(ig)
               gf%cc(ig) = fourpi/g2
            END DO
            IF (grid%have_g0) gf%cc(1) = 0.0_dp
         CASE (eri_operator_yukawa)
            a = eri_env%operator_parameter**2
            DO ig = grid%first_gne0, grid%ngpts_cut_local
               g2 = grid%gsq(ig)
               gf%cc(ig) = fourpi/(a+g2)
            END DO
            IF (grid%have_g0) gf%cc(1) = fourpi/a
         CASE (eri_operator_erf)
            CPABORT("")
         CASE (eri_operator_erfc)
            CPABORT("")
         CASE (eri_operator_gaussian)
            CPABORT("")
         CASE DEFAULT
            CPABORT("")
         END SELECT

      CASE (ANALYTIC0D)

         SELECT CASE (eri_env%operator)
         CASE (eri_operator_coulomb)
            rlength = green%radius
            DO ig = grid%first_gne0, grid%ngpts_cut_local
               g2 = grid%gsq(ig)
               gg = SQRT(g2)
               g3d = fourpi/g2
               gf%cc(ig) = g3d*(1.0_dp-COS(rlength*gg))
            END DO
            IF (grid%have_g0) gf%cc(1) = 0.5_dp*fourpi*rlength*rlength
         CASE (eri_operator_yukawa)
            CPABORT("")
         CASE (eri_operator_erf)
            CPABORT("")
         CASE (eri_operator_erfc)
            CPABORT("")
         CASE (eri_operator_gaussian)
            CPABORT("")
         CASE DEFAULT
            CPABORT("")
         END SELECT

      CASE DEFAULT
         CPABORT("")
      END SELECT

   END SUBROUTINE pw_eri_green_create

! **************************************************************************************************
!> \brief calculates combined index (ij)
!> \param i ...
!> \param j ...
!> \param n ...
!> \retval ij ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   FUNCTION get_index_ij(i, j, n) RESULT(ij)

      INTEGER                                            :: i, j, n, ij

      CPASSERT(i <= j)
      CPASSERT(i <= n)
      CPASSERT(j <= n)

      ij = (i-1)*n-((i-1)*(i-2))/2+(j-i+1)

      CPASSERT(ij <= (n*(n+1))/2)

   END FUNCTION get_index_ij
! **************************************************************************************************
!> \brief extracts indeces i and j from combined index ij
!> \param ij ...
!> \param n ...
!> \param i ...
!> \param j ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE get_index_pair(ij, n, i, j)
      INTEGER                                            :: ij, n, i, j

      INTEGER                                            :: m, m0

      m = MAX(ij/n, 1)
      DO i = m, n
         m0 = (i-1)*n-((i-1)*(i-2))/2
         j = ij-m0+i-1
         IF (j <= n) EXIT
      END DO

      CPASSERT(i > 0 .AND. i <= n)
      CPASSERT(j > 0 .AND. j <= n)
      CPASSERT(i <= j)

   END SUBROUTINE get_index_pair
! **************************************************************************************************
!> \brief calculates index range for processor in group mp_group
!> \param nindex ...
!> \param mp_group ...
!> \retval irange ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   FUNCTION get_irange_csr(nindex, mp_group) RESULT(irange)

      INTEGER                                            :: nindex, mp_group
      INTEGER, DIMENSION(2)                              :: irange

      INTEGER                                            :: numtask, taskid
      REAL(KIND=dp)                                      :: rat

      CALL mp_environ(numtask, taskid, mp_group)

      IF (numtask == 1 .AND. taskid == 0) THEN
         irange(1) = 1
         irange(2) = nindex
      ELSEIF (numtask >= nindex) THEN
         IF (taskid >= nindex) THEN
            irange(1) = 1
            irange(2) = 0
         ELSE
            irange(1) = taskid+1
            irange(2) = taskid+1
         END IF
      ELSE
         rat = REAL(nindex, KIND=dp)/REAL(numtask, KIND=dp)
         irange(1) = NINT(rat*taskid)+1
         irange(2) = NINT(rat*taskid+rat)
      END IF

   END FUNCTION get_irange_csr

! **************************************************************************************************
!> \brief Adds data for a new row to the csr matrix
!> \param csr_mat ...
!> \param nnz ...
!> \param rdat ...
!> \param rind ...
!> \param irow ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE update_csr_matrix(csr_mat, nnz, rdat, rind, irow)

      TYPE(csr_type), INTENT(INOUT)                      :: csr_mat
      INTEGER, INTENT(IN)                                :: nnz
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rdat
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rind
      INTEGER, INTENT(IN)                                :: irow

      INTEGER                                            :: k, nrow, nze, nze_new

      IF (irow /= 0) THEN
         CPASSERT(irow <= csr_mat%nrows_total)
         CPASSERT(irow > csr_mat%nrows_local)
         nze = csr_mat%nze_local
         nze_new = nze+nnz
         ! values
         CALL reallocate(csr_mat%nzval_local%r_dp, 1, nze_new)
         csr_mat%nzval_local%r_dp(nze+1:nze_new) = rdat(1:nnz)
         ! col indices
         CALL reallocate(csr_mat%colind_local, 1, nze_new)
         csr_mat%colind_local(nze+1:nze_new) = rind(1:nnz)
         ! rows
         nrow = csr_mat%nrows_local
         CALL reallocate(csr_mat%rowptr_local, 1, irow+1)
         csr_mat%rowptr_local(nrow+1:irow) = nze+1
         csr_mat%rowptr_local(irow+1) = nze_new+1
         ! nzerow
         CALL reallocate(csr_mat%nzerow_local, 1, irow)
         DO k = nrow+1, irow
            csr_mat%nzerow_local(k) = csr_mat%rowptr_local(k+1)-csr_mat%rowptr_local(k)
         ENDDO
         csr_mat%nrows_local = irow
         csr_mat%nze_local = csr_mat%nze_local+nnz
      END IF
      csr_mat%nze_total = csr_mat%nze_total+nnz
      csr_mat%has_indices = .TRUE.

   END SUBROUTINE update_csr_matrix

! **************************************************************************************************
!> \brief Computes and prints the active orbitals on Cube Files
!> \param input ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param mos ...
! **************************************************************************************************
   SUBROUTINE print_orbital_cubes(input, qs_env, mos)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_p_type), DIMENSION(:)                  :: mos

      CHARACTER(len=*), PARAMETER :: routineN = 'print_orbital_cubes', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filebody, filename, title
      INTEGER                                            :: i, imo, isp, nmo, str(3), unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: alist, blist, istride
      LOGICAL                                            :: do_mo, explicit_a, explicit_b
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: wf_g, wf_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CALL section_vals_val_get(input, "FILENAME", c_val=filebody)
      CALL section_vals_val_get(input, "STRIDE", i_vals=istride)
      IF (SIZE(istride) == 1) THEN
         str(1:3) = istride(1)
      ELSEIF (SIZE(istride) == 3) THEN
         str(1:3) = istride(1:3)
      ELSE
         CPABORT("STRIDE arguments inconsistent")
      END IF
      CALL section_vals_val_get(input, "ALIST", i_vals=alist, explicit=explicit_a)
      CALL section_vals_val_get(input, "BLIST", i_vals=blist, explicit=explicit_b)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      subsys=subsys, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      particle_set=particle_set, &
                      pw_env=pw_env)
      CALL qs_subsys_get(subsys, particles=particles)
      !
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_r%pw, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_g%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      !
      DO isp = 1, SIZE(mos)
         CALL get_mo_set(mo_set=mos(isp)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
         DO imo = 1, nmo
            IF (isp == 1 .AND. explicit_a) THEN
               IF (alist(1) == -1) THEN
                  do_mo = .TRUE.
               ELSE
                  do_mo = .FALSE.
                  DO i = 1, SIZE(alist)
                     IF (imo == alist(i)) do_mo = .TRUE.
                  END DO
               END IF
            ELSE IF (isp == 2 .AND. explicit_b) THEN
               IF (blist(1) == -1) THEN
                  do_mo = .TRUE.
               ELSE
                  do_mo = .FALSE.
                  DO i = 1, SIZE(blist)
                     IF (imo == blist(i)) do_mo = .TRUE.
                  END DO
               END IF
            ELSE
               do_mo = .TRUE.
            END IF
            IF (.NOT. do_mo) CYCLE
            CALL calculate_wavefunction(mo_coeff, imo, wf_r, wf_g, atomic_kind_set, &
                                        qs_kind_set, cell, dft_control, particle_set, pw_env)
            IF (para_env%ionode) THEN
               WRITE (filename, '(A,A1,I4.4,A1,I1.1,A)') TRIM(filebody), "_", imo, "_", isp, ".cube"
               CALL open_file(filename, unit_number=unit_nr, file_status="UNKNOWN", file_action="WRITE")
               WRITE (title, *) "Active Orbital ", imo, " spin ", isp
            ELSE
               unit_nr = -1
            END IF
            CALL cp_pw_to_cube(wf_r%pw, unit_nr, title, particles=particles, stride=istride)
            IF (para_env%ionode) THEN
               CALL close_file(unit_nr)
            END IF
         ENDDO
      ENDDO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g%pw)

   END SUBROUTINE print_orbital_cubes
! **************************************************************************************************
!> \brief Writes a FCIDUMP file
!> \param active_space_env ...
!> \param as_input ...
!> \param eref ...
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE fcidump(active_space_env, as_input, eref)

      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(section_vals_type), POINTER                   :: as_input
      REAL(KIND=dp), INTENT(IN)                          :: eref

      CHARACTER(len=*), PARAMETER :: routineN = 'fcidump', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, i1, i12, i12l, i2, i3, i34, i34l, i4, &
                                                            irange(2), irptr, isym, iw, ms2, &
                                                            nelec, nindex, norb, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: colind
      REAL(KIND=dp)                                      :: eeri, erint, esub, mval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: erival
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ks_mat, ks_ref, p_mat
      TYPE(cp_fm_type), POINTER                          :: matrix
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(csr_type), POINTER                            :: eri

      logger => cp_get_default_logger()
      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           as_input, "FCIDUMP"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, as_input, "FCIDUMP", &
                                   extension=".fcidump", file_status="REPLACE", file_action="WRITE", file_form="FORMATTED")
         !
         nspins = active_space_env%nspins
         IF (nspins == 1) THEN
            ms2 = 0 !active_space_env%multiplicity ?
            nelec = active_space_env%nactive
            CALL get_mo_set(active_space_env%mos_active(1)%mo_set, nmo=norb)
            IF (iw > 0) THEN
               WRITE (iw, "(A,A,I4,A,I4,A,I2,A)") "&FCI", " NORB=", norb, ",NELEC=", nelec, ",MS2=", ms2, ","
               isym = 1
               WRITE (iw, "(A,1000(I1,','))") "  ORBSYM=", (isym, i=1, norb)
               isym = 0
               WRITE (iw, "(A,I1,A)") "  ISYM=", isym, ","
               WRITE (iw, "(A)") " /"
            END IF
            !
            ! Loop over ERI, calculate subspace HF energy and Fock matrix
            !
            ! replicate KS, Core, and P matrices
            ALLOCATE (ks_mat(norb, norb), ks_ref(norb, norb), p_mat(norb, norb))
            ks_mat = 0.0_dp
            ks_ref = 0.0_dp
            p_mat = 0.0_dp
            DO i1 = 1, norb
               DO i2 = i1, norb
                  matrix => active_space_env%p_ref(1)%matrix
                  CALL cp_fm_get_element(matrix, i1, i2, mval)
                  p_mat(i1, i2) = mval
                  p_mat(i2, i1) = mval
                  matrix => active_space_env%ks_sub(1)%matrix
                  CALL cp_fm_get_element(matrix, i1, i2, mval)
                  ks_mat(i1, i2) = mval
                  ks_mat(i2, i1) = mval
               END DO
            END DO
            ! calculate energy
            eeri = 0.0_dp
            eri => active_space_env%eri%eri(1)%csr_mat
            nindex = (norb*(norb+1))/2
            irange = get_irange_csr(nindex, eri%mp_group)
            DO i1 = 1, norb
               DO i2 = i1, norb
                  i12 = get_index_ij(i1, i2, norb)
                  IF (i12 >= irange(1) .AND. i12 <= irange(2)) THEN
                     i12l = i12-irange(1)+1
                     irptr = eri%rowptr_local(i12l)-1
                     DO i34l = 1, eri%nzerow_local(i12l)
                        i34 = eri%colind_local(irptr+i34l)
                        CALL get_index_pair(i34, norb, i3, i4)
                        erint = eri%nzval_local%r_dp(irptr+i34l)
                        ! Coulomb
                        ks_ref(i1, i2) = ks_ref(i1, i2)+erint*p_mat(i3, i4)
                        IF (i3 /= i4) THEN
                           ks_ref(i1, i2) = ks_ref(i1, i2)+erint*p_mat(i3, i4)
                        END IF
                        IF (i12 /= i34) THEN
                           ks_ref(i3, i4) = ks_ref(i3, i4)+erint*p_mat(i1, i2)
                           IF (i1 /= i2) THEN
                              ks_ref(i3, i4) = ks_ref(i3, i4)+erint*p_mat(i1, i2)
                           END IF
                        END IF
                        ! Exchange
                        erint = -0.5_dp*erint
                        ks_ref(i1, i3) = ks_ref(i1, i3)+erint*p_mat(i2, i4)
                        IF (i1 /= i2) THEN
                           ks_ref(i2, i3) = ks_ref(i2, i3)+erint*p_mat(i1, i4)
                        END IF
                        IF (i3 /= i4) THEN
                           ks_ref(i1, i4) = ks_ref(i1, i4)+erint*p_mat(i2, i3)
                        END IF
                        IF (i1 /= i2 .AND. i3 /= i4) THEN
                           ks_ref(i2, i4) = ks_ref(i2, i4)+erint*p_mat(i1, i3)
                        END IF
                     END DO
                  END IF
               END DO
            END DO
            !
            DO i1 = 1, norb
               DO i2 = i1, norb
                  ks_ref(i1, i2) = ks_ref(i2, i1)
               END DO
            END DO
            CALL mp_sum(ks_ref, eri%mp_group)
            !
            eeri = 0.5_dp*SUM(ks_ref*p_mat)
            ks_mat(1:norb, 1:norb) = ks_mat(1:norb, 1:norb)+ks_ref(1:norb, 1:norb)
            esub = eref-SUM(ks_mat(1:norb, 1:norb)*p_mat(1:norb, 1:norb))+eeri
            !
            ! Print integrals: ERI
            eri => active_space_env%eri%eri(1)%csr_mat
            nindex = (norb*(norb+1))/2
            irange = get_irange_csr(nindex, eri%mp_group)
            ALLOCATE (erival(nindex), colind(nindex))
            DO i1 = 1, norb
               DO i2 = i1, norb
                  i12 = get_index_ij(i1, i2, norb)
                  IF (i12 >= irange(1) .AND. i12 <= irange(2)) THEN
                     i12l = i12-irange(1)+1
                     irptr = eri%rowptr_local(i12l)
                     nindex = eri%nzerow_local(i12l)
                     colind(1:nindex) = eri%colind_local(irptr:irptr+nindex-1)
                     erival(1:nindex) = eri%nzval_local%r_dp(irptr:irptr+nindex-1)
                  ELSE
                     erival = 0.0_dp
                     colind = 0
                     nindex = 0
                  END IF
                  CALL mp_sum(nindex, eri%mp_group)
                  CALL mp_sum(colind(1:nindex), eri%mp_group)
                  CALL mp_sum(erival(1:nindex), eri%mp_group)
                  CALL mp_sync(eri%mp_group)
                  IF (iw > 0) THEN
                     DO i34l = 1, nindex
                        i34 = colind(i34l)
                        erint = erival(i34l)
                        CALL get_index_pair(i34, norb, i3, i4)
                        WRITE (iw, "(ES23.16,4I4)") erint, i1, i2, i3, i4
                     END DO
                  END IF
               END DO
            END DO
            DEALLOCATE (erival, colind)
            ! Print integrals: Fij
            IF (iw > 0) THEN
               i3 = 0; i4 = 0
               DO i1 = 1, norb
                  DO i2 = i1, norb
                     WRITE (iw, "(ES23.16,4I4)") ks_mat(i1, i2), i1, i2, i3, i4
                  END DO
               END DO
            END IF
            ! Print energy
            i1 = 0; i2 = 0; i3 = 0; i4 = 0
            IF (iw > 0) WRITE (iw, "(ES23.16,4I4)") esub, i1, i2, i3, i4
            DEALLOCATE (ks_mat, ks_ref, p_mat)
         ELSE
            CPABORT("Option not yet available")
         END IF
         !
         CALL cp_print_key_finished_output(iw, logger, as_input, "FCIDUMP")
      END IF

   END SUBROUTINE fcidump
! **************************************************************************************************
!> \brief Creates a local basis
!> \param pro_basis_set ...
!> \param zval ...
!> \param ishell ...
!> \param nshell ...
!> \param lnam ...
!> \par History
!>      05.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE create_pro_basis(pro_basis_set, zval, ishell, nshell, lnam)
      TYPE(gto_basis_set_type), POINTER                  :: pro_basis_set
      INTEGER, INTENT(IN)                                :: zval, ishell
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nshell
      CHARACTER(len=*), DIMENSION(:), INTENT(IN)         :: lnam

      CHARACTER(len=*), PARAMETER :: routineN = 'create_pro_basis', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=6), DIMENSION(:), POINTER            :: sym
      INTEGER                                            :: i, l, nj
      INTEGER, DIMENSION(4, 7)                           :: ne
      INTEGER, DIMENSION(:), POINTER                     :: lq, nq
      REAL(KIND=dp), DIMENSION(:), POINTER               :: zet
      TYPE(sto_basis_set_type), POINTER                  :: sto_basis_set

      CPASSERT(.NOT. ASSOCIATED(pro_basis_set))
      NULLIFY (sto_basis_set)

      ! electronic configuration
      ne = 0
      DO l = 1, 4 !lq(1)+1
         nj = 2*(l-1)+1
         DO i = l, 7 ! nq(1)
            ne(l, i) = ptable(zval)%e_conv(l-1)-2*nj*(i-l)
            ne(l, i) = MAX(ne(l, i), 0)
            ne(l, i) = MIN(ne(l, i), 2*nj)
         END DO
      END DO
      ALLOCATE (nq(ishell), lq(ishell), zet(ishell), sym(ishell))
      DO i = 1, ishell
         nq(i) = nshell(i)
         SELECT CASE (lnam (i))
         CASE ('S', 's')
            lq(i) = 0
         CASE ('P', 'p')
            lq(i) = 1
         CASE ('D', 'd')
            lq(i) = 2
         CASE ('F', 'f')
            lq(i) = 3
         CASE DEFAULT
            CPABORT("Wrong l QN")
         END SELECT
         sym(i) = lnam(i)
         zet(i) = srules(zval, ne, nq(1), lq(1))
      END DO
      CALL allocate_sto_basis_set(sto_basis_set)
      CALL set_sto_basis_set(sto_basis_set, nshell=1, nq=nq, lq=lq, zet=zet, symbol=sym)
      CALL create_gto_from_sto_basis(sto_basis_set, pro_basis_set, 6)
      pro_basis_set%norm_type = 2
      CALL init_orb_basis_set(pro_basis_set)
      CALL deallocate_sto_basis_set(sto_basis_set)

   END SUBROUTINE create_pro_basis

END MODULE qs_active_space_methods
