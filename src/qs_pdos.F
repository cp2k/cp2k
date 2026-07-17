!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation and writing of projected density of  states
!>         The DOS is computed per angular momentum and per kind
!> \par History
!>      -
!> \author Marcella (29.02.2008,MK)
! **************************************************************************************************
MODULE qs_pdos
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: lowdin_kp_mo_coeff
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: nso,&
                                              nsoset
   USE orbital_symbols,                 ONLY: l_sym,&
                                              sgf_symbol
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_dos_utils,                    ONLY: &
        add_broadened_value, broadening_cutoff, broadening_function, dos_density_scale, &
        dos_energy_label, dos_energy_scale, dos_energy_unit_ev, dos_energy_zero_absolute, &
        dos_energy_zero_auto, dos_energy_zero_hoco, dos_energy_zero_label, &
        dos_resolve_energy_zero, write_broadening_info
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_scf_diagonalization,          ONLY: diag_kp_smat
   USE qs_scf_types,                    ONLY: qs_scf_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_pdos'

! **************************************************************************************************
   ! *** Public subroutines ***

   PUBLIC :: calculate_projected_dos, calculate_projected_dos_kp

   TYPE ldos_type
      INTEGER :: maxl = -1, nlist = -1
      LOGICAL :: separate_components = .FALSE.
      INTEGER, DIMENSION(:), POINTER           :: list_index => NULL()
      REAL(KIND=dp), DIMENSION(:, :), &
         POINTER                                :: pdos_array => NULL()
   END TYPE ldos_type

   TYPE r_ldos_type
      INTEGER :: nlist = -1, npoints = -1
      INTEGER, DIMENSION(:, :), POINTER :: index_grid_local => NULL()
      INTEGER, DIMENSION(:), POINTER           :: list_index => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER     :: x_range => NULL(), y_range => NULL(), z_range => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER     :: eval_range => NULL()
      REAL(KIND=dp), DIMENSION(:), &
         POINTER                                 :: pdos_array => NULL()
   END TYPE r_ldos_type

   TYPE ldos_p_type
      TYPE(ldos_type), POINTER :: ldos => NULL()
   END TYPE ldos_p_type

   TYPE r_ldos_p_type
      TYPE(r_ldos_type), POINTER :: ldos => NULL()
   END TYPE r_ldos_p_type
CONTAINS

! **************************************************************************************************
!> \brief   Compute and write projected density of states
!> \param mo_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param qs_env ...
!> \param dft_section ...
!> \param ispin ...
!> \param xas_mittle ...
!> \param external_matrix_shalf ...
!> \param unoccupied_orbs ...
!> \param unoccupied_evals ...
!> \param pdos_print_key ...
!> \param write_pdos ...
!> \param write_pdos_curve ...
!> \date    26.02.2008
!> \par History:
!>       - Added optional external matrix_shalf to avoid recomputing it (A. Bussy, 09.2019)
!> \par Variables
!>       -
!>       -
!> \author  MI
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE calculate_projected_dos(mo_set, atomic_kind_set, qs_kind_set, particle_set, qs_env, &
                                      dft_section, ispin, xas_mittle, external_matrix_shalf, &
                                      unoccupied_orbs, unoccupied_evals, pdos_print_key, write_pdos, write_pdos_curve)

      TYPE(mo_set_type), INTENT(IN)                      :: mo_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: dft_section
      INTEGER, INTENT(IN), OPTIONAL                      :: ispin
      CHARACTER(LEN=default_string_length), INTENT(IN), &
         OPTIONAL                                        :: xas_mittle
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL, TARGET     :: external_matrix_shalf, unoccupied_orbs
      TYPE(cp_1d_r_p_type), INTENT(IN), OPTIONAL, TARGET :: unoccupied_evals
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: pdos_print_key
      LOGICAL, INTENT(IN), OPTIONAL                      :: write_pdos, write_pdos_curve

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_projected_dos'

      CHARACTER(LEN=16)                                  :: energy_label, fmtstr2
      CHARACTER(LEN=27)                                  :: fmtstr1
      CHARACTER(LEN=32)                                  :: zero_label
      CHARACTER(LEN=6), ALLOCATABLE, DIMENSION(:, :, :)  :: tmp_str
      CHARACTER(LEN=default_string_length)               :: kind_name, my_act, my_mittle, my_pos, &
                                                            my_print_key, spin(2)
      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:)                       :: ldos_index, r_ldos_index
      INTEGER :: broaden_type, energy_unit, energy_zero, handle, homo, i, iatom, ikind, il, ildos, &
         im, imo, imo_ref, in_x, in_y, in_z, ir, irow, iset, isgf, ishell, iso, ispin_ref, &
         iterstep, iw, j, jx, jy, jz, k, lcomponent, lshell, maxl, maxlgto, my_spin, n_dependent, &
         n_r_ldos, n_rep, nao, natom, ncol_global, ndigits, nkind, nldos, nmo, nmo_ref, np_tot, &
         npoints, nrow_global, nset, nsgf, nvirt, out_each, output_unit, resolved_energy_zero
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: firstrow
      INTEGER, DIMENSION(:), POINTER                     :: list, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: bo, l
      LOGICAL :: append, calc_matsh, do_curve, do_ldos, do_r_ldos, do_virt, fractional_occupation, &
         ionode, separate_components, should_output, write_curve, write_pdos_file
      LOGICAL, DIMENSION(:, :), POINTER                  :: read_r
      REAL(KIND=dp) :: broaden_width, de, dh(3, 3), dvol, e_fermi, e_fermi_ref(2), energy_factor, &
         energy_ref, ev_factor, hoco, hoco_ref(2), r(3), r_vec(3), ratom(3), voigt_mixing
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, eval_ref, evals_virt, &
                                                            occ_ref, occupation_numbers
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vecbuffer
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: pdos_array
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: matrix_shalfc, matrix_work
      TYPE(cp_fm_type), POINTER                          :: matrix_shalf, mo_coeff, mo_virt
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: s_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(ldos_p_type), DIMENSION(:), POINTER           :: ldos_p
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos_ref
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_c1d_gs_type)                               :: wf_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: wf_r
      TYPE(r_ldos_p_type), DIMENSION(:), POINTER         :: r_ldos_p
      TYPE(section_vals_type), POINTER                   :: curve_section, ldos_section

      NULLIFY (logger, mos_ref, eval_ref, occ_ref)
      logger => cp_get_default_logger()
      ionode = logger%para_env%is_source()
      my_print_key = "PRINT%PDOS"
      IF (PRESENT(pdos_print_key)) my_print_key = TRIM(pdos_print_key)
      write_pdos_file = .TRUE.
      IF (PRESENT(write_pdos)) write_pdos_file = write_pdos
      curve_section => section_vals_get_subs_vals(dft_section, TRIM(my_print_key)//"%CURVE")
      CALL section_vals_get(curve_section, explicit=write_curve)
      IF (PRESENT(write_pdos_curve)) write_curve = write_pdos_curve
      should_output = BTEST(cp_print_key_should_output(logger%iter_info, dft_section, &
                                                       TRIM(my_print_key)), cp_p_file)
      output_unit = cp_logger_get_default_io_unit(logger)

      spin(1) = "ALPHA"
      spin(2) = "BETA"
      IF ((.NOT. should_output)) RETURN

      NULLIFY (context, s_matrix, orb_basis_set, para_env, pdos_array)
      NULLIFY (eigenvalues, fm_struct_tmp, mo_coeff, vecbuffer, mo_virt)
      NULLIFY (curve_section, ldos_section, list, cell, pw_env, auxbas_pw_pool, evals_virt)
      NULLIFY (occupation_numbers, ldos_p, r_ldos_p, dft_control, occupation_numbers)

      CALL timeset(routineN, handle)
      iterstep = logger%iter_info%iteration(logger%iter_info%n_rlevel)

      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT='(/,(T3,A,T61,I10))') &
         " Calculate PDOS at iteration step ", iterstep
      CALL get_qs_env(qs_env=qs_env, &
                      matrix_s=s_matrix, &
                      dft_control=dft_control)

      CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf, maxlgto=maxlgto)
      nkind = SIZE(atomic_kind_set)

      CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, homo=homo, nao=nao, nmo=nmo, &
                      mu=e_fermi)
      CALL cp_fm_get_info(mo_coeff, &
                          context=context, para_env=para_env, &
                          nrow_global=nrow_global, &
                          ncol_global=ncol_global)

      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%OUT_EACH_MO", i_val=out_each)
      IF (out_each == -1) out_each = nao + 1
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%DELTA_E", r_val=de)
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%CURVE%BROADEN%TYPE", i_val=broaden_type)
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%CURVE%BROADEN%WIDTH", r_val=broaden_width)
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%CURVE%BROADEN%VOIGT_MIXING", r_val=voigt_mixing)
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%NDIGITS", i_val=ndigits)
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%CURVE%ENERGY_UNIT", i_val=energy_unit)
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%CURVE%ENERGY_ZERO", i_val=energy_zero)
      ndigits = MIN(MAX(ndigits, 1), 10)
      IF (write_curve .AND. de <= 0.0_dp) THEN
         CPWARN("Broadened PDOS output requires DELTA_E > 0 and will be skipped")
         write_curve = .FALSE.
      END IF
      IF (write_curve .AND. broaden_width <= 0.0_dp) THEN
         CPWARN("Broadened PDOS output requires a finite WIDTH and will be skipped")
         write_curve = .FALSE.
      END IF
      do_curve = write_curve .AND. (broaden_width > 0.0_dp)
      IF (do_curve) de = MAX(de, 0.00001_dp)
      nvirt = 0
      NULLIFY (evals_virt)
      IF (PRESENT(unoccupied_orbs) .AND. PRESENT(unoccupied_evals)) THEN
         IF (ASSOCIATED(unoccupied_evals%array)) THEN
            nvirt = SIZE(unoccupied_evals%array)
            IF (nvirt > 0) THEN
               mo_virt => unoccupied_orbs
               evals_virt => unoccupied_evals%array
            END IF
         END IF
      END IF
      do_virt = (nvirt > 0)

      calc_matsh = .TRUE.
      IF (PRESENT(external_matrix_shalf)) calc_matsh = .FALSE.

      ! Create S^1/2 : from sparse to full matrix, if no external available
      IF (calc_matsh) THEN
         NULLIFY (matrix_shalf)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=context, &
                                  nrow_global=nrow_global, ncol_global=nrow_global)
         ALLOCATE (matrix_shalf)
         CALL cp_fm_create(matrix_shalf, fm_struct_tmp, name="matrix_shalf")
         CALL cp_fm_create(matrix_work, fm_struct_tmp, name="matrix_work")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL copy_dbcsr_to_fm(s_matrix(1)%matrix, matrix_shalf)
         CALL cp_fm_power(matrix_shalf, matrix_work, 0.5_dp, EPSILON(0.0_dp), n_dependent)
         CALL cp_fm_release(matrix_work)
      ELSE
         matrix_shalf => external_matrix_shalf
      END IF

      ! Multiply S^(1/2) time the mOS coefficients to get orthonormalized MOS
      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=context, &
                               nrow_global=nrow_global, ncol_global=ncol_global)
      CALL cp_fm_create(matrix_shalfc, fm_struct_tmp, name="matrix_shalfc")
      CALL parallel_gemm("N", "N", nrow_global, ncol_global, nrow_global, &
                         1.0_dp, matrix_shalf, mo_coeff, 0.0_dp, matrix_shalfc)
      CALL cp_fm_struct_release(fm_struct_tmp)

      IF (do_virt) THEN
         IF (output_unit > 0) WRITE (UNIT=output_unit, FMT='(/,(T3,A,T14,I10,T27,A))') &
            " Use ", nvirt, " additional unoccupied KS orbitals"
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=context, &
                                  nrow_global=nrow_global, ncol_global=nvirt)
         CALL cp_fm_create(matrix_work, fm_struct_tmp, name="matrix_shalfc")
         CALL parallel_gemm("N", "N", nrow_global, nvirt, nrow_global, &
                            1.0_dp, matrix_shalf, mo_virt, 0.0_dp, matrix_work)
         CALL cp_fm_struct_release(fm_struct_tmp)
      END IF

      IF (calc_matsh) THEN
         CALL cp_fm_release(matrix_shalf)
         DEALLOCATE (matrix_shalf)
      END IF
      ! Array to store the PDOS per kind and angular momentum
      do_ldos = .FALSE.
      ldos_section => section_vals_get_subs_vals(dft_section, TRIM(my_print_key)//"%LDOS")

      CALL section_vals_get(ldos_section, n_repetition=nldos)
      IF (nldos > 0) THEN
         IF (output_unit > 0) WRITE (UNIT=output_unit, FMT='(/,(T3,A,T61,I10))') &
            " Prepare the list of atoms for LDOS.   Number of lists: ", nldos
         do_ldos = .TRUE.
         ALLOCATE (ldos_p(nldos))
         ALLOCATE (ldos_index(nldos))
         DO ildos = 1, nldos
            WRITE (ldos_index(ildos), '(I0)') ildos
            ALLOCATE (ldos_p(ildos)%ldos)
            NULLIFY (ldos_p(ildos)%ldos%pdos_array)
            NULLIFY (ldos_p(ildos)%ldos%list_index)

            CALL section_vals_val_get(ldos_section, "LIST", i_rep_section=ildos, n_rep_val=n_rep)
            IF (n_rep > 0) THEN
               ldos_p(ildos)%ldos%nlist = 0
               DO ir = 1, n_rep
                  NULLIFY (list)
                  CALL section_vals_val_get(ldos_section, "LIST", i_rep_section=ildos, i_rep_val=ir, &
                                            i_vals=list)
                  IF (ASSOCIATED(list)) THEN
                     CALL reallocate(ldos_p(ildos)%ldos%list_index, 1, ldos_p(ildos)%ldos%nlist + SIZE(list))
                     DO i = 1, SIZE(list)
                        ldos_p(ildos)%ldos%list_index(i + ldos_p(ildos)%ldos%nlist) = list(i)
                     END DO
                     ldos_p(ildos)%ldos%nlist = ldos_p(ildos)%ldos%nlist + SIZE(list)
                  END IF
               END DO
            ELSE
               ! stop, LDOS without list of atoms is not implemented
            END IF

            IF (output_unit > 0) WRITE (UNIT=output_unit, FMT='((T10,A,T18,I6,T25,A,T36,I10,A))') &
               " List ", ildos, " contains ", ldos_p(ildos)%ldos%nlist, " atoms"
            CALL section_vals_val_get(ldos_section, "COMPONENTS", i_rep_section=ildos, &
                                      l_val=ldos_p(ildos)%ldos%separate_components)
            IF (ldos_p(ildos)%ldos%separate_components) THEN
               ALLOCATE (ldos_p(ildos)%ldos%pdos_array(nsoset(maxlgto), nmo + nvirt))
            ELSE
               ALLOCATE (ldos_p(ildos)%ldos%pdos_array(0:maxlgto, nmo + nvirt))
            END IF
            ldos_p(ildos)%ldos%pdos_array = 0.0_dp
            ldos_p(ildos)%ldos%maxl = -1

         END DO
      END IF

      do_r_ldos = .FALSE.
      ldos_section => section_vals_get_subs_vals(dft_section, TRIM(my_print_key)//"%R_LDOS")
      CALL section_vals_get(ldos_section, n_repetition=n_r_ldos)
      IF (n_r_ldos > 0) THEN
         do_r_ldos = .TRUE.
         IF (output_unit > 0) WRITE (UNIT=output_unit, FMT='(/,(T3,A,T61,I10))') &
            " Prepare the list of points for R_LDOS.   Number of lists: ", n_r_ldos
         ALLOCATE (r_ldos_p(n_r_ldos))
         ALLOCATE (r_ldos_index(n_r_ldos))
         CALL get_qs_env(qs_env=qs_env, &
                         cell=cell, &
                         dft_control=dft_control, &
                         pw_env=pw_env)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                         pw_pools=pw_pools)

         CALL auxbas_pw_pool%create_pw(wf_r)
         CALL auxbas_pw_pool%create_pw(wf_g)
         ALLOCATE (read_r(4, n_r_ldos))
         DO ildos = 1, n_r_ldos
            WRITE (r_ldos_index(ildos), '(I0)') ildos
            ALLOCATE (r_ldos_p(ildos)%ldos)
            NULLIFY (r_ldos_p(ildos)%ldos%pdos_array)
            NULLIFY (r_ldos_p(ildos)%ldos%list_index)

            CALL section_vals_val_get(ldos_section, "LIST", i_rep_section=ildos, n_rep_val=n_rep)
            IF (n_rep > 0) THEN
               r_ldos_p(ildos)%ldos%nlist = 0
               DO ir = 1, n_rep
                  NULLIFY (list)
                  CALL section_vals_val_get(ldos_section, "LIST", i_rep_section=ildos, i_rep_val=ir, &
                                            i_vals=list)
                  IF (ASSOCIATED(list)) THEN
                     CALL reallocate(r_ldos_p(ildos)%ldos%list_index, 1, r_ldos_p(ildos)%ldos%nlist + SIZE(list))
                     DO i = 1, SIZE(list)
                        r_ldos_p(ildos)%ldos%list_index(i + r_ldos_p(ildos)%ldos%nlist) = list(i)
                     END DO
                     r_ldos_p(ildos)%ldos%nlist = r_ldos_p(ildos)%ldos%nlist + SIZE(list)
                  END IF
               END DO
            ELSE
               ! stop, LDOS without list of atoms is not implemented
            END IF

            ALLOCATE (r_ldos_p(ildos)%ldos%pdos_array(nmo + nvirt))
            r_ldos_p(ildos)%ldos%pdos_array = 0.0_dp
            read_r(1:3, ildos) = .FALSE.
            CALL section_vals_val_get(ldos_section, "XRANGE", i_rep_section=ildos, explicit=read_r(1, ildos))
            IF (read_r(1, ildos)) THEN
               CALL section_vals_val_get(ldos_section, "XRANGE", i_rep_section=ildos, r_vals= &
                                         r_ldos_p(ildos)%ldos%x_range)
            ELSE
               ALLOCATE (r_ldos_p(ildos)%ldos%x_range(2))
               r_ldos_p(ildos)%ldos%x_range = 0.0_dp
            END IF
            CALL section_vals_val_get(ldos_section, "YRANGE", i_rep_section=ildos, explicit=read_r(2, ildos))
            IF (read_r(2, ildos)) THEN
               CALL section_vals_val_get(ldos_section, "YRANGE", i_rep_section=ildos, r_vals= &
                                         r_ldos_p(ildos)%ldos%y_range)
            ELSE
               ALLOCATE (r_ldos_p(ildos)%ldos%y_range(2))
               r_ldos_p(ildos)%ldos%y_range = 0.0_dp
            END IF
            CALL section_vals_val_get(ldos_section, "ZRANGE", i_rep_section=ildos, explicit=read_r(3, ildos))
            IF (read_r(3, ildos)) THEN
               CALL section_vals_val_get(ldos_section, "ZRANGE", i_rep_section=ildos, r_vals= &
                                         r_ldos_p(ildos)%ldos%z_range)
            ELSE
               ALLOCATE (r_ldos_p(ildos)%ldos%z_range(2))
               r_ldos_p(ildos)%ldos%z_range = 0.0_dp
            END IF

            CALL section_vals_val_get(ldos_section, "ERANGE", i_rep_section=ildos, explicit=read_r(4, ildos))
            IF (read_r(4, ildos)) THEN
               CALL section_vals_val_get(ldos_section, "ERANGE", i_rep_section=ildos, &
                                         r_vals=r_ldos_p(ildos)%ldos%eval_range)
            ELSE
               ALLOCATE (r_ldos_p(ildos)%ldos%eval_range(2))
               r_ldos_p(ildos)%ldos%eval_range(1) = -HUGE(0.0_dp)
               r_ldos_p(ildos)%ldos%eval_range(2) = +HUGE(0.0_dp)
            END IF

            bo => wf_r%pw_grid%bounds_local
            dh = wf_r%pw_grid%dh
            dvol = wf_r%pw_grid%dvol
            np_tot = wf_r%pw_grid%npts(1)*wf_r%pw_grid%npts(2)*wf_r%pw_grid%npts(3)
            ALLOCATE (r_ldos_p(ildos)%ldos%index_grid_local(3, np_tot))

            r_ldos_p(ildos)%ldos%npoints = 0
            DO jz = bo(1, 3), bo(2, 3)
            DO jy = bo(1, 2), bo(2, 2)
            DO jx = bo(1, 1), bo(2, 1)
               !compute the position of the grid point
               i = jx - wf_r%pw_grid%bounds(1, 1)
               j = jy - wf_r%pw_grid%bounds(1, 2)
               k = jz - wf_r%pw_grid%bounds(1, 3)
               r(3) = k*dh(3, 3) + j*dh(3, 2) + i*dh(3, 1)
               r(2) = k*dh(2, 3) + j*dh(2, 2) + i*dh(2, 1)
               r(1) = k*dh(1, 3) + j*dh(1, 2) + i*dh(1, 1)

               DO il = 1, r_ldos_p(ildos)%ldos%nlist
                  iatom = r_ldos_p(ildos)%ldos%list_index(il)
                  ratom = particle_set(iatom)%r
                  r_vec = pbc(ratom, r, cell)
                  IF (cell%orthorhombic) THEN
                     IF (cell%perd(1) == 0) r_vec(1) = MODULO(r_vec(1), cell%hmat(1, 1))
                     IF (cell%perd(2) == 0) r_vec(2) = MODULO(r_vec(2), cell%hmat(2, 2))
                     IF (cell%perd(3) == 0) r_vec(3) = MODULO(r_vec(3), cell%hmat(3, 3))
                  END IF

                  in_x = 0
                  in_y = 0
                  in_z = 0
                  IF (r_ldos_p(ildos)%ldos%x_range(1) /= 0.0_dp) THEN
                     IF (r_vec(1) > r_ldos_p(ildos)%ldos%x_range(1) .AND. &
                         r_vec(1) < r_ldos_p(ildos)%ldos%x_range(2)) THEN
                        in_x = 1
                     END IF
                  ELSE
                     in_x = 1
                  END IF
                  IF (r_ldos_p(ildos)%ldos%y_range(1) /= 0.0_dp) THEN
                     IF (r_vec(2) > r_ldos_p(ildos)%ldos%y_range(1) .AND. &
                         r_vec(2) < r_ldos_p(ildos)%ldos%y_range(2)) THEN
                        in_y = 1
                     END IF
                  ELSE
                     in_y = 1
                  END IF
                  IF (r_ldos_p(ildos)%ldos%z_range(1) /= 0.0_dp) THEN
                     IF (r_vec(3) > r_ldos_p(ildos)%ldos%z_range(1) .AND. &
                         r_vec(3) < r_ldos_p(ildos)%ldos%z_range(2)) THEN
                        in_z = 1
                     END IF
                  ELSE
                     in_z = 1
                  END IF
                  IF (in_x*in_y*in_z > 0) THEN
                     r_ldos_p(ildos)%ldos%npoints = r_ldos_p(ildos)%ldos%npoints + 1
                     r_ldos_p(ildos)%ldos%index_grid_local(1, r_ldos_p(ildos)%ldos%npoints) = jx
                     r_ldos_p(ildos)%ldos%index_grid_local(2, r_ldos_p(ildos)%ldos%npoints) = jy
                     r_ldos_p(ildos)%ldos%index_grid_local(3, r_ldos_p(ildos)%ldos%npoints) = jz
                     EXIT
                  END IF
               END DO
            END DO
            END DO
            END DO
            CALL reallocate(r_ldos_p(ildos)%ldos%index_grid_local, 1, 3, 1, r_ldos_p(ildos)%ldos%npoints)
            npoints = r_ldos_p(ildos)%ldos%npoints
            CALL para_env%sum(npoints)
            IF (output_unit > 0) WRITE (UNIT=output_unit, FMT='((T10,A,T18,I6,T25,A,T36,I10,A))') &
               " List ", ildos, " contains ", npoints, " grid points"
         END DO
      END IF

      IF (TRIM(my_print_key) == "PRINT%DOS") THEN
         CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%PDOS%COMPONENTS", l_val=separate_components)
      ELSE
         CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%COMPONENTS", l_val=separate_components)
      END IF
      IF (separate_components) THEN
         ALLOCATE (pdos_array(nsoset(maxlgto), nkind, nmo + nvirt))
      ELSE
         ALLOCATE (pdos_array(0:maxlgto, nkind, nmo + nvirt))
      END IF
      IF (do_virt) THEN
         ALLOCATE (eigenvalues(nmo + nvirt))
         eigenvalues(1:nmo) = mo_set%eigenvalues(1:nmo)
         eigenvalues(nmo + 1:nmo + nvirt) = evals_virt(1:nvirt)
         ALLOCATE (occupation_numbers(nmo + nvirt))
         occupation_numbers(:) = 0.0_dp
         occupation_numbers(1:nmo) = mo_set%occupation_numbers(1:nmo)
      ELSE
         eigenvalues => mo_set%eigenvalues
         occupation_numbers => mo_set%occupation_numbers
      END IF

      hoco = -HUGE(0.0_dp)
      fractional_occupation = .FALSE.
      DO imo = 1, nmo + nvirt
         IF (occupation_numbers(imo) > 1.0e-10_dp) hoco = MAX(hoco, eigenvalues(imo))
         IF (ABS(occupation_numbers(imo) - REAL(NINT(occupation_numbers(imo)), KIND=dp)) > &
             1.0e-8_dp) fractional_occupation = .TRUE.
      END DO
      IF (hoco < -0.5_dp*HUGE(0.0_dp)) hoco = e_fermi
      IF (PRESENT(ispin) .AND. dft_control%nspins == 2) THEN
         e_fermi_ref(:) = 0.0_dp
         hoco_ref(:) = -HUGE(0.0_dp)
         CALL get_qs_env(qs_env=qs_env, mos=mos_ref)
         IF (ASSOCIATED(mos_ref)) THEN
            DO ispin_ref = 1, dft_control%nspins
               CALL get_mo_set(mo_set=mos_ref(ispin_ref), nmo=nmo_ref, mu=e_fermi_ref(ispin_ref))
               eval_ref => mos_ref(ispin_ref)%eigenvalues
               occ_ref => mos_ref(ispin_ref)%occupation_numbers
               DO imo_ref = 1, nmo_ref
                  IF (occ_ref(imo_ref) > 1.0e-10_dp) THEN
                     hoco_ref(ispin_ref) = MAX(hoco_ref(ispin_ref), eval_ref(imo_ref))
                  END IF
                  IF (ABS(occ_ref(imo_ref) - REAL(NINT(occ_ref(imo_ref)), KIND=dp)) > &
                      1.0e-8_dp) fractional_occupation = .TRUE.
               END DO
               IF (hoco_ref(ispin_ref) < -0.5_dp*HUGE(0.0_dp)) hoco_ref(ispin_ref) = e_fermi_ref(ispin_ref)
            END DO
         ELSE
            e_fermi_ref(:) = e_fermi
            hoco_ref(:) = hoco
         END IF
      ELSE
         e_fermi_ref(:) = e_fermi
         hoco_ref(:) = hoco
      END IF
      resolved_energy_zero = dos_resolve_energy_zero(energy_zero, dft_control%smear, fractional_occupation)
      SELECT CASE (resolved_energy_zero)
      CASE (dos_energy_zero_absolute)
         energy_ref = 0.0_dp
      CASE (dos_energy_zero_hoco)
         energy_ref = MAXVAL(hoco_ref(1:dft_control%nspins))
      CASE DEFAULT
         energy_ref = MAXVAL(e_fermi_ref(1:dft_control%nspins))
      END SELECT
      energy_factor = dos_energy_scale(energy_unit)
      energy_label = dos_energy_label(energy_unit)
      zero_label = dos_energy_zero_label(resolved_energy_zero)
      IF (energy_zero == dos_energy_zero_auto) zero_label = "AUTO -> "//TRIM(zero_label)
      ev_factor = dos_energy_scale(dos_energy_unit_ev)

      pdos_array = 0.0_dp
      nao = mo_set%nao
      ALLOCATE (vecbuffer(1, nao))
      vecbuffer = 0.0_dp
      ALLOCATE (firstrow(natom))
      firstrow = 0

      !Adjust energy range for r_ldos
      DO ildos = 1, n_r_ldos
         IF (eigenvalues(1) > r_ldos_p(ildos)%ldos%eval_range(1)) THEN
            r_ldos_p(ildos)%ldos%eval_range(1) = eigenvalues(1)
         END IF
         IF (eigenvalues(nmo + nvirt) < r_ldos_p(ildos)%ldos%eval_range(2)) THEN
            r_ldos_p(ildos)%ldos%eval_range(2) = eigenvalues(nmo + nvirt)
         END IF
      END DO

      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT='(/,(T15,A))') &
         "---- PDOS: start iteration on the KS states --- "

      DO imo = 1, nmo + nvirt

         IF (output_unit > 0 .AND. MOD(imo, out_each) == 0) WRITE (UNIT=output_unit, FMT='((T20,A,I10))') &
            " KS state index : ", imo
         ! Extract the eigenvector from the distributed full matrix
         IF (imo > nmo) THEN
            CALL cp_fm_get_submatrix(matrix_work, vecbuffer, 1, imo - nmo, &
                                     nao, 1, transpose=.TRUE.)
         ELSE
            CALL cp_fm_get_submatrix(matrix_shalfc, vecbuffer, 1, imo, &
                                     nao, 1, transpose=.TRUE.)
         END IF

         ! Calculate the pdos for all the kinds
         irow = 1
         DO iatom = 1, natom
            firstrow(iatom) = irow
            NULLIFY (orb_basis_set)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                   nset=nset, &
                                   nshell=nshell, &
                                   l=l, maxl=maxl)
            IF (separate_components) THEN
               isgf = 1
               DO iset = 1, nset
                  DO ishell = 1, nshell(iset)
                     lshell = l(ishell, iset)
                     DO iso = 1, nso(lshell)
                        lcomponent = nsoset(lshell - 1) + iso
                        pdos_array(lcomponent, ikind, imo) = &
                           pdos_array(lcomponent, ikind, imo) + &
                           vecbuffer(1, irow)*vecbuffer(1, irow)
                        irow = irow + 1
                     END DO ! iso
                  END DO ! ishell
               END DO ! iset
            ELSE
               isgf = 1
               DO iset = 1, nset
                  DO ishell = 1, nshell(iset)
                     lshell = l(ishell, iset)
                     DO iso = 1, nso(lshell)
                        pdos_array(lshell, ikind, imo) = &
                           pdos_array(lshell, ikind, imo) + &
                           vecbuffer(1, irow)*vecbuffer(1, irow)
                        irow = irow + 1
                     END DO ! iso
                  END DO ! ishell
               END DO ! iset
            END IF
         END DO ! iatom

         ! Calculate the pdos for all the lists
         DO ildos = 1, nldos
            DO il = 1, ldos_p(ildos)%ldos%nlist
               iatom = ldos_p(ildos)%ldos%list_index(il)

               irow = firstrow(iatom)
               NULLIFY (orb_basis_set)
               CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)

               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      nset=nset, &
                                      nshell=nshell, &
                                      l=l, maxl=maxl)
               ldos_p(ildos)%ldos%maxl = MAX(ldos_p(ildos)%ldos%maxl, maxl)
               IF (ldos_p(ildos)%ldos%separate_components) THEN
                  isgf = 1
                  DO iset = 1, nset
                     DO ishell = 1, nshell(iset)
                        lshell = l(ishell, iset)
                        DO iso = 1, nso(lshell)
                           lcomponent = nsoset(lshell - 1) + iso
                           ldos_p(ildos)%ldos%pdos_array(lcomponent, imo) = &
                              ldos_p(ildos)%ldos%pdos_array(lcomponent, imo) + &
                              vecbuffer(1, irow)*vecbuffer(1, irow)
                           irow = irow + 1
                        END DO ! iso
                     END DO ! ishell
                  END DO ! iset
               ELSE
                  isgf = 1
                  DO iset = 1, nset
                     DO ishell = 1, nshell(iset)
                        lshell = l(ishell, iset)
                        DO iso = 1, nso(lshell)
                           ldos_p(ildos)%ldos%pdos_array(lshell, imo) = &
                              ldos_p(ildos)%ldos%pdos_array(lshell, imo) + &
                              vecbuffer(1, irow)*vecbuffer(1, irow)
                           irow = irow + 1
                        END DO ! iso
                     END DO ! ishell
                  END DO ! iset
               END IF
            END DO !il
         END DO !ildos

         ! Calculate the DOS projected in a given volume in real space
         DO ildos = 1, n_r_ldos
            IF (r_ldos_p(ildos)%ldos%eval_range(1) <= eigenvalues(imo) .AND. &
                r_ldos_p(ildos)%ldos%eval_range(2) >= eigenvalues(imo)) THEN

               IF (imo > nmo) THEN
                  CALL calculate_wavefunction(mo_virt, imo - nmo, &
                                              wf_r, wf_g, atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                              pw_env)
               ELSE
                  CALL calculate_wavefunction(mo_coeff, imo, &
                                              wf_r, wf_g, atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                              pw_env)
               END IF
               r_ldos_p(ildos)%ldos%pdos_array(imo) = 0.0_dp
               DO il = 1, r_ldos_p(ildos)%ldos%npoints
                  j = j + 1
                  jx = r_ldos_p(ildos)%ldos%index_grid_local(1, il)
                  jy = r_ldos_p(ildos)%ldos%index_grid_local(2, il)
                  jz = r_ldos_p(ildos)%ldos%index_grid_local(3, il)
                  r_ldos_p(ildos)%ldos%pdos_array(imo) = r_ldos_p(ildos)%ldos%pdos_array(imo) + &
                                                         wf_r%array(jx, jy, jz)*wf_r%array(jx, jy, jz)
               END DO
               r_ldos_p(ildos)%ldos%pdos_array(imo) = r_ldos_p(ildos)%ldos%pdos_array(imo)*dvol
            END IF
         END DO
      END DO ! imo

      CALL cp_fm_release(matrix_shalfc)
      DEALLOCATE (vecbuffer)

      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%APPEND", l_val=append)
      IF (append .AND. iterstep > 1) THEN
         my_pos = "APPEND"
      ELSE
         my_pos = "REWIND"
      END IF
      my_act = "WRITE"
      IF (write_pdos_file) THEN
      DO ikind = 1, nkind

         NULLIFY (orb_basis_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), name=kind_name)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_set, maxl=maxl)

         ! basis none has no associated maxl, and no pdos
         IF (maxl < 0) CYCLE

         IF (PRESENT(ispin)) THEN
            IF (PRESENT(xas_mittle)) THEN
               my_mittle = TRIM(xas_mittle)//TRIM(spin(ispin))//"_k"//TRIM(ADJUSTL(cp_to_string(ikind)))
            ELSE
               my_mittle = TRIM(spin(ispin))//"_k"//TRIM(ADJUSTL(cp_to_string(ikind)))
            END IF
            my_spin = ispin
         ELSE
            my_mittle = "k"//TRIM(ADJUSTL(cp_to_string(ikind)))
            my_spin = 1
         END IF

         IF (write_pdos_file .AND. do_curve) THEN
            IF (separate_components) THEN
               CALL write_broadened_pdos(logger, dft_section, TRIM(my_mittle)//"_curve", my_pos, my_act, iterstep, &
                                         e_fermi, hoco, energy_ref, TRIM(zero_label), &
                                         "Projected DOS for atomic kind "//TRIM(kind_name), &
                                         maxl, .TRUE., pdos_array(1:nsoset(maxl), ikind, :), &
                                         eigenvalues, nmo + nvirt, de, broaden_type, broaden_width, &
                                         voigt_mixing, ndigits, pdos_print_key=TRIM(my_print_key))
            ELSE
               CALL write_broadened_pdos(logger, dft_section, TRIM(my_mittle)//"_curve", my_pos, my_act, iterstep, &
                                         e_fermi, hoco, energy_ref, TRIM(zero_label), &
                                         "Projected DOS for atomic kind "//TRIM(kind_name), &
                                         maxl, .FALSE., pdos_array(0:maxl, ikind, :), &
                                         eigenvalues, nmo + nvirt, de, broaden_type, broaden_width, &
                                         voigt_mixing, ndigits, pdos_print_key=TRIM(my_print_key))
            END IF
         END IF

         IF (write_pdos_file) THEN
            iw = cp_print_key_unit_nr(logger, dft_section, TRIM(my_print_key), &
                                      extension=".pdos", file_position=my_pos, file_action=my_act, &
                                      file_form="FORMATTED", middle_name=TRIM(my_mittle))
            IF (iw > 0) THEN

               fmtstr1 = "(I8,2X,2F16.6,  (2X,F16.8))"
               fmtstr2 = "(A42,  (10X,A8))"
               IF (separate_components) THEN
                  WRITE (UNIT=fmtstr1(15:16), FMT="(I2)") nsoset(maxl) + 1
                  WRITE (UNIT=fmtstr2(6:7), FMT="(I2)") nsoset(maxl) + 1
               ELSE
                  WRITE (UNIT=fmtstr1(15:16), FMT="(I2)") maxl + 2
                  WRITE (UNIT=fmtstr2(6:7), FMT="(I2)") maxl + 2
               END IF

               WRITE (UNIT=iw, FMT="(A,I0)") &
                  "# Projected DOS for atomic kind "//TRIM(kind_name)//" at iteration step i = ", iterstep
               WRITE (UNIT=iw, FMT="(A,F12.6,A,F12.6,A)") &
                  "# E(Fermi) = ", e_fermi, " a.u. = ", e_fermi*ev_factor, " eV"
               WRITE (UNIT=iw, FMT="(A,F12.6,A,F12.6,A)") &
                  "# E(HOCO)  = ", hoco, " a.u. = ", hoco*ev_factor, " eV"
               IF (separate_components) THEN
                  ALLOCATE (tmp_str(0:0, 0:maxl, -maxl:maxl))
                  tmp_str = ""
                  DO j = 0, maxl
                     DO i = -j, j
                        tmp_str(0, j, i) = sgf_symbol(0, j, i)
                     END DO
                  END DO

                  WRITE (UNIT=iw, FMT=fmtstr2) &
                     "#          MO Energy[a.u.]      Occupation", "Total", &
                     ((TRIM(tmp_str(0, il, im)), im=-il, il), il=0, maxl)
                  DO imo = 1, nmo + nvirt
                     WRITE (UNIT=iw, FMT=fmtstr1) imo, eigenvalues(imo), &
                        occupation_numbers(imo), SUM(pdos_array(1:nsoset(maxl), ikind, imo)), &
                        (pdos_array(lshell, ikind, imo), lshell=1, nsoset(maxl))
                  END DO
                  DEALLOCATE (tmp_str)
               ELSE
                  WRITE (UNIT=iw, FMT=fmtstr2) &
                     "#          MO Energy[a.u.]      Occupation", "Total", &
                     (TRIM(l_sym(il)), il=0, maxl)
                  DO imo = 1, nmo + nvirt
                     WRITE (UNIT=iw, FMT=fmtstr1) imo, eigenvalues(imo), &
                        occupation_numbers(imo), SUM(pdos_array(0:maxl, ikind, imo)), &
                        (pdos_array(lshell, ikind, imo), lshell=0, maxl)
                  END DO
               END IF
            END IF
            CALL cp_print_key_finished_output(iw, logger, dft_section, &
                                              TRIM(my_print_key))
         END IF

      END DO ! ikind
      END IF

      ! write the pdos for the lists, each ona different file,
      ! the filenames are indexed with the list number
      DO ildos = 1, nldos
         ! basis none has no associated maxl, and no pdos
         IF (ldos_p(ildos)%ldos%maxl > 0) THEN

            IF (PRESENT(ispin)) THEN
               IF (PRESENT(xas_mittle)) THEN
                  my_mittle = TRIM(xas_mittle)//TRIM(spin(ispin))//"_list"//TRIM(ldos_index(ildos))
               ELSE
                  my_mittle = TRIM(spin(ispin))//"_list"//TRIM(ldos_index(ildos))
               END IF
               my_spin = ispin
            ELSE
               my_mittle = "list"//TRIM(ldos_index(ildos))
               my_spin = 1
            END IF

            IF (do_curve) THEN
               IF (ldos_p(ildos)%ldos%separate_components) THEN
                  CALL write_broadened_pdos(logger, dft_section, TRIM(my_mittle)//"_curve", my_pos, my_act, iterstep, &
                                            e_fermi, hoco, energy_ref, TRIM(zero_label), &
                                            "Projected DOS for atom list "//TRIM(ldos_index(ildos)), &
                                            ldos_p(ildos)%ldos%maxl, .TRUE., &
                                            ldos_p(ildos)%ldos%pdos_array(1:nsoset(ldos_p(ildos)%ldos%maxl), :), &
                                            eigenvalues, nmo + nvirt, de, broaden_type, broaden_width, &
                                            voigt_mixing, ndigits, pdos_print_key=TRIM(my_print_key))
               ELSE
                  CALL write_broadened_pdos(logger, dft_section, TRIM(my_mittle)//"_curve", my_pos, my_act, iterstep, &
                                            e_fermi, hoco, energy_ref, TRIM(zero_label), &
                                            "Projected DOS for atom list "//TRIM(ldos_index(ildos)), &
                                            ldos_p(ildos)%ldos%maxl, .FALSE., &
                                            ldos_p(ildos)%ldos%pdos_array(0:ldos_p(ildos)%ldos%maxl, :), &
                                            eigenvalues, nmo + nvirt, de, broaden_type, broaden_width, &
                                            voigt_mixing, ndigits, pdos_print_key=TRIM(my_print_key))
               END IF
            END IF

            iw = cp_print_key_unit_nr(logger, dft_section, TRIM(my_print_key), &
                                      extension=".pdos", file_position=my_pos, file_action=my_act, &
                                      file_form="FORMATTED", middle_name=TRIM(my_mittle))
            IF (iw > 0) THEN

               fmtstr1 = "(I8,2X,2F16.6,  (2X,F16.8))"
               fmtstr2 = "(A42,  (10X,A8))"
               IF (ldos_p(ildos)%ldos%separate_components) THEN
                  WRITE (UNIT=fmtstr1(15:16), FMT="(I2)") nsoset(ldos_p(ildos)%ldos%maxl) + 1
                  WRITE (UNIT=fmtstr2(6:7), FMT="(I2)") nsoset(ldos_p(ildos)%ldos%maxl) + 1
               ELSE
                  WRITE (UNIT=fmtstr1(15:16), FMT="(I2)") ldos_p(ildos)%ldos%maxl + 2
                  WRITE (UNIT=fmtstr2(6:7), FMT="(I2)") ldos_p(ildos)%ldos%maxl + 2
               END IF

               WRITE (UNIT=iw, FMT="(A,I0,A,I0,A,I0)") &
                  "# Projected DOS for list ", ildos, " of ", ldos_p(ildos)%ldos%nlist, &
                  " atoms, at iteration step i = ", iterstep
               WRITE (UNIT=iw, FMT="(A,F12.6,A,F12.6,A)") &
                  "# E(Fermi) = ", e_fermi, " a.u. = ", e_fermi*ev_factor, " eV"
               WRITE (UNIT=iw, FMT="(A,F12.6,A,F12.6,A)") &
                  "# E(HOCO)  = ", hoco, " a.u. = ", hoco*ev_factor, " eV"
               IF (ldos_p(ildos)%ldos%separate_components) THEN
                  ALLOCATE (tmp_str(0:0, 0:ldos_p(ildos)%ldos%maxl, -ldos_p(ildos)%ldos%maxl:ldos_p(ildos)%ldos%maxl))
                  tmp_str = ""
                  DO j = 0, ldos_p(ildos)%ldos%maxl
                     DO i = -j, j
                        tmp_str(0, j, i) = sgf_symbol(0, j, i)
                     END DO
                  END DO

                  WRITE (UNIT=iw, FMT=fmtstr2) &
                     "#          MO Energy[a.u.]      Occupation", "Total", &
                     ((TRIM(tmp_str(0, il, im)), im=-il, il), il=0, ldos_p(ildos)%ldos%maxl)
                  DO imo = 1, nmo + nvirt
                     WRITE (UNIT=iw, FMT=fmtstr1) imo, eigenvalues(imo), &
                        occupation_numbers(imo), &
                        SUM(ldos_p(ildos)%ldos%pdos_array(1:nsoset(ldos_p(ildos)%ldos%maxl), imo)), &
                        (ldos_p(ildos)%ldos%pdos_array(lshell, imo), &
                         lshell=1, nsoset(ldos_p(ildos)%ldos%maxl))
                  END DO
                  DEALLOCATE (tmp_str)
               ELSE
                  WRITE (UNIT=iw, FMT=fmtstr2) &
                     "#          MO Energy[a.u.]      Occupation", "Total", &
                     (TRIM(l_sym(il)), il=0, ldos_p(ildos)%ldos%maxl)
                  DO imo = 1, nmo + nvirt
                     WRITE (UNIT=iw, FMT=fmtstr1) imo, eigenvalues(imo), &
                        occupation_numbers(imo), &
                        SUM(ldos_p(ildos)%ldos%pdos_array(0:ldos_p(ildos)%ldos%maxl, imo)), &
                        (ldos_p(ildos)%ldos%pdos_array(lshell, imo), lshell=0, ldos_p(ildos)%ldos%maxl)
                  END DO
               END IF
            END IF
            CALL cp_print_key_finished_output(iw, logger, dft_section, &
                                              TRIM(my_print_key))
         END IF ! maxl>0
      END DO ! ildos

      ! write the pdos for the lists, each ona different file,
      ! the filenames are indexed with the list number
      DO ildos = 1, n_r_ldos

         npoints = r_ldos_p(ildos)%ldos%npoints
         CALL para_env%sum(npoints)
         CALL para_env%sum(np_tot)
         CALL para_env%sum(r_ldos_p(ildos)%ldos%pdos_array)
         IF (PRESENT(ispin)) THEN
            IF (PRESENT(xas_mittle)) THEN
               my_mittle = TRIM(xas_mittle)//TRIM(spin(ispin))//"_r_list"//TRIM(r_ldos_index(ildos))
            ELSE
               my_mittle = TRIM(spin(ispin))//"_r_list"//TRIM(r_ldos_index(ildos))
            END IF
            my_spin = ispin
         ELSE
            my_mittle = "r_list"//TRIM(r_ldos_index(ildos))
            my_spin = 1
         END IF

         iw = cp_print_key_unit_nr(logger, dft_section, TRIM(my_print_key), &
                                   extension=".pdos", file_position=my_pos, file_action=my_act, &
                                   file_form="FORMATTED", middle_name=TRIM(my_mittle))
         IF (iw > 0) THEN
            fmtstr1 = "(I8,2X,2F16.6,  (2X,F16.8))"
            fmtstr2 = "(A42,  (10X,A8))"

            WRITE (UNIT=iw, FMT="(A,I0,A,F12.6,F12.6,A)") &
               "# Projected DOS in real space, using ", npoints, &
               " points of the grid, and eval in the range", r_ldos_p(ildos)%ldos%eval_range(1:2), " Hartree"
            WRITE (UNIT=iw, FMT="(A,F12.6,A,F12.6,A)") &
               "# E(Fermi) = ", e_fermi, " a.u. = ", e_fermi*ev_factor, " eV"
            WRITE (UNIT=iw, FMT="(A,F12.6,A,F12.6,A)") &
               "# E(HOCO)  = ", hoco, " a.u. = ", hoco*ev_factor, " eV"
            WRITE (UNIT=iw, FMT="(A)") &
               "#          MO Energy[a.u.]      Occupation      LDOS"
            DO imo = 1, nmo + nvirt
               IF (r_ldos_p(ildos)%ldos%eval_range(1) <= eigenvalues(imo) .AND. &
                   r_ldos_p(ildos)%ldos%eval_range(2) >= eigenvalues(imo)) THEN
                  WRITE (UNIT=iw, FMT="(I8,2X,2F16.6,E20.10,E20.10)") imo, &
                     eigenvalues(imo), occupation_numbers(imo), &
                     r_ldos_p(ildos)%ldos%pdos_array(imo), r_ldos_p(ildos)%ldos%pdos_array(imo)*np_tot
               END IF
            END DO

         END IF
         CALL cp_print_key_finished_output(iw, logger, dft_section, &
                                           TRIM(my_print_key))
      END DO

      ! deallocate local variables
      DEALLOCATE (pdos_array)
      DEALLOCATE (firstrow)
      IF (do_ldos) THEN
         DO ildos = 1, nldos
            DEALLOCATE (ldos_p(ildos)%ldos%pdos_array)
            DEALLOCATE (ldos_p(ildos)%ldos%list_index)
            DEALLOCATE (ldos_p(ildos)%ldos)
         END DO
         DEALLOCATE (ldos_p)
         DEALLOCATE (ldos_index)
      END IF
      IF (do_r_ldos) THEN
         DO ildos = 1, n_r_ldos
            DEALLOCATE (r_ldos_p(ildos)%ldos%index_grid_local)
            DEALLOCATE (r_ldos_p(ildos)%ldos%pdos_array)
            DEALLOCATE (r_ldos_p(ildos)%ldos%list_index)
            IF (.NOT. read_r(1, ildos)) THEN
               DEALLOCATE (r_ldos_p(ildos)%ldos%x_range)
            END IF
            IF (.NOT. read_r(2, ildos)) THEN
               DEALLOCATE (r_ldos_p(ildos)%ldos%y_range)
            END IF
            IF (.NOT. read_r(3, ildos)) THEN
               DEALLOCATE (r_ldos_p(ildos)%ldos%z_range)
            END IF
            IF (.NOT. read_r(4, ildos)) THEN
               DEALLOCATE (r_ldos_p(ildos)%ldos%eval_range)
            END IF
            DEALLOCATE (r_ldos_p(ildos)%ldos)
         END DO
         DEALLOCATE (read_r)
         DEALLOCATE (r_ldos_p)
         DEALLOCATE (r_ldos_index)
         CALL auxbas_pw_pool%give_back_pw(wf_r)
         CALL auxbas_pw_pool%give_back_pw(wf_g)
      END IF
      IF (do_virt) THEN
         CALL cp_fm_release(matrix_work)
         DEALLOCATE (eigenvalues)
         DEALLOCATE (occupation_numbers)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_projected_dos

! **************************************************************************************************
!> \brief Compute and write broadened projected density of states for k-point calculations.
!> \param qs_env ...
!> \param dft_section ...
!> \param pdos_print_key ...
!> \param write_pdos ...
!> \param write_pdos_curve ...
! **************************************************************************************************
   SUBROUTINE calculate_projected_dos_kp(qs_env, dft_section, pdos_print_key, write_pdos, write_pdos_curve)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: dft_section
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: pdos_print_key
      LOGICAL, INTENT(IN), OPTIONAL                      :: write_pdos, write_pdos_curve

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_projected_dos_kp'

      CHARACTER(LEN=32)                                  :: zero_label
      CHARACTER(LEN=default_string_length)               :: kind_name, my_act, my_mittle, my_pos, &
                                                            my_print_key, spin(2)
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: zvecbuffer
      INTEGER :: broaden_type, energy_zero, fractional_occupation_int, handle, icomp, ik, ikind, &
         imo, ispin, iterstep, maxl, maxlgto, n_r_ldos, nao, ncomp, ndigits, nhist, nkind, nldos, &
         nmo_kp, nspins, output_unit, resolved_energy_zero
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_comp, ao_kind, ao_l, kind_maxl
      LOGICAL                                            :: append, fractional_occupation, &
                                                            separate_components, should_output, &
                                                            write_curve, write_pdos_file
      REAL(KIND=dp)                                      :: broaden_cutoff, broaden_width, de, e1, &
                                                            e2, e_fermi(2), emax, emin, &
                                                            energy_ref(2), hoco(2), voigt_mixing, &
                                                            wkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ao_weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: proj_weight, vecbuffer
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: pdos_curve
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_cfm_type)                                  :: cshalfc
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: shalfc
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: curve_section, ldos_section

      NULLIFY (logger, kpoints, dft_control, para_env, atomic_kind_set, qs_kind_set, particle_set)
      NULLIFY (matrix_s_kp, scf_env)
      NULLIFY (kp, mo_set, eigenvalues, fm_struct_tmp, orb_basis_set, curve_section, ldos_section)
      logger => cp_get_default_logger()
      my_print_key = "PRINT%PDOS"
      IF (PRESENT(pdos_print_key)) my_print_key = TRIM(pdos_print_key)
      write_pdos_file = .TRUE.
      IF (PRESENT(write_pdos)) write_pdos_file = write_pdos
      curve_section => section_vals_get_subs_vals(dft_section, TRIM(my_print_key)//"%CURVE")
      CALL section_vals_get(curve_section, explicit=write_curve)
      IF (PRESENT(write_pdos_curve)) write_curve = write_pdos_curve
      should_output = BTEST(cp_print_key_should_output(logger%iter_info, dft_section, &
                                                       TRIM(my_print_key)), cp_p_file)
      output_unit = cp_logger_get_default_io_unit(logger)
      IF ((.NOT. should_output)) RETURN

      CALL timeset(routineN, handle)
      iterstep = logger%iter_info%iteration(logger%iter_info%n_rlevel)

      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT='(/,(T3,A,T61,I10))') &
         " Calculate k-point PDOS at iteration step ", iterstep

      CALL get_qs_env(qs_env=qs_env, &
                      kpoints=kpoints, &
                      dft_control=dft_control, &
                      matrix_s_kp=matrix_s_kp, &
                      scf_env=scf_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set)
      para_env => kpoints%para_env_inter_kp
      nspins = dft_control%nspins
      nkind = SIZE(atomic_kind_set)
      CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto)
      IF (.NOT. ASSOCIATED(kpoints%kp_env)) THEN
         CPWARN("No local k points available for k-point PDOS")
         CALL timestop(handle)
         RETURN
      END IF
      IF (SIZE(kpoints%kp_env) == 0) THEN
         CPWARN("No local k points available for k-point PDOS")
         CALL timestop(handle)
         RETURN
      END IF

      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%DELTA_E", r_val=de)
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%APPEND", l_val=append)
      IF (TRIM(my_print_key) == "PRINT%DOS") THEN
         CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%PDOS%COMPONENTS", l_val=separate_components)
      ELSE
         CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%COMPONENTS", l_val=separate_components)
      END IF
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%NDIGITS", i_val=ndigits)
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%CURVE%BROADEN%TYPE", i_val=broaden_type)
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%CURVE%BROADEN%WIDTH", r_val=broaden_width)
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%CURVE%BROADEN%VOIGT_MIXING", r_val=voigt_mixing)
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%CURVE%ENERGY_ZERO", i_val=energy_zero)
      ndigits = MIN(MAX(ndigits, 1), 10)
      IF (write_curve .AND. de <= 0.0_dp) THEN
         CPWARN("Broadened k-point PDOS output requires DELTA_E > 0 and will be skipped")
         CALL timestop(handle)
         RETURN
      END IF
      IF (write_curve) de = MAX(de, 0.00001_dp)

      ldos_section => section_vals_get_subs_vals(dft_section, TRIM(my_print_key)//"%LDOS")
      CALL section_vals_get(ldos_section, n_repetition=nldos)
      ldos_section => section_vals_get_subs_vals(dft_section, TRIM(my_print_key)//"%R_LDOS")
      CALL section_vals_get(ldos_section, n_repetition=n_r_ldos)
      IF (nldos > 0 .OR. n_r_ldos > 0) THEN
         CPWARN("LDOS/R_LDOS are not implemented for k-point PDOS and will be ignored")
      END IF
      IF (write_pdos_file) THEN
         CPWARN("State-resolved k-point PDOS output is not implemented yet")
      END IF
      IF (.NOT. write_curve) THEN
         CALL timestop(handle)
         RETURN
      END IF
      IF (broaden_width <= 0.0_dp) THEN
         CPWARN("Broadened k-point PDOS output requires a finite WIDTH and will be skipped")
         CALL timestop(handle)
         RETURN
      END IF

      IF (separate_components) THEN
         ncomp = nsoset(maxlgto)
      ELSE
         ncomp = maxlgto + 1
      END IF
      ALLOCATE (kind_maxl(nkind))
      kind_maxl = -1
      DO ikind = 1, nkind
         NULLIFY (orb_basis_set)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_set, maxl=maxl)
         kind_maxl(ikind) = maxl
      END DO

      emin = HUGE(0.0_dp)
      emax = -HUGE(0.0_dp)
      e_fermi(:) = 0.0_dp
      hoco(:) = -HUGE(0.0_dp)
      fractional_occupation = .FALSE.
      IF (kpoints%nkp /= 0) THEN
         DO ik = 1, SIZE(kpoints%kp_env)
            kp => kpoints%kp_env(ik)%kpoint_env
            DO ispin = 1, nspins
               mo_set => kp%mos(1, ispin)
               CALL get_mo_set(mo_set=mo_set, nmo=nmo_kp, mu=e_fermi(ispin))
               eigenvalues => mo_set%eigenvalues
               occupation_numbers => mo_set%occupation_numbers
               DO imo = 1, nmo_kp
                  IF (occupation_numbers(imo) > 1.0e-10_dp) hoco(ispin) = MAX(hoco(ispin), eigenvalues(imo))
                  IF (ABS(occupation_numbers(imo) - REAL(NINT(occupation_numbers(imo)), KIND=dp)) > &
                      1.0e-8_dp) fractional_occupation = .TRUE.
               END DO
               e1 = MINVAL(eigenvalues(1:nmo_kp))
               e2 = MAXVAL(eigenvalues(1:nmo_kp))
               emin = MIN(emin, e1)
               emax = MAX(emax, e2)
            END DO
         END DO
      END IF
      CALL para_env%min(emin)
      CALL para_env%max(emax)
      CALL para_env%max(e_fermi)
      CALL para_env%max(hoco)
      fractional_occupation_int = MERGE(1, 0, fractional_occupation)
      CALL para_env%max(fractional_occupation_int)
      fractional_occupation = (fractional_occupation_int /= 0)
      DO ispin = 1, nspins
         IF (hoco(ispin) < -0.5_dp*HUGE(0.0_dp)) hoco(ispin) = e_fermi(ispin)
      END DO
      resolved_energy_zero = dos_resolve_energy_zero(energy_zero, dft_control%smear, fractional_occupation)
      SELECT CASE (resolved_energy_zero)
      CASE (dos_energy_zero_absolute)
         energy_ref(:) = 0.0_dp
      CASE (dos_energy_zero_hoco)
         energy_ref(:) = MAXVAL(hoco(1:nspins))
      CASE DEFAULT
         energy_ref(:) = MAXVAL(e_fermi(1:nspins))
      END SELECT
      zero_label = dos_energy_zero_label(resolved_energy_zero)
      IF (energy_zero == dos_energy_zero_auto) zero_label = "AUTO -> "//TRIM(zero_label)
      broaden_cutoff = broadening_cutoff(broaden_type, broaden_width)
      emin = emin - broaden_cutoff
      emax = emax + broaden_cutoff
      nhist = NINT((emax - emin)/de) + 1
      ALLOCATE (pdos_curve(nhist, ncomp, nkind, nspins))
      pdos_curve = 0.0_dp

      ! Ensure that S(k)^1/2 is available for the Lowdin projection.
      ! This is normally only constructed for Lowdin population/DFT+U paths.
      CALL diag_kp_smat(matrix_s_kp, kpoints, scf_env%scf_work1)

      ! Use the first local k point to construct the AO -> kind/l/component map.
      kp => kpoints%kp_env(1)%kpoint_env
      mo_set => kp%mos(1, 1)
      CALL get_mo_set(mo_set=mo_set, nao=nao)
      CALL build_pdos_ao_map(qs_kind_set, particle_set, nao, ao_kind, ao_l, ao_comp)
      ALLOCATE (ao_weight(nao), proj_weight(ncomp, nkind))
      ALLOCATE (vecbuffer(1, nao), zvecbuffer(1, nao))

      IF (kpoints%nkp /= 0) THEN
         DO ik = 1, SIZE(kpoints%kp_env)
            kp => kpoints%kp_env(ik)%kpoint_env
            wkp = kp%wkp
            DO ispin = 1, nspins
               mo_set => kp%mos(1, ispin)
               CALL get_mo_set(mo_set=mo_set, nao=nao, nmo=nmo_kp)
               eigenvalues => mo_set%eigenvalues
               CALL cp_fm_get_info(mo_set%mo_coeff, matrix_struct=fm_struct_tmp)
               IF (kpoints%use_real_wfn) THEN
                  CALL cp_fm_create(shalfc, fm_struct_tmp, name="shalfc")
                  CALL lowdin_kp_mo_coeff(kp, ispin, kpoints%use_real_wfn, shalfc=shalfc)
               ELSE
                  CALL cp_cfm_create(cshalfc, fm_struct_tmp, name="cshalfc")
                  CALL lowdin_kp_mo_coeff(kp, ispin, kpoints%use_real_wfn, cshalfc=cshalfc)
               END IF

               DO imo = 1, nmo_kp
                  IF (kpoints%use_real_wfn) THEN
                     CALL cp_fm_get_submatrix(shalfc, vecbuffer, 1, imo, nao, 1, transpose=.TRUE.)
                     ao_weight(:) = vecbuffer(1, 1:nao)**2
                  ELSE
                     CALL cp_cfm_get_submatrix(cshalfc, zvecbuffer, 1, imo, nao, 1, transpose=.TRUE.)
                     ao_weight(:) = REAL(CONJG(zvecbuffer(1, 1:nao))*zvecbuffer(1, 1:nao), KIND=dp)
                  END IF
                  proj_weight = 0.0_dp
                  CALL accumulate_pdos_weights(ao_weight, ao_kind, ao_l, ao_comp, &
                                               separate_components, proj_weight)
                  DO ikind = 1, nkind
                     IF (kind_maxl(ikind) < 0) CYCLE
                     IF (separate_components) THEN
                        DO icomp = 1, nsoset(kind_maxl(ikind))
                           CALL add_broadened_value(pdos_curve(:, icomp, ikind, ispin), &
                                                    emin, de, eigenvalues(imo), &
                                                    wkp*proj_weight(icomp, ikind), &
                                                    broaden_type, broaden_width, voigt_mixing)
                        END DO
                     ELSE
                        DO icomp = 1, kind_maxl(ikind) + 1
                           CALL add_broadened_value(pdos_curve(:, icomp, ikind, ispin), &
                                                    emin, de, eigenvalues(imo), &
                                                    wkp*proj_weight(icomp, ikind), &
                                                    broaden_type, broaden_width, voigt_mixing)
                        END DO
                     END IF
                  END DO
               END DO

               IF (kpoints%use_real_wfn) THEN
                  CALL cp_fm_release(shalfc)
               ELSE
                  CALL cp_cfm_release(cshalfc)
               END IF
            END DO
         END DO
      END IF
      CALL para_env%sum(pdos_curve)

      IF (append .AND. iterstep > 1) THEN
         my_pos = "APPEND"
      ELSE
         my_pos = "REWIND"
      END IF
      my_act = "WRITE"
      spin(1) = "ALPHA"
      spin(2) = "BETA"
      IF (write_pdos_file) THEN
      DO ikind = 1, nkind
         IF (kind_maxl(ikind) < 0) CYCLE
         CALL get_atomic_kind(atomic_kind_set(ikind), name=kind_name)
         DO ispin = 1, nspins
            IF (nspins == 2) THEN
               my_mittle = TRIM(spin(ispin))//"_k"//TRIM(ADJUSTL(cp_to_string(ikind)))
            ELSE
               my_mittle = "k"//TRIM(ADJUSTL(cp_to_string(ikind)))
            END IF
            CALL write_broadened_pdos_curve(logger, dft_section, TRIM(my_mittle)//"_curve", my_pos, my_act, &
                                            iterstep, e_fermi(ispin), hoco(ispin), energy_ref(ispin), &
                                            TRIM(zero_label), &
                                            "K-point projected DOS for atomic kind "//TRIM(kind_name), &
                                            kind_maxl(ikind), separate_components, &
                                            pdos_curve(:, :, ikind, ispin), emin, de, &
                                            broaden_type, broaden_width, voigt_mixing, ndigits, pdos_print_key=TRIM(my_print_key))
         END DO
      END DO
      END IF

      DEALLOCATE (ao_comp, ao_kind, ao_l, ao_weight, kind_maxl, pdos_curve, proj_weight, &
                  vecbuffer, zvecbuffer)

      CALL timestop(handle)

   END SUBROUTINE calculate_projected_dos_kp

! **************************************************************************************************
!> \brief Build AO mapping arrays for PDOS accumulation.
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param nao ...
!> \param ao_kind ...
!> \param ao_l ...
!> \param ao_comp ...
! **************************************************************************************************
   SUBROUTINE build_pdos_ao_map(qs_kind_set, particle_set, nao, ao_kind, ao_l, ao_comp)

      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: nao
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: ao_kind, ao_l, ao_comp

      INTEGER                                            :: iatom, ikind, irow, iset, ishell, iso, &
                                                            lshell, maxl, nset
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      ALLOCATE (ao_kind(nao), ao_l(nao), ao_comp(nao))
      irow = 0
      DO iatom = 1, SIZE(particle_set)
         NULLIFY (orb_basis_set)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                nset=nset, &
                                nshell=nshell, &
                                l=l, maxl=maxl)
         DO iset = 1, nset
            DO ishell = 1, nshell(iset)
               lshell = l(ishell, iset)
               DO iso = 1, nso(lshell)
                  irow = irow + 1
                  CPASSERT(irow <= nao)
                  ao_kind(irow) = ikind
                  ao_l(irow) = lshell
                  ao_comp(irow) = nsoset(lshell - 1) + iso
               END DO
            END DO
         END DO
      END DO
      CPASSERT(irow == nao)

   END SUBROUTINE build_pdos_ao_map

! **************************************************************************************************
!> \brief Accumulate AO weights into kind/l or kind/component projected weights.
!> \param ao_weight ...
!> \param ao_kind ...
!> \param ao_l ...
!> \param ao_comp ...
!> \param separate_components ...
!> \param proj_weight ...
! **************************************************************************************************
   SUBROUTINE accumulate_pdos_weights(ao_weight, ao_kind, ao_l, ao_comp, &
                                      separate_components, proj_weight)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: ao_weight
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ao_kind, ao_l, ao_comp
      LOGICAL, INTENT(IN)                                :: separate_components
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: proj_weight

      INTEGER                                            :: iao, icomp, ikind

      DO iao = 1, SIZE(ao_weight)
         ikind = ao_kind(iao)
         IF (separate_components) THEN
            icomp = ao_comp(iao)
         ELSE
            icomp = ao_l(iao) + 1
         END IF
         proj_weight(icomp, ikind) = proj_weight(icomp, ikind) + ao_weight(iao)
      END DO

   END SUBROUTINE accumulate_pdos_weights

! **************************************************************************************************
!> \brief Write a broadened k-point PDOS curve.
!> \param logger ...
!> \param dft_section ...
!> \param middle_name ...
!> \param file_position ...
!> \param file_action ...
!> \param iterstep ...
!> \param e_fermi ...
!> \param hoco ...
!> \param energy_ref ...
!> \param zero_label ...
!> \param title ...
!> \param maxl ...
!> \param separate_components ...
!> \param pdos_curve ...
!> \param emin ...
!> \param de ...
!> \param broaden_type ...
!> \param broaden_width ...
!> \param voigt_mixing ...
!> \param ndigits ...
!> \param pdos_print_key ...
! **************************************************************************************************
   SUBROUTINE write_broadened_pdos_curve(logger, dft_section, middle_name, file_position, &
                                         file_action, iterstep, e_fermi, hoco, energy_ref, zero_label, &
                                         title, maxl, &
                                         separate_components, pdos_curve, emin, de, &
                                         broaden_type, broaden_width, voigt_mixing, ndigits, pdos_print_key)

      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: dft_section
      CHARACTER(LEN=*), INTENT(IN)                       :: middle_name, file_position, file_action
      INTEGER, INTENT(IN)                                :: iterstep
      REAL(KIND=dp), INTENT(IN)                          :: e_fermi, hoco, energy_ref
      CHARACTER(LEN=*), INTENT(IN)                       :: zero_label, title
      INTEGER, INTENT(IN)                                :: maxl
      LOGICAL, INTENT(IN)                                :: separate_components
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: pdos_curve
      REAL(KIND=dp), INTENT(IN)                          :: emin, de
      INTEGER, INTENT(IN)                                :: broaden_type
      REAL(KIND=dp), INTENT(IN)                          :: broaden_width, voigt_mixing
      INTEGER, INTENT(IN)                                :: ndigits
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: pdos_print_key

      CHARACTER(LEN=16)                                  :: energy_label
      CHARACTER(LEN=20)                                  :: fmtstr_data
      CHARACTER(LEN=6), ALLOCATABLE, DIMENSION(:, :, :)  :: tmp_str
      CHARACTER(LEN=default_string_length)               :: my_print_key
      INTEGER                                            :: energy_unit, i, icomp, il, im, iw, &
                                                            ncomponents, nhist
      REAL(KIND=dp)                                      :: density_factor, energy_factor, &
                                                            ev_factor, eval

      my_print_key = "PRINT%PDOS"
      IF (PRESENT(pdos_print_key)) my_print_key = TRIM(pdos_print_key)

      nhist = SIZE(pdos_curve, 1)
      IF (separate_components) THEN
         ncomponents = nsoset(maxl)
      ELSE
         ncomponents = maxl + 1
      END IF
      CALL section_vals_val_get(dft_section, TRIM(my_print_key)//"%CURVE%ENERGY_UNIT", i_val=energy_unit)
      energy_factor = dos_energy_scale(energy_unit)
      density_factor = dos_density_scale(energy_unit)
      energy_label = dos_energy_label(energy_unit)
      ev_factor = dos_energy_scale(dos_energy_unit_ev)
      iw = cp_print_key_unit_nr(logger, dft_section, TRIM(my_print_key), &
                                extension=".pdos", file_position=file_position, file_action=file_action, &
                                file_form="FORMATTED", middle_name=TRIM(middle_name))
      IF (iw > 0) THEN
         WRITE (UNIT=iw, FMT="(A,I0)") "# "//TRIM(title)//" at iteration step i = ", iterstep
         WRITE (UNIT=iw, FMT="(A,F12.6,A,F12.6,A)") &
            "# E(Fermi) = ", e_fermi, " a.u. = ", e_fermi*ev_factor, " eV"
         WRITE (UNIT=iw, FMT="(A,F12.6,A,F12.6,A)") &
            "# E(HOCO)  = ", hoco, " a.u. = ", hoco*ev_factor, " eV"
         WRITE (UNIT=iw, FMT="(A,A)") "# Energy zero: ", TRIM(zero_label)
         CALL write_broadening_info(iw, broaden_type, broaden_width, voigt_mixing)
         WRITE (UNIT=iw, FMT="(A)", ADVANCE="NO") "# "//TRIM(energy_label)
         WRITE (UNIT=iw, FMT="(2X,A)", ADVANCE="NO") "total"
         IF (separate_components) THEN
            ALLOCATE (tmp_str(0:0, 0:maxl, -maxl:maxl))
            tmp_str = ""
            DO il = 0, maxl
               DO im = -il, il
                  tmp_str(0, il, im) = sgf_symbol(0, il, im)
                  WRITE (UNIT=iw, FMT="(2X,A)", ADVANCE="NO") TRIM(tmp_str(0, il, im))
               END DO
            END DO
            DEALLOCATE (tmp_str)
         ELSE
            DO il = 0, maxl
               WRITE (UNIT=iw, FMT="(2X,A)", ADVANCE="NO") TRIM(l_sym(il))
            END DO
         END IF
         WRITE (UNIT=iw, FMT="()")
         WRITE (UNIT=fmtstr_data, FMT="(A,I0,A)") "(2X,F20.", ndigits, ")"
         DO i = 1, nhist
            eval = (emin + (i - 1)*de - energy_ref)*energy_factor
            WRITE (UNIT=iw, FMT="(F15.8)", ADVANCE="NO") eval
            WRITE (UNIT=iw, FMT=fmtstr_data, ADVANCE="NO") SUM(pdos_curve(i, 1:ncomponents))*density_factor
            DO icomp = 1, ncomponents
               WRITE (UNIT=iw, FMT=fmtstr_data, ADVANCE="NO") pdos_curve(i, icomp)*density_factor
            END DO
            WRITE (UNIT=iw, FMT="()")
         END DO
      END IF
      CALL cp_print_key_finished_output(iw, logger, dft_section, TRIM(my_print_key))

   END SUBROUTINE write_broadened_pdos_curve

! **************************************************************************************************
!> \brief Write a broadened PDOS curve for a projected weight matrix.
!> \param logger ...
!> \param dft_section ...
!> \param middle_name ...
!> \param file_position ...
!> \param file_action ...
!> \param iterstep ...
!> \param e_fermi ...
!> \param hoco ...
!> \param energy_ref ...
!> \param zero_label ...
!> \param title ...
!> \param maxl ...
!> \param separate_components ...
!> \param weights ...
!> \param eigenvalues ...
!> \param nstates ...
!> \param de ...
!> \param broaden_type ...
!> \param broaden_width ...
!> \param voigt_mixing ...
!> \param ndigits ...
!> \param pdos_print_key ...
! **************************************************************************************************
   SUBROUTINE write_broadened_pdos(logger, dft_section, middle_name, file_position, file_action, &
                                   iterstep, e_fermi, hoco, energy_ref, zero_label, title, maxl, separate_components, weights, &
                                   eigenvalues, nstates, de, broaden_type, broaden_width, &
                                   voigt_mixing, ndigits, pdos_print_key)

      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: dft_section
      CHARACTER(LEN=*), INTENT(IN)                       :: middle_name, file_position, file_action
      INTEGER, INTENT(IN)                                :: iterstep
      REAL(KIND=dp), INTENT(IN)                          :: e_fermi, hoco, energy_ref
      CHARACTER(LEN=*), INTENT(IN)                       :: zero_label, title
      INTEGER, INTENT(IN)                                :: maxl
      LOGICAL, INTENT(IN)                                :: separate_components
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: weights
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues
      INTEGER, INTENT(IN)                                :: nstates
      REAL(KIND=dp), INTENT(IN)                          :: de
      INTEGER, INTENT(IN)                                :: broaden_type
      REAL(KIND=dp), INTENT(IN)                          :: broaden_width, voigt_mixing
      INTEGER, INTENT(IN)                                :: ndigits
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: pdos_print_key

      INTEGER                                            :: i, icomp, imo, ncomponents, nhist
      REAL(KIND=dp)                                      :: cutoff, emax, emin, eval, line_shape
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: pdos_curve

      IF (broaden_width <= 0.0_dp) RETURN

      ncomponents = SIZE(weights, 1)
      cutoff = broadening_cutoff(broaden_type, broaden_width)
      emin = MINVAL(eigenvalues(1:nstates)) - cutoff
      emax = MAXVAL(eigenvalues(1:nstates)) + cutoff
      nhist = NINT((emax - emin)/de) + 1
      ALLOCATE (pdos_curve(nhist, ncomponents))
      pdos_curve = 0.0_dp

      DO imo = 1, nstates
         DO i = MAX(1, FLOOR((eigenvalues(imo) - cutoff - emin)/de) + 1), &
            MIN(nhist, CEILING((eigenvalues(imo) + cutoff - emin)/de) + 1)
            eval = emin + (i - 1)*de
            line_shape = broadening_function(eval - eigenvalues(imo), broaden_type, broaden_width, &
                                             voigt_mixing)
            DO icomp = 1, ncomponents
               pdos_curve(i, icomp) = pdos_curve(i, icomp) + weights(icomp, imo)*line_shape
            END DO
         END DO
      END DO

      CALL write_broadened_pdos_curve(logger, dft_section, middle_name, file_position, file_action, &
                                      iterstep, e_fermi, hoco, energy_ref, zero_label, title, maxl, separate_components, &
                                      pdos_curve, emin, de, broaden_type, broaden_width, &
                                      voigt_mixing, ndigits, pdos_print_key=pdos_print_key)

      DEALLOCATE (pdos_curve)

   END SUBROUTINE write_broadened_pdos

END MODULE qs_pdos
