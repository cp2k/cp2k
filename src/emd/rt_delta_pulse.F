!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to apply a delta pulse for RTP and EMD
! **************************************************************************************************

MODULE rt_delta_pulse
   USE bibliography,                    ONLY: Mattiat2019,&
                                              Mattiat2022,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE commutator_rpnl,                 ONLY: build_com_mom_nl,&
                                              build_com_nl_mag
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              rtp_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_get_info, &
        dbcsr_init_p, dbcsr_p_type, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric, &
        dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_triangular_multiply,&
                                              cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert,&
                                              cp_fm_cholesky_reduce,&
                                              cp_fm_cholesky_restore
   USE cp_fm_diag,                      ONLY: cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: one,&
                                              twopi,&
                                              zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE moments_utils,                   ONLY: get_reference_point
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE qs_dftb_matrices,                ONLY: build_dftb_overlap
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: build_berry_moment_matrix,&
                                              build_local_magmom_matrix,&
                                              build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE rt_propagation_types,            ONLY: get_rtp,&
                                              rt_prop_create_mos,&
                                              rt_prop_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_delta_pulse'

   PUBLIC :: apply_delta_pulse

CONTAINS

! **************************************************************************************************
!> \brief Interface to call the delta pulse depending on the type of calculation.
!> \param qs_env ...
!> \param rtp ...
!> \param rtp_control ...
!> \author Update: Guillaume Le Breton (2023.01)
! **************************************************************************************************

   SUBROUTINE apply_delta_pulse(qs_env, rtp, rtp_control)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control

      CHARACTER(LEN=3), DIMENSION(3)                     :: rlab
      INTEGER                                            :: i, output_unit
      LOGICAL                                            :: my_apply_pulse, periodic_cell
      REAL(KIND=dp), DIMENSION(3)                        :: kvec
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new, mos_old
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(section_vals_type), POINTER                   :: input, rtp_section

      NULLIFY (logger, input, rtp_section)

      logger => cp_get_default_logger()
      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      input=input, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s)
      rtp_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION")
      output_unit = cp_print_key_unit_nr(logger, rtp_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")
      rlab = [CHARACTER(LEN=3) :: "X", "Y", "Z"]
      periodic_cell = ANY(cell%perd > 0)
      my_apply_pulse = .TRUE.
      CALL get_qs_env(qs_env, mos=mos)

      IF (rtp%linear_scaling) THEN
         IF (.NOT. ASSOCIATED(mos)) THEN
            CALL cp_warn(__LOCATION__, "Delta Pulse not implemented for Linear-Scaling based ground "// &
                         "state calculation. If you want to perform a Linear-Scaling RTP from a "// &
                         "Linear-Scaling GS calculation you can do the following: (i) LSCF froms "// &
                         "scratch, (ii) MO-based SCF (for 1 SCF loop for instance) with the LSCF "// &
                         "result as a restart and (iii) linear scaling RTP + delta kick (for 1 "// &
                         "SCF loop for instance).")
            my_apply_pulse = .FALSE.
         ELSE
            ! create temporary mos_old and mos_new to use delta kick routine designed for MOs-based RTP
            CALL rt_prop_create_mos(rtp, mos, qs_env%mpools, dft_control, &
                                    init_mos_old=.TRUE., init_mos_new=.TRUE., &
                                    init_mos_next=.FALSE., init_mos_admn=.FALSE.)
         END IF
      END IF

      IF (my_apply_pulse) THEN
         ! The amplitude of the perturbation for all the method, modulo some prefactor:
         kvec(:) = cell%h_inv(1, :)*rtp_control%delta_pulse_direction(1) + &
                   cell%h_inv(2, :)*rtp_control%delta_pulse_direction(2) + &
                   cell%h_inv(3, :)*rtp_control%delta_pulse_direction(3)
         kvec = kvec*twopi*rtp_control%delta_pulse_scale

         CALL get_rtp(rtp=rtp, mos_old=mos_old, mos_new=mos_new)
         IF (rtp_control%apply_delta_pulse) THEN
            IF (dft_control%qs_control%dftb) &
               CALL build_dftb_overlap(qs_env, 1, matrix_s)
            IF (rtp_control%periodic) THEN
               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(/,(T3,A,T40))") &
                     "An Electric Delta Kick within periodic condition is applied before running RTP.  "// &
                     "Its amplitude in atomic unit is:"
                  WRITE (output_unit, "(T3,3(A,A,E16.8,1X))") &
                     (TRIM(rlab(i)), "=", -kvec(i), i=1, 3)
               END IF
               CALL apply_delta_pulse_electric_periodic(qs_env, mos_old, mos_new, -kvec)
            ELSE
               IF (periodic_cell) THEN
                  CPWARN("This application of the delta pulse is not compatible with PBC!")
               END IF
               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(/,(T3,A,T40))") &
                     "An Electric Delta Kick within the length gauge is applied before running RTP.  "// &
                     "Its amplitude in atomic unit is:"
                  WRITE (output_unit, "(T3,3(A,A,E16.8,1X))") &
                     (TRIM(rlab(i)), "=", -kvec(i), i=1, 3)
               END IF
               CALL apply_delta_pulse_electric(qs_env, mos_old, mos_new, -kvec)
            END IF
         ELSE IF (rtp_control%apply_delta_pulse_mag) THEN
            IF (periodic_cell) THEN
               CPWARN("This application of the delta pulse is not compatible with PBC!")
            END IF
            ! The prefactor (strength of the magnetic field, should be divided by 2c)
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T40))") &
                  "A Magnetic Delta Kick is applied before running RTP.  "// &
                  "Its amplitude in atomic unit is:"
               WRITE (output_unit, "(T3,3(A,A,E16.8,1X))") &
                  (TRIM(rlab(i)), "=", -kvec(i)/2, i=1, 3)
            END IF
            CALL apply_delta_pulse_mag(qs_env, mos_old, mos_new, -kvec(:)/2)
         ELSE
            CPABORT("Code error: this case should not happen!")
         END IF
      END IF

   END SUBROUTINE apply_delta_pulse

! **************************************************************************************************
!> \brief uses perturbation theory to get the proper initial conditions
!>        The len_rep option is NOT compatible with periodic boundary conditions!
!> \param qs_env ...
!> \param mos_old ...
!> \param mos_new ...
!> \param kvec ...
!> \author Joost & Martin (2011)
! **************************************************************************************************

   SUBROUTINE apply_delta_pulse_electric_periodic(qs_env, mos_old, mos_new, kvec)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_old, mos_new
      REAL(KIND=dp), DIMENSION(3)                        :: kvec

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_delta_pulse_electric_periodic'

      INTEGER                                            :: handle, icol, idir, irow, ispin, nao, &
                                                            ncol_local, nmo, nrow_local, nvirt, &
                                                            reference
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: com_nl, len_rep, periodic_cell
      REAL(KIND=dp)                                      :: eps_ppnl, factor
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      REAL(KIND=dp), DIMENSION(3)                        :: rcc
      REAL(kind=dp), DIMENSION(:), POINTER               :: eigenvalues, ref_point
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, fm_struct_tmp
      TYPE(cp_fm_type)                                   :: eigenvectors, mat_ks, mat_tmp, momentum, &
                                                            S_chol, virtuals
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_r, matrix_rv, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      NULLIFY (cell, mos, rtp, matrix_s, matrix_ks, input, dft_control, particle_set, fm_struct)
      ! we need the overlap and ks matrix for a full diagionalization
      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      mos=mos, &
                      rtp=rtp, &
                      matrix_s=matrix_s, &
                      matrix_ks=matrix_ks, &
                      dft_control=dft_control, &
                      input=input, &
                      particle_set=particle_set)

      rtp_control => dft_control%rtp_control
      periodic_cell = ANY(cell%perd > 0)

      ! relevant input parameters
      com_nl = section_get_lval(section_vals=input, keyword_name="DFT%REAL_TIME_PROPAGATION%COM_NL")
      len_rep = section_get_lval(section_vals=input, keyword_name="DFT%REAL_TIME_PROPAGATION%LEN_REP")

      ! calculate non-local commutator if necessary
      IF (com_nl) THEN
         CALL cite_reference(Mattiat2019)
         NULLIFY (qs_kind_set, sab_orb, sap_ppnl)
         CALL get_qs_env(qs_env, &
                         sap_ppnl=sap_ppnl, &
                         sab_orb=sab_orb, &
                         qs_kind_set=qs_kind_set)
         eps_ppnl = dft_control%qs_control%eps_ppnl

         NULLIFY (matrix_rv)
         CALL dbcsr_allocate_matrix_set(matrix_rv, 3)
         DO idir = 1, 3
            CALL dbcsr_init_p(matrix_rv(idir)%matrix)
            CALL dbcsr_create(matrix_rv(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_rv(idir)%matrix, sab_orb)
            CALL dbcsr_set(matrix_rv(idir)%matrix, 0._dp)
         END DO
         CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, cell, matrix_rv=matrix_rv)
      END IF

      ! calculate dipole moment matrix if required, NOT for periodic boundary conditions!
      IF (len_rep) THEN
         CALL cite_reference(Mattiat2022)
         IF (periodic_cell) THEN
            CPWARN("This application of the delta pulse is not compatible with PBC!")
         END IF
         ! get reference point
         reference = section_get_ival(section_vals=input, &
                                      keyword_name="DFT%PRINT%MOMENTS%REFERENCE")
         NULLIFY (ref_point)
         CALL section_vals_val_get(input, "DFT%PRINT%MOMENTS%REF_POINT", r_vals=ref_point)
         CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)

         NULLIFY (sab_orb)
         CALL get_qs_env(qs_env, sab_orb=sab_orb)
         ! calculate dipole moment operator
         NULLIFY (matrix_r)
         CALL dbcsr_allocate_matrix_set(matrix_r, 3)
         DO idir = 1, 3
            CALL dbcsr_init_p(matrix_r(idir)%matrix)
            CALL dbcsr_create(matrix_r(idir)%matrix, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_symmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_r(idir)%matrix, sab_orb)
            CALL dbcsr_set(matrix_r(idir)%matrix, 0._dp)
         END DO
         CALL build_local_moment_matrix(qs_env, matrix_r, 1, rcc)
      END IF

      IF (rtp_control%velocity_gauge) THEN
         rtp_control%vec_pot = rtp_control%vec_pot + kvec
      END IF

      ! struct for fm matrices
      fm_struct => rtp%ao_ao_fmstruct

      ! create matrices and get Cholesky decomposition of S
      CALL cp_fm_create(mat_ks, matrix_struct=fm_struct, name="mat_ks")
      CALL cp_fm_create(eigenvectors, matrix_struct=fm_struct, name="eigenvectors")
      CALL cp_fm_create(S_chol, matrix_struct=fm_struct, name="S_chol")
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, S_chol)
      CALL cp_fm_cholesky_decompose(S_chol)

      ! get number of atomic orbitals
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

      DO ispin = 1, SIZE(matrix_ks)
         ! diagonalize KS matrix to get occ and virt mos
         ALLOCATE (eigenvalues(nao))
         CALL cp_fm_create(mat_tmp, matrix_struct=fm_struct, name="mat_tmp")
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, mat_ks)
         CALL cp_fm_cholesky_reduce(mat_ks, S_chol)
         CALL cp_fm_syevd(mat_ks, mat_tmp, eigenvalues)
         CALL cp_fm_cholesky_restore(mat_tmp, nao, S_chol, eigenvectors, "SOLVE")

         ! virtuals
         CALL get_mo_set(mo_set=mos(ispin), nmo=nmo)
         nvirt = nao - nmo
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=fm_struct%para_env, context=fm_struct%context, &
                                  nrow_global=nao, ncol_global=nvirt)
         CALL cp_fm_create(virtuals, matrix_struct=fm_struct_tmp, name="virtuals")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_to_fm(eigenvectors, virtuals, nvirt, nmo + 1, 1)

         ! occupied
         CALL cp_fm_to_fm(eigenvectors, mos_old(2*ispin - 1), nmo, 1, 1)

         CALL cp_fm_struct_create(fm_struct_tmp, para_env=fm_struct%para_env, context=fm_struct%context, &
                                  nrow_global=nvirt, ncol_global=nmo)
         CALL cp_fm_create(momentum, matrix_struct=fm_struct_tmp, name="momentum")
         CALL cp_fm_struct_release(fm_struct_tmp)

         ! the momentum operator (in a given direction)
         CALL cp_fm_set_all(mos_new(2*ispin - 1), 0.0_dp)

         DO idir = 1, 3
            factor = kvec(idir)
            IF (factor .NE. 0.0_dp) THEN
               IF (.NOT. len_rep) THEN
                  CALL cp_dbcsr_sm_fm_multiply(matrix_s(idir + 1)%matrix, mos_old(2*ispin - 1), &
                                               mos_old(2*ispin), ncol=nmo)
               ELSE
                  CALL cp_dbcsr_sm_fm_multiply(matrix_r(idir)%matrix, mos_old(2*ispin - 1), &
                                               mos_old(2*ispin), ncol=nmo)
               END IF

               CALL cp_fm_scale_and_add(1.0_dp, mos_new(2*ispin - 1), factor, mos_old(2*ispin))
               IF (com_nl) THEN
                  CALL cp_fm_set_all(mos_old(2*ispin), 0.0_dp)
                  CALL cp_dbcsr_sm_fm_multiply(matrix_rv(idir)%matrix, mos_old(2*ispin - 1), &
                                               mos_old(2*ispin), ncol=nmo)
                  CALL cp_fm_scale_and_add(1.0_dp, mos_new(2*ispin - 1), factor, mos_old(2*ispin))
               END IF
            END IF
         END DO

         CALL parallel_gemm('T', 'N', nvirt, nmo, nao, 1.0_dp, virtuals, mos_new(2*ispin - 1), 0.0_dp, momentum)

         ! the tricky bit ... rescale by the eigenvalue difference
         IF (.NOT. len_rep) THEN
            CALL cp_fm_get_info(momentum, nrow_local=nrow_local, ncol_local=ncol_local, &
                                row_indices=row_indices, col_indices=col_indices, local_data=local_data)
            DO icol = 1, ncol_local
               DO irow = 1, nrow_local
                  factor = 1/(eigenvalues(col_indices(icol)) - eigenvalues(nmo + row_indices(irow)))
                  local_data(irow, icol) = factor*local_data(irow, icol)
               END DO
            END DO
         END IF
         CALL cp_fm_release(mat_tmp)
         DEALLOCATE (eigenvalues)

         ! now obtain the initial condition in mos_old
         CALL cp_fm_to_fm(eigenvectors, mos_old(2*ispin - 1), nmo, 1, 1)
         CALL parallel_gemm("N", "N", nao, nmo, nvirt, 1.0_dp, virtuals, momentum, 0.0_dp, mos_old(2*ispin))

         CALL cp_fm_release(virtuals)
         CALL cp_fm_release(momentum)
      END DO

      ! release matrices
      CALL cp_fm_release(S_chol)
      CALL cp_fm_release(mat_ks)
      CALL cp_fm_release(eigenvectors)
      IF (com_nl) CALL dbcsr_deallocate_matrix_set(matrix_rv)
      IF (len_rep) CALL dbcsr_deallocate_matrix_set(matrix_r)

      ! orthonormalize afterwards
      CALL orthonormalize_complex_mos(qs_env, mos_old)

      CALL timestop(handle)

   END SUBROUTINE apply_delta_pulse_electric_periodic

! **************************************************************************************************
!> \brief applies exp(ikr) to the wavefunction.... stored in mos_old...
!> \param qs_env ...
!> \param mos_old ...
!> \param mos_new ...
!> \param kvec ...
!> \author Joost & Martin (2011)
! **************************************************************************************************

   SUBROUTINE apply_delta_pulse_electric(qs_env, mos_old, mos_new, kvec)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_old, mos_new
      REAL(KIND=dp), DIMENSION(3)                        :: kvec

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_delta_pulse_electric'

      INTEGER                                            :: handle, i, nao, nmo
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: S_inv_fm, tmp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control

      CALL timeset(routineN, handle)
      NULLIFY (cell, dft_control, matrix_s, mos, rtp, rtp_control)
      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      mos=mos, &
                      rtp=rtp)
      rtp_control => dft_control%rtp_control

      IF (rtp_control%velocity_gauge) THEN
         rtp_control%vec_pot = rtp_control%vec_pot + kvec
      END IF

      ! calculate exponentials (= Berry moments)
      NULLIFY (cosmat, sinmat)
      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_copy(cosmat, matrix_s(1)%matrix, 'COS MOM')
      CALL dbcsr_copy(sinmat, matrix_s(1)%matrix, 'SIN MOM')
      CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec)

      ! need inverse of overlap matrix
      CALL cp_fm_create(S_inv_fm, matrix_struct=rtp%ao_ao_fmstruct, name="S_inv_fm")
      CALL cp_fm_create(tmp, matrix_struct=rtp%ao_ao_fmstruct, name="tmp_mat")
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, S_inv_fm)
      CALL cp_fm_cholesky_decompose(S_inv_fm)
      CALL cp_fm_cholesky_invert(S_inv_fm)
      CALL cp_fm_upper_to_full(S_inv_fm, tmp)
      CALL cp_fm_release(tmp)

      DO i = 1, SIZE(mos)
         ! apply exponentials to mo coefficients
         CALL get_mo_set(mos(i), nao=nao, nmo=nmo)
         CALL cp_dbcsr_sm_fm_multiply(cosmat, mos(i)%mo_coeff, mos_new(2*i - 1), ncol=nmo)
         CALL cp_dbcsr_sm_fm_multiply(sinmat, mos(i)%mo_coeff, mos_new(2*i), ncol=nmo)

         CALL parallel_gemm("N", "N", nao, nmo, nao, 1.0_dp, S_inv_fm, mos_new(2*i - 1), 0.0_dp, mos_old(2*i - 1))
         CALL parallel_gemm("N", "N", nao, nmo, nao, 1.0_dp, S_inv_fm, mos_new(2*i), 0.0_dp, mos_old(2*i))
      END DO

      CALL cp_fm_release(S_inv_fm)
      CALL dbcsr_deallocate_matrix(cosmat)
      CALL dbcsr_deallocate_matrix(sinmat)

      ! orthonormalize afterwards
      CALL orthonormalize_complex_mos(qs_env, mos_old)

      CALL timestop(handle)

   END SUBROUTINE apply_delta_pulse_electric

! **************************************************************************************************
!> \brief apply magnetic delta pulse to linear order
!> \param qs_env ...
!> \param mos_old ...
!> \param mos_new ...
!> \param kvec ...
! **************************************************************************************************
   SUBROUTINE apply_delta_pulse_mag(qs_env, mos_old, mos_new, kvec)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_old, mos_new
      REAL(KIND=dp), DIMENSION(3)                        :: kvec

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_delta_pulse_mag'

      INTEGER                                            :: gauge_orig, handle, idir, ispin, nao, &
                                                            nmo, nrow_global, nvirt
      REAL(KIND=dp)                                      :: eps_ppnl, factor
      REAL(KIND=dp), DIMENSION(3)                        :: rcc
      REAL(kind=dp), DIMENSION(:), POINTER               :: eigenvalues, ref_point
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: eigenvectors, mat_ks, perturbation, &
                                                            S_chol, virtuals
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_mag, matrix_nl, &
                                                            matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      CALL cite_reference(Mattiat2022)

      NULLIFY (rtp, dft_control, matrix_ks, matrix_s, input, mos, cell, sab_orb, sab_all, sap_ppnl, &
               qs_kind_set, particle_set)

      CALL get_qs_env(qs_env, &
                      rtp=rtp, &
                      dft_control=dft_control, &
                      mos=mos, &
                      matrix_ks=matrix_ks, &
                      matrix_s=matrix_s, &
                      input=input, &
                      cell=cell, &
                      sab_orb=sab_orb, &
                      sab_all=sab_all, &
                      sap_ppnl=sap_ppnl)

      gauge_orig = section_get_ival(section_vals=input, &
                                    keyword_name="DFT%REAL_TIME_PROPAGATION%GAUGE_ORIG")
      NULLIFY (ref_point)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%GAUGE_ORIG_MANUAL", r_vals=ref_point)
      CALL get_reference_point(rcc, qs_env=qs_env, reference=gauge_orig, ref_point=ref_point)

      ! Create fm matrices
      CALL cp_fm_create(S_chol, matrix_struct=rtp%ao_ao_fmstruct, name='Cholesky S')
      CALL cp_fm_create(eigenvectors, matrix_struct=rtp%ao_ao_fmstruct, name="gs evecs fm")
      CALL cp_fm_create(mat_ks, matrix_struct=rtp%ao_ao_fmstruct, name='KS matrix')

      ! get nrows_global
      CALL cp_fm_get_info(mat_ks, nrow_global=nrow_global)

      ! cholesky decomposition of overlap matrix
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, S_chol)
      CALL cp_fm_cholesky_decompose(S_chol)

      ! initiate perturbation matrix
      NULLIFY (matrix_mag)
      CALL dbcsr_allocate_matrix_set(matrix_mag, 3)
      DO idir = 1, 3
         CALL dbcsr_init_p(matrix_mag(idir)%matrix)
         CALL dbcsr_create(matrix_mag(idir)%matrix, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_mag(idir)%matrix, sab_orb)
         CALL dbcsr_set(matrix_mag(idir)%matrix, 0._dp)
      END DO
      ! construct magnetic dipole moment matrix
      CALL build_local_magmom_matrix(qs_env, matrix_mag, 1, ref_point=rcc)

      ! work matrix for non-local potential part if necessary
      NULLIFY (matrix_nl)
      IF (ASSOCIATED(sap_ppnl)) THEN
         CALL dbcsr_allocate_matrix_set(matrix_nl, 3)
         DO idir = 1, 3
            CALL dbcsr_init_p(matrix_nl(idir)%matrix)
            CALL dbcsr_create(matrix_nl(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_nl(idir)%matrix, sab_orb)
            CALL dbcsr_set(matrix_nl(idir)%matrix, 0._dp)
         END DO
         ! construct non-local contribution
         CALL get_qs_env(qs_env, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set)
         eps_ppnl = dft_control%qs_control%eps_ppnl

         CALL build_com_nl_mag(qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, matrix_nl, rcc, cell)

         DO idir = 1, 3
            CALL dbcsr_add(matrix_mag(idir)%matrix, matrix_nl(idir)%matrix, -one, one)
         END DO

         CALL dbcsr_deallocate_matrix_set(matrix_nl)
      END IF

      DO ispin = 1, dft_control%nspins
         ! allocate eigenvalues
         NULLIFY (eigenvalues)
         ALLOCATE (eigenvalues(nrow_global))
         ! diagonalize KS matrix in AO basis using Cholesky decomp. of S
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, mat_ks)
         CALL cp_fm_cholesky_reduce(mat_ks, S_chol)
         CALL cp_fm_syevd(mat_ks, eigenvectors, eigenvalues)
         CALL cp_fm_triangular_multiply(S_chol, eigenvectors, invert_tr=.TRUE.)

         ! virtuals
         CALL get_mo_set(mo_set=mos(ispin), nao=nao, nmo=nmo)
         nvirt = nao - nmo
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=rtp%ao_ao_fmstruct%para_env, context=rtp%ao_ao_fmstruct%context, &
                                  nrow_global=nrow_global, ncol_global=nvirt)
         CALL cp_fm_create(virtuals, matrix_struct=fm_struct_tmp, name="virtuals")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_to_fm(eigenvectors, virtuals, nvirt, nmo + 1, 1)

         ! occupied
         CALL cp_fm_to_fm(eigenvectors, mos_old(2*ispin - 1), nmo, 1, 1)

         CALL cp_fm_struct_create(fm_struct_tmp, para_env=rtp%ao_ao_fmstruct%para_env, context=rtp%ao_ao_fmstruct%context, &
                                  nrow_global=nvirt, ncol_global=nmo)
         CALL cp_fm_create(perturbation, matrix_struct=fm_struct_tmp, name="perturbation")
         CALL cp_fm_struct_release(fm_struct_tmp)

         ! apply perturbation
         CALL cp_fm_set_all(mos_new(2*ispin - 1), 0.0_dp)

         DO idir = 1, 3
            factor = kvec(idir)
            IF (factor .NE. 0.0_dp) THEN
               CALL cp_dbcsr_sm_fm_multiply(matrix_mag(idir)%matrix, mos_old(2*ispin - 1), &
                                            mos_old(2*ispin), ncol=nmo)
               CALL cp_fm_scale_and_add(1.0_dp, mos_new(2*ispin - 1), factor, mos_old(2*ispin))
            END IF
         END DO

         CALL parallel_gemm('T', 'N', nvirt, nmo, nao, 1.0_dp, virtuals, mos_new(2*ispin - 1), 0.0_dp, perturbation)

         DEALLOCATE (eigenvalues)

         ! now obtain the initial condition in mos_old
         CALL cp_fm_to_fm(eigenvectors, mos_old(2*ispin - 1), nmo, 1, 1)
         CALL parallel_gemm("N", "N", nao, nmo, nvirt, 1.0_dp, virtuals, perturbation, 0.0_dp, mos_old(2*ispin))

         CALL cp_fm_release(virtuals)
         CALL cp_fm_release(perturbation)
      END DO

      ! deallocations
      CALL cp_fm_release(S_chol)
      CALL cp_fm_release(mat_ks)
      CALL cp_fm_release(eigenvectors)
      CALL dbcsr_deallocate_matrix_set(matrix_mag)

      ! orthonormalize afterwards
      CALL orthonormalize_complex_mos(qs_env, mos_old)

      CALL timestop(handle)

   END SUBROUTINE apply_delta_pulse_mag

! **************************************************************************************************
!> \brief orthonormalize complex mos, e. g. after non-unitary transformations using Löwdin's algorithm
!> \param qs_env ...
!> \param coeffs ...
! **************************************************************************************************
   SUBROUTINE orthonormalize_complex_mos(qs_env, coeffs)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: coeffs

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: eigenvalues_sqrt
      INTEGER                                            :: im, ispin, j, nao, nmo, nspins, re
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: oo_c, oo_v, oo_vt
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: oo_1, oo_2, S_fm, tmp
      TYPE(cp_fm_type), DIMENSION(2)                     :: coeffs_tmp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (para_env, blacs_env, dft_control, matrix_s, mos)
      CALL get_qs_env(qs_env, &
                      blacs_env=blacs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      mos=mos, &
                      para_env=para_env)
      nspins = dft_control%nspins
      CALL cp_fm_get_info(coeffs(1), nrow_global=nao)

      ! get overlap matrix
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nao, &
                               context=blacs_env, para_env=para_env)
      CALL cp_fm_create(S_fm, matrix_struct=fm_struct_tmp, name="overlap fm")
      CALL cp_fm_struct_release(fm_struct_tmp)
      ! copy overlap matrix
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, S_fm)

      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), nmo=nmo)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=nmo, ncol_global=nmo)
         CALL cp_fm_create(oo_1, matrix_struct=fm_struct_tmp, name="oo_1")
         CALL cp_fm_create(oo_2, matrix_struct=fm_struct_tmp, name="oo_2")
         CALL cp_fm_struct_release(fm_struct_tmp)

         CALL cp_fm_create(tmp, matrix_struct=coeffs(2*ispin - 1)%matrix_struct, name="tmp_mat")
         ! get the complex overlap matrix in MO basis
         ! x^T S x + y^T S y + i (-y^TS x+x^T S y)
         CALL parallel_gemm("N", "N", nao, nmo, nao, 1.0_dp, S_fm, coeffs(2*ispin - 1), 0.0_dp, tmp)
         CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, coeffs(2*ispin - 1), tmp, 0.0_dp, oo_1)
         CALL parallel_gemm("T", "N", nmo, nmo, nao, -1.0_dp, coeffs(2*ispin), tmp, 0.0_dp, oo_2)

         CALL parallel_gemm("N", "N", nao, nmo, nao, 1.0_dp, S_fm, coeffs(2*ispin), 0.0_dp, tmp)
         CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, coeffs(2*ispin), tmp, 1.0_dp, oo_1)
         CALL parallel_gemm("T", "N", nmo, nmo, nao, 1.0_dp, coeffs(2*ispin - 1), tmp, 1.0_dp, oo_2)
         CALL cp_fm_release(tmp)

         ! complex Löwdin
         CALL cp_cfm_create(oo_c, oo_1%matrix_struct)
         CALL cp_cfm_create(oo_v, oo_1%matrix_struct)
         CALL cp_cfm_create(oo_vt, oo_1%matrix_struct)
         oo_c%local_data = CMPLX(oo_1%local_data, oo_2%local_data, KIND=dp)

         ALLOCATE (eigenvalues(nmo))
         ALLOCATE (eigenvalues_sqrt(nmo))
         CALL cp_cfm_heevd(oo_c, oo_v, eigenvalues)
         eigenvalues_sqrt(:) = CMPLX(one/SQRT(eigenvalues(:)), zero, dp)
         CALL cp_cfm_to_cfm(oo_v, oo_vt)
         CALL cp_cfm_column_scale(oo_v, eigenvalues_sqrt)
         DEALLOCATE (eigenvalues)
         DEALLOCATE (eigenvalues_sqrt)
         CALL parallel_gemm('N', 'C', nmo, nmo, nmo, (1.0_dp, 0.0_dp), &
                            oo_v, oo_vt, (0.0_dp, 0.0_dp), oo_c)
         oo_1%local_data = REAL(oo_c%local_data, KIND=dp)
         oo_2%local_data = AIMAG(oo_c%local_data)
         CALL cp_cfm_release(oo_c)
         CALL cp_cfm_release(oo_v)
         CALL cp_cfm_release(oo_vt)

         ! transform coefficients accordingly
         DO j = 1, 2
            CALL cp_fm_create(coeffs_tmp(j), matrix_struct=coeffs(2*(ispin - 1) + j)%matrix_struct)
         END DO

         ! indices for coeffs_tmp
         re = 1
         im = 2
         CALL parallel_gemm("N", "N", nao, nmo, nmo, one, coeffs(2*ispin - 1), oo_1, zero, coeffs_tmp(re))
         CALL parallel_gemm("N", "N", nao, nmo, nmo, one, coeffs(2*ispin - 1), oo_2, zero, coeffs_tmp(im))

         CALL parallel_gemm("N", "N", nao, nmo, nmo, -one, coeffs(2*ispin), oo_2, zero, coeffs(2*ispin - 1))
         CALL cp_fm_scale_and_add(one, coeffs(2*ispin - 1), one, coeffs_tmp(re))

         CALL parallel_gemm("N", "N", nao, nmo, nmo, one, coeffs(2*ispin), oo_1, one, coeffs_tmp(im))
         CALL cp_fm_to_fm(coeffs_tmp(im), coeffs(2*ispin))

         DO j = 1, 2
            CALL cp_fm_release(coeffs_tmp(j))
         END DO
         CALL cp_fm_release(oo_1)
         CALL cp_fm_release(oo_2)
      END DO
      CALL cp_fm_release(S_fm)

   END SUBROUTINE orthonormalize_complex_mos

END MODULE rt_delta_pulse
