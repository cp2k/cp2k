!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines printing various data related to POD method
!> \author Z. Futera (03.2019)
! **************************************************************************************************
MODULE et_coupling_pod_print

   USE cp_blacs_calls,                  ONLY: cp_blacs_gridexit,&
                                              cp_blacs_gridinfo,&
                                              cp_blacs_gridinit
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_should_output
   USE et_coupling_pod_mo,              ONLY: pod_mo_sum_c2
   USE et_coupling_pod_types,           ONLY: et_cpl
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_recv,&
                                              mp_send
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: mo_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod_print'

   PUBLIC :: pod_print_couplings, &
             pod_print_state_fracs, &
             pod_print_states

CONTAINS

! **************************************************************************************************
!> \brief Print out specific MO coefficients
!> \param output_unit unit number of the open output stream
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param blk atomic-block ID
!> \param n_spins number of spin components
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_print_state_fracs(output_unit, qs_env, ec, blk, n_spins)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: blk, n_spins

      CHARACTER(len=*), PARAMETER :: routineN = 'pod_print_state_fracs', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: j, k, l, m, n, n_ao, n_mo
      INTEGER, DIMENSION(:), POINTER                     :: list_at, list_mo
      REAL(KIND=dp)                                      :: c1, c2
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mat_w
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: block_sec, print_sec

! Routine name for debug purposes
! Local variables

      ! Initialization
      NULLIFY (block_sec)
      NULLIFY (print_sec)
      NULLIFY (qs_kind_set)

      ! Atomic block data
      block_sec => section_vals_get_subs_vals(qs_env%input, &
                                              'PROPERTIES%ET_COUPLING%PROJECTION%BLOCK')

      print_sec => section_vals_get_subs_vals(block_sec, 'PRINT', i_rep_section=blk)

      ! List of atoms
      CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM', n_rep_val=n)

      IF (n > 0) THEN

         IF (output_unit > 0) &
            WRITE (output_unit, '(/,T3,A/)') 'Block state fractions:'

         ! Number of AO functions
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
         CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao)

         ! MOs in orthonormal basis set
         ALLOCATE (mat_w(n_spins))
         CPASSERT(ASSOCIATED(mat_w))
         DO j = 1, n_spins
            n_mo = ec%block(blk)%n_ao
            CALL cp_fm_create(matrix=mat_w(j)%matrix, &
                              matrix_struct=ec%block(blk)%mo(j)%mo_set%mo_coeff%matrix_struct, &
                              name='BLOCK MOs IN ORTHONORMAL BASIS SET')
            CALL cp_gemm("N", "N", n_ao, n_mo, n_ao, 1.0_dp, ec%tmatr, &
                         ec%block(blk)%mo(j)%mo_set%mo_coeff, 0.0_dp, mat_w(j)%matrix)
         END DO

         DO j = 1, n
            NULLIFY (list_at)
            CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM', &
                                      i_rep_val=j, i_vals=list_at)
            IF (ASSOCIATED(list_at)) THEN

               ! List of states
               CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM_STATE', n_rep_val=m)

               IF (m > 0) THEN

                  DO k = 1, m
                     NULLIFY (list_mo)
                     CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM_STATE', &
                                               i_rep_val=k, i_vals=list_mo)
                     IF (ASSOCIATED(list_mo)) THEN

                        IF (j > 1) THEN
                           IF (output_unit > 0) &
                              WRITE (output_unit, *)
                        END IF

                        DO l = 1, SIZE(list_mo)

                           IF (n_spins > 1) THEN
                              c1 = pod_mo_sum_c2(ec%block(blk)%atom, mat_w(1)%matrix, &
                                                 list_mo(l), list_at)
                              c2 = pod_mo_sum_c2(ec%block(blk)%atom, mat_w(2)%matrix, &
                                                 list_mo(l), list_at)
                              IF (output_unit > 0) &
                                 WRITE (output_unit, '(I5,A,I5,2F20.10)') j, ' /', list_mo(l), c1, c2
                           ELSE
                              c1 = pod_mo_sum_c2(ec%block(blk)%atom, mat_w(1)%matrix, &
                                                 list_mo(l), list_at)
                              IF (output_unit > 0) &
                                 WRITE (output_unit, '(I5,A,I5,F20.10)') j, ' /', list_mo(l), c1
                           END IF

                        END DO

                     END IF
                  END DO

               END IF

            END IF
         END DO

         ! Clean memory
         DO j = 1, n_spins
            CALL cp_fm_release(matrix=mat_w(j)%matrix)
         END DO
         IF (ASSOCIATED(mat_w)) &
            DEALLOCATE (mat_w)

      END IF

   END SUBROUTINE pod_print_state_fracs

! **************************************************************************************************
!> \brief Print out electronic states (MOs)
!> \param output_unit unit number of the open output stream
!> \param mo array of MO sets
!> \param n_spins number of spin components
!> \param label output label
!> \param mx_mo_a maximum number of alpha states to print out
!> \param mx_mo_b maximum number of beta states to print out
!> \param fermi print out Fermi level and number of electrons
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_print_states(output_unit, mo, n_spins, label, mx_mo_a, mx_mo_b, fermi)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo
      INTEGER                                            :: n_spins
      CHARACTER(LEN=*)                                   :: label
      INTEGER, OPTIONAL                                  :: mx_mo_a, mx_mo_b
      LOGICAL, OPTIONAL                                  :: fermi

      CHARACTER(len=*), PARAMETER :: routineN = 'pod_print_states', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, mx_a, mx_b, n
      LOGICAL                                            :: prnt_fm

! Routine name for debug purposes
! Local variables

      prnt_fm = .FALSE.
      IF (PRESENT(fermi)) &
         prnt_fm = fermi

      IF (output_unit > 0) THEN

         WRITE (output_unit, '(/,T3,A/)') 'State energies ('//TRIM(ADJUSTL(label))//'):'

         ! Spin-polarized calculation
         IF (n_spins > 1) THEN

            mx_a = mo(1)%mo_set%nmo
            IF (PRESENT(mx_mo_a)) &
               mx_a = MIN(mo(1)%mo_set%nmo, mx_mo_a)
            mx_b = mo(2)%mo_set%nmo
            IF (PRESENT(mx_mo_b)) &
               mx_b = MIN(mo(2)%mo_set%nmo, mx_mo_b)
            n = MAX(mx_a, mx_b)

            DO i = 1, n
               WRITE (output_unit, '(T3,I10)', ADVANCE='no') i
               IF (i <= mx_a) THEN
                  WRITE (output_unit, '(2F12.4)', ADVANCE='no') &
                     mo(1)%mo_set%occupation_numbers(i), mo(1)%mo_set%eigenvalues(i)
               ELSE
                  WRITE (output_unit, '(A)', ADVANCE='no') '                        '
               END IF
               WRITE (output_unit, '(A)', ADVANCE='no') '     '
               IF (i <= mx_b) THEN
                  WRITE (output_unit, '(2F12.4)') &
                     mo(2)%mo_set%occupation_numbers(i), mo(2)%mo_set%eigenvalues(i)
               ELSE
                  WRITE (output_unit, *)
               END IF
            END DO

            IF (prnt_fm) THEN
               WRITE (output_unit, '(/,T3,I10,F24.4,I10,F19.4)') &
                  mo(1)%mo_set%nelectron, mo(1)%mo_set%mu, &
                  mo(2)%mo_set%nelectron, mo(2)%mo_set%mu
            END IF

            ! Spin-restricted calculation
         ELSE

            mx_a = mo(1)%mo_set%nmo
            IF (PRESENT(mx_mo_a)) &
               mx_a = MIN(mo(1)%mo_set%nmo, mx_mo_a)

            DO i = 1, mx_a
               WRITE (output_unit, '(T3,I10,2F12.4)') &
                  i, mo(1)%mo_set%occupation_numbers(i), mo(1)%mo_set%eigenvalues(i)
            END DO

            IF (prnt_fm) THEN
               WRITE (output_unit, '(/,T3,I10,F24.4)') &
                  mo(1)%mo_set%nelectron, mo(1)%mo_set%mu
            END IF

         END IF

      END IF

   END SUBROUTINE pod_print_states

! **************************************************************************************************
!> \brief Print out donor-acceptor state couplings
!> \param output_unit unit number of the open output stream
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param n_states number of states
!> \param n_beta_states number of beta-spin states
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_print_couplings(output_unit, qs_env, ec, n_states, n_beta_states)

      IMPLICIT NONE

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: n_states
      INTEGER, OPTIONAL                                  :: n_beta_states

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_print_couplings', &
                                     routineP = moduleN//':'//routineN

      ! Local variables
      LOGICAL                                            :: do_print
      INTEGER                                            :: i, j, k, l, id
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: et_proj_sec, print_key
#if defined(__SCALAPACK)
      INTEGER                                            :: n, cntxt_loc, info, tg1, tg2
      INTEGER                                            :: ib, ix, ic, ir, ir_loc, ic_loc
      INTEGER                                            :: n_b_row, n_b_col
      INTEGER                                            :: n_p_row, n_p_col
      INTEGER                                            :: n_t_row, n_t_col
      INTEGER                                            :: p_row_i, p_col_i
      INTEGER                                            :: p_row_0, p_col_0
      INTEGER, DIMENSION(9)                              :: desc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: v1, v2
      REAL(KIND=dp), DIMENSION(:, :), POINTER             :: d1, d2
      TYPE(cp_fm_type), POINTER                          :: m1, m2

      INTEGER, EXTERNAL                                  :: numroc

      NULLIFY (m1)
      NULLIFY (m2)
      NULLIFY (d1)
      NULLIFY (d2)
      NULLIFY (v1)
      NULLIFY (v2)

      tg1 = 0
      tg2 = 0

#endif

      ! Initialization
      NULLIFY (et_proj_sec)
      NULLIFY (logger)

      ! Check whether to print coupling elements or not
      logger => cp_get_default_logger()
      et_proj_sec => section_vals_get_subs_vals(qs_env%input, &
                                                'PROPERTIES%ET_COUPLING%PROJECTION')
      print_key => section_vals_get_subs_vals(et_proj_sec, &
                                              'PRINT%COUPLING_ELEMENTS')
      do_print = BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)

      IF (do_print) THEN

         ! Header
         IF (output_unit > 0) &
            WRITE (output_unit, '(/,T3,A/)') 'Coupling elements [meV]:'

#if defined(__SCALAPACK)

         DO i = 1, ec%n_blocks
            id = 1
            DO j = i + 1, ec%n_blocks

               ! Matrix pointers
               m1 => ec%block(i)%hab(1, id)%matrix
               IF (SIZE(ec%block(i)%hab, 1) > 1) &
                  m2 => ec%block(i)%hab(2, id)%matrix
               ! Total number of rows/columns
               n_t_row = m1%matrix_struct%nrow_global
               n_t_col = m1%matrix_struct%ncol_global
               ! Number of rows/columns in Scalapack blocks
               n_b_row = m1%matrix_struct%nrow_block
               n_b_col = m1%matrix_struct%ncol_block

               ! Context and descriptor
               cntxt_loc = m1%matrix_struct%para_env%group
               CALL cp_blacs_gridinit(cntxt_loc, 'R', 1, m1%matrix_struct%para_env%num_pe)
               CALL cp_blacs_gridinfo(cntxt_loc, n_p_row, n_p_col, p_row_i, p_col_i)
               CALL descinit(desc, n_t_row, n_t_col, n_b_row, n_t_col, 0, 0, cntxt_loc, n_t_row, info)
               CPASSERT(info == 0)

               ! Number of rows on the actual processor
               n = numroc(n_t_row, n_b_row, p_row_i, 0, n_p_row)

               ! Local data arrays
               ALLOCATE (d1(MAX(1, n), n_t_col))
               ALLOCATE (v1(n_b_row*n_t_col))
               v1 = HUGE(1.0_dp)
               IF (SIZE(ec%block(i)%hab, 1) > 1) THEN
                  ALLOCATE (d2(MAX(1, n), n_t_col))
                  ALLOCATE (v2(n_b_row*n_t_col))
                  v2 = HUGE(1.0_dp)
               END IF

               ! Get the matrix elements
               CALL pdgemr2d(n_t_row, n_t_col, m1%local_data(1, 1), 1, 1, &
                             m1%matrix_struct%descriptor, &
                             d1(1, 1), 1, 1, desc, cntxt_loc)
               IF (SIZE(ec%block(i)%hab, 1) > 1) THEN
                  CALL pdgemr2d(n_t_row, n_t_col, m2%local_data(1, 1), 1, 1, &
                                m2%matrix_struct%descriptor, &
                                d2(1, 1), 1, 1, desc, cntxt_loc)
               END IF

               ! Row / column index
               k = 1
               l = 1

               ! Loop over blocks of rows
               DO ir = 1, n_t_row, MAX(n_b_row, 1)

                  ! Number of rows in the block
                  ib = MIN(n_b_row, n_t_row - ir + 1)

                  ! Local element position
                  CALL infog2l(ir, 1, desc, n_p_row, n_p_col, p_row_i, p_col_i, &
                               ir_loc, ic_loc, p_row_0, p_col_0)

                  ! Save block-row columns
                  IF (p_row_0 == p_row_i) THEN
                     DO ix = 1, ib
                        v1((ix - 1)*n_t_col + 1:n_t_col*ix) = d1(ir_loc + ix - 1, :)
                        IF (SIZE(ec%block(i)%hab, 1) > 1) &
                           v2((ix - 1)*n_t_col + 1:n_t_col*ix) = d2(ir_loc + ix - 1, :)
                     END DO
                  END IF

                  ! Transfer data
                  IF (p_row_0 /= 0) THEN
                     IF (p_row_0 == p_row_i) THEN
                        CALL mp_send(v1(:), 0, tg1, m1%matrix_struct%para_env%group)
                        IF (SIZE(ec%block(i)%hab, 1) > 1) &
                           CALL mp_send(v2(:), 0, tg2, m2%matrix_struct%para_env%group)
                     END IF
                     IF (p_row_i == 0) THEN
                        CALL mp_recv(v1(:), p_row_0, tg1, m1%matrix_struct%para_env%group)
                        IF (SIZE(ec%block(i)%hab, 1) > 1) &
                           CALL mp_recv(v2(:), p_row_0, tg2, m2%matrix_struct%para_env%group)
                     END IF
                  END IF

                  ! Output
                  IF (output_unit > 0) THEN
                     DO ix = 1, ib
                        DO ic = (ix - 1)*n_t_col + 1, n_t_col*ix

                           IF (k <= MIN(ec%block(i)%n_ao, n_states) .AND. &
                               l <= MIN(ec%block(j)%n_ao, n_states)) THEN

                              WRITE (output_unit, '(T3,I3,A,I4,A,I1,A,I4,A,E20.6)', ADVANCE='no') &
                                 i, "[", k, "] - ", j, "[", l, "] ", &
                                 v1(ic)*evolt*1000.0_dp

                              ! Second spin component
                              IF (SIZE(ec%block(i)%hab, 1) > 1 .AND. &
                                  k <= n_beta_states .AND. l <= n_beta_states) THEN
                                 WRITE (output_unit, '(E20.6)', ADVANCE='no') &
                                    v2(ic)*evolt*1000.0_dp
                              END IF

                              WRITE (output_unit, *)

                           END IF

                           ! Row / column IDs
                           l = l + 1
                           IF (l > n_t_col) THEN
                              l = 1
                              k = k + 1
                           END IF

                        END DO
                     END DO
                  END IF

               END DO

               ! clean memory
               CALL cp_blacs_gridexit(cntxt_loc)

               IF (ASSOCIATED(d1)) &
                  DEALLOCATE (d1)
               IF (ASSOCIATED(v1)) &
                  DEALLOCATE (v1)
               IF (SIZE(ec%block(i)%hab, 1) > 1) THEN
                  IF (ASSOCIATED(d2)) &
                     DEALLOCATE (d2)
                  IF (ASSOCIATED(v2)) &
                     DEALLOCATE (v2)
               END IF

               id = id + 1

            END DO
         END DO

#else

         IF (output_unit > 0) THEN

            DO i = 1, ec%n_blocks
               id = 1
               DO j = i + 1, ec%n_blocks

                  DO k = 1, MIN(ec%block(i)%n_ao, n_states)
                     DO l = 1, MIN(ec%block(j)%n_ao, n_states)

                        IF (SIZE(ec%block(i)%hab, 1) > 1) THEN

                           WRITE (output_unit, '(T3,I3,A,I4,A,I1,A,I4,A,E20.6)', ADVANCE='no') &
                              i, "[", k, "] - ", j, "[", l, "] ", &
                              ec%block(i)%hab(1, id)%matrix%local_data(k, l)*evolt*1000.0_dp
                           IF ((k <= n_beta_states) .AND. (l <= n_beta_states)) THEN
                              WRITE (output_unit, '(E20.6)') &
                                 ec%block(i)%hab(2, id)%matrix%local_data(k, l)*evolt*1000.0_dp
                           ELSE
                              WRITE (output_unit, *)
                           END IF

                        ELSE

                           WRITE (output_unit, '(T3,I3,A,I4,A,I1,A,I4,A,E20.6)') &
                              i, "[", k, "] - ", j, "[", l, "] ", &
                              ec%block(i)%hab(1, id)%matrix%local_data(k, l)*evolt*1000.0_dp

                        END IF

                     END DO
                  END DO

                  id = id + 1

               END DO
            END DO

         END IF

#endif

      END IF

   END SUBROUTINE pod_print_couplings

END MODULE et_coupling_pod_print
