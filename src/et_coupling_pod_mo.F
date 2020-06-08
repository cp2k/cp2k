!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines dealing with molecular orbitals localized by POD method
!> \author Z. Futera (03.2019)
! **************************************************************************************************
MODULE et_coupling_pod_mo

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type,&
                                              cp_fm_vectorssum
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE et_coupling_pod_types,           ONLY: et_cpl,&
                                              et_cpl_atom
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              mo_set_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod_mo'

   PUBLIC :: pod_mo_set, &
             pod_mo_save, &
             pod_mo_sum_c2

CONTAINS

! **************************************************************************************************
!> \brief Normalize set of MO vectors
!> \param qs_env QuickStep environment containing all system data
!> \param mo storage for the MO data set
!> \param n_ao number of AO basis functions
!> \param n_mo number of block states
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_mo_norm(qs_env, mo, n_ao, n_mo)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_type), POINTER                         :: mo
      INTEGER                                            :: n_ao, n_mo

      CHARACTER(len=*), PARAMETER :: routineN = 'pod_mo_norm', routineP = moduleN//':'//routineN

      REAL(KIND=dp), DIMENSION(:), POINTER               :: vec_t
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      TYPE(cp_fm_type), POINTER                          :: mat_sc, mat_t
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_s

! Routine name for debug purposes
! Local variables

! Routine name for debug purposes

      ! Initialization
      NULLIFY (blacs_env)
      NULLIFY (para_env)
      NULLIFY (fm_s)
      NULLIFY (mat_s)
      NULLIFY (mat_sc)
      NULLIFY (mat_t)
      NULLIFY (vec_t)

      ! Overlap matrix
      CALL get_qs_env(qs_env, matrix_s=mat_s)

      ! Calculate S*C product
      CALL cp_fm_create(matrix=mat_sc, matrix_struct=mo%mo_coeff%matrix_struct, &
                        name='S*C PRODUCT MATRIX')
      CALL cp_dbcsr_sm_fm_multiply(mat_s(1)%matrix, mo%mo_coeff, mat_sc, n_mo)

      ! Calculate C^T*S*C
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                               nrow_global=n_mo, ncol_global=n_mo)
      CALL cp_fm_create(matrix=mat_t, matrix_struct=fm_s, &
                        name='C^T*S*C OVERLAP PRODUCT MATRIX')
      CALL cp_gemm('T', 'N', n_mo, n_mo, n_ao, 1.0_dp, mo%mo_coeff, mat_sc, 0.0_dp, mat_t)

      ! Normalization
      ALLOCATE (vec_t(n_mo))
      CPASSERT(ASSOCIATED(vec_t))
      CALL cp_fm_vectorssum(mat_t, vec_t)
      vec_t = 1.0_dp/DSQRT(vec_t)
      CALL cp_fm_column_scale(mo%mo_coeff, vec_t)

      ! Clean memory
      CALL cp_fm_struct_release(fmstruct=fm_s)
      CALL cp_fm_release(matrix=mat_sc)
      CALL cp_fm_release(matrix=mat_t)
      IF (ASSOCIATED(vec_t)) &
         DEALLOCATE (vec_t)

   END SUBROUTINE pod_mo_norm

! **************************************************************************************************
!> \brief Transform block MO coefficients to original non-orthogonal basis set and save them
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param id block ID
!> \param mo storage for the MO data set
!> \param mat_u matrix of the block states
!> \param n_ao number of AO basis functions
!> \param n_mo number of block states
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_mo_set_coeffs(qs_env, ec, id, mo, mat_u, n_ao, n_mo)

      IMPLICIT NONE

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: id
      TYPE(mo_set_type), POINTER                         :: mo
      TYPE(cp_fm_type), POINTER                          :: mat_u
      INTEGER                                            :: n_ao
      INTEGER                                            :: n_mo

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_mo_set_coeffs', &
                                     routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: ir, ic, jr, jc, nr, nc
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), POINTER                          :: mat_w
      TYPE(cp_para_env_type), POINTER                    :: para_env
#if defined(__SCALAPACK)
      INTEGER                                            :: is_p_row, is_p_col
      INTEGER                                            :: iu_p_row, iu_p_col
      INTEGER                                            :: ns_p_rows, ns_p_cols
      INTEGER                                            :: nu_p_rows, nu_p_cols
      INTEGER                                            :: cs_row, cs_col, cu_row, cu_col
      INTEGER                                            :: cs_row_f, cs_col_f, cu_row_f, cu_col_f
      REAL(KIND=dp)                                      :: xu

      ! number of processors
      ns_p_rows = mo%mo_coeff%matrix_struct%context%num_pe(1)
      ns_p_cols = mo%mo_coeff%matrix_struct%context%num_pe(2)
      nu_p_rows = mat_u%matrix_struct%context%num_pe(1)
      nu_p_cols = mat_u%matrix_struct%context%num_pe(2)

      ! position of processors
      cs_row = mo%mo_coeff%matrix_struct%context%mepos(1)
      cs_col = mo%mo_coeff%matrix_struct%context%mepos(2)
      cu_row = mat_u%matrix_struct%context%mepos(1)
      cu_col = mat_u%matrix_struct%context%mepos(2)
#endif

      NULLIFY (blacs_env)
      NULLIFY (para_env)
      NULLIFY (mat_w)

      ! Working matrix
      CALL cp_fm_create(matrix=mat_w, matrix_struct=mo%mo_coeff%matrix_struct, &
                        name='BLOCK MO-TRANSFORMATION WORKING MATRIX')
      CALL cp_fm_set_all(mat_w, 0.0_dp)

      ! Matrix-element reordering
      nr = 1
      ! Rows
      DO ir = 1, ec%block(id)%n_atoms
         DO jr = 1, ec%block(id)%atom(ir)%n_ao
            ! Columns
            nc = 1
            DO ic = 1, ec%block(id)%n_atoms
               DO jc = 1, ec%block(id)%atom(ic)%n_ao
#if defined(__SCALAPACK)
                  CALL infog2l(ec%block(id)%atom(ir)%ao_pos + jr - 1, nc, &
                               mat_w%matrix_struct%descriptor, &
                               ns_p_rows, ns_p_cols, cs_row, cs_col, &
                               is_p_row, is_p_col, cs_row_f, cs_col_f)
                  CALL infog2l(nr, nc, mat_u%matrix_struct%descriptor, &
                               nu_p_rows, nu_p_cols, cu_row, cu_col, &
                               iu_p_row, iu_p_col, cu_row_f, cu_col_f)
                  ! Local U element
                  IF ((cu_row_f == cu_row) .AND. (cu_col_f == cu_col)) THEN
                     xu = mat_u%local_data(iu_p_row, iu_p_col)
                     CALL dgebs2d(mat_u%matrix_struct%context%group, 'All', ' ', 1, 1, xu, 1)
                     IF ((cs_row_f == cs_row) .AND. (cs_col_f == cs_col)) THEN
                        mat_w%local_data(is_p_row, is_p_col) = xu
                     END IF
                     ! Remote U element
                  ELSE
                     CALL dgebr2d(mat_u%matrix_struct%context%group, &
                                  'All', ' ', 1, 1, xu, 1, cu_row_f, cu_col_f)
                     IF ((cs_row_f == cs_row) .AND. (cs_col_f == cs_col)) THEN
                        mat_w%local_data(is_p_row, is_p_col) = xu
                     END IF
                  END IF
#else
                  mat_w%local_data(ec%block(id)%atom(ir)%ao_pos + jr - 1, nc) = &
                     mat_u%local_data(nr, nc)
#endif
                  nc = nc + 1
               END DO
            END DO
            nr = nr + 1
         END DO
      END DO

      ! Transformation to original non-orthogonal basis set
      CALL cp_gemm("N", "N", n_ao, n_mo, n_ao, 1.0_dp, ec%tmatf, mat_w, 0.0_dp, mo%mo_coeff)
      CALL pod_mo_norm(qs_env, mo, n_ao, n_mo)

      ! Clean memory
      CALL cp_fm_release(matrix=mat_w)

   END SUBROUTINE pod_mo_set_coeffs

! **************************************************************************************************
!> \brief Creates MO set corresponding to one atomic data block
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param id block ID
!> \param spin spin component
!> \param mat_u matrix of the block states
!> \param vec_e array of the block eigenvalues
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_mo_set(qs_env, ec, id, spin, mat_u, vec_e)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      INTEGER                                            :: id, spin
      TYPE(cp_fm_type), POINTER                          :: mat_u
      REAL(KIND=dp), DIMENSION(:), POINTER               :: vec_e

      CHARACTER(len=*), PARAMETER :: routineN = 'pod_mo_set', routineP = moduleN//':'//routineN

      INTEGER                                            :: n_ao, n_el, n_mo
      REAL(KIND=dp)                                      :: mx_occ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_s
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_cntrl
      TYPE(mo_set_type), POINTER                         :: mo
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_cntrl

! Routine name for debug purposes
! Local variables

      NULLIFY (blacs_env)
      NULLIFY (fm_s)
      NULLIFY (para_env)
      NULLIFY (dft_cntrl)
      NULLIFY (mo)
      NULLIFY (qs_kind_set)
      NULLIFY (scf_cntrl)

      ! Number of basis functions
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao)

      ! Number of states
      n_mo = mat_u%matrix_struct%nrow_global
      IF (n_mo /= mat_u%matrix_struct%ncol_global) &
         CPABORT('block state matrix is not square')
      IF (n_mo /= SIZE(vec_e)) &
         CPABORT('inconsistent number of states / energies')

      ! Maximal occupancy
      CALL get_qs_env(qs_env, dft_control=dft_cntrl)
      mx_occ = 2.0_dp
      IF (dft_cntrl%nspins > 1) &
         mx_occ = 1.0_dp

      ! Number of electrons
      n_el = ec%block(id)%n_electrons
      IF (dft_cntrl%nspins > 1) THEN
         n_el = n_el/2
         IF (MOD(ec%block(id)%n_electrons, 2) == 1) THEN
            IF (spin == 1) &
               n_el = n_el + 1
         END IF
      END IF

      ! Memory allocation
      NULLIFY (ec%block(id)%mo(spin)%mo_set)
      CALL allocate_mo_set(ec%block(id)%mo(spin)%mo_set, n_ao, n_mo, n_el, REAL(n_el, dp), mx_occ, 0.0_dp)
      mo => ec%block(id)%mo(spin)%mo_set

      ! State energies
      NULLIFY (mo%eigenvalues)
      ALLOCATE (mo%eigenvalues(n_mo))
      CPASSERT(ASSOCIATED(mo%eigenvalues))
      mo%eigenvalues = vec_e

      ! States coefficients
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL cp_fm_struct_create(fmstruct=fm_s, para_env=para_env, context=blacs_env, &
                               nrow_global=n_ao, ncol_global=n_mo)
      CALL cp_fm_create(matrix=mo%mo_coeff, matrix_struct=fm_s, name='BLOCK STATES')

      ! Transform MO coefficients to original non-orthogonal basis set
      CALL pod_mo_set_coeffs(qs_env, ec, id, mo, mat_u, n_ao, n_mo)

      ! Occupancies
      NULLIFY (mo%occupation_numbers)
      ALLOCATE (mo%occupation_numbers(n_mo))
      CPASSERT(ASSOCIATED(mo%occupation_numbers))
      mo%occupation_numbers = 0.0_dp

      IF (n_el > 0) THEN
         CALL get_qs_env(qs_env, scf_control=scf_cntrl)
         CALL set_mo_occupation(mo_set=mo, smear=scf_cntrl%smear)
      END IF

      ! Clean memory
      CALL cp_fm_struct_release(fmstruct=fm_s)

   END SUBROUTINE pod_mo_set

! **************************************************************************************************
!> \brief save given electronic state to cube files
!> \param qs_env QuickStep environment containing all system data
!> \param logger output logger
!> \param input input-file block print setting section
!> \param mo electronic states data
!> \param ib block ID
!> \param im state ID
!> \param is spin ID
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_mo_save(qs_env, logger, input, mo, ib, im, is)

      ! Routine arguments
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(mo_set_type), POINTER                         :: mo
      INTEGER                                            :: ib, im, is

      CHARACTER(len=*), PARAMETER :: routineN = 'pod_mo_save', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename
      CHARACTER(LEN=default_string_length)               :: title
      INTEGER                                            :: unit_nr
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: wf_g, wf_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys

! Routine name for debug purposes
! Local variables

      ! Initialization
      NULLIFY (particles)
      NULLIFY (subsys)

      NULLIFY (pw_env)
      NULLIFY (pw_pools)
      NULLIFY (auxbas_pw_pool)

      NULLIFY (atomic_kind_set)
      NULLIFY (cell)
      NULLIFY (dft_control)
      NULLIFY (particle_set)
      NULLIFY (qs_kind_set)

      ! Name of the cube file
      WRITE (filename, '(A7,I1.1,A1,I5.5,A1,I1.1)') 'POD-MO_', ib, '_', im, '_', is
      ! Open the file
      unit_nr = cp_print_key_unit_nr(logger, input, 'MO_CUBES', extension='.cube', &
                                     middle_name=TRIM(filename), file_position='REWIND', &
                                     log_filename=.FALSE.)
      ! Title of the file
      WRITE (title, *) 'WAVEFUNCTION ', im, ' block ', ib, ' spin ', is

      ! List of all atoms
      CALL get_qs_env(qs_env, subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)

      ! Grids for wavefunction
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_g%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      ! Calculate the grid values
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, dft_control=dft_control, particle_set=particle_set)
      CALL calculate_wavefunction(mo%mo_coeff, im, wf_r, wf_g, atomic_kind_set, &
                                  qs_kind_set, cell, dft_control, particle_set, pw_env)
      CALL cp_pw_to_cube(wf_r%pw, unit_nr, title, particles=particles, &
                         stride=section_get_ivals(input, 'MO_CUBES%STRIDE'))

      ! Close file
      CALL cp_print_key_finished_output(unit_nr, logger, input, 'MO_CUBES')

      ! Clean memory
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g%pw)

   END SUBROUTINE pod_mo_save

! **************************************************************************************************
!> \brief Return sum of selected squared MO coefficients
!> \param blk_at list of atoms in the block
!> \param mo array of MO sets
!> \param id state index
!> \param atom list of atoms for MO coefficient summing
!> \return ...
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   FUNCTION pod_mo_sum_c2(blk_at, mo, id, atom) RESULT(c2)

      IMPLICIT NONE

      ! Routine arguments
      TYPE(et_cpl_atom), DIMENSION(:), POINTER           :: blk_at
      TYPE(cp_fm_type), POINTER                          :: mo
      INTEGER                                            :: id
      INTEGER, DIMENSION(:), POINTER                     :: atom

      ! Returning value
      REAL(KIND=dp)                                      :: c2

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_mo_sum_c2', &
                                     routineP = moduleN//':'//routineN

      ! Local variables
      LOGICAL                                            :: found
      INTEGER                                            :: i, j, k
      REAL(KIND=dp)                                      :: c
#if defined(__SCALAPACK)
      INTEGER                                            :: c_row, c_col, c_row_f, c_col_f
      INTEGER                                            :: i_p_row, i_p_col
      INTEGER                                            :: n_p_rows, n_p_cols

      ! number of processors
      n_p_rows = mo%matrix_struct%context%num_pe(1)
      n_p_cols = mo%matrix_struct%context%num_pe(2)

      ! position of processors
      c_row = mo%matrix_struct%context%mepos(1)
      c_col = mo%matrix_struct%context%mepos(2)
#endif

      ! initialization
      c2 = 0.0d0

      ! selected atoms
      DO i = 1, SIZE(atom)

         ! find atomic function offset
         found = .FALSE.
         DO j = 1, SIZE(blk_at)
            IF (blk_at(j)%id == atom(i)) THEN
               found = .TRUE.
               EXIT
            END IF
         END DO

         IF (.NOT. found) &
            CPABORT('MO-fraction atom ID not defined in the block')

         ! sum MO coefficients from the atom
         DO k = 1, blk_at(j)%n_ao
#if defined(__SCALAPACK)
            c = 0.0_dp
            CALL infog2l(blk_at(j)%ao_pos + k - 1, id, mo%matrix_struct%descriptor, &
                         n_p_rows, n_p_cols, c_row, c_col, &
                         i_p_row, i_p_col, c_row_f, c_col_f)
            ! local element
            IF ((c_row_f == c_row) .AND. (c_col_f == c_col)) THEN
               c = mo%local_data(i_p_row, i_p_col)
               CALL dgebs2d(mo%matrix_struct%context%group, 'All', ' ', 1, 1, c, 1)
               ! remote element
            ELSE
               CALL dgebr2d(mo%matrix_struct%context%group, 'All', ' ', 1, 1, c, 1, c_row_f, c_col_f)
            END IF
#else
            c = mo%local_data(blk_at(j)%ao_pos + k - 1, id)
#endif
            c2 = c2 + c*c
         END DO

      END DO

   END FUNCTION pod_mo_sum_c2

END MODULE et_coupling_pod_mo
