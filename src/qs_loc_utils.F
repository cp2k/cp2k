!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Some utilities for the construction of
!>      the localization environment
!> \author MI (05-2005)
! *****************************************************************************
MODULE qs_loc_utils

  USE ai_moments,                      ONLY: contract_cossin,&
                                             cossin
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_1d_r_p_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_set,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
  USE cp_fm_diag,                      ONLY: choose_eigv_solver
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_element, cp_fm_get_info, cp_fm_get_submatrix, &
       cp_fm_maxabsval, cp_fm_p_type, cp_fm_release, cp_fm_set_all, &
       cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_type, cp_fm_write_unformatted
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE input_constants,                 ONLY: &
       do_loc_crazy, do_loc_direct, do_loc_jacobi, do_loc_l1_norm_sd, &
       do_loc_none, energy_loc_range, op_loc_berry, op_loc_boys, &
       op_loc_pipek, state_loc_all, state_loc_list, state_loc_none, &
       state_loc_range
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi,&
                                             twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                             localized_wfn_control_create,&
                                             localized_wfn_control_release,&
                                             localized_wfn_control_type,&
                                             qs_loc_env_new_type,&
                                             set_qs_loc_env
  USE qs_localization_methods,         ONLY: initialize_weights,&
                                             rotate_orbitals
  USE qs_mo_methods,                   ONLY: make_mo_eig
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_scf_types,                    ONLY: ot_method_nr
  USE scf_control_types,               ONLY: scf_control_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_loc_utils'

! *** Public ***
  PUBLIC :: jacobi_rotation_pipek, qs_loc_env_init, loc_write_restart,&
            retain_history, qs_loc_init,&
            set_loc_centers, set_loc_wfn_lists, qs_loc_control_init

CONTAINS

! *****************************************************************************
!> \brief copy old mos to new ones, allocating as necessary
!> \param mo_loc_history ...
!> \param mo_loc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE retain_history(mo_loc_history, mo_loc, error)

    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_loc_history, mo_loc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_history', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ncol_hist, ncol_loc
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.

    IF (.NOT. ASSOCIATED(mo_loc_history)) THEN
       ALLOCATE(mo_loc_history(SIZE(mo_loc)))
       DO i=1,SIZE(mo_loc_history)
          CALL cp_fm_create(mo_loc_history(i)%matrix, mo_loc(i)%matrix%matrix_struct,error=error)
       ENDDO
    ENDIF

    DO i=1,SIZE(mo_loc_history)
       CALL cp_fm_get_info(mo_loc_history(i)%matrix, ncol_global=ncol_hist, error=error)
       CALL cp_fm_get_info(mo_loc(i)%matrix,         ncol_global=ncol_loc, error=error)
       CPPrecondition(ncol_hist==ncol_loc,cp_failure_level,routineP,error,failure)
       CALL cp_fm_to_fm(mo_loc(i)%matrix, mo_loc_history(i)%matrix, error=error)
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE retain_history

! *****************************************************************************
!> \brief rotate the mo_new, so that the orbitals are as similar
!>        as possible to ones in mo_ref.
!> \param mo_new ...
!> \param mo_ref ...
!> \param matrix_S ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rotate_state_to_ref(mo_new, mo_ref, matrix_S, error)

    TYPE(cp_fm_type), POINTER                :: mo_new, mo_ref
    TYPE(cp_dbcsr_type), POINTER             :: matrix_S
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotate_state_to_ref', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ncol, ncol_ref, nrow
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: o1, o2, o3, o4, smo

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CALL cp_fm_get_info(mo_new, nrow_global=nrow, ncol_global=ncol, error=error)
    CALL cp_fm_get_info(mo_ref, ncol_global=ncol_ref, error=error)
    CPPrecondition(ncol==ncol_ref,cp_failure_level,routineP,error,failure)

    NULLIFY(fm_struct_tmp,o1,o2,o3,o4,smo)
    CALL cp_fm_create(smo,mo_ref%matrix_struct,error=error)

    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, &
            ncol_global=ncol, para_env=mo_new%matrix_struct%para_env, &
            context=mo_new%matrix_struct%context,error=error)
    CALL cp_fm_create(o1,fm_struct_tmp,error=error)
    CALL cp_fm_create(o2,fm_struct_tmp,error=error)
    CALL cp_fm_create(o3,fm_struct_tmp,error=error)
    CALL cp_fm_create(o4,fm_struct_tmp,error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    ! o1 = (mo_new)^T matrix_S mo_ref
    CALL cp_dbcsr_sm_fm_multiply(matrix_S,mo_ref,smo,ncol,error=error)
    CALL cp_gemm('T','N',ncol,ncol,nrow,1.0_dp, mo_new,smo,0.0_dp,o1,error=error)

    ! o2 = (o1^T o1)
    CALL cp_gemm('T','N',ncol,ncol,ncol,1.0_dp, o1,o1,0.0_dp,o2,error=error)

    ! o2 = (o1^T o1)^-1/2
    ALLOCATE(eigenvalues(ncol))
    CALL choose_eigv_solver(o2,o3,eigenvalues,error=error)
    CALL cp_fm_to_fm(o3,o4,error=error)
    eigenvalues(:) = 1.0_dp/SQRT(eigenvalues(:))
    CALL cp_fm_column_scale(o4,eigenvalues)
    CALL cp_gemm('N','T',ncol,ncol,ncol,1.0_dp,o3,o4,0.0_dp,o2,error=error)

    ! o3 = o1 (o1^T o1)^-1/2
    CALL cp_gemm('N','N',ncol,ncol,ncol,1.0_dp,o1,o2,0.0_dp,o3,error=error)

    ! mo_new o1 (o1^T o1)^-1/2
    CALL cp_gemm('N','N',nrow,ncol,ncol,1.0_dp,mo_new,o3,0.0_dp,smo,error=error)
    CALL cp_fm_to_fm(smo,mo_new,error=error)

    ! XXXXXXX testing
    ! CALL cp_gemm('N','T',ncol,ncol,ncol,1.0_dp,o3,o3,0.0_dp,o1,error=error)
    ! WRITE(*,*) o1%local_data
    ! CALL cp_gemm('T','N',ncol,ncol,ncol,1.0_dp,o3,o3,0.0_dp,o1,error=error)
    ! WRITE(*,*) o1%local_data

    CALL cp_fm_release(o1,error=error)
    CALL cp_fm_release(o2,error=error)
    CALL cp_fm_release(o3,error=error)
    CALL cp_fm_release(o4,error=error)
    CALL cp_fm_release(smo,error=error)

    CALL timestop(handle)

  END SUBROUTINE rotate_state_to_ref


! *****************************************************************************
!> \brief 2by2 rotation for the pipek operator
!>       in this case the z_ij numbers are reals
!> \param zij_fm_set ...
!> \param vectors ...
!> \param sweeps ...
!> \param error ...
!> \par History
!>       05-2005 created
!> \author MI
! *****************************************************************************
  SUBROUTINE jacobi_rotation_pipek(zij_fm_set, vectors, sweeps, error)

    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: zij_fm_set
    TYPE(cp_fm_type), POINTER                :: vectors
    INTEGER                                  :: sweeps
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'jacobi_rotation_pipek', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, istate, jstate, natom, &
                                                nstate
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aij, bij, ct, mii, mij, mjj, &
                                                ratio, st, theta, tolerance
    TYPE(cp_fm_type), POINTER                :: rmat

    failure = .FALSE.

    NULLIFY(rmat)

    CALL cp_fm_create(rmat,zij_fm_set(1,1)%matrix%matrix_struct,error=error)
    CALL cp_fm_set_all(rmat, 0.0_dp, 1.0_dp, error=error)

    CALL cp_fm_get_info ( rmat , nrow_global = nstate, error = error )
    tolerance = 1.0e10_dp
    sweeps = 0
    natom = SIZE(zij_fm_set,1)
    ! do jacobi sweeps until converged
    DO WHILE ( tolerance >= 1.0e-4_dp )
       sweeps = sweeps + 1
       DO istate = 1, nstate
         DO jstate = istate+1, nstate
           aij = 0.0_dp
           bij = 0.0_dp
           DO iatom = 1,natom
             CALL cp_fm_get_element(zij_fm_set(iatom,1)%matrix,istate,istate,mii)
             CALL cp_fm_get_element(zij_fm_set(iatom,1)%matrix,istate,jstate,mij)
             CALL cp_fm_get_element(zij_fm_set(iatom,1)%matrix,jstate,jstate,mjj)
             aij = aij +  mij*(mii - mjj)
             bij = bij +  mij*mij - 0.25_dp * (mii-mjj)*(mii-mjj)

           END DO
           IF (ABS (bij) > 1.E-10_dp) THEN
             ratio = -aij/bij
             theta = 0.25_dp * ATAN(ratio)
           ELSE
             bij = 0.0_dp
             theta = 0.0_dp
           END IF
           ! Check max or min
           ! To minimize the spread
           IF(theta > pi * 0.5_dp) THEN
             theta = theta - pi*0.25_dp
           ELSE IF(theta < - pi*0.5_dp) THEN
             theta =theta + pi*0.25_dp
           END IF

           ct = COS(theta)
           st = SIN(theta)
           CALL rotate_rmat_real(istate,jstate,st,ct,rmat,error=error)
           CALL rotate_zij_real(istate,jstate,st,ct,zij_fm_set,error=error)
         END DO
       END DO
       CALL check_tolerance_real(zij_fm_set,tolerance,error=error)
    END DO

    CALL rotate_orbitals ( rmat, vectors )
    CALL cp_fm_release ( rmat ,error=error)

  END SUBROUTINE jacobi_rotation_pipek

! *****************************************************************************
!> \brief ...
!> \param istate ...
!> \param jstate ...
!> \param st ...
!> \param ct ...
!> \param rmat ...
!> \param error ...
!> \par History
!>      04.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE rotate_rmat_real(istate,jstate,st,ct,rmat,error)

    INTEGER, INTENT(IN)                      :: istate, jstate
    REAL(dp), INTENT(IN)                     :: st, ct
    TYPE(cp_fm_type), POINTER                :: rmat
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: nstate
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desc
#endif
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL cp_fm_get_info(rmat, nrow_global=nstate, error=error)
#if defined(__SCALAPACK)
    desc(:) = rmat%matrix_struct%descriptor(:)
    CALL pzrot(nstate,rmat%local_data(1,1),1,istate,desc,1,&
         rmat%local_data(1,1),1,jstate,desc,1,ct,st)
#else
    CALL drot(nstate,rmat%local_data(1,istate),1,&
         rmat%local_data(1,jstate),1,ct,st)
#endif

  END SUBROUTINE rotate_rmat_real

! *****************************************************************************
!> \brief ...
!> \param istate ...
!> \param jstate ...
!> \param st ...
!> \param ct ...
!> \param zij ...
!> \param error ...
!> \par History
!>      04.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE rotate_zij_real(istate,jstate,st,ct,zij,error)

    INTEGER, INTENT(IN)                      :: istate, jstate
    REAL(dp), INTENT(IN)                     :: st, ct
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: zij
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: iatom, natom, nstate
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desc
#else
    INTEGER                                  :: stride
#endif
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL cp_fm_get_info(zij(1,1)%matrix,nrow_global = nstate, error=error)

    natom = SIZE(zij, 1)
    DO iatom = 1, natom
#if defined(__SCALAPACK)
      desc(:) = zij(iatom,1) % matrix %matrix_struct%descriptor(:)
      CALL pzrot(nstate,zij(iatom,1)%matrix%local_data(1,1),&
           1,istate,desc,1,zij(iatom,1)%matrix%local_data(1,1),&
           1,jstate,desc,nstate,ct,st)
      CALL pzrot(nstate,zij(iatom,1)%matrix%local_data(1,1),&
           istate,1,desc,1,zij(iatom,1)%matrix%local_data(1,1),&
           jstate,1,desc,nstate,ct,st)
#else
       CALL drot(nstate,zij(iatom,1)%matrix%local_data(1,istate),&
            1,zij(iatom,1)%matrix%local_data(1,jstate),1,ct,st)
       stride = SIZE(zij(iatom,1)%matrix%local_data,1)
       CALL drot(nstate,zij(iatom,1)%matrix%local_data(istate,1),&
            stride,zij(iatom,1)%matrix%local_data(jstate,1),stride,ct,st)
#endif
    END DO

  END SUBROUTINE rotate_zij_real

! *****************************************************************************
!> \brief ...
!> \param zij_fm_set ...
!> \param tolerance ...
!> \param error ...
!> \par History
!>      04.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE check_tolerance_real(zij_fm_set,tolerance,error)

    TYPE(cp_fm_p_type), DIMENSION(:, :)      :: zij_fm_set
    REAL(dp), INTENT(OUT)                    :: tolerance
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: iatom, istat, istate, jstate, &
                                                natom, ncol_local, &
                                                nrow_global, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    REAL(dp)                                 :: grad_ij, zii, zij, zjj
    REAL(dp), DIMENSION(:, :), POINTER       :: diag
    TYPE(cp_fm_type), POINTER                :: force

    failure = .FALSE.

    NULLIFY( force , diag, col_indices, row_indices )
    natom = SIZE(zij_fm_set,1)
    CALL cp_fm_get_info ( zij_fm_set(1,1)%matrix, nrow_local = nrow_local,  &
         ncol_local = ncol_local, nrow_global = nrow_global, &
         row_indices = row_indices, col_indices = col_indices, error = error )
    ALLOCATE(diag(nrow_global,natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,"check_tolerance_real",error,failure)

    CALL cp_fm_create(force , zij_fm_set(1,1)%matrix%matrix_struct,error=error)
    CALL cp_fm_set_all (force, 0._dp, error= error)

    DO iatom = 1,natom
      DO istate = 1,nrow_global
        CALL cp_fm_get_element(zij_fm_set(iatom,1)%matrix,istate,istate,diag(istate,iatom))
      END DO
    END DO

    DO istate = 1,nrow_local
      DO jstate = 1,ncol_local
        grad_ij = 0.0_dp
        DO iatom = 1,natom
          zii = diag(istate,iatom)
          zjj = diag(jstate,iatom)
          zij = zij_fm_set(iatom,1)%matrix%local_data(istate,jstate)
          grad_ij = grad_ij + 4.0_dp*zij*(zjj-zii)
        END DO
        force%local_data(istate,jstate) = grad_ij
      END DO
    END DO

    DEALLOCATE(diag,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,"check_tolerance_real",error,failure)

    CALL cp_fm_maxabsval ( force, tolerance, error = error )

    CALL cp_fm_release ( force ,error=error)

  END SUBROUTINE check_tolerance_real

! *****************************************************************************
!> \brief allocates the data, and initializes the operators
!> \param qs_loc_env new environment for the localization calculations
!> \param localized_wfn_control variables and directives for the localization
!> \param qs_env the qs_env in which the qs_env lives
!> \param myspin ...
!> \param do_localize ...
!> \param loc_coeff ...
!> \param mo_loc_history ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      04.2005 created [MI]
!> \author MI
!> \note
!>      similar to the old one, but not quite
! *****************************************************************************
  SUBROUTINE qs_loc_env_init(qs_loc_env,localized_wfn_control,qs_env,myspin,do_localize,&
             loc_coeff, mo_loc_history, error)

    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN), OPTIONAL            :: myspin
    LOGICAL, INTENT(IN), OPTIONAL            :: do_localize
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: loc_coeff, mo_loc_history
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_loc_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: dim_op, handle, i, iatom, imo, imoloc, ispin, istat, j, &
      l_spin, lb, my_nspins, nao, naosub, natoms, nmo, nmosub, nspins, &
      s_spin, ub
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: my_occ, occ_imo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupations
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: moloc_coeff
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mat_ptr, mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(mos,matrix_s,moloc_coeff,particle_set,para_env,cell,local_molecules,occupations,mat_ptr)
    IF (PRESENT(do_localize)) qs_loc_env%do_localize = do_localize
    IF (qs_loc_env%do_localize) THEN
       CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, cell=cell, &
            local_molecules=local_molecules,particle_set=particle_set, &
            para_env=para_env,mos=mos, error=error)
       nspins    = SIZE(mos,1)
       s_spin    = 1
       l_spin    = nspins
       my_nspins = nspins
       IF(PRESENT(myspin))  THEN
          s_spin = myspin
          l_spin = myspin
          my_nspins = 1
       END IF
       ALLOCATE (moloc_coeff(my_nspins), STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ispin = s_spin,l_spin
          NULLIFY(tmp_fm_struct,mo_coeff)
          CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,nao=nao,nmo=nmo)
          nmosub = localized_wfn_control%nloc_states(ispin)
          CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
               ncol_global=nmosub,para_env=para_env,context=mo_coeff%matrix_struct%context,error=error)
          CALL cp_fm_create (moloc_coeff(ispin)%matrix , tmp_fm_struct ,error=error)

          CALL cp_fm_get_info(moloc_coeff(ispin)%matrix, nrow_global=naosub, &
               ncol_global=nmosub,error=error)
          CPPostcondition(nao==naosub,cp_failure_level,routineP,error,failure)
          IF(localized_wfn_control%do_homo)THEN
             CPPostcondition(nmo>=nmosub,cp_failure_level,routineP,error,failure)
          ELSE
             CPPostcondition(nao-nmo>=nmosub,cp_failure_level,routineP,error,failure)
          END IF
          CALL cp_fm_set_all(moloc_coeff(ispin)%matrix,0.0_dp,error=error)
          CALL cp_fm_struct_release ( tmp_fm_struct ,error=error)
       END DO  ! ispin
       ! Copy the submatrix

       DO ispin =  s_spin,l_spin
          CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,&
               occupation_numbers=occupations, nao=nao, nmo=nmo)
          lb = localized_wfn_control%lu_bound_states(1,ispin)
          ub = localized_wfn_control%lu_bound_states(2,ispin)

          IF(PRESENT(loc_coeff)) THEN
            mat_ptr=> loc_coeff(ispin)%matrix
          ELSE
            mat_ptr=> mo_coeff
          END IF
          IF(localized_wfn_control%set_of_states == state_loc_list) THEN
             ALLOCATE(vecbuffer(1,nao),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             IF(localized_wfn_control%do_homo)THEN
                 my_occ =  occupations(localized_wfn_control%loc_states(1, ispin))
             END IF
             nmosub = SIZE(localized_wfn_control%loc_states,1)
             CPPostcondition(nmosub>0,cp_failure_level,routineP,error,failure)
             imoloc = 0
             DO i = lb,ub
                ! Get the index in the subset
                imoloc = imoloc + 1
                ! Get the index in the full set
                imo = localized_wfn_control%loc_states(i,ispin)
                IF(localized_wfn_control%do_homo)THEN
                   occ_imo = occupations(imo)
                   IF(ABS(occ_imo-my_occ) > localized_wfn_control% eps_occ) THEN
                      IF(localized_wfn_control%localization_method /= do_loc_none) &
                           CALL stop_program(routineN,moduleN,__LINE__,&
                                             "States with different occupations "//&
                                             "cannot be rotated together")
                   END IF
                END IF
                ! Take the imo vector from the full set and copy in the imoloc vector of the subset
                CALL cp_fm_get_submatrix(mat_ptr,vecbuffer,1,imo,&
                     nao,1,transpose=.TRUE.,error=error)
                CALL cp_fm_set_submatrix(moloc_coeff(ispin)%matrix,vecbuffer,1,imoloc,&
                     nao,1,transpose=.TRUE.,error=error)
             END DO
             DEALLOCATE(vecbuffer,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ELSE
             my_occ  = occupations( lb )
             occ_imo = occupations( ub )
             IF(ABS(occ_imo-my_occ) >  localized_wfn_control% eps_occ) THEN
                IF(localized_wfn_control%localization_method /= do_loc_none) &
                     CALL stop_program(routineN,moduleN,__LINE__,&
                                       "States with different occupations "//&
                                       "cannot be rotated together")
             END IF
             nmosub = localized_wfn_control%nloc_states(ispin)
             CALL cp_fm_to_fm(mat_ptr,moloc_coeff(ispin)%matrix,nmosub,lb,1)
          END IF

          ! we have the mo's to be localized now, see if we can rotate them according to the history
          ! only do that if we have a history of course. The history is filled
          IF (PRESENT(mo_loc_history)) THEN
            IF (localized_wfn_control%use_history .AND. ASSOCIATED(mo_loc_history)) THEN
               CALL rotate_state_to_ref(moloc_coeff(ispin)%matrix, &
                   mo_loc_history(ispin)%matrix, matrix_s(1)%matrix, &
                   error=error)
            ENDIF
          ENDIF

       END DO

       CALL set_qs_loc_env(qs_loc_env=qs_loc_env, cell=cell, local_molecules=local_molecules,&
            moloc_coeff=moloc_coeff, particle_set=particle_set,para_env=para_env,&
            localized_wfn_control=localized_wfn_control,error=error)

       ! Prepare the operators
       NULLIFY(tmp_fm_struct,mo_coeff)
       nmosub = MAXVAL(localized_wfn_control%nloc_states)
       CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmosub,&
            ncol_global=nmosub,para_env=para_env,context=mo_coeff%matrix_struct%context,error=error)

       IF ( localized_wfn_control%operator_type==op_loc_berry )  THEN
          IF(qs_loc_env%cell%orthorhombic) THEN
             dim_op = 3
          ELSE
             dim_op = 6
          END IF
          CALL set_qs_loc_env(qs_loc_env=qs_loc_env, dim_op=dim_op, error=error)
          ALLOCATE ( qs_loc_env % op_sm_set ( 2, dim_op ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO i=1, dim_op
             DO j=1, SIZE ( qs_loc_env % op_sm_set, 1 )
                NULLIFY(qs_loc_env%op_sm_set(j,i)%matrix)
                ALLOCATE(qs_loc_env%op_sm_set(j,i)%matrix)
                CALL cp_dbcsr_init(qs_loc_env%op_sm_set(j,i)%matrix, error=error)
                CALL cp_dbcsr_copy(qs_loc_env%op_sm_set(j,i)%matrix, matrix_s(1)%matrix,&
                     name="qs_loc_env%op_sm_"//TRIM(ADJUSTL(cp_to_string(j)))//"-"//TRIM(ADJUSTL(cp_to_string(i))),&
                     error=error)
                CALL cp_dbcsr_set(qs_loc_env%op_sm_set(j,i)%matrix,0.0_dp,error=error)
             ENDDO
          END DO

       ELSEIF ( localized_wfn_control%operator_type== op_loc_pipek) THEN
          natoms = SIZE ( qs_loc_env % particle_set, 1 )
          ALLOCATE ( qs_loc_env % op_fm_set ( natoms , 1 ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          CALL set_qs_loc_env(qs_loc_env=qs_loc_env, dim_op=natoms, error=error)
          DO ispin = 1, SIZE(qs_loc_env % op_fm_set ,2)
             CALL get_mo_set(mos(ispin)%mo_set,nmo=nmo)
             DO iatom = 1,natoms
                NULLIFY(qs_loc_env %op_fm_set(iatom,ispin)%matrix)

                CALL cp_fm_create (qs_loc_env % op_fm_set(iatom,ispin)%matrix , tmp_fm_struct ,error=error)

                CALL cp_fm_get_info(qs_loc_env %op_fm_set(iatom,ispin)%matrix, nrow_global=nmosub, &
                     error=error)
                CPPostcondition(nmo>=nmosub,cp_failure_level,routineP,error,failure)
                CALL cp_fm_set_all(qs_loc_env %op_fm_set(iatom,ispin)%matrix,0.0_dp,error=error)
             END DO  ! iatom
          END DO  ! ispin
       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Type of operator not implemented")
       ENDIF
       CALL cp_fm_struct_release ( tmp_fm_struct ,error=error)

       IF (  localized_wfn_control%operator_type==op_loc_berry ) THEN

          CALL initialize_weights ( qs_loc_env % cell, qs_loc_env % weights )

          CALL get_berry_operator ( qs_loc_env, qs_env , error=error)

       ELSEIF ( localized_wfn_control%operator_type== op_loc_pipek) THEN

          !!    here we don't have to do anything
          !!    CALL get_pipek_mezey_operator ( qs_loc_env, qs_env )

       ENDIF

       qs_loc_env%molecular_states = .FALSE.
       qs_loc_env%wannier_states = .FALSE.
    END IF
    CALL timestop(handle)

  END SUBROUTINE qs_loc_env_init

! *****************************************************************************
!> \brief Computes the Berry operator for periodic systems
!>       used to define the spread of the MOS
!>       Here the matrix elements of the type <mu|cos(kr)|nu> and  <mu|sin(kr)|nu>
!>       are computed, where mu and nu are the contracted basis functions.
!>       Namely the Berry operator is exp(ikr)
!>       k is defined somewhere
!>       the pair lists are exploited and sparse matrixes are constructed
!> \param qs_loc_env new environment for the localization calculations
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      04.2005 created [MI]
!> \author MI
!> \note
!>      The intgrals are computed analytically  using the primitives GTO
!>      The contraction is performed block-wise
! *****************************************************************************
  SUBROUTINE get_berry_operator ( qs_loc_env, qs_env, error )
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_berry_operator', &
      routineP = moduleN//':'//routineN

    INTEGER :: dim_op, handle, i, iatom, icol, ikind, inode, irow, iset, &
      istat, jatom, jkind, jset, last_jatom, ldab, ldsa, ldsb, ldwork, maxl, &
      ncoa, ncob, nkind, nrow, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(3)                    :: perd0
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found, new_atom_b
    REAL(KIND=dp)                            :: dab, kvec(3), rab2, &
                                                vector_k(3,6)
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cosab, rpgfa, rpgfb, sinab, &
                                                sphi_a, sphi_b, work, zeta, &
                                                zetb
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: op_cos, op_sin
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: op_sm_set
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b, &
                                                orb_basis_set
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY (qs_kind, qs_kind_set)
    NULLIFY ( cell, op_sm_set , orb_basis_set, particle_set)
    NULLIFY (sab_orb)
    NULLIFY (cosab, sinab, work)
    NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )

    CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set,&
                    particle_set=particle_set, sab_orb=sab_orb,&
                    error=error)

    nkind = SIZE(qs_kind_set)

    CALL get_qs_kind_set(qs_kind_set=qs_kind_set,&
                             maxco=ldwork, maxlgto=maxl )
    ldab = ldwork
    ALLOCATE(cosab(ldab,ldab))
    cosab = 0.0_dp
    ALLOCATE(sinab(ldab,ldab))
    sinab = 0.0_dp
    ALLOCATE(work(ldwork,ldwork))
    work = 0.0_dp

    CALL get_qs_loc_env(qs_loc_env=qs_loc_env, op_sm_set=op_sm_set, &
         cell=cell, dim_op=dim_op, error=error)

    ALLOCATE (op_cos(dim_op),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (op_sin(dim_op),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, dim_op
       NULLIFY (op_cos(i)%block)
       NULLIFY (op_sin(i)%block)
    END DO

    kvec = 0.0_dp
    vector_k = 0.0_dp
    vector_k(:,1) = twopi*cell%h_inv(1,:)
    vector_k(:,2) = twopi*cell%h_inv(2,:)
    vector_k(:,3) = twopi*cell%h_inv(3,:)
    vector_k(:,4) = twopi*(cell%h_inv(1,:)+cell%h_inv(2,:))
    vector_k(:,5) = twopi*(cell%h_inv(1,:)+cell%h_inv(3,:))
    vector_k(:,6) = twopi*(cell%h_inv(2,:)+cell%h_inv(3,:))


    ! This operator can be used only for periodic systems
    ! If an isolated system is used the periodicity is overimposed
    perd0(1:3) = cell%perd(1:3)
    cell%perd(1:3) = 1

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      qs_kind => qs_kind_set(ikind)
      CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       ra = pbc(particle_set(iatom)%r,cell)
       ! basis ikind
       first_sgfa   =>  basis_set_a%first_sgf
       la_max       =>  basis_set_a%lmax
       la_min       =>  basis_set_a%lmin
       npgfa        =>  basis_set_a%npgf
       nseta        =   basis_set_a%nset
       nsgfa        =>  basis_set_a%nsgf_set
       rpgfa        =>  basis_set_a%pgf_radius
       set_radius_a =>  basis_set_a%set_radius
       sphi_a       =>  basis_set_a%sphi
       zeta         =>  basis_set_a%zet
       ! basis jkind
       first_sgfb   =>  basis_set_b%first_sgf
       lb_max       =>  basis_set_b%lmax
       lb_min       =>  basis_set_b%lmin
       npgfb        =>  basis_set_b%npgf
       nsetb        =   basis_set_b%nset
       nsgfb        =>  basis_set_b%nsgf_set
       rpgfb        =>  basis_set_b%pgf_radius
       set_radius_b =>  basis_set_b%set_radius
       sphi_b       =>  basis_set_b%sphi
       zetb         =>  basis_set_b%zet

       ldsa = SIZE(sphi_a,1)
       ldsb = SIZE(sphi_b,1)
       IF(inode==1) last_jatom = 0

       rb = rab + ra

       IF (jatom /= last_jatom) THEN
         new_atom_b = .TRUE.
         last_jatom = jatom
       ELSE
         new_atom_b = .FALSE.
       END IF

       IF (new_atom_b) THEN
         IF (iatom <= jatom) THEN
           irow = iatom
           icol = jatom
         ELSE
           irow = jatom
           icol = iatom
         END IF

         DO i = 1,dim_op
            NULLIFY(op_cos(i)%block)
            CALL cp_dbcsr_get_block_p(matrix=op_sm_set(1,i)%matrix,&
                 row=irow,col=icol,block=op_cos(i)%block,found=found)
            NULLIFY(op_sin(i)%block)
            CALL cp_dbcsr_get_block_p(matrix=op_sm_set(2,i)%matrix,&
                 row=irow,col=icol,block=op_sin(i)%block,found=found)
         END DO
       END IF  ! new_atom_b

       rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       dab = SQRT(rab2)

       nrow = 0
       DO iset=1,nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

         DO jset=1,nsetb

           ncob = npgfb(jset)*ncoset(lb_max(jset))
           sgfb = first_sgfb(1,jset)

           IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN

!           *** Calculate the primitive overlap integrals ***
             DO i = 1,dim_op

               kvec(1:3) = vector_k(1:3,i)
               cosab = 0.0_dp
               sinab = 0.0_dp
               CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),&
                           la_min(iset),lb_max(jset),npgfb(jset),zetb(:,jset),&
                           rpgfb(:,jset),lb_min(jset),&
                           ra,rb,kvec,cosab,sinab)

               CALL contract_cossin(op_cos(i)%block,op_sin(i)%block,&
                          iatom,ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                          jatom,ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                          cosab,sinab,ldab,work,ldwork)
             END DO

           END IF  !  >= dab

         END DO ! jset

         nrow = nrow + ncoa

       END DO ! iset

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    ! Set back the correct periodicity
    cell%perd(1:3) =  perd0(1:3)

    DO i = 1,dim_op
       NULLIFY(op_cos(i)%block)
       NULLIFY(op_sin(i)%block)
    END DO
    DEALLOCATE(op_cos,op_sin, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(cosab,sinab,work,basis_set_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)
  END SUBROUTINE get_berry_operator

! *****************************************************************************
!> \brief ...
!> \param qs_loc_env ...
!> \param section ...
!> \param mo_array ...
!> \param coeff_localized ...
!> \param do_homo ...
!> \param evals ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE loc_write_restart(qs_loc_env,section,mo_array, coeff_localized,&
             do_homo, evals, error)
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: section
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: coeff_localized
    LOGICAL, INTENT(IN)                      :: do_homo
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: evals
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'loc_write_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: my_middle
    INTEGER                                  :: handle, ispin, max_block, &
                                                nao, nloc, nmo, output_unit, &
                                                rst_unit
    LOGICAL                                  :: failure
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY (logger)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    IF (qs_loc_env%do_localize) THEN

      print_key => section_vals_get_subs_vals(section,"LOC_RESTART",error=error)
      IF (BTEST(cp_print_key_should_output(logger%iter_info,&
                section,"LOC_RESTART",error=error),&
                cp_p_file)) THEN

         ! Open file
         rst_unit = -1

         IF(do_homo) THEN
            my_middle="LOC_HOMO"
         ELSE
            my_middle="LOC_LUMO"
         END IF

         rst_unit = cp_print_key_unit_nr(logger,section,"LOC_RESTART",&
                    extension=".wfn", file_status="REPLACE", file_action="WRITE",&
                    file_form="UNFORMATTED",middle_name=TRIM(my_middle), error=error)

         IF(rst_unit>0) filename = cp_print_key_generate_filename(logger,print_key,&
                    middle_name=TRIM(my_middle),extension=".wfn",&
                    my_local=.FALSE.,error=error)

         IF(output_unit>0) THEN
            WRITE (UNIT=output_unit,FMT="(/,T2,A, A/)")&
                 "LOCALIZATION| Write restart file for the localized MOS : ",&
                  TRIM(filename)
         END IF

         IF(rst_unit>0) THEN
           WRITE(rst_unit) qs_loc_env%localized_wfn_control%set_of_states
           WRITE(rst_unit) qs_loc_env%localized_wfn_control%lu_bound_states
           WRITE(rst_unit) qs_loc_env%localized_wfn_control%nloc_states
         END IF

         DO ispin = 1,SIZE(coeff_localized)
           mo_coeff => coeff_localized(ispin)%matrix
           CALL cp_fm_get_info(mo_coeff, nrow_global=nao,ncol_global=nmo,ncol_block=max_block,error=error)
           nloc=qs_loc_env%localized_wfn_control%nloc_states(ispin)
           IF(rst_unit>0) THEN
             WRITE(rst_unit) qs_loc_env%localized_wfn_control%loc_states(1:nloc,ispin)
             IF(do_homo) THEN
               WRITE (rst_unit) nmo,&
                     mo_array(ispin)%mo_set%homo,&
                     mo_array(ispin)%mo_set%lfomo,&
                     mo_array(ispin)%mo_set%nelectron
               WRITE (rst_unit) mo_array(ispin)%mo_set%eigenvalues(1:nmo),&
                     mo_array(ispin)%mo_set%occupation_numbers(1:nmo)
             ELSE
               WRITE (rst_unit) nmo
               WRITE (rst_unit) evals(ispin)%array(1:nmo)
             END IF
           END IF

           CALL cp_fm_write_unformatted(mo_coeff,rst_unit,error)

         END DO

         CALL cp_print_key_finished_output(rst_unit,logger,section,&
                    "LOC_RESTART", error=error)
      END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE  loc_write_restart


! *****************************************************************************
!> \brief ...
!> \param qs_loc_env ...
!> \param mos ...
!> \param mos_localized ...
!> \param section ...
!> \param section2 ...
!> \param para_env ...
!> \param do_homo ...
!> \param restart_found ...
!> \param evals ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE  loc_read_restart(qs_loc_env,mos,mos_localized,section, section2,para_env,&
              do_homo,restart_found, evals,error)

    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_localized
    TYPE(section_vals_type), POINTER         :: section, section2
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: do_homo
    LOGICAL, INTENT(INOUT)                   :: restart_found
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: evals
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'loc_read_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=25)                        :: fname_key
    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: my_middle
    INTEGER :: group, handle, homo_read, i, ispin, istat, lfomo_read, &
      max_nloc, n_rep_val, nao, nelectron_read, nloc, nmo, nmo_read, nspin, &
      output_unit, rst_unit, source
    LOGICAL                                  :: failure, file_exists
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eig_read, occ_read
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    failure = .FALSE.

    nspin = SIZE(mos_localized)
    nao = mos(1)%mo_set%nao
    rst_unit = -1

    group = para_env%group
    source = para_env%source
    output_unit = cp_print_key_unit_nr(logger,section2,&
             "PROGRAM_RUN_INFO", extension=".Log",error=error)

    IF(do_homo) THEN
      fname_key="LOCHOMO_RESTART_FILE_NAME"
    ELSE
      fname_key="LOCLUMO_RESTART_FILE_NAME"
      CALL cp_assert(PRESENT(evals),cp_failure_level,cp_assertion_failed,routineP,&
          "Missing argument to localize unoccupied states."//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
    END IF

    file_exists = .FALSE.
    CALL section_vals_val_get(section,fname_key,n_rep_val=n_rep_val,error=error)
    IF (n_rep_val>0) THEN
      CALL section_vals_val_get(section,fname_key,c_val=filename,error=error)
    ELSE

       print_key => section_vals_get_subs_vals(section2,"LOC_RESTART",error=error)
       IF(do_homo) THEN
          my_middle="LOC_HOMO"
       ELSE
          my_middle="LOC_LUMO"
       END IF
       filename = cp_print_key_generate_filename(logger,print_key,&
                  middle_name=TRIM(my_middle),extension=".wfn",&
                  my_local=.FALSE.,error=error)
    END IF

    IF (para_env%ionode)  INQUIRE(FILE=filename,exist=file_exists)

    IF(file_exists) THEN
      IF (para_env%ionode) THEN
        CALL open_file(file_name=filename,&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_status="OLD",&
                     unit_number=rst_unit)

         READ(rst_unit) qs_loc_env%localized_wfn_control%set_of_states
         READ(rst_unit) qs_loc_env%localized_wfn_control%lu_bound_states
         READ(rst_unit) qs_loc_env%localized_wfn_control%nloc_states
      END IF
    ELSE
       IF(output_unit > 0) WRITE(output_unit,"(/,T10,A)")&
               "Restart file not available filename=<"//TRIM(filename)//'>'
    END IF
    CALL mp_bcast(file_exists,source,group)

    IF(file_exists) THEN
      restart_found = .TRUE.


    CALL mp_bcast(qs_loc_env%localized_wfn_control%set_of_states,source,group)
    CALL mp_bcast(qs_loc_env%localized_wfn_control%lu_bound_states,source,group)
    CALL mp_bcast(qs_loc_env%localized_wfn_control%nloc_states,source,group)

    max_nloc = MAXVAL(qs_loc_env%localized_wfn_control%nloc_states(:))

    ALLOCATE (vecbuffer(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    IF(ASSOCIATED(qs_loc_env%localized_wfn_control%loc_states)) THEN
      DEALLOCATE(qs_loc_env%localized_wfn_control%loc_states,STAT=istat)
    END IF
    ALLOCATE(qs_loc_env%localized_wfn_control%loc_states(max_nloc,2),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    qs_loc_env%localized_wfn_control%loc_states = 0

    DO ispin=1,nspin
       IF(do_homo) THEN
         nmo=mos(ispin)%mo_set%nmo
       ELSE
         nmo = SIZE(evals(ispin)%array,1)
       END IF
       IF (para_env%ionode.AND.(nmo > 0)) THEN
         nloc = qs_loc_env%localized_wfn_control%nloc_states(ispin)
         READ(rst_unit) qs_loc_env%localized_wfn_control%loc_states(1:nloc,ispin)
         IF(do_homo) THEN
           READ (rst_unit) nmo_read, homo_read, lfomo_read, nelectron_read
           ALLOCATE(eig_read(nmo_read), occ_read(nmo_read), STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
           eig_read = 0.0_dp
           occ_read = 0.0_dp
           READ (rst_unit) eig_read(1:nmo_read), occ_read(1:nmo_read)
         ELSE
           READ (rst_unit) nmo_read
           ALLOCATE(eig_read(nmo_read), STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
           eig_read = 0.0_dp
           READ (rst_unit) eig_read(1:nmo_read)
         END IF
         CALL cp_assert((nmo_read >= nmo),cp_warning_level,cp_assertion_failed,routineP,&
              "The number of MOs on the restart unit is smaller than the number of "//&
              "the allocated MOs. "//&
              CPSourceFileRef,&
              only_ionode=.TRUE.)
         CALL cp_assert((nmo_read<=nmo),cp_warning_level,cp_assertion_failed,routineP,&
             "The number of MOs on the restart unit is greater than the number of "//&
             "the allocated MOs. The read MO set will be truncated!"//&
             CPSourceFileRef,&
             only_ionode=.TRUE.)

         nmo = MIN(nmo,nmo_read)
         IF(do_homo) THEN
           mos(ispin)%mo_set%eigenvalues(1:nmo) = eig_read(1:nmo)
           mos(ispin)%mo_set%occupation_numbers(1:nmo) = occ_read(1:nmo)
           DEALLOCATE(eig_read, occ_read, STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ELSE
           evals(ispin)%array(1:nmo) = eig_read(1:nmo)
           DEALLOCATE(eig_read, STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         END IF

      END IF
      IF(do_homo) THEN
        CALL mp_bcast(mos(ispin)%mo_set%eigenvalues,source,group)
        CALL mp_bcast(mos(ispin)%mo_set%occupation_numbers,source,group)
      ELSE
        CALL mp_bcast(evals(ispin)%array,source,group)
      END IF

      DO i = 1,nmo
          IF (para_env%ionode) THEN
             READ (rst_unit) vecbuffer
          ELSE
             vecbuffer(1,:) = 0.0_dp
          END IF
          CALL mp_bcast(vecbuffer,source,group)
          CALL cp_fm_set_submatrix(mos_localized(ispin)%matrix,&
                 vecbuffer,1,i,nao,1,transpose=.TRUE.,error=error)
      END DO
    END DO

    CALL mp_bcast(qs_loc_env%localized_wfn_control%loc_states,source,group)

    DEALLOCATE(vecbuffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END IF

    ! Close restart file
    IF (para_env%ionode) THEN
      IF(file_exists)  CALL close_file(unit_number=rst_unit)
    END IF

    CALL timestop(handle)

  END SUBROUTINE  loc_read_restart

! *****************************************************************************
!> \brief initializes everything needed for localization of the HOMOs
!> \param qs_loc_env ...
!> \param loc_section ...
!> \param do_homo ...
!> \param do_xas ...
!> \param nloc_xas ...
!> \param error ...
!> \par History
!>      2009 created
! *****************************************************************************
  SUBROUTINE qs_loc_control_init(qs_loc_env,loc_section,do_homo,do_xas,nloc_xas,error)

    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: loc_section
    LOGICAL, INTENT(IN)                      :: do_homo
    LOGICAL, INTENT(IN), OPTIONAL            :: do_xas
    INTEGER, INTENT(IN), OPTIONAL            :: nloc_xas
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_loc_control_init', &
      routineP = moduleN//':'//routineN

    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control

    NULLIFY(localized_wfn_control)

    CALL localized_wfn_control_create(localized_wfn_control,error=error)
    CALL set_qs_loc_env(qs_loc_env,localized_wfn_control=localized_wfn_control,error=error)
    CALL localized_wfn_control_release(localized_wfn_control,error=error)
    CALL get_qs_loc_env(qs_loc_env,localized_wfn_control=localized_wfn_control,error=error)
    localized_wfn_control%do_homo=do_homo
    CALL read_loc_section(localized_wfn_control,loc_section,&
         qs_loc_env%do_localize,do_xas,nloc_xas,error=error)

  END SUBROUTINE  qs_loc_control_init

! *****************************************************************************
!> \brief initializes everything needed for localization of the molecular orbitals
!> \param qs_env ...
!> \param qs_loc_env ...
!> \param localize_section ...
!> \param mos_localized ...
!> \param do_homo ...
!> \param do_mo_cubes ...
!> \param mo_loc_history ...
!> \param evals ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_loc_init(qs_env,qs_loc_env,localize_section,mos_localized, &
             do_homo,do_mo_cubes,mo_loc_history,evals,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: localize_section
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_localized
    LOGICAL, OPTIONAL                        :: do_homo, do_mo_cubes
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mo_loc_history
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: evals
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_loc_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, i, ilast_intocc, ilow, ispin, iup, n_mo(2), &
      n_mos(2), nao, nelectron, nmoloc(2), nspin, output_unit
    LOGICAL                                  :: my_do_homo, my_do_mo_cubes, &
                                                restart_found
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues, occupation
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, mo_derivs
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: loc_print_section

    CALL timeset(routineN,handle)

    NULLIFY(mos,mo_coeff,mo_eigenvalues,occupation,ks_rmpv,mo_derivs,scf_control,para_env)
    CALL get_qs_env(qs_env,&
                    mos=mos,&
                    matrix_ks=ks_rmpv,&
                    mo_derivs=mo_derivs,&
                    scf_control=scf_control,&
                    para_env=para_env,&
                    error=error)

    loc_print_section => section_vals_get_subs_vals(localize_section,"PRINT",error=error)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    nspin = SIZE(mos)
    IF(PRESENT(do_homo)) THEN 
       my_do_homo = do_homo
    ELSE
      my_do_homo = .TRUE.
    END IF
    IF(PRESENT(do_mo_cubes)) THEN 
       my_do_mo_cubes = do_mo_cubes
    ELSE
      my_do_mo_cubes = .FALSE.
    END IF
    restart_found = .FALSE.
    IF (qs_loc_env%do_localize) THEN
      ! Some setup for MOs to be localized
      CALL get_qs_loc_env(qs_loc_env,localized_wfn_control=localized_wfn_control,error=error)
      IF(localized_wfn_control%loc_restart) THEN
         CALL loc_read_restart(qs_loc_env,mos,mos_localized,localize_section,&
              loc_print_section,para_env,my_do_homo,restart_found,evals=evals,error=error)
         IF (output_unit>0) WRITE(output_unit,"(/,T2,A,A)") "LOCALIZATION| ",&
                    "   The orbitals to be localized are read from localization restart file."
         nmoloc = localized_wfn_control%nloc_states
      END IF
      IF(.NOT. restart_found) THEN
         nmoloc = 0
         DO ispin = 1,nspin
            CALL get_mo_set(mos(ispin)%mo_set,nmo=n_mo(ispin),nelectron=nelectron,homo=homo,nao=nao,&
                 mo_coeff=mo_coeff, eigenvalues=mo_eigenvalues, occupation_numbers=occupation)
            ! Get eigenstates (only needed if not already calculated before)
            IF ((.NOT.my_do_mo_cubes)&
!                  .OR. section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO",error=error)==0)&
                  .AND.my_do_homo .AND. qs_env%scf_env%method==ot_method_nr) THEN
               CALL make_mo_eig(mos,nspin,ks_rmpv,scf_control,mo_derivs,error=error)
            END IF
            IF(localized_wfn_control%set_of_states == state_loc_all.AND.my_do_homo) THEN
               nmoloc(ispin) = NINT(nelectron/occupation(1))
               IF(n_mo(ispin)>homo) THEN
                  DO i= nmoloc(ispin),1,-1
                     IF(occupation(1)-occupation(i) < localized_wfn_control%eps_occ) THEN
                        ilast_intocc = i
                        EXIT
                     END IF
                  END DO
               ELSE
                  ilast_intocc = nmoloc(ispin)
               END IF
               nmoloc(ispin) = ilast_intocc
               localized_wfn_control%lu_bound_states(1,ispin) = 1
               localized_wfn_control%lu_bound_states(2,ispin) = ilast_intocc
               IF(nmoloc(ispin)/=n_mo(ispin)) THEN
                  IF (output_unit>0) WRITE(output_unit,"(/,T2,A,I4,A,I6,A,F12.6,A,F12.6,A)") "LOCALIZATION| Spin ", &
                       ispin ," The first ",&
                       ilast_intocc," occupied orbitals are localized, corresponding to energies from ",&
                       mo_eigenvalues(1), " to ",  mo_eigenvalues(ilast_intocc) , " [a.u.]."
               END IF
            ELSE IF(localized_wfn_control%set_of_states == energy_loc_range.AND.my_do_homo) THEN
               ilow = 0
               iup = 0
               DO i = 1,n_mo(ispin)
                  IF(mo_eigenvalues(i)>=localized_wfn_control%lu_ene_bound(1) ) THEN
                     ilow = i
                     EXIT
                  END IF
               END DO
               DO i = n_mo(ispin),1,-1
                  IF(mo_eigenvalues(i)<=localized_wfn_control%lu_ene_bound(2) ) THEN
                     iup = i
                     EXIT
                  END IF
               END DO
               localized_wfn_control%lu_bound_states(1,ispin) = ilow
               localized_wfn_control%lu_bound_states(2,ispin) = iup
               localized_wfn_control%nloc_states(ispin) = iup-ilow +1
               nmoloc(ispin) = localized_wfn_control%nloc_states(ispin)
               IF(occupation(ilow)-occupation(iup) > localized_wfn_control%eps_occ) THEN
                  CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,&
                       "The selected energy range includes orbitals with different occupation number. "//&
                       " The localization procedure cannot be applied."//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
               END IF
               IF (output_unit>0) WRITE(output_unit,"(/,T2,A,I4,A,I6,A)") "LOCALIZATION| Spin ", ispin ," : ",&
                    nmoloc(ispin)," orbitals in the selected energy range are localized."
            ELSE
               nmoloc(ispin) = MIN(localized_wfn_control%nloc_states(1) ,n_mo(ispin))
               IF (output_unit>0.AND.my_do_homo) WRITE(output_unit,"(/,T2,A,I4,A,I6,A)") "LOCALIZATION| Spin ", ispin,&
                  " : ", nmoloc(ispin)," occupied orbitals are localized, as given in the input list."
               IF (output_unit>0.AND.(.NOT.my_do_homo)) WRITE(output_unit,"(/,T2,A,I4,A,I6,A)") "LOCALIZATION| Spin ",&
                   ispin, " : ", nmoloc(ispin)," unoccupied orbitals are localized, as given in the input list."
               IF(n_mo(ispin)>homo.AND.my_do_homo) THEN
                  ilow = localized_wfn_control%loc_states(1,ispin)
                  DO i = 2,nmoloc(ispin)
                    iup = localized_wfn_control%loc_states(i,ispin)
                    IF(ABS(occupation(ilow) - occupation(iup))>localized_wfn_control%eps_occ) THEN
                  ! write warning
                      CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                       "User requested the calculation of localized wavefunction from a subset of MOs, "//&
                       "including MOs with different occupations. Check the selected subset, "//&
                       " the electronic density is not invariant with "//&
                       "respect to rotations among orbitals with different occupation numbers!"//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
                    END IF
                  END DO
               END IF
            ENDIF
         ENDDO  ! ispin
         n_mos(:)=nao-n_mo(:)
         IF(my_do_homo)n_mos=n_mo
         CALL set_loc_wfn_lists(localized_wfn_control,nmoloc,n_mos,nspin,error=error)
      END IF
      CALL set_loc_centers(localized_wfn_control,nmoloc,nspin,error=error)
      IF(my_do_homo)THEN
         CALL qs_loc_env_init(qs_loc_env,localized_wfn_control,qs_env, &
                     loc_coeff=mos_localized,mo_loc_history=mo_loc_history,error=error)
      END IF
    ELSE
       ! Let's inform in case the section is not present in the input
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
            "User requested the calculation of the localized wavefunction but the section "//&
            "LOCALIZE (located in: DFT%LOCALIZE) was not specified. Localization "//&
            "will not be performed!"//&
            CPSourceFileRef,&
            only_ionode=.TRUE.)
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_loc_init

! *****************************************************************************
!> \brief read the controlparameter from input, using the new input scheme
!> \param localized_wfn_control ...
!> \param loc_section ...
!> \param localize ...
!> \param do_xas ...
!> \param nloc_xas ...
!> \param error ...
!> \par History
!>      05.2005 created [MI]
! *****************************************************************************
  SUBROUTINE read_loc_section(localized_wfn_control,loc_section,&
       localize,do_xas,nloc_xas,error)

    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(section_vals_type), POINTER         :: loc_section
    LOGICAL, INTENT(OUT)                     :: localize
    LOGICAL, INTENT(IN), OPTIONAL            :: do_xas
    INTEGER, INTENT(IN), OPTIONAL            :: nloc_xas
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_loc_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ind, ir, n_list, n_rep, &
                                                n_state, nline, output_unit, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: list, loc_list
    LOGICAL                                  :: failure, my_do_xas
    REAL(dp), POINTER                        :: ene(:)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: loc_print_section

    failure=.FALSE.
    my_do_xas=.FALSE.
    IF(PRESENT(do_xas)) THEN
       my_do_xas=do_xas
       CPPrecondition(PRESENT(nloc_xas),cp_failure_level,routineP,error,failure)
    END IF
    CPPrecondition(ASSOCIATED(loc_section),cp_failure_level,routineP,error,failure)
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    IF(.NOT. failure)  THEN
       CALL section_vals_val_get(loc_section,"_SECTION_PARAMETERS_",l_val=localize,error=error)
       IF(localize) THEN
          loc_print_section => section_vals_get_subs_vals(loc_section,"PRINT",error=error)
          NULLIFY (list)
          NULLIFY (loc_list)
          localized_wfn_control%lu_bound_states = 0
          localized_wfn_control%lu_ene_bound = 0.0_dp
          localized_wfn_control%nloc_states     = 0
          localized_wfn_control%set_of_states   = 0
          n_state = 0

          CALL section_vals_val_get(loc_section,"MAX_ITER",&
               i_val=localized_wfn_control%max_iter,error=error)
          CALL section_vals_val_get(loc_section,"MAX_CRAZY_ANGLE",&
               r_val=localized_wfn_control%max_crazy_angle,error=error)
          CALL section_vals_val_get(loc_section,"CRAZY_SCALE",&
               r_val=localized_wfn_control%crazy_scale,error=error)
          CALL section_vals_val_get(loc_section,"EPS_OCCUPATION",&
               r_val=localized_wfn_control%eps_occ,error=error)
          CALL section_vals_val_get(loc_section,"CRAZY_USE_DIAG",&
               l_val=localized_wfn_control%crazy_use_diag,error=error)
          CALL section_vals_val_get(loc_section,"OUT_ITER_EACH",&
               i_val=localized_wfn_control%out_each,error=error)
          CALL section_vals_val_get(loc_section,"EPS_LOCALIZATION",&
               r_val=localized_wfn_control%eps_localization,error=error)
          CALL section_vals_val_get(loc_section,"MIN_OR_MAX",&
               i_val=localized_wfn_control%min_or_max,error=error)
          CALL section_vals_val_get(loc_section,"JACOBI_FALLBACK",&
               l_val=localized_wfn_control%jacobi_fallback,error=error)
          CALL section_vals_val_get(loc_section,"METHOD",&
               i_val=localized_wfn_control%localization_method,error=error)
          CALL section_vals_val_get(loc_section,"OPERATOR",&
               i_val=localized_wfn_control%operator_type,error=error)
          CALL section_vals_val_get(loc_section,"RESTART",&
               l_val=localized_wfn_control%loc_restart,error=error)
          CALL section_vals_val_get(loc_section,"USE_HISTORY",&
               l_val=localized_wfn_control%use_history,error=error)

          IF(localized_wfn_control%do_homo)THEN
            ! List of States HOMO
            CALL section_vals_val_get(loc_section,"LIST",  n_rep_val=n_rep,error=error)
            IF(n_rep>0) THEN
              n_list = 0
              DO ir = 1,n_rep
                 NULLIFY(list)
                 CALL section_vals_val_get(loc_section,"LIST",i_rep_val=ir,i_vals=list,error=error)
                 IF(ASSOCIATED(list)) THEN
                    CALL reallocate(loc_list,1,n_list+SIZE(list))
                    DO i = 1, SIZE(list)
                       loc_list(n_list+i) = list(i)
                    END DO ! i
                    n_list = n_list +  SIZE(list)
                 END IF
              END DO ! ir
              IF(n_list/=0) THEN
                 localized_wfn_control%set_of_states=state_loc_list
                 ALLOCATE(localized_wfn_control%loc_states(n_list,2),STAT=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                 localized_wfn_control%loc_states = 0
                 localized_wfn_control%loc_states(:,1) = loc_list(:)
                 localized_wfn_control%loc_states(:,2) = loc_list(:)
                 localized_wfn_control%nloc_states(1) = n_list
                 IF(my_do_xas) THEN
                   localized_wfn_control%nloc_states(2) = 0
                   localized_wfn_control%loc_states(:,2) = 0
                 END IF
                 DEALLOCATE(loc_list, STAT=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
               END IF
            END IF

          ELSE
            ! List of States LUMO
            CALL section_vals_val_get(loc_section,"LIST_UNOCCUPIED",  n_rep_val=n_rep,error=error)
            IF(n_rep>0) THEN
              n_list = 0
              DO ir = 1,n_rep
                 NULLIFY(list)
                 CALL section_vals_val_get(loc_section,"LIST_UNOCCUPIED",i_rep_val=ir,i_vals=list,error=error)
                 IF(ASSOCIATED(list)) THEN
                    CALL reallocate(loc_list,1,n_list+SIZE(list))
                    DO i = 1, SIZE(list)
                       loc_list(n_list+i) = list(i)
                    END DO ! i
                    n_list = n_list +  SIZE(list)
                 END IF
              END DO ! ir
              IF(n_list/=0) THEN
                 localized_wfn_control%set_of_states=state_loc_list
                 ALLOCATE(localized_wfn_control%loc_states(n_list,2),STAT=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                 localized_wfn_control%loc_states = 0
                 localized_wfn_control%loc_states(:,1) = loc_list(:)
                 localized_wfn_control%loc_states(:,2) = loc_list(:)
                 localized_wfn_control%nloc_states(1) = n_list
                 DEALLOCATE(loc_list, STAT=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              END IF
            END IF
          END IF

          IF(localized_wfn_control%set_of_states == 0) THEN
            CALL section_vals_val_get(loc_section,"ENERGY_RANGE",r_vals=ene,error=error)
            IF(ene(1)/=ene(2)) THEN
               localized_wfn_control%set_of_states = energy_loc_range
               localized_wfn_control%lu_ene_bound(1) = ene(1)
               localized_wfn_control%lu_ene_bound(2) = ene(2)
            END IF
          END IF

          ! All States or XAS specific states
          IF(localized_wfn_control%set_of_states == 0)  THEN
             IF(my_do_xas) THEN
                localized_wfn_control%set_of_states = state_loc_range
                localized_wfn_control%nloc_states(1) =  nloc_xas
                localized_wfn_control%lu_bound_states(1,1) = 1
                localized_wfn_control%lu_bound_states(2,1) = nloc_xas
                localized_wfn_control%nloc_states(2) =  0
                localized_wfn_control%lu_bound_states(1,2) = 0
                localized_wfn_control%lu_bound_states(2,2) = 0
             ELSE
                localized_wfn_control%set_of_states = state_loc_all
             END IF
          END IF

          localized_wfn_control%print_centers = &
               BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
               "WANNIER_CENTERS",error=error),cp_p_file)
          localized_wfn_control%print_spreads = &
               BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
               "WANNIER_SPREADS",error=error),cp_p_file)
          localized_wfn_control%print_cubes = &
               BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
               "WANNIER_CUBES",error=error),cp_p_file)

          output_unit = cp_print_key_unit_nr(logger,loc_print_section,"PROGRAM_RUN_INFO",&
               extension=".Log",error=error)

          IF (output_unit > 0) THEN
             WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
                  "LOCALIZE| The spread relative to a set of orbitals is computed"

             SELECT CASE (localized_wfn_control%set_of_states)
             CASE (state_loc_all)
                WRITE (UNIT=output_unit,FMT="(T2,A)")&
                     "LOCALIZE| Orbitals to be localized: All the occupied"
                WRITE (UNIT=output_unit,FMT="(T2,A,F16.8)")&
                     "LOCALIZE| If fractional occupation, fully occupied MOs are those within occupation tolerance of ", &
                       localized_wfn_control%eps_occ
             CASE (state_loc_range)
                WRITE (UNIT=output_unit,FMT="(T2,A,T65,I8,A,I8)")&
                     "LOCALIZE| Orbitals to be localized: Those with index between ",&
                     localized_wfn_control%lu_bound_states(1,1), " and ",&
                     localized_wfn_control%lu_bound_states(2,1)
             CASE (state_loc_list)
                WRITE (UNIT=output_unit,FMT="(T2,A)")&
                     "LOCALIZE| Orbitals to be localized: Those with index in the following list"
                nline = localized_wfn_control%nloc_states(1)/10+1
                ind = 0
                DO i = 1,nline
                  IF(ind+10<localized_wfn_control%nloc_states(1)) THEN
                    WRITE (UNIT=output_unit,FMT="(T8,10I7)") localized_wfn_control%loc_states(ind+1:ind+10,1)
                    ind=ind+10
                  ELSE
                    WRITE (UNIT=output_unit,FMT="(T8,10I7)") &
                        localized_wfn_control%loc_states(ind+1:localized_wfn_control%nloc_states(1),1)
                    ind = localized_wfn_control%nloc_states(1)
                  END IF
                END DO
             CASE (energy_loc_range)
                WRITE (UNIT=output_unit,FMT="(T2,A,T65,/,f16.6,A,f16.6,A)")&
                   "LOCALIZE| Orbitals to be localized: Those with energy in the range between ",&
                   localized_wfn_control%lu_ene_bound(1), " and ",  localized_wfn_control%lu_ene_bound(2)," a.u."
             CASE DEFAULT
                WRITE (UNIT=output_unit,FMT="(T2,A)")&
                     "LOCALIZE| Orbitals to be localized: None "
             END SELECT

             SELECT CASE (localized_wfn_control%operator_type)
             CASE(op_loc_berry)
                WRITE (UNIT=output_unit,FMT="(T2,A)")&
                     "LOCALIZE| Spread defined by the Berry phase operator "
             CASE(op_loc_boys)
                WRITE (UNIT=output_unit,FMT="(T2,A)")&
                     "LOCALIZE| Spread defined by the Boys phase operator "
             CASE DEFAULT
                WRITE (UNIT=output_unit,FMT="(T2,A)")&
                     "LOCALIZE| Spread defined by the Pipek phase operator "
             END SELECT

             SELECT CASE (localized_wfn_control%localization_method)
             CASE (do_loc_jacobi)
                WRITE (UNIT=output_unit,FMT="(T2,A)")&
                     "LOCALIZE| The optimal unitary transformation is generated by the Jacobi algorithm"
             CASE (do_loc_crazy)
                WRITE (UNIT=output_unit,FMT="(T2,A)")&
                     "LOCALIZE| The optimal unitary transformation is generated by the Crazy angle algorithm"
                WRITE (UNIT=output_unit,FMT="(T2,A,F16.8)")&
                     "LOCALIZE| maximum angle: ",localized_wfn_control%max_crazy_angle
                WRITE (UNIT=output_unit,FMT="(T2,A,F16.8)")&
                     "LOCALIZE| scaling: ",localized_wfn_control%crazy_scale
                WRITE (UNIT=output_unit,FMT="(T2,A,L1)")&
                     "LOCALIZE| use diag:",localized_wfn_control%crazy_use_diag
             CASE (do_loc_direct)
                WRITE (UNIT=output_unit,FMT="(T2,A)")&
                     "LOCALIZE| The optimal unitary transformation is generated by the direct algorithm"
             CASE (do_loc_l1_norm_sd)
                WRITE (UNIT=output_unit,FMT="(T2,A)")&
                     "LOCALIZE| The optimal unitary transformation is generated by the "//&
                     " steepest descent algorithm applyed on an approximate l1 norm"
             CASE (do_loc_none)
                WRITE (UNIT=output_unit,FMT="(T2,A)")&
                     "LOCALIZE| No unitary transformation is applied"
             END SELECT

          END IF ! process has output_unit

          CALL cp_print_key_finished_output(output_unit,logger,loc_print_section,"PROGRAM_RUN_INFO",error=error)

       ELSE
          localized_wfn_control%localization_method = do_loc_none
          localized_wfn_control%localization_method = state_loc_none
          localized_wfn_control%print_centers  = .FALSE.
          localized_wfn_control%print_spreads  = .FALSE.
          localized_wfn_control%print_cubes  = .FALSE.
       END IF
    END IF

  END SUBROUTINE read_loc_section

! *****************************************************************************
!> \brief create the center and spread array and the file names for the output
!> \param localized_wfn_control ...
!> \param nmoloc ...
!> \param nspins ...
!> \param error ...
!> \par History
!>      04.2005 created [MI]
! *****************************************************************************
  SUBROUTINE set_loc_centers(localized_wfn_control,nmoloc,nspins,error)

    TYPE(localized_wfn_control_type)         :: localized_wfn_control
    INTEGER, DIMENSION(2), INTENT(IN)        :: nmoloc
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_loc_centers', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    DO ispin = 1,nspins
       ALLOCATE(localized_wfn_control%centers_set(ispin)%array(6,nmoloc(ispin)), STAT = istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       localized_wfn_control%centers_set(ispin)%array = 0.0_dp
    END DO

  END SUBROUTINE set_loc_centers

! *****************************************************************************
!> \brief create the lists of mos that are taken into account
!> \param localized_wfn_control ...
!> \param nmoloc ...
!> \param nmo ...
!> \param nspins ...
!> \param error ...
!> \par History
!>      04.2005 created [MI]
! *****************************************************************************
  SUBROUTINE set_loc_wfn_lists(localized_wfn_control,nmoloc,nmo,nspins,error)

    TYPE(localized_wfn_control_type)         :: localized_wfn_control
    INTEGER, DIMENSION(2), INTENT(IN)        :: nmoloc, nmo
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_loc_wfn_lists', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ispin, istat, max_iloc, &
                                                max_nmoloc, state
    LOGICAL                                  :: failure

    CALL timeset(routineN,state)

    failure = .FALSE.

    localized_wfn_control%nloc_states(1:2) = nmoloc(1:2)
    max_nmoloc = MAX(nmoloc(1),nmoloc(2))

    SELECT CASE(localized_wfn_control%set_of_states )
    CASE(state_loc_list)
       ! List
       CPPrecondition(ASSOCIATED(localized_wfn_control%loc_states),cp_failure_level,routineP,error,failure)
       DO ispin = 1,nspins
          localized_wfn_control%lu_bound_states(1,ispin) = 1
          localized_wfn_control%lu_bound_states(2,ispin) = nmoloc(ispin)
          IF(nmoloc(ispin)< 1) THEN
              localized_wfn_control%lu_bound_states(1,ispin) = 0
              localized_wfn_control%loc_states(:,ispin) = 0
          END IF
       END DO
    CASE(state_loc_range)
       ! Range
       ALLOCATE(localized_wfn_control%loc_states(max_nmoloc,2),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       localized_wfn_control%loc_states = 0
       DO ispin = 1,nspins
          localized_wfn_control%lu_bound_states(1,ispin) = &
               localized_wfn_control%lu_bound_states(1,1)
          localized_wfn_control%lu_bound_states(2,ispin) = &
               localized_wfn_control%lu_bound_states(1,1) +  nmoloc(ispin) -1
          max_iloc = localized_wfn_control%lu_bound_states(2,ispin)
          DO i = 1,nmoloc(ispin)
             localized_wfn_control%loc_states(i,ispin) = localized_wfn_control%lu_bound_states(1,ispin) + i -1
          END DO
          CPPostcondition(max_iloc<=nmo(ispin),cp_failure_level,routineP,error,failure)
       END DO
    CASE(energy_loc_range)
       ! Energy
       ALLOCATE(localized_wfn_control%loc_states(max_nmoloc,2),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       localized_wfn_control%loc_states = 0
       DO ispin = 1,nspins
          DO i = 1,nmoloc(ispin)
             localized_wfn_control%loc_states(i,ispin) = localized_wfn_control%lu_bound_states(1,ispin) + i -1
          END DO
       END DO
    CASE(state_loc_all)
       ! All
       ALLOCATE(localized_wfn_control%loc_states(max_nmoloc,2),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       localized_wfn_control%loc_states = 0

       DO ispin = 1,nspins
          localized_wfn_control%lu_bound_states(1,ispin) = 1
          localized_wfn_control%lu_bound_states(2,ispin) = nmoloc(ispin)
          IF(nmoloc(ispin)< 1) localized_wfn_control%lu_bound_states(1,ispin) = 0
          DO i = 1, nmoloc(ispin)
             localized_wfn_control%loc_states(i,ispin) = i
          END DO
       END DO
    END SELECT

    CALL timestop(state)

  END SUBROUTINE set_loc_wfn_lists


END MODULE qs_loc_utils

