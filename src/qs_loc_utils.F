!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Some utilities for the construction of
!>      the localization environment
!> \author MI (05-2005)
! *****************************************************************************
MODULE qs_loc_utils

  USE ai_moments,                      ONLY: contract_cossin,&
                                             cossin
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_1d_r_p_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_set,&
                                             dbcsr_obj_p_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_element, cp_fm_get_info, cp_fm_get_submatrix, &
       cp_fm_maxabsval, cp_fm_p_type, cp_fm_release, cp_fm_set_all, &
       cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_methods,                   ONLY: dbcsr_init
  USE dbcsr_operations,                ONLY: dbcsr_get_block_p
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE input_constants,                 ONLY: do_loc_none,&
                                             op_loc_berry,&
                                             op_loc_pipek,&
                                             state_loc_list
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi,&
                                             twopi
  USE message_passing,                 ONLY: mp_bcast
  USE orbital_pointers,                ONLY: nco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_loc_control,                  ONLY: localized_wfn_control_type
  USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                             qs_loc_env_new_type,&
                                             set_qs_loc_env
  USE qs_localization_methods,         ONLY: initialize_weights,&
                                             rotate_orbitals
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_loc_utils'

! *** Public ***
  PUBLIC ::  jacobi_rotation_pipek, qs_loc_env_init, loc_write_restart, loc_read_restart

CONTAINS

! *****************************************************************************
!> \brief 2by2 rotation for the pipek operator
!>       in this case the z_ij numbers are reals
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

    INTEGER                                  :: iatom, natom, nstate, stride
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9)                    :: desc
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
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      similar to the old one, but not quite
!> \par History
!>      04.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE qs_loc_env_init(qs_loc_env,localized_wfn_control,qs_env,myspin,do_localize,loc_coeff,error)

    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN), OPTIONAL            :: myspin
    LOGICAL, INTENT(IN), OPTIONAL            :: do_localize
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: loc_coeff
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
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: moloc_coeff
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mat_ptr, mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
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
            para_env=para_env,mos=mos,error=error)
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
                           CALL stop_program(routineP,"States with different occupations"//&
                           " cannot be rotated together")
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
                     CALL stop_program(routineP,"States with different occupations"//&
                     " cannot be rotated together")
             END IF
             nmosub = localized_wfn_control%nloc_states(ispin)
             CALL cp_fm_to_fm(mat_ptr,moloc_coeff(ispin)%matrix,nmosub,lb,1)
          END IF
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
                CALL dbcsr_init(qs_loc_env%op_sm_set(j,i)%matrix)
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
          CALL stop_program(routineP,"Type of operator not implemented")
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
!> \note
!>      The intgrals are computed analytically  using the primitives GTO
!>      The contraction is performed block-wise
!> \par History
!>      04.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE get_berry_operator ( qs_loc_env, qs_env, error )
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_berry_operator', &
      routineP = moduleN//':'//routineN

    INTEGER :: dim_op, handle, i, iab, iatom, icol, ikind, ilist, inode, &
      irow, iset, istat, ithread, j, jatom, jkind, jset, last_jatom, ldab, &
      ldsa, ldsb, ldwork, maxl, ncoa, ncob, nkind, nlist, nnode, nrow, nseta, &
      nsetb, nthread, sgfa, sgfb
    INTEGER, DIMENSION(3)                    :: perd0
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found, new_atom_b, tr
    REAL(KIND=dp)                            :: dab, kvec(3), rab2, &
                                                vector_k(3,6)
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cosab, rpgfa, rpgfb, sinab, &
                                                sphi_a, sphi_b, work, zeta, &
                                                zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: cosabt, sinabt, workt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: op_cos, op_sin
    TYPE(block_p_type), DIMENSION(:, :), &
      POINTER                                :: op_cost, op_sint
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dbcsr_obj_p_type), &
      DIMENSION(:, :), POINTER               :: op_sm_set
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY (atomic_kind, atomic_kind_set)
    NULLIFY ( cell, op_sm_set , orb_basis_set, particle_set)
    NULLIFY (sab_orb,sab_orb_neighbor_list,sab_orb_neighbor_list_local,&
             sab_orb_neighbor_node)
    NULLIFY (cosabt,sinabt, workt,cosab,sinab, work, distribution_2d)
    NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set, sab_orb=sab_orb,&
                    distribution_2d=distribution_2d,error=error)

    nkind = SIZE(atomic_kind_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=ldwork, maxlgto=maxl )
    ldab = ldwork
    nthread = 1
!$  nthread = omp_get_max_threads()
    ALLOCATE(cosabt(ldab,ldab,0:nthread-1))
    cosabt = 0.0_dp
    ALLOCATE(sinabt(ldab,ldab,0:nthread-1))
    sinabt = 0.0_dp
    ALLOCATE(workt(ldwork,ldwork,0:nthread-1))
    workt = 0.0_dp

    CALL get_qs_loc_env(qs_loc_env=qs_loc_env, op_sm_set=op_sm_set, &
         cell=cell, dim_op=dim_op, error=error)

    ALLOCATE (op_cost(dim_op,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (op_sint(dim_op,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, dim_op
      DO j=0,nthread-1
        NULLIFY (op_cost(i,j)%block)
        NULLIFY (op_sint(i,j)%block)
      END DO
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

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             pgf_radius=rpgfa,&
                             set_radius=set_radius_a,&
                             sphi=sphi_a,&
                             zet=zeta)
      ldsa = SIZE(sphi_a,1)
      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               first_sgf=first_sgfb,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               nsgf_set=nsgfb,&
                               pgf_radius=rpgfb,&
                               set_radius=set_radius_b,&
                               sphi=sphi_b,&
                               zet=zetb)

        iab = ikind + nkind*(jkind - 1)
        ldsb = SIZE(sphi_b,1)

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        CALL get_neighbor_list_set(neighbor_list_set=&
                                   sab_orb(iab)%neighbor_list_set,&
                                   nlist=nlist)

        NULLIFY ( sab_orb_neighbor_list )

!$OMP parallel &
!$OMP default(none) &
!$OMP private(ithread,cosab,sinab,work,ilist,iatom,jatom)&
!$OMP private(last_jatom,new_atom_b,istat,sgfa,sgfb,ra,rb)&
!$OMP private(nrow,iset,jset,ncoa,ncob,nnode,inode,irow,icol,found,tr)&
!$OMP private(op_cos,op_sin,rab,rab2,dab,i,kvec)&
!$OMP private(sab_orb_neighbor_node,sab_orb_neighbor_list_local)&
!$OMP shared(ikind,jkind,rpgfb,zetb,sphi_b,sphi_a,nsgfa,nsgfb,set_radius_b,ldab,ldwork)&
!$OMP shared(lb_max,npgfb,nsetb,iab,la_min,rpgfa,zeta)&
!$OMP shared(set_radius_a,nco,ncoset,npgfa,nseta,ldsa,ldsb)&
!$OMP shared(nkind,sab_orb_neighbor_list,nlist,workt,cosabt,sinabt,qs_env)&
!$OMP shared(op_cost,op_sint,dim_op,vector_k)&
!$OMP shared(first_sgfa,first_sgfb,sab_orb,lb_min,la_max,op_sm_set)&
!$OMP shared(atomic_kind_set,particle_set,cell)
        ithread = 0
!$      ithread = omp_get_thread_num()

        cosab => cosabt(:,:,ithread)
        sinab => sinabt(:,:,ithread)
        work  => workt(:,:,ithread)
        op_cos => op_cost(:,ithread)
        op_sin => op_sint(:,ithread)

!$OMP do
       DO ilist = 1,nlist

!$OMP critical (qs_core_neighbor_list)
          IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
            sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
          ELSE
            sab_orb_neighbor_list => next(sab_orb_neighbor_list)
          END IF
          sab_orb_neighbor_list_local => sab_orb_neighbor_list
!$OMP end critical (qs_core_neighbor_list)

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                 atom=iatom,&
                                 nnode=nnode)
          ra = pbc(particle_set(iatom)%r,cell)
          last_jatom = 0

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom, r=rab)

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
                 CALL dbcsr_get_block_p(matrix=op_sm_set(1,i)%matrix,&
                      row=irow,col=icol,block=op_cos(i)%block,tr=tr,found=found)
                 NULLIFY(op_sin(i)%block)
                 CALL dbcsr_get_block_p(matrix=op_sm_set(2,i)%matrix,&
                      row=irow,col=icol,block=op_sin(i)%block,tr=tr,found=found)
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

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO ! jatom => atom B

       END DO ! iatom => atom A

!$OMP end parallel

      END DO ! jkind

    END DO ! ikind

    ! Set back the correct periodicity
    cell%perd(1:3) =  perd0(1:3)

    DO i = 1,dim_op
      DO j = 0,nthread-1
         NULLIFY(op_cost(i,j)%block)
         NULLIFY(op_sint(i,j)%block)
      END DO
    END DO
    DEALLOCATE(op_cost,op_sint, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(cosabt,sinabt,workt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)
  END SUBROUTINE get_berry_operator

  SUBROUTINE loc_write_restart(qs_env,qs_loc_env,section,mo_array, coeff_localized, do_homo, evals, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
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
    INTEGER                                  :: handle, i, i_block, ispin, &
                                                istat, j, max_block, nao, &
                                                nloc, nmo, output_unit, &
                                                rst_unit
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer
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
                  filename
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
   
           ALLOCATE(vecbuffer(nao,max_block),STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 
           DO i=1,nmo,MAX(max_block,1)
               i_block=MIN(max_block,nmo-i+1)
               CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,i,nao,i_block,error=error)
               ! doing this in one write would increase efficiency, but breaks RESTART compatibility.
               ! to old ones, and in cases where max_block is different between runs, as might happen during
               ! restarts with a different number of CPUs
               DO j=1,i_block
                  IF (rst_unit>0) WRITE (rst_unit) vecbuffer(1:nao,j)
               ENDDO
   
           END DO
   
           DEALLOCATE (vecbuffer,STAT=istat)
           CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 
         END DO
   
   
         CALL cp_print_key_finished_output(rst_unit,logger,section,&
                    "LOC_RESTART", error=error)
      END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE  loc_write_restart


  SUBROUTINE  loc_read_restart(qs_loc_env,mos,mos_localized,section, section2,para_env,do_homo,evals,error)

    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_localized
    TYPE(section_vals_type), POINTER         :: section, section2
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: do_homo
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
      rst_unit, source
    LOGICAL                                  :: exist, failure
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

    IF(do_homo) THEN
      fname_key="LOCHOMO_RESTART_FILE_NAME"
    ELSE
      fname_key="LOCLUMO_RESTART_FILE_NAME"
    END IF

    exist = .FALSE. 
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

    INQUIRE(FILE=filename,exist=exist)

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

    ! Close restart file
    IF (para_env%ionode) CALL close_file(unit_number=rst_unit)

    CALL timestop(handle)

  END SUBROUTINE  loc_read_restart


END MODULE qs_loc_utils

