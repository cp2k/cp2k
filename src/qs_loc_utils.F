!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_loc_utils [1.0] *
!!
!!   NAME
!!     qs_loc_utils
!!
!!   FUNCTION
!!     Some utilities for the construction of
!!     the localization environment
!!
!!   AUTHORS
!!     MI (05-2005)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_loc_utils

! *****************************************************************************

  USE ai_moments,                      ONLY: contract_cossin,&
                                             cossin
  USE atomic_kind_types,               ONLY: atomic_kind_type, &
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_maxabsval,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp
  USE input_constants,                 ONLY: do_loc_min,&
                                             do_loc_max,&
                                             do_loc_jacobi,&
                                             do_loc_crazy,&
                                             do_loc_direct,&
                                             do_loc_none,&
                                             op_loc_berry,&
                                             op_loc_boys,&
                                             op_loc_pipek,&
                                             state_loc_all,&
                                             state_loc_upton,&
                                             state_loc_fromitoj,&
                                             state_loc_list
  USE mathconstants,                   ONLY: pi, twopi
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
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             block_p_type,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE cell_types,                 ONLY: cell_type,&
                                             pbc
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             get_matrix_info
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: module_name = "qs_loc_utils"


! *** Public ***
  PUBLIC ::  jacobi_rotation_pipek, qs_loc_env_init

! *****************************************************************************

CONTAINS


!******************************************************************************
!!****f* qs_loc_types/jacpbi_rotation_pipek [1.0] *
!!
!!   NAME
!!     jacobi_rotation_pipek
!!
!!   FUNCTION
!!      2by2 rotation for the pipek operator
!!      in this case the z_ij numbers are reals 
!!
!!   ARGUMENTS
!!
!!  
!!   AUTHOR
!!     MI
!!
!!   HISTORY
!!      05-2005 created
!!
!!*** **********************************************************************
  SUBROUTINE jacobi_rotation_pipek(zij_fm_set, vectors, sweeps, error)

    TYPE(cp_fm_p_type), DIMENSION(:,:),POINTER   :: zij_fm_set
    TYPE(cp_fm_type), POINTER                    :: vectors
    INTEGER                                      :: sweeps
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                   :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'jacobi_rotation_pipek', &
      routineP = module_name//':'//routine_name

    INTEGER :: iatom, istate, jstate, natom, nstate
     
    LOGICAL                                      :: failure
    REAL(KIND=dp)                                :: ct, st, theta, tolerance
    REAL(KIND=dp)                                :: mii, mij, mjj
    REAL(KIND=dp)                                :: aij, bij,ratio
    TYPE(cp_fm_type), POINTER                    :: rmat

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
    CALL cp_fm_release ( rmat )

  END SUBROUTINE jacobi_rotation_pipek

  SUBROUTINE rotate_rmat_real(istate,jstate,st,ct,rmat,error)

    INTEGER, INTENT(IN)              :: istate, jstate
    REAL(dp), INTENT(IN)             :: ct, st
    TYPE(cp_fm_type), POINTER        :: rmat
    INTEGER                          :: nstate
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                            :: error


#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desc
#endif

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

  SUBROUTINE rotate_zij_real(istate,jstate,st,ct,zij,error)

    INTEGER, INTENT(IN)                   :: istate, jstate
    REAL(dp), INTENT(IN)                  :: ct, st
    TYPE(cp_fm_p_type), DIMENSION(:,:), &
      POINTER                             :: zij
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                            :: error

    INTEGER :: iatom, natom, nstate, stride
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desc
#endif

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
       stride = size(zij(iatom,1)%matrix%local_data,1)
       CALL drot(nstate,zij(iatom,1)%matrix%local_data(istate,1),&
            stride,zij(iatom,1)%matrix%local_data(jstate,1),stride,ct,st)
#endif
    END DO
    
  END SUBROUTINE rotate_zij_real

  SUBROUTINE check_tolerance_real(zij_fm_set,tolerance,error)

    TYPE(cp_fm_p_type), DIMENSION(:,:)    :: zij_fm_set
    REAL(dp), INTENT(OUT)                  :: tolerance
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                            :: error
    
    INTEGER :: iatom, istat, istate, jstate, natom, ncol_local,&
               nrow_global, nrow_local 
    INTEGER, DIMENSION(:), POINTER        :: col_indices, row_indices
    LOGICAL                               :: failure
    REAL(dp)                              :: grad_ij, zii, zij, zjj
    REAL(dp), DIMENSION(:,:), POINTER     :: diag
    TYPE(cp_fm_type), POINTER             :: force

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

    Do iatom = 1,natom
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

    CALL cp_fm_release ( force ) 

  END SUBROUTINE check_tolerance_real
!******************************************************************************
!!****f* qs_loc_types/qs_loc_env_init [1.0] *
!!
!!   NAME
!!     qs_loc_env_init
!!
!!   FUNCTION
!!      allocates the data, and initializes the operators
!!
!!   NOTES
!!     similar to the old one, but not quite
!!
!!   ARGUMENTS
!!     - qs_loc_env: new environment for the localization calculations
!!     - localized_wfn_control: variables and directives for the localization 
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!*** **********************************************************************
  SUBROUTINE qs_loc_env_init(qs_loc_env,localized_wfn_control,qs_env,onlyspin1,error)

    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE ( localized_wfn_control_type ), &
      POINTER                                  :: localized_wfn_control
    TYPE(qs_environment_type), POINTER         :: qs_env
    LOGICAL, OPTIONAL                          :: onlyspin1
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'qs_loc_env_init', &
      routineP = module_name//':'//routine_name

    INTEGER  :: handle, i, iatom, imo, imoloc, isos, ispin, istat, j, &
                lb, my_nspins, nao, naosub, natoms,  nmo, nmosub, nspins, ub
    LOGICAL  :: failure, my_one_spin
    REAL(KIND=dp)                               :: my_occ, occ_imo
    REAL(KIND=dp), DIMENSION(:, :), POINTER     :: vecbuffer
    REAL(KIND=dp), DIMENSION(:), POINTER        :: occupations
    TYPE(cell_type), POINTER                    :: cell
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER   :: moloc_coeff
    TYPE(cp_fm_struct_type), POINTER            :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                   :: mo_coeff
    TYPE(cp_para_env_type), POINTER             :: para_env
    TYPE(distribution_1d_type), POINTER         :: local_molecules
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                   :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                   :: matrix_s

!   ---------------------------------------------------------------------------

    CALL timeset("qs_loc_env_init","I",'',handle)

    failure = .FALSE.
    my_one_spin = .FALSE.
    IF(PRESENT(onlyspin1)) my_one_spin = onlyspin1

    NULLIFY(mos,matrix_s,moloc_coeff,particle_set,para_env,cell,local_molecules,occupations)
  
    CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, cell=cell, &
                     local_molecules=local_molecules,&
                     particle_set=particle_set, &
                     para_env=para_env,&
                     mos=mos)

    nspins = SIZE(mos,1)
    my_nspins = nspins
    IF(my_one_spin) my_nspins = 1
    ALLOCATE (moloc_coeff(my_nspins), STAT=istat)
    DO ispin = 1,my_nspins
      NULLIFY(tmp_fm_struct,mo_coeff)
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,nao=nao,nmo=nmo)
      nmosub = localized_wfn_control%nloc_states(ispin)
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
           ncol_global=nmosub,para_env=para_env,context=mo_coeff%matrix_struct%context)

      CALL cp_fm_create (moloc_coeff(ispin)%matrix , tmp_fm_struct )

      CALL cp_fm_get_info(moloc_coeff(ispin)%matrix, nrow_global=naosub, &
           ncol_global=nmosub,error=error)
      CPPostcondition(nao==naosub,cp_failure_level,routineP,error,failure)
      CPPostcondition(nmo>=nmosub,cp_failure_level,routineP,error,failure)
      CALL cp_fm_set_all(moloc_coeff(ispin)%matrix,0.0_dp)
      CALL cp_fm_struct_release ( tmp_fm_struct )
    END DO  ! ispin
    ! Copy the submatrix

    DO ispin = 1,my_nspins
      CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,&
           occupation_numbers=occupations, nao=nao, nmo=nmo)
      lb = localized_wfn_control%lu_bound_states(1,ispin)
      ub = localized_wfn_control%lu_bound_states(2,ispin)

      IF(localized_wfn_control%set_of_states == state_loc_list) THEN
        ALLOCATE(vecbuffer(1,nao),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        nmosub = SIZE(localized_wfn_control%loc_states,1)
        CPPostcondition(nmosub>0,cp_failure_level,routineP,error,failure)
        imoloc = 0
        my_occ =  occupations(localized_wfn_control%loc_states(1))
        DO i = lb,ub
           ! Get the index in the subset
           imoloc = imoloc + 1
           ! Get the index in the full set
           imo = localized_wfn_control%loc_states(i)
           occ_imo = occupations(imo)
           IF(ABS(occ_imo-my_occ) > EPSILON(1.0_dp)) THEN
             IF(localized_wfn_control%localization_method /= do_loc_none) &
               CALL stop_program(routineP,"States with different occupations"//&
               " cannot be rotated together")
           END IF
           ! Take the imo vector from the full set and copy in the imoloc vector of the subset
           CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,imo,&
                nao,1,transpose=.TRUE.,error=error)
           CALL cp_fm_set_submatrix(moloc_coeff(ispin)%matrix,vecbuffer,1,imoloc,&
                nao,1,transpose=.TRUE.,error=error)
         END DO
         DEALLOCATE(vecbuffer,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ELSE
         my_occ = occupations( lb )
         occ_imo = occupations( ub )
         IF(ABS(occ_imo-my_occ) >  EPSILON(1.0_dp)) THEN
            IF(localized_wfn_control%localization_method /= do_loc_none) &
             CALL stop_program(routineP,"States with different occupations"//&
             " cannot be rotated together")
         END IF
         nmosub = localized_wfn_control%nloc_states(ispin)
         CALL cp_fm_to_fm(mo_coeff,moloc_coeff(ispin)%matrix,nmosub,lb,1)
       END IF
    END DO

    CALL set_qs_loc_env(qs_loc_env=qs_loc_env, cell=cell, local_molecules=local_molecules,&
                        moloc_coeff=moloc_coeff, particle_set=particle_set,para_env=para_env,&
                        localized_wfn_control=localized_wfn_control)

    ! Prepare the operators
    NULLIFY(tmp_fm_struct,mo_coeff)
    nmosub = MAXVAL(localized_wfn_control%nloc_states)
    CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff)
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmosub,&
         ncol_global=nmosub,para_env=para_env,context=mo_coeff%matrix_struct%context)

    IF ( localized_wfn_control%operator_type==op_loc_berry )  THEN
      ALLOCATE ( qs_loc_env % op_sm_set ( 2, 3 ), STAT = istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 
      DO i=1, SIZE ( qs_loc_env % op_sm_set, 2 )
         DO j=1, SIZE ( qs_loc_env % op_sm_set, 1 )
           NULLIFY(qs_loc_env%op_sm_set(j,i)%matrix)
           CALL replicate_matrix_structure(matrix_s(1)%matrix, &
                qs_loc_env%op_sm_set(j,i)%matrix,"qs_loc_env%op_sm_"//&
                TRIM(ADJUSTL(cp_to_string(j)))//"-"//TRIM(ADJUSTL(cp_to_string(i))))
           CALL set_matrix(qs_loc_env%op_sm_set(j,i)%matrix,0.0_dp)
         ENDDO
      END DO

      ALLOCATE ( qs_loc_env % op_fm_set( 2, 3 ), STAT = istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)


      DO i =   1, SIZE ( qs_loc_env % op_fm_set, 2 )
         DO j =   1, SIZE ( qs_loc_env % op_fm_set, 1 )
           NULLIFY(qs_loc_env % op_fm_set(j,i)%matrix)

           CALL cp_fm_create (qs_loc_env % op_fm_set(j,i)%matrix , tmp_fm_struct )

           CALL cp_fm_get_info(qs_loc_env %op_fm_set(j,i)%matrix, nrow_global=nmosub, &
               error=error)
           CPPostcondition(nmo>=nmosub,cp_failure_level,routineP,error,failure)
           CALL cp_fm_set_all(qs_loc_env %op_fm_set(j,i)%matrix,0.0_dp)
         END DO
      END DO
    ELSEIF ( localized_wfn_control%operator_type== op_loc_pipek) THEN
      natoms = SIZE ( qs_loc_env % particle_set, 1 )
      ALLOCATE ( qs_loc_env % op_fm_set ( natoms , 1 ), STAT = istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 
      DO ispin = 1, SIZE(qs_loc_env % op_fm_set ,2)
         CALL get_mo_set(mos(ispin)%mo_set,nmo=nmo)
         DO iatom = 1,natoms
           NULLIFY(qs_loc_env %op_fm_set(iatom,ispin)%matrix)

           CALL cp_fm_create (qs_loc_env % op_fm_set(iatom,ispin)%matrix , tmp_fm_struct )

           CALL cp_fm_get_info(qs_loc_env %op_fm_set(iatom,ispin)%matrix, nrow_global=nmosub, &
                error=error)
           CPPostcondition(nmo>=nmosub,cp_failure_level,routineP,error,failure)
           CALL cp_fm_set_all(qs_loc_env %op_fm_set(iatom,ispin)%matrix,0.0_dp)
         END DO  ! iatom
      END DO  ! ispin
    ELSE
      CALL stop_program(routineP,"Type of operator not implemented")
    ENDIF
    CALL cp_fm_struct_release ( tmp_fm_struct )

    IF (  localized_wfn_control%operator_type==op_loc_berry ) THEN 
 
      CALL initialize_weights ( qs_loc_env % cell, qs_loc_env % weights ) 
  
      CALL get_berry_operator ( qs_loc_env, qs_env , error=error)

    ELSEIF ( localized_wfn_control%operator_type== op_loc_pipek) THEN

!!    here we don't have to do anything
!!    CALL get_pipek_mezey_operator ( qs_loc_env, qs_env )

    ENDIF

    
    qs_loc_env%molecular_states = .FALSE.
  

! done
    qs_loc_env%is_initialized = 666
    CALL timestop(0.0_dp,handle) 

END SUBROUTINE qs_loc_env_init

!******************************************************************************
!!****f* qs_loc_types/get_berry_operator [1.0] *
!!
!!   NAME
!!     get_berry_operator
!!
!!   FUNCTION
!!      Computes the Berry operator for periodic systems
!!      used to define the spread of the MOS
!!      Here the matrix elements of the type <mu|cos(kr)|nu> and  <mu|sin(kr)|nu> 
!!      are computed, where mu and nu are the contracted basis functions.
!!      Namely the Berry operator is exp(ikr)
!!      k is defined somewhere
!!      the pair lists are exploited and sparse matrixes are constructed
!!
!!   NOTES
!!     The intgrals are computed analytically  using the primitives GTO
!!     The contraction is performed block-wise
!!
!!   ARGUMENTS
!!     - qs_loc_env: new environment for the localization calculations
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!*** **********************************************************************

  SUBROUTINE get_berry_operator ( qs_loc_env, qs_env, error ) 

! globals
    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE(qs_environment_type), POINTER         :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'get_berry_operator', &
      routineP = module_name//':'//routine_name


    INTEGER  :: handle, i, iab, iat, iatom, icol, ikind, ilist, inode, irow, &
                iset, istat, ithread, j, jat, jatom, jkind, jset, last_jatom, &
                ldab, ldsa, ldsb, ldwork, maxl, ncoa, ncob, nkind, nlist,& 
                nnode, nrow, nthread, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER             :: la_max, la_min, lb_max, lb_min
    INTEGER, DIMENSION(:), POINTER             :: npgfa, npgfb
    INTEGER, DIMENSION(:), POINTER             :: nsgfa, nsgfb
    INTEGER, DIMENSION(:,:), POINTER           :: first_sgfa,first_sgfb
    LOGICAL :: failure, new_atom_b, ortho
    REAL(KIND=dp)                              :: ooLx, ooLy, ooLz, rab2, dab
    REAL(KIND=dp)                              :: twopoLx, twopoLy, twopoLz
    REAL(KIND=dp), DIMENSION(3)                :: kvec, ra, rb, rab
    REAL(KIND=dp), DIMENSION(6)                :: weights
    REAL(KIND = dp), DIMENSION(:), POINTER     :: set_radius_a,set_radius_b
    REAL(KIND=dp), DIMENSION(:,:), POINTER     :: rpgfa, rpgfb, sphi_a,sphi_b, zeta, zetb 
    REAL(KIND=dp), DIMENSION(:,:), POINTER     :: cosab,sinab, work 
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER   :: cosabt,sinabt, workt
    TYPE(atomic_kind_type), DIMENSION(:),&
      POINTER                                  :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER            :: atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
          POINTER                              :: op_cos, op_sin
    TYPE(block_p_type), DIMENSION(:,:), &
      POINTER                                  :: op_cost, op_sint
    TYPE(cell_type), POINTER                   :: cell
    TYPE(gto_basis_set_type), POINTER          :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:),           POINTER          :: sab_orb
    TYPE(neighbor_list_type), POINTER          :: sab_orb_neighbor_list,&
                                                  sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER          :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                  :: particle_set
    TYPE (real_matrix_p_type), DIMENSION(:,:),&
      POINTER                                  :: op_sm_set

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num

    CALL timeset(routine_name,"I"," ",handle)

    failure = .FALSE.
    NULLIFY (atomic_kind, atomic_kind_set)
    NULLIFY ( cell, op_sm_set , orb_basis_set, particle_set)
    NULLIFY (sab_orb,sab_orb_neighbor_list,sab_orb_neighbor_list_local,&
             sab_orb_neighbor_node)
    NULLIFY (cosabt,sinabt, workt,cosab,sinab, work)
    NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )


    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set, sab_orb=sab_orb)

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

    ALLOCATE (op_cost(3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (op_sint(3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, 3
      DO j=0,nthread-1
        NULLIFY (op_cost(i,j)%block)
        NULLIFY (op_sint(i,j)%block)
      END DO
    END DO 
 
    CALL get_qs_loc_env(qs_loc_env=qs_loc_env,op_sm_set=op_sm_set, cell=cell)
    ooLx = 1.0_dp/SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) ) 
    ooLy = 1.0_dp/SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )
    ooLz = 1.0_dp/SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )
    twopoLx = twopi*ooLx
    twopoLy = twopi*ooLy
    twopoLz = twopi*ooLz
    ortho = cell%orthorhombic

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
!$OMP private(nrow,iset,jset,ncoa,ncob,nnode,inode,irow,icol)&
!$OMP private(op_cos,op_sin,rab,rab2,dab,i,kvec)&
!$OMP private(sab_orb_neighbor_node,sab_orb_neighbor_list_local)&
!$OMP shared(ikind,jkind,rpgfb,zetb,sphi_b,sphi_a,nsgfa,nsgfb,set_radius_b,ldab,ldwork)&
!$OMP shared(lb_max,npgfb,nsetb,iab,la_min,rpgfa,zeta)&
!$OMP shared(set_radius_a,nco,ncoset,npgfa,nseta,ldsa,ldsb)&
!$OMP shared(nkind,sab_orb_neighbor_list,nlist,workt,cosabt,sinabt,qs_env)&
!$OMP shared(op_cost,op_sint,ortho,twopoLx,twopoLy,twopoLz)&
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

!            rb = ra - rab
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

              DO i = 1,3
                 NULLIFY(op_cos(i)%block)
                 CALL get_block_node(matrix=op_sm_set(1,i)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=op_cos(i)%block)
                 NULLIFY(op_sin(i)%block)
                 CALL get_block_node(matrix=op_sm_set(2,i)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=op_sin(i)%block)
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
                  IF(ortho)THEN

                    ! x component
                    kvec = 0.0_dp
                    kvec(1) = twopoLx
                    cosab = 0.0_dp  
                    sinab = 0.0_dp  

                    CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                              lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                              ra,rb,kvec,cosab,sinab)

                    CALL contract_cossin(op_cos(1)%block,op_sin(1)%block,&
                               iatom,ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                               jatom,ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                               cosab,sinab,ldab,work,ldwork)

                    ! y component
                    kvec = 0.0_dp
                    kvec(2) = twopoLy
                    cosab = 0.0_dp  
                    sinab = 0.0_dp  


                    CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                              lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                              ra,rb,kvec,cosab,sinab)

                    CALL contract_cossin(op_cos(2)%block,op_sin(2)%block,&
                               iatom,ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                               jatom,ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                               cosab,sinab,ldab,work,ldwork)

                    ! z component
                    kvec = 0.0_dp
                    kvec(3) = twopoLz
                    cosab = 0.0_dp  
                    sinab = 0.0_dp  


                    CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                              lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                              ra,rb,kvec,cosab,sinab)

                    CALL contract_cossin(op_cos(3)%block,op_sin(3)%block,&
                               iatom,ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                               jatom,ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                               cosab,sinab,ldab,work,ldwork)

                  ELSE
                    CALL stop_program(routineP,&
                      "spread calculation with berry phase in non orthorombic cell NYI")
                  END IF  ! ortho

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

    DO i = 1,3
      DO j = 0,nthread-1
         NULLIFY(op_cost(i,j)%block)
         NULLIFY(op_sint(i,j)%block)
      END DO
    END DO 
    DEALLOCATE(op_cost,op_sint, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(cosabt,sinabt,workt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE get_berry_operator

!------------------------------------------------------------------------------

END MODULE qs_loc_utils

