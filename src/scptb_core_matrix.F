!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the Hamiltonian integral and Overlap matrices
!> <a|b> and <a|T+V|b> for SCP-TB methods
!> \author JGH
! *****************************************************************************
MODULE scptb_core_matrix
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_copy,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_init
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_allocate_matrix_set
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_output_handling,              ONLY: cp_print_key_should_output
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: use_orb_basis_set
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_core_hamiltonian,             ONLY: dump_info_core_hamiltonian
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kinetic,                      ONLY: build_kinetic_matrix
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: build_overlap_matrix
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scptb_types,                     ONLY: get_scptb_parameter,&
                                             scptb_parameter_p_type,&
                                             scptb_parameter_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scptb_core_matrix'

  PUBLIC :: build_scptb_core_matrix

CONTAINS

! *****************************************************************************
  SUBROUTINE build_scptb_core_matrix(qs_env,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_scptb_core_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, i, iatom, icol, icor, ikind, inode, &
      irow, istat, j, jatom, jkind, l, li, ll, maxder, natom, nder, nkind
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: defined, failure, found, &
                                                use_virial
    LOGICAL, DIMENSION(:), POINTER           :: scptb_defined
    REAL(KIND=dp)                            :: dr, f0, kh, kij
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: hmt, kk, umt
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ds_block, h_block, p_block, &
                                                s_block
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p, matrix_s, &
                                                matrix_t, matrix_w
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(scptb_parameter_p_type), &
      DIMENSION(:), POINTER                  :: scptb_kind_param
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind_a
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    NULLIFY ( rho,force,atomic_kind_set,sab_orb,matrix_p )

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    matrix_h=matrix_h,&
                    kinetic=matrix_t,&
                    matrix_s=matrix_s,&
                    particle_set=particle_set,&
                    sab_orb=sab_orb,&
                    virial=virial,&
                    para_env=para_env,error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    IF (calculate_forces) THEN
       nder = 1

       ALLOCATE (atom_of_kind(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

       CALL get_qs_env(qs_env=qs_env,force=force,matrix_w=matrix_w,error=error)

       CALL get_qs_env(qs_env=qs_env, rho=rho,error=error)
       matrix_p => rho%rho_ao
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                         alpha_scalar= 1.0_dp, beta_scalar=1.0_dp,error=error)
          CALL cp_dbcsr_add(matrix_p(2)%matrix, matrix_p(1)%matrix, &
                         alpha_scalar=-2.0_dp, beta_scalar=1.0_dp,error=error)
          CALL cp_dbcsr_add(matrix_w(1)%matrix, matrix_w(2)%matrix, &
                         alpha_scalar= 1.0_dp, beta_scalar=1.0_dp,error=error)
       END IF
    ELSE
       IF (cp_print_key_should_output(logger%iter_info,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/DERIVATIVES",error=error)/=0) THEN
          nder = 1
       ELSE
          nder = 0
       END IF
    END IF

    maxder = ncoset(nder)

    IF (calculate_forces) THEN
       ! S matrix
       CALL build_overlap_matrix(qs_env,nderivative=nder,matrix_s=matrix_s,&
                                 matrix_name="OVERLAP MATRIX",&
                                 basis_set_id_a=use_orb_basis_set,&
                                 basis_set_id_b=use_orb_basis_set, &
                                 sab_nl=sab_orb,calculate_forces=.TRUE.,&
                                 matrix_p=matrix_w(1)%matrix,&
                                 error=error)
       ! T matrix
       CALL build_kinetic_matrix(qs_env,nderivative=nder,matrix_t=matrix_t,&
                                 matrix_name="KINETIC ENERGY MATRIX",&
                                 basis_set_id=use_orb_basis_set,&
                                 sab_nl=sab_orb,calculate_forces=.TRUE.,&
                                 matrix_p=matrix_p(1)%matrix,&
                                 error=error)
    ELSE
       ! S matrix
       CALL build_overlap_matrix(qs_env,nderivative=nder,matrix_s=matrix_s,&
                                 matrix_name="OVERLAP MATRIX",&
                                 basis_set_id_a=use_orb_basis_set,&
                                 basis_set_id_b=use_orb_basis_set, &
                                 sab_nl=sab_orb,&
                                 error=error)
       ! T matrix
       CALL build_kinetic_matrix(qs_env,nderivative=nder,matrix_t=matrix_t,&
                                 matrix_name="KINETIC ENERGY MATRIX",&
                                 basis_set_id=use_orb_basis_set,&
                                 sab_nl=sab_orb,&
                                 error=error)
    END IF

    ! initialize H matrix
    CALL cp_dbcsr_allocate_matrix_set(matrix_h,1,error)
    ALLOCATE(matrix_h(1)%matrix,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(matrix_h(1)%matrix, error=error)
    CALL cp_dbcsr_copy(matrix_h(1)%matrix,matrix_t(1)%matrix,name="CORE HAMILTONIAN MATRIX",error=error)
    !The band energy elements will be added to the kinetic energy.

    ALLOCATE (scptb_kind_param(nkind),scptb_defined(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (umt(100,nkind),hmt(100,nkind),kk(nkind,nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    umt = 0._dp
    hmt = 0._dp
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,scptb_parameter=scptb_kind_a)
       scptb_kind_param(ikind)%scptb_param => scptb_kind_a
       CALL get_scptb_parameter(scptb_kind_a,defined=defined)
       scptb_defined(ikind)=defined
       li = 0
       DO l=0,scptb_kind_a%lmaxorb
          ll = ((l+1)*(l+2))/2
          DO i=1,scptb_kind_a%norb(l)
            hmt(li+1:li+ll,ikind) = scptb_kind_a%hcore(i,l)
            li = li + ll
          END DO
       END DO
    END DO
    !
    kk = 0.5_dp

    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,inode=inode,r=rij)
       IF (.NOT.scptb_defined(ikind)) CYCLE
       IF (.NOT.scptb_defined(jkind)) CYCLE

       kh = kk(ikind,jkind)

       NULLIFY (h_block)
       IF (iatom <= jatom) THEN
         irow = iatom
         icol = jatom
       ELSE
         irow = jatom
         icol = iatom
       END IF
       CALL cp_dbcsr_get_block_p(matrix_h(1)%matrix,irow,icol,h_block,found)
       CPPostcondition(ASSOCIATED(h_block),cp_failure_level,routineP,error,failure)

       dr = SUM(rij(:)**2)
       !actual core energy determination , eq. 2.17 in my notes ???
       IF (iatom == jatom .AND. dr < 0.0001_dp) THEN
          DO i=1,SIZE(h_block,1)
             h_block(i,i)=h_block(i,i)+hmt(i,ikind)
          END DO
       ELSE
          NULLIFY(s_block) !s_block is overlap [sub]matrix ???
          CALL cp_dbcsr_get_block_p(matrix_s(1)%matrix,irow,icol,s_block,found)
          CPPostcondition(ASSOCIATED(s_block),cp_failure_level,routineP,error,failure)
          IF ( irow == iatom ) THEN
            DO i=1,SIZE(h_block,1)
              DO j=1,SIZE(h_block,2)
                kij = kh*(hmt(i,ikind)+hmt(j,jkind))
                h_block(i,j)=h_block(i,j)+kij*s_block(i,j)
              END DO
            END DO
          ELSE
            DO i=1,SIZE(h_block,1)
              DO j=1,SIZE(h_block,2)
                kij = kh*(hmt(j,ikind)+hmt(i,jkind))
                h_block(i,j)=h_block(i,j)+kij*s_block(i,j)
              END DO
            END DO
          END IF
          IF(calculate_forces) THEN
            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)
            CALL cp_dbcsr_get_block_p(matrix_p(1)%matrix,irow,icol,p_block,found)
            CPPostcondition(ASSOCIATED(p_block),cp_failure_level,routineP,error,failure)
            DO icor=1,3
              force_ab(icor) = 0._dp
              CALL cp_dbcsr_get_block_p(matrix_s(icor+1)%matrix,irow,icol,ds_block,found)
              CPPostcondition(ASSOCIATED(ds_block),cp_failure_level,routineP,error,failure)
              IF ( irow == iatom ) THEN
                 DO i=1,SIZE(h_block,1)
                   DO j=1,SIZE(h_block,2)
                     kij = 2._dp*kh*(hmt(i,ikind)+hmt(j,jkind))
                     force_ab(icor)=force_ab(icor)+kij*ds_block(i,j)*p_block(i,j)
                   END DO
                 END DO
               ELSE
                 DO i=1,SIZE(h_block,1)
                   DO j=1,SIZE(h_block,2)
                     kij = 2._dp*kh*(hmt(j,ikind)+hmt(i,jkind))
                     force_ab(icor)=force_ab(icor)+kij*ds_block(i,j)*p_block(i,j)
                   END DO
                 END DO
               END IF
            END DO
            IF ( irow == iatom ) force_ab = -force_ab
            force(ikind)%all_potential(:,atom_a) = force(ikind)%all_potential(:,atom_a) - force_ab(:)
            force(jkind)%all_potential(:,atom_b) = force(jkind)%all_potential(:,atom_b) + force_ab(:)
            IF (use_virial) THEN
               f0 = -1.0_dp
               IF(iatom==jatom) f0=-0.5_dp
               CALL virial_pair_force ( virial%pv_virial, f0, force_ab, rij, error)
            END IF
          END IF

       END IF

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    DEALLOCATE (scptb_kind_param,scptb_defined,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (umt,hmt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF (calculate_forces) THEN
       DEALLOCATE (atom_of_kind,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       IF (SIZE(matrix_p) == 2) THEN
          CALL cp_dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                         alpha_scalar= 0.5_dp, beta_scalar=0.5_dp,error=error)
          CALL cp_dbcsr_add(matrix_p(2)%matrix, matrix_p(1)%matrix, &
                         alpha_scalar=-1.0_dp, beta_scalar=1.0_dp,error=error)
       END IF
    END IF

    ! *** Put the core Hamiltonian matrix in the QS environment ***
    CALL set_qs_env(qs_env=qs_env,&
                    matrix_s=matrix_s,&
                    kinetic=matrix_t,&
                    matrix_h=matrix_h,&
                    error=error)

    ! Print matrices if requested
    CALL dump_info_core_hamiltonian(matrix_s, matrix_t, matrix_h, &
         qs_env, calculate_forces, para_env, error)
    CALL timestop(handle)

  END SUBROUTINE build_scptb_core_matrix

END MODULE scptb_core_matrix

