!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004-2006 CP2K developers group
!-----------------------------------------------------------------------------!


!!****** cp2k/nddo_methods [1.1] *
!!
!!   FUNCTION
!!     Calculation of the Fock matrix for NDDO methods
!!
!!   AUTHOR
!!     JGH
!!
!!***
!******************************************************************************

MODULE nddo_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE semi_empirical_integrals,        ONLY: drotint,&
                                             drotnuc,&
                                             rotint,&
                                             rotnuc,&
                                             set_taper_fn
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE semi_empirical_utils,            ONLY: contract,&
                                             expand,&
                                             get_se_param,&
                                             get_se_type,&
                                             traceb
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             cp_sm_sm_trace,&
                                             create_dbm_from_sm,&
                                             deallocate_dbm_from_sm,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nddo_methods'

  PUBLIC :: build_nddo_fock_matrix

  INTERFACE fock1
    MODULE PROCEDURE rfock1, ufock1
  END INTERFACE

  INTERFACE fock2c
    MODULE PROCEDURE rfock2c, ufock2c
  END INTERFACE

  INTERFACE dfock2c
    MODULE PROCEDURE rdfock2c, udfock2c
  END INTERFACE

  INTERFACE fock2e
    MODULE PROCEDURE rfock2e, ufock2e
  END INTERFACE

  INTERFACE dfock2e
    MODULE PROCEDURE rdfock2e, udfock2e
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_nddo_fock_matrix (ks_env,qs_env,ks_matrix,rho,energy,&
               calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_nddo_fock_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, ncol_global, &
                                                nspins, output_unit
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers, &
                                                scaling_factor
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h
    TYPE(section_vals_type), POINTER         :: scf_section

!   ----------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    NULLIFY(matrix_h, dft_control, logger, scf_section)

    logger => cp_error_get_logger(error)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN

       CALL get_qs_env(qs_env=qs_env,&
            dft_control=dft_control,&
            matrix_h=matrix_h,&
            para_env=para_env)
 
       scf_section => section_vals_get_subs_vals(qs_env%input,"DFT%SCF",error=error)
       nspins=dft_control%nspins
       CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)
 
       DO ispin=1,nspins
         ! copy the core matrix into the fock matrix
         CALL transfer_matrix(matrix_h(1)%matrix,ks_matrix(ispin)%matrix)
       END DO
 
       CALL build_nddo_fock_matrix_exchange (ks_env,qs_env,ks_matrix,rho,energy,&
                calculate_forces,just_energy,error)
       CALL build_nddo_fock_matrix_coulomb (ks_env,qs_env,ks_matrix,rho,energy,&
               calculate_forces,just_energy,error)

       energy%mulliken = 0._dp
       energy%exc = 0._dp

       energy%total = energy%core + energy%core_overlap + 0.5_dp*energy%hartree + &
            energy%qmmm_el + energy%mulliken

       output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%DETAILED_ENERGY",&
            extension=".scfLog",error=error)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
               "Core Hamiltonian energy:                       ",energy%core,&
               "Two-electron integral energy:                  ",energy%hartree
          IF (qs_env%qmmm) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "QM/MM Electrostatic energy:                    ",energy%qmmm_el
          END IF
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%DETAILED_ENERGY", error=error)       
       ! here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C (plus occupation numbers)
       IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
          CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array)
          DO ispin=1,SIZE(mo_derivs)
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff, &
                             occupation_numbers=occupation_numbers )
             CALL cp_fm_get_info(mo_coeff,ncol_global=ncol_global)
             CALL cp_sm_fm_multiply(ks_matrix(ispin)%matrix,mo_coeff,&
                                          mo_derivs(ispin)%matrix, ncol_global)
             ALLOCATE(scaling_factor(SIZE(occupation_numbers)))
             scaling_factor=2.0_dp*occupation_numbers
             CALL cp_fm_column_scale(mo_derivs(ispin)%matrix,scaling_factor)
             DEALLOCATE(scaling_factor)
          ENDDO
       ENDIF

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_nddo_fock_matrix

! *****************************************************************************

  SUBROUTINE build_nddo_fock_matrix_coulomb (ks_env,qs_env,ks_matrix,rho,energy,&
               calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'build_nddo_fock_matrix_coulomb', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iab, iatom, icor, ikind, ilist, inode, &
      ispin, istat, itype, jatom, jkind, natom, natorb_a, natorb_b, nkind, &
      nlist, nnode, nspins
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: defined, failure, switch
    REAL(KIND=dp)                            :: delta, dr1, dr2, ecore2, &
                                                ecoul, enuc, enuclear, range, &
                                                rc
    REAL(KIND=dp), DIMENSION(10)             :: e1b, e2a, pvec
    REAL(KIND=dp), DIMENSION(10, 3)          :: de1b, de2a
    REAL(KIND=dp), DIMENSION(3)              :: denuc, force_ab, rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: ksa_block_a, ksa_block_b, &
      ksb_block_a, ksb_block_b, pa_block_a, pa_block_b, pb_block_a, pb_block_b
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_2c
    TYPE(neighbor_list_type), POINTER        :: sab_2c_neighbor_list, &
                                                sab_2c_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_2c_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: diagmat_ks, diagmat_p, &
                                                matrix_h, matrix_p
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b

!   ----------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(matrix_h,matrix_p,dft_control,cell,force,particle_set)
    NULLIFY(diagmat_ks,diagmat_p)

    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         matrix_h=matrix_h,&
         cell=cell,&
         para_env=para_env)

    ! set values for tapering function 
    rc = dft_control%qs_control%se_control%rc_coulomb
    range = dft_control%qs_control%se_control%rc_range
    CALL set_taper_fn (2._dp*rc,range)

    nspins=dft_control%nspins
    CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN

       CALL allocate_matrix_set(diagmat_ks,nspins)
       CALL allocate_matrix_set(diagmat_p,nspins)

       matrix_p => rho%rho_ao

       CALL get_qs_env(qs_env=qs_env,&
                       sab_2c=sab_2c,&
                       atomic_kind_set=atomic_kind_set)

       nkind = SIZE(atomic_kind_set)

       IF(calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env,&
                         particle_set=particle_set,&
                         force=force)
         natom = SIZE (particle_set)
         ALLOCATE (atom_of_kind(natom),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         delta = dft_control%qs_control%se_control%delta
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                  atom_of_kind=atom_of_kind)
       END IF

       DO ispin=1,nspins
          ! Allocate diagonal block matrices
          CALL create_dbm_from_sm(diagmat_p(ispin)%matrix,matrix_p(ispin)%matrix,&
                             "DIAGONAL BLOCK P MATRIX",.TRUE.,para_env)
          CALL create_dbm_from_sm(diagmat_ks(ispin)%matrix,ks_matrix(ispin)%matrix,&
                             "DIAGONAL BLOCK KS MATRIX",.FALSE.)
       END DO

       ecore2 = 0._dp
       enuclear = 0._dp
       itype = get_se_type(dft_control%qs_control%method)

       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,&
                            defined=defined,&
                            natorb=natorb_a)

          IF (.NOT.defined .OR. natorb_a < 1) CYCLE

          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                  se_parameter=se_kind_b)
             CALL get_se_param(se_kind_b,&
                               defined=defined,&
                               natorb=natorb_b)

             IF (.NOT.defined .OR. natorb_b < 1) CYCLE

             iab = ikind + nkind*(jkind - 1)
             IF (.NOT.ASSOCIATED(sab_2c(iab)%neighbor_list_set)) CYCLE
             CALL get_neighbor_list_set(neighbor_list_set=&
                                        sab_2c(iab)%neighbor_list_set,&
                                        nlist=nlist)
             NULLIFY(sab_2c_neighbor_list)

             DO ilist=1,nlist
               IF ( .NOT. ASSOCIATED(sab_2c_neighbor_list) ) THEN
                 sab_2c_neighbor_list => first_list(sab_2c(iab)%neighbor_list_set)
               ELSE
                 sab_2c_neighbor_list => next(sab_2c_neighbor_list)
               END IF
               sab_2c_neighbor_list_local => sab_2c_neighbor_list
               CALL get_neighbor_list(neighbor_list=sab_2c_neighbor_list_local,&
                                      atom=iatom,&
                                      nnode=nnode)

               sab_2c_neighbor_node => first_node(sab_2c_neighbor_list_local)

               CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
               IF ( nspins >= 1 ) THEN
                 CALL get_block_node(matrix=diagmat_p(1)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     BLOCK=pa_block_a)
                 CPPostcondition(ASSOCIATED(pa_block_a),cp_failure_level,routineP,error,failure)
                 CALL get_block_node(matrix=diagmat_ks(1)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     BLOCK=ksa_block_a)
                 CPPostcondition(ASSOCIATED(ksa_block_a),cp_failure_level,routineP,error,failure)
               END IF
               IF ( nspins >= 2 ) THEN
                 CALL get_block_node(matrix=diagmat_p(2)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     BLOCK=pa_block_b)
                 CPPostcondition(ASSOCIATED(pa_block_b),cp_failure_level,routineP,error,failure)
                 CALL get_block_node(matrix=diagmat_ks(2)%matrix,&
                                     block_row=iatom,&
                                     block_col=iatom,&
                                     BLOCK=ksa_block_b)
                 CPPostcondition(ASSOCIATED(ksa_block_b),cp_failure_level,routineP,error,failure)
               END IF

               DO inode=1,nnode
                 CALL get_neighbor_node(neighbor_node=sab_2c_neighbor_node,&
                                        neighbor=jatom,r=rij)

                 dr1 = SUM(rij(:)**2)
                 IF ( dr1 > 0._dp ) THEN

                   IF (iatom <= jatom) THEN
                     switch = .FALSE.
                   ELSE
                     switch = .TRUE.
                   END IF
                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                   IF ( nspins >= 1 ) THEN
                     CALL get_block_node(matrix=diagmat_p(1)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         BLOCK=pb_block_a)
                     CPPostcondition(ASSOCIATED(pb_block_a),cp_failure_level,routineP,error,failure)
                     CALL get_block_node(matrix=diagmat_ks(1)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         BLOCK=ksb_block_a)
                     CPPostcondition(ASSOCIATED(ksb_block_a),cp_failure_level,routineP,error,failure)
                   END IF
                   IF ( nspins >= 2 ) THEN
                     CALL get_block_node(matrix=diagmat_p(2)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         BLOCK=pb_block_b)
                     CPPostcondition(ASSOCIATED(pb_block_b),cp_failure_level,routineP,error,failure)
                     CALL get_block_node(matrix=diagmat_ks(2)%matrix,&
                                         block_row=jatom,&
                                         block_col=jatom,&
                                         BLOCK=ksb_block_b)
                     CPPostcondition(ASSOCIATED(ksb_block_b),cp_failure_level,routineP,error,failure)
                   END IF

                   SELECT CASE (dft_control%qs_control%method)
                     CASE ("EHT")
                     CASE ("MNDO","AM1","PM3","PDG")
                     IF ( .NOT. switch ) THEN
                       CALL rotnuc (se_kind_a,se_kind_b,rij,e1b,e2a,enuc,itype)
                     ELSE
                       CALL rotnuc (se_kind_b,se_kind_a,-rij,e2a,e1b,enuc,itype)
                     END IF
                     enuclear = enuclear + enuc
                     ! one-centre one-electron terms
                     CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                     IF ( nspins >= 1 ) THEN
                       ! contribution to iatom
                       IF ( SIZE(ksa_block_a,1) == 1 ) THEN
                         ksa_block_a(1,1)=ksa_block_a(1,1)+e1b(1)
                         ecore2=ecore2+e1b(1)*pa_block_a(1,1)
                       ELSE
                         CALL expand ( ksa_block_a, e1b )
                         ecore2=ecore2+traceb(pa_block_a,e1b)
                       END IF
                       ! contribution to jatom
                       IF ( SIZE(ksb_block_a,1) == 1 ) THEN
                         ksb_block_a(1,1)=ksb_block_a(1,1)+e2a(1)
                         ecore2=ecore2+e2a(1)*pb_block_a(1,1)
                       ELSE
                         CALL expand ( ksb_block_a, e2a )
                       ecore2=ecore2+traceb(pb_block_a,e2a)
                       END IF
                     END IF
                     IF ( nspins >= 2 ) THEN
                       ! contribution to iatom
                       IF ( SIZE(ksa_block_b,1) == 1 ) THEN
                         ksa_block_b(1,1)=ksa_block_b(1,1)+e1b(1)
                         ecore2=ecore2+e1b(1)*pa_block_b(1,1)
                       ELSE
                         CALL expand ( ksa_block_b, e1b )
                         ecore2=ecore2+traceb(pa_block_b,e1b)
                       END IF
                       ! contribution to jatom
                       IF ( SIZE(ksb_block_b,1) == 1 ) THEN
                         ksb_block_b(1,1)=ksb_block_b(1,1)+e2a(1)
                         ecore2=ecore2+e2a(1)*pb_block_b(1,1)
                       ELSE
                         CALL expand ( ksb_block_b, e2a )
                         ecore2=ecore2+traceb(pb_block_b,e2a)
                       END IF
                     END IF
                     IF(calculate_forces) THEN
                       atom_a = atom_of_kind(iatom)
                       atom_b = atom_of_kind(jatom)
                       IF ( .NOT. switch ) THEN
                         CALL drotnuc (se_kind_a,se_kind_b,rij,de1b,de2a,&
                                       denuc,itype,delta)
                       ELSE
                         CALL drotnuc (se_kind_b,se_kind_a,-rij,de2a,de1b,&
                                       denuc,itype,delta)
                       END IF
                       force_ab(1:3)=-denuc(1:3)
                       CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                       IF ( nspins >= 1 ) THEN
                           IF ( SIZE(ksa_block_a,1) == 1 ) THEN
                             force_ab(1:3)=force_ab(1:3)-de1b(1,1:3)*pa_block_a(1,1)
                           ELSE
                             CALL contract ( pa_block_a, pvec )
                             DO icor=1,3
                               force_ab(icor)=force_ab(icor)-&
                                              SUM(de1b(1:10,icor)*pvec(1:10))
                             END DO
                           END IF
                           IF ( SIZE(ksb_block_a,1) == 1 ) THEN
                             force_ab(1:3)=force_ab(1:3)-de2a(1,1:3)*pb_block_a(1,1)
                           ELSE
                             CALL contract ( pb_block_a, pvec )
                             DO icor=1,3
                               force_ab(icor)=force_ab(icor)-&
                                              SUM(de2a(1:10,icor)*pvec(1:10))
                             END DO
                           END IF
                       END IF
                       IF ( nspins >= 2 ) THEN
                           IF ( SIZE(ksa_block_b,1) == 1 ) THEN
                             force_ab(1:3)=force_ab(1:3)-de1b(1,1:3)*pa_block_b(1,1)
                           ELSE
                             CALL contract ( pa_block_b, pvec )
                             DO icor=1,3
                               force_ab(icor)=force_ab(icor)-&
                                              SUM(de1b(1:10,icor)*pvec(1:10))
                             END DO
                           END IF
                           IF ( SIZE(ksb_block_b,1) == 1 ) THEN
                             force_ab(1:3)=force_ab(1:3)-de2a(1,1:3)*pb_block_b(1,1)
                           ELSE
                             CALL contract ( pb_block_b, pvec )
                             DO icor=1,3
                               force_ab(icor)=force_ab(icor)-&
                                              SUM(de2a(1:10,icor)*pvec(1:10))
                             END DO
                           END IF
                       END IF
                       IF ( .NOT. switch ) force_ab = -force_ab
                          force(ikind)%all_potential(:,atom_a) =&
                            force(ikind)%all_potential(:,atom_a) - force_ab(:)
                          force(jkind)%all_potential(:,atom_b) =&
                            force(jkind)%all_potential(:,atom_b) + force_ab(:)
                     END IF

                     CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                     IF ( nspins == 1 ) THEN
                       IF ( .NOT. switch ) THEN
                         CALL fock2c(se_kind_a,se_kind_b,rij,pa_block_a,pb_block_a,&
                                    ksa_block_a,ksb_block_a,error)
                       ELSE
                         CALL fock2c(se_kind_b,se_kind_a,-rij,pb_block_a,pa_block_a,&
                                    ksb_block_a,ksa_block_a,error)
                       ENDIF
                     ELSE IF ( nspins == 2 ) THEN
                       IF ( .NOT. switch ) THEN
                         CALL fock2c(se_kind_a,se_kind_b,rij,pa_block_a,pb_block_a,&
                                    ksa_block_a,ksb_block_a,&
                                    pa_block_b,pb_block_b,&
                                    ksa_block_b,ksb_block_b,error)
                       ELSE
                         CALL fock2c(se_kind_b,se_kind_a,-rij,pb_block_a,pa_block_a,&
                                    ksb_block_a,ksa_block_a,&
                                    pb_block_b,pa_block_b,&
                                    ksb_block_b,ksa_block_b,error)
                       ENDIF
                     END IF

                     IF(calculate_forces) THEN
                       force_ab = 0._dp
                       CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                       IF ( nspins == 1 ) THEN
                         IF ( .NOT. switch ) THEN
                           CALL dfock2c(se_kind_a,se_kind_b,rij,pa_block_a,&
                                pb_block_a,force_ab,delta,error)
                         ELSE
                           CALL dfock2c(se_kind_b,se_kind_a,-rij,pb_block_a,&
                                pa_block_a,force_ab,delta,error)
                         ENDIF
                       ELSE IF ( nspins == 2 ) THEN
                         IF ( .NOT. switch ) THEN
                           CALL dfock2c(se_kind_a,se_kind_b,rij,pa_block_a,&
                                pb_block_a,pa_block_b,pb_block_b,force_ab,delta,error)
                         ELSE
                           CALL dfock2c(se_kind_b,se_kind_a,-rij,pb_block_a,&
                                pa_block_a,pb_block_b,pa_block_b,force_ab,delta,error)
                         ENDIF
                       END IF
                       IF ( switch ) force_ab = -force_ab
                       force(ikind)%rho_elec(:,atom_a) =&
                         force(ikind)%rho_elec(:,atom_a) - force_ab(:)
                       force(jkind)%rho_elec(:,atom_b) =&
                         force(jkind)%rho_elec(:,atom_b) + force_ab(:)
                     END IF

                   END SELECT

                 END IF

                 sab_2c_neighbor_node => next(sab_2c_neighbor_node)

               END DO ! inode => jatom(atom B)

             END DO ! ilist => iatom(atom A)

          END DO ! jkind

       END DO ! ikind

       DO ispin = 1, nspins
          CALL deallocate_dbm_from_sm(diagmat_p(ispin)%matrix,&
                                      matrix_p(1)%matrix,.FALSE.)
          CALL deallocate_dbm_from_sm(diagmat_ks(ispin)%matrix,&
                                      ks_matrix(ispin)%matrix,.TRUE.,para_env)
       END DO
       DEALLOCATE(diagmat_ks,diagmat_p)

       IF (calculate_forces) DEALLOCATE(atom_of_kind)

       CALL mp_sum(enuclear,para_env%group)
       CALL mp_sum(ecore2,para_env%group)

       energy%core_overlap = enuclear
       energy%hartree = ecore2 - energy%core
       DO ispin=1,nspins
          CALL cp_sm_sm_trace(ks_matrix(ispin)%matrix,&
                              matrix_p(ispin)%matrix,trace=ecoul,&
                              para_env=para_env,error=error)
          energy%hartree=energy%hartree+ecoul
       END DO

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_nddo_fock_matrix_coulomb

! *****************************************************************************

  SUBROUTINE build_nddo_fock_matrix_exchange (ks_env,qs_env,ks_matrix,rho,energy,&
               calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'build_nddo_fock_matrix_exchange', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iab, iatom, icol, ikind, ilist, inode, &
      irow, istat, jatom, jkind, natom, natorb_a, natorb_b, nkind, nlist, &
      nnode, nspins
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: defined, failure
    REAL(KIND=dp)                            :: delta, dr, gp2, gpp, gsp, &
                                                gss, hsp, range, rc
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ks_block_a, ks_block_b, &
                                                p_block_a, p_block_b
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b

!   ----------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(matrix_h,matrix_p,dft_control,cell,force,particle_set)

    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         matrix_h=matrix_h,&
         cell=cell,&
         para_env=para_env)

    ! set values for tapering function
    rc = dft_control%qs_control%se_control%rc_interaction
    range = dft_control%qs_control%se_control%rc_range
    CALL set_taper_fn (2._dp*rc,range)

    nspins=dft_control%nspins
    CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN

       matrix_p => rho%rho_ao

       CALL get_qs_env(qs_env=qs_env,&
                       sab_orb=sab_orb,&
                       atomic_kind_set=atomic_kind_set)

       nkind = SIZE(atomic_kind_set)

       IF(calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env,&
                         particle_set=particle_set,&
                         force=force)
         natom = SIZE (particle_set)
         ALLOCATE (atom_of_kind(natom),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         delta = dft_control%qs_control%se_control%delta
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                  atom_of_kind=atom_of_kind)
       END IF

       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,&
                            defined=defined,&
                            natorb=natorb_a,&
                            gss=gss,gsp=gsp,gpp=gpp,gp2=gp2,hsp=hsp)

          IF (.NOT.defined .OR. natorb_a < 1) CYCLE

          DO jkind=1,nkind
             atomic_kind => atomic_kind_set(jkind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                  se_parameter=se_kind_b)
             CALL get_se_param(se_kind_b,&
                               defined=defined,&
                               natorb=natorb_b)

             IF (.NOT.defined .OR. natorb_b < 1) CYCLE

             iab = ikind + nkind*(jkind - 1)
             IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE
             CALL get_neighbor_list_set(neighbor_list_set=&
                                        sab_orb(iab)%neighbor_list_set,&
                                        nlist=nlist)
             NULLIFY(sab_orb_neighbor_list)

             DO ilist=1,nlist
               IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
                 sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
               ELSE
                 sab_orb_neighbor_list => next(sab_orb_neighbor_list)
               END IF
               sab_orb_neighbor_list_local => sab_orb_neighbor_list
               CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                      atom=iatom,&
                                      nnode=nnode)

               sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

               DO inode=1,nnode
                 CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                        neighbor=jatom,r=rij)

                 IF (iatom <= jatom) THEN
                   irow = iatom
                   icol = jatom
                 ELSE
                   irow = jatom
                   icol = iatom
                 END IF
                 CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                 IF ( nspins >= 1 ) THEN
                   CALL get_block_node(matrix=ks_matrix(1)%matrix,&
                                       block_row=irow,&
                                       block_col=icol,&
                                       BLOCK=ks_block_a)
                   CPPostcondition(ASSOCIATED(ks_block_a),cp_failure_level,routineP,error,failure)
                   CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                       block_row=irow,&
                                       block_col=icol,&
                                       BLOCK=p_block_a)
                   CPPostcondition(ASSOCIATED(p_block_a),cp_failure_level,routineP,error,failure)
                 END IF
                 IF ( nspins >= 2 ) THEN
                   CALL get_block_node(matrix=ks_matrix(2)%matrix,&
                                       block_row=irow,&
                                       block_col=icol,&
                                       BLOCK=ks_block_b)
                   CPPostcondition(ASSOCIATED(ks_block_b),cp_failure_level,routineP,error,failure)
                   CALL get_block_node(matrix=matrix_p(2)%matrix,&
                                       block_row=irow,&
                                       block_col=icol,&
                                       BLOCK=p_block_b)
                   CPPostcondition(ASSOCIATED(p_block_b),cp_failure_level,routineP,error,failure)
                 END IF

                 dr = SUM(rij(:)**2)
                 IF ( iatom == jatom .AND. dr < 0.001_dp ) THEN

                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                   IF ( nspins == 1 ) THEN
                     CALL fock1(gss,gsp,gpp,gp2,hsp,p_block_a,ks_block_a,error)
                   ELSE IF ( nspins == 2 ) THEN
                     CALL fock1(gss,gsp,gpp,gp2,hsp,p_block_a,ks_block_a,&
                                p_block_b,ks_block_b,error)
                   END IF

                 ELSE

                   CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                   IF ( nspins == 1 ) THEN
                     IF ( irow == iatom ) THEN
                       CALL fock2e(se_kind_a,se_kind_b,rij,p_block_a,ks_block_a,error)
                     ELSE
                       CALL fock2e(se_kind_b,se_kind_a,-rij,p_block_a,ks_block_a,error)
                     ENDIF
                   ELSE IF ( nspins == 2 ) THEN
                     IF ( irow == iatom ) THEN
                       CALL fock2e(se_kind_a,se_kind_b,rij,p_block_a,ks_block_a,&
                                   p_block_b,ks_block_b,error)
                     ELSE
                       CALL fock2e(se_kind_b,se_kind_a,-rij,p_block_a,ks_block_a,&
                                   p_block_b,ks_block_b,error)
                     ENDIF
                   END IF

                   IF(calculate_forces) THEN
                     force_ab(1:3) = 0._dp
                     CPPrecondition(nspins<3,cp_failure_level,routineP,error,failure)
                     IF ( nspins == 1 ) THEN
                       IF ( irow == iatom ) THEN
                         CALL dfock2e(se_kind_a,se_kind_b,rij,p_block_a,force_ab,delta,error)
                       ELSE
                         CALL dfock2e(se_kind_b,se_kind_a,-rij,p_block_a,force_ab,delta,error)
                       ENDIF
                     ELSE IF ( nspins == 2 ) THEN
                       IF ( irow == iatom ) THEN
                         CALL dfock2e(se_kind_a,se_kind_b,rij,&
                              p_block_a,p_block_b,force_ab,delta,error)
                       ELSE
                         CALL dfock2e(se_kind_b,se_kind_a,-rij,&
                              p_block_a,p_block_b,force_ab,delta,error)
                       ENDIF
                     END IF
                     atom_a = atom_of_kind(iatom)
                     atom_b = atom_of_kind(jatom)
                     IF ( irow /= iatom ) force_ab = -force_ab
                     force(ikind)%rho_elec(:,atom_a) =&
                       force(ikind)%rho_elec(:,atom_a) - force_ab(:)
                     force(jkind)%rho_elec(:,atom_b) =&
                       force(jkind)%rho_elec(:,atom_b) + force_ab(:)
                   END IF

                 END IF

                 sab_orb_neighbor_node => next(sab_orb_neighbor_node)

               END DO ! inode => jatom(atom B)

             END DO ! ilist => iatom(atom A)

          END DO ! jkind

       END DO ! ikind

       IF (calculate_forces) DEALLOCATE(atom_of_kind)

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_nddo_fock_matrix_exchange

! *****************************************************************************

  SUBROUTINE rfock1 (gss,gsp,gpp,gp2,hsp,pmat,fmat,error)
    REAL(dp), INTENT(IN)                     :: gss, gsp, gpp, gp2, hsp
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pmat
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fmat
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rfock1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n
    LOGICAL                                  :: failure

    n = SIZE ( pmat,1 )

    CPPrecondition(n<=4,cp_failure_level,routineP,error,failure)
    IF ( n == 1 ) THEN
      fmat(1,1) = fmat(1,1) + 0.5_dp * pmat(1,1)*gss
    END IF
    IF ( n == 4 ) THEN
      fmat(1,1) = fmat(1,1) + 0.5_dp * pmat(1,1)*gss + &
         (pmat(2,2)+pmat(3,3)+pmat(4,4))*(gsp-0.5_dp*hsp)
      fmat(2,2) = fmat(2,2) + pmat(1,1)*(gsp-0.5_dp*hsp) + &
                              0.5_dp * pmat(2,2)*gpp + &
                              (pmat(3,3)+pmat(4,4))*gp2 - &
                              0.25_dp * (pmat(3,3)+pmat(4,4))*(gpp-gp2)
      fmat(3,3) = fmat(3,3) + pmat(1,1)*(gsp-0.5_dp*hsp) + &
                              0.5_dp * pmat(3,3)*gpp + &
                              (pmat(2,2)+pmat(4,4))*gp2 - &
                              0.25_dp * (pmat(2,2)+pmat(4,4))*(gpp-gp2)
      fmat(4,4) = fmat(4,4) + pmat(1,1)*(gsp-0.5_dp*hsp) + &
                              0.5_dp * pmat(4,4)*gpp + &
                              (pmat(2,2)+pmat(3,3))*gp2 - &
                              0.25_dp * (pmat(2,2)+pmat(3,3))*(gpp-gp2)

      fmat(1,2) = fmat(1,2) + 0.5_dp * pmat(1,2)*(3._dp*hsp-gsp)
      fmat(2,1) = fmat(1,2)
      fmat(1,3) = fmat(1,3) + 0.5_dp * pmat(1,3)*(3._dp*hsp-gsp)
      fmat(3,1) = fmat(1,3)
      fmat(1,4) = fmat(1,4) + 0.5_dp * pmat(1,4)*(3._dp*hsp-gsp)
      fmat(4,1) = fmat(1,4)

      fmat(2,3) = fmat(2,3) + pmat(2,3)*(gpp-gp2) - 0.25_dp * pmat(2,3)*(gpp+gp2)
      fmat(3,2) = fmat(2,3)
      fmat(2,4) = fmat(2,4) + pmat(2,4)*(gpp-gp2) - 0.25_dp * pmat(2,4)*(gpp+gp2)
      fmat(4,2) = fmat(2,4)
      fmat(3,4) = fmat(3,4) + pmat(3,4)*(gpp-gp2) - 0.25_dp * pmat(3,4)*(gpp+gp2)
      fmat(4,3) = fmat(3,4)
    END IF

  END SUBROUTINE rfock1

  SUBROUTINE ufock1 (gss,gsp,gpp,gp2,hsp,pmata,fmata,pmatb,fmatb,error)
    REAL(dp), INTENT(IN)                     :: gss, gsp, gpp, gp2, hsp
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pmata
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fmata
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pmatb
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fmatb
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ufock1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n
    LOGICAL                                  :: failure
    REAL(dp)                                 :: pss, pxx, pyy, pzz

    n = SIZE ( pmata,1 )

    CPPrecondition(n<=4,cp_failure_level,routineP,error,failure)
    IF ( n == 1 ) THEN
      fmata(1,1) = fmata(1,1) + pmatb(1,1)*gss
      fmatb(1,1) = fmatb(1,1) + pmata(1,1)*gss
    END IF
    IF ( n == 4 ) THEN
      pss = pmata(1,1)+pmatb(1,1)
      pxx = pmata(2,2)+pmatb(2,2)
      pyy = pmata(3,3)+pmatb(3,3)
      pzz = pmata(4,4)+pmatb(4,4)
      fmata(1,1) = fmata(1,1) + pmatb(1,1)*gss + &
         (pxx+pyy+pzz)*gsp - (pmata(2,2)+pmata(3,3)+pmata(4,4))*hsp
      fmatb(1,1) = fmatb(1,1) + pmata(1,1)*gss + &
         (pxx+pyy+pzz)*gsp - (pmatb(2,2)+pmatb(3,3)+pmatb(4,4))*hsp
      fmata(2,2) = fmata(2,2) + pmatb(2,2)*gpp - pmata(1,1)*hsp + pss*gsp + (pyy+pzz)*gp2 -&
                    0.5_dp*(pmata(3,3)+pmata(4,4))*(gpp-gp2)
      fmata(3,3) = fmata(3,3) + pmatb(3,3)*gpp - pmata(1,1)*hsp + pss*gsp + (pxx+pzz)*gp2 -&
                    0.5_dp*(pmata(2,2)+pmata(4,4))*(gpp-gp2)
      fmata(4,4) = fmata(4,4) + pmatb(4,4)*gpp - pmata(1,1)*hsp + pss*gsp + (pxx+pyy)*gp2 -&
                    0.5_dp*(pmata(2,2)+pmata(3,3))*(gpp-gp2)
      fmatb(2,2) = fmatb(2,2) + pmata(2,2)*gpp - pmatb(1,1)*hsp + pss*gsp + (pyy+pzz)*gp2 -&
                    0.5_dp*(pmatb(3,3)+pmatb(4,4))*(gpp-gp2)
      fmatb(3,3) = fmatb(3,3) + pmata(3,3)*gpp - pmatb(1,1)*hsp + pss*gsp + (pxx+pzz)*gp2 -&
                    0.5_dp*(pmatb(2,2)+pmatb(4,4))*(gpp-gp2)
      fmatb(4,4) = fmatb(4,4) + pmata(4,4)*gpp - pmatb(1,1)*hsp + pss*gsp + (pxx+pyy)*gp2 -&
                    0.5_dp*(pmatb(2,2)+pmatb(3,3))*(gpp-gp2)

      fmata(1,2) = fmata(1,2) + pmata(1,2)*(hsp-gsp) + 2._dp*pmatb(1,2)*hsp
      fmata(2,1) = fmata(1,2)
      fmata(1,3) = fmata(1,3) + pmata(1,3)*(hsp-gsp) + 2._dp*pmatb(1,3)*hsp
      fmata(3,1) = fmata(1,3)
      fmata(1,4) = fmata(1,4) + pmata(1,4)*(hsp-gsp) + 2._dp*pmatb(1,4)*hsp
      fmata(4,1) = fmata(1,4)
      fmatb(1,2) = fmatb(1,2) + pmatb(1,2)*(hsp-gsp) + 2._dp*pmata(1,2)*hsp
      fmatb(2,1) = fmatb(1,2)
      fmatb(1,3) = fmatb(1,3) + pmatb(1,3)*(hsp-gsp) + 2._dp*pmata(1,3)*hsp
      fmatb(3,1) = fmatb(1,3)
      fmatb(1,4) = fmatb(1,4) + pmatb(1,4)*(hsp-gsp) + 2._dp*pmata(1,4)*hsp
      fmatb(4,1) = fmatb(1,4)

      fmata(2,3) = fmata(2,3) + (pmata(2,3)+pmatb(2,3))*(gpp-gp2) &
                              - 0.5_dp * pmata(2,3)*(gpp+gp2)
      fmata(3,2) = fmata(2,3)
      fmata(2,4) = fmata(2,4) + (pmata(2,4)+pmatb(2,4))*(gpp-gp2) &
                              - 0.5_dp * pmata(2,4)*(gpp+gp2)
      fmata(4,2) = fmata(2,4)
      fmata(3,4) = fmata(3,4) + (pmata(3,4)+pmatb(3,4))*(gpp-gp2) &
                              - 0.5_dp * pmata(3,4)*(gpp+gp2)
      fmata(4,3) = fmata(3,4)
      fmatb(2,3) = fmatb(2,3) + (pmata(2,3)+pmatb(2,3))*(gpp-gp2) &
                              - 0.5_dp * pmatb(2,3)*(gpp+gp2)
      fmatb(3,2) = fmatb(2,3)
      fmatb(2,4) = fmatb(2,4) + (pmata(2,4)+pmatb(2,4))*(gpp-gp2) &
                              - 0.5_dp * pmatb(2,4)*(gpp+gp2)
      fmatb(4,2) = fmatb(2,4)
      fmatb(3,4) = fmatb(3,4) + (pmata(3,4)+pmatb(3,4))*(gpp-gp2) &
                              - 0.5_dp * pmatb(3,4)*(gpp+gp2)
      fmatb(4,3) = fmatb(3,4)

    END IF
   
  END SUBROUTINE ufock1

! *****************************************************************************

  SUBROUTINE rfock2c(sepa,sepb,rij,pa,pb,fa,fb,error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pa, pb
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fa, fb
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rfock2c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(10)             :: ca, cb, ppa, ppb
    REAL(KIND=dp), DIMENSION(100)            :: wint

    na = SIZE ( pa, 1 )
    nb = SIZE ( pb, 1 )

    CALL rotint (sepa,sepb,rij,wint)

    IF ( na==1 .AND. nb==1 ) THEN
      fa(1,1) = fa(1,1) + pb(1,1)*wint(1)
      fb(1,1) = fb(1,1) + pa(1,1)*wint(1)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb, ppb )
      fa(1,1) = fa(1,1) + SUM ( ppb(1:10)*wint(1:10) )
      ppa(1:10) = pa(1,1)*wint(1:10)
      CALL expand( fb, ppa )
    ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa, ppa )
      fb(1,1) = fb(1,1) + SUM ( ppa(1:10)*wint(1:10) )
      ppb(1:10) = pb(1,1)*wint(1:10)
      CALL expand( fa, ppb )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pa, ppa )
      CALL contract( pb, ppb )
      DO i=1,10
        ii = (i-1)*10 + 1
        ca(i) = SUM ( wint(ii:ii+9)*ppb(1:10) )
        cb(i) = SUM ( wint(i:90+i:10)*ppa(1:10) )
      END DO
      CALL expand( fa, ca )
      CALL expand( fb, cb )
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE rfock2c

  SUBROUTINE ufock2c(sepa,sepb,rij,pa_a,pb_a,fa_a,fb_a,&
                    pa_b,pb_b,fa_b,fb_b,error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pa_a, pb_a
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fa_a, fb_a
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pa_b, pb_b
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fa_b, fb_b
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ufock2c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ess
    REAL(KIND=dp), DIMENSION(10)             :: ca, cb, ppa_a, ppa_b, ppb_a, &
                                                ppb_b
    REAL(KIND=dp), DIMENSION(100)            :: wint

    na = SIZE ( pa_a, 1 )
    nb = SIZE ( pb_a, 1 )

    CALL rotint (sepa,sepb,rij,wint)

    IF ( na==1 .AND. nb==1 ) THEN
      fa_a(1,1) = fa_a(1,1) + (pb_a(1,1)+pb_b(1,1))*wint(1)
      fb_a(1,1) = fb_a(1,1) + (pa_a(1,1)+pa_b(1,1))*wint(1)
      fa_b(1,1) = fa_b(1,1) + (pb_a(1,1)+pb_b(1,1))*wint(1)
      fb_b(1,1) = fb_b(1,1) + (pa_a(1,1)+pa_b(1,1))*wint(1)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      ess = SUM ( (ppb_a(1:10)+ppb_b(1:10)) * wint(1:10) )
      fa_a(1,1) = fa_a(1,1) + ess
      fa_b(1,1) = fa_b(1,1) + ess
      ppa_a(1:10) = (pa_a(1,1)+pa_b(1,1))*wint(1:10)
      CALL expand( fb_a, ppa_a )
      CALL expand( fb_b, ppa_a )
    ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      ess = SUM ( (ppa_a(1:10)+ppa_b(1:10)) * wint(1:10) )
      fb_a(1,1) = fb_a(1,1) + ess
      fb_b(1,1) = fb_b(1,1) + ess
      ppb_a(1:10) = (pb_a(1,1)+pb_b(1,1))*wint(1:10)
      CALL expand( fa_a, ppb_a )
      CALL expand( fa_b, ppb_a )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      DO i=1,10
        ii = (i-1)*10 + 1
        ca(i) = SUM ( wint(ii:ii+9)*(ppb_a(1:10)+ppb_b(1:10)) )
        cb(i) = SUM ( wint(i:90+i:10)*(ppa_a(1:10)+ppa_b(1:10)) )
      END DO
      CALL expand( fa_a, ca )
      CALL expand( fa_b, ca )
      CALL expand( fb_a, cb )
      CALL expand( fb_b, cb )
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE ufock2c

! *****************************************************************************

  SUBROUTINE rdfock2c(sepa,sepb,rij,pa,pb,force,delta,error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pa, pb
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: force
    REAL(dp), INTENT(IN)                     :: delta
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rdfock2c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, j, na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(10)             :: ca, ppa, ppb
    REAL(KIND=dp), DIMENSION(100, 3)         :: dwint

    na = SIZE ( pa, 1 )
    nb = SIZE ( pb, 1 )

    CALL drotint (sepa,sepb,rij,dwint,delta)

    IF ( na==1 .AND. nb==1 ) THEN
      force(:) = force(:) + pa(1,1)*pb(1,1)*dwint(1,:)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb, ppb )
      DO i=1,3
        force(i) = force(i) + pa(1,1) * SUM ( ppb(1:10)*dwint(1:10,i) )
      END DO
   ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa, ppa )
      DO i=1,3
        force(i) = force(i) + pb(1,1) * SUM ( ppa(1:10)*dwint(1:10,i) )
      END DO
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pb, ppb )
      CALL contract( pa, ppa )
      DO j=1,3
        DO i=1,10
          ii = (i-1)*10 + 1
          ca(i) = SUM ( dwint(ii:ii+9,j)*ppb(1:10) )
        END DO
        force(j) = force(j) + SUM ( ca(1:10)*ppa(1:10) )
      END DO
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE rdfock2c

  SUBROUTINE udfock2c(sepa,sepb,rij,pa_a,pb_a,pa_b,pb_b,&
             force,delta,error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pa_a, pb_a, pa_b, pb_b
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: force
    REAL(dp), INTENT(IN)                     :: delta
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'udfock2c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, j, na, nb
    LOGICAL                                  :: failure
    REAL(dp)                                 :: ess, pta, ptb
    REAL(KIND=dp), DIMENSION(10)             :: ca, ppa_a, ppa_b, ppb_a, ppb_b
    REAL(KIND=dp), DIMENSION(100, 3)         :: dwint

    na = SIZE ( pa_a, 1 )
    nb = SIZE ( pb_a, 1 )

    CALL drotint (sepa,sepb,rij,dwint,delta)

    IF ( na==1 .AND. nb==1 ) THEN
      pta = pa_a(1,1)+pa_b(1,1)
      ptb = pb_a(1,1)+pb_b(1,1)
      force(:) = force(:) + pta*ptb*dwint(1,:)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      DO i=1,3
        ess = SUM ( (ppb_a(1:10)+ppb_b(1:10)) * dwint(1:10,i) )
        force(i) = force(i) + (pa_a(1,1)+pa_b(1,1)) * ess
      END DO
   ELSE IF ( na==4 .AND. nb==1 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      DO i=1,3
        ess = SUM ( (ppa_a(1:10)+ppa_b(1:10)) * dwint(1:10,i) )
        force(i) = force(i) + (pb_a(1,1)+pb_b(1,1)) * ess
      END DO
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL contract( pa_a, ppa_a )
      CALL contract( pa_b, ppa_b )
      CALL contract( pb_a, ppb_a )
      CALL contract( pb_b, ppb_b )
      DO j=1,3
        DO i=1,10
          ii = (i-1)*10 + 1
          ca(i) = SUM ( dwint(ii:ii+9,j)*(ppb_a(1:10)+ppb_b(1:10)) )
        END DO
        force(j) = force(j) + SUM ( ca(1:10)*(ppa_a(1:10)+ppa_b(1:10)) )
      END DO
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE udfock2c

  SUBROUTINE rfock2e(sepa,sepb,rij,pab,fab,error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fab
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rfock2e', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(100)            :: wint
    REAL(KIND=dp), DIMENSION(4, 4)           :: kab

    na = SIZE ( pab, 1 )
    nb = SIZE ( pab, 2 )

    CALL rotint (sepa,sepb,rij,wint)

    IF ( na==1 .AND. nb==1 ) THEN
      fab(1,1) = fab(1,1) - 0.5_dp * pab(1,1)*wint(1)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      fab(1,1) = fab(1,1) - 0.5_dp*( wint(1)*pab(1,1) +&
        wint(2)*pab(1,4) + wint(4)*pab(1,2) + wint(7)*pab(1,3) )
      fab(1,2) = fab(1,2) - 0.5_dp*( wint(4)*pab(1,1) +&
        wint(5)*pab(1,4) + wint(6)*pab(1,2) + wint(9)*pab(1,3) )
      fab(1,3) = fab(1,3) - 0.5_dp*( wint(7)*pab(1,1) +&
        wint(8)*pab(1,4) + wint(9)*pab(1,2) + wint(10)*pab(1,3) )
      fab(1,4) = fab(1,4) - 0.5_dp*( wint(2)*pab(1,1) +&
        wint(3)*pab(1,4) + wint(5)*pab(1,2) + wint(8)*pab(1,3) )
    ELSE IF ( na==4 .AND. nb==1 ) THEN
      fab(1,1) = fab(1,1) - 0.5_dp*( wint(1)*pab(1,1) +&
        wint(2)*pab(4,1) + wint(4)*pab(2,1) + wint(7)*pab(3,1) )
      fab(2,1) = fab(2,1) - 0.5_dp*( wint(4)*pab(1,1) +&
        wint(5)*pab(4,1) + wint(6)*pab(2,1) + wint(9)*pab(3,1) )
      fab(3,1) = fab(3,1) - 0.5_dp*( wint(7)*pab(1,1) +&
        wint(8)*pab(4,1) + wint(9)*pab(2,1) + wint(10)*pab(3,1) )
      fab(4,1) = fab(4,1) - 0.5_dp*( wint(2)*pab(1,1) +&
        wint(3)*pab(4,1) + wint(5)*pab(2,1) + wint(8)*pab(3,1) )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL exchange_sp ( kab, pab, wint )
      fab(1:4,1:4) = fab(1:4,1:4) - 0.5_dp*kab(1:4,1:4)
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE rfock2e

  SUBROUTINE ufock2e(sepa,sepb,rij,pab_a,fab_a,pab_b,fab_b,error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab_a
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fab_a
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab_b
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: fab_b
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ufock2e', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(100)            :: wint
    REAL(KIND=dp), DIMENSION(4, 4)           :: kab

    na = SIZE ( pab_a, 1 )
    nb = SIZE ( pab_a, 2 )

    CALL rotint (sepa,sepb,rij,wint)

    IF ( na==1 .AND. nb==1 ) THEN
      fab_a(1,1) = fab_a(1,1) - pab_a(1,1)*wint(1)
      fab_b(1,1) = fab_b(1,1) - pab_b(1,1)*wint(1)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      fab_a(1,1) = fab_a(1,1) - wint(1)*pab_a(1,1) -&
        wint(2)*pab_a(1,4) - wint(4)*pab_a(1,2) - wint(7)*pab_a(1,3)
      fab_a(1,2) = fab_a(1,2) - wint(4)*pab_a(1,1) -&
        wint(5)*pab_a(1,4) - wint(6)*pab_a(1,2) - wint(9)*pab_a(1,3)
      fab_a(1,3) = fab_a(1,3) - wint(7)*pab_a(1,1) -&
        wint(8)*pab_a(1,4) - wint(9)*pab_a(1,2) - wint(10)*pab_a(1,3)
      fab_a(1,4) = fab_a(1,4) - wint(2)*pab_a(1,1) -&
        wint(3)*pab_a(1,4) - wint(5)*pab_a(1,2) - wint(8)*pab_a(1,3)
      fab_b(1,1) = fab_b(1,1) - wint(1)*pab_b(1,1) -&
        wint(2)*pab_b(1,4) - wint(4)*pab_b(1,2) - wint(7)*pab_b(1,3)
      fab_b(1,2) = fab_b(1,2) - wint(4)*pab_b(1,1) -&
        wint(5)*pab_b(1,4) - wint(6)*pab_b(1,2) - wint(9)*pab_b(1,3)
      fab_b(1,3) = fab_b(1,3) - wint(7)*pab_b(1,1) -&
        wint(8)*pab_b(1,4) - wint(9)*pab_b(1,2) - wint(10)*pab_b(1,3)
      fab_b(1,4) = fab_b(1,4) - wint(2)*pab_b(1,1) -&
        wint(3)*pab_b(1,4) - wint(5)*pab_b(1,2) - wint(8)*pab_b(1,3)
    ELSE IF ( na==4 .AND. nb==1 ) THEN
      fab_a(1,1) = fab_a(1,1) - wint(1)*pab_a(1,1) -&
        wint(2)*pab_a(4,1) - wint(4)*pab_a(2,1) - wint(7)*pab_a(3,1)
      fab_a(2,1) = fab_a(2,1) - wint(4)*pab_a(1,1) -&
        wint(5)*pab_a(4,1) - wint(6)*pab_a(2,1) - wint(9)*pab_a(3,1)
      fab_a(3,1) = fab_a(3,1) - wint(7)*pab_a(1,1) -&
        wint(8)*pab_a(4,1) - wint(9)*pab_a(2,1) - wint(10)*pab_a(3,1)
      fab_a(4,1) = fab_a(4,1) - wint(2)*pab_a(1,1) -&
        wint(3)*pab_a(4,1) - wint(5)*pab_a(2,1) - wint(8)*pab_a(3,1)
      fab_b(1,1) = fab_b(1,1) - wint(1)*pab_b(1,1) -&
        wint(2)*pab_b(4,1) - wint(4)*pab_b(2,1) - wint(7)*pab_b(3,1)
      fab_b(2,1) = fab_b(2,1) - wint(4)*pab_b(1,1) -&
        wint(5)*pab_b(4,1) - wint(6)*pab_b(2,1) - wint(9)*pab_b(3,1)
      fab_b(3,1) = fab_b(3,1) - wint(7)*pab_b(1,1) -&
        wint(8)*pab_b(4,1) - wint(9)*pab_b(2,1) - wint(10)*pab_b(3,1)
      fab_b(4,1) = fab_b(4,1) - wint(2)*pab_b(1,1) -&
        wint(3)*pab_b(4,1) - wint(5)*pab_b(2,1) - wint(8)*pab_b(3,1)
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      CALL exchange_sp ( kab, pab_a, wint )
      fab_a(1:4,1:4) = fab_a(1:4,1:4) - kab(1:4,1:4)
      CALL exchange_sp ( kab, pab_b, wint )
      fab_b(1:4,1:4) = fab_b(1:4,1:4) - kab(1:4,1:4)
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE ufock2e

! *****************************************************************************

  SUBROUTINE rdfock2e(sepa,sepb,rij,pab,force,delta,error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: force
    REAL(dp), INTENT(IN)                     :: delta
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rdfock2e', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(100, 3)         :: dwint
    REAL(KIND=dp), DIMENSION(4, 4)           :: kab

    na = SIZE ( pab, 1 )
    nb = SIZE ( pab, 2 )

    CALL drotint (sepa,sepb,rij,dwint,delta)

    IF ( na==1 .AND. nb==1 ) THEN
      force(:) = force(:) - 0.5_dp*pab(1,1)*pab(1,1)*dwint(1,:)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      force(:) = force(:) - 0.5_dp*pab(1,1)*( dwint(1,:)*pab(1,1) +&
        dwint(2,:)*pab(1,4) + dwint(4,:)*pab(1,2) + dwint(7,:)*pab(1,3) )
      force(:) = force(:) - 0.5_dp*pab(1,2)*( dwint(4,:)*pab(1,1) +&
        dwint(5,:)*pab(1,4) + dwint(6,:)*pab(1,2) + dwint(9,:)*pab(1,3) )
      force(:) = force(:) - 0.5_dp*pab(1,3)*( dwint(7,:)*pab(1,1) +&
        dwint(8,:)*pab(1,4) + dwint(9,:)*pab(1,2) + dwint(10,:)*pab(1,3) )
      force(:) = force(:) - 0.5_dp*pab(1,4)*( dwint(2,:)*pab(1,1) +&
        dwint(3,:)*pab(1,4) + dwint(5,:)*pab(1,2) + dwint(8,:)*pab(1,3) )
   ELSE IF ( na==4 .AND. nb==1 ) THEN
      force(:) = force(:) - 0.5_dp*pab(1,1)*( dwint(1,:)*pab(1,1) +&
        dwint(2,:)*pab(4,1) + dwint(4,:)*pab(2,1) + dwint(7,:)*pab(3,1) )
      force(:) = force(:) - 0.5_dp*pab(2,1)*( dwint(4,:)*pab(1,1) +&
        dwint(5,:)*pab(4,1) + dwint(6,:)*pab(2,1) + dwint(9,:)*pab(3,1) )
      force(:) = force(:) - 0.5_dp*pab(3,1)*( dwint(7,:)*pab(1,1) +&
        dwint(8,:)*pab(4,1) + dwint(9,:)*pab(2,1) + dwint(10,:)*pab(3,1) )
      force(:) = force(:) - 0.5_dp*pab(4,1)*( dwint(2,:)*pab(1,1) +&
        dwint(3,:)*pab(4,1) + dwint(5,:)*pab(2,1) + dwint(8,:)*pab(3,1) )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      DO j=1,3
        CALL exchange_sp ( kab, pab, dwint(:,j) )
        force(j) = force(j) - 0.5_dp * SUM(pab*kab)
      END DO
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE rdfock2e

  SUBROUTINE udfock2e(sepa,sepb,rij,pab_a,pab_b,force,delta,error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepa, sepb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rij
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab_a, pab_b
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: force
    REAL(dp), INTENT(IN)                     :: delta
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'udfock2e', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, na, nb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(100, 3)         :: dwint
    REAL(KIND=dp), DIMENSION(4, 4)           :: kab

    na = SIZE ( pab_a, 1 )
    nb = SIZE ( pab_a, 2 )

    CALL drotint (sepa,sepb,rij,dwint,delta)

    IF ( na==1 .AND. nb==1 ) THEN
      force(:) = force(:) - pab_a(1,1)*pab_a(1,1)*dwint(1,:)
      force(:) = force(:) - pab_b(1,1)*pab_b(1,1)*dwint(1,:)
    ELSE IF ( na==1 .AND. nb==4 ) THEN
      force(:) = force(:) - pab_a(1,1)*( dwint(1,:)*pab_a(1,1) +&
        dwint(2,:)*pab_a(1,4) + dwint(4,:)*pab_a(1,2) + dwint(7,:)*pab_a(1,3) )
      force(:) = force(:) - pab_a(1,2)*( dwint(4,:)*pab_a(1,1) +&
        dwint(5,:)*pab_a(1,4) + dwint(6,:)*pab_a(1,2) + dwint(9,:)*pab_a(1,3) )
      force(:) = force(:) - pab_a(1,3)*( dwint(7,:)*pab_a(1,1) +&
        dwint(8,:)*pab_a(1,4) + dwint(9,:)*pab_a(1,2) + dwint(10,:)*pab_a(1,3) )
      force(:) = force(:) - pab_a(1,4)*( dwint(2,:)*pab_a(1,1) +&
        dwint(3,:)*pab_a(1,4) + dwint(5,:)*pab_a(1,2) + dwint(8,:)*pab_a(1,3) )
      force(:) = force(:) - pab_b(1,1)*( dwint(1,:)*pab_b(1,1) +&
        dwint(2,:)*pab_b(1,4) + dwint(4,:)*pab_b(1,2) + dwint(7,:)*pab_b(1,3) )
      force(:) = force(:) - pab_b(1,2)*( dwint(4,:)*pab_b(1,1) +&
        dwint(5,:)*pab_b(1,4) + dwint(6,:)*pab_b(1,2) + dwint(9,:)*pab_b(1,3) )
      force(:) = force(:) - pab_b(1,3)*( dwint(7,:)*pab_b(1,1) +&
        dwint(8,:)*pab_b(1,4) + dwint(9,:)*pab_b(1,2) + dwint(10,:)*pab_b(1,3) )
      force(:) = force(:) - pab_b(1,4)*( dwint(2,:)*pab_b(1,1) +&
        dwint(3,:)*pab_b(1,4) + dwint(5,:)*pab_b(1,2) + dwint(8,:)*pab_b(1,3) )
   ELSE IF ( na==4 .AND. nb==1 ) THEN
      force(:) = force(:) - pab_a(1,1)*( dwint(1,:)*pab_a(1,1) +&
        dwint(2,:)*pab_a(4,1) + dwint(4,:)*pab_a(2,1) + dwint(7,:)*pab_a(3,1) )
      force(:) = force(:) - pab_a(2,1)*( dwint(4,:)*pab_a(1,1) +&
        dwint(5,:)*pab_a(4,1) + dwint(6,:)*pab_a(2,1) + dwint(9,:)*pab_a(3,1) )
      force(:) = force(:) - pab_a(3,1)*( dwint(7,:)*pab_a(1,1) +&
        dwint(8,:)*pab_a(4,1) + dwint(9,:)*pab_a(2,1) + dwint(10,:)*pab_a(3,1) )
      force(:) = force(:) - pab_a(4,1)*( dwint(2,:)*pab_a(1,1) +&
        dwint(3,:)*pab_a(4,1) + dwint(5,:)*pab_a(2,1) + dwint(8,:)*pab_a(3,1) )
      force(:) = force(:) - pab_b(1,1)*( dwint(1,:)*pab_b(1,1) +&
        dwint(2,:)*pab_b(4,1) + dwint(4,:)*pab_b(2,1) + dwint(7,:)*pab_b(3,1) )
      force(:) = force(:) - pab_b(2,1)*( dwint(4,:)*pab_b(1,1) +&
        dwint(5,:)*pab_b(4,1) + dwint(6,:)*pab_b(2,1) + dwint(9,:)*pab_b(3,1) )
      force(:) = force(:) - pab_b(3,1)*( dwint(7,:)*pab_b(1,1) +&
        dwint(8,:)*pab_b(4,1) + dwint(9,:)*pab_b(2,1) + dwint(10,:)*pab_b(3,1) )
      force(:) = force(:) - pab_b(4,1)*( dwint(2,:)*pab_b(1,1) +&
        dwint(3,:)*pab_b(4,1) + dwint(5,:)*pab_b(2,1) + dwint(8,:)*pab_b(3,1) )
    ELSE IF ( na==4 .AND. nb==4 ) THEN
      DO j=1,3
        CALL exchange_sp ( kab, pab_a, dwint(:,j) )
        force(j) = force(j) - SUM(pab_a*kab)
        CALL exchange_sp ( kab, pab_b, dwint(:,j) )
        force(j) = force(j) - SUM(pab_b*kab)
      END DO
    ELSE
      CPPrecondition(na<=4,cp_failure_level,routineP,error,failure)
      CPPrecondition(nb<=4,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE udfock2e

! *****************************************************************************

  SUBROUTINE exchange_sp ( kab, pab, w )
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: kab
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab
    REAL(dp), DIMENSION(:), INTENT(IN)       :: w

    REAL(dp), DIMENSION(16)                  :: pk

    pk( 1) = pab(1,1)
    pk( 2) = pab(1,4)
    pk( 3) = pab(1,2)
    pk( 4) = pab(1,3)
    pk( 5) = pab(4,1)
    pk( 6) = pab(4,4)
    pk( 7) = pab(4,2)
    pk( 8) = pab(4,3)
    pk( 9) = pab(2,1)
    pk(10) = pab(2,4)
    pk(11) = pab(2,2)
    pk(12) = pab(2,3)
    pk(13) = pab(3,1)
    pk(14) = pab(3,4)
    pk(15) = pab(3,2)
    pk(16) = pab(3,3)

    kab(1,1)= &
      +pk( 1)*w(  1)+pk( 2)*w(  2)+pk( 3)*w(  4)+pk( 4)*w(  7) &
      +pk( 5)*w( 11)+pk( 6)*w( 12)+pk( 7)*w( 14)+pk( 8)*w( 17) &
      +pk( 9)*w( 31)+pk(10)*w( 32)+pk(11)*w( 34)+pk(12)*w( 37) &
      +pk(13)*w( 61)+pk(14)*w( 62)+pk(15)*w( 64)+pk(16)*w( 67)
    kab(1,4)= &
      +pk( 1)*w(  2)+pk( 2)*w(  3)+pk( 3)*w(  5)+pk( 4)*w(  8) &
      +pk( 5)*w( 12)+pk( 6)*w( 13)+pk( 7)*w( 15)+pk( 8)*w( 18) &
      +pk( 9)*w( 32)+pk(10)*w( 33)+pk(11)*w( 35)+pk(12)*w( 38) &
      +pk(13)*w( 62)+pk(14)*w( 63)+pk(15)*w( 65)+pk(16)*w( 68)
    kab(1,2)= &
      +pk( 1)*w(  4)+pk( 2)*w(  5)+pk( 3)*w(  6)+pk( 4)*w(  9) &
      +pk( 5)*w( 14)+pk( 6)*w( 15)+pk( 7)*w( 16)+pk( 8)*w( 19) &
      +pk( 9)*w( 34)+pk(10)*w( 35)+pk(11)*w( 36)+pk(12)*w( 39) &
      +pk(13)*w( 64)+pk(14)*w( 65)+pk(15)*w( 66)+pk(16)*w( 69)
    kab(1,3)= &
      +pk( 1)*w(  7)+pk( 2)*w(  8)+pk( 3)*w(  9)+pk( 4)*w( 10) &
      +pk( 5)*w( 17)+pk( 6)*w( 18)+pk( 7)*w( 19)+pk( 8)*w( 20) &
      +pk( 9)*w( 37)+pk(10)*w( 38)+pk(11)*w( 39)+pk(12)*w( 40) &
      +pk(13)*w( 67)+pk(14)*w( 68)+pk(15)*w( 69)+pk(16)*w( 70)
    kab(4,1)= &
      +pk( 1)*w( 11)+pk( 2)*w( 12)+pk( 3)*w( 14)+pk( 4)*w( 17) &
      +pk( 5)*w( 21)+pk( 6)*w( 22)+pk( 7)*w( 24)+pk( 8)*w( 27) &
      +pk( 9)*w( 41)+pk(10)*w( 42)+pk(11)*w( 44)+pk(12)*w( 47) &
      +pk(13)*w( 71)+pk(14)*w( 72)+pk(15)*w( 74)+pk(16)*w( 77)
    kab(4,4)= &
      +pk( 1)*w( 12)+pk( 2)*w( 13)+pk( 3)*w( 15)+pk( 4)*w( 18) &
      +pk( 5)*w( 22)+pk( 6)*w( 23)+pk( 7)*w( 25)+pk( 8)*w( 28) &
      +pk( 9)*w( 42)+pk(10)*w( 43)+pk(11)*w( 45)+pk(12)*w( 48) &
      +pk(13)*w( 72)+pk(14)*w( 73)+pk(15)*w( 75)+pk(16)*w( 78)
    kab(4,2)= &
      +pk( 1)*w( 14)+pk( 2)*w( 15)+pk( 3)*w( 16)+pk( 4)*w( 19) &
      +pk( 5)*w( 24)+pk( 6)*w( 25)+pk( 7)*w( 26)+pk( 8)*w( 29) &
      +pk( 9)*w( 44)+pk(10)*w( 45)+pk(11)*w( 46)+pk(12)*w( 49) &
      +pk(13)*w( 74)+pk(14)*w( 75)+pk(15)*w( 76)+pk(16)*w( 79)
    kab(4,3)= &
      +pk( 1)*w( 17)+pk( 2)*w( 18)+pk( 3)*w( 19)+pk( 4)*w( 20) &
      +pk( 5)*w( 27)+pk( 6)*w( 28)+pk( 7)*w( 29)+pk( 8)*w( 30) &
      +pk( 9)*w( 47)+pk(10)*w( 48)+pk(11)*w( 49)+pk(12)*w( 50) &
      +pk(13)*w( 77)+pk(14)*w( 78)+pk(15)*w( 79)+pk(16)*w( 80)
    kab(2,1)= &
      +pk( 1)*w( 31)+pk( 2)*w( 32)+pk( 3)*w( 34)+pk( 4)*w( 37) &
      +pk( 5)*w( 41)+pk( 6)*w( 42)+pk( 7)*w( 44)+pk( 8)*w( 47) &
      +pk( 9)*w( 51)+pk(10)*w( 52)+pk(11)*w( 54)+pk(12)*w( 57) &
      +pk(13)*w( 81)+pk(14)*w( 82)+pk(15)*w( 84)+pk(16)*w( 87)
    kab(2,4)= &
      +pk( 1)*w( 32)+pk( 2)*w( 33)+pk( 3)*w( 35)+pk( 4)*w( 38) &
      +pk( 5)*w( 42)+pk( 6)*w( 43)+pk( 7)*w( 45)+pk( 8)*w( 48) &
      +pk( 9)*w( 52)+pk(10)*w( 53)+pk(11)*w( 55)+pk(12)*w( 58) &
      +pk(13)*w( 82)+pk(14)*w( 83)+pk(15)*w( 85)+pk(16)*w( 88)
    kab(2,2)= &
      +pk( 1)*w( 34)+pk( 2)*w( 35)+pk( 3)*w( 36)+pk( 4)*w( 39) &
      +pk( 5)*w( 44)+pk( 6)*w( 45)+pk( 7)*w( 46)+pk( 8)*w( 49) &
      +pk( 9)*w( 54)+pk(10)*w( 55)+pk(11)*w( 56)+pk(12)*w( 59) &
      +pk(13)*w( 84)+pk(14)*w( 85)+pk(15)*w( 86)+pk(16)*w( 89)
    kab(2,3)= &
      +pk( 1)*w( 37)+pk( 2)*w( 38)+pk( 3)*w( 39)+pk( 4)*w( 40) &
      +pk( 5)*w( 47)+pk( 6)*w( 48)+pk( 7)*w( 49)+pk( 8)*w( 50) &
      +pk( 9)*w( 57)+pk(10)*w( 58)+pk(11)*w( 59)+pk(12)*w( 60) &
      +pk(13)*w( 87)+pk(14)*w( 88)+pk(15)*w( 89)+pk(16)*w( 90)
    kab(3,1)= &
      +pk( 1)*w( 61)+pk( 2)*w( 62)+pk( 3)*w( 64)+pk( 4)*w( 67) &
      +pk( 5)*w( 71)+pk( 6)*w( 72)+pk( 7)*w( 74)+pk( 8)*w( 77) &
      +pk( 9)*w( 81)+pk(10)*w( 82)+pk(11)*w( 84)+pk(12)*w( 87) &
      +pk(13)*w( 91)+pk(14)*w( 92)+pk(15)*w( 94)+pk(16)*w( 97)
    kab(3,4)= &
      +pk( 1)*w( 62)+pk( 2)*w( 63)+pk( 3)*w( 65)+pk( 4)*w( 68) &
      +pk( 5)*w( 72)+pk( 6)*w( 73)+pk( 7)*w( 75)+pk( 8)*w( 78) &
      +pk( 9)*w( 82)+pk(10)*w( 83)+pk(11)*w( 85)+pk(12)*w( 88) &
      +pk(13)*w( 92)+pk(14)*w( 93)+pk(15)*w( 95)+pk(16)*w( 98)
    kab(3,2)= &
      +pk( 1)*w( 64)+pk( 2)*w( 65)+pk( 3)*w( 66)+pk( 4)*w( 69) &
      +pk( 5)*w( 74)+pk( 6)*w( 75)+pk( 7)*w( 76)+pk( 8)*w( 79) &
      +pk( 9)*w( 84)+pk(10)*w( 85)+pk(11)*w( 86)+pk(12)*w( 89) &
      +pk(13)*w( 94)+pk(14)*w( 95)+pk(15)*w( 96)+pk(16)*w( 99)
    kab(3,3)= &
      +pk( 1)*w( 67)+pk( 2)*w( 68)+pk( 3)*w( 69)+pk( 4)*w( 70) &
      +pk( 5)*w( 77)+pk( 6)*w( 78)+pk( 7)*w( 79)+pk( 8)*w( 80) &
      +pk( 9)*w( 87)+pk(10)*w( 88)+pk(11)*w( 89)+pk(12)*w( 90) &
      +pk(13)*w( 97)+pk(14)*w( 98)+pk(15)*w( 99)+pk(16)*w(100)

  END SUBROUTINE exchange_sp

! *****************************************************************************

END MODULE nddo_methods

! *****************************************************************************

