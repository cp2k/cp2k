!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of Hamiltonian matrices in SCPTB
!> \author JGH
! *****************************************************************************
MODULE scptb_ks_matrix
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type,&
                                             scptb_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_iterator_blocks_left,&
                                             cp_dbcsr_iterator_next_block,&
                                             cp_dbcsr_iterator_start,&
                                             cp_dbcsr_iterator_stop,&
                                             cp_dbcsr_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE erf_fn,                          ONLY: erfc
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_sum
  USE mulliken,                        ONLY: mulliken_charges
  USE orbital_pointers,                ONLY: indso
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_transfer
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_collocate_density,            ONLY: calculate_scp_charge
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_integrate_potential,          ONLY: integrate_scp_rspace
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scptb_types,                     ONLY: &
       get_scptb_parameter, scp_vector_add, scp_vector_copy, &
       scp_vector_create, scp_vector_dot, scp_vector_release, &
       scp_vector_scale, scp_vector_set, scp_vector_type, &
       scptb_parameter_p_type, scptb_parameter_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scptb_ks_matrix'

  PUBLIC :: build_scptb_ks_matrix

CONTAINS

! *****************************************************************************
  SUBROUTINE build_scptb_ks_matrix (ks_env,qs_env,ks_matrix,rho,energy,&
               calculate_forces,just_energy,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_scptb_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ikind, &
                                                ispin, istat, natom, &
                                                natom_kind, nkind, nspins, &
                                                output_unit
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: zeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dqpair, mcharge
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation_numbers
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges, dmcharge
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_p, matrix_s, &
                                                mo_derivs
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind
    TYPE(section_vals_type), POINTER         :: scf_section

    CALL timeset(routineN,handle)

    NULLIFY(dft_control, logger, scf_section)
    NULLIFY(matrix_p, particle_set)

    logger => cp_error_get_logger(error)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN

       CALL get_qs_env(qs_env=qs_env,&
            dft_control=dft_control,&
            atomic_kind_set=atomic_kind_set,&
            matrix_h=matrix_h,&
            para_env=para_env,error=error)

       scf_section => section_vals_get_subs_vals(qs_env%input,"DFT%SCF",error=error)
       nspins=dft_control%nspins
       CPPrecondition(ASSOCIATED(matrix_h),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

       DO ispin=1,nspins
         ! copy the core matrix into the fock matrix
          CALL cp_dbcsr_copy(ks_matrix(ispin)%matrix,matrix_h(1)%matrix,error=error)
       END DO

      ! Mulliken charges
      CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,matrix_s=matrix_s,error=error)
      matrix_p => rho%rho_ao
      natom=SIZE(particle_set)
      ALLOCATE(charges(natom,nspins),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(mcharge(natom),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      IF(calculate_forces) THEN
         ALLOCATE(dmcharge(natom,3),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         CALL mulliken_charges(matrix_p,matrix_s,para_env,charges,dmcharge,error=error)
      ELSE
         NULLIFY(dmcharge)
         CALL mulliken_charges(matrix_p,matrix_s(1)%matrix,para_env,charges,error=error)
      END IF
      nkind = SIZE(atomic_kind_set)
      DO ikind=1,nkind
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,natom=natom_kind,scptb_parameter=scptb_kind)
         CALL get_scptb_parameter(scptb_kind,zeff=zeff)
         DO iatom=1,natom_kind
           atom_a = atomic_kind%atom_list(iatom)
           mcharge(atom_a) = zeff - SUM(charges(atom_a,1:nspins))
         END DO
      END DO
      DEALLOCATE(charges,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(dqpair(natom),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      CALL coulomb_pair(qs_env,dqpair,mcharge,calculate_forces,just_energy,error)
      IF(.NOT.just_energy) THEN
         CALL coulomb_pair_matrix(qs_env,ks_matrix,matrix_s,matrix_p,dqpair,dmcharge,calculate_forces,error)
      END IF
      CALL scp_coulomb(qs_env,ks_matrix,matrix_s,rho,mcharge,dmcharge,calculate_forces,just_energy,error)

      DEALLOCATE(mcharge,dqpair,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      IF(calculate_forces) THEN
         DEALLOCATE(dmcharge,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      energy%total = energy%core + energy%hartree + energy%qmmm_el + energy%repulsive + &
                     energy%core_self + energy%core_overlap + energy%dispersion

      output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%DETAILED_ENERGY",&
         extension=".scfLog",error=error)
      IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T9,A,T60,F20.10))")&
               "Repulsive pair potential energy:               ",energy%repulsive,&
               "Zeroth order Hamiltonian energy:               ",energy%core,&
               "Charge fluctuation (SCP) energy:               ",energy%hartree,&
               "Charge overlap correction energy:              ",energy%core_overlap,&
               "SCP self energy:                               ",energy%core_self,&
               "London dispersion energy:                      ",energy%dispersion
          IF (qs_env%qmmm) THEN
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "QM/MM Electrostatic energy:                    ",energy%qmmm_el
          END IF
      END IF
      CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
           "PRINT%DETAILED_ENERGY", error=error)
       ! here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C (plus occupation numbers)
      IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
         CALL get_qs_env(qs_env,mo_derivs=mo_derivs,mos=mo_array,error=error)
         DO ispin=1,SIZE(mo_derivs)
            CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                 mo_coeff_b=mo_coeff, occupation_numbers=occupation_numbers )
            IF(.NOT.mo_array(ispin)%mo_set%use_mo_coeff_b) THEN
               CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
            ENDIF
            CALL cp_dbcsr_multiply('n','n',1.0_dp,ks_matrix(ispin)%matrix,mo_coeff,&
                 0.0_dp,mo_derivs(ispin)%matrix, error=error)
         ENDDO
      ENDIF

    END IF

    CALL timestop(handle)

  END SUBROUTINE build_scptb_ks_matrix

! *****************************************************************************
  SUBROUTINE coulomb_pair(qs_env,dqpair,mcharge,calculate_forces,just_energy,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), DIMENSION(:), INTENT(out) :: dqpair
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: mcharge
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coulomb_pair', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, handle, &
                                                iatom, ikind, jatom, jkind, &
                                                natom, nkind, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: defined, failure, use_virial
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: scptb_defined
    REAL(KIND=dp)                            :: a, ab, b, dr1, dr2, eij, &
                                                eoverlap, fij
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sac_ae
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(scptb_parameter_p_type), &
      DIMENSION(:), POINTER                  :: scptb_kind_param
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind_a, scptb_kind_b
    TYPE(virial_type), POINTER               :: virial

    failure  = .FALSE.
    eoverlap = 0.0_dp

    CALL timeset(routineN,handle)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CALL get_qs_env(qs_env=qs_env,para_env=para_env,dft_control=dft_control,&
                    virial=virial,energy=energy,error=error)

    ! Parameters
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env=qs_env,sac_ae=sac_ae,atomic_kind_set=atomic_kind_set,error=error)

       nkind = SIZE(atomic_kind_set)
       IF(.NOT.just_energy) THEN
          dqpair = 0._dp
       END IF
       ! Possibly compute forces
       IF(calculate_forces) THEN
          CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,force=force,error=error)
          natom = SIZE (particle_set)
          ALLOCATE (atom_of_kind(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)
       END IF

       ALLOCATE (scptb_kind_param(nkind),scptb_defined(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,scptb_parameter=scptb_kind_a)
          scptb_kind_param(ikind)%scptb_param => scptb_kind_a
          CALL get_scptb_parameter(scptb_kind_a,defined=defined)
          scptb_defined(ikind)=defined
       END DO
       CALL neighbor_list_iterator_create(nl_iterator,sac_ae)
       DO WHILE (neighbor_list_iterate(nl_iterator)==0)
          CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,r=rij)
          IF (.NOT.scptb_defined(ikind)) CYCLE
          IF (.NOT.scptb_defined(jkind)) CYCLE
          scptb_kind_a => scptb_kind_param(ikind)%scptb_param
          scptb_kind_b => scptb_kind_param(jkind)%scptb_param
          dr2 = DOT_PRODUCT(rij,rij)
          IF ( dr2 > 0.00001_dp ) THEN
             ! charge overlap correction
             dr1 = SQRT(dr2)
             a = scptb_kind_a%a0
             b = scptb_kind_b%a0
             ab = SQRT(a*b/(a+b))
             eij = mcharge(iatom)*mcharge(jatom)/dr1 * erfc(ab*dr1)
             IF(.NOT.just_energy .AND. ABS(eij)>0.0_dp) THEN
               dqpair(iatom) = dqpair(iatom) + eij/mcharge(iatom)
               dqpair(jatom) = dqpair(jatom) + eij/mcharge(jatom)
             END IF
             IF(calculate_forces) THEN
                fij = eij/dr2 + 2._dp/SQRT(pi)*ab*mcharge(iatom)*mcharge(jatom)/dr2*EXP(-ab*dr2)
                force_ab = fij*rij

                atom_a = atom_of_kind(iatom)
                atom_b = atom_of_kind(jatom)

                ! Sum up force components
                force(ikind)%all_potential(1,atom_a) = force(ikind)%all_potential(1,atom_a) - force_ab(1)
                force(jkind)%all_potential(1,atom_b) = force(jkind)%all_potential(1,atom_b) + force_ab(1)

                force(ikind)%all_potential(2,atom_a) = force(ikind)%all_potential(2,atom_a) - force_ab(2)
                force(jkind)%all_potential(2,atom_b) = force(jkind)%all_potential(2,atom_b) + force_ab(2)

                force(ikind)%all_potential(3,atom_a) = force(ikind)%all_potential(3,atom_a) - force_ab(3)
                force(jkind)%all_potential(3,atom_b) = force(jkind)%all_potential(3,atom_b) + force_ab(3)
             END IF
          END IF
          eoverlap = eoverlap + eij
       END DO
       CALL neighbor_list_iterator_release(nl_iterator)

       DEALLOCATE (scptb_kind_param,scptb_defined,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       IF (calculate_forces) THEN
          DEALLOCATE(atom_of_kind,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       CALL mp_sum(eoverlap,para_env%group)
       energy%core_overlap = eoverlap

       IF(.NOT.just_energy) THEN
          CALL mp_sum(dqpair,para_env%group)
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE coulomb_pair

! *****************************************************************************
  SUBROUTINE coulomb_pair_matrix(qs_env, ks_matrix, matrix_s, matrix_p, dqpair,&
             dmcharge, calculate_forces, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_s, matrix_p
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: dqpair
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dmcharge
    LOGICAL, INTENT(IN)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coulomb_pair_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, atom_b, blk, i, &
                                                iatom, ikind, is, istat, &
                                                jatom, jkind, natom, nsize
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind, kind_of
    LOGICAL                                  :: failure, found, use_virial
    REAL(KIND=dp)                            :: dq, fi
    REAL(KIND=dp), DIMENSION(3)              :: fij, rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ksblock, sblock
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(virial_type), POINTER               :: virial

    nsize = SIZE(ks_matrix,1)
    IF ( calculate_forces ) THEN
      natom = SIZE(dmcharge,1)
      ALLOCATE (atom_of_kind(natom),kind_of(natom),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CALL get_qs_env(qs_env=qs_env,cell=cell,atomic_kind_set=atomic_kind_set,force=force,error=error)
      CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,virial=virial,error=error)
      use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               kind_of=kind_of,&
                               atom_of_kind=atom_of_kind)
    END IF

    CALL cp_dbcsr_iterator_start(iter,matrix_s(1)%matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
       dq = dqpair(iatom)+dqpair(jatom)
       DO is=1,nsize
          NULLIFY(ksblock)
          CALL cp_dbcsr_get_block_p(matrix=ks_matrix(is)%matrix,&
               row=iatom,col=jatom,block=ksblock,found=found)
          ksblock = ksblock + dq*sblock
       END DO
       IF ( calculate_forces ) THEN
            ikind  = kind_of(iatom)
            atom_a = atom_of_kind(iatom)
            jkind  = kind_of(jatom)
            atom_b = atom_of_kind(jatom)
            DO i=1,3
               fi = dqpair(jatom)*dmcharge(iatom,i) - dqpair(iatom)*dmcharge(jatom,i)
               force(ikind)%rho_elec(i,atom_a) = force(ikind)%rho_elec(i,atom_a) + fi
               force(jkind)%rho_elec(i,atom_b) = force(jkind)%rho_elec(i,atom_b) - fi
               fij(i) = fi
            END DO
            IF (use_virial) THEN
               rij = particle_set(iatom)%r - particle_set(jatom)%r
               rij = pbc(rij,cell)
               CALL virial_pair_force ( virial%pv_virial, -1._dp, fij, rij, error)
            END IF
       END IF
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    IF ( calculate_forces ) THEN
      DEALLOCATE (atom_of_kind,kind_of,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE coulomb_pair_matrix

! *****************************************************************************
  SUBROUTINE scp_coulomb(qs_env,ks_matrix,matrix_s,rho,mcharge,dmcharge,&
             calculate_forces,just_energy,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_s
    TYPE(qs_rho_type), POINTER               :: rho
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: mcharge
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: dmcharge
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_coulomb', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, iat, iatom, ibas, ikind, is, istat, jat, jatom, jkind, &
      kint, kmax, lmaxscp, nat, natom, nkind, norb, nsize
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, lscp, &
                                                natoms, nbasis
    LOGICAL                                  :: defined, failure, found
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: adef
    REAL(KIND=dp)                            :: ak, alpha, bk, chabs, dq, &
                                                epsabs, epsrel, ff, rhom1, &
                                                rhom2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: zeta
    REAL(KIND=dp), DIMENSION(0:3)            :: pol
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ksblock, sblock
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(scp_vector_type), POINTER           :: cp, cres, crhs, cvec, cw
    TYPE(scptb_control_type), POINTER        :: scptb_control
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
         energy=energy,para_env=para_env,error=error)

    energy%hartree = 0._dp
    energy%core_self = 0._dp

    chabs = SUM(ABS(mcharge))
    IF ( chabs > 1.e-8_dp ) THEN

       scptb_control => dft_control%qs_control%scptb_control

       NULLIFY(cvec,crhs,cp,cw,cres)

       epsrel = 1.e-8_dp

       ! Create scp vectors
       CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,error=error)
       nkind = SIZE(atomic_kind_set)
       ALLOCATE(nbasis(nkind),natoms(nkind),lscp(nkind),zeta(nkind),adef(nkind),stat=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,natom=nat,scptb_parameter=scptb_kind)
          CALL get_scptb_parameter(scptb_kind,defined=defined,lmaxscp=lmaxscp,ag=alpha)
          IF(defined) THEN
            norb = (lmaxscp+1)**2
          ELSE
            norb = 0
          END IF
          natoms(ikind) = nat
          nbasis(ikind) = norb
          lscp(ikind)   = lmaxscp
          zeta(ikind)   = alpha
          adef(ikind)   = defined
       END DO
       CALL scp_vector_create(cvec, nkind, natoms, nbasis, error)
       CALL scp_vector_set(cvec, 0._dp, error)
       ! set charge
       CALL scp_set_charge(cvec,mcharge,adef,natoms,atomic_kind_set)
       ! Calclulate scp vector
       IF (scptb_control%do_ewald) THEN
          kmax = MIN(SUM(nbasis),100) * 10
          CALL scp_vector_create(crhs, nkind, natoms, nbasis, error)
          CALL scp_vector_set(crhs, 0._dp, error)
          ! for ci=0, this creates the RHS of (Ax=b)
          CALL apply_scp_vector(cvec,crhs,qs_env,calculate_forces=.FALSE.,error=error)
          CALL scp_vector_scale(-1._dp, crhs, error)
          ! conjugate gradient optimization
          CALL scp_vector_create(cres, nkind, natoms, nbasis, error)
          CALL scp_vector_copy(crhs, cres, error)
          CALL scp_vector_create(cp, nkind, natoms, nbasis, error)
          CALL scp_vector_create(cw, nkind, natoms, nbasis, error)
          ! scp_dot is not using the c0 (charge) term
          rhom1 = scp_dot(cres,cres)
          epsabs = epsrel*SQRT(rhom1)
          DO kint=1,kmax
             IF(kint==1) THEN
                CALL scp_vector_copy(cres, cp, error)
             ELSE
                bk = rhom1/rhom2
                CALL scp_vector_scale(bk, cp, error)
                CALL scp_vector_add(1._dp, cres, cp, error)
             END IF
             ! charge term
             CALL scp_set_charge(cp,mcharge,adef,natoms,atomic_kind_set)
             ! Ax - b
             CALL scp_vector_set(cw, 0._dp, error)
             CALL apply_scp_vector(cp,cw,qs_env,calculate_forces=.FALSE.,error=error)
             ! add b to get Ax
             CALL scp_vector_add(1._dp, crhs, cw, error)
             ! add force from self term ci/ai
             DO ikind=1,nkind
               IF(adef(ikind)) THEN
                 atomic_kind => atomic_kind_set(ikind)
                 CALL get_atomic_kind(atomic_kind=atomic_kind,natom=nat,scptb_parameter=scptb_kind)
                 CALL get_scptb_parameter(scptb_kind,pol=pol)
                 DO iat=1,natoms(ikind)
                   DO ibas=1,nbasis(ikind)
                      cw%vector(ikind)%vmat(ibas,iat) = cw%vector(ikind)%vmat(ibas,iat) &
                         + cp%vector(ikind)%vmat(ibas,iat)/pol(indso(1,ibas))
                   END DO
                 END DO
               END IF
             END DO
             !
             ak = rhom1/scp_dot(cp,cw)
             CALL scp_vector_add(ak, cp, cvec, error)
             CALL scp_vector_add(-ak, cw, cres, error)
             rhom2 = rhom1
             rhom1 = scp_dot(cres,cres)
             IF ( SQRT(rhom1) < epsabs ) EXIT
          END DO
!deb
          WRITE(6,*) "SCP OPT   ",kint,SQRT(rhom1)
!deb
          ! Coulomb energy
          CALL scp_set_charge(cvec,mcharge,adef,natoms,atomic_kind_set)
          CALL scp_vector_set(cw, 0._dp, error)
          CALL apply_scp_vector(cvec,cw,qs_env,calculate_forces=.FALSE.,error=error)
          CALL scp_vector_dot(ak, cw, cvec, error)
          energy%hartree = 0.5_dp*ak
          ! release scratch vectors
          CALL scp_vector_release(crhs, error)
          CALL scp_vector_release(cres, error)
          CALL scp_vector_release(cp, error)
          CALL scp_vector_release(cw, error)
       ELSE
          !debug: do nothing for now
          ! Coulomb energy
          energy%hartree = 0._dp
       END IF
       ! Calculate self energy
       DO ikind=1,nkind
          IF(adef(ikind)) THEN
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind=atomic_kind,natom=nat,scptb_parameter=scptb_kind)
            CALL get_scptb_parameter(scptb_kind,pol=pol)
            DO iat=1,natoms(ikind)
               DO ibas=1,nbasis(ikind)
                  energy%core_self = energy%core_self + &
                     0.5_dp*(cvec%vector(ikind)%vmat(ibas,iat))**2/pol(indso(1,ibas))
               END DO
            END DO
          END IF
       END DO

       IF (.NOT.just_energy) THEN

          natom = SIZE(mcharge,1)
          ALLOCATE (atom_of_kind(natom),kind_of(natom),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)

          ! calculate KS matrix
          IF (scptb_control%do_ewald) THEN
             CALL scp_vector_create(cw, nkind, natoms, nbasis, error)
             CALL scp_set_charge(cvec,mcharge,adef,natoms,atomic_kind_set)
             CALL scp_vector_set(cw, 0._dp, error)
             CALL apply_scp_vector(cvec,cw,qs_env,calculate_forces=.FALSE.,error=error)

             nsize = SIZE(ks_matrix)
             CALL cp_dbcsr_iterator_start(iter,matrix_s(1)%matrix)
             DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
                CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sblock, blk)
                ikind = kind_of(iatom)
                iat   = atom_of_kind(iatom)
                jkind = kind_of(jatom)
                jat   = atom_of_kind(jatom)
                IF (iatom==jatom) THEN
                  ff = 0.5_dp
                ELSE
                  ff = 1.0_dp
                END IF
                dq = (cw%vector(ikind)%vmat(1,iat) + cw%vector(jkind)%vmat(1,jat))*ff
                DO is=1,nsize
                   NULLIFY(ksblock)
                   CALL cp_dbcsr_get_block_p(matrix=ks_matrix(is)%matrix,&
                        row=iatom,col=jatom,block=ksblock,found=found)
                   ksblock = ksblock + dq*sblock
                END DO
             END DO
             CALL cp_dbcsr_iterator_stop(iter)

             CALL scp_vector_release(cw, error)
          ELSE
             ! explicit integrals
          END IF
          IF (calculate_forces) THEN
             ! forces on atoms
          END IF
       END IF

       CALL scp_vector_release(cvec, error)
       DEALLOCATE(nbasis,natoms,lscp,zeta,adef,stat=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    END IF

  END SUBROUTINE scp_coulomb

  FUNCTION scp_dot(vec1,vec2) RESULT(res)
    TYPE(scp_vector_type)                    :: vec1, vec2
    REAL(KIND=dp)                            :: res

    INTEGER                                  :: i

    res = 0._dp
    DO i=1,SIZE(vec1%vector)
      res = res + SUM(vec1%vector(i)%vmat(2:,:)*vec2%vector(i)%vmat(2:,:))
    END DO

  END FUNCTION scp_dot

  SUBROUTINE scp_set_charge(vec,mcharge,adef,natoms,atomic_kind_set)
    TYPE(scp_vector_type)                    :: vec
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: mcharge
    LOGICAL, DIMENSION(:), INTENT(IN)        :: adef
    INTEGER, DIMENSION(:), INTENT(IN)        :: natoms
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

    INTEGER                                  :: iat, iatom, ikind, nkind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

! set charge

    nkind = SIZE(adef)
    DO ikind=1,nkind
       IF(adef(ikind)) THEN
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,atom_list=atom_list)
         DO iat=1,natoms(ikind)
            iatom = atom_list(iat)
            vec%vector(ikind)%vmat(1,iat) = mcharge(iatom)
         END DO
       END IF
    END DO

  END SUBROUTINE scp_set_charge

! *****************************************************************************

  SUBROUTINE apply_scp_vector(cin,cout,qs_env,calculate_forces,error)

    TYPE(scp_vector_type), POINTER           :: cin, cout
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_scp_vector', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: use_virial
    REAL(KIND=dp)                            :: ehartree
    REAL(KIND=dp), DIMENSION(3, 3)           :: h_stress
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: scp_pot, scp_rho, scp_rho_g
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(virial_type), POINTER               :: virial

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,virial=virial,para_env=para_env,error=error)
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,&
                    pw_pools=pw_pools,poisson_env=poisson_env,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                           scp_rho%pw,&
                           use_data=REALDATA3D,&
                           in_space=REALSPACE,&
                           error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                           scp_pot%pw, &
                           use_data=COMPLEXDATA1D,&
                           in_space=RECIPROCALSPACE,&
                           error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                           scp_rho_g%pw, &
                           use_data=COMPLEXDATA1D,&
                           in_space=RECIPROCALSPACE,&
                           error=error)

    CALL calculate_scp_charge(scp_rho,qs_env,cin,error)
    CALL pw_transfer(scp_rho%pw, scp_rho_g%pw, error=error)

    ! Getting the Hartree energy and Hartree potential.  Also getting the stress tensor
    ! from the Hartree term if needed. 
    IF (use_virial .AND. calculate_forces) THEN
       h_stress(:,:) = 0.0_dp
       CALL pw_poisson_solve(poisson_env,scp_rho_g%pw,ehartree,scp_pot%pw,h_stress=h_stress,error=error)
       virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe,dp)
    ELSE
       CALL pw_poisson_solve(poisson_env,scp_rho_g%pw,ehartree,scp_pot%pw,error=error)
    END IF
    CALL pw_pool_give_back_pw(auxbas_pw_pool,scp_rho_g%pw,error=error)

    CALL pw_transfer(scp_pot%pw, scp_rho%pw, error=error)
    CALL integrate_scp_rspace(scp_rho,qs_env,cout,calculate_forces,error)

    CALL pw_pool_give_back_pw(auxbas_pw_pool,scp_rho%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,scp_pot%pw,error=error)

  END SUBROUTINE apply_scp_vector

! *****************************************************************************

END MODULE scptb_ks_matrix

