!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Construction of the Exchange part of the Fock Matrix
!> \author Teodoro Laino [tlaino] (05.2009) - Split and module reorganization
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : d-orbitals
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Speed-up
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Periodic SE
! *****************************************************************************
MODULE se_fock_matrix_exchange
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             semi_empirical_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: do_multipole_none,&
                                             do_se_IS_kdso,&
                                             do_se_IS_kdso_d
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE se_fock_matrix_integrals,        ONLY: dfock2E,&
                                             fock1_2el,&
                                             fock2E
  USE semi_empirical_int_arrays,       ONLY: rij_threshold
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_type
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             se_int_control_type,&
                                             se_taper_type,&
                                             semi_empirical_p_type,&
                                             semi_empirical_type,&
                                             setup_se_int_control_type
  USE semi_empirical_utils,            ONLY: finalize_se_taper,&
                                             initialize_se_taper
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_fock_matrix_exchange'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module       = .FALSE.

  PUBLIC :: build_fock_matrix_exchange

CONTAINS


! *****************************************************************************
!> \brief Construction of the Exchange part of the Fock matrix
!> \param qs_env ...
!> \param ks_matrix ...
!> \param matrix_p ...
!> \param calculate_forces ...
!> \param store_int_env ...
!> \param error ...
!> \author JGH
! *****************************************************************************
  SUBROUTINE build_fock_matrix_exchange (qs_env, ks_matrix, matrix_p, calculate_forces,&
       store_int_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, matrix_p
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_fock_matrix_exchange', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iatom, icol, ikind, &
      integral_screening, irow, jatom, jkind, natom, natorb_a, nkind, nspins, &
      stat
    INTEGER, DIMENSION(2)                    :: size_p_block_a
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind
    LOGICAL                                  :: anag, check, defined, &
                                                failure, found, switch, &
                                                use_virial
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: se_defined
    REAL(KIND=dp)                            :: delta, dr
    REAL(KIND=dp), DIMENSION(3)              :: force_ab, rij
    REAL(KIND=dp), DIMENSION(45, 45)         :: p_block_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ks_block_a, ks_block_b, &
                                                p_block_a, p_block_b
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(se_int_control_type)                :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_p_type), &
      DIMENSION(:), POINTER                  :: se_kind_list
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b
    TYPE(virial_type), POINTER               :: virial

    failure=.FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(dft_control,cell,force,particle_set,se_control,se_taper)
    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,cell=cell,se_taper=se_taper,&
         para_env=para_env,virial=virial,error=error)

    CALL initialize_se_taper(se_taper,exchange=.TRUE.,error=error)
    se_control => dft_control%qs_control%se_control
    anag       =  se_control%analytical_gradients
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    nspins=dft_control%nspins

    CPPrecondition(ASSOCIATED(matrix_p),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN
       ! Identify proper integral screening (according user requests)
       integral_screening = se_control%integral_screening
       IF ((integral_screening==do_se_IS_kdso_d).AND.(.NOT.se_control%force_kdsod_EX)) THEN
          integral_screening = do_se_IS_kdso
       END IF
       CALL setup_se_int_control_type(se_int_control, shortrange=.FALSE.,&
            do_ewald_r3=.FALSE., do_ewald_gks=.FALSE., integral_screening=integral_screening,&
            max_multipole=do_multipole_none, pc_coulomb_int=.FALSE.)

       CALL get_qs_env(qs_env=qs_env,sab_orb=sab_orb,&
             atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,error=error)

       nkind = SIZE(atomic_kind_set)
       IF(calculate_forces) THEN
          CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,force=force,error=error)
          natom = SIZE (particle_set)
          ALLOCATE (atom_of_kind(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          delta = se_control%delta
          CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)
       END IF

     ALLOCATE (se_defined(nkind),se_kind_list(nkind),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO ikind=1,nkind
        CALL get_qs_kind(qs_kind_set(ikind),se_parameter=se_kind_a)
        se_kind_list(ikind)%se_param => se_kind_a
        CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a)
        se_defined(ikind) = (defined .AND. natorb_a >= 1)
     END DO

     CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
     DO WHILE (neighbor_list_iterate(nl_iterator)==0)
        CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,r=rij)
        IF (.NOT.se_defined(ikind)) CYCLE
        IF (.NOT.se_defined(jkind)) CYCLE
        se_kind_a => se_kind_list(ikind)%se_param
        se_kind_b => se_kind_list(jkind)%se_param

        IF (iatom <= jatom) THEN
           irow = iatom
           icol = jatom
           switch = .FALSE.
        ELSE
           irow = jatom
           icol = iatom
           switch = .TRUE.
        END IF
        ! Retrieve blocks for KS and P
        CALL cp_dbcsr_get_block_p(matrix=ks_matrix(1)%matrix,&
             row=irow,col=icol,BLOCK=ks_block_a,found=found)
        CPPostcondition(ASSOCIATED(ks_block_a),cp_failure_level,routineP,error,failure)
        CALL cp_dbcsr_get_block_p(matrix=matrix_p(1)%matrix,&
             row=irow,col=icol,BLOCK=p_block_a,found=found)
        CPPostcondition(ASSOCIATED(p_block_a),cp_failure_level,routineP,error,failure)
        size_p_block_a(1) = SIZE(p_block_a,1)
        size_p_block_a(2) = SIZE(p_block_a,2)
        p_block_tot(1:size_p_block_a(1),1:size_p_block_a(2)) = 2.0_dp * p_block_a

        ! Handle more configurations
        IF ( nspins == 2 ) THEN
           CALL cp_dbcsr_get_block_p(matrix=ks_matrix(2)%matrix,&
                row=irow,col=icol,BLOCK=ks_block_b,found=found)
           CPPostcondition(ASSOCIATED(ks_block_b),cp_failure_level,routineP,error,failure)
           CALL cp_dbcsr_get_block_p(matrix=matrix_p(2)%matrix,&
                row=irow,col=icol,BLOCK=p_block_b,found=found)
           CPPostcondition(ASSOCIATED(p_block_b),cp_failure_level,routineP,error,failure)
           check = (size_p_block_a(1)==SIZE(p_block_b,1)).AND.(size_p_block_a(2)==SIZE(p_block_b,2))
           CPPostcondition(check,cp_failure_level,routineP,error,failure)
           p_block_tot(1:SIZE(p_block_a,1),1:SIZE(p_block_a,2)) = p_block_a + p_block_b
        END IF

        dr = DOT_PRODUCT(rij,rij)
        IF ( iatom == jatom .AND. dr < rij_threshold ) THEN
           ! Once center - Two electron Terms
           IF      ( nspins == 1 ) THEN
              CALL fock1_2el(se_kind_a,p_block_tot,p_block_a,ks_block_a,factor=0.5_dp,error=error)
           ELSE IF ( nspins == 2 ) THEN
              CALL fock1_2el(se_kind_a,p_block_tot,p_block_a,ks_block_a,factor=1.0_dp,error=error)
              CALL fock1_2el(se_kind_a,p_block_tot,p_block_b,ks_block_b,factor=1.0_dp,error=error)
           END IF
        ELSE
           ! Exchange Terms
           IF      ( nspins == 1 ) THEN
              CALL fock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_a, ks_block_a,&
                   factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, &
                   store_int_env=store_int_env, error=error)
           ELSE IF ( nspins == 2 ) THEN
              CALL fock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_a, ks_block_a,&
                   factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, &
                   store_int_env=store_int_env, error=error)

              CALL fock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_b, ks_block_b,&
                   factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, &
                   store_int_env=store_int_env, error=error)
           END IF
           IF(calculate_forces) THEN
              atom_a = atom_of_kind(iatom)
              atom_b = atom_of_kind(jatom)
              force_ab = 0.0_dp
              IF      ( nspins == 1 ) THEN
                 CALL dfock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_a,&
                      factor=0.5_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab,&
                      delta=delta, error=error)
              ELSE IF ( nspins == 2 ) THEN
                 CALL dfock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_a,&
                      factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab,&
                      delta=delta, error=error)

                 CALL dfock2E(se_kind_a, se_kind_b, rij, switch, size_p_block_a, p_block_tot, p_block_b,&
                      factor=1.0_dp, anag=anag, se_int_control=se_int_control, se_taper=se_taper, force=force_ab,&
                      delta=delta, error=error)
              END IF
              IF (switch) THEN
                 force_ab(1) = -force_ab(1)
                 force_ab(2) = -force_ab(2)
                 force_ab(3) = -force_ab(3)
              END IF
              IF (use_virial) THEN
                 CALL virial_pair_force ( virial%pv_virial, -1.0_dp, force_ab, rij, error)
              END IF

              force(ikind)%rho_elec(1,atom_a) = force(ikind)%rho_elec(1,atom_a) - force_ab(1)
              force(jkind)%rho_elec(1,atom_b) = force(jkind)%rho_elec(1,atom_b) + force_ab(1)

              force(ikind)%rho_elec(2,atom_a) = force(ikind)%rho_elec(2,atom_a) - force_ab(2)
              force(jkind)%rho_elec(2,atom_b) = force(jkind)%rho_elec(2,atom_b) + force_ab(2)

              force(ikind)%rho_elec(3,atom_a) = force(ikind)%rho_elec(3,atom_a) - force_ab(3)
              force(jkind)%rho_elec(3,atom_b) = force(jkind)%rho_elec(3,atom_b) + force_ab(3)
           END IF
        END IF
     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     DEALLOCATE(se_kind_list,se_defined,stat=stat)
     CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

     IF (calculate_forces) THEN
        DEALLOCATE(atom_of_kind,stat=stat)
        CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
    END IF
    CALL finalize_se_taper(se_taper,error=error)

    CALL timestop(handle)

  END SUBROUTINE build_fock_matrix_exchange

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE build_fock_matrix_ph(qs_env, ks_matrix, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_fock_matrix_ph', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, icol, &
                                                ikind, irow, jatom, jkind, k, &
                                                natorb_a, nkind, nspins, stat
    LOGICAL                                  :: anag, defined, failure, &
                                                found, ifHi, ifHj, switch
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: se_defined
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ks_block_a, ks_block_b
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_p_type), &
      DIMENSION(:), POINTER                  :: se_kind_list
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b

    failure=.FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(dft_control, se_control)
    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
         para_env=para_env,error=error)

    se_control => dft_control%qs_control%se_control
    anag       =  se_control%analytical_gradients
    nspins=dft_control%nspins

    CPPrecondition(SIZE(ks_matrix)>0,cp_failure_level,routineP,error,failure)

    IF ( .NOT. failure ) THEN

       CALL get_qs_env(qs_env=qs_env,sab_orb=sab_orb,&
         atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,error=error)

       nkind = SIZE(atomic_kind_set)

     ALLOCATE (se_defined(nkind),se_kind_list(nkind),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO ikind=1,nkind
        CALL get_qs_kind(qs_kind_set(ikind), se_parameter=se_kind_a)
        se_kind_list(ikind)%se_param => se_kind_a
        CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a)
        se_defined(ikind) = (defined .AND. natorb_a >= 1)
     END DO

     CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
     DO WHILE (neighbor_list_iterate(nl_iterator)==0)
        CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom)
        IF (.NOT.se_defined(ikind)) CYCLE
        IF (.NOT.se_defined(jkind)) CYCLE
        se_kind_a => se_kind_list(ikind)%se_param
        se_kind_b => se_kind_list(jkind)%se_param

        IF (iatom <= jatom) THEN
           irow = iatom
           icol = jatom
           switch = .FALSE.
        ELSE
           irow = jatom
           icol = iatom
           switch = .TRUE.
        END IF
        ! Retrieve blocks for KS
        CALL cp_dbcsr_get_block_p(matrix=ks_matrix(1)%matrix,&
             row=irow,col=icol,BLOCK=ks_block_a,found=found)
        CPPostcondition(ASSOCIATED(ks_block_a),cp_failure_level,routineP,error,failure)

                   ifHi= (atomic_kind_set(ikind)%element_symbol=='H ')
                   ifHj= (atomic_kind_set(jkind)%element_symbol=='H ')

                   IF( (iatom == jatom).AND. ifHi )THEN
                     DO k=2,4
                       DO i=1,k-1
                         ks_block_a(i,k)= 0.0_dp
                         ks_block_a(k,i)= 0.0_dp
                       ENDDO
                     ENDDO
                   ENDIF

                   IF( ifHi )THEN
                     IF( irow==iatom )THEN
                       DO k=2,4
                         DO i=1,4
                           ks_block_a(k,i)= 0.0_dp
                         ENDDO
                       ENDDO
                     ELSE
                       DO k=2,4
                         DO i=1,4
                           ks_block_a(i,k)= 0.0_dp
                         ENDDO
                       ENDDO
                     ENDIF
                   ENDIF

                   IF( ifHj )THEN
                     IF( irow==iatom )THEN
                       DO k=2,4
                         DO i=1,4
                           ks_block_a(i,k)= 0.0_dp
                         ENDDO
                       ENDDO
                     ELSE
                       DO k=2,4
                         DO i=1,4
                           ks_block_a(k,i)= 0.0_dp
                         ENDDO
                       ENDDO
                     ENDIF
                   ENDIF

        ! Handle more configurations
        IF ( nspins == 2 ) THEN
           CALL cp_dbcsr_get_block_p(matrix=ks_matrix(2)%matrix,&
                row=irow,col=icol,BLOCK=ks_block_b,found=found)
           CPPostcondition(ASSOCIATED(ks_block_b),cp_failure_level,routineP,error,failure)
        END IF

     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     DEALLOCATE(se_kind_list,se_defined,stat=stat)
     CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    END IF

    CALL timestop(handle)

  END SUBROUTINE build_fock_matrix_ph

END MODULE se_fock_matrix_exchange

