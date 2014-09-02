!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief Perform the calculation of the dispersion contributions
!>        for the SCP method
!> \author CJM
!******************************************************************************
MODULE  scp_dispersion
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             semi_empirical_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: do_se_IS_slater
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
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
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE semi_empirical_int_arrays,       ONLY: se_orbital_pointer
  USE semi_empirical_integrals,        ONLY: drotint,&
                                             rotint
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
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Global parameters (only in this module)
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_dispersion'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

  ! Public Subroutine
  PUBLIC :: scp_nddo_dispersion


CONTAINS

! *****************************************************************************
!> \brief Analytic calculation of the dispersion energy and force
!>        for SCP using the QS environment
!> \param qs_env ...
!> \param calculate_forces ...
!> \param error ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE scp_nddo_dispersion( qs_env, calculate_forces, error )
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_nddo_dispersion', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_atom, iatom, ikind, &
                                                inode, j_atom, jatom, jkind, &
                                                natom, natorb_a, nkind, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    LOGICAL                                  :: anag, defined, failure, &
                                                forces, old_store_int, switch
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: se_defined
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: se_de
    REAL(kind=dp)                            :: dab, de_a, de_b, disp_factor, &
                                                e_scp, force_ab( 3 ), &
                                                rab( 3 ), rab2
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_se
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(se_int_control_type)                :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_p_type), &
      DIMENSION(:), POINTER                  :: se_kind_list
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(semi_empirical_type), POINTER       :: se_kind_a, se_kind_b

    CALL timeset(routineN,handle)

    failure = .FALSE.
    forces = .FALSE.
    IF ( PRESENT(calculate_forces) ) forces = calculate_forces

    e_scp = 0._dp

    NULLIFY( scp_env, se_taper, store_int_env)
    NULLIFY( energy, force, para_env, atomic_kind_set, qs_kind_set )
    NULLIFY( cell, sab_se, dft_control, se_control )

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
         qs_kind_set=qs_kind_set, cell=cell, scp_env=scp_env, se_taper=se_taper, &
         se_store_int_env=store_int_env, &
         para_env=para_env, particle_set=particle_set,&
         sab_se=sab_se, force=force, dft_control=dft_control, &
         error=error)

    CALL initialize_se_taper(se_taper,coulomb=.TRUE.,error=error)

    ! Parameters
    se_control => dft_control%qs_control%se_control
    anag       =  se_control%analytical_gradients
    CALL setup_se_int_control_type(se_int_control,do_ewald_r3=.FALSE.,&
         do_ewald_gks=.FALSE., integral_screening=do_se_IS_slater,&
         shortrange=.FALSE., max_multipole=se_control%max_multipole,&
         pc_coulomb_int=.FALSE.)

    CALL get_scp_env(scp_env=scp_env, energy=energy, error=error)

    nkind = SIZE(atomic_kind_set,1)
    natom = SIZE(particle_set,1)

    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

    ALLOCATE (se_defined(nkind),se_kind_list(nkind),se_de(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
       CALL get_qs_kind(qs_kind_set(ikind), scp_se_parameter=se_kind_a)
       se_kind_list(ikind)%se_param => se_kind_a
       CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a,de=de_a)
       se_defined(ikind) = (defined .AND. natorb_a >= 1)
       se_de(ikind) = de_a
    END DO

    CALL neighbor_list_iterator_create(nl_iterator,sab_se)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=iatom,jatom=jatom,inode=inode,r=rab)
       IF (.NOT.se_defined(ikind)) CYCLE
       IF (.NOT.se_defined(jkind)) CYCLE
       se_kind_a => se_kind_list(ikind)%se_param
       se_kind_b => se_kind_list(jkind)%se_param
       de_a = se_de(ikind)
       de_b = se_de(jkind)
       ! Computing the constants for the dispersion interaction
       IF ( (de_a + de_b) == 0._dp ) THEN
          disp_factor = 0._dp
       ELSE
          disp_factor = -0.25_dp * de_a * de_b / ( de_a + de_b )
       END IF
       ! this assumes MIC, but we have dispersion that would be long range?????
       IF ( iatom == jatom ) CYCLE

       rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       dab = SQRT ( rab2 )

       IF (iatom <= jatom) THEN
          switch = .FALSE.
       ELSE
          switch = .TRUE.
       END IF

       force_ab = 0.0_dp

       old_store_int = store_int_env%memory_parameter%do_all_on_the_fly
       store_int_env%memory_parameter%do_all_on_the_fly= .TRUE.

       CALL scp_disp_2C ( se_kind_a, se_kind_b, rab, switch, factor=disp_factor, &
            se_taper=se_taper, store_int_env=store_int_env,      &
            se_int_control=se_int_control, edisp=e_scp, error=error)
       IF ( forces ) THEN
          CALL dscp_disp_2C ( se_kind_a, se_kind_b, rab, switch, factor=disp_factor, &
               se_taper=se_taper, store_int_env=store_int_env,  &
               se_int_control=se_int_control, force=force_ab, error=error)
          i_atom = atom_of_kind(iatom)
          j_atom = atom_of_kind(jatom)
          ! Sum up force components
          force(ikind)%all_potential(1,i_atom) = force(ikind)%all_potential(1,i_atom) - force_ab(1)
          force(jkind)%all_potential(1,j_atom) = force(jkind)%all_potential(1,j_atom) + force_ab(1)

          force(ikind)%all_potential(2,i_atom) = force(ikind)%all_potential(2,i_atom) - force_ab(2)
          force(jkind)%all_potential(2,j_atom) = force(jkind)%all_potential(2,j_atom) + force_ab(2)

          force(ikind)%all_potential(3,i_atom) = force(ikind)%all_potential(3,i_atom) - force_ab(3)
          force(jkind)%all_potential(3,j_atom) = force(jkind)%all_potential(3,j_atom) + force_ab(3)
       END IF

       store_int_env%memory_parameter%do_all_on_the_fly = old_store_int
    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    CALL mp_sum(e_scp,para_env%group)
    energy % e_scp_disp =  e_scp

    DEALLOCATE(atom_of_kind,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (se_defined,se_kind_list,se_de,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL finalize_se_taper(se_taper,error=error)
    CALL timestop(handle)

  END SUBROUTINE scp_nddo_dispersion

! *****************************************************************************
!> \brief Construction of 2-center Coulomb Terms for SCP-NDDO dispersion
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param switch ...
!> \param factor ...
!> \param se_taper ...
!> \param store_int_env ...
!> \param se_int_control ...
!> \param edisp ...
!> \param error ...
!> \date 04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE scp_disp_2C(sepi, sepj, rij, switch, factor, se_taper, store_int_env,  &
                         se_int_control, edisp, error)

    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), INTENT(IN)                :: factor
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    REAL(KIND=dp), INTENT(INOUT)             :: edisp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_disp_2C', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aa, app_a, app_b, &
                                                ascp_a( 9, 9 ), &
                                                ascp_b( 9, 9 ), asp_a, asp_b, &
                                                ass_a, ass_b, bb, irij(3)
    REAL(KIND=dp), DIMENSION(2025)           :: w

    failure = .FALSE.
    CALL get_se_param(sepi, ass = ass_a, asp = asp_a, app = app_a )
    ascp_a( 1  ,  1) = ass_a
    ascp_a( 1  ,2:4) = asp_a
    ascp_a( 2:4,  1) = asp_a
    ascp_a( 2:4,2:4) = app_a
    CALL get_se_param(sepj, ass = ass_b, asp = asp_b, app = app_b )
    ascp_b( 1  ,  1) = ass_b
    ascp_b( 1  ,2:4) = asp_b
    ascp_b( 2:4,  1) = asp_b
    ascp_b( 2:4,2:4) = app_b
    ! Evaluate integrals
    IF (.NOT.switch) THEN
       CALL rotint (sepi,sepj, rij,w,anag=.TRUE.,se_int_control=se_int_control,se_taper=se_taper,&
            store_int_env=store_int_env, error=error)
    ELSE
       irij = -rij
       CALL rotint (sepj,sepi,irij,w,anag=.TRUE.,se_int_control=se_int_control,se_taper=se_taper,&
            store_int_env=store_int_env, error=error)
    END IF

    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                IF ( iL == 1 .AND.  jL == 1 ) CYCLE
                IF ( kL == 1 .AND.  lL == 1 ) CYCLE
                IF (switch) THEN
                  edisp = edisp + aa*bb*w(kr)*w(kr)*ascp_b(iL,jL)*ascp_a(kL,lL)*factor
                ELSE
                  edisp = edisp + aa*bb*w(kr)*w(kr)*ascp_a(iL,jL)*ascp_b(kL,lL)*factor
                END IF
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE scp_disp_2C

! *****************************************************************************
!> \brief Derivatives of 2-center Coulomb Terms for SCP-NDDO dispersion
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param switch ...
!> \param factor ...
!> \param se_taper ...
!> \param store_int_env ...
!> \param se_int_control ...
!> \param force ...
!> \param error ...
!> \date 04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE dscp_disp_2C(sepi, sepj, rij, switch, factor, &
                          se_taper, store_int_env, se_int_control, force, error)

    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), INTENT(IN)                :: factor
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dscp_disp_2C', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aa, app_a, app_b, &
                                                ascp_a( 9, 9 ), &
                                                ascp_b( 9, 9 ), asp_a, asp_b, &
                                                ass_a, ass_b, bb
    REAL(KIND=dp), DIMENSION(2025)           :: w
    REAL(KIND=dp), DIMENSION(3)              :: a, irij
    REAL(KIND=dp), DIMENSION(3, 2025)        :: dw

    failure = .FALSE.
    CALL get_se_param(sepi, ass = ass_a, asp = asp_a, app = app_a )
    ascp_a( 1  ,  1) = ass_a
    ascp_a( 1  ,2:4) = asp_a
    ascp_a( 2:4,  1) = asp_a
    ascp_a( 2:4,2:4) = app_a
    CALL get_se_param(sepj, ass = ass_b, asp = asp_b, app = app_b )
    ascp_b( 1  ,  1) = ass_b
    ascp_b( 1  ,2:4) = asp_b
    ascp_b( 2:4,  1) = asp_b
    ascp_b( 2:4,2:4) = app_b
    ! Evaluate integrals
    IF (.NOT.switch) THEN
       CALL rotint (sepi,sepj, rij,w,anag=.TRUE.,se_int_control=se_int_control,se_taper=se_taper,&
            store_int_env=store_int_env, error=error)
    ELSE
       irij = -rij
       CALL rotint (sepj,sepi,irij,w,anag=.TRUE.,se_int_control=se_int_control,se_taper=se_taper,&
            store_int_env=store_int_env, error=error)
    END IF
    ! Evaluate integrals' derivatives
    IF (.NOT.switch) THEN
       CALL drotint (sepi,sepj, rij,dw,delta=0.0_dp,anag=.TRUE.,se_int_control=se_int_control, &
                     se_taper=se_taper, error=error)
    ELSE
       irij = -rij
       CALL drotint (sepj,sepi,irij,dw,delta=0.0_dp,anag=.TRUE.,se_int_control=se_int_control, &
                     se_taper=se_taper, error=error)
    END IF

    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                IF ( iL == 1 .AND.  jL == 1 ) CYCLE
                IF ( kL == 1 .AND.  lL == 1 ) CYCLE
                a(1) = 2.0_dp*dw(1,kr)*w(kr)*factor*aa*bb
                a(2) = 2.0_dp*dw(2,kr)*w(kr)*factor*aa*bb
                a(3) = 2.0_dp*dw(3,kr)*w(kr)*factor*aa*bb
                IF (switch) THEN
                  force(1) = force(1) - a(1) * ascp_b(iL,jL)*ascp_a( kL,lL)
                  force(2) = force(2) - a(2) * ascp_b(iL,jL)*ascp_a( kL,lL)
                  force(3) = force(3) - a(3) * ascp_b(iL,jL)*ascp_a( kL,lL)
                ELSE
                  force(1) = force(1) + a(1) * ascp_a(iL,jL)*ascp_b( kL,lL)
                  force(2) = force(2) + a(2) * ascp_a(iL,jL)*ascp_b( kL,lL)
                  force(3) = force(3) + a(3) * ascp_a(iL,jL)*ascp_b( kL,lL)
                END IF
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE dscp_disp_2C

!******************************************************************************
END MODULE  scp_dispersion

