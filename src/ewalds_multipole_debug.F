!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Debug routines for multipoles
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param nonbond_env ...
!> \param cell ...
!> \param particle_set ...
!> \param local_particles ...
!> \param iw ...
!> \param debug_r_space ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! *****************************************************************************
SUBROUTINE debug_ewald_multipoles(ewald_env, ewald_pw, nonbond_env, cell, &
     particle_set, local_particles, iw, debug_r_space)
  USE cell_types,                      ONLY: cell_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: matvec_3x3
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type,&
                                             random_numbers
  USE particle_types,                  ONLY: particle_type
  USE ewalds_multipole,                ONLY: ewald_multipole_evaluate
#include "./common/cp_common_uses.f90"
  IMPLICIT NONE

  TYPE charge_mono_type
     REAL(KIND=dp), DIMENSION(:),&
          POINTER                          :: charge
     REAL(KIND=dp), DIMENSION(:,:),&
          POINTER                          :: pos
  END TYPE charge_mono_type
  TYPE multi_charge_type
     TYPE(charge_mono_type), DIMENSION(:),&
          POINTER                          :: charge_typ
  END TYPE multi_charge_type
  TYPE(ewald_environment_type), POINTER    :: ewald_env
  TYPE(ewald_pw_type), POINTER             :: ewald_pw
  TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
  TYPE(cell_type), POINTER                 :: cell
  TYPE(particle_type), DIMENSION(:), &
       POINTER                                :: particle_set
  TYPE(distribution_1d_type), POINTER      :: local_particles
  INTEGER, INTENT(IN)                      :: iw
  LOGICAL, INTENT(IN)                      :: debug_r_space

  CHARACTER(len=*), PARAMETER :: routineN = 'debug_ewald_multipoles', &
       routineP = "ewalds_multipole_debug"//':'//routineN

  INTEGER                                  :: nparticles
  LOGICAL                                  :: failure
  LOGICAL, DIMENSION(3)                    :: task
  REAL(KIND=dp)                            :: e_neut, e_self, g_energy, &
                                              r_energy, debug_energy
  REAL(KIND=dp), POINTER, DIMENSION(:)     :: charges
  REAL(KIND=dp), POINTER, &
       DIMENSION(:, :)                     :: dipoles, g_forces, g_pv, &
                                              r_forces, r_pv, e_field1,&
                                              e_field2
  REAL(KIND=dp), POINTER, &
       DIMENSION(:, :, :)                  :: quadrupoles
  TYPE(rng_stream_type), POINTER           :: random_stream
  TYPE(multi_charge_type), DIMENSION(:),&
       POINTER                             :: multipoles

  failure = .FALSE.
  NULLIFY(random_stream, multipoles, charges, dipoles, g_forces, g_pv,&
          r_forces, r_pv, e_field1, e_field2)
  CALL create_rng_stream(random_stream,name="DEBUG_EWALD_MULTIPOLE",&
       distribution_type=UNIFORM)
  ! check:  charge - charge
  task    = .FALSE.
  nparticles = SIZE(particle_set)

  ! Allocate charges, dipoles, quadrupoles
  ALLOCATE(charges(nparticles))
  ALLOCATE(dipoles(3,nparticles))
  ALLOCATE(quadrupoles(3,3,nparticles))

  ! Allocate arrays for forces
  ALLOCATE(r_forces(3,nparticles))
  ALLOCATE(g_forces(3,nparticles))
  ALLOCATE(e_field1(3,nparticles))
  ALLOCATE(e_field2(3,nparticles))
  ALLOCATE(g_pv(3,3))
  ALLOCATE(r_pv(3,3))

  ! Debug CHARGES-CHARGES
  task(1) = .TRUE.
  charges = 0.0_dp
  dipoles = 0.0_dp
  quadrupoles = 0.0_dp
  r_forces = 0.0_dp
  g_forces = 0.0_dp
  e_field1 = 0.0_dp
  e_field2 = 0.0_dp
  g_pv = 0.0_dp
  r_pv = 0.0_dp
  g_energy = 0.0_dp
  r_energy = 0.0_dp
  e_neut = 0.0_dp
  e_self = 0.0_dp

  CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "CHARGE", echarge=-1.0_dp, &
       random_stream=random_stream, charges=charges)
  CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "CHARGE", echarge=1.0_dp, &
       random_stream=random_stream, charges=charges)
  CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy,&
       debug_r_space)

  WRITE(*,*)"DEBUG ENERGY (CHARGE-CHARGE): ", debug_energy
  CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
       particle_set, local_particles, g_energy, r_energy, e_neut, e_self,&
       task, do_correction_bonded=.FALSE., do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE.,&
       charges=charges,dipoles=dipoles,quadrupoles=quadrupoles,forces_local=g_forces,&
       forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
  CALL release_multi_type(multipoles)


  ! Debug CHARGES-DIPOLES
  task(1) = .TRUE.
  task(2) = .TRUE.
  charges = 0.0_dp
  dipoles = 0.0_dp
  quadrupoles = 0.0_dp
  r_forces = 0.0_dp
  g_forces = 0.0_dp
  e_field1 = 0.0_dp
  e_field2 = 0.0_dp
  g_pv = 0.0_dp
  r_pv = 0.0_dp
  g_energy = 0.0_dp
  r_energy = 0.0_dp
  e_neut = 0.0_dp
  e_self = 0.0_dp

  CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "CHARGE", echarge=-1.0_dp, &
       random_stream=random_stream, charges=charges)
  CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "DIPOLE", echarge=0.5_dp, &
       random_stream=random_stream, dipoles=dipoles)
  WRITE(*,'("CHARGES",F15.9)')charges
  WRITE(*,'("DIPOLES",3F15.9)')dipoles
  CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy, &
       debug_r_space)

  WRITE(*,*)"DEBUG ENERGY (CHARGE-DIPOLE): ", debug_energy
  CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
       particle_set, local_particles, g_energy, r_energy, e_neut, e_self,&
       task,do_correction_bonded=.FALSE.,  do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE.,&
       charges=charges,dipoles=dipoles,quadrupoles=quadrupoles,forces_local=g_forces,&
       forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
  CALL release_multi_type(multipoles)

  ! Debug DIPOLES-DIPOLES
  task(2) = .TRUE.
  charges = 0.0_dp
  dipoles = 0.0_dp
  quadrupoles = 0.0_dp
  r_forces = 0.0_dp
  g_forces = 0.0_dp
  e_field1 = 0.0_dp
  e_field2 = 0.0_dp
  g_pv = 0.0_dp
  r_pv = 0.0_dp
  g_energy = 0.0_dp
  r_energy = 0.0_dp
  e_neut = 0.0_dp
  e_self = 0.0_dp

  CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "DIPOLE", echarge=10000.0_dp, &
       random_stream=random_stream,dipoles=dipoles)
  CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "DIPOLE", echarge=20000._dp, &
       random_stream=random_stream,dipoles=dipoles)
  WRITE(*,'("DIPOLES",3F15.9)')dipoles
  CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy, &
       debug_r_space)

  WRITE(*,*)"DEBUG ENERGY (DIPOLE-DIPOLE): ", debug_energy
  CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
       particle_set, local_particles, g_energy, r_energy, e_neut, e_self,&
       task,do_correction_bonded=.FALSE.,  do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE.,&
       charges=charges,dipoles=dipoles,quadrupoles=quadrupoles,forces_local=g_forces,&
       forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
  CALL release_multi_type(multipoles)

  ! Debug CHARGES-QUADRUPOLES
  task(1) = .TRUE.
  task(3) = .TRUE.
  charges = 0.0_dp
  dipoles = 0.0_dp
  quadrupoles = 0.0_dp
  r_forces = 0.0_dp
  g_forces = 0.0_dp
  e_field1 = 0.0_dp
  e_field2 = 0.0_dp
  g_pv = 0.0_dp
  r_pv = 0.0_dp
  g_energy = 0.0_dp
  r_energy = 0.0_dp
  e_neut = 0.0_dp
  e_self = 0.0_dp

  CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "CHARGE", echarge=-1.0_dp, &
       random_stream=random_stream,charges=charges)
  CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "QUADRUPOLE", echarge=10.0_dp, &
       random_stream=random_stream,quadrupoles=quadrupoles)
  WRITE(*,'("CHARGES",F15.9)')charges
  WRITE(*,'("QUADRUPOLES",9F15.9)')quadrupoles
  CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy, &
       debug_r_space)

  WRITE(*,*)"DEBUG ENERGY (CHARGE-QUADRUPOLE): ", debug_energy
  CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
       particle_set, local_particles, g_energy, r_energy, e_neut, e_self,&
       task,do_correction_bonded=.FALSE.,  do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE.,&
       charges=charges,dipoles=dipoles,quadrupoles=quadrupoles,forces_local=g_forces,&
       forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
  CALL release_multi_type(multipoles)

  ! Debug DIPOLES-QUADRUPOLES
  task(2) = .TRUE.
  task(3) = .TRUE.
  charges = 0.0_dp
  dipoles = 0.0_dp
  quadrupoles = 0.0_dp
  r_forces = 0.0_dp
  g_forces = 0.0_dp
  e_field1 = 0.0_dp
  e_field2 = 0.0_dp
  g_pv = 0.0_dp
  r_pv = 0.0_dp
  g_energy = 0.0_dp
  r_energy = 0.0_dp
  e_neut = 0.0_dp
  e_self = 0.0_dp

  CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "DIPOLE", echarge=10000.0_dp, &
       random_stream=random_stream,dipoles=dipoles)
  CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "QUADRUPOLE", echarge=10000.0_dp, &
       random_stream=random_stream,quadrupoles=quadrupoles)
  WRITE(*,'("DIPOLES",3F15.9)')dipoles
  WRITE(*,'("QUADRUPOLES",9F15.9)')quadrupoles
  CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy, &
       debug_r_space)

  WRITE(*,*)"DEBUG ENERGY (DIPOLE-QUADRUPOLE): ", debug_energy
  CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
       particle_set, local_particles, g_energy, r_energy, e_neut, e_self,&
       task,do_correction_bonded=.FALSE.,  do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE.,&
       charges=charges,dipoles=dipoles,quadrupoles=quadrupoles,forces_local=g_forces,&
       forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
  CALL release_multi_type(multipoles)

  ! Debug QUADRUPOLES-QUADRUPOLES
  task(3) = .TRUE.
  charges = 0.0_dp
  dipoles = 0.0_dp
  quadrupoles = 0.0_dp
  r_forces = 0.0_dp
  g_forces = 0.0_dp
  e_field1 = 0.0_dp
  e_field2 = 0.0_dp
  g_pv = 0.0_dp
  r_pv = 0.0_dp
  g_energy = 0.0_dp
  r_energy = 0.0_dp
  e_neut = 0.0_dp
  e_self = 0.0_dp

  CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "QUADRUPOLE", echarge=-20000.0_dp, &
       random_stream=random_stream,quadrupoles=quadrupoles)
  CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "QUADRUPOLE", echarge=10000.0_dp, &
       random_stream=random_stream,quadrupoles=quadrupoles)
  WRITE(*,'("QUADRUPOLES",9F15.9)')quadrupoles
  CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy, &
       debug_r_space)

  WRITE(*,*)"DEBUG ENERGY (QUADRUPOLE-QUADRUPOLE): ", debug_energy
  CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
       particle_set, local_particles, g_energy, r_energy, e_neut, e_self,&
       task,do_correction_bonded=.FALSE.,  do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE.,&
       charges=charges,dipoles=dipoles,quadrupoles=quadrupoles,forces_local=g_forces,&
       forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
  CALL release_multi_type(multipoles)


  DEALLOCATE(charges)
  DEALLOCATE(dipoles)
  DEALLOCATE(quadrupoles)
  DEALLOCATE(r_forces)
  DEALLOCATE(g_forces)
  DEALLOCATE(e_field1)
  DEALLOCATE(e_field2)
  DEALLOCATE(g_pv)
  DEALLOCATE(r_pv)
  CALL delete_rng_stream(random_stream)

CONTAINS
! *****************************************************************************
!> \brief  Debug routines for multipoles - low level - charge interactions
!> \param particle_set ...
!> \param cell ...
!> \param nonbond_env ...
!> \param multipoles ...
!> \param energy ...
!> \param debug_r_space ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! *****************************************************************************
  SUBROUTINE debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles,&
       energy, debug_r_space)
    IMPLICIT NONE
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
    TYPE(multi_charge_type), DIMENSION(:),&
         POINTER                             :: multipoles
    REAL(KIND=dp), INTENT(OUT)               :: energy
    LOGICAL, INTENT(IN)                      :: debug_r_space

    CHARACTER(len=*), PARAMETER :: routineN = 'debug_ewald_multipole_low', &
      routineP = "ewalds_multipole_debug"//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: atom_a, atom_b, iend, igrp, &
                                                ikind, ilist, ipair, istart, &
                                                jkind, nkinds, npairs,k,k1,l,l1,&
                                                ncells, icell, jcell, kcell
    INTEGER, DIMENSION(:, :), POINTER        :: list
    REAL(KIND=dp)                            :: r, rab2, rab2_max, fac_ij, q
    REAL(KIND=dp), DIMENSION(3)              :: cell_v, cvi, rab, rm, rab0
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(neighbor_kind_pairs_type), POINTER  :: neighbor_kind_pair
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc

    failure = .FALSE.
    energy  = 0.0_dp
    CALL fist_nonbond_env_get (nonbond_env, nonbonded=nonbonded, natom_types = nkinds,&
         r_last_update=r_last_update,r_last_update_pbc=r_last_update_pbc)
    rab2_max = HUGE(0.0_dp)
    IF (debug_r_space) THEN
       ! This debugs the real space part of the multipole Ewald summation scheme
       ! Starting the force loop
       Lists: DO ilist=1,nonbonded%nlists
          neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
          npairs=neighbor_kind_pair%npairs
          IF (npairs ==0) CYCLE
          list  => neighbor_kind_pair%list
          cvi   =  neighbor_kind_pair%cell_vector
          CALL matvec_3x3(cell_v, cell%hmat, cvi)
          Kind_Group_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
             istart  = neighbor_kind_pair%grp_kind_start(igrp)
             iend    = neighbor_kind_pair%grp_kind_end(igrp)
             ikind   = neighbor_kind_pair%ij_kind(1,igrp)
             jkind   = neighbor_kind_pair%ij_kind(2,igrp)
             Pairs: DO ipair = istart, iend
                fac_ij = 1.0_dp
                atom_a = list(1,ipair)
                atom_b = list(2,ipair)
                IF (atom_a==atom_b) fac_ij = 0.5_dp
                rab    = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
                rab    = rab + cell_v
                rab2   = rab(1)**2 + rab(2)**2 + rab(3)**2
                IF (rab2 <= rab2_max) THEN

                   DO k = 1, SIZE(multipoles(atom_a)%charge_typ)
                      DO k1 = 1, SIZE(multipoles(atom_a)%charge_typ(k)%charge)

                         DO l = 1, SIZE(multipoles(atom_b)%charge_typ)
                            DO l1 = 1, SIZE(multipoles(atom_b)%charge_typ(l)%charge)

                               rm = rab + multipoles(atom_b)%charge_typ(l)%pos(:,l1) - multipoles(atom_a)%charge_typ(k)%pos(:,k1)
                               r  = SQRT(DOT_PRODUCT(rm,rm))
                               q  = multipoles(atom_b)%charge_typ(l)%charge(l1) * multipoles(atom_a)%charge_typ(k)%charge(k1)
                               energy = energy + q / r * fac_ij
                            END DO
                         END DO

                      END DO
                   END DO

                END IF
             END DO Pairs
          END DO Kind_Group_Loop
       END DO Lists
    ELSE
       ncells = 6
       !Debugs the sum of real + space terms.. (Charge-Charge and Charge-Dipole should be anyway wrong but
       !all the other terms should be correct)
       DO atom_a = 1, SIZE(particle_set)
       DO atom_b = atom_a, SIZE(particle_set)
          fac_ij = 1.0_dp
          IF (atom_a==atom_b) fac_ij = 0.5_dp
          rab0   = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
          ! Loop over cells
          DO icell = -ncells,ncells
          DO jcell = -ncells,ncells
          DO kcell = -ncells,ncells
             cell_v= MATMUL(cell%hmat,REAL((/icell,jcell,kcell/),KIND=dp))
             IF (ALL(cell_v==0.0_dp).AND.(atom_a==atom_b)) CYCLE
             rab = rab0 + cell_v
             rab2   = rab(1)**2 + rab(2)**2 + rab(3)**2
             IF (rab2 <= rab2_max) THEN

                DO k = 1, SIZE(multipoles(atom_a)%charge_typ)
                   DO k1 = 1, SIZE(multipoles(atom_a)%charge_typ(k)%charge)

                      DO l = 1, SIZE(multipoles(atom_b)%charge_typ)
                         DO l1 = 1, SIZE(multipoles(atom_b)%charge_typ(l)%charge)

                            rm = rab + multipoles(atom_b)%charge_typ(l)%pos(:,l1) - multipoles(atom_a)%charge_typ(k)%pos(:,k1)
                            r  = SQRT(DOT_PRODUCT(rm,rm))
                            q  = multipoles(atom_b)%charge_typ(l)%charge(l1) * multipoles(atom_a)%charge_typ(k)%charge(k1)
                            energy = energy + q / r * fac_ij
                         END DO
                      END DO

                   END DO
                END DO

             END IF
          END DO
          END DO
          END DO
       END DO
       END DO
    END IF
  END SUBROUTINE debug_ewald_multipole_low

! *****************************************************************************
!> \brief  create multi_type for multipoles
!> \param multipoles ...
!> \param idim ...
!> \param istart ...
!> \param iend ...
!> \param label ...
!> \param echarge ...
!> \param random_stream ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! *****************************************************************************
  SUBROUTINE create_multi_type(multipoles, idim, istart, iend, label, echarge,&
       random_stream, charges, dipoles, quadrupoles)
    IMPLICIT NONE
    TYPE(multi_charge_type), DIMENSION(:),&
         POINTER                             :: multipoles
    INTEGER, INTENT(IN)                      :: idim, istart, iend
    CHARACTER(LEN=*), INTENT(IN)             :: label
    REAL(KIND=dp), INTENT(IN)                :: echarge
    TYPE(rng_stream_type), POINTER           :: random_stream
    REAL(KIND=dp), DIMENSION(:), POINTER,&
         OPTIONAL                            :: charges
    REAL(KIND=dp), DIMENSION(:,:), POINTER,&
         OPTIONAL                            :: dipoles
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER,&
         OPTIONAL                            :: quadrupoles

    CHARACTER(len=*), PARAMETER :: routineN = 'create_multi_type', &
      routineP = "ewalds_multipole_debug"//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: i, isize, k, l, m
    REAL(KIND=dp)                            :: dx, rvec(3), rvec1(3), rvec2(3), r2

    failure = .FALSE.
    IF (ASSOCIATED(multipoles)) THEN
       CPASSERT(SIZE(multipoles)==idim)
    ELSE
       ALLOCATE(multipoles(idim))
       DO i = 1, idim
          NULLIFY(multipoles(i)%charge_typ)
       END DO
    END IF
    DO i = istart, iend
       IF (ASSOCIATED(multipoles(i)%charge_typ)) THEN
          ! make a copy of the array and enlarge the array type by 1
          isize = SIZE(multipoles(i)%charge_typ)+1
       ELSE
          isize = 1
       END IF
       CALL reallocate_charge_type(multipoles(i)%charge_typ,1,isize)
       SELECT CASE(label)
       CASE("CHARGE")
          CPASSERT(PRESENT(charges))
          CPASSERT(ASSOCIATED(charges))
          ALLOCATE(multipoles(i)%charge_typ(isize)%charge(1))
          ALLOCATE(multipoles(i)%charge_typ(isize)%pos(3,1))

          multipoles(i)%charge_typ(isize)%charge(1)  = echarge
          multipoles(i)%charge_typ(isize)%pos(1:3,1) = 0.0_dp
          charges(i) = charges(i) + echarge
       CASE("DIPOLE")
          dx = 1.0E-4_dp
          CPASSERT(PRESENT(dipoles))
          CPASSERT(ASSOCIATED(dipoles))
          ALLOCATE(multipoles(i)%charge_typ(isize)%charge(2))
          ALLOCATE(multipoles(i)%charge_typ(isize)%pos(3,2))
          CALL random_numbers(rvec, random_stream)
          rvec = rvec/(2.0_dp*SQRT(DOT_PRODUCT(rvec,rvec)))*dx
          multipoles(i)%charge_typ(isize)%charge(1)  = echarge
          multipoles(i)%charge_typ(isize)%pos(1:3,1) = rvec
          multipoles(i)%charge_typ(isize)%charge(2)  =-echarge
          multipoles(i)%charge_typ(isize)%pos(1:3,2) =-rvec

          dipoles(:,i) = dipoles(:,i) + 2.0_dp*echarge*rvec
       CASE("QUADRUPOLE")
          dx = 1.0E-2_dp
          CPASSERT(PRESENT(quadrupoles))
          CPASSERT(ASSOCIATED(quadrupoles))
          ALLOCATE(multipoles(i)%charge_typ(isize)%charge(4))
          ALLOCATE(multipoles(i)%charge_typ(isize)%pos(3,4))
          CALL random_numbers(rvec1, random_stream)
          CALL random_numbers(rvec2, random_stream)
          rvec1 = rvec1/SQRT(DOT_PRODUCT(rvec1,rvec1))
          rvec2 = rvec2 - DOT_PRODUCT(rvec2,rvec1)*rvec1
          rvec2 = rvec2/SQRT(DOT_PRODUCT(rvec2,rvec2))
          !
          rvec1 = rvec1/2.0_dp*dx
          rvec2 = rvec2/2.0_dp*dx
          !       + (4)  ^      - (1)
          !              |rvec2
          !              |
          !              0------> rvec1
          !
          !
          !       - (3)         + (2)
          multipoles(i)%charge_typ(isize)%charge(1)  = -echarge
          multipoles(i)%charge_typ(isize)%pos(1:3,1) =  rvec1+rvec2
          multipoles(i)%charge_typ(isize)%charge(2)  =  echarge
          multipoles(i)%charge_typ(isize)%pos(1:3,2) =  rvec1-rvec2
          multipoles(i)%charge_typ(isize)%charge(3)  = -echarge
          multipoles(i)%charge_typ(isize)%pos(1:3,3) = -rvec1-rvec2
          multipoles(i)%charge_typ(isize)%charge(4)  =  echarge
          multipoles(i)%charge_typ(isize)%pos(1:3,4) = -rvec1+rvec2

          DO k = 1, 4
             r2 = DOT_PRODUCT(multipoles(i)%charge_typ(isize)%pos(:,k),multipoles(i)%charge_typ(isize)%pos(:,k))
             DO l = 1, 3
                DO m = 1, 3
                   quadrupoles(m,l,i) = quadrupoles(m,l,i) +  3.0_dp* 0.5_dp *multipoles(i)%charge_typ(isize)%charge(k) * &
                                                              multipoles(i)%charge_typ(isize)%pos(l,k)  * &
                                                              multipoles(i)%charge_typ(isize)%pos(m,k)
                   IF (m==l) quadrupoles(m,l,i) = quadrupoles(m,l,i) -  0.5_dp * multipoles(i)%charge_typ(isize)%charge(k) * r2
                END DO
             END DO
          END DO

       END SELECT
    END DO
  END SUBROUTINE create_multi_type

! *****************************************************************************
!> \brief  release multi_type for multipoles
!> \param multipoles ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! *****************************************************************************
  SUBROUTINE release_multi_type(multipoles)
    IMPLICIT NONE
    TYPE(multi_charge_type), DIMENSION(:),&
         POINTER                             :: multipoles

    CHARACTER(len=*), PARAMETER :: routineN = 'release_multi_type', &
      routineP = "ewalds_multipole_debug"//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: i, j

    failure = .FALSE.
    IF (ASSOCIATED(multipoles)) THEN
       DO i = 1, SIZE(multipoles)
          DO j = 1, SIZE(multipoles(i)%charge_typ)
             DEALLOCATE(multipoles(i)%charge_typ(j)%charge)
             DEALLOCATE(multipoles(i)%charge_typ(j)%pos)
          END DO
          DEALLOCATE(multipoles(i)%charge_typ)
       END DO
    END IF
  END SUBROUTINE release_multi_type

! *****************************************************************************
!> \brief  reallocates multi_type for multipoles
!> \param charge_typ ...
!> \param istart ...
!> \param iend ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! *****************************************************************************
  SUBROUTINE reallocate_charge_type(charge_typ, istart, iend)
    IMPLICIT NONE
    TYPE(charge_mono_type), DIMENSION(:),&
         POINTER                             :: charge_typ
    INTEGER, INTENT(IN)                      :: istart, iend

    CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_charge_type', &
      routineP = "ewalds_multipole_debug"//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: i, isize, jsize, jsize1, jsize2, j
        TYPE(charge_mono_type), DIMENSION(:),&
         POINTER                             :: charge_typ_bk

    failure = .FALSE.
    IF (ASSOCIATED(charge_typ)) THEN
       isize = SIZE(charge_typ)
       ALLOCATE(charge_typ_bk(1:isize))
       DO j = 1, isize
          jsize = SIZE(charge_typ(j)%charge)
          ALLOCATE(charge_typ_bk(j)%charge(jsize))
          jsize1 = SIZE(charge_typ(j)%pos,1)
          jsize2 = SIZE(charge_typ(j)%pos,2)
          ALLOCATE(charge_typ_bk(j)%pos(jsize1,jsize2))
          charge_typ_bk(j)%pos    = charge_typ(j)%pos
          charge_typ_bk(j)%charge = charge_typ(j)%charge
       END DO
       DO j = 1, SIZE(charge_typ)
          DEALLOCATE(charge_typ(j)%charge)
          DEALLOCATE(charge_typ(j)%pos)
       END DO
       DEALLOCATE(charge_typ)
       ! Reallocate
       ALLOCATE(charge_typ_bk(istart:iend))
       DO i = istart, isize
          jsize = SIZE(charge_typ_bk(j)%charge)
          ALLOCATE(charge_typ(j)%charge(jsize))
          jsize1 = SIZE(charge_typ_bk(j)%pos,1)
          jsize2 = SIZE(charge_typ_bk(j)%pos,2)
          ALLOCATE(charge_typ(j)%pos(jsize1,jsize2))
          charge_typ(j)%pos    = charge_typ_bk(j)%pos
          charge_typ(j)%charge = charge_typ_bk(j)%charge
       END DO
       DO j = 1, SIZE(charge_typ_bk)
          DEALLOCATE(charge_typ_bk(j)%charge)
          DEALLOCATE(charge_typ_bk(j)%pos)
       END DO
       DEALLOCATE(charge_typ_bk)
    ELSE
       ALLOCATE(charge_typ(istart:iend))
    END IF

   END SUBROUTINE reallocate_charge_type

END SUBROUTINE debug_ewald_multipoles

! *****************************************************************************
!> \brief  Routine to debug potential, field and electric field gradients
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param nonbond_env ...
!> \param cell ...
!> \param particle_set ...
!> \param local_particles ...
!> \param radii ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \param task ...
!> \param iw ...
!> \param atomic_kind_set ...
!> \param mm_section ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! *****************************************************************************
SUBROUTINE debug_ewald_multipoles_fields(ewald_env, ewald_pw, nonbond_env, cell,&
     particle_set, local_particles, radii, charges, dipoles, quadrupoles, task, iw, &
     atomic_kind_set, mm_section)
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE ewalds_multipole,                ONLY: ewald_multipole_evaluate
  USE fist_neighbor_list_control,      ONLY: list_control
#include "./common/cp_common_uses.f90"
  IMPLICIT NONE
  TYPE(ewald_environment_type), POINTER    :: ewald_env
  TYPE(ewald_pw_type), POINTER             :: ewald_pw
  TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
  TYPE(cell_type), POINTER                 :: cell
  TYPE(particle_type), POINTER             :: particle_set(:)
  TYPE(distribution_1d_type), POINTER      :: local_particles
  REAL(KIND=dp), DIMENSION(:), &
       POINTER, OPTIONAL                   :: radii, charges
  REAL(KIND=dp), DIMENSION(:, :), &
       POINTER, OPTIONAL                   :: dipoles
  REAL(KIND=dp), DIMENSION(:, :, :), &
       POINTER, OPTIONAL                   :: quadrupoles
  LOGICAL, DIMENSION(3), INTENT(IN)        :: task
  INTEGER, INTENT(IN)                      :: iw
  TYPE(atomic_kind_type), POINTER          :: atomic_kind_set(:)
  TYPE(section_vals_type), POINTER         :: mm_section

  CHARACTER(len=*), PARAMETER :: routineN = 'debug_ewald_multipoles_fields', &
      routineP = "ewalds_multipole_debug"//':'//routineN

  REAL(KIND=dp), POINTER, DIMENSION(:)     :: lcharges, efield0
  REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: ldipoles
  REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: lquadrupoles
  REAL(KIND=dp)                            :: dq, energy_local, energy_glob, e_neut,&
       e_self, pv_local(3,3), pv_glob(3,3), ene(2), pot, coord(3), efield1n(3), tot_ene,&
       o_tot_ene, enev(3,2), efield2n(3,3)
  REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE&
                                           :: forces_local, forces_glob, efield1, efield2
  INTEGER :: nparticles,i,k,j,nparticle_local, iparticle_kind
  TYPE(particle_type), POINTER, DIMENSION(:) :: shell_particle_set, core_particle_set
  TYPE(cp_logger_type), POINTER              :: logger

  NULLIFY(lcharges, ldipoles, lquadrupoles, shell_particle_set, core_particle_set)
  NULLIFY(logger)
  logger => cp_get_default_logger()

  nparticles = SIZE(particle_set)
  nparticle_local = 0
  DO iparticle_kind=1,SIZE(local_particles%n_el)
     nparticle_local = nparticle_local + local_particles%n_el(iparticle_kind)
  END DO
  ALLOCATE(lcharges(nparticles))
  ALLOCATE(forces_glob(3,nparticles))
  ALLOCATE(forces_local(3,nparticle_local))
  ALLOCATE(efield0(nparticles))
  ALLOCATE(efield1(3,nparticles))
  ALLOCATE(efield2(9,nparticles))
  forces_glob = 0.0_dp
  forces_local= 0.0_dp
  efield0     = 0.0_dp
  efield1     = 0.0_dp
  efield2     = 0.0_dp
  pv_local    = 0.0_dp
  pv_glob     = 0.0_dp
  energy_glob = 0.0_dp
  energy_local= 0.0_dp
  e_neut      = 0.0_dp
  e_self      = 0.0_dp
  CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, particle_set,&
       local_particles, energy_local, energy_glob, e_neut, e_self, task, .FALSE., .TRUE., .TRUE.,&
       .TRUE., radii, charges, dipoles, quadrupoles, forces_local, forces_glob, pv_local, pv_glob,&
       efield0, efield1, efield2, iw, do_debug=.FALSE.)
  o_tot_ene = energy_local + energy_glob + e_neut + e_self
  WRITE(*,*)"TOTAL ENERGY :: ========>",o_tot_ene
  ! Debug Potential
  dq      = 0.001_dp
  tot_ene = 0.0_dp
  DO i = 1, nparticles
     DO k = 1, 2
        lcharges    = charges
        lcharges(i) = charges(i) + (-1.0_dp)**k * dq
        forces_glob = 0.0_dp
        forces_local= 0.0_dp
        pv_local    = 0.0_dp
        pv_glob     = 0.0_dp
        energy_glob = 0.0_dp
        energy_local= 0.0_dp
        e_neut      = 0.0_dp
        e_self      = 0.0_dp
        CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, particle_set,&
             local_particles, energy_local, energy_glob, e_neut, e_self, task, .FALSE., .FALSE., .FALSE.,&
             .FALSE., radii, lcharges, dipoles, quadrupoles, iw=iw, do_debug=.FALSE.)
        ene(k) = energy_local + energy_glob + e_neut + e_self
     END DO
     pot = (ene(2)-ene(1))/(2.0_dp*dq)
     WRITE(*,'(A,I8,3(A,F15.9))')"POTENTIAL FOR ATOM: ",i, " NUMERICAL: ",pot, " ANALYTICAL: ",efield0(i),&
          " ERROR: ",pot-efield0(i)
     tot_ene = tot_ene + 0.5_dp*efield0(i)*charges(i)
  END DO
  WRITE(*,*)"ENERGIES: ",o_tot_ene, tot_ene, o_tot_ene-tot_ene
  WRITE(*,'(/,/,/)')
  ! Debug Field
  dq = 0.001_dp
  DO i = 1, nparticles
     coord = particle_set(i)%r
     DO j = 1, 3
        DO k = 1, 2
           particle_set(i)%r(j) = coord(j) + (-1.0_dp)**k * dq

           ! Rebuild neighbor lists
           CALL list_control ( atomic_kind_set, particle_set, local_particles, &
                cell, nonbond_env, logger%para_env, mm_section, &
                shell_particle_set, core_particle_set)

           forces_glob = 0.0_dp
           forces_local= 0.0_dp
           pv_local    = 0.0_dp
           pv_glob     = 0.0_dp
           energy_glob = 0.0_dp
           energy_local= 0.0_dp
           e_neut      = 0.0_dp
           e_self      = 0.0_dp
           efield0     = 0.0_dp
           CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, particle_set,&
                local_particles, energy_local, energy_glob, e_neut, e_self, task, .FALSE., .TRUE., .TRUE.,&
                .TRUE., radii, charges, dipoles, quadrupoles, forces_local, forces_glob, pv_local, pv_glob,&
                efield0, iw=iw, do_debug=.FALSE.)
           ene(k) = efield0(i)
           particle_set(i)%r(j) = coord(j)
        END DO
        efield1n(j) = -(ene(2)-ene(1))/(2.0_dp*dq)
     END DO
     WRITE(*,'(/,A,I8)')"FIELD FOR ATOM: ",i
     WRITE(*,'(A,3F15.9)')" NUMERICAL: ",efield1n, " ANALYTICAL: ",efield1(:,i),&
          " ERROR: ",efield1n-efield1(:,i)
     IF (task(2)) THEN
        tot_ene = tot_ene - 0.5_dp*DOT_PRODUCT(efield1(:,i),dipoles(:,i))
     END IF
  END DO
  WRITE(*,*)"ENERGIES: ",o_tot_ene, tot_ene, o_tot_ene-tot_ene

 ! Debug Field Gradient
  dq = 0.0001_dp
  DO i = 1, nparticles
     coord = particle_set(i)%r
     DO j = 1, 3
        DO k = 1, 2
           particle_set(i)%r(j) = coord(j) + (-1.0_dp)**k * dq

           ! Rebuild neighbor lists
           CALL list_control ( atomic_kind_set, particle_set, local_particles, &
                cell, nonbond_env, logger%para_env, mm_section, &
                shell_particle_set, core_particle_set)

           forces_glob = 0.0_dp
           forces_local= 0.0_dp
           pv_local    = 0.0_dp
           pv_glob     = 0.0_dp
           energy_glob = 0.0_dp
           energy_local= 0.0_dp
           e_neut      = 0.0_dp
           e_self      = 0.0_dp
           efield1     = 0.0_dp
           CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, particle_set,&
                local_particles, energy_local, energy_glob, e_neut, e_self, task,.FALSE., .TRUE., .TRUE.,&
                .TRUE., radii, charges, dipoles, quadrupoles, forces_local, forces_glob, pv_local, pv_glob,&
                efield1=efield1, iw=iw, do_debug=.FALSE.)
           enev(:,k) = efield1(:,i)
           particle_set(i)%r(j) = coord(j)
        END DO
        efield2n(:,j) = (enev(:,2)-enev(:,1))/(2.0_dp*dq)
     END DO
     WRITE(*,'(/,A,I8)')"FIELD GRADIENT FOR ATOM: ",i
     WRITE(*,'(A,9F15.9)')" NUMERICAL:  ",efield2n,&
                          " ANALYTICAL: ",efield2(:,i),&
                          " ERROR:      ",RESHAPE(efield2n,(/9/))-efield2(:,i)
  END DO
END SUBROUTINE debug_ewald_multipoles_fields

! *****************************************************************************
!> \brief  Routine to debug potential, field and electric field gradients
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param nonbond_env ...
!> \param cell ...
!> \param particle_set ...
!> \param local_particles ...
!> \param radii ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \param task ...
!> \param iw ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! *****************************************************************************
SUBROUTINE debug_ewald_multipoles_fields2(ewald_env, ewald_pw, nonbond_env, cell,&
     particle_set, local_particles, radii, charges, dipoles, quadrupoles, task, iw)
  USE cell_types,                      ONLY: cell_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_type
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                             neighbor_kind_pairs_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE ewalds_multipole,                ONLY: ewald_multipole_evaluate
  USE fist_neighbor_list_control,      ONLY: list_control
#include "./common/cp_common_uses.f90"
  IMPLICIT NONE
  TYPE(ewald_environment_type), POINTER    :: ewald_env
  TYPE(ewald_pw_type), POINTER             :: ewald_pw
  TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
  TYPE(cell_type), POINTER                 :: cell
  TYPE(particle_type), POINTER             :: particle_set(:)
  TYPE(distribution_1d_type), POINTER      :: local_particles
  REAL(KIND=dp), DIMENSION(:), &
       POINTER, OPTIONAL                   :: radii, charges
  REAL(KIND=dp), DIMENSION(:, :), &
       POINTER, OPTIONAL                   :: dipoles
  REAL(KIND=dp), DIMENSION(:, :, :), &
       POINTER, OPTIONAL                   :: quadrupoles
  LOGICAL, DIMENSION(3), INTENT(IN)        :: task
  INTEGER, INTENT(IN)                      :: iw

  CHARACTER(len=*), PARAMETER :: routineN = 'debug_ewald_multipoles_fields', &
      routineP = "ewalds_multipole_debug"//':'//routineN

  REAL(KIND=dp), POINTER, DIMENSION(:)     :: efield0
  REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: ldipoles
  REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: lquadrupoles
  REAL(KIND=dp)                            :: energy_local, energy_glob, e_neut,&
       e_self, pv_local(3,3), pv_glob(3,3), tot_ene, o_tot_ene, prod
  REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE&
                                           :: forces_local, forces_glob, efield1, efield2
  INTEGER :: nparticles,i,k,j,nparticle_local,iparticle_kind,ind
  TYPE(particle_type), POINTER, DIMENSION(:) :: shell_particle_set, core_particle_set
  TYPE(cp_logger_type), POINTER              :: logger

  NULLIFY(ldipoles, lquadrupoles, shell_particle_set, core_particle_set)
  NULLIFY(logger)
  logger => cp_get_default_logger()

  nparticles = SIZE(particle_set)
  nparticle_local = 0
  DO iparticle_kind=1,SIZE(local_particles%n_el)
     nparticle_local = nparticle_local + local_particles%n_el(iparticle_kind)
  END DO
  ALLOCATE(forces_glob(3,nparticles))
  ALLOCATE(forces_local(3,nparticle_local))
  ALLOCATE(efield0(nparticles))
  ALLOCATE(efield1(3,nparticles))
  ALLOCATE(efield2(9,nparticles))
  forces_glob = 0.0_dp
  forces_local= 0.0_dp
  efield0     = 0.0_dp
  efield1     = 0.0_dp
  efield2     = 0.0_dp
  pv_local    = 0.0_dp
  pv_glob     = 0.0_dp
  energy_glob = 0.0_dp
  energy_local= 0.0_dp
  e_neut      = 0.0_dp
  e_self      = 0.0_dp
  CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, particle_set,&
       local_particles, energy_local, energy_glob, e_neut, e_self, task,.FALSE.,.TRUE., .TRUE.,&
       .TRUE., radii, charges, dipoles, quadrupoles, forces_local, forces_glob, pv_local, pv_glob,&
       efield0, efield1, efield2, iw, do_debug=.FALSE.)
  o_tot_ene = energy_local + energy_glob + e_neut + e_self
  WRITE(*,*)"TOTAL ENERGY :: ========>",o_tot_ene

  ! Debug Potential
  tot_ene = 0.0_dp
  IF (task(1)) THEN
     DO i = 1, nparticles
        tot_ene = tot_ene + 0.5_dp*efield0(i)*charges(i)
     END DO
     WRITE(*,*)"ENERGIES: ",o_tot_ene, tot_ene, o_tot_ene-tot_ene
     WRITE(*,'(/,/,/)')
  END IF

  ! Debug Field
  IF (task(2)) THEN
     DO i = 1, nparticles
        tot_ene = tot_ene - 0.5_dp*DOT_PRODUCT(efield1(:,i),dipoles(:,i))
     END DO
     WRITE(*,*)"ENERGIES: ",o_tot_ene, tot_ene, o_tot_ene-tot_ene
     WRITE(*,'(/,/,/)')
  END IF

 ! Debug Field Gradient
  IF (task(3)) THEN
     DO i = 1, nparticles
        ind  = 0
        prod = 0.0_dp
        DO j = 1,3
           DO k = 1,3
              ind = ind + 1
              prod = prod + efield2(ind,i)*quadrupoles(j,k,i)
           END DO
        END DO
        tot_ene = tot_ene - 0.5_dp*(1.0_dp/3.0_dp)*prod
     END DO
     WRITE(*,*)"ENERGIES: ",o_tot_ene, tot_ene, o_tot_ene-tot_ene
     WRITE(*,'(/,/,/)')
  END IF

END SUBROUTINE debug_ewald_multipoles_fields2
