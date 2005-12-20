!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/KG_intra[1.0] *
!!
!!   NAME
!!     kg_intra
!!
!!   FUNCTION
!!     bonded_correction due to the non-interaction of gaussian charges on 
!!     the same molecule. Calculates the nuclear bonded_correction
!!     and the frozen_density bonded correction
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE kg_intra

! *****************************************************************************

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE erf_fn,                          ONLY: erf
  USE global_types,                    ONLY: global_environment_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_force_types,                  ONLY: kg_force_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE mathconstants,                   ONLY: oorootpi,&
                                             pi,&
                                             twopi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calculate_ebond_corr


! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE calculate_ebond_corr(kg_env,globenv,e_bc,calculate_forces)

!   Purpose: Calculate the overlap energy of the core charge distribution.

!   History: - Creation (24.10.2002)

!   ***************************************************************************


    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      pointer                             :: globenv
    REAL(KIND=dp), intent(out)               :: e_bc
    LOGICAL, INTENT(IN)                      :: calculate_forces

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE calculate_ebond_corr (MODULE kg_intra)"

    INTEGER :: atom_a, atom_b, handle, i, iat, iatom, iexl, ikind, ipgf, iset, &
      ishell, istat, j, jatom, jkind, jpgf, jset, jshell, group, nkind, &
      nexl, nparticle_local, npart, nseta, nsetb, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: npgfa, npgfb, nshella, nshellb
    LOGICAL                                  :: ionode
    REAL(KIND=dp) :: alpha_i, alpha_j, ann, &
      dij, eee, ene, enn, fee, fne, fnn, fscalar, gnorm,&
      idij, prefac, q_i, q_j, qnn, qnuc_i, qnuc_j,&
      rijsq, rootann, rootzee, rootzne, z_i, &
      z_j, zee, zne
    REAL(KIND=dp), DIMENSION(3)              :: ri, rj, rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcca, gccb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: box
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, DIMENSION ( : ), POINTER        :: list

!   ---------------------------------------------------------------------------
    CALL write_checkpoint_information("entering "//routine,globenv)

    CALL timeset("calculate_ebond_corr","I","",handle)

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
    prefac = (1.0_dp/pi)**1.5_dp

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=box,&
                    local_particles=local_particles,&
                    particle_set=particle_set) 

    IF(calculate_forces) THEN
      CALL get_kg_env(kg_env=kg_env,&
                    force=force) 
    END IF 
    nkind = SIZE ( atomic_kind_set )
    npart = SIZE ( particle_set )

!   *** Allocate work storage ***

    ALLOCATE (atom_of_kind(npart),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",npart*int_size)
    atom_of_kind(:) = 0

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                atom_of_kind=atom_of_kind)
    e_bc = 0.0_dp

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           alpha_core_charge=alpha_i,&
                           orb_basis_set=orb_basis_set,&                   
                           qeff=qnuc_i)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nshell=nshella,&
                             gcc=gcca,&
                             zet=zeta)
      nparticle_local = local_particles % n_el ( ikind )
      DO iatom = 1, nparticle_local
        i = local_particles % list ( ikind ) % array ( iatom )
        atom_a=atom_of_kind(i)
        ri ( : ) = particle_set ( i ) % r ( : )
        DO jkind = 1, nkind
          list => particle_set ( i ) % lpi % exclusion ( jkind ) % list_ei
          IF (.NOT.ASSOCIATED(list)) CYCLE
          nexl = SIZE ( list )
          DO iexl = 1, nexl
            j =  list ( iexl )
            IF ( j==i ) CYCLE
            rj ( : ) = particle_set ( j ) % r ( : )
            atom_b=atom_of_kind(j)
            atomic_kind => atomic_kind_set(jkind)
            rij(1) = ri(1) - rj(1)
            rij(2) = ri(2) - rj(2)
            rij(3) = ri(3) - rj(3)
            rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
            dij = SQRT ( rijsq )
            idij = 1.0_dp / dij
            CALL get_atomic_kind ( atomic_kind=atomic_kind, &
                 alpha_core_charge=alpha_j,&
                 orb_basis_set=orb_basis_set, qeff=qnuc_j)            
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                   npgf=npgfb,&
                                   nset=nsetb,&
                                   nshell=nshellb,&
                                   gcc=gccb,&
                                   zet=zetb)
            qnn = qnuc_i*qnuc_j
            ann = alpha_i*alpha_j/(alpha_i + alpha_j)
            rootann = SQRT(ann)
            fnn = 2.0_dp*oorootpi*qnn*rootann
! getting the bc potential due to nuclear-nuclear intramolecule interacton
! the factor of 1/2 comes from double counting in the exclusion list
            enn = qnn * idij * erf (rootann*dij)
            e_bc = e_bc - enn*0.5_dp
! subtracting the force from the total force
            fscalar = (enn - fnn*exp(-ann*rijsq))/rijsq
            IF (calculate_forces) THEN
              force(ikind)%f_bc(:,atom_a) = force(ikind)%f_bc(:,atom_a) - fscalar*rij(:)
            END IF

            DO iset = 1, nseta
              DO ishell =1, nshella(iset)
                DO ipgf= 1,npgfa(iset)

                  z_i=zeta(ipgf,iset)
                  gnorm = prefac*z_i**1.5_dp
                  q_i = gcca(ipgf,ishell,iset)/gnorm
                  zne = z_i*alpha_j/(alpha_j + z_i)
                  rootzne = SQRT(zne)
                  fne = 2.0_dp*oorootpi*qnuc_j*q_i*rootzne

! getting the bc potential due to nuclear-elec intramolecule interacton
! (n-e interactions are NOT double-counted!!!)
                  ene = q_i * qnuc_j * idij * erf (rootzne * dij)
                  e_bc = e_bc - ene
! subtracting the force from the total force
                  fscalar = (ene - fne*exp(-zne*rijsq))/rijsq

                  IF (calculate_forces) THEN
                    force(ikind)%f_bc(:,atom_a) = force(ikind)%f_bc(:,atom_a) - fscalar*rij(:)
                    force(jkind)%f_bc(:,atom_b) = force(jkind)%f_bc(:,atom_b) + fscalar*rij(:)
                  END IF
                  DO jset = 1, nsetb
                    DO jshell =1, nshellb(jset)
                      DO jpgf= 1,npgfb(jset)
                        z_j=zetb(jpgf,jset)
                        gnorm = prefac*z_j**1.5_dp
                        q_j = gccb(jpgf,jshell,jset)/gnorm
                        zee = z_i*z_j/(z_i+z_j)
                        rootzee=SQRT(zee)
                        fee = 2.0_dp*oorootpi*q_i*q_j*rootzee

! getting the bc potential due to elec-elec intramolecule interacton
! the factor of 1/2 comes from double counting in the exclusion list
                        eee = q_i * q_j * idij * erf (rootzee*dij)
                        e_bc = e_bc - eee*0.5_dp
! subtracting the force from the total force
                        fscalar = (eee - fee*exp(-zee*rijsq))/rijsq
                        IF (calculate_forces) THEN
                          force(ikind)%f_bc(:,atom_a) = force(ikind)%f_bc(:,atom_a) &
                                                        - fscalar*rij(:)
                        END IF
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

!   *** Release work storage ***

    DEALLOCATE (atom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")

    CALL mp_sum(e_bc,group)
    CALL write_checkpoint_information("leaving "//routine,globenv)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_ebond_corr  

! *****************************************************************************

END MODULE kg_intra
