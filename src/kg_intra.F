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
!!     gt (24.10.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE kg_intra

! *****************************************************************************

  USE atomic_kind_types,        ONLY: atomic_kind_type,&
                                      get_atomic_kind,&
                                      get_atomic_kind_set
  USE basis_set_types,          ONLY: get_gto_basis_set,&
                                      gto_basis_set_type
  USE erf_fn,                   ONLY: erf
  USE global_types,             ONLY : global_environment_type
  USE kg_environment_types,     ONLY: kg_environment_type, &
                                      get_kg_env
  USE kg_force_types,           ONLY: kg_force_type 
  USE kinds,                    ONLY: int_size,&
                                      dp,&
                                      dp_size
  USE mathconstants,            ONLY: oorootpi,twopi, pi
  USE molecule_types,           ONLY : molecule_structure_type, particle_node_type, &
                                       linklist_exclusion, linklist_atoms
  USE particle_types,           ONLY: particle_type
  USE simulation_cell,          ONLY: cell_type,&
                                      pbc
  USE termination,              ONLY: stop_memory
  USE timings,                  ONLY: timeset,&
                                      timestop

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calculate_ebond_corr


! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE calculate_ebond_corr(kg_env,globenv,e_bc,force)

!   Purpose: Calculate the overlap energy of the core charge distribution.

!   History: - Creation (24.10.2002)

!   ***************************************************************************

    IMPLICIT NONE

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(kg_environment_type), INTENT(INOUT)  :: kg_env
    REAL(KIND = dp), intent (out) :: e_bc
    TYPE(kg_force_type), DIMENSION(:), POINTER, OPTIONAL  :: force

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_ebond_corr (MODULE kg_intra)"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    TYPE(cell_type), POINTER        :: box
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    REAL(KIND = dp)                        :: ann,eee,ene,enn,&
                                       fee,fne,fnn,fscalar,gnorm, &
                                       dij,idij,prefac,q_i,q_j,qnn,rijsq,&
                                       rootann,rootzee,rootzne,&
                                       zne,zee,z_i,z_j             
    INTEGER                         :: atom_a,atom_b,handle, &
                                       iat,iatom,iset,ishell,ipgf,ikind,&
                                       imol,istat,j,jatom,jkind,jshell,jset,jpgf,&
                                       nseta,nsetb,maxatom,natom_a,nkind,nmol,npart,&
                                       output_unit
    LOGICAL                         :: ionode

    REAL(KIND = dp), DIMENSION(3) :: rij

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
    TYPE(particle_type), DIMENSION(:), POINTER    :: part
    INTEGER, DIMENSION(:), POINTER :: nshella, nshellb, npgfa, npgfb
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: gcca, gccb
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: zeta,zetb

    TYPE ( linklist_exclusion ), POINTER :: llex
    TYPE ( linklist_atoms ), POINTER :: llat
    REAL(KIND = dp), DIMENSION(:), ALLOCATABLE             :: alpha,qnuc
    INTEGER, DIMENSION(:), ALLOCATABLE              :: atom_of_kind

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_ebond_corr","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    prefac = (1.0_dp/pi)**1.5_dp

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    box=box,molecule=molecule,&
                    part=part) 

    nkind = SIZE ( atomic_kind_set )
    nmol = SIZE ( molecule )
    npart = SIZE ( part )

!   *** Allocate work storage ***

    ALLOCATE (alpha(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"alpha",nkind*dp_size)
    alpha(:) = 0.0_dp

    ALLOCATE (qnuc(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"qnuc",nkind*dp_size)
    qnuc(:) = 0.0_dp

    ALLOCATE (atom_of_kind(npart),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",npart*int_size)
    atom_of_kind(:) = 0

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                atom_of_kind=atom_of_kind)

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           alpha_core_charge=alpha(ikind),&
                           qeff=qnuc(ikind))

    END DO

    e_bc = 0.0_dp

    DO imol = 1, nmol
      llat => molecule ( imol ) % ll_atoms

      DO iat = 1, molecule ( imol ) % natoms_mol
        atom_a = llat % part % p % iatom         
        atomic_kind => llat % part % p % atomic_kind
        CALL get_atomic_kind ( atomic_kind=atomic_kind, &
             orb_basis_set=orb_basis_set, kind_number=ikind)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nshell=nshella,&
                             gcc=gcca,&
                             zet=zeta)
        iatom = atom_of_kind ( atom_a )
        llex => llat % part % ex
        DO j = 1, llat % part % nexcl
          rij(1) = llat % part % p % r(1) - llex % p % r(1)
          rij(2) = llat % part % p % r(2) - llex % p % r(2)
          rij(3) = llat % part % p % r(3) - llex % p % r(3)
          rij(:) = pbc(rij(:),box)
          rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
          dij = SQRT ( rijsq )
          idij = 1.0_dp / dij
          atom_b = llex % p % iatom
          atomic_kind => llex % p % atomic_kind
          CALL get_atomic_kind ( atomic_kind=atomic_kind, &
               orb_basis_set=orb_basis_set, kind_number=jkind)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 npgf=npgfb,&
                                 nset=nsetb,&
                                 nshell=nshellb,&
                                 gcc=gccb,&
                                 zet=zetb)
          jatom = atom_of_kind ( atom_b )
          qnn = qnuc(ikind)*qnuc(jkind)
          ann = alpha(ikind)*alpha(jkind)/(alpha(ikind) + alpha(jkind))
          rootann = SQRT(ann)
          fnn = 2.0_dp*oorootpi*qnn*rootann
! getting the bc potential due to nuclear-nuclear intramolecule interacton
! the factor of 1/2 comes from double counting in the exclusion list
          enn = qnn * idij * erf (rootann*dij)
          e_bc = e_bc - enn*0.5_dp
! subtracting the force from the total force           
          fscalar = (enn - fnn*exp(-ann*rijsq))/rijsq
          IF (PRESENT (force)) THEN
            force(ikind)%f_bc(:,iatom) = force(ikind)%f_bc(:,iatom) - fscalar*rij(:)
          END IF

          DO iset = 1, nseta
            DO ishell =1, nshella(iset)
              DO ipgf= 1,npgfa(iset)

                z_i=zeta(ipgf,iset)
                gnorm = prefac*z_i**1.5_dp
                q_i = gcca(ipgf,ishell,iset)/gnorm
                zne = z_i*alpha(jkind)/(alpha(jkind) + z_i)
                rootzne = SQRT(zne)
                fne = 2.0_dp*oorootpi*qnuc(jkind)*q_i*rootzne

! getting the bc potential due to nuclear-elec intramolecule interacton
! (n-e interactions are NOT double-counted!!!)
                ene = q_i * qnuc(jkind) * idij * erf (rootzne * dij)
                e_bc = e_bc - ene
! subtracting the force from the total force           
                fscalar = (ene - fne*exp(-zne*rijsq))/rijsq
 
                IF (PRESENT (force)) THEN
                  force(ikind)%f_bc(:,iatom) = force(ikind)%f_bc(:,iatom) - fscalar*rij(:)
                  force(jkind)%f_bc(:,jatom) = force(jkind)%f_bc(:,jatom) + fscalar*rij(:)
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
                      IF (PRESENT (force)) THEN
                        force(ikind)%f_bc(:,iatom) = force(ikind)%f_bc(:,iatom) - fscalar*rij(:)
                      END IF
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
          llex => llex % next
        END DO
        llat => llat % next
      END DO
    END DO

!   *** Release work storage ***

    DEALLOCATE (alpha,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"alpha")

    DEALLOCATE (qnuc,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"qeff")

    DEALLOCATE (atom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_ebond_corr  

! *****************************************************************************

  SUBROUTINE calculate_eself(kg_env,ecore_self)

!   Purpose: Calculate the self energy of the core charge distribution.

!   History: - Creation (14.10.2002) 

!   ***************************************************************************

    TYPE(kg_environment_type), INTENT(IN) :: kg_env
    REAL (KIND = dp), intent(OUT) :: ecore_self

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    REAL(KIND = dp)                        :: alpha_core_charge,qeff
    INTEGER                         :: ikind,natom

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set

!   ---------------------------------------------------------------------------

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set) 

    ecore_self = 0.0_dp

    DO ikind=1,SIZE(atomic_kind_set)

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom,&
                           qeff=qeff,&
                           alpha_core_charge=alpha_core_charge)

      ecore_self = ecore_self - REAL(natom*qeff**2,dp)*SQRT(alpha_core_charge)

    END DO

    ecore_self = ecore_self/SQRT(twopi)

  END SUBROUTINE calculate_eself

! *****************************************************************************

END MODULE kg_intra
