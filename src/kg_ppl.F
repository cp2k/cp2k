!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_ppl [1.1] *
!!
!!     b) Pseudopotential calculation (Goedecker, Teter and Hutter; GTH):
!!
!!        <a|V|b> = <a|(V(local) b>
!!
!!        <a|V(local)|b> = <a|-Z(eff)*erf(SQRT(2)*alpha*r)/r +
!!                            (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!!                             C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>
!!
!!        <a|V(non-local)|b> = O                            
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE kg_ppl               
  USE ai_overlap_ppl,                  ONLY: overlap_ppl
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             dbl,&
                                             dp_size
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset
  USE orbital_symbols,                 ONLY: cgf_symbol
  USE particle_types,                  ONLY: particle_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_force_types,                  ONLY: kg_force_type
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_ppl"



! *** Public subroutines ***

  PUBLIC :: calculate_rho0_ppl, calculate_drho_ppl

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE calculate_rho0_ppl(kg_env,globenv,energy,calculate_forces,force)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(kg_environment_type), POINTER        :: kg_env
    REAL(dbl), INTENT (out)                   :: energy
    LOGICAL, INTENT(IN)                       :: calculate_forces
    TYPE(kg_force_type), DIMENSION(:), POINTER, OPTIONAL            :: force

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_rho0_pp"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    TYPE(neighbor_list_type), POINTER               :: sac_ppl_neighbor_list
    TYPE(neighbor_node_type), POINTER               :: sac_ppl_neighbor_node

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sac_ppl
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set

    INTEGER                           :: ldsab,maxder,maxsgf,nder,nkind
    INTEGER  :: ab,ac,atom_a,atom_c,handle,&
                i,iatom,ikind,ilist,inode,ipgf,iset,ishell,istat,jkind,&
                katom,kkind,knode,kneigbor,l,lc_max,lc_min,maxco,maxl,&
                maxlgto,maxlppl,natom,natom_of_kind,ncoa,&
                nlist,nneighbor,nnode,nseta,sgfa

    INTEGER, DIMENSION(:), ALLOCATABLE    :: atom_of_kind,kind_of
    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,npgfa,nsgfa 
    INTEGER, DIMENSION(:), POINTER                :: atom_list
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa

    REAL(dbl) :: dab,dac,dbc,rab2,rac2,rbc2
    REAL(dbl), DIMENSION(1) :: zetc,zetb,rpgfb
    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: hab,pab,work

    REAL(dbl), DIMENSION(3)            :: force_a
    REAL(dbl), DIMENSION(:), POINTER   :: set_radius_a
    REAL(dbl), DIMENSION(:,:), POINTER :: rpgfa,sphi_a,zeta

    REAL(dbl), DIMENSION(3)                  :: rab,rac,rbc
    REAL(dbl), DIMENSION(:), ALLOCATABLE     :: alpha_ppl,ppl_radius

    LOGICAL  :: gth_potential_present, ppl_present

    TYPE cexp_ppl_type
      REAL(dbl), DIMENSION(:), POINTER :: c
    END TYPE cexp_ppl_type

    TYPE(cexp_ppl_type), DIMENSION(:), ALLOCATABLE :: cexp_ppl

    TYPE vppl_type
      REAL(dbl), DIMENSION(:), POINTER   :: r2
      INTEGER, DIMENSION(:), POINTER    :: neighbor
      REAL(dbl), DIMENSION(:,:), POINTER :: r
    END TYPE vppl_type

    TYPE(vppl_type), DIMENSION(:), ALLOCATABLE :: vppl

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    part=particle_set,&
                    sac_ppl=sac_ppl) 

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    zetb(:)=0._dbl
    rpgfb(:)=0._dbl
    energy=0._dbl

    IF (calculate_forces) THEN
      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",natom*int_size)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set,&
                             natom=natom_of_kind)
        IF (ASSOCIATED(orb_basis_set)) THEN
          IF (.NOT.ASSOCIATED(force(ikind)%f_gth_ppl)) THEN
            ALLOCATE (force(ikind)%f_gth_ppl(3,natom_of_kind),STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"force(ikind)%f_gth_ppl",&
                               3*natom_of_kind*dbl)
            END IF
          END IF
          force(ikind)%f_gth_ppl(:,:) = 0.0_dbl
        END IF
      END DO
      nder = 1
    END IF

    maxder = ncoset(nder)
    maxder = MAX(1,nder)


!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             maxco=maxco,&
                             maxlgto=maxlgto,&
                             maxlppl=maxlppl,&
                             maxsgf=maxsgf)

    maxl = MAX(maxlgto,maxlppl)

    CALL init_orbital_pointers(maxl+nder+1)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxder=maxder)      

    ldsab = MAX(maxco,maxsgf)

    ALLOCATE (hab(ldsab*maxder,1),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"hab",ldsab*maxder*dp_size)
    END IF
    hab(:,:) = 0.0_dbl

    ALLOCATE (work(ldsab*maxder,1),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"work",ldsab*ldsab*maxder*dp_size)
    END IF
    work(:,:) = 0.0_dbl

    IF (calculate_forces) THEN
      ALLOCATE (pab(maxco,1),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"pab",maxco*dp_size)
      END IF
      pab(:,:) = 0.0_dbl
    END IF

!   *** Load GTH pseudo potential data (local part -> PPL) ***

    ppl_present = ASSOCIATED(sac_ppl)

    IF (ppl_present) THEN

      ALLOCATE (alpha_ppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"alpha_ppl",nkind*dp_size)

      ALLOCATE (cexp_ppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"cexp_ppl",nkind*int_size)
      DO i=1,SIZE(cexp_ppl)
        NULLIFY(cexp_ppl(i)%c)
      END DO

      ALLOCATE (ppl_radius(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_radius",nkind*dp_size)

      DO ikind=1,nkind

        atomic_kind => atomic_kind_set(ikind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             gth_potential=gth_potential)

        IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             alpha_ppl=alpha_ppl(ikind),&
                             cexp_ppl=cexp_ppl(ikind)%c,&
                             ppl_radius=ppl_radius(ikind))
        END IF

      END DO

      ALLOCATE (vppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"vppl",nkind*int_size)
      DO jkind=1,SIZE(vppl)
        NULLIFY (vppl(jkind)%r2)
        NULLIFY (vppl(jkind)%neighbor)
        NULLIFY (vppl(jkind)%r)
      END DO

    END IF

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom_of_kind,&
                           atom_list=atom_list,&
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

      DO iatom = 1, natom_of_kind

        atom_a = atom_list(iatom)

        IF (ppl_present) THEN

!           *** Retrieve the data of the SAC_PPL neighbors of atom "iatom" ***
 
          DO kkind=1,nkind

            ac = ikind + (kkind - 1)*nkind
  
            IF (.NOT.ASSOCIATED(sac_ppl(ac)%neighbor_list_set)) CYCLE

            sac_ppl_neighbor_list =>&
            find_neighbor_list(neighbor_list_set=&
                                     sac_ppl(ac)%neighbor_list_set,&
                                     atom=atom_a)

            CALL get_neighbor_list(neighbor_list=sac_ppl_neighbor_list,&
                                       nnode=nnode)

            ALLOCATE (vppl(kkind)%r2(nnode),STAT=istat)
            IF (istat /= 0) THEN
               CALL stop_memory(routine,"vppl(kkind)%r2",nnode*dp_size)
            END IF

            ALLOCATE (vppl(kkind)%neighbor(nnode),STAT=istat)
            IF (istat /= 0) THEN
               CALL stop_memory(routine,"vppl(kkind)%neighbor",nnode*int_size)
            END IF

            ALLOCATE (vppl(kkind)%r(3,nnode),STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"vppl(kkind)%r",3*nnode*dp_size)
            END IF
  
            sac_ppl_neighbor_node => first_node(sac_ppl_neighbor_list)
  
            DO knode=1,nnode
              CALL get_neighbor_node(neighbor_node=sac_ppl_neighbor_node,&
                                     neighbor=vppl(kkind)%neighbor(knode),&
                                     r=vppl(kkind)%r(:,knode))
              vppl(kkind)%r2(knode) =&
                vppl(kkind)%r(1,knode)*vppl(kkind)%r(1,knode) +&
                vppl(kkind)%r(2,knode)*vppl(kkind)%r(2,knode) +&
                vppl(kkind)%r(3,knode)*vppl(kkind)%r(3,knode)
              sac_ppl_neighbor_node => next(sac_ppl_neighbor_node)
            END DO
  
          END DO

        END IF

        rab2 = 0._dbl
        dab = SQRT(rab2)

        DO iset=1,nseta

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)

          IF (calculate_forces) THEN

            DO i=1,nsgfa(iset)
              work(i,1)=1.d0                   
            ENDDO
         
  
            CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                      1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                      work(1,1),SIZE(work,1),&
                      0.0_dbl,pab(1,1),SIZE(pab,1))

          END IF


!      *** Loop over the GTH pseudo potential atoms (local part) ***

          IF (ppl_present) THEN

            DO kkind=1,nkind

              IF (.NOT.ASSOCIATED(vppl(kkind)%neighbor)) CYCLE

              DO knode=1,SIZE(vppl(kkind)%neighbor)

                katom = vppl(kkind)%neighbor(knode)
                rac(:) = vppl(kkind)%r(:,knode)
                rac2 = vppl(kkind)%r2(knode)
                dac = SQRT(rac2)
                IF (set_radius_a(iset) + ppl_radius(kkind) < dac) CYCLE

                hab(:,:) = 0.0_dbl
                rab(:) =0.0_dbl

                rbc(:) = rac(:) - rab(:)
                rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                dbc = SQRT(rbc2)

!             *** Calculate the GTH pseudo potential forces ***

                IF (calculate_forces) THEN

                  CALL overlap_ppl(&
                    la_max(iset),zeta(1:npgfa(iset),iset),&
                    rpgfa(1:npgfa(iset),iset),la_min(iset),&
                    0,zetb,&
                    rpgfb,0,&
                    cexp_ppl(kkind)%c,alpha_ppl(kkind),&
                    ppl_radius(kkind),&
                    rab,rab2,dab,&
                    rac,rac2,dac,&
                    rbc,rbc2,dbc,&
                    hab,nder,0,.FALSE.,&
                    pab,force_a)

!               *** The derivatives w.r.t. atomic center c are    ***
!               *** calculated using the translational invariance ***
!               *** of the first derivatives                      ***

                  atom_c = atom_of_kind(katom)

                  force(ikind)%f_gth_ppl(:,iatom) =&
                    force(ikind)%f_gth_ppl(:,iatom) + force_a(:)
                  force(kkind)%f_gth_ppl(:,atom_c) =&
                    force(kkind)%f_gth_ppl(:,atom_c) - force_a(:)

                ELSE

                  CALL overlap_ppl(&
                    la_max(iset),zeta(1:npgfa(iset),iset),&
                    rpgfa(1:npgfa(iset),iset),la_min(iset),&
                    0,zetb,rpgfb,0,&
                    cexp_ppl(kkind)%c,alpha_ppl(kkind),&
                    ppl_radius(kkind),&
                    rab,rab2,dab,&
                    rac,rac2,dac,&
                    rbc,rbc2,dbc,&
                    hab,0,0,.FALSE.)

                END IF

                CALL dgemm("T","N",nsgfa(iset),1,ncoa,&
                          1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                          hab(1,1),SIZE(hab,1),&
                          0.0_dbl,work(1,1),SIZE(work,1))

                DO i=1,nsgfa(iset)
                  energy = energy - work(i,1)
                END DO

              END DO

            END DO

          END IF

        END DO 

        IF (ppl_present) THEN

          DO kkind=1,nkind
            IF (ASSOCIATED(vppl(kkind)%r2)) THEN
              DEALLOCATE (vppl(kkind)%r2,STAT=istat)
              IF (istat /= 0) CALL stop_memory(routine,"vppl(kkind)%r2")
              DEALLOCATE (vppl(kkind)%neighbor,STAT=istat)
              IF (istat /= 0) CALL stop_memory(routine,"vppl(kkind)%neighbor")
              DEALLOCATE (vppl(kkind)%r,STAT=istat)
              IF (istat /= 0) CALL stop_memory(routine,"vppl(kkind)%r")
            END IF
          END DO

        END IF

      END DO ! iatom => atom A

    END DO ! ikind 

!   *** Release work storage ***

    IF (ppl_present) THEN

      DEALLOCATE (alpha_ppl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"alpha_ppl")

      DO ikind=1,nkind
        NULLIFY (cexp_ppl(ikind)%c)
      END DO
      DEALLOCATE (cexp_ppl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"cexp_ppl")

      DEALLOCATE (ppl_radius,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_radius")

      DO ikind=1,nkind
        IF (ASSOCIATED(vppl(ikind)%r2)) THEN
          DEALLOCATE (vppl(ikind)%r2,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"vppl(ikind)%r2")
          DEALLOCATE (vppl(ikind)%r,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"vppl(ikind)%r")
        END IF
      END DO
      DEALLOCATE (vppl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"vppl")

    END IF
 
    DEALLOCATE (hab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    IF (calculate_forces) THEN
      DEALLOCATE (atom_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")
      DEALLOCATE (pab,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pab")
    END IF

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_rho0_ppl
! *****************************************************************************

  SUBROUTINE calculate_drho_ppl(eigenvector,eigenforce,kg_env,globenv,energy,&
                               calculate_forces,force)

!   ***************************************************************************

    REAL(dbl), DIMENSION(:), POINTER            :: eigenvector
    REAL(dbl), DIMENSION(:), POINTER            :: eigenforce
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(kg_environment_type), POINTER        :: kg_env
    REAL(dbl), INTENT (out)                   :: energy
    LOGICAL, INTENT(IN)                       :: calculate_forces
    TYPE(kg_force_type), DIMENSION(:), POINTER, OPTIONAL            :: force

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_drho_pp"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: aux_basis_set

    TYPE(neighbor_list_type), POINTER               :: sac_ppl_neighbor_list
    TYPE(neighbor_node_type), POINTER               :: sac_ppl_neighbor_node

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sac_ppl_aux
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set

    INTEGER                           :: ldsab,maxder,maxsgf,nder,nkind
    INTEGER  :: ab,ac,atom_a,atom_c,handle,&
                i,iatom,ikind,ilist,inode,ipgf,iset,ishell,istat,jkind,&
                katom,kkind,knode,kneigbor,l,lc_max,lc_min,maxco,maxl,&
                maxlgto,maxlppl,natom,natom_of_kind,ncoa,&
                nlist,nneighbor,nnode,nseta,sgfa,offset

    INTEGER, DIMENSION(:), ALLOCATABLE    :: atom_of_kind,kind_of
    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,npgfa,nsgfa 
    INTEGER, DIMENSION(:), POINTER                :: atom_list
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa

    REAL(dbl) :: dab,dac,dbc,rab2,rac2,rbc2
    REAL(dbl), DIMENSION(1) :: zetc,zetb,rpgfb
    REAL(dbl), DIMENSION(:,:), ALLOCATABLE :: hab,pab,work

    REAL(dbl), DIMENSION(3)            :: force_a
    REAL(dbl), DIMENSION(:), POINTER   :: set_radius_a
    REAL(dbl), DIMENSION(:,:), POINTER :: rpgfa,sphi_a,zeta

    REAL(dbl), DIMENSION(3)                  :: rab,rac,rbc
    REAL(dbl), DIMENSION(:), ALLOCATABLE     :: alpha_ppl,ppl_radius

    LOGICAL  :: gth_potential_present, ppl_present

    TYPE cexp_ppl_type
      REAL(dbl), DIMENSION(:), POINTER :: c
    END TYPE cexp_ppl_type

    TYPE(cexp_ppl_type), DIMENSION(:), ALLOCATABLE :: cexp_ppl

    TYPE vppl_type
      REAL(dbl), DIMENSION(:), POINTER   :: r2
      INTEGER, DIMENSION(:), POINTER    :: neighbor
      REAL(dbl), DIMENSION(:,:), POINTER :: r
    END TYPE vppl_type

    TYPE(vppl_type), DIMENSION(:), ALLOCATABLE :: vppl

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    part=particle_set,&
                    sac_ppl_aux=sac_ppl_aux) 

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    zetb(:)=0._dbl
    rpgfb(:)=0._dbl
    energy=0._dbl

    IF (calculate_forces) THEN
      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",natom*int_size)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             aux_basis_set=aux_basis_set,&
                             natom=natom_of_kind)
        IF (ASSOCIATED(aux_basis_set)) THEN
          IF (.NOT.ASSOCIATED(force(ikind)%f_gth_ppl)) THEN
            ALLOCATE (force(ikind)%f_gth_ppl(3,natom_of_kind),STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"force(ikind)%f_gth_ppl",&
                               3*natom_of_kind*dbl)
            END IF
          END IF
        END IF
      END DO
      nder = 1
    END IF

    maxder = ncoset(nder)
    maxder = MAX(1,nder)


!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             maxco=maxco,&
                             maxlgto=maxlgto,&
                             maxlppl=maxlppl,&
                             maxsgf=maxsgf)

    maxl = MAX(maxlgto,maxlppl)

    CALL init_orbital_pointers(maxl+nder+1)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxder=maxder)      

    ldsab = MAX(maxco,maxsgf)

    ALLOCATE (hab(ldsab*maxder,1),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"hab",ldsab*maxder*dp_size)
    END IF
    hab(:,:) = 0.0_dbl

    ALLOCATE (work(ldsab*maxder,1),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"work",ldsab*ldsab*maxder*dp_size)
    END IF
    work(:,:) = 0.0_dbl

    IF (calculate_forces) THEN
      ALLOCATE (pab(maxco,1),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"pab",maxco*dp_size)
      END IF
      pab(:,:) = 0.0_dbl
    END IF

!   *** Load GTH pseudo potential data (local part -> PPL) ***

    ppl_present = ASSOCIATED(sac_ppl_aux)

    IF (ppl_present) THEN

      ALLOCATE (alpha_ppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"alpha_ppl",nkind*dp_size)

      ALLOCATE (cexp_ppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"cexp_ppl",nkind*int_size)
      DO i=1,SIZE(cexp_ppl)
        NULLIFY(cexp_ppl(i)%c)
      END DO

      ALLOCATE (ppl_radius(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_radius",nkind*dp_size)

      offset=0

      DO ikind=1,nkind

        atomic_kind => atomic_kind_set(ikind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             gth_potential=gth_potential)

        IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             alpha_ppl=alpha_ppl(ikind),&
                             cexp_ppl=cexp_ppl(ikind)%c,&
                             ppl_radius=ppl_radius(ikind))
        END IF

      END DO

      ALLOCATE (vppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"vppl",nkind*int_size)
      DO jkind=1,SIZE(vppl)
        NULLIFY (vppl(jkind)%r2)
        NULLIFY (vppl(jkind)%neighbor)
        NULLIFY (vppl(jkind)%r)
      END DO

    END IF

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=natom_of_kind,&
                           atom_list=atom_list,&
                           aux_basis_set=aux_basis_set)

      IF (.NOT.ASSOCIATED(aux_basis_set)) CYCLE

      CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
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

      DO iatom = 1, natom_of_kind

        atom_a = atom_list(iatom)

        IF (ppl_present) THEN

!           *** Retrieve the data of the SAC_PPL neighbors of atom "iatom" ***
 
          DO kkind=1,nkind

            ac = ikind + (kkind - 1)*nkind
  
            IF (.NOT.ASSOCIATED(sac_ppl_aux(ac)%neighbor_list_set)) CYCLE

            sac_ppl_neighbor_list =>&
            find_neighbor_list(neighbor_list_set=&
                                     sac_ppl_aux(ac)%neighbor_list_set,&
                                     atom=atom_a)

            CALL get_neighbor_list(neighbor_list=sac_ppl_neighbor_list,&
                                       nnode=nnode)

            ALLOCATE (vppl(kkind)%r2(nnode),STAT=istat)
            IF (istat /= 0) THEN
               CALL stop_memory(routine,"vppl(kkind)%r2",nnode*dp_size)
            END IF

            ALLOCATE (vppl(kkind)%neighbor(nnode),STAT=istat)
            IF (istat /= 0) THEN
               CALL stop_memory(routine,"vppl(kkind)%neighbor",nnode*int_size)
            END IF

            ALLOCATE (vppl(kkind)%r(3,nnode),STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"vppl(kkind)%r",3*nnode*dp_size)
            END IF
  
            sac_ppl_neighbor_node => first_node(sac_ppl_neighbor_list)
  
            DO knode=1,nnode
              CALL get_neighbor_node(neighbor_node=sac_ppl_neighbor_node,&
                                     neighbor=vppl(kkind)%neighbor(knode),&
                                     r=vppl(kkind)%r(:,knode))
              vppl(kkind)%r2(knode) =&
                vppl(kkind)%r(1,knode)*vppl(kkind)%r(1,knode) +&
                vppl(kkind)%r(2,knode)*vppl(kkind)%r(2,knode) +&
                vppl(kkind)%r(3,knode)*vppl(kkind)%r(3,knode)
              sac_ppl_neighbor_node => next(sac_ppl_neighbor_node)
            END DO
  
          END DO

        END IF

        rab2 = 0._dbl
        dab = SQRT(rab2)

        DO iset=1,nseta

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)

          IF (calculate_forces) THEN

            DO i=1,nsgfa(iset)
              work(i,1)=eigenvector(offset+i)  
            ENDDO
         
  
            CALL dgemm("N","N",ncoa,1,nsgfa(iset),&
                      1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                      work(1,1),SIZE(work,1),&
                      0.0_dbl,pab(1,1),SIZE(pab,1))

          END IF


!      *** Loop over the GTH pseudo potential atoms (local part) ***

          IF (ppl_present) THEN

            DO kkind=1,nkind

              IF (.NOT.ASSOCIATED(vppl(kkind)%neighbor)) CYCLE

              DO knode=1,SIZE(vppl(kkind)%neighbor)

                katom = vppl(kkind)%neighbor(knode)
                rac(:) = vppl(kkind)%r(:,knode)
                rac2 = vppl(kkind)%r2(knode)
                dac = SQRT(rac2)
                IF (set_radius_a(iset) + ppl_radius(kkind) < dac) CYCLE

                hab(:,:) = 0.0_dbl
                rab(:) =0.0_dbl

                rbc(:) = rac(:) - rab(:)
                rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                dbc = SQRT(rbc2)

!             *** Calculate the GTH pseudo potential forces ***

                IF (calculate_forces) THEN

                  CALL overlap_ppl(&
                    la_max(iset),zeta(1:npgfa(iset),iset),&
                    rpgfa(1:npgfa(iset),iset),la_min(iset),&
                    0,zetb,&
                    rpgfb,0,&
                    cexp_ppl(kkind)%c,alpha_ppl(kkind),&
                    ppl_radius(kkind),&
                    rab,rab2,dab,&
                    rac,rac2,dac,&
                    rbc,rbc2,dbc,&
                    hab,nder,0,.FALSE.,&
                    pab,force_a)

!               *** The derivatives w.r.t. atomic center c are    ***
!               *** calculated using the translational invariance ***
!               *** of the first derivatives                      ***

                  atom_c = atom_of_kind(katom)

                  force(ikind)%f_gth_ppl(:,iatom) =&
                    force(ikind)%f_gth_ppl(:,iatom) - force_a(:)
                  force(kkind)%f_gth_ppl(:,atom_c) =&
                    force(kkind)%f_gth_ppl(:,atom_c) + force_a(:)

                ELSE

                  CALL overlap_ppl(&
                    la_max(iset),zeta(1:npgfa(iset),iset),&
                    rpgfa(1:npgfa(iset),iset),la_min(iset),&
                    0,zetb,rpgfb,0,&
                    cexp_ppl(kkind)%c,alpha_ppl(kkind),&
                    ppl_radius(kkind),&
                    rab,rab2,dab,&
                    rac,rac2,dac,&
                    rbc,rbc2,dbc,&
                    hab,0,0,.FALSE.)

                END IF

                CALL dgemm("T","N",nsgfa(iset),1,ncoa,&
                          1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                          hab(1,1),SIZE(hab,1),&
                          0.0_dbl,work(1,1),SIZE(work,1))

                DO i=1,nsgfa(iset)
                  energy = energy + work(i,1)*eigenvector(offset+i)
                  eigenforce(i+offset) = eigenforce(i+offset) - work(i,1)
                END DO


              END DO

            END DO

          END IF

          offset=offset+nsgfa(iset)

        END DO 

        IF (ppl_present) THEN

          DO kkind=1,nkind
            IF (ASSOCIATED(vppl(kkind)%r2)) THEN
              DEALLOCATE (vppl(kkind)%r2,STAT=istat)
              IF (istat /= 0) CALL stop_memory(routine,"vppl(kkind)%r2")
              DEALLOCATE (vppl(kkind)%neighbor,STAT=istat)
              IF (istat /= 0) CALL stop_memory(routine,"vppl(kkind)%neighbor")
              DEALLOCATE (vppl(kkind)%r,STAT=istat)
              IF (istat /= 0) CALL stop_memory(routine,"vppl(kkind)%r")
            END IF
          END DO

        END IF

      END DO ! iatom => atom A

    END DO ! ikind 

!   *** Release work storage ***

    IF (ppl_present) THEN

      DEALLOCATE (alpha_ppl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"alpha_ppl")

      DO ikind=1,nkind
        NULLIFY (cexp_ppl(ikind)%c)
      END DO
      DEALLOCATE (cexp_ppl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"cexp_ppl")

      DEALLOCATE (ppl_radius,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_radius")

      DO ikind=1,nkind
        IF (ASSOCIATED(vppl(ikind)%r2)) THEN
          DEALLOCATE (vppl(ikind)%r2,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"vppl(ikind)%r2")
          DEALLOCATE (vppl(ikind)%r,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"vppl(ikind)%r")
        END IF
      END DO
      DEALLOCATE (vppl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"vppl")

    END IF
 
    DEALLOCATE (hab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    IF (calculate_forces) THEN
      DEALLOCATE (atom_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")
      DEALLOCATE (pab,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pab")
    END IF

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_drho_ppl

! *****************************************************************************

END MODULE kg_ppl
