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
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_distributed,&
                                             dyn_coeff_set_type,&
                                             dyn_coeff_type
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: global_environment_type
  USE kg_energy_types,                 ONLY: kg_energy_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kg_force_types,                  ONLY: kg_force_type
  USE kinds,                           ONLY: dbl,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
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

  SUBROUTINE calculate_rho0_ppl(kg_env,globenv,calculate_forces)

!   ***************************************************************************

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    LOGICAL, INTENT(IN)                      :: calculate_forces

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_rho0_pp", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER :: ac, atom_a, atom_c, group, handle, i, ialpha, iatom, ikind, &
      iparticle_local, iset, istat, jkind, katom, kkind, knode, ldai, ldsab, &
      maxco, maxder, maxl, maxlgto, maxlppl, maxsgf, natom, natom_of_kind, &
      ncexp, ncoa, nder, nkind, nnode, nparticle_local, nseta, output_unit, &
      sgfa
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: atom_list, la_max, la_min, &
                                                npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: elp_potential_present, &
                                                gth_potential_present, &
                                                ionode, ppl_present
    REAL(dbl)                                :: alfa, dab, dac, dbc, e_ppl, &
                                                rab2, rac2, rbc2
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: ppl_radius
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: hab, pab, work
    REAL(dbl), DIMENSION(1)                  :: rpgfb, zetb
    REAL(dbl), DIMENSION(3)                  :: force_a, rab, rac, rbc
    REAL(dbl), DIMENSION(:), POINTER         :: cexp, set_radius_a
    REAL(dbl), DIMENSION(:, :), POINTER      :: rpgfa, sphi_a, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sac_ppl
    TYPE(neighbor_list_type), POINTER        :: sac_ppl_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sac_ppl_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    TYPE cexp_ppl_type
      REAL(dbl), DIMENSION(:,:), POINTER :: c
    END TYPE cexp_ppl_type

    TYPE alpha_ppl_type
      REAL(dbl), DIMENSION(:), POINTER :: alpha 
    END TYPE alpha_ppl_type

    TYPE(cexp_ppl_type), DIMENSION(:), ALLOCATABLE :: cexp_ppl
    TYPE(alpha_ppl_type), DIMENSION(:), ALLOCATABLE :: alpha_ppl

    TYPE vppl_type
      REAL(dbl), DIMENSION(:), POINTER   :: r2
      INTEGER, DIMENSION(:), POINTER    :: neighbor
      REAL(dbl), DIMENSION(:,:), POINTER :: r
    END TYPE vppl_type

    TYPE(vppl_type), DIMENSION(:), ALLOCATABLE :: vppl

    REAL(dbl), DIMENSION(:,:,:,:), ALLOCATABLE :: ai_work

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)
    
    NULLIFY(atomic_kind_set,atomic_kind,cexp,elp_potential,gth_potential)
    NULLIFY(orb_basis_set,particle_set,sac_ppl,set_radius_a,rpgfa,sphi_a,zeta)
    NULLIFY(sac_ppl_neighbor_list,sac_ppl_neighbor_node,atom_list,first_sgfa)
    NULLIFY(la_max,la_min,npgfa,nsgfa,energy,force,local_particles)

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    local_particles=local_particles,&
                    sac_ppl=sac_ppl,&
                    energy=energy,&
                    force=force) 

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    zetb(:)=0._dbl
    rpgfb(:)=0._dbl
    e_ppl=0._dbl

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
          IF (.NOT.ASSOCIATED(force(ikind)%f_ppl)) THEN
            ALLOCATE (force(ikind)%f_ppl(3,natom_of_kind),STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"force(ikind)%f_ppl",&
                               3*natom_of_kind*dbl)
            END IF
          END IF
          force(ikind)%f_ppl(:,:) = 0.0_dbl
        END IF
      END DO
      nder = 1
    END IF

    maxder = ncoset(nder)
    maxder = MAX(1,nder)


!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             elp_potential_present=elp_potential_present,&
                             maxco=maxco,&
                             maxlgto=maxlgto,&
                             maxlppl=maxlppl,&
                             maxsgf=maxsgf)

    maxl = MAX(maxlgto,maxlppl)

    CALL init_orbital_pointers(maxl+nder+1)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxder=maxder)      

    ldai = ncoset(maxl+nder+1)
    ALLOCATE (ai_work(ldai,ldai,ncoset(maxlppl),ncoset(nder+1)),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "ai_work",ldai*ldai*ncoset(maxlppl)*ncoset(nder+1)*dp_size)
    END IF

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

!   *** Load pseudo potential data (local part -> PPL) ***

    ppl_present = ASSOCIATED(sac_ppl)

    IF (ppl_present) THEN

      ALLOCATE (alpha_ppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"alpha_ppl",nkind*dp_size)
      DO i=1,SIZE(alpha_ppl)
        NULLIFY(alpha_ppl(i)%alpha)
      END DO

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
                             elp_potential=elp_potential,&
                             gth_potential=gth_potential)

        IF (ASSOCIATED(elp_potential)) THEN
          CALL get_potential(potential=elp_potential,&
                             alpha_ppl=alpha_ppl(ikind)%alpha,&
                             cexp_ppl=cexp_ppl(ikind)%c,&
                             ppl_radius=ppl_radius(ikind))
        ELSE IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             alpha_ppl=alfa,& 
                             cexp_ppl=cexp,&
                             ppl_radius=ppl_radius(ikind))
          ALLOCATE(cexp_ppl(ikind)%c(1,SIZE(cexp)),STAT=istat) 
          IF (istat /= 0) CALL stop_memory(routine,"cexp_ppl%c",&
                          SIZE(cexp)*int_size)
          DO i=1,SIZE(cexp)
            cexp_ppl(ikind)%c(1,i)=cexp(i)
          END DO 
          ALLOCATE(alpha_ppl(ikind)%alpha(1),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"alpha_ppl%alpha",&
                          1*int_size)
          alpha_ppl(ikind)%alpha(1)=alfa    
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

      nparticle_local = local_particles%n_el(ikind)

      DO iparticle_local=1,nparticle_local
        iatom = local_particles%list(ikind)%array(iparticle_local)

        IF (ppl_present) THEN

!           *** Retrieve the data of the SAC_PPL neighbors of atom "iatom" ***
 
          DO kkind=1,nkind

            ac = ikind + (kkind - 1)*nkind
  
            IF (.NOT.ASSOCIATED(sac_ppl(ac)%neighbor_list_set)) CYCLE

            sac_ppl_neighbor_list =>&
            find_neighbor_list(neighbor_list_set=&
                                     sac_ppl(ac)%neighbor_list_set,&
                                     atom=iatom)

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

                DO ialpha=1,size(alpha_ppl(kkind)%alpha)
                ncexp=size(cexp_ppl(kkind)%c,2)

                IF (calculate_forces) THEN

                  CALL overlap_ppl(&
                    la_max(iset),la_min(iset),npgfa(iset),&
                    rpgfa(:,iset),zeta(:,iset),&
                    0,0,1,rpgfb,zetb,&
                    cexp_ppl(kkind)%c(ialpha,1:ncexp),&
                    alpha_ppl(kkind)%alpha(ialpha),&
                    ppl_radius(kkind),&
                    rab,dab,rac,dac,rbc,dbc,&
                    hab,nder,0,.FALSE.,ai_work,&
                    pab,force_a)


!               *** The derivatives w.r.t. atomic center c are    ***
!               *** calculated using the translational invariance ***
!               *** of the first derivatives                      ***

                  atom_c = atom_of_kind(katom)
                  atom_a = atom_of_kind(iatom)

                  force(ikind)%f_ppl(:,atom_a) =&
                    force(ikind)%f_ppl(:,atom_a) + force_a(:)
                  force(kkind)%f_ppl(:,atom_c) =&
                    force(kkind)%f_ppl(:,atom_c) - force_a(:)

                ELSE

                  CALL overlap_ppl(&
                    la_max(iset),la_min(iset),npgfa(iset),&
                    rpgfa(:,iset),zeta(:,iset),&
                    0,0,1,rpgfb,zetb,&
                    cexp_ppl(kkind)%c(ialpha,1:ncexp),&
                    alpha_ppl(kkind)%alpha(ialpha),&
                    ppl_radius(kkind),&
                    rab,dab,rac,dac,rbc,dbc,&
                    hab,0,0,.FALSE.,ai_work)

                END IF
                END DO 

                CALL dgemm("T","N",nsgfa(iset),1,ncoa,&
                          1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                          hab(1,1),SIZE(hab,1),&
                          0.0_dbl,work(1,1),SIZE(work,1))

                DO i=1,nsgfa(iset)
                  e_ppl = e_ppl - work(i,1)
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

      DO ikind=1,nkind
        NULLIFY (alpha_ppl(ikind)%alpha)
      END DO
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

    DEALLOCATE (ai_work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "ai_work")
 
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

    CALL mp_sum(e_ppl,group)

    energy%pseudo = e_ppl

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_rho0_ppl
! *****************************************************************************

  SUBROUTINE calculate_drho_ppl(kg_env,globenv,calculate_forces)

!   ***************************************************************************

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    LOGICAL, INTENT(IN)                      :: calculate_forces

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_drho_pp", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER :: ac, atom_a, atom_c, group, handle, i, ialpha, iatom, ikind, &
      iparticle_local, iset, istat, jkind, katom, kkind, knode, ldai, ldsab, &
      maxco, maxder, maxl, maxlgto, maxlppl, maxsgf, natom, natom_of_kind, &
      ncexp, ncoa, nder, nkind, nnode, nparticle_local, nseta, offset, &
      output_unit, sgfa
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: atom_list, la_max, la_min, &
                                                npgfa, nsgfa
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: elp_potential_present, &
                                                gth_potential_present, &
                                                ionode, ppl_present
    REAL(dbl)                                :: alfa, coef, dab, dac, dbc, &
                                                e_ppl, rab2, rac2, rbc2
    REAL(dbl), ALLOCATABLE, DIMENSION(:)     :: ppl_radius
    REAL(dbl), ALLOCATABLE, DIMENSION(:, :)  :: hab, pab, work
    REAL(dbl), DIMENSION(1)                  :: rpgfb, zetb
    REAL(dbl), DIMENSION(3)                  :: force_a, rab, rac, rbc
    REAL(dbl), DIMENSION(:), POINTER         :: cexp, set_radius_a
    REAL(dbl), DIMENSION(:, :), POINTER      :: rpgfa, sphi_a, zeta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: local_coeffs
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sac_ppl_aux
    TYPE(neighbor_list_type), POINTER        :: sac_ppl_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sac_ppl_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    TYPE cexp_ppl_type
      REAL(dbl), DIMENSION(:,:), POINTER :: c
    END TYPE cexp_ppl_type

    TYPE alpha_ppl_type
      REAL(dbl), DIMENSION(:), POINTER :: alpha 
    END TYPE alpha_ppl_type

    TYPE(cexp_ppl_type), DIMENSION(:), ALLOCATABLE :: cexp_ppl
    TYPE(alpha_ppl_type), DIMENSION(:), ALLOCATABLE :: alpha_ppl

    TYPE vppl_type
      REAL(dbl), DIMENSION(:), POINTER   :: r2
      INTEGER, DIMENSION(:), POINTER    :: neighbor
      REAL(dbl), DIMENSION(:,:), POINTER :: r
    END TYPE vppl_type

    TYPE(vppl_type), DIMENSION(:), ALLOCATABLE :: vppl

    REAL(dbl), DIMENSION(:,:,:,:), ALLOCATABLE :: ai_work

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)
    
    NULLIFY(atomic_kind_set,atomic_kind,cexp,elp_potential,gth_potential)
    NULLIFY(aux_basis_set,particle_set,sac_ppl_aux,set_radius_a,rpgfa,sphi_a,zeta)
    NULLIFY(sac_ppl_neighbor_list,sac_ppl_neighbor_node,atom_list,first_sgfa)
    NULLIFY(la_max,la_min,npgfa,nsgfa,energy,force,local_particles)
    NULLIFY(dyn_coeff_set,local_coeffs)

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dyn_coeff_set=dyn_coeff_set,&
                    particle_set=particle_set,&
                    local_particles=local_particles,&
                    sac_ppl_aux=sac_ppl_aux,&
                    energy=energy,&
                    force=force) 


    IF(dyn_coeff_set%distribution_method/=dyn_coeff_distributed) THEN
      CALL stop_program ( 'kg_ppl','replicated coefs not yet implemented')
    END IF

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    zetb(:)=0._dbl
    rpgfb(:)=0._dbl
    e_ppl=0._dbl

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
          IF (.NOT.ASSOCIATED(force(ikind)%f_ppl)) THEN
            ALLOCATE (force(ikind)%f_ppl(3,natom_of_kind),STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"force(ikind)%f_ppl",&
                               3*natom_of_kind*dbl)
            END IF
            force(ikind)%f_ppl(:,:) = 0.0_dbl
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
                             elp_potential_present=elp_potential_present,&
                             maxco=maxco,&
                             maxlgto=maxlgto,&
                             maxlppl=maxlppl,&
                             maxsgf=maxsgf)

    maxl = MAX(maxlgto,maxlppl)

    CALL init_orbital_pointers(maxl+nder+1)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxder=maxder)      

    ldai = ncoset(maxl+nder+1)
    ALLOCATE (ai_work(ldai,ldai,ncoset(maxlppl),ncoset(nder+1)),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "ai_work",ldai*ldai*ncoset(maxlppl)*ncoset(nder+1)*dp_size)
    END IF

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

!   *** Load pseudo potential data (local part -> PPL) ***

    ppl_present = ASSOCIATED(sac_ppl_aux)

    IF (ppl_present) THEN

      ALLOCATE (alpha_ppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"alpha_ppl",nkind*dp_size)
      DO i=1,SIZE(alpha_ppl)
        NULLIFY(alpha_ppl(i)%alpha)
      END DO

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
                             elp_potential=elp_potential,&
                             gth_potential=gth_potential)

        IF (ASSOCIATED(elp_potential)) THEN
          CALL get_potential(potential=elp_potential,&
                             alpha_ppl=alpha_ppl(ikind)%alpha,&
                             cexp_ppl=cexp_ppl(ikind)%c,&
                             ppl_radius=ppl_radius(ikind))
        ELSE IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             alpha_ppl=alfa,& 
                             cexp_ppl=cexp,&
                             ppl_radius=ppl_radius(ikind))
          ALLOCATE(cexp_ppl(ikind)%c(1,SIZE(cexp)),STAT=istat) 
          IF (istat /= 0) CALL stop_memory(routine,"cexp_ppl%c",&
                          SIZE(cexp)*int_size)
          DO i=1,SIZE(cexp)
            cexp_ppl(ikind)%c(1,i)=cexp(i)
          END DO 
          ALLOCATE(alpha_ppl(ikind)%alpha(1),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"alpha_ppl%alpha",&
                          1*int_size)
          alpha_ppl(ikind)%alpha(1)=alfa    
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

      local_coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
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

      nparticle_local = local_particles%n_el(ikind)

      DO iparticle_local=1,nparticle_local
        iatom = local_particles%list(ikind)%array(iparticle_local)

        IF (ppl_present) THEN

!           *** Retrieve the data of the SAC_PPL neighbors of atom "iatom" ***
 
          DO kkind=1,nkind

            ac = ikind + (kkind - 1)*nkind
  
            IF (.NOT.ASSOCIATED(sac_ppl_aux(ac)%neighbor_list_set)) CYCLE

            sac_ppl_neighbor_list =>&
            find_neighbor_list(neighbor_list_set=&
                                     sac_ppl_aux(ac)%neighbor_list_set,&
                                     atom=iatom)

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
        offset=0

        DO iset=1,nseta

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)

          IF (calculate_forces) THEN

            DO i=1,nsgfa(iset)
              work(i,1)=local_coeffs%pos(iparticle_local,i+offset)
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

                DO ialpha=1,size(alpha_ppl(kkind)%alpha)
                ncexp=size(cexp_ppl(kkind)%c,2)

                  IF (calculate_forces) THEN

                    CALL overlap_ppl(&
                      la_max(iset),la_min(iset),npgfa(iset),&
                      rpgfa(:,iset),zeta(:,iset),&
                      0,0,1,rpgfb,zetb,&
                      cexp_ppl(kkind)%c(ialpha,1:ncexp),&
                      alpha_ppl(kkind)%alpha(ialpha),&
                      ppl_radius(kkind),&
                      rab,dab,rac,dac,rbc,dbc,&
                      hab,nder,0,.FALSE.,ai_work,&
                      pab,force_a)

!               *** The derivatives w.r.t. atomic center c are    ***
!               *** calculated using the translational invariance ***
!               *** of the first derivatives                      ***

                    atom_c = atom_of_kind(katom)
                    atom_a = atom_of_kind(iatom)

                    force(ikind)%f_ppl(:,atom_a) =&
                      force(ikind)%f_ppl(:,atom_a) + force_a(:)
                    force(kkind)%f_ppl(:,atom_c) =&
                      force(kkind)%f_ppl(:,atom_c) - force_a(:)

                  ELSE

                    CALL overlap_ppl(&
                      la_max(iset),la_min(iset),npgfa(iset),&
                      rpgfa(:,iset),zeta(:,iset),&
                      0,0,1,rpgfb,zetb,&
                      cexp_ppl(kkind)%c(ialpha,1:ncexp),&
                      alpha_ppl(kkind)%alpha(ialpha),&
                      ppl_radius(kkind),&
                      rab,dab,rac,dac,rbc,dbc,&
                      hab,0,0,.FALSE.,ai_work)

                  END IF
                END DO 

                CALL dgemm("T","N",nsgfa(iset),1,ncoa,&
                          1.0_dbl,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                          hab(1,1),SIZE(hab,1),&
                          0.0_dbl,work(1,1),SIZE(work,1))


                DO i=1,nsgfa(iset)
                  coef=local_coeffs%pos(iparticle_local,i+offset)
                  e_ppl = e_ppl + work(i,1)*coef
                  local_coeffs%forces(iparticle_local,i+offset) = &
                      local_coeffs%forces(iparticle_local,i+offset) - work(i,1)
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

      DO ikind=1,nkind
        NULLIFY (alpha_ppl(ikind)%alpha)
      END DO
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
 
    DEALLOCATE (ai_work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "ai_work")
 
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

    CALL mp_sum(e_ppl,group)

    energy%ppseudo = e_ppl

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE calculate_drho_ppl
!*****************************************************************************

END MODULE kg_ppl
