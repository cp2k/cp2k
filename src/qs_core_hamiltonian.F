!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_core_hamiltonian [1.1] *
!!
!!   FUNCTION
!!     Calculation of the core Hamiltonian integral matrix <a|H|b> over
!!     Cartesian Gaussian-type functions.
!!
!!     <a|H|b> = <a|T|b> + <a|V|b>
!!
!!
!!     Kinetic energy:
!!
!!     <a|T|b> = <a|-nabla**2/2|b>
!!               \_______________/
!!                       |
!!                    kinetic
!!
!!     Nuclear potential energy:
!!
!!     a) Allelectron calculation:
!!
!!                         erfc(r)
!!        <a|V|b> = -Z*<a|---------|b>
!!                            r
!!
!!                         1 - erf(r)
!!                = -Z*<a|------------|b>
!!                             r
!!
!!                          1           erf(r)
!!                = -Z*(<a|---|b> - <a|--------|b>)
!!                          r             r
!!
!!                          1
!!                = -Z*(<a|---|b> - N*<ab||c>)
!!                          r
!!
!!                     -Z
!!                = <a|---|b> + Z*N*<ab||c>
!!                      r
!!                  \_______/       \_____/
!!                      |              |
!!                   nuclear        coulomb
!!
!!     b) Pseudopotential calculation (Goedecker, Teter and Hutter; GTH):
!!
!!        <a|V|b> = <a|(V(local) + V(non-local))|b>
!!
!!                = <a|(V(local)|b> + <a|V(non-local))|b>
!!
!!        <a|V(local)|b> = <a|-Z(eff)*erf(SQRT(2)*alpha*r)/r +
!!                            (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!!                             C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>
!!
!!        <a|V(non-local)|b> = <a|p(l,i)>*h(i,j)*<p(l,j)|b>
!!
!!   AUTHOR
!!     Matthias Krack (14.09.2000,21.03.02)
!!
!!   LITERATURE
!!     S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996)
!!     C. Hartwigsen, S. Goedecker and J. Hutter, Phys. Rev. B 58, 3641 (1998)
!!     M. Krack and M. Parrinello, Phys. Chem. Chem. Phys. 2, 2105 (2000)
!!     S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_core_hamiltonian

! *****************************************************************************

  USE kinds, ONLY: int_size,&
                   wp => dp,&
                   wp_size => dp_size

  USE ai_kinetic,               ONLY: kinetic
  USE ai_overlap,               ONLY: overlap
  USE ai_overlap_ppl,           ONLY: overlap_ppl
  USE atomic_kind_types,        ONLY: atomic_kind_type,&
                                      get_atomic_kind,&
                                      get_atomic_kind_set
  USE basis_set_types,          ONLY: get_gto_basis_set,&
                                      gto_basis_set_type
  USE external_potential_types, ONLY: get_potential,&
                                      gth_potential_type
  USE global_types,             ONLY: global_environment_type
  USE mathlib,                  ONLY: symmetrize_matrix
  USE message_passing,          ONLY: mp_sum
  USE orbital_pointers,         ONLY: nco,ncoset,nso
  USE particle_types,           ONLY: particle_type
  USE qs_environment_types,     ONLY: get_qs_env,&
                                      qs_environment_type,&
                                      set_qs_env
  USE qs_neighbor_list_types,   ONLY: find_neighbor_list,&
                                      first_list,&
                                      first_node,&
                                      get_neighbor_list,&
                                      get_neighbor_node,&
                                      neighbor_list_set_type,&
                                      neighbor_list_set_p_type,&
                                      neighbor_list_type,&
                                      neighbor_node_type,&
                                      next
  USE qs_overlap,               ONLY: write_sparse_matrix
  USE sparse_matrix_types,      ONLY: add_block_node,&
                                      allocate_matrix,&
                                      copy_sparse_to_full_matrix,&
                                      deallocate_matrix,&
                                      deallocate_matrix_row,&
                                      first_block_node,&
                                      get_block_node,&
                                      get_matrix_info,&
                                      next_block_node,&
                                      real_block_node_type,&
                                      real_matrix_p_type,&
                                      real_matrix_type,&
                                      replicate_matrix_structure
  USE termination,              ONLY: stop_memory
  USE timings,                  ONLY: timeset,timestop

  IMPLICIT NONE

  PRIVATE

  TYPE(atomic_kind_type), POINTER   :: atomic_kind
  TYPE(gth_potential_type), POINTER :: gth_potential
  TYPE(gto_basis_set_type), POINTER :: orb_basis_set
  TYPE(real_matrix_p_type), POINTER :: h,s
  INTEGER                           :: nkind

  TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
  TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_ppl,&
                                                           sca_ppnl,scb_ppnl
  TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set

  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: hab,sab,work

  REAL(wp), DIMENSION(:), POINTER   :: set_radius_a,set_radius_b
  INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,lb_max,lb_min,npgfa,&
                                       npgfb,nsgfa,nsgfb
  REAL(wp), DIMENSION(:,:), POINTER :: rpgfa,rpgfb,sphi_a,sphi_b,zeta,zetb
  INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa,first_sgfb

! *** Public subroutines ***

  PUBLIC :: build_core_hamiltonian_matrix

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_core_hamiltonian_matrix(qs_env,globenv)

!   Purpose: Calculation of the overlap matrix over Cartesian Gaussian
!            functions.

!   History: - Creation (11.03.2002, MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE build_core_hamiltonian_matrix (MODULE qs_core_hamiltonian)"

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: sab_orb_neighbor_list,&
                                         sac_ppl_neighbor_list
    TYPE(neighbor_node_type), POINTER :: sab_orb_neighbor_node,&
                                         sac_ppl_neighbor_node
    TYPE(real_matrix_p_type), POINTER :: t

    REAL(wp) :: dab,dac,dbc,rab2,rac2,rbc2
    INTEGER  :: ab,ac,first_setb,handle,iatom,ikind,inode,ipgf,iset,ishell,&
                istat,jatom,jkind,jpgf,jset,jshell,katom,kkind,maxco,maxfun,&
                maxppnl,maxsgf,ncoa,ncob,nneighbor,nnode,nseta,nsetb,sgfa,sgfb
    LOGICAL  :: gth_potential_present

    REAL(wp), DIMENSION(3)                :: rab,rac,rbc
    REAL(wp), DIMENSION(:), ALLOCATABLE   :: alpha_ppl,ppl_radius
    INTEGER, DIMENSION(:), ALLOCATABLE    :: sac_ppl_neighbors
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sac_ppl_r

    REAL(wp), DIMENSION(:,:), POINTER :: h_block,s_block,t_block

    TYPE cexp_ppl_type
      REAL(wp), DIMENSION(:), POINTER :: c
    END TYPE cexp_ppl_type

    TYPE(cexp_ppl_type), DIMENSION(:), ALLOCATABLE :: cexp_ppl

    TYPE vppl_type
      REAL(wp), DIMENSION(:), POINTER   :: r2
      REAL(wp), DIMENSION(:,:), POINTER :: r
    END TYPE vppl_type

    TYPE(vppl_type), DIMENSION(:), ALLOCATABLE :: vppl

!   ---------------------------------------------------------------------------

    CALL timeset("build_core_hamiltonian","I","",handle)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    h=h,&
                    s=s,&
                    sab_orb=sab_orb,&
                    sac_ppl=sac_ppl,&
                    sca_ppnl=sca_ppnl,&
                    scb_ppnl=scb_ppnl)

!   *** Allocate the kinetic energy integral ***
!   *** matrix (only need for printing)      ***

    IF (globenv%print%kinetic_energy_matrix) THEN
      NULLIFY (t)
      CALL replicate_matrix_structure(source=s,&
                                      target=t,&
                                      target_name="KINETIC ENERGY MATRIX")
    END IF

!   *** Allocate the core Hamiltonian matrix ***

    NULLIFY (h)
    CALL replicate_matrix_structure(source=s,&
                                    target=h,&
                                    target_name="CORE HAMILTONIAN MATRIX")

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             maxco=maxco,&
                             maxppnl=maxppnl,&
                             maxsgf=maxsgf)

    maxfun = MAX(maxco,maxppnl)

    ALLOCATE (hab(maxco,maxco),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"hab",maxco*maxco*wp_size)

    ALLOCATE (sab(maxfun,maxco),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"sab",maxfun*maxco*wp_size)

    ALLOCATE (work(maxfun,maxsgf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxfun*maxsgf*wp_size)

    nkind = SIZE(atomic_kind_set)

    IF (gth_potential_present) THEN

      ALLOCATE (alpha_ppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"alpha_ppl",nkind*wp_size)

      ALLOCATE (cexp_ppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"cexp_ppl",nkind*int_size)

      ALLOCATE (ppl_radius(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ppl_radius",nkind*wp_size)

      DO ikind=1,nkind

        atomic_kind => atomic_kind_set(ikind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             gth_potential=gth_potential)

        CALL get_potential(potential=gth_potential,&
                           alpha_ppl=alpha_ppl(ikind),&
                           cexp_ppl=cexp_ppl(ikind)%c,&
                           ppl_radius=ppl_radius(ikind))

      END DO

      ALLOCATE (vppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"vppl",nkind*int_size)

    END IF

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
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

      DO jkind=ikind,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               first_sgf=first_sgfb,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               nsgf_set=nsgfb,&
                               pgf_radius=rpgfb,&
                               set_radius=set_radius_b,&
                               sphi=sphi_b,&
                               zet=zetb)

        ab = ikind + jkind*(jkind - 1)/2

        sab_orb_neighbor_list => first_list(sab_orb(ab)%neighbor_list_set)

        DO WHILE (ASSOCIATED(sab_orb_neighbor_list))

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                                 atom=iatom)

!         *** Retrieve the data of the SAC_PPL neighbors of atom "iatom" ***

          IF (gth_potential_present) THEN

            DO kkind=1,nkind

              ac = ikind + (kkind - 1)*nkind

              IF (.NOT.ASSOCIATED(sac_ppl(ac)%neighbor_list_set)) THEN
                NULLIFY (vppl(kkind)%r2)
                NULLIFY (vppl(kkind)%r)
                CYCLE
              END IF

              sac_ppl_neighbor_list =>&
                find_neighbor_list(neighbor_list_set=&
                                     sac_ppl(ac)%neighbor_list_set,&
                                   atom=iatom)

              CALL get_neighbor_list(neighbor_list=sac_ppl_neighbor_list,&
                                     nnode=nnode)

              ALLOCATE (vppl(kkind)%r2(nnode),STAT=istat)
              IF (istat /= 0) THEN
                CALL stop_memory(routine,"vppl(kkind)%r2",nnode*wp_size)
              END IF

              ALLOCATE (vppl(kkind)%r(3,nnode),STAT=istat)
              IF (istat /= 0) THEN
                CALL stop_memory(routine,"vppl(kkind)%r",3*nnode*wp_size)
              END IF

              inode = 0

              sac_ppl_neighbor_node => first_node(sac_ppl_neighbor_list)

              DO WHILE (ASSOCIATED(sac_ppl_neighbor_node))
                inode = inode + 1
                CALL get_neighbor_node(neighbor_node=sac_ppl_neighbor_node,&
                                       r=vppl(kkind)%r(:,inode))
                vppl(kkind)%r2(inode) =&
                  vppl(kkind)%r(1,inode)*vppl(kkind)%r(1,inode) +&
                  vppl(kkind)%r(2,inode)*vppl(kkind)%r(2,inode) +&
                  vppl(kkind)%r(3,inode)*vppl(kkind)%r(3,inode)
                sac_ppl_neighbor_node => next(sac_ppl_neighbor_node)
              END DO

            END DO

          END IF

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

          DO WHILE (ASSOCIATED(sab_orb_neighbor_node))

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab(:))

            IF (iatom <= jatom) THEN
              CALL get_block_node(matrix=s%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=s_block)
            ELSE
              CALL get_block_node(matrix=s%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  block=s_block)
            END IF

!           *** Check, if the atomic block has to be ***
!           *** calculated by the current processor  ***

            IF (ASSOCIATED(s_block)) THEN

              IF (iatom <= jatom) THEN
                CALL get_block_node(matrix=h%matrix,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    block=h_block)
                IF (globenv%print%kinetic_energy_matrix) THEN
                  CALL get_block_node(matrix=t%matrix,&
                                      block_row=iatom,&
                                      block_col=jatom,&
                                      block=t_block)
                END IF
              ELSE
                CALL get_block_node(matrix=h%matrix,&
                                    block_row=jatom,&
                                    block_col=iatom,&
                                    block=h_block)
                IF (globenv%print%kinetic_energy_matrix) THEN
                  CALL get_block_node(matrix=t%matrix,&
                                      block_row=jatom,&
                                      block_col=iatom,&
                                      block=t_block)
                END IF
              END IF

              rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
              dab = SQRT(rab2)

              DO iset=1,nseta

                ncoa = npgfa(iset)*ncoset(la_max(iset))
                sgfa = first_sgfa(1,iset)

                IF (iatom == jatom) THEN
                  first_setb = iset
                ELSE
                  first_setb = 1
                END IF

                DO jset=first_setb,nsetb

                  IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1,jset)

!                 *** Calculate the primitive overlap integrals ***

                  CALL overlap(la_max(iset),zeta(1:npgfa(iset),iset),&
                               rpgfa(1:npgfa(iset),iset),la_min(iset),&
                               lb_max(jset),zetb(1:npgfb(jset),jset),&
                               rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                               rab,rab2,sab)

!                 *** Contraction step (overlap matrix) ***

                  CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                             1.0_wp,sab(1,1),SIZE(sab,1),&
                             sphi_b(1,sgfb),SIZE(sphi_b,1),&
                             0.0_wp,work(1,1),SIZE(work,1))

                  IF (iatom <= jatom) THEN
                    CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                               1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               work(1,1),SIZE(work,1),&
                               1.0_wp,s_block(sgfa,sgfb),SIZE(s_block,1))
                  ELSE
                    CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                               1.0_wp,work(1,1),SIZE(work,1),&
                               sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               1.0_wp,s_block(sgfb,sgfa),SIZE(s_block,1))
                  END IF

!                 *** Calculate the primitive kinetic energy integrals ***

                  CALL kinetic(la_max(iset),zeta(1:npgfa(iset),iset),&
                               rpgfa(1:npgfa(iset),iset),la_min(iset),&
                               lb_max(jset),zetb(1:npgfb(jset),jset),&
                               rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                               rab,rab2,sab,hab)

!                 *** Contraction step (kinetic energy integral ***
!                 *** matrix (only needed for printing)         ***

                  IF (globenv%print%kinetic_energy_matrix) THEN
                    CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                               1.0_wp,hab(1,1),SIZE(hab,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_wp,work(1,1),SIZE(work,1))
                    IF (iatom <= jatom) THEN
                      CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 work(1,1),SIZE(work,1),&
                                 1.0_wp,t_block(sgfa,sgfb),SIZE(t_block,1))
                    ELSE
                      CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                 1.0_wp,work(1,1),SIZE(work,1),&
                                 sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 1.0_wp,t_block(sgfb,sgfa),SIZE(t_block,1))
                    END IF
                  END IF

!                 *** Loop over the GTH pseudo potential atoms (local part) ***

                  DO kkind=1,nkind

                    IF (.NOT.ASSOCIATED(vppl(kkind)%r2)) CYCLE

                    DO katom=1,SIZE(vppl(kkind)%r2)

                      rac(:) = vppl(kkind)%r(:,katom)
                      rac2 = vppl(kkind)%r2(katom)
                      dac = SQRT(rac2)
                      IF (set_radius_a(iset) + ppl_radius(kkind) < dac) CYCLE

                      rbc(:) = rac(:) - rab(:)
                      rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                      dbc = SQRT(rbc2)
                      IF (set_radius_b(jset) + ppl_radius(kkind) < dbc) CYCLE

                      CALL overlap_ppl(la_max(iset),zeta(1:npgfa(iset),iset),&
                                       rpgfa(1:npgfa(iset),iset),la_min(iset),&
                                       lb_max(jset),zetb(1:npgfb(jset),jset),&
                                       rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                                       cexp_ppl(kkind)%c,alpha_ppl(kkind),&
                                       ppl_radius(kkind),&
                                       rab,rab2,dab,&
                                       rac,rac2,dac,&
                                       rbc,rbc2,dbc,hab)

                    END DO

                  END DO

!                 *** Contraction step (core Hamiltonian matrix) ***

                  CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                             1.0_wp,hab(1,1),SIZE(hab,1),&
                             sphi_b(1,sgfb),SIZE(sphi_b,1),&
                             0.0_wp,work(1,1),SIZE(work,1))

                  IF (iatom <= jatom) THEN
                    CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                               1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               work(1,1),SIZE(work,1),&
                               1.0_wp,h_block(sgfa,sgfb),SIZE(h_block,1))
                  ELSE
                    CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                               1.0_wp,work(1,1),SIZE(work,1),&
                               sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               1.0_wp,h_block(sgfb,sgfa),SIZE(h_block,1))
                  END IF

                END DO

              END DO

!             *** Symmetrize the diagonal blocks ***

              IF (iatom == jatom) THEN
                CALL symmetrize_matrix(s_block,"upper_to_lower")
                IF (globenv%print%kinetic_energy_matrix) THEN
                  CALL symmetrize_matrix(t_block,"upper_to_lower")
                END IF
                CALL symmetrize_matrix(h_block,"upper_to_lower")
              END IF

            END IF

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO

          IF (gth_potential_present) THEN
            DO kkind=1,nkind
              IF (ASSOCIATED(vppl(kkind)%r2)) THEN
                DEALLOCATE (vppl(kkind)%r2,STAT=istat)
                IF (istat /= 0) CALL stop_memory(routine,"vppl(kkind)%r2")
                DEALLOCATE (vppl(kkind)%r,STAT=istat)
                IF (istat /= 0) CALL stop_memory(routine,"vppl(kkind)%r")
              END IF
            END DO
          END IF

          sab_orb_neighbor_list => next(sab_orb_neighbor_list)

        END DO

      END DO

    END DO

!   *** Add the V(non-local) contribution of the GTH pseudo potential ***

    IF (gth_potential_present) CALL add_vppnl()

!   *** Put the core Hamiltonian matrix in the QS environment ***

    CALL set_qs_env(qs_env=qs_env,h=h)

!   *** Release work storage ***

    IF (gth_potential_present) THEN

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

    DEALLOCATE (sab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"sab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    CALL timestop(0.0_wp,handle)

!   *** Print the overlap integral matrix, if requested ***

    IF (globenv%print%overlap_matrix) THEN
      CALL write_sparse_matrix(s%matrix,4,6,qs_env,globenv)
    END IF

!   *** Print the kinetic energy integral matrix, if requested ***

    IF (globenv%print%kinetic_energy_matrix) THEN
      CALL write_sparse_matrix(t%matrix,4,6,qs_env,globenv)
      CALL deallocate_matrix(t)
    END IF

!   *** Print the core Hamiltonian matrix, if requested ***

    IF (globenv%print%core_hamiltonian_matrix) THEN
      CALL write_sparse_matrix(h%matrix,4,6,qs_env,globenv)
    END IF

  END SUBROUTINE build_core_hamiltonian_matrix

! *****************************************************************************

  SUBROUTINE add_vppnl()

!   Purpose: Add the non-local contribution of the GTH pseudo potential to the
!            core Hamiltonian matrix.

!   History: - Creation (24.10.2000, Matthias Krack)

!   ***************************************************************************

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE add_vppnl (MODULE qs_core_hamiltonian)"

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER   :: sca_neighbor_list,scb_neighbor_list
    TYPE(neighbor_node_type), POINTER   :: sca_neighbor_node,scb_neighbor_node
    TYPE(real_block_node_type), POINTER :: block_node_ca,block_node_cb
    TYPE(real_matrix_type), POINTER     :: sca,scb

    REAL(wp) :: dca,rca2
    INTEGER  :: atom_c,current_kind,first_prjc,handle,iatom,ikind,iset,istat,&
                jatom,kkind,l,lc_max,lc_min,lppnl,na,natom,nb,ncoa,ncoc,nppnl,&
                nprjc,nseta,nsgf,prjc,sgfa
    LOGICAL  :: ppnl_present

    REAL(wp), DIMENSION(1) :: rprjc,zetc
    INTEGER, DIMENSION(1)  :: first_prj,last_prj
    REAL(wp), DIMENSION(3) :: rca
    INTEGER, DIMENSION(3)  :: cell_c

    REAL(wp), DIMENSION(:), POINTER    :: alpha_ppnl
    INTEGER, DIMENSION(:), ALLOCATABLE :: kind_of
    INTEGER, DIMENSION(:), POINTER     :: nprj_ppnl
    REAL(wp), DIMENSION(:,:), POINTER  :: cprj,h_block,sca_block,scb_block,&
                                          vprj_ppnl

!   ---------------------------------------------------------------------------

    CALL timeset("add_vppnl","I","",handle)

    NULLIFY (sca)
    NULLIFY (scb)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom,&
                             nsgf=nsgf)

!   *** Allocate and initialise work storage ***

    ALLOCATE (kind_of(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_of",natom*int_size)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             kind_of=kind_of)

!   *** Loop over all GTH-PPNL operator atom kinds ***

    DO kkind=1,nkind

      atomic_kind => atomic_kind_set(kkind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           gth_potential=gth_potential)

      IF (.NOT.ASSOCIATED(gth_potential)) CYCLE

      CALL get_potential(potential=gth_potential,&
                         ppnl_present=ppnl_present)

      IF (.NOT.ppnl_present) CYCLE

      CALL get_potential(potential=gth_potential,&
                         alpha_ppnl=alpha_ppnl,&
                         cprj=cprj,&
                         lppnl=lppnl,&
                         nppnl=nppnl,&
                         nprj_ppnl=nprj_ppnl,&
                         ppnl_radius=rprjc(1),&
                         vprj_ppnl=vprj_ppnl)

      first_prj(1) = 1
      last_prj(1) = nppnl

!     *** Allocate the CA projector (ppnl) overlap matrices row ***

      CALL allocate_matrix(matrix=sca,&
                           nblock_row=1,&
                           nblock_col=natom,&
                           nrow=nppnl,&
                           ncol=nsgf,&
                           first_row=first_prj(:),&
                           last_row=last_prj(:),&
                           first_col=particle_set(:)%first_sgf,&
                           last_col=particle_set(:)%last_sgf,&
                           matrix_name="CA PROJECTOR (PPNL) OVERLAP MATRIX",&
                           matrix_symmetry="no symmetry")

!     *** Allocate the CB projector (ppnl) overlap matrices row ***

      CALL allocate_matrix(matrix=scb,&
                           nblock_row=1,&
                           nblock_col=natom,&
                           nrow=nppnl,&
                           ncol=nsgf,&
                           first_row=first_prj(:),&
                           last_row=last_prj(:),&
                           first_col=particle_set(:)%first_sgf,&
                           last_col=particle_set(:)%last_sgf,&
                           matrix_name="CB PROJECTOR (PPNL) OVERLAP MATRIX",&
                           matrix_symmetry="no symmetry")

!     *** Calculate the overlap integrals between the orbital ***
!     *** basis functions and the GTH(nl) projector functions ***

      sca_neighbor_list => first_list(sca_ppnl(kkind)%neighbor_list_set)
      scb_neighbor_list => first_list(scb_ppnl(kkind)%neighbor_list_set)

      DO WHILE (ASSOCIATED(sca_neighbor_list))

        current_kind = 0

        sca_neighbor_node => first_node(sca_neighbor_list)

        DO WHILE (ASSOCIATED(sca_neighbor_node))

          CALL get_neighbor_node(neighbor_node=sca_neighbor_node,&
                                 neighbor=iatom,&
                                 r=rca(:))

          ikind = kind_of(iatom)

          IF (ikind /= current_kind) THEN

            atomic_kind => atomic_kind_set(ikind)

            CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                 orb_basis_set=orb_basis_set)

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

          END IF

          rca2 = rca(1)*rca(1) + rca(2)*rca(2) + rca(3)*rca(3)
          dca = SQRT(rca2)

!         *** Create a new matrix block ***

          NULLIFY (sca_block)

          CALL add_block_node(matrix=sca,&
                              block_row=1,&
                              block_col=iatom,&
                              block=sca_block)

!         *** Calculate the current overlap matrix block ***

          prjc = 1

          DO l=0,lppnl

            nprjc = nprj_ppnl(l)*nco(l)

            IF (nprjc == 0) CYCLE

            lc_max = l + 2*(nprj_ppnl(l) - 1)
            lc_min = l
            zetc(1) = alpha_ppnl(l)
            ncoc = ncoset(lc_max)

            DO iset=1,nseta

              IF (rprjc(1) + set_radius_a(iset) < dca) CYCLE

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

!             *** Calculate the primitive overlap integrals ***

              CALL overlap(lc_max,zetc,rprjc,lc_min,&
                           la_max(iset),zeta(1:npgfa(iset),iset),&
                           rpgfa(1:npgfa(iset),iset),la_min(iset),&
                           rca,rca2,sab)

!             *** Contraction step (orbital basis functions) ***

              CALL dgemm("N","N",ncoc,nsgfa(iset),ncoa,&
                         1.0_wp,sab(1,1),SIZE(sab,1),&
                         sphi_a(1,sgfa),SIZE(sphi_a,1),&
                         0.0_wp,work(1,1),SIZE(work,1))

!             *** Contraction step (projector functions) ***

              CALL dgemm("T","N",nprjc,nsgfa(iset),ncoc,&
                         1.0_wp,cprj(1,prjc),SIZE(cprj,1),&
                         work(1,1),SIZE(work,1),&
                         1.0_wp,sca_block(prjc,sgfa),SIZE(sca_block,1))

            END DO

            prjc = prjc + nprjc

          END DO

          sca_neighbor_node => next(sca_neighbor_node)

        END DO

        current_kind = 0

        scb_neighbor_node => first_node(scb_neighbor_list)

        DO WHILE (ASSOCIATED(scb_neighbor_node))

          CALL get_neighbor_node(neighbor_node=scb_neighbor_node,&
                                 neighbor=iatom,&
                                 r=rca(:))

          ikind = kind_of(iatom)

          IF (ikind /= current_kind) THEN

            atomic_kind => atomic_kind_set(ikind)

            CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                 orb_basis_set=orb_basis_set)

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

          END IF

          rca2 = rca(1)*rca(1) + rca(2)*rca(2) + rca(3)*rca(3)
          dca = SQRT(rca2)

!         *** Create a new matrix block ***

          NULLIFY (scb_block)

          CALL add_block_node(matrix=scb,&
                              block_row=1,&
                              block_col=iatom,&
                              block=scb_block)

!         *** Calculate the current overlap matrix block ***

          prjc = 1

          DO l=0,lppnl

            nprjc = nprj_ppnl(l)*nco(l)

            IF (nprjc == 0) CYCLE

            lc_max = l + 2*(nprj_ppnl(l) - 1)
            lc_min = l
            zetc(1) = alpha_ppnl(l)
            ncoc = ncoset(lc_max)

            DO iset=1,nseta

              IF (rprjc(1) + set_radius_a(iset) < dca) CYCLE

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

!             *** Calculate the primitive overlap integrals ***

              CALL overlap(lc_max,zetc,rprjc,lc_min,&
                           la_max(iset),zeta(1:npgfa(iset),iset),&
                           rpgfa(1:npgfa(iset),iset),la_min(iset),&
                           rca,rca2,sab)

!             *** Contraction step (orbital basis functions) ***

              CALL dgemm("N","N",ncoc,nsgfa(iset),ncoa,&
                         1.0_wp,sab(1,1),SIZE(sab,1),&
                         sphi_a(1,sgfa),SIZE(sphi_a,1),&
                         0.0_wp,work(1,1),SIZE(work,1))

!             *** Contraction step (projector functions) ***

              CALL dgemm("T","N",nprjc,nsgfa(iset),ncoc,&
                         1.0_wp,cprj(1,prjc),SIZE(cprj,1),&
                         work(1,1),SIZE(work,1),&
                         1.0_wp,scb_block(prjc,sgfa),SIZE(scb_block,1))

            END DO

            prjc = prjc + nprjc

          END DO

          scb_neighbor_node => next(scb_neighbor_node)

        END DO

!       *** Add the V_PPNL contribution to the core Hamiltonian matrix ***

        block_node_ca => first_block_node(sca,1)

        DO WHILE (ASSOCIATED(block_node_ca))

          CALL get_block_node(block_node=block_node_ca,&
                              block_col=iatom,&
                              block=sca_block)

          nprjc = SIZE(sca_block,1)
          na = SIZE(sca_block,2)

          CALL dgemm("T","N",nprjc,na,nprjc,&
                     1.0_wp,vprj_ppnl(1,1),nprjc,&
                     sca_block(1,1),nprjc,&
                     0.0_wp,work(1,1),SIZE(work,1))

          block_node_cb => block_node_ca

          DO WHILE (ASSOCIATED(block_node_cb))

            CALL get_block_node(block_node=block_node_cb,&
                                block_col=jatom,&
                                block=scb_block)

            nb = SIZE(scb_block,2)

            IF (jatom > iatom) THEN
              CALL get_block_node(matrix=h%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=h_block)
              IF (.NOT.ASSOCIATED(h_block)) THEN
                CALL add_block_node(matrix=h%matrix,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    block=h_block)
              END IF
              CALL dgemm("T","N",na,nb,nprjc,&
                         1.0_wp,work(1,1),SIZE(work,1),&
                         scb_block(1,1),nprjc,&
                         1.0_wp,h_block(1,1),na)
            ELSE
              CALL get_block_node(matrix=h%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  block=h_block)
              IF (.NOT.ASSOCIATED(h_block)) THEN
                CALL add_block_node(matrix=h%matrix,&
                                    block_row=jatom,&
                                    block_col=iatom,&
                                    block=h_block)
              END IF
              CALL dgemm("T","N",nb,na,nprjc,&
                         1.0_wp,scb_block(1,1),nprjc,&
                         work(1,1),SIZE(work,1),&
                         1.0_wp,h_block(1,1),nb)
            END IF

            block_node_cb => next_block_node(block_node_cb)

          END DO

          block_node_cb => first_block_node(scb,1)

          DO WHILE (ASSOCIATED(block_node_cb))

            CALL get_block_node(block_node=block_node_cb,&
                                block_col=jatom,&
                                block=scb_block)

            IF (jatom >= iatom) THEN
              nb = SIZE(scb_block,2)
              CALL get_block_node(matrix=h%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=h_block)
              IF (.NOT.ASSOCIATED(h_block)) THEN
                CALL add_block_node(matrix=h%matrix,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    block=h_block)
              END IF
              CALL dgemm("T","N",na,nb,nprjc,&
                         1.0_wp,work(1,1),SIZE(work,1),&
                         scb_block(1,1),nprjc,&
                         1.0_wp,h_block(1,1),na)
            END IF

            block_node_cb => next_block_node(block_node_cb)

          END DO

          block_node_ca => next_block_node(block_node_ca)

        END DO

        CALL deallocate_matrix_row(sca,1)
        CALL deallocate_matrix_row(scb,1)

        sca_neighbor_list => next(sca_neighbor_list)
        scb_neighbor_list => next(scb_neighbor_list)

      END DO

      CALL deallocate_matrix(sca)
      CALL deallocate_matrix(scb)

    END DO

!   *** Release work storage ***

    DEALLOCATE (kind_of,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_of")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE add_vppnl

! *****************************************************************************

END MODULE qs_core_hamiltonian
