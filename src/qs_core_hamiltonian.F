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
  USE ai_overlap_new,           ONLY: overlap
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
  USE orbital_pointers,         ONLY: indco,init_orbital_pointers,nco,ncoset
  USE orbital_symbols,          ONLY: cgf_symbol
  USE particle_types,           ONLY: particle_type
  USE qs_environment_types,     ONLY: get_qs_env,&
                                      qs_environment_type,&
                                      set_qs_env
  USE qs_force_types,           ONLY: qs_force_type
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
                                      allocate_matrix_set,&
                                      copy_sparse_to_full_matrix,&
                                      deallocate_matrix,&
                                      deallocate_matrix_set,&
                                      deallocate_matrix_row,&
                                      first_block_node,&
                                      get_block_node,&
                                      get_matrix_info,&
                                      next_block_node,&
                                      real_block_node_type,&
                                      real_matrix_p_type,&
                                      real_matrix_type,&
                                      replicate_matrix_structure,&
                                      set_matrix
  USE string_utilities,         ONLY: compress,&
                                      uppercase
  USE termination,              ONLY: stop_memory,&
                                      stop_program
  USE timings,                  ONLY: timeset,timestop

  IMPLICIT NONE

  PRIVATE

! *** Global types (in this module) ***

  TYPE block_p_type
    REAL(wp), DIMENSION(:,:), POINTER :: block
  END TYPE block_p_type

! *** Global variables (in this module) ***

  TYPE(atomic_kind_type), POINTER   :: atomic_kind
  TYPE(gth_potential_type), POINTER :: gth_potential
  TYPE(gto_basis_set_type), POINTER :: orb_basis_set
  CHARACTER(LEN=12)                 :: cgfsym
  CHARACTER(LEN=60)                 :: name
  INTEGER                           :: ldsab,maxder,maxl,maxsgf,nder,nkind


  TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
  TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_ppl,&
                                                           sca_ppnl,scb_ppnl
  TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
  TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
  TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: h,s,t

  INTEGER, DIMENSION(:), ALLOCATABLE    :: atom_of_kind,kind_of
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: hab,pab,sab,vdab,vadb,work

  REAL(wp), DIMENSION(3)            :: force_a,force_b
  REAL(wp), DIMENSION(:), POINTER   :: set_radius_a,set_radius_b
  INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,lb_max,lb_min,npgfa,&
                                       npgfb,nsgfa,nsgfb
  REAL(wp), DIMENSION(:,:), POINTER :: h_block,p_block,rpgfa,rpgfb,sphi_a,&
                                       sphi_b,w_block,zeta,zetb
  INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa,first_sgfb

! *** Public subroutines ***

  PUBLIC :: build_core_hamiltonian_matrix

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_core_hamiltonian_matrix(qs_env,globenv,calculate_forces,p,w)

!   Purpose: Calculation of the overlap matrix over Cartesian Gaussian
!            functions.

!   History: - Creation (11.03.2002, MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env
    LOGICAL, INTENT(IN)                       :: calculate_forces
    TYPE(real_matrix_type), OPTIONAL, POINTER :: p,w

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE build_core_hamiltonian_matrix (MODULE qs_core_hamiltonian)"

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: sab_orb_neighbor_list,&
                                         sac_ppl_neighbor_list
    TYPE(neighbor_node_type), POINTER :: sab_orb_neighbor_node,&
                                         sac_ppl_neighbor_node

    REAL(wp) :: dab,dac,dbc,f,f0,rab2,rac2,rbc2
    INTEGER  :: ab,ac,atom_a,atom_b,atom_c,da_max,db_max,first_col,first_setb,&
                handle,i,iatom,ikind,inode,ipgf,iset,ishell,istat,jatom,jkind,&
                jpgf,jset,jshell,katom,kkind,maxco,maxdco,maxl,maxppnl,natom,&
                natom_of_kind,ncoa,ncob,nneighbor,nnode,nseta,nsetb,sgfa,sgfb
    LOGICAL  :: gth_potential_present,ppnl_present,&
                return_s_derivatives,return_t_derivatives

    REAL(wp), DIMENSION(3)                  :: rab,rac,rbc
    REAL(wp), DIMENSION(:), ALLOCATABLE     :: alpha_ppl,ppl_radius
    INTEGER, DIMENSION(:), ALLOCATABLE      :: sac_ppl_neighbors
    REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: sac_ppl_r
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: sdab

    TYPE(block_p_type), DIMENSION(:), ALLOCATABLE :: sint,tint

    TYPE cexp_ppl_type
      REAL(wp), DIMENSION(:), POINTER :: c
    END TYPE cexp_ppl_type

    TYPE(cexp_ppl_type), DIMENSION(:), ALLOCATABLE :: cexp_ppl

    TYPE vppl_type
      REAL(wp), DIMENSION(:), POINTER   :: r2
      INTEGER, DIMENSION(:), POINTER    :: neighbor
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
                    scb_ppnl=scb_ppnl,&
                    force=force)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    IF (calculate_forces) THEN
      IF (.NOT.PRESENT(p)) THEN
        CALL stop_program(routine,"The optional argument p "//&
                                  "(density matrix) has to be specified "//&
                                  "for the force calculation")
      END IF
      IF (.NOT.PRESENT(w)) THEN
        CALL stop_program(routine,"The optional argument w "//&
                                  "(W matrix) has to be specified "//&
                                  "for the force calculation")
      END IF
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
          IF (.NOT.ASSOCIATED(force(ikind)%overlap)) THEN
            ALLOCATE (force(ikind)%overlap(3,natom_of_kind),STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"force(ikind)%overlap",&
                               3*natom_of_kind*wp_size)
            END IF
          END IF
          force(ikind)%overlap(:,:) = 0.0_wp
          IF (.NOT.ASSOCIATED(force(ikind)%kinetic)) THEN
            ALLOCATE (force(ikind)%kinetic(3,natom_of_kind),STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"force(ikind)%kinetic",&
                               3*natom_of_kind*wp_size)
            END IF
          END IF
          force(ikind)%kinetic(:,:) = 0.0_wp
          IF (.NOT.ASSOCIATED(force(ikind)%gth_ppl)) THEN
            ALLOCATE (force(ikind)%gth_ppl(3,natom_of_kind),STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"force(ikind)%gth_ppl",&
                               3*natom_of_kind*wp_size)
            END IF
          END IF
          force(ikind)%gth_ppl(:,:) = 0.0_wp
          IF (.NOT.ASSOCIATED(force(ikind)%gth_ppnl)) THEN
            ALLOCATE (force(ikind)%gth_ppnl(3,natom_of_kind),STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine,"force(ikind)%gth_ppnl",&
                               3*natom_of_kind*wp_size)
            END IF
          END IF
          force(ikind)%gth_ppnl(:,:) = 0.0_wp
        END IF
      END DO
      nder = 1
    END IF

    IF (globenv%print%derivatives) nder = 1

    maxder = ncoset(nder)

    return_s_derivatives = (globenv%print%overlap_matrix.AND.&
                            globenv%print%derivatives)

!   *** Initialise the overlap matrix ***

    CALL set_matrix(s(1)%matrix,0.0_wp)

!   *** Allocate the kinetic energy integral ***
!   *** matrix (only needed for printing)    ***

    IF (globenv%print%kinetic_energy_matrix) THEN

      IF (globenv%print%derivatives) THEN
        CALL allocate_matrix_set(t,maxder)
        return_t_derivatives = .TRUE.
      ELSE
        CALL allocate_matrix_set(t,1)
      END IF

      CALL replicate_matrix_structure(source=s(1)%matrix,&
                                      target=t(1)%matrix,&
                                      target_name="KINETIC ENERGY MATRIX")

      DO i=2,SIZE(t)
        cgfsym = cgf_symbol(1,indco(1:3,i))
        name = cgfsym(4:)//" DERIVATIVE OF THE KINETIC ENERGY MATRIX"
        CALL compress(name)
        CALL uppercase(name)
        CALL replicate_matrix_structure(source=t(1)%matrix,&
                                        target=t(i)%matrix,&
                                        target_name=name,&
                                        target_symmetry="antisymmetric")
      END DO

    END IF

!   *** Allocate the core Hamiltonian matrix ***

    CALL allocate_matrix_set(h,1)

    CALL replicate_matrix_structure(source=s(1)%matrix,&
                                    target=h(1)%matrix,&
                                    target_name="CORE HAMILTONIAN MATRIX")

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             maxco=maxco,&
                             maxl=maxl,&
                             maxppnl=maxppnl,&
                             maxsgf=maxsgf)

    CALL init_orbital_pointers(maxl+nder+1)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxdco,&
                             maxder=MAX(1,nder))

    ldsab = MAX(maxco,maxppnl)

    ALLOCATE (hab(maxco,maxco*maxder),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"hab",maxco*maxco*maxder*wp_size)
    END IF
    hab(:,:) = 0.0_wp

    ALLOCATE (sab(ldsab,ldsab*maxder),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"sab",ldsab*ldsab*maxder*wp_size)
    END IF
    sab(:,:) = 0.0_wp

    ALLOCATE (sdab(maxdco,maxco,4),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"sdab",maxdco*maxco*4*wp_size)
    END IF
    sdab(:,:,:) = 0.0_wp

    ALLOCATE (work(ldsab,ldsab*maxder),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"work",ldsab*ldsab*maxder*wp_size)
    END IF
    work(:,:) = 0.0_wp

    IF (calculate_forces) THEN
      ALLOCATE (pab(maxco,maxco),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"pab",maxco*maxco*wp_size)
      END IF
      pab(:,:) = 0.0_wp
    END IF

    ALLOCATE (sint(SIZE(s)),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,"sint",SIZE(s)*int_size)
    END IF
    DO i=1,SIZE(sint)
      NULLIFY (sint(i)%block)
    END DO

    IF (globenv%print%kinetic_energy_matrix) THEN
      ALLOCATE (tint(SIZE(t)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"sint",SIZE(t)*int_size)
      DO i=1,SIZE(tint)
        NULLIFY (tint(i)%block)
      END DO
    END IF

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

        IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             alpha_ppl=alpha_ppl(ikind),&
                             cexp_ppl=cexp_ppl(ikind)%c,&
                             ppl_radius=ppl_radius(ikind))
        END IF

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

          IF (calculate_forces) atom_a = atom_of_kind(iatom)

!         *** Retrieve the data of the SAC_PPL neighbors of atom "iatom" ***

          IF (gth_potential_present) THEN

            DO kkind=1,nkind

              ac = ikind + (kkind - 1)*nkind

              IF (.NOT.ASSOCIATED(sac_ppl(ac)%neighbor_list_set)) THEN
                NULLIFY (vppl(kkind)%r2)
                NULLIFY (vppl(kkind)%neighbor)
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
              ALLOCATE (vppl(kkind)%neighbor(nnode),STAT=istat)
              IF (istat /= 0) THEN
                CALL stop_memory(routine,"vppl(kkind)%neighbor",nnode*int_size)
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
                                       neighbor=vppl(kkind)%neighbor(inode),&
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

            IF (calculate_forces) atom_b = atom_of_kind(jatom)

!           *** Use the symmetry of the first derivatives ***

            IF (iatom == jatom) THEN
              f0 = 1.0_wp
            ELSE
              f0 = 2.0_wp
            END IF

            IF (iatom <= jatom) THEN
              CALL get_block_node(matrix=s(1)%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=sint(1)%block)
            ELSE
              CALL get_block_node(matrix=s(1)%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  block=sint(1)%block)
            END IF

!           *** Check, if the atomic block has to be ***
!           *** calculated by the current processor  ***

            IF (ASSOCIATED(sint(1)%block)) THEN

!             *** Get the pointers to all needed matrix blocks ***

              DO i=2,SIZE(s)
                IF (iatom <= jatom) THEN
                  CALL get_block_node(matrix=s(i)%matrix,&
                                      block_row=iatom,&
                                      block_col=jatom,&
                                      block=sint(i)%block)
                ELSE
                  CALL get_block_node(matrix=s(i)%matrix,&
                                      block_row=jatom,&
                                      block_col=iatom,&
                                      block=sint(i)%block)
                END IF
              END DO

              IF (globenv%print%kinetic_energy_matrix) THEN
                DO i=1,SIZE(t)
                  IF (iatom <= jatom) THEN
                    CALL get_block_node(matrix=t(i)%matrix,&
                                        block_row=iatom,&
                                        block_col=jatom,&
                                        block=tint(i)%block)
                  ELSE
                    CALL get_block_node(matrix=t(i)%matrix,&
                                        block_row=jatom,&
                                        block_col=iatom,&
                                        block=tint(i)%block)
                  END IF
                END DO
              END IF

              IF (iatom <= jatom) THEN
                CALL get_block_node(matrix=h(1)%matrix,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    block=h_block)
                IF (calculate_forces) THEN
                  CALL get_block_node(matrix=p,&
                                      block_row=iatom,&
                                      block_col=jatom,&
                                      block=p_block)
                  CALL get_block_node(matrix=w,&
                                      block_row=iatom,&
                                      block_col=jatom,&
                                      block=w_block)
                END IF
              ELSE
                CALL get_block_node(matrix=h(1)%matrix,&
                                    block_row=jatom,&
                                    block_col=iatom,&
                                    block=h_block)
                IF (calculate_forces) THEN
                  CALL get_block_node(matrix=p,&
                                      block_row=jatom,&
                                      block_col=iatom,&
                                      block=p_block)
                  CALL get_block_node(matrix=w,&
                                      block_row=jatom,&
                                      block_col=iatom,&
                                      block=w_block)
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

!                 *** No derivatives are needed for the diagonal blocks ***
!                 *** of the two-center integrals except the overlap    ***
!                 *** integrals which are needed for the calculation of ***
!                 *** the kinetic energy integrals                      ***

                  IF (iatom == jatom) THEN
                    da_max = 0
                  ELSE
                    da_max = nder
                  END IF

                  IF (calculate_forces.AND.(iatom /= jatom)) THEN

!                   *** Decontract W matrix block ***

                    IF (iatom <= jatom) THEN
                      CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 w_block(sgfa,sgfb),SIZE(w_block,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                    ELSE
                      CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 w_block(sgfb,sgfa),SIZE(w_block,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                    END IF
                    CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                               1.0_wp,work(1,1),SIZE(work,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_wp,pab(1,1),SIZE(pab,1))

!                   *** Calculate the primitive overlap integrals ***
!                   *** and the corresponding force contribution  ***

                    CALL overlap(la_max(iset),zeta(1:npgfa(iset),iset),&
                                 rpgfa(1:npgfa(iset),iset),la_min(iset),&
                                 lb_max(jset),zetb(1:npgfb(jset),jset),&
                                 rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                                 rab,rab2,dab,sab,da_max,return_s_derivatives,&
                                 sdab,pab,force_a)

                    force(ikind)%overlap(:,atom_a) =&
                      force(ikind)%overlap(:,atom_a) - f0*force_a(:)
                    force(jkind)%overlap(:,atom_b) =&
                      force(jkind)%overlap(:,atom_b) + f0*force_a(:)

                  ELSE

!                   *** Calculate the primitive overlap integrals ***

                    CALL overlap(la_max(iset),zeta(1:npgfa(iset),iset),&
                                 rpgfa(1:npgfa(iset),iset),la_min(iset),&
                                 lb_max(jset),zetb(1:npgfb(jset),jset),&
                                 rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                                 rab,rab2,dab,sab,da_max,return_s_derivatives,&
                                 sdab)

                  END IF

!                 *** Contraction step (overlap matrix) ***

                  CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                             1.0_wp,sab(1,1),SIZE(sab,1),&
                             sphi_b(1,sgfb),SIZE(sphi_b,1),&
                             0.0_wp,work(1,1),SIZE(work,1))
                  IF (iatom <= jatom) THEN
                    CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                               1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               work(1,1),SIZE(work,1),&
                               1.0_wp,sint(1)%block(sgfa,sgfb),&
                               SIZE(sint(1)%block,1))
                  ELSE
                    CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                               1.0_wp,work(1,1),SIZE(work,1),&
                               sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               1.0_wp,sint(1)%block(sgfb,sgfa),&
                               SIZE(sint(1)%block,1))
                  END IF

!                 *** Contraction step (derivatives of the overlap matrix) ***

                  IF (iatom /= jatom) THEN
                    DO i=2,SIZE(s)
                      first_col = (i - 1)*SIZE(sab,1) + 1
                      CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                                 1.0_wp,sab(1,first_col),SIZE(sab,1),&
                                 sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                      IF (iatom <= jatom) THEN
                        CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                   1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                   work(1,1),SIZE(work,1),&
                                   1.0_wp,sint(i)%block(sgfa,sgfb),&
                                   SIZE(sint(i)%block,1))
                      ELSE
                        CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                   -1.0_wp,work(1,1),SIZE(work,1),&
                                   sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                   1.0_wp,sint(i)%block(sgfb,sgfa),&
                                   SIZE(sint(i)%block,1))
                      END IF
                    END DO
                  END IF

!                 *** Decontract density matrix block ***

                  IF (calculate_forces) THEN
                    IF (iatom <= jatom) THEN
                      CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 p_block(sgfa,sgfb),SIZE(p_block,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                    ELSE
                      CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 p_block(sgfb,sgfa),SIZE(p_block,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                    END IF
                    CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                               1.0_wp,work(1,1),SIZE(work,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_wp,pab(1,1),SIZE(pab,1))
                  END IF

                  IF (calculate_forces.AND.(iatom /= jatom)) THEN

!                   *** Calculate the primitive kinetic energy integrals ***
!                   *** and the corresponding force contribution         ***

                    CALL kinetic(la_max(iset),zeta(1:npgfa(iset),iset),&
                                 rpgfa(1:npgfa(iset),iset),la_min(iset),&
                                 lb_max(jset),zetb(1:npgfb(jset),jset),&
                                 rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                                 rab,rab2,dab,sdab,hab,da_max,&
                                 return_t_derivatives,pab,force_a)

                    force(ikind)%kinetic(:,atom_a) =&
                      force(ikind)%kinetic(:,atom_a) + f0*force_a(:)
                    force(jkind)%kinetic(:,atom_b) =&
                      force(jkind)%kinetic(:,atom_b) - f0*force_a(:)

                  ELSE

!                   *** Calculate the primitive kinetic energy integrals ***

                    CALL kinetic(la_max(iset),zeta(1:npgfa(iset),iset),&
                                 rpgfa(1:npgfa(iset),iset),la_min(iset),&
                                 lb_max(jset),zetb(1:npgfb(jset),jset),&
                                 rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                                 rab,rab2,dab,sdab,hab,da_max,&
                                 return_t_derivatives)

                  END IF

!                 *** Contraction step (kinetic energy integral matrix ***
!                 *** and its derivatives (only needed for printing)   ***

                  IF (globenv%print%kinetic_energy_matrix) THEN
                    CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                               1.0_wp,hab(1,1),SIZE(hab,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_wp,work(1,1),SIZE(work,1))
                    IF (iatom <= jatom) THEN
                      CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 work(1,1),SIZE(work,1),&
                                 1.0_wp,tint(1)%block(sgfa,sgfb),&
                                 SIZE(tint(1)%block,1))
                    ELSE
                      CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                 1.0_wp,work(1,1),SIZE(work,1),&
                                 sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 1.0_wp,tint(1)%block(sgfb,sgfa),&
                                 SIZE(tint(1)%block,1))
                    END IF
                    IF (iatom /= jatom) THEN
                      DO i=2,SIZE(t)
                        first_col = (i - 1)*SIZE(hab,1) + 1
                        CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                                   1.0_wp,hab(1,first_col),SIZE(hab,1),&
                                   sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                   0.0_wp,work(1,1),SIZE(work,1))
                        IF (iatom <= jatom) THEN
                          CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                     1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                     work(1,1),SIZE(work,1),&
                                     1.0_wp,tint(i)%block(sgfa,sgfb),&
                                     SIZE(tint(i)%block,1))
                        ELSE
                          CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                     -1.0_wp,work(1,1),SIZE(work,1),&
                                     sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                     1.0_wp,tint(i)%block(sgfb,sgfa),&
                                     SIZE(tint(i)%block,1))
                        END IF
                      END DO
                    END IF
                  END IF

!                 *** Loop over the GTH pseudo potential atoms (local part) ***

                  DO kkind=1,nkind

                    IF (.NOT.ASSOCIATED(vppl(kkind)%neighbor)) CYCLE

                    DO inode=1,SIZE(vppl(kkind)%neighbor)

                      katom = vppl(kkind)%neighbor(inode)
                      rac(:) = vppl(kkind)%r(:,inode)
                      rac2 = vppl(kkind)%r2(inode)
                      dac = SQRT(rac2)
                      IF (set_radius_a(iset) + ppl_radius(kkind) < dac) CYCLE

                      rbc(:) = rac(:) - rab(:)
                      rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                      dbc = SQRT(rbc2)
                      IF (set_radius_b(jset) + ppl_radius(kkind) < dbc) CYCLE

!                     *** Calculate the GTH pseudo potential forces ***

                      IF (calculate_forces) THEN

                        CALL overlap_ppl(&
                          la_max(iset),zeta(1:npgfa(iset),iset),&
                          rpgfa(1:npgfa(iset),iset),la_min(iset),&
                          lb_max(jset),zetb(1:npgfb(jset),jset),&
                          rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                          cexp_ppl(kkind)%c,alpha_ppl(kkind),&
                          ppl_radius(kkind),&
                          rab,rab2,dab,&
                          rac,rac2,dac,&
                          rbc,rbc2,dbc,&
                          hab,nder,nder,.FALSE.,&
                          pab,force_a,force_b)

!                       *** The derivatives w.r.t. atomic center c are    ***
!                       *** calculated using the translational invariance ***
!                       *** of the first derivatives                      ***

                        atom_c = atom_of_kind(katom)

                        force(ikind)%gth_ppl(:,atom_a) =&
                          force(ikind)%gth_ppl(:,atom_a) + f0*force_a(:)
                        force(kkind)%gth_ppl(:,atom_c) =&
                          force(kkind)%gth_ppl(:,atom_c) - f0*force_a(:)

                        force(jkind)%gth_ppl(:,atom_b) =&
                          force(jkind)%gth_ppl(:,atom_b) + f0*force_b(:)
                        force(kkind)%gth_ppl(:,atom_c) =&
                          force(kkind)%gth_ppl(:,atom_c) - f0*force_b(:)

                      ELSE

                        CALL overlap_ppl(&
                          la_max(iset),zeta(1:npgfa(iset),iset),&
                          rpgfa(1:npgfa(iset),iset),la_min(iset),&
                          lb_max(jset),zetb(1:npgfb(jset),jset),&
                          rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                          cexp_ppl(kkind)%c,alpha_ppl(kkind),&
                          ppl_radius(kkind),&
                          rab,rab2,dab,&
                          rac,rac2,dac,&
                          rbc,rbc2,dbc,&
                          hab,0,0,.FALSE.)

                      END IF

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
                CALL symmetrize_matrix(sint(1)%block,"upper_to_lower")
                DO i=2,SIZE(s)
                  sint(i)%block(:,:) = 0.0_wp
                END DO
                IF (globenv%print%kinetic_energy_matrix) THEN
                  CALL symmetrize_matrix(tint(1)%block,"upper_to_lower")
                  DO i=2,SIZE(t)
                    tint(i)%block(:,:) = 0.0_wp
                  END DO
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

    IF (gth_potential_present) CALL add_vppnl(calculate_forces,p)

!   *** Put the core Hamiltonian matrix in the QS environment ***

    CALL set_qs_env(qs_env=qs_env,s=s,h=h)

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

    DEALLOCATE (sdab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"sdab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work")

    IF (calculate_forces) THEN
      DEALLOCATE (atom_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")
      DEALLOCATE (pab,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"pab")
    END IF

    DEALLOCATE (sint,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"sint")

    IF (globenv%print%kinetic_energy_matrix) THEN
      DEALLOCATE (tint,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"tint")
    END IF

    CALL timestop(0.0_wp,handle)

!   *** Print the overlap integral matrix, if requested ***

    IF (globenv%print%overlap_matrix) THEN
      CALL write_sparse_matrix(s(1)%matrix,4,6,qs_env,globenv)
      IF (globenv%print%derivatives) THEN
        DO i=2,SIZE(s)
          CALL write_sparse_matrix(s(i)%matrix,4,6,qs_env,globenv)
        END DO
      END IF
    END IF

!   *** Print the kinetic energy integral matrix, if requested ***

    IF (globenv%print%kinetic_energy_matrix) THEN
      DO i=1,SIZE(t)
        CALL write_sparse_matrix(t(i)%matrix,4,6,qs_env,globenv)
      END DO
    END IF

!   *** Print the core Hamiltonian matrix, if requested ***

    IF (globenv%print%core_hamiltonian_matrix) THEN
      CALL write_sparse_matrix(h(1)%matrix,4,6,qs_env,globenv)
    END IF

  END SUBROUTINE build_core_hamiltonian_matrix

! *****************************************************************************

  SUBROUTINE add_vppnl(calculate_forces,p)

!   Purpose: Add the non-local contribution of the GTH pseudo potential to the
!            core Hamiltonian matrix.

!   History: - Creation (24.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(real_matrix_type), OPTIONAL, POINTER :: p
    LOGICAL, INTENT(IN)                       :: calculate_forces

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE add_vppnl (MODULE qs_core_hamiltonian)"

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: sca_neighbor_list,scb_neighbor_list
    TYPE(neighbor_node_type), POINTER :: sca_neighbor_node,scb_neighbor_node

    REAL(wp) :: dca,f,f0,rca2
    INTEGER  :: atom_a,atom_b,atom_c,current_kind,first_col,first_prjc,handle,&
                i,iatom,ikind,iset,istat,jatom,jkind,katom,kkind,l,lc_max,&
                lc_min,lppnl,na,natom,nb,ncoa,ncoc,nppnl,nprjc,nseta,nsgf,&
                prjc,sgfa
    LOGICAL  :: ppnl_present

    REAL(wp), DIMENSION(1) :: rprjc,zetc
    INTEGER, DIMENSION(1)  :: first_prj,last_prj
    REAL(wp), DIMENSION(3) :: rca
    INTEGER, DIMENSION(3)  :: cell_c

    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: sca,scb
    TYPE(block_p_type), DIMENSION(:), ALLOCATABLE   :: scaint,scbint

    REAL(wp), DIMENSION(:), POINTER       :: alpha_ppnl
    INTEGER, DIMENSION(:), POINTER        :: nprj_ppnl
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: v_block
    REAL(wp), DIMENSION(:,:), POINTER     :: cprj,h_block,sca_block,scb_block,&
                                             vprj_ppnl

!   *** Local data types ***

    TYPE real_block_node_p_type
      TYPE(real_block_node_type), POINTER :: block_node
    END TYPE real_block_node_p_type
    TYPE(real_block_node_type), POINTER :: block_node_dum

    TYPE(real_block_node_p_type), DIMENSION(:), ALLOCATABLE :: ca,cb

!   ---------------------------------------------------------------------------

    CALL timeset("add_vppnl","I","",handle)

!   *** Allocate and initialise work storage ***

    NULLIFY (sca,scb)
    CALL allocate_matrix_set(sca,maxder)
    CALL allocate_matrix_set(scb,maxder)

    IF (calculate_forces) THEN
      ALLOCATE (v_block(maxsgf,maxsgf),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"v_block",maxsgf*maxsgf*wp_size)
      v_block(:,:) = 0.0_wp
    END IF

    ALLOCATE (scaint(maxder),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scaint",maxder*int_size)
    DO i=1,SIZE(scaint)
      NULLIFY (scaint(i)%block)
    END DO

    ALLOCATE (scbint(maxder),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scbint",maxder*int_size)
    DO i=1,SIZE(scbint)
      NULLIFY (scbint(i)%block)
    END DO

    ALLOCATE (ca(maxder),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ca",maxder*int_size)
    DO i=1,SIZE(ca)
      NULLIFY (ca(i)%block_node)
    END DO

    ALLOCATE (cb(maxder),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"cb",maxder*int_size)
    DO i=1,SIZE(cb)
      NULLIFY (cb(i)%block_node)
    END DO

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom,&
                             nsgf=nsgf)

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

!     *** Allocate an CA projector (ppnl) overlap matrix row ***

      CALL allocate_matrix(matrix=sca(1)%matrix,&
                           nblock_row=1,&
                           nblock_col=natom,&
                           nrow=nppnl,&
                           ncol=nsgf,&
                           first_row=first_prj(:),&
                           last_row=last_prj(:),&
                           first_col=particle_set(:)%first_sgf,&
                           last_col=particle_set(:)%last_sgf,&
                           matrix_name="SCA GTH-PPNL OVERLAP MATRIX ROW",&
                           matrix_symmetry="no symmetry")

      DO i=2,maxder
        cgfsym = cgf_symbol(1,indco(1:3,i))
        name = cgfsym(4:)//" DERIVATIVE OF THE SCA GTH-PPNL OVERLAP MATRIX ROW"
        CALL compress(name)
        CALL uppercase(name)
        CALL replicate_matrix_structure(source=sca(1)%matrix,&
                                        target=sca(i)%matrix,&
                                        target_name=name,&
                                        target_symmetry="no symmetry")
      END DO

!     *** Allocate the CB projector (ppnl) overlap matrices row ***

      CALL allocate_matrix(matrix=scb(1)%matrix,&
                           nblock_row=1,&
                           nblock_col=natom,&
                           nrow=nppnl,&
                           ncol=nsgf,&
                           first_row=first_prj(:),&
                           last_row=last_prj(:),&
                           first_col=particle_set(:)%first_sgf,&
                           last_col=particle_set(:)%last_sgf,&
                           matrix_name="SCB GTH-PPNL OVERLAP MATRIX ROW",&
                           matrix_symmetry="no symmetry")

      DO i=2,maxder
        cgfsym = cgf_symbol(1,indco(1:3,i))
        name = cgfsym(4:)//" DERIVATIVE OF THE SCB GTH-PPNL OVERLAP MATRIX ROW"
        CALL compress(name)
        CALL uppercase(name)
        CALL replicate_matrix_structure(source=scb(1)%matrix,&
                                        target=scb(i)%matrix,&
                                        target_name=name,&
                                        target_symmetry="no symmetry")
      END DO

!     *** Calculate the overlap integrals between the orbital ***
!     *** basis functions and the GTH(nl) projector functions ***

      sca_neighbor_list => first_list(sca_ppnl(kkind)%neighbor_list_set)
      scb_neighbor_list => first_list(scb_ppnl(kkind)%neighbor_list_set)

      DO WHILE (ASSOCIATED(sca_neighbor_list))

        current_kind = 0

        CALL get_neighbor_list(neighbor_list=sca_neighbor_list,&
                               atom=katom)

        IF (calculate_forces) atom_c = atom_of_kind(katom)

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

          DO i=1,maxder
            NULLIFY (scaint(i)%block)
            CALL add_block_node(matrix=sca(i)%matrix,&
                                block_row=1,&
                                block_col=iatom,&
                                block=scaint(i)%block)
          END DO

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
                           rca,rca2,dca,sab,nder,.TRUE.)

!             *** Contraction step (orbital basis functions) ***

              CALL dgemm("N","N",ncoc,nsgfa(iset),ncoa,&
                         1.0_wp,sab(1,1),SIZE(sab,1),&
                         sphi_a(1,sgfa),SIZE(sphi_a,1),&
                         0.0_wp,work(1,1),SIZE(work,1))

!             *** Contraction step (projector functions) ***

              CALL dgemm("T","N",nprjc,nsgfa(iset),ncoc,&
                         1.0_wp,cprj(1,prjc),SIZE(cprj,1),&
                         work(1,1),SIZE(work,1),&
                         1.0_wp,scaint(1)%block(prjc,sgfa),&
                         SIZE(scaint(1)%block,1))

!             *** Contraction step (derivatives) ***

              DO i=2,maxder
                first_col = (i - 1)*SIZE(sab,1) + 1
                CALL dgemm("N","N",ncoc,nsgfa(iset),ncoa,&
                           -1.0_wp,sab(1,first_col),SIZE(sab,1),&
                           sphi_a(1,sgfa),SIZE(sphi_a,1),&
                           0.0_wp,work(1,1),SIZE(work,1))
                CALL dgemm("T","N",nprjc,nsgfa(iset),ncoc,&
                           1.0_wp,cprj(1,prjc),SIZE(cprj,1),&
                           work(1,1),SIZE(work,1),&
                           1.0_wp,scaint(i)%block(prjc,sgfa),&
                           SIZE(scaint(i)%block,1))
              END DO

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

          DO i=1,maxder
            NULLIFY (scbint(i)%block)
            CALL add_block_node(matrix=scb(i)%matrix,&
                                block_row=1,&
                                block_col=iatom,&
                                block=scbint(i)%block)
          END DO

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
                           rca,rca2,dca,sab,nder,.TRUE.)

!             *** Contraction step (orbital basis functions) ***

              CALL dgemm("N","N",ncoc,nsgfa(iset),ncoa,&
                         1.0_wp,sab(1,1),SIZE(sab,1),&
                         sphi_a(1,sgfa),SIZE(sphi_a,1),&
                         0.0_wp,work(1,1),SIZE(work,1))

!             *** Contraction step (projector functions) ***

              CALL dgemm("T","N",nprjc,nsgfa(iset),ncoc,&
                         1.0_wp,cprj(1,prjc),SIZE(cprj,1),&
                         work(1,1),SIZE(work,1),&
                         1.0_wp,scbint(1)%block(prjc,sgfa),&
                         SIZE(scbint(1)%block,1))

!             *** Contraction step (derivatives) ***

              DO i=2,maxder
                first_col = (i - 1)*SIZE(sab,1) + 1
                CALL dgemm("N","N",ncoc,nsgfa(iset),ncoa,&
                           -1.0_wp,sab(1,first_col),SIZE(sab,1),&
                           sphi_a(1,sgfa),SIZE(sphi_a,1),&
                           0.0_wp,work(1,1),SIZE(work,1))
                CALL dgemm("T","N",nprjc,nsgfa(iset),ncoc,&
                           1.0_wp,cprj(1,prjc),SIZE(cprj,1),&
                           work(1,1),SIZE(work,1),&
                           1.0_wp,scbint(i)%block(prjc,sgfa),&
                           SIZE(scbint(i)%block,1))
              END DO

            END DO

            prjc = prjc + nprjc

          END DO

          scb_neighbor_node => next(scb_neighbor_node)

        END DO

!       *** Add the V_PPNL contribution to the core Hamiltonian matrix ***

        DO i=1,maxder
          ca(i)%block_node => first_block_node(sca(i)%matrix,1)
        END DO

        DO WHILE (ASSOCIATED(ca(1)%block_node))

          block_node_dum=>ca(1)%block_node
          CALL get_block_node(block_node=block_node_dum,&
                              block_col=iatom,&
                              block=scaint(1)%block)
          DO i=2,maxder
            block_node_dum=>ca(i)%block_node
            CALL get_block_node(block_node=block_node_dum,&
                                block=scaint(i)%block)
          END DO

          IF (calculate_forces) atom_a = atom_of_kind(iatom)
          ikind = kind_of(iatom)
          nprjc = SIZE(scaint(1)%block,1)
          na = SIZE(scaint(1)%block,2)

          DO i=1,maxder
            first_col = (i - 1)*SIZE(work,1) + 1
            CALL dgemm("T","N",nprjc,na,nprjc,&
                       1.0_wp,vprj_ppnl(1,1),nprjc,&
                       scaint(i)%block(1,1),nprjc,&
                       0.0_wp,work(1,first_col),SIZE(work,1))
          END DO

          DO i=1,maxder
            cb(i)%block_node => ca(i)%block_node
          END DO

          DO WHILE (ASSOCIATED(cb(1)%block_node))
            
            block_node_dum => cb(1)%block_node
            CALL get_block_node(block_node=block_node_dum,&
                                block_col=jatom,&
                                block=scbint(1)%block)
            DO i=2,maxder
              block_node_dum => cb(i)%block_node
              CALL get_block_node(block_node=block_node_dum,&
                                  block=scbint(i)%block)
            END DO

            IF (calculate_forces) atom_b = atom_of_kind(jatom)
            jkind = kind_of(jatom)
            nb = SIZE(scbint(1)%block,2)

            IF (jatom > iatom) THEN
              CALL get_block_node(matrix=h(1)%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=h_block)
              IF (.NOT.ASSOCIATED(h_block)) THEN
                CALL add_block_node(matrix=h(1)%matrix,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    block=h_block)
              END IF
              CALL dgemm("T","N",na,nb,nprjc,&
                         1.0_wp,work(1,1),SIZE(work,1),&
                         scbint(1)%block(1,1),nprjc,&
                         1.0_wp,h_block(1,1),na)
              IF (calculate_forces) THEN
                CALL get_block_node(matrix=p,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    block=p_block)
                DO i=1,3
                  first_col = i*SIZE(work,1) + 1
                  CALL dgemm("T","N",na,nb,nprjc,&
                             1.0_wp,work(1,first_col),SIZE(work,1),&
                             scbint(1)%block(1,1),nprjc,&
                             0.0_wp,v_block(1,1),SIZE(v_block,1))
                  f = 2.0_wp*trace(na,nb,p_block,v_block)
                  force(ikind)%gth_ppnl(i,atom_a) =&
                    force(ikind)%gth_ppnl(i,atom_a) + f
                  force(kkind)%gth_ppnl(i,atom_c) =&
                    force(kkind)%gth_ppnl(i,atom_c) - f
                  CALL dgemm("T","N",na,nb,nprjc,&
                             1.0_wp,work(1,1),SIZE(work,1),&
                             scbint(i+1)%block(1,1),nprjc,&
                             0.0_wp,v_block(1,1),SIZE(v_block,1))
                  f = 2.0_wp*trace(na,nb,p_block,v_block)
                  force(jkind)%gth_ppnl(i,atom_b) =&
                    force(jkind)%gth_ppnl(i,atom_b) + f
                  force(kkind)%gth_ppnl(i,atom_c) =&
                    force(kkind)%gth_ppnl(i,atom_c) - f
                END DO
              END IF
            ELSE
              CALL get_block_node(matrix=h(1)%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  block=h_block)
              IF (.NOT.ASSOCIATED(h_block)) THEN
                CALL add_block_node(matrix=h(1)%matrix,&
                                    block_row=jatom,&
                                    block_col=iatom,&
                                    block=h_block)
              END IF
              CALL dgemm("T","N",nb,na,nprjc,&
                         1.0_wp,scbint(1)%block(1,1),nprjc,&
                         work(1,1),SIZE(work,1),&
                         1.0_wp,h_block(1,1),nb)
              IF (calculate_forces) THEN
                IF (iatom == jatom) THEN
                  f0 = 1.0_wp
                ELSE
                  f0 = 2.0_wp
                END IF
                CALL get_block_node(matrix=p,&
                                    block_row=jatom,&
                                    block_col=iatom,&
                                    block=p_block)
                DO i=1,3
                  first_col = i*SIZE(work,1) + 1
                  CALL dgemm("T","N",nb,na,nprjc,&
                             1.0_wp,scbint(1)%block(1,1),nprjc,&
                             work(1,first_col),SIZE(work,1),&
                             0.0_wp,v_block(1,1),SIZE(v_block,1))
                  f = f0*trace(nb,na,p_block,v_block)
                  force(ikind)%gth_ppnl(i,atom_a) =&
                    force(ikind)%gth_ppnl(i,atom_a) + f
                  force(kkind)%gth_ppnl(i,atom_c) =&
                    force(kkind)%gth_ppnl(i,atom_c) - f
                  CALL dgemm("T","N",nb,na,nprjc,&
                             1.0_wp,scbint(i+1)%block(1,1),nprjc,&
                             work(1,1),SIZE(work,1),&
                             0.0_wp,v_block(1,1),SIZE(v_block,1))
                  f = f0*trace(nb,na,p_block,v_block)
                  force(jkind)%gth_ppnl(i,atom_b) =&
                    force(jkind)%gth_ppnl(i,atom_b) + f
                  force(kkind)%gth_ppnl(i,atom_c) =&
                    force(kkind)%gth_ppnl(i,atom_c) - f
                END DO
              END IF
            END IF

            DO i=1,maxder
              block_node_dum => cb(i)%block_node
              cb(i)%block_node => next_block_node(block_node_dum)
            END DO

          END DO

          DO i=1,maxder
            cb(i)%block_node => first_block_node(scb(i)%matrix,1)
          END DO

          DO WHILE (ASSOCIATED(cb(1)%block_node))

            block_node_dum => cb(1)%block_node
            CALL get_block_node(block_node=block_node_dum,&
                                block_col=jatom,&
                                block=scbint(1)%block)
            DO i=2,maxder
              block_node_dum => cb(i)%block_node
              CALL get_block_node(block_node=block_node_dum,&
                                  block=scbint(i)%block)
            END DO

            IF (calculate_forces) atom_b = atom_of_kind(jatom)
            jkind = kind_of(jatom)

            IF (jatom >= iatom) THEN
              nb = SIZE(scbint(1)%block,2)
              CALL get_block_node(matrix=h(1)%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=h_block)
              IF (.NOT.ASSOCIATED(h_block)) THEN
                CALL add_block_node(matrix=h(1)%matrix,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    block=h_block)
              END IF
              CALL dgemm("T","N",na,nb,nprjc,&
                         1.0_wp,work(1,1),SIZE(work,1),&
                         scbint(1)%block(1,1),nprjc,&
                         1.0_wp,h_block(1,1),na)
              IF (calculate_forces) THEN
                IF (iatom == jatom) THEN
                  f0 = 1.0_wp
                ELSE
                  f0 = 2.0_wp
                END IF
                CALL get_block_node(matrix=p,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    block=p_block)
                DO i=1,3
                  first_col = i*SIZE(work,1) + 1
                  CALL dgemm("T","N",na,nb,nprjc,&
                             1.0_wp,work(1,first_col),SIZE(work,1),&
                             scbint(1)%block(1,1),nprjc,&
                             0.0_wp,v_block(1,1),SIZE(v_block,1))
                  f = f0*trace(na,nb,p_block,v_block)
                  force(ikind)%gth_ppnl(i,atom_a) =&
                    force(ikind)%gth_ppnl(i,atom_a) + f
                  force(kkind)%gth_ppnl(i,atom_c) =&
                    force(kkind)%gth_ppnl(i,atom_c) - f
                  CALL dgemm("T","N",na,nb,nprjc,&
                             1.0_wp,work(1,1),SIZE(work,1),&
                             scbint(i+1)%block(1,1),nprjc,&
                             0.0_wp,v_block(1,1),SIZE(v_block,1))
                  f = f0*trace(na,nb,p_block,v_block)
                  force(jkind)%gth_ppnl(i,atom_b) =&
                    force(jkind)%gth_ppnl(i,atom_b) + f
                  force(kkind)%gth_ppnl(i,atom_c) =&
                    force(kkind)%gth_ppnl(i,atom_c) - f
                END DO
              END IF
            END IF

            DO i=1,maxder
              block_node_dum => cb(i)%block_node
              cb(i)%block_node => next_block_node(block_node_dum)
            END DO

          END DO

          DO i=1,maxder
            block_node_dum => ca(i)%block_node
            ca(i)%block_node => next_block_node(block_node_dum)
          END DO

        END DO

        DO i=1,maxder
          CALL deallocate_matrix_row(sca(i)%matrix,1)
          CALL deallocate_matrix_row(scb(i)%matrix,1)
        END DO

        sca_neighbor_list => next(sca_neighbor_list)
        scb_neighbor_list => next(scb_neighbor_list)

      END DO

      DO i=1,maxder
        CALL deallocate_matrix(sca(i)%matrix)
        CALL deallocate_matrix(scb(i)%matrix)
      END DO

    END DO

!   *** Release work storage ***

    CALL deallocate_matrix_set(sca)
    CALL deallocate_matrix_set(scb)

    IF (calculate_forces) THEN
      DEALLOCATE (v_block,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"vab")
    END IF

    DEALLOCATE (scaint,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scaint")

    DEALLOCATE (scbint,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scbint")

    DEALLOCATE (ca,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"ca")

    DEALLOCATE (cb,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"cb")

    DEALLOCATE (kind_of,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_of")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE add_vppnl

! *****************************************************************************

  FUNCTION trace(nrow,ncol,a,b) RESULT(tr)

    REAL(wp), DIMENSION(:,:), INTENT(IN) :: a,b
    INTEGER, INTENT(IN)                  :: ncol,nrow

    REAL(wp) :: tr

    INTEGER :: i,j

!   ---------------------------------------------------------------------------

    tr = 0.0_wp

    DO j=1,ncol
      DO i=1,nrow
        tr = tr + a(i,j)*b(i,j)
      END DO
    END DO

  END FUNCTION trace

! *****************************************************************************

END MODULE qs_core_hamiltonian
