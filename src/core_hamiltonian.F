!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/core_hamiltonian [1.0] *
!!
!!   NAME
!!     core_hamiltonian
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
!!     Matthias Krack (14.09.2000)
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

MODULE core_hamiltonian

! *****************************************************************************

! Index:

! SUBROUTINE build_core_hamiltonian_matrix(globenv)
! SUBROUTINE distribute_s_matrix()
! SUBROUTINE build_h_matrix()
! SUBROUTINE add_ppnl_matrix()

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE atoms,                    ONLY: atom_info,natom,nsgf
  USE atomic_kinds,             ONLY: all_potential_present,&
                                      gth_potential_present,&
                                      kind_info,nkind
  USE basis_set_types,          ONLY: maxco,maxsgf
  USE global_types,             ONLY: global_environment_type
  USE matrix_types,             ONLY: add_block_node,&
                                      allocate_matrix,&
                                      deallocate_matrix,&
                                      deallocate_matrix_row,&
                                      first_block_node,&
                                      get_block_node,&
                                      get_matrix_info,&
                                      next_block_node,&
                                      real_block_node_type,&
                                      real_matrix_set_type,&
                                      replicate_matrix_structure
  USE memory_utilities,         ONLY: reallocate
  USE message_passing,          ONLY: mp_sum
  USE method_specifications,    ONLY: allchem,gpw,maxder
  USE neighbor_list_types,      ONLY: extract_neighbor_list,&
                                      find_neighbor_list,&
                                      first_neighbor_list,&
                                      first_neighbor_node,&
                                      get_neighbor_list,&
                                      get_neighbor_list_set_info,&
                                      get_neighbor_node,&
                                      neighbor_list_type,&
                                      neighbor_node_type,&
                                      next_neighbor_list,&
                                      next_neighbor_node
  USE nl,                       ONLY: orb_neighbor_list_set,&
                                      verfc_neighbor_list_set,&
                                      vppl_neighbor_list_set,&
                                      vppnl_ca_neighbor_list_set
  USE om_utilities,             ONLY: write_spherical_matrix
  USE orbital_pointers,         ONLY: nco,ncoset
  USE external_potential_types, ONLY: maxlppnl
  USE termination,              ONLY: stop_memory
  USE timings,                  ONLY: timeset,timestop

  IMPLICIT NONE

  PRIVATE

! *** Global variables ***

  TYPE(real_matrix_set_type) :: h,s,sca_prj_ppnl,scb_prj_ppnl,t
  LOGICAL                    :: init_pointers_done = .FALSE.

  REAL(wp), DIMENSION(:,:), POINTER :: hab,sab,work

! *** Public variables ***

  PUBLIC :: h,s

! *** Public subroutines ***

  PUBLIC :: build_core_hamiltonian_matrix

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_core_hamiltonian_matrix(globenv)

!   Purpose: Calculation of the core Hamiltonian integral matrix.

!   History: - Creation (14.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE build_core_hamiltonian_matrix (MODULE core_hamiltonian)"

!   *** Local variables ***

    INTEGER :: handle,istat,maxfun

!   ---------------------------------------------------------------------------

    CALL timeset("build_core_hamiltonian_matrix","I","",handle)

!   *** Initialize the matrix pointers ***

    IF (.NOT.init_pointers_done) THEN
      NULLIFY (h%matrix)
      NULLIFY (s%matrix)
      NULLIFY (t%matrix)
      NULLIFY (sca_prj_ppnl%matrix)
      NULLIFY (scb_prj_ppnl%matrix)
      init_pointers_done = .TRUE.
    END IF

!   *** Allocate work storage ***

    maxfun = MAX(maxco,ncoset(maxlppnl))

    NULLIFY (hab,sab,work)

    hab => reallocate(hab,1,maxco,1,maxco)
    sab => reallocate(sab,1,maxfun,1,maxco)
    work => reallocate(work,1,maxfun,1,maxsgf)

!   *** Allocate the overlap matrix and distribute the atomic blocks ***

    CALL distribute_s_matrix(globenv)

!   *** Allocate the core Hamiltonian matrix ***

    CALL replicate_matrix_structure(source=s%matrix,&
                                    target=h%matrix,&
                                    target_name="CORE HAMILTONIAN MATRIX")

    IF (globenv%print%kinetic_energy_matrix) THEN
      CALL replicate_matrix_structure(source=s%matrix,&
                                      target=t%matrix,&
                                      target_name="KINETIC ENERGY MATRIX")
    END IF

!   *** Calculate the overlap and the core Hamiltonian matrix ***

    CALL build_h_matrix(globenv)

!   *** Printing ***

    IF (globenv%print%overlap_matrix) THEN
      CALL write_spherical_matrix(s%matrix,4,6,globenv)
    END IF

    IF (globenv%print%kinetic_energy_matrix) THEN
      CALL write_spherical_matrix(t%matrix,4,6,globenv)
    END IF

    IF (globenv%print%core_hamiltonian_matrix) THEN
      CALL write_spherical_matrix(h%matrix,4,6,globenv)
    END IF

!   *** Release work storage ***

    IF (ASSOCIATED(hab)) THEN
      DEALLOCATE (hab,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"hab")
      END IF
    END IF

    IF (ASSOCIATED(sab)) THEN
      DEALLOCATE (sab,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"sab")
      END IF
    END IF

    IF (ASSOCIATED(work)) THEN
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"work")
      END IF
    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE build_core_hamiltonian_matrix

! *****************************************************************************

  SUBROUTINE distribute_s_matrix(globenv)

!   Purpose: Distribute the overlap matrix.

!   History: - Creation (26.06.2000, Matthias Krack)

!   ***************************************************************************

!   ij: The element (i,j) of an upper triangular matrix mapped to the
!       corresponding index of a packed vector.

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE distribute_s_matrix (MODULE core_hamiltonian)"

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: neighbor_list
    TYPE(neighbor_node_type), POINTER :: neighbor_node

    INTEGER :: group,iatom,ij,ikind,ipe,istat,jatom,jkind,jpe,mype,&
               nblock_full_matrix,npe,nsgfa,nsgfb,output_unit,sum_nblock_pe
    LOGICAL :: ionode

    INTEGER, DIMENSION(natom*(natom+1)/2)  :: nblock,nelement
    INTEGER, DIMENSION(0:globenv%num_pe-1) :: nblock_pe,nelement_pe

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos
    npe = globenv%num_pe
    output_unit = globenv%scr

    nblock(:) = 0
    nelement(:) = 0

!   *** Allocate the overlap matrix ***

    CALL allocate_matrix(matrix=s%matrix,&
                         nblock_row=natom,&
                         nblock_col=natom,&
                         nrow=nsgf,&
                         ncol=nsgf,&
                         first_row=atom_info(:)%first_sgf,&
                         last_row=atom_info(:)%last_sgf,&
                         first_col=atom_info(:)%first_sgf,&
                         last_col=atom_info(:)%last_sgf,&
                         matrix_name="OVERLAP MATRIX",&
                         matrix_symmetry="symmetric")

!   *** Loop over all pairs of adjacent atoms and distribute work ***

    neighbor_list => first_neighbor_list(orb_neighbor_list_set)

    DO WHILE (ASSOCIATED(neighbor_list))

      CALL get_neighbor_list(neighbor_list=neighbor_list,&
                             atom=iatom)

      ikind = atom_info(iatom)%kind
      nsgfa = kind_info(ikind)%orb_basis_set%nsgf

      neighbor_node => first_neighbor_node(neighbor_list)

      DO WHILE (ASSOCIATED(neighbor_node))

        CALL get_neighbor_node(neighbor_node=neighbor_node,&
                               neighbor=jatom)

        jkind = atom_info(jatom)%kind
        nsgfb = kind_info(jkind)%orb_basis_set%nsgf

!       *** Process only the upper triangle matrix ***

        IF (iatom <= jatom) THEN
          ij = iatom + jatom*(jatom - 1)/2
          nblock(ij) = nblock(ij) + 1
          nelement(ij) = nelement(ij) + nsgfa*nsgfb
        END IF

        neighbor_node => next_neighbor_node(neighbor_node)

      END DO

      neighbor_list => next_neighbor_list(neighbor_list)

    END DO

!   *** Distribute the atom blocks ***

    nblock_pe(:) = 0
    nelement_pe(:) = 0

    DO iatom=1,natom
      DO jatom=iatom,natom

        ij = iatom + jatom*(jatom - 1)/2

        IF (nelement(ij) > 0) THEN

          ipe = 0

          DO jpe=0,npe-1
            IF (nelement_pe(jpe) < nelement_pe(ipe)) ipe = jpe
          END DO

          nblock_pe(ipe) = nblock_pe(ipe) + nblock(ij)
          nelement_pe(ipe) = nelement_pe(ipe) + nelement(ij)

          IF (ipe == mype) CALL add_block_node(matrix=s%matrix,&
                                               block_row=iatom,&
                                               block_col=jatom)

        END IF

      END DO
    END DO

!   *** Print the distribution of the overlap matrix ***

    IF (globenv%print%distribution) THEN

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T2,A,/,/,T3,A,/,/,T5,A,/,/,(I6,8X,I8,8X,I10))")&
          "DISTRIBUTION OF THE OVERLAP MATRIX ELEMENTS",&
          "Image atoms included:",&
          "PE   Matrix blocks   Matrix elements",&
          (ipe,nblock_pe(ipe),nelement_pe(ipe),ipe=0,npe-1)
        WRITE (UNIT=output_unit,FMT="(/,T4,A3,8X,I8,8X,I10)")&
          "Sum",SUM(nblock_pe),SUM(nelement_pe)
      END IF

      nblock_pe(:) = 0
      nelement_pe(:) = 0

      CALL get_matrix_info(matrix=s%matrix,&
                           nblock_allocated=nblock_pe(mype),&
                           nelement_allocated=nelement_pe(mype))

      CALL mp_sum(nblock_pe,group)
      CALL mp_sum(nelement_pe,group)

      nblock_full_matrix = natom*(natom + 1)/2
      sum_nblock_pe = SUM(nblock_pe)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T3,A,/,/,T5,A,/,/,(I6,8X,I8,8X,I10))")&
          "Image atoms not included:",&
          "PE   Matrix blocks   Matrix elements",&
          (ipe,nblock_pe(ipe),nelement_pe(ipe),ipe=0,npe-1)
        WRITE (UNIT=output_unit,FMT="(/,T4,A3,8X,I8,8X,I10)")&
          "Sum",sum_nblock_pe,SUM(nelement_pe)
        WRITE (UNIT=output_unit,FMT="(/,T4,A3,8X,I8,A,F5.1,A)")&
          " of",nblock_full_matrix," blocks in the full matrix (",&
          100.0_wp*REAL(sum_nblock_pe,wp)/REAL(nblock_full_matrix,wp),&
          " % occupation)"
      END IF

    END IF

  END SUBROUTINE distribute_s_matrix

! *****************************************************************************

  SUBROUTINE build_h_matrix(globenv)

!   Purpose: Calculate the core Hamiltonian integral matrix.

!   History: - Creation (14.09.2000, Matthias Krack)

!   ***************************************************************************

    USE ai_coulomb,     ONLY: coulomb3
    USE ai_kinetic,     ONLY: kinetic
    USE ai_nuclear,     ONLY: nuclear
    USE ai_overlap,     ONLY: overlap
    USE ai_overlap_ppl, ONLY: overlap_ppl
    USE ai_verfc,       ONLY: verfc
    USE mathlib,        ONLY: symmetrize_matrix

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE build_h_matrix (MODULE core_hamiltonian)"

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: orb_neighbor_list,&
                                         verfc_neighbor_list,&
                                         vppl_neighbor_list

    REAL(wp) :: dab,dac,dbc,rab2,rac2,rbc2,radius_a,radius_b,radius_c,zetc
    INTEGER  :: first_sgfa,first_sgfb,i,&
                iatom,ikind,ineighbor,ipgf,iset,ishell,istat,&
                jatom,jkind,jneighbor,jpgf,jset,jshell,&
                katom,kkind,kkind_atom,&
                la_max,la_min,last_sgfa,last_sgfb,lb_max,lb_min,&
                nsgfa,nsgfb,ncoa,ncob,npgfa,npgfb,nseta,nsetb,nshella,nshellb

    REAL(wp), DIMENSION(1)            :: gccc
    REAL(wp), DIMENSION(3)            :: rab,rac,rbc
    INTEGER, DIMENSION(:), POINTER    :: orb_neighbors,verfc_neighbors,&
                                         vppl_neighbors
    REAL(wp), DIMENSION(:), POINTER   :: orb_r2,rpgfa,rpgfb,verfc_r2,vppl_r2,&
                                         zeta,zetb
    REAL(wp), DIMENSION(:,:), POINTER :: h_block,orb_r,s_block,t_block,&
                                         verfc_r,vppl_r

!   ---------------------------------------------------------------------------

    NULLIFY (orb_neighbor_list,orb_neighbors,orb_r,orb_r2)
    NULLIFY (verfc_neighbor_list,verfc_neighbors,verfc_r,verfc_r2)
    NULLIFY (vppl_neighbor_list,vppl_neighbors,vppl_r,vppl_r2)

!   *** Loop over all atoms ***

    orb_neighbor_list => first_neighbor_list(orb_neighbor_list_set)

    DO WHILE (ASSOCIATED(orb_neighbor_list))

      CALL get_neighbor_list(orb_neighbor_list,iatom)

      CALL extract_neighbor_list(orb_neighbor_list,orb_r2,orb_neighbors,orb_r)

      IF (gpw.AND.all_potential_present) THEN
        verfc_neighbor_list => find_neighbor_list(verfc_neighbor_list_set,&
                                                  iatom)
        IF (ASSOCIATED(verfc_neighbor_list)) THEN
          CALL extract_neighbor_list(verfc_neighbor_list,verfc_r2,&
                                     verfc_neighbors,verfc_r)
        END IF
      END IF

      IF (gth_potential_present) THEN
        vppl_neighbor_list =>&
          find_neighbor_list(vppl_neighbor_list_set,iatom)
        IF (ASSOCIATED(vppl_neighbor_list)) THEN
          CALL extract_neighbor_list(vppl_neighbor_list,vppl_r2,&
                                     vppl_neighbors,vppl_r)
        END IF
      END IF

      ikind = atom_info(iatom)%kind
      nseta = kind_info(ikind)%orb_basis_set%nset

!     *** Loop over all neighbor atoms of the current atom "iatom" ***

      DO ineighbor=1,SIZE(orb_neighbors)

        jatom = orb_neighbors(ineighbor)

!       *** Skip all atomic blocks in the lower triangle matrix ***

        IF (iatom > jatom) CYCLE

        rab2 = orb_r2(ineighbor)
        rab(:) = orb_r(:,ineighbor)
        dab = SQRT(rab2)

        CALL get_block_node(matrix=s%matrix,&
                            block_row=iatom,&
                            block_col=jatom,&
                            block=s_block)

!       *** Check, if the atomic block has to be ***
!       *** calculated by the current processor  ***

        IF (ASSOCIATED(s_block)) THEN

          IF (globenv%print%kinetic_energy_matrix) THEN
            CALL get_block_node(matrix=t%matrix,&
                                block_row=iatom,&
                                block_col=jatom,&
                                block=t_block)
          END IF

          CALL get_block_node(matrix=h%matrix,&
                              block_row=iatom,&
                              block_col=jatom,&
                              block=h_block)

          jkind = atom_info(jatom)%kind
          nsetb = kind_info(jkind)%orb_basis_set%nset

          DO iset=1,nseta

            la_max = kind_info(ikind)%orb_basis_set%lmax(iset)
            la_min = kind_info(ikind)%orb_basis_set%lmin(iset)

            npgfa = kind_info(ikind)%orb_basis_set%npgf(iset)
            nshella = kind_info(ikind)%orb_basis_set%nshell(iset)
            rpgfa => kind_info(ikind)%orb_basis_set%pgf_radius(1:npgfa,iset)
            zeta => kind_info(ikind)%orb_basis_set%zet(1:npgfa,iset)

            first_sgfa = kind_info(ikind)%orb_basis_set%first_sgf(1,iset)
            last_sgfa = kind_info(ikind)%orb_basis_set%last_sgf(nshella,iset)
            nsgfa = last_sgfa - first_sgfa + 1
            ncoa = npgfa*ncoset(la_max)

            IF (iatom == jatom) THEN
              i = iset
            ELSE
              i = 1
            END IF

            radius_a = kind_info(ikind)%orb_basis_set%set_radius(iset)

            DO jset=i,nsetb

              radius_b = kind_info(jkind)%orb_basis_set%set_radius(jset)

              IF (radius_a + radius_b < dab) CYCLE

              lb_max = kind_info(jkind)%orb_basis_set%lmax(jset)
              lb_min = kind_info(jkind)%orb_basis_set%lmin(jset)

              npgfb = kind_info(jkind)%orb_basis_set%npgf(jset)
              nshellb = kind_info(jkind)%orb_basis_set%nshell(jset)
              rpgfb => kind_info(jkind)%orb_basis_set%pgf_radius(1:npgfb,jset)
              zetb => kind_info(jkind)%orb_basis_set%zet(1:npgfb,jset)

              first_sgfb = kind_info(jkind)%orb_basis_set%first_sgf(1,jset)
              last_sgfb = kind_info(jkind)%orb_basis_set%last_sgf(nshellb,jset)
              nsgfb = last_sgfb - first_sgfb + 1
              ncob = npgfb*ncoset(lb_max)

!             *** Calculate the primitive overlap integrals ***

              CALL overlap(la_max,zeta,rpgfa,la_min,&
                           lb_max,zetb,rpgfb,lb_min,&
                           rab,rab2,sab)

!             *** Contraction step (overlap matrix) ***

              CALL dgemm("N","N",ncoa,nsgfb,ncob,1.0_wp,sab(1,1),SIZE(sab,1),&
                         kind_info(jkind)%orb_basis_set%sphi(1,first_sgfb),&
                         SIZE(kind_info(jkind)%orb_basis_set%sphi,1),0.0_wp,&
                         work(1,1),SIZE(work,1))

              CALL dgemm("T","N",nsgfa,nsgfb,ncoa,1.0_wp,&
                         kind_info(ikind)%orb_basis_set%sphi(1,first_sgfa),&
                         SIZE(kind_info(ikind)%orb_basis_set%sphi,1),&
                         work(1,1),SIZE(work,1),1.0_wp,&
                         s_block(first_sgfa,first_sgfb),SIZE(s_block,1))

!             *** Calculate the primitive kinetic energy integrals ***

              CALL kinetic(la_max,zeta,rpgfa,la_min,&
                           lb_max,zetb,rpgfb,lb_min,&
                           rab,rab2,sab,hab)

!             *** Contraction step (kinetic energy matrix, only for output) ***

              IF (globenv%print%kinetic_energy_matrix) THEN
                CALL dgemm("N","N",ncoa,nsgfb,ncob,1.0_wp,hab(1,1),SIZE(hab,1),&
                           kind_info(jkind)%orb_basis_set%sphi(1,first_sgfb),&
                           SIZE(kind_info(jkind)%orb_basis_set%sphi,1),0.0_wp,&
                           work(1,1),SIZE(work,1))
                CALL dgemm("T","N",nsgfa,nsgfb,ncoa,1.0_wp,&
                           kind_info(ikind)%orb_basis_set%sphi(1,first_sgfa),&
                           SIZE(kind_info(ikind)%orb_basis_set%sphi,1),&
                           work(1,1),SIZE(work,1),1.0_wp,&
                           t_block(first_sgfa,first_sgfb),SIZE(t_block,1))
              END IF

!             *** Calculate the nuclear attraction contribution ***

              IF (allchem) THEN

!               *** All-electron atoms ***

                DO kkind=1,nkind
                  IF (ASSOCIATED(kind_info(kkind)%all_potential)) THEN
                    DO kkind_atom=1,kind_info(kkind)%natom
                      katom = kind_info(kkind)%atom_list(kkind_atom)
                      rac(:) = atom_info(katom)%r(:) - atom_info(iatom)%r(:)
                      CALL nuclear(la_max,zeta,rpgfa,la_min,&
                                   lb_max,zetb,rpgfb,lb_min,&
                                   REAL(kind_info(kkind)%z,wp),&
                                   rab,rab2,rac,hab)
                    END DO
                  END IF
                END DO

!               *** GTH pseudopotential atoms ***

                IF (ASSOCIATED(vppl_neighbor_list)) THEN

                  DO jneighbor=1,SIZE(vppl_neighbors)

                    katom = vppl_neighbors(jneighbor)
                    kkind = atom_info(katom)%kind

                    radius_c = kind_info(kkind)%gth_potential%ppl_radius

                    rac2 = vppl_r2(jneighbor)
                    rac(:) = vppl_r(:,jneighbor)
                    dac = SQRT(rac2)

                    IF (radius_a + radius_c < dac) CYCLE

                    rbc(:) = rac(:) - rab(:)
                    rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                    dbc = SQRT(rbc2)

                    IF (radius_b + radius_c < dbc) CYCLE

                    zetc = kind_info(kkind)%gth_potential%alpha_ppl
                    gccc(1) = -kind_info(kkind)%gth_potential%cerf_ppl

                    CALL coulomb3(la_max,zeta,rpgfa,la_min,&
                                  lb_max,zetb,rpgfb,lb_min,&
                                  0,zetc,radius_c,0,gccc,&
                                  rab,rab2,rac,rac2,rbc,rbc2,hab)

                  END DO

                END IF

              ELSE

!               *** All-electron atoms ***

                IF (ASSOCIATED(verfc_neighbor_list)) THEN

                  DO jneighbor=1,SIZE(verfc_neighbors)

                    katom = verfc_neighbors(jneighbor)
                    kkind = atom_info(katom)%kind

                    radius_c = kind_info(kkind)%core_charge_radius

                    rac2 = verfc_r2(jneighbor)
                    rac(:) = verfc_r(:,jneighbor)
                    dac = SQRT(rac2)

                    IF (radius_a + radius_c < dac) CYCLE

                    rbc(:) = rac(:) - rab(:)
                    rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                    dbc = SQRT(rbc2)

                    IF (radius_b + radius_c < dbc) CYCLE

                    zetc = kind_info(kkind)%alpha_core_charge

                    CALL verfc(la_max,zeta,rpgfa,la_min,&
                               lb_max,zetb,rpgfb,lb_min,&
                               zetc,radius_c,&
                               REAL(kind_info(kkind)%z,wp),&
                               kind_info(kkind)%ccore_charge,&
                               rab,rab2,rac,rac2,rbc,rbc2,hab)

                  END DO

                END IF

              END IF

!             *** GTH pseudopotential atoms ***

              IF (ASSOCIATED(vppl_neighbor_list)) THEN

                DO jneighbor=1,SIZE(vppl_neighbors)

                  katom = vppl_neighbors(jneighbor)
                  kkind = atom_info(katom)%kind

                  radius_c = kind_info(kkind)%gth_potential%ppl_radius

                  rac2 = vppl_r2(jneighbor)
                  rac(:) = vppl_r(:,jneighbor)
                  dac = SQRT(rac2)

                  IF (radius_a + radius_c < dac) CYCLE

                  rbc(:) = rac(:) - rab(:)
                  rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                  dbc = SQRT(rbc2)

                  IF (radius_b + radius_c < dbc) CYCLE

                  zetc = kind_info(kkind)%gth_potential%alpha_ppl

                  CALL overlap_ppl(la_max,zeta,rpgfa,la_min,&
                                   lb_max,zetb,rpgfb,lb_min,&
                                   kind_info(kkind)%gth_potential%cexp_ppl,&
                                   zetc,radius_c,&
                                   rab,rab2,rac,rac2,rbc,rbc2,hab)
                END DO

              END IF

!             *** Contraction step (core Hamiltonian matrix) ***

              CALL dgemm("N","N",ncoa,nsgfb,ncob,1.0_wp,hab(1,1),SIZE(hab,1),&
                         kind_info(jkind)%orb_basis_set%sphi(1,first_sgfb),&
                         SIZE(kind_info(jkind)%orb_basis_set%sphi,1),0.0_wp,&
                         work(1,1),SIZE(work,1))

              CALL dgemm("T","N",nsgfa,nsgfb,ncoa,1.0_wp,&
                         kind_info(ikind)%orb_basis_set%sphi(1,first_sgfa),&
                         SIZE(kind_info(ikind)%orb_basis_set%sphi,1),&
                         work(1,1),SIZE(work,1),1.0_wp,&
                         h_block(first_sgfa,first_sgfb),SIZE(h_block,1))

            END DO

          END DO

!         *** Symmetrize the diagonal blocks ***

          IF (iatom == jatom) THEN
            CALL symmetrize_matrix(s_block,"upper_to_lower")
            IF (globenv%print%kinetic_energy_matrix) THEN
              CALL symmetrize_matrix(t_block,"upper_to_lower")
            END IF
            CALL symmetrize_matrix(h_block,"upper_to_lower")
          END IF

        END IF

      END DO

      orb_neighbor_list => next_neighbor_list(orb_neighbor_list)

    END DO

!   *** Add the V(non-local) contribution of the GTH pseudopotential ***

    IF (gth_potential_present) CALL add_ppnl_matrix()

!   *** Release work storage ***

    IF (ASSOCIATED(orb_neighbors)) THEN
      DEALLOCATE (orb_neighbors,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"orb_neighbors")
      END IF
    END IF

    IF (ASSOCIATED(orb_r)) THEN
      DEALLOCATE (orb_r,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"orb_r")
      END IF
    END IF

    IF (ASSOCIATED(orb_r2)) THEN
      DEALLOCATE (orb_r2,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"orb_r2")
      END IF
    END IF

    IF (ASSOCIATED(verfc_neighbors)) THEN
      DEALLOCATE (verfc_neighbors,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"verfc_neighbors")
      END IF
    END IF

    IF (ASSOCIATED(verfc_r)) THEN
      DEALLOCATE (verfc_r,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"verfc_r")
      END IF
    END IF

    IF (ASSOCIATED(verfc_r2)) THEN
      DEALLOCATE (verfc_r2,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"verfc_r2")
      END IF
    END IF

    IF (ASSOCIATED(vppl_neighbors)) THEN
      DEALLOCATE (vppl_neighbors,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"vppl_neighbors")
      END IF
    END IF

    IF (ASSOCIATED(vppl_r)) THEN
      DEALLOCATE (vppl_r,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"vppl_r")
      END IF
    END IF

    IF (ASSOCIATED(vppl_r2)) THEN
      DEALLOCATE (vppl_r2,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"vppl_r2")
      END IF
    END IF

  END SUBROUTINE build_h_matrix

! *****************************************************************************

  SUBROUTINE add_ppnl_matrix()

!   Purpose: Add the non-local contribution of the GTH pseudopotential to the
!            core Hamiltonian matrix.

!   History: - Creation (24.10.2000, Matthias Krack)

!   ***************************************************************************

    USE ai_overlap, ONLY: overlap

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE add_ppnl_matrix (MODULE core_hamiltonian)"

    TYPE(neighbor_list_type), POINTER   :: neighbor_list
    TYPE(real_block_node_type), POINTER :: block_node_ca,block_node_cb

    REAL(wp) :: dca,rca2,radius_a
    INTEGER  :: first_sgfa,first_prjc,iatom,iatom_prj_ppnl,ikind,ineighbor,&
                iset,istat,jatom,katom,kkind,l,la_max,la_min,last_sgfa,lc_max,&
                lc_min,natom_prj_ppnl,nsgfa,nsgfb,ncoa,ncoc,npgfa,nprj_ppnl,&
                nprjc,nseta,nshella

    REAL(wp), DIMENSION(1)            :: rprjc,zetc
    REAL(wp), DIMENSION(3)            :: rca
    REAL(wp), DIMENSION(:), POINTER   :: r2,rpgfa,zeta
    INTEGER, DIMENSION(:), POINTER    :: first_prj_ppnl,last_prj_ppnl,&
                                         neighbors
    REAL(wp), DIMENSION(:,:), POINTER :: r,h_block,s_block,sca,scb
    INTEGER, DIMENSION(:,:), POINTER  :: cell

!   ---------------------------------------------------------------------------

    NULLIFY (cell,first_prj_ppnl,last_prj_ppnl,neighbors,r,r2)

!   *** Get the number of operator atoms ***

    CALL get_neighbor_list_set_info(neighbor_list_set=&
                                      vppnl_ca_neighbor_list_set,&
                                    nlist=natom_prj_ppnl)

!   *** Prepare the allocation of the projector (ppnl) overlap matrix ***

    first_prj_ppnl => reallocate(first_prj_ppnl,1,natom_prj_ppnl)
    last_prj_ppnl => reallocate(last_prj_ppnl,1,natom_prj_ppnl)

    iatom_prj_ppnl = 0
    nprj_ppnl = 0

    neighbor_list => first_neighbor_list(vppnl_ca_neighbor_list_set)

    DO WHILE (ASSOCIATED(neighbor_list))
      iatom_prj_ppnl = iatom_prj_ppnl + 1
      CALL get_neighbor_list(neighbor_list=neighbor_list,atom=katom)
      kkind = atom_info(katom)%kind
      first_prj_ppnl(iatom_prj_ppnl) = nprj_ppnl + 1
      nprj_ppnl = nprj_ppnl + kind_info(kkind)%gth_potential%nppnl
      last_prj_ppnl(iatom_prj_ppnl) = nprj_ppnl
      neighbor_list => next_neighbor_list(neighbor_list)
    END DO

!   *** Allocate the projector (ppnl) overlap matrices ***

    CALL allocate_matrix(matrix=sca_prj_ppnl%matrix,&
                         nblock_row=natom_prj_ppnl,&
                         nblock_col=natom,&
                         nrow=nprj_ppnl,&
                         ncol=nsgf,&
                         first_row=first_prj_ppnl(:),&
                         last_row=last_prj_ppnl(:),&
                         first_col=atom_info(:)%first_sgf,&
                         last_col=atom_info(:)%last_sgf,&
                         matrix_name="CA PROJECTOR (PPNL) OVERLAP MATRIX",&
                         matrix_symmetry="no symmetry")

    CALL allocate_matrix(matrix=scb_prj_ppnl%matrix,&
                         nblock_row=natom_prj_ppnl,&
                         nblock_col=natom,&
                         nrow=nprj_ppnl,&
                         ncol=nsgf,&
                         first_row=first_prj_ppnl(:),&
                         last_row=last_prj_ppnl(:),&
                         first_col=atom_info(:)%first_sgf,&
                         last_col=atom_info(:)%last_sgf,&
                         matrix_name="CB PROJECTOR (PPNL) OVERLAP MATRIX",&
                         matrix_symmetry="no symmetry")

!   *** Calculate the overlap integrals between the orbital ***
!   *** basis functions and the GTH(nl) projector functions ***

    iatom_prj_ppnl = 0

    neighbor_list => first_neighbor_list(vppnl_ca_neighbor_list_set)

    DO WHILE (ASSOCIATED(neighbor_list))

      iatom_prj_ppnl = iatom_prj_ppnl + 1

!     *** Load the current neighbor list ***

      CALL get_neighbor_list(neighbor_list=neighbor_list,&
                             atom=katom)

      CALL extract_neighbor_list(neighbor_list,r2,neighbors,r,cell)

      kkind = atom_info(katom)%kind
      rprjc(1) = kind_info(kkind)%gth_potential%ppnl_radius

!     *** Loop over all neighbor atoms of the current operator atom "katom" ***

      DO ineighbor=1,SIZE(neighbors)

        iatom = neighbors(ineighbor)
        ikind = atom_info(iatom)%kind
        nseta = kind_info(ikind)%orb_basis_set%nset

!       *** Create a new matrix block ***

        NULLIFY (s_block)

        IF ((cell(1,ineighbor) == 0).AND.&
            (cell(2,ineighbor) == 0).AND.&
            (cell(3,ineighbor) == 0)) THEN
          CALL add_block_node(matrix=sca_prj_ppnl%matrix,&
                              block_row=iatom_prj_ppnl,&
                              block_col=iatom,&
                              block=s_block)
        ELSE
          CALL add_block_node(matrix=scb_prj_ppnl%matrix,&
                              block_row=iatom_prj_ppnl,&
                              block_col=iatom,&
                              block=s_block)
        END IF

!       *** Calculate the current overlap matrix block ***

        rca2 = r2(ineighbor)
        rca(:) = r(:,ineighbor)

        dca = SQRT(rca2)

        first_prjc = 1

        DO l=0,kind_info(kkind)%gth_potential%lppnl

          nprjc = kind_info(kkind)%gth_potential%nprj_ppnl(l)*nco(l)

          IF (nprjc == 0) CYCLE

          lc_max = l + 2*(kind_info(kkind)%gth_potential%nprj_ppnl(l) - 1)
          lc_min = l
          zetc(1) = kind_info(kkind)%gth_potential%alpha_ppnl(l)
          ncoc = ncoset(lc_max)

          DO iset=1,nseta

            radius_a = kind_info(ikind)%orb_basis_set%set_radius(iset)

            IF (rprjc(1) + radius_a < dca) CYCLE

            la_max = kind_info(ikind)%orb_basis_set%lmax(iset)
            la_min = kind_info(ikind)%orb_basis_set%lmin(iset)

            npgfa = kind_info(ikind)%orb_basis_set%npgf(iset)
            nshella = kind_info(ikind)%orb_basis_set%nshell(iset)
            rpgfa => kind_info(ikind)%orb_basis_set%pgf_radius(1:npgfa,iset)
            zeta => kind_info(ikind)%orb_basis_set%zet(1:npgfa,iset)

            first_sgfa = kind_info(ikind)%orb_basis_set%first_sgf(1,iset)
            last_sgfa = kind_info(ikind)%orb_basis_set%last_sgf(nshella,iset)
            nsgfa = last_sgfa - first_sgfa + 1
            ncoa = npgfa*ncoset(la_max)

!           *** Calculate the primitive overlap integrals ***

            CALL overlap(lc_max,zetc,rprjc,lc_min,&
                         la_max,zeta,rpgfa,la_min,&
                         rca,rca2,sab)

!           *** Contraction step (orbital basis functions) ***

            CALL dgemm("N","N",ncoc,nsgfa,ncoa,1.0_wp,sab(1,1),SIZE(sab,1),&
                       kind_info(ikind)%orb_basis_set%sphi(1,first_sgfa),&
                       SIZE(kind_info(ikind)%orb_basis_set%sphi,1),0.0_wp,&
                       work(1,1),SIZE(work,1))

!           *** Contraction step (projector functions) ***

            CALL dgemm("T","N",nprjc,nsgfa,ncoc,1.0_wp,&
                       kind_info(kkind)%gth_potential%cprj(1,first_prjc),&
                       SIZE(kind_info(kkind)%gth_potential%cprj,1),&
                       work(1,1),SIZE(work,1),1.0_wp,&
                       s_block(first_prjc,first_sgfa),SIZE(s_block,1))

          END DO

          first_prjc = first_prjc + nprjc

        END DO

      END DO

!     *** Add the V(non-local) contribution to the core Hamiltonian matrix ***

      block_node_ca => first_block_node(sca_prj_ppnl%matrix,iatom_prj_ppnl)

      DO WHILE (ASSOCIATED(block_node_ca))

        CALL get_block_node(block_node=block_node_ca,&
                            block_col=iatom,&
                            block=sca)

        nprjc = SIZE(sca,1)
        nsgfa = SIZE(sca,2)

        CALL dgemm("T","N",nprjc,nsgfa,nprjc,1.0_wp,&
                   kind_info(kkind)%gth_potential%vprj_ppnl(1,1),nprjc,&
                   sca(1,1),nprjc,0.0_wp,work(1,1),SIZE(work,1))

        block_node_cb => block_node_ca

        DO WHILE (ASSOCIATED(block_node_cb))

          CALL get_block_node(block_node=block_node_cb,&
                              block_col=jatom,&
                              block=scb)

          nsgfb = SIZE(scb,2)

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
            CALL dgemm("T","N",nsgfa,nsgfb,nprjc,1.0_wp,work(1,1),&
                       SIZE(work,1),scb(1,1),nprjc,1.0_wp,h_block(1,1),nsgfa)
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
            CALL dgemm("T","N",nsgfb,nsgfa,nprjc,1.0_wp,scb(1,1),nprjc,&
                       work(1,1),SIZE(work,1),1.0_wp,h_block(1,1),nsgfb)
          END IF

          block_node_cb => next_block_node(block_node_cb)

        END DO

        block_node_cb => first_block_node(scb_prj_ppnl%matrix,iatom_prj_ppnl)

        DO WHILE (ASSOCIATED(block_node_cb))

          CALL get_block_node(block_node=block_node_cb,&
                              block_col=jatom,&
                              block=scb)

          IF (jatom >= iatom) THEN
            nsgfb = SIZE(scb,2)
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
            CALL dgemm("T","N",nsgfa,nsgfb,nprjc,1.0_wp,work(1,1),&
                       SIZE(work,1),scb(1,1),nprjc,1.0_wp,h_block(1,1),nsgfa)
          END IF

          block_node_cb => next_block_node(block_node_cb)

        END DO

        block_node_ca => next_block_node(block_node_ca)

      END DO

!     *** Release work storage ***

      CALL deallocate_matrix_row(sca_prj_ppnl%matrix,iatom_prj_ppnl)
      CALL deallocate_matrix_row(scb_prj_ppnl%matrix,iatom_prj_ppnl)

      neighbor_list => next_neighbor_list(neighbor_list)

    END DO

!   *** Release work storage ***

    CALL deallocate_matrix(sca_prj_ppnl%matrix)
    CALL deallocate_matrix(scb_prj_ppnl%matrix)

    IF (ASSOCIATED(r2)) THEN
      DEALLOCATE (r2,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"r2")
      END IF
    END IF

    IF (ASSOCIATED(neighbors)) THEN
      DEALLOCATE (neighbors,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"neighbors")
      END IF
    END IF

    IF (ASSOCIATED(r)) THEN
      DEALLOCATE (r,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"r")
      END IF
    END IF

    IF (ASSOCIATED(cell)) THEN
      DEALLOCATE (cell,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"cell")
      END IF
    END IF

    IF (ASSOCIATED(first_prj_ppnl)) THEN
      DEALLOCATE (first_prj_ppnl,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"first_prj_ppnl")
      END IF
    END IF

    IF (ASSOCIATED(last_prj_ppnl)) THEN
      DEALLOCATE (last_prj_ppnl,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,"last_prj_ppnl")
      END IF
    END IF

  END SUBROUTINE add_ppnl_matrix

! *****************************************************************************

END MODULE core_hamiltonian
