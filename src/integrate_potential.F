!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/integrate_potential [1.0] *
!!
!!   NAME
!!     integrate_potential
!!
!!   FUNCTION
!!     Build up the plane wave density by collocating the primitive Gaussian
!!     functions (pgf).
!!
!!   AUTHOR
!!     Matthias Krack (03.04.2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE integrate_potential

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE atomic_kinds,        ONLY: kind_info,nkind
  USE atoms,               ONLY: atom_info
  USE basis_set_types,     ONLY: gto_basis_set_type,maxco,maxsgf,maxsgf_set
  USE coefficient_types,   ONLY: coeff_type
  USE mathlib,             ONLY: symmetrize_matrix
  USE matrix_types,        ONLY: get_block_node,&
                                 real_matrix_set_type
  USE memory_utilities,    ONLY: reallocate
  USE neighbor_list_types, ONLY: extract_neighbor_list,&
                                 find_neighbor_list,&
                                 first_neighbor_list,&
                                 first_neighbor_node,&
                                 get_neighbor_list,&
                                 get_neighbor_node,&
                                 neighbor_list_type,&
                                 neighbor_node_type,&
                                 next_neighbor_list,&
                                 next_neighbor_node
  USE nl,                  ONLY: orb_neighbor_list_set
  USE orbital_pointers,    ONLY: coset,indco,ncoset
  USE pw_types,            ONLY: pw_type
  USE timings,             ONLY: timeset,timestop

  IMPLICIT NONE

  PRIVATE

  INTEGER :: npgf_product = 0

  PUBLIC :: integrate_v_rspace

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE integrate_v_rspace(v_rspace,h,nproduct)

    TYPE(coeff_type), INTENT(IN)              :: v_rspace
    TYPE(real_matrix_set_type), INTENT(INOUT) :: h
    INTEGER, OPTIONAL, INTENT(OUT)            :: nproduct

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: basis_set_a,basis_set_b
    TYPE(neighbor_list_type), POINTER :: orb_neighbor_list

    REAL(wp) :: dab,rab2,radius_kind_a,radius_kind_b,radius_set_a,radius_set_b
    INTEGER  :: first_isgf,first_jset,first_jsgf,first_sgfa,&
                first_sgfb,handle,iatom,ico,ikind,ineighbor,ipgf,iset,isgf,&
                isgfa,isgfb,ishell,istat,jatom,jco,jkind,jneighbor,jpgf,jset,&
                jsgf,jshell,la_max,la_min,last_sgfa,last_sgfb,lb_max,lb_min,&
                na1,na2,nb1,nb2,nsgfa,nsgfb,ncoa,ncob,npgfa,npgfb,nseta,nsetb,&
                nshella,nshellb

    REAL(wp), DIMENSION(3) :: ra,rab

    INTEGER, DIMENSION(:), POINTER    :: orb_neighbors
    REAL(wp), DIMENSION(:), POINTER   :: orb_r2,rpgfa,rpgfb,zeta,zetb
    REAL(wp), DIMENSION(:,:), POINTER :: h_block,hab,orb_r,work

!   ---------------------------------------------------------------------------

    CALL timeset("integrate_v_rspace","I","",handle)

    npgf_product = 0

!   *** Allocate work storage ***

    NULLIFY (hab,work)

    hab => reallocate(hab,1,maxco,1,maxco)
    work => reallocate(work,1,maxco,1,maxsgf_set)

    NULLIFY (orb_neighbor_list,orb_neighbors,orb_r,orb_r2)

!   *** Loop over all atoms ***

    orb_neighbor_list => first_neighbor_list(orb_neighbor_list_set)

    DO WHILE (ASSOCIATED(orb_neighbor_list))

      CALL get_neighbor_list(orb_neighbor_list,iatom)

      CALL extract_neighbor_list(orb_neighbor_list,orb_r2,orb_neighbors,orb_r)

      ikind = atom_info(iatom)%kind
      first_isgf = atom_info(iatom)%first_sgf
      ra(:) = atom_info(iatom)%r_pbc(:)
      basis_set_a => kind_info(ikind)%orb_basis_set
      nseta = basis_set_a%nset
      radius_kind_a = basis_set_a%kind_radius

!     *** Loop over all neighbor atoms of the current atom "iatom" ***

      DO ineighbor=1,SIZE(orb_neighbors)

        jatom = orb_neighbors(ineighbor)

!       *** Skip all atomic blocks in the lower triangle matrix ***

        IF (iatom > jatom) CYCLE

        jkind = atom_info(jatom)%kind
        first_jsgf = atom_info(jatom)%first_sgf
        rab2 = orb_r2(ineighbor)
        rab(:) = orb_r(:,ineighbor)
        dab = SQRT(rab2)
        basis_set_b => kind_info(jkind)%orb_basis_set
        nsetb = basis_set_b%nset
        radius_kind_b = basis_set_b%kind_radius

        CALL get_block_node(matrix=h%matrix,&
                            block_row=iatom,&
                            block_col=jatom,&
                            block=h_block)

        DO iset=1,nseta

          radius_set_a = basis_set_a%set_radius(iset)

          IF (radius_set_a + radius_kind_b < dab) CYCLE

          la_max = basis_set_a%lmax(iset)
          la_min = basis_set_a%lmin(iset)
          npgfa = basis_set_a%npgf(iset)
          nshella = basis_set_a%nshell(iset)
          rpgfa => basis_set_a%pgf_radius(1:npgfa,iset)
          zeta => basis_set_a%zet(1:npgfa,iset)

          first_sgfa = basis_set_a%first_sgf(1,iset)
          last_sgfa = basis_set_a%last_sgf(nshella,iset)
          nsgfa = last_sgfa - first_sgfa + 1
          ncoa = npgfa*ncoset(la_max)

          IF (iatom == jatom) THEN
            first_jset = iset
          ELSE
            first_jset = 1
          END IF

          DO jset=first_jset,nsetb

            radius_set_b = basis_set_b%set_radius(jset)

            IF (radius_set_a + radius_set_b < dab) CYCLE

            lb_max = basis_set_b%lmax(jset)
            lb_min = basis_set_b%lmin(jset)
            npgfb = basis_set_b%npgf(jset)
            nshellb = basis_set_b%nshell(jset)
            rpgfb => basis_set_b%pgf_radius(1:npgfb,jset)
            zetb => basis_set_b%zet(1:npgfb,jset)
            first_sgfb = basis_set_b%first_sgf(1,jset)
            last_sgfb = basis_set_b%last_sgf(nshellb,jset)
            nsgfb = last_sgfb - first_sgfb + 1
            ncob = npgfb*ncoset(lb_max)

            DO jco=1,ncob
              DO ico=1,ncoa
                hab(ico,jco) = 0.0_wp
              END DO
            END DO

            DO ipgf=1,npgfa

              IF (rpgfa(ipgf) + radius_set_b < dab) CYCLE

              na1 = (ipgf - 1)*ncoset(la_max) + 1
              na2 = ipgf*ncoset(la_max)

              DO jpgf=1,npgfb

                IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) CYCLE

                nb1 = (jpgf - 1)*ncoset(lb_max) + 1
                nb2 = jpgf*ncoset(lb_max)

                CALL integrate_pgf_product_rspace(la_max,zeta(ipgf),la_min,&
                                                  lb_max,zetb(jpgf),lb_min,&
                                                  ra,rab,rab2,v_rspace%pw,&
                                                  hab(na1:na2,nb1:nb2))

              END DO

            END DO

!           *** Contract the operator matrix set block ***

            CALL dgemm("N","N",ncoa,nsgfb,ncob,1.0_wp,hab(1,1),SIZE(hab,1),&
                       basis_set_b%sphi(1,first_sgfb),&
                       SIZE(basis_set_b%sphi,1),0.0_wp,work(1,1),SIZE(work,1))
            CALL dgemm("T","N",nsgfa,nsgfb,ncoa,1.0_wp,&
                       basis_set_a%sphi(1,first_sgfa),&
                       SIZE(basis_set_a%sphi,1),work(1,1),SIZE(work,1),1.0_wp,&
                       h_block(first_sgfa,first_sgfb),SIZE(h_block,1))

          END DO

        END DO

        IF (iatom == jatom) CALL symmetrize_matrix(h_block,"upper_to_lower")

      END DO

      orb_neighbor_list => next_neighbor_list(orb_neighbor_list)

    END DO

    IF (PRESENT(nproduct)) nproduct = npgf_product

    DEALLOCATE (hab,work)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE integrate_v_rspace

! *****************************************************************************

  SUBROUTINE integrate_pgf_product_rspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,pw,hab)

    USE interactions, ONLY: eps_rho_rspace,exp_radius

    TYPE(pw_type), TARGET, INTENT(IN)  :: pw
    REAL(wp), INTENT(IN)               :: rab2,zeta,zetb
    INTEGER, INTENT(IN)                :: la_max,la_min,lb_max,lb_min
    REAL(wp), DIMENSION(3), INTENT(IN) :: ra,rab

    REAL(wp), DIMENSION(ncoset(la_max),ncoset(lb_max)), INTENT(INOUT) :: hab

!   *** Local variables ***

    REAL(wp) :: f,prefactor,radius,yz,z,zetp
    INTEGER  :: ax,ay,az,bx,by,bz,handle,i,ico,ig,j,jco,jg,k,kg,la,lb,&
                lb_cube_min,ub_cube_max

    REAL(wp), DIMENSION(3) :: rap,rbp,roffset,rp
    INTEGER, DIMENSION(3)  :: cubecenter,lb_cube,ub_cube

    REAL(wp), DIMENSION(:), POINTER       :: dr,rag,rbg,rpg
    INTEGER, DIMENSION(:), POINTER        :: ng
    INTEGER, DIMENSION(:,:), POINTER      :: map
    REAL(wp), DIMENSION(:,:,:), POINTER   :: grid
    REAL(wp), DIMENSION(:,:,:,:), POINTER :: cubeaxis

!   ---------------------------------------------------------------------------

    CALL timeset("integrate_pgf_product_rspace","I","",handle)

    zetp = zeta + zetb
    f = zetb/zetp
    prefactor = EXP(-zeta*f*rab2)

    radius = exp_radius(0,zetp,eps_rho_rspace,prefactor)

    IF (radius == 0.0_wp) THEN
      CALL timestop(0.0_wp,handle)
      RETURN
    END IF

    npgf_product = npgf_product + (ncoset(la_max) - ncoset(la_min-1))*&
                                  (ncoset(lb_max) - ncoset(lb_min-1))

    dr => pw%pw_grid%dr(:)
    ng => pw%pw_grid%npts(:)
    grid => pw%cr3d(:,:,:)

    rap(:) = f*rab(:)
    rbp(:) = rap(:) - rab(:)
    rp(:) = ra(:) + rap(:)
    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:) = rp(:) - REAL(cubecenter(:),wp)*dr(:)
    lb_cube(:) = FLOOR(-radius/dr(:))
    ub_cube(:) = CEILING(radius/dr(:))

!   *** Calculate the Gaussian function values for each cube axis ***

    lb_cube_min = MINVAL(lb_cube(:))
    ub_cube_max = MAXVAL(ub_cube(:))

    NULLIFY (cubeaxis,map,rag,rbg,rpg)

    cubeaxis => reallocate(cubeaxis,lb_cube_min,ub_cube_max,1,3,0,la_max,&
                           0,lb_max)
    map => reallocate(map,lb_cube_min,ub_cube_max,1,3)
    rag => reallocate(rag,lb_cube_min,ub_cube_max)
    rbg => reallocate(rbg,lb_cube_min,ub_cube_max)
    rpg => reallocate(rpg,lb_cube_min,ub_cube_max)

    DO i=1,3

      DO ig=lb_cube(i),ub_cube(i)
        map(ig,i) = MODULO(cubecenter(i) + ig,ng(i)) + 1
        rpg(ig) = REAL(ig,wp)*dr(i) - roffset(i)
        cubeaxis(ig,i,0,0) = EXP(-zetp*rpg(ig)**2)
      END DO

      IF (la_max > 0) THEN
        DO ig=lb_cube(i),ub_cube(i)
          rag(ig) = rap(i) + rpg(ig)
          cubeaxis(ig,i,1,0) = rag(ig)*cubeaxis(ig,i,0,0)
        END DO
        DO la=2,la_max
          DO ig=lb_cube(i),ub_cube(i)
            cubeaxis(ig,i,la,0) = rag(ig)*cubeaxis(ig,i,la-1,0)
          END DO
        END DO
        IF (lb_max > 0) THEN
          DO ig=lb_cube(i),ub_cube(i)
            rbg(ig) = rag(ig) - rab(i)
          END DO
          DO la=0,la_max
            DO lb=1,lb_max
              DO ig=lb_cube(i),ub_cube(i)
                cubeaxis(ig,i,la,lb) = rbg(ig)*cubeaxis(ig,i,la,lb-1)
              END DO
            END DO
          END DO
        END IF
      ELSE
        IF (lb_max > 0) THEN
          DO ig=lb_cube(i),ub_cube(i)
            rbg(ig) = rbp(i) + rpg(ig)
            cubeaxis(ig,i,0,1) = rbg(ig)*cubeaxis(ig,i,0,0)
          END DO
          DO lb=2,lb_max
            DO ig=lb_cube(i),ub_cube(i)
              cubeaxis(ig,i,0,lb) = rbg(ig)*cubeaxis(ig,i,0,lb-1)
            END DO
          END DO
        END IF
      END IF

    END DO

!   *** Add the current Gaussian function product to grid ***

    DO ico=ncoset(la_min-1)+1,ncoset(la_max)

      ax = indco(1,ico)
      ay = indco(2,ico)
      az = indco(3,ico)

      DO jco=ncoset(lb_min-1)+1,ncoset(lb_max)

        bx = indco(1,jco)
        by = indco(2,jco)
        bz = indco(3,jco)

        DO kg=lb_cube(3),ub_cube(3)
          k = map(kg,3)
          z = prefactor*cubeaxis(kg,3,az,bz)
          IF (ABS(z) < eps_rho_rspace) CYCLE
          DO jg=lb_cube(2),ub_cube(2)
            j = map(jg,2)
            yz = cubeaxis(jg,2,ay,by)*z
            IF (ABS(yz) < eps_rho_rspace) CYCLE
            DO ig=lb_cube(1),ub_cube(1)
              i = map(ig,1)
              hab(ico,jco) = hab(ico,jco) + grid(i,j,k)*cubeaxis(ig,1,ax,bx)*yz
            END DO
          END DO
        END DO

      END DO

    END DO

    DEALLOCATE (cubeaxis,map,rag,rbg,rpg)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE integrate_pgf_product_rspace

! *****************************************************************************

END MODULE integrate_potential
