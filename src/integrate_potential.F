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
!!     Joost VandeVondele (02.2002)
!!           1) rewrote collocate_pgf for increased accuracy and speed
!!           2) collocate_core hack for PGI compiler
!!           3) added multiple grid feature
!!
!!   SOURCE
!******************************************************************************
!#include "f_hpm.h"
MODULE integrate_potential

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE atomic_kinds,        ONLY: kind_info,nkind
  USE atoms,               ONLY: atom_info
  USE basis_set_types,     ONLY: gto_basis_set_type,maxco,maxsgf,maxsgf_set
  USE coefficient_types,   ONLY: coeff_type,coeff_transform_space,coeff_copy
  USE gaussian_gridlevels, ONLY: gaussian_gridlevel,gridlevel_info_type
  USE global_types,        ONLY: global_environment_type
  USE mathlib,             ONLY: symmetrize_matrix
  USE matrix_types,        ONLY: get_block_node,&
                                 real_matrix_set_type
  USE memory_utilities,    ONLY: reallocate
  USE message_passing,     ONLY: mp_sum
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
  USE cube_utils,          ONLY: cube_info_type,return_cube
  USE l_utils,             ONLY: l_info_type,return_l_info

  IMPLICIT NONE

  PRIVATE

  INTEGER :: npgf_product = 0

  PUBLIC :: integrate_v_rspace

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE integrate_v_rspace(p,v_rspace,v_gspace,gridlevel_info,v_rspace_old,  &
                                   cube_info,l_info,rebuild_gvg,gvg,h,nproduct, & 
                                   globenv)

    TYPE(coeff_type), DIMENSION(:), POINTER  :: v_rspace,v_gspace
    TYPE(coeff_type), INTENT(INOUT)           :: v_rspace_old ! old will contain new
    TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
    TYPE(gridlevel_info_type), INTENT(IN) :: gridlevel_info

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_set_type), INTENT(INOUT) :: gvg,h
    INTEGER, INTENT(OUT)                      :: nproduct
    TYPE(real_matrix_set_type), INTENT(IN)    :: p
    LOGICAL                                   :: rebuild_gvg             
    TYPE(l_info_type), INTENT(IN)             :: l_info

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: basis_set_a,basis_set_b
    TYPE(neighbor_list_type), POINTER :: orb_neighbor_list

    REAL(wp) :: dab,rab2,radius_kind_a,radius_kind_b,radius_set_a,radius_set_b,scaling_factor
    INTEGER  :: first_isgf,first_jset,first_jsgf,first_sgfa,&
                first_sgfb,handle,iatom,ico,ikind,ineighbor,ipgf,iset,isgf,&
                isgfa,isgfb,ishell,istat,jatom,jco,jkind,jneighbor,jpgf,jset,&
                jsgf,jshell,la_max,la_min,last_sgfa,last_sgfb,lb_max,lb_min,&
                na1,na2,nb1,nb2,nsgfa,nsgfb,ncoa,ncob,npgfa,npgfb,nseta,nsetb,&
                nshella,nshellb,iatom_old,jatom_old,igrid_level

    REAL(wp), DIMENSION(3) :: ra,rab

    INTEGER, DIMENSION(:), POINTER    :: orb_neighbors
    REAL(wp), DIMENSION(:), POINTER   :: orb_r2,rpgfa,rpgfb,zeta,zetb
    REAL(wp), DIMENSION(:,:), POINTER :: gvg_block,gvgab,orb_r,s_block,&
                                         work,p_block,pab,h_block
    REAL(wp) :: maxvnew,maxvdiff,vratio,old_gvg_factor
!   ---------------------------------------------------------------------------

    CALL timeset("integrate_v_rspace","I","",handle)

    npgf_product = 0

!   *** get difference potential ***
    maxvnew=MAXVAL(ABS(v_rspace(1)%pw%cr3d))
    if (.not.rebuild_gvg) then
      v_rspace(1)%pw%cr3d=v_rspace(1)%pw%cr3d-v_rspace_old%pw%cr3d
    endif
    maxvdiff=MAXVAL(ABS(v_rspace(1)%pw%cr3d))
    vratio=maxvdiff/maxvnew 

!   *** get the potential on the subgrids in real space, via fft
    IF (gridlevel_info%ngrid_levels .gt. 1) then
       CALL coeff_transform_space(v_rspace(1),v_gspace(1)) ! we have the gspace potential
       DO igrid_level=2,gridlevel_info%ngrid_levels
          CALL coeff_copy(v_gspace(1),v_gspace(igrid_level)) 
          CALL coeff_transform_space(v_gspace(igrid_level),v_rspace(igrid_level))
          !*** multiply by the grid volume element ratio ...
          scaling_factor=v_rspace(igrid_level)%pw%pw_grid%dvol/v_rspace(1)%pw%pw_grid%dvol
          v_rspace(igrid_level)%pw%cr3d=v_rspace(igrid_level)%pw%cr3d*scaling_factor
       ENDDO
    ENDIF

!   *** Allocate work storage ***

    NULLIFY (pab,gvgab,work)

    gvgab => reallocate(gvgab,1,maxco,1,maxco)
    pab => reallocate(pab,1,maxco,1,maxco)
    work => reallocate(work,1,maxco,1,maxsgf_set)

!   *** incremental stuff ! goal is to add gvg to h or zero it
    iatom_old=0
    jatom_old=0
    NULLIFY (orb_neighbor_list,orb_neighbors,orb_r,orb_r2)
    orb_neighbor_list => first_neighbor_list(orb_neighbor_list_set)
    DO WHILE (ASSOCIATED(orb_neighbor_list))
      CALL get_neighbor_list(orb_neighbor_list,iatom)
      CALL extract_neighbor_list(orb_neighbor_list,orb_r2,orb_neighbors,orb_r)
      DO ineighbor=1,SIZE(orb_neighbors)
        jatom = orb_neighbors(ineighbor)
        IF (iatom > jatom) CYCLE

        ! each iatom,jatom pair should appear only once (e.g. image atoms)
        IF (iatom .eq. iatom_old .and. jatom .eq. jatom_old) cycle 
        iatom_old=iatom
        jatom_old=jatom

        CALL get_block_node(matrix=gvg%matrix,&
                            block_row=iatom,&
                            block_col=jatom,&
                            block=gvg_block)

        CALL get_block_node(matrix=h%matrix,&
                            block_row=iatom,&
                            block_col=jatom,&
                            block=h_block)
        IF (.NOT.ASSOCIATED(gvg_block) .and. .NOT.ASSOCIATED(h_block)) CYCLE

        if (rebuild_gvg) then
           gvg_block(:,:)=0.0_wp
        else
           h_block(:,:)=h_block(:,:)+gvg_block(:,:)
        endif

      END DO
      orb_neighbor_list => next_neighbor_list(orb_neighbor_list)
    END DO

!   *** Loop over all atoms ***

    NULLIFY (orb_neighbor_list,orb_neighbors,orb_r,orb_r2)

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

        CALL get_block_node(matrix=gvg%matrix,&
                            block_row=iatom,&
                            block_col=jatom,&
                            block=gvg_block)
        IF (.NOT.ASSOCIATED(gvg_block)) CYCLE

        CALL get_block_node(matrix=p%matrix,&
                            block_row=iatom,&
                            block_col=jatom,&
                            block=p_block)
        IF (.NOT.ASSOCIATED(p_block)) CYCLE

        CALL get_block_node(matrix=h%matrix,&
                            block_row=iatom,&
                            block_col=jatom,&
                            block=h_block)
        IF (.NOT.ASSOCIATED(h_block)) CYCLE

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

!           *** Decontract the density matrix set block ***

            CALL dgemm("N","N",ncoa,nsgfb,nsgfa,1.0_wp,&
                       basis_set_a%sphi(1,first_sgfa),&
                       SIZE(basis_set_a%sphi,1),&
                       p_block(first_sgfa,first_sgfb),SIZE(p_block,1),&
                       0.0_wp,work(1,1),SIZE(work,1))
            CALL dgemm("N","T",ncoa,ncob,nsgfb,1.0_wp,&
                       work(1,1),SIZE(work,1),basis_set_b%sphi(1,first_sgfb),&
                       SIZE(basis_set_b%sphi,1),0.0_wp,pab(1,1),SIZE(pab,1))


            DO jco=1,ncob
              DO ico=1,ncoa
                gvgab(ico,jco) = 0.0_wp
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

                igrid_level=gaussian_gridlevel(gridlevel_info,zeta(ipgf)+zetb(jpgf))
                CALL integrate_pgf_product_rspace(la_max,zeta(ipgf),la_min,&
                                                  lb_max,zetb(jpgf),lb_min,&
                                                  ra,rab,rab2,v_rspace(igrid_level)%pw,&
                                                  cube_info(igrid_level),l_info,&
                                                  gvgab(na1:na2,nb1:nb2),&
                                                  pab(na1:na2,nb1:nb2),vratio,rebuild_gvg)

              END DO

            END DO

!           *** Contract the operator matrix set block ***

            CALL dgemm("N","N",ncoa,nsgfb,ncob,1.0_wp,gvgab(1,1),SIZE(gvgab,1),&
                       basis_set_b%sphi(1,first_sgfb),&
                       SIZE(basis_set_b%sphi,1),0.0_wp,work(1,1),SIZE(work,1))
            CALL dgemm("T","N",nsgfa,nsgfb,ncoa,1.0_wp,&
                       basis_set_a%sphi(1,first_sgfa),&
                       SIZE(basis_set_a%sphi,1),work(1,1),&
                       SIZE(work,1),1.0_wp,&
                       gvg_block(first_sgfa,first_sgfb),SIZE(gvg_block,1))

            CALL dgemm("T","N",nsgfa,nsgfb,ncoa,1.0_wp,&
                       basis_set_a%sphi(1,first_sgfa),&
                       SIZE(basis_set_a%sphi,1),work(1,1),&
                       SIZE(work,1),1.0_wp,&
                       h_block(first_sgfa,first_sgfb),SIZE(h_block,1))

          END DO

        END DO

        IF (iatom == jatom) CALL symmetrize_matrix(h_block,"upper_to_lower")

      END DO

      orb_neighbor_list => next_neighbor_list(orb_neighbor_list)

    END DO

!   *** Release work storage ***

    DEALLOCATE (pab,gvgab,work)

    CALL mp_sum(npgf_product,globenv%group)

    nproduct = npgf_product

! reconstruct the potential
    if (.not.rebuild_gvg) then
      v_rspace(1)%pw%cr3d=v_rspace(1)%pw%cr3d+v_rspace_old%pw%cr3d
    endif
    v_rspace_old%pw%cr3d=v_rspace(1)%pw%cr3d

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE integrate_v_rspace

! *****************************************************************************

    SUBROUTINE integrate_pgf_product_rspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,pw,cube_info,l_info,hab,pab,&
                                          vratio,rebuild_gvg)

    USE interactions, ONLY: eps_gvg,exp_radius,exp_radius_very_extended

    TYPE(pw_type), TARGET, INTENT(IN)  :: pw
    REAL(wp), INTENT(IN)               :: rab2,zeta,zetb,vratio
    INTEGER, INTENT(IN)                :: la_max,la_min,lb_max,lb_min
    REAL(wp), DIMENSION(3), INTENT(IN) :: ra,rab
    TYPE(cube_info_type),INTENT(IN)    :: cube_info
    TYPE(l_info_type),INTENT(IN)       :: l_info
    LOGICAL                            :: rebuild_gvg

    REAL(wp), DIMENSION(ncoset(la_max),ncoset(lb_max)), INTENT(INOUT) :: hab
    REAL(wp), DIMENSION(ncoset(la_max),ncoset(lb_max)), INTENT(IN) :: pab

!   *** Local variables ***

    REAL(wp) :: f,pij,prefactor,radius,yz,z,zetp
    INTEGER  :: ax,ay,az,bx,by,bz,handle,i,ico,ig,j,jco,jg,k,kg,la,lb,&
                lb_cube_min,ub_cube_max,gridbounds(2,3)

    REAL(wp), DIMENSION(3) :: dr,rap,rbp,roffset,rp,rb
    INTEGER, DIMENSION(3)  :: cubecenter,lb_cube,ng,ub_cube,ub_grid

    INTEGER, DIMENSION(:,:), POINTER      :: map
    REAL(wp), DIMENSION(:,:,:), POINTER   :: grid
    INTEGER lxyz_max,lxy_max,lx_max,lx,lxb,lya,lyb,lza,lzb,lz,ly,icoef,lxy
    INTEGER coef_max,lxa,lxyz,l
    REAL(wp) a,b,binomial_k_lxa,binomial_l_lxb,zbp,zap,za,zb,ya,yb,yap,ybp,rpg,s,pg
    REAL(wp)  :: xa,xb,xap,xbp,radius2,maxpab
    REAL(wp) :: dx2,dy2,dz2,dxi,dyi,dzi,cutoff
    INTEGER  :: kgmax,kgmin,jgmax,jgmin,igmax,igmin,j2,jg2,lar(3),lbr(3), &
                start,length,offset


    INTEGER, pointer, dimension(:) :: ly_max,lz_max
    REAL(wp), allocatable :: pzyx(:), alpha(:,:)
    REAL(wp), allocatable :: polz(:,:),dpz(:,:),poly(:,:),dpy(:,:),polx(:,:)
    REAL(wp) , allocatable :: dpx(:,:)
    INTEGER, pointer, dimension(:)  :: sphere_bounds

!   ---------------------------------------------------------------------------

    coef_max=la_max+lb_max+1
    SELECT CASE(coef_max)
    CASE(1)
      CALL timeset("integrate_pgf_rspace_1","I","",handle)
    CASE(2)
      CALL timeset("integrate_pgf_rspace_2","I","",handle)
    CASE(3)
      CALL timeset("integrate_pgf_rspace_3","I","",handle)
    CASE(4)
      CALL timeset("integrate_pgf_rspace_4","I","",handle)
    CASE(5)
      CALL timeset("integrate_pgf_rspace_5","I","",handle)
    CASE DEFAULT
      CALL timeset("integrate_pgf_rspace_SLOW","I","",handle)
    END SELECT
    zetp = zeta + zetb
    f = zetb/zetp
    prefactor = EXP(-zeta*f*rab2)
!   *** position of the gaussian product
    rap(:) = f*rab(:)
    rbp(:) = rap(:) - rab(:)
    rp(:) = ra(:) + rap(:)  ! this is the gaussian center in real coordinates
    rb(:) = ra(:) + rab(:)


    ! the cutoff could be vratio instead of 1.0_wp if not rebuild_gvg
    ! in this case the minimization goes fine up to the point where the energy is 
    ! noisy
    cutoff=1.0_wp
    radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,ra,rb,rp,zetp,&
                                eps_gvg,prefactor*vratio,cutoff)

    IF (radius == 0.0_wp) THEN
      CALL timestop(0.0_wp,handle)
      RETURN
    END IF

    
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)
    CALL return_l_info(l_info,la_min,la_max,lb_min,lb_max,lx_max,lxy_max, &
                         lxyz_max,ly_max,lz_max)


!   *** irrelevant ? ***
    npgf_product = npgf_product + (ncoset(la_max) - ncoset(la_min-1))*&
                                  (ncoset(lb_max) - ncoset(lb_min-1))

!   *** properties of the grid ***
    dr(:) = pw%pw_grid%dr(:)
    ng(:) = pw%pw_grid%npts(:)
    ub_grid(:) = pw%pw_grid%bounds(2,:)
    grid => pw%cr3d(:,:,:)
    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:) = rp(:) - REAL(cubecenter(:),wp)*dr(:)
    lb_cube_min = MINVAL(lb_cube(:))
    ub_cube_max = MAXVAL(ub_cube(:))

    NULLIFY(map)
    map => reallocate(map,lb_cube_min,ub_cube_max,1,3)
    allocate(pzyx(lxyz_max), &
             polz(lxyz_max,lb_cube(3):ub_cube(3)), &
             dpz(0:lb_max,0:la_max), &
             poly(lxy_max,lb_cube(2):ub_cube(2)), &
             dpy(0:lb_max,0:la_max), &
             alpha(coef_max,lx_max), &
             polx(coef_max,lb_cube(1):ub_cube(1)))

!   *** a mapping so that the ig corresponds to the right grid point, also with pbc
    DO i=1,3
        start=lb_cube(i)
        DO
         offset=MODULO(cubecenter(i)+start,ng(i))+1-start
         length=MIN(ub_cube(i),ng(i)-offset)-start
         DO ig=start,start+length
          map(ig,i) = ig+offset
         END DO
         if (start+length.ge.ub_cube(i)) EXIT
         start=start+length+1
        END DO
    ENDDO

!   *** initialise the pol x,y,z terms
    DO ig=lb_cube(3),ub_cube(3)
      lxyz=0
      rpg = REAL(ig,wp)*dr(3) - roffset(3)
      zap = EXP(-zetp*rpg**2)*prefactor  ! take into account the prefactor here
      za  = rpg + rap(3)
      zb  = za  - rab(3)
      DO lza=0,la_max
       zbp=1.0_wp
       DO lzb=0,lb_max
          dpz(lzb,lza)=zap*zbp
          zbp=zbp*zb
       ENDDO
       zap=zap*za
      ENDDO
      DO lxa=0,la_max
      DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
          DO lza=max(la_min-lxa-lya,0),la_max-lxa-lya
          DO lzb=max(lb_min-lxb-lyb,0),lb_max-lxb-lyb
             lxyz=lxyz+1
             polz(lxyz,ig)=dpz(lzb,lza)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
      ENDDO
      ENDDO 
    ENDDO

    DO ig=lb_cube(2),ub_cube(2)
      rpg = REAL(ig,wp)*dr(2) - roffset(2)
      yap = EXP(-zetp*rpg**2)
      ya  = rpg + rap(2)
      yb  = ya  - rab(2)
      DO lya=0,la_max
       ybp=1.0_wp
       DO lyb=0,lb_max
          dpy(lyb,lya)=yap*ybp
          ybp=ybp*yb
       ENDDO
       yap=yap*ya
      ENDDO

      lxy=0
      DO lxa=0,la_max
      DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
             lxy=lxy+1
             poly(lxy,ig)=dpy(lyb,lya)
       ENDDO
       ENDDO
      ENDDO
      ENDDO
    ENDDO

 
!   *** make the alpha matrix ***
    alpha(:,:)=0.0_wp
    lx=0
    DO lxa=0,la_max
    DO lxb=0,lb_max
       lx=lx+1
       binomial_k_lxa=1.0_wp
       a=1.0_wp
       DO k=0,lxa
        binomial_l_lxb=1.0_wp
        b=1.0_wp
        DO l=0,lxb
           alpha(lxa-l+lxb-k+1,lx)=alpha(lxa-l+lxb-k+1,lx)+ &
                             binomial_k_lxa*binomial_l_lxb*a*b
           binomial_l_lxb=binomial_l_lxb*real(lxb-l,wp)/real(l+1,wp)
           b=b*(rp(1)-(ra(1)+rab(1)))
        ENDDO
        binomial_k_lxa=binomial_k_lxa*real(lxa-k,wp)/real(k+1,wp)
        a=a*(rp(1)-ra(1))
       ENDDO
    ENDDO
    ENDDO


    DO ig=lb_cube(1),ub_cube(1)
      rpg = REAL(ig,wp)*dr(1) - roffset(1)
      pg  = EXP(-zetp*rpg**2)
      DO icoef=1,coef_max
         polx(icoef,ig)=pg
         pg=pg*(rpg)
      ENDDO
    ENDDO

!   *** do the loop over the grid

    call integrate_core(polx,poly,polz,grid,alpha, &
                        lx_max,lxy_max,lxyz_max,coef_max,lb_cube,ub_cube, &
                        ly_max,lz_max,gridbounds,lb_cube_min,ub_cube_max,map,hab, &
                        la_min,la_max,lb_min,lb_max,sphere_bounds)

    DEALLOCATE (map,pzyx,alpha,polz,dpz,dpy,poly,polx)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE integrate_pgf_product_rspace

  SUBROUTINE integrate_core(polx,poly,polz,grid,alpha,lx_max,lxy_max,lxyz_max,&
                            coef_max,lb_cube,ub_cube,ly_max,lz_max,gridbounds,&
                            lb_cube_min,ub_cube_max,map,hab,la_min,la_max,&
                            lb_min,lb_max,sphere_bounds)
    integer,  INTENT(IN) :: sphere_bounds(:)
    integer,  INTENT(IN) :: lx_max,lxy_max,lxyz_max,coef_max,&
                            lb_cube(3),ub_cube(3),la_min,la_max,lb_min,lb_max
    integer,  INTENT(IN) :: ly_max(lx_max),lz_max(lxy_max)
    integer,  INTENT(IN) :: lb_cube_min,ub_cube_max,map(lb_cube_min:ub_cube_max,1:3)
    integer,  INTENT(IN) :: gridbounds(2,3)
    real(wp), INTENT(IN) :: polz(lxyz_max,lb_cube(3):ub_cube(3))
    real(wp), INTENT(IN) :: poly(lxy_max,lb_cube(2):ub_cube(2))
    real(wp), INTENT(IN) :: polx(coef_max,lb_cube(1):ub_cube(1))
    real(wp), INTENT(IN) :: alpha(coef_max,lx_max)
! help the compiler not to copy the array, but allow optimization on other machines
#ifdef __PGI
    real(wp), DIMENSION(:,:,:), POINTER :: grid
#else
    real(wp), INTENT(IN) :: grid(gridbounds(1,1):gridbounds(2,1),&
                            gridbounds(1,2):gridbounds(2,2),&
                            gridbounds(1,3):gridbounds(2,3)) 
#endif
    REAL(wp), DIMENSION(ncoset(la_max),ncoset(lb_max)), INTENT(INOUT) :: hab
    real(wp) coef(coef_max,4)
    real(wp) pyx(4,lxy_max),pzyx(lxyz_max),radius2,dz2,dy2

    integer kg,kgmax,k,jgmax,jgmin,jg,j,j2,igmax,ig,i,kgmin,igmin,k2
    integer lx,ly,lz,lxy,lxyz,icoef,lxa,lxb,lya,lyb,lza,lzb,ico,jco,i2,sci
    real(wp) s1,s2,s3,s4,s5,s6,s7,s8,s9,sa,g1,g2
    real(wp) s1k,s2k,s3k,s4k,s5k,s6k,s7k,s8k,s9k,sak,g1k,g2k
    real(wp) r1,r2,r3,r4,r5,r6

!    call f_hpmstart(21,"gvg")
  
    pzyx(:)=0.0_wp
    sci=1
    kgmin=sphere_bounds(sci)
    sci=sci+1
    kgmax=sphere_bounds(sci)
    sci=sci+1

    DO kg=kgmin,kgmax,2
       k=map(kg,3)
       k2=map(kg+1,3)
       pyx(:,:)=0.0_wp
       jgmin=sphere_bounds(sci)
       sci=sci+1
       jgmax=sphere_bounds(sci)
       sci=sci+1
       DO jg=jgmin,jgmax,2
          j=map(jg,2)
          j2=map(jg+1,2)
          igmin=sphere_bounds(sci)
          sci=sci+1
          igmax=sphere_bounds(sci)
          sci=sci+1
          SELECT CASE (coef_max)
          CASE (1)
             s1=0.0_wp
             s2=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             DO ig=igmin,igmax,2
                i =map(ig,1)
                i2 =map(ig+1,1)
                r1=polx(1,ig)
                r2=polx(1,ig+1)
                s1=s1+grid(i,j,k)     *r1
                s2=s2+grid(i,j2,k)    *r1
                s1k=s1k+grid(i,j,k2)  *r1
                s2k=s2k+grid(i,j2,k2) *r1
                s3=s3+grid(i2,j,k)     *r2
                s4=s4+grid(i2,j2,k)    *r2
                s3k=s3k+grid(i2,j,k2)  *r2
                s4k=s4k+grid(i2,j2,k2) *r2
             END DO
             r3=alpha(1,1)
             r1=poly(1,jg)*r3
             r2=poly(1,jg+1)*r3
             pyx(1,1)=pyx(1,1)+(s1+s3)*r1
             pyx(2,1)=pyx(2,1)+(s2+s4)*r2
             pyx(3,1)=pyx(3,1)+(s1k+s3k)*r1
             pyx(4,1)=pyx(4,1)+(s2k+s4k)*r2
          CASE (2)
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             DO ig=igmin,igmax
                i=map(ig,1)
                r1=grid(i,j,k)
                r2=grid(i,j2,k)
                r3=grid(i,j,k2)
                r4=grid(i,j2,k2)
                r5=polx(1,ig)
                r6=polx(2,ig)
                s1=s1+r1*r5
                s3=s3+r1*r6
                s2=s2+r2*r5
                s4=s4+r2*r6
                s1k=s1k+r3*r5
                s3k=s3k+r3*r6
                s2k=s2k+r4*r5
                s4k=s4k+r4*r6
             END DO
             lxy=0
             DO lx=1,lx_max
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              g1=r1*s1+r2*s3
              g2=r1*s2+r2*s4
              g1k=r1*s1k+r2*s3k
              g2k=r1*s2k+r2*s4k
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=poly(lxy,jg)
                 r2=poly(lxy,jg+1)
                 pyx(1,lxy)=pyx(1,lxy)+g1*r1
                 pyx(2,lxy)=pyx(2,lxy)+g2*r2
                 pyx(3,lxy)=pyx(3,lxy)+g1k*r1
                 pyx(4,lxy)=pyx(4,lxy)+g2k*r2
              ENDDO
             ENDDO
          CASE (3)
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s5=0.0_wp
             s6=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             s5k=0.0_wp
             s6k=0.0_wp
             DO ig=igmin,igmax
                i=map(ig,1)
                r1=grid(i,j,k)
                r2=grid(i,j2,k)
                r5=polx(1,ig)
                s1=s1+r1*r5
                s2=s2+r2*r5
                r3=grid(i,j,k2)
                r4=grid(i,j2,k2)
                s1k=s1k+r3*r5
                s2k=s2k+r4*r5
                r5=polx(2,ig)
                s3=s3+r1*r5
                s4=s4+r2*r5
                s3k=s3k+r3*r5
                s4k=s4k+r4*r5
                r5=polx(3,ig)
                s5=s5+r1*r5
                s6=s6+r2*r5
                s5k=s5k+r3*r5
                s6k=s6k+r4*r5
             END DO
             lxy=0
             DO lx=1,lx_max
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              r3=alpha(3,lx)
              g1=r1*s1+r2*s3+r3*s5
              g2=r1*s2+r2*s4+r3*s6
              g1k=r1*s1k+r2*s3k+r3*s5k
              g2k=r1*s2k+r2*s4k+r3*s6k
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=poly(lxy,jg)
                 r2=poly(lxy,jg+1)
                 pyx(1,lxy)=pyx(1,lxy)+g1*r1
                 pyx(2,lxy)=pyx(2,lxy)+g2*r2
                 pyx(3,lxy)=pyx(3,lxy)+g1k*r1
                 pyx(4,lxy)=pyx(4,lxy)+g2k*r2
              ENDDO
             ENDDO
          CASE (4)
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s5=0.0_wp
             s6=0.0_wp
             s7=0.0_wp
             s8=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             s5k=0.0_wp
             s6k=0.0_wp
             s7k=0.0_wp
             s8k=0.0_wp
             DO ig=igmin,igmax
                i=map(ig,1)
                s1=s1+grid(i,j ,k)*polx(1,ig)
                s2=s2+grid(i,j2,k)*polx(1,ig)
                s3=s3+grid(i,j ,k)*polx(2,ig)
                s4=s4+grid(i,j2,k)*polx(2,ig)
                s5=s5+grid(i,j ,k)*polx(3,ig)
                s6=s6+grid(i,j2,k)*polx(3,ig)
                s7=s7+grid(i,j ,k)*polx(4,ig)
                s8=s8+grid(i,j2,k)*polx(4,ig)
                s1k=s1k+grid(i,j ,k2)*polx(1,ig)
                s2k=s2k+grid(i,j2,k2)*polx(1,ig)
                s3k=s3k+grid(i,j ,k2)*polx(2,ig)
                s4k=s4k+grid(i,j2,k2)*polx(2,ig)
                s5k=s5k+grid(i,j ,k2)*polx(3,ig)
                s6k=s6k+grid(i,j2,k2)*polx(3,ig)
                s7k=s7k+grid(i,j ,k2)*polx(4,ig)
                s8k=s8k+grid(i,j2,k2)*polx(4,ig)
             END DO
             lxy=0
             DO lx=1,lx_max
              g1=alpha(1,lx)*s1+alpha(2,lx)*s3+alpha(3,lx)*s5+alpha(4,lx)*s7
              g2=alpha(1,lx)*s2+alpha(2,lx)*s4+alpha(3,lx)*s6+alpha(4,lx)*s8
              g1k=alpha(1,lx)*s1k+alpha(2,lx)*s3k+alpha(3,lx)*s5k+alpha(4,lx)*s7k
              g2k=alpha(1,lx)*s2k+alpha(2,lx)*s4k+alpha(3,lx)*s6k+alpha(4,lx)*s8k
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 pyx(1,lxy)=pyx(1,lxy)+g1*poly(lxy,jg)
                 pyx(2,lxy)=pyx(2,lxy)+g2*poly(lxy,jg+1)
                 pyx(3,lxy)=pyx(3,lxy)+g1k*poly(lxy,jg)
                 pyx(4,lxy)=pyx(4,lxy)+g2k*poly(lxy,jg+1)
              ENDDO
             ENDDO
           CASE (5)
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s5=0.0_wp
             s6=0.0_wp
             s7=0.0_wp
             s8=0.0_wp
             s9=0.0_wp
             sa=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             s5k=0.0_wp
             s6k=0.0_wp
             s7k=0.0_wp
             s8k=0.0_wp
             s9k=0.0_wp
             sak=0.0_wp
             DO ig=igmin,igmax
                i=map(ig,1)
                g1=grid(i,j,k)
                g2=grid(i,j2,k)
                g1k=grid(i,j,k2)
                g2k=grid(i,j2,k2)
                r1=polx(1,ig)
                r2=polx(2,ig)
                r3=polx(3,ig)
                r4=polx(4,ig)
                r5=polx(5,ig)
                s1=s1+g1*r1
                s2=s2+g2*r1
                s3=s3+g1*r2
                s4=s4+g2*r2
                s5=s5+g1*r3
                s6=s6+g2*r3
                s7=s7+g1*r4
                s8=s8+g2*r4
                s9=s9+g1*r5
                sa=sa+g2*r5
                s1k=s1k+g1k*r1
                s2k=s2k+g2k*r1
                s3k=s3k+g1k*r2
                s4k=s4k+g2k*r2
                s5k=s5k+g1k*r3
                s6k=s6k+g2k*r3
                s7k=s7k+g1k*r4
                s8k=s8k+g2k*r4
                s9k=s9k+g1k*r5
                sak=sak+g2k*r5
             END DO
             lxy=0
             DO lx=1,lx_max
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              r3=alpha(3,lx)
              r4=alpha(4,lx)
              r5=alpha(5,lx)
              g1=r1*s1
              g2=r1*s2
              g1k=r1*s1k
              g2k=r1*s2k
              g1=g1+r2*s3
              g2=g2+r2*s4
              g1k=g1k+r2*s3k
              g2k=g2k+r2*s4k
              g1=g1+r3*s5
              g2=g2+r3*s6
              g1k=g1k+r3*s5k
              g2k=g2k+r3*s6k
              g1=g1+r4*s7
              g2=g2+r4*s8
              g1k=g1k+r4*s7k
              g2k=g2k+r4*s8k
              g1=g1+r5*s9
              g2=g2+r5*sa
              g1k=g1k+r5*s9k
              g2k=g2k+r5*sak
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=poly(lxy,jg)
                 r2=poly(lxy,jg+1)
                 pyx(1,lxy)=pyx(1,lxy)+g1*r1
                 pyx(2,lxy)=pyx(2,lxy)+g2*r2
                 pyx(3,lxy)=pyx(3,lxy)+g1k*r1
                 pyx(4,lxy)=pyx(4,lxy)+g2k*r2
              ENDDO
             ENDDO
          CASE DEFAULT
             coef(:,:)=0.0_wp
             DO ig=igmin,igmax
                i=map(ig,1)
                DO icoef=1,coef_max  
                  coef(icoef,1)=coef(icoef,1)+grid(i,j,k)*polx(icoef,ig)
                  coef(icoef,2)=coef(icoef,2)+grid(i,j2,k)*polx(icoef,ig)
                  coef(icoef,3)=coef(icoef,3)+grid(i,j,k2)*polx(icoef,ig)
                  coef(icoef,4)=coef(icoef,4)+grid(i,j2,k2)*polx(icoef,ig)
                ENDDO
             END DO
             lxy=0
             DO lx=1,lx_max
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO icoef=1,coef_max
                  g1=g1+alpha(icoef,lx)*coef(icoef,1)   
                  g2=g2+alpha(icoef,lx)*coef(icoef,2)   
                  g1k=g1k+alpha(icoef,lx)*coef(icoef,3)   
                  g2k=g2k+alpha(icoef,lx)*coef(icoef,4)   
              ENDDO
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 pyx(1,lxy)=pyx(1,lxy)+g1*poly(lxy,jg)
                 pyx(2,lxy)=pyx(2,lxy)+g2*poly(lxy,jg+1)
                 pyx(3,lxy)=pyx(3,lxy)+g1k*poly(lxy,jg)
                 pyx(4,lxy)=pyx(4,lxy)+g2k*poly(lxy,jg+1)
              ENDDO
             ENDDO
          END SELECT
       END DO ! loop over y, jg
       lxyz=0
       DO lxy=1,lxy_max
          s1=pyx(1,lxy)+pyx(2,lxy)
          s1k=pyx(3,lxy)+pyx(4,lxy)
          DO lz=1,lz_max(lxy)
             lxyz=lxyz+1
             pzyx(lxyz)=pzyx(lxyz)+s1*polz(lxyz,kg)+s1k*polz(lxyz,kg+1)
          ENDDO
       ENDDO
    END DO ! loop over z, kg

    lxyz=0
    DO lxa=0,la_max
    DO lxb=0,lb_max
       DO lya=0,la_max-lxa
       DO lyb=0,lb_max-lxb
          DO lza=max(la_min-lxa-lya,0),la_max-lxa-lya
          DO lzb=max(lb_min-lxb-lyb,0),lb_max-lxb-lyb
             lxyz=lxyz+1
             ico=coset(lxa,lya,lza)
             jco=coset(lxb,lyb,lzb)
             hab(ico,jco)=hab(ico,jco)+pzyx(lxyz)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO

!    call f_hpmstop(21)

  END SUBROUTINE


END MODULE integrate_potential
