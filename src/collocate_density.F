!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/collocate_density [1.0] *
!!
!!   NAME
!!     collocate_density
!!
!!   FUNCTION
!!     Calculate the plane wave density by collocating the primitive Gaussian
!!     functions (pgf).
!!
!!   AUTHOR
!!     Matthias Krack (03.04.2001)
!!     1) Joost VandeVondele (01.2002)
!!
!!   MODIFICATION HISTORY
!!     1) rewrote collocate for increased accuracy and speed
!!     2) introduced the PGI hack for increased speed with that compiler (22.02.02)
!!     3) Added Multiple Grid feature
!!     4) new way to get over the grid (01.03.02)
!!     5) removed timing calls since they were getting expensive
!!
!!   SOURCE
!******************************************************************************
!#include "f_hpm.h"
MODULE collocate_density

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE atomic_kinds,        ONLY: kind_info,nkind
  USE atoms,               ONLY: atom_info
  USE basis_set_types,     ONLY: gto_basis_set_type,maxco,maxsgf,maxsgf_set
  USE coefficient_types,   ONLY: coeff_type,coeff_transform_space,coeff_add,coeff_zero,coeff_copy
  USE global_types,        ONLY: global_environment_type
  USE gaussian_gridlevels, ONLY: gaussian_gridlevel,gridlevel_info_type
  USE interactions,        ONLY: eps_rho_gspace,eps_rho_rspace,exp_radius,&
                                 exp_radius_very_extended
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
  USE pw_types,            ONLY: pw_type, pw_copy, pw_allocate, pw_deallocate, COMPLEXDATA3D
  USE timings,             ONLY: timeset,timestop
  USE cube_utils,          ONLY: return_cube,cube_info_type
  USE l_utils,          ONLY: return_l_info,l_info_type
  USE termination,         ONLY: stop_program

  IMPLICIT NONE

  PRIVATE

  INTEGER :: npgf_product = 0

! *** Public subroutines ***

  PUBLIC :: calculate_rho_core,&
            calculate_rho_elec

! *** Public functions ***

  PUBLIC :: calculate_total_rho

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE calculate_rho_core(rho,cube_info,l_info,total_rho,globenv)

    TYPE(coeff_type), INTENT(INOUT)           :: rho
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(cube_info_type), INTENT(IN)          :: cube_info
    TYPE(l_info_type), INTENT(IN)          :: l_info
    REAL(wp), INTENT(OUT)                     :: total_rho

!   *** Local variables ***

    REAL(wp) :: alpha
    INTEGER  :: handle,iatom,ikind

    REAL(wp), DIMENSION(3)   :: ra
    REAL(wp), DIMENSION(1,1) :: pab

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_rho_core","I","",handle)

    DO ikind=1,nkind

      pab(1,1) = kind_info(ikind)%ccore_charge
      alpha = kind_info(ikind)%alpha_core_charge

      DO iatom=1,kind_info(ikind)%natom

        ra(:) = atom_info(kind_info(ikind)%atom_list(iatom))%r_pbc(:)

        IF (ASSOCIATED(rho%pw%cc3d)) THEN
          CALL collocate_pgf_product_gspace(0,alpha,0,0,0.0_wp,0,ra,&
                                            (/0.0_wp,0.0_wp,0.0_wp/),&
                                            0.0_wp,-1.0_wp,pab,rho%pw)
        ELSE IF (ASSOCIATED(rho%pw%cr3d)) THEN
          CALL collocate_pgf_product_rspace(0,alpha,0,0,0.0_wp,0,ra,&
                                            (/0.0_wp,0.0_wp,0.0_wp/),&
                                            0.0_wp,-1.0_wp,pab,rho%pw,cube_info,l_info)

        END IF

      END DO

    END DO

!   IF (ASSOCIATED(rho%pw%cc3d)) THEN
!     CALL mp_sum(rho%pw%cc3d,globenv%group)
!   ELSE IF (ASSOCIATED(rho%pw%cr3d)) THEN
!     CALL mp_sum(rho%pw%cr3d,globenv%group)
!   END IF

    total_rho = calculate_total_rho(rho)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_rho_core

! *****************************************************************************

  SUBROUTINE calculate_rho_elec(p,p_old,rebuild_rho,rho,rho_rspace_old,rho_gspace,&
       gridlevel_info, &
       cube_info,l_info, &
       total_rho,nproduct,globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(coeff_type), pointer, DIMENSION(:)   :: rho
    TYPE(coeff_type), pointer, DIMENSION(:)   :: rho_rspace_old
    TYPE(coeff_type), pointer, DIMENSION(:)   :: rho_gspace
    TYPE(gridlevel_info_type), INTENT(IN)     :: gridlevel_info
    TYPE(cube_info_type), pointer, DIMENSION(:)       :: cube_info
    TYPE(l_info_type), INTENT(IN)             :: l_info
    TYPE(real_matrix_set_type), INTENT(IN)    :: p,p_old
    LOGICAL, INTENT(IN)                       :: rebuild_rho
    REAL(wp), INTENT(OUT)                     :: total_rho
    INTEGER, INTENT(OUT)                      :: nproduct

    !   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: basis_set_a, basis_set_b
    TYPE(neighbor_list_type), POINTER :: orb_neighbor_list

    REAL(wp) :: radius_kind_a, radius_set_a, rab2

    INTEGER  :: handle, iatom, ineighbor, jatom, nseta
         
    REAL(wp), DIMENSION(3) :: ra, rab

    INTEGER, DIMENSION(:), POINTER    :: orb_neighbors
    INTEGER                           :: igrid_level
    REAL(wp), DIMENSION(:), POINTER   :: orb_r2
    REAL(wp), DIMENSION(:,:), POINTER :: orb_r, pab, pab_old, work, p_old_block

    !   ---------------------------------------------------------------------------

    CALL timeset("calculate_rho_elec","I","",handle)

    !   get the right density to start with

    DO igrid_level=1,gridlevel_info%ngrid_levels
       if (rebuild_rho) then
          CALL coeff_zero(rho(igrid_level))
       else
          CALL coeff_copy(rho_rspace_old(igrid_level),rho(igrid_level))
       endif
    ENDDO

    npgf_product = 0

    !   *** Allocate work storage ***

    NULLIFY (pab,work)
    pab => reallocate(pab,1,maxco,1,maxco)
    work => reallocate(work,1,maxco,1,maxsgf_set)

    if (.not.rebuild_rho) then 
       NULLIFY (pab_old)
       pab_old => reallocate(pab_old,1,maxco,1,maxco)
    endif

    NULLIFY (orb_neighbor_list,orb_neighbors,orb_r,orb_r2)

    !   *** Loop over all atoms ***

    orb_neighbor_list => first_neighbor_list(orb_neighbor_list_set)

    DO WHILE (ASSOCIATED(orb_neighbor_list))

       CALL get_neighbor_list(orb_neighbor_list,iatom)

       CALL extract_neighbor_list(orb_neighbor_list,orb_r2,orb_neighbors,orb_r)

       ra(:) = atom_info(iatom)%r_pbc(:)
       basis_set_a => kind_info(atom_info(iatom)%kind)%orb_basis_set
       nseta = basis_set_a%nset
       radius_kind_a = basis_set_a%kind_radius


       !     *** Loop over all neighbor atoms of the current atom "iatom" ***
!$omp parallel do &
!$omp default ( private ) &
!$omp shared ( p, p_old, pab, pab_old, work ) &
       DO ineighbor=1,SIZE(orb_neighbors)

          jatom = orb_neighbors(ineighbor)

          ! *** Skip all atomic blocks in the lower triangle matrix ***
          IF (iatom > jatom) CYCLE

          rab(:) = orb_r(:,ineighbor)
          rab2 = orb_r2(ineighbor)

          CALL ineighbor_loop_work ( &
               ineighbor, iatom, jatom, &
               p, p_old, rho, rebuild_rho, &
               basis_set_a, nseta, radius_kind_a, &
               ra, rab, rab2, &
               gridlevel_info, cube_info, l_info, &
               pab, pab_old, work )


       END DO
!$omp end parallel do

       orb_neighbor_list => next_neighbor_list(orb_neighbor_list)

    ENDDO


    !   *** storage ***
    DO igrid_level=1,gridlevel_info%ngrid_levels
       CALL coeff_copy(rho(igrid_level),rho_rspace_old(igrid_level))
    ENDDO

    !   *** get the density on the highest cutoff grid


    if ((gridlevel_info%ngrid_levels.gt.1) .or. (rho_gspace(1)%pw%in_use/=COMPLEXDATA3D )) then
       CALL coeff_transform_space(rho(1),rho_gspace(1))
    endif
    do igrid_level=2,gridlevel_info%ngrid_levels
       CALL coeff_transform_space(rho(igrid_level),rho_gspace(igrid_level))
       CALL coeff_add(rho_gspace(1),rho_gspace(igrid_level),rho_gspace(1))
    enddo
    ! the latter is needed for a consistent density in rho and g space
    if ((gridlevel_info%ngrid_levels.gt.1) .or. (rho_gspace(1)%pw%in_use/=COMPLEXDATA3D )) then
       CALL coeff_transform_space(rho_gspace(1),rho(1))
    endif

    !   *** Release work storage ***

    DEALLOCATE (pab,work)

    if (.not.rebuild_rho) then
       deallocate(pab_old) 
    endif

    CALL mp_sum(npgf_product,globenv%group)

    nproduct = npgf_product

    IF (ASSOCIATED(rho(1)%pw%cc3d)) THEN
       CALL mp_sum(rho(1)%pw%cc3d,globenv%group)
    ELSE IF (ASSOCIATED(rho(1)%pw%cr3d)) THEN
       CALL mp_sum(rho(1)%pw%cr3d,globenv%group)
    END IF

    total_rho = calculate_total_rho(rho(1))

    CALL timestop(0.0_wp,handle)


  END SUBROUTINE calculate_rho_elec

! *****************************************************************************

  FUNCTION calculate_total_rho(rho) RESULT(total_rho)

    TYPE(coeff_type), TARGET, INTENT(IN) :: rho

    REAL(wp) :: total_rho

!   *** Local variables ***

    INTEGER, DIMENSION(:), POINTER :: lb_grid

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(rho%pw%cc3d)) THEN
       lb_grid => rho%pw%pw_grid%bounds(1,:)
       total_rho = -rho%pw%pw_grid%vol*rho%pw%cc3d(lb_grid(1),&
                                                   lb_grid(2),&
                                                   lb_grid(3))
    ELSE IF (ASSOCIATED(rho%pw%cr3d)) THEN
      total_rho = -rho%pw%pw_grid%dvol*SUM(rho%pw%cr3d)
    ELSE IF (ASSOCIATED(rho%pw%cc)) THEN
      total_rho = -rho%pw%pw_grid%vol*rho%pw%cc(1)
    ELSE
      CALL stop_program("calculate_total_rho","not defined")
    END IF

  END FUNCTION calculate_total_rho

! *****************************************************************************

  SUBROUTINE collocate_pgf_product_gspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,scale,pab,pw)

    USE mathconstants, ONLY: pi,twopi

    TYPE(pw_type), TARGET, INTENT(INOUT) :: pw
    REAL(wp), INTENT(IN)                 :: rab2,scale,zeta,zetb
    INTEGER, INTENT(IN)                  :: la_max,la_min,lb_max,lb_min
    REAL(wp), DIMENSION(3), INTENT(IN)   :: ra,rab

    REAL(wp), DIMENSION(ncoset(la_max),ncoset(lb_max)), INTENT(IN) :: pab

!   *** Local variables ***

    COMPLEX(wp) :: yz,z
    REAL(wp)    :: f,fa,fb,pij,prefactor,rzetp,twozetp,zetp
    INTEGER     :: ax,ay,az,bx,by,bz,handle,ico,i,ig,ig2,jco,jg,kg,la,lb,&
                   lb_cube_min,ub_cube_max

    COMPLEX(wp), DIMENSION(3) :: phasefactor
    REAL(wp), DIMENSION(3)    :: dg,expfactor,fap,fbp,rap,rbp,rp
    INTEGER, DIMENSION(3)     :: lb_cube,ub_cube

    COMPLEX(wp), DIMENSION(:), POINTER :: rag,rbg
    REAL(wp), DIMENSION(:), POINTER    :: g

    COMPLEX(wp), DIMENSION(:,:,:,:), POINTER :: cubeaxis

!   ---------------------------------------------------------------------------

    CALL timeset("collocate_pgf_product_gspace","I","",handle)

    npgf_product = npgf_product + (ncoset(la_max) - ncoset(la_min-1))*&
                                  (ncoset(lb_max) - ncoset(lb_min-1))

    dg(:) = twopi/(pw%pw_grid%npts(:)*pw%pw_grid%dr(:))

    zetp = zeta + zetb
    rzetp = 1.0_wp/zetp
    f = zetb*rzetp
    rap(:) = f*rab(:)
    rbp(:) = rap(:) - rab(:)
    rp(:) = ra(:) + rap(:)
    twozetp = 2.0_wp*zetp
    fap(:) = twozetp*rap(:)
    fbp(:) = twozetp*rbp(:)

    prefactor = scale*SQRT((pi*rzetp)**3)*EXP(-zeta*f*rab2)
    phasefactor(:) = EXP(CMPLX(0.0_wp,-rp(:)*dg(:),wp))
    expfactor(:) = EXP(-0.25*rzetp*dg(:)*dg(:))

    lb_cube(:) = pw%pw_grid%bounds(1,:)
    ub_cube(:) = pw%pw_grid%bounds(2,:)

    lb_cube_min = MINVAL(lb_cube(:))
    ub_cube_max = MAXVAL(ub_cube(:))

    NULLIFY (cubeaxis,g,rag,rbg)

    cubeaxis => reallocate(cubeaxis,lb_cube_min,ub_cube_max,1,3,0,la_max,&
                           0,lb_max)
    g => reallocate(g,lb_cube_min,ub_cube_max)
    rag => reallocate(rag,lb_cube_min,ub_cube_max)
    rbg => reallocate(rbg,lb_cube_min,ub_cube_max)

    DO i=1,3

      DO ig=lb_cube(i),ub_cube(i)
        ig2 = ig*ig
        cubeaxis(ig,i,0,0) = expfactor(i)**ig2*phasefactor(i)**ig
      END DO

      IF (la_max > 0) THEN
        DO ig=lb_cube(i),ub_cube(i)
          g(ig) = REAL(ig,wp)*dg(i)
          rag(ig) = CMPLX(fap(i),-g(ig),wp)
          cubeaxis(ig,i,1,0) = rag(ig)*cubeaxis(ig,i,0,0)
        END DO
        DO la=2,la_max
          fa = REAL(la-1,wp)*twozetp
          DO ig=lb_cube(i),ub_cube(i)
            cubeaxis(ig,i,la,0) = rag(ig)*cubeaxis(ig,i,la-1,0) +&
                                  fa*cubeaxis(ig,i,la-2,0)
          END DO
        END DO
        IF (lb_max > 0) THEN
          fa = twozetp
          DO ig=lb_cube(i),ub_cube(i)
            rbg(ig) = CMPLX(fbp(i),-g(ig),wp)
            cubeaxis(ig,i,0,1) = rbg(ig)*cubeaxis(ig,i,0,0)
            cubeaxis(ig,i,1,1) = rbg(ig)*cubeaxis(ig,i,1,0) +&
                                 fa*cubeaxis(ig,i,0,0)
          END DO
          DO lb=2,lb_max
            fb = REAL(lb-1,wp)*twozetp
            DO ig=lb_cube(i),ub_cube(i)
              cubeaxis(ig,i,0,lb) = rbg(ig)*cubeaxis(ig,i,0,lb-1) +&
                                    fb*cubeaxis(ig,i,0,lb-2)
              cubeaxis(ig,i,1,lb) = rbg(ig)*cubeaxis(ig,i,1,lb-1) +&
                                    fb*cubeaxis(ig,i,1,lb-2) +&
                                    fa*cubeaxis(ig,i,0,lb-1)
            END DO
          END DO
          DO la=2,la_max
            fa = REAL(la,wp)*twozetp
            DO lb=1,lb_max
              fb = REAL(lb-1,wp)*twozetp
              DO ig=lb_cube(i),ub_cube(i)
                cubeaxis(ig,i,la,lb) = rbg(ig)*cubeaxis(ig,i,la,lb-1) +&
                                       fb*cubeaxis(ig,i,la,lb-2) +&
                                       fa*cubeaxis(ig,i,la-1,lb-1)
              END DO
            END DO
          END DO
        END IF
      ELSE
        IF (lb_max > 0) THEN
          DO ig=lb_cube(i),ub_cube(i)
            g(ig) = REAL(ig,wp)*dg(i)
            rbg(ig) = CMPLX(fbp(i),-g(ig),wp)
            cubeaxis(ig,i,0,1) = rbg(ig)*cubeaxis(ig,i,0,0)
          END DO
          DO lb=2,lb_max
            fb = REAL(lb-1,wp)*twozetp
            DO ig=lb_cube(i),ub_cube(i)
              cubeaxis(ig,i,0,lb) = rbg(ig)*cubeaxis(ig,i,0,lb-1) +&
                                    fb*cubeaxis(ig,i,0,lb-2)
            END DO
          END DO
        END IF
      END IF

    END DO

    DO la=0,la_max
      DO lb=0,lb_max
        IF (la + lb == 0) CYCLE
        fa = (1.0_wp/twozetp)**(la + lb)
        DO i=1,3
          DO ig=lb_cube(i),ub_cube(i)
            cubeaxis(ig,i,la,lb) = fa*cubeaxis(ig,i,la,lb)
          END DO
        END DO
      END DO
    END DO

!   *** Add the current primitive Gaussian function product to grid ***

    DO ico=ncoset(la_min-1)+1,ncoset(la_max)

      ax = indco(1,ico)
      ay = indco(2,ico)
      az = indco(3,ico)

      DO jco=ncoset(lb_min-1)+1,ncoset(lb_max)

        pij = prefactor*pab(ico,jco)

        IF (ABS(pij) < eps_rho_gspace) CYCLE

        bx = indco(1,jco)
        by = indco(2,jco)
        bz = indco(3,jco)

        DO kg=lb_cube(3),ub_cube(3)
          z = pij*cubeaxis(kg,3,az,bz)
          IF (ABS(z) < eps_rho_gspace) CYCLE
          DO jg=lb_cube(2),ub_cube(2)
            yz = cubeaxis(jg,2,ay,by)*z
            IF (ABS(yz) < eps_rho_gspace) CYCLE
            DO ig=lb_cube(1),ub_cube(1)
              pw%cc3d(ig,jg,kg) = pw%cc3d(ig,jg,kg) + cubeaxis(ig,1,ax,bx)*yz
            END DO
          END DO
        END DO

      END DO

    END DO

    DEALLOCATE (cubeaxis,g,rag,rbg)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE collocate_pgf_product_gspace

! *****************************************************************************

  SUBROUTINE collocate_pgf_product_rspace(la_max,zeta,la_min,&
                                          lb_max,zetb,lb_min,&
                                          ra,rab,rab2,scale,pab,pw,cube_info,l_info)

    TYPE(pw_type), TARGET, INTENT(INOUT) :: pw
    TYPE(cube_info_type), INTENT(IN)     :: cube_info
    TYPE(l_info_type), INTENT(IN)        :: l_info
    REAL(wp), INTENT(IN)                 :: rab2,scale,zeta,zetb
    INTEGER, INTENT(IN)                  :: la_max,la_min,lb_max,lb_min
    REAL(wp), DIMENSION(3), INTENT(IN)   :: ra,rab

    REAL(wp), DIMENSION(ncoset(la_max),ncoset(lb_max)), INTENT(IN) :: pab

!   *** Local variables ***

    REAL(wp) :: f,prefactor,radius,zetp
    INTEGER  :: i,ico,ig,j,jco,jg,k,kg

    REAL(wp), DIMENSION(3) :: dr,rap,rbp,roffset,rp,rb
    INTEGER, DIMENSION(3)  :: cubecenter,lb_cube,ng,ub_cube,ub_grid

    INTEGER, DIMENSION(:,:), POINTER      :: map
    REAL(wp), DIMENSION(:,:,:), POINTER   :: grid
    INTEGER lxyz_max,lxy_max,lx_max,lx,lxb,lya,lyb,lza,lzb,lz,ly,icoef,lxy
    INTEGER coef_max,lxa,lxyz,l,gridbounds(2,3)
    INTEGER, pointer,dimension(:) :: ly_max,lz_max
    REAL(wp) a,b,binomial_k_lxa,binomial_l_lxb,zbp,zap,za,zb,ya,yb,yap,ybp,rpg,s,pg

    REAL(wp), pointer :: dpx(:,:)
    REAL(wp), pointer :: pzyx(:), alpha(:,:)
    REAL(wp), pointer :: polz(:,:),dpz(:,:),poly(:,:),dpy(:,:),polx(:,:)
    INTEGER, pointer, dimension(:)  :: sphere_bounds
    REAL(wp)  :: xa,xb,xap,xbp,polypart,cutoff
    INTEGER  :: start,length,offset,cmax

!   ---------------------------------------------------------------------------

    coef_max=la_max+lb_max+1

    zetp      = zeta + zetb
    f         = zetb/zetp
    prefactor = scale*EXP(-zeta*f*rab2)
    rap(:)    = f*rab(:)
    rbp(:)    = rap(:) - rab(:)
    rp(:)     = ra(:) + rap(:)  
    rb(:)     = ra(:)+rab(:)
    cutoff    = 0.0_wp

    radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,ra,rb,rp, &
                   zetp,eps_rho_rspace,prefactor,cutoff)

    IF (radius .eq. 0.0_wp ) THEN
      RETURN
    END IF

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

!   *** irrelevant ? ***
    npgf_product = npgf_product + (ncoset(la_max) - ncoset(la_min-1))*&
                                  (ncoset(lb_max) - ncoset(lb_min-1))

!   *** get the sub grid properties for the given radius ***
    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)

!   *** get the l_info logic and arrays ***
    CALL return_l_info(l_info,la_min,la_max,lb_min,lb_max,lx_max, &
                       lxy_max,lxyz_max,ly_max,lz_max, &
                       map,polx,poly,polz,dpy,dpz,alpha,pzyx,cmax)

!   *** position of the gaussian product

    cubecenter(:) = FLOOR(rp(:)/dr(:))
    roffset(:)    = rp(:) - REAL(cubecenter(:),wp)*dr(:)

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

!   *** initialise the p terms and loop logic
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
             pzyx(lxyz)=prefactor*pab(ico,jco)
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO

 
!   *** initialise the pol x,y,z terms
    DO ig=lb_cube(3),ub_cube(3)
      lxyz=0
      rpg = REAL(ig,wp)*dr(3) - roffset(3)
      zap = EXP(-zetp*rpg**2)
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
        a=a*(-ra(1)+rp(1))
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

    CALL collocate_core(pzyx,polx,poly,polz,grid,alpha,lx_max,lxy_max,&
                lxyz_max,coef_max, &
                cmax,ly_max,lz_max,gridbounds,&
                map,sphere_bounds)

  END SUBROUTINE collocate_pgf_product_rspace

! ***********************************************************************************************
  SUBROUTINE collocate_core(pzyx,polx,   &
                poly,polz,grid,alpha,lx_max,lxy_max,lxyz_max,coef_max, &
                cmax,ly_max,lz_max,gridbounds,&
                map,sphere_bounds)
  IMPLICIT NONE
    integer, INTENT(IN) ::  lx_max,lxy_max,lxyz_max,coef_max,cmax
    integer, INTENT(IN) ::  ly_max(lx_max),lz_max(lxy_max)
    integer, INTENT(IN) ::  map(-cmax:cmax,1:3)
    integer, INTENT(IN) ::  gridbounds(2,3)
    integer, INTENT(IN) :: sphere_bounds(:)
    real(wp), INTENT(IN) :: pzyx(lxyz_max)
    real(wp), INTENT(IN) :: polz(lxyz_max,-cmax:cmax)
    real(wp), INTENT(IN) :: poly(lxy_max,-cmax:cmax)
    real(wp), INTENT(IN) :: polx(coef_max,-cmax:cmax)
    real(wp), INTENT(IN) :: alpha(coef_max,lx_max)
    real(wp) ::  coef(coef_max,4)
    real(wp) ::  pyx(2,lxy_max)
! help the compiler not to copy the array, but allow optimization on other machines
#ifdef __PGI
    real(wp), DIMENSION(:,:,:), POINTER :: grid
#else
    real(wp), INTENT(INOUT) :: grid(gridbounds(1,1):gridbounds(2,1),gridbounds(1,2):gridbounds(2,2),gridbounds(1,3):gridbounds(2,3))
#endif

    integer kg,kgmax,k,jgmax,jgmin,jg,j,j2,igmax,ig,i,kgmin,igmin,k2,i2,jg2,kg2
    integer lx,ly,lz,lxy,lxyz,icoef,sci
    real(wp) s1,s2,s3,s4,s5,s6,s7,s8,s9,sa,g1,g2,g3
    real(wp) s1k,s2k,s3k,s4k,s5k,s6k,s7k,s8k,s9k,sak,g1k,g2k
    real(wp) r1,r2,r3,r4,r5,r6,r7,r8

!    call f_hpmstart(11,"collo")
    sci=1
    kgmin=sphere_bounds(sci)
    sci=sci+1
    DO kg=kgmin,0
       kg2=1-kg
       k=map(kg,3)
       k2=map(kg2,3)
       lxyz=0
       DO lxy=1,lxy_max
          pyx(1,lxy)=0.0_wp
          pyx(2,lxy)=0.0_wp
          DO lz=1,lz_max(lxy)
             lxyz=lxyz+1
             pyx(1,lxy)=pyx(1,lxy)+pzyx(lxyz)*polz(lxyz,kg)
             pyx(2,lxy)=pyx(2,lxy)+pzyx(lxyz)*polz(lxyz,kg2)
          ENDDO
       ENDDO
       jgmin=sphere_bounds(sci)
       sci=sci+1
       DO jg=jgmin,0
          jg2=1-jg
          j=map(jg,2)
          j2=map(jg2,2)
          igmin=sphere_bounds(sci)
          sci=sci+1
          igmax=1-igmin
          SELECT CASE (coef_max)
          CASE(1)
             s1=0.0_wp
             s2=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             lxy=0
             DO ly=1,ly_max(1)
                lxy=lxy+1
                s1=s1+pyx(1,lxy)*poly(lxy,jg)
                s2=s2+pyx(1,lxy)*poly(lxy,jg2)
                s1k=s1k+pyx(2,lxy)*poly(lxy,jg)
                s2k=s2k+pyx(2,lxy)*poly(lxy,jg2)
             ENDDO
             s1=alpha(1,1)*s1
             s2=alpha(1,1)*s2
             s1k=alpha(1,1)*s1k
             s2k=alpha(1,1)*s2k
             DO ig=igmin,igmax,2
                i=map(ig,1)
                i2=map(ig+1,1)
                g1=polx(1,ig)
                g2=polx(1,ig+1)
                r1=grid(i,j,k)
                r2=grid(i,j2,k)
                r3=grid(i,j,k2)
                r4=grid(i,j2,k2)
                r1= r1 + s1*g1
                r2= r2 + s2*g1
                r3= r3 + s1k*g1
                r4= r4 + s2k*g1
                r5=grid(i2,j,k)
                r6=grid(i2,j2,k)
                r7=grid(i2,j,k2)
                r8=grid(i2,j2,k2)
                grid(i,j,k)   = r1 
                grid(i,j2,k)  = r2
                grid(i,j,k2)  = r3
                grid(i,j2,k2) = r4 
                r5 = r5 + s1*g2
                r6 = r6 + s2*g2
                r7 = r7 + s1k*g2
                r8 = r8 + s2k*g2
                grid(i2,j,k)   = r5 
                grid(i2,j2,k)  = r6 
                grid(i2,j,k2)  = r7
                grid(i2,j2,k2) = r8 
             END DO
          CASE(2)
             lxy=0
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             DO lx=1,2
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=pyx(1,lxy)
                 r2=pyx(2,lxy)
                 r3=poly(lxy,jg)
                 r4=poly(lxy,jg2)
                 g1=g1+r1*r3
                 g2=g2+r1*r4
                 g1k=g1k+r2*r3
                 g2k=g2k+r2*r4
              ENDDO
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              s1=s1+r1*g1
              s2=s2+r2*g1
              s3=s3+r1*g2
              s4=s4+r2*g2
              s1k=s1k+r1*g1k
              s2k=s2k+r2*g1k
              s3k=s3k+r1*g2k
              s4k=s4k+r2*g2k
             ENDDO
             DO ig=igmin,igmax,2
                i=map(ig,1)
                i2=map(ig+1,1)
                g1=polx(1,ig)
                g2=polx(1,ig+1)
                r1=grid(i,j,k)+ s1*g1
                r2=grid(i,j2,k)+ s3*g1
                r3=grid(i,j,k2)+ s1k*g1
                r4=grid(i,j2,k2)+ s3k*g1
                r5=grid(i2,j,k)+ s1*g2
                r6=grid(i2,j2,k)+ s3*g2
                r7=grid(i2,j,k2)+ s1k*g2
                r8=grid(i2,j2,k2)+ s3k*g2
                g1=polx(2,ig)
                g2=polx(2,ig+1)
                grid(i,j,k) = r1+ s2*g1
                grid(i,j2,k) = r2+ s4*g1
                grid(i,j,k2) = r3+ s2k*g1
                grid(i,j2,k2) = r4+ s4k*g1
                grid(i2,j,k) = r5+ s2*g2
                grid(i2,j2,k) = r6+ s4*g2
                grid(i2,j,k2) = r7+ s2k*g2
                grid(i2,j2,k2) = r8+ s4k*g2
             END DO
          CASE(3)
             lxy=0
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
             DO lx=1,lx_max
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=pyx(1,lxy)
                 r2=pyx(2,lxy)
                 r3=poly(lxy,jg)
                 r4=poly(lxy,jg2)
                 g1=g1+r1*r3
                 g2=g2+r1*r4
                 g1k=g1k+r2*r3
                 g2k=g2k+r2*r4
              ENDDO
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              r3=alpha(3,lx)
              s1=s1+r1*g1
              s2=s2+r2*g1
              s3=s3+r3*g1
              s4=s4+r1*g2
              s5=s5+r2*g2
              s6=s6+r3*g2
              s1k=s1k+r1*g1k
              s2k=s2k+r2*g1k
              s3k=s3k+r3*g1k
              s4k=s4k+r1*g2k
              s5k=s5k+r2*g2k
              s6k=s6k+r3*g2k
             ENDDO
             DO ig=igmin,igmax,2
                i=map(ig,1)
                i2=map(ig+1,1)
                g1=polx(1,ig)
                g2=polx(1,ig+1)
                r1=grid(i,j,k)+ s1*g1
                r2=grid(i,j2,k)+ s4*g1
                r3=grid(i,j,k2)+ s1k*g1
                r4=grid(i,j2,k2)+ s4k*g1
                r5=grid(i2,j,k)+ s1*g2
                r6=grid(i2,j2,k)+ s4*g2
                r7=grid(i2,j,k2)+ s1k*g2
                r8=grid(i2,j2,k2)+ s4k*g2
                g1=polx(2,ig)
                g2=polx(2,ig+1)
                r1=r1+s2*g1
                r2=r2+s5*g1
                r3=r3+s2k*g1
                r4=r4+s5k*g1
                r5=r5+s2*g2
                r6=r6+s5*g2
                r7=r7+s2k*g2
                r8=r8+s5k*g2
                g1=polx(3,ig)
                g2=polx(3,ig+1)
                grid(i,j,k) = r1 + s3*g1
                grid(i,j2,k)= r2 + s6*g1
                grid(i,j,k2) = r3 + s3k*g1
                grid(i,j2,k2) = r4 + s6k*g1
                grid(i2,j,k) = r5 + s3*g2
                grid(i2,j2,k)= r6 + s6*g2
                grid(i2,j,k2) = r7 + s3k*g2
                grid(i2,j2,k2) = r8 + s6k*g2
             END DO
          CASE(4)
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
             lxy=0
             DO lx=1,lx_max
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 r1=pyx(1,lxy)
                 r2=pyx(2,lxy)
                 r3=poly(lxy,jg)
                 r4=poly(lxy,jg2)
                 g1=g1+r1*r3
                 g2=g2+r1*r4
                 g1k=g1k+r2*r3
                 g2k=g2k+r2*r4
              ENDDO
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              r3=alpha(3,lx)
              r4=alpha(4,lx)
              s1=s1+r1*g1
              s2=s2+r2*g1
              s3=s3+r3*g1
              s4=s4+r4*g1
              s5=s5+r1*g2
              s6=s6+r2*g2
              s7=s7+r3*g2
              s8=s8+r4*g2
              s1k=s1k+r1*g1k
              s2k=s2k+r2*g1k
              s3k=s3k+r3*g1k
              s4k=s4k+r4*g1k
              s5k=s5k+r1*g2k
              s6k=s6k+r2*g2k
              s7k=s7k+r3*g2k
              s8k=s8k+r4*g2k
             ENDDO
             DO ig=igmin,igmax,2
                i=map(ig,1)
                i2=map(ig+1,1)
                g1=polx(1,ig)
                g2=polx(1,ig+1)
                r1=grid(i,j,k)+ s1*g1
                r2=grid(i,j2,k)+ s5*g1
                r3=grid(i,j,k2)+ s1k*g1
                r4=grid(i,j2,k2)+ s5k*g1
                r5=grid(i2,j,k)+ s1*g2
                r6=grid(i2,j2,k)+ s5*g2
                r7=grid(i2,j,k2)+ s1k*g2
                r8=grid(i2,j2,k2)+ s5k*g2
                g1=polx(2,ig)
                g2=polx(2,ig+1)
                r1=r1+s2*g1
                r2=r2+s6*g1
                r3=r3+s2k*g1
                r4=r4+s6k*g1
                r5=r5+s2*g2
                r6=r6+s6*g2
                r7=r7+s2k*g2
                r8=r8+s6k*g2
                g1=polx(3,ig)
                g2=polx(3,ig+1)
                r1=r1+s3*g1
                r2=r2+s7*g1
                r3=r3+s3k*g1
                r4=r4+s7k*g1
                r5=r5+s3*g2
                r6=r6+s7*g2
                r7=r7+s3k*g2
                r8=r8+s7k*g2
                g1=polx(4,ig)
                g2=polx(4,ig+1)
                grid(i,j,k) = r1 + s4*g1
                grid(i,j2,k) = r2+ s8*g1
                grid(i,j,k2) = r3+ s4k*g1
                grid(i,j2,k2) = r4+ s8k*g1
                grid(i2,j,k) = r5 + s4*g2
                grid(i2,j2,k) = r6+ s8*g2
                grid(i2,j,k2) = r7+ s4k*g2
                grid(i2,j2,k2) = r8+ s8k*g2
             END DO
          CASE(5)
             lxy=0
             s1=0.0_wp
             s2=0.0_wp
             s3=0.0_wp
             s4=0.0_wp
             s9=0.0_wp
             s5=0.0_wp
             s6=0.0_wp
             s7=0.0_wp
             s8=0.0_wp
             sa=0.0_wp
             s1k=0.0_wp
             s2k=0.0_wp
             s3k=0.0_wp
             s4k=0.0_wp
             s9k=0.0_wp
             s5k=0.0_wp
             s6k=0.0_wp
             s7k=0.0_wp
             s8k=0.0_wp
             sak=0.0_wp
             DO lx=1,lx_max
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO ly=1,ly_max(lx)
                 lxy=lxy+1 
                 r1=pyx(1,lxy)
                 r2=pyx(2,lxy)
                 r3=poly(lxy,jg)
                 r4=poly(lxy,jg2)
                 g1=g1+r1*r3
                 g2=g2+r1*r4
                 g1k=g1k+r2*r3
                 g2k=g2k+r2*r4
              ENDDO
              r1=alpha(1,lx)
              r2=alpha(2,lx)
              r3=alpha(3,lx)
              r4=alpha(4,lx)
              r5=alpha(5,lx)
              s1=s1+r1*g1
              s2=s2+r2*g1
              s3=s3+r3*g1
              s4=s4+r4*g1
              s9=s9+r5*g1
              s5=s5+r1*g2
              s6=s6+r2*g2
              s7=s7+r3*g2
              s8=s8+r4*g2
              sa=sa+r5*g2
              s1k=s1k+r1*g1k
              s2k=s2k+r2*g1k
              s3k=s3k+r3*g1k
              s4k=s4k+r4*g1k
              s9k=s9k+r5*g1k
              s5k=s5k+r1*g2k
              s6k=s6k+r2*g2k
              s7k=s7k+r3*g2k
              s8k=s8k+r4*g2k
              sak=sak+r5*g2k
             ENDDO
             DO ig=igmin,igmax,2
                i=map(ig,1)
                i2=map(ig+1,1)
                g1=polx(1,ig)
                g2=polx(1,ig+1)
                r1= grid(i,j,k) + s1*g1
                r2=grid(i,j2,k) + s5*g1
                r3=grid(i,j,k2) + s1k*g1
                r4=grid(i,j2,k2)+ s5k*g1
                r5= grid(i2,j,k) + s1*g2
                r6=grid(i2,j2,k) + s5*g2
                r7=grid(i2,j,k2) + s1k*g2
                r8=grid(i2,j2,k2)+ s5k*g2
                g1=polx(2,ig)
                g2=polx(2,ig+1)
                r1=r1+s2*g1
                r2=r2+s6*g1
                r3=r3+s2k*g1
                r4=r4+s6k*g1
                r5=r5+s2*g2
                r6=r6+s6*g2
                r7=r7+s2k*g2
                r8=r8+s6k*g2
                g1=polx(3,ig)
                g2=polx(3,ig+1)
                r1=r1+s3*g1
                r2=r2+s7*g1
                r3=r3+s3k*g1
                r4=r4+s7k*g1
                r5=r5+s3*g2
                r6=r6+s7*g2
                r7=r7+s3k*g2
                r8=r8+s7k*g2
                g1=polx(4,ig)
                g2=polx(4,ig+1)
                r1=r1+s4*g1
                r2=r2+s8*g1
                r3=r3+s4k*g1
                r4=r4+s8k*g1
                r5=r5+s4*g2
                r6=r6+s8*g2
                r7=r7+s4k*g2
                r8=r8+s8k*g2
                g1=polx(5,ig)
                g2=polx(5,ig+1)
                grid(i,j,k) = r1 + s9*g1
                grid(i,j2,k)= r2 + sa*g1
                grid(i,j,k2)= r3 + s9k*g1
                grid(i,j2,k2)= r4+ sak*g1
                grid(i2,j,k) = r5 + s9*g2
                grid(i2,j2,k)= r6 + sa*g2
                grid(i2,j,k2)= r7 + s9k*g2
                grid(i2,j2,k2)= r8+ sak*g2
             END DO

          CASE DEFAULT
             coef(:,:)=0.0
             lxy=0
             DO lx=1,lx_max
              g1=0.0_wp
              g2=0.0_wp
              g1k=0.0_wp
              g2k=0.0_wp
              DO ly=1,ly_max(lx)
                 lxy=lxy+1
                 g1=g1+pyx(1,lxy)*poly(lxy,jg)
                 g2=g2+pyx(1,lxy)*poly(lxy,jg2)
                 g1k=g1k+pyx(2,lxy)*poly(lxy,jg)
                 g2k=g2k+pyx(2,lxy)*poly(lxy,jg2)
              ENDDO
              DO icoef=1,coef_max
                 coef(icoef,1)=coef(icoef,1)+alpha(icoef,lx)*g1
                 coef(icoef,2)=coef(icoef,2)+alpha(icoef,lx)*g2
                 coef(icoef,3)=coef(icoef,3)+alpha(icoef,lx)*g1k
                 coef(icoef,4)=coef(icoef,4)+alpha(icoef,lx)*g2k
              ENDDO
             ENDDO
             DO ig=igmin,igmax
                i=map(ig,1)
                s1=0.0_wp
                s2=0.0_wp
                s1k=0.0_wp
                s2k=0.0_wp
                DO icoef=1,coef_max
                  s1=s1+coef(icoef,1)*polx(icoef,ig)
                  s2=s2+coef(icoef,2)*polx(icoef,ig)
                  s1k=s1k+coef(icoef,3)*polx(icoef,ig)
                  s2k=s2k+coef(icoef,4)*polx(icoef,ig)
                ENDDO
                grid(i,j,k) = grid(i,j,k)     + s1
                grid(i,j2,k) = grid(i,j2,k)   + s2
                grid(i,j,k2) = grid(i,j,k2)   + s1k
                grid(i,j2,k2) = grid(i,j2,k2) + s2k
             END DO
          END SELECT
       END DO
    END DO

!    call f_hpmstop(11)
  END SUBROUTINE

! *****************************************************************************

  SUBROUTINE ineighbor_loop_work ( &
       ineighbor, iatom, jatom, &
       p, p_old, rho, rebuild_rho, &
       basis_set_a, nseta, radius_kind_a, &
       ra, rab, rab2, &
       gridlevel_info, cube_info, l_info, &
       pab, pab_old, work )

    IMPLICIT NONE

    ! ARGUMENTS
    INTEGER, INTENT(IN) :: ineighbor, iatom, jatom, nseta

    LOGICAL, INTENT(IN) :: rebuild_rho

    REAL(wp), INTENT(IN) :: radius_kind_a, rab2
    REAL(wp), DIMENSION(3), INTENT(IN) :: ra, rab
    REAL(wp), DIMENSION(:,:), POINTER :: pab, pab_old, work ! INTENT(INOUT)

    TYPE(coeff_type), POINTER, DIMENSION(:) :: rho ! INTENT(INOUT)
    TYPE(gridlevel_info_type), INTENT(IN) :: gridlevel_info
    TYPE(gto_basis_set_type), POINTER :: basis_set_a ! INTENT(IN)
    TYPE(real_matrix_set_type), INTENT(IN) :: p, p_old
    TYPE(cube_info_type), POINTER, DIMENSION(:) :: cube_info
    TYPE(l_info_type), INTENT(IN) :: l_info

    ! LOCALS
    INTEGER :: jpgf, jset, &
         first_jset, first_sgfa, first_jpgf, &
         first_sgfb, ipgf, iset, &
         la_max, la_min, last_sgfa, last_sgfb, lb_max, lb_min, &
         na1, na2, nb1, nb2, nsgfa, nsgfb, ncoa, ncob, npgfa, npgfb, &
         nsetb, nshella, nshellb, &
         igrid_level

    REAL(wp) :: scale, dab, radius_set_a, radius_kind_b, radius_set_b
    REAL(wp), DIMENSION(:), POINTER :: rpgfa, rpgfb, zeta, zetb
    REAL(wp), DIMENSION(:,:), POINTER :: p_block, p_old_block

    TYPE(gto_basis_set_type), POINTER :: basis_set_b

!   ---------------------------------------------------------------------------

    dab = SQRT(rab2)
    basis_set_b => kind_info(atom_info(jatom)%kind)%orb_basis_set
    nsetb = basis_set_b%nset
    radius_kind_b = basis_set_b%kind_radius

    CALL get_block_node(matrix=p%matrix,&
         block_row=iatom,&
         block_col=jatom,&
         block=p_block)


!!!    IF (.NOT.ASSOCIATED(p_block)) CYCLE
    IF (.NOT.ASSOCIATED(p_block)) RETURN

    IF (rebuild_rho) then
!!!       IF (MAXVAL(ABS(p_block(:,:))) < 1.0E-18_wp) CYCLE
       IF (MAXVAL(ABS(p_block(:,:))) < 1.0E-18_wp) RETURN
    ELSE
       CALL get_block_node(matrix=p_old%matrix,&
            block_row=iatom,&
            block_col=jatom,&
            block=p_old_block)
!!!       IF (MAXVAL(ABS(p_block(:,:)-p_old_block(:,:))) < 1.0E-18_wp) CYCLE
       IF (MAXVAL(ABS(p_block(:,:)-p_old_block(:,:))) < 1.0E-18_wp) RETURN
    ENDIF


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

          ! *** Decontract the density matrix set block ***

          CALL dgemm("N","N",ncoa,nsgfb,nsgfa,1.0_wp,&
               basis_set_a%sphi(1,first_sgfa),&
               SIZE(basis_set_a%sphi,1),&
               p_block(first_sgfa,first_sgfb),SIZE(p_block,1),&
               0.0_wp,work(1,1),SIZE(work,1))
          CALL dgemm("N","T",ncoa,ncob,nsgfb,1.0_wp,&
               work(1,1),SIZE(work,1),basis_set_b%sphi(1,first_sgfb),&
               SIZE(basis_set_b%sphi,1),0.0_wp,pab(1,1),SIZE(pab,1))

          if (.not.rebuild_rho) then
             CALL dgemm("N","N",ncoa,nsgfb,nsgfa,1.0_wp,&
                  basis_set_a%sphi(1,first_sgfa),&
                  SIZE(basis_set_a%sphi,1),&
                  p_old_block(first_sgfa,first_sgfb),SIZE(p_old_block,1),&
                  0.0_wp,work(1,1),SIZE(work,1))
             CALL dgemm("N","T",ncoa,ncob,nsgfb,1.0_wp,&
                  work(1,1),SIZE(work,1),basis_set_b%sphi(1,first_sgfb),&
                  SIZE(basis_set_b%sphi,1),0.0_wp,pab_old(1,1),SIZE(pab_old,1))

             pab(:,:)=pab(:,:)-pab_old(:,:)
          endif


          DO ipgf=1,npgfa

             IF (rpgfa(ipgf) + radius_set_b < dab) CYCLE

             na1 = (ipgf - 1)*ncoset(la_max) + 1
             na2 = ipgf*ncoset(la_max)

             IF ((iatom == jatom).AND.(iset == jset)) THEN
                first_jpgf = ipgf
             ELSE
                first_jpgf = 1
             END IF

             DO jpgf=first_jpgf,npgfb

                IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) CYCLE

                nb1 = (jpgf - 1)*ncoset(lb_max) + 1
                nb2 = jpgf*ncoset(lb_max)

                IF ((iatom == jatom).AND.&
                     (iset == jset).AND.&
                     (ipgf == jpgf)) THEN
                   scale = 1.0_wp
                ELSE
                   scale = 2.0_wp
                END IF
                ! at this point calling this with gspace might be rather wrong
                IF (ASSOCIATED(rho(1)%pw%cc3d)) THEN
                   CALL collocate_pgf_product_gspace(la_max,zeta(ipgf),la_min,&
                        lb_max,zetb(jpgf),lb_min,&
                        ra,rab,rab2,scale,&
                        pab(na1:na2,nb1:nb2),&
                        rho(1)%pw)
                ELSE IF (ASSOCIATED(rho(1)%pw%cr3d)) THEN
                   igrid_level=gaussian_gridlevel(gridlevel_info,zeta(ipgf)+zetb(jpgf))
                   CALL collocate_pgf_product_rspace(la_max,zeta(ipgf),la_min,&
                        lb_max,zetb(jpgf),lb_min,&
                        ra,rab,rab2,scale,&
                        pab(na1:na2,nb1:nb2),&
                        rho(igrid_level)%pw, &
                        cube_info(igrid_level),l_info)
                END IF

             END DO ! jpgf

          END DO ! ipgf

       END DO ! jset

    END DO ! iset
          
  END SUBROUTINE ineighbor_loop_work

! *****************************************************************************

END MODULE collocate_density
