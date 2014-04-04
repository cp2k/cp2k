!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE qs_integrate_potential_low
  USE cell_types,                      ONLY: cell_type
  USE cube_utils,                      ONLY: compute_cube_center,&
                                             cube_info_type,&
                                             return_cube,&
                                             return_cube_nonortho
  USE d3_poly,                         ONLY: poly_d32cp2k
  USE gauss_colloc,                    ONLY: integrateGaussFull
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants,                   ONLY: fac
  USE orbital_pointers,                ONLY: coset,&
                                             current_maxl,&
                                             ncoset
  USE qs_interactions,                 ONLY: exp_radius_very_extended
  USE realspace_grid_types,            ONLY: realspace_grid_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integrate_potential_low'

! *** Don't include this routines directly, use the interface to
! *** qs_integrate_potential

  PUBLIC :: integrate_pgf_product_rspace

CONTAINS

! *****************************************************************************
!> \brief low level function to compute matrix elements of primitive gaussian functions
! *****************************************************************************
    SUBROUTINE integrate_pgf_product_rspace(la_max,zeta,la_min,&
                                            lb_max,zetb,lb_min,&
                                            ra,rab,rab2,rsgrid,cell,&
                                            cube_info,hab,pab,o1,o2,&
                                            eps_gvg_rspace,&
                                            calculate_forces,hdab,hadb,force_a,force_b,&
                                            compute_tau,map_consistent,&
                                            collocate_rho0,rpgf0_s,use_virial,my_virial_a,&
                                            my_virial_b,a_hdab,use_subpatch,subpatch_pattern,error)

    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), INTENT(IN)                :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), INTENT(IN)                :: zetb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rab
    REAL(KIND=dp), INTENT(IN)                :: rab2
    TYPE(realspace_grid_type), POINTER       :: rsgrid
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: hab
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: pab
    INTEGER, INTENT(IN)                      :: o1, o2
    REAL(KIND=dp), INTENT(IN)                :: eps_gvg_rspace
    LOGICAL, INTENT(IN)                      :: calculate_forces
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: hdab, hadb
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT), OPTIONAL                :: force_a, force_b
    LOGICAL, INTENT(IN), OPTIONAL            :: compute_tau, map_consistent, &
                                                collocate_rho0
    REAL(dp), INTENT(IN), OPTIONAL           :: rpgf0_s
    LOGICAL, INTENT(IN), OPTIONAL            :: use_virial
    REAL(KIND=dp), DIMENSION(3,3), OPTIONAL  :: my_virial_a, my_virial_b
    REAL(KIND=dp), DIMENSION(:,:,:,:), OPTIONAL, POINTER  :: a_hdab
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, OPTIONAL                        :: use_subpatch
    INTEGER(KIND=int_8), INTENT(IN), OPTIONAL :: subpatch_pattern

    CHARACTER(len=*), PARAMETER :: routineN = 'integrate_pgf_product_rspace', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, by, bz, cmax, coef_max, gridbounds(2,3), i, &
      ico, icoef, ig, j, jco, k, l, la, la_max_local, la_min_local, lb, &
      lb_cube_min, lb_max_local, lb_min_local, length, lx, lx_max, lxa, lxb, lxy, &
      lxy_max, lxyz, lxyz_max, lya, lyb, lza, lzb, offset, start, ub_cube_max, handle
    INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ng, &
                                                ub_cube
    INTEGER, DIMENSION(:), POINTER           :: ly_max, lz_max, sphere_bounds
    LOGICAL                                  :: my_collocate_rho0, &
                                                my_compute_tau, &
                                                my_map_consistent, &
                                                my_use_virial,&
                                                subpatch_integrate
    REAL(KIND=dp) :: a, axpm0, b, binomial_k_lxa, binomial_l_lxb, cutoff, &
      der_a(3), der_b(3), exp_x0, exp_x1, exp_x2, f, ftza, ftzb, pabval, pg, &
      prefactor, radius, rpg, ya, yap, yb, ybp, za, zap, zb, zbp, zetp
    REAL(KIND=dp), DIMENSION(3)              :: dr, rap, rb, rbp, roffset, rp
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid

    INTEGER :: lxp,lyp,lzp,lp,iaxis
    INTEGER,       ALLOCATABLE, DIMENSION(:,:) :: map
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: alpha
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xyz
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xyt
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coef_xtt
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: coef_ttz
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: coef_tyz

    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: pol_z
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: pol_y
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: pol_x
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:) :: vab

    REAL(KIND=dp) :: t_exp_1,t_exp_2,t_exp_min_1,t_exp_min_2,t_exp_plus_1,t_exp_plus_2
    LOGICAL  :: failure


    failure = .FALSE.
    subpatch_integrate = .FALSE.

    IF(PRESENT(use_subpatch)) THEN
       IF(use_subpatch)THEN
         subpatch_integrate = .TRUE.
         CPPrecondition(PRESENT(subpatch_pattern),cp_failure_level,routineP,error,failure)
       ENDIF
    ENDIF

    IF (PRESENT(use_virial)) THEN
       my_use_virial=use_virial
    ELSE
       my_use_virial=.FALSE.
    ENDIF

    ! my_compute_tau defaults to .FALSE.
    ! IF (.true.) it will compute 0.5 * (nabla x_a).(v(r) nabla x_b)
    IF (PRESENT(compute_tau)) THEN
       my_compute_tau=compute_tau
    ELSE
       my_compute_tau=.FALSE.
    ENDIF

    ! use identical radii for integrate and collocate ?
    IF (PRESENT(map_consistent)) THEN
       my_map_consistent=map_consistent
    ELSE
       my_map_consistent=.FALSE.
    ENDIF

    IF (PRESENT(collocate_rho0).AND.PRESENT(rpgf0_s)) THEN
       my_collocate_rho0=collocate_rho0
    ELSE
       my_collocate_rho0=.FALSE.
    END IF

    IF (calculate_forces) THEN
      la_max_local=la_max+1  ! needed for the derivative of the gaussian, unimportant which one
      la_min_local=MAX(la_min-1,0) ! just in case the la_min,lb_min is not zero
      lb_min_local=MAX(lb_min-1,0)
      lb_max_local=lb_max
      IF (my_use_virial) THEN
         la_max_local=la_max_local+1
         lb_max_local=lb_max_local+1
      ENDIF
    ELSE
      la_max_local=la_max
      la_min_local=la_min
      lb_min_local=lb_min
      lb_max_local=lb_max
    END IF

    IF (my_compute_tau) THEN
      la_max_local=la_max_local+1
      lb_max_local=lb_max_local+1
      la_min_local=MAX(la_min_local-1,0)
      lb_min_local=MAX(lb_min_local-1,0)
    ENDIF

    coef_max=la_max_local+lb_max_local+1
    zetp = zeta + zetb
    f = zetb/zetp
    prefactor = EXP(-zeta*f*rab2)
!   *** position of the gaussian product
    rap(:) = f*rab(:)
    rbp(:) = rap(:) - rab(:)
    rp(:) = ra(:) + rap(:)  ! this is the gaussian center in real coordinates
    rb(:) = ra(:) + rab(:)

    IF (my_map_consistent) THEN ! still assumes that eps_gvg_rspace=eps_rho_rspace
       cutoff=1.0_dp
       radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
               zetp=zetp,eps=eps_gvg_rspace,prefactor=prefactor,cutoff=cutoff)
    ELSE IF (my_collocate_rho0) THEN
       cutoff    = 0.0_dp
       prefactor = 1.0_dp
       radius = rpgf0_s
    ELSE
       cutoff=1.0_dp
       IF (PRESENT(pab)) THEN
          radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,pab,o1,o2,ra,rb,rp,&
                                       zetp,eps_gvg_rspace,prefactor,cutoff)
       ELSE
          radius=exp_radius_very_extended(la_min,la_max,lb_min,lb_max,ra=ra,rb=rb,rp=rp,&
               zetp=zetp,eps=eps_gvg_rspace,prefactor=prefactor,cutoff=cutoff)
       ENDIF
    ENDIF

    IF (radius == 0.0_dp) THEN
       RETURN
    ENDIF

    ng(:) = rsgrid%desc%npts(:)
    grid => rsgrid%r(:,:,:)
    ALLOCATE(vab(ncoset(la_max_local),ncoset(lb_max_local)))
    vab=0.0_dp

    ! the likely call is integrate_ortho
    IF (.NOT. subpatch_integrate) THEN
        IF (rsgrid%desc%orthorhombic ) THEN
          CALL integrate_ortho()
        ELSE
          CALL integrate_general_wings()
        END IF
    ELSE
        CALL integrate_general_subpatch()
    END IF

!   *** vab contains all the information needed to find the elements of hab
!   *** and optionally of derivatives of these elements

    ftza = 2.0_dp*zeta
    ftzb = 2.0_dp*zetb

    DO la=la_min,la_max
      DO ax=0,la
        DO ay=0,la-ax
          az = la - ax - ay
          ico=coset(ax,ay,az)
          DO lb=lb_min,lb_max
            DO bx=0,lb
              DO by=0,lb-bx
                bz = lb - bx - by
                jco=coset(bx,by,bz)
                IF (.NOT.my_compute_tau) THEN
                    axpm0 = vab(ico,jco)
                ELSE
                    axpm0 =  0.5_dp * ( ax * bx * vab(coset(MAX(ax-1,0),ay,az),coset(MAX(bx-1,0),by,bz)) +  &
                                        ay * by * vab(coset(ax,MAX(ay-1,0),az),coset(bx,MAX(by-1,0),bz)) +  &
                                        az * bz * vab(coset(ax,ay,MAX(az-1,0)),coset(bx,by,MAX(bz-1,0)))  &
                                        - ftza * bx * vab(coset(ax+1,ay,az),coset(MAX(bx-1,0),by,bz))  &
                                        - ftza * by * vab(coset(ax,ay+1,az),coset(bx,MAX(by-1,0),bz))  &
                                        - ftza * bz * vab(coset(ax,ay,az+1),coset(bx,by,MAX(bz-1,0)))  &
                                        - ax * ftzb * vab(coset(MAX(ax-1,0),ay,az),coset(bx+1,by,bz))  &
                                        - ay * ftzb * vab(coset(ax,MAX(ay-1,0),az),coset(bx,by+1,bz))  &
                                        - az * ftzb * vab(coset(ax,ay,MAX(az-1,0)),coset(bx,by,bz+1)) +  &
                                        ftza * ftzb * vab(coset(ax+1,ay,az),coset(bx+1,by,bz)) + &
                                        ftza * ftzb * vab(coset(ax,ay+1,az),coset(bx,by+1,bz)) + &
                                        ftza * ftzb * vab(coset(ax,ay,az+1),coset(bx,by,bz+1)) )
                ENDIF
                hab(o1+ico,o2+jco) = hab(o1+ico,o2+jco) + axpm0
                IF (calculate_forces .AND. PRESENT(force_a)) THEN
                  IF (my_compute_tau) THEN
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*ax * bx
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,MAX(ax-1,0),ay,az,MAX(bx-1,0),by,bz,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*ay * by
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,MAX(ay-1,0),az,bx,MAX(by-1,0),bz,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*az * bz
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,MAX(az-1,0),bx,by,MAX(bz-1,0),vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ftza * bx )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax+1,ay,az,MAX(bx-1,0),by,bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ftza * by )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay+1,az,bx,MAX(by-1,0),bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ftza * bz  )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,az+1,bx,by,MAX(bz-1,0) ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ax * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,MAX(ax-1,0),ay,az,bx+1,by,bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- ay * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,MAX(ay-1,0),az,bx,by+1,bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(- az * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,MAX(az-1,0),bx,by,bz+1 ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(ftza * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax+1,ay,az,bx+1,by,bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(ftza * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay+1,az,bx,by+1,bz ,vab)
                     pabval=pab(o1+ico,o2+jco)*0.5_dp*(ftza * ftzb )
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,az+1,bx,by,bz+1 ,vab)
                  ELSE
                     pabval=pab(o1+ico,o2+jco)
                     CALL force_update(force_a,force_b,rab,pabval,ftza,ftzb,ax,ay,az,bx,by,bz,vab)
                     IF (my_use_virial) THEN
                       CALL virial_update(my_virial_a,my_virial_b,rab,pabval,ftza,ftzb,ax,ay,az,bx,by,bz,vab)
                     ENDIF
                  ENDIF
                END IF
                IF (calculate_forces .AND. PRESENT(hdab)) THEN
                  der_a(1:3) = 0.0_dp
                  der_b(1:3) = 0.0_dp
                  CALL hab_derivatives(der_a,der_b,rab,ftza,ftzb,ax,ay,az,bx,by,bz,vab)
                  hdab(1:3,o1+ico,o2+jco) = der_a(1:3)
                  hadb(1:3,o1+ico,o2+jco) = der_b(1:3)
                  pabval=1.0_dp
                  IF (my_use_virial .AND. PRESENT(a_hdab)) THEN
                    my_virial_a =0.0_dp
                    my_virial_b =0.0_dp
                    CALL virial_update(my_virial_a,my_virial_b,rab,pabval,ftza,ftzb,ax,ay,az,bx,by,bz,vab)
                    DO j=1,3
                      a_hdab(1:3,j,o1+ico,o2+jco) = a_hdab(1:3,j,o1+ico,o2+jco) + my_virial_a(1:3,j)
                    END DO
                  ENDIF
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(vab)


  CONTAINS

! *****************************************************************************
   SUBROUTINE integrate_ortho()

    CALL return_cube(cube_info,radius,lb_cube,ub_cube,sphere_bounds)
    cmax=MAXVAL(ub_cube)

    dr(1) = rsgrid%desc%dh(1,1)
    dr(2) = rsgrid%desc%dh(2,2)
    dr(3) = rsgrid%desc%dh(3,3)

    gridbounds(1,1)=LBOUND(GRID,1)
    gridbounds(2,1)=UBOUND(GRID,1)
    gridbounds(1,2)=LBOUND(GRID,2)
    gridbounds(2,2)=UBOUND(GRID,2)
    gridbounds(1,3)=LBOUND(GRID,3)
    gridbounds(2,3)=UBOUND(GRID,3)

    CALL compute_cube_center(cubecenter,rsgrid%desc,zeta,zetb,ra,rab)
    roffset(:) = rp(:) - REAL(cubecenter(:),dp)*dr(:)
    lb_cube_min = MINVAL(lb_cube(:))
    ub_cube_max = MAXVAL(ub_cube(:))

!   *** a mapping so that the ig corresponds to the right grid point, also with pbc
    ALLOCATE(map(-cmax:cmax,3))
    DO i=1,3
      IF ( rsgrid % desc % perd ( i ) == 1 ) THEN
        start=lb_cube(i)
        DO
         offset=MODULO(cubecenter(i)+start,ng(i))+1-start
         length=MIN(ub_cube(i),ng(i)-offset)-start
         DO ig=start,start+length
            map(ig,i) = ig+offset
         END DO
         IF (start+length.GE.ub_cube(i)) EXIT
         start=start+length+1
        END DO
      ELSE
        ! this takes partial grid + border regions into account
        offset=MODULO(cubecenter(i)+lb_cube(i)+rsgrid%desc%lb(i)-rsgrid%lb_local(i),ng(i))+1-lb_cube(i)
        ! check for out of bounds
        IF (ub_cube(i)+offset>UBOUND(grid,i).OR.lb_cube(i)+offset<LBOUND(grid,i)) THEN
           CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        ENDIF
        DO ig=lb_cube(i),ub_cube(i)
           map(ig,i) = ig+offset
        END DO
      END IF
    ENDDO

    lp=la_max_local+lb_max_local
    ALLOCATE(coef_xyz(((lp+1)*(lp+2)*(lp+3))/6))
    ALLOCATE(pol_z(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_y(1:2,0:lp,-cmax:0))
    ALLOCATE(pol_x(0:lp,-cmax:cmax))

#include "prep.f90"

    CALL call_integrate()


    CALL call_to_xyz_to_vab(prefactor, coef_xyz, lp, la_max_local, lb_max_local, &
                                  rp, ra, rab, vab, coset, la_min_local, & 
                                  lb_min_local, current_maxl, ncoset(la_max_local), ncoset(lb_max_local))
    
    DEALLOCATE(coef_xyz)
    DEALLOCATE(pol_z)
    DEALLOCATE(pol_y)
    DEALLOCATE(pol_x)
    DEALLOCATE(map)
    END SUBROUTINE integrate_ortho

! *****************************************************************************
    SUBROUTINE call_integrate


      SELECT CASE(lp)
      CASE(0)
      CALL integrate_core_0(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                             map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
      CASE(1)
      CALL integrate_core_1(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                             map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
      CASE(2)
      CALL integrate_core_2(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                             map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
      CASE(3)
      CALL integrate_core_3(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                             map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
      CASE(4)
      CALL integrate_core_4(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                             map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
      CASE(5)
      CALL integrate_core_5(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                             map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
      CASE(6)
      CALL integrate_core_6(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                             map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
      CASE(7)
      CALL integrate_core_7(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                             map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
      CASE(8)
      CALL integrate_core_8(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                             map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
      CASE(9)
      CALL integrate_core_9(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                             map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
      CASE DEFAULT
      CALL integrate_core_default(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax),pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                             map(-cmax,1),sphere_bounds(1),lp,cmax,gridbounds(1,1))
      END SELECT



    END SUBROUTINE call_integrate

       

! *****************************************************************************
! general optimized routine, not called for orthogonal cells
! *****************************************************************************
    SUBROUTINE integrate_general_opt()
    INTEGER :: i, i_index, il, ilx, ily, ilz, index_max(3), index_min(3), &
      ismax, ismin, j, j_index, jl, jlx, jly, jlz, k, k_index, kl, klx, kly, &
      klz, lpx, lpy, lpz, lx, ly, lz, offset(3)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: grid_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: coef_map
    REAL(KIND=dp) :: a, b, c, d, di, dip, dj, djp, dk, dkp, exp0i, exp1i, &
      exp2i, gp(3), gridval, hmatgrid(3,3), pointj(3), pointk(3), v(3)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: coef_ijk
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: hmatgridp

!
! transform P_{lxp,lyp,lzp} into a P_{lip,ljp,lkp} such that
! sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-x_p)**lxp (y-y_p)**lyp (z-z_p)**lzp =
! sum_{lip,ljp,lkp} P_{lip,ljp,lkp} (i-i_p)**lip (j-j_p)**ljp (k-k_p)**lkp
!

      lp=la_max_local+lb_max_local
      ALLOCATE(coef_xyz(((lp+1)*(lp+2)*(lp+3))/6))
      ALLOCATE(coef_ijk(((lp+1)*(lp+2)*(lp+3))/6))
      ALLOCATE(coef_xyt(((lp+1)*(lp+2))/2))
      ALLOCATE(coef_xtt(0:lp))

      ! aux mapping array to simplify life
      ALLOCATE(coef_map(0:lp,0:lp,0:lp))
      coef_map=HUGE(coef_map)
      lxyz=0
      DO lzp=0,lp
      DO lyp=0,lp-lzp
      DO lxp=0,lp-lzp-lyp
          lxyz=lxyz+1
          coef_map(lxp,lyp,lzp)=lxyz
      ENDDO
      ENDDO
      ENDDO

      ! cell hmat in grid points
      hmatgrid(:,1)=cell%hmat(:,1)/ng(1)
      hmatgrid(:,2)=cell%hmat(:,2)/ng(2)
      hmatgrid(:,3)=cell%hmat(:,3)/ng(3)

      ! center in grid coords
      gp=MATMUL(cell%h_inv,rp)*ng

      ! added bt matt
      cubecenter(:) = FLOOR(gp)

      !t2=nanotime_ia32()
      !write(6,*) t2-t1
      !t1=nanotime_ia32()

      CALL return_cube_nonortho(cube_info,radius,index_min,index_max,rp)

      offset(:)=MODULO(index_min(:)+rsgrid%desc%lb(:)-rsgrid%lb_local(:),ng(:))+1

      ALLOCATE(grid_map(index_min(1):index_max(1)))
      DO i=index_min(1),index_max(1)
         grid_map(i)=MODULO(i,ng(1))+1
         IF (rsgrid % desc % perd ( 1 )==1) THEN
            grid_map(i)=MODULO(i,ng(1))+1
         ELSE
            grid_map(i)=i-index_min(1)+offset(1)
         ENDIF
      ENDDO


      coef_ijk=0.0_dp

      ! go over the grid, but cycle if the point is not within the radius
      DO k=index_min(3),index_max(3)
        dk=k-gp(3)
        pointk=hmatgrid(:,3)*dk

        ! allow for generalised rs grids
        IF (rsgrid % desc % perd ( 3 )==1) THEN
           k_index=MODULO(k,ng(3))+1
        ELSE
           k_index=k-index_min(3)+offset(3)
        ENDIF

        coef_xyt=0.0_dp

        DO j=index_min(2),index_max(2)
          dj=j-gp(2)
          pointj=pointk+hmatgrid(:,2)*dj
          IF (rsgrid % desc % perd ( 2 )==1) THEN
             j_index=MODULO(j,ng(2))+1
          ELSE
             j_index=j-index_min(2)+offset(2)
          ENDIF

          coef_xtt=0.0_dp

          ! find bounds for the inner loop
          ! based on a quadratic equation in i
          ! a*i**2+b*i+c=radius**2
          v=pointj-gp(1)*hmatgrid(:,1)
          a=DOT_PRODUCT(hmatgrid(:,1),hmatgrid(:,1))
          b=2*DOT_PRODUCT(v,hmatgrid(:,1))
          c=DOT_PRODUCT(v,v)
          d=b*b-4*a*(c-radius**2)

          IF (d<0) THEN
              CYCLE
          ELSE
              d=SQRT(d)
              ismin=CEILING((-b-d)/(2*a))
              ismax=FLOOR((-b+d)/(2*a))
          ENDIF
          ! prepare for computing -zetp*rsq
          a=-zetp*a
          b=-zetp*b
          c=-zetp*c
          i=ismin-1
          exp2i=EXP((a*i+b)*i+c)
          exp1i=EXP(2*a*i+a+b)
          exp0i=EXP(2*a)

          coef_xtt=0.0_dp

          DO i=ismin,ismax
             di=i-gp(1)

             exp2i=exp2i*exp1i
             exp1i=exp1i*exp0i

             i_index=grid_map(i)
             gridval=grid(i_index,j_index,k_index)*exp2i

             dip=1.0_dp
             DO il=0,lp
                coef_xtt(il)=coef_xtt(il)+gridval*dip
                dip=dip*di
             ENDDO
          ENDDO

          lxy=0
          djp=1.0_dp
          DO jl=0,lp
            DO il=0,lp-jl
               lxy=lxy+1
               coef_xyt(lxy)=coef_xyt(lxy)+coef_xtt(il)*djp
            ENDDO
            djp=djp*dj
          ENDDO

        ENDDO

        lxyz = 0
        dkp=1.0_dp
        DO kl=0,lp
           lxy=0
           DO jl=0,lp-kl
              DO il=0,lp-kl-jl
                 lxyz=lxyz+1 ; lxy=lxy+1
                 coef_ijk(lxyz)=coef_ijk(lxyz)+dkp*coef_xyt(lxy)
              ENDDO
              lxy=lxy+kl
           ENDDO
           dkp=dkp*dk
        ENDDO

      ENDDO

      ! transform using multinomials
      ALLOCATE(hmatgridp(3,3,0:lp))
      hmatgridp(:,:,0)=1.0_dp
      DO k=1,lp
         hmatgridp(:,:,k)=hmatgridp(:,:,k-1)*hmatgrid(:,:)
      ENDDO

      coef_xyz=0.0_dp
      lpx=lp
      DO klx=0,lpx
      DO jlx=0,lpx-klx
      DO ilx=0,lpx-klx-jlx
         lx=ilx+jlx+klx
         lpy=lp-lx
         DO kly=0,lpy
         DO jly=0,lpy-kly
         DO ily=0,lpy-kly-jly
            ly=ily+jly+kly
            lpz=lp-lx-ly
            DO klz=0,lpz
            DO jlz=0,lpz-klz
            DO ilz=0,lpz-klz-jlz
               lz=ilz+jlz+klz

               il=ilx+ily+ilz
               jl=jlx+jly+jlz
               kl=klx+kly+klz
               coef_xyz(coef_map(lx,ly,lz))=coef_xyz(coef_map(lx,ly,lz))+ coef_ijk(coef_map(il,jl,kl))* &
                                            hmatgridp(1,1,ilx) * hmatgridp(1,2,jlx) * hmatgridp(1,3,klx) * &
                                            hmatgridp(2,1,ily) * hmatgridp(2,2,jly) * hmatgridp(2,3,kly) * &
                                            hmatgridp(3,1,ilz) * hmatgridp(3,2,jlz) * hmatgridp(3,3,klz) * &
                                            fac(lx)*fac(ly)*fac(lz)/ &
                        (fac(ilx)*fac(ily)*fac(ilz)*fac(jlx)*fac(jly)*fac(jlz)*fac(klx)*fac(kly)*fac(klz))
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL call_to_xyz_to_vab(prefactor, coef_xyz, lp, la_max_local, lb_max_local, &
                                  rp, ra, rab, vab, coset, la_min_local, & 
                                  lb_min_local, current_maxl, ncoset(la_max_local), & 
                                  ncoset(lb_max_local))


      ! deallocation needed to pass around a pgi bug..
      DEALLOCATE(hmatgridp)
      DEALLOCATE(grid_map)
      DEALLOCATE(coef_map)
      DEALLOCATE(coef_xtt)
      DEALLOCATE(coef_xyt)
      DEALLOCATE(coef_ijk)
      DEALLOCATE(coef_xyz)
    END SUBROUTINE integrate_general_opt

! *****************************************************************************

    SUBROUTINE integrate_general_subpatch
    INTEGER                                  :: stat
    INTEGER, DIMENSION(2, 3)                 :: local_b
    INTEGER, DIMENSION(3)                    :: local_s, periodic
    REAL(dp), DIMENSION((&
      la_max_local+lb_max_local+1)*(&
      la_max_local+lb_max_local+2)*(&
      la_max_local+lb_max_local+3)/6)        :: poly_d3

        periodic=1 ! cell%perd
        lp=la_max_local+lb_max_local
        local_b(1,:)=rsgrid%lb_real-rsgrid%desc%lb
        local_b(2,:)=rsgrid%ub_real-rsgrid%desc%lb
        local_s=rsgrid%lb_real-rsgrid%lb_local
        IF (BTEST(subpatch_pattern,0)) local_b(1,1)=local_b(1,1)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,1)) local_b(2,1)=local_b(2,1)+rsgrid%desc%border
        IF (BTEST(subpatch_pattern,2)) local_b(1,2)=local_b(1,2)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,3)) local_b(2,2)=local_b(2,2)+rsgrid%desc%border
        IF (BTEST(subpatch_pattern,4)) local_b(1,3)=local_b(1,3)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,5)) local_b(2,3)=local_b(2,3)+rsgrid%desc%border
        IF (BTEST(subpatch_pattern,0)) local_s(1)=local_s(1)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,2)) local_s(2)=local_s(2)-rsgrid%desc%border
        IF (BTEST(subpatch_pattern,4)) local_s(3)=local_s(3)-rsgrid%desc%border
        CALL integrateGaussFull(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rp,max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b,local_shift=local_s,&
            error=error,scale=rsgrid%desc%ngpts/ABS(cell%deth))
        ! defaults: local_shift=(/0,0,0/),poly_shift=(/0.0_dp,0.0_dp,0.0_dp/),scale=1.0_dp,
        ALLOCATE(coef_xyz(((lp+1)*(lp+2)*(lp+3))/6),stat=stat)
        CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
        CALL poly_d32cp2k(coef_xyz,lp,poly_d3,error)
        CALL call_to_xyz_to_vab(prefactor, coef_xyz, lp, la_max_local, lb_max_local, &
                                  rp, ra, rab, vab, coset, la_min_local, lb_min_local, &
                                  current_maxl, ncoset(la_max_local), & 
                                  ncoset(lb_max_local))

        DEALLOCATE(coef_xyz,stat=stat)
        CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    END SUBROUTINE

    SUBROUTINE integrate_general_wings()
    INTEGER                                  :: periodic(3), stat
    INTEGER, DIMENSION(2, 3)                 :: local_b
    REAL(dp), DIMENSION((&
      la_max_local+lb_max_local+1)*(&
      la_max_local+lb_max_local+2)*(&
      la_max_local+lb_max_local+3)/6)        :: poly_d3
    REAL(dp), DIMENSION(3)                   :: local_shift, rShifted

        periodic=1 ! cell%perd
        local_b(1,:)=0
        local_b(2,:)=MIN(rsgrid%desc%npts-1,rsgrid%ub_local-rsgrid%lb_local)
        local_shift=REAL(rsgrid%desc%lb-rsgrid%lb_local,dp)/REAL(rsgrid%desc%npts,dp)
        rShifted(1)=rp(1)+cell%hmat(1,1)*local_shift(1)&
             +cell%hmat(1,2)*local_shift(2)&
             +cell%hmat(1,3)*local_shift(3)
        rShifted(2)=rp(2)+cell%hmat(2,1)*local_shift(1)&
             +cell%hmat(2,2)*local_shift(2)&
             +cell%hmat(2,3)*local_shift(3)
        rShifted(3)=rp(3)+cell%hmat(3,1)*local_shift(1)&
             +cell%hmat(3,2)*local_shift(2)&
             +cell%hmat(3,3)*local_shift(3)
        lp=la_max_local+lb_max_local
        CALL integrateGaussFull(h=cell%hmat,h_inv=cell%h_inv,&
            grid=grid,poly=poly_d3,alphai=zetp,posi=rShifted,&
            max_r2=radius*radius,&
            periodic=periodic,gdim=ng,local_bounds=local_b,&
            error=error,scale=rsgrid%desc%ngpts/ABS(cell%deth))
        ! defaults: local_shift=(/0,0,0/),poly_shift=(/0.0_dp,0.0_dp,0.0_dp/),scale=1.0_dp,
        ALLOCATE(coef_xyz(((lp+1)*(lp+2)*(lp+3))/6),stat=stat)
        CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
        CALL poly_d32cp2k(coef_xyz,lp,poly_d3,error)
        CALL call_to_xyz_to_vab(prefactor, coef_xyz, lp, la_max_local, lb_max_local, &
                                  rp, ra, rab, vab, coset, la_min_local, lb_min_local, &
                                current_maxl, ncoset(la_max_local), ncoset(lb_max_local))


        DEALLOCATE(coef_xyz,stat=stat)
        CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    END SUBROUTINE

! *****************************************************************************
    SUBROUTINE integrate_general()
    INTEGER                                  :: i, index_max(3), &
                                                index_min(3), ipoint(3), j, k
    REAL(KIND=dp)                            :: gridval, point(3)

      CALL return_cube_nonortho(cube_info,radius,index_min,index_max,rp)

      ! go over the grid, but cycle if the point is not within the radius
      DO k=index_min(3),index_max(3)
      DO j=index_min(2),index_max(2)
      DO i=index_min(1),index_max(1)
         ! point in real space
         point=MATMUL(cell%hmat,REAL((/i,j,k/),KIND=dp)/ng)
         ! skip if outside of the sphere
         IF (SUM((point-rp)**2)>radius**2) CYCLE
         ! point on the grid (including pbc)
         ipoint=MODULO((/i,j,k/),ng)+1
         ! integrate on the grid
         gridval=grid(ipoint(1),ipoint(2),ipoint(3))
         CALL primitive_integrate(point,gridval)
      ENDDO
      ENDDO
      ENDDO
    END SUBROUTINE integrate_general

! *****************************************************************************
    SUBROUTINE primitive_integrate(point,gridval)
    REAL(KIND=dp)                            :: point(3), gridval

    REAL(KIND=dp)                            :: dra(3), drap(3), drb(3), &
                                                drbp(3), myexp

       myexp=EXP(-zetp*SUM((point-rp)**2))*prefactor*gridval
        dra=point-ra
        drb=point-rb
        drap(1)=1.0_dp
        DO lxa=0,la_max_local
        drbp(1)=1.0_dp
        DO lxb=0,lb_max_local
           drap(2)=1.0_dp
           DO lya=0,la_max_local-lxa
           drbp(2)=1.0_dp
           DO lyb=0,lb_max_local-lxb
              drap(3)=1.0_dp
              DO lza=1,MAX(la_min_local-lxa-lya,0)
                 drap(3)=drap(3)*dra(3)
              ENDDO
              DO lza=MAX(la_min_local-lxa-lya,0),la_max_local-lxa-lya
              drbp(3)=1.0_dp
              DO lzb=1,MAX(lb_min_local-lxb-lyb,0)
                 drbp(3)=drbp(3)*drb(3)
              ENDDO
              DO lzb=MAX(lb_min_local-lxb-lyb,0),lb_max_local-lxb-lyb
                ico=coset(lxa,lya,lza)
                jco=coset(lxb,lyb,lzb)
                vab(ico,jco)=vab(ico,jco)+myexp*PRODUCT(drap)*PRODUCT(drbp)
                drbp(3)=drbp(3)*drb(3)
              ENDDO
              drap(3)=drap(3)*dra(3)
              ENDDO
           drbp(2)=drbp(2)*drb(2)
           ENDDO
           drap(2)=drap(2)*dra(2)
           ENDDO
        drbp(1)=drbp(1)*drb(1)
        ENDDO
        drap(1)=drap(1)*dra(1)
        ENDDO

    END SUBROUTINE


  END SUBROUTINE integrate_pgf_product_rspace

! *****************************************************************************
!> \brief given a set of matrix elements, perform the correct contraction to obtain the virial
! *****************************************************************************
  SUBROUTINE virial_update(my_virial_a,my_virial_b,rab,pab,&
                           ftza,ftzb,ax,ay,az,bx,by,bz,vab)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(INOUT)                          :: my_virial_a, my_virial_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: pab, ftza, ftzb
    INTEGER, INTENT(IN)                      :: ax, ay, az, bx, by, bz
    REAL(KIND=dp)                            :: vab(:,:)

    my_virial_a(1,1) = my_virial_a(1,1) &
        + pab*ftza*vab(coset(ax+2,ay,az),coset(bx,by,bz)) &
        - pab*REAL(ax,dp)*vab(coset(MAX(0,ax-1)+1,ay,az),coset(bx,by,bz))
    my_virial_a(1,2) = my_virial_a(1,2) &
        + pab*ftza*vab(coset(ax+1,ay+1,az),coset(bx,by,bz)) &
        - pab*REAL(ax,dp)*vab(coset(MAX(0,ax-1),ay+1,az),coset(bx,by,bz))
    my_virial_a(1,3) = my_virial_a(1,3) &
        + pab*ftza*vab(coset(ax+1,ay,az+1),coset(bx,by,bz)) &
        - pab*REAL(ax,dp)*vab(coset(MAX(0,ax-1),ay,az+1),coset(bx,by,bz))
    my_virial_a(2,1) = my_virial_a(2,1) &
        + pab*ftza*vab(coset(ax+1,ay+1,az),coset(bx,by,bz)) &
        - pab*REAL(ay,dp)*vab(coset(ax+1,MAX(0,ay-1),az),coset(bx,by,bz))
    my_virial_a(2,2) = my_virial_a(2,2) &
        + pab*ftza*vab(coset(ax,ay+2,az),coset(bx,by,bz)) &
        - pab*REAL(ay,dp)*vab(coset(ax,MAX(0,ay-1)+1,az),coset(bx,by,bz))
    my_virial_a(2,3) = my_virial_a(2,3) &
        + pab*ftza*vab(coset(ax,ay+1,az+1),coset(bx,by,bz)) &
        - pab*REAL(ay,dp)*vab(coset(ax,MAX(0,ay-1),az+1),coset(bx,by,bz))
    my_virial_a(3,1) = my_virial_a(3,1) &
        + pab*ftza*vab(coset(ax+1,ay,az+1),coset(bx,by,bz)) &
        - pab*REAL(az,dp)*vab(coset(ax+1,ay,MAX(0,az-1)),coset(bx,by,bz))
    my_virial_a(3,2) = my_virial_a(3,2) &
        + pab*ftza*vab(coset(ax,ay+1,az+1),coset(bx,by,bz)) &
        - pab*REAL(az,dp)*vab(coset(ax,ay+1,MAX(0,az-1)),coset(bx,by,bz))
    my_virial_a(3,3) = my_virial_a(3,3) &
        + pab*ftza*vab(coset(ax,ay,az+2),coset(bx,by,bz)) &
        - pab*REAL(az,dp)*vab(coset(ax,ay,MAX(0,az-1)+1),coset(bx,by,bz))

    my_virial_b(1,1) = my_virial_b(1,1) + pab*ftzb* ( &
          vab(coset(ax+2,ay,az),coset(bx,by,bz)) &
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(1) &
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(1) &
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(1)*rab(1) ) &
        - pab*REAL(bx,dp)*vab(coset(ax,ay,az),coset(MAX(0,bx-1)+1,by,bz))
    my_virial_b(1,2) = my_virial_b(1,2) + pab*ftzb* ( &
          vab(coset(ax+1,ay+1,az),coset(bx,by,bz)) &
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(1) &
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(2) &
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(1)*rab(2) ) &
        - pab*REAL(bx,dp)*vab(coset(ax,ay,az),coset(MAX(0,bx-1),by+1,bz))
    my_virial_b(1,3) = my_virial_b(1,3) + pab*ftzb* ( &
          vab(coset(ax+1,ay,az+1),coset(bx,by,bz)) &
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(1) &
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(3) &
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(1)*rab(3) ) &
        - pab*REAL(bx,dp)*vab(coset(ax,ay,az),coset(MAX(0,bx-1),by,bz+1))
    my_virial_b(2,1) = my_virial_b(2,1) + pab*ftzb* ( &
          vab(coset(ax+1,ay+1,az),coset(bx,by,bz)) &
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(2) &
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(1) &
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(2)*rab(1) ) &
        - pab*REAL(by,dp)*vab(coset(ax,ay,az),coset(bx+1,MAX(0,by-1),bz))
    my_virial_b(2,2) = my_virial_b(2,2) + pab*ftzb* ( &
          vab(coset(ax,ay+2,az),coset(bx,by,bz)) &
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(2) &
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(2) &
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(2)*rab(2) ) &
        - pab*REAL(by,dp)*vab(coset(ax,ay,az),coset(bx,MAX(0,by-1)+1,bz))
    my_virial_b(2,3) = my_virial_b(2,3) + pab*ftzb* ( &
          vab(coset(ax,ay+1,az+1),coset(bx,by,bz)) &
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(2) &
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(3) &
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(2)*rab(3) ) &
        - pab*REAL(by,dp)*vab(coset(ax,ay,az),coset(bx,MAX(0,by-1),bz+1))
    my_virial_b(3,1) = my_virial_b(3,1) + pab*ftzb* ( &
          vab(coset(ax+1,ay,az+1),coset(bx,by,bz)) &
        - vab(coset(ax+1,ay,az),coset(bx,by,bz))*rab(3) &
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(1) &
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(3)*rab(1) ) &
        - pab*REAL(bz,dp)*vab(coset(ax,ay,az),coset(bx+1,by,MAX(0,bz-1)))
    my_virial_b(3,2) = my_virial_b(3,2) + pab*ftzb* ( &
          vab(coset(ax,ay+1,az+1),coset(bx,by,bz)) &
        - vab(coset(ax,ay+1,az),coset(bx,by,bz))*rab(3) &
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(2) &
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(3)*rab(2) ) &
        - pab*REAL(bz,dp)*vab(coset(ax,ay,az),coset(bx,by+1,MAX(0,bz-1)))
    my_virial_b(3,3) = my_virial_b(3,3) + pab*ftzb* ( &
          vab(coset(ax,ay,az+2),coset(bx,by,bz)) &
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(3) &
        - vab(coset(ax,ay,az+1),coset(bx,by,bz))*rab(3) &
        + vab(coset(ax,ay,az),coset(bx,by,bz))*rab(3)*rab(3) ) &
        - pab*REAL(bz,dp)*vab(coset(ax,ay,az),coset(bx,by,MAX(0,bz-1)+1))

  END SUBROUTINE virial_update

! *****************************************************************************
!> \brief given a bunch of matrix elements, performe the right contractions to obtain the forces
! *****************************************************************************
  SUBROUTINE force_update(force_a,force_b,rab,pab,ftza,ftzb,ax,ay,az,bx,by,bz,vab)
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force_a, force_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: pab, ftza, ftzb
    INTEGER, INTENT(IN)                      :: ax, ay, az, bx, by, bz
    REAL(KIND=dp)                            :: vab(:,:)

    REAL(KIND=dp)                            :: axm1, axp1, axpm0, aym1, &
                                                ayp1, azm1, azp1, bxm1, bym1, &
                                                bzm1

    axpm0 = vab(coset(ax,ay,az),coset(bx,by,bz))
    axp1=vab(coset(ax+1,ay,az),coset(bx,by,bz))
    axm1=vab(coset(MAX(0,ax-1),ay,az),coset(bx,by,bz))
    ayp1=vab(coset(ax,ay+1,az),coset(bx,by,bz))
    aym1=vab(coset(ax,MAX(0,ay-1),az),coset(bx,by,bz))
    azp1=vab(coset(ax,ay,az+1),coset(bx,by,bz))
    azm1=vab(coset(ax,ay,MAX(0,az-1)),coset(bx,by,bz))
    bxm1=vab(coset(ax,ay,az),coset(MAX(0,bx-1),by,bz))
    bym1=vab(coset(ax,ay,az),coset(bx,MAX(0,by-1),bz))
    bzm1=vab(coset(ax,ay,az),coset(bx,by,MAX(0,bz-1)))
    force_a(1) = force_a(1) + pab*(ftza*axp1 - REAL(ax,dp)* axm1)
    force_a(2) = force_a(2) + pab*(ftza*ayp1 - REAL(ay,dp)* aym1)
    force_a(3) = force_a(3) + pab*(ftza*azp1 - REAL(az,dp)* azm1)
    force_b(1) = force_b(1) + pab*(ftzb*(axp1 - rab(1)*axpm0) - REAL(bx,dp)* bxm1)
    force_b(2) = force_b(2) + pab*(ftzb*(ayp1 - rab(2)*axpm0) - REAL(by,dp)* bym1)
    force_b(3) = force_b(3) + pab*(ftzb*(azp1 - rab(3)*axpm0) - REAL(bz,dp)* bzm1)

  END SUBROUTINE force_update

! *****************************************************************************
!> \brief given a bunch of matrix elements perform the right contractions to obtain the
!>      derivatives of the hab matirx
! *****************************************************************************
  SUBROUTINE hab_derivatives(der_a,der_b,rab,ftza,ftzb,ax,ay,az,bx,by,bz,vab)
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: der_a, der_b
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: ftza, ftzb
    INTEGER, INTENT(IN)                      :: ax, ay, az, bx, by, bz
    REAL(KIND=dp)                            :: vab(:,:)

    REAL(KIND=dp)                            :: axm1, axp1, axpm0, aym1, &
                                                ayp1, azm1, azp1, bxm1, bym1, &
                                                bzm1

    axpm0 = vab(coset(ax,ay,az),coset(bx,by,bz))
    axp1=vab(coset(ax+1,ay,az),coset(bx,by,bz))
    axm1=vab(coset(MAX(0,ax-1),ay,az),coset(bx,by,bz))
    ayp1=vab(coset(ax,ay+1,az),coset(bx,by,bz))
    aym1=vab(coset(ax,MAX(0,ay-1),az),coset(bx,by,bz))
    azp1=vab(coset(ax,ay,az+1),coset(bx,by,bz))
    azm1=vab(coset(ax,ay,MAX(0,az-1)),coset(bx,by,bz))
    bxm1=vab(coset(ax,ay,az),coset(MAX(0,bx-1),by,bz))
    bym1=vab(coset(ax,ay,az),coset(bx,MAX(0,by-1),bz))
    bzm1=vab(coset(ax,ay,az),coset(bx,by,MAX(0,bz-1)))
    der_a(1) =  (ftza*axp1 - REAL(ax,dp)* axm1)
    der_a(2) =  (ftza*ayp1 - REAL(ay,dp)* aym1)
    der_a(3) =  (ftza*azp1 - REAL(az,dp)* azm1)
    der_b(1) =  (ftzb*(axp1 - rab(1)*axpm0) - REAL(bx,dp)* bxm1)
    der_b(2) =  (ftzb*(ayp1 - rab(2)*axpm0) - REAL(by,dp)* bym1)
    der_b(3) =  (ftzb*(azp1 - rab(3)*axpm0) - REAL(bz,dp)* bzm1)

  END SUBROUTINE hab_derivatives

END MODULE qs_integrate_potential_low
