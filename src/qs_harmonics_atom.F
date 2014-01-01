!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE qs_harmonics_atom

  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE lebedev,                         ONLY: lebedev_grid
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indco,&
                                             indso,&
                                             nco,&
                                             ncoset,&
                                             nso,&
                                             nsoset,&
                                             soset
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE spherical_harmonics,             ONLY: y_lm
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_harmonics_atom'

  TYPE harmonics_atom_type
     INTEGER                                :: max_s_harm,llmax, &
                                               max_iso_not0, &
                                               dmax_iso_not0, &
                                               damax_iso_not0, &
                                               ngrid
     REAL(dp), DIMENSION(:,:),    POINTER   :: a, slm
     REAL(dp), DIMENSION(:,:,:),  POINTER   :: dslm, dslm_dxyz
     REAL(dp), DIMENSION(:,:,:),  POINTER   :: my_CG, my_dCG
     REAL(dp), DIMENSION(:,:,:,:),  POINTER :: my_CG_dxyz
     REAL(dp), DIMENSION(:,:,:,:),  POINTER :: my_CG_dxyz_asym
     REAL(dp), DIMENSION(:), POINTER        :: slm_int

  END TYPE harmonics_atom_type

  PUBLIC :: allocate_harmonics_atom, &
            create_harmonics_atom, &
            deallocate_harmonics_atom,&
            get_none0_cg_list

  PUBLIC :: harmonics_atom_type, get_maxl_CG

  INTERFACE get_none0_cg_list
     MODULE PROCEDURE get_none0_cg_list3,get_none0_cg_list4
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief   Allocate a spherical harmonics set for the atom grid.
!> \version 1.0
! *****************************************************************************
  SUBROUTINE allocate_harmonics_atom(harmonics,error)

    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_harmonics_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(harmonics)) CALL deallocate_harmonics_atom(harmonics,error)

    ALLOCATE (harmonics,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    harmonics%max_s_harm=0
    harmonics%llmax=0
    harmonics%max_iso_not0=0
    harmonics%dmax_iso_not0=0
    harmonics%damax_iso_not0=0
    harmonics%ngrid=0

    NULLIFY(harmonics%slm)
    NULLIFY(harmonics%dslm)
    NULLIFY(harmonics%dslm_dxyz)
    NULLIFY(harmonics%slm_int)
    NULLIFY(harmonics%my_CG)
    NULLIFY(harmonics%my_dCG)
    NULLIFY(harmonics%my_CG_dxyz)
    NULLIFY(harmonics%my_CG_dxyz_asym)
    NULLIFY(harmonics%a)

  END SUBROUTINE allocate_harmonics_atom

! *****************************************************************************
!> \brief   Deallocate the spherical harmonics set for the atom grid.
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_harmonics_atom(harmonics,error)

    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_harmonics_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(harmonics)) THEN

       IF(ASSOCIATED(harmonics%slm)) THEN
          DEALLOCATE (harmonics%slm,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       IF(ASSOCIATED(harmonics%dslm)) THEN
          DEALLOCATE (harmonics%dslm,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       IF(ASSOCIATED(harmonics%dslm_dxyz)) THEN
          DEALLOCATE (harmonics%dslm_dxyz,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       IF(ASSOCIATED(harmonics%slm_int)) THEN
          DEALLOCATE (harmonics%slm_int,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       IF(ASSOCIATED(harmonics%my_CG)) THEN
          DEALLOCATE (harmonics%my_CG,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       IF(ASSOCIATED(harmonics%my_dCG)) THEN
          DEALLOCATE (harmonics%my_dCG,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       IF(ASSOCIATED(harmonics%my_CG_dxyz)) THEN
          DEALLOCATE (harmonics%my_CG_dxyz,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       IF(ASSOCIATED(harmonics%my_CG_dxyz_asym)) THEN
          DEALLOCATE (harmonics%my_CG_dxyz_asym,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       IF(ASSOCIATED(harmonics%a)) THEN
          DEALLOCATE (harmonics%a,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ENDIF

       DEALLOCATE (harmonics,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer harmonics is not associated and "//&
            "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_harmonics_atom

! *****************************************************************************
  SUBROUTINE create_harmonics_atom(harmonics,my_CG,na,&
       llmax,maxs,max_s_harm,ll,wa,azi,cos_azi,pol,cos_pol,cotan_azi,error)

    TYPE(harmonics_atom_type), POINTER       :: harmonics
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    INTEGER, INTENT(IN)                      :: na, llmax, maxs, max_s_harm, &
                                                ll
    REAL(dp), DIMENSION(:), INTENT(IN)       :: wa, azi, cos_azi, pol, &
                                                cos_pol, cotan_azi
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_harmonics_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ia, ic, is, is1, &
                                                is2, iso, iso1, iso2, iso_mm, &
                                                iso_mm1, iso_mp1, l, l1, l2, &
                                                lx, ly, lz, m, m1, m2, n
    LOGICAL                                  :: failure
    REAL(dp)                                 :: abs_m, cosmp1, drx, dry, drz, &
                                                int1, int2, int3, l_r, m_r, &
                                                rx, ry, rz, sinmp1
    REAL(dp), DIMENSION(:), POINTER          :: slm_int, y
    REAL(dp), DIMENSION(:, :), POINTER       :: dc, dy, slm
    REAL(dp), DIMENSION(:, :, :), POINTER    :: dslm_dxyz

    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY (y,dy,slm,dslm_dxyz,dc)

    CPPrecondition(ASSOCIATED(harmonics),cp_failure_level,routineP,error,failure)

    harmonics%max_s_harm = max_s_harm
    harmonics%llmax = llmax
    harmonics%ngrid = na

    NULLIFY(harmonics%my_CG,harmonics%my_CG_dxyz)
    CALL reallocate(harmonics%my_CG,1,maxs,1,maxs,1,max_s_harm)
    CALL reallocate(harmonics%my_CG_dxyz,1,3,1,maxs,1,maxs,1,max_s_harm)
    CALL reallocate(harmonics%my_CG_dxyz_asym,1,3,1,maxs,1,maxs,1,max_s_harm)

    DO i = 1,max_s_harm
       DO is1 = 1, maxs
          harmonics%my_CG(1:maxs,is1,i) =  my_CG(1:maxs,is1,i)
       END DO
    END DO

    !   *** allocate and calculate the spherical harmonics LM for this grid ***
    !       and their derivatives
    CALL reallocate(harmonics%slm,1,na,1,max_s_harm)
    CALL reallocate(harmonics%dslm,1,2,1,na,1,max_s_harm)
    CALL reallocate(harmonics%dslm_dxyz,1,3,1,na,1,max_s_harm)
    CALL reallocate(harmonics%a,1,3,1,na)
    CALL reallocate(y,1,na)
    CALL reallocate(dy,1,3,1,na)
    CALL reallocate(dc,1,nco(llmax),1,3)
    CALL reallocate(harmonics%slm_int,1,max_s_harm)

    slm => harmonics%slm
    dslm_dxyz => harmonics%dslm_dxyz
    slm_int => harmonics%slm_int
    slm_int = 0.0_dp

    DO iso = 1,max_s_harm

       l = indso(1,iso)
       m = indso(2,iso)

       CALL y_lm(lebedev_grid(ll)%r,y,l,m)
       slm(1:na,iso) = y(1:na)

       DO i = 1,na
          slm_int(iso) = slm_int(iso) + slm(i,iso)*wa(i)
       END DO ! i

    END DO  ! iso

    DO i = 1,3
       harmonics%a(i,1:na) = lebedev_grid(ll)%r(i,1:na)
    END DO

    DO ia = 1,na
       DO l = 0,indso(1,max_s_harm)
          DO ic = 1,nco(l)
             lx = indco(1,ic+ncoset(l-1))
             ly = indco(2,ic+ncoset(l-1))
             lz = indco(3,ic+ncoset(l-1))

             IF (lx == 0) THEN
                rx = 1.0_dp
                drx = 0.0_dp
             ELSE IF (lx == 1) THEN
                rx = lebedev_grid(ll)%r(1,ia)
                drx = 1.0_dp
             ELSE
                rx = lebedev_grid(ll)%r(1,ia)**lx
                drx = REAL(lx,dp)*lebedev_grid(ll)%r(1,ia)**(lx - 1)
             END IF
             IF (ly == 0) THEN
                ry = 1.0_dp
                dry = 0.0_dp
             ELSE IF (ly == 1) THEN
                ry = lebedev_grid(ll)%r(2,ia)
                dry = 1.0_dp
             ELSE
                ry = lebedev_grid(ll)%r(2,ia)**ly
                dry = REAL(ly,dp)*lebedev_grid(ll)%r(2,ia)**(ly - 1)
             END IF
             IF (lz == 0) THEN
                rz = 1.0_dp
                drz = 0.0_dp
             ELSE IF (lz == 1) THEN
                rz = lebedev_grid(ll)%r(3,ia)
                drz = 1.0_dp
             ELSE
                rz = lebedev_grid(ll)%r(3,ia)**lz
                drz = REAL(lz,dp)*lebedev_grid(ll)%r(3,ia)**(lz - 1)
             END IF
             dc(ic,1) = drx*ry*rz
             dc(ic,2) = rx*dry*rz
             dc(ic,3) = rx*ry*drz
          END DO
          n = nsoset(l-1)
          DO is = 1,nso(l)
             iso = is + n
             dslm_dxyz(1:3,ia,iso) = 0.0_dp
             DO ic=1,nco(l)
                dslm_dxyz(1,ia,iso) = dslm_dxyz(1,ia,iso) +&
                     orbtramat(l)%slm(is,ic)*dc(ic,1)
                dslm_dxyz(2,ia,iso) = dslm_dxyz(2,ia,iso) +&
                     orbtramat(l)%slm(is,ic)*dc(ic,2)
                dslm_dxyz(3,ia,iso) = dslm_dxyz(3,ia,iso) +&
                     orbtramat(l)%slm(is,ic)*dc(ic,3)
             END DO
          END DO
       END DO  ! l
    END DO

    ! Expansion coefficients of the cartesian derivatives
    ! of the product of two harmonics :
    ! d(Y(l1m1) * Y(l2m2))/dx ; d(Y(l1m1) * Y(l2m2))/dy ; d(Y(l1m1) * Y(l2m2))/dz
    harmonics%my_CG_dxyz = 0.0_dp
    DO iso1 = 1,maxs
       DO iso2 = 1,maxs
          int1 = 0.0_dp
          int2 = 0.0_dp
          int3 = 0.0_dp
          DO ia = 1,na
             int1 = int1 + wa(ia)* &
                  (dslm_dxyz(1,ia,iso1)*slm(ia,iso2)+slm(ia,iso1)*dslm_dxyz(1,ia,iso2))
             int2 = int2 + wa(ia)* &
                  (dslm_dxyz(2,ia,iso1)*slm(ia,iso2)+slm(ia,iso1)*dslm_dxyz(2,ia,iso2))
             int3 = int3 + wa(ia)* &
                  (dslm_dxyz(3,ia,iso1)*slm(ia,iso2)+slm(ia,iso1)*dslm_dxyz(3,ia,iso2))
          END DO
          DO iso = 1,max_s_harm
             rx = 0.0_dp
             drx = 0.0_dp
             ry = 0.0_dp
             dry = 0.0_dp
             rz = 0.0_dp
             drz = 0.0_dp
             DO ia = 1,na
                rx = rx + wa(ia)*slm(ia,iso)*&
                     (dslm_dxyz(1,ia,iso1)*slm(ia,iso2)+slm(ia,iso1)*dslm_dxyz(1,ia,iso2))
                drx = drx + wa(ia)*slm(ia,iso)
                ry = ry + wa(ia)*slm(ia,iso)*&
                     (dslm_dxyz(2,ia,iso1)*slm(ia,iso2)+slm(ia,iso1)*dslm_dxyz(2,ia,iso2))
                dry = dry + wa(ia)*slm(ia,iso)
                rz = rz + wa(ia)*slm(ia,iso)*&
                     (dslm_dxyz(3,ia,iso1)*slm(ia,iso2)+slm(ia,iso1)*dslm_dxyz(3,ia,iso2))
                drz = drz + wa(ia)*slm(ia,iso)
             END DO

             harmonics%my_CG_dxyz(1,iso1,iso2,iso)=rx

             harmonics%my_CG_dxyz(2,iso1,iso2,iso)=ry

             harmonics%my_CG_dxyz(3,iso1,iso2,iso)=rz

          END DO
       END DO
    END DO

    ! Expansion coefficients of the cartesian of the combinations
    ! Y(l1m1) * d(Y(l2m2))/dx -  d(Y(l1m1))/dx * Y(l2m2)
    ! Y(l1m1) * d(Y(l2m2))/dy -  d(Y(l1m1))/dy * Y(l2m2)
    ! Y(l1m1) * d(Y(l2m2))/dz -  d(Y(l1m1))/dz * Y(l2m2)
    harmonics%my_CG_dxyz_asym = 0.0_dp
    DO iso1 = 1,maxs
       DO iso2 = 1,maxs
          DO iso = 1,max_s_harm
             rx = 0.0_dp
             drx = 0.0_dp
             ry = 0.0_dp
             dry = 0.0_dp
             rz = 0.0_dp
             drz = 0.0_dp
             DO ia = 1,na
                rx = rx + wa(ia)*slm(ia,iso)*&
                     (-dslm_dxyz(1,ia,iso1)*slm(ia,iso2)+&
                     slm(ia,iso1)*dslm_dxyz(1,ia,iso2))
                drx = drx + wa(ia)*slm(ia,iso)
                ry = ry + wa(ia)*slm(ia,iso)*&
                     (-dslm_dxyz(2,ia,iso1)*slm(ia,iso2)+&
                     slm(ia,iso1)*dslm_dxyz(2,ia,iso2))
                dry = dry + wa(ia)*slm(ia,iso)
                rz = rz + wa(ia)*slm(ia,iso)*&
                     (-dslm_dxyz(3,ia,iso1)*slm(ia,iso2)+&
                     slm(ia,iso1)*dslm_dxyz(3,ia,iso2))
                drz = drz + wa(ia)*slm(ia,iso)
             END DO

             harmonics%my_CG_dxyz_asym(1,iso1,iso2,iso)=rx

             harmonics%my_CG_dxyz_asym(2,iso1,iso2,iso)=ry

             harmonics%my_CG_dxyz_asym(3,iso1,iso2,iso)=rz

          END DO  ! iso
       END DO  ! iso2
    END DO  ! iso1

    ! Calculate the derivatives of the harmonics with respect of the 2 angles
    ! the first angle is acos(lebedev_grid(ll)%r(3))
    ! the second angle is atan(lebedev_grid(ll)%r(2)/lebedev_grid(ll)%r(1))
    DO iso = 1,max_s_harm

       l = indso(1,iso)
       m = indso(2,iso)
       l_r = REAL(l,dp)
       m_r = REAL(m,dp)

       !  Derivative with respect to the polar angle
       IF(m==0) THEN
          harmonics%dslm(2,1:na,iso) = 0.0_dp
       ELSE
          iso_mm = soset(l,-m)
          harmonics%dslm(2,1:na,iso) = - m_r * harmonics%slm(1:na,iso_mm)
       END IF

       ! Derivative with respect to the azimuthal angle
       IF(m > 0) THEN

          DO ia = 1,na

             harmonics%dslm(1,ia,iso) = m_r * cotan_azi(ia)*harmonics%slm(ia,iso)

             IF(m /= l ) THEN
                iso_mp1 = soset(l,m+1)

                cosmp1 = COS( (m_r+1) *pol(ia) )
                IF (ABS(cosmp1) > EPSILON(1.0_dp ) ) THEN
                   harmonics%dslm(1,ia,iso) = harmonics%dslm(1,ia,iso) +&
                        COS( m_r * pol(ia) )*harmonics%slm(ia,iso_mp1)/ &
                        cosmp1 * SQRT( (l_r-m_r) * (l_r+m_r+1._dp) )
                END IF

             END IF
          END DO

       ELSEIF(m == 0) THEN
          IF(l == 0) THEN
             harmonics%dslm(1,1:na,iso) = 0.0_dp
          ELSE
             iso_mp1 = soset(l,1)
             iso_mm1 = soset(l,-1)
             DO ia = 1,na

                IF(ABS(cos_pol(ia)) > EPSILON(1.0_dp ) ) THEN
                   harmonics%dslm(1,ia,iso) = SQRT(0.5_dp * l_r * (l_r+1._dp) )*&
                        harmonics%slm(ia,iso_mp1)/cos_pol(ia)
                END IF
             END DO
          END IF
       ELSE

          abs_m = ABS(m_r)
          DO ia = 1,na
             harmonics%dslm(1,ia,iso) = abs_m * cotan_azi(ia)*&
                  harmonics%slm(ia,iso)
             IF(m /= -l) THEN
                iso_mm1 = soset(l,m-1)

                sinmp1 =  SIN( (abs_m+1._dp)*pol(ia) )
                IF(ABS(sinmp1) > EPSILON(1.0_dp) ) THEN
                   harmonics%dslm(1,ia,iso) = harmonics%dslm(1,ia,iso) +&
                        SIN( abs_m*pol(ia) ) * harmonics%slm(ia,iso_mm1)/&
                        sinmp1 * SQRT( (l_r-abs_m) * (l_r+abs_m+1._dp) )
                END IF
             END IF
          END DO

       END IF

    END DO

    ! expansion coefficients of product of theta derivatives (dslm(1...)) in
    ! spherical harmonics (used for tau functionals)
    CALL reallocate(harmonics%my_dCG,1,maxs,1,maxs,1,max_s_harm)
    DO is1=1,maxs
      l1 = indso(1,is1)
      m1 = indso(2,is1)
      DO is2=1,maxs
        l2 = indso(1,is2)
        m2 = indso(2,is2)
        DO iso=1,max_s_harm
          l = indso(1,iso)
          m = indso(2,iso)
          CALL y_lm(lebedev_grid(ll)%r,y,l,m)
          y(1:na) = y(1:na)*lebedev_grid(ll)%w(1:na)*harmonics%dslm(1,1:na,is1)*harmonics%dslm(1,1:na,is2)
          harmonics%my_dCG(is1,is2,iso) = SUM(y(1:na))
        END DO
      END DO
    END DO

    DEALLOCATE(y,dy,dc)

    CALL timestop(handle)

  END SUBROUTINE  create_harmonics_atom

! *****************************************************************************

  SUBROUTINE get_maxl_CG(harmonics,orb_basis,my_CG,llmax,max_s_harm,error)

    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    INTEGER, INTENT(IN)                      :: llmax, max_s_harm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_maxl_CG', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: damax_iso_not0, &
                                                dmax_iso_not0, handle, is1, &
                                                is2, itmp, max_iso_not0, nset
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL timeset(routineN,handle)

    CPPrecondition(ASSOCIATED(harmonics),cp_failure_level,routineP,error,failure)

    CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,nset=nset)

    !   *** Assign indexes for the non null CG coefficients ***
    max_iso_not0 = 0
    dmax_iso_not0 = 0
    damax_iso_not0 =  0
    DO is1 = 1,nset
       DO is2 = 1,nset
          CALL get_none0_cg_list(harmonics%my_CG,&
                                 lmin(is1),lmax(is1),lmin(is2),lmax(is2),&
                                 max_s_harm,llmax,max_iso_not0=itmp,error=error)
          max_iso_not0 = MAX(max_iso_not0,itmp)
          CALL get_none0_cg_list(harmonics%my_dCG,&
                                 lmin(is1),lmax(is1),lmin(is2),lmax(is2),&
                                 max_s_harm,llmax,max_iso_not0=itmp,error=error)
          dmax_iso_not0 = MAX(dmax_iso_not0,itmp)
          CALL get_none0_cg_list(harmonics%my_CG_dxyz_asym,&
                                 lmin(is1),lmax(is1),lmin(is2),lmax(is2),&
                                 max_s_harm,llmax,max_iso_not0=itmp,error=error)
          damax_iso_not0 = MAX(damax_iso_not0,itmp)
       END DO ! is2
    END DO ! is1
    harmonics%max_iso_not0 = max_iso_not0
    harmonics%dmax_iso_not0 = dmax_iso_not0
    harmonics%damax_iso_not0 = damax_iso_not0

    CALL timestop(handle)

  END SUBROUTINE  get_maxl_CG

  SUBROUTINE get_none0_cg_list4(cgc,lmin1,lmax1,lmin2,lmax2,max_s_harm,llmax,&
                                list,n_list,max_iso_not0,error)

    REAL(dp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: cgc
    INTEGER, INTENT(IN)                      :: lmin1, lmax1, lmin2, lmax2, &
                                                max_s_harm, llmax
    INTEGER, DIMENSION(:, :, :), &
      INTENT(OUT), OPTIONAL                  :: list
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: n_list
    INTEGER, INTENT(OUT)                     :: max_iso_not0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_none0_cg_list4', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iso, iso1, iso2, l1, l2, nlist
    LOGICAL                                  :: failure

    CPPrecondition(nsoset(lmax1).LE.SIZE(cgc,2),cp_failure_level,routineP,error,failure)
    CPPrecondition(nsoset(lmax2).LE.SIZE(cgc,3),cp_failure_level,routineP,error,failure)
    CPPrecondition(max_s_harm.LE.SIZE(cgc,4),cp_failure_level,routineP,error,failure)
    IF(PRESENT(n_list).AND.PRESENT(list)) THEN
       CPPrecondition(max_s_harm.LE.SIZE(list,3),cp_failure_level,routineP,error,failure)
    ENDIF
    max_iso_not0 = 0
    IF(PRESENT(n_list).AND.PRESENT(list)) n_list = 0
    DO iso = 1,max_s_harm
       nlist = 0
       DO l1 = lmin1,lmax1
          DO iso1 = nsoset(l1-1)+1,nsoset(l1)
             DO l2 = lmin2,lmax2
                IF(l1+l2 > llmax) CYCLE
                DO iso2 = nsoset(l2-1)+1,nsoset(l2)
                   IF(ABS(cgc(1,iso1,iso2,iso))+&
                      ABS(cgc(2,iso1,iso2,iso))+&
                      ABS(cgc(3,iso1,iso2,iso)) > 1.E-8_dp ) THEN
                      nlist = nlist + 1
                      IF(PRESENT(n_list).AND.PRESENT(list)) THEN
                         list(1,nlist,iso) = iso1
                         list(2,nlist,iso) = iso2
                      ENDIF
                      max_iso_not0 = MAX(max_iso_not0,iso)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       IF(PRESENT(n_list).AND.PRESENT(list)) n_list(iso) = nlist
    ENDDO
  END SUBROUTINE get_none0_cg_list4

  SUBROUTINE get_none0_cg_list3(cgc,lmin1,lmax1,lmin2,lmax2,max_s_harm,llmax,&
                                list,n_list,max_iso_not0,error)

    REAL(dp), DIMENSION(:, :, :), INTENT(IN) :: cgc
    INTEGER, INTENT(IN)                      :: lmin1, lmax1, lmin2, lmax2, &
                                                max_s_harm, llmax
    INTEGER, DIMENSION(:, :, :), &
      INTENT(OUT), OPTIONAL                  :: list
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: n_list
    INTEGER, INTENT(OUT)                     :: max_iso_not0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_none0_cg_list3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iso, iso1, iso2, l1, l2, nlist
    LOGICAL                                  :: failure

    CPPrecondition(nsoset(lmax1).LE.SIZE(cgc,1),cp_failure_level,routineP,error,failure)
    CPPrecondition(nsoset(lmax2).LE.SIZE(cgc,2),cp_failure_level,routineP,error,failure)
    CPPrecondition(max_s_harm.LE.SIZE(cgc,3),cp_failure_level,routineP,error,failure)
    IF(PRESENT(n_list).AND.PRESENT(list)) THEN
       CPPrecondition(max_s_harm.LE.SIZE(list,3),cp_failure_level,routineP,error,failure)
    ENDIF
    max_iso_not0 = 0
    IF(PRESENT(n_list).AND.PRESENT(list)) n_list = 0
    DO iso = 1,max_s_harm
       nlist = 0
       DO l1 = lmin1,lmax1
          DO iso1 = nsoset(l1-1)+1,nsoset(l1)
             DO l2 = lmin2, lmax2
                IF(l1+l2 > llmax) CYCLE
                DO iso2 = nsoset(l2-1)+1,nsoset(l2)
                   IF(ABS(cgc(iso1,iso2,iso)) > 1.E-8_dp ) THEN
                      nlist = nlist + 1
                      IF(PRESENT(n_list).AND.PRESENT(list)) THEN
                         list(1,nlist,iso) = iso1
                         list(2,nlist,iso) = iso2
                      ENDIF
                      max_iso_not0 = MAX(max_iso_not0,iso)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       IF(PRESENT(n_list).AND.PRESENT(list)) n_list(iso) = nlist
    ENDDO
  END SUBROUTINE get_none0_cg_list3

END MODULE qs_harmonics_atom
