!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE qs_rho0_types

  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE erf_fn,                          ONLY: erf
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE mathconstants,                   ONLY: fourpi,&
                                             pi,&
                                             rootpi
  USE memory_utilities,                ONLY: reallocate
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_rho_atom_types,               ONLY: rho_atom_coeff
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE whittaker,                       ONLY: whittaker_c0a,&
                                             whittaker_ci
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho0_types'

! *** Define multipole type ***

! *****************************************************************************
  TYPE  mpole_rho_atom
    REAL(dp), DIMENSION(:), POINTER             ::  Qlm_h, &
                                                    Qlm_s, &
                                                    Qlm_tot,&
                                                    Qlm_car
    REAL(dp)                                    ::  Qlm_z
    REAL(dp), DIMENSION(2)                      ::  Q0
  END TYPE mpole_rho_atom

! *****************************************************************************
  TYPE mpole_gau_overlap
    REAL(dp), DIMENSION(:,:,:), POINTER         :: Qlm_gg
    REAL(dp), DIMENSION(:,:), POINTER           :: g0_h, vg0_h
    REAL(dp)                                    :: rpgf0_h, rpgf0_s
  END TYPE mpole_gau_overlap

! *****************************************************************************
  TYPE rho0_mpole_type
    TYPE(mpole_rho_atom), DIMENSION(:), POINTER  :: mp_rho
    TYPE(mpole_gau_overlap), DIMENSION(:), &
                                      POINTER   :: mp_gau
    REAL(dp)                                    :: zet0_h,&
                                                   total_rho0_h
    REAL(dp)                                    :: max_rpgf0_s
    REAL(dp), DIMENSION(:), POINTER             :: norm_g0l_h
    INTEGER, DIMENSION(:),  POINTER             :: lmax0_kind
    INTEGER                                     :: lmax_0,igrid_zet0_s
    TYPE(pw_p_type), POINTER                    :: rho0_s_rs,&
                                                   rho0_s_gs
  END TYPE rho0_mpole_type

! *****************************************************************************
  TYPE rho0_atom_type
    TYPE(rho_atom_coeff), POINTER               :: rho0_rad_h, &
                                                   vrho0_rad_h
  END TYPE rho0_atom_type

! Public Types

  PUBLIC :: mpole_rho_atom, mpole_gau_overlap, &
            rho0_atom_type, rho0_mpole_type

! Public Subroutine

  PUBLIC :: allocate_multipoles, allocate_rho0_mpole, &
            allocate_rho0_atom, allocate_rho0_atom_rad, &
            deallocate_rho0_atom, deallocate_rho0_mpole, &
            calculate_g0, get_rho0_mpole, initialize_mpole_rho, &
            write_rho0_info

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param mp_rho ...
!> \param natom ...
!> \param mp_gau ...
!> \param nkind ...
! *****************************************************************************
  SUBROUTINE allocate_multipoles(mp_rho,natom,mp_gau,nkind)

    TYPE(mpole_rho_atom), DIMENSION(:), &
      POINTER                                :: mp_rho
    INTEGER, INTENT(IN)                      :: natom
    TYPE(mpole_gau_overlap), DIMENSION(:), &
      POINTER                                :: mp_gau
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_multipoles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, ikind, istat

    IF(ASSOCIATED(mp_rho)) THEN
      CALL deallocate_mpole_rho(mp_rho)
    END IF

    ALLOCATE (mp_rho(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                            "mp_rho",natom*int_size)

    DO iat = 1,natom
      NULLIFY(mp_rho(iat)%Qlm_h)
      NULLIFY(mp_rho(iat)%Qlm_s)
      NULLIFY(mp_rho(iat)%Qlm_tot)
      NULLIFY(mp_rho(iat)%Qlm_car)
    END DO

    IF(ASSOCIATED(mp_gau)) THEN
      CALL deallocate_mpole_gau(mp_gau)
    END IF

    ALLOCATE (mp_gau(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                            "mp_gau",nkind*int_size)

    DO ikind = 1,nkind
      NULLIFY(mp_gau(ikind)%Qlm_gg)
      NULLIFY(mp_gau(ikind)%g0_h)
      NULLIFY(mp_gau(ikind)%vg0_h)
      mp_gau(ikind)%rpgf0_h = 0.0_dp
      mp_gau(ikind)%rpgf0_s = 0.0_dp
    END DO

  END SUBROUTINE allocate_multipoles

! *****************************************************************************
!> \brief ...
!> \param rho0_set ...
!> \param natom ...
! *****************************************************************************
  SUBROUTINE allocate_rho0_atom(rho0_set,natom)

    TYPE(rho0_atom_type), DIMENSION(:), &
      POINTER                                :: rho0_set
    INTEGER                                  :: natom

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_rho0_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, istat

    IF(ASSOCIATED(rho0_set)) THEN
      CALL deallocate_rho0_atom(rho0_set)
    END IF

    ALLOCATE (rho0_set(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                            "rho0_atom_set",natom*int_size)

    DO iat = 1,natom
      NULLIFY(rho0_set(iat)%rho0_rad_h)
      NULLIFY(rho0_set(iat)%vrho0_rad_h)
    ENDDO

  END SUBROUTINE allocate_rho0_atom

! *****************************************************************************
!> \brief ...
!> \param rho0_atom ...
!> \param nr ...
!> \param nchannels ...
! *****************************************************************************
  SUBROUTINE allocate_rho0_atom_rad(rho0_atom,nr,nchannels)

    TYPE(rho0_atom_type)                     :: rho0_atom
    INTEGER                                  :: nr, nchannels

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_rho0_atom_rad', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    ALLOCATE(rho0_atom%rho0_rad_h,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                            "rho0_rad_h",int_size)

    NULLIFY(rho0_atom%rho0_rad_h%r_coef)
    ALLOCATE (rho0_atom%rho0_rad_h%r_coef(1:nr,1:nchannels),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                            "rho0_rad_h",int_size)
    rho0_atom%rho0_rad_h%r_coef = 0.0_dp

    ALLOCATE(rho0_atom%vrho0_rad_h,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                            "vrho0_rad_h",int_size)

    NULLIFY(rho0_atom%vrho0_rad_h%r_coef)
    ALLOCATE (rho0_atom%vrho0_rad_h%r_coef(1:nr,1:nchannels),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                            "vrho0_rad_h",int_size)
    rho0_atom%vrho0_rad_h%r_coef = 0.0_dp

  END SUBROUTINE allocate_rho0_atom_rad

! *****************************************************************************
!> \brief ...
!> \param rho0 ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE allocate_rho0_mpole(rho0,error)

    TYPE(rho0_mpole_type), POINTER           :: rho0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_rho0_mpole', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF(ASSOCIATED(rho0)) THEN
      CALL deallocate_rho0_mpole(rho0,error=error)
    END IF

    ALLOCATE (rho0,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                            "rho0_mpole",int_size)

    NULLIFY(rho0%mp_rho)
    NULLIFY(rho0%mp_gau)
    NULLIFY(rho0%norm_g0l_h)
    NULLIFY(rho0%lmax0_kind)
    NULLIFY(rho0%rho0_s_rs)
    NULLIFY(rho0%rho0_s_gs)

  END SUBROUTINE allocate_rho0_mpole

! *****************************************************************************
!> \brief ...
!> \param rho0_mpole ...
!> \param grid_atom ...
!> \param ik ...
! *****************************************************************************
  SUBROUTINE calculate_g0(rho0_mpole,grid_atom,ik)

    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(grid_atom_type), POINTER            :: grid_atom
    INTEGER                                  :: ik

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_g0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ir, istat, l, lmax, nr
    REAL(dp)                                 :: c1, prefactor, root_z_h, z_h
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: erf_z_h, gexp, gh_tmp, int1, &
                                                int2

    nr = grid_atom%nr
    lmax = rho0_mpole%lmax0_kind(ik)
    z_h = rho0_mpole%zet0_h
    root_z_h = SQRT(z_h)

!   Allocate g0
    CALL reallocate(rho0_mpole%mp_gau(ik)%g0_h,1,nr,0,lmax)
    CALL reallocate(rho0_mpole%mp_gau(ik)%vg0_h,1,nr,0,lmax)

    ALLOCATE(gexp(nr),gh_tmp(nr),erf_z_h(nr),int1(nr),int2(nr),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                              "gexp,gh_tmp,erf_z_h,in1,int2",6*nr*dp_size)

    gh_tmp(1:nr) = EXP(-z_h*grid_atom%rad2(1:nr))

    DO ir = 1,nr
      erf_z_h(ir) = erf(grid_atom%rad(ir)*root_z_h)
    END DO

    DO ir =1,nr
      IF(gh_tmp(ir) < 1.0E-30_dp) gh_tmp(ir) = 0.0_dp
    END DO

    gexp(1:nr) = gh_tmp(1:nr)
    rho0_mpole%mp_gau(ik)%g0_h(1:nr,0) = gh_tmp(1:nr)* &
                                         rho0_mpole%norm_g0l_h(0)
    CALL whittaker_c0a(int1,grid_atom%rad,gh_tmp,erf_z_h,z_h,0,0,nr)
    CALL whittaker_ci(int2,grid_atom%rad,gh_tmp,z_h,0,nr)

    prefactor = fourpi*rho0_mpole%norm_g0l_h(0)

    c1 = SQRT(pi*pi*pi/(z_h*z_h*z_h))*rho0_mpole%norm_g0l_h(0)

    DO ir = 1,nr
      rho0_mpole%mp_gau(ik)%vg0_h(ir,0) = c1*erf_z_h(ir)*grid_atom%oorad2l(ir,1)
    END DO

    DO l = 1,lmax
      gh_tmp(1:nr) = gh_tmp(1:nr)*grid_atom%rad(1:nr)
      rho0_mpole%mp_gau(ik)%g0_h(1:nr,l) = gh_tmp(1:nr)* &
                                           rho0_mpole%norm_g0l_h(l)

      prefactor = fourpi/(2.0_dp*l+1.0_dp)*rho0_mpole%norm_g0l_h(l)
      CALL whittaker_c0a(int1,grid_atom%rad,gexp,erf_z_h,z_h,l,l,nr)
      DO ir = 1,nr
        rho0_mpole%mp_gau(ik)%vg0_h(ir,l) = prefactor*(int1(ir)+&
                                      int2(ir)*grid_atom%rad2l(ir,l))
      END DO

    END DO  ! l

    DEALLOCATE (gexp,erf_z_h,gh_tmp,int1,int2,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "gexp,gh_tmp,erf_z_h,int1,int2")
  END SUBROUTINE calculate_g0

! *****************************************************************************
!> \brief ...
!> \param mp_gau ...
! *****************************************************************************
  SUBROUTINE deallocate_mpole_gau(mp_gau)

    TYPE(mpole_gau_overlap), DIMENSION(:), &
      POINTER                                :: mp_gau

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_mpole_gau', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, nkind

    istat = 0
    nkind = SIZE(mp_gau)

    DO ikind = 1,nkind

      IF(ASSOCIATED(mp_gau(ikind)%Qlm_gg)) THEN
        DEALLOCATE(mp_gau(ikind)%Qlm_gg,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                      "Qlm_gg")
      END IF

      DEALLOCATE(mp_gau(ikind)%g0_h, STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                      "mp_gau(ikind)%g0_h")

      DEALLOCATE(mp_gau(ikind)%vg0_h, STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                      "mp_gau(ikind)%vg0_h")
    END DO

    DEALLOCATE(mp_gau, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                      "mp_gau")

  END SUBROUTINE deallocate_mpole_gau

! *****************************************************************************
!> \brief ...
!> \param mp_rho ...
! *****************************************************************************
  SUBROUTINE deallocate_mpole_rho(mp_rho)

    TYPE(mpole_rho_atom), DIMENSION(:), &
      POINTER                                :: mp_rho

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_mpole_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, istat, natom

    natom = SIZE(mp_rho)

    DO iat = 1,natom
      DEALLOCATE(mp_rho(iat)%Qlm_h,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                      "Qlm_h")
      DEALLOCATE(mp_rho(iat)%Qlm_s,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                      "Qlm_s")
      DEALLOCATE(mp_rho(iat)%Qlm_tot,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                      "Qlm_tot")
      DEALLOCATE(mp_rho(iat)%Qlm_car,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                      "Qlm_car")
    END DO

    DEALLOCATE(mp_rho,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                      "mp_rho")

  END SUBROUTINE deallocate_mpole_rho

! *****************************************************************************
!> \brief ...
!> \param rho0_atom_set ...
! *****************************************************************************
  SUBROUTINE deallocate_rho0_atom(rho0_atom_set)

    TYPE(rho0_atom_type), DIMENSION(:), &
      POINTER                                :: rho0_atom_set

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_rho0_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, istat, natom

    IF (ASSOCIATED(rho0_atom_set)) THEN

      natom = SIZE(rho0_atom_set)

      DO iat = 1,natom
        IF(ASSOCIATED(rho0_atom_set(iat)%rho0_rad_h)) THEN
          DEALLOCATE (rho0_atom_set(iat)%rho0_rad_h%r_coef)
          DEALLOCATE (rho0_atom_set(iat)%rho0_rad_h)
        ENDIF
        IF(ASSOCIATED(rho0_atom_set(iat)%vrho0_rad_h)) THEN
          DEALLOCATE (rho0_atom_set(iat)%vrho0_rad_h%r_coef)
          DEALLOCATE (rho0_atom_set(iat)%vrho0_rad_h)
        ENDIF
      ENDDO

      DEALLOCATE (rho0_atom_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "rho0_atom_set")
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer rho0_atom_set is not associated and "//&
                        "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_rho0_atom
! *****************************************************************************
!> \brief ...
!> \param rho0 ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE deallocate_rho0_mpole(rho0,error)

    TYPE(rho0_mpole_type), POINTER           :: rho0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_rho0_mpole', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF (ASSOCIATED(rho0)) THEN

      IF (ASSOCIATED(rho0%mp_gau)) CALL deallocate_mpole_gau(rho0%mp_gau)

      IF (ASSOCIATED(rho0%mp_rho)) CALL deallocate_mpole_rho(rho0%mp_rho)

      IF (ASSOCIATED(rho0%lmax0_kind)) THEN
         DEALLOCATE(rho0%lmax0_kind,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                         "rho0%lmax0_kind")
      END IF

      IF (ASSOCIATED(rho0%norm_g0l_h)) THEN
         DEALLOCATE(rho0%norm_g0l_h, STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                         "rho0%norm_g0l_h")
      END IF

      IF (ASSOCIATED(rho0%rho0_s_rs)) THEN
          CALL pw_release(rho0%rho0_s_rs%pw,error=error)
          DEALLOCATE(rho0%rho0_s_rs, STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                         "rho0%rho0_s_rs")
      ENDIF

      IF (ASSOCIATED(rho0%rho0_s_gs)) THEN
          CALL pw_release(rho0%rho0_s_gs%pw,error=error)
          DEALLOCATE(rho0%rho0_s_gs, STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                         "rho0%rho0_s_gs")

      ENDIF
      DEALLOCATE (rho0,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "rho0_mpole")
    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer rho0 is not associated and "//&
                        "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_rho0_mpole

! *****************************************************************************
!> \brief ...
!> \param rho0_mpole ...
!> \param g0_h ...
!> \param vg0_h ...
!> \param iat ...
!> \param ikind ...
!> \param lmax_0 ...
!> \param l0_ikind ...
!> \param mp_gau_ikind ...
!> \param mp_rho ...
!> \param norm_g0l_h ...
!> \param Qlm_gg ...
!> \param Qlm_car ...
!> \param Qlm_tot ...
!> \param zet0_h ...
!> \param igrid_zet0_s ...
!> \param rpgf0_h ...
!> \param rpgf0_s ...
!> \param max_rpgf0_s ...
!> \param rho0_s_rs ...
!> \param rho0_s_gs ...
! *****************************************************************************
  SUBROUTINE get_rho0_mpole(rho0_mpole, g0_h, vg0_h, iat, ikind, lmax_0, l0_ikind,&
                            mp_gau_ikind, mp_rho, norm_g0l_h, &
                            Qlm_gg, Qlm_car,  Qlm_tot,  &
                            zet0_h, igrid_zet0_s, rpgf0_h,rpgf0_s, &
                            max_rpgf0_s, rho0_s_rs, rho0_s_gs)

    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: g0_h, vg0_h
    INTEGER, INTENT(IN), OPTIONAL            :: iat, ikind
    INTEGER, INTENT(OUT), OPTIONAL           :: lmax_0, l0_ikind
    TYPE(mpole_gau_overlap), OPTIONAL, &
      POINTER                                :: mp_gau_ikind
    TYPE(mpole_rho_atom), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mp_rho
    REAL(dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: norm_g0l_h
    REAL(dp), DIMENSION(:, :, :), OPTIONAL, &
      POINTER                                :: Qlm_gg
    REAL(dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: Qlm_car, Qlm_tot
    REAL(dp), INTENT(OUT), OPTIONAL          :: zet0_h
    INTEGER, INTENT(OUT), OPTIONAL           :: igrid_zet0_s
    REAL(dp), INTENT(OUT), OPTIONAL          :: rpgf0_h, rpgf0_s, max_rpgf0_s
    TYPE(pw_p_type), OPTIONAL, POINTER       :: rho0_s_rs, rho0_s_gs

    CHARACTER(len=*), PARAMETER :: routineN = 'get_rho0_mpole', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(rho0_mpole)) THEN

      IF(PRESENT(lmax_0)) lmax_0 = rho0_mpole%lmax_0
      IF(PRESENT(mp_rho)) mp_rho => rho0_mpole%mp_rho
      IF(PRESENT(norm_g0l_h)) norm_g0l_h => rho0_mpole%norm_g0l_h
      IF(PRESENT(zet0_h)) zet0_h = rho0_mpole%zet0_h
      IF(PRESENT(igrid_zet0_s)) igrid_zet0_s = rho0_mpole%igrid_zet0_s
      IF(PRESENT(max_rpgf0_s)) max_rpgf0_s = rho0_mpole%max_rpgf0_s
      IF(PRESENT(rho0_s_rs))    rho0_s_rs => rho0_mpole%rho0_s_rs
      IF(PRESENT(rho0_s_gs))    rho0_s_gs => rho0_mpole%rho0_s_gs

      IF(PRESENT(ikind)) THEN
        IF(PRESENT(l0_ikind))     l0_ikind = rho0_mpole%lmax0_kind(ikind)
        IF(PRESENT(mp_gau_ikind)) mp_gau_ikind => rho0_mpole%mp_gau(ikind)
        IF(PRESENT(g0_h))         g0_h => rho0_mpole%mp_gau(ikind)%g0_h
        IF(PRESENT(vg0_h))        vg0_h => rho0_mpole%mp_gau(ikind)%vg0_h
        IF(PRESENT(Qlm_gg))       Qlm_gg => rho0_mpole%mp_gau(ikind)%Qlm_gg
        IF(PRESENT(rpgf0_h))      rpgf0_h = rho0_mpole%mp_gau(ikind)%rpgf0_h
        IF(PRESENT(rpgf0_s))      rpgf0_s = rho0_mpole%mp_gau(ikind)%rpgf0_s
      END IF
      IF(PRESENT(iat)) THEN
        IF(PRESENT(Qlm_car)) Qlm_car => rho0_mpole%mp_rho(iat)%Qlm_car
        IF(PRESENT(Qlm_tot)) Qlm_tot => rho0_mpole%mp_rho(iat)%Qlm_tot
      END IF

    ELSE
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "The pointer rho0_mpole is not associated")
    END IF

  END SUBROUTINE get_rho0_mpole

! *****************************************************************************
!> \brief ...
!> \param mp_rho ...
!> \param nchan_s ...
!> \param nchan_c ...
!> \param zeff ...
!> \param tddft ...
! *****************************************************************************
  SUBROUTINE initialize_mpole_rho(mp_rho,nchan_s,nchan_c,zeff,tddft)

    TYPE(mpole_rho_atom)                     :: mp_rho
    INTEGER, INTENT(IN)                      :: nchan_s, nchan_c
    REAL(KIND=dp), INTENT(IN)                :: zeff
    LOGICAL, OPTIONAL                        :: tddft

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_mpole_rho', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: my_tddft

    my_tddft = .FALSE.
    IF (PRESENT(tddft)) my_tddft = tddft

   CALL reallocate(mp_rho%Qlm_h,1,nchan_s)
   CALL reallocate(mp_rho%Qlm_s,1,nchan_s)
   CALL reallocate(mp_rho%Qlm_tot,1,nchan_s)
   CALL reallocate(mp_rho%Qlm_car,1,nchan_c)

   mp_rho%Qlm_h = 0.0_dp
   mp_rho%Qlm_s = 0.0_dp
   mp_rho%Qlm_tot = 0.0_dp
   mp_rho%Qlm_car = 0.0_dp
   IF (.NOT.my_tddft) THEN
      mp_rho%Qlm_z = -2.0_dp*rootpi*Zeff
   ELSE
      mp_rho%Qlm_z = 0.0_dp
   END IF
   mp_rho%Q0 = 0.0_dp

  END SUBROUTINE initialize_mpole_rho

! *****************************************************************************
!> \brief ...
!> \param rho0_mpole ...
!> \param unit_str ...
!> \param output_unit ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_rho0_info(rho0_mpole,unit_str,output_unit,error)

    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    CHARACTER(LEN=*), INTENT(IN)             :: unit_str
    INTEGER, INTENT(in)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_rho0_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, l, nkind
    REAL(dp)                                 :: conv

    IF (ASSOCIATED(rho0_mpole)) THEN
      conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)

      WRITE (UNIT=output_unit,FMT="(/,T5,A,/)") &
               "*** Compensation density charges data set ***"
      WRITE (UNIT=output_unit,FMT="(T2,A,T35,f16.10)")&
            "- Rho0 exponent :",rho0_mpole%zet0_h
      WRITE (UNIT=output_unit,FMT="(T2,A,T35,I5)")&
            "- Global max l :",rho0_mpole%lmax_0

      WRITE (UNIT=output_unit,FMT="(T2,A)")&
            "- Normalization constants for g0"
      DO l = 0,rho0_mpole%lmax_0
        WRITE (UNIT=output_unit,FMT="(T20,A,T31,I2,T38,A,f15.5)")&
            "ang. mom.= ", l, " hard= ", rho0_mpole%norm_g0l_h(l)
      END DO

      nkind = SIZE(rho0_mpole%lmax0_kind,1)
      DO ikind = 1, nkind
        WRITE (UNIT=output_unit,FMT="(/,T2,A,T55,I2)")&
            "- rho0 max L and radii in "//TRIM(unit_str)//&
            " for the atom kind :", ikind

        WRITE (UNIT=output_unit,FMT="(T2,T20,A,T55,I5)")&
            "=> l max  :",rho0_mpole%lmax0_kind(ikind)

        WRITE (UNIT=output_unit,FMT="(T2,T20,A,T55,f20.10)")&
            "=> max radius of g0: ",&
            rho0_mpole%mp_gau(ikind)%rpgf0_h*conv
      END DO  ! ikind

    ELSE
      WRITE(UNIT=output_unit,FMT="(/,T5,A,/)") &
           ' WARNING: I cannot print rho0, it is not associated'
    END IF

  END SUBROUTINE write_rho0_info
END MODULE qs_rho0_types
