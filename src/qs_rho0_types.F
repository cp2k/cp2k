!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_rho0_types
!!
!!   NAME
!!     qs_multipole_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   SOURCE
!******************************************************************************

MODULE qs_rho0_types

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE mathconstants,                   ONLY: rootpi
  USE memory_utilities,                ONLY: reallocate
  USE qs_rho_atom_types,               ONLY: rho_atom_coeff
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_multipole_types"

! *** Define multipole type ***

  TYPE  mpole_rho_atom
    REAL(dp), DIMENSION(:), POINTER ::  Qlm_h,Qlm_s,Qlm_tot
    REAL(dp)                        ::  Qlm_z
  END TYPE mpole_rho_atom

  TYPE mpole_gau_overlap
    REAL(dp), DIMENSION(:,:,:), POINTER :: Qlm_gg 
  END TYPE mpole_gau_overlap

  TYPE rho0_mpole_type
    TYPE(mpole_rho_atom), DIMENSION(:), POINTER  :: mp_rho
    TYPE(mpole_gau_overlap), DIMENSION(:), &
                                      POINTER   :: mp_gau
    REAL(dp)                                    :: zet0_h,zet0_s
    REAL(dp), DIMENSION(:), POINTER             :: norm_g0l_h, &
                                                   norm_g0l_s
    INTEGER                                     :: lmax_0
    INTEGER, DIMENSION(:),  POINTER             :: lmax0_kind
  END TYPE rho0_mpole_type

  TYPE rho0_atom_type
    TYPE(rho_atom_coeff), POINTER               :: rho0_rad_h, &
                                                   rho0_rad_s
  END TYPE rho0_atom_type

! Public Types

  PUBLIC :: mpole_rho_atom, mpole_gau_overlap, &
            rho0_atom_type, rho0_mpole_type

! Public Subroutine
 
  PUBLIC :: allocate_multipoles, allocate_rho0_mpole, &
            allocate_rho0_atom, deallocate_rho0_atom, &
            allocate_rho0_atom_rad, &
            initialize_mpole_rho, & 
            deallocate_rho0_mpole, write_rho0_info

!******************************************************************************
 
CONTAINS
 
!******************************************************************************
  SUBROUTINE allocate_multipoles(mp_rho,natom,mp_gau,nkind)

    TYPE(mpole_rho_atom), DIMENSION(:),    POINTER :: mp_rho
    TYPE(mpole_gau_overlap), DIMENSION(:), POINTER :: mp_gau
    INTEGER, INTENT(IN)                            :: natom,nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_multipoles" 

    INTEGER                                        :: iat, ikind, istat

!   ---------------------------------------------------------------------------

    IF(ASSOCIATED(mp_rho)) THEN
      CALL deallocate_mpole_rho(mp_rho)
    END IF

    ALLOCATE (mp_rho(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                            "allocate_multipoles",natom*int_size)

    DO iat = 1,natom
      NULLIFY(mp_rho(iat)%Qlm_h)
      NULLIFY(mp_rho(iat)%Qlm_s)
      NULLIFY(mp_rho(iat)%Qlm_tot)
    END DO

    IF(ASSOCIATED(mp_rho)) THEN
      CALL deallocate_mpole_gau(mp_gau)
    END IF

    ALLOCATE (mp_gau(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                            "allocate_multipoles",nkind*int_size)

    DO ikind = 1,nkind
      NULLIFY(mp_gau(ikind)%Qlm_gg)
    END DO

  END SUBROUTINE allocate_multipoles

! ******************************************************************************
  SUBROUTINE allocate_rho0_atom(rho0_set,natom)

    TYPE(rho0_atom_type), DIMENSION(:), &
                                   POINTER :: rho0_set
    INTEGER                                :: natom
    
    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_rho0_atom"

    INTEGER                                :: iat, istat

!   ---------------------------------------------------------------------------

    IF(ASSOCIATED(rho0_set)) THEN
      CALL deallocate_rho0_atom(rho0_set)
    END IF 

    ALLOCATE (rho0_set(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                            "allocate_rho0_atom",natom*int_size)

    DO iat = 1,natom 
      NULLIFY(rho0_set(iat)%rho0_rad_h) 
      NULLIFY(rho0_set(iat)%rho0_rad_s)
!      NULLIFY(rho0_set(iat)%drho0_rad_h)
!      NULLIFY(rho0_set(iat)%drho0_rad_s)
    ENDDO

  END SUBROUTINE allocate_rho0_atom

! ******************************************************************************

  SUBROUTINE allocate_rho0_atom_rad(rho0_atom,nr,nchannels)

    TYPE(rho0_atom_type)                :: rho0_atom
    INTEGER                             :: nr,nchannels
    
    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_rho0_atom_rad"

    INTEGER                             :: iat

!   ---------------------------------------------------------------------------

    ALLOCATE (rho0_atom%rho0_rad_h%r_coef(1:nr,1:nchannels))
    rho0_atom%rho0_rad_h%r_coef = 0.0_dp
    ALLOCATE (rho0_atom%rho0_rad_s%r_coef(1:nr,1:nchannels))
    rho0_atom%rho0_rad_s%r_coef = 0.0_dp

  END SUBROUTINE allocate_rho0_atom_rad
! ******************************************************************************

  SUBROUTINE allocate_rho0_mpole(rho0)

    TYPE(rho0_mpole_type), POINTER    :: rho0

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_rho0_mpole"

    INTEGER                     :: istat

!   ---------------------------------------------------------------------------

    IF(ASSOCIATED(rho0)) THEN
      CALL deallocate_rho0_mpole(rho0)
    END IF

    ALLOCATE (rho0,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                            "allocate_rho0_mpole",int_size)

    NULLIFY(rho0%mp_rho)
    NULLIFY(rho0%mp_gau)
    NULLIFY(rho0%norm_g0l_h)
    NULLIFY(rho0%norm_g0l_s)
    NULLIFY(rho0%lmax0_kind)

  END SUBROUTINE allocate_rho0_mpole

! ******************************************************************************

  SUBROUTINE deallocate_mpole_gau(mp_gau)

    TYPE(mpole_gau_overlap), DIMENSION(:), POINTER :: mp_gau

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_mpole_gau"

    INTEGER                     :: ikind,istat,nkind

!   ---------------------------------------------------------------------------

    nkind = SIZE(mp_gau)

    DO ikind = 1,nkind
      DEALLOCATE(mp_gau(ikind)%Qlm_gg,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                      "Qlm_gg")
    END DO

    DEALLOCATE(mp_gau, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                      "mp_gau")

  END SUBROUTINE deallocate_mpole_gau

! ******************************************************************************

  SUBROUTINE deallocate_mpole_rho(mp_rho)

    TYPE(mpole_rho_atom), DIMENSION(:), POINTER :: mp_rho

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_mpole_rho"

    INTEGER                     :: iat,istat,natom

!   ---------------------------------------------------------------------------

    natom = SIZE(mp_rho)

    DO iat = 1,natom
      DEALLOCATE(mp_rho(iat)%Qlm_h,STAT=istat)
      DEALLOCATE(mp_rho(iat)%Qlm_s,STAT=istat) 
      DEALLOCATE(mp_rho(iat)%Qlm_tot,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                      "Qlm_h,Qlm_s,Qlm_tot")
    END DO

    DEALLOCATE(mp_rho,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                      "mp_rho")

  END SUBROUTINE deallocate_mpole_rho

! ******************************************************************************

  SUBROUTINE deallocate_rho0_atom(rho0_atom_set)

    TYPE(rho0_atom_type), DIMENSION(:), &
                                   POINTER :: rho0_atom_set
    
    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_rho0_atom"

    INTEGER                                :: iat, istat, natom

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(rho0_atom_set)) THEN

      natom = SIZE(rho0_atom_set)

      DO iat = 1,natom
        IF(ASSOCIATED(rho0_atom_set(iat)%rho0_rad_h)) THEN
          DEALLOCATE (rho0_atom_set(iat)%rho0_rad_h%r_coef)
          DEALLOCATE (rho0_atom_set(iat)%rho0_rad_h)
        ENDIF
        IF(ASSOCIATED(rho0_atom_set(iat)%rho0_rad_s)) THEN
          DEALLOCATE (rho0_atom_set(iat)%rho0_rad_s%r_coef)
          DEALLOCATE (rho0_atom_set(iat)%rho0_rad_s)
        ENDIF
      ENDDO

      DEALLOCATE (rho0_atom_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "rho0_atom_set")
    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer rho0_atom_set is not associated and "//&
                        "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_rho0_atom
! ******************************************************************************
  SUBROUTINE deallocate_rho0_mpole(rho0)

    TYPE(rho0_mpole_type), POINTER    :: rho0

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_rho0_mpole"

    INTEGER                     :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(rho0)) THEN

      IF (ASSOCIATED(rho0%mp_gau)) CALL deallocate_mpole_gau(rho0%mp_gau)

      IF (ASSOCIATED(rho0%mp_rho)) CALL deallocate_mpole_rho(rho0%mp_rho)

      IF (ASSOCIATED(rho0%lmax0_kind)) THEN
         DEALLOCATE(rho0%lmax0_kind,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                         "rho0%lmax0_kind")
      END IF

      IF (ASSOCIATED(rho0%norm_g0l_h)) THEN
         DEALLOCATE(rho0%norm_g0l_h, STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                         "rho0%norm_g0l_h")
      END IF

      IF (ASSOCIATED(rho0%norm_g0l_s)) THEN
         DEALLOCATE(rho0%norm_g0l_s, STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                         "rho0%norm_g0l_s")
      END IF

      DEALLOCATE (rho0,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "rho0_mpole")
    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer rho0 is not associated and "//&
                        "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_rho0_mpole

! ******************************************************************************
  SUBROUTINE initialize_mpole_rho(mp_rho,nchannels,zeff)

    TYPE(mpole_rho_atom)         :: mp_rho
    INTEGER , INTENT(IN)         :: nchannels , zeff

    CHARACTER(LEN=*), PARAMETER :: routine_name = "initialize_mpole_rho"

!   ---------------------------------------------------------------------------

   CALL reallocate(mp_rho%Qlm_h,1,nchannels)
   CALL reallocate(mp_rho%Qlm_s,1,nchannels)
   CALL reallocate(mp_rho%Qlm_tot,1,nchannels)

   mp_rho%Qlm_h = 0.0_dp
   mp_rho%Qlm_s = 0.0_dp
   mp_rho%Qlm_tot = 0.0_dp
   mp_rho%Qlm_z = -2._dp*ROOTPI*Zeff


  END SUBROUTINE initialize_mpole_rho
! ******************************************************************************

  SUBROUTINE write_rho0_info(rho0_mpole,globenv)

    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_rho0_info"

    INTEGER                                  :: ikind,l,nkind,output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (ASSOCIATED(rho0_mpole)) THEN
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,T5,A,/)") &
               "*** Compensation density charges data set ***"
      WRITE (UNIT=output_unit,FMT="(T2,A,T35,f8.3)")&
            "- Hard exponent :",rho0_mpole%zet0_h
      WRITE (UNIT=output_unit,FMT="(T2,A,T35,f8.3)")&
            "- Soft exponent :",rho0_mpole%zet0_s
      WRITE (UNIT=output_unit,FMT="(T2,A,T35,I5)")&
            "- Global max l :",rho0_mpole%lmax_0
      WRITE (UNIT=output_unit,FMT="(T2,A,T25,A4,T35,I2,T40,I5)")&
            "- l max per kind :","kind", 1,rho0_mpole%lmax0_kind(1)
      nkind = SIZE(rho0_mpole%lmax0_kind)
      DO ikind = 1,nkind
        WRITE (UNIT=output_unit,FMT="(T25,A4,T35,I2,T40,I5)")&
                      "kind",ikind,rho0_mpole%lmax0_kind(ikind)  
      END DO
      WRITE (UNIT=output_unit,FMT="(T2,A)")&
            "- Normalization constants for soft and hard g0"
      DO l = 0,rho0_mpole%lmax_0
        WRITE (UNIT=output_unit,FMT="(T20,A,T35,I2,T38,A,f10.5,A,f10.5)")&
            "ang. mom.= ", l, " hard= ", rho0_mpole%norm_g0l_h(l),&
                              " soft= ", rho0_mpole%norm_g0l_s(l)
      END DO
    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
        "The pointer rho0_mpole or rho0_atom_set is not associated",globenv)
    END IF

  END SUBROUTINE write_rho0_info
! ******************************************************************************
END MODULE qs_rho0_types
