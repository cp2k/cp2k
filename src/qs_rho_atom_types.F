!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_rho_atom_types [1.0] *
!!
!!   NAME
!!     qs_rho_atom_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     
!!
!!   SOURCE
!******************************************************************************

MODULE qs_rho_atom_types

  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indso,&
                                             ncoset
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE util,                            ONLY: locate,&
                                             sort

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_rho_atom_types"

! *** Define the rho_atom_type ***

  TYPE rho_atom_coeff
    REAL(dp), DIMENSION(:,:), POINTER :: r_coef
  END TYPE rho_atom_coeff


  TYPE rho_atom_type
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER  :: cpc_h, cpc_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER  :: rho_rad_h,&
                                           rho_rad_s,&
                                           vrho_rad_h,&
                                           vrho_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:),&
                               POINTER  :: drho_rad_h,&
                                           drho_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:,:),&
                               POINTER  :: rho_rad_h_d,&
                                           rho_rad_s_d

    INTEGER                             :: rhoa_of_atom
    REAL(dp)                            :: exc_h,&
                                           exc_s
    TYPE(rho_atom_coeff),DIMENSION(:),&
                              POINTER   :: ga_Vxc_gb_h,&
                                           ga_Vxc_gb_s
  END TYPE rho_atom_type


! *** Public subroutines ***

  PUBLIC :: allocate_rho_atom_set,&
            deallocate_rho_atom_set,&
            get_rho_atom


! *** Public data types ***

  PUBLIC :: rho_atom_coeff,rho_atom_type

! *****************************************************************************

 CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_rho_atom_set(rho_atom_set,natom)

!   ***************************************************************************

    TYPE(rho_atom_type), DIMENSION(:), & 
      POINTER                          :: rho_atom_set
    INTEGER, INTENT(IN)                :: natom

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_rho_atom_set"

    INTEGER                            :: iat, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(rho_atom_set)) THEN
      CALL deallocate_rho_atom_set(rho_atom_set)
    END IF

    ALLOCATE (rho_atom_set(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                            "allocate_rho_atom_set",natom*int_size)

    DO iat = 1,natom

      NULLIFY (rho_atom_set(iat)%cpc_h)
      NULLIFY (rho_atom_set(iat)%cpc_s)
      NULLIFY (rho_atom_set(iat)%rho_rad_h)
      NULLIFY (rho_atom_set(iat)%rho_rad_s)
      NULLIFY (rho_atom_set(iat)%drho_rad_h)
      NULLIFY (rho_atom_set(iat)%drho_rad_s)
      NULLIFY (rho_atom_set(iat)%rho_rad_h_d)
      NULLIFY (rho_atom_set(iat)%rho_rad_s_d)
      NULLIFY (rho_atom_set(iat)%vrho_rad_h)
      NULLIFY (rho_atom_set(iat)%vrho_rad_s)
      NULLIFY (rho_atom_set(iat)%ga_Vxc_gb_h)
      NULLIFY (rho_atom_set(iat)%ga_Vxc_gb_s)

    ENDDO

  END SUBROUTINE allocate_rho_atom_set

! *****************************************************************************

  SUBROUTINE deallocate_rho_atom_set(rho_atom_set) 

    TYPE(rho_atom_type), DIMENSION(:), &   
      POINTER                          :: rho_atom_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_rho_atom_set"

    INTEGER                            :: i, iat, istat,j, n, natom

! *****************************************************************************

    IF (ASSOCIATED(rho_atom_set)) THEN

      natom = SIZE(rho_atom_set)

      DO iat = 1,natom
        IF(ASSOCIATED(rho_atom_set(iat)%cpc_h)) THEN
          IF(ASSOCIATED(rho_atom_set(iat)%cpc_h(1)%r_coef)) THEN
            n = SIZE(rho_atom_set(iat)%cpc_h,1)
            DO i = 1,n 
              deallocate(rho_atom_set(iat)%cpc_h(i)%r_coef)
              deallocate(rho_atom_set(iat)%cpc_s(i)%r_coef)
              deallocate(rho_atom_set(iat)%ga_Vxc_gb_h(i)%r_coef)
              deallocate(rho_atom_set(iat)%ga_Vxc_gb_s(i)%r_coef)
              DEALLOCATE(rho_atom_set(iat)%drho_rad_h(i)%r_coef)
              DEALLOCATE(rho_atom_set(iat)%drho_rad_s(i)%r_coef)
              DO j = 1,3
                DEALLOCATE(rho_atom_set(iat)%rho_rad_h_d(j,i)%r_coef)
                DEALLOCATE(rho_atom_set(iat)%rho_rad_s_d(j,i)%r_coef)
              END DO
            ENDDO
          END IF 
          DEALLOCATE (rho_atom_set(iat)%cpc_h)
          DEALLOCATE (rho_atom_set(iat)%cpc_s)
          DEALLOCATE (rho_atom_set(iat)%drho_rad_h)
          DEALLOCATE (rho_atom_set(iat)%drho_rad_s)
          DEALLOCATE (rho_atom_set(iat)%rho_rad_h_d)
          DEALLOCATE (rho_atom_set(iat)%rho_rad_s_d)
          DEALLOCATE (rho_atom_set(iat)%ga_Vxc_gb_h)
          DEALLOCATE (rho_atom_set(iat)%ga_Vxc_gb_s)

        ENDIF

        IF(ASSOCIATED(rho_atom_set(iat)%rho_rad_h)) THEN
          n = SIZE(rho_atom_set(iat)%rho_rad_h)
          DO i = 1,n 
            deallocate(rho_atom_set(iat)%rho_rad_h(i)%r_coef)
          ENDDO
          DEALLOCATE (rho_atom_set(iat)%rho_rad_h)
        ENDIF

        IF(ASSOCIATED(rho_atom_set(iat)%rho_rad_s)) THEN
          n = SIZE(rho_atom_set(iat)%rho_rad_s)
          DO i = 1,n 
            deallocate(rho_atom_set(iat)%rho_rad_s(i)%r_coef)
          ENDDO
          DEALLOCATE (rho_atom_set(iat)%rho_rad_s)
        ENDIF

        IF(ASSOCIATED(rho_atom_set(iat)%vrho_rad_h)) THEN
          n = SIZE(rho_atom_set(iat)%vrho_rad_h)
          DO i = 1,n 
            deallocate(rho_atom_set(iat)%vrho_rad_h(i)%r_coef)
          ENDDO
          DEALLOCATE (rho_atom_set(iat)%vrho_rad_h)
        ENDIF

        IF(ASSOCIATED(rho_atom_set(iat)%vrho_rad_s)) THEN
          n = SIZE(rho_atom_set(iat)%vrho_rad_s)
          DO i = 1,n 
            deallocate(rho_atom_set(iat)%vrho_rad_s(i)%r_coef)
          ENDDO
          DEALLOCATE (rho_atom_set(iat)%vrho_rad_s)
        ENDIF

      ENDDO

      DEALLOCATE (rho_atom_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "rho_atom_set")

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer rho_atom_set is not associated and "//&
                        "cannot be deallocated")

    END IF

  END SUBROUTINE deallocate_rho_atom_set

! *****************************************************************************

  SUBROUTINE get_rho_atom(rho_atom,cpc_h,cpc_s,rho_rad_h,rho_rad_s,&
                          drho_rad_h,drho_rad_s,vrho_rad_h,vrho_rad_s,&
                          ga_Vxc_gb_h,ga_Vxc_gb_s,rho_rad_h_d,rho_rad_s_d)

    TYPE(rho_atom_type), POINTER        :: rho_atom
    TYPE(rho_atom_coeff), DIMENSION(:), &
                     POINTER, OPTIONAL  :: cpc_h, cpc_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
                      POINTER, OPTIONAL :: rho_rad_h,&
                                           rho_rad_s,&
                                           vrho_rad_h,&
                                           vrho_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
                      POINTER, OPTIONAL :: drho_rad_h,&
                                           drho_rad_s
    TYPE(rho_atom_coeff), DIMENSION(:,:), &
                      POINTER, OPTIONAL :: rho_rad_h_d,&
                                           rho_rad_s_d

    TYPE(rho_atom_coeff),DIMENSION(:),POINTER, &
                               OPTIONAL :: ga_Vxc_gb_h,&
                                           ga_Vxc_gb_s

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_rho_atom"

    INTEGER                                  :: iset

! *****************************************************************************
    IF (ASSOCIATED(rho_atom)) THEN
      IF (PRESENT(cpc_h)) cpc_h => rho_atom%cpc_h
      IF (PRESENT(cpc_s)) cpc_s => rho_atom%cpc_s
      IF (PRESENT(rho_rad_h)) rho_rad_h => rho_atom%rho_rad_h
      IF (PRESENT(rho_rad_s)) rho_rad_s => rho_atom%rho_rad_s
      IF (PRESENT(drho_rad_h)) drho_rad_h => rho_atom%drho_rad_h
      IF (PRESENT(drho_rad_s)) drho_rad_s => rho_atom%drho_rad_s
      IF (PRESENT(rho_rad_h_d)) rho_rad_h_d => rho_atom%rho_rad_h_d
      IF (PRESENT(rho_rad_s_d)) rho_rad_s_d => rho_atom%rho_rad_s_d
      IF (PRESENT(vrho_rad_h)) vrho_rad_h => rho_atom%vrho_rad_h
      IF (PRESENT(vrho_rad_s)) vrho_rad_s => rho_atom%vrho_rad_s
      IF (PRESENT(ga_Vxc_gb_h)) ga_Vxc_gb_h => rho_atom%ga_Vxc_gb_h
      IF (PRESENT(ga_Vxc_gb_s)) ga_Vxc_gb_s => rho_atom%ga_Vxc_gb_s
    ELSE
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer rho_atom is not associated")
    END IF

  END SUBROUTINE  get_rho_atom

!****************************************************************************
END MODULE qs_rho_atom_types
