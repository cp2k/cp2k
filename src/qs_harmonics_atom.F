!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_harmonics_atom [1.0] *
!!
!!   NAME
!!     qs_harmonics_atom
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     
!!
!!   SOURCE
!******************************************************************************

MODULE qs_harmonics_atom

  USE basis_set_types,                 ONLY: get_gto_basis_set, &
                                             gto_basis_set_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE lebedev,                         ONLY: lebedev_grid

  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset
  USE spherical_harmonics,             ONLY: clebsch_gordon,&
                                             clebsch_gordon_init,&
                                             y_lm
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_harmonics_atom"

  TYPE harmonics_atom_type

     INTEGER                         :: max_s_harm, maxcg, ngrid
     INTEGER, DIMENSION(:,:,:), &
                           POINTER   :: ncg
     INTEGER, DIMENSION(:,:,:,:), &
                           POINTER   :: ind_cg_a,ind_cg_b
     REAL(dp), DIMENSION(:,:), &
                             POINTER :: slm
     REAL(dp), DIMENSION(:,:,:), &
                             POINTER :: dslm
     REAL(dp), DIMENSION(:,:,:), &
                             POINTER :: my_CG

  END TYPE harmonics_atom_type


! *** Public subroutines ***

  PUBLIC :: allocate_harmonics_atom, &
            create_harmonics_atom, &
            deallocate_harmonics_atom


! *** Public data types ***

  PUBLIC :: harmonics_atom_type


!!***
! *****************************************************************************

 CONTAINS

! *****************************************************************************


  SUBROUTINE allocate_harmonics_atom(harmonics)

!   Purpose: Allocate a spherical harmonics set for the atom grid.

!   History: - Creation ()

!   ***************************************************************************

    TYPE(harmonics_atom_type), POINTER        ::  harmonics

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_harmonics_atom"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(harmonics)) CALL deallocate_harmonics_atom(harmonics)

    ALLOCATE (harmonics,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics",0)

    NULLIFY(harmonics%ind_cg_a)
    NULLIFY(harmonics%ind_cg_b)
    NULLIFY(harmonics%slm)
    NULLIFY(harmonics%dslm)
    NULLIFY(harmonics%ncg)
    NULLIFY(harmonics%my_CG)

  END SUBROUTINE allocate_harmonics_atom

! *****************************************************************************


  SUBROUTINE deallocate_harmonics_atom(harmonics)

!   Purpose: Deallocate the spherical harmonics set for the atom grid.

!   History: - Creation ()

!   ***************************************************************************

    TYPE(harmonics_atom_type), POINTER        ::  harmonics

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_harmonics_atom"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(harmonics)) THEN 

      IF(ASSOCIATED(harmonics%ind_cg_a)) THEN
        DEALLOCATE (harmonics%ind_cg_a,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%ind_cg_a")
      ENDIF
      IF(ASSOCIATED(harmonics%ind_cg_b)) THEN
        DEALLOCATE (harmonics%ind_cg_b,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%ind_cg_b")
      ENDIF
      IF(ASSOCIATED(harmonics%slm)) THEN
        DEALLOCATE (harmonics%slm,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%slm")
      ENDIF
      IF(ASSOCIATED(harmonics%dslm)) THEN
        DEALLOCATE (harmonics%dslm,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%dslm")
      ENDIF
      IF(ASSOCIATED(harmonics%ncg)) THEN
        DEALLOCATE (harmonics%ncg,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%ncg")
      ENDIF
      IF(ASSOCIATED(harmonics%my_CG)) THEN
        DEALLOCATE (harmonics%my_CG,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%my_CG")
      ENDIF
      DEALLOCATE (harmonics,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics")
    ELSE
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer harmonics is not associated and "//&
                        "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_harmonics_atom

! *****************************************************************************

  SUBROUTINE create_harmonics_atom(harmonics,orb_basis,my_CG,na,&
                                   llmax,max_s_harm,ll,globenv)

    TYPE(harmonics_atom_type), POINTER     :: harmonics
    TYPE(gto_basis_set_type), POINTER      :: orb_basis
    REAL(dp), DIMENSION(:,:,:), POINTER    :: my_CG
    INTEGER                                :: ll,na,llmax,max_s_harm 
    TYPE(global_environment_type)          :: globenv

    CHARACTER(len=*), PARAMETER :: routine_name = 'create_harmonics_atom', &
      routineP = module_name//':'//routine_name

    REAL(dp), DIMENSION(:), POINTER      :: y 
    REAL(dp), DIMENSION(:,:), POINTER    :: slm
    REAL(dp), DIMENSION(:,:,:), POINTER  :: dslm
    INTEGER, DIMENSION(:), POINTER       :: lmax, lmin
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind_cg_a,ind_cg_b

    INTEGER   :: i, incg, is1, is2, iso, iso1, iso2, l, l1, l2, &
                 m, mm, maxcg, maxl, maxs, nset

! *****************************************************************************

    NULLIFY (y,slm,dslm)

    IF (ASSOCIATED(harmonics)) THEN
      CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin, &
                           maxl=maxl,nset=nset)

      harmonics%max_s_harm = max_s_harm
      harmonics%ngrid = na
      maxs = nsoset(maxl)

      NULLIFY(harmonics%my_CG)
      CALL reallocate(harmonics%my_CG,1,maxs,1,maxs,1,max_s_harm)
      DO i = 1,max_s_harm
        DO is1 = 1, maxs
          harmonics%my_CG(1:maxs,is1,i) =  my_CG(1:maxs,is1,i)
        END DO
      END Do

!   *** allocate and calculate the spherical harmonics LM for this grid ***
!       and their derivatives 
      call reallocate(harmonics%slm,1,na,1,max_s_harm)
      call reallocate(harmonics%dslm,1,2,1,na,1,max_s_harm)
      call reallocate(harmonics%ncg,1,max_s_harm,1,nset,1,nset)
      call reallocate(y,1,na)

      slm => harmonics%slm
      dslm => harmonics%dslm

      DO iso = 1,max_s_harm

        l = indso(1,iso)
        m = indso(2,iso)

        CALL y_lm(lebedev_grid(ll)%r,y,l,m)
        slm(1:na,iso) = y(1:na)

        DO i =1,2
!!        TO DO
!!        calculate the derivatives of the spherical harmonics
!!        with respect to teta (i=1) and to phi (i=2)
!!        where   cos(teta)          = lebedev_grid(ll)%r(3)
!!                cos(phi)*sin(teta) = lebedev_grid(ll)%r(1)
!!                sin(phi)*sin(teta) = lebedev_grid(ll)%r(2)
!!        CALL dy_lm(lebedev_grid(ll)%r,dy,l,m,i)
          dslm(i,1:na,iso) =   y(1:na)  !dydteta(1:na) 
        END DO 
      END DO

!   *** Assign indexes for the non null CG coefficients ***      
      mm = nsoset(maxl)*nsoset(maxl)
      ALLOCATE (ind_cg_a(mm))
      ALLOCATE (ind_cg_b(mm))

      NULLIFY(harmonics%ind_cg_a,harmonics%ind_cg_b)
            CALL reallocate(harmonics%ind_cg_a,&
                        1,mm,1,max_s_harm,1,nset,1,nset)
            CALL reallocate(harmonics%ind_cg_b,&
                        1,mm,1,max_s_harm,1,nset,1,nset)
      maxcg = 0
      DO is1 = 1,nset
        DO is2 = 1,nset
          DO iso = 1,max_s_harm
            incg = 0
            ind_cg_a = 0
            ind_cg_b = 0
            DO l1 = lmin(is1),lmax(is1)
              DO iso1 = nsoset(l1-1)+1,nsoset(l1)
                DO l2 = lmin(is2), lmax(is2)
                  IF(l1+l2 > llmax) CYCLE
                  DO iso2 = nsoset(l2-1)+1,nsoset(l2)
                    IF(abs(my_CG(iso1,iso2,iso)) .LT. 1.E-8_dp) CYCLE
                    incg = incg + 1
                    ind_cg_a(incg) = iso1
                    ind_cg_b(incg) = iso2
                  END DO
                END DO
              END DO
            END DO
            maxcg = max(maxcg,incg)

            harmonics%ind_cg_a(1:incg,iso,is1,is2) = &
                                    ind_cg_a(1:incg)
            harmonics%ind_cg_b(1:incg,iso,is1,is2) = &
                                    ind_cg_b(1:incg)
            harmonics%ncg(iso,is1,is2) = incg 
          END DO
        END DO ! is2
      END DO ! is1
      harmonics%maxcg =  maxcg
      DEALLOCATE(ind_cg_a)
      DEALLOCATE(ind_cg_b)

    ELSE
      CALL stop_program(routine_name,module_name,__LINE__,&
               "The pointer harmonics is not associated",&
                globenv)
    ENDIF

  END SUBROUTINE  create_harmonics_atom
  
!*******************************************************************************
END MODULE qs_harmonics_atom
