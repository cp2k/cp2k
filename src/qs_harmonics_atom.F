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

  USE mathconstants,                   ONLY: degree, pi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset,&
                                             soset
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

     INTEGER                              :: max_s_harm, &
                                             max_iso_not0, &
                                             maxcg, ngrid
     INTEGER, DIMENSION(:,:,:), POINTER   :: ncg
     INTEGER, DIMENSION(:,:,:,:), POINTER :: ind_cg_a,ind_cg_b
     REAL(dp), DIMENSION(:,:),    POINTER :: slm
     REAL(dp), DIMENSION(:,:,:),  POINTER :: dslm
     REAL(dp), DIMENSION(:,:,:),  POINTER :: my_CG
     REAL(dp), DIMENSION(:), POINTER      :: slm_int

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
    NULLIFY(harmonics%slm_int)
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

      IF(ASSOCIATED(harmonics%slm_int)) THEN
        DEALLOCATE (harmonics%slm_int,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%slm_int")
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
                                   llmax,max_s_harm,ll,wa,&
                                   ang1,cos_ang2,cotan_ang2,globenv)

    TYPE(harmonics_atom_type), POINTER     :: harmonics
    TYPE(gto_basis_set_type), POINTER      :: orb_basis
    REAL(dp), DIMENSION(:,:,:), POINTER    :: my_CG
    INTEGER,  INTENT(IN)                   :: ll,na,llmax,max_s_harm 
    TYPE(global_environment_type)          :: globenv
    REAL(dp), DIMENSION(:), INTENT(IN)     :: wa
    REAL(dp), DIMENSION(:), INTENT(IN)     :: ang1
    REAL(dp), DIMENSION(:), INTENT(IN)     :: cos_ang2, cotan_ang2

    CHARACTER(len=*), PARAMETER :: routine_name = 'create_harmonics_atom', &
      routineP = module_name//':'//routine_name

    REAL(dp), DIMENSION(:), POINTER      :: y , slm_int
    REAL(dp), DIMENSION(:,:), POINTER    :: slm
    REAL(dp), DIMENSION(:,:,:), POINTER  :: dslm
    INTEGER, DIMENSION(:), POINTER       :: lmax, lmin
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind_cg_a,ind_cg_b

    INTEGER   :: i, ia, ia_bar, incg, is1, is2, iso, iso_mm, &
                 iso_mm1, iso_mp1, iso1, iso2, ja, l, l1, l2, &
                 m, mm, maxcg, maxl, maxs, nset,max_iso_not0
    REAL(dp)  :: abs_m, ang1_bar, l_r, m_r
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
      CALL reallocate(harmonics%slm_int,1,max_s_harm)

      slm => harmonics%slm
      dslm => harmonics%dslm
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

! Calculate the derivatives of theharminics with respect of the 2 angles
! the first angle is acos(lebedev_grid(ll)%r(3))
! the second angle is atan(lebedev_grid(ll)%r(2)/lebedev_grid(ll)%r(1))

!dbg
!      DO iso = 1,max_s_harm
!        l = indso(1,iso)
!        m = indso(2,iso)
!        l_r = REAL(l,dp)
!        m_r = REAL(m,dp)
!        DO ia = 1,na
!          IF(ABS(COS( (m_r+1._dp) * ang1(ia) )) < 1.E-8_dp) THEN
!            DO ja = 1,na
!              IF (ABS(COS( (m_r+1._dp) * ang1(ja))) > 1.E-8_dp .AND. &
!                ABS(cos_ang2(ia)-cos_ang2(ja)) < 1.E-8_dp) THEN
!                ang1_bar = ang1(ja)
!                ia_bar = ja
!    WRITE(*,'(4I3,4f20.10)') l, m+1, ia, ja, &
!            ang1(ia)*degree,ang1(ja)*degree,cos_ang2(ia),cos_ang2(ja)
!                GOTO 10
!              END IF
!            END DO
!              WRITE(*,'(3I3)') l, m+1, ia
!              STOP 'non trovato'
!10          CONTINUE
!          END IF
!        END DO
!      END DO
!  STOP 'ia ja'
!dbg

      DO iso = 1,max_s_harm

        l = indso(1,iso)
        m = indso(2,iso)
        l_r = REAL(l,dp)
        m_r = REAL(m,dp)


!  First Angle
        IF(m==0) THEN
          harmonics%dslm(1,1:na,iso) = 0.0_dp
        ELSE
          iso_mm = soset(l,-m)
          harmonics%dslm(1,1:na,iso) = - m_r * harmonics%slm(1:na,iso_mm)
        END IF

! Second Angle
        IF(m > 0) THEN

          harmonics%dslm(2,1:na,iso) = m_r * cotan_ang2(1:na)*&
                                       harmonics%slm(1:na,iso)
          IF(m /= l) THEN
            iso_mp1 = soset(l,m+1)
            DO ia = 1,na 
              IF(ABS(lebedev_grid(ll)%r(2,ia)) < 1.E-8_dp ) THEN
                harmonics%dslm(2,ia,iso) = 0.0_dp
              ELSE IF(ABS(COS( (m_r+1._dp) * ang1(ia) )) < 1.E-8_dp) THEN

!               Find a point where COS((m_r+1)*ang1) is different from zero
                DO ja = 1,na
                  IF (ABS(COS( (m_r+1._dp) * ang1(ja))) > 1.E-8_dp .AND. &
                       ABS(cos_ang2(ia)-cos_ang2(ja)) < 1.E-8_dp) THEN
                    ang1_bar = ang1(ja)
                    ia_bar = ja
                    EXIT
                  END IF
                END DO 


                harmonics%dslm(2,ia,iso) = harmonics%dslm(2,ia,iso) +&
                    COS( m_r * ang1(ia) )/COS( (m_r+1._dp) * ang1_bar ) * &
                    SQRT( (l_r-m_r) * (l_r+m_r+1._dp) ) * harmonics%slm(ia_bar,iso_mp1)
              ELSE
                harmonics%dslm(2,ia,iso) = harmonics%dslm(2,ia,iso) +&
                    COS( m_r * ang1(ia) )/COS( (m_r+1._dp) * ang1(ia) )*&
                    SQRT( (l_r-m_r) * (l_r+m_r+1._dp) ) * harmonics%slm(ia,iso_mp1)
              END IF
            END DO
          END IF

        ELSEIF(m == 0) THEN
          IF(l == 0) THEN
            harmonics%dslm(2,1:na,iso) = 0.0_dp
          ELSE 
            iso_mp1 = soset(l,1)
            DO ia = 1,na
              IF(ABS(lebedev_grid(ll)%r(2,ia)) < 1.E-8_dp ) THEN
                harmonics%dslm(2,ia,iso) = 0.0_dp
              ELSE IF(ABS(COS( ang1(ia))) < 1.E-8_dp) THEN

!               Find a point where COS((m_r+1)*ang1) is different from zero
                DO ja = 1,na
                  IF (ABS(COS( ang1(ja))) > 1.E-8_dp .AND. &
                       ABS(cos_ang2(ia)-cos_ang2(ja)) < 1.E-8_dp) THEN
                    ang1_bar = ang1(ja)
                    ia_bar = ja
                    EXIT
                  END IF
                END DO 

                harmonics%dslm(2,ia,iso) = SQRT(0.5_dp * l_r * (l_r+1._dp) )*&
                                  harmonics%slm(ia_bar,iso_mp1)/COS(ang1_bar)
              ELSE
                harmonics%dslm(2,ia,iso) = SQRT(0.5_dp * l_r * (l_r+1._dp) )*&
                                  harmonics%slm(ia,iso_mp1)/COS(ang1(ia))
              END IF
            END DO
          END IF
        ELSE
          abs_m = ABS(m_r) 
          harmonics%dslm(2,1:na,iso) = abs_m * cotan_ang2(1:na)*&
                                      harmonics%slm(1:na,iso)
          IF(m /= -l) THEN
            iso_mm1 = soset(l,m-1)
            DO ia = 1,na
              IF(ABS(lebedev_grid(ll)%r(2,ia)) < 1.E-8_dp .OR.&
                                  ABS(SIN( ang1(ia) )) < 1.E-8_dp) THEN 
                 harmonics%dslm(2,ia,iso) = 0.0_dp
              ELSE IF(ABS(SIN( abs_m * ang1(ia) )) < 1.E-8_dp) THEN
!                 Nothing
              ELSE IF(ABS(SIN((abs_m+1.0_dp) * ang1(ia) )) < 1.E-8_dp) THEN
!               Find a point where SIN((m_r+1)*ang1) is different from zero
                DO ja = 1,na
                  IF (ABS(SIN(  (abs_m+1.0_dp) * ang1(ja))) > 1.E-8_dp .AND. &
                       ABS(cos_ang2(ia)-cos_ang2(ja)) < 1.E-8_dp) THEN
                    ang1_bar = ang1(ja)
                    ia_bar = ja
                    EXIT
                  END IF
                END DO 
                harmonics%dslm(2,ia,iso) = harmonics%dslm(2,ia,iso) +&
                   SIN( abs_m * ang1(ia) )/SIN( (abs_m+1.0_dp) * ang1_bar )*&
                   SQRT( (l_r-abs_m) * (l_r+abs_m+1._dp) )* harmonics%slm(ia_bar,iso_mm1)
              ELSE
                harmonics%dslm(2,ia,iso) = harmonics%dslm(2,ia,iso) +&
                   SIN( abs_m * ang1(ia) )/SIN( (abs_m+1.0_dp) * ang1(ia) )*&
                   SQRT( (l_r-abs_m) * (l_r+abs_m+1._dp) )* harmonics%slm(ia,iso_mm1)
              END IF
            END DO
          END IF

        END IF

      END DO 

!  STOP 'dslm'
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
      max_iso_not0 = 0
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
                    max_iso_not0 = MAX(max_iso_not0,iso)
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
      harmonics%max_iso_not0 = max_iso_not0
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
