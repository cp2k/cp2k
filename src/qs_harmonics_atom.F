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
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE lebedev,                         ONLY: lebedev_grid

  USE mathconstants,                   ONLY: degree, pi,twopi,fourpi,fac,dfac
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset,nso,&
                                             soset,& 
                                             ncoset,nco,indco
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE spherical_harmonics,             ONLY: clebsch_gordon,&
                                             clebsch_gordon_init,&
                                             y_lm,dPof1,legendre!,drvy_lm
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_harmonics_atom"

  TYPE harmonics_atom_type

     INTEGER                              :: max_s_harm, &
                                             max_iso_not0, &
                                             maxcg, maxcg_d, ngrid
     INTEGER, DIMENSION(:,:,:), POINTER   :: ncg, ncg_dxyz
     INTEGER, DIMENSION(:,:,:,:), POINTER :: ind_cg_a,ind_cg_b
     INTEGER, DIMENSION(:,:,:,:), POINTER :: ind_cg_a_dxyz,ind_cg_b_dxyz
     REAL(dp), DIMENSION(:,:),    POINTER :: a, slm
     REAL(dp), DIMENSION(:,:,:),  POINTER :: dslm, dslm_dxyz
     REAL(dp), DIMENSION(:,:,:),  POINTER :: my_CG
     REAL(dp), DIMENSION(:,:,:,:),  POINTER :: my_CG_dxyz
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
    NULLIFY(harmonics%ind_cg_a_dxyz)
    NULLIFY(harmonics%ind_cg_b_dxyz)
    NULLIFY(harmonics%slm)
    NULLIFY(harmonics%dslm)
    NULLIFY(harmonics%dslm_dxyz)
    NULLIFY(harmonics%slm_int)
    NULLIFY(harmonics%ncg)
    NULLIFY(harmonics%ncg_dxyz)
    NULLIFY(harmonics%my_CG)
    NULLIFY(harmonics%my_CG_dxyz)
    NULLIFY(harmonics%a)

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

      IF(ASSOCIATED(harmonics%ind_cg_a_dxyz)) THEN
        DEALLOCATE (harmonics%ind_cg_a_dxyz,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%ind_cg_a_dxyz")
      ENDIF

      IF(ASSOCIATED(harmonics%ind_cg_b_dxyz)) THEN
        DEALLOCATE (harmonics%ind_cg_b_dxyz,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%ind_cg_b_dxyz")
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

      IF(ASSOCIATED(harmonics%dslm_dxyz)) THEN
        DEALLOCATE (harmonics%dslm_dxyz,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%dslm_dxyz")
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
      IF(ASSOCIATED(harmonics%ncg_dxyz)) THEN
        DEALLOCATE (harmonics%ncg_dxyz,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%ncg_dxyz")
      ENDIF

      IF(ASSOCIATED(harmonics%my_CG)) THEN
        DEALLOCATE (harmonics%my_CG,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%my_CG")
      ENDIF

      IF(ASSOCIATED(harmonics%my_CG_dxyz)) THEN
        DEALLOCATE (harmonics%my_CG_dxyz,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%my_CG_dxyz")
      ENDIF

      IF(ASSOCIATED(harmonics%a)) THEN
        DEALLOCATE (harmonics%a,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "harmonics%a")
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
                                   azi,cos_azi,pol,cos_pol,cotan_pol)

    TYPE(harmonics_atom_type), POINTER     :: harmonics
    TYPE(gto_basis_set_type), POINTER      :: orb_basis
    REAL(dp), DIMENSION(:,:,:), POINTER    :: my_CG
    INTEGER,  INTENT(IN)                   :: ll,na,llmax,max_s_harm 
    REAL(dp), DIMENSION(:), INTENT(IN)     :: wa
    REAL(dp), DIMENSION(:), INTENT(IN)     :: azi, pol, cos_azi
    REAL(dp), DIMENSION(:), INTENT(IN)     :: cos_pol, cotan_pol

    CHARACTER(len=*), PARAMETER :: routine_name = 'create_harmonics_atom', &
      routineP = module_name//':'//routine_name

    REAL(dp), DIMENSION(:), POINTER      :: y , slm_int
    REAL(dp), DIMENSION(:,:), POINTER    :: dy, slm,dc
    REAL(dp), DIMENSION(:,:,:), POINTER  :: dslm, dslm_dxyz
    INTEGER, DIMENSION(:), POINTER       :: lmax, lmin
    INTEGER, DIMENSION(:), ALLOCATABLE   :: ind_cg_a,ind_cg_b
    INTEGER, DIMENSION(:), ALLOCATABLE   :: ind_cg_a_d,ind_cg_b_d

    INTEGER   :: i, ia, ia_bar, incg, incg_d, is1, is2, iso, iso_mm, &
                 iso_mm1, iso_mp1, iso1, iso2, ja, l, l1, l2, &
                 m, mm, maxcg, maxcg_d, maxl, maxs, nset,max_iso_not0,&
                 lx, ly, lz, n, ic, is
    REAL(dp)  :: abs_m, azi_bar, int, int1, int2, int3, l_r, m_r, rx, ry, rz, drx, dry, drz
    REAL(dp), DIMENSION(:,:,:),ALLOCATABLE :: ss
!dbg
    CHARACTER(len=80) :: findex,filename
    REAL(dp) :: fact,oorootfourpi,ooroottwopi,plm
    REAL(dp), DIMENSION(:,:),ALLOCATABLE :: my_slm 
! *****************************************************************************

    NULLIFY (y,dy,slm,dslm,dslm_dxyz,dc)

    IF (ASSOCIATED(harmonics)) THEN
      CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin, &
                           maxl=maxl,nset=nset)

      harmonics%max_s_harm = max_s_harm
      harmonics%ngrid = na
      maxs = nsoset(maxl)

      NULLIFY(harmonics%my_CG,harmonics%my_CG_dxyz)
      CALL reallocate(harmonics%my_CG,1,maxs,1,maxs,1,max_s_harm)
      CALL reallocate(harmonics%my_CG_dxyz,1,3,1,maxs,1,maxs,1,max_s_harm)
      DO i = 1,max_s_harm
        DO is1 = 1, maxs
          harmonics%my_CG(1:maxs,is1,i) =  my_CG(1:maxs,is1,i)
        END DO
      END Do

!   *** allocate and calculate the spherical harmonics LM for this grid ***
!       and their derivatives 
      call reallocate(harmonics%slm,1,na,1,max_s_harm)
      call reallocate(harmonics%dslm,1,2,1,na,1,max_s_harm)
      CALL reallocate(harmonics%dslm_dxyz,1,3,1,na,1,max_s_harm)
      call reallocate(harmonics%ncg,1,max_s_harm,1,nset,1,nset)
      call reallocate(harmonics%ncg_dxyz,1,max_s_harm,1,nset,1,nset)
      CALL reallocate(harmonics%a,1,3,1,na)
      call reallocate(y,1,na)
      CALL reallocate(dy,1,3,1,na)
      CALL reallocate(dc,1,nco(6),1,3)
      CALL reallocate(harmonics%slm_int,1,max_s_harm)

      slm => harmonics%slm
      dslm => harmonics%dslm
      dslm_dxyz => harmonics%dslm_dxyz
      slm_int => harmonics%slm_int
      slm_int = 0.0_dp

      DO iso = 1,max_s_harm

        l = indso(1,iso)
        m = indso(2,iso)

        CALL y_lm(lebedev_grid(ll)%r,y,l,m)
        slm(1:na,iso) = y(1:na)
!        CALL drvy_lm(lebedev_grid(ll)%r,dy,l,m)

        DO i = 1,na
!           dslm_dxyz(1:3,i,iso) = dy(1:3,i)
           slm_int(iso) = slm_int(iso) + slm(i,iso)*wa(i)
        END DO ! i

      END DO  ! iso

      DO i = 1,3
        harmonics%a(i,1:na) = lebedev_grid(ll)%r(i,1:na)
      END DO

! dbg test slm
  allocate(my_slm(na,max_s_harm))
  my_slm = 0.0_dp
  oorootfourpi = 1._dp/sqrt(fourpi)
  ooroottwopi = 1._dp/sqrt(twopi)
  do iso = 1,9
    l = indso(1,iso)
    m = indso(2,iso)

   write(findex,'(I4)') 1000+iso
    filename =  'slm_'//trim(findex)//'.dat'
!    open(10,file=filename,status="unknown",form="formatted")

    fact = sqrt((2._dp*real(l)+1._dp)*fac(l-abs(m))/fac(l+abs(m)))
    if(m>0) then
      fact = fact*ooroottwopi
      do ia = 1,na
        plm = legendre (cos_azi(ia),l,m)
        my_slm(ia,iso) = fact*plm*cos(real(m)*pol(ia)) 
!        write(10,'(I4,5f12.6,2f20.10)') ia, harmonics%a(1,ia) ,&
!                                            harmonics%a(2,ia) ,&
!                                            harmonics%a(3,ia) ,&
!                                            azi(ia)*degree,pol(ia)*degree,&
!                                            slm(ia,iso),my_slm(ia,iso)
      end do
    elseif(m==0) then
      fact = fact*oorootfourpi
      do ia = 1,na
        plm = legendre (cos_azi(ia),l,0)
        my_slm(ia,iso) = fact*plm
!        write(10,'(I4,5f12.6,2f20.10)') ia, harmonics%a(1,ia) ,&
!                                            harmonics%a(2,ia) ,&
!                                            harmonics%a(3,ia) ,&
!                                            azi(ia)*degree,pol(ia)*degree,&
!                                            slm(ia,iso),my_slm(ia,iso)
      end do
    else
      fact = fact*ooroottwopi
      do ia = 1,na
        plm = legendre (cos_azi(ia),l,abs(m))
        my_slm(ia,iso) = fact*plm*sin(real(abs(m))*pol(ia)) 
!        write(10,'(I4,5f12.6,2f20.10)') ia, harmonics%a(1,ia) ,&
!                                            harmonics%a(2,ia) ,&
!                                            harmonics%a(3,ia) ,&
!                                            azi(ia)*degree,pol(ia)*degree,&
!                                            slm(ia,iso),my_slm(ia,iso)
      end do
    end if
!    close(10)
  
  end do ! iso
!  stop 'check slm'
! dbg 

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

!  WRITE(*,*) indso(1,12), indso(2,12) 
!  DO i =1, na
!   WRITE(*,'(I4,4f20.10)')  i, slm(i,7),dslm_dxyz(1,i,7),dslm_dxyz(2,i,7),dslm_dxyz(3,i,7)
!  end do

! STOP 'dslm'
! Calculate the derivatives of theharminics with respect of the 2 angles
! the first angle is acos(lebedev_grid(ll)%r(3))
! the second angle is atan(lebedev_grid(ll)%r(2)/lebedev_grid(ll)%r(1))


      DO iso = 1,max_s_harm

        l = indso(1,iso)
        m = indso(2,iso)
        l_r = REAL(l,dp)
        m_r = REAL(m,dp)


!  Derivative with respect to the azimuthal angle
        IF(m==0) THEN
          harmonics%dslm(2,1:na,iso) = 0.0_dp
        ELSE
          iso_mm = soset(l,-m)
          harmonics%dslm(2,1:na,iso) = - m_r * harmonics%slm(1:na,iso_mm)
        END IF

! Derivative with respect to the polar angle
        IF(m > 0) THEN

          DO ia = 1,na 

            IF (m == 1 .AND. ABS(cos_pol(ia)-1.0_dp) < 1.E-8_dp ) THEN
              harmonics%dslm(1,ia,iso) = -SQRT( (2.0_dp*l_r) / (2.0_dp*pi*l_r*(l_r+1.0_dp) ) )*&
                             COS( azi(ia) )*SIGN( 1.0_dp,cos_pol(ia) ) * dPof1( cos_pol(ia),l )
            ELSEIF( ABS(cos_pol(ia)-1.0_dp) < 1.E-8_dp ) THEN
              harmonics%dslm(1,ia,iso) = 0.0_dp
            ELSE
              harmonics%dslm(1,ia,iso) = m_r * cotan_pol(ia)*harmonics%slm(ia,iso)
            END IF

            IF(m /= l) THEN
              iso_mp1 = soset(l,m+1)
              iso_mm1 = soset(l,-m-1)

              harmonics%dslm(1,ia,iso) = harmonics%dslm(1,ia,iso) +&
                    (COS( azi(ia) )*harmonics%slm(ia,iso_mp1) +&
                     SIN( azi(ia) )*harmonics%slm(ia,iso_mm1))*&
                    SQRT( (l_r-m_r) * (l_r+m_r+1._dp) )

            END IF
          END DO

        ELSEIF(m == 0) THEN
          IF(l == 0) THEN
            harmonics%dslm(1,1:na,iso) = 0.0_dp
          ELSE 
            iso_mp1 = soset(l,1)
            iso_mm1 = soset(l,-1)
            DO ia = 1,na

              IF(ABS(COS( azi(ia))) < 1.E-8_dp) THEN
                harmonics%dslm(1,ia,iso) = SQRT(0.5_dp * l_r * (l_r+1._dp) )*&
                                  harmonics%slm(ia,iso_mm1)*SIN(azi(ia))
              ELSE
                harmonics%dslm(1,ia,iso) = SQRT(0.5_dp * l_r * (l_r+1._dp) )*&
                                  harmonics%slm(ia,iso_mp1)/COS(azi(ia))
              END IF
            END DO
          END IF
        ELSE

          abs_m = ABS(m_r) 
          DO ia = 1,na 
            IF (ABS(m) == 1 .AND. ABS(cos_pol(ia)-1.0_dp) < 1.E-8_dp ) THEN
              harmonics%dslm(1,ia,iso) = -SQRT( (2.0_dp*l_r) / (2.0_dp*pi*l_r*( l_r+1.0_dp ) ))*&
                           SIN( azi(ia) )*SIGN( 1.0_dp,cos_pol(ia) ) * dPof1( cos_pol(ia),l )
            ELSEIF( ABS(cos_pol(ia)-1.0_dp) < 1.E-8_dp ) THEN
              harmonics%dslm(1,ia,iso) = 0.0_dp
            ELSE
              harmonics%dslm(1,ia,iso) = abs_m * cotan_pol(ia)*&
                                           harmonics%slm(ia,iso)
            END IF
            IF(m /= -l) THEN
              iso_mm1 = soset(l,m-1)
              iso_mm1 = soset(l,-m+1)

              harmonics%dslm(1,ia,iso) = harmonics%dslm(1,ia,iso) +&
                   (COS(  azi(ia) ) * harmonics%slm(ia,iso_mm1) -&
                    SIN( azi(ia) )* harmonics%slm(ia,iso_mp1) ) *&
                   SQRT( (l_r-abs_m) * (l_r+abs_m+1._dp) )

            END IF
          END DO

        END IF

      END DO 

!      DO iso = 1,max_s_harm
!        l = indso(1,iso)
!        m = indso(2,iso)
!        DO ia = 1,na
!          WRITE(25,'(2I4,2f12.5,3f20.16)') l,m,ACOS(cos_pol(ia)),azi(ia),&
!                    harmonics%slm(ia,iso),&
!                    harmonics%dslm(1,ia,iso),harmonics%dslm(2,ia,iso)
!        END DO 
!      END DO

!  STOP 'dslm'
!   *** Assign indexes for the non null CG coefficients ***      
      mm = nsoset(maxl)*nsoset(maxl)
      ALLOCATE (ind_cg_a(mm),ind_cg_a_d(mm))
      ALLOCATE (ind_cg_b(mm),ind_cg_b_d(mm))

      NULLIFY(harmonics%ind_cg_a,harmonics%ind_cg_b)
      NULLIFY(harmonics%ind_cg_a_dxyz,harmonics%ind_cg_b_dxyz)
      CALL reallocate(harmonics%ind_cg_a,&
                        1,mm,1,max_s_harm,1,nset,1,nset)
      CALL reallocate(harmonics%ind_cg_b,&
                        1,mm,1,max_s_harm,1,nset,1,nset)
      CALL reallocate(harmonics%ind_cg_a_dxyz,&
                        1,mm,1,max_s_harm,1,nset,1,nset)
      CALL reallocate(harmonics%ind_cg_b_dxyz,&
                        1,mm,1,max_s_harm,1,nset,1,nset)

      maxcg = 0
      maxcg_d = 0
      max_iso_not0 = 0
      DO is1 = 1,nset
        DO is2 = 1,nset
          DO iso = 1,max_s_harm
            incg = 0
            incg_d = 0
            ind_cg_a = 0
            ind_cg_b = 0
            ind_cg_a_d = 0
            ind_cg_b_d = 0
            DO l1 = lmin(is1),lmax(is1)
              DO iso1 = nsoset(l1-1)+1,nsoset(l1)
                DO l2 = lmin(is2), lmax(is2)
                  IF(l1+l2 > llmax) CYCLE
                  DO iso2 = nsoset(l2-1)+1,nsoset(l2)
                    IF(abs(my_CG(iso1,iso2,iso)) > 1.E-8_dp) THEN
                      incg = incg + 1
                      ind_cg_a(incg) = iso1
                      ind_cg_b(incg) = iso2
                      max_iso_not0 = MAX(max_iso_not0,iso)
                    END IF
                    rx = abs(harmonics%my_CG_dxyz(1,iso1,iso2,iso))+&
                         abs(harmonics%my_CG_dxyz(2,iso1,iso2,iso))+&
                         abs(harmonics%my_CG_dxyz(3,iso1,iso2,iso))
                    IF(rx > 1.E-8_dp) THEN
                      incg_d = incg_d + 1
                      ind_cg_a_d(incg_d) = iso1
                      ind_cg_b_d(incg_d) = iso2
                    END IF
                  END DO
                END DO
              END DO
            END DO
            maxcg = max(maxcg,incg)
            maxcg_d = MAX(maxcg_d,incg_d)

            harmonics%ind_cg_a(1:incg,iso,is1,is2) = &
                                    ind_cg_a(1:incg)
            harmonics%ind_cg_b(1:incg,iso,is1,is2) = &
                                    ind_cg_b(1:incg)
            harmonics%ncg(iso,is1,is2) = incg 

            harmonics%ind_cg_a_dxyz(1:incg_d,iso,is1,is2) = &
                                            ind_cg_a_d(1:incg_d)
            harmonics%ind_cg_b_dxyz(1:incg_d,iso,is1,is2) = &
                                            ind_cg_b_d(1:incg_d)
            harmonics%ncg_dxyz(iso,is1,is2) = incg_d

          END DO
        END DO ! is2
      END DO ! is1
      harmonics%maxcg =  maxcg
      harmonics%maxcg_d =  maxcg_d
      harmonics%max_iso_not0 = max_iso_not0
      DEALLOCATE(ind_cg_a)
      DEALLOCATE(ind_cg_b)
      DEALLOCATE(ind_cg_a_d)
      DEALLOCATE(ind_cg_b_d)


!dbg
!      ALLOCATE(ss(na,max_s_harm,max_s_harm))
!      ss = 0.0_dp
!      DO is1 = 1,nset
!        DO is2 = 1,nset
!          DO iso = 1,max_iso_not0
!            DO i = 1,harmonics%ncg_dxyz(iso,is1,is2)
!              iso1 = harmonics%ind_cg_a_dxyz(i,iso,is1,is2)
!              iso2 = harmonics%ind_cg_b_dxyz(i,iso,is1,is2)
!              l1 = indso(1,iso1) + indso(1,iso2)
!   DO iso1 = 1,maxs
!   DO iso2 = 1,maxs
!              ss(1:na,iso1,iso2) = ss(1:na,iso1,iso2)+&
!                        harmonics%my_CG_dxyz(2,iso1,iso2,iso)*slm(1:na,iso)
!              ss(1:na,iso1,iso2) = ss(1:na,iso1,iso2)+&
!                            my_CG(iso1,iso2,iso)*slm(1:na,iso)
!            END DO
!   END DO
!   END DO
!          END DO
!        END DO
!      END DO 
!          iso1 = 6
!          iso2 = 8
!          DO ia = 1,na
!            rx = dslm_dxyz(2,ia,iso1)*slm(ia,iso2)+slm(ia,iso1)*dslm_dxyz(2,ia,iso2)
!             rx = slm(ia,iso1)*slm(ia,iso2)
!            WRITE(*,*) ia,rx,ss(ia,iso1,iso2)
!          END DO
!  DEALLOCATE(ss)
!  stop 'ss'


    ELSE
      CALL stop_program(routine_name,module_name,__LINE__,&
               "The pointer harmonics is not associated")
    ENDIF

  END SUBROUTINE  create_harmonics_atom
  
!*******************************************************************************
END MODULE qs_harmonics_atom
