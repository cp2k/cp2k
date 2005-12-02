!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/soft_basis_set [1.0] *
!!
!!   NAME
!!     soft_basis_set
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MI (08.01.2004)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE soft_basis_set

! *****************************************************************************

  USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                             deallocate_gto_basis_set,&
                                             init_cphi_and_sphi,&
                                             get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: gapw_control_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: dfac,pi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: coset,init_orbital_pointers,&
                                             indco,nco,ncoset,nso
  USE orbital_symbols,                 ONLY: cgf_symbol,&
                                             sgf_symbol
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE qs_util,                         ONLY: exp_radius
  USE termination,                     ONLY: stop_memory,&
                                             stop_program


  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "soft_basis_set"

  INTEGER, PARAMETER :: max_name_length = 60

! *** Public subroutines ***

  PUBLIC :: create_soft_basis

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE create_soft_basis(orb_basis,soft_basis,gapw,rc,paw_atom)
   
!   Purpose: create the soft basis from a GTO basis 

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER   :: orb_basis,soft_basis
    TYPE(gapw_control_type),  POINTER   :: gapw
    REAL(dp), INTENT(IN)                :: rc
    LOGICAL, INTENT(INOUT)              :: paw_atom

    CHARACTER(LEN=*), PARAMETER :: routine = "create_soft_basis"

    CHARACTER(LEN=default_string_length)      :: bsname
    INTEGER :: i, ico, ikind, ipgf, iset, ishell, istat, lshell, m, maxco, maxl, &
      maxpgf, maxshell, ncgf, ncgf_set, nkind, nmin, nset, nsgf, &
      ipgf_s, iset_s, lshell_old, maxpgf_s, maxshell_s, nset_s, &
      icgf,jpgf,lx,ly,lz
    INTEGER, DIMENSION(:), POINTER      :: lmax,lmin,npgf,nshell
    INTEGER, DIMENSION(:,:), POINTER    :: l,n
    REAL(dp)                            :: eps_fit, expzet, fnorm, & 
                                           gcca, gccb, minzet, & 
                                           prefac, radius, &
                                           zeta,zetb
    REAL(dp), DIMENSION(:,:), POINTER   :: zet
    REAL(dp), DIMENSION(:,:,:), POINTER :: gcc
    INTEGER, DIMENSION(:),ALLOCATABLE   :: iset_s2h

!   ---------------------------------------------------------------------------

    NULLIFY(gcc,l,lmax,lmin,n,npgf,nshell,zet)

    eps_fit=gapw%eps_fit

    CALL get_gto_basis_set(gto_basis_set=orb_basis,name=bsname,&
                           maxpgf=maxpgf,maxshell=maxshell,nset=nset)

    paw_atom = .FALSE.
    soft_basis%name = TRIM(bsname)//"_soft"

    CALL reallocate(npgf,1,nset)
    CALL reallocate(nshell,1,nset)
    CALL reallocate(lmax,1,nset)
    CALL reallocate(lmin,1,nset)


    CALL reallocate(n,1,maxshell,1,nset)
    CALL reallocate(l,1,maxshell,1,nset)

    CALL reallocate(zet,1,maxpgf,1,nset)
    CALL reallocate(gcc,1,maxpgf,1,maxshell,1,nset)

    ALLOCATE(iset_s2h(nset))
    iset_s2h = 0
    iset_s = 0
    maxpgf_s =0
    maxshell_s = 0

    DO iset=1,nset                            ! iset
       minzet = orb_basis%zet(orb_basis%npgf(iset),iset)
       DO ipgf = orb_basis%npgf(iset)-1,1,-1
         IF (orb_basis%zet(ipgf,iset) < minzet) THEN
            minzet = orb_basis%zet(ipgf,iset)
         ENDIF
       ENDDO
       radius = exp_radius(orb_basis%lmax(iset),&
                           minzet,eps_fit,1.0_dp)

!      The soft basis contains this set
       iset_s = iset_s + 1
       nshell(iset_s) = orb_basis%nshell(iset)
       lmax(iset_s)   = orb_basis%lmax(iset)
       lmin(iset_s)   = orb_basis%lmin(iset)

       iset_s2h(iset_s) = iset

       DO ishell=1,nshell(iset_s)
          n(ishell,iset_s) =  orb_basis%n(ishell,iset)
          l(ishell,iset_s) =  orb_basis%l(ishell,iset)
       ENDDO
    
       IF ( nshell(iset_s) > maxshell_s) THEN
          maxshell_s = nshell(iset_s)
       ENDIF

       if (radius < rc ) then
!      The soft basis does not contain this set
!      For the moment I keep the set as a dummy set
!      with no exponents, in order to have the right number of contractions
!      In a second time it can be taken away, by creating a pointer
!      which connects the remaining sets to the right contraction index
         paw_atom = .true.
         npgf(iset_s) = 0 
         cycle
       endif

       ipgf_s = 0
       DO ipgf = 1,orb_basis%npgf(iset)         ! ipgf 
         IF (orb_basis%zet(ipgf,iset) > 100.0_dp) THEN
!          The soft basis does not contain this exponent
           paw_atom = .true.
           cycle
         ENDIF

         radius = exp_radius(orb_basis%lmax(iset),&
                             orb_basis%zet(ipgf,iset),&
                             eps_fit,1.0_dp)
         IF (radius < rc ) THEN
!           The soft basis does not contain this exponent
            paw_atom = .true.
            cycle
         ENDIF

!        The soft basis contains this exponent
         ipgf_s = ipgf_s + 1
         zet(ipgf_s,iset_s) = orb_basis%zet(ipgf,iset)

         lshell_old = orb_basis%l(1,iset)
         radius = exp_radius(lshell_old,&
                             zet(ipgf_s,iset_s),&
                             eps_fit,1.0_dp)

         DO ishell = 1,nshell(iset_s)           ! ishell
            lshell = orb_basis%l(ishell,iset)      
            IF (lshell == lshell_old) THEN
            ELSE
              lshell_old = lshell
              radius = exp_radius(lshell_old,&
                          zet(ipgf_s,iset_s),&
                               eps_fit,1.0_dp)
            ENDIF
            IF (radius < rc ) THEN
               gcc(ipgf_s,ishell,iset_s) = 0.0_dp
               paw_atom = .true.
            ELSE
               gcc(ipgf_s,ishell,iset_s) = orb_basis%gcc(ipgf,ishell,iset)
            ENDIF
              
          ENDDO                                  ! ishell
       ENDDO                                    ! ipgf
       npgf(iset_s) = ipgf_s
       IF ( ipgf_s > maxpgf_s) THEN
          maxpgf_s = ipgf_s
       ENDIF
    ENDDO                                      ! iset    
    nset_s = iset_s
 
    IF (.not. paw_atom) THEN
!      CALL  deallocate_gto_basis_set(soft_basis)
        DEALLOCATE (soft_basis,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "soft_basis_set")
      soft_basis => orb_basis
    ELSE
       soft_basis%nset = nset_s 
       CALL reallocate(soft_basis%lmax,1,nset_s)
       CALL reallocate(soft_basis%lmin,1,nset_s)
       CALL reallocate(soft_basis%npgf,1,nset_s)
       CALL reallocate(soft_basis%nshell,1,nset_s)
       CALL reallocate(soft_basis%n,1,maxshell_s,1,nset_s)
       CALL reallocate(soft_basis%l,1,maxshell_s,1,nset_s)
       CALL reallocate(soft_basis%zet,1,maxpgf_s,1,nset_s)
       CALL reallocate(soft_basis%gcc,1,maxpgf_s,1,maxshell_s,1,nset_s)

!      *** Copy the basis set information into the data structure ***

       DO iset=1,nset_s
          soft_basis%lmax(iset) = lmax(iset)
          soft_basis%lmin(iset) = lmin(iset)
          soft_basis%npgf(iset) = npgf(iset)
          soft_basis%nshell(iset) = nshell(iset)
          DO ishell=1,nshell(iset)
            soft_basis%n(ishell,iset) = n(ishell,iset)
            soft_basis%l(ishell,iset) = l(ishell,iset)
            DO ipgf=1,npgf(iset)
              soft_basis%gcc(ipgf,ishell,iset) = gcc(ipgf,ishell,iset)
            END DO
          END DO
          DO ipgf=1,npgf(iset)
            soft_basis%zet(ipgf,iset) = zet(ipgf,iset)
          END DO
       END DO

!      *** Initialise the depending soft_basis pointers ***

       CALL reallocate(soft_basis%set_radius,1,nset_s)
       CALL reallocate(soft_basis%pgf_radius,1,maxpgf_s,1,nset_s)
       CALL reallocate(soft_basis%first_cgf,1,maxshell_s,1,nset_s)
       CALL reallocate(soft_basis%first_sgf,1,maxshell_s,1,nset_s)
       CALL reallocate(soft_basis%last_cgf,1,maxshell_s,1,nset_s)
       CALL reallocate(soft_basis%last_sgf,1,maxshell_s,1,nset_s)
       CALL reallocate(soft_basis%ncgf_set,1,nset_s)
       CALL reallocate(soft_basis%nsgf_set,1,nset_s)

       maxco = 0
       ncgf = 0
       nsgf = 0

       DO iset=1,nset_s
          soft_basis%ncgf_set(iset) = 0
          soft_basis%nsgf_set(iset) = 0
          DO ishell=1,nshell(iset)
            lshell = soft_basis%l(ishell,iset)
            soft_basis%first_cgf(ishell,iset) = ncgf + 1
            ncgf = ncgf + nco(lshell)
            soft_basis%last_cgf(ishell,iset) = ncgf
            soft_basis%ncgf_set(iset) =&
            soft_basis%ncgf_set(iset) + nco(lshell)
            soft_basis%first_sgf(ishell,iset) = nsgf + 1
            nsgf = nsgf + nso(lshell)
            soft_basis%last_sgf(ishell,iset) = nsgf
            soft_basis%nsgf_set(iset) =&
            soft_basis%nsgf_set(iset) + nso(lshell)
          END DO
          maxco = MAX(maxco,npgf(iset)*ncoset(lmax(iset)))
       END DO
       soft_basis%ncgf = ncgf
       soft_basis%nsgf = nsgf

       CALL reallocate(soft_basis%cphi,1,maxco,1,ncgf)
       CALL reallocate(soft_basis%sphi,1,maxco,1,nsgf)
       CALL reallocate(soft_basis%lx,1,ncgf)
       CALL reallocate(soft_basis%ly,1,ncgf)
       CALL reallocate(soft_basis%lz,1,ncgf)
       CALL reallocate(soft_basis%m,1,nsgf)
       CALL reallocate(soft_basis%norm_cgf,1,ncgf)
       ALLOCATE (soft_basis%cgf_symbol(ncgf),STAT=istat)
       IF (istat /= 0) THEN
         CALL stop_memory(routine,module_name,__LINE__,&
                         "soft_basis_set%cgf_symbol",12*ncgf)
       END IF
       ALLOCATE (soft_basis%sgf_symbol(nsgf),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routine,module_name,__LINE__,&
                          "soft_basis_set%sgf_symbol",6*nsgf)
       END IF

       ncgf = 0
       nsgf = 0

       DO iset=1,nset_s
          DO ishell=1,nshell(iset)
            lshell = soft_basis%l(ishell,iset)
            DO ico=ncoset(lshell-1)+1,ncoset(lshell)
              ncgf = ncgf + 1
              soft_basis%lx(ncgf) = indco(1,ico)
              soft_basis%ly(ncgf) = indco(2,ico)
              soft_basis%lz(ncgf) = indco(3,ico)
              soft_basis%cgf_symbol(ncgf) =&
              cgf_symbol(n(ishell,iset),(/soft_basis%lx(ncgf),&
                                         soft_basis%ly(ncgf),&
                                         soft_basis%lz(ncgf)/))
            END DO
            DO m=-lshell,lshell
              nsgf = nsgf + 1
              soft_basis%m(nsgf) = m
              soft_basis%sgf_symbol(nsgf) =&
              sgf_symbol(n(ishell,iset),lshell,m)
            END DO
          END DO
       ENDDO

!      *** Normalization factor of the contracted Gaussians *** 
       soft_basis%norm_type = orb_basis%norm_type

       SELECT CASE (soft_basis%norm_type)
       CASE ( 0 )
         !   *** No normalisation requested ***
       CASE ( 1 ) 
         DO iset_s = 1 , nset_s            ! iset_s
            iset = iset_s2h(iset_s)
            DO ishell = 1,nshell(iset_s)
              DO icgf=soft_basis%first_cgf(ishell,iset_s),&
                soft_basis%last_cgf(ishell,iset_s)
                soft_basis%norm_cgf(icgf) = 1.0_dp
              END DO
            END DO                          ! ishell
         END DO
       CASE ( 2 )
         DO iset_s = 1 , nset_s            ! iset_s
            iset = iset_s2h(iset_s)
            DO ishell = 1,nshell(iset_s)    ! ishell
              lshell = l(ishell,iset_s)
              expzet = 0.5_dp*REAL(2*lshell+3,dp)
              fnorm  = 0.0_dp
              DO ipgf=1,orb_basis%npgf(iset)
                gcca = orb_basis%gcc(ipgf,ishell,iset)
                zeta = orb_basis%zet(ipgf,iset)
                DO jpgf=1,orb_basis%npgf(iset)
                  gccb = orb_basis%gcc(jpgf,ishell,iset)
                  zetb = orb_basis%zet(jpgf,iset)
                  fnorm = fnorm + gcca*gccb/(zeta + zetb)**expzet
                END DO
              END DO
              fnorm = 0.5_dp**lshell*pi**1.5_dp*fnorm
              DO icgf=soft_basis%first_cgf(ishell,iset_s),&
                soft_basis%last_cgf(ishell,iset_s)
                lx = soft_basis%lx(icgf)
                ly = soft_basis%ly(icgf)
                lz = soft_basis%lz(icgf)
                prefac = dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1)
                soft_basis%norm_cgf(icgf) = 1.0_dp/SQRT(prefac*fnorm)
              END DO
            END DO                          ! ishell
         END DO                            ! iset_s
       END SELECT 
!      *** Initialize the transformation matrices ***
       CALL init_cphi_and_sphi(soft_basis)          
     ENDIF 

     DEALLOCATE (gcc,l,lmax,lmin,n,npgf,nshell,zet,iset_s2h,STAT=istat)
     IF (istat /= 0) THEN
        CALL stop_memory(routine,module_name,__LINE__,&
                        "gcc,l,lmax,lmin,n,npgf,nshell,zet")
     END IF
     
  END SUBROUTINE create_soft_basis

! *****************************************************************************

END MODULE soft_basis_set
