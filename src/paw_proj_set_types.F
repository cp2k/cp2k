!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/paw_proj_set_types [1.0] *
!!
!!   NAME
!!   paw_proj_set_types 
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

MODULE paw_proj_set_types

! *****************************************************************************

  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: gapw_control_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE mathconstants,                   ONLY: dfac,pi,rootpi,oorootpi
  USE mathlib,                         ONLY: invert_matrix
  USE memory_utilities,                ONLY: reallocate
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE orbital_pointers,                ONLY: coset,&
                                             indco,nco,ncoset,nso,&
                                             nsoset, indso
  USE qs_util,                         ONLY: gauss_exponent,&
                                             exp_radius
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
 
! *****************************************************************************

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "paw_proj_set_types"

  INTEGER, PARAMETER :: max_name_length = 60

! *** Define the projectors types ***

  TYPE paw_proj_set_type
    INTEGER                                  :: maxl,ncgauprj,nsgauprj
    INTEGER, DIMENSION(:), POINTER           :: nprj       ! 0:maxl
    INTEGER, DIMENSION(:), POINTER           :: lx,ly,lz   ! ncgauprj
    INTEGER, DIMENSION(:), POINTER           :: ll, m      ! nsgauprj
    INTEGER, DIMENSION(:), POINTER           :: first_prj,last_prj ! 0:maxl 
    INTEGER, DIMENSION(:), POINTER           :: first_prjs ! 0:maxl 
    REAL(dp)                                 :: rcprj
    REAL(dp), DIMENSION(:), POINTER          :: zisomin,zprjisomin
    REAL(dp), DIMENSION(:,:), POINTER        :: zetprj     ! maxnprj,0:maxl
    REAL(dp), DIMENSION(:,:), POINTER        :: rzetprj    ! maxnprj,0:maxl
    REAL(dp), DIMENSION(:,:), POINTER        :: cprj       ! ncgauprj, maxco*nset
    REAL(dp), DIMENSION(:,:), POINTER        :: cprj_s     ! nsgauprj, maxso*nset
    REAL(dp), DIMENSION(:,:), POINTER        :: local_oce_cphi_h,local_oce_cphi_s ! maxco,ncgf
    REAL(dp), DIMENSION(:,:), POINTER        :: local_oce_sphi_h,local_oce_sphi_s ! maxco,nsgf
    REAL(dp), DIMENSION(:,:,:,:), POINTER    :: gccprj     ! maxnprj,maxnpgf,0:maxl,nset
    LOGICAL, DIMENSION(:,:), POINTER         :: isoprj     ! maxnprj,0:maxl

  END TYPE paw_proj_set_type



! *** Public subroutines ***

  PUBLIC ::  allocate_paw_proj_set,&
             deallocate_paw_proj_set,&
             get_paw_proj_set,&
             projectors,&
             set_paw_proj_set

! *** Public data types ***

  PUBLIC ::  paw_proj_set_type

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_paw_proj_set(paw_proj_set)

!  Purpose: Allocate projector type for GAPW 

! *****************************************************************************

   TYPE(paw_proj_set_type), POINTER :: paw_proj_set

    CHARACTER(LEN=*), PARAMETER :: routine =" allocate_paw_proj_set"

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(paw_proj_set)) CALL deallocate_paw_proj_set(paw_proj_set)

    ALLOCATE (paw_proj_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,"paw_proj_set",0)

    NULLIFY (paw_proj_set%nprj)
    NULLIFY (paw_proj_set%lx)
    NULLIFY (paw_proj_set%ly)
    NULLIFY (paw_proj_set%lz)
    NULLIFY (paw_proj_set%ll)
    NULLIFY (paw_proj_set%m)
    NULLIFY (paw_proj_set%first_prj)
    NULLIFY (paw_proj_set%last_prj)
    NULLIFY (paw_proj_set%first_prjs)

    NULLIFY (paw_proj_set%zisomin) 
    NULLIFY (paw_proj_set%zprjisomin)
    NULLIFY (paw_proj_set%zetprj)
    NULLIFY (paw_proj_set%cprj)
    NULLIFY (paw_proj_set%cprj_s)
    NULLIFY (paw_proj_set%local_oce_cphi_h)
    NULLIFY (paw_proj_set%local_oce_cphi_s)
    NULLIFY (paw_proj_set%local_oce_sphi_h)
    NULLIFY (paw_proj_set%local_oce_sphi_s)
    NULLIFY (paw_proj_set%gccprj)
    NULLIFY (paw_proj_set%rzetprj)

    NULLIFY (paw_proj_set%isoprj)

  END SUBROUTINE allocate_paw_proj_set

! *****************************************************************************

  SUBROUTINE deallocate_paw_proj_set(paw_proj_set)

!   Purpose: Deallocate a projector-type set data set.

!   ***************************************************************************

    TYPE(paw_proj_set_type), POINTER :: paw_proj_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine = " deallocate_paw_proj_set "

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    DEALLOCATE (paw_proj_set%zisomin,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%zisomin")

    DEALLOCATE (paw_proj_set%zprjisomin,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%zprjisomin")

    DEALLOCATE (paw_proj_set%nprj,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%nprj")

    DEALLOCATE (paw_proj_set%lx,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%lx")

    DEALLOCATE (paw_proj_set%ly,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%ly")

    DEALLOCATE (paw_proj_set%lz,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%lz")

    DEALLOCATE (paw_proj_set%first_prj,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%first_prj")

    DEALLOCATE (paw_proj_set%last_prj,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%last_prj")

    DEALLOCATE (paw_proj_set%first_prjs,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%first_prjs")

    DEALLOCATE (paw_proj_set%zetprj,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%zetprj")

    DEALLOCATE (paw_proj_set%cprj,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%cprj")

    DEALLOCATE (paw_proj_set%cprj_s,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%cprj_s")

    DEALLOCATE (paw_proj_set%local_oce_cphi_h,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%local_oce_cphi_h")

    DEALLOCATE (paw_proj_set%local_oce_cphi_s,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%local_oce_cphi_s")

    DEALLOCATE (paw_proj_set%local_oce_sphi_h,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%local_oce_sphi_h")

    DEALLOCATE (paw_proj_set%local_oce_sphi_s,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%local_oce_sphi_s")

    DEALLOCATE (paw_proj_set%gccprj,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%gccprj")
 
    DEALLOCATE (paw_proj_set%isoprj,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%isoprj")
 
    DEALLOCATE (paw_proj_set%rzetprj,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set%rzetprj")

    DEALLOCATE (paw_proj_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,module_name,__LINE__,&
                                     "paw_proj_set")

  END SUBROUTINE deallocate_paw_proj_set

! *****************************************************************************

  SUBROUTINE projectors(paw_proj,orb_basis,rc,gapw,globenv)

!   Purpose: initialize the projector-type set data set.

!   ***************************************************************************

    TYPE(paw_proj_set_type),  POINTER         :: paw_proj
    TYPE(gto_basis_set_type), POINTER         :: orb_basis
    TYPE(gapw_control_type),  POINTER         :: gapw
    TYPE(global_environment_type), INTENT(IN) :: globenv
    REAL(dp)                                  :: rc

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine = " projectors "

!   *** Local variables ***

    REAL(dp)                                  :: eps_fit,eps_iso

!   ***************************************************************************

    eps_fit=gapw%eps_fit
    eps_iso=gapw%eps_iso

    CALL build_projector(paw_proj,orb_basis,eps_fit,eps_iso,rc,globenv)

  END SUBROUTINE  projectors

! *****************************************************************************

  SUBROUTINE build_projector(paw_proj,orb_basis,&
                             eps_fit,eps_iso,rc,globenv)

!   Purpose: initialize the projector-type set data set.

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER :: orb_basis
    TYPE(paw_proj_set_type),  POINTER :: paw_proj
    REAL(dp), INTENT(IN)              :: eps_fit,eps_iso,rc
    TYPE(global_environment_type), INTENT(IN) :: globenv


    CHARACTER(LEN=*), PARAMETER :: routine = "  build_projectors "

    CHARACTER(LEN=default_string_length)   :: bsname
    INTEGER  :: icgf, icgfmin, icgfmax, ico, icomin, icomax, il, ip, ipgf, iprj, &
      iprjfirst, isgfmin, iset, ishell, jp, kp, lshell, m, maxco, maxl, maxnprj, maxpgf, & 
      maxso, iso, is1, is2, lx, ly, lz,  &
      n, nc, ncgauprj, ncgf, nisop, np, npgfg, ns, nset, nsgf, output_unit
    INTEGER, DIMENSION(:),    POINTER :: lmin,lmax,npgf,nshell
    INTEGER, DIMENSION(:,:),  POINTER :: first_cgf,first_sgf,last_cgf,l
    LOGICAL,  DIMENSION(:),   POINTER :: isoprj
    REAL(dp)                          :: expzet, error, prefac, radius,&
                                         x, zetmin, zetval
    REAL(dp), DIMENSION(:),   POINTER :: zet,zetp
    REAL(dp), DIMENSION(:,:), POINTER :: gcc,smat,work
    REAL(dp), DIMENSION(:,:), POINTER :: cphi

    INTEGER :: mc, ms, ic, iprjs, nsgauprj
    REAL(dp), DIMENSION(:,:), POINTER ::  sphi_s,sphi_h

!   ***************************************************************************

    NULLIFY (first_cgf,first_sgf,last_cgf,isoprj,gcc,l)
    NULLIFY (cphi,lmax,lmin,npgf,nshell,zet,zetp,smat,work)
    NULLIFY (sphi_s,sphi_h)

!   ***************************************************************************

    CALL get_gto_basis_set(gto_basis_set=orb_basis,name=bsname,&
                           ncgf=ncgf,nset=nset,nsgf=nsgf,&
                           lmax=lmax,lmin=lmin,npgf=npgf,&
                           nshell=nshell,cphi=cphi,&
                           first_cgf=first_cgf,first_sgf=first_sgf,&
                           l=l,last_cgf=last_cgf,&
                           maxco=maxco,maxso=maxso,maxl=maxl,maxpgf=maxpgf)

    paw_proj%maxl = maxl
    CALL reallocate(paw_proj%zisomin,0,maxl)
    CALL reallocate(paw_proj%zprjisomin,0,maxl)
    CALL reallocate(paw_proj%nprj,0,maxl)

    IF (globenv%ionode .AND. globenv%print%basis_set) THEN
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
             "Projectors for the basis functions of"//TRIM(bsname)
    ENDIF

!   *** Determine the max # of projectors for an ang. mom. q.n. ***
    maxnprj = 0
    DO lshell = 0,maxl     ! lshell
      np = 0
      DO iset = 1,nset
        IF (lshell >= lmin(iset) .AND. lshell <= lmax(iset))THEN
          np = np + npgf(iset)
        ENDIF
      ENDDO
      maxnprj = MAX(maxnprj,np)
      
      paw_proj%nprj(lshell) = np 

    ENDDO                  ! lshell


!   *** Allocate exponents and coefficients ***
    CALL reallocate(paw_proj%zetprj,1,maxnprj,0,maxl)
    CALL reallocate(paw_proj%gccprj,1,maxnprj,1,maxpgf,0,maxl,1,nset)
    CALL reallocate(paw_proj%rzetprj,1,maxnprj,0,maxl)
    allocate(paw_proj%isoprj(maxnprj,0:maxl))
    paw_proj%isoprj = .FALSE.


!   *** Generate the projetor basis for each ang. mom. q.n. ***
    DO lshell = 0,maxl    ! lshell

      np = paw_proj%nprj(lshell) 

      allocate(isoprj(np))
      isoprj = .FALSE.
      CALL reallocate(zet,1,np)
      CALL reallocate(zetp,1,np)
      CALL reallocate(gcc,1,np,1,np)
      CALL reallocate(smat,1,np,1,np)
      CALL reallocate(work,1,np,1,np)

      npgfg = 0
!     *** Collect all the exponent which contribute to lshell ***
      DO iset = 1,nset     ! iset
        IF (lshell >= lmin(iset) .AND. lshell <= lmax(iset))THEN
          zet(npgfg+1:npgfg+npgf(iset)) = orb_basis%zet(1:npgf(iset),iset)
          npgfg = npgfg + npgf(iset)
        ENDIF
      ENDDO                ! iset

!     *** Smallest exp. due to eps_iso: concerned as an isolated projector ***
      paw_proj%zisomin(lshell) = gauss_exponent(lshell,rc,eps_iso,1.0_dp)

      isoprj(1:np) = (zet(1:np) >= paw_proj%zisomin(lshell))

!     *** Smallest exp. due to eps_fit: where to start geometric progression ***
      zetmin = gauss_exponent(lshell,rc,eps_fit,1.0_dp)

!     *** Generate the projectors by the geometric progression ***
      x = 2.0_dp

      zetval = zetmin
      DO ip = np,1, -1
        IF(.NOT. isoprj(ip)) THEN
          zetp(ip) = zetval
          zetval = x*zetval
        ENDIF
      ENDDO

      paw_proj%zprjisomin(lshell) = zetval

      nisop = 0
      DO ip = np,1, -1
        IF(isoprj(ip)) THEN
          zetp(ip) = zetval
          zetval = x*zetval
          nisop = nisop + 1
        ENDIF
      ENDDO


!     *** Build the overlap matrix: <projector|primitive> ***
      prefac = 0.5_dp**(lshell + 2)*rootpi*dfac(2*lshell + 1)
      expzet = REAL(lshell,dp) + 1.5_dp

      DO ip = 1,np
        IF (isoprj(ip)) THEN
          smat(ip,ip) = prefac/(zetp(ip) + zet(ip))**expzet
        ELSE
          DO jp = 1,np
            IF (.not. isoprj(jp)) THEN
              smat(ip,jp) = prefac/(zetp(ip) + zet(jp))**expzet
            ENDIF
           ENDDO
        ENDIF
      ENDDO

!     *** Compute inverse of the transpose ***
      CALL invert_matrix(smat,gcc,error,"T")

!     *** Compute the orthogonality error ***

!     *** Set the coefficient of the isolated projectors to 0 ***
      DO ip = 1, np
        IF(isoprj(ip)) THEN
          gcc(ip,ip) = 0.0_dp
        ENDIF
      ENDDO

!     *** Transfer data from local to global variables ***

      paw_proj%zetprj(1:np,lshell) = zetp(1:np)
      paw_proj%isoprj(1:np,lshell) = isoprj(1:np)

      npgfg = 0
      DO iset = 1,nset     ! iset
        IF (lshell >= lmin(iset) .AND. lshell <= lmax(iset))THEN
          paw_proj%gccprj(1:np,1:npgf(iset),lshell,iset) = &
               gcc(1:np,npgfg+1:npgfg+npgf(iset))
          npgfg = npgfg + npgf(iset)
        ENDIF
      ENDDO                ! iset

!     *** Print exponents and coefficients of the projectors ***
      IF (globenv%ionode .AND. globenv%print%basis_set) THEN
        output_unit = globenv%scr

        WRITE (UNIT=output_unit,FMT="(/,/,T2,A,I2)")&
         "Built projector for angular momentum quantum number l= ", lshell
        WRITE (UNIT=output_unit,FMT="(T2,A,I2)")&
         "Number of isolated projectors = ", nisop
        DO iset = 1,nset     ! iset
          IF (lshell >= lmin(iset) .AND. lshell <= lmax(iset))THEN
            WRITE (UNIT=output_unit,FMT="(/,T2,A,I5,/,/,T4,A5,(T9,4f15.6))")&
                   "Set ",iset,"exp prj: ", &
                   (paw_proj%zetprj(ip,lshell),ip=1,np)
            DO jp = 1,npgf(iset)
              WRITE (UNIT=output_unit,FMT="(/,T4,A5,F15.6,/,T4,A5,(t9,4E15.6))") &
                "exp gto: ",orb_basis%zet(jp,iset),&
                "coeff.:  ",(paw_proj%gccprj(ip,jp,&
                             lshell,iset),ip=1,np)
            ENDDO
          ENDIF
        ENDDO                ! iset


      ENDIF

!     *** Release the working storage for the current value lshell ***
      NULLIFY(isoprj,gcc,zet,zetp,smat,work)

    ENDDO                 ! lshell

!   *** Count primitives basis functions for the projectors
    paw_proj%ncgauprj = 0
    paw_proj%nsgauprj = 0
    DO lshell = 0,maxl
      paw_proj%ncgauprj = paw_proj%ncgauprj +nco(lshell)*paw_proj%nprj(lshell)
      paw_proj%nsgauprj = paw_proj%nsgauprj +nso(lshell)*paw_proj%nprj(lshell)
    ENDDO

    ncgauprj = paw_proj%ncgauprj
    nsgauprj = paw_proj%nsgauprj
    CALL reallocate(paw_proj%cprj,1,ncgauprj,1,maxso*nset)
    CALL reallocate(paw_proj%lx,1,ncgauprj)
    CALL reallocate(paw_proj%ly,1,ncgauprj)
    CALL reallocate(paw_proj%lz,1,ncgauprj)
    CALL reallocate(paw_proj%first_prj,0,maxl)
    CALL reallocate(paw_proj%last_prj, 0,maxl)
    CALL reallocate(paw_proj%cprj_s,1,nsgauprj,1,maxso*nset)
    CALL reallocate(paw_proj%ll,1,nsgauprj)
    CALL reallocate(paw_proj%m,1,nsgauprj)
    CALL reallocate(paw_proj%first_prjs,0,maxl)

!    NULLIFY (matl)
!    CALL reallocate(matl,1,ncoset(maxl),1,ncoset(maxl),1,1)

    ncgauprj = 0
    nsgauprj  = 0
    DO lshell = 0,maxl
      np = paw_proj%nprj(lshell)
      paw_proj%first_prj(lshell) = ncgauprj + 1
      paw_proj%first_prjs(lshell) = nsgauprj + 1
      paw_proj%last_prj(lshell) = ncgauprj + nco(lshell)*np
      DO ip = 1,np
        DO ico = ncoset(lshell-1) + 1, ncoset(lshell)
          ncgauprj = ncgauprj + 1
          paw_proj%lx(ncgauprj) = indco(1,ico)
          paw_proj%ly(ncgauprj) = indco(2,ico)
          paw_proj%lz(ncgauprj) = indco(3,ico)
        ENDDO   ! ico
        DO iso = nsoset(lshell - 1) + 1 , nsoset(lshell)
          nsgauprj = nsgauprj + 1
          paw_proj%ll(nsgauprj) = indso(1,iso)
          paw_proj%m(nsgauprj)  = indso(2,iso)
        ENDDO
      ENDDO     ! ip
    ENDDO       ! lshell

    ms = 0
    DO iset=1,orb_basis%nset
      ns = nsoset(lmax(iset)) 
      DO lshell = lmin(iset),lmax(iset)
        iprjfirst =  paw_proj%first_prjs(lshell)
        np = paw_proj%nprj(lshell)
        DO ipgf = 1,npgf(iset)
          DO ip = 1, np
            DO il = 1, nso(lshell) 
              iprjs =  iprjfirst - 1 + il + (ip-1)* nso(lshell)    
              iso =  nsoset(lshell-1)+1+(lshell+paw_proj%m(iprjs))

              iso = iso + (ipgf-1)*ns + ms
              paw_proj%cprj_s(iprjs,iso) = & 
                      paw_proj%gccprj(ip,ipgf,lshell,iset)
            ENDDO        ! iprjs
          ENDDO          ! ip
        ENDDO            ! ipgf
      ENDDO              ! lshell
      ms = ms + maxso
    ENDDO                ! iset


    ms = 0
    DO iset=1,orb_basis%nset
      ns = nsoset(lmax(iset))
      DO ipgf = 1,npgf(iset)
        DO iso = nsoset(lmin(iset)-1)+1+ns*(ipgf-1)+ms,ns*ipgf+ms 
          DO lshell = 0,maxl
            np = paw_proj%nprj(lshell)
            DO ip = 1, np
              jp = paw_proj%first_prj(lshell) + (ip-1)* nco(lshell)-1
              kp = paw_proj%first_prjs(lshell)+ (ip-1)* nso(lshell)-1
              DO il = 1,nco(lshell)
                lx = indco(1,il+ncoset(lshell-1))
                ly = indco(2,il+ncoset(lshell-1))
                lz = indco(3,il+ncoset(lshell-1))
                iprj = jp + il
                DO ic = 1,nso(lshell)
                  iprjs = kp + ic
                  paw_proj%cprj(iprj,iso) = paw_proj%cprj(iprj,iso) +&
                  orbtramat(lshell)%s2c(ic,il)*paw_proj%cprj_s(iprjs,iso)*&
                  sqrt(dfac(2*lshell+1)/(4.0_dp*pi)/&
                  dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)) 
                ENDDO ! ic
              ENDDO ! il
            ENDDO ! ip
          ENDDO  ! lshell
        ENDDO  ! iso
      ENDDO  ! ipgf
      ms = ms + maxso
    ENDDO ! iset

!   local coefficients for the one center expansions : oce
!   the coefficients are calculated for the full and soft expansions
    CALL reallocate(paw_proj%local_oce_cphi_h, 1,maxco,1,ncgf)
    CALL reallocate(paw_proj%local_oce_sphi_h, 1,maxso,1,nsgf)
    CALL reallocate(sphi_h, 1,maxco,1,nsgf)

    CALL reallocate(paw_proj%local_oce_cphi_s, 1,maxco,1,ncgf)
    CALL reallocate(paw_proj%local_oce_sphi_s, 1,maxso,1,nsgf)
    CALL reallocate(sphi_s, 1,maxco,1,nsgf)


!   *** Cartesian ***
    DO iset = 1,nset
      n = ncoset(lmax(iset))
      DO ipgf = 1,npgf(iset)
        DO ishell = 1,nshell(iset)
          lshell = l(ishell,iset)
          icomin = ncoset(lshell-1) + 1 + n*(ipgf - 1)
          icomax = ncoset(lshell) + n*(ipgf - 1)
          icgfmin = first_cgf(ishell,iset)
          icgfmax = last_cgf(ishell,iset)
          radius = exp_radius(lshell,orb_basis%zet(ipgf,iset),&
                              eps_fit,1.0_dp)
          DO icgf = icgfmin,icgfmax 
            paw_proj%local_oce_cphi_h(icomin:icomax,icgf) = &
                              cphi(icomin:icomax,icgf)
            IF (radius < rc ) THEN
              paw_proj%local_oce_cphi_s(icomin:icomax,icgf)=0.0_dp
            ELSE
              paw_proj%local_oce_cphi_s(icomin:icomax,icgf) = &
                              cphi(icomin:icomax,icgf)
            ENDIF
          END DO
        ENDDO        ! ishell
      ENDDO          ! ipgf 
    ENDDO            ! iset


!   *** Spherical***
    DO ic = 1,maxco
    DO iset=1,nset
      DO ishell=1,nshell(iset)
        lshell = l(ishell,iset)
        icgfmin = first_cgf(ishell,iset)
        isgfmin = first_sgf(ishell,iset)
        nc = nco(lshell)
        ns = nso(lshell)
        DO iso = 1,ns
          is1 = isgfmin + iso - 1
          DO ico = 1,nc
            is2 = icgfmin + ico - 1
            lx = indco(1,ico+ncoset(lshell-1))
            ly = indco(2,ico+ncoset(lshell-1))
            lz = indco(3,ico+ncoset(lshell-1))  

            sphi_h(ic,is1) =  sphi_h(ic,is1) + &
               paw_proj%local_oce_cphi_h(ic,is2) * &
               orbtramat(lshell)%c2s(iso,ico) ! * &
!               sqrt(4.0_dp*pi/dfac(2*lshell+1)*&
!              dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1))
 
            sphi_s(ic,is1) =  sphi_s(ic,is1) + &
               paw_proj%local_oce_cphi_s(ic,is2) * &
               orbtramat(lshell)%c2s(iso,ico)  !* &
!               sqrt(4.0_dp*pi/dfac(2*lshell+1)*&
!               dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1))

          ENDDO  ! ico
        ENDDO  ! iso
      END DO  ! ishell
    END DO  ! iset
    ENDDO  ! ic


    DO iset=1,nset
      ns = nsoset(lmax(iset))
      nc = ncoset(lmax(iset))
      isgfmin = first_sgf(1,iset)
      ms = orb_basis%nsgf_set(iset)
      DO is1 = isgfmin,isgfmin+ms-1
        DO ipgf = 1,npgf(iset)
          DO lshell = lmin(iset),lmax(iset)
            iso = nsoset(lshell - 1)  + ns*(ipgf-1) 
            ico = ncoset(lshell - 1)  + nc*(ipgf-1) 
            DO is2 = 1,nso(lshell)
              DO ic = 1,nco(lshell)
                lx = indco(1,ic+ncoset(lshell-1))
                ly = indco(2,ic+ncoset(lshell-1))
                lz = indco(3,ic+ncoset(lshell-1))
                paw_proj%local_oce_sphi_h(iso+is2,is1) = &
                    paw_proj%local_oce_sphi_h(iso+is2,is1) + & 
                    orbtramat(lshell)%c2s(is2,ic)*sphi_h(ico+ic,is1)*&
                    sqrt(4.0_dp*pi/dfac(2*lshell+1)*&
                    dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1))

                paw_proj%local_oce_sphi_s(iso+is2,is1) = &
                    paw_proj%local_oce_sphi_s(iso+is2,is1) + &
                    orbtramat(lshell)%c2s(is2,ic)*sphi_s(ico+ic,is1)*&
                    sqrt(4.0_dp*pi/dfac(2*lshell+1)*&
                    dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1))

              ENDDO  ! ic2
            ENDDO  ! is2
          ENDDO  ! lshell 
        ENDDO  ! ipgf
      ENDDO  ! is1
    ENDDO  ! iset

    NULLIFY(sphi_s,sphi_h)

  END SUBROUTINE  build_projector

! *****************************************************************************

  SUBROUTINE get_paw_proj_set(paw_proj_set,cprj,cprj_s,&
                              first_prj,first_prjs,last_prj, & 
                              local_oce_sphi_h,local_oce_sphi_s,&
                              maxl,ncgauprj,nsgauprj,nprj,&
                              rcprj,rzetprj,zisomin,zetprj)

!   Purpose: Get informations about a paw projectors set.

!   ***************************************************************************

    TYPE(paw_proj_set_type),                POINTER :: paw_proj_set
    INTEGER, OPTIONAL, INTENT(OUT)                  :: maxl, ncgauprj,nsgauprj 
    INTEGER, DIMENSION(:), OPTIONAL, POINTER        :: nprj
    INTEGER, DIMENSION(:), OPTIONAL, POINTER        :: first_prj,first_prjs,last_prj
    REAL(dp), OPTIONAL, INTENT(OUT)                 :: rcprj
    REAL(dp), DIMENSION(:), OPTIONAL, POINTER       :: zisomin
    REAL(dp), DIMENSION(:,:), OPTIONAL, POINTER     :: zetprj 
    REAL(dp), DIMENSION(:,:), OPTIONAL, POINTER     :: rzetprj 
    REAL(dp), DIMENSION(:,:), OPTIONAL, POINTER     :: cprj     
    REAL(dp), DIMENSION(:,:), OPTIONAL, POINTER     :: cprj_s

    REAL(dp), DIMENSION(:,:),OPTIONAL, POINTER     :: local_oce_sphi_h,local_oce_sphi_s 


    CHARACTER(LEN=*), PARAMETER :: routine = " get_paw_proj_set"


!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(paw_proj_set)) THEN
      IF (PRESENT(cprj)) cprj => paw_proj_set%cprj
      IF (PRESENT(cprj_s)) cprj_s => paw_proj_set%cprj_s
      IF (PRESENT(local_oce_sphi_h)) local_oce_sphi_h => paw_proj_set%local_oce_sphi_h
      IF (PRESENT(local_oce_sphi_s)) local_oce_sphi_s => paw_proj_set%local_oce_sphi_s
      IF (PRESENT(first_prj)) first_prj => paw_proj_set%first_prj
      IF (PRESENT(last_prj)) last_prj => paw_proj_set%last_prj
      IF (PRESENT(first_prjs)) first_prjs => paw_proj_set%first_prjs
      IF (PRESENT(maxl)) maxl = paw_proj_set%maxl
      IF (PRESENT(ncgauprj )) ncgauprj  = paw_proj_set%ncgauprj
      IF (PRESENT(nsgauprj )) nsgauprj  = paw_proj_set%nsgauprj
      IF (PRESENT(nprj))   nprj => paw_proj_set%nprj
      IF (PRESENT(rcprj))  rcprj = paw_proj_set%rcprj
      IF (PRESENT(rzetprj)) rzetprj => paw_proj_set%rzetprj 
      IF (PRESENT(zisomin)) zisomin => paw_proj_set%zisomin
      IF (PRESENT(zetprj)) zetprj => paw_proj_set%zetprj

    ELSE
      CALL stop_program(routine,module_name,__LINE__,&
                       "The pointer  paw_proj_set is not associated")
    ENDIF

  END SUBROUTINE get_paw_proj_set

! *****************************************************************************

  SUBROUTINE set_paw_proj_set(paw_proj_set, rzetprj, rcprj)

!   Purpose: Set informations about a paw projectors set.

!   ***************************************************************************

    TYPE(paw_proj_set_type),            POINTER :: paw_proj_set
    REAL(dp), DIMENSION(:,:), OPTIONAL, POINTER :: rzetprj 
    REAL(dp), INTENT(IN), OPTIONAL              :: rcprj

    CHARACTER(LEN=*), PARAMETER :: routine = " set_paw_proj_set"

!   ***************************************************************************

    IF (ASSOCIATED(paw_proj_set)) THEN
      IF (PRESENT(rzetprj)) paw_proj_set%rzetprj(:,0:) = rzetprj(:,0:)
      IF (PRESENT(rcprj))  paw_proj_set%rcprj = rcprj
    ELSE
      CALL stop_program(routine,module_name,__LINE__,&
                        "The pointer paw_proj_set is not associated")
    END IF

  END SUBROUTINE set_paw_proj_set
! *****************************************************************************

END MODULE paw_proj_set_types

! *****************************************************************************
