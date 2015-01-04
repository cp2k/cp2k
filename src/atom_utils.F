!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Some basic routines for atomic calculations
!> \author  jgh
!> \date    01.04.2008
!> \version 1.0
!>
! *****************************************************************************
MODULE atom_utils
  USE atom_types,                      ONLY: CGTO_BASIS,&
                                             GTO_BASIS,&
                                             atom_basis_type,&
                                             atom_gthpot_type,&
                                             atom_state,&
                                             atom_type,&
                                             eri
  USE basis_set_types,                 ONLY: srules
  USE cp_files,                        ONLY: close_file,&
                                             get_unit_number,&
                                             open_file
  USE erf_fn,                          ONLY: erf
  USE input_constants,                 ONLY: do_rhf_atom,&
                                             do_rks_atom,&
                                             do_rohf_atom,&
                                             do_uhf_atom,&
                                             do_uks_atom
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE lapack
  USE mathconstants,                   ONLY: dfac,&
                                             fac,&
                                             fourpi,&
                                             maxfac,&
                                             rootpi
  USE periodic_table,                  ONLY: nelem,&
                                             ptable
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE string_utilities,                ONLY: uppercase
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_utils'

  PUBLIC :: atom_set_occupation, get_maxl_occ, get_maxn_occ
  PUBLIC :: atom_denmat, atom_density, atom_core_density
  PUBLIC :: integrate_grid, atom_trace, atom_solve
  PUBLIC :: coulomb_potential_numeric, coulomb_potential_analytic
  PUBLIC :: exchange_numeric, exchange_semi_analytic
  PUBLIC :: numpot_matrix, ceri_contract, eeri_contract, err_matrix
  PUBLIC :: slater_density, wigner_slater_functional, atom_local_potential
  PUBLIC :: atom_orbital_charge, atom_orbital_nodes, atom_consistent_method
  PUBLIC :: atom_orbital_max, atom_wfnr0
! ZMP added public subroutines
  PUBLIC :: atom_read_external_density
  PUBLIC :: atom_read_external_vxc
  PUBLIC :: atom_read_zmp_restart
  PUBLIC :: atom_write_zmp_restart

!-----------------------------------------------------------------------------!

  INTERFACE integrate_grid
     MODULE PROCEDURE integrate_grid_function1, &
                      integrate_grid_function2, &
                      integrate_grid_function3
  END INTERFACE


! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param ostring ...
!> \param occupation ...
!> \param wfnocc ...
!> \param multiplicity ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_set_occupation(ostring,occupation,wfnocc,multiplicity,error)
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: ostring
    REAL(Kind=dp), DIMENSION(0:3, 10)        :: occupation, wfnocc
    INTEGER, INTENT(OUT), OPTIONAL           :: multiplicity
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_set_occupation', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2)                         :: elem
    CHARACTER(LEN=default_string_length)     :: pstring
    INTEGER                                  :: i, i1, i2, ielem, is, jd, jf, &
                                                jp, js, k, l, mult, n, no
    LOGICAL                                  :: failure
    REAL(Kind=dp)                            :: e0, el, oo

    failure = .FALSE.

    occupation = 0._dp

    CPPrecondition(ASSOCIATED(ostring),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ostring) > 0,cp_failure_level,routineP,error,failure)

    no = SIZE(ostring)

    is = 1
    ! look for multiplicity
    mult = -1 !not specified
    IF ( is <= no ) THEN
      IF ( INDEX(ostring(is),"(") /= 0 ) THEN
        i1 = INDEX(ostring(is),"(")
        i2 = INDEX(ostring(is),")")
        CPAssert((i2-i1-1 > 0).AND.(i2-i1-1 < 3),cp_failure_level,routineP,error,failure)
        elem = ostring(is)(i1+1:i2-1)
        IF ( INDEX(elem,"HS") /= 0 ) THEN
          mult=-2 !High spin
        ELSE IF ( INDEX(elem,"LS") /= 0 ) THEN
          mult=-3 !Low spin
        ELSE
          READ(elem,*) mult
        END IF
        is=is+1
      END IF
    END IF

    IF ( is <= no ) THEN
      IF ( INDEX(ostring(is),"CORE") /= 0 ) is = is+1   !Pseudopotential detected
    END IF
    IF ( is <= no ) THEN
      IF ( INDEX(ostring(is),"none") /= 0 ) is = is+1   !no electrons, used with CORE
    END IF
    IF ( is <= no ) THEN
      IF ( INDEX(ostring(is),"[") /= 0 ) THEN
        ! core occupation from element [XX]
        i1 = INDEX(ostring(is),"[")
        i2 = INDEX(ostring(is),"]")
        CPAssert((i2-i1-1 > 0).AND.(i2-i1-1 < 3),cp_failure_level,routineP,error,failure)
        elem = ostring(is)(i1+1:i2-1)
        ielem=0
        DO k=1,nelem
          IF(elem==ptable(k)%symbol) THEN
            ielem=k
            EXIT
          END IF
        END DO
        CPAssert(ielem /= 0,cp_failure_level,routineP,error,failure)
          DO l=0,3
          el = 2._dp * (2._dp*REAL(l,dp) + 1._dp)
          e0 = ptable(ielem)%e_conv(l)
          DO k=1,10
            occupation(l,k) = MIN(el,e0)
            e0 = e0 - el
            IF ( e0 <= 0._dp ) EXIT
          END DO
        END DO
        is = is+1
      END IF

    END IF

    DO i=is,no
      pstring = ostring(i)
      CALL uppercase(pstring)
      js = INDEX(pstring,"S")
      jp = INDEX(pstring,"P")
      jd = INDEX(pstring,"D")
      jf = INDEX(pstring,"F")
      CPAssert(js+jp+jd+jf > 0,cp_failure_level,routineP,error,failure)
      IF ( js > 0 ) THEN
        CPAssert(jp+jd+jf == 0,cp_failure_level,routineP,error,failure)
        READ(pstring(1:js-1),*) n
        READ(pstring(js+1:),*) oo
        CPAssert(n > 0,cp_failure_level,routineP,error,failure)
        CPAssert(oo >= 0._dp,cp_failure_level,routineP,error,failure)
        CPAssert(occupation(0,n)==0,cp_failure_level,routineP,error,failure)
        occupation(0,n) = oo
      END IF
      IF ( jp > 0 ) THEN
        CPAssert(js+jd+jf == 0,cp_failure_level,routineP,error,failure)
        READ(pstring(1:jp-1),*) n
        READ(pstring(jp+1:),*) oo
        CPAssert(n > 1,cp_failure_level,routineP,error,failure)
        CPAssert(oo >= 0._dp,cp_failure_level,routineP,error,failure)
        CPAssert(occupation(1,n-1)==0,cp_failure_level,routineP,error,failure)
        occupation(1,n-1) = oo
      END IF
      IF ( jd > 0 ) THEN
        CPAssert(js+jp+jf == 0,cp_failure_level,routineP,error,failure)
        READ(pstring(1:jd-1),*) n
        READ(pstring(jd+1:),*) oo
        CPAssert(n > 2,cp_failure_level,routineP,error,failure)
        CPAssert(oo >= 0._dp,cp_failure_level,routineP,error,failure)
        CPAssert(occupation(2,n-2)==0,cp_failure_level,routineP,error,failure)
        occupation(2,n-2) = oo
      END IF
      IF ( jf > 0 ) THEN
        CPAssert(js+jp+jd == 0,cp_failure_level,routineP,error,failure)
        READ(pstring(1:jf-1),*) n
        READ(pstring(jf+1:),*) oo
        CPAssert(n > 3,cp_failure_level,routineP,error,failure)
        CPAssert(oo >= 0._dp,cp_failure_level,routineP,error,failure)
        CPAssert(occupation(3,n-3)==0,cp_failure_level,routineP,error,failure)
        occupation(3,n-3) = oo
      END IF

    END DO

    wfnocc = 0._dp
    DO l=0,3
      k = 0
      DO i=1,10
        IF ( occupation(l,i) /= 0._dp ) THEN
          k = k + 1
          wfnocc(l,k) = occupation(l,i)
        END IF
      END DO
    END DO

    !Check for consistency with multiplicity
    IF (mult /= -1) THEN
      ! count open shells
      js=0
      DO l=0,3
        k=2*(2*l+1)
        DO i=1,10
          IF(wfnocc(l,i) /= 0._dp .AND. wfnocc(l,i) /= REAL(k,dp)) THEN
            js = js+1
            i1 = l
            i2 = i
          END IF
        END DO
      END DO

      IF (js==0 .AND. mult==-2) mult=1
      IF (js==0 .AND. mult==-3) mult=1
      IF (js==0) THEN
        CPAssert(mult==1,cp_failure_level,routineP,error,failure)
      END IF
      IF (js==1) THEN
        l=i1
        i=i2
        k = NINT(wfnocc(l,i))
        IF(k > (2*l+1)) k=2*(2*l+1)-k
        IF(mult==-2) mult=k+1
        IF(mult==-3) mult=MOD(k,2)+1
        CPAssert(MOD(k+1-mult,2)==0,cp_failure_level,routineP,error,failure)
      END IF
      IF (js>1 .AND. mult/=-2) THEN
        CPAssert(mult==-2,cp_failure_level,routineP,error,failure)
      END IF

    END IF

    IF(PRESENT(multiplicity)) multiplicity = mult

  END SUBROUTINE atom_set_occupation

! *****************************************************************************
!> \brief ...
!> \param occupation ...
!> \retval maxl ...
! *****************************************************************************
  FUNCTION get_maxl_occ(occupation) RESULT(maxl)
    REAL(Kind=dp), DIMENSION(0:3, 10)        :: occupation
    INTEGER                                  :: maxl

    INTEGER                                  :: l

    maxl = 0
    DO l=0,3
      IF(SUM(occupation(l,:)) /= 0._dp) maxl=l
    END DO

  END FUNCTION get_maxl_occ

! *****************************************************************************
!> \brief ...
!> \param occupation ...
!> \retval maxn ...
! *****************************************************************************
  FUNCTION get_maxn_occ(occupation) RESULT(maxn)
    REAL(Kind=dp), DIMENSION(0:3, 10)        :: occupation
    INTEGER, DIMENSION(0:3)                  :: maxn

    INTEGER                                  :: k, l

    maxn = 0
    DO l=0,3
      DO k=1,10
        IF(occupation(l,k) /= 0._dp) maxn(l)=maxn(l)+1
      END DO
    END DO

  END FUNCTION get_maxn_occ

! *****************************************************************************
!> \brief ...
!> \param pmat ...
!> \param wfn ...
!> \param nbas ...
!> \param occ ...
!> \param maxl ...
!> \param maxn ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_denmat(pmat,wfn,nbas,occ,maxl,maxn,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: pmat
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: wfn
    INTEGER, DIMENSION(0:3), INTENT(IN)      :: nbas
    REAL(KIND=dp), DIMENSION(0:, :), &
      INTENT(IN)                             :: occ
    INTEGER, INTENT(IN)                      :: maxl
    INTEGER, DIMENSION(0:3), INTENT(IN)      :: maxn
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_denmat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, l, n

    pmat = 0._dp
    n = SIZE(wfn,2)
    DO l=0,maxl
      DO i=1,MIN(n,maxn(l))
        DO k=1,nbas(l)
          DO j=1,nbas(l)
            pmat(j,k,l) =  pmat(j,k,l) + occ(l,i)*wfn(j,i,l)*wfn(k,i,l)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE atom_denmat
! *****************************************************************************
!> \brief ...
!> \param density ...
!> \param pmat ...
!> \param basis ...
!> \param maxl ...
!> \param typ ...
!> \param rr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_density(density,pmat,basis,maxl,typ,rr,error)
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: density
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: pmat
    TYPE(atom_basis_type), INTENT(IN)        :: basis
    INTEGER, INTENT(IN)                      :: maxl
    CHARACTER(LEN=*), OPTIONAL               :: typ
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: rr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_density', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3)                         :: my_typ
    INTEGER                                  :: i, j, l, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ff

    my_typ = "RHO"
    failure = .FALSE.
    IF(PRESENT(typ)) my_typ = typ(1:3)
    IF(my_typ == "KIN") THEN
      CPPrecondition(PRESENT(rr),cp_failure_level,routineP,error,failure)
    END IF

    density=0._dp
    DO l=0,maxl
      n=basis%nbas(l)
      DO i=1,n
        DO j=i,n
          ff=pmat(i,j,l)
          IF(i/=j) ff=2._dp*pmat(i,j,l)
          IF(my_typ == "RHO") THEN
            density(:) = density(:) + ff*basis%bf(:,i,l)*basis%bf(:,j,l)
          ELSE IF(my_typ == "DER") THEN
            density(:) = density(:) + ff*basis%dbf(:,i,l)*basis%bf(:,j,l) &
                                    + ff*basis%bf(:,i,l)*basis%dbf(:,j,l)
          ELSE IF(my_typ == "KIN") THEN
            density(:) = density(:) + 0.5_dp*ff*(&
              basis%dbf(:,i,l)*basis%dbf(:,j,l) + &
              REAL(l*(l+1),dp)*basis%bf(:,i,l)*basis%bf(:,j,l)/rr(:))
          ELSE
            CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
        END DO
      END DO
    END DO
    ! this factor from the product of two spherical harmonics
    density=density/fourpi

  END SUBROUTINE atom_density
! *****************************************************************************
!> \brief   ZMP subroutine to write external restart file 
!>          
!> \param atom ...
!> \date    07.10.2013
!> \author  D. Varsano [daniele.varsano@nano.cnr.it]
!> \version 1.0
! *****************************************************************************
  SUBROUTINE atom_write_zmp_restart(atom)
    TYPE(atom_type), INTENT(IN)              :: atom

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_write_zmp_restart', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: extunit, i, k, l, n

    extunit=get_unit_number()
    CALL open_file(file_name=atom%zmp_restart_file,file_status="UNKNOWN",&
                   file_form="FORMATTED", file_action="WRITE",&
                   unit_number=extunit)

    n = SIZE(atom%orbitals%wfn,2)
    WRITE(extunit,*) atom%basis%nbas
    DO l=0,atom%state%maxl_occ
      DO i=1,MIN(n,atom%state%maxn_occ(l))
        DO k=1,atom%basis%nbas(l)
          WRITE(extunit,*) atom%orbitals%wfn(k,i,l)
        END DO
      END DO
    END DO

    CALL close_file(unit_number=extunit)

  END SUBROUTINE atom_write_zmp_restart

! *****************************************************************************
!> \brief   ZMP subroutine to read external restart file 
!>          
!> \param atom ...
!> \param doguess ...
!> \param iw ...
!> \date    07.10.2013
!> \author  D. Varsano [daniele.varsano@nano.cnr.it]
!> \version 1.0
! *****************************************************************************
  SUBROUTINE atom_read_zmp_restart(atom,doguess,iw)
    TYPE(atom_type), INTENT(INOUT)           :: atom
    LOGICAL, INTENT(INOUT)                   :: doguess
    INTEGER, INTENT(IN)                      :: iw

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_read_zmp_restart', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: er, extunit, i, k, l, maxl, n
    INTEGER, DIMENSION(0:3)                  :: maxn, nbas

    INQUIRE(file=atom%zmp_restart_file,exist=atom%doread)

    IF (atom%doread) THEN
      WRITE(iw,fmt="(' ZMP       | Restart file found')")
      extunit=get_unit_number()
      CALL open_file(file_name=atom%zmp_restart_file,file_status="OLD",&
                     file_form="FORMATTED", file_action="READ",&
                     unit_number=extunit)

      READ(extunit,*,IOSTAT=er) nbas

      IF ( er .NE. 0 ) THEN
        WRITE(iw,fmt="(' ZMP       | ERROR! Restart file unreadable')")
        WRITE(iw,fmt="(' ZMP       | ERROR! Starting ZMP calculation form initial atomic guess')") 
        doguess=.TRUE.
        atom%doread =.FALSE.
      ELSE IF ( nbas(1).NE.atom%basis%nbas(1) ) THEN
        WRITE(iw,fmt="(' ZMP       | ERROR! Restart file contains a different basis set')")
        WRITE(iw,fmt="(' ZMP       | ERROR! Starting ZMP calculation form initial atomic guess')")
        doguess=.TRUE.
        atom%doread =.FALSE.
      ELSE
        nbas=atom%basis%nbas
        maxl=atom%state%maxl_occ
        maxn=atom%state%maxn_occ
        n = SIZE(atom%orbitals%wfn,2)
        DO l=0,atom%state%maxl_occ
          DO i=1,MIN(n,atom%state%maxn_occ(l))
            DO k=1,atom%basis%nbas(l)
              READ(extunit,*) atom%orbitals%wfn(k,i,l)
            END DO
          END DO
        END DO
        doguess=.FALSE.
      END IF
      CALL close_file(unit_number=extunit)
    ELSE
      WRITE(iw,fmt="(' ZMP       | WARNING! Restart file not found')")
      WRITE(iw,fmt="(' ZMP       | WARNING! Starting ZMP calculation form initial atomic guess')")
    ENDIF
  END SUBROUTINE atom_read_zmp_restart
! *****************************************************************************
!> \brief   ZMP subroutine to read external density from linear grid of 
!>          density matrix. 
!> \param density ...
!> \param atom ...
!> \param iw ...
!> \date    07.10.2013
!> \author  D. Varsano [daniele.varsano@nano.cnr.it]
!> \version 1.0
! *****************************************************************************
  SUBROUTINE atom_read_external_density(density,atom,iw)
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: density
    TYPE(atom_type), INTENT(INOUT)           :: atom
    INTEGER, INTENT(IN)                      :: iw

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_read_external_density', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: filename
    INTEGER                                  :: extunit, ir, j, k, l, &
                                                maxl_occ, maxnbas, nbas, nr
    LOGICAL                                  :: ldm
    REAL(KIND=dp)                            :: rr
    REAL(KIND=dp), ALLOCATABLE               :: pmatread(:,:,:)
    TYPE(cp_error_type)                      :: error

    filename=atom%ext_file
    ldm=atom%dm
    extunit=get_unit_number()

    CALL open_file(file_name=filename,file_status="OLD",&
                   file_form="FORMATTED", file_action="READ",&
                   unit_number=extunit)

    IF (.not.ldm) THEN
      READ(extunit,*) nr

      IF (nr.ne.atom%basis%grid%nr) THEN
        IF (iw > 0 ) WRITE(iw,fmt="(' ZMP       | ERROR! External grid dimension ',I4,' differs from internal grid ',I4)")&
                                                                                    nr, atom%basis%grid%nr
        IF (iw > 0 ) WRITE(iw,fmt="(' ZMP       | ERROR! Stopping ZMP calculation')")
        STOP
      ENDIF

      DO ir=1,nr
        READ(extunit,*) rr,density(ir)
        IF (ABS(rr-atom%basis%grid%rad(ir)).GT. atom%zmpgrid_tol) THEN
          IF (iw > 0 ) WRITE(iw,fmt="(' ZMP       | ERROR! Grid points do not coincide: ')")
          IF (iw > 0 ) WRITE(iw,fmt='(" ZMP       |",T20,"R_out[bohr]",T36,"R_in[bohr]",T61,"R_diff[bohr]")')
          IF (iw > 0 ) WRITE(iw,fmt='(" ZMP       |",T14,E24.15,T39,E24.15,T64,E24.15)') &
                     rr, atom%basis%grid%rad(ir), ABS(rr-atom%basis%grid%rad(ir))
          STOP
        ENDIF
      ENDDO
      CALL close_file(unit_number=extunit)
    ELSE
      READ(extunit,*) maxl_occ
      maxnbas=MAXVAL(atom%basis%nbas)
      ALLOCATE(pmatread(maxnbas,maxnbas,0:maxl_occ))
      pmatread=0.0
      DO l=0,maxl_occ
        nbas=atom%basis%nbas(l)
        READ(extunit,*) ! Read empty line 
        DO k=1,nbas
          READ(extunit,*) (pmatread(j,k,l),j=1,k)
          DO j=1,k
            pmatread(k,j,l)=pmatread(j,k,l)
          ENDDO
        ENDDO
      ENDDO

      CALL close_file(unit_number=extunit)
  
      CALL atom_density(density,pmatread,atom%basis,maxl_occ,typ="RHO",error=error)

      extunit=get_unit_number()
      CALL open_file(file_name="rho_target.dat",file_status="UNKNOWN",&
                     file_form="FORMATTED", file_action="WRITE", unit_number=extunit)

      IF (iw > 0 ) WRITE(iw,fmt="(' ZMP       | Writing target density from density matrix')")

      WRITE(extunit,fmt='("# Target density built from density matrix : ",A20)') filename
      WRITE(extunit,fmt='("#",T10,"R[bohr]",T36,"Rho[au]")')

      nr=atom%basis%grid%nr

      DO ir=1,nr
        WRITE(extunit,fmt='(T1,E24.15,T26,E24.15)') &
                            atom%basis%grid%rad(ir), density(ir) 
      ENDDO
      DEALLOCATE(pmatread)

      CALL close_file(unit_number=extunit) 

    ENDIF

  END SUBROUTINE atom_read_external_density

! *****************************************************************************
!> \brief ZMP subroutine to read external v_xc in the atomic code
!>
!> \param vxc ...
!> \param atom ...
!> \param iw ...
!> \author D. Varsano [daniele.varsano@nano.cnr.it]
! *****************************************************************************
  SUBROUTINE atom_read_external_vxc(vxc,atom,iw)
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: vxc
    TYPE(atom_type), INTENT(INOUT)           :: atom
    INTEGER, INTENT(IN)                      :: iw

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_read_external_vxc', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: adum, filename
    INTEGER                                  :: extunit, ir, nr
    REAL(KIND=dp)                            :: rr

    filename=atom%ext_vxc_file
    extunit=get_unit_number()

    CALL open_file(file_name=filename,file_status="OLD",&
                   file_form="FORMATTED", file_action="READ",&
                   unit_number=extunit)

    READ(extunit,*)
    READ(extunit,*) adum, nr

    IF (nr.ne.atom%basis%grid%nr) THEN
      IF (iw > 0 ) WRITE(iw,fmt="(' ZMP       | ERROR! External grid dimension ',I4,' differs from internal grid ',I4)")&
                                                                                  nr, atom%basis%grid%nr
      IF (iw > 0 ) WRITE(iw,fmt="(' ZMP       | ERROR! Stopping ZMP calculation')")
      STOP
    ENDIF
    DO ir=1,nr
      READ(extunit,*) rr,vxc(ir)
      IF (ABS(rr-atom%basis%grid%rad(ir)).GT. atom%zmpvxcgrid_tol) THEN
        IF (iw > 0 ) WRITE(iw,fmt="(' ZMP       | ERROR! Grid points do not coincide: ')")
        IF (iw > 0 ) WRITE(iw,fmt='(" ZMP       |",T20,"R_out[bohr]",T36,"R_in[bohr]",T61,"R_diff[bohr]")')
        IF (iw > 0 ) WRITE(iw,fmt='(" ZMP       |",T14,E24.15,T39,E24.15,T64,E24.15)')  &
                      rr, atom%basis%grid%rad(ir), ABS(rr-atom%basis%grid%rad(ir))
        STOP
      ENDIF
    ENDDO

  END SUBROUTINE atom_read_external_vxc

! *****************************************************************************
!> \brief ...
!> \param charge ...
!> \param wfn ...
!> \param rcov ...
!> \param l ...
!> \param basis ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_orbital_charge(charge,wfn,rcov,l,basis,error)
    REAL(KIND=dp), INTENT(OUT)               :: charge
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: wfn
    REAL(KIND=dp), INTENT(IN)                :: rcov
    INTEGER, INTENT(IN)                      :: l
    TYPE(atom_basis_type), INTENT(IN)        :: basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_orbital_charge', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, j, m, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: den

    failure = .FALSE.
    charge=0._dp
    m=SIZE(basis%bf,1)
    ALLOCATE(den(m),STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    n=basis%nbas(l)
    den=0._dp
    DO i=1,n
      DO j=1,n
        ff=wfn(i)*wfn(j)
        den(1:m) = den(1:m) + ff*basis%bf(1:m,i,l)*basis%bf(1:m,j,l)
      END DO
    END DO
    DO i=1,m
      IF(basis%grid%rad(i) > rcov) den(i)=0._dp
    END DO
    charge = SUM ( den(1:m) * basis%grid%wr(1:m) )
    DEALLOCATE(den,STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE atom_orbital_charge
! *****************************************************************************
!> \brief ...
!> \param corden ...
!> \param gthpot ...
!> \param typ ...
!> \param rr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_core_density(corden,gthpot,typ,rr,error)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: corden
    TYPE(atom_gthpot_type), INTENT(IN)       :: gthpot
    CHARACTER(LEN=*), OPTIONAL               :: typ
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_core_density', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3)                         :: my_typ
    INTEGER                                  :: i, ierr, j, m, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: fe, rc

    my_typ = "RHO"
    failure = .FALSE.
    IF(PRESENT(typ)) my_typ = typ(1:3)

    IF(gthpot%nlcc) THEN
      m = SIZE(corden)
      ALLOCATE(fe(m),rc(m),STAT=ierr)
      CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      n=gthpot%nexp_nlcc
      DO i=1,n
        a=gthpot%alpha_nlcc(i)
        ! note that all terms are computed with rc, not rr
        rc(:)=rr(:)/a
        fe(:)=EXP(-0.5_dp*rc(:)*rc(:))
        DO j=1,gthpot%nct_nlcc(i)
          IF(my_typ == "RHO") THEN
            corden(:) = corden(:) + fe(:)*rc**(2*j-2)*gthpot%cval_nlcc(j,i)
          ELSE IF(my_typ == "DER") THEN
            corden(:) = corden(:) - fe(:)*rc**(2*j-1)*gthpot%cval_nlcc(j,i)/a
            IF(j>1) THEN
              corden(:) = corden(:) + REAL(2*j-2,dp)*fe(:)*rc**(2*j-3)*gthpot%cval_nlcc(j,i)/a
            END IF
          ELSE
            CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
        END DO
      END DO
      DEALLOCATE(fe,rc,STAT=ierr)
      CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE atom_core_density
! *****************************************************************************
!> \brief ...
!> \param locpot ...
!> \param gthpot ...
!> \param rr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_local_potential(locpot,gthpot,rr,error)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: locpot
    TYPE(atom_gthpot_type), INTENT(IN)       :: gthpot
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_local_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, j, m, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: fe, rc

    failure = .FALSE.
    m = SIZE(locpot)
    ALLOCATE(fe(m),rc(m),STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    rc(:)=rr(:)/gthpot%rc
    DO i=1,m
       locpot(i) = -gthpot%zion*erf(rc(i)/SQRT(2._dp))/rr(i)
    END DO
    n=gthpot%ncl
    fe(:)=EXP(-0.5_dp*rc(:)*rc(:))
    DO i=1,n
       locpot(:) = locpot(:) + fe(:)*rc**(2*i-2)*gthpot%cl(i)
    END DO
    IF(gthpot%lpotextended) THEN
       DO j=1,gthpot%nexp_lpot
          a = gthpot%alpha_lpot(j) 
          rc(:)=rr(:)/a
          fe(:)=EXP(-0.5_dp*rc(:)*rc(:))
          n = gthpot%nct_lpot(j)
          DO i=1,n
             locpot(:) = locpot(:) + fe(:)*rc**(2*i-2)*gthpot%cval_lpot(i,j)
          END DO
       END DO
    END IF
    DEALLOCATE(fe,rc,STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE atom_local_potential
! *****************************************************************************
!> \brief ...
!> \param rmax ...
!> \param wfn ...
!> \param rcov ...
!> \param l ...
!> \param basis ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_orbital_max(rmax,wfn,rcov,l,basis,error)
    REAL(KIND=dp), INTENT(OUT)               :: rmax
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: wfn
    REAL(KIND=dp), INTENT(IN)                :: rcov
    INTEGER, INTENT(IN)                      :: l
    TYPE(atom_basis_type), INTENT(IN)        :: basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_orbital_max', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, m, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dorb

    failure = .FALSE.
    m=SIZE(basis%bf,1)
    ALLOCATE(dorb(m),STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    n=basis%nbas(l)
    dorb=0._dp
    DO i=1,n
      ff=wfn(i)
      dorb(1:m) = dorb(1:m) + ff*basis%dbf(1:m,i,l)
    END DO
    rmax=-1._dp
    DO i=1,m-1
      IF(basis%grid%rad(i) < 2*rcov) THEN
        IF (dorb(i)*dorb(i+1) < 0._dp) THEN
          rmax=MAX(rmax,basis%grid%rad(i))
        END IF
      END IF
    END DO
    DEALLOCATE(dorb,STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE atom_orbital_max
! *****************************************************************************
!> \brief ...
!> \param node ...
!> \param wfn ...
!> \param rcov ...
!> \param l ...
!> \param basis ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_orbital_nodes(node,wfn,rcov,l,basis,error)
    INTEGER, INTENT(OUT)                     :: node
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: wfn
    REAL(KIND=dp), INTENT(IN)                :: rcov
    INTEGER, INTENT(IN)                      :: l
    TYPE(atom_basis_type), INTENT(IN)        :: basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_orbital_nodes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, m, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: orb

    failure = .FALSE.
    node=0
    m=SIZE(basis%bf,1)
    ALLOCATE(orb(m),STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    n=basis%nbas(l)
    orb=0._dp
    DO i=1,n
      ff=wfn(i)
      orb(1:m) = orb(1:m) + ff*basis%bf(1:m,i,l)
    END DO
    DO i=1,m-1
      IF(basis%grid%rad(i) < rcov) THEN
        IF (orb(i)*orb(i+1) < 0._dp) node=node+1
      END IF
    END DO
    DEALLOCATE(orb,STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE atom_orbital_nodes
! *****************************************************************************
!> \brief ...
!> \param value ...
!> \param wfn ...
!> \param basis ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_wfnr0(value,wfn,basis,error)
    REAL(KIND=dp), INTENT(OUT)               :: value
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: wfn
    TYPE(atom_basis_type), INTENT(IN)        :: basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_wfnr0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, m, n

    value=0._dp
    m=MAXVAL(MINLOC(basis%grid%rad))
    n=basis%nbas(0)
    DO i=1,n
      value = value + wfn(i)*basis%bf(m,i,0)
    END DO
  END SUBROUTINE atom_wfnr0
! *****************************************************************************
!> \brief ...
!> \param hmat ...
!> \param umat ...
!> \param orb ...
!> \param ener ...
!> \param nb ...
!> \param nv ...
!> \param maxl ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_solve(hmat,umat,orb,ener,nb,nv,maxl,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: hmat, umat
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: orb
    REAL(KIND=dp), DIMENSION(:, 0:), &
      INTENT(INOUT)                          :: ener
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nb, nv
    INTEGER, INTENT(IN)                      :: maxl
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_solve', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, info, l, lwork, m, n
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: w, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: a

    failure = .FALSE.
    CPPrecondition(ALL(nb >= nv),cp_failure_level,routineP,error,failure)

    orb = 0._dp
    DO l=0,maxl
      n = nb(l)
      m = nv(l)
      IF ( n > 0 .AND. m > 0 ) THEN
        lwork = 10*m
        ALLOCATE(a(n,n),w(n),work(lwork),STAT=ierr)
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        a(1:m,1:m) = MATMUL(TRANSPOSE(umat(1:n,1:m,l)),MATMUL(hmat(1:n,1:n,l),umat(1:n,1:m,l)))
        CALL lapack_ssyev ( "V", "U", m, a(1:m,1:m), m, w(1:m), work, lwork, info )
        a(1:n,1:m) = MATMUL(umat(1:n,1:m,l),a(1:m,1:m))

        m = MIN(m,SIZE(orb,2))
        orb(1:n,1:m,l) = a(1:n,1:m)
        ener(1:m,l) = w(1:m)

        DEALLOCATE(a,w,work,STAT=ierr)
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      END IF
    END DO

  END SUBROUTINE atom_solve
! *****************************************************************************
!> \brief ...
!> \param fun ...
!> \param deps ...
!> \retval nc ...
! *****************************************************************************
  FUNCTION prune_grid ( fun, deps ) RESULT(nc)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fun
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: deps
    INTEGER                                  :: nc

    INTEGER                                  :: i, nr
    REAL(KIND=dp)                            :: meps

    meps = 1.e-14_dp
    IF ( PRESENT(deps) ) meps=deps

    nr = SIZE(fun)
    nc = 0
    DO i=nr,1,-1
      IF (ABS(fun(i)) > meps) THEN
        nc = i
        EXIT
      END IF
    END DO

  END FUNCTION prune_grid
! *****************************************************************************
!> \brief ...
!> \param fun ...
!> \param grid ...
!> \retval integral ...
! *****************************************************************************
  FUNCTION integrate_grid_function1 ( fun, grid ) RESULT(integral)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fun
    TYPE(grid_atom_type)                     :: grid
    REAL(KIND=dp)                            :: integral

    INTEGER                                  :: nc

    nc = SIZE(fun)
    integral = SUM ( fun(1:nc) * grid%wr(1:nc) )

  END FUNCTION integrate_grid_function1

! *****************************************************************************
!> \brief ...
!> \param fun1 ...
!> \param fun2 ...
!> \param grid ...
!> \retval integral ...
! *****************************************************************************
  FUNCTION integrate_grid_function2 ( fun1, fun2, grid ) RESULT(integral)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fun1, fun2
    TYPE(grid_atom_type)                     :: grid
    REAL(KIND=dp)                            :: integral

    INTEGER                                  :: nc

    nc = MIN(SIZE(fun1),SIZE(fun2))
    integral = SUM ( fun1(1:nc) * fun2(1:nc) * grid%wr(1:nc) )

  END FUNCTION integrate_grid_function2

! *****************************************************************************
!> \brief ...
!> \param fun1 ...
!> \param fun2 ...
!> \param fun3 ...
!> \param grid ...
!> \retval integral ...
! *****************************************************************************
  FUNCTION integrate_grid_function3 ( fun1, fun2, fun3, grid ) RESULT(integral)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fun1, fun2, fun3
    TYPE(grid_atom_type)                     :: grid
    REAL(KIND=dp)                            :: integral

    INTEGER                                  :: nc

    nc = MIN(SIZE(fun1),SIZE(fun2),SIZE(fun3))
    integral = SUM ( fun1(1:nc) * fun2(1:nc) * fun3(1:nc) * grid%wr(1:nc) )

  END FUNCTION integrate_grid_function3
! *****************************************************************************
!> \brief ...
!> \param cpot ...
!> \param density ...
!> \param grid ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE coulomb_potential_numeric ( cpot, density, grid, error )
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: cpot
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: density
    TYPE(grid_atom_type), INTENT(IN)         :: grid
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coulomb_potential_numeric', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nc
    LOGICAL                                  :: failure
    REAL(dp)                                 :: int1, int2
    REAL(dp), DIMENSION(:), POINTER          :: r, wr

    failure = .FALSE.
    nc = MIN(SIZE(cpot),SIZE(density))
    r => grid%rad
    wr => grid%wr

    int1 = fourpi * integrate_grid ( density, grid )
    int2 = 0._dp
    cpot(nc:) = int1/r(nc:)

    ! test that grid is decreasing
    CPPostcondition(r(1)>r(nc),cp_failure_level,routineP,error,failure)
    DO i = 1, nc
      cpot(i) = int1/r(i) + int2
      int1 = int1 - fourpi * density(i)*wr(i)
      int2 = int2 + fourpi * density(i)*wr(i)/r(i)
    END DO

  END SUBROUTINE coulomb_potential_numeric

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param cpot ...
!> \param pmat ...
!> \param basis ...
!> \param grid ...
!> \param maxl ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE coulomb_potential_analytic ( cpot, pmat, basis, grid, maxl, error )
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: cpot
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: pmat
    TYPE(atom_basis_type), INTENT(IN)        :: basis
    TYPE(grid_atom_type)                     :: grid
    INTEGER, INTENT(IN)                      :: maxl
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coulomb_potential_analytic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, j, k, l, m, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, b, ff, oab, sab
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: erfa, expa, z
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: unp

    failure = .FALSE.
    m = SIZE(cpot)
    ALLOCATE ( erfa(1:m), expa(1:m), z(1:m), STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

    cpot = 0._dp

    DO l = 0, maxl
      IF ( MAXVAL(ABS(pmat(:,:,l))) < 1.e-14_dp ) CYCLE
      SELECT CASE ( basis%basis_type )
        CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        CASE ( GTO_BASIS )
          DO i = 1, basis%nbas(l)
            DO j = i, basis%nbas(l)
              IF ( ABS(pmat(i,j,l)) < 1.e-14_dp ) CYCLE
              ff = pmat(i,j,l)
              IF (i /= j) ff = 2._dp*ff
              a = basis%am(i,l)
              b = basis%am(j,l)
              sab = SQRT(a+b)
              oab = rootpi/(a+b)**(l+1.5_dp)*ff
              z(:) = sab*grid%rad(:)
              DO k=1,SIZE(erfa)
                erfa(k)= oab*erf(z(k))/grid%rad(k)
              END DO
              expa(:)= EXP(-z(:)**2)*ff/(a+b)**(l+1)
              SELECT CASE (l)
                CASE DEFAULT
                  CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                CASE (0)
                  cpot(:) = cpot(:) + 0.25_dp*erfa(:)
                CASE (1)
                  cpot(:) = cpot(:) + 0.375_dp*erfa(:)-0.25_dp*expa(:)
                CASE (2)
                  cpot(:) = cpot(:) + 0.9375_dp*erfa(:)-expa(:)*(0.875_dp+0.25_dp*z(:)**2)
                CASE (3)
                  cpot(:) = cpot(:) + 3.28125_dp*erfa(:)-expa(:)*(3.5625_dp+1.375_dp*z(:)**2+0.25*z(:)**4)
              END SELECT
            END DO
          END DO
        CASE ( CGTO_BASIS )
          n = basis%nprim(l)
          m = basis%nbas(l)
          ALLOCATE ( unp(n,n), STAT=ierr )
          CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

          unp(1:n,1:n) = MATMUL ( MATMUL(basis%cm(1:n,1:m,l),pmat(1:m,1:m,l)),&
                                  TRANSPOSE(basis%cm(1:n,1:m,l)))
          DO i = 1, basis%nprim(l)
            DO j = i, basis%nprim(l)
              IF ( ABS(unp(i,j)) < 1.e-14_dp ) CYCLE
              ff = unp(i,j)
              IF (i /= j) ff = 2._dp*ff
              a = basis%am(i,l)
              b = basis%am(j,l)
              sab = SQRT(a+b)
              oab = rootpi/(a+b)**(l+1.5_dp)*ff
              z(:) = sab*grid%rad(:)
              DO k=1,SIZE(erfa)
                erfa(k)= oab*erf(z(k))/grid%rad(k)
              END DO
              expa(:)= EXP(-z(:)**2)*ff/(a+b)**(l+1)
              SELECT CASE (l)
                CASE DEFAULT
                  CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                CASE (0)
                  cpot(:) = cpot(:) + 0.25_dp*erfa(:)
                CASE (1)
                  cpot(:) = cpot(:) + 0.375_dp*erfa(:)-0.25_dp*expa(:)
                CASE (2)
                  cpot(:) = cpot(:) + 0.9375_dp*erfa(:)-expa(:)*(0.875_dp+0.25_dp*z(:)**2)
                CASE (3)
                  cpot(:) = cpot(:) + 3.28125_dp*erfa(:)-expa(:)*(3.5625_dp+1.375_dp*z(:)**2+0.25*z(:)**4)
              END SELECT
            END DO
          END DO

          DEALLOCATE ( unp, STAT=ierr )
          CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      END SELECT
    END DO
    DEALLOCATE ( erfa, expa, z, STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE coulomb_potential_analytic

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param kmat ...
!> \param state ...
!> \param occ ...
!> \param wfn ...
!> \param basis ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE exchange_numeric ( kmat, state, occ, wfn, basis, error )
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: kmat
    TYPE(atom_state), INTENT(IN)             :: state
    REAL(KIND=dp), DIMENSION(0:, :), &
      INTENT(IN)                             :: occ
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: wfn
    TYPE(atom_basis_type), INTENT(IN)        :: basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'exchange_numeric', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ia, ib, ierr, k, lad, lbc, &
                                                lh, ll, nbas, norb, nr, nu
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: almn
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cpot, nai, nbi
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: orb
    REAL(KIND=dp), DIMENSION(0:maxfac)       :: arho

    failure = .FALSE.
    arho = 0._dp
    DO ll=0,maxfac,2
      lh=ll/2
      arho(ll)=fac(ll)/fac(lh)**2
    END DO

    kmat = 0._dp

    nr = basis%grid%nr
    ALLOCATE ( nai(nr), nbi(nr), cpot(nr), STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

    DO lad = 0, state%maxl_calc
      DO lbc = 0 , state%maxl_occ
        norb = state%maxn_occ(lbc)
        nbas = basis%nbas(lbc)
        ! calculate orbitals for angmom lbc
        ALLOCATE ( orb(nr,norb), STAT=ierr )
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        orb = 0._dp
        DO i=1,norb
          DO k=1,nbas
            orb(:,i) = orb(:,i) + wfn(k,i,lbc)*basis%bf(:,k,lbc)
          END DO
        END DO
        DO nu = ABS(lad-lbc),lad+lbc,2
          almn = arho(-lad+lbc+nu)*arho(lad-lbc+nu)*arho(lad+lbc-nu)/(REAL(lad+lbc+nu+1,dp)*arho(lad+lbc+nu))
          almn = -0.5_dp*almn

          DO ia=1,basis%nbas(lad)
            DO i=1,norb
              nai(:) = orb(:,i)*basis%bf(:,ia,lad)
              CALL potential_numeric(cpot,nai,nu,basis%grid,error)
              DO ib=1,basis%nbas(lad)
                kmat(ia,ib,lad) = kmat(ia,ib,lad) + almn*occ(lbc,i)* &
                  integrate_grid(cpot,orb(:,i),basis%bf(:,ib,lad),basis%grid)
              END DO
            END DO
          END DO

        END DO
        DEALLOCATE ( orb, STAT=ierr )
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      END DO
    END DO

    DEALLOCATE ( nai, nbi, cpot, STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE exchange_numeric

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param kmat ...
!> \param state ...
!> \param occ ...
!> \param wfn ...
!> \param basis ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE exchange_semi_analytic ( kmat, state, occ, wfn, basis, error )
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: kmat
    TYPE(atom_state), INTENT(IN)             :: state
    REAL(KIND=dp), DIMENSION(0:, :), &
      INTENT(IN)                             :: occ
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: wfn
    TYPE(atom_basis_type), INTENT(IN)        :: basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'exchange_semi_analytic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ia, ib, ierr, k, lad, lbc, &
                                                lh, ll, nbas, norb, nr, nu
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: almn
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cpot, nai, nbi
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: orb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: pot
    REAL(KIND=dp), DIMENSION(0:maxfac)       :: arho

    failure = .FALSE.
    arho = 0._dp
    DO ll=0,maxfac,2
      lh=ll/2
      arho(ll)=fac(ll)/fac(lh)**2
    END DO

    kmat = 0._dp

    nr = basis%grid%nr
    nbas = MAXVAL(basis%nbas)
    ALLOCATE ( pot(nr,nbas,nbas), STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( nai(nr), nbi(nr), cpot(nr), STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

    DO lad = 0, state%maxl_calc
      DO lbc = 0 , state%maxl_occ
        norb = state%maxn_occ(lbc)
        nbas = basis%nbas(lbc)
        ! calculate orbitals for angmom lbc
        ALLOCATE ( orb(nr,norb), STAT=ierr )
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        orb = 0._dp
        DO i=1,norb
          DO k=1,nbas
            orb(:,i) = orb(:,i) + wfn(k,i,lbc)*basis%bf(:,k,lbc)
          END DO
        END DO
        DO nu = ABS(lad-lbc),lad+lbc,2
          almn = arho(-lad+lbc+nu)*arho(lad-lbc+nu)*arho(lad+lbc-nu)/(REAL(lad+lbc+nu+1,dp)*arho(lad+lbc+nu))
          almn = -0.5_dp*almn
          ! calculate potential for basis function pair (lad,lbc)
          pot = 0._dp
          CALL potential_analytic ( pot, lad, lbc, nu, basis, error )
          DO ia=1,basis%nbas(lad)
            DO i=1,norb
              cpot = 0._dp
              DO k=1,nbas
                cpot(:) = cpot(:) + pot(:,ia,k)*wfn(k,i,lbc)
              END DO
              DO ib=1,basis%nbas(lad)
                kmat(ia,ib,lad) = kmat(ia,ib,lad) + almn*occ(lbc,i)* &
                  integrate_grid(cpot,orb(:,i),basis%bf(:,ib,lad),basis%grid)
              END DO
            END DO
          END DO
        END DO
        DEALLOCATE ( orb, STAT=ierr )
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      END DO
    END DO

    DEALLOCATE ( pot, STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( nai, nbi, cpot, STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE exchange_semi_analytic

! *****************************************************************************
!> \brief ...
!> \param cpot ...
!> \param density ...
!> \param nu ...
!> \param grid ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE potential_numeric ( cpot, density, nu, grid, error )
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: cpot
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: density
    INTEGER, INTENT(IN)                      :: nu
    TYPE(grid_atom_type), INTENT(IN)         :: grid
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'potential_numeric', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nc
    LOGICAL                                  :: failure
    REAL(dp)                                 :: int1, int2
    REAL(dp), DIMENSION(:), POINTER          :: r, wr

    failure = .FALSE.
    nc = MIN(SIZE(cpot),SIZE(density))
    r => grid%rad
    wr => grid%wr

    int1 = integrate_grid ( density, r**nu, grid )
    int2 = 0._dp
    cpot(nc:) = int1/r(nc:)**(nu+1)

    ! test that grid is decreasing
    CPPostcondition(r(1)>r(nc),cp_failure_level,routineP,error,failure)
    DO i = 1, nc
      cpot(i) = int1/r(i)**(nu+1) + int2*r(i)**nu
      int1 = int1 - r(i)**(nu)*density(i)*wr(i)
      int2 = int2 + density(i)*wr(i)/r(i)**(nu+1)
    END DO

  END SUBROUTINE potential_numeric
! *****************************************************************************
!> \brief ...
!> \param cpot ...
!> \param la ...
!> \param lb ...
!> \param nu ...
!> \param basis ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE potential_analytic ( cpot, la, lb, nu, basis, error )
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: cpot
    INTEGER, INTENT(IN)                      :: la, lb, nu
    TYPE(atom_basis_type), INTENT(IN)        :: basis
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'potential_analytic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, j, k, l, ll, m
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, b, oab, sab
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: erfa, expa, z

    failure = .FALSE.
    m = SIZE(cpot,1)
    ALLOCATE ( erfa(1:m), expa(1:m), z(1:m), STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

    ll = la+lb

    cpot = 0._dp

    SELECT CASE ( basis%basis_type )
      CASE DEFAULT
        CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
      CASE ( GTO_BASIS )
        DO i = 1, basis%nbas(la)
          DO j = 1, basis%nbas(lb)
            a = basis%am(i,la)
            b = basis%am(j,lb)
            sab = SQRT(a+b)
            oab = dfac(ll+nu+1)*rootpi/(sab**(ll+2))/(2._dp**((ll+nu)/2+2))
            z(:) = sab*basis%grid%rad(:)
            DO k=1,SIZE(erfa)
              erfa(k)= oab*erf(z(k))/z(k)**(nu+1)
            END DO
            cpot(:,i,j) = erfa(:)
            expa(:)= EXP(-z(:)**2)/(sab**(ll+2))/(2._dp**((ll+nu)/2+2))
            SELECT CASE (ll)
              CASE DEFAULT
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
              CASE (0)
                CPPostcondition(nu==0,cp_failure_level,routineP,error,failure)
              CASE (1)
                CPPostcondition(nu==1,cp_failure_level,routineP,error,failure)
                cpot(:,i,j) = cpot(:,i,j)-6._dp*expa(:)/z(:)
              CASE (2)
                SELECT CASE (nu)
                  CASE DEFAULT
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                  CASE (0)
                    cpot(:,i,j) = cpot(:,i,j)-2._dp*expa(:)
                  CASE (2)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(20._dp+30._dp/z(:)**2)
                END SELECT
              CASE (3)
                SELECT CASE (nu)
                  CASE DEFAULT
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                  CASE (1)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(12._dp*z(:)+30._dp/z(:))
                  CASE (3)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(56._dp*z(:)+140._dp/z(:)+210._dp*z(:)**3)
                END SELECT
              CASE (4)
                SELECT CASE (nu)
                  CASE DEFAULT
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                  CASE (0)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(4._dp*z(:)**2+14._dp)
                  CASE (2)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(40._dp*z(:)**2+140._dp+210._dp/z(:)**2)
                  CASE (4)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(144._dp*z(:)**2+504._dp+1260._dp/z(:)**2+1890._dp/z(:)**4)
                END SELECT
              CASE (5)
                SELECT CASE (nu)
                  CASE DEFAULT
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                  CASE (1)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(24._dp*z(:)**3+108._dp*z(:)+210._dp/z(:))
                  CASE (3)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(112._dp*z(:)**3+504._dp*z(:)+1260._dp/z(:)+1890._dp/z(:)**3)
                  CASE (5)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(352._dp*z(:)**3+1584._dp*z(:)+5544._dp/z(:)+&
                                  13860._dp/z(:)**3+20790._dp/z(:)**5)
                END SELECT
              CASE (6)
                SELECT CASE (nu)
                  CASE DEFAULT
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                  CASE (0)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(8._dp*z(:)**4+44._dp*z(:)**2+114._dp)
                  CASE (2)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(80._dp*z(:)**4+440._dp*z(:)**2+1260._dp+1896._dp/z(:)**2)
                  CASE (4)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(288._dp*z(:)**4+1584._dp*z(:)**2+5544._dp+&
                                  13860._dp/z(:)**2+20790._dp/z(:)**4)
                  CASE (6)
                    cpot(:,i,j) = cpot(:,i,j)-expa(:)*(832._dp*z(:)**4+4576._dp*z(:)**2+20592._dp+&
                                  72072._dp/z(:)**2+180180._dp/z(:)**4+270270._dp/z(:)**6)
                END SELECT
            END SELECT
          END DO
        END DO
      CASE ( CGTO_BASIS )
        DO i = 1, basis%nprim(la)
          DO j = 1, basis%nprim(lb)
            a = basis%am(i,la)
            b = basis%am(j,lb)
            sab = SQRT(a+b)
            oab = dfac(ll+nu+1)*rootpi/sab**(ll+2)/2._dp**((ll+nu)/2+2)
            z(:) = sab*basis%grid%rad(:)
            DO k=1,SIZE(erfa)
              erfa(k)= oab*erf(z(k))/z(k)**(nu+1)
            END DO
            expa(:)= EXP(-z(:)**2)/sab**(ll+2)/2._dp**((ll+nu)/2+2)
            SELECT CASE (ll)
              CASE DEFAULT
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
              CASE (0)
                CPPostcondition(nu==0,cp_failure_level,routineP,error,failure)
              CASE (1)
                CPPostcondition(nu==1,cp_failure_level,routineP,error,failure)
                erfa(:) = erfa(:)-6._dp*expa(:)/z(:)
              CASE (2)
                SELECT CASE (nu)
                  CASE DEFAULT
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                  CASE (0)
                    erfa(:) = erfa(:)-2._dp*expa(:)
                  CASE (2)
                    erfa(:) = erfa(:)-expa(:)*(20._dp+30._dp/z(:)**2)
                END SELECT
              CASE (3)
                SELECT CASE (nu)
                  CASE DEFAULT
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                  CASE (1)
                    erfa(:) = erfa(:)-expa(:)*(12._dp*z(:)+30._dp/z(:))
                  CASE (3)
                    erfa(:) = erfa(:)-expa(:)*(56._dp*z(:)+140._dp/z(:)+210._dp*z(:)**3)
                END SELECT
              CASE (4)
                SELECT CASE (nu)
                  CASE DEFAULT
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                  CASE (0)
                    erfa(:) = erfa(:)-expa(:)*(4._dp*z(:)**2+14._dp)
                  CASE (2)
                    erfa(:) = erfa(:)-expa(:)*(40._dp*z(:)**2+140._dp+210._dp/z(:)**2)
                  CASE (4)
                    erfa(:) = erfa(:)-expa(:)*(144._dp*z(:)**2+504._dp+1260._dp/z(:)**2+1890._dp/z(:)**4)
                END SELECT
              CASE (5)
                SELECT CASE (nu)
                  CASE DEFAULT
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                  CASE (1)
                    erfa(:) = erfa(:)-expa(:)*(24._dp*z(:)**3+108._dp*z(:)+210._dp/z(:))
                  CASE (3)
                    erfa(:) = erfa(:)-expa(:)*(112._dp*z(:)**3+504._dp*z(:)+1260._dp/z(:)+1890._dp/z(:)**3)
                  CASE (5)
                    erfa(:) = erfa(:)-expa(:)*(352._dp*z(:)**3+1584._dp*z(:)+5544._dp/z(:)+&
                                  13860._dp/z(:)**3+20790._dp/z(:)**5)
                END SELECT
              CASE (6)
                SELECT CASE (nu)
                  CASE DEFAULT
                    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                  CASE (0)
                    erfa(:) = erfa(:)-expa(:)*(8._dp*z(:)**4+44._dp*z(:)**2+114._dp)
                  CASE (2)
                    erfa(:) = erfa(:)-expa(:)*(80._dp*z(:)**4+440._dp*z(:)**2+1260._dp+1896._dp/z(:)**2)
                  CASE (4)
                    erfa(:) = erfa(:)-expa(:)*(288._dp*z(:)**4+1584._dp*z(:)**2+5544._dp+&
                                  13860._dp/z(:)**2+20790._dp/z(:)**4)
                  CASE (6)
                    erfa(:) = erfa(:)-expa(:)*(832._dp*z(:)**4+4576._dp*z(:)**2+20592._dp+&
                                  72072._dp/z(:)**2+180180._dp/z(:)**4+270270._dp/z(:)**6)
                END SELECT
            END SELECT
            DO k = 1, basis%nbas(la)
              DO l = 1, basis%nbas(lb)
                cpot(:,k,l) = cpot(:,k,l) + erfa(:)*basis%cm(i,k,la)*basis%cm(j,l,lb)
              END DO
            END DO
          END DO
        END DO

    END SELECT

    DEALLOCATE ( erfa, expa, z, STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE potential_analytic
! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param opmat ...
!> \param pmat ...
!> \param error ...
!> \retval trace ...
! *****************************************************************************
  FUNCTION atom_trace(opmat,pmat,error) RESULT(trace)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: opmat, pmat
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    REAL(KIND=dp)                            :: trace

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_trace', &
      routineP = moduleN//':'//routineN

    trace = SUM(opmat(:,:,:)*pmat(:,:,:))

  END FUNCTION atom_trace

! *****************************************************************************
!> \brief ...
!> \param imat ...
!> \param cpot ...
!> \param basis ...
!> \param derivatives ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE numpot_matrix(imat,cpot,basis,derivatives,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: imat
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: cpot
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    INTEGER, INTENT(IN)                      :: derivatives
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'numpot_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, l, n
    LOGICAL                                  :: failure

    failure = .FALSE.

    SELECT CASE (derivatives)
      CASE (0)
        DO l=0,3
          n=basis%nbas(l)
          DO i=1,n
            DO j=i,n
              imat(i,j,l) = imat(i,j,l) + &
                  integrate_grid ( cpot, basis%bf(:,i,l), basis%bf(:,j,l), basis%grid )
              imat(j,i,l) = imat(i,j,l)
            END DO
          END DO
        END DO
      CASE (1)
        DO l=0,3
          n=basis%nbas(l)
          DO i=1,n
            DO j=i,n
              imat(i,j,l) = imat(i,j,l) + &
                  integrate_grid ( cpot, basis%dbf(:,i,l), basis%bf(:,j,l), basis%grid )
              imat(i,j,l) = imat(i,j,l) + &
                  integrate_grid ( cpot, basis%bf(:,i,l), basis%dbf(:,j,l), basis%grid )
              imat(j,i,l) = imat(i,j,l)
            END DO
          END DO
        END DO
      CASE (2)
        DO l=0,3
          n=basis%nbas(l)
          DO i=1,n
            DO j=i,n
              imat(i,j,l) = imat(i,j,l) + &
                 integrate_grid ( cpot, basis%dbf(:,i,l), basis%dbf(:,j,l), basis%grid )
              imat(j,i,l) = imat(i,j,l)
            END DO
          END DO
        END DO
      CASE DEFAULT
        CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE numpot_matrix
! *****************************************************************************
!> \brief ...
!> \param jmat ...
!> \param erint ...
!> \param pmat ...
!> \param nsize ...
!> \param all_nu ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE ceri_contract(jmat,erint,pmat,nsize,all_nu,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: jmat
    TYPE(eri), DIMENSION(:), INTENT(IN)      :: erint
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: pmat
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nsize
    LOGICAL, INTENT(IN), OPTIONAL            :: all_nu
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ceri_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, ij1, ij2, j1, j2, l1, &
                                                l2, ll, n1, n2
    LOGICAL                                  :: failure, have_all_nu
    REAL(KIND=dp)                            :: eint, f1, f2

    failure = .FALSE.
    IF ( PRESENT(all_nu) ) THEN
      have_all_nu = all_nu
    ELSE
      have_all_nu = .FALSE.
    END IF

    jmat(:,:,:) = 0._dp
    ll = 0
    DO l1=0,3
      n1 = nsize(l1)
      DO l2=0,l1
        n2 = nsize(l2)
        ll = ll + 1
        ij1=0
        DO i1=1,n1
          DO j1=i1,n1
            ij1=ij1+1
            f1 = 1._dp
            IF(i1/=j1) f1 = 2._dp
            ij2=0
            DO i2=1,n2
              DO j2=i2,n2
                ij2=ij2+1
                f2 = 1._dp
                IF(i2/=j2) f2 = 2._dp
                eint = erint(ll)%int(ij1,ij2)
                IF ( l1==l2 ) THEN
                  jmat(i1,j1,l1) = jmat(i1,j1,l1) + f2*pmat(i2,j2,l2)*eint
                ELSE
                  jmat(i1,j1,l1) = jmat(i1,j1,l1) + f2*pmat(i2,j2,l2)*eint
                  jmat(i2,j2,l2) = jmat(i2,j2,l2) + f1*pmat(i1,j1,l1)*eint
                END IF
              END DO
            END DO
          END DO
        END DO
        IF ( have_all_nu ) THEN
          ! skip integral blocks with nu/=0
          ll = ll + l2
        END IF
      END DO
    END DO
    DO l1=0,3
      n1 = nsize(l1)
      DO i1=1,n1
        DO j1=i1,n1
          jmat(j1,i1,l1) = jmat(i1,j1,l1)
        END DO
      END DO
    END DO

  END SUBROUTINE ceri_contract
! *****************************************************************************
!> \brief ...
!> \param kmat ...
!> \param erint ...
!> \param pmat ...
!> \param nsize ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE eeri_contract(kmat,erint,pmat,nsize,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: kmat
    TYPE(eri), DIMENSION(:), INTENT(IN)      :: erint
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: pmat
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nsize
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eeri_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, ij1, ij2, j1, j2, l1, &
                                                l2, lh, ll, n1, n2, nu
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: almn, eint, f1, f2
    REAL(KIND=dp), DIMENSION(0:maxfac)       :: arho

    failure = .FALSE.

    arho = 0._dp
    DO ll=0,maxfac,2
      lh=ll/2
      arho(ll)=fac(ll)/fac(lh)**2
    END DO

    kmat(:,:,:) = 0._dp
    ll = 0
    DO l1=0,3
      n1 = nsize(l1)
      DO l2=0,l1
        n2 = nsize(l2)
        DO nu = ABS(l1-l2),l1+l2,2
          almn = arho(-l1+l2+nu)*arho(l1-l2+nu)*arho(l1+l2-nu)/(REAL(l1+l2+nu+1,dp)*arho(l1+l2+nu))
          almn = -0.5_dp*almn
          ll = ll + 1
          ij1=0
          DO i1=1,n1
            DO j1=i1,n1
              ij1=ij1+1
              f1 = 1._dp
              IF(i1/=j1) f1 = 2._dp
              ij2=0
              DO i2=1,n2
                DO j2=i2,n2
                  ij2=ij2+1
                  f2 = 1._dp
                  IF(i2/=j2) f2 = 2._dp
                  eint = erint(ll)%int(ij1,ij2)
                  IF ( l1==l2 ) THEN
                    kmat(i1,j1,l1) = kmat(i1,j1,l1) + f2*almn*pmat(i2,j2,l2)*eint
                  ELSE
                    kmat(i1,j1,l1) = kmat(i1,j1,l1) + f2*almn*pmat(i2,j2,l2)*eint
                    kmat(i2,j2,l2) = kmat(i2,j2,l2) + f1*almn*pmat(i1,j1,l1)*eint
                  END IF
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    DO l1=0,3
      n1 = nsize(l1)
      DO i1=1,n1
        DO j1=i1,n1
          kmat(j1,i1,l1) = kmat(i1,j1,l1)
        END DO
      END DO
    END DO

  END SUBROUTINE eeri_contract
! *****************************************************************************
!> \brief ...
!> \param emat ...
!> \param demax ...
!> \param kmat ...
!> \param pmat ...
!> \param umat ...
!> \param upmat ...
!> \param nval ...
!> \param nbs ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE err_matrix(emat,demax,kmat,pmat,umat,upmat,nval,nbs,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(OUT)                            :: emat
    REAL(KIND=dp), INTENT(OUT)               :: demax
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: kmat, pmat, umat, upmat
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nval, nbs
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'err_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, l, m, n
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: tkmat, tpmat

    failure = .FALSE.
    emat = 0._dp
    DO l=0,3
      n = nval(l)
      m = nbs(l)
      IF (m > 0) THEN
        ALLOCATE(tkmat(1:m,1:m),tpmat(1:m,1:m),STAT=ierr)
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        tkmat = 0._dp
        tpmat = 0._dp
        tkmat(1:m,1:m) = MATMUL(TRANSPOSE(umat(1:n,1:m,l)),MATMUL(kmat(1:n,1:n,l),umat(1:n,1:m,l)))
        tpmat(1:m,1:m) = MATMUL(TRANSPOSE(umat(1:n,1:m,l)),MATMUL(pmat(1:n,1:n,l),umat(1:n,1:m,l)))
        tpmat(1:m,1:m) = MATMUL(upmat(1:m,1:m,l),MATMUL(tpmat(1:m,1:m),upmat(1:m,1:m,l)))

        emat(1:m,1:m,l) = MATMUL(tkmat(1:m,1:m),tpmat(1:m,1:m)) - MATMUL(tpmat(1:m,1:m),tkmat(1:m,1:m))

        DEALLOCATE(tkmat,tpmat,STAT=ierr)
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      END IF
    END DO
    demax = MAXVAL(ABS(emat))

  END SUBROUTINE err_matrix
! *****************************************************************************
!> \brief ...
!> \param density1 ...
!> \param density2 ...
!> \param zcore ...
!> \param state ...
!> \param grid ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE slater_density(density1,density2,zcore,state,grid,error)
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: density1, density2
    INTEGER, INTENT(IN)                      :: zcore
    TYPE(atom_state), INTENT(IN)             :: state
    TYPE(grid_atom_type), INTENT(IN)         :: grid
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'slater_density', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, l, n
    INTEGER, DIMENSION(4, 7)                 :: ne
    REAL(KIND=dp)                            :: a, pf

    density1 = 0._dp
    density2 = 0._dp
    ne(1:4,1:7) = NINT(state%core(0:3,1:7) + state%occ(0:3,1:7))
    DO l=0,state%maxl_occ
      DO i=1,SIZE(state%occ,2)
        IF ( state%occ(l,i) > 0._dp ) THEN
          n = i+l
          a = srules(zcore,ne,n,l)
          pf = 1._dp/SQRT(fac(2*n)) * (2._dp*a)**(n+0.5_dp)
          IF (state%multiplicity == -1) THEN
            density1(:) = density1(:)+state%occ(l,i)/fourpi*(grid%rad(:)**(n-1)*EXP(-a*grid%rad(:))*pf)**2
          ELSE
            density1(:) = density1(:)+state%occa(l,i)/fourpi*(grid%rad(:)**(n-1)*EXP(-a*grid%rad(:))*pf)**2
            density2(:) = density2(:)+state%occb(l,i)/fourpi*(grid%rad(:)**(n-1)*EXP(-a*grid%rad(:))*pf)**2
          END IF
        END IF
      END DO
    END DO

  END SUBROUTINE slater_density
! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param vxc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE wigner_slater_functional(rho,vxc,error)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rho
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: vxc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wigner_slater_functional', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: ec, ex, rs, vc, vx

    vxc = 0._dp
    DO i=1,SIZE(rho)
      IF ( rho(i) > 1.e-20_dp ) THEN
         ex = -0.7385588_dp * rho(i)**0.333333333_dp
         vx = 1.333333333_dp * ex
         rs = (3._dp/fourpi/rho(i))**0.333333333_dp
         ec = -0.88_dp/(rs+7.8_dp)
         vc = ec*(1._dp + rs/(3._dp*(rs+7.8_dp)))
         vxc(i) = vx + vc
      END IF
    END DO

  END SUBROUTINE wigner_slater_functional
! *****************************************************************************
!> \brief ...
!> \param method ...
!> \param multiplicity ...
!> \retval consistent ...
! *****************************************************************************
  FUNCTION atom_consistent_method(method,multiplicity) RESULT(consistent)
    INTEGER, INTENT(IN)                      :: method, multiplicity
    LOGICAL                                  :: consistent

    SELECT CASE (method)
      CASE DEFAULT
        consistent=.FALSE.
      CASE (do_rks_atom)
        consistent=(multiplicity==-1)
      CASE (do_uks_atom)
        consistent=(multiplicity/=-1)
      CASE (do_rhf_atom)
        consistent=(multiplicity==-1)
      CASE (do_uhf_atom)
        consistent=(multiplicity/=-1)
      CASE (do_rohf_atom)
        consistent=.FALSE.
    END SELECT

  END FUNCTION atom_consistent_method
! *****************************************************************************

END MODULE atom_utils
