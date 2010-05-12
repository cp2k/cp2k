!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that build the integrals of the Vxc potential calculated
!>      for the atomic code
! *****************************************************************************
MODULE atom_output
  USE atom_types,                      ONLY: &
       CGTO_BASIS, GTH_PSEUDO, GTO_BASIS, NO_PSEUDO, NUM_BASIS, STO_BASIS, &
       atom_basis_type, atom_potential_type, atom_state, atom_type
  USE atom_utils,                      ONLY: get_maxl_occ,&
                                             get_maxn_occ
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_dkh0_atom, do_dkh1_atom, do_dkh2_atom, do_dkh3_atom, do_dkh4_atom, &
       do_dkh5_atom, do_nonrel_atom, do_rhf_atom, do_rks_atom, do_rohf_atom, &
       do_uhf_atom, do_uks_atom, do_zoramp_atom, xc_none
  USE input_cp2k_check,                ONLY: xc_functionals_expand
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: evolt
  USE xc_derivatives,                  ONLY: xc_functional_get_info
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_output'

  PUBLIC :: atom_print_state, atom_print_energies, atom_print_iteration, &
            atom_print_basis, atom_print_method, atom_print_info, atom_print_potential

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE atom_print_info(zval,info,iw,error)
    INTEGER, INTENT(IN)                      :: zval
    CHARACTER(len=*), INTENT(IN)             :: info
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_print_info', &
      routineP = moduleN//':'//routineN

    WRITE(iw,'(/," ",A,T40,A," [",A,"]",T62,"Atomic number:",T78,I3,/)') &
      ADJUSTL(TRIM(info)),TRIM(ptable(zval)%name),TRIM(ptable(zval)%symbol),zval

  END SUBROUTINE atom_print_info
! *****************************************************************************
  SUBROUTINE atom_print_state(state,iw,error)
    TYPE(atom_state)                         :: state
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_print_state', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1), DIMENSION(0:3)         :: label = (/ "S","P","D","F" /)
    INTEGER                                  :: j, l, mc, mlc, mlo, mm(0:3), &
                                                mo
    LOGICAL                                  :: failure = .FALSE.

    WRITE(iw,'(/,T2,A)') "Electronic structure"
    WRITE(iw,'(T5,A,T71,F10.2)') "Total number of core electrons",SUM(state%core)
    WRITE(iw,'(T5,A,T71,F10.2)') "Total number of valence electrons",SUM(state%occ)
    WRITE(iw,'(T5,A,T71,F10.2)') "Total number of electrons",SUM(state%occ+state%core)
    SELECT CASE (state%multiplicity)
      CASE (-1)
        WRITE(iw,'(T5,A,T68,A)') "Multiplicity","not specified"
      CASE (-2)
        WRITE(iw,'(T5,A,T72,A)') "Multiplicity","high spin"
      CASE (-3)
        WRITE(iw,'(T5,A,T73,A)') "Multiplicity","low spin"
      CASE (1)
        WRITE(iw,'(T5,A,T74,A)') "Multiplicity","singlet"
      CASE (2)
        WRITE(iw,'(T5,A,T74,A)') "Multiplicity","doublet"
      CASE (3)
        WRITE(iw,'(T5,A,T74,A)') "Multiplicity","triplet"
      CASE (4)
        WRITE(iw,'(T5,A,T74,A)') "Multiplicity","quartet"
      CASE (5)
        WRITE(iw,'(T5,A,T74,A)') "Multiplicity","quintet"
      CASE (6)
        WRITE(iw,'(T5,A,T75,A)') "Multiplicity","sextet"
      CASE (7)
        WRITE(iw,'(T5,A,T75,A)') "Multiplicity","septet"
      CASE DEFAULT
    END SELECT

    mlo = get_maxl_occ(state%occ)
    mlc = get_maxl_occ(state%core)
    mm  = get_maxn_occ(state%core)

    IF (state%multiplicity == -1) THEN
      DO l=0,MAX(mlo,mlc)
        mo = state%maxn_occ(l)
        IF ( SUM(state%core(l,:)) == 0 ) THEN
           WRITE(iw,'(A5,T10,10F6.2)') label(l),(state%occ(l,j),j=1,mo)
        ELSE
           mc = mm(l)
           CPPrecondition(SUM(state%occ(l,1:mc))==0,cp_failure_level,routineP,error,failure)
           WRITE(iw,ADVANCE="no",FMT='(A5,T9,A1,10F6.2)') label(l),"[",(state%core(l,j),j=1,mc)
           WRITE(iw,FMT='(A1,F5.2,10F6.2)') "]",(state%occ(l,j),j=mc+1,mc+mo)
        END IF
      END DO
    ELSE
      WRITE(iw,'(T5,A)') "Alpha Electrons"
      DO l=0,MAX(mlo,mlc)
        mo = state%maxn_occ(l)
        IF ( SUM(state%core(l,:)) == 0 ) THEN
           WRITE(iw,'(A5,T10,10F6.2)') label(l),(state%occa(l,j),j=1,mo)
        ELSE
           mc = mm(l)
           WRITE(iw,ADVANCE="no",FMT='(A5,T9,A1,10F6.2)') label(l),"[",(0.5_dp*state%core(l,j),j=1,mc)
           WRITE(iw,FMT='(A1,F5.2,10F6.2)') "]",(state%occa(l,j),j=1,mo)
        END IF
      END DO
      WRITE(iw,'(T5,A)') "Beta Electrons"
      DO l=0,MAX(mlo,mlc)
        mo = state%maxn_occ(l)
        IF ( SUM(state%core(l,:)) == 0 ) THEN
           WRITE(iw,'(A5,T10,10F6.2)') label(l),(state%occb(l,j),j=1,mo)
        ELSE
           mc = mm(l)
           WRITE(iw,ADVANCE="no",FMT='(A5,T9,A1,10F6.2)') label(l),"[",(0.5_dp*state%core(l,j),j=1,mc)
           WRITE(iw,FMT='(A1,F5.2,10F6.2)') "]",(state%occb(l,j),j=1,mo)
        END IF
      END DO
    END IF
    WRITE(iw,*)

  END SUBROUTINE atom_print_state
! *****************************************************************************
  SUBROUTINE atom_print_energies(atom,iw,error)
    TYPE(atom_type)                          :: atom
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_print_energies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, l, n

    WRITE(iw,'(/,A,T36,A,T61,F20.12)') " Energy components [Hartree]",&
                                   "    Total Energy ::",atom%energy%etot
    WRITE(iw,'(T36,A,T61,F20.12)') "     Band Energy ::",atom%energy%eband
    WRITE(iw,'(T36,A,T61,F20.12)') "  Kinetic Energy ::",atom%energy%ekin
    WRITE(iw,'(T36,A,T61,F20.12)') "Potential Energy ::",atom%energy%epot
    IF (atom%energy%ekin /= 0.0_dp) THEN
       WRITE(iw,'(T36,A,T61,F20.12)') "   Virial (-V/T) ::",-atom%energy%epot/atom%energy%ekin
    END IF
    WRITE(iw,'(T36,A,T61,F20.12)') "     Core Energy ::",atom%energy%ecore
    IF(atom%energy%exc /= 0._dp) &
      WRITE(iw,'(T36,A,T61,F20.12)') "       XC Energy ::",atom%energy%exc
    WRITE(iw,'(T36,A,T61,F20.12)') "  Coulomb Energy ::",atom%energy%ecoulomb
    IF(atom%energy%eexchange /= 0._dp) &
      WRITE(iw,'(T34,A,T61,F20.12)') "HF Exchange Energy ::",atom%energy%eexchange
    IF ( atom%potential%ppot_type /= NO_PSEUDO ) THEN
       WRITE(iw,'(T20,A,T61,F20.12)') "    Total Pseudopotential Energy ::",atom%energy%epseudo
       WRITE(iw,'(T20,A,T61,F20.12)') "    Local Pseudopotential Energy ::",atom%energy%eploc
       IF(atom%energy%elsd /= 0._dp) &
         WRITE(iw,'(T20,A,T61,F20.12)') "     Local Spin-potential Energy ::",atom%energy%elsd
       WRITE(iw,'(T20,A,T61,F20.12)') " Nonlocal Pseudopotential Energy ::",atom%energy%epnl
    END IF
    IF ( atom%potential%confinement ) THEN
       WRITE(iw,'(T36,A,T61,F20.12)') "     Confinement ::",atom%energy%econfinement
    END IF

    IF ( atom%state%multiplicity == -1 ) THEN
       WRITE(iw,'(/,A,T20,A,T30,A,T36,A,T49,A,T71,A,/)') " Orbital energies",&
            "State","L","Occupation","Energy[a.u.]","Energy[eV]"
       DO l=0,atom%state%maxl_calc
         n = atom%state%maxn_calc(l)
         DO i=1,n 
           WRITE(iw,'(T23,I2,T30,I1,T36,F10.3,T46,F15.6,T66,F15.6)') &
              i,l,atom%state%occupation(l,i),atom%orbitals%ener(i,l),atom%orbitals%ener(i,l)*evolt
         END DO
         IF (n>0) WRITE(iw,*)
       END DO
    ELSE
       WRITE(iw,'(/,A,T20,A,T30,A,T36,A,T42,A,T55,A,T71,A,/)') " Orbital energies",&
            "State","Spin","L","Occupation","Energy[a.u.]","Energy[eV]"
       DO l=0,atom%state%maxl_calc
         n = atom%state%maxn_calc(l)
         DO i=1,n 
           WRITE(iw,'(T23,I2,T29,A,T36,I1,T42,F10.3,T52,F15.6,T68,F13.6)') &
              i,"alpha",l,atom%state%occa(l,i),atom%orbitals%enera(i,l),atom%orbitals%enera(i,l)*evolt
         END DO
         DO i=1,n 
           WRITE(iw,'(T23,I2,T29,A,T36,I1,T42,F10.3,T52,F15.6,T68,F13.6)') &
              i," beta",l,atom%state%occb(l,i),atom%orbitals%enerb(i,l),atom%orbitals%enerb(i,l)*evolt
         END DO
         IF (n>0) WRITE(iw,*)
       END DO
    END IF

  END SUBROUTINE atom_print_energies
! *****************************************************************************
  SUBROUTINE atom_print_iteration(iter,deps,etot,iw,error)
    INTEGER, INTENT(IN)                      :: iter
    REAL(dp), INTENT(IN)                     :: deps, etot
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_print_iteration', &
      routineP = moduleN//':'//routineN

    IF (iter==1) THEN
      WRITE(iw,'(/," ",79("*"),/,T19,A,T38,A,T70,A,/," ",79("*"))') &
           "Iteration","Convergence","Energy [au]"
    END IF
    WRITE(iw,'(T20,i8,T34,G14.6,T61,F20.12)') iter,deps,etot

  END SUBROUTINE atom_print_iteration
! *****************************************************************************
  SUBROUTINE atom_print_basis(atom_basis,iw,title,error)
    TYPE(atom_basis_type)                    :: atom_basis
    INTEGER, INTENT(IN)                      :: iw
    CHARACTER(len=*)                         :: title
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_print_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, l
    LOGICAL                                  :: failure = .FALSE.

    WRITE(iw,'(/,A)') TRIM(title)
    SELECT CASE (atom_basis%basis_type)
      CASE (GTO_BASIS)
        IF ( atom_basis%geometrical ) THEN
          WRITE(iw,'(/," ",21("*"),A,22("*"))') " Geometrical Gaussian Type Orbitals "
          WRITE(iw,'(A,F15.8,T41,A,F15.8)') " Initial exponent: ",atom_basis%aval,&
                " Proportionality factor: ",atom_basis%cval
        ELSE
          WRITE(iw,'(/," ",21("*"),A,21("*"))') " Uncontracted Gaussian Type Orbitals "
        END IF
        IF ( atom_basis%nbas(0) > 0 ) THEN
          WRITE(iw,'(/,T2,A,(T30,I5,T51,F30.8))') &
           "s Exponents: ",(i,atom_basis%am(i,0),i=1,atom_basis%nbas(0))
        END IF
        IF ( atom_basis%nbas(1) > 0 ) THEN
          WRITE(iw,'(/,T2,A,(T30,I5,T51,F30.8))') &
           "p Exponents: ",(i,atom_basis%am(i,1),i=1,atom_basis%nbas(1))
        END IF
        IF ( atom_basis%nbas(2) > 0 ) THEN
          WRITE(iw,'(/,T2,A,(T30,I5,T51,F30.8))') &
           "d Exponents: ",(i,atom_basis%am(i,2),i=1,atom_basis%nbas(2))
        END IF
        IF ( atom_basis%nbas(3) > 0 ) THEN
          WRITE(iw,'(/,T2,A,(T30,I5,T51,F30.8))') &
           "f Exponents: ",(i,atom_basis%am(i,3),i=1,atom_basis%nbas(3))
        END IF
        WRITE(iw,'(" ",79("*"))')
      CASE (CGTO_BASIS)
        WRITE(iw,'(/," ",22("*"),A,22("*"))') " Contracted Gaussian Type Orbitals "
        DO l=0,3
          IF ( atom_basis%nbas(l) > 0 ) THEN
            IF (l==0) WRITE(iw,'(A)') " s Functions"
            IF (l==1) WRITE(iw,'(A)') " p Functions"
            IF (l==2) WRITE(iw,'(A)') " d Functions"
            IF (l==3) WRITE(iw,'(A)') " f Functions"
            DO i=1,atom_basis%nprim(l)
              WRITE(iw,'(F15.6,5(T21,6F10.6,/))') &
                  atom_basis%am(i,l),(atom_basis%cm(i,j,l),j=1,atom_basis%nbas(l))
            END DO
          END IF
        END DO
        WRITE(iw,'(" ",79("*"))')
      CASE (STO_BASIS)
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (NUM_BASIS)
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
    END SELECT

  END SUBROUTINE atom_print_basis
! *****************************************************************************
  SUBROUTINE atom_print_method(atom,iw,error)
    TYPE(atom_type)                          :: atom
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_print_method', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=160)                       :: reference
    CHARACTER(len=20)                        :: tmpStr
    INTEGER                                  :: ifun, il, meth, myfun, reltyp
    LOGICAL                                  :: failure = .FALSE., lsd
    TYPE(section_vals_type), POINTER         :: xc_fun, xc_fun_section, &
                                                xc_section

    NULLIFY(xc_fun,xc_fun_section,xc_section)

    meth = atom%method_type

    xc_section => atom%xc_section
    xc_fun_section => section_vals_get_subs_vals(xc_section,"XC_FUNCTIONAL",error=error)
    SELECT CASE (meth)
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rks_atom)
        CALL section_vals_val_get(xc_fun_section,"_SECTION_PARAMETERS_",i_val=myfun,error=error)
      CASE (do_uks_atom)
        CALL section_vals_val_get(xc_fun_section,"_SECTION_PARAMETERS_",i_val=myfun,error=error)
      CASE (do_rhf_atom)
        myfun = xc_none
      CASE (do_uhf_atom)
        myfun = xc_none
      CASE (do_rohf_atom)
        myfun = xc_none
    END SELECT

    SELECT CASE (meth)
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rks_atom)
        IF (iw > 0) WRITE (iw, fmt="(/,' METHOD    | Restricted Kohn-Sham Calculation')")
      CASE (do_uks_atom)
        IF (iw > 0) WRITE (iw, fmt="(/,' METHOD    | Unrestricted Kohn-Sham Calculation')")
      CASE (do_rhf_atom)
        IF (iw > 0) WRITE (iw, fmt="(/,' METHOD    | Restricted Hartree-Fock Calculation')")
      CASE (do_uhf_atom)
        IF (iw > 0) WRITE (iw, fmt="(/,' METHOD    | Unrestricted Hartree-Fock Calculation')")
      CASE (do_rohf_atom)
        IF (iw > 0) WRITE (iw, fmt="(/,' METHOD    | Restricted Open-Shell Kohn-Sham Calculation')")
    END SELECT

    IF (atom%pp_calc) THEN
      IF (iw > 0) WRITE (iw, fmt="(' METHOD    | Nonrelativistic Calculation')")
    ELSE
      reltyp = atom%relativistic

      SELECT CASE (reltyp)
        CASE DEFAULT
          CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
        CASE (do_nonrel_atom)
          IF (iw > 0) WRITE (iw, fmt="(' METHOD    | Nonrelativistic Calculation')")
        CASE (do_zoramp_atom)
          IF (iw > 0) WRITE (iw, fmt="(' METHOD    | Relativistic Calculation using ZORA(MP)')")
        CASE (do_dkh0_atom)
          IF (iw > 0) WRITE (iw, fmt="(' METHOD    | Relativistic Calculation using Douglas-Kroll 0th order')")
          IF (iw > 0) WRITE (iw, fmt="(' METHOD    | Relativistic Calculation using kietic energy scaling')")
        CASE (do_dkh1_atom)
          IF (iw > 0) WRITE (iw, fmt="(' METHOD    | Relativistic Calculation using Douglas-Kroll 1st order')")
          IF (iw > 0) WRITE (iw, fmt="(' METHOD    | Relativistic Calculation using Foldy-Wouthuysen transformation')")
        CASE (do_dkh2_atom)
          IF (iw > 0) WRITE (iw, fmt="(' METHOD    | Relativistic Calculation using Douglas-Kroll 2nd order')")
        CASE (do_dkh3_atom)
          IF (iw > 0) WRITE (iw, fmt="(' METHOD    | Relativistic Calculation using Douglas-Kroll 3rd order')")
        CASE (do_dkh4_atom)
          IF (iw > 0) WRITE (iw, fmt="(' METHOD    | Relativistic Calculation using Douglas-Kroll 4th order')")
        CASE (do_dkh5_atom)
          IF (iw > 0) WRITE (iw, fmt="(' METHOD    | Relativistic Calculation using Douglas-Kroll 5th order')")
      END SELECT
    END IF

    lsd = (meth==do_uks_atom)

    IF (myfun/=xc_none) THEN
      CALL section_vals_val_get(xc_section,"FUNCTIONAL_ROUTINE",c_val=tmpStr,error=error)
      IF (iw > 0) WRITE (iw, fmt="(' FUNCTIONAL| ROUTINE=',a)") TRIM(tmpStr)
      CALL xc_functionals_expand(xc_fun_section,xc_section,error=error)
      IF (iw > 0) THEN
        ifun=0
        DO
          ifun=ifun+1
          xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
          IF (.NOT.ASSOCIATED(xc_fun)) EXIT
          CALL xc_functional_get_info(xc_fun, lsd=lsd, reference=reference, error=error)
          WRITE (iw, fmt="(' FUNCTIONAL| ',a,':')") TRIM(xc_fun%section%name)
          il = LEN_TRIM(reference)
          IF (il < 67) THEN
            WRITE (iw, fmt="(' FUNCTIONAL| ',a)") reference(1:il)
          ELSE IF (il < 134) THEN
            WRITE (iw, fmt="(' FUNCTIONAL| ',a)") reference(1:67)
            WRITE (iw, fmt="(' FUNCTIONAL| ',a)") reference(68:il)
          ELSE
            WRITE (iw, fmt="(' FUNCTIONAL| ',a)") reference(1:67)
            WRITE (iw, fmt="(' FUNCTIONAL| ',a)") reference(68:134)
            WRITE (iw, fmt="(' FUNCTIONAL| ',a)") reference(134:il)
          END IF
        END DO
      END IF
    ELSE
      IF(iw>0) WRITE (iw, fmt="(' FUNCTIONAL| NO EXCHANGE-CORRELATION FUNCTIONAL USED.')")
    END IF

  END SUBROUTINE atom_print_method
! *****************************************************************************
  SUBROUTINE atom_print_potential(potential,iw,error)
    TYPE(atom_potential_type)                :: potential
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_print_potential', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=60)                        :: pline
    INTEGER                                  :: i, j, k, l
    LOGICAL                                  :: failure = .FALSE.

    SELECT CASE (potential%ppot_type)
      CASE (NO_PSEUDO)
        WRITE(iw,'(/," ",28("*"),A,27("*"))') " All Electron Potential "
      CASE (GTH_PSEUDO)
        WRITE(iw,'(/," ",29("*"),A,29("*"))') " GTH Pseudopotential "
        WRITE(iw,'(T10,A,T76,F5.1)') " Core Charge ",potential%gth_pot%zion
        WRITE(iw,'(T10,A,T66,F15.6)') " Rc ",potential%gth_pot%rc
        WRITE(pline,'(5F12.6)') (potential%gth_pot%cl(i),i=1,potential%gth_pot%ncl)
        WRITE(iw,'(T10,A,T21,A60)') " C1 C2 ... ",ADJUSTR(pline)
        IF(potential%gth_pot%lpotextended) THEN
          DO k=1,potential%gth_pot%nexp_lpot
            WRITE(iw,'(T10,A,F10.6,T38,A,4F10.6)') " LPot: rc=",potential%gth_pot%alpha_lpot(k),&
              "CX=",(potential%gth_pot%cval_lpot(i,k),i=1,potential%gth_pot%nct_lpot(k))
          END DO
        END IF
        IF(potential%gth_pot%nlcc) THEN
          DO k=1,potential%gth_pot%nexp_nlcc
            WRITE(iw,'(T10,A,F10.6,T38,A,4F10.6)') " LSDPot: rc=",potential%gth_pot%alpha_nlcc(k),&
              "CX=",(potential%gth_pot%cval_nlcc(i,k),i=1,potential%gth_pot%nct_nlcc(k))
          END DO
        END IF
        IF(potential%gth_pot%lsdpot) THEN
          DO k=1,potential%gth_pot%nexp_lsd
            WRITE(iw,'(T10,A,F10.6,T38,A,4F10.6)') " LSDPot: rc=",potential%gth_pot%alpha_lsd(k),&
              "CX=",(potential%gth_pot%cval_lsd(i,k),i=1,potential%gth_pot%nct_lsd(k))
          END DO
        END IF
        DO l=0,3
          IF (potential%gth_pot%nl(l) > 0) THEN
             WRITE(iw,'(T10,A,T76,I5)') " Angular momentum ",l
             WRITE(iw,'(T10,A,T66,F15.6)') " Rcnl ",potential%gth_pot%rcnl(l)
             WRITE(iw,'(T10,A,T76,I5)') " Nl ",potential%gth_pot%nl(l)
             WRITE(pline,'(5F12.6)') (potential%gth_pot%hnl(1,j,l),j=1,potential%gth_pot%nl(l))
             WRITE(iw,'(T10,A,T21,A60)') " Hnl ",ADJUSTR(pline)
             DO i=2,potential%gth_pot%nl(l)
               WRITE(pline,'(T21,5F12.6)') (potential%gth_pot%hnl(i,j,l),j=i,potential%gth_pot%nl(l))
               WRITE(iw,'(T21,A60)') ADJUSTR(pline)
             END DO
          END IF
        END DO
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
    END SELECT
    IF (potential%confinement) THEN
      WRITE(iw,'(/,T10,A,T51,F12.6," * (R /",F6.2,")**",I2)') &
        " Confinement Potential ",potential%acon,potential%rcon,potential%ncon
    END IF
    WRITE(iw,'(" ",79("*"))')

  END SUBROUTINE atom_print_potential
! *****************************************************************************

END MODULE atom_output

! *****************************************************************************
