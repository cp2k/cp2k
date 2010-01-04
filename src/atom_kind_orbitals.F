!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculate the orbitals for a given atomic kind type       
! *****************************************************************************
MODULE atom_kind_orbitals
  USE atom_electronic_structure,       ONLY: calculate_atom
  USE atom_operators,                  ONLY: atom_int_release,&
                                             atom_int_setup,&
                                             atom_ppint_release,&
                                             atom_ppint_setup,&
                                             atom_relint_release,&
                                             atom_relint_setup
  USE atom_types,                      ONLY: &
       CGTO_BASIS, atom_basis_type, atom_integrals, atom_orbitals, &
       atom_potential_type, atom_type, create_atom_orbs, create_atom_type, &
       release_atom_basis, release_atom_potential, release_atom_type, set_atom
  USE atom_utils,                      ONLY: get_maxl_occ,&
                                             get_maxn_occ
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE input_constants,                 ONLY: do_analytic,&
                                             do_gapw_log,&
                                             do_nonrel_atom,&
                                             do_numeric,&
                                             do_rks_atom,&
                                             gth_pseudo,&
                                             no_pseudo
  USE kinds,                           ONLY: dp
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: bohr
  USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                             create_grid_atom,&
                                             grid_atom_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_kind_orbitals'

  PUBLIC :: calculate_atomic_orbitals

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE calculate_atomic_orbitals (atomic_kind,iunit,pmat,ispin,nspin,error)
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    INTEGER, INTENT(IN), OPTIONAL            :: iunit
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: pmat
    INTEGER, INTENT(IN)                      :: ispin, nspin
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_atomic_orbitals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, ii, ipgf, j, k, k1, &
                                                k2, l, ll, lm, m, mb, mo, n, &
                                                ne, ngp, nj, nn, nr, ns, &
                                                nset, nsgf, quadtype, z
    INTEGER, DIMENSION(0:3)                  :: econfx
    INTEGER, DIMENSION(0:3, 10)              :: ncalc, ncore, nelem
    INTEGER, DIMENSION(0:3, 100)             :: set_index, shell_index
    INTEGER, DIMENSION(:), POINTER           :: econf, lmax, lmin, npgf, &
                                                nppnl, nshell, ppeconf
    INTEGER, DIMENSION(:, :), POINTER        :: addel, first_sgf, laddel, &
                                                last_sgf, ls, naddel
    LOGICAL                                  :: bs_occupation, &
                                                failure = .FALSE.
    REAL(KIND=dp)                            :: ac, al, ear, rk, scal, zeff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ap, ce
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc, hp
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atom_basis_type), POINTER           :: basis
    TYPE(atom_integrals), POINTER            :: integrals
    TYPE(atom_orbitals), POINTER             :: orbitals
    TYPE(atom_potential_type), POINTER       :: potential
    TYPE(atom_type), POINTER                 :: atom
    TYPE(grid_atom_type), POINTER            :: grid
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

    NULLIFY(atom)
    CALL create_atom_type(atom,error)

    CALL get_atomic_kind(atomic_kind,z=z,zeff=zeff)
    NULLIFY(all_potential,gth_potential,orb_basis_set)
    CALL get_atomic_kind(atomic_kind=atomic_kind,&
                         all_potential=all_potential,&
                         gth_potential=gth_potential,&
                         orb_basis_set=orb_basis_set,&
                         bs_occupation=bs_occupation,&
                         addel=addel,laddel=laddel,naddel=naddel)
    CPPostcondition(ASSOCIATED(orb_basis_set), cp_failure_level, routineP, error, failure)

    atom%z = z
    CALL set_atom(atom,&
                  pp_calc=ASSOCIATED(gth_potential),&
                  method_type=do_rks_atom,&
                  relativistic=do_nonrel_atom,&
                  coulomb_integral_type=do_numeric,&
                  exchange_integral_type=do_numeric,&
                  error=error)

    ALLOCATE (potential,basis,integrals,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    IF ( ASSOCIATED(gth_potential) ) THEN
      potential%ppot_type=gth_pseudo
      potential%confinement=.TRUE. 
      potential%acon=0.1_dp
      potential%rcon=2._dp*ptable(z)%vdw_radius*bohr
      potential%ncon=2

      CALL get_potential(gth_potential,&
                         zeff=zeff,&
                         elec_conf=ppeconf,&
                         alpha_core_charge=ac,&
                         nexp_ppl=ne,&
                         cexp_ppl=ce,&
                         lppnl=lm,&
                         nprj_ppnl=nppnl,&
                         alpha_ppnl=ap,&
                         hprj_ppnl=hp)

      potential%gth_pot%zion      = zeff
      potential%gth_pot%rc        = SQRT(0.5_dp/ac)
      potential%gth_pot%ncl       = ne
      potential%gth_pot%cl(:)     = 0._dp
      IF (ac > 0._dp) THEN
        DO i=1,ne
          potential%gth_pot%cl(i) = ce(i)/(2._dp*ac)**(i-1) 
        END DO
      END IF
      potential%gth_pot%nl(:)     = 0 
      potential%gth_pot%rcnl(:)   = 0._dp
      potential%gth_pot%hnl(:,:,:)= 0._dp
      DO l=0,lm
        n = nppnl(l)
        potential%gth_pot%nl(l)     = n
        potential%gth_pot%rcnl(l)   = SQRT(0.5_dp/ap(l))
        potential%gth_pot%hnl(1:n,1:n,l)= hp(1:n,1:n,l)
      END DO

      CALL set_atom(atom,zcore=NINT(zeff),potential=potential,error=error)
    ELSE
      potential%ppot_type=no_pseudo
      potential%confinement=.FALSE.
      CALL set_atom(atom,zcore=z,potential=potential,error=error)
    END IF

    CALL get_gto_basis_set(orb_basis_set,&
      nset=nset,nshell=nshell,npgf=npgf,lmin=lmin,lmax=lmax,l=ls,nsgf=nsgf,zet=zet,gcc=gcc,&
      first_sgf=first_sgf,last_sgf=last_sgf)

    NULLIFY(grid)
    ngp = 400
    quadtype = do_gapw_log 
    CALL allocate_grid_atom(grid,error)
    CALL create_grid_atom(grid,ngp,1,1,quadtype)
    grid%nr = ngp
    basis%grid => grid

    NULLIFY(basis%am,basis%cm,basis%as,basis%bf,basis%dbf)
    basis%basis_type = CGTO_BASIS
    basis%eps_eig = 1.e-12_dp

    ! fill in the basis data structures
    set_index = 0
    shell_index = 0
    basis%nprim = 0
    basis%nbas = 0
    DO i=1,nset
      DO j=lmin(i),MIN(lmax(i),3)
        basis%nprim(j) = basis%nprim(j) + npgf(i)
      END DO
      DO j=1,nshell(i)
        l = ls(j,i)
        IF ( l <= 3 ) THEN
          basis%nbas(l) = basis%nbas(l)+1
          k = basis%nbas(l)
          CPPostcondition(k<=100, cp_failure_level, routineP, error, failure)
          set_index(l,k) = i
          shell_index(l,k) = j
        END IF
      END DO
    END DO

    nj = MAXVAL(basis%nprim)
    ns = MAXVAL(basis%nbas)
    ALLOCATE (basis%am(nj,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    basis%am = 0._dp
    ALLOCATE (basis%cm(nj,ns,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    basis%cm = 0._dp
    DO j=0,3
      nj = 0
      ns = 0
      DO i=1,nset
        IF (j >= lmin(i) .AND. j <= lmax(i)) THEN
          DO ipgf = 1,npgf(i)
            basis%am(nj+ipgf,j) = zet(ipgf,i)
          END DO
          DO ii=1,nshell(i)
            IF ( ls(ii,i) == j ) THEN
              ns = ns + 1
              DO ipgf=1,npgf(i)
                basis%cm(nj+ipgf,ns,j) = gcc(ipgf,ii,i)
              END DO
            END IF
          END DO
          nj = nj + npgf(i)
        END IF
      END DO
    END DO

    ! initialize basis function on a radial grid
    nr = basis%grid%nr
    m  = MAXVAL(basis%nbas)
    ALLOCATE (basis%bf(nr,m,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    ALLOCATE (basis%dbf(nr,m,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    basis%bf =  0._dp
    basis%dbf = 0._dp
    DO l=0,3
      DO i=1,basis%nprim(l)
        al  = basis%am(i,l)
        DO k=1,nr
          rk  = basis%grid%rad(k)
          ear = EXP(-al*basis%grid%rad(k)**2)
          DO j=1,basis%nbas(l)
            basis%bf(k,j,l) = basis%bf(k,j,l) + rk**l * ear*basis%cm(i,j,l)
            basis%dbf(k,j,l) = basis%dbf(k,j,l) &
              + ( REAL(l,dp)*rk**(l-1) - 2._dp*al*rk**(l+1) ) * ear*basis%cm(i,j,l)
          END DO
        END DO
      END DO
    END DO

    CALL set_atom(atom,basis=basis,error=error)

    ! optimization defaults
    atom%optimization%damping    = 0.2_dp
    atom%optimization%eps_scf    = 1.e-6_dp
    atom%optimization%eps_diis   = 100._dp
    atom%optimization%max_iter   = 50
    atom%optimization%n_diis     = 5

    ! electronic state
    nelem = 0
    ncore = 0
    ncalc = 0
    NULLIFY(econf)
    IF ( ASSOCIATED(gth_potential) ) THEN
      CALL get_potential(gth_potential,elec_conf=econf)
      econfx = 0
      econfx(0:SIZE(econf)-1) = econf
      IF ( SUM(econf) >= 0 ) THEN
        DO l=0,3
          ll = 2*(2*l+1)
          nn = ptable(z)%e_conv(l)-econfx(l)
          ii = 0
          DO
            ii = ii + 1
            IF(nn <= ll) THEN
              ncore(l,ii) = nn
              EXIT
            ELSE
              ncore(l,ii) = ll
              nn = nn - ll
            END IF
          END DO
        END DO
        DO l=0,3
          ll = 2*(2*l+1)
          nn = ptable(z)%e_conv(l)
          ii = 0
          DO
            ii = ii + 1
            IF(nn <= ll) THEN
              nelem(l,ii) = nn
              EXIT
            ELSE
              nelem(l,ii) = ll
              nn = nn - ll
            END IF
          END DO
        END DO
        ncalc = nelem - ncore
      ELSE
        ncore = 0
        ncalc = 0
        DO l=0,3
          ll = 2*(2*l+1)
          nn = ABS(econfx(l))
          ii = 0
          DO
            ii = ii + 1
            IF(nn <= ll) THEN
              ncalc(l,ii) = -nn
              EXIT
            ELSE
              ncalc(l,ii) = -ll
              nn = nn - ll
            END IF
          END DO
        END DO
      END IF
    ELSE
      DO l=0,3
        ll = 2*(2*l+1)
        nn = ptable(z)%e_conv(l)
        ii = 0
        DO
          ii = ii + 1
          IF(nn <= ll) THEN
            nelem(l,ii) = nn
            EXIT
          ELSE
            nelem(l,ii) = ll
            nn = nn - ll
          END IF
        END DO
      END DO
      ncalc = nelem - ncore
    END IF

    IF(bs_occupation) THEN
      ! readjust the occupation number of the atomic orbitals
      ! according to the changes required from input in order to bias the initial guess
      DO i=1,SIZE(addel,1)
        ne=addel(i,ispin)
        l=laddel(i,ispin)
        nn=naddel(i,ispin)-l
        IF(ne/=0) THEN
          IF(nn==0) THEN
            DO ii=SIZE(nelem,2),1,-1
              IF(ncalc(l,ii)>0) THEN
                 IF((ncalc(l,ii)+ne) < 2*(2*l+1)+1) THEN
                   ncalc(l,ii) = ncalc(l,ii)+ne
                   nn = ii
                 ELSE
                   ncalc(l,ii+1) = ncalc(l,ii+1)+ne
                   nn = ii + 1 
                 END IF
                 EXIT
              ELSE IF (ii==1) THEN
                 ncalc(l,ii)=ncalc(l,ii)+ne
                 nn = ii 
              END IF
            END DO
          ELSE
            ncalc(l,nn) = ncalc(l,nn) + ne
          END IF
          IF( ncalc(l,nn)<0 ) THEN
             ncalc(l,nn) = 0
          END IF
        END IF
      END DO
    END IF

    IF ( atomic_kind%ghost ) THEN
      nelem = 0
      ncore = 0
      ncalc = 0
    END IF

    ALLOCATE (atom%state,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    atom%state%core = 0._dp
    atom%state%core(0:3,1:7) = REAL(ncore(0:3,1:7),dp)
    atom%state%occ = 0._dp
    atom%state%occ(0:3,1:7)  = REAL(ncalc(0:3,1:7),dp)
    atom%state%occupation = 0._dp
    atom%state%multiplicity = -1
    DO l=0,3
      k = 0
      DO i=1,7
        IF ( ncalc(l,i) > 0 ) THEN
          k = k + 1
          atom%state%occupation(l,k) = REAL(ncalc(l,i),dp)
        END IF
      END DO
    END DO


    atom%state%maxl_occ = get_maxl_occ(atom%state%occupation)
    atom%state%maxn_occ = get_maxn_occ(atom%state%occupation)
    atom%state%maxl_calc = atom%state%maxl_occ
    atom%state%maxn_calc = atom%basis%nbas

    ! calculate integrals
    ! general integrals
    CALL atom_int_setup(integrals,basis,potential=atom%potential,&
                        eri_coulomb=(atom%coulomb_integral_type==do_analytic),&
                        eri_exchange=(atom%exchange_integral_type==do_analytic),error=error)
    ! potential
    CALL atom_ppint_setup(integrals,basis,potential=atom%potential,error=error)
    ! relativistic correction terms
    NULLIFY(integrals%tzora,integrals%hdkh)
    CALL atom_relint_setup(integrals,basis,atom%relativistic,zcore=REAL(atom%zcore,dp),error=error)
    CALL set_atom(atom,integrals=integrals,error=error)

    NULLIFY(orbitals)
    mo = MAXVAL(atom%state%maxn_calc)
    mb = MAXVAL(atom%basis%nbas)
    CALL create_atom_orbs(orbitals,mb,mo,error)
    CALL set_atom(atom,orbitals=orbitals,error=error)

    IF(PRESENT(iunit)) THEN
      CALL calculate_atom(atom,iunit,error=error)
    ELSE
      CALL calculate_atom(atom,-1,error=error)
    END IF

    IF (PRESENT(pmat)) THEN
      ! recover density matrix in CP2K/GPW order and normalization
      IF(ASSOCIATED(pmat)) THEN
        DEALLOCATE (pmat,STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      END IF
      ALLOCATE (pmat(nsgf,nsgf),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      pmat = 0._dp
      DO l=0,3
         ll = 2*l
         DO k1=1,atom%basis%nbas(l)
           DO k2=1,atom%basis%nbas(l)
             scal=SQRT(atom%integrals%ovlp(k1,k1,l)*atom%integrals%ovlp(k2,k2,l))/REAL(2*l+1,KIND=dp)
             i=first_sgf(shell_index(l,k1),set_index(l,k1))
             j=first_sgf(shell_index(l,k2),set_index(l,k2))
             DO m=0,ll
               pmat(i+m,j+m) = atom%orbitals%pmat(k1,k2,l)*scal
             END DO
           END DO
         END DO
      ENDDO
    END IF

    ! clean up
    CALL atom_int_release(integrals,error)
    CALL atom_ppint_release(integrals,error)
    CALL atom_relint_release(integrals,error)
    CALL release_atom_basis(basis,error)
    CALL release_atom_potential(potential,error)
    CALL release_atom_type(atom,error)

    DEALLOCATE (potential,basis,integrals,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE calculate_atomic_orbitals

! *****************************************************************************

END MODULE atom_kind_orbitals
