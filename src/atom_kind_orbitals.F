!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculate the orbitals for a given atomic kind type
! *****************************************************************************
MODULE atom_kind_orbitals
  USE ai_onecenter,                    ONLY: sg_erfc
  USE atom_electronic_structure,       ONLY: calculate_atom
  USE atom_fit,                        ONLY: atom_fit_density
  USE atom_operators,                  ONLY: atom_int_release,&
                                             atom_int_setup,&
                                             atom_ppint_release,&
                                             atom_ppint_setup,&
                                             atom_relint_release,&
                                             atom_relint_setup
  USE atom_types,                      ONLY: &
       CGTO_BASIS, Clementi_geobas, GTO_BASIS, atom_basis_type, &
       atom_integrals, atom_orbitals, atom_potential_type, atom_type, &
       create_atom_orbs, create_atom_type, release_atom_basis, &
       release_atom_potential, release_atom_type, set_atom
  USE atom_utils,                      ONLY: get_maxl_occ,&
                                             get_maxn_occ
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE input_constants,                 ONLY: &
       do_analytic, do_dkh0_atom, do_dkh1_atom, do_dkh2_atom, do_dkh3_atom, &
       do_dkh4_atom, do_dkh5_atom, do_gapw_log, do_nonrel_atom, do_numeric, &
       do_rks_atom, do_sczoramp_atom, do_zoramp_atom, gth_pseudo, no_pseudo, &
       rel_dkh, rel_none, rel_sczora_mp, rel_zora, rel_zora_full, rel_zora_mp
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: dfac,&
                                             pi
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: bohr
  USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                             create_grid_atom,&
                                             grid_atom_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             init_atom_electronic_state,&
                                             qs_kind_type
  USE rel_control_types,               ONLY: rel_control_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_kind_orbitals'

  PUBLIC :: calculate_atomic_orbitals, calculate_atomic_density, &
            calculate_atomic_relkin

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param atomic_kind ...
!> \param qs_kind ...
!> \param iunit ...
!> \param pmat ...
!> \param ispin ...
!> \param confine ...
!> \param xc_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_atomic_orbitals (atomic_kind,qs_kind,iunit,pmat,ispin,confine,&
                                        xc_section,error)
    TYPE(atomic_kind_type), INTENT(IN)       :: atomic_kind
    TYPE(qs_kind_type), INTENT(IN)           :: qs_kind
    INTEGER, INTENT(IN), OPTIONAL            :: iunit
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: pmat
    INTEGER, INTENT(IN)                      :: ispin
    LOGICAL, INTENT(IN), OPTIONAL            :: confine
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: xc_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_atomic_orbitals', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, ierr, ii, ipgf, j, k, k1, k2, l, ll, lm, m, mb, mo, n, ne, &
      nexp_lpot, nexp_lsd, nexp_nlcc, ngp, nj, nr, ns, nset, nsgf, quadtype, z
    INTEGER, DIMENSION(0:3, 10)              :: ncalc, ncore, nelem
    INTEGER, DIMENSION(0:3, 100)             :: set_index, shell_index
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, nct_lpot, &
                                                nct_lsd, nct_nlcc, npgf, &
                                                nppnl, nshell, ppeconf
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgf, last_sgf, ls
    LOGICAL                                  :: failure, lpot_present, &
                                                lsd_present, nlcc_present
    REAL(KIND=dp)                            :: ac, al, ear, rk, scal, zeff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_lpot, alpha_lsd, &
                                                alpha_nlcc, ap, ce
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cval_lpot, cval_lsd, &
                                                cval_nlcc, zet
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

!INTEGER, DIMENSION(0:3)                  :: econfx
!INTEGER, DIMENSION(:), POINTER           :: econf
!INTEGER, DIMENSION(:, :), POINTER        :: addel, laddel, naddel 
!LOGICAL                                  :: bs_occupation

    failure = .FALSE.
    NULLIFY(atom)
    CALL create_atom_type(atom,error)

    IF(PRESENT(xc_section)) THEN
       atom%xc_section => xc_section
    ELSE
       NULLIFY(atom%xc_section)
    END IF

    CALL get_atomic_kind(atomic_kind,z=z)
    NULLIFY(all_potential,gth_potential,orb_basis_set)
    CALL get_qs_kind(qs_kind,zeff=zeff,&
                     basis_set=orb_basis_set,&
                     all_potential=all_potential,&
                     gth_potential=gth_potential, error=error)

    !!RZK CALL get_qs_kind(qs_kind,&
    !!RZK                  bs_occupation=bs_occupation,&
    !!RZK                  addel=addel,laddel=laddel,naddel=naddel)

    CPPostcondition(ASSOCIATED(orb_basis_set), cp_failure_level, routineP, error, failure)

    atom%z = z
    CALL set_atom(atom,&
                  pp_calc=ASSOCIATED(gth_potential),&
                  do_zmp=.FALSE., &
                  doread=.FALSE., &
                  read_vxc=.FALSE., &
                  method_type=do_rks_atom,&
                  relativistic=do_nonrel_atom,&
                  coulomb_integral_type=do_numeric,&
                  exchange_integral_type=do_numeric,&
                  error=error)

    ALLOCATE (potential,basis,integrals,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    IF ( ASSOCIATED(gth_potential) ) THEN
      potential%ppot_type=gth_pseudo
      IF ( PRESENT(confine) ) THEN
         potential%confinement=confine
      ELSE
         potential%confinement=.TRUE.
      END IF
      potential%acon=0.1_dp
      potential%rcon=2.0_dp*ptable(z)%vdw_radius*bohr
      potential%scon=2.0_dp

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
      !extended type
      potential%gth_pot%lpotextended = .FALSE.
      potential%gth_pot%lsdpot = .FALSE.
      potential%gth_pot%nlcc = .FALSE.
      potential%gth_pot%nexp_lpot = 0
      potential%gth_pot%nexp_lsd = 0
      potential%gth_pot%nexp_nlcc = 0
      CALL get_potential(gth_potential,&
                         lpot_present=lpot_present,&
                         lsd_present=lsd_present,&
                         nlcc_present=nlcc_present)
      IF (lpot_present) THEN
        CALL get_potential(gth_potential,&
                           nexp_lpot=nexp_lpot,&
                           alpha_lpot=alpha_lpot,&
                           nct_lpot=nct_lpot,&
                           cval_lpot=cval_lpot)
        potential%gth_pot%lpotextended = .TRUE.
        potential%gth_pot%nexp_lpot = nexp_lpot
        potential%gth_pot%alpha_lpot(1:nexp_lpot) = SQRT(0.5_dp/alpha_lpot(1:nexp_lpot))
        potential%gth_pot%nct_lpot(1:nexp_lpot) = nct_lpot(1:nexp_lpot)
        DO j=1,nexp_lpot
          ac = alpha_lpot(j)
          DO i=1,4
            potential%gth_pot%cval_lpot(i,j) = cval_lpot(i,j)/(2._dp*ac)**(i-1)
          END DO
        END DO
      END IF
      IF (lsd_present) THEN
        CALL get_potential(gth_potential,&
                           nexp_lsd=nexp_lsd,&
                           alpha_lsd=alpha_lsd,&
                           nct_lsd=nct_lsd,&
                           cval_lsd=cval_lsd)
        potential%gth_pot%lsdpot = .TRUE.
        potential%gth_pot%nexp_lsd = nexp_lsd
        potential%gth_pot%alpha_lsd(1:nexp_lsd) = SQRT(0.5_dp/alpha_lsd(1:nexp_lsd))
        potential%gth_pot%nct_lsd(1:nexp_lsd) = nct_lsd(1:nexp_lsd)
        DO j=1,nexp_lpot
          ac = alpha_lsd(j)
          DO i=1,4
            potential%gth_pot%cval_lsd(i,j) = cval_lsd(i,j)/(2._dp*ac)**(i-1)
          END DO
        END DO
      END IF

      ! nonlocal part
      potential%gth_pot%nl(:)     = 0
      potential%gth_pot%rcnl(:)   = 0._dp
      potential%gth_pot%hnl(:,:,:)= 0._dp
      DO l=0,lm
        n = nppnl(l)
        potential%gth_pot%nl(l)     = n
        potential%gth_pot%rcnl(l)   = SQRT(0.5_dp/ap(l))
        potential%gth_pot%hnl(1:n,1:n,l)= hp(1:n,1:n,l)
      END DO

      IF (nlcc_present) THEN
        CALL get_potential(gth_potential,&
                           nexp_nlcc=nexp_nlcc,&
                           alpha_nlcc=alpha_nlcc,&
                           nct_nlcc=nct_nlcc,&
                           cval_nlcc=cval_nlcc)
        potential%gth_pot%nlcc = .TRUE.
        potential%gth_pot%nexp_nlcc = nexp_nlcc
        potential%gth_pot%alpha_nlcc(1:nexp_nlcc) = alpha_nlcc(1:nexp_nlcc)
        potential%gth_pot%nct_nlcc(1:nexp_nlcc) = nct_nlcc(1:nexp_nlcc)
        potential%gth_pot%cval_nlcc(1:4,1:nexp_nlcc) = cval_nlcc(1:4,1:nexp_nlcc)
      END IF

      CALL set_atom(atom,zcore=NINT(zeff),potential=potential,error=error)
    ELSE
      potential%ppot_type=no_pseudo
      IF ( PRESENT(confine) ) THEN
         potential%confinement=confine
      ELSE
         potential%confinement=.FALSE.
      END IF
      potential%acon=0.1_dp
      potential%rcon=2.0_dp*ptable(z)%vdw_radius*bohr
      potential%scon=2.0_dp
      CALL set_atom(atom,zcore=z,potential=potential,error=error)
    END IF

    CALL get_gto_basis_set(orb_basis_set,&
      nset=nset,nshell=nshell,npgf=npgf,lmin=lmin,lmax=lmax,l=ls,nsgf=nsgf,zet=zet,gcc=gcc,&
      first_sgf=first_sgf,last_sgf=last_sgf)

    NULLIFY(grid)
    ngp = 400
    quadtype = do_gapw_log
    CALL allocate_grid_atom(grid,error)
    CALL create_grid_atom(grid,ngp,1,1,0,quadtype)
    grid%nr = ngp
    basis%grid => grid

    NULLIFY(basis%am,basis%cm,basis%as,basis%ns,basis%bf,basis%dbf)
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

    ! set up the electronic state
    CALL init_atom_electronic_state (atomic_kind=atomic_kind, &
                                     qs_kind=qs_kind, &
                                     ispin=ispin, & 
                                     ncalc=ncalc, &
                                     ncore=ncore, &
                                     nelem=nelem, &
                                     error=error)

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
    atom%state%maxn_calc = atom%state%maxn_occ

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

! *****************************************************************************
!> \brief ...
!> \param density ...
!> \param atomic_kind ...
!> \param qs_kind ...
!> \param ngto ...
!> \param iunit ...
!> \param allelectron ...
!> \param confine ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_atomic_density (density,atomic_kind,qs_kind,ngto,iunit,allelectron,confine,error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: density
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(qs_kind_type), POINTER              :: qs_kind
    INTEGER, INTENT(IN)                      :: ngto
    INTEGER, INTENT(IN), OPTIONAL            :: iunit
    LOGICAL, INTENT(IN), OPTIONAL            :: allelectron, confine
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_atomic_density', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: num_gto = 40

    INTEGER :: i, ierr, ii, j, k, l, ll, lm, m, mb, mo, n, ne, nexp_lpot, &
      nexp_lsd, nexp_nlcc, ngp, nn, nr, quadtype, relativistic, z
    INTEGER, DIMENSION(0:3)                  :: econfx, starti
    INTEGER, DIMENSION(0:3, 10)              :: ncalc, ncore, nelem
    INTEGER, DIMENSION(:), POINTER           :: econf, nct_lpot, nct_lsd, &
                                                nct_nlcc, nppnl, ppeconf
    LOGICAL                                  :: failure, lpot_present, &
                                                lsd_present, nlcc_present
    REAL(KIND=dp)                            :: ac, al, aval, cc, cval, ear, &
                                                rk, xx, zeff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: alpha_lpot, alpha_lsd, &
                                                alpha_nlcc, ap, ce
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cval_lpot, cval_lsd, cval_nlcc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: hp
    REAL(KIND=dp), DIMENSION(num_gto+2)      :: results
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atom_basis_type), POINTER           :: basis
    TYPE(atom_integrals), POINTER            :: integrals
    TYPE(atom_orbitals), POINTER             :: orbitals
    TYPE(atom_potential_type), POINTER       :: potential
    TYPE(atom_type), POINTER                 :: atom
    TYPE(grid_atom_type), POINTER            :: grid
    TYPE(gth_potential_type), POINTER        :: gth_potential

    failure = .FALSE.
    NULLIFY(atom)
    CALL create_atom_type(atom,error)

    CALL get_atomic_kind(atomic_kind,z=z)
    NULLIFY(all_potential, gth_potential)
    CALL get_qs_kind(qs_kind,zeff=zeff,&
                     all_potential=all_potential,&
                     gth_potential=gth_potential, error=error)

    IF(PRESENT(allelectron)) THEN
       IF(allelectron) THEN
          NULLIFY(gth_potential)
          zeff = z
       END IF
    END IF

    CPPrecondition(ngto<=num_gto, cp_failure_level, routineP, error, failure)

    IF ( ASSOCIATED(gth_potential) ) THEN
       ! PP calculation are non-relativistic
       relativistic = do_nonrel_atom
    ELSE
       ! AE calculations use DKH2
       relativistic = do_dkh2_atom
    END IF

    atom%z = z
    CALL set_atom(atom,&
                  pp_calc=ASSOCIATED(gth_potential),&
                  method_type=do_rks_atom,&
                  relativistic=relativistic,&
                  coulomb_integral_type=do_numeric,&
                  exchange_integral_type=do_numeric,&
                  error=error)

    ALLOCATE (potential,basis,integrals,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    IF ( ASSOCIATED(gth_potential) ) THEN
      potential%ppot_type=gth_pseudo
      IF ( PRESENT(confine) ) THEN
         potential%confinement=confine
      ELSE
         potential%confinement=.TRUE.
      END IF
      potential%acon=0.1_dp
      potential%rcon=2.0_dp*ptable(z)%vdw_radius*bohr
      potential%scon=2.0_dp

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
      !extended type
      potential%gth_pot%lpotextended = .FALSE.
      potential%gth_pot%lsdpot = .FALSE.
      potential%gth_pot%nlcc = .FALSE.
      potential%gth_pot%nexp_lpot = 0
      potential%gth_pot%nexp_lsd = 0
      potential%gth_pot%nexp_nlcc = 0
      CALL get_potential(gth_potential,&
                         lpot_present=lpot_present,&
                         lsd_present=lsd_present,&
                         nlcc_present=nlcc_present)
      IF (lpot_present) THEN
        CALL get_potential(gth_potential,&
                           nexp_lpot=nexp_lpot,&
                           alpha_lpot=alpha_lpot,&
                           nct_lpot=nct_lpot,&
                           cval_lpot=cval_lpot)
        potential%gth_pot%lpotextended = .TRUE.
        potential%gth_pot%nexp_lpot = nexp_lpot
        potential%gth_pot%alpha_lpot(1:nexp_lpot) = SQRT(0.5_dp/alpha_lpot(1:nexp_lpot))
        potential%gth_pot%nct_lpot(1:nexp_lpot) = nct_lpot(1:nexp_lpot)
        DO j=1,nexp_lpot
          ac = alpha_lpot(j)
          DO i=1,4
            potential%gth_pot%cval_lpot(i,j) = cval_lpot(i,j)/(2._dp*ac)**(i-1)
          END DO
        END DO
      END IF
      IF (lsd_present) THEN
        CALL get_potential(gth_potential,&
                           nexp_lsd=nexp_lsd,&
                           alpha_lsd=alpha_lsd,&
                           nct_lsd=nct_lsd,&
                           cval_lsd=cval_lsd)
        potential%gth_pot%lsdpot = .TRUE.
        potential%gth_pot%nexp_lsd = nexp_lsd
        potential%gth_pot%alpha_lsd(1:nexp_lsd) = SQRT(0.5_dp/alpha_lsd(1:nexp_lsd))
        potential%gth_pot%nct_lsd(1:nexp_lsd) = nct_lsd(1:nexp_lsd)
        DO j=1,nexp_lpot
          ac = alpha_lsd(j)
          DO i=1,4
            potential%gth_pot%cval_lsd(i,j) = cval_lsd(i,j)/(2._dp*ac)**(i-1)
          END DO
        END DO
      END IF

      ! nonlocal part
      potential%gth_pot%nl(:)     = 0
      potential%gth_pot%rcnl(:)   = 0._dp
      potential%gth_pot%hnl(:,:,:)= 0._dp
      DO l=0,lm
        n = nppnl(l)
        potential%gth_pot%nl(l)     = n
        potential%gth_pot%rcnl(l)   = SQRT(0.5_dp/ap(l))
        potential%gth_pot%hnl(1:n,1:n,l)= hp(1:n,1:n,l)
      END DO

      IF (nlcc_present) THEN
        CALL get_potential(gth_potential,&
                           nexp_nlcc=nexp_nlcc,&
                           alpha_nlcc=alpha_nlcc,&
                           nct_nlcc=nct_nlcc,&
                           cval_nlcc=cval_nlcc)
        potential%gth_pot%nlcc = .TRUE.
        potential%gth_pot%nexp_nlcc = nexp_nlcc
        potential%gth_pot%alpha_nlcc(1:nexp_nlcc) = alpha_nlcc(1:nexp_nlcc)
        potential%gth_pot%nct_nlcc(1:nexp_nlcc) = nct_nlcc(1:nexp_nlcc)
        potential%gth_pot%cval_nlcc(1:4,1:nexp_nlcc) = cval_nlcc(1:4,1:nexp_nlcc)
      END IF

      CALL set_atom(atom,zcore=NINT(zeff),potential=potential,error=error)
    ELSE
      potential%ppot_type=no_pseudo
      IF ( PRESENT(confine) ) THEN
         potential%confinement=confine
      ELSE
         potential%confinement=.FALSE.
      END IF
      potential%acon=0.1_dp
      potential%rcon=2.0_dp*ptable(z)%vdw_radius*bohr
      potential%scon=2.0_dp
      CALL set_atom(atom,zcore=z,potential=potential,error=error)
    END IF

    ! atomic grid
    NULLIFY(grid)
    ngp = 400
    quadtype = do_gapw_log
    CALL allocate_grid_atom(grid,error)
    CALL create_grid_atom(grid,ngp,1,1,0,quadtype)
    grid%nr = ngp
    basis%grid => grid

    NULLIFY(basis%am,basis%cm,basis%as,basis%ns,basis%bf,basis%dbf)

    ! fill in the basis data structures
    basis%eps_eig = 1.e-12_dp
    basis%basis_type = GTO_BASIS
    CALL Clementi_geobas(z,cval,aval,basis%nbas,starti,error)
    basis%nprim = basis%nbas
    m  = MAXVAL(basis%nbas)
    ALLOCATE (basis%am(m,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    basis%am = 0._dp
    DO l=0,3
      DO i=1,basis%nbas(l)
        ll = i - 1 + starti(l)
        basis%am(i,l) = aval * cval**(ll)
      END DO
    END DO

    basis%geometrical = .TRUE.
    basis%aval = aval
    basis%cval = cval
    basis%start = starti

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
      DO i=1,basis%nbas(l)
        al  = basis%am(i,l)
        DO k=1,nr
          rk  = basis%grid%rad(k)
          ear = EXP(-al*basis%grid%rad(k)**2)
          basis%bf(k,i,l) = rk**l * ear
          basis%dbf(k,i,l) = ( REAL(l,dp)*rk**(l-1) - 2._dp*al*rk**(l+1) ) * ear
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

    IF ( qs_kind%ghost ) THEN
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
    atom%state%maxn_calc = atom%state%maxn_occ

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
      CALL atom_fit_density (atom,ngto,0,iunit,results=results,error=error)
    ELSE
      CALL calculate_atom(atom,-1,error=error)
      CALL atom_fit_density (atom,ngto,0,-1,results=results,error=error)
    END IF

    xx = results(1)
    cc = results(2)
    DO i=1,ngto
       density(i,1) = xx * cc**i
       density(i,2) = results(2+i)
    END DO

    ! clean up
    CALL atom_int_release(integrals,error)
    CALL atom_ppint_release(integrals,error)
    CALL atom_relint_release(integrals,error)
    CALL release_atom_basis(basis,error)
    CALL release_atom_potential(potential,error)
    CALL release_atom_type(atom,error)

    DEALLOCATE (potential,basis,integrals,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE calculate_atomic_density

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param atomic_kind ...
!> \param qs_kind ...
!> \param rel_control ...
!> \param rtmat ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_atomic_relkin (atomic_kind,qs_kind,rel_control,rtmat,error)
    TYPE(atomic_kind_type), INTENT(IN)       :: atomic_kind
    TYPE(qs_kind_type), INTENT(IN)           :: qs_kind
    TYPE(rel_control_type), POINTER          :: rel_control
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rtmat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_atomic_relkin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, ii, ipgf, j, k, k1, &
                                                k2, l, ll, m, n, ngp, nj, nn, &
                                                nr, ns, nset, nsgf, quadtype, &
                                                relativistic, z
    INTEGER, DIMENSION(0:3, 10)              :: ncalc, ncore, nelem
    INTEGER, DIMENSION(0:3, 100)             :: set_index, shell_index
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgf, last_sgf, ls
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: al, alpha, ear, prefac, rk, &
                                                zeff
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: omat
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atom_basis_type), POINTER           :: basis
    TYPE(atom_integrals), POINTER            :: integrals
    TYPE(atom_potential_type), POINTER       :: potential
    TYPE(atom_type), POINTER                 :: atom
    TYPE(grid_atom_type), POINTER            :: grid
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

    failure = .FALSE.

    IF(rel_control%rel_method == rel_none) RETURN

    NULLIFY(all_potential,orb_basis_set)
    CALL get_qs_kind(qs_kind, basis_set=orb_basis_set, all_potential=all_potential, error=error)

    CPPostcondition(ASSOCIATED(orb_basis_set), cp_failure_level, routineP, error, failure)

    IF ( ASSOCIATED(all_potential) ) THEN
       ! only all electron atoms will get the relativistic correction

       CALL get_atomic_kind(atomic_kind,z=z)
       CALL get_qs_kind(qs_kind,zeff=zeff,error=error)
       NULLIFY(atom)
       CALL create_atom_type(atom,error)
       NULLIFY(atom%xc_section)
       NULLIFY(atom%orbitals)
       atom%z = z
       alpha = SQRT(all_potential%alpha_core_charge)

       ! set the method flag
       SELECT CASE (rel_control%rel_method)
          CASE DEFAULT
             CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
          CASE (rel_dkh)
             SELECT CASE (rel_control%rel_DKH_order)
                CASE DEFAULT
                   CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
                CASE (0)
                   relativistic = do_dkh0_atom
                CASE (1)
                   relativistic = do_dkh1_atom
                CASE (2)
                   relativistic = do_dkh2_atom
                CASE (3)
                   relativistic = do_dkh3_atom
                CASE (4)
                   relativistic = do_dkh4_atom
                CASE (5)
                   relativistic = do_dkh5_atom
             END SELECT
          CASE (rel_zora)
             SELECT CASE (rel_control%rel_zora_type)
                CASE DEFAULT
                   CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
                CASE (rel_zora_full)
                   CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
                CASE (rel_zora_mp)
                   relativistic = do_zoramp_atom
                CASE (rel_sczora_mp)
                   relativistic = do_sczoramp_atom
             END SELECT
       END SELECT

       CALL set_atom(atom,&
                  pp_calc=.FALSE.,&
                  method_type=do_rks_atom,&
                  relativistic=relativistic,&
                  coulomb_integral_type=do_numeric,&
                  exchange_integral_type=do_numeric,&
                  error=error)

       ALLOCATE (potential,basis,integrals,STAT=ierr)
       CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

       potential%ppot_type=no_pseudo
       CALL set_atom(atom,zcore=z,potential=potential,error=error)

       CALL get_gto_basis_set(orb_basis_set,&
          nset=nset,nshell=nshell,npgf=npgf,lmin=lmin,lmax=lmax,l=ls,nsgf=nsgf,zet=zet,gcc=gcc,&
          first_sgf=first_sgf,last_sgf=last_sgf)

       NULLIFY(grid)
       ngp = 400
       quadtype = do_gapw_log
       CALL allocate_grid_atom(grid,error)
       CALL create_grid_atom(grid,ngp,1,1,0,quadtype)
       grid%nr = ngp
       basis%grid => grid

       NULLIFY(basis%am,basis%cm,basis%as,basis%ns,basis%bf,basis%dbf)
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

       ! Normalization as used in the atomic code
       ! We have to undo the Quickstep normalization
       DO j=0,3
          prefac =  2.0_dp * SQRT( pi/dfac(2*j+1) )
          DO ipgf=1,basis%nprim(j)
             DO ii=1,basis%nbas(j)
                basis%cm(ipgf,ii,j) = prefac * basis%cm(ipgf,ii,j)
             END DO
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

       IF ( qs_kind%ghost ) THEN
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
       atom%state%maxn_calc = atom%state%maxn_occ

       ! calculate integrals
       ! general integrals
       CALL atom_int_setup(integrals,basis,error=error)
       ! potential
       CALL atom_ppint_setup(integrals,basis,potential=atom%potential,error=error)
       ! relativistic correction terms
       NULLIFY(integrals%tzora,integrals%hdkh)
       CALL atom_relint_setup(integrals,basis,atom%relativistic,zcore=REAL(atom%zcore,dp),&
                              alpha=alpha,error=error)
       CALL set_atom(atom,integrals=integrals,error=error)

       ! for DKH we need erfc integrals to correct non-relativistic
       integrals%core = 0.0_dp
       DO l=0,3
          n = integrals%n(l)
          m = basis%nprim(l)
          ALLOCATE (omat(m,m),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

          CALL sg_erfc ( omat(1:m,1:m), l, alpha, basis%am(1:m,l), basis%am(1:m,l) )
          integrals%core(1:n,1:n,l) = MATMUL(TRANSPOSE(basis%cm(1:m,1:n,l)),&
                                             MATMUL(omat(1:m,1:m),basis%cm(1:m,1:n,l)))

          DEALLOCATE (omat,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
       END DO

       ! recover relativistic kinetic matrix in CP2K/GPW order and normalization
       IF(ASSOCIATED(rtmat)) THEN
          DEALLOCATE (rtmat,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
       END IF
       ALLOCATE (rtmat(nsgf,nsgf),STAT=ierr)
       CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
       rtmat = 0._dp
       DO l=0,3
          ll = 2*l
          DO k1=1,basis%nbas(l)
             DO k2=1,basis%nbas(l)
                i=first_sgf(shell_index(l,k1),set_index(l,k1))
                j=first_sgf(shell_index(l,k2),set_index(l,k2))
                SELECT CASE (atom%relativistic)
                   CASE DEFAULT
                      CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
                   CASE (do_zoramp_atom,do_sczoramp_atom)
                      DO m=0,ll
                         rtmat(i+m,j+m) = integrals%tzora(k1,k2,l)
                      END DO
                   CASE (do_dkh0_atom,do_dkh1_atom,do_dkh2_atom,do_dkh3_atom,do_dkh4_atom,do_dkh5_atom)
                      DO m=0,ll
                         rtmat(i+m,j+m) = integrals%hdkh(k1,k2,l)-integrals%kin(k1,k2,l)+atom%zcore*integrals%core(k1,k2,l)
                      END DO
                END SELECT
             END DO
          END DO
       ENDDO
       DO k1=1,nsgf
          DO k2=k1,nsgf
             rtmat(k1,k2) = 0.5_dp*(rtmat(k1,k2) + rtmat(k2,k1))
             rtmat(k2,k1) = rtmat(k1,k2)
          END DO
       END DO

       ! clean up
       CALL atom_int_release(integrals,error)
       CALL atom_ppint_release(integrals,error)
       CALL atom_relint_release(integrals,error)
       CALL release_atom_basis(basis,error)
       CALL release_atom_potential(potential,error)
       CALL release_atom_type(atom,error)

       DEALLOCATE (potential,basis,integrals,STAT=ierr)
       CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    ELSE

       IF(ASSOCIATED(rtmat)) THEN
          DEALLOCATE (rtmat,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
       END IF
       NULLIFY(rtmat)

    END IF

  END SUBROUTINE calculate_atomic_relkin

! *****************************************************************************

END MODULE atom_kind_orbitals
