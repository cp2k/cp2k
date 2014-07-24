!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Define the atom type and its sub types
!> \author  jgh
!> \date    03.03.2008
!> \version 1.0
!>
! *****************************************************************************
MODULE atom_types
  USE cp_linked_list_val,              ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object,&
                                             parser_search_string,&
                                             parser_test_next_token
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE f77_blas
  USE input_constants,                 ONLY: contracted_gto,&
                                             gaussian,&
                                             geometrical_gto,&
                                             numerical,&
                                             slater
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_list_get,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: dfac,&
                                             fac,&
                                             pi
  USE periodic_table,                  ONLY: ptable
  USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                             create_grid_atom,&
                                             deallocate_grid_atom,&
                                             grid_atom_type
  USE string_utilities,                ONLY: remove_word,&
                                             uppercase
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_types'

  INTEGER, PARAMETER                                 :: GTO_BASIS=100,&
                                                        CGTO_BASIS=101,&
                                                        STO_BASIS=102,&
                                                        NUM_BASIS=103
  INTEGER, PARAMETER                                 :: NO_PSEUDO=0,&
                                                        GTH_PSEUDO=1

!> Provides all information about a basis set
! *****************************************************************************
  TYPE atom_basis_type
    INTEGER                                       :: basis_type
    INTEGER, DIMENSION(0:3)                       :: nbas
    INTEGER, DIMENSION(0:3)                       :: nprim
    REAL(KIND=dp),DIMENSION(:,:),POINTER          :: am         !GTO exponents
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: cm         !Contraction coeffs
    REAL(KIND=dp),DIMENSION(:,:),POINTER          :: as         !STO exponents
    INTEGER,DIMENSION(:,:),POINTER                :: ns         !STO n-quantum numbers
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: bf         !num. bsf
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: dbf        !derivatives (num)
    REAL(KIND=dp)                                 :: eps_eig
    TYPE(grid_atom_type),POINTER                  :: grid
    LOGICAL                                       :: geometrical
    REAL(KIND=dp)                                 :: aval, cval
    INTEGER, DIMENSION(0:3)                       :: start
  END TYPE atom_basis_type

!> Provides all information about a pseudopotential
! *****************************************************************************
  TYPE atom_gthpot_type
    CHARACTER (LEN=2)                             :: symbol
    CHARACTER (LEN=default_string_length)         :: pname
    INTEGER, DIMENSION(0:3)                       :: econf
    REAL(dp)                                      :: zion
    REAL(dp)                                      :: rc
    INTEGER                                       :: ncl
    REAL(dp), DIMENSION(5)                        :: cl
    INTEGER, DIMENSION(0:3)                       :: nl
    REAL(dp), DIMENSION(0:3)                      :: rcnl
    REAL(dp), DIMENSION(4,4,0:3)                  :: hnl
    ! type extensions
    ! NLCC
    LOGICAL                                    :: nlcc
    INTEGER                                    :: nexp_nlcc
    REAL(KIND = dp), DIMENSION(10)             :: alpha_nlcc
    INTEGER, DIMENSION(10)                     :: nct_nlcc
    REAL(KIND = dp), DIMENSION(4,10)           :: cval_nlcc
    ! LSD potential
    LOGICAL                                    :: lsdpot
    INTEGER                                    :: nexp_lsd
    REAL(KIND = dp), DIMENSION(10)             :: alpha_lsd
    INTEGER, DIMENSION(10)                     :: nct_lsd
    REAL(KIND = dp), DIMENSION(4,10)           :: cval_lsd
    ! extended local potential
    LOGICAL                                    :: lpotextended
    INTEGER                                    :: nexp_lpot
    REAL(KIND = dp), DIMENSION(10)             :: alpha_lpot
    INTEGER, DIMENSION(10)                     :: nct_lpot
    REAL(KIND = dp), DIMENSION(4,10)           :: cval_lpot
  END TYPE atom_gthpot_type

  TYPE atom_potential_type
    INTEGER                                       :: ppot_type
    LOGICAL                                       :: confinement
    REAL(dp)                                      :: acon
    REAL(dp)                                      :: rcon
    REAL(dp)                                      :: scon
    TYPE(atom_gthpot_type)                        :: gth_pot
  END TYPE atom_potential_type

!> Provides all information on states and occupation
! *****************************************************************************
  TYPE atom_state
    REAL(KIND=dp),DIMENSION(0:3,10)               :: occ
    REAL(KIND=dp),DIMENSION(0:3,10)               :: core
    REAL(KIND=dp),DIMENSION(0:3,10)               :: occupation
    INTEGER                                       :: maxl_occ
    INTEGER,DIMENSION(0:3)                        :: maxn_occ
    INTEGER                                       :: maxl_calc
    INTEGER,DIMENSION(0:3)                        :: maxn_calc
    INTEGER                                       :: multiplicity
    REAL(KIND=dp),DIMENSION(0:3,10)               :: occa, occb
  END TYPE atom_state

!> Holds atomic integrals
! *****************************************************************************
  TYPE eri
    REAL(KIND=dp),DIMENSION(:,:),POINTER          :: int
  END TYPE eri

  TYPE atom_integrals
    INTEGER                                       :: status=0
    INTEGER                                       :: ppstat=0
    LOGICAL                                       :: eri_coulomb
    LOGICAL                                       :: eri_exchange
    LOGICAL                                       :: all_nu
    INTEGER, DIMENSION(0:3)                       :: n, nne
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: ovlp, kin, core, clsd
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: utrans, uptrans
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: hnl
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: conf
    TYPE(eri),DIMENSION(20)                       :: ceri
    TYPE(eri),DIMENSION(20)                       :: eeri
    INTEGER                                       :: dkhstat=0
    INTEGER                                       :: zorastat=0
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: tzora
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: hdkh
  END TYPE atom_integrals

!> Holds atomic orbitals and energies
! *****************************************************************************
  TYPE atom_orbitals
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: wfn,wfna,wfnb
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: pmat,pmata,pmatb
    REAL(KIND=dp),DIMENSION(:,:),POINTER          :: ener,enera,enerb
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: refene,refchg,refnod
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: wrefene,wrefchg,wrefnod
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: crefene,crefchg,crefnod
    REAL(KIND=dp),DIMENSION(:,:),POINTER          :: wpsir0
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: rcmax
  END TYPE atom_orbitals

!> Operator matrices
! *****************************************************************************
  TYPE opmat_type
    INTEGER, DIMENSION(0:3)                       :: n
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: op
  END TYPE opmat_type

!> Operator grids
! *****************************************************************************
  TYPE opgrid_type
    REAL(KIND=dp),DIMENSION(:),POINTER            :: op
    TYPE(grid_atom_type),POINTER                  :: grid
  END TYPE opgrid_type

!> All energies
! *****************************************************************************
  TYPE atom_energy_type
    REAL(KIND=dp)                                 :: etot
    REAL(KIND=dp)                                 :: eband
    REAL(KIND=dp)                                 :: ekin
    REAL(KIND=dp)                                 :: epot
    REAL(KIND=dp)                                 :: ecore
    REAL(KIND=dp)                                 :: elsd
    REAL(KIND=dp)                                 :: epseudo
    REAL(KIND=dp)                                 :: eploc
    REAL(KIND=dp)                                 :: epnl
    REAL(KIND=dp)                                 :: exc
    REAL(KIND=dp)                                 :: ecoulomb
    REAL(KIND=dp)                                 :: eexchange
    REAL(KIND=dp)                                 :: econfinement
  END TYPE atom_energy_type

!> Information on optimization procedure
! *****************************************************************************
  TYPE atom_optimization_type
    REAL(KIND=dp)                                 :: damping
    REAL(KIND=dp)                                 :: eps_scf
    REAL(KIND=dp)                                 :: eps_diis
    INTEGER                                       :: max_iter
    INTEGER                                       :: n_diis
  END TYPE atom_optimization_type

!> Provides all information about an atomic kind
! *****************************************************************************
  TYPE atom_type
    INTEGER                                       :: z
    INTEGER                                       :: zcore
    LOGICAL                                       :: pp_calc
! ZMP adding in type some variables
    LOGICAL                                       :: do_zmp, doread, read_vxc, dm
    CHARACTER(LEN=default_string_length)          :: ext_file, ext_vxc_file, &
                                                     zmp_restart_file
    INTEGER                                       :: method_type
    INTEGER                                       :: relativistic
    INTEGER                                       :: coulomb_integral_type
    INTEGER                                       :: exchange_integral_type
! ZMP
    REAL(KIND=dp)                                 :: lambda
    REAL(KIND=dp)                                 :: rho_diff_integral
    REAL(KIND=dp)                                 :: weight, zmpgrid_tol, zmpvxcgrid_tol
    TYPE(atom_basis_type), POINTER                :: basis
    TYPE(atom_potential_type), POINTER            :: potential
    TYPE(atom_state), POINTER                     :: state
    TYPE(atom_integrals), POINTER                 :: integrals
    TYPE(atom_orbitals), POINTER                  :: orbitals
    TYPE(atom_energy_type)                        :: energy
    TYPE(atom_optimization_type)                  :: optimization
    TYPE(section_vals_type), POINTER              :: xc_section, zmp_section
  END TYPE atom_type
! *****************************************************************************
  TYPE atom_p_type
    TYPE(atom_type), POINTER                      :: atom
  END TYPE atom_p_type

  PUBLIC :: atom_p_type, atom_type, atom_basis_type, atom_state, atom_integrals
  PUBLIC :: atom_orbitals, eri, atom_potential_type, atom_gthpot_type
  PUBLIC :: atom_optimization_type
  PUBLIC :: create_atom_type, release_atom_type, set_atom
  PUBLIC :: create_atom_orbs
  PUBLIC :: init_atom_basis, release_atom_basis
  PUBLIC :: init_atom_potential, release_atom_potential
  PUBLIC :: read_atom_opt_section
  PUBLIC :: Clementi_geobas
  PUBLIC :: GTO_BASIS, CGTO_BASIS, STO_BASIS, NUM_BASIS
  PUBLIC :: opmat_type, create_opmat, release_opmat
  PUBLIC :: opgrid_type, create_opgrid, release_opgrid
  PUBLIC :: NO_PSEUDO, GTH_PSEUDO

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief Initialize the basis for the atomic code
!> \param basis ...
!> \param basis_section ...
!> \param zval ...
!> \param btyp ...
!> \param error ...
!> \note  Highly accurate relativistic universal Gaussian basis set: Dirac-Fock-Coulomb calculations
!>        for atomic systems up to nobelium
!>        J. Chem. Phys. 101, 6829 (1994); DOI:10.1063/1.468311
!>        G. L. Malli and A. B. F. Da Silva
!>        Department of Chemistry, Simon Fraser University, Burnaby, B.C., Canada
!>        Yasuyuki Ishikawa
!>        Department of Chemistry, University of Puerto Rico, San Juan, Puerto Rico
!>
!>        A universal Gaussian basis set is developed that leads to relativistic Dirac-Fock SCF energies
!>        of comparable accuracy as that obtained by the accurate numerical finite-difference method
!>        (GRASP2 package) [J. Phys. B 25, 1 (1992)]. The Gaussian-type functions of our universal basis
!>        set satisfy the relativistic boundary conditions associated with the finite nuclear model for a
!>        finite speed of light and conform to the so-called kinetic balance at the nonrelativistic limit.
!>        We attribute the exceptionally high accuracy obtained in our calculations to the fact that the
!>        representation of the relativistic dynamics of an electron in a spherical ball finite nucleus
!>        near the origin in terms of our universal Gaussian basis set is as accurate as that provided by
!>        the numerical finite-difference method. Results of the Dirac-Fock-Coulomb energies for a number
!>        of atoms up to No (Z=102) and some negative ions are presented and compared with the recent
!>        results obtained with the numerical finite-difference method and geometrical Gaussian basis sets
!>        by Parpia, Mohanty, and Clementi [J. Phys. B 25, 1 (1992)]. The accuracy of our calculations is
!>        estimated to be within a few parts in 109 for all the atomic systems studied.
! *****************************************************************************
  SUBROUTINE init_atom_basis(basis,basis_section,zval,btyp,error)
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    TYPE(section_vals_type), POINTER         :: basis_section
    INTEGER, INTENT(IN)                      :: zval
    CHARACTER(LEN=2)                         :: btyp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_atom_basis', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nua = 40, nup = 16
    REAL(KIND=dp), DIMENSION(nua), PARAMETER :: ugbs = (/0.007299_dp,&
      0.013705_dp, 0.025733_dp, 0.048316_dp, 0.090718_dp, 0.170333_dp,&
      0.319819_dp, 0.600496_dp, 1.127497_dp, 2.117000_dp, 3.974902_dp,&
      7.463317_dp, 14.013204_dp, 26.311339_dp, 49.402449_dp, 92.758561_dp,&
      174.164456_dp, 327.013024_dp, 614.003114_dp, 1152.858743_dp,&
      2164.619772_dp, 4064.312984_dp, 7631.197056_dp, 14328.416324_dp,&
      26903.186074_dp, 50513.706789_dp, 94845.070265_dp, 178082.107320_dp,&
      334368.848683_dp, 627814.487663_dp, 1178791.123851_dp, 2213310.684886_dp&
      ,4155735.557141_dp,7802853.046713_dp, 14650719.428954_dp,&
      27508345.793637_dp, 51649961.080194_dp, 96978513.342764_dp,&
      182087882.613702_dp, 341890134.751331_dp /)

    CHARACTER(LEN=default_string_length)     :: basis_fn, basis_name
    INTEGER                                  :: basistype, i, ierr, j, k, l, &
                                                ll, m, ngp, nl, nr, nu, &
                                                quadtype
    INTEGER, DIMENSION(0:3)                  :: starti
    INTEGER, DIMENSION(:), POINTER           :: nqm, num_gto, num_slater, &
                                                sindex
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: al, amax, aval, cval, ear, &
                                                pf, rk
    REAL(KIND=dp), DIMENSION(:), POINTER     :: expo
    TYPE(section_vals_type), POINTER         :: gto_basis_section

!   btyp = AE : standard all-electron basis
!   btyp = PP : standard pseudopotential basis
!   btyp = AA : high accuracy all-electron basis
!   btyp = AP : high accuracy pseudopotential basis

    failure = .FALSE.
    NULLIFY(basis%am,basis%cm,basis%as,basis%ns,basis%bf,basis%dbf)
    ! get information on quadrature type and number of grid points
    ! allocate and initialize the atomic grid
    CALL allocate_grid_atom(basis%grid,error)
    CALL section_vals_val_get(basis_section,"QUADRATURE",i_val=quadtype,error=error)
    CALL section_vals_val_get(basis_section,"GRID_POINTS",i_val=ngp,error=error)
    CALL cp_assert(ngp > 0,cp_failure_level,cp_assertion_failed,routineP,&
                   "# point radial grid < 0",error,failure)
    CALL create_grid_atom(basis%grid,ngp,1,1,0,quadtype)
    basis%grid%nr = ngp
    basis%geometrical = .FALSE.
    basis%aval  = 0._dp
    basis%cval  = 0._dp
    basis%start = 0

    CALL section_vals_val_get(basis_section,"BASIS_TYPE",i_val=basistype,error=error)
    CALL section_vals_val_get(basis_section,"EPS_EIGENVALUE",r_val=basis%eps_eig,error=error)
    SELECT CASE (basistype)
      CASE DEFAULT
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      CASE (gaussian)
         basis%basis_type = GTO_BASIS
         NULLIFY(num_gto)
         CALL section_vals_val_get(basis_section,"NUM_GTO",i_vals=num_gto,error=error)
         IF ( num_gto(1) < 1 ) THEN
           ! use default basis
           IF ( btyp == "AE" ) THEN
             nu = nua
           ELSEIF ( btyp == "PP" ) THEN
             nu = nup
           ELSE
             nu = nua
           ENDIF
           basis%nbas = nu
           basis%nprim = nu
           ALLOCATE (basis%am(nu,0:3),STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
           basis%am(1:nu,0) = ugbs(1:nu)
           basis%am(1:nu,1) = ugbs(1:nu)
           basis%am(1:nu,2) = ugbs(1:nu)
           basis%am(1:nu,3) = ugbs(1:nu)
         ELSE
           basis%nbas = 0
           DO i=1,SIZE(num_gto)
             basis%nbas(i-1) = num_gto(i)
           END DO
           basis%nprim = basis%nbas
           m  = MAXVAL(basis%nbas)
           ALLOCATE (basis%am(m,0:3),STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
           basis%am = 0._dp
           IF ( basis%nbas(0) > 0 ) THEN
             NULLIFY(expo)
             CALL section_vals_val_get(basis_section,"S_EXPONENTS",r_vals=expo,error=error)
             CPPostcondition(SIZE(expo)>=basis%nbas(0), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(0)
               basis%am(i,0) = expo(i)
             END DO
           END IF
           IF ( basis%nbas(1) > 0 ) THEN
             NULLIFY(expo)
             CALL section_vals_val_get(basis_section,"P_EXPONENTS",r_vals=expo,error=error)
             CPPostcondition(SIZE(expo)>=basis%nbas(1), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(1)
               basis%am(i,1) = expo(i)
             END DO
           END IF
           IF ( basis%nbas(2) > 0 ) THEN
             NULLIFY(expo)
             CALL section_vals_val_get(basis_section,"D_EXPONENTS",r_vals=expo,error=error)
             CPPostcondition(SIZE(expo)>=basis%nbas(2), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(2)
               basis%am(i,2) = expo(i)
             END DO
           END IF
           IF ( basis%nbas(3) > 0 ) THEN
             NULLIFY(expo)
             CALL section_vals_val_get(basis_section,"F_EXPONENTS",r_vals=expo,error=error)
             CPPostcondition(SIZE(expo)>=basis%nbas(3), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(3)
               basis%am(i,3) = expo(i)
             END DO
           END IF
         END IF
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
      CASE (geometrical_gto)
         basis%basis_type = GTO_BASIS
         NULLIFY(num_gto)
         CALL section_vals_val_get(basis_section,"NUM_GTO",i_vals=num_gto,error=error)
         IF ( num_gto(1) < 1 ) THEN
           IF ( btyp == "AE" ) THEN
             ! use the Clementi extra large basis
             CALL Clementi_geobas(zval,cval,aval,basis%nbas,starti,error)
           ELSEIF ( btyp == "PP" ) THEN
             ! use the Clementi extra large basis
             CALL Clementi_geobas(zval,cval,aval,basis%nbas,starti,error)
           ELSEIF ( btyp == "AA" ) THEN
             CALL Clementi_geobas(zval,cval,aval,basis%nbas,starti,error)
             amax = cval**(basis%nbas(0)-1)
             basis%nbas(0) = NINT((LOG(amax)/LOG(1.6_dp)))
             cval = 1.6_dp
             starti = 0
             basis%nbas(1) = basis%nbas(0) - 4
             basis%nbas(2) = basis%nbas(0) - 8
             basis%nbas(3) = basis%nbas(0) -12
           ELSEIF ( btyp == "AP" ) THEN
             CALL Clementi_geobas(zval,cval,aval,basis%nbas,starti,error)
             amax = 500._dp/aval
             basis%nbas = NINT((LOG(amax)/LOG(1.6_dp)))
             cval = 1.6_dp
             starti = 0
           ELSE
             ! use the Clementi extra large basis
             CALL Clementi_geobas(zval,cval,aval,basis%nbas,starti,error)
           ENDIF
           basis%nprim = basis%nbas
         ELSE
           basis%nbas = 0
           DO i=1,SIZE(num_gto)
             basis%nbas(i-1) = num_gto(i)
           END DO
           basis%nprim = basis%nbas
           NULLIFY(sindex)
           CALL section_vals_val_get(basis_section,"START_INDEX",i_vals=sindex,error=error)
           starti = 0
           DO i=1,SIZE(sindex)
             starti(i-1) = sindex(i)
             CPPostcondition(sindex(i)>=0, cp_failure_level, routineP, error, failure)
           END DO
           CALL section_vals_val_get(basis_section,"GEOMETRICAL_FACTOR",r_val=cval,error=error)
           CALL section_vals_val_get(basis_section,"GEO_START_VALUE",r_val=aval,error=error)
         END IF
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
      CASE (contracted_gto)
         basis%basis_type = CGTO_BASIS
         CALL section_vals_val_get(basis_section,"BASIS_SET_FILE_NAME",c_val=basis_fn,error=error)
         CALL section_vals_val_get(basis_section,"BASIS_SET",c_val=basis_name,error=error)
         gto_basis_section => section_vals_get_subs_vals(basis_section,"BASIS",error=error)
         CALL read_basis_set(ptable(zval)%symbol,basis,basis_name,basis_fn,&
                             gto_basis_section,error)

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
      CASE (slater)
         basis%basis_type = STO_BASIS
         NULLIFY(num_slater)
         CALL section_vals_val_get(basis_section,"NUM_SLATER",i_vals=num_slater,error=error)
         IF ( num_slater(1) < 1 ) THEN
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
         ELSE
           basis%nbas = 0
           DO i=1,SIZE(num_slater)
             basis%nbas(i-1) = num_slater(i)
           END DO
           basis%nprim = basis%nbas
           m  = MAXVAL(basis%nbas)
           ALLOCATE (basis%as(m,0:3),basis%ns(m,0:3),STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
           basis%as = 0._dp
           basis%ns = 0
           IF ( basis%nbas(0) > 0 ) THEN
             NULLIFY(expo)
             CALL section_vals_val_get(basis_section,"S_EXPONENTS",r_vals=expo,error=error)
             CPPostcondition(SIZE(expo)>=basis%nbas(0), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(0)
               basis%as(i,0) = expo(i)
             END DO
             NULLIFY(nqm)
             CALL section_vals_val_get(basis_section,"S_QUANTUM_NUMBERS",i_vals=nqm,error=error)
             CPPostcondition(SIZE(nqm)>=basis%nbas(0), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(0)
               basis%ns(i,0) = nqm(i)
             END DO
           END IF
           IF ( basis%nbas(1) > 0 ) THEN
             NULLIFY(expo)
             CALL section_vals_val_get(basis_section,"P_EXPONENTS",r_vals=expo,error=error)
             CPPostcondition(SIZE(expo)>=basis%nbas(1), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(1)
               basis%as(i,1) = expo(i)
             END DO
             NULLIFY(nqm)
             CALL section_vals_val_get(basis_section,"P_QUANTUM_NUMBERS",i_vals=nqm,error=error)
             CPPostcondition(SIZE(nqm)>=basis%nbas(1), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(1)
               basis%ns(i,1) = nqm(i)
             END DO
           END IF
           IF ( basis%nbas(2) > 0 ) THEN
             NULLIFY(expo)
             CALL section_vals_val_get(basis_section,"D_EXPONENTS",r_vals=expo,error=error)
             CPPostcondition(SIZE(expo)>=basis%nbas(2), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(2)
               basis%as(i,2) = expo(i)
             END DO
             NULLIFY(nqm)
             CALL section_vals_val_get(basis_section,"D_QUANTUM_NUMBERS",i_vals=nqm,error=error)
             CPPostcondition(SIZE(nqm)>=basis%nbas(2), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(2)
               basis%ns(i,2) = nqm(i)
             END DO
           END IF
           IF ( basis%nbas(3) > 0 ) THEN
             NULLIFY(expo)
             CALL section_vals_val_get(basis_section,"F_EXPONENTS",r_vals=expo,error=error)
             CPPostcondition(SIZE(expo)>=basis%nbas(3), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(3)
               basis%as(i,3) = expo(i)
             END DO
             NULLIFY(nqm)
             CALL section_vals_val_get(basis_section,"F_QUANTUM_NUMBERS",i_vals=nqm,error=error)
             CPPostcondition(SIZE(nqm)>=basis%nbas(3), cp_failure_level, routineP, error, failure)
             DO i=1,basis%nbas(3)
               basis%ns(i,3) = nqm(i)
             END DO
           END IF
         END IF
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
             al  = basis%as(i,l)
             nl  = basis%ns(i,l)
             pf  = (2._dp*al)**nl * SQRT(2._dp*al/fac(2*nl))
             DO k=1,nr
               rk  = basis%grid%rad(k)
               ear = EXP(-al*rk)
               basis%bf(k,i,l) = pf * rk**(nl-1) * ear
               basis%dbf(k,i,l) = pf * ( REAL(nl-1,dp)/rk - al ) * rk**(nl-1) * ear
             END DO
           END DO
         END DO
      CASE (numerical)
         basis%basis_type = NUM_BASIS
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE init_atom_basis

! *****************************************************************************
!> \brief ...
!> \param basis ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE release_atom_basis(basis,error)
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_atom_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF(ASSOCIATED(basis%am)) THEN
       DEALLOCATE (basis%am,STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF
    IF(ASSOCIATED(basis%cm)) THEN
       DEALLOCATE (basis%cm,STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF
    IF(ASSOCIATED(basis%as)) THEN
       DEALLOCATE (basis%as,STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF
    IF(ASSOCIATED(basis%ns)) THEN
       DEALLOCATE (basis%ns,STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF
    IF(ASSOCIATED(basis%bf)) THEN
       DEALLOCATE (basis%bf,STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF
    IF(ASSOCIATED(basis%dbf)) THEN
       DEALLOCATE (basis%dbf,STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF

    CALL deallocate_grid_atom(basis%grid,error)

  END SUBROUTINE release_atom_basis
! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param atom ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_atom_type(atom,error)
    TYPE(atom_type), POINTER                 :: atom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_atom_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPAssert(.NOT.ASSOCIATED(atom),cp_failure_level,routineP,error,failure)

    ALLOCATE(atom,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    NULLIFY(atom%zmp_section)
    NULLIFY(atom%xc_section)
    atom%do_zmp = .FALSE.  
    atom%doread = .FALSE. 
    atom%read_vxc = .FALSE.
    atom%dm = .FALSE. 

  END SUBROUTINE create_atom_type

! *****************************************************************************
!> \brief ...
!> \param atom ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE release_atom_type(atom,error)
    TYPE(atom_type), POINTER                 :: atom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_atom_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPAssert(ASSOCIATED(atom),cp_failure_level,routineP,error,failure)

    NULLIFY(atom%basis)
    NULLIFY(atom%integrals)
    IF(ASSOCIATED(atom%state)) THEN
      DEALLOCATE(atom%state,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(atom%state)
    END IF
    IF(ASSOCIATED(atom%orbitals)) THEN
      CALL release_atom_orbs(atom%orbitals,error)
    END IF

    DEALLOCATE(atom,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    NULLIFY(atom)

  END SUBROUTINE release_atom_type

! ZMP adding input variables in subroutine do_zmp,doread,read_vxc,method_type
! *****************************************************************************
!> \brief ...
!> \param atom ...
!> \param basis ...
!> \param state ...
!> \param integrals ...
!> \param orbitals ...
!> \param potential ...
!> \param zcore ...
!> \param pp_calc ...
!> \param do_zmp ...
!> \param doread ...
!> \param read_vxc ...
!> \param method_type ...
!> \param relativistic ...
!> \param coulomb_integral_type ...
!> \param exchange_integral_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE set_atom(atom,basis,state,integrals,orbitals,potential,zcore,pp_calc,do_zmp,doread,&
    read_vxc,method_type,relativistic,coulomb_integral_type,exchange_integral_type,error)
    TYPE(atom_type), POINTER                 :: atom
    TYPE(atom_basis_type), OPTIONAL, POINTER :: basis
    TYPE(atom_state), OPTIONAL, POINTER      :: state
    TYPE(atom_integrals), OPTIONAL, POINTER  :: integrals
    TYPE(atom_orbitals), OPTIONAL, POINTER   :: orbitals
    TYPE(atom_potential_type), OPTIONAL, &
      POINTER                                :: potential
    INTEGER, INTENT(IN), OPTIONAL            :: zcore
    LOGICAL, INTENT(IN), OPTIONAL            :: pp_calc, do_zmp, doread, &
                                                read_vxc
    INTEGER, INTENT(IN), OPTIONAL            :: method_type, relativistic, &
                                                coulomb_integral_type, &
                                                exchange_integral_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_atom', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    CPAssert(ASSOCIATED(atom),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
      IF(PRESENT(basis)) atom%basis => basis
      IF(PRESENT(state)) atom%state => state
      IF(PRESENT(integrals)) atom%integrals => integrals
      IF(PRESENT(orbitals)) atom%orbitals => orbitals
      IF(PRESENT(potential)) atom%potential => potential
      IF(PRESENT(zcore)) atom%zcore = zcore
      IF(PRESENT(pp_calc)) atom%pp_calc = pp_calc
! ZMP assigning variable values if present
      IF(PRESENT(do_zmp)) atom%do_zmp = do_zmp
      IF(PRESENT(doread)) atom%doread = doread
      IF(PRESENT(read_vxc)) atom%read_vxc = read_vxc
 
      IF(PRESENT(method_type)) atom%method_type = method_type
      IF(PRESENT(relativistic)) atom%relativistic = relativistic
      IF(PRESENT(coulomb_integral_type)) atom%coulomb_integral_type = coulomb_integral_type
      IF(PRESENT(exchange_integral_type)) atom%exchange_integral_type = exchange_integral_type
    END IF

  END SUBROUTINE set_atom

! *****************************************************************************
!> \brief ...
!> \param orbs ...
!> \param mbas ...
!> \param mo ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_atom_orbs(orbs,mbas,mo,error)
    TYPE(atom_orbitals), POINTER             :: orbs
    INTEGER, INTENT(IN)                      :: mbas, mo
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_atom_orbs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPAssert(.NOT.ASSOCIATED(orbs),cp_failure_level,routineP,error,failure)

    ALLOCATE(orbs,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    ALLOCATE(orbs%wfn(mbas,mo,0:3),orbs%wfna(mbas,mo,0:3),orbs%wfnb(mbas,mo,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    orbs%wfn = 0._dp
    orbs%wfna = 0._dp
    orbs%wfnb = 0._dp

    ALLOCATE(orbs%pmat(mbas,mbas,0:3),orbs%pmata(mbas,mbas,0:3),orbs%pmatb(mbas,mbas,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    orbs%pmat = 0._dp
    orbs%pmata = 0._dp
    orbs%pmatb = 0._dp

    ALLOCATE(orbs%ener(mo,0:3),orbs%enera(mo,0:3),orbs%enerb(mo,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    orbs%ener = 0._dp
    orbs%enera = 0._dp
    orbs%enerb = 0._dp

    ALLOCATE(orbs%refene(mo,0:3,2),orbs%refchg(mo,0:3,2),orbs%refnod(mo,0:3,2),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    orbs%refene = 0._dp
    orbs%refchg = 0._dp
    orbs%refnod = 0._dp
    ALLOCATE(orbs%wrefene(mo,0:3,2),orbs%wrefchg(mo,0:3,2),orbs%wrefnod(mo,0:3,2),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    orbs%wrefene = 0._dp
    orbs%wrefchg = 0._dp
    orbs%wrefnod = 0._dp
    ALLOCATE(orbs%crefene(mo,0:3,2),orbs%crefchg(mo,0:3,2),orbs%crefnod(mo,0:3,2),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    orbs%crefene = 0._dp
    orbs%crefchg = 0._dp
    orbs%crefnod = 0._dp
    ALLOCATE(orbs%rcmax(mo,0:3,2),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    orbs%rcmax = 0._dp
    ALLOCATE(orbs%wpsir0(mo,2),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    orbs%wpsir0 = 0._dp

  END SUBROUTINE create_atom_orbs

! *****************************************************************************
!> \brief ...
!> \param orbs ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE release_atom_orbs(orbs,error)
    TYPE(atom_orbitals), POINTER             :: orbs
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_atom_orbs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPAssert(ASSOCIATED(orbs),cp_failure_level,routineP,error,failure)

    IF(ASSOCIATED(orbs%wfn)) THEN
      DEALLOCATE(orbs%wfn,orbs%wfna,orbs%wfnb,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%wfn,orbs%wfna,orbs%wfnb)
    END IF
    IF(ASSOCIATED(orbs%pmat)) THEN
      DEALLOCATE(orbs%pmat,orbs%pmata,orbs%pmatb,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%pmat,orbs%pmata,orbs%pmatb)
    END IF
    IF(ASSOCIATED(orbs%ener)) THEN
      DEALLOCATE(orbs%ener,orbs%enera,orbs%enerb,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%ener,orbs%enera,orbs%enerb)
    END IF
    IF(ASSOCIATED(orbs%refene)) THEN
      DEALLOCATE(orbs%refene,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%refene)
    END IF
    IF(ASSOCIATED(orbs%refchg)) THEN
      DEALLOCATE(orbs%refchg,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%refchg)
    END IF
    IF(ASSOCIATED(orbs%refnod)) THEN
      DEALLOCATE(orbs%refnod,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%refnod)
    END IF
    IF(ASSOCIATED(orbs%wrefene)) THEN
      DEALLOCATE(orbs%wrefene,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%wrefene)
    END IF
    IF(ASSOCIATED(orbs%wrefchg)) THEN
      DEALLOCATE(orbs%wrefchg,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%wrefchg)
    END IF
    IF(ASSOCIATED(orbs%wrefnod)) THEN
      DEALLOCATE(orbs%wrefnod,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%wrefnod)
    END IF
    IF(ASSOCIATED(orbs%crefene)) THEN
      DEALLOCATE(orbs%crefene,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%crefene)
    END IF
    IF(ASSOCIATED(orbs%crefchg)) THEN
      DEALLOCATE(orbs%crefchg,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%crefchg)
    END IF
    IF(ASSOCIATED(orbs%crefnod)) THEN
      DEALLOCATE(orbs%crefnod,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%crefnod)
    END IF
    IF(ASSOCIATED(orbs%rcmax)) THEN
      DEALLOCATE(orbs%rcmax,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%rcmax)
    END IF
    IF(ASSOCIATED(orbs%wpsir0)) THEN
      DEALLOCATE(orbs%wpsir0,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%wpsir0)
    END IF

    DEALLOCATE(orbs,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    NULLIFY(orbs)

  END SUBROUTINE release_atom_orbs

! *****************************************************************************
!> \brief ...
!> \param opmat ...
!> \param n ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_opmat(opmat,n,error)
    TYPE(opmat_type), POINTER                :: opmat
    INTEGER, DIMENSION(0:3), INTENT(IN)      :: n
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_opmat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, m
    LOGICAL                                  :: failure

    failure = .FALSE.
    m=MAXVAL(n)

    CPPrecondition(.NOT.ASSOCIATED(opmat), cp_failure_level, routineP, error, failure)

    ALLOCATE(opmat,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    opmat%n = n
    ALLOCATE(opmat%op(m,m,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    opmat%op = 0._dp

  END SUBROUTINE create_opmat

! *****************************************************************************
!> \brief ...
!> \param opmat ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE release_opmat(opmat,error)
    TYPE(opmat_type), POINTER                :: opmat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_opmat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(opmat), cp_failure_level, routineP, error, failure)

    opmat%n = 0
    DEALLOCATE(opmat%op,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DEALLOCATE(opmat,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE release_opmat

! *****************************************************************************
!> \brief ...
!> \param opgrid ...
!> \param grid ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_opgrid(opgrid,grid,error)
    TYPE(opgrid_type), POINTER               :: opgrid
    TYPE(grid_atom_type), POINTER            :: grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_opgrid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, nr
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(opgrid), cp_failure_level, routineP, error, failure)

    ALLOCATE(opgrid,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    opgrid%grid => grid

    nr = grid%nr

    ALLOCATE(opgrid%op(nr),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    opgrid%op = 0._dp

  END SUBROUTINE create_opgrid

! *****************************************************************************
!> \brief ...
!> \param opgrid ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE release_opgrid(opgrid,error)
    TYPE(opgrid_type), POINTER               :: opgrid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_opgrid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(opgrid), cp_failure_level, routineP, error, failure)

    NULLIFY(opgrid%grid)
    DEALLOCATE(opgrid%op,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DEALLOCATE(opgrid,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE release_opgrid

! *****************************************************************************
!> \brief ...
!> \param zval ...
!> \param cval ...
!> \param aval ...
!> \param ngto ...
!> \param ival ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE Clementi_geobas(zval,cval,aval,ngto,ival,error)
    INTEGER, INTENT(IN)                      :: zval
    REAL(dp), INTENT(OUT)                    :: cval, aval
    INTEGER, DIMENSION(0:3), INTENT(OUT)     :: ngto, ival
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'Clementi_geobas', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    ngto = 0
    ival = 0
    cval = 0._dp
    aval = 0._dp

    SELECT CASE (zval)
      CASE DEFAULT
        CPPrecondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (1)   ! this is from the general geometrical basis and extended
        cval = 2.0_dp
        aval = 0.016_dp
        ngto(0) = 20
      CASE (2)
        cval = 2.14774520_dp
        aval = 0.04850670_dp
        ngto(0) = 20
      CASE (3)
        cval = 2.08932430_dp
        aval = 0.02031060_dp
        ngto(0) = 23
      CASE (4)
        cval = 2.09753060_dp
        aval = 0.03207070_dp
        ngto(0) = 23
      CASE (5)
        cval = 2.10343410_dp
        aval = 0.03591970_dp
        ngto(0) = 23
        ngto(1) = 16
      CASE (6)
        cval = 2.10662820_dp
        aval = 0.05292410_dp
        ngto(0) = 23
        ngto(1) = 16
      CASE (7)
        cval = 2.13743840_dp
        aval = 0.06291970_dp
        ngto(0) = 23
        ngto(1) = 16
      CASE (8)
        cval = 2.08687310_dp
        aval = 0.08350860_dp
        ngto(0) = 23
        ngto(1) = 16
      CASE (9)
        cval = 2.12318180_dp
        aval = 0.09899170_dp
        ngto(0) = 23
        ngto(1) = 16
      CASE (10)
        cval = 2.13164810_dp
        aval = 0.11485350_dp
        ngto(0) = 23
        ngto(1) = 16
      CASE (11)
        cval = 2.11413310_dp
        aval = 0.00922630_dp
        ngto(0) = 26
        ngto(1) = 16
        ival(1) = 4
      CASE (12)
        cval = 2.12183620_dp
        aval = 0.01215850_dp
        ngto(0) = 26
        ngto(1) = 16
        ival(1) = 4
      CASE (13)
        cval = 2.06073230_dp
        aval = 0.01449350_dp
        ngto(0) = 26
        ngto(1) = 20
        ival(0) = 1
      CASE (14)
        cval = 2.08563660_dp
        aval = 0.01861460_dp
        ngto(0) = 26
        ngto(1) = 20
        ival(0) = 1
      CASE (15)
        cval = 2.04879270_dp
        aval = 0.02147790_dp
        ngto(0) = 26
        ngto(1) = 20
        ival(0) = 1
      CASE (16)
        cval = 2.06216660_dp
        aval = 0.01978920_dp
        ngto(0) = 26
        ngto(1) = 20
        ival(0) = 1
      CASE (17)
        cval = 2.04628670_dp
        aval = 0.02451470_dp
        ngto(0) = 26
        ngto(1) = 20
        ival(0) = 1
      CASE (18)
        cval = 2.08675200_dp
        aval = 0.02635040_dp
        ngto(0) = 26
        ngto(1) = 20
        ival(0) = 1
      CASE (19)
        cval = 2.02715220_dp
        aval = 0.01822040_dp
        ngto(0) = 29
        ngto(1) = 20
        ival(1) = 2
      CASE (20)
        cval = 2.01465650_dp
        aval = 0.01646570_dp
        ngto(0) = 29
        ngto(1) = 20
        ival(1) = 2
      CASE (21)
        cval = 2.01605240_dp
        aval = 0.01254190_dp
        ngto(0) = 30
        ngto(1) = 20
        ngto(2) = 18
        ival(1) = 2
      CASE (22)
        cval = 2.01800000_dp
        aval = 0.01195490_dp
        ngto(0) = 30
        ngto(1) = 21
        ngto(2) = 17
        ival(1) = 2
        ival(2) = 1
      CASE (23)
        cval = 1.98803560_dp
        aval = 0.02492140_dp
        ngto(0) = 30
        ngto(1) = 21
        ngto(2) = 17
        ival(1) = 2
        ival(2) = 1
      CASE (24)
        cval = 1.98984000_dp
        aval = 0.02568400_dp
        ngto(0) = 30
        ngto(1) = 21
        ngto(2) = 17
        ival(1) = 2
        ival(2) = 1
      CASE (25)
        cval = 2.01694380_dp
        aval = 0.02664480_dp
        ngto(0) = 30
        ngto(1) = 21
        ngto(2) = 17
        ival(1) = 2
        ival(2) = 1
      CASE (26)
        cval = 2.01824090_dp
        aval = 0.01355000_dp
        ngto(0) = 30
        ngto(1) = 21
        ngto(2) = 17
        ival(1) = 2
        ival(2) = 1
      CASE (27)
        cval = 1.98359400_dp
        aval = 0.01702210_dp
        ngto(0) = 30
        ngto(1) = 21
        ngto(2) = 17
        ival(1) = 2
        ival(2) = 2
      CASE (28)
        cval = 1.96797340_dp
        aval = 0.02163180_dp
        ngto(0) = 30
        ngto(1) = 22
        ngto(2) = 17
        ival(1) = 3
        ival(2) = 2
      CASE (29)
        cval = 1.98955180_dp
        aval = 0.02304480_dp
        ngto(0) = 30
        ngto(1) = 20
        ngto(2) = 17
        ival(1) = 3
        ival(2) = 2
      CASE (30)
        cval = 1.98074320_dp
        aval = 0.02754320_dp
        ngto(0) = 30
        ngto(1) = 21
        ngto(2) = 17
        ival(1) = 3
        ival(2) = 2
      CASE (31)
        cval = 2.00551070_dp
        aval = 0.02005530_dp
        ngto(0) = 30
        ngto(1) = 23
        ngto(2) = 17
        ival(0) = 1
        ival(2) = 2
      CASE (32)
        cval = 2.00000030_dp
        aval = 0.02003000_dp
        ngto(0) = 30
        ngto(1) = 24
        ngto(2) = 17
        ival(0) = 1
        ival(2) = 2
      CASE (33)
        cval = 2.00609100_dp
        aval = 0.02055620_dp
        ngto(0) = 30
        ngto(1) = 23
        ngto(2) = 17
        ival(0) = 1
        ival(2) = 2
      CASE (34)
        cval = 2.00701000_dp
        aval = 0.02230400_dp
        ngto(0) = 30
        ngto(1) = 24
        ngto(2) = 17
        ival(0) = 1
        ival(2) = 2
      CASE (35)
        cval = 2.01508710_dp
        aval = 0.02685790_dp
        ngto(0) = 30
        ngto(1) = 24
        ngto(2) = 17
        ival(0) = 1
        ival(2) = 2
      CASE (36)
        cval = 2.01960430_dp
        aval = 0.02960430_dp
        ngto(0) = 30
        ngto(1) = 24
        ngto(2) = 17
        ival(0) = 1
        ival(2) = 2
      CASE (37)
        cval = 2.00031000_dp
        aval = 0.00768400_dp
        ngto(0) = 32
        ngto(1) = 25
        ngto(2) = 17
        ival(0) = 1
        ival(1) = 1
        ival(2) = 4
      CASE (38)
        cval = 1.99563960_dp
        aval = 0.01401940_dp
        ngto(0) = 33
        ngto(1) = 24
        ngto(2) = 17
        ival(1) = 1
        ival(2) = 4
      CASE (39)
        cval = 1.98971210_dp
        aval = 0.01558470_dp
        ngto(0) = 33
        ngto(1) = 24
        ngto(2) = 20
        ival(1) = 1
      CASE (40)
        cval = 1.97976190_dp
        aval = 0.01705520_dp
        ngto(0) = 33
        ngto(1) = 24
        ngto(2) = 20
        ival(1) = 1
      CASE (41)
        cval = 1.97989290_dp
        aval = 0.01527040_dp
        ngto(0) = 33
        ngto(1) = 24
        ngto(2) = 20
        ival(1) = 1
      CASE (42)
        cval = 1.97909240_dp
        aval = 0.01879720_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 20
        ival(1) = 1
      CASE (43)
        cval = 1.98508430_dp
        aval = 0.01497550_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 20
        ival(1) = 2
        ival(2) = 1
      CASE (44)
        cval = 1.98515010_dp
        aval = 0.01856670_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 20
        ival(1) = 2
        ival(2) = 1
      CASE (45)
        cval = 1.98502970_dp
        aval = 0.01487000_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 20
        ival(1) = 2
        ival(2) = 1
      CASE (46)
        cval = 1.97672850_dp
        aval = 0.01762500_dp
        ngto(0) = 30
        ngto(1) = 24
        ngto(2) = 20
        ival(0) = 2
        ival(1) = 2
        ival(2) = 1
      CASE (47)
        cval = 1.97862730_dp
        aval = 0.01863310_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 20
        ival(1) = 2
        ival(2) = 1
      CASE (48)
        cval = 1.97990020_dp
        aval = 0.01347150_dp
        ngto(0) = 33
        ngto(1) = 24
        ngto(2) = 20
        ival(1) = 2
        ival(2) = 2
      CASE (49)
        cval = 1.97979410_dp
        aval = 0.00890265_dp
        ngto(0) = 33
        ngto(1) = 27
        ngto(2) = 20
        ival(0) = 2
        ival(2) = 2
      CASE (50)
        cval = 1.98001000_dp
        aval = 0.00895215_dp
        ngto(0) = 33
        ngto(1) = 27
        ngto(2) = 20
        ival(0) = 2
        ival(2) = 2
      CASE (51)
        cval = 1.97979980_dp
        aval = 0.01490290_dp
        ngto(0) = 33
        ngto(1) = 26
        ngto(2) = 20
        ival(1) = 1
        ival(2) = 2
      CASE (52)
        cval = 1.98009310_dp
        aval = 0.01490390_dp
        ngto(0) = 33
        ngto(1) = 26
        ngto(2) = 20
        ival(1) = 1
        ival(2) = 2
      CASE (53)
        cval = 1.97794750_dp
        aval = 0.01425880_dp
        ngto(0) = 33
        ngto(1) = 26
        ngto(2) = 20
        ival(0) = 2
        ival(1) = 1
        ival(2) = 2
      CASE (54)
        cval = 1.97784450_dp
        aval = 0.01430130_dp
        ngto(0) = 33
        ngto(1) = 26
        ngto(2) = 20
        ival(0) = 2
        ival(1) = 1
        ival(2) = 2
      CASE (55)
        cval = 1.97784450_dp
        aval = 0.00499318_dp
        ngto(0) = 32
        ngto(1) = 25
        ngto(2) = 17
        ival(0) = 1
        ival(1) = 3
        ival(2) = 6
      CASE (56)
        cval = 1.97764820_dp
        aval = 0.00500392_dp
        ngto(0) = 32
        ngto(1) = 25
        ngto(2) = 17
        ival(0) = 1
        ival(1) = 3
        ival(2) = 6
      CASE (57)
        cval = 1.97765150_dp
        aval = 0.00557083_dp
        ngto(0) = 32
        ngto(1) = 25
        ngto(2) = 20
        ival(0) = 1
        ival(1) = 3
        ival(2) = 3
      CASE (58)
        cval = 1.97768750_dp
        aval = 0.00547531_dp
        ngto(0) = 32
        ngto(1) = 25
        ngto(2) = 20
        ngto(3) = 16
        ival(0) = 1
        ival(1) = 3
        ival(2) = 3
        ival(3) = 3
      CASE (59)
        cval = 1.96986600_dp
        aval = 0.00813143_dp
        ngto(0) = 32
        ngto(1) = 25
        ngto(2) = 17
        ngto(3) = 16
        ival(0) = 1
        ival(1) = 3
        ival(2) = 6
        ival(3) = 4
      CASE (60)
        cval = 1.97765720_dp
        aval = 0.00489201_dp
        ngto(0) = 32
        ngto(1) = 25
        ngto(2) = 17
        ngto(3) = 16
        ival(0) = 1
        ival(1) = 3
        ival(2) = 6
        ival(3) = 4
      CASE (61)
        cval = 1.97768120_dp
        aval = 0.00499000_dp
        ngto(0) = 32
        ngto(1) = 25
        ngto(2) = 17
        ngto(3) = 16
        ival(0) = 1
        ival(1) = 3
        ival(2) = 6
        ival(3) = 4
      CASE (62)
        cval = 1.97745700_dp
        aval = 0.00615587_dp
        ngto(0) = 32
        ngto(1) = 25
        ngto(2) = 17
        ngto(3) = 16
        ival(0) = 1
        ival(1) = 3
        ival(2) = 6
        ival(3) = 4
      CASE (63)
        cval = 1.97570240_dp
        aval = 0.00769959_dp
        ngto(0) = 32
        ngto(1) = 25
        ngto(2) = 17
        ngto(3) = 16
        ival(0) = 1
        ival(1) = 3
        ival(2) = 6
        ival(3) = 4
      CASE (64)
        cval = 1.97629350_dp
        aval = 0.00706610_dp
        ngto(0) = 32
        ngto(1) = 25
        ngto(2) = 20
        ngto(3) = 16
        ival(0) = 1
        ival(1) = 3
        ival(2) = 3
        ival(3) = 4
      CASE (65)
        cval = 1.96900000_dp
        aval = 0.01019150_dp
        ngto(0) = 32
        ngto(1) = 26
        ngto(2) = 18
        ngto(3) = 16
        ival(0) = 1
        ival(1) = 3
        ival(2) = 6
        ival(3) = 4
      CASE (66)
        cval = 1.97350000_dp
        aval = 0.01334320_dp
        ngto(0) = 33
        ngto(1) = 26
        ngto(2) = 18
        ngto(3) = 16
        ival(0) = 1
        ival(1) = 3
        ival(2) = 6
        ival(3) = 4
      CASE (67)
        cval = 1.97493000_dp
        aval = 0.01331360_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 17
        ngto(3) = 14
        ival(1) = 2
        ival(2) = 5
        ival(3) = 4
      CASE (68)
        cval = 1.97597670_dp
        aval = 0.01434040_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 17
        ngto(3) = 14
        ival(0) = 0
        ival(1) = 2
        ival(2) = 5
        ival(3) = 4
      CASE (69)
        cval = 1.97809240_dp
        aval = 0.01529430_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 17
        ngto(3) = 14
        ival(0) = 0
        ival(1) = 2
        ival(2) = 5
        ival(3) = 4
      CASE (70)
        cval = 1.97644360_dp
        aval = 0.01312770_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 17
        ngto(3) = 14
        ival(0) = 0
        ival(1) = 2
        ival(2) = 5
        ival(3) = 4
      CASE (71)
        cval = 1.96998000_dp
        aval = 0.01745150_dp
        ngto(0) = 31
        ngto(1) = 24
        ngto(2) = 20
        ngto(3) = 14
        ival(0) = 1
        ival(1) = 2
        ival(2) = 2
        ival(3) = 4
      CASE (72)
        cval = 1.97223830_dp
        aval = 0.01639750_dp
        ngto(0) = 31
        ngto(1) = 24
        ngto(2) = 20
        ngto(3) = 14
        ival(0) = 1
        ival(1) = 2
        ival(2) = 2
        ival(3) = 4
      CASE (73)
        cval = 1.97462110_dp
        aval = 0.01603680_dp
        ngto(0) = 31
        ngto(1) = 24
        ngto(2) = 20
        ngto(3) = 14
        ival(0) = 1
        ival(1) = 2
        ival(2) = 2
        ival(3) = 4
      CASE (74)
        cval = 1.97756000_dp
        aval = 0.02030570_dp
        ngto(0) = 31
        ngto(1) = 24
        ngto(2) = 20
        ngto(3) = 14
        ival(0) = 1
        ival(1) = 2
        ival(2) = 2
        ival(3) = 4
      CASE (75)
        cval = 1.97645760_dp
        aval = 0.02057180_dp
        ngto(0) = 31
        ngto(1) = 24
        ngto(2) = 20
        ngto(3) = 14
        ival(0) = 1
        ival(1) = 2
        ival(2) = 2
        ival(3) = 4
      CASE (76)
        cval = 1.97725820_dp
        aval = 0.02058210_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 20
        ngto(3) = 15
        ival(0) = 0
        ival(1) = 2
        ival(2) = 2
        ival(3) = 4
      CASE (77)
        cval = 1.97749380_dp
        aval = 0.02219380_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 20
        ngto(3) = 15
        ival(0) = 0
        ival(1) = 2
        ival(2) = 2
        ival(3) = 4
      CASE (78)
        cval = 1.97946280_dp
        aval = 0.02216280_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 20
        ngto(3) = 15
        ival(0) = 0
        ival(1) = 2
        ival(2) = 2
        ival(3) = 4
      CASE (79)
        cval = 1.97852130_dp
        aval = 0.02168500_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 20
        ngto(3) = 15
        ival(0) = 0
        ival(1) = 2
        ival(2) = 2
        ival(3) = 4
      CASE (80)
        cval = 1.98045190_dp
        aval = 0.02177860_dp
        ngto(0) = 32
        ngto(1) = 24
        ngto(2) = 20
        ngto(3) = 15
        ival(0) = 0
        ival(1) = 2
        ival(2) = 2
        ival(3) = 4
      CASE (81)
        cval = 1.97000000_dp
        aval = 0.02275000_dp
        ngto(0) = 31
        ngto(1) = 25
        ngto(2) = 18
        ngto(3) = 13
        ival(0) = 1
        ival(1) = 0
        ival(2) = 3
        ival(3) = 6
      CASE (82)
        cval = 1.97713580_dp
        aval = 0.02317030_dp
        ngto(0) = 31
        ngto(1) = 27
        ngto(2) = 18
        ngto(3) = 13
        ival(0) = 1
        ival(1) = 0
        ival(2) = 3
        ival(3) = 6
      CASE (83)
        cval = 1.97537880_dp
        aval = 0.02672860_dp
        ngto(0) = 32
        ngto(1) = 27
        ngto(2) = 17
        ngto(3) = 13
        ival(0) = 1
        ival(1) = 0
        ival(2) = 3
        ival(3) = 6
      CASE (84)
        cval = 1.97545360_dp
        aval = 0.02745360_dp
        ngto(0) = 31
        ngto(1) = 27
        ngto(2) = 17
        ngto(3) = 13
        ival(0) = 1
        ival(1) = 0
        ival(2) = 3
        ival(3) = 6
      CASE (85)
        cval = 1.97338370_dp
        aval = 0.02616310_dp
        ngto(0) = 32
        ngto(1) = 27
        ngto(2) = 19
        ngto(3) = 13
        ival(0) = 1
        ival(1) = 0
        ival(2) = 3
        ival(3) = 6
      CASE (86)
        cval = 1.97294240_dp
        aval = 0.02429220_dp
        ngto(0) = 32
        ngto(1) = 27
        ngto(2) = 19
        ngto(3) = 13
        ival(0) = 1
        ival(1) = 0
        ival(2) = 3
        ival(3) = 6
      CASE (87:106) ! these numbers are a educated guess
        cval = 1.98000000_dp
        aval = 0.01400000_dp
        ngto(0) = 34
        ngto(1) = 28
        ngto(2) = 20
        ngto(3) = 15
        ival(0) = 0
        ival(1) = 0
        ival(2) = 3
        ival(3) = 6
    END SELECT

  END SUBROUTINE Clementi_geobas
! *****************************************************************************
!> \brief ...
!> \param element_symbol ...
!> \param basis ...
!> \param basis_set_name ...
!> \param basis_set_file ...
!> \param basis_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE read_basis_set(element_symbol,basis,basis_set_name,basis_set_file,&
                            basis_section,error)

    CHARACTER(LEN=*), INTENT(IN)             :: element_symbol
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    CHARACTER(LEN=*), INTENT(IN)             :: basis_set_name, basis_set_file
    TYPE(section_vals_type), POINTER         :: basis_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_basis_set', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: maxpri = 40, maxset = 20

    CHARACTER(len=20*default_string_length)  :: line_att
    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(LEN=LEN(basis_set_name))       :: bsname
    CHARACTER(LEN=LEN(basis_set_name)+2)     :: bsname2
    CHARACTER(LEN=LEN(element_symbol))       :: symbol
    CHARACTER(LEN=LEN(element_symbol)+2)     :: symbol2
    INTEGER                                  :: i, ierr, ii, ipgf, irep, &
                                                iset, ishell, j, k, lshell, &
                                                nj, nmin, ns, nset, strlen1, &
                                                strlen2
    INTEGER, DIMENSION(maxpri, maxset)       :: l
    INTEGER, DIMENSION(maxset)               :: lmax, lmin, n, npgf, nshell
    LOGICAL                                  :: failure, found, is_ok, match, &
                                                read_from_input
    REAL(dp)                                 :: expzet, gcca, prefac, zeta
    REAL(dp), &
      DIMENSION(maxpri, maxpri, maxset)      :: gcc
    REAL(dp), DIMENSION(maxpri, maxset)      :: zet
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(val_type), POINTER                  :: val

    failure = .FALSE.

    bsname = basis_set_name
    symbol = element_symbol
    irep   = 0

    nset = 0
    lmin = 0
    lmax = 0
    npgf = 0
    n = 0
    l = 0
    zet = 0._dp
    gcc = 0._dp

    read_from_input = .FALSE.
    CALL section_vals_get(basis_section,explicit=read_from_input, error=error)
    IF (read_from_input) THEN
      NULLIFY(list,val)
      CALL section_vals_list_get(basis_section,"_DEFAULT_KEYWORD_",list=list,error=error)
      CALL uppercase(symbol)
      CALL uppercase(bsname)
      is_ok=cp_sll_val_next(list,val,error=error)
      CPPrecondition(is_ok, cp_failure_level, routineP, error, failure)
      CALL val_get(val,c_val=line_att,error=error)
      READ(line_att,*) nset
      CPPrecondition(nset <= maxset, cp_failure_level, routineP, error, failure)
      DO iset=1,nset
        is_ok=cp_sll_val_next(list,val,error=error)
        CPPrecondition(is_ok, cp_failure_level, routineP, error, failure)
        CALL val_get(val,c_val=line_att,error=error)
        READ(line_att,*) n(iset)
        CALL remove_word(line_att)
        READ(line_att,*) lmin(iset)
        CALL remove_word(line_att)
        READ(line_att,*) lmax(iset)
        CALL remove_word(line_att)
        READ(line_att,*) npgf(iset)
        CALL remove_word(line_att)
        CPPrecondition(npgf(iset) <= maxpri, cp_failure_level, routineP, error, failure)
        nshell(iset) = 0
        DO lshell=lmin(iset),lmax(iset)
          nmin = n(iset) + lshell - lmin(iset)
          READ(line_att,*) ishell
          CALL remove_word(line_att)
          nshell(iset) = nshell(iset) + ishell
          DO i=1,ishell
            l(nshell(iset)-ishell+i,iset) = lshell
          END DO
        END DO
        CPPrecondition(LEN_TRIM(line_att)==0, cp_failure_level, routineP, error, failure)
        DO ipgf=1,npgf(iset)
          is_ok=cp_sll_val_next(list,val,error=error)
          CPPrecondition(is_ok, cp_failure_level, routineP, error, failure)
          CALL val_get(val,c_val=line_att,error=error)
          READ(line_att,*) zet(ipgf,iset),(gcc(ipgf,ishell,iset),ishell=1,nshell(iset))
        END DO
      END DO
    ELSE
      NULLIFY(parser)
      CALL parser_create(parser,basis_set_file,error=error)
      ! Search for the requested basis set in the basis set file
      ! until the basis set is found or the end of file is reached
      search_loop: DO
        CALL parser_search_string(parser,TRIM(bsname),.TRUE.,found,line,error=error)
        IF (found) THEN
          CALL uppercase(symbol)
          CALL uppercase(bsname)
          match = .FALSE.
          CALL uppercase(line)
          ! Check both the element symbol and the basis set name
          line2 = " "//line//" "
          symbol2 = " "//TRIM(symbol)//" "
          bsname2 = " "//TRIM(bsname)//" "
          strlen1 = LEN_TRIM(symbol2) + 1
          strlen2 = LEN_TRIM(bsname2) + 1

          IF ( (INDEX(line2,symbol2(:strlen1)) > 0).AND.&
               (INDEX(line2,bsname2(:strlen2)) > 0) ) match = .TRUE.

          IF (match) THEN
            ! Read the basis set information
            CALL parser_get_object(parser,nset,newline=.TRUE.,error=error)
            CPPrecondition(nset <= maxset, cp_failure_level, routineP, error, failure)
            DO iset=1,nset
              CALL parser_get_object(parser,n(iset),newline=.TRUE.,error=error)
              CALL parser_get_object(parser,lmin(iset),error=error)
              CALL parser_get_object(parser,lmax(iset),error=error)
              CALL parser_get_object(parser,npgf(iset),error=error)
              CPPrecondition(npgf(iset) <= maxpri, cp_failure_level, routineP, error, failure)
              nshell(iset) = 0
              DO lshell=lmin(iset),lmax(iset)
                nmin = n(iset) + lshell - lmin(iset)
                CALL parser_get_object(parser,ishell,error=error)
                nshell(iset) = nshell(iset) + ishell
                DO i=1,ishell
                  l(nshell(iset)-ishell+i,iset) = lshell
                END DO
              END DO
              DO ipgf=1,npgf(iset)
                CALL parser_get_object(parser,zet(ipgf,iset),newline=.TRUE.,error=error)
                DO ishell=1,nshell(iset)
                  CALL parser_get_object(parser,gcc(ipgf,ishell,iset),error=error)
                END DO
              END DO
            END DO

            EXIT search_loop

          END IF
        ELSE
          ! Stop program, if the end of file is reached
          CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
        END IF

      END DO search_loop

      CALL parser_release(parser,error=error)
    END IF

    ! fill in the basis data structures
    basis%nprim = 0
    basis%nbas = 0
    DO i=1,nset
      DO j=lmin(i),lmax(i)
        basis%nprim(j) = basis%nprim(j) + npgf(i)
      END DO
      DO j=1,nshell(i)
        k = l(j,i)
        basis%nbas(k) = basis%nbas(k)+1
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
            IF ( l(ii,i) == j ) THEN
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

    ! Normalization
    DO j=0,3
      expzet = 0.25_dp*REAL(2*j + 3,dp)
      prefac = SQRT( SQRT(pi)/2._dp**(j+2)*dfac(2*j+1) )
      DO ipgf=1,basis%nprim(j)
        DO ii=1,basis%nbas(j)
          gcca = basis%cm(ipgf,ii,j)
          zeta = 2._dp*basis%am(ipgf,j)
          basis%cm(ipgf,ii,j) = zeta**expzet*gcca/prefac
        END DO
      END DO
    END DO

  END SUBROUTINE read_basis_set

! *****************************************************************************
!> \brief ...
!> \param optimization ...
!> \param opt_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE read_atom_opt_section(optimization,opt_section,error)
    TYPE(atom_optimization_type), &
      INTENT(INOUT)                          :: optimization
    TYPE(section_vals_type), POINTER         :: opt_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_atom_opt_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: miter, ndiis
    REAL(KIND=dp)                            :: damp, eps_diis, eps_scf

    CALL section_vals_val_get(opt_section,"MAX_ITER",i_val=miter,error=error)
    CALL section_vals_val_get(opt_section,"EPS_SCF",r_val=eps_scf,error=error)
    CALL section_vals_val_get(opt_section,"N_DIIS",i_val=ndiis,error=error)
    CALL section_vals_val_get(opt_section,"EPS_DIIS",r_val=eps_diis,error=error)
    CALL section_vals_val_get(opt_section,"DAMPING",r_val=damp,error=error)

    optimization%max_iter = miter
    optimization%eps_scf  = eps_scf
    optimization%n_diis   = ndiis
    optimization%eps_diis = eps_diis
    optimization%damping  = damp

  END SUBROUTINE read_atom_opt_section
! *****************************************************************************
!> \brief ...
!> \param potential ...
!> \param potential_section ...
!> \param zval ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE init_atom_potential(potential,potential_section,zval,error)
    TYPE(atom_potential_type), INTENT(INOUT) :: potential
    TYPE(section_vals_type), POINTER         :: potential_section
    INTEGER, INTENT(IN)                      :: zval
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_atom_potential', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: pseudo_fn, pseudo_name
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:), POINTER          :: convals
    TYPE(section_vals_type), POINTER         :: gth_potential_section

    failure = .FALSE.

    IF ( zval > 0 ) THEN
      CALL section_vals_val_get(potential_section,"PSEUDO_TYPE",i_val=potential%ppot_type,error=error)

      SELECT CASE (potential%ppot_type)
        CASE (gth_pseudo)
           CALL section_vals_val_get(potential_section,"POTENTIAL_FILE_NAME",c_val=pseudo_fn,error=error)
           CALL section_vals_val_get(potential_section,"POTENTIAL_NAME",c_val=pseudo_name,error=error)
           gth_potential_section => section_vals_get_subs_vals(potential_section,"GTH_POTENTIAL",error=error)
           CALL read_gth_potential(ptable(zval)%symbol,potential%gth_pot,&
                 pseudo_name,pseudo_fn,gth_potential_section,error)
        CASE (no_pseudo)
          ! do nothing
        CASE DEFAULT
          CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      END SELECT
    ELSE
      potential%ppot_type = no_pseudo
    END IF

    ! confinement
    NULLIFY(convals)
    CALL section_vals_val_get(potential_section,"CONFINEMENT",r_vals=convals,error=error)
    IF ( SIZE (convals) >= 1 ) THEN
      IF ( convals(1) > 0.0_dp ) THEN
        potential%confinement = .TRUE.
        potential%acon = convals(1)
        IF ( SIZE (convals) >= 2 ) THEN
          potential%rcon = convals(2)
        ELSE
          potential%rcon = 4.0_dp
        END IF
        IF ( SIZE (convals) >= 3 ) THEN
          potential%scon = convals(3)
        ELSE
          potential%scon = 2.0_dp
        END IF
      ELSE
        potential%confinement = .FALSE.
      END IF
    ELSE
      potential%confinement = .FALSE.
    END IF

  END SUBROUTINE init_atom_potential
! *****************************************************************************
!> \brief ...
!> \param potential ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE release_atom_potential(potential,error)
    TYPE(atom_potential_type), INTENT(INOUT) :: potential
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_atom_potential', &
      routineP = moduleN//':'//routineN

    potential%confinement = .FALSE.

  END SUBROUTINE release_atom_potential
! *****************************************************************************
!> \brief ...
!> \param element_symbol ...
!> \param potential ...
!> \param pseudo_name ...
!> \param pseudo_file ...
!> \param potential_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE read_gth_potential(element_symbol,potential,pseudo_name,pseudo_file,&
                            potential_section,error)

    CHARACTER(LEN=*), INTENT(IN)             :: element_symbol
    TYPE(atom_gthpot_type), INTENT(INOUT)    :: potential
    CHARACTER(LEN=*), INTENT(IN)             :: pseudo_name, pseudo_file
    TYPE(section_vals_type), POINTER         :: potential_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_gth_potential', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(len=5*default_string_length)   :: line_att
    CHARACTER(LEN=LEN(element_symbol))       :: symbol
    CHARACTER(LEN=LEN(element_symbol)+2)     :: symbol2
    CHARACTER(LEN=LEN(pseudo_name))          :: apname
    CHARACTER(LEN=LEN(pseudo_name)+2)        :: apname2
    INTEGER                                  :: i, ic, ipot, j, l, nlmax, &
                                                strlen1, strlen2
    INTEGER, DIMENSION(0:4)                  :: elec_conf
    LOGICAL                                  :: failure, found, is_ok, match, &
                                                read_from_input
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(val_type), POINTER                  :: val

    failure = .FALSE.

    elec_conf = 0

    apname = pseudo_name
    symbol = element_symbol

    potential%symbol = symbol
    potential%pname = apname
    potential%econf = 0
    potential%rc = 0._dp
    potential%ncl = 0
    potential%cl = 0._dp
    potential%nl = 0
    potential%rcnl = 0._dp
    potential%hnl = 0._dp

    potential%lpotextended = .FALSE.
    potential%lsdpot = .FALSE.
    potential%nlcc = .FALSE.
    potential%nexp_lpot = 0
    potential%nexp_lsd = 0
    potential%nexp_nlcc = 0

    read_from_input = .FALSE.
    CALL section_vals_get(potential_section,explicit=read_from_input, error=error)
    IF (read_from_input) THEN
      CALL section_vals_list_get(potential_section,"_DEFAULT_KEYWORD_",list=list,error=error)
      CALL uppercase(symbol)
      CALL uppercase(apname)
      ! Read the electronic configuration, not used here
      l = 0
      is_ok=cp_sll_val_next(list,val,error=error)
      CPPostcondition(is_ok, cp_failure_level, routineP, error, failure)
      CALL val_get(val,c_val=line_att,error=error)
      READ(line_att,*) elec_conf(l)
      CALL remove_word(line_att)
      DO WHILE (LEN_TRIM(line_att) /= 0)
        l = l + 1
        READ(line_att,*)elec_conf(l)
        CALL remove_word(line_att)
      END DO
      potential%econf(0:3) = elec_conf(0:3)
      potential%zion = REAL ( SUM(elec_conf), dp )
      ! Read r(loc) to define the exponent of the core charge
      is_ok=cp_sll_val_next(list,val,error=error)
      CPPostcondition(is_ok, cp_failure_level, routineP, error, failure)
      CALL val_get(val,c_val=line_att,error=error)
      READ(line_att,*) potential%rc
      CALL remove_word(line_att)
      ! Read the parameters for the local part of the GTH pseudopotential (ppl)
      READ(line_att,*) potential%ncl
      CALL remove_word(line_att)
      DO i=1,potential%ncl
        READ(line_att,*) potential%cl(i)
        CALL remove_word(line_att)
      END DO
      ! Check for the next entry: LPOT, NLCC, LSD, or ppnl
      DO
        is_ok=cp_sll_val_next(list,val,error=error)
        CPPostcondition(is_ok, cp_failure_level, routineP, error, failure)
        CALL val_get(val,c_val=line_att,error=error)
        IF(INDEX(line_att,"LPOT") /= 0) THEN
          potential%lpotextended = .TRUE.
          CALL remove_word(line_att)
          READ(line_att,*) potential%nexp_lpot
          DO ipot=1,potential%nexp_lpot
             is_ok=cp_sll_val_next(list,val,error=error)
             CPPostcondition(is_ok, cp_failure_level, routineP, error, failure)
             CALL val_get(val,c_val=line_att,error=error)
             READ(line_att,*) potential%alpha_lpot(ipot)
             CALL remove_word(line_att)
             READ(line_att,*) potential%nct_lpot(ipot)
             CALL remove_word(line_att)
             DO ic=1,potential%nct_lpot(ipot)
               READ(line_att,*) potential%cval_lpot(ic,ipot)
               CALL remove_word(line_att)
             END DO
          END DO
        ELSEIF(INDEX(line_att,"NLCC") /= 0) THEN
          potential%nlcc = .TRUE.
          CALL remove_word(line_att)
          READ(line_att,*) potential%nexp_nlcc
          DO ipot=1,potential%nexp_nlcc
             is_ok=cp_sll_val_next(list,val,error=error)
             CPPostcondition(is_ok, cp_failure_level, routineP, error, failure)
             CALL val_get(val,c_val=line_att,error=error)
             READ(line_att,*) potential%alpha_nlcc(ipot)
             CALL remove_word(line_att)
             READ(line_att,*) potential%nct_nlcc(ipot)
             CALL remove_word(line_att)
             DO ic=1,potential%nct_nlcc(ipot)
               READ(line_att,*) potential%cval_nlcc(ic,ipot)
               !make cp2k compatible with bigdft 
               potential%cval_nlcc(ic,ipot)=potential%cval_nlcc(ic,ipot)/(4.0_dp*pi)
               CALL remove_word(line_att)
             END DO
          END DO
        ELSEIF(INDEX(line_att,"LSD") /= 0) THEN
          potential%lsdpot = .TRUE.
          CALL remove_word(line_att)
          READ(line_att,*) potential%nexp_lsd
          DO ipot=1,potential%nexp_lsd
             is_ok=cp_sll_val_next(list,val,error=error)
             CPPostcondition(is_ok, cp_failure_level, routineP, error, failure)
             CALL val_get(val,c_val=line_att,error=error)
             READ(line_att,*) potential%alpha_lsd(ipot)
             CALL remove_word(line_att)
             READ(line_att,*) potential%nct_lsd(ipot)
             CALL remove_word(line_att)
             DO ic=1,potential%nct_lsd(ipot)
               READ(line_att,*) potential%cval_lsd(ic,ipot)
               CALL remove_word(line_att)
             END DO
          END DO
        ELSE
          EXIT
        END IF
      END DO
      ! Read the parameters for the non-local part of the GTH pseudopotential (ppnl)
      READ(line_att,*) nlmax
      CALL remove_word(line_att)
      IF (nlmax > 0) THEN
        ! Load the parameter for nlmax non-local projectors
        DO l=0,nlmax-1
          is_ok=cp_sll_val_next(list,val,error=error)
          CPPostcondition(is_ok, cp_failure_level, routineP, error, failure)
          CALL val_get(val,c_val=line_att,error=error)
          READ(line_att,*) potential%rcnl(l)
          CALL remove_word(line_att)
          READ(line_att,*) potential%nl(l)
          CALL remove_word(line_att)
          DO i=1,potential%nl(l)
            IF (i == 1) THEN
              READ(line_att,*) potential%hnl(1,1,l)
              CALL remove_word(line_att)
            ELSE
              CPPostcondition(LEN_TRIM(line_att)==0, cp_failure_level, routineP, error, failure)
              is_ok=cp_sll_val_next(list,val,error=error)
              CPPostcondition(is_ok, cp_failure_level, routineP, error, failure)
              CALL val_get(val,c_val=line_att,error=error)
              READ(line_att,*) potential%hnl(i,i,l)
              CALL remove_word(line_att)
            END IF
            DO j=i+1,potential%nl(l)
              READ(line_att,*) potential%hnl(i,j,l)
              potential%hnl(j,i,l) = potential%hnl(i,j,l)
              CALL remove_word(line_att)
            END DO
          END DO
          CPPostcondition(LEN_TRIM(line_att)==0, cp_failure_level, routineP, error, failure)
        END DO
      END IF
    ELSE
      NULLIFY(parser)
      CALL parser_create(parser,pseudo_file,error=error)

      search_loop: DO
        CALL parser_search_string(parser,TRIM(apname),.TRUE.,found,line,error=error)
        IF (found) THEN
          CALL uppercase(symbol)
          CALL uppercase(apname)
          ! Check both the element symbol and the atomic potential name
          match = .FALSE.
          CALL uppercase(line)
          line2 = " "//line//" "
          symbol2 = " "//TRIM(symbol)//" "
          apname2 = " "//TRIM(apname)//" "
          strlen1 = LEN_TRIM(symbol2) + 1
          strlen2 = LEN_TRIM(apname2) + 1

          IF ( (INDEX(line2,symbol2(:strlen1)) > 0).AND.&
                  (INDEX(line2,apname2(:strlen2)) > 0) ) match = .TRUE.

          IF (match) THEN
             ! Read the electronic configuration
             l = 0
             CALL parser_get_object(parser,elec_conf(l),newline=.TRUE.,error=error)
             DO WHILE (parser_test_next_token(parser,error=error) == "INT")
               l = l + 1
               CALL parser_get_object(parser,elec_conf(l),error=error)
             END DO
             potential%econf(0:3) = elec_conf(0:3)
             potential%zion = REAL ( SUM(elec_conf), dp )
             ! Read r(loc) to define the exponent of the core charge
             CALL parser_get_object(parser,potential%rc,newline=.TRUE.,error=error)
             ! Read the parameters for the local part of the GTH pseudopotential (ppl)
             CALL parser_get_object(parser,potential%ncl,error=error)
             DO i=1,potential%ncl
               CALL parser_get_object(parser,potential%cl(i),error=error)
             END DO
             ! Extended type input
             DO
               CALL parser_get_next_line(parser,1,error=error)
               IF(parser_test_next_token(parser,error=error) == "INT") THEN
                 EXIT
               ELSEIF(parser_test_next_token(parser,error=error) == "STR") THEN
                 CALL parser_get_object(parser,line,error=error)
                 IF(INDEX(LINE,"LPOT") /= 0) THEN
                   ! local potential
                   potential%lpotextended = .TRUE.
                   CALL parser_get_object(parser,potential%nexp_lpot,error=error)
                   DO ipot=1,potential%nexp_lpot
                     CALL parser_get_object(parser,potential%alpha_lpot(ipot),newline=.TRUE.,error=error)
                     CALL parser_get_object(parser,potential%nct_lpot(ipot),error=error)
                     DO ic=1,potential%nct_lpot(ipot)
                       CALL parser_get_object(parser,potential%cval_lpot(ic,ipot),error=error)
                     END DO
                   END DO
                 ELSEIF(INDEX(LINE,"NLCC") /= 0) THEN
                   ! NLCC
                   potential%nlcc = .TRUE.
                   CALL parser_get_object(parser,potential%nexp_nlcc,error=error)
                   DO ipot=1,potential%nexp_nlcc
                     CALL parser_get_object(parser,potential%alpha_nlcc(ipot),newline=.TRUE.,error=error)
                     CALL parser_get_object(parser,potential%nct_nlcc(ipot),error=error)
                     DO ic=1,potential%nct_nlcc(ipot)
                       CALL parser_get_object(parser,potential%cval_nlcc(ic,ipot),error=error)
                       !make cp2k compatible with bigdft 
                       potential%cval_nlcc(ic,ipot)=potential%cval_nlcc(ic,ipot)/(4.0_dp*pi)
                     END DO
                   END DO
                 ELSEIF(INDEX(LINE,"LSD") /= 0) THEN
                   ! LSD potential
                   potential%lsdpot = .TRUE.
                   CALL parser_get_object(parser,potential%nexp_lsd,error=error)
                   DO ipot=1,potential%nexp_lsd
                     CALL parser_get_object(parser,potential%alpha_lsd(ipot),newline=.TRUE.,error=error)
                     CALL parser_get_object(parser,potential%nct_lsd(ipot),error=error)
                     DO ic=1,potential%nct_lsd(ipot)
                       CALL parser_get_object(parser,potential%cval_lsd(ic,ipot),error=error)
                     END DO
                   END DO
                 ELSE
                   CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
                 END IF
               ELSE
                 CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
               END IF
             END DO
             ! Read the parameters for the non-local part of the GTH pseudopotential (ppnl)
             CALL parser_get_object(parser,nlmax,error=error)
             IF (nlmax > 0) THEN
                ! Load the parameter for n non-local projectors
                DO l=0,nlmax-1
                   CALL parser_get_object(parser,potential%rcnl(l),newline=.TRUE.,error=error)
                   CALL parser_get_object(parser,potential%nl(l),error=error)
                   DO i=1,potential%nl(l)
                      IF (i == 1) THEN
                        CALL parser_get_object(parser,potential%hnl(i,i,l),error=error)
                      ELSE
                        CALL parser_get_object(parser,potential%hnl(i,i,l),newline=.TRUE.,error=error)
                      END IF
                      DO j=i+1,potential%nl(l)
                        CALL parser_get_object(parser,potential%hnl(i,j,l),error=error)
                        potential%hnl(j,i,l) = potential%hnl(i,j,l)
                      END DO
                   END DO
                END DO
             END IF
             EXIT search_loop
          END IF
        ELSE
          ! Stop program, if the end of file is reached
          CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
        END IF

      END DO search_loop

      CALL parser_release(parser,error=error)
    END IF

  END SUBROUTINE read_gth_potential
! *****************************************************************************

END MODULE atom_types
