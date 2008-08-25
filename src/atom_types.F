!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Define the atom type and its sub types
!> \author  jgh
!> \date    03.03.2008
!> \version 1.0
!>
! *****************************************************************************
MODULE atom_types

  USE input_constants,                 ONLY: contracted_gto,&
                                             gaussian,&
                                             geometrical_gto,&
                                             numerical,&
                                             slater
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                             create_grid_atom,&
                                             deallocate_grid_atom,&
                                             grid_atom_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_types'

  INTEGER, PARAMETER                                 :: GTO_BASIS=100,&
                                                        CGTO_BASIS=101,&
                                                        STO_BASIS=102,&
                                                        NUM_BASIS=103

!> Provides all information about a basis set
! *****************************************************************************
  TYPE atom_basis_type
    INTEGER                                       :: basis_type
    INTEGER, DIMENSION(0:3)                       :: nbas
    INTEGER                                       :: ngto
    INTEGER, DIMENSION(0:3)                       :: nprim
    REAL(KIND=dp),DIMENSION(:),POINTER            :: al         !GTO exponents
    REAL(KIND=dp),DIMENSION(:,:),POINTER          :: am         !CGTO exponents
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: cm         !Contraction coeffs
    REAL(KIND=dp),DIMENSION(:),POINTER            :: as         !STO exponents
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: bf         !num. bsf
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: dbf        !derivatives (num)
    TYPE(grid_atom_type),POINTER                  :: grid
  END TYPE atom_basis_type

!> Provides all information on states and occupation
! *****************************************************************************
  TYPE atom_state
    REAL(KIND=dp),DIMENSION(0:3,10)               :: occ
    REAL(KIND=dp),DIMENSION(0:3,10)               :: core
    INTEGER                                       :: maxl_occ
    INTEGER,DIMENSION(0:3)                        :: maxn_occ
    INTEGER                                       :: maxl_calc
    INTEGER,DIMENSION(0:3)                        :: maxn_calc
  END TYPE atom_state

!> Holds atomic integrals
! *****************************************************************************
  TYPE eri
    REAL(KIND=dp),DIMENSION(:,:),POINTER          :: int
  END TYPE eri

  TYPE atom_integrals
    INTEGER                                       :: status
    LOGICAL                                       :: eri_coulomb
    LOGICAL                                       :: eri_exchange
    LOGICAL                                       :: all_nu
    INTEGER, DIMENSION(0:3)                       :: n
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: ovlp, kin, core, proj
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: conf
    TYPE(eri),DIMENSION(20)                       :: ceri
    TYPE(eri),DIMENSION(20)                       :: eeri
  END TYPE atom_integrals

!> Holds atomic orbitals and energies
! *****************************************************************************
  TYPE atom_orbitals 
    INTEGER, DIMENSION(0:3)                       :: n
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: wfn
    REAL(KIND=dp),DIMENSION(:,:,:),POINTER        :: pmat
    REAL(KIND=dp),DIMENSION(:,:),POINTER          :: ener
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
    REAL(KIND=dp)                                 :: ecore
    REAL(KIND=dp)                                 :: epseudo
    REAL(KIND=dp)                                 :: eploc
    REAL(KIND=dp)                                 :: epnl
    REAL(KIND=dp)                                 :: exc
    REAL(KIND=dp)                                 :: ecoulomb
    REAL(KIND=dp)                                 :: eexchange
  END TYPE atom_energy_type

!> Provides all information about an atomic kind
! *****************************************************************************
  TYPE atom_type
    INTEGER                                       :: z
    INTEGER                                       :: zcore
    INTEGER                                       :: ppot_type
    TYPE(atom_basis_type), POINTER                :: basis
    TYPE(atom_state), POINTER                     :: state
    TYPE(atom_integrals), POINTER                 :: integrals
    TYPE(atom_orbitals), POINTER                  :: orbitals
    TYPE(atom_energy_type)                        :: energy
  END TYPE atom_type
! *****************************************************************************
  TYPE atom_p_type
    TYPE(atom_type), POINTER                      :: atom
  END TYPE atom_p_type

  PUBLIC :: atom_p_type, atom_type, atom_basis_type, atom_state, atom_integrals
  PUBLIC :: atom_orbitals, atom_energy_type, eri
  PUBLIC :: create_atom_type, release_atom_type, set_atom
  PUBLIC :: create_atom_orbs, release_atom_orbs
  PUBLIC :: init_atom_basis, release_atom_basis
  PUBLIC :: GTO_BASIS, CGTO_BASIS, STO_BASIS, NUM_BASIS
  PUBLIC :: opmat_type, create_opmat, release_opmat
  PUBLIC :: opgrid_type, create_opgrid, release_opgrid

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE init_atom_basis(basis,basis_section,error)
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    TYPE(section_vals_type), POINTER         :: basis_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_atom_basis', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nu = 40
    REAL(KIND=dp), DIMENSION(nu), PARAMETER :: ugbs = (/0.007299_dp, &
      0.013705_dp, 0.025733_dp, 0.048316_dp, 0.090718_dp, 0.170333_dp,&
      0.319819_dp, 0.600496_dp, 1.127497_dp, 2.117000_dp, 3.974902_dp, &
      1.463317_dp,14.013204_dp, 26.311339_dp, 49.402449_dp, 92.758561_dp, &
      174.164456_dp, 327.013024_dp,614.003114_dp, 1152.858743_dp, &
      2164.619772_dp, 4064.312984_dp, 7631.197056_dp,14328.416324_dp, &
      26903.186074_dp, 50513.706789_dp, 94845.070265_dp, 178082.107320_dp,&
      334368.848683_dp, 627814.487663_dp, 1178791.123851_dp, 2213310.684886_dp&
      , 4155735.557141_dp,7802853.046713_dp, 14650719.428954_dp, &
      27508345.793637_dp, 51649961.080194_dp, 96978513.342764_dp,&
      182087882.613702_dp, 341890134.751331_dp /)

    INTEGER                                  :: basistype, i, ierr, k, l, m, &
                                                ngp, nr, num_gto, quadtype
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: al, ear, rk

!   INTEGER, PARAMETER                       :: nu = 4
!   REAL(KIND=dp), DIMENSION(nu), PARAMETER :: ugbs = (/0.100000_dp, 0.50000_dp, 1.500000_dp, 3.00000_dp/)
!   INTEGER, PARAMETER                       :: nu = 2
!   REAL(KIND=dp), DIMENSION(nu), PARAMETER :: ugbs = (/0.50000_dp, 1.500000_dp/)
!   INTEGER, PARAMETER                       :: nu = 1
!   REAL(KIND=dp), DIMENSION(nu), PARAMETER :: ugbs = (/0.30000_dp/)
!
! Highly accurate relativistic universal Gaussian basis set: Dirac-Fock-Coulomb calculations 
! for atomic systems up to nobelium
! J. Chem. Phys. 101, 6829 (1994); DOI:10.1063/1.468311
! G. L. Malli and A. B. F. Da Silva
! Department of Chemistry, Simon Fraser University, Burnaby, B.C., Canada
! Yasuyuki Ishikawa
! Department of Chemistry, University of Puerto Rico, San Juan, Puerto Rico
!
! A universal Gaussian basis set is developed that leads to relativistic Dirac-Fock SCF energies 
! of comparable accuracy as that obtained by the accurate numerical finite-difference method 
! (GRASP2 package) [J. Phys. B 25, 1 (1992)]. The Gaussian-type functions of our universal basis 
! set satisfy the relativistic boundary conditions associated with the finite nuclear model for a 
! finite speed of light and conform to the so-called kinetic balance at the nonrelativistic limit. 
! We attribute the exceptionally high accuracy obtained in our calculations to the fact that the 
! representation of the relativistic dynamics of an electron in a spherical ball finite nucleus 
! near the origin in terms of our universal Gaussian basis set is as accurate as that provided by 
! the numerical finite-difference method. Results of the Dirac-Fock-Coulomb energies for a number 
! of atoms up to No (Z=102) and some negative ions are presented and compared with the recent 
! results obtained with the numerical finite-difference method and geometrical Gaussian basis sets 
! by Parpia, Mohanty, and Clementi [J. Phys. B 25, 1 (1992)]. The accuracy of our calculations is 
! estimated to be within a few parts in 109 for all the atomic systems studied.
!
!

    failure = .FALSE.

    NULLIFY(basis%al,basis%am,basis%cm,basis%as,basis%bf,basis%dbf)
    ! get information on quadrature type and number of grid points
    ! allocate and initialize the atomic grid 
    CALL allocate_grid_atom(basis%grid,error)
    CALL section_vals_val_get(basis_section,"QUADRATURE",i_val=quadtype,error=error)
    CALL section_vals_val_get(basis_section,"GRID_POINTS",i_val=ngp,error=error)
    CALL cp_assert(ngp > 0,cp_failure_level,cp_assertion_failed,routineP,&
                   "# point radial grid < 0",error,failure)
    CALL create_grid_atom(basis%grid,ngp,1,1,quadtype)
    basis%grid%nr = ngp

    CALL section_vals_val_get(basis_section,"BASIS_TYPE",i_val=basistype,error=error)
    SELECT CASE (basistype)
      CASE DEFAULT
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      CASE (gaussian)
         basis%basis_type = GTO_BASIS
         CALL section_vals_val_get(basis_section,"NUM_GTO",i_val=num_gto,error=error)
         IF ( num_gto < 1 ) THEN
           ! use default basis
           basis % nbas = nu
           basis % ngto = nu
           basis % nprim = 0
           ALLOCATE (basis%al(nu),STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
           basis%al(1:nu) = ugbs(1:nu)
         ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
         END IF
         ! initialize basis function on a radial grid
         nr = basis%grid%nr
         m  = MAXVAL(basis%nbas)
         ALLOCATE (basis%bf(nr,m,0:3),STAT=ierr)
         CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
         ALLOCATE (basis%dbf(nr,m,0:3),STAT=ierr)
         CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
         basis%bf =  0._dp
         DO l=0,3
           DO i=1,basis%nbas(l)
             al  = basis%al(i)
             DO k=1,nr
               rk  = basis%grid%rad(k)
               ear = EXP(-basis%al(i)*basis%grid%rad(k)**2)
               basis%bf(k,i,l) = rk**l * ear
               basis%dbf(k,i,l) = ( REAL(l,dp)*rk**(l-1) - 2._dp*al*rk**(l+1) ) * ear
             END DO
           END DO
         END DO
      CASE (geometrical_gto)
         basis%basis_type = GTO_BASIS
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      CASE (contracted_gto)
         basis%basis_type = CGTO_BASIS
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      CASE (slater)
         basis%basis_type = STO_BASIS
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      CASE (numerical)
         basis%basis_type = NUM_BASIS
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE init_atom_basis

  SUBROUTINE release_atom_basis(basis,error)
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_atom_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF(ASSOCIATED(basis%al)) THEN
       DEALLOCATE (basis%al,STAT=stat)
       CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF
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

    atom%ppot_type = 0

  END SUBROUTINE create_atom_type

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

  SUBROUTINE set_atom(atom,basis,state,integrals,orbitals,error)
    TYPE(atom_type), POINTER                 :: atom
    TYPE(atom_basis_type), INTENT(IN), &
      OPTIONAL, TARGET                       :: basis
    TYPE(atom_state), INTENT(IN), OPTIONAL, &
      TARGET                                 :: state
    TYPE(atom_integrals), INTENT(IN), &
      OPTIONAL, TARGET                       :: integrals
    TYPE(atom_orbitals), INTENT(IN), &
      OPTIONAL, TARGET                       :: orbitals
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
    END IF

  END SUBROUTINE set_atom

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
    
    ALLOCATE(orbs%wfn(mbas,mo,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    orbs%wfn = 0._dp
    
    ALLOCATE(orbs%pmat(mbas,mbas,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    orbs%pmat = 0._dp

    ALLOCATE(orbs%ener(mo,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    orbs%ener = 0._dp

  END SUBROUTINE create_atom_orbs
  
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
      DEALLOCATE(orbs%wfn,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%wfn)
    END IF
    IF(ASSOCIATED(orbs%pmat)) THEN
      DEALLOCATE(orbs%pmat,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%pmat)
    END IF
    IF(ASSOCIATED(orbs%ener)) THEN
      DEALLOCATE(orbs%ener,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      NULLIFY(orbs%ener)
    END IF

    DEALLOCATE(orbs,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    NULLIFY(orbs)

  END SUBROUTINE release_atom_orbs

! *****************************************************************************
  SUBROUTINE create_opmat(opmat,n,error)
    TYPE(opmat_type), POINTER                :: opmat
    INTEGER, DIMENSION(0:3), INTENT(IN)      :: n
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_opmat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, m
    LOGICAL                                  :: failure = .FALSE.

    m=MAXVAL(n)

    CPPrecondition(.NOT.ASSOCIATED(opmat), cp_failure_level, routineP, error, failure)

    ALLOCATE(opmat,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    opmat%n = n
    ALLOCATE(opmat%op(m,m,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    opmat%op = 0._dp

  END SUBROUTINE create_opmat

  SUBROUTINE release_opmat(opmat,error)
    TYPE(opmat_type), POINTER                :: opmat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_opmat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure = .FALSE.

    CPPrecondition(ASSOCIATED(opmat), cp_failure_level, routineP, error, failure)

    opmat%n = 0
    DEALLOCATE(opmat%op,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DEALLOCATE(opmat,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE release_opmat

! *****************************************************************************
  SUBROUTINE create_opgrid(opgrid,grid,error)
    TYPE(opgrid_type), POINTER               :: opgrid
    TYPE(grid_atom_type), POINTER            :: grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_opgrid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, nr
    LOGICAL                                  :: failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(opgrid), cp_failure_level, routineP, error, failure)

    ALLOCATE(opgrid,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    opgrid%grid => grid

    nr = grid%nr

    ALLOCATE(opgrid%op(nr),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    opgrid%op = 0._dp

  END SUBROUTINE create_opgrid

  SUBROUTINE release_opgrid(opgrid,error)
    TYPE(opgrid_type), POINTER               :: opgrid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_opgrid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure = .FALSE.

    CPPrecondition(ASSOCIATED(opgrid), cp_failure_level, routineP, error, failure)

    NULLIFY(opgrid%grid)
    DEALLOCATE(opgrid%op,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DEALLOCATE(opgrid,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE release_opgrid

! *****************************************************************************

END MODULE atom_types
