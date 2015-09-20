!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - 02.2004 flexible normalization of basis sets [jgh]
!>      - 07.2014 Add a set of contraction coefficient that only work on active 
!>                functions
!> \author Matthias Krack (04.07.2000)
! *****************************************************************************
MODULE basis_set_types

  USE ai_coulomb,                      ONLY: coulomb2
  USE bibliography,                    ONLY: VandeVondele2007,&
                                             cite_reference
  USE cp_linked_list_val,              ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_object,&
                                             parser_search_string
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE input_section_types,             ONLY: section_vals_list_get,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: dfac,&
                                             pi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: coset,&
                                             indco,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset,&
                                             nso,&
                                             nsoset
  USE orbital_symbols,                 ONLY: cgf_symbol,&
                                             sgf_symbol
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE sto_ng,                          ONLY: get_sto_ng
  USE string_utilities,                ONLY: integer_to_string,&
                                             remove_word,&
                                             uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'basis_set_types'

! *****************************************************************************
  ! Define the Gaussian-type orbital basis set type

  TYPE gto_basis_set_type
     !MK PRIVATE
     CHARACTER(LEN=default_string_length)       :: name
     REAL(KIND = dp)                            :: kind_radius
     REAL(KIND = dp)                            :: short_kind_radius
     INTEGER                                    :: norm_type = -1
     INTEGER                                    :: ncgf,nset,nsgf
     CHARACTER(LEN=12), DIMENSION(:), POINTER   :: cgf_symbol
     CHARACTER(LEN=6), DIMENSION(:), POINTER    :: sgf_symbol
     REAL(KIND = dp), DIMENSION(:), POINTER     :: norm_cgf,set_radius
     INTEGER, DIMENSION(:), POINTER             :: lmax,lmin,lx,ly,lz,m,ncgf_set,&
          npgf,nsgf_set,nshell
     REAL(KIND = dp), DIMENSION(:,:), POINTER   :: cphi,pgf_radius,sphi,scon,zet
     INTEGER, DIMENSION(:,:), POINTER           :: first_cgf,first_sgf,l,&
          last_cgf,last_sgf,n
     REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: gcc
  END TYPE gto_basis_set_type

  TYPE gto_basis_set_p_type
     TYPE(gto_basis_set_type), POINTER :: gto_basis_set
  END TYPE gto_basis_set_p_type

! *****************************************************************************
  ! Define the Slater-type orbital basis set type

  TYPE sto_basis_set_type
     PRIVATE
     CHARACTER(LEN=default_string_length)       :: name
     INTEGER                                    :: nshell
     CHARACTER(LEN=6), DIMENSION(:), POINTER    :: symbol
     INTEGER, DIMENSION(:), POINTER             :: nq,lq
     REAL(KIND = dp), DIMENSION(:), POINTER     :: zet
  END TYPE sto_basis_set_type

  TYPE sto_basis_set_p_type
     TYPE(sto_basis_set_type), POINTER :: sto_basis_set
  END TYPE sto_basis_set_p_type

! *****************************************************************************
  INTERFACE read_gto_basis_set
     MODULE PROCEDURE read_gto_basis_set1,read_gto_basis_set2
  END INTERFACE
! *****************************************************************************


  ! Public subroutines
  PUBLIC :: allocate_gto_basis_set,&
            deallocate_gto_basis_set,&
            get_gto_basis_set,&
            init_aux_basis_set,&
            init_cphi_and_sphi,&
            init_orb_basis_set,&
            read_gto_basis_set,&
            copy_gto_basis_set,&
            combine_basis_sets,&
            set_gto_basis_set,&
            write_gto_basis_set,&
            write_orb_basis_set

  PUBLIC :: allocate_sto_basis_set,&
            create_gto_from_sto_basis,&
            deallocate_sto_basis_set,&
            set_sto_basis_set,&
            srules

  ! Public data types
  PUBLIC :: gto_basis_set_p_type,&
            gto_basis_set_type,&
            sto_basis_set_type

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
! *****************************************************************************
  SUBROUTINE allocate_gto_basis_set(gto_basis_set)

    ! Allocate a Gaussian-type orbital (GTO) basis set data set.

    ! - Creation (26.10.2000,MK)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_gto_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL deallocate_gto_basis_set(gto_basis_set)

    ALLOCATE (gto_basis_set,STAT=stat)
    CPASSERT(stat==0)

    NULLIFY (gto_basis_set%cgf_symbol)
    NULLIFY (gto_basis_set%first_cgf)
    NULLIFY (gto_basis_set%first_sgf)
    NULLIFY (gto_basis_set%gcc)
    NULLIFY (gto_basis_set%l)
    NULLIFY (gto_basis_set%last_cgf)
    NULLIFY (gto_basis_set%last_sgf)
    NULLIFY (gto_basis_set%lmax)
    NULLIFY (gto_basis_set%lmin)
    NULLIFY (gto_basis_set%lx)
    NULLIFY (gto_basis_set%ly)
    NULLIFY (gto_basis_set%lz)
    NULLIFY (gto_basis_set%m)
    NULLIFY (gto_basis_set%n)
    NULLIFY (gto_basis_set%ncgf_set)
    NULLIFY (gto_basis_set%norm_cgf)
    NULLIFY (gto_basis_set%npgf)
    NULLIFY (gto_basis_set%nsgf_set)
    NULLIFY (gto_basis_set%nshell)
    NULLIFY (gto_basis_set%pgf_radius)
    NULLIFY (gto_basis_set%cphi)
    NULLIFY (gto_basis_set%sphi)
    NULLIFY (gto_basis_set%scon)
    NULLIFY (gto_basis_set%set_radius)
    NULLIFY (gto_basis_set%sgf_symbol)
    NULLIFY (gto_basis_set%zet)
    gto_basis_set%short_kind_radius = 0.0_dp

  END SUBROUTINE allocate_gto_basis_set

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
! *****************************************************************************
  SUBROUTINE deallocate_gto_basis_set(gto_basis_set)

    ! Deallocate a Gaussian-type orbital (GTO) basis set data set.

    ! - Creation (03.11.2000,MK)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_gto_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(gto_basis_set))  THEN
       IF (ASSOCIATED(gto_basis_set%cgf_symbol)) THEN
          DEALLOCATE (gto_basis_set%cgf_symbol,STAT=stat)
          CPASSERT(stat==0)
       ENDIF
       IF (ASSOCIATED(gto_basis_set%sgf_symbol)) THEN
          DEALLOCATE (gto_basis_set%sgf_symbol,STAT=stat)
          CPASSERT(stat==0)
       ENDIF
       DEALLOCATE (gto_basis_set%norm_cgf,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%set_radius,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%lmax,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%lmin,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%lx,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%ly,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%lz,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%m,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%ncgf_set,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%npgf,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%nsgf_set,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%nshell,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%cphi,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%pgf_radius,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%sphi,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%scon,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%zet,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%first_cgf,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%first_sgf,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%l,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%last_cgf,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%last_sgf,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%n,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set%gcc,STAT=stat)
       CPASSERT(stat==0)
       DEALLOCATE (gto_basis_set,STAT=stat)
       CPASSERT(stat==0)
    END IF
  END SUBROUTINE deallocate_gto_basis_set

! *****************************************************************************
!> \brief ...
!> \param basis_set_in ...
!> \param basis_set_out ...
! *****************************************************************************
  SUBROUTINE copy_gto_basis_set(basis_set_in, basis_set_out)

    ! Copy a Gaussian-type orbital (GTO) basis set data set.

    TYPE(gto_basis_set_type), POINTER        :: basis_set_in, basis_set_out

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_gto_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: maxco, maxpgf, maxshell, &
                                                ncgf, nset, nsgf, stat
    LOGICAL                                  :: failure
    TYPE(gto_basis_set_type), POINTER        :: bin, bout

    failure = .FALSE.

    CPASSERT(ASSOCIATED(basis_set_in))

    CALL allocate_gto_basis_set(basis_set_out)

    bin => basis_set_in
    bout => basis_set_out

    bout%name = bin%name
    bout%kind_radius = bin%kind_radius
    bout%norm_type = bin%norm_type
    bout%nset = bin%nset
    bout%ncgf = bin%ncgf
    bout%nsgf = bin%nsgf
    nset = bin%nset
    ncgf = bin%ncgf
    nsgf = bin%nsgf
    ALLOCATE (bout%cgf_symbol(ncgf),STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE (bout%sgf_symbol(nsgf),STAT=stat)
    CPASSERT(stat==0)
    bout%cgf_symbol = bin%cgf_symbol
    bout%sgf_symbol = bin%sgf_symbol
    ALLOCATE (bout%norm_cgf(ncgf),STAT=stat)
    CPASSERT(stat==0)
    bout%norm_cgf = bin%norm_cgf
    ALLOCATE (bout%set_radius(nset),STAT=stat)
    CPASSERT(stat==0)
    bout%set_radius = bin%set_radius
    ALLOCATE (bout%lmax(nset),bout%lmin(nset),bout%npgf(nset),bout%nshell(nset),STAT=stat)
    CPASSERT(stat==0)
    bout%lmax = bin%lmax
    bout%lmin = bin%lmin
    bout%npgf = bin%npgf
    bout%nshell = bin%nshell
    ALLOCATE (bout%lx(ncgf),bout%ly(ncgf),bout%lz(ncgf),bout%m(nsgf),STAT=stat)
    CPASSERT(stat==0)
    bout%lx = bin%lx
    bout%ly = bin%ly
    bout%lz = bin%lz
    bout%m = bin%m
    ALLOCATE (bout%ncgf_set(nset),bout%nsgf_set(nset),STAT=stat)
    CPASSERT(stat==0)
    bout%ncgf_set = bin%ncgf_set
    bout%nsgf_set = bin%nsgf_set
    maxco = SIZE(bin%cphi,1)
    ALLOCATE (bout%cphi(maxco,ncgf),bout%sphi(maxco,nsgf),bout%scon(maxco,nsgf),STAT=stat)
    CPASSERT(stat==0)
    bout%cphi = bin%cphi
    bout%sphi = bin%sphi
    bout%scon = bin%scon
    maxpgf = MAXVAL(bin%npgf)
    ALLOCATE (bout%pgf_radius(maxpgf,nset),bout%zet(maxpgf,nset),STAT=stat)
    CPASSERT(stat==0)
    bout%pgf_radius = bin%pgf_radius
    bout%zet = bin%zet
    maxshell = MAXVAL(bin%nshell)
    ALLOCATE (bout%first_cgf(maxshell,nset),bout%first_sgf(maxshell,nset),STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE (bout%last_cgf(maxshell,nset),bout%last_sgf(maxshell,nset),STAT=stat)
    CPASSERT(stat==0)
    bout%first_cgf = bin%first_cgf
    bout%first_sgf = bin%first_sgf
    bout%last_cgf = bin%last_cgf
    bout%last_sgf = bin%last_sgf
    ALLOCATE (bout%n(maxshell,nset),bout%l(maxshell,nset),STAT=stat)
    CPASSERT(stat==0)
    bout%n = bin%n
    bout%l = bin%l
    ALLOCATE (bout%gcc(maxpgf,maxshell,nset),STAT=stat)
    CPASSERT(stat==0)
    bout%gcc = bin%gcc

  END SUBROUTINE copy_gto_basis_set

! *****************************************************************************
!> \brief ...
!> \param basis_set ...
!> \param basis_set_add ...
! *****************************************************************************
  SUBROUTINE combine_basis_sets(basis_set, basis_set_add)

    ! Combine two Gaussian-type orbital (GTO) basis sets.

    TYPE(gto_basis_set_type), POINTER        :: basis_set, basis_set_add

    CHARACTER(len=*), PARAMETER :: routineN = 'combine_basis_sets', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=12), DIMENSION(:), POINTER :: cgf_symbol
    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: sgf_symbol
    INTEGER :: iset, ishell, lshell, maxco, maxpgf, maxshell, nc, ncgf, &
      ncgfn, ncgfo, ns, nset, nsetn, nseto, nsgf, nsgfn, nsgfo, stat
    LOGICAL                                  :: failure
    TYPE(gto_basis_set_type), POINTER        :: bad, bout

    failure = .FALSE.

    CPASSERT(ASSOCIATED(basis_set))
    CPASSERT(ASSOCIATED(basis_set_add))

    bad => basis_set_add
    bout => basis_set

    bout%name = bout%name // bad%name
    bout%nset = bout%nset + bad%nset
    bout%ncgf = bout%ncgf + bad%ncgf
    bout%nsgf = bout%nsgf + bad%nsgf
    nset = bout%nset
    ncgf = bout%ncgf
    nsgf = bout%nsgf

    nsetn = bad%nset
    nseto = nset - nsetn
    CALL reallocate(bout%set_radius,1,nset) ! to be defined later
    CALL reallocate(bout%lmax,1,nset)
    CALL reallocate(bout%lmin,1,nset)
    CALL reallocate(bout%npgf,1,nset)
    CALL reallocate(bout%nshell,1,nset)
    bout%lmax(nseto+1:nset) = bad%lmax(1:nsetn)
    bout%lmin(nseto+1:nset) = bad%lmin(1:nsetn)
    bout%npgf(nseto+1:nset) = bad%npgf(1:nsetn)
    bout%nshell(nseto+1:nset) = bad%nshell(1:nsetn)
    CALL reallocate(bout%ncgf_set,1,nset)
    CALL reallocate(bout%nsgf_set,1,nset)
    bout%ncgf_set(nseto+1:nset) = bad%ncgf_set(1:nsetn)
    bout%nsgf_set(nseto+1:nset) = bad%nsgf_set(1:nsetn)

    nsgfn = bad%nsgf
    nsgfo = nsgf - nsgfn
    ncgfn = bad%ncgf
    ncgfo = ncgf - ncgfn

    ALLOCATE (cgf_symbol(ncgf),sgf_symbol(nsgf),STAT=stat)
    CPASSERT(stat==0)
    cgf_symbol(1:ncgfo) = bout%cgf_symbol(1:ncgfo)
    cgf_symbol(ncgfo+1:ncgf) = bad%cgf_symbol(1:ncgfn)
    sgf_symbol(1:nsgfo) = bout%sgf_symbol(1:nsgfo)
    sgf_symbol(nsgfo+1:nsgf) = bad%sgf_symbol(1:nsgfn)
    DEALLOCATE (bout%cgf_symbol,bout%sgf_symbol,STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE (bout%cgf_symbol(ncgf),bout%sgf_symbol(nsgf),STAT=stat)
    CPASSERT(stat==0)
    bout%cgf_symbol = cgf_symbol
    bout%sgf_symbol = sgf_symbol
    DEALLOCATE (cgf_symbol,sgf_symbol,STAT=stat)
    CPASSERT(stat==0)

    CALL reallocate(bout%lx,1,ncgf)
    CALL reallocate(bout%ly,1,ncgf)
    CALL reallocate(bout%lz,1,ncgf)
    CALL reallocate(bout%m,1,nsgf)
    bout%lx(ncgfo+1:ncgf) = bad%lx(1:ncgfn)
    bout%ly(ncgfo+1:ncgf) = bad%ly(1:ncgfn)
    bout%lz(ncgfo+1:ncgf) = bad%lz(1:ncgfn)
    bout%m(nsgfo+1:nsgf) = bad%m(1:nsgfn)

    maxpgf = MAXVAL(bout%npgf)
    CALL reallocate(bout%zet,1,maxpgf,1,nset)
    bout%zet(:,nseto+1:nset) = bad%zet(:,1:nsetn)

    maxshell = MAXVAL(bout%nshell)
    CALL reallocate(bout%l,1,maxshell,1,nset)
    CALL reallocate(bout%n,1,maxshell,1,nset)
    bout%l(:,nseto+1:nset) = bad%l(:,1:nsetn)
    bout%n(:,nseto+1:nset) = bad%n(:,1:nsetn)

    CALL reallocate(bout%first_cgf,1,maxshell,1,nset)
    CALL reallocate(bout%first_sgf,1,maxshell,1,nset)
    CALL reallocate(bout%last_cgf,1,maxshell,1,nset)
    CALL reallocate(bout%last_sgf,1,maxshell,1,nset)
    nc = 0
    ns = 0
    DO iset=1,nset
       DO ishell=1,bout%nshell(iset)
          lshell = bout%l(ishell,iset)
          bout%first_cgf(ishell,iset) = nc + 1
          nc = nc + nco(lshell)
          bout%last_cgf(ishell,iset) = nc
          bout%first_sgf(ishell,iset) = ns + 1
          ns = ns + nso(lshell)
          bout%last_sgf(ishell,iset) = ns
       END DO
    END DO

    CALL reallocate(bout%gcc,1,maxpgf,1,maxshell,1,nset)
    bout%gcc(:,:,nseto+1:nset) = bad%gcc(:,:,1:nsetn)

    ! these arrays are determined later using initialization calls
    CALL reallocate(bout%norm_cgf,1,ncgf)
    maxco = MAX(SIZE(bout%cphi,1),SIZE(bad%cphi,1))
    CALL reallocate(bout%cphi,1,maxco,1,ncgf)
    CALL reallocate(bout%sphi,1,maxco,1,nsgf)
    CALL reallocate(bout%scon,1,maxco,1,nsgf)
    CALL reallocate(bout%pgf_radius,1,maxpgf,1,nset)

  END SUBROUTINE combine_basis_sets

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
!> \param name ...
!> \param norm_type ...
!> \param kind_radius ...
!> \param ncgf ...
!> \param nset ...
!> \param nsgf ...
!> \param cgf_symbol ...
!> \param sgf_symbol ...
!> \param norm_cgf ...
!> \param set_radius ...
!> \param lmax ...
!> \param lmin ...
!> \param lx ...
!> \param ly ...
!> \param lz ...
!> \param m ...
!> \param ncgf_set ...
!> \param npgf ...
!> \param nsgf_set ...
!> \param nshell ...
!> \param cphi ...
!> \param pgf_radius ...
!> \param sphi ...
!> \param scon ...
!> \param zet ...
!> \param first_cgf ...
!> \param first_sgf ...
!> \param l ...
!> \param last_cgf ...
!> \param last_sgf ...
!> \param n ...
!> \param gcc ...
!> \param maxco ...
!> \param maxl ...
!> \param maxpgf ...
!> \param maxsgf_set ...
!> \param maxshell ...
!> \param maxso ...
!> \param nco_sum ...
!> \param npgf_sum ...
!> \param nshell_sum ...
!> \param maxder ...
!> \param short_kind_radius ...
! *****************************************************************************
  SUBROUTINE get_gto_basis_set(gto_basis_set,name,norm_type,kind_radius,ncgf,&
       nset,nsgf,cgf_symbol,sgf_symbol,norm_cgf,set_radius,lmax,lmin,lx,ly,lz,&
       m,ncgf_set,npgf,nsgf_set,nshell,cphi,pgf_radius,sphi,scon,zet,first_cgf,first_sgf,l,&
       last_cgf,last_sgf,n,gcc,maxco,maxl,maxpgf,maxsgf_set,maxshell,maxso,nco_sum,&
       npgf_sum,nshell_sum,maxder,short_kind_radius)

    ! Get informations about a Gaussian-type orbital (GTO) basis set.

    ! - Creation (10.01.2002,MK)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    INTEGER, INTENT(OUT), OPTIONAL           :: norm_type
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: kind_radius
    INTEGER, INTENT(OUT), OPTIONAL           :: ncgf, nset, nsgf
    CHARACTER(LEN=12), DIMENSION(:), &
      OPTIONAL, POINTER                      :: cgf_symbol
    CHARACTER(LEN=6), DIMENSION(:), &
      OPTIONAL, POINTER                      :: sgf_symbol
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: norm_cgf, set_radius
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: lmax, lmin, lx, ly, lz, m, &
                                                ncgf_set, npgf, nsgf_set, &
                                                nshell
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: cphi, pgf_radius, sphi, scon, &
                                                zet
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: first_cgf, first_sgf, l, &
                                                last_cgf, last_sgf, n
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: gcc
    INTEGER, INTENT(OUT), OPTIONAL           :: maxco, maxl, maxpgf, &
                                                maxsgf_set, maxshell, maxso, &
                                                nco_sum, npgf_sum, nshell_sum
    INTEGER, INTENT(IN), OPTIONAL            :: maxder
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: short_kind_radius

    CHARACTER(len=*), PARAMETER :: routineN = 'get_gto_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iset, nder

! -------------------------------------------------------------------------

    IF (ASSOCIATED(gto_basis_set)) THEN
       IF (PRESENT(name)) name = gto_basis_set%name
       IF (PRESENT(norm_type)) norm_type = gto_basis_set%norm_type
       IF (PRESENT(kind_radius)) kind_radius = gto_basis_set%kind_radius
       IF (PRESENT(short_kind_radius)) short_kind_radius = gto_basis_set%short_kind_radius
       IF (PRESENT(ncgf)) ncgf = gto_basis_set%ncgf
       IF (PRESENT(nset)) nset = gto_basis_set%nset
       IF (PRESENT(nsgf)) nsgf = gto_basis_set%nsgf
       IF (PRESENT(cgf_symbol)) cgf_symbol => gto_basis_set%cgf_symbol
       IF (PRESENT(sgf_symbol)) sgf_symbol => gto_basis_set%sgf_symbol
       IF (PRESENT(norm_cgf)) norm_cgf => gto_basis_set%norm_cgf
       IF (PRESENT(set_radius)) set_radius => gto_basis_set%set_radius
       IF (PRESENT(lmax)) lmax => gto_basis_set%lmax
       IF (PRESENT(lmin)) lmin => gto_basis_set%lmin
       IF (PRESENT(lx)) lx => gto_basis_set%lx
       IF (PRESENT(ly)) ly => gto_basis_set%ly
       IF (PRESENT(lz)) lz => gto_basis_set%lz
       IF (PRESENT(m)) m => gto_basis_set%m
       IF (PRESENT(ncgf_set)) ncgf_set => gto_basis_set%ncgf_set
       IF (PRESENT(npgf)) npgf => gto_basis_set%npgf
       IF (PRESENT(nsgf_set)) nsgf_set => gto_basis_set%nsgf_set
       IF (PRESENT(nshell)) nshell => gto_basis_set%nshell
       IF (PRESENT(cphi)) cphi => gto_basis_set%cphi
       IF (PRESENT(pgf_radius)) pgf_radius => gto_basis_set%pgf_radius
       IF (PRESENT(sphi)) sphi => gto_basis_set%sphi
       IF (PRESENT(scon)) scon => gto_basis_set%scon
       IF (PRESENT(zet)) zet => gto_basis_set%zet
       IF (PRESENT(first_cgf)) first_cgf => gto_basis_set%first_cgf
       IF (PRESENT(first_sgf)) first_sgf => gto_basis_set%first_sgf
       IF (PRESENT(l)) l => gto_basis_set%l
       IF (PRESENT(last_cgf)) last_cgf => gto_basis_set%last_cgf
       IF (PRESENT(last_sgf)) last_sgf => gto_basis_set%last_sgf
       IF (PRESENT(n)) n => gto_basis_set%n
       IF (PRESENT(gcc)) gcc => gto_basis_set%gcc
       IF (PRESENT(maxco)) THEN
          maxco = 0
          IF (PRESENT(maxder)) THEN
             nder = maxder
          ELSE
             nder = 0
          END IF
          DO iset=1,gto_basis_set%nset
             maxco = MAX(maxco,gto_basis_set%npgf(iset)*&
                  ncoset(gto_basis_set%lmax(iset)+nder))
          END DO
       END IF
       IF (PRESENT(maxl)) THEN
          maxl = -1
          DO iset=1,gto_basis_set%nset
             maxl = MAX(maxl,gto_basis_set%lmax(iset))
          END DO
       END IF
       IF (PRESENT(maxpgf)) THEN
          maxpgf = 0
          DO iset=1,gto_basis_set%nset
             maxpgf = MAX(maxpgf,gto_basis_set%npgf(iset))
          END DO
       END IF
       IF (PRESENT(maxsgf_set)) THEN
          maxsgf_set = 0
          DO iset=1,gto_basis_set%nset
             maxsgf_set = MAX(maxsgf_set,gto_basis_set%nsgf_set(iset))
          END DO
       END IF
       IF (PRESENT(maxshell)) THEN ! MAXVAL on structure component avoided
          maxshell = 0
          DO iset=1,gto_basis_set%nset
             maxshell = MAX(maxshell,gto_basis_set%nshell(iset))
          END DO
       END IF
       IF (PRESENT(maxso)) THEN
          maxso = 0
          DO iset=1,gto_basis_set%nset
             maxso = MAX(maxso,gto_basis_set%npgf(iset)*&
                  nsoset(gto_basis_set%lmax(iset)))
          END DO
       END IF

       IF (PRESENT(nco_sum)) THEN
          nco_sum = 0
          DO iset=1,gto_basis_set%nset
             nco_sum = nco_sum + gto_basis_set%npgf(iset)*&
                  ncoset(gto_basis_set%lmax(iset))
          END DO
       END IF
       IF (PRESENT(npgf_sum)) npgf_sum = SUM(gto_basis_set%npgf)
       IF (PRESENT(nshell_sum)) nshell_sum = SUM(gto_basis_set%nshell)
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer gto_basis_set is not associated")
    END IF

  END SUBROUTINE get_gto_basis_set

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
! *****************************************************************************
  SUBROUTINE init_aux_basis_set(gto_basis_set)

    ! Initialise a Gaussian-type orbital (GTO) basis set data set.

    ! - Creation (06.12.2000,MK)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(len=*), PARAMETER :: routineN = 'init_aux_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

! -------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(gto_basis_set)) RETURN

    CALL timeset(routineN,handle)

    SELECT CASE (gto_basis_set%norm_type)
    CASE ( 0 )
       ! No normalisation requested
    CASE ( 1 )
       CALL init_norm_cgf_aux_2(gto_basis_set)
    CASE ( 2 )
       ! WARNING this was never tested
       CALL init_norm_cgf_aux(gto_basis_set)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
            "Normalization method not specified")
    END SELECT

    ! Initialise the transformation matrices "pgf" -> "cgf"
    CALL init_cphi_and_sphi(gto_basis_set)
  
    CALL timestop(handle)

  END SUBROUTINE init_aux_basis_set

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
! *****************************************************************************
  SUBROUTINE init_cphi_and_sphi(gto_basis_set)

    ! Initialise the matrices for the transformation of primitive Cartesian
    ! Gaussian-type functions to contracted Cartesian (cphi) and spherical
    ! (sphi) Gaussian-type functions.

    ! - Creation (20.09.2000,MK)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(len=*), PARAMETER :: routineN = 'init_cphi_and_sphi', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_cgf, first_sgf, handle, icgf, ico, ipgf, iset, ishell, &
      l, last_sgf, lmax, lmin, n, n1, n2, ncgf, nn, nn1, nn2, npgf, nsgf

! -------------------------------------------------------------------------
! Build the Cartesian transformation matrix "cphi"

    CALL timeset(routineN,handle)

    DO iset=1,gto_basis_set%nset
       n = ncoset(gto_basis_set%lmax(iset))
       DO ishell=1,gto_basis_set%nshell(iset)
          DO icgf=gto_basis_set%first_cgf(ishell,iset),&
               gto_basis_set%last_cgf(ishell,iset)
             ico = coset(gto_basis_set%lx(icgf),&
                  gto_basis_set%ly(icgf),&
                  gto_basis_set%lz(icgf))
             DO ipgf=1,gto_basis_set%npgf(iset)
                gto_basis_set%cphi(ico,icgf) = gto_basis_set%norm_cgf(icgf)*&
                     gto_basis_set%gcc(ipgf,ishell,iset)
                ico = ico + n
             END DO
          END DO
       END DO
    END DO

    ! Build the spherical transformation matrix "sphi"

    n = SIZE(gto_basis_set%cphi,1)

    IF(n.GT.0) THEN
       DO iset=1,gto_basis_set%nset
          DO ishell=1,gto_basis_set%nshell(iset)
             l = gto_basis_set%l(ishell,iset)
             first_cgf = gto_basis_set%first_cgf(ishell,iset)
             first_sgf = gto_basis_set%first_sgf(ishell,iset)
             ncgf = nco(l)
             nsgf = nso(l)
             CALL dgemm("N","T",n,nsgf,ncgf,&
                  1.0_dp,gto_basis_set%cphi(1,first_cgf),n,&
                  orbtramat(l)%c2s(1,1),nsgf,&
                  0.0_dp,gto_basis_set%sphi(1,first_sgf),n)
          END DO
       END DO
    ENDIF

    ! Build the reduced transformation matrix "scon"
    ! This matrix transforms from cartesian primitifs to spherical contracted functions
    ! "scon" only includes primitifs (lmin -> lmax), whereas sphi is (0 -> lmax)

    gto_basis_set%scon = 0.0_dp
    IF(n.GT.0) THEN
       DO iset=1,gto_basis_set%nset
          lmin = gto_basis_set%lmin(iset)
          lmax = gto_basis_set%lmax(iset)
          npgf = gto_basis_set%npgf(iset)
          nn = ncoset(lmax)-ncoset(lmin-1)
          DO ishell=1,gto_basis_set%nshell(iset)
             first_sgf = gto_basis_set%first_sgf(ishell,iset)
             last_sgf = gto_basis_set%last_sgf(ishell,iset)
             DO ipgf=1,npgf
                nn1 = (ipgf-1) * ncoset(lmax) + ncoset(lmin-1) + 1
                nn2 = ipgf * ncoset(lmax)
                n1 = (ipgf-1) * nn + 1
                n2 = ipgf * nn
                gto_basis_set%scon(n1:n2,first_sgf:last_sgf) = gto_basis_set%sphi(nn1:nn2,first_sgf:last_sgf)
             END DO
          END DO
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE init_cphi_and_sphi

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
! *****************************************************************************
  SUBROUTINE init_norm_cgf_aux(gto_basis_set)

    ! Initialise the normalization factors of the contracted Cartesian Gaussian
    ! functions, if the Gaussian functions represent charge distributions.

    ! - Creation (07.12.2000,MK)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(len=*), PARAMETER :: routineN = 'init_norm_cgf_aux', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icgf, ico, ipgf, iset, &
                                                ishell, jco, jpgf, ll, lmax, &
                                                lmin, lx, ly, lz, n, npgfa, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: fnorm, gcca, gccb
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ff
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: gaa
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: vv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rpgfa, zeta

! -------------------------------------------------------------------------

    failure = .FALSE.
    n = 0
    ll = 0
    DO iset=1,gto_basis_set%nset
       n = MAX(n,gto_basis_set%npgf(iset)*ncoset(gto_basis_set%lmax(iset)))
       ll = MAX(ll,gto_basis_set%lmax(iset))
    END DO

    ALLOCATE (gaa(n,n),STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE (vv(ncoset(ll),ncoset(ll),ll+ll+1),STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE (ff(0:ll+ll),STAT=stat)
    CPASSERT(stat==0)

    DO iset=1,gto_basis_set%nset
       lmax = gto_basis_set%lmax(iset)
       lmin = gto_basis_set%lmin(iset)
       n = ncoset(lmax)
       npgfa = gto_basis_set%npgf(iset)
       rpgfa => gto_basis_set%pgf_radius(1:npgfa,iset)
       zeta => gto_basis_set%zet(1:npgfa,iset)
       CALL coulomb2(lmax,npgfa,zeta,rpgfa,lmin,&
            lmax,npgfa,zeta,rpgfa,lmin,&
            (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,gaa,vv,ff(0:))
       DO ishell=1,gto_basis_set%nshell(iset)
          DO icgf=gto_basis_set%first_cgf(ishell,iset),&
               gto_basis_set%last_cgf(ishell,iset)
             lx = gto_basis_set%lx(icgf)
             ly = gto_basis_set%ly(icgf)
             lz = gto_basis_set%lz(icgf)
             ico = coset(lx,ly,lz)
             fnorm = 0.0_dp
             DO ipgf=1,npgfa
                gcca = gto_basis_set%gcc(ipgf,ishell,iset)
                jco = coset(lx,ly,lz)
                DO jpgf=1,npgfa
                   gccb = gto_basis_set%gcc(jpgf,ishell,iset)
                   fnorm = fnorm + gcca*gccb*gaa(ico,jco)
                   jco = jco + n
                END DO
                ico = ico + n
             END DO
             gto_basis_set%norm_cgf(icgf) = 1.0_dp/SQRT(fnorm)
          END DO
       END DO
    END DO

    DEALLOCATE (vv,ff, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE (gaa,STAT=stat)
    CPASSERT(stat==0)

  END SUBROUTINE init_norm_cgf_aux

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
! *****************************************************************************
  SUBROUTINE init_norm_cgf_aux_2(gto_basis_set)

    ! Initialise the normalization factors of the auxiliary Cartesian Gaussian
    ! functions (Kim-Gordon polarization basis) Norm = 1.

    ! - Creation (07.12.2000,GT)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    INTEGER                                  :: icgf, iset, ishell

! -------------------------------------------------------------------------

    DO iset=1,gto_basis_set%nset
       DO ishell=1,gto_basis_set%nshell(iset)
          DO icgf=gto_basis_set%first_cgf(ishell,iset),&
               gto_basis_set%last_cgf(ishell,iset)
             gto_basis_set%norm_cgf(icgf) = 1.0_dp
          END DO
       END DO
    END DO

  END SUBROUTINE init_norm_cgf_aux_2

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
! *****************************************************************************
  SUBROUTINE init_norm_cgf_orb(gto_basis_set)

    ! Initialise the normalization factors of the contracted Cartesian Gaussian
    ! functions.

    ! - Creation (14.04.2000,MK)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    INTEGER                                  :: icgf, ipgf, iset, ishell, &
                                                jpgf, l, lx, ly, lz
    REAL(KIND=dp)                            :: expzet, fnorm, gcca, gccb, &
                                                prefac, zeta, zetb

! -------------------------------------------------------------------------

    DO iset=1,gto_basis_set%nset
       DO ishell=1,gto_basis_set%nshell(iset)

          l = gto_basis_set%l(ishell,iset)

          expzet = 0.5_dp*REAL(2*l + 3,dp)

          fnorm = 0.0_dp

          DO ipgf=1,gto_basis_set%npgf(iset)
             gcca = gto_basis_set%gcc(ipgf,ishell,iset)
             zeta = gto_basis_set%zet(ipgf,iset)
             DO jpgf=1,gto_basis_set%npgf(iset)
                gccb = gto_basis_set%gcc(jpgf,ishell,iset)
                zetb = gto_basis_set%zet(jpgf,iset)
                fnorm = fnorm + gcca*gccb/(zeta + zetb)**expzet
             END DO
          END DO

          fnorm = 0.5_dp**l*pi**1.5_dp*fnorm

          DO icgf=gto_basis_set%first_cgf(ishell,iset),&
               gto_basis_set%last_cgf(ishell,iset)
             lx = gto_basis_set%lx(icgf)
             ly = gto_basis_set%ly(icgf)
             lz = gto_basis_set%lz(icgf)
             prefac = dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1)
             gto_basis_set%norm_cgf(icgf) = 1.0_dp/SQRT(prefac*fnorm)
          END DO

       END DO
    END DO

  END SUBROUTINE init_norm_cgf_orb

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
! *****************************************************************************
  SUBROUTINE init_norm_cgf_orb_den(gto_basis_set)

    ! Initialise the normalization factors of the contracted Cartesian Gaussian
    ! functions used for frozen density representation.

    ! - Creation (21.09.2002,GT)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    INTEGER                                  :: icgf, ipgf, iset, ishell, l
    REAL(KIND=dp)                            :: expzet, gcca, prefac, zeta

! -------------------------------------------------------------------------

    DO iset=1,gto_basis_set%nset
       DO ishell=1,gto_basis_set%nshell(iset)
          l = gto_basis_set%l(ishell,iset)
          expzet = 0.5_dp*REAL(2*l + 3,dp)
          prefac = (1.0_dp/pi)**1.5_dp
          DO ipgf=1,gto_basis_set%npgf(iset)
             gcca = gto_basis_set%gcc(ipgf,ishell,iset)
             zeta = gto_basis_set%zet(ipgf,iset)
             gto_basis_set%gcc(ipgf,ishell,iset) = prefac*zeta**expzet*gcca
          END DO
          DO icgf=gto_basis_set%first_cgf(ishell,iset),&
               gto_basis_set%last_cgf(ishell,iset)
             gto_basis_set%norm_cgf(icgf) = 1.0_dp
          END DO
       END DO
    END DO

  END SUBROUTINE init_norm_cgf_orb_den

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
! *****************************************************************************
  SUBROUTINE init_orb_basis_set(gto_basis_set)

    ! Initialise a Gaussian-type orbital (GTO) basis set data set.

    ! - Creation (26.10.2000,MK)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(len=*), PARAMETER :: routineN = 'init_orb_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

! -------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(gto_basis_set)) RETURN

    CALL timeset(routineN,handle)

    SELECT CASE (gto_basis_set%norm_type)
    CASE ( 0 )
       ! No normalisation requested
    CASE ( 1 )
       CALL init_norm_cgf_orb_den(gto_basis_set)
    CASE ( 2 )
       ! Normalise the primitive Gaussian functions
       CALL normalise_gcc_orb(gto_basis_set)
       ! Compute the normalization factors of the contracted Gaussian-type
       ! functions
       CALL init_norm_cgf_orb(gto_basis_set)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
            "Normalization method not specified")
    END SELECT

    ! Initialise the transformation matrices "pgf" -> "cgf"

    CALL init_cphi_and_sphi(gto_basis_set)

    CALL timestop(handle)

  END SUBROUTINE init_orb_basis_set

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
! *****************************************************************************
  SUBROUTINE normalise_gcc_orb(gto_basis_set)

    ! Normalise the primitive Cartesian Gaussian functions. The normalization
    ! factor is included in the Gaussian contraction coefficients.

    ! - Creation (20.08.1999,MK)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    INTEGER                                  :: ipgf, iset, ishell, l
    REAL(KIND=dp)                            :: expzet, gcca, prefac, zeta

! -------------------------------------------------------------------------

    DO iset=1,gto_basis_set%nset
       DO ishell=1,gto_basis_set%nshell(iset)
          l = gto_basis_set%l(ishell,iset)
          expzet = 0.25_dp*REAL(2*l + 3,dp)
          prefac = 2.0_dp**l*(2.0_dp/pi)**0.75_dp
          DO ipgf=1,gto_basis_set%npgf(iset)
             gcca = gto_basis_set%gcc(ipgf,ishell,iset)
             zeta = gto_basis_set%zet(ipgf,iset)
             gto_basis_set%gcc(ipgf,ishell,iset) = prefac*zeta**expzet*gcca
          END DO
       END DO
    END DO

  END SUBROUTINE normalise_gcc_orb

! *****************************************************************************
!> \brief ...
!> \param element_symbol ...
!> \param basis_set_name ...
!> \param gto_basis_set ...
!> \param para_env ...
!> \param dft_section ...
! *****************************************************************************
  SUBROUTINE read_gto_basis_set1(element_symbol,basis_set_name,gto_basis_set,&
       para_env,dft_section)

    ! Read a Gaussian-type orbital (GTO) basis set from the database file.

    ! - Creation (13.04.2000,MK)

    CHARACTER(LEN=*), INTENT(IN)             :: element_symbol, basis_set_name
    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: dft_section

    CHARACTER(len=*), PARAMETER :: routineN = 'read_gto_basis_set1', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(len=default_path_length)       :: basis_set_file_name, tmp
    CHARACTER(LEN=default_path_length), &
      DIMENSION(:), POINTER                  :: cbasis
    CHARACTER(LEN=LEN(basis_set_name))       :: bsname
    CHARACTER(LEN=LEN(basis_set_name)+2)     :: bsname2
    CHARACTER(LEN=LEN(element_symbol))       :: symbol
    CHARACTER(LEN=LEN(element_symbol)+2)     :: symbol2
    INTEGER :: i, ibasis, ico, ipgf, irep, iset, ishell, lshell, m, maxco, &
      maxl, maxpgf, maxshell, nbasis, ncgf, nmin, nset, nsgf, stat, strlen1, &
      strlen2
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, n
    LOGICAL                                  :: basis_found, failure, found, &
                                                match
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc
    TYPE(cp_parser_type), POINTER            :: parser

    failure = .FALSE.
    line = ""
    line2 = ""
    symbol = ""
    symbol2 = ""
    bsname = ""
    bsname2 = ""

    nbasis = 1

    gto_basis_set%name = basis_set_name
    CALL section_vals_val_get(dft_section,"BASIS_SET_FILE_NAME",&
         n_rep_val=nbasis)
    ALLOCATE(cbasis(nbasis))
    DO ibasis = 1,nbasis
       CALL section_vals_val_get(dft_section,"BASIS_SET_FILE_NAME",&
            i_rep_val=ibasis,c_val=cbasis(ibasis))
       NULLIFY(parser)
       basis_set_file_name = cbasis(ibasis)
       tmp=basis_set_file_name
       CALL uppercase(tmp)
       IF (INDEX(tmp,"MOLOPT").NE.0) CALL cite_reference(VandeVondele2007)
    END DO

    ! Search for the requested basis set in the basis set file
    ! until the basis set is found or the end of file is reached

    basis_found = .FALSE.
    basis_loop:DO ibasis = 1,nbasis
      IF( basis_found ) EXIT basis_loop
      NULLIFY(parser)
      basis_set_file_name = cbasis(ibasis)
      CALL parser_create(parser,basis_set_file_name,para_env=para_env)

      bsname = basis_set_name
      symbol = element_symbol
      irep   = 0

      tmp=basis_set_name
      CALL uppercase(tmp)
      IF (INDEX(tmp,"MOLOPT").NE.0) CALL cite_reference(VandeVondele2007)

      nset = 0
      maxshell=0
      maxpgf=0
      maxco=0
      ncgf=0
      nsgf=0
      gto_basis_set%nset = nset
      gto_basis_set%ncgf = ncgf
      gto_basis_set%nsgf = nsgf
      CALL reallocate(gto_basis_set%lmax,1,nset)
      CALL reallocate(gto_basis_set%lmin,1,nset)
      CALL reallocate(gto_basis_set%npgf,1,nset)
      CALL reallocate(gto_basis_set%nshell,1,nset)
      CALL reallocate(gto_basis_set%n,1,maxshell,1,nset)
      CALL reallocate(gto_basis_set%l,1,maxshell,1,nset)
      CALL reallocate(gto_basis_set%zet,1,maxpgf,1,nset)
      CALL reallocate(gto_basis_set%gcc,1,maxpgf,1,maxshell,1,nset)
      CALL reallocate(gto_basis_set%set_radius,1,nset)
      CALL reallocate(gto_basis_set%pgf_radius,1,maxpgf,1,nset)
      CALL reallocate(gto_basis_set%first_cgf,1,maxshell,1,nset)
      CALL reallocate(gto_basis_set%first_sgf,1,maxshell,1,nset)
      CALL reallocate(gto_basis_set%last_cgf,1,maxshell,1,nset)
      CALL reallocate(gto_basis_set%last_sgf,1,maxshell,1,nset)
      CALL reallocate(gto_basis_set%ncgf_set,1,nset)
      CALL reallocate(gto_basis_set%nsgf_set,1,nset)
      CALL reallocate(gto_basis_set%cphi,1,maxco,1,ncgf)
      CALL reallocate(gto_basis_set%sphi,1,maxco,1,nsgf)
      CALL reallocate(gto_basis_set%scon,1,maxco,1,nsgf)
      CALL reallocate(gto_basis_set%lx,1,ncgf)
      CALL reallocate(gto_basis_set%ly,1,ncgf)
      CALL reallocate(gto_basis_set%lz,1,ncgf)
      CALL reallocate(gto_basis_set%m,1,nsgf)
      CALL reallocate(gto_basis_set%norm_cgf,1,ncgf)

      IF (tmp.NE."NONE") THEN
        search_loop: DO

           CALL parser_search_string(parser,TRIM(bsname),.TRUE.,found,line)
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
                 NULLIFY (gcc,l,lmax,lmin,n,npgf,nshell,zet)
                 ! Read the basis set information
                 CALL parser_get_object(parser,nset,newline=.TRUE.)

                 CALL reallocate(npgf,1,nset)
                 CALL reallocate(nshell,1,nset)
                 CALL reallocate(lmax,1,nset)
                 CALL reallocate(lmin,1,nset)
                 CALL reallocate(n,1,1,1,nset)

                 maxl = 0
                 maxpgf = 0
                 maxshell = 0

                 DO iset=1,nset
                    CALL parser_get_object(parser,n(1,iset),newline=.TRUE.)
                    CALL parser_get_object(parser,lmin(iset))
                    CALL parser_get_object(parser,lmax(iset))
                    CALL parser_get_object(parser,npgf(iset))
                    maxl = MAX(maxl,lmax(iset))
                    IF (npgf(iset) > maxpgf) THEN
                       maxpgf = npgf(iset)
                       CALL reallocate(zet,1,maxpgf,1,nset)
                       CALL reallocate(gcc,1,maxpgf,1,maxshell,1,nset)
                    END IF
                    nshell(iset) = 0
                    DO lshell=lmin(iset),lmax(iset)
                       nmin = n(1,iset) + lshell - lmin(iset)
                       CALL parser_get_object(parser,ishell)
                       nshell(iset) = nshell(iset) + ishell
                       IF (nshell(iset) > maxshell) THEN
                          maxshell = nshell(iset)
                          CALL reallocate(n,1,maxshell,1,nset)
                          CALL reallocate(l,1,maxshell,1,nset)
                          CALL reallocate(gcc,1,maxpgf,1,maxshell,1,nset)
                       END IF
                       DO i=1,ishell
                          n(nshell(iset)-ishell+i,iset) = nmin + i - 1
                          l(nshell(iset)-ishell+i,iset) = lshell
                       END DO
                    END DO
                    DO ipgf=1,npgf(iset)
                       CALL parser_get_object(parser,zet(ipgf,iset),newline=.TRUE.)
                       DO ishell=1,nshell(iset)
                          CALL parser_get_object(parser,gcc(ipgf,ishell,iset))
                       END DO
                    END DO
                 END DO

                 ! Maximum angular momentum quantum number of the atomic kind

                 CALL init_orbital_pointers(maxl)

                 ! Allocate the global variables

                 gto_basis_set%nset = nset
                 CALL reallocate(gto_basis_set%lmax,1,nset)
                 CALL reallocate(gto_basis_set%lmin,1,nset)
                 CALL reallocate(gto_basis_set%npgf,1,nset)
                 CALL reallocate(gto_basis_set%nshell,1,nset)
                 CALL reallocate(gto_basis_set%n,1,maxshell,1,nset)
                 CALL reallocate(gto_basis_set%l,1,maxshell,1,nset)
                 CALL reallocate(gto_basis_set%zet,1,maxpgf,1,nset)
                 CALL reallocate(gto_basis_set%gcc,1,maxpgf,1,maxshell,1,nset)

                 ! Copy the basis set information into the data structure

                 DO iset=1,nset
                    gto_basis_set%lmax(iset) = lmax(iset)
                    gto_basis_set%lmin(iset) = lmin(iset)
                    gto_basis_set%npgf(iset) = npgf(iset)
                    gto_basis_set%nshell(iset) = nshell(iset)
                    DO ishell=1,nshell(iset)
                       gto_basis_set%n(ishell,iset) = n(ishell,iset)
                       gto_basis_set%l(ishell,iset) = l(ishell,iset)
                       DO ipgf=1,npgf(iset)
                          gto_basis_set%gcc(ipgf,ishell,iset) = gcc(ipgf,ishell,iset)
                       END DO
                    END DO
                    DO ipgf=1,npgf(iset)
                       gto_basis_set%zet(ipgf,iset) = zet(ipgf,iset)
                    END DO
                 END DO

                 ! Initialise the depending atomic kind information

                 CALL reallocate(gto_basis_set%set_radius,1,nset)
                 CALL reallocate(gto_basis_set%pgf_radius,1,maxpgf,1,nset)
                 CALL reallocate(gto_basis_set%first_cgf,1,maxshell,1,nset)
                 CALL reallocate(gto_basis_set%first_sgf,1,maxshell,1,nset)
                 CALL reallocate(gto_basis_set%last_cgf,1,maxshell,1,nset)
                 CALL reallocate(gto_basis_set%last_sgf,1,maxshell,1,nset)
                 CALL reallocate(gto_basis_set%ncgf_set,1,nset)
                 CALL reallocate(gto_basis_set%nsgf_set,1,nset)

                 maxco = 0
                 ncgf = 0
                 nsgf = 0

                 DO iset=1,nset
                    gto_basis_set%ncgf_set(iset) = 0
                    gto_basis_set%nsgf_set(iset) = 0
                    DO ishell=1,nshell(iset)
                       lshell = gto_basis_set%l(ishell,iset)
                       gto_basis_set%first_cgf(ishell,iset) = ncgf + 1
                       ncgf = ncgf + nco(lshell)
                       gto_basis_set%last_cgf(ishell,iset) = ncgf
                       gto_basis_set%ncgf_set(iset) =&
                            gto_basis_set%ncgf_set(iset) + nco(lshell)
                       gto_basis_set%first_sgf(ishell,iset) = nsgf + 1
                       nsgf = nsgf + nso(lshell)
                       gto_basis_set%last_sgf(ishell,iset) = nsgf
                       gto_basis_set%nsgf_set(iset) =&
                            gto_basis_set%nsgf_set(iset) + nso(lshell)
                    END DO
                    maxco = MAX(maxco,npgf(iset)*ncoset(lmax(iset)))
                 END DO

                 gto_basis_set%ncgf = ncgf
                 gto_basis_set%nsgf = nsgf

                 CALL reallocate(gto_basis_set%cphi,1,maxco,1,ncgf)
                 CALL reallocate(gto_basis_set%sphi,1,maxco,1,nsgf)
                 CALL reallocate(gto_basis_set%scon,1,maxco,1,nsgf)
                 CALL reallocate(gto_basis_set%lx,1,ncgf)
                 CALL reallocate(gto_basis_set%ly,1,ncgf)
                 CALL reallocate(gto_basis_set%lz,1,ncgf)
                 CALL reallocate(gto_basis_set%m,1,nsgf)
                 CALL reallocate(gto_basis_set%norm_cgf,1,ncgf)
                 ALLOCATE (gto_basis_set%cgf_symbol(ncgf),STAT=stat)
                 CPASSERT(stat==0)

                 ALLOCATE (gto_basis_set%sgf_symbol(nsgf),STAT=stat)
                 CPASSERT(stat==0)

                 ncgf = 0
                 nsgf = 0

                 DO iset=1,nset
                    DO ishell=1,nshell(iset)
                       lshell = gto_basis_set%l(ishell,iset)
                       DO ico=ncoset(lshell-1)+1,ncoset(lshell)
                          ncgf = ncgf + 1
                          gto_basis_set%lx(ncgf) = indco(1,ico)
                          gto_basis_set%ly(ncgf) = indco(2,ico)
                          gto_basis_set%lz(ncgf) = indco(3,ico)
                          gto_basis_set%cgf_symbol(ncgf) =&
                               cgf_symbol(n(ishell,iset),(/gto_basis_set%lx(ncgf),&
                               gto_basis_set%ly(ncgf),&
                               gto_basis_set%lz(ncgf)/))
                       END DO
                       DO m=-lshell,lshell
                          nsgf = nsgf + 1
                          gto_basis_set%m(nsgf) = m
                          gto_basis_set%sgf_symbol(nsgf) =&
                               sgf_symbol(n(ishell,iset),lshell,m)
                       END DO
                    END DO
                 END DO

                 DEALLOCATE (gcc,l,lmax,lmin,n,npgf,nshell,zet,STAT=stat)
                 CPASSERT(stat==0)

                 basis_found = .TRUE.
                 EXIT search_loop
              END IF
           ELSE
             EXIT search_loop
           END IF
        END DO search_loop
      ELSE
        match=.FALSE.
        ALLOCATE (gto_basis_set%cgf_symbol(ncgf),STAT=stat)
        CPASSERT(stat==0)
        ALLOCATE (gto_basis_set%sgf_symbol(nsgf),STAT=stat)
        CPASSERT(stat==0)
      ENDIF

      CALL parser_release(parser)

    END DO basis_loop

    IF (tmp .NE. "NONE") THEN
      IF( .NOT. basis_found ) THEN
        basis_set_file_name = ""
        DO ibasis = 1,nbasis
          basis_set_file_name = TRIM(basis_set_file_name)//"<"//TRIM(cbasis(ibasis))//"> "
        END DO
        CALL stop_program(routineN,moduleN,__LINE__,&
             "The requested basis set <"//TRIM(bsname)//&
             "> for element <"//TRIM(symbol)//"> was not "//&
             "found in the basis set files "//&
             TRIM(basis_set_file_name),para_env)
      END IF
    END IF
    DEALLOCATE(cbasis)
     
  END SUBROUTINE read_gto_basis_set1

! *****************************************************************************
!> \brief ...
!> \param element_symbol ...
!> \param basis_type ...
!> \param gto_basis_set ...
!> \param basis_section ...
!> \param irep ...
! *****************************************************************************
  SUBROUTINE read_gto_basis_set2(element_symbol,basis_type,gto_basis_set,&
       basis_section,irep)

    ! Read a Gaussian-type orbital (GTO) basis set from the database file.

    ! - Creation (13.04.2000,MK)

    CHARACTER(LEN=*), INTENT(IN)             :: element_symbol
    CHARACTER(LEN=*), INTENT(INOUT)          :: basis_type
    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: basis_section
    INTEGER                                  :: irep

    CHARACTER(len=*), PARAMETER :: routineN = 'read_gto_basis_set2', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=20*default_string_length)  :: line_att
    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(LEN=default_path_length)       :: bsname, bsname2
    CHARACTER(LEN=LEN(element_symbol))       :: symbol
    CHARACTER(LEN=LEN(element_symbol)+2)     :: symbol2
    INTEGER                                  :: i, ico, ipgf, iset, ishell, &
                                                lshell, m, maxco, maxl, &
                                                maxpgf, maxshell, ncgf, nmin, &
                                                nset, nsgf, stat
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, n
    LOGICAL                                  :: failure, is_ok
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(val_type), POINTER                  :: val

    failure = .FALSE.
    line = ""
    line2 = ""
    symbol = ""
    symbol2 = ""
    bsname = ""
    bsname2 = ""

    gto_basis_set%name = " "

    bsname = " "
    symbol = element_symbol

    nset = 0
    maxshell=0
    maxpgf=0
    maxco=0
    ncgf=0
    nsgf=0
    gto_basis_set%nset = nset
    gto_basis_set%ncgf = ncgf
    gto_basis_set%nsgf = nsgf
    CALL reallocate(gto_basis_set%lmax,1,nset)
    CALL reallocate(gto_basis_set%lmin,1,nset)
    CALL reallocate(gto_basis_set%npgf,1,nset)
    CALL reallocate(gto_basis_set%nshell,1,nset)
    CALL reallocate(gto_basis_set%n,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%l,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%zet,1,maxpgf,1,nset)
    CALL reallocate(gto_basis_set%gcc,1,maxpgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%set_radius,1,nset)
    CALL reallocate(gto_basis_set%pgf_radius,1,maxpgf,1,nset)
    CALL reallocate(gto_basis_set%first_cgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%first_sgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%last_cgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%last_sgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%ncgf_set,1,nset)
    CALL reallocate(gto_basis_set%nsgf_set,1,nset)
    CALL reallocate(gto_basis_set%cphi,1,maxco,1,ncgf)
    CALL reallocate(gto_basis_set%sphi,1,maxco,1,nsgf)
    CALL reallocate(gto_basis_set%scon,1,maxco,1,nsgf)
    CALL reallocate(gto_basis_set%lx,1,ncgf)
    CALL reallocate(gto_basis_set%ly,1,ncgf)
    CALL reallocate(gto_basis_set%lz,1,ncgf)
    CALL reallocate(gto_basis_set%m,1,nsgf)
    CALL reallocate(gto_basis_set%norm_cgf,1,ncgf)

    basis_type=""
    CALL section_vals_val_get(basis_section,"_SECTION_PARAMETERS_",i_rep_val=irep,c_val=basis_type)
    IF(basis_type=="Orbital") basis_type="ORB"

    NULLIFY(list,val)
    CALL section_vals_list_get(basis_section,"_DEFAULT_KEYWORD_",i_rep_section=irep,list=list)
    CALL uppercase(symbol)
    CALL uppercase(bsname)

    NULLIFY (gcc,l,lmax,lmin,n,npgf,nshell,zet)
    ! Read the basis set information
    is_ok=cp_sll_val_next(list,val)
    IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
        __LINE__,&
        "Error reading the Basis set from input file!!")
    CALL val_get(val,c_val=line_att)
    READ(line_att,*)nset

    CALL reallocate(npgf,1,nset)
    CALL reallocate(nshell,1,nset)
    CALL reallocate(lmax,1,nset)
    CALL reallocate(lmin,1,nset)
    CALL reallocate(n,1,1,1,nset)

    maxl = 0
    maxpgf = 0
    maxshell = 0

    DO iset=1,nset
       is_ok=cp_sll_val_next(list,val)
       IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
          __LINE__,&
          "Error reading the Basis set from input file!!")
       CALL val_get(val,c_val=line_att)
       READ(line_att,*)n(1,iset)
       CALL remove_word(line_att)
       READ(line_att,*)lmin(iset)
       CALL remove_word(line_att)
       READ(line_att,*)lmax(iset)
       CALL remove_word(line_att)
       READ(line_att,*)npgf(iset)
       CALL remove_word(line_att)
       maxl = MAX(maxl,lmax(iset))
       IF (npgf(iset) > maxpgf) THEN
          maxpgf = npgf(iset)
          CALL reallocate(zet,1,maxpgf,1,nset)
          CALL reallocate(gcc,1,maxpgf,1,maxshell,1,nset)
       END IF
       nshell(iset) = 0
       DO lshell=lmin(iset),lmax(iset)
          nmin = n(1,iset) + lshell - lmin(iset)
          READ(line_att,*)ishell
          CALL remove_word(line_att)
          nshell(iset) = nshell(iset) + ishell
          IF (nshell(iset) > maxshell) THEN
             maxshell = nshell(iset)
             CALL reallocate(n,1,maxshell,1,nset)
             CALL reallocate(l,1,maxshell,1,nset)
             CALL reallocate(gcc,1,maxpgf,1,maxshell,1,nset)
          END IF
          DO i=1,ishell
             n(nshell(iset)-ishell+i,iset) = nmin + i - 1
             l(nshell(iset)-ishell+i,iset) = lshell
          END DO
       END DO
       IF (LEN_TRIM(line_att)/=0) CALL stop_program(routineN,moduleN,&
           __LINE__,&
           "Error reading the Basis from input file!!")
       DO ipgf=1,npgf(iset)
          is_ok=cp_sll_val_next(list,val)
          IF (.NOT.is_ok) CALL stop_program(routineN,moduleN,&
             __LINE__,&
             "Error reading the Basis set from input file!!")
          CALL val_get(val,c_val=line_att)
          READ(line_att,*)zet(ipgf,iset),(gcc(ipgf,ishell,iset),ishell=1,nshell(iset))
       END DO
    END DO

    ! Maximum angular momentum quantum number of the atomic kind

    CALL init_orbital_pointers(maxl)

    ! Allocate the global variables

    gto_basis_set%nset = nset
    CALL reallocate(gto_basis_set%lmax,1,nset)
    CALL reallocate(gto_basis_set%lmin,1,nset)
    CALL reallocate(gto_basis_set%npgf,1,nset)
    CALL reallocate(gto_basis_set%nshell,1,nset)
    CALL reallocate(gto_basis_set%n,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%l,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%zet,1,maxpgf,1,nset)
    CALL reallocate(gto_basis_set%gcc,1,maxpgf,1,maxshell,1,nset)

    ! Copy the basis set information into the data structure

    DO iset=1,nset
       gto_basis_set%lmax(iset) = lmax(iset)
       gto_basis_set%lmin(iset) = lmin(iset)
       gto_basis_set%npgf(iset) = npgf(iset)
       gto_basis_set%nshell(iset) = nshell(iset)
       DO ishell=1,nshell(iset)
          gto_basis_set%n(ishell,iset) = n(ishell,iset)
          gto_basis_set%l(ishell,iset) = l(ishell,iset)
          DO ipgf=1,npgf(iset)
             gto_basis_set%gcc(ipgf,ishell,iset) = gcc(ipgf,ishell,iset)
          END DO
       END DO
       DO ipgf=1,npgf(iset)
          gto_basis_set%zet(ipgf,iset) = zet(ipgf,iset)
       END DO
    END DO

    ! Initialise the depending atomic kind information

    CALL reallocate(gto_basis_set%set_radius,1,nset)
    CALL reallocate(gto_basis_set%pgf_radius,1,maxpgf,1,nset)
    CALL reallocate(gto_basis_set%first_cgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%first_sgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%last_cgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%last_sgf,1,maxshell,1,nset)
    CALL reallocate(gto_basis_set%ncgf_set,1,nset)
    CALL reallocate(gto_basis_set%nsgf_set,1,nset)

    maxco = 0
    ncgf = 0
    nsgf = 0

    DO iset=1,nset
       gto_basis_set%ncgf_set(iset) = 0
       gto_basis_set%nsgf_set(iset) = 0
       DO ishell=1,nshell(iset)
          lshell = gto_basis_set%l(ishell,iset)
          gto_basis_set%first_cgf(ishell,iset) = ncgf + 1
          ncgf = ncgf + nco(lshell)
          gto_basis_set%last_cgf(ishell,iset) = ncgf
          gto_basis_set%ncgf_set(iset) =&
               gto_basis_set%ncgf_set(iset) + nco(lshell)
          gto_basis_set%first_sgf(ishell,iset) = nsgf + 1
          nsgf = nsgf + nso(lshell)
          gto_basis_set%last_sgf(ishell,iset) = nsgf
          gto_basis_set%nsgf_set(iset) =&
               gto_basis_set%nsgf_set(iset) + nso(lshell)
       END DO
       maxco = MAX(maxco,npgf(iset)*ncoset(lmax(iset)))
    END DO

    gto_basis_set%ncgf = ncgf
    gto_basis_set%nsgf = nsgf

    CALL reallocate(gto_basis_set%cphi,1,maxco,1,ncgf)
    CALL reallocate(gto_basis_set%sphi,1,maxco,1,nsgf)
    CALL reallocate(gto_basis_set%scon,1,maxco,1,nsgf)
    CALL reallocate(gto_basis_set%lx,1,ncgf)
    CALL reallocate(gto_basis_set%ly,1,ncgf)
    CALL reallocate(gto_basis_set%lz,1,ncgf)
    CALL reallocate(gto_basis_set%m,1,nsgf)
    CALL reallocate(gto_basis_set%norm_cgf,1,ncgf)
    ALLOCATE (gto_basis_set%cgf_symbol(ncgf),STAT=stat)
    CPASSERT(stat==0)

    ALLOCATE (gto_basis_set%sgf_symbol(nsgf),STAT=stat)
    CPASSERT(stat==0)

    ncgf = 0
    nsgf = 0

    DO iset=1,nset
       DO ishell=1,nshell(iset)
          lshell = gto_basis_set%l(ishell,iset)
          DO ico=ncoset(lshell-1)+1,ncoset(lshell)
             ncgf = ncgf + 1
             gto_basis_set%lx(ncgf) = indco(1,ico)
             gto_basis_set%ly(ncgf) = indco(2,ico)
             gto_basis_set%lz(ncgf) = indco(3,ico)
             gto_basis_set%cgf_symbol(ncgf) =&
                  cgf_symbol(n(ishell,iset),(/gto_basis_set%lx(ncgf),&
                  gto_basis_set%ly(ncgf),&
                  gto_basis_set%lz(ncgf)/))
          END DO
          DO m=-lshell,lshell
             nsgf = nsgf + 1
             gto_basis_set%m(nsgf) = m
             gto_basis_set%sgf_symbol(nsgf) =&
                  sgf_symbol(n(ishell,iset),lshell,m)
          END DO
       END DO
    END DO

    DEALLOCATE (gcc,l,lmax,lmin,n,npgf,nshell,zet,STAT=stat)
    CPASSERT(stat==0)

  END SUBROUTINE read_gto_basis_set2

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
!> \param name ...
!> \param norm_type ...
!> \param kind_radius ...
!> \param ncgf ...
!> \param nset ...
!> \param nsgf ...
!> \param cgf_symbol ...
!> \param sgf_symbol ...
!> \param norm_cgf ...
!> \param set_radius ...
!> \param lmax ...
!> \param lmin ...
!> \param lx ...
!> \param ly ...
!> \param lz ...
!> \param m ...
!> \param ncgf_set ...
!> \param npgf ...
!> \param nsgf_set ...
!> \param nshell ...
!> \param cphi ...
!> \param pgf_radius ...
!> \param sphi ...
!> \param scon ...
!> \param zet ...
!> \param first_cgf ...
!> \param first_sgf ...
!> \param l ...
!> \param last_cgf ...
!> \param last_sgf ...
!> \param n ...
!> \param gcc ...
!> \param short_kind_radius ...
! *****************************************************************************
  SUBROUTINE set_gto_basis_set(gto_basis_set,name,norm_type,kind_radius,ncgf,&
       nset,nsgf,cgf_symbol,sgf_symbol,norm_cgf,set_radius,lmax,&
       lmin,lx,ly,lz,m,ncgf_set,npgf,nsgf_set,nshell,&
       cphi,pgf_radius,sphi,scon,zet,first_cgf,first_sgf,l,&
       last_cgf,last_sgf,n,gcc,short_kind_radius)

    ! Set the components of Gaussian-type orbital (GTO) basis set data set.

    ! - Creation (10.01.2002,MK)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name
    INTEGER, INTENT(IN), OPTIONAL            :: norm_type
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: kind_radius
    INTEGER, INTENT(IN), OPTIONAL            :: ncgf, nset, nsgf
    CHARACTER(LEN=12), DIMENSION(:), &
      OPTIONAL, POINTER                      :: cgf_symbol
    CHARACTER(LEN=6), DIMENSION(:), &
      OPTIONAL, POINTER                      :: sgf_symbol
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: norm_cgf, set_radius
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: lmax, lmin, lx, ly, lz, m, &
                                                ncgf_set, npgf, nsgf_set, &
                                                nshell
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: cphi, pgf_radius, sphi, scon, &
                                                zet
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: first_cgf, first_sgf, l, &
                                                last_cgf, last_sgf, n
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: gcc
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: short_kind_radius

    CHARACTER(len=*), PARAMETER :: routineN = 'set_gto_basis_set', &
      routineP = moduleN//':'//routineN

! -------------------------------------------------------------------------

    IF (ASSOCIATED(gto_basis_set)) THEN
       IF (PRESENT(name)) gto_basis_set%name = name
       IF (PRESENT(norm_type)) gto_basis_set%norm_type = norm_type
       IF (PRESENT(kind_radius)) gto_basis_set%kind_radius = kind_radius
       IF (PRESENT(short_kind_radius)) gto_basis_set%short_kind_radius = short_kind_radius
       IF (PRESENT(ncgf)) gto_basis_set%ncgf = ncgf
       IF (PRESENT(nset)) gto_basis_set%nset = nset
       IF (PRESENT(nsgf)) gto_basis_set%nsgf = nsgf
       IF (PRESENT(cgf_symbol)) gto_basis_set%cgf_symbol(:) = cgf_symbol(:)
       IF (PRESENT(sgf_symbol)) gto_basis_set%sgf_symbol(:) = sgf_symbol(:)
       IF (PRESENT(norm_cgf)) gto_basis_set%norm_cgf(:) = norm_cgf(:)
       IF (PRESENT(set_radius)) gto_basis_set%set_radius(:) = set_radius(:)
       IF (PRESENT(lmax)) gto_basis_set%lmax(:) = lmax(:)
       IF (PRESENT(lmin)) gto_basis_set%lmin(:) = lmin(:)
       IF (PRESENT(lx)) gto_basis_set%lx(:) = lx(:)
       IF (PRESENT(ly)) gto_basis_set%ly(:) = ly(:)
       IF (PRESENT(lz)) gto_basis_set%lz(:) = lz(:)
       IF (PRESENT(m)) gto_basis_set%m(:) = m(:)
       IF (PRESENT(ncgf_set)) gto_basis_set%ncgf_set(:) = ncgf_set(:)
       IF (PRESENT(npgf)) gto_basis_set%npgf(:) = npgf(:)
       IF (PRESENT(nsgf_set)) gto_basis_set%nsgf_set(:) = nsgf_set(:)
       IF (PRESENT(nshell)) gto_basis_set%nshell(:) = nshell(:)
       IF (PRESENT(cphi)) gto_basis_set%cphi(:,:) = cphi(:,:)
       IF (PRESENT(pgf_radius)) gto_basis_set%pgf_radius(:,:) = pgf_radius(:,:)
       IF (PRESENT(sphi)) gto_basis_set%sphi(:,:) = sphi(:,:)
       IF (PRESENT(scon)) gto_basis_set%scon(:,:) = scon(:,:)
       IF (PRESENT(zet)) gto_basis_set%zet(:,:) = zet(:,:)
       IF (PRESENT(first_cgf)) gto_basis_set%first_cgf(:,:) = first_cgf(:,:)
       IF (PRESENT(first_sgf)) gto_basis_set%first_sgf(:,:) = first_sgf(:,:)
       IF (PRESENT(l)) l(:,:) = gto_basis_set%l(:,:)
       IF (PRESENT(last_cgf)) gto_basis_set%last_cgf(:,:) = last_cgf(:,:)
       IF (PRESENT(last_sgf)) gto_basis_set%last_sgf(:,:) = last_sgf(:,:)
       IF (PRESENT(n)) gto_basis_set%n(:,:) = n(:,:)
       IF (PRESENT(gcc)) gto_basis_set%gcc(:,:,:) = gcc(:,:,:)
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer gto_basis_set is not associated")
    END IF

  END SUBROUTINE set_gto_basis_set

! *****************************************************************************
!> \brief ...
!> \param gto_basis_set ...
!> \param output_unit ...
!> \param header ...
! *****************************************************************************
  SUBROUTINE write_gto_basis_set(gto_basis_set,output_unit,header)

    ! Write a Gaussian-type orbital (GTO) basis set data set to the output
    ! unit.

    ! - Creation (09.01.2002,MK)

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set
    INTEGER, INTENT(in)                      :: output_unit
    CHARACTER(len=*), OPTIONAL               :: header

    INTEGER                                  :: ipgf, iset, ishell

! -------------------------------------------------------------------------

    IF (ASSOCIATED(gto_basis_set).AND.(output_unit>0)) THEN

       IF(PRESENT(header)) THEN
         WRITE (UNIT=output_unit,FMT="(/,T6,A,T41,A40)")&
              TRIM(header),TRIM(gto_basis_set%name)
       END IF

       WRITE (UNIT=output_unit,FMT="(/,(T8,A,T71,I10))")&
            "Number of orbital shell sets:            ",&
            gto_basis_set%nset,&
            "Number of orbital shells:                ",&
            SUM(gto_basis_set%nshell(:)),&
            "Number of primitive Cartesian functions: ",&
            SUM(gto_basis_set%npgf(:)),&
            "Number of Cartesian basis functions:     ",&
            gto_basis_set%ncgf,&
            "Number of spherical basis functions:     ",&
            gto_basis_set%nsgf,&
            "Norm type:                               ",&
            gto_basis_set%norm_type

       WRITE (UNIT=output_unit,FMT="(/,T6,A,T41,A40,/,/,T25,A)")&
            "GTO basis set information for",TRIM(gto_basis_set%name),&
            "Set   Shell     n   l            Exponent    Coefficient"

       DO iset=1,gto_basis_set%nset
          WRITE (UNIT=output_unit,FMT="(A)") ""
          DO ishell=1,gto_basis_set%nshell(iset)
             WRITE (UNIT=output_unit,&
                  FMT="(T25,I3,4X,I4,4X,I2,2X,I2,(T51,2F15.6))")&
                  iset,ishell,&
                  gto_basis_set%n(ishell,iset),&
                  gto_basis_set%l(ishell,iset),&
                  (gto_basis_set%zet(ipgf,iset),&
                  gto_basis_set%gcc(ipgf,ishell,iset),&
                  ipgf=1,gto_basis_set%npgf(iset))
          END DO
       END DO

    END IF

  END SUBROUTINE write_gto_basis_set

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param orb_basis_set ...
!> \param output_unit ...
!> \param header ...
! *****************************************************************************
  SUBROUTINE write_orb_basis_set(orb_basis_set,output_unit,header)

    ! Write a Gaussian-type orbital (GTO) basis set data set to the output
    ! unit.

    ! - Creation (09.01.2002,MK)

    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    INTEGER, INTENT(in)                      :: output_unit
    CHARACTER(len=*), OPTIONAL               :: header

    INTEGER                                  :: icgf, ico, ipgf, iset, ishell

    IF (ASSOCIATED(orb_basis_set).AND.(output_unit>0)) THEN
       IF(PRESENT(header)) THEN
         WRITE (UNIT=output_unit,FMT="(/,T6,A,T41,A40)")&
              TRIM(header),TRIM(orb_basis_set%name)
       END IF

       WRITE (UNIT=output_unit,FMT="(/,(T8,A,T71,I10))")&
            "Number of orbital shell sets:            ",&
            orb_basis_set%nset,&
            "Number of orbital shells:                ",&
            SUM(orb_basis_set%nshell(:)),&
            "Number of primitive Cartesian functions: ",&
            SUM(orb_basis_set%npgf(:)),&
            "Number of Cartesian basis functions:     ",&
            orb_basis_set%ncgf,&
            "Number of spherical basis functions:     ",&
            orb_basis_set%nsgf,&
            "Norm type:                               ",&
            orb_basis_set%norm_type

       WRITE (UNIT=output_unit,FMT="(/,T8,A,/,/,T25,A)")&
            "Normalised Cartesian orbitals:",&
            "Set   Shell   Orbital            Exponent    Coefficient"

       icgf = 0

       DO iset=1,orb_basis_set%nset
          DO ishell=1,orb_basis_set%nshell(iset)
             WRITE (UNIT=output_unit,FMT="(A)") ""
             DO ico=1,nco(orb_basis_set%l(ishell,iset))
                icgf = icgf + 1
                WRITE (UNIT=output_unit,&
                     FMT="(T25,I3,4X,I4,3X,A12,(T51,2F15.6))")&
                     iset,ishell,orb_basis_set%cgf_symbol(icgf),&
                     (orb_basis_set%zet(ipgf,iset),&
                     orb_basis_set%norm_cgf(icgf)*&
                     orb_basis_set%gcc(ipgf,ishell,iset),&
                     ipgf=1,orb_basis_set%npgf(iset))
             END DO
          END DO
       END DO
    END IF

  END SUBROUTINE write_orb_basis_set

! *****************************************************************************
!> \brief ...
!> \param sto_basis_set ...
! *****************************************************************************
  SUBROUTINE allocate_sto_basis_set(sto_basis_set)

    TYPE(sto_basis_set_type), POINTER        :: sto_basis_set

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_sto_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

! -------------------------------------------------------------------------

    failure = .FALSE.
    CALL deallocate_sto_basis_set(sto_basis_set)

    ALLOCATE (sto_basis_set,STAT=stat)
    CPASSERT(stat==0)

    sto_basis_set%name="NONAME"
    NULLIFY (sto_basis_set%symbol)
    NULLIFY (sto_basis_set%nq)
    NULLIFY (sto_basis_set%lq)
    NULLIFY (sto_basis_set%zet)

  END SUBROUTINE allocate_sto_basis_set

! *****************************************************************************
!> \brief ...
!> \param sto_basis_set ...
! *****************************************************************************
  SUBROUTINE deallocate_sto_basis_set(sto_basis_set)

    TYPE(sto_basis_set_type), POINTER        :: sto_basis_set

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_sto_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

! -------------------------------------------------------------------------

    failure = .FALSE.
    IF (ASSOCIATED(sto_basis_set)) THEN
       IF (ASSOCIATED(sto_basis_set%symbol)) THEN
          DEALLOCATE (sto_basis_set%symbol,STAT=stat)
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED(sto_basis_set%nq)) THEN
          DEALLOCATE (sto_basis_set%nq,STAT=stat)
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED(sto_basis_set%lq)) THEN
          DEALLOCATE (sto_basis_set%lq,STAT=stat)
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED(sto_basis_set%zet)) THEN
          DEALLOCATE (sto_basis_set%zet,STAT=stat)
          CPASSERT(stat==0)
       END IF

       DEALLOCATE (sto_basis_set,STAT=stat)
       CPASSERT(stat==0)
    END IF
  END SUBROUTINE deallocate_sto_basis_set

! *****************************************************************************
!> \brief ...
!> \param sto_basis_set ...
!> \param name ...
!> \param nshell ...
!> \param symbol ...
!> \param nq ...
!> \param lq ...
!> \param zet ...
!> \param maxlq ...
!> \param numsto ...
! *****************************************************************************
  SUBROUTINE get_sto_basis_set(sto_basis_set,name,nshell,symbol,nq,lq,zet,maxlq,numsto)

    TYPE(sto_basis_set_type), POINTER        :: sto_basis_set
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    INTEGER, INTENT(OUT), OPTIONAL           :: nshell
    CHARACTER(LEN=6), DIMENSION(:), &
      OPTIONAL, POINTER                      :: symbol
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: nq, lq
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: zet
    INTEGER, INTENT(OUT), OPTIONAL           :: maxlq, numsto

    CHARACTER(len=*), PARAMETER :: routineN = 'get_sto_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iset

! -------------------------------------------------------------------------

    IF (ASSOCIATED(sto_basis_set)) THEN
       IF (PRESENT(name)) name = sto_basis_set%name
       IF (PRESENT(nshell)) nshell = sto_basis_set%nshell
       IF (PRESENT(symbol)) symbol => sto_basis_set%symbol
       IF (PRESENT(nq)) nq => sto_basis_set%nq
       IF (PRESENT(lq)) lq => sto_basis_set%lq
       IF (PRESENT(zet)) zet => sto_basis_set%zet
       IF (PRESENT(maxlq)) THEN
          maxlq = MAXVAL ( sto_basis_set%lq (1:sto_basis_set%nshell) )
       END IF
       IF (PRESENT(numsto)) THEN
          numsto = 0
          DO iset=1,sto_basis_set%nshell
             numsto = numsto + 2*sto_basis_set%lq(iset)+1
          END DO
       END IF
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer sto_basis_set is not associated")
    END IF

  END SUBROUTINE get_sto_basis_set

! *****************************************************************************
!> \brief ...
!> \param sto_basis_set ...
!> \param name ...
!> \param nshell ...
!> \param symbol ...
!> \param nq ...
!> \param lq ...
!> \param zet ...
! *****************************************************************************
  SUBROUTINE set_sto_basis_set(sto_basis_set,name,nshell,symbol,nq,lq,zet)

    TYPE(sto_basis_set_type), POINTER        :: sto_basis_set
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name
    INTEGER, INTENT(IN), OPTIONAL            :: nshell
    CHARACTER(LEN=6), DIMENSION(:), &
      OPTIONAL, POINTER                      :: symbol
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: nq, lq
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: zet

    CHARACTER(len=*), PARAMETER :: routineN = 'set_sto_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ns

! -------------------------------------------------------------------------

    IF (ASSOCIATED(sto_basis_set)) THEN
       IF (PRESENT(name)) sto_basis_set%name = name
       IF (PRESENT(nshell)) sto_basis_set%nshell = nshell
       IF (PRESENT(symbol)) THEN
          ns = SIZE(symbol)
          IF (ASSOCIATED(sto_basis_set%symbol)) DEALLOCATE (sto_basis_set%symbol)
          ALLOCATE (sto_basis_set%symbol(1:ns))
          sto_basis_set%symbol(:) = symbol(:)
       END IF
       IF (PRESENT(nq)) THEN
          ns = SIZE(nq)
          CALL reallocate(sto_basis_set%nq,1,ns)
          sto_basis_set%nq = nq(:)
       END IF
       IF (PRESENT(lq)) THEN
          ns = SIZE(lq)
          CALL reallocate(sto_basis_set%lq,1,ns)
          sto_basis_set%lq = lq(:)
       END IF
       IF (PRESENT(zet)) THEN
          ns = SIZE(zet)
          CALL reallocate(sto_basis_set%zet,1,ns)
          sto_basis_set%zet = zet(:)
       END IF
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer sto_basis_set is not associated")
    END IF

  END SUBROUTINE set_sto_basis_set

! *****************************************************************************
!> \brief ...
!> \param sto_basis_set ...
!> \param gto_basis_set ...
!> \param ngauss ...
! *****************************************************************************
  SUBROUTINE create_gto_from_sto_basis(sto_basis_set,gto_basis_set,ngauss)

    TYPE(sto_basis_set_type), POINTER        :: sto_basis_set
    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set
    INTEGER, OPTIONAL                        :: ngauss

    CHARACTER(len=*), PARAMETER :: routineN = 'create_gto_from_sto_basis', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: maxng = 6

    CHARACTER(LEN=default_string_length)     :: name, sng
    INTEGER                                  :: ico, ipgf, iset, lshell, m, &
                                                maxco, maxl, ncgf, ng, np, &
                                                nset, nsgf, nshell, stat
    INTEGER, DIMENSION(:), POINTER           :: lq, nq
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: zet
    REAL(KIND=dp), DIMENSION(maxng)          :: gcc, zetg

    failure = .FALSE.
    ng = 6
    IF (PRESENT(ngauss)) ng=ngauss
    IF (ng > maxng) CALL stop_program(routineN,moduleN,__LINE__,&
         "Too many Gaussian primitives requested")

    CALL allocate_gto_basis_set(gto_basis_set)

    CALL get_sto_basis_set(sto_basis_set,name=name,nshell=nshell,nq=nq,&
         lq=lq,zet=zet)

    maxl = MAXVAL(lq)
    CALL init_orbital_pointers(maxl)

    CALL integer_to_string(ng,sng)
    gto_basis_set%name = TRIM(name)//"_STO-"//TRIM(sng)//"G"

    nset = nshell
    gto_basis_set%nset = nset
    CALL reallocate(gto_basis_set%lmax,1,nset)
    CALL reallocate(gto_basis_set%lmin,1,nset)
    CALL reallocate(gto_basis_set%npgf,1,nset)
    CALL reallocate(gto_basis_set%nshell,1,nset)
    CALL reallocate(gto_basis_set%n,1,1,1,nset)
    CALL reallocate(gto_basis_set%l,1,1,1,nset)
    CALL reallocate(gto_basis_set%zet,1,ng,1,nset)
    CALL reallocate(gto_basis_set%gcc,1,ng,1,1,1,nset)

    DO iset=1,nset
       CALL get_sto_ng ( zet(iset), ng, nq(iset), lq(iset), zetg, gcc)
       gto_basis_set%lmax(iset) = lq(iset)
       gto_basis_set%lmin(iset) = lq(iset)
       gto_basis_set%npgf(iset) = ng
       gto_basis_set%nshell(iset) = 1
       gto_basis_set%n(1,iset) = lq(iset)+1
       gto_basis_set%l(1,iset) = lq(iset)
       DO ipgf=1,ng
          gto_basis_set%gcc(ipgf,1,iset) = gcc(ipgf)
          gto_basis_set%zet(ipgf,iset) = zetg(ipgf)
       END DO
    END DO

    CALL reallocate(gto_basis_set%set_radius,1,nset)
    CALL reallocate(gto_basis_set%pgf_radius,1,ng,1,nset)
    CALL reallocate(gto_basis_set%first_cgf,1,1,1,nset)
    CALL reallocate(gto_basis_set%first_sgf,1,1,1,nset)
    CALL reallocate(gto_basis_set%last_cgf,1,1,1,nset)
    CALL reallocate(gto_basis_set%last_sgf,1,1,1,nset)
    CALL reallocate(gto_basis_set%ncgf_set,1,nset)
    CALL reallocate(gto_basis_set%nsgf_set,1,nset)

    maxco = 0
    ncgf = 0
    nsgf = 0

    DO iset=1,nset
       gto_basis_set%ncgf_set(iset) = 0
       gto_basis_set%nsgf_set(iset) = 0
       lshell = gto_basis_set%l(1,iset)
       gto_basis_set%first_cgf(1,iset) = ncgf + 1
       ncgf = ncgf + nco(lshell)
       gto_basis_set%last_cgf(1,iset) = ncgf
       gto_basis_set%ncgf_set(iset) =&
            gto_basis_set%ncgf_set(iset) + nco(lshell)
       gto_basis_set%first_sgf(1,iset) = nsgf + 1
       nsgf = nsgf + nso(lshell)
       gto_basis_set%last_sgf(1,iset) = nsgf
       gto_basis_set%nsgf_set(iset) =&
            gto_basis_set%nsgf_set(iset) + nso(lshell)
       maxco = MAX(maxco,ng*ncoset(lshell))
    END DO

    gto_basis_set%ncgf = ncgf
    gto_basis_set%nsgf = nsgf

    CALL reallocate(gto_basis_set%cphi,1,maxco,1,ncgf)
    CALL reallocate(gto_basis_set%sphi,1,maxco,1,nsgf)
    CALL reallocate(gto_basis_set%scon,1,maxco,1,nsgf)
    CALL reallocate(gto_basis_set%lx,1,ncgf)
    CALL reallocate(gto_basis_set%ly,1,ncgf)
    CALL reallocate(gto_basis_set%lz,1,ncgf)
    CALL reallocate(gto_basis_set%m,1,nsgf)
    CALL reallocate(gto_basis_set%norm_cgf,1,ncgf)
    ALLOCATE (gto_basis_set%cgf_symbol(ncgf),STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE (gto_basis_set%sgf_symbol(nsgf),STAT=stat)
    CPASSERT(stat==0)

    ncgf = 0
    nsgf = 0

    DO iset=1,nset
       lshell = gto_basis_set%l(1,iset)
       np = lshell + 1
       DO ico=ncoset(lshell-1)+1,ncoset(lshell)
          ncgf = ncgf + 1
          gto_basis_set%lx(ncgf) = indco(1,ico)
          gto_basis_set%ly(ncgf) = indco(2,ico)
          gto_basis_set%lz(ncgf) = indco(3,ico)
          gto_basis_set%cgf_symbol(ncgf) =&
               cgf_symbol(np,(/gto_basis_set%lx(ncgf),&
               gto_basis_set%ly(ncgf),&
               gto_basis_set%lz(ncgf)/))
       END DO
       DO m=-lshell,lshell
          nsgf = nsgf + 1
          gto_basis_set%m(nsgf) = m
          gto_basis_set%sgf_symbol(nsgf) = sgf_symbol(np,lshell,m)
       END DO
    END DO

    gto_basis_set%norm_type = -1

  END SUBROUTINE create_gto_from_sto_basis

! *****************************************************************************
!> \brief ...
!> \param z ...
!> \param ne ...
!> \param n ...
!> \param l ...
!> \retval srules ...
! *****************************************************************************
  FUNCTION srules(z,ne,n,l)
    ! Slater rules
    INTEGER                                  :: z
    INTEGER, DIMENSION(4, 7)                 :: ne
    INTEGER                                  :: n, l
    REAL(dp)                                 :: srules

    REAL(dp), DIMENSION(7), PARAMETER :: &
      xns = (/ 1.0_dp,2.0_dp,3.0_dp,3.7_dp,4.0_dp,4.2_dp,4.4_dp /)

    INTEGER                                  :: i, l1, l2, m, m1, m2
    REAL(dp)                                 :: s

    s = 0.0_dp
    ! The complete shell
    l1=l+1
    IF(l1 == 1) l2=2
    IF(l1 == 2) l2=1
    IF(l1 == 3) l2=4
    IF(l1 == 4) l2=3
    ! Rule a) no contribution from shells further out
    ! Rule b) 0.35 (1s 0.3) from each other electron in the same shell
    IF(n == 1) THEN
       m=ne(1,1)
       s=s+0.3_dp*REAL(m-1,dp)
    ELSE
       m=ne(l1,n)+ne(l2,n)
       s=s+0.35_dp*REAL(m-1,dp)
    END IF
    ! Rule c) if (s,p) shell 0.85 from each electron with n-1, and 1.0
    !      from all electrons further in
    IF(l1+l2 == 3) THEN
       IF(n > 1) THEN
          m1=ne(1,n-1)+ne(2,n-1)+ne(3,n-1)+ne(4,n-1)
          m2=0
          DO i=1,n-2
             m2=m2+ne(1,i)+ne(2,i)+ne(3,i)+ne(4,I)
          END DO
          s=s+0.85_dp*REAL(m1,dp)+1._dp*REAL(m2,dp)
       END IF
    ELSE
       ! Rule d) if (d,f) shell 1.0 from each electron inside
       m=0
       DO i=1,n-1
          m=m+ne(1,i)+ne(2,i)+ne(3,i)+ne(4,i)
       END DO
       s=s+1._dp*REAL(m,dp)
    END IF
    ! Slater exponent is (Z-S)/NS
    srules = (REAL(z,dp) - s)/xns(n)
  END FUNCTION srules

! *****************************************************************************

END MODULE basis_set_types
