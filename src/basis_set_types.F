!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/basis_set_types [1.0] *
!!
!!   NAME
!!     basis_set_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     Matthias Krack (04.07.2000)
!!
!!   MODIFICATION HISTORY
!!     - 02.2004 flexible normalization of basis sets [jgh]
!!
!!***
!******************************************************************************

MODULE basis_set_types

  USE ai_coulomb,                      ONLY: coulomb2
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size
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
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             search_object,&
                                             start_parser
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module) ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "basis_set_types)"

! *** Define the Gaussian-type orbital basis set type ***

  TYPE gto_basis_set_type
!MK PRIVATE
    CHARACTER(LEN=default_string_length)       :: name
    REAL(KIND = dp)                            :: kind_radius
    INTEGER                                    :: norm_type
    INTEGER                                    :: ncgf,nset,nsgf
    CHARACTER(LEN=12), DIMENSION(:), POINTER   :: cgf_symbol
    CHARACTER(LEN=6), DIMENSION(:), POINTER    :: sgf_symbol
    REAL(KIND = dp), DIMENSION(:), POINTER     :: norm_cgf,set_radius
    INTEGER, DIMENSION(:), POINTER             :: lmax,lmin,lx,ly,lz,m,ncgf_set,&
                                                  npgf,nsgf_set,nshell
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: cphi,pgf_radius,sphi,zet
    INTEGER, DIMENSION(:,:), POINTER           :: first_cgf,first_sgf,l,&
                                                  last_cgf,last_sgf,n
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: gcc
  END TYPE gto_basis_set_type

! *** Public subroutines ***

  PUBLIC :: allocate_gto_basis_set,&
            deallocate_gto_basis_set,&
            get_gto_basis_set,&
            init_aux_basis_set,&
            init_cphi_and_sphi,&
            init_orb_basis_set,&
            read_gto_basis_set,&
            set_gto_basis_set,&
            write_aux_basis_set,&
            write_gto_basis_set,&
            write_orb_basis_set

! *** Public data types ***

  PUBLIC :: gto_basis_set_type

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_gto_basis_set(gto_basis_set)

!   Purpose: Allocate a Gaussian-type orbital (GTO) basis set data set.

!   History: - Creation (26.10.2000,MK)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_gto_basis_set"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(gto_basis_set)) CALL deallocate_gto_basis_set(gto_basis_set)

    ALLOCATE (gto_basis_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set",0)

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
    NULLIFY (gto_basis_set%set_radius)
    NULLIFY (gto_basis_set%sgf_symbol)
    NULLIFY (gto_basis_set%zet)

  END SUBROUTINE allocate_gto_basis_set

! *****************************************************************************

  SUBROUTINE deallocate_gto_basis_set(gto_basis_set)

!   Purpose: Deallocate a Gaussian-type orbital (GTO) basis set data set.

!   History: - Creation (03.11.2000,MK)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_gto_basis_set"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    DEALLOCATE (gto_basis_set%cgf_symbol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%cgf_symbol")

    DEALLOCATE (gto_basis_set%sgf_symbol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%sgf_symbol")

    DEALLOCATE (gto_basis_set%norm_cgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%norm_cgf")

    DEALLOCATE (gto_basis_set%set_radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%set_radius")

    DEALLOCATE (gto_basis_set%lmax,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%lmax")

    DEALLOCATE (gto_basis_set%lmin,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%lmin")

    DEALLOCATE (gto_basis_set%lx,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%lx")

    DEALLOCATE (gto_basis_set%ly,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%ly")

    DEALLOCATE (gto_basis_set%lz,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%lz")

    DEALLOCATE (gto_basis_set%m,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%m")

    DEALLOCATE (gto_basis_set%ncgf_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%ncgf_set")

    DEALLOCATE (gto_basis_set%npgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%npgf")

    DEALLOCATE (gto_basis_set%nsgf_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%nsgf_set")

    DEALLOCATE (gto_basis_set%nshell,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%nshell")

    DEALLOCATE (gto_basis_set%cphi,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%cphi")

    DEALLOCATE (gto_basis_set%pgf_radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%pgf_radius")

    DEALLOCATE (gto_basis_set%sphi,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%sphi")

    DEALLOCATE (gto_basis_set%zet,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%zet")

    DEALLOCATE (gto_basis_set%first_cgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%first_cgf")

    DEALLOCATE (gto_basis_set%first_sgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%first_sgf")

    DEALLOCATE (gto_basis_set%l,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%l")

    DEALLOCATE (gto_basis_set%last_cgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%last_cgf")

    DEALLOCATE (gto_basis_set%last_sgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%last_sgf")

    DEALLOCATE (gto_basis_set%n,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%n")

    DEALLOCATE (gto_basis_set%gcc,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set%gcc")

    DEALLOCATE (gto_basis_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gto_basis_set")

  END SUBROUTINE deallocate_gto_basis_set

! *****************************************************************************

  SUBROUTINE get_gto_basis_set(gto_basis_set,name,norm_type,kind_radius,ncgf,&
                    nset,nsgf,cgf_symbol,sgf_symbol,norm_cgf,set_radius,lmax,&
                    lmin,lx,ly,lz,m,ncgf_set,npgf,nsgf_set,nshell,&
                    cphi,pgf_radius,sphi,zet,first_cgf,first_sgf,l,&
                    last_cgf,last_sgf,n,gcc,maxco,maxl,maxpgf,&
                    maxsgf_set,maxshell,maxso,nco_sum,npgf_sum,&
                    nshell_sum,maxder)

!   Purpose: Get informations about a Gaussian-type orbital (GTO) basis set.

!   History: - Creation (10.01.2002,MK)

!   ***************************************************************************

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
      OPTIONAL, POINTER                      :: cphi, pgf_radius, sphi, zet
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: first_cgf, first_sgf, l, &
                                                last_cgf, last_sgf, n
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: gcc
    INTEGER, INTENT(OUT), OPTIONAL           :: maxco, maxl, maxpgf, &
                                                maxsgf_set, maxshell, maxso,&
                                                nco_sum, npgf_sum, nshell_sum
    INTEGER, INTENT(IN), OPTIONAL            :: maxder

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_gto_basis_set"

    INTEGER                                  :: iset, nder

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(gto_basis_set)) THEN
      IF (PRESENT(name)) name = gto_basis_set%name
      IF (PRESENT(norm_type)) norm_type = gto_basis_set%norm_type
      IF (PRESENT(kind_radius)) kind_radius = gto_basis_set%kind_radius
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
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer gto_basis_set is not associated")
    END IF

  END SUBROUTINE get_gto_basis_set

! *****************************************************************************

  SUBROUTINE init_aux_basis_set(gto_basis_set)

!   Purpose: Initialise a Gaussian-type orbital (GTO) basis set data set.

!   History: - Creation (06.12.2000,MK)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_aux_basis_set"

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(gto_basis_set)) RETURN

    SELECT CASE (gto_basis_set%norm_type)
    CASE ( 0 )
      !   *** No normalisation requested ***
    CASE ( 1 )
      CALL init_norm_cgf_aux_2(gto_basis_set)
    CASE ( 2 )
! WARNING this was never tested
      CALL init_norm_cgf_aux(gto_basis_set)
    CASE DEFAULT
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "Normalization method not specified")
    END SELECT

!   *** Initialise the transformation matrices "pgf" -> "cgf" ***

    CALL init_cphi_and_sphi(gto_basis_set)

  END SUBROUTINE init_aux_basis_set

! *****************************************************************************

  SUBROUTINE init_cphi_and_sphi(gto_basis_set)

!   Purpose: Initialise the matrices for the transformation of primitive
!            Cartesian Gaussian-type functions to contracted Cartesian
!            (cphi) and spherical (sphi) Gaussian-type functions.

!   History: - Creation (20.09.2000,MK)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    INTEGER                                  :: first_cgf, first_sgf, icgf, &
                                                ico, ipgf, iset, ishell, l, &
                                                n, ncgf, nsgf

!   ---------------------------------------------------------------------------
!   *** Build the Cartesian transformation matrix "cphi" ***

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

!   *** Build the spherical transformation matrix "sphi" ***

    n = SIZE(gto_basis_set%cphi,1)

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

  END SUBROUTINE init_cphi_and_sphi

! *****************************************************************************

  SUBROUTINE init_norm_cgf_aux(gto_basis_set)

!   Purpose: Initialise the normalization factors of the contracted Cartesian
!            Gaussian functions.

!   History: - Creation (07.12.2000,MK)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_norm_cgf_aux"

    INTEGER                                  :: icgf, ico, ipgf, iset, &
                                                ishell, istat, jco, jpgf, &
                                                lmax, lmin, lx, ly, lz, n, &
                                                npgfa
    REAL(KIND=dp)                            :: fnorm, gcca, gccb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: gaa
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rpgfa, zeta

!   ---------------------------------------------------------------------------

    n = 0

    DO iset=1,gto_basis_set%nset
      n = MAX(n,gto_basis_set%npgf(iset)*ncoset(gto_basis_set%lmax(iset)))
    END DO

    ALLOCATE (gaa(n,n),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gaa",n*n*dp_size)

    DO iset=1,gto_basis_set%nset
      lmax = gto_basis_set%lmax(iset)
      lmin = gto_basis_set%lmin(iset)
      n = ncoset(lmax)
      npgfa = gto_basis_set%npgf(iset)
      rpgfa => gto_basis_set%pgf_radius(1:npgfa,iset)
      zeta => gto_basis_set%zet(1:npgfa,iset)
      CALL coulomb2(lmax,zeta,rpgfa,lmin,&
                    lmax,zeta,rpgfa,lmin,&
                    (/0.0_dp,0.0_dp,0.0_dp/),0.0_dp,gaa)
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

    DEALLOCATE (gaa,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gaa")

  END SUBROUTINE init_norm_cgf_aux

! *****************************************************************************

  SUBROUTINE init_norm_cgf_aux_2(gto_basis_set)

!   Purpose: Initialise the normalization factors of the auxiliary Cartesian
!            Gaussian functions (kim-gordon polarization basis) Norm=1.

!   History: - Creation (07.12.2000,GT)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_norm_cgf_aux_2"

    INTEGER                                  :: icgf, iset, ishell

!   ---------------------------------------------------------------------------

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

  SUBROUTINE init_norm_cgf_orb(gto_basis_set)

!   Purpose: Initialise the normalization factors of the contracted Cartesian
!            Gaussian functions.

!   History: - Creation (14.04.2000,MK)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    INTEGER                                  :: icgf, ipgf, iset, ishell, &
                                                jpgf, l, lx, ly, lz
    REAL(KIND=dp)                            :: expzet, fnorm, gcca, gccb, &
                                                prefac, zeta, zetb

!   ---------------------------------------------------------------------------

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

  SUBROUTINE init_norm_cgf_orb_den(gto_basis_set)

!   Purpose: Initialise the normalization factors of the contracted Cartesian
!            Gaussian functions used for frozen density representation

!   History: - Creation (21.09.2002,GT)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_norm_cgf_orb_den"

    INTEGER                                  :: icgf, ipgf, iset, ishell, l
    REAL(KIND=dp)                            :: expzet, gcca, prefac, zeta

!   ---------------------------------------------------------------------------

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

  SUBROUTINE init_orb_basis_set(gto_basis_set)

!   Purpose: Initialise a Gaussian-type orbital (GTO) basis set data set.

!   History: - Creation (26.10.2000,MK)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_orb_basis_set"

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(gto_basis_set)) RETURN

    SELECT CASE (gto_basis_set%norm_type)
    CASE ( 0 )
      !   *** No normalisation requested ***
    CASE ( 1 )
      CALL init_norm_cgf_orb_den(gto_basis_set)
    CASE ( 2 )
      !   *** Normalise the primitive Gaussian functions ***
      CALL normalise_gcc_orb(gto_basis_set)
      !   *** Compute the normalization factors of   ***
      !   *** the contracted Gaussian-type functions ***
      CALL init_norm_cgf_orb(gto_basis_set)
    CASE DEFAULT
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "Normalization method not specified")
    END SELECT 

!   *** Initialise the transformation matrices "pgf" -> "cgf" ***
    CALL init_cphi_and_sphi(gto_basis_set)

  END SUBROUTINE init_orb_basis_set

! *****************************************************************************

  SUBROUTINE normalise_gcc_orb(gto_basis_set)

!   Purpose: Normalise the primitive Cartesian Gaussian functions. The
!            normalization factor is included in the Gaussian contraction
!            coefficients.

!   History: - Creation (20.08.1999,MK)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set

    INTEGER                                  :: ipgf, iset, ishell, l
    REAL(KIND=dp)                            :: expzet, gcca, prefac, zeta

!   ---------------------------------------------------------------------------

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

  SUBROUTINE read_gto_basis_set(element_symbol,basis_set_name,gto_basis_set,&
                                globenv)

!   Purpose: Read a Gaussian-type orbital (GTO) basis set from the database
!            file.

!   History: - Creation (13.04.2000,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: element_symbol, basis_set_name
    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_gto_basis_set"

    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(LEN=LEN(basis_set_name))       :: bsname
    CHARACTER(LEN=LEN(basis_set_name)+2)     :: bsname2
    CHARACTER(LEN=LEN(element_symbol))       :: symbol
    CHARACTER(LEN=LEN(element_symbol)+2)     :: symbol2
    INTEGER :: i, ico, ipgf, iset, ishell, istat, lshell, m, maxco, maxl, &
      maxpgf, maxshell, ncgf, nmin, nset, nsgf, strlen1, strlen2
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, n
    LOGICAL                                  :: found
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zet
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: gcc

!   ---------------------------------------------------------------------------

    line2 = ""
    symbol2 = ""
    bsname2 = ""

    gto_basis_set%name = basis_set_name

    CALL start_parser(globenv%basis_set_file_name,globenv)

!   *** Search for the requested basis set in the basis set file   ***
!   *** until the basis set is found or the end of file is reached ***

    bsname = basis_set_name
    symbol = element_symbol

    search_loop: DO

      CALL search_object(TRIM(bsname),.TRUE.,found,line)

      IF (found) THEN

        CALL uppercase(line)
        CALL uppercase(symbol)
        CALL uppercase(bsname)

!       *** Check both the element symbol and the basis set name ***

        line2 = " "//line//" "
        symbol2 = " "//TRIM(symbol)//" "
        bsname2 = " "//TRIM(bsname)//" "
        strlen1 = LEN_TRIM(symbol2) + 1
        strlen2 = LEN_TRIM(bsname2) + 1

        IF ((INDEX(line2,symbol2(:strlen1)) > 0).AND.&
            (INDEX(line2,bsname2(:strlen2)) > 0)) THEN

          NULLIFY (gcc,l,lmax,lmin,n,npgf,nshell,zet)

!         *** Read the basis set information ***

          CALL read_object(nset,newline=.TRUE.)

          CALL reallocate(npgf,1,nset)
          CALL reallocate(nshell,1,nset)
          CALL reallocate(lmax,1,nset)
          CALL reallocate(lmin,1,nset)
          CALL reallocate(n,1,1,1,nset)

          maxl = 0
          maxpgf = 0
          maxshell = 0

          DO iset=1,nset
            CALL read_object(n(1,iset),newline=.TRUE.)
            CALL read_object(lmin(iset))
            CALL read_object(lmax(iset))
            maxl = MAX(maxl,lmax(iset))
            CALL read_object(npgf(iset))
            IF (npgf(iset) > maxpgf) THEN
              maxpgf = npgf(iset)
              CALL reallocate(zet,1,maxpgf,1,nset)
              CALL reallocate(gcc,1,maxpgf,1,maxshell,1,nset)
            END IF
            nshell(iset) = 0
            DO lshell=lmin(iset),lmax(iset)
              nmin = n(1,iset) + lshell - lmin(iset)
              CALL read_object(ishell)
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
              CALL read_object(zet(ipgf,iset),newline=.TRUE.)
              DO ishell=1,nshell(iset)
                CALL read_object(gcc(ipgf,ishell,iset))
              END DO
            END DO
          END DO

!         *** Maximum angular momentum quantum number of the atomic kind ***

          CALL init_orbital_pointers(maxl)

!         *** Allocate the global variables ***

          gto_basis_set%nset = nset
          CALL reallocate(gto_basis_set%lmax,1,nset)
          CALL reallocate(gto_basis_set%lmin,1,nset)
          CALL reallocate(gto_basis_set%npgf,1,nset)
          CALL reallocate(gto_basis_set%nshell,1,nset)
          CALL reallocate(gto_basis_set%n,1,maxshell,1,nset)
          CALL reallocate(gto_basis_set%l,1,maxshell,1,nset)
          CALL reallocate(gto_basis_set%zet,1,maxpgf,1,nset)
          CALL reallocate(gto_basis_set%gcc,1,maxpgf,1,maxshell,1,nset)

!         *** Copy the basis set information into the data structure ***

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

!         *** Initialise the depending atomic kind information ***

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
          CALL reallocate(gto_basis_set%lx,1,ncgf)
          CALL reallocate(gto_basis_set%ly,1,ncgf)
          CALL reallocate(gto_basis_set%lz,1,ncgf)
          CALL reallocate(gto_basis_set%m,1,nsgf)
          CALL reallocate(gto_basis_set%norm_cgf,1,ncgf)
          ALLOCATE (gto_basis_set%cgf_symbol(ncgf),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine_name,module_name,__LINE__,&
                             "gto_basis_set%cgf_symbol",12*ncgf)
          END IF
          ALLOCATE (gto_basis_set%sgf_symbol(nsgf),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine_name,module_name,__LINE__,&
                             "gto_basis_set%sgf_symbol",6*nsgf)
          END IF

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

          DEALLOCATE (gcc,l,lmax,lmin,n,npgf,nshell,zet,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine_name,module_name,__LINE__,&
                             "gcc,l,lmax,lmin,n,npgf,nshell,zet")
          END IF

          EXIT search_loop

        END IF

      ELSE

!       *** Stop program, if the end of file is reached ***

        CALL stop_program(routine_name,module_name,__LINE__,&
                          "The requested basis set <"//TRIM(bsname)//&
                          "> for element <"//TRIM(symbol)//"> was not "//&
                          "found in the basis set file <"//&
                          TRIM(globenv%basis_set_file_name)//">",globenv)

      END IF

    END DO search_loop

    CALL finish_parser()

  END SUBROUTINE read_gto_basis_set

! *****************************************************************************

  SUBROUTINE set_gto_basis_set(gto_basis_set,name,norm_type,kind_radius,ncgf,&
                     nset,nsgf,cgf_symbol,sgf_symbol,norm_cgf,set_radius,lmax,&
                     lmin,lx,ly,lz,m,ncgf_set,npgf,nsgf_set,nshell,&
                     cphi,pgf_radius,sphi,zet,first_cgf,first_sgf,l,&
                     last_cgf,last_sgf,n,gcc)

!   Purpose: Set the components of Gaussian-type orbital (GTO) basis set data
!            set.

!   History: - Creation (10.01.2002,MK)

!   ***************************************************************************

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
      OPTIONAL, POINTER                      :: cphi, pgf_radius, sphi, zet
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: first_cgf, first_sgf, l, &
                                                last_cgf, last_sgf, n
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: gcc

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_gto_basis_set"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(gto_basis_set)) THEN
      IF (PRESENT(name)) gto_basis_set%name = name
      IF (PRESENT(norm_type)) gto_basis_set%norm_type = norm_type
      IF (PRESENT(kind_radius)) gto_basis_set%kind_radius = kind_radius
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
      IF (PRESENT(zet)) gto_basis_set%zet(:,:) = zet(:,:)
      IF (PRESENT(first_cgf)) gto_basis_set%first_cgf(:,:) = first_cgf(:,:)
      IF (PRESENT(first_sgf)) gto_basis_set%first_sgf(:,:) = first_sgf(:,:)
      IF (PRESENT(l)) l(:,:) = gto_basis_set%l(:,:)
      IF (PRESENT(last_cgf)) gto_basis_set%last_cgf(:,:) = last_cgf(:,:)
      IF (PRESENT(last_sgf)) gto_basis_set%last_sgf(:,:) = last_sgf(:,:)
      IF (PRESENT(n)) gto_basis_set%n(:,:) = n(:,:)
      IF (PRESENT(gcc)) gto_basis_set%gcc(:,:,:) = gcc(:,:,:)
    ELSE
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer gto_basis_set is not associated")
    END IF

  END SUBROUTINE set_gto_basis_set

! *****************************************************************************

  SUBROUTINE write_aux_basis_set(aux_basis_set,globenv)

!   Purpose: Write an auxiliary Gaussian-type orbital (GTO) basis set data set
!            to the output unit.

!   History: - Creation (09.01.2002,MK)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: icgf, ico, ipgf, iset, &
                                                ishell, output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (ASSOCIATED(aux_basis_set).AND.globenv%print%basis_set) THEN

      output_unit = globenv%scr

      WRITE (UNIT=output_unit,FMT="(/,T4,A,/,/,T6,A)")&
        "- Auxiliary basis set information",&
        "Auxiliary basis set name: "//TRIM(aux_basis_set%name)

      WRITE (UNIT=output_unit,FMT="(/,(T6,A,I6))")&
        "Number of orbital shell sets:            ",&
        aux_basis_set%nset,&
        "Number of orbital shells:                ",&
        SUM(aux_basis_set%nshell(:)),&
        "Number of primitive Cartesian functions: ",&
        SUM(aux_basis_set%npgf(:)),&
        "Number of Cartesian basis functions:     ",&
        aux_basis_set%ncgf,&
        "Number of spherical basis functions:     ",&
        aux_basis_set%nsgf

      WRITE (UNIT=output_unit,FMT="(/,T6,A,/,/,T6,A)")&
        "Normalised Cartesian auxiliary functions",&
        "Set   Shell   Function           Exponent    Coefficient"

      icgf = 0

      DO iset=1,aux_basis_set%nset
        DO ishell=1,aux_basis_set%nshell(iset)
          WRITE (UNIT=output_unit,FMT="(A)") ""
          DO ico=1,nco(aux_basis_set%l(ishell,iset))
            icgf = icgf + 1
            WRITE (UNIT=output_unit,&
                   FMT="(T6,I3,4X,I4,3X,A12,(T32,2F15.6))")&
              iset,ishell,&
              aux_basis_set%cgf_symbol(icgf),&
              (aux_basis_set%zet(ipgf,iset),&
               aux_basis_set%norm_cgf(icgf)*&
               aux_basis_set%gcc(ipgf,ishell,iset),&
               ipgf=1,aux_basis_set%npgf(iset))
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE write_aux_basis_set

! *****************************************************************************

  SUBROUTINE write_gto_basis_set(gto_basis_set,globenv)

!   Purpose: Write a Gaussian-type orbital (GTO) basis set data set to the
!            output unit.

!   History: - Creation (09.01.2002,MK)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: gto_basis_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: ipgf, iset, ishell, &
                                                output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (ASSOCIATED(gto_basis_set).AND.globenv%print%basis_set) THEN

      output_unit = globenv%scr

      WRITE (UNIT=output_unit,FMT="(/,T4,A,2(/,/,T6,A))")&
        "- GTO basis set information",&
        "GTO basis set name: "//TRIM(gto_basis_set%name),&
        "Set   Shell     n   l            Exponent    Coefficient"

      DO iset=1,gto_basis_set%nset
        WRITE (UNIT=output_unit,FMT="(A)") ""
        DO ishell=1,gto_basis_set%nshell(iset)
          WRITE (UNIT=output_unit,&
                 FMT="(T6,I3,4X,I4,4X,I2,2X,I2,(T32,2F15.6))")&
            iset,ishell,&
            gto_basis_set%n(ishell,iset),&
            gto_basis_set%l(ishell,iset),&
            (gto_basis_set%zet(ipgf,iset),gto_basis_set%gcc(ipgf,ishell,iset),&
             ipgf=1,gto_basis_set%npgf(iset))
        END DO
      END DO

    END IF

  END SUBROUTINE write_gto_basis_set

! *****************************************************************************

  SUBROUTINE write_orb_basis_set(orb_basis_set,globenv)

!   Purpose: Write a Gaussian-type orbital (GTO) basis set data set to the
!            output unit.

!   History: - Creation (09.01.2002,MK)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: icgf, ico, ipgf, iset, &
                                                ishell, output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (ASSOCIATED(orb_basis_set).AND.globenv%print%basis_set) THEN

      output_unit = globenv%scr

      WRITE (UNIT=output_unit,FMT="(/,T4,A,/,/,T6,A)")&
        "- Orbital basis set information",&
        "Orbital basis set name: "//TRIM(orb_basis_set%name)

      WRITE (UNIT=output_unit,FMT="(/,(T6,A,I6))")&
        "Number of orbital shell sets:            ",&
        orb_basis_set%nset,&
        "Number of orbital shells:                ",&
        SUM(orb_basis_set%nshell(:)),&
        "Number of primitive Cartesian functions: ",&
        SUM(orb_basis_set%npgf(:)),&
        "Number of Cartesian basis functions:     ",&
        orb_basis_set%ncgf,&
        "Number of spherical basis functions:     ",&
        orb_basis_set%nsgf

      WRITE (UNIT=output_unit,FMT="(/,T6,A,/,/,T6,A)")&
        "Normalised Cartesian orbitals",&
          "Set   Shell   Orbital            Exponent    Coefficient"

      icgf = 0

      DO iset=1,orb_basis_set%nset
        DO ishell=1,orb_basis_set%nshell(iset)
          WRITE (UNIT=output_unit,FMT="(A)") ""
          DO ico=1,nco(orb_basis_set%l(ishell,iset))
            icgf = icgf + 1
            WRITE (UNIT=output_unit,FMT="(T6,I3,4X,I4,3X,A12,(T32,2F15.6))")&
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

END MODULE basis_set_types
