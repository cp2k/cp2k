!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE basis_set_types

! *****************************************************************************

! Index:

! SUBROUTINE allocate_gto_basis_set(gto_basis_set)
! SUBROUTINE deallocate_gto_basis_set(gto_basis_set)
! SUBROUTINE init_cphi_and_sphi(gto_basis_set)
! SUBROUTINE init_norm_cgf(gto_basis_set)
! SUBROUTINE init_aux_basis_set(gto_basis_set)
! SUBROUTINE init_orb_basis_set(gto_basis_set)
! SUBROUTINE normalize_gcc_aux(gto_basis_set)
! SUBROUTINE normalize_gcc_orb(gto_basis_set)
! SUBROUTINE read_gto_basis_set(element_symbol,basis_set_name,gto_basis_set,&
!                               globenv)
! SUBROUTINE write_gto_basis_set(gto_basis_set,lunit)

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  IMPLICIT NONE

  PRIVATE

! *** Define the Gaussian-type orbital basis set type ***

  TYPE gto_basis_set_type
    REAL(wp)                                 :: kind_radius
    INTEGER                                  :: ncgf,nset,nsgf
    CHARACTER(LEN=12), DIMENSION(:), POINTER :: cgf_symbol
    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: sgf_symbol
    REAL(wp), DIMENSION(:), POINTER          :: norm_cgf,set_radius
    INTEGER, DIMENSION(:), POINTER           :: lmax,lmin,lx,ly,lz,m,npgf,&
                                                nshell
    REAL(wp), DIMENSION(:,:), POINTER        :: cphi,pgf_radius,sphi,zet
    INTEGER, DIMENSION(:,:), POINTER         :: first_cgf,first_sgf,l,&
                                                last_cgf,last_sgf,n
    REAL(wp), DIMENSION(:,:,:), POINTER      :: gcc
  END TYPE gto_basis_set_type

  INTEGER :: maxlcgf = 0

! *** Public subroutines ***

  PUBLIC :: allocate_basis_set,&
            deallocate_basis_set,&
            init_aux_basis_set,&
            init_orb_basis_set,&
            read_basis_set,&
            write_basis_set

! *** Public data types ***

  PUBLIC :: gto_basis_set_type

! *** Public variables ***

  PUBLIC :: maxlcgf

! *****************************************************************************

  INTERFACE allocate_basis_set
    MODULE PROCEDURE allocate_gto_basis_set
  END INTERFACE

  INTERFACE deallocate_basis_set
    MODULE PROCEDURE deallocate_gto_basis_set
  END INTERFACE

  INTERFACE read_basis_set
    MODULE PROCEDURE read_gto_basis_set
  END INTERFACE

  INTERFACE write_basis_set
    MODULE PROCEDURE write_gto_basis_set
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_gto_basis_set(gto_basis_set)

!   Purpose: Allocate a Gaussian-type orbital basis set data structure.

!   History: - Creation (26.10.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(gto_basis_set_type), POINTER :: gto_basis_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_gto_basis_set (MODULE basis_set_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(gto_basis_set)) CALL deallocate_gto_basis_set(gto_basis_set)

    ALLOCATE (gto_basis_set,STAT=istat)

    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set",0)

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
    NULLIFY (gto_basis_set%norm_cgf)
    NULLIFY (gto_basis_set%npgf)
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

!   Purpose: Deallocate a Gaussian-type orbital basis set data structure.

!   History: - Creation (03.11.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(gto_basis_set_type), POINTER :: gto_basis_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_gto_basis_set (MODULE basis_set_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    DEALLOCATE (gto_basis_set%cgf_symbol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%cgf_symbol")

    DEALLOCATE (gto_basis_set%sgf_symbol,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%sgf_symbol")

    DEALLOCATE (gto_basis_set%norm_cgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%norm_cgf")

    DEALLOCATE (gto_basis_set%set_radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%set_radius")

    DEALLOCATE (gto_basis_set%lmax,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%lmax")

    DEALLOCATE (gto_basis_set%lmin,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%lmin")

    DEALLOCATE (gto_basis_set%lx,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%lx")

    DEALLOCATE (gto_basis_set%ly,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%ly")

    DEALLOCATE (gto_basis_set%lz,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%lz")

    DEALLOCATE (gto_basis_set%m,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%m")

    DEALLOCATE (gto_basis_set%npgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%npgf")

    DEALLOCATE (gto_basis_set%nshell,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%nshell")

    DEALLOCATE (gto_basis_set%cphi,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%cphi")

    DEALLOCATE (gto_basis_set%pgf_radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%pgf_radius")

    DEALLOCATE (gto_basis_set%sphi,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%sphi")

    DEALLOCATE (gto_basis_set%zet,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%zet")

    DEALLOCATE (gto_basis_set%first_cgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%first_cgf")

    DEALLOCATE (gto_basis_set%first_sgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%first_sgf")

    DEALLOCATE (gto_basis_set%l,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%l")

    DEALLOCATE (gto_basis_set%last_cgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%last_cgf")

    DEALLOCATE (gto_basis_set%last_sgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%last_sgf")

    DEALLOCATE (gto_basis_set%n,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%n")

    DEALLOCATE (gto_basis_set%gcc,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set%gcc")

    DEALLOCATE (gto_basis_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"gto_basis_set")

  END SUBROUTINE deallocate_gto_basis_set

! *****************************************************************************

  SUBROUTINE init_aux_basis_set(gto_basis_set)

!   Purpose: Initialize a Gaussian-type orbital basis set data structure.

!   History: - Creation (06.12.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER :: gto_basis_set

!   ---------------------------------------------------------------------------

!   *** Compute the normalization factors of   ***
!   *** the contracted Gaussian-type functions ***

    CALL init_norm_cgf_aux(gto_basis_set)

!   *** Initialize the transformation matrices "pgf" -> "cgf" ***

    CALL init_cphi_and_sphi(gto_basis_set)

  END SUBROUTINE init_aux_basis_set

! *****************************************************************************

  SUBROUTINE init_cphi_and_sphi(gto_basis_set)

!   Purpose: Initialize the matrices for the transformation of primitive
!            Cartesian Gaussian-type functions to contracted Cartesian
!            (cphi) and spherical (sphi) Gaussian-type functions.

!   History: - Creation (20.09.2000, Matthias Krack)

!   ***************************************************************************

    USE orbital_pointers,                ONLY: coset,nco,ncoset,nso
    USE orbital_transformation_matrices, ONLY: orbtramat

    TYPE(gto_basis_set_type), POINTER :: gto_basis_set

!   *** Local variables ***

    INTEGER :: first_cgf,first_sgf,icgf,ico,ipgf,iset,ishell,l,n,ncgf,nsgf

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
                   1.0_wp,gto_basis_set%cphi(1,first_cgf),n,&
                   orbtramat(l)%c2s(1,1),nsgf,&
                   0.0_wp,gto_basis_set%sphi(1,first_sgf),n)
      END DO
    END DO

  END SUBROUTINE init_cphi_and_sphi

! *****************************************************************************

  SUBROUTINE init_norm_cgf_aux(gto_basis_set)

!   Purpose: Initialize the normalization factors of the contracted Cartesian
!            Gaussian functions.

!   History: - Creation (07.12.2000, Matthias Krack)

!   ***************************************************************************

    USE ai_coulomb,       ONLY: coulomb2
    USE memory_utilities, ONLY: reallocate
    USE orbital_pointers, ONLY: coset,ncoset

    TYPE(gto_basis_set_type), POINTER :: gto_basis_set

!   *** Local variables ***

    REAL(wp) :: fnorm,gcca,gccb
    INTEGER  :: ico,icgf,ipgf,iset,ishell,jco,jpgf,lmax,lmin,lx,ly,lz,n,npgfa

    REAL(wp), DIMENSION(:), POINTER   :: rpgfa,zeta
    REAL(wp), DIMENSION(:,:), POINTER :: gaa

!   ---------------------------------------------------------------------------

    n = 0

    DO iset=1,gto_basis_set%nset
      n = MAX(n,gto_basis_set%npgf(iset)*ncoset(gto_basis_set%lmax(iset)))
    END DO

    NULLIFY (gaa)

    gaa => reallocate(gaa,1,n,1,n)

    DO iset=1,gto_basis_set%nset
      lmax = gto_basis_set%lmax(iset)
      lmin = gto_basis_set%lmin(iset)
      n = ncoset(lmax)
      npgfa = gto_basis_set%npgf(iset)
      rpgfa => gto_basis_set%pgf_radius(1:npgfa,iset)
      zeta => gto_basis_set%zet(1:npgfa,iset)
      CALL coulomb2(lmax,zeta,rpgfa,lmin,&
                    lmax,zeta,rpgfa,lmin,&
                    (/0.0_wp,0.0_wp,0.0_wp/),0.0_wp,gaa)
      DO ishell=1,gto_basis_set%nshell(iset)
        DO icgf=gto_basis_set%first_cgf(ishell,iset),&
                gto_basis_set%last_cgf(ishell,iset)
          lx = gto_basis_set%lx(icgf)
          ly = gto_basis_set%ly(icgf)
          lz = gto_basis_set%lz(icgf)
          ico = coset(lx,ly,lz)
          fnorm = 0.0_wp
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
          gto_basis_set%norm_cgf(icgf) = 1.0_wp/SQRT(fnorm)
        END DO
      END DO
    END DO

  END SUBROUTINE init_norm_cgf_aux

! *****************************************************************************

  SUBROUTINE init_norm_cgf_orb(gto_basis_set)

!   Purpose: Initialize the normalization factors of the contracted Cartesian
!            Gaussian functions.

!   History: - Creation (14.04.2000, Matthias Krack)

!   ***************************************************************************

    USE mathconstants, ONLY: dfac,pi

    TYPE(gto_basis_set_type), POINTER :: gto_basis_set

!   *** Local variables ***

    REAL(wp) :: expzet,fnorm,gcca,gccb,prefac,zeta,zetb
    INTEGER  :: icgf,ipgf,iset,ishell,jpgf,l,lx,ly,lz

!   ---------------------------------------------------------------------------

    DO iset=1,gto_basis_set%nset
      DO ishell=1,gto_basis_set%nshell(iset)

        l = gto_basis_set%l(ishell,iset)

        expzet = 0.5_wp*REAL(2*l + 3,wp)

        fnorm = 0.0_wp

        DO ipgf=1,gto_basis_set%npgf(iset)
          gcca = gto_basis_set%gcc(ipgf,ishell,iset)
          zeta = gto_basis_set%zet(ipgf,iset)
          DO jpgf=1,gto_basis_set%npgf(iset)
            gccb = gto_basis_set%gcc(jpgf,ishell,iset)
            zetb = gto_basis_set%zet(jpgf,iset)
            fnorm = fnorm + gcca*gccb/(zeta + zetb)**expzet
          END DO
        END DO

        fnorm = 0.5_wp**l*pi**1.5_wp*fnorm

        DO icgf=gto_basis_set%first_cgf(ishell,iset),&
                gto_basis_set%last_cgf(ishell,iset)
          lx = gto_basis_set%lx(icgf)
          ly = gto_basis_set%ly(icgf)
          lz = gto_basis_set%lz(icgf)
          prefac = dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1)
          gto_basis_set%norm_cgf(icgf) = 1.0_wp/SQRT(prefac*fnorm)
        END DO

      END DO
    END DO

  END SUBROUTINE init_norm_cgf_orb

! *****************************************************************************

  SUBROUTINE init_orb_basis_set(gto_basis_set)

!   Purpose: Initialize a Gaussian-type orbital basis set data structure.

!   History: - Creation (26.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER :: gto_basis_set

!   ---------------------------------------------------------------------------

!   *** Normalize the primitive Gaussian functions ***

    CALL normalize_gcc_orb(gto_basis_set)

!   *** Compute the normalization factors of   ***
!   *** the contracted Gaussian-type functions ***

    CALL init_norm_cgf_orb(gto_basis_set)

!   *** Initialize the transformation matrices "pgf" -> "cgf" ***

    CALL init_cphi_and_sphi(gto_basis_set)

  END SUBROUTINE init_orb_basis_set

! *****************************************************************************

  SUBROUTINE normalize_gcc_orb(gto_basis_set)

!   Purpose: Normalize the primitive Cartesian Gaussian functions. The
!            normalization factor is included in the Gaussian contraction
!            coefficients.

!   History: - Creation (20.08.1999, Matthias Krack)

!   ***************************************************************************

    USE mathconstants, ONLY: pi

    TYPE(gto_basis_set_type), POINTER :: gto_basis_set

!   *** Local variables ***

    REAL(wp) :: expzet,gcca,prefac,zeta
    INTEGER  :: ipgf,iset,ishell,l

!   ---------------------------------------------------------------------------

    DO iset=1,gto_basis_set%nset
      DO ishell=1,gto_basis_set%nshell(iset)
        l = gto_basis_set%l(ishell,iset)
        expzet = 0.25_wp*REAL(2*l + 3,wp)
        prefac = 2.0_wp**l*(2.0_wp/pi)**0.75_wp
        DO ipgf=1,gto_basis_set%npgf(iset)
          gcca = gto_basis_set%gcc(ipgf,ishell,iset)
          zeta = gto_basis_set%zet(ipgf,iset)
          gto_basis_set%gcc(ipgf,ishell,iset) = prefac*zeta**expzet*gcca
        END DO
      END DO
    END DO

  END SUBROUTINE normalize_gcc_orb

! *****************************************************************************

  SUBROUTINE read_gto_basis_set(element_symbol,basis_set_name,gto_basis_set,&
                                globenv)

!   Purpose: Read a Gaussian-type orbital basis set from the database file.

!   History: - Creation (13.04.2000, Matthias Krack)

!   ***************************************************************************

    USE global_types,          ONLY: global_environment_type
    USE input_utilities,       ONLY: finish_input_session,&
                                     read_object,&
                                     search,&
                                     start_input_session
    USE memory_utilities,      ONLY: reallocate
    USE method_specifications, ONLY: maxder
    USE orbital_pointers,      ONLY: init_orbital_pointers,nco,ncoset,nso
    USE orbital_symbols,       ONLY: cgf_symbol,sgf_symbol
    USE string_utilities,      ONLY: uppercase
    USE termination,           ONLY: stop_memory,stop_program

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(gto_basis_set_type), POINTER         :: gto_basis_set
    CHARACTER(LEN=*), INTENT(IN)              :: basis_set_name,element_symbol

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_gto_basis_set (MODULE basis_set_types)"

!   *** Local variables ***

    CHARACTER(LEN=LEN(basis_set_name)) :: bsname
    CHARACTER(LEN=240)                 :: line
    CHARACTER(LEN=LEN(element_symbol)) :: symbol

    INTEGER :: i,ipgf,iset,ishell,istat,lshell,lx,ly,lz,m,maxco,maxl,maxpgf,&
               maxshell,ncgf,ncgf_set,nmin,nset,nsgf
    LOGICAL :: found

    INTEGER, DIMENSION(:), POINTER      :: lmax,lmin,npgf,nshell
    REAL(wp), DIMENSION(:,:), POINTER   :: zet
    INTEGER, DIMENSION(:,:), POINTER    :: l,n
    REAL(wp), DIMENSION(:,:,:), POINTER :: gcc

!   ---------------------------------------------------------------------------

    CALL start_input_session(globenv%basis_set_file_name,globenv)

!   *** Search for the requested basis set in the basis set file   ***
!   *** until the basis set is found or the end of file is reached ***

    bsname = basis_set_name
    symbol = element_symbol

    search_loop: DO

      CALL search(bsname,.TRUE.,found,line)

      IF (found) THEN

        CALL uppercase(line)
        CALL uppercase(symbol)
        CALL uppercase(bsname)

        IF ((INDEX(" "//line//" "," "//TRIM(symbol)//" ") > 0).AND.&
            (INDEX(" "//line//" "," "//TRIM(bsname)//" ") > 0)) THEN

          NULLIFY (gcc,l,lmax,lmin,n,npgf,nshell,zet)

!         *** Read the basis set information ***

          CALL read_object(nset,skip_lines=1)

          npgf => reallocate(npgf,1,nset)
          nshell => reallocate(nshell,1,nset)
          lmax => reallocate(lmax,1,nset)
          lmin => reallocate(lmin,1,nset)
          n => reallocate(n,1,1,1,nset)

          maxl = 0
          maxpgf = 0
          maxshell = 0

          DO iset=1,nset
            CALL read_object(n(1,iset),skip_lines=1)
            CALL read_object(lmin(iset))
            CALL read_object(lmax(iset))
            maxl = MAX(maxl,lmax(iset))
            CALL read_object(npgf(iset))
            IF (npgf(iset) > maxpgf) THEN
              maxpgf = npgf(iset)
              zet => reallocate(zet,1,maxpgf,1,nset)
              gcc => reallocate(gcc,1,maxpgf,1,maxshell,1,nset)
            END IF
            nshell(iset) = 0
            DO lshell=lmin(iset),lmax(iset)
              nmin = n(1,iset) + lshell - lmin(iset)
              CALL read_object(ishell)
              nshell(iset) = nshell(iset) + ishell
              IF (nshell(iset) > maxshell) THEN
                maxshell = nshell(iset)
                n => reallocate(n,1,maxshell,1,nset)
                l => reallocate(l,1,maxshell,1,nset)
                gcc => reallocate(gcc,1,maxpgf,1,maxshell,1,nset)
              END IF
              DO i=1,ishell
                n(nshell(iset)-ishell+i,iset) = nmin + i - 1
                l(nshell(iset)-ishell+i,iset) = lshell
              END DO
            END DO
            DO ipgf=1,npgf(iset)
              CALL read_object(zet(ipgf,iset),skip_lines=1)
              DO ishell=1,nshell(iset)
                CALL read_object(gcc(ipgf,ishell,iset))
              END DO
            END DO
          END DO

!         *** Maximum angular momentum quantum number of the atomic kind ***

          maxlcgf = MAX(maxlcgf,maxl)

          CALL init_orbital_pointers(maxlcgf+maxder)

!         *** Allocate the global variables ***

          gto_basis_set%nset = nset
          gto_basis_set%lmax => reallocate(gto_basis_set%lmax,1,nset)
          gto_basis_set%lmin => reallocate(gto_basis_set%lmin,1,nset)
          gto_basis_set%npgf => reallocate(gto_basis_set%npgf,1,nset)
          gto_basis_set%nshell => reallocate(gto_basis_set%nshell,1,nset)
          gto_basis_set%n => reallocate(gto_basis_set%n,1,maxshell,1,nset)
          gto_basis_set%l => reallocate(gto_basis_set%l,1,maxshell,1,nset)
          gto_basis_set%zet => reallocate(gto_basis_set%zet,1,maxpgf,1,nset)
          gto_basis_set%gcc =>&
            reallocate(gto_basis_set%gcc,1,maxpgf,1,maxshell,1,nset)

!         *** Copy the basis set information into the global variables ***

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

!         *** Initialize the depending atomic kind information ***

          gto_basis_set%set_radius =>&
            reallocate(gto_basis_set%set_radius,1,nset)
          gto_basis_set%pgf_radius =>&
            reallocate(gto_basis_set%pgf_radius,1,maxpgf,1,nset)
          gto_basis_set%first_cgf =>&
            reallocate(gto_basis_set%first_cgf,1,maxshell,1,nset)
          gto_basis_set%first_sgf =>&
            reallocate(gto_basis_set%first_sgf,1,maxshell,1,nset)
          gto_basis_set%last_cgf =>&
            reallocate(gto_basis_set%last_cgf,1,maxshell,1,nset)
          gto_basis_set%last_sgf =>&
            reallocate(gto_basis_set%last_sgf,1,maxshell,1,nset)

          maxco = 0
          ncgf = 0
          nsgf = 0

          DO iset=1,nset
            DO ishell=1,nshell(iset)
              lshell = gto_basis_set%l(ishell,iset)
              gto_basis_set%first_cgf(ishell,iset) = ncgf + 1
              ncgf = ncgf + nco(lshell)
              gto_basis_set%last_cgf(ishell,iset) = ncgf
              gto_basis_set%first_sgf(ishell,iset) = nsgf + 1
              nsgf = nsgf + nso(lshell)
              gto_basis_set%last_sgf(ishell,iset) = nsgf
            END DO
            maxco = MAX(maxco,npgf(iset)*ncoset(lmax(iset)))
          END DO

          gto_basis_set%ncgf = ncgf
          gto_basis_set%nsgf = nsgf

          gto_basis_set%cphi => reallocate(gto_basis_set%cphi,1,maxco,1,ncgf)
          gto_basis_set%sphi => reallocate(gto_basis_set%sphi,1,maxco,1,nsgf)
          gto_basis_set%lx => reallocate(gto_basis_set%lx,1,ncgf)
          gto_basis_set%ly => reallocate(gto_basis_set%ly,1,ncgf)
          gto_basis_set%lz => reallocate(gto_basis_set%lz,1,ncgf)
          gto_basis_set%m => reallocate(gto_basis_set%m,1,nsgf)
          gto_basis_set%norm_cgf => reallocate(gto_basis_set%norm_cgf,1,ncgf)
          ALLOCATE (gto_basis_set%cgf_symbol(ncgf),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"gto_basis_set%cgf_symbol",12*ncgf)
          END IF
          ALLOCATE (gto_basis_set%sgf_symbol(nsgf),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine,"gto_basis_set%sgf_symbol",6*nsgf)
          END IF

          ncgf = 0
          nsgf = 0

          DO iset=1,nset
            DO ishell=1,nshell(iset)
              lshell = gto_basis_set%l(ishell,iset)
              DO lx=lshell,0,-1
                DO ly=lshell-lx,0,-1
                  lz = lshell - lx - ly
                  ncgf = ncgf + 1
                  gto_basis_set%lx(ncgf) = lx
                  gto_basis_set%ly(ncgf) = ly
                  gto_basis_set%lz(ncgf) = lz
                  gto_basis_set%cgf_symbol(ncgf) =&
                    cgf_symbol(n(ishell,iset),(/lx,ly,lz/))
                END DO
              END DO
              DO m=-lshell,lshell
                nsgf = nsgf + 1
                gto_basis_set%m(nsgf) = m
                gto_basis_set%sgf_symbol(nsgf) =&
                  sgf_symbol(n(ishell,iset),lshell,m)
              END DO
            END DO
          END DO

          DEALLOCATE (gcc,l,lmax,lmin,n,npgf,nshell,zet)

          EXIT search_loop

        END IF

      ELSE

!       *** Stop program, if the end of file is reached ***

        CALL stop_program(routine,&
                          "The requested basis set <"//&
                          TRIM(basis_set_name)//&
                          "> was not found in the basis set file <"//&
                          TRIM(globenv%basis_set_file_name)//">")

      END IF

    END DO search_loop

    CALL finish_input_session

  END SUBROUTINE read_gto_basis_set

! *****************************************************************************

  SUBROUTINE write_gto_basis_set(gto_basis_set,lunit)

!   Purpose: Write the Gaussian-type orbital basis set data read from the
!            database file to the logical unit number lunit.

!   History: - Creation (06.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(gto_basis_set_type), POINTER :: gto_basis_set
    INTEGER, INTENT(IN)               :: lunit

!   *** Local variables ***

    INTEGER :: ipgf,iset,ishell

!   ---------------------------------------------------------------------------

    DO iset=1,gto_basis_set%nset
      WRITE (lunit,"(A)") ""
      DO ishell=1,gto_basis_set%nshell(iset)
        WRITE (lunit,"(T6,I3,4X,I4,4X,I2,2X,I2,(T32,2F15.6))")&
          iset,ishell,&
          gto_basis_set%n(ishell,iset),&
          gto_basis_set%l(ishell,iset),&
          (gto_basis_set%zet(ipgf,iset),gto_basis_set%gcc(ipgf,ishell,iset),&
           ipgf=1,gto_basis_set%npgf(iset))
      END DO
    END DO

  END SUBROUTINE write_gto_basis_set

! *****************************************************************************

END MODULE basis_set_types
